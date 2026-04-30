#include "config.h"
#include "fragmenter.h"
#include "matrix.h"
#include "reference.h"
#include "simulator.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

double estimate_empirical_fraction(const std::vector<OffsetEntry> &source_offsets,
                                   const std::vector<OffsetEntry> &target_offsets) {
    if (target_offsets.empty()) {
        return 0.0;
    }
    if (source_offsets.empty()) {
        return 1.0;
    }

    std::unordered_map<std::string, std::size_t> target_bins_by_name;
    std::size_t total_target_bins = 0;
    for (const auto &offset : target_offsets) {
        const std::size_t span = offset.end_bin > offset.start_bin ? offset.end_bin - offset.start_bin : 0;
        target_bins_by_name[offset.contig] = span;
        total_target_bins += span;
    }
    if (total_target_bins == 0) {
        return 0.0;
    }

    std::size_t matched_bins = 0;
    for (const auto &offset : source_offsets) {
        const auto it = target_bins_by_name.find(offset.contig);
        if (it != target_bins_by_name.end()) {
            matched_bins += it->second;
        }
    }

    const double fraction = static_cast<double>(matched_bins) / static_cast<double>(total_target_bins);
    return std::max(0.0, std::min(1.0, fraction));
}

struct ContactMatrixStats {
    std::size_t cis_contacts = 0;
    std::size_t trans_contacts = 0;
    double cis_weight = 0.0;
    double trans_weight = 0.0;
};

ContactMatrixStats summarize_matrix(const ContactMatrix &matrix,
                                    const std::vector<OffsetEntry> &offsets) {
    std::vector<int> bin_to_contig(matrix.bin_count, -1);
    for (std::size_t i = 0; i < offsets.size(); ++i) {
        const std::size_t end_bin = std::min(offsets[i].end_bin, matrix.bin_count);
        for (std::size_t bin = offsets[i].start_bin; bin < end_bin; ++bin) {
            bin_to_contig[bin] = static_cast<int>(i);
        }
    }

    ContactMatrixStats stats;
    for (const auto &contact : matrix.contacts) {
        if (contact.bin1 >= matrix.bin_count || contact.bin2 >= matrix.bin_count ||
            !std::isfinite(contact.weight) || contact.weight <= 0.0) {
            continue;
        }
        const int contig1 = bin_to_contig[contact.bin1];
        const int contig2 = bin_to_contig[contact.bin2];
        if (contig1 >= 0 && contig1 == contig2) {
            ++stats.cis_contacts;
            stats.cis_weight += contact.weight;
        } else {
            ++stats.trans_contacts;
            stats.trans_weight += contact.weight;
        }
    }
    return stats;
}

void print_usage() {
    std::cerr 
        << "Usage:\n"
        << "  hicreate ref.fa 10000 [options]\n"
        << "  hicreate --reference ref.fa --bin-size 10000 [options]\n\n"
        << "Required positional arguments:\n"
        << "  ref.fa                 Reference FASTA. No --reference flag is required.\n"
        << "  10000                  Genomic bin size. No --bin-size flag is required.\n\n"
        << "Long-form compatibility:\n"
        << "  -r, --reference FILE   Reference FASTA.\n"
        << "  -b, --bin-size N       Genomic bin size used by the contact matrix.\n\n"
        << "Read count options (choose one):\n"
        << "  -c, --coverage X       Target genome depth. Read pairs are computed as\n"
        << "                         ceil(X * reference_bases / 300) for PE150 reads.\n"
        << "      --depth X          Alias for --coverage.\n"
        << "  -p, --pairs N          Exact number of 150 bp paired-end read pairs.\n"
        << "                         Default: 100000 when --coverage is omitted.\n"
        << "                         --coverage and --pairs cannot be used together.\n\n"
        << "Output and reproducibility:\n"
        << "  -o, --output-prefix PREFIX\n"
        << "                         Prefix for output FASTQ files. Default: sim\n"
        << "  -s, --seed N           Random seed. Default: 1\n"
        << "  -j, --threads N        Worker threads for read generation. Default: 1\n"
        << "                         Use 0 to auto-detect hardware threads.\n\n"
        << "Reference digestion:\n"
        << "  -e, --enzyme-site SEQ  Restriction enzyme motif. Default: AAGCTT\n"
        << "                         Common cut offsets are recognized for HindIII (AAGCTT)\n"
        << "                         and DpnII/MboI (GATC).\n\n"
        << "Optional input matrix:\n"
        << "  -m, --matrix FILE      Sparse or dense Hi-C contact matrix. If omitted,\n"
        << "                         a synthetic matrix is generated from the reference.\n"
        << "  -f, --offset FILE      Optional matrix bin-to-contig mapping.\n"
        << "                         Format: contig start_bin end_bin\n"
        << "  Sparse matrix format:  bin1 bin2 value\n"
        << "  Dense matrix format:   headerless square numeric matrix\n\n"
        << "Synthetic contact model:\n"
        << "  -t, --trans-ratio X    Fraction of trans-chromosomal contact mass.\n"
        << "                         Default: species-aware auto, human/auto uses 0.12\n"
        << "  --synthetic-contacts N Number of sparse contacts in synthetic matrix. Default: auto\n"
        << "  --cis-decay-alpha X    Cis distance-decay exponent. Default: 1.0\n"
        << "  --min-cis-distance-bins N\n"
        << "                         Minimum cis bin separation for synthetic contacts.\n"
        << "                         Default: 0, preserving strong same-bin contacts.\n"
        << "  --max-cis-distance-bins N\n"
        << "                         Maximum cis separation sampled for synthetic contacts.\n"
        << "                         Default: 200\n"
        << "  -S, --species-model NAME\n"
        << "                         Preset: auto, generic_plant, human, arabidopsis,\n"
        << "                         rice, maize, wheat, barley. Default: auto\n"
        << "  -A, --arrangement-model M\n"
        << "                         auto, territory, rabl, rosette, nonrabl. Default: auto\n"
        << "  -T, --trans-model M    auto, territory, random, telomere, centromere,\n"
        << "                         compartment, hubs.\n"
        << "                         Default: auto\n"
        << "  --trans-hotspots N     Number of hub bins for --trans-model hubs. Default: 8\n"
        << "  --collision-randomness X\n"
        << "                         Mix between random collision and arrangement effects.\n"
        << "                         Range: 0..1, default: 0.35\n\n"
        << "Output:\n"
        << "  PREFIX_R1.fastq and PREFIX_R2.fastq, always 150 bp paired-end reads.\n";
}

Config parse_args(int argc, char **argv) {
    Config cfg;
    std::vector<std::string> positional_args;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto require_value = [&](const std::string &name) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value for " + name);
            }
            return argv[++i];
        };

        if (arg == "--reference" || arg == "-r") {
            cfg.reference_path = require_value(arg);
        } else if (arg == "--matrix" || arg == "-m") {
            cfg.matrix_path = require_value(arg);
        } else if (arg == "--output-prefix" || arg == "-o") {
            cfg.output_prefix = require_value(arg);
        } else if (arg == "--offset" || arg == "-f") {
            cfg.offset_path = require_value(arg);
        } else if (arg == "--enzyme-site" || arg == "-e") {
            cfg.enzyme_site = require_value(arg);
        } else if (arg == "--bin-size" || arg == "-b") {
            cfg.bin_size = std::stoull(require_value(arg));
        } else if (arg == "--read-length") {
            cfg.read_length = std::stoull(require_value(arg));
        } else if (arg == "--pairs" || arg == "-p") {
            cfg.pair_count = std::stoull(require_value(arg));
            cfg.pair_count_explicit = true;
        } else if (arg == "--coverage" || arg == "--depth" || arg == "-c") {
            cfg.coverage_depth = std::stod(require_value(arg));
        } else if (arg == "--seed" || arg == "-s") {
            cfg.seed = std::stoull(require_value(arg));
        } else if (arg == "--threads" || arg == "-j") {
            cfg.thread_count = std::stoull(require_value(arg));
        } else if (arg == "--trans-ratio" || arg == "-t") {
            cfg.trans_ratio = std::stod(require_value(arg));
            cfg.trans_ratio_explicit = true;
        } else if (arg == "--synthetic-contacts") {
            cfg.synthetic_contact_count = std::stoull(require_value(arg));
        } else if (arg == "--cis-decay-alpha") {
            cfg.cis_decay_alpha = std::stod(require_value(arg));
        } else if (arg == "--min-cis-distance-bins") {
            cfg.min_cis_distance_bins = std::stoull(require_value(arg));
        } else if (arg == "--max-cis-distance-bins") {
            cfg.max_cis_distance_bins = std::stoull(require_value(arg));
        } else if (arg == "--species-model" || arg == "-S") {
            cfg.species_model = require_value(arg);
        } else if (arg == "--arrangement-model" || arg == "-A") {
            cfg.arrangement_model = require_value(arg);
        } else if (arg == "--trans-model" || arg == "-T") {
            cfg.trans_model = require_value(arg);
        } else if (arg == "--trans-hotspots") {
            cfg.trans_hotspots = std::stoull(require_value(arg));
        } else if (arg == "--collision-randomness") {
            cfg.collision_randomness = std::stod(require_value(arg));
        } else if (arg == "--help" || arg == "-h") {
            print_usage();
            std::exit(0);
        } else if (!arg.empty() && arg[0] != '-') {
            positional_args.push_back(arg);
        } else {
            throw std::runtime_error("Unknown argument: " + arg);
        }
    }

    if (!positional_args.empty()) {
        if (positional_args.size() > 2) {
            throw std::runtime_error("Too many positional arguments. Use: hicreate ref.fa bin_size [options].");
        }
        for (const auto &value : positional_args) {
            if (cfg.reference_path.empty()) {
                cfg.reference_path = value;
                continue;
            }
            if (cfg.bin_size == 0) {
                cfg.bin_size = std::stoull(value);
                continue;
            }
            throw std::runtime_error("Unexpected positional argument: " + value);
        }
    }

    if (cfg.reference_path.empty() || cfg.bin_size == 0) {
        throw std::runtime_error("Missing required arguments.");
    }
    if (cfg.output_prefix.empty()) {
        throw std::runtime_error("--output-prefix must not be empty.");
    }
    if (cfg.read_length == 0) {
        throw std::runtime_error("--read-length must be positive.");
    }
    if (cfg.read_length != 150) {
        throw std::runtime_error("Only 150 bp paired-end reads are supported.");
    }
    if (cfg.coverage_depth < 0.0) {
        throw std::runtime_error("--coverage must be non-negative.");
    }
    if (cfg.coverage_depth > 0.0 && cfg.pair_count_explicit) {
        throw std::runtime_error("Use either --coverage or --pairs, not both.");
    }
    if (cfg.coverage_depth == 0.0 && cfg.pair_count == 0) {
        throw std::runtime_error("--pairs must be positive when --coverage is omitted.");
    }
    if (cfg.thread_count == 0) {
        cfg.thread_count = std::max<unsigned int>(1, std::thread::hardware_concurrency());
    }
    if (cfg.trans_ratio < 0.0 || cfg.trans_ratio > 1.0) {
        throw std::runtime_error("--trans-ratio must be within [0, 1].");
    }
    if (cfg.cis_decay_alpha <= 0.0) {
        throw std::runtime_error("--cis-decay-alpha must be positive.");
    }
    if (cfg.min_cis_distance_bins > cfg.max_cis_distance_bins) {
        throw std::runtime_error("--min-cis-distance-bins cannot exceed --max-cis-distance-bins.");
    }
    if (cfg.collision_randomness < 0.0 || cfg.collision_randomness > 1.0) {
        throw std::runtime_error("--collision-randomness must be within [0, 1].");
    }

    return cfg;
}

}  // namespace

int main(int argc, char **argv) {
    try {
        Config cfg = parse_args(argc, argv);
        std::cerr << "Loading reference: " << cfg.reference_path << '\n';
        const ReferenceGenome reference = load_reference_fasta(cfg.reference_path);
        std::cerr << "Reference loaded: " << reference.contigs.size()
                  << " contigs, " << reference.total_length() << " bp\n";
        if (cfg.read_length >= reference.total_length()) {
            throw std::runtime_error("Read length must be shorter than the total reference length.");
        }
        if (cfg.coverage_depth > 0.0) {
            const double requested_pairs =
                std::ceil(cfg.coverage_depth * static_cast<double>(reference.total_length()) /
                          static_cast<double>(2 * cfg.read_length));
            cfg.pair_count = std::max<std::size_t>(1, static_cast<std::size_t>(requested_pairs));
        }

        const std::size_t total_bins = reference.total_bins(cfg.bin_size);
        if (total_bins == 0) {
            throw std::runtime_error("Reference produced zero bins.");
        }
        std::cerr << "Global bins: " << total_bins << '\n';

        const std::vector<OffsetEntry> reference_offsets = build_reference_offsets(reference, cfg.bin_size);
        std::vector<OffsetEntry> source_offsets;
        if (!cfg.offset_path.empty()) {
            source_offsets = load_offsets(cfg.offset_path);
        }

        ContactMatrix matrix;
        SyntheticModelOptions model_options;
        model_options.trans_ratio = cfg.trans_ratio;
        model_options.trans_ratio_explicit = cfg.trans_ratio_explicit;
        model_options.cis_decay_alpha = cfg.cis_decay_alpha;
        model_options.min_cis_distance_bins = cfg.min_cis_distance_bins;
        model_options.max_cis_distance_bins = cfg.max_cis_distance_bins;
        model_options.species_model = cfg.species_model;
        model_options.arrangement_model = cfg.arrangement_model;
        model_options.trans_model = cfg.trans_model;
        model_options.trans_hotspots = cfg.trans_hotspots;
        model_options.collision_randomness = cfg.collision_randomness;
        model_options.seed = cfg.seed;

        if (!cfg.matrix_path.empty()) {
            std::cerr << "Loading and remapping input matrix: " << cfg.matrix_path << '\n';
            const ContactMatrix source_matrix = load_matrix(cfg.matrix_path, 0);
            const ContactMatrix remapped_matrix =
                remap_matrix_to_reference(source_matrix, source_offsets, reference_offsets);
            const double empirical_fraction = source_offsets.empty()
                                                  ? 1.0
                                                  : std::max(0.15, estimate_empirical_fraction(source_offsets, reference_offsets));

            if (empirical_fraction < 0.999) {
                std::size_t synthetic_contacts = cfg.synthetic_contact_count;
                if (synthetic_contacts == 0) {
                    synthetic_contacts = std::max<std::size_t>(
                        20000, std::min<std::size_t>(2000000, source_matrix.contacts.size() * 8 + total_bins * 10));
                }
                const ContactMatrix synthetic_fill =
                    generate_synthetic_matrix(total_bins, reference_offsets, synthetic_contacts, model_options);
                matrix = blend_contact_matrices(remapped_matrix, synthetic_fill, empirical_fraction);
            } else {
                matrix = remapped_matrix;
            }

            matrix = apply_trans_ratio(matrix, reference_offsets, cfg.trans_ratio, cfg.seed);
        } else {
            std::size_t synthetic_contacts = cfg.synthetic_contact_count;
            if (synthetic_contacts == 0) {
                synthetic_contacts = std::max<std::size_t>(
                    20000, std::min<std::size_t>(2000000, total_bins * 10));
            }
            std::cerr << "Generating synthetic contact matrix with "
                      << synthetic_contacts << " sampled contacts...\n";
            matrix = generate_synthetic_matrix(total_bins,
                                              reference_offsets,
                                              synthetic_contacts,
                                              model_options);
        }
        std::cerr << "Contact matrix ready: " << matrix.contacts.size() << " sparse contacts\n";
        const ContactMatrixStats matrix_stats = summarize_matrix(matrix, reference_offsets);
        const double matrix_weight_total = matrix_stats.cis_weight + matrix_stats.trans_weight;
        if (matrix_weight_total > std::numeric_limits<double>::min()) {
            std::cerr << "Matrix cis/trans: "
                      << matrix_stats.cis_contacts << " cis contacts, "
                      << matrix_stats.trans_contacts << " trans contacts; "
                      << "weight fractions cis="
                      << (matrix_stats.cis_weight / matrix_weight_total)
                      << ", trans="
                      << (matrix_stats.trans_weight / matrix_weight_total)
                      << '\n';
        }
        std::cerr << "Streaming " << cfg.pair_count << " read pairs to FASTQ with "
                  << cfg.thread_count << " worker thread(s)...\n";
        PairedReadWriter writer(cfg);
        write_hic_reads(cfg, reference, reference_offsets, matrix, writer);

        std::cout << "Contigs: " << reference.contigs.size() << '\n'
                  << "Reference length: " << reference.total_length() << '\n'
                  << "Global bins: " << total_bins << '\n'
                  << "Contacts loaded: " << matrix.contacts.size() << '\n'
                  << "Read pairs: " << writer.count() << '\n'
                  << "Coverage: " << (static_cast<double>(writer.count() * 2 * cfg.read_length) /
                                      static_cast<double>(reference.total_length())) << "x\n"
                  << "Read length: " << cfg.read_length << '\n'
                  << "Output FASTQ R1: " << writer.read1_path() << '\n'
                  << "Output FASTQ R2: " << writer.read2_path() << '\n';
    } catch (const std::exception &ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        print_usage();
        return 1;
    }

    return 0;
}
