#include "config.h"
#include "fragmenter.h"
#include "matrix.h"
#include "reference.h"
#include "simulator.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <stdexcept>
#include <string>

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

void print_usage() {
    std::cerr 
        << "Usage:\n"
        << "  hicreate --reference ref.fa [--matrix matrix.tsv] --bin-size 10000\n"
        << "           [--coverage 30] [--pairs 100000] [--output-prefix sim]\n"
        << "           [--offset offsets.tsv] [--enzyme-site AAGCTT]\n"
        << "           [--seed 123]\n"
        << "           [--trans-ratio 0.10] [--synthetic-contacts 200000]\n"
        << "           [--cis-decay-alpha 1.0] [--max-cis-distance-bins 200]\n"
        << "           [--species-model generic_plant] [--arrangement-model auto]\n"
        << "           [--trans-model auto] [--trans-hotspots 8]\n"
        << "           [--collision-randomness 0.35]\n\n"
        << "Output reads:\n"
        << "  Always writes 150 bp paired-end FASTQ.\n"
        << "  --coverage X sets read pairs to ceil(X * reference_bases / 300).\n"
        << "  --pairs is the exact number of read pairs when --coverage is omitted.\n\n"
        << "Input matrix:\n"
        << "  Sparse: bin1 bin2 value\n"
        << "  Dense: headerless square numeric matrix\n\n"
        << "Offset format:\n"
        << "  contig start_bin end_bin\n";
}

Config parse_args(int argc, char **argv) {
    Config cfg;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto require_value = [&](const std::string &name) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value for " + name);
            }
            return argv[++i];
        };

        if (arg == "--reference") {
            cfg.reference_path = require_value(arg);
        } else if (arg == "--matrix") {
            cfg.matrix_path = require_value(arg);
        } else if (arg == "--output-prefix") {
            cfg.output_prefix = require_value(arg);
        } else if (arg == "--offset") {
            cfg.offset_path = require_value(arg);
        } else if (arg == "--enzyme-site") {
            cfg.enzyme_site = require_value(arg);
        } else if (arg == "--bin-size") {
            cfg.bin_size = std::stoull(require_value(arg));
        } else if (arg == "--read-length") {
            cfg.read_length = std::stoull(require_value(arg));
        } else if (arg == "--pairs") {
            cfg.pair_count = std::stoull(require_value(arg));
            cfg.pair_count_explicit = true;
        } else if (arg == "--coverage" || arg == "--depth") {
            cfg.coverage_depth = std::stod(require_value(arg));
        } else if (arg == "--seed") {
            cfg.seed = std::stoull(require_value(arg));
        } else if (arg == "--trans-ratio") {
            cfg.trans_ratio = std::stod(require_value(arg));
        } else if (arg == "--synthetic-contacts") {
            cfg.synthetic_contact_count = std::stoull(require_value(arg));
        } else if (arg == "--cis-decay-alpha") {
            cfg.cis_decay_alpha = std::stod(require_value(arg));
        } else if (arg == "--max-cis-distance-bins") {
            cfg.max_cis_distance_bins = std::stoull(require_value(arg));
        } else if (arg == "--species-model") {
            cfg.species_model = require_value(arg);
        } else if (arg == "--arrangement-model") {
            cfg.arrangement_model = require_value(arg);
        } else if (arg == "--trans-model") {
            cfg.trans_model = require_value(arg);
        } else if (arg == "--trans-hotspots") {
            cfg.trans_hotspots = std::stoull(require_value(arg));
        } else if (arg == "--collision-randomness") {
            cfg.collision_randomness = std::stod(require_value(arg));
        } else if (arg == "--help" || arg == "-h") {
            print_usage();
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown argument: " + arg);
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
    if (cfg.trans_ratio < 0.0 || cfg.trans_ratio > 1.0) {
        throw std::runtime_error("--trans-ratio must be within [0, 1].");
    }
    if (cfg.cis_decay_alpha <= 0.0) {
        throw std::runtime_error("--cis-decay-alpha must be positive.");
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
        const ReferenceGenome reference = load_reference_fasta(cfg.reference_path);
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

        const std::vector<OffsetEntry> reference_offsets = build_reference_offsets(reference, cfg.bin_size);
        std::vector<OffsetEntry> source_offsets;
        if (!cfg.offset_path.empty()) {
            source_offsets = load_offsets(cfg.offset_path);
        }

        ContactMatrix matrix;
        SyntheticModelOptions model_options;
        model_options.trans_ratio = cfg.trans_ratio;
        model_options.cis_decay_alpha = cfg.cis_decay_alpha;
        model_options.max_cis_distance_bins = cfg.max_cis_distance_bins;
        model_options.species_model = cfg.species_model;
        model_options.arrangement_model = cfg.arrangement_model;
        model_options.trans_model = cfg.trans_model;
        model_options.trans_hotspots = cfg.trans_hotspots;
        model_options.collision_randomness = cfg.collision_randomness;
        model_options.seed = cfg.seed;

        if (!cfg.matrix_path.empty()) {
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
                    20000, std::min<std::size_t>(2000000, cfg.pair_count * 20));
            }
            matrix = generate_synthetic_matrix(total_bins,
                                              reference_offsets,
                                              synthetic_contacts,
                                              model_options);
        }
        const auto read_templates = create_hic_read_templates(cfg, reference, reference_offsets, matrix);
        simulate_paired_reads(cfg, read_templates);

        std::cout << "Contigs: " << reference.contigs.size() << '\n'
                  << "Reference length: " << reference.total_length() << '\n'
                  << "Global bins: " << total_bins << '\n'
                  << "Contacts loaded: " << matrix.contacts.size() << '\n'
                  << "Read pairs: " << read_templates.size() << '\n'
                  << "Coverage: " << (static_cast<double>(read_templates.size() * 2 * cfg.read_length) /
                                      static_cast<double>(reference.total_length())) << "x\n"
                  << "Read length: " << cfg.read_length << '\n';
    } catch (const std::exception &ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        print_usage();
        return 1;
    }

    return 0;
}
