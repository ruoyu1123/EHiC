#include "config.h"
#include "cifi.h"
#include "fragmenter.h"
#include "matrix.h"
#include "porec.h"
#include "reference.h"
#include "simulator.h"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

std::string to_lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return value;
}

AssayType parse_assay_type(const std::string &value) {
    const std::string assay = to_lower(value);
    if (assay == "hic" || assay == "hi-c") {
        return AssayType::Hic;
    }
    if (assay == "porec" || assay == "pore-c") {
        return AssayType::PoreC;
    }
    if (assay == "cifi" || assay == "ci-fi") {
        return AssayType::CiFi;
    }
    throw std::runtime_error("--assay must be one of: hic, porec, cifi.");
}

std::string assay_type_to_string(AssayType assay) {
    switch (assay) {
        case AssayType::Hic: return "hic";
        case AssayType::PoreC: return "porec";
        case AssayType::CiFi: return "cifi";
    }
    return "unknown";
}

MatrixFormat parse_matrix_format(const std::string &value) {
    const std::string format = to_lower(value);
    if (format == "auto") {
        return MatrixFormat::Auto;
    }
    if (format == "sparse") {
        return MatrixFormat::Sparse;
    }
    if (format == "dense") {
        return MatrixFormat::Dense;
    }
    if (format == "binary" || format == "binary-sparse" || format == "hcm") {
        return MatrixFormat::BinarySparse;
    }
    throw std::runtime_error("--matrix-format must be one of: auto, sparse, dense, binary.");
}

std::string matrix_format_to_string(MatrixFormat format) {
    switch (format) {
        case MatrixFormat::Auto:
            return "auto";
        case MatrixFormat::Sparse:
            return "sparse";
        case MatrixFormat::Dense:
            return "dense";
        case MatrixFormat::BinarySparse:
            return "binary";
    }
    return "unknown";
}

struct EmpiricalModelSpec {
    std::string name;
    std::string matrix_path;
    std::string offset_path;
    MatrixFormat matrix_format = MatrixFormat::Auto;
};

std::vector<std::string> split_manifest_tokens(const std::string &line) {
    std::string normalized;
    normalized.reserve(line.size());
    for (char ch : line) {
        normalized.push_back(ch == '=' || ch == ',' ? ' ' : ch);
    }
    std::istringstream stream(normalized);
    std::vector<std::string> tokens;
    std::string token;
    while (stream >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

EmpiricalModelSpec load_empirical_model_spec(const std::string &model_name,
                                             const std::string &model_dir) {
    namespace fs = std::filesystem;
    fs::path model_path(model_name);
    if (!fs::exists(model_path)) {
        model_path = fs::path(model_dir) / model_name;
    }
    if (!fs::exists(model_path)) {
        throw std::runtime_error("Empirical model not found: " + model_name);
    }

    fs::path manifest_path = fs::is_directory(model_path)
                                 ? model_path / "manifest.tsv"
                                 : model_path;
    if (!fs::exists(manifest_path)) {
        throw std::runtime_error("Empirical model manifest not found: " + manifest_path.string());
    }

    EmpiricalModelSpec spec;
    spec.name = model_path.stem().string();
    std::ifstream in(manifest_path);
    if (!in) {
        throw std::runtime_error("Failed to open empirical model manifest: " + manifest_path.string());
    }

    std::string line;
    while (std::getline(in, line)) {
        const auto tokens = split_manifest_tokens(line);
        if (tokens.empty() || tokens[0].empty() || tokens[0][0] == '#') {
            continue;
        }
        if (tokens.size() < 2) {
            throw std::runtime_error("Model manifest rows must contain key and value.");
        }
        const std::string key = to_lower(tokens[0]);
        const std::string value = tokens[1];
        if (key == "name") {
            spec.name = value;
        } else if (key == "matrix") {
            spec.matrix_path = (manifest_path.parent_path() / value).string();
        } else if (key == "offset") {
            spec.offset_path = (manifest_path.parent_path() / value).string();
        } else if (key == "format" || key == "matrix-format") {
            spec.matrix_format = parse_matrix_format(value);
        }
    }

    if (spec.matrix_path.empty() || spec.offset_path.empty()) {
        throw std::runtime_error("Empirical model manifest requires matrix and offset entries.");
    }
    return spec;
}

std::string infer_empirical_model(const Config &cfg) {
    namespace fs = std::filesystem;
    if (!cfg.empirical_model.empty()) {
        return cfg.empirical_model;
    }

    const std::string species = to_lower(cfg.species_model.empty() ? "auto" : cfg.species_model);
    if (species == "auto" || species == "human" || species == "homo_sapiens" ||
        species == "mammal" || species == "chm13") {
        return "human_cell_40mb";
    }

    const fs::path direct_model = fs::path(cfg.model_dir) / cfg.species_model;
    if (fs::exists(direct_model / "manifest.tsv") || fs::exists(direct_model)) {
        return cfg.species_model;
    }

    throw std::runtime_error("No empirical matrix model preset is linked to --species-model " +
                             cfg.species_model +
                             ". Provide --empirical-model or add models/<species>/manifest.tsv.");
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

ContactMatrix align_matrix_to_reference(const ContactMatrix &source_matrix,
                                        const std::vector<OffsetEntry> &source_offsets,
                                        const std::vector<OffsetEntry> &reference_offsets,
                                        std::uint64_t seed,
                                        const std::string &label,
                                        bool force_contig_reuse) {
    const std::size_t reference_bin_count =
        reference_offsets.empty() ? 0 : reference_offsets.back().end_bin;
    std::cerr << "Preparing " << label
              << ": source bins=" << source_matrix.bin_count
              << ", source offsets=" << source_offsets.size()
              << ", target offsets=" << reference_offsets.size()
              << ", target bins=" << reference_bin_count << '\n';

    if (source_offsets.empty()) {
        std::cerr << "No offset provided for " << label
                  << "; falling back to global bin scaling.\n";
        std::vector<std::string> remap_warnings;
        ContactMatrix remapped = remap_matrix_to_reference(source_matrix,
                                                           source_offsets,
                                                           reference_offsets,
                                                           seed,
                                                           &remap_warnings,
                                                           false);
        for (const auto &warning : remap_warnings) {
            std::cerr << "Warning: " << warning << '\n';
        }
        return remapped;
    }

    if (source_offsets.size() < reference_offsets.size() && !force_contig_reuse) {
        throw std::runtime_error(label + " offset has " +
                                 std::to_string(source_offsets.size()) +
                                 " contigs, but the target reference has " +
                                 std::to_string(reference_offsets.size()) +
                                 ". Refusing to reuse source contigs because trans contacts "
                                 "between duplicated templates cannot be inferred. Provide an "
                                 "offset with at least as many contigs as the target reference, "
                                 "or explicitly enable --force-contig-reuse.");
    }
    if (source_offsets.size() < reference_offsets.size()) {
        std::cerr << "Warning: forcing source contig reuse for " << label
                  << "; duplicate-target trans blocks will be synthesized from real source "
                     "trans contacts and normalized.\n";
    }

    std::string offset_match_reason;
    if (offsets_match_reference_exactly(source_offsets, reference_offsets, &offset_match_reason)) {
        if (source_matrix.bin_count != reference_offsets.back().end_bin) {
            std::cerr << "Using exact offset match for " << label
                      << ", but matrix bin_count (" << source_matrix.bin_count
                      << ") differs from target bin count (" << reference_offsets.back().end_bin
                      << "); remapping to reconcile bin count.\n";
        } else {
            std::cerr << "Offset layout for " << label
                      << " matches target reference exactly; skipping resize/remap.\n";
            return source_matrix;
        }
    } else {
        std::cerr << "Resizing/remapping " << label << " because "
                  << offset_match_reason << ".\n";
    }

    std::vector<std::string> remap_warnings;
    ContactMatrix remapped = remap_matrix_to_reference(source_matrix,
                                                       source_offsets,
                                                       reference_offsets,
                                                       seed,
                                                       &remap_warnings,
                                                       force_contig_reuse);
    std::cerr << "Completed resize/remap for " << label
              << ": output bins=" << remapped.bin_count
              << ", sparse contacts=" << remapped.contacts.size() << '\n';
    for (const auto &warning : remap_warnings) {
        std::cerr << "Warning: " << warning << '\n';
    }
    return remapped;
}

void print_usage() {
    std::cerr 
        << "Usage:\n"
        << "  hicreate ref.fa 10000 [options]\n"
        << "  hicreate --reference ref.fa --bin-size 10000 [options]\n\n"
        << "Assay selection:\n"
        << "      --assay TYPE       hic, porec, or cifi. Default: hic\n"
        << "                         Hi-C writes PE150 FASTQ; Pore-C and CiFi write\n"
        << "                         single-molecule concatemer FASTQ plus truth TSV.\n\n"
        << "Required positional arguments:\n"
        << "  ref.fa                 Reference FASTA. No --reference flag is required.\n"
        << "  10000                  Genomic bin size. No --bin-size flag is required.\n\n"
        << "Long-form compatibility:\n"
        << "  -r, --reference FILE   Reference FASTA.\n"
        << "  -b, --bin-size N       Genomic bin size used by the contact matrix.\n\n"
        << "Read count options (choose one):\n"
        << "  -c, --coverage X       Target genome depth. Read pairs are computed as\n"
        << "                         ceil(X * reference_bases / 300) for PE150 reads.\n"
        << "                         Long-read counts use the configured mean length.\n"
        << "      --depth X          Alias for --coverage.\n"
        << "  -p, --pairs N          Exact number of 150 bp paired-end read pairs.\n"
        << "                         Default: 100000 when --coverage is omitted.\n"
        << "      --reads N          Exact number of Pore-C/CiFi long reads.\n"
        << "                         Default: 10000 when --coverage is omitted.\n"
        << "                         --coverage and exact count options cannot be combined.\n\n"
        << "Output and reproducibility:\n"
        << "  -o, --output-prefix PREFIX\n"
        << "                         Prefix for output FASTQ files. Default: sim\n"
        << "  -s, --seed N           Random seed. Default: 1\n"
        << "  -j, --threads N        Worker threads for read generation. Default: 1\n"
        << "                         Use 0 to auto-detect hardware threads.\n\n"
        << "Reference digestion:\n"
        << "  -e, --enzyme-site SEQ  Restriction enzyme motif. Defaults: Hi-C AAGCTT,\n"
        << "                         Pore-C/CiFi GATC.\n"
        << "                         Common cut offsets are recognized for HindIII (AAGCTT)\n"
        << "                         and DpnII/MboI (GATC).\n\n"
        << "Long-read options (Pore-C/CiFi only):\n"
        << "      --long-read-mean N Mean read length override.\n"
        << "      --long-read-min N  Minimum read length override.\n"
        << "      --long-read-max N  Maximum read length override.\n"
        << "      --long-read-sigma X\n"
        << "                         Lognormal read-length sigma override.\n"
        << "      --long-read-qv X   Mean platform QV override.\n"
        << "      --max-segments N   Maximum concatemer segments per read. Default: 256\n\n"
        << "Matrix and model options:\n"
        << "  -m, --matrix FILE      Sparse or dense Hi-C contact matrix. If omitted,\n"
        << "                         the empirical model linked to --species-model is used.\n"
        << "      --matrix-format F  Matrix parser: auto, sparse, dense, or binary. Default: auto\n"
        << "      --empirical-model NAME|MANIFEST\n"
        << "                         Matrix model from --model-dir. Without --matrix, use it\n"
        << "                         as the full matrix. With --matrix, keep input cis and\n"
        << "                         replace trans with this model after remapping.\n"
        << "      --model-dir DIR    Empirical model directory. Default: models\n"
        << "  -S, --species-model NAME\n"
        << "                         Empirical preset selector. Default: auto\n"
        << "                         Installed default: auto/human/chm13 -> human_cell_40mb\n"
        << "  -f, --offset FILE      Optional matrix bin-to-contig mapping.\n"
        << "                         Formats: contig start_bin end_bin, or\n"
        << "                                  name bin_offset length bin_num. Source offsets\n"
        << "                         must contain at least as many contigs as the reference.\n"
        << "      --force-contig-reuse\n"
        << "                         Allow fewer source contigs. Duplicate-target trans blocks\n"
        << "                         are copied from real source trans blocks and normalized.\n"
        << "  Sparse matrix format:  bin1 bin2 value\n"
        << "  Dense matrix format:   headerless square numeric matrix\n\n"
        << "Trans contact control:\n"
        << "  -t, --trans-ratio X    Fraction of trans-chromosomal contact mass.\n"
        << "                         Applied when an empirical model is used for the full\n"
        << "                         matrix or for trans replacement.\n\n"
        << "Output:\n"
        << "  Hi-C:  PREFIX_R1.fastq and PREFIX_R2.fastq (150 bp paired-end)\n"
        << "  Pore-C: PREFIX_porec.fastq and PREFIX_porec_truth.tsv\n"
        << "  CiFi:   PREFIX_cifi.fastq and PREFIX_cifi_truth.tsv\n";
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

        if (arg == "--assay") {
            cfg.assay = parse_assay_type(require_value(arg));
        } else if (arg == "--reference" || arg == "-r") {
            cfg.reference_path = require_value(arg);
        } else if (arg == "--matrix" || arg == "-m") {
            cfg.matrix_path = require_value(arg);
        } else if (arg == "--matrix-format") {
            cfg.matrix_format = parse_matrix_format(require_value(arg));
        } else if (arg == "--empirical-model" || arg == "--model") {
            cfg.empirical_model = require_value(arg);
        } else if (arg == "--model-dir") {
            cfg.model_dir = require_value(arg);
        } else if (arg == "--output-prefix" || arg == "-o") {
            cfg.output_prefix = require_value(arg);
        } else if (arg == "--offset" || arg == "-f") {
            cfg.offset_path = require_value(arg);
        } else if (arg == "--enzyme-site" || arg == "-e") {
            cfg.enzyme_site = require_value(arg);
            cfg.enzyme_site_explicit = true;
        } else if (arg == "--bin-size" || arg == "-b") {
            cfg.bin_size = std::stoull(require_value(arg));
        } else if (arg == "--read-length") {
            cfg.read_length = std::stoull(require_value(arg));
            cfg.read_length_explicit = true;
        } else if (arg == "--pairs" || arg == "-p") {
            cfg.pair_count = std::stoull(require_value(arg));
            cfg.pair_count_explicit = true;
            cfg.pairs_option_used = true;
        } else if (arg == "--reads") {
            cfg.pair_count = std::stoull(require_value(arg));
            cfg.pair_count_explicit = true;
            cfg.reads_option_used = true;
        } else if (arg == "--coverage" || arg == "--depth" || arg == "-c") {
            cfg.coverage_depth = std::stod(require_value(arg));
        } else if (arg == "--seed" || arg == "-s") {
            cfg.seed = std::stoull(require_value(arg));
        } else if (arg == "--threads" || arg == "-j") {
            cfg.thread_count = std::stoull(require_value(arg));
        } else if (arg == "--trans-ratio" || arg == "-t") {
            cfg.trans_ratio = std::stod(require_value(arg));
            cfg.trans_ratio_explicit = true;
        } else if (arg == "--force-contig-reuse") {
            cfg.force_contig_reuse = true;
        } else if (arg == "--long-read-mean") {
            cfg.long_read_mean = std::stoull(require_value(arg));
            cfg.long_read_options_used = true;
            if (cfg.long_read_mean == 0) {
                throw std::runtime_error("--long-read-mean must be positive.");
            }
        } else if (arg == "--long-read-min") {
            cfg.long_read_min = std::stoull(require_value(arg));
            cfg.long_read_options_used = true;
            if (cfg.long_read_min == 0) {
                throw std::runtime_error("--long-read-min must be positive.");
            }
        } else if (arg == "--long-read-max") {
            cfg.long_read_max = std::stoull(require_value(arg));
            cfg.long_read_options_used = true;
            if (cfg.long_read_max == 0) {
                throw std::runtime_error("--long-read-max must be positive.");
            }
        } else if (arg == "--long-read-sigma") {
            cfg.long_read_sigma = std::stod(require_value(arg));
            cfg.long_read_options_used = true;
            cfg.long_read_sigma_explicit = true;
        } else if (arg == "--long-read-qv") {
            cfg.long_read_qv = std::stod(require_value(arg));
            cfg.long_read_options_used = true;
            cfg.long_read_qv_explicit = true;
        } else if (arg == "--max-segments") {
            cfg.max_concatemer_segments = std::stoull(require_value(arg));
            cfg.long_read_options_used = true;
        } else if (arg == "--species-model" || arg == "-S") {
            cfg.species_model = require_value(arg);
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
    if (cfg.pairs_option_used && cfg.reads_option_used) {
        throw std::runtime_error("Use either --pairs or --reads, not both.");
    }
    if (cfg.assay == AssayType::Hic && cfg.reads_option_used) {
        throw std::runtime_error("--reads is only valid with --assay porec or cifi.");
    }
    if (cfg.assay != AssayType::Hic && cfg.pairs_option_used) {
        throw std::runtime_error("--pairs is only valid with --assay hic; use --reads.");
    }
    if (cfg.assay != AssayType::Hic && !cfg.pair_count_explicit) {
        cfg.pair_count = 10000;
    }
    if (cfg.assay != AssayType::Hic && !cfg.enzyme_site_explicit) {
        cfg.enzyme_site = "GATC";
    }
    if (cfg.assay != AssayType::Hic && cfg.read_length_explicit) {
        throw std::runtime_error(
            "--read-length is a Hi-C option; use --long-read-mean for Pore-C/CiFi.");
    }
    if (cfg.read_length == 0) {
        throw std::runtime_error("--read-length must be positive.");
    }
    if (cfg.assay == AssayType::Hic && cfg.read_length != 150) {
        throw std::runtime_error("Only 150 bp paired-end reads are supported.");
    }
    if (!std::isfinite(cfg.coverage_depth) || cfg.coverage_depth < 0.0) {
        throw std::runtime_error("--coverage must be non-negative.");
    }
    if (cfg.coverage_depth > 0.0 && cfg.pair_count_explicit) {
        throw std::runtime_error("Use --coverage or an exact --pairs/--reads count, not both.");
    }
    if (cfg.coverage_depth == 0.0 && cfg.pair_count == 0) {
        throw std::runtime_error("--pairs/--reads must be positive when --coverage is omitted.");
    }
    if (cfg.thread_count == 0) {
        cfg.thread_count = std::max<unsigned int>(1, std::thread::hardware_concurrency());
    }
    if (!std::isfinite(cfg.trans_ratio) || cfg.trans_ratio < 0.0 || cfg.trans_ratio > 1.0) {
        throw std::runtime_error("--trans-ratio must be within [0, 1].");
    }
    if (cfg.assay == AssayType::Hic && cfg.long_read_options_used) {
        throw std::runtime_error("Long-read options require --assay porec or cifi.");
    }
    if (!std::isfinite(cfg.long_read_sigma) || cfg.long_read_sigma < 0.0 ||
        cfg.long_read_sigma > 3.0) {
        throw std::runtime_error("--long-read-sigma must be within [0, 3].");
    }
    if (!std::isfinite(cfg.long_read_qv) || cfg.long_read_qv < 0.0 ||
        cfg.long_read_qv > 41.0) {
        throw std::runtime_error("--long-read-qv must be within (0, 41].");
    }
    if (cfg.long_read_qv_explicit && cfg.long_read_qv <= 0.0) {
        throw std::runtime_error("--long-read-qv must be within (0, 41].");
    }
    if (cfg.max_concatemer_segments < 2 || cfg.max_concatemer_segments > 10000) {
        throw std::runtime_error("--max-segments must be within [2, 10000].");
    }
    constexpr std::size_t maximum_long_read_length = 10000000;
    if (cfg.long_read_mean > maximum_long_read_length ||
        cfg.long_read_min > maximum_long_read_length ||
        cfg.long_read_max > maximum_long_read_length) {
        throw std::runtime_error("Long-read length overrides cannot exceed 10000000 bases.");
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
        if (cfg.assay == AssayType::Hic && cfg.read_length >= reference.total_length()) {
            throw std::runtime_error("Read length must be shorter than the total reference length.");
        }
        if (cfg.assay == AssayType::Hic && cfg.coverage_depth > 0.0) {
            const double requested_pairs =
                std::ceil(cfg.coverage_depth * static_cast<double>(reference.total_length()) /
                          static_cast<double>(2 * cfg.read_length));
            if (!std::isfinite(requested_pairs) ||
                requested_pairs > static_cast<double>(std::numeric_limits<std::size_t>::max())) {
                throw std::runtime_error("Requested Hi-C coverage is too large.");
            }
            cfg.pair_count = std::max<std::size_t>(1, static_cast<std::size_t>(requested_pairs));
        }

        LongReadProfile long_read_profile;
        if (cfg.assay == AssayType::PoreC) {
            long_read_profile = porec_profile(cfg);
        } else if (cfg.assay == AssayType::CiFi) {
            long_read_profile = cifi_profile(cfg);
        }

        const std::size_t total_bins = reference.total_bins(cfg.bin_size);
        if (total_bins == 0) {
            throw std::runtime_error("Reference produced zero bins.");
        }
        std::cerr << "Global bins: " << total_bins << '\n';
        std::cerr << "Run configuration:\n"
                  << "  assay: " << assay_type_to_string(cfg.assay) << '\n'
                  << "  output prefix: " << cfg.output_prefix << '\n'
                  << "  bin size: " << cfg.bin_size << '\n'
                  << "  enzyme site: " << cfg.enzyme_site << '\n'
                  << "  threads: " << cfg.thread_count << '\n'
                  << "  seed: " << cfg.seed << '\n'
                  << "  matrix path: " << (cfg.matrix_path.empty() ? "<none>" : cfg.matrix_path) << '\n'
                  << "  matrix format: " << matrix_format_to_string(cfg.matrix_format) << '\n'
                  << "  offset path: " << (cfg.offset_path.empty() ? "<none>" : cfg.offset_path) << '\n'
                  << "  empirical model: " << (cfg.empirical_model.empty() ? "<auto/none>" : cfg.empirical_model) << '\n'
                  << "  species model: " << cfg.species_model << '\n'
                  << "  model dir: " << cfg.model_dir << '\n'
                  << "  force contig reuse: " << (cfg.force_contig_reuse ? "yes" : "no") << '\n'
                  << "  trans ratio: ";
        if (cfg.trans_ratio_explicit) {
            std::cerr << cfg.trans_ratio << " (explicit)\n";
        } else {
            std::cerr << cfg.trans_ratio << " (default / used only when applicable)\n";
        }
        if (cfg.assay == AssayType::Hic) {
            std::cerr << "  read length: " << cfg.read_length << '\n';
            if (cfg.coverage_depth > 0.0) {
                std::cerr << "  coverage: " << cfg.coverage_depth << "x\n"
                          << "  read pairs: " << cfg.pair_count << " (derived from coverage)\n";
            } else {
                std::cerr << "  coverage: <none>\n"
                          << "  read pairs: " << cfg.pair_count;
                if (cfg.pair_count_explicit) {
                    std::cerr << " (explicit)\n";
                } else {
                    std::cerr << " (default)\n";
                }
            }
        } else {
            std::cerr << "  long-read mean/min/max: " << long_read_profile.mean_length << '/'
                      << long_read_profile.min_length << '/' << long_read_profile.max_length << '\n'
                      << "  long-read lognormal sigma: " << long_read_profile.length_log_sigma << '\n'
                      << "  long-read mean QV: " << long_read_profile.mean_qv << '\n'
                      << "  maximum concatemer segments: " << cfg.max_concatemer_segments << '\n';
            if (cfg.coverage_depth > 0.0) {
                std::cerr << "  requested coverage: " << cfg.coverage_depth
                          << "x (read count derived from mean read length)\n";
            } else {
                std::cerr << "  coverage: <none>\n"
                          << "  long reads: " << cfg.pair_count
                          << (cfg.pair_count_explicit ? " (explicit)\n" : " (default)\n");
            }
        }

        const std::vector<OffsetEntry> reference_offsets = build_reference_offsets(reference, cfg.bin_size);
        std::vector<OffsetEntry> source_offsets;
        if (!cfg.offset_path.empty()) {
            source_offsets = load_offsets(cfg.offset_path);
        }

        ContactMatrix matrix;
        if (!cfg.matrix_path.empty()) {
            std::cerr << "Loading input matrix: " << cfg.matrix_path << '\n';
            const ContactMatrix source_matrix = load_matrix(cfg.matrix_path, 0, cfg.matrix_format);
            matrix = align_matrix_to_reference(source_matrix,
                                              source_offsets,
                                              reference_offsets,
                                              cfg.seed,
                                              "input matrix",
                                              cfg.force_contig_reuse);

            if (!cfg.empirical_model.empty()) {
                const EmpiricalModelSpec empirical_spec =
                    load_empirical_model_spec(cfg.empirical_model, cfg.model_dir);
                std::cerr << "Loading trans replacement model: " << empirical_spec.name << '\n'
                          << "  matrix: " << empirical_spec.matrix_path << '\n'
                          << "  offset: " << empirical_spec.offset_path << '\n';
                const ContactMatrix empirical_source =
                    load_matrix(empirical_spec.matrix_path, 0, empirical_spec.matrix_format);
                const std::vector<OffsetEntry> empirical_offsets =
                    load_offsets(empirical_spec.offset_path);
                const ContactMatrix empirical_matrix =
                    align_matrix_to_reference(empirical_source,
                                              empirical_offsets,
                                              reference_offsets,
                                              cfg.seed + 1,
                                              "empirical trans model",
                                              cfg.force_contig_reuse);

                std::cerr << "Replacing input trans contacts with empirical model trans contacts...\n";
                matrix = replace_trans_contacts(matrix,
                                                empirical_matrix,
                                                reference_offsets,
                                                cfg.trans_ratio,
                                                cfg.trans_ratio_explicit);
            } else if (cfg.trans_ratio_explicit) {
                std::cerr << "Ignoring --trans-ratio because no trans replacement model was provided.\n";
            }
        } else {
            const std::string selected_model = infer_empirical_model(cfg);
            const EmpiricalModelSpec empirical_spec =
                load_empirical_model_spec(selected_model, cfg.model_dir);
            std::cerr << "Loading empirical matrix model: " << empirical_spec.name << '\n'
                      << "  matrix: " << empirical_spec.matrix_path << '\n'
                      << "  offset: " << empirical_spec.offset_path << '\n';
            const ContactMatrix empirical_source =
                load_matrix(empirical_spec.matrix_path, 0, empirical_spec.matrix_format);
            const std::vector<OffsetEntry> empirical_offsets =
                load_offsets(empirical_spec.offset_path);
            matrix = align_matrix_to_reference(empirical_source,
                                               empirical_offsets,
                                               reference_offsets,
                                               cfg.seed,
                                               "empirical matrix model",
                                               cfg.force_contig_reuse);
            std::cerr << "Using empirical model as the full contact matrix.\n";
            if (cfg.trans_ratio_explicit) {
                matrix = apply_trans_ratio(matrix, reference_offsets, cfg.trans_ratio);
            }
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
        if (cfg.assay == AssayType::Hic) {
            std::cerr << "Streaming " << cfg.pair_count << " read pairs to FASTQ with "
                      << cfg.thread_count << " worker thread(s)...\n";
            PairedReadWriter writer(cfg);
            write_hic_reads(cfg, reference, reference_offsets, matrix, writer);

            std::cout << "Assay: hic\n"
                      << "Contigs: " << reference.contigs.size() << '\n'
                      << "Reference length: " << reference.total_length() << '\n'
                      << "Global bins: " << total_bins << '\n'
                      << "Contacts loaded: " << matrix.contacts.size() << '\n'
                      << "Read pairs: " << writer.count() << '\n'
                      << "Coverage: "
                      << (static_cast<double>(writer.count() * 2 * cfg.read_length) /
                          static_cast<double>(reference.total_length())) << "x\n"
                      << "Read length: " << cfg.read_length << '\n'
                      << "Output FASTQ R1: " << writer.read1_path() << '\n'
                      << "Output FASTQ R2: " << writer.read2_path() << '\n';
        } else {
            const LongReadSimulationResult result =
                cfg.assay == AssayType::PoreC
                    ? write_porec_reads(cfg, reference, reference_offsets, matrix)
                    : write_cifi_reads(cfg, reference, reference_offsets, matrix);
            const double actual_coverage =
                static_cast<double>(result.sequenced_bases) /
                static_cast<double>(reference.total_length());
            const double mean_read_length =
                result.read_count == 0
                    ? 0.0
                    : static_cast<double>(result.sequenced_bases) / result.read_count;
            const double mean_segments =
                result.read_count == 0
                    ? 0.0
                    : static_cast<double>(result.concatemer_segments) / result.read_count;
            std::cout << "Assay: " << assay_type_to_string(cfg.assay) << '\n'
                      << "Contigs: " << reference.contigs.size() << '\n'
                      << "Reference length: " << reference.total_length() << '\n'
                      << "Global bins: " << total_bins << '\n'
                      << "Contacts loaded: " << matrix.contacts.size() << '\n'
                      << "Long reads: " << result.read_count << '\n'
                      << "Sequenced bases: " << result.sequenced_bases << '\n'
                      << "Coverage: " << actual_coverage << "x\n"
                      << "Mean read length: " << mean_read_length << '\n'
                      << "Mean concatemer segments: " << mean_segments << '\n'
                      << "Segment-limited reads: " << result.segment_limited_reads << '\n'
                      << "Output FASTQ: " << result.fastq_path << '\n'
                      << "Output truth TSV: " << result.truth_path << '\n';
        }
    } catch (const std::exception &ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        print_usage();
        return 1;
    }

    return 0;
}
