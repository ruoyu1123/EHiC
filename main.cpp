#include "config.h"
#include "fragmenter.h"
#include "matrix.h"
#include "reference.h"
#include "simulator.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void print_usage() {
    std::cerr 
        << "Usage:\n"
        << "  hicreate --reference ref.fa [--matrix matrix.tsv] --bin-size 10000\n"
        << "           --read-length 150 --pairs 100000 --output-prefix sim\n"
        << "           [--offset offsets.tsv] [--enzyme-site AAGCTT] [--skip-art]\n"
        << "           [--seed 123] [--insert-mean 150] [--insert-std 25]\n"
        << "           [--trans-ratio 0.10] [--synthetic-contacts 200000]\n"
        << "           [--cis-decay-alpha 1.0] [--max-cis-distance-bins 200]\n"
        << "           [--species-model generic_plant] [--arrangement-model auto]\n"
        << "           [--collision-randomness 0.35]\n\n"
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
        } else if (arg == "--seed") {
            cfg.seed = std::stoull(require_value(arg));
        } else if (arg == "--insert-mean") {
            cfg.insert_mean = std::stod(require_value(arg));
        } else if (arg == "--insert-std") {
            cfg.insert_std = std::stod(require_value(arg));
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
        } else if (arg == "--collision-randomness") {
            cfg.collision_randomness = std::stod(require_value(arg));
        } else if (arg == "--skip-art") {
            cfg.skip_art = true;
        } else if (arg == "--help" || arg == "-h") {
            print_usage();
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown argument: " + arg);
        }
    }

    if (cfg.reference_path.empty() || cfg.output_prefix.empty() || cfg.bin_size == 0 ||
        cfg.read_length == 0 || cfg.pair_count == 0) {
        throw std::runtime_error("Missing required arguments.");
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
        const Config cfg = parse_args(argc, argv);
        const ReferenceGenome reference = load_reference_fasta(cfg.reference_path);
        if (cfg.read_length >= reference.total_length()) {
            throw std::runtime_error("Read length must be shorter than the total reference length.");
        }

        const std::size_t total_bins = reference.total_bins(cfg.bin_size);
        if (total_bins == 0) {
            throw std::runtime_error("Reference produced zero bins.");
        }

        std::vector<OffsetEntry> offsets;
        if (!cfg.offset_path.empty()) {
            offsets = load_offsets(cfg.offset_path);
        } else {
            std::size_t start_bin = 0;
            for (const auto &contig : reference.contigs) {
                const std::size_t contig_bins = (contig.sequence.size() + cfg.bin_size - 1) / cfg.bin_size;
                offsets.push_back(OffsetEntry{contig.name, start_bin, start_bin + contig_bins});
                start_bin += contig_bins;
            }
        }

        ContactMatrix matrix;
        if (!cfg.matrix_path.empty()) {
            matrix = load_matrix(cfg.matrix_path, total_bins);
            matrix = apply_trans_ratio(matrix, offsets, cfg.trans_ratio, cfg.seed);
        } else {
            std::size_t synthetic_contacts = cfg.synthetic_contact_count;
            if (synthetic_contacts == 0) {
                synthetic_contacts = std::max<std::size_t>(
                    20000, std::min<std::size_t>(2000000, cfg.pair_count * 20));
            }
            SyntheticModelOptions model_options;
            model_options.trans_ratio = cfg.trans_ratio;
            model_options.cis_decay_alpha = cfg.cis_decay_alpha;
            model_options.max_cis_distance_bins = cfg.max_cis_distance_bins;
            model_options.species_model = cfg.species_model;
            model_options.arrangement_model = cfg.arrangement_model;
            model_options.collision_randomness = cfg.collision_randomness;
            model_options.seed = cfg.seed;
            matrix = generate_synthetic_matrix(total_bins,
                                              offsets,
                                              synthetic_contacts,
                                              model_options);
        }
        const auto ligation_products = create_ligation_products(cfg, reference, offsets, matrix);
        simulate_reads_with_artillumina(cfg, ligation_products);

        std::size_t total_library_bases = 0;
        for (const auto &product : ligation_products) {
            total_library_bases += product.sequence.size();
        }

        std::cout << "Contigs: " << reference.contigs.size() << '\n'
                  << "Reference length: " << reference.total_length() << '\n'
                  << "Global bins: " << total_bins << '\n'
                  << "Contacts loaded: " << matrix.contacts.size() << '\n'
                  << "Ligation products: " << ligation_products.size() << '\n'
                  << "Library bases: " << total_library_bases << '\n'
                  << "Output FASTA: " << cfg.output_prefix << "_fragments.fa\n";
    } catch (const std::exception &ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        print_usage();
        return 1;
    }

    return 0;
}
