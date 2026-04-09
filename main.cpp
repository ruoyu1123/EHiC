#include "config.h"
#include "fragmenter.h"
#include "matrix.h"
#include "reference.h"
#include "simulator.h"
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void print_usage() {
    std::cerr
        << "Usage:\n"
        << "  hicreate --reference ref.fa --matrix matrix.tsv --bin-size 10000\n"
        << "           --read-length 150 --pairs 100000 --output-prefix sim\n"
        << "           [--offset offsets.tsv] [--enzyme-site AAGCTT] [--skip-art]\n"
        << "           [--seed 123] [--insert-mean 150] [--insert-std 25]\n\n"
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
        } else if (arg == "--skip-art") {
            cfg.skip_art = true;
        } else if (arg == "--help" || arg == "-h") {
            print_usage();
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown argument: " + arg);
        }
    }

    if (cfg.reference_path.empty() || cfg.matrix_path.empty() ||
        cfg.output_prefix.empty() || cfg.bin_size == 0 ||
        cfg.read_length == 0 || cfg.pair_count == 0) {
        throw std::runtime_error("Missing required arguments.");
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

        const ContactMatrix matrix = load_matrix(cfg.matrix_path, total_bins);
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
