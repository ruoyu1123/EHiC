#include "simulator.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void write_wrapped_sequence(std::ofstream &out, const std::string &sequence) {
    constexpr std::size_t line_width = 80;
    for (std::size_t pos = 0; pos < sequence.size(); pos += line_width) {
        out << sequence.substr(pos, line_width) << '\n';
    }
}

}  // namespace

void simulate_reads_with_artillumina(const Config &cfg, const std::vector<LigationProduct> &products) {
    const std::string fasta_path = cfg.output_prefix + "_fragments.fa";
    std::ofstream fasta_out(fasta_path);
    if (!fasta_out) {
        throw std::runtime_error("Failed to open FASTA output file: " + fasta_path);
    }

    std::size_t total_bases = 0;
    for (const auto &product : products) {
        fasta_out << '>' << product.name << '\n';
        write_wrapped_sequence(fasta_out, product.sequence);
        total_bases += product.sequence.size();
    }
    fasta_out.close();

    if (cfg.skip_art) {
        std::cout << "Skipping ART simulation because --skip-art was requested.\n";
        return;
    }

    if (total_bases == 0) {
        throw std::runtime_error("No ligation products were generated.");
    }

    const double coverage =
        static_cast<double>(cfg.pair_count * 2 * cfg.read_length) / static_cast<double>(total_bases);

    const std::string command = "art_illumina -i " + fasta_path +
                                " -o " + cfg.output_prefix +
                                " -l " + std::to_string(cfg.read_length) +
                                " -f " + std::to_string(coverage) +
                                " -m " + std::to_string(static_cast<int>(cfg.insert_mean)) +
                                " -s " + std::to_string(static_cast<int>(cfg.insert_std)) +
                                " -p -na";

    std::cout << "Running: " << command << '\n';
    const int result = std::system(command.c_str());
    if (result != 0) {
        std::cerr << "Warning: art_illumina failed or is unavailable (exit code " << result
                  << "). FASTA library was still written to " << fasta_path << ".\n";
    }
}
