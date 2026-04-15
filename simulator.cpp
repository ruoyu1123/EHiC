#include "simulator.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cctype>

namespace {

void write_wrapped_sequence(std::ofstream &out, const std::string &sequence) {
    constexpr std::size_t line_width = 80;
    for (std::size_t pos = 0; pos < sequence.size(); pos += line_width) {
        out << sequence.substr(pos, line_width) << '\n';
    }
}

char complement(char base) {
    switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

std::string reverse_complement(const std::string &sequence) {
    std::string rc;
    rc.reserve(sequence.size());
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        rc.push_back(complement(*it));
    }
    return rc;
}

std::string quote_arg(const std::string &value) {
    std::string quoted = "\"";
    for (char ch : value) {
        if (ch == '"') {
            quoted += "\\\"";
        } else {
            quoted.push_back(ch);
        }
    }
    quoted.push_back('"');
    return quoted;
}

std::string pad_or_trim_read(const std::string &sequence, std::size_t read_length) {
    if (sequence.size() >= read_length) {
        return sequence.substr(0, read_length);
    }
    return sequence + std::string(read_length - sequence.size(), 'N');
}

void write_fastq_record(std::ofstream &out,
                        const std::string &name,
                        const std::string &sequence) {
    out << '@' << name << '\n'
        << sequence << '\n'
        << "+\n"
        << std::string(sequence.size(), 'I') << '\n';
}

void write_builtin_fastq_pairs(const Config &cfg, const std::vector<LigationProduct> &products) {
    const std::string read1_path = cfg.output_prefix + "1.fq";
    const std::string read2_path = cfg.output_prefix + "2.fq";
    std::ofstream read1_out(read1_path);
    std::ofstream read2_out(read2_path);
    if (!read1_out || !read2_out) {
        throw std::runtime_error("Failed to open FASTQ output files for built-in read simulation.");
    }

    for (std::size_t i = 0; i < products.size(); ++i) {
        const auto &product = products[i];
        const std::string read1 = pad_or_trim_read(product.sequence, cfg.read_length);
        const std::size_t tail_start = product.sequence.size() > cfg.read_length
                                           ? product.sequence.size() - cfg.read_length
                                           : 0;
        const std::string read2_template = product.sequence.substr(tail_start);
        const std::string read2 = pad_or_trim_read(reverse_complement(read2_template), cfg.read_length);
        const std::string name = product.name + "/" + std::to_string(i + 1);
        write_fastq_record(read1_out, name, read1);
        write_fastq_record(read2_out, name, read2);
    }

    std::cout << "Wrote built-in paired FASTQ reads: "
              << read1_path << " and " << read2_path << '\n';
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

    const std::string command = "art_illumina -i " + quote_arg(fasta_path) +
                                " -o " + quote_arg(cfg.output_prefix) +
                                " -l " + std::to_string(cfg.read_length) +
                                " -f " + std::to_string(coverage) +
                                " -m " + std::to_string(static_cast<int>(cfg.insert_mean)) +
                                " -s " + std::to_string(static_cast<int>(cfg.insert_std)) +
                                " -p -na";

    std::cout << "Running: " << command << '\n';
    const int result = std::system(command.c_str());
    if (result != 0) {
        std::cerr << "Warning: art_illumina failed or is unavailable (exit code " << result
                  << "). Falling back to built-in FASTQ generation.\n";
        write_builtin_fastq_pairs(cfg, products);
    }
}
