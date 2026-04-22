#include "simulator.h"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

namespace {

char mutate_base(char base, std::mt19937_64 &rng) {
    const char upper = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    if (upper != 'A' && upper != 'C' && upper != 'G' && upper != 'T') {
        return 'N';
    }

    char alternatives[3];
    int index = 0;
    for (char candidate : {'A', 'C', 'G', 'T'}) {
        if (candidate != upper) {
            alternatives[index++] = candidate;
        }
    }
    std::uniform_int_distribution<int> dist(0, 2);
    return alternatives[dist(rng)];
}

int illumina_like_quality(std::size_t position, std::size_t read_length, std::mt19937_64 &rng) {
    const double x = read_length <= 1
                         ? 0.0
                         : static_cast<double>(position) / static_cast<double>(read_length - 1);
    const double mean_q = 36.5 - 7.5 * std::pow(x, 1.7);
    const double stddev_q = 1.2 + 2.0 * x;
    std::normal_distribution<double> q_dist(mean_q, stddev_q);
    const int q = static_cast<int>(std::lround(q_dist(rng)));
    return std::max(2, std::min(41, q));
}

std::string normalize_template(const std::string &sequence, std::size_t read_length) {
    std::string normalized;
    normalized.reserve(read_length);
    for (std::size_t i = 0; i < read_length; ++i) {
        const char base = i < sequence.size()
                              ? static_cast<char>(std::toupper(static_cast<unsigned char>(sequence[i])))
                              : 'N';
        switch (base) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'N':
                normalized.push_back(base);
                break;
            default:
                normalized.push_back('N');
                break;
        }
    }
    return normalized;
}

struct SimulatedRead {
    std::string sequence;
    std::string qualities;
};

SimulatedRead simulate_read(const std::string &read_template,
                            std::size_t read_length,
                            std::mt19937_64 &rng) {
    std::string sequence = normalize_template(read_template, read_length);
    std::string qualities;
    qualities.reserve(read_length);

    std::uniform_real_distribution<double> unit(0.0, 1.0);
    for (std::size_t i = 0; i < read_length; ++i) {
        const int q = illumina_like_quality(i, read_length, rng);
        qualities.push_back(static_cast<char>(q + 33));

        if (sequence[i] == 'N') {
            continue;
        }

        const double error_probability = std::pow(10.0, -static_cast<double>(q) / 10.0);
        if (unit(rng) < error_probability) {
            sequence[i] = mutate_base(sequence[i], rng);
        }
    }

    return SimulatedRead{sequence, qualities};
}

void write_fastq_record(std::ofstream &out,
                        const std::string &name,
                        const SimulatedRead &read) {
    out << '@' << name << '\n'
        << read.sequence << '\n'
        << "+\n"
        << read.qualities << '\n';
}

}  // namespace

void simulate_paired_reads(const Config &cfg, const std::vector<ReadPairTemplate> &read_templates) {
    if (read_templates.empty()) {
        throw std::runtime_error("No read templates were generated.");
    }
    if (cfg.read_length != 150) {
        throw std::runtime_error("Built-in read simulation currently supports only 150 bp paired-end reads.");
    }

    const std::string read1_path = cfg.output_prefix + "_R1.fastq";
    const std::string read2_path = cfg.output_prefix + "_R2.fastq";
    std::ofstream read1_out(read1_path);
    std::ofstream read2_out(read2_path);
    if (!read1_out || !read2_out) {
        throw std::runtime_error("Failed to open paired FASTQ output files.");
    }

    std::mt19937_64 rng(cfg.seed + 2027);
    for (std::size_t i = 0; i < read_templates.size(); ++i) {
        const auto &templ = read_templates[i];
        const SimulatedRead read1 = simulate_read(templ.read1, cfg.read_length, rng);
        const SimulatedRead read2 = simulate_read(templ.read2, cfg.read_length, rng);
        write_fastq_record(read1_out, templ.name + "/1", read1);
        write_fastq_record(read2_out, templ.name + "/2", read2);
    }

    std::cout << "Output read pairs: " << read_templates.size() << '\n'
              << "Output FASTQ R1: " << read1_path << '\n'
              << "Output FASTQ R2: " << read2_path << '\n';
}
