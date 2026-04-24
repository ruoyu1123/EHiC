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

void append_fastq_record(std::string &out,
                         const std::string &name,
                         const SimulatedRead &read) {
    out.push_back('@');
    out += name;
    out.push_back('\n');
    out += read.sequence;
    out += "\n+\n";
    out += read.qualities;
    out.push_back('\n');
}

}  // namespace

void append_simulated_fastq_pair(const Config &cfg,
                                 const ReadPairTemplate &read_template,
                                 std::mt19937_64 &rng,
                                 std::string &read1_out,
                                 std::string &read2_out) {
    const SimulatedRead read1 = simulate_read(read_template.read1, cfg.read_length, rng);
    const SimulatedRead read2 = simulate_read(read_template.read2, cfg.read_length, rng);
    append_fastq_record(read1_out, read_template.name + "/1", read1);
    append_fastq_record(read2_out, read_template.name + "/2", read2);
}

PairedReadWriter::PairedReadWriter(const Config &cfg)
    : cfg_(cfg),
      read1_path_(cfg.output_prefix + "_R1.fastq"),
      read2_path_(cfg.output_prefix + "_R2.fastq"),
      read1_out_(read1_path_),
      read2_out_(read2_path_),
      rng_(cfg.seed + 2027) {
    if (cfg_.read_length != 150) {
        throw std::runtime_error("Built-in read simulation currently supports only 150 bp paired-end reads.");
    }
    if (!read1_out_ || !read2_out_) {
        throw std::runtime_error("Failed to open paired FASTQ output files.");
    }
}

void PairedReadWriter::write_template(const ReadPairTemplate &read_template) {
    const SimulatedRead read1 = simulate_read(read_template.read1, cfg_.read_length, rng_);
    const SimulatedRead read2 = simulate_read(read_template.read2, cfg_.read_length, rng_);
    write_fastq_record(read1_out_, read_template.name + "/1", read1);
    write_fastq_record(read2_out_, read_template.name + "/2", read2);
    const std::size_t previous_count = count_;
    ++count_;
    if (count_ / 1000000 > previous_count / 1000000) {
        std::cerr << "Wrote " << count_ << " read pairs...\n";
    }
}

void PairedReadWriter::write_block(const FastqBlock &block) {
    read1_out_ << block.read1;
    read2_out_ << block.read2;
    const std::size_t previous_count = count_;
    count_ += block.pair_count;
    if (count_ / 1000000 > previous_count / 1000000) {
        std::cerr << "Wrote " << count_ << " read pairs...\n";
    }
}

std::size_t PairedReadWriter::count() const {
    return count_;
}

const std::string &PairedReadWriter::read1_path() const {
    return read1_path_;
}

const std::string &PairedReadWriter::read2_path() const {
    return read2_path_;
}
