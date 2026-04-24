#pragma once

#include "config.h"
#include <cstddef>
#include <fstream>
#include <random>
#include <string>

#include "fragmenter.h"

struct FastqBlock {
    std::size_t pair_count = 0;
    std::string read1;
    std::string read2;
};

void append_simulated_fastq_pair(const Config &cfg,
                                 const ReadPairTemplate &read_template,
                                 std::mt19937_64 &rng,
                                 std::string &read1_out,
                                 std::string &read2_out);

class PairedReadWriter {
public:
    explicit PairedReadWriter(const Config &cfg);
    void write_template(const ReadPairTemplate &read_template);
    void write_block(const FastqBlock &block);
    std::size_t count() const;
    const std::string &read1_path() const;
    const std::string &read2_path() const;

private:
    const Config &cfg_;
    std::string read1_path_;
    std::string read2_path_;
    std::ofstream read1_out_;
    std::ofstream read2_out_;
    std::mt19937_64 rng_;
    std::size_t count_ = 0;
};
