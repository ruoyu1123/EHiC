#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "matrix.h"

struct Contig {
    std::string name;
    std::string sequence;
};

struct ReferenceGenome {
    std::vector<Contig> contigs;

    std::size_t total_length() const;
    std::size_t total_bins(std::size_t bin_size) const;
};

ReferenceGenome load_reference_fasta(const std::string &path);
std::vector<OffsetEntry> build_reference_offsets(const ReferenceGenome &reference, std::size_t bin_size);
