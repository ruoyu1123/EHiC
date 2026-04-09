#pragma once

#include <cstddef>
#include <string>
#include <vector>

struct OffsetEntry {
    std::string contig;
    std::size_t start_bin = 0;
    std::size_t end_bin = 0;
};

struct Contact {
    std::size_t bin1 = 0;
    std::size_t bin2 = 0;
    double weight = 0.0;
};

struct ContactMatrix {
    std::size_t bin_count = 0;
    std::vector<Contact> contacts;
};

std::vector<OffsetEntry> load_offsets(const std::string &path);
ContactMatrix load_matrix(const std::string &path, std::size_t expected_bin_count);
