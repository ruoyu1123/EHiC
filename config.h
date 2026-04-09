#pragma once

#include <cstdint>
#include <string>

struct Config {
    std::string reference_path;
    std::string matrix_path;
    std::string offset_path;
    std::string output_prefix;
    std::string enzyme_site = "AAGCTT";
    std::size_t bin_size = 0;
    std::size_t read_length = 0;
    std::size_t pair_count = 0;
    std::uint64_t seed = 1;
    double insert_mean = 150.0;
    double insert_std = 25.0;
    bool skip_art = false;
};
