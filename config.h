#pragma once

#include <cstdint>
#include "matrix.h"
#include <string>

struct Config {
    std::string reference_path;
    std::string matrix_path;
    std::string offset_path;
    std::string output_prefix = "sim";
    std::string enzyme_site = "AAGCTT";
    MatrixFormat matrix_format = MatrixFormat::Auto;
    std::string empirical_model;
    std::string model_dir = "models";
    std::size_t bin_size = 0;
    std::size_t read_length = 150;
    std::size_t pair_count = 100000;
    std::size_t thread_count = 1;
    bool pair_count_explicit = false;
    double coverage_depth = 0.0;
    std::uint64_t seed = 1;
    double trans_ratio = 0.10;
    bool trans_ratio_explicit = false;
    bool force_contig_reuse = false;
    std::string species_model = "auto";
};
