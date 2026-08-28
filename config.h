#pragma once

#include <cstdint>
#include "matrix.h"
#include <string>

enum class AssayType {
    Hic,
    PoreC,
    CiFi
};

struct Config {
    AssayType assay = AssayType::Hic;
    std::string reference_path;
    std::string matrix_path;
    std::string offset_path;
    std::string output_prefix = "sim";
    std::string enzyme_site = "AAGCTT";
    bool enzyme_site_explicit = false;
    MatrixFormat matrix_format = MatrixFormat::Auto;
    std::string empirical_model;
    std::string model_dir = "models";
    std::size_t bin_size = 0;
    std::size_t read_length = 150;
    bool read_length_explicit = false;
    std::size_t pair_count = 100000;
    std::size_t thread_count = 1;
    bool pair_count_explicit = false;
    bool pairs_option_used = false;
    bool reads_option_used = false;
    double coverage_depth = 0.0;
    std::uint64_t seed = 1;
    double trans_ratio = 0.10;
    bool trans_ratio_explicit = false;
    bool force_contig_reuse = false;
    std::string species_model = "auto";
    std::size_t long_read_mean = 0;
    std::size_t long_read_min = 0;
    std::size_t long_read_max = 0;
    double long_read_sigma = 0.0;
    double long_read_qv = 0.0;
    std::size_t max_concatemer_segments = 256;
    bool long_read_options_used = false;
    bool long_read_sigma_explicit = false;
    bool long_read_qv_explicit = false;
};
