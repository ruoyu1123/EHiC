#pragma once

#include <cstdint>
#include <string>

struct Config {
    std::string reference_path;
    std::string matrix_path;
    std::string offset_path;
    std::string output_prefix = "sim";
    std::string enzyme_site = "AAGCTT";
    std::size_t bin_size = 0;
    std::size_t read_length = 150;
    std::size_t pair_count = 100000;
    std::size_t thread_count = 1;
    bool pair_count_explicit = false;
    double coverage_depth = 0.0;
    std::uint64_t seed = 1;
    double trans_ratio = 0.10;
    bool trans_ratio_explicit = false;
    std::size_t synthetic_contact_count = 0;
    double cis_decay_alpha = 1.0;
    std::size_t max_cis_distance_bins = 200;
    std::string species_model = "auto";
    std::string arrangement_model = "auto";
    std::string trans_model = "auto";
    std::size_t trans_hotspots = 8;
    double collision_randomness = 0.35;
};
