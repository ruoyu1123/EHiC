#pragma once

#include <cstddef>
#include <cstdint>
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

struct SyntheticModelOptions {
    double trans_ratio = 0.10;
    double cis_decay_alpha = 1.0;
    std::size_t max_cis_distance_bins = 200;
    std::string species_model = "generic_plant";
    std::string arrangement_model = "auto";
    double collision_randomness = 0.35;
    std::uint64_t seed = 1;
};

std::vector<OffsetEntry> load_offsets(const std::string &path);
ContactMatrix load_matrix(const std::string &path, std::size_t expected_bin_count);
ContactMatrix generate_synthetic_matrix(std::size_t bin_count,
                                        const std::vector<OffsetEntry> &offsets,
                                        std::size_t synthetic_contact_count,
                                        const SyntheticModelOptions &options);
ContactMatrix apply_trans_ratio(const ContactMatrix &matrix,
                                const std::vector<OffsetEntry> &offsets,
                                double target_trans_ratio,
                                std::uint64_t seed);
