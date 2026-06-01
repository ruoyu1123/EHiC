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

enum class MatrixFormat {
    Auto,
    Sparse,
    Dense,
    BinarySparse
};

std::vector<OffsetEntry> load_offsets(const std::string &path);
ContactMatrix load_matrix(const std::string &path,
                          std::size_t expected_bin_count,
                          MatrixFormat format = MatrixFormat::Auto);
ContactMatrix remap_matrix_to_reference(const ContactMatrix &source_matrix,
                                        const std::vector<OffsetEntry> &source_offsets,
                                        const std::vector<OffsetEntry> &target_offsets);
ContactMatrix apply_trans_ratio(const ContactMatrix &matrix,
                                const std::vector<OffsetEntry> &offsets,
                                double target_trans_ratio);
ContactMatrix replace_trans_contacts(const ContactMatrix &base,
                                     const ContactMatrix &trans_template,
                                     const std::vector<OffsetEntry> &offsets,
                                     double target_trans_ratio,
                                     bool target_trans_ratio_explicit);
ContactMatrix blend_contact_matrices(const ContactMatrix &primary,
                                     const ContactMatrix &secondary,
                                     double primary_fraction);
