#include "matrix.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace {

struct ContactKey {
    std::size_t bin1 = 0;
    std::size_t bin2 = 0;

    bool operator==(const ContactKey &other) const {
        return bin1 == other.bin1 && bin2 == other.bin2;
    }
};

struct ContactKeyHash {
    std::size_t operator()(const ContactKey &key) const {
        const std::size_t left = std::hash<std::size_t>{}(key.bin1);
        const std::size_t right = std::hash<std::size_t>{}(key.bin2);
        return left ^ (right + 0x9e3779b97f4a7c15ULL + (left << 6U) + (left >> 2U));
    }
};

using ContactWeightMap = std::unordered_map<ContactKey, double, ContactKeyHash>;

std::string to_lower(std::string s);

std::string trim(const std::string &s) {
    const auto begin = s.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) {
        return "";
    }
    const auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(begin, end - begin + 1);
}

std::string normalize_delimiters(std::string line) {
    std::string normalized;
    normalized.reserve(line.size());
    for (std::size_t i = 0; i < line.size(); ++i) {
        if (line[i] == '`' && i + 1 < line.size() && line[i + 1] == 't') {
            normalized.push_back('\t');
            ++i;
        } else if (line[i] == ',') {
            normalized.push_back(' ');
        } else {
            normalized.push_back(line[i]);
        }
    }
    return normalized;
}

std::vector<std::string> split_tokens(const std::string &line) {
    std::istringstream stream(normalize_delimiters(line));
    std::vector<std::string> tokens;
    std::string token;
    while (stream >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

std::string token_lower(const std::string &token) {
    return to_lower(trim(token));
}

bool is_sparse_header(const std::vector<std::string> &tokens) {
    if (tokens.size() < 3) {
        return false;
    }
    return token_lower(tokens[0]) == "bin1" &&
           token_lower(tokens[1]) == "bin2" &&
           (token_lower(tokens[2]) == "value" ||
            token_lower(tokens[2]) == "weight" ||
            token_lower(tokens[2]) == "count");
}

bool is_offset_header(const std::vector<std::string> &tokens) {
    if (tokens.size() < 3) {
        return false;
    }
    return token_lower(tokens[0]) == "contig" &&
           token_lower(tokens[1]) == "start_bin" &&
           token_lower(tokens[2]) == "end_bin";
}

bool has_binary_matrix_suffix(const std::string &path) {
    const std::string lower = to_lower(path);
    return lower.size() >= 4 &&
           (lower.substr(lower.size() - 4) == ".hcm" ||
            lower.substr(lower.size() - 5) == ".hcmx");
}

std::uint64_t read_u64(std::ifstream &in, const std::string &label) {
    std::uint64_t value = 0;
    in.read(reinterpret_cast<char *>(&value), sizeof(value));
    if (!in) {
        throw std::runtime_error("Failed to read binary matrix " + label + ".");
    }
    return value;
}

double read_f64(std::ifstream &in, const std::string &label) {
    double value = 0.0;
    in.read(reinterpret_cast<char *>(&value), sizeof(value));
    if (!in) {
        throw std::runtime_error("Failed to read binary matrix " + label + ".");
    }
    return value;
}

ContactMatrix load_binary_sparse_matrix(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open binary matrix file: " + path);
    }

    char magic[8] = {};
    in.read(magic, sizeof(magic));
    const char expected[8] = {'H', 'C', 'M', 'O', 'D', '1', '\0', '\0'};
    if (!in || !std::equal(std::begin(magic), std::end(magic), std::begin(expected))) {
        throw std::runtime_error("Binary matrix has invalid magic header; expected HCMOD1.");
    }

    const std::uint64_t raw_bin_count = read_u64(in, "bin_count");
    const std::uint64_t raw_contact_count = read_u64(in, "contact_count");
    if (raw_bin_count > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()) ||
        raw_contact_count > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::runtime_error("Binary matrix is too large for this platform.");
    }

    ContactMatrix matrix;
    matrix.bin_count = static_cast<std::size_t>(raw_bin_count);
    matrix.contacts.reserve(static_cast<std::size_t>(raw_contact_count));
    for (std::uint64_t i = 0; i < raw_contact_count; ++i) {
        const std::uint64_t raw_bin1 = read_u64(in, "bin1");
        const std::uint64_t raw_bin2 = read_u64(in, "bin2");
        const double weight = read_f64(in, "weight");
        if (raw_bin1 >= raw_bin_count || raw_bin2 >= raw_bin_count) {
            throw std::runtime_error("Binary matrix contact bin exceeds bin_count.");
        }
        if (std::isfinite(weight) && weight > 0.0) {
            matrix.contacts.push_back(Contact{
                static_cast<std::size_t>(raw_bin1),
                static_cast<std::size_t>(raw_bin2),
                weight
            });
        }
    }
    if (matrix.contacts.empty()) {
        throw std::runtime_error("Binary sparse matrix does not contain any positive contacts.");
    }
    return matrix;
}

ContactMatrix dense_to_sparse(const std::vector<std::vector<double>> &rows,
                              std::size_t expected_bin_count) {
    if (rows.empty() || rows.front().empty()) {
        throw std::runtime_error("Matrix file is empty.");
    }

    const std::size_t row_count = rows.size();
    const std::size_t col_count = rows.front().size();
    if (row_count != col_count) {
        throw std::runtime_error("Dense matrix must be square for Hi-C contacts.");
    }

    const std::size_t final_size = expected_bin_count == 0 ? row_count : expected_bin_count;
    std::vector<double> resized(final_size * final_size, 0.0);

    auto at = [&](std::size_t r, std::size_t c) -> double & {
        return resized[r * final_size + c];
    };

    auto sample = [&](double row_pos, double col_pos) {
        const auto clamp_index = [](double x, std::size_t max_index) {
            return std::max(0.0, std::min(x, static_cast<double>(max_index)));
        };

        row_pos = clamp_index(row_pos, row_count - 1);
        col_pos = clamp_index(col_pos, col_count - 1);

        const std::size_t r0 = static_cast<std::size_t>(std::floor(row_pos));
        const std::size_t c0 = static_cast<std::size_t>(std::floor(col_pos));
        const std::size_t r1 = std::min(r0 + 1, row_count - 1);
        const std::size_t c1 = std::min(c0 + 1, col_count - 1);
        const double fr = row_pos - static_cast<double>(r0);
        const double fc = col_pos - static_cast<double>(c0);

        const double v00 = rows[r0][c0];
        const double v01 = rows[r0][c1];
        const double v10 = rows[r1][c0];
        const double v11 = rows[r1][c1];

        const double top = v00 * (1.0 - fc) + v01 * fc;
        const double bottom = v10 * (1.0 - fc) + v11 * fc;
        return top * (1.0 - fr) + bottom * fr;
    };

    if (final_size == 1) {
        at(0, 0) = std::max(0.0, rows[0][0]);
    } else {
        const double row_scale =
            row_count > 1 ? static_cast<double>(row_count - 1) / static_cast<double>(final_size - 1) : 0.0;
        const double col_scale =
            col_count > 1 ? static_cast<double>(col_count - 1) / static_cast<double>(final_size - 1) : 0.0;

        for (std::size_t r = 0; r < final_size; ++r) {
            for (std::size_t c = 0; c < final_size; ++c) {
                at(r, c) = std::max(0.0, sample(static_cast<double>(r) * row_scale,
                                                static_cast<double>(c) * col_scale));
            }
        }
    }

    ContactMatrix matrix;
    matrix.bin_count = final_size;
    for (std::size_t r = 0; r < final_size; ++r) {
        for (std::size_t c = r; c < final_size; ++c) {
            const double weight = 0.5 * (at(r, c) + at(c, r));
            if (weight > 0.0) {
                matrix.contacts.push_back(Contact{r, c, weight});
            }
        }
    }
    return matrix;
}

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

std::vector<int> build_bin_to_contig(std::size_t bin_count, const std::vector<OffsetEntry> &offsets) {
    std::vector<int> bin_to_contig(bin_count, -1);
    for (std::size_t i = 0; i < offsets.size(); ++i) {
        const auto &offset = offsets[i];
        if (offset.start_bin >= offset.end_bin || offset.end_bin > bin_count) {
            throw std::runtime_error("Offset bins are invalid for matrix bin count.");
        }
        for (std::size_t bin = offset.start_bin; bin < offset.end_bin; ++bin) {
            bin_to_contig[bin] = static_cast<int>(i);
        }
    }
    return bin_to_contig;
}

std::size_t total_offset_bins(const std::vector<OffsetEntry> &offsets) {
    std::size_t total = 0;
    for (const auto &offset : offsets) {
        if (offset.end_bin > offset.start_bin) {
            total += offset.end_bin - offset.start_bin;
        }
    }
    return total;
}

double sum_contact_weights(const ContactMatrix &matrix) {
    double total = 0.0;
    for (const auto &contact : matrix.contacts) {
        if (std::isfinite(contact.weight) && contact.weight > 0.0) {
            total += contact.weight;
        }
    }
    return total;
}

std::vector<int> build_contig_lookup_by_bin(std::size_t bin_count, const std::vector<OffsetEntry> &offsets) {
    std::vector<int> bin_to_contig(bin_count, -1);
    for (std::size_t i = 0; i < offsets.size(); ++i) {
        const auto &offset = offsets[i];
        if (offset.start_bin >= offset.end_bin || offset.start_bin >= bin_count) {
            continue;
        }
        const std::size_t end_bin = std::min(offset.end_bin, bin_count);
        for (std::size_t bin = offset.start_bin; bin < end_bin; ++bin) {
            if (bin_to_contig[bin] == -1) {
                bin_to_contig[bin] = static_cast<int>(i);
            }
        }
    }
    return bin_to_contig;
}

bool contact_is_cis(const Contact &contact, const std::vector<int> &bin_to_contig) {
    if (contact.bin1 >= bin_to_contig.size() || contact.bin2 >= bin_to_contig.size()) {
        return false;
    }
    const int contig1 = bin_to_contig[contact.bin1];
    const int contig2 = bin_to_contig[contact.bin2];
    return contig1 >= 0 && contig1 == contig2;
}

std::size_t map_bin_by_scale(std::size_t source_bin, std::size_t source_bin_count, std::size_t target_bin_count) {
    if (target_bin_count == 0) {
        throw std::runtime_error("Target bin count must be positive.");
    }
    if (target_bin_count == 1 || source_bin_count <= 1) {
        return 0;
    }

    const double scaled =
        (static_cast<double>(source_bin) + 0.5) * static_cast<double>(target_bin_count) /
        static_cast<double>(source_bin_count);
    const std::size_t mapped = static_cast<std::size_t>(std::floor(scaled));
    return std::min(target_bin_count - 1, mapped);
}

std::size_t map_bin_within_contig(std::size_t source_bin,
                                  const OffsetEntry &source_offset,
                                  const OffsetEntry &target_offset) {
    const std::size_t source_span = source_offset.end_bin - source_offset.start_bin;
    const std::size_t target_span = target_offset.end_bin - target_offset.start_bin;
    if (target_span == 0) {
        throw std::runtime_error("Target offset span must be positive.");
    }
    if (target_span == 1 || source_span <= 1) {
        return target_offset.start_bin;
    }

    const std::size_t local_source = source_bin - source_offset.start_bin;
    const double scaled =
        (static_cast<double>(local_source) + 0.5) * static_cast<double>(target_span) /
        static_cast<double>(source_span);
    const std::size_t local_target = static_cast<std::size_t>(std::floor(scaled));
    return target_offset.start_bin + std::min(target_span - 1, local_target);
}

ContactKey contact_key(std::size_t bin1, std::size_t bin2) {
    const std::size_t lo = std::min(bin1, bin2);
    const std::size_t hi = std::max(bin1, bin2);
    return ContactKey{lo, hi};
}

void add_contact_weight(ContactWeightMap &weights,
                        std::size_t bin1,
                        std::size_t bin2,
                        double value) {
    if (!std::isfinite(value) || value <= 0.0) {
        return;
    }
    weights[contact_key(bin1, bin2)] += value;
}

ContactMatrix finalize_sparse_contacts(std::size_t bin_count,
                                       const ContactWeightMap &weights) {
    ContactMatrix matrix;
    matrix.bin_count = bin_count;
    matrix.contacts.reserve(weights.size());
    for (const auto &[key, weight] : weights) {
        if (!std::isfinite(weight) || weight <= 0.0) {
            continue;
        }
        const std::size_t bin1 = key.bin1;
        const std::size_t bin2 = key.bin2;
        if (bin1 >= bin_count || bin2 >= bin_count) {
            continue;
        }
        matrix.contacts.push_back(Contact{bin1, bin2, weight});
    }

    if (matrix.contacts.empty()) {
        throw std::runtime_error("No positive contacts remain in matrix.");
    }

    return matrix;
}

}  // namespace

std::vector<OffsetEntry> load_offsets(const std::string &path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open offset file: " + path);
    }

    std::vector<OffsetEntry> offsets;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') {
            continue;
        }

        const auto tokens = split_tokens(line);
        if (tokens.size() < 3) {
            throw std::runtime_error("Offset rows must contain: contig start_bin end_bin");
        }
        if (is_offset_header(tokens)) {
            continue;
        }

        OffsetEntry entry;
        entry.contig = tokens[0];
        entry.start_bin = std::stoull(tokens[1]);
        entry.end_bin = std::stoull(tokens[2]);
        if (entry.end_bin <= entry.start_bin) {
            throw std::runtime_error("Offset end_bin must be greater than start_bin.");
        }
        offsets.push_back(entry);
    }

    if (offsets.empty()) {
        throw std::runtime_error("Offset file is empty.");
    }

    std::sort(offsets.begin(), offsets.end(),
              [](const OffsetEntry &lhs, const OffsetEntry &rhs) { return lhs.start_bin < rhs.start_bin; });
    for (std::size_t i = 1; i < offsets.size(); ++i) {
        if (offsets[i].start_bin < offsets[i - 1].end_bin) {
            throw std::runtime_error("Offset intervals overlap.");
        }
    }

    return offsets;
}

ContactMatrix load_matrix(const std::string &path,
                          std::size_t expected_bin_count,
                          MatrixFormat format) {
    if (format == MatrixFormat::BinarySparse ||
        (format == MatrixFormat::Auto && has_binary_matrix_suffix(path))) {
        ContactMatrix matrix = load_binary_sparse_matrix(path);
        if (expected_bin_count != 0 && matrix.bin_count != expected_bin_count) {
            throw std::runtime_error("Binary sparse matrix bin_count does not match expected_bin_count.");
        }
        return matrix;
    }

    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open matrix file: " + path);
    }

    std::string line;
    std::vector<std::vector<double>> dense_rows;
    ContactWeightMap sparse_weights;
    std::size_t detected_columns = 0;
    bool sparse_mode = format == MatrixFormat::Sparse;
    bool mode_decided = format != MatrixFormat::Auto;
    std::size_t max_bin = 0;

    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') {
            continue;
        }

        const auto tokens = split_tokens(line);
        if (tokens.empty()) {
            continue;
        }
        if (is_sparse_header(tokens)) {
            if (format == MatrixFormat::Dense) {
                throw std::runtime_error("Sparse matrix header found while --matrix-format dense was requested.");
            }
            sparse_mode = true;
            mode_decided = true;
            continue;
        }

        if (!mode_decided) {
            sparse_mode = tokens.size() == 3;
            mode_decided = true;
        }

        if (sparse_mode) {
            if (tokens.size() != 3) {
                throw std::runtime_error("Sparse matrix rows must contain exactly 3 columns: bin1 bin2 value");
            }

            const std::size_t bin1 = std::stoull(tokens[0]);
            const std::size_t bin2 = std::stoull(tokens[1]);
            const double weight = std::stod(tokens[2]);
            if (!std::isfinite(weight) || weight <= 0.0) {
                continue;
            }

            add_contact_weight(sparse_weights, bin1, bin2, weight);
            max_bin = std::max(max_bin, std::max(bin1, bin2));
        } else {
            if (detected_columns == 0) {
                detected_columns = tokens.size();
            } else if (tokens.size() != detected_columns) {
                throw std::runtime_error("Dense matrix rows have inconsistent column counts.");
            }

            std::vector<double> row;
            row.reserve(tokens.size());
            for (const auto &token : tokens) {
                const double value = std::stod(token);
                if (!std::isfinite(value)) {
                    throw std::runtime_error("Matrix contains non-finite values.");
                }
                row.push_back(value);
            }
            dense_rows.push_back(std::move(row));
        }
    }

    if (!mode_decided) {
        throw std::runtime_error("Matrix file is empty.");
    }

    if (!sparse_mode) {
        return dense_to_sparse(dense_rows, expected_bin_count);
    }

    ContactMatrix matrix;
    matrix.bin_count = max_bin + 1;
    for (const auto &[key, weight] : sparse_weights) {
        matrix.contacts.push_back(Contact{key.bin1, key.bin2, weight});
    }

    if (matrix.contacts.empty()) {
        throw std::runtime_error("Sparse matrix does not contain any positive contacts.");
    }

    return matrix;
}

ContactMatrix remap_matrix_to_reference(const ContactMatrix &source_matrix,
                                        const std::vector<OffsetEntry> &source_offsets,
                                        const std::vector<OffsetEntry> &target_offsets) {
    const std::size_t target_bin_count = total_offset_bins(target_offsets);
    if (target_bin_count == 0) {
        throw std::runtime_error("Reference offsets are empty.");
    }
    if (source_matrix.bin_count == 0 || source_matrix.contacts.empty()) {
        return ContactMatrix{target_bin_count, {}};
    }

    const std::vector<int> source_bin_to_contig =
        build_contig_lookup_by_bin(source_matrix.bin_count, source_offsets);
    std::unordered_map<std::string, OffsetEntry> target_by_name;
    target_by_name.reserve(target_offsets.size());
    for (const auto &offset : target_offsets) {
        target_by_name[offset.contig] = offset;
    }

    ContactWeightMap remapped_weights;
    remapped_weights.reserve(source_matrix.contacts.size() * 2);

    for (const auto &contact : source_matrix.contacts) {
        if (contact.weight <= 0.0 || !std::isfinite(contact.weight)) {
            continue;
        }
        if (contact.bin1 >= source_matrix.bin_count || contact.bin2 >= source_matrix.bin_count) {
            continue;
        }

        std::size_t mapped1 = 0;
        std::size_t mapped2 = 0;
        bool mapped_ok = false;

        const int source_contig1 = contact.bin1 < source_bin_to_contig.size() ? source_bin_to_contig[contact.bin1] : -1;
        const int source_contig2 = contact.bin2 < source_bin_to_contig.size() ? source_bin_to_contig[contact.bin2] : -1;

        if (!source_offsets.empty() && source_contig1 >= 0 && source_contig2 >= 0) {
            const OffsetEntry &src_offset1 = source_offsets[static_cast<std::size_t>(source_contig1)];
            const OffsetEntry &src_offset2 = source_offsets[static_cast<std::size_t>(source_contig2)];
            const auto target_it1 = target_by_name.find(src_offset1.contig);
            const auto target_it2 = target_by_name.find(src_offset2.contig);
            if (target_it1 != target_by_name.end() && target_it2 != target_by_name.end()) {
                mapped1 = map_bin_within_contig(contact.bin1, src_offset1, target_it1->second);
                mapped2 = map_bin_within_contig(contact.bin2, src_offset2, target_it2->second);
                mapped_ok = true;
            }
        }

        if (!mapped_ok) {
            mapped1 = map_bin_by_scale(contact.bin1, source_matrix.bin_count, target_bin_count);
            mapped2 = map_bin_by_scale(contact.bin2, source_matrix.bin_count, target_bin_count);
        }

        add_contact_weight(remapped_weights, mapped1, mapped2, contact.weight);
    }

    if (remapped_weights.empty()) {
        return ContactMatrix{target_bin_count, {}};
    }

    return finalize_sparse_contacts(target_bin_count, remapped_weights);
}

ContactMatrix blend_contact_matrices(const ContactMatrix &primary,
                                     const ContactMatrix &secondary,
                                     double primary_fraction) {
    primary_fraction = std::max(0.0, std::min(1.0, primary_fraction));
    const std::size_t bin_count = std::max(primary.bin_count, secondary.bin_count);
    if (bin_count == 0) {
        return ContactMatrix{};
    }

    const double primary_total = sum_contact_weights(primary);
    const double secondary_total = sum_contact_weights(secondary);
    ContactWeightMap combined;
    combined.reserve(primary.contacts.size() + secondary.contacts.size());

    if (primary_total > 0.0 && primary_fraction > 0.0) {
        const double scale = primary_fraction / primary_total;
        for (const auto &contact : primary.contacts) {
            if (contact.bin1 < bin_count && contact.bin2 < bin_count) {
                add_contact_weight(combined, contact.bin1, contact.bin2, contact.weight * scale);
            }
        }
    }

    if (secondary_total > 0.0 && primary_fraction < 1.0) {
        const double scale = (1.0 - primary_fraction) / secondary_total;
        for (const auto &contact : secondary.contacts) {
            if (contact.bin1 < bin_count && contact.bin2 < bin_count) {
                add_contact_weight(combined, contact.bin1, contact.bin2, contact.weight * scale);
            }
        }
    }

    if (combined.empty()) {
        if (primary_total > 0.0) {
            return primary;
        }
        return secondary;
    }

    return finalize_sparse_contacts(bin_count, combined);
}

ContactMatrix apply_trans_ratio(const ContactMatrix &matrix,
                                const std::vector<OffsetEntry> &offsets,
                                double target_trans_ratio) {
    if (target_trans_ratio < 0.0 || target_trans_ratio > 1.0) {
        throw std::runtime_error("target_trans_ratio must be within [0, 1].");
    }

    if (matrix.bin_count == 0 || matrix.contacts.empty()) {
        return matrix;
    }

    const std::vector<int> bin_to_contig = build_bin_to_contig(matrix.bin_count, offsets);
    ContactWeightMap combined;
    combined.reserve(matrix.contacts.size() * 2);

    double cis_sum = 0.0;
    double trans_sum = 0.0;
    for (const auto &contact : matrix.contacts) {
        if (contact.bin1 >= matrix.bin_count || contact.bin2 >= matrix.bin_count || contact.weight <= 0.0) {
            continue;
        }
        const int contig1 = bin_to_contig[contact.bin1];
        const int contig2 = bin_to_contig[contact.bin2];
        if (contig1 >= 0 && contig1 == contig2) {
            cis_sum += contact.weight;
        } else {
            trans_sum += contact.weight;
        }
    }

    const double total_sum = cis_sum + trans_sum;
    if (total_sum <= std::numeric_limits<double>::min()) {
        return matrix;
    }

    const double desired_trans = total_sum * target_trans_ratio;
    const double desired_cis = total_sum - desired_trans;
    if (desired_trans > 0.0 && trans_sum <= std::numeric_limits<double>::min()) {
        throw std::runtime_error("Cannot apply --trans-ratio because the matrix has no trans contacts. "
                                 "Use an empirical model with trans contacts.");
    }
    if (desired_cis > 0.0 && cis_sum <= std::numeric_limits<double>::min()) {
        throw std::runtime_error("Cannot apply --trans-ratio because the matrix has no cis contacts.");
    }
    const double cis_scale = cis_sum > 0.0 ? desired_cis / cis_sum : 0.0;
    const double trans_scale = trans_sum > 0.0 ? desired_trans / trans_sum : 0.0;

    for (const auto &contact : matrix.contacts) {
        if (contact.bin1 >= matrix.bin_count || contact.bin2 >= matrix.bin_count || contact.weight <= 0.0) {
            continue;
        }
        const int contig1 = bin_to_contig[contact.bin1];
        const int contig2 = bin_to_contig[contact.bin2];
        const bool is_cis = contig1 >= 0 && contig1 == contig2;
        const double scaled_weight = contact.weight * (is_cis ? cis_scale : trans_scale);
        add_contact_weight(combined, contact.bin1, contact.bin2, scaled_weight);
    }

    return finalize_sparse_contacts(matrix.bin_count, combined);
}

ContactMatrix replace_trans_contacts(const ContactMatrix &base,
                                     const ContactMatrix &trans_template,
                                     const std::vector<OffsetEntry> &offsets,
                                     double target_trans_ratio,
                                     bool target_trans_ratio_explicit) {
    if (base.bin_count == 0 || base.contacts.empty()) {
        return base;
    }
    if (target_trans_ratio < 0.0 || target_trans_ratio > 1.0) {
        throw std::runtime_error("target_trans_ratio must be within [0, 1].");
    }
    if (offsets.size() < 2) {
        return base;
    }
    if (trans_template.bin_count != base.bin_count) {
        throw std::runtime_error("Trans template must be remapped to the same bin count before replacement.");
    }

    const std::vector<int> bin_to_contig = build_bin_to_contig(base.bin_count, offsets);
    ContactWeightMap combined;
    combined.reserve(base.contacts.size() + trans_template.contacts.size());

    double base_cis_sum = 0.0;
    double base_trans_sum = 0.0;
    for (const auto &contact : base.contacts) {
        if (contact.bin1 >= base.bin_count || contact.bin2 >= base.bin_count ||
            !std::isfinite(contact.weight) || contact.weight <= 0.0) {
            continue;
        }
        if (contact_is_cis(contact, bin_to_contig)) {
            add_contact_weight(combined, contact.bin1, contact.bin2, contact.weight);
            base_cis_sum += contact.weight;
        } else {
            base_trans_sum += contact.weight;
        }
    }

    double template_trans_sum = 0.0;
    for (const auto &contact : trans_template.contacts) {
        if (contact.bin1 >= trans_template.bin_count || contact.bin2 >= trans_template.bin_count ||
            !std::isfinite(contact.weight) || contact.weight <= 0.0) {
            continue;
        }
        if (!contact_is_cis(contact, bin_to_contig)) {
            template_trans_sum += contact.weight;
        }
    }

    if (template_trans_sum <= std::numeric_limits<double>::min()) {
        throw std::runtime_error("Empirical trans model does not contain positive trans contacts after remapping.");
    }

    double target_trans_sum = base_trans_sum;
    if (target_trans_ratio_explicit) {
        if (target_trans_ratio >= 1.0) {
            target_trans_sum = base_cis_sum > 0.0 ? base_cis_sum : template_trans_sum;
        } else {
            target_trans_sum = base_cis_sum * target_trans_ratio / (1.0 - target_trans_ratio);
        }
    } else if (target_trans_sum <= std::numeric_limits<double>::min()) {
        const double template_cis_sum = sum_contact_weights(trans_template) - template_trans_sum;
        if (template_cis_sum > std::numeric_limits<double>::min() && base_cis_sum > 0.0) {
            target_trans_sum = base_cis_sum * template_trans_sum / template_cis_sum;
        } else {
            target_trans_sum = template_trans_sum;
        }
    }

    if (target_trans_sum <= std::numeric_limits<double>::min()) {
        return finalize_sparse_contacts(base.bin_count, combined);
    }

    const double trans_scale = target_trans_sum / template_trans_sum;
    for (const auto &contact : trans_template.contacts) {
        if (contact.bin1 >= trans_template.bin_count || contact.bin2 >= trans_template.bin_count ||
            !std::isfinite(contact.weight) || contact.weight <= 0.0) {
            continue;
        }
        if (!contact_is_cis(contact, bin_to_contig)) {
            add_contact_weight(combined, contact.bin1, contact.bin2, contact.weight * trans_scale);
        }
    }

    return finalize_sparse_contacts(base.bin_count, combined);
}
