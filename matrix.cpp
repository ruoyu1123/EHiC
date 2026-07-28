#include "matrix.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

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

struct WeightedBin {
    std::size_t bin = 0;
    double fraction = 0.0;
};

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
    const std::string first = token_lower(tokens[0]);
    const std::string second = token_lower(tokens[1]);
    const std::string third = token_lower(tokens[2]);
    if (first == "contig" && second == "start_bin" && third == "end_bin") {
        return true;
    }
    return tokens.size() >= 4 &&
           (first == "name" || first == "contig") &&
           second == "bin_offset" &&
           third == "length" &&
           token_lower(tokens[3]) == "bin_num";
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

struct ContactMassSummary {
    double cis = 0.0;
    double trans = 0.0;
    double total = 0.0;
    double trans_ratio = 0.0;
};

ContactMassSummary summarize_contact_mass(const ContactMatrix &matrix,
                                          const std::vector<int> &bin_to_contig) {
    ContactMassSummary summary;
    for (const auto &contact : matrix.contacts) {
        if (contact.bin1 >= matrix.bin_count || contact.bin2 >= matrix.bin_count ||
            !std::isfinite(contact.weight) || contact.weight <= 0.0) {
            continue;
        }
        const int contig1 = contact.bin1 < bin_to_contig.size() ? bin_to_contig[contact.bin1] : -1;
        const int contig2 = contact.bin2 < bin_to_contig.size() ? bin_to_contig[contact.bin2] : -1;
        if (contig1 >= 0 && contig1 == contig2) {
            summary.cis += contact.weight;
        } else {
            summary.trans += contact.weight;
        }
    }
    summary.total = summary.cis + summary.trans;
    if (summary.total > std::numeric_limits<double>::min()) {
        summary.trans_ratio = summary.trans / summary.total;
    }
    return summary;
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

std::vector<WeightedBin> build_weighted_bin_mapping(double mapped_center,
                                                    std::size_t target_start_bin,
                                                    std::size_t target_span) {
    if (target_span == 0) {
        throw std::runtime_error("Target bin span must be positive.");
    }
    if (target_span == 1) {
        return std::vector<WeightedBin>{{target_start_bin, 1.0}};
    }

    const double clamped_center =
        std::max(0.0, std::min(mapped_center, static_cast<double>(target_span - 1)));
    const std::size_t lower_local = static_cast<std::size_t>(std::floor(clamped_center));
    const double upper_fraction = clamped_center - static_cast<double>(lower_local);
    std::vector<WeightedBin> mapped_bins;
    mapped_bins.push_back(WeightedBin{target_start_bin + lower_local, 1.0 - upper_fraction});
    if (upper_fraction > std::numeric_limits<double>::epsilon() && lower_local + 1 < target_span) {
        mapped_bins.push_back(WeightedBin{target_start_bin + lower_local + 1, upper_fraction});
    } else {
        mapped_bins.front().fraction = 1.0;
    }
    return mapped_bins;
}

std::vector<WeightedBin> map_bin_by_scale(std::size_t source_bin,
                                          std::size_t source_bin_count,
                                          std::size_t target_bin_count) {
    if (target_bin_count == 0) {
        throw std::runtime_error("Target bin count must be positive.");
    }
    if (target_bin_count == 1 || source_bin_count <= 1) {
        return std::vector<WeightedBin>{{0, 1.0}};
    }

    const double scaled_center =
        (static_cast<double>(source_bin) + 0.5) * static_cast<double>(target_bin_count) /
        static_cast<double>(source_bin_count) - 0.5;
    return build_weighted_bin_mapping(scaled_center, 0, target_bin_count);
}

std::vector<WeightedBin> map_bin_within_contig(std::size_t source_bin,
                                               const OffsetEntry &source_offset,
                                               const OffsetEntry &target_offset) {
    const std::size_t source_span = source_offset.end_bin - source_offset.start_bin;
    const std::size_t target_span = target_offset.end_bin - target_offset.start_bin;
    if (target_span == 0) {
        throw std::runtime_error("Target offset span must be positive.");
    }
    if (target_span == 1 || source_span <= 1) {
        return std::vector<WeightedBin>{{target_offset.start_bin, 1.0}};
    }

    const std::size_t local_source = source_bin - source_offset.start_bin;
    const double scaled_center =
        (static_cast<double>(local_source) + 0.5) * static_cast<double>(target_span) /
        static_cast<double>(source_span) - 0.5;
    return build_weighted_bin_mapping(scaled_center, target_offset.start_bin, target_span);
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
            throw std::runtime_error("Offset rows must contain either: contig start_bin end_bin "
                                     "or name bin_offset length bin_num");
        }
        if (is_offset_header(tokens)) {
            continue;
        }

        OffsetEntry entry;
        entry.contig = tokens[0];
        entry.start_bin = std::stoull(tokens[1]);
        if (tokens.size() >= 4) {
            const std::size_t bin_num = std::stoull(tokens[3]);
            if (bin_num > std::numeric_limits<std::size_t>::max() - entry.start_bin) {
                throw std::runtime_error("Offset bin_offset + bin_num overflows size_t.");
            }
            entry.end_bin = entry.start_bin + bin_num;
        } else {
            entry.end_bin = std::stoull(tokens[2]);
        }
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
    std::unordered_set<std::string> contig_names;
    contig_names.reserve(offsets.size());
    for (std::size_t i = 0; i < offsets.size(); ++i) {
        if (!contig_names.insert(offsets[i].contig).second) {
            throw std::runtime_error("Offset contig names must be unique: " + offsets[i].contig);
        }
        if (i == 0) {
            continue;
        }
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

bool offsets_match_reference_exactly(const std::vector<OffsetEntry> &source_offsets,
                                     const std::vector<OffsetEntry> &target_offsets,
                                     std::string *reason) {
    auto set_reason = [&](const std::string &message) {
        if (reason != nullptr) {
            *reason = message;
        }
    };

    if (source_offsets.empty()) {
        set_reason("source offsets are empty");
        return false;
    }

    if (source_offsets.size() != target_offsets.size()) {
        set_reason("contig count differs between matrix offsets (" +
                   std::to_string(source_offsets.size()) + ") and target reference (" +
                   std::to_string(target_offsets.size()) + ")");
        return false;
    }

    for (std::size_t i = 0; i < source_offsets.size(); ++i) {
        const OffsetEntry &source = source_offsets[i];
        const OffsetEntry &target = target_offsets[i];
        if (source.contig != target.contig) {
            set_reason("contig name/order differs at index " + std::to_string(i) +
                       ": matrix has '" + source.contig + "', target has '" + target.contig + "'");
            return false;
        }
        if (source.start_bin != target.start_bin || source.end_bin != target.end_bin) {
            set_reason("bin range differs for contig '" + source.contig + "': matrix has [" +
                       std::to_string(source.start_bin) + ", " + std::to_string(source.end_bin) +
                       "), target has [" + std::to_string(target.start_bin) + ", " +
                       std::to_string(target.end_bin) + ")");
            return false;
        }
    }

    set_reason("offsets match target reference exactly");
    return true;
}

ContactMatrix remap_matrix_to_reference(const ContactMatrix &source_matrix,
                                        const std::vector<OffsetEntry> &source_offsets,
                                        const std::vector<OffsetEntry> &target_offsets,
                                        std::uint64_t seed,
                                        std::vector<std::string> *warnings,
                                        bool force_contig_reuse) {
    const std::size_t target_bin_count = total_offset_bins(target_offsets);
    if (target_bin_count == 0) {
        throw std::runtime_error("Reference offsets are empty.");
    }
    if (source_matrix.bin_count == 0 || source_matrix.contacts.empty()) {
        return ContactMatrix{target_bin_count, {}};
    }

    const std::vector<int> source_bin_to_contig =
        build_contig_lookup_by_bin(source_matrix.bin_count, source_offsets);
    std::unordered_map<std::string, std::size_t> source_index_by_name;
    source_index_by_name.reserve(source_offsets.size());
    for (std::size_t i = 0; i < source_offsets.size(); ++i) {
        source_index_by_name[source_offsets[i].contig] = i;
    }

    std::vector<bool> source_has_trans(source_offsets.size(), false);
    for (const auto &contact : source_matrix.contacts) {
        if (contact.bin1 >= source_bin_to_contig.size() || contact.bin2 >= source_bin_to_contig.size() ||
            !std::isfinite(contact.weight) || contact.weight <= 0.0) {
            continue;
        }
        const int contig1 = source_bin_to_contig[contact.bin1];
        const int contig2 = source_bin_to_contig[contact.bin2];
        if (contig1 >= 0 && contig2 >= 0 && contig1 != contig2) {
            source_has_trans[static_cast<std::size_t>(contig1)] = true;
            source_has_trans[static_cast<std::size_t>(contig2)] = true;
        }
    }

    std::vector<std::vector<std::size_t>> targets_by_source(source_offsets.size());
    if (!source_offsets.empty()) {
        std::vector<bool> source_assigned(source_offsets.size(), false);
        std::vector<bool> target_assigned(target_offsets.size(), false);

        // Prefer same-name contigs, then pair the remaining entries by offset order.
        for (std::size_t target_index = 0; target_index < target_offsets.size(); ++target_index) {
            const auto source_it = source_index_by_name.find(target_offsets[target_index].contig);
            if (source_it != source_index_by_name.end() && !source_assigned[source_it->second]) {
                targets_by_source[source_it->second].push_back(target_index);
                source_assigned[source_it->second] = true;
                target_assigned[target_index] = true;
            }
        }
        std::size_t next_source = 0;
        std::vector<std::size_t> reusable_sources;
        for (std::size_t source_index = 0; source_index < source_offsets.size(); ++source_index) {
            if (source_has_trans[source_index]) {
                reusable_sources.push_back(source_index);
            }
        }
        std::mt19937_64 reuse_rng(seed);
        std::shuffle(reusable_sources.begin(), reusable_sources.end(), reuse_rng);
        std::size_t reuse_cursor = 0;
        for (std::size_t target_index = 0; target_index < target_offsets.size(); ++target_index) {
            if (target_assigned[target_index]) {
                continue;
            }
            while (next_source < source_assigned.size() && source_assigned[next_source]) {
                ++next_source;
            }
            if (next_source == source_assigned.size()) {
                if (!force_contig_reuse) {
                    throw std::runtime_error("Not enough source offset contigs to map the target reference.");
                }
                if (reusable_sources.empty()) {
                    throw std::runtime_error("--force-contig-reuse requires at least two source contigs with "
                                             "positive trans contacts; no reusable trans template is available.");
                }
                const std::size_t reused_source =
                    reusable_sources[reuse_cursor++ % reusable_sources.size()];
                targets_by_source[reused_source].push_back(target_index);
                target_assigned[target_index] = true;
                if (warnings != nullptr) {
                    warnings->push_back("Target contig '" + target_offsets[target_index].contig +
                                        "' reuses source contig '" +
                                        source_offsets[reused_source].contig +
                                        "' under --force-contig-reuse.");
                }
                continue;
            }
            targets_by_source[next_source].push_back(target_index);
            source_assigned[next_source] = true;
            if (warnings != nullptr) {
                warnings->push_back("Target contig '" + target_offsets[target_index].contig +
                                    "' has no matching source name; mapping by offset order to source contig '" +
                                    source_offsets[next_source].contig + "'.");
            }
        }
    }
    ContactWeightMap remapped_weights;
    remapped_weights.reserve(source_matrix.contacts.size() * 2);
    bool warned_dropped_contacts = false;
    bool warned_global_scaling = false;

    for (const auto &contact : source_matrix.contacts) {
        if (contact.weight <= 0.0 || !std::isfinite(contact.weight)) {
            continue;
        }
        if (contact.bin1 >= source_matrix.bin_count || contact.bin2 >= source_matrix.bin_count) {
            continue;
        }

        bool mapped_ok = false;

        const int source_contig1 = contact.bin1 < source_bin_to_contig.size() ? source_bin_to_contig[contact.bin1] : -1;
        const int source_contig2 = contact.bin2 < source_bin_to_contig.size() ? source_bin_to_contig[contact.bin2] : -1;

        if (!source_offsets.empty() && source_contig1 >= 0 && source_contig2 >= 0) {
            const std::size_t source_index1 = static_cast<std::size_t>(source_contig1);
            const std::size_t source_index2 = static_cast<std::size_t>(source_contig2);
            const OffsetEntry &src_offset1 = source_offsets[source_index1];
            const OffsetEntry &src_offset2 = source_offsets[source_index2];
            const auto &target_indices1 = targets_by_source[source_index1];
            const auto &target_indices2 = targets_by_source[source_index2];

            if (source_index1 == source_index2) {
                if (!target_indices1.empty()) {
                    const double split_weight =
                        contact.weight / static_cast<double>(target_indices1.size());
                    for (const std::size_t target_index : target_indices1) {
                        const OffsetEntry &target_offset = target_offsets[target_index];
                        const std::vector<WeightedBin> mapped_bins1 =
                            map_bin_within_contig(contact.bin1, src_offset1, target_offset);
                        const std::vector<WeightedBin> mapped_bins2 =
                            map_bin_within_contig(contact.bin2, src_offset2, target_offset);
                        for (const auto &mapped_bin1 : mapped_bins1) {
                            for (const auto &mapped_bin2 : mapped_bins2) {
                                add_contact_weight(remapped_weights,
                                                   mapped_bin1.bin,
                                                   mapped_bin2.bin,
                                                   split_weight * mapped_bin1.fraction * mapped_bin2.fraction);
                            }
                        }
                    }
                    mapped_ok = true;
                }
            } else if (!target_indices1.empty() && !target_indices2.empty()) {
                const double split_weight =
                    contact.weight /
                    static_cast<double>(target_indices1.size() * target_indices2.size());
                for (const std::size_t target_index1 : target_indices1) {
                    for (const std::size_t target_index2 : target_indices2) {
                        const std::vector<WeightedBin> mapped_bins1 =
                            map_bin_within_contig(contact.bin1, src_offset1, target_offsets[target_index1]);
                        const std::vector<WeightedBin> mapped_bins2 =
                            map_bin_within_contig(contact.bin2, src_offset2, target_offsets[target_index2]);
                        for (const auto &mapped_bin1 : mapped_bins1) {
                            for (const auto &mapped_bin2 : mapped_bins2) {
                                add_contact_weight(remapped_weights,
                                                   mapped_bin1.bin,
                                                   mapped_bin2.bin,
                                                   split_weight * mapped_bin1.fraction * mapped_bin2.fraction);
                            }
                        }
                    }
                }
                mapped_ok = true;
            }

            if (!mapped_ok) {
                if (warnings != nullptr && !warned_dropped_contacts) {
                    warnings->push_back("Dropping contacts from source contigs that are not mapped to the target reference.");
                    warned_dropped_contacts = true;
                }
                continue;
            }
        }

        if (!mapped_ok) {
            const std::vector<WeightedBin> mapped_bins1 =
                map_bin_by_scale(contact.bin1, source_matrix.bin_count, target_bin_count);
            const std::vector<WeightedBin> mapped_bins2 =
                map_bin_by_scale(contact.bin2, source_matrix.bin_count, target_bin_count);
            if (!source_offsets.empty() && warnings != nullptr && !warned_global_scaling) {
                warnings->push_back("Falling back to global bin scaling for some contacts because they could not be mapped through contig offsets.");
                warned_global_scaling = true;
            }
            for (const auto &mapped_bin1 : mapped_bins1) {
                for (const auto &mapped_bin2 : mapped_bins2) {
                    add_contact_weight(remapped_weights,
                                       mapped_bin1.bin,
                                       mapped_bin2.bin,
                                       contact.weight * mapped_bin1.fraction * mapped_bin2.fraction);
                }
            }
        }
    }

    if (force_contig_reuse) {
        const std::vector<int> target_bin_to_contig =
            build_contig_lookup_by_bin(target_bin_count, target_offsets);
        const auto target_is_trans = [&](const ContactKey &key) {
            const int contig1 = key.bin1 < target_bin_to_contig.size()
                                    ? target_bin_to_contig[key.bin1]
                                    : -1;
            const int contig2 = key.bin2 < target_bin_to_contig.size()
                                    ? target_bin_to_contig[key.bin2]
                                    : -1;
            return contig1 >= 0 && contig2 >= 0 && contig1 != contig2;
        };

        double baseline_trans_sum = 0.0;
        for (const auto &[key, weight] : remapped_weights) {
            if (target_is_trans(key) && std::isfinite(weight) && weight > 0.0) {
                baseline_trans_sum += weight;
            }
        }

        std::unordered_map<ContactKey, std::vector<const Contact *>, ContactKeyHash>
            source_trans_by_pair;
        for (const auto &contact : source_matrix.contacts) {
            if (contact.bin1 >= source_bin_to_contig.size() || contact.bin2 >= source_bin_to_contig.size() ||
                !std::isfinite(contact.weight) || contact.weight <= 0.0) {
                continue;
            }
            const int contig1 = source_bin_to_contig[contact.bin1];
            const int contig2 = source_bin_to_contig[contact.bin2];
            if (contig1 >= 0 && contig2 >= 0 && contig1 != contig2) {
                source_trans_by_pair[contact_key(static_cast<std::size_t>(contig1),
                                                 static_cast<std::size_t>(contig2))]
                    .push_back(&contact);
            }
        }

        std::mt19937_64 donor_rng(seed + 32452843ULL);
        std::size_t synthetic_pair_count = 0;
        std::size_t logged_pair_count = 0;
        for (std::size_t source_index = 0; source_index < targets_by_source.size(); ++source_index) {
            const auto &target_indices = targets_by_source[source_index];
            if (target_indices.size() < 2) {
                continue;
            }

            std::vector<std::size_t> donor_sources;
            for (std::size_t donor_index = 0; donor_index < source_offsets.size(); ++donor_index) {
                if (donor_index != source_index &&
                    source_trans_by_pair.find(contact_key(source_index, donor_index)) !=
                        source_trans_by_pair.end()) {
                    donor_sources.push_back(donor_index);
                }
            }
            if (donor_sources.empty()) {
                throw std::runtime_error("Reused source contig '" + source_offsets[source_index].contig +
                                         "' has no positive trans block with another source contig.");
            }
            std::shuffle(donor_sources.begin(), donor_sources.end(), donor_rng);

            std::size_t donor_cursor = 0;
            for (std::size_t left_pos = 0; left_pos < target_indices.size(); ++left_pos) {
                for (std::size_t right_pos = left_pos + 1; right_pos < target_indices.size(); ++right_pos) {
                    const std::size_t donor_index =
                        donor_sources[donor_cursor++ % donor_sources.size()];
                    const auto pair_it =
                        source_trans_by_pair.find(contact_key(source_index, donor_index));
                    const OffsetEntry &source_offset = source_offsets[source_index];
                    const OffsetEntry &donor_offset = source_offsets[donor_index];
                    const OffsetEntry &left_target = target_offsets[target_indices[left_pos]];
                    const OffsetEntry &right_target = target_offsets[target_indices[right_pos]];
                    const double donor_pair_scale =
                        1.0 / static_cast<double>(target_indices.size() *
                                                  targets_by_source[donor_index].size());

                    for (const Contact *contact : pair_it->second) {
                        const int contact_contig1 = source_bin_to_contig[contact->bin1];
                        const bool source_is_first =
                            contact_contig1 == static_cast<int>(source_index);
                        const std::size_t source_bin =
                            source_is_first ? contact->bin1 : contact->bin2;
                        const std::size_t donor_bin =
                            source_is_first ? contact->bin2 : contact->bin1;
                        const auto mapped_source_bins =
                            map_bin_within_contig(source_bin, source_offset, left_target);
                        const auto mapped_donor_bins =
                            map_bin_within_contig(donor_bin, donor_offset, right_target);
                        for (const auto &mapped_source : mapped_source_bins) {
                            for (const auto &mapped_donor : mapped_donor_bins) {
                                add_contact_weight(remapped_weights,
                                                   mapped_source.bin,
                                                   mapped_donor.bin,
                                                   contact->weight * donor_pair_scale *
                                                       mapped_source.fraction * mapped_donor.fraction);
                            }
                        }
                    }

                    ++synthetic_pair_count;
                    if (warnings != nullptr && logged_pair_count < 20) {
                        warnings->push_back("Synthetic trans block '" + left_target.contig +
                                            "' <-> '" + right_target.contig + "' uses source trans '" +
                                            source_offsets[source_index].contig + "' <-> '" +
                                            source_offsets[donor_index].contig + "'.");
                        ++logged_pair_count;
                    }
                }
            }
        }

        if (synthetic_pair_count > 0) {
            if (baseline_trans_sum <= std::numeric_limits<double>::min()) {
                throw std::runtime_error("Cannot normalize forced contig reuse because the remapped matrix "
                                         "has no positive baseline trans mass.");
            }
            double expanded_trans_sum = 0.0;
            for (const auto &[key, weight] : remapped_weights) {
                if (target_is_trans(key) && std::isfinite(weight) && weight > 0.0) {
                    expanded_trans_sum += weight;
                }
            }
            const double trans_scale = baseline_trans_sum / expanded_trans_sum;
            for (auto &[key, weight] : remapped_weights) {
                if (target_is_trans(key)) {
                    weight *= trans_scale;
                }
            }
            if (warnings != nullptr) {
                if (synthetic_pair_count > logged_pair_count) {
                    warnings->push_back("Additional synthetic trans blocks omitted from detailed log: " +
                                        std::to_string(synthetic_pair_count - logged_pair_count) + ".");
                }
                warnings->push_back("Forced contig reuse synthesized " +
                                    std::to_string(synthetic_pair_count) +
                                    " duplicate-target trans block(s); trans mass normalized from " +
                                    std::to_string(expanded_trans_sum) + " to " +
                                    std::to_string(baseline_trans_sum) + ".");
            }
        }
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
    std::cerr << "Applying trans ratio override: current cis=" << cis_sum
              << ", current trans=" << trans_sum
              << ", target trans ratio=" << target_trans_ratio
              << ", cis scale=" << cis_scale
              << ", trans scale=" << trans_scale << '\n';

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

    const ContactMassSummary base_summary = summarize_contact_mass(base, bin_to_contig);
    const ContactMassSummary template_summary = summarize_contact_mass(trans_template, bin_to_contig);

    double base_cis_sum = 0.0;
    for (const auto &contact : base.contacts) {
        if (contact.bin1 >= base.bin_count || contact.bin2 >= base.bin_count ||
            !std::isfinite(contact.weight) || contact.weight <= 0.0) {
            continue;
        }
        if (contact_is_cis(contact, bin_to_contig)) {
            add_contact_weight(combined, contact.bin1, contact.bin2, contact.weight);
            base_cis_sum += contact.weight;
        }
    }

    if (template_summary.trans <= std::numeric_limits<double>::min()) {
        throw std::runtime_error("Empirical trans model does not contain positive trans contacts after remapping.");
    }

    const double template_trans_ratio = template_summary.trans_ratio;
    const double template_cis_ratio = 1.0 - template_trans_ratio;
    if (template_trans_ratio <= std::numeric_limits<double>::min()) {
        throw std::runtime_error("Empirical trans model has no trans ratio after remapping.");
    }

    double target_trans_sum = base_summary.trans;
    if (target_trans_ratio_explicit) {
        if (target_trans_ratio >= 1.0) {
            target_trans_sum = base_cis_sum > 0.0 ? base_cis_sum : template_summary.trans;
        } else {
            target_trans_sum = base_cis_sum * target_trans_ratio / (1.0 - target_trans_ratio);
        }
    } else if (target_trans_sum <= std::numeric_limits<double>::min()) {
        if (base_cis_sum > 0.0 && template_cis_ratio > std::numeric_limits<double>::min()) {
            target_trans_sum = base_cis_sum * template_trans_ratio / template_cis_ratio;
        } else {
            target_trans_sum = template_summary.trans;
        }
    }

    if (target_trans_sum <= std::numeric_limits<double>::min()) {
        return finalize_sparse_contacts(base.bin_count, combined);
    }

    const double trans_scale = target_trans_sum / template_summary.trans;
    std::cerr << "Trans replacement summary: base cis=" << base_summary.cis
              << ", base trans=" << base_summary.trans
              << ", base trans ratio=" << base_summary.trans_ratio
              << ", template cis=" << template_summary.cis
              << ", template trans=" << template_summary.trans
              << ", template trans ratio=" << template_trans_ratio
              << ", target trans sum=" << target_trans_sum
              << ", trans scale=" << trans_scale;
    if (target_trans_ratio_explicit) {
        std::cerr << ", source=explicit --trans-ratio (" << target_trans_ratio << ")";
    } else if (base_summary.trans > std::numeric_limits<double>::min()) {
        std::cerr << ", source=input matrix trans mass";
    } else {
        std::cerr << ", source=template trans ratio fallback";
    }
    std::cerr << '\n';

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
