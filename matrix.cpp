#include "matrix.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace {

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

ContactMatrix load_matrix(const std::string &path, std::size_t expected_bin_count) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open matrix file: " + path);
    }

    std::string line;
    std::vector<std::vector<double>> dense_rows;
    std::unordered_map<std::uint64_t, double> sparse_weights;
    std::size_t detected_columns = 0;
    bool sparse_mode = false;
    bool mode_decided = false;
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

            const std::size_t lo = std::min(bin1, bin2);
            const std::size_t hi = std::max(bin1, bin2);
            const std::uint64_t key =
                (static_cast<std::uint64_t>(lo) << 32U) | static_cast<std::uint64_t>(hi);
            sparse_weights[key] += weight;
            max_bin = std::max(max_bin, hi);
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
    matrix.bin_count = expected_bin_count == 0 ? max_bin + 1 : expected_bin_count;
    for (const auto &[key, weight] : sparse_weights) {
        const std::size_t bin1 = static_cast<std::size_t>(key >> 32U);
        const std::size_t bin2 = static_cast<std::size_t>(key & 0xffffffffULL);
        if (bin1 >= matrix.bin_count || bin2 >= matrix.bin_count) {
            throw std::runtime_error("Sparse matrix bin index exceeds expected bin count.");
        }
        matrix.contacts.push_back(Contact{bin1, bin2, weight});
    }

    if (matrix.contacts.empty()) {
        throw std::runtime_error("Sparse matrix does not contain any positive contacts.");
    }

    return matrix;
}
