#include "matrix.h"
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <random>
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

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

enum class ArrangementMode {
    Territory,
    Rabl,
    Rosette,
    NonRabl
};

struct SpeciesProfile {
    ArrangementMode arrangement = ArrangementMode::Territory;
    double default_trans_ratio = 0.10;
    double default_cis_decay_alpha = 1.0;
    double default_collision_randomness = 0.35;
    double rabl_arm_bias = 0.0;
    double rosette_center_bias = 0.0;
    double subtelomere_trans_bias = 0.0;
    double trans_distance_gamma = 1.2;
};

struct ResolvedSyntheticOptions {
    ArrangementMode arrangement = ArrangementMode::Territory;
    double trans_ratio = 0.10;
    double cis_decay_alpha = 1.0;
    double collision_randomness = 0.35;
    std::size_t max_cis_distance_bins = 200;
    double rabl_arm_bias = 0.0;
    double rosette_center_bias = 0.0;
    double subtelomere_trans_bias = 0.0;
    double trans_distance_gamma = 1.2;
    std::uint64_t seed = 1;
};

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

SpeciesProfile species_profile_for(const std::string &species_model) {
    const std::string species = to_lower(species_model);
    if (species == "arabidopsis" || species == "athaliana" || species == "arabidopsis_thaliana") {
        return SpeciesProfile{
            ArrangementMode::Rosette, 0.06, 1.08, 0.30, 0.0, 0.35, 0.05, 1.1
        };
    }
    if (species == "wheat" || species == "triticum" || species == "triticum_aestivum") {
        return SpeciesProfile{
            ArrangementMode::Rabl, 0.14, 0.95, 0.28, 0.28, 0.0, 0.45, 1.35
        };
    }
    if (species == "barley" || species == "hordeum" || species == "hordeum_vulgare") {
        return SpeciesProfile{
            ArrangementMode::Rabl, 0.12, 0.95, 0.30, 0.24, 0.0, 0.40, 1.30
        };
    }
    if (species == "rice" || species == "oryza" || species == "oryza_sativa") {
        return SpeciesProfile{
            ArrangementMode::NonRabl, 0.08, 1.03, 0.36, 0.04, 0.08, 0.10, 1.2
        };
    }
    if (species == "maize" || species == "zea" || species == "zea_mays") {
        return SpeciesProfile{
            ArrangementMode::Rabl, 0.11, 1.0, 0.32, 0.16, 0.0, 0.26, 1.25
        };
    }
    return SpeciesProfile{};
}

ArrangementMode parse_arrangement(const std::string &arrangement_model, ArrangementMode fallback) {
    const std::string model = to_lower(arrangement_model);
    if (model.empty() || model == "auto") {
        return fallback;
    }
    if (model == "territory") {
        return ArrangementMode::Territory;
    }
    if (model == "rabl") {
        return ArrangementMode::Rabl;
    }
    if (model == "rosette") {
        return ArrangementMode::Rosette;
    }
    if (model == "nonrabl" || model == "non-rabl") {
        return ArrangementMode::NonRabl;
    }
    throw std::runtime_error("Unknown arrangement model: " + arrangement_model);
}

ResolvedSyntheticOptions resolve_synthetic_options(const SyntheticModelOptions &options) {
    if (options.trans_ratio < 0.0 || options.trans_ratio > 1.0) {
        throw std::runtime_error("trans_ratio must be within [0, 1].");
    }
    if (!std::isfinite(options.cis_decay_alpha) || options.cis_decay_alpha <= 0.0) {
        throw std::runtime_error("cis_decay_alpha must be a positive finite value.");
    }
    if (options.collision_randomness < 0.0 || options.collision_randomness > 1.0) {
        throw std::runtime_error("collision_randomness must be within [0, 1].");
    }

    const SpeciesProfile profile = species_profile_for(options.species_model);
    ResolvedSyntheticOptions resolved;
    resolved.arrangement = parse_arrangement(options.arrangement_model, profile.arrangement);
    resolved.trans_ratio =
        std::abs(options.trans_ratio - 0.10) < 1e-12 ? profile.default_trans_ratio : options.trans_ratio;
    resolved.cis_decay_alpha =
        std::abs(options.cis_decay_alpha - 1.0) < 1e-12 ? profile.default_cis_decay_alpha : options.cis_decay_alpha;
    resolved.collision_randomness =
        std::abs(options.collision_randomness - 0.35) < 1e-12 ? profile.default_collision_randomness : options.collision_randomness;
    resolved.max_cis_distance_bins = options.max_cis_distance_bins;
    resolved.rabl_arm_bias = profile.rabl_arm_bias;
    resolved.rosette_center_bias = profile.rosette_center_bias;
    resolved.subtelomere_trans_bias = profile.subtelomere_trans_bias;
    resolved.trans_distance_gamma = profile.trans_distance_gamma;
    resolved.seed = options.seed;
    return resolved;
}

double squared_distance(const Vec3 &a, const Vec3 &b) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

std::vector<Vec3> build_chromosome_centers(const std::vector<OffsetEntry> &offsets,
                                           ArrangementMode arrangement,
                                           std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::uniform_real_distribution<double> u01(0.0, 1.0);

    std::vector<Vec3> centers(offsets.size());
    const double pi = 3.14159265358979323846;
    for (std::size_t i = 0; i < offsets.size(); ++i) {
        if (arrangement == ArrangementMode::Rabl) {
            const double theta = 2.0 * pi * (static_cast<double>(i) + 0.5) / static_cast<double>(offsets.size());
            centers[i] = Vec3{0.25 * u(rng), std::cos(theta), 0.55 * std::sin(theta)};
            continue;
        }
        if (arrangement == ArrangementMode::Rosette) {
            const double theta = 2.0 * pi * (static_cast<double>(i) + 0.5) / static_cast<double>(offsets.size());
            const double radius = 0.75 + 0.15 * u01(rng);
            centers[i] = Vec3{radius * std::cos(theta), radius * std::sin(theta), 0.20 * u(rng)};
            continue;
        }
        const double theta = 2.0 * pi * u01(rng);
        const double phi = std::acos(std::max(-1.0, std::min(1.0, u(rng))));
        const double radius = std::cbrt(u01(rng));
        centers[i] = Vec3{
            radius * std::sin(phi) * std::cos(theta),
            radius * std::sin(phi) * std::sin(theta),
            radius * std::cos(phi)
        };
    }
    return centers;
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

std::uint64_t contact_key(std::size_t bin1, std::size_t bin2) {
    const std::size_t lo = std::min(bin1, bin2);
    const std::size_t hi = std::max(bin1, bin2);
    return (static_cast<std::uint64_t>(lo) << 32U) | static_cast<std::uint64_t>(hi);
}

void add_contact_weight(std::unordered_map<std::uint64_t, double> &weights,
                        std::size_t bin1,
                        std::size_t bin2,
                        double value) {
    if (!std::isfinite(value) || value <= 0.0) {
        return;
    }
    weights[contact_key(bin1, bin2)] += value;
}

ContactMatrix finalize_sparse_contacts(std::size_t bin_count,
                                       const std::unordered_map<std::uint64_t, double> &weights) {
    ContactMatrix matrix;
    matrix.bin_count = bin_count;
    matrix.contacts.reserve(weights.size());
    for (const auto &[key, weight] : weights) {
        if (!std::isfinite(weight) || weight <= 0.0) {
            continue;
        }
        const std::size_t bin1 = static_cast<std::size_t>(key >> 32U);
        const std::size_t bin2 = static_cast<std::size_t>(key & 0xffffffffULL);
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

std::pair<std::size_t, std::size_t> sample_trans_contact(const std::vector<OffsetEntry> &offsets,
                                                         const std::vector<std::array<std::size_t, 2>> &contig_pairs,
                                                         std::discrete_distribution<std::size_t> &pair_dist,
                                                         const ResolvedSyntheticOptions &options,
                                                         std::mt19937_64 &rng) {
    if (offsets.size() < 2) {
        throw std::runtime_error("Cannot sample trans contacts with fewer than two contigs.");
    }

    const auto pair = contig_pairs[pair_dist(rng)];
    const std::size_t left_contig = pair[0];
    const std::size_t right_contig = pair[1];

    const OffsetEntry &left_offset = offsets[left_contig];
    const OffsetEntry &right_offset = offsets[right_contig];
    const auto sample_from_offset = [&](const OffsetEntry &offset) -> std::size_t {
        const std::size_t span = offset.end_bin - offset.start_bin;
        if (span == 1) {
            return offset.start_bin;
        }
        std::uniform_int_distribution<std::size_t> uniform_dist(offset.start_bin, offset.end_bin - 1);
        if (options.arrangement == ArrangementMode::Rabl && options.subtelomere_trans_bias > 0.0) {
            if (std::bernoulli_distribution(options.subtelomere_trans_bias)(rng)) {
                const std::size_t window = std::max<std::size_t>(1, span / 6);
                std::uniform_int_distribution<std::size_t> edge_dist(0, window - 1);
                if (std::bernoulli_distribution(0.5)(rng)) {
                    return offset.start_bin + edge_dist(rng);
                }
                return (offset.end_bin - 1) - edge_dist(rng);
            }
        }
        if (options.arrangement == ArrangementMode::Rosette && options.rosette_center_bias > 0.0) {
            if (std::bernoulli_distribution(options.rosette_center_bias * 0.5)(rng)) {
                const std::size_t center = offset.start_bin + span / 2;
                const std::size_t half_window = std::max<std::size_t>(1, span / 8);
                const std::size_t lo = center > half_window ? center - half_window : offset.start_bin;
                const std::size_t hi = std::min(offset.end_bin - 1, center + half_window);
                std::uniform_int_distribution<std::size_t> center_dist(lo, hi);
                return center_dist(rng);
            }
        }
        return uniform_dist(rng);
    };

    return {sample_from_offset(left_offset), sample_from_offset(right_offset)};
}

std::pair<std::size_t, std::size_t> sample_cis_contact(const OffsetEntry &offset,
                                                       std::discrete_distribution<std::size_t> &distance_dist,
                                                       const ResolvedSyntheticOptions &options,
                                                       std::mt19937_64 &rng) {
    const std::size_t contig_bins = offset.end_bin - offset.start_bin;
    if (contig_bins == 0) {
        throw std::runtime_error("Cannot sample cis contacts on an empty contig range.");
    }
    if (contig_bins == 1) {
        return {offset.start_bin, offset.start_bin};
    }

    if (options.arrangement == ArrangementMode::Rabl &&
        options.rabl_arm_bias > 0.0 &&
        std::bernoulli_distribution(options.rabl_arm_bias)(rng)) {
        const std::size_t mid = contig_bins / 2;
        if (mid > 0 && contig_bins - mid > 0) {
            std::uniform_int_distribution<std::size_t> left_half(0, mid - 1);
            std::uniform_int_distribution<std::size_t> right_half(mid, contig_bins - 1);
            return {offset.start_bin + left_half(rng), offset.start_bin + right_half(rng)};
        }
    }

    if (options.arrangement == ArrangementMode::Rosette &&
        options.rosette_center_bias > 0.0 &&
        std::bernoulli_distribution(options.rosette_center_bias)(rng)) {
        const std::size_t center = contig_bins / 2;
        const std::size_t half_window = std::max<std::size_t>(1, contig_bins / 8);
        const std::size_t lo = center > half_window ? center - half_window : 0;
        const std::size_t hi = std::min(contig_bins - 1, center + half_window);
        std::uniform_int_distribution<std::size_t> center_dist(lo, hi);
        return {offset.start_bin + center_dist(rng), offset.start_bin + center_dist(rng)};
    }

    const std::size_t distance = distance_dist(rng);

    std::uniform_int_distribution<std::size_t> start_dist(0, contig_bins - 1 - distance);
    const std::size_t left_local = start_dist(rng);
    const std::size_t right_local = left_local + distance;
    return {offset.start_bin + left_local, offset.start_bin + right_local};
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
    matrix.bin_count = max_bin + 1;
    for (const auto &[key, weight] : sparse_weights) {
        const std::size_t bin1 = static_cast<std::size_t>(key >> 32U);
        const std::size_t bin2 = static_cast<std::size_t>(key & 0xffffffffULL);
        matrix.contacts.push_back(Contact{bin1, bin2, weight});
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

    std::unordered_map<std::uint64_t, double> remapped_weights;
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
    std::unordered_map<std::uint64_t, double> combined;
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

ContactMatrix generate_synthetic_matrix(std::size_t bin_count,
                                        const std::vector<OffsetEntry> &offsets,
                                        std::size_t synthetic_contact_count,
                                        const SyntheticModelOptions &options) {
    if (bin_count == 0) {
        throw std::runtime_error("Synthetic matrix generation requires positive bin_count.");
    }
    if (offsets.empty()) {
        throw std::runtime_error("Synthetic matrix generation requires contig offsets.");
    }
    const ResolvedSyntheticOptions resolved = resolve_synthetic_options(options);

    std::vector<OffsetEntry> valid_offsets;
    valid_offsets.reserve(offsets.size());
    for (const auto &offset : offsets) {
        if (offset.start_bin < offset.end_bin && offset.end_bin <= bin_count) {
            valid_offsets.push_back(offset);
        }
    }
    if (valid_offsets.empty()) {
        throw std::runtime_error("No valid contig offsets are available for synthetic matrix generation.");
    }

    if (synthetic_contact_count == 0) {
        const std::size_t auto_count = std::max<std::size_t>(10000, std::min<std::size_t>(2000000, bin_count * 100));
        synthetic_contact_count = auto_count;
    }

    std::vector<double> cis_contig_weights;
    cis_contig_weights.reserve(valid_offsets.size());
    for (const auto &offset : valid_offsets) {
        const double bins = static_cast<double>(offset.end_bin - offset.start_bin);
        cis_contig_weights.push_back(std::max(0.0, bins));
    }
    std::discrete_distribution<std::size_t> cis_contig_dist(cis_contig_weights.begin(), cis_contig_weights.end());

    std::vector<std::discrete_distribution<std::size_t>> cis_distance_dists;
    cis_distance_dists.reserve(valid_offsets.size());
    for (const auto &offset : valid_offsets) {
        const std::size_t contig_bins = offset.end_bin - offset.start_bin;
        const std::size_t max_distance =
            contig_bins > 1 ? std::min(resolved.max_cis_distance_bins, contig_bins - 1) : 0;
        std::vector<double> distance_weights(max_distance + 1, 1.0);
        for (std::size_t d = 0; d <= max_distance; ++d) {
            distance_weights[d] = 1.0 / std::pow(static_cast<double>(d + 1), resolved.cis_decay_alpha);
        }
        cis_distance_dists.emplace_back(distance_weights.begin(), distance_weights.end());
    }

    std::mt19937_64 rng(resolved.seed + 311);
    std::bernoulli_distribution use_trans(resolved.trans_ratio);

    const std::vector<Vec3> centers = build_chromosome_centers(valid_offsets, resolved.arrangement, rng);
    std::vector<std::array<std::size_t, 2>> trans_pairs;
    std::vector<double> trans_pair_weights;
    for (std::size_t i = 0; i < valid_offsets.size(); ++i) {
        for (std::size_t j = i + 1; j < valid_offsets.size(); ++j) {
            const double bins_i = static_cast<double>(valid_offsets[i].end_bin - valid_offsets[i].start_bin);
            const double bins_j = static_cast<double>(valid_offsets[j].end_bin - valid_offsets[j].start_bin);
            const double d = std::sqrt(squared_distance(centers[i], centers[j]));
            const double arranged = 1.0 / std::pow(d + 0.25, resolved.trans_distance_gamma);
            const double collision = 1.0;
            const double mixed = resolved.collision_randomness * collision +
                                 (1.0 - resolved.collision_randomness) * arranged;
            trans_pairs.push_back({i, j});
            trans_pair_weights.push_back(std::max(0.0, bins_i * bins_j * mixed));
        }
    }
    std::discrete_distribution<std::size_t> trans_pair_dist(trans_pair_weights.begin(), trans_pair_weights.end());

    std::unordered_map<std::uint64_t, double> weights;
    weights.reserve(synthetic_contact_count * 2);

    for (std::size_t i = 0; i < synthetic_contact_count; ++i) {
        bool do_trans = use_trans(rng) && valid_offsets.size() > 1;
        if (do_trans) {
            const auto [left, right] =
                sample_trans_contact(valid_offsets, trans_pairs, trans_pair_dist, resolved, rng);
            add_contact_weight(weights, left, right, 1.0);
            continue;
        }

        const std::size_t contig_index = cis_contig_dist(rng);
        const auto [left, right] = sample_cis_contact(valid_offsets[contig_index],
                                                      cis_distance_dists[contig_index],
                                                      resolved,
                                                      rng);
        add_contact_weight(weights, left, right, 1.0);
    }

    return finalize_sparse_contacts(bin_count, weights);
}

ContactMatrix apply_trans_ratio(const ContactMatrix &matrix,
                                const std::vector<OffsetEntry> &offsets,
                                double target_trans_ratio,
                                std::uint64_t seed) {
    if (target_trans_ratio < 0.0 || target_trans_ratio > 1.0) {
        throw std::runtime_error("target_trans_ratio must be within [0, 1].");
    }

    if (matrix.bin_count == 0 || matrix.contacts.empty()) {
        return matrix;
    }

    const std::vector<int> bin_to_contig = build_bin_to_contig(matrix.bin_count, offsets);
    std::unordered_map<std::uint64_t, double> combined;
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

    if (desired_trans > 0.0 && trans_sum == 0.0 && offsets.size() > 1) {
        std::mt19937_64 rng(seed + 7331);
        std::vector<double> trans_contig_weights;
        trans_contig_weights.reserve(offsets.size());
        for (const auto &offset : offsets) {
            const double bins = static_cast<double>(offset.end_bin - offset.start_bin);
            trans_contig_weights.push_back(std::max(0.0, bins));
        }
        std::vector<std::array<std::size_t, 2>> trans_pairs;
        std::vector<double> trans_pair_weights;
        for (std::size_t i = 0; i < offsets.size(); ++i) {
            for (std::size_t j = i + 1; j < offsets.size(); ++j) {
                trans_pairs.push_back({i, j});
                trans_pair_weights.push_back(std::max(0.0, trans_contig_weights[i] * trans_contig_weights[j]));
            }
        }
        std::discrete_distribution<std::size_t> trans_pair_dist(trans_pair_weights.begin(), trans_pair_weights.end());
        ResolvedSyntheticOptions trans_fill_options;
        trans_fill_options.arrangement = ArrangementMode::Territory;
        const std::size_t trans_contacts_to_add =
            std::max<std::size_t>(1000, std::min<std::size_t>(200000, static_cast<std::size_t>(matrix.contacts.size() / 5 + 1)));
        const double per_contact = desired_trans / static_cast<double>(trans_contacts_to_add);
        for (std::size_t i = 0; i < trans_contacts_to_add; ++i) {
            const auto [left, right] =
                sample_trans_contact(offsets, trans_pairs, trans_pair_dist, trans_fill_options, rng);
            add_contact_weight(combined, left, right, per_contact);
        }
    }

    if (desired_cis > 0.0 && cis_sum == 0.0 && !offsets.empty()) {
        std::mt19937_64 rng(seed + 9393);
        const std::size_t cis_contacts_to_add =
            std::max<std::size_t>(1000, std::min<std::size_t>(200000, static_cast<std::size_t>(matrix.contacts.size() / 5 + 1)));
        const double per_contact = desired_cis / static_cast<double>(cis_contacts_to_add);

        std::vector<double> contig_weights;
        contig_weights.reserve(offsets.size());
        for (const auto &offset : offsets) {
            const double bins = static_cast<double>(offset.end_bin - offset.start_bin);
            contig_weights.push_back(std::max(0.0, bins));
        }
        std::discrete_distribution<std::size_t> contig_dist(contig_weights.begin(), contig_weights.end());

        for (std::size_t i = 0; i < cis_contacts_to_add; ++i) {
            const OffsetEntry &offset = offsets[contig_dist(rng)];
            if (offset.start_bin >= offset.end_bin) {
                continue;
            }
            std::uniform_int_distribution<std::size_t> left_dist(offset.start_bin, offset.end_bin - 1);
            std::uniform_int_distribution<std::size_t> right_dist(offset.start_bin, offset.end_bin - 1);
            add_contact_weight(combined, left_dist(rng), right_dist(rng), per_contact);
        }
    }

    return finalize_sparse_contacts(matrix.bin_count, combined);
}
