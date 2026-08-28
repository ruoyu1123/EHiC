#include "long_read_common.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <condition_variable>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <random>
#include <stdexcept>
#include <thread>
#include <unordered_map>

namespace {

struct RestrictionFragment {
    std::size_t id = 0;
    std::size_t contig_index = 0;
    std::size_t start = 0;
    std::size_t end = 0;
};

struct FragmentIndex {
    std::vector<RestrictionFragment> fragments;
    std::vector<std::vector<std::size_t>> bin_to_fragments;
};

struct WeightedNeighbor {
    std::size_t bin = 0;
    double cumulative_weight = 0.0;
};

struct ContactSamplingIndex {
    std::vector<Contact> contacts;
    std::vector<double> cumulative_contact_weights;
    double total_contact_weight = 0.0;
    std::vector<std::vector<WeightedNeighbor>> neighbors;
};

struct TruthSegment {
    std::string contig;
    std::size_t start = 0;
    std::size_t end = 0;
    char strand = '+';
    std::size_t bin = 0;
    std::size_t template_start = 0;
    std::size_t template_end = 0;
};

struct ConcatemerTemplate {
    std::string sequence;
    std::vector<TruthSegment> segments;
    bool hit_segment_limit = false;
};

struct SimulatedSequence {
    std::string sequence;
    std::string qualities;
};

struct LongReadBlock {
    std::size_t read_count = 0;
    std::size_t sequenced_bases = 0;
    std::size_t concatemer_segments = 0;
    std::size_t segment_limited_reads = 0;
    std::string fastq;
    std::string truth;
};

std::string normalize_site(std::string site) {
    for (char &base : site) {
        base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
        if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
            throw std::runtime_error("Long-read simulation requires an A/C/G/T enzyme motif.");
        }
    }
    if (site.empty()) {
        throw std::runtime_error("Enzyme site cannot be empty.");
    }
    return site;
}

std::size_t infer_cut_offset(const std::string &site) {
    if (site == "AAGCTT") {
        return 1;
    }
    if (site == "GATC") {
        return 0;
    }
    return site.size() / 2;
}

std::string build_ligation_junction(const std::string &site, std::size_t cut_offset) {
    return site.substr(0, site.size() - cut_offset) + site.substr(cut_offset);
}

char complement(char base) {
    switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

std::string reverse_complement(const std::string &sequence) {
    std::string result;
    result.reserve(sequence.size());
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        result.push_back(complement(*it));
    }
    return result;
}

std::vector<RestrictionFragment> digest_contig(const Contig &contig,
                                               std::size_t contig_index,
                                               const std::string &site,
                                               std::size_t cut_offset,
                                               std::size_t &next_id) {
    std::vector<RestrictionFragment> fragments;
    std::size_t start = 0;
    std::size_t search_from = 0;
    while (true) {
        const std::size_t position = contig.sequence.find(site, search_from);
        if (position == std::string::npos) {
            break;
        }
        const std::size_t cut = position + cut_offset;
        if (cut > start) {
            fragments.push_back(RestrictionFragment{next_id++, contig_index, start, cut});
        }
        start = cut;
        search_from = position + 1;
    }
    if (start < contig.sequence.size()) {
        fragments.push_back(
            RestrictionFragment{next_id++, contig_index, start, contig.sequence.size()});
    }
    if (fragments.empty()) {
        fragments.push_back(
            RestrictionFragment{next_id++, contig_index, 0, contig.sequence.size()});
    }
    return fragments;
}

FragmentIndex build_fragment_index(const Config &cfg,
                                   const ReferenceGenome &reference,
                                   const std::vector<OffsetEntry> &offsets,
                                   std::size_t matrix_bin_count) {
    std::unordered_map<std::string, std::size_t> contig_lookup;
    for (std::size_t i = 0; i < reference.contigs.size(); ++i) {
        contig_lookup.emplace(reference.contigs[i].name, i);
    }

    FragmentIndex index;
    index.bin_to_fragments.resize(matrix_bin_count);
    const std::string site = normalize_site(cfg.enzyme_site);
    const std::size_t cut_offset = infer_cut_offset(site);
    std::size_t next_id = 0;

    for (const auto &offset : offsets) {
        const auto contig_it = contig_lookup.find(offset.contig);
        if (contig_it == contig_lookup.end()) {
            throw std::runtime_error("Offset contig not found in FASTA: " + offset.contig);
        }
        const std::size_t contig_index = contig_it->second;
        const auto &contig = reference.contigs[contig_index];
        const std::size_t expected_bins =
            (contig.sequence.size() + cfg.bin_size - 1) / cfg.bin_size;
        if (offset.end_bin - offset.start_bin != expected_bins) {
            throw std::runtime_error("Offset bins do not match contig length for " + offset.contig);
        }

        auto fragments = digest_contig(contig, contig_index, site, cut_offset, next_id);
        for (const auto &fragment : fragments) {
            const std::size_t first_local_bin = fragment.start / cfg.bin_size;
            const std::size_t last_local_bin = (fragment.end - 1) / cfg.bin_size;
            for (std::size_t local_bin = first_local_bin; local_bin <= last_local_bin; ++local_bin) {
                const std::size_t global_bin = offset.start_bin + local_bin;
                if (global_bin >= index.bin_to_fragments.size()) {
                    throw std::runtime_error("Restriction fragment bin exceeds matrix bin count.");
                }
                index.bin_to_fragments[global_bin].push_back(fragment.id);
            }
            index.fragments.push_back(fragment);
        }
    }

    if (index.fragments.empty()) {
        throw std::runtime_error("Restriction digest produced no fragments.");
    }
    return index;
}

ContactSamplingIndex build_contact_index(const ContactMatrix &matrix,
                                         const FragmentIndex &fragments) {
    ContactSamplingIndex index;
    index.neighbors.resize(matrix.bin_count);
    for (const auto &contact : matrix.contacts) {
        const bool valid = contact.bin1 < fragments.bin_to_fragments.size() &&
                           contact.bin2 < fragments.bin_to_fragments.size() &&
                           !fragments.bin_to_fragments[contact.bin1].empty() &&
                           !fragments.bin_to_fragments[contact.bin2].empty() &&
                           std::isfinite(contact.weight) && contact.weight > 0.0;
        if (!valid) {
            continue;
        }
        index.contacts.push_back(contact);
        index.total_contact_weight += contact.weight;
        index.cumulative_contact_weights.push_back(index.total_contact_weight);

        auto &forward = index.neighbors[contact.bin1];
        const double forward_total = forward.empty() ? 0.0 : forward.back().cumulative_weight;
        forward.push_back(WeightedNeighbor{contact.bin2, forward_total + contact.weight});
        if (contact.bin1 != contact.bin2) {
            auto &reverse = index.neighbors[contact.bin2];
            const double reverse_total = reverse.empty() ? 0.0 : reverse.back().cumulative_weight;
            reverse.push_back(WeightedNeighbor{contact.bin1, reverse_total + contact.weight});
        }
    }
    if (index.contacts.empty()) {
        throw std::runtime_error("Contact matrix does not overlap restriction fragments.");
    }
    if (!std::isfinite(index.total_contact_weight) || index.total_contact_weight <= 0.0) {
        throw std::runtime_error("Usable contact matrix weights overflowed or summed to zero.");
    }
    return index;
}

template <typename Item, typename WeightAccessor>
std::size_t sample_cumulative_index(const std::vector<Item> &items,
                                    double total,
                                    WeightAccessor cumulative_weight,
                                    std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> distribution(0.0, total);
    const double value = distribution(rng);
    const auto it = std::upper_bound(
        items.begin(), items.end(), value,
        [&](double needle, const Item &item) { return needle < cumulative_weight(item); });
    if (it == items.end()) {
        return items.size() - 1;
    }
    return static_cast<std::size_t>(std::distance(items.begin(), it));
}

const Contact &sample_contact(const ContactSamplingIndex &index, std::mt19937_64 &rng) {
    const std::size_t sampled = sample_cumulative_index(
        index.cumulative_contact_weights,
        index.total_contact_weight,
        [](double value) { return value; },
        rng);
    return index.contacts[sampled];
}

std::size_t sample_neighbor(const ContactSamplingIndex &index,
                            std::size_t bin,
                            std::mt19937_64 &rng) {
    if (bin >= index.neighbors.size() || index.neighbors[bin].empty()) {
        const Contact &fallback = sample_contact(index, rng);
        return std::uniform_int_distribution<int>(0, 1)(rng) == 0
                   ? fallback.bin1
                   : fallback.bin2;
    }
    const auto &neighbors = index.neighbors[bin];
    const std::size_t sampled = sample_cumulative_index(
        neighbors,
        neighbors.back().cumulative_weight,
        [](const WeightedNeighbor &neighbor) { return neighbor.cumulative_weight; },
        rng);
    return neighbors[sampled].bin;
}

const RestrictionFragment &sample_fragment(const FragmentIndex &index,
                                           std::size_t bin,
                                           std::size_t previous_fragment,
                                           std::mt19937_64 &rng) {
    if (bin >= index.bin_to_fragments.size() || index.bin_to_fragments[bin].empty()) {
        throw std::runtime_error("Sampled contact bin has no restriction fragments.");
    }
    const auto &candidates = index.bin_to_fragments[bin];
    std::uniform_int_distribution<std::size_t> distribution(0, candidates.size() - 1);
    std::size_t fragment_id = candidates[distribution(rng)];
    if (candidates.size() > 1) {
        for (int retry = 0; retry < 12 && fragment_id == previous_fragment; ++retry) {
            fragment_id = candidates[distribution(rng)];
        }
    }
    return index.fragments[fragment_id];
}

std::size_t sample_target_length(const LongReadProfile &profile, std::mt19937_64 &rng) {
    const double sigma = profile.length_log_sigma;
    const double mu = std::log(static_cast<double>(profile.mean_length)) - sigma * sigma / 2.0;
    std::lognormal_distribution<double> distribution(mu, sigma);
    const double sampled = distribution(rng);
    const double bounded = std::max(
        static_cast<double>(profile.min_length),
        std::min(sampled, static_cast<double>(profile.max_length)));
    return static_cast<std::size_t>(std::llround(bounded));
}

std::string fragment_sequence(const ReferenceGenome &reference,
                              const RestrictionFragment &fragment,
                              bool reverse,
                              std::size_t take,
                              TruthSegment &truth) {
    const auto &contig = reference.contigs[fragment.contig_index];
    take = std::min(take, fragment.end - fragment.start);
    truth.contig = contig.name;
    truth.strand = reverse ? '-' : '+';
    if (!reverse) {
        truth.start = fragment.start;
        truth.end = fragment.start + take;
        return contig.sequence.substr(truth.start, take);
    }
    truth.start = fragment.end - take;
    truth.end = fragment.end;
    return reverse_complement(contig.sequence.substr(truth.start, take));
}

ConcatemerTemplate build_concatemer(const Config &cfg,
                                    const LongReadProfile &profile,
                                    const ReferenceGenome &reference,
                                    const FragmentIndex &fragments,
                                    const ContactSamplingIndex &contacts,
                                    const std::string &junction,
                                    std::mt19937_64 &rng) {
    const std::size_t target_length = sample_target_length(profile, rng);
    const Contact &first_contact = sample_contact(contacts, rng);
    std::vector<std::size_t> bins{first_contact.bin1, first_contact.bin2};
    if (std::uniform_int_distribution<int>(0, 1)(rng) == 1) {
        std::swap(bins[0], bins[1]);
    }

    ConcatemerTemplate result;
    result.sequence.reserve(target_length);
    result.segments.reserve(std::min<std::size_t>(cfg.max_concatemer_segments, 32));
    std::size_t previous_fragment = std::numeric_limits<std::size_t>::max();
    std::size_t current_bin = bins.front();
    const std::size_t anchor_bin = bins.front();
    const std::size_t minimum_terminal_segment =
        std::min<std::size_t>(500, std::max<std::size_t>(50, target_length / 20));

    for (std::size_t segment_index = 0;
         segment_index < cfg.max_concatemer_segments && result.sequence.size() < target_length;
         ++segment_index) {
        if (segment_index < bins.size()) {
            current_bin = bins[segment_index];
        } else {
            // A pairwise matrix cannot identify a unique higher-order joint distribution.
            // Sampling around one anchor avoids amplifying a single trans step into an
            // artificial chromosome-to-chromosome random walk.
            current_bin = sample_neighbor(contacts, anchor_bin, rng);
        }

        const RestrictionFragment &fragment =
            sample_fragment(fragments, current_bin, previous_fragment, rng);
        previous_fragment = fragment.id;

        if (!result.segments.empty()) {
            const std::size_t remaining = target_length - result.sequence.size();
            if (remaining <= junction.size()) {
                break;
            }
            result.sequence += junction;
        }

        std::size_t available = target_length - result.sequence.size();
        if (segment_index == 0 && available > minimum_terminal_segment) {
            available -= minimum_terminal_segment;
        }
        const std::size_t take =
            std::min(fragment.end - fragment.start, std::max<std::size_t>(1, available));
        TruthSegment truth;
        truth.bin = current_bin;
        truth.template_start = result.sequence.size();
        const bool reverse = std::uniform_int_distribution<int>(0, 1)(rng) == 1;
        result.sequence += fragment_sequence(reference, fragment, reverse, take, truth);
        truth.template_end = result.sequence.size();
        result.segments.push_back(std::move(truth));

        if (result.segments.size() >= 2 && result.sequence.size() >= target_length) {
            break;
        }
    }

    if (result.segments.size() < 2) {
        throw std::runtime_error(
            "Could not construct a multi-fragment concatemer from the contact matrix.");
    }
    result.hit_segment_limit =
        result.segments.size() == cfg.max_concatemer_segments &&
        result.sequence.size() < target_length;
    return result;
}

char mutate_base(char base, std::mt19937_64 &rng) {
    const char upper = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    const char alphabet[] = {'A', 'C', 'G', 'T'};
    if (upper != 'A' && upper != 'C' && upper != 'G' && upper != 'T') {
        return alphabet[std::uniform_int_distribution<int>(0, 3)(rng)];
    }
    char replacement = upper;
    while (replacement == upper) {
        replacement = alphabet[std::uniform_int_distribution<int>(0, 3)(rng)];
    }
    return replacement;
}

char random_base(std::mt19937_64 &rng) {
    static const char alphabet[] = {'A', 'C', 'G', 'T'};
    return alphabet[std::uniform_int_distribution<int>(0, 3)(rng)];
}

SimulatedSequence add_platform_errors(const std::string &input,
                                      const LongReadProfile &profile,
                                      std::mt19937_64 &rng) {
    SimulatedSequence result;
    result.sequence.reserve(input.size() + input.size() / 20);
    result.qualities.reserve(result.sequence.capacity());
    std::normal_distribution<double> quality_distribution(profile.mean_qv, profile.qv_sd);
    std::uniform_real_distribution<double> unit(0.0, 1.0);

    for (char base : input) {
        const int q = std::max(2, std::min(41,
            static_cast<int>(std::lround(quality_distribution(rng)))));
        const double error_probability = std::pow(10.0, -static_cast<double>(q) / 10.0);
        const double event = unit(rng);
        if (event < error_probability * profile.deletion_fraction) {
            continue;
        }

        char emitted = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
        if (event < error_probability *
                        (profile.deletion_fraction + profile.substitution_fraction)) {
            emitted = mutate_base(emitted, rng);
        } else if (emitted != 'A' && emitted != 'C' && emitted != 'G' && emitted != 'T' &&
                   emitted != 'N') {
            emitted = 'N';
        }
        result.sequence.push_back(emitted);
        result.qualities.push_back(static_cast<char>(q + 33));

        if (unit(rng) < error_probability * profile.insertion_fraction) {
            result.sequence.push_back(random_base(rng));
            result.qualities.push_back(static_cast<char>(q + 33));
        }
    }
    return result;
}

void append_truth_rows(std::string &output,
                       const std::string &read_name,
                       const std::vector<TruthSegment> &segments) {
    for (std::size_t i = 0; i < segments.size(); ++i) {
        const auto &segment = segments[i];
        output += read_name + '\t' + std::to_string(i + 1) + '\t' + segment.contig + '\t' +
                  std::to_string(segment.start) + '\t' + std::to_string(segment.end) + '\t' +
                  segment.strand + '\t' + std::to_string(segment.bin) + '\t' +
                  std::to_string(segment.template_start) + '\t' +
                  std::to_string(segment.template_end) + '\n';
    }
}

LongReadBlock make_block(const Config &cfg,
                         const LongReadProfile &profile,
                         const ReferenceGenome &reference,
                         const FragmentIndex &fragments,
                         const ContactSamplingIndex &contacts,
                         const std::string &junction,
                         std::size_t block_index,
                         std::size_t begin_read,
                         std::size_t end_read) {
    LongReadBlock block;
    block.read_count = end_read - begin_read;
    block.fastq.reserve(block.read_count * (profile.mean_length * 2 + 128));
    std::mt19937_64 rng(cfg.seed + 130363ULL * (block_index + 1));
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::vector<ConcatemerTemplate> templates;
    templates.reserve(block.read_count);

    for (std::size_t read_index = begin_read; read_index < end_read; ++read_index) {
        const bool duplicate_template = !templates.empty() &&
            unit(rng) < profile.template_duplicate_rate;
        std::size_t duplicate_source = 0;
        ConcatemerTemplate read_template;
        if (duplicate_template) {
            duplicate_source = std::uniform_int_distribution<std::size_t>(
                0, templates.size() - 1)(rng);
            read_template = templates[duplicate_source];
        } else {
            read_template =
                build_concatemer(cfg, profile, reference, fragments, contacts, junction, rng);
        }
        templates.push_back(read_template);
        SimulatedSequence read = add_platform_errors(read_template.sequence, profile, rng);
        std::string name = "hicreate_" + profile.assay_name + '_' +
                           std::to_string(read_index + 1) + "_segments" +
                           std::to_string(read_template.segments.size());
        if (duplicate_template) {
            name += "_dup" + std::to_string(begin_read + duplicate_source + 1);
        }
        block.fastq += '@' + name + '\n' + read.sequence + "\n+\n" + read.qualities + '\n';
        append_truth_rows(block.truth, name, read_template.segments);
        block.sequenced_bases += read.sequence.size();
        block.concatemer_segments += read_template.segments.size();
        if (read_template.hit_segment_limit) {
            ++block.segment_limited_reads;
        }
    }
    return block;
}

std::size_t requested_read_count(const Config &cfg,
                                 const ReferenceGenome &reference,
                                 const LongReadProfile &profile) {
    if (cfg.coverage_depth <= 0.0) {
        return cfg.pair_count;
    }
    const long double requested_bases =
        static_cast<long double>(cfg.coverage_depth) * reference.total_length();
    const long double reads =
        std::ceil(requested_bases / static_cast<long double>(profile.mean_length));
    if (!std::isfinite(static_cast<double>(reads)) ||
        reads > static_cast<long double>(std::numeric_limits<std::size_t>::max())) {
        throw std::runtime_error("Requested long-read coverage is too large.");
    }
    return std::max<std::size_t>(1, static_cast<std::size_t>(reads));
}

}  // namespace

LongReadProfile apply_long_read_overrides(const Config &cfg, LongReadProfile profile) {
    if (cfg.long_read_mean > 0) {
        profile.mean_length = cfg.long_read_mean;
    }
    if (cfg.long_read_min > 0) {
        profile.min_length = cfg.long_read_min;
    }
    if (cfg.long_read_max > 0) {
        profile.max_length = cfg.long_read_max;
    }
    if (cfg.long_read_sigma_explicit) {
        profile.length_log_sigma = cfg.long_read_sigma;
    }
    if (cfg.long_read_qv_explicit) {
        profile.mean_qv = cfg.long_read_qv;
    }
    if (profile.min_length > profile.mean_length || profile.mean_length > profile.max_length) {
        throw std::runtime_error(
            "Long-read lengths must satisfy min <= mean <= max after applying defaults.");
    }
    return profile;
}

LongReadSimulationResult write_long_contact_reads(
    const Config &cfg,
    const ReferenceGenome &reference,
    const std::vector<OffsetEntry> &offsets,
    const ContactMatrix &matrix,
    const LongReadProfile &profile) {
    const std::string site = normalize_site(cfg.enzyme_site);
    const std::string junction = build_ligation_junction(site, infer_cut_offset(site));
    const FragmentIndex fragments =
        build_fragment_index(cfg, reference, offsets, matrix.bin_count);
    const ContactSamplingIndex contacts = build_contact_index(matrix, fragments);
    const std::size_t read_count = requested_read_count(cfg, reference, profile);
    constexpr std::size_t target_block_bytes = 8 * 1024 * 1024;
    const std::size_t estimated_record_bytes = profile.mean_length * 2 + 256;
    const std::size_t reads_per_block = std::max<std::size_t>(
        1, std::min<std::size_t>(128, target_block_bytes / estimated_record_bytes));
    const std::size_t total_blocks = (read_count - 1) / reads_per_block + 1;

    LongReadSimulationResult result;
    result.fastq_path = cfg.output_prefix + '_' + profile.assay_name + ".fastq";
    result.truth_path = cfg.output_prefix + '_' + profile.assay_name + "_truth.tsv";
    std::ofstream fastq_out(result.fastq_path);
    std::ofstream truth_out(result.truth_path);
    if (!fastq_out || !truth_out) {
        throw std::runtime_error("Failed to open long-read FASTQ or truth output.");
    }
    truth_out << "read_name\tsegment_order\tcontig\tstart\tend\tstrand\tmatrix_bin"
                 "\ttemplate_start\ttemplate_end\n";

    const std::size_t worker_count =
        std::max<std::size_t>(1, std::min(cfg.thread_count, total_blocks));
    std::cerr << "Restriction digest for " << profile.assay_name << ": "
              << fragments.fragments.size() << " fragments; "
              << contacts.contacts.size() << " usable matrix contacts\n"
              << "Streaming " << read_count << ' ' << profile.assay_name
              << " long reads with " << worker_count << " worker thread(s)...\n";

    const auto write_block = [&](const LongReadBlock &block) {
        fastq_out << block.fastq;
        truth_out << block.truth;
        if (!fastq_out || !truth_out) {
            throw std::runtime_error("Failed while writing long-read output.");
        }
        result.read_count += block.read_count;
        result.sequenced_bases += block.sequenced_bases;
        result.concatemer_segments += block.concatemer_segments;
        result.segment_limited_reads += block.segment_limited_reads;
    };

    if (worker_count == 1 || total_blocks < 2) {
        for (std::size_t block_index = 0; block_index < total_blocks; ++block_index) {
            const std::size_t begin = block_index * reads_per_block;
            const std::size_t end = std::min(read_count, begin + reads_per_block);
            write_block(make_block(cfg, profile, reference, fragments, contacts, junction,
                                   block_index, begin, end));
        }
    } else {
        const std::size_t max_ready_blocks = worker_count * 2 + 2;
        std::atomic<std::size_t> next_block{0};
        std::atomic<std::size_t> finished_workers{0};
        std::atomic<bool> stop_workers{false};
        std::mutex mutex;
        std::condition_variable cv;
        std::map<std::size_t, LongReadBlock> ready_blocks;
        std::exception_ptr worker_error;
        std::exception_ptr writer_error;
        std::size_t next_to_write = 0;

        const auto worker = [&]() {
            try {
                while (!stop_workers.load()) {
                    const std::size_t block_index = next_block.fetch_add(1);
                    if (block_index >= total_blocks) {
                        break;
                    }
                    const std::size_t begin = block_index * reads_per_block;
                    const std::size_t end = std::min(read_count, begin + reads_per_block);
                    LongReadBlock block = make_block(cfg, profile, reference, fragments, contacts,
                                                     junction, block_index, begin, end);
                    {
                        std::unique_lock<std::mutex> lock(mutex);
                        cv.wait(lock, [&]() {
                            return stop_workers.load() ||
                                   ready_blocks.size() < max_ready_blocks ||
                                   block_index == next_to_write;
                        });
                        if (stop_workers.load()) {
                            break;
                        }
                        ready_blocks.emplace(block_index, std::move(block));
                    }
                    cv.notify_one();
                }
            } catch (...) {
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    if (!worker_error) {
                        worker_error = std::current_exception();
                    }
                    stop_workers.store(true);
                }
                cv.notify_all();
            }
            finished_workers.fetch_add(1);
            cv.notify_all();
        };

        std::vector<std::thread> workers;
        workers.reserve(worker_count);
        for (std::size_t i = 0; i < worker_count; ++i) {
            workers.emplace_back(worker);
        }

        try {
            while (next_to_write < total_blocks) {
                LongReadBlock block;
                {
                    std::unique_lock<std::mutex> lock(mutex);
                    cv.wait(lock, [&]() {
                        return ready_blocks.find(next_to_write) != ready_blocks.end() ||
                               worker_error || finished_workers.load() == worker_count;
                    });
                    auto it = ready_blocks.find(next_to_write);
                    if (it == ready_blocks.end()) {
                        break;
                    }
                    block = std::move(it->second);
                    ready_blocks.erase(it);
                }
                cv.notify_all();
                write_block(block);
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    ++next_to_write;
                }
                cv.notify_all();
            }
        } catch (...) {
            writer_error = std::current_exception();
        }

        stop_workers.store(true);
        cv.notify_all();
        for (auto &thread : workers) {
            if (thread.joinable()) {
                thread.join();
            }
        }
        if (worker_error) {
            std::rethrow_exception(worker_error);
        }
        if (writer_error) {
            std::rethrow_exception(writer_error);
        }
        if (next_to_write != total_blocks) {
            throw std::runtime_error("Long-read generation stopped before all blocks were written.");
        }
    }

    fastq_out.flush();
    truth_out.flush();
    if (!fastq_out || !truth_out) {
        throw std::runtime_error("Failed to finalize long-read output.");
    }
    if (result.segment_limited_reads > 0) {
        std::cerr << "Warning: " << result.segment_limited_reads
                  << " long reads reached --max-segments before their sampled target length; "
                     "increase --max-segments if this is not intended.\n";
    }
    return result;
}
