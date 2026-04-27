#include "fragmenter.h"
#include "simulator.h"
#include <algorithm>
#include <atomic>
#include <cctype>
#include <condition_variable>
#include <exception>
#include <cmath>
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

std::size_t infer_cut_offset(const std::string &site) {
    if (site == "AAGCTT") {
        return 1;  // HindIII: A|AGCTT
    }
    if (site == "GATC") {
        return 0;  // DpnII/MboI: ^GATC
    }
    return site.size() / 2;
}

std::string build_ligation_junction(const std::string &site, std::size_t cut_offset) {
    if (site.empty() || cut_offset > site.size()) {
        throw std::runtime_error("Invalid enzyme site or cut offset.");
    }
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
    std::string rc;
    rc.reserve(sequence.size());
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        rc.push_back(complement(*it));
    }
    return rc;
}

std::string pad_to_length(std::string sequence, std::size_t length) {
    if (sequence.size() > length) {
        sequence.resize(length);
    }
    if (sequence.size() < length) {
        sequence.append(length - sequence.size(), 'N');
    }
    return sequence;
}

std::string reference_slice(const ReferenceGenome &reference,
                            const RestrictionFragment &fragment,
                            std::size_t start,
                            std::size_t end) {
    const auto &sequence = reference.contigs[fragment.contig_index].sequence;
    if (start > end || start < fragment.start || end > fragment.end || end > sequence.size()) {
        throw std::runtime_error("Invalid restriction fragment slice.");
    }
    return sequence.substr(start, end - start);
}

void append_prefix(std::string &out, const std::string &sequence, std::size_t target_length) {
    if (out.size() >= target_length) {
        return;
    }
    const std::size_t remaining = target_length - out.size();
    out.append(sequence, 0, std::min(remaining, sequence.size()));
}

struct FragmentEndTemplate {
    std::string outward;
};

FragmentEndTemplate get_fragment_end_template(const ReferenceGenome &reference,
                                              const RestrictionFragment &fragment,
                                              bool use_right_end,
                                              std::size_t read_length) {
    const std::size_t fragment_length = fragment.end - fragment.start;
    const std::size_t take = std::min(fragment_length, read_length);
    if (!use_right_end) {
        return FragmentEndTemplate{
            reference_slice(reference, fragment, fragment.start, fragment.start + take)
        };
    }

    return FragmentEndTemplate{
        reverse_complement(reference_slice(reference, fragment, fragment.end - take, fragment.end))
    };
}

std::vector<RestrictionFragment> digest_contig(const Contig &contig,
                                               std::size_t contig_index,
                                               const std::string &site,
                                               std::size_t cut_offset,
                                               std::size_t &next_id) {
    if (site.empty()) {
        throw std::runtime_error("Enzyme site cannot be empty.");
    }

    std::vector<RestrictionFragment> fragments;
    std::size_t start = 0;
    std::size_t pos = 0;
    while (true) {
        pos = contig.sequence.find(site, pos);
        if (pos == std::string::npos) {
            break;
        }

        const std::size_t cut = pos + cut_offset;
        if (cut > start) {
            fragments.push_back(RestrictionFragment{
                next_id++, contig_index, start, cut});
        }
        start = cut;
        pos = pos + 1;
    }

    if (start < contig.sequence.size()) {
        fragments.push_back(RestrictionFragment{
            next_id++, contig_index, start, contig.sequence.size()});
    }

    if (fragments.empty()) {
        fragments.push_back(RestrictionFragment{next_id++, contig_index, 0, contig.sequence.size()});
    }

    return fragments;
}

std::vector<std::size_t> bins_for_fragment(const RestrictionFragment &fragment, std::size_t bin_size) {
    const std::size_t start_bin = fragment.start / bin_size;
    const std::size_t end_bin = (fragment.end - 1) / bin_size;
    std::vector<std::size_t> bins;
    bins.reserve(end_bin - start_bin + 1);
    for (std::size_t bin = start_bin; bin <= end_bin; ++bin) {
        bins.push_back(bin);
    }
    return bins;
}

struct FragmentSamplingIndex {
    std::vector<RestrictionFragment> fragments;
    std::vector<std::vector<std::size_t>> bin_to_fragments;
};

FragmentSamplingIndex build_fragment_sampling_index(const Config &cfg,
                                                    const ReferenceGenome &reference,
                                                    const std::vector<OffsetEntry> &offsets,
                                                    std::size_t matrix_bin_count) {
    std::unordered_map<std::string, std::size_t> contig_lookup;
    for (std::size_t i = 0; i < reference.contigs.size(); ++i) {
        contig_lookup.emplace(reference.contigs[i].name, i);
    }

    FragmentSamplingIndex index;
    index.bin_to_fragments.resize(matrix_bin_count);

    std::size_t next_fragment_id = 0;
    const std::size_t cut_offset = infer_cut_offset(cfg.enzyme_site);
    for (const auto &offset : offsets) {
        const auto it = contig_lookup.find(offset.contig);
        if (it == contig_lookup.end()) {
            throw std::runtime_error("Offset contig not found in FASTA: " + offset.contig);
        }

        const std::size_t contig_index = it->second;
        const std::size_t contig_bins =
            (reference.contigs[contig_index].sequence.size() + cfg.bin_size - 1) / cfg.bin_size;
        if (offset.end_bin - offset.start_bin != contig_bins) {
            throw std::runtime_error("Offset bins do not match contig length for " + offset.contig);
        }

        auto contig_fragments =
            digest_contig(reference.contigs[contig_index], contig_index, cfg.enzyme_site, cut_offset, next_fragment_id);
        for (const auto &fragment : contig_fragments) {
            const auto local_bins = bins_for_fragment(fragment, cfg.bin_size);
            for (const std::size_t local_bin : local_bins) {
                const std::size_t global_bin = offset.start_bin + local_bin;
                if (global_bin >= index.bin_to_fragments.size()) {
                    throw std::runtime_error("Fragment bin exceeds matrix bin count.");
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

std::vector<double> build_contact_sampling_weights(const ContactMatrix &matrix,
                                                   const std::vector<std::vector<std::size_t>> &bin_to_fragments) {
    std::vector<double> weights;
    weights.reserve(matrix.contacts.size());
    double total_weight = 0.0;
    for (const auto &contact : matrix.contacts) {
        const bool valid_bin1 = contact.bin1 < bin_to_fragments.size() && !bin_to_fragments[contact.bin1].empty();
        const bool valid_bin2 = contact.bin2 < bin_to_fragments.size() && !bin_to_fragments[contact.bin2].empty();
        const double weight = valid_bin1 && valid_bin2 ? contact.weight : 0.0;
        weights.push_back(weight);
        total_weight += weight;
    }

    if (total_weight <= 0.0) {
        throw std::runtime_error("Contact matrix does not overlap any valid reference bins.");
    }

    return weights;
}

struct CumulativeWeights {
    std::vector<double> cumulative;
    double total = 0.0;
};

CumulativeWeights build_cumulative_weights(const std::vector<double> &weights) {
    CumulativeWeights result;
    result.cumulative.reserve(weights.size());
    for (const double weight : weights) {
        if (std::isfinite(weight) && weight > 0.0) {
            result.total += weight;
        }
        result.cumulative.push_back(result.total);
    }
    if (result.total <= 0.0) {
        throw std::runtime_error("Contact matrix does not contain positive usable weights.");
    }
    return result;
}

std::size_t sample_contact_index(const CumulativeWeights &weights, std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> dist(0.0, weights.total);
    const double value = dist(rng);
    const auto it = std::upper_bound(weights.cumulative.begin(), weights.cumulative.end(), value);
    if (it == weights.cumulative.end()) {
        return weights.cumulative.size() - 1;
    }
    return static_cast<std::size_t>(std::distance(weights.cumulative.begin(), it));
}

ReadPairTemplate make_read_template(const std::string &name,
                                    const FragmentEndTemplate &left,
                                    const std::string &ligation_junction,
                                    const FragmentEndTemplate &right,
                                    std::size_t read_length) {
    std::string read1_template;
    read1_template.reserve(read_length);
    append_prefix(read1_template, left.outward, read_length);
    append_prefix(read1_template, ligation_junction, read_length);
    append_prefix(read1_template, right.outward, read_length);

    std::string read2_template;
    read2_template.reserve(read_length);
    append_prefix(read2_template, right.outward, read_length);
    append_prefix(read2_template, ligation_junction, read_length);
    append_prefix(read2_template, left.outward, read_length);

    return ReadPairTemplate{
        name,
        pad_to_length(read1_template, read_length),
        pad_to_length(read2_template, read_length)
    };
}

ReadPairTemplate sample_one_read_template(const Config &cfg,
                                          const ReferenceGenome &reference,
                                          const FragmentSamplingIndex &index,
                                          const ContactMatrix &matrix,
                                          const std::string &ligation_junction,
                                          const CumulativeWeights &contact_weights,
                                          std::mt19937_64 &rng,
                                          std::size_t pair_index) {
    const Contact &contact = matrix.contacts[sample_contact_index(contact_weights, rng)];
    const auto &left_candidates = index.bin_to_fragments[contact.bin1];
    const auto &right_candidates = index.bin_to_fragments[contact.bin2];
    std::uniform_int_distribution<std::size_t> left_dist(0, left_candidates.size() - 1);
    std::uniform_int_distribution<std::size_t> right_dist(0, right_candidates.size() - 1);

    const RestrictionFragment &left = index.fragments[left_candidates[left_dist(rng)]];
    std::size_t right_index = right_candidates[right_dist(rng)];
    if (right_candidates.size() > 1 || left_candidates.size() > 1) {
        for (int retry = 0; retry < 8 && right_index == left.id; ++retry) {
            right_index = right_candidates[right_dist(rng)];
        }
    }
    const RestrictionFragment &right = index.fragments[right_index];

    const bool use_left_right_end = std::uniform_int_distribution<int>(0, 1)(rng) == 1;
    const bool use_right_right_end = std::uniform_int_distribution<int>(0, 1)(rng) == 1;
    const FragmentEndTemplate left_ends =
        get_fragment_end_template(reference, left, use_left_right_end, cfg.read_length);
    const FragmentEndTemplate right_ends =
        get_fragment_end_template(reference, right, use_right_right_end, cfg.read_length);

    const std::string name = "hicreate_" + std::to_string(pair_index + 1) +
                             "_bin" + std::to_string(contact.bin1) +
                             "_bin" + std::to_string(contact.bin2);
    return make_read_template(name, left_ends, ligation_junction, right_ends, cfg.read_length);
}

template <typename EmitReadPair>
void sample_hic_read_templates(const Config &cfg,
                               const ReferenceGenome &reference,
                               const std::vector<OffsetEntry> &offsets,
                               const ContactMatrix &matrix,
                               EmitReadPair emit) {
    const auto index = build_fragment_sampling_index(cfg, reference, offsets, matrix.bin_count);
    const std::size_t cut_offset = infer_cut_offset(cfg.enzyme_site);
    const std::string ligation_junction = build_ligation_junction(cfg.enzyme_site, cut_offset);
    const std::vector<double> weights = build_contact_sampling_weights(matrix, index.bin_to_fragments);
    const CumulativeWeights contact_weights = build_cumulative_weights(weights);

    std::mt19937_64 rng(cfg.seed + 1);

    for (std::size_t i = 0; i < cfg.pair_count; ++i) {
        emit(sample_one_read_template(cfg, reference, index, matrix, ligation_junction, contact_weights, rng, i));
    }
}

}  // namespace

void write_hic_reads(const Config &cfg,
                     const ReferenceGenome &reference,
                     const std::vector<OffsetEntry> &offsets,
                     const ContactMatrix &matrix,
                     PairedReadWriter &writer) {
    const auto index = build_fragment_sampling_index(cfg, reference, offsets, matrix.bin_count);
    const std::size_t cut_offset = infer_cut_offset(cfg.enzyme_site);
    const std::string ligation_junction = build_ligation_junction(cfg.enzyme_site, cut_offset);
    const std::vector<double> weights = build_contact_sampling_weights(matrix, index.bin_to_fragments);
    const CumulativeWeights contact_weights = build_cumulative_weights(weights);

    const std::size_t worker_count = std::max<std::size_t>(1, cfg.thread_count);
    if (worker_count == 1 || cfg.pair_count < 20000) {
        std::mt19937_64 rng(cfg.seed + 1);
        for (std::size_t i = 0; i < cfg.pair_count; ++i) {
            writer.write_template(
                sample_one_read_template(cfg, reference, index, matrix, ligation_junction, contact_weights, rng, i));
        }
        return;
    }

    constexpr std::size_t pairs_per_block = 8192;
    const std::size_t total_blocks = (cfg.pair_count + pairs_per_block - 1) / pairs_per_block;
    const std::size_t max_ready_blocks = worker_count * 2 + 2;
    std::atomic<std::size_t> next_block{0};
    std::atomic<std::size_t> finished_workers{0};
    std::atomic<bool> stop_workers{false};
    std::mutex mutex;
    std::condition_variable cv;
    std::map<std::size_t, FastqBlock> ready_blocks;
    std::exception_ptr worker_error;

    const auto make_block = [&](std::size_t block_index) {
        const std::size_t begin_pair = block_index * pairs_per_block;
        const std::size_t end_pair = std::min(cfg.pair_count, begin_pair + pairs_per_block);

        FastqBlock block;
        block.pair_count = end_pair - begin_pair;
        block.read1.reserve(block.pair_count * 360);
        block.read2.reserve(block.pair_count * 360);

        std::mt19937_64 rng(cfg.seed + 104729ULL * (block_index + 1));
        for (std::size_t pair_index = begin_pair; pair_index < end_pair; ++pair_index) {
            ReadPairTemplate read_pair =
                sample_one_read_template(cfg, reference, index, matrix, ligation_junction, contact_weights, rng, pair_index);
            append_simulated_fastq_pair(cfg, read_pair, rng, block.read1, block.read2);
        }
        return block;
    };

    const auto worker = [&]() {
        try {
            while (!stop_workers.load()) {
                const std::size_t block_index = next_block.fetch_add(1);
                if (block_index >= total_blocks) {
                    break;
                }
                FastqBlock block = make_block(block_index);
                {
                    std::unique_lock<std::mutex> lock(mutex);
                    cv.wait(lock, [&]() {
                        return stop_workers.load() || ready_blocks.size() < max_ready_blocks;
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

    std::size_t next_to_write = 0;
    while (next_to_write < total_blocks) {
        FastqBlock block;
        {
            std::unique_lock<std::mutex> lock(mutex);
            cv.wait(lock, [&]() {
                return ready_blocks.find(next_to_write) != ready_blocks.end() ||
                       worker_error ||
                       finished_workers.load() == worker_count;
            });

            auto it = ready_blocks.find(next_to_write);
            if (it == ready_blocks.end()) {
                break;
            }
            block = std::move(it->second);
            ready_blocks.erase(it);
        }
        cv.notify_all();
        writer.write_block(block);
        ++next_to_write;
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
    if (next_to_write != total_blocks) {
        throw std::runtime_error("Read generation stopped before all FASTQ blocks were written.");
    }
}
