#include "fragmenter.h"
#include <algorithm>
#include <cctype>
#include <random>
#include <stdexcept>
#include <unordered_map>

namespace {

struct RestrictionFragment {
    std::size_t id = 0;
    std::size_t contig_index = 0;
    std::size_t start = 0;
    std::size_t end = 0;
    std::string sequence;
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
                next_id++, contig_index, start, cut, contig.sequence.substr(start, cut - start)});
        }
        start = cut;
        pos = pos + 1;
    }

    if (start < contig.sequence.size()) {
        fragments.push_back(RestrictionFragment{
            next_id++, contig_index, start, contig.sequence.size(),
            contig.sequence.substr(start, contig.sequence.size() - start)});
    }

    if (fragments.empty()) {
        fragments.push_back(RestrictionFragment{next_id++, contig_index, 0, contig.sequence.size(), contig.sequence});
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

}  // namespace

std::vector<LigationProduct> create_ligation_products(const Config &cfg,
                                                      const ReferenceGenome &reference,
                                                      const std::vector<OffsetEntry> &offsets,
                                                      const ContactMatrix &matrix) {
    std::unordered_map<std::string, std::size_t> contig_lookup;
    for (std::size_t i = 0; i < reference.contigs.size(); ++i) {
        contig_lookup.emplace(reference.contigs[i].name, i);
    }

    std::vector<RestrictionFragment> fragments;
    std::vector<std::vector<std::size_t>> bin_to_fragments(matrix.bin_count);
    std::size_t next_fragment_id = 0;
    const std::size_t cut_offset = infer_cut_offset(cfg.enzyme_site);
    const std::string ligation_junction = build_ligation_junction(cfg.enzyme_site, cut_offset);

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
                if (global_bin >= bin_to_fragments.size()) {
                    throw std::runtime_error("Fragment bin exceeds matrix bin count.");
                }
                bin_to_fragments[global_bin].push_back(fragment.id);
            }
            fragments.push_back(fragment);
        }
    }

    if (fragments.empty()) {
        throw std::runtime_error("Restriction digest produced no fragments.");
    }

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

    std::discrete_distribution<std::size_t> contact_dist(weights.begin(), weights.end());
    std::mt19937_64 rng(cfg.seed + 1);
    std::vector<LigationProduct> products;
    products.reserve(cfg.pair_count);

    for (std::size_t i = 0; i < cfg.pair_count; ++i) {
        const Contact &contact = matrix.contacts[contact_dist(rng)];
        const auto &left_candidates = bin_to_fragments[contact.bin1];
        const auto &right_candidates = bin_to_fragments[contact.bin2];
        std::uniform_int_distribution<std::size_t> left_dist(0, left_candidates.size() - 1);
        std::uniform_int_distribution<std::size_t> right_dist(0, right_candidates.size() - 1);

        const RestrictionFragment &left = fragments[left_candidates[left_dist(rng)]];
        std::size_t right_index = right_candidates[right_dist(rng)];
        if (right_candidates.size() > 1 || left_candidates.size() > 1) {
            for (int retry = 0; retry < 8 && right_index == left.id; ++retry) {
                right_index = right_candidates[right_dist(rng)];
            }
        }
        const RestrictionFragment &right = fragments[right_index];

        std::string left_seq = left.sequence;
        std::string right_seq = right.sequence;
        if (std::uniform_int_distribution<int>(0, 1)(rng) == 1) {
            left_seq = reverse_complement(left_seq);
        }
        if (std::uniform_int_distribution<int>(0, 1)(rng) == 1) {
            right_seq = reverse_complement(right_seq);
        }

        LigationProduct product;
        product.name = "ligation_" + std::to_string(i + 1) +
                       "_bin" + std::to_string(contact.bin1) +
                       "_bin" + std::to_string(contact.bin2);
        product.sequence = left_seq + ligation_junction + right_seq;
        products.push_back(std::move(product));
    }

    return products;
}
