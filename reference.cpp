#include "reference.h"
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

std::string normalize_contig_name(const std::string &header_text) {
    std::istringstream stream(header_text);
    std::string name;
    stream >> name;
    return name;
}

}

std::size_t ReferenceGenome::total_length() const {
    std::size_t total = 0;
    for (const auto &contig : contigs) {
        total += contig.sequence.size();
    }
    return total;
}

std::size_t ReferenceGenome::total_bins(std::size_t bin_size) const {
    if (bin_size == 0) {
        return 0;
    }

    std::size_t total = 0;
    for (const auto &contig : contigs) {
        total += (contig.sequence.size() + bin_size - 1) / bin_size;
    }
    return total;
}

ReferenceGenome load_reference_fasta(const std::string &path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open reference FASTA: " + path);
    }

    ReferenceGenome genome;
    std::string line;
    Contig current;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '>') {
            if (!current.name.empty()) {
                genome.contigs.push_back(current);
                current = Contig{};
            }
            current.name = normalize_contig_name(line.substr(1));
            if (current.name.empty()) {
                throw std::runtime_error("FASTA header cannot be empty.");
            }
            continue;
        }

        if (current.name.empty()) {
            current.name = "contig_1";
        }

        for (char ch : line) {
            const char upper = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            if (!std::isspace(static_cast<unsigned char>(upper))) {
                current.sequence.push_back(upper);
            }
        }
    }

    if (!current.name.empty()) {
        genome.contigs.push_back(current);
    }

    if (genome.contigs.empty()) {
        throw std::runtime_error("Reference FASTA does not contain sequence.");
    }

    for (const auto &contig : genome.contigs) {
        if (contig.sequence.empty()) {
            throw std::runtime_error("Reference FASTA contains an empty contig: " + contig.name);
        }
    }

    return genome;
}

std::vector<OffsetEntry> build_reference_offsets(const ReferenceGenome &reference, std::size_t bin_size) {
    if (bin_size == 0) {
        throw std::runtime_error("Bin size must be positive when building reference offsets.");
    }

    std::vector<OffsetEntry> offsets;
    offsets.reserve(reference.contigs.size());

    std::size_t start_bin = 0;
    for (const auto &contig : reference.contigs) {
        const std::size_t contig_bins = (contig.sequence.size() + bin_size - 1) / bin_size;
        offsets.push_back(OffsetEntry{contig.name, start_bin, start_bin + contig_bins});
        start_bin += contig_bins;
    }

    return offsets;
}
