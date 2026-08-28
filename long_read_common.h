#pragma once

#include "config.h"
#include "matrix.h"
#include "reference.h"

#include <cstddef>
#include <string>
#include <vector>

struct LongReadProfile {
    std::string assay_name;
    std::size_t mean_length = 0;
    std::size_t min_length = 0;
    std::size_t max_length = 0;
    double length_log_sigma = 0.0;
    double mean_qv = 0.0;
    double qv_sd = 0.0;
    double substitution_fraction = 0.0;
    double deletion_fraction = 0.0;
    double insertion_fraction = 0.0;
    double template_duplicate_rate = 0.0;
};

struct LongReadSimulationResult {
    std::string fastq_path;
    std::string truth_path;
    std::size_t read_count = 0;
    std::size_t sequenced_bases = 0;
    std::size_t concatemer_segments = 0;
    std::size_t segment_limited_reads = 0;
};

LongReadProfile apply_long_read_overrides(const Config &cfg, LongReadProfile profile);

LongReadSimulationResult write_long_contact_reads(
    const Config &cfg,
    const ReferenceGenome &reference,
    const std::vector<OffsetEntry> &offsets,
    const ContactMatrix &matrix,
    const LongReadProfile &profile);
