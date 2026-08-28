#pragma once

#include "long_read_common.h"

LongReadProfile cifi_profile(const Config &cfg);

LongReadSimulationResult write_cifi_reads(
    const Config &cfg,
    const ReferenceGenome &reference,
    const std::vector<OffsetEntry> &offsets,
    const ContactMatrix &matrix);
