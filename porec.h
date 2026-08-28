#pragma once

#include "long_read_common.h"

LongReadProfile porec_profile(const Config &cfg);

LongReadSimulationResult write_porec_reads(
    const Config &cfg,
    const ReferenceGenome &reference,
    const std::vector<OffsetEntry> &offsets,
    const ContactMatrix &matrix);
