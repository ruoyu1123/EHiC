#pragma once

#include "config.h"
#include "matrix.h"
#include "reference.h"
#include <string>
#include <vector>

struct ReadPairTemplate {
    std::string name;
    std::string read1;
    std::string read2;
};

std::vector<ReadPairTemplate> create_hic_read_templates(const Config &cfg,
                                                        const ReferenceGenome &reference,
                                                        const std::vector<OffsetEntry> &offsets,
                                                        const ContactMatrix &matrix);
