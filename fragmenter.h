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

class PairedReadWriter;

void write_hic_reads(const Config &cfg,
                     const ReferenceGenome &reference,
                     const std::vector<OffsetEntry> &offsets,
                     const ContactMatrix &matrix,
                     PairedReadWriter &writer);
