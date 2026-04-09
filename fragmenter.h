#pragma once

#include "config.h"
#include "matrix.h"
#include "reference.h"
#include <string>
#include <vector>

struct LigationProduct {
    std::string name;
    std::string sequence;
};

std::vector<LigationProduct> create_ligation_products(const Config &cfg,
                                                      const ReferenceGenome &reference,
                                                      const std::vector<OffsetEntry> &offsets,
                                                      const ContactMatrix &matrix);
