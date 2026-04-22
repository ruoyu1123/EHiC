#pragma once

#include "config.h"
#include "fragmenter.h"
#include <vector>

void simulate_paired_reads(const Config &cfg, const std::vector<ReadPairTemplate> &read_templates);
