#include "porec.h"

LongReadProfile porec_profile(const Config &cfg) {
    LongReadProfile profile;
    profile.assay_name = "porec";
    profile.mean_length = 10000;
    profile.min_length = 3000;
    profile.max_length = 100000;
    profile.length_log_sigma = 0.80;
    profile.mean_qv = 13.0;
    profile.qv_sd = 2.0;
    profile.substitution_fraction = 0.25;
    profile.deletion_fraction = 0.45;
    profile.insertion_fraction = 0.30;
    return apply_long_read_overrides(cfg, profile);
}

LongReadSimulationResult write_porec_reads(
    const Config &cfg,
    const ReferenceGenome &reference,
    const std::vector<OffsetEntry> &offsets,
    const ContactMatrix &matrix) {
    return write_long_contact_reads(cfg, reference, offsets, matrix, porec_profile(cfg));
}
