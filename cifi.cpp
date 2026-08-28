#include "cifi.h"

LongReadProfile cifi_profile(const Config &cfg) {
    LongReadProfile profile;
    profile.assay_name = "cifi";
    profile.mean_length = 9350;
    profile.min_length = 5000;
    profile.max_length = 25000;
    profile.length_log_sigma = 0.65;
    profile.mean_qv = 38.0;
    profile.qv_sd = 1.5;
    profile.substitution_fraction = 0.55;
    profile.deletion_fraction = 0.25;
    profile.insertion_fraction = 0.20;
    profile.template_duplicate_rate = 0.018;
    return apply_long_read_overrides(cfg, profile);
}

LongReadSimulationResult write_cifi_reads(
    const Config &cfg,
    const ReferenceGenome &reference,
    const std::vector<OffsetEntry> &offsets,
    const ContactMatrix &matrix) {
    return write_long_contact_reads(cfg, reference, offsets, matrix, cifi_profile(cfg));
}
