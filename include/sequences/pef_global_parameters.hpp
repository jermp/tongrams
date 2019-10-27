#pragma once

namespace pef {
struct pef_global_parameters {
    pef_global_parameters()
        : ef_log_sampling0(9)
        , ef_log_sampling1(8)
        , rb_log_rank1_sampling(9)
        , rb_log_sampling1(8) {}

    uint8_t ef_log_sampling0;
    uint8_t ef_log_sampling1;
    uint8_t rb_log_rank1_sampling;
    uint8_t rb_log_sampling1;
};
}  // namespace pef
