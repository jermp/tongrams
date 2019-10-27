#pragma once

#include "utils/util.hpp"
#include "vectors/bit_vector.hpp"

namespace pef {

// note: n can be 0
void write_gamma(tongrams::bit_vector_builder& bvb, uint64_t n) {
    uint64_t nn = n + 1;
    uint64_t l = tongrams::util::msb(nn);
    uint64_t hb = uint64_t(1) << l;
    bvb.append_bits(hb, l + 1);
    bvb.append_bits(nn ^ hb, l);
}

void write_gamma_nonzero(tongrams::bit_vector_builder& bvb, uint64_t n) {
    assert(n > 0);
    write_gamma(bvb, n - 1);
}

uint64_t read_gamma(tongrams::bits_iterator<tongrams::bit_vector>& it) {
    uint64_t l = it.skip_zeros();
    return (it.get_bits(l) | (uint64_t(1) << l)) - 1;
}

uint64_t read_gamma_nonzero(tongrams::bits_iterator<tongrams::bit_vector>& it) {
    return read_gamma(it) + 1;
}

void write_delta(tongrams::bit_vector_builder& bvb, uint64_t n) {
    uint64_t nn = n + 1;
    uint64_t l = tongrams::util::msb(nn);
    uint64_t hb = uint64_t(1) << l;
    write_gamma(bvb, l);
    bvb.append_bits(nn ^ hb, l);
}

uint64_t read_delta(tongrams::bits_iterator<tongrams::bit_vector>& it) {
    uint64_t l = read_gamma(it);
    return (it.get_bits(l) | (uint64_t(1) << l)) - 1;
}

}  // namespace pef
