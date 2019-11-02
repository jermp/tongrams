#pragma once

#include <iostream>
#include <vector>
#include <fstream>
#include <ctime>
#include <sys/time.h>
#include <type_traits>
#include <cassert>
#include <locale>
#include <string.h>

#include <xmmintrin.h>
#if TONGRAMS_USE_POPCNT
#include <smmintrin.h>
#endif

#include "utils/binary_header.hpp"
#include "../external/essentials/include/essentials.hpp"

#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#define NOINLINE __attribute__((noinline))
#define ALWAYSINLINE __attribute__((always_inline))

namespace tongrams {

namespace bits {
// NOTE: for union punning
typedef union {
    float float_value;
    uint64_t uint64_value;
} reinterpret;

// pack 2 floats into 1 unit64_t
void pack(uint64_t& packed, float prob, float backoff) {
    reinterpret reint;
    reint.float_value = backoff;
    packed = reint.uint64_value << 32;
    reint.float_value = prob;
    packed |= reint.uint64_value;
}

inline void unpack(uint64_t packed, float& prob, float& backoff) {
    const static uint64_t mask = (uint64_t(1) << 32) - 1;
    reinterpret reint;
    reint.uint64_value = packed & mask;
    prob = reint.float_value;
    reint.uint64_value = packed >> 32;
    backoff = reint.float_value;
}
}  // namespace bits

namespace global {
static const uint8_t null_header = uint8_t(-1);

static const uint8_t max_order = 8;
static const uint8_t max_remapping_order = 2;
static const uint64_t not_found = uint64_t(-1);

static const float default_unk_prob = -100.0;
static const uint8_t default_probs_quantization_bits = 8;
static const uint8_t default_backoffs_quantization_bits = 8;
}  // namespace global

namespace style {
std::ostream& bold(std::ostream& os) {
    return os << "\e[1m";
}

std::ostream& underline(std::ostream& os) {
    return os << "\e[4m";
}

std::ostream& off(std::ostream& os) {
    return os << "\e[0m";
}
}  // namespace style

namespace bytes {
// if s is formatted as: X separator Y
// the bytes dedicated to X are from begin to pos,
// the ones dedicated to Y  from pos + 1 to s.size()
// (excluding the null terminator)
// WARNING: unsafe version that assumes separator
// is always present
byte_range split_upon(std::string const& s, const char separator) {
    const uint8_t* begin = reinterpret_cast<uint8_t const*>(s.c_str());
    auto pos = begin;
    for (; *pos != separator; ++pos)
        ;
    return {begin, pos};
}

// safe version of split_upon:
// used when separator could be missing
byte_range split_upon_check_end(std::string const& s, const char separator) {
    const uint8_t* begin = reinterpret_cast<uint8_t const*>(s.c_str());
    // exclude null terminator
    const uint8_t* end = begin + s.size();
    auto pos = begin;
    for (; pos != end; ++pos) {
        if (*pos == separator) {
            break;
        }
    }
    return {begin, pos};
}

byte_range predecessor(byte_range range) {
    // WARNING: unsafe version that assumes
    // whitespace is always present, use
    // commented code below for the safe variant
    auto pos = range.second - 1;
    while (*pos-- != ' ')
        ;
    return {range.first, pos + 1};
    // for (auto pos = range.second - 1;
    //           pos != range.first; --pos) {
    //     if (*pos == ' ') {
    //         return {range.first, pos};
    //     }
    // }
    // throw std::runtime_error("whitespace not found");
}

// return all but the first token
byte_range suffix(byte_range range) {
    auto pos = range.first;
    while (*pos++ != ' ')
        ;
    return {pos, range.second};
}

// returns the word back to the specified amount
// examples:
//     back_to(X Y Z W, 1) = Z
//     back_to(X Y Z W, 2) = Y
//     back_to(X Y Z W, 3) = X
byte_range back_to(byte_range range, uint32_t back_offset) {
    // WARNING: unsafe version that assumes
    // whitespace is always present
    auto pos = range.second - 1;
    uint32_t spaces = 0;
    auto end = pos;
    while (pos != range.first) {
        if (*pos == ' ') {
            ++spaces;
            if (spaces == back_offset) {
                end = pos;
            }
            if (spaces == back_offset + 1) {
                break;
            }
        }
        --pos;
    }
    return {pos + 1, end};
}

// returns the word from the beginning to the specified amount
// examples:
//     back_to(X Y Z W, 1) = Y
//     back_to(X Y Z W, 2) = Z
//     back_to(X Y Z W, 3) = W
byte_range to(byte_range range, uint32_t to_offset) {
    // WARNING: unsafe version that assumes
    // whitespace is always present
    auto pos = range.first;
    uint32_t spaces = 0;
    auto begin = pos;
    while (pos != range.second) {
        if (*pos == ' ') {
            ++spaces;
            if (spaces == to_offset) {
                begin = pos + 1;
            }
            if (spaces == to_offset + 1) {
                break;
            }
        }
        ++pos;
    }

    return {begin, pos};
}

byte_range next(byte_range range) {
    // WARNING: unsafe version that assumes
    // whitespace and next always present
    auto pos = range.second + 1;
    while (*pos++ != ' ')
        ;
    return {range.second + 1, pos - 1};
}

byte_range prev(byte_range range) {
    // WARNING: unsafe version that assumes
    // whitespace and next always present
    auto pos = range.first;
    --pos;  // skip space
    --pos;  // first non-whitespace character
    while (*pos-- != ' ')
        ;
    return {pos + 2, range.first - 1};
}

bool equal_bytes(byte_range x, byte_range y) {
    auto len = x.second - x.first;
    if (len == (y.second - y.first)) {
        return memcmp(x.first, y.first, len) == 0;
    }
    return false;
}
}  // namespace bytes

namespace tables {
const uint8_t select_in_byte[2048] = {
    8, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3,
    0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0,
    1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1,
    0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0,
    2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2,
    0, 1, 0, 7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0,
    1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1,
    0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0,
    3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0,
    1, 0, 2, 0, 1, 0, 8, 8, 8, 1, 8, 2, 2, 1, 8, 3, 3, 1, 3, 2, 2, 1, 8, 4, 4,
    1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1,
    3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 6, 6, 1, 6,
    2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2,
    2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2,
    1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 7, 7, 1, 7, 2, 2, 1, 7, 3, 3, 1, 3, 2, 2, 1,
    7, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 7, 5, 5, 1, 5, 2, 2, 1, 5,
    3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 7, 6,
    6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3,
    1, 3, 2, 2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1,
    4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 8, 8, 8, 8, 8, 8, 2, 8, 8, 8, 3, 8,
    3, 3, 2, 8, 8, 8, 4, 8, 4, 4, 2, 8, 4, 4, 3, 4, 3, 3, 2, 8, 8, 8, 5, 8, 5,
    5, 2, 8, 5, 5, 3, 5, 3, 3, 2, 8, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3,
    2, 8, 8, 8, 6, 8, 6, 6, 2, 8, 6, 6, 3, 6, 3, 3, 2, 8, 6, 6, 4, 6, 4, 4, 2,
    6, 4, 4, 3, 4, 3, 3, 2, 8, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2, 6,
    5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2, 8, 8, 8, 7, 8, 7, 7, 2, 8, 7,
    7, 3, 7, 3, 3, 2, 8, 7, 7, 4, 7, 4, 4, 2, 7, 4, 4, 3, 4, 3, 3, 2, 8, 7, 7,
    5, 7, 5, 5, 2, 7, 5, 5, 3, 5, 3, 3, 2, 7, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3,
    4, 3, 3, 2, 8, 7, 7, 6, 7, 6, 6, 2, 7, 6, 6, 3, 6, 3, 3, 2, 7, 6, 6, 4, 6,
    4, 4, 2, 6, 4, 4, 3, 4, 3, 3, 2, 7, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3,
    3, 2, 6, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 3, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 4, 8, 4, 4, 3,
    8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 5, 8, 5, 5, 3, 8, 8, 8, 5, 8, 5, 5, 4, 8,
    5, 5, 4, 5, 4, 4, 3, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 3, 8, 8,
    8, 6, 8, 6, 6, 4, 8, 6, 6, 4, 6, 4, 4, 3, 8, 8, 8, 6, 8, 6, 6, 5, 8, 6, 6,
    5, 6, 5, 5, 3, 8, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 8, 8, 8, 8,
    8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 3, 8, 8, 8, 7, 8, 7, 7, 4, 8, 7, 7, 4, 7,
    4, 4, 3, 8, 8, 8, 7, 8, 7, 7, 5, 8, 7, 7, 5, 7, 5, 5, 3, 8, 7, 7, 5, 7, 5,
    5, 4, 7, 5, 5, 4, 5, 4, 4, 3, 8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6,
    3, 8, 7, 7, 6, 7, 6, 6, 4, 7, 6, 6, 4, 6, 4, 4, 3, 8, 7, 7, 6, 7, 6, 6, 5,
    7, 6, 6, 5, 6, 5, 5, 3, 7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8,
    8, 8, 8, 8, 5, 8, 8, 8, 5, 8, 5, 5, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 4, 8, 8, 8, 8, 8,
    8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 5, 8, 8, 8, 6, 8, 6, 6, 5, 8, 6, 6, 5, 6, 5,
    5, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8,
    7, 8, 8, 8, 7, 8, 7, 7, 4, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 5,
    8, 8, 8, 7, 8, 7, 7, 5, 8, 7, 7, 5, 7, 5, 5, 4, 8, 8, 8, 8, 8, 8, 8, 7, 8,
    8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 4, 8, 8,
    8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 5, 8, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6,
    5, 6, 5, 5, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    6, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 6,
    8, 8, 8, 6, 8, 6, 6, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 5, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7,
    8, 7, 7, 6, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 7, 8,
    7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7};
}

namespace building_util {
void check_order(uint8_t order) {
    if (order < 2) {
        throw std::invalid_argument("order must be at least 2");
    }

    if (order > global::max_order) {
        throw std::invalid_argument("tongrams max supported order is " +
                                    std::to_string(global::max_order));
    }
}

void check_remapping_order(uint8_t order) {
    if (order > global::max_remapping_order) {
        throw std::invalid_argument(
            "tongrams max supported remapping order is " +
            std::to_string(global::max_remapping_order));
    }
}

void check_unk_logprob(float unk_logprob) {
    if (unk_logprob >= 0.0) {
        throw std::invalid_argument("log probability of <unk> must be < 0.0");
    }
}

void check_num_threads(uint64_t num_threads) {
    if (!num_threads) {
        throw std::invalid_argument("number of threads must be > 0");
    }
}

bool request_help(int argc, char** argv) {
    for (int i = 0; i < argc; ++i) {
        if (argv[i] == std::string("--help")) {
            return true;
        }
    }
    return false;
}

void display_legend() {
    std::cout << "Legend:\n"
              << "  " << style::bold << "bold     " << style::off
              << "  means mandatory argument\n"
              << "  " << style::underline << "underline" << style::off
              << "  means replaceable argument\n"
              << "  []         means optional argument" << std::endl
              << std::endl;
}

void print_general_params() {
    std::cout << "\t[--u " << style::underline << "value" << style::off << "]\n"
              << "\t[--p " << style::underline << "value" << style::off << "]\n"
              << "\t[--b " << style::underline << "value" << style::off << "]\n"
              << "\t[--arpa " << style::underline << "arpa_filename"
              << style::off << "]\n"
              << "\t[--out " << style::underline << "output_filename"
              << style::off << "]" << std::endl;
}

void print_general_info() {
    std::cout << "-------------------------------------------------------------"
                 "------------------------------------------\n"
              << style::bold << style::underline << "order" << style::off
              << " defines the maximum-order of the grams.\n"
              << style::bold << style::underline << "value_type" << style::off
              << " is either 'count' or 'prob_backoff'.\n"
              << "'--u' specifies the log10 probability of the unknown word.\n"
              << "'--p' specifies the number of quantization bits for "
                 "probability values (default is 8).\n"
              << "'--b' specifies the number of quantization bits for backoff "
                 "values (default is 8).\n"
              << "'--arpa' specifies the input arpa file.\n"
              << "'--out' specifies the output binary filename. If omitted it "
                 "will be given the name of the model type.\n"
              << "If 'count' is specified as " << style::bold
              << style::underline << "value_type" << style::off
              << ", then '--arpa' option will be ignored.\n"
              << "If '--p' and '--b' are only valid with 'prob_backoff' "
              << style::bold << style::underline << "value_type" << style::off
              << " specified." << std::endl;
}

void unknown_type(std::string const& type) {
    std::cerr << "Error: unknown type "
              << "'" << type << "'" << std::endl;
}

bool is_empty(std::ifstream& is) {
    return is.peek() == std::ifstream::traits_type::eof();
}

void write(std::ofstream& os, std::string const& line) {
    os.write(line.data(), line.size() * sizeof(char));
}
}  // namespace building_util

namespace util {

void not_found(std::string const& what) {
    throw std::runtime_error(("'" + what + "' not found.").c_str());
}

inline static uint32_t toul(byte_range const& br) {
    return std::strtoul(reinterpret_cast<const char*>(br.first), nullptr, 10);
}

inline static uint64_t toull(byte_range const& br) {
    return std::strtoull(reinterpret_cast<const char*>(br.first), nullptr, 10);
}

inline static uint64_t toull(const char* s) {
    return std::strtoull(s, nullptr, 10);
}

inline double get_time_usecs() {
    timeval tv;
    gettimeofday(&tv, NULL);
    return double(tv.tv_sec) * 1000000 + double(tv.tv_usec);
}

template <typename T>
inline void prefetch(T const* ptr) {
    _mm_prefetch((const char*)ptr, _MM_HINT_T0);
}

void check(uint64_t index, uint64_t got, uint64_t expected,
           std::string const& what) {
    if (got != expected) {
        std::cout << "Error at " << index << ":\n\t"
                  << "got " + what + " " << got << ", but "
                  << "expected " + what + " " << expected << std::endl;
        std::abort();
    }
}

void input_filename(const char* input_dir, uint8_t order,
                    std::string& filename) {
    filename = std::string(input_dir) + "/" + std::to_string(order) +
               "-grams.sorted.gz";
}

void check_filename(std::string const& filename) {
    std::ifstream is(filename.c_str());
    if (!is.good()) {
        std::cerr << "Expected file: '" << filename
                  << "', but it does not exist." << std::endl;
        std::abort();
    }
}

template <typename RandomAccessSequence, typename Adaptor>
bool binary_search(RandomAccessSequence const& s, uint64_t n, uint64_t x,
                   uint64_t& rank, Adaptor adaptor) {
    uint64_t lo = 0, hi = n - 1;
    while (lo <= hi) {
        uint64_t mid = lo + ((hi - lo) >> 1);
        uint64_t mid_v = adaptor.first(s, mid);
        if (mid_v == x) {
            rank = adaptor.second(s, mid);
            return true;
        }
        if (mid_v > x) {
            hi = mid - 1;
        } else {
            lo = mid + 1;
        }
    }
    return false;
}

template <typename T>
void save(uint8_t header, T const& data_structure,
          char const* output_filename) {
    if (output_filename == nullptr) {
        throw std::runtime_error(
            "You must specify the name of the output file.");
    }
    std::ofstream os(output_filename, std::ios::binary);
    essentials::save_pod(os, header);
    data_structure.save(os);
    os.close();
}

template <typename T>
size_t load(T& data_structure, std::string const& binary_filename) {
    std::ifstream is(binary_filename, std::ios::binary);
    if (!is.good()) {
        throw std::runtime_error(
            "Error in opening binary file, it may not exist or be malformed.");
    }
    uint8_t header = 0;
    essentials::load_pod(is, header);
    (void)header;  // skip header
    data_structure.load(is);
    size_t bytes = (size_t)is.tellg();
    is.close();
    return bytes;
}

std::string get_model_type(std::string const& binary_filename) {
    std::ifstream is(binary_filename, std::ios::binary);
    if (!is.good()) {
        throw std::runtime_error(
            "Error in opening binary file, it may not exist or be malformed.");
    }
    uint8_t header = 0;
    essentials::load_pod(is, header);
    binary_header bin_header;
    static constexpr bool verbose = true;
    auto model_string_type = bin_header.parse(header, verbose);
    is.close();
    return model_string_type;
}

inline uint8_t msb(uint64_t x) {
    assert(x);
    unsigned long ret = -1U;
    if (x) {
        ret = (unsigned long)(63 - __builtin_clzll(x));
    }
    return (uint8_t)ret;
}

inline bool bsr64(unsigned long* const index, const uint64_t mask) {
    if (mask) {
        *index = (unsigned long)(63 - __builtin_clzll(mask));
        return true;
    } else {
        return false;
    }
}

inline uint8_t msb(uint64_t x, unsigned long& ret) {
    return bsr64(&ret, x);
}

inline uint8_t lsb(uint64_t x, unsigned long& ret) {
    if (x) {
        ret = (unsigned long)__builtin_ctzll(x);
        return true;
    }
    return false;
}

inline uint8_t lsb(uint64_t x) {
    assert(x);
    unsigned long ret = -1U;
    lsb(x, ret);
    return (uint8_t)ret;
}

inline uint64_t ceil_log2(const uint64_t x) {
    return (x > 1) ? msb(x - 1) + 1 : 0;
}

inline uint64_t floor_log2(const uint64_t x) {
    return (x > 1) ? msb(x) : 0;
}

static const uint64_t ones_step_4 = 0x1111111111111111ULL;
static const uint64_t ones_step_8 = 0x0101010101010101ULL;
static const uint64_t msbs_step_8 = 0x80ULL * ones_step_8;

inline uint64_t byte_counts(uint64_t x) {
    x = x - ((x & 0xa * ones_step_4) >> 1);
    x = (x & 3 * ones_step_4) + ((x >> 2) & 3 * ones_step_4);
    x = (x + (x >> 4)) & 0x0f * ones_step_8;
    return x;
}

inline uint64_t bytes_sum(uint64_t x) {
    return x * ones_step_8 >> 56;
}

inline uint64_t popcount(uint64_t x) {
#if TONGRAMS_USE_POPCNT
    return uint64_t(_mm_popcnt_u64(x));
#else
    return bytes_sum(byte_counts(x));
#endif
}

// this is the select-in-word algorithm presented in
// "A Fast x86 Implementation of Select" by
// P. Pandey, M. A. Bender, and R. Johnson
// the algorithm uses only four x86 machine instructions,
// two of which were introduced in Intel’s Haswell CPUs in 2013
// source: https://github.com/splatlab/rankselect/blob/master/popcount.h
inline uint64_t select64_pdep_tzcnt(uint64_t x, const uint64_t k) {
    uint64_t i = 1ULL << k;
    asm("pdep %[x], %[mask], %[x]" : [ x ] "+r"(x) : [ mask ] "r"(i));
    asm("tzcnt %[bit], %[index]" : [ index ] "=r"(i) : [ bit ] "g"(x) : "cc");
    return i;
}

inline uint64_t select_in_word(const uint64_t x, const uint64_t k) {
    assert(k < popcount(x));
#if TONGRAMS_USE_PDEP
    return select64_pdep_tzcnt(x, k);
#else
    uint64_t byte_sums = byte_counts(x) * ones_step_8;
    const uint64_t k_step_8 = k * ones_step_8;
    const uint64_t geq_k_step_8 =
        (((k_step_8 | msbs_step_8) - byte_sums) & msbs_step_8);
    const uint64_t place = popcount(geq_k_step_8) * 8;
    const uint64_t byte_rank =
        k - (((byte_sums << 8) >> place) & uint64_t(0xFF));
    return place +
           tables::select_in_byte[((x >> place) & 0xFF) | (byte_rank << 8)];
#endif
}

template <typename IntType1, typename IntType2>
inline IntType1 ceil_div(IntType1 dividend, IntType2 divisor) {
    IntType1 d = IntType1(divisor);
    return IntType1(dividend + d - 1) / d;
}
}  // namespace util
}  // namespace tongrams
