#pragma once

#include <vector>

namespace tongrams {

struct pointer_range {
    uint64_t begin;
    uint64_t end;
};

enum data_structure_type {
    hash = 0,      // minimal perfect hash (MPH)
    ef_trie = 1,   // Elias-Fano trie
    pef_trie = 2,  // partitioned Elias-Fano trie
};

enum value_type { count = 0, prob_backoff = 1, none = 2 };

enum ranks_type {
    IC = 0,    // indexed codewords
    PSEF = 1,  // prefix-sums + Elias-Fano
    PSPEF = 2  // prefix-sums + partitioned Elias-Fano
};

typedef std::pair<uint64_t, uint64_t> uint64_pair;

typedef std::tuple<uint64_t, uint64_t, uint64_t> uint64_triplet;

typedef std::vector<uint64_pair> pairs_vector;

typedef std::pair<uint8_t const*, uint8_t const*> byte_range;

struct count_record {
    count_record() : count(0) {}

    count_record(byte_range g, uint64_t c) : gram(g), count(c) {}

    byte_range gram;
    uint64_t count;
};

struct prob_backoff_record {
    prob_backoff_record() : prob(0.0), backoff(0.0) {}

    prob_backoff_record(byte_range g, float p, float b)
        : gram(g), prob(p), backoff(b) {}

    byte_range gram;
    float prob, backoff;
};

struct identity_adaptor {
    byte_range operator()(byte_range br) const {
        return br;
    }
};

// NOTE: throughout the code we do NOT include the null terminator '\0'
// when building a minimal perfect hash function over a set of strings,
// therefore the following adaptor should be used for clients willing
// to use std::string
struct stl_string_adaptor {
    byte_range operator()(std::string const& s) const {
        const uint8_t* buf = reinterpret_cast<uint8_t const*>(s.c_str());
        const uint8_t* end = buf + s.size();  // exclude the null terminator
        return {buf, end};
    }
};

// NOTE: we could use sizeof to work with any kind of unsigned ints
struct uint64_adaptor {
    byte_range operator()(uint64_t const& x) const {
        const uint8_t* buf = reinterpret_cast<const uint8_t*>(&x);
        return {buf, buf + 8};
    }
};

}  // namespace tongrams
