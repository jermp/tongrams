#pragma once

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/cat.hpp>

#include "../external/emphf/base_hash.hpp"
#include "utils/mph_tables.hpp"
#include "mph_count_lm.hpp"
#include "mph_prob_lm.hpp"
#include "trie_count_lm.hpp"
#include "trie_prob_lm.hpp"
#include "mappers.hpp"
#include "vectors/hash_compact_vector.hpp"
#include "sequences/ef_sequence.hpp"
#include "sequences/sequence_collections.hpp"
#include "sequences/indexed_codewords_sequence.hpp"
#include "sequences/pointer_sequence.hpp"
#include "sequences/fast_ef_sequence.hpp"
#include "sequences/prefix_summed_sequence.hpp"
#include "sequences/uniform_pef_sequence.hpp"

namespace tongrams {

#define TONGRAMS_MPH_COUNT_TYPE(HASH_KEY_BITS)                 \
    mph_count_lm<sequence_collection,                          \
                 hash_compact_vector<uint##HASH_KEY_BITS##_t>, \
                 emphf::jenkins##HASH_KEY_BITS##_hasher>

#define TONGRAMS_MPH_PROB_TYPE(HASH_KEY_BITS)                 \
    mph_prob_lm<quantized_sequence_collection,                \
                hash_compact_vector<uint##HASH_KEY_BITS##_t>, \
                emphf::jenkins##HASH_KEY_BITS##_hasher>

typedef TONGRAMS_MPH_COUNT_TYPE(32) mph32_count_lm;
typedef TONGRAMS_MPH_COUNT_TYPE(64) mph64_count_lm;
typedef TONGRAMS_MPH_PROB_TYPE(32) mph32_prob_lm;
typedef TONGRAMS_MPH_PROB_TYPE(64) mph64_prob_lm;

#define TONGRAMS_TRIE_COUNT_TYPE(MAPPER, COUNT_RANKS, GRAM_SEQUENCE_TYPE) \
    trie_count_lm<single_valued_mpht64, MAPPER, sequence_collection,      \
                  COUNT_RANKS, GRAM_SEQUENCE_TYPE, ef_sequence>

typedef TONGRAMS_TRIE_COUNT_TYPE(identity_mapper, indexed_codewords_sequence,
                                 fast_ef_sequence) ef_trie_IC_ranks_count_lm;
typedef TONGRAMS_TRIE_COUNT_TYPE(identity_mapper,
                                 prefix_summed_sequence<ef_sequence>,
                                 fast_ef_sequence) ef_trie_PSEF_ranks_count_lm;
typedef TONGRAMS_TRIE_COUNT_TYPE(
    identity_mapper, prefix_summed_sequence<pef::uniform_pef_sequence>,
    fast_ef_sequence) ef_trie_PSPEF_ranks_count_lm;
typedef TONGRAMS_TRIE_COUNT_TYPE(sorted_array_mapper,
                                 indexed_codewords_sequence,
                                 fast_ef_sequence) ef_rtrie_IC_ranks_count_lm;
typedef TONGRAMS_TRIE_COUNT_TYPE(sorted_array_mapper,
                                 prefix_summed_sequence<ef_sequence>,
                                 fast_ef_sequence) ef_rtrie_PSEF_ranks_count_lm;
typedef TONGRAMS_TRIE_COUNT_TYPE(
    sorted_array_mapper, prefix_summed_sequence<pef::uniform_pef_sequence>,
    fast_ef_sequence) ef_rtrie_PSPEF_ranks_count_lm;

typedef TONGRAMS_TRIE_COUNT_TYPE(identity_mapper, indexed_codewords_sequence,
                                 pef::uniform_pef_sequence)
    pef_trie_IC_ranks_count_lm;
typedef TONGRAMS_TRIE_COUNT_TYPE(
    identity_mapper, prefix_summed_sequence<ef_sequence>,
    pef::uniform_pef_sequence) pef_trie_PSEF_ranks_count_lm;
typedef TONGRAMS_TRIE_COUNT_TYPE(
    identity_mapper, prefix_summed_sequence<pef::uniform_pef_sequence>,
    pef::uniform_pef_sequence) pef_trie_PSPEF_ranks_count_lm;
typedef TONGRAMS_TRIE_COUNT_TYPE(
    sorted_array_mapper, indexed_codewords_sequence,
    pef::uniform_pef_sequence) pef_rtrie_IC_ranks_count_lm;
typedef TONGRAMS_TRIE_COUNT_TYPE(
    sorted_array_mapper, prefix_summed_sequence<ef_sequence>,
    pef::uniform_pef_sequence) pef_rtrie_PSEF_ranks_count_lm;
typedef TONGRAMS_TRIE_COUNT_TYPE(
    sorted_array_mapper, prefix_summed_sequence<pef::uniform_pef_sequence>,
    pef::uniform_pef_sequence) pef_rtrie_PSPEF_ranks_count_lm;

#define TONGRAMS_TRIE_PROB_TYPE(MAPPER, GRAM_SEQUENCE_TYPE)                   \
    trie_prob_lm<double_valued_mpht64, MAPPER, quantized_sequence_collection, \
                 compact_vector, GRAM_SEQUENCE_TYPE, ef_sequence>

typedef TONGRAMS_TRIE_PROB_TYPE(identity_mapper,
                                fast_ef_sequence) ef_trie_prob_lm;
typedef TONGRAMS_TRIE_PROB_TYPE(identity_mapper,
                                pef::uniform_pef_sequence) pef_trie_prob_lm;
typedef TONGRAMS_TRIE_PROB_TYPE(sorted_array_mapper,
                                fast_ef_sequence) ef_rtrie_prob_lm;
typedef TONGRAMS_TRIE_PROB_TYPE(sorted_array_mapper,
                                pef::uniform_pef_sequence) pef_rtrie_prob_lm;

// for print_stats.cpp
#define TONGRAMS_TYPES                                                      \
    (mph32_count_lm)(mph64_count_lm)(ef_trie_IC_ranks_count_lm)(            \
        ef_trie_PSEF_ranks_count_lm)(ef_trie_PSPEF_ranks_count_lm)(         \
        ef_rtrie_IC_ranks_count_lm)(ef_rtrie_PSEF_ranks_count_lm)(          \
        ef_rtrie_PSPEF_ranks_count_lm)(pef_trie_IC_ranks_count_lm)(         \
        pef_trie_PSEF_ranks_count_lm)(pef_trie_PSPEF_ranks_count_lm)(       \
        pef_rtrie_IC_ranks_count_lm)(pef_rtrie_PSEF_ranks_count_lm)(        \
        pef_rtrie_PSPEF_ranks_count_lm)(ef_trie_prob_lm)(pef_trie_prob_lm)( \
        ef_rtrie_prob_lm)(pef_rtrie_prob_lm)(mph32_prob_lm)(mph64_prob_lm)

// for check_count_model.cpp
//     lookup_perf_test.cpp
#define TONGRAMS_COUNT_TYPES                                          \
    (mph32_count_lm)(mph64_count_lm)(ef_trie_IC_ranks_count_lm)(      \
        ef_trie_PSEF_ranks_count_lm)(ef_trie_PSPEF_ranks_count_lm)(   \
        ef_rtrie_IC_ranks_count_lm)(ef_rtrie_PSEF_ranks_count_lm)(    \
        ef_rtrie_PSPEF_ranks_count_lm)(pef_trie_IC_ranks_count_lm)(   \
        pef_trie_PSEF_ranks_count_lm)(pef_trie_PSPEF_ranks_count_lm)( \
        pef_rtrie_IC_ranks_count_lm)(pef_rtrie_PSEF_ranks_count_lm)(  \
        pef_rtrie_PSPEF_ranks_count_lm)

// for build_mph_lm.cpp
#define TONGRAMS_HASH_COUNT_TYPES (mph32_count_lm)(mph64_count_lm)

// for build_mph_lm.cpp
#define TONGRAMS_HASH_PROB_TYPES (mph32_prob_lm)(mph64_prob_lm)

// for build_trie_lm.cpp
#define TONGRAMS_TRIE_COUNT_TYPES                                     \
    (ef_trie_IC_ranks_count_lm)(ef_trie_PSEF_ranks_count_lm)(         \
        ef_trie_PSPEF_ranks_count_lm)(ef_rtrie_IC_ranks_count_lm)(    \
        ef_rtrie_PSEF_ranks_count_lm)(ef_rtrie_PSPEF_ranks_count_lm)( \
        pef_trie_IC_ranks_count_lm)(pef_trie_PSEF_ranks_count_lm)(    \
        pef_trie_PSPEF_ranks_count_lm)(pef_rtrie_IC_ranks_count_lm)(  \
        pef_rtrie_PSEF_ranks_count_lm)(pef_rtrie_PSPEF_ranks_count_lm)

// for build_trie_lm.cpp
#define TONGRAMS_TRIE_PROB_TYPES \
    (ef_trie_prob_lm)(pef_trie_prob_lm)(ef_rtrie_prob_lm)(pef_rtrie_prob_lm)

// for score.cpp
#define TONGRAMS_SCORE_TYPES                                                  \
    (ef_trie_prob_lm)(pef_trie_prob_lm)(ef_rtrie_prob_lm)(pef_rtrie_prob_lm)( \
        mph32_prob_lm)(mph64_prob_lm)

}  // namespace tongrams