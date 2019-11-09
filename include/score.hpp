#pragma once

#include "trie_prob_lm.hpp"
#include "mph_prob_lm.hpp"

namespace tongrams {

template <typename Vocabulary, typename Mapper, typename Values, typename Ranks,
          typename Grams, typename Pointers>
float trie_prob_lm<Vocabulary, Mapper, Values, Ranks, Grams, Pointers>::score(
    state_type& state, byte_range const word, bool& is_OOV) {
    uint64_pair word_id = m_vocab.lookup_pair(word, identity_adaptor());
    state.add_word(word_id.first);

    uint8_t longest_matching_history_len = 0;
    uint64_t order_m1 = 1;
    float prob = 0.0;

    // STEP (1): determine longest matching history
    if (word_id.first != global::not_found) {
        float backoff;
        bits::unpack(word_id.second, prob, backoff);
        state.add_backoff(backoff);

        if (backoff) {
            longest_matching_history_len = 1;
        }

        auto words_rbegin = state.words.rbegin();
        ++words_rbegin;  // skip just added word id

        // needed for remapping
        uint64_t prev_id = word_id.first;
        uint64_t prev_prev_id = prev_id;

        auto r = m_arrays[0].range(word_id.first);

        for (; order_m1 <= state.length; ++order_m1, ++words_rbegin) {
            state.advance();

            if (r.end - r.begin == 0) {
                // no extension to the left, i.e.,
                // no successors in reversed trie
                break;
            }

            uint64_t id = *words_rbegin;

            if (Mapper::context_remapping && order_m1 > m_remapping_order) {
                id = m_mapper.map_id(
                    prev_id,
                    prev_prev_id,  // pass the two parent ids for remapping
                    id, &m_arrays.front(), m_remapping_order);
            }

            uint64_t pos = m_arrays[order_m1].position(r, id);
            if (pos == global::not_found) {
                break;
            }

            uint64_t probs_quantization_bits =
                m_probs_averages.quantization_bits(order_m1 - 1);
            uint64_t mask = (uint64_t(1) << probs_quantization_bits) - 1;
            uint64_t prob_backoff_rank =
                m_arrays[order_m1].prob_backoff_rank(pos);
            uint64_t prob_rank = prob_backoff_rank & mask;
            uint64_t backoff_rank =
                prob_backoff_rank >> probs_quantization_bits;
            prob = m_probs_averages.access(order_m1 - 1, prob_rank);

            if (order_m1 != order() - 1) {
                backoff =
                    m_backoffs_averages.access(order_m1 - 1, backoff_rank);
                state.add_backoff(backoff);
                r = m_arrays[order_m1].range(pos);
                if (backoff) {
                    longest_matching_history_len = order_m1 + 1;
                }
            }

            prev_prev_id = prev_id;
            prev_id = *words_rbegin;
        }

    } else {  // unseen word
        ++state.OOVs;
        is_OOV = true;
        prob = m_unk_prob;
        state.add_backoff(0.0);
    }

    // STEP (2): add backoff weights
    // if we encountered unseen ngrams during STEP (1)
    for (uint64_t i = order_m1 - 1; i < state.length; ++i) {
        prob += state.backoff(i);
    }

    state.length = longest_matching_history_len;
    state.finalize();
    assert(prob < 0.0);
    return prob;
}

}  // namespace tongrams