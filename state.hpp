#pragma once

#include "utils/util.hpp"
#include "utils/circular_buffer.hpp"

namespace tongrams
{
    // NOTE:
    // the template parameter Word is uint64_t for a trie model (ID of a word);
    // it is a uint8_t const* for a hash model (pointer to the last character of the word).
    template<typename Word>
    struct prob_model_state {
        prob_model_state(uint8_t max_context_length)
            : length(0)
            , words(max_context_length)
            , OOVs(0)
            , m_pos(0)
        {
            std::fill(m_curr_backoffs,
                      m_curr_backoffs + max_context_length - 1, 1.0);
            std::fill(m_prev_backoffs,
                      m_prev_backoffs + max_context_length - 1, 1.0);
        }

        inline void advance() {
            ++m_pos;
        }

        void finalize() {
            m_pos = 0;
            std::swap_ranges(std::begin(m_curr_backoffs),
                             std::end(m_curr_backoffs),
                             std::begin(m_prev_backoffs));
        }

        inline void add_word(Word word) {
            words.push_back(word);
        }

        inline void add_backoff(float backoff) {
            m_curr_backoffs[m_pos] = backoff;
        }

        inline float backoff(size_t i) const {
            return m_prev_backoffs[i];
        }

        inline void init() {
            words.init();
            std::fill(m_prev_backoffs,
                      m_prev_backoffs + length, 0.0);
            length = 0;
            OOVs = 0;
            m_pos = 0;
        }

        uint8_t length;
        circular_buffer<Word> words;
        uint64_t OOVs;

    private:
        uint64_t m_pos;
        float m_curr_backoffs[global::max_order];
        float m_prev_backoffs[global::max_order];
    };
}
