#pragma once

#include <unistd.h>

#include "utils/mph_tables.hpp"
#include "utils/iterators.hpp"
#include "state.hpp"

namespace tongrams {

template <typename Values, typename KeyRankSequence, typename BaseHasher>
struct mph_prob_lm {
    typedef single_valued_mpht<KeyRankSequence, BaseHasher> hash_table;

    struct builder {
        builder() {}

        builder(const char* arpa_filename, uint8_t order, float unk_prob,
                uint8_t probs_quantization_bits,
                uint8_t backoffs_quantization_bits)
            : m_arpa_filename(arpa_filename)
            , m_order(order)
            , m_unk_prob(unk_prob) {
            building_util::check_order(m_order);
            m_tables.reserve(m_order);

            typename Values::builder probs_builder(m_order - 1);
            typename Values::builder backoffs_builder(m_order - 2);

            arpa_parser ap(arpa_filename);
            std::vector<uint64_t> counts;
            ap.read_header(counts);

            if (m_order > counts.size()) {
                throw std::invalid_argument(
                    "specified order exceeds arpa file order");
            }

            for (uint8_t order = 1; order <= m_order; ++order) {
                std::vector<float> probs;
                std::vector<float> backoffs;
                uint64_t n = counts[order - 1];
                probs.reserve(n);
                backoffs.reserve(n);

                ap.read_values(order, n, probs, backoffs);
                assert(probs.size() == n);

                if (order !=
                    1) {  // need scan uni-grams anyway to set arpa offsets
                    probs_builder.build_probs_sequence(probs,
                                                       probs_quantization_bits);
                    if (order != m_order) {
                        backoffs_builder.build_backoffs_sequence(
                            backoffs, backoffs_quantization_bits);
                    }
                }
            }

            ap.parse_eof();
            std::vector<uint64_t> arpa_offsets = ap.offsets();

            util::logger("Building vocabulary");
            build_vocabulary(arpa_offsets.front());

            size_t available_ram =
                sysconf(_SC_PAGESIZE) * sysconf(_SC_PHYS_PAGES);
            for (uint8_t order = 2; order <= m_order; ++order) {
                util::logger("Building " + std::to_string(order) + "-grams");
                arpa_iterator it(m_arpa_filename, order,
                                 arpa_offsets[order - 1]);
                uint64_t n = it.num_grams();
                grams_probs_pool pool(n, available_ram);

                for (uint64_t i = 0; i < n; ++i) {
                    auto const& tuple = it.next();
                    pool.append(tuple);
                }

                std::vector<byte_range> bytes;
                bytes.reserve(n);
                auto& pool_index = pool.index();

                compact_vector::builder cvb(
                    n, probs_quantization_bits +
                           (order != m_order ? backoffs_quantization_bits : 0));

                for (auto const& record : pool_index) {
                    bytes.push_back(record.gram);
                    float prob = record.prob;
                    float backoff = record.backoff;
                    // store interleaved ranks
                    uint64_t prob_rank = probs_builder.rank(order - 2, prob, 0);
                    uint64_t packed = prob_rank;
                    if (order != m_order) {
                        uint64_t backoff_rank =
                            backoffs_builder.rank(order - 2, backoff, 1);
                        packed |= backoff_rank << probs_quantization_bits;
                    }
                    cvb.push_back(packed);
                }

                m_tables.emplace_back(bytes, compact_vector(cvb),
                                      identity_adaptor());
            }

            probs_builder.build(m_probs_averages);
            backoffs_builder.build(m_backoffs_averages);
        }

        void build(mph_prob_lm<Values, KeyRankSequence, BaseHasher>& mph) {
            mph.m_order = m_order;
            mph.m_unk_prob = m_unk_prob;
            mph.m_probs_averages.swap(m_probs_averages);
            mph.m_backoffs_averages.swap(m_backoffs_averages);
            mph.m_tables.resize(m_order);
            for (uint8_t i = 0; i < m_order; ++i) {
                m_tables[i].build(mph.m_tables[i]);
            }
            builder().swap(*this);
        }

        void swap(builder& other) {
            std::swap(m_order, other.m_order);
            std::swap(m_unk_prob, other.m_unk_prob);
            m_probs_averages.swap(other.m_probs_averages);
            m_backoffs_averages.swap(other.m_backoffs_averages);
            m_tables.swap(other.m_tables);
        }

    private:
        char const* m_arpa_filename;
        uint8_t m_order;
        float m_unk_prob;
        Values m_probs_averages;
        Values m_backoffs_averages;
        std::vector<typename hash_table::builder> m_tables;

        void build_vocabulary(uint64_t unigrams_arpa_offset) {
            arpa_iterator it(m_arpa_filename, 1, unigrams_arpa_offset);
            uint64_t n = it.num_grams();

            size_t available_ram =
                sysconf(_SC_PAGESIZE) * sysconf(_SC_PHYS_PAGES);
            grams_probs_pool pool(n, available_ram);

            for (uint64_t i = 0; i < n; ++i) {
                auto const& record = it.next();
                auto gram = record.gram;
                // check if unk is present and NOT already set
                if (std::string(gram.first, gram.second) == "<unk>" &&
                    m_unk_prob == global::default_unk_prob) {
                    m_unk_prob = record.prob;
                    std::cout << "<unk> probability found to be: " << m_unk_prob
                              << std::endl;
                }
                pool.append(record);
            }

            std::vector<byte_range> bytes;
            bytes.reserve(n);
            auto& pool_index = pool.index();
            compact_vector::builder
                // unigrams' values are not quantized
                values_cvb(n, 64);

            for (auto const& record : pool_index) {
                bytes.push_back(record.gram);
                float prob = record.prob;
                float backoff = record.backoff;
                uint64_t packed = 0;
                bits::pack(packed, prob, backoff);
                values_cvb.push_back(packed);
            }

            m_tables.emplace_back(bytes, compact_vector(values_cvb),
                                  identity_adaptor());
        }
    };

    mph_prob_lm() : m_order(0), m_unk_prob(global::default_unk_prob) {}

    typedef prob_model_state<uint8_t const*> state_type;

    state_type state() {
        return state_type(order());
    }

    void score(state_type& state, byte_range const& word, bool& is_OOV,
               float& prob) {
        uint64_t value = m_tables[0].lookup(word, identity_adaptor());
        state.add_word(word.first);  // save beginning of word

        uint8_t longest_matching_history_len = 0;
        uint64_t order_m1 = 1;

        // STEP (1): determine longest matching history
        if (value != global::not_found) {
            float backoff;
            bits::unpack(value, prob, backoff);
            state.add_backoff(backoff);

            if (backoff) {
                longest_matching_history_len = 1;
            }

            auto words_rbegin = state.words.rbegin();
            ++words_rbegin;  // skip just added word

            for (; order_m1 <= state.length; ++order_m1, ++words_rbegin) {
                state.advance();

                auto prev_word_begin = *words_rbegin;
                byte_range gram(prev_word_begin, word.second);

                uint64_t rank =
                    m_tables[order_m1].lookup(gram, identity_adaptor());
                if (rank == global::not_found) {
                    break;
                }

                if (order_m1 != order() - 1) {
                    uint64_t probs_quantization_bits =
                        m_probs_averages.quantization_bits(order_m1 - 1);
                    uint64_t mask =
                        (uint64_t(1) << probs_quantization_bits) - 1;
                    uint64_t prob_rank = rank & mask;
                    uint64_t backoff_rank = rank >> probs_quantization_bits;
                    prob = m_probs_averages.access(order_m1 - 1, prob_rank);
                    backoff =
                        m_backoffs_averages.access(order_m1 - 1, backoff_rank);
                    state.add_backoff(backoff);
                    if (backoff) {
                        longest_matching_history_len = order_m1 + 1;
                    }

                } else {
                    prob = m_probs_averages.access(order_m1 - 1, rank);
                }
            }

        } else {  // unseen word
            ++state.OOVs;
            is_OOV = true;
            prob = m_unk_prob;
            state.add_backoff(0.0);
        }

        // if we encountered unseen ngrams during STEP (1)
        for (uint64_t i = order_m1 - 1; i < state.length; ++i) {
            prob += state.backoff(i);
        }

        state.length = longest_matching_history_len;
        state.finalize();
        assert(prob < 0.0);
    }

    void print_stats(size_t bytes) const;

    uint64_t order() const {
        return uint64_t(m_order);
    }

    size_t size() const {
        size_t size = 0;
        for (auto const& t : m_tables) {
            size += t.size();
        }
        return size;
    }

    void save(std::ostream& os) const {
        util::save_pod(os, &m_order);
        util::save_pod(os, &m_unk_prob);
        m_probs_averages.save(os);
        m_backoffs_averages.save(os);
        for (auto const& t : m_tables) {
            t.save(os);
        }
    }

    void load(std::istream& is) {
        util::load_pod(is, &m_order);
        util::load_pod(is, &m_unk_prob);
        m_probs_averages.load(is, m_order - 1);
        m_backoffs_averages.load(is, m_order - 2);
        m_tables.resize(m_order);
        for (auto& t : m_tables) {
            t.load(is);
        }
    }

private:
    uint8_t m_order;
    float m_unk_prob;
    Values m_probs_averages;
    Values m_backoffs_averages;
    std::vector<hash_table> m_tables;
};
}  // namespace tongrams
