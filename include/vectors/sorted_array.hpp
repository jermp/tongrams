#pragma once

#include "../sequences/pointer_sequence.hpp"
#include "../utils/util.hpp"

namespace tongrams {
template <typename Grams, typename Ranks, typename Pointers>
struct sorted_array {
    struct estimation_builder {
        estimation_builder() {}

        template <typename T>
        void build_word_ids(uint8_t order, sorted_array& sa, T& partition) {
            sa.m_grams.build(word_ids.begin(), word_ids.size(), partition,
                             order);
            compact_vector::builder().swap(word_ids);
        }

        void build_probs_backoffs_ranks(sorted_array& sa) {
            probs_backoffs_ranks.build(sa.m_probs_backoffs_ranks);
            compact_vector::builder().swap(probs_backoffs_ranks);
        }

        void build_pointers(sorted_array& sa) {
            sa.m_pointers.build(pointers);
            compact_vector::builder().swap(pointers);
        }

        compact_vector::builder word_ids;
        compact_vector::builder probs_backoffs_ranks;
        compact_vector::builder pointers;
    };

    struct builder {
        builder() {}

        builder(uint64_t num_grams, uint64_t max_gram_id,
                uint64_t max_count_rank, uint8_t quantization_bits)
            : m_size(num_grams)
            , m_grams(num_grams, util::ceil_log2(max_gram_id + 1))
            , m_counts_ranks(num_grams, util::ceil_log2(max_count_rank + 1))
            , m_probs_backoffs_ranks(num_grams, quantization_bits) {}

        void add_gram(uint64_t id) {
            m_grams.push_back(id);
        }

        void add_count_rank(uint64_t rank) {
            m_counts_ranks.push_back(rank);
        }

        // prob and backoff ranks are stored interleaved
        void add_prob_backoff_rank(uint64_t rank) {
            m_probs_backoffs_ranks.push_back(rank);
        }

        template <typename T>
        void build(sorted_array& sa, T& pointers, uint8_t order, int value_t) {
            compact_vector grams(m_grams);
            sa.m_grams.build(grams.begin(), grams.size(), pointers, order);

            switch (value_t) {
                case value_type::count:
                    build_counts_ranks(sa, order);
                    break;
                case value_type::prob_backoff:
                    build_probs_backoffs_ranks(sa, order);
                    break;
                case value_type::none:
                    break;
                default:
                    assert(false);
            }

            builder().swap(*this);
        }

        void build_counts_ranks(sorted_array& sa, uint8_t order) {
            compact_vector counts_ranks(m_counts_ranks);
            sa.m_counts_ranks.build(counts_ranks.begin(), counts_ranks.size(),
                                    order);
            compact_vector::builder().swap(m_counts_ranks);
        }

        void build_probs_backoffs_ranks(sorted_array& sa, uint8_t order) {
            compact_vector probs_backoffs_ranks(m_probs_backoffs_ranks);
            sa.m_probs_backoffs_ranks.build(probs_backoffs_ranks.begin(),
                                            probs_backoffs_ranks.size(), order);
            compact_vector::builder().swap(m_probs_backoffs_ranks);
        }

        template <typename T>
        static void build_pointers(sorted_array& sa, T& pointers) {
            sa.m_pointers.build(pointers);
        }

        void swap(builder& other) {
            std::swap(m_size, other.m_size);
            m_grams.swap(other.m_grams);
            m_counts_ranks.swap(other.m_counts_ranks);
            m_probs_backoffs_ranks.swap(other.m_probs_backoffs_ranks);
        }

    private:
        uint64_t m_size;
        compact_vector::builder m_grams;
        compact_vector::builder m_counts_ranks;
        compact_vector::builder m_probs_backoffs_ranks;
    };

    sorted_array() {}

    sorted_array(uint64_t size) : m_size(size) {}

    inline pointer_range range(uint64_t pos) {
        assert(pos < size());
        return m_pointers[pos];
    }

    inline uint64_t next(pointer_range& r, uint64_t id) {
        uint64_t pos = position(r, id);
        if (pos == global::not_found) {
            return global::not_found;
        }
        r = range(pos);
        return pos;
    }

    inline uint64_t count_rank(uint64_t pos) {
        assert(pos < size());
        return m_counts_ranks[pos];
    }

    inline uint64_t prob_backoff_rank(uint64_t pos) const {
        assert(pos < size());
        return m_probs_backoffs_ranks.access(pos);
    }

    inline uint64_t position(pointer_range r, uint64_t id) {
        uint64_t pos = 0;
        m_grams.find(r, id, &pos);
        return pos;
    }

    Grams* grams() {
        return &m_grams;
    }

    Ranks const* counts_ranks() const {
        return &m_counts_ranks;
    }

    Ranks const& probs_backoffs_ranks() const {
        return m_probs_backoffs_ranks;
    }

    Ranks const& ranks() const {
        return m_counts_ranks;
    }

    pointer_sequence<Pointers> const* ptrs() const {
        return &m_pointers;
    }

    uint64_t size() const {
        return m_size;
    }

    void print_stats(uint8_t order,
                     pointer_sequence<Pointers> const* ranges) const;

    uint64_t grams_bytes() const {
        return m_grams.bytes();
    }

    uint64_t counts_ranks_bytes() const {
        return m_counts_ranks.bytes();
    }

    uint64_t probs_backoffs_ranks_bytes() const {
        return m_probs_backoffs_ranks.bytes();
    }

    uint64_t pointers_bytes() const {
        return m_pointers.bytes();
    }

    void save(std::ostream& os, uint8_t order, int value_t) const {
        essentials::save_pod(os, m_size);

        if (order != 1) {
            m_grams.save(os);
        }

        switch (value_t) {
            case value_type::count:
                m_counts_ranks.save(os);
                break;
            case value_type::prob_backoff:
                m_probs_backoffs_ranks.save(os);
                break;
            case value_type::none:
                break;
            default:
                assert(false);
        }

        m_pointers.save(os);
    }

    void load(std::istream& is, uint8_t order, int value_t) {
        essentials::load_pod(is, m_size);

        if (order != 1) {
            m_grams.load(is);
        }

        switch (value_t) {
            case value_type::count:
                m_counts_ranks.load(is);
                break;
            case value_type::prob_backoff:
                m_probs_backoffs_ranks.load(is);
                break;
            case value_type::none:
                break;
            default:
                assert(false);
        }

        m_pointers.load(is);
    }

private:
    uint64_t m_size;
    Grams m_grams;
    Ranks m_counts_ranks;
    Ranks m_probs_backoffs_ranks;
    pointer_sequence<Pointers> m_pointers;
};
}  // namespace tongrams
