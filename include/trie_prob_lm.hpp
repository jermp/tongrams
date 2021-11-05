#pragma once

#include "utils/util.hpp"
#include "state.hpp"
#include "utils/iterators.hpp"
#include "vectors/sorted_array.hpp"
#include "../external/essentials/include/essentials.hpp"

namespace tongrams {

template <typename Vocabulary, typename Mapper, typename Values, typename Ranks,
          typename Grams, typename Pointers>
struct trie_prob_lm {
    typedef sorted_array<Grams, Ranks, Pointers> sorted_array_type;

    struct estimation_builder;

    struct builder {
        builder() {}

        builder(const char* arpa_filename, uint8_t order,
                uint8_t remapping_order, float unk_prob,
                uint8_t probs_quantization_bits,
                uint8_t backoffs_quantization_bits)
            : m_arpa_filename(arpa_filename)
            , m_order(order)
            , m_remapping_order(remapping_order)
            , m_unk_prob(unk_prob) {
            building_util::check_order(m_order);
            building_util::check_remapping_order(m_remapping_order);
            m_arrays.reserve(m_order);

            typename Values::builder probs_builder(m_order - 1);
            typename Values::builder backoffs_builder(m_order - 2);

            arpa_parser ap(m_arpa_filename);
            std::vector<uint64_t> counts;
            ap.read_header(counts);

            if (m_order > counts.size()) {
                throw std::invalid_argument(
                    "specified order exceeds arpa file order");
            }

            for (uint8_t ord = 1; ord <= m_order; ++ord) {
                std::vector<float> probs;
                std::vector<float> backoffs;
                uint64_t n = counts[ord - 1];
                probs.reserve(n);
                backoffs.reserve(n);

                m_arrays.push_back(sorted_array_type(n));
                essentials::logger("Reading " + std::to_string(ord) +
                                   "-grams probs/backoffs");
                ap.read_values(ord, n, probs, backoffs);
                assert(probs.size() == n);

                if (ord != 1) {  // unigrams are NOT quantized
                    probs_builder.build_probs_sequence(probs,
                                                       probs_quantization_bits);
                    if (ord != m_order) {
                        backoffs_builder.build_backoffs_sequence(
                            backoffs, backoffs_quantization_bits);
                    }
                }
            }

            ap.parse_eof();
            std::vector<uint64_t> arpa_offsets = ap.offsets();

            essentials::logger("Building vocabulary");
            build_vocabulary(arpa_offsets.front());

            for (uint8_t ord = 2; ord <= m_order; ++ord) {
                std::string order_grams(std::to_string(ord) + "-grams");
                uint64_t n = counts[ord - 1];

                typename sorted_array_type::builder sa_builder(
                    n,
                    m_vocab.size(),            // max_gram_id
                    0,                         // max_count_rank not used
                    probs_quantization_bits +  // quantization_bits
                        (ord != m_order ? backoffs_quantization_bits : 0));
                uint64_t num_pointers = counts[ord - 2] + 1;

                // NOTE: we could use this to save pointers' space at building
                // time compact_vector::builder pointers(num_pointers,
                // util::ceil_log2(n + 1));
                std::vector<uint64_t> pointers;
                pointers.reserve(num_pointers);

                essentials::logger("Building " + order_grams);
                build_ngrams(ord, arpa_offsets[ord - 1], arpa_offsets[ord - 2],
                             probs_builder, backoffs_builder, pointers,
                             sa_builder);
                assert(pointers.back() == n);
                assert(pointers.size() == num_pointers);
                essentials::logger("Writing " + order_grams);
                sa_builder.build(m_arrays[ord - 1], pointers, ord,
                                 value_type::prob_backoff);
                essentials::logger("Writing pointers");
                sorted_array_type::builder::build_pointers(m_arrays[ord - 2],
                                                           pointers);
            }

            probs_builder.build(m_probs_averages);
            backoffs_builder.build(m_backoffs_averages);
        }

        void build(trie_prob_lm& trie) {
            trie.m_order = m_order;
            trie.m_remapping_order = m_remapping_order;
            trie.m_unk_prob = m_unk_prob;
            trie.m_probs_averages.swap(m_probs_averages);
            trie.m_backoffs_averages.swap(m_backoffs_averages);
            trie.m_vocab.swap(m_vocab);
            trie.m_arrays.swap(m_arrays);
            builder().swap(*this);
        }

        void swap(builder& other) {
            std::swap(m_order, other.m_order);
            std::swap(m_remapping_order, other.m_remapping_order);
            std::swap(m_unk_prob, other.m_unk_prob);
            m_probs_averages.swap(other.m_probs_averages);
            m_backoffs_averages.swap(other.m_backoffs_averages);
            m_vocab.swap(other.m_vocab);
            m_arrays.swap(other.m_arrays);
        }

    private:
        char const* m_arpa_filename;
        uint8_t m_order;
        uint8_t m_remapping_order;
        float m_unk_prob;
        Mapper m_mapper;
        Values m_probs_averages;
        Values m_backoffs_averages;
        Vocabulary m_vocab;
        std::vector<sorted_array_type> m_arrays;

        void build_vocabulary(uint64_t unigrams_arpa_offset) {
            arpa_iterator it(m_arpa_filename, 1, unigrams_arpa_offset);
            uint64_t n = it.num_grams();

            strings_pool unigrams_pool;
            unigrams_pool.reserve(
                (sysconf(_SC_PAGESIZE) * sysconf(_SC_PHYS_PAGES)) * 0.8);
            std::vector<size_t> offsets;
            offsets.reserve(n);
            offsets.push_back(0);

            // unigrams' values are not quantized
            compact_vector::builder values_cvb(n, 64);

            for (uint64_t i = 0; i < n; ++i) {
                auto const& record = it.next();
                auto gram = record.gram;
                float prob = record.prob;
                float backoff = record.backoff;

                // check if <unk> is present and NOT already set
                if (std::string(gram.first, gram.second) == "<unk>" &&
                    m_unk_prob == global::default_unk_prob) {
                    m_unk_prob = prob;
                    std::cout << "<unk> probability found to be: " << m_unk_prob
                              << std::endl;
                }

                unigrams_pool.append(gram);
                offsets.push_back(unigrams_pool.bytes());
                uint64_t packed = 0;
                bits::pack(packed, prob, backoff);
                values_cvb.push_back(packed);
            }

            std::vector<byte_range> bytes;
            bytes.reserve(n);

            auto base_addr = unigrams_pool.base_addr();
            for (uint64_t i = 0; i < offsets.size() - 1; ++i) {
                bytes.emplace_back(unigrams_pool.get_bytes(
                    base_addr, offsets[i], offsets[i + 1]));
            }

            compact_vector::builder ids_cvb(n, util::ceil_log2(n + 1));
            for (uint64_t id = 0; id < n; ++id) {
                ids_cvb.push_back(id);
            }

            // NOTE:
            // build vocabulary excluding null terminators
            // from unigrams strings so that we can lookup
            // for any substring of a n-gram
            // without allocating a std::string
            m_vocab.build(
                bytes,
                compact_vector(),  // empty vector to use default hash-keys
                compact_vector(ids_cvb), compact_vector(values_cvb),
                identity_adaptor());
        }

        template <typename T>
        void build_ngrams(uint8_t order, uint64_t curr_order_arpa_offset,
                          uint64_t prev_order_arpa_offset,
                          typename Values::builder const& probs_builder,
                          typename Values::builder const& backoffs_builder,
                          T& pointers,
                          typename sorted_array_type::builder& sa_builder) {
            pointers.push_back(0);
            identity_adaptor adaptor;

            uint64_t pos = 0;

            arpa_iterator curr_it(m_arpa_filename, order,
                                  curr_order_arpa_offset);
            arpa_iterator prev_it(m_arpa_filename, order - 1,
                                  prev_order_arpa_offset);
            uint64_t curr_n = curr_it.num_grams();
            uint64_t prev_n = prev_it.num_grams();
            auto prev_record = prev_it.next();
            uint8_t probs_quantization_bits =
                probs_builder.quantization_bits(order - 2);

            uint64_t j = 0;
            for (uint64_t i = 0; i < curr_n; ++i) {
                auto const& record = curr_it.next();
                auto gram = record.gram;

                // NOTE:
                // in a BACKWARD trie, 'pattern' is the suffix of 'gram'
                // and 'token' is the first token of 'gram'
                byte_range pattern = bytes::suffix(gram);
                byte_range token(gram.first, pattern.first - 1);

                while (j != prev_n  // NOTE:
                                    // this test is here only to
                                    // guarantee termination in
                                    // case of wrong data:
                                    // 'pattern' should ALWAYS
                                    // be found within previous order grams
                       && !bytes::equal_bytes(pattern, prev_record.gram)) {
                    pointers.push_back(pos);
                    prev_record = prev_it.next();
                    ++j;
                }

                // check correctness of arpa file
                if (j == prev_n) {
                    std::cerr << "arpa file contains wrong data:\n";
                    std::cerr
                        << "\t'" << std::string(pattern.first, pattern.second)
                        << "'"
                        << " should have been found within previous order grams"
                        << std::endl;
                    exit(1);
                }

                ++pos;

                uint64_pair token_id = m_vocab.lookup_pair(token, adaptor);

                // apply remapping if Mapper has to
                if (Mapper::context_remapping &&
                    order > m_remapping_order + 1) {
                    token_id.first =
                        m_mapper.map_id(gram, token_id.first, &m_vocab,
                                        &m_arrays.front(), m_remapping_order,
                                        true);  // BACKWARD trie
                }

                sa_builder.add_gram(token_id.first);

                uint64_t prob_rank =
                    probs_builder.rank(order - 2, record.prob, 0);
                uint64_t rank = prob_rank;
                if (order != m_order) {
                    uint64_t backoff_rank =
                        backoffs_builder.rank(order - 2, record.backoff, 1);
                    // store probs/backoffs ranks interleaved
                    rank |= backoff_rank << probs_quantization_bits;
                }

                sa_builder.add_prob_backoff_rank(rank);
            }

            // set remaining pointers (if any)
            for (; j != prev_n; ++j) {
                pointers.push_back(pos);
            }
        }
    };

    trie_prob_lm()
        : m_order(0)
        , m_remapping_order(0)
        , m_unk_prob(global::default_unk_prob) {}

    typedef prob_model_state<uint64_t> state_type;

    state_type state() {
        return state_type(order());
    }

    float score(state_type& state, byte_range const word, bool& is_OOV);

    inline uint64_t order() const {
        return uint64_t(m_order);
    }

    uint64_t remapping_order() const {
        return uint64_t(m_remapping_order);
    }

    void print_stats(size_t bytes) const;

    uint64_t size() const {
        uint64_t size = 0;
        for (auto const& a : m_arrays) {
            size += a.size();
        }
        return size;
    }

    void save(std::ostream& os) const {
        essentials::save_pod(os, m_order);
        essentials::save_pod(os, m_remapping_order);
        essentials::save_pod(os, m_unk_prob);
        m_probs_averages.save(os);
        m_backoffs_averages.save(os);
        m_vocab.save(os);
        m_arrays.front().save(os, 1, value_type::none);
        for (uint8_t order_m1 = 1; order_m1 < m_order; ++order_m1) {
            m_arrays[order_m1].save(os, order_m1 + 1, value_type::prob_backoff);
        }
    }

    void load(std::istream& is) {
        essentials::load_pod(is, m_order);
        essentials::load_pod(is, m_remapping_order);
        essentials::load_pod(is, m_unk_prob);
        m_probs_averages.load(is, m_order - 1);
        m_backoffs_averages.load(is, m_order - 2);
        m_vocab.load(is);
        m_arrays.resize(m_order);
        m_arrays.front().load(is, 1, value_type::none);
        for (uint8_t order_m1 = 1; order_m1 < m_order; ++order_m1) {
            m_arrays[order_m1].load(is, order_m1 + 1, value_type::prob_backoff);
        }
    }

private:
    uint8_t m_order;
    uint8_t m_remapping_order;
    float m_unk_prob;
    Mapper m_mapper;
    Values m_probs_averages;
    Values m_backoffs_averages;
    Vocabulary m_vocab;
    std::vector<sorted_array_type> m_arrays;
};

}  // namespace tongrams
