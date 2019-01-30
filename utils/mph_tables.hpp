#pragma once

#include <vector>
#include <numeric>

#include "mphf.hpp"
#include "util.hpp"

#include "../emphf/common.hpp"
#include "../emphf/base_hash.hpp"
#include "../emphf/mmap_memory_model.hpp"
#include "../emphf/hypergraph_sorter_scan.hpp"

#include "../vectors/compact_vector.hpp"
#include "../vectors/hash_compact_vector.hpp"
#include "../vectors/compact_triplets_vector.hpp"

namespace tongrams {

    template<typename KeyValueSequence,
             typename BaseHasher>
    struct single_valued_mpht {
        typedef tongrams::mphf<BaseHasher> hash_function;

        struct builder {
            builder()
            {}

            template<typename T, typename Adaptor>
            builder(std::vector<T> const& ngrams,
                    compact_vector const& values,
                    Adaptor adaptor)
                : m_data_builder(values.size(), values.width())
            {
                assert(ngrams.size() == values.size());
                using namespace emphf;
                size_t n = ngrams.size();

                typedef hypergraph_sorter_scan<uint32_t, mmap_memory_model> hs32_t;
                typedef hypergraph_sorter_scan<uint64_t, mmap_memory_model> hs64_t;
                size_t max_nodes = (size_t(std::ceil(double(n) * 1.23)) + 2) / 3 * 3;
                if (max_nodes >= uint64_t(1) << 32) {
                    hs64_t sorter;
                    hash_function(sorter, n, ngrams, adaptor).swap(m_h);
                } else {
                    hs32_t sorter;
                    hash_function(sorter, n, ngrams, adaptor).swap(m_h);
                }

                auto it = ngrams.begin();
                for (auto value: values) {
                    auto hashes = m_h.hashes(*it, adaptor);
                    uint64_t key = m_h.mix_hashes(hashes);
                    uint64_t pos = m_h.lookup(hashes);
                    m_data_builder.set(pos, key, value);
                    ++it;
                }
            }

            void build(single_valued_mpht<KeyValueSequence,
                                          BaseHasher>& mpht)
            {
                mpht.m_h.swap(m_h);
                mpht.m_data.build(m_data_builder);
                builder().swap(*this);
            }

            void swap(builder& other) {
                m_h.swap(other.m_h);
                m_data_builder.swap(other.m_data_builder);
            }

        private:
            hash_function m_h;
            typename KeyValueSequence::builder m_data_builder;
        };

        single_valued_mpht()
        {}

        single_valued_mpht(single_valued_mpht::builder& in) {
            in.build(*this);
            single_valued_mpht::builder().swap(in);
        }

        template<typename T, typename Adaptor>
        uint64_t lookup(T gram, Adaptor adaptor) const {
            auto hashes = m_h.hashes(gram, adaptor);
            uint64_t key = m_h.mix_hashes(hashes);
            uint64_t pos = m_h.lookup(hashes);
            auto const& p = m_data[pos];
            return p.first == key ? p.second : global::not_found;
        }

        size_t size() const {
            return m_h.size();
        }

        hash_function const& h() const {
            return m_h;
        }

        size_t data_bytes() const {
            return m_data.bytes();
        }

        void swap(single_valued_mpht& other) {
            m_h.swap(other.m_h);
            m_data.swap(other.m_data);
        }

        void save(std::ostream& os) const {
            m_h.save(os);
            m_data.save(os);
        }

        void load(std::istream& is) {
            m_h.load(is);
            m_data.load(is);
        }

    private:
        hash_function m_h;
        KeyValueSequence m_data;
    };

    typedef single_valued_mpht<
                               hash_compact_vector<uint32_t>,
                               emphf::jenkins32_hasher
                              >
                              single_valued_mpht32;

    typedef single_valued_mpht<
                               hash_compact_vector<uint64_t>,
                               emphf::jenkins64_hasher
                              >
                              single_valued_mpht64;

    template<typename BaseHasher>
    struct double_valued_mpht {

        typedef tongrams::mphf<BaseHasher> hash_function;

        double_valued_mpht()
        {}

        template <typename T, typename Adaptor>
        double_valued_mpht(std::vector<T> const& ngrams,
                           compact_vector const& keys,
                           compact_vector const& values1,
                           compact_vector const& values2,
                           Adaptor adaptor)
        {
            build(ngrams, keys, values1, values2, adaptor);
        }

        template <typename T, typename Adaptor>
        void build(std::vector<T> const& ngrams,
                   compact_vector const& keys,
                   compact_vector const& values1,
                   compact_vector const& values2,
                   Adaptor adaptor)
        {
            assert(ngrams.size() == values1.size());
            assert(ngrams.size() == values2.size());

            using namespace emphf;
            size_t n = ngrams.size();

            typedef hypergraph_sorter_scan<uint32_t, mmap_memory_model> hs32_t;
            typedef hypergraph_sorter_scan<uint64_t, mmap_memory_model> hs64_t;
            size_t max_nodes = (size_t(std::ceil(double(n) * 1.23)) + 2) / 3 * 3;
            if (max_nodes >= uint64_t(1) << 32) {
                hs64_t sorter;
                hash_function(sorter, n, ngrams, adaptor).swap(m_h);
            } else {
                hs32_t sorter;
                hash_function(sorter, n, ngrams, adaptor).swap(m_h);
            }


            compact_triplets_vector::builder data_builder(n,
                                                          keys.size() ? keys.width() : 64,
                                                          values1.width(),
                                                          values2.width());
            auto it_keys = keys.begin();
            auto it_values1 = values1.begin();
            auto it_values2 = values2.begin();
            for (auto gram: ngrams) {
                auto hashes = m_h.hashes(gram, adaptor);
                // NOTE:
                // use 64-bit hashes as keys if not provided by client
                // trie_prob_lm uses default 64-bit hashes
                uint64_t key = keys.size() ? *it_keys : m_h.mix_hashes(hashes);
                uint64_t pos = m_h.lookup(hashes);
                uint64_t value = *it_values2;
                compact_triplets_vector::value_type triple(key, *it_values1, value);
                data_builder.set(pos, triple);
                ++it_keys;
                ++it_values1;
                ++it_values2;
            }
            data_builder.build(m_data);
        }

        // always compare with default hash keys
        template<typename T, typename Adaptor>
        uint64_pair lookup_pair(T gram, Adaptor adaptor) const {
            auto hashes = m_h.hashes(gram, adaptor);
            uint64_t key = m_h.mix_hashes(hashes);
            uint64_t pos = m_h.lookup(hashes);
            auto const& triple = m_data[pos];
            uint64_pair values;
            if (std::get<0>(triple) == key) {
                values.first = std::get<1>(triple);
                values.second = std::get<2>(triple);
            } else {
                // if key not found, just assign default values
                values.first = global::not_found;
                values.second = global::not_found;
            }
            return values;
        }

        // compare with the passed key
        // template<typename T, typename Adaptor>
        // uint64_pair lookup_pair(T gram, uint64_t key, Adaptor adaptor) const {
        //     auto hashes = m_h.hashes(gram, adaptor);
        //     uint64_t pos = m_h.lookup(hashes);
        //     auto const& triple = m_data[pos];
        //     uint64_pair values;
        //     if (std::get<0>(triple) == key) {
        //         values.first = std::get<1>(triple);
        //         values.second = std::get<2>(triple);
        //     } else {
        //         // if key not found, just assign default values
        //         values.first = global::not_found;
        //         values.second = global::not_found;
        //     }
        //     return values;
        // }

        // just return the first field, comparing with default hash key
        template<typename T, typename Adaptor>
        uint64_t lookup(T gram, Adaptor adaptor) const {
            auto hashes = m_h.hashes(gram, adaptor);
            uint64_t key = m_h.mix_hashes(hashes);
            uint64_t pos = m_h.lookup(hashes);
            auto const& triple = m_data[pos];
            return std::get<0>(triple) == key
                 ? std::get<1>(triple) : global::not_found;
        }

        size_t size() const {
            return m_h.size();
        }

        hash_function const& h() const {
            return m_h;
        }

        size_t data_bytes() const {
            return m_data.bytes();
        }

        void swap(double_valued_mpht& other) {
            m_h.swap(other.m_h);
            m_data.swap(other.m_data);
        }

        void save(std::ostream& os) const {
            m_h.save(os);
            m_data.save(os);
        }

        void load(std::istream& is) {
            m_h.load(is);
            m_data.load(is);
        }

    private:
        hash_function m_h;
        compact_triplets_vector m_data;
    };

    typedef double_valued_mpht<emphf::jenkins32_hasher> double_valued_mpht32;
    typedef double_valued_mpht<emphf::jenkins64_hasher> double_valued_mpht64;

    // NOTE: used by sequences/fast_ef_sequence.hpp
    template<typename UintValueType1,
             typename UintValueType2,
             typename BaseHasher = emphf::jenkins64_hasher>
    struct uint_mpht {

        typedef tongrams::mphf<BaseHasher> hash_function;

        uint_mpht()
        {}

        template<typename Adaptor>
        uint_mpht(std::vector<UintValueType1> const& from,
                 std::vector<UintValueType2> const& to,
                 Adaptor adaptor)
        {
            build(from, to, adaptor);
        }

        template<typename Adaptor>
        void build(std::vector<UintValueType1> const& from,
                   std::vector<UintValueType2> const& to,
                   Adaptor adaptor)
        {
            using namespace emphf;
            size_t n = from.size();

            typedef hypergraph_sorter_scan<uint32_t, mmap_memory_model> hs32_t;
            typedef hypergraph_sorter_scan<uint64_t, mmap_memory_model> hs64_t;
            size_t max_nodes = (size_t(std::ceil(double(n) * 1.23)) + 2) / 3 * 3;
            if (max_nodes >= uint64_t(1) << 32) {
                hs64_t sorter;
                hash_function(sorter, n, from, adaptor).swap(m_h);
            } else {
                hs32_t sorter;
                hash_function(sorter, n, from, adaptor).swap(m_h);
            }

            auto it = from.begin();
            compact_vector::builder cvb(to.size(), util::ceil_log2(to.back() + 1));
            for (auto t: to) {
                auto hashes = m_h.hashes(*it, adaptor);
                uint64_t pos = m_h.lookup(hashes);
                cvb.set(pos, t);
                ++it;
            }
            m_values.build(cvb);
        }

        template<typename Adaptor>
        UintValueType2 lookup(UintValueType1 x,
                              Adaptor adaptor) const {
            auto hashes = m_h.hashes(x, adaptor);
            uint64_t pos = m_h.lookup(hashes);
            return m_values[pos];
        }

        size_t size() const {
            return m_h.size();
        }

        hash_function const& h() const {
            return m_h;
        }

        // does NOT include bytes
        // for hash function itself
        size_t data_bytes() const {
            return m_values.bytes();
        }

        void swap(uint_mpht& other) {
            m_h.swap(other.m_h);
            m_values.swap(other.m_values);
        }

        void save(std::ostream& os) const {
            m_h.save(os);
            m_values.save(os);
        }

        void load(std::istream& is) {
            m_h.load(is);
            m_values.load(is);
        }

    private:
        hash_function m_h;
        compact_vector m_values;
    };
}
