#pragma once

#include <iostream>
#include <vector>
#include <cassert>

#include "../utils/util.hpp"

namespace tongrams
{
    template<typename HashType>
    struct hash_compact_vector
    {
        typedef HashType hash_t;
        typedef std::pair<hash_t, uint64_t> key_value_pair;
        static const uint64_t hash_bits = sizeof(hash_t) * 8;

        struct builder
        {
            builder()
                : m_size(0)
                , m_width(0)
                , m_mask(0)
            {}

            builder(uint64_t n, uint64_t w)
                : m_size(n)
                , m_width(w)
                , m_mask(-(w == hash_bits || w == 64) | ((uint64_t(1) << w) - 1))
                , m_bits(util::words_for<hash_t>(m_size * (hash_bits + m_width)), 0)
            {
                if (!m_width) {
                    std::cerr << "Error: value width must be non zero." << std::endl;
                    std::abort();
                }

                if (m_width > 64) {
                    std::cerr << "Error: value width must be <= 64." << std::endl;
                    std::abort();
                }
            }

            void set(uint64_t i, hash_t k, uint64_t v)
            {
                assert(i < m_size);
                hash_t pos = i * (hash_bits + m_width);
                hash_t block = pos / hash_bits;
                hash_t shift = pos & (hash_bits - 1);

                if (shift) {
                    // writing of key
                    uint64_t hash_key_mask = (hash_t(1) << shift) - 1;
                    m_bits[block] &= hash_key_mask;
                    m_bits[block] |= k << shift;
                    hash_t res_shift = hash_bits - shift;           
                    m_bits[block + 1] &= ~hash_key_mask;
                    m_bits[block + 1] |= k >> res_shift;

                    // writing of value
                    m_bits[block + 1] &= ~(m_mask << shift);
                    m_bits[block + 1] |= v << shift;

                    if (m_width > res_shift) {
                        m_bits[block + 2] &= ~(m_mask >> res_shift);
                        m_bits[block + 2] |= v >> res_shift;
                        if (hash_bits == 32 && m_width > res_shift + 32) {
                            // NOTE (1): this case is never happening when hash_bits == 64
                            uint64_t res_shift2 = m_width - 32 - res_shift;
                            m_bits[block + 3] &= ~(m_mask >> res_shift2);
                            m_bits[block + 3] |= v >> (32 + res_shift);
                        }
                    }
                } else {
                    m_bits[block] = k;
                    m_bits[block + 1] &= ~m_mask;
                    m_bits[block + 1] |= v;
                    if (hash_bits == 32 && m_width > 32) { // see NOTE (1)
                        uint64_t res_shift = m_width - 32;
                        m_bits[block + 2] &= ~(m_mask >> res_shift);
                        m_bits[block + 2] |= v >> 32;
                    }
                }
            }

            void swap(hash_compact_vector::builder& other) {
                std::swap(m_size, other.m_size);
                std::swap(m_width, other.m_width);
                std::swap(m_mask, other.m_mask);
                m_bits.swap(other.m_bits);
            }

            uint64_t size() const {
                return m_size;
            }

            uint64_t width() const {
                return m_width;
            }

            uint64_t mask() const {
                return m_mask;
            }

            std::vector<hash_t>& bits() {
                return m_bits;
            }

        private:
            uint64_t m_size;
            uint64_t m_width;
            uint64_t m_mask;
            std::vector<hash_t> m_bits;
        };

        hash_compact_vector()
        {}

        hash_compact_vector(hash_compact_vector::builder& in) {
            build(in);
        }

        void build(hash_compact_vector::builder& in) {
            m_size = in.size();
            m_width = in.width();
            m_mask = in.mask();
            m_bits.swap(in.bits());
            hash_compact_vector::builder().swap(in);
        }

        key_value_pair operator[](uint64_t i) const
        {
            hash_t pos = i * (hash_bits + m_width);
            hash_t block = pos / hash_bits;
            hash_t shift = pos & (hash_bits - 1);

            hash_t k = 0;
            uint64_t v = 0;
            if (shift) {
                hash_t res_shift = hash_bits - shift;
                k = (m_bits[block] >> shift)
                  | (m_bits[block + 1] << res_shift);
                v = (m_bits[block + 1] >> shift) & m_mask;

                if (m_width > res_shift) {
                    uint64_t res_shift2 = m_width - res_shift;
                    v |= (m_bits[block + 2] & ((uint64_t(1) << res_shift2) - 1)) << res_shift;
                    if (hash_bits == 32 && m_width > res_shift + 32) { // see NOTE (1)
                        uint64_t res_shift2 = m_width - 32 - res_shift;
                        v |= (m_bits[block + 3] & ((uint64_t(1) << res_shift2) - 1))
                          << (32 + res_shift);
                    }
                }
            } else {
                k = m_bits[block];
                v = m_bits[block + 1] & m_mask;
                if (hash_bits == 32 && m_width > 32) { // see NOTE (1)
                    uint64_t res_shift = m_width - 32;
                    v |= (m_bits[block + 2] & ((uint64_t(1) << res_shift) - 1)) << 32;
                }
            }

            return {k, v};
        }

        size_t bytes() const {
            return sizeof(m_size)
                 + sizeof(m_width)
                 + sizeof(m_mask)
                 + m_bits.size() * sizeof(hash_t);
        }

        void swap(hash_compact_vector& other) {
            std::swap(m_size, other.m_size);
            std::swap(m_width, other.m_width);
            std::swap(m_mask, other.m_mask);
            m_bits.swap(other.m_bits);
        }

        void save(std::ostream& os) const {
            util::save_pod(os, &m_size);
            util::save_pod(os, &m_width);
            util::save_pod(os, &m_mask);
            util::save_vec(os, m_bits);
        }

        void load(std::istream& is) {
            util::load_pod(is, &m_size);
            util::load_pod(is, &m_width);
            util::load_pod(is, &m_mask);
            util::load_vec(is, m_bits);
        }

    private:
        uint64_t m_size;
        uint64_t m_width;
        uint64_t m_mask;
        std::vector<hash_t> m_bits;
    };
}
