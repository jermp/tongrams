#pragma once

#include <cassert>

#include "../utils/util.hpp"
#include "bit_vector.hpp"

namespace tongrams
{
    struct compact_triplets_vector
    {
        typedef std::tuple<uint64_t, uint64_t, uint64_t> value_type;

        struct builder
        {
            builder()
                : m_size(0)
                , m_width1(0)
                , m_width2(0)
                , m_width3(0)
            {}

            builder(uint64_t n, uint64_t w1, uint64_t w2, uint64_t w3)
                : m_size(n)
                , m_width1(width(w1))
                , m_width2(width(w2))
                , m_width3(width(w3))
                , m_bits(m_size * (m_width1 + m_width2 + m_width3), 0)
            {
                check_width(m_width1);
                check_width(m_width2);
                check_width(m_width3);
            }

            uint64_t width(uint64_t w) {
                return !w ? w + 1 : w;
            }

            uint64_t mask(uint64_t w) {
                return -(w == 64) | ((uint64_t(1) << w) - 1);
            }

            void check_width(uint64_t w) {
                if (w > 64) {
                    std::cerr << "Error: width must be <= 64, "
                              << "got " << w << std::endl;
                    std::terminate();
                }
            }

            inline void set(uint64_t i, value_type triplet) {
                uint64_t pos = i * (m_width1 + m_width2 + m_width3);
                m_bits.set_bits(pos, std::get<0>(triplet), m_width1);
                pos += m_width1;
                m_bits.set_bits(pos, std::get<1>(triplet), m_width2);
                pos += m_width2;
                m_bits.set_bits(pos, std::get<2>(triplet), m_width3);
            }

            void build(compact_triplets_vector& in) {
                in.m_size = m_size;
                in.m_width1 = m_width1;
                in.m_width2 = m_width2;
                in.m_width3 = m_width3;
                in.m_bits.build(&m_bits);
                builder().swap(*this);
            }

            void swap(compact_triplets_vector::builder& other) {
                std::swap(m_size, other.m_size);
                std::swap(m_width1, other.m_width1);
                std::swap(m_width2, other.m_width2);
                std::swap(m_width3, other.m_width3);
                m_bits.swap(other.m_bits);
            }

            uint64_t size() const {
                return m_size;
            }

        private:
            uint64_t m_size;
            uint64_t m_width1;
            uint64_t m_width2;
            uint64_t m_width3;
            bit_vector_builder m_bits;
        };

        compact_triplets_vector()
        {}

        inline value_type operator[](uint64_t pos) const {
            assert(pos < m_size);
            value_type triplet;
            bits_iterator<bit_vector> it(m_bits, pos * (m_width1 + m_width2 + m_width3));
            std::get<0>(triplet) = it.get_bits(m_width1);
            std::get<1>(triplet) = it.get_bits(m_width2);
            std::get<2>(triplet) = it.get_bits(m_width3);
            return triplet;
        }

        size_t bytes() const {
            return sizeof(m_size)
                 + sizeof(m_width1)
                 + sizeof(m_width2)
                 + sizeof(m_width3)
                 + m_bits.bytes();
        }

        void swap(compact_triplets_vector& other) {
            std::swap(m_size, other.m_size);
            std::swap(m_width1, other.m_width1);
            std::swap(m_width2, other.m_width2);
            std::swap(m_width3, other.m_width3);
            m_bits.swap(other.m_bits);
        }

        void save(std::ostream& os) const {
            util::save_pod(os, &m_size);
            util::save_pod(os, &m_width1);
            util::save_pod(os, &m_width2);
            util::save_pod(os, &m_width3);
            m_bits.save(os);
        }

        void load(std::istream& is) {
            util::load_pod(is, &m_size);
            util::load_pod(is, &m_width1);
            util::load_pod(is, &m_width2);
            util::load_pod(is, &m_width3);
            m_bits.load(is);
        }

    private:
        uint64_t m_size;
        uint64_t m_width1;
        uint64_t m_width2;
        uint64_t m_width3;
        bit_vector m_bits;
    };
}
