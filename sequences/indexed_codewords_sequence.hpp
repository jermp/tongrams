#pragma once

#include "../utils/util.hpp"
#include "../vectors/bit_vector.hpp"
#include "darray.hpp"

namespace tongrams
{
    struct indexed_codewords_sequence
    {
        indexed_codewords_sequence()
            : m_size(0)
        {}

        template<typename Iterator>
        void build(Iterator begin, uint64_t n, uint8_t /*order*/)
        {
            m_size = n;
            auto start = begin;
            uint64_t bits = 0;
            for (uint64_t i = 0; i < n; ++i, ++start) {
                bits += util::floor_log2(*start + 2);
            }
            bit_vector_builder bvb_index(bits + 1);
            bit_vector_builder bvb_codewords(bits);

            uint64_t pos = 0;
            for (uint64_t i = 0; i < n; ++i, ++begin) {
                auto v = *begin;
                uint64_t len = util::floor_log2(v + 2);
                assert(len <= 64);
                uint64_t cw = v + 2 - (uint64_t(1) << len);
                bvb_codewords.set_bits(pos, cw, len);
                bvb_index.set(pos, 1);
                pos += len;
            }
            // NOTE: store a last 1 to delimit last codeword:
            // avoid test for last codeword in operator[]
            assert(pos == bits);
            bvb_index.set(pos, 1);

            bit_vector(&bvb_codewords).swap(m_codewords);
            bit_vector(&bvb_index).swap(m_index);
            darray1(m_index).swap(m_index_d1);
        }

        inline uint64_t operator[](uint64_t i) const {
            uint64_t pos = m_index_d1.select(m_index, i);
            assert(pos < m_index.size());
            bit_vector::unary_iterator e(m_index, pos + 1);
            uint64_t next = e.next();
            uint64_t len = next - pos;
            assert(len <= 64);
            uint64_t cw = m_codewords.get_bits(pos, len);
            uint64_t value = cw + (uint64_t(1) << len) - 2;
            return value;
        }

        struct iterator
        {
            iterator(indexed_codewords_sequence const* seq, uint64_t pos)
                : m_seq(seq)
                , m_pos(pos)
                , m_len(0)
            {
                m_it = bit_vector::unary_iterator(m_seq->m_index, pos);
                m_pos = m_it.next();
            }

            uint64_t operator*() {
                m_len = m_it.next() - m_pos;
                assert(m_len);
                uint64_t cw = m_seq->m_codewords.get_bits(m_pos, m_len);
                return cw + (uint64_t(1) << m_len) - 2;
            }

            void operator++() {
                m_pos += m_len;
            }

            bool operator==(iterator const& other) const {
                return m_pos == other.m_pos;
            }

            bool operator!=(iterator const& other) const {
                return !(*this == other);
            }

        private:
            indexed_codewords_sequence const* m_seq;
            uint64_t m_pos;
            uint64_t m_len;
            bit_vector::unary_iterator m_it;
        };

        iterator begin() const {
            return iterator(this, 0);
        }

        iterator end() const {
            return iterator(this, m_codewords.size());
        }

        uint64_t size() const {
            return m_size;
        }

        void save(std::ostream& os) const {
            util::save_pod(os, &m_size);
            m_codewords.save(os);
            m_index.save(os);
            m_index_d1.save(os);
        }

        void load(std::istream& is) {
            util::load_pod(is, &m_size);
            m_codewords.load(is);
            m_index.load(is);
            m_index_d1.load(is);
        }

        uint64_t bytes() const {
            return sizeof(m_size)
                 + m_codewords.bytes()
                 + m_index.bytes()
                 + m_index_d1.bytes();
        }

    private:
        uint64_t m_size;
        bit_vector m_codewords;
        bit_vector m_index;
        darray1 m_index_d1;
    };
}
