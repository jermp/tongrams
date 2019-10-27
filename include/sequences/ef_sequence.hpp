#pragma once

#include "sequences/darray.hpp"
#include "vectors/bit_vector.hpp"

namespace tongrams {

struct ef_sequence {
    ef_sequence() : m_l(0), m_size(0) {}

    template <typename Iterator>
    ef_sequence(Iterator begin, uint64_t n, uint64_t u,
                bool index_on_zeros = false) {
        build(begin, n, u, index_on_zeros);
    }

    template <typename Iterator>
    void build(Iterator begin, uint64_t n, uint64_t u,
               bool index_on_zeros = false) {
        m_size = n;
        m_l = uint8_t((n && u / n) ? util::msb(u / n) : 0);
        bit_vector_builder bvb_high_bits(n + (u >> m_l) + 1);
        bit_vector_builder bvb_low_bits;
        bvb_low_bits.reserve(n * m_l);

        uint64_t low_mask = (uint64_t(1) << m_l) - 1;
        uint64_t last = 0;
        for (size_t i = 0; i < n; ++i, ++begin) {
            auto v = *begin;
            if (i && v < last) {
                std::cout << "at pos: " << i << "/" << n << std::endl;
                std::cout << "v = " << v << "; last = " << last << std::endl;
                throw std::runtime_error("sequence is not sorted.");
            }
            if (m_l) {
                bvb_low_bits.append_bits(v & low_mask, m_l);
            }
            bvb_high_bits.set((v >> m_l) + i, 1);
            last = v;
        }

        bit_vector(&bvb_high_bits).swap(m_high_bits);
        bit_vector(&bvb_low_bits).swap(m_low_bits);
        darray1(m_high_bits).swap(m_high_bits_d1);
        if (index_on_zeros) {
            darray0(m_high_bits).swap(m_high_bits_d0);
        }
    }

    inline uint64_t operator[](uint64_t i) {
        assert(i < size());
        return ((m_high_bits_d1.select(m_high_bits, i) - i) << m_l) |
               m_low_bits.get_bits(i * m_l, m_l);
    }

    inline uint64_t num_ones() const {
        return m_high_bits_d1.num_positions();
    }

    struct iterator {
        iterator(ef_sequence const& ef, uint64_t i = 0)
            : m_ef(&ef), m_i(i), m_l(ef.m_l) {
            m_low_mask = (uint64_t(1) << m_l) - 1;
            m_low_buf = 0;
            if (m_l) {
                m_chunks_in_word = 64 / m_l;
                m_chunks_avail = 0;
            } else {
                m_chunks_in_word = 0;
                m_chunks_avail = m_ef->num_ones();
            }

            if (!m_ef->num_ones())
                return;
            uint64_t pos = m_ef->m_high_bits_d1.select(m_ef->m_high_bits, m_i);
            m_high_enum = bit_vector::unary_iterator(m_ef->m_high_bits, pos);
            assert(m_l < 64);
        }

        uint64_t next() {
            if (!m_chunks_avail--) {
                m_low_buf = m_ef->m_low_bits.get_word64(m_i * m_l);
                m_chunks_avail = m_chunks_in_word - 1;
            }

            uint64_t high = m_high_enum.next();
            assert(high == m_ef->m_high_bits_d1.select(m_ef->m_high_bits, m_i));
            uint64_t low = m_low_buf & m_low_mask;
            uint64_t ret = (((high - m_i) << m_l) | low);
            ++m_i;
            m_low_buf >>= m_l;

            return ret;
        }

    private:
        ef_sequence const* m_ef;
        uint64_t m_i;
        uint64_t m_l;
        bit_vector::unary_iterator m_high_enum;
        uint64_t m_low_buf;
        uint64_t m_low_mask;
        uint64_t m_chunks_in_word;
        uint64_t m_chunks_avail;
    };

    inline uint64_pair pair(uint64_t i) {
        return {operator[](i), operator[](i + 1)};
    }

    uint64_t size() const {
        return m_size;
    }

    uint64_t universe() {
        return operator[](m_size - 1);
    }

    uint64_t bytes() const {
        return sizeof(m_l) + m_high_bits.bytes() + m_high_bits_d1.bytes() +
               m_high_bits_d0.bytes() + m_low_bits.bytes() + sizeof(m_size);
    }

    void swap(ef_sequence& other) {
        std::swap(other.m_size, m_size);
        other.m_high_bits.swap(m_high_bits);
        other.m_high_bits_d1.swap(m_high_bits_d1);
        other.m_high_bits_d0.swap(m_high_bits_d0);
        other.m_low_bits.swap(m_low_bits);
        std::swap(other.m_l, m_l);
    }

    void save(std::ostream& os) const {
        util::save_pod(os, &m_l);
        m_high_bits.save(os);
        m_high_bits_d1.save(os);
        m_high_bits_d0.save(os);
        m_low_bits.save(os);
        util::save_pod(os, &m_size);
    }

    void load(std::istream& is) {
        util::load_pod(is, &m_l);
        m_high_bits.load(is);
        m_high_bits_d1.load(is);
        m_high_bits_d0.load(is);
        m_low_bits.load(is);
        util::load_pod(is, &m_size);
    }

private:
    uint8_t m_l;
    bit_vector m_high_bits;
    darray1 m_high_bits_d1;
    darray0 m_high_bits_d0;
    bit_vector m_low_bits;
    uint64_t m_size;
};
}  // namespace tongrams
