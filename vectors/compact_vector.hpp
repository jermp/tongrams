#pragma once

#include "../utils/util.hpp"

namespace tongrams
{
    struct compact_vector
    {
        template<typename Data>
        struct iterator
            : public std::iterator<std::input_iterator_tag, uint64_t>
        {
            iterator(Data const* data, uint64_t i = 0)
                : m_i(i)
                , m_cur_val(0)
                , m_cur_block((i * data->m_width) >> 6)
                , m_cur_shift((i * data->m_width) & 63)
                , m_data(data)
            {}

            uint64_t operator*() {
                read();
                return m_cur_val;
            }

            iterator& operator++() {
                ++m_i;
                return *this;
            }

            bool operator==(iterator const& other) const {
                return m_i == other.m_i;
            }

            bool operator!=(iterator const& other) const {
                return !(*this == other);
            }

        private:
            uint64_t m_i;
            uint64_t m_cur_val;
            uint64_t m_cur_block;
            int64_t m_cur_shift;
            Data const* m_data;

            void read()
            {
                if (m_cur_shift + m_data->m_width <= 64) {
                    m_cur_val = m_data->m_bits[m_cur_block]
                             >> m_cur_shift & m_data->m_mask;
                } else {
                    uint64_t res_shift = 64 - m_cur_shift;
                    m_cur_val = (m_data->m_bits[m_cur_block] >> m_cur_shift)
                              | (m_data->m_bits[m_cur_block + 1] << res_shift & m_data->m_mask);
                    ++m_cur_block;
                    m_cur_shift = -res_shift;
                }

                m_cur_shift += m_data->m_width;

                if (m_cur_shift == 64) {
                    m_cur_shift = 0;
                    ++m_cur_block;
                }
            }
        };

        struct builder
        {
            builder(uint64_t n = 0, uint64_t w = 0)
                : m_size(n)
                , m_width(!w ? w + 1 : w)
                , m_mask(-(w == 64) | ((1ULL << w) - 1))
                , m_back(0)
                , m_cur_block(0)
                , m_cur_shift(0)
                , m_bits(util::words_for(m_size * m_width), 0)
            {
                if (m_width > 64) {
                    std::cerr << "Error: width must be <= 64."
                              << std::endl;
                    std::terminate();
                }
            }

            template<typename Iterator>
            builder(Iterator begin, uint64_t n, uint64_t w)
                : builder(n, w)
            {
                fill(begin, n);
            }

            template<typename Iterator>
            void fill(Iterator begin, uint64_t n)
            {
                if (!m_width) {
                    throw std::runtime_error("width must be greater than 0");
                }

                for (uint64_t i = 0; i < n; ++i, ++begin) {
                    push_back(*begin);
                }
            }

            inline void set(uint64_t i, uint64_t v)
            {
                assert(m_width);
                assert(i < m_size);
                if (i == m_size - 1) {
                    m_back = v;
                }

                uint64_t pos = i * m_width;
                uint64_t block = pos >> 6;
                uint64_t shift = pos & 63;

                m_bits[block] &= ~(m_mask << shift);
                m_bits[block] |= v << shift;

                uint64_t res_shift = 64 - shift;
                if (res_shift < m_width) {
                    m_bits[block + 1] &= ~(m_mask >> res_shift);
                    m_bits[block + 1] |= v >> res_shift;
                }
            }

            inline void push_back(uint64_t v)
            {
                assert(m_width);
                m_back = v;
                m_bits[m_cur_block] &= ~(m_mask << m_cur_shift);
                m_bits[m_cur_block] |= v << m_cur_shift;

                uint64_t res_shift = 64 - m_cur_shift;
                if (res_shift < m_width) {
                    ++m_cur_block;
                    m_bits[m_cur_block] &= ~(m_mask >> res_shift);
                    m_bits[m_cur_block] |= v >> res_shift;
                    m_cur_shift = -res_shift;
                }

                m_cur_shift += m_width;

                if (m_cur_shift == 64) {
                    m_cur_shift = 0;
                    ++m_cur_block;
                }
            }

            friend struct iterator<builder>;

            typedef iterator<builder> iterator_type;

            iterator_type begin() const {
                return iterator_type(this); 
            }

            iterator_type end() const {
                return iterator_type(this, m_size);
            }

            void swap(compact_vector::builder& other)
            {
                std::swap(m_size, other.m_size);
                std::swap(m_width, other.m_width);
                std::swap(m_mask, other.m_mask);
                std::swap(m_cur_block, other.m_cur_block);
                std::swap(m_cur_shift, other.m_cur_shift);
                m_bits.swap(other.m_bits);
            }

            uint64_t back() const {
                return m_back;
            }

            uint64_t size() const {
                return m_size;
            }

            uint64_t width() const {
                return m_width;
            }

            std::vector<uint64_t>& bits() {
                return m_bits;
            }

        private:
            uint64_t m_size;
            uint64_t m_width;
            uint64_t m_mask;
            uint64_t m_back;
            uint64_t m_cur_block;
            int64_t m_cur_shift;
            std::vector<uint64_t> m_bits;
        };

        compact_vector()
            : m_size(0)
            , m_width(0)
            , m_mask(0)
        {}

        compact_vector(compact_vector::builder& in) {
            build(in);
        }

        void build(compact_vector::builder& in) {
            m_size = in.size();
            m_width = in.width();
            m_mask = -(m_width == 64) | ((1ULL << m_width) - 1);
            m_bits.swap(in.bits());
        }

        template<typename Iterator>
        void build(Iterator begin, size_t n, uint8_t /*order*/) {
            uint64_t max = 0;
            auto it = begin;
            for (uint64_t i = 0; i < n; ++i, ++it) {
                uint64_t v = *it;
                if (v > max) {
                    max = v;
                }
            }
            compact_vector::builder cvb(n, util::ceil_log2(max + 1));
            cvb.fill(begin, n);
            build(cvb);
        }

        inline uint64_t operator[](uint64_t i) const {
            assert(i < size());
            uint64_t pos = i * m_width;
            uint64_t block = pos >> 6;
            uint64_t shift = pos & 63;
            return shift + m_width <= 64
                 ? m_bits[block] >> shift & m_mask
                 :(m_bits[block] >> shift) | (m_bits[block + 1] << (64 - shift) & m_mask);
        }

        // it retrieves at least 57 bits
        inline uint64_t access(uint64_t i) const {
            assert(i < size());
            uint64_t pos = i * m_width;
            const char* ptr = reinterpret_cast<const char*>(m_bits.data());
            return (*(reinterpret_cast<uint64_t const*>(ptr + (pos >> 3))) >> (pos & 7)) & m_mask;
        }

        inline void prefetch(size_t i) const {
            util::prefetch(m_bits.data() + i);
        }

        uint64_t back() const {
            return operator[](size() - 1);
        }

        uint64_t size() const {
            return m_size;
        }

        uint64_t width() const {
            return m_width;
        }

        typedef iterator<compact_vector> iterator_type;

        iterator_type begin(uint64_t pos = 0) const {
            return iterator_type(this, pos);
        }

        iterator_type end() const {
            return iterator_type(this, m_size);
        }

        std::vector<uint64_t> const& bits() const {
            return m_bits;
        }

        size_t bytes() const {
            return m_bits.size() * 8
                 + sizeof(m_size)
                 + sizeof(m_width)
                 + sizeof(m_mask);
        }

        void swap(compact_vector& other) {
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
        std::vector<uint64_t> m_bits;
    };
}
