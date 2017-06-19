#pragma once

#include "../vectors/bit_vector.hpp"
#include "../utils/util.hpp"

namespace tongrams
{
    namespace detail
    {
        template<typename WordGetter>
        struct darray
        {
            darray()
                : m_positions()
            {}

            darray(bit_vector const& bv)
                : m_positions()
            {
                std::vector<uint64_t> const& data = bv.data();
                std::vector<uint64_t> cur_block_positions;
                std::vector<int64_t> block_inventory;
                std::vector<uint16_t> subblock_inventory;
                std::vector<uint64_t> overflow_positions;

                for (size_t word_idx = 0; word_idx < data.size(); ++word_idx)
                {
                    size_t cur_pos = word_idx << 6;
                    uint64_t cur_word = WordGetter()(data, word_idx);
                    unsigned long l;
                    while (util::lsb(cur_word, l)) {
                        cur_pos += l;
                        cur_word >>= l;
                        if (cur_pos >= bv.size()) break;

                        cur_block_positions.push_back(cur_pos);

                        if (cur_block_positions.size() == block_size) {
                            flush_cur_block(cur_block_positions, block_inventory, subblock_inventory, overflow_positions);
                        }

                        // can't do >>= l + 1, can be 64
                        cur_word >>= 1;
                        cur_pos += 1;
                        m_positions += 1;
                    }
                }
                if (cur_block_positions.size()) {
                    flush_cur_block(cur_block_positions, block_inventory, subblock_inventory, overflow_positions);
                }
                m_block_inventory.swap(block_inventory);
                m_subblock_inventory.swap(subblock_inventory);
                m_overflow_positions.swap(overflow_positions);
            }

            void swap(darray& other) {
                std::swap(other.m_positions, m_positions);
                m_block_inventory.swap(other.m_block_inventory);
                m_subblock_inventory.swap(other.m_subblock_inventory);
                m_overflow_positions.swap(other.m_overflow_positions);
            }

            inline uint64_t select(bit_vector const& bv, uint64_t idx) const
            {
                assert(idx < num_positions());
                uint64_t block = idx / block_size;
                int64_t block_pos = m_block_inventory[block];
                if (block_pos < 0) { // sparse super-block
                    uint64_t overflow_pos = uint64_t(-block_pos - 1);
                    return m_overflow_positions[overflow_pos + (idx & (block_size - 1))];
                }

                size_t subblock = idx / subblock_size;
                size_t start_pos = uint64_t(block_pos) + m_subblock_inventory[subblock];
                size_t reminder = idx & (subblock_size - 1);
                std::vector<uint64_t> const& data = bv.data();

                if (!reminder) {
                    return start_pos;
                } else {
                    size_t word_idx = start_pos >> 6;
                    size_t word_shift = start_pos & 63;
                    uint64_t word = WordGetter()(data, word_idx) & (uint64_t(-1) << word_shift);

                    while (true) {
                        size_t popcnt = util::popcount(word);
                        if (reminder < popcnt) break;
                        reminder -= popcnt;
                        word = WordGetter()(data, ++word_idx);
                    }

                    return (word_idx << 6) + util::select_in_word(word, reminder);
                }
            }

            inline uint64_t num_positions() const {
                return m_positions;
            }

            uint64_t bytes() const {
                return sizeof(m_positions)
                     + m_block_inventory.size() * sizeof(m_block_inventory[0])
                     + m_subblock_inventory.size() * sizeof(m_subblock_inventory[0])
                     + m_overflow_positions.size() * sizeof(m_overflow_positions[0]);
            }

            void save(std::ostream& os) const {
                util::save_pod(os, &m_positions);
                util::save_vec(os, m_block_inventory);
                util::save_vec(os, m_subblock_inventory);
                util::save_vec(os, m_overflow_positions);
            }

            void load(std::istream& is) {
                util::load_pod(is, &m_positions);
                util::load_vec(is, m_block_inventory);
                util::load_vec(is, m_subblock_inventory);
                util::load_vec(is, m_overflow_positions);
            }

        protected:
            static void flush_cur_block(std::vector<uint64_t>& cur_block_positions, std::vector<int64_t>& block_inventory,
                                        std::vector<uint16_t>& subblock_inventory, std::vector<uint64_t>& overflow_positions)
            {
                if (cur_block_positions.back() - cur_block_positions.front() < max_in_block_distance) {
                    block_inventory.push_back(int64_t(cur_block_positions.front()));
                    for (size_t i = 0; i < cur_block_positions.size(); i += subblock_size) {
                        subblock_inventory.push_back(uint16_t(cur_block_positions[i] - cur_block_positions.front()));
                    }
                } else {
                    block_inventory.push_back(-int64_t(overflow_positions.size()) - 1);
                    for (size_t i = 0; i < cur_block_positions.size(); ++i) {
                        overflow_positions.push_back(cur_block_positions[i]);
                    }
                    for (size_t i = 0; i < cur_block_positions.size(); i += subblock_size) {
                        subblock_inventory.push_back(uint16_t(-1));
                    }
                }
                cur_block_positions.clear();
            }

            static const size_t block_size = 1024;
            static const size_t subblock_size = 32;
            static const size_t max_in_block_distance = 1 << 16;

            size_t m_positions;
            std::vector<int64_t> m_block_inventory;
            std::vector<uint16_t> m_subblock_inventory;
            std::vector<uint64_t> m_overflow_positions;
        };

        struct identity_getter {
            uint64_t operator()(std::vector<uint64_t> const& data, size_t idx) const {
                return data[idx];
            }
        };

        struct negating_getter {
            uint64_t operator()(std::vector<uint64_t> const& data, size_t idx) const {
                return ~data[idx];
            }
        };
    }

    typedef detail::darray<detail::identity_getter> darray1;
    typedef detail::darray<detail::negating_getter> darray0;
}
