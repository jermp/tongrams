#pragma once

#include "util_types.hpp"

#include <tuple>

namespace tongrams
{
    struct strings_pool
    {
        strings_pool()
        {}

        void reserve(size_t num_bytes = 0) {
            m_data.reserve(num_bytes);
        }

        void append(byte_range br) {
            m_data.insert(m_data.end(),
                          br.first, br.second);
        }

        uint8_t const* base_addr() const {
            return reinterpret_cast<uint8_t const*>(m_data.data());
        }

        size_t bytes() const {
            return m_data.size();
        }

        inline byte_range get_bytes(uint8_t const* base_addr,
                                    size_t begin, size_t end) const {
            return {base_addr + begin, base_addr + end};
        }

    private:
        std::vector<uint8_t> m_data;
    };

    struct grams_probs_pool
    {
        grams_probs_pool(size_t num_bytes = 0)
            : m_max_bytes(num_bytes)
        {
            m_strings_pool.reserve(num_bytes);
            m_base_addr =
                reinterpret_cast<uint8_t const*>(m_strings_pool.data());
        }

        grams_probs_pool(size_t num_index_entries, size_t num_bytes)
            : grams_probs_pool(num_bytes)
        {
            m_index.reserve(num_index_entries);
        }

        bool append(prob_backoff_record const& record)
        {
            auto begin = m_strings_pool.size();
            auto gram = record.gram;
            size_t gram_bytes = gram.second - gram.first;

            if (gram_bytes) {
                if (m_strings_pool.size() + gram_bytes > m_max_bytes) {
                    return false;
                }

                m_strings_pool.insert(m_strings_pool.end(),
                                      gram.first, gram.second);
                auto end = m_strings_pool.size();

                m_index.emplace_back(byte_range(m_base_addr + begin,
                                                m_base_addr + end),
                                                record.prob,
                                                record.backoff);

                if (m_strings_pool.size() % GB == 0) {
                    util::logger("Loaded " + std::to_string(m_strings_pool.size()) + " bytes");
                }
            }

            return true;
        }

        void clear() {
            m_index.clear();
            m_strings_pool.clear();
        }

        std::vector<prob_backoff_record>& index() {
            return m_index;
        }

    private:
        size_t m_max_bytes;
        uint8_t const* m_base_addr;
        std::vector<prob_backoff_record> m_index;
        std::vector<uint8_t> m_strings_pool;
    };

    struct grams_counts_pool
    {
        grams_counts_pool(size_t num_bytes = 0)
            : m_max_bytes(num_bytes)
        {
            m_strings_pool.reserve(num_bytes);
            m_base_addr =
                reinterpret_cast<uint8_t const*>(m_strings_pool.data());
        }

        grams_counts_pool(size_t num_index_entries, size_t num_bytes)
            : grams_counts_pool(num_bytes)
        {
            m_index.reserve(num_index_entries);
        }

        bool append(count_record const& record)
        {
            auto begin = m_strings_pool.size();
            auto gram = record.gram;
            size_t gram_bytes = gram.second - gram.first;

            if (gram_bytes) {
                if (m_strings_pool.size() + gram_bytes > m_max_bytes) {
                    return false;
                }

                m_strings_pool.insert(m_strings_pool.end(),
                                      gram.first, gram.second);
                auto end = m_strings_pool.size();

                m_index.emplace_back(byte_range(m_base_addr + begin,
                                                m_base_addr + end),
                                                record.count);

                if (m_strings_pool.size() % GB == 0) {
                    util::logger("Loaded " + std::to_string(m_strings_pool.size()) + " bytes");
                }
            }

            return true;
        }

        template<typename Parser>
        void load_from(const char* filename) {
            clear();
            Parser parser(filename);
            uint64_t index_entries = parser.num_lines();
            m_index.reserve(index_entries);
            for (auto const& l: parser) {
                if (!append(count_record(l.gram, l.count))) {
                    throw std::runtime_error("max available memory pool excedeed");
                }
            }
        }

        void clear() {
            m_index.clear();
            m_strings_pool.clear();
        }

        std::vector<count_record>& index() {
            return m_index;
        }

    private:
        size_t m_max_bytes;
        uint8_t const* m_base_addr;
        std::vector<count_record> m_index;
        std::vector<uint8_t> m_strings_pool;
    };
}
