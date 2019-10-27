#pragma once

#include <boost/iostreams/device/mapped_file.hpp>
#include <sys/mman.h>

#include "utils/parsers.hpp"
#include "utils/util_types.hpp"

namespace tongrams {

struct text_lines {
    text_lines(const char* filename) : m_pos(0), m_num_words(0), m_eol(false) {
        m_file.open(filename);
        if (!m_file.is_open()) {
            throw std::runtime_error("Error opening file");
        }

        m_data = (uint8_t const*)m_file.data();
        m_size = m_file.size() / sizeof(m_data[0]);

        auto ret = posix_madvise((void*)m_data, m_size, POSIX_MADV_SEQUENTIAL);
        if (ret) {
            std::cerr << "Error calling madvice: " << errno << std::endl;
        }
    }

    void next_word(byte_range& word) {
        uint64_t pos = m_pos;
        for (; m_data[pos] != ' '; ++pos) {
            if (m_data[pos] == '\n') {
                m_eol = true;
                break;
            }
        }
        ++m_num_words;
        word.first = &m_data[m_pos];
        word.second = &m_data[pos];
        m_pos = pos + 1;
    }

    void begin_line() {
        m_eol = false;
    }

    bool end_of_line() const {
        return m_eol;
    }

    bool end_of_file() const {
        return m_pos == m_size;
    }

    uint64_t num_words() const {
        return m_num_words;
    }

private:
    boost::iostreams::mapped_file_source m_file;
    uint8_t const* m_data;
    size_t m_size;
    uint64_t m_pos;
    uint64_t m_num_words;
    bool m_eol;
};

struct arpa_iterator {
    arpa_iterator(char const* arpa_filename, uint8_t order, uint64_t offset)
        : m_num_grams(0), m_arpa_parser(arpa_filename) {
        std::vector<uint64_t> counts;
        m_arpa_parser.read_header(counts);
        m_arpa_parser.skip_to(offset);
        m_arpa_parser.m_order = order;
        m_num_grams = counts[order - 1];
        if (order > m_arpa_parser.m_max_order) {
            throw std::invalid_argument("order must be less than max_order");
        }
    }

    uint64_t num_grams() const {
        return m_num_grams;
    }

    prob_backoff_record next() {
        m_arpa_parser.read_line();
        return parse_prob_backoff_line(m_arpa_parser.m_cur_line);
    }

private:
    uint64_t m_num_grams;
    arpa_parser m_arpa_parser;
};

struct forward_byte_range_iterator {
    void init(byte_range const& range) {
        m_cur_pos = range.first;
        m_begin = range.first;
        m_end = range.second;
    }

    uint64_t spaces() const {
        return std::count(m_begin, m_end, ' ');
    }

    byte_range next() {
        auto pos = m_cur_pos;
        for (; pos != m_end; ++pos) {
            if (*pos == ' ' || *pos == '\n') {
                break;
            }
        }
        auto br = byte_range(m_cur_pos, pos);
        m_cur_pos = pos + 1;
        return br;
    }

private:
    uint8_t const* m_cur_pos;
    uint8_t const* m_begin;
    uint8_t const* m_end;
};

struct backward_byte_range_iterator {
    void init(byte_range const& range) {
        m_begin = range.first;
        m_end = range.second;
        m_cur_pos = m_end;
    }

    uint64_t spaces() const {
        return std::count(m_begin, m_end, ' ');
    }

    byte_range next() {
        auto pos = m_cur_pos;
        --pos;
        for (; pos != m_begin; --pos) {
            if (*pos == ' ') {
                break;
            }
        }
        auto br = byte_range(pos != m_begin ? pos + 1 : m_begin, m_cur_pos);
        m_cur_pos = pos;
        return br;
    }

private:
    uint8_t const* m_cur_pos;
    uint8_t const* m_begin;
    uint8_t const* m_end;
};

}  // namespace tongrams
