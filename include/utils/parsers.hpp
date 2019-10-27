#pragma once

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "utils/util.hpp"
#include "utils/pools.hpp"
#include "vectors/compact_vector.hpp"

namespace tongrams {

count_record parse_count_line(std::string const& line) {
    auto br = bytes::split_upon_check_end(line, '\t');
    byte_range gram(br.first, br.second);
    count_record record;
    record.gram = gram;

    // if parsed end of line, count will default to 0
    if (br.second ==
        reinterpret_cast<uint8_t const*>(line.c_str()) + line.size()) {
        return record;
    }

    auto pos = reinterpret_cast<const char*>(br.second) + 1;
    uint64_t count = util::toull(pos);
    record.count = count;
    return record;
}

prob_backoff_record parse_prob_backoff_line(std::string const& line) {
    prob_backoff_record record;
    char* end;
    float prob = std::strtod(line.data(), &end);
    if (prob > 0.0) {
        std::cerr << "Warning: positive log10 probability detected."
                  << " This will be mapped to 0." << std::endl;
        prob = 0.0;
    }
    record.prob = prob;

    const uint8_t* buf = reinterpret_cast<uint8_t const*>(line.c_str());
    byte_range br(reinterpret_cast<const uint8_t*>(end) + 1, buf + line.size());

    auto pos = br.first;
    for (; pos != br.second; ++pos) {
        if (*pos == '\t') {
            break;
        }
    }
    record.gram = byte_range(br.first, pos);
    // if pos is equal to the end of the line
    // then backoff is missing and defaults to 0.0
    if (pos != br.second) {
        record.backoff =
            std::strtod(reinterpret_cast<const char*>(pos), nullptr);
    }

    return record;
}

struct prob_backoff_line_handler {
    static const int value_t = value_type::prob_backoff;

    static prob_backoff_record parse_line(std::string const& line) {
        return parse_prob_backoff_line(line);
    }

    static void format_line(prob_backoff_record const& record,
                            std::string& formatted_line) {
        formatted_line =
            std::to_string(record.prob) + "\t" +
            std::string(record.gram.first, record.gram.second)
            // omit 0.0 backoffs
            + (record.backoff ? "\t" + std::to_string(record.backoff) : "");
    }
};

struct count_line_handler {
    static const int value_t = value_type::count;

    static count_record parse_line(std::string const& line) {
        return parse_count_line(line);
    }

    static void format_line(count_record const& record,
                            std::string& formatted_line) {
        formatted_line = std::string(record.gram.first, record.gram.second)
                         // omit 0 counts?
                         + "\t" + std::to_string(record.count);
    }
};

struct arpa_parser {
    friend struct arpa_iterator;

    arpa_parser() : m_cur_line_num(0) {}

    arpa_parser(char const* arpa_filename)
        : m_is(std::ifstream(arpa_filename)), m_cur_line_num(0) {
        if (!m_is.good()) {
            throw std::runtime_error(
                "error in opening arpa file, it may not exist.");
        }
    }

    void read_header(std::vector<uint64_t>& counts) {
        skip_blanklines();
        expect("\\data\\");

        uint32_t i = 1;
        read_line();
        do {
            auto br = bytes::split_upon(m_cur_line, '=');
            uint64_t count =
                util::toull(reinterpret_cast<const char*>(br.second) + 1);

            expect("ngram " + std::to_string(i) + "=" + std::to_string(count));
            counts.push_back(count);
            ++i;
        } while (read_line());

        counts.shrink_to_fit();
        m_max_order = counts.size();

        skip_blanklines();
    }

    void read_values(uint8_t order, uint64_t num_ngrams,
                     std::vector<float>& probs, std::vector<float>& backoffs) {
        assert(order);
        m_order = order;
        skip_blanklines();

        expect("\\" + std::to_string(order) + "-grams:");
        m_offsets.push_back(m_is.tellg());

        for (uint64_t i = 0; i < num_ngrams; ++i) {
            read_line();
            parse_line();
            auto prob = m_cur_record.prob;
            probs.push_back(prob);
            auto backoff = m_cur_record.backoff;
            if (order == 1) {  // unigrams are not quantized
                backoffs.push_back(backoff);
            } else {
                if (backoff) {  // discard 0.0 backoffs for quantization
                    backoffs.push_back(backoff);
                }
            }
        }

        read_line();
    }

    bool read_line() {
        std::getline(m_is, m_cur_line);
        ++m_cur_line_num;
        if (m_cur_line_num % 100000000 == 0) {
            util::logger("Processed " + std::to_string(m_cur_line_num) +
                         " lines");
        }
        return !m_cur_line.empty();
    }

    void parse_line() {
        m_cur_record = parse_prob_backoff_line(m_cur_line);
    }

    template <typename Pool>
    bool parse_line(Pool& pool) {
        m_cur_record = parse_prob_backoff_line(m_cur_line);
        return pool.append(m_cur_record);
    }

    void skip_blanklines() {
        while (m_cur_line.empty()) {
            read_line();
        }
    }

    void parse_eof() {
        read_line();
        expect("\\end\\");
    }

    std::vector<uint64_t> offsets() const {
        return m_offsets;
    }

    prob_backoff_record parsed_line() const {
        return m_cur_record;
    }

    void expect(std::string const& s) {
        if (m_cur_line != s) {
            std::cerr << "Error during parsing arpa file at "
                      << "line " << m_cur_line_num << ": "
                      << "expected '" << s << "' but got '" << m_cur_line << "'"
                      << std::endl;
            std::abort();
        }
    }

private:
    std::ifstream m_is;
    std::string m_cur_line;
    prob_backoff_record m_cur_record;
    uint64_t m_cur_line_num;
    std::vector<uint64_t> m_offsets;
    uint8_t m_order;
    uint8_t m_max_order;

    void skip_to(uint64_t offset) {
        m_is.seekg(offset, m_is.beg);
    }
};

template <typename Buffer>
struct count_records_iterator {
    count_records_iterator(Buffer* buf, size_t line_num)
        : m_cur_line_num(line_num), m_buf(buf) {
        if (m_buf) {
            std::getline(*m_buf, m_cur_line);
        }
    }

    count_record operator*() {
        return parse_count_line(m_cur_line);
    }

    count_records_iterator& operator++() {
        std::getline(*m_buf, m_cur_line);
        ++m_cur_line_num;
        if (m_cur_line_num % 100000000 == 0) {
            util::logger("Processed " + std::to_string(m_cur_line_num) +
                         " lines");
        }
        return *this;
    }

    bool operator==(count_records_iterator const& other) const {
        return m_cur_line_num == other.m_cur_line_num;
    }

    bool operator!=(count_records_iterator const& other) const {
        return !(*this == other);
    }

private:
    size_t m_cur_line_num;
    std::string m_cur_line;
    Buffer* m_buf;
};

struct grams_gzparser {
    typedef count_records_iterator<boost::iostreams::filtering_istream>
        iterator;

    grams_gzparser(char const* grams_filename)
        : m_is(std::ifstream(grams_filename,
                             std::ios_base::in | std::ios_base::binary))
        , m_num_lines(0) {
        if (!m_is.good()) {
            throw std::runtime_error(
                "error in opening grams file, it may not exist.");
        }

        m_fi.push(boost::iostreams::gzip_decompressor());
        m_fi.push(m_is);

        std::string line;
        std::getline(m_fi, line);
        if (line.empty()) {
            throw std::runtime_error(
                "first line must be non-empty and contain the number of "
                "lines.");
        }

        try {
            m_num_lines = std::stoull(line);
            if (!m_num_lines) {
                throw std::runtime_error("number of lines must not be 0.");
            }
        } catch (std::invalid_argument& e) {
            std::cerr << e.what() << std::endl;
            std::cerr << "firt line must contain the number of lines."
                      << std::endl;
            exit(1);
        }
    }

    uint64_t num_lines() {
        return m_num_lines;
    }

    iterator begin() {
        return iterator(&m_fi, 0);
    }

    iterator end() {
        return iterator(nullptr, m_num_lines);
    }

private:
    std::ifstream m_is;
    boost::iostreams::filtering_istream m_fi;
    uint64_t m_num_lines;
};

// NOTE: non-gzipped version
// struct grams_parser
// {
//     typedef count_records_iterator<std::ifstream> iterator;

//     grams_parser(char const* grams_filename)
//         : m_is(std::ifstream(grams_filename))
//         , m_num_lines(0)
//     {
//         if (!m_is.good()) {
//             throw std::runtime_error("error in opening grams file, it may not
//             exist.");
//         }

//         std::string line;
//         std::getline(m_is, line);
//         if (line.empty()) {
//             throw std::runtime_error("first line must be non-empty and
//             contain the number of lines.");
//         }
//         m_num_lines = std::stoull(line);
//         if (!m_num_lines) {
//             throw std::runtime_error("number of lines must not be 0.");
//         }
//     }

//     uint64_t num_lines() {
//         return m_num_lines;
//     }

//     iterator begin() {
//         return iterator(&m_is, 0);
//     }

//     iterator end() {
//         return iterator(nullptr, m_num_lines);
//     }

// private:
//     std::ifstream m_is;
//     uint64_t m_num_lines;
// };
}  // namespace tongrams
