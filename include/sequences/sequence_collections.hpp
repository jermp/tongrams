#pragma once

#include <algorithm>
#include <functional>
#include <unordered_map>
#include <unordered_set>

#include "utils/util.hpp"

namespace tongrams {

struct quantized_sequence_collection {
    struct builder {
        builder(size_t n = 0) {
            m_sequences.reserve(n);
        }

        void resize(size_t n) {
            m_quantization_bits.resize(n);
            m_sequences.resize(n);
        }

        void add_sequence(uint8_t order, uint8_t quantization_bits,
                          std::vector<float>& averages) {
            check_quantization_bits(quantization_bits, averages.size());
            m_quantization_bits[order - 1] = quantization_bits;
            m_sequences[order - 1].swap(averages);
        }

        uint64_t rank(uint8_t order_m1, float value, size_t reserved) const {
            // NOTE: in order to return the same results as KenLM
            if (reserved == 1  // value is a backoff
                && value == 0.0) {
                return 0;
            }

            assert(order_m1 < m_sequences.size());
            auto const& data = m_sequences[order_m1];

            auto above =
                std::lower_bound(data.begin() + reserved, data.end(), value);
            if (above == data.begin() + reserved) {
                return reserved;
            }
            if (above == data.end()) {
                return data.end() - data.begin() - 1;
            }
            auto prev_above = above;
            --prev_above;
            auto to_ret =
                above - data.begin() - (value - *prev_above < *above - value);
            assert(size_t(to_ret) < data.size());

            return to_ret;
        }

        void make_bins(std::vector<float>& sequence,
                       std::vector<float>& averages, size_t num_bins) {
            std::vector<float>::const_iterator start = sequence.begin(), finish;
            for (uint32_t i = 0; i < num_bins; ++i, start = finish) {
                finish = sequence.begin() +
                         ((sequence.size() * static_cast<uint64_t>(i + 1)) /
                          num_bins);
                if (finish == start) {
                    averages.push_back(
                        i ? averages.back()
                          : -std::numeric_limits<float>::infinity());
                } else {
                    averages.push_back(std::accumulate(start, finish, 0.0) /
                                       static_cast<float>(finish - start));
                }
            }
        }

        void check_quantization_bits(uint8_t quantization_bits, size_t n) {
            if (quantization_bits < 2 || quantization_bits > 32) {
                throw std::invalid_argument(
                    "quantization bits must be in [2, 32]");
            }
            uint64_t num_bins = uint64_t(1) << quantization_bits;
            if (n && n < num_bins) {
                std::cerr
                    << "too many quantization bits used: these must be in [2, "
                    << util::ceil_log2(n) << "]" << std::endl;
                std::abort();
            }
        }

        void build_backoffs_sequence(std::vector<float>& sequence,
                                     uint8_t quantization_bits) {
            check_quantization_bits(quantization_bits, sequence.size());
            m_quantization_bits.push_back(quantization_bits);
            std::sort(sequence.begin(), sequence.end());
            std::vector<float> averages;
            uint64_t num_bins = uint64_t(1) << quantization_bits;
            averages.reserve(num_bins);
            averages.push_back(0.0);
            make_bins(sequence, averages, num_bins - 1);
            assert(averages.size() == num_bins);
            m_sequences.push_back(std::move(averages));
        }

        void build_probs_sequence(std::vector<float>& sequence,
                                  uint8_t quantization_bits) {
            check_quantization_bits(quantization_bits, sequence.size());
            m_quantization_bits.push_back(quantization_bits);
            std::sort(sequence.begin(), sequence.end());
            std::vector<float> averages;
            uint64_t num_bins = uint64_t(1) << quantization_bits;
            averages.reserve(num_bins);
            make_bins(sequence, averages, num_bins);
            assert(averages.size() == num_bins);
            m_sequences.push_back(std::move(averages));
        }

        uint64_t quantization_bits(uint64_t order_m1) const {
            assert(order_m1 < m_quantization_bits.size());
            return m_quantization_bits[order_m1];
        }

        void swap(builder& other) {
            m_quantization_bits.swap(other.m_quantization_bits);
            m_sequences.swap(other.m_sequences);
        }

        void build(quantized_sequence_collection& qsc) {
            qsc.m_quantization_bits.swap(m_quantization_bits);
            qsc.m_sequences.swap(m_sequences);
            builder().swap(*this);
        }

    private:
        // NOTE: use a vector since we could specify
        // different quantization bits for different orders
        // (not yet supported)
        std::vector<uint8_t> m_quantization_bits;
        std::vector<std::vector<float>> m_sequences;
    };

    quantized_sequence_collection() {}

    inline float access(uint64_t order_m1, uint64_t i) const {
        assert(order_m1 < m_sequences.size());
        assert(i < m_sequences[order_m1].size());
        return m_sequences[order_m1][i];
    }

    size_t bytes() const {
        size_t bytes = m_quantization_bits.size() * sizeof(uint8_t);
        for (auto const& s : m_sequences) {
            bytes += s.size() * sizeof(float);
        }
        return bytes;
    }

    size_t size(uint64_t order_m1) const {
        return m_sequences[order_m1].size();
    }

    uint8_t quantization_bits(uint64_t order_m1) const {
        return m_quantization_bits[order_m1];
    }

    void swap(quantized_sequence_collection& other) {
        m_quantization_bits.swap(other.m_quantization_bits);
        m_sequences.swap(other.m_sequences);
    }

    void save(std::ostream& os) const {
        os.write(
            reinterpret_cast<char const*>(m_quantization_bits.data()),
            (std::streamsize)(sizeof(uint8_t) * m_quantization_bits.size()));
        for (auto const& s : m_sequences) {
            essentials::save_vec(os, s);
        }
    }

    void load(std::istream& is, uint64_t order) {
        m_quantization_bits.resize(order);
        is.read(reinterpret_cast<char*>(m_quantization_bits.data()),
                (std::streamsize)(sizeof(uint8_t) * order));
        m_sequences.resize(order);
        for (auto& s : m_sequences) {
            essentials::load_vec(is, s);
        }
    }

private:
    std::vector<uint8_t> m_quantization_bits;
    std::vector<std::vector<float>> m_sequences;
};

template <typename CountType>
struct byte_aligned_sequence_collection {
    struct builder {
        struct adaptor {
            uint64_t first(pairs_vector const& s, uint64_t i) const {
                return s[i].first;
            }

            uint64_t second(pairs_vector const& s, uint64_t i) const {
                return s[i].second;
            }
        };

        builder(size_t n = 0) {
            m_pairs_sequences.reserve(n);
            m_sequences.reserve(n);
        }

        template <typename Iterator>
        void build_sequence(Iterator begin, uint64_t n) {
            if (n) {
                std::unordered_map<uint64_t, uint64_t> x;
                auto end = x.end();
                uint64_t max = 0;
                for (uint64_t i = 0; i < n; ++i, ++begin) {
                    uint64_t v = *begin;
                    if (v > max) {
                        max = v;
                    }
                    if (x.find(v) != end) {
                        ++x[v];
                    } else {
                        x.emplace(v, 1);
                    }
                }

                uint64_t num_distinct_values = x.size();
                pairs_vector sorted_x;
                sorted_x.reserve(num_distinct_values);
                for (auto v : x) {
                    sorted_x.emplace_back(v.first, v.second);
                }

                // sort on frequency of counters
                std::sort(sorted_x.begin(), sorted_x.end(),
                          [&](uint64_pair const& x, uint64_pair const& y) {
                              return x.second > y.second;
                          });

                // assign ranks
                std::vector<CountType> distinct_values;
                distinct_values.reserve(num_distinct_values);
                for (uint64_t i = 0; i < num_distinct_values; ++i) {
                    auto& p = sorted_x[i];
                    distinct_values.push_back(p.first);
                    p.second = i;
                }
                m_sequences.push_back(std::move(distinct_values));

                // sort on values to enable binary search
                std::sort(sorted_x.begin(), sorted_x.end(),
                          [&](uint64_pair const& x, uint64_pair const& y) {
                              return x.first < y.first;
                          });
                m_pairs_sequences.push_back(std::move(sorted_x));
            } else {  // push empty sequences
                m_sequences.push_back(std::vector<CountType>());
                m_pairs_sequences.push_back(pairs_vector());
            }
        }

        uint64_t rank(uint8_t order_m1, uint64_t value) const {
            assert(order_m1 < m_pairs_sequences.size());
            auto const& distinct_values = m_pairs_sequences[order_m1];
            uint64_t rank = 0;
            if (!util::binary_search(distinct_values, distinct_values.size(),
                                     value, rank, adaptor())) {
                throw std::runtime_error("value not found");
            }
            return rank;
        }

        size_t size(uint64_t order_m1) const {
            return m_sequences[order_m1].size();
        }

        void swap(builder& other) {
            m_sequences.swap(other.m_sequences);
        }

        void build(byte_aligned_sequence_collection& sc) {
            sc.m_sequences.swap(m_sequences);
            builder().swap(*this);
        }

    private:
        std::vector<pairs_vector> m_pairs_sequences;
        std::vector<std::vector<CountType>> m_sequences;
    };

    byte_aligned_sequence_collection() {}

    inline uint64_t access(uint64_t order_m1, uint64_t i) const {
        assert(order_m1 < m_sequences.size());
        return m_sequences[order_m1][i];
    }

    size_t bytes() const {
        size_t bytes = 0;
        for (auto const& s : m_sequences) {
            bytes += s.size() * sizeof(CountType);
        }
        return bytes;
    }

    size_t size(uint64_t order_m1) const {
        return m_sequences[order_m1].size();
    }

    void swap(byte_aligned_sequence_collection& other) {
        m_sequences.swap(other.m_sequences);
    }

    void save(std::ostream& os) const {
        for (auto const& s : m_sequences) {
            essentials::save_vec(os, s);
        }
    }

    void load(std::istream& is, uint8_t order) {
        m_sequences.resize(order);
        for (auto& s : m_sequences) {
            essentials::load_vec(is, s);
        }
    }

private:
    std::vector<std::vector<CountType>> m_sequences;
};

struct sequence_collection {
    struct builder {
        struct adaptor {
            uint64_t first(pairs_vector const& s, uint64_t i) const {
                return s[i].first;
            }

            uint64_t second(pairs_vector const& s, uint64_t i) const {
                return s[i].second;
            }
        };

        builder(size_t n = 0) {
            m_pairs_sequences.reserve(n);
            m_sequences.reserve(n);
        }

        template <typename Iterator>
        void build_sequence(Iterator begin, uint64_t n) {
            if (n) {
                std::unordered_map<uint64_t, uint64_t> x;
                auto end = x.end();
                uint64_t max = 0;
                for (uint64_t i = 0; i < n; ++i, ++begin) {
                    uint64_t v = *begin;
                    if (v > max) {
                        max = v;
                    }
                    if (x.find(v) != end) {
                        ++x[v];
                    } else {
                        x.emplace(v, 1);
                    }
                }

                uint64_t m = x.size();
                pairs_vector sorted_x;
                sorted_x.reserve(m);
                for (auto v : x) {
                    sorted_x.emplace_back(v.first, v.second);
                }

                // sort on frequency of counters
                std::sort(sorted_x.begin(), sorted_x.end(),
                          [&](uint64_pair const& x, uint64_pair const& y) {
                              return x.second > y.second;
                          });

                // assign ranks
                compact_vector::builder cvb(m, util::ceil_log2(max + 1));
                for (uint64_t i = 0; i < m; ++i) {
                    auto& p = sorted_x[i];
                    cvb.push_back(p.first);
                    p.second = i;
                }
                m_sequences.emplace_back(cvb);

                // sort on values to enable binary search
                std::sort(sorted_x.begin(), sorted_x.end(),
                          [&](uint64_pair const& x, uint64_pair const& y) {
                              return x.first < y.first;
                          });
                m_pairs_sequences.push_back(std::move(sorted_x));
            } else {  // push empty sequences
                m_sequences.push_back(compact_vector());
                m_pairs_sequences.push_back(pairs_vector());
            }
        }

        uint64_t rank(uint8_t order_m1, uint64_t value) const {
            assert(order_m1 < m_pairs_sequences.size());
            auto const& distinct_values = m_pairs_sequences[order_m1];
            uint64_t rank = 0;
            if (!util::binary_search(distinct_values, distinct_values.size(),
                                     value, rank, adaptor())) {
                throw std::runtime_error("value not found");
            }
            return rank;
        }

        size_t size(uint64_t order_m1) const {
            return m_sequences[order_m1].size();
        }

        void swap(builder& other) {
            m_sequences.swap(other.m_sequences);
        }

        void build(sequence_collection& sc) {
            sc.m_sequences.swap(m_sequences);
            builder().swap(*this);
        }

    private:
        std::vector<pairs_vector> m_pairs_sequences;
        std::vector<compact_vector> m_sequences;
    };

    sequence_collection() {}

    inline uint64_t access(uint64_t order_m1, uint64_t i) const {
        assert(order_m1 < m_sequences.size());
        return m_sequences[order_m1][i];
    }

    size_t bytes() const {
        size_t bytes = 0;
        for (auto const& s : m_sequences) {
            bytes += s.bytes();
        }
        return bytes;
    }

    size_t size(uint64_t order_m1) const {
        return m_sequences[order_m1].size();
    }

    void swap(sequence_collection& other) {
        m_sequences.swap(other.m_sequences);
    }

    void save(std::ostream& os) const {
        for (auto const& s : m_sequences) {
            s.save(os);
        }
    }

    void load(std::istream& is, uint8_t order) {
        m_sequences.resize(order);
        for (auto& s : m_sequences) {
            s.load(is);
        }
    }

private:
    std::vector<compact_vector> m_sequences;
};

}  // namespace tongrams
