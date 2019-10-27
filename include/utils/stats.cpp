#pragma once

#include <unordered_map>

#include "mph_count_lm.hpp"
#include "mph_prob_lm.hpp"
#include "trie_count_lm.hpp"
#include "trie_prob_lm.hpp"

namespace tongrams {

template <typename Values, typename KeyRankSequence, typename BaseHasher>
void mph_count_lm<Values, KeyRankSequence, BaseHasher>::print_stats(
    size_t bytes) const {
    util::logger("========= MPH_COUNT_LM statistics =========");
    uint64_t num_grams = size();
    std::cout << "order: " << order() << "\n";
    std::cout << "num. of grams: " << num_grams << "\n";
    std::cout << "tot. bytes: " << bytes << "\n";
    std::cout << "bytes per gram: " << double(bytes) / num_grams << "\n";
    std::cout << "unique values bytes: " << m_distinct_counts.bytes() << "\n";
    uint64_t i = 1;
    size_t data_bytes = m_distinct_counts.bytes();
    for (auto const& t : m_tables) {
        std::cout << i << "-grams stats:\n";
        std::cout << "\tdistinct frequency counters: "
                  << m_distinct_counts.size(i - 1) << "\n";
        size_t x = t.data_bytes();
        std::cout << "\tdata table bytes: " << x << std::endl;
        std::cout << "\t(does NOT include hash function bytes)\n";
        data_bytes += x;
        ++i;
    }

    uint64_t hash_key_bytes = KeyRankSequence::hash_bits / 8;
    std::cout << "hash keys bytes: " << hash_key_bytes * num_grams << " ("
              << hash_key_bytes * num_grams * 100.0 / bytes << "%)\n"
              << "\tper gram: " << double(hash_key_bytes) / num_grams
              << std::endl;

    uint64_t counts_bytes = data_bytes - hash_key_bytes * num_grams;
    std::cout << "counts bytes: " << counts_bytes << " ("
              << counts_bytes * 100.0 / bytes << "%)\n"
              << "\tper gram: " << double(counts_bytes) / num_grams
              << std::endl;

    uint64_t hash_function_bytes = bytes - data_bytes;
    std::cout << "hash functions bytes: " << hash_function_bytes << " ("
              << hash_function_bytes * 100.0 / bytes << "%)\n"
              << "\tper gram: " << double(hash_function_bytes) / num_grams
              << std::endl;

    std::cout << "hash functions bits per gram: "
              << hash_function_bytes * 8.0 / num_grams << std::endl;
}

template <typename Values, typename KeyRankSequence, typename BaseHasher>
void mph_prob_lm<Values, KeyRankSequence, BaseHasher>::print_stats(
    size_t bytes) const {
    util::logger("========= MPH_PROB_LM statistics =========");
    uint64_t num_grams = size();
    std::cout << "order: " << order() << "\n";
    std::cout << "num. of grams: " << num_grams << "\n";
    std::cout << "tot. bytes: " << bytes << "\n";
    std::cout << "bytes per gram: " << double(bytes) / num_grams << "\n";
    std::cout << "quantized probs bytes: " << m_probs_averages.bytes() << "\n";
    std::cout << "quantized backoffs bytes: " << m_backoffs_averages.bytes()
              << "\n";
    uint64_t i = 1;
    size_t data_bytes = m_probs_averages.bytes() + m_backoffs_averages.bytes();
    for (auto const& t : m_tables) {
        std::cout << i << "-grams stats:\n";
        size_t x = t.data_bytes();
        std::cout << "\tdata table bytes: " << x << std::endl;
        std::cout << "\t(does NOT include hash function bytes)\n";
        data_bytes += x;
        ++i;
    }

    uint64_t hash_key_bytes = KeyRankSequence::hash_bits / 8;
    std::cout << "hash keys bytes: " << hash_key_bytes * num_grams << " ("
              << hash_key_bytes * num_grams * 100.0 / bytes << "%)\n"
              << "\tper gram: " << double(hash_key_bytes) / num_grams
              << std::endl;

    uint64_t ranks_bytes = data_bytes - hash_key_bytes * num_grams;
    std::cout << "ranks bytes: " << ranks_bytes << " ("
              << ranks_bytes * 100.0 / bytes << "%)\n"
              << "\tper gram: " << double(ranks_bytes) / num_grams << std::endl;

    uint64_t hash_function_bytes = bytes - data_bytes;
    std::cout << "hash functions bytes: " << hash_function_bytes << " ("
              << hash_function_bytes * 100.0 / bytes << "%)\n"
              << "\tper gram: " << double(hash_function_bytes) / num_grams
              << std::endl;

    std::cout << "hash functions bits per gram: "
              << hash_function_bytes * 8.0 / num_grams << std::endl;
}

struct topk_queue {
    topk_queue(size_t size) {
        m_data.reserve(size);
    }

    void add(uint64_t x) {
        if (m_data.size() < k()) {
            m_data.push_back(x);
            std::push_heap(m_data.begin(), m_data.end(),
                           std::greater<uint64_t>());
            return;
        }

        if (x > m_data.front()) {
            std::pop_heap(m_data.begin(), m_data.end(),
                          std::greater<uint64_t>());
            m_data.back() = x;
            std::push_heap(m_data.begin(), m_data.end(),
                           std::greater<uint64_t>());
        }
    }

    size_t k() const {
        return m_data.capacity();
    }

    void finalize() {
        std::sort(m_data.begin(), m_data.end(), std::greater<uint64_t>());
    }

    void print() const {
        for (auto v : m_data) {
            std::cout << v << " ";
        }
        std::cout << std::endl;
    }

private:
    std::vector<uint64_t> m_data;
};

template <typename Grams, typename Ranks, typename Pointers>
void sorted_array<Grams, Ranks, Pointers>::print_stats(
    uint8_t order, pointer_sequence<Pointers> const* ranges) const {
    uint64_t n = size();
    std::cout << "\tnum grams: " << n << "\n";
    if (m_grams.size()) {
        uint64_t u = m_grams.universe();
        std::cout << "\tuniverse: " << u << "\n";
        double avg_gap = double(u) / n;
        std::cout << "\tavg gap: " << avg_gap << "\n";
        // std::cout << "\tEF_single bytes: "
        //           << uint64_t(double(n * (util::ceil_log2(avg_gap) + 2)) / 8)
        //           << "\n";
    }

    auto m = m_pointers.size();
    if (m) {
        uint64_t sum = 0;
        uint64_t max = 0;
        uint64_t min = uint64_t(-1);
        uint64_t non_null_ranges = 0;

        {
            topk_queue top_k(5);
            auto it = m_pointers.begin();
            uint64_t begin = it.next();
            std::cout << "\tnum ptrs.: " << m << "\n";
            for (uint64_t i = 0; i < m - 1; ++i) {
                uint64_t end = it.next();
                uint64_t range = end - begin;
                if (range) {
                    if (range < min) {
                        min = range;
                    }
                    if (range > max) {
                        max = range;
                    }
                    sum += range;
                    ++non_null_ranges;
                    top_k.add(range);
                }
                begin = end;
            }
            std::cout << "\tnon-empty ranges: " << non_null_ranges << "\n";
            std::cout << "\tdensity: " << non_null_ranges * 100.0 / m << "%\n";
            std::cout << "\tmin range: " << min << "\n";
            std::cout << "\tmax range: " << max << "\n";
            std::cout << "\tavg range: " << double(sum) / non_null_ranges
                      << std::endl;
            std::cout << "\ttop-" << top_k.k() << " ranges: ";
            top_k.finalize();
            top_k.print();
        }

        // {
        //     auto it = m_pointers.begin();
        //     uint64_t begin = it.next();
        //     double percs[5] = {0.01, 0.1, 1, 10, 20};
        //     uint64_t leq_range_sizes[5] = {0, 0, 0, 0, 0};
        //     uint64_t geq_range_sizes[5] = {0, 0, 0, 0, 0};
        //     uint64_t leq_ints[5] = {0, 0, 0, 0, 0};
        //     uint64_t geq_ints[5] = {0, 0, 0, 0, 0};
        //     for (uint64_t i = 0; i < m - 1; ++i) {
        //         uint64_t end = it.next();
        //         uint64_t range = end - begin;
        //         if (range) {
        //             for (uint32_t i = 0; i < 5; ++i) {
        //                 if (range <= max * percs[i] / 100) {
        //                     ++leq_range_sizes[i];
        //                     leq_ints[i] += range;
        //                 } else {
        //                     ++geq_range_sizes[i];
        //                     geq_ints[i] += range;
        //                 }
        //             }
        //         }
        //         begin = end;
        //     }
        //     for (uint32_t i = 0; i < 5; ++i) {
        //         std::cout <<
        //         "\t========================================================================="
        //         << std::endl; std::cout << "\tnum of ranges whose length is
        //         <= " << percs[i]
        //                   << "% of max range (" << max * percs[i] / 100 << ")
        //                   = "
        //                   << leq_range_sizes[i] << " (" << leq_range_sizes[i]
        //                   * 100.0 / non_null_ranges << "%)" << std::endl;
        //         std::cout << "\tnumber of ints belonging to such ranges = "
        //         << leq_ints[i]
        //                   << " (" << leq_ints[i] * 100.0 /
        //                   next_order_num_grams << "%)" << std::endl;
        //         std::cout << "\tnum of ranges whose length is > " << percs[i]
        //                   << "% of max range (" << max * percs[i] / 100 << ")
        //                   = "
        //                   << geq_range_sizes[i] << " (" << geq_range_sizes[i]
        //                   * 100.0 / non_null_ranges << "%)" << std::endl;
        //         std::cout << "\tnumber of ints belonging to such ranges = "
        //         << geq_ints[i]
        //                   << " (" << geq_ints[i] * 100.0 /
        //                   next_order_num_grams << "%)" << std::endl;
        //     }
        // }
    }

    if (order == 2 || order == 3) {
        std::unordered_map<uint32_t, uint32_t> occs;

        auto pointers_it = ranges->begin();
        uint64_t m = ranges->size();
        uint64_t begin = pointers_it.next();
        uint64_t prev_upper = 0;
        auto grams_it = m_grams.begin();
        auto occs_end = occs.end();

        double H_1 = 0.0;

        for (uint64_t i = 0; i < m - 1; ++i) {
            uint64_t end = pointers_it.next();
            uint64_t range = end - begin;
            if (range) {
                H_1 += range * std::log2(range);
            }

            for (uint64_t j = 0; j < range; ++j) {
                uint64_t v = grams_it.next();
                uint64_t original_v = v - prev_upper;
                if (occs.find(original_v) != occs_end) {
                    ++occs[original_v];
                } else {
                    occs.emplace(original_v, 1);
                }

                if (j == range - 1) {
                    prev_upper = v;
                }
            }

            begin = end;
        }

        H_1 /= n;
        double H_0 = 0.0;
        for (auto const& pair : occs) {
            uint32_t occ = pair.second;
            H_0 += occ * std::log2(double(n) / occ);
        }
        H_0 /= n;

        if (order == 2) {
            std::cout << "\tH_0 = " << H_0 / 8 << " bytes" << std::endl;
            std::cout << "\tH_1 = " << H_1 / 8 << " bytes" << std::endl;
        }

        if (order == 3) {
            std::cout << "\tH_2 = " << H_1 / 8 << " bytes" << std::endl;
        }
    }
}

template <typename Vocabulary, typename Mapper, typename Values, typename Ranks,
          typename Grams, typename Pointers>
void trie_count_lm<Vocabulary, Mapper, Values, Ranks, Grams,
                   Pointers>::print_stats(size_t bytes) const {
    util::logger("========= TRIE_COUNT_LM statistics =========");
    uint64_t i = 1;
    uint64_t grams_bytes = 0;
    uint64_t counts_ranks_bytes = 0;
    uint64_t pointers_bytes = 0;

    for (auto const& a : m_arrays) {
        std::cout << i << "-grams bytes:\n";
        uint64_t x = a.grams_bytes();
        uint64_t y = a.counts_ranks_bytes();
        uint64_t z = a.pointers_bytes();
        uint64_t n = a.size();

        if (i != 1) {
            std::cout << "\tgrams: " << x << " (" << double(x) / n
                      << " per gram)\n";
        }

        std::cout << "\tranks: " << y << " (" << double(y) / n
                  << " per gram)\n";

        if (i != order()) {
            std::cout << "\tpointers: " << z << " (" << double(z) / n
                      << " per gram)" << std::endl;
        }

        ++i;
        grams_bytes += x;
        counts_ranks_bytes += y;
        pointers_bytes += z;
    }

    uint64_t num_grams = size();
    std::cout << "order: " << order() << "\n";
    std::cout << "num. of grams: " << num_grams << "\n";
    std::cout << "tot. bytes: " << bytes << "\n";
    std::cout << "bytes per gram: " << double(bytes) / num_grams << "\n";
    std::cout << "vocabulary data bytes: " << m_vocab.data_bytes() << "\n";
    std::cout << "\t(does NOT include hash function bytes)\n";
    std::cout << "unique values bytes: " << m_distinct_counts.bytes() << "\n";

    std::cout << "tot. grams bytes: " << grams_bytes << " ("
              << grams_bytes * 100.0 / bytes << "%)\n"
              << "\tper gram: " << double(grams_bytes) / num_grams << "\n";

    std::cout << "tot. ranks bytes: " << counts_ranks_bytes << " ("
              << counts_ranks_bytes * 100.0 / bytes << "%)\n"
              << "\tper gram: " << double(counts_ranks_bytes) / num_grams
              << "\n";

    std::cout << "tot. pointers bytes: " << pointers_bytes << " ("
              << pointers_bytes * 100.0 / bytes << "%)\n"
              << "\tper gram: " << double(pointers_bytes) / num_grams
              << std::endl;

    std::cout << "number of unique values:\n";
    for (uint8_t i = 0; i < m_order; ++i) {
        std::cout << "\t" << i + 1 << "-grams: " << m_distinct_counts.size(i)
                  << std::endl;
    }
    std::cout << std::endl;

    std::cout << "trie level statistics:\n";
    for (uint8_t i = 0; i < m_order; ++i) {
        uint8_t order = i + 1;
        std::cout << "\t" << uint32_t(order) << "-grams:\n";
        m_arrays[i].print_stats(order, i ? m_arrays[i - 1].ptrs() : nullptr);
        std::cout << std::endl;
    }
}

template <typename Vocabulary, typename Mapper, typename Values, typename Ranks,
          typename Grams, typename Pointers>
void trie_prob_lm<Vocabulary, Mapper, Values, Ranks, Grams,
                  Pointers>::print_stats(size_t bytes) const {
    util::logger("========= TRIE_PROB_LM statistics =========");
    uint64_t i = 1;
    uint64_t grams_bytes = 0;
    uint64_t probs_backoffs_ranks_bytes = 0;
    uint64_t pointers_bytes = 0;
    for (auto const& a : m_arrays) {
        std::cout << i << "-grams bytes:\n";
        uint64_t x = a.grams_bytes();
        uint64_t y = a.probs_backoffs_ranks_bytes();
        uint64_t z = a.pointers_bytes();
        std::cout << "\tgrams: " << x << "\n";
        std::cout << "\tprob/backoff ranks: " << y << "\n";
        std::cout << "\tpointers: " << z << std::endl;
        ++i;
        grams_bytes += x;
        probs_backoffs_ranks_bytes += y;
        pointers_bytes += z;
    }

    uint64_t num_grams = size();
    std::cout << "order: " << order() << "\n";
    std::cout << "num. of grams: " << num_grams << "\n";
    std::cout << "tot. bytes: " << bytes << "\n";
    std::cout << "bytes per gram: " << double(bytes) / num_grams << "\n";
    std::cout << "vocabulary data bytes: " << m_vocab.data_bytes() << "\n";
    std::cout << "\t(does NOT include hash function bytes)\n";
    std::cout << "probs averages bytes: " << m_probs_averages.bytes() << "\n";
    std::cout << "probs backoffs bytes: " << m_backoffs_averages.bytes()
              << "\n";

    std::cout << "tot. grams bytes: " << grams_bytes << " ("
              << grams_bytes * 100.0 / bytes << "%)\n"
              << "\tper gram: " << double(grams_bytes) / num_grams << "\n";

    std::cout << "tot. probs/backoffs bytes: " << probs_backoffs_ranks_bytes
              << " (" << probs_backoffs_ranks_bytes * 100.0 / bytes << "%)\n"
              << "\tper gram: "
              << double(probs_backoffs_ranks_bytes) / num_grams << "\n";

    std::cout << "tot. pointers bytes: " << pointers_bytes << " ("
              << pointers_bytes * 100.0 / bytes << "%)\n"
              << "\tper gram: " << double(pointers_bytes) / num_grams
              << std::endl;

    std::cout << "trie level statistics:\n";
    for (uint8_t i = 0; i < m_order; ++i) {
        uint8_t order = i + 1;
        std::cout << "\t" << uint32_t(order) << "-grams:\n";
        m_arrays[i].print_stats(order, i ? m_arrays[i - 1].ptrs() : nullptr);
        std::cout << std::endl;
    }
}

}  // namespace tongrams
