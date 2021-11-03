#pragma once

#include <unistd.h>

#include "utils/mph_tables.hpp"
#include "utils/parsers.hpp"
#include "utils/util.hpp"

namespace tongrams {

template <typename Values, typename KeyRankSequence, typename BaseHasher>
struct mph_count_lm {
    typedef single_valued_mpht<KeyRankSequence, BaseHasher> hash_table_type;

    mph_count_lm() : m_order(0) {}

    mph_count_lm(const char* input_dir, uint8_t order) : m_order(order) {
        building_util::check_order(m_order);
        m_tables.reserve(m_order);

        typename Values::builder counts_builder(m_order);

        for (uint8_t ord = 1; ord <= m_order; ++ord) {
            std::string filename;
            util::input_filename(input_dir, ord, filename);
            util::check_filename(filename);
            grams_gzparser gp(filename.c_str());

            std::vector<uint64_t> counts;
            counts.reserve(gp.num_lines());

            essentials::logger("Reading " + std::to_string(ord) +
                               "-grams counts");
            for (auto const& l : gp) {
                counts_builder.eat_value(l.count);
            }

            counts_builder.build_sequence();
        }

        size_t available_ram = sysconf(_SC_PAGESIZE) * sysconf(_SC_PHYS_PAGES);
        for (uint8_t ord = 1; ord <= m_order; ++ord) {
            essentials::logger("Building " + std::to_string(ord) + "-grams");
            grams_counts_pool pool(available_ram * 0.8);
            std::string filename;
            util::input_filename(input_dir, ord, filename);
            pool.load_from<grams_gzparser>(filename.c_str());

            auto& pool_index = pool.index();
            uint64_t n = pool_index.size();
            compact_vector::builder counts_ranks_cvb(
                n, util::ceil_log2(counts_builder.size(ord - 1) + 1));

            std::vector<byte_range> byte_ranges;
            byte_ranges.reserve(n);
            for (auto const& record : pool_index) {
                byte_ranges.push_back(record.gram);
                uint64_t rank = counts_builder.rank(ord - 1, record.count);
                counts_ranks_cvb.push_back(rank);
            }

            typename hash_table_type::builder builder(
                byte_ranges, compact_vector(counts_ranks_cvb),
                identity_adaptor());
            m_tables.emplace_back(builder);
        }

        counts_builder.build(m_distinct_counts);
    }

    template <typename T, typename Adaptor>
    uint64_t lookup(T gram, Adaptor adaptor) const {
        byte_range br = adaptor(gram);
        auto order = std::count(br.first, br.second, ' ');  // order minus 1
        assert(order < m_order);
        uint64_t rank = m_tables[order].lookup(gram, adaptor);
        return m_distinct_counts.access(order, rank);
    }

    void print_stats(size_t bytes) const;

    uint64_t order() const {
        return uint64_t(m_order);
    }

    size_t size() const {
        size_t size = 0;
        for (auto const& t : m_tables) {
            size += t.size();
        }
        return size;
    }

    void save(std::ostream& os) const {
        essentials::save_pod(os, m_order);
        m_distinct_counts.save(os);
        for (auto const& t : m_tables) {
            t.save(os);
        }
    }

    void load(std::istream& is) {
        essentials::load_pod(is, m_order);
        m_distinct_counts.load(is, m_order);
        m_tables.resize(m_order);
        for (auto& t : m_tables) {
            t.load(is);
        }
    }

private:
    uint8_t m_order;
    Values m_distinct_counts;
    std::vector<hash_table_type> m_tables;
};

}  // namespace tongrams
