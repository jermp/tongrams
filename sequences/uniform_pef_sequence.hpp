#pragma once

#include <stdexcept>

#include "integer_codes.hpp"
#include "compact_elias_fano.hpp"

#include "../utils/util.hpp"
#include "../vectors/bit_vector.hpp"
#include "../vectors/compact_vector.hpp"

namespace pef
{
    struct uniform_pef_sequence
    {
        uniform_pef_sequence()
            : m_size(0)
            , m_universe(0)
            , m_partitions(0)
            , m_log_partition_size(0)
        {}

        template<typename Iterator>
        void build(Iterator grams_begin, uint64_t num_grams,
                   std::vector<uint64_t>& pointers, uint8_t order)
        {
            std::vector<uint64_t> values;
            values.reserve(num_grams);
            uint64_t prev_upper = 0;
            auto pointers_it = pointers.begin();
            uint64_t start = *pointers_it;
            ++pointers_it;
            uint64_t end = *pointers_it;
            uint64_t run = end - start;
            uint64_t within = 0;
            for (uint64_t i = 0; i < num_grams; ++i, ++grams_begin) {
                if (within == run) {
                    within = 0;
                    do {
                        start = end;
                        ++pointers_it;
                        end = *pointers_it;
                        run = end - start;
                    } while (!run);
                    prev_upper = values.size()
                               ? values.back() : 0;
                }
                uint64_t v = *grams_begin;
                values.push_back(v + prev_upper);
                ++within;
            }
            assert(values.size() == num_grams);
            write(values.begin(), values.back(),
                  values.size(), order);
        }

        template<typename Iterator>
        void build(Iterator begin, uint64_t n,
                   uint64_t universe, uint8_t order)
        {
            write(begin, universe, n, order);
        }

        template <typename Iterator>
        void write(Iterator begin,
                   uint64_t universe,
                   uint64_t n, uint8_t order)
        {
            assert(n > 0);

            m_size = n;
            m_universe = universe;

            if (order <= 2) {
                m_log_partition_size = 6;
            } else {
               m_log_partition_size = 7;
            }
            
            pef_global_parameters params;
            uint64_t partition_size = uint64_t(1) << m_log_partition_size;
            size_t partitions = tongrams::util::ceil_div(n, partition_size);

            m_partitions = partitions;

            tongrams::bit_vector_builder bvb;
            std::vector<uint64_t> cur_partition;

            uint64_t cur_base = 0;
            if (partitions == 1) {
                cur_base = *begin;
                Iterator it = begin;

                for (size_t i = 0; i < n; ++i, ++it) {
                    cur_partition.push_back(*it - cur_base);
                }

                uint64_t universe_bits = tongrams::util::ceil_log2(universe);
                bvb.append_bits(cur_base, universe_bits);
                // write universe only if non-singleton and not tight
                if (n > 1) {
                    if (cur_base + cur_partition.back() + 1 == universe) {
                        // tight universe
                        write_delta(bvb, 0);
                    } else {
                        write_delta(bvb, cur_partition.back());
                    }
                }

                compact_elias_fano::write(bvb, cur_partition.begin(),
                                          cur_partition.back() + 1,
                                          cur_partition.size(), params);
            } else {

                tongrams::bit_vector_builder bv_sequences;
                std::vector<uint64_t> endpoints;
                std::vector<uint64_t> upper_bounds;

                uint64_t cur_i = 0;
                Iterator it = begin;
                cur_base = *begin;
                upper_bounds.push_back(cur_base);

                for (size_t p = 0; p < partitions; ++p) {
                    cur_partition.clear();
                    uint64_t value = 0;
                    for (; cur_i < ((p + 1) * partition_size) && cur_i < n;
                         ++cur_i, ++it) {
                        value = *it;
                        cur_partition.push_back(value - cur_base);
                    }

                    assert(cur_partition.size() <= partition_size);
                    assert((p == partitions - 1)
                           || cur_partition.size() == partition_size);

                    uint64_t upper_bound = value;
                    assert(cur_partition.size() > 0);

                    compact_elias_fano::write(bv_sequences, cur_partition.begin(),
                                              cur_partition.back() + 1,
                                              cur_partition.size(), params);
                    endpoints.push_back(bv_sequences.size());
                    upper_bounds.push_back(upper_bound);
                    cur_base = upper_bound;
                }


                tongrams::compact_vector::builder
                    ub_cv_builder(upper_bounds.size(),
                                  tongrams::util::ceil_log2(upper_bounds.back() + 1));
                for (auto u: upper_bounds) {
                    ub_cv_builder.push_back(u);
                }
                m_upper_bounds.build(ub_cv_builder);

                uint64_t endpoint_bits = tongrams::util::ceil_log2(bv_sequences.size() + 1);
                write_gamma(bvb, endpoint_bits);
                for (uint64_t p = 0; p < endpoints.size() - 1; ++p) {
                    bvb.append_bits(endpoints[p], endpoint_bits);
                }

                bvb.append(bv_sequences);
                m_data.build(&bvb);
            }

            // init enumerator to map ids needed by pef_rtrie
            e.init(m_data, m_upper_bounds, m_size,
                   m_universe, m_partitions,
                   m_log_partition_size);
        }

        inline uint64_t operator[](uint64_t position) {
            return e.move(position).second;
        }

        void find(tongrams::pointer_range const& r, uint64_t id, uint64_t* pos)
        {
            assert(r.end > r.begin);
            assert(r.end <= size());

            if (!id) {
                *pos = r.begin;
                return;
            }

            uint64_t partition_begin = r.begin >> m_log_partition_size;
            uint64_t partition_end = r.end >> m_log_partition_size;

            e.switch_partition(partition_begin);

            uint64_t prev_upper = 0;
            if (LIKELY(r.begin)) {
                prev_upper = operator[](r.begin - 1);
            }
            
            id += prev_upper;
            auto pos_value = e.next_geq(id, r.end, partition_end);
            if (pos_value.second == id) {
                *pos = pos_value.first;
                return;
            }
            *pos = tongrams::global::not_found;
        }

        struct enumerator
        {
            typedef std::pair<uint64_t, uint64_t> value_type; // (position, value)

            enumerator()
            {}

            void init(tongrams::bit_vector const& bv,
                      tongrams::compact_vector const& upper_bounds,
                      uint64_t n, uint64_t universe, uint64_t partitions,
                      uint8_t log_partition_size)
            {
                m_partitions = partitions;
                m_size = n;
                m_universe = universe;
                m_bv = &bv;
                m_upper_bounds = &upper_bounds;
                m_log_partition_size = log_partition_size;

                pef_global_parameters params;
                tongrams::bits_iterator<tongrams::bit_vector> it(bv);
                if (m_partitions == 1) {
                    m_cur_partition = 0;
                    m_cur_begin = 0;
                    m_cur_end = n;

                    uint64_t universe_bits = tongrams::util::ceil_log2(universe);
                    m_cur_base = it.get_bits(universe_bits);
                    auto ub = 0;
                    if (n > 1) {
                        uint64_t universe_delta = read_delta(it);
                        ub = universe_delta ? universe_delta : (universe - m_cur_base - 1);
                    }

                    m_partition_enum =
                        compact_elias_fano::enumerator
                        (*m_bv, it.position(), ub + 1, n, params);

                    m_cur_upper_bound = m_cur_base + ub;
                } else {
                    m_endpoint_bits = read_gamma(it);
                    uint64_t cur_offset = it.position();
                    m_endpoints_offset = cur_offset;
                    uint64_t endpoints_size = m_endpoint_bits * (m_partitions - 1);
                    cur_offset += endpoints_size;
                    m_sequences_offset = cur_offset;
                }

                m_position = 0;
                m_first = true; // for next()
                slow_move();
            }

            const uint64_t linear_scan_threshold = 8;

            value_type ALWAYSINLINE move(uint64_t position)
            {
                assert(position <= size());
                m_position = position;
                if (m_position >= m_cur_begin && m_position < m_cur_end) {
                    uint64_t val = m_cur_base + m_partition_enum.move(m_position - m_cur_begin).second;
                    return value_type(m_position, val);
                }
                return slow_move();
            }

            value_type ALWAYSINLINE next_geq(uint64_t lower_bound, uint64_t range_end, uint64_t partition_end)
            {
                if (LIKELY(lower_bound >= m_cur_base && lower_bound <= m_cur_upper_bound)) {
                    auto val = m_partition_enum.next_geq(lower_bound - m_cur_base);
                    m_position = m_cur_begin + val.first;
                    return value_type(m_position, m_cur_base + val.second);
                }

                // bounds checking
                if (m_cur_partition > partition_end) { // out of bounds form the right
                    return value_type(m_position, tongrams::global::not_found);
                }

                if (lower_bound < m_cur_base) { // out of bounds form the left
                    return value_type(m_position, tongrams::global::not_found);
                }
                
                return slow_next_geq(lower_bound, range_end, partition_end);
            }

            uint64_t size() const {
                return m_size;
            }

            template<typename T>
            uint64_t bsearch(T const* vec, uint64_t lower_bound,
                             uint64_t partition_begin,
                             uint64_t partition_end)
            {
                // optimized small skips with linear scan
                if (LIKELY(partition_end - partition_begin <= linear_scan_threshold)) {
                    uint64_t id = scan(vec, partition_begin, partition_end, lower_bound);
                    return id;
                }

                uint64_t partition_id = partition_begin;
                uint64_t lo = partition_begin;
                uint64_t hi = partition_end;

                while (lo <= hi) {
                    uint64_t mid = (lo + hi) / 2;
                    uint64_t mid_value = vec->access(mid);
                    
                    if (mid_value > lower_bound) {
                        hi = mid != 0 ? mid - 1 : 0;
                        if (lower_bound > vec->access(hi)) {
                            return mid;
                        }
                    } else if (mid_value < lower_bound) {
                        lo = mid + 1;
                        if (lower_bound < (lo != vec->size() - 1 ? vec->access(lo) : vec->back())) {
                            partition_id = lo;
                        }
                    } else {
                        return mid;
                    }

                    if (hi - lo <= linear_scan_threshold) {
                        return scan(vec, lo, hi, lower_bound);
                    }
                }

                return partition_id;
            }

            template<typename T>
            uint64_t inline scan(T const* vec,
                                 uint64_t lo, uint64_t hi,
                                 uint64_t lower_bound) {
                while (lo <= hi) {
                    if (vec->access(lo) >= lower_bound) {
                        break;
                    }
                    ++lo;
                }
                return lo;
            }

            value_type NOINLINE slow_move()
            {
                if (m_position == size()) {
                    if (m_partitions > 1) {
                        switch_partition(m_partitions - 1);
                    }
                    m_partition_enum.move(m_partition_enum.size());
                    return value_type(m_position, m_universe);
                }
                uint64_t partition = m_position >> m_log_partition_size;
                switch_partition(partition);
                uint64_t val = m_cur_base + m_partition_enum.move(m_position - m_cur_begin).second;
                return value_type(m_position, val);
            }

            value_type NOINLINE slow_next_geq(uint64_t lower_bound, uint64_t range_end, uint64_t partition_end)
            {
                if (m_partitions == 1) {
                    if (lower_bound < m_cur_base) {
                        return move(0);
                    } else {
                        return move(size());
                    }
                }

                uint64_t partition_id = bsearch(m_upper_bounds, lower_bound, m_cur_partition + 1, partition_end + 1);

                if (partition_id == 0) {
                    return move(0);
                }

                if (partition_id == m_upper_bounds->size()) {
                    return move(size());
                }

                switch_partition(partition_id - 1);
                return next_geq(lower_bound, range_end, partition_end);
            }

            // used for linear scan in stats.cpp
            uint64_t next() {
                if (UNLIKELY(m_first)) {
                    uint64_t offset = m_partition_enum.move(0).second;
                    m_first = false;
                    return m_cur_base + offset;
                }

                ++m_position;
                if (LIKELY(m_position < m_cur_end)) {
                    uint64_t offset = m_partition_enum.next().second;
                    return m_cur_base + offset;
                }
                return slow_next();
            }

            // used for linear scan in stats.cpp
            uint64_t slow_next() {
                if (UNLIKELY(m_position == m_size)) {
                    assert(m_cur_partition == m_partitions - 1);
                    auto val = m_partition_enum.next();
                    assert(val.first == m_partition_enum.size()); (void)val;
                    return m_universe;
                }
                switch_partition(m_cur_partition + 1);
                uint64_t val = m_cur_base + m_partition_enum.move(0).second;
                return val;
            }

            void switch_partition(uint64_t partition)
            {
                assert(m_partitions > 1);

                uint64_t endpoint = partition
                    ? m_bv->get_bits(m_endpoints_offset +
                                     (partition - 1) * m_endpoint_bits,
                                     m_endpoint_bits)
                    : 0;

                m_cur_partition_begin = m_sequences_offset + endpoint;
                tongrams::util::prefetch(m_bv->data().data() + m_cur_partition_begin / 64);

                m_cur_partition = partition;
                m_cur_begin = partition << m_log_partition_size;
                m_cur_end = std::min(size(), (partition + 1) << m_log_partition_size);

                m_cur_upper_bound = m_upper_bounds->access(partition + 1);
                m_cur_base = m_upper_bounds->access(partition);

                m_partition_enum = compact_elias_fano::enumerator
                    (*m_bv, m_cur_partition_begin,
                     m_cur_upper_bound - m_cur_base + 1,
                     m_cur_end - m_cur_begin, m_params);
            }

            uint8_t m_log_partition_size;

            pef_global_parameters m_params;
            uint64_t m_partitions;
            uint64_t m_endpoints_offset;
            uint64_t m_endpoint_bits;
            uint64_t m_sequences_offset;
            uint64_t m_size;
            uint64_t m_universe;

            uint64_t m_position;
            uint64_t m_cur_partition_begin;
            uint64_t m_cur_partition;
            uint64_t m_cur_begin;
            uint64_t m_cur_end;
            uint64_t m_cur_base;
            uint64_t m_cur_upper_bound;

            tongrams::bit_vector const* m_bv;
            compact_elias_fano::enumerator m_partition_enum;
            tongrams::compact_vector const* m_upper_bounds;

            bool m_first;
        };

        // used for linear scan in stats.cpp
        enumerator begin() const {
            enumerator e;
            e.init(m_data, m_upper_bounds, m_size,
                   m_universe, m_partitions,
                   m_log_partition_size);
            return e;
        }

        uint64_t size() const {
            return m_size;
        }

        uint64_t universe() const {
            return m_universe;
        }

        uint64_t num_partitions() const {
            return m_partitions;
        }

        uint64_t bytes() const {
            return sizeof(m_size)
                 + sizeof(m_universe)
                 + sizeof(m_partitions)
                 + m_upper_bounds.bytes()
                 + m_data.bytes()
                 + sizeof(m_log_partition_size);
        }

        void save(std::ostream& os) const {
            tongrams::util::save_pod(os, &m_size);
            tongrams::util::save_pod(os, &m_universe);
            tongrams::util::save_pod(os, &m_partitions);
            m_upper_bounds.save(os);
            m_data.save(os);
            tongrams::util::save_pod(os, &m_log_partition_size);
        }

        void load(std::istream& is) {
            tongrams::util::load_pod(is, &m_size);
            tongrams::util::load_pod(is, &m_universe);
            tongrams::util::load_pod(is, &m_partitions);
            m_upper_bounds.load(is);
            m_data.load(is);
            tongrams::util::load_pod(is, &m_log_partition_size);

            if (m_size) {
                e.init(m_data, m_upper_bounds, m_size,
                       m_universe, m_partitions,
                       m_log_partition_size);
            }
        }

    private:
        uint64_t m_size;
        uint64_t m_universe;
        uint64_t m_partitions;
        pef_global_parameters m_params;
        tongrams::compact_vector m_upper_bounds;
        tongrams::bit_vector m_data;
        uint8_t m_log_partition_size;
        enumerator e;
    };
}
