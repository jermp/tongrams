#pragma once

#include "utils/util.hpp"
#include "sequences/darray.hpp"
#include "vectors/bit_vector.hpp"
#include "../external/emphf/common.hpp"
#include "../external/emphf/mmap_memory_model.hpp"
#include "../external/emphf/hypergraph_sorter_scan.hpp"
#include "../external/emphf/base_hash.hpp"
#include "utils/mphf.hpp"
#include "vectors/compact_vector.hpp"
#include "utils/mph_tables.hpp"

#include <numeric>

namespace tongrams {

struct fast_ef_sequence {
    typedef uint64_t sample_t;
    static const uint64_t sampling_threshold = 128;
    static const uint64_t log2_sampling_threshold = 7;
    static const uint64_t linear_scan_threshold = 64;

    fast_ef_sequence() : m_size(0), m_l(0) {}

    template <typename Iterator>
    void build(Iterator grams_begin, uint64_t num_grams,
               std::vector<uint64_t> const& pointers, uint8_t /*order*/) {
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
                prev_upper = values.size() ? values.back() : 0;
            }
            uint64_t v = *grams_begin;
            values.push_back(v + prev_upper);
            ++within;
        }
        assert(values.size() == num_grams);
        build(values.begin(), values.size(), values.back(), pointers);
    }

    inline uint64_t operator[](uint64_t i) const {
        return ((m_high_bits_d1.select(m_high_bits, i) - i) << m_l) |
               m_low_bits.get_bits(i * m_l, m_l);
    }

    inline uint64_t num_ones() const {
        return m_high_bits_d1.num_positions();
    }

    void find(pointer_range const& r, uint64_t id, uint64_t* pos) {
        assert(r.end > r.begin);
        assert(r.end <= size());

        // NOTE: switch to this code if samplings ARE NOT to be used
        // uint64_t prev_upper = r.begin
        //                     ? operator[](r.begin - 1) : 0;
        // id += prev_upper;
        // uint64_t lo = r.begin;
        // uint64_t hi = r.end;

        // while (lo <= hi) {
        //     if (LIKELY(hi - lo <= linear_scan_threshold)) {
        //         scan(lo, hi, id, pos);
        //         return;
        //     }
        //     uint64_t mid = (lo + hi) >> 1;
        //     uint64_t v = operator[](mid);
        //     if (v == id) {
        //         *pos = mid;
        //         return;
        //     }
        //     if (v > id) {
        //         // NOTE: not mid-1 since we are searching in [lo, hi)
        //         hi = mid;
        //     } else {
        //         lo = mid + 1;
        //     }

        //     // NOTE: check the invariant that id must be searched in [lo, hi)
        //     // notice that the assert must be commented when searching for
        //     // an id that could NOT be present, as it is for the perplexity
        //     benchmark
        //     // assert(id >= operator[](lo)
        //     //     && id <= operator[](hi - 1) // NOTE: not id <
        //     operator[](hi)
        //     //                                 // because the sequence could
        //     NOT be strictly increasing
        //     //                                 // and it could be that
        //     operator[](hi) == operator[](hi - 1)
        //     //     );
        // }
        // *pos = global::not_found;

        // NOTE: switch to this code if samplings ARE to be used
        size_t lo = r.begin;
        size_t hi = r.end;
        size_t run = hi - lo;

        if (run >= sampling_threshold) {
            size_t base = m_offsets.lookup(lo, uint64_adaptor());
            uint64_t prev_upper = m_samplings[base];
            id += prev_upper;
            uint64_t lower_bound = prev_upper;
            uint64_t tree_height =
                util::ceil_log2(run) - log2_sampling_threshold;

            for (size_t i = 1, h = 0; h != tree_height; ++h) {
                sample_t sample = m_samplings[base + i];
                size_t mid = (lo + hi) >> 1;
                assert(sample == operator[](mid));
                if (id == sample) {
                    *pos = mid;
                    return;
                }
                if (id < sample) {
                    hi = mid;
                    i = (i << 1);
                } else {
                    lo = mid + 1;
                    lower_bound = sample;
                    i = (i << 1) + 1;
                }
            }

            bsearch_scan(lo, hi, id, pos, lower_bound);
            return;
        }

        uint64_t prev_upper = lo ? operator[](lo - 1) : 0;
        bsearch_scan(lo, hi, id + prev_upper, pos, prev_upper);
    }

    struct iterator {
        iterator(bit_vector const& high_bits, bit_vector const& low_bits,
                 darray1 const& high_bits_d1, uint8_t l, uint64_t i = 0)
            : m_high_bits(&high_bits)
            , m_low_bits(&low_bits)
            , m_high_bits_d1(&high_bits_d1)
            , m_l(l)
            , m_i(i) {
            m_low_mask = (uint64_t(1) << m_l) - 1;
            m_low_buf = 0;
            if (m_l) {
                m_chunks_in_word = 64 / m_l;
                m_chunks_avail = 0;
            } else {
                m_chunks_in_word = 0;
                m_chunks_avail = m_high_bits_d1->num_positions();
            }

            if (!(m_high_bits_d1->num_positions())) {
                return;
            }
            uint64_t pos = m_high_bits_d1->select(*m_high_bits, m_i);
            m_high_enum = bit_vector::unary_iterator(*m_high_bits, pos);
            assert(m_l < 64);
        }

        uint64_t next() {
            if (!m_chunks_avail--) {
                m_low_buf = m_low_bits->get_word64(m_i * m_l);
                m_chunks_avail = m_chunks_in_word - 1;
            }

            uint64_t high = m_high_enum.next();
            assert(high == m_high_bits_d1->select(*m_high_bits, m_i));
            uint64_t low = m_low_buf & m_low_mask;
            uint64_t ret = (((high - m_i) << m_l) | low);
            ++m_i;
            m_low_buf >>= m_l;

            return ret;
        }

    private:
        bit_vector const* m_high_bits;
        bit_vector const* m_low_bits;
        darray1 const* m_high_bits_d1;
        uint64_t m_l;
        uint64_t m_i;
        bit_vector::unary_iterator m_high_enum;
        uint64_t m_low_buf;
        uint64_t m_low_mask;
        uint64_t m_chunks_in_word;
        uint64_t m_chunks_avail;
    };

    // used for linear scan in stats.cpp
    iterator begin() const {
        return iterator(m_high_bits, m_low_bits, m_high_bits_d1, m_l);
    }

    uint64_t size() const {
        return m_size;
    }

    uint64_t universe() const {
        return operator[](size() - 1);
    }

    uint64_t bytes() const {
        return sizeof(m_size) + m_high_bits.bytes() + m_high_bits_d1.bytes() +
               m_low_bits.bytes() + sizeof(m_l) + samplings_bytes();
    }

    uint64_t samplings_bytes() const {
        return m_samplings.size() * sizeof(sample_t) + m_offsets.data_bytes();
    }

    void swap(fast_ef_sequence& other) {
        std::swap(other.m_size, m_size);
        other.m_offsets.swap(m_offsets);
        other.m_samplings.swap(m_samplings);
        other.m_high_bits.swap(m_high_bits);
        other.m_high_bits_d1.swap(m_high_bits_d1);
        other.m_low_bits.swap(m_low_bits);
        std::swap(other.m_l, m_l);
    }

    void save(std::ostream& os) const {
        essentials::save_pod(os, m_size);
        m_offsets.save(os);
        essentials::save_vec(os, m_samplings);
        m_high_bits.save(os);
        m_high_bits_d1.save(os);
        m_low_bits.save(os);
        essentials::save_pod(os, m_l);
    }

    void load(std::istream& is) {
        essentials::load_pod(is, m_size);
        m_offsets.load(is);
        essentials::load_vec(is, m_samplings);
        m_high_bits.load(is);
        m_high_bits_d1.load(is);
        m_low_bits.load(is);
        essentials::load_pod(is, m_l);
    }

private:
    uint64_t m_size;
    uint_mpht<uint64_t, uint64_t> m_offsets;
    std::vector<sample_t> m_samplings;
    bit_vector m_high_bits;
    darray1 m_high_bits_d1;
    bit_vector m_low_bits;
    uint8_t m_l;

    // iterator-based scan
    void scan(uint64_t lo, uint64_t hi, uint64_t id, uint64_t* pos) const {
        iterator it(m_high_bits, m_low_bits, m_high_bits_d1, m_l, lo);
        size_t run = hi - lo;
        for (uint64_t i = 0; i < run; ++i) {
            auto value = it.next();
            if (value == id) {
                *pos = lo + i;
                return;
            }
        }
        *pos = global::not_found;
    }

    // optimized scan
    void scan(uint64_t lo, uint64_t hi, uint64_t id, uint64_t* pos,
              uint64_t lower_bound) const {
        // STEP (1): check high part only
        uint64_t begin = lo                       // # of 1's
                         + (lower_bound >> m_l);  // # of 0's

        bit_vector::unary_iterator it(m_high_bits, begin);
        uint64_t high_id = id >> m_l;
        uint64_t v = it.next();

        while (v - lo < high_id) {
            ++lo;
            if (lo == hi) {
                *pos = global::not_found;
                return;
            }
            v = it.next();
        }

        assert(v - lo >= high_id);

        if (v - lo > high_id) {
            *pos = global::not_found;
            return;
        }

        assert(v - lo == high_id);

        // STEP (2): check low part only
        uint64_t chunks_in_word;
        uint64_t chunks_avail;
        if (m_l) {
            chunks_in_word = 64 / m_l;
            chunks_avail = 0;
        } else {
            chunks_in_word = 0;
            chunks_avail = num_ones();
        }
        uint64_t low_mask = (uint64_t(1) << m_l) - 1;
        uint64_t low_id = id & low_mask;
        uint64_t word = 0;  // just to silence warning

        it.skip0(1);
        uint64_t run =
            std::min(it.position() - v  // # of integers sharing the same high
                                        // part identified during STEP (1)
                     ,
                     hi - lo);

        for (size_t i = 0; i < run; ++i, ++lo) {
            if (!chunks_avail--) {
                word = m_low_bits.get_word64(lo * m_l);
                chunks_avail = chunks_in_word - 1;
            }
            uint64_t cur_low = word & low_mask;

            if (cur_low == low_id) {
                *pos = lo;
                assert(operator[](lo) == id);
                return;
            }

            if (cur_low > low_id) {
                break;
            }

            word >>= m_l;
        }

        *pos = global::not_found;
    }

    void bsearch_scan(uint64_t lo, uint64_t hi, uint64_t id, uint64_t* pos,
                      uint64_t lower_bound) const {
        while (hi - lo > linear_scan_threshold) {
            size_t mid = (lo + hi) >> 1;
            uint64_t v = operator[](mid);
            if (id == v) {
                *pos = mid;
                return;
            }
            if (id < v) {
                hi = mid;
            } else {
                lo = mid + 1;
                lower_bound = v;
            }
        }

        return scan(lo, hi, id, pos, lower_bound);
    }

    template <typename T>
    struct deque_t {
        deque_t() : m_low(0) {}

        void reserve(size_t n) {
            m_data.reserve(n);
        }

        T& pop_front() {
            return m_data[m_low++];
        }

        size_t size() const {
            return m_data.size() - m_low;
        }

        void push_back(T const& obj) {
            m_data.push_back(obj);
        }

    private:
        size_t m_low;
        std::vector<T> m_data;
    };

    template <typename RandomAccessIterator>
    static void fill_samplings(uint64_t lo, uint64_t hi, uint64_t tree_height,
                               RandomAccessIterator it,
                               std::vector<sample_t>& samplings) {
        typedef std::pair<uint64_t, uint64_t> range_t;
        range_t root(lo, hi);
        deque_t<range_t> ranges;  // BFS-order
        ranges.reserve((size_t(1) << tree_height) - 1);
        ranges.push_back(root);
        for (uint64_t h = 0; h != tree_height; ++h) {
            uint64_t nodes = uint64_t(1) << h;
            for (uint64_t i = 0; i < nodes; ++i) {
                auto r = ranges.pop_front();
                uint64_t l = r.first;
                uint64_t h = r.second;
                size_t mid = (l + h) >> 1;
                samplings.push_back(it[mid]);
                range_t left(l, mid);
                range_t right(mid + 1, h);
                ranges.push_back(left);
                ranges.push_back(right);
            }
        }
    }

    template <typename RandomAccessIterator>
    void build(RandomAccessIterator begin, uint64_t n, uint64_t u,
               std::vector<uint64_t> const& pointers) {
        {
            std::vector<uint64_t> from;
            std::vector<uint64_t> to;
            std::vector<sample_t> samplings;

            uint64_t ptr_begin = pointers.front();
            for (uint64_t i = 1; i < pointers.size(); ++i) {
                uint64_t ptr_end = pointers[i];
                uint64_t range = ptr_end - ptr_begin;
                if (range >= sampling_threshold) {
                    from.push_back(ptr_begin);
                    to.push_back(samplings.size());
                    uint64_t tree_height =
                        util::ceil_log2(range) - log2_sampling_threshold;
                    // push previous range upperbound
                    samplings.push_back(ptr_begin ? begin[ptr_begin - 1] : 0);
                    fill_samplings(ptr_begin, ptr_end, tree_height, begin,
                                   samplings);
                }
                ptr_begin = ptr_end;
            }
            m_samplings.swap(samplings);

            if (from.size()) {
                m_offsets.build(from, to, uint64_adaptor());
            }
        }

        m_size = n;
        m_l = uint8_t((n && u / n) ? util::msb(u / n) : 0);
        bit_vector_builder bvb_high_bits(n + (u >> m_l) + 1);
        bit_vector_builder bvb_low_bits;
        bvb_low_bits.reserve(n * m_l);

        uint64_t low_mask = (uint64_t(1) << m_l) - 1;
        uint64_t last = 0;
        for (size_t i = 0; i < n; ++i, ++begin) {
            uint64_t v = *begin;
            if (i && v < last) {  // diagnostic purpose only
                std::cerr << "error at position " << i << "/" << n << std::endl;
                std::cerr << v << " < " << last << std::endl;
                throw std::runtime_error("sequence is not sorted");
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
    }
};
}  // namespace tongrams
