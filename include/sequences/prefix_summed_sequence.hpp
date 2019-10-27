#pragma once

namespace tongrams {

template <typename Sequence>
struct prefix_summed_sequence {
    prefix_summed_sequence() : m_size(0) {}

    template <typename Iterator>
    void build(Iterator begin, uint64_t n, uint8_t order) {
        m_size = n;
        std::vector<uint64_t> prefix_sums;
        prefix_sums.reserve(n);
        uint64_t last = 0;
        for (uint64_t i = 0; i < n; ++i, ++begin) {
            uint64_t value = *begin + last;
            prefix_sums.push_back(value);
            last = value;
        }
        m_sequence.build(prefix_sums.begin(), n, prefix_sums.back(), order);
    }

    inline uint64_t operator[](uint64_t i) {
        return m_sequence[i] - (i ? m_sequence[i - 1] : 0);
    }

    uint64_t size() const {
        return m_size;
    }

    uint64_t bytes() const {
        return sizeof(m_size) + m_sequence.bytes();
    }

    void swap(prefix_summed_sequence& other) {
        std::swap(other.m_size, m_size);
        other.m_sequence.swap(m_sequence);
    }

    void save(std::ostream& os) const {
        util::save_pod(os, &m_size);
        m_sequence.save(os);
    }

    void load(std::istream& is) {
        util::load_pod(is, &m_size);
        m_sequence.load(is);
    }

private:
    uint64_t m_size;
    Sequence m_sequence;
};

}  // namespace tongrams
