#pragma once

namespace tongrams {

template <typename Pointers>
struct pointer_sequence {
    pointer_sequence() {}

    template <typename T>
    void build(T& pointers) {
        m_pointers.build(pointers.begin(), pointers.size(), pointers.back());
    }

    inline pointer_range operator[](uint64_t i) {
        auto p = m_pointers.pair(i);
        return {p.first, p.second};
    }

    uint64_t size() const {
        return m_pointers.size();
    }

    uint64_t universe() {
        return m_pointers.universe();
    }

    void save(std::ostream& os) const {
        m_pointers.save(os);
    }

    void load(std::istream& is) {
        m_pointers.load(is);
    }

    uint64_t bytes() const {
        return m_pointers.bytes();
    }

    // proxy of Pointer::iterator,
    // used in sorted_array.print_stats()
    struct iterator {
        iterator(Pointers const& ptrs)
            : m_it(ptrs)  // always begins at 0
        {}

        uint64_t next() {
            return m_it.next();
        }

    private:
        typename Pointers::iterator m_it;
    };

    iterator begin() const {
        return iterator(m_pointers);
    }

private:
    Pointers m_pointers;
};

}  // namespace tongrams
