#pragma once

namespace tongrams
{
    template<typename T>
    struct circular_buffer {
        circular_buffer(size_t capacity)
            : m_capacity(capacity)
            , m_begin(0)
            , m_end(0)
        {
            m_data.resize(capacity);
        }

        void init() {
            m_begin = 0;
            m_end = 0;
        }

        inline void push_back(T x) {
            m_data[m_end] = x;

            ++m_end;
            fall_back(m_end);

            if (m_end == m_begin) {
                ++m_begin;
                fall_back(m_begin);
            }
        }

        struct reverse_iterator {
            reverse_iterator(circular_buffer<T> const* cb)
                : m_cb(cb)
                , m_pos(cb->m_end)
            {
                // one step back since
                // m_end is one past the end
                this->operator++();
            }

            inline T operator*() const {
                return m_cb->m_data[m_pos];
            }

            inline void operator++() {
                if (!m_pos) {
                    m_pos = m_cb->m_capacity;
                }
                --m_pos;
            }

        private:
            circular_buffer const* m_cb;
            uint64_t m_pos;
        };

        reverse_iterator rbegin() const {
            return reverse_iterator(this);
        }

    private:
        uint64_t m_capacity;
        uint64_t m_begin;
        uint64_t m_end;
        std::vector<T> m_data;

        void fall_back(uint64_t& x) const {
            if (x == m_capacity) {
                x = 0;
            }
        }
    };
}
