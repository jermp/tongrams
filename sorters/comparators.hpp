#pragma once

#include "../utils/iterators.hpp"

namespace tongrams {

    template<typename Vocabulary,
             typename Record,
             typename Iterator>
    struct comparator {
        comparator(Vocabulary const& vocab)
            : m_vocab(&vocab)
        {}

        bool operator()(Record const& x, Record const& y)
        {
            m_x_it.init(x.gram);
            m_y_it.init(y.gram);
            uint32_t order = m_x_it.spaces() + 1;

            for (uint32_t i = 0; i < order; ++i)
            {
                auto br_x = m_x_it.next();
                auto br_y = m_y_it.next();

                // empty strings
                if (br_x.first == br_x.second
                 || br_y.first == br_y.second) {
                    return false;
                }

                if (!bytes::equal_bytes(br_x, br_y))
                {
                    uint64_t id_x = m_vocab->lookup(br_x, identity_adaptor());
                    if (id_x == global::not_found) {
                        return false;
                    }

                    uint64_t id_y = m_vocab->lookup(br_y, identity_adaptor());
                    if (id_y == global::not_found) {
                        return false;
                    }

                    return id_x < id_y;
                }
            }

            return false;
        }

    private:
        Vocabulary const* m_vocab;
        Iterator m_x_it, m_y_it;
    };

    #define prefix_order_comparator(Vocabulary, Record) comparator<                             \
                                                                   Vocabulary,                  \
                                                                   Record,                      \
                                                                   forward_byte_range_iterator  \
                                                                  >                             \

    #define suffix_order_comparator(Vocabulary, Record) comparator<                             \
                                                                   Vocabulary,                  \
                                                                   Record,                      \
                                                                   backward_byte_range_iterator \
                                                                  >
}
