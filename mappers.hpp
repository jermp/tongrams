#pragma once

namespace tongrams
{
    struct identity_mapper
    {
        static const bool context_remapping = false;

        template<typename SortedArray> static inline
        uint64_t map_id(uint64_t /*parent_id*/,
                        uint64_t /*parent_parent_id*/,
                        uint64_t id,
                        SortedArray* /*unigrams*/,
                        uint8_t /*remapping_order*/)
        {
            return id;
        }

        template<typename Vocabulary,
                 typename SortedArray> static inline
        uint64_t map_id(byte_range /*gram*/, uint64_t id,
                        Vocabulary const* /*vocab*/,
                        SortedArray const* /*unigrams*/,
                        uint8_t /*remapping_order*/,
                        bool /*reversed*/)
        {
            return id;
        }

        template<typename Vocabulary,
                 typename SortedArray> static inline
        uint64_t map_query(byte_range range, uint64_t* word_ids,
                           Vocabulary const* vocab,
                           SortedArray const* /*unigrams*/,
                           uint8_t /*remapping_order*/)
        {
            uint8_t const* pos = range.first;
            uint8_t const* prev_pos = pos;
            identity_adaptor adaptor;
            uint64_t order_m1 = 0; // order minus 1
            for (; pos != range.second; ++pos) {
                // assume words separated by whitespaces
                if (*pos == ' ') {

                    // the token to lookup could be a whitespace ' '
                    if (pos == range.first) {
                        ++pos;
                        break;
                    }

                    byte_range br(prev_pos, pos);
                    uint64_t id = 0;
                    vocab->lookup(br, id, adaptor);
                    word_ids[order_m1] = id;
                    prev_pos = pos + 1;
                    ++order_m1;
                }

            }

            uint64_t id = 0;
            byte_range br(prev_pos, pos);
            vocab->lookup(br, id, adaptor);
            word_ids[order_m1] = id;
            return order_m1;
        }
    };

    struct sorted_array_mapper
    {
        static const bool context_remapping = true;

        // NOTE: always works in BACKWARD direction
        // because it is used by trie_prob_lm::score()
        template<typename SortedArray> static inline
        uint64_t map_id(uint64_t parent_id,
                        uint64_t parent_parent_id,
                        uint64_t id,
                        SortedArray* unigrams,
                        uint8_t remapping_order)
        {
            pointer_range r;
            if (remapping_order == 1) {
                r = unigrams->range(parent_id);
            } else {
                r = unigrams->range(parent_parent_id);
                (unigrams + 1)->next(r, parent_id);
            }
            return map_id(r, id, unigrams + remapping_order);
        }

        template<typename Vocabulary,
                 typename SortedArray> static inline
        uint64_t map_id(byte_range gram, uint64_t id,
                        Vocabulary const* vocab,
                        SortedArray* unigrams,
                        uint8_t remapping_order,
                        bool reversed) // controls the order of remapping
                                       // example, assuming reversed is false:
                                       // given 'X Y Z'
                                       // order = 1 --> Z is remapped as the position it occupies within the successors of 'Y'
                                       // order = 2 --> Z is remapped as the position it occupies within the successors of 'X Y'
        {
            byte_range begin;
            uint64_t parent_id = 0;
            if (reversed) {
                begin = bytes::to(gram, remapping_order);
            } else {
                begin = bytes::back_to(gram, remapping_order);
            }
            vocab->lookup(begin, parent_id, identity_adaptor());
            auto r = unigrams->range(parent_id);

            if (remapping_order == 2) {
                byte_range next_word;
                uint64_t next_id = 0;
                if (reversed) {
                    next_word = bytes::prev(begin);
                } else {
                    next_word = bytes::next(begin);
                }
                vocab->lookup(next_word, next_id, identity_adaptor());
                (unigrams + 1)->next(r, next_id);
            }
            return map_id(r, id, unigrams + remapping_order);
        }

        template<typename SortedArray> static inline
        uint64_t map_id(pointer_range r, uint64_t id,
                        SortedArray* ngrams)
        {
            uint64_t pos = 0;
            auto grams = ngrams->grams();
            grams->find(r, id, &pos);
            if (pos == global::not_found) {
                return global::not_found;
            }
            assert(pos >= r.begin && pos < r.end);
            return pos - r.begin;
        }

        // NOTE: only works in FORWARD direction
        // as it is used in trie_count_lm::lookup()
        template<typename Vocabulary,
                 typename SortedArray> static inline
        uint64_t map_query(byte_range range, uint64_t* word_ids,
                           Vocabulary const* vocab,
                           SortedArray* unigrams,
                           uint8_t remapping_order)
        {
            uint8_t const* pos = range.first;
            uint8_t const* prev_pos = pos;
            identity_adaptor adaptor;
            uint64_t order_m1 = 0; // order minus 1

            // parent_ids are NOT remapped
            // keep parent_ids aligned with word_ids
            static uint64_t parent_ids[global::max_order];

            for (; pos != range.second; ++pos) {
                // assume words separated by whitespaces
                if (*pos == ' ') {

                    // the token to lookup could be a whitespace ' '
                    if (pos == range.first) {
                        ++pos;
                        break;
                    }

                    byte_range br(prev_pos, pos);
                    uint64_t id = 0;
                    vocab->lookup(br, id, adaptor);
                    word_ids[order_m1] = id;
                    prev_pos = pos + 1;
                    ++order_m1;
                    parent_ids[order_m1] = id;

                    if (order_m1 == 2) { // 3-4-5-gram case
                        ++pos;

                        for (; pos != range.second; ++pos) {
                            if (*pos == ' ') {
                                uint64_t id = 0;
                                uint64_t mapped_id = 0;
                                byte_range br(prev_pos, pos);
                                vocab->lookup(br, id, adaptor);

                                if (order_m1 == 2 && remapping_order == 2) {
                                    mapped_id = id;
                                } else {
                                    mapped_id = map_id(id, parent_ids,
                                                       unigrams, order_m1,
                                                       remapping_order);
                                }

                                word_ids[order_m1] = mapped_id;
                                prev_pos = pos + 1;
                                ++order_m1;
                                parent_ids[order_m1] = id;
                            }
                        }

                        uint64_t mapped_id = 0;
                        byte_range br(prev_pos, pos);
                        vocab->lookup(br, mapped_id, adaptor);

                        // if we have a trigram and remapping order is 2
                        // just pick vocab id
                        if (order_m1 != 2 || remapping_order != 2) {
                            mapped_id = map_id(mapped_id, parent_ids,
                                               unigrams, order_m1,
                                               remapping_order);
                        }

                        word_ids[order_m1] = mapped_id;
                        return order_m1;
                    }
                }
            }

            // unigram or bigram case
            uint64_t id = 0;
            byte_range br(prev_pos, pos);
            vocab->lookup(br, id, adaptor);
            word_ids[order_m1] = id;
            return order_m1;
        }

    private:
        template<typename SortedArray> static inline
        uint64_t map_id(uint64_t id, uint64_t* parent_ids,
                        SortedArray* unigrams,
                        uint8_t order_m1,
                        uint8_t remapping_order)
        {
            uint64_t parent_id = parent_ids[order_m1 - remapping_order + 1];
            auto r = unigrams->range(parent_id);
            if (remapping_order == 2) {
                (unigrams + 1)->next(r, parent_ids[order_m1]);
            }
            return map_id(r, id, unigrams + remapping_order);
        }
    };
}
