#include <iostream>

#include "utils/util.hpp"
#include "sequences/uniform_pef_sequence.hpp"
#include "../external/essentials/include/essentials.hpp"
#include "../external/cmd_line_parser/include/parser.hpp"

int main(int argc, char** argv) {
    using namespace tongrams;
    cmd_line_parser::parser parser(argc, argv);
    parser.add("num_of_values", "Number of values.");
    parser.add("max_range_len", "Maximum range length.");
    if (!parser.parse()) return 1;

    uint64_t n = parser.get<uint64_t>("num_of_values");
    uint64_t max_range_len = parser.get<uint64_t>("max_range_len");
    if (n == 0 or max_range_len == 0) {
        std::cerr << "Arguments must be both non zero." << std::endl;
        std::terminate();
    }

    std::random_device rd;
    std::mt19937 rng(rd());

    {
        essentials::uniform_int_rng<uint64_t> value_distr(
            0,  // sequence is not strictly increasing
            100, essentials::get_random_seed());

        std::vector<uint64_t> values;
        values.reserve(n);
        uint64_t last_value = 0;
        for (uint64_t i = 0; i < n; ++i) {
            values.push_back(last_value + value_distr.gen());
            last_value = values.back();
        }
        assert(values.size() == n);

        pef::uniform_pef_sequence seq;
        essentials::logger("Building sequence");
        seq.build(values.begin(), values.size(), values.back(),
                  2);  // partition size is 64

        essentials::logger("Testing iterator");
        auto it = seq.begin();
        for (uint64_t i = 0; i < n; ++i) {
            util::check(i, it.next(), values[i], "value");
        }
        essentials::logger("OK");

        essentials::logger("Testing uniform_pef_sequence::operator[]()");
        for (uint64_t i = 0; i < n; ++i) {
            util::check(i, seq[i], values[i], "value");
        }
        essentials::logger("OK");
    }

    std::uniform_int_distribution<uint64_t> range_distr(1,  // not empty ranges
                                                        max_range_len);

    std::vector<uint64_t> pointers;
    pointers.push_back(0);
    std::vector<pointer_range> pointer_ranges;
    uint64_t last_offset = 0;
    for (uint64_t offset = 0; offset < n;) {
        uint64_t range = range_distr(rng);
        offset += range;
        if (offset > n) {
            offset = n;
        }
        pointers.push_back(offset);
        pointer_ranges.push_back({last_offset, offset});
        last_offset = offset;
    }
    assert(pointers.size() == pointer_ranges.size() + 1);
    std::cout << "number of created ranges: " << pointer_ranges.size()
              << std::endl;

    essentials::uniform_int_rng<uint64_t> value_distr(
        1,  // elements are distinct within a range
        50, essentials::get_random_seed());

    std::vector<uint64_t> values;
    values.reserve(n);
    for (auto const& ptr_range : pointer_ranges) {
        uint64_t last_value = 0;
        values.push_back(
            last_value);  // make first element of each range equal to 0
        for (uint64_t k = ptr_range.begin + 1; k < ptr_range.end; ++k) {
            values.push_back(last_value + value_distr.gen());
            last_value = values.back();
        }
    }
    assert(values.size() == n);

    pef::uniform_pef_sequence seq;
    essentials::logger("Building sequence");
    seq.build(values.begin(), values.size(), pointers,
              2);  // partition size is 64

    essentials::logger("Testing iterator");
    auto it = seq.begin();
    uint64_t j = 0;
    for (uint64_t i = 0; i < n; ++i) {
        auto ptr_range = pointer_ranges[j];
        if (i >= ptr_range.end) {
            ++j;
            ptr_range = pointer_ranges[j];
        }
        util::check(i, it.next() - (ptr_range.begin ? seq[ptr_range.begin] : 0),
                    values[i], "value");
    }
    assert(j == pointer_ranges.size() - 1);
    essentials::logger("OK");

    essentials::logger("Testing uniform_pef_sequence::find()");
    j = 0;
    for (uint64_t i = 0; i < n; ++i) {
        auto ptr_range = pointer_ranges[j];
        if (i >= ptr_range.end) {
            ++j;
            ptr_range = pointer_ranges[j];
        }
        uint64_t pos = 0;
        seq.find(ptr_range, values[i], &pos);
        util::check(i, pos, i, "position");
    }
    assert(j == pointer_ranges.size() - 1);
    essentials::logger("OK");

    return 0;
}
