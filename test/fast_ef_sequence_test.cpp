#include <numeric>
#include <random>

#include "../utils/util.hpp"
#include "../sequences/fast_ef_sequence.hpp"

int main(int argc, char** argv)
{
    using namespace tongrams;
    if (argc < 3 || building_util::request_help(argc, argv)) {
        building_util::display_legend();
        std::cout << "Usage " << argv[0] << ":\n"
                  << "\t" << style::bold << style::underline << "num_of_values" << style::off << "\n"
                  << "\t" << style::bold << style::underline << "max_range_len" << style::off
                  << std::endl;
        return 1;
    }

    uint64_t n = util::toull(argv[1]);
    if (!n) {
        std::cerr << "Error: number of values must be non-zero."
                  << std::endl;
        return 1;
    }

    uint64_t max_range_len = util::toull(argv[2]);
    if (!max_range_len) {
        std::cerr << "Error: max_range_len must be non-zero."
                  << std::endl;
        return 1;
    }

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<uint64_t>
        range_distr(1, // not empty ranges
                    max_range_len);
    std::uniform_int_distribution<uint64_t>
        value_distr(1, // elements are distinct within a range
                    50);

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
    std::cout << "number of created ranges: "
              << pointer_ranges.size() << std::endl;

    std::vector<uint64_t> values;
    values.reserve(n);
    for (auto const& ptr_range: pointer_ranges) {
        uint64_t last_value = 0;
        values.push_back(last_value); // make first element of each range equal to 0
        for (uint64_t k = ptr_range.begin + 1; k < ptr_range.end; ++k) {
            values.push_back(last_value + value_distr(rng));
            last_value = values.back();
        }
    }   
    assert(values.size() == n);

    fast_ef_sequence seq;
    util::logger("Building sequence");
    seq.build(values.begin(), values.size(), pointers,
              0); // ngram order is not used

    util::logger("Testing iterator");
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
    util::logger("OK");

    util::logger("Testing fast_ef_sequence::find()");
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
    util::logger("OK");

    return 0;
}
