#include <random>

#include "utils/util.hpp"
#include "vectors/compact_vector.hpp"
#include "../external/essentials/include/essentials.hpp"
#include "../external/cmd_line_parser/include/parser.hpp"

int main(int argc, char** argv) {
    using namespace tongrams;
    cmd_line_parser::parser parser(argc, argv);
    parser.add("num_of_values", "Number of values.");
    parser.add("bits_per_value", "Bits per value.");
    if (!parser.parse()) return 1;

    uint64_t n = parser.get<uint64_t>("num_of_values");
    uint64_t w = parser.get<uint64_t>("bits_per_value");
    if (w > 64) {
        std::cerr << "w must be < 64" << std::endl;
        return 1;
    }
    if (n == 0 or w == 0) {
        std::cerr << "Arguments must be both non zero." << std::endl;
        return 1;
    }

    std::vector<uint64_t> v;
    v.reserve(n);

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<uint64_t> distr(0, std::pow(2, w - 1));
    compact_vector::builder cvb(n, w);

    for (uint64_t i = 0; i < n; ++i) {
        uint64_t x = distr(rng);
        v.push_back(x);
        cvb.push_back(x);
    }
    compact_vector values(cvb);

    essentials::logger("Checking random access");
    uint32_t i = 0;
    for (; i < n; ++i) {
        util::check(i, values[i], v[i], "value");
    }

    essentials::logger("Checking size");
    util::check(0, values.size(), n, "size");

    essentials::logger("Checking iterator");
    auto it = v.begin();
    i = 0;
    for (auto val : values) {
        util::check(i++, val, *it++, "value");
    }

    essentials::logger("Checking sequential filler");

    double tick = util::get_time_usecs();
    for (uint32_t j = 0; j < 10; ++j) {
        compact_vector::builder cvb(v.begin(), n, w);
        values.build(cvb);
        i = 0;
        it = v.begin();
        for (auto val : values) {
            util::check(i++, val, *it++, "value");
        }
    }
    double elapsed = (util::get_time_usecs() - tick) / 1000000;
    std::cout << "\ttook " << elapsed << " [secs] for 10 iterations"
              << std::endl;

    essentials::logger("Writing to disk");
    util::save(global::null_header, values, "./tmp.out");

    essentials::logger("Loading from disk");
    compact_vector loaded_values;
    size_t file_size = util::load(loaded_values, "./tmp.out");
    essentials::logger("read " + std::to_string(file_size) + " bytes");

    essentials::logger("Checking loaded values");
    for (i = 0; i < loaded_values.size(); ++i) {
        util::check(i, loaded_values[i], v[i], "value");
    }
    essentials::logger("OK");

    essentials::logger("Checking loaded iterator values");
    it = v.begin();
    i = 0;
    for (auto val : loaded_values) {
        util::check(i++, val, *it++, "value");
    }
    essentials::logger("OK");
    std::remove("./tmp.out");

    return 0;
}
