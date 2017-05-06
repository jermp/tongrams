#include <random>
#include <cmath>

#include "../utils/util.hpp"
#include "../vectors/hash_compact_vector.hpp"

using namespace tongrams;

template<typename HashType>
void perf_test(uint64_t n, uint64_t w)
{
    typedef HashType hash_t;
    std::vector<hash_t> keys;
    std::vector<uint64_t> values;
    keys.reserve(n);
    values.reserve(n);

    typename hash_compact_vector<hash_t>::builder hcvb(n, w);
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_int_distribution<uint64_t>
        distr(0, std::pow(2, w - 1));

    util::logger("Generating random (key, value) pairs");
    for (uint64_t i = 0; i < n; ++i) {
        hash_t k = distr(eng);
        uint64_t v = distr(eng);
        keys.push_back(k);
        values.push_back(v);
        hcvb.set(i, k, v);
    }

    hash_compact_vector<hash_t> cv(hcvb);
    util::logger("Checking hash_compact_vector random access");
    for (uint64_t i = 0; i < n; ++i) {
        auto got = cv[i];
        util::check(i, got.first, keys[i], "key");
        util::check(i, got.second, values[i], "value");
    }
    util::logger("OK");

    util::logger("Writing to disk");
    util::save(global::null_header, cv, "./tmp.out");

    util::logger("Loading from disk");
    hash_compact_vector<hash_t> loaded;
    size_t file_size = util::load(loaded, "./tmp.out");
    util::logger("read " + std::to_string(file_size) + " bytes");

    util::logger("Checking loaded values");
    for (uint64_t i = 0; i < n; ++i) {
        auto got = loaded[i];
        util::check(i, got.first, keys[i], "key");
        util::check(i, got.second, values[i], "value");
    }
    util::logger("OK");
    std::remove("./tmp.out");
}

int main(int argc, char** argv)
{
    if (argc < 4 || building_util::request_help(argc, argv)) {
        building_util::display_legend();
        std::cout << "Usage " << argv[0] << ":\n"
                  << "\t" << style::bold << style::underline << "num_of_values" << style::off << "\n"
                  << "\t" << style::bold << style::underline << "bits_per_hash" << style::off << "\n"
                  << "\t" << style::bold << style::underline << "bits_per_value" << style::off
                  << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        std::cout << style::bold << style::underline << "bits_per_hash" << style::off << " must be 32 or 64."
                  << std::endl;
        return 1;
    }

    uint64_t n = util::toull(argv[1]);
    if (!n) {
        std::cerr << "Error: number of values must be non-zero."
                  << std::endl;
        std::terminate();
    }
    
    uint64_t hash_width = util::toull(argv[2]);
    uint64_t w = util::toull(argv[3]);

    if (hash_width == 32) {
        perf_test<uint32_t>(n, w);
    } else if (hash_width == 64) {
        perf_test<uint64_t>(n, w);
    } else {
        std::cerr << "Error: hash width must be 32 or 64."
                  << std::endl;
        std::terminate();
    }

    return 0;
}
