#include <random>

#include "../utils/util.hpp"
#include "../vectors/compact_vector.hpp"

int main(int argc, char** argv)
{
    using namespace tongrams;
    if (argc < 3 || building_util::request_help(argc, argv)) {
        building_util::display_legend();
        std::cout << "Usage " << argv[0] << ":\n"
                  << "\t" << style::bold << style::underline << "num_of_values" << style::off << "\n"
                  << "\t" << style::bold << style::underline << "bits_per_value" << style::off
                  << std::endl;
        return 1;
    }

    uint64_t n = util::toull(argv[1]);
    if (!n) {
        std::cerr << "Error: number of values must be non-zero."
                  << std::endl;
        std::terminate();
    }

    uint64_t w = util::toull(argv[2]);

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

    util::logger("Checking random access");
    uint32_t i = 0;
    for (; i < n; ++i) {
        util::check(i, values[i], v[i], "value");
    }

    util::logger("Checking size");
    util::check(0, values.size(), n, "size");

    util::logger("Checking iterator");
    auto it = v.begin();
    i = 0;
    for (auto val: values) {
        util::check(i++, val, *it++, "value");
    }

    util::logger("Checking sequential filler");

    double tick = util::get_time_usecs();
    for (uint32_t j = 0; j < 10; ++j)
    {
        compact_vector::builder cvb(v.begin(), n, w);
        values.build(cvb);
        i = 0;
        it = v.begin();
        for (auto val: values) {
            util::check(i++, val, *it++, "value");
        }
    }
    double elapsed = (util::get_time_usecs() - tick) / 1000000;
    std::cout << "\ttook " << elapsed << " [secs] for 10 iterations" << std::endl;

    util::logger("Writing to disk");
    util::save(global::null_header, values, "./tmp.out");

    util::logger("Loading from disk");
    compact_vector loaded_values;
    size_t file_size = util::load(loaded_values, "./tmp.out");
    util::logger("read " + std::to_string(file_size) + " bytes");

    util::logger("Checking loaded values");
    for (i = 0; i < loaded_values.size(); ++i) {
        util::check(i, loaded_values[i], v[i], "value");
    }
    util::logger("OK");

    util::logger("Checking loaded iterator values");
    it = v.begin();
    i = 0;
    for (auto val: loaded_values) {
        util::check(i++, val, *it++, "value");
    }
    util::logger("OK");
    std::remove("./tmp.out");

    return 0;
}
