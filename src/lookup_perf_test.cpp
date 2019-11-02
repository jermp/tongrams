#include <iostream>

#include "lm_types.hpp"
#include "utils/util.hpp"
#include "utils/pools.hpp"
#include "../external/essentials/include/essentials.hpp"
#include "../external/cmd_line_parser/include/parser.hpp"

using namespace tongrams;

template <typename Model>
void perf_test(std::string const& index_filename,
               std::string const& query_filename, uint32_t runs) {
    strings_pool sp;
    std::vector<size_t> offsets;
    offsets.push_back(0);

    essentials::logger("Loading strings in memory for faster lookup");
    {
        emphf::file_lines lines(query_filename.c_str());
        for (auto& l : lines) {
            auto br = bytes::split_upon_check_end(l, '\t');
            sp.append(br);
            offsets.push_back(sp.bytes());
        }
    }

    size_t test_strings = offsets.size() - 1;
    identity_adaptor adaptor;

    Model model;
    essentials::logger("Loading data structure");
    size_t file_size = util::load(model, index_filename);
    std::cout << "\tTotal bytes: " << file_size << "\n";
    std::cout << "\tTotal ngrams: " << model.size() << "\n";
    std::cout << "\tBytes per gram: " << double(file_size) / model.size()
              << std::endl;

    uint8_t const* base_addr = sp.base_addr();

    essentials::logger("Performing lookups");
    std::vector<double> query_times;
    query_times.reserve(runs - 1);

    for (size_t run = 0; run < runs; ++run) {
        auto tick = util::get_time_usecs();
        for (size_t i = 0; i < test_strings; ++i) {
            byte_range br = sp.get_bytes(base_addr, offsets[i], offsets[i + 1]);
            uint64_t count = model.lookup(br, adaptor);
            essentials::do_not_optimize_away(count);
        }
        double elapsed = double(util::get_time_usecs() - tick);
        if (run) {  // first run is not timed
            query_times.push_back(elapsed);
        }
    }

    if (false) {
        for (auto t : query_times) {
            std::cout << (t / 1000) << std::endl;
        }
    } else {
        double avg =
            std::accumulate(query_times.begin(), query_times.end(), 0.0) /
            query_times.size();
        std::cout << "\tMean per run: " << avg / 1000000 << " [sec]"
                  << std::endl;
        std::cout << "\tMean per query: " << avg / test_strings << " [musec]"
                  << std::endl;
    }
}

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Index filename.");
    parser.add("query_filename", "Query filename.");
    parser.add("runs",
               "Number of runs for the benchmark. Must be greater than 1.");
    if (!parser.parse()) return 1;

    auto index_filename = parser.get<std::string>("index_filename");
    auto query_filename = parser.get<std::string>("query_filename");
    auto runs = parser.get<uint32_t>("runs");

    if (runs < 2) {
        std::cerr << "Error: number of runs must be greater than 1."
                  << std::endl;
        return 1;
    }

    auto model_string_type = util::get_model_type(index_filename);

    if (false) {
#define LOOP_BODY(R, DATA, T)                              \
    }                                                      \
    else if (model_string_type == BOOST_PP_STRINGIZE(T)) { \
        perf_test<T>(index_filename, query_filename, runs);

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, SXLM_COUNT_TYPES);
#undef LOOP_BODY
    } else {
        std::cerr << "Error: lookup() not supported with type "
                  << "'" << model_string_type << "'." << std::endl;
    }

    return 0;
}
