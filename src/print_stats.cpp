#include <iostream>

#include "lm_types.hpp"
#include "utils/util.hpp"
#include "utils/stats.cpp"
#include "../external/essentials/include/essentials.hpp"
#include "../external/cmd_line_parser/include/parser.hpp"

using namespace tongrams;

template <typename T>
void print_stats(std::string const& index_filename) {
    T model;
    essentials::logger("Loading data structure");
    size_t bytes = util::load(model, index_filename);
    model.print_stats(bytes);
}

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Index filename.");
    if (!parser.parse()) return 1;

    auto index_filename = parser.get<std::string>("index_filename");
    auto model_string_type = util::get_model_type(index_filename);

    if (false) {
#define LOOP_BODY(R, DATA, T)                              \
    }                                                      \
    else if (model_string_type == BOOST_PP_STRINGIZE(T)) { \
        print_stats<T>(index_filename);

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, SXLM_TYPES);
#undef LOOP_BODY
    } else {
        building_util::unknown_type(model_string_type);
    }

    return 0;
}
