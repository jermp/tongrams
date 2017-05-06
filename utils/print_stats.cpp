#include <iostream>

#include "../lm_types.hpp"
#include "util.hpp"
#include "stats.cpp"

using namespace tongrams;

template<typename T>
void print_stats(const char* binary_filename) {
    T model;
    util::logger("Loading data structure");
    size_t bytes = util::load(model, binary_filename);
    model.print_stats(bytes);
}

int main(int argc, char** argv)
{
    if (argc < 2 || building_util::request_help(argc, argv)) {
        building_util::display_legend();
        std::cerr << "Usage " << argv[0] << ":\n"
                  << "\t" << style::bold << style::underline << "binary_filename" << style::off
                  << std::endl;
        return 1;
    }

    const char* binary_filename = argv[1];
    auto model_string_type = util::get_model_type(binary_filename);

    if (false) {
#define LOOP_BODY(R, DATA, T)                                   \
    } else if (model_string_type == BOOST_PP_STRINGIZE(T)) {    \
        print_stats<T>(binary_filename);                        \

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, SXLM_TYPES);
#undef LOOP_BODY
    } else {
        building_util::unknown_type(model_string_type);
    }

    return 0;
}
