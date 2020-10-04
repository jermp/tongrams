#include <iostream>

#include "utils/util.hpp"
#include "lm_types.hpp"
#include "../external/essentials/include/essentials.hpp"
#include "../external/cmd_line_parser/include/parser.hpp"

using namespace tongrams;

template <typename Model>
void check_model(Model& model, std::string const& input_folder) {
    identity_adaptor adaptor;
    for (uint8_t order = 1; order <= model.order(); ++order) {
        std::string str_order = std::to_string(order);
        tongrams::grams_gzparser grams_parser(
            (input_folder + "/" + str_order + "-grams.sorted.gz").c_str());
        essentials::logger("Checking " + str_order + "-grams");
        uint64_t i = 0;
        for (auto const& l : grams_parser) {
            uint64_t count = model.lookup(l.gram, adaptor);
            if (count == global::not_found) {
                std::cout << int(order) << "-gram = '"
                          << std::string(l.gram.first, l.gram.second)
                          << "' not found at line " << i << std::endl;
            } else {
                util::check(i, count, l.count, "value");
            }
            ++i;
        }
        essentials::logger("OK");
    }
}

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("binary_filename", "Binary filename.");
    parser.add("input_folder", "Input folder.");
    if (!parser.parse()) return 1;

    auto binary_filename = parser.get<std::string>("binary_filename");
    auto input_folder = parser.get<std::string>("input_folder");
    auto model_string_type = util::get_model_type(binary_filename.c_str());

    if (false) {
#define LOOP_BODY(R, DATA, T)                                                 \
    }                                                                         \
    else if (model_string_type == BOOST_PP_STRINGIZE(T)) {                    \
        T model;                                                              \
        essentials::logger("Loading data structure");                         \
        size_t file_size = util::load(model, binary_filename.c_str());        \
        std::cout << "\tTotal bytes: " << file_size << "\n";                  \
        std::cout << "\tTotal ngrams: " << model.size() << "\n";              \
        std::cout << "\tBytes per gram: " << double(file_size) / model.size() \
                  << std::endl;                                               \
        check_model<T>(model, input_folder);

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, TONGRAMS_COUNT_TYPES);
#undef LOOP_BODY
    } else {
        std::cerr << "Error: check not supported with type "
                  << "'" << model_string_type << "'." << std::endl;
    }

    return 0;
}
