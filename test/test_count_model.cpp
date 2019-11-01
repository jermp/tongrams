#include <iostream>

#include "utils/util.hpp"
#include "lm_types.hpp"
#include "../external/essentials/include/essentials.hpp"

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
    if (argc < 3 || building_util::request_help(argc, argv)) {
        building_util::display_legend();
        std::cout << "Usage " << argv[0] << ":\n"
                  << "\t" << style::bold << style::underline
                  << "binary_filename" << style::off << "\t" << style::bold
                  << style::underline << "input_folder" << style::off
                  << std::endl;
        std::cout << "-----------------------------------------------\n";
        return 1;
    }

    const char* binary_filename = argv[1];
    std::string input_folder = argv[2];

    auto model_string_type = util::get_model_type(binary_filename);

    if (false) {
#define LOOP_BODY(R, DATA, T)                                                 \
    }                                                                         \
    else if (model_string_type == BOOST_PP_STRINGIZE(T)) {                    \
        T model;                                                              \
        essentials::logger("Loading data structure");                         \
        size_t file_size = util::load(model, binary_filename);                \
        std::cout << "\tTotal bytes: " << file_size << "\n";                  \
        std::cout << "\tTotal ngrams: " << model.size() << "\n";              \
        std::cout << "\tBytes per gram: " << double(file_size) / model.size() \
                  << std::endl;                                               \
        check_model<T>(model, input_folder);

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, SXLM_COUNT_TYPES);
#undef LOOP_BODY
    } else {
        std::cerr << "Error: check not supported with type "
                  << "'" << model_string_type << "'." << std::endl;
    }

    return 0;
}
