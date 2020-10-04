#include <iostream>

#include "utils/util.hpp"
#include "lm_types.hpp"
#include "../external/cmd_line_parser/include/parser.hpp"

int main(int argc, char** argv) {
    using namespace tongrams;
    cmd_line_parser::parser parser(argc, argv);

    parser.add("order", "Language model order. It must be > 0 and <= " +
                            std::to_string(global::max_order) + ".");
    parser.add("hash_key_bytes",
               "Number of bytes for hash keys: either 4 or 8.");
    parser.add("value_type",
               "Value type. It must be either 'count' or 'prob_backoff'.");
    parser.add("dir",
               "Input directory for n-gram counts. Valid if 'count' value "
               "type is specified.",
               "--dir", false);
    parser.add(
        "arpa",
        "Input ARPA filename. Valid if 'prob_backoff' value type is specified.",
        "--arpa", false);
    parser.add("p",
               "Probability quantization bits. Valid if 'prob_backoff' value "
               "type is specified.",
               "--p", false);
    parser.add("b",
               "Backoff quantization bits. Valid if 'prob_backoff' value "
               "type is specified.",
               "--b", false);
    parser.add("unk",
               "log10 probability for the unknown <unk> word. Valid if "
               "'prob_backoff' value "
               "type is specified.",
               "--u", false);
    parser.add("out", "Output filename.", "--out", false);

    if (!parser.parse()) return 1;

    auto order = parser.get<uint64_t>("order");
    building_util::check_order(order);
    auto hash_key_bytes = parser.get<uint64_t>("hash_key_bytes");
    auto value_type = parser.get<std::string>("value_type");

    if (hash_key_bytes != 4 and hash_key_bytes != 8) {
        std::cerr << "Error: invalid number of bytes for hash keys.\n"
                  << "It must be 4 or 8." << std::endl;
        return 1;
    }

    binary_header bin_header;
    bin_header.data_structure_t = data_structure_type::hash;
    bin_header.hash_key_bytes = hash_key_bytes;

    if (value_type == "count") {
        bin_header.value_t = value_type::count;
    } else if (value_type == "prob_backoff") {
        bin_header.value_t = value_type::prob_backoff;
    }
    if (binary_header::is_invalid(bin_header.value_t)) {
        std::cerr << "Error: invalid data type.\n"
                  << "Either 'count' or 'prob_backoff' must be specified."
                  << std::endl;
        return 1;
    }

    float unk_prob = global::default_unk_prob;
    uint32_t probs_quantization_bits = global::default_probs_quantization_bits;
    uint32_t backoffs_quantization_bits =
        global::default_backoffs_quantization_bits;
    const char* input_dir = ".";
    const char* arpa_filename = nullptr;
    const char* output_filename = nullptr;

    if (parser.parsed("unk")) {
        unk_prob = parser.get<float>("unk");
        building_util::check_unk_logprob(unk_prob);
        std::cout << "substituting <unk> probability with: " << unk_prob
                  << std::endl;
    }
    if (parser.parsed("p")) {
        probs_quantization_bits = parser.get<uint32_t>("p");
    }
    if (parser.parsed("b")) {
        backoffs_quantization_bits = parser.get<uint32_t>("b");
    }

    auto dir = parser.get<std::string>("dir");
    if (parser.parsed("dir")) {
        input_dir = dir.c_str();
    }

    auto arpa = parser.get<std::string>("arpa");
    if (parser.parsed("arpa")) {
        arpa_filename = arpa.c_str();
    }

    uint8_t header = bin_header.get();
    auto model_string_type = bin_header.parse(header);

    auto default_out_filename = model_string_type + std::string(".out");
    output_filename = default_out_filename.c_str();
    auto out = parser.get<std::string>("out");
    if (parser.parsed("out")) {
        output_filename = out.c_str();
    }

    if (bin_header.value_t == value_type::count and arpa_filename != nullptr) {
        std::cerr << "warning: option '--arpa' ignored with data type 'count' "
                     "specified."
                  << std::endl;
    }

    if (bin_header.value_t == value_type::count) {
        if (false) {
#define LOOP_BODY(R, DATA, T)                              \
    }                                                      \
    else if (model_string_type == BOOST_PP_STRINGIZE(T)) { \
        T model(input_dir, order);                         \
        util::save(header, model, output_filename);

            BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, TONGRAMS_HASH_COUNT_TYPES);
#undef LOOP_BODY
        } else {
            building_util::unknown_type(model_string_type);
        }
    } else {  // bin_header.value_t == value_type::prob_backoff
        if (false) {
#define LOOP_BODY(R, DATA, T)                              \
    }                                                      \
    else if (model_string_type == BOOST_PP_STRINGIZE(T)) { \
        T::builder builder(arpa_filename, order, unk_prob, \
                           probs_quantization_bits,        \
                           backoffs_quantization_bits);    \
        T model;                                           \
        builder.build(model);                              \
        util::save(header, model, output_filename);

            BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, TONGRAMS_HASH_PROB_TYPES);
#undef LOOP_BODY
        } else {
            building_util::unknown_type(model_string_type);
        }
    }

    return 0;
}
