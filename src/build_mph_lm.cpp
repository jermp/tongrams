#include <iostream>

#include "utils/util.hpp"
#include "lm_types.hpp"

int main(int argc, char** argv) {
    using namespace tongrams;
    if (argc < 5 || building_util::request_help(argc, argv)) {
        building_util::display_legend();
        std::cout << "Usage " << argv[0] << ":\n"
                  << "\t" << style::bold << style::underline << "order"
                  << style::off << "\n"
                  << "\t" << style::bold << style::underline << "hash_key_bytes"
                  << style::off << "\n"
                  << "\t" << style::bold << style::underline << "count_bytes"
                  << style::off << "\n"
                  << "\t" << style::bold << style::underline << "value_type"
                  << style::off << "\n"
                  << "\t[--dir " << style::underline << "input_dir"
                  << style::off << "]" << std::endl;
        building_util::print_general_params();
        building_util::print_general_info();
        std::cout << style::underline << "input_dir" << style::off
                  << " is the directory containing the N-gram counts files. "
                  << "If omitted is assumed to be the current directory."
                  << std::endl;
        std::cout << style::bold << style::underline << "hash_key_bytes"
                  << style::off
                  << " specifies the number of bytes used for the hash keys: "
                     "it must be 4 or 8."
                  << std::endl;
        std::cout << style::bold << style::underline << "count_bytes"
                  << style::off
                  << " specifies the number of bytes used for the distinct "
                     "counts: it must be 4 or 8."
                  << std::endl;
        return 1;
    }

    uint64_t order = std::stoull(argv[1]);
    building_util::check_order(order);
    uint64_t hash_key_bytes = std::stoull(argv[2]);
    uint64_t count_bytes = std::stoull(argv[3]);
    float unk_prob = global::default_unk_prob;
    uint8_t probs_quantization_bits = global::default_probs_quantization_bits;
    uint8_t backoffs_quantization_bits =
        global::default_backoffs_quantization_bits;
    char const* input_dir = "./";
    char const* arpa_filename = nullptr;
    char const* output_filename = nullptr;

    if (hash_key_bytes != 4 && hash_key_bytes != 8) {
        std::cerr << "Error: invalid number of bytes for hash keys.\n"
                  << "It must be 4 or 8." << std::endl;
        return 1;
    }

    if (count_bytes != 4 && count_bytes != 8) {
        std::cerr << "Error: invalid number of bytes for distinct counts.\n"
                  << "It must be 4 or 8." << std::endl;
        return 1;
    }

    binary_header bin_header;
    bin_header.data_structure_t = data_structure_type::hash;
    bin_header.hash_key_bytes = hash_key_bytes;
    bin_header.hash_count_bytes = count_bytes;

    for (int i = 4; i < argc; ++i) {
        if (argv[i] == std::string("count")) {
            bin_header.value_t = value_type::count;
        } else if (argv[i] == std::string("prob_backoff")) {
            bin_header.value_t = value_type::prob_backoff;
        } else if (argv[i] == std::string("--u")) {
            unk_prob = std::stof(argv[++i]);
            building_util::check_unk_logprob(unk_prob);
            std::cout << "substituting <unk> probability with: " << unk_prob
                      << std::endl;
        } else if (argv[i] == std::string("--p")) {
            probs_quantization_bits = std::stoull(argv[++i]);
        } else if (argv[i] == std::string("--b")) {
            backoffs_quantization_bits = std::stoull(argv[++i]);
        } else if (argv[i] == std::string("--arpa")) {
            arpa_filename = argv[++i];
        } else if (argv[i] == std::string("--dir")) {
            input_dir = argv[++i];
        } else if (argv[i] == std::string("--out")) {
            output_filename = argv[++i];
        } else {
            std::cerr << "unknown option: '" << argv[i] << "'" << std::endl;
            return 1;
        }
    }

    if (binary_header::is_invalid(bin_header.value_t)) {
        std::cerr << "Error: invalid data type.\n"
                  << "Either 'count' or 'prob_backoff' must be specified."
                  << std::endl;
        return 1;
    }

    if (bin_header.value_t == value_type::count && arpa_filename != nullptr) {
        std::cerr << "warning: option '--arpa' ignored with data type 'count' "
                     "specified."
                  << std::endl
                  << std::endl;
    }

    uint8_t header = bin_header.get();
    auto model_string_type = bin_header.parse(header);

    auto out_filename = model_string_type;
    if (output_filename == nullptr) {  // assign default output filename
        out_filename += std::string(".out");
        output_filename = out_filename.c_str();
    }

    if (bin_header.value_t == value_type::count) {
        if (false) {
#define LOOP_BODY(R, DATA, T)                              \
    }                                                      \
    else if (model_string_type == BOOST_PP_STRINGIZE(T)) { \
        T model(input_dir, order);                         \
        util::save(header, model, output_filename);

            BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, SXLM_HASH_COUNT_TYPES);
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

            BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, SXLM_HASH_PROB_TYPES);
#undef LOOP_BODY
        } else {
            building_util::unknown_type(model_string_type);
        }
    }

    return 0;
}
