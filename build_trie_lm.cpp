#include <iostream>

#include "utils/util.hpp"
#include "lm_types.hpp"

int main(int argc, char** argv)
{
    using namespace tongrams;
    if (argc < 4 || building_util::request_help(argc, argv)) {
        building_util::display_legend();
        std::cout << "Usage " << argv[0] << ":\n"
                  << "\t" << style::bold << style::underline << "data_structure_type" << style::off << "\n"
                  << "\t" << style::bold << style::underline << "order" << style::off << "\n"
                  << "\t" << style::bold << style::underline << "value_type" << style::off << "\n"
                  << "\t[--dir " << style::underline << "input_dir" << style::off << "]\n"
                  << "\t[--remapping " << style::underline << "order" << style::off << "]\n"
                  << "\t[--ranks " << style::underline << "type" << style::off << "]"
                  << std::endl;
        building_util::print_general_params();
        building_util::print_general_info();
        std::cout << style::underline << "input_dir" << style::off << " is the directory containing the N-gram counts files. "
                  << "If omitted is assumed to be the current directory."
                  << std::endl;
        std::cout << "Remapping " << style::underline << "order" << style::off << " must be an integer in [0, 3). "
                  << "If omitted is given a default value of 0."
                  << std::endl;
        std::cout << "Ranks " << style::underline << "type" << style::off << " must be either 'IC', 'PSEF' or 'PSPEF'. "
                  << "If omitted is assigned 'IC' by default.\n"
                  << "See also the file 'utils/util_types.hpp'."
                  << std::endl;
        return 1;
    }
    
    std::string data_structure_t = argv[1];
    uint64_t order = std::stoull(argv[2]);
    building_util::check_order(order);
    float unk_prob = global::default_unk_prob;
    uint8_t probs_quantization_bits = global::default_probs_quantization_bits;
    uint8_t backoffs_quantization_bits = global::default_backoffs_quantization_bits;
    uint64_t remapping_order = 0;
    char const* input_dir = "./";
    char const* arpa_filename = nullptr;
    char const* output_filename = nullptr;

    binary_header bin_header;
    bin_header.remapping_order = 0;        // default remapping order
    bin_header.ranks_t = ranks_type::IC;   // default ranks_type

    if (data_structure_t == std::string("ef_trie")) {
        bin_header.data_structure_t = data_structure_type::ef_trie;
    }
    if (data_structure_t == std::string("pef_trie")) {
        bin_header.data_structure_t = data_structure_type::pef_trie;
    }

    if (binary_header::is_invalid(bin_header.data_structure_t)) {
        std::cerr << "Error: invalid data structure type.\n"
                  << "Either 'ef_trie' or 'pef_trie' must be specified."
                  << std::endl;
        return 1;
    }

    for (int i = 3; i < argc; ++i) {
        if (argv[i] == std::string("count")) {
            bin_header.value_t = value_type::count;
        } else if (argv[i] == std::string("prob_backoff")) {
            bin_header.value_t = value_type::prob_backoff;
        } else if (argv[i] == std::string("--p")) {
            probs_quantization_bits = std::stoull(argv[++i]);
        } else if (argv[i] == std::string("--b")) {
            backoffs_quantization_bits = std::stoull(argv[++i]);
        } else if (argv[i] == std::string("--u")) {
            unk_prob = std::stof(argv[++i]);
            building_util::check_unk_logprob(unk_prob);
            std::cout << "sustituting <unk> probability with: "
                      << unk_prob << std::endl;
        } else if (argv[i] == std::string("--remapping")) {
            remapping_order = std::stoull(argv[++i]);
            building_util::check_remapping_order(remapping_order);
            bin_header.remapping_order = remapping_order;
        } else if (argv[i] == std::string("--ranks")) {
            ++i;            
            if (argv[i] == std::string("IC")) {
                bin_header.ranks_t = ranks_type::IC;
            }
            if (argv[i] == std::string("PSEF")) {
                bin_header.ranks_t = ranks_type::PSEF;
            }
            if (argv[i] == std::string("PSPEF")) {
                bin_header.ranks_t = ranks_type::PSPEF;
            }
        } else if (argv[i] == std::string("--arpa")) {
            arpa_filename = argv[++i];
        } else if (argv[i] == std::string("--dir")) {
            input_dir = argv[++i];
        } else if (argv[i] == std::string("--out")) {
            output_filename = argv[++i];
        } else {
            std::cerr << "Unknown option: '"
                      << argv[i] << "'" << std::endl;
            return 1;
        }
    }

    if (binary_header::is_invalid(bin_header.value_t)) {
        std::cerr << "Error: invalid data type.\n"
                  << "Either 'count' or 'prob' must be specified."
                  << std::endl;
        return 1;
    }

    if (bin_header.value_t == value_type::count && arpa_filename != nullptr) {
        std::cerr << "warning: option '--arpa' ignored with data type 'count' specified."
                  << std::endl << std::endl;    
    }

    uint8_t header = bin_header.get();
    auto model_string_type = bin_header.parse(header);

    auto out_filename = model_string_type;
    if (output_filename == nullptr) { // assign default output filename
        out_filename += std::string(".out");
        output_filename = out_filename.c_str();
    }

    if (bin_header.value_t == value_type::count) {
        if (false) {
    #define LOOP_BODY(R, DATA, T)                                       \
        } else if (model_string_type == BOOST_PP_STRINGIZE(T)) {        \
                                                                        \
            T::builder builder(input_dir, order, remapping_order);      \
            T model;                                                    \
            builder.build(model);                                       \
            util::save(header, model, output_filename);                 \

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, SXLM_TRIE_COUNT_TYPES);
    #undef LOOP_BODY
        } else {
            building_util::unknown_type(model_string_type);
        }
    } else { // bin_header.value_t == value_type::prob_backoff
        if (false) {
    #define LOOP_BODY(R, DATA, T)                                   \
        } else if (model_string_type == BOOST_PP_STRINGIZE(T)) {    \
                                                                    \
            T::builder builder(arpa_filename,                       \
                               order,                               \
                               remapping_order,                     \
                               unk_prob,                            \
                               probs_quantization_bits,             \
                               backoffs_quantization_bits);         \
            T model;                                                \
            builder.build(model);                                   \
            util::save(header, model, output_filename);             \

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, SXLM_TRIE_PROB_TYPES);
    #undef LOOP_BODY
        } else {
            building_util::unknown_type(model_string_type);
        }
    }

    return 0;
}
