#include <iostream>
#include <limits>
#include <cmath>

#include "utils/util.hpp"
#include "utils/iterators.hpp"
#include "lm_types.hpp"
#include "../external/essentials/include/essentials.hpp"

using namespace tongrams;

template <typename Model>
void score_corpus(const char* binary_filename, const char* corpus_filename) {
    Model model;
    essentials::logger("Loading data structure");
    util::load(model, binary_filename);

    text_lines corpus(corpus_filename);

    uint64_t tot_OOVs = 0;
    uint64_t corpus_tokens = 0;
    uint64_t corpus_sentences = 0;
    float tot_log10_prob = 0.0;
    float tot_log10_prob_only_OOVs = 0.0;

    float sentence_log10_prob;
    bool is_OOV;
    float log10_prob = 0.0;

    auto state = model.state();
    byte_range word;

    essentials::logger("Scoring");
    auto tick = util::get_time_usecs();

    while (!corpus.end_of_file()) {  // assume one sentence per line

        state.init();

        sentence_log10_prob = 0.0;
        is_OOV = false;

        // std::cerr << "{";

        corpus.begin_line();
        while (!corpus.end_of_line()) {
            corpus.next_word(word);
            model.score(state, word, is_OOV, log10_prob);

            // std::cerr << "\"word\" : \"" << std::string(word.first,
            // word.second) << "\", "
            //        << "\"log10_prob\" : " << log10_prob << "\n";

            sentence_log10_prob += log10_prob;
            if (is_OOV) {
                tot_log10_prob_only_OOVs += log10_prob;
                is_OOV = false;
            }
        }

        // std::cerr << "\"total\" : " << sentence_log10_prob << ", "
        //           << "\"OOVs\" : " << state.OOVs << "}" << std::endl;

        tot_OOVs += state.OOVs;
        tot_log10_prob += sentence_log10_prob;
        ++corpus_sentences;
    }

    double elapsed = double(util::get_time_usecs() - tick);

    corpus_tokens = corpus.num_words();
    std::cout.precision(8);
    std::cout << "tot_log10_prob = " << tot_log10_prob << std::endl;
    std::cout << "tot_log10_prob_only_OOVs = " << tot_log10_prob_only_OOVs
              << std::endl;
    std::cout << "perplexity including OOVs = "
              << pow(10.0,
                     -(tot_log10_prob / static_cast<double>(corpus_tokens)))
              << std::endl;
    if (corpus_tokens - tot_OOVs) {
        std::cout << "perplexity excluding OOVs = "
                  << pow(10.0,
                         -((tot_log10_prob - tot_log10_prob_only_OOVs) /
                           (static_cast<double>(corpus_tokens - tot_OOVs))))
                  << std::endl;
    }
    std::cout << "OOVs = " << tot_OOVs << "\n"
              << "corpus tokens = " << corpus_tokens << "\n"
              << "corpus sentences = " << corpus_sentences << "\n"
              << "elapsed time: " << elapsed / 1000000 << " [sec]" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 3 || building_util::request_help(argc, argv)) {
        building_util::display_legend();
        std::cout << "Usage " << argv[0] << ":\n"
                  << "\t" << style::bold << style::underline
                  << "binary_filename" << style::off << "\n"
                  << "\t" << style::bold << style::underline
                  << "corpus_filename" << style::off << std::endl;
        return 1;
    }

    const char* binary_filename = argv[1];
    const char* corpus_filename = argv[2];

    auto model_string_type = util::get_model_type(binary_filename);

    if (false) {
#define LOOP_BODY(R, DATA, T)                              \
    }                                                      \
    else if (model_string_type == BOOST_PP_STRINGIZE(T)) { \
        score_corpus<T>(binary_filename, corpus_filename);

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, SXLM_SCORE_TYPES);
#undef LOOP_BODY
    } else {
        std::cerr << "Error: score() not supported with type "
                  << "'" << model_string_type << "'." << std::endl;
    }

    return 0;
}
