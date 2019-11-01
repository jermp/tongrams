#include <unistd.h>

#include "sorters/sorter.hpp"
#include "sorters/comparators.hpp"
#include "lm_types.hpp"
#include "utils/parsers.hpp"
#include "utils/mph_tables.hpp"
#include "utils/pools.hpp"
#include "../external/essentials/include/essentials.hpp"

using namespace tongrams;

void build_vocabulary(grams_counts_pool& pool, single_valued_mpht64& vocab) {
    auto const& index = pool.index();
    auto n = index.size();
    std::vector<byte_range> bytes;
    bytes.reserve(n);
    compact_vector::builder ids(n, util::ceil_log2(n + 1));
    uint64_t id = 0;
    for (auto const& record : index) {
        bytes.push_back(record.gram);
        ids.push_back(id++);
    }

    // NOTE: build vocabulary excluding null terminators
    // from unigrams strings so that we can lookup
    // for any substring of a n-gram
    // without allocating a std::string
    single_valued_mpht64::builder builder(bytes, compact_vector(ids),
                                          identity_adaptor());
    builder.build(vocab);
}

int main(int argc, char** argv) {
    if (argc < 4 || building_util::request_help(argc, argv)) {
        building_util::display_legend();
        std::cerr << "Usage " << argv[0] << ":\n"
                  << "\t" << style::bold << style::underline
                  << "ngrams_filename" << style::off << "\n"
                  << "\t" << style::bold << style::underline << "vocab_filename"
                  << style::off << "\n"
                  << "\t" << style::bold << style::underline
                  << "output_filename" << style::off << "\n"
                  << "\t[--t " << style::underline << "tmp_dir" << style::off
                  << "]\n"
                  << "\t[--ram " << style::underline << "percentage"
                  << style::off << "]" << std::endl;
        std::cerr << "---------------------------------------------------------"
                     "-------\n"
                  << style::bold << style::underline << "tmp_dir" << style::off
                  << " is the directory for temporaries.\n"
                  << "If omitted is assumed to be the current directory.\n"
                  << "RAM percentage is expressed as real in (0.0, 100.0]."
                  << std::endl;
        return 1;
    }

    const char* ngrams_filename = argv[1];
    const char* vocab_filename = argv[2];
    const char* output_filename = argv[3];
    std::string default_tmp_dir("./");
    std::string tmp_dir = default_tmp_dir;

    size_t available_ram = sysconf(_SC_PAGESIZE) * sysconf(_SC_PHYS_PAGES);
    size_t ram_percentage = available_ram;
    double perc = 100.0;

    for (int i = 4; i < argc; ++i) {
        if (argv[i] == std::string("--ram")) {
            perc = std::stod(argv[++i]);
            if (perc <= 0.0 || perc > 100.0) {
                std::cerr << "percentage must be a vaue within (0.0, 100.0]"
                          << std::endl;
                return 1;
            }
        } else if (argv[i] == std::string("--t")) {
            tmp_dir = std::string(argv[++i]);
            essentials::create_directory(tmp_dir);
        } else {
            std::cerr << "unknown option: '" << argv[i] << "'" << std::endl;
            return 1;
        }
    }

    ram_percentage *= perc / 100;
    std::cout << "Sorting with " << perc << "\% of available RAM"
              << " (" << ram_percentage << "/" << available_ram << ")"
              << std::endl;

    single_valued_mpht64 vocab;
    {
        // assume unigrams fit in memory
        grams_counts_pool unigrams_pool(available_ram);
        unigrams_pool.load_from<grams_gzparser>(vocab_filename);
        essentials::logger("Building vocabulary");
        build_vocabulary(unigrams_pool, vocab);
    }

    grams_gzparser parser(ngrams_filename);

    auto n = parser.num_lines();
    grams_counts_pool gp(n, ram_percentage);

    auto begin = parser.begin();
    auto const end = parser.end();

    // ngrams are sorted in PREFIX order
    typedef prefix_order_comparator(single_valued_mpht64, count_record)
        comparator_type;
    comparator_type cmp(vocab);
    {
        sorter<comparator_type, count_line_handler> sorter(
            n, cmp, output_filename, tmp_dir);

        for (uint64_t i = 0; i < n - 1;) {
            auto const& l = *begin;
            gp.append(count_record(l.gram, l.count));
            ++begin;
            while (begin != end) {
                auto const& l = *begin;
                if (gp.append(count_record(l.gram, l.count))) {
                    ++i;
                    ++begin;
                } else {
                    break;
                }
            }

            auto& grams_index = gp.index();
            sorter.sort(grams_index.begin(), grams_index.end());
            gp.clear();
            ++i;
        }
    }

    if (tmp_dir != default_tmp_dir) {
        if (!essentials::remove_directory(tmp_dir)) {
            std::cerr << "directory '" << tmp_dir << "' not removed"
                      << std::endl;
        }
    }

    return 0;
}
