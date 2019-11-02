#include <unistd.h>

#include "sorters/sorter.hpp"
#include "sorters/comparators.hpp"
#include "lm_types.hpp"
#include "utils/parsers.hpp"
#include "utils/mph_tables.hpp"
#include "utils/pools.hpp"
#include "../external/essentials/include/essentials.hpp"
#include "../external/cmd_line_parser/include/parser.hpp"

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
    single_valued_mpht64::builder builder(bytes, compact_vector(ids),
                                          identity_adaptor());
    builder.build(vocab);
}

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("ngrams_filename", "Input filename to sort.");
    parser.add("vocab_filename", "Vocabulary filename.");
    parser.add("output_filename", "Output filename.");
    parser.add("tmp_dir", "Temporary directory for sorting.", "--tmp", false);
    parser.add("ram", "Percentage of RAM to use. It must be in (0,100].",
               "--ram", false);
    if (!parser.parse()) return 1;

    auto ngrams_filename = parser.get<std::string>("ngrams_filename");
    auto vocab_filename = parser.get<std::string>("vocab_filename");
    auto output_filename = parser.get<std::string>("output_filename");

    std::string default_tmp_dir("./");
    std::string tmp_dir = default_tmp_dir;
    if (parser.parsed("tmp_dir")) {
        tmp_dir = parser.get<std::string>("tmp_dir");
        essentials::create_directory(tmp_dir);
    }

    size_t available_ram = sysconf(_SC_PAGESIZE) * sysconf(_SC_PHYS_PAGES);
    size_t ram_percentage = available_ram;
    double perc = 100.0;
    if (parser.parsed("ram")) perc = parser.get<double>("ram");
    ram_percentage *= perc / 100;
    std::cout << "Sorting with " << perc << "\% of available RAM"
              << " (" << ram_percentage << "/" << available_ram << ")"
              << std::endl;

    single_valued_mpht64 vocab;
    {
        // assume unigrams fit in memory
        grams_counts_pool unigrams_pool(available_ram);
        unigrams_pool.load_from<grams_gzparser>(vocab_filename.c_str());
        essentials::logger("Building vocabulary");
        build_vocabulary(unigrams_pool, vocab);
    }

    grams_gzparser input(ngrams_filename.c_str());
    auto n = input.num_lines();
    grams_counts_pool gp(n, ram_percentage);
    auto begin = input.begin();
    auto const end = input.end();

    typedef prefix_order_comparator(single_valued_mpht64, count_record)
        comparator_type;  // NOTE: prefix order
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
