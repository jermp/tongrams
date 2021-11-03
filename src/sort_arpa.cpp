#include <unistd.h>

#include "sorters/sorter.hpp"
#include "sorters/sorter_common.hpp"
#include "lm_types.hpp"
#include "utils/parsers.hpp"
#include "../external/essentials/include/essentials.hpp"
#include "../external/cmd_line_parser/include/parser.hpp"

using namespace tongrams;

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("order", "n-gram order. Must be larger than 0.");
    parser.add("arpa_filename", "ARPA filename.");
    parser.add("vocab_filename", "Vocabulary filename.");
    parser.add("output_filename", "Output filename.");
    parser.add("tmp_dir", "Temporary directory for sorting.", "--tmp", false);
    parser.add("ram", "Percentage of RAM to use. It must be in (0,100].",
               "--ram", false);
    if (!parser.parse()) return 1;

    auto order = parser.get<uint32_t>("order");
    auto arpa_filename = parser.get<std::string>("arpa_filename");
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

    {
        std::vector<uint64_t> counts;
        arpa_parser ap(arpa_filename.c_str());
        ap.read_header(counts);
        if (order == 0 or order > counts.size()) {
            std::cerr << "invalid specified order" << std::endl;
            return 1;
        }

        single_valued_mpht64 vocab;
        essentials::logger("Building vocabulary");
        build_vocabulary(vocab_filename.c_str(), vocab, available_ram * 0.8);

        ap.read_line();

        // skip to specified order
        assert(order > 0);
        for (uint32_t i = 0; i != order - 1; ++i) {
            while (ap.read_line())
                ;
            ap.read_line();
            ap.read_line();
        }

        auto n = counts[order - 1];
        grams_probs_pool pool(n, ram_percentage);

        // NOTE: SUFFIX order
        typedef suffix_order_comparator(single_valued_mpht64,
                                        prob_backoff_record) comparator_type;
        comparator_type cmp(vocab);
        sorter<comparator_type, prob_backoff_line_handler> sorter(
            n, cmp, output_filename, tmp_dir);

        {  // write ARPA header
            std::ofstream os(output_filename);
            std::string header("\\" + std::to_string(order) + "-grams:\n");
            os.write(header.data(), header.size() * sizeof(char));
            os.close();
        }

        for (uint64_t i = 0; i < n - 1;) {
            ap.parse_line(pool);
            while (i != n - 1) {
                if (!ap.read_line()) {
                    std::cerr << "skipping empty line at: " << i << "/" << n - 1
                              << " "
                              << "(arpa file is malformed)" << std::endl;
                } else {
                    if (ap.parse_line(pool)) {
                        ++i;
                    } else {
                        break;
                    }
                }
            }

            auto& grams_index = pool.index();
            sorter.sort(grams_index.begin(), grams_index.end());
            pool.clear();
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
