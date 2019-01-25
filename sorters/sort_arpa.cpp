#include <unistd.h>

#include "sorter.hpp"
#include "comparators.hpp"
#include "../lm_types.hpp"
#include "../utils/parsers.hpp"
#include "../utils/mph_tables.hpp"
#include "../utils/pools.hpp"

using namespace tongrams;
size_t available_ram;

void build_vocabulary(const char* vocab_filename, single_valued_mpht64& vocab)
{
    // assume unigrams fit in memory
    grams_counts_pool unigrams(available_ram);
    unigrams.load_from<grams_gzparser>(vocab_filename);
    auto& unigrams_pool_index = unigrams.index();
    uint64_t n = unigrams_pool_index.size();

    std::vector<byte_range> bytes;
    bytes.reserve(n);
    for (auto const& record: unigrams_pool_index) {
        bytes.push_back(record.gram);
    }

    compact_vector::builder cvb(n, util::ceil_log2(n + 1));
    for (uint64_t id = 0; id < n; ++id) {
        cvb.push_back(id);
    }

    // NOTE: build vocabulary excluding null terminators
    // from unigrams strings so that we can lookup
    // for any substring of a n-gram
    // without allocating a std::string
    single_valued_mpht64::builder
        builder(bytes,
                compact_vector(cvb),
                identity_adaptor());
    builder.build(vocab);
}

// sort always in SUFFIX order
int main(int argc, char** argv)
{
    if (argc < 5 || building_util::request_help(argc, argv)) {
        building_util::display_legend();
        std::cerr << "Usage " << argv[0] << ":\n"
                  << "\t" << style::bold << style::underline << "order" << style::off << "\n"
                  << "\t" << style::bold << style::underline << "arpa_filename" << style::off << "\n"
                  << "\t" << style::bold << style::underline << "vocab_filename" << style::off << "\n"
                  << "\t" << style::bold << style::underline << "output_filename" << style::off << "\n"
                  << "\t[--t " << style::underline << "tmp_dir" << style::off << "]\n"
                  << "\t[--ram " << style::underline << "percentage" << style::off << "]" << std::endl;
        std::cerr << "----------------------------------------------------------------\n"
                  << style::bold << style::underline << "tmp_dir" << style::off
                  << " is the directory for temporaries.\n"
                  << "If omitted is assumed to be the current directory.\n"
                  << "RAM percentage is expressed as real in (0.0, 100.0]."
                  << std::endl;
        return 1;
    }

    uint32_t order = std::atoi(argv[1]);
    const char* arpa_filename = argv[2];
    const char* vocab_filename = argv[3];
    const char* output_filename = argv[4];
    std::string tmp_dir("./");

    available_ram = sysconf(_SC_PAGESIZE)
                  * sysconf(_SC_PHYS_PAGES);
    size_t ram_percentage = available_ram;
    double perc = 100.0;

    for (int i = 5; i < argc; ++i) {
        if (argv[i] == std::string("--ram")) {
            perc = std::stod(argv[++i]);
            if (perc <= 0.0 || perc > 100.0) {
                std::cerr << "percentage must be a vaue within (0.0, 100.0]"
                          << std::endl;
                return 1;
            }
        } else if (argv[i] == std::string("--t")) {
            tmp_dir = std::string(argv[++i]);
            building_util::create_directory(tmp_dir);
        } else {
            std::cerr << "unknown option: '" << argv[i] << "'"
                      << std::endl;
            return 1;
        }
    }

    ram_percentage *= perc / 100;
    std::cout << "Sorting with " << perc << "\% of available RAM"
              << " (" << ram_percentage << "/" << available_ram << ")"
              << std::endl;

    std::vector<uint64_t> counts;
    arpa_parser ap(arpa_filename);
    ap.read_header(counts);
    if (!order || order > counts.size()) {
        std::cerr << "invalid specified order"
                  << std::endl;
        return 1;
    }

    single_valued_mpht64 vocab;
    util::logger("Building vocabulary");
    build_vocabulary(vocab_filename, vocab);

    ap.read_line();

    // skip to specified order
    for (uint8_t i = 0; i < order - 1; ++i) {
        while (ap.read_line());
        ap.read_line();
        ap.read_line();
    }

    auto n = counts[order - 1];
    grams_probs_pool pool(n, ram_percentage);

    // ngrams are sorted in SUFFIX order
    typedef suffix_order_comparator(single_valued_mpht64,
                                    prob_backoff_record) comparator_type;
    comparator_type cmp(vocab);
    sorter<
           comparator_type,
           prob_backoff_line_handler
          > sorter(n, cmp, output_filename, tmp_dir);

    {   // write ARPA header
        std::ofstream os(output_filename);
        std::string header("\\" + std::to_string(order) + "-grams:\n");
        os.write(header.data(), header.size() * sizeof(char));
        os.close();
    }

    for (uint64_t i = 0; i < n - 1;)
    {
        ap.parse_line(pool);
        while (i != n - 1) {
            if (!ap.read_line()) {
                std::cerr << "skipping empty line at: " << i << "/"
                          << n - 1 << " "
                          << "(arpa file is malformed)"
                          << std::endl;
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

    return 0;
}
