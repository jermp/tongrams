#pragma once

#include "utils/util_types.hpp"

namespace tongrams {

namespace version {
static const uint8_t version_number = 10;  // xy read as 'x.y'

std::string lib_version() {
    std::string v = std::to_string(version_number / 10) + "." +
                    std::to_string(version_number % 10);
    return v;
}
}  // namespace version

struct binary_header {
    /*
                            2 bits           2 bits      1 bit      2 bits
                      ---------------------------------------------------------
        trie_count    |   ranks_type |remapping_order|value_t|data_structure_t|
                      ---------------------------------------------------------
                  ------------------------------------------
        trie_prob |remapping_order|value_t|data_structure_t|
                  ------------------------------------------


                           1 bits      1 bit      2 bits
                      -----------------------------------------
        hash_count    |hash_key_bytes|value_t|data_structure_t|
                      -----------------------------------------
                    -----------------------------------------
        hash_prob   |hash_key_bytes|value_t|data_structure_t|
                    -----------------------------------------
    */

    static const int invalid = -1;

    binary_header()
        : data_structure_t(invalid)
        , value_t(invalid)
        , remapping_order(invalid)
        , hash_key_bytes(invalid)
        , ranks_t(invalid) {}

    static bool is_invalid(int param) {
        return param == invalid;
    }

    static void check_is_valid(int param) {
        if (is_invalid(param)) {
            throw std::runtime_error("Error: invalid header param");
        }
    }

    uint8_t get() {
        uint8_t header = 0;
        uint8_t position = 0;

        check_is_valid(data_structure_t);
        header |= data_structure_t;
        position += 2;

        check_is_valid(value_t);
        header |= value_t << position;
        ++position;

        if (data_structure_t == data_structure_type::hash) {
            check_is_valid(hash_key_bytes);
            header |= (hash_key_bytes / 4 - 1) << position;
        } else {
            check_is_valid(remapping_order);
            header |= remapping_order << position;
            if (value_t == value_type::count) {
                position += 2;
                check_is_valid(ranks_t);
                header |= ranks_t << position;
            }
        }

        return header;
    }

    std::string parse(uint8_t header, bool verbose = false) {
        if (verbose) {
            std::cout << "==== tongrams binary format ====\n"
                      << "library version: " << version::lib_version() << "\n";
        }

        std::string model_string_type = "";
        int data_structure_t = header & 3;
        switch (data_structure_t) {
            case data_structure_type::hash:
                model_string_type += "mph";
                break;
            case data_structure_type::ef_trie:
                model_string_type += "ef_trie";
                break;
            case data_structure_type::pef_trie:
                model_string_type += "pef_trie";
                break;
            default:
                assert(false);
        }

        header >>= 2;
        int value_t = header & 1;
        header >>= 1;

        if (verbose) {
            std::cout << "data structure type: " << model_string_type << "\n";
        }

        if (data_structure_t == data_structure_type::hash) {
            int hash_key_bytes = ((header & 1) + 1) * 4;
            model_string_type += hash_key_bytes == 4 ? "32" : "64";
            if (verbose) {
                std::cout << "hash_key_bytes: " << hash_key_bytes << "\n";
            }
        } else {
            int remapping_order = header & 3;
            header >>= 2;
            if (remapping_order) {
                if (data_structure_t == data_structure_type::ef_trie) {
                    model_string_type = "ef_rtrie";
                } else {
                    model_string_type = "pef_rtrie";
                }
            }

            if (verbose) {
                std::cout << "remapping order: " << remapping_order << "\n";
            }

            if (value_t == value_type::count) {
                int ranks_t = header & 3;
                switch (ranks_t) {
                    case ranks_type::IC:
                        model_string_type += "_IC_ranks";
                        if (verbose) {
                            std::cout << "ranks type: IC"
                                      << "\n";
                        }
                        break;
                    case ranks_type::PSEF:
                        model_string_type += "_PSEF_ranks";
                        if (verbose) {
                            std::cout << "ranks type: PSEF"
                                      << "\n";
                        }
                        break;
                    case ranks_type::PSPEF:
                        model_string_type += "_PSPEF_ranks";
                        if (verbose) {
                            std::cout << "ranks type: PSPEF"
                                      << "\n";
                        }
                        break;
                    default:
                        assert(false);
                }
            }
        }

        switch (value_t) {
            case value_type::count:
                model_string_type += "_count_lm";
                if (verbose) {
                    std::cout << "value type: count"
                              << "\n";
                }
                break;
            case value_type::prob_backoff:
                model_string_type += "_prob_lm";
                if (verbose) {
                    std::cout << "value type: prob_backoff"
                              << "\n";
                }
                break;
            default:
                assert(false);
        }

        if (verbose) {
            std::cout << "================================" << std::endl;
        }

        return model_string_type;
    }

    int data_structure_t;
    int value_t;
    int remapping_order;
    int hash_key_bytes;
    int ranks_t;
};

}  // namespace tongrams
