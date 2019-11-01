`tongrams` - Tons of *N*-Grams
----------
`tongrams` is a C++ library implementing the compressed data structures described in the paper [*Efficient Data Structures for Massive N-Gram Datasets*](http://pages.di.unipi.it/pibiri/papers/SIGIR17.pdf), by Giulio Ermanno Pibiri and Rossano Venturini, published in ACM SIGIR 2017 [1]. The proposed data structures can be used to map *N*-grams to their corresponding (integer) frequency counts or to (floating point) probabilities and backoffs for backoff-interpolated [Knenser-Ney](https://en.wikipedia.org/wiki/Kneser%E2%80%93Ney_smoothing) models.

The library features a compressed trie data structure in which *N*-grams are assigned integer identifiers (IDs) and compressed with *Elias-Fano* (Subsection 3.1 of [1]) as to support efficient searches within compressed space. The *context-based remapping* of such identifiers (Subsection 3.2 of [1]) permits to encode a word following a context of fixed length _k_, i.e., its preceding _k_ words, with an integer whose value is bounded by the number of words that follow such context and _not_ by the size of the whole vocabulary (number of uni-grams).
Additionally to the trie data structure, the library allows to build models based on *minimal perfect hashing* (MPH), for constant-time retrieval (Section 4 of [1]).

When used to store frequency counts, the data structures support a `lookup()` operation that returns the number of occurrences of the specified *N*-gram. Differently, when used to store probabilities and backoffs, the data structures implement a `score()` function that, given a text as input, computes the [perplexity](https://en.wikipedia.org/wiki/Perplexity) score of the text.

This guide is meant to provide a brief overview of the library and to illustrate its functionalities through some examples.
##### Table of contents
* [Building the code](#building-the-code)
* [Input data format](#input-data-format)
* [Building the data structures](#building-the-data-structures)
* [Tests](#tests)
* [Benchmarks](#benchmarks)
* [Statistics](#statistics)
* [Authors](#authors)
* [Bibliography](#bibliography)

Building the code
-----------------
The code has been tested on Linux Ubuntu with `gcc` 5.4.1, 7.3.0 and 8.3.0; Mac OS X El Capitan with `clang` 7.3.0; Mac OS X Mojave with `clang` 10.0.0.

The following dependencies are needed for the build: [`CMake`](https://cmake.org) and [`Boost`](https://www.boost.org).

If you have cloned the repository
without `--recursive`, you will need to perform the following commands before
building:

    $ git submodule init
    $ git submodule update

To build the code on Unix systems (see file `CMakeLists.txt` for the used compilation flags), it is sufficient to do the following.

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make

You can enable parallel compilation by specifying some jobs: `make -j4`.

For best of performace, compile as follows.

    $ cmake .. -DCMAKE_BUILD_TYPE=Release  -DTONGRAMS_USE_SANITIZERS=OFF -DEMPHF_USE_POPCOUNT=ON -DTONGRAMS_USE_POPCNT=ON -DTONGRAMS_USE_PDEP=ON
    $ make

For a debug environment, compile as follows instead.

    $ cmake .. -DCMAKE_BUILD_TYPE=Debug -DTONGRAMS_USE_SANITIZERS=ON
    $ make

Unless otherwise specified, for the rest of this guide we assume that we type the terminal commands of the following examples from the created directory `build`.

Input data format
-----------------
The *N*-gram counts files follow the [Google format](http://storage.googleapis.com/books/ngrams/books/datasetsv2.html), i.e., one separate file for each distinct value of *N* (order) listing one gram per row. We enrich this format with a file header indicating the total number of *N*-grams in the file (rows):

    <total_number_of_rows>
    <gram1> <TAB> <count1>
    <gram2> <TAB> <count2>
    <gram3> <TAB> <count3>
    ...

Such *N* files must be named according to the following convention: `<order>-grams`, where `<order>` is a placeholder for the value of *N*. The files can be left unsorted if only MPH-based models have to be built, whereas these must be sorted in *prefix order* for trie-based data structures, *according to the chosen vocabulary mapping*, which should be represented by the uni-gram file (see Subsection 3.1 of [1]). Compressing the input files with standard utilities, such as `gzip`, is highly recommended.
The utility `sort_grams` can be used to sort the *N*-gram counts files in prefix order.
In conclusion, the data structures storing frequency counts are built from a directory containing the files
* `1-grams.sorted.gz`
* `2-grams.sorted.gz`
* `3-grams.sorted.gz`
* ...

formatted as explained above.

The file listing *N*-gram probabilities and backoffs is conform to, instead, the [ARPA file format](http://www.speech.sri.com/projects/srilm/manpages/ngram-format.5.html).
The *N*-grams in the ARPA file must be sorted in *suffix order* in order to build the reversed trie data structure.
The utility `sort_arpa` can be used for that purpose.

The directory `test_data` contains:
* all *N*-gram counts files (for a total of 252,550 *N*-grams), for *N* going from 1 to 5, extracted from the Agner Fog's manual *Optimizing software in C++*, sorted in prefix order and compressed with `gzip`;
* the query file `queries.random.5K` comprising 5,000 *N*-grams (1,000 for each order and drawn at random);
* the ARPA file `arpa` which lists all *N*-grams sorted in suffix order as to build backward tries efficiently;
* the `sample_text` query file (6,075 sentence for a total of 153,583 words) used for the perplexity benchmark; its companion `sample_text.LESSER` file includes just the first 10 sentences.

For the following examples, we assume to work with the sample data contained in `test_data`.

Building the data structures
----------------------------
The two executables `build_trie_lm` and `build_mph_lm` must be used to build trie-based and (minimal perfect) hash-based language models respectively. By running the executables without any arguments, or with the flag `--help`, the detailed list of expected input parameters is shown, along with their meanings. For example, by doing

    $ ./build_trie_lm --help

the following message will show up (plus additional details, here not reported for brevity):

    $ Usage ./build_trie_lm:
        data_structure_type
        order
        value_type
        [--dir input_dir]
        [--remapping order]
        [--ranks type]
        [--u value]
        [--p value]
        [--b value]
        [--arpa arpa_filename]
        [--out output_filename]

For the available `data_structure_type`(s), `value_type`(s) and `rank_type`(s) see also the file `utils/util_types.hpp`.
We now show some examples.

##### Example 1.
The command

    $ ./build_trie_lm ef_trie 5 count --dir ../test_data --out ef_trie.count.bin

builds an Elias-Fano trie
* of order 5;
* that stores frequency counts;
* from the *N*-gram counts files contained in the directory `test_data`;
* with no context-based remapping (default);
* whose counts ranks are encoded with the indexed codewords (IC) technique (default);
* that is serialized to the binary file `ef_trie.count.bin`.

##### Example 2.
The command

    $ ./build_trie_lm pef_trie 5 count --dir ../test_data --remapping 1 --ranks PSEF  --out pef_trie.count.out

builds a partitioned Elias-Fano trie
* of order 5;
* that stores frequency counts;
* from the *N*-gram counts files contained in the directory `test_data`;
* with context-based remapping of order 1;
* whose counts ranks are encoded with prefix sums (PS) + Elias-Fano (EF);
* that is serialized to the binary file `pef_trie.count.out`.

##### Example 3.
The command

    $ ./build_trie_lm ef_trie 5 prob_backoff --remapping 2 --u -20.0 --p 8 --b 8 --arpa ../test_data/arpa --out ef_trie.prob_backoff.bin

builds an Elias-Fano trie
* of order 5;
* that stores probabilities and backoffs;
* with context-based remapping of order 2;
* with `<unk>` probability of -20.0 and using 8 bits for quantizing probabilities (`--p`) and backoffs (`--b`);
* from the arpa file named `arpa`;
* that is serialized to the binary file `ef_trie.prob_backoff.bin`.

##### Example 4.
The command

    $ ./build_mph_lm 5 8 4 count --dir ../test_data --out hash.bin

builds a MPH-based model
* of order 5;
* that uses 8 bytes per hash key and 4 bytes per unique count;
* that stores frequency counts;
* from the *N*-gram counts files contained in the directory `test_data`;
* that is serialized to the binary file `hash.bin`.

Tests
-----
The `test` directory contains the unit tests of some of the fundamental building blocks used by the implemented data structures. As usual, running the executables without any arguments will show the list of their expected input parameters.
Examples:

    $ ./compact_vector_test 10000 13
    $ ./fast_ef_sequence_test 1000000 128

The directory also contains the unit test for the data structures storing frequency counts, named `check_count_model`, which validates the implementation by checking that each count stored in the data structure is the same as the one provided in the input files from which the data structure was previously built.
Example:

    $ ./check_count_model count_data_structure.bin ../test_data

where `count_data_structure.bin` is the name of the data structure binary file and `test_data` is the name of the folder containing the input *N*-gram counts files.

Benchmarks
----------
For the examples in this section, we used a desktop machine running Mac OS X Mojave, equipped with a 2.3 GHz Intel Core i5 processor (referred to as *Desktop Mac*). The code was compiled with Apple LLVM version 10.0.0 `clang` *with* all optimizations (see section [Building the code](#building-the-code)).
We additionally replicate some experiments with an Intel(R) Core(TM) i9-9900K CPU @ 3.60 GHz, under Ubuntu 19.04, 64 bits (referred to as *Server Linux*). In this case the code was compiled with `gcc` 8.3.0.

For a data structure storing frequency counts, we can test the speed of lookup queries by using the benchmark program `lookup_perf_test`.
In the following example, we show how to build and benchmark three different data structures: **EF-Trie** with no remapping, **EF-RTrie** with remapping order 1 and **PEF-RTrie** with remapping order 2 (we use the same names for the data structures as presented in [1]). Each experiment is repeated 1,000 times over the test query file `queries.random.5K`. The benchmark program `lookup_perf_test` will show mean time per run and mean time per query (along with the total number of *N*-grams, total bytes of the data structure and bytes per *N*-gram).

    $ ./build_trie_lm ef_trie 5 count --dir ../test_data --out ef_trie.bin
    $ ./lookup_perf_test ef_trie.bin ../test_data/queries.random.5K 1000

    $ ./build_trie_lm ef_trie 5 count --remapping 1 --dir ../test_data --out ef_trie.r1.bin
    $ ./lookup_perf_test ef_trie.r1.bin ../test_data/queries.random.5K 1000

    $ ./build_trie_lm pef_trie 5 count --remapping 2 --dir ../test_data --out pef_trie.r2.bin
    $ ./lookup_perf_test pef_trie.r2.bin ../test_data/queries.random.5K 1000

The results of this (micro) benchmark are summarized in the following table.

|**Data structure** |**Remapping order** | **Bytes x gram**  | **µs x query** - Desktop Mac| **µs x query** - Server Linux|
|-------------------|:------------------:|-------------------|--------------|------------|
|EF-Trie            |                   0|2.40               |      0.435   |      0.316 |
|EF-RTrie           |                   1|1.93 (**-19.7%**)  |      0.583   |      0.428 |
|PEF-RTrie          |                   2|1.75 (**-26.9%**)  |      0.595   |      0.427 |

For a data structure storing probabilities and backoffs, we can instead test the speed of scoring a text file by using the benchmark program `score`. A complete example follows.

    $ ./build_trie_lm ef_trie 5 prob_backoff --u -10.0 --p 8 --b 8 --arpa ../test_data/arpa --out ef_trie.prob_backoff.8.8.bin
    $ ./score ef_trie.prob_backoff.8.8.bin ../test_data/sample_text

The first command will build the data structure, the second one will score the text file `sample_text` contained in `test_data`. The input text file must contain one sentence per line, with words separated by spaces. During the scoring of the file, we do not wrap each sentence with markers `<s>` and `</s>`.

An examplar output could be (OOV stands for *Out Of Vocabulary*):

	perplexity including OOVs = 493720.19
	perplexity excluding OOVs = 1094.2574
	OOVs = 55868
	corpus tokens = 153583
	corpus sentences = 6075
	elapsed time: 0.037301 [sec]

Statistics
----------
The executable `print_stats` can be used to gather useful statistics regarding the space usage of the various data structure components (e.g., gram-ID and pointer sequences for tries), as well as structual properties of the indexed *N*-gram dataset (e.g., number of unique counts, min/max range lengths, average gap of gram-ID sequences, ecc.).

As an example, the following command:

    $ ./print_stats data_structure.bin

will show the statistics for the data structure serialized to the file `data_structure.bin`.

Authors
-------
* [Giulio Ermanno Pibiri](http://pages.di.unipi.it/pibiri/), <giulio.pibiri@di.unipi.it>
* [Rossano Venturini](http://pages.di.unipi.it/rossano/), <rossano.venturini@unipi.it>

Bibliography
------------
* [1] Giulio Ermanno Pibiri and Rossano Venturini, *Efficient Data Structures for Massive N-Gram Datasets*. In the Proceedings of the 40-th ACM Conference on Research and Development in Information Retrieval (SIGIR 2017).
