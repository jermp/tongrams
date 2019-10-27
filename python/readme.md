To install the python wrapper simply run

    $ bash build_install.sh

To test it, first create some indexes.

For example, run the following commands from within a `build` directory
where the library was compiler.

    $ ./build_trie_lm pef_trie 5 count --dir ../test_data/ --ranks PSEF --out pef_trie.count.bin
    $ ./build_trie_lm pef_trie 5 prob_backoff --arpa ../test_data/arpa --out pef_trie.prob_backoff.bin

And then execute the example tests.

    $ python test_trie_count.py
    $ python test_trie_prob.py
