To install the python wrapper simply run

    $ build_install.sh

To test it, first create an index. For example, run the following command from the `build/` directory

    $ ./build_trie_lm pef_trie 5 count --dir ../test_data/ --ranks PSEF --out ../data_pef_trie.count.out

And then

    $ python test_after_install.py

The python wrapper should be generalized to support all data structures (currently supports only the `pef_trie` with `PSEF` ranks).

