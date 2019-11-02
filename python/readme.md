PyTongrams
----------

To install the python wrapper simply run

    $ bash install.sh

To use `tongrams` from Python, first create some indexes.

For example, run the following commands from within a `build` directory
where the library was compiled.

    $ ./build_trie ef_trie 5 count --dir ../test_data --out ef_trie.count.bin
    $ ./build_trie ef_trie 5 prob_backoff --u -10.0 --arpa ../test_data/arpa --out ef_trie.prob_backoff.bin

And then execute the example

    $ python3 example.py

### Basic usage

```python
import tongrams
model = tongrams.ProbModel('../build/ef_trie.prob_backoff.bin')
sentence = "Having a little flexibility on that issue would go a long way to putting together a final package ."
score = model.score_sentence(sentence, len(sentence))
print("Score: " + str(score))
```