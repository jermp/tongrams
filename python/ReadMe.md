PyTongrams
----------

PyTongrams is a Python wrapper around the C++ library Tongrams.
The wrapper is realized using the standard approach of [extending Python with C++ modules](https://docs.python.org/3/extending/extending.html), via the `Python.h` API.

### Installation

To install the wrapper just run

    bash install.sh

### Example

To use `tongrams` from Python, first create some indexes.

For example, run the following commands from within a `build` directory
where the library was compiled.

    ./build_trie ef_trie 5 count --dir ../test_data --out ef_trie.count.bin
    ./build_trie ef_trie 5 prob_backoff --u -10.0 --arpa ../test_data/arpa --out ef_trie.prob_backoff.bin

And then execute the example

    python3 example.py

### Basic usage

```python
import tongrams
model = tongrams.ProbModel('../build/ef_trie.prob_backoff.bin')
sentence = "Having a little flexibility on that issue would go a long way to putting together a final package ."
score = model.score_sentence(sentence, len(sentence))
print("Score: " + str(score))
```