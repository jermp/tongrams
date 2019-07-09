#include <Python.h>

#include "lm_types.hpp"

#include "../utils/util.hpp"
#include "../utils/pools.hpp"

using namespace tongrams;
using namespace std;


// https://docs.python.org/3/extending/newtypes_tutorial.html

static PyObject * tongram_load(PyObject *self, PyObject *args){
    const char *binary_filename;
    pef_trie_PSEF_ranks_count_lm * model_p;

    model_p = new pef_trie_PSEF_ranks_count_lm;
    if (!PyArg_ParseTuple(args, "s", &binary_filename)) return NULL;

    auto model_string_type = util::get_model_type(binary_filename);

    std::cout << "tongram model_string_type: " << model_string_type << "\n";
    util::logger("Loading data structure");
    size_t file_size = util::load(*model_p, binary_filename);

    std::cout << "\tTotal bytes: " << file_size << "\n";
    std::cout << "\tTotal ngrams: " << model_p->size() << "\n";
    std::cout << "\tBytes per gram: " << double(file_size) / model_p->size() << std::endl;

    return Py_BuildValue("k", model_p); // unsigned long. https://docs.python.org/3/c-api/arg.html
}

// ef_trie.prob_backoff.bin
static PyObject * tongram_load_ef_trie_prob_backoff(PyObject *self, PyObject *args){
    const char *binary_filename;
    // ef_trie_prob_lm
    //pef_trie_PSEF_ranks_count_lm * model_p;
    ef_trie_prob_lm * model_p;

    model_p = new ef_trie_prob_lm;
    if (!PyArg_ParseTuple(args, "s", &binary_filename)) return NULL;
//
    auto model_string_type = util::get_model_type(binary_filename);
//
    std::cout << "tongram model_string_type: " << model_string_type << "\n";
    util::logger("Loading data structure");
    size_t file_size = util::load(*model_p, binary_filename);
//
    std::cout << "\tTotal bytes: " << file_size << "\n";
    std::cout << "\tTotal ngrams: " << model_p->size() << "\n";
    std::cout << "\tBytes per gram: " << double(file_size) / model_p->size() << std::endl;

    return Py_BuildValue("k", model_p); // unsigned long. https://docs.python.org/3/c-api/arg.html
}


static PyObject * tongram_lookup(PyObject *self, PyObject *args){
    const char *ngrams_space_separated;
    pef_trie_PSEF_ranks_count_lm *pointer_to_model;

    unsigned long pointer_as_long;
    if (!PyArg_ParseTuple(args, "ls", &pointer_as_long, &ngrams_space_separated))  return NULL; // https://docs.python.org/3/extending/extending.html
    pointer_to_model = (pef_trie_PSEF_ranks_count_lm *)pointer_as_long;

    stl_string_adaptor adaptor;

	uint64_t value1 = pointer_to_model->lookup(ngrams_space_separated, adaptor);

    // cout << "lookup value: " << value1 << std::endl;

	if (value1 == global::not_found) {
	    return Py_BuildValue("K", 0); // unsigned long long // https://docs.python.org/3/c-api/arg.html
	} else {
	    return Py_BuildValue("K", value1); // unsigned long long
	}
}

// _ef_trie_prob_backoff
static PyObject * tongram_score_ef_trie_prob_backoff(PyObject *self, PyObject *args){
    const char *ngrams_space_separated;

    ef_trie_prob_lm  *pointer_to_model;

    unsigned long pointer_as_long;
    if (!PyArg_ParseTuple(args, "ls", &pointer_as_long, &ngrams_space_separated))  return NULL; // https://docs.python.org/3/extending/extending.html
    pointer_to_model = (ef_trie_prob_lm *)pointer_as_long;

    float log10_prob = 0.0;
    float sentence_log10_prob;
    bool is_OOV = false;
    auto state = pointer_to_model->state();

    byte_range word;

    string_text_lines corpus(ngrams_space_separated);

    state.init();

    sentence_log10_prob = 0.0;

    while (!corpus.end_of_line()) {

        corpus.next_word(word);
        pointer_to_model->score(state, word, is_OOV, log10_prob);

//         std::cout << "\"word\" : \"" << std::string(word.first, word.second) << "\", "
//                << "\"log10_prob\" : " << log10_prob << "\n";

        sentence_log10_prob += log10_prob;
    }

//    std::cout << "score value: " << sentence_log10_prob << std::endl;

    return Py_BuildValue("f", sentence_log10_prob); // unsigned long long
}

static PyMethodDef TongramMethods[] = {

	{"load",                         tongram_load                       , METH_VARARGS, "Load a binary tongram model."},
	{"load_ef_trie_prob_backoff",    tongram_load_ef_trie_prob_backoff  , METH_VARARGS, "Load a binary tongram model - _ef_trie_prob_backoff."},
	{"lookup",                       tongram_lookup                     , METH_VARARGS, "Look up a space separated ngram."},
	{"score_ef_trie_prob_backoff",  tongram_score_ef_trie_prob_backoff, METH_VARARGS, "Look up a space separated ngram - _ef_trie_prob_backoff."},
    
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC inittongram(void) {
    (void) Py_InitModule("tongram", TongramMethods);
}

int main(int argc, char *argv[]) {
    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName(argv[0]);

    /* Initialize the Python interpreter.  Required. */
    Py_Initialize();

    /* Add a static module */
    inittongram();
}
