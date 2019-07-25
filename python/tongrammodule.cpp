#include <Python.h>
#include "lm_types.hpp"
#include "../utils/util.hpp"
#include "../utils/pools.hpp"
#include <sstream>

using namespace tongrams;
using namespace std;

template<typename Model>
static PyObject * tongram_load_by_type(PyObject *self, PyObject *args, std::string model_string_type){

    const char *binary_filename;
    if (!PyArg_ParseTuple(args, "s", &binary_filename)) return NULL;

    Model * model_p = new Model;

    std::cout << "Loading data structure type: " << model_string_type << "\n";
    size_t file_size = util::load(*model_p, binary_filename);

    std::cout << "\tTotal bytes: " << file_size << "\n";
    std::cout << "\tTotal ngrams: " << model_p->size() << "\n";
    std::cout << "\tBytes per gram: " << double(file_size) / model_p->size() << std::endl;

    return Py_BuildValue("ks", model_p, model_string_type.c_str()); // unsigned long. https://docs.python.org/3/c-api/arg.html
} // EOF tongram_load_by_type


static PyObject * tongram_load(PyObject *self, PyObject *args){
    const char *binary_filename;

    if (!PyArg_ParseTuple(args, "s", &binary_filename)) return NULL;

    std::string model_string_type = util::get_model_type(binary_filename);

    if (false) {
#define LOOP_BODY(R, DATA, T)                                   \
    } else if (model_string_type == BOOST_PP_STRINGIZE(T)) {    \
        return tongram_load_by_type<T>(self, args, model_string_type);      \

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, SXLM_SCORE_TYPES);
#undef LOOP_BODY
    }

    if (false) {
#define LOOP_BODY(R, DATA, T)                                   \
    } else if (model_string_type == BOOST_PP_STRINGIZE(T)) {    \
        return tongram_load_by_type<T>(self, args, model_string_type);      \

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, SXLM_TRIE_COUNT_TYPES);
#undef LOOP_BODY
    }

    std::ostringstream oss;
    oss << "Error: tongram_load() not supported with type " << "'" << model_string_type << "' ["<< binary_filename <<"].";

    PyErr_SetString(PyExc_RuntimeError, oss.str().c_str());
    return NULL;

} // EOF tongram_load


template<typename Model>
static PyObject * tongram_lookup_by_type(PyObject *self, unsigned long pointer_as_long, const char * ngrams_space_separated){

    Model * pointer_to_model = (Model *)pointer_as_long;
    stl_string_adaptor adaptor;

	uint64_t value1 = pointer_to_model->lookup(ngrams_space_separated, adaptor);

	if (value1 == global::not_found) {
	    return Py_BuildValue("K", 0); // unsigned long long // https://docs.python.org/3/c-api/arg.html
	} else {
	    return Py_BuildValue("K", value1); // unsigned long long
	}

} // EOF tongram_lookup_by_type


static PyObject * tongram_lookup(PyObject *self, PyObject *args){
    const    char *  ngrams_space_separated;
    unsigned long    pointer_as_long;
    const    char *  model_string_type;

    if (!PyArg_ParseTuple(args, "(ls)s", &pointer_as_long, &model_string_type, &ngrams_space_separated))  return NULL; // https://docs.python.org/3/extending/extending.html

    std::string tmp_string = model_string_type;

    if (false) {
#define LOOP_BODY(R, DATA, T)                                   \
    } else if (tmp_string == BOOST_PP_STRINGIZE(T)) {    \
        return tongram_lookup_by_type<T>(self, pointer_as_long, ngrams_space_separated);      \

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, SXLM_TRIE_COUNT_TYPES);
#undef LOOP_BODY
    }

    std::ostringstream oss;
    oss << "Error: tongram_lookup() not supported with type " << "'" << model_string_type;

    PyErr_SetString(PyExc_RuntimeError, oss.str().c_str());
    return NULL;

} // EOF tongram_score

template<typename Model>
static PyObject * tongram_score_by_type(PyObject *self, unsigned long pointer_as_long, const char * ngrams_space_separated){

    Model * pointer_to_model = (Model *)pointer_as_long;

    float log10_prob          = 0.0;
    float sentence_log10_prob = 0.0;
    bool  is_OOV              = false;
    auto  state               = pointer_to_model->state();

    byte_range word;

    string_text_lines corpus(ngrams_space_separated);
    state.init();

    while (!corpus.end_of_line()) {
        corpus.next_word(word);
        pointer_to_model->score(state, word, is_OOV, log10_prob);
        std::cout << "\"Word\" : \"" << std::string(word.first, word.second) << "\", " << "\"log10_prob\" : " << log10_prob << "\n";
        sentence_log10_prob += log10_prob;
    }

    return Py_BuildValue("f", sentence_log10_prob); // unsigned long long
} // EOF tongram_score_by_type


static PyObject * tongram_score(PyObject *self, PyObject *args){
    const    char *  ngrams_space_separated;
    unsigned long    pointer_as_long;
    const    char *  model_string_type;

    if (!PyArg_ParseTuple(args, "(ls)s", &pointer_as_long, &model_string_type, &ngrams_space_separated))  return NULL; // https://docs.python.org/3/extending/extending.html

    std::string tmp_string = model_string_type;

    if (false) {
#define LOOP_BODY(R, DATA, T)                                   \
    } else if (tmp_string == BOOST_PP_STRINGIZE(T)) {    \
        return tongram_score_by_type<T>(self, pointer_as_long, ngrams_space_separated);      \

        BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, SXLM_SCORE_TYPES);
#undef LOOP_BODY
    }

    std::ostringstream oss;
    oss << "Error: tongram_score() not supported with type " << "'" << model_string_type;

    PyErr_SetString(PyExc_RuntimeError, oss.str().c_str());
    return NULL;

} // EOF tongram_score

static PyMethodDef TongramMethods[] = {

	{"load",                         tongram_load                       , METH_VARARGS, "Load a binary tongram model."},
	{"lookup",                       tongram_lookup                     , METH_VARARGS, "Look up a space separated ngram."},
	{"score",                        tongram_score                      , METH_VARARGS, "score a space separated ngram"},
    
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

