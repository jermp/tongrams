#include <Python.h>

#include "lm_types.hpp"
#include "score.hpp"
#include "utils/util.hpp"

using namespace tongrams;

typedef ef_trie_IC_ranks_count_lm count_lm_type;
typedef ef_trie_prob_lm prob_lm_type;

template <typename Model>
struct PyModel {
    PyObject_HEAD Model* model;
};

typedef PyModel<count_lm_type> PyCountLM;
typedef PyModel<prob_lm_type> PyProbLM;

static PyModuleDef tongramsmodule = {
    PyModuleDef_HEAD_INIT,
    "tongrams",
    "Python wrapper for the C++ library tongrams.",
    -1,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL};

template <typename Model>
static int load_from_binary(PyModel<Model>* self, PyObject* args,
                            PyObject* kwds) {
    const char* binary_filename = nullptr;
    if (!PyArg_ParseTuple(args, "s", &binary_filename)) return 0;
    self->model = new Model;
    std::string model_string_type = util::get_model_type(binary_filename);
    std::cout << "Loading data structure type: " << model_string_type << "\n";
    size_t file_size = util::load(*(self->model), binary_filename);
    std::cout << "\tTotal bytes: " << file_size << "\n";
    std::cout << "\tTotal ngrams: " << self->model->size() << "\n";
    std::cout << "\tBytes per gram: " << double(file_size) / self->model->size()
              << std::endl;
    return 0;
}

static int PyCountLM_init(PyCountLM* self, PyObject* args, PyObject* kwds) {
    return load_from_binary<count_lm_type>(self, args, kwds);
}

static int PyProbLM_init(PyProbLM* self, PyObject* args, PyObject* kwds) {
    return load_from_binary<prob_lm_type>(self, args, kwds);
}

static void PyCountLM_dealloc(PyCountLM* self) {
    delete self->model;
    Py_TYPE(self)->tp_free(self);
}

static void PyProbLM_dealloc(PyProbLM* self) {
    delete self->model;
    Py_TYPE(self)->tp_free(self);
}

static PyObject* PyCountLM_lookup(PyCountLM* self, PyObject* args) {
    const char* ngram = nullptr;
    if (!PyArg_ParseTuple(args, "s", &ngram)) return 0;
    stl_string_adaptor adaptor;
    uint64_t count = self->model->lookup(ngram, adaptor);
    if (count == global::not_found) {
        return Py_BuildValue("K", 0);
    } else {
        return Py_BuildValue("K", count);
    }
}

static PyObject* PyProbLM_score_corpus(PyProbLM* self, PyObject* args) {
    const char* corpus_filename = nullptr;
    if (!PyArg_ParseTuple(args, "s", &corpus_filename)) return 0;
    text_lines corpus(corpus_filename);
    auto state = self->model->state();
    float tot_log10_prob = 0.0;
    while (!corpus.end_of_file()) {
        state.init();
        float sentence_log10_prob = 0.0;
        corpus.begin_line();
        while (!corpus.end_of_line()) {
            auto word = corpus.next_word();
            bool is_OOV = false;
            float log10_prob = self->model->score(state, word, is_OOV);
            sentence_log10_prob += log10_prob;
        }
        tot_log10_prob += sentence_log10_prob;
    }
    return Py_BuildValue("fk", tot_log10_prob, corpus.num_words());
}

static PyObject* PyProbLM_score_sentence(PyProbLM* self, PyObject* args) {
    const char* sentence = nullptr;
    uint64_t size = 0;
    if (!PyArg_ParseTuple(args, "sk", &sentence, &size)) return 0;
    uint8_t const* ptr = reinterpret_cast<uint8_t const*>(sentence);
    forward_byte_range_iterator it;
    it.init({ptr, ptr + size});
    auto state = self->model->state();
    state.init();
    float sentence_log10_prob = 0.0;
    uint64_t words = 0;
    while (it.has_next()) {
        auto word = it.next();
        bool is_OOV = false;
        float log10_prob = self->model->score(state, word, is_OOV);
        sentence_log10_prob += log10_prob;
        ++words;
    }
    return Py_BuildValue("fk", sentence_log10_prob, words);
}

static PyMethodDef PyCountLM_methods[] = {
    {"lookup", (PyCFunction)PyCountLM_lookup, METH_VARARGS,
     "Lookup a space-separated n-gram and returns its count."},
    {NULL} /* Sentinel */
};

static PyMethodDef PyProbLM_methods[] = {
    {"score_sentence", (PyCFunction)PyProbLM_score_sentence, METH_VARARGS,
     "Score a space-separated sentence its total log_10(prob) score."},
    {"score_corpus", (PyCFunction)PyProbLM_score_corpus, METH_VARARGS,
     "Score a space-separated text and returns its total log_10(prob) score. "
     "Assume "
     "one sentence per line."},
    {NULL} /* Sentinel */
};

static PyTypeObject PyCountLMType = {
    PyVarObject_HEAD_INIT(NULL, 0) "tongrams.CountModel" /* tp_name */
};

static PyTypeObject PyProbLMType = {
    PyVarObject_HEAD_INIT(NULL, 0) "tongrams.ProbModel" /* tp_name */
};

PyMODINIT_FUNC PyInit_tongrams(void) {
    PyObject* m;

    PyCountLMType.tp_new = PyType_GenericNew;
    PyCountLMType.tp_basicsize = sizeof(PyCountLM);
    PyCountLMType.tp_dealloc = (destructor)PyCountLM_dealloc;
    PyCountLMType.tp_flags = Py_TPFLAGS_DEFAULT;
    PyCountLMType.tp_doc = "CountModel class";
    PyCountLMType.tp_methods = PyCountLM_methods;
    PyCountLMType.tp_init = (initproc)PyCountLM_init;

    PyProbLMType.tp_new = PyType_GenericNew;
    PyProbLMType.tp_basicsize = sizeof(PyProbLM);
    PyProbLMType.tp_dealloc = (destructor)PyProbLM_dealloc;
    PyProbLMType.tp_flags = Py_TPFLAGS_DEFAULT;
    PyProbLMType.tp_doc = "ProbModel class";
    PyProbLMType.tp_methods = PyProbLM_methods;
    PyProbLMType.tp_init = (initproc)PyProbLM_init;

    if (PyType_Ready(&PyCountLMType) < 0) return NULL;
    if (PyType_Ready(&PyProbLMType) < 0) return NULL;

    m = PyModule_Create(&tongramsmodule);
    if (m == NULL) return NULL;

    Py_INCREF(&PyCountLMType);
    PyModule_AddObject(
        m, "CountModel",
        (PyObject*)&PyCountLMType);  // Add CountModel object to the module

    Py_INCREF(&PyProbLMType);
    PyModule_AddObject(
        m, "ProbModel",
        (PyObject*)&PyProbLMType);  // Add ProbModel object to the module

    return m;
}
