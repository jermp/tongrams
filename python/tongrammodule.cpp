#include <Python.h>

#include "lm_types.hpp"

#include "../utils/util.hpp"
#include "../utils/pools.hpp"

using namespace tongrams;
using namespace std;

pef_trie_PSEF_ranks_count_lm model;


static PyObject * tongram_load(PyObject *self, PyObject *args){
    const char *binary_filename;

    if (!PyArg_ParseTuple(args, "s", &binary_filename)) return NULL;

    auto model_string_type = util::get_model_type(binary_filename);

	std::cout << "tongram model_string_type: " << model_string_type << "\n";
	util::logger("Loading data structure");
	size_t file_size = util::load(model, binary_filename);

	std::cout << "\tTotal bytes: " << file_size << "\n";
	std::cout << "\tTotal ngrams: " << model.size() << "\n";
	std::cout << "\tBytes per gram: " << double(file_size) / model.size() << std::endl;

    return Py_BuildValue("i", 0);
}

static PyObject * tongram_lookup(PyObject *self, PyObject *args){
    const char *ngrams_space_separated;

    if (!PyArg_ParseTuple(args, "s", &ngrams_space_separated)) return NULL;

    stl_string_adaptor adaptor;
	uint64_t value1 = model.lookup(ngrams_space_separated, adaptor);

    return Py_BuildValue("i", value1);
}

static PyMethodDef TongramMethods[] = {

	{"load",    tongram_load  , METH_VARARGS, "Load a binary tongram model."},
	{"lookup",  tongram_lookup, METH_VARARGS, "Look up a space separated ngram."},
    
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
