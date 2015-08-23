#include "hw.h"
#include <Python.h>
#include <boost/python.hpp>

std::string heyworld(std::string dude)//Implementation
{
    return "Hello " + dude + "!";
}

static PyObject* helloworld(PyObject* self)
{
    return Py_BuildValue("s", "Hello, Python extensions!!");
}

static boost::python::object heyworld_wrapper(PyObject * self, PyObject * args)
{
  std::string input;
  std::string result;

  // parse arguments
  if (!PyArg_ParseTuple(args, "s", &input)) {
    return NULL;
  }

  // build the resulting string into a Python object.
  boost::python::object ret = heyworld(input);
  // free(result);

  return ret;
}

static char helloworld_docs[] =
    "helloworld( ): Any message you want to put here!!\nheyworld( str ): greet a value!!!";

static PyMethodDef helloworld_functions[] = {
    {"helloworld", (PyCFunction)helloworld, 
     METH_NOARGS, helloworld_docs},
    { "hello", heyworld_wrapper, METH_VARARGS, "Say hello" },
    { NULL, NULL }
};

void inithelloworld(void)
{
    Py_InitModule3("helloworld", helloworld_functions,
                   "Extension module example!");
}