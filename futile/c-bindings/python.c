#include "config.h"

#ifdef HAVE_PYTHON
#include <Python.h>
#endif

#include <futile.h>

static void  INThandler(int sig)
{
  exit(0);
}

void FC_FUNC_(f_python_initialize, F_PYTHON_INITIALIZE)()
{
#ifdef HAVE_PYTHON
  PyObject *m, *sys, *p, *path;

  if (Py_IsInitialized())
    return;

  Py_Initialize();

  sys = PyImport_AddModule("sys");
  p = PyString_FromString(PYTHON_EXECDIR);
  path = PyObject_GetAttrString(sys, "path");
  PyList_Append(path, p);
  Py_DECREF(path);
  Py_DECREF(p);

  m = PyImport_ImportModule("futile");
  if (m == NULL)
    PyErr_Print();
  else
    PyObject_SetAttrString(PyImport_AddModule("__main__"), "futile", m);
#endif
  signal(SIGINT, INThandler);
}

void FC_FUNC_(f_python_finalize, F_PYTHON_FINALIZE)()
{
#ifdef HAVE_PYTHON
  Py_Finalize();
#endif
}

static char* f2c(const char *fbuf, int ln)
{
  char *buf;

  buf = malloc(sizeof(char) * (ln + 1));
  memcpy(buf, fbuf, sizeof(char) * ln);
  buf[ln] = '\0';
  return buf;
}

void FC_FUNC_(f_python_add_object, F_PYTHON_ADD_OBJECT)(const char *obj_id, const char *varname, void *add, int ln_obj_id, int ln_varname)
{
#ifdef HAVE_PYTHON
  char *varid;
  PyObject *pymain, *futile, *fobj, *var, *args;

  pymain = PyImport_AddModule("__main__");
  futile = PyImport_AddModule("futile");

  fobj = PyObject_GetAttrString(futile, "FObject");
  
  if (fobj != NULL)
    {
      args = Py_BuildValue("s#O&", obj_id, ln_obj_id, PyLong_FromVoidPtr, add);
      var = PyObject_CallObject(fobj, args);
      Py_DECREF(args);
      Py_DECREF(fobj);

      varid = f2c(varname, ln_varname);
      PyObject_SetAttrString(pymain, varid, var);
      free(varid);
      Py_DECREF(var);
    }
#endif
}

void FC_FUNC_(f_python_execute, F_PYTHON_EXECUTE)(const char *script, int ln_script)
{
  char *data;

  data = f2c(script, ln_script);

#ifdef HAVE_PYTHON
  PyRun_SimpleString(data);
#endif

  free(data);
}
