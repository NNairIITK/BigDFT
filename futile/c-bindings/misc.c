#include <config.h>

#include <string.h>
#include <stdio.h>

#include "futile.h"

void FC_FUNC_(f_lib_initialize, F_LIB_INITIALIZE)(void);
void FC_FUNC_(f_lib_finalize, F_LIB_FINALIZE)(void);
void FC_FUNC_(f_object_get_method, F_OBJECT_GET_METHOD)(const char *obj_id, const char *meth_id, int *n_args, void **callback, int ln_obj_id, int ln_meth_id);

void futile_initialize()
{
  FC_FUNC_(f_lib_initialize, F_LIB_INITIALIZE)();
}
void futile_finalize()
{
  FC_FUNC_(f_lib_finalize, F_LIB_FINALIZE)();
}

gboolean futile_object_get_method(FutileMethod *meth,
                                  const char *obj_id, const char *meth_id)
{
  void *callback;

  FC_FUNC_(f_object_get_method, F_OBJECT_GET_METHOD)(obj_id, meth_id, (int*)&meth->n_args, &callback, strlen(obj_id), strlen(meth_id));
  if (!callback)
    return FALSE;
  
  meth->callback = (FutileMethodFortranFunc)callback;

  memset(meth->args, '\0', sizeof(meth->args));
  meth->n_strs = 0;
  memset(meth->strlens, '\0', sizeof(meth->strlens));
  meth->n_arrs = 0;
  memset(meth->arrays, '\0', sizeof(meth->arrays));
  return TRUE;
}

void futile_object_method_add_arg(FutileMethod *meth, void *arg)
{
  unsigned int i;

  for (i = 0; i < FUTILE_METHOD_ARG_MAX; i++)
    if (!meth->args[i])
      {
        meth->args[i] = arg;
        return;
      }
}
void futile_object_method_add_arg_str(FutileMethod *meth, char *arg, int ln)
{
  futile_object_method_add_arg(meth, arg);
  meth->strlens[meth->n_strs++] = ln;
}
void futile_object_method_add_arg_arr(FutileMethod *meth, void *arg,
                                      FutileNumeric type, size_t size)
{
  FutileArray arr = {meth->n_arrs, type, size};

  meth->arrays[meth->n_arrs++] = arr;
  futile_object_method_add_arg(meth, arg);
  futile_object_method_add_arg(meth, &meth->arrays[meth->n_arrs - 1].size);
}

void futile_object_method_execute(FutileMethod *meth)
{
  switch (meth->n_args + meth->n_arrs)
    {
    case (0):
      switch (meth->n_strs)
        {
        case (0):
          if (meth->callback) meth->callback();
          break;
        case (1):
          if (meth->callback) meth->callback(meth->strlens[0]);
          break;
        case (2):
          if (meth->callback) meth->callback(meth->strlens[0], meth->strlens[1]);
          break;
        default:
          break;
        }
      break;
    case (1):
      switch (meth->n_strs)
        {
        case (0):
          if (meth->callback) meth->callback(meth->args[0]);
          break;
        case (1):
          if (meth->callback) meth->callback(meth->args[0], meth->strlens[0]);
          break;
        case (2):
          if (meth->callback) meth->callback(meth->args[0], meth->strlens[0], meth->strlens[1]);
          break;
        default:
          break;
        }
      break;
    case (2):
      switch (meth->n_strs)
        {
        case (0):
          if (meth->callback) meth->callback(meth->args[0], meth->args[1]);
          break;
        case (1):
          if (meth->callback) meth->callback(meth->args[0], meth->args[1], meth->strlens[0]);
          break;
        case (2):
          if (meth->callback) meth->callback(meth->args[0], meth->args[1], meth->strlens[0], meth->strlens[1]);
          break;
        default:
          break;
        }
      break;
    case (3):
      switch (meth->n_strs)
        {
        case (0):
          if (meth->callback) meth->callback(meth->args[0], meth->args[1], meth->args[2]);
          break;
        case (1):
          if (meth->callback) meth->callback(meth->args[0], meth->args[1], meth->args[2], meth->strlens[0]);
          break;
        case (2):
          if (meth->callback) meth->callback(meth->args[0], meth->args[1], meth->args[2], meth->strlens[0], meth->strlens[1]);
          break;
        default:
          break;
        }
      break;
    case (4):
      switch (meth->n_strs)
        {
        case (0):
          if (meth->callback) meth->callback(meth->args[0], meth->args[1], meth->args[2], meth->args[3]);
          break;
        case (1):
          if (meth->callback) meth->callback(meth->args[0], meth->args[1], meth->args[2], meth->args[3], meth->strlens[0]);
          break;
        case (2):
          if (meth->callback) meth->callback(meth->args[0], meth->args[1], meth->args[2], meth->args[3], meth->strlens[0], meth->strlens[1]);
          break;
        default:
          break;
        }
      break;
    default:
      fprintf(stderr, "Implement %d arguments.\n", meth->n_args + meth->n_arrs);
      break;
    }
}
