#ifndef FUTILE_H
#define FUTILE_H

#include "futile_cst.h"
#include "dict.h"

void futile_initialize();
void futile_finalize();

#define FUTILE_F_POINTER_SIZE 64

#define FUTILE_TRANSFER_NONE      0
#define FUTILE_TRANSFER_OWNERSHIP 1

#define FUTILE_METHOD_ARG_MAX 32

typedef enum _FutileNumeric FutileNumeric;
enum _FutileNumeric
  {
    FUTILE_INTEGER_4,
    FUTILE_REAL_8
  };

typedef union
{
  int ival;
  double dval;
} FutileNumericValue;

typedef void (*FutileMethodFortranFunc)();

typedef struct _FutileArray FutileArray;
struct _FutileArray
{
  int iarg;
  FutileNumeric type;
  size_t size;
  void *data;
};

typedef struct _FutileMethod FutileMethod;
struct _FutileMethod
{
  unsigned int n_args;
  unsigned int n_strs;
  unsigned int n_arrs;
  int isfunc;
  FutileMethodFortranFunc callback;
  void *args[FUTILE_METHOD_ARG_MAX];
  int strlens[FUTILE_METHOD_ARG_MAX];
  FutileArray arrays[FUTILE_METHOD_ARG_MAX];
};

gboolean futile_object_get_method(FutileMethod *meth,
                                  const char *obj_id, const char *meth_id);
void futile_object_method_add_arg(FutileMethod *meth, void *arg);
void futile_object_method_add_arg_str(FutileMethod *meth, char *arg, int ln);
void futile_object_method_add_arg_arr(FutileMethod *meth, void *arg,
                                      FutileNumeric type, size_t size,
                                      gboolean transfer);
void futile_object_method_execute(FutileMethod *meth);
void* futile_object_method_get_arg_arr(FutileMethod *meth, unsigned int i);
void futile_object_method_clean(FutileMethod *meth);

#endif
