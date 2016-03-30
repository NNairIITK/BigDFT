#include <config.h>

void futile_initialize()
{
  FC_FUNC_(f_lib_initialize, F_LIB_INITIALIZE)();
}
void futile_finalize()
{
  FC_FUNC_(f_lib_finalize, F_LIB_FINALIZE)();
}
