#include <config.h>

#include "addresses.h"
void FC_FUNC_(callable_void, CALLABLE_VOID)(FFunc_void *func)
{
  if (func)
    (*func)();
}

void FC_FUNC(getlongaddress, GETLONGADDRESS)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  //printf("\n test long address = %p %lli\n", (void*)ptr,*address);
  return;
}

