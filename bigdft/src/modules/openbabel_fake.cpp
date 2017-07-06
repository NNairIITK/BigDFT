#include <config.h>
#include <stdlib.h>
extern "C" {
#include <futile.h>
}

extern "C" void FC_FUNC_(openbabel_load, OPENBABEL_LOAD)(f90_dictionary_pointer *dict_posinp,
                                              const char *filename, unsigned int *flen)
{
}

extern "C" void FC_FUNC_(openbabel_dump, OPENBABEL_DUMP)(f90_dictionary_pointer *dict_posinp,
                                              f90_dictionary_pointer *dict_types,
                                              const char *filename, unsigned int *flen)
{
}
