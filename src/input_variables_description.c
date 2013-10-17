#include <config.h>

#define _GNU_SOURCE

#include <fcntl.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char bigdft_input_vars_def[] = 
#include "input_variables_definition-inc.h"
  ;

// function to recuperate the address and the size of the definition of the input variables
void FC_FUNC(getstaticinputdef, GETSTATICINPUTDEF)(void **pt, long *pt_len) {
  *pt = (void*)bigdft_input_vars_def;
  *pt_len = strlen(bigdft_input_vars_def);
}
