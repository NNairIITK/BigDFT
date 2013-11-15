#include <config.h>

#define _GNU_SOURCE

#include <fcntl.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char bigdft_input_vars_def[] = 
#include "input_variables_definition-inc.h"
  ;

// function to recuperate the size of the definition of the input variables
void FC_FUNC(getinputdefsize, GETINPUTDEFSIZE)(long long int *pt_len) {
  *pt_len = strlen(bigdft_input_vars_def);
}

void FC_FUNC(getinputdef, GETINPUTDEF)(char *to)
{
  memcpy(to,bigdft_input_vars_def, sizeof(char) * strlen(bigdft_input_vars_def));
}


// function to recuperate the address and the size of the definition of the input variables
void FC_FUNC(getstaticinputdef, GETSTATICINPUTDEF)( long long int *pt, long long int *pt_len) {
  //  *pt = (void*)bigdft_input_vars_def; // it was void **
  *pt = (long long int) bigdft_input_vars_def;
  *pt_len = strlen(bigdft_input_vars_def);
}
