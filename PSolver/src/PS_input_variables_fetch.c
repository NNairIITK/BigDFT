#include <config.h>

#define _GNU_SOURCE

#include <fcntl.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char ps_input_vars_def[] = 
#include "PS_input_variables_definition-inc.h"
  ;

// function to recuperate the size of the definition of the input variables
void FC_FUNC(getpsinputdefsize, GETPSINPUTDEFSIZE)(int *pt_len) {
  *pt_len = strlen(ps_input_vars_def);
}

void FC_FUNC(getpsinputdef, GETPSINPUTDEF)(char *to)
{
  memcpy(to,ps_input_vars_def, sizeof(char) * strlen(ps_input_vars_def));
}

// function to recuperate the address and the size of the definition of the input variables
void FC_FUNC(getstaticpsinputdef, GETSTATICPSINPUTDEF)( long long int *pt, int *pt_len) {
  //  *pt = (void*)bigdft_input_vars_def; // it was void **
  *pt = (long long int) ps_input_vars_def;
  *pt_len = strlen(ps_input_vars_def);
}
