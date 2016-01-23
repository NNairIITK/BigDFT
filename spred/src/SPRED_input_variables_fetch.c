#include <config.h>

#define _GNU_SOURCE

#include <fcntl.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char spred_input_vars_def[] = 
#include "PS_input_variables_definition-inc.h"
  ;

// function to recuperate the size of the definition of the input variables
void FC_FUNC(getspredinputdefsize, GETPSINPUTDEFSIZE)(int *pt_len) {
  *pt_len = strlen(spred_input_vars_def);
}

void FC_FUNC(getspredinputdef, GETPSINPUTDEF)(char *to)
{
  memcpy(to,spred_input_vars_def, sizeof(char) * strlen(spred_input_vars_def));
}
