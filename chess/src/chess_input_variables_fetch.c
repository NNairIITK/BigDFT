#include <config.h>

#define _GNU_SOURCE

#include <fcntl.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char chess_input_vars_def[] = 
#include "chess_input_variables_definition-inc.h"
  ;

// function to recuperate the size of the definition of the input variables
void FC_FUNC(getchessinputdefsize, GETCHESSINPUTDEFSIZE)(int *pt_len) {
  *pt_len = strlen(chess_input_vars_def);
}

void FC_FUNC(getchessinputdef, GETCHESSINPUTDEF)(char *to)
{
  memcpy(to,chess_input_vars_def, sizeof(char) * strlen(chess_input_vars_def));
}
