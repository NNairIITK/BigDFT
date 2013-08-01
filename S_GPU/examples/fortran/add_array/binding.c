#include <stdio.h>

#include "type.h"
#include "s_gpu.h"

extern void callback_mul(void *);

void binding_callback_mull_(double **a, double **b, double **c, int *N, sg_stream_ptr_t *stream_ptr)
{
  param_t param;
  param.a_d = *a;
  param.b_d = *b;
  param.c_d = *c;
  param.N = *N;


  printf(" *N %i\n",param.N);
 // callback_mull(&param);
  sg_gpu_send_calc(&callback_mul,&param,sizeof(param_t),*stream_ptr);

}

