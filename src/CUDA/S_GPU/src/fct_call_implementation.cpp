#include "fct_call_implementation.h"
#include "cudafct.h"

void fct_call_trsf_CPU_GPU::operator()(int tID)
{
 
  c_cuda_gpu_send_pi(dest,src,mem_size);
 
}
void fct_call_trsf_GPU_CPU::operator()(int tID)
{

  c_cuda_gpu_recv_pi(dest,src,mem_size);
 
}

void fct_call_memcpy::operator()(int tID)
{
  memcpy(dest,src, mem_size);
}


//================= fct_call_calc_generic ======================

fct_call_calc_generic::~fct_call_calc_generic()
{
  free(param_local); 
}


void fct_call_calc_generic::initParam_local(void *dest, const void *src, size_t param_size)
{
  param_local = malloc(param_size);
  memcpy(dest,src,param_size);
}

void fct_call_calc_generic::operator()(int tID)
{
  (*f_call)(param_local);
}
