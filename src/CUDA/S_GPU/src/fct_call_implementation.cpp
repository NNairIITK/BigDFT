#include <iostream>

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



void fct_call_calc_generic::malloc_and_copy(void* src,size_t size)
{
  local_param =  malloc(size);
  memcpy(local_param,src,size);
}

fct_call_calc_generic::~fct_call_calc_generic()
{
  free(local_param);
}

void fct_call_calc_generic::operator()(int tID)
{
  (*f_call)(local_param);
}







