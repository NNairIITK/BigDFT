#ifndef __checkcuda__
#define __checkcuda__

#include "dynamicGPU/exceptions.h"

template<class E>
inline void check_cuda_error()
{
  cudaError_t error = cudaGetLastError();


  check<E>(error != cudaSuccess,cudaGetErrorString(error));
  /*  if(error != cudaSuccess) 
    {
      fprintf(stderr,"ERROR: %s: %s\n", message, cudaGetErrorString(error) );
     
      }*/
}


#endif
