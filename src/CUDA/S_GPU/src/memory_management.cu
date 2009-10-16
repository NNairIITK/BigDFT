#include <iostream>
#include <stdio.h>
#include "exceptions.h"  

//functions for allocate / deallocate / send memory to GPU


#define CUERR { cudaError_t err; \
 if ((err = cudaGetLastError()) != cudaSuccess) { \
 printf("CUDA error: %s, line %d\n", cudaGetErrorString(err), __LINE__); }}




extern "C" 
void sg_cpu_pinned_alloc__(void **CPU_pointer,
			   int *nsize, //memory size
			   int *precision,
			   int *ierr) // error code, 1 if failure
{
  unsigned int mem_size = (*nsize)*(*precision);
  *ierr=0;
  
  try
    {
      check<cuda_error>(cudaMallocHost(CPU_pointer, mem_size ) != cudaSuccess,"CPU pinned allocation",__FILE__,__LINE__);
    }
  
  catch(std::exception &e)
    {
      
      std::cerr << "*** CUDA ERROR DETECTED" << std::endl;
      std::cerr << "ERROR MESSAGE : " << e.what() << std::endl;
      *ierr = 1;
    }
 
}


extern "C" 
void sg_cpu_pinned_free__(void **CPU_pointer,
			  int *ierr) // error code, 1 if failure
{
  *ierr=0;

  try
    {
      check<cuda_error>(cudaFreeHost(*CPU_pointer) != cudaSuccess,"cudaFreeHost",__FILE__,__LINE__);
    }
  
  catch(std::exception &e)
    {
      
      std::cerr << "*** CUDA ERROR DETECTED" << std::endl;
      std::cerr << "ERROR MESSAGE : " << e.what() << std::endl;
      *ierr = 1;
    }


  CUERR;
}






extern "C" 
void sg_gpu_alloc__(void **GPU_pointer, // pointer indicating the GPU address
		    int *nsize, //memory size
		    int *precision,
		    int *ierr) // error code, 1 if failure
{

  unsigned int mem_size = (*nsize)*(*precision);

  
  *ierr=0;


 try
    {
      check<cuda_error>(cudaMalloc( GPU_pointer, mem_size) != 0,"GPU allocation",__FILE__,__LINE__);
    }
  
 catch(std::exception &e)
    {
      
      std::cerr << "*** CUDA ERROR DETECTED" << std::endl;
      std::cerr << "ERROR MESSAGE : " << e.what() << std::endl;
      *ierr = 1;
    }
}


extern "C" 
void sg_gpu_free__(void **GPU_pointer, // pointer indicating the GPU address
		   int *ierr) // error code, 1 if failure
{

  *ierr=0;

 try
    {
      check<cuda_error>(cudaFree(*GPU_pointer) != 0,"CUDA free",__FILE__,__LINE__);
    }
  
  catch(std::exception e)
    {
      
      std::cerr << "*** CUDA ERROR DETECTED" << std::endl;
      std::cerr << "ERROR MESSAGE : " << e.what() << std::endl;
      *ierr = 1;
    }
}



extern "C"
void sg_gpu_imm_send__(void **GPU_pointer,
		       void *CPU_pointer, 
		       int *nsize,
		       int *precision,
		       int *ierr)
{
  unsigned int mem_size = (*nsize)*(*precision);


  *ierr=0;
  std::cout << "SEND memsi " << mem_size << ", prec " << *precision << std::endl;


 try
    {
      check<cuda_error>(cudaMemcpy(*GPU_pointer, CPU_pointer, mem_size, cudaMemcpyHostToDevice)  != 0,"copy host to device",__FILE__,__LINE__);
    }
  
  catch(std::exception &e)
    {
      
      std::cerr << "*** CUDA ERROR DETECTED" << std::endl;
      std::cerr << "ERROR MESSAGE : " << e.what() << std::endl;
      *ierr = 1;
    }

}



extern "C" 
void sg_gpu_imm_recv__(void *CPU_pointer, 
		       void **GPU_pointer,
		       int *nsize,
		       int *precision,
		       int *ierr)
{

  unsigned int mem_size = (*nsize)*(*precision);

 
  *ierr=0;


try
    {
      check<cuda_error>(cudaMemcpy(CPU_pointer,*GPU_pointer, mem_size, cudaMemcpyDeviceToHost),"copy device to host",__FILE__,__LINE__);
    }
  
  catch(std::exception& e)
    {
      
      std::cerr << "*** CUDA ERROR DETECTED" << std::endl;
      std::cerr << "ERROR MESSAGE : " << e.what() << std::endl;
      *ierr = 1;
    }



}


