#include <stdio.h>


//functions for allocate / dealocate / send memory to GPU


#define CUERR { cudaError_t err; \
 if ((err = cudaGetLastError()) != cudaSuccess) { \
 printf("CUDA error: %s, line %d\n", cudaGetErrorString(err), __LINE__); }}

#define GPU_SIMPLE 1
#define GPU_DOUBLE 2

static int getPrecisionSize();
static short gpu_precision;







extern "C"
void set_gpu_simple__()
{
  gpu_precision = GPU_SIMPLE;
}

extern "C"
void set_gpu_double__()
{
  gpu_precision = GPU_DOUBLE;
}




extern "C" 
void gpu_allocate__(int *nsize, //memory size
		    float **GPU_pointer, // pointer indicating the GPU address
		    int *ierr) // error code, 1 if failure

		    
{

  unsigned int mem_size = (*nsize)*getPrecisionSize();


  //allocate memory on GPU, return error code in case of problems
  *ierr=0;
  if(cudaMalloc( (void**) (GPU_pointer), mem_size) != 0)
    {
      printf("GPU allocation error \n");
      *ierr=1;
      return;
    }
}

extern "C" 
void gpu_deallocate__(float **GPU_pointer, // pointer indicating the GPU address
		      int *ierr) // error code, 1 if failure
{
  //deallocate memory on GPU, return error code in case of problems
  *ierr=0;
  if(cudaFree(*GPU_pointer) != 0)
    {
      CUERR
      printf("GPU deallocation error \n");
      *ierr=1;
      return;
    }
}


//Temporary send-receive operations, displacements to be added (other routines?)


extern "C"
void gpu_send__(int *nsize,
		float *CPU_pointer, 
		float **GPU_pointer,
		int *ierr)
{

  unsigned int mem_size = (*nsize)*getPrecisionSize();

  //copy V to GPU
  *ierr=0;
  if(cudaMemcpy(*GPU_pointer, CPU_pointer, mem_size, cudaMemcpyHostToDevice)  != 0)
    {
      printf("HostToDevice Memcpy error \n");
      *ierr=1;
      return;
    }

}

extern "C" 
void gpu_receive__(int *nsize,
		   float *CPU_pointer, 
		   float **GPU_pointer,
		   int *ierr)
{

  unsigned int mem_size = (*nsize)*getPrecisionSize();

 
  *ierr=0;
  if(cudaMemcpy(CPU_pointer,*GPU_pointer, mem_size, cudaMemcpyDeviceToHost)  != 0)
    {
      CUERR;
      printf("DeviceToHost Memcpy error \n");
      printf(" %i \n",mem_size);
      *ierr=1;
      return;
    }

}

//private fct
static int getPrecisionSize()
{
  if(gpu_precision == GPU_DOUBLE)
    return sizeof(double);
  else
    return sizeof(float);
}
