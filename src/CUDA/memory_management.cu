#include <stdio.h>


//functions for allocate / deallocate / send memory to GPU


#define CUERR { cudaError_t err; \
 if ((err = cudaGetLastError()) != cudaSuccess) { \
 printf("CUDA error: %s, line %d\n", cudaGetErrorString(err), __LINE__); }}

#define GPU_SIMPLE 1
#define GPU_DOUBLE 2

unsigned int getPrecisionSize();
short gpu_precision;


int memoryGPUSum = 0;
int memoryPISum = 0;


// functions to set memory operations precision
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
void cpu_pinned_allocation__(int *nsize, //memory size
			     void **CPU_pointer,
			     int *ierr) // error code, 1 if failure

		    
{
  unsigned int mem_size = (*nsize)*getPrecisionSize();
  *ierr=0;

  if(cudaMallocHost(CPU_pointer, mem_size ) != cudaSuccess)
    {
      printf("CPU pinned allocation error \n");
      *ierr=1;
      return;
    }

  // memoryPISum += mem_size;
  // printf("CurMem PI = %i \n",memoryPISum);
  
}


extern "C" 
void cpu_pinned_deallocation__(void **CPU_pointer,
			       int *ierr) // error code, 1 if failure

		    
{
 
  *ierr=0;

  if(cudaFreeHost(*CPU_pointer) != cudaSuccess)
    {
      printf("CPU pinned dealocation error \n");
      *ierr=1;
      return;
    }

  CUERR;
}






extern "C" 
void gpu_allocate__(int *nsize, //memory size
		    void **GPU_pointer, // pointer indicating the GPU address
		    int *ierr) // error code, 1 if failure

		    
{

  unsigned int mem_size = (*nsize)*getPrecisionSize();


  //allocate memory on GPU, return error code in case of problems
  *ierr=0;
  if(cudaMalloc( GPU_pointer, mem_size) != 0)
    {
      CUERR
      printf("GPU allocation error \n");
      *ierr=1;
      return;
    }

  // memoryGPUSum += mem_size;
  //  printf("CurMem PI = %i \n",memoryGPUSum);


}

extern "C" 
void gpu_int_allocate__(int *nsize, //memory size
			void **GPU_pointer, // pointer indicating the GPU address
			int *ierr) // error code, 1 if failure

		    
{
 
  unsigned int mem_size = (*nsize)*sizeof(int);


  //allocate memory on GPU, return error code in case of problems
  *ierr=0;
  if(cudaMalloc( GPU_pointer, mem_size) != 0)
    {
      printf("GPU allocation error \n");
      *ierr=1;
      return;
    }
}


extern "C" 
void gpu_deallocate__(void **GPU_pointer, // pointer indicating the GPU address
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
		void *CPU_pointer, 
		void **GPU_pointer,
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
void gpu_int_send__(int *nsize,
		    int *CPU_pointer, 
		    void **GPU_pointer,
		    int *ierr)
{

  unsigned int mem_size = (*nsize)*sizeof(int);

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
		   void *CPU_pointer, 
		   void **GPU_pointer,
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
unsigned int getPrecisionSize()
{
  if(gpu_precision == GPU_DOUBLE)
    return sizeof(double);
  else
    return sizeof(float);
}
