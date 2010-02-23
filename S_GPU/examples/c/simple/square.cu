// example1.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
//#include <cuda.h>
#include "s_gpu.h"
#include <mpi.h>


typedef struct param_t{
    float *a;
    int N;
} param_t;

// Kernel that executes on the CUDA device
__global__ void square_array(float *a, int N)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx<N) a[idx] = a[idx] * a[idx];
}

void callback_square(void *_param) {
  param_t *param = ((param_t*)(_param));
  float *a_d = param->a;
  int N = param->N;
  // Do calculation on device:
  int block_size = 4;
  int n_blocks = N/block_size + (N%block_size == 0 ? 0:1);

  square_array <<< n_blocks, block_size >>> (a_d, N);

}
  
void callback_beg(int i) {
  printf(">> beg on card %d\n", i);
}
void callback_end(int i) {
  printf("<< end on card %d\n", i);
}

// main routine that executes on the host
int main(int narg, char *args[]){
  int gpushare, usegpu;
  sg_init(&gpushare, &usegpu, 0);
  tracer_callbacks_t tracs;
	memset ( &tracs, NULL, sizeof ( tracer_callbacks_t ) );
  tracs.calc_beg = &callback_beg;
  tracs.calc_end = &callback_end;
  add_tracer_callbacks(&tracs);

  MPI_Init(&narg,&args);

  float *a_host, *a_pinned, *a_device;  // Pointer to host & device arrays
  const int N = 10;  // Number of elements in arrays
  size_t size = N * sizeof(float);



  a_host = (float *)malloc(size);        // Allocate array on host
  // Initialize host array and copy it to CUDA device
  for (int i=0; i<N; i++) a_host[i] = (float)i;
  
  
  sg_cpu_malloc_pinned((void**)&a_pinned, size);
  sg_gpu_malloc((void**)&a_device, size);

  sg_stream_ptr_t stream_ptr = sg_create_stream();
// ca ne peut pas marcher, car ensuite dans send_mem, les @ passées ne sont pas encore allouées (le malloc n'a pas encore eu lieu)
// pour résoudre ce problème, il faudrait que la fonction send_mem (et recv_mem) prennent en paramètre des ** et non des *.
// bref, ça change pas mal les choses.
//  sg_malloc_pinned((void**)&a_pinned, size, stream_ptr);
//  sg_malloc_gpu((void**)&a_device, size, stream_ptr);
  
  sg_gpu_send_mem(a_device, a_host, a_pinned, size, stream_ptr);

  param_t param;
  param.a = a_device;
  param.N = N;
  sg_gpu_send_calc(&callback_square,&param,sizeof(param_t),stream_ptr);

  // Retrieve result from device and store it in host array
  sg_gpu_recv_mem(a_host, a_device, a_pinned, size, stream_ptr);

  sg_exec_all_streams();

  // Print results
  for (int i=0; i<N; i++) printf("%d %f\n", i, a_host[i]);
  
  // Cleanup
  free(a_host); cudaFree(a_device); cudaFree(a_pinned);




  MPI_Finalize();
}

