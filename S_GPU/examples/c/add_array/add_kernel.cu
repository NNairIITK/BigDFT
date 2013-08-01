#include "type.h"

// Kernel that executes on the CUDA device
__global__ void add_arrays ( double *a, double *b, double *c, int N, int tNum )
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if ( idx < N ) c[idx] = a[idx] + b[idx];

//  int posBlock = (N/blockDim.x ) * blockIdx.x;
//  int posThread = ((N/blockDim.x)/tNum) * threadIdx.x;
//  for (int i=0;i<(N/blockDim.x)/tNum;i++)
//    {
//      c[i + posBlock + posThread ] = a[i + posBlock + posThread] * b[i + posBlock + posThread];
//    }
}

extern "C"
{
	void callback_add ( void * );
}

void callback_add ( void *_param )
{
	param_t *param = ( ( param_t* ) ( _param ) );
	double *a_d = param->a_d;
	double *b_d = param->b_d;
	double *c_d = param->c_d;
	int N = param->N;
	// Do calculation on device:

	int block_size = 64;
	int n_blocks = N / block_size + ( N % block_size == 0 ? 0 : 1 );
	add_arrays <<< n_blocks, block_size >>> ( a_d, b_d, c_d, N, 0 );

//  const int tnum = 64;
//  dim3 threads(tnum,1);
//  dim3 grid((N)/tnum, 1);
//
//  add_arrays <<< grid, threads >>> (a_d, b_d, c_d, N, tnum);
}
