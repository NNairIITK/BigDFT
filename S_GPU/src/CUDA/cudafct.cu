/**
 * @file cudafct.cu
 * @author Matthieu Ospici
 * 
 * @brief
 * cudafct is wrapping cuda function calls. This wrapping adds error checking when
 * calling these cuda functions.
 * 
 * @section LICENSE
 * 
 * Copyright (C) 2010 BULL LIG CEA-INAC UJF
 *
 * This file is part of S_GPU library.
 * 
 * S_GPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * S_GPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with S_GPU.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include <stdio.h>
#include <cuda.h>

#define CUERR(fct_name) { cudaError_t err; \
 if ((err = cudaGetLastError()) != cudaSuccess) { \
 printf("CUDA error: %s, in file %s:%d\n", cudaGetErrorString(err), __FILE__, __LINE__); } \
 printf ( "**** ERROR *** : %s\n", fct_name );}


extern "C"
int c_cudaGetDeviceCount ( int *s_gpuCount )
{
	return cudaGetDeviceCount ( s_gpuCount );
}

extern "C"
int c_cudaSetDevice ( int device )
{
	return cudaSetDevice ( device );
}

extern "C"
int c_cudaMalloc ( void** p, unsigned int memsize )
{
	int ret = cudaMalloc ( p, memsize );
	if ( ret != cudaSuccess )
	{
		CUERR("c_cuda_malloc");
	}
	return ret;
}

extern "C"
int c_cudaMallocHost ( void** p, unsigned int memsize )
{
	int ret = cudaMallocHost ( p, memsize );
	if ( ret != cudaSuccess )
	{
		CUERR("c_cudaMallocHost");
	}
	return ret;
}

//run with fstrict-aliasing
typedef union conv_u
{
	float **p;
	void **v;
} conv_u_t;

extern "C"
int c_cudaMallocHostFloat ( float **p, unsigned int memsize )
{
	conv_u_t tmp;
	tmp.p = p;
	return cudaMallocHost ( tmp.v, memsize );
}

extern "C"
int c_cuda_get_device ( int *dev )
{
	return cudaGetDevice ( dev );
}

extern "C"
int c_cuda_gpu_send_pi ( void *dest, const void *src,  size_t memByte )
{
	// std::cout << "cudaSucess : " << cudaSuccess
	//	    << "cudaErrorInvalidValue : " << cudaErrorInvalidValue
	//	    << "cudaErrorInvalidDevicePointer : " << cudaErrorInvalidDevicePointer
	//	    << "cudaErrorInvalidMemcpyDirection : " << cudaErrorInvalidMemcpyDirection<< std::endl;
	int ret =  cudaMemcpyAsync ( dest, src, memByte, cudaMemcpyHostToDevice, 0 );
	if ( ret != cudaSuccess )
	{
		CUERR("c_cuda_gpu_send_pi");
	}
  // the asynchronous copy is synchronized here
	if ( cudaStreamSynchronize ( 0 ) != cudaSuccess )
	{
		CUERR("c_cuda_gpu_send_pi stream");
	}
	return ret;
}

extern "C"
int c_cuda_gpu_recv_pi ( void *dest, const void *src,  size_t memByte )
{
	int ret = cudaMemcpyAsync ( dest, src, memByte, cudaMemcpyDeviceToHost, 0 );
	if ( ret != cudaSuccess )
	{
		CUERR("c_cuda_gpu_recv_pi");
	}
	if ( cudaStreamSynchronize ( 0 ) != cudaSuccess )
	{
		CUERR("c_cuda_gpu_recv_pi stream");
	}
	return ret;
}

extern "C"
int c_cuda_setdevice ( int device )
{
	return cudaSetDevice ( device );
}

/*int c_cuda_setdevice_ctx(int device)
{

  CUcontext pCtx;

  CUdevice dev;
  if( cuDeviceGet(&dev, device) != CUDA_SUCCESS)
    {
      CUERR;
      std::cout << "**** ERROR *** : cudeviceget" << std::endl;
    }

  if(cuCtxCreate(&pCtx, 0, dev)  != CUDA_SUCCESS)
    {
      CUERR;
      std::cout << "**** ERROR *** : cuctxcreate" << std::endl;
    }



  return 0;
}

int c_cuda_ctxpopcur(void *contex)
{
  if(cuCtxPopCurrent((CUcontext*)contex) != CUDA_SUCCESS)
    {
      CUERR;
      std::cout << "**** ERROR *** : cuctxcreate" << std::endl;
    }

  return 0;
  }

int c_cuda_ctxpushcur(void *contex)
{
  if(cuCtxPushCurrent(*(CUcontext*)contex) !=CUDA_SUCCESS)
    {
      CUERR;
      std::cout << "**** ERROR *** : cuctxcreate" << std::endl;
    }
 return 0;
 }*/
