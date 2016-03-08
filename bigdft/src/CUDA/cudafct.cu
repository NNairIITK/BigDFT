/*!
 @author
    Copyright (C) 2010-2011 BigDFT group (LG)
    This file is distributed under the terms of the
    GNU General Public License, see ~/COPYING file
    or http://www.gnu.org/copyleft/gpl.txt .
    For the list of contributors, see ~/AUTHORS 

*/

#include <stdio.h>
#include <cuda.h>

#define CUERR { cudaError_t err; \
 if ((err = cudaGetLastError()) != cudaSuccess) { \
 printf("CUDA error: %s, line %d\n", cudaGetErrorString(err), __LINE__); }}

int c_cudaGetDeviceCount(int *s_gpuCount)
{
  return cudaGetDeviceCount(s_gpuCount);
}


int c_cudaSetDevice(int device)
{
  return cudaSetDevice(device);
}


int c_cudaMalloc(void** p,unsigned int memsize)
{
  int ret = cudaMalloc(p,memsize);
    if(ret != cudaSuccess)
    {
      CUERR;
      printf("**** ERROR *** : c_cuda_malloc\n");
    }
  return ret;
 
}


int c_cudaMallocHost(void** p,unsigned int memsize)
{
  int ret = cudaMallocHost(p,memsize);
    if(ret != cudaSuccess)
    {
      CUERR;
      printf("**** ERROR *** : c_cuda_mallocHost\n");
    }
  return ret;
  
}


//run with fstrict-aliasing
typedef union conv_u
{
  float **p;
  void **v;
} conv_u_t;

int c_cudaMallocHost(float **p,unsigned int memsize)
{
  conv_u_t tmp;
  tmp.p = p;

  return cudaMallocHost(tmp.v,memsize);
}

int c_cuda_get_device(int *dev)
{
  return cudaGetDevice(dev);
}

int c_cuda_gpu_send_pi(void *dest, const void *src,  size_t memByte)
{ 
  // std::cout << "cudaSucess : " << cudaSuccess 
  //	    << "cudaErrorInvalidValue : " << cudaErrorInvalidValue
  //	    << "cudaErrorInvalidDevicePointer : " << cudaErrorInvalidDevicePointer
  //	    << "cudaErrorInvalidMemcpyDirection : " << cudaErrorInvalidMemcpyDirection<< std::endl;
  int ret =  cudaMemcpyAsync(dest, src, memByte, cudaMemcpyHostToDevice,0);
  if(ret != cudaSuccess)
    {
      CUERR;
      printf("**** ERROR *** : c_cuda_gpu_send_pi\n");
    }

  if(cudaStreamSynchronize(0) != cudaSuccess)
    {
      CUERR;
      printf("**** ERROR *** : c_cuda_gpu_recv_pi STREAM\n");
    }

  return ret;

}

int c_cuda_gpu_recv_pi(void *dest, const void *src,  size_t memByte)
{ 
  int ret = cudaMemcpyAsync(dest, src, memByte, cudaMemcpyDeviceToHost,0);
  if(ret != cudaSuccess)
    {
      CUERR;
      printf("**** ERROR *** : c_cuda_gpu_recv_pi\n");
    }

    if(cudaStreamSynchronize(0) != cudaSuccess)
    {
      CUERR;
      printf("**** ERROR *** : c_cuda_gpu_recv_pi STREAM\n");
      }
  return ret;

}

int c_cuda_setdevice(int device)
{
  return cudaSetDevice(device);
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
/****/
