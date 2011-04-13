//! @file
//!  Check streams (OpenCL)
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include "OpenCL_wrappers.h"
#include "Stream.h"
#include "MagicFilter.h"
#define SIZE_I 128
#define NB_STREAM 8

inline void magicfilter_generic_stream(cl_kernel kernel, ocl_stream stream, cl_uint n,cl_uint ndat, cl_mem psi, cl_mem out){
    cl_int ciErrNum;
    int FILTER_WIDTH=16;
    assert(n>=FILTER_WIDTH);
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;
    cl_uint i = 0;
    clSetKernelArg(kernel, i++,sizeof(n), (void*)&n);
    clSetKernelArg(kernel, i++,sizeof(ndat), (void*)&ndat);
    clSetKernelArg(kernel, i++,sizeof(psi), (void*)&psi);
    clSetKernelArg(kernel, i++,sizeof(out), (void*)&out);
    clSetKernelArg(kernel, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,n), shrRoundUp(block_size_j,ndat)};
    ciErrNum = oclEnstreamNDRangeKernel  (stream, kernel, 2, NULL, globalWorkSize, localWorkSize);
    oclErrorCheck(ciErrNum,"Failed to enqueue magic filter kernel!");
}

int main() {
  int i,j;
  size_t size;
  cl_uint n, ndat;
  cl_context context;
  cl_command_queue queue;
  cl_platform_id platform_id;
  clGetPlatformIDs(1, &platform_id, NULL);
  cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform_id, 0 };

  cl_int ciErrNum = CL_SUCCESS;
  context = clCreateContextFromType(properties, CL_DEVICE_TYPE_GPU, NULL, NULL, &ciErrNum);
  oclErrorCheck(ciErrNum,"Failed to create GPU context!");

  size_t nContextDescriptorSize;
  clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, 0, &nContextDescriptorSize);
  cl_device_id * aDevices = (cl_device_id *) malloc(nContextDescriptorSize);
  clGetContextInfo(context, CL_CONTEXT_DEVICES, nContextDescriptorSize, aDevices, 0);
  queue = clCreateCommandQueue(context, aDevices[0], 0, &ciErrNum);
  oclErrorCheck(ciErrNum,"Failed to create command queue!");

  ciErrNum = oclInitStreams(context);
  oclErrorCheck(ciErrNum,"Failed to init streams!");
  build_magicfilter_programs(&context);
  struct bigdft_kernels kernels;
  create_magicfilter_kernels(&kernels);

  double * data[NB_STREAM];
  double * results[NB_STREAM];
  cl_mem input[NB_STREAM];
  cl_mem output[NB_STREAM];
  ocl_stream streams[NB_STREAM];
  

  size = sizeof(double) * SIZE_I * SIZE_I * SIZE_I;
  for(i=0; i<NB_STREAM; i++) {
    data[i] = (double *)malloc(size);
    results[i] = (double *)malloc(size);
    input[i] = clCreateBuffer( context, CL_MEM_READ_WRITE, size, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create read buffer!");
    output[i] = clCreateBuffer( context, CL_MEM_READ_WRITE, size, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create write buffer!");
    streams[i] = oclCreateStream(queue, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create stream!");
    for(j=0; j < SIZE_I * SIZE_I * SIZE_I; j++)
      data[i][j] = (double)rand() / (double)RAND_MAX;
  }

  printf("Enstreaming writes...\n");
  for(i=0; i<NB_STREAM; i++) {
    ciErrNum = oclEnstreamWriteBuffer(streams[i], input[i], 0, size, data[i]);
    oclErrorCheck(ciErrNum,"Failed to enstream write buffer!");
  }
  printf("Enstreaming kernels...\n");
  for(i=0; i<NB_STREAM; i++) {
    for(j=0; j<500; j++){
       magicfilter_generic_stream(kernels.magicfilter1d_kernel_d, streams[i], SIZE_I, SIZE_I * SIZE_I, input[i], output[i]);
       magicfilter_generic_stream(kernels.magicfilter1d_kernel_d, streams[i], SIZE_I, SIZE_I * SIZE_I, output[i], input[i]);
    }
  }
  printf("Enstreaming reads...\n");
  for(i=0; i<NB_STREAM; i++) {
    ciErrNum = oclEnstreamReadBuffer(streams[i], output[i], 0, size, results[i]);
    oclErrorCheck(ciErrNum,"Failed to enstream read buffer!");
  }
  printf("Waiting for kernels to finish...\n");
  for(i=0; i<NB_STREAM; i++) {
    ciErrNum = oclStreamFinish(streams[i]);
    oclErrorCheck(ciErrNum,"Failed to finish stream!");
  }
  printf("Streams finished.\n");
  printf("Enqueuing writes...\n");
  for(i=0; i<NB_STREAM; i++) {
    ciErrNum = clEnqueueWriteBuffer(queue, input[i], CL_FALSE, 0, size, data[i],0,NULL,NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue write buffer!");
  }
  printf("Enqueuing kernels...\n");
  n = SIZE_I;
  ndat = SIZE_I * SIZE_I;
  for(i=0; i<NB_STREAM; i++) {
    for(j=0; j<500; j++){
       magicfilter_generic(kernels.magicfilter1d_kernel_d, queue, &n, &ndat, &(input[i]), &(output[i]));
       magicfilter_generic(kernels.magicfilter1d_kernel_d, queue, &n, &ndat, &(output[i]), &(input[i]));
    }
    ciErrNum = clFlush(queue);
    oclErrorCheck(ciErrNum,"Failed to flush queue!");
  }
  printf("Enqueuing reads...\n");
  for(i=0; i<NB_STREAM; i++) {
    ciErrNum = clEnqueueReadBuffer(queue, output[i], CL_FALSE, 0, size, results[i],0,NULL,NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue read buffer!");
  }
  printf("Waiting for kernels to finish...\n");
  ciErrNum = clFinish(queue);
  oclErrorCheck(ciErrNum,"Failed to finish queue!");
  printf("Queue finished.\n");
  for(i=0; i<NB_STREAM; i++) { 
    ciErrNum = oclReleaseStream(streams[i]);
    oclErrorCheck(ciErrNum,"Failed to release stream!");
    ciErrNum = clReleaseMemObject(input[i]);
    oclErrorCheck(ciErrNum,"Failed to release buffer!");
    ciErrNum = clReleaseMemObject(output[i]);
    oclErrorCheck(ciErrNum,"Failed to release buffer!");
  }
  ciErrNum = oclEndStreams();
  oclErrorCheck(ciErrNum,"Failed to end streams!");
  clean_magicfilter_kernels(&kernels);
  clean_magicfilter_programs();
  ciErrNum = clReleaseCommandQueue(queue);
  oclErrorCheck(ciErrNum,"Failed to release command queue!");
  ciErrNum = clReleaseContext(context);
  oclErrorCheck(ciErrNum,"Failed to release context!");
  return 0;
}

