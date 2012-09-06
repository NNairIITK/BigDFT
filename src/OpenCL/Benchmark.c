//! @file
//!  Those kernels are benchmark, they try to focus on
//!  1 performance aspect each.
//!  benchmark_mops measures the number of memory copy
//!  that can be processed. It generates n*8 double memory
//!  copy.
//!  benchmark_flops mesures the number of mad that can be
//!  expected from the device. It generates 16*128*n mads,
//!  or 16*128*n*2 flop.
//!  transpose measures the performances of the transposition
//!  method used in almost every kernel.
//! 
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include "OpenCL_wrappers.h"
#include "Benchmark_Generator.h"

cl_program benchmarkProgram;

void create_benchmark_kernels(struct bigdft_kernels * kernels){
    cl_int ciErrNum = CL_SUCCESS;
    kernels->benchmark_mops_kernel_d=clCreateKernel(benchmarkProgram,"benchmark_mopsKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create benchmark_mopsKernel_d kernel!");
    kernels->benchmark_flops_kernel_d=clCreateKernel(benchmarkProgram,"benchmark_flopsKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create benchmark_flopsKernel_d kernel!");
    kernels->transpose_kernel_d=clCreateKernel(benchmarkProgram,"transposeKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create transposeKernel_d kernel!");
    kernels->notranspose_kernel_d=clCreateKernel(benchmarkProgram,"notransposeKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create notransposeKernel_d kernel!");
}

void build_benchmark_programs(bigdft_context * context){
    struct bigdft_device_infos infos;
    get_context_devices_infos(context, &infos);
    cl_int ciErrNum=CL_SUCCESS;
    char * code = generate_benchmark_program(&infos);
    benchmarkProgram = clCreateProgramWithSource((*context)->context,1,(const char**) &code, NULL, &ciErrNum);
    free(code);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(benchmarkProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to build benchmark program!\n",ciErrNum);
        char cBuildLog[10240];
        clGetProgramBuildInfo(benchmarkProgram, oclGetFirstDev((*context)->context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
}

inline void benchmark_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n,cl_mem *in,cl_mem *out){
    cl_int ciErrNum;
    int FILTER_WIDTH=64;
    assert(*n>=FILTER_WIDTH);
    size_t block_size_i=FILTER_WIDTH;
    cl_uint i = 0;
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*in), (void*)in);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    size_t localWorkSize[] = { block_size_i };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n)};
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue benchmark kernel!");
}

inline void transpose_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    int FILTER_WIDTH=16;
    assert(*n>=FILTER_WIDTH);
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;
    cl_uint i = 0;
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(cl_double)*block_size_j*(block_size_i+1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue transpose kernel!");
}

void FC_FUNC_(benchmark_flops_d,BENCHMARK_FLOPS_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *out){
    benchmark_generic((*command_queue)->kernels.benchmark_flops_kernel_d, (*command_queue)->command_queue, n, in, out);
}

void FC_FUNC_(benchmark_mops_d,BENCHMARK_MOPS_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *out){
    benchmark_generic((*command_queue)->kernels.benchmark_mops_kernel_d, (*command_queue)->command_queue, n, in, out);
}

void FC_FUNC_(transpose_d,TRANSPOSE_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    transpose_generic((*command_queue)->kernels.transpose_kernel_d, (*command_queue)->command_queue, n, ndat, psi, out);
}

void FC_FUNC_(notranspose_d,NOTRANSPOSE_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    transpose_generic((*command_queue)->kernels.notranspose_kernel_d, (*command_queue)->command_queue, n, ndat, psi, out);
}

void clean_benchmark_kernels(struct bigdft_kernels * kernels){
  cl_int ciErrNum;
  ciErrNum = clReleaseKernel(kernels->benchmark_flops_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->benchmark_mops_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->transpose_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->notranspose_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
}

void clean_benchmark_programs(){
  cl_int ciErrNum;
  ciErrNum = clReleaseProgram(benchmarkProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
}
