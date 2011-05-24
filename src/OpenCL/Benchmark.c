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

char * benchmark_program = "\
#ifdef cl_khr_fp64\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#elif defined (cl_amd_fp64)\n\
#pragma OPENCL EXTENSION cl_amd_fp64: enable \n\
#endif\n\
#define NB_ITER 16\n\
__kernel void benchmark_mopsKernel_d(uint n, __global const double *in, __global double *out){\n\
size_t i = get_global_id(0);\n\
i = get_group_id(0) == get_num_groups(0) - 1 ? i - ( get_global_size(0) - n ) : i;\n\
out[i] = in[i];\n\
out[i] = in[i];\n\
out[i] = in[i];\n\
out[i] = in[i];\n\
out[i] = in[i];\n\
out[i] = in[i];\n\
out[i] = in[i];\n\
out[i] = in[i];\n\
};\n\
__kernel void benchmark_flopsKernel_d(uint n, __global const double *in, __global double *out){\n\
size_t i = get_global_id(0);\n\
double a = in[i]*1.15;\n\
double b = in[i]*1.16;\n\
i = get_group_id(0) == get_num_groups(0) - 1 ? i - ( get_global_size(0) - n ) : i;\n\
int j=0;\n\
for(j=0;j<NB_ITER;j++){\n\
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
}\n\
out[i]=a+b;\n\
};\n\
#define FILTER_WIDTH 16\n\
__kernel void transposeKernel_d(uint n, uint ndat, __global const double *psi, __global double *out, __local double *tmp ) {\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1);\n\
ptrdiff_t igt = get_group_id(0);\n\
ptrdiff_t jgt = get_group_id(1);\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
igt = ig - i2 + j2;\n\
jgt = jg - j2 + i2;\n\
//If I'm on the outside, select a border element to load\n\
tmp[i2 * (FILTER_WIDTH + 1) + j2] = psi[jgt + igt * ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
out[(jg*n+ig)]=tmp[j2 * (FILTER_WIDTH + 1) + i2];\n\
};\n\
__kernel void notransposeKernel_d(uint n, uint ndat, __global const double *psi, __global double *out, __local double *tmp ) {\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1);\n\
ptrdiff_t igt = get_group_id(0);\n\
ptrdiff_t jgt = get_group_id(1);\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
//igt = ig - i2 + j2;\n\
//jgt = jg - j2 + i2;\n\
//If I'm on the outside, select a border element to load\n\
//tmp[i2 * (FILTER_WIDTH + 1) + j2] = psi[jgt + igt * ndat];\n\
//barrier(CLK_LOCAL_MEM_FENCE);\n\
\
out[jg*n+ig]=psi[jg*n+ig];//tmp[j2 * (FILTER_WIDTH + 1) + i2];\n\
};\n\
";

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

void build_benchmark_programs(cl_context * context){
    cl_int ciErrNum = CL_SUCCESS;
    benchmarkProgram = clCreateProgramWithSource(*context,1,(const char**) &benchmark_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(benchmarkProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to build benchmark program!\n",ciErrNum);
        char cBuildLog[10240];
        clGetProgramBuildInfo(benchmarkProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
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
