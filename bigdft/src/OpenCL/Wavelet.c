//! @file
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include "Wavelet.h"
#include "OpenCL_wrappers.h"
#include "Wavelet_Generator.h"


inline void ana_block_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out){
    cl_int ciErrNum;
    int FILTER_WIDTH = 16;
    int ELEM_PER_THREAD=2;
    assert(*n>=FILTER_WIDTH);
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH/ELEM_PER_THREAD;
    cl_uint i = 0;
    clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    clSetKernelArg(kernel, i++,sizeof(double)*block_size_j*ELEM_PER_THREAD*(block_size_i*2+FILTER_WIDTH + 1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j*ELEM_PER_THREAD,*ndat)*block_size_j/FILTER_WIDTH};
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue analysis kernel!");

}

inline void ana_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out){
    cl_int ciErrNum;
    int FILTER_WIDTH = 16;
    assert(*n>=FILTER_WIDTH);
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;
    cl_uint i = 0;
    clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue analysis kernel!");

}

inline void syn_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out) {
    cl_int ciErrNum;
    int FILTER_WIDTH = 8;
    int SIZE_I = 2*FILTER_WIDTH;
    assert(*n>=SIZE_I);
    size_t block_size_i=SIZE_I, block_size_j=SIZE_I;
    cl_uint i = 0;
    clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue synthesis kernel!");
}

//cl_program anaProgram;
//cl_program synProgram;

void create_wavelet_kernels(bigdft_context * context, struct bigdft_kernels * kernels) {
    cl_int ciErrNum = CL_SUCCESS;
    kernels->anashrink1d_kernel_d=clCreateKernel((*context)->anaProgram,"anashrink1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create anashrink1dKernel_d kernel!");
    kernels->ana1d_kernel_d=clCreateKernel((*context)->anaProgram,"ana1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create ana1dKernel_d kernel!");
    kernels->ana1d_block_kernel_d=clCreateKernel((*context)->anaProgram,"ana1d_blockKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create ana1d_blockKernel_d kernel!");
    kernels->syngrow1d_kernel_d=clCreateKernel((*context)->synProgram,"syngrow1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create syngrow1dKernel_d kernel!");
    kernels->syn1d_kernel_d=clCreateKernel((*context)->synProgram,"syn1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create syn1dKernel_d kernel!");
}

void build_wavelet_programs(bigdft_context * context){
    struct bigdft_device_infos infos;
    get_context_devices_infos(context, &infos);
    cl_int ciErrNum = CL_SUCCESS;
    char * code = generate_ana_program(&infos);
    (*context)->anaProgram = clCreateProgramWithSource((*context)->context, 1, (const char**) &code, NULL, &ciErrNum);
    free(code);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram((*context)->anaProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build ana program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo((*context)->anaProgram, oclGetFirstDev((*context)->context), CL_PROGRAM_BUILD_LOG, sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    code = generate_syn_program(&infos);
    (*context)->synProgram = clCreateProgramWithSource((*context)->context, 1, (const char**) &code, NULL, &ciErrNum);
    free(code);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram((*context)->synProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build syn program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo((*context)->synProgram, oclGetFirstDev((*context)->context), CL_PROGRAM_BUILD_LOG, sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
}

void FC_FUNC_(anashrink1d_d,ANASHRINK1D_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
  ana_generic((*command_queue)->kernels.anashrink1d_kernel_d, (*command_queue)->command_queue, n, ndat, psi, out);
}

void FC_FUNC_(ana1d_d,ANA1D_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
  ana_generic((*command_queue)->kernels.ana1d_kernel_d, (*command_queue)->command_queue, n, ndat, psi, out);
}

void FC_FUNC_(ana1d_block_d,ANA1D_BLOCK_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
  ana_block_generic((*command_queue)->kernels.ana1d_block_kernel_d, (*command_queue)->command_queue, n, ndat, psi, out);
}

void FC_FUNC_(ana_d_generic,ANA_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *psi, cl_mem *out){
  cl_uint n, ndat;
  cl_uint n1 = dimensions[0];
  cl_uint n2 = dimensions[1];
  cl_uint n3 = dimensions[2];
  if( !periodic[0] ) n1 += 7;
  if( !periodic[1] ) n2 += 7;
  if( !periodic[2] ) n3 += 7;
  ndat = n2 * n1 * 4;
  if( periodic[2] ) {
    n = n3;
    ana_generic((*command_queue)->kernels.ana1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, psi, out);
  } else {
    n3 -= 7;
    n = n3;
    ana_generic((*command_queue)->kernels.anashrink1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, psi, out);
  }
  ndat = n1 * n3 * 4;
  if( periodic[1] ) {
    n = n2;
    ana_generic((*command_queue)->kernels.ana1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, out, tmp);
  } else {
    n2 -= 7;
    n = n2;
    ana_generic((*command_queue)->kernels.anashrink1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, out, tmp);
  }
  ndat = n2 * n3 * 4;
  if( periodic[0] ) {
    n = n1;
    ana_generic((*command_queue)->kernels.ana1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, tmp, out);
  } else {
    n1 -= 7;
    n = n1;
    ana_generic((*command_queue)->kernels.anashrink1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, tmp, out);
  }
}

void FC_FUNC_(ana_self_d_generic,ANA_SELF_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *psi, cl_mem *out){
  cl_uint n, ndat;
  cl_uint n1 = dimensions[0];
  cl_uint n2 = dimensions[1];
  cl_uint n3 = dimensions[2];
  if( !periodic[0] ) n1 += 7;
  if( !periodic[1] ) n2 += 7;
  if( !periodic[2] ) n3 += 7;
  ndat = n2 * n1 * 4;
  if( periodic[2] ) {
    n = n3;
    ana_generic((*command_queue)->kernels.ana1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, psi, out);
  } else {
    n3 -= 7;
    n = n3;
    ana_generic((*command_queue)->kernels.anashrink1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, psi, out);
  }
  ndat = n1 * n3 * 4;
  if( periodic[1] ) {
    n = n2;
    ana_generic((*command_queue)->kernels.ana1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, out, psi);
  } else {
    n2 -= 7;
    n = n2;
    ana_generic((*command_queue)->kernels.anashrink1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, out, psi);
  }
  ndat = n2 * n3 * 4;
  if( periodic[0] ) {
    n = n1;
    ana_generic((*command_queue)->kernels.ana1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, psi, out);
  } else {
    n1 -= 7;
    n = n1;
    ana_generic((*command_queue)->kernels.anashrink1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, psi, out);
  }
}

void FC_FUNC_(ana_d,ANA_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out){
  cl_uint n1 = dimensions[0];
  cl_uint n2 = dimensions[1];
  cl_uint n3 = dimensions[2];
  cl_uint ndat = n2 * n1 * 4;
  ana_generic((*command_queue)->kernels.ana1d_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, out);
  ndat = n1 * n3 * 4;
  ana_generic((*command_queue)->kernels.ana1d_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, tmp);
  ndat = n2 * n3 * 4;
  ana_generic((*command_queue)->kernels.ana1d_kernel_d, (*command_queue)->command_queue, &n1, &ndat, tmp, out);
}

void FC_FUNC_(ana_block_d,ANA_BLOCK_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out){
  cl_uint n1 = dimensions[0];
  cl_uint n2 = dimensions[1];
  cl_uint n3 = dimensions[2];
  cl_uint ndat = n2 * n1 * 4;
  ana_block_generic((*command_queue)->kernels.ana1d_block_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, out);
  ndat = n1 * n3 * 4;
  ana_block_generic((*command_queue)->kernels.ana1d_block_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, tmp);
  ndat = n2 * n3 * 4;
  ana_block_generic((*command_queue)->kernels.ana1d_block_kernel_d, (*command_queue)->command_queue, &n1, &ndat, tmp, out);
}

void FC_FUNC_(ana_self_d,ANA_SELF_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *psi, cl_mem *out){
  cl_uint n1 = dimensions[0];
  cl_uint n2 = dimensions[1];
  cl_uint n3 = dimensions[2];
  cl_uint ndat = n2 * n1 * 4;
  ana_generic((*command_queue)->kernels.ana1d_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, out);
  ndat = n1 * n3 * 4;
  ana_generic((*command_queue)->kernels.ana1d_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, psi);
  ndat = n2 * n3 * 4;
  ana_generic((*command_queue)->kernels.ana1d_kernel_d, (*command_queue)->command_queue, &n1, &ndat, psi, out);
}

void FC_FUNC_(syngrow1d_d,SYNGROW1D_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_uint n1 = *n+7;
    syn_generic((*command_queue)->kernels.syngrow1d_kernel_d, (*command_queue)->command_queue, &n1, ndat, psi, out);
}

void FC_FUNC_(syn1d_d,SYN1D_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_uint *ndat,cl_mem *psi,cl_mem *out){
    syn_generic((*command_queue)->kernels.syn1d_kernel_d, (*command_queue)->command_queue, n, ndat, psi, out);
}

void FC_FUNC_(syn_d_generic,SYN_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *psi, cl_mem *out){
  cl_uint n, ndat;
  cl_uint n1 = dimensions[0];
  cl_uint n2 = dimensions[1];
  cl_uint n3 = dimensions[2];
  ndat = n2 * n1 * 4;
  if( periodic[2] ) {
    n = n3;
    syn_generic((*command_queue)->kernels.syn1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, psi, out);
  } else {
    n3 += 7;
    n = n3;
    syn_generic((*command_queue)->kernels.syngrow1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, psi, out);
  }
  ndat = n1 * n3 * 4;
  if( periodic[1] ) {
    n = n2;
    syn_generic((*command_queue)->kernels.syn1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, out, tmp);
  } else {
    n2 += 7;
    n = n2;
    syn_generic((*command_queue)->kernels.syngrow1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, out, tmp);
  }
  ndat = n2 * n3 * 4;
  if( periodic[0] ) {
    n = n1;
    syn_generic((*command_queue)->kernels.syn1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, tmp, out);
  } else {
    n1 += 7;
    n = n1;
    syn_generic((*command_queue)->kernels.syngrow1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, tmp, out);
  }
}

void FC_FUNC_(syn_self_d_generic,SYN_SELF_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *psi, cl_mem *out){
  cl_uint n, ndat;
  cl_uint n1 = dimensions[0];
  cl_uint n2 = dimensions[1];
  cl_uint n3 = dimensions[2];
  ndat = n2 * n1 * 4;
  if( periodic[2] ) {
    n = n3;
    syn_generic((*command_queue)->kernels.syn1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, psi, out);
  } else {
    n3 += 7;
    n = n3;
    syn_generic((*command_queue)->kernels.syngrow1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, psi, out);
  }
  ndat = n1 * n3 * 4;
  if( periodic[1] ) {
    n = n2;
    syn_generic((*command_queue)->kernels.syn1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, out, psi);
  } else {
    n2 += 7;
    n = n2;
    syn_generic((*command_queue)->kernels.syngrow1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, out, psi);
  }
  ndat = n2 * n3 * 4;
  if( periodic[0] ) {
    n = n1;
    syn_generic((*command_queue)->kernels.syn1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, psi, out);
  } else {
    n1 += 7;
    n = n1;
    syn_generic((*command_queue)->kernels.syngrow1d_kernel_d, (*command_queue)->command_queue, &n, &ndat, psi, out);
  }
}

void FC_FUNC_(syn_d,SYN_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out){
  cl_uint n1 = dimensions[0];
  cl_uint n2 = dimensions[1];
  cl_uint n3 = dimensions[2];
  cl_uint ndat = n2 * n1 * 4;
  syn_generic((*command_queue)->kernels.syn1d_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, out);
  ndat = n1 * n3 * 4;
  syn_generic((*command_queue)->kernels.syn1d_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, tmp);
  ndat = n2 * n3 * 4;
  syn_generic((*command_queue)->kernels.syn1d_kernel_d, (*command_queue)->command_queue, &n1, &ndat, tmp, out);
}
void FC_FUNC_(syn_self_d,SYN_SELF_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *psi, cl_mem *out){
  cl_uint n1 = dimensions[0];
  cl_uint n2 = dimensions[1];
  cl_uint n3 = dimensions[2];
  cl_uint ndat = n2 * n1 * 4;
  syn_generic((*command_queue)->kernels.syn1d_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, out);
  ndat = n1 * n3 * 4;
  syn_generic((*command_queue)->kernels.syn1d_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, psi);
  ndat = n2 * n3 * 4;
  syn_generic((*command_queue)->kernels.syn1d_kernel_d, (*command_queue)->command_queue, &n1, &ndat, psi, out);
}

void clean_wavelet_kernels(struct bigdft_kernels * kernels){
  cl_int ciErrNum;
  ciErrNum = clReleaseKernel(kernels->ana1d_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->ana1d_block_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->anashrink1d_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->syn1d_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->syngrow1d_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
}

void clean_wavelet_programs(bigdft_context * context){
  cl_int ciErrNum;
  ciErrNum = clReleaseProgram((*context)->anaProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
  ciErrNum = clReleaseProgram((*context)->synProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
}
