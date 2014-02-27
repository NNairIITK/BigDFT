//! @file
/*!
  Magicfilter kernels share commons features. They are composed
  of one convolution and a transposition. Each work item is
  responsible for processing one element of the result matrix.
  The convolution filter is 16 elemnts long, and buffers are 33*16
  elements long to avoid bank conflicts. Local work size is 16*16
  elements, so each work item is responsible for loading 2 elements.
  The first kernel is commented.
*/
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include "MagicFilter.h"
#include "OpenCL_wrappers.h"
#include "MagicFilter_Generator.h"

inline void magicfilter_block_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    int FILTER_WIDTH=16;
    int ELEM_PER_THREAD=2;
    assert(*n>=FILTER_WIDTH);
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH/ELEM_PER_THREAD;
    cl_uint i = 0;
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j*ELEM_PER_THREAD,*ndat)*block_size_j/FILTER_WIDTH};
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue magic filter kernel!");
}

inline void magicfilter_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    int FILTER_WIDTH=16;
    assert(*n>=FILTER_WIDTH);
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;
    cl_uint i = 0;
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue magic filter kernel!");
}

inline void magicfilter_pot_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *pot, cl_mem *out) {
    cl_int ciErrNum;
    int FILTER_WIDTH = 16;
    assert(*n>=FILTER_WIDTH);
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;
    cl_uint i = 0;
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*pot), (void*)pot);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue magic filter pot kernel!");
}

cl_program magicfilterProgram;

void create_magicfilter_kernels(struct bigdft_kernels * kernels){
    cl_int ciErrNum = CL_SUCCESS;
    kernels->magicfiltergrow1d_kernel_d=clCreateKernel(magicfilterProgram,"magicfiltergrow1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfiltergrow1dKernel_d kernel!");
    kernels->magicfiltergrow1d_den_kernel_d=clCreateKernel(magicfilterProgram,"magicfiltergrow1d_denKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfiltergrow1d_denKernel_d kernel!");
    kernels->magicfiltershrink1d_kernel_d=clCreateKernel(magicfilterProgram,"magicfiltershrink1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfiltershrink1dKernel_d kernel!");
    kernels->magicfiltergrow1d_pot_kernel_d=clCreateKernel(magicfilterProgram,"magicfiltergrow1d_potKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfiltergrow1d_potKernel_d kernel!");
    kernels->magicfilter1d_kernel_d=clCreateKernel(magicfilterProgram,"magicfilter1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfilter1dKernel_d kernel!");
    kernels->magicfilter1d_den_kernel_d=clCreateKernel(magicfilterProgram,"magicfilter1d_denKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfilter1d_denKernel_d kernel!");
    kernels->magicfilter1d_pot_kernel_d=clCreateKernel(magicfilterProgram,"magicfilter1d_potKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfilter1d_potKernel_d kernel!");
    kernels->magicfilter1d_t_kernel_d=clCreateKernel(magicfilterProgram,"magicfilter1d_tKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfilter1d_tKernel_d kernel!");
    kernels->magicfilter1d_straight_kernel_d=clCreateKernel(magicfilterProgram,"magicfilter1d_straightKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfilter1d_straightKernel_d kernel!");
    kernels->magicfilter1d_block_kernel_d=clCreateKernel(magicfilterProgram,"magicfilter1d_blockKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfilter1d_blockKernel_d kernel!");
}

void build_magicfilter_programs(bigdft_context * context){
    struct bigdft_device_infos infos;
    get_context_devices_infos(context, &infos);
    cl_int ciErrNum=CL_SUCCESS;
    char * code = generate_magicfilter_program(&infos);

    magicfilterProgram = clCreateProgramWithSource((*context)->context, 1, (const char**) &code, NULL, &ciErrNum);
    free(code);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(magicfilterProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to build magicfilter program!\n", ciErrNum);
        char cBuildLog[10240];
        clGetProgramBuildInfo(magicfilterProgram, oclGetFirstDev((*context)->context), CL_PROGRAM_BUILD_LOG, sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr, "%s\n", cBuildLog);
        exit(1);
    }
}

void FC_FUNC_(magicfiltershrink1d_d,MAGICFILTERSHRINK1D_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    magicfilter_generic((*command_queue)->kernels.magicfiltershrink1d_kernel_d,(*command_queue)->command_queue,n,ndat,psi,out);
}

void FC_FUNC_(magicfiltergrow1d_d,MAGICFILTERGROW1D_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_uint n1 = *n + 15;
    magicfilter_generic((*command_queue)->kernels.magicfiltergrow1d_kernel_d, (*command_queue)->command_queue, &n1, ndat, psi, out);
}

void FC_FUNC_(magicfilter1d_d,MAGICFILTER1D_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue, n, ndat, psi, out);
}

void FC_FUNC_(magicfilter1d_straight_d,MAGICFILTER1D_STRAIGHT_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_straight_kernel_d, (*command_queue)->command_queue, n, ndat, psi, out);
}

void FC_FUNC_(magicfilter1d_block_d,MAGICFILTER1D_BLOCK_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    magicfilter_block_generic((*command_queue)->kernels.magicfilter1d_block_kernel_d, (*command_queue)->command_queue, n, ndat, psi, out);
}

void FC_FUNC_(magicfilter1d_pot_d,MAGICFILTER1D_POT_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *pot, cl_mem *out){
    magicfilter_pot_generic((*command_queue)->kernels.magicfilter1d_pot_kernel_d, (*command_queue)->command_queue, n, ndat, psi, pot, out);
}

void FC_FUNC_(magicfilter1d_t_d,MAGICFILTER1D_T_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_t_kernel_d,(*command_queue)->command_queue,n,ndat,psi,out);
}

void FC_FUNC_(magicfilter_n_self_d,MAGICFILTER_N_SELF_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *psi, cl_mem *out){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, out);
    ndat = n1 * n3;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, psi);
    ndat = n2 * n3;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue, &n1, &ndat, psi, out);
}

void FC_FUNC_(magicfilter_n_d,MAGICFILTER_N_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, out);
    ndat = n1 * n3;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, tmp);
    ndat = n2 * n3;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue, &n1, &ndat, tmp, out);
}

void FC_FUNC_(magicfilter_n_straight_d,MAGICFILTER_N_STRAIGHT_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n2 * n3;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_straight_kernel_d, (*command_queue)->command_queue, &n1, &ndat, psi, out);
    ndat = n1 * n3;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_straight_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, tmp);
    ndat = n1 * n2;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_straight_kernel_d, (*command_queue)->command_queue, &n3, &ndat, tmp, out);
}

void FC_FUNC_(magicfilter_n_block_d,MAGICFILTER_N_BLOCK_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    magicfilter_block_generic((*command_queue)->kernels.magicfilter1d_block_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, out);
    ndat = n1 * n3;
    magicfilter_block_generic((*command_queue)->kernels.magicfilter1d_block_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, tmp);
    ndat = n2 * n3;
    magicfilter_block_generic((*command_queue)->kernels.magicfilter1d_block_kernel_d, (*command_queue)->command_queue, &n1, &ndat, tmp, out);
}

void FC_FUNC_(magicfilter_den_d,MAGICFILTER_DEN_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out){
    cl_uint n1 = dimensions[0]*2;
    cl_uint n2 = dimensions[1]*2;
    cl_uint n3 = dimensions[2]*2;
    cl_uint ndat = n1 * n2;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, out);
    ndat = n1 * n3;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, tmp);
    ndat = n2 * n3;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_den_kernel_d, (*command_queue)->command_queue, &n1, &ndat, tmp, out);
}

void FC_FUNC_(magicfilter_den_d_generic,MAGICFILTER_DEN_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *psi, cl_mem *out){
    cl_uint n1 = dimensions[0]*2;
    cl_uint n2 = dimensions[1]*2;
    cl_uint n3 = dimensions[2]*2;
    if( !periodic[0] ) n1 += 14;
    if( !periodic[1] ) n2 += 14;
    if( !periodic[2] ) n3 += 14;
    cl_uint ndat = n1 * n2;
    if( periodic[2] ) {
      magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue,  &n3, &ndat, psi, out);
    } else {
      n3 += 15;
      magicfilter_generic((*command_queue)->kernels.magicfiltergrow1d_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, out);
    }
    ndat = n1 * n3;
    if( periodic[1] ) {
      magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue,  &n2, &ndat, out, tmp);
    } else {
      n2 += 15;
      magicfilter_generic((*command_queue)->kernels.magicfiltergrow1d_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, tmp);
    }
    ndat = n2 * n3;
    if( periodic[0] ) {
      magicfilter_generic((*command_queue)->kernels.magicfilter1d_den_kernel_d, (*command_queue)->command_queue, &n1, &ndat, tmp, out);
    } else {
      n1 += 15;
      magicfilter_generic((*command_queue)->kernels.magicfiltergrow1d_den_kernel_d, (*command_queue)->command_queue, &n1, &ndat, tmp, out);
    }
}

void FC_FUNC_(magicfilter_t_self_d,MAGICFILTER_T_SELF_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *psi, cl_mem *out){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_t_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, out);
    ndat = n1 * n3;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_t_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, psi);
    ndat = n2 * n3;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_t_kernel_d, (*command_queue)->command_queue, &n1, &ndat, psi, out);
}

void FC_FUNC_(magicfilter_t_d,MAGICFILTER_T_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_t_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, out);
    ndat = n1 * n3;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_t_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, tmp);
    ndat = n2 * n3;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_t_kernel_d, (*command_queue)->command_queue, &n1, &ndat, tmp, out);
}

void FC_FUNC_(magic_filter_3d_generic,MAGIC_FILTER_3D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *tmp_dot, cl_mem *psi, cl_mem *out) {
  cl_uint ndat;
  cl_uint n1, n2, n3;
  n1 = dimensions[0] * 2;
  n2 = dimensions[1] * 2;
  n3 = dimensions[2] * 2;
  if( !periodic[0] ) n1 += 14;
  if( !periodic[1] ) n2 += 14;
  if( !periodic[2] ) n3 += 14;
  ndat = n1 * n2;
  if( periodic[2] ) {
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue,  &n3, &ndat, psi, tmp);
  } else {
    n3 += 15;
    magicfilter_generic((*command_queue)->kernels.magicfiltergrow1d_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, tmp);
  }
  ndat = n1 * n3;
  if( periodic[1] ) {
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue,  &n2, &ndat, tmp, tmp_dot);
  } else {
    n2 += 15;
    magicfilter_generic((*command_queue)->kernels.magicfiltergrow1d_kernel_d, (*command_queue)->command_queue, &n2, &ndat, tmp, tmp_dot);
  }
  ndat =  n2 * n3;
  if( periodic[0] ) {
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue,  &n1, &ndat, tmp_dot, out);
  } else {
    n1 += 15;
    magicfilter_generic((*command_queue)->kernels.magicfiltergrow1d_kernel_d, (*command_queue)->command_queue, &n1, &ndat, tmp_dot, out);
  }
}

void FC_FUNC_(magic_filter_t_3d_generic,MAGIC_FILTER_T_3D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *tmp_dot, cl_mem *psi, cl_mem *out) {
  cl_uint ndat;
  cl_uint n1, n2, n3;
  n1 = dimensions[0] * 2;
  n2 = dimensions[1] * 2;
  n3 = dimensions[2] * 2;
  if( !periodic[0] ) n1 += 14 + 15;
  if( !periodic[1] ) n2 += 14 + 15;
  if( !periodic[2] ) n3 += 14 + 15;
  ndat = n1 * n2;
  if( periodic[2] ) {
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_t_kernel_d, (*command_queue)->command_queue,  &n3, &ndat, psi, tmp);
  } else {
    n3 -= 15;
    magicfilter_generic((*command_queue)->kernels.magicfiltershrink1d_kernel_d, (*command_queue)->command_queue,  &n3, &ndat, psi, tmp);
  }
  ndat = n1 * n3;
  if( periodic[1] ) {
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_t_kernel_d, (*command_queue)->command_queue,  &n2, &ndat, tmp, tmp_dot);
  } else {
    n2 -= 15;
    magicfilter_generic((*command_queue)->kernels.magicfiltershrink1d_kernel_d, (*command_queue)->command_queue,  &n2, &ndat, tmp, tmp_dot);
  }
  ndat = n2 * n3;
  if( periodic[0] ) {
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_t_kernel_d, (*command_queue)->command_queue,  &n1, &ndat, tmp_dot, out);
  } else {
    n1 -= 15;
    magicfilter_generic((*command_queue)->kernels.magicfiltershrink1d_kernel_d, (*command_queue)->command_queue,  &n1, &ndat, tmp_dot, out);
  }
}

void FC_FUNC_(potential_application_d_generic,POTENTIAL_APPLICATION_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *tmp_dot, cl_mem *psi, cl_mem *out, cl_mem *pot, double *epot) {
  cl_uint ndat;
  cl_uint n1, n2, n3;
  n1 = dimensions[0] * 2;
  n2 = dimensions[1] * 2;
  n3 = dimensions[2] * 2;
  if( !periodic[0] ) n1 += 14;
  if( !periodic[1] ) n2 += 14;
  if( !periodic[2] ) n3 += 14;
  ndat = n1 * n2;
  if( periodic[2] ) {
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue,  &n3, &ndat, psi, tmp);
  } else {
    n3 += 15;
    magicfilter_generic((*command_queue)->kernels.magicfiltergrow1d_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, tmp);
  }
  ndat = n1 * n3;
  if( periodic[1] ) {
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue,  &n2, &ndat, tmp, out);
  } else {
    n2 += 15;
    magicfilter_generic((*command_queue)->kernels.magicfiltergrow1d_kernel_d, (*command_queue)->command_queue, &n2, &ndat, tmp, out);
  }
  ndat =  n2 * n3;
  if( periodic[0] ) {
    magicfilter_pot_generic((*command_queue)->kernels.magicfilter1d_pot_kernel_d, (*command_queue)->command_queue,  &n1, &ndat, out, pot, tmp);
  } else {
    n1 += 15;
    magicfilter_pot_generic((*command_queue)->kernels.magicfiltergrow1d_pot_kernel_d, (*command_queue)->command_queue, &n1, &ndat, out, pot, tmp);
  }
  ndat = n1 * n2;
  if( periodic[2] ) {
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_t_kernel_d, (*command_queue)->command_queue,  &n3, &ndat, tmp, out);
  } else {
    n3 -= 15;
    magicfilter_generic((*command_queue)->kernels.magicfiltershrink1d_kernel_d, (*command_queue)->command_queue,  &n3, &ndat, tmp, out);
  }
  ndat = n1 * n3;
  if( periodic[1] ) {
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_t_kernel_d, (*command_queue)->command_queue,  &n2, &ndat, out, tmp);
  } else {
    n2 -= 15;
    magicfilter_generic((*command_queue)->kernels.magicfiltershrink1d_kernel_d, (*command_queue)->command_queue,  &n2, &ndat, out, tmp);
  }
  ndat = n2 * n3;
  if( periodic[0] ) {
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_t_kernel_d, (*command_queue)->command_queue,  &n1, &ndat, tmp, out);
  } else {
    n1 -= 15;
    magicfilter_generic((*command_queue)->kernels.magicfiltershrink1d_kernel_d, (*command_queue)->command_queue,  &n1, &ndat, tmp, out);
  }
  ndat = n1*n2*n3;
  FC_FUNC_(dot_d_async, DOT_D_ASYNC)(command_queue, &ndat, psi, out, tmp, tmp_dot, epot);
}

void FC_FUNC_(potential_application_d,POTENTIAL_APPLICATION_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out, cl_mem *pot) {
    cl_uint n1 = dimensions[0] * 2;
    cl_uint n2 = dimensions[1] * 2;
    cl_uint n3 = dimensions[2] * 2;
    cl_uint ndat = n1 * n2;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, tmp);
    ndat = n1 * n3;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_kernel_d, (*command_queue)->command_queue, &n2, &ndat, tmp, out);
    ndat = n2 * n3;
    magicfilter_pot_generic((*command_queue)->kernels.magicfilter1d_pot_kernel_d, (*command_queue)->command_queue, &n1, &ndat, out, pot, tmp);
    ndat = n1 * n2;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_t_kernel_d, (*command_queue)->command_queue, &n3, &ndat, tmp, out);
    ndat = n1 * n3;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_t_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, tmp);
    ndat = n2 * n3;
    magicfilter_generic((*command_queue)->kernels.magicfilter1d_t_kernel_d, (*command_queue)->command_queue, &n1, &ndat, tmp, out);
}

void clean_magicfilter_kernels(struct bigdft_kernels * kernels){
  cl_int ciErrNum;
  ciErrNum = clReleaseKernel(kernels->magicfilter1d_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->magicfilter1d_den_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->magicfilter1d_pot_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->magicfilter1d_t_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->magicfiltershrink1d_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->magicfiltergrow1d_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->magicfiltergrow1d_den_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->magicfiltergrow1d_pot_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->magicfilter1d_straight_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->magicfilter1d_block_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
}

void clean_magicfilter_programs(){
  cl_int ciErrNum;
  ciErrNum = clReleaseProgram(magicfilterProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
}
