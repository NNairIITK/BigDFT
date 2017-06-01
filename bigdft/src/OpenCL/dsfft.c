#include <string.h>
#include "OpenCL_wrappers.h"
#include "fft_generator.h"


cl_uint use_constant_memory=1;

inline void fft_generated_generic(cl_kernel kernel, bigdft_command_queue command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out,cl_mem *cosi){
  cl_int ciErrNum;
  size_t block_size_i, block_size_j, elem_per_thread;
  fft_compute_sizes(&(command_queue->device_infos), *n, &block_size_i, &block_size_j, &elem_per_thread);
  cl_uint i = 0;
  ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  ciErrNum = clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
  ciErrNum = clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
  ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
  if(!use_constant_memory){
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*cosi), (void*)cosi);
  }
  size_t localWorkSize[] = { block_size_i,block_size_j };
//  printf("%lu %lu\n",block_size_i,block_size_j);
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
  ciErrNum = clEnqueueNDRangeKernel(command_queue->command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue fft kernel!");
}

inline void fft_k_generated_generic(cl_kernel kernel, bigdft_command_queue command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out,cl_mem *k,cl_mem *cosi){
  cl_int ciErrNum;
  size_t block_size_i, block_size_j, elem_per_thread;
  fft_compute_sizes(&(command_queue->device_infos), *n, &block_size_i, &block_size_j, &elem_per_thread);
  cl_uint i = 0;
  ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  ciErrNum = clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
  ciErrNum = clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
  ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
  ciErrNum = clSetKernelArg(kernel, i++,sizeof(*k), (void*)k);
  if(!use_constant_memory){
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*cosi), (void*)cosi);
  }
  size_t localWorkSize[] = { block_size_i,block_size_j };
//    printf("%lu %lu\n",block_size_i,block_size_j);
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
  ciErrNum = clEnqueueNDRangeKernel(command_queue->command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue fft kernel!");
}


void FC_FUNC_(customize_fft,CUSTOMIZE_FFT)(bigdft_context * context, cl_uint *dimensions) {
  assert(dimensions[0]<=1024);
  assert(dimensions[1]<=1024);
  assert(dimensions[2]<=1024);
  (*context)->fft_size[0] = dimensions[0];
  (*context)->fft_size[1] = dimensions[1];
  (*context)->fft_size[2] = dimensions[2];
}

//cl_mem cossind0;
//cl_mem cossind1;
//cl_mem cossind2;

void FC_FUNC_(fft1d_d,FFT1D_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    if((*command_queue)->context->fft_size[0] == *n)
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d0_d, *command_queue, n, ndat, psi, out, &((*command_queue)->context->cossind0));
    else if((*command_queue)->context->fft_size[1] == *n)
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d1_d, *command_queue, n, ndat, psi, out, &((*command_queue)->context->cossind1));
    else if((*command_queue)->context->fft_size[2] == *n)
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d2_d, *command_queue, n, ndat, psi, out, &((*command_queue)->context->cossind2));
}

void FC_FUNC_(fft1d_r_d,FFT1D_R_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    if((*command_queue)->context->fft_size[0] == *n)
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d0_r_d, *command_queue, n, ndat, psi, out, &((*command_queue)->context->cossind0));
    else if((*command_queue)->context->fft_size[1] == *n)
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d1_r_d, *command_queue, n, ndat, psi, out, &((*command_queue)->context->cossind1));
    else if((*command_queue)->context->fft_size[2] == *n)
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d2_r_d, *command_queue, n, ndat, psi, out, &((*command_queue)->context->cossind2));
}

void FC_FUNC_(fft3d_d,FFT1D_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,cl_mem *psi,cl_mem *out,cl_mem *tmp){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    assert((*command_queue)->context->fft_size[2]==n3);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d2_d, *command_queue, &n3, &ndat, psi, out, &((*command_queue)->context->cossind2));
    ndat = n1 * n3;
    assert((*command_queue)->context->fft_size[1]==n2);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d1_d, *command_queue, &n2, &ndat, out, tmp, &((*command_queue)->context->cossind1));
    ndat = n2 * n3;
    assert((*command_queue)->context->fft_size[0]==n1);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d0_d, *command_queue, &n1, &ndat, tmp, out, &((*command_queue)->context->cossind0));
}

void FC_FUNC_(fft3d_k_r2c_d,FFT1D_K_R2C_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,cl_mem *psi,cl_mem *out,cl_mem *tmp,cl_mem *k){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    assert((*command_queue)->context->fft_size[2]==n3);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d2_r2c_d, *command_queue, &n3, &ndat, psi, out, &((*command_queue)->context->cossind2));
    ndat = n1 * n3;
    assert((*command_queue)->context->fft_size[1]==n2);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d1_d, *command_queue, &n2, &ndat, out, tmp, &((*command_queue)->context->cossind1));
    ndat = n2 * n3;
    assert((*command_queue)->context->fft_size[0]==n1);
    fft_k_generated_generic((*command_queue)->kernels.fft_kernel_k_d0_d, *command_queue, &n1, &ndat, tmp, out, k, &((*command_queue)->context->cossind0));
}


void FC_FUNC_(fft3d_r_d,FFT1D_R_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,cl_mem *psi,cl_mem *out,cl_mem *tmp){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    assert((*command_queue)->context->fft_size[2]==n3);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d2_r_d, *command_queue, &n3, &ndat, psi, out, &((*command_queue)->context->cossind2));
    ndat = n1 * n3;
    assert((*command_queue)->context->fft_size[1]==n2);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d1_r_d, *command_queue, &n2, &ndat, out, tmp, &((*command_queue)->context->cossind1));
    ndat = n2 * n3;
    assert((*command_queue)->context->fft_size[0]==n1);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d0_r_d, *command_queue, &n1, &ndat, tmp, out, &((*command_queue)->context->cossind0));
}

void FC_FUNC_(fft3d_r_c2r_d,FFT1D_R_C2R_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,cl_mem *psi,cl_mem *out,cl_mem *tmp){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    assert((*command_queue)->context->fft_size[2]==n3);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d2_r_d, *command_queue, &n3, &ndat, psi, out, &((*command_queue)->context->cossind2));
    ndat = n1 * n3;
    assert((*command_queue)->context->fft_size[1]==n2);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d1_r_d, *command_queue, &n2, &ndat, out, tmp, &((*command_queue)->context->cossind1));
    ndat = n2 * n3;
    assert((*command_queue)->context->fft_size[0]==n1);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d0_r_c2r_d, *command_queue, &n1, &ndat, tmp, out, &((*command_queue)->context->cossind0));
}

void FC_FUNC_(fft3d_k_r2c_d_generic,FFT1D_K_R2C_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *psi,cl_mem *out,cl_mem *tmp,cl_mem *k){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    if( periodic[2] ){
      assert((*command_queue)->context->fft_size[2]==n3);
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d2_r2c_d, *command_queue, &n3, &ndat, psi, out, &((*command_queue)->context->cossind2));
    } else {
      n3 *= 2;
      assert((*command_queue)->context->fft_size[2]==n3);
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d2_r2c_f_d, *command_queue, &n3, &ndat, psi, out, &((*command_queue)->context->cossind2));
    }
    ndat = n1 * n3;
    if( periodic[1] ){
      assert((*command_queue)->context->fft_size[1]==n2);
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d1_d, *command_queue, &n2, &ndat, out, tmp, &((*command_queue)->context->cossind1));
    } else {
      n2 *= 2;
      assert((*command_queue)->context->fft_size[1]==n2);
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d1_f_d, *command_queue, &n2, &ndat, out, tmp, &((*command_queue)->context->cossind1));
    }
    ndat = n2 * n3;
    if( periodic[0] ){
      assert((*command_queue)->context->fft_size[0]==n1);
      fft_k_generated_generic((*command_queue)->kernels.fft_kernel_k_d0_d, *command_queue, &n1, &ndat, tmp, out, k, &((*command_queue)->context->cossind0));
    } else {
      n1 *= 2;
      assert((*command_queue)->context->fft_size[0]==n1);
      fft_k_generated_generic((*command_queue)->kernels.fft_kernel_k_d0_f_d, *command_queue, &n1, &ndat, tmp, out, k, &((*command_queue)->context->cossind0));
    }
}

void FC_FUNC_(fft3d_r_c2r_d_generic,FFT1D_R_C2R_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *psi,cl_mem *out,cl_mem *tmp){
    cl_uint n1 = dimensions[0]*2;
    cl_uint n2 = dimensions[1]*2;
    cl_uint n3 = dimensions[2]*2;
    if( periodic[0] ) n1 *= 2;
    if( periodic[1] ) n2 *= 2;
    if( periodic[2] ) n3 *= 2;

    cl_uint ndat = n1 * n2;

    assert((*command_queue)->context->fft_size[2]==n3);
    if( periodic[2] )
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d2_r_d, *command_queue, &n3, &ndat, psi, out, &((*command_queue)->context->cossind2));
    else {
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d2_r_f_d, *command_queue, &n3, &ndat, psi, out, &((*command_queue)->context->cossind2));
      n3 /= 2;
    }
    ndat = n1 * n3;
    assert((*command_queue)->context->fft_size[1]==n2);
    if( periodic[1] )
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d1_r_d, *command_queue, &n2, &ndat, out, tmp, &((*command_queue)->context->cossind1));
    else {
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d1_r_f_d, *command_queue, &n2, &ndat, out, tmp, &((*command_queue)->context->cossind1));
      n2 /= 2;
    }
    ndat = n2 * n3;
    assert((*command_queue)->context->fft_size[0]==n1);
    if( periodic[0] )
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d0_r_c2r_d, *command_queue, &n1, &ndat, tmp, out, &((*command_queue)->context->cossind0));
    else {
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d0_r_c2r_f_d, *command_queue, &n1, &ndat, tmp, out, &((*command_queue)->context->cossind0));
      n1 /= 2;
    }
}


//cl_program fftProgramd0;
//cl_program fftProgramd1;
//cl_program fftProgramd2;

//cl_uint fft_size[3]={0,0,0};

void create_fft_kernels(bigdft_context * context, struct bigdft_kernels * kernels){
  cl_int ciErrNum = CL_SUCCESS;
  char kernel_name[256];
  if((*context)->fft_size[0]!=0){
    sprintf(kernel_name,"fftKernel_%u_d",(*context)->fft_size[0]);
    kernels->fft_kernel_d0_d=clCreateKernel((*context)->fftProgramd0,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d0_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r2c_d",(*context)->fft_size[0]);
    kernels->fft_kernel_d0_r2c_d=clCreateKernel((*context)->fftProgramd0,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d0_r2c_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r_d",(*context)->fft_size[0]);
    kernels->fft_kernel_d0_r_d=clCreateKernel((*context)->fftProgramd0,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d0_r_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r_c2r_d",(*context)->fft_size[0]);
    kernels->fft_kernel_d0_r_c2r_d=clCreateKernel((*context)->fftProgramd0,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d0_r_c2r_d kernel!");
    sprintf(kernel_name,"fftKernel_k_%u_d",(*context)->fft_size[0]);
    kernels->fft_kernel_k_d0_d=clCreateKernel((*context)->fftProgramd0,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_k_d0_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_f_d",(*context)->fft_size[0]);
    kernels->fft_kernel_d0_f_d=clCreateKernel((*context)->fftProgramd0,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d0_f_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r2c_f_d",(*context)->fft_size[0]);
    kernels->fft_kernel_d0_r2c_f_d=clCreateKernel((*context)->fftProgramd0,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d0_r2c_f_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r_f_d",(*context)->fft_size[0]);
    kernels->fft_kernel_d0_r_f_d=clCreateKernel((*context)->fftProgramd0,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d0_r_f_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r_c2r_f_d",(*context)->fft_size[0]);
    kernels->fft_kernel_d0_r_c2r_f_d=clCreateKernel((*context)->fftProgramd0,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d0_r_c2r_f_d kernel!");
    sprintf(kernel_name,"fftKernel_k_%u_f_d",(*context)->fft_size[0]);
    kernels->fft_kernel_k_d0_f_d=clCreateKernel((*context)->fftProgramd0,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_k_d0_f_d kernel!");
  }
  if((*context)->fft_size[1]!=0){
    sprintf(kernel_name,"fftKernel_%u_d",(*context)->fft_size[1]);
    kernels->fft_kernel_d1_d=clCreateKernel((*context)->fftProgramd1,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d1_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r2c_d",(*context)->fft_size[1]);
    kernels->fft_kernel_d1_r2c_d=clCreateKernel((*context)->fftProgramd1,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d1_r2c_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r_d",(*context)->fft_size[1]);
    kernels->fft_kernel_d1_r_d=clCreateKernel((*context)->fftProgramd1,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d1_r_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r_c2r_d",(*context)->fft_size[1]);
    kernels->fft_kernel_d1_r_c2r_d=clCreateKernel((*context)->fftProgramd1,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d1_r_c2r_d kernel!");
    sprintf(kernel_name,"fftKernel_k_%u_d",(*context)->fft_size[1]);
    kernels->fft_kernel_k_d1_d=clCreateKernel((*context)->fftProgramd1,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_k_d1_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_f_d",(*context)->fft_size[1]);
    kernels->fft_kernel_d1_f_d=clCreateKernel((*context)->fftProgramd1,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d1_f_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r2c_f_d",(*context)->fft_size[1]);
    kernels->fft_kernel_d1_r2c_f_d=clCreateKernel((*context)->fftProgramd1,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d1_r2c_f_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r_f_d",(*context)->fft_size[1]);
    kernels->fft_kernel_d1_r_f_d=clCreateKernel((*context)->fftProgramd1,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d1_r_f_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r_c2r_f_d",(*context)->fft_size[1]);
    kernels->fft_kernel_d1_r_c2r_f_d=clCreateKernel((*context)->fftProgramd1,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d1_r_c2r_f_d kernel!");
    sprintf(kernel_name,"fftKernel_k_%u_f_d",(*context)->fft_size[1]);
    kernels->fft_kernel_k_d1_f_d=clCreateKernel((*context)->fftProgramd1,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_k_d1_f_d kernel!");
  }
  if((*context)->fft_size[2]!=0){
    sprintf(kernel_name,"fftKernel_%u_d",(*context)->fft_size[2]);
    kernels->fft_kernel_d2_d=clCreateKernel((*context)->fftProgramd2,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d2_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r2c_d",(*context)->fft_size[2]);
    kernels->fft_kernel_d2_r2c_d=clCreateKernel((*context)->fftProgramd2,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d2_r2c_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r_d",(*context)->fft_size[2]);
    kernels->fft_kernel_d2_r_d=clCreateKernel((*context)->fftProgramd2,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d2_r_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r_c2r_d",(*context)->fft_size[2]);
    kernels->fft_kernel_d2_r_c2r_d=clCreateKernel((*context)->fftProgramd2,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d2_r_c2r_d kernel!");
    sprintf(kernel_name,"fftKernel_k_%u_d",(*context)->fft_size[2]);
    kernels->fft_kernel_k_d2_d=clCreateKernel((*context)->fftProgramd2,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_k_d2_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_f_d",(*context)->fft_size[2]);
    kernels->fft_kernel_d2_f_d=clCreateKernel((*context)->fftProgramd2,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d2_f_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r2c_f_d",(*context)->fft_size[2]);
    kernels->fft_kernel_d2_r2c_f_d=clCreateKernel((*context)->fftProgramd2,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d2_r2c_f_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r_f_d",(*context)->fft_size[2]);
    kernels->fft_kernel_d2_r_f_d=clCreateKernel((*context)->fftProgramd2,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d2_r_f_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r_c2r_f_d",(*context)->fft_size[2]);
    kernels->fft_kernel_d2_r_c2r_f_d=clCreateKernel((*context)->fftProgramd2,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d2_r_c2r_f_d kernel!");
    sprintf(kernel_name,"fftKernel_k_%u_f_d",(*context)->fft_size[2]);
    kernels->fft_kernel_k_d2_f_d=clCreateKernel((*context)->fftProgramd2,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_k_d2_f_d kernel!");
  }
}


void build_fft_programs(bigdft_context * context){
    struct bigdft_device_infos infos;
    get_context_devices_infos(context, &infos);
    cl_int ciErrNum = CL_SUCCESS;
    cl_image_format format = { CL_RGBA, CL_UNSIGNED_INT32 };
    fft_code * c;
    if((*context)->fft_size[0]!=0){
      c = generate_fft_program((*context)->fft_size[0],&infos);
      //printf("%s\n",c->code);
      (*context)->fftProgramd0 = clCreateProgramWithSource((*context)->context, 1, (const char**) &(c->code), NULL, &ciErrNum);
      oclErrorCheck(ciErrNum,"Failed to create programd0!");
      ciErrNum = clBuildProgram((*context)->fftProgramd0, 0, NULL, "-cl-mad-enable", NULL, NULL);
      if (ciErrNum != CL_SUCCESS)
      {
          fprintf(stderr,"Error %d: Failed to build fft program d0!\n",ciErrNum);
          char cBuildLog[10240];
          clGetProgramBuildInfo((*context)->fftProgramd0, oclGetFirstDev((*context)->context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	  fprintf(stderr,"%s\n",cBuildLog);
          exit(1);
      }
      if(!use_constant_memory) {
#ifdef CL_VERSION_1_2
        if( compare_opencl_version((*context)->PLATFORM_VERSION, opencl_version_1_2) >= 0 ) {
          cl_image_desc desc;
          memset(&desc,0,sizeof(cl_image_desc));
          desc.image_type = CL_MEM_OBJECT_IMAGE1D;
          desc.image_width = (*context)->fft_size[0];
          desc.buffer = NULL;
          (*context)->cossind0 = clCreateImage((*context)->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR , &format, &desc, c->cossin, &ciErrNum);
        } else
#endif
          (*context)->cossind0 = clCreateImage2D((*context)->context,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR , &format,(*context)->fft_size[0],1,0,c->cossin,&ciErrNum);
        if (ciErrNum != CL_SUCCESS)
          fprintf(stderr,"Error %d: Failed to allocate image buffer (*context)->cossind0!\n",ciErrNum);
      }
      if(c->cossin) free(c->cossin);
      free(c->code);
      free(c);
    }
    if((*context)->fft_size[1]!=0){
      c = generate_fft_program((*context)->fft_size[1], &infos);
      //printf("%s\n",c->code);
      (*context)->fftProgramd1 = clCreateProgramWithSource((*context)->context, 1, (const char**) &(c->code), NULL, &ciErrNum);
      oclErrorCheck(ciErrNum,"Failed to create programd1!");
      ciErrNum = clBuildProgram((*context)->fftProgramd1, 0, NULL, "-cl-mad-enable", NULL, NULL);
      if (ciErrNum != CL_SUCCESS)
      {
          fprintf(stderr,"Error %d: Failed to build fft program d1!\n",ciErrNum);
          char cBuildLog[10240];
          clGetProgramBuildInfo((*context)->fftProgramd1, oclGetFirstDev((*context)->context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	  fprintf(stderr,"%s\n",cBuildLog);
          exit(1);
      }
      if(!use_constant_memory) {
#ifdef CL_VERSION_1_2
        if( compare_opencl_version((*context)->PLATFORM_VERSION, opencl_version_1_2) >= 0 ) {
          cl_image_desc desc;
          memset(&desc,0,sizeof(cl_image_desc));
          desc.image_type = CL_MEM_OBJECT_IMAGE1D;
          desc.image_width = (*context)->fft_size[1];
          desc.buffer = NULL;
          (*context)->cossind1 = clCreateImage((*context)->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR , &format, &desc, c->cossin, &ciErrNum);
        } else
#endif
          (*context)->cossind1 = clCreateImage2D((*context)->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR , &format,(*context)->fft_size[1],1,0,c->cossin,&ciErrNum);
        if (ciErrNum != CL_SUCCESS)
          fprintf(stderr,"Error %d: Failed to allocate image buffer (*context)->cossind1!\n",ciErrNum);
      }
      if(c->cossin) free(c->cossin);
      free(c->code);
      free(c);
    }
    if((*context)->fft_size[2]!=0){
      c = generate_fft_program((*context)->fft_size[2], &infos);
      //printf("%s\n",c->code);
      (*context)->fftProgramd2 = clCreateProgramWithSource((*context)->context,1,(const char**) &(c->code), NULL, &ciErrNum);
      oclErrorCheck(ciErrNum,"Failed to create programd1!");
      ciErrNum = clBuildProgram((*context)->fftProgramd2, 0, NULL, "-cl-mad-enable", NULL, NULL);
      if (ciErrNum != CL_SUCCESS)
      {
          fprintf(stderr,"Error %d: Failed to build fft program d2!\n",ciErrNum);
          char cBuildLog[10240];
          clGetProgramBuildInfo((*context)->fftProgramd2, oclGetFirstDev((*context)->context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	  fprintf(stderr,"%s\n",cBuildLog);
          exit(1);
      }
      if(!use_constant_memory) {
#ifdef CL_VERSION_1_2
        if( compare_opencl_version((*context)->PLATFORM_VERSION, opencl_version_1_2) >= 0 ) {
          cl_image_desc desc;
          memset(&desc,0,sizeof(cl_image_desc));
          desc.image_type = CL_MEM_OBJECT_IMAGE1D;
          desc.image_width = (*context)->fft_size[2];
          desc.buffer = NULL;
          (*context)->cossind2 = clCreateImage((*context)->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR , &format, &desc, c->cossin, &ciErrNum);
        } else
#endif
          (*context)->cossind2 = clCreateImage2D((*context)->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR , &format,(*context)->fft_size[2],1,0,c->cossin,&ciErrNum);
        if (ciErrNum != CL_SUCCESS)
          fprintf(stderr,"Error %d: Failed to allocate image buffer (*context)->cossind2!\n",ciErrNum);
      }
      if(c->cossin) free(c->cossin);
      free(c->code);
      free(c);
    }
}

void clean_fft_kernels(bigdft_context * context, struct bigdft_kernels * kernels){
  cl_int ciErrNum;
  if((*context)->fft_size[0]!=0){
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d0_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d0_r2c_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d0_r_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d0_r_c2r_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_k_d0_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d0_f_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d0_r2c_f_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d0_r_f_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d0_r_c2r_f_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_k_d0_f_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
  }
  if((*context)->fft_size[1]!=0){
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d1_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d1_r2c_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d1_r_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d1_r_c2r_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_k_d1_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d1_f_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d1_r2c_f_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d1_r_f_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d1_r_c2r_f_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_k_d1_f_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
  }
  if((*context)->fft_size[2]!=0){
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d2_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d2_r2c_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d2_r_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d2_r_c2r_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_k_d2_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d2_f_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d2_r2c_f_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d2_r_f_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d2_r_c2r_f_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_k_d2_f_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
  }
}

void clean_fft_programs(bigdft_context * context){
  cl_int ciErrNum;
  if((*context)->fft_size[0]!=0){
    if(!use_constant_memory){
      ciErrNum = clReleaseMemObject ((*context)->cossind0);
      oclErrorCheck(ciErrNum,"Failed to release buffer!");
    }
    ciErrNum = clReleaseProgram((*context)->fftProgramd0);
    oclErrorCheck(ciErrNum,"Failed to release program!");
  }
  if((*context)->fft_size[1]!=0){
    if(!use_constant_memory){
      ciErrNum = clReleaseMemObject ((*context)->cossind1);
      oclErrorCheck(ciErrNum,"Failed to release buffer!");
    }
    ciErrNum = clReleaseProgram((*context)->fftProgramd1);
    oclErrorCheck(ciErrNum,"Failed to release program!");
  }
  if((*context)->fft_size[2]!=0){
    if(!use_constant_memory){
      ciErrNum = clReleaseMemObject ((*context)->cossind2);
      oclErrorCheck(ciErrNum,"Failed to release buffer!");
    }
    ciErrNum = clReleaseProgram((*context)->fftProgramd2);
    oclErrorCheck(ciErrNum,"Failed to release program!");
  }
}
