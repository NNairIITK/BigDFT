#include "OpenCL_wrappers.h"

char * dgemm_program="\
//size is supposed to be 16*16\n\
#define BUFFER_SIZE 16\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
__kernel void gemmKernel_d( uint m, uint n, uint k, double alpha, __global const double *a, __global const double *b, double beta, __global double * c, __local double *tmp1, __local double *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  size_t gi = get_group_id(0);\n\
  size_t gj = get_group_id(1);\n\
  size_t si = get_local_size(0);\n\
  size_t sj = get_local_size(1);\n\
  \n\
  size_t index = 0;\n\
  size_t sumi;\n\
  double result = 0.0;\n\
  while( index < k) {\n\
    //load first matrix in tmp1\n\
    tmp1[j*(BUFFER_SIZE) + i] = (jg < m && (index + i) <  k) ? a[jg*k + index+i] : 0.0;\n\
    //load second matrix in tmp2\n\
    tmp2[j*(BUFFER_SIZE) + i] = (ig < n && (index + j) < k) ? b[(index+j)*n + ig] : 0.0;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    //compute partial result (only half bank conlicts removed)\n\
    #pragma unroll\n\
    for(sumi=0; sumi<BUFFER_SIZE; sumi++)\n\
      result += tmp1[j*(BUFFER_SIZE) + sumi] * tmp2[ sumi*(BUFFER_SIZE) + i];\n\
    index += BUFFER_SIZE;\n\
  }\n\
  if(ig < n && jg < m){\n\
    if( beta != 0.0 )\n\
      c[jg*n + ig] = alpha * result + beta * c[jg*n + ig];\n\
    else\n\
      c[jg*n + ig] = alpha * result;\n\
  }\n\
}\n\
__kernel void gemmKernel_d_opti( uint m, uint n, uint k, double alpha, __global const double *a, __global const double *b, double beta, __global double * c, __local double *tmp1, __local double *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  size_t gi = get_group_id(0);\n\
  size_t gj = get_group_id(1);\n\
  size_t si = get_local_size(0);\n\
  size_t sj = get_local_size(1);\n\
  size_t jt = (j + i)%BUFFER_SIZE;\n\
  size_t jgt = jg - j + jt;\n\
  \n\
  size_t index = 0;\n\
  size_t sumi;\n\
  double result = 0.0;\n\
  while( index < k) {\n\
    //load first matrix in tmp1\n\
    tmp1[j*(BUFFER_SIZE+1) + i] = (jg < m && (index + i) < k) ? a[jg*k + index+i] : 0.0;\n\
    //load second matrix in tmp2\n\
    tmp2[j*(BUFFER_SIZE) + i] = (ig < n && (index + j) < k) ? b[(index+j)*n + ig] : 0.0;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(sumi=0; sumi<BUFFER_SIZE; sumi++)\n\
      result += tmp1[jt*(BUFFER_SIZE+1) + sumi] * tmp2[sumi*(BUFFER_SIZE) + i];\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    index += BUFFER_SIZE;\n\
  }\n\
  if(ig < n && jgt < m){\n\
    if( beta != 0.0 )\n\
      c[jgt*n + ig] = alpha * result + beta * c[jgt*n + ig];\n\
    else\n\
      c[jgt*n + ig] = alpha * result;\n\
  }\n\
}\n\
";

char * reduction_program="\
//group_size is supposed to be 512\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
__kernel void reductionKernel_d( uint n, __global const double *x, __global double *y, __local double *tmp ) {\n\
  size_t i = get_local_id(0);\n\
  size_t g = get_group_id(0)*1024+i;\n\
  if(g<n) {\n\
    tmp[i] = x[g];\n\
  } else {\n\
    tmp[i] = 0.0;\n\
  }\n\
  if(g+512<n) {\n\
    tmp[i+512] = x[g+512];\n\
  } else {\n\
    tmp[i+512] = 0.0;\n\
  }\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  tmp[i] = tmp[i] + tmp[i+512];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<256 )\n\
    tmp[i] = tmp[i] + tmp[i+256];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<128 )\n\
    tmp[i] = tmp[i] + tmp[i+128];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<64 )\n\
    tmp[i] = tmp[i] + tmp[i+64];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<32 )\n\
    tmp[i] = tmp[i] + tmp[i+32];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<16 )\n\
    tmp[i] = tmp[i] + tmp[i+16];\n\
  if( i<8 )\n\
    tmp[i] = tmp[i] + tmp[i+8];\n\
  if( i<4 )\n\
    tmp[i] = tmp[i] + tmp[i+4];\n\
  if( i<2 )\n\
    tmp[i] = tmp[i] + tmp[i+2];\n\
  if( i==0 )\n\
    y[get_group_id(0)] = tmp[0]+tmp[1];\n\
}\n\
__kernel void reduction_dotKernel_d( uint n, __global const double *x, __global double *y, __local double *tmp ) {\n\
  size_t i = get_local_id(0);\n\
  size_t g = get_group_id(0)*1024+i;\n\
  double tt;\n\
  if(g<n) {\n\
    tt = x[g];\n\
    tmp[i] = tt*tt;\n\
  } else {\n\
    tmp[i] = 0.0;\n\
  }\n\
  if(g+512<n) {\n\
    tt = x[g+512];\n\
    tmp[i+512] = tt*tt;\n\
  } else {\n\
    tmp[i+512] = 0.0;\n\
  }\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  tmp[i] = tmp[i] + tmp[i+512];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<256 )\n\
    tmp[i] = tmp[i] + tmp[i+256];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<128 )\n\
    tmp[i] = tmp[i] + tmp[i+128];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<64 )\n\
    tmp[i] = tmp[i] + tmp[i+64];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<32 )\n\
    tmp[i] = tmp[i] + tmp[i+32];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<16 )\n\
    tmp[i] = tmp[i] + tmp[i+16];\n\
  if( i<8 )\n\
    tmp[i] = tmp[i] + tmp[i+8];\n\
  if( i<4 )\n\
    tmp[i] = tmp[i] + tmp[i+4];\n\
  if( i<2 )\n\
    tmp[i] = tmp[i] + tmp[i+2];\n\
  if( i==0 )\n\
    y[get_group_id(0)] = tmp[0]+tmp[1];\n\
}\n\
__kernel void dotKernel_d( uint n, __global const double *x, __global double *y, __global double *z, __local double *tmp ) {\n\
  size_t i = get_local_id(0);\n\
  size_t g = get_group_id(0)*1024+i;\n\
  if(g<n)\n\
    tmp[i] = x[g]*y[g];\n\
  else\n\
    tmp[i] = 0.0;\n\
  if(g+512<n)\n\
    tmp[i+512] = x[g+512]*y[g+512];\n\
  else\n\
    tmp[i+512] = 0.0;\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  tmp[i] = tmp[i] + tmp[i+512];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<256 )\n\
    tmp[i] = tmp[i] + tmp[i+256];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<128 )\n\
    tmp[i] = tmp[i] + tmp[i+128];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<64 )\n\
    tmp[i] = tmp[i] + tmp[i+64];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<32 )\n\
    tmp[i] = tmp[i] + tmp[i+32];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<16 )\n\
    tmp[i] = tmp[i] + tmp[i+16];\n\
  if( i<8 )\n\
    tmp[i] = tmp[i] + tmp[i+8];\n\
  if( i<4 )\n\
    tmp[i] = tmp[i] + tmp[i+4];\n\
  if( i<2 )\n\
    tmp[i] = tmp[i] + tmp[i+2];\n\
  if( i==0 )\n\
    z[get_group_id(0)] = tmp[0]+tmp[1];\n\
}\n\
__kernel void axpyKernel_d( uint n, double alpha, __global const double *x, __global const double *y, __global double *out) {\n\
  size_t ig = get_global_id(0);\n\
  if( ig < n)\n\
    out[ig] = y[ig] + alpha * x[ig];\n\
}\n\
__kernel void axpy_offsetKernel_d( uint n, double alpha, uint offset_x, __global const double *x, uint offset_y, __global double *y, uint offset_out, __global double *out) {\n\
  size_t ig = get_global_id(0);\n\
  if( ig < n)\n\
    out[ig+offset_out] = y[ig+offset_y] + alpha * x[ig+offset_x];\n\
}\n\
__kernel void scalKernel_d( uint n, double alpha, __global const double *x, __global double *y) {\n\
  size_t ig = get_global_id(0);\n\
  if( ig < n)\n\
    y[ig] = alpha * x[ig];\n\
}\n\
__kernel void copyKernel_d( uint n, __global const double *x, __global double *y) {\n\
  size_t ig = get_global_id(0);\n\
  if( ig < n)\n\
    y[ig] = x[ig];\n\
}\n\
__kernel void setKernel_d( uint n, const double val, __global double *x) {\n\
  size_t ig = get_global_id(0);\n\
  if( ig < n)\n\
    x[ig] = val;\n\
}\n\
\n\
";

void inline gemm_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *m, cl_uint *n, cl_uint *k, double *alpha, cl_mem *a, cl_mem *b, double *beta, cl_mem *c) {
  cl_int ciErrNum;
  size_t block_size_i=16;
  size_t block_size_j=16;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*m), (void*)m);
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*k), (void*)k);
  clSetKernelArg(kernel, i++,sizeof(*alpha), (void*)alpha);
  clSetKernelArg(kernel, i++,sizeof(*a), (void*)a);
  clSetKernelArg(kernel, i++,sizeof(*b), (void*)b);
  clSetKernelArg(kernel, i++,sizeof(*beta), (void*)beta);
  clSetKernelArg(kernel, i++,sizeof(*c), (void*)c);
  clSetKernelArg(kernel, i++,sizeof(double)*(block_size_i+1)*block_size_j, NULL);
  clSetKernelArg(kernel, i++,sizeof(double)*(block_size_i+1)*block_size_j, NULL);
  size_t localWorkSize[] = { block_size_i, block_size_j };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*m)  };
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue gemm kernel!");
}

void inline set_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, double *val, cl_mem *x) {
  cl_int ciErrNum;
  size_t block_size_i=64;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*val), (void*)val);
  clSetKernelArg(kernel, i++,sizeof(*x), (void*)x);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n) };
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue set kernel!");
}

void inline dot_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, cl_mem *x, cl_mem *y, cl_mem *out) {
  cl_int ciErrNum;
  size_t block_size_i=512;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*x), (void*)x);
  clSetKernelArg(kernel, i++,sizeof(*y), (void*)y);
  clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
  clSetKernelArg(kernel, i++,sizeof(double)*block_size_i*2, NULL);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n) };
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue dot kernel!");
}

void inline copy_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *out) {
  cl_int ciErrNum;
  size_t block_size_i=64;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*in), (void*)in);
  clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n) };
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue copy kernel!");
}

void inline scal_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, double *alpha, cl_mem *in, cl_mem *out) {
  cl_int ciErrNum;
  size_t block_size_i=64;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*alpha), (void*)alpha);
  clSetKernelArg(kernel, i++,sizeof(*in), (void*)in);
  clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n) };
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue scal kernel!");
}

void inline axpy_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, double *alpha, cl_mem *x, cl_mem *y, cl_mem *out) {
  cl_int ciErrNum;
  size_t block_size_i=64;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*alpha), (void*)alpha);
  clSetKernelArg(kernel, i++,sizeof(*x), (void*)x);
  clSetKernelArg(kernel, i++,sizeof(*y), (void*)y);
  clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n) };
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue axpy kernel!");
}

void inline axpy_offset_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, double *alpha, 
                                                                                   cl_uint *offset_x, cl_mem *x,
                                                                                   cl_uint *offset_y, cl_mem *y,
                                                                                   cl_uint *offset_out, cl_mem *out) {
  cl_int ciErrNum;
  size_t block_size_i=64;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*alpha), (void*)alpha);
  clSetKernelArg(kernel, i++,sizeof(*offset_x), (void*)offset_x);
  clSetKernelArg(kernel, i++,sizeof(*x), (void*)x);
  clSetKernelArg(kernel, i++,sizeof(*offset_y), (void*)offset_y);
  clSetKernelArg(kernel, i++,sizeof(*y), (void*)y);
  clSetKernelArg(kernel, i++,sizeof(*offset_out), (void*)offset_out);
  clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n) };
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue axpy kernel!");
}

void inline reduction_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *out) {
  cl_int ciErrNum;
  size_t block_size_i=512;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
  clSetKernelArg(kernel, i++,sizeof(*in), (void*)in);
  clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
  clSetKernelArg(kernel, i++,sizeof(double)*block_size_i*2, NULL);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i*2,*ndat)/2 };
  ciErrNum = clEnqueueNDRangeKernel(*command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue reduction kernel!");
}

void FC_FUNC_(set_d,SET_D)(cl_command_queue *command_queue, cl_uint *n, double *val, cl_mem *x){
  if(*n==0) return;
  set_generic(set_kernel_d, command_queue, n, val, x);
}

void FC_FUNC_(copy_d,COPY_D)(cl_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *out){
  if(*n==0) return;
  copy_generic(copy_kernel_d, command_queue, n, in, out);
}

void FC_FUNC_(scal_self_d,SCAL_SELF_D)(cl_command_queue *command_queue, cl_uint *n, double *alpha, cl_mem *inout){
  if(*n==0)
    return;
  scal_generic(scal_kernel_d, command_queue, n, alpha, inout, inout);
}

void FC_FUNC_(scal_d,SCAL_D)(cl_command_queue *command_queue, cl_uint *n, double *alpha, cl_mem *in, cl_mem *out){
  if(*n==0)
    return;
  scal_generic(scal_kernel_d, command_queue, n, alpha, in, out);
}

void FC_FUNC_(axpy_self_d,AXPY_SELF_D)(cl_command_queue *command_queue, cl_uint *n, double *alpha, cl_mem *in, cl_mem *inout){
  if(*n==0)
    return;
  axpy_generic(axpy_kernel_d, command_queue, n, alpha, in, inout, inout);
}

void FC_FUNC_(axpy_d,AXPY_D)(cl_command_queue *command_queue, cl_uint *n, double *alpha, cl_mem *x, cl_mem *y, cl_mem *z){
  if(*n==0)
    return;
  axpy_generic(axpy_kernel_d, command_queue, n, alpha, x, y, z);
}

void FC_FUNC_(axpy_offset_d,AXPY_OFFSET_D)(cl_command_queue *command_queue, cl_uint *n, double *alpha,
                                                                            cl_uint *offset_x, cl_mem *x,
                                                                            cl_uint *offset_y, cl_mem *y,
                                                                            cl_uint *offset_z, cl_mem *z){
  if(*n==0)
    return;
  axpy_offset_generic(axpy_offset_kernel_d, command_queue, n, alpha, offset_x, x, offset_y, y, offset_z, z);
}

void FC_FUNC_(axpy_offset_self_d,AXPY_OFFSET_SELF_D)(cl_command_queue *command_queue, cl_uint *n, double *alpha,
                                                                            cl_uint *offset_x, cl_mem *x,
                                                                            cl_uint *offset_y, cl_mem *y){
  if(*n==0)
    return;
  axpy_offset_generic(axpy_offset_kernel_d, command_queue, n, alpha, offset_x, x, offset_y, y, offset_y, y);
}

void FC_FUNC_(asum_self_d,ASUM_SELF_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work, double *out) {
  if(*ndat==0){
    *out = 0.0;
    return;
  }
  cl_uint n = *ndat;
  cl_mem *input = in;
  cl_mem *output = work;
  cl_mem *tmp;
  do {
    reduction_generic(reduction_kernel_d, command_queue, &n, input, output);
    tmp = input;
    input = output;
    output = tmp;
    n = shrRoundUp(1024,n)/1024;
  } while(n>1);
  clEnqueueReadBuffer(*command_queue, *input, CL_TRUE, 0, sizeof(double), out, 0, NULL, NULL);
}

void FC_FUNC_(nrm2sq_self_d,NRM2SQ_SELF_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work, double *out) {
  if(*ndat==0){
   *out = 0.0;
   return;
  }
  cl_uint n = *ndat;
  cl_mem *input = in;
  cl_mem *output = work;
  cl_mem *tmp;
  reduction_generic(reduction_dot_kernel_d, command_queue, &n, input, output);
  input = work;
  output = in;
  n = shrRoundUp(1024,n)/1024;
  if(n>1) {
    do {
      reduction_generic(reduction_kernel_d, command_queue, &n, input, output);
      tmp = input;
      input = output;
      output = tmp;
      n = shrRoundUp(1024,n)/1024;
    } while(n>1);
  }
  clEnqueueReadBuffer(*command_queue, *input, CL_TRUE, 0, sizeof(double), out, 0, NULL, NULL);
}

void FC_FUNC_(gemm_d,GEMM_D)(cl_command_queue *command_queue,  cl_uint *m, cl_uint *n, cl_uint *k, double *alpha, cl_mem *a, cl_mem *b, double *beta, cl_mem *c) {
  gemm_generic(gemm_kernel_d, command_queue, m, n, k, alpha, a, b, beta, c);
}

void FC_FUNC_(gemm_d_opti,GEMM_D_OPTI)(cl_command_queue *command_queue,  cl_uint *m, cl_uint *n, cl_uint *k, double *alpha, cl_mem *a, cl_mem *b, double *beta, cl_mem *c) {
  gemm_generic(gemm_kernel_d_opti, command_queue, m, n, k, alpha, a, b, beta, c);
}

void FC_FUNC_(asum_d,ASUM_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work1, cl_mem *work2, double *out) {
  if(*ndat==0){
   *out = 0.0;
   return;
  }
  cl_uint n = *ndat;
  cl_mem *input = in;
  cl_mem *output = work1;
  cl_mem *tmp;
  reduction_generic(reduction_kernel_d, command_queue, &n, input, output);
  input = work1;
  output = work2;
  n = shrRoundUp(1024,n)/1024;
  if(n>1) {
    do {
      reduction_generic(reduction_kernel_d, command_queue, &n, input, output);
      tmp = input;
      input = output;
      output = tmp;
      n = shrRoundUp(1024,n)/1024;
    } while(n>1);
  }
  clEnqueueReadBuffer(*command_queue, *input, CL_TRUE, 0, sizeof(double), out, 0, NULL, NULL);
}

void FC_FUNC_(nrm2sq_d,NRM2SQ_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work1, cl_mem *work2, double *out) {
  if(*ndat==0){
   *out = 0.0;
   return;
  }
  cl_uint n = *ndat;
  reduction_generic(reduction_dot_kernel_d, command_queue, &n, in, work1);
  cl_mem *input = work1;
  cl_mem *output = work2;
  cl_mem *tmp;
  n = shrRoundUp(1024,n)/1024;
  if(n>1) {
    do {
      reduction_generic(reduction_kernel_d, command_queue, &n, input, output);
      tmp = input;
      input = output;
      output = tmp;
      n = shrRoundUp(1024,n)/1024;
    } while(n>1);
  }
  clEnqueueReadBuffer(*command_queue, *input, CL_TRUE, 0, sizeof(double), out, 0, NULL, NULL);
}

void FC_FUNC_(dot_d,DOT_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *x, cl_mem *y, cl_mem *work1, cl_mem *work2, double *out) {
  if(*ndat==0){
   *out = 0.0;
   return;
  }
  cl_uint n = *ndat;
  dot_generic(dot_kernel_d, command_queue, &n, x, y, work1);
  cl_mem *input=work1;
  cl_mem *output=work2;
  cl_mem *tmp;
  n = shrRoundUp(1024,n)/1024;
  if(n>1) {
    do {
      reduction_generic(reduction_kernel_d, command_queue, &n, input, output);
      tmp = input;
      input = output;
      output = tmp;
      n = shrRoundUp(1024,n)/1024;
    } while(n>1);
  }
  clEnqueueReadBuffer(*command_queue, *input, CL_TRUE, 0, sizeof(double), out, 0, NULL, NULL);
}

cl_kernel reduction_kernel_d;
cl_kernel reduction_dot_kernel_d;
cl_kernel axpy_kernel_d;
cl_kernel axpy_offset_kernel_d;
cl_kernel scal_kernel_d;
cl_kernel copy_kernel_d;
cl_kernel dot_kernel_d;
cl_kernel set_kernel_d;
cl_kernel gemm_kernel_d;
cl_kernel gemm_kernel_d_opti;
cl_program reductionProgram;
cl_program dgemmProgram;

void create_reduction_kernels(){
    cl_int ciErrNum = CL_SUCCESS;
    axpy_offset_kernel_d=clCreateKernel(reductionProgram,"axpy_offsetKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    axpy_kernel_d=clCreateKernel(reductionProgram,"axpyKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    scal_kernel_d=clCreateKernel(reductionProgram,"scalKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    reduction_kernel_d=clCreateKernel(reductionProgram,"reductionKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    reduction_dot_kernel_d=clCreateKernel(reductionProgram,"reduction_dotKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    copy_kernel_d=clCreateKernel(reductionProgram,"copyKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    dot_kernel_d=clCreateKernel(reductionProgram,"dotKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    set_kernel_d=clCreateKernel(reductionProgram,"setKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    gemm_kernel_d=clCreateKernel(dgemmProgram,"gemmKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    gemm_kernel_d_opti=clCreateKernel(dgemmProgram,"gemmKernel_d_opti",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
}

void build_reduction_programs(cl_context * context){
    cl_int ciErrNum = CL_SUCCESS;
    reductionProgram = clCreateProgramWithSource(*context,1,(const char**) &reduction_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(reductionProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build reduction program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(reductionProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
        fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    dgemmProgram = clCreateProgramWithSource(*context,1,(const char**) &dgemm_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(dgemmProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build dgemm program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(dgemmProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
        fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
}

void clean_reduction_kernels(){
  cl_int ciErrNum;
  ciErrNum = clReleaseKernel(reduction_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(reduction_dot_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(axpy_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(axpy_offset_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(scal_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(copy_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(dot_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(set_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(gemm_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(gemm_kernel_d_opti);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseProgram(reductionProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
  ciErrNum = clReleaseProgram(dgemmProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
}
