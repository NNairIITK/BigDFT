//! @file
/*!
  Matrix multiplication kernels written here share a common basic design.
  Each work item computes 1 (or 2 elements in bloc design) element of
  the result matrix. Local work size is 16*16 elements. Both arrays are 
  loaded by block of 16*16 into shared memory. Buffers are big enough
  so data can be padded to 17*16. Loading both matrix allows uniform
  performance regarding transposition an conjugation.
  As computation is performed by block, at each iteration a work item is
  responsible for loading 2 elements, one from each source matrix.
  Elements out of bouns are replaced by 0.
  The first kernel is commented, other kernels are similar.
  gemmsy is also commented.
  Matrix are stored in column.
*/
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include "OpenCL_wrappers.h"
#include "Reduction_Generator.h"

char * dgemm_program="\
//size is supposed to be 16*16\n\
#ifdef cl_khr_fp64\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#elif defined (cl_amd_fp64)\n\
#pragma OPENCL EXTENSION cl_amd_fp64: enable \n\
#endif\n\
#define BUFFER_SIZE 16\n\
#define ELEM_PER_THREAD 2\n\
#define init16(array) \
array[0]=0.0;\
array[1]=0.0;\
array[2]=0.0;\
array[3]=0.0;\
array[4]=0.0;\
array[5]=0.0;\
array[6]=0.0;\
array[7]=0.0;\
array[8]=0.0;\
array[9]=0.0;\
array[10]=0.0;\
array[11]=0.0;\
array[12]=0.0;\
array[13]=0.0;\
array[14]=0.0;\
array[15]=0.0;\n\
#define axpy16(a_c,a_a,a_b,o) \
a_c[0] += a_a * a_b[0*16+o];\
a_c[1] += a_a * a_b[1*16+o];\
a_c[2] += a_a * a_b[2*16+o];\
a_c[3] += a_a * a_b[3*16+o];\
a_c[4] += a_a * a_b[4*16+o];\
a_c[5] += a_a * a_b[5*16+o];\
a_c[6] += a_a * a_b[6*16+o];\
a_c[7] += a_a * a_b[7*16+o];\
a_c[8] += a_a * a_b[8*16+o];\
a_c[9] += a_a * a_b[9*16+o];\
a_c[10] += a_a * a_b[10*16+o];\
a_c[11] += a_a * a_b[11*16+o];\
a_c[12] += a_a * a_b[12*16+o];\
a_c[13] += a_a * a_b[13*16+o];\
a_c[14] += a_a * a_b[14*16+o];\
a_c[15] += a_a * a_b[15*16+o];\n\
__kernel void gemm_volkovKernel_d( uint m, uint n, uint k, double alpha, __global const double *a, uint lda, __global const double *b, uint ldb, double beta, __global double * c, uint ldc, __local double *tmp){\n\
  double a_t;\n\
  double c_t[16];\n\
  size_t i = get_local_id(0);\n\
//  size_t j = get_local_id(1);\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1)*16;\n\
  size_t index = 0;\n\
  double result = 0.0;\n\
  bool condm = ig < m;\n\
  init16(c_t);\n\
  size_t jt = i/16+jg;\n\
  b+=jt*ldb+i%16;\n\
  a+=ig;\n\
  while( index < k - (k%16)) {\n\
    size_t it = i%16+index;\n\
    tmp[i] = (jt < n && it < k) ? b[0] : 0.0;\n\
    tmp[i+64] = (jt+4 < n && it < k) ? b[4*ldb] : 0.0;\n\
    tmp[i+128] = (jt+8 < n && it < k) ? b[8*ldb] : 0.0;\n\
    tmp[i+192] = (jt+12 < n && it < k) ? b[12*ldb] : 0.0;\n\
    b+=16;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    if(condm){\n\
    a_t = a[(index)*lda];\n\
    axpy16(c_t,a_t,tmp,0);\n\
    a_t = a[(index+1)*lda];\n\
    axpy16(c_t,a_t,tmp,1);\n\
    a_t = a[(index+2)*lda];\n\
    axpy16(c_t,a_t,tmp,2);\n\
    a_t = a[(index+3)*lda];\n\
    axpy16(c_t,a_t,tmp,3);\n\
    a_t = a[(index+4)*lda];\n\
    axpy16(c_t,a_t,tmp,4);\n\
    a_t = a[(index+5)*lda];\n\
    axpy16(c_t,a_t,tmp,5);\n\
    a_t = a[(index+6)*lda];\n\
    axpy16(c_t,a_t,tmp,6);\n\
    a_t = a[(index+7)*lda];\n\
    axpy16(c_t,a_t,tmp,7);\n\
    a_t = a[(index+8)*lda];\n\
    axpy16(c_t,a_t,tmp,8);\n\
    a_t = a[(index+9)*lda];\n\
    axpy16(c_t,a_t,tmp,9);\n\
    a_t = a[(index+10)*lda];\n\
    axpy16(c_t,a_t,tmp,10);\n\
    a_t = a[(index+11)*lda];\n\
    axpy16(c_t,a_t,tmp,11);\n\
    a_t = a[(index+12)*lda];\n\
    axpy16(c_t,a_t,tmp,12);\n\
    a_t = a[(index+13)*lda];\n\
    axpy16(c_t,a_t,tmp,13);\n\
    a_t = a[(index+14)*lda];\n\
    axpy16(c_t,a_t,tmp,14);\n\
    a_t = a[(index+15)*lda];\n\
    axpy16(c_t,a_t,tmp,15);\n\
    }\n\
    index+=16;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
  }\n\
  if( index < k ) {\n\
    size_t it = i%16+index;\n\
    tmp[i] = (jt < n && it < k) ? b[0] : 0.0;\n\
    tmp[i+64] = (jt+4 < n && it < k) ? b[4*ldb] : 0.0;\n\
    tmp[i+128] = (jt+8 < n && it < k) ? b[8*ldb] : 0.0;\n\
    tmp[i+192] = (jt+12 < n && it < k) ? b[12*ldb] : 0.0;\n\
    b+=16;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    if(condm) {\n\
    int offset=0;\n\
    switch(k-index){\n\
    case 15:\n\
    a_t = a[(index++)*lda];\n\
    axpy16(c_t,a_t,tmp,offset);\n\
    offset++;\n\
    case 14:\n\
    a_t = a[(index++)*lda];\n\
    axpy16(c_t,a_t,tmp,offset);\n\
    offset++;\n\
    case 13:\n\
    a_t = a[(index++)*lda];\n\
    axpy16(c_t,a_t,tmp,offset);\n\
    offset++;\n\
    case 12:\n\
    a_t = a[(index++)*lda];\n\
    axpy16(c_t,a_t,tmp,offset);\n\
    offset++;\n\
    case 11:\n\
    a_t = a[(index++)*lda];\n\
    axpy16(c_t,a_t,tmp,offset);\n\
    offset++;\n\
    case 10:\n\
    a_t = a[(index++)*lda];\n\
    axpy16(c_t,a_t,tmp,offset);\n\
    offset++;\n\
    case 9:\n\
    a_t = a[(index++)*lda];\n\
    axpy16(c_t,a_t,tmp,offset);\n\
    offset++;\n\
    case 8:\n\
    a_t = a[(index++)*lda];\n\
    axpy16(c_t,a_t,tmp,offset);\n\
    offset++;\n\
    case 7:\n\
    a_t = a[(index++)*lda];\n\
    axpy16(c_t,a_t,tmp,offset);\n\
    offset++;\n\
    case 6:\n\
    a_t = a[(index++)*lda];\n\
    axpy16(c_t,a_t,tmp,offset);\n\
    offset++;\n\
    case 5:\n\
    a_t = a[(index++)*lda];\n\
    axpy16(c_t,a_t,tmp,offset);\n\
    offset++;\n\
    case 4:\n\
    a_t = a[(index++)*lda];\n\
    axpy16(c_t,a_t,tmp,offset);\n\
    case 3:\n\
    a_t = a[(index++)*lda];\n\
    axpy16(c_t,a_t,tmp,offset);\n\
    offset++;\n\
    case 2:\n\
    a_t = a[(index++)*lda];\n\
    axpy16(c_t,a_t,tmp,offset);\n\
    offset++;\n\
    case 1:\n\
    a_t = a[(index)*lda];\n\
    axpy16(c_t,a_t,tmp,offset);\n\
    }\n\
    }\n\
  }\n\
  (condm && jg < n) ? c[jg*ldc + ig] = alpha * c_t[0] + beta * c[jg*ldc + ig]:0.0;\n\
  jg++;\n\
  (condm && jg < n) ? c[jg*ldc + ig] = alpha * c_t[1] + beta * c[jg*ldc + ig]:0.0;\n\
  jg++;\n\
  (condm && jg < n) ? c[jg*ldc + ig] = alpha * c_t[2] + beta * c[jg*ldc + ig]:0.0;\n\
  jg++;\n\
  (condm && jg < n) ? c[jg*ldc + ig] = alpha * c_t[3] + beta * c[jg*ldc + ig]:0.0;\n\
  jg++;\n\
  (condm && jg < n) ? c[jg*ldc + ig] = alpha * c_t[4] + beta * c[jg*ldc + ig]:0.0;\n\
  jg++;\n\
  (condm && jg < n) ? c[jg*ldc + ig] = alpha * c_t[5] + beta * c[jg*ldc + ig]:0.0;\n\
  jg++;\n\
  (condm && jg < n) ? c[jg*ldc + ig] = alpha * c_t[6] + beta * c[jg*ldc + ig]:0.0;\n\
  jg++;\n\
  (condm && jg < n) ? c[jg*ldc + ig] = alpha * c_t[7] + beta * c[jg*ldc + ig]:0.0;\n\
  jg++;\n\
  (condm && jg < n) ? c[jg*ldc + ig] = alpha * c_t[8] + beta * c[jg*ldc + ig]:0.0;\n\
  jg++;\n\
  (condm && jg < n) ? c[jg*ldc + ig] = alpha * c_t[9] + beta * c[jg*ldc + ig]:0.0;\n\
  jg++;\n\
  (condm && jg < n) ? c[jg*ldc + ig] = alpha * c_t[10] + beta * c[jg*ldc + ig]:0.0;\n\
  jg++;\n\
  (condm && jg < n) ? c[jg*ldc + ig] = alpha * c_t[11] + beta * c[jg*ldc + ig]:0.0;\n\
  jg++;\n\
  (condm && jg < n) ? c[jg*ldc + ig] = alpha * c_t[12] + beta * c[jg*ldc + ig]:0.0;\n\
  jg++;\n\
  (condm && jg < n) ? c[jg*ldc + ig] = alpha * c_t[13] + beta * c[jg*ldc + ig]:0.0;\n\
  jg++;\n\
  (condm && jg < n) ? c[jg*ldc + ig] = alpha * c_t[14] + beta * c[jg*ldc + ig]:0.0;\n\
  jg++;\n\
  (condm && jg < n) ? c[jg*ldc + ig] = alpha * c_t[15] + beta * c[jg*ldc + ig]:0.0;\n\
}\n\
__kernel void gemmKernel_d( uint m, uint n, uint k, double alpha, __global const double *a, uint lda, __global const double *b, uint ldb, double beta, __global double * c, uint ldc, __local double *tmp1, __local double *tmp2){\n\
  //get our position in the local workgroup\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  //get our position in the result matrix\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  \n\
  size_t index = 0;\n\
  double result = 0.0;\n\
  //for each block of 16*16 elements in the line of matrix a and column of matrix b\n\
  while( index < k) {\n\
    //load first matrix element in tmp1, 0.0 if out of bound.\n\
    tmp1[j*(BUFFER_SIZE) + i] = (ig < m && (index + j) < k) ? a[(index+j)*lda + ig] : 0.0;\n\
    //load second matrix element in tmp2, 0.0 if out of bound.\n\
    tmp2[j*(BUFFER_SIZE) + i] = (jg < n && (index + i) < k) ? b[(jg)*ldb + index+i] : 0.0;\n\
    //wait for buffer to be full\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    //compute partial sum, iterating over a line of tmp1 and a column of tmp2.\n\
    #pragma unroll\n\
    for(size_t sumi=0; sumi<BUFFER_SIZE; sumi++)\n\
      result += tmp1[sumi*(BUFFER_SIZE) + i] * tmp2[ j*(BUFFER_SIZE) + sumi];\n\
    //wait for buffer to be fully read.\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    index += BUFFER_SIZE;\n\
  }\n\
  //write output result if we are inside the bounds.\n\
  if(ig < m && jg < n)\n\
    c[jg*ldc + ig] = alpha * result + beta * c[jg*ldc + ig];\n\
}\n\
__kernel void gemmKernel_d_ta( uint m, uint n, uint k, double alpha, __global const double *a, uint lda, __global const double *b, uint ldb, double beta, __global double * c, uint ldc, __local double *tmp1, __local double *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  size_t igt = ig - i + j;\n\
  \n\
  size_t index = 0;\n\
  double result = 0.0;\n\
  while( index < k) {\n\
    tmp1[i*(BUFFER_SIZE+1) + j] = (igt < m && (index + i) <  k) ? a[(igt)*lda + index+i] : 0.0;\n\
    tmp2[j*(BUFFER_SIZE) + i] = (jg < n && (index + i) < k) ? b[(jg)*ldb + index+i] : 0.0;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t sumi=0; sumi<BUFFER_SIZE; sumi++)\n\
      result += tmp1[sumi*(BUFFER_SIZE+1) + i] * tmp2[ j*(BUFFER_SIZE) + sumi];\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    index += BUFFER_SIZE;\n\
  }\n\
  if(ig < m && jg < n)\n\
    c[jg*ldc + ig] = alpha * result + beta * c[jg*ldc + ig];\n\
}\n\
__kernel void gemm_blockKernel_d( uint m, uint n, uint k, double alpha, __global const double *a, uint lda, __global const double *b, uint ldb, double beta, __global double * c, uint ldc, __local double *tmp1, __local double *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1)*ELEM_PER_THREAD;\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1)*ELEM_PER_THREAD;\n\
  \n\
  size_t index = 0;\n\
  double result_1 = 0.0;\n\
  double result_2 = 0.0;\n\
  a += ig+j*lda;\n\
  b+= i+jg*ldb;\n\
  while( index < k) {\n\
    __local double *tmp2_t = tmp2 + j*BUFFER_SIZE;\n\
    tmp1[j*(BUFFER_SIZE) + i] = (ig < m && (index + j) < k) ? *a : 0.0;\n\
    tmp1[(j+1)*(BUFFER_SIZE) + i] = (ig < m && (index + j+1) < k) ? a[lda] : 0.0;\n\
    tmp2_t[i] = (jg < n && (index + i) < k) ? *b : 0.0;\n\
    tmp2_t[BUFFER_SIZE + i] = ((jg+1) < n && (index + i) < k) ? b[ldb] : 0.0;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    __local double *tmp1_t = tmp1 + i;\n\
    for(size_t sumi=0; sumi<BUFFER_SIZE; sumi++){\n\
      result_1 += *tmp1_t * *tmp2_t;\n\
      result_2 += *tmp1_t * tmp2_t[BUFFER_SIZE];\n\
      tmp1_t+=BUFFER_SIZE;\n\
      tmp2_t++;\n\
    }\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    index += BUFFER_SIZE;\n\
    b += BUFFER_SIZE;\n\
    a += BUFFER_SIZE*lda;\n\
  }\n\
  c += jg*ldc + ig;\n\
  if(ig < m && jg < n){\n\
    *c = alpha * result_1 + beta * *c;\n\
    c += ldc;\n\
    *c = alpha * result_2 + beta * *c;\n\
  }\n\
}\n\
//blocks are numbered in the following way :\n\
// 0 1 3 6 10\n\
//   2 4 7 11\n\
//     5 8 12\n\
//       9 13\n\
//         14\n\
//position in the matrix is obtained via solving ig * ( ig + 1) / 2 = index\n\
//and rounding to the apropriate integer.\n\
//non numbered blocks are obtained by transposing numbered blocks.\n\
__kernel void gemmsyKernel_d( uint m, uint k, double alpha, __global const double *a, uint lda, __global const double *b, uint ldb, double beta, __global double * c, uint ldc, __local double *tmp1, __local double *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  size_t ig = get_group_id(0);\n\
  ig = ceil(-0.500001+0.5*sqrt(9.0+8*ig))-1;\n\
  size_t jg = get_group_id(0) - ig * (ig + 1) / 2;\n\
  ig = ig * BUFFER_SIZE + i;\n\
  jg = jg * BUFFER_SIZE + j;\n\
  \n\
  size_t index = 0;\n\
  double result = 0.0;\n\
  while( index < k) {\n\
    tmp1[j*(BUFFER_SIZE) + i] = (ig < m && (index + j) < k) ? a[(index+j)*lda + ig] : 0.0;\n\
    tmp2[j*(BUFFER_SIZE) + i] = (jg < m && (index + i) < k) ? b[(jg)*ldb + index+i] : 0.0;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t sumi=0; sumi<BUFFER_SIZE; sumi++)\n\
      result += tmp1[sumi*(BUFFER_SIZE) + i] * tmp2[ j*(BUFFER_SIZE) + sumi];\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    index += BUFFER_SIZE;\n\
  }\n\
  result = alpha * result + beta * c[jg*ldc + ig];\n\
  tmp1[j*(BUFFER_SIZE+1) + i] = result;\n\
  size_t igt = ig - i + j;\n\
  size_t jgt = jg - j + i;\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if(ig < m && jg < m)\n\
    c[jg*ldc + ig] = result;\n\
  if(igt < m && jgt < m)\n\
    c[igt*ldc + jgt] = tmp1[i*(BUFFER_SIZE+1) + j];\n\
}\n\
__kernel void gemmKernel_d_tb( uint m, uint n, uint k, double alpha, __global const double *a, uint lda, __global const double *b, uint ldb, double beta, __global double * c, uint ldc, __local double *tmp1, __local double *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  size_t jgt = jg - j + i;\n\
  \n\
  size_t index = 0;\n\
  double result = 0.0;\n\
  while( index < k) {\n\
    tmp1[j*(BUFFER_SIZE) + i] = (ig < m && (index + j) <  k) ? a[(index+j)*lda + ig] : 0.0;\n\
    tmp2[i*(BUFFER_SIZE+1) + j] = (jgt < n && (index + j) < k) ? b[(index+j)*ldb + jgt] : 0.0;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t sumi=0; sumi<BUFFER_SIZE; sumi++)\n\
      result += tmp1[sumi*(BUFFER_SIZE) + i] * tmp2[ j*(BUFFER_SIZE+1) + sumi];\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    index += BUFFER_SIZE;\n\
  }\n\
  if(ig < m && jg < n)\n\
    c[jg*ldc + ig] = alpha * result + beta * c[jg*ldc + ig];\n\
}\n\
__kernel void gemmKernel_d_tatb( uint m, uint n, uint k, double alpha, __global const double *a, uint lda, __global const double *b, uint ldb, double beta, __global double * c, uint ldc, __local double *tmp1, __local double *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  size_t jgt = jg - j + i;\n\
  size_t igt = ig - i + j;\n\
  \n\
  size_t index = 0;\n\
  double result = 0.0;\n\
  while( index < k) {\n\
    tmp1[i*(BUFFER_SIZE+1) + j] = (igt < m && (index + i) <  k) ? a[(igt)*lda + index+i] : 0.0;\n\
    tmp2[i*(BUFFER_SIZE+1) + j] = (jgt < n && (index + j) < k) ? b[(index+j)*ldb + jgt] : 0.0;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t sumi=0; sumi<BUFFER_SIZE; sumi++)\n\
      result += tmp1[sumi*(BUFFER_SIZE+1) + i] * tmp2[ j*(BUFFER_SIZE+1) + sumi];\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    index += BUFFER_SIZE;\n\
  }\n\
  if(ig < m && jg < n)\n\
    c[jg*ldc + ig] = alpha * result + beta * c[jg*ldc + ig];\n\
}\n\
__kernel __attribute__((reqd_work_group_size(16, 16, 1))) __attribute__((vec_type_hint(double2))) void gemmKernel_z( const uint m, const uint n, const uint k, const double2 alpha, __global const double2 *a, const uint lda, __global const double2 *b, const uint ldb, const double2 beta, __global double2 * c, const uint ldc, __local double2 *tmp1, __local double2 *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  \n\
  double2 result __attribute__ ((aligned (16)));\n\
  double2 a_t __attribute__ ((aligned (16)));\n\
  double2 b_t __attribute__ ((aligned (16)));\n\
  result = (double2)(0.0, 0.0);\n\
  size_t index = 0;\n\
  tmp2 += j*(BUFFER_SIZE);\n\
  tmp1 += i;\n\
  a += j*lda + ig;\n\
  b += jg*ldb + i;\n\
  c += jg*ldc + ig;\n\
  bool condm = ig < m;\n\
  bool condn = jg < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    tmp1[j*(BUFFER_SIZE)] = (condm && (index + j) < k) ? a[index*lda] : (double2)(0.0, 0.0);\n\
    tmp2[i] = (condn && (index + i) < k) ? b[index] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*BUFFER_SIZE];\n\
      b_t = tmp2[k];\n\
      result.x += a_t.x * b_t.x;\n\
      result.x -= a_t.y * b_t.y;\n\
      result.y += a_t.x * b_t.y;\n\
      result.y += a_t.y * b_t.x;\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x -= alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x -= beta.y*c_t.y;\n\
    final_result.y += beta.x*c_t.y;\n\
    final_result.y += beta.y*c_t.x;\n\
    c[0] = final_result;\n\
  }\n\
}\n\
__kernel __attribute__((reqd_work_group_size(16, 16, 1))) __attribute__((vec_type_hint(double2))) void gemmKernel_z_ta( const uint m, const uint n, const uint k, const double2 alpha, __global const double2 *a, const uint lda, __global const double2 *b, const uint ldb, const double2 beta, __global double2 * c, const uint ldc, __local double2 *tmp1, __local double2 *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  size_t igt = ig - i + j;\n\
  \n\
  double2 result __attribute__ ((aligned (16)));\n\
  double2 a_t __attribute__ ((aligned (16)));\n\
  double2 b_t __attribute__ ((aligned (16)));\n\
  result = (double2)(0.0, 0.0);\n\
  size_t index = 0;\n\
  tmp2 += j*(BUFFER_SIZE);\n\
  a += igt*lda + i;\n\
  b += jg*ldb + i;\n\
  c += jg*ldc + ig;\n\
  bool condm = igt < m;\n\
  bool condn = jg < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    tmp1[i*(BUFFER_SIZE+1)+j] = (condm && (index + i) < k) ? a[index] : (double2)(0.0, 0.0);\n\
    tmp2[i] = (condn && (index + i) < k) ? b[index] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*(BUFFER_SIZE+1) + i];\n\
      b_t = tmp2[k];\n\
      result.x += a_t.x * b_t.x;\n\
      result.x -= a_t.y * b_t.y;\n\
      result.y += a_t.x * b_t.y;\n\
      result.y += a_t.y * b_t.x;\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x -= alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x -= beta.y*c_t.y;\n\
    final_result.y += beta.x*c_t.y;\n\
    final_result.y += beta.y*c_t.x;\n\
    c[0] = final_result;\n\
  }\n\
}\n\
__kernel __attribute__((reqd_work_group_size(16, 16, 1))) __attribute__((vec_type_hint(double2))) void gemmKernel_z_tb( const uint m, const uint n, const uint k, const double2 alpha, __global const double2 *a, const uint lda, __global const double2 *b, const uint ldb, const double2 beta, __global double2 * c, const uint ldc, __local double2 *tmp1, __local double2 *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  size_t jgt = jg - j + i;\n\
  \n\
  double2 result __attribute__ ((aligned (16)));\n\
  double2 a_t __attribute__ ((aligned (16)));\n\
  double2 b_t __attribute__ ((aligned (16)));\n\
  result = (double2)(0.0, 0.0);\n\
  size_t index = 0;\n\
  tmp1 += i;\n\
  a += j*lda + ig;\n\
  b += j*ldb + jgt;\n\
  c += jg*ldc + ig;\n\
  bool condm = ig < m;\n\
  bool condn = jgt < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    tmp1[j*(BUFFER_SIZE)] = (condm && (index + j) < k) ? a[index*lda] : (double2)(0.0, 0.0);\n\
    tmp2[i*(BUFFER_SIZE+1) + j] = (condn && (index + j) < k) ? b[index*ldb] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*BUFFER_SIZE];\n\
      b_t = tmp2[j*(BUFFER_SIZE+1) + k];\n\
      result.x += a_t.x * b_t.x;\n\
      result.x -= a_t.y * b_t.y;\n\
      result.y += a_t.x * b_t.y;\n\
      result.y += a_t.y * b_t.x;\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x -= alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x -= beta.y*c_t.y;\n\
    final_result.y += beta.x*c_t.y;\n\
    final_result.y += beta.y*c_t.x;\n\
    c[0] = final_result;\n\
  }\n\
}\n\
__kernel __attribute__((reqd_work_group_size(16, 16, 1))) __attribute__((vec_type_hint(double2))) void gemmKernel_z_tatb( const uint m, const uint n, const uint k, const double2 alpha, __global const double2 *a, const uint lda, __global const double2 *b, const uint ldb, const double2 beta, __global double2 * c, const uint ldc, __local double2 *tmp1, __local double2 *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  size_t jgt = jg - j + i;\n\
  size_t igt = ig - i + j;\n\
  \n\
  double2 result __attribute__ ((aligned (16)));\n\
  double2 a_t __attribute__ ((aligned (16)));\n\
  double2 b_t __attribute__ ((aligned (16)));\n\
  result = (double2)(0.0, 0.0);\n\
  size_t index = 0;\n\
  a += igt*lda + i;\n\
  b += j*ldb + jgt;\n\
  c += jg*ldc + ig;\n\
  bool condm = igt < m;\n\
  bool condn = jgt < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    tmp1[i*(BUFFER_SIZE+1) + j] = (condm && (index + i) < k) ? a[index] : (double2)(0.0, 0.0);\n\
    tmp2[i*(BUFFER_SIZE+1) + j] = (condn && (index + j) < k) ? b[index*ldb] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*(BUFFER_SIZE+1) + i];\n\
      b_t = tmp2[j*(BUFFER_SIZE+1) + k];\n\
      result.x += a_t.x * b_t.x;\n\
      result.x -= a_t.y * b_t.y;\n\
      result.y += a_t.x * b_t.y;\n\
      result.y += a_t.y * b_t.x;\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x -= alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x -= beta.y*c_t.y;\n\
    final_result.y += beta.x*c_t.y;\n\
    final_result.y += beta.y*c_t.x;\n\
    c[0] = final_result;\n\
  }\n\
}\n\
__kernel __attribute__((reqd_work_group_size(16, 16, 1))) __attribute__((vec_type_hint(double2))) void gemmKernel_z_ca( const uint m, const uint n, const uint k, const double2 alpha, __global const double2 *a, const uint lda, __global const double2 *b, const uint ldb, const double2 beta, __global double2 * c, const uint ldc, __local double2 *tmp1, __local double2 *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  size_t igt = ig - i + j;\n\
  \n\
  double2 result __attribute__ ((aligned (16)));\n\
  double2 a_t __attribute__ ((aligned (16)));\n\
  double2 b_t __attribute__ ((aligned (16)));\n\
  result = (double2)(0.0, 0.0);\n\
  size_t index = 0;\n\
  tmp2 += j*(BUFFER_SIZE);\n\
  a += igt*lda + i;\n\
  b += jg*ldb + i;\n\
  c += jg*ldc + ig;\n\
  bool condm = igt < m;\n\
  bool condn = jg < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    tmp1[i*(BUFFER_SIZE+1)+j] = (condm && (index + i) < k) ? a[index] : (double2)(0.0, 0.0);\n\
    tmp2[i] = (condn && (index + i) < k) ? b[index] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*(BUFFER_SIZE+1) + i];\n\
      b_t = tmp2[k];\n\
      result.x += a_t.x * b_t.x;\n\
      result.x += a_t.y * b_t.y;\n\
      result.y += a_t.x * b_t.y;\n\
      result.y -= a_t.y * b_t.x;\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x -= alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x -= beta.y*c_t.y;\n\
    final_result.y += beta.x*c_t.y;\n\
    final_result.y += beta.y*c_t.x;\n\
    c[0] = final_result;\n\
  }\n\
}\n\
__kernel __attribute__((reqd_work_group_size(16, 16, 1))) __attribute__((vec_type_hint(double2))) void gemmKernel_z_cb( const uint m, const uint n, const uint k, const double2 alpha, __global const double2 *a, const uint lda, __global const double2 *b, const uint ldb, const double2 beta, __global double2 * c, const uint ldc, __local double2 *tmp1, __local double2 *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  size_t jgt = jg - j + i;\n\
  \n\
  double2 result __attribute__ ((aligned (16)));\n\
  double2 a_t __attribute__ ((aligned (16)));\n\
  double2 b_t __attribute__ ((aligned (16)));\n\
  result = (double2)(0.0, 0.0);\n\
  size_t index = 0;\n\
  tmp1 += i;\n\
  a += j*lda + ig;\n\
  b += j*ldb + jgt;\n\
  c += jg*ldc + ig;\n\
  bool condm = ig < m;\n\
  bool condn = jgt < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    tmp1[j*(BUFFER_SIZE)] = (condm && (index + j) < k) ? a[index*lda] : (double2)(0.0, 0.0);\n\
    tmp2[i*(BUFFER_SIZE+1) + j] = (condn && (index + j) < k) ? (double2)(1,-1) * b[index*ldb] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*BUFFER_SIZE];\n\
      b_t = tmp2[j*(BUFFER_SIZE+1) + k];\n\
      result.x += a_t.x * b_t.x;\n\
      result.x -= a_t.y * b_t.y;\n\
      result.y += a_t.x * b_t.y;\n\
      result.y += a_t.y * b_t.x;\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x -= alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x -= beta.y*c_t.y;\n\
    final_result.y += beta.x*c_t.y;\n\
    final_result.y += beta.y*c_t.x;\n\
    c[0] = final_result;\n\
  }\n\
}\n\
__kernel __attribute__((reqd_work_group_size(16, 16, 1))) __attribute__((vec_type_hint(double2))) void gemmKernel_z_cacb( const uint m, const uint n, const uint k, const double2 alpha, __global const double2 *a, const uint lda, __global const double2 *b, const uint ldb, const double2 beta, __global double2 * c, const uint ldc, __local double2 *tmp1, __local double2 *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  size_t jgt = jg - j + i;\n\
  size_t igt = ig - i + j;\n\
  \n\
  double2 result __attribute__ ((aligned (16)));\n\
  double2 a_t __attribute__ ((aligned (16)));\n\
  double2 b_t __attribute__ ((aligned (16)));\n\
  result = (double2)(0.0, 0.0);\n\
  size_t index = 0;\n\
  a += igt*lda + i;\n\
  b += j*ldb + jgt;\n\
  c += jg*ldc + ig;\n\
  bool condm = igt < m;\n\
  bool condn = jgt < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    tmp1[i*(BUFFER_SIZE+1) + j] = (condm && (index + i) < k) ? a[index] : (double2)(0.0, 0.0);\n\
    tmp2[i*(BUFFER_SIZE+1) + j] = (condn && (index + j) < k) ? (double2)(1,-1) * b[index*ldb] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*(BUFFER_SIZE+1) + i];\n\
      b_t = tmp2[j*(BUFFER_SIZE+1) + k];\n\
      result.x += a_t.x * b_t.x;\n\
      result.x += a_t.y * b_t.y;\n\
      result.y += a_t.x * b_t.y;\n\
      result.y -= a_t.y * b_t.x;\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x -= alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x -= beta.y*c_t.y;\n\
    final_result.y += beta.x*c_t.y;\n\
    final_result.y += beta.y*c_t.x;\n\
    c[0] = final_result;\n\
  }\n\
}\n\
__kernel __attribute__((reqd_work_group_size(16, 16, 1))) __attribute__((vec_type_hint(double2))) void gemmKernel_z_catb( const uint m, const uint n, const uint k, const double2 alpha, __global const double2 *a, const uint lda, __global const double2 *b, const uint ldb, const double2 beta, __global double2 * c, const uint ldc, __local double2 *tmp1, __local double2 *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  size_t jgt = jg - j + i;\n\
  size_t igt = ig - i + j;\n\
  \n\
  double2 result __attribute__ ((aligned (16)));\n\
  double2 a_t __attribute__ ((aligned (16)));\n\
  double2 b_t __attribute__ ((aligned (16)));\n\
  result = (double2)(0.0, 0.0);\n\
  size_t index = 0;\n\
  a += igt*lda + i;\n\
  b += j*ldb + jgt;\n\
  c += jg*ldc + ig;\n\
  bool condm = igt < m;\n\
  bool condn = jgt < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    tmp1[i*(BUFFER_SIZE+1) + j] = (condm && (index + i) < k) ? a[index] : (double2)(0.0, 0.0);\n\
    tmp2[i*(BUFFER_SIZE+1) + j] = (condn && (index + j) < k) ? b[index*ldb] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*(BUFFER_SIZE+1) + i];\n\
      b_t = tmp2[j*(BUFFER_SIZE+1) + k];\n\
      result.x += a_t.x * b_t.x;\n\
      result.x += a_t.y * b_t.y;\n\
      result.y += a_t.x * b_t.y;\n\
      result.y -= a_t.y * b_t.x;\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x -= alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x -= beta.y*c_t.y;\n\
    final_result.y += beta.x*c_t.y;\n\
    final_result.y += beta.y*c_t.x;\n\
    c[0] = final_result;\n\
  }\n\
}\n\
__kernel __attribute__((reqd_work_group_size(16, 16, 1))) __attribute__((vec_type_hint(double2))) void gemmKernel_z_tacb( const uint m, const uint n, const uint k, const double2 alpha, __global const double2 *a, const uint lda, __global const double2 *b, const uint ldb, const double2 beta, __global double2 * c, const uint ldc, __local double2 *tmp1, __local double2 *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  size_t jgt = jg - j + i;\n\
  size_t igt = ig - i + j;\n\
  \n\
  double2 result __attribute__ ((aligned (16)));\n\
  double2 a_t __attribute__ ((aligned (16)));\n\
  double2 b_t __attribute__ ((aligned (16)));\n\
  result = (double2)(0.0, 0.0);\n\
  size_t index = 0;\n\
  a += igt*lda + i;\n\
  b += j*ldb + jgt;\n\
  c += jg*ldc + ig;\n\
  bool condm = igt < m;\n\
  bool condn = jgt < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    tmp1[i*(BUFFER_SIZE+1) + j] = (condm && (index + i) < k) ? a[index] : (double2)(0.0, 0.0);\n\
    tmp2[i*(BUFFER_SIZE+1) + j] = (condn && (index + j) < k) ? (double2)(1,-1) * b[index*ldb] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*(BUFFER_SIZE+1) + i];\n\
      b_t = tmp2[j*(BUFFER_SIZE+1) + k];\n\
      result.x += a_t.x * b_t.x;\n\
      result.x -= a_t.y * b_t.y;\n\
      result.y += a_t.x * b_t.y;\n\
      result.y += a_t.y * b_t.x;\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x -= alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x -= beta.y*c_t.y;\n\
    final_result.y += beta.x*c_t.y;\n\
    final_result.y += beta.y*c_t.x;\n\
    c[0] = final_result;\n\
  }\n\
}\n\
";



void inline zgemm_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *m, cl_uint *n, cl_uint *k, cl_double2 *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double2 *beta, cl_mem *c, cl_uint *ldc) {
  cl_int ciErrNum;
  size_t block_size_i=16;
  size_t block_size_j=16;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*m), (void*)m);
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*k), (void*)k);
  clSetKernelArg(kernel, i++,sizeof(*alpha), (void*)alpha);
  clSetKernelArg(kernel, i++,sizeof(*a), (void*)a);
  clSetKernelArg(kernel, i++,sizeof(*lda), (void*)lda);
  clSetKernelArg(kernel, i++,sizeof(*b), (void*)b);
  clSetKernelArg(kernel, i++,sizeof(*ldb), (void*)ldb);
  clSetKernelArg(kernel, i++,sizeof(*beta), (void*)beta);
  clSetKernelArg(kernel, i++,sizeof(*c), (void*)c);
  clSetKernelArg(kernel, i++,sizeof(*ldc), (void*)ldc);
  clSetKernelArg(kernel, i++,sizeof(cl_double2)*(block_size_i+1)*block_size_j, NULL);
  clSetKernelArg(kernel, i++,sizeof(cl_double2)*(block_size_i+1)*block_size_j, NULL);
  size_t localWorkSize[] = { block_size_i, block_size_j };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*m), shrRoundUp(block_size_j,*n)  };
  ciErrNum = clEnqueueNDRangeKernel(command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue zgemm kernel!");
}
void inline gemm_block_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *m, cl_uint *n, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc) {
  cl_int ciErrNum;
  int ELEM_PER_THREAD=2;
  size_t block_size_i=16;
  size_t block_size_j=16/ELEM_PER_THREAD;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*m), (void*)m);
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*k), (void*)k);
  clSetKernelArg(kernel, i++,sizeof(*alpha), (void*)alpha);
  clSetKernelArg(kernel, i++,sizeof(*a), (void*)a);
  clSetKernelArg(kernel, i++,sizeof(*lda), (void*)lda);
  clSetKernelArg(kernel, i++,sizeof(*b), (void*)b);
  clSetKernelArg(kernel, i++,sizeof(*ldb), (void*)ldb);
  clSetKernelArg(kernel, i++,sizeof(*beta), (void*)beta);
  clSetKernelArg(kernel, i++,sizeof(*c), (void*)c);
  clSetKernelArg(kernel, i++,sizeof(*ldc), (void*)ldc);
  clSetKernelArg(kernel, i++,sizeof(cl_double)*(block_size_i+1)*block_size_j*ELEM_PER_THREAD, NULL);
  clSetKernelArg(kernel, i++,sizeof(cl_double)*(block_size_i+1)*block_size_j*ELEM_PER_THREAD, NULL);
  size_t localWorkSize[] = { block_size_i, block_size_j };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*m), shrRoundUp(block_size_j*ELEM_PER_THREAD,*n)*block_size_j/16  };
  ciErrNum = clEnqueueNDRangeKernel(command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue gemm_block kernel!");
}
void inline gemm_volkov_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *m, cl_uint *n, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc) {
  cl_int ciErrNum;
  int ELEM_PER_THREAD=16;
  size_t block_size_i=64;
  size_t block_size_j=16/ELEM_PER_THREAD;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*m), (void*)m);
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*k), (void*)k);
  clSetKernelArg(kernel, i++,sizeof(*alpha), (void*)alpha);
  clSetKernelArg(kernel, i++,sizeof(*a), (void*)a);
  clSetKernelArg(kernel, i++,sizeof(*lda), (void*)lda);
  clSetKernelArg(kernel, i++,sizeof(*b), (void*)b);
  clSetKernelArg(kernel, i++,sizeof(*ldb), (void*)ldb);
  clSetKernelArg(kernel, i++,sizeof(*beta), (void*)beta);
  clSetKernelArg(kernel, i++,sizeof(*c), (void*)c);
  clSetKernelArg(kernel, i++,sizeof(*ldc), (void*)ldc);
  clSetKernelArg(kernel, i++,sizeof(cl_double)*(block_size_i/4+1)*block_size_j*ELEM_PER_THREAD, NULL);
  size_t localWorkSize[] = { block_size_i, block_size_j };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*m), shrRoundUp(block_size_j*ELEM_PER_THREAD,*n)*block_size_j/16  };
  ciErrNum = clEnqueueNDRangeKernel(command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue gemm_volkov kernel!");
}
void inline gemm_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *m, cl_uint *n, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc) {
  cl_int ciErrNum;
  size_t block_size_i=16;
  size_t block_size_j=16;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*m), (void*)m);
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*k), (void*)k);
  clSetKernelArg(kernel, i++,sizeof(*alpha), (void*)alpha);
  clSetKernelArg(kernel, i++,sizeof(*a), (void*)a);
  clSetKernelArg(kernel, i++,sizeof(*lda), (void*)lda);
  clSetKernelArg(kernel, i++,sizeof(*b), (void*)b);
  clSetKernelArg(kernel, i++,sizeof(*ldb), (void*)ldb);
  clSetKernelArg(kernel, i++,sizeof(*beta), (void*)beta);
  clSetKernelArg(kernel, i++,sizeof(*c), (void*)c);
  clSetKernelArg(kernel, i++,sizeof(*ldc), (void*)ldc);
  clSetKernelArg(kernel, i++,sizeof(cl_double)*(block_size_i+1)*block_size_j, NULL);
  clSetKernelArg(kernel, i++,sizeof(cl_double)*(block_size_i+1)*block_size_j, NULL);
  size_t localWorkSize[] = { block_size_i, block_size_j };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*m), shrRoundUp(block_size_j,*n)  };
  ciErrNum = clEnqueueNDRangeKernel(command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue gemm kernel!");
}
void inline gemmsy_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *m, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc) {
  cl_int ciErrNum;
  size_t block_size_i=16;
  size_t block_size_j=16;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*m), (void*)m);
  clSetKernelArg(kernel, i++,sizeof(*k), (void*)k);
  clSetKernelArg(kernel, i++,sizeof(*alpha), (void*)alpha);
  clSetKernelArg(kernel, i++,sizeof(*a), (void*)a);
  clSetKernelArg(kernel, i++,sizeof(*lda), (void*)lda);
  clSetKernelArg(kernel, i++,sizeof(*b), (void*)b);
  clSetKernelArg(kernel, i++,sizeof(*ldb), (void*)ldb);
  clSetKernelArg(kernel, i++,sizeof(*beta), (void*)beta);
  clSetKernelArg(kernel, i++,sizeof(*c), (void*)c);
  clSetKernelArg(kernel, i++,sizeof(*ldc), (void*)ldc);
  clSetKernelArg(kernel, i++,sizeof(cl_double)*(block_size_i+1)*block_size_j, NULL);
  clSetKernelArg(kernel, i++,sizeof(cl_double)*(block_size_i+1)*block_size_j, NULL);
  size_t localWorkSize[] = { block_size_i, block_size_j };
  size_t group_number = shrRoundUp(block_size_i,*m)/block_size_i;
  group_number = group_number * ( group_number + 1 ) / 2;
  size_t globalWorkSize[] ={ group_number * block_size_i, block_size_j  };
  ciErrNum = clEnqueueNDRangeKernel(command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue gemmsy kernel!");
}
void inline set_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n, cl_double *val, cl_mem *x) {
  cl_int ciErrNum;
  size_t block_size_i=64;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*val), (void*)val);
  clSetKernelArg(kernel, i++,sizeof(*x), (void*)x);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n) };
  ciErrNum = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue set kernel!");
}

void inline dot_generic(cl_kernel kernel, bigdft_command_queue command_queue, cl_uint *n, cl_mem *x, cl_mem *y, cl_mem *out) {
  cl_int ciErrNum;
  size_t block_size_i=command_queue->device_infos.MAX_WORK_GROUP_SIZE;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*x), (void*)x);
  clSetKernelArg(kernel, i++,sizeof(*y), (void*)y);
  clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
  clSetKernelArg(kernel, i++,sizeof(double)*block_size_i*2, NULL);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n) };
  ciErrNum = clEnqueueNDRangeKernel(command_queue->command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue dot kernel!");
}

void inline copy_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n, cl_mem *in, cl_mem *out) {
  cl_int ciErrNum;
  size_t block_size_i=64;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*in), (void*)in);
  clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n) };
  ciErrNum = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue copy kernel!");
}

void inline scal_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n, cl_double *alpha, cl_mem *in, cl_mem *out) {
  cl_int ciErrNum;
  size_t block_size_i=64;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*alpha), (void*)alpha);
  clSetKernelArg(kernel, i++,sizeof(*in), (void*)in);
  clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n) };
  ciErrNum = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue scal kernel!");
}

void inline axpy_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n, cl_double *alpha, cl_mem *x, cl_mem *y, cl_mem *out) {
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
  ciErrNum = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue axpy kernel!");
}

void inline axpy_offset_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n, cl_double *alpha, 
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
  ciErrNum = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue axpy kernel!");
}

void inline reduction_generic(cl_kernel kernel, bigdft_command_queue command_queue, cl_uint *ndat, cl_mem *in, cl_mem *out) {
  cl_int ciErrNum;
  size_t block_size_i=command_queue->device_infos.MAX_WORK_GROUP_SIZE;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
  clSetKernelArg(kernel, i++,sizeof(*in), (void*)in);
  clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
  clSetKernelArg(kernel, i++,sizeof(cl_double)*block_size_i*2, NULL);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i*2,*ndat)/2 };
  ciErrNum = clEnqueueNDRangeKernel(command_queue->command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue reduction kernel!");
}

cl_program reductionProgram;
cl_program dgemmProgram;

void FC_FUNC_(set_d,SET_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_double *val, cl_mem *x){
  if(*n==0) return;
  set_generic((*command_queue)->kernels.set_kernel_d, (*command_queue)->command_queue, n, val, x);
}

void FC_FUNC_(copy_d,COPY_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *out){
  if(*n==0) return;
  copy_generic((*command_queue)->kernels.copy_kernel_d, (*command_queue)->command_queue, n, in, out);
}

void FC_FUNC_(scal_self_d,SCAL_SELF_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *inout){
  if(*n==0)
    return;
  scal_generic((*command_queue)->kernels.scal_kernel_d, (*command_queue)->command_queue, n, alpha, inout, inout);
}

void FC_FUNC_(scal_d,SCAL_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *in, cl_mem *out){
  if(*n==0)
    return;
  scal_generic((*command_queue)->kernels.scal_kernel_d, (*command_queue)->command_queue, n, alpha, in, out);
}

void FC_FUNC_(axpy_self_d,AXPY_SELF_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *in, cl_mem *inout){
  if(*n==0)
    return;
  axpy_generic((*command_queue)->kernels.axpy_kernel_d, (*command_queue)->command_queue, n, alpha, in, inout, inout);
}

void FC_FUNC_(axpy_d,AXPY_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *x, cl_mem *y, cl_mem *z){
  if(*n==0)
    return;
  axpy_generic((*command_queue)->kernels.axpy_kernel_d, (*command_queue)->command_queue, n, alpha, x, y, z);
}

void FC_FUNC_(axpy_offset_d,AXPY_OFFSET_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_double *alpha,
                                                                            cl_uint *offset_x, cl_mem *x,
                                                                            cl_uint *offset_y, cl_mem *y,
                                                                            cl_uint *offset_z, cl_mem *z){
  if(*n==0)
    return;
  axpy_offset_generic((*command_queue)->kernels.axpy_offset_kernel_d, (*command_queue)->command_queue, n, alpha, offset_x, x, offset_y, y, offset_z, z);
}

void FC_FUNC_(axpy_offset_self_d,AXPY_OFFSET_SELF_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_double *alpha,
                                                                            cl_uint *offset_x, cl_mem *x,
                                                                            cl_uint *offset_y, cl_mem *y){
  if(*n==0)
    return;
  axpy_offset_generic((*command_queue)->kernels.axpy_offset_kernel_d, (*command_queue)->command_queue, n, alpha, offset_x, x, offset_y, y, offset_y, y);
}


void FC_FUNC_(gemm_volkov_d,GEMM_VOLKOV_D)(bigdft_command_queue *command_queue, char *transa, char *transb, cl_uint *m, cl_uint *n, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc) {
  if( *transa == 't' || *transa == 'c' || *transa == 'T' || *transa == 'C' ) {
    if ( *transb == 't' || *transb == 'c' || *transb == 'T' || *transb == 'C' ) {
      gemm_generic((*command_queue)->kernels.gemm_kernel_d_tatb, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
      gemm_generic((*command_queue)->kernels.gemm_kernel_d_ta, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
  } else {
    if ( *transb == 't' || *transb == 'c' || *transb == 'T' || *transb == 'C' ) {
      gemm_generic((*command_queue)->kernels.gemm_kernel_d_tb, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
      gemm_volkov_generic((*command_queue)->kernels.gemm_volkov_kernel_d, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
  }
}
void FC_FUNC_(gemm_d,GEMM_D)(bigdft_command_queue *command_queue, char *transa, char *transb, cl_uint *m, cl_uint *n, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc) {
  if( *transa == 't' || *transa == 'c' || *transa == 'T' || *transa == 'C' ) {
    if ( *transb == 't' || *transb == 'c' || *transb == 'T' || *transb == 'C' ) {
      gemm_generic((*command_queue)->kernels.gemm_kernel_d_tatb, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
      gemm_generic((*command_queue)->kernels.gemm_kernel_d_ta, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
  } else {
    if ( *transb == 't' || *transb == 'c' || *transb == 'T' || *transb == 'C' ) {
      gemm_generic((*command_queue)->kernels.gemm_kernel_d_tb, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
      gemm_generic((*command_queue)->kernels.gemm_kernel_d, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
  }
}
void FC_FUNC_(gemm_block_d,GEMM_BLOCK_D)(bigdft_command_queue *command_queue, char *transa, char *transb, cl_uint *m, cl_uint *n, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc) {
  if( *transa == 't' || *transa == 'c' || *transa == 'T' || *transa == 'C' ) {
    if ( *transb == 't' || *transb == 'c' || *transb == 'T' || *transb == 'C' ) {
      gemm_generic((*command_queue)->kernels.gemm_kernel_d_tatb, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
      gemm_generic((*command_queue)->kernels.gemm_kernel_d_ta, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
  } else {
    if ( *transb == 't' || *transb == 'c' || *transb == 'T' || *transb == 'C' ) {
      gemm_generic((*command_queue)->kernels.gemm_kernel_d_tb, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
      gemm_block_generic((*command_queue)->kernels.gemm_block_kernel_d, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
  }
}

void FC_FUNC_(gemmsy_d,GEMMSY_D)(bigdft_command_queue *command_queue, char *transa, char *transb, cl_uint *m, cl_uint *n, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc) {
  assert(m == n);
  if( *transa == 't' || *transa == 'c' || *transa == 'T' || *transa == 'C' ) {
    if ( *transb == 't' || *transb == 'c' || *transb == 'T' || *transb == 'C' ) {
      gemm_generic((*command_queue)->kernels.gemm_kernel_d_tatb, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
      gemm_generic((*command_queue)->kernels.gemm_kernel_d_ta, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
  } else {
    if ( *transb == 't' || *transb == 'c' || *transb == 'T' || *transb == 'C' ) {
      gemm_generic((*command_queue)->kernels.gemm_kernel_d_tb, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
      gemmsy_generic((*command_queue)->kernels.gemmsy_kernel_d, (*command_queue)->command_queue, m, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
  }
}

void FC_FUNC_(gemm_z,GEMM_Z)(bigdft_command_queue *command_queue, char *transa, char *transb, cl_uint *m, cl_uint *n, cl_uint *k, cl_double2 *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double2 *beta, cl_mem *c, cl_uint *ldc) {
  if( *transa == 't' || *transa == 'T' ) {
    if ( *transb == 't' || *transb == 'T' ) {
      zgemm_generic((*command_queue)->kernels.gemm_kernel_z_tatb, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if ( *transb == 'c' || *transb == 'C' ) {
      zgemm_generic((*command_queue)->kernels.gemm_kernel_z_tacb, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
      zgemm_generic((*command_queue)->kernels.gemm_kernel_z_ta, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
  } else if ( *transa == 'c' || *transa == 'C' ) {
    if ( *transb == 't' || *transb == 'T' ) {
      zgemm_generic((*command_queue)->kernels.gemm_kernel_z_catb, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if ( *transb == 'c' || *transb == 'C' ) {
      zgemm_generic((*command_queue)->kernels.gemm_kernel_z_cacb, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
      zgemm_generic((*command_queue)->kernels.gemm_kernel_z_ca, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
  } else {
    if ( *transb == 't' || *transb == 'T' ) {
      zgemm_generic((*command_queue)->kernels.gemm_kernel_z_tb, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if ( *transb == 'c' || *transb == 'C' ) {
      zgemm_generic((*command_queue)->kernels.gemm_kernel_z_cb, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
      zgemm_generic((*command_queue)->kernels.gemm_kernel_z, (*command_queue)->command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
  }
}

void FC_FUNC_(asum_self_d,ASUM_SELF_D)(bigdft_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work, cl_double *out) {
  if(*ndat==0){
    *out = 0.0;
    return;
  }
  size_t max_wgs = (*command_queue)->device_infos.MAX_WORK_GROUP_SIZE*2;
  cl_uint n = *ndat;
  cl_mem *input = in;
  cl_mem *output = work;
  cl_mem *tmp;
  do {
    reduction_generic((*command_queue)->kernels.reduction_kernel_d, *command_queue, &n, input, output);
    tmp = input;
    input = output;
    output = tmp;
    n = shrRoundUp(max_wgs,n)/max_wgs;
  } while(n>1);
  clEnqueueReadBuffer((*command_queue)->command_queue, *input, CL_TRUE, 0, sizeof(cl_double), out, 0, NULL, NULL);
}

void FC_FUNC_(nrm2sq_self_d,NRM2SQ_SELF_D)(bigdft_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work, cl_double *out) {
  if(*ndat==0){
   *out = 0.0;
   return;
  }
  size_t max_wgs = (*command_queue)->device_infos.MAX_WORK_GROUP_SIZE*2;
  cl_uint n = *ndat;
  cl_mem *input = in;
  cl_mem *output = work;
  cl_mem *tmp;
  reduction_generic((*command_queue)->kernels.reduction_dot_kernel_d, *command_queue, &n, input, output);
  input = work;
  output = in;
  n = shrRoundUp(max_wgs,n)/max_wgs;
  if(n>1) {
    do {
      reduction_generic((*command_queue)->kernels.reduction_kernel_d, *command_queue, &n, input, output);
      tmp = input;
      input = output;
      output = tmp;
      n = shrRoundUp(max_wgs,n)/max_wgs;
    } while(n>1);
  }
  clEnqueueReadBuffer((*command_queue)->command_queue, *input, CL_TRUE, 0, sizeof(cl_double), out, 0, NULL, NULL);
}

void FC_FUNC_(asum_d,ASUM_D)(bigdft_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work1, cl_mem *work2, cl_double *out) {
  if(*ndat==0){
   *out = 0.0;
   return;
  }
  size_t max_wgs = (*command_queue)->device_infos.MAX_WORK_GROUP_SIZE*2;
  cl_uint n = *ndat;
  cl_mem *input = in;
  cl_mem *output = work1;
  cl_mem *tmp;
  reduction_generic((*command_queue)->kernels.reduction_kernel_d, *command_queue, &n, input, output);
  input = work1;
  output = work2;
  n = shrRoundUp(max_wgs,n)/max_wgs;
  if(n>1) {
    do {
      reduction_generic((*command_queue)->kernels.reduction_kernel_d, *command_queue, &n, input, output);
      tmp = input;
      input = output;
      output = tmp;
      n = shrRoundUp(max_wgs,n)/max_wgs;
    } while(n>1);
  }
  clEnqueueReadBuffer((*command_queue)->command_queue, *input, CL_TRUE, 0, sizeof(cl_double), out, 0, NULL, NULL);
}

void FC_FUNC_(nrm2sq_d,NRM2SQ_D)(bigdft_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work1, cl_mem *work2, cl_double *out) {
  if(*ndat==0){
   *out = 0.0;
   return;
  }
  size_t max_wgs = (*command_queue)->device_infos.MAX_WORK_GROUP_SIZE*2;
  cl_uint n = *ndat;
  reduction_generic((*command_queue)->kernels.reduction_dot_kernel_d, *command_queue, &n, in, work1);
  cl_mem *input = work1;
  cl_mem *output = work2;
  cl_mem *tmp;
  n = shrRoundUp(max_wgs,n)/max_wgs;
  if(n>1) {
    do {
      reduction_generic((*command_queue)->kernels.reduction_kernel_d, *command_queue, &n, input, output);
      tmp = input;
      input = output;
      output = tmp;
      n = shrRoundUp(max_wgs,n)/max_wgs;
    } while(n>1);
  }
  clEnqueueReadBuffer((*command_queue)->command_queue, *input, CL_TRUE, 0, sizeof(cl_double), out, 0, NULL, NULL);
}

void FC_FUNC_(dot_d,DOT_D)(bigdft_command_queue *command_queue, cl_uint *ndat, cl_mem *x, cl_mem *y, cl_mem *work1, cl_mem *work2, cl_double *out) {
  if(*ndat==0){
   *out = 0.0;
   return;
  }
  size_t max_wgs = (*command_queue)->device_infos.MAX_WORK_GROUP_SIZE*2;
  cl_uint n = *ndat;
  dot_generic((*command_queue)->kernels.dot_kernel_d, *command_queue, &n, x, y, work1);
  cl_mem *input=work1;
  cl_mem *output=work2;
  cl_mem *tmp;
  n = shrRoundUp(max_wgs,n)/max_wgs;
  if(n>1) {
    do {
      reduction_generic((*command_queue)->kernels.reduction_kernel_d, *command_queue, &n, input, output);
      tmp = input;
      input = output;
      output = tmp;
      n = shrRoundUp(max_wgs,n)/max_wgs;
    } while(n>1);
  }
  clEnqueueReadBuffer((*command_queue)->command_queue, *input, CL_TRUE, 0, sizeof(cl_double), out, 0, NULL, NULL);
}

void FC_FUNC_(dot_d_async,DOT_D_ASYNC)(bigdft_command_queue *command_queue, cl_uint *ndat, cl_mem *x, cl_mem *y, cl_mem *work1, cl_mem *work2, cl_double *out) {
  if(*ndat==0){
   *out = 0.0;
   return;
  }
  size_t max_wgs = (*command_queue)->device_infos.MAX_WORK_GROUP_SIZE*2;
  cl_uint n = *ndat;
  dot_generic((*command_queue)->kernels.dot_kernel_d, *command_queue, &n, x, y, work1);
  cl_mem *input=work1;
  cl_mem *output=work2;
  cl_mem *tmp;
  n = shrRoundUp(max_wgs,n)/max_wgs;
  if(n>1) {
    do {
      reduction_generic((*command_queue)->kernels.reduction_kernel_d, *command_queue, &n, input, output);
      tmp = input;
      input = output;
      output = tmp;
      n = shrRoundUp(max_wgs,n)/max_wgs;
    } while(n>1);
  }
  clEnqueueReadBuffer((*command_queue)->command_queue, *input, CL_FALSE, 0, sizeof(cl_double), out, 0, NULL, NULL);
}


void create_reduction_kernels(struct bigdft_kernels * kernels){
    cl_int ciErrNum = CL_SUCCESS;
    kernels->axpy_offset_kernel_d=clCreateKernel(reductionProgram,"axpy_offsetKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->axpy_kernel_d=clCreateKernel(reductionProgram,"axpyKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->scal_kernel_d=clCreateKernel(reductionProgram,"scalKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->reduction_kernel_d=clCreateKernel(reductionProgram,"reductionKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->reduction_dot_kernel_d=clCreateKernel(reductionProgram,"reduction_dotKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->copy_kernel_d=clCreateKernel(reductionProgram,"copyKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->dot_kernel_d=clCreateKernel(reductionProgram,"dotKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->set_kernel_d=clCreateKernel(reductionProgram,"setKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->gemm_block_kernel_d=clCreateKernel(dgemmProgram,"gemm_blockKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->gemmsy_kernel_d=clCreateKernel(dgemmProgram,"gemmsyKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->gemm_volkov_kernel_d=clCreateKernel(dgemmProgram,"gemm_volkovKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->gemm_kernel_d=clCreateKernel(dgemmProgram,"gemmKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->gemm_kernel_z=clCreateKernel(dgemmProgram,"gemmKernel_z",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->gemm_kernel_d_tb=clCreateKernel(dgemmProgram,"gemmKernel_d_tb",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->gemm_kernel_z_tb=clCreateKernel(dgemmProgram,"gemmKernel_z_tb",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->gemm_kernel_z_cb=clCreateKernel(dgemmProgram,"gemmKernel_z_cb",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->gemm_kernel_d_ta=clCreateKernel(dgemmProgram,"gemmKernel_d_ta",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->gemm_kernel_z_ta=clCreateKernel(dgemmProgram,"gemmKernel_z_ta",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->gemm_kernel_z_ca=clCreateKernel(dgemmProgram,"gemmKernel_z_ca",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->gemm_kernel_d_tatb=clCreateKernel(dgemmProgram,"gemmKernel_d_tatb",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->gemm_kernel_z_tatb=clCreateKernel(dgemmProgram,"gemmKernel_z_tatb",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->gemm_kernel_z_catb=clCreateKernel(dgemmProgram,"gemmKernel_z_catb",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->gemm_kernel_z_tacb=clCreateKernel(dgemmProgram,"gemmKernel_z_tacb",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    kernels->gemm_kernel_z_cacb=clCreateKernel(dgemmProgram,"gemmKernel_z_cacb",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
}

void build_reduction_programs(bigdft_context * context){
    struct bigdft_device_infos infos;
    get_context_devices_infos(context, &infos);
    cl_int ciErrNum = CL_SUCCESS;
    char * code = generate_reduction_program(&infos);
    reductionProgram = clCreateProgramWithSource((*context)->context, 1, (const char**) &(code), NULL, &ciErrNum);
    free(code);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(reductionProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build reduction program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(reductionProgram, oclGetFirstDev((*context)->context), CL_PROGRAM_BUILD_LOG, sizeof(cBuildLog), cBuildLog, NULL );
        fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    dgemmProgram = clCreateProgramWithSource((*context)->context, 1, (const char**) &dgemm_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(dgemmProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build dgemm program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(dgemmProgram, oclGetFirstDev((*context)->context), CL_PROGRAM_BUILD_LOG, sizeof(cBuildLog), cBuildLog, NULL );
        fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
}

void clean_reduction_kernels(struct bigdft_kernels * kernels){
  cl_int ciErrNum;
  ciErrNum = clReleaseKernel(kernels->reduction_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->reduction_dot_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->axpy_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->axpy_offset_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->scal_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->copy_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->dot_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->set_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->gemm_block_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->gemmsy_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->gemm_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->gemm_volkov_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->gemm_kernel_z);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->gemm_kernel_d_tb);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->gemm_kernel_z_tb);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->gemm_kernel_z_cb);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->gemm_kernel_d_ta);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->gemm_kernel_z_ta);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->gemm_kernel_z_ca);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->gemm_kernel_d_tatb);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->gemm_kernel_z_tatb);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->gemm_kernel_z_catb);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->gemm_kernel_z_tacb);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->gemm_kernel_z_cacb);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
}
void clean_reduction_programs(){
  cl_int ciErrNum;
  ciErrNum = clReleaseProgram(reductionProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
  ciErrNum = clReleaseProgram(dgemmProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
}
