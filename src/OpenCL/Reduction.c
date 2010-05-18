#include "OpenCL_wrappers.h"

char * dgemm_program="\
//size is supposed to be 16*16\n\
#define BUFFER_SIZE 16\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
__kernel void gemmKernel_d( uint m, uint n, uint k, double alpha, __global const double *a, uint lda, __global const double *b, uint ldb, double beta, __global double * c, uint ldc, __local double *tmp1, __local double *tmp2){\n\
  size_t i = get_local_id(0);\n\
  size_t j = get_local_id(1);\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  \n\
  size_t index = 0;\n\
  double result = 0.0;\n\
  while( index < k) {\n\
    //load first matrix in tmp1\n\
    tmp1[j*(BUFFER_SIZE) + i] = (ig < m && (index + j) < k) ? a[(index+j)*lda + ig] : 0.0;\n\
    //load second matrix in tmp2\n\
    tmp2[j*(BUFFER_SIZE) + i] = (jg < n && (index + i) < k) ? b[(jg)*ldb + index+i] : 0.0;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t sumi=0; sumi<BUFFER_SIZE; sumi++)\n\
      result += tmp1[sumi*(BUFFER_SIZE) + i] * tmp2[ j*(BUFFER_SIZE) + sumi];\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    index += BUFFER_SIZE;\n\
  }\n\
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
    //load first matrix in tmp1\n\
    tmp1[i*(BUFFER_SIZE+1) + j] = (igt < m && (index + i) <  k) ? a[(igt)*lda + index+i] : 0.0;\n\
    //load second matrix in tmp2\n\
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
    //load first matrix in tmp1\n\
    tmp1[j*(BUFFER_SIZE) + i] = (ig < m && (index + j) <  k) ? a[(index+j)*lda + ig] : 0.0;\n\
    //load second matrix in tmp2\n\
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
    //load first matrix in tmp1\n\
    tmp1[i*(BUFFER_SIZE+1) + j] = (igt < m && (index + i) <  k) ? a[(igt)*lda + index+i] : 0.0;\n\
    //load second matrix in tmp2\n\
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
    //load first matrix in tmp1\n\
    tmp1[j*(BUFFER_SIZE)] = (condm && (index + j) < k) ? a[index*lda] : (double2)(0.0, 0.0);\n\
    //load second matrix in tmp2\n\
    tmp2[i] = (condn && (index + i) < k) ? b[index] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*BUFFER_SIZE];\n\
      b_t = tmp2[k];\n\
      result.x = mad(a_t.x,b_t.x,mad(-a_t.y,b_t.y,result.x));\n\
      result.y = mad(a_t.x,b_t.y,mad(a_t.y,b_t.x,result.y));\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x +=-alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x +=-beta.y*c_t.y;\n\
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
//  tmp1 += i;\n\
  a += igt*lda + i;\n\
  b += jg*ldb + i;\n\
  c += jg*ldc + ig;\n\
  bool condm = igt < m;\n\
  bool condn = jg < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    //load first matrix in tmp1\n\
    tmp1[i*(BUFFER_SIZE+1)+j] = (condm && (index + i) < k) ? a[index] : (double2)(0.0, 0.0);\n\
    //load second matrix in tmp2\n\
    tmp2[i] = (condn && (index + i) < k) ? b[index] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*(BUFFER_SIZE+1) + i];\n\
      b_t = tmp2[k];\n\
      result.x = mad(a_t.x,b_t.x,mad(-a_t.y,b_t.y,result.x));\n\
      result.y = mad(a_t.x,b_t.y,mad(a_t.y,b_t.x,result.y));\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x +=-alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x +=-beta.y*c_t.y;\n\
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
//  tmp2 += j*(BUFFER_SIZE);\n\
  tmp1 += i;\n\
  a += j*lda + ig;\n\
  b += j*ldb + jgt;\n\
  c += jg*ldc + ig;\n\
  bool condm = ig < m;\n\
  bool condn = jgt < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    //load first matrix in tmp1\n\
    tmp1[j*(BUFFER_SIZE)] = (condm && (index + j) < k) ? a[index*lda] : (double2)(0.0, 0.0);\n\
    //load second matrix in tmp2\n\
    tmp2[i*(BUFFER_SIZE+1) + j] = (condn && (index + j) < k) ? b[index*ldb] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*BUFFER_SIZE];\n\
      b_t = tmp2[j*(BUFFER_SIZE+1) + k];\n\
      result.x = mad(a_t.x,b_t.x,mad(-a_t.y,b_t.y,result.x));\n\
      result.y = mad(a_t.x,b_t.y,mad(a_t.y,b_t.x,result.y));\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x +=-alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x +=-beta.y*c_t.y;\n\
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
//  tmp2 += j*(BUFFER_SIZE);\n\
//  tmp1 += i;\n\
  a += igt*lda + i;\n\
  b += j*ldb + jgt;\n\
  c += jg*ldc + ig;\n\
  bool condm = igt < m;\n\
  bool condn = jgt < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    //load first matrix in tmp1\n\
    tmp1[i*(BUFFER_SIZE+1) + j] = (condm && (index + i) < k) ? a[index] : (double2)(0.0, 0.0);\n\
    //load second matrix in tmp2\n\
    tmp2[i*(BUFFER_SIZE+1) + j] = (condn && (index + j) < k) ? b[index*ldb] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*(BUFFER_SIZE+1) + i];\n\
      b_t = tmp2[j*(BUFFER_SIZE+1) + k];\n\
      result.x = mad(a_t.x,b_t.x,mad(-a_t.y,b_t.y,result.x));\n\
      result.y = mad(a_t.x,b_t.y,mad(a_t.y,b_t.x,result.y));\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x +=-alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x +=-beta.y*c_t.y;\n\
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
//  tmp1 += i;\n\
  a += igt*lda + i;\n\
  b += jg*ldb + i;\n\
  c += jg*ldc + ig;\n\
  bool condm = igt < m;\n\
  bool condn = jg < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    //load first matrix in tmp1\n\
    tmp1[i*(BUFFER_SIZE+1)+j] = (condm && (index + i) < k) ? a[index] : (double2)(0.0, 0.0);\n\
    //load second matrix in tmp2\n\
    tmp2[i] = (condn && (index + i) < k) ? b[index] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*(BUFFER_SIZE+1) + i];\n\
      b_t = tmp2[k];\n\
      result.x = mad(a_t.x,b_t.x,mad(a_t.y,b_t.y,result.x));\n\
      result.y = mad(a_t.x,b_t.y,mad(-a_t.y,b_t.x,result.y));\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x +=-alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x +=-beta.y*c_t.y;\n\
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
//  tmp2 += j*(BUFFER_SIZE);\n\
  tmp1 += i;\n\
  a += j*lda + ig;\n\
  b += j*ldb + jgt;\n\
  c += jg*ldc + ig;\n\
  bool condm = ig < m;\n\
  bool condn = jgt < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    //load first matrix in tmp1\n\
    tmp1[j*(BUFFER_SIZE)] = (condm && (index + j) < k) ? a[index*lda] : (double2)(0.0, 0.0);\n\
    //load second matrix in tmp2\n\
    tmp2[i*(BUFFER_SIZE+1) + j] = (condn && (index + j) < k) ? b[index*ldb] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*BUFFER_SIZE];\n\
      b_t = tmp2[j*(BUFFER_SIZE+1) + k];\n\
      result.x = mad(a_t.x,b_t.x,mad(-a_t.y,-b_t.y,result.x));\n\
      result.y = mad(a_t.x,-b_t.y,mad(a_t.y,b_t.x,result.y));\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x +=-alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x +=-beta.y*c_t.y;\n\
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
//  tmp2 += j*(BUFFER_SIZE);\n\
//  tmp1 += i;\n\
  a += igt*lda + i;\n\
  b += j*ldb + jgt;\n\
  c += jg*ldc + ig;\n\
  bool condm = igt < m;\n\
  bool condn = jgt < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    //load first matrix in tmp1\n\
    tmp1[i*(BUFFER_SIZE+1) + j] = (condm && (index + i) < k) ? a[index] : (double2)(0.0, 0.0);\n\
    //load second matrix in tmp2\n\
    tmp2[i*(BUFFER_SIZE+1) + j] = (condn && (index + j) < k) ? b[index*ldb] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*(BUFFER_SIZE+1) + i];\n\
      b_t = tmp2[j*(BUFFER_SIZE+1) + k];\n\
      result.x = mad(a_t.x,b_t.x,mad(a_t.y,-b_t.y,result.x));\n\
      result.y = mad(a_t.x,-b_t.y,mad(-a_t.y,b_t.x,result.y));\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x +=-alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x +=-beta.y*c_t.y;\n\
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
//  tmp2 += j*(BUFFER_SIZE);\n\
//  tmp1 += i;\n\
  a += igt*lda + i;\n\
  b += j*ldb + jgt;\n\
  c += jg*ldc + ig;\n\
  bool condm = igt < m;\n\
  bool condn = jgt < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    //load first matrix in tmp1\n\
    tmp1[i*(BUFFER_SIZE+1) + j] = (condm && (index + i) < k) ? a[index] : (double2)(0.0, 0.0);\n\
    //load second matrix in tmp2\n\
    tmp2[i*(BUFFER_SIZE+1) + j] = (condn && (index + j) < k) ? b[index*ldb] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*(BUFFER_SIZE+1) + i];\n\
      b_t = tmp2[j*(BUFFER_SIZE+1) + k];\n\
      result.x = mad(a_t.x,b_t.x,mad(a_t.y,b_t.y,result.x));\n\
      result.y = mad(a_t.x,b_t.y,mad(-a_t.y,b_t.x,result.y));\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x +=-alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x +=-beta.y*c_t.y;\n\
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
//  tmp2 += j*(BUFFER_SIZE);\n\
//  tmp1 += i;\n\
  a += igt*lda + i;\n\
  b += j*ldb + jgt;\n\
  c += jg*ldc + ig;\n\
  bool condm = igt < m;\n\
  bool condn = jgt < n;\n\
  while( index < k) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    //load first matrix in tmp1\n\
    tmp1[i*(BUFFER_SIZE+1) + j] = (condm && (index + i) < k) ? a[index] : (double2)(0.0, 0.0);\n\
    //load second matrix in tmp2\n\
    tmp2[i*(BUFFER_SIZE+1) + j] = (condn && (index + j) < k) ? b[index*ldb] : (double2)(0.0, 0.0);\n\
    index += BUFFER_SIZE;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    #pragma unroll\n\
    for(size_t k=0; k<BUFFER_SIZE; k++){\n\
      a_t = tmp1[k*(BUFFER_SIZE+1) + i];\n\
      b_t = tmp2[j*(BUFFER_SIZE+1) + k];\n\
      result.x = mad(a_t.x,b_t.x,mad(-a_t.y,-b_t.y,result.x));\n\
      result.y = mad(a_t.x,-b_t.y,mad(a_t.y,b_t.x,result.y));\n\
    }\n\
  }\n\
  if(condm && condn){\n\
    double2 final_result = (double2)(0.0, 0.0);\n\
    final_result.x += alpha.x*result.x;\n\
    final_result.x +=-alpha.y*result.y;\n\
    final_result.y += alpha.x*result.y;\n\
    final_result.y += alpha.y*result.x;\n\
    double2 c_t = c[0];\n\
    final_result.x += beta.x*c_t.x;\n\
    final_result.x +=-beta.y*c_t.y;\n\
    final_result.y += beta.x*c_t.y;\n\
    final_result.y += beta.y*c_t.x;\n\
    c[0] = final_result;\n\
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

void inline zgemm_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *m, cl_uint *n, cl_uint *k, cl_double2 *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double2 *beta, cl_mem *c, cl_uint *ldc) {
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
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue zgemm kernel!");
}

void inline gemm_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *m, cl_uint *n, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc) {
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
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue gemm kernel!");
}

void inline set_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, cl_double *val, cl_mem *x) {
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

void inline scal_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *in, cl_mem *out) {
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

void inline axpy_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *x, cl_mem *y, cl_mem *out) {
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

void inline axpy_offset_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, cl_double *alpha, 
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
  clSetKernelArg(kernel, i++,sizeof(cl_double)*block_size_i*2, NULL);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i*2,*ndat)/2 };
  ciErrNum = clEnqueueNDRangeKernel(*command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue reduction kernel!");
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
cl_kernel gemm_kernel_z;
cl_kernel gemm_kernel_z_ta;
cl_kernel gemm_kernel_z_tb;
cl_kernel gemm_kernel_z_tatb;
cl_kernel gemm_kernel_z_ca;
cl_kernel gemm_kernel_z_cb;
cl_kernel gemm_kernel_z_cacb;
cl_kernel gemm_kernel_z_tacb;
cl_kernel gemm_kernel_z_catb;
cl_kernel gemm_kernel_d_tb;
cl_kernel gemm_kernel_d_ta;
cl_kernel gemm_kernel_d_tatb;
cl_program reductionProgram;
cl_program dgemmProgram;

void FC_FUNC_(set_d,SET_D)(cl_command_queue *command_queue, cl_uint *n, cl_double *val, cl_mem *x){
  if(*n==0) return;
  set_generic(set_kernel_d, command_queue, n, val, x);
}

void FC_FUNC_(copy_d,COPY_D)(cl_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *out){
  if(*n==0) return;
  copy_generic(copy_kernel_d, command_queue, n, in, out);
}

void FC_FUNC_(scal_self_d,SCAL_SELF_D)(cl_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *inout){
  if(*n==0)
    return;
  scal_generic(scal_kernel_d, command_queue, n, alpha, inout, inout);
}

void FC_FUNC_(scal_d,SCAL_D)(cl_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *in, cl_mem *out){
  if(*n==0)
    return;
  scal_generic(scal_kernel_d, command_queue, n, alpha, in, out);
}

void FC_FUNC_(axpy_self_d,AXPY_SELF_D)(cl_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *in, cl_mem *inout){
  if(*n==0)
    return;
  axpy_generic(axpy_kernel_d, command_queue, n, alpha, in, inout, inout);
}

void FC_FUNC_(axpy_d,AXPY_D)(cl_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *x, cl_mem *y, cl_mem *z){
  if(*n==0)
    return;
  axpy_generic(axpy_kernel_d, command_queue, n, alpha, x, y, z);
}

void FC_FUNC_(axpy_offset_d,AXPY_OFFSET_D)(cl_command_queue *command_queue, cl_uint *n, cl_double *alpha,
                                                                            cl_uint *offset_x, cl_mem *x,
                                                                            cl_uint *offset_y, cl_mem *y,
                                                                            cl_uint *offset_z, cl_mem *z){
  if(*n==0)
    return;
  axpy_offset_generic(axpy_offset_kernel_d, command_queue, n, alpha, offset_x, x, offset_y, y, offset_z, z);
}

void FC_FUNC_(axpy_offset_self_d,AXPY_OFFSET_SELF_D)(cl_command_queue *command_queue, cl_uint *n, cl_double *alpha,
                                                                            cl_uint *offset_x, cl_mem *x,
                                                                            cl_uint *offset_y, cl_mem *y){
  if(*n==0)
    return;
  axpy_offset_generic(axpy_offset_kernel_d, command_queue, n, alpha, offset_x, x, offset_y, y, offset_y, y);
}

void FC_FUNC_(asum_self_d,ASUM_SELF_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work, cl_double *out) {
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
  clEnqueueReadBuffer(*command_queue, *input, CL_TRUE, 0, sizeof(cl_double), out, 0, NULL, NULL);
}

void FC_FUNC_(nrm2sq_self_d,NRM2SQ_SELF_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work, cl_double *out) {
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
  clEnqueueReadBuffer(*command_queue, *input, CL_TRUE, 0, sizeof(cl_double), out, 0, NULL, NULL);
}

void FC_FUNC_(gemm_d,GEMM_D)(cl_command_queue *command_queue, char *transa, char *transb, cl_uint *m, cl_uint *n, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc) {
  if( *transa == 'n' && *transb == 'n' )
    gemm_generic(gemm_kernel_d, command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
  else if ( *transa == 'n' && ( *transb == 't' || *transb == 'c' ) )
    gemm_generic(gemm_kernel_d_tb, command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
  else if ( *transb == 'n' && ( *transa == 't' || *transa == 'c' ) )
    gemm_generic(gemm_kernel_d_ta, command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
  else
    gemm_generic(gemm_kernel_d_tatb, command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}
void FC_FUNC_(gemm_z,GEMM_Z)(cl_command_queue *command_queue, char *transa, char *transb, cl_uint *m, cl_uint *n, cl_uint *k, cl_double2 *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double2 *beta, cl_mem *c, cl_uint *ldc) {
  if( *transa == 't' ) {
    if ( *transb == 't' ) {
      zgemm_generic(gemm_kernel_z_tatb, command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if ( *transb == 'c' ) {
      zgemm_generic(gemm_kernel_z_tacb, command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
      zgemm_generic(gemm_kernel_z_ta, command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
  } else if ( *transa == 'c' ) {
    if ( *transb == 't' ) {
      zgemm_generic(gemm_kernel_z_catb, command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if ( *transb == 'c' ) {
      zgemm_generic(gemm_kernel_z_cacb, command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
      zgemm_generic(gemm_kernel_z_ca, command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
  } else {
    if ( *transb == 't' ) {
      zgemm_generic(gemm_kernel_z_tb, command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else if ( *transb == 'c' ) {
      zgemm_generic(gemm_kernel_z_cb, command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
      zgemm_generic(gemm_kernel_z, command_queue, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
  }
}

void FC_FUNC_(asum_d,ASUM_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work1, cl_mem *work2, cl_double *out) {
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
  clEnqueueReadBuffer(*command_queue, *input, CL_TRUE, 0, sizeof(cl_double), out, 0, NULL, NULL);
}

void FC_FUNC_(nrm2sq_d,NRM2SQ_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work1, cl_mem *work2, cl_double *out) {
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
  clEnqueueReadBuffer(*command_queue, *input, CL_TRUE, 0, sizeof(cl_double), out, 0, NULL, NULL);
}

void FC_FUNC_(dot_d,DOT_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *x, cl_mem *y, cl_mem *work1, cl_mem *work2, cl_double *out) {
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
  clEnqueueReadBuffer(*command_queue, *input, CL_TRUE, 0, sizeof(cl_double), out, 0, NULL, NULL);
}

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
    gemm_kernel_z=clCreateKernel(dgemmProgram,"gemmKernel_z",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    gemm_kernel_d_tb=clCreateKernel(dgemmProgram,"gemmKernel_d_tb",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    gemm_kernel_z_tb=clCreateKernel(dgemmProgram,"gemmKernel_z_tb",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    gemm_kernel_z_cb=clCreateKernel(dgemmProgram,"gemmKernel_z_cb",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    gemm_kernel_d_ta=clCreateKernel(dgemmProgram,"gemmKernel_d_ta",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    gemm_kernel_z_ta=clCreateKernel(dgemmProgram,"gemmKernel_z_ta",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    gemm_kernel_z_ca=clCreateKernel(dgemmProgram,"gemmKernel_z_ca",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    gemm_kernel_d_tatb=clCreateKernel(dgemmProgram,"gemmKernel_d_tatb",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    gemm_kernel_z_tatb=clCreateKernel(dgemmProgram,"gemmKernel_z_tatb",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    gemm_kernel_z_catb=clCreateKernel(dgemmProgram,"gemmKernel_z_catb",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    gemm_kernel_z_tacb=clCreateKernel(dgemmProgram,"gemmKernel_z_tacb",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    gemm_kernel_z_cacb=clCreateKernel(dgemmProgram,"gemmKernel_z_cacb",&ciErrNum);
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
  ciErrNum = clReleaseKernel(gemm_kernel_z);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(gemm_kernel_d_tb);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(gemm_kernel_z_tb);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(gemm_kernel_d_ta);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(gemm_kernel_z_ta);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(gemm_kernel_d_tatb);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(gemm_kernel_z_tatb);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseProgram(reductionProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
  ciErrNum = clReleaseProgram(dgemmProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
}
