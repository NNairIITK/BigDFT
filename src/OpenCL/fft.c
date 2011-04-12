#include "OpenCL_wrappers.h"
#include "fft_generator.h"

char * fft_program="\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#define TWOPI 3.14159265358979323846264338327950288419716939937510*2\n\
#define LINE_NUMBER 16\n\
#define BUFFER_DEPTH LINE_NUMBER+1\n\
#define FFT_LENGTH 16\n\
#define radix2m(il,jl,N,A,B,in,out) \
{ \
  double2 tmp,val,w;\
  int a,b,p,r;\
  b = jl / (2*A);\
  r = jl % (2*A);\
/*  p = r / A;*/\
  a = r % A;\
  val = in[A*b+a][il];\
  tmp.x = val.x;\
  tmp.y = val.y;\
  val = in[(N/2)+A*b+a][il];\
  w.x = cos((TWOPI/(A*2))*r);\
  w.y = sin((TWOPI/(A*2))*r);\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  out[jl][il]=tmp;\
} \n\
#define radix3m(il,jl,N,A,B,in,out) \
{ \
  double2 tmp,val,w;\
  int a,b,p,r;\
  b = jl / (3*A);\
  r = jl % (3*A);\
  a = r % A;\
  val = in[A*b+a][il];\
  tmp.x = val.x;\
  tmp.y = val.y;\
  val = in[(N/3)+A*b+a][il];\
  w.x = cos((TWOPI/(A*3))*r);\
  w.y = sin((TWOPI/(A*3))*r);\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/3)*2+A*b+a][il];\
  w.x = cos((TWOPI/(A*3))*r*2);\
  w.y = sin((TWOPI/(A*3))*r*2);\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  out[jl][il]=tmp;\
} \n\
#define radix5m(il,jl,N,A,B,in,out) \
{ \
  double2 tmp,val,w;\
  int a,b,p,r;\
  b = jl / (5*A);\
  r = jl % (5*A);\
  a = r % A;\
  val = in[A*b+a][il];\
  tmp.x = val.x;\
  tmp.y = val.y;\
  val = in[(N/5)+A*b+a][il];\
  w.x = cos((TWOPI/(A*5))*r);\
  w.y = sin((TWOPI/(A*5))*r);\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/5)*2+A*b+a][il];\
  w.x = cos((TWOPI/(A*5))*r*2);\
  w.y = sin((TWOPI/(A*5))*r*2);\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/5)*3+A*b+a][il];\
  w.x = cos((TWOPI/(A*5))*r*3);\
  w.y = sin((TWOPI/(A*5))*r*3);\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/5)*4+A*b+a][il];\
  w.x = cos((TWOPI/(A*5))*r*4);\
  w.y = sin((TWOPI/(A*5))*r*4);\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  out[jl][il]=tmp;\
} \n\
inline void radix2(size_t il, size_t jl, uint N, uint A, uint B, __local double2 in[FFT_LENGTH][BUFFER_DEPTH], __local double2 out[FFT_LENGTH][BUFFER_DEPTH]) {\n\
   double2 tmp, val,w;\n\
   int a,b,r;\n\
   b = jl / (2*A);\n\
   r = jl % (2*A);\n\
   a = r % A;\n\
\n\
  val = in[A*b+a][il];\n\
  tmp.x = val.x;\n\
  tmp.y = val.y;\n\
  val = in[(N/2)+A*b+a][il];\n\
  w.x = cos((TWOPI/(A*2))*r);\n\
  w.y = sin((TWOPI/(A*2))*r);\n\
  tmp.x += val.x * w.x;\n\
  tmp.x += val.y * w.y;\n\
  tmp.y -= val.x * w.y;\n\
  tmp.y += val.y * w.x;\n\
  out[jl][il]=tmp;\n\
   \n\
}\n\
\n\
inline void radix4(size_t il, size_t jl, uint N, uint A, uint B, __local double2 in[FFT_LENGTH][BUFFER_DEPTH], __local double2 out[FFT_LENGTH][BUFFER_DEPTH]) {\n\
   unsigned int i;\n\
   double2 tmp=(double2)(0.0, 0.0);\n\
   double2 val;\n\
   double wx,wy;\n\
   int a,b,p,r;\n\
   b = jl/(4*A);\n\
   r = (jl - b*4*A);\n\
   p = r/A;\n\
   a = r - A*p;\n\
\n\
   for(i=0; i<4; i++){\n\
     val = in[(N/4)*i+A*b+a][il];\n\
     wx = cos(TWOPI*r*i/((double)A*4));\n\
     wy = -sin(TWOPI*r*i/((double)A*4));\n\
     tmp.x += val.x * wx - val.y * wy;\n\
     tmp.y += val.x * wy + val.y * wx;\n\
   }\n\
\n\
   out[jl][il]=tmp;\n\
   \n\
}\n\
inline void radix16(size_t il, size_t jl, uint N, uint A, uint B, __local double2 in[FFT_LENGTH][BUFFER_DEPTH], __local double2 out[FFT_LENGTH][BUFFER_DEPTH]) {\n\
   unsigned int i;\n\
   double2 tmp=(double2)(0.0, 0.0);\n\
   double2 val;\n\
   double wx,wy;\n\
   int a,b,p,r;\n\
   b = jl/(16*A);\n\
   r = (jl - b*16*A);\n\
   p = r/A;\n\
   a = r - A*p;\n\
\n\
   for(i=0; i<16; i++){\n\
     val = in[(N/16)*i+A*b+a][il];\n\
     wx = cos(TWOPI*r*i/((double)A*16));\n\
     wy = -sin(TWOPI*r*i/((double)A*16));\n\
     tmp.x += val.x * wx - val.y * wy;\n\
     tmp.y += val.x * wy + val.y * wx;\n\
   }\n\
\n\
   out[jl][il]=tmp;\n\
   \n\
}\n\
\n\
__kernel void fftKernel_test_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out){\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 (*tmp_in)[BUFFER_DEPTH];\n\
__local double2 (*tmp_out)[BUFFER_DEPTH];\n\
__local double2 (*tmp)[BUFFER_DEPTH];\n\
\n\
  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t jlt = il;\n\
  size_t ilt = jl;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + il;\n\
  tmp1[jl][il] = jl < n ? psi[jgt + ( jl ) * ndat] : 0.0;\n\
  \n\
  unsigned int A,B;\n\
  unsigned int j;\n\
  A=1;\n\
  B=FFT_LENGTH;\n\
  tmp_in=tmp1;\n\
  tmp_out=tmp2;\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  #pragma unroll\n\
  for(j=0;j<4;j++){\n\
    B /= 2;\n\
    radix2(il, jl, n, A, B, tmp_in, tmp_out);\n\
//    tmp_out[il][jl]=tmp_in[il][jl];\n\
    tmp = tmp_in; tmp_in = tmp_out; tmp_out = tmp;\n\
    A *= 2;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
  }\n\
/*    radix2m(il, jl, n, 1, 8, tmp_in, tmp_out);\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    radix2m(il, jl, n, 2, 4, tmp_out, tmp_in);\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    radix2m(il, jl, n, 4, 2, tmp_in, tmp_out);\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    radix2m(il, jl, n, 8, 1, tmp_out, tmp_in);\n\
    barrier(CLK_LOCAL_MEM_FENCE);*/\n\
/*  for(j=0;j<2;j++){\n\
    B /= 4;\n\
    radix4(il, jl, n, A, B, tmp_in, tmp_out);\n\
    tmp = tmp_in; tmp_in = tmp_out; tmp_out = tmp;\n\
    A *= 4;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
  }*/\n\
/*  for(j=0;j<1;j++){\n\
    B /= 16;\n\
    radix16(il, jl, n, A, B, tmp_in, tmp_out);\n\
    tmp = tmp_in; tmp_in = tmp_out; tmp_out = tmp;\n\
    A *= 16;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
  }*/\n\
\n\
  if(il<n)\n\
    out[jg*n+il] = tmp_in[il][jl];\n\
}\n\
inline void radix2_2(size_t il, size_t jl, size_t offset, uint N, uint A, uint B, __local double2 in[FFT_LENGTH][BUFFER_DEPTH], __local double2 out[FFT_LENGTH][BUFFER_DEPTH]) {\n\
   double2 tmp,tmp2, val,w;\n\
   int a,b,r,r2;\n\
   b = jl / (2*A);\n\
   r = jl % (2*A);\n\
   r2 = (jl+offset) % (2*A);\n\
   a = r % A;\n\
\n\
  val = in[A*b+a][il];\n\
  tmp2.x = tmp.x = val.x;\n\
  tmp2.y = tmp.y = val.y;\n\
  val = in[(N/2)+A*b+a][il];\n\
  w.x = cos((TWOPI/(A*2))*r);\n\
  w.y = sin((TWOPI/(A*2))*r);\n\
  tmp.x += val.x * w.x;\n\
  tmp.x += val.y * w.y;\n\
  tmp.y -= val.x * w.y;\n\
  tmp.y += val.y * w.x;\n\
  w.x = cos((TWOPI/(A*2))*r2);\n\
  w.y = sin((TWOPI/(A*2))*r2);\n\
  tmp2.x += val.x * w.x;\n\
  tmp2.x += val.y * w.y;\n\
  tmp2.y -= val.x * w.y;\n\
  tmp2.y += val.y * w.x;\n\
  out[jl][il]=tmp;\n\
  out[jl+offset][il]=tmp2;\n\
   \n\
}\n\
__kernel void fftKernel_test_2_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out){\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 (*tmp_in)[BUFFER_DEPTH];\n\
__local double2 (*tmp_out)[BUFFER_DEPTH];\n\
__local double2 (*tmp)[BUFFER_DEPTH];\n\
\n\
  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1)*2;\n\
  size_t jl2 = jl+1;\n\
  size_t jg = get_global_id(1)*2;\n\
  size_t jlt = il;\n\
  size_t ilt = jl;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1)*2 - ndat ) : jg;\n\
  jgt = jg - jl + il;\n\
  tmp1[jl][il] = jl < n ? psi[jgt + ( jl ) * ndat] : 0.0;\n\
  tmp1[jl2][il] = jl2 < n ? psi[jgt + ( jl2 ) * ndat] : 0.0;\n\
  \n\
  unsigned int A,B;\n\
  unsigned int j;\n\
  A=1;\n\
  B=FFT_LENGTH;\n\
  tmp_in=tmp1;\n\
  tmp_out=tmp2;\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
//  #pragma unroll\n\
//  for(j=0;j<4;j++){\n\
    radix2_2(il, jl, 1, n, 1, 8, tmp_in, tmp_out);\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    radix2_2(il, (jl/4)*4+((jl/2)%2), 2, n, 2, 4, tmp_out, tmp_in);\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    radix2_2(il, (jl/8)*8+((jl/2)%4), 4, n, 4, 2, tmp_in, tmp_out);\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    radix2_2(il, jl/2, 8, n, 8, 1, tmp_out, tmp_in);\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
//  }\n\
\n\
  if(il<n){\n\
    out[jg*n+il] = tmp_in[il][jl];\n\
    out[(jg+1)*n+il] = tmp_in[il][jl2];\n\
  }\n\
}\n\
#undef LINE_NUMBER\n\
#define LINE_NUMBER 8\n\
#undef BUFFER_DEPTH\n\
#define BUFFER_DEPTH LINE_NUMBER+1\n\
__kernel void fftKernel_16_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out){\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 (*tmp_in)[BUFFER_DEPTH];\n\
__local double2 (*tmp_out)[BUFFER_DEPTH];\n\
__local double2 (*tmp)[BUFFER_DEPTH];\n\
\n\
  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t ilt = jl+(il/LINE_NUMBER)*LINE_NUMBER;\n\
  size_t jlt = il%LINE_NUMBER;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + jlt;\n\
  tmp1[ilt][jlt] = ilt < 16 ? psi[jgt + ( ilt ) * ndat] : 0.0;\n\
  \n\
  unsigned int A,B;\n\
  unsigned int j;\n\
  A=1;\n\
  B=16;\n\
  tmp_in=tmp1;\n\
  tmp_out=tmp2;\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  #pragma unroll\n\
  for(j=0;j<4;j++){\n\
    B /= 2;\n\
    radix2m(jlt, ilt, 16, A, B, tmp_in, tmp_out);\n\
    tmp = tmp_in; tmp_in = tmp_out; tmp_out = tmp;\n\
    A *= 2;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
  }\n\
\n\
  if(il<16)\n\
    out[jg*n+il] = tmp_in[il][jl];\n\
}\n\
__kernel void fftKernel_15_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out){\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
\n\
  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t ilt = jl+(il/LINE_NUMBER)*LINE_NUMBER;\n\
  size_t jlt = il%LINE_NUMBER;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + jlt;\n\
  tmp1[ilt][jlt] = ilt < 15 ? psi[jgt + ( ilt ) * ndat] : 0.0;\n\
  \n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix3m(jlt, ilt, 15, 1, 5, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix5m(jlt, ilt, 15, 3, 1, tmp2, tmp1);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
  if(il<15)\n\
    out[jg*n+il] = tmp1[il][jl];\n\
}\n\
#undef FFT_LENGTH\n\
#define FFT_LENGTH 32\n\
#undef LINE_NUMBER\n\
#define LINE_NUMBER 4\n\
#undef BUFFER_DEPTH\n\
#define BUFFER_DEPTH LINE_NUMBER+1\n\
__kernel void fftKernel_18_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out){\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
\n\
  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t ilt = jl+(il/LINE_NUMBER)*LINE_NUMBER;\n\
  size_t jlt = il%LINE_NUMBER;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + jlt;\n\
  tmp1[ilt][jlt] = ilt < 18 ? psi[jgt + ( ilt ) * ndat] : 0.0;\n\
  \n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, 18, 1, 9, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix3m(jlt, ilt, 18, 2, 3, tmp2, tmp1);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix3m(jlt, ilt, 18, 6, 1, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
  if(il<18)\n\
    out[jg*n+il] = tmp2[il][jl];\n\
}\n\
__kernel void fftKernel_20_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out){\n\
  const uint fftl=20;\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
\n\
  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t ilt = jl+(il/LINE_NUMBER)*LINE_NUMBER;\n\
  size_t jlt = il%LINE_NUMBER;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + jlt;\n\
  tmp1[ilt][jlt] = ilt < fftl ? psi[jgt + ( ilt ) * ndat] : 0.0;\n\
  \n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 1, 10, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 2, 5, tmp2, tmp1);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix5m(jlt, ilt, fftl, 4, 1, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
  if(il<fftl)\n\
    out[jg*n+il] = tmp2[il][jl];\n\
}\n\
__kernel void fftKernel_24_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out){\n\
  const uint fftl=24;\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
\n\
  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t ilt = jl+(il/LINE_NUMBER)*LINE_NUMBER;\n\
  size_t jlt = il%LINE_NUMBER;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + jlt;\n\
  tmp1[ilt][jlt] = ilt < fftl ? psi[jgt + ( ilt ) * ndat] : 0.0;\n\
  \n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 1, 12, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 2, 6, tmp2, tmp1);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 4, 3, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix3m(jlt, ilt, fftl, 8, 1, tmp2, tmp1);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
  if(il<fftl)\n\
    out[jg*n+il] = tmp1[il][jl];\n\
}\n\
__kernel void fftKernel_30_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out){\n\
  const uint fftl=30;\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
\n\
  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t ilt = jl+(il/LINE_NUMBER)*LINE_NUMBER;\n\
  size_t jlt = il%LINE_NUMBER;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + jlt;\n\
  tmp1[ilt][jlt] = ilt < fftl ? psi[jgt + ( ilt ) * ndat] : 0.0;\n\
  \n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 1, 15, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix3m(jlt, ilt, fftl, 2, 5, tmp2, tmp1);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix5m(jlt, ilt, fftl, 6, 1, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
  if(il<fftl)\n\
    out[jg*n+il] = tmp2[il][jl];\n\
}\n\
__kernel void fftKernel_32_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out){\n\
  const uint fftl=32;\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
\n\
  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t ilt = jl+(il/LINE_NUMBER)*LINE_NUMBER;\n\
  size_t jlt = il%LINE_NUMBER;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + jlt;\n\
  tmp1[ilt][jlt] = ilt < fftl ? psi[jgt + ( ilt ) * ndat] : 0.0;\n\
  \n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 1, 16, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 2, 8, tmp2, tmp1);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 4, 4, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 8, 2, tmp2, tmp1);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 16, 1, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
  if(il<fftl)\n\
    out[jg*n+il] = tmp2[il][jl];\n\
}\n\
#undef FFT_LENGTH\n\
#define FFT_LENGTH 48\n\
#undef LINE_NUMBER\n\
#define LINE_NUMBER 2\n\
#undef BUFFER_DEPTH\n\
#define BUFFER_DEPTH LINE_NUMBER+1\n\
__kernel void fftKernel_48_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out){\n\
  const uint fftl=48;\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
\n\
  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t ilt = jl+(il/LINE_NUMBER)*LINE_NUMBER;\n\
  size_t jlt = il%LINE_NUMBER;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + jlt;\n\
  tmp1[ilt][jlt] = ilt < fftl ? psi[jgt + ( ilt ) * ndat] : 0.0;\n\
  \n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 1, 24, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 2, 12, tmp2, tmp1);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 4, 6, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 8, 3, tmp2, tmp1);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix3m(jlt, ilt, fftl, 16, 1, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
  if(il<fftl)\n\
    out[jg*n+il] = tmp2[il][jl];\n\
}\n\
#undef FFT_LENGTH\n\
#define FFT_LENGTH 64\n\
#undef LINE_NUMBER\n\
#define LINE_NUMBER 2\n\
#undef BUFFER_DEPTH\n\
#define BUFFER_DEPTH LINE_NUMBER+1\n\
__kernel void fftKernel_64_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out){\n\
  const uint fftl=64;\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
\n\
  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t ilt = jl+(il/LINE_NUMBER)*LINE_NUMBER;\n\
  size_t jlt = il%LINE_NUMBER;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + jlt;\n\
  tmp1[ilt][jlt] = ilt < fftl ? psi[jgt + ( ilt ) * ndat] : 0.0;\n\
  \n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 1, 32, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 2, 16, tmp2, tmp1);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 4, 8, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 8, 4, tmp2, tmp1);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 16, 2, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 32, 1, tmp2, tmp1);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
  if(il<fftl)\n\
    out[jg*n+il] = tmp1[il][jl];\n\
}\n\
#undef FFT_LENGTH\n\
#define FFT_LENGTH 128\n\
#undef LINE_NUMBER\n\
#define LINE_NUMBER 1\n\
#undef BUFFER_DEPTH\n\
#define BUFFER_DEPTH LINE_NUMBER\n\
__kernel void fftKernel_128_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out){\n\
  const uint fftl=128;\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
\n\
  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t ilt = il;//jl+(il/LINE_NUMBER)*LINE_NUMBER;\n\
  size_t jlt = 0;//il%LINE_NUMBER;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + jlt;\n\
  tmp1[ilt][jlt] = ilt < fftl ? psi[jgt + ( ilt ) * ndat] : 0.0;\n\
  \n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 1, 64, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 2, 32, tmp2, tmp1);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 4, 16, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 8, 8, tmp2, tmp1);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 16, 4, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 32, 2, tmp2, tmp1);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  radix2m(jlt, ilt, fftl, 64, 1, tmp1, tmp2);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
  if(il<fftl)\n\
    out[jg*n+il] = tmp2[il][jl];\n\
}\n\
";


inline void fft_test_2_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    int FFT_LENGTH=16;
    int LINE_NUMBER=16;
    assert(*n==FFT_LENGTH);
    size_t block_size_i=FFT_LENGTH, block_size_j=LINE_NUMBER/2;
    cl_uint i = 0;
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j*2,*ndat)/2};
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue fft kernel!");
}


inline void fft_test_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    int FFT_LENGTH=16;
    int LINE_NUMBER=16;
    assert(*n==FFT_LENGTH);
    size_t block_size_i=FFT_LENGTH, block_size_j=LINE_NUMBER;
    cl_uint i = 0;
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue fft kernel!");
}

inline void fft_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    int FFT_LENGTH=128;
    int LINE_NUMBER=1;
//    assert(*n==FFT_LENGTH);
    size_t block_size_i=FFT_LENGTH, block_size_j=LINE_NUMBER;
    cl_uint i = 0;
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue fft kernel!");
}

inline void fft_generated_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    cl_uint shared_size_used=0;
    if(*n <= 64)
       shared_size_used=512;
    else if(*n <= 256)
       shared_size_used=1024;
    int FFT_LENGTH=(*n / 16) * 16 + (*n % 16 ? 16 : 0);
    int LINE_LIMIT=(shared_size_used/FFT_LENGTH) & (~3);
    int LINE_NUMBER = 1;
    while( LINE_LIMIT /= 2 )
      LINE_NUMBER *= 2;
    if( LINE_NUMBER > FFT_LENGTH )
       LINE_NUMBER = FFT_LENGTH;
    if( LINE_NUMBER < 1 )
       LINE_NUMBER=1;
    size_t block_size_i=FFT_LENGTH, block_size_j=LINE_NUMBER;
    cl_uint i = 0;
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue fft kernel!");
}

void FC_FUNC_(customize_fft,CUSTOMIZE_FFT)(cl_uint *fft) {
  fft_size = *fft;
}

void FC_FUNC_(fft1d_d,FFT1D_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    if(fft_size!=0)
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d, (*command_queue)->command_queue, n, ndat, psi, out);
    else
      fft_generic((*command_queue)->kernels.fft_kernel_d, (*command_queue)->command_queue, n, ndat, psi, out);
}

void FC_FUNC_(fft3d_d,FFT1D_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,cl_mem *psi,cl_mem *out,cl_mem *tmp){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    if(fft_size!=0)
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, out);
    else
      fft_generic((*command_queue)->kernels.fft_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, out);
    ndat = n1 * n3;
    if(fft_size!=0)
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, tmp);
    else
      fft_generic((*command_queue)->kernels.fft_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, tmp);
    ndat = n2 * n3;
    if(fft_size!=0)
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d, (*command_queue)->command_queue, &n1, &ndat, tmp, out);
    else
      fft_generic((*command_queue)->kernels.fft_kernel_d, (*command_queue)->command_queue, &n1, &ndat, tmp, out);
}

cl_program fftProgram;

cl_uint fft_size=0;

void create_fft_kernels(struct bigdft_kernels * kernels){
    cl_int ciErrNum = CL_SUCCESS;
    kernels->fft_kernel_d=clCreateKernel(fftProgram,"fftKernel_128_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d kernel!");
}

void create_fft_generated_kernels(struct bigdft_kernels * kernels, cl_uint fft){
    char kernel_name[256];
    sprintf(kernel_name,"fftKernel_%u_d",fft);
    cl_int ciErrNum = CL_SUCCESS;
    kernels->fft_kernel_d=clCreateKernel(fftProgram,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d kernel!");
}


void build_fft_programs(cl_context * context){
    fft_size=0;
    cl_int ciErrNum = CL_SUCCESS;
    fftProgram = clCreateProgramWithSource(*context,1,(const char**) &fft_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(fftProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to build fft program!\n",ciErrNum);
        char cBuildLog[10240];
        clGetProgramBuildInfo(fftProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
}

void build_fft_generated_programs(cl_context * context, cl_uint fft){
    cl_int ciErrNum = CL_SUCCESS;
    char * code;
    code = generate_fft_program(fft);
    printf("%s",code);
    fftProgram = clCreateProgramWithSource(*context,1,(const char**) &code, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(fftProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to build fft program!\n",ciErrNum);
        char cBuildLog[10240];
        clGetProgramBuildInfo(fftProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
}


void clean_fft_kernels(struct bigdft_kernels * kernels){
  cl_int ciErrNum;
  ciErrNum = clReleaseKernel(kernels->fft_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
}

void clean_fft_programs(){
  cl_int ciErrNum;
  ciErrNum = clReleaseProgram(fftProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
}
