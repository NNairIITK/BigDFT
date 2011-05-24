//! @file
/*!
  Kernels in the ana_program peform a wavelet analysis,
  which can be reduced to a convolution and a transposition.
  In the kernels, each work item processes 2 elements, and
  loads 3 elements.
  The filter macro is optimized to reduce register usage
  in order to imporve occupancy.
  Buffer is of size 16*49.
  Memory pattern somewhat ressembles this :
  ///////   ///////
  \\\\\\\   ///////
  ///////   \\\\\\\
  \\\\\\\   \\\\\\\
  Even line centered convolutions go int the upper part
  of the resulting matrix, while odd line centered convolutions
  go to the lower part of the resulting matrix.
  Size  of the data is 2*n * ndat.
*/
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include "Wavelet.h"
#include "OpenCL_wrappers.h"

char * ana_program="\
#ifdef cl_khr_fp64\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#elif defined (cl_amd_fp64)\n\
#pragma OPENCL EXTENSION cl_amd_fp64: enable \n\
#endif\n\
#define FILTER_WIDTH 16\n\
#define filter(ci,di,tmp)\
ci += tmp[14] * -0.00030292051472413308126;\
di += tmp[14] *  0.00054213233180001068935;\
ci += tmp[ 1] * -0.00054213233180001068935;\
di += tmp[ 1] * -0.00030292051472413308126;\
ci += tmp[15] *  0.0018899503327676891843;\
di += tmp[15] * -0.0033824159510050025955;\
ci += tmp[ 0] * -0.0033824159510050025955;\
di += tmp[ 0] * -0.0018899503327676891843;\
ci += tmp[12] *  0.0038087520138944894631;\
di += tmp[12] * -0.0076074873249766081919;\
ci += tmp[ 3] *  0.0076074873249766081919;\
di += tmp[ 3] *  0.0038087520138944894631;\
ci += tmp[13] * -0.014952258337062199118;\
di += tmp[13] *  0.031695087811525991431;\
ci += tmp[ 2] *  0.031695087811525991431;\
di += tmp[ 2] *  0.014952258337062199118;\
ci += tmp[10] * -0.027219029917103486322;\
di += tmp[10] *  0.061273359067811077843;\
ci += tmp[ 5] * -0.061273359067811077843;\
di += tmp[ 5] * -0.027219029917103486322;\
ci += tmp[11] *  0.049137179673730286787;\
di += tmp[11] * -0.14329423835127266284;\
ci += tmp[ 4] * -0.14329423835127266284;\
di += tmp[ 4] * -0.049137179673730286787;\
ci += tmp[ 9] * -0.051945838107881800736;\
di += tmp[ 9] *  0.48135965125905339159;\
ci += tmp[ 6] *  0.48135965125905339159;\
di += tmp[ 6] *  0.051945838107881800736;\
ci += tmp[ 8] *  0.36444189483617893676;\
di += tmp[ 8] * -0.77718575169962802862;\
ci += tmp[ 7] *  0.77718575169962802862;\
di += tmp[ 7] *  0.36444189483617893676;\n\
//periodic boundary condition of the filter\n\
__kernel void ana1dKernel_d(uint n, uint ndat, __global const double *psi, __global double *out){\n\
__local double tmp1[FILTER_WIDTH*(3*FILTER_WIDTH+1)];\n\
__local double *tmp = &tmp1[0];\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1);\n\
ptrdiff_t igt = get_group_id(0);\n\
ptrdiff_t jgt = get_group_id(1);\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
igt = (ig - i2) * 2 + j2 - FILTER_WIDTH/2;\n\
jgt = jg - j2 + i2;\n\
//if we are outside the matrix, load elements wrapping aroud.\n\
if ( igt < 0 ) \n\
  tmp[i2 * (3 * FILTER_WIDTH + 1) + j2] = psi[jgt + ( 2*n + igt ) * ndat];\n\
else \n\
  tmp[i2 * (3 * FILTER_WIDTH + 1) + j2] = psi[jgt + igt * ndat];\n\
igt += FILTER_WIDTH;\n\
tmp[i2 * (3 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = psi[jgt + igt * ndat];\n\
igt += FILTER_WIDTH;\n\
if ( igt >= 2*n ) \n\
  tmp[i2 * (3 * FILTER_WIDTH + 1) + j2 + 2*FILTER_WIDTH] = psi[jgt + ( igt - 2*n ) * ndat];\n\
else\n\
  tmp[i2 * (3 * FILTER_WIDTH + 1) + j2 + 2*FILTER_WIDTH] = psi[jgt +  igt * ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
double ci = 0.0;\n\
double di = 0.0;\n\
tmp = tmp + j2*(3*FILTER_WIDTH+1) + 2*i2 + 1;\n\
filter(ci,di,tmp);\n\
\
out[(jg*(2*n)+ig)]=ci;\n\
out[(jg*(2*n)+ig+n)]=di;\n\
};\n\
//non periodic boundary condition version of the filter\n\
//output data is shrinked in regard of the input data\n\
__kernel void anashrink1dKernel_d(uint n, uint ndat, __global const double *psi, __global double *out){\n\
__local double tmp1[FILTER_WIDTH*(3*FILTER_WIDTH+1)];\n\
__local double *tmp = &tmp1[0];\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1);\n\
ptrdiff_t igt = get_group_id(0);\n\
ptrdiff_t jgt = get_group_id(1);\n\
ptrdiff_t jb;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
igt = 2*(ig - i2) + j2;\n\
jgt = jg - j2 + i2;\n\
psi += 7*ndat;\n\
igt -= 7;\n\
//shrinking, so no elements nedd to be loaded outside the matrix.\n\
tmp[i2 * (3 * FILTER_WIDTH + 1) + j2]=psi[jgt+igt*ndat];\n\
igt += FILTER_WIDTH - 1;\n\
tmp[i2 * (3 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH - 1]=psi[jgt+igt*ndat];\n\
igt += FILTER_WIDTH - 1;\n\
tmp[i2 * (3 * FILTER_WIDTH + 1) + j2 + 2 * FILTER_WIDTH - 2]=psi[jgt+igt*ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
double ci = 0.0;\n\
double di = 0.0;\n\
tmp = tmp + j2*(3*FILTER_WIDTH+1) + 2*i2;\n\
filter(ci,di,tmp);\n\
\
out[(jg*(2*n)+ig)]=ci;\n\
out[(jg*(2*n)+ig+n)]=di;\n\
};\n\
#define ELEM_PER_THREAD 2\n\
__kernel void ana1d_blockKernel_d(uint n, uint ndat, __global const double *psi, __global double *out, __local double *tmp){\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1)*ELEM_PER_THREAD;\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1)*ELEM_PER_THREAD;\n\
ptrdiff_t igt = get_group_id(0);\n\
ptrdiff_t jgt = get_group_id(1);\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1)*ELEM_PER_THREAD - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
igt = (ig - i2) * 2 + j2 - FILTER_WIDTH/2;\n\
jgt = jg - j2 + i2;\n\
__local double * tmp_1 = tmp + i2 * (3 * FILTER_WIDTH + 1) + j2;\n\
psi += jgt;\n\
if ( igt < 0 ) \n\
  *tmp_1++ = psi[( 2*n + igt++ ) * ndat];\n\
else \n\
  *tmp_1++ = psi[igt++ * ndat];\n\
if ( igt < 0 ) \n\
  *tmp_1++ = psi[( 2*n + igt++ ) * ndat];\n\
else \n\
  *tmp_1++ = psi[igt++ * ndat];\n\
igt += FILTER_WIDTH-ELEM_PER_THREAD;\n\
tmp_1 += FILTER_WIDTH-ELEM_PER_THREAD;\n\
*tmp_1++ = psi[igt++ * ndat];\n\
*tmp_1++ = psi[igt++ * ndat];\n\
igt += FILTER_WIDTH-ELEM_PER_THREAD;\n\
tmp_1 += FILTER_WIDTH-ELEM_PER_THREAD;\n\
if ( igt >= 2*n ) \n\
  *tmp_1++ = psi[(igt++ - 2*n) * ndat];\n\
else\n\
  *tmp_1++ = psi[igt++ * ndat];\n\
if ( igt >= 2*n ) \n\
  *tmp_1++ = psi[(igt++ - 2*n) * ndat];\n\
else\n\
  *tmp_1++ = psi[igt++ * ndat];\n\
\
double ci_1 = 0.0;\n\
double di_1 = 0.0;\n\
double ci_2 = 0.0;\n\
double di_2 = 0.0;\n\
tmp_1 = tmp + j2*(3*FILTER_WIDTH+1) + 2*i2 + 1;\n\
__local double * tmp_2 = tmp_1 + (3*FILTER_WIDTH+1);\n\
out += (jg*(2*n)+ig);\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
filter(ci_1,di_1,tmp_1);\n\
filter(ci_2,di_2,tmp_2);\n\
\
*out = ci_1;\n\
out += n;\n\
*out = di_1;\n\
out += n;\n\
*out = ci_2;\n\
out += n;\n\
*out = di_2;\n\
};\n\
";


/*
  Synthesis is the reciprocal operation of analysis.
  Each work item processes 2 elements.
  Buffer is of size 16*49.
  memory pattern is obviously reversed :
  
  ///////    ///////
  ///////    \\\\\\\
  \\\\\\\    ///////
  \\\\\\\    \\\\\\\
  As buffer is split into 2 zones, each of length:
  4+16+4, some work items load 4 elements, while other
  load 2 elements. The performance inpact of this is minimal.
  Size  of the data is 2*n * ndat.
*/
char * syn_program="\
#ifdef cl_khr_fp64\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#elif defined (cl_amd_fp64)\n\
#pragma OPENCL EXTENSION cl_amd_fp64: enable \n\
#endif\n\
#define FILTER_WIDTH 8\n\
#define SIZE_I 16\n\
__kernel void syn1dKernel_d(uint n, uint ndat, __global const double *psi, __global double *out){\n\
__local double tmp1[SIZE_I*(2*FILTER_WIDTH+2*SIZE_I+1)];\n\
__local double *tmp = &tmp1[0];\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1);\n\
ptrdiff_t igt = get_group_id(0);\n\
ptrdiff_t jgt = get_group_id(1);\n\
size_t ioff;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
igt = ig - i2 + j2 - FILTER_WIDTH/2;\n\
jgt = jg - j2 + i2;\n\
//If I'm on the outside, select a border element to load\n\
ioff = i2*(2*FILTER_WIDTH+2*SIZE_I+1) + j2;\n\
if( igt < 0 ) {\n\
  tmp[ioff] = psi[jgt+(n+igt)*ndat];\n\
  tmp[ioff+FILTER_WIDTH+SIZE_I] = psi[jgt+(n+igt+n)*ndat];\n\
} else {\n\
  tmp[ioff] = psi[jgt+igt*ndat];\n\
  tmp[ioff+FILTER_WIDTH+SIZE_I] = psi[jgt+(igt+n)*ndat];\n\
}\n\
igt += SIZE_I;\n\
//only half the work items load these elements\n\
if( j2 < SIZE_I - FILTER_WIDTH){\n\
  if ( igt >=n ) {\n\
    tmp[ioff+SIZE_I] = psi[jgt+(igt-n)*ndat];\n\
    tmp[ioff+FILTER_WIDTH+2*SIZE_I] = psi[jgt+(igt-n+n)*ndat];\n\
  } else {\n\
    tmp[ioff+SIZE_I] = psi[jgt+igt*ndat];\n\
    tmp[ioff+FILTER_WIDTH+2*SIZE_I] = psi[jgt+(igt+n)*ndat];\n\
  }\n\
}\n\
tmp += j2*(2*FILTER_WIDTH+2*SIZE_I+1) + FILTER_WIDTH/2+i2;\n\
double se = 0.0;\n\
double so = 0.0;\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
se += tmp[-3] * -0.014952258337062199118;\n\
so += tmp[-3] * -0.00030292051472413308126;\n\
se += tmp[FILTER_WIDTH+SIZE_I + 3] * -0.00030292051472413308126;\n\
so += tmp[FILTER_WIDTH+SIZE_I + 3] *  0.014952258337062199118;\n\
se += tmp[-2] *  0.049137179673730286787;\n\
so += tmp[-2] *  0.0038087520138944894631;\n\
se += tmp[FILTER_WIDTH+SIZE_I + 2] *  0.0038087520138944894631;\n\
so += tmp[FILTER_WIDTH+SIZE_I + 2] * -0.049137179673730286787;\n\
se += tmp[-1] * -0.051945838107881800736;\n\
so += tmp[-1] * -0.027219029917103486322;\n\
se += tmp[FILTER_WIDTH+SIZE_I + 1] * -0.027219029917103486322;\n\
so += tmp[FILTER_WIDTH+SIZE_I + 1] *  0.051945838107881800736;\n\
se += tmp[ 0] *  0.77718575169962802862;\n\
so += tmp[ 0] *  0.36444189483617893676;\n\
se += tmp[FILTER_WIDTH+SIZE_I + 0] *  0.36444189483617893676;\n\
so += tmp[FILTER_WIDTH+SIZE_I + 0] * -0.77718575169962802862;\n\
se += tmp[ 1] * -0.061273359067811077843;\n\
so += tmp[ 1] *  0.48135965125905339159;\n\
se += tmp[FILTER_WIDTH+SIZE_I - 1] *  0.48135965125905339159;\n\
so += tmp[FILTER_WIDTH+SIZE_I - 1] *  0.061273359067811077843;\n\
se += tmp[ 2] *  0.0076074873249766081919;\n\
so += tmp[ 2] * -0.14329423835127266284;\n\
se += tmp[FILTER_WIDTH+SIZE_I - 2] * -0.14329423835127266284;\n\
so += tmp[FILTER_WIDTH+SIZE_I - 2] * -0.0076074873249766081919;\n\
se += tmp[ 3] * -0.00054213233180001068935;\n\
so += tmp[ 3] *  0.031695087811525991431;\n\
se += tmp[FILTER_WIDTH+SIZE_I - 3] *  0.031695087811525991431;\n\
so += tmp[FILTER_WIDTH+SIZE_I - 3] *  0.00054213233180001068935;\n\
se += tmp[-4] *  0.0018899503327676891843;\n\
so += tmp[ 4] * -0.0033824159510050025955;\n\
se += tmp[FILTER_WIDTH+SIZE_I - 4] * -0.0033824159510050025955;\n\
so += tmp[FILTER_WIDTH+SIZE_I + 4] * -0.0018899503327676891843;\n\
\
out[jg*(2*n)+ig*2]=se;\n\
out[jg*(2*n)+ig*2+1]=so;\n\
};\n\
__kernel void syngrow1dKernel_d(uint n, uint ndat, __global const double *psi, __global double *out){\n\
__local double tmp1[SIZE_I*(2*FILTER_WIDTH+2*SIZE_I+1)];\n\
__local double *tmp = &tmp1[0];\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1);\n\
ptrdiff_t igt = get_group_id(0);\n\
ptrdiff_t jgt = get_group_id(1);\n\
ptrdiff_t jb;\n\
size_t ioff;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
igt = ig - i2 + j2 - FILTER_WIDTH/2 - 4;\n\
jgt = jg - j2 + i2;\n\
//If I'm on the outside, select a border element to load\n\
ioff = i2*(2*FILTER_WIDTH+2*SIZE_I+1) + j2;\n\
if( igt < 0 ) {\n\
  tmp[ioff] = 0.0;\n\
  tmp[ioff+FILTER_WIDTH+SIZE_I] = 0.0;\n\
} else {\n\
  tmp[ioff] = psi[jgt+igt*ndat];\n\
  tmp[ioff+FILTER_WIDTH+SIZE_I] = psi[jgt+(igt+n-7)*ndat];\n\
}\n\
igt += SIZE_I;\n\
if( j2 < SIZE_I - FILTER_WIDTH){\n\
  if ( igt >=n-7 ) {\n\
    tmp[ioff+SIZE_I] = 0.0;\n\
    tmp[ioff+FILTER_WIDTH+2*SIZE_I] = 0.0;\n\
  } else {\n\
    tmp[ioff+SIZE_I] = psi[jgt+igt*ndat];\n\
    tmp[ioff+FILTER_WIDTH+2*SIZE_I] = psi[jgt+(igt+n-7)*ndat];\n\
  }\n\
}\n\
tmp += j2*(2*FILTER_WIDTH+2*SIZE_I+1) + FILTER_WIDTH/2+i2;\n\
double se = 0.0;\n\
double so = 0.0;\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
so += tmp[FILTER_WIDTH+SIZE_I + 4] * -0.0018899503327676891843;\n\
se += tmp[FILTER_WIDTH+SIZE_I + 4] * -0.00030292051472413308126;\n\
so += tmp[-3] * -0.00030292051472413308126;\n\
se += tmp[-3] *  0.0018899503327676891843;\n\
so += tmp[FILTER_WIDTH+SIZE_I + 3] *  0.014952258337062199118;\n\
se += tmp[FILTER_WIDTH+SIZE_I + 3] *  0.0038087520138944894631;\n\
so += tmp[-2] *  0.0038087520138944894631;\n\
se += tmp[-2] * -0.014952258337062199118;\n\
so += tmp[FILTER_WIDTH+SIZE_I + 2] * -0.049137179673730286787;\n\
se += tmp[FILTER_WIDTH+SIZE_I + 2] * -0.027219029917103486322;\n\
so += tmp[-1] * -0.027219029917103486322;\n\
se += tmp[-1] *  0.049137179673730286787;\n\
so += tmp[FILTER_WIDTH+SIZE_I + 1] *  0.051945838107881800736;\n\
se += tmp[FILTER_WIDTH+SIZE_I + 1] *  0.36444189483617893676;\n\
so += tmp[ 0] *  0.36444189483617893676;\n\
se += tmp[ 0] * -0.051945838107881800736;\n\
so += tmp[FILTER_WIDTH+SIZE_I + 0] * -0.77718575169962802862;\n\
se += tmp[FILTER_WIDTH+SIZE_I + 0] *  0.48135965125905339159;\n\
so += tmp[ 1] *  0.48135965125905339159;\n\
se += tmp[ 1] *  0.77718575169962802862;\n\
so += tmp[FILTER_WIDTH+SIZE_I - 1] *  0.061273359067811077843;\n\
se += tmp[FILTER_WIDTH+SIZE_I - 1] * -0.14329423835127266284;\n\
so += tmp[ 2] * -0.14329423835127266284;\n\
se += tmp[ 2] * -0.061273359067811077843;\n\
so += tmp[FILTER_WIDTH+SIZE_I - 2] * -0.0076074873249766081919;\n\
se += tmp[FILTER_WIDTH+SIZE_I - 2] *  0.031695087811525991431;\n\
so += tmp[ 3] *  0.031695087811525991431;\n\
se += tmp[ 3] *  0.0076074873249766081919;\n\
so += tmp[FILTER_WIDTH+SIZE_I - 3] *  0.00054213233180001068935;\n\
se += tmp[FILTER_WIDTH+SIZE_I - 3] * -0.0033824159510050025955;\n\
so += tmp[ 4] * -0.0033824159510050025955;\n\
se += tmp[ 4] * -0.00054213233180001068935;\n\
\
out[jg*2*n+2*ig]=so;\n\
out[jg*2*n+2*ig+1]=se;\n\
};\n\
";

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

cl_program anaProgram;
cl_program synProgram;

void create_wavelet_kernels(struct bigdft_kernels * kernels) {
    cl_int ciErrNum = CL_SUCCESS;
    kernels->anashrink1d_kernel_d=clCreateKernel(anaProgram,"anashrink1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create anashrink1dKernel_d kernel!");
    kernels->ana1d_kernel_d=clCreateKernel(anaProgram,"ana1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create ana1dKernel_d kernel!");
    kernels->ana1d_block_kernel_d=clCreateKernel(anaProgram,"ana1d_blockKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create ana1d_blockKernel_d kernel!");
    kernels->syngrow1d_kernel_d=clCreateKernel(synProgram,"syngrow1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create syngrow1dKernel_d kernel!");
    kernels->syn1d_kernel_d=clCreateKernel(synProgram,"syn1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create syn1dKernel_d kernel!");
}

void build_wavelet_programs(cl_context * context){
    cl_int ciErrNum = CL_SUCCESS;
    anaProgram = clCreateProgramWithSource(*context,1,(const char**) &ana_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(anaProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build ana program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(anaProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    synProgram = clCreateProgramWithSource(*context,1,(const char**) &syn_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(synProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build syn program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(synProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
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

void clean_wavelet_programs(){
  cl_int ciErrNum;
  ciErrNum = clReleaseProgram(anaProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
  ciErrNum = clReleaseProgram(synProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
}
