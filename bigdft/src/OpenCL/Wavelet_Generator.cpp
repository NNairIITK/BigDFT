#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include "Wavelet_Generator.h"
#include "OpenCL_wrappers.h"

#define FILTER_WIDTH 16
static void generate_ana_header(std::stringstream &program){
  program<<"#ifdef cl_khr_fp64\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#elif defined (cl_amd_fp64)\n\
#pragma OPENCL EXTENSION cl_amd_fp64: enable \n\
#endif\n\
#define FILTER_WIDTH "<<FILTER_WIDTH<<"\n";
}

static void generate_ana_filter(std::stringstream &program){
  program<<"#define filter(ci,di,tmp)\
ci = mad(tmp[14], -0.00030292051472413308126, ci);\
di = mad(tmp[14],  0.00054213233180001068935, di);\
ci = mad(tmp[ 1], -0.00054213233180001068935, ci);\
di = mad(tmp[ 1], -0.00030292051472413308126, di);\
ci = mad(tmp[15],  0.0018899503327676891843 , ci);\
di = mad(tmp[15], -0.0033824159510050025955 , di);\
ci = mad(tmp[ 0], -0.0033824159510050025955 , ci);\
di = mad(tmp[ 0], -0.0018899503327676891843 , di);\
ci = mad(tmp[12],  0.0038087520138944894631 , ci);\
di = mad(tmp[12], -0.0076074873249766081919 , di);\
ci = mad(tmp[ 3],  0.0076074873249766081919 , ci);\
di = mad(tmp[ 3],  0.0038087520138944894631 , di);\
ci = mad(tmp[13], -0.014952258337062199118  , ci);\
di = mad(tmp[13],  0.031695087811525991431  , di);\
ci = mad(tmp[ 2],  0.031695087811525991431  , ci);\
di = mad(tmp[ 2],  0.014952258337062199118  , di);\
ci = mad(tmp[10], -0.027219029917103486322  , ci);\
di = mad(tmp[10],  0.061273359067811077843  , di);\
ci = mad(tmp[ 5], -0.061273359067811077843  , ci);\
di = mad(tmp[ 5], -0.027219029917103486322  , di);\
ci = mad(tmp[11],  0.049137179673730286787  , ci);\
di = mad(tmp[11], -0.14329423835127266284   , di);\
ci = mad(tmp[ 4], -0.14329423835127266284   , ci);\
di = mad(tmp[ 4], -0.049137179673730286787  , di);\
ci = mad(tmp[ 9], -0.051945838107881800736  , ci);\
di = mad(tmp[ 9],  0.48135965125905339159   , di);\
ci = mad(tmp[ 6],  0.48135965125905339159   , ci);\
di = mad(tmp[ 6],  0.051945838107881800736  , di);\
ci = mad(tmp[ 8],  0.36444189483617893676   , ci);\
di = mad(tmp[ 8], -0.77718575169962802862   , di);\
ci = mad(tmp[ 7],  0.77718575169962802862   , ci);\
di = mad(tmp[ 7],  0.36444189483617893676   , di);\n\
#define filter_vector2(tmp) \
double2 ci = (double2)(0.0, 0.0);\
double2 di = (double2)(0.0, 0.0);\
__local double2 *tmp2= (__local double2 *)tmp;\
ci = mad( *tmp2,   (double2)(-0.0033824159510050025955  ,-0.00054213233180001068935), ci);\
di = mad( *tmp2++, (double2)(-0.0018899503327676891843  ,-0.00030292051472413308126), di);\
ci = mad( *tmp2,   (double2)( 0.031695087811525991431   , 0.0076074873249766081919) , ci);\
di = mad( *tmp2++, (double2)( 0.014952258337062199118   , 0.0038087520138944894631) , di);\
ci = mad( *tmp2,   (double2)(-0.14329423835127266284    ,-0.061273359067811077843)  , ci);\
di = mad( *tmp2++, (double2)(-0.049137179673730286787   ,-0.027219029917103486322)  , di);\
ci = mad( *tmp2,   (double2)( 0.48135965125905339159    , 0.77718575169962802862)   , ci);\
di = mad( *tmp2++, (double2)( 0.051945838107881800736   , 0.36444189483617893676)   , di);\
ci = mad( *tmp2,   (double2)( 0.36444189483617893676    ,-0.051945838107881800736)  , ci);\
di = mad( *tmp2++, (double2)(-0.77718575169962802862    , 0.48135965125905339159)   , di);\
ci = mad( *tmp2,   (double2)(-0.027219029917103486322   , 0.049137179673730286787)  , ci);\
di = mad( *tmp2++, (double2)( 0.061273359067811077843   ,-0.14329423835127266284)   , di);\
ci = mad( *tmp2,   (double2)( 0.0038087520138944894631  ,-0.014952258337062199118)  , ci);\
di = mad( *tmp2++, (double2)(-0.0076074873249766081919  , 0.031695087811525991431)  , di);\
ci = mad( *tmp2,   (double2)(-0.00030292051472413308126 , 0.0018899503327676891843) , ci);\
di = mad( *tmp2++, (double2)( 0.00054213233180001068935 ,-0.0033824159510050025955) , di);\n";
}

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
static void generate_ana1dKernel(std::stringstream &program, struct bigdft_device_infos * infos){
  program<<"//periodic boundary condition of the filter\n\
__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<", 1))) void ana1dKernel_d(uint n, uint ndat, __global const double * restrict psi, __global double * restrict out){\n\
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
barrier(CLK_LOCAL_MEM_FENCE);\n";
  if(strncmp(infos->NAME,"Cayman",strlen("Cayman"))==0){
    program<<"tmp = tmp + j2*(3*FILTER_WIDTH+1) + 2*i2 + 1;\n\
filter_vector2(tmp);\n\
\
out[(jg*(2*n)+ig)]=ci.x+ci.y;\n\
out[(jg*(2*n)+ig+n)]=di.x+di.y;\n\
};\n";
  } else {
    program<<"double ci = 0.0;\n\
double di = 0.0;\n\
tmp = tmp + j2*(3*FILTER_WIDTH+1) + 2*i2 + 1;\n\
filter(ci,di,tmp);\n\
\
out[(jg*(2*n)+ig)]=ci;\n\
out[(jg*(2*n)+ig+n)]=di;\n\
};\n";
  }
}

static void generate_anashrink1dKernel(std::stringstream &program, struct bigdft_device_infos * infos){
  program<<"//non periodic boundary condition version of the filter\n\
//output data is shrinked in regard of the input data\n\
__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<", 1))) void anashrink1dKernel_d(uint n, uint ndat, __global const double * restrict psi, __global double * restrict out){\n\
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
barrier(CLK_LOCAL_MEM_FENCE);\n";
  if(strncmp(infos->NAME,"Cayman",strlen("Cayman"))==0){
    program<<"tmp = tmp + j2*(3*FILTER_WIDTH+1) + 2*i2;\n\
filter_vector2(tmp);\n\
\
out[(jg*(2*n)+ig)]=ci.x+ci.y;\n\
out[(jg*(2*n)+ig+n)]=di.x+di.y;\n\
};\n";
  } else {
    program<<"double ci = 0.0;\n\
double di = 0.0;\n\
tmp = tmp + j2*(3*FILTER_WIDTH+1) + 2*i2;\n\
filter(ci,di,tmp);\n\
\
out[(jg*(2*n)+ig)]=ci;\n\
out[(jg*(2*n)+ig+n)]=di;\n\
};\n";
  }
}

static void generate_ana1d_blockKernel(std::stringstream &program){
  program<<"#define ELEM_PER_THREAD 2\n\
__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<"/ELEM_PER_THREAD, 1))) void ana1d_blockKernel_d(uint n, uint ndat, __global const double *psi, __global double * restrict out, __local double * restrict tmp){\n\
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
};\n";
}

static void generate_syn_header(std::stringstream &program){
  program<<"#ifdef cl_khr_fp64\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#elif defined (cl_amd_fp64)\n\
#pragma OPENCL EXTENSION cl_amd_fp64: enable \n\
#endif\n\
#define FILTER_WIDTH 8\n\
#define SIZE_I 16\n";
}

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
static void generate_syn_filter(std::stringstream &program){
  program<<"#define filter_vector2(tmp1,tmp2) \
double2 se = (double2)(0.0, 0.0);\
double2 so = (double2)(0.0, 0.0);\
se = mad((double2)(*tmp1,*(tmp1+7)), (double2)(0.0018899503327676891843,-0.00054213233180001068935),  se);\
__local double2 *tmp1b = (__local double2 *)(tmp1+1);\
so = mad(*tmp1b,   (double2)(-0.00030292051472413308126, 0.0038087520138944894631) ,so);\
se = mad(*tmp1b++, (double2)(-0.014952258337062199118,   0.049137179673730286787)  ,se);\
so = mad(*tmp1b,   (double2)(-0.027219029917103486322,   0.36444189483617893676)   ,so);\
se = mad(*tmp1b++, (double2)(-0.051945838107881800736,   0.77718575169962802862)   ,se);\
so = mad(*tmp1b,   (double2)( 0.48135965125905339159,   -0.14329423835127266284)   ,so);\
se = mad(*tmp1b++, (double2)(-0.061273359067811077843,   0.0076074873249766081919) ,se);\
so = mad(*tmp1b,   (double2)( 0.031695087811525991431,  -0.0033824159510050025955) ,so);\
se = mad((double2)(*tmp2,*(tmp2+7)), (double2)(-0.0033824159510050025955, -0.00030292051472413308126),  se);\
__local double2 *tmp2b = (__local double2 *)(tmp2+1);\
so = mad(*tmp2b,   (double2)( 0.00054213233180001068935,-0.0076074873249766081919) ,so);\
se = mad(*tmp2b++, (double2)( 0.031695087811525991431,  -0.14329423835127266284)   ,se);\
so = mad(*tmp2b,   (double2)( 0.061273359067811077843,  -0.77718575169962802862)   ,so);\
se = mad(*tmp2b++, (double2)( 0.48135965125905339159,    0.36444189483617893676)   ,se);\
so = mad(*tmp2b,   (double2)( 0.051945838107881800736,  -0.049137179673730286787)  ,so);\
se = mad(*tmp2b++, (double2)(-0.027219029917103486322,   0.0038087520138944894631) ,se);\
so = mad(*tmp2b,   (double2)( 0.014952258337062199118,  -0.0018899503327676891843) ,so);\n\
#define filter(tmp1,tmp2) \
double se = 0.0;\
double so = 0.0;\
__local double *tmp1b = tmp1;\
__local double *tmp2b = tmp2+8;\
se = mad(*tmp1b++,  0.0018899503327676891843,  se);\
se = mad(*tmp1b,   -0.014952258337062199118,   se);\
so = mad(*tmp1b++, -0.00030292051472413308126, so);\
so = mad(*tmp2b--, -0.0018899503327676891843,  so);\
se = mad(*tmp2b,   -0.00030292051472413308126, se);\
so = mad(*tmp2b--,  0.014952258337062199118,   so);\
se = mad(*tmp1b,    0.049137179673730286787,   se);\
so = mad(*tmp1b++,  0.0038087520138944894631,  so);\
se = mad(*tmp2b,    0.0038087520138944894631,  se);\
so = mad(*tmp2b--, -0.049137179673730286787,   so);\
se = mad(*tmp1b,   -0.051945838107881800736,   se);\
so = mad(*tmp1b++, -0.027219029917103486322,   so);\
se = mad(*tmp2b,   -0.027219029917103486322,   se);\
so = mad(*tmp2b--,  0.051945838107881800736,   so);\
se = mad(*tmp1b,    0.77718575169962802862,    se);\
so = mad(*tmp1b++,  0.36444189483617893676,    so);\
se = mad(*tmp2b,    0.36444189483617893676,    se);\
so = mad(*tmp2b--, -0.77718575169962802862,    so);\
se = mad(*tmp1b,   -0.061273359067811077843,   se);\
so = mad(*tmp1b++,  0.48135965125905339159,    so);\
se = mad(*tmp2b,    0.48135965125905339159,    se);\
so = mad(*tmp2b--,  0.061273359067811077843,   so);\
se = mad(*tmp1b,    0.0076074873249766081919,  se);\
so = mad(*tmp1b++, -0.14329423835127266284,    so);\
se = mad(*tmp2b,   -0.14329423835127266284,    se);\
so = mad(*tmp2b--, -0.0076074873249766081919,  so);\
se = mad(*tmp1b,   -0.00054213233180001068935, se);\
so = mad(*tmp1b++,  0.031695087811525991431,   so);\
se = mad(*tmp2b,    0.031695087811525991431,   se);\
so = mad(*tmp2b--,  0.00054213233180001068935, so);\
so = mad(*tmp1b++, -0.0033824159510050025955,  so);\
se = mad(*tmp2b--, -0.0033824159510050025955,  se);\n\
#define filter_grow_vector2(tmp1,tmp2) \
double2 so = (double2)(0.0, 0.0);\
double2 se = (double2)(0.0, 0.0);\
__local double2 *tmp1b = (__local double2 *)tmp1;\
__local double2 *tmp2b = (__local double2 *)tmp2;\
so = mad(*tmp1b,   (double2)(-0.00030292051472413308126, 0.0038087520138944894631) , so);\
se = mad(*tmp1b++, (double2)( 0.0018899503327676891843, -0.014952258337062199118)  , se);\
so = mad(*tmp1b,   (double2)(-0.027219029917103486322,   0.36444189483617893676)   , so);\
se = mad(*tmp1b++, (double2)( 0.049137179673730286787,  -0.051945838107881800736)  , se);\
so = mad(*tmp1b,   (double2)( 0.48135965125905339159,   -0.14329423835127266284)   , so);\
se = mad(*tmp1b++, (double2)( 0.77718575169962802862,   -0.061273359067811077843)  , se);\
so = mad(*tmp1b,   (double2)( 0.031695087811525991431,  -0.0033824159510050025955) , so);\
se = mad(*tmp1b++, (double2)( 0.0076074873249766081919, -0.00054213233180001068935), se);\
so = mad(*tmp2b,   (double2)( 0.00054213233180001068935,-0.0076074873249766081919) , so);\
se = mad(*tmp2b++, (double2)(-0.0033824159510050025955,  0.031695087811525991431)  , se);\
so = mad(*tmp2b,   (double2)( 0.061273359067811077843,  -0.77718575169962802862)   , so);\
se = mad(*tmp2b++, (double2)(-0.14329423835127266284,    0.48135965125905339159)   , se);\
so = mad(*tmp2b,   (double2)( 0.051945838107881800736,  -0.049137179673730286787)  , so);\
se = mad(*tmp2b++, (double2)( 0.36444189483617893676,   -0.027219029917103486322)  , se);\
so = mad(*tmp2b,   (double2)( 0.014952258337062199118,  -0.0018899503327676891843) , so);\
se = mad(*tmp2b++, (double2)( 0.0038087520138944894631, -0.00030292051472413308126), se);\n\
#define filter_grow(tmp1,tmp2) \
double se = 0.0;\
double so = 0.0;\
__local double *tmp1b = tmp1;\
__local double *tmp2b = tmp2+7;\
so = mad(*tmp2b,   -0.0018899503327676891843,  so);\
se = mad(*tmp2b--, -0.00030292051472413308126, se);\
so = mad(*tmp1b,   -0.00030292051472413308126, so);\
se = mad(*tmp1b++,  0.0018899503327676891843,  se);\
so = mad(*tmp2b,    0.014952258337062199118,   so);\
se = mad(*tmp2b--,  0.0038087520138944894631,  se);\
so = mad(*tmp1b,    0.0038087520138944894631,  so);\
se = mad(*tmp1b++, -0.014952258337062199118,   se);\
so = mad(*tmp2b,   -0.049137179673730286787,   so);\
se = mad(*tmp2b--, -0.027219029917103486322,   se);\
so = mad(*tmp1b,   -0.027219029917103486322,   so);\
se = mad(*tmp1b++,  0.049137179673730286787,   se);\
so = mad(*tmp2b,    0.051945838107881800736,   so);\
se = mad(*tmp2b--,  0.36444189483617893676,    se);\
so = mad(*tmp1b,    0.36444189483617893676,    so);\
se = mad(*tmp1b++, -0.051945838107881800736,   se);\
so = mad(*tmp2b,   -0.77718575169962802862,    so);\
se = mad(*tmp2b--,  0.48135965125905339159,    se);\
so = mad(*tmp1b,    0.48135965125905339159,    so);\
se = mad(*tmp1b++,  0.77718575169962802862,    se);\
so = mad(*tmp2b,    0.061273359067811077843,   so);\
se = mad(*tmp2b--, -0.14329423835127266284,    se);\
so = mad(*tmp1b,   -0.14329423835127266284,    so);\
se = mad(*tmp1b++, -0.061273359067811077843,   se);\
so = mad(*tmp2b,   -0.0076074873249766081919,  so);\
se = mad(*tmp2b--,  0.031695087811525991431,   se);\
so = mad(*tmp1b,    0.031695087811525991431,   so);\
se = mad(*tmp1b++,  0.0076074873249766081919,  se);\
so = mad(*tmp2b,    0.00054213233180001068935, so);\
se = mad(*tmp2b--, -0.0033824159510050025955,  se);\
so = mad(*tmp1b,   -0.0033824159510050025955,  so);\
se = mad(*tmp1b++, -0.00054213233180001068935, se);\n\
";
}
static void generate_syn1dKernel(std::stringstream &program, struct bigdft_device_infos * infos){
  program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<", 1))) void syn1dKernel_d(uint n, uint ndat, __global const double * restrict psi, __global double * restrict out){\n\
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
barrier(CLK_LOCAL_MEM_FENCE);\n";
  if(strncmp(infos->NAME,"Cayman",strlen("Cayman"))==0){
    program<<"filter_vector2(&tmp[-4],&tmp[FILTER_WIDTH+SIZE_I - 4]);\n\
out[jg*(2*n)+ig*2]=se.x+se.y;\n\
out[jg*(2*n)+ig*2+1]=so.x+so.y;\n\
};\n";
  } else {
    program<<"filter(&tmp[-4],&tmp[FILTER_WIDTH+SIZE_I - 4]);\n\
out[jg*(2*n)+ig*2]=se;\n\
out[jg*(2*n)+ig*2+1]=so;\n\
};\n";
  }
}

static void generate_syngrow1dKernel(std::stringstream &program, struct bigdft_device_infos * infos){
  program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<", 1))) void syngrow1dKernel_d(uint n, uint ndat, __global const double * restrict psi, __global double * restrict out){\n\
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
barrier(CLK_LOCAL_MEM_FENCE);\n";
  if(strncmp(infos->NAME,"Cayman",strlen("Cayman"))==0){
    program<<"filter_grow_vector2(&tmp[-3],&tmp[FILTER_WIDTH+SIZE_I - 3]);\n\
out[jg*2*n+2*ig]=so.x+so.y;\n\
out[jg*2*n+2*ig+1]=se.x+se.y;\n\
};\n";
  } else {
    program<<"filter_grow(&tmp[-3],&tmp[FILTER_WIDTH+SIZE_I - 3]);\n\
out[jg*2*n+2*ig]=so;\n\
out[jg*2*n+2*ig+1]=se;\n\
};\n";
  }
}

extern "C" char* generate_ana_program(struct bigdft_device_infos * infos){
  char * output;
  std::stringstream program;

  generate_ana_header(program);
  generate_ana_filter(program);
  generate_ana1dKernel(program,infos);
  generate_anashrink1dKernel(program,infos);
  generate_ana1d_blockKernel(program);

  output = (char *)malloc((program.str().size()+1)*sizeof(char));
  strcpy(output, program.str().c_str());
  return output;
}

extern "C" char* generate_syn_program(struct bigdft_device_infos * infos){
  char * output;
  std::stringstream program;

  generate_syn_header(program);
  generate_syn_filter(program);
  generate_syn1dKernel(program,infos);
  generate_syngrow1dKernel(program,infos);
  
  output = (char *)malloc((program.str().size()+1)*sizeof(char));
  strcpy(output, program.str().c_str());
  return output;
}

