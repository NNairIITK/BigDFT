#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include "MagicFilter_Generator.h"
#include "OpenCL_wrappers.h"

#define FILTER_WIDTH 16
static void generate_header(std::stringstream &program){
  program<<"#ifdef cl_khr_fp64\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#elif defined (cl_amd_fp64)\n\
#pragma OPENCL EXTENSION cl_amd_fp64: enable \n\
#endif\n\
#define FILTER_WIDTH "<<FILTER_WIDTH<<"\n\
#define FILT0   8.4334247333529341094733325815816e-7\n\
#define FILT1  -0.1290557201342060969516786758559028e-4\n\
#define FILT2   0.8762984476210559564689161894116397e-4\n\
#define FILT3  -0.30158038132690463167163703826169879e-3\n\
#define FILT4   0.174723713672993903449447812749852942e-2\n\
#define FILT5  -0.942047030201080385922711540948195075e-2\n\
#define FILT6   0.2373821463724942397566389712597274535e-1\n\
#define FILT7   0.612625895831207982195380597e-1\n\
#define FILT8   0.9940415697834003993178616713\n\
#define FILT9  -0.604895289196983516002834636e-1\n\
#define FILT10 -0.2103025160930381434955489412839065067e-1\n\
#define FILT11  0.1337263414854794752733423467013220997e-1\n\
#define FILT12 -0.344128144493493857280881509686821861e-2\n\
#define FILT13  0.49443227688689919192282259476750972e-3\n\
#define FILT14 -0.5185986881173432922848639136911487e-4\n\
#define FILT15  2.72734492911979659657715313017228e-6\n\
#define FILTER_WIDTH 16\n\
#define filter_nomad(tt,tmp) \
tt += *tmp++ * FILT0;\
tt += *tmp++ * FILT1;\
tt += *tmp++ * FILT2;\
tt += *tmp++ * FILT3;\
tt += *tmp++ * FILT4;\
tt += *tmp++ * FILT5;\
tt += *tmp++ * FILT6;\
tt += *tmp++ * FILT7;\
tt += *tmp++ * FILT8;\
tt += *tmp++ * FILT9;\
tt += *tmp++ * FILT10;\
tt += *tmp++ * FILT11;\
tt += *tmp++ * FILT12;\
tt += *tmp++ * FILT13;\
tt += *tmp++ * FILT14;\
tt += *tmp++ * FILT15;\n\
#define filter(tmp) \
double tt = 0.0;\
tt = mad(*tmp++, FILT0, tt);\
tt = mad(*tmp++, FILT1, tt);\
tt = mad(*tmp++, FILT2, tt);\
tt = mad(*tmp++, FILT3, tt);\
tt = mad(*tmp++, FILT4, tt);\
tt = mad(*tmp++, FILT5, tt);\
tt = mad(*tmp++, FILT6, tt);\
tt = mad(*tmp++, FILT7, tt);\
tt = mad(*tmp++, FILT8, tt);\
tt = mad(*tmp++, FILT9, tt);\
tt = mad(*tmp++, FILT10, tt);\
tt = mad(*tmp++, FILT11, tt);\
tt = mad(*tmp++, FILT12, tt);\
tt = mad(*tmp++, FILT13, tt);\
tt = mad(*tmp++, FILT14, tt);\
tt = mad(*tmp++, FILT15, tt);\n\
#define filterp(tt,tmp) \
tt = mad(*tmp++, FILT0, tt);\
tt = mad(*tmp++, FILT1, tt);\
tt = mad(*tmp++, FILT2, tt);\
tt = mad(*tmp++, FILT3, tt);\
tt = mad(*tmp++, FILT4, tt);\
tt = mad(*tmp++, FILT5, tt);\
tt = mad(*tmp++, FILT6, tt);\
tt = mad(*tmp++, FILT7, tt);\
tt = mad(*tmp++, FILT8, tt);\
tt = mad(*tmp++, FILT9, tt);\
tt = mad(*tmp++, FILT10, tt);\
tt = mad(*tmp++, FILT11, tt);\
tt = mad(*tmp++, FILT12, tt);\
tt = mad(*tmp++, FILT13, tt);\
tt = mad(*tmp++, FILT14, tt);\
tt = mad(*tmp++, FILT15, tt);\n\
#define filter_vector2(tmp) \
double2 tt = (double2)(0.0, 0.0);\
__local double2 *tmp2= (__local double2 *)tmp;\
tt = mad(*tmp2++, (double2)(FILT0,FILT1), tt);\
tt = mad(*tmp2++, (double2)(FILT2,FILT3), tt);\
tt = mad(*tmp2++, (double2)(FILT4,FILT5), tt);\
tt = mad(*tmp2++, (double2)(FILT6,FILT7), tt);\
tt = mad(*tmp2++, (double2)(FILT8,FILT9), tt);\
tt = mad(*tmp2++, (double2)(FILT10,FILT11), tt);\
tt = mad(*tmp2++, (double2)(FILT12,FILT13), tt);\
tt = mad(*tmp2++, (double2)(FILT14,FILT15), tt);\n\
#define filter_vector4(tt,tmp) \
tt = mad(*tmp++, (double4)(FILT0,FILT1,FILT2,FILT3), tt);\
tt = mad(*tmp++, (double4)(FILT4,FILT5,FILT6,FILT7), tt);\
tt = mad(*tmp++, (double4)(FILT8,FILT9,FILT10,FILT11), tt);\
tt = mad(*tmp++, (double4)(FILT12,FILT13,FILT14,FILT15), tt);\n\
#define filter_reverse_nomad(tt,tmp) \
tt += *tmp++ *  FILT15;\
tt += *tmp++ *  FILT14;\
tt += *tmp++ *  FILT13;\
tt += *tmp++ *  FILT12;\
tt += *tmp++ *  FILT11;\
tt += *tmp++ *  FILT10;\
tt += *tmp++ *  FILT9;\
tt += *tmp++ *  FILT8;\
tt += *tmp++ *  FILT7;\
tt += *tmp++ *  FILT6;\
tt += *tmp++ *  FILT5;\
tt += *tmp++ *  FILT4;\
tt += *tmp++ *  FILT3;\
tt += *tmp++ *  FILT2;\
tt += *tmp++ *  FILT1;\
tt += *tmp++ *  FILT0;\n\
#define filter_reverse(tmp) \
double tt = 0.0;\
tt = mad(*tmp++, FILT15, tt);\
tt = mad(*tmp++, FILT14, tt);\
tt = mad(*tmp++, FILT13, tt);\
tt = mad(*tmp++, FILT12, tt);\
tt = mad(*tmp++, FILT11, tt);\
tt = mad(*tmp++, FILT10, tt);\
tt = mad(*tmp++, FILT9, tt);\
tt = mad(*tmp++, FILT8, tt);\
tt = mad(*tmp++, FILT7, tt);\
tt = mad(*tmp++, FILT6, tt);\
tt = mad(*tmp++, FILT5, tt);\
tt = mad(*tmp++, FILT4, tt);\
tt = mad(*tmp++, FILT3, tt);\
tt = mad(*tmp++, FILT2, tt);\
tt = mad(*tmp++, FILT1, tt);\
tt = mad(*tmp++, FILT0, tt);\n\
#define filter_reverse_vector2(tmp) \
double2 tt = (double2)(0.0, 0.0);\
__local double2 *tmp2= (__local double2 *)tmp;\
tt = mad(*tmp2++, (double2)(FILT15,FILT14), tt);\
tt = mad(*tmp2++, (double2)(FILT13,FILT12), tt);\
tt = mad(*tmp2++, (double2)(FILT11,FILT10), tt);\
tt = mad(*tmp2++, (double2)(FILT9,FILT8), tt);\
tt = mad(*tmp2++, (double2)(FILT7,FILT6), tt);\
tt = mad(*tmp2++, (double2)(FILT5,FILT4), tt);\
tt = mad(*tmp2++, (double2)(FILT3,FILT2), tt);\
tt = mad(*tmp2++, (double2)(FILT1,FILT0), tt);\n\
";
}

static void generate_magicfilter1dKernel(std::stringstream &program, struct bigdft_device_infos * infos){
  program<<"//n is supposed to be greater or equal than get_local_size(0)\n\
//this filter is for periodic boundary conditions\n\
__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<", 1))) void magicfilter1dKernel_d(uint n, uint ndat, __global const double * restrict psi, __global double * restrict out){\n\
__local double tmp1[FILTER_WIDTH*(2*FILTER_WIDTH+1)];\n\
__local double *tmp = &tmp1[0];\n\
//get our position in the local work group\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
//get our position in the result matrix\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1);\n\
//get our group number\n\
ptrdiff_t igt = get_group_id(0);\n\
ptrdiff_t jgt = get_group_id(1);\n\
//if data are ill dimentioned border blocks recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
//transpose indexes in the work group in order to read transposed data\n\
igt = ig - i2 + j2 - FILTER_WIDTH/2;\n\
jgt = jg - j2 + i2;\n\
//if we are on the outside, select a border element to load, wrapping around\n\
//we will be loading 2 elements each\n\
if ( igt < 0 ) \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = psi[jgt + ( n + igt ) * ndat];\n\
else \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = psi[jgt + igt * ndat];\n\
igt += FILTER_WIDTH;\n\
if ( igt >= n ) \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = psi[jgt + ( igt - n ) * ndat];\n\
else\n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = psi[jgt +  igt * ndat];\n\
//rest position in the buffer to first element involved in the convolution\n\
tmp += j2*(2*FILTER_WIDTH+1) + i2;\n\
//wait for buffer to be full\n\
barrier(CLK_LOCAL_MEM_FENCE);\n";
  if(strncmp(infos->NAME,"Cayman",strlen("Cayman"))==0){
    program<<"//apply filter\n\
filter_vector2(tmp);\n\
//store the result\n\
out[(jg*n+ig)]=tt.x+tt.y;\n\
};\n";
  } else {
    program<<"//apply filter\n\
filter(tmp);\n\
//store the result\n\
out[(jg*n+ig)]=tt;\n\
};\n";
  }
}

static void generate_magicfilter1d_tKernel(std::stringstream &program, struct bigdft_device_infos * infos){
program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<", 1))) void magicfilter1d_tKernel_d(uint n, uint ndat, __global const double * restrict psi, __global double * restrict out){\n\
__local double tmp1[FILTER_WIDTH*(2*FILTER_WIDTH+1)];\n\
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
igt = ig - i2 + j2 - FILTER_WIDTH/2;\n\
jgt = jg - j2 + i2;\n\
//If I'm on the outside, select a border element to load\n\
if ( igt < 0 ) \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = psi[jgt + ( n + igt ) * ndat];\n\
else \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = psi[jgt + igt * ndat];\n\
igt += FILTER_WIDTH;\n\
if ( igt >= n ) \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = psi[jgt + ( igt - n ) * ndat];\n\
else\n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = psi[jgt +  igt * ndat];\n\
tmp += j2*(2*FILTER_WIDTH+1) + i2 + 1;\n\
barrier(CLK_LOCAL_MEM_FENCE);\n";
  if(strncmp(infos->NAME,"Cayman",strlen("Cayman"))==0){
    program<<"filter_reverse_vector2(tmp);\n\
out[(jg*n+ig)]=tt.x+tt.y;\n\
};\n";
  } else {
    program<<"filter_reverse(tmp);\n\
out[(jg*n+ig)]=tt;\n\
};\n";
  }
}

static void generate_magicfiltergrow1dKernel(std::stringstream &program, struct bigdft_device_infos * infos){
  program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<", 1))) void magicfiltergrow1dKernel_d(uint n, uint ndat, __global const double * restrict psi, __global double * restrict out){\n\
__local double tmp1[FILTER_WIDTH*(2*FILTER_WIDTH+1)];\n\
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
igt = ig - i2 + j2 - FILTER_WIDTH/2 - 7;\n\
jgt = jg - j2 + i2;\n\
//If I'm on the outside, select a border element to load\n\
if ( igt < 0 ) \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = 0.0;\n\
else \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = psi[jgt + igt * ndat];\n\
igt += FILTER_WIDTH;\n\
if ( igt >= n - 15 ) \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = 0.0;\n\
else\n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = psi[jgt +  igt * ndat];\n\
tmp += j2*(2*FILTER_WIDTH+1) + i2;\n\
barrier(CLK_LOCAL_MEM_FENCE);\n";
  if(strncmp(infos->NAME,"Cayman",strlen("Cayman"))==0){
  program<<"filter_vector2(tmp);\n\
out[(jg*n+ig)]=tt.x+tt.y;\n\
};\n";
  } else {
  program<<"filter(tmp);\n\
out[(jg*n+ig)]=tt;\n\
};\n";
  }
}

static void generate_magicfiltergrow1d_denKernel(std::stringstream &program, struct bigdft_device_infos * infos){
  program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<", 1))) void magicfiltergrow1d_denKernel_d(uint n, uint ndat, __global const double * restrict psi, __global double * restrict out){\n\
__local double tmp1[FILTER_WIDTH*(2*FILTER_WIDTH+1)];\n\
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
igt = ig - i2 + j2 - FILTER_WIDTH/2 - 7;\n\
jgt = jg - j2 + i2;\n\
//If I'm on the outside, select a border element to load\n\
if ( igt < 0 ) \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = 0.0;\n\
else \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = psi[jgt + igt * ndat];\n\
igt += FILTER_WIDTH;\n\
if ( igt >= n - 15 ) \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = 0.0;\n\
else\n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = psi[jgt +  igt * ndat];\n\
tmp += j2*(2*FILTER_WIDTH+1) + i2;\n\
barrier(CLK_LOCAL_MEM_FENCE);\n";
  if(strncmp(infos->NAME,"Cayman",strlen("Cayman"))==0){
    program<<"filter_vector2(tmp);\n\
tt.x += tt.y;\n\
out[(jg*n+ig)]=tt.x*tt.x;\n\
};\n";
  } else {
    program<<"filter(tmp);\n\
out[(jg*n+ig)]=tt*tt;\n\
};\n";
  }
}

static void generate_magicfiltergrow1d_potKernel(std::stringstream &program, struct bigdft_device_infos * infos){
  program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<", 1))) void magicfiltergrow1d_potKernel_d(uint n, uint ndat, __global const double *psi, __global const double * restrict pot, __global double * restrict out){\n\
__local double tmp1[FILTER_WIDTH*(2*FILTER_WIDTH+1)];\n\
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
igt = ig - i2 + j2 - FILTER_WIDTH/2 - 7;\n\
jgt = jg - j2 + i2;\n\
//If I'm on the outside, select a border element to load\n\
if ( igt < 0 ) \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = 0.0;\n\
else \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = psi[jgt + igt * ndat];\n\
igt += FILTER_WIDTH;\n\
if ( igt >= n - 15 ) \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = 0.0;\n\
else\n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = psi[jgt +  igt * ndat];\n\
tmp += j2*(2*FILTER_WIDTH+1) + i2;\n\
barrier(CLK_LOCAL_MEM_FENCE);\n";
  if(strncmp(infos->NAME,"Cayman",strlen("Cayman"))==0){
    program<<"filter_vector2(tmp);\n\
out[(jg*n+ig)]=(tt.x+tt.y)*pot[(jg*n+ig)];\n\
};\n";
  } else {
    program<<"filter(tmp);\n\
out[(jg*n+ig)]=tt*pot[(jg*n+ig)];\n\
};\n";
  }
}

static void generate_magicfiltershrink1dKernel(std::stringstream &program, struct bigdft_device_infos * infos){
  program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<", 1))) void magicfiltershrink1dKernel_d(uint n, uint ndat, __global const double * restrict psi, __global double * restrict out){\n\
__local double tmp1[FILTER_WIDTH*(2*FILTER_WIDTH+1)];\n\
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
igt = ig - i2 + j2;\n\
jgt = jg - j2 + i2;\n\
igt -= FILTER_WIDTH/2;\n\
psi = psi + 8 * ndat;\n\
//If I'm on the outside, select a border element to load\n\
tmp[i2 * (2 * FILTER_WIDTH + 1) + j2]=psi[jgt + igt * ndat];\n\
igt += FILTER_WIDTH - 1;\n\
tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH - 1]=psi[jgt + igt * ndat];\n\
tmp += j2*(2*FILTER_WIDTH+1) + i2;\n\
barrier(CLK_LOCAL_MEM_FENCE);\n";
  if(strncmp(infos->NAME,"Cayman",strlen("Cayman"))==0){
    program<<"filter_reverse_vector2(tmp);\n\
out[(jg*n+ig)]=tt.x+tt.y;\n\
};\n";
  } else {
    program<<"filter_reverse(tmp);\n\
out[(jg*n+ig)]=tt;\n\
};\n";
  }
}

static void generate_magicfilter1d_potKernel(std::stringstream &program, struct bigdft_device_infos * infos){
  program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<", 1))) void magicfilter1d_potKernel_d(uint n, uint ndat, __global const double *psi, __global double * restrict pot, __global double * restrict out){\n\
__local double tmp1[FILTER_WIDTH*(2*FILTER_WIDTH+1)];\n\
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
igt = ig - i2 + j2;\n\
jgt = jg - j2 + i2;\n\
igt -= FILTER_WIDTH/2;\n\
//If I'm on the outside, select a border element to load\n\
if ( igt < 0 ) \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = psi[jgt + ( n + igt ) * ndat];\n\
else \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = psi[jgt + igt * ndat];\n\
igt += FILTER_WIDTH;\n\
if ( igt >= n ) \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = psi[jgt + ( igt - n ) * ndat];\n\
else\n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = psi[jgt +  igt * ndat];\n\
tmp += j2*(2*FILTER_WIDTH+1) + i2;\n\
barrier(CLK_LOCAL_MEM_FENCE);\n";
  if(strncmp(infos->NAME,"Cayman",strlen("Cayman"))==0){
    program<<"filter_vector2(tmp);\n\
out[(jg*n+ig)]=(tt.x+tt.y)*pot[jg*n+ig];\n\
};\n";
  } else {
    program<<"filter(tmp);\n\
out[(jg*n+ig)]=tt*pot[jg*n+ig];\n\
};\n";
  }
}

static void generate_magicfilter1d_blockKernel(std::stringstream &program, struct bigdft_device_infos * infos){
program<<"#define ELEM_PER_THREAD 2\n\
__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<"/ELEM_PER_THREAD, 1))) void magicfilter1d_blockKernel_d(uint n, uint ndat, __global const double * restrict psi, __global double * restrict out){\n\
__local double tmp1[FILTER_WIDTH*(2*FILTER_WIDTH+1)];\n\
__local double *tmp = &tmp1[0];\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1)*ELEM_PER_THREAD;\n\
size_t i2 = get_local_id(0);\n\
size_t j2 = get_local_id(1)*ELEM_PER_THREAD;\n\
ptrdiff_t igt = get_group_id(0);\n\
ptrdiff_t jgt = get_group_id(1);\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1)*ELEM_PER_THREAD - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
igt = ig - i2 + j2 - FILTER_WIDTH/2;\n\
jgt = jg - j2 + i2;\n\
//If I'm on the outside, select a border element to load\n\
__local double * tmp_1 = tmp + i2 * (2 * FILTER_WIDTH + 1) + j2;\n\
psi += jgt;\n\
if ( igt < 0 )\n\
  *tmp_1++ = psi[(n + igt++) * ndat];\n\
else \n\
  *tmp_1++ = psi[igt++ * ndat];\n\
if ( igt < 0 )\n\
  *tmp_1++ = psi[(n + igt++) * ndat];\n\
else \n\
  *tmp_1++ = psi[igt++ * ndat];\n\
igt += FILTER_WIDTH-ELEM_PER_THREAD;\n\
tmp_1 += FILTER_WIDTH-ELEM_PER_THREAD;\n\
if ( igt >= n )\n\
  *tmp_1++ = psi[(igt++ - n) * ndat];\n\
else\n\
  *tmp_1++ = psi[igt++ * ndat];\n\
if ( igt >= n )\n\
  *tmp_1++ = psi[(igt++ - n) * ndat];\n\
else\n\
  *tmp_1++ = psi[igt++ * ndat];\n\
\
tmp_1 = tmp + j2*(2*FILTER_WIDTH+1) + i2;\n\
__local double * tmp_2 = tmp_1 + (2*FILTER_WIDTH+1);\n\
out += jg*n + ig;\n\
double tt_1 = 0.0;\n\
double tt_2 = 0.0;\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
filterp(tt_1,tmp_1);\n\
filterp(tt_2,tmp_2);\n\
*out = tt_1;\n\
out += n;\n\
*out = tt_2;\n\
};\n";
}

static void generate_magicfilter1d_straightKernel(std::stringstream &program, struct bigdft_device_infos * infos){
program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<", 1))) void magicfilter1d_straightKernel_d(uint n, uint ndat, __global const double * restrict psi, __global double * restrict out){\n\
__local double tmp1[FILTER_WIDTH*(2*FILTER_WIDTH+1)];\n\
__local double *tmp = &tmp1[0];\n\
ptrdiff_t ig = get_global_id(0);\n\
ptrdiff_t jg = get_global_id(1);\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1);\n\
ptrdiff_t igt = get_group_id(0);\n\
ptrdiff_t jgt = get_group_id(1);\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
igt = ig - i2 + j2;\n\
ig -= FILTER_WIDTH/2;\n\
jgt = jg - j2 + i2;\n\
//If I'm on the outside, select a border element to load\n\
if ( ig < 0 ) \n\
  tmp[j2 * (2 * FILTER_WIDTH + 1) + i2] = psi[jg*n + n + ig];\n\
else \n\
  tmp[j2 * (2 * FILTER_WIDTH + 1) + i2] = psi[jg*n + ig];\n\
ig += FILTER_WIDTH;\n\
if ( ig >= n ) \n\
  tmp[j2 * (2 * FILTER_WIDTH + 1) + i2 + FILTER_WIDTH] = psi[jg*n + ig - n];\n\
else\n\
  tmp[j2 * (2 * FILTER_WIDTH + 1) + i2 + FILTER_WIDTH] = psi[jg*n + ig];\n\
tmp += i2*(2*FILTER_WIDTH+1) + j2;\n\
barrier(CLK_LOCAL_MEM_FENCE);\n";
  if(strncmp(infos->NAME,"Cayman",strlen("Cayman"))==0){
    program<<"filter_vector2(tmp);\n\
out[(igt*ndat+jgt)]=tt.x+tt.y;\n\
};\n";
  } else {
    program<<"filter(tmp);\n\
out[(igt*ndat+jgt)]=tt;\n\
};\n";
  }
}

static void generate_magicfilter1d_denKernel(std::stringstream &program, struct bigdft_device_infos * infos){
program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<", 1))) void magicfilter1d_denKernel_d(uint n, uint ndat, __global const double * restrict psi, __global double * restrict out){\n\
__local double tmp1[FILTER_WIDTH*(2*FILTER_WIDTH+1)];\n\
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
igt = ig - i2 + j2 - FILTER_WIDTH/2;\n\
jgt = jg - j2 + i2;\n\
//If I'm on the outside, select a border element to load\n\
if ( igt < 0 ) \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = psi[jgt + ( n + igt ) * ndat];\n\
else \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = psi[jgt + igt * ndat];\n\
igt += FILTER_WIDTH;\n\
if ( igt >= n ) \n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = psi[jgt + ( igt - n ) * ndat];\n\
else\n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = psi[jgt +  igt * ndat];\n\
tmp += j2*(2*FILTER_WIDTH+1) + i2;\n\
barrier(CLK_LOCAL_MEM_FENCE);\n";
  if(strncmp(infos->NAME,"Cayman",strlen("Cayman"))==0){
    program<<"filter_vector2(tmp);\n\
tt.x+=tt.y;\n\
out[(jg*n+ig)]=tt.x*tt.x;\n\
};\n";
  } else {
    program<<"filter(tmp);\n\
out[(jg*n+ig)]=tt*tt;\n\
};\n";
  }
}


extern "C" char* generate_magicfilter_program(struct bigdft_device_infos * infos){
  char * output;
  std::stringstream program;

  generate_header(program);
  generate_magicfilter1dKernel(program,infos);
  generate_magicfiltergrow1dKernel(program,infos);
  generate_magicfiltergrow1d_denKernel(program,infos);
  generate_magicfiltergrow1d_potKernel(program,infos);
  generate_magicfiltershrink1dKernel(program,infos);
  generate_magicfilter1d_potKernel(program,infos);
  generate_magicfilter1d_blockKernel(program,infos);
  generate_magicfilter1d_straightKernel(program,infos);
  generate_magicfilter1d_denKernel(program,infos);
  generate_magicfilter1d_tKernel(program,infos);

  output = (char *)malloc((program.str().size()+1)*sizeof(char));
  strcpy(output, program.str().c_str());
  return output;
}


