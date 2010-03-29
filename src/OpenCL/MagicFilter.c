#include "MagicFilter.h"
#include "OpenCL_wrappers.h"

char * magicfilter_program="\
#define FILTER_WIDTH 16\n\
//n is supposed to be greater or equal than get_local_size(0)\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
__kernel void magicfiltergrow1dKernel_d(uint n, uint ndat, __global const double *psi, __global double *out, __local double tmp[]){\n\
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
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
double tt = 0.0;\n\
__local double * tmp_o = tmp + j2*(2*FILTER_WIDTH+1) + FILTER_WIDTH/2+i2;\n\
tt += tmp_o[-8] *  8.4334247333529341094733325815816e-7\n\
    + tmp_o[+7] *  2.72734492911979659657715313017228e-6\n\
    + tmp_o[-7] * -0.1290557201342060969516786758559028e-4\n\
    + tmp_o[+6] * -0.5185986881173432922848639136911487e-4\n\
    + tmp_o[-6] *  0.8762984476210559564689161894116397e-4\n\
    + tmp_o[-5] * -0.30158038132690463167163703826169879e-3\n\
    + tmp_o[+5] *  0.49443227688689919192282259476750972e-3\n\
    + tmp_o[-4] *  0.174723713672993903449447812749852942e-2;\n\
tt += tmp_o[+4] * -0.344128144493493857280881509686821861e-2\n\
    + tmp_o[-3] * -0.942047030201080385922711540948195075e-2\n\
    + tmp_o[+3] *  0.1337263414854794752733423467013220997e-1\n\
    + tmp_o[+2] * -0.2103025160930381434955489412839065067e-1\n\
    + tmp_o[-2] *  0.2373821463724942397566389712597274535e-1\n\
    + tmp_o[+1] * -0.604895289196983516002834636e-1\n\
    + tmp_o[-1] *  0.612625895831207982195380597e-1\n\
    + tmp_o[+0] *  0.9940415697834003993178616713;\n\
out[(jg*n+ig)]=tt;\n\
};\n\
__kernel void magicfiltergrow1d_potKernel_d(uint n, uint ndat, __global const double *psi, __global const double *pot, __global double *out, __local double tmp[]){\n\
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
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
double tt = 0.0;\n\
__local double * tmp_o = tmp + j2*(2*FILTER_WIDTH+1) + FILTER_WIDTH/2+i2;\n\
tt += tmp_o[-8] *  8.4334247333529341094733325815816e-7\n\
    + tmp_o[+7] *  2.72734492911979659657715313017228e-6\n\
    + tmp_o[-7] * -0.1290557201342060969516786758559028e-4\n\
    + tmp_o[+6] * -0.5185986881173432922848639136911487e-4\n\
    + tmp_o[-6] *  0.8762984476210559564689161894116397e-4\n\
    + tmp_o[-5] * -0.30158038132690463167163703826169879e-3\n\
    + tmp_o[+5] *  0.49443227688689919192282259476750972e-3\n\
    + tmp_o[-4] *  0.174723713672993903449447812749852942e-2;\n\
tt += tmp_o[+4] * -0.344128144493493857280881509686821861e-2\n\
    + tmp_o[-3] * -0.942047030201080385922711540948195075e-2\n\
    + tmp_o[+3] *  0.1337263414854794752733423467013220997e-1\n\
    + tmp_o[+2] * -0.2103025160930381434955489412839065067e-1\n\
    + tmp_o[-2] *  0.2373821463724942397566389712597274535e-1\n\
    + tmp_o[+1] * -0.604895289196983516002834636e-1\n\
    + tmp_o[-1] *  0.612625895831207982195380597e-1\n\
    + tmp_o[+0] *  0.9940415697834003993178616713;\n\
out[(jg*n+ig)]=tt*pot[(jg*n+ig)];\n\
};\n\
__kernel void magicfiltershrink1dKernel_d(uint n, uint ndat, __global const double *psi, __global double *out, __local double tmp[]){\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1);\n\
size_t igt = get_group_id(0);\n\
size_t jgt = get_group_id(1);\n\
size_t jb;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
igt = ig - i2 + j2;\n\
jgt = jg - j2 + i2;\n\
psi = psi + 8 * ndat;\n\
//If I'm on the outside, select a border element to load\n\
if(j2 < FILTER_WIDTH/2) {\n\
    jb = igt - FILTER_WIDTH/2; \n\
    tmp[i2 * (2 * FILTER_WIDTH + 1) + j2]=psi[jgt + jb * ndat];\n\
  }\n\
if (j2 >= FILTER_WIDTH/2+1) {\n\
    jb = igt + FILTER_WIDTH/2-1;\n\
    tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH - 1]=psi[jgt + jb * ndat];\n\
  }\n\
//check boundaries\
//Load the element I am to calculate\n\
tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH/2]=psi[jgt + igt * ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
double tt = 0.0;\n\
__local double * tmp_o = tmp + j2*(2*FILTER_WIDTH+1) + FILTER_WIDTH/2+i2;\n\
tt += tmp_o[+7] *  8.4334247333529341094733325815816e-7\n\
    + tmp_o[-8] *  2.72734492911979659657715313017228e-6\n\
    + tmp_o[+6] * -0.1290557201342060969516786758559028e-4\n\
    + tmp_o[-7] * -0.5185986881173432922848639136911487e-4\n\
    + tmp_o[+5] *  0.8762984476210559564689161894116397e-4\n\
    + tmp_o[+4] * -0.30158038132690463167163703826169879e-3\n\
    + tmp_o[-6] *  0.49443227688689919192282259476750972e-3\n\
    + tmp_o[+3] *  0.174723713672993903449447812749852942e-2;\n\
tt += tmp_o[-5] * -0.344128144493493857280881509686821861e-2\n\
    + tmp_o[+2] * -0.942047030201080385922711540948195075e-2\n\
    + tmp_o[-4] *  0.1337263414854794752733423467013220997e-1\n\
    + tmp_o[-3] * -0.2103025160930381434955489412839065067e-1\n\
    + tmp_o[+1] *  0.2373821463724942397566389712597274535e-1\n\
    + tmp_o[-2] * -0.604895289196983516002834636e-1\n\
    + tmp_o[+0] *  0.612625895831207982195380597e-1\n\
    + tmp_o[-1] *  0.9940415697834003993178616713;\n\
out[(jg*n+ig)]=tt;\n\
};\n\
__kernel void magicfilter1d_potKernel_d(uint n, uint ndat, __global const double *psi, __global double *pot, __global double *out, __local double tmp[]){\n\
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
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
double tt = 0.0;\n\
__local double * tmp_o = tmp + j2*(2*FILTER_WIDTH+1) + FILTER_WIDTH/2+i2;\n\
tt += tmp_o[-8] *  8.4334247333529341094733325815816e-7\n\
    + tmp_o[+7] *  2.72734492911979659657715313017228e-6\n\
    + tmp_o[-7] * -0.1290557201342060969516786758559028e-4\n\
    + tmp_o[+6] * -0.5185986881173432922848639136911487e-4\n\
    + tmp_o[-6] *  0.8762984476210559564689161894116397e-4\n\
    + tmp_o[-5] * -0.30158038132690463167163703826169879e-3\n\
    + tmp_o[+5] *  0.49443227688689919192282259476750972e-3\n\
    + tmp_o[-4] *  0.174723713672993903449447812749852942e-2;\n\
tt += tmp_o[+4] * -0.344128144493493857280881509686821861e-2\n\
    + tmp_o[-3] * -0.942047030201080385922711540948195075e-2\n\
    + tmp_o[+3] *  0.1337263414854794752733423467013220997e-1\n\
    + tmp_o[+2] * -0.2103025160930381434955489412839065067e-1\n\
    + tmp_o[-2] *  0.2373821463724942397566389712597274535e-1\n\
    + tmp_o[+1] * -0.604895289196983516002834636e-1\n\
    + tmp_o[-1] *  0.612625895831207982195380597e-1\n\
    + tmp_o[+0] *  0.9940415697834003993178616713;\n\
out[(jg*n+ig)]=tt*pot[jg*n+ig];\n\
};\n\
__kernel void magicfilter1dKernel_d(uint n, uint ndat, __global const double *psi, __global double *out, __local double tmp[]){\n\
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
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
double tt = 0.0;\n\
__local double * tmp_o = tmp + j2*(2*FILTER_WIDTH+1) + FILTER_WIDTH/2+i2;\n\
tt += tmp_o[-8] *  8.4334247333529341094733325815816e-7\n\
    + tmp_o[+7] *  2.72734492911979659657715313017228e-6\n\
    + tmp_o[-7] * -0.1290557201342060969516786758559028e-4\n\
    + tmp_o[+6] * -0.5185986881173432922848639136911487e-4\n\
    + tmp_o[-6] *  0.8762984476210559564689161894116397e-4\n\
    + tmp_o[-5] * -0.30158038132690463167163703826169879e-3\n\
    + tmp_o[+5] *  0.49443227688689919192282259476750972e-3\n\
    + tmp_o[-4] *  0.174723713672993903449447812749852942e-2;\n\
tt += tmp_o[+4] * -0.344128144493493857280881509686821861e-2\n\
    + tmp_o[-3] * -0.942047030201080385922711540948195075e-2\n\
    + tmp_o[+3] *  0.1337263414854794752733423467013220997e-1\n\
    + tmp_o[+2] * -0.2103025160930381434955489412839065067e-1\n\
    + tmp_o[-2] *  0.2373821463724942397566389712597274535e-1\n\
    + tmp_o[+1] * -0.604895289196983516002834636e-1\n\
    + tmp_o[-1] *  0.612625895831207982195380597e-1\n\
    + tmp_o[+0] *  0.9940415697834003993178616713;\n\
out[(jg*n+ig)]=tt;\n\
};\n\
__kernel void magicfilter1d_denKernel_d(uint n, uint ndat, __global const double *psi, __global double *out, __local double tmp[]){\n\
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
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
double tt = 0.0;\n\
__local double * tmp_o = tmp + j2*(2*FILTER_WIDTH+1) + FILTER_WIDTH/2+i2;\n\
tt += tmp_o[-8] *  8.4334247333529341094733325815816e-7\n\
    + tmp_o[+7] *  2.72734492911979659657715313017228e-6\n\
    + tmp_o[-7] * -0.1290557201342060969516786758559028e-4\n\
    + tmp_o[+6] * -0.5185986881173432922848639136911487e-4\n\
    + tmp_o[-6] *  0.8762984476210559564689161894116397e-4\n\
    + tmp_o[-5] * -0.30158038132690463167163703826169879e-3\n\
    + tmp_o[+5] *  0.49443227688689919192282259476750972e-3\n\
    + tmp_o[-4] *  0.174723713672993903449447812749852942e-2;\n\
tt += tmp_o[+4] * -0.344128144493493857280881509686821861e-2\n\
    + tmp_o[-3] * -0.942047030201080385922711540948195075e-2\n\
    + tmp_o[+3] *  0.1337263414854794752733423467013220997e-1\n\
    + tmp_o[+2] * -0.2103025160930381434955489412839065067e-1\n\
    + tmp_o[-2] *  0.2373821463724942397566389712597274535e-1\n\
    + tmp_o[+1] * -0.604895289196983516002834636e-1\n\
    + tmp_o[-1] *  0.612625895831207982195380597e-1\n\
    + tmp_o[+0] *  0.9940415697834003993178616713;\n\
out[(jg*n+ig)]=tt*tt;\n\
};\n\
__kernel void magicfilter1d_tKernel_d(uint n, uint ndat, __global const double *psi, __global double *out, __local double tmp[]){\n\
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
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
double tt = 0.0;\n\
__local double * tmp_o = tmp + j2*(2*FILTER_WIDTH+1) + FILTER_WIDTH/2+i2;\n\
tt += tmp_o[+8] *  8.4334247333529341094733325815816e-7\n\
    + tmp_o[-7] *  2.72734492911979659657715313017228e-6\n\
    + tmp_o[+7] * -0.1290557201342060969516786758559028e-4\n\
    + tmp_o[-6] * -0.5185986881173432922848639136911487e-4\n\
    + tmp_o[+6] *  0.8762984476210559564689161894116397e-4\n\
    + tmp_o[+5] * -0.30158038132690463167163703826169879e-3\n\
    + tmp_o[-5] *  0.49443227688689919192282259476750972e-3\n\
    + tmp_o[+4] *  0.174723713672993903449447812749852942e-2;\n\
tt += tmp_o[-4] * -0.344128144493493857280881509686821861e-2\n\
    + tmp_o[+3] * -0.942047030201080385922711540948195075e-2\n\
    + tmp_o[-3] *  0.1337263414854794752733423467013220997e-1\n\
    + tmp_o[-2] * -0.2103025160930381434955489412839065067e-1\n\
    + tmp_o[+2] *  0.2373821463724942397566389712597274535e-1\n\
    + tmp_o[-1] * -0.604895289196983516002834636e-1\n\
    + tmp_o[+1] *  0.612625895831207982195380597e-1\n\
    + tmp_o[+0] *  0.9940415697834003993178616713;\n\
out[(jg*n+ig)]=tt;\n\
};\n\
\n\
";

inline void magicfilter_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    int FILTER_WIDTH=16;
    assert(*n>=FILTER_WIDTH);
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;
    cl_uint i = 0;
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue magic filter kernel!");
}

inline void magicfilter_pot_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *pot, cl_mem *out) {
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
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue magic filter pot kernel!");
}

cl_kernel magicfilter1d_kernel_d;
cl_kernel magicfilter1d_den_kernel_d;
cl_kernel magicfilter1d_pot_kernel_d;
cl_kernel magicfilter1d_t_kernel_d;
cl_kernel magicfiltershrink1d_kernel_d;
cl_kernel magicfiltergrow1d_kernel_d;
cl_kernel magicfiltergrow1d_pot_kernel_d;
cl_program magicfilterProgram;

void create_magicfilter_kernels(){
    cl_int ciErrNum = CL_SUCCESS;
    magicfiltergrow1d_kernel_d=clCreateKernel(magicfilterProgram,"magicfiltergrow1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfiltergrow1dKernel_d kernel!");
    magicfiltershrink1d_kernel_d=clCreateKernel(magicfilterProgram,"magicfiltershrink1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfiltershrink1dKernel_d kernel!");
    magicfiltergrow1d_pot_kernel_d=clCreateKernel(magicfilterProgram,"magicfiltergrow1d_potKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfiltergrow1d_potKernel_d kernel!");
    magicfilter1d_kernel_d=clCreateKernel(magicfilterProgram,"magicfilter1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfilter1dKernel_d kernel!");
    magicfilter1d_den_kernel_d=clCreateKernel(magicfilterProgram,"magicfilter1d_denKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfilter1d_denKernel_d kernel!");
    magicfilter1d_pot_kernel_d=clCreateKernel(magicfilterProgram,"magicfilter1d_potKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfilter1d_potKernel_d kernel!");
    magicfilter1d_t_kernel_d=clCreateKernel(magicfilterProgram,"magicfilter1d_tKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create magicfilter1d_tKernel_d kernel!");
}

void build_magicfilter_programs(cl_context * context){
    cl_int ciErrNum = CL_SUCCESS;
    magicfilterProgram = clCreateProgramWithSource(*context,1,(const char**) &magicfilter_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(magicfilterProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to build magicfilter program!\n",ciErrNum);
        char cBuildLog[10240];
        clGetProgramBuildInfo(magicfilterProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
}

void FC_FUNC_(magicfiltershrink1d_d,MAGICFILTERSHRINK1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    magicfilter_generic(magicfiltershrink1d_kernel_d,command_queue,n,ndat,psi,out);
}

void FC_FUNC_(magicfiltergrow1d_d,MAGICFILTERGROW1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_uint n1 = *n + 15;
    magicfilter_generic(magicfiltergrow1d_kernel_d, command_queue, &n1, ndat, psi, out);
}

void FC_FUNC_(magicfilter1d_d,MAGICFILTER1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, n, ndat, psi, out);
}

void FC_FUNC_(magicfilter1d_pot_d,MAGICFILTER1D_POT_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *pot, cl_mem *out){
    magicfilter_pot_generic(magicfilter1d_pot_kernel_d, command_queue, n, ndat, psi, pot, out);
}


void FC_FUNC_(magicfilter1d_t_d,MAGICFILTER1D_T_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    magicfilter_generic(magicfilter1d_t_kernel_d,command_queue,n,ndat,psi,out);
}

void FC_FUNC_(magicfilter_n_self_d,MAGICFILTER_N_SELF_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *psi, cl_mem *out){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, &n3, &ndat, psi, out);
    ndat = n1 * n3;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, &n2, &ndat, out, psi);
    ndat = n2 * n3;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, &n1, &ndat, psi, out);
}

void FC_FUNC_(magicfilter_n_d,MAGICFILTER_N_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, &n3, &ndat, psi, out);
    ndat = n1 * n3;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, &n2, &ndat, out, tmp);
    ndat = n2 * n3;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, &n1, &ndat, tmp, out);
}

void FC_FUNC_(magicfilter_den_d,MAGICFILTER_DEN_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out){
    cl_uint n1 = dimensions[0]*2;
    cl_uint n2 = dimensions[1]*2;
    cl_uint n3 = dimensions[2]*2;
    cl_uint ndat = n1 * n2;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, &n3, &ndat, psi, out);
    ndat = n1 * n3;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, &n2, &ndat, out, tmp);
    ndat = n2 * n3;
    magicfilter_generic(magicfilter1d_den_kernel_d, command_queue, &n1, &ndat, tmp, out);
}

void FC_FUNC_(magicfilter_t_self_d,MAGICFILTER_T_SELF_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *psi, cl_mem *out){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, &n3, &ndat, psi, out);
    ndat = n1 * n3;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, &n2, &ndat, out, psi);
    ndat = n2 * n3;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, &n1, &ndat, psi, out);
}

void FC_FUNC_(magicfilter_t_d,MAGICFILTER_T_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, &n3, &ndat, psi, out);
    ndat = n1 * n3;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, &n2, &ndat, out, tmp);
    ndat = n2 * n3;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, &n1, &ndat, tmp, out);
}

void FC_FUNC_(potential_application_d_generic,POTENTIAL_APPLICATION_D_GENERIC)(cl_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *psi, cl_mem *out, cl_mem *pot) {
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
    magicfilter_generic(magicfilter1d_kernel_d, command_queue,  &n3, &ndat, psi, tmp);
  } else {
    n3 += 15;
    magicfilter_generic(magicfiltergrow1d_kernel_d, command_queue, &n3, &ndat, psi, tmp);
  }
  ndat = n1 * n3;
  if( periodic[1] ) {
    magicfilter_generic(magicfilter1d_kernel_d, command_queue,  &n2, &ndat, tmp, out);
  } else {
    n2 += 15;
    magicfilter_generic(magicfiltergrow1d_kernel_d, command_queue, &n2, &ndat, tmp, out);
  }
  ndat =  n2 * n3;
  if( periodic[0] ) {
    magicfilter_pot_generic(magicfilter1d_pot_kernel_d, command_queue,  &n1, &ndat, out, pot, tmp);
  } else {
    n1 += 15;
    magicfilter_pot_generic(magicfiltergrow1d_pot_kernel_d, command_queue, &n1, &ndat, out, pot, tmp);
  }
  ndat = n1 * n2;
  if( periodic[2] ) {
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue,  &n3, &ndat, tmp, out);
  } else {
    n3 -= 15;
    magicfilter_generic(magicfiltershrink1d_kernel_d, command_queue,  &n3, &ndat, tmp, out);
  }
  ndat = n1 * n3;
  if( periodic[1] ) {
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue,  &n2, &ndat, out, tmp);
  } else {
    n2 -= 15;
    magicfilter_generic(magicfiltershrink1d_kernel_d, command_queue,  &n2, &ndat, out, tmp);
  }
  ndat = n2 * n3;
  if( periodic[0] ) {
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue,  &n1, &ndat, tmp, out);
  } else {
    n1 -= 15;
    magicfilter_generic(magicfiltershrink1d_kernel_d, command_queue,  &n1, &ndat, tmp, out);
  }
}

void FC_FUNC_(potential_application_d,POTENTIAL_APPLICATION_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out, cl_mem *pot) {
    cl_uint n1 = dimensions[0] * 2;
    cl_uint n2 = dimensions[1] * 2;
    cl_uint n3 = dimensions[2] * 2;
    cl_uint ndat = n1 * n2;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, &n3, &ndat, psi, tmp);
    ndat = n1 * n3;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, &n2, &ndat, tmp, out);
    ndat = n2 * n3;
    magicfilter1d_pot_d_(command_queue, &n1, &ndat, out, pot, tmp);
    ndat = n1 * n2;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, &n3, &ndat, tmp, out);
    ndat = n1 * n3;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, &n2, &ndat, out, tmp);
    ndat = n2 * n3;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, &n1, &ndat, tmp, out);
}

void clean_magicfilter_kernels(){
  cl_int ciErrNum;
  ciErrNum = clReleaseKernel(magicfilter1d_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(magicfilter1d_den_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(magicfilter1d_pot_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(magicfilter1d_t_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(magicfiltershrink1d_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(magicfiltergrow1d_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(magicfiltergrow1d_pot_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseProgram(magicfilterProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
}
