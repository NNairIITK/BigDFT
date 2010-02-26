#ifndef MAGICFILTER_H
#define MAGICFILTER_H

char * magicfilter1d_program="\
#define FILTER_WIDTH 16\n\
//n is supposed to be greater or equal than get_local_size(0)\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
__kernel void magicfiltergrow1dKernel_d(size_t n, size_t ndat, __global const double *psi, __global double *out, __local double tmp[]){\n\
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
//If I'm on the outside, select a border element to load\n\
if(j2 < FILTER_WIDTH/2)\n\
  { if (igt < FILTER_WIDTH/2)\n\
      { tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = 0.0; }\n\
    else {\n\
      jb = igt - FILTER_WIDTH/2 - 7;\n\
      if (jb < 0) { \n\
        tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = 0.0;\n\
      } else {\n\
        tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = psi[jgt+jb*ndat];\n\
      }\n\
    }\n\
  }\n\
if (j2 >= FILTER_WIDTH/2)\n\
  { if (igt >= n - FILTER_WIDTH/2)\n\
      { tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = 0.0; }\n\
    else {\n\
      jb = igt + FILTER_WIDTH/2 - 7;\n\
      if (igt + FILTER_WIDTH/2  >= n - 8) {//if (jb >= n-8 ) {\n\
        tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = 0.0;\n\
      } else { \n\
        tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = psi[jgt + jb * ndat];\n\
      }\n\
    }\n\
  }\n\
//check boundaries\
//Load the element I am to calculate\n\
if(igt < 7 || igt >= n - 8) {\n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH/2]=0.0;\n\
} else {\n\
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH/2]=psi[jgt + (igt - 7) * ndat];\n\
}\n\
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
__kernel void magicfiltershrink1dKernel_d(size_t n, size_t ndat, __global const double *psi, __global double *out, __local double tmp[]){\n\
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
__kernel void magicfilter1dKernel_d_ref(size_t n, size_t ndat, __global const double *psi, __global double *out, __local double tmp[]){\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
size_t i2 = get_local_id(0);\n\
size_t j2 = get_local_id(1);\n\
size_t is = get_local_size(0);\n\
size_t ib;\n\
size_t it;\n\
size_t base_i;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
if( jg >= ndat ) return;\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
__local double * tmp_o = tmp + j2*(is+FILTER_WIDTH);\n\
base_i = FILTER_WIDTH/2+i2;\n\
//If I'm on the outside, select a border element to load\n\
if(i2 < FILTER_WIDTH/2)\n\
  { it = i2;\n\
    if (ig < FILTER_WIDTH/2)\n\
      { ib =  n + i2 - FILTER_WIDTH/2; }\n\
    else { ib = ig - FILTER_WIDTH/2; }\n\
    tmp_o[it]=psi[jg+ib*ndat];\n\
  }\n\
if (i2 >= (is - FILTER_WIDTH/2) || (ig >= n - FILTER_WIDTH/2))\n\
  { it = i2 + FILTER_WIDTH;\n\
    if (ig >= n - FILTER_WIDTH/2)\n\
      { ib = ig - n + FILTER_WIDTH/2; }\n\
    else { ib = ig + FILTER_WIDTH/2; }\n\
    tmp_o[it]=psi[jg+ib*ndat];\n\
  }\n\
//check boundaries\
//Load the element I am to calculate\n\
tmp_o[base_i]=psi[jg+ig*ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
double tt = 0.0;\n\
tmp_o = tmp_o + base_i;\n\
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
__kernel void magicfilter1dKernel_d(size_t n, size_t ndat, __global const double *psi, __global double *out, __local double tmp[]){\n\
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
//If I'm on the outside, select a border element to load\n\
if(j2 < FILTER_WIDTH/2)\n\
  { if (igt < FILTER_WIDTH/2)\n\
      { jb =  n + j2 - FILTER_WIDTH/2; }\n\
    else { jb = igt - FILTER_WIDTH/2; }\n\
    tmp[i2 * (2 * FILTER_WIDTH + 1) + j2]=psi[jgt + jb * ndat];\n\
  }\n\
if (j2 >= FILTER_WIDTH/2)\n\
  { if (igt >= n - FILTER_WIDTH/2)\n\
      { jb = igt - n + FILTER_WIDTH/2; }\n\
    else { jb = igt + FILTER_WIDTH/2; }\n\
    tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH]=psi[jgt + jb * ndat];\n\
  }\n\
//check boundaries\
//Load the element I am to calculate\n\
tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH/2]=psi[jgt + igt * ndat];\n\
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
\n\
__kernel void magicfilter1dKernel_l(size_t n, size_t ndat, __global const float *psi, __global float *out, __local float *tmp){\n\
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
//If I'm on the outside, select a border element to load\n\
if(j2 < FILTER_WIDTH/2)\n\
  { if (igt < FILTER_WIDTH/2)\n\
      { jb =  n + j2 - FILTER_WIDTH/2; }\n\
    else { jb = igt - FILTER_WIDTH/2; }\n\
    tmp[i2 * (2 * FILTER_WIDTH + 1) + j2]=psi[jgt + jb * ndat];\n\
  }\n\
if (j2 >= FILTER_WIDTH/2)\n\
  { if (igt >= n - FILTER_WIDTH/2)\n\
      { jb = igt - n + FILTER_WIDTH/2; }\n\
    else { jb = igt + FILTER_WIDTH/2; }\n\
    tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH]=psi[jgt + jb * ndat];\n\
  }\n\
//check boundaries\
//Load the element I am to calculate\n\
tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH/2]=psi[jgt + igt*ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
float tt = 0.0;\n\
__local float * tmp_o = tmp + j2*(2*FILTER_WIDTH+1) + FILTER_WIDTH/2+i2;\n\
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
\n\
//straight version\n\
__kernel void magicfilter1dKernel_s_l(size_t n, size_t ndat, __global const float *psi, __global float *out, __local float tmp[]){\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
size_t i2 = get_local_id(0);\n\
size_t j2 = get_local_id(1);\n\
size_t is = get_local_size(0);\n\
size_t ib;\n\
size_t it;\n\
size_t base_i;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
if( jg >= ndat ) return;\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
__local float * tmp_o = tmp + j2*(is+FILTER_WIDTH);\n\
base_i = FILTER_WIDTH/2+i2;\n\
//If I'm on the outside, select a border element to load\n\
if(i2 < FILTER_WIDTH/2)\n\
  { it = i2;\n\
    if (ig < FILTER_WIDTH/2)\n\
      { ib =  n + i2 - FILTER_WIDTH/2; }\n\
    else { ib = ig - FILTER_WIDTH/2; }\n\
    tmp_o[it]=psi[jg*n+ib];\n\
  }\n\
if (i2 >= (is - FILTER_WIDTH/2) || (ig >= n - FILTER_WIDTH/2))\n\
  { it = i2 + FILTER_WIDTH;\n\
    if (ig >= n - FILTER_WIDTH/2)\n\
      { ib = ig - n + FILTER_WIDTH/2; }\n\
    else { ib = ig + FILTER_WIDTH/2; }\n\
    tmp_o[it]=psi[jg*n+ib];\n\
  }\n\
//check boundaries\
//Load the element I am to calculate\n\
tmp_o[base_i]=psi[jg*n+ig];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
float tt = 0.0;\n\
tmp_o = tmp_o + base_i;\n\
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
\n\
//filter multiplied by a potential\n\
__kernel void magicfilter1d_potKernel_l(size_t n, size_t ndat, __global const float *psi, float pot, __global float *out, __local float tmp[]){\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
size_t i2 = get_local_id(0);\n\
size_t j2 = get_local_id(1);\n\
size_t is = get_local_size(0);\n\
size_t ib;\n\
size_t it;\n\
size_t base_i;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
if( jg >= ndat ) return;\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
__local float * tmp_o = tmp + j2*(is+FILTER_WIDTH);\n\
base_i = FILTER_WIDTH/2+i2;\n\
//If I'm on the outside, select a border element to load\n\
if(i2 < FILTER_WIDTH/2)\n\
  { it = i2;\n\
    if (ig < FILTER_WIDTH/2)\n\
      { ib =  n + i2 - FILTER_WIDTH/2; }\n\
    else { ib = ig - FILTER_WIDTH/2; }\n\
    tmp_o[it]=psi[jg+ib*ndat];\n\
  }\n\
if (i2 >= (is - FILTER_WIDTH/2) || (ig >= n - FILTER_WIDTH/2))\n\
  { it = i2 + FILTER_WIDTH;\n\
    if (ig >= n - FILTER_WIDTH/2)\n\
      { ib = ig - n + FILTER_WIDTH/2; }\n\
    else { ib = ig + FILTER_WIDTH/2; }\n\
    tmp_o[it]=psi[jg+ib*ndat];\n\
  }\n\
//check boundaries\
//Load the element I am to calculate\n\
tmp_o[base_i]=psi[jg+ig*ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
float tt = 0.0;\n\
tmp_o = tmp_o + base_i;\n\
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
out[(jg*n+ig)]=tt*pot;\n\
};\
\n\
//transposed kernel\n\
__kernel void magicfilter1d_tKernel_l(size_t n, size_t ndat, __global const float *psi, __global float *out, __local float tmp[]){\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
size_t i2 = get_local_id(0);\n\
size_t j2 = get_local_id(1);\n\
size_t is = get_local_size(0);\n\
size_t ib;\n\
size_t it;\n\
size_t base_i;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
if( jg >= ndat ) return;\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
__local float * tmp_o = tmp + j2*(is+FILTER_WIDTH);\n\
base_i = FILTER_WIDTH/2+i2;\n\
//If I'm on the outside, select a border element to load\n\
if(i2 < FILTER_WIDTH/2)\n\
  { it = i2;\n\
    if (ig < FILTER_WIDTH/2)\n\
      { ib =  n + i2 - FILTER_WIDTH/2; }\n\
    else { ib = ig - FILTER_WIDTH/2; }\n\
    tmp_o[it]=psi[jg+ib*ndat];\n\
  }\n\
if (i2 >= (is - FILTER_WIDTH/2) || (ig >= n - FILTER_WIDTH/2))\n\
  { it = i2 + FILTER_WIDTH;\n\
    if (ig >= n - FILTER_WIDTH/2)\n\
      { ib = ig - n + FILTER_WIDTH/2; }\n\
    else { ib = ig + FILTER_WIDTH/2; }\n\
    tmp_o[it]=psi[jg+ib*ndat];\n\
  }\n\
//check boundaries\
//Load the element I am to calculate\n\
tmp_o[base_i]=psi[jg+ig*ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
float tt = 0.0;\n\
tmp_o = tmp_o + base_i;\n\
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
\n\
//filter and square of the psi function\n\
__kernel void magicfilter1d_denKernel_l(size_t n, size_t ndat, __global const float *psi, __global float *out, __local float tmp[]){\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
size_t i2 = get_local_id(0);\n\
size_t j2 = get_local_id(1);\n\
size_t is = get_local_size(0);\n\
size_t ib;\n\
size_t it;\n\
size_t base_i;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
if( jg >= ndat ) return;\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
__local float * tmp_o = tmp + j2*(is+FILTER_WIDTH);\n\
base_i = FILTER_WIDTH/2+i2;\n\
//If I'm on the outside, select a border element to load\n\
if(i2 < FILTER_WIDTH/2)\n\
  { it = i2;\n\
    if (ig < FILTER_WIDTH/2)\n\
      { ib =  n + i2 - FILTER_WIDTH/2; }\n\
    else { ib = ig - FILTER_WIDTH/2; }\n\
    tmp_o[it]=psi[jg+ib*ndat];\n\
  }\n\
if (i2 >= (is - FILTER_WIDTH/2) || (ig >= n - FILTER_WIDTH/2))\n\
  { it = i2 + FILTER_WIDTH;\n\
    if (ig >= n - FILTER_WIDTH/2)\n\
      { ib = ig - n + FILTER_WIDTH/2; }\n\
    else { ib = ig + FILTER_WIDTH/2; }\n\
    tmp_o[it]=psi[jg+ib*ndat];\n\
  }\n\
//check boundaries\
//Load the element I am to calculate\n\
tmp_o[base_i]=psi[jg+ig*ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
float tt = 0.0;\n\
tmp_o = tmp_o + base_i;\n\
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
";

#endif
