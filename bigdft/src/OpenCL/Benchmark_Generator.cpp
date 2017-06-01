#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include "Benchmark_Generator.h"
#include "OpenCL_wrappers.h"

#define FILTER_WIDTH 16
static void generate_header(std::stringstream &program){
  program<<"#ifdef cl_khr_fp64\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#elif defined (cl_amd_fp64)\n\
#pragma OPENCL EXTENSION cl_amd_fp64: enable \n\
#endif\n\
#define FILTER_WIDTH "<<FILTER_WIDTH<<"\n";
}

static void generate_mopsKernel(std::stringstream &program){
  program<<"#define NB_ITER 16\n\
__kernel void benchmark_mopsKernel_d(uint n, __global const double *in, __global double *out){\n\
size_t i = get_global_id(0);\n\
i = get_group_id(0) == get_num_groups(0) - 1 ? i - ( get_global_size(0) - n ) : i;\n\
out[i] = in[i];\n\
out[i] = in[i];\n\
out[i] = in[i];\n\
out[i] = in[i];\n\
out[i] = in[i];\n\
out[i] = in[i];\n\
out[i] = in[i];\n\
out[i] = in[i];\n\
};\n";
}

static void generate_flopsKernel(std::stringstream &program){
  program<<"__kernel void benchmark_flopsKernel_d(uint n, __global const double *in, __global double *out){\n\
size_t i = get_global_id(0);\n\
double a = in[i]*1.15;\n\
double b = in[i]*1.16;\n\
i = get_group_id(0) == get_num_groups(0) - 1 ? i - ( get_global_size(0) - n ) : i;\n\
int j=0;\n\
for(j=0;j<NB_ITER;j++){\n\
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
     a = b * a + b; \
     b = a * b + a; \
}\n\
out[i]=a+b;\n\
};\n";
}

static void generate_transposeKernel(std::stringstream &program){
  program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<", 1))) void transposeKernel_d(uint n, uint ndat, __global const double *psi, __global double *out, __local double *tmp ) {\n\
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
tmp[i2 * (FILTER_WIDTH + 1) + j2] = psi[jgt + igt * ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
out[(jg*n+ig)]=tmp[j2 * (FILTER_WIDTH + 1) + i2];\n\
};\n";
}

static void generate_notransposeKernel(std::stringstream &program){
  program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<FILTER_WIDTH<<", 1))) void notransposeKernel_d(uint n, uint ndat, __global const double *psi, __global double *out, __local double *tmp ) {\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1);\n\
ptrdiff_t igt = get_group_id(0);\n\
ptrdiff_t jgt = get_group_id(1);\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
//igt = ig - i2 + j2;\n\
//jgt = jg - j2 + i2;\n\
//If I'm on the outside, select a border element to load\n\
//tmp[i2 * (FILTER_WIDTH + 1) + j2] = psi[jgt + igt * ndat];\n\
//barrier(CLK_LOCAL_MEM_FENCE);\n\
\
out[jg*n+ig]=psi[jg*n+ig];//tmp[j2 * (FILTER_WIDTH + 1) + i2];\n\
};\n";
}

extern "C" char* generate_benchmark_program(struct bigdft_device_infos * infos){
  char * output;
  std::stringstream program;

  generate_header(program);
  generate_mopsKernel(program);
  generate_flopsKernel(program);
  generate_transposeKernel(program);
  generate_notransposeKernel(program);

  output = (char *)malloc((program.str().size()+1)*sizeof(char));
  strcpy(output, program.str().c_str());
  return output;
}

