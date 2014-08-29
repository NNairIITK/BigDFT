#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include "Kinetic_Generator.h"
#include "OpenCL_wrappers.h"

#define FILTER_WIDTH 32
static void generate_header(std::stringstream &program){
  program<<"#ifdef cl_khr_fp64\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#elif defined (cl_amd_fp64)\n\
#pragma OPENCL EXTENSION cl_amd_fp64: enable \n\
#endif\n\
#define FILTER_WIDTH "<<FILTER_WIDTH<<"\n\
#define FILT1_14 -6.924474940639200152025730585882e-18\n\
#define FILT1_13  2.70800493626319438269856689037647576e-13\n\
#define FILT1_12 -5.813879830282540547959250667e-11\n\
#define FILT1_11 -1.05857055496741470373494132287e-8\n\
#define FILT1_10 -3.7230763047369275848791496973044e-7\n\
#define FILT1_9   2.0904234952920365957922889447361e-6\n\
#define FILT1_8  -0.2398228524507599670405555359023135e-4\n\
#define FILT1_7   0.45167920287502235349480037639758496e-3\n\
#define FILT1_6  -0.409765689342633823899327051188315485e-2\n\
#define FILT1_5   0.02207029188482255523789911295638968409e0\n\
#define FILT1_4  -0.0822663999742123340987663521e0\n\
#define FILT1_3   0.2371780582153805636239247476e0\n\
#define FILT1_2  -0.6156141465570069496314853949e0\n\
#define FILT1_1   2.2191465938911163898794546405e0\n\
#define FILT1_0  -3.5536922899131901941296809374e0\n";
}

static void generate_filters(std::stringstream &program){
  program<<"#define filter1(tt,tmp) \
tt = mad(tmp[14] + tmp[-14], FILT1_14, tt);\
tt = mad(tmp[13] + tmp[-13], FILT1_13, tt);\
tt = mad(tmp[12] + tmp[-12], FILT1_12, tt);\
tt = mad(tmp[11] + tmp[-11], FILT1_11, tt);\
tt = mad(tmp[10] + tmp[-10], FILT1_10, tt);\
tt = mad(tmp[ 9] + tmp[ -9], FILT1_9 , tt);\
tt = mad(tmp[ 8] + tmp[ -8], FILT1_8 , tt);\
tt = mad(tmp[ 7] + tmp[ -7], FILT1_7 , tt);\
tt = mad(tmp[ 6] + tmp[ -6], FILT1_6 , tt);\
tt = mad(tmp[ 5] + tmp[ -5], FILT1_5 , tt);\
tt = mad(tmp[ 4] + tmp[ -4], FILT1_4 , tt);\
tt = mad(tmp[ 3] + tmp[ -3], FILT1_3 , tt);\
tt = mad(tmp[ 2] + tmp[ -2], FILT1_2 , tt);\
tt = mad(tmp[ 1] + tmp[ -1], FILT1_1 , tt);\
tt = mad(tmp[ 0]           , FILT1_0 , tt);\n";
}

static void generate_kinetic1dKernel(std::stringstream &program, struct bigdft_device_infos * infos){
  size_t max_wgs = infos->MAX_WORK_GROUP_SIZE;
  size_t buff_l = max_wgs / FILTER_WIDTH;
  if(buff_l > 16)
    buff_l = 16;
  program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<buff_l<<", 1))) void kinetic1dKernel_d(uint n, uint ndat, double scale, __global const double *  restrict x_in, __global double *  restrict x_out, __global const double *  restrict y_in, __global double *  restrict y_out, __local double * tmp, __local double * tmp2) {\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  const size_t i2 = get_local_id(0);\n\
  const size_t j2 = get_local_id(1);\n\
  ptrdiff_t igt = get_group_id(0);\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  size_t jb;\n\
  const size_t j2t = "<<FILTER_WIDTH/buff_l<<"*j2 + i2/"<<buff_l<<";\n\
  const size_t i2t = i2 - "<<buff_l<<" * (i2 / "<<buff_l<<");\n\
  //if data are ill dimentioned last block recomputes part of the data\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
  igt = ig - i2 + j2t;\n\
  jgt = jg - j2 + i2t;\n\
  tmp2[i2t * (32 + 1) + j2t] = y_in[jgt+igt*ndat];\n\
  igt -= FILTER_WIDTH/2;\n\
  if ( igt < 0 ) \n\
    tmp[i2t * (2 * FILTER_WIDTH + 1) + j2t] = x_in[jgt + ( n + igt ) * ndat];\n\
  else \n\
    tmp[i2t * (2 * FILTER_WIDTH + 1) + j2t] = x_in[jgt + igt * ndat];\n\
  igt += FILTER_WIDTH;\n\
  if ( igt >= n ) \n\
    tmp[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = x_in[jgt + ( igt - n ) * ndat];\n\
  else\n\
    tmp[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = x_in[jgt +  igt * ndat];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  __local double * tmp_o = tmp + j2*(2*32 + 1) + FILTER_WIDTH/2+i2;\n\
  double conv = 0.0;\n\
  filter1(conv, tmp_o);\n\
  //y_out[jg*n + ig] = y_in[jg + ig*ndat] - scale * conv;\n\
  y_out[jg*n + ig] = tmp2[j2*(32+1) + i2] + scale * conv;\n\
  x_out[jg*n + ig] = tmp_o[0];\n\
};\n";
}

static void generate_kinetic1d_fKernel(std::stringstream &program, struct bigdft_device_infos * infos){
  size_t max_wgs = infos->MAX_WORK_GROUP_SIZE;
  size_t buff_l = max_wgs / FILTER_WIDTH;
  if(buff_l > 16)
    buff_l = 16;
  program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<buff_l<<", 1))) void kinetic1d_fKernel_d(uint n, uint ndat, double scale, __global const double *  restrict x_in, __global double *  restrict x_out, __global const double *  restrict y_in, __global double *  restrict y_out, __local double * tmp, __local double * tmp2) {\n\
  size_t ig = get_global_id(0);\n\
  size_t jg = get_global_id(1);\n\
  const size_t i2 = get_local_id(0);\n\
  const size_t j2 = get_local_id(1);\n\
  ptrdiff_t igt = get_group_id(0);\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  size_t jb;\n\
  const size_t j2t = "<<FILTER_WIDTH/buff_l<<"*j2 + i2/"<<buff_l<<";\n\
  const size_t i2t = i2 - "<<buff_l<<" * (i2 / "<<buff_l<<");\n\
  //if data are ill dimentioned last block recomputes part of the data\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
  igt = ig - i2 + j2t;\n\
  jgt = jg - j2 + i2t;\n\
  tmp2[i2t * (32 + 1) + j2t] = y_in[jgt+igt*ndat];\n\
  igt -= FILTER_WIDTH/2;\n\
  if ( igt < 0 ) \n\
    tmp[i2t * (2 * FILTER_WIDTH + 1) + j2t] = 0.0;\n\
  else \n\
    tmp[i2t * (2 * FILTER_WIDTH + 1) + j2t] = x_in[jgt + igt * ndat];\n\
  igt += FILTER_WIDTH;\n\
  if ( igt >= n ) \n\
    tmp[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = 0.0;\n\
  else\n\
    tmp[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = x_in[jgt +  igt * ndat];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  __local double * tmp_o = tmp + j2*(2*32 + 1) + FILTER_WIDTH/2+i2;\n\
  double conv = 0.0;\n\
  filter1(conv, tmp_o);\n\
  //y_out[jg*n + ig] = y_in[jg + ig*ndat] - scale * conv;\n\
  y_out[jg*n + ig] = tmp2[j2*(32+1) + i2] + scale * conv;\n\
  x_out[jg*n + ig] = tmp_o[0];\n\
};\n";
}

extern "C" char* generate_kinetic_program(struct bigdft_device_infos * infos){
  char * output;
  std::stringstream program;

  generate_header(program);
  generate_filters(program);
  generate_kinetic1dKernel(program, infos);
  generate_kinetic1d_fKernel(program, infos);

  output = (char *)malloc((program.str().size()+1)*sizeof(char));
  strcpy(output, program.str().c_str());
  return output;
}

