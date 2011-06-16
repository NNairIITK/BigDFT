#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include "Reduction_Generator.h"
#include "OpenCL_wrappers.h"

#define FILTER_WIDTH 32
static void generate_header(std::stringstream &program){
  program<<"#ifdef cl_khr_fp64\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#elif defined (cl_amd_fp64)\n\
#pragma OPENCL EXTENSION cl_amd_fp64: enable \n\
#endif\n\
#define FILTER_WIDTH "<<FILTER_WIDTH<<"\n";
}

static void generate_kinetic1dKernel(std::stringstream &program, struct bigdft_device_infos * infos){
  size_t max_wgs = infos->MAX_WORK_GROUP_SIZE;
  size_t buff_l = max_wgs / FILTER_WIDTH;
  if(buff_l > 16)
    buff_l = 16;
  program<<"__kernel void kinetic1dKernel_d(uint n, uint ndat, double scale, __global const double * x_in, __global double * x_out, __global const double * y_in, __global double * y_out, __local double * tmp, __local double * tmp2) {\n\
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
  conv += (tmp_o[14] + tmp_o[-14]) * -6.924474940639200152025730585882e-18;\n\
  conv += (tmp_o[13] + tmp_o[-13]) *  2.70800493626319438269856689037647576e-13;\n\
  conv += (tmp_o[12] + tmp_o[-12]) * -5.813879830282540547959250667e-11;\n\
  conv += (tmp_o[11] + tmp_o[-11]) * -1.05857055496741470373494132287e-8;\n\
  conv += (tmp_o[10] + tmp_o[-10]) * -3.7230763047369275848791496973044e-7;\n\
  conv += (tmp_o[ 9] + tmp_o[ -9]) *  2.0904234952920365957922889447361e-6;\n\
  conv += (tmp_o[ 8] + tmp_o[ -8]) * -0.2398228524507599670405555359023135e-4;\n\
  conv += (tmp_o[ 7] + tmp_o[ -7]) *  0.45167920287502235349480037639758496e-3;\n\
  conv += (tmp_o[ 6] + tmp_o[ -6]) * -0.409765689342633823899327051188315485e-2;\n\
  conv += (tmp_o[ 5] + tmp_o[ -5]) *  0.02207029188482255523789911295638968409e0;\n\
  conv += (tmp_o[ 4] + tmp_o[ -4]) * -0.0822663999742123340987663521e0;\n\
  conv += (tmp_o[ 3] + tmp_o[ -3]) *  0.2371780582153805636239247476e0;\n\
  conv += (tmp_o[ 2] + tmp_o[ -2]) * -0.6156141465570069496314853949e0;\n\
  conv += (tmp_o[ 1] + tmp_o[ -1]) *  2.2191465938911163898794546405e0;\n\
  conv +=  tmp_o[ 0]             * -3.5536922899131901941296809374e0;\n\
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
  program<<"__kernel void kinetic1d_fKernel_d(uint n, uint ndat, double scale, __global const double * x_in, __global double * x_out, __global const double * y_in, __global double * y_out, __local double * tmp, __local double * tmp2) {\n\
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
  conv += (tmp_o[14] + tmp_o[-14]) * -6.924474940639200152025730585882e-18;\n\
  conv += (tmp_o[13] + tmp_o[-13]) *  2.70800493626319438269856689037647576e-13;\n\
  conv += (tmp_o[12] + tmp_o[-12]) * -5.813879830282540547959250667e-11;\n\
  conv += (tmp_o[11] + tmp_o[-11]) * -1.05857055496741470373494132287e-8;\n\
  conv += (tmp_o[10] + tmp_o[-10]) * -3.7230763047369275848791496973044e-7;\n\
  conv += (tmp_o[ 9] + tmp_o[ -9]) *  2.0904234952920365957922889447361e-6;\n\
  conv += (tmp_o[ 8] + tmp_o[ -8]) * -0.2398228524507599670405555359023135e-4;\n\
  conv += (tmp_o[ 7] + tmp_o[ -7]) *  0.45167920287502235349480037639758496e-3;\n\
  conv += (tmp_o[ 6] + tmp_o[ -6]) * -0.409765689342633823899327051188315485e-2;\n\
  conv += (tmp_o[ 5] + tmp_o[ -5]) *  0.02207029188482255523789911295638968409e0;\n\
  conv += (tmp_o[ 4] + tmp_o[ -4]) * -0.0822663999742123340987663521e0;\n\
  conv += (tmp_o[ 3] + tmp_o[ -3]) *  0.2371780582153805636239247476e0;\n\
  conv += (tmp_o[ 2] + tmp_o[ -2]) * -0.6156141465570069496314853949e0;\n\
  conv += (tmp_o[ 1] + tmp_o[ -1]) *  2.2191465938911163898794546405e0;\n\
  conv +=  tmp_o[ 0]             * -3.5536922899131901941296809374e0;\n\
  //y_out[jg*n + ig] = y_in[jg + ig*ndat] - scale * conv;\n\
  y_out[jg*n + ig] = tmp2[j2*(32+1) + i2] + scale * conv;\n\
  x_out[jg*n + ig] = tmp_o[0];\n\
};\n";
}

extern "C" char* generate_kinetic_program(struct bigdft_device_infos * infos){
  char * output;
  std::stringstream program;

  generate_header(program);
  generate_kinetic1dKernel(program, infos);
  generate_kinetic1d_fKernel(program, infos);

  output = (char *)malloc((program.str().size()+1)*sizeof(char));
  strcpy(output, program.str().c_str());
  return output;
}

