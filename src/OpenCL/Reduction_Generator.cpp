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

static void generate_header(std::stringstream &program){
  program<<"#ifdef cl_khr_fp64\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#elif defined (cl_amd_fp64)\n\
#pragma OPENCL EXTENSION cl_amd_fp64: enable \n\
#endif\n";
}

static void generate_reductionKernel(std::stringstream &program,struct bigdft_device_infos * infos){
  size_t max_wgs = infos->MAX_WORK_GROUP_SIZE;
  cl_device_type device_type = infos->DEVICE_TYPE;
  size_t cutoff;
  if( infos->DEVICE_TYPE == CL_DEVICE_TYPE_CPU )
    cutoff = 4;
  else
    cutoff = 64;
  program<<"__kernel void reductionKernel_d( uint n, __global const double *x, __global double *y, __local volatile double *tmp ) {\n\
  //get our position in the local buffer\n\
  size_t i = get_local_id(0);\n\
  //get our position in the input data\n\
  size_t g = get_group_id(0)*"<<max_wgs*2<<"+i;\n\
  //copy our 2 data elements in the buffer, or store 0 if out of bounds.\n\
  if(g<n)\n\
    tmp[i] = x[g];\n\
  else\n\
    tmp[i] = 0.0;\n\
  if(g+"<<max_wgs<<"<n)\n\
    tmp[i+"<<max_wgs<<"] = x[g+"<<max_wgs<<"];\n\
  else\n\
    tmp[i+"<<max_wgs<<"] = 0.0;\n\
  //wait for buffer to be full\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  tmp[i] = tmp[i] + tmp[i+"<<max_wgs<<"];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n";
  do {
    max_wgs /= 2;
    program<<"  if( i<"<<max_wgs<<" )\n\
    tmp[i] = tmp[i] + tmp[i+"<<max_wgs<<"];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n";
  } while(max_wgs >= cutoff);
  while(max_wgs >= 4) {
    max_wgs /= 2;
    program<<"  if( i<"<<max_wgs<<" )\n\
    tmp[i] = tmp[i] + tmp[i+"<<max_wgs<<"];\n";
  }
  program<<"  if( i==0 )\n\
    y[get_group_id(0)] = tmp[0]+tmp[1];\n\
}\n";
}

static void generate_reduction_dotKernel(std::stringstream &program,struct bigdft_device_infos * infos){
  size_t max_wgs = infos->MAX_WORK_GROUP_SIZE;
  size_t cutoff;
  if( infos->DEVICE_TYPE == CL_DEVICE_TYPE_CPU )
    cutoff = 4;
  else
    cutoff = 64;
  program<<"__kernel void reduction_dotKernel_d( uint n, __global const double *x, __global double *y, __local volatile double *tmp ) {\n\
  size_t i = get_local_id(0);\n\
  size_t g = get_group_id(0)*"<<max_wgs*2<<"+i;\n\
  double tt;\n\
  if(g<n) {\n\
    tt = x[g];\n\
    tmp[i] = tt*tt;\n\
  } else\n\
    tmp[i] = 0.0;\n\
  if(g+"<<max_wgs<<"<n) {\n\
    tt = x[g+"<<max_wgs<<"];\n\
    tmp[i+"<<max_wgs<<"] = tt*tt;\n\
  } else\n\
    tmp[i+"<<max_wgs<<"] = 0.0;\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  tmp[i] = tmp[i] + tmp[i+"<<max_wgs<<"];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n";
  do {
    max_wgs /= 2;
    program<<"  if( i<"<<max_wgs<<" )\n\
    tmp[i] = tmp[i] + tmp[i+"<<max_wgs<<"];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n";
  } while(max_wgs >= cutoff);
  while(max_wgs >= 4) {
    max_wgs /= 2;
    program<<"  if( i<"<<max_wgs<<" )\n\
    tmp[i] = tmp[i] + tmp[i+"<<max_wgs<<"];\n";
  }
  program<<"  if( i==0 )\n\
    y[get_group_id(0)] = tmp[0]+tmp[1];\n\
}\n";
}

static void generate_dotKernel(std::stringstream &program,struct bigdft_device_infos * infos){
  size_t max_wgs = infos->MAX_WORK_GROUP_SIZE;
  size_t cutoff;
  if( infos->DEVICE_TYPE == CL_DEVICE_TYPE_CPU )
    cutoff = 4;
  else
    cutoff = 64;
  program<<"__kernel void dotKernel_d( uint n, __global const double *x, __global const double *y, __global double *z, __local volatile double *tmp ) {\n\
  size_t i = get_local_id(0);\n\
  size_t g = get_group_id(0)*"<<max_wgs*2<<"+i;\n\
  if(g<n)\n\
    tmp[i] = x[g]*y[g];\n\
  else\n\
    tmp[i] = 0.0;\n\
  if(g+"<<max_wgs<<"<n)\n\
    tmp[i+"<<max_wgs<<"] = x[g+"<<max_wgs<<"]*y[g+"<<max_wgs<<"];\n\
  else\n\
    tmp[i+"<<max_wgs<<"] = 0.0;\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  tmp[i] = tmp[i] + tmp[i+"<<max_wgs<<"];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n";
  do {
    max_wgs /= 2;
    program<<"  if( i<"<<max_wgs<<" )\n\
    tmp[i] = tmp[i] + tmp[i+"<<max_wgs<<"];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n";
  } while(max_wgs >= cutoff);
  while(max_wgs >= 4) {
    max_wgs /= 2;
    program<<"  if( i<"<<max_wgs<<" )\n\
    tmp[i] = tmp[i] + tmp[i+"<<max_wgs<<"];\n";
  }
  program<<"  if( i==0 )\n\
    z[get_group_id(0)] = tmp[0]+tmp[1];\n\
}\n";
}

static void generate_axpyKernel(std::stringstream &program){
  program<<"__kernel void axpyKernel_d( uint n, double alpha, __global const double *x, __global const double *y, __global double *out) {\n\
  size_t ig = get_global_id(0);\n\
  if( ig < n)\n\
    out[ig] = y[ig] + alpha * x[ig];\n\
}\n";
}

static void generate_axpy_offsetKernel(std::stringstream &program){
  program<<"__kernel void axpy_offsetKernel_d( uint n, double alpha, uint offset_x, __global const double *x, uint offset_y, __global double *y, uint offset_out, __global double *out) {\n\
  size_t ig = get_global_id(0);\n\
  if( ig < n)\n\
    out[ig+offset_out] = y[ig+offset_y] + alpha * x[ig+offset_x];\n\
}\n";
}

static void generate_scalKernel(std::stringstream &program){
  program<<"__kernel void scalKernel_d( uint n, double alpha, __global const double *x, __global double *y) {\n\
  size_t ig = get_global_id(0);\n\
  if( ig < n)\n\
    y[ig] = alpha * x[ig];\n\
}\n";
}

static void generate_copyKernel(std::stringstream &program){
  program<<"__kernel void copyKernel_d( uint n, __global const double *x, __global double *y) {\n\
  size_t ig = get_global_id(0);\n\
  if( ig < n)\n\
    y[ig] = x[ig];\n\
}\n";
}

static void generate_setKernel(std::stringstream &program){
  program<<"__kernel void setKernel_d( uint n, const double val, __global double *x) {\n\
  size_t ig = get_global_id(0);\n\
  if( ig < n)\n\
    x[ig] = val;\n\
}\n";
}

extern "C" char* generate_reduction_program(struct bigdft_device_infos * infos){
  char * output;
  std::stringstream program;

  generate_header(program);

  generate_reductionKernel(program,infos);
  generate_reduction_dotKernel(program,infos);
  generate_dotKernel(program,infos);
  generate_axpyKernel(program);
  generate_axpy_offsetKernel(program);
  generate_scalKernel(program);
  generate_copyKernel(program);
  generate_setKernel(program);

  output = (char *)malloc((program.str().size()+1)*sizeof(char));
  strcpy(output, program.str().c_str());
  return output;
}

