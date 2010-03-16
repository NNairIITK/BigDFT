#include "OpenCL_wrappers.h"

char * reduction_program="\
//group_size is supposed to be 512\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
__kernel void reductionKernel_d( size_t n, __global const double *x, __global double *y, __local double *tmp ) {\n\
  size_t i = get_local_id(0);\n\
  size_t g = get_group_id(0)*1024+i;\n\
  if(g<n) {\n\
    tmp[i] = x[g];\n\
  } else {\n\
    tmp[i] = 0.0;\n\
  }\n\
  if(g+512<n) {\n\
    tmp[i+512] = x[g+512];\n\
  } else {\n\
    tmp[i+512] = 0.0;\n\
  }\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  tmp[i] = tmp[i] + tmp[i+512];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<256 )\n\
    tmp[i] = tmp[i] + tmp[i+256];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<128 )\n\
    tmp[i] = tmp[i] + tmp[i+128];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<64 )\n\
    tmp[i] = tmp[i] + tmp[i+64];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<32 )\n\
    tmp[i] = tmp[i] + tmp[i+32];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<16 )\n\
    tmp[i] = tmp[i] + tmp[i+16];\n\
  if( i<8 )\n\
    tmp[i] = tmp[i] + tmp[i+8];\n\
  if( i<4 )\n\
    tmp[i] = tmp[i] + tmp[i+4];\n\
  if( i<2 )\n\
    tmp[i] = tmp[i] + tmp[i+2];\n\
  if( i==0 )\n\
    y[get_group_id(0)] = tmp[0]+tmp[1];\n\
}\n\
__kernel void reduction_dotKernel_d( size_t n, __global const double *x, __global double *y, __local double *tmp ) {\n\
  size_t i = get_local_id(0);\n\
  size_t g = get_group_id(0)*1024+i;\n\
  double tt;\n\
  if(g<n) {\n\
    tt = x[g];\n\
    tmp[i] = tt*tt;\n\
  } else {\n\
    tmp[i] = 0.0;\n\
  }\n\
  if(g+512<n) {\n\
    tt = x[g+512];\n\
    tmp[i+512] = tt*tt;\n\
  } else {\n\
    tmp[i+512] = 0.0;\n\
  }\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  tmp[i] = tmp[i] + tmp[i+512];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<256 )\n\
    tmp[i] = tmp[i] + tmp[i+256];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<128 )\n\
    tmp[i] = tmp[i] + tmp[i+128];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<64 )\n\
    tmp[i] = tmp[i] + tmp[i+64];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<32 )\n\
    tmp[i] = tmp[i] + tmp[i+32];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if( i<16 )\n\
    tmp[i] = tmp[i] + tmp[i+16];\n\
  if( i<8 )\n\
    tmp[i] = tmp[i] + tmp[i+8];\n\
  if( i<4 )\n\
    tmp[i] = tmp[i] + tmp[i+4];\n\
  if( i<2 )\n\
    tmp[i] = tmp[i] + tmp[i+2];\n\
  if( i==0 )\n\
    y[get_group_id(0)] = tmp[0]+tmp[1];\n\
}\n\
";

void inline reduction_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *out) {
  cl_int ciErrNum;
  size_t block_size_i=512;
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
  clSetKernelArg(kernel, i++,sizeof(*in), (void*)in);
  clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
  clSetKernelArg(kernel, i++,sizeof(double)*block_size_i*2, NULL);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i*2,*ndat)/2 };
  ciErrNum = clEnqueueNDRangeKernel(*command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue reduction kernel!");
}

cl_kernel reduction_kernel_d;
cl_kernel reduction_dot_kernel_d;

void FC_FUNC_(reduction_self_d,REDUCTION_SELF_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work, double *out) {
  assert(*ndat>0);
  cl_uint n = *ndat;
  cl_mem *input = in;
  cl_mem *output = work;
  cl_mem *tmp;
  do {
    reduction_generic(reduction_kernel_d, command_queue, &n, input, output);
    tmp = input;
    input = output;
    output = tmp;
    n = shrRoundUp(1024,n)/1024;
  } while(n>1);
  clEnqueueReadBuffer(*command_queue, *input, CL_TRUE, 0, sizeof(double), out, 0, NULL, NULL);
}

void FC_FUNC_(reduction_dot_self_d,REDUCTION_SELF_DOT_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work, double *out) {
  assert(*ndat>0);
  cl_uint n = *ndat;
  cl_mem *input = in;
  cl_mem *output = work;
  cl_mem *tmp;
  reduction_generic(reduction_dot_kernel_d, command_queue, &n, input, output);
  input = work;
  output = in;
  n = shrRoundUp(1024,n)/1024;
  if(n>1) {
    do {
      reduction_generic(reduction_kernel_d, command_queue, &n, input, output);
      tmp = input;
      input = output;
      output = tmp;
      n = shrRoundUp(1024,n)/1024;
    } while(n>1);
  }
  clEnqueueReadBuffer(*command_queue, *input, CL_TRUE, 0, sizeof(double), out, 0, NULL, NULL);
}

void FC_FUNC_(reduction_d,REDUCTION_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work1, cl_mem *work2, double *out) {
  assert(*ndat>0);
  cl_uint n = *ndat;
  cl_mem *input = in;
  cl_mem *output = work1;
  cl_mem *tmp;
  reduction_generic(reduction_kernel_d, command_queue, &n, input, output);
  input = work1;
  output = work2;
  n = shrRoundUp(1024,n)/1024;
  if(n>1) {
    do {
      reduction_generic(reduction_kernel_d, command_queue, &n, input, output);
      tmp = input;
      input = output;
      output = tmp;
      n = shrRoundUp(1024,n)/1024;
    } while(n>1);
  }
  clEnqueueReadBuffer(*command_queue, *input, CL_TRUE, 0, sizeof(double), out, 0, NULL, NULL);
}

void FC_FUNC_(reduction_dot_d,REDUCTION_DOT_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work1, cl_mem *work2, double *out) {
  assert(*ndat>0);
  cl_uint n = *ndat;
  cl_mem *input = in;
  cl_mem *output = work1;
  cl_mem *tmp;
  reduction_generic(reduction_dot_kernel_d, command_queue, &n, input, output);
  input = work1;
  output = work2;
  n = shrRoundUp(1024,n)/1024;
  if(n>1) {
    do {
      reduction_generic(reduction_kernel_d, command_queue, &n, input, output);
      tmp = input;
      input = output;
      output = tmp;
      n = shrRoundUp(1024,n)/1024;
    } while(n>1);
  }
  clEnqueueReadBuffer(*command_queue, *input, CL_TRUE, 0, sizeof(double), out, 0, NULL, NULL);
}

void build_reduction_kernels(cl_context * context){
    cl_int ciErrNum = CL_SUCCESS;

    cl_program reductionProgram = clCreateProgramWithSource(*context,1,(const char**) &reduction_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(reductionProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build reduction program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(reductionProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
        fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    reduction_kernel_d=clCreateKernel(reductionProgram,"reductionKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    reduction_dot_kernel_d=clCreateKernel(reductionProgram,"reduction_dotKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = clReleaseProgram(reductionProgram);
    oclErrorCheck(ciErrNum,"Failed to release program!");
}

void clean_reduction_kernels(){
  clReleaseKernel(reduction_kernel_d);
  clReleaseKernel(reduction_dot_kernel_d);
}
