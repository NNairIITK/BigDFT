#include "OpenCL_wrappers.h"
#include "Initialize.h"

char * initialize_program="\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
__kernel void c_initializeKernel_d(uint n, __global const double * x_in, __global double * y_in, double c) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
y_in[ig] = x_in[ig] * c;\n\
};\n\
__kernel void v_initializeKernel_d(uint n, __global double * y_in, double v) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
y_in[ig] = v;\n\
};\n\
__kernel void p_initializeKernel_d(uint n, __global const double * x, __global double * y) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
y[ig] = x[ig] * x[ig];\n\
};\n\
";

cl_kernel c_initialize_kernel_d;
cl_kernel v_initialize_kernel_d;
cl_kernel p_initialize_kernel_d;
cl_program initializeProgram;

void inline p_initialize_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *out) {
  cl_int ciErrNum;
  size_t block_size_i=64;
  assert(*ndat>=block_size_i);
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
  clSetKernelArg(kernel, i++,sizeof(*in), (void*)in);
  clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*ndat) };
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue p_initialize kernel!");
}

void inline c_initialize_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *inout, double *c) {
  cl_int ciErrNum;
  size_t block_size_i=64;
  assert(*ndat>=block_size_i);
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
  clSetKernelArg(kernel, i++,sizeof(*in), (void*)in);
  clSetKernelArg(kernel, i++,sizeof(*inout), (void*)inout);
  clSetKernelArg(kernel, i++,sizeof(*c), (void*)c);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*ndat) };
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue c_initialize kernel!");
}

void inline v_initialize_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *ndat, cl_mem *out, double *c) {
  cl_int ciErrNum;
  size_t block_size_i=64;
  assert(*ndat>=block_size_i);
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
  clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
  clSetKernelArg(kernel, i++,sizeof(*c), (void*)c);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*ndat) };
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue v_initialize kernel!");
}

void create_initialize_kernels(){
    cl_int ciErrNum = CL_SUCCESS;
    c_initialize_kernel_d=clCreateKernel(initializeProgram,"c_initializeKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create c_initializeKernel_d kernel!");
    v_initialize_kernel_d=clCreateKernel(initializeProgram,"v_initializeKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create v_initializeKernel_d kernel!");
    p_initialize_kernel_d=clCreateKernel(initializeProgram,"p_initializeKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create p_initializeKernel_d kernel!");
}

void build_initialize_programs(cl_context * context){
    cl_int ciErrNum = CL_SUCCESS;

    initializeProgram = clCreateProgramWithSource(*context,1,(const char**) &initialize_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(initializeProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build c_initialize program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(initializeProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
}

void clean_initialize_kernels(){
  cl_int ciErrNum;
  ciErrNum = clReleaseKernel(c_initialize_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(v_initialize_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(p_initialize_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseProgram(initializeProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
}