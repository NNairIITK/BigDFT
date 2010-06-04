#include "OpenCL_wrappers.h"

char * benchmark_program = "\
#pragma OPENCL EXTENSION cl_khr_fp64: enable\n\
#define NB_ITER 16\n\
__kernel void benchmark_flopsKernel_d(uint n, __global const double *in, __global double *out){\n\
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
};\n\
";

cl_kernel benchmark_flops_kernel_d;
cl_program benchmarkProgram;

void create_benchmark_kernels(){
    cl_int ciErrNum = CL_SUCCESS;
    benchmark_flops_kernel_d=clCreateKernel(benchmarkProgram,"benchmark_flopsKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create benchmark_flopsKernel_d kernel!");
}

void build_benchmark_programs(cl_context * context){
    cl_int ciErrNum = CL_SUCCESS;
    benchmarkProgram = clCreateProgramWithSource(*context,1,(const char**) &benchmark_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(benchmarkProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to build benchmark program!\n",ciErrNum);
        char cBuildLog[10240];
        clGetProgramBuildInfo(benchmarkProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
}

inline void benchmark_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n,cl_mem *in,cl_mem *out){
    cl_int ciErrNum;
    int FILTER_WIDTH=32;
    assert(*n>=FILTER_WIDTH);
    size_t block_size_i=FILTER_WIDTH;
    cl_uint i = 0;
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*in), (void*)in);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    size_t localWorkSize[] = { block_size_i };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue benchmark kernel!");
}

void FC_FUNC_(benchmark_flops_d,BENCHMARK_FLOPS_D)(cl_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *out){
    benchmark_generic(benchmark_flops_kernel_d, command_queue, n, in, out);
}

void clean_benchmark_kernels(){
  cl_int ciErrNum;
  ciErrNum = clReleaseKernel(benchmark_flops_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseProgram(benchmarkProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
}
