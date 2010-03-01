#include "OpenCL_wrappers.h"
#include "Initialize.h"

cl_kernel c_initialize_kernel_d;
cl_kernel v_initialize_kernel_d;


void build_initialize_kernels(cl_context * context){
    cl_int ciErrNum = CL_SUCCESS;

    cl_program c_initializeProgram = clCreateProgramWithSource(*context,1,(const char**) &c_initialize_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(c_initializeProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build c_initialize program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(c_initializeProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    c_initialize_kernel_d=clCreateKernel(c_initializeProgram,"c_initializeKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    v_initialize_kernel_d=clCreateKernel(c_initializeProgram,"v_initializeKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = clReleaseProgram(c_initializeProgram);
    oclErrorCheck(ciErrNum,"Failed to release program!");
}

void clean_initialize_kernels(){
  clReleaseKernel(c_initialize_kernel_d);
  clReleaseKernel(v_initialize_kernel_d);
}
