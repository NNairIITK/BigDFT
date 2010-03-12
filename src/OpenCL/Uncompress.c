#include "OpenCL_wrappers.h"
#include "Uncompress.h"

cl_kernel uncompress_coarse_kernel_d;
cl_kernel uncompress_fine_kernel_d;
cl_kernel compress_coarse_kernel_d;
cl_kernel compress_fine_kernel_d;

void build_uncompress_kernels(cl_context * context){
    cl_int ciErrNum = CL_SUCCESS;
    cl_program uncompressProgram = clCreateProgramWithSource(*context,1,(const char**) &uncompress_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(uncompressProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build uncompress program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(uncompressProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    uncompress_coarse_kernel_d=clCreateKernel(uncompressProgram,"uncompress_coarseKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    uncompress_fine_kernel_d=clCreateKernel(uncompressProgram,"uncompress_fineKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = clReleaseProgram(uncompressProgram);
    oclErrorCheck(ciErrNum,"Failed to release program!");

    cl_program compressProgram = clCreateProgramWithSource(*context,1,(const char**) &compress_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(compressProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build compress program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(compressProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    compress_coarse_kernel_d=clCreateKernel(compressProgram,"compress_coarseKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    compress_fine_kernel_d=clCreateKernel(compressProgram,"compress_fineKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = clReleaseProgram(compressProgram);
    oclErrorCheck(ciErrNum,"Failed to release program!");
}

void FC_FUNC_(uncompress_d,UNCOMPRESS_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,
                                       cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                       cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
				       cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi_out) {
    cl_uint full_size = *n1 * *n2 * *n3 * 8;
    double init = 0.0;
    v_initialize_generic(v_initialize_kernel_d, command_queue, &full_size, psi_out, &init);
    uncompress_coarse_generic(uncompress_coarse_kernel_d, command_queue, n1, n2, n3, nseg_c, nvctr_c, keyg_c, keyv_c, psi_c, psi_out);
    if(nvctr_f == 0) return;
    uncompress_fine_generic(uncompress_fine_kernel_d, command_queue, n1, n2, n3, nseg_f, nvctr_f, keyg_f, keyv_f, psi_f, psi_out);
}

void FC_FUNC_(compress_d,COMPRESS_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,
                                     cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                     cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                     cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi) {
    compress_coarse_generic(compress_coarse_kernel_d, command_queue, n1, n2, n3, nseg_c, nvctr_c, keyg_c, keyv_c, psi_c, psi);
    if(nvctr_f == 0) return;
    compress_fine_generic(compress_fine_kernel_d, command_queue, n1, n2, n3, nseg_f, nvctr_f, keyg_f, keyv_f, psi_f, psi);
}

void clean_uncompress_kernels(){
  clReleaseKernel(uncompress_coarse_kernel_d);
  clReleaseKernel(uncompress_fine_kernel_d);
  clReleaseKernel(compress_coarse_kernel_d);
  clReleaseKernel(compress_fine_kernel_d);
}
