#include "Wavelet.h"

cl_kernel ana1d_kernel_d;
cl_kernel anashrink1d_kernel_d;
cl_kernel syn1d_kernel_d;
cl_kernel syngrow1d_kernel_d;


void build_wavelet_kernels(cl_context * context){
    cl_int ciErrNum = CL_SUCCESS;
    cl_program ana1dProgram = clCreateProgramWithSource(*context,1,(const char**) &ana1d_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(ana1dProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build ana1d program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(ana1dProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    anashrink1d_kernel_d=clCreateKernel(ana1dProgram,"anashrink1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    ana1d_kernel_d=clCreateKernel(ana1dProgram,"ana1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = clReleaseProgram(ana1dProgram);
    oclErrorCheck(ciErrNum,"Failed to release program!");

    cl_program syn1dProgram = clCreateProgramWithSource(*context,1,(const char**) &syn1d_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(syn1dProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build syn1d program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(syn1dProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    syngrow1d_kernel_d=clCreateKernel(syn1dProgram,"syngrow1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    syn1d_kernel_d=clCreateKernel(syn1dProgram,"syn1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = clReleaseProgram(syn1dProgram);
    oclErrorCheck(ciErrNum,"Failed to release program!");

}

void FC_FUNC_(anashrink1d_d,ANASHRINK1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
  ana_generic(anashrink1d_kernel_d, command_queue, n, ndat, psi, out);
}

void FC_FUNC_(ana1d_d,ANA1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
  ana_generic(ana1d_kernel_d, command_queue, n, ndat, psi, out);
}

void FC_FUNC_(ana_d,ANA_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,cl_mem *tmp,cl_mem *psi,cl_mem *out){
    cl_uint ndat = *n2 * *n1 * 4;
    ana_generic(ana1d_kernel_d, command_queue, n3, &ndat, psi, out);
    ndat = *n1 * *n3 * 4;
    ana_generic(ana1d_kernel_d, command_queue, n2, &ndat, out, tmp);
    ndat = *n2 * *n3 * 4;
    ana_generic(ana1d_kernel_d, command_queue, n1, &ndat, tmp, out);
}

void FC_FUNC_(ana_self_d,ANA_SELF_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,cl_mem *psi,cl_mem *out){
    cl_uint ndat = *n2 * *n1 * 4;
    ana_generic(ana1d_kernel_d, command_queue, n3, &ndat, psi, out);
    ndat = *n1 * *n3 * 4;
    ana_generic(ana1d_kernel_d, command_queue, n2, &ndat, out, psi);
    ndat = *n2 * *n3 * 4;
    ana_generic(ana1d_kernel_d, command_queue, n1, &ndat, psi, out);
}

void FC_FUNC_(syngrow1d_d,SYNGROW1D_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_uint n1 = *n+7;
    syn_generic(syngrow1d_kernel_d, command_queue, &n1, ndat, psi, out);
}

void FC_FUNC_(syn1d_d,SYN1D_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat,cl_mem *psi,cl_mem *out){
    syn_generic(syn1d_kernel_d, command_queue, n, ndat, psi, out);
}

void FC_FUNC_(syn_d,SYN_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,cl_mem *tmp, cl_mem *psi, cl_mem *out){
    cl_uint ndat = *n2 * *n1 * 4;
    syn_generic(syn1d_kernel_d, command_queue, n3, &ndat, psi, out);
    ndat = *n1 * *n3 * 4;
    syn_generic(syn1d_kernel_d, command_queue, n2, &ndat, out, tmp);
    ndat = *n2 * *n3 * 4;
    syn_generic(syn1d_kernel_d, command_queue, n1, &ndat, tmp, out);
}
void FC_FUNC_(syn_self_d,SYN_SELF_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3, cl_mem *psi, cl_mem *out){
    cl_uint ndat = *n2 * *n1 * 4;
    syn_generic(syn1d_kernel_d, command_queue, n3, &ndat, psi, out);
    ndat = *n1 * *n3 * 4;
    syn_generic(syn1d_kernel_d, command_queue, n2, &ndat, out, psi);
    ndat = *n2 * *n3 * 4;
    syn_generic(syn1d_kernel_d, command_queue, n1, &ndat, psi, out);
}

void clean_wavelet_kernels(){
  clReleaseKernel(ana1d_kernel_d);
  clReleaseKernel(anashrink1d_kernel_d);
  clReleaseKernel(syn1d_kernel_d);
  clReleaseKernel(syngrow1d_kernel_d);
}
