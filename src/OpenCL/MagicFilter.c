#include "MagicFilter.h"

cl_kernel magicfilter1d_kernel_d;
cl_kernel magicfilter1d_pot_kernel_d;
cl_kernel magicfilter1d_t_kernel_d;
cl_kernel magicfiltershrink1d_kernel_d;
cl_kernel magicfiltergrow1d_kernel_d;

void build_magicfilter_kernels(cl_context * context){
    cl_int ciErrNum = CL_SUCCESS;
    cl_program magicfilter1dProgram = clCreateProgramWithSource(*context,1,(const char**) &magicfilter1d_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(magicfilter1dProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to build magicfilter1d program!\n",ciErrNum);
        char cBuildLog[10240];
        clGetProgramBuildInfo(magicfilter1dProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    magicfiltergrow1d_kernel_d=clCreateKernel(magicfilter1dProgram,"magicfiltergrow1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    magicfiltershrink1d_kernel_d=clCreateKernel(magicfilter1dProgram,"magicfiltershrink1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    magicfilter1d_kernel_d=clCreateKernel(magicfilter1dProgram,"magicfilter1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    magicfilter1d_pot_kernel_d=clCreateKernel(magicfilter1dProgram,"magicfilter1d_potKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    magicfilter1d_t_kernel_d=clCreateKernel(magicfilter1dProgram,"magicfilter1d_tKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = clReleaseProgram(magicfilter1dProgram);
    oclErrorCheck(ciErrNum,"Failed to release program!");
}

void FC_FUNC_(magicfiltershrink1d_d,MAGICFILTERSHRINK1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    magicfilter_generic(magicfiltershrink1d_kernel_d,command_queue,n,ndat,psi,out);
}

void FC_FUNC_(magicfiltergrow1d_d,MAGICFILTERGROW1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_uint n1 = *n + 15;
    magicfilter_generic(magicfiltergrow1d_kernel_d,command_queue,&n1,ndat,psi,out);
}

void FC_FUNC_(magicfilter1d_d,MAGICFILTER1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    magicfilter_generic(magicfilter1d_kernel_d,command_queue,n,ndat,psi,out);
}

void FC_FUNC_(magicfilter1d_pot_d,MAGICFILTER1D_POT_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *pot, cl_mem *out){
    cl_int ciErrNum;
    cl_event e;
    int FILTER_WIDTH = 16;
    assert(*n>=FILTER_WIDTH);
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;
    cl_uint i = 0;
    ciErrNum = clSetKernelArg(magicfilter1d_pot_kernel_d, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(magicfilter1d_pot_kernel_d, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_pot_kernel_d, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(magicfilter1d_pot_kernel_d, i++,sizeof(*pot), (void*)pot);
    ciErrNum = clSetKernelArg(magicfilter1d_pot_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_pot_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_pot_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &e);
#if PROFILING
    event ev;
    ev.e = e;
    ev.comment = __func__;
    addToEventList(ev);
#endif
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_pot_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }
}


void FC_FUNC_(magicfilter1d_t_d,MAGICFILTER1D_T_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    magicfilter_generic(magicfilter1d_t_kernel_d,command_queue,n,ndat,psi,out);
}

void FC_FUNC_(magicfilter_n_self_d,MAGICFILTER_N_SELF_D)(cl_command_queue *command_queue, cl_uint *n1,cl_uint *n2,cl_uint *n3,cl_mem *psi,cl_mem *out){
    cl_uint ndat = *n1 * *n2;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, n3, &ndat, psi, out);
    ndat = *n1 * *n3;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, n2, &ndat, out, psi);
    ndat = *n2 * *n3;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, n1, &ndat, psi, out);
}

void FC_FUNC_(magicfilter_n_d,MAGICFILTER_N_D)(cl_command_queue *command_queue, cl_uint *n1,cl_uint *n2,cl_uint *n3,cl_mem *tmp,cl_mem *psi,cl_mem *out){
    cl_uint ndat = *n1 * *n2;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, n3, &ndat, psi, out);
    ndat = *n1 * *n3;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, n2, &ndat, out, tmp);
    ndat = *n2 * *n3;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, n1, &ndat, tmp, out);
}

void FC_FUNC_(magicfilter_t_self_d,MAGICFILTER_T_SELF_D)(cl_command_queue *command_queue, cl_uint *n1,cl_uint *n2,cl_uint *n3,cl_mem *psi,cl_mem *out){
    cl_uint ndat = *n1 * *n2;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, n3, &ndat, psi, out);
    ndat = *n1 * *n3;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, n2, &ndat, out, psi);
    ndat = *n2 * *n3;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, n1, &ndat, psi, out);
}

void FC_FUNC_(magicfilter_t_d,MAGICFILTER_T_D)(cl_command_queue *command_queue, cl_uint *n1,cl_uint *n2,cl_uint *n3,cl_mem *tmp,cl_mem *psi,cl_mem *out){
    cl_uint ndat = *n1 * *n2;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, n3, &ndat, psi, out);
    ndat = *n1 * *n3;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, n2, &ndat, out, tmp);
    ndat = *n2 * *n3;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, n1, &ndat, tmp, out);
}

void FC_FUNC_(potential_application_d,POTENTIAL_APPLICATION_D)(cl_command_queue *command_queue, cl_uint *n1,cl_uint *n2,cl_uint *n3, cl_mem *tmp, cl_mem *psi, cl_mem *out, cl_mem *pot) {
    cl_uint n = *n3 * 2;
    cl_uint ndat = *n1 * *n2 * 4;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, &n, &ndat, psi, tmp);
    n = *n2 * 2;
    ndat = *n1 * *n3 * 4;
    magicfilter_generic(magicfilter1d_kernel_d, command_queue, &n, &ndat, tmp, out);
    n = *n1 * 2;
    ndat = *n2 * *n3 * 4;
    magicfilter1d_pot_d_(command_queue, &n, &ndat, out, pot, tmp);
    n = *n3 * 2;
    ndat = *n1 * *n2 * 4;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, &n, &ndat, tmp, out);
    n = *n2 * 2;
    ndat = *n1 * *n3 * 4;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, &n, &ndat, out, tmp);
    n = *n1 * 2;
    ndat = *n2 * *n3 * 4;
    magicfilter_generic(magicfilter1d_t_kernel_d, command_queue, &n, &ndat, tmp, out);
}

void clean_magicfilter_kernels(){
  clReleaseKernel(magicfilter1d_kernel_d);
  clReleaseKernel(magicfilter1d_pot_kernel_d);
  clReleaseKernel(magicfilter1d_t_kernel_d);
  clReleaseKernel(magicfiltershrink1d_kernel_d);
  clReleaseKernel(magicfiltergrow1d_kernel_d);
}
