#include "OpenCL_wrappers.h"
#include "Wavelet.h"

cl_kernel ana1d_kernel_l;
cl_kernel syn1d_kernel_l;


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
    ana1d_kernel_l=clCreateKernel(ana1dProgram,"ana1dKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");

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
    syn1d_kernel_l=clCreateKernel(syn1dProgram,"syn1dKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");

}


void FC_FUNC_(ana1d_l,ANA1D_L)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %lu, dimension dat: %lu, psi: %p, out: %p\n",*command_queue, (long unsigned)*n, (long unsigned)*ndat, *psi, *out);
#endif
    int FILTER_WIDTH = 16;
    if(*n<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=256/FILTER_WIDTH;
    while (*n > block_size_i >= 1 && block_size_j > 8)
        { block_size_i *= 2; block_size_j /= 2;}
    cl_uint i = 0;
    clSetKernelArg(ana1d_kernel_l, i++,sizeof(*n), (void*)n);
    clSetKernelArg(ana1d_kernel_l, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(ana1d_kernel_l, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(ana1d_kernel_l, i++,sizeof(*out), (void*)out);
    clSetKernelArg(ana1d_kernel_l, i++,sizeof(float)*block_size_j*(block_size_i*2+FILTER_WIDTH), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, ana1d_kernel_l, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue ana1d_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }

}

void FC_FUNC_(syn1d_l,SYN1D_L)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %lu, dimension dat: %lu, psi: %p, out: %p\n",*command_queue, (long unsigned)*n, (long unsigned)*ndat, *psi, *out);
#endif
    int FILTER_WIDTH = 16;
    if(*n<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=256/FILTER_WIDTH;
    while (*n > block_size_i >= 1 && block_size_j > 8)
        { block_size_i *= 2; block_size_j /= 2;}
    cl_uint i = 0;
    clSetKernelArg(syn1d_kernel_l, i++,sizeof(*n), (void*)n);
    clSetKernelArg(syn1d_kernel_l, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(syn1d_kernel_l, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(syn1d_kernel_l, i++,sizeof(*out), (void*)out);
    clSetKernelArg(syn1d_kernel_l, i++,sizeof(float)*block_size_j*(block_size_i+FILTER_WIDTH), NULL);
    clSetKernelArg(syn1d_kernel_l, i++,sizeof(float)*block_size_j*(block_size_i+FILTER_WIDTH), NULL);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, syn1d_kernel_l, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {         
        fprintf(stderr,"Error %d: Failed to enqueue syn1d_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }  
}


