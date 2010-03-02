#include "OpenCL_wrappers.h"
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
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %lu, dimension dat: %lu, psi: %p, out: %p\n",*command_queue, (long unsigned)*n, (long unsigned)*ndat, *psi, *out);
#endif
    int FILTER_WIDTH = 16;
    if(*n<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small : n = %d!\n", __func__, __FILE__,*n); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;
    cl_uint i = 0;
    clSetKernelArg(anashrink1d_kernel_d, i++,sizeof(*n), (void*)n);
    clSetKernelArg(anashrink1d_kernel_d, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(anashrink1d_kernel_d, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(anashrink1d_kernel_d, i++,sizeof(*out), (void*)out);
    clSetKernelArg(anashrink1d_kernel_d, i++,sizeof(double)*block_size_j*(block_size_i*2+FILTER_WIDTH + 1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, anashrink1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue anashrink1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }

}

void FC_FUNC_(ana1d_d,ANA1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %lu, dimension dat: %lu, psi: %p, out: %p\n",*command_queue, (long unsigned)*n, (long unsigned)*ndat, *psi, *out);
#endif
    int FILTER_WIDTH = 16;
    if(*n<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;
    cl_uint i = 0;
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*n), (void*)n);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*out), (void*)out);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(double)*block_size_j*(block_size_i*2+FILTER_WIDTH + 1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, ana1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue ana1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }

}

void FC_FUNC_(ana_d,ANA_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,cl_mem *tmp,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    int FILTER_WIDTH = 16;
    if(*n1<FILTER_WIDTH || *n2<FILTER_WIDTH || *n3<FILTER_WIDTH) {fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;
    cl_uint ndat = *n2 * *n1 * 4;
    cl_uint i = 0;
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*n3), (void*)n3);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*out), (void*)out);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(double)*block_size_j*(block_size_i*2+FILTER_WIDTH + 1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n3), shrRoundUp(block_size_j,ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, ana1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue ana1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }
    ndat = *n1 * *n3 * 4;
    i = 0;
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*n2), (void*)n2);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*out), (void*)out);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*tmp), (void*)tmp);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(double)*block_size_j*(block_size_i*2+FILTER_WIDTH + 1), 0);
    localWorkSize[0] =  block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,*n2); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, ana1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue ana1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }
    ndat = *n2 * *n3 * 4;
    i = 0;
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*n1), (void*)n1);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*tmp), (void*)tmp);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*out), (void*)out);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(double)*block_size_j*(block_size_i*2+FILTER_WIDTH + 1), 0);
    localWorkSize[0] =  block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,*n1); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, ana1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue ana1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }
}

void FC_FUNC_(ana_self_d,ANA_SELF_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    int FILTER_WIDTH = 16;
    if(*n1<FILTER_WIDTH || *n2<FILTER_WIDTH || *n3<FILTER_WIDTH) {fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;
    cl_uint ndat = *n2 * *n1 * 4;
    cl_uint i = 0;
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*n3), (void*)n3);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*out), (void*)out);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(double)*block_size_j*(block_size_i*2+FILTER_WIDTH + 1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n3), shrRoundUp(block_size_j,ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, ana1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue ana1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }
    ndat = *n1 * *n3 * 4;
    i = 0;
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*n2), (void*)n2);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*out), (void*)out);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(double)*block_size_j*(block_size_i*2+FILTER_WIDTH + 1), 0);
    localWorkSize[0] =  block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,*n2); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, ana1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue ana1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }
    ndat = *n2 * *n3 * 4;
    i = 0;
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*n1), (void*)n1);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(*out), (void*)out);
    clSetKernelArg(ana1d_kernel_d, i++,sizeof(double)*block_size_j*(block_size_i*2+FILTER_WIDTH + 1), 0);
    localWorkSize[0] =  block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,*n1); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, ana1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue ana1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }
}

void FC_FUNC_(syngrow1d_d,SYNGROW1D_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %lu, dimension dat: %lu, psi: %p, out: %p\n",*command_queue, (long unsigned)*n, (long unsigned)*ndat, *psi, *out);
#endif
    int FILTER_WIDTH = 8;
    cl_uint n1 = *n+7;
    if(n1<16) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=16, block_size_j=16;
    cl_uint i = 0;
    clSetKernelArg(syngrow1d_kernel_d, i++,sizeof(n1), (void*)&n1);
    clSetKernelArg(syngrow1d_kernel_d, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(syngrow1d_kernel_d, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(syngrow1d_kernel_d, i++,sizeof(*out), (void*)out);
    clSetKernelArg(syngrow1d_kernel_d, i++,sizeof(double)*block_size_j*(block_size_i*2+FILTER_WIDTH*2+1), NULL);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,n1), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, syngrow1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {         
        fprintf(stderr,"Error %d: Failed to enqueue syngrow1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }  
}

void FC_FUNC_(syn1d_d,SYN1D_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %lu, dimension dat: %lu, psi: %p, out: %p\n",*command_queue, (long unsigned)*n, (long unsigned)*ndat, *psi, *out);
#endif
    int FILTER_WIDTH = 8;
    if(*n<16) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=16, block_size_j=16;
    cl_uint i = 0;
    clSetKernelArg(syn1d_kernel_d, i++,sizeof(*n), (void*)n);
    clSetKernelArg(syn1d_kernel_d, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(syn1d_kernel_d, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(syn1d_kernel_d, i++,sizeof(*out), (void*)out);
    clSetKernelArg(syn1d_kernel_d, i++,sizeof(double)*block_size_j*(block_size_i*2+FILTER_WIDTH*2+1), NULL);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, syn1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {         
        fprintf(stderr,"Error %d: Failed to enqueue syn1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }  
}

void clean_wavelet_kernels(){
  clReleaseKernel(ana1d_kernel_d);
  clReleaseKernel(anashrink1d_kernel_d);
  clReleaseKernel(syn1d_kernel_d);
  clReleaseKernel(syngrow1d_kernel_d);
}
