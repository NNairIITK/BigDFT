#include "OpenCL_wrappers.h"
#include "Kinetic.h"

cl_kernel kinetic1d_kernel_l;

void build_kinetic_kernels(cl_context * context){
    cl_int ciErrNum=CL_SUCCESS;

    cl_program kinetic1dProgram = clCreateProgramWithSource(*context,1,(const char**) &kinetic1d_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(kinetic1dProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build kinetic1d program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(kinetic1dProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    kinetic1d_kernel_l=clCreateKernel(kinetic1dProgram,"kinetic1dKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = clReleaseProgram(kinetic1dProgram);
    oclErrorCheck(ciErrNum,"Failed to release program!");
}

void FC_FUNC_(kinetic1d_l,KINETIC1D_L)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, float *h, float*c, cl_mem *x, cl_mem *y, cl_mem *workx, cl_mem *worky,float *ekin){
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %lu, dimension dat: %lu, h: %f, c: %f, x: %p, workx: %p, y: %p, worky: %p\n",*command_queue, (long unsigned)*n, (long unsigned)*ndat, *h, *c, *x, *workx, *y, *worky);
#endif
    int FILTER_WIDTH = 30;
    if(*n<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=32, block_size_j=256/32;
    while (*n > block_size_i >= 1 && block_size_j > 4)
        { block_size_i *= 2; block_size_j /= 2;}
    cl_uint i = 0;
    clSetKernelArg(c_initialize_kernel_l, i++,sizeof(*n), (void*)n);
    clSetKernelArg(c_initialize_kernel_l, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(c_initialize_kernel_l, i++,sizeof(*x), (void*)x);
    clSetKernelArg(c_initialize_kernel_l, i++,sizeof(*worky), (void*)worky);
    clSetKernelArg(c_initialize_kernel_l, i++,sizeof(*c), (void*)c);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, c_initialize_kernel_l, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue c_initialize_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }
    i = 0;
    float scale = 0.5 / ( *h * (*h) );
    clSetKernelArg(kinetic1d_kernel_l, i++,sizeof(*n), (void*)n);
    clSetKernelArg(kinetic1d_kernel_l, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(kinetic1d_kernel_l, i++,sizeof(scale), (void*)&scale);
    clSetKernelArg(kinetic1d_kernel_l, i++,sizeof(*x), (void*)x);
    clSetKernelArg(kinetic1d_kernel_l, i++,sizeof(*workx), (void*)workx);
    clSetKernelArg(kinetic1d_kernel_l, i++,sizeof(*worky), (void*)worky);
    clSetKernelArg(kinetic1d_kernel_l, i++,sizeof(*y), (void*)y);
    clSetKernelArg(kinetic1d_kernel_l, i++,sizeof(float)*block_size_j*(block_size_i+FILTER_WIDTH), NULL);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kinetic1d_kernel_l, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue kinetic1d_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }
}
