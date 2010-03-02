#include "OpenCL_wrappers.h"
#include "Kinetic.h"
#include "Kinetic_k.h"

cl_kernel kinetic1d_kernel_d;
cl_kernel kinetic_k1d_kernel_d;

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
    kinetic1d_kernel_d=clCreateKernel(kinetic1dProgram,"kinetic1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = clReleaseProgram(kinetic1dProgram);
    oclErrorCheck(ciErrNum,"Failed to release program!");

    ciErrNum = CL_SUCCESS; 
    cl_program kinetic_k1dProgram = clCreateProgramWithSource(*context,1,(const char**) &kinetic_k1d_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(kinetic_k1dProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build kinetic_k1d program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(kinetic_k1dProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    kinetic_k1d_kernel_d=clCreateKernel(kinetic_k1dProgram,"kinetic_k1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = clReleaseProgram(kinetic_k1dProgram);
    oclErrorCheck(ciErrNum,"Failed to release program!");
}

void FC_FUNC_(kinetic_k_d,KINETIC_K_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3, double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y, double * c_in,  double *k1, double *k2, double *k3) {
  cl_int ciErrNum;
  double c = *c_in  + .5 * (*k1 * *k1 + *k2 * *k2 + *k3 * *k3);
  cl_uint ng = *n1 * *n2 * *n3 * 2;
  cl_uint i = 0;
  size_t block_size_i1=64;

  clSetKernelArg(c_initialize_kernel_d, i++,sizeof(ng), (void*)&ng);
  clSetKernelArg(c_initialize_kernel_d, i++,sizeof(*x), (void*)x);
  clSetKernelArg(c_initialize_kernel_d, i++,sizeof(*work_y), (void*)work_y);
  clSetKernelArg(c_initialize_kernel_d, i++,sizeof(c), (void*)&c);
  size_t localWorkSize1[] = { block_size_i1 };
  size_t globalWorkSize1[] ={ shrRoundUp(block_size_i1,ng) };
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, c_initialize_kernel_d, 1, NULL, globalWorkSize1, localWorkSize1, 0, NULL, NULL);
  if (ciErrNum != CL_SUCCESS)
  {
    fprintf(stderr,"Error %d: Failed to enqueue c_initialize_d kernel!\n",ciErrNum);
    exit(1);
  }
  double un = 1.0;
  i = 0;
  clSetKernelArg(c_initialize_kernel_d, i++,sizeof(ng), (void*)&ng);
  clSetKernelArg(c_initialize_kernel_d, i++,sizeof(*x), (void*)x);
  clSetKernelArg(c_initialize_kernel_d, i++,sizeof(*work_x), (void*)work_x);
  clSetKernelArg(c_initialize_kernel_d, i++,sizeof(un), (void*)&un);
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, c_initialize_kernel_d, 1, NULL, globalWorkSize1, localWorkSize1, 0, NULL, NULL);
  if (ciErrNum != CL_SUCCESS)
  {
    fprintf(stderr,"Error %d: Failed to enqueue c_initialize_d kernel!\n",ciErrNum);
    exit(1);
  }
  int FILTER_WIDTH = 32;
  size_t block_size_i2=32;
  size_t block_size_j2=8;
  i = 0;
  double scale_1 = -0.5 / ( h[0] * h[0] );
  double scale_2 = *k1 / h[0];
//  printf("%lf %lf\n",scale_1,scale_2);
  ng = *n2 * *n3;
  size_t localWorkSize2[] = { block_size_i2, block_size_j2 };
  size_t globalWorkSize2[] ={ shrRoundUp(block_size_i2,*n1), shrRoundUp(block_size_j2,ng) };
//__kernel void kinetic_k1dKernel_d(size_t n, size_t ndat, double scale_1, double scale_2, __global const double * x_in, __global double * x, __global const double * y_in, __global double * y, __local double * tmp, __local double * tmp_y )
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(*n1), (void*)n1);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(ng), (void*)&ng);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(scale_1), (void*)&scale_1);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(scale_2), (void*)&scale_2);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(*work_x), (void*)work_x);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(*x), (void*)x);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(*work_y), (void*)work_y);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(*y), (void*)y);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(double)*block_size_j2*2*(block_size_i2+FILTER_WIDTH+1), NULL);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(double)*block_size_j2*2*(block_size_i2+1), NULL);
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kinetic_k1d_kernel_d, 2, NULL, globalWorkSize2, localWorkSize2, 0, NULL, NULL);
  if (ciErrNum != CL_SUCCESS)
  {
      fprintf(stderr,"Error %d: Failed to enqueue kinetic_k1d_d 1 kernel!\n",ciErrNum);
      fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize2[0],(long unsigned)globalWorkSize2[1]);
      fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize2[0],(long unsigned)localWorkSize2[1]);
      exit(1);
  }
  i = 0;
  scale_1 = -0.5 / ( h[1] * h[1] );
  scale_2 = *k2 / h[1];
  ng = *n1 * *n3;
  localWorkSize2[0] =  block_size_i2; localWorkSize2[1] = block_size_j2 ;
  globalWorkSize2[0] = shrRoundUp(block_size_i2,*n2); globalWorkSize2[1] = shrRoundUp(block_size_j2,ng) ;
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(*n2), (void*)n2);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(ng), (void*)&ng);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(scale_1), (void*)&scale_1);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(scale_2), (void*)&scale_2);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(*x), (void*)x);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(*work_x), (void*)work_x);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(*y), (void*)y);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(*work_y), (void*)work_y);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(double)*block_size_j2*2*(block_size_i2+FILTER_WIDTH+1), NULL);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(double)*block_size_j2*2*(block_size_i2+1), NULL);
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kinetic_k1d_kernel_d, 2, NULL, globalWorkSize2, localWorkSize2, 0, NULL, NULL);
  if (ciErrNum != CL_SUCCESS)
  {
      fprintf(stderr,"Error %d: Failed to enqueue kinetic_k1d_d 2 kernel!\n",ciErrNum);
      fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize2[0],(long unsigned)globalWorkSize2[1]);
      fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize2[0],(long unsigned)localWorkSize2[1]);
      exit(1);
  }
  i = 0;
  scale_1 = -0.5 / ( h[2] * h[2] );
  scale_2 = *k3 / h[2];
  ng = *n1 * *n2;
  localWorkSize2[0] =  block_size_i2; localWorkSize2[1] = block_size_j2 ;
  globalWorkSize2[0] = shrRoundUp(block_size_i2,*n3); globalWorkSize2[1] = shrRoundUp(block_size_j2,ng) ;
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(*n3), (void*)n3);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(ng), (void*)&ng);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(scale_1), (void*)&scale_1);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(scale_2), (void*)&scale_2);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(*work_x), (void*)work_x);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(*x), (void*)x);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(*work_y), (void*)work_y);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(*y), (void*)y);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(double)*block_size_j2*2*(block_size_i2+FILTER_WIDTH+1), NULL);
  clSetKernelArg(kinetic_k1d_kernel_d, i++,sizeof(double)*block_size_j2*2*(block_size_i2+1), NULL);
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kinetic_k1d_kernel_d, 2, NULL, globalWorkSize2, localWorkSize2, 0, NULL, NULL);
  if (ciErrNum != CL_SUCCESS)
  {
      fprintf(stderr,"Error %d: Failed to enqueue kinetic_k1d_d 3 kernel!\n",ciErrNum);
      fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize2[0],(long unsigned)globalWorkSize2[1]);
      fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize2[0],(long unsigned)localWorkSize2[1]);
      exit(1);
  }

}

void FC_FUNC_(kinetic1d_d,KINETIC1D_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, double *h, double*c, cl_mem *x, cl_mem *y, cl_mem *workx, cl_mem *worky,double *ekin){
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %lu, dimension dat: %lu, h: %f, c: %f, x: %p, workx: %p, y: %p, worky: %p\n",*command_queue, (long unsigned)*n, (long unsigned)*ndat, *h, *c, *x, *workx, *y, *worky);
#endif
    int FILTER_WIDTH = 30;
    if(*n<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i1=64;
    cl_uint ng = *n * *ndat;
    cl_uint i = 0;
    clSetKernelArg(c_initialize_kernel_d, i++,sizeof(ng), (void*)&ng);
    clSetKernelArg(c_initialize_kernel_d, i++,sizeof(*x), (void*)x);
    clSetKernelArg(c_initialize_kernel_d, i++,sizeof(*worky), (void*)worky);
    clSetKernelArg(c_initialize_kernel_d, i++,sizeof(*c), (void*)c);
    size_t localWorkSize1[] = { block_size_i1 };
    size_t globalWorkSize1[] ={ shrRoundUp(block_size_i1,ng) };
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, c_initialize_kernel_d, 1, NULL, globalWorkSize1, localWorkSize1, 0, NULL, NULL);
    
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue c_initialize_d kernel!\n",ciErrNum);
        exit(1);
    }
    i = 0;
    double scale = 0.5 / ( *h * (*h) );
    size_t block_size_i2=32;
    size_t block_size_j2=16;
    size_t localWorkSize2[] = { block_size_i2, block_size_j2 };
    size_t globalWorkSize2[] ={ shrRoundUp(block_size_i2,*n), shrRoundUp(block_size_j2,*ndat) };

    clSetKernelArg(kinetic1d_kernel_d, i++,sizeof(*n), (void*)n);
    clSetKernelArg(kinetic1d_kernel_d, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(kinetic1d_kernel_d, i++,sizeof(scale), (void*)&scale);
    clSetKernelArg(kinetic1d_kernel_d, i++,sizeof(*x), (void*)x);
    clSetKernelArg(kinetic1d_kernel_d, i++,sizeof(*workx), (void*)workx);
    clSetKernelArg(kinetic1d_kernel_d, i++,sizeof(*worky), (void*)worky);
    clSetKernelArg(kinetic1d_kernel_d, i++,sizeof(*y), (void*)y);
    clSetKernelArg(kinetic1d_kernel_d, i++,sizeof(double)*block_size_j2*(block_size_i2+FILTER_WIDTH+3), NULL);
    clSetKernelArg(kinetic1d_kernel_d, i++,sizeof(double)*block_size_j2*(block_size_i2+1), NULL);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kinetic1d_kernel_d, 2, NULL, globalWorkSize2, localWorkSize2, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue kinetic1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize2[0],(long unsigned)globalWorkSize2[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize2[0],(long unsigned)localWorkSize2[1]);
        exit(1);
    }
}


void clean_kinetic_kernels(){
  clReleaseKernel(kinetic1d_kernel_d);
  clReleaseKernel(kinetic_k1d_kernel_d);
}
