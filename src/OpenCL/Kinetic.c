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
  double c = *c_in  + .5 * (*k1 * *k1 + *k2 * *k2 + *k3 * *k3);
  cl_uint ng = *n1 * *n2 * *n3 * 2;
  c_initialize_generic(c_initialize_kernel_d, command_queue, &ng, x, work_y, &c);  
  c = 1.0;
  c_initialize_generic(c_initialize_kernel_d, command_queue, &ng, x, work_x, &c);  
  double scale_1 = -0.5 / ( h[2] * h[2] );
  double scale_2 = *k3 / h[2];
  ng = *n2 * *n1;
  kinetic_k_generic(kinetic_k1d_kernel_d, command_queue, n3, &ng, &scale_1, &scale_2, work_x, x, work_y, y);
  scale_1 = -0.5 / ( h[1] * h[1] );
  scale_2 = *k2 / h[1];
  ng = *n1 * *n3;
  kinetic_k_generic(kinetic_k1d_kernel_d, command_queue, n2, &ng, &scale_1, &scale_2, x, work_x, y, work_y);
  scale_1 = -0.5 / ( h[0] * h[0] );
  scale_2 = *k1 / h[0];
  ng = *n3 * *n2;
  kinetic_k_generic(kinetic_k1d_kernel_d, command_queue, n1, &ng, &scale_1, &scale_2, work_x, x, work_y, y);
}

void FC_FUNC_(kinetic_stable_d,KINETIC_STABLE_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3, double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y, cl_mem *tmp_x, cl_mem *tmp_y) {
  double scale = -0.5 / ( h[2] * h[2] );
  cl_uint ng = *n2 * *n1 * 4;
  cl_uint n = 2 * *n3;
  kinetic_generic(kinetic1d_kernel_d, command_queue, &n, &ng, &scale, x, work_x, y, work_y);
  scale = -0.5 / ( h[1] * h[1] );
  n = 2 * *n2;
  ng = *n1 * *n3 * 4;
  kinetic_generic(kinetic1d_kernel_d, command_queue, &n, &ng, &scale, work_x, tmp_x, work_y, tmp_y);
  scale = -0.5 / ( h[0] * h[0] );
  n = 2 * *n1;
  ng = *n3 * *n2 * 4;
  kinetic_generic(kinetic1d_kernel_d, command_queue, &n, &ng, &scale, tmp_x, work_x, tmp_y, work_y);
}

void FC_FUNC_(kinetic_d,KINETIC_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3, double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y) {
  double scale = -0.5 / ( h[2] * h[2] );
  cl_uint ng = *n2 * *n1 * 4;
  cl_uint n = 2 * *n3;
  kinetic_generic(kinetic1d_kernel_d, command_queue, &n, &ng, &scale, x, work_x, y, work_y);
  scale = -0.5 / ( h[1] * h[1] );
  n = 2 * *n2;
  ng = *n1 * *n3 * 4;
  kinetic_generic(kinetic1d_kernel_d, command_queue, &n, &ng, &scale, work_x, x, work_y, y);
  scale = -0.5 / ( h[0] * h[0] );
  n = 2 * *n1;
  ng = *n3 * *n2 * 4;
  kinetic_generic(kinetic1d_kernel_d, command_queue, &n, &ng, &scale, x, work_x, y, work_y);
}

void FC_FUNC_(kinetic1d_d,KINETIC1D_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, double *h, double*c, cl_mem *x, cl_mem *y, cl_mem *workx, cl_mem *worky,double *ekin){
  cl_uint ng = *n * *ndat;
  c_initialize_generic(c_initialize_kernel_d, command_queue, &ng, x, worky, c);
  double scale = - 0.5 / ( *h * *h );
  kinetic_generic(kinetic1d_kernel_d, command_queue, n, ndat, &scale, x, workx, worky, y);
}

void clean_kinetic_kernels(){
  clReleaseKernel(kinetic1d_kernel_d);
  clReleaseKernel(kinetic_k1d_kernel_d);
}
