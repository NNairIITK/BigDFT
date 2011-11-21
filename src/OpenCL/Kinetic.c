//! @file
//!  Kinetic kernel in OpenCL
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include "Initialize.h"
#include "Kinetic.h"
#include "Kinetic_k.h"
#include "OpenCL_wrappers.h"
#include "Kinetic_Generator.h"
#include "Kinetic_k_Generator.h"

inline void kinetic_generic(cl_kernel kernel, bigdft_command_queue command_queue, cl_uint *n, cl_uint *ndat, double *scale, cl_mem *x_in, cl_mem *x_out, cl_mem *y_in, cl_mem *y_out) {
  int FILTER_WIDTH = 32;
  cl_int ciErrNum;
  assert(*n>=FILTER_WIDTH);
  size_t block_size_i=FILTER_WIDTH;
  size_t block_size_j=command_queue->device_infos.MAX_WORK_GROUP_SIZE/FILTER_WIDTH;
  if(block_size_j>16)
    block_size_j=16;
  size_t localWorkSize[] = { block_size_i, block_size_j };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat) };
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
  clSetKernelArg(kernel, i++,sizeof(*scale), (void*)scale);
  clSetKernelArg(kernel, i++,sizeof(*x_in), (void*)x_in);
  clSetKernelArg(kernel, i++,sizeof(*x_out), (void*)x_out);
  clSetKernelArg(kernel, i++,sizeof(*y_in), (void*)y_in);
  clSetKernelArg(kernel, i++,sizeof(*y_out), (void*)y_out);
  clSetKernelArg(kernel, i++,sizeof(double)*block_size_j*(block_size_i+FILTER_WIDTH+1), NULL);
  clSetKernelArg(kernel, i++,sizeof(double)*block_size_j*(block_size_i+1), NULL);
  ciErrNum = clEnqueueNDRangeKernel  (command_queue->command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue kinetic kernel!");
}


inline void kinetic_k_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n, cl_uint *ndat, double *scale1, double *scale2, cl_mem *x_in, cl_mem *x_out, cl_mem *y_in, cl_mem *y_out) {
  int FILTER_WIDTH = 32;
  cl_int ciErrNum;
  assert(*n>=FILTER_WIDTH);
  size_t block_size_i=FILTER_WIDTH;
  size_t block_size_j=8;
  size_t localWorkSize[] = { block_size_i, block_size_j };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat) };
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
  clSetKernelArg(kernel, i++,sizeof(*scale1), (void*)scale1);
  clSetKernelArg(kernel, i++,sizeof(*scale2), (void*)scale2);
  clSetKernelArg(kernel, i++,sizeof(*x_in), (void*)x_in);
  clSetKernelArg(kernel, i++,sizeof(*x_out), (void*)x_out);
  clSetKernelArg(kernel, i++,sizeof(*y_in), (void*)y_in);
  clSetKernelArg(kernel, i++,sizeof(*y_out), (void*)y_out);
  ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue kinetic_k kernel!");
} 

inline void kinetic_k_generic_2(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n, cl_uint *ndat, double *scale1, double *scale2, cl_mem *x_in_r, cl_mem *x_in_i, cl_mem *x_out_r, cl_mem *x_out_i, cl_mem *y_in_r, cl_mem *y_in_i, cl_mem *y_out_r, cl_mem *y_out_i) {
  int FILTER_WIDTH = 32;
  cl_int ciErrNum;
  assert(*n>=FILTER_WIDTH);
  size_t block_size_i=FILTER_WIDTH;
  size_t block_size_j=8;
  size_t localWorkSize[] = { block_size_i, block_size_j };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat) };
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
  clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
  clSetKernelArg(kernel, i++,sizeof(*scale1), (void*)scale1);
  clSetKernelArg(kernel, i++,sizeof(*scale2), (void*)scale2);
  clSetKernelArg(kernel, i++,sizeof(*x_in_r), (void*)x_in_r);
  clSetKernelArg(kernel, i++,sizeof(*x_in_i), (void*)x_in_i);
  clSetKernelArg(kernel, i++,sizeof(*x_out_r), (void*)x_out_r);
  clSetKernelArg(kernel, i++,sizeof(*x_out_i), (void*)x_out_i);
  clSetKernelArg(kernel, i++,sizeof(*y_in_r), (void*)y_in_r);
  clSetKernelArg(kernel, i++,sizeof(*y_in_i), (void*)y_in_i);
  clSetKernelArg(kernel, i++,sizeof(*y_out_r), (void*)y_out_r);
  clSetKernelArg(kernel, i++,sizeof(*y_out_i), (void*)y_out_i);
  ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue kinetic_k_2 kernel!");
} 



cl_program kineticProgram;
cl_program kinetic_kProgram;

void create_kinetic_kernels(struct bigdft_kernels * kernels) {
    cl_int ciErrNum=CL_SUCCESS;
    kernels->kinetic1d_kernel_d=clCreateKernel(kineticProgram,"kinetic1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kinetic1dKernel_d kernel!");
    kernels->kinetic1d_f_kernel_d=clCreateKernel(kineticProgram,"kinetic1d_fKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kinetic1d_fKernel_d kernel!");
    kernels->kinetic_k1d_kernel_d=clCreateKernel(kinetic_kProgram,"kinetic_k1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kinetic_k1dKernel_d kernel!");
    kernels->kinetic_k1d_kernel_d_2=clCreateKernel(kinetic_kProgram,"kinetic_k1dKernel_d_2",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kinetic_k1dKernel_d_2 kernel!");
    kernels->kinetic_k1d_f_kernel_d_2=clCreateKernel(kinetic_kProgram,"kinetic_k1d_fKernel_d_2",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kinetic_k1d_fKernel_d_2 kernel!");
}

/*static void print_program_binary(cl_program prog){
    cl_int ciErrNum;
    cl_uint dev_num;
    ciErrNum = clGetProgramInfo(kinetic_kProgram, CL_PROGRAM_NUM_DEVICES, sizeof(cl_uint), &dev_num, NULL);
    oclErrorCheck(ciErrNum,"Failed to get num devices!");
    size_t *prog_sizes=(size_t *)malloc(dev_num*sizeof(cl_uint));
    ciErrNum = clGetProgramInfo(kinetic_kProgram, CL_PROGRAM_BINARY_SIZES , sizeof(size_t)*dev_num, prog_sizes, NULL);
    oclErrorCheck(ciErrNum,"Failed to get program binary sizes!");
    unsigned char ** binaries = (unsigned char **)malloc(dev_num * sizeof(unsigned char *));
    unsigned int i;
    size_t prog_sum_size = 0;
    for(i=0; i<dev_num; i++){
      binaries[i] = (unsigned char *)malloc(prog_sizes[i]*sizeof(unsigned char));
      prog_sum_size += prog_sizes[i];
    }
    ciErrNum = clGetProgramInfo(kinetic_kProgram, CL_PROGRAM_BINARIES , prog_sum_size, binaries, NULL);
    oclErrorCheck(ciErrNum,"Failed to get program binaries!");
    for(i=0; i<dev_num; i++){
      printf("%.*s\n",(int)prog_sizes[i],binaries[i]);
      free(binaries[i]);
    }
    free(binaries);
}*/

void build_kinetic_programs(cl_context * context){
    struct bigdft_device_infos infos;
    get_context_devices_infos(context, &infos);
    cl_int ciErrNum=CL_SUCCESS;
    char * code = generate_kinetic_program(&infos);
    kineticProgram = clCreateProgramWithSource(*context,1,(const char**) &code, NULL, &ciErrNum);
    free(code);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(kineticProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build kinetic program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(kineticProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    code = generate_kinetic_k_program(&infos);
    kinetic_kProgram = clCreateProgramWithSource(*context,1,(const char**) &code, NULL, &ciErrNum);
    free(code);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(kinetic_kProgram, 0, NULL, "", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build kinetic_k program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(kinetic_kProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
//    print_program_binary(kinetic_kProgram);
}

void FC_FUNC_(kinetic_k_d,KINETIC_K_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y, double * c_in,  double *k) {
  double c = *c_in  + .5 * (k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
  cl_uint n1 = dimensions[0];
  cl_uint n2 = dimensions[1];
  cl_uint n3 = dimensions[2];
  cl_uint ng = n1 * n2 * n3 * 2;
  c_initialize_generic((*command_queue)->kernels.c_initialize_kernel_d, (*command_queue)->command_queue, &ng, x, work_y, &c);  
  c = 1.0;
  c_initialize_generic((*command_queue)->kernels.c_initialize_kernel_d, (*command_queue)->command_queue, &ng, x, work_x, &c);  
  double scale_1 = -0.5 / ( h[2] * h[2] );
  double scale_2 = k[2] / h[2];
  ng = n2 * n1;
  kinetic_k_generic((*command_queue)->kernels.kinetic_k1d_kernel_d, (*command_queue)->command_queue, &n3, &ng, &scale_1, &scale_2, work_x, x, work_y, y);
  scale_1 = -0.5 / ( h[1] * h[1] );
  scale_2 = k[1] / h[1];
  ng = n1 * n3;
  kinetic_k_generic((*command_queue)->kernels.kinetic_k1d_kernel_d, (*command_queue)->command_queue, &n2, &ng, &scale_1, &scale_2, x, work_x, y, work_y);
  scale_1 = -0.5 / ( h[0] * h[0] );
  scale_2 = k[0] / h[0];
  ng = n3 * n2;
  kinetic_k_generic((*command_queue)->kernels.kinetic_k1d_kernel_d, (*command_queue)->command_queue, &n1, &ng, &scale_1, &scale_2, work_x, x, work_y, y);
}

void FC_FUNC_(kinetic_stable_d,KINETIC_STABLE_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y, cl_mem *tmp_x, cl_mem *tmp_y) {
  cl_uint n1 = dimensions[0] * 2;
  cl_uint n2 = dimensions[1] * 2;
  cl_uint n3 = dimensions[2] * 2;
  double scale = -0.5 / ( h[2] * h[2] );
  cl_uint ng = n2 * n1;
  kinetic_generic((*command_queue)->kernels.kinetic1d_kernel_d, (*command_queue), &n3, &ng, &scale, x, work_x, y, work_y);
  scale = -0.5 / ( h[1] * h[1] );
  ng = n1 * n3;
  kinetic_generic((*command_queue)->kernels.kinetic1d_kernel_d, (*command_queue), &n2, &ng, &scale, work_x, tmp_x, work_y, tmp_y);
  scale = -0.5 / ( h[0] * h[0] );
  ng = n3 * n2;
  kinetic_generic((*command_queue)->kernels.kinetic1d_kernel_d, (*command_queue), &n1, &ng, &scale, tmp_x, work_x, tmp_y, work_y);
}

void FC_FUNC_(kinetic_k_d_generic,KINETIC_K_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, double *h, double *k, cl_mem *x_r, cl_mem *x_i, cl_mem *y_r, cl_mem *y_i, cl_mem *work_x_r, cl_mem *work_x_i, cl_mem *work_y_r, cl_mem *work_y_i){
  double scale_1, scale_2;
  cl_uint ng;
  cl_uint n1 = dimensions[0] * 2;             
  cl_uint n2 = dimensions[1] * 2;
  cl_uint n3 = dimensions[2] * 2;
  if( !periodic[0] ) n1 += 2*7;
  if( !periodic[1] ) n2 += 2*7;
  if( !periodic[2] ) n3 += 2*7;
  scale_1 = -0.5 / ( h[2] * h[2] );
  scale_2 = k[2] / h[2];
  ng = n1 * n2;
  if( periodic[2] ) {
    kinetic_k_generic_2((*command_queue)->kernels.kinetic_k1d_kernel_d_2, (*command_queue)->command_queue, &n3, &ng, &scale_1, &scale_2, x_r, x_i, work_x_r, work_x_i, y_r, y_i, work_y_r, work_y_i);
  } else {
    kinetic_k_generic_2((*command_queue)->kernels.kinetic_k1d_f_kernel_d_2, (*command_queue)->command_queue, &n3, &ng, &scale_1, &scale_2, x_r, x_i, work_x_r, work_x_i, y_r, y_i, work_y_r, work_y_i);
  }
  scale_1 = -0.5 / ( h[1] * h[1] );
  scale_2 = k[1] / h[1];
  ng = n1 * n3;
  if( periodic[1] ) {
    kinetic_k_generic_2((*command_queue)->kernels.kinetic_k1d_kernel_d_2, (*command_queue)->command_queue, &n2, &ng, &scale_1, &scale_2, work_x_r, work_x_i, x_r, x_i, work_y_r, work_y_i, y_r, y_i);
  } else {
    kinetic_k_generic_2((*command_queue)->kernels.kinetic_k1d_f_kernel_d_2, (*command_queue)->command_queue, &n2, &ng, &scale_1, &scale_2, work_x_r, work_x_i, x_r, x_i, work_y_r, work_y_i, y_r, y_i);
  }
  scale_1 = -0.5 / ( h[0] * h[0] );
  scale_2 = k[0] / h[0];
  ng = n2 * n3;
  if( periodic[0] ) {
    kinetic_k_generic_2((*command_queue)->kernels.kinetic_k1d_kernel_d_2, (*command_queue)->command_queue, &n1, &ng, &scale_1, &scale_2, x_r, x_i, work_x_r, work_x_i, y_r, y_i, work_y_r, work_y_i);
  } else {
    kinetic_k_generic_2((*command_queue)->kernels.kinetic_k1d_f_kernel_d_2, (*command_queue)->command_queue, &n1, &ng, &scale_1, &scale_2, x_r, x_i, work_x_r, work_x_i, y_r, y_i, work_y_r, work_y_i);
  }
}
 

void FC_FUNC_(kinetic_d_generic,KINETIC_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y) {
  double scale;
  cl_uint ng;
  cl_uint n1 = dimensions[0] * 2;             
  cl_uint n2 = dimensions[1] * 2;
  cl_uint n3 = dimensions[2] * 2;
  if( !periodic[0] ) n1 += 2*7;
  if( !periodic[1] ) n2 += 2*7;
  if( !periodic[2] ) n3 += 2*7;
  scale = -0.5 / ( h[2] * h[2] );
  ng = n1 * n2;
  if( periodic[2] ) {
    kinetic_generic((*command_queue)->kernels.kinetic1d_kernel_d, (*command_queue), &n3, &ng, &scale, x, work_x, y, work_y);
  } else {
    kinetic_generic((*command_queue)->kernels.kinetic1d_f_kernel_d, (*command_queue), &n3, &ng, &scale, x, work_x, y, work_y);
  }
  scale = -0.5 / ( h[1] * h[1] );
  ng = n1 * n3;
  if( periodic[1] ) {
    kinetic_generic((*command_queue)->kernels.kinetic1d_kernel_d, (*command_queue), &n2, &ng, &scale, work_x, x, work_y, y);
  } else {
    kinetic_generic((*command_queue)->kernels.kinetic1d_f_kernel_d, (*command_queue), &n2, &ng, &scale, work_x, x, work_y, y);
  }
  scale = -0.5 / ( h[0] * h[0] );
  ng = n2 * n3;
  if( periodic[0] ) {
    kinetic_generic((*command_queue)->kernels.kinetic1d_kernel_d, (*command_queue), &n1, &ng, &scale, x, work_x, y, work_y);
  } else {
    kinetic_generic((*command_queue)->kernels.kinetic1d_f_kernel_d, (*command_queue), &n1, &ng, &scale, x, work_x, y, work_y);
  }
}

void FC_FUNC_(kinetic_d,KINETIC_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y) {
  cl_uint n1 = dimensions[0] * 2;             
  cl_uint n2 = dimensions[1] * 2;
  cl_uint n3 = dimensions[2] * 2;
  double scale = -0.5 / ( h[2] * h[2] );
  cl_uint ng = n2 * n1;
  kinetic_generic((*command_queue)->kernels.kinetic1d_kernel_d, (*command_queue), &n3, &ng, &scale, x, work_x, y, work_y);
  scale = -0.5 / ( h[1] * h[1] );
  ng = n1 * n3;
  kinetic_generic((*command_queue)->kernels.kinetic1d_kernel_d, (*command_queue), &n2, &ng, &scale, work_x, x, work_y, y);
  scale = -0.5 / ( h[0] * h[0] );
  ng = n3 * n2;
  kinetic_generic((*command_queue)->kernels.kinetic1d_kernel_d, (*command_queue), &n1, &ng, &scale, x, work_x, y, work_y);
}

void FC_FUNC_(kinetic1d_d,KINETIC1D_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_uint *ndat, double *h, double*c, cl_mem *x, cl_mem *y, cl_mem *workx, cl_mem *worky,double *ekin){
  cl_uint ng = *n * *ndat;
  c_initialize_generic((*command_queue)->kernels.c_initialize_kernel_d, (*command_queue)->command_queue, &ng, x, worky, c);
  double scale = - 0.5 / ( *h * *h );
  kinetic_generic((*command_queue)->kernels.kinetic1d_kernel_d, (*command_queue), n, ndat, &scale, x, workx, worky, y);
  *ekin = 0.0;
}

void clean_kinetic_kernels(struct bigdft_kernels * kernels){
  cl_int ciErrNum;
  ciErrNum = clReleaseKernel(kernels->kinetic1d_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->kinetic1d_f_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->kinetic_k1d_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->kinetic_k1d_kernel_d_2);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
  ciErrNum = clReleaseKernel(kernels->kinetic_k1d_f_kernel_d_2);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
}

void clean_kinetic_programs(){
  cl_int ciErrNum;
  ciErrNum = clReleaseProgram(kineticProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
  ciErrNum = clReleaseProgram(kinetic_kProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
}
