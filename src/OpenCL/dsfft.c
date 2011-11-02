#include "OpenCL_wrappers.h"
#include "fft_generator.h"

cl_uint use_constant_memory=1;

inline void fft_generated_generic(cl_kernel kernel, bigdft_command_queue command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out,cl_mem *cosi){
    cl_int ciErrNum;
    cl_uint shared_size_used=0;
    if(*n <= 64)
       shared_size_used=512;
    else if(*n <= 256)
       shared_size_used=1024;
    int FFT_LENGTH = (*n / 16) * 16 + (*n % 16 ? 16 : 0);
    int LINE_LIMIT = ( (shared_size_used/FFT_LENGTH) & (~3) ) <= command_queue->device_infos.MAX_WORK_GROUP_SIZE/FFT_LENGTH ? (shared_size_used/FFT_LENGTH) & (~3) : command_queue->device_infos.MAX_WORK_GROUP_SIZE/FFT_LENGTH ;
    int LINE_NUMBER = 1;
    while( LINE_LIMIT /= 2 )
      LINE_NUMBER *= 2;
    if( LINE_NUMBER > FFT_LENGTH )
       LINE_NUMBER = FFT_LENGTH;
    if( LINE_NUMBER < 1 )
       LINE_NUMBER=1;
    size_t block_size_i=FFT_LENGTH, block_size_j=LINE_NUMBER;
    cl_uint i = 0;
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    if(!use_constant_memory){
      ciErrNum = clSetKernelArg(kernel, i++,sizeof(*cosi), (void*)cosi);
    }
    size_t localWorkSize[] = { block_size_i,block_size_j };
//    printf("%lu %lu\n",block_size_i,block_size_j);
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel(command_queue->command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue fft kernel!");
}

inline void fft_k_generated_generic(cl_kernel kernel, bigdft_command_queue command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out,cl_mem *k,cl_mem *cosi){
    cl_int ciErrNum;
    cl_uint shared_size_used=0;
    if(*n <= 64)
       shared_size_used=512;
    else if(*n <= 256)
       shared_size_used=1024;
    int FFT_LENGTH = (*n / 16) * 16 + (*n % 16 ? 16 : 0);
    int LINE_LIMIT = ( (shared_size_used/FFT_LENGTH) & (~3) ) <= command_queue->device_infos.MAX_WORK_GROUP_SIZE/FFT_LENGTH ? (shared_size_used/FFT_LENGTH) & (~3) : command_queue->device_infos.MAX_WORK_GROUP_SIZE/FFT_LENGTH ;
    int LINE_NUMBER = 1;
    while( LINE_LIMIT /= 2 )
      LINE_NUMBER *= 2;
    if( LINE_NUMBER > FFT_LENGTH )
       LINE_NUMBER = FFT_LENGTH;
    if( LINE_NUMBER < 1 )
       LINE_NUMBER=1;
    size_t block_size_i=FFT_LENGTH, block_size_j=LINE_NUMBER;
    cl_uint i = 0;
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*k), (void*)k);
    if(!use_constant_memory){
      ciErrNum = clSetKernelArg(kernel, i++,sizeof(*cosi), (void*)cosi);
    }
    size_t localWorkSize[] = { block_size_i,block_size_j };
//    printf("%lu %lu\n",block_size_i,block_size_j);
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel(command_queue->command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue fft kernel!");
}


void FC_FUNC_(customize_fft,CUSTOMIZE_FFT)(cl_uint *dimensions) {
  assert(dimensions[0]<=256);
  assert(dimensions[1]<=256);
  assert(dimensions[2]<=256);
  fft_size[0] = dimensions[0];
  fft_size[1] = dimensions[1];
  fft_size[2] = dimensions[2];
}

cl_mem cossind0;
cl_mem cossind1;
cl_mem cossind2;

void FC_FUNC_(fft1d_d,FFT1D_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    if(fft_size[0] == *n)
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d0_d, *command_queue, n, ndat, psi, out, &cossind0);
    else if(fft_size[1] == *n)
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d1_d, *command_queue, n, ndat, psi, out, &cossind1);
    else if(fft_size[2] == *n)
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d2_d, *command_queue, n, ndat, psi, out, &cossind2);
}

void FC_FUNC_(fft1d_r_d,FFT1D_R_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    if(fft_size[0] == *n)
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d0_r_d, *command_queue, n, ndat, psi, out, &cossind0);
    else if(fft_size[1] == *n)
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d1_r_d, *command_queue, n, ndat, psi, out, &cossind1);
    else if(fft_size[2] == *n)
      fft_generated_generic((*command_queue)->kernels.fft_kernel_d2_r_d, *command_queue, n, ndat, psi, out, &cossind2);
}

void FC_FUNC_(fft3d_d,FFT1D_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,cl_mem *psi,cl_mem *out,cl_mem *tmp){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    assert(fft_size[2]==n3);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d2_d, *command_queue, &n3, &ndat, psi, out, &cossind2);
    ndat = n1 * n3;
    assert(fft_size[1]==n2);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d1_d, *command_queue, &n2, &ndat, out, tmp, &cossind1);
    ndat = n2 * n3;
    assert(fft_size[0]==n1);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d0_d, *command_queue, &n1, &ndat, tmp, out, &cossind0);
}

void FC_FUNC_(fft3d_k_d,FFT1D_K_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,cl_mem *psi,cl_mem *out,cl_mem *tmp,cl_mem *k){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    assert(fft_size[2]==n3);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d2_d, *command_queue, &n3, &ndat, psi, out, &cossind2);
    ndat = n1 * n3;
    assert(fft_size[1]==n2);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d1_d, *command_queue, &n2, &ndat, out, tmp, &cossind1);
    ndat = n2 * n3;
    assert(fft_size[0]==n1);
    fft_k_generated_generic((*command_queue)->kernels.fft_kernel_k_d0_d, *command_queue, &n1, &ndat, tmp, out, k, &cossind0);
}


void FC_FUNC_(fft3d_r_d,FFT1D_R_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,cl_mem *psi,cl_mem *out,cl_mem *tmp){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    assert(fft_size[2]==n3);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d2_r_d, *command_queue, &n3, &ndat, psi, out, &cossind2);
    ndat = n1 * n3;
    assert(fft_size[1]==n2);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d1_r_d, *command_queue, &n2, &ndat, out, tmp, &cossind1);
    ndat = n2 * n3;
    assert(fft_size[0]==n1);
    fft_generated_generic((*command_queue)->kernels.fft_kernel_d0_r_d, *command_queue, &n1, &ndat, tmp, out, &cossind0);
}

cl_program fftProgramd0;
cl_program fftProgramd1;
cl_program fftProgramd2;

cl_uint fft_size[3]={0,0,0};

void create_fft_kernels(struct bigdft_kernels * kernels){
  cl_int ciErrNum = CL_SUCCESS;
  char kernel_name[256];
  if(fft_size[0]!=0){
    sprintf(kernel_name,"fftKernel_%u_d",fft_size[0]);
    kernels->fft_kernel_d0_d=clCreateKernel(fftProgramd0,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d0_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r_d",fft_size[0]);
    kernels->fft_kernel_d0_r_d=clCreateKernel(fftProgramd0,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d0_r_d kernel!");
    sprintf(kernel_name,"fftKernel_k_%u_d",fft_size[0]);
    kernels->fft_kernel_k_d0_d=clCreateKernel(fftProgramd0,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_k_d0_d kernel!");
  }
  if(fft_size[1]!=0){
    sprintf(kernel_name,"fftKernel_%u_d",fft_size[1]);
    kernels->fft_kernel_d1_d=clCreateKernel(fftProgramd1,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d1_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r_d",fft_size[1]);
    kernels->fft_kernel_d1_r_d=clCreateKernel(fftProgramd1,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d1_r_d kernel!");
    sprintf(kernel_name,"fftKernel_k_%u_d",fft_size[1]);
    kernels->fft_kernel_k_d1_d=clCreateKernel(fftProgramd1,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_k_d1_d kernel!");
  }
  if(fft_size[2]!=0){
    sprintf(kernel_name,"fftKernel_%u_d",fft_size[2]);
    kernels->fft_kernel_d2_d=clCreateKernel(fftProgramd2,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d2_d kernel!");
    sprintf(kernel_name,"fftKernel_%u_r_d",fft_size[2]);
    kernels->fft_kernel_d2_r_d=clCreateKernel(fftProgramd2,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d2_r_d kernel!");
    sprintf(kernel_name,"fftKernel_k_%u_d",fft_size[2]);
    kernels->fft_kernel_k_d2_d=clCreateKernel(fftProgramd2,kernel_name,&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_k_d2_d kernel!");
  }
}


void build_fft_programs(cl_context * context){
    struct bigdft_device_infos infos;
    get_context_devices_infos(context, &infos);
    cl_int ciErrNum = CL_SUCCESS;
    cl_image_format format = { CL_RGBA, CL_UNSIGNED_INT32 };
    fft_code * c;
    if(fft_size[0]!=0){
      c = generate_fft_program(fft_size[0],&infos);
      printf("%s\n",c->code);
      fftProgramd0 = clCreateProgramWithSource(*context,1,(const char**) &(c->code), NULL, &ciErrNum);
      oclErrorCheck(ciErrNum,"Failed to create programd0!");
      ciErrNum = clBuildProgram(fftProgramd0, 0, NULL, "-cl-mad-enable", NULL, NULL);
      if (ciErrNum != CL_SUCCESS)
      {
          fprintf(stderr,"Error %d: Failed to build fft program d0!\n",ciErrNum);
          char cBuildLog[10240];
          clGetProgramBuildInfo(fftProgramd0, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	  fprintf(stderr,"%s\n",cBuildLog);
          exit(1);
      }
      if(!use_constant_memory) {
        cossind0 = clCreateImage2D(*context,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR , &format,fft_size[0],1,0,c->cossin,&ciErrNum);
        if (ciErrNum != CL_SUCCESS)
          fprintf(stderr,"Error %d: Failed to allocate image buffer cossind0!\n",ciErrNum);
      }
      if(c->cossin) free(c->cossin);
      free(c->code);
      free(c);
    }
    if(fft_size[1]!=0){
      c = generate_fft_program(fft_size[1], &infos);
      printf("%s\n",c->code);
      fftProgramd1 = clCreateProgramWithSource(*context,1,(const char**) &(c->code), NULL, &ciErrNum);
      oclErrorCheck(ciErrNum,"Failed to create programd1!");
      ciErrNum = clBuildProgram(fftProgramd1, 0, NULL, "-cl-mad-enable", NULL, NULL);
      if (ciErrNum != CL_SUCCESS)
      {
          fprintf(stderr,"Error %d: Failed to build fft program d1!\n",ciErrNum);
          char cBuildLog[10240];
          clGetProgramBuildInfo(fftProgramd1, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	  fprintf(stderr,"%s\n",cBuildLog);
          exit(1);
      }
      if(!use_constant_memory) {
        cossind1 = clCreateImage2D(*context,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR , &format,fft_size[1],1,0,c->cossin,&ciErrNum);
        if (ciErrNum != CL_SUCCESS)
          fprintf(stderr,"Error %d: Failed to allocate image buffer cossind1!\n",ciErrNum);
      }
      if(c->cossin) free(c->cossin);
      free(c->code);
      free(c);
    }
    if(fft_size[2]!=0){
      c = generate_fft_program(fft_size[2], &infos);
      printf("%s\n",c->code);
      fftProgramd2 = clCreateProgramWithSource(*context,1,(const char**) &(c->code), NULL, &ciErrNum);
      oclErrorCheck(ciErrNum,"Failed to create programd1!");
      ciErrNum = clBuildProgram(fftProgramd2, 0, NULL, "-cl-mad-enable", NULL, NULL);
      if (ciErrNum != CL_SUCCESS)
      {
          fprintf(stderr,"Error %d: Failed to build fft program d2!\n",ciErrNum);
          char cBuildLog[10240];
          clGetProgramBuildInfo(fftProgramd2, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	  fprintf(stderr,"%s\n",cBuildLog);
          exit(1);
      }
      if(!use_constant_memory) {
        cossind2 = clCreateImage2D(*context,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR , &format,fft_size[2],1,0,c->cossin,&ciErrNum);
        if (ciErrNum != CL_SUCCESS)
          fprintf(stderr,"Error %d: Failed to allocate image buffer cossind2!\n",ciErrNum);
      }
      if(c->cossin) free(c->cossin);
      free(c->code);
      free(c);
    }
}

void clean_fft_kernels(struct bigdft_kernels * kernels){
  cl_int ciErrNum;
  if(fft_size[0]!=0){
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d0_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d0_r_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_k_d0_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
  }
  if(fft_size[1]!=0){
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d1_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d1_r_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_k_d1_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
  }
  if(fft_size[2]!=0){
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d2_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_d2_r_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
    ciErrNum = clReleaseKernel(kernels->fft_kernel_k_d2_d);
    oclErrorCheck(ciErrNum,"Failed to release kernel!");
  }
}

void clean_fft_programs(){
  cl_int ciErrNum;
  if(fft_size[0]!=0){
    if(!use_constant_memory){
      ciErrNum = clReleaseMemObject (cossind0);
      oclErrorCheck(ciErrNum,"Failed to release buffer!");
    }
    ciErrNum = clReleaseProgram(fftProgramd0);
    oclErrorCheck(ciErrNum,"Failed to release program!");
  }
  if(fft_size[1]!=0){
    if(!use_constant_memory){
      ciErrNum = clReleaseMemObject (cossind1);
      oclErrorCheck(ciErrNum,"Failed to release buffer!");
    }
    ciErrNum = clReleaseProgram(fftProgramd1);
    oclErrorCheck(ciErrNum,"Failed to release program!");
  }
  if(fft_size[2]!=0){
    if(!use_constant_memory){
      ciErrNum = clReleaseMemObject (cossind2);
      oclErrorCheck(ciErrNum,"Failed to release buffer!");
    }
    ciErrNum = clReleaseProgram(fftProgramd2);
    oclErrorCheck(ciErrNum,"Failed to release program!");
  }
}
