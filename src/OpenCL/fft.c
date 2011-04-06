#include "OpenCL_wrappers.h"

char * fft_program="\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#define TWOPI 3.1415926535898*2\n\
#define LINE_NUMBER 16\n\
#define BUFFER_DEPTH LINE_NUMBER+1\n\
#define FFT_LENGTH 16\n\
#define radix2m(il,jl,N,A,B,in,out) \
{ \
  double2 tmp,val,w;\
  int a,b,p,r;\
  b = jl / (2*A);\
  r = jl % (2*A);\
/*  p = r / A;*/\
  a = r % A;\
  val = in[A*b+a][il];\
  tmp.x = val.x;\
  tmp.y = val.y;\
  val = in[(N/2)+A*b+a][il];\
  w.x = cos((TWOPI/(A*2))*r);\
  w.y = sin((TWOPI/(A*2))*r);\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  out[jl][il]=tmp;\
} \n\
inline void radix2(size_t il, size_t jl, uint N, uint A, uint B, __local double2 in[FFT_LENGTH][BUFFER_DEPTH], __local double2 out[FFT_LENGTH][BUFFER_DEPTH]) {\n\
   double2 tmp, val,w;\n\
   int a,b,r;\n\
   b = jl / (2*A);\n\
   r = jl % (2*A);\n\
   a = r % A;\n\
\n\
  val = in[A*b+a][il];\n\
  tmp.x = val.x;\n\
  tmp.y = val.y;\n\
  val = in[(N/2)+A*b+a][il];\n\
  w.x = cos((TWOPI/(A*2))*r);\n\
  w.y = sin((TWOPI/(A*2))*r);\n\
  tmp.x += val.x * w.x;\n\
  tmp.x += val.y * w.y;\n\
  tmp.y -= val.x * w.y;\n\
  tmp.y += val.y * w.x;\n\
  out[jl][il]=tmp;\n\
   \n\
}\n\
\n\
inline void radix4(size_t il, size_t jl, uint N, uint A, uint B, __local double2 in[FFT_LENGTH][BUFFER_DEPTH], __local double2 out[FFT_LENGTH][BUFFER_DEPTH]) {\n\
   unsigned int i;\n\
   double2 tmp=(double2)(0.0, 0.0);\n\
   double2 val;\n\
   double wx,wy;\n\
   int a,b,p,r;\n\
   b = jl/(4*A);\n\
   r = (jl - b*4*A);\n\
   p = r/A;\n\
   a = r - A*p;\n\
\n\
   for(i=0; i<4; i++){\n\
     val = in[(N/4)*i+A*b+a][il];\n\
     wx = cos(TWOPI*r*i/((double)A*4));\n\
     wy = -sin(TWOPI*r*i/((double)A*4));\n\
     tmp.x += val.x * wx - val.y * wy;\n\
     tmp.y += val.x * wy + val.y * wx;\n\
   }\n\
\n\
   out[jl][il]=tmp;\n\
   \n\
}\n\
inline void radix16(size_t il, size_t jl, uint N, uint A, uint B, __local double2 in[FFT_LENGTH][BUFFER_DEPTH], __local double2 out[FFT_LENGTH][BUFFER_DEPTH]) {\n\
   unsigned int i;\n\
   double2 tmp=(double2)(0.0, 0.0);\n\
   double2 val;\n\
   double wx,wy;\n\
   int a,b,p,r;\n\
   b = jl/(16*A);\n\
   r = (jl - b*16*A);\n\
   p = r/A;\n\
   a = r - A*p;\n\
\n\
   for(i=0; i<16; i++){\n\
     val = in[(N/16)*i+A*b+a][il];\n\
     wx = cos(TWOPI*r*i/((double)A*16));\n\
     wy = -sin(TWOPI*r*i/((double)A*16));\n\
     tmp.x += val.x * wx - val.y * wy;\n\
     tmp.y += val.x * wy + val.y * wx;\n\
   }\n\
\n\
   out[jl][il]=tmp;\n\
   \n\
}\n\
\n\
__kernel void fftKernel_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out){\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 (*tmp_in)[BUFFER_DEPTH];\n\
__local double2 (*tmp_out)[BUFFER_DEPTH];\n\
__local double2 (*tmp)[BUFFER_DEPTH];\n\
\n\
  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t jlt = il;\n\
  size_t ilt = jl;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + il;\n\
  tmp1[jl][il] = jl < n ? psi[jgt + ( jl ) * ndat] : 0.0;\n\
  \n\
  unsigned int A,B;\n\
  unsigned int j;\n\
  A=1;\n\
  B=FFT_LENGTH;\n\
  tmp_in=tmp1;\n\
  tmp_out=tmp2;\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  #pragma unroll\n\
  for(j=0;j<4;j++){\n\
    B /= 2;\n\
    radix2(il, jl, n, A, B, tmp_in, tmp_out);\n\
//    tmp_out[il][jl]=tmp_in[il][jl];\n\
    tmp = tmp_in; tmp_in = tmp_out; tmp_out = tmp;\n\
    A *= 2;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
  }\n\
/*    radix2m(il, jl, n, 1, 8, tmp_in, tmp_out);\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    radix2m(il, jl, n, 2, 4, tmp_out, tmp_in);\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    radix2m(il, jl, n, 4, 2, tmp_in, tmp_out);\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    radix2m(il, jl, n, 8, 1, tmp_out, tmp_in);\n\
    barrier(CLK_LOCAL_MEM_FENCE);*/\n\
/*  for(j=0;j<2;j++){\n\
    B /= 4;\n\
    radix4(il, jl, n, A, B, tmp_in, tmp_out);\n\
    tmp = tmp_in; tmp_in = tmp_out; tmp_out = tmp;\n\
    A *= 4;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
  }*/\n\
/*  for(j=0;j<1;j++){\n\
    B /= 16;\n\
    radix16(il, jl, n, A, B, tmp_in, tmp_out);\n\
    tmp = tmp_in; tmp_in = tmp_out; tmp_out = tmp;\n\
    A *= 16;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
  }*/\n\
\n\
  if(il<n)\n\
    out[jg*n+il] = tmp_in[il][jl];\n\
}";


inline void fft_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
#define FFT_LENGTH 16
#define BUFFER_DEPTH 16
    assert(*n==FFT_LENGTH);
    size_t block_size_i=FFT_LENGTH, block_size_j=BUFFER_DEPTH;
    cl_uint i = 0;
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue fft kernel!");
}

void FC_FUNC_(fft1d_d,FFT1D_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    fft_generic((*command_queue)->kernels.fft_kernel_d, (*command_queue)->command_queue, n, ndat, psi, out);
}

void FC_FUNC_(fft3d_d,FFT1D_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,cl_mem *psi,cl_mem *out,cl_mem *tmp){
    cl_uint n1 = dimensions[0];
    cl_uint n2 = dimensions[1];
    cl_uint n3 = dimensions[2];
    cl_uint ndat = n1 * n2;
    fft_generic((*command_queue)->kernels.fft_kernel_d, (*command_queue)->command_queue, &n3, &ndat, psi, out);
    ndat = n1 * n3;
    fft_generic((*command_queue)->kernels.fft_kernel_d, (*command_queue)->command_queue, &n2, &ndat, out, tmp);
    ndat = n2 * n3;
    fft_generic((*command_queue)->kernels.fft_kernel_d, (*command_queue)->command_queue, &n1, &ndat, tmp, out);
}

cl_program fftProgram;

void create_fft_kernels(struct bigdft_kernels * kernels){
    cl_int ciErrNum = CL_SUCCESS;
    kernels->fft_kernel_d=clCreateKernel(fftProgram,"fftKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create fftKernel_d kernel!");
}

void build_fft_programs(cl_context * context){
    cl_int ciErrNum = CL_SUCCESS;
    fftProgram = clCreateProgramWithSource(*context,1,(const char**) &fft_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(fftProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to build fft program!\n",ciErrNum);
        char cBuildLog[10240];
        clGetProgramBuildInfo(fftProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
}

void clean_fft_kernels(struct bigdft_kernels * kernels){
  cl_int ciErrNum;
  ciErrNum = clReleaseKernel(kernels->fft_kernel_d);
  oclErrorCheck(ciErrNum,"Failed to release kernel!");
}

void clean_fft_programs(){
  cl_int ciErrNum;
  ciErrNum = clReleaseProgram(fftProgram);
  oclErrorCheck(ciErrNum,"Failed to release program!");
}
