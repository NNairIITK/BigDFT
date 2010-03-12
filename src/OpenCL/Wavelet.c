#include "Wavelet.h"
#include "OpenCL_wrappers.h"

char * ana1d_program="\
#define FILTER_WIDTH 16\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
__kernel void anashrink1dKernel_d(size_t n, size_t ndat, __global const double *psi, __global double *out, __local double tmp[]){\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1);\n\
size_t igt = get_group_id(0);\n\
size_t jgt = get_group_id(1);\n\
size_t jb;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
igt = ig - i2 + j2;\n\
jgt = jg - j2 + i2;\n\
psi += 7*ndat;\n\
//If I'm on the outside, select a border element to load\n\
if (j2 < 7) {\n\
    jb =  2 * igt - j2 - 7;\n\
    tmp[i2 * (3 * FILTER_WIDTH + 1) + j2]=psi[jgt+jb*ndat];\n\
  }\n\
if (j2 >= 9) {\n\
    jb = 2 * igt + FILTER_WIDTH - j2 + 7;\n\
    tmp[i2 * (3 * FILTER_WIDTH + 1) + j2 + 2 * FILTER_WIDTH - 2]=psi[jgt+jb*ndat];\n\
  }\n\
//Load the elements I am to calculate\n\
tmp[i2 * (3 * FILTER_WIDTH + 1) + 2 * j2 + 7]=psi[jgt+igt*2*ndat];\n\
tmp[i2 * (3 * FILTER_WIDTH + 1) + 2 * j2 + 7 + 1]=psi[jgt+(igt*2+1)*ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
double ci = 0.0;\n\
double di = 0.0;\n\
tmp = tmp + j2*(3*FILTER_WIDTH+1) + 7 + 2*i2;\n\
ci += tmp[ 7] * -0.00030292051472413308126;\n\
di += tmp[ 7] *  0.00054213233180001068935;\n\
di += tmp[-6] * -0.00030292051472413308126;\n\
ci += tmp[-6] * -0.00054213233180001068935;\n\
ci += tmp[ 8] *  0.0018899503327676891843;\n\
di += tmp[ 8] * -0.0033824159510050025955;\n\
di += tmp[-7] * -0.0018899503327676891843;\n\
ci += tmp[-7] * -0.0033824159510050025955;\n\
ci += tmp[ 5] *  0.0038087520138944894631;\n\
di += tmp[ 5] * -0.0076074873249766081919;\n\
di += tmp[-4] *  0.0038087520138944894631;\n\
ci += tmp[-4] *  0.0076074873249766081919;\n\
ci += tmp[ 6] * -0.014952258337062199118;\n\
di += tmp[ 6] *  0.031695087811525991431;\n\
di += tmp[-5] *  0.014952258337062199118;\n\
ci += tmp[-5] *  0.031695087811525991431;\n\
ci += tmp[ 3] * -0.027219029917103486322;\n\
di += tmp[ 3] *  0.061273359067811077843;\n\
di += tmp[-2] * -0.027219029917103486322;\n\
ci += tmp[-2] * -0.061273359067811077843;\n\
di += tmp[-3] * -0.049137179673730286787;\n\
ci += tmp[-3] * -0.14329423835127266284;\n\
ci += tmp[ 2] * -0.051945838107881800736;\n\
di += tmp[ 2] *  0.48135965125905339159;\n\
di += tmp[-1] *  0.051945838107881800736;\n\
ci += tmp[-1] *  0.48135965125905339159;\n\
ci += tmp[ 4] *  0.049137179673730286787;\n\
di += tmp[ 4] * -0.14329423835127266284;\n\
ci += tmp[ 1] *  0.36444189483617893676;\n\
di += tmp[ 1] * -0.77718575169962802862;\n\
di += tmp[ 0] *  0.36444189483617893676;\n\
ci += tmp[ 0] *  0.77718575169962802862;\n\
\
\
out[(jg*(2*n)+ig)]=ci;\n\
out[(jg*(2*n)+ig+n)]=di;\n\
};\n\
__kernel void ana1dKernel_d(size_t n, size_t ndat, __global const double *psi, __global double *out, __local double tmp[]){\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1);\n\
ptrdiff_t igt = get_group_id(0);\n\
ptrdiff_t jgt = get_group_id(1);\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
igt = (ig - i2) * 2 + j2 - FILTER_WIDTH/2;\n\
jgt = jg - j2 + i2;\n\
if ( igt < 0 ) \n\
  tmp[i2 * (3 * FILTER_WIDTH + 1) + j2] = psi[jgt + ( 2*n + igt ) * ndat];\n\
else \n\
  tmp[i2 * (3 * FILTER_WIDTH + 1) + j2] = psi[jgt + igt * ndat];\n\
igt += FILTER_WIDTH;\n\
tmp[i2 * (3 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = psi[jgt + igt * ndat];\n\
igt += FILTER_WIDTH;\n\
if ( igt >= 2*n ) \n\
  tmp[i2 * (3 * FILTER_WIDTH + 1) + j2 + 2*FILTER_WIDTH] = psi[jgt + ( igt - 2*n ) * ndat];\n\
else\n\
  tmp[i2 * (3 * FILTER_WIDTH + 1) + j2 + 2*FILTER_WIDTH] = psi[jgt +  igt * ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
double ci = 0.0;\n\
double di = 0.0;\n\
tmp = tmp + j2*(3*FILTER_WIDTH+1) + FILTER_WIDTH/2+2*i2;\n\
ci += tmp[ 7] * -0.00030292051472413308126;\n\
di += tmp[ 7] *  0.00054213233180001068935;\n\
di += tmp[-6] * -0.00030292051472413308126;\n\
ci += tmp[-6] * -0.00054213233180001068935;\n\
ci += tmp[ 8] *  0.0018899503327676891843;\n\
di += tmp[ 8] * -0.0033824159510050025955;\n\
di += tmp[-7] * -0.0018899503327676891843;\n\
ci += tmp[-7] * -0.0033824159510050025955;\n\
ci += tmp[ 5] *  0.0038087520138944894631;\n\
di += tmp[ 5] * -0.0076074873249766081919;\n\
di += tmp[-4] *  0.0038087520138944894631;\n\
ci += tmp[-4] *  0.0076074873249766081919;\n\
ci += tmp[ 6] * -0.014952258337062199118;\n\
di += tmp[ 6] *  0.031695087811525991431;\n\
di += tmp[-5] *  0.014952258337062199118;\n\
ci += tmp[-5] *  0.031695087811525991431;\n\
ci += tmp[ 3] * -0.027219029917103486322;\n\
di += tmp[ 3] *  0.061273359067811077843;\n\
di += tmp[-2] * -0.027219029917103486322;\n\
ci += tmp[-2] * -0.061273359067811077843;\n\
di += tmp[-3] * -0.049137179673730286787;\n\
ci += tmp[-3] * -0.14329423835127266284;\n\
ci += tmp[ 2] * -0.051945838107881800736;\n\
di += tmp[ 2] *  0.48135965125905339159;\n\
di += tmp[-1] *  0.051945838107881800736;\n\
ci += tmp[-1] *  0.48135965125905339159;\n\
ci += tmp[ 4] *  0.049137179673730286787;\n\
di += tmp[ 4] * -0.14329423835127266284;\n\
ci += tmp[ 1] *  0.36444189483617893676;\n\
di += tmp[ 1] * -0.77718575169962802862;\n\
di += tmp[ 0] *  0.36444189483617893676;\n\
ci += tmp[ 0] *  0.77718575169962802862;\n\
\
\
out[(jg*(2*n)+ig)]=ci;\n\
out[(jg*(2*n)+ig+n)]=di;\n\
};\n\
";

char * syn1d_program="\
#define FILTER_WIDTH 8\n\
#define SIZE_I 16\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
__kernel void syngrow1dKernel_d(size_t n, size_t ndat, __global const double *psi, __global double *out, __local double tmp_1[]){\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
{\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1);\n\
ptrdiff_t igt = get_group_id(0);\n\
ptrdiff_t jgt = get_group_id(1);\n\
ptrdiff_t jb;\n\
size_t ioff;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
igt = ig - i2 + j2 - FILTER_WIDTH/2 - 4;\n\
jgt = jg - j2 + i2;\n\
//If I'm on the outside, select a border element to load\n\
ioff = i2*(2*FILTER_WIDTH+2*SIZE_I+1) + j2;\n\
if( igt < 0 ) {\n\
  tmp_1[ioff] = 0.0;\n\
  tmp_1[ioff+FILTER_WIDTH+SIZE_I] = 0.0;\n\
} else {\n\
  tmp_1[ioff] = psi[jgt+igt*ndat];\n\
  tmp_1[ioff+FILTER_WIDTH+SIZE_I] = psi[jgt+(igt+n-7)*ndat];\n\
}\n\
igt += SIZE_I;\n\
if( j2 < SIZE_I - FILTER_WIDTH){\n\
  if ( igt >=n-7 ) {\n\
    tmp_1[ioff+SIZE_I] = 0.0;\n\
    tmp_1[ioff+FILTER_WIDTH+2*SIZE_I] = 0.0;\n\
  } else {\n\
    tmp_1[ioff+SIZE_I] = psi[jgt+igt*ndat];\n\
    tmp_1[ioff+FILTER_WIDTH+2*SIZE_I] = psi[jgt+(igt+n-7)*ndat];\n\
  }\n\
}\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
ioff = j2*(2*FILTER_WIDTH+2*SIZE_I+1) + FILTER_WIDTH/2+i2;\n\
tmp_1 = tmp_1 + ioff;\n\
}\n\
double se = 0.0;\n\
double so = 0.0;\n\
se += tmp_1[FILTER_WIDTH+SIZE_I + 4] * -0.00030292051472413308126;\n\
so += tmp_1[FILTER_WIDTH+SIZE_I + 3] *  0.014952258337062199118;\n\
se += tmp_1[-2] * -0.014952258337062199118;\n\
so += tmp_1[-3] * -0.00030292051472413308126;\n\
se += tmp_1[ 4] * -0.00054213233180001068935;\n\
so += tmp_1[ 3] *  0.031695087811525991431;\n\
se += tmp_1[FILTER_WIDTH+SIZE_I - 2] *  0.031695087811525991431;\n\
so += tmp_1[FILTER_WIDTH+SIZE_I - 3] *  0.00054213233180001068935;\n\
se += tmp_1[FILTER_WIDTH+SIZE_I + 3] *  0.0038087520138944894631;\n\
so += tmp_1[FILTER_WIDTH+SIZE_I + 2] * -0.049137179673730286787;\n\
se += tmp_1[-1] *  0.049137179673730286787;\n\
so += tmp_1[-2] *  0.0038087520138944894631;\n\
se += tmp_1[FILTER_WIDTH+SIZE_I - 1] * -0.14329423835127266284;\n\
so += tmp_1[FILTER_WIDTH+SIZE_I - 2] * -0.0076074873249766081919;\n\
se += tmp_1[ 3] *  0.0076074873249766081919;\n\
so += tmp_1[ 2] * -0.14329423835127266284;\n\
se += tmp_1[FILTER_WIDTH+SIZE_I + 2] * -0.027219029917103486322;\n\
so += tmp_1[FILTER_WIDTH+SIZE_I + 1] *  0.051945838107881800736;\n\
se += tmp_1[ 0] * -0.051945838107881800736;\n\
so += tmp_1[-1] * -0.027219029917103486322;\n\
se += tmp_1[ 2] * -0.061273359067811077843;\n\
so += tmp_1[ 1] *  0.48135965125905339159;\n\
se += tmp_1[FILTER_WIDTH+SIZE_I - 0] *  0.48135965125905339159;\n\
so += tmp_1[FILTER_WIDTH+SIZE_I - 1] *  0.061273359067811077843;\n\
se += tmp_1[FILTER_WIDTH+SIZE_I + 1] *  0.36444189483617893676;\n\
so += tmp_1[FILTER_WIDTH+SIZE_I + 0] * -0.77718575169962802862;\n\
se += tmp_1[ 1] *  0.77718575169962802862;\n\
so += tmp_1[ 0] *  0.36444189483617893676;\n\
so += tmp_1[FILTER_WIDTH+SIZE_I + 4] * -0.0018899503327676891843;\n\
so += tmp_1[ 4] * -0.0033824159510050025955;\n\
se += tmp_1[-3] *  0.0018899503327676891843;\n\
se += tmp_1[FILTER_WIDTH+SIZE_I - 3] * -0.0033824159510050025955;\n\
\
out[jg*2*n+2*ig]=so;\n\
out[jg*2*n+2*ig+1]=se;\n\
};\n\
__kernel void syn1dKernel_d(size_t n, size_t ndat, __global const double *psi, __global double *out, __local double tmp_1[]){\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
{\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1);\n\
ptrdiff_t igt = get_group_id(0);\n\
ptrdiff_t jgt = get_group_id(1);\n\
size_t ioff;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
igt = ig - i2 + j2 - FILTER_WIDTH/2;\n\
jgt = jg - j2 + i2;\n\
//If I'm on the outside, select a border element to load\n\
ioff = i2*(2*FILTER_WIDTH+2*SIZE_I+1) + j2;\n\
if( igt < 0 ) {\n\
  tmp_1[ioff] = psi[jgt+(n+igt)*ndat];\n\
  tmp_1[ioff+FILTER_WIDTH+SIZE_I] = psi[jgt+(n+igt+n)*ndat];\n\
} else {\n\
  tmp_1[ioff] = psi[jgt+igt*ndat];\n\
  tmp_1[ioff+FILTER_WIDTH+SIZE_I] = psi[jgt+(igt+n)*ndat];\n\
}\n\
igt += SIZE_I;\n\
if( j2 < SIZE_I - FILTER_WIDTH){\n\
  if ( igt >=n ) {\n\
    tmp_1[ioff+SIZE_I] = psi[jgt+(igt-n)*ndat];\n\
    tmp_1[ioff+FILTER_WIDTH+2*SIZE_I] = psi[jgt+(igt-n+n)*ndat];\n\
  } else {\n\
    tmp_1[ioff+SIZE_I] = psi[jgt+igt*ndat];\n\
    tmp_1[ioff+FILTER_WIDTH+2*SIZE_I] = psi[jgt+(igt+n)*ndat];\n\
  }\n\
}\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
ioff = j2*(2*FILTER_WIDTH+2*SIZE_I+1) + FILTER_WIDTH/2+i2;\n\
tmp_1 = tmp_1 + ioff;\n\
}\n\
double se = 0.0;\n\
double so = 0.0;\n\
se += tmp_1[FILTER_WIDTH+SIZE_I + 3] * -0.00030292051472413308126;\n\
so += tmp_1[FILTER_WIDTH+SIZE_I + 3] *  0.014952258337062199118;\n\
se += tmp_1[-3] * -0.014952258337062199118;\n\
so += tmp_1[-3] * -0.00030292051472413308126;\n\
se += tmp_1[ 3] * -0.00054213233180001068935;\n\
so += tmp_1[ 3] *  0.031695087811525991431;\n\
se += tmp_1[FILTER_WIDTH+SIZE_I - 3] *  0.031695087811525991431;\n\
so += tmp_1[FILTER_WIDTH+SIZE_I - 3] *  0.00054213233180001068935;\n\
se += tmp_1[FILTER_WIDTH+SIZE_I + 2] *  0.0038087520138944894631;\n\
so += tmp_1[FILTER_WIDTH+SIZE_I + 2] * -0.049137179673730286787;\n\
se += tmp_1[-2] *  0.049137179673730286787;\n\
so += tmp_1[-2] *  0.0038087520138944894631;\n\
se += tmp_1[FILTER_WIDTH+SIZE_I - 2] * -0.14329423835127266284;\n\
so += tmp_1[FILTER_WIDTH+SIZE_I - 2] * -0.0076074873249766081919;\n\
se += tmp_1[ 2] *  0.0076074873249766081919;\n\
so += tmp_1[ 2] * -0.14329423835127266284;\n\
se += tmp_1[FILTER_WIDTH+SIZE_I + 1] * -0.027219029917103486322;\n\
so += tmp_1[FILTER_WIDTH+SIZE_I + 1] *  0.051945838107881800736;\n\
se += tmp_1[-1] * -0.051945838107881800736;\n\
so += tmp_1[-1] * -0.027219029917103486322;\n\
se += tmp_1[ 1] * -0.061273359067811077843;\n\
so += tmp_1[ 1] *  0.48135965125905339159;\n\
se += tmp_1[FILTER_WIDTH+SIZE_I - 1] *  0.48135965125905339159;\n\
so += tmp_1[FILTER_WIDTH+SIZE_I - 1] *  0.061273359067811077843;\n\
se += tmp_1[FILTER_WIDTH+SIZE_I + 0] *  0.36444189483617893676;\n\
so += tmp_1[FILTER_WIDTH+SIZE_I + 0] * -0.77718575169962802862;\n\
se += tmp_1[ 0] *  0.77718575169962802862;\n\
so += tmp_1[ 0] *  0.36444189483617893676;\n\
so += tmp_1[FILTER_WIDTH+SIZE_I + 4] * -0.0018899503327676891843;\n\
so += tmp_1[ 4] * -0.0033824159510050025955;\n\
se += tmp_1[-4] *  0.0018899503327676891843;\n\
se += tmp_1[FILTER_WIDTH+SIZE_I - 4] * -0.0033824159510050025955;\n\
\
out[jg*(2*n)+ig*2]=se;\n\
out[jg*(2*n)+ig*2+1]=so;\n\
};\n\
";

inline void ana_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out){
    cl_int ciErrNum;
    int FILTER_WIDTH = 16;
    assert(*n>=FILTER_WIDTH);
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;
    cl_uint i = 0;
    clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    clSetKernelArg(kernel, i++,sizeof(double)*block_size_j*(block_size_i*2+FILTER_WIDTH + 1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue analysis kernel!");

}

inline void syn_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out) {
    cl_int ciErrNum;
    int FILTER_WIDTH = 8;
    int SIZE_I = 16;
    assert(*n>=SIZE_I);
    size_t block_size_i=SIZE_I, block_size_j=SIZE_I;
    cl_uint i = 0;
    clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    clSetKernelArg(kernel, i++,sizeof(double)*block_size_j*(block_size_i*2+FILTER_WIDTH*2+1), NULL);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue synthesis kernel!");
}

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
