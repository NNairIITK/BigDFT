#include <CL/cl.h>
#include <stdio.h>
#include <stdlib.h>
#include <config.h>
#include <math.h>

#define DEBUG 0

char * magicfilter1d_program="\
#define FILTER_WIDTH 16\n\
//n is supposed to be greater or equal than get_local_size(0)\n\
__kernel void magicfilter1dKernel_l(size_t n, size_t ndat, __global const float *psi, __global float *out, __local float tmp[]){\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
size_t i2 = get_local_id(0);\n\
size_t j2 = get_local_id(1);\n\
size_t is = get_local_size(0);\n\
size_t ib;\n\
size_t it;\n\
size_t base_i;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
if( jg >= ndat ) return;\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
__local float * tmp_o = tmp + j2*(is+FILTER_WIDTH);\n\
base_i = FILTER_WIDTH/2+i2;\n\
//If I'm on the outside, select a border element to load\n\
if(i2 < FILTER_WIDTH/2)\n\
  { it = i2;\n\
    if (ig < FILTER_WIDTH/2)\n\
      { ib =  n + i2 - FILTER_WIDTH/2; }\n\
    else { ib = ig - FILTER_WIDTH/2; }\n\
    tmp_o[it]=psi[jg+ib*ndat];\n\
  }\n\
if (i2 >= (is - FILTER_WIDTH/2) || (ig >= n - FILTER_WIDTH/2))\n\
  { it = i2 + FILTER_WIDTH;\n\
    if (ig >= n - FILTER_WIDTH/2)\n\
      { ib = ig - n + FILTER_WIDTH/2; }\n\
    else { ib = ig + FILTER_WIDTH/2; }\n\
    tmp_o[it]=psi[jg+ib*ndat];\n\
  }\n\
//check boundaries\
//Load the element I am to calculate\n\
tmp_o[base_i]=psi[jg+ig*ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
float tt = 0.0;\n\
tmp_o = tmp_o + base_i;\n\
tt += tmp_o[-8] *  8.4334247333529341094733325815816e-7\n\
    + tmp_o[+7] *  2.72734492911979659657715313017228e-6\n\
    + tmp_o[-7] * -0.1290557201342060969516786758559028e-4\n\
    + tmp_o[+6] * -0.5185986881173432922848639136911487e-4\n\
    + tmp_o[-6] *  0.8762984476210559564689161894116397e-4\n\
    + tmp_o[-5] * -0.30158038132690463167163703826169879e-3\n\
    + tmp_o[+5] *  0.49443227688689919192282259476750972e-3\n\
    + tmp_o[-4] *  0.174723713672993903449447812749852942e-2;\n\
tt += tmp_o[+4] * -0.344128144493493857280881509686821861e-2\n\
    + tmp_o[-3] * -0.942047030201080385922711540948195075e-2\n\
    + tmp_o[+3] *  0.1337263414854794752733423467013220997e-1\n\
    + tmp_o[+2] * -0.2103025160930381434955489412839065067e-1\n\
    + tmp_o[-2] *  0.2373821463724942397566389712597274535e-1\n\
    + tmp_o[+1] * -0.604895289196983516002834636e-1\n\
    + tmp_o[-1] *  0.612625895831207982195380597e-1\n\
    + tmp_o[+0] *  0.9940415697834003993178616713;\n\
out[(jg*n+ig)]=tt;\n\
};\
";

char * ana1d_program="\
#define FILTER_WIDTH 16\n\
__kernel void ana1dKernel_l(size_t n, size_t ndat, __global const float *psi, __global float *out, __local float tmp[]){\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
size_t i2 = get_local_id(0);\n\
size_t j2 = get_local_id(1);\n\
size_t is = get_local_size(0);\n\
size_t ib;\n\
size_t it;\n\
size_t base_i;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
if( jg >= ndat ) return;\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
__local float * tmp_o = tmp + j2*(2*is+FILTER_WIDTH);\n\
base_i = FILTER_WIDTH/2+2*i2;\n\
//If I'm on the outside, select a border element to load\n\
if(i2 < FILTER_WIDTH/2)\n\
  { it = i2;\n\
    if (ig < FILTER_WIDTH/2)\n\
      { ib = 2 * n  - ( FILTER_WIDTH/2 - i2 ); }\n\
    else { ib =  2 * ig - i2 - FILTER_WIDTH/2; }// 2 * ( ig - i2 ) - ( FILTER_WIDTH/2 - i2 )\n\
    tmp_o[it]=psi[jg+ib*ndat];\n\
  }\n\
if (i2 >= (is - FILTER_WIDTH/2) || (ig >= n - FILTER_WIDTH/2))\n\
  { it = is + i2 + FILTER_WIDTH;\n\
    if (ig >= n - FILTER_WIDTH/2)\n\
      { ib = ig - n + FILTER_WIDTH/2; }\n\
    else { ib = 2 * ig +  is - i2 + FILTER_WIDTH/2; }// 2 * ( ig + is - i2 ) + ( FILTER_WIDTH/2 - ( is - i2 ) )\n\
    tmp_o[it]=psi[jg+ib*ndat];\n\
  }\n\
//Load the elements I am to calculate\n\
tmp_o[base_i]=psi[jg+ig*2*ndat];\n\
tmp_o[base_i+1]=psi[jg+(ig*2+1)*ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
float ci = 0.0;\n\
float di = 0.0;\n\
tmp_o = tmp_o + base_i;\n\
ci += tmp_o[ 7] * -0.00030292051472413308126;\n\
di += tmp_o[ 7] *  0.00054213233180001068935;\n\
di += tmp_o[-6] * -0.00030292051472413308126;\n\
ci += tmp_o[-6] * -0.00054213233180001068935;\n\
ci += tmp_o[ 8] *  0.0018899503327676891843;\n\
di += tmp_o[ 8] * -0.0033824159510050025955;\n\
di += tmp_o[-7] * -0.0018899503327676891843;\n\
ci += tmp_o[-7] * -0.0033824159510050025955;\n\
ci += tmp_o[ 5] *  0.0038087520138944894631;\n\
di += tmp_o[ 5] * -0.0076074873249766081919;\n\
di += tmp_o[-4] *  0.0038087520138944894631;\n\
ci += tmp_o[-4] *  0.0076074873249766081919;\n\
ci += tmp_o[ 6] * -0.014952258337062199118;\n\
di += tmp_o[ 6] *  0.031695087811525991431;\n\
di += tmp_o[-5] *  0.014952258337062199118;\n\
ci += tmp_o[-5] *  0.031695087811525991431;\n\
ci += tmp_o[ 3] * -0.027219029917103486322;\n\
di += tmp_o[ 3] *  0.061273359067811077843;\n\
di += tmp_o[-2] * -0.027219029917103486322;\n\
ci += tmp_o[-2] * -0.061273359067811077843;\n\
di += tmp_o[-3] * -0.049137179673730286787;\n\
ci += tmp_o[-3] * -0.14329423835127266284;\n\
ci += tmp_o[ 2] * -0.051945838107881800736;\n\
di += tmp_o[ 2] *  0.48135965125905339159;\n\
di += tmp_o[-1] *  0.051945838107881800736;\n\
ci += tmp_o[-1] *  0.48135965125905339159;\n\
ci += tmp_o[ 4] *  0.049137179673730286787;\n\
di += tmp_o[ 4] * -0.14329423835127266284;\n\
ci += tmp_o[ 1] *  0.36444189483617893676;\n\
di += tmp_o[ 1] * -0.77718575169962802862;\n\
di += tmp_o[ 0] *  0.36444189483617893676;\n\
ci += tmp_o[ 0] *  0.77718575169962802862;\n\
\
\
out[(jg*(2*n)+ig)]=ci;\n\
out[(jg*(2*n)+ig+n)]=di;\n\
};\n\
";

char * syn1d_program="\
#define FILTER_WIDTH 16\n\
__kernel void syn1dKernel_l(size_t n, size_t ndat, __global const float *psi, __global float *out, __local float tmp_1[], __local float tmp_2[]){\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
size_t i2 = get_local_id(0);\n\
size_t j2 = get_local_id(1);\n\
size_t is = get_local_size(0);\n\
size_t ib;\n\
size_t it;\n\
size_t base_i;\n\
if( jg >= ndat ) return;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
base_i = FILTER_WIDTH/2+i2;\n\
__local float * tmp_o_1 = tmp_1 + j2*(is+FILTER_WIDTH);\n\
__local float * tmp_o_2 = tmp_2 + j2*(is+FILTER_WIDTH);\n\
//If I'm on the outside, select a border element to load\n\
if(i2 < FILTER_WIDTH/2)\n\
  { it = i2;\n\
    if (ig < FILTER_WIDTH/2)\n\
      { ib = n - ( FILTER_WIDTH/2 - i2 ); }\n\
    else { ib = ig - FILTER_WIDTH/2; }\n\
    tmp_o_1[it]=psi[jg+ib*ndat];\n\
    tmp_o_2[it]=psi[jg+(ib+n)*ndat];\n\
  }\n\
if (i2 >= (is - FILTER_WIDTH/2) || (ig >= n - FILTER_WIDTH/2))\n\
  { it = i2 + FILTER_WIDTH;\n\
    if (ig >= n - FILTER_WIDTH/2)\n\
      { ib = ig - n + FILTER_WIDTH/2; }\n\
    else { ib = ig + FILTER_WIDTH/2; }\n\
    tmp_o_1[it]=psi[jg+ib*ndat];\n\
    tmp_o_2[it]=psi[jg+(ib+n)*ndat];\n\
  }\n\
//Load the elements I am to calculate\n\
tmp_o_1 = tmp_o_1 + base_i;\n\
tmp_o_2 = tmp_o_2 + base_i;\n\
tmp_o_1[0]=psi[jg+ig*ndat];\n\
tmp_o_2[0]=psi[jg+(ig+n)*ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
float se = 0.0;\n\
float so = 0.0;\n\
se += tmp_o_2[ 3] * -0.00030292051472413308126;\n\
so += tmp_o_2[ 3] *  0.014952258337062199118;\n\
se += tmp_o_1[ 3] * -0.00054213233180001068935;\n\
so += tmp_o_1[ 3] *  0.031695087811525991431;\n\
se += tmp_o_1[-4] *  0.0018899503327676891843;\n\
se += tmp_o_2[-4] * -0.0033824159510050025955;\n\
so += tmp_o_2[ 4] * -0.0018899503327676891843;\n\
so += tmp_o_1[ 4] * -0.0033824159510050025955;\n\
se += tmp_o_2[ 2] *  0.0038087520138944894631;\n\
so += tmp_o_2[ 2] * -0.049137179673730286787;\n\
se += tmp_o_1[ 2] *  0.0076074873249766081919;\n\
so += tmp_o_1[ 2] * -0.14329423835127266284;\n\
se += tmp_o_1[-3] * -0.014952258337062199118;\n\
so += tmp_o_1[-3] * -0.00030292051472413308126;\n\
se += tmp_o_2[ 1] * -0.027219029917103486322;\n\
so += tmp_o_2[ 1] *  0.051945838107881800736;\n\
se += tmp_o_2[-3] *  0.031695087811525991431;\n\
so += tmp_o_2[-3] *  0.00054213233180001068935;\n\
se += tmp_o_1[-2] *  0.049137179673730286787;\n\
so += tmp_o_1[-2] *  0.0038087520138944894631;\n\
se += tmp_o_1[-1] * -0.051945838107881800736;\n\
so += tmp_o_1[-1] * -0.027219029917103486322;\n\
se += tmp_o_1[ 1] * -0.061273359067811077843;\n\
so += tmp_o_1[ 1] *  0.48135965125905339159;\n\
se += tmp_o_2[-2] * -0.14329423835127266284;\n\
so += tmp_o_2[-2] * -0.0076074873249766081919;\n\
se += tmp_o_2[-1] *  0.48135965125905339159;\n\
so += tmp_o_2[-1] *  0.061273359067811077843;\n\
se += tmp_o_2[ 0] *  0.36444189483617893676;\n\
so += tmp_o_2[ 0] * -0.77718575169962802862;\n\
se += tmp_o_1[ 0] *  0.77718575169962802862;\n\
so += tmp_o_1[ 0] *  0.36444189483617893676;\n\
\
\
out[(jg*(2*n)+ig*2)]=se;\n\
out[(jg*(2*n)+ig*2+1)]=so;\n\
};\n\
";

char * c_initialize_program="\
__kernel void c_initializeKernel_l(int n, int ndat, __global const float * x_in, __global float * y_in, float c) {\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
if( jg >= ndat ) return;\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
size_t pos = jg+ig*ndat;\n\
y_in[pos] = x_in[pos] * c;\n\
};\n\
";

char * kinetic1d_program="\
#define FILTER_WIDTH 30\n\
__kernel void kinetic1dKernel_l(int n, int ndat, float scale, __global const float * x_in, __global float * x_out, __global const float * y_in, __global float * y_out, __local float * tmp) {\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
size_t i2 = get_local_id(0);\n\
size_t j2 = get_local_id(1);\n\
size_t is = get_local_size(0);\n\
size_t ib;\n\
size_t it;\n\
size_t base_i;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
if( jg >= ndat ) return;\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
__local float * tmp_o = tmp + j2*(is+FILTER_WIDTH);\n\
base_i = FILTER_WIDTH/2+i2;\n\
//If I'm on the outside, select a border element to load\n\
if(i2 < FILTER_WIDTH/2)\n\
  { it = i2;\n\
    if (ig < FILTER_WIDTH/2)\n\
      { ib =  n + i2 - FILTER_WIDTH/2; }\n\
    else { ib = ig - FILTER_WIDTH/2; }\n\
    tmp_o[it]=x_in[jg+ib*ndat];\n\
  }\n\
if (i2 >= (is - FILTER_WIDTH/2) || (ig >= n - FILTER_WIDTH/2))\n\
  { it = i2 + FILTER_WIDTH;\n\
    if (ig >= n - FILTER_WIDTH/2)\n\
      { ib = ig - n + FILTER_WIDTH/2; }\n\
    else { ib = ig + FILTER_WIDTH/2; }\n\
    tmp_o[it]=x_in[jg+ib*ndat];\n\
  }\n\
//check boundaries\
//Load the element I am to calculate\n\
tmp_o[base_i]=x_in[jg+ig*ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
tmp_o = tmp_o + base_i;\n\
float conv = 0.0;\n\
conv += (tmp_o[14] + tmp_o[-14]) * -6.924474940639200152025730585882e-18;\n\
conv += (tmp_o[13] + tmp_o[-13]) *  2.70800493626319438269856689037647576e-13;\n\
conv += (tmp_o[12] + tmp_o[-12]) * -5.813879830282540547959250667e-11;\n\
conv += (tmp_o[11] + tmp_o[-11]) * -1.05857055496741470373494132287e-8;\n\
conv += (tmp_o[10] + tmp_o[-10]) * -3.7230763047369275848791496973044e-7;\n\
conv += (tmp_o[ 9] + tmp_o[ -9]) *  2.0904234952920365957922889447361e-6;\n\
conv += (tmp_o[ 8] + tmp_o[ -8]) * -0.2398228524507599670405555359023135e-4;\n\
conv += (tmp_o[ 7] + tmp_o[ -7]) *  0.45167920287502235349480037639758496e-3;\n\
conv += (tmp_o[ 6] + tmp_o[ -6]) * -0.409765689342633823899327051188315485e-2;\n\
conv += (tmp_o[ 5] + tmp_o[ -5]) *  0.02207029188482255523789911295638968409e0;\n\
conv += (tmp_o[ 4] + tmp_o[ -4]) * -0.0822663999742123340987663521e0;\n\
conv += (tmp_o[ 3] + tmp_o[ -3]) *  0.2371780582153805636239247476e0;\n\
conv += (tmp_o[ 2] + tmp_o[ -2]) * -0.6156141465570069496314853949e0;\n\
conv += (tmp_o[ 1] + tmp_o[ -1]) *  2.2191465938911163898794546405e0;\n\
conv +=  tmp_o[ 0]             * -3.5536922899131901941296809374e0;\n\
y_out[jg*n + ig] = y_in[jg + ig*ndat] - scale * conv;\n\
x_out[jg*n + ig] = tmp[0];\n\
};\n\
";

#define oclErrorCheck(errorCode,message) if(errorCode!=CL_SUCCESS) { fprintf(stderr,"Error(%i) (%s: %s): %s\n", errorCode,__FILE__,__func__,message);exit(1);} 


cl_kernel magicfilter1d_kernel_l;
cl_kernel ana1d_kernel_l;
cl_kernel syn1d_kernel_l;
cl_kernel kinetic1d_kernel_l;
cl_kernel c_initialize_kernel_l;

cl_device_id oclGetFirstDev(cl_context cxGPUContext)
{
    size_t szParmDataBytes;
    cl_device_id* cdDevices;

    // get the list of GPU devices associated with context
    clGetContextInfo(cxGPUContext, CL_CONTEXT_DEVICES, 0, NULL, &szParmDataBytes);
    cdDevices = (cl_device_id*) malloc(szParmDataBytes);

    clGetContextInfo(cxGPUContext, CL_CONTEXT_DEVICES, szParmDataBytes, cdDevices, NULL);

    cl_device_id first = cdDevices[0];
    free(cdDevices);

    return first;
}

void FC_FUNC_(ocl_build_kernels,OCL_BUILD_KERNELS)(cl_context * context) {
    cl_int ciErrNum = CL_SUCCESS;
    cl_program magicfilter1dProgram = clCreateProgramWithSource(*context,1,(const char**) &magicfilter1d_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(magicfilter1dProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build magicfilter1d program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(magicfilter1dProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    magicfilter1d_kernel_l=clCreateKernel(magicfilter1dProgram,"magicfilter1dKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");

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

    cl_program c_initializeProgram = clCreateProgramWithSource(*context,1,(const char**) &c_initialize_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(c_initializeProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build c_initialize program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(c_initializeProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    c_initialize_kernel_l=clCreateKernel(c_initializeProgram,"c_initializeKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");

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
 
}



void FC_FUNC_(ocl_create_gpu_context,OCL_CREATE_GPU_CONTEXT)(cl_context * context) {
    cl_int ciErrNum = CL_SUCCESS;
    *context = clCreateContextFromType(0, CL_DEVICE_TYPE_GPU, NULL, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p\n",*context);
#endif
    oclErrorCheck(ciErrNum,"Failed to create GPU context!");
}

void FC_FUNC_(ocl_create_cpu_context,OCL_CREATE_CPU_CONTEXT)(cl_context * context) {
    cl_int ciErrNum = CL_SUCCESS;
    *context = clCreateContextFromType(0, CL_DEVICE_TYPE_CPU, NULL, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p\n",*context);
#endif
    oclErrorCheck(ciErrNum,"Failed to create CPU context!");
}

void FC_FUNC_(ocl_create_read_buffer,OCL_CREATE_READ_BUFFER)(cl_context *context, size_t *size, cl_mem *buff_ptr) {
    cl_int ciErrNum = CL_SUCCESS;
    *buff_ptr = clCreateBuffer( *context, CL_MEM_READ_ONLY, *size, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, memory address: %p, size: %d\n",*context,*buff_ptr,*size);
#endif
    oclErrorCheck(ciErrNum,"Failed to create read buffer!");
}

void FC_FUNC_(ocl_create_read_write_buffer,OCL_CREATE_READ_WRITE_BUFFER)(cl_context *context, size_t *size, cl_mem *buff_ptr) {
    cl_int ciErrNum = CL_SUCCESS;
    *buff_ptr = clCreateBuffer( *context, CL_MEM_READ_WRITE, *size, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, memory address: %p, size: %d\n",*context,*buff_ptr,*size);
#endif
    oclErrorCheck(ciErrNum,"Failed to create read_write buffer!");
}

void FC_FUNC_(ocl_create_read_buffer_and_copy,OCL_CREATE_READ_BUFFER_AND_COPY)(cl_context *context, size_t *size, void *host_ptr, cl_mem *buff_ptr) {
    cl_int ciErrNum = CL_SUCCESS;
    *buff_ptr = clCreateBuffer( *context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, *size, host_ptr, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create initialized read buffer!");
}

void FC_FUNC_(ocl_create_write_buffer,OCL_CREATE_WRITE_BUFFER)(cl_context *context, size_t *size, cl_mem *buff_ptr) {
    cl_int ciErrNum = CL_SUCCESS;
    *buff_ptr = clCreateBuffer( *context, CL_MEM_WRITE_ONLY, *size, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, memory address: %p, size: %d\n",*context,*buff_ptr,*size);
#endif
    oclErrorCheck(ciErrNum,"Failed to create write buffer!");
}

void FC_FUNC_(ocl_release_mem_object,OCL_RELEASE_MEM_OBJECT)(cl_mem *buff_ptr) {
    cl_int ciErrNum = clReleaseMemObject( *buff_ptr);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("memory address: %p\n",*buff_ptr);
#endif
    oclErrorCheck(ciErrNum,"Failed to release buffer!");
}

void FC_FUNC_(ocl_enqueue_read_buffer,OCL_ENQUEUE_READ_BUFFER)(cl_command_queue *command_queue, cl_mem *buffer, size_t *size, void *ptr){
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, memory address: %p, size: %d, target: %p\n",*command_queue,*buffer,*size, ptr);
#endif
    cl_int ciErrNum = clEnqueueReadBuffer( *command_queue, *buffer, CL_TRUE, 0, *size, ptr, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue read buffer!");
}

void FC_FUNC_(ocl_enqueue_write_buffer,OCL_ENQUEUE_WRITE_BUFFER)(cl_command_queue *command_queue, cl_mem *buffer, size_t *size,	const void *ptr){
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, memory address: %p, size: %d, source: %p\n",*command_queue,*buffer,*size, ptr);
#endif
    cl_int ciErrNum = clEnqueueWriteBuffer( *command_queue, *buffer, CL_TRUE, 0, *size, ptr, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue write buffer!");
}

void FC_FUNC_(ocl_create_command_queue,OCL_CREATE_COMMAND_QUEUE)(cl_command_queue *hCmdQueue, cl_context *context){
    size_t nContextDescriptorSize;
    cl_int ciErrNum;
    clGetContextInfo(*context, CL_CONTEXT_DEVICES, 0, 0, &nContextDescriptorSize);
    cl_device_id * aDevices = (cl_device_id *) malloc(nContextDescriptorSize);
    clGetContextInfo(*context, CL_CONTEXT_DEVICES, nContextDescriptorSize, aDevices, 0);
    // create a command queue for first device the context reported
    *hCmdQueue = clCreateCommandQueue(*context, aDevices[0], 0, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, command queue: %p\n",*context, *hCmdQueue);
#endif
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to create command queue!\n");
        exit(1);
    }
}


size_t shrRoundUp(size_t group_size, size_t global_size)
{
    size_t r = global_size % group_size;
    if(r == 0)
    {
        return global_size;
    } else
    {
        return global_size + group_size - r;
    }
}


void magicfilter1dKernelCheck(size_t n, size_t ndat, double *psi, double *out){
double filter[]={8.4334247333529341094733325815816e-7,
                -0.1290557201342060969516786758559028e-4,
                0.8762984476210559564689161894116397e-4,
                -0.30158038132690463167163703826169879e-3,
                0.174723713672993903449447812749852942e-2,
                -0.942047030201080385922711540948195075e-2,
                0.2373821463724942397566389712597274535e-1,
                0.612625895831207982195380597e-1,
                0.9940415697834003993178616713,
                -0.604895289196983516002834636e-1,
                -0.2103025160930381434955489412839065067e-1,
                0.1337263414854794752733423467013220997e-1,
                -0.344128144493493857280881509686821861e-2,
                0.49443227688689919192282259476750972e-3,
                -0.5185986881173432922848639136911487e-4,
                2.72734492911979659657715313017228e-6};
double *filt = filter + 8;
int i,j,l,k;
for(j = 0; j < ndat; j++) {
    for( i = 0; i < n; i++) {
        double tt = 0.0;
        for( l = -8; l <= 7; l++) {
            k = (i+l)%n;
/*            printf("%d %lf\n", k, psi[k*ndat + j]);*/
            tt = tt + psi[k*ndat + j] * filt[l];
        }
/*        printf("%lf\n",tt);*/
        out[i + j*n] = tt;
/*        return;*/
    }

}
}

void FC_FUNC_(magicfilter1d_check,MAGICFILTER1D_CHECK)(size_t *n,size_t *ndat,void *psi,void *out){

	magicfilter1dKernelCheck(*n, *ndat, psi, out);

}

#define BLOCK_SIZE_I 64
#define BLOCK_SIZE_J 4

void FC_FUNC_(magicfilter1d_l,MAGICFILTER1D_L)(cl_command_queue *command_queue, size_t *n,size_t *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %d, dimension dat: %d, psi: %p, out: %p\n",*command_queue, *n, *ndat, *psi, *out);
#endif
    int FILTER_WIDTH = 16;
    if(*n<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=256/FILTER_WIDTH;
    while (*n > block_size_i >= 1 && block_size_j > 4)
	{ block_size_i *= 2; block_size_j /= 2;}

    cl_uint i = 0;
    clSetKernelArg(magicfilter1d_kernel_l, i++,sizeof(*n), (void*)n);
    clSetKernelArg(magicfilter1d_kernel_l, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(magicfilter1d_kernel_l, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(magicfilter1d_kernel_l, i++,sizeof(*out), (void*)out);
    clSetKernelArg(magicfilter1d_kernel_l, i++,sizeof(float)*block_size_j*(block_size_i+FILTER_WIDTH), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_kernel_l, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %d, %d}\n",globalWorkSize[0],globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %d, %d}\n",localWorkSize[0],localWorkSize[1]);
        exit(1);
    }   
}

void FC_FUNC_(ana1d_l,ANA1D_L)(cl_command_queue *command_queue, size_t *n,size_t *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %d, dimension dat: %d, psi: %p, out: %p\n",*command_queue, *n, *ndat, *psi, *out);
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
        fprintf(stderr,"globalWorkSize = { %d, %d}\n",globalWorkSize[0],globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %d, %d}\n",localWorkSize[0],localWorkSize[1]);
        exit(1);
    }

}

void FC_FUNC_(syn1d_l,SYN1D_L)(cl_command_queue *command_queue, size_t *n,size_t *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %d, dimension dat: %d, psi: %p, out: %p\n",*command_queue, *n, *ndat, *psi, *out);
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
        fprintf(stderr,"globalWorkSize = { %d, %d}\n",globalWorkSize[0],globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %d, %d}\n",localWorkSize[0],localWorkSize[1]);
        exit(1);
    }  
}

void FC_FUNC_(kinetic1d_l,KINETIC1D_L)(cl_command_queue *command_queue, size_t *n,size_t *ndat, float *h, float*c, cl_mem *x, cl_mem *y, cl_mem *workx, cl_mem *worky,float *ekin){
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %d, dimension dat: %d, h: %f, c: %f, x: %p, workx: %p, y: %p, worky: %p\n",*command_queue, *n, *ndat, *h, *c, *x, *workx, *y, *worky);
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
        fprintf(stderr,"globalWorkSize = { %d, %d}\n",globalWorkSize[0],globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %d, %d}\n",localWorkSize[0],localWorkSize[1]);
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
        fprintf(stderr,"globalWorkSize = { %d, %d}\n",globalWorkSize[0],globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %d, %d}\n",localWorkSize[0],localWorkSize[1]);
        exit(1);
    }
}

void FC_FUNC_(ocl_finish,OCL_FINISH)(cl_command_queue *command_queue){
    cl_int ciErrNum;
    ciErrNum = clFinish(*command_queue);
    oclErrorCheck(ciErrNum,"Failed to finish!");
}

void FC_FUNC_(ocl_enqueue_barrier,OCL_ENQUEUE_BARRIER)(cl_command_queue *command_queue){
    cl_int ciErrNum;
    ciErrNum = clEnqueueBarrier(*command_queue);
    oclErrorCheck(ciErrNum,"Failed to enqueue barrier!");
}
