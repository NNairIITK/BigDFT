#include "Kinetic.h"
#include "Kinetic_k.h"
#include "OpenCL_wrappers.h"

char * kinetic1d_program="\
#define FILTER_WIDTH 32\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
__kernel void kinetic1dKernel_d(size_t n, size_t ndat, double scale, __global const double * x_in, __global double * x_out, __global const double * y_in, __global double * y_out, __local double * tmp, __local double * tmp2) {\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1);\n\
ptrdiff_t igt = get_group_id(0);\n\
ptrdiff_t jgt = get_group_id(1);\n\
size_t jb;\n\
const size_t j2t = 2*j2 + i2/16;\n\
const size_t i2t = i2 - 16 * (i2 / 16);\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
igt = ig - i2 + j2t;\n\
jgt = jg - j2 + i2t;\n\
tmp2[i2t * (32 + 1) + j2t] = y_in[jgt+igt*ndat];\n\
igt -= FILTER_WIDTH/2;\n\
if ( igt < 0 ) \n\
  tmp[i2t * (2 * FILTER_WIDTH + 1) + j2t] = x_in[jgt + ( n + igt ) * ndat];\n\
else \n\
  tmp[i2t * (2 * FILTER_WIDTH + 1) + j2t] = x_in[jgt + igt * ndat];\n\
igt += FILTER_WIDTH;\n\
if ( igt >= n ) \n\
  tmp[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = x_in[jgt + ( igt - n ) * ndat];\n\
else\n\
  tmp[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = x_in[jgt +  igt * ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
__local double * tmp_o = tmp + j2*(2*32 + 1) + FILTER_WIDTH/2+i2;\n\
double conv = 0.0;\n\
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
//y_out[jg*n + ig] = y_in[jg + ig*ndat] - scale * conv;\n\
y_out[jg*n + ig] = tmp2[j2*(32+1) + i2] + scale * conv;\n\
x_out[jg*n + ig] = tmp_o[0];\n\
};\n\
\n\
";

inline void kinetic_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, double *scale, cl_mem *x_in, cl_mem *x_out, cl_mem *y_in, cl_mem *y_out) {
  int FILTER_WIDTH = 32;
  cl_int ciErrNum;
  assert(*n>=FILTER_WIDTH);
  size_t block_size_i=FILTER_WIDTH;
  size_t block_size_j=16;
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
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue kinetic kernel!");
}

char * kinetic_k1d_program="\
#define FILTER_WIDTH 32\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
__kernel void kinetic_k1dKernel_d(size_t n, size_t ndat, double scale_1, double scale_2, __global const double * x_in, __global double * x, __global const double * y_in, __global double * y, __local double * tmp, __local double * tmp_y ) {\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
const size_t i2 = get_local_id(0);\n\
const size_t j2 = get_local_id(1);\n\
ptrdiff_t igt = get_group_id(0);\n\
ptrdiff_t jgt = get_group_id(1);\n\
const size_t j2t = 4*j2 + i2/8;\n\
const size_t i2t = i2 - 8 * (i2 / 8);\n\
//if data are ill dimentioned last block recomputes part of the data\n\
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
igt = ig - i2 + j2t;\n\
jgt = jg - j2 + i2t;\n\
tmp_y[i2t * (2 * FILTER_WIDTH + 1) + 2 * j2t] = y_in[2*jgt+2*igt*ndat];\n\
tmp_y[i2t * (2 * FILTER_WIDTH + 1) + 2 * j2t + 1] = y_in[2*jgt+2*igt*ndat+1];\n\
igt -= FILTER_WIDTH/2;\n\
if ( igt < 0 ) {\n\
  tmp[i2t * (4 * FILTER_WIDTH + 1) + 2 * j2t] = x_in[2*jgt + 2 * ( n + igt ) * ndat];\n\
  tmp[i2t * (4 * FILTER_WIDTH + 1) + 2 * j2t + 1] = x_in[2*jgt + 2 * ( n + igt ) * ndat + 1];\n\
} else {\n\
  tmp[i2t * (4 * FILTER_WIDTH + 1) + 2 * j2t] = x_in[2*jgt + 2 * igt * ndat];\n\
  tmp[i2t * (4 * FILTER_WIDTH + 1) + 2 * j2t + 1] = x_in[2*jgt + 2 * igt * ndat+1];\n\
}\n\
igt += FILTER_WIDTH;\n\
if ( igt >= n ) {\n\
  tmp[i2t * (4 * FILTER_WIDTH + 1) + 2 * j2t + 2 * FILTER_WIDTH] = x_in[2*jgt + 2 * ( igt - n ) * ndat];\n\
  tmp[i2t * (4 * FILTER_WIDTH + 1) + 2 * j2t + 2 * FILTER_WIDTH + 1] = x_in[2*jgt + 2 * ( igt - n ) * ndat + 1];\n\
} else {\n\
  tmp[i2t * (4 * FILTER_WIDTH + 1) + 2 * j2t + 2 * FILTER_WIDTH] = x_in[2*jgt +  2 * igt * ndat];\n\
  tmp[i2t * (4 * FILTER_WIDTH + 1) + 2 * j2t + 2 * FILTER_WIDTH + 1] = x_in[2*jgt +  2 * igt * ndat + 1];\n\
}\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
__local double * tmp_o = tmp + j2*(4 * FILTER_WIDTH + 1) + FILTER_WIDTH + 2*i2;\n\
double tt1_1=0;\n\
tt1_1 += (tmp_o[2*14] + tmp_o[2*-14]) * -6.924474940639200152025730585882e-18;\n\
tt1_1 += (tmp_o[2*13] + tmp_o[2*-13]) *  2.70800493626319438269856689037647576e-13;\n\
tt1_1 += (tmp_o[2*12] + tmp_o[2*-12]) * -5.813879830282540547959250667e-11;\n\
tt1_1 += (tmp_o[2*11] + tmp_o[2*-11]) * -1.05857055496741470373494132287e-8;\n\
tt1_1 += (tmp_o[2*10] + tmp_o[2*-10]) * -3.7230763047369275848791496973044e-7;\n\
tt1_1 += (tmp_o[2* 9] + tmp_o[2* -9]) *  2.0904234952920365957922889447361e-6;\n\
tt1_1 += (tmp_o[2* 8] + tmp_o[2* -8]) * -0.2398228524507599670405555359023135e-4;\n\
tt1_1 += (tmp_o[2* 7] + tmp_o[2* -7]) *  0.45167920287502235349480037639758496e-3;\n\
tt1_1 += (tmp_o[2* 6] + tmp_o[2* -6]) * -0.409765689342633823899327051188315485e-2;\n\
tt1_1 += (tmp_o[2* 5] + tmp_o[2* -5]) *  0.02207029188482255523789911295638968409e0;\n\
tt1_1 += (tmp_o[2* 4] + tmp_o[2* -4]) * -0.0822663999742123340987663521e0;\n\
tt1_1 += (tmp_o[2* 3] + tmp_o[2* -3]) *  0.2371780582153805636239247476e0;\n\
tt1_1 += (tmp_o[2* 2] + tmp_o[2* -2]) * -0.6156141465570069496314853949e0;\n\
tt1_1 += (tmp_o[2* 1] + tmp_o[2* -1]) *  2.2191465938911163898794546405e0;\n\
tt1_1 +=  tmp_o[2* 0]                 * -3.5536922899131901941296809374e0;\n\
double tt1_2=0;\n\
tt1_2 += (tmp_o[2*14+1] + tmp_o[2*-14+1]) * -6.924474940639200152025730585882e-18;\n\
tt1_2 += (tmp_o[2*13+1] + tmp_o[2*-13+1]) *  2.70800493626319438269856689037647576e-13;\n\
tt1_2 += (tmp_o[2*12+1] + tmp_o[2*-12+1]) * -5.813879830282540547959250667e-11;\n\
tt1_2 += (tmp_o[2*11+1] + tmp_o[2*-11+1]) * -1.05857055496741470373494132287e-8;\n\
tt1_2 += (tmp_o[2*10+1] + tmp_o[2*-10+1]) * -3.7230763047369275848791496973044e-7;\n\
tt1_2 += (tmp_o[2* 9+1] + tmp_o[2* -9+1]) *  2.0904234952920365957922889447361e-6;\n\
tt1_2 += (tmp_o[2* 8+1] + tmp_o[2* -8+1]) * -0.2398228524507599670405555359023135e-4;\n\
tt1_2 += (tmp_o[2* 7+1] + tmp_o[2* -7+1]) *  0.45167920287502235349480037639758496e-3;\n\
tt1_2 += (tmp_o[2* 6+1] + tmp_o[2* -6+1]) * -0.409765689342633823899327051188315485e-2;\n\
tt1_2 += (tmp_o[2* 5+1] + tmp_o[2* -5+1]) *  0.02207029188482255523789911295638968409e0;\n\
tt1_2 += (tmp_o[2* 4+1] + tmp_o[2* -4+1]) * -0.0822663999742123340987663521e0;\n\
tt1_2 += (tmp_o[2* 3+1] + tmp_o[2* -3+1]) *  0.2371780582153805636239247476e0;\n\
tt1_2 += (tmp_o[2* 2+1] + tmp_o[2* -2+1]) * -0.6156141465570069496314853949e0;\n\
tt1_2 += (tmp_o[2* 1+1] + tmp_o[2* -1+1]) *  2.2191465938911163898794546405e0;\n\
tt1_2 +=  tmp_o[2* 0+1]                   * -3.5536922899131901941296809374e0;\n\
double tt2_1=0;\n\
tt2_1 += (tmp_o[2*14] - tmp_o[2*-14]) * -1.585464751677102510097179e-19;\n\
tt2_1 += (tmp_o[2*13] - tmp_o[2*-13]) *  1.240078536096648534547439e-14;\n\
tt2_1 += (tmp_o[2*12] - tmp_o[2*-12]) * -7.252206916665149851135592e-13;\n\
tt2_1 += (tmp_o[2*11] - tmp_o[2*-11]) * -9.697184925637300947553069e-10;\n\
tt2_1 += (tmp_o[2*10] - tmp_o[2*-10]) * -7.207948238588481597101904e-8;\n\
tt2_1 += (tmp_o[2* 9] - tmp_o[2* -9]) *  3.993810456408053712133667e-8;\n\
tt2_1 += (tmp_o[2* 8] - tmp_o[2* -8]) *  2.451992111053665419191564e-7;\n\
tt2_1 += (tmp_o[2* 7] - tmp_o[2* -7]) *  0.00007667706908380351933901775e0;\n\
tt2_1 += (tmp_o[2* 6] - tmp_o[2* -6]) * -0.001031530213375445369097965e0;\n\
tt2_1 += (tmp_o[2* 5] - tmp_o[2* -5]) *  0.006958379116450707495020408e0;\n\
tt2_1 += (tmp_o[2* 4] - tmp_o[2* -4]) * -0.03129014783948023634381564e0;\n\
tt2_1 += (tmp_o[2* 3] - tmp_o[2* -3]) *  0.1063640682894442760934532e0;\n\
tt2_1 += (tmp_o[2* 2] - tmp_o[2* -2]) * -0.3032593514765938346887962e0;\n\
tt2_1 += (tmp_o[2* 1] - tmp_o[2* -1]) *  0.8834460460908270942785856e0;\n\
double tt2_2=0;\n\
tt2_2 += (tmp_o[2*14+1] - tmp_o[2*-14+1]) * -1.585464751677102510097179e-19;\n\
tt2_2 += (tmp_o[2*13+1] - tmp_o[2*-13+1]) *  1.240078536096648534547439e-14;\n\
tt2_2 += (tmp_o[2*12+1] - tmp_o[2*-12+1]) * -7.252206916665149851135592e-13;\n\
tt2_2 += (tmp_o[2*11+1] - tmp_o[2*-11+1]) * -9.697184925637300947553069e-10;\n\
tt2_2 += (tmp_o[2*10+1] - tmp_o[2*-10+1]) * -7.207948238588481597101904e-8;\n\
tt2_2 += (tmp_o[2* 9+1] - tmp_o[2* -9+1]) *  3.993810456408053712133667e-8;\n\
tt2_2 += (tmp_o[2* 8+1] - tmp_o[2* -8+1]) *  2.451992111053665419191564e-7;\n\
tt2_2 += (tmp_o[2* 7+1] - tmp_o[2* -7+1]) *  0.00007667706908380351933901775e0;\n\
tt2_2 += (tmp_o[2* 6+1] - tmp_o[2* -6+1]) * -0.001031530213375445369097965e0;\n\
tt2_2 += (tmp_o[2* 5+1] - tmp_o[2* -5+1]) *  0.006958379116450707495020408e0;\n\
tt2_2 += (tmp_o[2* 4+1] - tmp_o[2* -4+1]) * -0.03129014783948023634381564e0;\n\
tt2_2 += (tmp_o[2* 3+1] - tmp_o[2* -3+1]) *  0.1063640682894442760934532e0;\n\
tt2_2 += (tmp_o[2* 2+1] - tmp_o[2* -2+1]) * -0.3032593514765938346887962e0;\n\
tt2_2 += (tmp_o[2* 1+1] - tmp_o[2* -1+1]) *  0.8834460460908270942785856e0;\n\
y[jg*2*n + 2*ig] = tmp_y[j2*(2*FILTER_WIDTH+1) + 2*i2] + tt1_1 * scale_1 - tt2_2 * scale_2;\n\
y[jg*2*n + 2*ig + 1] = tmp_y[j2*(2*FILTER_WIDTH+1) + 2*i2 + 1] + tt1_2 * scale_1 + tt2_1 * scale_2;\n\
x[jg*2*n + 2*ig] = tmp_o[0];\n\
x[jg*2*n + 2*ig + 1] = tmp_o[1];\n\
}\n\
";

inline void kinetic_k_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, double *scale1, double *scale2, cl_mem *x_in, cl_mem *x_out, cl_mem *y_in, cl_mem *y_out) {
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
  clSetKernelArg(kernel, i++,sizeof(double)*block_size_j*2*(block_size_i+FILTER_WIDTH+1), NULL);
  clSetKernelArg(kernel, i++,sizeof(double)*block_size_j*2*(block_size_i+1), NULL);
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue kinetic_k kernel!");
} 


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
