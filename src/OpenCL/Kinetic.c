#include "Initialize.h"
#include "Kinetic.h"
#include "Kinetic_k.h"
#include "OpenCL_wrappers.h"
#include "Kinetic_Generator.h"

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

char * kinetic_k_program="\
#ifdef cl_khr_fp64\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#elif defined (cl_amd_fp64)\n\
#pragma OPENCL EXTENSION cl_amd_fp64: enable \n\
#endif\n\
#define FILTER_WIDTH 32\n\
__kernel void kinetic_k1dKernel_d(uint n, uint ndat, double scale_1, double scale_2, __global const double * x_in, __global double * x, __global const double * y_in, __global double * y, __local double * tmp, __local double * tmp_y ) {\n\
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
__kernel void kinetic_k1dKernel_d_2(uint n, uint ndat, double scale_1, double scale_2, __global const double * x_in_r, __global const double * x_in_i, __global double * x_r, __global double * x_i, __global const double * y_in_r, __global const double * y_in_i, __global double * y_r, __global double * y_i, __local double * tmp_r, __local double * tmp_i, __local double * tmp_y_r, __local double * tmp_y_i ) {\n\
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
tmp_y_r[i2t * (FILTER_WIDTH + 1) + j2t] = y_in_r[jgt+igt*ndat];\n\
tmp_y_i[i2t * (FILTER_WIDTH + 1) + j2t] = y_in_i[jgt+igt*ndat];\n\
igt -= FILTER_WIDTH/2;\n\
if ( igt < 0 ) {\n\
  tmp_r[i2t * (2 * FILTER_WIDTH + 1) + j2t] = x_in_r[jgt + ( n + igt ) * ndat];\n\
  tmp_i[i2t * (2 * FILTER_WIDTH + 1) + j2t] = x_in_i[jgt + ( n + igt ) * ndat];\n\
} else {\n\
  tmp_r[i2t * (2 * FILTER_WIDTH + 1) + j2t] = x_in_r[jgt + igt * ndat];\n\
  tmp_i[i2t * (2 * FILTER_WIDTH + 1) + j2t] = x_in_i[jgt + igt * ndat];\n\
}\n\
igt += FILTER_WIDTH;\n\
if ( igt >= n ) {\n\
  tmp_r[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = x_in_r[jgt + ( igt - n ) * ndat];\n\
  tmp_i[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = x_in_i[jgt + ( igt - n ) * ndat];\n\
} else {\n\
  tmp_r[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = x_in_r[jgt +  igt * ndat];\n\
  tmp_i[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = x_in_i[jgt +  igt * ndat];\n\
}\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
__local double * tmp_o_r = tmp_r + j2*(2 * FILTER_WIDTH + 1) + FILTER_WIDTH/2 + i2;\n\
__local double * tmp_o_i = tmp_i + j2*(2 * FILTER_WIDTH + 1) + FILTER_WIDTH/2 + i2;\n\
double tt1_1=0;\n\
tt1_1 += (tmp_o_r[14] + tmp_o_r[-14]) * -6.924474940639200152025730585882e-18;\n\
tt1_1 += (tmp_o_r[13] + tmp_o_r[-13]) *  2.70800493626319438269856689037647576e-13;\n\
tt1_1 += (tmp_o_r[12] + tmp_o_r[-12]) * -5.813879830282540547959250667e-11;\n\
tt1_1 += (tmp_o_r[11] + tmp_o_r[-11]) * -1.05857055496741470373494132287e-8;\n\
tt1_1 += (tmp_o_r[10] + tmp_o_r[-10]) * -3.7230763047369275848791496973044e-7;\n\
tt1_1 += (tmp_o_r[ 9] + tmp_o_r[ -9]) *  2.0904234952920365957922889447361e-6;\n\
tt1_1 += (tmp_o_r[ 8] + tmp_o_r[ -8]) * -0.2398228524507599670405555359023135e-4;\n\
tt1_1 += (tmp_o_r[ 7] + tmp_o_r[ -7]) *  0.45167920287502235349480037639758496e-3;\n\
tt1_1 += (tmp_o_r[ 6] + tmp_o_r[ -6]) * -0.409765689342633823899327051188315485e-2;\n\
tt1_1 += (tmp_o_r[ 5] + tmp_o_r[ -5]) *  0.02207029188482255523789911295638968409e0;\n\
tt1_1 += (tmp_o_r[ 4] + tmp_o_r[ -4]) * -0.0822663999742123340987663521e0;\n\
tt1_1 += (tmp_o_r[ 3] + tmp_o_r[ -3]) *  0.2371780582153805636239247476e0;\n\
tt1_1 += (tmp_o_r[ 2] + tmp_o_r[ -2]) * -0.6156141465570069496314853949e0;\n\
tt1_1 += (tmp_o_r[ 1] + tmp_o_r[ -1]) *  2.2191465938911163898794546405e0;\n\
tt1_1 +=  tmp_o_r[ 0]                 * -3.5536922899131901941296809374e0;\n\
double tt1_2=0;\n\
tt1_2 += (tmp_o_i[14] + tmp_o_i[-14]) * -6.924474940639200152025730585882e-18;\n\
tt1_2 += (tmp_o_i[13] + tmp_o_i[-13]) *  2.70800493626319438269856689037647576e-13;\n\
tt1_2 += (tmp_o_i[12] + tmp_o_i[-12]) * -5.813879830282540547959250667e-11;\n\
tt1_2 += (tmp_o_i[11] + tmp_o_i[-11]) * -1.05857055496741470373494132287e-8;\n\
tt1_2 += (tmp_o_i[10] + tmp_o_i[-10]) * -3.7230763047369275848791496973044e-7;\n\
tt1_2 += (tmp_o_i[ 9] + tmp_o_i[ -9]) *  2.0904234952920365957922889447361e-6;\n\
tt1_2 += (tmp_o_i[ 8] + tmp_o_i[ -8]) * -0.2398228524507599670405555359023135e-4;\n\
tt1_2 += (tmp_o_i[ 7] + tmp_o_i[ -7]) *  0.45167920287502235349480037639758496e-3;\n\
tt1_2 += (tmp_o_i[ 6] + tmp_o_i[ -6]) * -0.409765689342633823899327051188315485e-2;\n\
tt1_2 += (tmp_o_i[ 5] + tmp_o_i[ -5]) *  0.02207029188482255523789911295638968409e0;\n\
tt1_2 += (tmp_o_i[ 4] + tmp_o_i[ -4]) * -0.0822663999742123340987663521e0;\n\
tt1_2 += (tmp_o_i[ 3] + tmp_o_i[ -3]) *  0.2371780582153805636239247476e0;\n\
tt1_2 += (tmp_o_i[ 2] + tmp_o_i[ -2]) * -0.6156141465570069496314853949e0;\n\
tt1_2 += (tmp_o_i[ 1] + tmp_o_i[ -1]) *  2.2191465938911163898794546405e0;\n\
tt1_2 +=  tmp_o_i[ 0]                   * -3.5536922899131901941296809374e0;\n\
double tt2_1=0;\n\
tt2_1 += (tmp_o_r[14] - tmp_o_r[-14]) * -1.585464751677102510097179e-19;\n\
tt2_1 += (tmp_o_r[13] - tmp_o_r[-13]) *  1.240078536096648534547439e-14;\n\
tt2_1 += (tmp_o_r[12] - tmp_o_r[-12]) * -7.252206916665149851135592e-13;\n\
tt2_1 += (tmp_o_r[11] - tmp_o_r[-11]) * -9.697184925637300947553069e-10;\n\
tt2_1 += (tmp_o_r[10] - tmp_o_r[-10]) * -7.207948238588481597101904e-8;\n\
tt2_1 += (tmp_o_r[ 9] - tmp_o_r[ -9]) *  3.993810456408053712133667e-8;\n\
tt2_1 += (tmp_o_r[ 8] - tmp_o_r[ -8]) *  2.451992111053665419191564e-7;\n\
tt2_1 += (tmp_o_r[ 7] - tmp_o_r[ -7]) *  0.00007667706908380351933901775e0;\n\
tt2_1 += (tmp_o_r[ 6] - tmp_o_r[ -6]) * -0.001031530213375445369097965e0;\n\
tt2_1 += (tmp_o_r[ 5] - tmp_o_r[ -5]) *  0.006958379116450707495020408e0;\n\
tt2_1 += (tmp_o_r[ 4] - tmp_o_r[ -4]) * -0.03129014783948023634381564e0;\n\
tt2_1 += (tmp_o_r[ 3] - tmp_o_r[ -3]) *  0.1063640682894442760934532e0;\n\
tt2_1 += (tmp_o_r[ 2] - tmp_o_r[ -2]) * -0.3032593514765938346887962e0;\n\
tt2_1 += (tmp_o_r[ 1] - tmp_o_r[ -1]) *  0.8834460460908270942785856e0;\n\
double tt2_2=0;\n\
tt2_2 += (tmp_o_i[14] - tmp_o_i[-14]) * -1.585464751677102510097179e-19;\n\
tt2_2 += (tmp_o_i[13] - tmp_o_i[-13]) *  1.240078536096648534547439e-14;\n\
tt2_2 += (tmp_o_i[12] - tmp_o_i[-12]) * -7.252206916665149851135592e-13;\n\
tt2_2 += (tmp_o_i[11] - tmp_o_i[-11]) * -9.697184925637300947553069e-10;\n\
tt2_2 += (tmp_o_i[10] - tmp_o_i[-10]) * -7.207948238588481597101904e-8;\n\
tt2_2 += (tmp_o_i[ 9] - tmp_o_i[ -9]) *  3.993810456408053712133667e-8;\n\
tt2_2 += (tmp_o_i[ 8] - tmp_o_i[ -8]) *  2.451992111053665419191564e-7;\n\
tt2_2 += (tmp_o_i[ 7] - tmp_o_i[ -7]) *  0.00007667706908380351933901775e0;\n\
tt2_2 += (tmp_o_i[ 6] - tmp_o_i[ -6]) * -0.001031530213375445369097965e0;\n\
tt2_2 += (tmp_o_i[ 5] - tmp_o_i[ -5]) *  0.006958379116450707495020408e0;\n\
tt2_2 += (tmp_o_i[ 4] - tmp_o_i[ -4]) * -0.03129014783948023634381564e0;\n\
tt2_2 += (tmp_o_i[ 3] - tmp_o_i[ -3]) *  0.1063640682894442760934532e0;\n\
tt2_2 += (tmp_o_i[ 2] - tmp_o_i[ -2]) * -0.3032593514765938346887962e0;\n\
tt2_2 += (tmp_o_i[ 1] - tmp_o_i[ -1]) *  0.8834460460908270942785856e0;\n\
y_r[jg*n + ig] = tmp_y_r[j2*(FILTER_WIDTH+1) + i2] + tt1_1 * scale_1 - tt2_2 * scale_2;\n\
y_i[jg*n + ig] = tmp_y_i[j2*(FILTER_WIDTH+1) + i2] + tt1_2 * scale_1 + tt2_1 * scale_2;\n\
x_r[jg*n + ig] = tmp_o_r[0];\n\
x_i[jg*n + ig] = tmp_o_i[0];\n\
}\n\
__kernel void kinetic_k1d_fKernel_d_2(uint n, uint ndat, double scale_1, double scale_2, __global const double * x_in_r, __global const double * x_in_i, __global double * x_r, __global double * x_i, __global const double * y_in_r, __global const double * y_in_i, __global double * y_r, __global double * y_i, __local double * tmp_r, __local double * tmp_i, __local double * tmp_y_r, __local double * tmp_y_i ) {\n\
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
tmp_y_r[i2t * (FILTER_WIDTH + 1) + j2t] = y_in_r[jgt+igt*ndat];\n\
tmp_y_i[i2t * (FILTER_WIDTH + 1) + j2t] = y_in_i[jgt+igt*ndat];\n\
igt -= FILTER_WIDTH/2;\n\
if ( igt < 0 ) {\n\
  tmp_r[i2t * (2 * FILTER_WIDTH + 1) + j2t] = 0.0;\n\
  tmp_i[i2t * (2 * FILTER_WIDTH + 1) + j2t] = 0.0;\n\
} else {\n\
  tmp_r[i2t * (2 * FILTER_WIDTH + 1) + j2t] = x_in_r[jgt + igt * ndat];\n\
  tmp_i[i2t * (2 * FILTER_WIDTH + 1) + j2t] = x_in_i[jgt + igt * ndat];\n\
}\n\
igt += FILTER_WIDTH;\n\
if ( igt >= n ) {\n\
  tmp_r[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = 0.0;\n\
  tmp_i[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = 0.0;\n\
} else {\n\
  tmp_r[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = x_in_r[jgt +  igt * ndat];\n\
  tmp_i[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = x_in_i[jgt +  igt * ndat];\n\
}\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
__local double * tmp_o_r = tmp_r + j2*(2 * FILTER_WIDTH + 1) + FILTER_WIDTH/2 + i2;\n\
__local double * tmp_o_i = tmp_i + j2*(2 * FILTER_WIDTH + 1) + FILTER_WIDTH/2 + i2;\n\
double tt1_1=0;\n\
tt1_1 += (tmp_o_r[14] + tmp_o_r[-14]) * -6.924474940639200152025730585882e-18;\n\
tt1_1 += (tmp_o_r[13] + tmp_o_r[-13]) *  2.70800493626319438269856689037647576e-13;\n\
tt1_1 += (tmp_o_r[12] + tmp_o_r[-12]) * -5.813879830282540547959250667e-11;\n\
tt1_1 += (tmp_o_r[11] + tmp_o_r[-11]) * -1.05857055496741470373494132287e-8;\n\
tt1_1 += (tmp_o_r[10] + tmp_o_r[-10]) * -3.7230763047369275848791496973044e-7;\n\
tt1_1 += (tmp_o_r[ 9] + tmp_o_r[ -9]) *  2.0904234952920365957922889447361e-6;\n\
tt1_1 += (tmp_o_r[ 8] + tmp_o_r[ -8]) * -0.2398228524507599670405555359023135e-4;\n\
tt1_1 += (tmp_o_r[ 7] + tmp_o_r[ -7]) *  0.45167920287502235349480037639758496e-3;\n\
tt1_1 += (tmp_o_r[ 6] + tmp_o_r[ -6]) * -0.409765689342633823899327051188315485e-2;\n\
tt1_1 += (tmp_o_r[ 5] + tmp_o_r[ -5]) *  0.02207029188482255523789911295638968409e0;\n\
tt1_1 += (tmp_o_r[ 4] + tmp_o_r[ -4]) * -0.0822663999742123340987663521e0;\n\
tt1_1 += (tmp_o_r[ 3] + tmp_o_r[ -3]) *  0.2371780582153805636239247476e0;\n\
tt1_1 += (tmp_o_r[ 2] + tmp_o_r[ -2]) * -0.6156141465570069496314853949e0;\n\
tt1_1 += (tmp_o_r[ 1] + tmp_o_r[ -1]) *  2.2191465938911163898794546405e0;\n\
tt1_1 +=  tmp_o_r[ 0]                 * -3.5536922899131901941296809374e0;\n\
double tt1_2=0;\n\
tt1_2 += (tmp_o_i[14] + tmp_o_i[-14]) * -6.924474940639200152025730585882e-18;\n\
tt1_2 += (tmp_o_i[13] + tmp_o_i[-13]) *  2.70800493626319438269856689037647576e-13;\n\
tt1_2 += (tmp_o_i[12] + tmp_o_i[-12]) * -5.813879830282540547959250667e-11;\n\
tt1_2 += (tmp_o_i[11] + tmp_o_i[-11]) * -1.05857055496741470373494132287e-8;\n\
tt1_2 += (tmp_o_i[10] + tmp_o_i[-10]) * -3.7230763047369275848791496973044e-7;\n\
tt1_2 += (tmp_o_i[ 9] + tmp_o_i[ -9]) *  2.0904234952920365957922889447361e-6;\n\
tt1_2 += (tmp_o_i[ 8] + tmp_o_i[ -8]) * -0.2398228524507599670405555359023135e-4;\n\
tt1_2 += (tmp_o_i[ 7] + tmp_o_i[ -7]) *  0.45167920287502235349480037639758496e-3;\n\
tt1_2 += (tmp_o_i[ 6] + tmp_o_i[ -6]) * -0.409765689342633823899327051188315485e-2;\n\
tt1_2 += (tmp_o_i[ 5] + tmp_o_i[ -5]) *  0.02207029188482255523789911295638968409e0;\n\
tt1_2 += (tmp_o_i[ 4] + tmp_o_i[ -4]) * -0.0822663999742123340987663521e0;\n\
tt1_2 += (tmp_o_i[ 3] + tmp_o_i[ -3]) *  0.2371780582153805636239247476e0;\n\
tt1_2 += (tmp_o_i[ 2] + tmp_o_i[ -2]) * -0.6156141465570069496314853949e0;\n\
tt1_2 += (tmp_o_i[ 1] + tmp_o_i[ -1]) *  2.2191465938911163898794546405e0;\n\
tt1_2 +=  tmp_o_i[ 0]                   * -3.5536922899131901941296809374e0;\n\
double tt2_1=0;\n\
tt2_1 += (tmp_o_r[14] - tmp_o_r[-14]) * -1.585464751677102510097179e-19;\n\
tt2_1 += (tmp_o_r[13] - tmp_o_r[-13]) *  1.240078536096648534547439e-14;\n\
tt2_1 += (tmp_o_r[12] - tmp_o_r[-12]) * -7.252206916665149851135592e-13;\n\
tt2_1 += (tmp_o_r[11] - tmp_o_r[-11]) * -9.697184925637300947553069e-10;\n\
tt2_1 += (tmp_o_r[10] - tmp_o_r[-10]) * -7.207948238588481597101904e-8;\n\
tt2_1 += (tmp_o_r[ 9] - tmp_o_r[ -9]) *  3.993810456408053712133667e-8;\n\
tt2_1 += (tmp_o_r[ 8] - tmp_o_r[ -8]) *  2.451992111053665419191564e-7;\n\
tt2_1 += (tmp_o_r[ 7] - tmp_o_r[ -7]) *  0.00007667706908380351933901775e0;\n\
tt2_1 += (tmp_o_r[ 6] - tmp_o_r[ -6]) * -0.001031530213375445369097965e0;\n\
tt2_1 += (tmp_o_r[ 5] - tmp_o_r[ -5]) *  0.006958379116450707495020408e0;\n\
tt2_1 += (tmp_o_r[ 4] - tmp_o_r[ -4]) * -0.03129014783948023634381564e0;\n\
tt2_1 += (tmp_o_r[ 3] - tmp_o_r[ -3]) *  0.1063640682894442760934532e0;\n\
tt2_1 += (tmp_o_r[ 2] - tmp_o_r[ -2]) * -0.3032593514765938346887962e0;\n\
tt2_1 += (tmp_o_r[ 1] - tmp_o_r[ -1]) *  0.8834460460908270942785856e0;\n\
double tt2_2=0;\n\
tt2_2 += (tmp_o_i[14] - tmp_o_i[-14]) * -1.585464751677102510097179e-19;\n\
tt2_2 += (tmp_o_i[13] - tmp_o_i[-13]) *  1.240078536096648534547439e-14;\n\
tt2_2 += (tmp_o_i[12] - tmp_o_i[-12]) * -7.252206916665149851135592e-13;\n\
tt2_2 += (tmp_o_i[11] - tmp_o_i[-11]) * -9.697184925637300947553069e-10;\n\
tt2_2 += (tmp_o_i[10] - tmp_o_i[-10]) * -7.207948238588481597101904e-8;\n\
tt2_2 += (tmp_o_i[ 9] - tmp_o_i[ -9]) *  3.993810456408053712133667e-8;\n\
tt2_2 += (tmp_o_i[ 8] - tmp_o_i[ -8]) *  2.451992111053665419191564e-7;\n\
tt2_2 += (tmp_o_i[ 7] - tmp_o_i[ -7]) *  0.00007667706908380351933901775e0;\n\
tt2_2 += (tmp_o_i[ 6] - tmp_o_i[ -6]) * -0.001031530213375445369097965e0;\n\
tt2_2 += (tmp_o_i[ 5] - tmp_o_i[ -5]) *  0.006958379116450707495020408e0;\n\
tt2_2 += (tmp_o_i[ 4] - tmp_o_i[ -4]) * -0.03129014783948023634381564e0;\n\
tt2_2 += (tmp_o_i[ 3] - tmp_o_i[ -3]) *  0.1063640682894442760934532e0;\n\
tt2_2 += (tmp_o_i[ 2] - tmp_o_i[ -2]) * -0.3032593514765938346887962e0;\n\
tt2_2 += (tmp_o_i[ 1] - tmp_o_i[ -1]) *  0.8834460460908270942785856e0;\n\
y_r[jg*n + ig] = tmp_y_r[j2*(FILTER_WIDTH+1) + i2] + tt1_1 * scale_1 - tt2_2 * scale_2;\n\
y_i[jg*n + ig] = tmp_y_i[j2*(FILTER_WIDTH+1) + i2] + tt1_2 * scale_1 + tt2_1 * scale_2;\n\
x_r[jg*n + ig] = tmp_o_r[0];\n\
x_i[jg*n + ig] = tmp_o_i[0];\n\
}\n\
";

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
  clSetKernelArg(kernel, i++,sizeof(double)*block_size_j*2*(block_size_i+FILTER_WIDTH+1), NULL);
  clSetKernelArg(kernel, i++,sizeof(double)*block_size_j*2*(block_size_i+1), NULL);
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
  clSetKernelArg(kernel, i++,sizeof(double)*block_size_j*(block_size_i+FILTER_WIDTH+1), NULL);
  clSetKernelArg(kernel, i++,sizeof(double)*block_size_j*(block_size_i+FILTER_WIDTH+1), NULL);
  clSetKernelArg(kernel, i++,sizeof(double)*block_size_j*(block_size_i+1), NULL);
  clSetKernelArg(kernel, i++,sizeof(double)*block_size_j*(block_size_i+1), NULL);
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
    kinetic_kProgram = clCreateProgramWithSource(*context,1,(const char**) &kinetic_k_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(kinetic_kProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build kinetic_k program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(kinetic_kProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
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
