#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include "Kinetic_k_Generator.h"
#include "OpenCL_wrappers.h"

#define FILTER_WIDTH 32
#define BUFFER_WIDTH 8
static void generate_header(std::stringstream &program){
  program<<"#ifdef cl_khr_fp64\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#elif defined (cl_amd_fp64)\n\
#pragma OPENCL EXTENSION cl_amd_fp64: enable \n\
#endif\n\
#define FILTER_WIDTH "<<FILTER_WIDTH<<"\n\
#define FILT1_14 -6.924474940639200152025730585882e-18\n\
#define FILT1_13  2.70800493626319438269856689037647576e-13\n\
#define FILT1_12 -5.813879830282540547959250667e-11\n\
#define FILT1_11 -1.05857055496741470373494132287e-8\n\
#define FILT1_10 -3.7230763047369275848791496973044e-7\n\
#define FILT1_9   2.0904234952920365957922889447361e-6\n\
#define FILT1_8  -0.2398228524507599670405555359023135e-4\n\
#define FILT1_7   0.45167920287502235349480037639758496e-3\n\
#define FILT1_6  -0.409765689342633823899327051188315485e-2\n\
#define FILT1_5   0.02207029188482255523789911295638968409e0\n\
#define FILT1_4  -0.0822663999742123340987663521e0\n\
#define FILT1_3   0.2371780582153805636239247476e0\n\
#define FILT1_2  -0.6156141465570069496314853949e0\n\
#define FILT1_1   2.2191465938911163898794546405e0\n\
#define FILT1_0  -3.5536922899131901941296809374e0\n\
#define FILT2_14 -1.585464751677102510097179e-19\n\
#define FILT2_13  1.240078536096648534547439e-14\n\
#define FILT2_12 -7.252206916665149851135592e-13\n\
#define FILT2_11 -9.697184925637300947553069e-10\n\
#define FILT2_10 -7.207948238588481597101904e-8\n\
#define FILT2_9   3.993810456408053712133667e-8\n\
#define FILT2_8   2.451992111053665419191564e-7\n\
#define FILT2_7   0.00007667706908380351933901775e0\n\
#define FILT2_6  -0.001031530213375445369097965e0\n\
#define FILT2_5   0.006958379116450707495020408e0\n\
#define FILT2_4  -0.03129014783948023634381564e0\n\
#define FILT2_3   0.1063640682894442760934532e0\n\
#define FILT2_2  -0.3032593514765938346887962e0\n\
#define FILT2_1   0.8834460460908270942785856e0\n";
}

static void generate_interleaved_filters(std::stringstream &program){
  program<<"#define filter1(tt,tmp) \
tt = mad(tmp[2*14] + tmp[2*-14], FILT1_14, tt);\
tt = mad(tmp[2*13] + tmp[2*-13], FILT1_13, tt);\
tt = mad(tmp[2*12] + tmp[2*-12], FILT1_12, tt);\
tt = mad(tmp[2*11] + tmp[2*-11], FILT1_11, tt);\
tt = mad(tmp[2*10] + tmp[2*-10], FILT1_10, tt);\
tt = mad(tmp[2* 9] + tmp[2* -9], FILT1_9 , tt);\
tt = mad(tmp[2* 8] + tmp[2* -8], FILT1_8 , tt);\
tt = mad(tmp[2* 7] + tmp[2* -7], FILT1_7 , tt);\
tt = mad(tmp[2* 6] + tmp[2* -6], FILT1_6 , tt);\
tt = mad(tmp[2* 5] + tmp[2* -5], FILT1_5 , tt);\
tt = mad(tmp[2* 4] + tmp[2* -4], FILT1_4 , tt);\
tt = mad(tmp[2* 3] + tmp[2* -3], FILT1_3 , tt);\
tt = mad(tmp[2* 2] + tmp[2* -2], FILT1_2 , tt);\
tt = mad(tmp[2* 1] + tmp[2* -1], FILT1_1 , tt);\
tt = mad(tmp[2* 0]             , FILT1_0 , tt);\n\
#define filter2(tt,tmp) \
tt = mad(tmp[2*14] - tmp[2*-14], FILT2_14, tt);\
tt = mad(tmp[2*13] - tmp[2*-13], FILT2_13, tt);\
tt = mad(tmp[2*12] - tmp[2*-12], FILT2_12, tt);\
tt = mad(tmp[2*11] - tmp[2*-11], FILT2_11, tt);\
tt = mad(tmp[2*10] - tmp[2*-10], FILT2_10, tt);\
tt = mad(tmp[2* 9] - tmp[2* -9], FILT2_9,  tt);\
tt = mad(tmp[2* 8] - tmp[2* -8], FILT2_8,  tt);\
tt = mad(tmp[2* 7] - tmp[2* -7], FILT2_7,  tt);\
tt = mad(tmp[2* 6] - tmp[2* -6], FILT2_6,  tt);\
tt = mad(tmp[2* 5] - tmp[2* -5], FILT2_5,  tt);\
tt = mad(tmp[2* 4] - tmp[2* -4], FILT2_4,  tt);\
tt = mad(tmp[2* 3] - tmp[2* -3], FILT2_3,  tt);\
tt = mad(tmp[2* 2] - tmp[2* -2], FILT2_2,  tt);\
tt = mad(tmp[2* 1] - tmp[2* -1], FILT2_1,  tt);\n";
}

static void generate_kinetic_k1dKernel(std::stringstream &program){
  program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<BUFFER_WIDTH<<", 1))) void kinetic_k1dKernel_d(uint n, uint ndat, double scale_1, double scale_2, __global const double *  restrict x_in, __global double *  restrict x, __global const double *  restrict y_in, __global double *  restrict y) {\n\
__local double tmp_r["<<BUFFER_WIDTH<<"*(2*"<<FILTER_WIDTH<<"+1)];\n\
__local double tmp_i["<<BUFFER_WIDTH<<"*(2*"<<FILTER_WIDTH<<"+1)];\n\
__local double tmp_y_r["<<BUFFER_WIDTH<<"*("<<FILTER_WIDTH<<"+1)];\n\
__local double tmp_y_i["<<BUFFER_WIDTH<<"*("<<FILTER_WIDTH<<"+1)];\n\
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
tmp_y_r[i2t * (FILTER_WIDTH + 1) + j2t] = y_in[2*jgt+2*igt*ndat];\n\
tmp_y_i[i2t * (FILTER_WIDTH + 1) + j2t] = y_in[2*jgt+2*igt*ndat+1];\n\
igt -= FILTER_WIDTH/2;\n\
if ( igt < 0 ) {\n\
  tmp_r[i2t * (2 * FILTER_WIDTH + 1) + j2t] = x_in[2*jgt + 2 * ( n + igt ) * ndat];\n\
  tmp_i[i2t * (2 * FILTER_WIDTH + 1) + j2t] = x_in[2*jgt + 2 * ( n + igt ) * ndat + 1];\n\
} else {\n\
  tmp_r[i2t * (2 * FILTER_WIDTH + 1) + j2t] = x_in[2*jgt + 2 * igt * ndat];\n\
  tmp_i[i2t * (2 * FILTER_WIDTH + 1) + j2t] = x_in[2*jgt + 2 * igt * ndat+1];\n\
}\n\
igt += FILTER_WIDTH;\n\
if ( igt >= n ) {\n\
  tmp_r[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = x_in[2*jgt + 2 * ( igt - n ) * ndat];\n\
  tmp_i[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = x_in[2*jgt + 2 * ( igt - n ) * ndat + 1];\n\
} else {\n\
  tmp_r[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = x_in[2*jgt +  2 * igt * ndat];\n\
  tmp_i[i2t * (2 * FILTER_WIDTH + 1) + j2t + FILTER_WIDTH] = x_in[2*jgt +  2 * igt * ndat + 1];\n\
}\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
__local const double * restrict tmp_o_r = tmp_r + j2*(2 * FILTER_WIDTH + 1) + FILTER_WIDTH/2 + i2;\n\
__local const double * restrict tmp_o_i = tmp_i + j2*(2 * FILTER_WIDTH + 1) + FILTER_WIDTH/2 + i2;\n\
double tt1_1=0.0;\n\
filter1(tt1_1,tmp_o_r);\n\
double tt2_2=0.0;\n\
filter2(tt2_2,tmp_o_i);\n\
double tt1_2=0.0;\n\
filter1(tt1_2,tmp_o_i);\n\
double tt2_1=0.0;\n\
filter2(tt2_1,tmp_o_r);\n\
y[jg*2*n + 2*ig] = tmp_y_r[j2*(FILTER_WIDTH+1) + i2] + tt1_1 * scale_1 - tt2_2 * scale_2;\n\
y[jg*2*n + 2*ig + 1] = tmp_y_i[j2*(FILTER_WIDTH+1) + i2] + tt1_2 * scale_1 + tt2_1 * scale_2;\n\
x[jg*2*n + 2*ig] = tmp_o_r[0];\n\
x[jg*2*n + 2*ig + 1] = tmp_o_i[0];\n\
}\n";
}

static void generate_filters_nomad(std::stringstream &program){
  program<<"#undef filter1\n\
#define filter1(tt,tmp) \
tt += (tmp[14] + tmp[-14])*FILT1_14;\
tt += (tmp[13] + tmp[-13])*FILT1_13;\
tt += (tmp[12] + tmp[-12])*FILT1_12;\
tt += (tmp[11] + tmp[-11])*FILT1_11;\
tt += (tmp[10] + tmp[-10])*FILT1_10;\
tt += (tmp[ 9] + tmp[ -9])*FILT1_9;\
tt += (tmp[ 8] + tmp[ -8])*FILT1_8;\
tt += (tmp[ 7] + tmp[ -7])*FILT1_7;\
tt += (tmp[ 6] + tmp[ -6])*FILT1_6;\
tt += (tmp[ 5] + tmp[ -5])*FILT1_5;\
tt += (tmp[ 4] + tmp[ -4])*FILT1_4;\
tt += (tmp[ 3] + tmp[ -3])*FILT1_3;\
tt += (tmp[ 2] + tmp[ -2])*FILT1_2;\
tt += (tmp[ 1] + tmp[ -1])*FILT1_1;\
tt += (tmp[ 0]           )*FILT1_0;\n\
#undef filter2\n\
#define filter2(tt,tmp) \
tt += (tmp[14] - tmp[-14])*FILT2_14;\
tt += (tmp[13] - tmp[-13])*FILT2_13;\
tt += (tmp[12] - tmp[-12])*FILT2_12;\
tt += (tmp[11] - tmp[-11])*FILT2_11;\
tt += (tmp[10] - tmp[-10])*FILT2_10;\
tt += (tmp[ 9] - tmp[ -9])*FILT2_9;\
tt += (tmp[ 8] - tmp[ -8])*FILT2_8;\
tt += (tmp[ 7] - tmp[ -7])*FILT2_7;\
tt += (tmp[ 6] - tmp[ -6])*FILT2_6;\
tt += (tmp[ 5] - tmp[ -5])*FILT2_5;\
tt += (tmp[ 4] - tmp[ -4])*FILT2_4;\
tt += (tmp[ 3] - tmp[ -3])*FILT2_3;\
tt += (tmp[ 2] - tmp[ -2])*FILT2_2;\
tt += (tmp[ 1] - tmp[ -1])*FILT2_1;\n";
}

static void generate_filters(std::stringstream &program){
  program<<"#undef filter1\n\
#define filter1(tt,tmp) \
tt = mad(tmp[14] + tmp[-14], FILT1_14, tt);\
tt = mad(tmp[13] + tmp[-13], FILT1_13, tt);\
tt = mad(tmp[12] + tmp[-12], FILT1_12, tt);\
tt = mad(tmp[11] + tmp[-11], FILT1_11, tt);\
tt = mad(tmp[10] + tmp[-10], FILT1_10, tt);\
tt = mad(tmp[ 9] + tmp[ -9], FILT1_9,  tt);\
tt = mad(tmp[ 8] + tmp[ -8], FILT1_8,  tt);\
tt = mad(tmp[ 7] + tmp[ -7], FILT1_7,  tt);\
tt = mad(tmp[ 6] + tmp[ -6], FILT1_6,  tt);\
tt = mad(tmp[ 5] + tmp[ -5], FILT1_5,  tt);\
tt = mad(tmp[ 4] + tmp[ -4], FILT1_4,  tt);\
tt = mad(tmp[ 3] + tmp[ -3], FILT1_3,  tt);\
tt = mad(tmp[ 2] + tmp[ -2], FILT1_2,  tt);\
tt = mad(tmp[ 1] + tmp[ -1], FILT1_1,  tt);\
tt = mad(tmp[ 0]           , FILT1_0,  tt);\n\
#undef filter2\n\
#define filter2(tt,tmp) \
tt = mad(tmp[14] - tmp[-14], FILT2_14, tt);\
tt = mad(tmp[13] - tmp[-13], FILT2_13, tt);\
tt = mad(tmp[12] - tmp[-12], FILT2_12, tt);\
tt = mad(tmp[11] - tmp[-11], FILT2_11, tt);\
tt = mad(tmp[10] - tmp[-10], FILT2_10, tt);\
tt = mad(tmp[ 9] - tmp[ -9], FILT2_9,  tt);\
tt = mad(tmp[ 8] - tmp[ -8], FILT2_8,  tt);\
tt = mad(tmp[ 7] - tmp[ -7], FILT2_7,  tt);\
tt = mad(tmp[ 6] - tmp[ -6], FILT2_6,  tt);\
tt = mad(tmp[ 5] - tmp[ -5], FILT2_5,  tt);\
tt = mad(tmp[ 4] - tmp[ -4], FILT2_4,  tt);\
tt = mad(tmp[ 3] - tmp[ -3], FILT2_3,  tt);\
tt = mad(tmp[ 2] - tmp[ -2], FILT2_2,  tt);\
tt = mad(tmp[ 1] - tmp[ -1], FILT2_1,  tt);\n";
}

static void generate_kinetic_k1dKernel_2(std::stringstream &program){
  program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<BUFFER_WIDTH<<", 1))) void kinetic_k1dKernel_d_2(uint n, uint ndat, double scale_1, double scale_2, __global const double *  restrict x_in_r, __global const double *  restrict x_in_i, __global double *  restrict x_r, __global double *  restrict x_i, __global const double *  restrict y_in_r, __global const double *  restrict y_in_i, __global double *  restrict y_r, __global double *  restrict y_i ) {\n\
__local double tmp_r["<<BUFFER_WIDTH<<"*(2*"<<FILTER_WIDTH<<"+1)];\n\
__local double tmp_i["<<BUFFER_WIDTH<<"*(2*"<<FILTER_WIDTH<<"+1)];\n\
__local double tmp_y_r["<<BUFFER_WIDTH<<"*("<<FILTER_WIDTH<<"+1)];\n\
__local double tmp_y_i["<<BUFFER_WIDTH<<"*("<<FILTER_WIDTH<<"+1)];\n\
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
filter1(tt1_1,tmp_o_r);\n\
double tt1_2=0;\n\
filter1(tt1_2,tmp_o_i);\n\
double tt2_1=0;\n\
filter2(tt2_1,tmp_o_r);\n\
double tt2_2=0;\n\
filter2(tt2_2,tmp_o_i);\n\
y_r[jg*n + ig] = tmp_y_r[j2*(FILTER_WIDTH+1) + i2] + tt1_1 * scale_1 - tt2_2 * scale_2;\n\
y_i[jg*n + ig] = tmp_y_i[j2*(FILTER_WIDTH+1) + i2] + tt1_2 * scale_1 + tt2_1 * scale_2;\n\
x_r[jg*n + ig] = tmp_o_r[0];\n\
x_i[jg*n + ig] = tmp_o_i[0];\n\
}\n";
}

static void generate_kinetic_k1d_fKernel_2(std::stringstream &program){
  program<<"__kernel __attribute__((reqd_work_group_size("<<FILTER_WIDTH<<","<<BUFFER_WIDTH<<", 1))) void kinetic_k1d_fKernel_d_2(uint n, uint ndat, double scale_1, double scale_2, __global const double *  restrict x_in_r, __global const double *  restrict x_in_i, __global double *  restrict x_r, __global double *  restrict x_i, __global const double *  restrict y_in_r, __global const double *  restrict y_in_i, __global double *  restrict y_r, __global double *  restrict y_i ) {\n\
__local double tmp_r["<<BUFFER_WIDTH<<"*(2*"<<FILTER_WIDTH<<"+1)];\n\
__local double tmp_i["<<BUFFER_WIDTH<<"*(2*"<<FILTER_WIDTH<<"+1)];\n\
__local double tmp_y_r["<<BUFFER_WIDTH<<"*("<<FILTER_WIDTH<<"+1)];\n\
__local double tmp_y_i["<<BUFFER_WIDTH<<"*("<<FILTER_WIDTH<<"+1)];\n\
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
filter1(tt1_1,tmp_o_r);\n\
double tt1_2=0;\n\
filter1(tt1_2,tmp_o_i);\n\
double tt2_1=0;\n\
filter2(tt2_1,tmp_o_r);\n\
double tt2_2=0;\n\
filter2(tt2_2,tmp_o_i);\n\
y_r[jg*n + ig] = tmp_y_r[j2*(FILTER_WIDTH+1) + i2] + tt1_1 * scale_1 - tt2_2 * scale_2;\n\
y_i[jg*n + ig] = tmp_y_i[j2*(FILTER_WIDTH+1) + i2] + tt1_2 * scale_1 + tt2_1 * scale_2;\n\
x_r[jg*n + ig] = tmp_o_r[0];\n\
x_i[jg*n + ig] = tmp_o_i[0];\n\
}\n";
}




extern "C" char* generate_kinetic_k_program(struct bigdft_device_infos * infos){
  char * output;
  std::stringstream program;

  generate_header(program);
  generate_filters(program);
  generate_kinetic_k1dKernel(program);
  generate_kinetic_k1dKernel_2(program);
  generate_kinetic_k1d_fKernel_2(program);

  output = (char *)malloc((program.str().size()+1)*sizeof(char));
  strcpy(output, program.str().c_str());
  return output;
}

