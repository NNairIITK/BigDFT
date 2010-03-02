#ifndef KINETIC_K_H
#define KINETIC_K_H

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

#endif
