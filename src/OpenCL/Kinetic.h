#ifndef KINETIC_H
#define KINETIC_H

char * kinetic1d_program="\
#define FILTER_WIDTH 30\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
__kernel void kinetic1dKernel_d(size_t n, size_t ndat, double scale, __global const double * x_in, __global double * x_out, __global const double * y_in, __global double * y_out, __local double * tmp) {\n\
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
__local double * tmp_o = tmp + j2*(is+FILTER_WIDTH);\n\
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
y_out[jg*n + ig] = y_in[jg + ig*ndat] - scale * conv;\n\
x_out[jg*n + ig] = tmp[0];\n\
};\n\
\n\
__kernel void kinetic1dKernel_l(size_t n, size_t ndat, float scale, __global const float * x_in, __global float * x_out, __global const float * y_in, __global float * y_out, __local float * tmp) {\n\
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

#endif
