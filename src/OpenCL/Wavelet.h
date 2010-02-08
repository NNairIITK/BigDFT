#ifndef WAVELET_H
#define WAVELET_H

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


#endif
