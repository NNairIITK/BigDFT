#ifndef INITIALIZE_H
#define INITIALIZE_H

char * c_initialize_program="\
__kernel void c_initializeKernel_l(size_t n, size_t ndat, __global const float * x_in, __global float * y_in, float c) {\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
if( jg >= ndat ) return;\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
size_t pos = jg+ig*ndat;\n\
y_in[pos] = x_in[pos] * c;\n\
};\n\
__kernel void v_initializeKernel_l(size_t n, __global float * y_in, float v) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
y_in[ig] = v;\n\
};\n\
";

#endif
