#ifndef INITIALIZE_H
#define INITIALIZE_H

char * c_initialize_program="\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
__kernel void c_initializeKernel_d(size_t n, __global const double * x_in, __global double * y_in, double c) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
y_in[ig] = x_in[ig] * c;\n\
};\n\
__kernel void v_initializeKernel_d(size_t n, __global double * y_in, double v) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
y_in[ig] = v;\n\
};\n\
";

#endif
