#ifndef KINETIC_K_H
#define KINETIC_K_H
#include <CL/cl.h>

inline void kinetic_k_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n, cl_uint *ndat, double *scale1, double *scale2, cl_mem *x_in, cl_mem *x_out, cl_mem *y_in, cl_mem *y_out);

#endif
