#ifndef KINETIC_H
#define KINETIC_H
#include <CL/cl.h>
#include <OpenCL_wrappers.h>

inline void kinetic_generic(cl_kernel kernel, bigdft_command_queue command_queue, cl_uint *n, cl_uint *ndat, double *scale, cl_mem *x_in, cl_mem *x_out, cl_mem *y_in, cl_mem *y_out);

#endif
