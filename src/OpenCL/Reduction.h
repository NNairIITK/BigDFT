#ifndef REDUCTION_H
#define REDUCTION_H
#include <CL/cl.h>

void inline reduction_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *out);
void inline scal_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, double *alpha, cl_mem *in, cl_mem *inout);
void inline axpy_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, double *alpha, cl_mem *in, cl_mem *inout);

#endif
