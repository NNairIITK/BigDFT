#ifndef INITIALIZE_H
#define INITIALIZE_H
#include <CL/cl.h>

void c_initialize_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *ndat, cl_mem *in, cl_mem *inout, double *c);
void v_initialize_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *ndat, cl_mem *out, double *c);
void p_initialize_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *ndat, cl_mem *in, cl_mem *out);

#endif
