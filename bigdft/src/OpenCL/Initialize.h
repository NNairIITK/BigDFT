#ifndef INITIALIZE_H
#define INITIALIZE_H
#include <CL/cl.h>

void inline c_initialize_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *ndat, cl_mem *in, cl_mem *inout, double *c);
void inline v_initialize_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *ndat, cl_mem *out, double *c);
void inline p_initialize_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *ndat, cl_mem *in, cl_mem *out);

#endif
