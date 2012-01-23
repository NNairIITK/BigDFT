#ifndef MAGICFILTER_H
#define MAGICFILTER_H
#include <CL/cl.h>


inline void magicfilter_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);

#endif
