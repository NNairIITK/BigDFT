#ifndef REDUCTION_H
#define REDUCTION_H
#include <CL/cl.h>

void inline reduction_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *out);
void inline scal_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, double *alpha, cl_mem *in, cl_mem *inout);
void inline axpy_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, double *alpha, cl_mem *in, cl_mem *inout);
void inline axpy_offset_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, double *alpha, 
                                                                                   cl_uint *offset_x, cl_mem *x,
                                                                                   cl_uint *offset_y, cl_mem *y,
                                                                                   cl_uint *offset_out, cl_mem *out);
void inline dot_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, cl_mem *x, cl_mem *y, cl_mem *out);
void inline copy_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *out);
void inline set_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n, double *val, cl_mem *x);
void inline gemm_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *m, cl_uint *n, cl_uint *k, double *alpha, cl_mem *a, cl_mem *b, double *beta, cl_mem *c);

#endif
