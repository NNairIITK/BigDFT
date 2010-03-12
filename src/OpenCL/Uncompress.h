#ifndef UNCOMPRESS_H
#define UNCOMPRESS_H
#include <CL/cl.h>

inline void uncompress_coarse_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,
                                    cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                    cl_mem *psi_c, cl_mem * psi_out);

inline void uncompress_fine_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,
                                    cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                    cl_mem *psi_f, cl_mem * psi_out);

inline void compress_coarse_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,
                                    cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                    cl_mem *psi_c, cl_mem * psi);

inline void compress_fine_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,
                                    cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                    cl_mem *psi_f, cl_mem * psi);

#endif
