#ifndef UNCOMPRESS_H
#define UNCOMPRESS_H
#include <CL/cl.h>

inline void un_compress_coarse_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *dimensions,
                                    cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                    cl_mem *psi_c, cl_mem * psi);

inline void un_compress_fine_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *dimensions,
                                    cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                    cl_mem *psi_f, cl_mem * psi);

inline void scale_psi_coarse_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *nvctr_c, double *h, double *c, cl_mem *psi_c);
inline void scale_psi_fine_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *nvctr_f, double *h, double *c, cl_mem *psi_f);
inline void un_compress_scale_coarse_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *dimensions, double *h, double *c,
                                    cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                    cl_mem *psi_c, cl_mem * psi);
inline void un_compress_scale_fine_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *dimensions, double *h, double *c,
                                    cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                    cl_mem *psi_f, cl_mem * psi);
#endif
