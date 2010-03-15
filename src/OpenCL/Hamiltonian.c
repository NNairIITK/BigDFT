#include "OpenCL_wrappers.h"

void FC_FUNC_(ocl_fulllocham,OCL_FULLLOCHAM)(cl_command_queue *command_queue,
                                          cl_uint *n1, cl_uint *n2, cl_uint *n3,
                                          double *h,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *pot, 
                                          cl_mem *psi, cl_mem *out,
                                          cl_mem *work, cl_mem *kinres,
                                          double *epot,double *ekinpot) {
  cl_uint dimensions[] = { *n1, *n2, *n3 };
  uncompress_d_(command_queue, dimensions,
                               nseg_c, nvctr_c, keyg_c, keyv_c,
                               nseg_f, nvctr_f, keyg_f, keyv_f,
                               psi_c, psi_f, out);
  syn_self_d_(command_queue, n1, n2, n3, out, psi);
  potential_application_d_(command_queue, n1, n2, n3, work, psi, out, pot);
  kinetic_d_(command_queue, n1, n2, n3, h, psi, out, work, kinres);
  ana_self_d_(command_queue, n1, n2, n3, kinres, psi);
  compress_d_(command_queue, dimensions,
	      nseg_c, nvctr_c, keyg_c, keyv_c,
	      nseg_f, nvctr_f, keyg_f, keyv_f,
	      psi_c, psi_f, psi);

  *ekinpot = 0.0;
}

void FC_FUNC_(ocl_fulllocham_generic,OCL_FULLLOCHAM_GENERIC)(cl_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_uint *periodic,
                                          double *h,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *pot, 
                                          cl_mem *psi, cl_mem *out,
                                          cl_mem *work, cl_mem *kinres,
                                          double *epot,double *ekinpot) {
  uncompress_d_(command_queue, dimensions,
                               nseg_c, nvctr_c, keyg_c, keyv_c,
                               nseg_f, nvctr_f, keyg_f, keyv_f,
                               psi_c, psi_f, out);
  syn_self_d_generic_(command_queue, dimensions, periodic, out, psi);
  potential_application_d_generic_(command_queue, dimensions, periodic, work, psi, out, pot);
  kinetic_d_generic_(command_queue, dimensions, periodic, h, psi, out, work, kinres);
  ana_self_d_generic_(command_queue, dimensions, periodic, kinres, psi);
  compress_d_(command_queue, dimensions,
	      nseg_c, nvctr_c, keyg_c, keyv_c,
	      nseg_f, nvctr_f, keyg_f, keyv_f,
	      psi_c, psi_f, psi);

  *ekinpot = 0.0;
}
