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

  uncompress_d_(command_queue, n1, n2, n3,
                               nseg_c, nvctr_c, keyg_c, keyv_c,
                               nseg_f, nvctr_f, keyg_f, keyv_f,
                               psi_c, psi_f, psi);
  analysis_d_(command_queue, n1, n2, n3, psi, out);
  potential_application_d_(command_queue, n1, n2, n3, work, out, psi, pot);
  kinetic_d_(command_queue, n1, n2, n3, h, out, psi, work, kinres);
  synthesis_d_(command_queue, n1, n2, n3, kinres, out);
  compress_d_(command_queue, n1, n2, n3,
                             nseg_c, nvctr_c, keyg_c, keyv_c,
                             nseg_f, nvctr_f, keyg_f, keyv_f,
                             psi_c, psi_f, out);
  *ekinpot = 0.0;
}
