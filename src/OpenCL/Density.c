#include "OpenCL_wrappers.h"


void FC_FUNC_(ocl_locden,OCL_LOCDEN)(cl_command_queue *command_queue,
                                     cl_uint *dimensions,
                                     double *hfac,
                                     cl_uint *iaddjmp,
                                     cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                     cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                     cl_mem *psi_c, cl_mem *psi_f,
                                     cl_mem *psi, cl_mem *out, cl_mem *work,
                                     cl_mem *pot)
{
   uncompress_d_(command_queue, dimensions,
                 nseg_c, nvctr_c, keyg_c, keyv_c,
                 nseg_f, nvctr_f, keyg_f, keyv_f,
                 psi_c, psi_f, out);
   syn_self_d_(command_queue, dimensions, out, psi);
   magicfilter_den_d_(command_queue, dimensions, work, psi, out);
   cl_uint offset=0;
   cl_uint ndat = 8 * dimensions[0] * dimensions[1] * dimensions[2];
   axpy_offset_self_d_(command_queue, &ndat, hfac, &offset, out, iaddjmp, pot);
}
 
