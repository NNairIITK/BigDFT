#include "OpenCL_wrappers.h"

void apply_hp(bigdft_command_queue *command_queue,
              cl_uint *dimensions,
              double *h,
              double *c,
              cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
              cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
              cl_mem *psi_c, cl_mem *psi_f,
              cl_mem *psi_c_out, cl_mem *psi_f_out,
              cl_mem *psi, cl_mem *out,
              cl_mem *work, cl_mem *kinres) {
  uncompress_scale_d_(command_queue, dimensions, h, c,
                      nseg_c, nvctr_c, keyg_c, keyv_c,
                      nseg_f, nvctr_f, keyg_f, keyv_f,
                      psi_c, psi_f, out);
  syn_self_d_(command_queue, dimensions, out, psi);
  cl_uint n = dimensions[0] * dimensions[1] * dimensions[2] * 8;
  scal_d_(command_queue, &n, c, psi, out);
  kinetic_d_(command_queue, dimensions, h, psi, out, work, kinres);
  ana_self_d_(command_queue, dimensions, kinres, out);
  compress_scale_d_(command_queue, dimensions, h, c,
                      nseg_c, nvctr_c, keyg_c, keyv_c,
                      nseg_f, nvctr_f, keyg_f, keyv_f,
                      psi_c_out, psi_f_out, out);
}

void FC_FUNC_(ocl_preconditioner,OCL_PRECONDITIONER)(bigdft_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          double *h,
                                          double *c,
                                          cl_uint *ncong,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *psi_c_r, cl_mem *psi_f_r,
                                          cl_mem *psi_c_b, cl_mem *psi_f_b,
                                          cl_mem *psi_c_d, cl_mem *psi_f_d,
                                          cl_mem *work1, cl_mem *work2, cl_mem *work3, cl_mem *work4) {
  cl_uint nvctr_f_7 = (*nvctr_f)*7;

  apply_hp(command_queue, dimensions, h, c,
           nseg_c, nvctr_c, keyg_c, keyv_c,
           nseg_f, nvctr_f, keyg_f, keyv_f,
           psi_c, psi_f,
           psi_c_d, psi_f_d,
           work1, work2, work3, work4);

  double alpha, beta, val;
  val = -1.0;
  axpy_d_(command_queue, nvctr_c, &val, psi_c_d, psi_c_b, psi_c_r);
  axpy_d_(command_queue, &nvctr_f_7, &val, psi_f_d, psi_f_b, psi_f_r);
  copy_d_(command_queue, nvctr_c, psi_c_r, psi_c_d);
  copy_d_(command_queue, &nvctr_f_7, psi_f_r, psi_f_d);
//    return;
  double norm_sq_c;
  double norm_sq_f;
  nrm2sq_d_(command_queue, nvctr_c, psi_c_r, work1, work2, &norm_sq_c);
  nrm2sq_d_(command_queue, &nvctr_f_7, psi_f_r, work1, work2, &norm_sq_f);
  double norm_sq = norm_sq_c + norm_sq_f;
  cl_uint i;
  for(i=0; i<*ncong; i++) {
    apply_hp(command_queue, dimensions, h, c,
             nseg_c, nvctr_c, keyg_c, keyv_c,
             nseg_f, nvctr_f, keyg_f, keyv_f, 
             psi_c_d, psi_f_d,
             psi_c_b, psi_f_b,
             work1, work2, work3, work4);
    dot_d_(command_queue, nvctr_c, psi_c_d, psi_c_b, work1, work2, &norm_sq_c);
    dot_d_(command_queue, &nvctr_f_7, psi_f_d, psi_f_b, work1, work2, &norm_sq_f);
    alpha = norm_sq / ( norm_sq_c + norm_sq_f);
    val = alpha;
    axpy_self_d_(command_queue, nvctr_c, &val, psi_c_d, psi_c);
    axpy_self_d_(command_queue, &nvctr_f_7, &val, psi_f_d, psi_f);
    if (i != *ncong-1) {
      val = -alpha;
      axpy_self_d_(command_queue, nvctr_c, &val, psi_c_b, psi_c_r);
      axpy_self_d_(command_queue, &nvctr_f_7, &val, psi_f_b, psi_f_r);
      nrm2sq_d_(command_queue, nvctr_c, psi_c_r, work1, work2, &norm_sq_c);
      nrm2sq_d_(command_queue, &nvctr_f_7, psi_f_r, work1, work2, &norm_sq_f);
      beta = (norm_sq_c + norm_sq_f)/norm_sq;
      val = beta;
      
      axpy_d_(command_queue, nvctr_c, &val, psi_c_d, psi_c_r, psi_c_d);
      axpy_d_(command_queue, &nvctr_f_7, &val, psi_f_d, psi_f_r, psi_f_d);
      norm_sq = norm_sq_c + norm_sq_f;
    }
  }
  scale_psi_d_(command_queue, nvctr_c, nvctr_f, h, c, psi_c,  psi_f);
}

void apply_hp_generic(bigdft_command_queue *command_queue,
              cl_uint *dimensions,
              cl_uint *periodic,
              double *h,
              double *c,
              cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
              cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
              cl_mem *psi_c, cl_mem *psi_f,
              cl_mem *psi_c_out, cl_mem *psi_f_out,
              cl_mem *psi, cl_mem *out,
              cl_mem *work, cl_mem *kinres) {
  uncompress_scale_d_(command_queue, dimensions, h, c,
                      nseg_c, nvctr_c, keyg_c, keyv_c,
                      nseg_f, nvctr_f, keyg_f, keyv_f,
                      psi_c, psi_f, out);
  syn_self_d_generic_(command_queue, dimensions, periodic, out, psi);
  cl_uint n = (dimensions[0] + (periodic[0]?0:7)) * (dimensions[1] + (periodic[1]?0:7)) * (dimensions[2] + (periodic[2]?0:7)) * 8;
  scal_d_(command_queue, &n, c, psi, out);
  kinetic_d_generic_(command_queue, dimensions, periodic, h, psi, out, work, kinres);
  ana_self_d_generic_(command_queue, dimensions, periodic, kinres, out);
  compress_scale_d_(command_queue, dimensions, h, c,
                      nseg_c, nvctr_c, keyg_c, keyv_c,
                      nseg_f, nvctr_f, keyg_f, keyv_f,
                      psi_c_out, psi_f_out, out);
}


void FC_FUNC_(ocl_preconditioner_generic,OCL_PRECONDITIONER_GENERIC)(bigdft_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_uint *periodic,
                                          double *h,
                                          double *c,
                                          cl_uint *ncong,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *psi_c_r, cl_mem *psi_f_r,
                                          cl_mem *psi_c_b, cl_mem *psi_f_b,
                                          cl_mem *psi_c_d, cl_mem *psi_f_d,
                                          cl_mem *work1, cl_mem *work2, cl_mem *work3, cl_mem *work4) {
  cl_uint nvctr_f_7 = (*nvctr_f)*7;

  apply_hp_generic(command_queue, dimensions, periodic, h, c,
           nseg_c, nvctr_c, keyg_c, keyv_c,
           nseg_f, nvctr_f, keyg_f, keyv_f,
           psi_c, psi_f,
           psi_c_d, psi_f_d,
           work1, work2, work3, work4);

  double alpha, beta, val;
  val = -1.0;
  axpy_d_(command_queue, nvctr_c, &val, psi_c_d, psi_c_b, psi_c_r);
  axpy_d_(command_queue, &nvctr_f_7, &val, psi_f_d, psi_f_b, psi_f_r);
  copy_d_(command_queue, nvctr_c, psi_c_r, psi_c_d);
  copy_d_(command_queue, &nvctr_f_7, psi_f_r, psi_f_d);
//    return;
  double norm_sq_c;
  double norm_sq_f;
  nrm2sq_d_(command_queue, nvctr_c, psi_c_r, work1, work2, &norm_sq_c);
  nrm2sq_d_(command_queue, &nvctr_f_7, psi_f_r, work1, work2, &norm_sq_f);
  double norm_sq = norm_sq_c + norm_sq_f;
  cl_uint i;
  for(i=0; i<*ncong; i++) {
    apply_hp_generic(command_queue, dimensions, periodic, h, c,
             nseg_c, nvctr_c, keyg_c, keyv_c,
             nseg_f, nvctr_f, keyg_f, keyv_f, 
             psi_c_d, psi_f_d,
             psi_c_b, psi_f_b,
             work1, work2, work3, work4);
    dot_d_(command_queue, nvctr_c, psi_c_d, psi_c_b, work1, work2, &norm_sq_c);
    dot_d_(command_queue, &nvctr_f_7, psi_f_d, psi_f_b, work1, work2, &norm_sq_f);
    alpha = norm_sq / ( norm_sq_c + norm_sq_f);
    val = alpha;
    axpy_self_d_(command_queue, nvctr_c, &val, psi_c_d, psi_c);
    axpy_self_d_(command_queue, &nvctr_f_7, &val, psi_f_d, psi_f);
    if (i != *ncong-1) {
      val = -alpha;
      axpy_self_d_(command_queue, nvctr_c, &val, psi_c_b, psi_c_r);
      axpy_self_d_(command_queue, &nvctr_f_7, &val, psi_f_b, psi_f_r);
      nrm2sq_d_(command_queue, nvctr_c, psi_c_r, work1, work2, &norm_sq_c);
      nrm2sq_d_(command_queue, &nvctr_f_7, psi_f_r, work1, work2, &norm_sq_f);
      beta = (norm_sq_c + norm_sq_f)/norm_sq;
      val = beta;
      
      axpy_d_(command_queue, nvctr_c, &val, psi_c_d, psi_c_r, psi_c_d);
      axpy_d_(command_queue, &nvctr_f_7, &val, psi_f_d, psi_f_r, psi_f_d);
      norm_sq = norm_sq_c + norm_sq_f;
    }
  }
  scale_psi_d_(command_queue, nvctr_c, nvctr_f, h, c, psi_c,  psi_f);
}

