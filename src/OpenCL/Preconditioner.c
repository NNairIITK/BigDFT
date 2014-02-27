//! @file
//!  OpenCl wrappers for hamiltonian applications
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include "OpenCL_wrappers.h"

static void apply_hp(bigdft_command_queue *command_queue,
              cl_uint *dimensions,
              double *h,
              double *c,
              cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
              cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
              cl_mem *psi_c, cl_mem *psi_f,
              cl_mem *psi_c_out, cl_mem *psi_f_out,
              cl_mem *psi, cl_mem *out,
              cl_mem *work, cl_mem *kinres) {
  FC_FUNC_(uncompress_scale_d, UNCOMPRESS_SCALE_D)(command_queue, dimensions, h, c,
                      nseg_c, nvctr_c, keyg_c, keyv_c,
                      nseg_f, nvctr_f, keyg_f, keyv_f,
                      psi_c, psi_f, out);
  FC_FUNC_(syn_self_d, SYN_SELF_D)(command_queue, dimensions, out, psi);
  cl_uint n = dimensions[0] * dimensions[1] * dimensions[2] * 8;
  FC_FUNC_(scal_d, SCAL_D)(command_queue, &n, c, psi, out);
  FC_FUNC_(kinetic_d, KINETIC_D)(command_queue, dimensions, h, psi, out, work, kinres);
  FC_FUNC_(ana_self_d, ANA_SELF_D)(command_queue, dimensions, kinres, out);
  FC_FUNC_(compress_scale_d, COMPRESS_SCALE_D)(command_queue, dimensions, h, c,
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
  FC_FUNC_(axpy_d, AXPY_D)(command_queue, nvctr_c, &val, psi_c_d, psi_c_b, psi_c_r);
  FC_FUNC_(axpy_d, AXPY_D)(command_queue, &nvctr_f_7, &val, psi_f_d, psi_f_b, psi_f_r);
  FC_FUNC_(copy_d, COPY_D)(command_queue, nvctr_c, psi_c_r, psi_c_d);
  FC_FUNC_(copy_d, COPY_D)(command_queue, &nvctr_f_7, psi_f_r, psi_f_d);
//    return;
  double norm_sq_c;
  double norm_sq_f;
  FC_FUNC_(nrm2sq_d, NRM2SQ_D)(command_queue, nvctr_c, psi_c_r, work1, work2, &norm_sq_c);
  FC_FUNC_(nrm2sq_d, NRM2SQ_D)(command_queue, &nvctr_f_7, psi_f_r, work1, work2, &norm_sq_f);
  double norm_sq = norm_sq_c + norm_sq_f;
  cl_uint i;
  for(i=0; i<*ncong; i++) {
    apply_hp(command_queue, dimensions, h, c,
             nseg_c, nvctr_c, keyg_c, keyv_c,
             nseg_f, nvctr_f, keyg_f, keyv_f, 
             psi_c_d, psi_f_d,
             psi_c_b, psi_f_b,
             work1, work2, work3, work4);
    FC_FUNC_(dot_d, DOT_D)(command_queue, nvctr_c, psi_c_d, psi_c_b, work1, work2, &norm_sq_c);
    FC_FUNC_(dot_d, DOT_D)(command_queue, &nvctr_f_7, psi_f_d, psi_f_b, work1, work2, &norm_sq_f);
    alpha = norm_sq / ( norm_sq_c + norm_sq_f);
    val = alpha;
    FC_FUNC_(axpy_self_d, AXPY_SELF_D)(command_queue, nvctr_c, &val, psi_c_d, psi_c);
    FC_FUNC_(axpy_self_d, AXPY_SELF_D)(command_queue, &nvctr_f_7, &val, psi_f_d, psi_f);
    if (i != *ncong-1) {
      val = -alpha;
      FC_FUNC_(axpy_self_d, AXPY_SELF_D)(command_queue, nvctr_c, &val, psi_c_b, psi_c_r);
      FC_FUNC_(axpy_self_d, AXPY_SELF_D)(command_queue, &nvctr_f_7, &val, psi_f_b, psi_f_r);
      FC_FUNC_(nrm2sq_d, NRM2SQ_D)(command_queue, nvctr_c, psi_c_r, work1, work2, &norm_sq_c);
      FC_FUNC_(nrm2sq_d, NRM2SQ_D)(command_queue, &nvctr_f_7, psi_f_r, work1, work2, &norm_sq_f);
      beta = (norm_sq_c + norm_sq_f)/norm_sq;
      val = beta;
      
      FC_FUNC_(axpy_d, AXPY_D)(command_queue, nvctr_c, &val, psi_c_d, psi_c_r, psi_c_d);
      FC_FUNC_(axpy_d, AXPY_D)(command_queue, &nvctr_f_7, &val, psi_f_d, psi_f_r, psi_f_d);
      norm_sq = norm_sq_c + norm_sq_f;
    }
  }
  FC_FUNC_(scale_psi_d, SCALE_PSI_D)(command_queue, nvctr_c, nvctr_f, h, c, psi_c,  psi_f);
}

static void apply_hp_generic(bigdft_command_queue *command_queue,
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
  FC_FUNC_(uncompress_scale_d, UNCOMPRESS_SCALE_D)(command_queue, dimensions, h, c,
                      nseg_c, nvctr_c, keyg_c, keyv_c,
                      nseg_f, nvctr_f, keyg_f, keyv_f,
                      psi_c, psi_f, out);
  FC_FUNC_(syn_self_d_generic, SYN_SELF_D_GENERIC)(command_queue, dimensions, periodic, out, psi);
  cl_uint n = (dimensions[0] + (periodic[0]?0:7)) * (dimensions[1] + (periodic[1]?0:7)) * (dimensions[2] + (periodic[2]?0:7)) * 8;
  FC_FUNC_(scal_d, SCAL_D)(command_queue, &n, c, psi, out);
  FC_FUNC_(kinetic_d_generic, KINETIC_D_GENERIC)(command_queue, dimensions, periodic, h, psi, out, work, kinres);
  FC_FUNC_(ana_self_d_generic, ANA_SELF_D_GENERIC)(command_queue, dimensions, periodic, kinres, out);
  FC_FUNC_(compress_scale_d, COMPRESS_SCALE_D)(command_queue, dimensions, h, c,
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
  FC_FUNC_(axpy_d, AXPY_D)(command_queue, nvctr_c, &val, psi_c_d, psi_c_b, psi_c_r);
  FC_FUNC_(axpy_d, AXPY_D)(command_queue, &nvctr_f_7, &val, psi_f_d, psi_f_b, psi_f_r);
  FC_FUNC_(copy_d, COPY_D)(command_queue, nvctr_c, psi_c_r, psi_c_d);
  FC_FUNC_(copy_d, COPY_D)(command_queue, &nvctr_f_7, psi_f_r, psi_f_d);
//    return;
  double norm_sq_c;
  double norm_sq_f;
  FC_FUNC_(nrm2sq_d, NRM2SQ_D)(command_queue, nvctr_c, psi_c_r, work1, work2, &norm_sq_c);
  FC_FUNC_(nrm2sq_d, NRM2SQ_D)(command_queue, &nvctr_f_7, psi_f_r, work1, work2, &norm_sq_f);
  double norm_sq = norm_sq_c + norm_sq_f;
  cl_uint i;
  for(i=0; i<*ncong; i++) {
    apply_hp_generic(command_queue, dimensions, periodic, h, c,
             nseg_c, nvctr_c, keyg_c, keyv_c,
             nseg_f, nvctr_f, keyg_f, keyv_f, 
             psi_c_d, psi_f_d,
             psi_c_b, psi_f_b,
             work1, work2, work3, work4);
    FC_FUNC_(dot_d, DOT_D)(command_queue, nvctr_c, psi_c_d, psi_c_b, work1, work2, &norm_sq_c);
    FC_FUNC_(dot_d, DOT_D)(command_queue, &nvctr_f_7, psi_f_d, psi_f_b, work1, work2, &norm_sq_f);
    alpha = norm_sq / ( norm_sq_c + norm_sq_f);
    val = alpha;
    FC_FUNC_(axpy_self_d, AXPY_SELF_D)(command_queue, nvctr_c, &val, psi_c_d, psi_c);
    FC_FUNC_(axpy_self_d, AXPY_SELF_D)(command_queue, &nvctr_f_7, &val, psi_f_d, psi_f);
    if (i != *ncong-1) {
      val = -alpha;
      FC_FUNC_(axpy_self_d, AXPY_SELF_D)(command_queue, nvctr_c, &val, psi_c_b, psi_c_r);
      FC_FUNC_(axpy_self_d, AXPY_SELF_D)(command_queue, &nvctr_f_7, &val, psi_f_b, psi_f_r);
      FC_FUNC_(nrm2sq_d, NRM2SQ_D)(command_queue, nvctr_c, psi_c_r, work1, work2, &norm_sq_c);
      FC_FUNC_(nrm2sq_d, NRM2SQ_D)(command_queue, &nvctr_f_7, psi_f_r, work1, work2, &norm_sq_f);
      beta = (norm_sq_c + norm_sq_f)/norm_sq;
      val = beta;
      
      FC_FUNC_(axpy_d, AXPY_D)(command_queue, nvctr_c, &val, psi_c_d, psi_c_r, psi_c_d);
      FC_FUNC_(axpy_d, AXPY_D)(command_queue, &nvctr_f_7, &val, psi_f_d, psi_f_r, psi_f_d);
      norm_sq = norm_sq_c + norm_sq_f;
    }
  }
  FC_FUNC_(scale_psi_d, SCALE_PSI_D)(command_queue, nvctr_c, nvctr_f, h, c, psi_c,  psi_f);
}

static void apply_hp_generic_k(bigdft_command_queue *command_queue,
              cl_uint *dimensions,
              cl_uint *periodic,
              double *h,
              double *k,
              double *c,
              cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
              cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
              cl_mem *psi_c_r, cl_mem *psi_f_r,
              cl_mem *psi_c_i, cl_mem *psi_f_i,
              cl_mem *psi_c_out_r, cl_mem *psi_f_out_r,
              cl_mem *psi_c_out_i, cl_mem *psi_f_out_i,
              cl_mem *psi_r, cl_mem *out_r, cl_mem *work_r, cl_mem *kinres_r,
              cl_mem *psi_i, cl_mem *out_i, cl_mem *work_i, cl_mem *kinres_i,
              cl_uint *nspinor) {
  FC_FUNC_(uncompress_scale_d, UNCOMPRESS_SCALE_D)(command_queue, dimensions, h, c,
                      nseg_c, nvctr_c, keyg_c, keyv_c,
                      nseg_f, nvctr_f, keyg_f, keyv_f,
                      psi_c_r, psi_f_r, out_r);
  FC_FUNC_(syn_self_d_generic, SYN_SELF_D_GENERIC)(command_queue, dimensions, periodic, out_r, psi_r);
  if( *nspinor == 2 ) {
    FC_FUNC_(uncompress_scale_d, UNCOMPRESS_SCALE_D)(command_queue, dimensions, h, c,
                        nseg_c, nvctr_c, keyg_c, keyv_c,
                        nseg_f, nvctr_f, keyg_f, keyv_f,
                        psi_c_i, psi_f_i, out_i);
    FC_FUNC_(syn_self_d_generic, SYN_SELF_D_GENERIC)(command_queue, dimensions, periodic, out_i, psi_i);
  }
  cl_uint n = (dimensions[0] + (periodic[0]?0:7)) * (dimensions[1] + (periodic[1]?0:7)) * (dimensions[2] + (periodic[2]?0:7)) * 8;
  if( *nspinor == 2 ) {
    cl_double c2 = .5 * (k[0] * k[0] + k[1] * k[1] + k[2] * k[2])+*c;
    FC_FUNC_(scal_d, SCAL_D)(command_queue, &n, &c2, psi_r, out_r);
    FC_FUNC_(scal_d, SCAL_D)(command_queue, &n, &c2, psi_i, out_i);
    FC_FUNC_(kinetic_k_d_generic, KINETIC_K_D_GENERIC)(command_queue, dimensions, periodic, h, k, psi_r, psi_i, out_r, out_i, work_r, work_i, kinres_r, kinres_i);
  } else {
    FC_FUNC_(scal_d, SCAL_D)(command_queue, &n, c, psi_r, out_r);
    FC_FUNC_(kinetic_d_generic, KINETIC_D_GENERIC)(command_queue, dimensions, periodic, h, psi_r, out_r, work_r, kinres_r);
  }
  FC_FUNC_(ana_self_d_generic, ANA_SELF_D_GENERIC)(command_queue, dimensions, periodic, kinres_r, out_r);
  FC_FUNC_(compress_scale_d, COMPRESS_SCALE_D)(command_queue, dimensions, h, c,
                    nseg_c, nvctr_c, keyg_c, keyv_c,
                    nseg_f, nvctr_f, keyg_f, keyv_f,
                    psi_c_out_r, psi_f_out_r, out_r);
  if( *nspinor == 2 ) {
    FC_FUNC_(ana_self_d_generic, ANA_SELF_D_GENERIC)(command_queue, dimensions, periodic, kinres_i, out_i);
    FC_FUNC_(compress_scale_d, COMPRESS_SCALE_D)(command_queue, dimensions, h, c,
                      nseg_c, nvctr_c, keyg_c, keyv_c,
                      nseg_f, nvctr_f, keyg_f, keyv_f,
                      psi_c_out_i, psi_f_out_i, out_i);
  }
}


void FC_FUNC_(ocl_preconditioner_generic_k,OCL_PRECONDITIONER_GENERIC_K)(bigdft_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_uint *periodic,
                                          double *h,
                                          double *k,
                                          double *c,
                                          cl_uint *ncong,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c_r, cl_mem *psi_f_r,
                                          cl_mem *psi_c_i, cl_mem *psi_f_i,
                                          cl_mem *psi_c_r_r, cl_mem *psi_f_r_r,
                                          cl_mem *psi_c_r_i, cl_mem *psi_f_r_i,
                                          cl_mem *psi_c_b_r, cl_mem *psi_f_b_r,
                                          cl_mem *psi_c_b_i, cl_mem *psi_f_b_i,
                                          cl_mem *psi_c_d_r, cl_mem *psi_f_d_r,
                                          cl_mem *psi_c_d_i, cl_mem *psi_f_d_i,
                                          cl_mem *work1_r, cl_mem *work2_r, cl_mem *work3_r, cl_mem *work4_r,
                                          cl_mem *work1_i, cl_mem *work2_i, cl_mem *work3_i, cl_mem *work4_i,
                                          cl_uint *nspinor, double *norm_sq_cf) {
  cl_uint nvctr_f_7 = (*nvctr_f)*7;

  apply_hp_generic_k(command_queue, dimensions, periodic, h, k, c,
           nseg_c, nvctr_c, keyg_c, keyv_c,
           nseg_f, nvctr_f, keyg_f, keyv_f,
           psi_c_r, psi_f_r,
           psi_c_i, psi_f_i,
           psi_c_d_r, psi_f_d_r,
           psi_c_d_i, psi_f_d_i,
           work1_r, work2_r, work3_r, work4_r,
           work1_i, work2_i, work3_i, work4_i,
           nspinor);

  
  double alpha, beta, val;
  double norm_sq, norm_sq2;
  val = -1.0;
  FC_FUNC_(axpy_d, AXPY_D)(command_queue, nvctr_c, &val, psi_c_d_r, psi_c_b_r, psi_c_r_r);
  FC_FUNC_(axpy_d, AXPY_D)(command_queue, &nvctr_f_7, &val, psi_f_d_r, psi_f_b_r, psi_f_r_r);
  FC_FUNC_(copy_d, COPY_D)(command_queue, nvctr_c, psi_c_r_r, psi_c_d_r);
  FC_FUNC_(copy_d, COPY_D)(command_queue, &nvctr_f_7, psi_f_r_r, psi_f_d_r);
  FC_FUNC_(nrm2sq_d, NRM2SQ_D)(command_queue, nvctr_c, psi_c_r_r, work1_r, work2_r, &norm_sq_cf[0]);
  FC_FUNC_(nrm2sq_d, NRM2SQ_D)(command_queue, &nvctr_f_7, psi_f_r_r, work1_r, work2_r, &norm_sq_cf[1]);
  norm_sq = norm_sq_cf[0] + norm_sq_cf[1];
  if(*nspinor == 2) {
    FC_FUNC_(axpy_d, AXPY_D)(command_queue, nvctr_c, &val, psi_c_d_i, psi_c_b_i, psi_c_r_i);
    FC_FUNC_(axpy_d, AXPY_D)(command_queue, &nvctr_f_7, &val, psi_f_d_i, psi_f_b_i, psi_f_r_i);
    FC_FUNC_(copy_d, COPY_D)(command_queue, nvctr_c, psi_c_r_i, psi_c_d_i);
    FC_FUNC_(copy_d, COPY_D)(command_queue, &nvctr_f_7, psi_f_r_i, psi_f_d_i);
    FC_FUNC_(nrm2sq_d, NRM2SQ_D)(command_queue, nvctr_c, psi_c_r_i, work1_i, work2_i, &norm_sq_cf[0]);
    FC_FUNC_(nrm2sq_d, NRM2SQ_D)(command_queue, &nvctr_f_7, psi_f_r_i, work1_i, work2_i, &norm_sq_cf[1]);
    norm_sq += norm_sq_cf[0] + norm_sq_cf[1];
  }
  cl_uint i;
  for(i=0; i<*ncong; i++) {
    apply_hp_generic_k(command_queue, dimensions, periodic, h, k, c,
             nseg_c, nvctr_c, keyg_c, keyv_c,
             nseg_f, nvctr_f, keyg_f, keyv_f, 
             psi_c_d_r, psi_f_d_r,
             psi_c_d_i, psi_f_d_i,
             psi_c_b_r, psi_f_b_r,
             psi_c_b_i, psi_f_b_i,
             work1_r, work2_r, work3_r, work4_r,
             work1_i, work2_i, work3_i, work4_i,
             nspinor);
    FC_FUNC_(dot_d, DOT_D)(command_queue, nvctr_c, psi_c_d_r, psi_c_b_r, work1_r, work2_r, &norm_sq_cf[0]);
    FC_FUNC_(dot_d, DOT_D)(command_queue, &nvctr_f_7, psi_f_d_r, psi_f_b_r, work1_r, work2_r, &norm_sq_cf[1]);
    norm_sq2 = norm_sq_cf[0] + norm_sq_cf[1];
    if(*nspinor == 2) {
      FC_FUNC_(dot_d, DOT_D)(command_queue, nvctr_c, psi_c_d_i, psi_c_b_i, work1_i, work2_i, &norm_sq_cf[0]);
      FC_FUNC_(dot_d, DOT_D)(command_queue, &nvctr_f_7, psi_f_d_i, psi_f_b_i, work1_i, work2_i, &norm_sq_cf[1]);
      norm_sq2 += norm_sq_cf[0] + norm_sq_cf[1];
    }
    alpha = norm_sq / norm_sq2;
    val = alpha;
    FC_FUNC_(axpy_self_d, AXPY_SELF_D)(command_queue, nvctr_c, &val, psi_c_d_r, psi_c_r);
    FC_FUNC_(axpy_self_d, AXPY_SELF_D)(command_queue, &nvctr_f_7, &val, psi_f_d_r, psi_f_r);
    if(*nspinor == 2) {
      FC_FUNC_(axpy_self_d, AXPY_SELF_D)(command_queue, nvctr_c, &val, psi_c_d_i, psi_c_i);
      FC_FUNC_(axpy_self_d, AXPY_SELF_D)(command_queue, &nvctr_f_7, &val, psi_f_d_i, psi_f_i);
    }
    if (i != *ncong-1) {
      val = -alpha;
      FC_FUNC_(axpy_self_d, AXPY_SELF_D)(command_queue, nvctr_c, &val, psi_c_b_r, psi_c_r_r);
      FC_FUNC_(axpy_self_d, AXPY_SELF_D)(command_queue, &nvctr_f_7, &val, psi_f_b_r, psi_f_r_r);
      FC_FUNC_(nrm2sq_d, NRM2SQ_D)(command_queue, nvctr_c, psi_c_r_r, work1_r, work2_r, &norm_sq_cf[0]);
      FC_FUNC_(nrm2sq_d, NRM2SQ_D)(command_queue, &nvctr_f_7, psi_f_r_r, work1_r, work2_r, &norm_sq_cf[1]);
      norm_sq2 = norm_sq_cf[0] + norm_sq_cf[1];
      if(*nspinor == 2) {
        FC_FUNC_(axpy_self_d, AXPY_SELF_D)(command_queue, nvctr_c, &val, psi_c_b_i, psi_c_r_i);
        FC_FUNC_(axpy_self_d, AXPY_SELF_D)(command_queue, &nvctr_f_7, &val, psi_f_b_i, psi_f_r_i);
        FC_FUNC_(nrm2sq_d, NRM2SQ_D)(command_queue, nvctr_c, psi_c_r_i, work1_i, work2_i, &norm_sq_cf[0]);
        FC_FUNC_(nrm2sq_d, NRM2SQ_D)(command_queue, &nvctr_f_7, psi_f_r_i, work1_i, work2_i, &norm_sq_cf[1]);
        norm_sq2 += norm_sq_cf[0] + norm_sq_cf[1];
      }

      beta = norm_sq2 / norm_sq;
      val = beta;
      
      FC_FUNC_(axpy_d, AXPY_D)(command_queue, nvctr_c, &val, psi_c_d_r, psi_c_r_r, psi_c_d_r);
      FC_FUNC_(axpy_d, AXPY_D)(command_queue, &nvctr_f_7, &val, psi_f_d_r, psi_f_r_r, psi_f_d_r);
      if(*nspinor == 2) {
        FC_FUNC_(axpy_d, AXPY_D)(command_queue, nvctr_c, &val, psi_c_d_i, psi_c_r_i, psi_c_d_i);
        FC_FUNC_(axpy_d, AXPY_D)(command_queue, &nvctr_f_7, &val, psi_f_d_i, psi_f_r_i, psi_f_d_i);
      }
      norm_sq = norm_sq2;
    }
  }
  FC_FUNC_(scale_psi_d, SCALE_PSI_D)(command_queue, nvctr_c, nvctr_f, h, c, psi_c_r,  psi_f_r);
  if(*nspinor == 2) {
    FC_FUNC_(scale_psi_d, SCALE_PSI_D)(command_queue, nvctr_c, nvctr_f, h, c, psi_c_i,  psi_f_i);
  }
}

