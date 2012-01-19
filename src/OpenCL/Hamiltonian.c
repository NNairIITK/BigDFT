//! @file
//!  Wrappers for OpenCL
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include "OpenCL_wrappers.h"

void FC_FUNC_(ocl_fulllocham,OCL_FULLLOCHAM)(bigdft_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_double *h,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *pot, 
                                          cl_mem *psi, cl_mem *out,
                                          cl_mem *work, cl_mem *kinres,
                                          cl_double *epot,cl_double *ekinpot) {
  uncompress_d_(command_queue, dimensions,
                               nseg_c, nvctr_c, keyg_c, keyv_c,
                               nseg_f, nvctr_f, keyg_f, keyv_f,
                               psi_c, psi_f, out);
  syn_self_d_(command_queue, dimensions, out, psi);
  potential_application_d_(command_queue, dimensions, work, psi, out, pot);
  kinetic_d_(command_queue, dimensions, h, psi, out, work, kinres);
  ana_self_d_(command_queue, dimensions, kinres, psi);
  compress_d_(command_queue, dimensions,
	      nseg_c, nvctr_c, keyg_c, keyv_c,
	      nseg_f, nvctr_f, keyg_f, keyv_f,
	      psi_c, psi_f, psi);

  *ekinpot = 0.0;
  *epot = 0.0;
}

void FC_FUNC_(ocl_fulllocham_generic,OCL_FULLLOCHAM_GENERIC)(bigdft_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_uint *periodic,
                                          cl_double *h,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *pot, 
                                          cl_mem *psi, cl_mem *out,
                                          cl_mem *work, cl_mem *kinres,
                                          cl_double *epot,cl_double *ekinpot) {
  uncompress_d_(command_queue, dimensions,
                               nseg_c, nvctr_c, keyg_c, keyv_c,
                               nseg_f, nvctr_f, keyg_f, keyv_f,
                               psi_c, psi_f, out);
  syn_self_d_generic_(command_queue, dimensions, periodic, out, psi);
  potential_application_d_generic_(command_queue, dimensions, periodic, work, kinres, psi, out, pot, epot);
  kinetic_d_generic_(command_queue, dimensions, periodic, h, psi, out, work, kinres);
  cl_uint n1 = dimensions[0] * 2;             
  cl_uint n2 = dimensions[1] * 2;
  cl_uint n3 = dimensions[2] * 2;
  if( !periodic[0] ) n1 += 2*7;
  if( !periodic[1] ) n2 += 2*7;
  if( !periodic[2] ) n3 += 2*7;
  cl_uint ndat = n1*n2*n3;
  dot_d_async_(command_queue,&ndat, work, kinres, psi, out, ekinpot);
  ana_self_d_generic_(command_queue, dimensions, periodic, kinres, psi);
  compress_d_(command_queue, dimensions,
	      nseg_c, nvctr_c, keyg_c, keyv_c,
	      nseg_f, nvctr_f, keyg_f, keyv_f,
	      psi_c, psi_f, psi);
}

void FC_FUNC_(ocl_isf_to_daub,OCL_ISF_TO_DAUB)(bigdft_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_uint *periodic,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *psi, cl_mem *out,
                                          cl_mem *work, cl_mem *kinres) {
  magic_filter_t_3d_generic_(command_queue, dimensions, periodic, work, kinres, psi, out);
  ana_self_d_generic_(command_queue, dimensions, periodic, out, psi);
  compress_d_(command_queue, dimensions,
	      nseg_c, nvctr_c, keyg_c, keyv_c,
	      nseg_f, nvctr_f, keyg_f, keyv_f,
	      psi_c, psi_f, psi);
}
void FC_FUNC_(ocl_daub_to_isf,OCL_DAUB_TO_ISF)(bigdft_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_uint *periodic,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *psi, cl_mem *out,
                                          cl_mem *work, cl_mem *kinres) {
  uncompress_d_(command_queue, dimensions,
                               nseg_c, nvctr_c, keyg_c, keyv_c,
                               nseg_f, nvctr_f, keyg_f, keyv_f,
                               psi_c, psi_f, out);
  syn_self_d_generic_(command_queue, dimensions, periodic, out, psi);
  magic_filter_3d_generic_(command_queue, dimensions, periodic, work, kinres, psi, out);
}

void FC_FUNC_(ocl_fulllocham_generic_k,OCL_FULLLOCHAM_GENERIC_K)(bigdft_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_uint *periodic,
                                          cl_double *h,
                                          cl_double *k,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c_r, cl_mem *psi_f_r,
                                          cl_mem *psi_c_i, cl_mem *psi_f_i,
                                          cl_mem *pot, 
                                          cl_mem *psi_r, cl_mem *out_r, cl_mem *work_r,
                                          cl_mem *psi_i, cl_mem *out_i, cl_mem *work_i,
                                          cl_mem *kinres_r,
                                          cl_mem *kinres_i,
                                          cl_uint *nspinor,
                                          cl_double *epot,cl_double *ekinpot){
//                                          double *epot_i,double *ekinpot_i) {
  epot[1] = 0.0;
  ekinpot[1] = 0.0;
  uncompress_d_(command_queue, dimensions,
                               nseg_c, nvctr_c, keyg_c, keyv_c,
                               nseg_f, nvctr_f, keyg_f, keyv_f,
                               psi_c_r, psi_f_r, out_r);
  syn_self_d_generic_(command_queue, dimensions, periodic, out_r, psi_r);
  potential_application_d_generic_(command_queue, dimensions, periodic, work_r, kinres_r, psi_r, out_r, pot, &epot[0]);
  if(*nspinor==2){
    uncompress_d_(command_queue, dimensions,
                               nseg_c, nvctr_c, keyg_c, keyv_c,
                               nseg_f, nvctr_f, keyg_f, keyv_f,
                               psi_c_i, psi_f_i, out_i);
    syn_self_d_generic_(command_queue, dimensions, periodic, out_i, psi_i);
    potential_application_d_generic_(command_queue, dimensions, periodic, work_i, kinres_i, psi_i, out_i, pot, &epot[1]);
  }


  cl_uint n1 = dimensions[0] * 2;             
  cl_uint n2 = dimensions[1] * 2;
  cl_uint n3 = dimensions[2] * 2;
  if( !periodic[0] ) n1 += 2*7;
  if( !periodic[1] ) n2 += 2*7;
  if( !periodic[2] ) n3 += 2*7;
  cl_uint ndat = n1*n2*n3;

  if(*nspinor==2) {
    double c = .5 * (k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
    axpy_self_d_(command_queue, &ndat, &c, psi_r, out_r);
    axpy_self_d_(command_queue, &ndat, &c, psi_i, out_i);
    kinetic_k_d_generic_(command_queue, dimensions, periodic, h, k, psi_r, psi_i, out_r, out_i, work_r, work_i, kinres_r, kinres_i);
  } else {
    kinetic_d_generic_(command_queue, dimensions, periodic, h, psi_r, out_r, work_r, kinres_r);
  }

  dot_d_async_(command_queue,&ndat, work_r, kinres_r, psi_r, out_r, &ekinpot[0]);
  ana_self_d_generic_(command_queue, dimensions, periodic, kinres_r, psi_r);
  compress_d_(command_queue, dimensions,
	      nseg_c, nvctr_c, keyg_c, keyv_c,
	      nseg_f, nvctr_f, keyg_f, keyv_f,
	      psi_c_r, psi_f_r, psi_r);
  if(*nspinor==2){
    dot_d_async_(command_queue,&ndat, work_i, kinres_i, psi_i, out_i, &ekinpot[1]);
    ana_self_d_generic_(command_queue, dimensions, periodic, kinres_i, psi_i);
    compress_d_(command_queue, dimensions,
	      nseg_c, nvctr_c, keyg_c, keyv_c,
	      nseg_f, nvctr_f, keyg_f, keyv_f,
	      psi_c_i, psi_f_i, psi_i);

  }
}
