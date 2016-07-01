//! @file
//!  Libray for bench convolutions (OpenCL)
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include <stdlib.h>
#include <assert.h>
#include "bench_lib.h"

void init_random(cl_double * data, size_t size) {
  size_t i;
  for(i=0; i<size; i++)
    data[i] = (cl_double)rand() / (cl_double)RAND_MAX;
}

bigdft_context context;
bigdft_command_queue queue;

void init_key(cl_int *keyg, cl_int *keyv, cl_int nseg, cl_int *nvctr_cf) {
  int i,j;
  cl_uint tmp=1;
  for(j=0;j<nseg;j++)
    for(i=0;i<nseg;i++){
      keyg[nseg*2*j+2*i] = j*nseg*nseg + i*nseg + 6;
      keyg[nseg*2*j+2*i+1] = j*nseg*nseg + i*nseg + nseg - 5;
      keyv[nseg*j+i] = tmp;
      tmp += nseg - 10;
  //    printf("%d %d : %d %d %d\n",j,i,keyg[nseg*2*j+i],keyg[nseg*2*j+i+1],keyv[nseg*j+i]);
  }
  *nvctr_cf = keyv[nseg*nseg-1] + nseg - 10 -1;
}

void bench_magicfilter1d(cl_uint n1, cl_uint n2, cl_uint n3, cl_double * in, cl_double * out) {
  cl_mem psi_GPU, work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size = n*ndat*sizeof(cl_double);

  ocl_create_write_buffer_(&context, &size, &psi_GPU);
  ocl_create_read_buffer_(&context, &size, &work_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size, in);
  magicfilter1d_d_(&queue,&n,&ndat,&work_GPU,&psi_GPU);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &psi_GPU, &size, out);
  ocl_release_mem_object_(&psi_GPU);
  ocl_release_mem_object_(&work_GPU);
}
void bench_magicfilter1d_straight(cl_uint n1, cl_uint n2, cl_uint n3, cl_double * in, cl_double * out) {
  cl_mem psi_GPU,work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size = n*ndat*sizeof(cl_double);

  ocl_create_write_buffer_(&context, &size, &psi_GPU);
  ocl_create_read_buffer_(&context, &size, &work_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size, in);
  magicfilter1d_straight_d_(&queue,&n,&ndat,&work_GPU,&psi_GPU);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &psi_GPU, &size, out);
  ocl_release_mem_object_(&psi_GPU);
  ocl_release_mem_object_(&work_GPU);
}
void bench_magicfilter1d_block(cl_uint n1, cl_uint n2, cl_uint n3, cl_double * in, cl_double * out) {
  cl_mem psi_GPU,work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size = n*ndat*sizeof(cl_double);

  ocl_create_write_buffer_(&context, &size, &psi_GPU);
  ocl_create_read_buffer_(&context, &size, &work_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size, in);
  magicfilter1d_block_d_(&queue,&n,&ndat,&work_GPU,&psi_GPU);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &psi_GPU, &size, out);
  ocl_release_mem_object_(&psi_GPU);
  ocl_release_mem_object_(&work_GPU);
}

void bench_magicfiltershrink1d(cl_uint n1, cl_uint n2, cl_uint n3, cl_double * in, cl_double * out) {
  cl_mem psi_GPU,work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size_i = n*ndat*sizeof(cl_double);
  cl_uint size_o = (n-15)*ndat*sizeof(cl_double);
  n = n1 - 15;

  ocl_create_write_buffer_(&context, &size_o, &psi_GPU);
  ocl_create_read_buffer_(&context, &size_i, &work_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size_i, in);
  magicfiltershrink1d_d_(&queue,&n,&ndat,&work_GPU,&psi_GPU);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &psi_GPU, &size_o, out);
  ocl_release_mem_object_(&psi_GPU);
  ocl_release_mem_object_(&work_GPU);
}

void bench_magicfiltergrow1d(cl_uint n1, cl_uint n2, cl_uint n3, cl_double * in, cl_double * out) {
  cl_mem psi_GPU,work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size_o = n*ndat*sizeof(cl_double);
  cl_uint size_i = (n-15)*ndat*sizeof(cl_double);
  n = n1 - 15;

  ocl_create_write_buffer_(&context, &size_o, &psi_GPU);
  ocl_create_read_buffer_(&context, &size_i, &work_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size_i, in);
  magicfiltergrow1d_d_(&queue,&n,&ndat,&work_GPU,&psi_GPU);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &psi_GPU, &size_o, out);
  ocl_release_mem_object_(&psi_GPU);
  ocl_release_mem_object_(&work_GPU);
}

void bench_kinetic1d(cl_uint n1, cl_uint n2, cl_uint n3, cl_double * in, cl_double * out) {
  cl_mem psi_GPU,work_GPU,work2_GPU,v_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size = n*ndat*sizeof(cl_double);
  cl_double hx = 0.1;
  cl_double c = 0.0;
  cl_double ekinGPU;

  ocl_create_write_buffer_(&context, &size, &psi_GPU);
  ocl_create_read_buffer_(&context, &size, &work_GPU);
  ocl_create_write_buffer_(&context, &size, &work2_GPU);
  ocl_create_read_write_buffer_(&context, &size, &v_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size, in);
  kinetic1d_d_(&queue, &n, &ndat, &hx, &c, &work_GPU, &psi_GPU, &work2_GPU, &v_GPU, &ekinGPU);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &psi_GPU, &size, out);
  ocl_release_mem_object_(&psi_GPU);
  ocl_release_mem_object_(&work_GPU);
  ocl_release_mem_object_(&work2_GPU);
  ocl_release_mem_object_(&v_GPU);
}

void bench_ana1d(cl_uint n1, cl_uint n2, cl_uint n3, cl_double * in, cl_double * out) {
  cl_mem psi_GPU,work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size = n*ndat*sizeof(cl_double);
  n = n1/2;

  ocl_create_write_buffer_(&context, &size, &psi_GPU);
  ocl_create_read_buffer_(&context, &size, &work_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size, in);
  ana1d_d_(&queue, &n, &ndat, &work_GPU, &psi_GPU);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &psi_GPU, &size, out);
  ocl_release_mem_object_(&psi_GPU);
  ocl_release_mem_object_(&work_GPU);
}

void bench_ana1d_block(cl_uint n1, cl_uint n2, cl_uint n3, cl_double * in, cl_double * out) {
  cl_mem psi_GPU,work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size = n*ndat*sizeof(cl_double);
  n = n1/2;

  ocl_create_write_buffer_(&context, &size, &psi_GPU);
  ocl_create_read_buffer_(&context, &size, &work_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size, in);
  ana1d_block_d_(&queue, &n, &ndat, &work_GPU, &psi_GPU);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &psi_GPU, &size, out);
  ocl_release_mem_object_(&psi_GPU);
  ocl_release_mem_object_(&work_GPU);
}


void bench_anashrink1d(cl_uint n1, cl_uint n2, cl_uint n3, cl_double * in, cl_double * out) {
  cl_mem psi_GPU,work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size_o = n*ndat*sizeof(cl_double);
  cl_uint size_f = (n-14)*ndat*sizeof(cl_double);
  n = n1/2-7;


  ocl_create_write_buffer_(&context, &size_f, &psi_GPU);
  ocl_create_read_buffer_(&context, &size_o, &work_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size_o, in);
  anashrink1d_d_(&queue, &n, &ndat, &work_GPU, &psi_GPU);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &psi_GPU, &size_f, out);
  ocl_release_mem_object_(&psi_GPU);
  ocl_release_mem_object_(&work_GPU);

}

void bench_syn1d(cl_uint n1, cl_uint n2, cl_uint n3, cl_double * in, cl_double * out) {
  cl_mem psi_GPU,work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size = n*ndat*sizeof(cl_double);
  n = n1/2;

  ocl_create_write_buffer_(&context, &size, &psi_GPU);
  ocl_create_read_buffer_(&context, &size, &work_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size, in);
  syn1d_d_(&queue, &n,&ndat,&work_GPU,&psi_GPU);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &psi_GPU, &size, out);
  ocl_release_mem_object_(&psi_GPU);
  ocl_release_mem_object_(&work_GPU);

}

void bench_syngrow1d(cl_uint n1, cl_uint n2, cl_uint n3, cl_double * in, cl_double * out) {
  cl_mem psi_GPU,work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size_f = n*ndat*sizeof(cl_double);
  cl_uint size_o = (n-14)*ndat*sizeof(cl_double);
  n = n1/2-7;

  ocl_create_write_buffer_(&context, &size_f, &psi_GPU);
  ocl_create_read_buffer_(&context, &size_o, &work_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size_o, in);
  syngrow1d_d_(&queue, &n, &ndat, &work_GPU, &psi_GPU);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &psi_GPU, &size_f, out);
  ocl_release_mem_object_(&psi_GPU);
  ocl_release_mem_object_(&work_GPU);
}

void bench_gemm(cl_uint n1, cl_uint n2, cl_uint n3, cl_double * in1, cl_double * in2, cl_double * out){
  cl_mem a, b, c;
  cl_double alpha = 1.2;
  cl_double beta = 1.3;
  cl_uint m = n1;
  cl_uint n = n2*n3;
  cl_uint k = n1;
  cl_uint size_a = m * k * sizeof(cl_double);
  cl_uint size_b = n * k * sizeof(cl_double);
  cl_uint size_c = n * m * sizeof(cl_double);
  ocl_create_write_buffer_(&context, &size_c, &c);
  ocl_create_read_buffer_(&context, &size_b, &b);
  ocl_create_read_buffer_(&context, &size_a, &a);
  ocl_enqueue_write_buffer_(&queue, &a, &size_a, in1);
  ocl_enqueue_write_buffer_(&queue, &b, &size_b, in2);
  char transa = 'n';
  char transb = 'n';
  gemm_d_(&queue, &transa, &transb, &m, &n, &k, &alpha, &a, &m, &b, &k, &beta, &c, &m);
  gemm_block_d_(&queue, &transa, &transb, &m, &n, &k, &alpha, &a, &m, &b, &k, &beta, &c, &m);
  transa = 'n';
  transb = 't';
  gemm_d_(&queue, &transa, &transb, &m, &n, &k, &alpha, &a, &m, &b, &n, &beta, &c, &m);
  transa = 't';
  transb = 'n';
  gemm_d_(&queue, &transa, &transb, &m, &n, &k, &alpha, &a, &k, &b, &k, &beta, &c, &m);
  transa = 't';
  transb = 't';
  gemm_d_(&queue, &transa, &transb, &m, &n, &k, &alpha, &a, &k, &b, &n, &beta, &c, &m);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &c, &size_c, out);
  ocl_release_mem_object_(&a);
  ocl_release_mem_object_(&b);
  ocl_release_mem_object_(&c);
}
void bench_zgemm(cl_uint n1, cl_uint n2, cl_uint n3, cl_double * in1, cl_double * in2, cl_double * out){
  cl_mem a, b, c;
  cl_double2 alpha;
  alpha.x = 1.2;
  alpha.y = 1.1;
  cl_double2 beta;
  beta.x = 1.3;
  beta.y = 1.4;
  cl_uint m = n1/2;
  cl_uint n = n2*n3;
  cl_uint k = n1/2;
  cl_uint size_a = m * k * 2 * sizeof(cl_double);
  cl_uint size_b = n * k * 2 * sizeof(cl_double);
  cl_uint size_c = n * m * 2 * sizeof(cl_double);
  ocl_create_write_buffer_(&context, &size_c, &c);
  ocl_create_read_buffer_(&context, &size_b, &b);
  ocl_create_read_buffer_(&context, &size_a, &a);
  ocl_enqueue_write_buffer_(&queue, &a, &size_a, in1);
  ocl_enqueue_write_buffer_(&queue, &b, &size_b, in2);
  char transa = 'n';
  char transb = 'n';
  gemm_z_(&queue, &transa, &transb, &m, &n, &k, &alpha, &a, &m, &b, &k, &beta, &c, &m);
  transa = 'n';
  transb = 't';
  gemm_z_(&queue, &transa, &transb, &m, &n, &k, &alpha, &a, &m, &b, &n, &beta, &c, &m);
  transa = 't';
  transb = 'n';
  gemm_z_(&queue, &transa, &transb, &m, &n, &k, &alpha, &a, &k, &b, &k, &beta, &c, &m);
  transa = 't';
  transb = 't';
  gemm_z_(&queue, &transa, &transb, &m, &n, &k, &alpha, &a, &k, &b, &n, &beta, &c, &m);
  transa = 'n';
  transb = 'c';
  gemm_z_(&queue, &transa, &transb, &m, &n, &k, &alpha, &a, &m, &b, &n, &beta, &c, &m);
  transa = 'c';
  transb = 'n';
  gemm_z_(&queue, &transa, &transb, &m, &n, &k, &alpha, &a, &k, &b, &k, &beta, &c, &m);
  transa = 'c';
  transb = 't';
  gemm_z_(&queue, &transa, &transb, &m, &n, &k, &alpha, &a, &k, &b, &n, &beta, &c, &m);
  transa = 'c';
  transb = 'c';
  gemm_z_(&queue, &transa, &transb, &m, &n, &k, &alpha, &a, &k, &b, &n, &beta, &c, &m);
  transa = 't';
  transb = 'c';
  gemm_z_(&queue, &transa, &transb, &m, &n, &k, &alpha, &a, &k, &b, &n, &beta, &c, &m);

  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &c, &size_c, out);
  ocl_release_mem_object_(&a);
  ocl_release_mem_object_(&b);
  ocl_release_mem_object_(&c);
}
void bench_uncompress(cl_uint n1, cl_uint n2, cl_uint n3, cl_uint nseg, cl_uint nvctr_cf, cl_uint * keyg, cl_uint * keyv, cl_double * psi_in, cl_double * psi_out) {
  cl_mem work_GPU,psi_c_GPU,psi_f_GPU,keyg_GPU, keyv_GPU ;
  cl_uint size;

  size = nvctr_cf*sizeof(cl_double);
  ocl_create_read_buffer_(&context, &size, &psi_c_GPU);
  ocl_enqueue_write_buffer_(&queue, &psi_c_GPU, &size, psi_in);
  size = nvctr_cf*sizeof(cl_double)*7;
  ocl_create_read_buffer_(&context, &size, &psi_f_GPU);
  ocl_enqueue_write_buffer_(&queue, &psi_f_GPU, &size, psi_in + nvctr_cf);
  size = nseg*sizeof(cl_uint)*2;
  ocl_create_read_buffer_(&context, &size, &keyg_GPU);
  ocl_enqueue_write_buffer_(&queue, &keyg_GPU, &size, keyg);
  size = nseg*sizeof(cl_uint);
  ocl_create_read_buffer_(&context, &size, &keyv_GPU);
  ocl_enqueue_write_buffer_(&queue, &keyv_GPU, &size, keyv);
  size = n1*n2*n3*sizeof(cl_double)*8;
  ocl_create_write_buffer_(&context, &size, &work_GPU);
  cl_uint dimensions[] = { n1, n2, n3};
  uncompress_d_(&queue , dimensions,
                &nseg, &nvctr_cf, &keyg_GPU, &keyv_GPU,
                &nseg, &nvctr_cf, &keyg_GPU, &keyv_GPU,
                &psi_c_GPU, &psi_f_GPU, &work_GPU);
  ocl_finish_(&queue);
  size = n1*n2*n3*sizeof(cl_double)*8;
  ocl_enqueue_read_buffer_(&queue, &work_GPU, &size, psi_out);
  ocl_release_mem_object_(&psi_c_GPU);
  ocl_release_mem_object_(&psi_f_GPU);
  ocl_release_mem_object_(&keyg_GPU);
  ocl_release_mem_object_(&keyv_GPU);
  ocl_release_mem_object_(&work_GPU);

  size = nvctr_cf*sizeof(cl_double);
  ocl_create_write_buffer_(&context, &size, &psi_c_GPU);
  size = nvctr_cf*sizeof(cl_double)*7;
  ocl_create_write_buffer_(&context, &size, &psi_f_GPU);
  size = nseg*sizeof(cl_uint)*2;
  ocl_create_read_buffer_(&context, &size, &keyg_GPU);
  ocl_enqueue_write_buffer_(&queue, &keyg_GPU, &size, keyg);
  size = nseg*sizeof(cl_uint);
  ocl_create_read_buffer_(&context, &size, &keyv_GPU);
  ocl_enqueue_write_buffer_(&queue, &keyv_GPU, &size, keyv);
  size = n1*n2*n3*sizeof(cl_double)*8;
  ocl_create_read_buffer_(&context, &size, &work_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size, psi_out);
  compress_d_(&queue , dimensions,
              &nseg, &nvctr_cf, &keyg_GPU, &keyv_GPU,
              &nseg, &nvctr_cf, &keyg_GPU, &keyv_GPU,
              &psi_c_GPU, &psi_f_GPU, &work_GPU);
  ocl_finish_(&queue);
  size = nvctr_cf*sizeof(cl_double);
  ocl_enqueue_read_buffer_(&queue, &psi_c_GPU, &size, psi_in);
  size = nvctr_cf*sizeof(cl_double)*7;
  ocl_enqueue_read_buffer_(&queue, &psi_f_GPU, &size, psi_in + nvctr_cf);
  ocl_release_mem_object_(&psi_c_GPU);
  ocl_release_mem_object_(&psi_f_GPU);
  ocl_release_mem_object_(&keyg_GPU);
  ocl_release_mem_object_(&keyv_GPU);
  ocl_release_mem_object_(&work_GPU);


}
