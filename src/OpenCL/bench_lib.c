#include <stdlib.h>
#include <assert.h>
#include "bench_lib.h"

void init_random(double * data, size_t size) {
  size_t i;
  for(i=0; i<size; i++)
    data[i] = (double)rand() / (double)RAND_MAX;
}

cl_context context;
cl_command_queue queue;

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

void bench_magicfilter1d(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out) {
  cl_mem psi_GPU,work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size = n*ndat*sizeof(double);

  ocl_create_write_buffer_(&context, &size, &psi_GPU);
  ocl_create_read_buffer_(&context, &size, &work_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size, in);
  magicfilter1d_d_(&queue,&n,&ndat,&work_GPU,&psi_GPU);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &psi_GPU, &size, out);
  ocl_release_mem_object_(&psi_GPU);
  ocl_release_mem_object_(&work_GPU);
}

void bench_magicfiltershrink1d(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out) {
  cl_mem psi_GPU,work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size_i = n*ndat*sizeof(double);
  cl_uint size_o = (n-15)*ndat*sizeof(double);
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

void bench_magicfiltergrow1d(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out) {
  cl_mem psi_GPU,work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size_o = n*ndat*sizeof(double);
  cl_uint size_i = (n-15)*ndat*sizeof(double);
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

void bench_kinetic1d(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out) {
  cl_mem psi_GPU,work_GPU,work2_GPU,v_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size = n*ndat*sizeof(double);
  double hx = 0.1;
  double c = 0.0;
  double ekinGPU;

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

void bench_ana1d(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out) {
  cl_mem psi_GPU,work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size = n*ndat*sizeof(double);
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

void bench_anashrink1d(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out) {
  cl_mem psi_GPU,work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size_o = n*ndat*sizeof(double);
  cl_uint size_f = (n-14)*ndat*sizeof(double);
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

void bench_syn1d(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out) {
  cl_mem psi_GPU,work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size = n*ndat*sizeof(double);
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

void bench_syngrow1d(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out) {
  cl_mem psi_GPU,work_GPU;
  cl_uint n = n1;
  cl_uint ndat = n2*n3;
  cl_uint size_f = n*ndat*sizeof(double);
  cl_uint size_o = (n-14)*ndat*sizeof(double);
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

void bench_uncompress(cl_uint n1, cl_uint n2, cl_uint n3, cl_uint nseg, cl_uint nvctr_cf, cl_uint * keyg, cl_uint * keyv, double * psi_in, double * psi_out) {
  cl_mem work_GPU,psi_c_GPU,psi_f_GPU,keyg_GPU, keyv_GPU ;
  cl_uint size;

  size = nvctr_cf*sizeof(double);
  ocl_create_read_buffer_(&context, &size, &psi_c_GPU);
  ocl_enqueue_write_buffer_(&queue, &psi_c_GPU, &size, psi_in);
  size = nvctr_cf*sizeof(double)*7;
  ocl_create_read_buffer_(&context, &size, &psi_f_GPU);
  ocl_enqueue_write_buffer_(&queue, &psi_f_GPU, &size, psi_in + nvctr_cf);
  size = nseg*sizeof(cl_uint)*2;
  ocl_create_read_buffer_(&context, &size, &keyg_GPU);
  ocl_enqueue_write_buffer_(&queue, &keyg_GPU, &size, keyg);
  size = nseg*sizeof(cl_uint);
  ocl_create_read_buffer_(&context, &size, &keyv_GPU);
  ocl_enqueue_write_buffer_(&queue, &keyv_GPU, &size, keyv);
  size = n1*n2*n3*sizeof(double)*8;
  ocl_create_write_buffer_(&context, &size, &work_GPU);
  cl_uint dimensions[] = { n1, n2, n3};
  uncompress_d_(&queue , dimensions,
                &nseg, &nvctr_cf, &keyg_GPU, &keyv_GPU,
                &nseg, &nvctr_cf, &keyg_GPU, &keyv_GPU,
                &psi_c_GPU, &psi_f_GPU, &work_GPU);
  ocl_finish_(&queue);
  size = n1*n2*n3*sizeof(double)*8;
  ocl_enqueue_read_buffer_(&queue, &work_GPU, &size, psi_out);
  ocl_release_mem_object_(&psi_c_GPU);
  ocl_release_mem_object_(&psi_f_GPU);
  ocl_release_mem_object_(&keyg_GPU);
  ocl_release_mem_object_(&keyv_GPU);
  ocl_release_mem_object_(&work_GPU);

  size = nvctr_cf*sizeof(double);
  ocl_create_write_buffer_(&context, &size, &psi_c_GPU);
  size = nvctr_cf*sizeof(double)*7;
  ocl_create_write_buffer_(&context, &size, &psi_f_GPU);
  size = nseg*sizeof(cl_uint)*2;
  ocl_create_read_buffer_(&context, &size, &keyg_GPU);
  ocl_enqueue_write_buffer_(&queue, &keyg_GPU, &size, keyg);
  size = nseg*sizeof(cl_uint);
  ocl_create_read_buffer_(&context, &size, &keyv_GPU);
  ocl_enqueue_write_buffer_(&queue, &keyv_GPU, &size, keyv);
  size = n1*n2*n3*sizeof(double)*8;
  ocl_create_read_buffer_(&context, &size, &work_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size, psi_out);
  compress_d_(&queue , dimensions,
              &nseg, &nvctr_cf, &keyg_GPU, &keyv_GPU,
              &nseg, &nvctr_cf, &keyg_GPU, &keyv_GPU,
              &psi_c_GPU, &psi_f_GPU, &work_GPU);
  ocl_finish_(&queue);
  size = nvctr_cf*sizeof(double);
  ocl_enqueue_read_buffer_(&queue, &psi_c_GPU, &size, psi_in);
  size = nvctr_cf*sizeof(double)*7;
  ocl_enqueue_read_buffer_(&queue, &psi_f_GPU, &size, psi_in + nvctr_cf);
  ocl_release_mem_object_(&psi_c_GPU);
  ocl_release_mem_object_(&psi_f_GPU);
  ocl_release_mem_object_(&keyg_GPU);
  ocl_release_mem_object_(&keyv_GPU);
  ocl_release_mem_object_(&work_GPU);


}
