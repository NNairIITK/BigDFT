//! @file
//!  Bench zgemm with OpenCL
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include "bench_lib.h"


#define MAX_N1 4096
#define MAX_N2 4096
#define MAX_N3 16
#define N1_STEP 4096
#define N2_STEP 4096
#define N3_STEP 16
#define N1_URANGE 0 
#define N2_URANGE 0 
#define N3_URANGE 0
#define N1_USTEP 1
#define N2_USTEP 1
#define N3_USTEP 1
#define NB_ITER 100

void bench_zgemm_simple(cl_uint n1, cl_uint n2, cl_uint n3, cl_double * in1, cl_double * in2, cl_double * out){
  cl_mem a, b, c;
  cl_double2 alpha;
  alpha.x = 1.2;
  alpha.y = 1.1;
  cl_double2 beta;
  beta.x = 1.3;
  beta.y = 1.4;
  cl_uint m = n1;
  cl_uint n = n2;
  cl_uint k = n3;
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
  int i;
  for(i=0;i<NB_ITER;i++)
    gemm_z_(&queue, &transa, &transb, &m, &n, &k, &alpha, &a, &m, &b, &k, &beta, &c, &m);

  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &c, &size_c, out);
  ocl_release_mem_object_(&a);
  ocl_release_mem_object_(&b);
  ocl_release_mem_object_(&c);
}

int main(){
  cl_uint n1,n2,n3;
  cl_uint un1,un2,un3;
  double *in, *out;
  cl_uint size_out = (MAX_N1+N1_URANGE)*(MAX_N2+N2_URANGE)*2;
  cl_uint size_in  = (MAX_N1+N1_URANGE)*(MAX_N3+N3_URANGE)*16*2;
  cl_ulong t0,t1;
  cl_uint device_number;



//  in = (double*) malloc(size*sizeof(double));
//  out = (double*) malloc(size*sizeof(double));


  ocl_create_gpu_context_(&context,&device_number);
  ocl_build_programs_(&context);
  ocl_create_command_queue_(&queue,&context);
  init_event_list_();
 
  cl_mem cmPinnedBufIn = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, size_in*sizeof(double), NULL, NULL);
  cl_mem cmPinnedBufOut = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, size_out*sizeof(double), NULL, NULL);
  in = (double *)clEnqueueMapBuffer(queue->command_queue, cmPinnedBufIn, CL_TRUE, CL_MAP_WRITE, 0, size_in*sizeof(double), 0, NULL, NULL, NULL);
  out = (double *)clEnqueueMapBuffer(queue->command_queue, cmPinnedBufOut, CL_TRUE, CL_MAP_READ, 0, size_out*sizeof(double), 0, NULL, NULL, NULL);
  init_random(in, size_in);

  for( n1 = N1_STEP; n1 <= MAX_N1; n1 += N1_STEP ){
    for( un1 = n1 - N1_URANGE; un1 <= n1 + N1_URANGE; un1 += N1_USTEP){
      for( n2 = N2_STEP; n2 <= MAX_N2; n2 += N2_STEP ){
        for( un2 = n2 - N2_URANGE; un2 <= n2 + N2_URANGE; un2 += N2_USTEP){
          for( n3 = N3_STEP; n3 <= MAX_N3; n3 += N3_STEP ){
            for( un3 = n3 - N3_URANGE; un3 <= n3 + N3_URANGE; un3 += N3_USTEP){
              nanosec_(&t0);
              bench_zgemm_simple(un1,un2,un3,in,in,out);
              nanosec_(&t1);
              printf("%lu usec, (1 run : %lu usec)\n",(t1-t0)/1000,(t1-t0)/(1000*NB_ITER));
              printf("ops : %lu\n",((cl_ulong)un1)*un2*un3*4*2*NB_ITER);
              printf("Gflops : %lf\n",(double)(((cl_ulong)un1)*un2*un3*4*2*NB_ITER)/(double)(t1-t0));
            }
          }
        }
      }
    }
  }

  ocl_release_mem_object_(&cmPinnedBufIn);
  ocl_release_mem_object_(&cmPinnedBufOut);
  print_event_list_();
  ocl_clean_command_queue_(&queue);
  ocl_clean_(&context);
  return 0;
}
