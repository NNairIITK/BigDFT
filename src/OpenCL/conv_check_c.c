#include "OpenCL_wrappers.h"
#include <stdlib.h>
#include <assert.h>

#define SIZE_X 64
#define SIZE_Y 2048


void init_random(double * data, size_t size) {
  size_t i;
  for(i=0; i<size; i++)
    data[i] = (double)rand() / (double)RAND_MAX;
}

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

int main(){
  cl_context context;
  cl_command_queue queue;
  cl_mem psi_GPU,work_GPU,work2_GPU,v_GPU,psi_c_GPU,psi_f_GPU,keyg_GPU, keyv_GPU ;
  double *in, *out;
  double *psi;
  double *psi_in, *psi_out;
  int i,j,k;
  cl_int size = SIZE_X * SIZE_Y * sizeof(double);
  cl_int size_x = SIZE_X;
  cl_int size_x2 = SIZE_X/2;
  cl_int size_y = SIZE_Y;
  cl_int nseg = SIZE_X * SIZE_X;
  cl_int *keyg, *keyv;
  cl_int nvctr_cf;
  keyg = (cl_uint*)malloc(SIZE_X * SIZE_X * 2 * sizeof(cl_uint));
  //printf("%p\n",keyg);
  keyv = (cl_uint*)malloc(SIZE_X * SIZE_X * sizeof(cl_uint));
  //printf("%p\n",keyv);
  init_key( keyg, keyv, SIZE_X, &nvctr_cf);
  psi = (double *)malloc(nvctr_cf * 8 * sizeof(double));
  psi = (double *)malloc(nvctr_cf * 8 * sizeof(double));
  //printf("%p\n",psi);
  for(i=0;i<nvctr_cf;i++){
    psi[i] = (double)i+1;   
  }
  for(i=0;i<nvctr_cf;i++){
    for(j=0;j<7;j++){
      psi[nvctr_cf+7*i+j]= i+1 + 0.1 * (j+1);
      assert(nvctr_cf+7*i+j < nvctr_cf * 8); 
    }
  }
  for(i=0;i<nvctr_cf*8;i++){
  //   printf("%d %lf\n",i,psi[i]);
  }
  psi_in = (double *)malloc(SIZE_X*SIZE_X*SIZE_X* sizeof(double)*8);
  psi_out = (double *)malloc(SIZE_X*SIZE_X*SIZE_X* sizeof(double)*8);
  printf("nvctr_cf : %d\n",nvctr_cf);

  in = (double*) malloc(size);
  out = (double*) malloc(size);
  double hx = 0.1;
  double c = 0.0;
  double ekinGPU;
  init_random(in, size_x * size_y);

  ocl_create_gpu_context_(&context);
  ocl_create_command_queue_(&queue,&context);
  ocl_build_kernels_(&context);
  init_event_list_();
/*

  ocl_create_write_buffer_(&context, &size, &psi_GPU);
  ocl_create_read_buffer_(&context, &size, &work_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size, in);
  magicfilter1d_d_(&queue,&size_x,&size_y,&work_GPU,&psi_GPU);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &psi_GPU, &size, out);
  ocl_release_mem_object_(&psi_GPU);
  ocl_release_mem_object_(&work_GPU);

  ocl_create_write_buffer_(&context, &size, &psi_GPU);
  ocl_create_read_buffer_(&context, &size, &work_GPU);
  ocl_create_write_buffer_(&context, &size, &work2_GPU);
  ocl_create_read_write_buffer_(&context, &size, &v_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size, in);
  kinetic1d_d_(&queue,&size_x,&size_y,&hx,&c,&work_GPU,&psi_GPU,&work2_GPU,&v_GPU,&ekinGPU);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &psi_GPU, &size, out);
  ocl_release_mem_object_(&psi_GPU);
  ocl_release_mem_object_(&work_GPU);
  ocl_release_mem_object_(&work2_GPU);
  ocl_release_mem_object_(&v_GPU);

  ocl_create_write_buffer_(&context, &size, &psi_GPU);
  ocl_create_read_buffer_(&context, &size, &work_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size, in);
  ana1d_d_(&queue,&size_x2,&size_y,&work_GPU,&psi_GPU);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &psi_GPU, &size, out);
  ocl_release_mem_object_(&psi_GPU);
  ocl_release_mem_object_(&work_GPU);

  ocl_create_write_buffer_(&context, &size, &psi_GPU);
  ocl_create_read_buffer_(&context, &size, &work_GPU);
  ocl_enqueue_write_buffer_(&queue, &work_GPU, &size, in);
  syn1d_d_(&queue, &size_x2,&size_y,&work_GPU,&psi_GPU);
  ocl_finish_(&queue);
  ocl_enqueue_read_buffer_(&queue, &psi_GPU, &size, out);
  ocl_release_mem_object_(&psi_GPU);
  ocl_release_mem_object_(&work_GPU);
*/  
  cl_uint data_size=4;
  size = nvctr_cf*sizeof(double);
  ocl_create_read_buffer_(&context, &size, &psi_c_GPU);
  ocl_enqueue_write_buffer_(&queue, &psi_c_GPU, &size, psi);
//  ocl_enqueue_write_buffer_d_(&queue, &psi_c_GPU, &size, psi, &data_size);
  size = nvctr_cf*sizeof(double)*7;
  ocl_create_read_buffer_(&context, &size, &psi_f_GPU);
  ocl_enqueue_write_buffer_(&queue, &psi_f_GPU, &size, psi + nvctr_cf);
//  ocl_enqueue_write_buffer_d_(&queue, &psi_f_GPU, &size, psi + nvctr_cf, &data_size);
  size = SIZE_X*SIZE_X*sizeof(cl_uint)*2;
  ocl_create_read_buffer_(&context, &size, &keyg_GPU);
//  ocl_enqueue_write_buffer_(&queue, &keyg_GPU, &size, keyg);
  ocl_enqueue_write_buffer_(&queue, &keyg_GPU, &size, keyg);
  size = SIZE_X*SIZE_X*sizeof(cl_uint);
  ocl_create_read_buffer_(&context, &size, &keyv_GPU);
  ocl_enqueue_write_buffer_(&queue, &keyv_GPU, &size, keyv);
  size = SIZE_X*SIZE_X*SIZE_X*sizeof(double)*8;
  ocl_create_write_buffer_(&context, &size, &work_GPU);
  uncompress_d_(&queue , &size_x, &size_x, &size_x,
                &nseg, &nvctr_cf, &keyg_GPU, &keyv_GPU,
                &nseg, &nvctr_cf, &keyg_GPU, &keyv_GPU,
                &psi_c_GPU, &psi_f_GPU, &work_GPU);
  ocl_finish_(&queue);
  size = SIZE_X*SIZE_X*SIZE_X*sizeof(double)*8;
  ocl_enqueue_read_buffer_(&queue, &work_GPU, &size, psi_out);
  ocl_release_mem_object_(&psi_c_GPU);
  ocl_release_mem_object_(&psi_f_GPU);
  ocl_release_mem_object_(&keyg_GPU);
  ocl_release_mem_object_(&keyv_GPU);
  ocl_release_mem_object_(&work_GPU);

for(i=0;i<2*size_x;i++)
  for(j=0;j<2*size_x;j++)
    for(k=0;k<2*size_x;k++){
//      printf("%4d %4d %4d", i+1,j+1,k+1);
//      printf(" %6.1lf",psi_out[(i*size_x*2 +j)*size_x*2 +k]);
//      printf("\n");
    }
  print_event_list_();
  ocl_clean_(&queue, &context);
}
