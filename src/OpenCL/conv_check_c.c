#include "OpenCL_wrappers.h"
#include <stdlib.h>

#define SIZE_X 64
#define SIZE_Y 2048


void init_random(double * data, size_t size) {
  size_t i;
  for(i=0; i<size; i++)
    data[i] = (double)rand() / (double)RAND_MAX;
}

int main(){
  cl_context context;
  cl_command_queue queue;
  cl_mem psi_GPU,work_GPU,work2_GPU,v_GPU;
  double * in, *out;
  int i;
  cl_uint size = SIZE_X * SIZE_Y * sizeof(double);
  cl_uint size_x = SIZE_X;
  cl_uint size_x2 = SIZE_X/2;
  cl_uint size_y = SIZE_Y;
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



  print_event_list_();
  ocl_clean_(&queue, &context);
}
