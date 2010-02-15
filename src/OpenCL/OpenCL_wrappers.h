#ifndef OPENCL_WRAPPERS_H
#define OPENCL_WRAPPERS_H

#include <CL/cl.h>
#include <stdio.h>
#include <stdlib.h>
#include <config.h>
#include <math.h>

#define DEBUG 0
#define PROFILING 1

#define oclErrorCheck(errorCode,message) if(errorCode!=CL_SUCCESS) { fprintf(stderr,"Error(%i) (%s: %s): %s\n", errorCode,__FILE__,__func__,message);exit(1);} 

extern cl_kernel magicfilter1d_kernel_l;
extern cl_kernel magicfilter1d_kernel_d;
extern cl_kernel magicfilter1d_kernel_s_l;
extern cl_kernel magicfilter1d_t_kernel_l;
extern cl_kernel magicfilter1d_pot_kernel_l;
extern cl_kernel magicfilter1d_den_kernel_l;
extern cl_kernel ana1d_kernel_d;
extern cl_kernel syn1d_kernel_d;
extern cl_kernel ana1d_kernel_l;
extern cl_kernel syn1d_kernel_l;
extern cl_kernel kinetic1d_kernel_l;
extern cl_kernel kinetic1d_kernel_d;
extern cl_kernel c_initialize_kernel_l;
extern cl_kernel c_initialize_kernel_d;
extern cl_kernel v_initialize_kernel_l;
extern cl_kernel v_initialize_kernel_d;
extern cl_kernel uncompress_coarse_kernel_l;
extern cl_kernel uncompress_fine_kernel_l;
extern cl_kernel compress_coarse_kernel_l;
extern cl_kernel compress_fine_kernel_l;
extern cl_kernel uncompress_coarse_kernel_d;
extern cl_kernel uncompress_fine_kernel_d;
extern cl_kernel compress_coarse_kernel_d;
extern cl_kernel compress_fine_kernel_d;

void build_magicfilter_kernels(cl_context * context);
void clean_magicfilter_kernels();
void build_kinetic_kernels(cl_context * context);
void clean_kinetic_kernels();
void build_wavelet_kernels(cl_context * context);
void clean_wavelet_kernels();
void build_uncompress_kernels(cl_context * context);
void clean_uncompress_kernels();
void build_initialize_kernels(cl_context * context);
void clean_initialize_kernels();

cl_device_id oclGetFirstDev(cl_context cxGPUContext);
size_t shrRoundUp(size_t group_size, size_t global_size);

typedef struct {
	cl_event e;
	char *comment;
} event;

void FC_FUNC(init_event_list,INIT_EVENT_LIST)();
int addToEventList (event ev);
extern event * event_list;
extern size_t event_number;
void FC_FUNC(print_event_list,PRINT_EVENT_LIST)();
void FC_FUNC_(ocl_build_kernels,OCL_BUILD_KERNELS)(cl_context * context);
void FC_FUNC_(ocl_create_gpu_context,OCL_CREATE_GPU_CONTEXT)(cl_context * context);
void FC_FUNC_(ocl_create_cpu_context,OCL_CREATE_CPU_CONTEXT)(cl_context * context);
void FC_FUNC_(ocl_create_read_buffer,OCL_CREATE_READ_BUFFER)(cl_context *context, cl_uint *size, cl_mem *buff_ptr);
void FC_FUNC_(ocl_create_read_write_buffer,OCL_CREATE_READ_WRITE_BUFFER)(cl_context *context, cl_uint *size, cl_mem *buff_ptr);
void FC_FUNC_(ocl_create_read_buffer_and_copy,OCL_CREATE_READ_BUFFER_AND_COPY)(cl_context *context, cl_uint *size, void *host_ptr, cl_mem *buff_ptr);
void FC_FUNC_(ocl_create_write_buffer,OCL_CREATE_WRITE_BUFFER)(cl_context *context, cl_uint *size, cl_mem *buff_ptr);
void FC_FUNC_(ocl_release_mem_object,OCL_RELEASE_MEM_OBJECT)(cl_mem *buff_ptr);
void FC_FUNC_(ocl_enqueue_read_buffer,OCL_ENQUEUE_READ_BUFFER)(cl_command_queue *command_queue, cl_mem *buffer, cl_uint *size, void *ptr);
void FC_FUNC_(ocl_enqueue_write_buffer,OCL_ENQUEUE_WRITE_BUFFER)(cl_command_queue *command_queue, cl_mem *buffer, cl_uint *size, const void *ptr);
void FC_FUNC_(ocl_create_command_queue,OCL_CREATE_COMMAND_QUEUE)(cl_command_queue *hCmdQueue, cl_context *context);
void FC_FUNC_(ocl_finish,OCL_FINISH)(cl_command_queue *command_queue);
void FC_FUNC_(ocl_enqueue_barrier,OCL_ENQUEUE_BARRIER)(cl_command_queue *command_queue);
void FC_FUNC_(ocl_clean,OCL_CLEAN)(cl_command_queue *command_queue, cl_context *context);
#endif
