//! @file
//!  OpenCL wrappers
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include "OpenCL_wrappers.h"


//void FC_FUNC_(rdtsc,RDTSC)(cl_ulong * t){
//  rdtscll(*t);
//}

cl_device_id oclGetFirstDev(cl_context cxGPUContext)
{
    size_t szParmDataBytes;
    cl_device_id* cdDevices;

    // get the list of GPU devices associated with context
    clGetContextInfo(cxGPUContext, CL_CONTEXT_DEVICES, 0, NULL, &szParmDataBytes);
    cdDevices = (cl_device_id*) malloc(szParmDataBytes);

    clGetContextInfo(cxGPUContext, CL_CONTEXT_DEVICES, szParmDataBytes, cdDevices, NULL);

    cl_device_id first = cdDevices[0];
    free(cdDevices);

    return first;
}

void FC_FUNC_(ocl_build_programs,OCL_BUILD_PROGRAMS)(cl_context * context) {
    build_magicfilter_programs(context);
    build_benchmark_programs(context);
    build_kinetic_programs(context);
    build_wavelet_programs(context);
    build_uncompress_programs(context);
    build_initialize_programs(context);
    build_reduction_programs(context);
    build_fft_programs(context);
}

void create_kernels(struct bigdft_kernels *kernels){
    create_magicfilter_kernels(kernels);
    create_benchmark_kernels(kernels);
    create_kinetic_kernels(kernels);
    create_wavelet_kernels(kernels);
    create_uncompress_kernels(kernels);
    create_initialize_kernels(kernels);
    create_reduction_kernels(kernels);
    create_fft_kernels(kernels);
}


// WARNING : devices are supposed to be uniform in a context
void get_context_devices_infos(cl_context * context, struct bigdft_device_infos * infos){
    cl_uint device_number;

#ifdef CL_VERSION_1_1
    clGetContextInfo(*context, CL_CONTEXT_NUM_DEVICES, sizeof(device_number), &device_number, NULL);
#else
    size_t nContextDescriptorSize;
    clGetContextInfo(*context, CL_CONTEXT_DEVICES, 0, 0, &nContextDescriptorSize);
    device_number = nContextDescriptorSize/sizeof(cl_device_id);
#endif
    cl_device_id * aDevices = (cl_device_id *) malloc(sizeof(cl_device_id)*device_number);
    clGetContextInfo(*context, CL_CONTEXT_DEVICES, sizeof(cl_device_id)*device_number, aDevices, 0);

    get_device_infos(aDevices[0], infos);

    free(aDevices);
}

void get_device_infos(cl_device_id device, struct bigdft_device_infos * infos){
    clGetDeviceInfo(device, CL_DEVICE_TYPE, sizeof(infos->DEVICE_TYPE), &(infos->DEVICE_TYPE), NULL);
    clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(infos->LOCAL_MEM_SIZE), &(infos->LOCAL_MEM_SIZE), NULL);
    clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(infos->MAX_WORK_GROUP_SIZE), &(infos->MAX_WORK_GROUP_SIZE), NULL);
    clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(infos->MAX_COMPUTE_UNITS), &(infos->MAX_COMPUTE_UNITS), NULL);
    clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(infos->NAME), infos->NAME, NULL);
}

void FC_FUNC_(ocl_create_gpu_context,OCL_CREATE_GPU_CONTEXT)(cl_context * context) {
    cl_int ciErrNum = CL_SUCCESS;
    cl_platform_id platform_id;
    clGetPlatformIDs(1, &platform_id, NULL);
    cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform_id, 0 };
    *context = clCreateContextFromType( properties , CL_DEVICE_TYPE_GPU, NULL, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p\n",*context);
#endif
    oclErrorCheck(ciErrNum,"Failed to create GPU context!");
}

void FC_FUNC_(ocl_create_cpu_context,OCL_CREATE_CPU_CONTEXT)(cl_context * context) {
    cl_int ciErrNum = CL_SUCCESS;
    cl_platform_id platform_id;
    clGetPlatformIDs(1, &platform_id, NULL);
    cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform_id, 0 };
    *context = clCreateContextFromType( properties, CL_DEVICE_TYPE_CPU, NULL, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p\n",*context);
#endif
    oclErrorCheck(ciErrNum,"Failed to create CPU context!");
}

void FC_FUNC_(ocl_create_read_buffer,OCL_CREATE_READ_BUFFER)(cl_context *context, cl_uint *size, cl_mem *buff_ptr) {
    cl_int ciErrNum = CL_SUCCESS;
    *buff_ptr = clCreateBuffer( *context, CL_MEM_READ_ONLY, *size, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, memory address: %p, size: %lu\n",*context,*buff_ptr,(long unsigned)*size);
#endif
    oclErrorCheck(ciErrNum,"Failed to create read buffer!");
}

void FC_FUNC_(ocl_create_read_write_buffer,OCL_CREATE_READ_WRITE_BUFFER)(cl_context *context, cl_uint *size, cl_mem *buff_ptr) {
    cl_int ciErrNum = CL_SUCCESS;
    *buff_ptr = clCreateBuffer( *context, CL_MEM_READ_WRITE, *size, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, memory address: %p, size: %lu\n",*context,*buff_ptr,(long unsigned)*size);
#endif
    oclErrorCheck(ciErrNum,"Failed to create read_write buffer!");
}

void FC_FUNC_(ocl_create_read_buffer_and_copy,OCL_CREATE_READ_BUFFER_AND_COPY)(cl_context *context, cl_uint *size, void *host_ptr, cl_mem *buff_ptr) {
    cl_int ciErrNum = CL_SUCCESS;
    *buff_ptr = clCreateBuffer( *context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, *size, host_ptr, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create initialized read buffer!");
}

void FC_FUNC_(ocl_create_write_buffer,OCL_CREATE_WRITE_BUFFER)(cl_context *context, cl_uint *size, cl_mem *buff_ptr) {
    cl_int ciErrNum = CL_SUCCESS;
    *buff_ptr = clCreateBuffer( *context, CL_MEM_WRITE_ONLY, *size, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, memory address: %p, size: %lu\n",*context,*buff_ptr,(long unsigned)*size);
#endif
    oclErrorCheck(ciErrNum,"Failed to create write buffer!");
}

void FC_FUNC_(ocl_release_mem_object,OCL_RELEASE_MEM_OBJECT)(cl_mem *buff_ptr) {
    cl_int ciErrNum = clReleaseMemObject( *buff_ptr);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("memory address: %p\n",*buff_ptr);
#endif
    oclErrorCheck(ciErrNum,"Failed to release buffer!");
}

void FC_FUNC_(ocl_enqueue_read_buffer,OCL_ENQUEUE_READ_BUFFER)(bigdft_command_queue *command_queue, cl_mem *buffer, cl_uint *size, void *ptr){
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, memory address: %p, size: %lu, target: %p\n",(*command_queue)->command_queue,*buffer,(long unsigned)*size, ptr);
#endif
    cl_int ciErrNum = clEnqueueReadBuffer( (*command_queue)->command_queue, *buffer, CL_TRUE, 0, *size, ptr, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue read buffer!");
}

void FC_FUNC_(ocl_enqueue_write_buffer,OCL_ENQUEUE_WRITE_BUFFER)(bigdft_command_queue *command_queue, cl_mem *buffer, cl_uint *size, const void *ptr){
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, memory address: %p, size: %lu, source: %p\n",(*command_queue)->command_queue,*buffer,(long unsigned)*size, ptr);
#endif
    cl_int ciErrNum = clEnqueueWriteBuffer( (*command_queue)->command_queue, *buffer, CL_TRUE, 0, *size, ptr, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue write buffer!");
}

void FC_FUNC_(ocl_enqueue_read_buffer_async,OCL_ENQUEUE_READ_BUFFER_ASYNC)(bigdft_command_queue *command_queue, cl_mem *buffer, cl_uint *size, void *ptr){
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, memory address: %p, size: %lu, target: %p\n",(*command_queue)->command_queue,*buffer,(long unsigned)*size, ptr);
#endif
    cl_int ciErrNum = clEnqueueReadBuffer( (*command_queue)->command_queue, *buffer, CL_FALSE, 0, *size, ptr, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue read buffer!");
}

void FC_FUNC_(ocl_enqueue_write_buffer_async,OCL_ENQUEUE_WRITE_BUFFER_ASYNC)(bigdft_command_queue *command_queue, cl_mem *buffer, cl_uint *size, const void *ptr){
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, memory address: %p, size: %lu, source: %p\n",(*command_queue)->command_queue,*buffer,(long unsigned)*size, ptr);
#endif
    cl_int ciErrNum = clEnqueueWriteBuffer( (*command_queue)->command_queue, *buffer, CL_FALSE, 0, *size, ptr, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue write buffer!");
}


void FC_FUNC_(ocl_create_command_queue,OCL_CREATE_COMMAND_QUEUE)(bigdft_command_queue *command_queue, cl_context *context){
    cl_int ciErrNum;
    cl_uint device_number;
    *command_queue = (struct _bigdft_command_queue *)malloc(sizeof(struct _bigdft_command_queue));
    if(*command_queue == NULL) {
      fprintf(stderr,"Error: Failed to create command queue (out of memory)!\n");
      exit(1);
    }
#if __OPENCL_VERSION__ <= CL_VERSION_1_0
    size_t nContextDescriptorSize; 
    clGetContextInfo(*context, CL_CONTEXT_DEVICES, 0, 0, &nContextDescriptorSize);
    device_number = nContextDescriptorSize/sizeof(cl_device_id);
#else
    clGetContextInfo(*context, CL_CONTEXT_NUM_DEVICES, sizeof(device_number), &device_number, NULL);
#endif
    cl_device_id * aDevices = (cl_device_id *) malloc(sizeof(cl_device_id)*device_number);
    clGetContextInfo(*context, CL_CONTEXT_DEVICES, sizeof(cl_device_id)*device_number, aDevices, 0);
#if PROFILING
    (*command_queue)->command_queue = clCreateCommandQueue(*context, aDevices[0], CL_QUEUE_PROFILING_ENABLE, &ciErrNum);
#else
    (*command_queue)->command_queue = clCreateCommandQueue(*context, aDevices[0], 0, &ciErrNum);
#endif
    get_device_infos(aDevices[0], &((*command_queue)->device_infos));
    free(aDevices);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, command queue: %p\n",*context, (*command_queue)->command_queue);
#endif
    oclErrorCheck(ciErrNum,"Failed to create command queue!");
    create_kernels(&((*command_queue)->kernels));
}

void FC_FUNC_(ocl_create_command_queue_id,OCL_CREATE_COMMAND_QUEUE_ID)(bigdft_command_queue *command_queue, cl_context *context, cl_uint *index){
    cl_int ciErrNum;
    cl_uint device_number;
    *command_queue = (struct _bigdft_command_queue *)malloc(sizeof(struct _bigdft_command_queue));
    if(*command_queue == NULL) {
      fprintf(stderr,"Error: Failed to create command queue (out of memory)!\n");
      exit(1);
    }
#if __OPENCL_VERSION__ <= CL_VERSION_1_0
    size_t nContextDescriptorSize; 
    clGetContextInfo(*context, CL_CONTEXT_DEVICES, 0, 0, &nContextDescriptorSize);
    device_number = nContextDescriptorSize/sizeof(cl_device_id);
#else
    clGetContextInfo(*context, CL_CONTEXT_NUM_DEVICES, sizeof(device_number), &device_number, NULL);
#endif
    cl_device_id * aDevices = (cl_device_id *) malloc(sizeof(cl_device_id)*device_number);
    clGetContextInfo(*context, CL_CONTEXT_DEVICES, sizeof(cl_device_id)*device_number, aDevices, 0);
#if PROFILING
    (*command_queue)->command_queue = clCreateCommandQueue(*context, aDevices[*index % device_number], CL_QUEUE_PROFILING_ENABLE, &ciErrNum);
#else
    (*command_queue)->command_queue = clCreateCommandQueue(*context, aDevices[*index % device_number], 0, &ciErrNum);
    /*printf("Queue created index : %d, gpu chosen :%d, gpu number : %d\n", *index, *index % device_number, device_number);*/ 
#endif
    get_device_infos(aDevices[*index % device_number], &((*command_queue)->device_infos));
    free(aDevices);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, command queue: %p\n",*context, (*command_queue)->command_queue);
#endif
    oclErrorCheck(ciErrNum,"Failed to create command queue!");
    create_kernels(&((*command_queue)->kernels));
}

size_t shrRoundUp(size_t group_size, size_t global_size)
{
    size_t r = global_size % group_size;
    if(r == 0)
    {
        return global_size;
    } else
    {
        return global_size + group_size - r;
    }
}




void FC_FUNC_(ocl_finish,OCL_FINISH)(bigdft_command_queue *command_queue){
    cl_int ciErrNum;
    ciErrNum = clFinish((*command_queue)->command_queue);
    oclErrorCheck(ciErrNum,"Failed to finish!");
}

void FC_FUNC_(ocl_enqueue_barrier,OCL_ENQUEUE_BARRIER)(bigdft_command_queue *command_queue){
    cl_int ciErrNum;
#ifdef CL_VERSION_1_2
    ciErrNum = clEnqueueBarrierWithWaitList((*command_queue)->command_queue, 0, NULL, NULL);
#else
    ciErrNum = clEnqueueBarrier((*command_queue)->command_queue);
#endif
    oclErrorCheck(ciErrNum,"Failed to enqueue barrier!");
}

void FC_FUNC_(ocl_clean_command_queue,OCL_CLEAN_COMMAND_QUEUE)(bigdft_command_queue *command_queue) {
  clean_magicfilter_kernels(&((*command_queue)->kernels));
  clean_benchmark_kernels(&((*command_queue)->kernels));
  clean_kinetic_kernels(&((*command_queue)->kernels));
  clean_initialize_kernels(&((*command_queue)->kernels));
  clean_wavelet_kernels(&((*command_queue)->kernels));
  clean_uncompress_kernels(&((*command_queue)->kernels));
  clean_reduction_kernels(&((*command_queue)->kernels));
  clean_fft_kernels(&((*command_queue)->kernels));
  clReleaseCommandQueue((*command_queue)->command_queue);
  free(*command_queue);
}

void FC_FUNC_(ocl_clean,OCL_CLEAN)(cl_context *context){
  size_t i;
  clean_magicfilter_programs();
  clean_benchmark_programs();
  clean_kinetic_programs();
  clean_initialize_programs();
  clean_wavelet_programs();
  clean_uncompress_programs();
  clean_reduction_programs();
  clean_fft_programs();
  for(i=0;i<event_number;i++){
    clReleaseEvent(event_list[i].e);
  }
  clReleaseContext(*context);
}
