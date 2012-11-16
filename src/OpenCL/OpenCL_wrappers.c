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

 
struct _opencl_version opencl_version_1_0 = {1,0};
struct _opencl_version opencl_version_1_1 = {1,1};
struct _opencl_version opencl_version_1_2 = {1,2};

cl_int compare_opencl_version(struct _opencl_version v1, struct _opencl_version v2) {
  if(v1.major > v2.major)
    return 1;
  if(v1.major < v2.major)
    return -1;
  if(v1.minor > v2.minor)
    return 1;
  if(v1.minor < v2.minor)
    return -1;
  return 0;
}

static void get_platform_version(cl_platform_id platform_id, struct _opencl_version * version) {
    cl_int ciErrNum = CL_SUCCESS;
    size_t cl_platform_version_size;
    ciErrNum = clGetPlatformInfo(platform_id, CL_PLATFORM_VERSION, 0, NULL, &cl_platform_version_size);
    oclErrorCheck(ciErrNum,"Failed to get CL_PLATFORM_VERSION size!");
    char * cl_platform_version;
    cl_platform_version = (char *)malloc( cl_platform_version_size );
    if(cl_platform_version == NULL) {
      fprintf(stderr,"Error: Failed to create string (out of memory)!\n");
      exit(1);
    }
    ciErrNum = clGetPlatformInfo(platform_id, CL_PLATFORM_VERSION, cl_platform_version_size, cl_platform_version, NULL);
    oclErrorCheck(ciErrNum,"Failed to get CL_PLATFORM_VERSION!");
    //OpenCL<space><major_version.minor_version><space><platform-specific information>
    char minor[2], major[2];
    major[0] = cl_platform_version[7];
    major[1] = 0;
    minor[0] = cl_platform_version[9];
    minor[1] = 0;
    version->major = atoi( major );
    version->major = atoi( minor );
    free(cl_platform_version);
}

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

void FC_FUNC_(ocl_build_programs,OCL_BUILD_PROGRAMS)(bigdft_context * context) {
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
void get_context_devices_infos(bigdft_context * context, struct bigdft_device_infos * infos){
    cl_uint device_number;

#ifdef CL_VERSION_1_1
    if( compare_opencl_version((*context)->PLATFORM_VERSION, opencl_version_1_1) >= 0 )
      clGetContextInfo((*context)->context, CL_CONTEXT_NUM_DEVICES, sizeof(device_number), &device_number, NULL);
    else
#endif
    {
      size_t nContextDescriptorSize;
      clGetContextInfo((*context)->context, CL_CONTEXT_DEVICES, 0, 0, &nContextDescriptorSize);
      device_number = nContextDescriptorSize/sizeof(cl_device_id);
    }
    cl_device_id * aDevices = (cl_device_id *) malloc(sizeof(cl_device_id)*device_number);
    clGetContextInfo((*context)->context, CL_CONTEXT_DEVICES, sizeof(cl_device_id)*device_number, aDevices, 0);

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

void FC_FUNC_(ocl_create_gpu_context,OCL_CREATE_GPU_CONTEXT)(bigdft_context * context,cl_uint *device_number) {
    cl_int ciErrNum = CL_SUCCESS;
    cl_platform_id *platform_ids;
    cl_uint num_platforms;
    clGetPlatformIDs(0, NULL, &num_platforms);
    //printf("num_platforms: %d\n",num_platforms);
       if(num_platforms == 0) {
      fprintf(stderr,"No OpenCL platform available!\n");
      exit(1);
    }
    platform_ids = (cl_platform_id *)malloc(num_platforms * sizeof(cl_platform_id));
    clGetPlatformIDs(num_platforms, platform_ids, NULL);
    cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform_ids[0], 0 };

    *context = (struct _bigdft_context *)malloc(sizeof(struct _bigdft_context));
    if(*context == NULL) {
      fprintf(stderr,"Error: Failed to create context (out of memory)!\n");
      exit(1);
    }
    (*context)->context = clCreateContextFromType( properties , CL_DEVICE_TYPE_GPU, NULL, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p\n",*context);
#endif
    oclErrorCheck(ciErrNum,"Failed to create GPU context!");
    
    get_platform_version(platform_ids[0], &((*context)->PLATFORM_VERSION));
    //getting the number of devices available in the context (devices which are of DEVICE_TYPE_GPU of platform platform_ids[0])
#ifdef CL_VERSION_1_1
    if( compare_opencl_version((*context)->PLATFORM_VERSION, opencl_version_1_1) >= 0 )
      clGetContextInfo((*context)->context, CL_CONTEXT_NUM_DEVICES, sizeof(*device_number), device_number, NULL);
    else
#endif
    {
      size_t nContextDescriptorSize;
      clGetContextInfo((*context)->context, CL_CONTEXT_DEVICES, 0, 0, &nContextDescriptorSize);
      *device_number = nContextDescriptorSize/sizeof(cl_device_id);
    }
    //printf("num_devices: %d\n",*device_number);
}

void FC_FUNC_(ocl_create_cpu_context,OCL_CREATE_CPU_CONTEXT)(bigdft_context * context) {
    cl_int ciErrNum = CL_SUCCESS;
    cl_platform_id platform_id;
    cl_uint platform_number;
    clGetPlatformIDs(1, &platform_id, &platform_number);
    if(platform_number == 0) {
      fprintf(stderr,"No OpenCL platform available!\n");
      exit(1);
    }

    cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform_id, 0 };
    *context = (struct _bigdft_context *)malloc(sizeof(struct _bigdft_context));
    if(*context == NULL) {
      fprintf(stderr,"Error: Failed to create context (out of memory)!\n");
      exit(1);
    }
    (*context)->context = clCreateContextFromType( properties , CL_DEVICE_TYPE_CPU, NULL, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p\n",*context);
#endif
    oclErrorCheck(ciErrNum,"Failed to create CPU context!");

    get_platform_version(platform_id, &((*context)->PLATFORM_VERSION));
}

void FC_FUNC_(ocl_create_read_buffer,OCL_CREATE_READ_BUFFER)(bigdft_context *context, cl_uint *size, cl_mem *buff_ptr) {
    cl_int ciErrNum = CL_SUCCESS;
    *buff_ptr = clCreateBuffer( (*context)->context, CL_MEM_READ_ONLY, *size, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, memory address: %p, size: %lu\n",*context,*buff_ptr,(long unsigned)*size);
#endif
    oclErrorCheck(ciErrNum,"Failed to create read buffer!");
}

void FC_FUNC_(ocl_create_read_write_buffer,OCL_CREATE_READ_WRITE_BUFFER)(bigdft_context *context, cl_uint *size, cl_mem *buff_ptr) {
    cl_int ciErrNum = CL_SUCCESS;
    *buff_ptr = clCreateBuffer( (*context)->context, CL_MEM_READ_WRITE, *size, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, memory address: %p, size: %lu\n",*context,*buff_ptr,(long unsigned)*size);
#endif
    oclErrorCheck(ciErrNum,"Failed to create read_write buffer!");
}

void FC_FUNC_(ocl_create_read_buffer_and_copy,OCL_CREATE_READ_BUFFER_AND_COPY)(bigdft_context *context, cl_uint *size, void *host_ptr, cl_mem *buff_ptr) {
    cl_int ciErrNum = CL_SUCCESS;
    *buff_ptr = clCreateBuffer( (*context)->context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, *size, host_ptr, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create initialized read buffer!");
}

void FC_FUNC_(ocl_create_write_buffer,OCL_CREATE_WRITE_BUFFER)(bigdft_context *context, cl_uint *size, cl_mem *buff_ptr) {
    cl_int ciErrNum = CL_SUCCESS;
    *buff_ptr = clCreateBuffer( (*context)->context, CL_MEM_WRITE_ONLY, *size, NULL, &ciErrNum);
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


void FC_FUNC_(ocl_create_command_queue,OCL_CREATE_COMMAND_QUEUE)(bigdft_command_queue *command_queue, bigdft_context *context){
    cl_int ciErrNum;
    cl_uint device_number;
    *command_queue = (struct _bigdft_command_queue *)malloc(sizeof(struct _bigdft_command_queue));
    if(*command_queue == NULL) {
      fprintf(stderr,"Error: Failed to create command queue (out of memory)!\n");
      exit(1);
    }
#ifdef CL_VERSION_1_1
    if( compare_opencl_version((*context)->PLATFORM_VERSION, opencl_version_1_1) >= 0 )
      clGetContextInfo((*context)->context, CL_CONTEXT_NUM_DEVICES, sizeof(device_number), &device_number, NULL);
    else
#endif
    {
      size_t nContextDescriptorSize; 
      clGetContextInfo((*context)->context, CL_CONTEXT_DEVICES, 0, 0, &nContextDescriptorSize);
      device_number = nContextDescriptorSize/sizeof(cl_device_id);
    }
    cl_device_id * aDevices = (cl_device_id *) malloc(sizeof(cl_device_id)*device_number);
    clGetContextInfo((*context)->context, CL_CONTEXT_DEVICES, sizeof(cl_device_id)*device_number, aDevices, 0);
#if PROFILING
    (*command_queue)->command_queue = clCreateCommandQueue((*context)->context, aDevices[0], CL_QUEUE_PROFILING_ENABLE, &ciErrNum);
#else
    (*command_queue)->command_queue = clCreateCommandQueue((*context)->context, aDevices[0], 0, &ciErrNum);
#endif
    oclErrorCheck(ciErrNum,"Failed to create command queue!");
    (*command_queue)->PLATFORM_VERSION = (*context)->PLATFORM_VERSION;
    get_device_infos(aDevices[0], &((*command_queue)->device_infos));
    free(aDevices);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, command queue: %p\n",*context, (*command_queue)->command_queue);
#endif
    create_kernels(&((*command_queue)->kernels));
}

void FC_FUNC_(ocl_create_command_queue_id,OCL_CREATE_COMMAND_QUEUE_ID)(bigdft_command_queue *command_queue, bigdft_context *context, cl_uint *index){
    cl_int ciErrNum;
    cl_uint device_number;
    *command_queue = (struct _bigdft_command_queue *)malloc(sizeof(struct _bigdft_command_queue));
    if(*command_queue == NULL) {
      fprintf(stderr,"Error: Failed to create command queue (out of memory)!\n");
      exit(1);
    }
#ifdef CL_VERSION_1_1
    if( compare_opencl_version((*context)->PLATFORM_VERSION, opencl_version_1_1) >= 0 )
      clGetContextInfo((*context)->context, CL_CONTEXT_NUM_DEVICES, sizeof(device_number), &device_number, NULL);
    else
#endif
    {
      size_t nContextDescriptorSize; 
      clGetContextInfo((*context)->context, CL_CONTEXT_DEVICES, 0, 0, &nContextDescriptorSize);
      device_number = nContextDescriptorSize/sizeof(cl_device_id);
    }
    cl_device_id * aDevices = (cl_device_id *) malloc(sizeof(cl_device_id)*device_number);
    clGetContextInfo((*context)->context, CL_CONTEXT_DEVICES, sizeof(cl_device_id)*device_number, aDevices, 0);
#if PROFILING
    (*command_queue)->command_queue = clCreateCommandQueue((*context)->context, aDevices[*index % device_number], CL_QUEUE_PROFILING_ENABLE, &ciErrNum);
#else
    (*command_queue)->command_queue = clCreateCommandQueue((*context)->context, aDevices[*index % device_number], 0, &ciErrNum);
    /*printf("Queue created index : %d, gpu chosen :%d, gpu number : %d\n", *index, *index % device_number, device_number);*/ 
#endif
    oclErrorCheck(ciErrNum,"Failed to create command queue!");
    (*command_queue)->PLATFORM_VERSION = (*context)->PLATFORM_VERSION;
    get_device_infos(aDevices[*index % device_number], &((*command_queue)->device_infos));
    free(aDevices);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, command queue: %p\n",*context, (*command_queue)->command_queue);
#endif
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
    if( compare_opencl_version((*command_queue)->PLATFORM_VERSION, opencl_version_1_2) >= 0 )
      ciErrNum = clEnqueueBarrierWithWaitList((*command_queue)->command_queue, 0, NULL, NULL);
    else
#endif
      ciErrNum = clEnqueueBarrier((*command_queue)->command_queue);
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

void FC_FUNC_(ocl_clean,OCL_CLEAN)(bigdft_context *context){
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
  clReleaseContext((*context)->context);
  free(*context);
}
