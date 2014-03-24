#ifndef OPENCL_WRAPPERS_H
#define OPENCL_WRAPPERS_H
#define CL_USE_DEPRECATED_OPENCL_1_0_APIS
#define CL_USE_DEPRECATED_OPENCL_1_1_APIS
#include <CL/opencl.h>
#include <stdio.h>
#include <stdlib.h>
#include <config.h>
#include <math.h>
#include <assert.h>
//#include <Tool.h>
//#include <time.h>

/** @file OpenCL_wrappers.h
 *  @brief Contains global declarations and fortran bindings for OpenCL convolutions.
 *  Warning : every fortran visible procedure is passing its argument per address.
 *  Warning : every floating point data is double precision.
 *  Warning : aliasing of buffer is not supported.
 */

/** Activate debugging info. */
#define DEBUG 0
/** Activate profiling info. */
#define PROFILING 0
#define ERR_STRING(n,s) switch(n)\
    { case  0:  s="CL_SUCCESS"; break;\
      case -1:  s="CL_DEVICE_NOT_FOUND"; break;\
      case -2:  s="CL_DEVICE_NOT_AVAILABLE"; break;\
      case -3:  s="CL_COMPILER_NOT_AVAILABLE"; break;\
      case -4:  s="CL_MEM_OBJECT_ALLOCATION_FAILURE"; break;\
      case -5:  s="CL_OUT_OF_RESOURCES"; break;\
      case -6:  s="CL_OUT_OF_HOST_MEMORY"; break;\
      case -7:  s="CL_PROFILING_INFO_NOT_AVAILABLE"; break;\
      case -8:  s="CL_MEM_COPY_OVERLAP"; break;\
      case -9:  s="CL_IMAGE_FORMAT_MISMATCH"; break;\
      case -10: s="CL_IMAGE_FORMAT_NOT_SUPPORTED"; break;\
      case -11: s="CL_BUILD_PROGRAM_FAILURE"; break;\
      case -12: s="CL_MAP_FAILURE"; break;\
      case -13: s="CL_MISALIGNED_SUB_BUFFER_OFFSET"; break;\
      case -14: s="CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST"; break;\
      case -15: s="CL_COMPILE_PROGRAM_FAILURE"; break;\
      case -16: s="CL_LINKER_NOT_AVAILABLE"; break;\
      case -17: s="CL_LINK_PROGRAM_FAILURE"; break;\
      case -18: s="CL_DEVICE_PARTITION_FAILED"; break;\
      case -19: s="CL_KERNEL_ARG_INFO_NOT_AVAILABLE"; break;\
      case -30: s="CL_INVALID_VALUE"; break;\
      case -31: s="CL_INVALID_DEVICE_TYPE"; break;\
      case -32: s="CL_INVALID_PLATFORM"; break;\
      case -33: s="CL_INVALID_DEVICE"; break;\
      case -34: s="CL_INVALID_CONTEXT"; break;\
      case -35: s="CL_INVALID_QUEUE_PROPERTIES"; break;\
      case -36: s="CL_INVALID_COMMAND_QUEUE"; break;\
      case -37: s="CL_INVALID_HOST_PTR"; break;\
      case -38: s="CL_INVALID_MEM_OBJECT"; break;\
      case -39: s="CL_INVALID_IMAGE_FORMAT_DESCRIPTOR"; break;\
      case -40: s="CL_INVALID_IMAGE_SIZE"; break;\
      case -41: s="CL_INVALID_SAMPLER"; break;\
      case -42: s="CL_INVALID_BINARY"; break;\
      case -43: s="CL_INVALID_BUILD_OPTIONS"; break;\
      case -44: s="CL_INVALID_PROGRAM"; break;\
      case -45: s="CL_INVALID_PROGRAM_EXECUTABLE"; break;\
      case -46: s="CL_INVALID_KERNEL_NAME"; break;\
      case -47: s="CL_INVALID_KERNEL_DEFINITION"; break;\
      case -48: s="CL_INVALID_KERNEL"; break;\
      case -49: s="CL_INVALID_ARG_INDEX"; break;\
      case -50: s="CL_INVALID_ARG_VALUE"; break;\
      case -51: s="CL_INVALID_ARG_SIZE"; break;\
      case -52: s="CL_INVALID_KERNEL_ARGS"; break;\
      case -53: s="CL_INVALID_WORK_DIMENSION"; break;\
      case -54: s="CL_INVALID_WORK_GROUP_SIZE"; break;\
      case -55: s="CL_INVALID_WORK_ITEM_SIZE"; break;\
      case -56: s="CL_INVALID_GLOBAL_OFFSET"; break;\
      case -57: s="CL_INVALID_EVENT_WAIT_LIST"; break;\
      case -58: s="CL_INVALID_EVENT"; break;\
      case -59: s="CL_INVALID_OPERATION"; break;\
      case -60: s="CL_INVALID_GL_OBJECT"; break;\
      case -61: s="CL_INVALID_BUFFER_SIZE"; break;\
      case -62: s="CL_INVALID_MIP_LEVEL"; break;\
      case -63: s="CL_INVALID_GLOBAL_WORK_SIZE"; break;\
      case -64: s="CL_INVALID_PROPERTY"; break;\
      case -65: s="CL_INVALID_IMAGE_DESCRIPTOR"; break;\
      case -66: s="CL_INVALID_COMPILER_OPTIONS"; break;\
      case -67: s="CL_INVALID_LINKER_OPTIONS"; break;\
      case -68: s="CL_INVALID_DEVICE_PARTITION_COUNT"; break;\
      case -1001: s="CL_PLATFORM_NOT_FOUND_KHR"; break;\
      default: s="UNKNOWN_ERROR";\
      }

#define oclErrorCheck(errorCode,message) if(errorCode!=CL_SUCCESS) { char *s; ERR_STRING(errorCode,s) ; fprintf(stderr,"Error(%i)(%s)  (%s: %s): %s\n", errorCode, s,__FILE__,__func__,message); fflush(NULL); exit(1);}


struct bigdft_kernels { 
  cl_kernel c_initialize_kernel_d;
  cl_kernel v_initialize_kernel_d;
  cl_kernel p_initialize_kernel_d;
  cl_kernel kinetic1d_kernel_d;
  cl_kernel kinetic1d_f_kernel_d;
  cl_kernel kinetic_k1d_kernel_d;
  cl_kernel kinetic_k1d_kernel_d_2;
  cl_kernel kinetic_k1d_f_kernel_d_2;
  cl_kernel magicfilter1d_kernel_d;
  cl_kernel magicfilter1d_straight_kernel_d;
  cl_kernel magicfilter1d_block_kernel_d;
  cl_kernel magicfilter1d_den_kernel_d;
  cl_kernel magicfilter1d_pot_kernel_d;
  cl_kernel magicfilter1d_t_kernel_d;
  cl_kernel magicfiltershrink1d_kernel_d;
  cl_kernel magicfiltergrow1d_kernel_d;
  cl_kernel magicfiltergrow1d_den_kernel_d;
  cl_kernel magicfiltergrow1d_pot_kernel_d;
  cl_kernel reduction_kernel_d;
  cl_kernel reduction_dot_kernel_d;
  cl_kernel axpy_kernel_d;
  cl_kernel axpy_offset_kernel_d;
  cl_kernel scal_kernel_d;
  cl_kernel copy_kernel_d;
  cl_kernel dot_kernel_d;
  cl_kernel set_kernel_d;
  cl_kernel void_kernel;
  cl_kernel uncompress_coarse_kernel_d;
  cl_kernel uncompress_fine_kernel_d;
  cl_kernel uncompress_scale_coarse_kernel_d;
  cl_kernel uncompress_scale_fine_kernel_d;
  cl_kernel compress_coarse_kernel_d;
  cl_kernel compress_fine_kernel_d;
  cl_kernel compress_scale_coarse_kernel_d;
  cl_kernel compress_scale_fine_kernel_d;
  cl_kernel scale_psi_fine_kernel_d;
  cl_kernel scale_psi_coarse_kernel_d;
  cl_kernel ana1d_kernel_d;
  cl_kernel ana1d_block_kernel_d;
  cl_kernel anashrink1d_kernel_d;
  cl_kernel syn1d_kernel_d;
  cl_kernel syngrow1d_kernel_d;
  cl_kernel gemm_kernel_d;
  cl_kernel gemm_volkov_kernel_d;
  cl_kernel gemm_kernel_d_tb;
  cl_kernel gemm_kernel_d_ta;
  cl_kernel gemm_kernel_d_tatb;
  cl_kernel gemm_block_kernel_d;
  cl_kernel gemmsy_kernel_d;
  cl_kernel gemm_kernel_z;
  cl_kernel gemm_kernel_z_tb;
  cl_kernel gemm_kernel_z_cb;
  cl_kernel gemm_kernel_z_ta;
  cl_kernel gemm_kernel_z_ca;
  cl_kernel gemm_kernel_z_tatb;
  cl_kernel gemm_kernel_z_catb;
  cl_kernel gemm_kernel_z_tacb;
  cl_kernel gemm_kernel_z_cacb;
  cl_kernel benchmark_flops_kernel_d;
  cl_kernel benchmark_mops_kernel_d;
  cl_kernel transpose_kernel_d;
  cl_kernel notranspose_kernel_d;
  cl_kernel fft_kernel_d0_d;
  cl_kernel fft_kernel_d1_d;
  cl_kernel fft_kernel_d2_d;
  cl_kernel fft_kernel_d0_r2c_d;
  cl_kernel fft_kernel_d1_r2c_d;
  cl_kernel fft_kernel_d2_r2c_d;
  cl_kernel fft_kernel_d0_r_d;
  cl_kernel fft_kernel_d1_r_d;
  cl_kernel fft_kernel_d2_r_d;
  cl_kernel fft_kernel_d0_r_c2r_d;
  cl_kernel fft_kernel_d1_r_c2r_d;
  cl_kernel fft_kernel_d2_r_c2r_d;
  cl_kernel fft_kernel_k_d0_d;
  cl_kernel fft_kernel_k_d1_d;
  cl_kernel fft_kernel_k_d2_d;
  cl_kernel fft_kernel_d0_f_d;
  cl_kernel fft_kernel_d1_f_d;
  cl_kernel fft_kernel_d2_f_d;
  cl_kernel fft_kernel_d0_r2c_f_d;
  cl_kernel fft_kernel_d1_r2c_f_d;
  cl_kernel fft_kernel_d2_r2c_f_d;
  cl_kernel fft_kernel_d0_r_f_d;
  cl_kernel fft_kernel_d1_r_f_d;
  cl_kernel fft_kernel_d2_r_f_d;
  cl_kernel fft_kernel_d0_r_c2r_f_d;
  cl_kernel fft_kernel_d1_r_c2r_f_d;
  cl_kernel fft_kernel_d2_r_c2r_f_d;
  cl_kernel fft_kernel_k_d0_f_d;
  cl_kernel fft_kernel_k_d1_f_d;
  cl_kernel fft_kernel_k_d2_f_d;
};

struct bigdft_device_infos {
  cl_device_type DEVICE_TYPE;
  size_t MAX_WORK_GROUP_SIZE;
  cl_ulong LOCAL_MEM_SIZE;
  cl_uint MAX_COMPUTE_UNITS;
  char NAME[1024];
};

struct _opencl_version {
  cl_uint minor;
  cl_uint major;
}; 

/** Structure associating an OpenCL event with a comment, for profiling purpose. */
typedef struct {
	cl_event e;
	char *comment;
} event;

struct _bigdft_context {
  cl_context context;
  struct _opencl_version PLATFORM_VERSION;
  cl_program benchmarkProgram;
  cl_program fftProgramd0;
  cl_program fftProgramd1;
  cl_program fftProgramd2;
  cl_program initializeProgram;
  cl_program kineticProgram;
  cl_program kinetic_kProgram;
  cl_program magicfilterProgram;
  cl_program reductionProgram;
  cl_program dgemmProgram;
  cl_program uncompressProgram;
  cl_program compressProgram;
  cl_program anaProgram;
  cl_program synProgram;
  cl_mem cossind0;
  cl_mem cossind1;
  cl_mem cossind2;
  cl_uint fft_size[3];
  event * event_list;
  size_t event_number;
  size_t event_allocated;
};

typedef struct _bigdft_context * bigdft_context;

struct _bigdft_command_queue {
  bigdft_context context;
  struct bigdft_kernels kernels;
  struct bigdft_device_infos device_infos;
  cl_command_queue command_queue;
  struct _opencl_version PLATFORM_VERSION;
};

extern struct _opencl_version opencl_version_1_0;
extern struct _opencl_version opencl_version_1_1;
extern struct _opencl_version opencl_version_1_2;

typedef struct _bigdft_command_queue * bigdft_command_queue;

void FC_FUNC_(customize_fft,CUSTOMIZE_FFT)(bigdft_context * context, cl_uint *dimensions);

cl_int compare_opencl_version(struct _opencl_version v1, struct _opencl_version v2);
/** Recovers device info used by BigDFT code generator. */
void get_context_devices_infos(bigdft_context * context, struct bigdft_device_infos * infos);
void get_device_infos(cl_device_id device, struct bigdft_device_infos * infos);
/** Creates all bigdft kernels*/
void create_kernels(bigdft_context * context, struct bigdft_kernels *kernels);
/** Creates magicfilter kernels. to be called after building the magicfilter programs. */
void create_magicfilter_kernels(bigdft_context * context, struct bigdft_kernels * kernels);
void create_benchmark_kernels(bigdft_context * context, struct bigdft_kernels * kernels);
void create_kinetic_kernels(bigdft_context * context, struct bigdft_kernels * kernels);
void create_wavelet_kernels(bigdft_context * context, struct bigdft_kernels * kernels);
void create_uncompress_kernels(bigdft_context * context, struct bigdft_kernels * kernels);
void create_initialize_kernels(bigdft_context * context, struct bigdft_kernels * kernels);
void create_reduction_kernels(bigdft_context * context, struct bigdft_kernels * kernels);
void create_fft_kernels(bigdft_context * context, struct bigdft_kernels * kernels);
/** Compiles magicfilter programs in the given context. */
void build_magicfilter_programs(bigdft_context * context);
void build_reduction_programs(bigdft_context * context);
void build_benchmark_programs(bigdft_context * context);
void build_kinetic_programs(bigdft_context * context);
void build_wavelet_programs(bigdft_context * context);
void build_uncompress_programs(bigdft_context * context);
void build_initialize_programs(bigdft_context * context);
void build_fft_programs(bigdft_context * context);
/** Releases magicfilter kernels. */
void clean_magicfilter_kernels(struct bigdft_kernels * kernels);
void clean_benchmark_kernels(struct bigdft_kernels * kernels);
void clean_kinetic_kernels(struct bigdft_kernels * kernels);
void clean_wavelet_kernels(struct bigdft_kernels * kernels);
void clean_uncompress_kernels(struct bigdft_kernels * kernels);
void clean_initialize_kernels(struct bigdft_kernels * kernels);
void clean_reduction_kernels(struct bigdft_kernels * kernels);
void clean_fft_kernels(bigdft_context * context, struct bigdft_kernels * kernels);
/** Releases magicfilter programs. */
void clean_magicfilter_programs(bigdft_context * context);
void clean_benchmark_programs(bigdft_context * context);
void clean_kinetic_programs(bigdft_context * context);
void clean_initialize_programs(bigdft_context * context);
void clean_wavelet_programs(bigdft_context * context);
void clean_uncompress_programs(bigdft_context * context);
void clean_reduction_programs(bigdft_context * context);
void clean_fft_programs(bigdft_context * context);

/** Returns the first device available in a given context. */
cl_device_id oclGetFirstDev(cl_context cxGPUContext);

/** Returns the next integer that is equal or greater than global_size and a multiple of group_size. */
size_t shrRoundUp(size_t group_size, size_t global_size);


/** Adds an event to the event list. */
int addToEventList(bigdft_context * context, event ev);

/** Reads the processor time stamp counter. */
void FC_FUNC_(rdtsc,RDTSC)(cl_ulong * t);
/** Return the real-time clock time in nanosecond since the epoch. */
void FC_FUNC_(nanosec,NANOSEC)(cl_ulong * t);

/** Initializes the event list. For profiling purpose. */
void FC_FUNC_(init_event_list,INIT_EVENT_LIST)(bigdft_context * context);
/** Prints the event list. */
void FC_FUNC_(print_event_list,PRINT_EVENT_LIST)(bigdft_context * context);
/** Buids and create the OpenCL kernel int the given context. */
void FC_FUNC_(ocl_build_programs,OCL_BUILD_PROGRAMS)(bigdft_context * context);
/** Creates a context containing devices of the type specified from the chosen platform*/
void FC_FUNC_(ocl_create_context,OCL_CREATE_CONTEXT)(bigdft_context * context, const char * platform, const char * devices, cl_int *device_type, cl_uint *device_number);
/** Creates a context containing all GPUs from the default platform */
void FC_FUNC_(ocl_create_gpu_context,OCL_CREATE_GPU_CONTEXT)(bigdft_context * context, cl_uint *device_number);
/** Creates a context containing all CPUs from the default platform */
void FC_FUNC_(ocl_create_cpu_context,OCL_CREATE_CPU_CONTEXT)(bigdft_context * context);
/** Creates a OpenCL read only buffer.
 *  @param context where the buffer is created.
 *  @param size of the buffer.
 *  @param buff_ptr return value : a buffer object reference.
 */
void FC_FUNC_(ocl_create_read_buffer,OCL_CREATE_READ_BUFFER)(bigdft_context *context, cl_uint *size, cl_mem *buff_ptr);
/** Creates an OpenCL buffer.
 *  @param context where the buffer is created.
 *  @param size of the buffer.
 *  @param buff_ptr return value : a buffer object reference.
 */
void FC_FUNC_(ocl_create_read_write_buffer,OCL_CREATE_READ_WRITE_BUFFER)(bigdft_context *context, cl_uint *size, cl_mem *buff_ptr);
void FC_FUNC_(ocl_create_read_buffer_and_copy,OCL_CREATE_READ_BUFFER_AND_COPY)(bigdft_context *context, cl_uint *size, void *host_ptr, cl_mem *buff_ptr);
/** Creates a OpenCL write only buffer.
 *  @param context where the buffer is created.
 *  @param size of the buffer.
 *  @param buff_ptr return value : a buffer object reference.
 */
void FC_FUNC_(ocl_create_write_buffer,OCL_CREATE_WRITE_BUFFER)(bigdft_context *context, cl_uint *size, cl_mem *buff_ptr);
/** Releases an OpenCL buffer. */
void FC_FUNC_(ocl_release_mem_object,OCL_RELEASE_MEM_OBJECT)(cl_mem *buff_ptr);
/** Copies data from an OpenCL buffer to Host memory.
 *  @param command_queue a pointer to the command queue used to make the copy.
 *  @param buffer to copy data from.
 *  @param size of the data to copy.
 *  @param host_ptr to copy the data to.
 */
void FC_FUNC_(ocl_enqueue_read_buffer,OCL_ENQUEUE_READ_BUFFER)(bigdft_command_queue *command_queue, cl_mem *buffer, cl_uint *size, void *host_ptr);
/** Copies data from Host memory to an OpenCL buffer.
 *  @param command_queue a pointer to the command queue used to make the copy.
 *  @param buffer to copy data to.
 *  @param size of the data to copy.
 *  @param host_ptr to copy the data from.
 */
void FC_FUNC_(ocl_enqueue_write_buffer,OCL_ENQUEUE_WRITE_BUFFER)(bigdft_command_queue *command_queue, cl_mem *buffer, cl_uint *size, const void *host_ptr);
/** Copies data from an OpenCL buffer to Host memory asynchronously.
 *  @param command_queue a pointer to the command queue used to make the copy.
 *  @param buffer to copy data from.
 *  @param size of the data to copy.
 *  @param host_ptr to copy the data to.
 */
void FC_FUNC_(ocl_enqueue_read_buffer_async,OCL_ENQUEUE_READ_BUFFER_ASYNC)(bigdft_command_queue *command_queue, cl_mem *buffer, cl_uint *size, void *host_ptr);
/** Copies data from Host memory to an OpenCL buffer asynchronously.
 *  @param command_queue a pointer to the command queue used to make the copy.
 *  @param buffer to copy data to.
 *  @param size of the data to copy.
 *  @param host_ptr to copy the data from.
 */
void FC_FUNC_(ocl_enqueue_write_buffer_async,OCL_ENQUEUE_WRITE_BUFFER_ASYNC)(bigdft_command_queue *command_queue, cl_mem *buffer, cl_uint *size, const void *host_ptr);
/** Creates a command queue in the given context, associating it to the first device in the context */
void FC_FUNC_(ocl_create_command_queue,OCL_CREATE_COMMAND_QUEUE)(bigdft_command_queue *command_queue, bigdft_context *context);
/** Creates a command queue in the given context, associating it to the device specified by index modulo the number of device. */
void FC_FUNC_(ocl_create_command_queue_id,OCL_CREATE_COMMAND_QUEUE_ID)(bigdft_command_queue *command_queue, bigdft_context *context, cl_uint *index);
/** Waits for all commands in a queue to complete. */
void FC_FUNC_(ocl_finish,OCL_FINISH)(bigdft_command_queue *command_queue);
/** Enqueues a barrier in a queue. Commands enqueued after the barrier will wait
 *  for commands enqueued before the barrier to be processed before being sent to the device. */
void FC_FUNC_(ocl_enqueue_barrier,OCL_ENQUEUE_BARRIER)(bigdft_command_queue *command_queue);
/** Releases a command queue and the associated kernels. */
void FC_FUNC_(ocl_clean_command_queue,OCL_CLEAN_COMMAND_QUEUE)(bigdft_command_queue *command_queue);
/** Releases the context, and beforehand releases the programs. */
void FC_FUNC_(ocl_clean,OCL_CLEAN)(bigdft_context *context);

/** Performs the one dimensional wavelet analysis and transposition with periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * (2 * n) * sizeof(double), stored in column major order.
 *  @param out output buffer of size (2 * n) * ndat * sizeof(double), stored in column major order.
 */
void FC_FUNC_(ana1d_d,ANA1D_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out);
/** Performs the one dimensional wavelet analysis and transposition with open boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * (2 * n + 14) * sizeof(double), stored in column major order.
 *  @param out output buffer of size (2 * n) * ndat * sizeof(double), stored in column major order.
 */
void FC_FUNC_(anashrink1d_d,ANASHRINK1D_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out);
/** Slightly more performing version of anashrink1d_d. @see anashrink1d_d. */
void FC_FUNC_(ana1d_block_d,ANA1D_BLOCK_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
/** Performs the three-dimensional wavelet analysis with periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param tmp temporary buffer to store intermediate results. Dimensions : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double).
 *  @param psi input buffer of dimension : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 *  @param out output buffer of dimensions : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 */
void FC_FUNC_(ana_d,ANA_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out);
/** Slightly more performing version of ana_d. @see ana_d. */
void FC_FUNC_(ana_block_d,ANA_BLOCK_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out);
/** Performs the three-dimensional wavelet analysis with periodic or non periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param periodic periodicity of the convolution. Vector of three value, one for each dimension. Non zero means periodic.
 *  @param tmp temporary buffer to store intermediate results. Must be of at least (2 * dimensions[0] + (periodic[0]?0:14)) * (2 * dimensions[1] + (periodic[1]?0:14)) * (2 * dimensions[2]) * sizeof(double) in size.
 *  @param psi input buffer of dimension : (2 * dimensions[0] + (periodic[0]?0:14)) * (2 * dimensions[1] + (periodic[1]?0:14)) * (2 * dimensions[2] + (periodic[2]?0:14)) * sizeof(double). Stored in collumn major order.
 *  @param out output buffer of dimensions : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 */
void FC_FUNC_(ana_d_generic,ANA_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *psi, cl_mem *out);
/** Version of ana_d without the temporary buffer, psi is erased during the computation. @see ana_d. */
void FC_FUNC_(ana_self_d,ANA_SELF_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *psi, cl_mem *out);
/** Version of ana_sef_d without the temporary buffer, psi is erased during the computation. @see ana_d_generic. */
void FC_FUNC_(ana_self_d_generic,ANA_SELF_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *psi, cl_mem *out);
/** Performs the one dimensional wavelet synthesis and transposition with periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * (2 * n) * sizeof(double), stored in column major order.
 *  @param out output buffer of size (2 * n) * ndat * sizeof(double), stored in column major order.
 */
void FC_FUNC_(syn1d_d,SYN1D_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out);
/** Performs the one dimensional wavelet analysis and transposition with open boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * (2 * n) * sizeof(double), stored in column major order.
 *  @param out output buffer of size (2 * n + 14) * ndat * sizeof(double), stored in column major order.
 */
void FC_FUNC_(syngrow1d_d,SYNGROW1D_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out);
/** Performs the three-dimensional wavelet synthesis with periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param tmp temporary buffer to store intermediate results. Dimensions : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double).
 *  @param psi input buffer of dimension : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 *  @param out output buffer of dimensions : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 */
void FC_FUNC_(syn_d,SYN_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out);
/** Performs the three-dimensional wavelet analysis with periodic or non periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param periodic periodicity of the convolution. Vector of three value, one for each dimension. Non zero means periodic.
 *  @param tmp temporary buffer to store intermediate results. Must be of at least (2 * dimensions[0]) * (2 * dimensions[1] + (periodic[1]?0:14)) * (2 * dimensions[2] + (periodic[1]?0:14)) * sizeof(double) in size.
 *  @param psi input buffer of dimension : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 *  @param out output buffer of dimensions : (2 * dimensions[0] + (periodic[0]?0:14)) * (2 * dimensions[1] + (periodic[1]?0:14)) * (2 * dimensions[2] + (periodic[2]?0:14)) * sizeof(double). Stored in column major order.
 */
void FC_FUNC_(syn_d_generic,SYN_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *psi, cl_mem *out);
/** Version of syn_d without the temporary buffer, psi is erased during the computation. @see syn_d. */
void FC_FUNC_(syn_self_d,SYN_SELF_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *psi, cl_mem *out);
/** Version of syn_sef_d without the temporary buffer, psi is erased during the computation. @see syn_d_generic. */
void FC_FUNC_(syn_self_d_generic,SYN_SELF_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *psi, cl_mem *out);

/** Performs the one dimensional magicfilter and transposition with periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * n * sizeof(double), stored in collumn major order.
 *  @param out output buffer of size n * ndat * sizeof(double), stored in collumn major order.
 */
void FC_FUNC_(magicfilter1d_d,MAGICFILTER1D_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
/** Performs the one dimensional magicfilter and transposition with open boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * n * sizeof(double), stored in column major order.
 *  @param out output buffer of size (n + 15) * ndat * sizeof(double), stored in column major order.
 */
void FC_FUNC_(magicfiltergrow1d_d,MAGICFILTERGROW1D_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
/** Performs the one dimensional reciprocal magicfilter and transposition with periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * n * sizeof(double), stored in collumn major order.
 *  @param out output buffer of size n * ndat * sizeof(double), stored in collumn major order.
 */
void FC_FUNC_(magicfilter1d_t_d,MAGICFILTER1D_T_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
/** Performs the one dimensional reciprocal magicfilter and transposition with open boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * (n +15) * sizeof(double), stored in column major order.
 *  @param out output buffer of size n * ndat * sizeof(double), stored in column major order.
 */
void FC_FUNC_(magicfiltershrink1d_d,MAGICFILTERSHRINK1D_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
/** Performs the one dimensional magicfilter and transposition with periodic boundary conditions. Storage of matrix is changed with respect to magicfilter1d_d.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size n * ndat * sizeof(double), stored in collumn major order.
 *  @param out output buffer of size ndat * n * sizeof(double), stored in collumn major order.
 */
void FC_FUNC_(magicfilter1d_straight_d,MAGICFILTER1D_STRAIGHT_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
/** Slightly more performing version of magicfilter1d_d. @see magicfilter1d_d. */
void FC_FUNC_(magicfilter1d_block_d,MAGICFILTER1D_BLOCK_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
/** Performs the one dimensional magicfilter and transposition with periodic boundary conditions and multiplies by a potential.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * n * sizeof(double), stored in collumn major order.
 *  @param pot potential applied during the convolution. Size is n * ndat * sizeof(double), stored in collumn major order.
 *  @param out output buffer of size n * ndat * sizeof(double), stored in collumn major order.
 */
void FC_FUNC_(magicfilter1d_pot_d,MAGICFILTER1D_POT_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *pot, cl_mem *out);
/** Performs the three-dimensional magicfilter with periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param tmp temporary buffer to store intermediate results. Must be of at least dimensions[0] * dimensions[1] * dimensions[2] * sizeof(double) in size.
 *  @param psi input buffer of dimension : dimensions[0] * dimensions[1] * dimensions[2] * sizeof(double). Stored in column major order.
 *  @param out output buffer of dimensions : dimensions[0] * dimensions[1] * dimensions[2] * sizeof(double). Stored in column major order.
 */
void FC_FUNC_(magicfilter_n_d,MAGICFILTER_N_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out);
/** Version of magicfilter_n_d without the temporary buffer, psi is erased during the computation. @see magicfilter_n_d. */
void FC_FUNC_(magicfilter_n_self_d,MAGICFILTER_N_SELF_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *psi, cl_mem *out);
/** Slightly more performing version of magicfilter_n_d. @see magicfilter_n_d. */
void FC_FUNC_(magicfilter_n_straight_d,MAGICFILTER_N_STRAIGHT_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out);
/** Slightly more performing version of magicfilter_n_d. @see magicfilter_n_d. */
void FC_FUNC_(magicfilter_n_block_d,MAGICFILTER_N_BLOCK_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out);
/** Performs the three-dimensional magicfilter with periodic boundary conditions, and squares the values to compute the density.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param tmp temporary buffer to store intermediate results. Must be of at least dimensions[0] * dimensions[1] * dimensions[2] * sizeof(double) in size.
 *  @param psi input buffer of dimension : dimensions[0] * dimensions[1] * dimensions[2] * sizeof(double). Stored in column major order.
 *  @param out output buffer of dimensions : dimensions[0] * dimensions[1] * dimensions[2] * sizeof(double). Stored in column major order.
 */
void FC_FUNC_(magicfilter_den_d,MAGICFILTER_DEN_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out);
/** Performs the three-dimensional magicfilter with periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param tmp temporary buffer to store intermediate results. Must be of at least dimensions[0] * dimensions[1] * dimensions[2] * sizeof(double) in size.
 *  @param psi input buffer of dimension : dimensions[0] * dimensions[1] * dimensions[2] * sizeof(double). Stored in column major order.
 *  @param out output buffer of dimensions : dimensions[0] * dimensions[1] * dimensions[2] * sizeof(double). Stored in column major order.
 */
void FC_FUNC_(magicfilter_den_d_generic,MAGICFILTER_DEN_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *psi, cl_mem *out);
/** Performs the three-dimensional magicfilter with periodic or non-periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param periodic periodicity of the convolution. Vector of three value, one for each dimension. Non zero means periodic.
 *  @param tmp temporary buffer to store intermediate results. Must be of at least (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double) in size.
 *  @param psi input buffer of dimension : (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double). Stored in column major order.
 *  @param out output buffer of dimensions : (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double). Stored in column major order.
 */
void FC_FUNC_(magicfilter_t_d,MAGICFILTER_T_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out);
/** Version of magicfilter_t_d without the temporary buffer, psi is erased during the computation. @see magicfilter_t_d. */
void FC_FUNC_(magicfilter_t_self_d,MAGICFILTER_T_SELF_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *psi, cl_mem *out);
/** Performs the three-dimensional magicfilter, applies the potential then applies the reciprocal three dimension magicfilter. With periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param tmp temporary buffer to store intermediate results. Must be of at least dimensions[0] * dimensions[1] * dimensions[2] * sizeof(double) in size.
 *  @param psi input buffer of dimension : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 *  @param out output buffer of dimensions : (2 * dimensions[0]) * ( 2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 *  @param pot potential applied, buffer of dimensions : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 */
void FC_FUNC_(potential_application_d,POTENTIAL_APPLICATION_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out, cl_mem *pot);
/** Performs the three-dimensional magicfilter, applies the potential then applies the reciprocal three dimension magicfilter. The potential energy is also computed. With periodic or non periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param periodic periodicity of the convolution. Vector of three value, one for each dimension. Non zero means periodic.
 *  @param tmp temporary buffer to store intermediate results. Must be of at least (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double) in size.
 *  @param tmp_dot temporary buffer to store intermediate results. Must be of at least (2 * dimensions[0] + (periodic[0]?0:14)) * (2 * dimensions[1] + (periodic[1]?0:14)) * (2 * dimensions[2] + (periodic[2]?0:14)) * sizeof(double) in size.
 *  @param psi input buffer of dimension : (2 * dimensions[0] + (periodic[0]?0:14)) * (2 * dimensions[1] + (periodic[1]?0:14)) * (2 * dimensions[2] + (periodic[2]?0:14)) * sizeof(double). Stored in column major order.
 *  @param out output buffer of dimensions : (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double). Stored in column major order.
 *  @param pot potential applied, buffer of dimensions : (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double). Stored in column major order.
 *  @param epot potential energy computed.
 */
void FC_FUNC_(potential_application_d_generic,POTENTIAL_APPLICATION_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *tmp_dot, cl_mem *psi, cl_mem *out, cl_mem *pot, cl_double *epot);

/** Benchmark to evaluate the throughput of the copy mechanism used in the convolutions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the first dimension.
 *  @param ndat size of the second dimension.
 *  @param psi input buffer of size ndat * n * sizeof(double), stored in column major order.
 *  @param out output buffer of size ndat * n * sizeof(double), stored in column major order.
 */
void FC_FUNC_(notranspose_d,NOTRANSPOSE_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
/** Benchmark to evaluate the throughput of the transposition mechanism used in the convolutions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the first dimension.
 *  @param ndat size of the second dimension.
 *  @param psi input buffer of size ndat * n * sizeof(double), stored in column major order.
 *  @param out output buffer of size n * ndat * sizeof(double), stored in column major order.
 */
void FC_FUNC_(transpose_d,TRANSPOSE_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
/** Benchmark to evaluate the throughput of the OpenCL device in FLOPS, each element processed generates 4096 FLOP.
 *  @param command_queue used to process the convolution.
 *  @param n number of elements.
 *  @param in input buffer of size n * sizeof(double), stored in column major order.
 *  @param out output buffer of size n * sizeof(double), stored in column major order.
 */
void FC_FUNC_(benchmark_flops_d,BENCHMARK_FLOPS_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *out);
/** Benchmark to evaluate the throughput of the OpenCL device global memory in MOPS, each element processed generates 8 read and 8 write, which are coalesced.
 *  @param command_queue used to process the convolution.
 *  @param n number of elements.
 *  @param in input buffer of size n * sizeof(double), stored in column major order.
 *  @param out output buffer of size n * sizeof(double), stored in column major order.
 */
void FC_FUNC_(benchmark_mops_d,BENCHMARK_MOPS_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *out);

/** Performs the one dimensional kinetic filter and transposition with periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param h hgrid along the dimension processed.
 *  @param c scaling factor.
 *  @param x input buffer of size ndat * n * sizeof(double), stored in column major order.
 *  @param y output buffer of size n * ndat * sizeof(double), stored in column major order.
 *  @param workx output buffer of size n * ndat * sizeof(double), stored in column major order. Transposition of x.
 *  @param work_y temporary buffer used to store intermediate results. Size ndat * n * sizeof(double), stored in column major order.
 *  @param ekin dummy argument. Will be used to compute the kinetic energy.
 */
void FC_FUNC_(kinetic1d_d,KINETIC1D_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_double *h, cl_double *c, cl_mem *x, cl_mem *y, cl_mem *workx, cl_mem *worky, cl_double *ekin);
/** Performs the three dimensional Kinetic filter with periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param h hgrid along the three dimensions. Vector of three values.
 *  @param x input buffer of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 *  @param y input buffer of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 *  @param work_x work buffer of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 *  @param work_y output buffer of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order. work_y = y + kinetic(x).
 */
void FC_FUNC_(kinetic_d,KINETIC_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y);
/** Performs the three dimensional Kinetic filter with periodic or non periodic boundary conditions. Input arrays are lost.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param h hgrid along the three dimensions. Vector of three values.
 *  @param x input buffer of size (2 * dimensions[0] + (periodic[0]?0:14)) * (2 * dimensions[1] + (periodic[1]?0:14)) * (2 * dimensions[2] + (periodic[2]?0:14)) * sizeof(double). Stored in column major order.
 *  @param y input buffer of size (2 * dimensions[0] + (periodic[0]?0:14)) * (2 * dimensions[1] + (periodic[1]?0:14)) * (2 * dimensions[2] + (periodic[2]?0:14)) * sizeof(double). Stored in column major order.
 *  @param work_x work buffer of size (2 * dimensions[0] + (periodic[0]?0:14)) * (2 * dimensions[1] + (periodic[1]?0:14)) * (2 * dimensions[2] + (periodic[2]?0:14)) * sizeof(double). Stored in column major order.
 *  @param work_y output buffer of size (2 * dimensions[0] + (periodic[0]?0:14)) * (2 * dimensions[1] + (periodic[1]?0:14)) * (2 * dimensions[2] + (periodic[2]?0:14)) * sizeof(double). Stored in column major order. work_y = y + kinetic(x).
 */
void FC_FUNC_(kinetic_d_generic,KINETIC_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y);
/** Version of kinetic_d using two temporary buffer to avoid erasing input arrays. @see kinetic_d. */
void FC_FUNC_(kinetic_stable_d,KINETIC_STABLE_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y, cl_mem *tmp_x, cl_mem *tmp_y);
/** Performs the three dimensional Kinetic filter with periodic boundary conditions on K point data.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param h hgrid along the three dimensions. Vector of three values.
 *  @param x input buffer of size dimensions[0] * dimensions[1] * dimensions[2] * sizeof(complex double). Stored in column major order.
 *  @param y input and output buffer of size dimensions[0] * dimensions[1] * dimensions[2] * sizeof(complex double). Stored in column major order.
 *  @param work_x work buffer of size dimensions[0] * dimensions[1] * dimensions[2] * sizeof(complex double). Stored in column major order.
 *  @param work_y work buffer of size dimensions[0] * dimensions[1] * dimensions[2] * sizeof(complex double). Stored in column major order. work_y = y + kinetic(x).
 *  @param c_in constant affecting the scaling factor.
 *  @param k point coordinates. Verctor of three values.
 */
void FC_FUNC_(kinetic_k_d,KINETIC_K_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y, cl_double * c_in,  cl_double *k);
void FC_FUNC_(kinetic_k_d_generic,KINETIC_K_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, double *h, double *k, cl_mem *x_r, cl_mem *x_i, cl_mem *y_r, cl_mem *y_i, cl_mem *work_x_r, cl_mem *work_x_i, cl_mem *work_y_r, cl_mem *work_y_i);
/** Computes the sum of the components of a vector.
 *  @param command_queue used to process the data.
 *  @param n number of element of the vector.
 *  @param in buffer of size n * sizeof(double), containing the input data.
 *  @param work1 temporary buffer of size n * sizeof(double).
 *  @param work2 temporary buffer of size n * sizeof(double).
 *  @param out pointer to the resulting sum.
 */
void FC_FUNC_(asum_d,ASUM_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *work1, cl_mem *work2, cl_double *out);
/** Version of asum_d without the second temporary buffer, erasing the input. @see asum_d. */
void FC_FUNC_(asum_self_d,ASUM_SELF_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *work, cl_double *out);
/** Computes the squared norm 2 of a vector.
 *  @param command_queue used to process the data.
 *  @param n number of element of the vector.
 *  @param in buffer of size n * sizeof(double), containing the input data.
 *  @param work1 temporary buffer of size n * sizeof(double).
 *  @param work2 temporary buffer of size n * sizeof(double).
 *  @param out pointer to the resulting sum.
 */
void FC_FUNC_(nrm2sq_d,NRM2SQ_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *work1, cl_mem *work2, cl_double *out);
/** Version of nrm2sq_d without the second temporary buffer, erasing the input. @see nrm2sq_d. */
void FC_FUNC_(nrm2sq_self_d,NRM2SQ_SELF_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *work, cl_double *out);
/** Computes z = alpha * x + y for vectors.
 *  @param command_queue used to process the data.
 *  @param n number of element of the vectors.
 *  @param alpha scaling coefficient.
 *  @param x input buffer, vector of size n * sizeof(double).
 *  @param y input buffer, vector of size n * sizeof(double).
 *  @param z output buffer, vector of size n * sizeof(double).
 */
void FC_FUNC_(axpy_d,AXPY_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *x, cl_mem *y, cl_mem *z);
/** Computes y = alpha * x + y for vectors.
 *  @param command_queue used to process the data.
 *  @param n number of element of the vectors.
 *  @param alpha scaling coefficient.
 *  @param x input buffer, vector of size n * sizeof(double).
 *  @param y input and output buffer, vector of size n * sizeof(double).
 */
void FC_FUNC_(axpy_self_d,AXPY_SELF_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *x, cl_mem *y);
/** Computes y = alpha * x for vectors.
 *  @param command_queue used to process the data.
 *  @param n number of element of the vectors.
 *  @param alpha scaling coefficient.
 *  @param x input buffer, vector of size n * sizeof(double).
 *  @param y output buffer, vector of size n * sizeof(double).
 */
void FC_FUNC_(scal_d,SCAL_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *x, cl_mem *y);
/** Computes x = alpha * x for vectors.
 *  @param command_queue used to process the data.
 *  @param n number of element of the vectors.
 *  @param alpha scaling coefficient.
 *  @param x input and output buffer, vector of size ndat * sizeof(double).
 */
void FC_FUNC_(scal_self_d,SCAL_SELF_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *x);
/** Computes the dot product of 2 vectors.
 *  @param command_queue used to process the data.
 *  @param n number of element of the vectors.
 *  @param x input buffer, vector of size n * sizeof(double).
 *  @param y input buffer, vector of size n * sizeof(double).
 *  @param work1 temporary buffer of size n * sizeof(double).
 *  @param work2 temporary buffer of size n * sizeof(double).
 *  @param out result.
 */
void FC_FUNC_(dot_d,DOT_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_mem *x, cl_mem *y, cl_mem *work1, cl_mem *work2, cl_double *out);
/** Computes the dot product of 2 vectors asynchronously.
 *  @param command_queue used to process the data.
 *  @param n number of element of the vectors.
 *  @param x input buffer, vector of size n * sizeof(double).
 *  @param y input buffer, vector of size n * sizeof(double).
 *  @param work1 temporary buffer of size n * sizeof(double).
 *  @param work2 temporary buffer of size n * sizeof(double).
 *  @param out result.
 */
void FC_FUNC_(dot_d_async,DOT_D_ASYNC)(bigdft_command_queue *command_queue, cl_uint *n, cl_mem *x, cl_mem *y, cl_mem *work1, cl_mem *work2, cl_double *out);
/** Initializes every component of a vector to a given value.
 *  @param command_queue used to process the data.
 *  @param n number of element of the vectors.
 *  @param value used to initialize the vector.
 *  @param x buffer to initialize, of size n * sizeof(double).
 */
void FC_FUNC_(set_d,SET_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_double *value, cl_mem *x);
/** Computex y = x for vectors.
 *  @param command_queue used to process the data.
 *  @param n number of element of the vectors.
 *  @param x input buffer, vector of size n * sizeof(double).
 *  @param y output buffer, vector of size n * sizeof(double).
 */
void FC_FUNC_(copy_d,COPY_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_mem *x, cl_mem *y);
/** Computes z[offset_z:offset_z+n) = alpha * x[offset_x:offset_x+n) + y[offset_y:offset_y+x)..
 *  @param command_queue used to process the data.
 *  @param n number of element to process.
 *  @param alpha scaling coefficient.
 *  @param offset_x offset in the vector x.
 *  @param x input buffer, vector of size n * sizeof(double).
 *  @param offset_y offset in the vector y.
 *  @param y input buffer, vector of size n * sizeof(double).
 *  @param offset_z offset in the vector z.
 *  @param z output buffer, vector of size n * sizeof(double).
 */
void FC_FUNC_(axpy_offset_d,AXPY_OFFSET_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_double *alpha,
                                                                            cl_uint *offset_x, cl_mem *x,
                                                                            cl_uint *offset_y, cl_mem *y,
                                                                            cl_uint *offset_z, cl_mem *z);
/** Computes y[offset_y:offset_y+n) = alpha * x[offset_x:offset_x+n) + y[offset_y:offset_y+x)..
 *  @param command_queue used to process the data.
 *  @param n number of element to process.
 *  @param alpha scaling coefficient.
 *  @param offset_x offset in the vector x.
 *  @param x input buffer, vector of size n * sizeof(double).
 *  @param offset_y offset in the vector y.
 *  @param y input buffer, vector of size n * sizeof(double).
 */
void FC_FUNC_(axpy_offset_self_d,AXPY_OFFSET_SELF_D)(bigdft_command_queue *command_queue, cl_uint *n, cl_double *alpha,
                                                                            cl_uint *offset_x, cl_mem *x,
                                                                            cl_uint *offset_y, cl_mem *y);
/** Computes the multiplication of 2 matrix.
 *  Usage is identical to the BLAS routine DGEMM.
 *  Refer to the BLAS documentation for the meaning of the
 *  different parameters, in respect of the different
 *  transposition and conjugation possible.
 *  @param command_queue used to process the data.
 */
void FC_FUNC_(gemm_d,GEMM_D)(bigdft_command_queue *command_queue, char *transa, char *transb, cl_uint *m, cl_uint *n, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc);
void FC_FUNC_(gemm_volkov_d,GEMM_VOLKOV_D)(bigdft_command_queue *command_queue, char *transa, char *transb, cl_uint *m, cl_uint *n, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc);
/** Computes the multiplication of 2 matrix, knowing the result to be a symmetric matrix.
 *  m and n have to be identical.
 *  Usage is identical to the BLAS routine DGEMM.
 *  Refer to the BLAS documentation for the meaning of the
 *  different parameters, in respect of the different
 *  transposition and conjugation possible.
 *  @param command_queue used to process the data.
 */
void FC_FUNC_(gemmsy_d,GEMMSY_D)(bigdft_command_queue *command_queue, char *transa, char *transb, cl_uint *m, cl_uint *n, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc);
/** Slightly more performing version of gemm_d. @see gemm_d. */
void FC_FUNC_(gemm_block_d,GEMM_BLOCK_D)(bigdft_command_queue *command_queue, char *transa, char *transb, cl_uint *m, cl_uint *n, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc);
/** Computes the multiplication of 2 complex matrix.
 *  Usage is identical to the BLAS routine ZGEMM.
 *  Refer to the BLAS documentation for the meaning of the
 *  different parameters, in respect of the different
 *  transposition and conjugation possible.
 *  @param command_queue used to process the data.
 */
void FC_FUNC_(gemm_z,GEMM_Z)(bigdft_command_queue *command_queue, char *transa, char *transb, cl_uint *m, cl_uint *n, cl_uint *k, cl_double2 *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double2 *beta, cl_mem *c, cl_uint *ldc);

/** Uncompresses a wave function using BigDFT sparse wave function compression.
 *  @param command_queue used to process the data.
 *  @param dimensions of the output data, vector of 3 values.
 *  @param nseg_c number of segment of coarse data.
 *  @param nvctr_c number of point of coarse data.
 *  @param keyg_c array of size 2 * nseg_c * sizeof(uint), representing the beginning and end of coarse segments in the output data.
 *  @param keyv_c array of size nseg_c * sizeof(uint), representing the beginning of coarse segments in the input data.
 *  @param nseg_f number of segment of fine data.
 *  @param nvctr_f number of point of fine data.
 *  @param keyg_f array of size 2 * nseg_f * sizeof(uint), representing the beginning and end of fine segments in the output data.
 *  @param keyv_f array of size nseg_f * sizeof(uint), representing the beginning of fine segments in the input data.
 *  @param psi_c array of size nvctr_c * sizeof(double), containing coarse input data.
 *  @param psi_f array of size nvctr_f * 7 * sizeof(double), containing fine input data.
 *  @param psi_out array of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) containing output data.
 */
void FC_FUNC_(uncompress_d,UNCOMPRESS_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,
                                       cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                       cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                       cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi_out);
/** Compresses a wave function using BigDFT sparse wave function compression.
 *  @param command_queue used to process the data.
 *  @param dimensions of the input data, vector of 3 values.
 *  @param nseg_c number of segment of coarse data.
 *  @param nvctr_c number of point of coarse data.
 *  @param keyg_c array of size 2 * nseg_c * sizeof(uint), representing the beginning and end of coarse segments in the input data.
 *  @param keyv_c array of size nseg_c * sizeof(uint), representing the beginning of coarse segments in the output data.
 *  @param nseg_f number of segment of fine data.
 *  @param nvctr_f number of point of fine data.
 *  @param keyg_f array of size 2 * nseg_f * sizeof(uint), representing the beginning and end of fine segments in the input data.
 *  @param keyv_f array of size nseg_f * sizeof(uint), representing the beginning of fine segments in the output data.
 *  @param psi_c array of size nvctr_c * sizeof(double), containing coarse output data.
 *  @param psi_f array of size nvctr_f * 7 * sizeof(double), containing fine output data.
 *  @param psi_out array of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) containing input data.
 */
void FC_FUNC_(compress_d,COMPRESS_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,
                                     cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                     cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                     cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi);
/** Scales the wavefunction in compressed form.
 *  @param command_queue used to process the data.
 *  @param nvctr_c number of point of coarse data.
 *  @param nvctr_f number of point of fine data.
 *  @param h hgrid of the system, vector of 3 values.
 *  @param c constant part of the scaling.
 *  @param psi_c array of size nvctr_c * sizeof(double), containing coarse input and output data.
 *  @param psi_f array of size nvctr_f * 7 * sizeof(double), containing fine input and output data.
 */
void FC_FUNC_(scale_psi_d,SCALE_PSI_D)(bigdft_command_queue *command_queue, cl_uint *nvctr_c, cl_uint *nvctr_f, cl_double *h, cl_double *c, cl_mem *psi_c,  cl_mem *psi_f);
/** Uncompresses and scales a wave function using BigDFT sparse wave function compression.
 *  @param command_queue used to process the data.
 *  @param dimensions of the output data, vector of 3 values.
 *  @param h hgrid of the system, vector of 3 values.
 *  @param c constant part of the scaling.
 *  @param nseg_c number of segment of coarse data.
 *  @param nvctr_c number of point of coarse data.
 *  @param keyg_c array of size 2 * nseg_c * sizeof(uint), representing the beginning and end of coarse segments in the output data.
 *  @param keyv_c array of size nseg_c * sizeof(uint), representing the beginning of coarse segments in the input data.
 *  @param nseg_f number of segment of fine data.
 *  @param nvctr_f number of point of fine data.
 *  @param keyg_f array of size 2 * nseg_f * sizeof(uint), representing the beginning and end of fine segments in the output data.
 *  @param keyv_f array of size nseg_f * sizeof(uint), representing the beginning of fine segments in the input data.
 *  @param psi_c array of size nvctr_c * sizeof(double), containing coarse input data.
 *  @param psi_f array of size nvctr_f * 7 * sizeof(double), containing fine input data.
 *  @param psi_out array of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) containing output data.
 */
void FC_FUNC_(uncompress_scale_d,UNCOMPRESS_SCALE_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_double *h, cl_double *c,
                                       cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                       cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                       cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi_out);
/** Compresses and scales a wave function using BigDFT sparse wave function compression.
 *  @param command_queue used to process the data.
 *  @param dimensions of the input data, vector of 3 values.
 *  @param h hgrid of the system, vector of 3 values.
 *  @param c constant part of the scaling.
 *  @param nseg_c number of segment of coarse data.
 *  @param nvctr_c number of point of coarse data.
 *  @param keyg_c array of size 2 * nseg_c * sizeof(uint), representing the beginning and end of coarse segments in the input data.
 *  @param keyv_c array of size nseg_c * sizeof(uint), representing the beginning of coarse segments in the output data.
 *  @param nseg_f number of segment of fine data.
 *  @param nvctr_f number of point of fine data.
 *  @param keyg_f array of size 2 * nseg_f * sizeof(uint), representing the beginning and end of fine segments in the input data.
 *  @param keyv_f array of size nseg_f * sizeof(uint), representing the beginning of fine segments in the output data.
 *  @param psi_c array of size nvctr_c * sizeof(double), containing coarse output data.
 *  @param psi_f array of size nvctr_f * 7 * sizeof(double), containing fine output data.
 *  @param psi_out array of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) containing input data.
 */
void FC_FUNC_(compress_scale_d,COMPRESS_SCALE_D)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_double *h, cl_double *c,
                                     cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                     cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                     cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi);

/** Performs a full local hamiltonian on a compressed wave function.
 *  @param command_queue used to process the data.
 *  @param dimensions of the input data, vector of 3 values.
 *  @param periodic periodicity of the convolution. Vector of three value, one for each dimension. Non zero means periodic.
 *  @param h hgrid of the system, vector of 3 values.
 *  @param nseg_c number of segment of coarse data.
 *  @param nvctr_c number of point of coarse data.
 *  @param keyg_c array of size 2 * nseg_c * sizeof(uint), representing the beginning and end of coarse segments in the compressed data.
 *  @param keyv_c array of size nseg_c * sizeof(uint), representing the beginning of coarse segments in the uncompressed data.
 *  @param nseg_f number of segment of fine data.
 *  @param nvctr_f number of point of fine data.
 *  @param keyg_f array of size 2 * nseg_f * sizeof(uint), representing the beginning and end of fine segments in the compressed data.
 *  @param keyv_f array of size nseg_f * sizeof(uint), representing the beginning of fine segments in the uncompressed data.
 *  @param psi_c array of size nvctr_c * sizeof(double), containing coarse input and output data.
 *  @param psi_f array of size nvctr_f * 7 * sizeof(double), containing fine input and output data.
 *  @param pot array of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double) containing potential input data.
 *  @param work1 temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param work2 temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param work3 temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param work4 temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param epot potential energy of the system.
 *  @param ekinpot kinetic energy of the system.
*/
void FC_FUNC_(ocl_fulllocham_generic,OCL_FULLLOCHAM_GENERIC)(bigdft_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_uint *periodic,
                                          cl_double *h,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *pot,
                                          cl_mem *work1, cl_mem *work2,
                                          cl_mem *work3, cl_mem *work4,
                                          cl_double *epot, cl_double *ekinpot);

/** Performs a full local hamiltonian on a compressed wave function.
 *  @param command_queue used to process the data.
 *  @param dimensions of the input data, vector of 3 values.
 *  @param periodic periodicity of the convolution. Vector of three value, one for each dimension. Non zero means periodic.
 *  @param h hgrid of the system, vector of 3 values.
 *  @param nseg_c number of segment of coarse data.
 *  @param nvctr_c number of point of coarse data.
 *  @param keyg_c array of size 2 * nseg_c * sizeof(uint), representing the beginning and end of coarse segments in the compressed data.
 *  @param keyv_c array of size nseg_c * sizeof(uint), representing the beginning of coarse segments in the uncompressed data.
 *  @param nseg_f number of segment of fine data.
 *  @param nvctr_f number of point of fine data.
 *  @param keyg_f array of size 2 * nseg_f * sizeof(uint), representing the beginning and end of fine segments in the compressed data.
 *  @param keyv_f array of size nseg_f * sizeof(uint), representing the beginning of fine segments in the uncompressed data.
 *  @param psi_c_r array of size nvctr_c * sizeof(double), containing coarse input and output real data.
 *  @param psi_f_r array of size nvctr_f * 7 * sizeof(double), containing fine input and output real data.
 *  @param psi_c_i array of size nvctr_c * sizeof(double), containing coarse input and output imaginary data.
 *  @param psi_f_i array of size nvctr_f * 7 * sizeof(double), containing fine input and output imaginary data.
 *  @param pot array of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double) containing potential input data.
 *  @param psi_r temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param out_r temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param work_r temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param psi_i temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param out_i temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param work_i temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param kinres_r temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param kinres_i temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param nspinor identifies the number of internal components used in the computation (1: real only, 2: both, 4: unsupported).
 *  @param epot potential energy of the system.
 *  @param ekinpot kinetic energy of the system.
*/
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
                                          cl_double *epot, cl_double *ekinpot);

/** Performs a full local hamiltonian in periodic boundary conditions.
 *  @param command_queue used to process the data.
 *  @param dimensions of the input data, vector of 3 values.
 *  @param periodic periodicity of the convolution. Vector of three value, one for each dimension. Non zero means periodic.
 *  @param h hgrid of the system, vector of 3 values.
 *  @param k k-points values in the reduced Brillouin zone.
 *  @param nseg_c number of segment of coarse data.
 *  @param nvctr_c number of point of coarse data.
 *  @param keyg_c array of size 2 * nseg_c * sizeof(uint), representing the beginning and end of coarse segments in the compressed data.
 *  @param keyv_c array of size nseg_c * sizeof(uint), representing the beginning of coarse segments in the uncompressed data.
 *  @param nseg_f number of segment of fine data.
 *  @param nvctr_f number of point of fine data.
 *  @param keyg_f array of size 2 * nseg_f * sizeof(uint), representing the beginning and end of fine segments in the compressed data.
 *  @param keyv_f array of size nseg_f * sizeof(uint), representing the beginning of fine segments in the uncompressed data.
 *  @param psi_c array of size nvctr_c * sizeof(double), containing coarse input data.
 *  @param psi_f array of size nvctr_f * 7 * sizeof(double), containing fine input data.
 *  @param pot array of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double) containing potential input data.
 *  @param work1 temporary buffer of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double).
 *  @param work2 temporary buffer of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double).
 *  @param work3 temporary buffer of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double).
 *  @param work4 temporary buffer of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double).
 *  @param epot potential energy of the system.
 *  @param ekinpot kinetic energy of the system.
*/
void FC_FUNC_(ocl_fulllocham,OCL_FULLLOCHAM)(bigdft_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_double *h,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *pot,
                                          cl_mem *psi, cl_mem *out,
                                          cl_mem *work, cl_mem *kinres,
                                          cl_double *epot, cl_double *ekinpot);

/** Computes the local partial density of a wave function in periodic boundary condition.
 *  @param command_queue used to process the data.
 *  @param dimensions of the input data, vector of 3 values.
 *  @param hfac
 *  @apram iaddjmp
 *  @param nseg_c number of segment of coarse data.
 *  @param nvctr_c number of point of coarse data.
 *  @param keyg_c array of size 2 * nseg_c * sizeof(uint), representing the beginning and end of coarse segments in the compressed data.
 *  @param keyv_c array of size nseg_c * sizeof(uint), representing the beginning of coarse segments in the uncompressed data.
 *  @param nseg_f number of segment of fine data.
 *  @param nvctr_f number of point of fine data.
 *  @param keyg_f array of size 2 * nseg_f * sizeof(uint), representing the beginning and end of fine segments in the compressed data.
 *  @param keyv_f array of size nseg_f * sizeof(uint), representing the beginning of fine segments in the uncompressed data.
 *  @param psi_c array of size nvctr_c * sizeof(double), containing coarse input data.
 *  @param psi_f array of size nvctr_f * 7 * sizeof(double), containing fine input data.
 *  @param psi temporary buffer of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double).
 *  @param out temporary buffer of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double).
 *  @param work temporary buffer of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double).
 *  @param pot array of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double) containing potential output data.
*/
void FC_FUNC_(ocl_locden,OCL_LOCDEN)(bigdft_command_queue *command_queue,
                                     cl_uint *dimensions,
                                     cl_double *hfac,
                                     cl_uint *iaddjmp,
                                     cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                     cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                     cl_mem *psi_c, cl_mem *psi_f,
                                     cl_mem *psi, cl_mem *out, cl_mem *work,
                                     cl_mem *pot);

/** Computes the local partial density of a wave function.
 *  @param command_queue used to process the data.
 *  @param dimensions of the input data, vector of 3 values.
 *  @param periodic periodicity of the convolution. Vector of three value, one for each dimension. Non zero means periodic.
 *  @param hfac
 *  @param nseg_c number of segment of coarse data.
 *  @param nvctr_c number of point of coarse data.
 *  @param keyg_c array of size 2 * nseg_c * sizeof(uint), representing the beginning and end of coarse segments in the compressed data.
 *  @param keyv_c array of size nseg_c * sizeof(uint), representing the beginning of coarse segments in the uncompressed data.
 *  @param nseg_f number of segment of fine data.
 *  @param nvctr_f number of point of fine data.
 *  @param keyg_f array of size 2 * nseg_f * sizeof(uint), representing the beginning and end of fine segments in the compressed data.
 *  @param keyv_f array of size nseg_f * sizeof(uint), representing the beginning of fine segments in the uncompressed data.
 *  @param psi_c array of size nvctr_c * sizeof(double), containing coarse input data.
 *  @param psi_f array of size nvctr_f * 7 * sizeof(double), containing fine input data.
 *  @param psi temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param out temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param work temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param pot array of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double) containing potential output data.
*/
void FC_FUNC_(ocl_locden_generic,OCL_LOCDEN_GENERIC)(bigdft_command_queue *command_queue,
                                     cl_uint *dimensions,
                                     cl_uint *periodic,
                                     cl_double *hfac,
                                     cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                     cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                     cl_mem *psi_c, cl_mem *psi_f,
                                     cl_mem *psi, cl_mem *out, cl_mem *work,
                                     cl_mem *pot);

/** Preconditions the data for further hamiltonian minimizing, in periodic boundary conditions.
 *  @param command_queue used to process the data.
 *  @param dimensions of the input data, vector of 3 values.
 *  @param h hgrid of the system, vector of 3 values.
 *  @param c scaling factor.
 *  @param ncong number of iterations of the preconditioner.
 *  @param nseg_c number of segment of coarse data.
 *  @param nvctr_c number of point of coarse data.
 *  @param keyg_c array of size 2 * nseg_c * sizeof(uint), representing the beginning and end of coarse segments in the compressed data.
 *  @param keyv_c array of size nseg_c * sizeof(uint), representing the beginning of coarse segments in the uncompressed data.
 *  @param nseg_f number of segment of fine data.
 *  @param nvctr_f number of point of fine data.
 *  @param keyg_f array of size 2 * nseg_f * sizeof(uint), representing the beginning and end of fine segments in the compressed data.
 *  @param keyv_f array of size nseg_f * sizeof(uint), representing the beginning of fine segments in the uncompressed data.
 *  @param psi_c array of size nvctr_c * sizeof(double), containing coarse input and output data.
 *  @param psi_f array of size nvctr_f * 7 * sizeof(double), containing fine input and output data.
 *  @param psi_c_r temporary buffer of size nvctr_c * sizeof(double).
 *  @param psi_f_r temporary buffer of size nvctr_f * 7 * sizeof(double).
 *  @param psi_c_b temporary buffer of size nvctr_c * sizeof(double).
 *  @param psi_f_b temporary buffer of size nvctr_f * 7 * sizeof(double).
 *  @param psi_c_d temporary buffer of size nvctr_c * sizeof(double).
 *  @param psi_f_d temporary buffer of size nvctr_f * 7 * sizeof(double).
 *  @param work1 temporary buffer of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double).
 *  @param work2 temporary buffer of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double).
 *  @param work3 temporary buffer of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double).
 *  @param work4 temporary buffer of size (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double).
 */
void FC_FUNC_(ocl_preconditioner,OCL_PRECONDITIONER)(bigdft_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_double *h,
                                          cl_double *c,
                                          cl_uint *ncong,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *psi_c_r, cl_mem *psi_f_r,
                                          cl_mem *psi_c_b, cl_mem *psi_f_b,
                                          cl_mem *psi_c_d, cl_mem *psi_f_d,
                                          cl_mem *work1, cl_mem *work2, cl_mem *work3, cl_mem *work4);

/** Preconditions the data for further hamiltonian minimizing.
 *  @param command_queue used to process the data.
 *  @param dimensions of the input data, vector of 3 values.
 *  @param periodic periodicity of the convolution. Vector of three value, one for each dimension. Non zero means periodic.
 *  @param h hgrid of the system, vector of 3 values.
 *  @param c scaling factor.
 *  @param ncong number of iterations of the preconditioner.
 *  @param nseg_c number of segment of coarse data.
 *  @param nvctr_c number of point of coarse data.
 *  @param keyg_c array of size 2 * nseg_c * sizeof(uint), representing the beginning and end of coarse segments in the compressed data.
 *  @param keyv_c array of size nseg_c * sizeof(uint), representing the beginning of coarse segments in the uncompressed data.
 *  @param nseg_f number of segment of fine data.
 *  @param nvctr_f number of point of fine data.
 *  @param keyg_f array of size 2 * nseg_f * sizeof(uint), representing the beginning and end of fine segments in the compressed data.
 *  @param keyv_f array of size nseg_f * sizeof(uint), representing the beginning of fine segments in the uncompressed data.
 *  @param psi_c array of size nvctr_c * sizeof(double), containing coarse input and output data.
 *  @param psi_f array of size nvctr_f * 7 * sizeof(double), containing fine input and output data.
 *  @param psi_c_r temporary buffer of size nvctr_c * sizeof(double).
 *  @param psi_f_r temporary buffer of size nvctr_f * 7 * sizeof(double).
 *  @param psi_c_b temporary buffer of size nvctr_c * sizeof(double).
 *  @param psi_f_b temporary buffer of size nvctr_f * 7 * sizeof(double).
 *  @param psi_c_d temporary buffer of size nvctr_c * sizeof(double).
 *  @param psi_f_d temporary buffer of size nvctr_f * 7 * sizeof(double).
 *  @param work1 temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param work2 temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param work3 temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param work4 temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 */
void FC_FUNC_(ocl_preconditioner_generic,OCL_PRECONDITIONER_GENERIC)(bigdft_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_uint *periodic,
                                          cl_double *h,
                                          cl_double *c,
                                          cl_uint *ncong,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *psi_c_r, cl_mem *psi_f_r,
                                          cl_mem *psi_c_b, cl_mem *psi_f_b,
                                          cl_mem *psi_c_d, cl_mem *psi_f_d,
                                          cl_mem *work1, cl_mem *work2, cl_mem *work3, cl_mem *work4);

/** Preconditions the data for further hamiltonian minimizing.
 *  @param command_queue used to process the data.
 *  @param dimensions of the input data, vector of 3 values.
 *  @param periodic periodicity of the convolution. Vector of three value, one for each dimension. Non zero means periodic.
 *  @param h hgrid of the system, vector of 3 values.
 *  @param k k-points values in the reduced Brillouin zone.
 *  @param c scaling factor.
 *  @param ncong number of iterations of the preconditioner.
 *  @param nseg_c number of segment of coarse data.
 *  @param nvctr_c number of point of coarse data.
 *  @param keyg_c array of size 2 * nseg_c * sizeof(uint), representing the beginning and end of coarse segments in the compressed data.
 *  @param keyv_c array of size nseg_c * sizeof(uint), representing the beginning of coarse segments in the uncompressed data.
 *  @param nseg_f number of segment of fine data.
 *  @param nvctr_f number of point of fine data.
 *  @param keyg_f array of size 2 * nseg_f * sizeof(uint), representing the beginning and end of fine segments in the compressed data.
 *  @param keyv_f array of size nseg_f * sizeof(uint), representing the beginning of fine segments in the uncompressed data.
 *  @param psi_c_r array of size nvctr_c * sizeof(double), containing coarse input and output data.
 *  @param psi_f_r array of size nvctr_f * 7 * sizeof(double), containing fine input and output data.
 *  @param psi_c_i array of size nvctr_c * sizeof(double), containing coarse input and output data.
 *  @param psi_f_i array of size nvctr_f * 7 * sizeof(double), containing fine input and output data.
 *  @param psi_c_r_r temporary buffer of size nvctr_c * sizeof(double).
 *  @param psi_f_r_r temporary buffer of size nvctr_f * 7 * sizeof(double).
 *  @param psi_c_r_i temporary buffer of size nvctr_c * sizeof(double).
 *  @param psi_f_r_i temporary buffer of size nvctr_f * 7 * sizeof(double).
 *  @param psi_c_b_r temporary buffer of size nvctr_c * sizeof(double).
 *  @param psi_f_b_r temporary buffer of size nvctr_f * 7 * sizeof(double).
 *  @param psi_c_b_i temporary buffer of size nvctr_c * sizeof(double).
 *  @param psi_f_b_i temporary buffer of size nvctr_f * 7 * sizeof(double).
 *  @param psi_c_d_r temporary buffer of size nvctr_c * sizeof(double).
 *  @param psi_f_d_r temporary buffer of size nvctr_f * 7 * sizeof(double).
 *  @param psi_c_d_i temporary buffer of size nvctr_c * sizeof(double).
 *  @param psi_f_d_i temporary buffer of size nvctr_f * 7 * sizeof(double).
 *  @param work1_r temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param work2_r temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param work3_r temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param work4_r temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param work1_i temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param work2_i temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param work3_i temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param work4_i temporary buffer of size (2 * dimensions[0] + (periodic[0]?0:14+15)) * (2 * dimensions[1] + (periodic[1]?0:14+15)) * (2 * dimensions[2] + (periodic[2]?0:14+15)) * sizeof(double).
 *  @param nspinor identifies the number of internal components used in the computation (1: real only, 2: both, 4: unsupported).
 */
void FC_FUNC_(ocl_preconditioner_generic_k,OCL_PRECONDITIONER_GENERIC_K)(bigdft_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_uint *periodic,
                                          cl_double *h,
                                          cl_double *k,
                                          cl_double *c,
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
                                          cl_uint *nspinor, double *norm_sq_cf);

void FC_FUNC_(ocl_isf_to_daub,OCL_ISF_TO_DAUB)(bigdft_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_uint *periodic,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *psi, cl_mem *out,
                                          cl_mem *work, cl_mem *kinres);
void FC_FUNC_(ocl_daub_to_isf,OCL_DAUB_TO_ISF)(bigdft_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_uint *periodic,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *psi, cl_mem *out,
                                          cl_mem *work, cl_mem *kinres);
void FC_FUNC_(magic_filter_t_3d_generic,MAGIC_FILTER_T_3D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *tmp_dot, cl_mem *psi, cl_mem *out);
void FC_FUNC_(magic_filter_3d_generic,MAGIC_FILTER_3D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *tmp_dot, cl_mem *psi, cl_mem *out);
void FC_FUNC_(fft1d_d,FFT1D_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
void FC_FUNC_(fft1d_r_d,FFT1D_R_D)(bigdft_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
void FC_FUNC_(fft3d_d,FFT3D_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,cl_mem *psi,cl_mem *out,cl_mem *tmp);
void FC_FUNC_(fft3d_k_r2c_d,FFT3D_K_R2C_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,cl_mem *psi,cl_mem *out,cl_mem *tmp,cl_mem *k);
void FC_FUNC_(fft3d_r_d,FFT3D_R_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,cl_mem *psi,cl_mem *out,cl_mem *tmp);
void FC_FUNC_(fft3d_r_c2r_d,FFT1D_R_C2R_D)(bigdft_command_queue *command_queue, cl_uint *dimensions,cl_mem *psi,cl_mem *out,cl_mem *tmp);
void FC_FUNC_(fft3d_k_r2c_d_generic,FFT1D_K_R2C_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *psi,cl_mem *out,cl_mem *tmp,cl_mem *k);
void FC_FUNC_(fft3d_r_c2r_d_generic,FFT1D_R_C2R_D_GENERIC)(bigdft_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *psi,cl_mem *out,cl_mem *tmp);
void FC_FUNC_(ocl_create_read_buffer_host,OCL_CREATE_READ_BUFFER_HOST)(bigdft_context *context, cl_uint *size, void *host_ptr, cl_mem *buff_ptr );
void FC_FUNC_(ocl_create_write_buffer_host,OCL_CREATE_WRITE_BUFFER_HOST)(bigdft_context *context, cl_uint *size, void *host_ptr, cl_mem *buff_ptr );
void FC_FUNC_(ocl_create_read_write_buffer_host,OCL_CREATE_READ_WRITE_BUFFER_HOST)(bigdft_context *context, cl_uint *size, void *host_ptr, cl_mem *buff_ptr );
void FC_FUNC_(ocl_map_read_buffer,OCL_MAP_READ_BUFFER)(bigdft_command_queue *command_queue, cl_mem *buff_ptr, cl_uint *offset, cl_uint *size);
void FC_FUNC_(ocl_map_write_buffer,OCL_MAP_WRITE_BUFFER)(bigdft_command_queue *command_queue, cl_mem *buff_ptr, cl_uint *offset, cl_uint *size);
void FC_FUNC_(ocl_map_read_write_buffer,OCL_MAP_READ_WRITE_BUFFER)(bigdft_command_queue *command_queue, cl_mem *buff_ptr, cl_uint *offset, cl_uint *size);
void FC_FUNC_(ocl_map_read_buffer_async,OCL_MAP_READ_BUFFER_ASYNC)(bigdft_command_queue *command_queue, cl_mem *buff_ptr, cl_uint *offset, cl_uint *size);
void FC_FUNC_(ocl_map_write_buffer_async,OCL_MAP_WRITE_BUFFER_ASYNC)(bigdft_command_queue *command_queue, cl_mem *buff_ptr, cl_uint *offset, cl_uint *size);
void FC_FUNC_(ocl_map_read_write_buffer_async,OCL_MAP_READ_WRITE_BUFFER_ASYNC)(bigdft_command_queue *command_queue, cl_mem *buff_ptr, cl_uint *offset, cl_uint *size);
void FC_FUNC_(ocl_unmap_mem_object,OCL_UNMAP_MEM_OBJECT)(bigdft_command_queue *command_queue, cl_mem *buff_ptr, void *host_ptr);
void FC_FUNC_(ocl_pin_read_write_buffer,OCL_PIN_READ_WRITE_BUFFER)(bigdft_context *context, bigdft_command_queue *command_queue, cl_uint *size, void *host_ptr, cl_mem *buff_ptr );
void FC_FUNC_(ocl_pin_write_buffer,OCL_PIN_WRITE_BUFFER)(bigdft_context *context, bigdft_command_queue *command_queue, cl_uint *size, void *host_ptr, cl_mem *buff_ptr );
void FC_FUNC_(ocl_pin_read_buffer,OCL_PIN_READ_BUFFER)(bigdft_context *context, bigdft_command_queue *command_queue, cl_uint *size, void *host_ptr, cl_mem *buff_ptr );
void FC_FUNC_(ocl_pin_read_write_buffer_async,OCL_PIN_READ_WRITE_BUFFER_ASYNC)(bigdft_context *context, bigdft_command_queue *command_queue, cl_uint *size, void *host_ptr, cl_mem *buff_ptr );
void FC_FUNC_(ocl_pin_write_buffer_async,OCL_PIN_WRITE_BUFFER_ASYNC)(bigdft_context *context, bigdft_command_queue *command_queue, cl_uint *size, void *host_ptr, cl_mem *buff_ptr );
void FC_FUNC_(ocl_pin_read_buffer_async,OCL_PIN_READ_BUFFER_ASYNC)(bigdft_context *context, bigdft_command_queue *command_queue, cl_uint *size, void *host_ptr, cl_mem *buff_ptr );

#endif
