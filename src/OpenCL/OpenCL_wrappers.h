#ifndef OPENCL_WRAPPERS_H
#define OPENCL_WRAPPERS_H
#include <CL/cl.h>
#include <stdio.h>
#include <stdlib.h>
#include <config.h>
#include <math.h>
#include <assert.h>
#include <Tool.h>
#include <time.h>

/** @file OpenCL_wrappers.h
 *  @brief Contains global declarations and fortran bindings for OpenCL convolutions.
 *  Warning : every fortran visible procedure is passing its argument per address.
 *  Warning : every floating point data is double precision.
 */

/** Activate debugging info. */
#define DEBUG 0
/** Activate profiling info. */
#define PROFILING 0

#define oclErrorCheck(errorCode,message) if(errorCode!=CL_SUCCESS) { fprintf(stderr,"Error(%i) (%s: %s): %s\n", errorCode,__FILE__,__func__,message);exit(1);} 

extern cl_kernel c_initialize_kernel_d;
extern cl_kernel v_initialize_kernel_d;
extern cl_kernel p_initialize_kernel_d;
extern cl_kernel kinetic1d_kernel_d;
extern cl_kernel kinetic1d_f_kernel_d;
extern cl_kernel kinetic_k1d_kernel_d;
extern cl_kernel magicfilter1d_kernel_d;
extern cl_kernel magicfilter1d_straight_kernel_d;
extern cl_kernel magicfilter1d_block_kernel_d;
extern cl_kernel magicfilter1d_den_kernel_d;
extern cl_kernel magicfilter1d_pot_kernel_d;
extern cl_kernel magicfilter1d_t_kernel_d;
extern cl_kernel magicfiltershrink1d_kernel_d;
extern cl_kernel magicfiltergrow1d_kernel_d;
extern cl_kernel magicfiltergrow1d_pot_kernel_d;
extern cl_kernel reduction_kernel_d;
extern cl_kernel reduction_dot_kernel_d;
extern cl_kernel axpy_kernel_d;
extern cl_kernel axpy_offset_kernel_d;
extern cl_kernel scal_kernel_d;
extern cl_kernel copy_kernel_d;
extern cl_kernel dot_kernel_d;
extern cl_kernel set_kernel_d;
extern cl_kernel void_kernel;
extern cl_kernel uncompress_coarse_kernel_d;
extern cl_kernel uncompress_fine_kernel_d;
extern cl_kernel uncompress_scale_coarse_kernel_d;
extern cl_kernel uncompress_scale_fine_kernel_d;
extern cl_kernel compress_coarse_kernel_d;
extern cl_kernel compress_fine_kernel_d;
extern cl_kernel compress_scale_coarse_kernel_d;
extern cl_kernel compress_scale_fine_kernel_d;
extern cl_kernel scale_psi_fine_kernel_d;
extern cl_kernel scale_psi_coarse_kernel_d;
extern cl_kernel ana1d_kernel_d;
extern cl_kernel ana1d_block_kernel_d;
extern cl_kernel anashrink1d_kernel_d;
extern cl_kernel syn1d_kernel_d;
extern cl_kernel syngrow1d_kernel_d;
extern cl_kernel gemm_kernel_d;
extern cl_kernel gemm_kernel_d_tb;
extern cl_kernel gemm_kernel_d_ta;
extern cl_kernel gemm_kernel_d_tatb;
extern cl_kernel benchmark_flops_kernel_d;
extern cl_kernel benchmark_mops_kernel_d;
extern cl_program benchmarkProgram;

/** Creates magicfilter kernels. to be called after building the magicfilter programs. */
void create_magicfilter_kernels();
/** Compiles magicfilter programs in the given context. */
void build_magicfilter_programs(cl_context * context);
/** Releases magicfilter kernels. */
void clean_magicfilter_kernels();
void create_benchmark_kernels();
void build_benchmark_programs(cl_context * context);
void clean_benchmark_kernels();
void create_kinetic_kernels();
void build_kinetic_programs(cl_context * context);
void clean_kinetic_kernels();
void create_wavelet_kernels();
void build_wavelet_programs(cl_context * context);
void clean_wavelet_kernels();
void create_uncompress_kernels();
void build_uncompress_programs(cl_context * context);
void clean_uncompress_kernels();
void create_initialize_kernels();
void build_initialize_programs(cl_context * context);
void clean_initialize_kernels();
void create_reduction_kernels();
void build_reduction_programs(cl_context * context);
void clean_reduction_kernels();

/** Returns the first device available in a given context. */
cl_device_id oclGetFirstDev(cl_context cxGPUContext);

/** Returns the next integer that is equal or greater than global_size and a multiple of group_size. */
size_t shrRoundUp(size_t group_size, size_t global_size);

/** Structure associating an OpenCL event with a comment, for profiling purpose. */
typedef struct {
	cl_event e;
	char *comment;
} event;

/** Adds an event to the global event list. */
int addToEventList (event ev);
/** The global event list. */
extern event * event_list;
/** The number of event in the event_list. */
extern size_t event_number;

/** Reads the processor time stamp counter. */
void FC_FUNC_(rdtsc,RDTSC)(cl_ulong * t);
/** Return the real-time clock time in nanosecond since the epoch. */
void FC_FUNC_(nanosec,NANOSEC)(cl_ulong * t);

/** Initializes the event list. For profiling purpose. */
void FC_FUNC_(init_event_list,INIT_EVENT_LIST)();
/** Prints the event list. */
void FC_FUNC_(print_event_list,PRINT_EVENT_LIST)();
/** Buids and create the OpenCL kernel int the given context. */
void FC_FUNC_(ocl_build_kernels,OCL_BUILD_KERNELS)(cl_context * context);
/** Creates a context containing all GPUs from the default platform */
void FC_FUNC_(ocl_create_gpu_context,OCL_CREATE_GPU_CONTEXT)(cl_context * context);
/** Creates a context containing all CPUs from the default platform */
void FC_FUNC_(ocl_create_cpu_context,OCL_CREATE_CPU_CONTEXT)(cl_context * context);
/** Creates a OpenCL read only buffer.
 *  @param context where the buffer is created.
 *  @param size of the buffer.
 *  @param buff_ptr return value : a buffer object reference.
 */
void FC_FUNC_(ocl_create_read_buffer,OCL_CREATE_READ_BUFFER)(cl_context *context, cl_uint *size, cl_mem *buff_ptr);
/** Creates an OpenCL buffer.
 *  @param context where the buffer is created.
 *  @param size of the buffer.
 *  @param buff_ptr return value : a buffer object reference.
 */
void FC_FUNC_(ocl_create_read_write_buffer,OCL_CREATE_READ_WRITE_BUFFER)(cl_context *context, cl_uint *size, cl_mem *buff_ptr);
void FC_FUNC_(ocl_create_read_buffer_and_copy,OCL_CREATE_READ_BUFFER_AND_COPY)(cl_context *context, cl_uint *size, void *host_ptr, cl_mem *buff_ptr);
/** Creates a OpenCL write only buffer.
 *  @param context where the buffer is created.
 *  @param size of the buffer.
 *  @param buff_ptr return value : a buffer object reference.
 */
void FC_FUNC_(ocl_create_write_buffer,OCL_CREATE_WRITE_BUFFER)(cl_context *context, cl_uint *size, cl_mem *buff_ptr);
/** Releases an OpenCL buffer. */
void FC_FUNC_(ocl_release_mem_object,OCL_RELEASE_MEM_OBJECT)(cl_mem *buff_ptr);
/** Copies data from an OpenCL buffer to Host memory.
 *  @param command_queue a pointer to the command queue used to make the copy.
 *  @param buffer to copy data from.
 *  @param size of the data to copy.
 *  @param host_ptr to copy the data to.
 */
void FC_FUNC_(ocl_enqueue_read_buffer,OCL_ENQUEUE_READ_BUFFER)(cl_command_queue *command_queue, cl_mem *buffer, cl_uint *size, void *host_ptr);
/** Copies data from Host memory to an OpenCL buffer.
 *  @param command_queue a pointer to the command queue used to make the copy.
 *  @param buffer to copy data to.
 *  @param size of the data to copy.
 *  @param host_ptr to copy the data from.
 */
void FC_FUNC_(ocl_enqueue_write_buffer,OCL_ENQUEUE_WRITE_BUFFER)(cl_command_queue *command_queue, cl_mem *buffer, cl_uint *size, const void *host_ptr);
/** Creates a command queue in the given context, associating it to the first device in the context */
void FC_FUNC_(ocl_create_command_queue,OCL_CREATE_COMMAND_QUEUE)(cl_command_queue *command_queue, cl_context *context);
/** Creates a command queue in the given context, associating it to the device specified by index modulo the number of device. */
void FC_FUNC_(ocl_create_command_queue_id,OCL_CREATE_COMMAND_QUEUE_ID)(cl_command_queue *command_queue, cl_context *context, cl_uint *index);
/** Waits for all commands in a queue to complete. */
void FC_FUNC_(ocl_finish,OCL_FINISH)(cl_command_queue *command_queue);
/** Enqueues a barrier in a queue. Commands enqueued after the barrier will wait
 *  for commands enqueued before the barrier to be processed before being sent to the device. */
void FC_FUNC_(ocl_enqueue_barrier,OCL_ENQUEUE_BARRIER)(cl_command_queue *command_queue);
/** Releases the command queue and the context, and beforehand releases the program and kernels. */
void FC_FUNC_(ocl_clean,OCL_CLEAN)(cl_command_queue *command_queue, cl_context *context);

/** Performs the one dimensional wavelet analysis and transposition with periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * (2 * n) * sizeof(double), stored in column major order.
 *  @param out output buffer of size (2 * n) * ndat * sizeof(double), stored in column major order.
 */
void FC_FUNC_(ana1d_d,ANA1D_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out);
/** Performs the one dimensional wavelet analysis and transposition with open boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * (2 * n + 14) * sizeof(double), stored in column major order.
 *  @param out output buffer of size (2 * n) * ndat * sizeof(double), stored in column major order.
 */
void FC_FUNC_(anashrink1d_d,ANASHRINK1D_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out);
/** Slightly more performing version of anashrink1d_d. @see anashrink1d_d. */
void FC_FUNC_(ana1d_block_d,ANA1D_BLOCK_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
/** Performs the three-dimensional wavelet analysis with periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param tmp temporary buffer to store intermediate results. Dimensions : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double).
 *  @param psi input buffer of dimension : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 *  @param out output buffer of dimensions : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 */
void FC_FUNC_(ana_d,ANA_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out);
/** Slightly more performing version of ana_d. @see ana_d. */
void FC_FUNC_(ana_block_d,ANA_BLOCK_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out);
/** Performs the three-dimensional wavelet analysis with periodic or non periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param periodic periodicity of the convolution. Vector of three value, one for each dimension. Non zero means periodic.
 *  @param tmp temporary buffer to store intermediate results. Must be of at least (2 * dimensions[0] + (periodic[0]?0:14)) * (2 * dimensions[1] + (periodic[1]?0:14)) * (2 * dimensions[2]) * sizeof(double) in size.
 *  @param psi input buffer of dimension : (2 * dimensions[0] + (periodic[0]?0:14)) * (2 * dimensions[1] + (periodic[1]?0:14)) * (2 * dimensions[2] + (periodic[2]?0:14)) * sizeof(double). Stored in collumn major order.
 *  @param out output buffer of dimensions : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 */
void FC_FUNC_(ana_d_generic,ANA_D_GENERIC)(cl_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *psi, cl_mem *out);
/** Version of ana_d without the temporary buffer, psi is erased during the computation. @see ana_d. */
void FC_FUNC_(ana_self_d,ANA_SELF_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *psi, cl_mem *out);
/** Version of ana_sef_d without the temporary buffer, psi is erased during the computation. @see ana_d_generic. */
void FC_FUNC_(ana_self_d_generic,ANA_SELF_D_GENERIC)(cl_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *psi, cl_mem *out);
/** Performs the one dimensional wavelet synthesis and transposition with periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * (2 * n) * sizeof(double), stored in column major order.
 *  @param out output buffer of size (2 * n) * ndat * sizeof(double), stored in column major order.
 */
void FC_FUNC_(syn1d_d,SYN1D_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out);
/** Performs the one dimensional wavelet analysis and transposition with open boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * (2 * n) * sizeof(double), stored in column major order.
 *  @param out output buffer of size (2 * n + 14) * ndat * sizeof(double), stored in column major order.
 */
void FC_FUNC_(syngrow1d_d,SYNGROW1D_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out);
/** Performs the three-dimensional wavelet synthesis with periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param tmp temporary buffer to store intermediate results. Dimensions : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double).
 *  @param psi input buffer of dimension : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 *  @param out output buffer of dimensions : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 */
void FC_FUNC_(syn_d,SYN_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out);
/** Performs the three-dimensional wavelet analysis with periodic or non periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param dimensions of the input data. Vector of three values, one for each dimension.
 *  @param periodic periodicity of the convolution. Vector of three value, one for each dimension. Non zero means periodic.
 *  @param tmp temporary buffer to store intermediate results. Must be of at least (2 * dimensions[0]) * (2 * dimensions[1] + (periodic[1]?0:14)) * (2 * dimensions[2] + (periodic[1]?0:14)) * sizeof(double) in size.
 *  @param psi input buffer of dimension : (2 * dimensions[0]) * (2 * dimensions[1]) * (2 * dimensions[2]) * sizeof(double). Stored in column major order.
 *  @param out output buffer of dimensions : (2 * dimensions[0] + (periodic[0]?0:14)) * (2 * dimensions[1] + (periodic[1]?0:14)) * (2 * dimensions[2] + (periodic[2]?0:14)) * sizeof(double). Stored in column major order.
 */
void FC_FUNC_(syn_d_generic,SYN_D_GENERIC)(cl_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *psi, cl_mem *out);
/** Version of syn_d without the temporary buffer, psi is erased during the computation. @see syn_d. */
void FC_FUNC_(syn_self_d,SYN_SELF_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *psi, cl_mem *out);
/** Version of syn_sef_d without the temporary buffer, psi is erased during the computation. @see syn_d_generic. */
void FC_FUNC_(syn_self_d_generic,SYN_SELF_D_GENERIC)(cl_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *psi, cl_mem *out);

/** Performs the one dimensional magicfilter and transposition with periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * n * sizeof(double), stored in collumn major order.
 *  @param out output buffer of size n * ndat * sizeof(double), stored in collumn major order.
 */
void FC_FUNC_(magicfilter1d_d,MAGICFILTER1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
/** Performs the one dimensional magicfilter and transposition with open boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * n * sizeof(double), stored in column major order.
 *  @param out output buffer of size (n + 15) * ndat * sizeof(double), stored in column major order.
 */
void FC_FUNC_(magicfiltergrow1d_d,MAGICFILTERGROW1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
/** Performs the one dimensional reciprocal magicfilter and transposition with periodic boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * n * sizeof(double), stored in collumn major order.
 *  @param out output buffer of size n * ndat * sizeof(double), stored in collumn major order.
 */
void FC_FUNC_(magicfilter1d_t_d,MAGICFILTER1D_T_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
/** Performs the one dimensional reciprocal magicfilter and transposition with open boundary conditions.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * (n +15) * sizeof(double), stored in column major order.
 *  @param out output buffer of size n * ndat * sizeof(double), stored in column major order.
 */
void FC_FUNC_(magicfiltershrink1d_d,MAGICFILTERSHRINK1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
/** Performs the one dimensional magicfilter and transposition with periodic boundary conditions. Storage of matrix is changed with respect to magicfilter1d_d.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size n * ndat * sizeof(double), stored in collumn major order.
 *  @param out output buffer of size ndat * n * sizeof(double), stored in collumn major order.
 */
void FC_FUNC_(magicfilter1d_straight_d,MAGICFILTER1D_STRAIGHT_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
/** Slightly more performing version of magicfilter1d_d. @see magicfilter1d_d. */
void FC_FUNC_(magicfilter1d_block_d,MAGICFILTER1D_BLOCK_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
/** Performs the one dimensional magicfilter and transposition with periodic boundary conditions and multiplies by a potential.
 *  @param command_queue used to process the convolution.
 *  @param n size of the dimension to process the convolution.
 *  @param ndat size of the other dimension.
 *  @param psi input buffer of size ndat * n * sizeof(double), stored in collumn major order.
 *  @param pot potential applied during the convolution. Size is n * ndat * sizeof(double), stored in collumn major order.
 *  @param out output buffer of size n * ndat * sizeof(double), stored in collumn major order.
 */
void FC_FUNC_(magicfilter1d_pot_d,MAGICFILTER1D_POT_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *pot, cl_mem *out);
void FC_FUNC_(magicfilter_n_self_d,MAGICFILTER_N_SELF_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *psi, cl_mem *out);
void FC_FUNC_(magicfilter_n_d,MAGICFILTER_N_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out);
void FC_FUNC_(magicfilter_n_straight_d,MAGICFILTER_N_STRAIGHT_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out);
void FC_FUNC_(magicfilter_n_block_d,MAGICFILTER_N_BLOCK_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out);
void FC_FUNC_(magicfilter_den_d,MAGICFILTER_DEN_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out);
void FC_FUNC_(magicfilter_t_self_d,MAGICFILTER_T_SELF_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *psi, cl_mem *out);
void FC_FUNC_(magicfilter_t_d,MAGICFILTER_T_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out);
void FC_FUNC_(potential_application_d,POTENTIAL_APPLICATION_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_mem *tmp, cl_mem *psi, cl_mem *out, cl_mem *pot);
void FC_FUNC_(potential_application_d_generic,POTENTIAL_APPLICATION_D_GENERIC)(cl_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *tmp_dot, cl_mem *psi, cl_mem *out, cl_mem *pot, cl_double *epot);

void FC_FUNC_(transpose_d,TRANSPOSE_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
void FC_FUNC_(benchmark_flops_d,BENCHMARK_FLOPS_D)(cl_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *out);
void FC_FUNC_(benchmark_mops_d,BENCHMARK_MOPS_D)(cl_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *out);

void FC_FUNC_(kinetic_k_d,KINETIC_K_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y, cl_double * c_in,  cl_double *k);
void FC_FUNC_(kinetic_stable_d,KINETIC_STABLE_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y, cl_mem *tmp_x, cl_mem *tmp_y);
void FC_FUNC_(kinetic_d,KINETIC_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y);
void FC_FUNC_(kinetic_d_generic,KINETIC_D_GENERIC)(cl_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y);
void FC_FUNC_(kinetic1d_d,KINETIC1D_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_double *h, cl_double*c, cl_mem *x, cl_mem *y, cl_mem *workx, cl_mem *worky, cl_double *ekin);

void FC_FUNC_(asum_self_d,ASUM_SELF_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work, cl_double *out);
void FC_FUNC_(nrm2sq_self_d,NRM2SQ_SELF_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work, cl_double *out);
void FC_FUNC_(asum_d,ASUM_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work1, cl_mem *work2, cl_double *out);
void FC_FUNC_(nrm2sq_d,NRM2SQ_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *work1, cl_mem *work2, cl_double *out);
void FC_FUNC_(axpy_self_d,AXPY_SELF_D)(cl_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *in, cl_mem *inout);
void FC_FUNC_(axpy_d,AXPY_D)(cl_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *x, cl_mem *y, cl_mem *out);
void FC_FUNC_(scal_d,SCAL_D)(cl_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *in, cl_mem *out);
void FC_FUNC_(scal_self_d,SCAL_SELF_D)(cl_command_queue *command_queue, cl_uint *n, cl_double *alpha, cl_mem *inout);
void FC_FUNC_(dot_d,DOT_D)(cl_command_queue *command_queue, cl_uint *ndat, cl_mem *x, cl_mem *y, cl_mem *work1, cl_mem *work2, cl_double *out);
void FC_FUNC_(set_d,SET_D)(cl_command_queue *command_queue, cl_uint *n, cl_double *val, cl_mem *x);
void FC_FUNC_(copy_d,COPY_D)(cl_command_queue *command_queue, cl_uint *n, cl_mem *in, cl_mem *out);
void FC_FUNC_(axpy_offset_self_d,AXPY_OFFSET_SELF_D)(cl_command_queue *command_queue, cl_uint *n, cl_double *alpha,
                                                                            cl_uint *offset_x, cl_mem *x,
                                                                            cl_uint *offset_y, cl_mem *y);
void FC_FUNC_(axpy_offset_d,AXPY_OFFSET_D)(cl_command_queue *command_queue, cl_uint *n, cl_double *alpha,
                                                                            cl_uint *offset_x, cl_mem *x,
                                                                            cl_uint *offset_y, cl_mem *y,
                                                                            cl_uint *offset_z, cl_mem *z);
void FC_FUNC_(gemm_d,GEMM_D)(cl_command_queue *command_queue, char *transa, char *transb, cl_uint *m, cl_uint *n, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc);
void FC_FUNC_(gemmsy_d,GEMMSY_D)(cl_command_queue *command_queue, char *transa, char *transb, cl_uint *m, cl_uint *n, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc);
void FC_FUNC_(gemm_block_d,GEMM_BLOCK_D)(cl_command_queue *command_queue, char *transa, char *transb, cl_uint *m, cl_uint *n, cl_uint *k, cl_double *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double *beta, cl_mem *c, cl_uint *ldc);
void FC_FUNC_(gemm_z,GEMM_Z)(cl_command_queue *command_queue, char *transa, char *transb, cl_uint *m, cl_uint *n, cl_uint *k, cl_double2 *alpha, cl_mem *a, cl_uint *lda, cl_mem *b, cl_uint *ldb, cl_double2 *beta, cl_mem *c, cl_uint *ldc);

void FC_FUNC_(uncompress_d,UNCOMPRESS_D)(cl_command_queue *command_queue, cl_uint *dimensions,
                                       cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                       cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                       cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi_out);
void FC_FUNC_(compress_d,COMPRESS_D)(cl_command_queue *command_queue, cl_uint *dimensions,
                                     cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                     cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                     cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi);
void FC_FUNC_(scale_psi_d,SCALE_PSI_D)(cl_command_queue *command_queue, cl_uint *nvctr_c, cl_uint *nvctr_f, cl_double *h, cl_double *c, cl_mem *psi_c,  cl_mem *psi_f);
void FC_FUNC_(uncompress_scale_d,UNCOMPRESS_SCALE_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_double *h, cl_double *c,
                                       cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                       cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                       cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi_out);
void FC_FUNC_(compress_scale_d,COMPRESS_SCALE_D)(cl_command_queue *command_queue, cl_uint *dimensions, cl_double *h, cl_double *c,
                                     cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                     cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                     cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi);

void FC_FUNC_(ocl_fulllocham,OCL_FULLLOCHAM)(cl_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_double *h,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *pot,
                                          cl_mem *psi, cl_mem *out,
                                          cl_mem *work, cl_mem *kinres,
                                          cl_double *epot, cl_double *ekinpot);

void FC_FUNC_(ocl_fulllocham_generic,OCL_FULLLOCHAM_GENERIC)(cl_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_uint *periodic,
                                          cl_double *h,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *pot,
                                          cl_mem *psi, cl_mem *out,
                                          cl_mem *work, cl_mem *kinres,
                                          cl_double *epot, cl_double *ekinpot);

void FC_FUNC_(ocl_preconditioner,OCL_PRECONDITIONER)(cl_command_queue *command_queue,
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
#endif
