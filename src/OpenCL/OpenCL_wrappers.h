#ifndef OPENCL_WRAPPERS_H
#define OPENCL_WRAPPERS_H

#include <CL/cl.h>
#include <stdio.h>
#include <stdlib.h>
#include <config.h>
#include <math.h>
#include <assert.h>

#define DEBUG 0
#define PROFILING 0

#define oclErrorCheck(errorCode,message) if(errorCode!=CL_SUCCESS) { fprintf(stderr,"Error(%i) (%s: %s): %s\n", errorCode,__FILE__,__func__,message);exit(1);} 

extern cl_kernel magicfilter1d_kernel_d;
extern cl_kernel ana1d_kernel_d;
extern cl_kernel syn1d_kernel_d;
extern cl_kernel kinetic1d_kernel_d;
extern cl_kernel c_initialize_kernel_d;
extern cl_kernel v_initialize_kernel_d;
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

int addToEventList (event ev);
extern event * event_list;
extern size_t event_number;

void FC_FUNC(init_event_list,INIT_EVENT_LIST)();
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

void FC_FUNC_(anashrink1d_d,ANASHRINK1D_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out);
void FC_FUNC_(ana1d_d,ANA1D_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out);
void FC_FUNC_(ana_d,ANA_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,cl_mem *tmp,cl_mem *psi, cl_mem *out);
void FC_FUNC_(ana_d_generic,ANA_D_GENERIC)(cl_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *psi, cl_mem *out);
void FC_FUNC_(ana_self_d,ANA_SELF_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,cl_mem *psi,cl_mem *out);
void FC_FUNC_(ana_self_d_generic,ANA_SELF_D_GENERIC)(cl_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *psi, cl_mem *out);
void FC_FUNC_(syngrow1d_d,SYNGROW1D_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out);
void FC_FUNC_(syn1d_d,SYN1D_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *out);
void FC_FUNC_(syn_d,SYN_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3, cl_mem *tmp, cl_mem *psi, cl_mem *out);
void FC_FUNC_(syn_d_generic,SYN_D_GENERIC)(cl_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *psi, cl_mem *out);
void FC_FUNC_(syn_self_d,SYN_SELF_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3, cl_mem *psi, cl_mem *out);
void FC_FUNC_(syn_self_d_generic,SYN_SELF_D_GENERIC)(cl_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *psi, cl_mem *out);

void FC_FUNC_(magicfiltershrink1d_d,MAGICFILTERSHRINK1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
void FC_FUNC_(magicfiltergrow1d_d,MAGICFILTERGROW1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
void FC_FUNC_(magicfilter1d_d,MAGICFILTER1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
void FC_FUNC_(magicfilter1d_pot_d,MAGICFILTER1D_POT_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, cl_mem *psi, cl_mem *pot, cl_mem *out);
void FC_FUNC_(magicfilter1d_t_d,MAGICFILTER1D_T_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out);
void FC_FUNC_(magicfilter_n_self_d,MAGICFILTER_N_SELF_D)(cl_command_queue *command_queue, cl_uint *n1,cl_uint *n2,cl_uint *n3,cl_mem *psi,cl_mem *out);
void FC_FUNC_(magicfilter_n_d,MAGICFILTER_N_D)(cl_command_queue *command_queue, cl_uint *n1,cl_uint *n2,cl_uint *n3,cl_mem *tmp,cl_mem *psi,cl_mem *out);
void FC_FUNC_(magicfilter_t_self_d,MAGICFILTER_T_SELF_D)(cl_command_queue *command_queue, cl_uint *n1,cl_uint *n2,cl_uint *n3,cl_mem *psi,cl_mem *out);
void FC_FUNC_(magicfilter_t_d,MAGICFILTER_T_D)(cl_command_queue *command_queue, cl_uint *n1,cl_uint *n2,cl_uint *n3,cl_mem *tmp,cl_mem *psi,cl_mem *out);
void FC_FUNC_(potential_application_d,POTENTIAL_APPLICATION_D)(cl_command_queue *command_queue, cl_uint *n1,cl_uint *n2,cl_uint *n3, cl_mem *tmp, cl_mem *psi, cl_mem *out, cl_mem *pot);
void FC_FUNC_(potential_application_d_generic,POTENTIAL_APPLICATION_D_GENERIC)(cl_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, cl_mem *tmp, cl_mem *psi, cl_mem *out, cl_mem *pot);

void FC_FUNC_(kinetic_k_d,KINETIC_K_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3, double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y, double * c_in,  double *k1, double *k2, double *k3);
void FC_FUNC_(kinetic_stable_d,KINETIC_STABLE_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3, double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y, cl_mem *tmp_x, cl_mem *tmp_y);
void FC_FUNC_(kinetic_d,KINETIC_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3, double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y);
void FC_FUNC_(kinetic_d_generic,KINETIC_D_GENERIC)(cl_command_queue *command_queue, cl_uint *dimensions, cl_uint *periodic, double *h, cl_mem *x, cl_mem *y, cl_mem *work_x, cl_mem *work_y);
void FC_FUNC_(kinetic1d_d,KINETIC1D_D)(cl_command_queue *command_queue, cl_uint *n, cl_uint *ndat, double *h, double*c, cl_mem *x, cl_mem *y, cl_mem *workx, cl_mem *worky,double *ekin);


void FC_FUNC_(uncompress_d,UNCOMPRESS_D)(cl_command_queue *command_queue, cl_uint *dimensions,
                                       cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                       cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                       cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi_out);
void FC_FUNC_(compress_d,COMPRESS_D)(cl_command_queue *command_queue, cl_uint *dimensions,
                                     cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                     cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                     cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi);

void FC_FUNC_(ocl_fulllocham,OCL_FULLLOCHAM)(cl_command_queue *command_queue,
                                          cl_uint *n1, cl_uint *n2, cl_uint *n3,
                                          double *h,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *pot,
                                          cl_mem *psi, cl_mem *out,
                                          cl_mem *work, cl_mem *kinres,
                                          double *epot,double *ekinpot);

void FC_FUNC_(ocl_fulllocham_generic,OCL_FULLLOCHAM_GENERIC)(cl_command_queue *command_queue,
                                          cl_uint *dimensions,
                                          cl_uint *periodic,
                                          double *h,
                                          cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c,
                                          cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                          cl_mem *psi_c, cl_mem *psi_f,
                                          cl_mem *pot,
                                          cl_mem *psi, cl_mem *out,
                                          cl_mem *work, cl_mem *kinres,
                                          double *epot,double *ekinpot);

#endif
