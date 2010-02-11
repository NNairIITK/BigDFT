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
extern cl_kernel magicfilter1d_t_kernel_l;
extern cl_kernel magicfilter1d_pot_kernel_l;
extern cl_kernel magicfilter1d_den_kernel_l;
extern cl_kernel ana1d_kernel_l;
extern cl_kernel syn1d_kernel_l;
extern cl_kernel kinetic1d_kernel_l;
extern cl_kernel c_initialize_kernel_l;
extern cl_kernel v_initialize_kernel_l;
extern cl_kernel uncompress_coarse_kernel_l;
extern cl_kernel uncompress_fine_kernel_l;
extern cl_kernel compress_coarse_kernel_l;
extern cl_kernel compress_fine_kernel_l;

void build_magicfilter_kernels(cl_context * context);
void build_kinetic_kernels(cl_context * context);
void build_wavelet_kernels(cl_context * context);
void build_uncompress_kernels(cl_context * context);
void build_initialize_kernels(cl_context * context);

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

#endif
