#include <CL/cl.h>
#include <stdio.h>
#include <stdlib.h>
#include <config.h>
#include <math.h>
#include "MagicFilter.h"
#include "Wavelet.h"
#include "Kinetic.h"
#include "Uncompress.h"

#define DEBUG 0
char * c_initialize_program="\
__kernel void c_initializeKernel_l(size_t n, size_t ndat, __global const float * x_in, __global float * y_in, float c) {\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
if( jg >= ndat ) return;\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
size_t pos = jg+ig*ndat;\n\
y_in[pos] = x_in[pos] * c;\n\
};\n\
__kernel void v_initializeKernel_l(size_t n, __global float * y_in, float v) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
y_in[ig] = v;\n\
};\n\
";

#define oclErrorCheck(errorCode,message) if(errorCode!=CL_SUCCESS) { fprintf(stderr,"Error(%i) (%s: %s): %s\n", errorCode,__FILE__,__func__,message);exit(1);} 


cl_kernel magicfilter1d_kernel_l;
cl_kernel ana1d_kernel_l;
cl_kernel syn1d_kernel_l;
cl_kernel kinetic1d_kernel_l;
cl_kernel c_initialize_kernel_l;
cl_kernel v_initialize_kernel_l;
cl_kernel uncompress_coarse_kernel_l;
cl_kernel uncompress_fine_kernel_l;
cl_kernel compress_coarse_kernel_l;
cl_kernel compress_fine_kernel_l;

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

void FC_FUNC_(ocl_build_kernels,OCL_BUILD_KERNELS)(cl_context * context) {
    cl_int ciErrNum = CL_SUCCESS;
    cl_program magicfilter1dProgram = clCreateProgramWithSource(*context,1,(const char**) &magicfilter1d_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(magicfilter1dProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build magicfilter1d program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(magicfilter1dProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    magicfilter1d_kernel_l=clCreateKernel(magicfilter1dProgram,"magicfilter1dKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");

    cl_program ana1dProgram = clCreateProgramWithSource(*context,1,(const char**) &ana1d_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(ana1dProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build ana1d program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(ana1dProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    ana1d_kernel_l=clCreateKernel(ana1dProgram,"ana1dKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");

    cl_program syn1dProgram = clCreateProgramWithSource(*context,1,(const char**) &syn1d_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(syn1dProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build syn1d program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(syn1dProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    syn1d_kernel_l=clCreateKernel(syn1dProgram,"syn1dKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");

    cl_program c_initializeProgram = clCreateProgramWithSource(*context,1,(const char**) &c_initialize_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(c_initializeProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build c_initialize program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(c_initializeProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    c_initialize_kernel_l=clCreateKernel(c_initializeProgram,"c_initializeKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    v_initialize_kernel_l=clCreateKernel(c_initializeProgram,"v_initializeKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");

    cl_program kinetic1dProgram = clCreateProgramWithSource(*context,1,(const char**) &kinetic1d_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(kinetic1dProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build kinetic1d program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(kinetic1dProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    kinetic1d_kernel_l=clCreateKernel(kinetic1dProgram,"kinetic1dKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");

    cl_program uncompressProgram = clCreateProgramWithSource(*context,1,(const char**) &uncompress_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(uncompressProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build uncompress program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(uncompressProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    uncompress_coarse_kernel_l=clCreateKernel(uncompressProgram,"uncompress_coarseKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    uncompress_fine_kernel_l=clCreateKernel(uncompressProgram,"uncompress_fineKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");

    cl_program compressProgram = clCreateProgramWithSource(*context,1,(const char**) &compress_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(compressProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build compress program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(compressProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    compress_coarse_kernel_l=clCreateKernel(compressProgram,"compress_coarseKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    compress_fine_kernel_l=clCreateKernel(compressProgram,"compress_fineKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");


}



void FC_FUNC_(ocl_create_gpu_context,OCL_CREATE_GPU_CONTEXT)(cl_context * context) {
    cl_int ciErrNum = CL_SUCCESS;
    *context = clCreateContextFromType(0, CL_DEVICE_TYPE_GPU, NULL, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p\n",*context);
#endif
    oclErrorCheck(ciErrNum,"Failed to create GPU context!");
}

void FC_FUNC_(ocl_create_cpu_context,OCL_CREATE_CPU_CONTEXT)(cl_context * context) {
    cl_int ciErrNum = CL_SUCCESS;
    *context = clCreateContextFromType(0, CL_DEVICE_TYPE_CPU, NULL, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p\n",*context);
#endif
    oclErrorCheck(ciErrNum,"Failed to create CPU context!");
}

void FC_FUNC_(ocl_create_read_buffer,OCL_CREATE_READ_BUFFER)(cl_context *context, size_t *size, cl_mem *buff_ptr) {
    cl_int ciErrNum = CL_SUCCESS;
    *buff_ptr = clCreateBuffer( *context, CL_MEM_READ_ONLY, *size, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, memory address: %p, size: %d\n",*context,*buff_ptr,*size);
#endif
    oclErrorCheck(ciErrNum,"Failed to create read buffer!");
}

void FC_FUNC_(ocl_create_read_write_buffer,OCL_CREATE_READ_WRITE_BUFFER)(cl_context *context, size_t *size, cl_mem *buff_ptr) {
    cl_int ciErrNum = CL_SUCCESS;
    *buff_ptr = clCreateBuffer( *context, CL_MEM_READ_WRITE, *size, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, memory address: %p, size: %d\n",*context,*buff_ptr,*size);
#endif
    oclErrorCheck(ciErrNum,"Failed to create read_write buffer!");
}

void FC_FUNC_(ocl_create_read_buffer_and_copy,OCL_CREATE_READ_BUFFER_AND_COPY)(cl_context *context, size_t *size, void *host_ptr, cl_mem *buff_ptr) {
    cl_int ciErrNum = CL_SUCCESS;
    *buff_ptr = clCreateBuffer( *context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, *size, host_ptr, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create initialized read buffer!");
}

void FC_FUNC_(ocl_create_write_buffer,OCL_CREATE_WRITE_BUFFER)(cl_context *context, size_t *size, cl_mem *buff_ptr) {
    cl_int ciErrNum = CL_SUCCESS;
    *buff_ptr = clCreateBuffer( *context, CL_MEM_WRITE_ONLY, *size, NULL, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, memory address: %p, size: %d\n",*context,*buff_ptr,*size);
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

void FC_FUNC_(ocl_enqueue_read_buffer,OCL_ENQUEUE_READ_BUFFER)(cl_command_queue *command_queue, cl_mem *buffer, size_t *size, void *ptr){
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, memory address: %p, size: %d, target: %p\n",*command_queue,*buffer,*size, ptr);
#endif
    cl_int ciErrNum = clEnqueueReadBuffer( *command_queue, *buffer, CL_TRUE, 0, *size, ptr, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue read buffer!");
}

void FC_FUNC_(ocl_enqueue_write_buffer,OCL_ENQUEUE_WRITE_BUFFER)(cl_command_queue *command_queue, cl_mem *buffer, size_t *size,	const void *ptr){
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, memory address: %p, size: %d, source: %p\n",*command_queue,*buffer,*size, ptr);
#endif
    cl_int ciErrNum = clEnqueueWriteBuffer( *command_queue, *buffer, CL_TRUE, 0, *size, ptr, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue write buffer!");
}

void FC_FUNC_(ocl_create_command_queue,OCL_CREATE_COMMAND_QUEUE)(cl_command_queue *hCmdQueue, cl_context *context){
    size_t nContextDescriptorSize;
    cl_int ciErrNum;
    clGetContextInfo(*context, CL_CONTEXT_DEVICES, 0, 0, &nContextDescriptorSize);
    cl_device_id * aDevices = (cl_device_id *) malloc(nContextDescriptorSize);
    clGetContextInfo(*context, CL_CONTEXT_DEVICES, nContextDescriptorSize, aDevices, 0);
    // create a command queue for first device the context reported
    *hCmdQueue = clCreateCommandQueue(*context, aDevices[0], 0, &ciErrNum);
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("contexte address: %p, command queue: %p\n",*context, *hCmdQueue);
#endif
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to create command queue!\n");
        exit(1);
    }
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


void magicfilter1dKernelCheck(size_t n, size_t ndat, double *psi, double *out){
double filter[]={8.4334247333529341094733325815816e-7,
                -0.1290557201342060969516786758559028e-4,
                0.8762984476210559564689161894116397e-4,
                -0.30158038132690463167163703826169879e-3,
                0.174723713672993903449447812749852942e-2,
                -0.942047030201080385922711540948195075e-2,
                0.2373821463724942397566389712597274535e-1,
                0.612625895831207982195380597e-1,
                0.9940415697834003993178616713,
                -0.604895289196983516002834636e-1,
                -0.2103025160930381434955489412839065067e-1,
                0.1337263414854794752733423467013220997e-1,
                -0.344128144493493857280881509686821861e-2,
                0.49443227688689919192282259476750972e-3,
                -0.5185986881173432922848639136911487e-4,
                2.72734492911979659657715313017228e-6};
double *filt = filter + 8;
int i,j,l,k;
for(j = 0; j < ndat; j++) {
    for( i = 0; i < n; i++) {
        double tt = 0.0;
        for( l = -8; l <= 7; l++) {
            k = (i+l)%n;
/*            printf("%d %lf\n", k, psi[k*ndat + j]);*/
            tt = tt + psi[k*ndat + j] * filt[l];
        }
/*        printf("%lf\n",tt);*/
        out[i + j*n] = tt;
/*        return;*/
    }

}
}

void FC_FUNC_(magicfilter1d_check,MAGICFILTER1D_CHECK)(size_t *n,size_t *ndat,void *psi,void *out){

	magicfilter1dKernelCheck(*n, *ndat, psi, out);

}

#define BLOCK_SIZE_I 64
#define BLOCK_SIZE_J 4

void FC_FUNC_(magicfilter1d_l,MAGICFILTER1D_L)(cl_command_queue *command_queue, size_t *n,size_t *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %d, dimension dat: %d, psi: %p, out: %p\n",*command_queue, *n, *ndat, *psi, *out);
#endif
    int FILTER_WIDTH = 16;
    if(*n<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=256/FILTER_WIDTH;
    while (*n > block_size_i >= 1 && block_size_j > 4)
	{ block_size_i *= 2; block_size_j /= 2;}

    cl_uint i = 0;
    clSetKernelArg(magicfilter1d_kernel_l, i++,sizeof(*n), (void*)n);
    clSetKernelArg(magicfilter1d_kernel_l, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(magicfilter1d_kernel_l, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(magicfilter1d_kernel_l, i++,sizeof(*out), (void*)out);
    clSetKernelArg(magicfilter1d_kernel_l, i++,sizeof(float)*block_size_j*(block_size_i+FILTER_WIDTH), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_kernel_l, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %d, %d}\n",globalWorkSize[0],globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %d, %d}\n",localWorkSize[0],localWorkSize[1]);
        exit(1);
    }   
}

void FC_FUNC_(ana1d_l,ANA1D_L)(cl_command_queue *command_queue, size_t *n,size_t *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %d, dimension dat: %d, psi: %p, out: %p\n",*command_queue, *n, *ndat, *psi, *out);
#endif
    int FILTER_WIDTH = 16;
    if(*n<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=256/FILTER_WIDTH;
    while (*n > block_size_i >= 1 && block_size_j > 8)
        { block_size_i *= 2; block_size_j /= 2;}
    cl_uint i = 0;
    clSetKernelArg(ana1d_kernel_l, i++,sizeof(*n), (void*)n);
    clSetKernelArg(ana1d_kernel_l, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(ana1d_kernel_l, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(ana1d_kernel_l, i++,sizeof(*out), (void*)out);
    clSetKernelArg(ana1d_kernel_l, i++,sizeof(float)*block_size_j*(block_size_i*2+FILTER_WIDTH), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, ana1d_kernel_l, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue ana1d_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %d, %d}\n",globalWorkSize[0],globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %d, %d}\n",localWorkSize[0],localWorkSize[1]);
        exit(1);
    }

}

void FC_FUNC_(syn1d_l,SYN1D_L)(cl_command_queue *command_queue, size_t *n,size_t *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %d, dimension dat: %d, psi: %p, out: %p\n",*command_queue, *n, *ndat, *psi, *out);
#endif
    int FILTER_WIDTH = 16;
    if(*n<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=256/FILTER_WIDTH;
    while (*n > block_size_i >= 1 && block_size_j > 8)
        { block_size_i *= 2; block_size_j /= 2;}
    cl_uint i = 0;
    clSetKernelArg(syn1d_kernel_l, i++,sizeof(*n), (void*)n);
    clSetKernelArg(syn1d_kernel_l, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(syn1d_kernel_l, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(syn1d_kernel_l, i++,sizeof(*out), (void*)out);
    clSetKernelArg(syn1d_kernel_l, i++,sizeof(float)*block_size_j*(block_size_i+FILTER_WIDTH), NULL);
    clSetKernelArg(syn1d_kernel_l, i++,sizeof(float)*block_size_j*(block_size_i+FILTER_WIDTH), NULL);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, syn1d_kernel_l, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {         
        fprintf(stderr,"Error %d: Failed to enqueue syn1d_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %d, %d}\n",globalWorkSize[0],globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %d, %d}\n",localWorkSize[0],localWorkSize[1]);
        exit(1);
    }  
}

void FC_FUNC_(kinetic1d_l,KINETIC1D_L)(cl_command_queue *command_queue, size_t *n,size_t *ndat, float *h, float*c, cl_mem *x, cl_mem *y, cl_mem *workx, cl_mem *worky,float *ekin){
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %d, dimension dat: %d, h: %f, c: %f, x: %p, workx: %p, y: %p, worky: %p\n",*command_queue, *n, *ndat, *h, *c, *x, *workx, *y, *worky);
#endif
    int FILTER_WIDTH = 30;
    if(*n<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=32, block_size_j=256/32;
    while (*n > block_size_i >= 1 && block_size_j > 4)
        { block_size_i *= 2; block_size_j /= 2;}
    cl_uint i = 0;
    clSetKernelArg(c_initialize_kernel_l, i++,sizeof(*n), (void*)n);
    clSetKernelArg(c_initialize_kernel_l, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(c_initialize_kernel_l, i++,sizeof(*x), (void*)x);
    clSetKernelArg(c_initialize_kernel_l, i++,sizeof(*worky), (void*)worky);
    clSetKernelArg(c_initialize_kernel_l, i++,sizeof(*c), (void*)c);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, c_initialize_kernel_l, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue c_initialize_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %d, %d}\n",globalWorkSize[0],globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %d, %d}\n",localWorkSize[0],localWorkSize[1]);
        exit(1);
    }
    i = 0;
    float scale = 0.5 / ( *h * (*h) );
    clSetKernelArg(kinetic1d_kernel_l, i++,sizeof(*n), (void*)n);
    clSetKernelArg(kinetic1d_kernel_l, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(kinetic1d_kernel_l, i++,sizeof(scale), (void*)&scale);
    clSetKernelArg(kinetic1d_kernel_l, i++,sizeof(*x), (void*)x);
    clSetKernelArg(kinetic1d_kernel_l, i++,sizeof(*workx), (void*)workx);
    clSetKernelArg(kinetic1d_kernel_l, i++,sizeof(*worky), (void*)worky);
    clSetKernelArg(kinetic1d_kernel_l, i++,sizeof(*y), (void*)y);
    clSetKernelArg(kinetic1d_kernel_l, i++,sizeof(float)*block_size_j*(block_size_i+FILTER_WIDTH), NULL);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kinetic1d_kernel_l, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue kinetic1d_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %d, %d}\n",globalWorkSize[0],globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %d, %d}\n",localWorkSize[0],localWorkSize[1]);
        exit(1);
    }
}

void FC_FUNC_(uncompress_l,UNCOMPRESS_L)(cl_command_queue *command_queue, size_t *n1, size_t *n2, size_t *n3,
                                       size_t *nseg_c, size_t *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                       size_t *nseg_f, size_t *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
				       cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi_out) {
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n1: %d, dimension n2: %d, dimension n3: %d,\n\
nseg_c: %d, nvctr_c: %d, keyg_c: %p, keyv_c: %p,\n\
nseg_f: %d, nvctr_f: %d, keyg_f: %p, keyv_f: %p,\n\
psi_c: %p, psi_f: %p, psi_out: %p\n",
*command_queue, *n1, *n2, *n3,
*nseg_c, *nvctr_c, *keyg_c, *keyv_c,
*nseg_f, *nvctr_f, *keyg_f, *keyv_f,
*psi_c, *psi_f, *psi_out);
#endif
    size_t block_size_i=64;
    size_t full_size = *n1 * *n2 * *n3 * 8;
    float init = 0.0;
    while ( block_size_i > full_size ) block_size_i /= 2;
    cl_uint i = 0;
    clSetKernelArg(v_initialize_kernel_l, i++, sizeof(full_size), (void*)&full_size);
    clSetKernelArg(v_initialize_kernel_l, i++, sizeof(*psi_out), (void*)psi_out);
    clSetKernelArg(v_initialize_kernel_l, i++, sizeof(init), (void*)&init);
    size_t localWorkSize[] = { block_size_i };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i, full_size) };
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, v_initialize_kernel_l, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue v_initialize_l kernel!\n",ciErrNum);
        exit(1);
    } 

    block_size_i=64;
    while ( block_size_i > *nvctr_c ) block_size_i /= 2;
    i = 0;
    clSetKernelArg(uncompress_coarse_kernel_l, i++, sizeof(*n1), (void*)n1);
    clSetKernelArg(uncompress_coarse_kernel_l, i++, sizeof(*n2), (void*)n2);
    clSetKernelArg(uncompress_coarse_kernel_l, i++, sizeof(*n3), (void*)n3);
    clSetKernelArg(uncompress_coarse_kernel_l, i++, sizeof(*nseg_c), (void*)nseg_c);
    clSetKernelArg(uncompress_coarse_kernel_l, i++, sizeof(*nvctr_c), (void*)nvctr_c);
    clSetKernelArg(uncompress_coarse_kernel_l, i++, sizeof(*keyg_c), (void*)keyg_c);
    clSetKernelArg(uncompress_coarse_kernel_l, i++, sizeof(*keyv_c), (void*)keyv_c);
    clSetKernelArg(uncompress_coarse_kernel_l, i++, sizeof(*psi_c), (void*)psi_c);
    clSetKernelArg(uncompress_coarse_kernel_l, i++, sizeof(*psi_out), (void*)psi_out);
    localWorkSize[0]= block_size_i ;
    globalWorkSize[0] = shrRoundUp(block_size_i, *nvctr_c) ;
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, uncompress_coarse_kernel_l, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue uncompress_coarse_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %d }\n",globalWorkSize[0]);
        fprintf(stderr,"localWorkSize = { %d }\n",localWorkSize[0]);
        exit(1);
    }

    block_size_i=64;
    while ( block_size_i > *nvctr_f ) block_size_i /= 2;
    clSetKernelArg(uncompress_fine_kernel_l, i++, sizeof(*n1), (void*)n1);
    clSetKernelArg(uncompress_fine_kernel_l, i++, sizeof(*n2), (void*)n2);
    clSetKernelArg(uncompress_fine_kernel_l, i++, sizeof(*n3), (void*)n3);
    clSetKernelArg(uncompress_fine_kernel_l, i++, sizeof(*nseg_f), (void*)nseg_f);
    clSetKernelArg(uncompress_fine_kernel_l, i++, sizeof(*nvctr_f), (void*)nvctr_f);
    clSetKernelArg(uncompress_fine_kernel_l, i++, sizeof(*keyg_f), (void*)keyg_f);
    clSetKernelArg(uncompress_fine_kernel_l, i++, sizeof(*keyv_f), (void*)keyv_f);
    clSetKernelArg(uncompress_fine_kernel_l, i++, sizeof(*psi_f), (void*)psi_f);
    clSetKernelArg(uncompress_fine_kernel_l, i++, sizeof(*psi_out), (void*)psi_out);
    localWorkSize[0]= block_size_i ;
    globalWorkSize[0] = shrRoundUp(block_size_i, *nvctr_f) ;
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, uncompress_fine_kernel_l, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue uncompress_fine_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %d }\n",globalWorkSize[0]);
        fprintf(stderr,"localWorkSize = { %d }\n",localWorkSize[0]);
        exit(1);
    }

}

void FC_FUNC_(compress_l,COMPRESS_L)(cl_command_queue *command_queue, size_t *n1, size_t *n2, size_t *n3,
                                     size_t *nseg_c, size_t *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                     size_t *nseg_f, size_t *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                     cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi) {
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n1: %d, dimension n2: %d, dimension n3: %d,\n\
nseg_c: %d, nvctr_c: %d, keyg_c: %p, keyv_c: %p,\n\
nseg_f: %d, nvctr_f: %d, keyg_f: %p, keyv_f: %p,\n\
psi_c: %p, psi_f: %p, psi: %p\n",
*command_queue, *n1, *n2, *n3,
*nseg_c, *nvctr_c, *keyg_c, *keyv_c,
*nseg_f, *nvctr_f, *keyg_f, *keyv_f,
*psi_c, *psi_f, *psi);
#endif
    size_t block_size_i=64;
    size_t full_size = *n1 * *n2 * *n3 * 8;
    float init = 0.0;
    while ( block_size_i > *nvctr_c ) block_size_i /= 2;
    cl_uint i = 0;
    clSetKernelArg(compress_coarse_kernel_l, i++, sizeof(*n1), (void*)n1);
    clSetKernelArg(compress_coarse_kernel_l, i++, sizeof(*n2), (void*)n2);
    clSetKernelArg(compress_coarse_kernel_l, i++, sizeof(*n3), (void*)n3);
    clSetKernelArg(compress_coarse_kernel_l, i++, sizeof(*nseg_c), (void*)nseg_c);
    clSetKernelArg(compress_coarse_kernel_l, i++, sizeof(*nvctr_c), (void*)nvctr_c);
    clSetKernelArg(compress_coarse_kernel_l, i++, sizeof(*keyg_c), (void*)keyg_c);
    clSetKernelArg(compress_coarse_kernel_l, i++, sizeof(*keyv_c), (void*)keyv_c);
    clSetKernelArg(compress_coarse_kernel_l, i++, sizeof(*psi_c), (void*)psi_c);
    clSetKernelArg(compress_coarse_kernel_l, i++, sizeof(*psi), (void*)psi);
    size_t localWorkSize[]= { block_size_i };
    size_t globalWorkSize[]={ shrRoundUp(block_size_i, *nvctr_c) } ;
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, compress_coarse_kernel_l, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue compress_coarse_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %d }\n",globalWorkSize[0]);
        fprintf(stderr,"localWorkSize = { %d }\n",localWorkSize[0]);
        exit(1);
    }

    block_size_i=64;
    while ( block_size_i > *nvctr_f ) block_size_i /= 2;
    clSetKernelArg(compress_fine_kernel_l, i++, sizeof(*n1), (void*)n1);
    clSetKernelArg(compress_fine_kernel_l, i++, sizeof(*n2), (void*)n2);
    clSetKernelArg(compress_fine_kernel_l, i++, sizeof(*n3), (void*)n3);
    clSetKernelArg(compress_fine_kernel_l, i++, sizeof(*nseg_f), (void*)nseg_f);
    clSetKernelArg(compress_fine_kernel_l, i++, sizeof(*nvctr_f), (void*)nvctr_f);
    clSetKernelArg(compress_fine_kernel_l, i++, sizeof(*keyg_f), (void*)keyg_f);
    clSetKernelArg(compress_fine_kernel_l, i++, sizeof(*keyv_f), (void*)keyv_f);
    clSetKernelArg(compress_fine_kernel_l, i++, sizeof(*psi_f), (void*)psi_f);
    clSetKernelArg(compress_fine_kernel_l, i++, sizeof(*psi), (void*)psi);
    localWorkSize[0]= block_size_i ;
    globalWorkSize[0] = shrRoundUp(block_size_i, *nvctr_f) ;
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, compress_fine_kernel_l, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue compress_fine_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %d }\n",globalWorkSize[0]);
        fprintf(stderr,"localWorkSize = { %d }\n",localWorkSize[0]);
        exit(1);
    }

}

void FC_FUNC_(ocl_finish,OCL_FINISH)(cl_command_queue *command_queue){
    cl_int ciErrNum;
    ciErrNum = clFinish(*command_queue);
    oclErrorCheck(ciErrNum,"Failed to finish!");
}

void FC_FUNC_(ocl_enqueue_barrier,OCL_ENQUEUE_BARRIER)(cl_command_queue *command_queue){
    cl_int ciErrNum;
    ciErrNum = clEnqueueBarrier(*command_queue);
    oclErrorCheck(ciErrNum,"Failed to enqueue barrier!");
}
