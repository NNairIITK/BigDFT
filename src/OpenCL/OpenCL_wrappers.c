#include <CL/cl.h>
#include <stdio.h>
#include <stdlib.h>
#include <config.h>
#include <math.h>

#define DEBUG 0

char * magicfilter1_program="\
__kernel void magicfilter1dKernel_l(size_t n, size_t ndat, __global const float *psi, __global float *out){\n\
const float filter[]={8.4334247333529341094733325815816e-7,\
		-0.1290557201342060969516786758559028e-4,\
		0.8762984476210559564689161894116397e-4,\
                -0.30158038132690463167163703826169879e-3,\
                0.174723713672993903449447812749852942e-2,\
                -0.942047030201080385922711540948195075e-2,\
                0.2373821463724942397566389712597274535e-1,\
                0.612625895831207982195380597e-1,\
                0.9940415697834003993178616713,\
                -0.604895289196983516002834636e-1,\
                -0.2103025160930381434955489412839065067e-1,\
                0.1337263414854794752733423467013220997e-1,\
                -0.344128144493493857280881509686821861e-2,\
                0.49443227688689919192282259476750972e-3,\
                -0.5185986881173432922848639136911487e-4,\
                2.72734492911979659657715313017228e-6};\
float *filt = filter + 8;\n\
size_t i = get_global_id(0);\n\
size_t j = get_global_id(1);\n\
float tt = 0.0;\n\
for( int l=-8; l<=7; l++){\n\
    int k = (i+l)%n;\n\
    tt = tt + psi[j + k*ndat] * filt[l];\n\
}\n\
out[j*n+i]=tt;\n\
};\n\
";
char * magicfilter1_program_optim="\
#define BLOCK_SIZE_I 16\n\
#define BLOCK_SIZE_J 8\n\
__kernel void magicfilter1dKernel_l_optim(size_t n, size_t ndat, __global const float *psi, __global float *out, __local float tmp[BLOCK_SIZE_J][3*BLOCK_SIZE_I]){\n\
size_t i = get_group_id(0);\n\
size_t j = get_group_id(1);\n\
size_t i2 = get_local_id(0);\n\
size_t j2 = get_local_id(1);\n\
size_t is = get_local_size(0);\n\
size_t js = get_local_size(1);\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
size_t ib;\n\
//load first matrix of size BLOCK_SIZE_I*BLOCK_SIZE_J\n\
if(i==0) ib = get_num_groups(0)-1;\n\
else ib = i-1;\n\
tmp[j2][i2]=psi[jg+(ib*is+i2)*ndat];\n\
//load second matrix\n\
tmp[j2][is+i2]=psi[jg+ig*ndat];\n\
if(i==get_num_groups(0)-1) ib = 0;\n\
else ib = i+1;\n\
tmp[j2][is*2+i2]=psi[jg+(ib*is+i2)*ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
//float *filt = filter + 8;\n\
float tt = 0.0;\n\
size_t base_i = i2+is;\n\
//for( int l=-8,m=7; l<0; l++,m--){\n\
//    tt += tmp[j2][base_i+l] * filt[l];\n\
//    tt += tmp[j2][base_i+m] * filt[m];\n\
//}\n\
tt += tmp[j2][base_i-8] * 8.4334247333529341094733325815816e-7;\n\
tt += tmp[j2][base_i+7] * 2.72734492911979659657715313017228e-6;\n\
tt += tmp[j2][base_i-7] * -0.1290557201342060969516786758559028e-4;\n\
tt += tmp[j2][base_i+6] * -0.5185986881173432922848639136911487e-4;\n\
tt += tmp[j2][base_i-6] * 0.8762984476210559564689161894116397e-4;\n\
tt += tmp[j2][base_i-5] * -0.30158038132690463167163703826169879e-3;\n\
tt += tmp[j2][base_i+5] * 0.49443227688689919192282259476750972e-3;\n\
tt += tmp[j2][base_i-4] * 0.174723713672993903449447812749852942e-2;\n\
tt += tmp[j2][base_i+4] * -0.344128144493493857280881509686821861e-2;\n\
tt += tmp[j2][base_i-3] * -0.942047030201080385922711540948195075e-2;\n\
tt += tmp[j2][base_i+3] * 0.1337263414854794752733423467013220997e-1;\n\
tt += tmp[j2][base_i+2] * -0.2103025160930381434955489412839065067e-1;\n\
tt += tmp[j2][base_i-2] * 0.2373821463724942397566389712597274535e-1;\n\
tt += tmp[j2][base_i+1] * -0.604895289196983516002834636e-1;\n\
tt += tmp[j2][base_i-1] * 0.612625895831207982195380597e-1;\n\
tt += tmp[j2][base_i] * 0.9940415697834003993178616713;\n\
out[jg*n+ig]=tt;\n\
};\
";

char * magicfilter1_program_optim_new="\
#define FILTER_WIDTH 16\n\
//n is supposed to be greater or equal than get_local_size(0)\n\
__kernel void magicfilter1dKernel_l_optim_new(size_t n, size_t ndat, __global const float *psi, __global float *out, __local float tmp[]){\n\
size_t ig = get_global_id(0);\n\
size_t jg = get_global_id(1);\n\
size_t i2 = get_local_id(0);\n\
size_t j2 = get_local_id(1);\n\
size_t is = get_local_size(0);\n\
size_t ib;\n\
size_t it;\n\
size_t base_i;\n\
//if data are ill dimentioned last block recomputes part of the data\n\
if( jg >= ndat ) return;\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
__local float * tmp_o = tmp + j2*(is+FILTER_WIDTH);\n\
base_i = FILTER_WIDTH/2+i2;\n\
//If I'm on the outside, select a border element to load\n\
if(i2 < FILTER_WIDTH/2)\n\
  { it = i2;\n\
    if (ig < FILTER_WIDTH/2)\n\
      { ib =  n + i2 - FILTER_WIDTH/2; }\n\
    else { ib = ig - FILTER_WIDTH/2; }\n\
    tmp_o[it]=psi[jg+ib*ndat];\n\
  }\n\
if (i2 >= (is - FILTER_WIDTH/2) || (ig >= n - FILTER_WIDTH/2))\n\
  { it = i2 + FILTER_WIDTH;\n\
    if (ig >= n - FILTER_WIDTH/2)\n\
      { ib = ig - n + FILTER_WIDTH/2; }\n\
    else { ib = ig + FILTER_WIDTH/2; }\n\
    tmp_o[it]=psi[jg+ib*ndat];\n\
  }\n\
//check boundaries\
//Load the element I am to calculate\n\
tmp_o[base_i]=psi[jg+ig*ndat];\n\
barrier(CLK_LOCAL_MEM_FENCE);\n\
\
float tt = 0.0;\n\
tmp_o = tmp_o + base_i;\n\
tt += tmp_o[-8] * 8.4334247333529341094733325815816e-7;\n\
tt += tmp_o[+7] * 2.72734492911979659657715313017228e-6;\n\
tt += tmp_o[-7] * -0.1290557201342060969516786758559028e-4;\n\
tt += tmp_o[+6] * -0.5185986881173432922848639136911487e-4;\n\
tt += tmp_o[-6] * 0.8762984476210559564689161894116397e-4;\n\
tt += tmp_o[-5] * -0.30158038132690463167163703826169879e-3;\n\
tt += tmp_o[+5] * 0.49443227688689919192282259476750972e-3;\n\
tt += tmp_o[-4] * 0.174723713672993903449447812749852942e-2;\n\
tt += tmp_o[+4] * -0.344128144493493857280881509686821861e-2;\n\
tt += tmp_o[-3] * -0.942047030201080385922711540948195075e-2;\n\
tt += tmp_o[+3] * 0.1337263414854794752733423467013220997e-1;\n\
tt += tmp_o[+2] * -0.2103025160930381434955489412839065067e-1;\n\
tt += tmp_o[-2] * 0.2373821463724942397566389712597274535e-1;\n\
tt += tmp_o[+1] * -0.604895289196983516002834636e-1;\n\
tt += tmp_o[-1] * 0.612625895831207982195380597e-1;\n\
tt += tmp_o[+0] * 0.9940415697834003993178616713;\n\
out[(jg*n+ig)]=tt;\n\
};\
";

#define oclErrorCheck(errorCode,message) if(errorCode!=CL_SUCCESS) { fprintf(stderr,"Error(%i) (%s: %s): %s\n", errorCode,__FILE__,__func__,message);exit(1);} 


cl_kernel magicfilter1d_kernel_l;
cl_kernel magicfilter1d_kernel_l_optim;
cl_kernel magicfilter1d_kernel_l_optim_new;

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
    cl_program magicfilter1dProgram = clCreateProgramWithSource(*context,1,(const char**) &magicfilter1_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(magicfilter1dProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(magicfilter1dProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    magicfilter1d_kernel_l=clCreateKernel(magicfilter1dProgram,"magicfilter1dKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    
    cl_program magicfilter1dProgram_optim = clCreateProgramWithSource(*context,1,(const char**) &magicfilter1_program_optim, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(magicfilter1dProgram_optim, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(magicfilter1dProgram_optim, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    magicfilter1d_kernel_l_optim=clCreateKernel(magicfilter1dProgram_optim,"magicfilter1dKernel_l_optim",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create optim kernel!");


    ciErrNum = CL_SUCCESS;
    cl_program magicfilter1dProgram_optim_new = clCreateProgramWithSource(*context,1,(const char**) &magicfilter1_program_optim_new, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(magicfilter1dProgram_optim_new, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error: Failed to build program!\n");
        char cBuildLog[10240];
        clGetProgramBuildInfo(magicfilter1dProgram_optim_new, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    magicfilter1d_kernel_l_optim_new=clCreateKernel(magicfilter1dProgram_optim_new,"magicfilter1dKernel_l_optim_new",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create new optim kernel!");
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
#define FILTER_WIDTH 16

void FC_FUNC_(magicfilter1d_l,MAGICFILTER1D_L)(cl_command_queue *command_queue, size_t *n,size_t *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %d, dimension dat: %d, psi: %p, out: %p\n",*command_queue, *n, *ndat, *psi, *out);
#endif
    cl_uint i = 0;
    clSetKernelArg(magicfilter1d_kernel_l, i++,sizeof(*n), (void*)n);
    clSetKernelArg(magicfilter1d_kernel_l, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(magicfilter1d_kernel_l, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(magicfilter1d_kernel_l, i++,sizeof(*out), (void*)out);
    size_t localWorkSize[] = { BLOCK_SIZE_I,BLOCK_SIZE_J };
    size_t globalWorkSize[] ={ shrRoundUp(BLOCK_SIZE_I,*n), shrRoundUp(BLOCK_SIZE_J,*ndat)};

    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_kernel_l, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue kernel!");
}

void FC_FUNC_(magicfilter1d_l_optim,MAGICFILTER1D_L_OPTIM)(cl_command_queue *command_queue, size_t *n,size_t *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %d, dimension dat: %d, psi: %p, out: %p\n",*command_queue, *n, *ndat, *psi, *out);
#endif
    cl_uint i = 0;
    clSetKernelArg(magicfilter1d_kernel_l_optim, i++,sizeof(*n), (void*)n);
    clSetKernelArg(magicfilter1d_kernel_l_optim, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(magicfilter1d_kernel_l_optim, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(magicfilter1d_kernel_l_optim, i++,sizeof(*out), (void*)out);
    clSetKernelArg(magicfilter1d_kernel_l_optim, i++,sizeof(float)*BLOCK_SIZE_J*BLOCK_SIZE_I*3, 0);
    size_t localWorkSize[] = { BLOCK_SIZE_I,BLOCK_SIZE_J };
    size_t globalWorkSize[] ={ shrRoundUp(BLOCK_SIZE_I,*n), shrRoundUp(BLOCK_SIZE_J,*ndat)};

    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_kernel_l_optim, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue optim kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %d, %d}\n",globalWorkSize[0],globalWorkSize[1]);
        exit(1);
    }   
}


void FC_FUNC_(magicfilter1d_l_optim_new,MAGICFILTER1D_L_OPTIM_NEW)(cl_command_queue *command_queue, size_t *n,size_t *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %d, dimension dat: %d, psi: %p, out: %p\n",*command_queue, *n, *ndat, *psi, *out);
#endif
    if(*n<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=256/FILTER_WIDTH;
    while (*n/block_size_i >= 1 && block_size_j > 1)
	{ block_size_i *= 2; block_size_j /= 2;}
    //size_t nbblock_i = ceil((float)(*n)/(float)block_size_i);
    //block_size_i = ( ceil((float)*n/(float)nbblock_i) );
    block_size_j = 256/block_size_i;

    block_size_i = 64;
    block_size_j = 4;
	

    cl_uint i = 0;
    clSetKernelArg(magicfilter1d_kernel_l_optim_new, i++,sizeof(*n), (void*)n);
    clSetKernelArg(magicfilter1d_kernel_l_optim_new, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(magicfilter1d_kernel_l_optim_new, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(magicfilter1d_kernel_l_optim_new, i++,sizeof(*out), (void*)out);
    clSetKernelArg(magicfilter1d_kernel_l_optim_new, i++,sizeof(float)*block_size_j*(block_size_i+FILTER_WIDTH), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
//        fprintf(stderr,"globalWorkSize = { %d, %d}\n",globalWorkSize[0],globalWorkSize[1]);
//        fprintf(stderr,"localWorkSize = { %d, %d}\n",localWorkSize[0],localWorkSize[1]);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_kernel_l_optim_new, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue optim kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %d, %d}\n",globalWorkSize[0],globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %d, %d}\n",localWorkSize[0],localWorkSize[1]);
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
