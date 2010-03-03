#include "OpenCL_wrappers.h"
#include "MagicFilter.h"

cl_kernel magicfilter1d_kernel_d;
cl_kernel magicfilter1d_pot_kernel_d;
cl_kernel magicfilter1d_t_kernel_d;
cl_kernel magicfiltershrink1d_kernel_d;
cl_kernel magicfiltergrow1d_kernel_d;

void build_magicfilter_kernels(cl_context * context){
    cl_int ciErrNum = CL_SUCCESS;
    cl_program magicfilter1dProgram = clCreateProgramWithSource(*context,1,(const char**) &magicfilter1d_program, NULL, &ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create program!");
    ciErrNum = clBuildProgram(magicfilter1dProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to build magicfilter1d program!\n",ciErrNum);
        char cBuildLog[10240];
        clGetProgramBuildInfo(magicfilter1dProgram, oclGetFirstDev(*context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
	fprintf(stderr,"%s\n",cBuildLog);
        exit(1);
    }
    ciErrNum = CL_SUCCESS;
    magicfiltergrow1d_kernel_d=clCreateKernel(magicfilter1dProgram,"magicfiltergrow1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    magicfiltershrink1d_kernel_d=clCreateKernel(magicfilter1dProgram,"magicfiltershrink1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    magicfilter1d_kernel_d=clCreateKernel(magicfilter1dProgram,"magicfilter1dKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    magicfilter1d_pot_kernel_d=clCreateKernel(magicfilter1dProgram,"magicfilter1d_potKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    magicfilter1d_t_kernel_d=clCreateKernel(magicfilter1dProgram,"magicfilter1d_tKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = clReleaseProgram(magicfilter1dProgram);
    oclErrorCheck(ciErrNum,"Failed to release program!");
}


void magicfilter1dKernelCheck(cl_uint n, cl_uint ndat, double *psi, double *out){
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

void FC_FUNC_(magicfilter1d_check,MAGICFILTER1D_CHECK)(cl_uint *n,cl_uint *ndat,void *psi,void *out){

	magicfilter1dKernelCheck(*n, *ndat, psi, out);

}

void FC_FUNC_(magicfiltershrink1d_d,MAGICFILTERSHRINK1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    cl_event e;
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %lu, dimension dat: %lu, psi: %p, out: %p\n",*command_queue, (long unsigned)*n, (long unsigned)*ndat, *psi, *out);
#endif
    int FILTER_WIDTH = 16;
    if(*n<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;

    cl_uint i = 0;
    ciErrNum = clSetKernelArg(magicfiltershrink1d_kernel_d, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(magicfiltershrink1d_kernel_d, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(magicfiltershrink1d_kernel_d, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(magicfiltershrink1d_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfiltershrink1d_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfiltershrink1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &e);
#if PROFILING
    event ev;
    ev.e = e;
    ev.comment = __func__;
    addToEventList(ev);
#endif
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfiltershrink1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
}
void FC_FUNC_(magicfiltergrow1d_d,MAGICFILTERGROW1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    cl_event e;
    cl_uint n1 = *n + 15;
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %lu, dimension dat: %lu, psi: %p, out: %p\n",*command_queue, (long unsigned)*n, (long unsigned)*ndat, *psi, *out);
#endif
    int FILTER_WIDTH = 16;
    if(n1<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;

    cl_uint i = 0;
    ciErrNum = clSetKernelArg(magicfiltergrow1d_kernel_d, i++,sizeof(n1), (void*)&n1);
    ciErrNum = clSetKernelArg(magicfiltergrow1d_kernel_d, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(magicfiltergrow1d_kernel_d, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(magicfiltergrow1d_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfiltergrow1d_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,n1), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfiltergrow1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &e);
#if PROFILING
    event ev;
    ev.e = e;
    ev.comment = __func__;
    addToEventList(ev);
#endif
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfiltergrow1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
}
void FC_FUNC_(magicfilter1d_d,MAGICFILTER1D_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    cl_event e;
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %lu, dimension dat: %lu, psi: %p, out: %p\n",*command_queue, (long unsigned)*n, (long unsigned)*ndat, *psi, *out);
#endif
    int FILTER_WIDTH = 16;
    if(*n<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;

    cl_uint i = 0;
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &e);
#if PROFILING
    event ev;
    ev.e = e;
    ev.comment = __func__;
    addToEventList(ev);
#endif
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
}

void FC_FUNC_(magicfilter1d_t_d,MAGICFILTER1D_T_D)(cl_command_queue *command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    cl_event e;
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %lu, dimension dat: %lu, psi: %p, out: %p\n",*command_queue, (long unsigned)*n, (long unsigned)*ndat, *psi, *out);
#endif
    int FILTER_WIDTH = 16;
    if(*n<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;

    cl_uint i = 0;
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_t_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &e);
#if PROFILING
    event ev;
    ev.e = e;
    ev.comment = __func__;
    addToEventList(ev);
#endif
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_t_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
}
void FC_FUNC_(magicfilter_n_self_d,MAGICFILTER_N_SELF_D)(cl_command_queue *command_queue, cl_uint *n1,cl_uint *n2,cl_uint *n3,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    int FILTER_WIDTH = 16;
    if(*n1<FILTER_WIDTH || *n2<FILTER_WIDTH || *n3<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;

    cl_uint i = 0;
    cl_uint ndat = *n1 * *n2;
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*n3), (void*)n3);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n3), shrRoundUp(block_size_j,ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
    i = 0;
    ndat = *n1 * *n3;
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*n2), (void*)n2);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    localWorkSize[0] = block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,*n2); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
    i = 0;
    ndat = *n2 * *n3;
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*n1), (void*)n1);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    localWorkSize[0] = block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,*n1); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
}

void FC_FUNC_(magicfilter_n_d,MAGICFILTER_N_D)(cl_command_queue *command_queue, cl_uint *n1,cl_uint *n2,cl_uint *n3,cl_mem *tmp,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    int FILTER_WIDTH = 16;
    if(*n1<FILTER_WIDTH || *n2<FILTER_WIDTH || *n3<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;

    cl_uint i = 0;
    cl_uint ndat = *n1 * *n2;
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*n3), (void*)n3);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n3), shrRoundUp(block_size_j,ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
    i = 0;
    ndat = *n1 * *n3;
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*n2), (void*)n2);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*tmp), (void*)tmp);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    localWorkSize[0] = block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,*n2); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
    i = 0;
    ndat = *n2 * *n3;
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*n1), (void*)n1);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*tmp), (void*)tmp);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    localWorkSize[0] = block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,*n1); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
}

void FC_FUNC_(magicfilter_t_self_d,MAGICFILTER_T_SELF_D)(cl_command_queue *command_queue, cl_uint *n1,cl_uint *n2,cl_uint *n3,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    int FILTER_WIDTH = 16;
    if(*n1<FILTER_WIDTH || *n2<FILTER_WIDTH || *n3<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;

    cl_uint i = 0;
    cl_uint ndat = *n1 * *n2;
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*n3), (void*)n3);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n3), shrRoundUp(block_size_j,ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_t_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_t_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
    i = 0;
    ndat = *n1 * *n3;
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*n2), (void*)n2);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    localWorkSize[0] = block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,*n2); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_t_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_t_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
    i = 0;
    ndat = *n2 * *n3;
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*n1), (void*)n1);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    localWorkSize[0] = block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,*n1); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_t_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_t_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
}

void FC_FUNC_(magicfilter_t_d,MAGICFILTER_T_D)(cl_command_queue *command_queue, cl_uint *n1,cl_uint *n2,cl_uint *n3,cl_mem *tmp,cl_mem *psi,cl_mem *out){
    cl_int ciErrNum;
    int FILTER_WIDTH = 16;
    if(*n1<FILTER_WIDTH || *n2<FILTER_WIDTH || *n3<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;

    cl_uint i = 0;
    cl_uint ndat = *n1 * *n2;
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*n3), (void*)n3);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n3), shrRoundUp(block_size_j,ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_t_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_t_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
    i = 0;
    ndat = *n1 * *n3;
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*n2), (void*)n2);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*tmp), (void*)tmp);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    localWorkSize[0] = block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,*n2); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_t_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_t_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
    i = 0;
    ndat = *n2 * *n3;
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*n1), (void*)n1);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*tmp), (void*)tmp);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    localWorkSize[0] = block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,*n1); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_t_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_t_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
}

void FC_FUNC_(potential_application_d,POTENTIAL_APPLICATION_D)(cl_command_queue *command_queue, cl_uint *n1,cl_uint *n2,cl_uint *n3, cl_mem *tmp, cl_mem *psi, cl_mem *out, cl_mem *pot) {
    cl_int ciErrNum;
    int FILTER_WIDTH = 16;
    if(*n1*2<FILTER_WIDTH || *n2*2<FILTER_WIDTH || *n3*2<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=FILTER_WIDTH;

    cl_uint i = 0;
    cl_uint n = *n3 * 2;
    cl_uint ndat = *n1 * *n2 * 4;
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(n), (void*)&n);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*tmp), (void*)tmp);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,n), shrRoundUp(block_size_j,ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
    i = 0;
    n = *n2 * 2;
    ndat = *n1 * *n3 * 4;
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(n), (void*)&n);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*tmp), (void*)tmp);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    localWorkSize[0] = block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,n); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
    i = 0;
    n = *n1 * 2;
    ndat = *n2 * *n3 * 4;
    ciErrNum = clSetKernelArg(magicfilter1d_pot_kernel_d, i++,sizeof(n), (void*)&n);
    ciErrNum = clSetKernelArg(magicfilter1d_pot_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_pot_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_pot_kernel_d, i++,sizeof(*pot), (void*)pot);
    ciErrNum = clSetKernelArg(magicfilter1d_pot_kernel_d, i++,sizeof(*tmp), (void*)tmp);
    ciErrNum = clSetKernelArg(magicfilter1d_pot_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    localWorkSize[0] = block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,n); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_pot_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
    i = 0;
    n = *n3 * 2;
    ndat = *n1 * *n2 * 4;
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(n), (void*)&n);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*tmp), (void*)tmp);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    localWorkSize[0] = block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,n); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_t_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_t_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
    i = 0;
    n = *n2 * 2;
    ndat = *n1 * *n3 * 4;
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(n), (void*)&n);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*tmp), (void*)tmp);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    localWorkSize[0] = block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,n); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_t_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_t_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }   
    i = 0;
    n = *n1 * 2;
    ndat = *n2 * *n3 * 4;
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(n), (void*)&n);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(ndat), (void*)&ndat);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*tmp), (void*)tmp);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(*out), (void*)out);
    ciErrNum = clSetKernelArg(magicfilter1d_t_kernel_d, i++,sizeof(cl_double)*block_size_j*(block_size_i+FILTER_WIDTH+1), 0);
    localWorkSize[0] = block_size_i; localWorkSize[1] = block_size_j;
    globalWorkSize[0] = shrRoundUp(block_size_i,n); globalWorkSize[1] = shrRoundUp(block_size_j,ndat);
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_t_kernel_d, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_t_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu, %lu}\n",(long unsigned)globalWorkSize[0],(long unsigned)globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %lu, %lu}\n",(long unsigned)localWorkSize[0],(long unsigned)localWorkSize[1]);
        exit(1);
    }

}


void clean_magicfilter_kernels(){
  clReleaseKernel(magicfilter1d_kernel_d);
  clReleaseKernel(magicfilter1d_pot_kernel_d);
  clReleaseKernel(magicfilter1d_t_kernel_d);
  clReleaseKernel(magicfiltershrink1d_kernel_d);
  clReleaseKernel(magicfiltergrow1d_kernel_d);
}
