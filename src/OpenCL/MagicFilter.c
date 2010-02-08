#include "OpenCL_wrappers.h"
#include "MagicFilter.h"

cl_kernel magicfilter1d_kernel_l;
cl_kernel magicfilter1d_t_kernel_l;
cl_kernel magicfilter1d_pot_kernel_l;
cl_kernel magicfilter1d_den_kernel_l;

void build_magicfilter_kernels(cl_context * context){
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
    ciErrNum = CL_SUCCESS;
    magicfilter1d_t_kernel_l=clCreateKernel(magicfilter1dProgram,"magicfilter1d_tKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    magicfilter1d_pot_kernel_l=clCreateKernel(magicfilter1dProgram,"magicfilter1d_potKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    magicfilter1d_den_kernel_l=clCreateKernel(magicfilter1dProgram,"magicfilter1d_denKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
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

void FC_FUNC_(magicfilter1d_t_l,MAGICFILTER1D_T_L)(cl_command_queue *command_queue, size_t *n,size_t *ndat,cl_mem *psi,cl_mem *out){
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
    clSetKernelArg(magicfilter1d_t_kernel_l, i++,sizeof(*n), (void*)n);
    clSetKernelArg(magicfilter1d_t_kernel_l, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(magicfilter1d_t_kernel_l, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(magicfilter1d_t_kernel_l, i++,sizeof(*out), (void*)out);
    clSetKernelArg(magicfilter1d_t_kernel_l, i++,sizeof(float)*block_size_j*(block_size_i+FILTER_WIDTH), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_t_kernel_l, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_t_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %d, %d}\n",globalWorkSize[0],globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %d, %d}\n",localWorkSize[0],localWorkSize[1]);
        exit(1);
    }   
}

void FC_FUNC_(magicfilter1d_den_l,MAGICFILTER1D_DEN_L)(cl_command_queue *command_queue, size_t *n,size_t *ndat,cl_mem *psi,cl_mem *out){
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
    clSetKernelArg(magicfilter1d_den_kernel_l, i++,sizeof(*n), (void*)n);
    clSetKernelArg(magicfilter1d_den_kernel_l, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(magicfilter1d_den_kernel_l, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(magicfilter1d_den_kernel_l, i++,sizeof(*out), (void*)out);
    clSetKernelArg(magicfilter1d_den_kernel_l, i++,sizeof(float)*block_size_j*(block_size_i+FILTER_WIDTH), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_den_kernel_l, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_den_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %d, %d}\n",globalWorkSize[0],globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %d, %d}\n",localWorkSize[0],localWorkSize[1]);
        exit(1);
    }   
}

void FC_FUNC_(magicfilter1d_pot_l,MAGICFILTER1D_POT_L)(cl_command_queue *command_queue, size_t *n,size_t *ndat, cl_mem *psi, float *p, cl_mem *out){
    cl_int ciErrNum;
#if DEBUG
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n: %d, dimension dat: %d, psi: %p, pot: %f, out: %p\n",*command_queue, *n, *ndat, *psi, *p, *out);
#endif
    int FILTER_WIDTH = 16;
    if(*n<FILTER_WIDTH) { fprintf(stderr,"%s %s : matrix is too small!\n", __func__, __FILE__); exit(1);}
    size_t block_size_i=FILTER_WIDTH, block_size_j=256/FILTER_WIDTH;
    while (*n > block_size_i >= 1 && block_size_j > 4)
	{ block_size_i *= 2; block_size_j /= 2;}

    cl_uint i = 0;
    clSetKernelArg(magicfilter1d_pot_kernel_l, i++,sizeof(*n), (void*)n);
    clSetKernelArg(magicfilter1d_pot_kernel_l, i++,sizeof(*ndat), (void*)ndat);
    clSetKernelArg(magicfilter1d_pot_kernel_l, i++,sizeof(*psi), (void*)psi);
    clSetKernelArg(magicfilter1d_pot_kernel_l, i++,sizeof(*p), (void*)p);
    clSetKernelArg(magicfilter1d_pot_kernel_l, i++,sizeof(*out), (void*)out);
    clSetKernelArg(magicfilter1d_pot_kernel_l, i++,sizeof(float)*block_size_j*(block_size_i+FILTER_WIDTH), 0);
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*n), shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, magicfilter1d_pot_kernel_l, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue magicfilter1d_pot_l kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %d, %d}\n",globalWorkSize[0],globalWorkSize[1]);
        fprintf(stderr,"localWorkSize = { %d, %d}\n",localWorkSize[0],localWorkSize[1]);
        exit(1);
    }   
}



