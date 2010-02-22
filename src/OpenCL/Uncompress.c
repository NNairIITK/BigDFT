#include "OpenCL_wrappers.h"
#include "Uncompress.h"

cl_kernel uncompress_coarse_kernel_l;
cl_kernel uncompress_fine_kernel_l;
cl_kernel compress_coarse_kernel_l;
cl_kernel compress_fine_kernel_l;
cl_kernel uncompress_coarse_kernel_d;
cl_kernel uncompress_fine_kernel_d;
cl_kernel compress_coarse_kernel_d;
cl_kernel compress_fine_kernel_d;

void build_uncompress_kernels(cl_context * context){
    cl_int ciErrNum = CL_SUCCESS;
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
    uncompress_coarse_kernel_d=clCreateKernel(uncompressProgram,"uncompress_coarseKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    uncompress_fine_kernel_d=clCreateKernel(uncompressProgram,"uncompress_fineKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    uncompress_coarse_kernel_l=clCreateKernel(uncompressProgram,"uncompress_coarseKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    uncompress_fine_kernel_l=clCreateKernel(uncompressProgram,"uncompress_fineKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = clReleaseProgram(uncompressProgram);
    oclErrorCheck(ciErrNum,"Failed to release program!");

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
//    FILE * f = fopen("compress.bin","w");
//    size_t prog_size;
//    clGetProgramInfo( compressProgram,CL_PROGRAM_BINARY_SIZES,sizeof(prog_size),&prog_size,NULL);
//    printf("binary_size %lu\n",(long unsigned int) prog_size);
//    char * bin[1];
//    bin[0] = (char*)malloc(prog_size*sizeof(char));
//    clGetProgramInfo( compressProgram,CL_PROGRAM_BINARIES,prog_size,bin,NULL);
//    fwrite(bin[0],sizeof(char),prog_size,f);
//    fclose(f);
    ciErrNum = CL_SUCCESS;
    compress_coarse_kernel_d=clCreateKernel(compressProgram,"compress_coarseKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    compress_fine_kernel_d=clCreateKernel(compressProgram,"compress_fineKernel_d",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    compress_coarse_kernel_l=clCreateKernel(compressProgram,"compress_coarseKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = CL_SUCCESS;
    compress_fine_kernel_l=clCreateKernel(compressProgram,"compress_fineKernel_l",&ciErrNum);
    oclErrorCheck(ciErrNum,"Failed to create kernel!");
    ciErrNum = clReleaseProgram(compressProgram);
    oclErrorCheck(ciErrNum,"Failed to release program!");
}

void FC_FUNC_(uncompress_d,UNCOMPRESS_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,
                                       cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                       cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
				       cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi_out) {
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n1: %lu, dimension n2: %lu, dimension n3: %lu,\n\
nseg_c: %lu, nvctr_c: %lu, keyg_c: %p, keyv_c: %p,\n\
nseg_f: %lu, nvctr_f: %lu, keyg_f: %p, keyv_f: %p,\n\
psi_c: %p, psi_f: %p, psi_out: %p\n",
*command_queue, (long unsigned)*n1, (long unsigned)*n2, (long unsigned)*n3,
(long unsigned)*nseg_c, (long unsigned)*nvctr_c, *keyg_c, *keyv_c,
(long unsigned)*nseg_f, (long unsigned)*nvctr_f, *keyg_f, *keyv_f,
*psi_c, *psi_f, *psi_out);
#endif
    size_t block_size_i=64;
    cl_uint full_size = *n1 * *n2 * *n3 * 8;
    double init = 0.0;
    while ( block_size_i > full_size ) block_size_i /= 2;
    cl_uint i = 0;
    clSetKernelArg(v_initialize_kernel_d, i++, sizeof(full_size), (void*)&full_size);
    clSetKernelArg(v_initialize_kernel_d, i++, sizeof(*psi_out), (void*)psi_out);
    clSetKernelArg(v_initialize_kernel_d, i++, sizeof(init), (void*)&init);
    size_t localWorkSize[] = { block_size_i };
    size_t globalWorkSize[] ={ shrRoundUp(block_size_i, full_size) };
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, v_initialize_kernel_d, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue v_initialize_d kernel!\n",ciErrNum);
        exit(1);
    } 

    block_size_i=64;
    while ( block_size_i > *nvctr_c ) block_size_i /= 2;
    i = 0;
    clSetKernelArg(uncompress_coarse_kernel_d, i++, sizeof(*n1), (void*)n1);
    clSetKernelArg(uncompress_coarse_kernel_d, i++, sizeof(*n2), (void*)n2);
    clSetKernelArg(uncompress_coarse_kernel_d, i++, sizeof(*n3), (void*)n3);
    clSetKernelArg(uncompress_coarse_kernel_d, i++, sizeof(*nseg_c), (void*)nseg_c);
    clSetKernelArg(uncompress_coarse_kernel_d, i++, sizeof(*nvctr_c), (void*)nvctr_c);
    clSetKernelArg(uncompress_coarse_kernel_d, i++, sizeof(*keyg_c), (void*)keyg_c);
    clSetKernelArg(uncompress_coarse_kernel_d, i++, sizeof(*keyv_c), (void*)keyv_c);
    clSetKernelArg(uncompress_coarse_kernel_d, i++, sizeof(*psi_c), (void*)psi_c);
    clSetKernelArg(uncompress_coarse_kernel_d, i++, sizeof(*psi_out), (void*)psi_out);
    localWorkSize[0]= block_size_i ;
    globalWorkSize[0] = shrRoundUp(block_size_i, *nvctr_c) ;
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, uncompress_coarse_kernel_d, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue uncompress_coarse_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu }\n",(long unsigned)globalWorkSize[0]);
        fprintf(stderr,"localWorkSize = { %lu }\n",(long unsigned)localWorkSize[0]);
        exit(1);
    }

    block_size_i=64;
//    while ( block_size_i > *nvctr_f ) block_size_i /= 2;
    i=0;
    clSetKernelArg(uncompress_fine_kernel_d, i++, sizeof(*n1), (void*)n1);
    clSetKernelArg(uncompress_fine_kernel_d, i++, sizeof(*n2), (void*)n2);
    clSetKernelArg(uncompress_fine_kernel_d, i++, sizeof(*n3), (void*)n3);
    clSetKernelArg(uncompress_fine_kernel_d, i++, sizeof(*nseg_f), (void*)nseg_f);
    clSetKernelArg(uncompress_fine_kernel_d, i++, sizeof(*nvctr_f), (void*)nvctr_f);
    clSetKernelArg(uncompress_fine_kernel_d, i++, sizeof(*keyg_f), (void*)keyg_f);
    clSetKernelArg(uncompress_fine_kernel_d, i++, sizeof(*keyv_f), (void*)keyv_f);
    clSetKernelArg(uncompress_fine_kernel_d, i++, sizeof(*psi_f), (void*)psi_f);
    clSetKernelArg(uncompress_fine_kernel_d, i++, sizeof(*psi_out), (void*)psi_out);
    clSetKernelArg(uncompress_fine_kernel_d, i++, sizeof(double)*block_size_i*7, NULL);
    localWorkSize[0]= block_size_i ;
    globalWorkSize[0] = shrRoundUp(block_size_i, *nvctr_f) ;
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, uncompress_fine_kernel_d, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue uncompress_fine_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu }\n",(long unsigned)globalWorkSize[0]);
        fprintf(stderr,"localWorkSize = { %lu }\n",(long unsigned)localWorkSize[0]);
        exit(1);
    }

}

void FC_FUNC_(uncompress_l,UNCOMPRESS_L)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,
                                       cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                       cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
				       cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi_out) {
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n1: %lu, dimension n2: %lu, dimension n3: %lu,\n\
nseg_c: %lu, nvctr_c: %lu, keyg_c: %p, keyv_c: %p,\n\
nseg_f: %lu, nvctr_f: %lu, keyg_f: %p, keyv_f: %p,\n\
psi_c: %p, psi_f: %p, psi_out: %p\n",
*command_queue, (long unsigned)*n1, (long unsigned)*n2, (long unsigned)*n3,
(long unsigned)*nseg_c, (long unsigned)*nvctr_c, *keyg_c, *keyv_c,
(long unsigned)*nseg_f, (long unsigned)*nvctr_f, *keyg_f, *keyv_f,
*psi_c, *psi_f, *psi_out);
#endif
    size_t block_size_i=64;
    cl_uint full_size = *n1 * *n2 * *n3 * 8;
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
        fprintf(stderr,"globalWorkSize = { %lu }\n",(long unsigned)globalWorkSize[0]);
        fprintf(stderr,"localWorkSize = { %lu }\n",(long unsigned)localWorkSize[0]);
        exit(1);
    }


    block_size_i=64;
    while ( block_size_i > *nvctr_f ) block_size_i /= 2;
    i=0;
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
        fprintf(stderr,"globalWorkSize = { %lu }\n",(long unsigned)globalWorkSize[0]);
        fprintf(stderr,"localWorkSize = { %lu }\n",(long unsigned)localWorkSize[0]);
        exit(1);
    }

}

void FC_FUNC_(compress_d,COMPRESS_D)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,
                                     cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                     cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                     cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi) {
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n1: %lu, dimension n2: %lu, dimension n3: %lu,\n\
nseg_c: %lu, nvctr_c: %lu, keyg_c: %p, keyv_c: %p,\n\
nseg_f: %lu, nvctr_f: %lu, keyg_f: %p, keyv_f: %p,\n\
psi_c: %p, psi_f: %p, psi: %p\n",
*command_queue, (long unsigned)*n1, (long unsigned)*n2, (long unsigned)*n3,
(long unsigned)*nseg_c, (long unsigned)*nvctr_c, *keyg_c, *keyv_c,
(long unsigned)*nseg_f, (long unsigned)*nvctr_f, *keyg_f, *keyv_f,
*psi_c, *psi_f, *psi);
#endif
    size_t block_size_i=64;
    while ( block_size_i > *nvctr_c ) block_size_i /= 2;
    cl_uint i = 0;
    clSetKernelArg(compress_coarse_kernel_d, i++, sizeof(*n1), (void*)n1);
    clSetKernelArg(compress_coarse_kernel_d, i++, sizeof(*n2), (void*)n2);
    clSetKernelArg(compress_coarse_kernel_d, i++, sizeof(*n3), (void*)n3);
    clSetKernelArg(compress_coarse_kernel_d, i++, sizeof(*nseg_c), (void*)nseg_c);
    clSetKernelArg(compress_coarse_kernel_d, i++, sizeof(*nvctr_c), (void*)nvctr_c);
    clSetKernelArg(compress_coarse_kernel_d, i++, sizeof(*keyg_c), (void*)keyg_c);
    clSetKernelArg(compress_coarse_kernel_d, i++, sizeof(*keyv_c), (void*)keyv_c);
    clSetKernelArg(compress_coarse_kernel_d, i++, sizeof(*psi_c), (void*)psi_c);
    clSetKernelArg(compress_coarse_kernel_d, i++, sizeof(*psi), (void*)psi);
    size_t localWorkSize[]= { block_size_i };
    size_t globalWorkSize[]={ shrRoundUp(block_size_i, *nvctr_c) } ;
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, compress_coarse_kernel_d, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue compress_coarse_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu }\n",(long unsigned)globalWorkSize[0]);
        fprintf(stderr,"localWorkSize = { %lu }\n",(long unsigned)localWorkSize[0]);
        exit(1);
    }

    i = 0;
    block_size_i=64;
//    while ( block_size_i > *nvctr_f ) block_size_i /= 2;
    clSetKernelArg(compress_fine_kernel_d, i++, sizeof(*n1), (void*)n1);
    clSetKernelArg(compress_fine_kernel_d, i++, sizeof(*n2), (void*)n2);
    clSetKernelArg(compress_fine_kernel_d, i++, sizeof(*n3), (void*)n3);
    clSetKernelArg(compress_fine_kernel_d, i++, sizeof(*nseg_f), (void*)nseg_f);
    clSetKernelArg(compress_fine_kernel_d, i++, sizeof(*nvctr_f), (void*)nvctr_f);
    clSetKernelArg(compress_fine_kernel_d, i++, sizeof(*keyg_f), (void*)keyg_f);
    clSetKernelArg(compress_fine_kernel_d, i++, sizeof(*keyv_f), (void*)keyv_f);
    clSetKernelArg(compress_fine_kernel_d, i++, sizeof(*psi_f), (void*)psi_f);
    clSetKernelArg(compress_fine_kernel_d, i++, sizeof(*psi), (void*)psi);
    clSetKernelArg(compress_fine_kernel_d, i++, sizeof(double)*block_size_i*7, NULL);
    localWorkSize[0]= block_size_i ;
    globalWorkSize[0] = shrRoundUp(block_size_i, *nvctr_f) ;
    ciErrNum = clEnqueueNDRangeKernel  (*command_queue, compress_fine_kernel_d, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        fprintf(stderr,"Error %d: Failed to enqueue compress_fine_d kernel!\n",ciErrNum);
        fprintf(stderr,"globalWorkSize = { %lu }\n",(long unsigned)globalWorkSize[0]);
        fprintf(stderr,"localWorkSize = { %lu }\n",(long unsigned)localWorkSize[0]);
        exit(1);
    }

}

void FC_FUNC_(compress_l,COMPRESS_L)(cl_command_queue *command_queue, cl_uint *n1, cl_uint *n2, cl_uint *n3,
                                     cl_uint *nseg_c, cl_uint *nvctr_c, cl_mem *keyg_c, cl_mem *keyv_c, 
                                     cl_uint *nseg_f, cl_uint *nvctr_f, cl_mem *keyg_f, cl_mem *keyv_f,
                                     cl_mem *psi_c, cl_mem *psi_f, cl_mem * psi) {
    cl_int ciErrNum;
#if DEBUG     
    printf("%s %s\n", __func__, __FILE__);
    printf("command queue: %p, dimension n1: %lu, dimension n2: %lu, dimension n3: %lu,\n\
nseg_c: %lu, nvctr_c: %lu, keyg_c: %p, keyv_c: %p,\n\
nseg_f: %lu, nvctr_f: %lu, keyg_f: %p, keyv_f: %p,\n\
psi_c: %p, psi_f: %p, psi: %p\n",
*command_queue, (long unsigned)*n1, (long unsigned)*n2, (long unsigned)*n3,
(long unsigned)*nseg_c, (long unsigned)*nvctr_c, *keyg_c, *keyv_c,
(long unsigned)*nseg_f, (long unsigned)*nvctr_f, *keyg_f, *keyv_f,
*psi_c, *psi_f, *psi);
#endif
    size_t block_size_i=64;
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
        fprintf(stderr,"globalWorkSize = { %lu }\n",(long unsigned)globalWorkSize[0]);
        fprintf(stderr,"localWorkSize = { %lu }\n",(long unsigned)localWorkSize[0]);
        exit(1);
    }
    i=0;
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
        fprintf(stderr,"globalWorkSize = { %lu }\n",(long unsigned)globalWorkSize[0]);
        fprintf(stderr,"localWorkSize = { %lu }\n",(long unsigned)localWorkSize[0]);
        exit(1);
    }

}

void clean_uncompress_kernels(){
  clReleaseKernel(uncompress_coarse_kernel_l);
  clReleaseKernel(uncompress_fine_kernel_l);
  clReleaseKernel(compress_coarse_kernel_l);
  clReleaseKernel(compress_fine_kernel_l);
  clReleaseKernel(uncompress_coarse_kernel_d);
  clReleaseKernel(uncompress_fine_kernel_d);
  clReleaseKernel(compress_coarse_kernel_d);
  clReleaseKernel(compress_fine_kernel_d);
}
