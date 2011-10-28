#include "OpenCL_wrappers.h"
#include "fft_generator.h"


inline void fft_generated_no_shared_generic(cl_kernel kernel, cl_command_queue command_queue, cl_uint *n,cl_uint *ndat,cl_mem *psi,cl_mem *out,cl_mem *cosi){
    cl_int ciErrNum;
    size_t block_size_i=1, block_size_j=64;
    cl_uint i = 0;
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*n), (void*)n);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*psi), (void*)psi);
    ciErrNum = clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
    if(!use_constant_memory){
      ciErrNum = clSetKernelArg(kernel, i++,sizeof(*cosi), (void*)cosi);
    }
    size_t localWorkSize[] = { block_size_i,block_size_j };
    size_t globalWorkSize[] ={ *n/ *n, shrRoundUp(block_size_j,*ndat)};
    ciErrNum = clEnqueueNDRangeKernel  (command_queue, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
    oclErrorCheck(ciErrNum,"Failed to enqueue fft kernel!");
}



