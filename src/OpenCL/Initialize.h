#ifndef INITIALIZE_H
#define INITIALIZE_H

char * c_initialize_program="\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
__kernel void c_initializeKernel_d(size_t n, __global const double * x_in, __global double * y_in, double c) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
y_in[ig] = x_in[ig] * c;\n\
};\n\
__kernel void v_initializeKernel_d(size_t n, __global double * y_in, double v) {\n\
size_t ig = get_global_id(0);\n\
ig = get_group_id(0) == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;\n\
y_in[ig] = v;\n\
};\n\
";

void inline c_initialize_generic(cl_kernel kernel, cl_command_queue *command_queue, cl_uint *ndat, cl_mem *in, cl_mem *out, double *c) {
  cl_int ciErrNum;
  size_t block_size_i=64;
  assert(*ndat>=block_size_i);
  cl_uint i=0;
  clSetKernelArg(kernel, i++,sizeof(*ndat), (void*)ndat);
  clSetKernelArg(kernel, i++,sizeof(*in), (void*)in);
  clSetKernelArg(kernel, i++,sizeof(*out), (void*)out);
  clSetKernelArg(kernel, i++,sizeof(*c), (void*)c);
  size_t localWorkSize[] = { block_size_i };
  size_t globalWorkSize[] ={ shrRoundUp(block_size_i,*ndat) };
  ciErrNum = clEnqueueNDRangeKernel  (*command_queue, kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
  oclErrorCheck(ciErrNum,"Failed to enqueue c_initialize kernel!");
}

#endif
