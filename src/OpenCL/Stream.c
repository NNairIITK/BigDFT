//! @file
//!  Stream for OpenCL (??)
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include "Stream.h"
#include <stdio.h>
#include "OpenCL_wrappers.h"

cl_kernel void_kernel;

cl_int oclInitStreams(cl_context context){
  cl_int ciErrNum = CL_SUCCESS;
  const char * stream_program ="__kernel void voidKernel(void) {return;};";
  cl_program streamProgram = clCreateProgramWithSource(context,1,(const char**) &stream_program, NULL, &ciErrNum);
  if( ciErrNum != CL_SUCCESS)
    return ciErrNum;
  ciErrNum = clBuildProgram(streamProgram, 0, NULL, "-cl-mad-enable", NULL, NULL);
  if (ciErrNum != CL_SUCCESS) {
    char cBuildLog[10240];
    clGetProgramBuildInfo(streamProgram, oclGetFirstDev(context), CL_PROGRAM_BUILD_LOG,sizeof(cBuildLog), cBuildLog, NULL );
    fprintf(stderr,"%s\n",cBuildLog);
    return ciErrNum;
  }
  ciErrNum = CL_SUCCESS;
  void_kernel=clCreateKernel(streamProgram,"voidKernel",&ciErrNum);
  if (ciErrNum != CL_SUCCESS)
    return ciErrNum;
  ciErrNum = clReleaseProgram(streamProgram);
  if (ciErrNum != CL_SUCCESS)
    return ciErrNum;
  return CL_SUCCESS;
}

cl_int oclEndStreams(){
  return clReleaseKernel(void_kernel);
}

ocl_stream oclCreateStream(cl_command_queue command_queue, cl_int *errcode_ret) {
  ocl_stream s = (struct _ocl_stream *)malloc(sizeof(struct _ocl_stream));
  if(s == NULL) {
    *errcode_ret = CL_OUT_OF_HOST_MEMORY;
    return NULL;
  }
  size_t localWorksize[] = {32};
  size_t globalWorkSize[] = {32};
  *errcode_ret = clEnqueueNDRangeKernel(command_queue, void_kernel, 1, NULL, globalWorkSize, localWorksize, 0, NULL, &(s->event));
  if(*errcode_ret != CL_SUCCESS) {
    free(s);
    return NULL;
  }
  *errcode_ret = clRetainCommandQueue(command_queue);
  if(*errcode_ret != CL_SUCCESS) {
    free(s);
    return NULL;
  }
  s->command_queue = command_queue;
  s->reference_counter = 1;
  *errcode_ret = CL_SUCCESS;
  return s;
}

cl_int oclReleaseStream( ocl_stream stream) {
  cl_int errcode_ret;
  errcode_ret = clReleaseCommandQueue(stream->command_queue);
  if( errcode_ret != CL_SUCCESS)
    return errcode_ret;
  errcode_ret = clReleaseEvent(stream->event);
  if( errcode_ret != CL_SUCCESS)
    return errcode_ret;
  stream->reference_counter--;
  if( stream->reference_counter==0 )
    free(stream);
  return CL_SUCCESS;
}

cl_int oclEnstreamNDRangeKernel(ocl_stream stream, 
                                     cl_kernel kernel, 
                                     cl_uint work_dim, 
                                     const size_t *global_work_offset,
                                     const size_t *global_work_size,
                                     const size_t *local_work_size)
{
  cl_int errcode_ret;
  cl_event event;
  errcode_ret = clEnqueueNDRangeKernel(stream->command_queue, kernel, work_dim, global_work_offset, global_work_size, local_work_size, 1, &(stream->event), &event);
  if( errcode_ret != CL_SUCCESS)
    return errcode_ret;
  errcode_ret = clReleaseEvent(stream->event);
  if( errcode_ret != CL_SUCCESS)
    return errcode_ret;
  stream->event = event;
  return CL_SUCCESS;
}

cl_int oclStreamFlush(ocl_stream stream) {
  return clFlush(stream->command_queue);
}

cl_int oclStreamFinish(ocl_stream stream) {
  return clWaitForEvents(1, &(stream->event));
}

cl_int oclEnstreamWriteBuffer(ocl_stream stream, cl_mem buffer, size_t offset, size_t cb, const void *ptr) {
  cl_int errcode_ret;
  cl_event event;
  errcode_ret = clEnqueueWriteBuffer(stream->command_queue, buffer, CL_FALSE, offset, cb, ptr, 1, &(stream->event), &(event));
  if( errcode_ret != CL_SUCCESS)
    return errcode_ret;
  errcode_ret = clReleaseEvent(stream->event);
  if( errcode_ret != CL_SUCCESS)
    return errcode_ret;
  stream->event = event;
  return CL_SUCCESS;
}

cl_int oclEnstreamReadBuffer(ocl_stream stream, cl_mem buffer,  size_t offset, size_t cb, void *ptr) {
  cl_int errcode_ret;
  cl_event event;
  errcode_ret = clEnqueueReadBuffer(stream->command_queue, buffer, CL_FALSE, offset, cb, ptr, 1, &(stream->event), &(event));
  if( errcode_ret != CL_SUCCESS)
    return errcode_ret;
  errcode_ret = clReleaseEvent(stream->event);
  if( errcode_ret != CL_SUCCESS)
    return errcode_ret;
  stream->event = event;
  return CL_SUCCESS;
}
