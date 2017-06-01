#ifndef STREAM_H
#define STREAM_H
#include <CL/cl.h>
#include <stdlib.h>

struct _ocl_stream {
  cl_int reference_counter;
  cl_event event;
  cl_command_queue command_queue;
};
typedef struct _ocl_stream * ocl_stream;


cl_int oclInitStreams(cl_context context);
cl_int oclEndStreams();
ocl_stream oclCreateStream(cl_command_queue command_queue, cl_int *errcode_ret);
cl_int oclReleaseStream(ocl_stream stream);
cl_int oclEnstreamNDRangeKernel(ocl_stream stream, 
                                cl_kernel kernel, 
                                cl_uint work_dim, 
                                const size_t *global_work_offset,
                                const size_t *global_work_size,
                                const size_t *local_work_size);
cl_int oclEnstreamWriteBuffer(ocl_stream stream, cl_mem buffer, size_t offset, size_t cb, const void *ptr);
cl_int oclEnstreamReadBuffer(ocl_stream stream, cl_mem buffer,  size_t offset, size_t cb, void *ptr);
cl_int oclStreamFinish(ocl_stream stream);
cl_int oclStreamFlush(ocl_stream stream);

#endif
