#ifndef __sgpuh__
#define __sgpuh__
#include <stddef.h>
#include "sg_common_def.h"

extern "C"
{

  int sg_init(bool *GPUshare, bool *useGPU,int iproc);

  sg_stream_ptr sg_create_stream();
  
  void sg_delete_stream(sg_stream_ptr);
  
  void sg_launch_all_streams();

  int sg_mem_copy(void *dest, const void *src, size_t mem_size, sg_stream_ptr stream);


  int sg_gpu_send_arr(GPU_ptr dest, const void *src, size_t mem_size, sg_stream_ptr stream);
  
  int sg_gpu_receiv_arr(void *dest, const GPU_ptr src, size_t mem_size, sg_stream_ptr stream);

  
  int sg_calc(sg_callback_ptr f_call, void *param,  size_t param_size,sg_stream_ptr stream);
}

#endif
