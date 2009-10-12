#include "s_gpu.h"

void sg_init__(bool *GPUshare, bool *useGPU,int *iproc,int *error)
{
  *error = sg_init(GPUshare,useGPU,*iproc);
}


//=========================

void sg_gpu_send_arr(GPU_ptr *dest_GPU_pointer,
		     void **src_CPU_pointer, 
		     int *nsize,
		     int *precision,
		     sg_stream_ptr *stream,
		     int *ierr)
{
  
  *ierr = sg_gpu_send_arr(*dest_GPU_pointer, *src_CPU_pointer, *nsize*(*precision), *stream);
}

//=========================


void sg_gpu_receiv_arr(void **dest_CPU_pointer,
		       GPU_ptr *src_GPU_pointer, 
		       int *nsize,
		       int *precision,
		       sg_stream_ptr *stream,
		       int *ierr)
{
  *ierr = sg_gpu_receiv_arr(*dest_CPU_pointer, *src_GPU_pointer, *nsize*(*precision), *stream);
}


//=========================


void sg_memcpy_f_to_c(void **dest,
		      void *srcFortran,
		      int *nsize,
		      int *precision,
		      sg_stream_ptr *stream,
		      int *ierr)
{
  *ierr = sg_mem_copy(*dest, srcFortran,  *nsize*(*precision),  *stream);
}

//=========================

void sg_memcpy_c_to_f(void *destFortran,
		      void **src,
		      int *nsize,
		      int *precision,
		      sg_stream_ptr *stream,
		      int *ierr)
{
  *ierr = sg_mem_copy(destFortran, *src,  *nsize*(*precision),  *stream);
}

//=========================


void sg_calc(sg_callback_ptr *f_call, void *param, int * size_param, sg_stream_ptr *stream, int *ierr)
{


  *ierr = sg_calc(*f_call, param,  *size_param, *stream);

}
