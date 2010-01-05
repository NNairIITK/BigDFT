#include "s_gpu.h"
#include <config.h>

void FC_FUNC(sg_init,SG_INIT)(int *GPUshare, int *useGPU,int *iproc,int *error)
{
  *error = sg_init(GPUshare,useGPU,*iproc);
}

//=========================

void FC_FUNC(sg_end,SG_END)()
{
  sg_end();
}


//=========================


void FC_FUNC(sg_create_stream,SG_CREATE_STREAM)(sg_stream_ptr *stream)
{
  *stream = sg_create_stream();
}

//=========================

void FC_FUNC(sg_launch_all_streams,SG_LAUNCH_ALL_STREAMS)()
{
  sg_launch_all_streams();
}
//=========================

void FC_FUNC(sg_gpu_pi_send,SG_GPU_PI_SEND)(GPU_ptr *dest_GPU_pointer,
		      void **src_CPU_pointer, 
		      int *nsize,
		      int *precision,
		      sg_stream_ptr *stream,
		      int *ierr)
{
  
  *ierr = sg_gpu_send_arr(*dest_GPU_pointer, *src_CPU_pointer, *nsize*(*precision), *stream);
}

//=========================


void FC_FUNC(sg_gpu_pi_recv,SG_GPU_PI_RECV)(void **dest_CPU_pointer,
		       GPU_ptr *src_GPU_pointer, 
		       int *nsize,
		       int *precision,
		       sg_stream_ptr *stream,
		       int *ierr)
{
   *ierr = sg_gpu_receiv_arr(*dest_CPU_pointer, *src_GPU_pointer, *nsize*(*precision), *stream);


}


//=========================


void FC_FUNC(sg_memcpy_f_to_c,SG_MEMCPY_F_TO_C)(void **dest,
			void *srcFortran,
			int *nsize,
			int *precision,
			sg_stream_ptr *stream,
			int *ierr)
{


    *ierr = sg_mem_copy(*dest, srcFortran,  *nsize*(*precision),  *stream);
 

}

//=========================

void FC_FUNC(sg_memcpy_c_to_f,SG_MEMCPY_C_TO_F)(void *destFortran,
			void **src,
			int *nsize,
			int *precision,
			sg_stream_ptr *stream,
			int *ierr)
{


      *ierr = sg_mem_copy(destFortran, *src,  *nsize*(*precision),  *stream);





 

 



}

//=========================


/*void FC_FUNC(sg_calc,SG_CALC)(sg_callback_ptr *f_call, void **param, sg_stream_ptr *stream, int *ierr)
{


  *ierr = sg_calc(*f_call, *param, *stream);

  }*/
