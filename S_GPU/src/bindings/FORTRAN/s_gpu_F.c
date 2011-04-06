/**
 * @file s_gpu_F.cpp
 * @author Matthieu Ospici
 * 
 * @brief
 * S_GPU is a library written in C++. In order to used by Fortran code, API function
 * names must rewritten with underscore and double underscores added at the end.
 * 
 * @section LICENSE
 * 
 * Copyright (C) 2010 BULL LIG CEA-INAC UJF
 *
 * This file is part of S_GPU library.
 * 
 * S_GPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * S_GPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with S_GPU.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include "s_gpu.h"
#include <stdio.h>

void sg_init__ ( int *GPUshare, int *useGPU, int *iproc, int *nproc_node, int *error )
{
	*error = sg_init ( GPUshare, useGPU, *iproc, *nproc_node );
}

//=========================

void sg_end__()
{
	sg_end();
}


//=========================


void sg_create_stream__ ( sg_stream_ptr_t *stream )
{
	*stream = sg_create_stream();
}

//=========================

void sg_exec_all_streams__()
{
	sg_exec_all_streams();
}
//=========================

void sg_gpu_send_pimem__ ( GPU_ptr_t *dest_GPU_pointer,
                           pinned_mem_ptr_t *src_CPU_pointer,
                           int *nsize,
                           int *precision,
                           sg_stream_ptr_t *stream,
                           int *ierr )
{

	*ierr = sg_gpu_send_pimem ( *dest_GPU_pointer, *src_CPU_pointer, *nsize * ( *precision ), *stream );
}

void sg_gpu_send_mem__( GPU_ptr_t *dest_GPU_pointer,
			void *src_CPU_pointer,
			pinned_mem_ptr_t *buff_pi,
			int *nsize,
			int *precision,
			sg_stream_ptr_t *stream,
			int *ierr )
{

   *ierr = sg_gpu_send_mem ( *dest_GPU_pointer, src_CPU_pointer, *buff_pi,*nsize * ( *precision ), *stream );

 }


//=========================


void sg_gpu_recv_pimem__ ( pinned_mem_ptr_t *dest_CPU_pointer,
                           GPU_ptr_t *src_GPU_pointer,
                           int *nsize,
                           int *precision,
                           sg_stream_ptr_t *stream,
                           int *ierr )
{
	*ierr = sg_gpu_recv_pimem ( *dest_CPU_pointer, *src_GPU_pointer, *nsize * ( *precision ), *stream );


}


void sg_gpu_recv_mem__ ( void *dest_CPU_pointer,
			 GPU_ptr_t *src_GPU_pointer,
			 pinned_mem_ptr_t *buff_pi,
			 int *nsize,
			 int *precision,
			 sg_stream_ptr_t *stream,
			 int *ierr )
{
  *ierr = sg_gpu_recv_mem (dest_CPU_pointer, *src_GPU_pointer,*buff_pi, *nsize * ( *precision ), *stream );


}



//=========================


void sg_memcpy_f_to_c__ ( void **dest,
                          void *srcFortran,
                          int *nsize,
                          int *precision,
                          sg_stream_ptr_t *stream,
                          int *ierr )
{


	*ierr = sg_mem_copy ( *dest, srcFortran,  *nsize * ( *precision ),  *stream );


}

//=========================

void sg_memcpy_c_to_f__ ( void *destFortran,
                          void **src,
                          int *nsize,
                          int *precision,
                          sg_stream_ptr_t *stream,
                          int *ierr )
{


	*ierr = sg_mem_copy ( destFortran, *src,  *nsize * ( *precision ),  *stream );

}

//=========================



void sg_cpu_malloc_pinned__(pinned_mem_ptr_t *pi_mem,
			    int *nsize, 
			    int *precision, 	     
			    int *ierr) 
{
  unsigned int mem_size = (*nsize)*(*precision);
  *ierr = sg_cpu_malloc_pinned(pi_mem,mem_size);
}
			  


void sg_cpu_free_pinned__(pinned_mem_ptr_t *pi_mem,int *ierr)
{
 *ierr = sg_cpu_free_pinned(*pi_mem);
}
  

void sg_gpu_malloc__(GPU_ptr_t *gpu_mem,
		    int *nsize, 
		    int *precision, 	     
		    int *ierr) 

{
    unsigned int mem_size = (*nsize)*(*precision);
    *ierr = sg_gpu_malloc(gpu_mem,mem_size);
}


void sg_gpu_free__ ( GPU_ptr_t *gpu_ptr,int *ierr)
{
*ierr = sg_gpu_free(*gpu_ptr);
}



void sg_send_mem_instantaneously__ ( GPU_ptr_t *GPU_pointer,
				     void *CPU_pointer,
				     int *nsize,
				     int *precision,
				     int *ierr )
{
  unsigned int mem_size = (*nsize)*(*precision);

  *ierr = sg_send_mem_instantaneously(*GPU_pointer, CPU_pointer,mem_size);
}


void sg_recv_mem_instantaneously__ ( void *dest,
				     GPU_ptr_t *gpu_ptr_src,
				     int *nsize,
				     int *precision,
				     int *ierr )
{
  unsigned int mem_size = (*nsize)*(*precision);
  *ierr = sg_recv_mem_instantaneously(dest, *gpu_ptr_src,mem_size);
}


// ================  one underscore support

void sg_init_ ( int *GPUshare, int *useGPU, int *iproc, int *nproc_node, int *error )
{
	*error = sg_init ( GPUshare, useGPU, *iproc, *nproc_node );
}

//=========================

void sg_end_()
{
	sg_end();
}


//=========================


void sg_create_stream_ ( sg_stream_ptr_t *stream )
{
	*stream = sg_create_stream();
}

//=========================

void sg_exec_all_streams_()
{
	sg_exec_all_streams();
}
//=========================

void sg_gpu_send_pimem_ ( GPU_ptr_t *dest_GPU_pointer,
                           pinned_mem_ptr_t *src_CPU_pointer,
                           int *nsize,
                           int *precision,
                           sg_stream_ptr_t *stream,
                           int *ierr )
{

	*ierr = sg_gpu_send_pimem ( *dest_GPU_pointer, *src_CPU_pointer, *nsize * ( *precision ), *stream );
}

void sg_gpu_send_mem_( GPU_ptr_t *dest_GPU_pointer,
			void *src_CPU_pointer,
			pinned_mem_ptr_t *buff_pi,
			int *nsize,
			int *precision,
			sg_stream_ptr_t *stream,
			int *ierr )
{

   *ierr = sg_gpu_send_mem ( *dest_GPU_pointer, src_CPU_pointer, *buff_pi,*nsize * ( *precision ), *stream );

 }


//=========================


void sg_gpu_recv_pimem_ ( pinned_mem_ptr_t *dest_CPU_pointer,
                           GPU_ptr_t *src_GPU_pointer,
                           int *nsize,
                           int *precision,
                           sg_stream_ptr_t *stream,
                           int *ierr )
{
	*ierr = sg_gpu_recv_pimem ( *dest_CPU_pointer, *src_GPU_pointer, *nsize * ( *precision ), *stream );


}


void sg_gpu_recv_mem_ ( void *dest_CPU_pointer,
			GPU_ptr_t *src_GPU_pointer,
			 pinned_mem_ptr_t *buff_pi,
			 int *nsize,
			 int *precision,
			 sg_stream_ptr_t *stream,
			 int *ierr )
{
  *ierr = sg_gpu_recv_mem (dest_CPU_pointer, *src_GPU_pointer,*buff_pi, *nsize * ( *precision ), *stream );


}



//=========================


void sg_memcpy_f_to_c_ ( void **dest,
                          void *srcFortran,
                          int *nsize,
                          int *precision,
                          sg_stream_ptr_t *stream,
                          int *ierr )
{


	*ierr = sg_mem_copy ( *dest, srcFortran,  *nsize * ( *precision ),  *stream );


}

//=========================

void sg_memcpy_c_to_f_ ( void *destFortran,
                          void **src,
                          int *nsize,
                          int *precision,
                          sg_stream_ptr_t *stream,
                          int *ierr )
{


	*ierr = sg_mem_copy ( destFortran, *src,  *nsize * ( *precision ),  *stream );

}

//=========================



void sg_cpu_malloc_pinned_(pinned_mem_ptr_t *pi_mem,
			    int *nsize, 
			    int *precision, 	     
			    int *ierr) 
{
  unsigned int mem_size = (*nsize)*(*precision);
  *ierr = sg_cpu_malloc_pinned(pi_mem,mem_size);
}
			  


void sg_cpu_free_pinned_(pinned_mem_ptr_t *pi_mem,int *ierr)
{
 *ierr = sg_cpu_free_pinned(*pi_mem);
}
  

void sg_gpu_malloc_(GPU_ptr_t *gpu_mem,
		    int *nsize, 
		    int *precision, 	     
		    int *ierr) 

{
    unsigned int mem_size = (*nsize)*(*precision);
    *ierr = sg_gpu_malloc(gpu_mem,mem_size);
}


void sg_gpu_free_ ( GPU_ptr_t *gpu_ptr,int *ierr)
{
*ierr = sg_gpu_free(*gpu_ptr);
}



void sg_send_mem_instantaneously_ ( GPU_ptr_t *GPU_pointer,
				     void *CPU_pointer,
				     int *nsize,
				     int *precision,
				     int *ierr )
{
  unsigned int mem_size = (*nsize)*(*precision);

  *ierr = sg_send_mem_instantaneously(*GPU_pointer, CPU_pointer,mem_size);
}


void sg_recv_mem_instantaneously_ ( void *dest,
				     GPU_ptr_t *gpu_ptr_src,
				     int *nsize,
				     int *precision,
				     int *ierr )
{
  unsigned int mem_size = (*nsize)*(*precision);
  *ierr = sg_recv_mem_instantaneously(dest, *gpu_ptr_src,mem_size);
}
