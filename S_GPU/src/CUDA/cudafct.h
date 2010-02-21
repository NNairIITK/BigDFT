/**
 * @file cudafct.h
 * @author Matthieu Ospici
 * 
 * @brief
 * cudafct is wrapping cuda function calls. This wrapping adds error checking when
 * calling these cuda functions.
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

#ifndef   	CUDAFCT_H_
#define   	CUDAFCT_H_

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

int c_cudaGetDeviceCount ( int *s_gpuCount );
int c_cudaSetDevice ( int device );

int c_cudaMalloc ( void**, size_t memsize );

int c_cudaMallocHost ( void**, size_t memsize );
int c_cudaMallocHostFloat ( float**, size_t memsize );

int c_cuda_gpu_send_pi ( void *dest, const void *src,  size_t memByte );

int c_cuda_gpu_recv_pi ( void *dest, const void *src,  size_t memByte );

int c_cuda_setdevice ( int device );
int c_cuda_get_device ( int *dev );

int c_cuda_setdevice_ctx ( int device );
int c_cuda_ctxpopcur ( void *contex );
int c_cuda_ctxpushcur ( void *contex );
#ifdef __cplusplus
}
#endif

#endif 	    /* !CUDAFCT_H_ */
