/**
 * @file cuda_memory_management.cu
 * @author Matthieu Ospici
 * 
 * @brief
 * cuda_memory_management is wrapping cuda function calls for memory management.
 * This wrapping adds error checking when calling these cuda functions.
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

#include <iostream>
#include <stdio.h>

#include "../exceptions.h"

//functions for allocate / deallocate / send memory to GPU
#define CUERR { cudaError_t err; \
 if ((err = cudaGetLastError()) != cudaSuccess) { \
 printf("CUDA error: %s, in file %s:%d\n", cudaGetErrorString(err), __FILE__, __LINE__); }}


void cuda_pinned_malloc ( void **CPU_pointer,
			 size_t nsize) throw (cuda_error)
{
  check<cuda_error> ( cudaMallocHost ( CPU_pointer, nsize ) != cudaSuccess, "CPU pinned allocation", __FILE__, __LINE__ );
}


void cuda_pinned_free ( void *CPU_pointer) throw (cuda_error)
{
  check<cuda_error> ( cudaFreeHost ( CPU_pointer ) != cudaSuccess, "cudaFreeHost", __FILE__, __LINE__ );
  CUERR;
}


void cuda_malloc( void **GPU_pointer, // pointer indicating the GPU address
		  size_t nsize ) throw (cuda_error)
{
  check<cuda_error> ( cudaMalloc ( GPU_pointer, nsize ) != 0, "GPU allocation", __FILE__, __LINE__ );
}


void cuda_free ( void *GPU_pointer) throw (cuda_error)
{
  check<cuda_error> ( cudaFree (GPU_pointer ) != 0, "CUDA free", __FILE__, __LINE__ );
}


void cuda_send ( void *GPU_pointer,
		 void *CPU_pointer,
		 size_t nsize ) throw (cuda_error)
{
  check<cuda_error> ( cudaMemcpy (GPU_pointer, CPU_pointer, nsize, cudaMemcpyHostToDevice )  != 0, "copy host to device", __FILE__, __LINE__ );
}


void cuda_recv ( void *CPU_pointer,
		 void *GPU_pointer,
		 size_t nsize ) throw (cuda_error)
{
  check<cuda_error> ( cudaMemcpy ( CPU_pointer, GPU_pointer, nsize, cudaMemcpyDeviceToHost ), "copy device to host", __FILE__, __LINE__ );
}
