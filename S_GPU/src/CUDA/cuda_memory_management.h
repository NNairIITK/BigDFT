/**
 * @file cuda_memory_management.h
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

#ifndef CUDA_MEMMGMH
#define CUDA_MEMMGMH

#ifdef HAVE_CONFIG
#include <config.h>
#endif

void cuda_pinned_malloc ( void **CPU_pointer,
			  size_t nsize) throw (cuda_error);

void cuda_pinned_free ( void *CPU_pointer) throw (cuda_error);

void cuda_malloc( void **GPU_pointer, // pointer indicating the GPU address
		  size_t nsize ) throw (cuda_error);

void cuda_free ( void *GPU_pointer) throw (cuda_error);

void cuda_send ( void *GPU_pointer,
		 void *CPU_pointer,
		 size_t nsize ) throw (cuda_error);

void cuda_recv ( void *CPU_pointer,
		 void *GPU_pointer,
		 size_t nsize ) throw (cuda_error);

#endif
