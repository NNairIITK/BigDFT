/**
 * @file   sg_common_def.h
 * @author Matthieu Ospici
 * @date   Fri Jan  8 12:21:13 2010
 * 
 * @brief  Common definition for the S_GPU API
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

#ifndef SG_COMMON_DEF_H_
#define SG_COMMON_DEF_H_

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif
  /** 
   * Callback type used to send kernels to GPU
   * 
   * @typedef sg_callback_ptr_t 
   * 
    */
	typedef void ( *sg_callback_ptr_t ) ( void* );

 /** 
   * A stream pointer
   * 
   * @typedef sg_stream_ptr_t 
   * 
    */
	typedef void *sg_stream_ptr_t;
 /** 
   * A pointer to a GPU memory area
   * 
   * @typedef GPU_ptr_t
   * 
    */
	typedef void *GPU_ptr_t;
 /** 
   * A pointer to a pinned memory area
   * 
   * @typedef pinned_mem_ptr_t
   * 
    */
	typedef void *pinned_mem_ptr_t;

	typedef void ( *sg_callback_tracer_ptr_t ) ( int );


	typedef struct tracer_cbs
	{
		sg_callback_tracer_ptr_t calc_beg;
		sg_callback_tracer_ptr_t calc_end;

		sg_callback_tracer_ptr_t mem_gpu_beg;
		sg_callback_tracer_ptr_t mem_gpu_end;

		sg_callback_tracer_ptr_t memcpy_beg;
		sg_callback_tracer_ptr_t memcpy_end;
	} tracer_callbacks_t;

#ifdef __cplusplus
}
#endif

#endif
