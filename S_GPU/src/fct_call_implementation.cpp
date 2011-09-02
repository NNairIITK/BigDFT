/**
 * @file fct_call_implementation.cpp
 * @author Matthieu Ospici
 * 
 * @brief
 * Implementation of the extended functions.
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
#ifdef VERBOSE_S_GPU
#include <stdio.h>
#endif

#include "fct_call_implementation.h"
#include "CUDA/cudafct.h"

void fct_call_malloc_pinned::operator() ( int tID )
{
#ifdef VERBOSE_S_GPU
  std::cout << ">> start in " << __FILE__ << ":" << "malloc_pinned" << std::endl;
#endif
  c_cudaMallocHost(dest, size);
#ifdef VERBOSE_S_GPU
  std::cout << "<< end in " << __FILE__ << ":" << "malloc_pinned" << std::endl;
#endif
}

void fct_call_malloc_gpu::operator() ( int tID )
{
#ifdef VERBOSE_S_GPU
  std::cout << ">> start in " << __FILE__ << ":" << "malloc_gpu" << std::endl;
#endif
  c_cudaMalloc(dest, size);
#ifdef VERBOSE_S_GPU
  std::cout << "<< end in " << __FILE__ << ":" << "malloc_gpu" << std::endl;
#endif
}

void fct_call_trsf_CPU_GPU::operator() ( int curGPU )
{
#ifdef TRACE_S_GPU
	if ( tcal->mem_gpu_beg )
		tcal->mem_gpu_beg ( curGPU );
#endif

#ifdef VERBOSE_S_GPU
  std::cout << ">> start in " << __FILE__ << ":" << "trsf_cpu_gpu" << std::endl;
#endif
	c_cuda_gpu_send_pi ( dest, src, mem_size );
#ifdef VERBOSE_S_GPU
  std::cout << "<< end in " << __FILE__ << ":" << "trsf_cpu_gpu" << std::endl;
#endif

#ifdef TRACE_S_GPU
	if ( tcal->mem_gpu_end )
		tcal->mem_gpu_end ( curGPU );
#endif
}

void fct_call_trsf_GPU_CPU::operator() ( int curGPU )
{
#ifdef TRACE_S_GPU
	if ( tcal->mem_gpu_beg )
		tcal->mem_gpu_beg ( curGPU );
#endif

#ifdef VERBOSE_S_GPU
  std::cout << ">> start in " << __FILE__ << ":" << "trsf_gpu_cpu" << std::endl;
#endif
	c_cuda_gpu_recv_pi ( dest, src, mem_size );
#ifdef VERBOSE_S_GPU
  std::cout << "<< end in " << __FILE__ << ":" << "trsf_gpu_cpu" << std::endl;
#endif

#ifdef TRACE_S_GPU
	if ( tcal->mem_gpu_end )
		tcal->mem_gpu_end ( curGPU );
#endif
}

void fct_call_memcpy::operator() ( int curGPU )
{
#ifdef TRACE_S_GPU
	if ( tcal->memcpy_beg )
		tcal->memcpy_beg ( curGPU );
#endif

#ifdef VERBOSE_S_GPU
  std::cout << ">> start in " << __FILE__ << ":" << "memcpy" << std::endl;
#endif
	memcpy ( dest, src, mem_size );
#ifdef VERBOSE_S_GPU
  std::cout << "<< end in " << __FILE__ << ":" << "memcpy" << std::endl;
#endif

#ifdef TRACE_S_GPU
	if ( tcal->memcpy_end )
		tcal->memcpy_end ( curGPU );
#endif
}

//================= fct_call_calc_generic ======================

void fct_call_calc_generic::malloc_and_copy ( void* src, size_t size )
{
#ifdef VERBOSE_S_GPU
  std::cout << ">> start in " << __FILE__ << ":" << "malloc_and_copy" << std::endl;
#endif
	local_param =  malloc ( size );
	memcpy ( local_param, src, size );
#ifdef VERBOSE_S_GPU
  std::cout << "<< end in " << __FILE__ << ":" << "malloc_and_copy" << std::endl;
#endif
}

fct_call_calc_generic::~fct_call_calc_generic()
{
#ifdef VERBOSE_S_GPU
  std::cout << ">> start in " << __FILE__ << ":" << "~calc_generic" << std::endl;
#endif
	free ( local_param );
#ifdef VERBOSE_S_GPU
  std::cout << "<< end in " << __FILE__ << ":" << "~calc_generic" << std::endl;
#endif
}

void fct_call_calc_generic::operator() ( int curGPU )
{
#ifdef TRACE_S_GPU
	if ( tcal->calc_beg )
		tcal->calc_beg ( curGPU );
#endif

#ifdef VERBOSE_S_GPU
  std::cout << ">> start in " << __FILE__ << ":" << "calc_generic" << std::endl;
#endif
	( *f_call ) ( local_param );
#ifdef VERBOSE_S_GPU
  std::cout << "<< end in " << __FILE__ << ":" << "calc_generic" << std::endl;
#endif

#ifdef TRACE_S_GPU
	if ( tcal->calc_end )
		tcal->calc_end ( curGPU );
#endif
}
