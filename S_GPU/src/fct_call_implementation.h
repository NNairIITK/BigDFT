/**
 * @file fct_call_implementation.h
 * @author Matthieu Ospici
 * 
 * @brief
 * Here are extended all the functions from fct_call.
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


#ifndef FCT_CALL_IMPLEMENTATION_H_
#define FCT_CALL_IMPLEMENTATION_H_

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include "fct_call.h"
#include "init_network.h"
#include "sg_common_def.h"

class fct_call_malloc_gpu : public fct_call
{
public:
#ifdef TRACE_S_GPU
	fct_call_malloc_gpu ( void **dest_, size_t size_, tracer_callbacks_t *_tcal )
			: fct_call (_tcal), dest ( dest_ ), size ( size_ ) {}
#else
	fct_call_malloc_gpu ( void **dest_, size_t size_ )
			: dest ( dest_ ), size ( size_ ) {}
#endif
	virtual void operator() ( int );
	virtual ~fct_call_malloc_gpu() {};

private:
	void **dest;
	size_t size;
};

class fct_call_malloc_pinned : public fct_call
{
public:
#ifdef TRACE_S_GPU
	fct_call_malloc_pinned ( void **dest_, size_t size_, tracer_callbacks_t *_tcal )
			: fct_call (_tcal), dest ( dest_ ), size ( size_ ) {}
#else
	fct_call_malloc_pinned ( void **dest_, size_t size_ )
			: dest ( dest_ ), size ( size_ ) {}
#endif
	virtual void operator() ( int );
	virtual ~fct_call_malloc_pinned() {};

private:
	void **dest;
	size_t size;
};

class fct_call_trsf_CPU_GPU : public fct_call
{
public:
#ifdef TRACE_S_GPU
	fct_call_trsf_CPU_GPU ( const void *src_, GPU_ptr_t dest_, size_t mem_size_, tracer_callbacks_t *_tcal )
			: fct_call ( _tcal ), src ( src_ ), dest ( dest_ ), mem_size ( mem_size_ ) {}
#else
	fct_call_trsf_CPU_GPU ( const void *src_, GPU_ptr_t dest_, size_t mem_size_ )
			: src ( src_ ), dest ( dest_ ), mem_size ( mem_size_ ) {}
#endif
	virtual void operator() ( int );
	virtual ~fct_call_trsf_CPU_GPU() {};

private:
	const void *src;
	GPU_ptr_t dest;
	size_t mem_size;
};

class fct_call_trsf_GPU_CPU : public fct_call
{
public:
#ifdef TRACE_S_GPU
	fct_call_trsf_GPU_CPU ( const GPU_ptr_t src_, void* dest_, size_t mem_size_, tracer_callbacks_t *_tcal )
			: fct_call ( _tcal ), src ( src_ ), dest ( dest_ ), mem_size ( mem_size_ ) {}
#else
	fct_call_trsf_GPU_CPU ( const GPU_ptr_t src_, void* dest_, size_t mem_size_ )
			: src ( src_ ), dest ( dest_ ), mem_size ( mem_size_ ) {}
#endif
	virtual void operator() ( int );
	virtual ~fct_call_trsf_GPU_CPU() {};

private:
	const GPU_ptr_t src;
	void *dest;
	size_t mem_size;
};

class fct_call_memcpy : public fct_call
{
public:
#ifdef TRACE_S_GPU
	fct_call_memcpy ( const void *src_, void* dest_, size_t mem_size_, tracer_callbacks_t *_tcal )
			: fct_call ( _tcal ), src ( src_ ), dest ( dest_ ), mem_size ( mem_size_ ) {}
#else
	fct_call_memcpy ( const void *src_, void* dest_, size_t mem_size_ )
			: src ( src_ ), dest ( dest_ ), mem_size ( mem_size_ ) {}
#endif
	virtual void operator() ( int );
	virtual ~fct_call_memcpy() {};

private:
	const void *src;
	void *dest;
	size_t mem_size;
};

class fct_call_calc_generic : public fct_call
{
public:
#ifdef TRACE_S_GPU
	fct_call_calc_generic ( sg_callback_ptr_t f_call_, void *param_, size_t size_param_, tracer_callbacks_t *_tcal )
			: fct_call ( _tcal ), f_call ( f_call_ ), size_param ( size_param_ )
#else
	fct_call_calc_generic ( sg_callback_ptr_t f_call_, void *param_, size_t size_param_ )
			: f_call ( f_call_ ), size_param ( size_param_ )
#endif
	{
		malloc_and_copy ( param_, size_param_ );
	}
	virtual void operator() ( int );
	virtual ~fct_call_calc_generic();

private:
	void malloc_and_copy ( void* src, size_t size );
	sg_callback_ptr_t f_call;
	void *local_param;
	size_t size_param;
};
#endif
