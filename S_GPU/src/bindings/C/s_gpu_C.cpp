/**
 * @file s_gpu_C.cpp
 * @author Matthieu Ospici
 * 
 * @brief
 * S_GPU is a library written in C++. In order to be used by C code, API functions
 * must be implemented with a extern "C" preamble.
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
#include <vector>
#include <sstream>
#include <s_gpu.h>



#include "exceptions.h"
#include "init_network.h"
#include "read_conf_file.h"
#include "localqueu.h"
#include "manage_cpu_affinity.h"
#include "set_repartition.h"
#include "manage_global_var.h"
#include "fct_call_implementation.h"


#include "CUDA/cuda_memory_management.h"



global_gpu_attach *g_gpu_attach = NULL;

local_network *l = NULL;
localqueu *locq = NULL;

sem_unix *sem_gpu_CALC;
sem_unix *sem_gpu_TRSF;

#ifdef TRACE_S_GPU
tracer_callbacks_t *tracer_callbacks;
#endif

extern "C"
int sg_init ( int *GPUshare, int *useGPU, int iproc, int nproc_node )
{
	*useGPU = 0;
	try
	{
		g_gpu_attach = new global_gpu_attach();
		const char *NAME_FILE = "GPU.config";
		readConfFile read_conf ( NAME_FILE );
		int mpi_tasks_per_node, num_GPU;
		int use_shared;
		//  int iconv_param,iblas_param;
		//read file
		//read_conf.get ( "MPI_TASKS_PER_NODE", &mpi_tasks_per_node );
        mpi_tasks_per_node = nproc_node;
		read_conf.get ( "NUM_GPU", &num_GPU );
		if (num_GPU <  mpi_tasks_per_node){
		  use_shared=1;
		    }
		else
		  {
		  use_shared=0;
		  }
		//read_conf.get ( "USE_SHARED", &use_shared );
		//  read_conf.get("USE_GPU_BLAS",&iblas_param);
		//  read_conf.get("USE_GPU_CONV",&iconv_param);
		manage_cpu_affinity mca ( iproc );
		for ( int i = 0; i < num_GPU; ++i )
		{
			std::ostringstream iss;
			std::string aff;
			iss << "GPU_CPUS_AFF" << "_" << i;
			read_conf.get ( iss.str(), aff );
			mca.add_connexion ( gpu_cpus_connexion ( aff ) );
		}
		set_repartition *set_r;
		if ( use_shared == 1 )
		{
			set_r = new set_repartition_shared ( mpi_tasks_per_node, num_GPU, iproc, g_gpu_attach );
			*GPUshare = 1;
		}
		else
		{
			set_r = new set_repartition_static ( mpi_tasks_per_node, num_GPU, iproc, g_gpu_attach );
			*GPUshare = 0;
		}
		//init node
		l = new local_network ( mpi_tasks_per_node, num_GPU, mca, set_r, iproc );
		if ( iproc == 0 )
			std::cout << "Check card on all nodes...." << std::endl;
		//disable GPU for tasks that not need it
		if ( g_gpu_attach->getIsAttached() )
		{
			//check the card precision, in order to detect error
			//call a fortran function in check_card/check_init.f90
			// checker::runTestOne(); //check only if the card has one GPU...
			*useGPU = 1;
		}
		//print repartition affinity
		if ( iproc == 0 )
			mca.print_affinity_matrix();
		delete set_r; //ugly, to change...
		locq = new localqueu();
		sem_gpu_CALC = l->getSemCalc();
		sem_gpu_TRSF = l->getSemTrsf();
#ifdef TRACE_S_GPU
		//initialization of tracer_callback struct
		tracer_callbacks = new tracer_callbacks_t;
		memset ( tracer_callbacks, NULL, sizeof ( tracer_callbacks_t ) );
#endif
	}

	catch ( synchronization_error& se )
	{
		std::cerr << "*** ERROR(s) DETECTED AT THE INITIALIZATION OF THE INTER-NODE COMMUNICATION SYSTEM ***" << std::endl;
		std::cerr << "ERROR MESSAGE : " << se.what() << std::endl;
		return 1;
	}

	catch ( inter_node_communication_error& ie )
	{
		std::cerr << "*** ERROR(s) DETECTED AT THE INITIALIZATION OF THE INTER-NODE COMMUNICATION SYSTEM ***" << std::endl;
		std::cerr << "ERROR MESSAGE : " << ie.what() << std::endl;
		return 1;
	}


	catch ( read_not_found& re )
	{
		std::cerr << "*** ERROR : INVALID CONFIG FILE. You have to set the number of mpi tasks per node and the number of GPU to use per node ***" << std::endl;
		std::cerr << "Missing information : " << re.what() << std::endl;
		return 1;
	}

	catch ( file_not_found& fe )
	{
		std::cerr << "*** ERROR : CONFIG FILE NOT FOUND" << std::endl;
		std::cerr << "File not found : " << fe.what() << std::endl;
		return 1;
	}

	catch ( check_calc_error& cce )
	{
		std::cerr << "*** ERROR : HARDWARE PROBLEME ON A CARD" << std::endl;
		std::cerr << "We have send calculations to a card and the result was bad. *** Hostname " << cce.what() << "***" << std::endl;
		return 1;
	}

	catch ( std::exception& e )
	{
		std::cerr << "*** ERROR(s) DETECTED AT THE INITIALIZATION OF THE INTER-NODE COMMUNICATION SYSTEM ***" << std::endl;
		std::cerr << "ERROR MESSAGE : " << e.what() << std::endl;
		return 1;
	}

	catch ( ... )
	{
		std::cerr << "** Unexpected exception " << std::endl;
		return 1;
	}

	return 0;
	// std::ostringstream ostr;
	// ostr << "trace_" << iproc;
	// tracer = new trace_exec(ostr.str(),false);
}


extern "C"
void sg_end()
{
	delete g_gpu_attach;
	delete locq;
	delete l;
#ifdef TRACE_S_GPU
	delete tracer_callbacks;
#endif
}

extern "C"
sg_stream_ptr_t sg_create_stream()
{
	gpu_stream *new_stream = new gpu_stream ( l->getCurrGPU() );
	locq->addStream ( new_stream );
	return ( sg_stream_ptr_t ) new_stream;
}

extern "C"
void sg_exec_all_streams()
{
	l->messageLoopNetwork ( *locq );
	//now we can remove the empty streams
	locq->removeStreams();
}

extern "C"
int sg_gpu_send_pimem ( GPU_ptr_t dest, const pinned_mem_ptr_t src, size_t mem_size, sg_stream_ptr_t stream )
{
	fct_call_trsf_CPU_GPU *trsfCPU_GPU =
#ifdef TRACE_S_GPU
	  new fct_call_trsf_CPU_GPU ( src, dest, mem_size, tracer_callbacks );
#else
	  new fct_call_trsf_CPU_GPU ( src, dest, mem_size );
#endif
	( ( gpu_stream* ) stream )->addOp ( trsfCPU_GPU, TRANSF );
	return 0;
}

extern "C"
int sg_gpu_send_mem ( GPU_ptr_t dest, const void* src, pinned_mem_ptr_t buff_pi, size_t mem_size, sg_stream_ptr_t stream )
{
  int err;
  try
    {
      void **loc_buff_pi;
    
      //if pinned pointer is 0, allocate pinned memory
      if(buff_pi == 0)
	{
	  cuda_pinned_malloc (loc_buff_pi, mem_size);
	}
      else
	loc_buff_pi = &buff_pi;

      
      
      //first copy from standard mem to pinned mem
      //no exception can be throwed on this 2 folowing functions, so no memory leak can appear...
      sg_mem_copy ( *loc_buff_pi, src,  mem_size, stream ); 
      sg_gpu_send_pimem ( dest, *loc_buff_pi, mem_size, stream );
      
      //delete pinned if allocated...
      if(buff_pi == 0)
	{
	  cuda_pinned_free(*loc_buff_pi);
	}
    }
  catch(std::exception &e)
    {
      std::cerr << "*** CUDA ERROR DETECTED" << std::endl;
      std::cerr << "ERROR MESSAGE : " << e.what() << std::endl;
      err = 1;
    }

  return err;

  
}

extern "C"
int sg_gpu_recv_pimem ( pinned_mem_ptr_t dest, const GPU_ptr_t src, size_t mem_size, sg_stream_ptr_t stream )
{
	fct_call_trsf_GPU_CPU *trsfGPU_CPU =
#ifdef TRACE_S_GPU
	  new fct_call_trsf_GPU_CPU ( src, dest, mem_size, tracer_callbacks );
#else
	  new fct_call_trsf_GPU_CPU ( src, dest, mem_size );
#endif
	( ( gpu_stream* ) stream )->addOp ( trsfGPU_CPU, TRANSF );
	return 0;
}

extern "C"
int sg_gpu_recv_mem ( void *dest, const GPU_ptr_t src, pinned_mem_ptr_t buff_pi, size_t mem_size, sg_stream_ptr_t stream )
{
	//copy to pinned mem
	sg_gpu_recv_pimem ( buff_pi,  src, mem_size, stream );
	//copy from pinned to standard
	sg_mem_copy ( dest,buff_pi, mem_size, stream );
	return 0;
}

extern "C"
int sg_mem_copy ( void *dest, const void *src, size_t mem_size, sg_stream_ptr_t stream )
{
	fct_call_memcpy *trsf_memcpy =
#ifdef TRACE_S_GPU
	  new fct_call_memcpy ( src, dest, mem_size, tracer_callbacks );
#else
	  new fct_call_memcpy ( src, dest, mem_size );
#endif
	( ( gpu_stream* ) stream )->addOp ( trsf_memcpy, TRANSF );
	return 0;
}

extern "C"
int sg_malloc_gpu ( void **dest, size_t size, sg_stream_ptr_t stream )
{
	fct_call_malloc_gpu *malloc_gpu =
#ifdef TRACE_S_GPU
	  new fct_call_malloc_gpu ( dest, size, tracer_callbacks );
#else
	  new fct_call_malloc_gpu ( dest, size );
#endif
	( ( gpu_stream* ) stream )->addOp ( malloc_gpu, TRANSF );
	return 0;
}

extern "C"
int sg_cpu_malloc_pinned(pinned_mem_ptr_t *pi_mem,  size_t size)
{
  int err = 0;
  try
    {
      cuda_pinned_malloc(pi_mem,size);
    }
  catch(std::exception &e)
    {
      std::cerr << "*** CUDA ERROR DETECTED" << std::endl;
      std::cerr << "ERROR MESSAGE : " << e.what() << std::endl;
      err = 1;
    }

  return err;
}

extern "C"
int sg_cpu_free_pinned(pinned_mem_ptr_t pi_mem)
{
  int err = 0;
  try
    {
      cuda_pinned_free( pi_mem);
    }
  catch(std::exception &e)
    {
      std::cerr << "*** CUDA ERROR DETECTED" << std::endl;
      std::cerr << "ERROR MESSAGE : " << e.what() << std::endl;
	err = 1;
    }
  
  return err;
}


extern "C"
int sg_gpu_malloc ( GPU_ptr_t *gpu_ptr, size_t size)
{
  int err = 0;
  try
    {
	cuda_malloc(gpu_ptr,size);
    }
  catch(std::exception &e)
    {
      std::cerr << "*** CUDA ERROR DETECTED" << std::endl;
	std::cerr << "ERROR MESSAGE : " << e.what() << std::endl;
	err = 1;
    }

  return err;
}

extern "C"
int sg_gpu_free ( GPU_ptr_t gpu_ptr)
{
  int err = 0;
  try
    {
      cuda_free(gpu_ptr);
    }
  catch(std::exception &e)
    {
      std::cerr << "*** CUDA ERROR DETECTED" << std::endl;
      std::cerr << "ERROR MESSAGE : " << e.what() << std::endl;
      err = 1;
    }
  
  return err;
}

extern "C"
int sg_send_mem_instantaneously(GPU_ptr_t gpu_ptr_dest, void *src,size_t size)
{
  int err = 0;
  try
    {
      cuda_send(gpu_ptr_dest,src,size);
    }
  catch(std::exception &e)
    {
      std::cerr << "*** CUDA ERROR DETECTED" << std::endl;
      std::cerr << "ERROR MESSAGE : " << e.what() << std::endl;
      err = 1;
    }
  
  return err;
}

extern "C"
int sg_recv_mem_instantaneously(void *dest,GPU_ptr_t gpu_ptr_src,size_t size)
{
  int err = 0;
  try
    {
      cuda_recv(dest,gpu_ptr_src,size);
    }
  catch(std::exception &e)
    {
      std::cerr << "*** CUDA ERROR DETECTED" << std::endl;
      std::cerr << "ERROR MESSAGE : " << e.what() << std::endl;
      err = 1;
    }
  
  return err;
}


extern "C"
int sg_gpu_send_calc ( sg_callback_ptr_t f_call, void *param, size_t size_param, sg_stream_ptr_t stream )
{
	fct_call_calc_generic *calc_generic =
#ifdef TRACE_S_GPU
	  new fct_call_calc_generic ( f_call, param, size_param, tracer_callbacks );
#else
	  new fct_call_calc_generic ( f_call, param, size_param );
#endif
	( ( gpu_stream* ) stream )->addOp ( calc_generic, CALC );
	return 0;
	}








#ifdef TRACE_S_GPU
extern "C"
void add_tracer_callbacks ( tracer_callbacks_t *tcal )
{
	//copy user struct to the s_gpu struct
	memcpy ( tracer_callbacks, tcal, sizeof ( tracer_callbacks_t ) );
//  tracer_callbacks = tcal;
//  tracer_callbacks->calc_beg(1);
}
#endif
