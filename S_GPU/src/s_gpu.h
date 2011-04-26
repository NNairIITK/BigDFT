/**
 * @file   s_gpu.h
 * @author Matthieu Ospici
 * @date   Fri Jan  8 12:21:13 2010
 * @brief  C interface for the S_GPU library. 
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
 * 
 */
#ifndef S_GPU_H_
#define S_GPU_H_

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include <stddef.h>

#include "sg_common_def.h"

/**
 * @mainpage
 *<b> Overview of S GPU</b>


S_GPU is a software stack between the application and the GPU itself. Each process runs its operations on the GPUs using the functions provided by S_GPU. If multiple GPUs are available, S_GPU will automatically assign to each process one GPU, without any programmer intervention. 

The purpose of S_GPU is to provide a kind of GPU virtualization. When  some processes are launch (typically MPI processes), each process see only one @b private GPU. Thus, all processes can do their operations on the GPU regardless the other processes. The S_GPU library will do the necessary works in order to share the GPU(s) between the CPUs.

S_GPU is built around the concept of stream. A S_GPU stream can be viewed as a list of dependencies. If we put on the stream the three operations (A, B, C), this indicates that the operation B requires the operation A is finished for being executed, the same is true for C, which requires B.


During a call to a S_GPU function, as a memory transfer, the command is put into a stream and it is not immediately executed on the GPU. It is possible to stack any number of operation into a stream and to have any number of streams. If you want to do one hundred times the operations:<em> memory transfer CPU → GPU, computation, memory transfer GPU → CPU</em>, you can create one hundred streams, each containing the three commands <em>memory transfer CPU → GPU, computing, memory transfer GPU → CPU.</em>


The S_GPU API provides the necessary functions to do these operations. A stream is created with the @ref sg_create_stream function. The  @ref sg_gpu_send_mem, @ref sg_gpu_recv_mem, @ref sg_gpu_send_calc, @ref sg_mem_copy, @ref sg_gpu_send_pimem and @ref sg_gpu_recv_pimem are used to managed streamed memory transfers and calculation. When one of these function is called, the operation is put on the stream, but the function code is executed only when the streams are launched by calling @ref sg_exec_all_streams function.


Some functions (@ref sg_cpu_malloc_pinned, @ref sg_cpu_free_pinned,  @ref sg_gpu_malloc,  @ref sg_gpu_free, @ref sg_send_mem_instantaneously and @ref sg_recv_mem_instantaneously)  does not use stream. It means that their code are executed immediatly, without the stream system. They can be usefull for certain case. For example, when only some byte of data has to be transfered; or @ref sg_cpu_malloc_pinned, when you want to create the pinned buffer needed by the @ref sg_gpu_send_mem and @ref sg_gpu_recv_mem functions.

Two examples, well documented, are provided in the @b example/ directory. The  first one is coded in FORTRAN, the second one in C.

<b>Setting the S_GPU environment</b>

The file GPU.config is the config file for S_GPU. It must be put on the directory where the application using S_GPU is launch and contains several variables

 * USE_SHARED=[0 or 1]
Set if the repartition is shared (1) or static (0)

 * MPI_TASKS_PER_NODE=[n]
Define the number of mpi tasks per node

 * NUM_GPU=[n]
Define the number of GPU device per node

 * GPU_CPUS_AFF_[n]=[x1],[x2],[x3]...
Define the affinity of the processes among the GPUs. It must define as much affinities as GPU devices. eg. if the computer has 2 GPUs and the program run with 8 MPI tasks, the definition could looks like:
GPU_CPUS_AFF_0=0,1,2,3
GPU_CPUS_AFF_1=4,5,6,7

<b> Fortran notes </b>

In fortran, in order to call kernel, you have to call an extra function coded in C. In the fortran example, the binding function is on the binding.c files. It is explained by the fact that the callback function called for the calculation has a C structure in parameter.

!! Be careful, fortan parameter are passed by pointer. !!


 */


#ifdef __cplusplus
extern "C"
{
#endif
  /** 
   * Initialize the S_GPU library. Must be called before any other S_GPU functions
   * 
   * @param GPUshare [out] set to 1 if GPU sharing is enabled
   * @param useGPU [out] set to 1 if GPU should be used
   * @param iproc [in] MPI task number (usefull for debugging)
   * 
   * @return 0 OK, 1 error

   * @b FORTRAN @b NOTES

  In fortran this function has an additional argument ierr at the end of the argument list. ierr is an integer and has the same meaning as the return value of the routine in C. 

   */
  int sg_init ( int *GPUshare, int *useGPU, int iproc, int nproc_node );

  /** 
   * Freed all memory and data structures used by sgpu.
   * 
  
   */

	void sg_end();
  /** 
   * Create one stream
   * 
   * 
   * @return Pointer to the stream created

   * @b FORTRAN @b NOTES

  In fortran this function has an additional argument : stream. stream must be a least a 8 byte data structure (for example a double precision number)

   */
	sg_stream_ptr_t sg_create_stream();

  /** 
   * Free a created stream
   * 
   * @param stream_ptr A pointer to a valid stream
   */
	void sg_delete_stream ( sg_stream_ptr_t stream_ptr);

  /** 
   * Execute ALL the streams created before. After the execution, streams are automaticaly removed
   * 
   */
	void sg_exec_all_streams();
  /** 
   *  Copy  mem_size bytes from memory area src to memory
       area dest.  Function executed on a stream.

   * @param dest pointer to destination area
   * @param src pointer to source area
   * @param mem_size byte to copy
   * @param stream stream to put this function
   * 
   * @return 0 OK, 1 error
   */
	int sg_mem_copy ( void *dest, const void *src, size_t mem_size, sg_stream_ptr_t stream );


  /** 
   * Allocate size byte of pinned memory. The memory is localized on the host.
   * This function is executed immediatly, not on a stream
   * 
   * @param pi_mem [out] pointer to the pinned area allocated
   * @param size byte to allocate
   * 
   * @return 0 OK, 1 error


   * @b FORTRAN @b NOTES

  In fortran this function has an additional argument ierr at the end of the argument list. ierr is an integer and has the same meaning as the return value of the routine in C. 

   */
  int sg_cpu_malloc_pinned(pinned_mem_ptr_t *pi_mem,  size_t size);

  /** 
   * Free a pined memory area
   * 
   * @param pi_mem the pointer
   * 
   * @return 0 OK, 1 error
   */
  int sg_cpu_free_pinned(pinned_mem_ptr_t pi_mem);
  
  /** 
   * Alloc size byte of memory on the GPU
   * 
   * @param gpu_ptr A GPU pointer to the area allocated
   * @param size byte to allocate
   * 
   * @return 0 OK, 1 error

   * @b FORTRAN @b NOTES

  In fortran this function has an additional argument ierr at the end of the argument list. ierr is an integer and has the same meaning as the return value of the routine in C. 


   */
  int sg_gpu_malloc ( GPU_ptr_t *gpu_ptr, size_t size);

  /** 
   * Free a GPU pointer
   * 
   * @param gpu_ptr the pointer
   * 
   * @return 0 OK, 1 error

   * @b FORTRAN @b NOTES

  In fortran this function has an additional argument ierr at the end of the argument list. ierr is an integer and has the same meaning as the return value of the routine in C. 
   */
  int sg_gpu_free ( GPU_ptr_t gpu_ptr);

 
  /** 
   * Send mem_size byte of pinned memory to the GPU. This function must be put on a stream
   * 
   * @param dest pointer to destination area (GPU memory)
   * @param src  pointer to source area (CPU PINNED memory)
   * @param mem_size byte to copy
   * @param stream stream to put this function
   * 
   * @return 0 OK, 1 error

   * @b FORTRAN @b NOTES

  In fortran this function has an additional argument ierr at the end of the argument list. ierr is an integer and has the same meaning as the return value of the routine in C. 

   */
  int sg_gpu_send_pimem ( GPU_ptr_t dest,
			  const pinned_mem_ptr_t src,
			  size_t mem_size,
			  sg_stream_ptr_t stream );
  /** 
   * Send mem_size byte of  memory to the GPU. This memory is automaticaly copied to buff_pi and then copied to the GPU

 This function must be put on a stream
   * 
   * @param dest pointer to destination area (GPU memory)
   * @param src pointer to source area (CPU  memory)
   * @param buff_pi pointer to a PINNED memory buffer 
   * @param mem_size  byte to copy
   * @param stream stream to put this function
   * 
   * @return 0 OK, 1 error


   * @b FORTRAN @b NOTES

  In fortran this function has an additional argument ierr at the end of the argument list. ierr is an integer and has the same meaning as the return value of the routine in C. 


   */
	int sg_gpu_send_mem ( GPU_ptr_t dest,
	                      const void* src,
	                      pinned_mem_ptr_t buff_pi,
	                      size_t mem_size,
	                      sg_stream_ptr_t stream );

  /** 
   * Receiv mem_size byte of  PINNED memory from the GPU. This function must be put on a stream
   * 
   * @param dest  pointer to destination area (CPU PINNED memory)
   * @param src pointer to source area (GPU memory)
   * @param mem_size byte to copy
   * @param stream stream to put this function
   * 

   * @b FORTRAN @b NOTES

  In fortran this function has an additional argument ierr at the end of the argument list. ierr is an integer and has the same meaning as the return value of the routine in C. 

   * @return 
   */
	int sg_gpu_recv_pimem ( pinned_mem_ptr_t dest,
	                        const GPU_ptr_t src,
	                        size_t mem_size,
	                        sg_stream_ptr_t stream );
  /** 
   * Receiv mem_size byte of  memory from the GPU. This memory is automaticaly copied to buff_pi and then copied to dest
   * 
   * @param dest pointer to destination area (CPU memory)
   * @param src pointer to source area (GPU memory)
   * @param buff_pi pointer to a PINNED memory buffer 
   * @param mem_size byte to copy
   * @param stream stream to put this function
   * 
 * @return 

   * @b FORTRAN @b NOTES

  In fortran this function has an additional argument ierr at the end of the argument list. ierr is an integer and has the same meaning as the return value of the routine in C. 

  
   */
	int sg_gpu_recv_mem ( void *dest,
	                      const GPU_ptr_t src,
	                      pinned_mem_ptr_t buff_pi,
	                      size_t mem_size,
	                      sg_stream_ptr_t stream );


  /** 
   * Send a compute kernel to GPU
   * 
   * @param f_call a callback to a function designed to run on the GPU
   * @param param parameter read by the callback
   * @param param_size size in byte of the parameter structure
   * @param stream stream to put this function
   * 

   * @b FORTRAN @b NOTES

This function cannot be called directly in fortran, you have to provide a C function in order to bind s_gpu and your fortran code. See the fortran example for more details

   * @return 
   */
	int sg_gpu_send_calc ( sg_callback_ptr_t f_call,
	                       void *param,
	                       size_t param_size,
	                       sg_stream_ptr_t stream );

  /** 
   * Send size byte of memory from src to the GPU pointer gpu_ptr_dest
   * This function does not use stream, it perfoms the memory transfers instantaneously and return when the transfer is achieved.
   * 
   * @param gpu_ptr_dest destination (GPU memory)
   * @param src  source (CPU memory)
   * @param size size to copy in byte
   * 
   * @return 

  * @b FORTRAN @b NOTES

  In fortran this function has an additional argument ierr at the end of the argument list. ierr is an integer and has the same meaning as the return value of the routine in C. 

   */
  int sg_send_mem_instantaneously(GPU_ptr_t gpu_ptr_dest,
				  void *src,
				  size_t size);

  /** 
   * Receiv size byte of memory from  the GPU pointer gpu_ptr_src to dest.
   * This function does not use stream, it perfoms the memory transfers instantaneously and return when the transfer is achieved.
   * 
   * @param dest destination (CPU memory)
   * @param gpu_ptr_src source (GPU memory)
   * @param size size to copy in byte
   * 
   * @return 

  * @b FORTRAN @b NOTES

  In fortran this function has an additional argument ierr at the end of the argument list. ierr is an integer and has the same meaning as the return value of the routine in C. 

   */
  int sg_recv_mem_instantaneously(void *dest,
				  GPU_ptr_t gpu_ptr_src,
				  size_t size);

	//tracer function
	void add_tracer_callbacks ( tracer_callbacks_t *tcal );

#ifdef __cplusplus
}
#endif

#endif
