#include <stdio.h>
#include <cutil.h>
//#include <multithreading.h>
#include <pthread.h>
#include <semaphore.h>
#include <sched.h>

#include <time.h>
#include <sys/time.h>

#include "convolution.h"

//#include "conv_pot_shared_kernel.cu"
//#include "deffct.h" //convseria

#define NUM_MULTIPROCESSORS 30 // valid for GTX280
#define BLOCK_SIZE 256 //must be at least 64


__device__ void copy_rectangle(unsigned int ndat, //total number of lines
			       unsigned int num_elem, //copied by each thread
			       unsigned int num_lines, //lines treated 
			       unsigned int dim_rectangle
			       const unsigned int tid, //thread id (1-dim)
			       float *psi_in, //input array (global)
			       float *psi_sh) //copy array (shared)
{
  int tidk,nx_sh,ndatx_sh

  for(int k=0; k < num_elem-1; ++k)
    {
      tidk=num_elem*tid+k

      //n axis position of the given thread in shared memory
      int nx_sh = tidk/num_lines;
      //ndat axis position of the given thread in shared memory
      // to be controlled if it is faster than tidk % num_lines
      int ndatx_sh = tidk - nx_sh*num_lines;

      //each thread copies the corresponding point in shared memory

      //bank conflicts for psi_sh: for each half-warp there is 
      //a linear addressing with stride num_lines (which MUST be odd
      //or equal to the number of banks)

      //coalesced accesses for psi_in: the access are completely
      //uncoalesced, must pass through texture mapping
      if(tidk <= dim_rectangle) //the warp-consistency should be checked separately
	{
	  psi_sh[tidk]=psi_in[nx_sh*ndat+ndatx_sh];
	}
  
    }

} 

__global__ void conv_lg(unsigned int n1,unsigned int n2,float *t_out,float *t_in,int nf)
{
  __shared__ float psi_sh[dim_shared];

  const unsigned int tid = threadIdx.x; //ID of current thread
  const unsigned int bid = blockIdx.x;

  //If the total number of element to be convoluted
  //is lower that NUM_MULTIPROCESSORS*BLOCK_SIZE the GPU strategy is
  //not anymore convenient

  //starting point of n axis in the given block (global memory)
  int n0_gb = tba;

  //starting point of ndat axis in the given block (global memory)
  int ndat0_gb = tba;

  int istart_gb = n0_gb*ndat + ndat0_gb

  //here we can start the for loop
  //we can copy NUM_THREADS*num_elem elements
  //so num_lines is NUM_THREADS/num_columns
  copy_rectangle(ndat,num_elem,num_lines,numthreads*num_elem,tid,
		 psi_in + istart_gb,psi_sh + num_lines*lowfil);

  //continue copying up to fill what it is lacking

  //then come back in the 

  //after the main copy we should copy the data for periodicity
  //the copy should be taken at the edge or on the other side
  //depending on the value of the starting point

  //number of other elements to be copied (right part)
  int nright = lupfil*num_lines ;
  
  //number of elements which still have to be copied in the right part
  int nrightstill = - numthreads*num_elem (num_lines+lupfil)*num_columns;

  //number of other elements to be copied (left part)
  int nleft = lowfil*num_lines;

  //end memory transfer
   __syncthreads();
 
  //then we perform convolutions, switching the values of shared
  //memory array


  //and then we copy back to global memory, in a coalesced way
  //beware of the problem of alignment for the array
  

   psi_out[iout_gb+tid]=psi_sh[iout_sh+tid];




  //copy into memory

  const unsigned int  step = (gridDim.y == 1) ? 0 : blockIdx.y/(gridDim.y - 1);
  const unsigned int TO_COPY =  param.fetchTabs.tabSeq[step][thid/16];    
  const int offset = param.fetchTabs.currPosTab[step][thid/16]; // collums offset (per thread)         
  const unsigned int TO_COPY_CALC =  param.calcTabs.tabSeq[step][thid/16];       
  const unsigned int offset_calc = param.calcTabs.currPosTab[step][thid/16]; // calc collums offset (per thread)


  int baseOffset;
      
  if(gridDim.y == 1)
    baseOffset = 0; //only one block
  else
    {
      baseOffset = param.sizeLineFirstBlocks * blockIdx.y;
    }
  
  
  
  for(int i=0,mod = baseOffset + offset + param.lowfil; i < TO_COPY; ++i)
    {
      if(mod >= param.SIZE_SHARED_2)
	{ 
	  mod = mod - param.SIZE_SHARED_2;
	}
      else if(mod < 0)
	{
	  mod = (param.SIZE_SHARED_2 + mod);
	}
      
      shared_temp[( i + offset)*16 + lineNumber] = t_in[(lineNumberFetch + SIZE_SHARED_1 *(bidx)) + mod*n1];
      
      
      ++mod;
    }
  
  
  __syncthreads();
  
  
  for(int i=0;i < TO_COPY_CALC;++i)
    {
      register float tmp = 0;
      
      #pragma unroll 20
      for(unsigned int j=0 ;j < nf;++j)
	{
	  
	  
	  tmp += shared_temp[lineNumber + (i + offset_calc + j)*16] * param.f_fct[j];

	}   
      t_out[(lineNumberFetch + SIZE_SHARED_1 *(bidx))*n2 +(i+offset_calc + baseOffset)] = tmp;
    }

}

