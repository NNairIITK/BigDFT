/****u* CUDA/conv_shared_multi_kernel.cu
**
** 
** AUTHOR
**  Matthieu Ospici
**
** SOURCE
*/

#ifndef _CONV_SHARED_KERNEL_MULTI_
#define _CONV_SHARED_KERNEL_MULTI_
#include <math.h>

//#include "cuda_runtime.h"
#include "structUtil.h"
///////////////////////////////////////////////////////////////////////////////
//f*g convolution



///////////////////////////////////////////////////////////////////////////////

//24-bit multiplication is faster on G80,
//but we must be sure to multiply integers
//only within [-8M, 8M - 1] range
#define IMUL(a, b) __mul24(a, b)



__device__ inline int threadToTab(int thid, const int* tab)
{
  return tab[thid/16];
}

__device__ inline unsigned int offsetToTab(int thid, const unsigned int* tab)
{
  return tab[thid/16];
}



__constant__ param_t param;


__global__ void conv_shared_multi(unsigned int n1,unsigned int n2,float *t_out,float *t_in,int nf)
{
  __shared__ float shared_temp[SIZE_SHARED_TOTAL];




  const unsigned int thid = threadIdx.x; //ID of current thread
  const unsigned int bidx = blockIdx.x;
  const unsigned int SIZE_SHARED_1 = param.SIZE_SHARED_1; 
  const unsigned int lineNumber = thid % 16;  


  int lineNumberFetch = lineNumber;

  if(bidx == gridDim.x - 1)
    {
      lineNumberFetch = (int)lineNumber - (int)param.lineLastBlock;

      if(lineNumberFetch <= 0)
	lineNumberFetch += (int)param.lineLastBlock;
      else
	lineNumberFetch = 0;
    }


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

#endif 
/****/
