/*! CUDA/anasyn.cu
 @author
    Copyright (C) 2010-2011 BigDFT group (LG)
    This file is distributed under the terms of the
    GNU General Public License, see ~/COPYING file
    or http://www.gnu.org/copyleft/gpl.txt .
    For the list of contributors, see ~/AUTHORS 
*/
  
#include <stdio.h>
#include "commonDef.h"
#include "structDef_anasyn.h"
#include "GPUparameters.h"

//__constant__ parGPU_t par[3];


#include "kernels_anasyn.hcu"

template<typename T>
int wavana1d(int n,int ndat,
	   T *psi, T *out)
{

  //create the parameters
  parGPU_t parCPU[1];

  //calculate the number of threads and blocks
  //unsigned int num_lines = min(16,ndat); //hard coded for the moment
  unsigned int numBlocks,linecuts,num_halfwarps;

  //calculate the parameters in constant memory for each of the 1D convolution
  //define the number of threads and blocks according to parameter definitions
  GPUParameters<T>(&parCPU[0],&num_halfwarps,n,ndat,2,LOWFILWT,LUPFILWT,&linecuts,&numBlocks);
  dim3  grid3(linecuts,  numBlocks, 1);  
  dim3  threads3(HALF_WARP_SIZE, num_halfwarps , 1);

  //send them to constant memory, once and for all
  if(cudaMemcpyToSymbol(*par,&parCPU, sizeof(parGPU_t)) != 0)
    {
      printf("MemcpyToSymbol error\n");

      return 1;
    }

  waveletanalysis<T> <<< grid3, threads3 >>>(n,ndat,psi,out,0);
  cudaThreadSynchronize();

  return 0;

}

template<typename T>
int wavsyn1d(int n,int ndat,
	   T *psi, T *out)
{

  //create the parameters
  parGPU_t parCPU[1];

  //calculate the number of threads and blocks
  //unsigned int num_lines = min(16,ndat); //hard coded for the moment
  unsigned int numBlocks,linecuts,num_halfwarps;

  //calculate the parameters in constant memory for each of the 1D convolution
  //define the number of threads and blocks according to parameter definitions
  GPUParameters<T>(&parCPU[0],&num_halfwarps,n,ndat,2,LOWFILWT,LUPFILWT,&linecuts,&numBlocks);
  dim3  grid3(linecuts,  numBlocks, 1);  
  dim3  threads3(HALF_WARP_SIZE, num_halfwarps , 1);

  //send them to constant memory, once and for all
  if(cudaMemcpyToSymbol(*par,&parCPU, sizeof(parGPU_t)) != 0)
    {
      printf("MemcpyToSymbol error\n");

      return 1;
    }

  waveletsynthesis<T> <<< grid3, threads3 >>>(n,ndat,psi,out,0);
  cudaThreadSynchronize();

  return 0;

}



template<typename T>
int wavana(int n1,int n2, int n3,
	   T *psi, T *out)
{

  //create the parameters
  parGPU_t parCPU[3];

  //calculate the number of threads and blocks
  //unsigned int num_lines = min(16,ndat); //hard coded for the moment
  unsigned int numBlocks,linecuts,num_halfwarps;

  //calculate the parameters in constant memory for each of the 1D convolution
  //define the number of threads and blocks according to parameter definitions
  GPUParameters<T>(&parCPU[2],&num_halfwarps,n3,n1*n2,2,LOWFILWT,LUPFILWT,&linecuts,&numBlocks);
  dim3  grid3(linecuts,  numBlocks, 1);  
  dim3  threads3(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n3,n1*n2);

  GPUParameters<T>(&parCPU[1],&num_halfwarps,n2,n1*n3,2,LOWFILWT,LUPFILWT,&linecuts,&numBlocks);
  dim3  grid2(linecuts,  numBlocks, 1);  
  dim3  threads2(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n2,n1*n3);

  GPUParameters<T>(&parCPU[0],&num_halfwarps,n1,n2*n3,2,LOWFILWT,LUPFILWT,&linecuts,&numBlocks);
  dim3  grid1(linecuts,  numBlocks, 1);  
  dim3  threads1(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n1,n3*n2);

  //send them to constant memory, once and for all
  if(cudaMemcpyToSymbol(*par,&parCPU, 3*sizeof(parGPU_t)) != 0)
    {
      printf("MemcpyToSymbol error\n");

      return 1;
    }

  waveletanalysis<T> <<< grid3, threads3 >>>(n3,n1*n2,psi,out,2);
  cudaThreadSynchronize();

  waveletanalysis<T> <<< grid2, threads2 >>>(n2,n1*n3,out,psi,1);
  cudaThreadSynchronize();

  waveletanalysis<T> <<< grid1, threads1 >>>(n1,n2*n3,psi,out,0);
  cudaThreadSynchronize();


  return 0;

}

//MUST CONTINUE WITH #D ANALYSIS SYNTHESIS VERSION

/****/

