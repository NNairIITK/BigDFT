/*!
 @author
    Copyright (C) 2010-2011 BigDFT group (LG)
    This file is distributed under the terms of the
    GNU General Public License, see ~/COPYING file
    or http://www.gnu.org/copyleft/gpl.txt .
    For the list of contributors, see ~/AUTHORS 

*/
  
#include <stdio.h>
#include "locpot.h" //function declaration
#include "commonDef.h"
#include "structDef_locpot.h"

//__constant__ parGPU_t par[3]; //must be defined here because it is used by both by this file and kernels_locpot.cu

#include "kernels_locpot.hcu"






template<typename T>
int magicfilterpot(int n1,int n2, int n3,
		   T *psi,
		   T *work,
		   T *pot,
		   T *epot)
{

  //create the parameters
  parGPU_t parCPU[3];

  //calculate the number of threads and blocks
  //unsigned int num_lines = min(16,ndat); //hard coded for the moment
  unsigned int numBlocks,linecuts,num_halfwarps;

  //calculate the parameters in constant memory for each of the 1D convolution
  //define the number of threads and blocks according to parameter definitions
  GPUParameters<T>(&parCPU[2],&num_halfwarps,n3,n1*n2,1,LOWFILMF,LUPFILMF,&linecuts,&numBlocks);
  dim3  grid3(linecuts,  numBlocks, 1);  
  dim3  threads3(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n3,n1*n2);

  GPUParameters<T>(&parCPU[1],&num_halfwarps,n2,n1*n3,1,LOWFILMF,LUPFILMF,&linecuts,&numBlocks);
  dim3  grid2(linecuts,  numBlocks, 1);  
  dim3  threads2(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n2,n1*n3);

  GPUParameters<T>(&parCPU[0],&num_halfwarps,n1,n2*n3,1,LOWFILMF,LUPFILMF,&linecuts,&numBlocks);
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

  //magicfilter1d <<< grid1, threads1 >>>(n1,n2*n3,psi,work,0);
  //cudaThreadSynchronize();


  //direct MF calculation

  //bind the texture reference to the input array
  //cudaBindTexture(NULL,psi_tex,GPU_idata,n*ndat*sizeof(float));

  //launch the kernel grid
  //conv1d_stride <<< grid1, threads1 >>>(n,ndat, GPU_odata);
  magicfilter1d<T> <<< grid3, threads3 >>>(n3,n1*n2,psi,work,2);
  //unbind the texture
  //cudaUnbindTexture(psi_tex);

  cudaThreadSynchronize();

  magicfilter1d<T> <<< grid2, threads2 >>>(n2,n1*n3,work,psi,1);
  cudaThreadSynchronize();

  magicfilter1d_pot<T> <<< grid1, threads1 >>>(n1,n2*n3,psi,pot,work,0);
  cudaThreadSynchronize();

  //here one should combine psi and work to calculate the potential
  //energy
  //reducearrays<T>(n1,n2*n3,psi,work,epot);
  //cudaThreadSynchronize(); //can be removed since the arrays are not overwritten

  //reverse MF calculation
  magicfilter1d_t<T> <<< grid3, threads3 >>>(n3,n1*n2,work,psi,2);
  cudaThreadSynchronize();

  magicfilter1d_t<T> <<< grid2, threads2 >>>(n2,n1*n3,psi,work,1);
  cudaThreadSynchronize();

  magicfilter1d_t<T> <<< grid1, threads1 >>>(n1,n2*n3,work,psi,0);
  cudaThreadSynchronize();

  return 0;

}

template<typename T>
int mf1d(int ndat, int n3,
	 T *psi,
	 T *out)
{

  //create the parameters
  parGPU_t parCPU[1];

  //calculate the number of threads and blocks
  //unsigned int num_lines = min(16,ndat); //hard coded for the moment
  unsigned int numBlocks,linecuts,num_halfwarps;

  //calculate the parameters in constant memory for each of the 1D convolution
  //define the number of threads and blocks according to parameter definitions
  GPUParameters<T>(&parCPU[0],&num_halfwarps,n3,ndat,1,LOWFILMF,LUPFILMF,&linecuts,&numBlocks);
  dim3  grid3(linecuts,  numBlocks, 1);  
  dim3  threads3(HALF_WARP_SIZE, num_halfwarps , 1);

  if(cudaMemcpyToSymbol(*par,&parCPU, sizeof(parGPU_t)) != 0)
    {
      printf("MemcpyToSymbol error\n");

      return 1;
    }

  magicfilter1d<T> <<< grid3, threads3 >>>(n3,ndat,psi,out,0);
  cudaThreadSynchronize();

  return 0;

}
/****/
