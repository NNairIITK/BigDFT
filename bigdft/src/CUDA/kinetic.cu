/*!
 @author
    Copyright (C) 2010-2011 BigDFT group (LG)
    This file is distributed under the terms of the
    GNU General Public License, see ~/COPYING file
    or http://www.gnu.org/copyleft/gpl.txt .
    For the list of contributors, see ~/AUTHORS 

*/
  
#include <stdio.h>

#include "kinetic.h"

#include "commonDef.h"
#include "structDef_kinetic.h"
#include "GPUparameters.h"


//__constant__ parGPU_t par[3];


#include "kernels_kinetic.hcu"



template<typename T>
int kineticfilter(int n1,int n2, int n3,
		  T h1,T h2,T h3,T c,
		  T *x,
		  T *workx,
		  T *y,
		  T *worky,
		  T *ekin)
{

  //create the parameters
  parGPU_t parCPU[3];

  //calculate the number of threads and blocks
  unsigned int numBlocks,linecuts,num_halfwarps;

  //calculate the parameters in constant memory for each of the 1D convolution
  //define the number of threads and blocks according to parameter definitions
  GPUParameters<T>(&parCPU[2],&num_halfwarps,n3,n1*n2,1,LOWFILK,LUPFILK,&linecuts,&numBlocks);
  dim3  grid3(linecuts,  numBlocks, 1);  
  dim3  threads3(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n3,n1*n2);

  GPUParameters<T>(&parCPU[1],&num_halfwarps,n2,n1*n3,1,LOWFILK,LUPFILK,&linecuts,&numBlocks);
  dim3  grid2(linecuts,  numBlocks, 1);  
  dim3  threads2(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n2,n1*n3);

  GPUParameters<T>(&parCPU[0],&num_halfwarps,n1,n2*n3,1,LOWFILK,LUPFILK,&linecuts,&numBlocks);
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


  //here the worky array should be initialised to c*x
  c_initialize<T> <<< grid3, threads3 >>>(n3,n1*n2,x,worky,c,2);
  cudaThreadSynchronize();


  //define the scale factor to be applied to the convolution
  T scale=0.5/(h3*h3);

  kinetic1d<T> <<< grid3, threads3 >>>(n3,n1*n2,scale,x,workx,worky,y,2);
  cudaThreadSynchronize();

  scale=0.5/(h2*h2);
  kinetic1d<T> <<< grid2, threads2 >>>(n2,n1*n3,scale,workx,x,y,worky,1);
  cudaThreadSynchronize();

  scale=0.5/(h1*h1);
  kinetic1d<T> <<< grid1, threads1 >>>(n1,n2*n3,scale,x,workx,worky,y,0);
  cudaThreadSynchronize();

  //then calculate the kinetic energy
  // reducearrays<T>(n1,n2*n3,x,y,ekin);
  //  cudaThreadSynchronize();

  return 0;

}

template<typename T>
int k1d(int ndat, int n,
	T h,T c,
	T *x,
	T *workx,
	T *y,
	T *worky,
	T *ekin)
{

  //create the parameters
  parGPU_t parCPU[1];

  //calculate the number of threads and blocks
  unsigned int numBlocks,linecuts,num_halfwarps;

  //calculate the parameters in constant memory for each of the 1D convolution
  //define the number of threads and blocks according to parameter definitions
  GPUParameters<T>(&parCPU[0],&num_halfwarps,n,ndat,1,LOWFILK,LUPFILK,&linecuts,&numBlocks);
  dim3  grid3(linecuts,  numBlocks, 1);  
  dim3  threads3(HALF_WARP_SIZE, num_halfwarps , 1);

  //send them to constant memory, once and for all
  if(cudaMemcpyToSymbol(*par,&parCPU, sizeof(parGPU_t)) != 0)
    {
      printf("MemcpyToSymbol error\n");

      return 1;
    }


  //define the scale factor to be applied to the convolution
  T scale=0.5/(h*h);

  //here the worky array should be initialised to c*x
  c_initialize<T> <<< grid3, threads3 >>>(n,ndat,x,worky,c,0);
  cudaThreadSynchronize();

  kinetic1d<T> <<< grid3, threads3 >>>(n,ndat,scale,x,workx,worky,y,0);
  cudaThreadSynchronize();

  //then calculate the kinetic energy
  //  reducearrays<T>(n,ndat,x,y,ekin);
  //  cudaThreadSynchronize();
  return 0;

}
/****/
