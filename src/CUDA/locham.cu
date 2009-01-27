/****u* CUDA/anasyn.cu
**
** 
** AUTHOR
**  Luigi Genovese
**
** SOURCE
*/
  
#include <stdio.h>
#include "commonDef.h"
#include "reduction.h"
#include "GPUparameters.h"

__constant__ parGPU_t par[9];

//global variables for launching different grids

unsigned int numBlocksWTz,linecutsWTz,num_halfwarpsWTz;
unsigned int numBlocksWTy,linecutsWTy,num_halfwarpsWTy;
unsigned int numBlocksWTx,linecutsWTx,num_halfwarpsWTx;
unsigned int numBlocksKz,linecutsKz,num_halfwarpsKz;
unsigned int numBlocksKy,linecutsKy,num_halfwarpsKy;
unsigned int numBlocksKx,linecutsKx,num_halfwarpsKx;
unsigned int numBlocksMFz,linecutsMFz,num_halfwarpsMFz;
unsigned int numBlocksMFy,linecutsMFy,num_halfwarpsMFy;
unsigned int numBlocksMFx,linecutsMFx,num_halfwarpsMFx;


#include "kernels_anasyn.hcu"
#include "kernels_locpot.hcu"
#include "kernels_kinetic.hcu"
#include "kernels_compress.hcu"


template<typename T>
int copygpulochamparameters(int n1,int n2, int n3)
{

  //create the parameters
  parGPU_t parCPU[9];

  //calculate the parameters in constant memory for each of the 1D convolution
  //define the number of threads and blocks according to parameter
  //definitions

  //parameters for the wavelet transformation
  GPUParameters<T>(&parCPU[0],&num_halfwarpsWTz,n3,4*n1*n2,2,LOWFILWT,LUPFILWT,&linecutsWTz,&numBlocksWTz);
  GPUParameters<T>(&parCPU[1],&num_halfwarpsWTy,n2,4*n1*n3,2,LOWFILWT,LUPFILWT,&linecutsWTy,&numBlocksWTy);
  GPUParameters<T>(&parCPU[2],&num_halfwarpsWTx,n1,4*n2*n3,2,LOWFILWT,LUPFILWT,&linecutsWTx,&numBlocksWTx);

  //parameters for the kinetic operator
  GPUParameters<T>(&parCPU[3],&num_halfwarpsKz,2*n3,4*n1*n2,1,LOWFILK,LUPFILK,&linecutsKz,&numBlocksKz);
  GPUParameters<T>(&parCPU[4],&num_halfwarpsKy,2*n2,4*n1*n3,1,LOWFILK,LUPFILK,&linecutsKy,&numBlocksKy);
  GPUParameters<T>(&parCPU[5],&num_halfwarpsKx,2*n1,4*n2*n3,1,LOWFILK,LUPFILK,&linecutsKx,&numBlocksKx);

  //parameters for the magic filters
  GPUParameters<T>(&parCPU[6],&num_halfwarpsMFz,2*n3,4*n1*n2,1,8,7,&linecutsMFz,&numBlocksMFz);
  GPUParameters<T>(&parCPU[7],&num_halfwarpsMFy,2*n2,4*n1*n3,1,8,7,&linecutsMFy,&numBlocksMFy);
  GPUParameters<T>(&parCPU[8],&num_halfwarpsMFx,2*n1,4*n2*n3,1,8,7,&linecutsMFx,&numBlocksMFx);

  //send them to constant memory, once and for all
  if(cudaMemcpyToSymbol(*par,&parCPU, 9*sizeof(parGPU_t)) != 0)
    {
      printf("MemcpyToSymbol error\n");

      return 1;
    }
  return 0;
}

extern "C" 
void creategpuparameters_(int *n1,int *n2, int *n3)
{

  
  if(copygpulochamparameters<double>(*n1+1,*n2+1,*n3+1)!= 0) 
    {
      printf("ERROR: GPU creategpuparameters\n ");
      return;
    } 
  return; 
}


template<typename T>
int completelocalhamiltonian(int n1,int n2, int n3,
			     T h1,T h2,T h3,
			     T *psi,T *out,T *pot,			     
			     T *work,
			     T *work2,
			     T *epot,T *ekinpot)
{

  //create the parameters
  parGPU_t parCPU[9];

  //calculate the number of threads and blocks
  //unsigned int num_lines = min(16,ndat); //hard coded for the moment
  unsigned int numBlocks,linecuts,num_halfwarps;

  //calculate the parameters in constant memory for each of the 1D convolution
  //define the number of threads and blocks according to parameter
  //definitions

  //parameters for the wavelet transformation
  GPUParameters<T>(&parCPU[0],&num_halfwarps,n3,4*n1*n2,2,4,4,&linecuts,&numBlocks);
  dim3  gridWT3(linecuts,  numBlocks, 1);  
  dim3  threadsWT3(HALF_WARP_SIZE, num_halfwarps , 1);
  GPUParameters<T>(&parCPU[1],&num_halfwarps,n2,4*n1*n3,2,4,4,&linecuts,&numBlocks);
  dim3  gridWT2(linecuts,  numBlocks, 1);  
  dim3  threadsWT2(HALF_WARP_SIZE, num_halfwarps , 1);
  GPUParameters<T>(&parCPU[2],&num_halfwarps,n1,4*n2*n3,2,4,4,&linecuts,&numBlocks);
  dim3  gridWT1(linecuts,  numBlocks, 1);  
  dim3  threadsWT1(HALF_WARP_SIZE, num_halfwarps , 1);

  //parameters for the kinetic operator
  GPUParameters<T>(&parCPU[3],&num_halfwarps,2*n3,4*n1*n2,1,14,14,&linecuts,&numBlocks);
  dim3  gridK3(linecuts,  numBlocks, 1);  
  dim3  threadsK3(HALF_WARP_SIZE, num_halfwarps , 1);
  GPUParameters<T>(&parCPU[4],&num_halfwarps,2*n2,4*n1*n3,1,14,14,&linecuts,&numBlocks);
  dim3  gridK2(linecuts,  numBlocks, 1);  
  dim3  threadsK2(HALF_WARP_SIZE, num_halfwarps , 1);
  GPUParameters<T>(&parCPU[5],&num_halfwarps,2*n1,4*n2*n3,1,14,14,&linecuts,&numBlocks);
  dim3  gridK1(linecuts,  numBlocks, 1);  
  dim3  threadsK1(HALF_WARP_SIZE, num_halfwarps , 1);

  //parameters for the magic filters
  GPUParameters<T>(&parCPU[6],&num_halfwarps,2*n3,4*n1*n2,1,8,7,&linecuts,&numBlocks);
  dim3  gridMF3(linecuts,  numBlocks, 1);  
  dim3  threadsMF3(HALF_WARP_SIZE, num_halfwarps , 1);
  GPUParameters<T>(&parCPU[7],&num_halfwarps,2*n2,4*n1*n3,1,8,7,&linecuts,&numBlocks);
  dim3  gridMF2(linecuts,  numBlocks, 1);  
  dim3  threadsMF2(HALF_WARP_SIZE, num_halfwarps , 1);
  GPUParameters<T>(&parCPU[8],&num_halfwarps,2*n1,4*n2*n3,1,8,7,&linecuts,&numBlocks);
  dim3  gridMF1(linecuts,  numBlocks, 1);  
  dim3  threadsMF1(HALF_WARP_SIZE, num_halfwarps , 1);

  //send them to constant memory, once and for all
  if(cudaMemcpyToSymbol(*par,&parCPU, 9*sizeof(parGPU_t)) != 0)
    {
      printf("MemcpyToSymbol error\n");

      return 1;
    }

  //start the hamiltonian calculation
  
  //wavelet synthesis
  waveletsynthesis<T> <<< gridWT3, threadsWT3 >>>(n3,4*n1*n2,psi,work,0);
  cudaThreadSynchronize();
  
  waveletsynthesis<T> <<< gridWT2, threadsWT2 >>>(n2,4*n1*n3,work,work2,1);
  cudaThreadSynchronize();

  waveletsynthesis<T> <<< gridWT1, threadsWT1 >>>(n1,4*n2*n3,work2,psi,2);
  cudaThreadSynchronize();


  //then calculate the potential energy
  magicfilter1d<T> <<< gridMF3, threadsMF3 >>>(2*n3,4*n1*n2,psi,work,6);
  cudaThreadSynchronize();

  magicfilter1d<T> <<< gridMF2, threadsMF2 >>>(2*n2,4*n1*n3,work,work2,7);
  cudaThreadSynchronize();

  //use other work array and save psi for kinetic operator
  magicfilter1d_pot<T> <<< gridMF1, threadsMF1 >>>(2*n1,4*n2*n3,work2,pot,work,8);
  cudaThreadSynchronize();

  //calculate the potential energy
  //reducearrays<T>(2*n1,4*n2*n3,work2,work,epot);
  //cudaThreadSynchronize();

  //reverse MF calculation
  magicfilter1d_t<T> <<< gridMF3, threadsMF3 >>>(2*n3,4*n1*n2,work,work2,6);
  cudaThreadSynchronize();

  magicfilter1d_t<T> <<< gridMF2, threadsMF2 >>>(2*n2,4*n1*n3,work2,work,7);
  cudaThreadSynchronize();

  magicfilter1d_t<T> <<< gridMF1, threadsMF1 >>>(2*n1,4*n2*n3,work,work2,8);
  cudaThreadSynchronize();


  //calculate the kinetic operator and add that to the potential
  //energy
  //use work2 inside the worky array to add the result to the vpsi
  //array

  //here the worky array should be initialised to c*x
  //c_initialize<T> <<< gridK3, threadsK3 >>>(2*n3,4*n1*n2,psi,work2,0.,3);
  //cudaThreadSynchronize();


  //define the scale factor to be applied to the convolution
  T scale=0.5/(h3*h3);
  kinetic1d<T> <<< gridK3, threadsK3 >>>(2*n3,4*n1*n2,scale,psi,work,work2,out,3);
  cudaThreadSynchronize();

  scale=0.5/(h2*h2);
  kinetic1d<T> <<< gridK2, threadsK2 >>>(2*n2,4*n1*n3,scale,work,psi,out,work2,4);
  cudaThreadSynchronize();

  scale=0.5/(h1*h1);
  kinetic1d<T> <<< gridK1, threadsK1 >>>(2*n1,4*n2*n3,scale,psi,work,work2,out,5);
  cudaThreadSynchronize();

  //calculate the reduction of the final array, not correct for ekin
  //reducearrays<T>(2*n1,4*n2*n3,work,out,ekinpot);
  //cudaThreadSynchronize();

  //wavelet analysis
  waveletanalysis<T> <<< gridWT3, threadsWT3 >>>(n3,4*n1*n2,out,work,0);
  cudaThreadSynchronize();

  waveletanalysis<T> <<< gridWT2, threadsWT2 >>>(n2,4*n1*n3,work,work2,1);
  cudaThreadSynchronize();

  waveletanalysis<T> <<< gridWT1, threadsWT1 >>>(n1,4*n2*n3,work2,out,2);
  cudaThreadSynchronize();

  return 0;

}


template<typename T>
int fulllocalhamiltonian(int n1,int n2, int n3,
			 T h1,T h2,T h3,
			 T *psiw,T *hpsiw,T *pot,int *keys,
			 T *psi,T *out, 
			 T *work,
			 T *work2,
			 T *epot,T *ekinpot)
{

  
  dim3  gridWT3(linecutsWTz,  numBlocksWTz, 1);  
  dim3  threadsWT3(HALF_WARP_SIZE, num_halfwarpsWTz , 1);
  dim3  gridWT2(linecutsWTy,  numBlocksWTy, 1);  
  dim3  threadsWT2(HALF_WARP_SIZE, num_halfwarpsWTy , 1);
  dim3  gridWT1(linecutsWTx,  numBlocksWTx, 1);  
  dim3  threadsWT1(HALF_WARP_SIZE, num_halfwarpsWTx , 1);

  dim3  gridK3(linecutsKz,  numBlocksKz, 1);  
  dim3  threadsK3(HALF_WARP_SIZE, num_halfwarpsKz , 1);
  dim3  gridK2(linecutsKy,  numBlocksKy, 1);  
  dim3  threadsK2(HALF_WARP_SIZE, num_halfwarpsKy , 1);
  dim3  gridK1(linecutsKx,  numBlocksKx, 1);  
  dim3  threadsK1(HALF_WARP_SIZE, num_halfwarpsKx , 1);

  dim3  gridMF3(linecutsMFz,  numBlocksMFz, 1);  
  dim3  threadsMF3(HALF_WARP_SIZE, num_halfwarpsMFz , 1);
  dim3  gridMF2(linecutsMFy,  numBlocksMFy, 1);  
  dim3  threadsMF2(HALF_WARP_SIZE, num_halfwarpsMFy , 1);
  dim3  gridMF1(linecutsMFx,  numBlocksMFx, 1);  
  dim3  threadsMF1(HALF_WARP_SIZE, num_halfwarpsMFx , 1);

  dim3  gridC(nblocksC, 1, 1);  
  dim3  threadsC(ELEMS_BLOCK, nseg_blockC , 1);

  //uncompression
  //set the value of the psig array to zero
  cudaMemset((void*) psi,0,8*n1*n2*n3*sizeof(T));

  uncompresscoarsefine<T> <<< gridC, threadsC >>>(n1,n2,n3,psiw,psi,keys);
  cudaThreadSynchronize();


  //start the hamiltonian calculation
  
  //wavelet synthesis
  waveletsynthesis<T> <<< gridWT3, threadsWT3 >>>(n3,4*n1*n2,psi,work,0);
  cudaThreadSynchronize();
  
  waveletsynthesis<T> <<< gridWT2, threadsWT2 >>>(n2,4*n1*n3,work,work2,1);
  cudaThreadSynchronize();

  waveletsynthesis<T> <<< gridWT1, threadsWT1 >>>(n1,4*n2*n3,work2,psi,2);
  cudaThreadSynchronize();

  //then calculate the potential energy
  magicfilter1d<T> <<< gridMF3, threadsMF3 >>>(2*n3,4*n1*n2,psi,work,6);
  cudaThreadSynchronize();

  magicfilter1d<T> <<< gridMF2, threadsMF2 >>>(2*n2,4*n1*n3,work,work2,7);
  cudaThreadSynchronize();

  //use other work array and save psi for kinetic operator
  magicfilter1d_pot<T> <<< gridMF1, threadsMF1 >>>(2*n1,4*n2*n3,work2,pot,work,8);
  cudaThreadSynchronize();

  //calculate the potential energy
  //reducearrays<T>(2*n1,4*n2*n3,work2,work,epot);
  //cudaThreadSynchronize();

  //reverse MF calculation
  magicfilter1d_t<T> <<< gridMF3, threadsMF3 >>>(2*n3,4*n1*n2,work,work2,6);
  cudaThreadSynchronize();

  magicfilter1d_t<T> <<< gridMF2, threadsMF2 >>>(2*n2,4*n1*n3,work2,work,7);
  cudaThreadSynchronize();

  magicfilter1d_t<T> <<< gridMF1, threadsMF1 >>>(2*n1,4*n2*n3,work,work2,8);
  cudaThreadSynchronize();


  //calculate the kinetic operator and add that to the potential
  //energy
  //use work2 inside the worky array to add the result to the vpsi
  //array


  //here the worky array should be initialised to c*x
  //c_initialize<T> <<< gridK3, threadsK3 >>>(2*n3,4*n1*n2,psi,work2,0.,3);
  //cudaThreadSynchronize();


  //define the scale factor to be applied to the convolution
  T scale=0.5/(h3*h3);
  kinetic1d<T> <<< gridK3, threadsK3 >>>(2*n3,4*n1*n2,scale,psi,work,work2,out,3);
  cudaThreadSynchronize();

  scale=0.5/(h2*h2);
  kinetic1d<T> <<< gridK2, threadsK2 >>>(2*n2,4*n1*n3,scale,work,psi,out,work2,4);
  cudaThreadSynchronize();

  scale=0.5/(h1*h1);
  kinetic1d<T> <<< gridK1, threadsK1 >>>(2*n1,4*n2*n3,scale,psi,work,work2,out,5);
  cudaThreadSynchronize();


  //calculate the reduction of the final array, not correct for ekin
  //reducearrays<T>(2*n1,4*n2*n3,work,out,ekinpot);
  //cudaThreadSynchronize();

  //wavelet analysis
  waveletanalysis<T> <<< gridWT3, threadsWT3 >>>(n3,4*n1*n2,out,work,0);
  cudaThreadSynchronize();

  waveletanalysis<T> <<< gridWT2, threadsWT2 >>>(n2,4*n1*n3,work,work2,1);
  cudaThreadSynchronize();

  waveletanalysis<T> <<< gridWT1, threadsWT1 >>>(n1,4*n2*n3,work2,out,2);
  cudaThreadSynchronize();

  //recompress
  compresscoarsefine<T> <<< gridC, threadsC >>>(n1,n2,n3,out,hpsiw,keys);
  cudaThreadSynchronize();


  return 0;

}


extern "C" 
void gpulocham_(int *n1,int *n2, int *n3,
		double *h1,double *h2,double *h3,
		double **psi,double **out,double **pot,			     
		double **work,
		double **work2,
		double *epot,double *ekinpot)
{

  
  if(completelocalhamiltonian<double>(*n1+1,*n2+1,*n3+1,
				      *h1,*h2,*h3,
				      *psi,*out,*pot,*work,*work2,
				      epot,ekinpot)!= 0) 
    {
      printf("ERROR: GPU localhamiltonian\n ");
      return;
    } 
  return; 
}


/****/


extern "C" 
void gpufulllocham_(int *n1,int *n2, int *n3,
		    double *h1,double *h2,double *h3,
		    double **psiw,double **hpsiw,double **pot,int **keys, 
		    double **psi,double **out,
		    double **work,
		    double **work2,
		    double *epot,double *ekinpot)
{

  
  if(fulllocalhamiltonian<double>(*n1+1,*n2+1,*n3+1,
				  *h1,*h2,*h3,
				  *psiw, *hpsiw, *pot, *keys,*psi,*out,*work,*work2,
				  epot,ekinpot)!= 0) 
    {
      printf("ERROR: GPU fulllocalhamiltonian\n ");
      return;
    } 
  return; 
}

extern "C" 
void ana1d_(int *n,int *ndat,
	    double **psi,double **out) 

{

  
  if(wavana1d<double>(*n+1,*ndat,
		      *psi,
		      *out) != 0) 
    {
      printf("ERROR: GPU waveletanalysis\n ");
      return;
    } 
  return; 
}

extern "C" 
void syn1d_(int *n,int *ndat,
	    double **psi,double **out) 

{

  
  if(wavsyn1d<double>(*n+1,*ndat,
		      *psi,
		      *out) != 0) 
    {
      printf("ERROR: GPU waveletanalysis\n ");
      return;
    } 
  return; 
}


extern "C" 
void localpotential_(int *n1,
		     int *n2,
		     int *n3,
		     float **psi,
		     float **work,
		     float **pot,
		     float *epot) 

{  
  if(magicfilterpot<float>(*n1+1,*n2+1,*n3+1,
		    *psi, 
		    *work,
		    *pot,
		    epot) != 0) 
    {
      printf("ERROR: GPU magicfilterpot\n ");
      return;
    } 
  return; 
}

extern "C" 
void localpotentiald_(int *n1,
		     int *n2,
		     int *n3,
		     double **psi,
		     double **work,
		     double **pot,
		     double *epot) 

{  
  if(magicfilterpot<double>(*n1+1,*n2+1,*n3+1,
		    *psi, 
		    *work,
		    *pot,
		    epot) != 0) 
    {
      printf("ERROR: GPU magicfilterpot\n ");
      return;
    } 
  return; 
}

extern "C" 
void magicfilter1d_(int *n,int *ndat,
		    double **psi,
		    double **out) 

{  
  if(mf1d<double>(*ndat,*n+1,
		  *psi, 
		  *out) != 0) 
    {
      printf("ERROR: GPU magicfilterpot\n ");
      return;
    } 
  return; 
}


extern "C"
void uncompressgpu_(int* n1, int* n2,int *n3 ,
		    double **psi, double **psig, int **keys)
{
  if(uncompressgpu<double>(*n1+1,*n2+1,*n3+1,*psi,*psig,*keys) != 0)
    {
      printf("MemcpyToSymbol error\n");

      return;
    }

}

extern "C"
void compressgpu_(int* n1, int* n2,int *n3 ,
		    double **psig, double **psi, int **keys)
{
  if(compressgpu<double>(*n1+1,*n2+1,*n3+1,*psig,*psi,*keys) != 0)
    {
      printf("MemcpyToSymbol error\n");

      return;
    }

}

extern "C"
void copynseginfo_(int *nblocks, int *nseg_block)
{
  nseg_blockC = *nseg_block;  
  nblocksC = *nblocks;

  //printf("nsegC %i, nblC %i,\n",nseg_blockC,nblocksC);
}


extern "C"
void readkeysinfo_(int *elems_block, int *nblocks_max)
{
  *elems_block=ELEMS_BLOCK;
  *nblocks_max=NBLOCKS_MAX;
}

