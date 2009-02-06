/****u* CUDA/anasyn.cu
**
** 
** AUTHOR
**  Luigi Genovese
**
** SOURCE
*/
  
#include <stdio.h>
#include <cublas.h>
#include "commonDef.h"
#include "reduction.h"
#include "GPUparameters.h"


//#include "locpot.h"  //for mf1d and magicfilterpot fcts
//#include "anasyn.h" //wavana1d wavsyn1d
//#include "compress.h" // uncompressgpu compressgpu
__constant__ parGPU_t par[9];

//global variables for launching different grids
static unsigned int numBlocksWTz,linecutsWTz,num_halfwarpsWTz;
static unsigned int numBlocksWTy,linecutsWTy,num_halfwarpsWTy;
static unsigned int numBlocksWTx,linecutsWTx,num_halfwarpsWTx;
static unsigned int numBlocksKz,linecutsKz,num_halfwarpsKz;
static unsigned int numBlocksKy,linecutsKy,num_halfwarpsKy;
static unsigned int numBlocksKx,linecutsKx,num_halfwarpsKx;
static unsigned int numBlocksMFz,linecutsMFz,num_halfwarpsMFz;
static unsigned int numBlocksMFy,linecutsMFy,num_halfwarpsMFy;
static unsigned int numBlocksMFx,linecutsMFx,num_halfwarpsMFx;
static unsigned int nseg_blockC,nblocksC;



#include "anasyn.cu"
#include "locpot.cu"
#include "kinetic.cu"
#include "compress.cu"



//template specialization for simple or double
template<typename T>
T do_blasDot(int size, T *tab1, T *tab2)
{
  return -1; //if this code is expended is a bug because T is neither double nor simple
}

template<>
float do_blasDot<float>(int size, float *tab1, float *tab2)
{
    return cublasSdot(size, 
		     tab1, 1,
		     tab2, 1);
}

template<>
double do_blasDot(int size, double *tab1, double *tab2)
{

  double ret = cublasDdot(size, 
			  tab1, 1,
			  tab2, 1);


 
  /*  if(cublasGetError() == CUBLAS_STATUS_NOT_INITIALIZED)
    {
      printf("CUBLAS_STATUS_NOT_INITIALIZED\n");
    }
  else if(cublasGetError() == CUBLAS_STATUS_ALLOC_FAILED)
      {
	printf("CUBLAS_STATUS_ALLOC_FAILED\n");
      }
  else if(cublasGetError() == CUBLAS_STATUS_ARCH_MISMATCH)
	{
	  printf("CUBLAS_STATUS_ARCH_MISMATCH\n");
	}
  else if(cublasGetError() == CUBLAS_STATUS_EXECUTION_FAILED)
	  {
	    printf("CUBLAS_STATUS_EXECUTION_FAILED\n");
	    }*/

   
  return ret;
}
// end template specialisation



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
  GPUParameters<T>(&parCPU[6],&num_halfwarpsMFz,2*n3,4*n1*n2,1,LOWFILMF,LUPFILMF,&linecutsMFz,&numBlocksMFz);
  GPUParameters<T>(&parCPU[7],&num_halfwarpsMFy,2*n2,4*n1*n3,1,LOWFILMF,LUPFILMF,&linecutsMFy,&numBlocksMFy);
  GPUParameters<T>(&parCPU[8],&num_halfwarpsMFx,2*n1,4*n2*n3,1,LOWFILMF,LUPFILMF,&linecutsMFx,&numBlocksMFx);

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
int waveletstosf(int n1,int n2, int n3,
		 dim3 gridWT1, dim3 gridWT2, dim3 gridWT3,
		 dim3 threadsWT1, dim3 threadsWT2, dim3 threadsWT3,
		 T *in,T *out)
{

  //wavelet synthesis
  waveletsynthesis<T> <<< gridWT3, threadsWT3 >>>(n3,4*n1*n2,in,out,0);
  cudaThreadSynchronize();
  
  waveletsynthesis<T> <<< gridWT2, threadsWT2 >>>(n2,4*n1*n3,out,in,1);
  cudaThreadSynchronize();

  waveletsynthesis<T> <<< gridWT1, threadsWT1 >>>(n1,4*n2*n3,in,out,2);
  cudaThreadSynchronize();

  return 0;
}

template<typename T>
int sftowavelets(int n1,int n2, int n3,
		 dim3 gridWT1, dim3 gridWT2, dim3 gridWT3,
		 dim3 threadsWT1, dim3 threadsWT2, dim3 threadsWT3,
		 T *in,T *out)
{

  //wavelet analysis
  waveletanalysis<T> <<< gridWT3, threadsWT3 >>>(n3,4*n1*n2,in,out,0);
  cudaThreadSynchronize();

  waveletanalysis<T> <<< gridWT2, threadsWT2 >>>(n2,4*n1*n3,out,in,1);
  cudaThreadSynchronize();

  waveletanalysis<T> <<< gridWT1, threadsWT1 >>>(n1,4*n2*n3,in,out,2);
  cudaThreadSynchronize();

  return 0;
}



template<typename T>
int potentialapplication(int n1,int n2, int n3,
			 dim3 gridMF1, dim3 gridMF2, dim3 gridMF3,
			 dim3 threadsMF1, dim3 threadsMF2, dim3 threadsMF3,
			 T *psi,T *out,T *pot,T *epot,
			 T *work)
{

  //calculate the MF transformation
  magicfilter1d<T> <<< gridMF3, threadsMF3 >>>(2*n3,4*n1*n2,psi,work,6);
  cudaThreadSynchronize();

  magicfilter1d<T> <<< gridMF2, threadsMF2 >>>(2*n2,4*n1*n3,work,out,7);
  cudaThreadSynchronize();

  magicfilter1d_pot<T> <<< gridMF1, threadsMF1 >>>(2*n1,4*n2*n3,out,pot,work,8);
  cudaThreadSynchronize();

  //reverse MF calculation
  magicfilter1d_t<T> <<< gridMF3, threadsMF3 >>>(2*n3,4*n1*n2,work,out,6);
  cudaThreadSynchronize();

  magicfilter1d_t<T> <<< gridMF2, threadsMF2 >>>(2*n2,4*n1*n3,out,work,7);
  cudaThreadSynchronize();

  magicfilter1d_t<T> <<< gridMF1, threadsMF1 >>>(2*n1,4*n2*n3,work,out,8);
  cudaThreadSynchronize();

  //calculate potential energy
  *epot = do_blasDot(8*n1*n2*n3, psi, out);
 
  return 0;

}

template<typename T>
int kineticapplication(int n1,int n2, int n3,
		       T h1,T h2,T h3,
		       dim3 gridK1, dim3 gridK2, dim3 gridK3,
		       dim3 threadsK1, dim3 threadsK2, dim3 threadsK3,
		       T *psi,T *vpsi,T *out,T *ekinpot,
		       T *work)
{
  //define the scale factor to be applied to the convolution
  T scale=0.5/(h3*h3);
  kinetic1d<T> <<< gridK3, threadsK3 >>>(2*n3,4*n1*n2,scale,psi,work,vpsi,out,3);
  cudaThreadSynchronize();

  scale=0.5/(h2*h2);
  kinetic1d<T> <<< gridK2, threadsK2 >>>(2*n2,4*n1*n3,scale,work,psi,out,vpsi,4);
  cudaThreadSynchronize();

  scale=0.5/(h1*h1);
  kinetic1d<T> <<< gridK1, threadsK1 >>>(2*n1,4*n2*n3,scale,psi,work,vpsi,out,5);
  cudaThreadSynchronize();

  //calculate potential energy
  *ekinpot = do_blasDot(8*n1*n2*n3, out, work);

  return 0;
}

template<typename T>
int completelocalhamiltonian(int n1,int n2, int n3,
			     T h1,T h2,T h3,
			     T *psi,T *out,T *pot,			     
			     T *psifscf,
			     T *vpsifscf,
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

  waveletstosf<T>(n1,n2,n3,gridWT1,gridWT2,gridWT3,
		  threadsWT1,threadsWT2,threadsWT3,
		  psi,psifscf);


  potentialapplication<T>(n1,n2,n3,gridMF1,gridMF2,gridMF3,
			  threadsMF1,threadsMF2,threadsMF3,
			  psifscf,vpsifscf,pot,epot,
			  out); //work array

  kineticapplication<T>(n1,n2,n3,h1,h2,h3,
			gridK1,gridK2,gridK3,
			threadsK1,threadsK2,threadsK3,
			psifscf,vpsifscf,psi,ekinpot,
			out); //work array

  *ekinpot -= *epot;

  sftowavelets<T>(n1,n2,n3,gridWT1,gridWT2,gridWT3,
		  threadsWT1,threadsWT2,threadsWT3,
		  psi,out);


  return 0;

}


template<typename T>
int fulllocalhamiltonian(int n1,int n2, int n3,
			 T h1,T h2,T h3,
			 T *psiw,T *hpsiw,T *pot,int *keys,
			 T *psi,T *out, 
			 T *psifscf,
			 T *vpsifscf,
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
  cudaThreadSynchronize();

  uncompresscoarsefine<T> <<< gridC, threadsC >>>(n1,n2,n3,psiw,psi,keys);
  cudaThreadSynchronize();
  

  //cudaMemset((void*) psi,1,8*n1*n2*n3*sizeof(T));
  //cudaThreadSynchronize();


  //start the hamiltonian calculation
  waveletstosf<T>(n1,n2,n3,gridWT1,gridWT2,gridWT3,
		  threadsWT1,threadsWT2,threadsWT3,
		  psi,psifscf);


  potentialapplication<T>(n1,n2,n3,gridMF1,gridMF2,gridMF3,
			  threadsMF1,threadsMF2,threadsMF3,
			  psifscf,vpsifscf,pot,epot,
			  out); //work array

  //here the worky array should be initialised to c*x
  //c_initialize<T> <<< gridK3, threadsK3 >>>(2*n3,4*n1*n2,psi,work2,0.,3);
  //cudaThreadSynchronize();

  kineticapplication<T>(n1,n2,n3,h1,h2,h3,
			gridK1,gridK2,gridK3,
			threadsK1,threadsK2,threadsK3,
			psifscf,vpsifscf,psi,ekinpot,
			out); //work array

  *ekinpot -= *epot;

  sftowavelets<T>(n1,n2,n3,gridWT1,gridWT2,gridWT3,
		  threadsWT1,threadsWT2,threadsWT3,
		  psi,out);


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

extern "C" 
void kineticterm_(int *n1,int *n2,int *n3,
		  float *hx,float *hy,float *hz,float *c,
		  float **x,float **y,float **workx,float **worky,
		  float *ekin) 

{

  
  if(kineticfilter<float>(*n1+1,*n2+1,*n3+1,
			  *hx,*hy,*hz,*c,
			  *x,*workx,*y,*worky,
			  ekin) != 0)
    {
      printf("ERROR: GPU kineticfilter\n ");
      return;
    } 
  return; 
}

extern "C" 
void kinetictermd_(int *n1,int *n2,int *n3,
		   double *hx,double *hy,double *hz,double *c,
		   double **x,double **y,double **workx,double **worky,
		   double *ekin) 

{

  
  if(kineticfilter<double>(*n1+1,*n2+1,*n3+1,
			  *hx,*hy,*hz,*c,
			  *x,*workx,*y,*worky,
			  ekin) != 0)
    {
      printf("ERROR: GPU kineticfilter\n ");
      return;
    } 
  return; 
  }

extern "C" 
void kinetic1d_(int *n,int *ndat,
		double *h,double *c,
		double **x,double **y,double **workx,double **worky,
		double *ekin) 

{

  
  if(k1d<double>(*ndat,*n+1,
		 *h,*c,
		 *x,*workx,*y,*worky,ekin) != 0)
    {
      printf("ERROR: GPU kineticfilter\n ");
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





