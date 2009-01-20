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
#include "structDef_anasyn.h"

__constant__ parGPU_t par[3];

#include "kernels_anasyn.hcu"

// parameters to be used for calculating the convolution
template<typename T>
void WTParameters(parGPU_t* par,
		  unsigned int* num_halfwarps,
		  int n,
		  int ndat,
		  unsigned int* linecuts,
		  unsigned int* num_blocks)

{

  //number of total allowed elements of a input line
  unsigned int num_elem_tot=MAX_SHARED_SIZE/sizeof(T)/NUM_LINES/2; //between1024and64
  
  //number of elements of the output
  unsigned int num_elem_max=min(num_elem_tot-LOWFILWT-LUPFILWT-1,n); //between 1008 and 48 for 16-fil

  //number of pieces in which a line is divided
  //if the line is too small and not a multiple of ElementsPerHalfWarp
  //divide the line in two
  *linecuts=
    (n <= num_elem_max && n % HW_ELEM !=0 ? 2 : (n-1)/num_elem_max+1);

  //number of blocks in ndat direction
  *num_blocks=((ndat-1)/NUM_LINES + 1);

  //number of elements treated by each block 
  //this may pose problems for values of n dimensions less than 48
  //when n is not a multiple of ElementsPerHalfWarp
  //to be controlled for wavelet transform
  par->ElementsPerBlock = 
    min(HW_ELEM*(((n-1)/(*linecuts))/HW_ELEM+1),n);

  int halfwarps=16;
  //calculate the maximum number of halfwarps (between 4 and 16)
  for(int i =3; i>=0; --i)
    {
      if(par->ElementsPerBlock/HW_ELEM >= 1 << i)
	{
	  halfwarps = 1 << i;
	  break;
	}
    }

  *num_halfwarps = halfwarps;


  //printf("num_elem_tot %i,num_elem_max %i,linecuts %i,num_blocks %i,elemperBL %i, halfwarps %i,\n",
  //num_elem_tot,num_elem_max,*linecuts,*num_blocks, par->ElementsPerBlock,halfwarps);


  for(int j=0;j < HALF_WARP_SIZE ; ++j)
    {
      par->thline[j]= j & (NUM_LINES - 1); //num_lines always a power of two 
      par->thelem[j]= j / NUM_LINES; 
    }

  //define the sequences of the number of elements
  correctSequence(halfwarps,par->ElementsPerBlock/HW_ELEM,par->hwelem_calc);

  correctSequence(halfwarps,(par->ElementsPerBlock+LOWFILWT+LUPFILWT+1)/HW_ELEM,
		  par->hwelem_copy);

  //define the offsets
  for(int j=0,pos_calc=0,pos_copy=0;j < halfwarps ; ++j)
    {
      par->hwoffset_calc[j]=pos_calc;
      par->hwoffset_copy[j]=pos_copy;
      pos_calc+=HW_ELEM*par->hwelem_calc[j];
      pos_copy+=HW_ELEM*par->hwelem_copy[j];
      //printf("j %i,offset_calc %i,offset_copy %i,elem_calc %i,elem_copy %i \n",
      //j,par->hwoffset_calc[j],par->hwoffset_copy[j],
      //par->hwelem_calc[j],par->hwelem_copy[j]);

    }
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
  WTParameters<T>(&parCPU[2],&num_halfwarps,n3,n1*n2,&linecuts,&numBlocks);
  dim3  grid3(linecuts,  numBlocks, 1);  
  dim3  threads3(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n3,n1*n2);

  WTParameters<T>(&parCPU[1],&num_halfwarps,n2,n1*n3,&linecuts,&numBlocks);
  dim3  grid2(linecuts,  numBlocks, 1);  
  dim3  threads2(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n2,n1*n3);

  WTParameters<T>(&parCPU[0],&num_halfwarps,n1,n2*n3,&linecuts,&numBlocks);
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
  WTParameters<T>(&parCPU[0],&num_halfwarps,n,ndat,&linecuts,&numBlocks);
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
  WTParameters<T>(&parCPU[0],&num_halfwarps,n,ndat,&linecuts,&numBlocks);
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


//MUST CONTINUE WITH #D ANALYSIS SYNTHESIS VERSION

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
void somenamed_(int *n1,int *n2,int *n3,
		     double **psi,double **out) 

{

  
  if(wavana<double>(*n1+1,*n2+1,*n3+1,
	       *psi,
	       *out) != 0) 
    {
      printf("ERROR: GPU waveletanalysis\n ");
      return;
    } 
  return; 
}


/****/

