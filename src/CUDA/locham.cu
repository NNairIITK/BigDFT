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

__constant__ parGPU_t par[9];

#include "kernels_anasyn.hcu"
#include "kernels_locpot.hcu"
#include "kernels_kinetic.hcu"

// parameters to be used for calculating the convolution
template<typename T>
void GPUParameters(parGPU_t* par,
		   unsigned int* num_halfwarps,
		   int n,
		   int ndat,
		   int nfac,
		   int lowfil,int lupfil,
		   unsigned int* linecuts,
		   unsigned int* num_blocks)

{

  //number of total allowed elements of a input line
  //nfac added for wavelet transformation routines
  unsigned int num_elem_tot=MAX_SHARED_SIZE/sizeof(T)/NUM_LINES/nfac; //between1024and64
  
  //number of elements of the output
  unsigned int num_elem_max=min(num_elem_tot-lowfil-lupfil-1,n); //between 1008 and 48 for 16-fil

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

  correctSequence(halfwarps,(par->ElementsPerBlock+lowfil+lupfil+1)/HW_ELEM,
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
  reducearrays<T>(2*n1,4*n2*n3,work2,work,epot);
  cudaThreadSynchronize();

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
  reducearrays<T>(2*n1,4*n2*n3,work,out,ekinpot);
  cudaThreadSynchronize();

  //wavelet analysis
  waveletanalysis<T> <<< gridWT3, threadsWT3 >>>(n3,4*n1*n2,out,work,0);
  cudaThreadSynchronize();

  waveletanalysis<T> <<< gridWT2, threadsWT2 >>>(n2,4*n1*n3,work,work2,1);
  cudaThreadSynchronize();

  waveletanalysis<T> <<< gridWT1, threadsWT1 >>>(n1,4*n2*n3,work2,out,2);
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

