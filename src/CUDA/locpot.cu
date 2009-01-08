/****u* CUDA/1Dconv_new.cu
**
** 
** AUTHOR
**  Luigi Genovese
**
** SOURCE
*/
  
#include <stdio.h>
#include "locpot.h" //function declaration
#include "structDef_locpot.h"
#include "commonDef.h"


__constant__ parMF_t par[3]; //must be defined here because it is used by both by this file and kernels_locpot.cu

#include "reduction.hcu"
#include "kernels_locpot.hcu"

#include "reduction.h"





extern "C" 
void localpotential_(int *n1,
		     int *n2,
		     int *n3,
		     float **psi,
		     float **work,
		     float **pot,
		     float *epot) 

{  
  if(magicfilterpot(*n1+1,*n2+1,*n3+1,
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






// Magic Filter parameters to be used for calculating the convolution
void MFParameters(parMF_t* par,
		  unsigned int* num_halfwarps,
		  //unsigned int num_lines,
		  int n,
		  int ndat,
		  unsigned int* linecuts,
		  unsigned int* num_blocks)

{

  //number of total allowed elements of a input line
  unsigned int num_elem_tot=MAX_SHARED_SIZE/sizeof(float)/NUM_LINES; //between1024and64
  
  //number of elements of the output
  unsigned int num_elem_max=min(num_elem_tot-LOWFIL-LUPFIL-1,n); //between 1008 and 48 for 16-fil

  //number of pieces in which a line is divided
  //if the line is too small and not a multiple of ElementsPerHalfWarp
  //divide the line in two
  *linecuts=
    (n <= num_elem_max && n % HW_ELEM !=0 ? 2 : (n-1)/num_elem_max+1);

  //number of blocks in ndat direction
  *num_blocks=((ndat-1)/NUM_LINES + 1);

  //printf("num_elem_tot %i,num_elem_max %i,linecuts %i,num_blocks %i,elemperHW %i \n",
  //num_elem_tot,num_elem_max,*linecuts,*num_blocks, par -> ElementsPerHalfWarp);

  //number of elements treated by each block 
  //this may pose problems for values of n dimensions less than 48
  //when n is not a multiple of ElementsPerHalfWarp
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

  for(int j=0;j < HALF_WARP_SIZE ; ++j)
    {
      par->thline[j]= j & (NUM_LINES - 1); //num_lines always a power of two 
      par->thelem[j]= j / NUM_LINES; 
    }

  //define the sequences of the number of elements
  correctSequence(halfwarps,par->ElementsPerBlock/HW_ELEM,par->hwelem_calc);

  correctSequence(halfwarps,(par->ElementsPerBlock+LOWFIL+LUPFIL+1)/HW_ELEM,
		  par->hwelem_copy);

  //define the offsets
  for(int j=0,pos_calc=0,pos_copy=0;j < halfwarps ; ++j)
    {
      par->hwoffset_calc[j]=pos_calc;
      par->hwoffset_copy[j]=pos_copy;
      pos_calc+=HW_ELEM*par->hwelem_calc[j];
      pos_copy+=HW_ELEM*par->hwelem_copy[j];
    }
 
  //printf("ElementsPerBlock %i,HalfWarpCalculatedElements %i,HalfWarpCopiedElements %i,LastHalfWarpCalcElements %i, LastHalfWarpCopiedElements %i \n",

}










int magicfilterpot(int n1,int n2, int n3,
		   float *psi,
		   float *work,
		   float *pot,
		   float *epot)
{

  //create the parameters
  parMF_t parCPU[3];

  //calculate the number of threads and blocks
  //unsigned int num_lines = min(16,ndat); //hard coded for the moment
  unsigned int numBlocks,linecuts,num_halfwarps;

  //calculate the parameters in constant memory for each of the 1D convolution
  //define the number of threads and blocks according to parameter definitions
  MFParameters(&parCPU[2],&num_halfwarps,n3,n1*n2,&linecuts,&numBlocks);
  dim3  grid3(linecuts,  numBlocks, 1);  
  dim3  threads3(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n3,n1*n2);

  MFParameters(&parCPU[1],&num_halfwarps,n2,n1*n3,&linecuts,&numBlocks);
  dim3  grid2(linecuts,  numBlocks, 1);  
  dim3  threads2(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n2,n1*n3);

  MFParameters(&parCPU[0],&num_halfwarps,n1,n2*n3,&linecuts,&numBlocks);
  dim3  grid1(linecuts,  numBlocks, 1);  
  dim3  threads1(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n1,n3*n2);

  //send them to constant memory, once and for all
  if(cudaMemcpyToSymbol(*par,&parCPU, 3*sizeof(parMF_t)) != 0)
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
  magicfilter1d <<< grid3, threads3 >>>(n3,n1*n2,psi,work,2);
  //unbind the texture
  //cudaUnbindTexture(psi_tex);

  cudaThreadSynchronize();

  magicfilter1d <<< grid2, threads2 >>>(n2,n1*n3,work,psi,1);
  cudaThreadSynchronize();

  magicfilter1d_pot <<< grid1, threads1 >>>(n1,n2*n3,psi,pot,work,0);
  cudaThreadSynchronize();

  //here one should combine psi and work to calculate the potential
  //energy
  reducearrays(n1,n2*n3,psi,work,epot);
  cudaThreadSynchronize(); //can be removed since the arrays are not overwritten

  //reverse MF calculation
  magicfilter1d_t <<< grid3, threads3 >>>(n3,n1*n2,work,psi,2);
  cudaThreadSynchronize();

  magicfilter1d_t <<< grid2, threads2 >>>(n2,n1*n3,psi,work,1);
  cudaThreadSynchronize();

  magicfilter1d_t <<< grid1, threads1 >>>(n1,n2*n3,work,psi,0);
  cudaThreadSynchronize();

  return 0;

}


