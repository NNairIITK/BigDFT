/****u* CUDA/kinetic.cu
**
** 
** AUTHOR
**  Luigi Genovese
**
** SOURCE
*/
  
#include <stdio.h>

#include "kinetic.h"

#include "commonDef.h"


#include "reduction.h"


#include "structDef_kinetic.h"



__constant__ parGPU_t par[3];


#include "kernels_kinetic.hcu"

// Kinetic parameters to be used for calculating the convolution
template<typename T>
void KParameters(parGPU_t* par,
		 unsigned int* num_halfwarps,
		 int n,
		 int ndat,
		 unsigned int* linecuts,
		 unsigned int* num_blocks)

{

  //number of total allowed elements of a input line
   int num_elem_tot = MAX_SHARED_SIZE/sizeof(T)/NUM_LINES; //between1024and64
  
  //number of elements of the output
   int num_elem_max = min(num_elem_tot-LOWFILK-LUPFILK-1,n); //between 996 and 35 

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
    min(HW_ELEM*(((n-1)/(int)(*linecuts))/HW_ELEM+1),n);

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

  correctSequence(halfwarps,(par->ElementsPerBlock+LOWFILK+LUPFILK+1)/HW_ELEM,
		  par->hwelem_copy);

  //define the offsets
  for(int j=0,pos_calc=0,pos_copy=0;j < halfwarps ; ++j)
    {
      par->hwoffset_calc[j]=pos_calc;
      par->hwoffset_copy[j]=pos_copy;
      pos_calc+=HW_ELEM*par->hwelem_calc[j];
      pos_copy+=HW_ELEM*par->hwelem_copy[j];
    }
 
}



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
  KParameters<T>(&parCPU[2],&num_halfwarps,n3,n1*n2,&linecuts,&numBlocks);
  dim3  grid3(linecuts,  numBlocks, 1);  
  dim3  threads3(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n3,n1*n2);

  KParameters<T>(&parCPU[1],&num_halfwarps,n2,n1*n3,&linecuts,&numBlocks);
  dim3  grid2(linecuts,  numBlocks, 1);  
  dim3  threads2(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n2,n1*n3);

  KParameters<T>(&parCPU[0],&num_halfwarps,n1,n2*n3,&linecuts,&numBlocks);
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
  reducearrays<T>(n1,n2*n3,x,y,ekin);
  cudaThreadSynchronize();

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
  KParameters<T>(&parCPU[0],&num_halfwarps,n,ndat,&linecuts,&numBlocks);
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
  reducearrays<T>(n,ndat,x,y,ekin);
  cudaThreadSynchronize();
  return 0;

}

/****/



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
