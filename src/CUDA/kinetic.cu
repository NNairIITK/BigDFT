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
/*
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
*/

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
