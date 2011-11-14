#ifndef __GPUparameters__
#define __GPUparameters__

//filters for second derivative in daubechies basis
#define A2  3.55369228991319019
#define B2  24.8758460293923314

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

  //extra terms which should be added for the copy
  unsigned int nextra=HW_ELEM*((lowfil+lupfil)/HW_ELEM+1);

  //number of total allowed elements of a input line
  //nfac added for wavelet transformation routines
  unsigned int num_elem_tot=MAX_SHARED_SIZE/sizeof(T)/NUM_LINES/nfac; //between1024and64
  
  //number of elements of the output
  unsigned int num_elem_max=min(num_elem_tot-nextra,n); //between 1008 and 48 for 16-fil

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

  correctSequence(halfwarps,par->ElementsPerBlock/HW_ELEM
		  +nextra/HW_ELEM, //exceeding term for the copy
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

// parameters for the preconditioning
template<typename T>
void GPUprecondparameters(T *GPUscal,T hx,T hy, T hz)
{
  T h1,h2,h3;

  h1= ((T)0.5)/(hx*hx);
  h2= ((T)0.5)/(hy*hy);
  h3= ((T)0.5)/(hz*hz);

  GPUscal[0]=A2*h1+A2*h2+A2*h3;
  GPUscal[1]=B2*h1+A2*h2+A2*h3;
  GPUscal[2]=A2*h1+B2*h2+A2*h3;
  GPUscal[3]=B2*h1+B2*h2+A2*h3;
  GPUscal[4]=A2*h1+A2*h2+B2*h3;
  GPUscal[5]=B2*h1+A2*h2+B2*h3;
  GPUscal[6]=A2*h1+B2*h2+B2*h3;
  GPUscal[7]=B2*h1+B2*h2+B2*h3;

}

#endif
