/****u* CUDA/1Dconv_new.cu
**
** 
** AUTHOR
**  Luigi Genovese
**
** SOURCE
*/

#include <stdio.h>
#include <cutil.h>
#include <pthread.h>
#include <semaphore.h>
#include <sched.h>

#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)

//maximum size of the shared memory array
//conceived for maximize occupancy on a hardware of compute
//capability 1.2 and higher (1024 threads at same time on a given multiprocessor)
#define MAX_SHARED_SIZE 3072 //16*256 4 kB (should be =~ 3.9 kB, try also 3072)
#define HALF_WARP_SIZE 16 // for all architectures
#define NUM_LINES 16 
#define HW_ELEM 1 //this is HALF_WARP_SIZE/NUM_LINES

//parameter related to the Magic Filter convolution
//lowfil + lupfil + 1  must be a multiple of 16
#define LOWFIL 8
#define LUPFIL 7

//convolution filters
#define MFIL0   8.4334247333529341094733325815816e-7f
#define MFIL1  -0.1290557201342060969516786758559028e-4f
#define MFIL2   0.8762984476210559564689161894116397e-4f
#define MFIL3  -0.30158038132690463167163703826169879e-3f
#define MFIL4   0.174723713672993903449447812749852942e-2f
#define MFIL5  -0.942047030201080385922711540948195075e-2f
#define MFIL6   0.2373821463724942397566389712597274535e-1f
#define MFIL7   0.612625895831207982195380597e-1f
#define MFIL8   0.9940415697834003993178616713f
#define MFIL9  -0.604895289196983516002834636e-1f
#define MFIL10 -0.2103025160930381434955489412839065067e-1f
#define MFIL11  0.1337263414854794752733423467013220997e-1f
#define MFIL12 -0.344128144493493857280881509686821861e-2f
#define MFIL13  0.49443227688689919192282259476750972e-3f
#define MFIL14 -0.5185986881173432922848639136911487e-4f
#define MFIL15  2.72734492911979659657715313017228e-6f

#include "reduction.h"
 

typedef struct  _parMF
{
  unsigned int ElementsPerBlock;

  int thline[HALF_WARP_SIZE]; //line considered by a thread within the half-warp
  int thelem[HALF_WARP_SIZE]; //elements considered by a thread within the half-warp
  int hwelem_calc[16]; //maximum number of half warps
  int hwelem_copy[16]; //maximum number of half-warps
  int hwoffset_calc[16]; //maximum number of half warps
  int hwoffset_copy[16]; //maximum number of half-warps

} parMF_t;

__constant__ parMF_t par[3];

//declare the texture for binding the input psi
//texture<float> psi_tex;

void correctSequence(int thds,int elem,int * tab);

int magicfilterpot(int n1,int n2, int n3,
		   float *psi,
		   float *work,
		   float *pot,
		   float epot);



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

//uniformise the tabular of the number of elements treated by each
//thread (counterpart of uniformiseTab)
void correctSequence(int thds,int elem,int * tab)
{
  //put to zero all the values;
  for(int j=0;j< elem; ++j)
    {
      tab[j]=0;
    }

  //then start to fill consecutively until reaching of the end
  //if elem > thds no element will be zero
  //this is the most balanced choice
  for(int i=0;i< elem; ++i)
    {
      tab[i % thds]+=1;
    }
}



//1D convolution of multiple lines in the same block
//__global__ void magicfilter1d(int n,int ndat, float *psi_out) //for textures
__global__ void magicfilter1d(int n,int ndat, float *psi_in, float *psi_out,int idim)
{

  //line treated by the given block
  unsigned int lineOffset = min(blockIdx.y*NUM_LINES,ndat-NUM_LINES);
  //starting element treated by the block
  unsigned int elemOffset = min(blockIdx.x*par[idim].ElementsPerBlock,n-par[idim].ElementsPerBlock);

  //half-warp id
  const unsigned int hwid = threadIdx.y;
  //tid within the HW
  const unsigned int tid_hw = threadIdx.x;

  //shared memory array
  __shared__ float psi_sh[MAX_SHARED_SIZE/sizeof(float)];

  //line treated by the given thread in ndat axis
  //which is the input base element
  unsigned int BaseElem = par[idim].thline[tid_hw] + lineOffset;
  //write data in shared memory
  //element treated by the given thread in n-axis
  unsigned int thelem = par[idim].thelem[tid_hw] + par[idim].hwoffset_copy[hwid];

  unsigned int ShBaseElem = tid_hw + NUM_LINES*par[idim].hwoffset_copy[hwid];

  int epsilon,npos;

  //NOTE: it is assumed that for non-first segments the starting
  //points is far enough for the filter to be contained
  //and the same for non-last segments.
  //in other terms: lenght of the line is always bigger than
  //max(lowfil,lupfil)

  for(int i=0,ipos=elemOffset-LOWFIL+thelem;i < par[idim].hwelem_copy[hwid] ; ++i)
    {
      epsilon=(ipos < 0 ? -1 : ipos/n);
      npos=ipos-epsilon*n;
      psi_sh[ShBaseElem]=psi_in[BaseElem+ndat*npos];
      //psi_sh[ShBaseElem]=tex1Dfetch(psi_tex,BaseElem+ndat*npos);

      ShBaseElem += HALF_WARP_SIZE;
      ipos += HW_ELEM;
      
    }

  //end shared memory copy
  __syncthreads();

  //element treated by the given thread in n-axis
  thelem = par[idim].thelem[tid_hw] + par[idim].hwoffset_calc[hwid];
  //base element for the given thread in shared memory
  ShBaseElem = tid_hw + NUM_LINES*par[idim].hwoffset_calc[hwid];

  //output base element, from the input one
  BaseElem =  n*BaseElem+ thelem + elemOffset;

  //perform convolution in shared memory 
  //each thread calculate a number of elements, identical for each
  //half-warp
  //#pragma unroll 5 (to be tested if it is important)

  for(int i=0;i < par[idim].hwelem_calc[hwid]; ++i)
    {
      //values of the convolution
      register float conv = 
	//hand-unrolled loop (16 elements for this filter)
	//order changed for increasing the precision
	MFIL0 *psi_sh[ShBaseElem               ] +
	MFIL15*psi_sh[ShBaseElem + 15*NUM_LINES] +
	MFIL1 *psi_sh[ShBaseElem +   NUM_LINES ] +
	MFIL14*psi_sh[ShBaseElem + 14*NUM_LINES] +
	MFIL2 *psi_sh[ShBaseElem + 2*NUM_LINES ] +
	MFIL13*psi_sh[ShBaseElem + 13*NUM_LINES] +
	MFIL3 *psi_sh[ShBaseElem + 3*NUM_LINES ] +
	MFIL12*psi_sh[ShBaseElem + 12*NUM_LINES] +
	MFIL4 *psi_sh[ShBaseElem + 4*NUM_LINES ] +
	MFIL11*psi_sh[ShBaseElem + 11*NUM_LINES] +
	MFIL5 *psi_sh[ShBaseElem + 5*NUM_LINES ] +
	MFIL10*psi_sh[ShBaseElem + 10*NUM_LINES] +
	MFIL6 *psi_sh[ShBaseElem + 6*NUM_LINES ] +
	MFIL9 *psi_sh[ShBaseElem + 9*NUM_LINES ] +
	MFIL7 *psi_sh[ShBaseElem + 7*NUM_LINES ] +
	MFIL8 *psi_sh[ShBaseElem + 8*NUM_LINES ] ;

      psi_out[BaseElem]=conv;
      //psi_sh[ShBaseElem+LOWFIL*par[idim].LinesPerBlock]; //for testing only

      ShBaseElem += HALF_WARP_SIZE;
      BaseElem += HW_ELEM;
      
    }

}

//1D convolution of multiple lines in the same block
//multiplies by the potential and calculate the potential energy
//__global__ void magicfilter1d_pot(int n,int ndat, float *psi_out)
__global__ void magicfilter1d_pot(int n,int ndat, float *psi_in, 
				  float *pot, float *psi_out,int idim)
{

  //line treated by the given block
  unsigned int lineOffset = min(blockIdx.y*NUM_LINES,ndat-NUM_LINES);
  //starting element treated by the block
  unsigned int elemOffset = min(blockIdx.x*par[idim].ElementsPerBlock,n-par[idim].ElementsPerBlock);

  //half-warp id
  const unsigned int hwid = threadIdx.y;
  //tid within the HW
  const unsigned int tid_hw = threadIdx.x;

  //shared memory array
  __shared__ float psi_sh[MAX_SHARED_SIZE/sizeof(float)];

  //line treated by the given thread in ndat axis
  //which is the input base element
  unsigned int BaseElem = par[idim].thline[tid_hw] + lineOffset;
  //write data in shared memory
  //element treated by the given thread in n-axis
  unsigned int thelem = par[idim].thelem[tid_hw] + par[idim].hwoffset_copy[hwid];

  unsigned int ShBaseElem = tid_hw + NUM_LINES*par[idim].hwoffset_copy[hwid];

  int epsilon,npos;

  //NOTE: it is assumed that for non-first segments the starting
  //points is far enough for the filter to be contained
  //and the same for non-last segments.
  //in other terms: lenght of the line is always bigger than
  //max(lowfil,lupfil)

  for(int i=0,ipos=elemOffset-LOWFIL+thelem;i < par[idim].hwelem_copy[hwid] ; ++i)
    {
      //control flag for periodic boundary conditions
      epsilon=(ipos < 0 ? -1 : ipos/n);
      npos=ipos-epsilon*n;

      psi_sh[ShBaseElem]=psi_in[BaseElem+ndat*npos];
      //psi_sh[ShBaseElem]=tex1Dfetch(psi_tex,BaseElem+ndat*npos);

      ShBaseElem += HALF_WARP_SIZE;
      ipos += HW_ELEM;
    }

  //end shared memory copy
  __syncthreads();

  //element treated by the given thread in n-axis
  thelem = par[idim].thelem[tid_hw] + par[idim].hwoffset_calc[hwid];
  //base element for the given thread in shared memory
  ShBaseElem = tid_hw + NUM_LINES*par[idim].hwoffset_calc[hwid];

  //output base element, from the input one
  BaseElem =  n*BaseElem+ thelem + elemOffset;

  //limit element for which the block treats unique elements

  //perform convolution in shared memory 
  //each thread calculate a number of elements, identical for each
  //half-warp

  for(int i=0;i < par[idim].hwelem_calc[hwid]; ++i)
    {
      //values of the convolution
      register float conv = 
	//hand-unrolled loop (16 elements for this filter)
	//order changed for increasing the precision
	MFIL0 *psi_sh[ShBaseElem               ] +
	MFIL15*psi_sh[ShBaseElem + 15*NUM_LINES] +
	MFIL1 *psi_sh[ShBaseElem +   NUM_LINES ] +
	MFIL14*psi_sh[ShBaseElem + 14*NUM_LINES] +
	MFIL2 *psi_sh[ShBaseElem + 2*NUM_LINES ] +
	MFIL13*psi_sh[ShBaseElem + 13*NUM_LINES] +
	MFIL3 *psi_sh[ShBaseElem + 3*NUM_LINES ] +
	MFIL12*psi_sh[ShBaseElem + 12*NUM_LINES] +
	MFIL4 *psi_sh[ShBaseElem + 4*NUM_LINES ] +
	MFIL11*psi_sh[ShBaseElem + 11*NUM_LINES] +
	MFIL5 *psi_sh[ShBaseElem + 5*NUM_LINES ] +
	MFIL10*psi_sh[ShBaseElem + 10*NUM_LINES] +
	MFIL6 *psi_sh[ShBaseElem + 6*NUM_LINES ] +
	MFIL9 *psi_sh[ShBaseElem + 9*NUM_LINES ] +
	MFIL7 *psi_sh[ShBaseElem + 7*NUM_LINES ] +
	MFIL8 *psi_sh[ShBaseElem + 8*NUM_LINES ] ;

      //register float v=tex1Dfetch(pot_tex,BaseElem);

      psi_out[BaseElem]=conv*pot[BaseElem];

      ShBaseElem += HALF_WARP_SIZE;
      BaseElem += HW_ELEM;
      
    }
 
}

//transposed convolution
__global__ void magicfilter1d_t(int n,int ndat, float *psi_in, float *psi_out,int idim)
{

  //line treated by the given block
  unsigned int lineOffset = min(blockIdx.y*NUM_LINES,ndat-NUM_LINES);
  //starting element treated by the block
  unsigned int elemOffset = 
    min(blockIdx.x*par[idim].ElementsPerBlock,n-par[idim].ElementsPerBlock);

  //half-warp id
  const unsigned int hwid = threadIdx.y;
  //tid within the HW
  const unsigned int tid_hw = threadIdx.x;

  //shared memory array
  __shared__ float psi_sh[MAX_SHARED_SIZE/sizeof(float)];

  //line treated by the given thread in ndat axis
  //which is the input base element
  unsigned int BaseElem = par[idim].thline[tid_hw] + lineOffset;
  //write data in shared memory
  //element treated by the given thread in n-axis
  unsigned int thelem = par[idim].thelem[tid_hw] + par[idim].hwoffset_copy[hwid];

  unsigned int ShBaseElem = tid_hw + NUM_LINES*par[idim].hwoffset_copy[hwid];

  int epsilon,npos;

  //NOTE: it is assumed that for non-first segments the starting
  //points is far enough for the filter to be contained
  //and the same for non-last segments.
  //in other terms: lenght of the line is always bigger than
  //max(lowfil,lupfil)

  for(int i=0,ipos=elemOffset-LUPFIL+thelem;i < par[idim].hwelem_copy[hwid] ; ++i)
    {
      epsilon=(ipos < 0 ? -1 : ipos/n);
      npos=ipos-epsilon*n;
      psi_sh[ShBaseElem]=psi_in[BaseElem+ndat*npos];
      //psi_sh[ShBaseElem]=tex1Dfetch(psi_tex,BaseElem+ndat*npos);

      ShBaseElem += HALF_WARP_SIZE;
      ipos += HW_ELEM;
      
    }

  //end shared memory copy
  __syncthreads();

  //element treated by the given thread in n-axis
  thelem = par[idim].thelem[tid_hw] + par[idim].hwoffset_calc[hwid];
  //base element for the given thread in shared memory
  ShBaseElem = tid_hw + NUM_LINES*par[idim].hwoffset_calc[hwid];

  //output base element, from the input one
  BaseElem =  n*BaseElem+ thelem + elemOffset;

  //perform convolution in shared memory 
  //each thread calculate a number of elements, identical for each
  //half-warp
  //#pragma unroll 5 (to be tested if it is important)

  for(int i=0;i < par[idim].hwelem_calc[hwid]; ++i)
    {
      //values of the convolution
      register float conv = 
	//hand-unrolled loop (16 elements for this filter)
	//order changed for increasing the precision
	MFIL15*psi_sh[ShBaseElem               ] +
	MFIL0 *psi_sh[ShBaseElem + 15*NUM_LINES] +
	MFIL14*psi_sh[ShBaseElem +   NUM_LINES ] +
	MFIL1 *psi_sh[ShBaseElem + 14*NUM_LINES] +
	MFIL13*psi_sh[ShBaseElem + 2*NUM_LINES ] +
	MFIL2 *psi_sh[ShBaseElem + 13*NUM_LINES] +
	MFIL12*psi_sh[ShBaseElem + 3*NUM_LINES ] +
	MFIL3 *psi_sh[ShBaseElem + 12*NUM_LINES] +
	MFIL11*psi_sh[ShBaseElem + 4*NUM_LINES ] +
	MFIL4 *psi_sh[ShBaseElem + 11*NUM_LINES] +
	MFIL10*psi_sh[ShBaseElem + 5*NUM_LINES ] +
	MFIL5 *psi_sh[ShBaseElem + 10*NUM_LINES] +
	MFIL9 *psi_sh[ShBaseElem + 6*NUM_LINES ] +
	MFIL6 *psi_sh[ShBaseElem + 9*NUM_LINES ] +
	MFIL8 *psi_sh[ShBaseElem + 7*NUM_LINES ] +
	MFIL7 *psi_sh[ShBaseElem + 8*NUM_LINES ] ;

      psi_out[BaseElem]=conv;
      //psi_sh[ShBaseElem+LOWFIL*par[idim].LinesPerBlock]; //for testing only

      ShBaseElem += HALF_WARP_SIZE;
      BaseElem += HW_ELEM;
      
    }

 
}

extern "C" 
void localpotential_(int *n1,int *n2,int *n3,
		     float **psi,float **work,float **pot,
		     float *epot) 

{

  
  if(magicfilterpot(*n1+1,*n2+1,*n3+1,
	       *psi,
	       *work,
	       *pot,
	       *epot) != 0)
    {
      printf("ERROR: GPU magicfilterpot\n ");
      return;
    } 
  return; 
}


int magicfilterpot(int n1,int n2, int n3,
		   float *psi,
		   float *work,
		   float *pot,
		   float epot)
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
  epot=reducearrays(n1,n2*n3,psi,work);

  //reverse MF calculation
  magicfilter1d_t <<< grid3, threads3 >>>(n3,n1*n2,work,psi,2);
  cudaThreadSynchronize();

  magicfilter1d_t <<< grid2, threads2 >>>(n2,n1*n3,psi,work,1);
  cudaThreadSynchronize();

  magicfilter1d_t <<< grid1, threads1 >>>(n1,n2*n3,work,psi,0);
  cudaThreadSynchronize();

  return 0;

}

/****/
