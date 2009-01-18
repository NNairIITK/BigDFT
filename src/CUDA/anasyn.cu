/****u* CUDA/anasyn.cu
**
** 
** AUTHOR
**  Luigi Genovese
**
** SOURCE
*/
  
#include <stdio.h>
//#include <pthread.h>
//#include <semaphore.h>
//#include <sched.h>

#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)

//maximum size of the shared memory array
//conceived for maximize occupancy on a hardware of compute
//capability 1.2 and higher (1024 threads at same time on a given multiprocessor)
#define MAX_SHARED_SIZE 3072 //16*256 4 kB (should be =~ 3.9 kB, try also 3072)
#define HALF_WARP_SIZE 16 // for all architectures
#define NUM_LINES 16 
#define HW_ELEM 1 //this is HALF_WARP_SIZE/NUM_LINES

//parameter related to convolution
#define LOWSYN 4
#define LUPSYN 4

#define LOWANA 7
#define LUPANA 8

/*
data ch  /  0.e0_wp            , -0.0033824159510050025955_wp, & 
  -0.00054213233180001068935_wp, 0.031695087811525991431_wp, & 
  0.0076074873249766081919_wp  , -0.14329423835127266284_wp, & 
  -0.061273359067811077843_wp  , 0.48135965125905339159_wp,  & 
  0.77718575169962802862_wp    ,0.36444189483617893676_wp, &
  -0.051945838107881800736_wp  ,-0.027219029917103486322_wp, &
  0.049137179673730286787_wp   ,0.0038087520138944894631_wp, &
  -0.014952258337062199118_wp  ,-0.00030292051472413308126_wp, &
  0.0018899503327676891843_wp  , 0.e0_wp /

  data cg  / 0.e0_wp           , -0.0018899503327676891843_wp, &
  -0.00030292051472413308126_wp, 0.014952258337062199118_wp, &
  0.0038087520138944894631_wp  , -0.049137179673730286787_wp, &
  -0.027219029917103486322_wp  , 0.051945838107881800736_wp, &
  0.36444189483617893676_wp    , -0.77718575169962802862_wp, &
  0.48135965125905339159_wp    , 0.061273359067811077843_wp, &
  -0.14329423835127266284_wp   , -0.0076074873249766081919_wp, &
  0.031695087811525991431_wp   , 0.00054213233180001068935_wp, &
  -0.0033824159510050025955_wp , 0.e0_wp /
*/

  //wavelet filters reordered for matching convolution 

  //even part of h filters, reverted
#define WTFIL0    0.0018899503327676891843f
#define WTFIL1   -0.014952258337062199118f
#define WTFIL2    0.049137179673730286787f
#define WTFIL3   -0.051945838107881800736f
#define WTFIL4    0.77718575169962802862f
#define WTFIL5   -0.061273359067811077843f
#define WTFIL6    0.0076074873249766081919f
#define WTFIL7   -0.00054213233180001068935f
  //even part of g filters, reverted (or odd part of h filters)
#define WTFIL8   -0.0033824159510050025955f 
#define WTFIL9    0.031695087811525991431f   
#define WTFIL10  -0.14329423835127266284f   
#define WTFIL11   0.48135965125905339159f    
#define WTFIL12   0.36444189483617893676f    
#define WTFIL13  -0.027219029917103486322f  
#define WTFIL14   0.0038087520138944894631f  
#define WTFIL15  -0.00030292051472413308126f



  do j=1,ndat

    do i=0,n
      se=0.e0_wp
      so=0.e0_wp
      do l=-4,4
	k=modulo(i-l,n+1)
     se=se+ch(2*l  )*x(  k,j)+cg(2*l  )*x(n+1+k  ,j)
     so=so+ch(2*l+1)*x(  k,j)+cg(2*l+1)*x(n+1+k  ,j)
     enddo
        y(j,2*i  )=se
        y(j,2*i+1)=so
	enddo

	enddo

typedef struct  _parWT
{
  unsigned int ElementsPerBlock;

  int thline[HALF_WARP_SIZE]; //line considered by a thread within the half-warp
  int thelem[HALF_WARP_SIZE]; //elements considered by a thread within the half-warp
  int hwelem_calc[16]; //maximum number of half warps
  int hwelem_copy[16]; //maximum number of half-warps
  int hwoffset_calc[16]; //maximum number of half warps
  int hwoffset_copy[16]; //maximum number of half-warps

} parWT_t;

__constant__ parWT_t par[3];

void correctSequence(int thds,int elem,int * tab);

int waveletransform(int n1,int n2, int n3,
		   float *psi,
		   float *work);

// Magic Filter parameters to be used for calculating the convolution
void WTParameters(parMF_t* par,
		  unsigned int* num_halfwarps,
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
}

//1D convolution of multiple lines in the same block
__global__ void waveletsynthesis(int n,int ndat, float *psi_in, float *psi_out,int idim)
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
  __shared__ float psi_sh1[MAX_SHARED_SIZE/sizeof(float)/2];
  __shared__ float psi_sh2[MAX_SHARED_SIZE/sizeof(float)/2];

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
      psi_sh1[ShBaseElem]=psi_in[BaseElem+ndat*npos];
      psi_sh2[ShBaseElem]=psi_in[BaseElem+ndat*(npos+n)];

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
  BaseElem =  (2*n+1)*BaseElem+ 2*(thelem + elemOffset);

  //perform convolution in shared memory 
  //each thread calculate a number of elements, identical for each
  //half-warp
  //#pragma unroll 5 (to be tested if it is important)

  for(int i=0;i < par[idim].hwelem_calc[hwid]; ++i)
    {
      //values of the convolution, even components
      register float conv = 
	//hand-unrolled loop 
	//order changed for increasing the precision
	WTFIL7 *psi_sh1[ShBaseElem + 7*NUM_LINES] +
	WTFIL15*psi_sh2[ShBaseElem + 7*NUM_LINES] +
	WTFIL0 *psi_sh1[ShBaseElem              ] +
	WTFIL8 *psi_sh2[ShBaseElem              ] +
	WTFIL6 *psi_sh1[ShBaseElem + 6*NUM_LINES] +
	WTFIL14*psi_sh2[ShBaseElem + 6*NUM_LINES] +
	WTFIL1 *psi_sh1[ShBaseElem + 1*NUM_LINES] +
	WTFIL9 *psi_sh2[ShBaseElem + 1*NUM_LINES] +
	WTFIL5 *psi_sh1[ShBaseElem + 5*NUM_LINES] +
	WTFIL13*psi_sh2[ShBaseElem + 5*NUM_LINES] +
	WTFIL2 *psi_sh1[ShBaseElem + 2*NUM_LINES] +
	WTFIL10*psi_sh2[ShBaseElem + 2*NUM_LINES] +
	WTFIL4 *psi_sh1[ShBaseElem + 4*NUM_LINES] +
	WTFIL12*psi_sh2[ShBaseElem + 4*NUM_LINES] +
	WTFIL3 *psi_sh1[ShBaseElem + 3*NUM_LINES] +
	WTFIL11*psi_sh2[ShBaseElem + 3*NUM_LINES] ;

      psi_out[BaseElem]=conv;


      //values of the convolution, odd components
      conv = 
	//hand-unrolled loop 
	//order changed for increasing the precision
	WTFIL15*psi_sh1[ShBaseElem + 1*NUM_LINES] -
	WTFIL7 *psi_sh2[ShBaseElem + 1*NUM_LINES] +
	WTFIL8 *psi_sh1[ShBaseElem + 8*NUM_LINES] -
	WTFIL0 *psi_sh2[ShBaseElem + 8*NUM_LINES] +
	WTFIL14*psi_sh1[ShBaseElem + 2*NUM_LINES] -
	WTFIL6 *psi_sh2[ShBaseElem + 2*NUM_LINES] +
	WTFIL9 *psi_sh1[ShBaseElem + 7*NUM_LINES] -
	WTFIL1 *psi_sh2[ShBaseElem + 7*NUM_LINES] +
	WTFIL13*psi_sh1[ShBaseElem + 3*NUM_LINES] -
	WTFIL5 *psi_sh2[ShBaseElem + 3*NUM_LINES] +
	WTFIL10*psi_sh1[ShBaseElem + 6*NUM_LINES] -
	WTFIL2 *psi_sh2[ShBaseElem + 6*NUM_LINES] +
	WTFIL12*psi_sh1[ShBaseElem + 4*NUM_LINES] -
	WTFIL4 *psi_sh2[ShBaseElem + 4*NUM_LINES] +
	WTFIL11*psi_sh1[ShBaseElem + 5*NUM_LINES] -
	WTFIL3 *psi_sh2[ShBaseElem + 5*NUM_LINES] ;

      psi_out[BaseElem+1]=conv;

      //psi_sh[ShBaseElem+LOWFIL*par[idim].LinesPerBlock]; //for testing only

      ShBaseElem += HALF_WARP_SIZE;
      BaseElem += 2*HW_ELEM;
      
    }

}

//transposed convolution
__global__ void waveletanalysis(int n,int ndat, float *psi_in, float *psi_out,int idim)
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

  unsigned int ShBaseElem = tid_hw + 2*NUM_LINES*par[idim].hwoffset_copy[hwid];

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
      psi_sh[ShBaseElem]=psi_in[BaseElem+ndat*2*npos];
      psi_sh[ShBaseElem+1]=psi_in[BaseElem+ndat*(2*npos+1)];
      ShBaseElem += 2*HALF_WARP_SIZE;
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

      //values of the convolution, even components
      register float conv = 
	//hand-unrolled loop 
	//order changed for increasing the precision
	WTFIL7 *psi_sh[ShBaseElem + 2*NUM_LINES] +
	WTFIL15*psi_sh[ShBaseElem +15*NUM_LINES] +
	WTFIL0 *psi_sh[ShBaseElem +16*NUM_LINES] +
	WTFIL8 *psi_sh[ShBaseElem + 1*NUM_LINES] +
	WTFIL6 *psi_sh[ShBaseElem + 4*NUM_LINES] +
	WTFIL14*psi_sh[ShBaseElem +13*NUM_LINES] +
	WTFIL1 *psi_sh[ShBaseElem +14*NUM_LINES] +
	WTFIL9 *psi_sh[ShBaseElem + 3*NUM_LINES] +
	WTFIL5 *psi_sh[ShBaseElem + 6*NUM_LINES] +
	WTFIL13*psi_sh[ShBaseElem +11*NUM_LINES] +
	WTFIL2 *psi_sh[ShBaseElem +12*NUM_LINES] +
	WTFIL10*psi_sh[ShBaseElem + 5*NUM_LINES] +
	WTFIL4 *psi_sh[ShBaseElem + 8*NUM_LINES] +
	WTFIL12*psi_sh[ShBaseElem + 9*NUM_LINES] +
	WTFIL3 *psi_sh[ShBaseElem +10*NUM_LINES] +
	WTFIL11*psi_sh[ShBaseElem + 7*NUM_LINES] ;

      psi_out[BaseElem]=conv;



      ///HERE WE SHOULD CHANGE ACCORDINGLY TO THE FILTERS
      //values of the convolution, odd components
      conv = 
	//hand-unrolled loop 
	//order changed for increasing the precision
	WTFIL15*psi_sh[ShBaseElem + 2*NUM_LINES] -
	WTFIL7 *psi_sh[ShBaseElem +15*NUM_LINES] +
	WTFIL8 *psi_sh[ShBaseElem +16*NUM_LINES] -
	WTFIL0 *psi_sh[ShBaseElem + 1*NUM_LINES] +
	WTFIL14*psi_sh[ShBaseElem + 4*NUM_LINES] -
	WTFIL6 *psi_sh[ShBaseElem +13*NUM_LINES] +
	WTFIL9 *psi_sh[ShBaseElem +14*NUM_LINES] -
	WTFIL1 *psi_sh[ShBaseElem + 3*NUM_LINES] +
	WTFIL13*psi_sh[ShBaseElem + 6*NUM_LINES] -
	WTFIL5 *psi_sh[ShBaseElem +11*NUM_LINES] +
	WTFIL10*psi_sh[ShBaseElem +12*NUM_LINES] -
	WTFIL2 *psi_sh[ShBaseElem + 5*NUM_LINES] +
	WTFIL12*psi_sh[ShBaseElem + 8*NUM_LINES] -
	WTFIL4 *psi_sh[ShBaseElem + 9*NUM_LINES] +
	WTFIL11*psi_sh[ShBaseElem +10*NUM_LINES] -
	WTFIL3 *psi_sh[ShBaseElem + 7*NUM_LINES] ;

      psi_out[BaseElem+n]=conv;

      //psi_sh[ShBaseElem+LOWFIL*par[idim].LinesPerBlock]; //for testing only

      ShBaseElem += HALF_WARP_SIZE;
      BaseElem += HW_ELEM;
      
    }

 
}

//MUST CONTINUE WITH #D ANALYSIS SYNTHESIS VERSION

extern "C" 
void localpotential_(int *n1,int *n2,int *n3,
		     float **psi,float **work,float **pot,
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


int magicfilterpot(int n1,int n2, int n3,
		   float *psi,
		   float *work)
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

/****/

