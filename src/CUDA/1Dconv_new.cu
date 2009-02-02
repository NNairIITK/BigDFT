/****u* CUDA/1Dconv_new.cu
**
** AUTHOR
**  Luigi Genovese
**
** SOURCE
*/

#include <stdio.h>
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

//error handling in CUDA
#define CUERR { cudaError_t err; \
 if ((err = cudaGetLastError()) != cudaSuccess) { \
 printf("CUDA error: %s, line %d\n", cudaGetErrorString(err), __LINE__); }}


typedef struct  _par
{
  //unsigned int LinesPerBlock;
  unsigned int ElementsPerBlock;
  //unsigned int ElementsPerHalfWarp;

  int lowfil, lupfil; //structure of f_fct
  int thline[HALF_WARP_SIZE]; //line considered by a thread within the half-warp
  int thelem[HALF_WARP_SIZE]; //elements considered by a thread within the half-warp
  int hwelem_calc[16]; //maximum number of half warps
  int hwelem_copy[16]; //maximum number of half-warps
  int hwoffset_calc[16]; //maximum number of half warps
  int hwoffset_copy[16]; //maximum number of half-warps
  
  //float fil[16];

} par_t;

__constant__ par_t par;

//declare the texture for binding the input psi
texture<float,2, cudaReadModeElementType> psi_tex;

#include "kernel_tex.cu"

int dogenconv(int ndat,
	      int n, 
	      float *GPU_idata,
	      float *GPU_odata,
	      int lowfil,
	      int lupfil);

int dogenconv_n(int ndat,
	      int n, 
	      float *GPU_idata,
	      float *GPU_odata,
	      int lowfil,
	      int lupfil);

void correctSequence(int thds,int elem,int * tab);

//create the parameters to be used for calculating the convolution
//with a given stride
void constantParameters(par_t* par,
			unsigned int* num_halfwarps,
			//unsigned int num_lines,
			int n,
			int ndat,
			int lowfil, //lowfil + lupfil must be a multiple of 16
			int lupfil,
			unsigned int* linecuts,
			unsigned int* num_blocks)

{

  //number of lines treated by each block
  //par->LinesPerBlock = num_lines;

  //number of total allowed elements of a input line
  unsigned int num_elem_tot=MAX_SHARED_SIZE/sizeof(float)/NUM_LINES; //between1024and64
  //unsigned int num_elem_tot=MAX_SHARED_SIZE/sizeof(float)/num_lines; //between 1024 and 64
  
  //number of elements of the output
  unsigned int num_elem_max=min(num_elem_tot-lowfil-lupfil-1,n); //between 1008 and 48 for 16-fil

  //number of elements treated by the single half-warp
  //par -> ElementsPerHalfWarp = HALF_WARP_SIZE/num_lines; //it is assumed they are multiples

  //number of pieces in which a line is divided
  //if the line is too small and not a multiple of ElementsPerHalfWarp
  //divide the line in two
  *linecuts=
    //(n <= num_elem_max && n % par -> ElementsPerHalfWarp !=0 ? 2 : (n-1)/num_elem_max+1);
    (n <= num_elem_max && n % HW_ELEM !=0 ? 2 : (n-1)/num_elem_max+1);

  //number of blocks in ndat direction
  //*num_blocks=((ndat-1)/num_lines + 1);
  *num_blocks=((ndat-1)/NUM_LINES + 1);

  //printf("num_elem_tot %i,num_elem_max %i,linecuts %i,num_blocks %i,elemperHW %i \n",
  //num_elem_tot,num_elem_max,*linecuts,*num_blocks, par -> ElementsPerHalfWarp);

  //number of elements treated by each block 
  //this may pose problems for values of n dimensions less than 48
  //when n is not a multiple of ElementsPerHalfWarp
  par->ElementsPerBlock = 
    //min(par->ElementsPerHalfWarp*(((n-1)/(*linecuts))/par->ElementsPerHalfWarp+1),n);
    min(HW_ELEM*(((n-1)/(*linecuts))/HW_ELEM+1),n);

  int halfwarps=16;
  //calculate the maximum number of halfwarps (between 4 and 16)
  for(int i =3; i>=0; --i)
    {
      //if(par->ElementsPerBlock/par->ElementsPerHalfWarp >= 1 << i)
      if(par->ElementsPerBlock/HW_ELEM >= 1 << i)
	{
	  halfwarps = 1 << i;
	  break;
	}
    }

  *num_halfwarps = halfwarps;

  for(int j=0;j < HALF_WARP_SIZE ; ++j)
    {
      //par->thline[j]= j & (num_lines - 1); //num_lines always a power of two 
      //par->thelem[j]= j / num_lines; 

      par->thline[j]= j & (NUM_LINES - 1); //num_lines always a power of two 
      par->thelem[j]= j / NUM_LINES; 
    }

  //define the sequences of the number of elements
  correctSequence(halfwarps,par->ElementsPerBlock/HW_ELEM,par->hwelem_calc);

  correctSequence(halfwarps,(par->ElementsPerBlock+lowfil+lupfil+1)/HW_ELEM,
		  par->hwelem_copy);

  //correctSequence(halfwarps,par->ElementsPerBlock/par->ElementsPerHalfWarp,
  //par->hwelem_calc);

  //correctSequence(halfwarps,(par->ElementsPerBlock+lowfil+lupfil+1)/par->ElementsPerHalfWarp,
  //par->hwelem_copy);


  //define the offsets
  for(int j=0,pos_calc=0,pos_copy=0;j < halfwarps ; ++j)
    {
      par->hwoffset_calc[j]=pos_calc;
      par->hwoffset_copy[j]=pos_copy;
      pos_calc+=HW_ELEM*par->hwelem_calc[j];
      pos_copy+=HW_ELEM*par->hwelem_copy[j];
      //pos_calc+=par->ElementsPerHalfWarp*par->hwelem_calc[j];
      //pos_copy+=par->ElementsPerHalfWarp*par->hwelem_copy[j];

    }
 
  //lowfil and lupfil parameters
  par->lowfil = lowfil;
  par->lupfil = lupfil;

  //printf("ElementsPerBlock %i,HalfWarpCalculatedElements %i,HalfWarpCopiedElements %i,LastHalfWarpCalcElements %i, LastHalfWarpCopiedElements %i \n",
  //par->ElementsPerBlock,par->hwelem_calc[0],par->hwelem_copy[0],
  //par->hwelem_calc[halfwarps-1],par->hwelem_copy[halfwarps-1]);

  /*
   //filter values for this convolution, hard coded
  par->fil[0] = 8.4334247333529341094733325815816e-7f;
  par->fil[1] =-0.1290557201342060969516786758559028e-4f;
  par->fil[2] = 0.8762984476210559564689161894116397e-4f;
  par->fil[3] =-0.30158038132690463167163703826169879e-3f;
  par->fil[4] = 0.174723713672993903449447812749852942e-2f;
  par->fil[5] =-0.942047030201080385922711540948195075e-2f;
  par->fil[6] = 0.2373821463724942397566389712597274535e-1f;
  par->fil[7] = 0.612625895831207982195380597e-1f;
  par->fil[8] = 0.9940415697834003993178616713f;
  par->fil[9] =-0.604895289196983516002834636e-1f;
  par->fil[10]=-0.2103025160930381434955489412839065067e-1f;
  par->fil[11]= 0.1337263414854794752733423467013220997e-1f;
  par->fil[12]=-0.344128144493493857280881509686821861e-2f;
  par->fil[13]= 0.49443227688689919192282259476750972e-3f;
  par->fil[14]=-0.5185986881173432922848639136911487e-4f;
  par->fil[15]= 2.72734492911979659657715313017228e-6f;
  */

}

/* //already defined in locpot
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
*/


//1D convolution of multiple lines in the same block
__global__ void conv1d_stride(int n,int ndat, float *psi_in, float *psi_out)
{

  //line treated by the given block
  unsigned int lineOffset = min(blockIdx.y*NUM_LINES,ndat-NUM_LINES);
  //starting element treated by the block
  unsigned int elemOffset = min(blockIdx.x*par.ElementsPerBlock,n-par.ElementsPerBlock);

  //line treated by the given block
  //unsigned int lineOffset = min(blockIdx.y*par.LinesPerBlock,ndat-par.LinesPerBlock);
  //starting element treated by the block
  //unsigned int elemOffset = min(blockIdx.x*par.ElementsPerBlock,n-par.ElementsPerBlock);

  //half-warp id
  const unsigned int hwid = threadIdx.y;
  //tid within the HW
  const unsigned int tid_hw = threadIdx.x;

  //shared memory array
  __shared__ float psi_sh[MAX_SHARED_SIZE/sizeof(float)];

  //line treated by the given thread in ndat axis
  //which is the input base element
  unsigned int BaseElem = par.thline[tid_hw] + lineOffset;
  //write data in shared memory
  //element treated by the given thread in n-axis
  unsigned int thelem = par.thelem[tid_hw] + par.hwoffset_copy[hwid];

  unsigned int ShBaseElem = tid_hw + NUM_LINES*par.hwoffset_copy[hwid];
  //unsigned int ShBaseElem = tid_hw + par.LinesPerBlock*par.hwoffset_copy[hwid];

  int epsilon,npos;

  //NOTE: it is assumed that for non-first segments the starting
  //points is far enough for the filter to be contained
  //and the same for non-last segments.
  //in other terms: lenght of the line is always bigger than
  //max(lowfil,lupfil)

  for(int i=0,ipos=elemOffset-par.lowfil+thelem;i < par.hwelem_copy[hwid] ; ++i)
    {
      epsilon=(ipos < 0 ? -1 : ipos/n);
      npos=ipos-epsilon*n;
      psi_sh[ShBaseElem]=psi_in[BaseElem+ndat*npos];
      //psi_sh[ShBaseElem]=tex1Dfetch(psi_tex,BaseElem+ndat*npos);

      ShBaseElem += HALF_WARP_SIZE;
      ipos += HW_ELEM;
      //ipos += par.ElementsPerHalfWarp;
      
    }

  //end shared memory copy
  __syncthreads();

  //element treated by the given thread in n-axis
  thelem = par.thelem[tid_hw] + par.hwoffset_calc[hwid];
  //base element for the given thread in shared memory
  ShBaseElem = tid_hw + NUM_LINES*par.hwoffset_calc[hwid];
  //ShBaseElem = tid_hw + par.LinesPerBlock*par.hwoffset_calc[hwid];

  //output base element, from the input one
  BaseElem =  n*BaseElem+ thelem + elemOffset;

  //perform convolution in shared memory 
  //each thread calculate a number of elements, identical for each
  //half-warp
  //#pragma unroll 5 (to be tested if it is important)

  for(int i=0;i < par.hwelem_calc[hwid]; ++i)
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

	/*
	par.fil[0]*psi_sh[ShBaseElem               ] +
	par.fil[15]*psi_sh[ShBaseElem + 15*NUM_LINES] +
	par.fil[1]*psi_sh[ShBaseElem +   NUM_LINES ] +
	par.fil[14]*psi_sh[ShBaseElem + 14*NUM_LINES] +
	par.fil[2]*psi_sh[ShBaseElem + 2*NUM_LINES ] +
	par.fil[13]*psi_sh[ShBaseElem + 13*NUM_LINES] +
	par.fil[3]*psi_sh[ShBaseElem + 3*NUM_LINES ] +
	par.fil[12]*psi_sh[ShBaseElem + 12*NUM_LINES] +
	par.fil[4]*psi_sh[ShBaseElem + 4*NUM_LINES ] +
	par.fil[11]*psi_sh[ShBaseElem + 11*NUM_LINES] +
	par.fil[5]*psi_sh[ShBaseElem + 5*NUM_LINES ] +
	par.fil[10]*psi_sh[ShBaseElem + 10*NUM_LINES] +
	par.fil[6]*psi_sh[ShBaseElem + 6*NUM_LINES ] +
	par.fil[9]*psi_sh[ShBaseElem + 9*NUM_LINES ] +
	par.fil[7]*psi_sh[ShBaseElem + 7*NUM_LINES ] +
	par.fil[8]*psi_sh[ShBaseElem + 8*NUM_LINES ] ;

      
	par.fil[0]*psi_sh[ShBaseElem               ] +
	par.fil[15]*psi_sh[ShBaseElem + 15*par.LinesPerBlock] +
	par.fil[1]*psi_sh[ShBaseElem +   par.LinesPerBlock ] +
	par.fil[14]*psi_sh[ShBaseElem + 14*par.LinesPerBlock] +
	par.fil[2]*psi_sh[ShBaseElem + 2*par.LinesPerBlock ] +
	par.fil[13]*psi_sh[ShBaseElem + 13*par.LinesPerBlock] +
	par.fil[3]*psi_sh[ShBaseElem + 3*par.LinesPerBlock ] +
	par.fil[12]*psi_sh[ShBaseElem + 12*par.LinesPerBlock] +
	par.fil[4]*psi_sh[ShBaseElem + 4*par.LinesPerBlock ] +
	par.fil[11]*psi_sh[ShBaseElem + 11*par.LinesPerBlock] +
	par.fil[5]*psi_sh[ShBaseElem + 5*par.LinesPerBlock ] +
	par.fil[10]*psi_sh[ShBaseElem + 10*par.LinesPerBlock] +
	par.fil[6]*psi_sh[ShBaseElem + 6*par.LinesPerBlock ] +
	par.fil[9]*psi_sh[ShBaseElem + 9*par.LinesPerBlock ] +
	par.fil[7]*psi_sh[ShBaseElem + 7*par.LinesPerBlock ] +
	par.fil[8]*psi_sh[ShBaseElem + 8*par.LinesPerBlock ] ;


	par.fil[0]*psi_sh[ShBaseElem               ] +
	par.fil[1]*psi_sh[ShBaseElem +   par.LinesPerBlock ] +
	par.fil[2]*psi_sh[ShBaseElem + 2*par.LinesPerBlock ] +
	par.fil[3]*psi_sh[ShBaseElem + 3*par.LinesPerBlock ] +
	par.fil[4]*psi_sh[ShBaseElem + 4*par.LinesPerBlock ] +
	par.fil[5]*psi_sh[ShBaseElem + 5*par.LinesPerBlock ] +
	par.fil[6]*psi_sh[ShBaseElem + 6*par.LinesPerBlock ] +
	par.fil[7]*psi_sh[ShBaseElem + 7*par.LinesPerBlock ] +
	par.fil[8]*psi_sh[ShBaseElem + 8*par.LinesPerBlock ] +
	par.fil[9]*psi_sh[ShBaseElem + 9*par.LinesPerBlock ] +
	par.fil[10]*psi_sh[ShBaseElem + 10*par.LinesPerBlock] +
	par.fil[11]*psi_sh[ShBaseElem + 11*par.LinesPerBlock] +
	par.fil[12]*psi_sh[ShBaseElem + 12*par.LinesPerBlock] +
	par.fil[13]*psi_sh[ShBaseElem + 13*par.LinesPerBlock] +
	par.fil[14]*psi_sh[ShBaseElem + 14*par.LinesPerBlock] +
	par.fil[15]*psi_sh[ShBaseElem + 15*par.LinesPerBlock];
      */

      psi_out[BaseElem]=conv;
      //psi_sh[ShBaseElem+par.lowfil*par.LinesPerBlock]; //for testing only

      ShBaseElem += HALF_WARP_SIZE;
      BaseElem += HW_ELEM;
      //BaseElem += par.ElementsPerHalfWarp;

      
    }

 
}


//assume psi_in is a 2d texture
__global__ void conv1d_stride_tex(int n,int ndat,float *psi_out) 
{

  //line treated by the given block
  unsigned int lineOffset = min(blockIdx.y*NUM_LINES,ndat-NUM_LINES);
  //starting element treated by the block
  unsigned int elemOffset = min(blockIdx.x*par.ElementsPerBlock,n-par.ElementsPerBlock);
  //half-warp id
  const unsigned int hwid = threadIdx.y;
  //tid within the HW
  const unsigned int tid_hw = threadIdx.x;

  //shared memory array
  __shared__ float psi_sh[MAX_SHARED_SIZE/sizeof(float)];

  //line treated by the given thread in ndat axis
  //which is the input base element
  unsigned int BaseElem = par.thline[tid_hw] + lineOffset;
  //write data in shared memory
  //element treated by the given thread in n-axis
  //const unsigned int thelem = par.thelem[tid_hw] + par.hwoffset_copy[hwid]
  //unsigned int thelem2 = (unsigned int) par.thelem[tid_hw] + par.hwoffset_copy[hwid];
  unsigned int thelem2 = par.hwoffset_copy[hwid];

  unsigned int ShBaseElem = tid_hw + NUM_LINES*par.hwoffset_copy[hwid];

  //int epsilon,npos;
  //float x,y;

  //NOTE: it is assumed that for non-first segments the starting
  //points is far enough for the filter to be contained
  //and the same for non-last segments.
  //in other terms: lenght of the line is always bigger than
  //max(lowfil,lupfil)

  for(int i=0,ipos=elemOffset-par.lowfil+thelem2;i < par.hwelem_copy[hwid] ; ++i)
    {
      //epsilon=(ipos < 0 ? -1 : ipos/n);
      //npos=ipos-epsilon*n;
      //psi_sh[ShBaseElem]=psi_in[BaseElem+ndat*npos];
      //x=(float) (npos)/;
      //y=(float)npos +0.5f; 
      //x=((float) (BaseElem) + 0.5f )/((float) (ndat));
      //y=((float) (npos) + 0.5f )/((float) (n));
      psi_sh[ShBaseElem]=12.f;//tex2D(psi_tex,x,y);
      //CUERR;
      ShBaseElem += HALF_WARP_SIZE;
      ipos += HW_ELEM;
      
    }

  //end shared memory copy
  __syncthreads();

  //element treated by the given thread in n-axis
  unsigned int thelem = par.thelem[tid_hw] + par.hwoffset_calc[hwid];
  //base element for the given thread in shared memory
  ShBaseElem = tid_hw + NUM_LINES*par.hwoffset_calc[hwid];
  //ShBaseElem = tid_hw + par.LinesPerBlock*par.hwoffset_calc[hwid];

  //output base element, from the input one
  BaseElem =  n*BaseElem + thelem + elemOffset;

  //perform convolution in shared memory 
  //each thread calculate a number of elements, identical for each
  //half-warp
  //#pragma unroll 5 (to be tested if it is important)

  for(int i=0;i < par.hwelem_calc[hwid]; ++i)
    {
    //values of the convolution
    /* register float conv = 
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
	MFIL8 *psi_sh[ShBaseElem + 8*NUM_LINES ] ; */

    psi_out[BaseElem]=//=conv;
	psi_sh[ShBaseElem+par.lowfil*NUM_LINES]; //for testing only

    ShBaseElem += HALF_WARP_SIZE;
    BaseElem += HW_ELEM;
    //BaseElem += par.ElementsPerHalfWarp;

      
    }

 
}


//interface, only the 1d convolution
extern "C" 
void g1dconv_(int *n, 
	      int *ndat, 
	      float **data_in, 
	      float **data_out, 
	      float *filters, 
	      int *lowfil, 
	      int *lupfil)
{

  const int n1 = *ndat;
  const int n2 = *n+1;

  
  if(dogenconv_n(n1,
	       n2, 
	       *data_in,
	       *data_out,
	       *lowfil,
	       *lupfil) != 0)
    {
      return;
    } 
  return; 
}


int dogenconv(int ndat,
	      int n, 
	      float *GPU_idata,
	      float *GPU_odata,
	      int lowfil,
	      int lupfil)
{

  //create the parameters
  par_t parCPU;
  //cudaChannelFormatDesc channelDesc = 
  //cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

  cudaArray* psiCA;

  //allocate the cuda array for the Psi wavefunction
  //create the channel format descriptor
  //cudaChannelFormatDesc channelDesc=
  // cudaCreateChannelDesc<float>();

  if(cudaMallocArray(&psiCA,&psi_tex.channelDesc,ndat,n) != 0)
    {
      printf("GPU allocation error \n");
      return 1;
    }


  //copy the input array to the cuda Array
  if(cudaMemcpyToArray(psiCA,0,0,GPU_idata, ndat*n*sizeof(float),cudaMemcpyDeviceToDevice)  != 0)
    {
      printf("DeviceToDevice Memcpy error \n");
      return 1;
    }


  //calculate the number of threads and blocks
  //unsigned int num_lines = min(16,ndat); //hard coded for the moment
  unsigned int numBlocks,linecuts,num_halfwarps;
  //int tex_offset;
  //size_t offset;

  constantParameters(&parCPU,&num_halfwarps,n,ndat,lowfil,lupfil,
		     &linecuts,&numBlocks);
  //constantParameters(&parCPU,&num_halfwarps,num_lines,n,ndat,lowfil,lupfil,
  //&linecuts,&numBlocks);



  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i\n",linecuts,numBlocks,num_halfwarps);

  //send them to constant memory
  if(cudaMemcpyToSymbol(par,&parCPU, sizeof(par_t)) != 0)
    {
      printf("MemcpyToSymbol error\n");

      return 1;
    }
 
  //define the number of threads and blocks according to parameter definitions
  dim3  grid1(linecuts,  numBlocks, 1);  
  dim3  threads1(HALF_WARP_SIZE, num_halfwarps , 1);

  //set texture parameters
  psi_tex.addressMode[0] = cudaAddressModeWrap;
  psi_tex.addressMode[1] = cudaAddressModeWrap;
  psi_tex.filterMode = cudaFilterModePoint;
  psi_tex.normalized = 1;


  //bind the texture reference to the CUDA array
  cudaBindTextureToArray(psi_tex,psiCA,psi_tex.channelDesc);
  CUERR;

  //element offset for reading from the texture
  //tex_offset = offset/sizeof(float);
  
  //printf(" offset %i\n",tex_offset); 
  //launch the kernel grid
  //conv1d_stride_tex <<< grid1, threads1 >>>(n,ndat, GPU_odata);
  //conv1d_stride <<< grid1, threads1 >>>(n,ndat, GPU_idata, GPU_odata);

  conv1d_stride_tex_m<<< grid1, threads1 >>>(n,ndat, GPU_idata, GPU_odata);
  cudaThreadSynchronize();

  //unbind the texture
  cudaUnbindTexture(psi_tex);
  CUERR;
  //free the CUDA Array
  //CUDA_SAFE_CALL(cudaFreeArray(psiCA));;
  cudaFreeArray(psiCA);;

  //cudaThreadSynchronize();

  return 0;

}

int dogenconv_n(int ndat,
	      int n, 
	      float *GPU_idata,
	      float *GPU_odata,
	      int lowfil,
	      int lupfil)
{

  //create the parameters
  par_t parCPU;

  //calculate the number of threads and blocks
  unsigned int numBlocks,linecuts,num_halfwarps;

  constantParameters(&parCPU,&num_halfwarps,n,ndat,lowfil,lupfil,
		     &linecuts,&numBlocks);

  //send them to constant memory
  if(cudaMemcpyToSymbol(par,&parCPU, sizeof(par_t)) != 0)
    {
      printf("MemcpyToSymbol error\n");

      return 1;
    }
 
  //define the number of threads and blocks according to parameter definitions
  dim3  grid1(linecuts,  numBlocks, 1);  
  dim3  threads1(HALF_WARP_SIZE, num_halfwarps , 1);

  //launch the kernel grid
  conv1d_stride <<< grid1, threads1 >>>(n,ndat, GPU_idata, GPU_odata);
  cudaThreadSynchronize();
  
  return 0;

}

/****/
