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
#define MAX_SHARED_SIZE 3840 //16*240 =~ 3.9 kB
#define HALF_WARP_SIZE 16 // for all architectures

typedef struct   _par
{
  unsigned int LinesPerBlock;
  unsigned int ElementsPerBlock;
  unsigned int HalfWarpCalculatedElements;
  unsigned int HalfWarpCopiedElements;
  unsigned int LastHalfWarpCopiedElements;

  int lowfil, lupfil; //structure of f_fct
  float fil[16];

} par_t;

__constant__ par_t par;

int dogenconv(int ndat,
	      int n, 
	      float *GPU_idata,
	      float *GPU_odata,
	      int lowfil,
	      int lupfil);

//create the parameters to be used for calculating the convolution
//with a given stride
void constantParameters(par_t* par,
			unsigned int num_threads,
			unsigned int num_lines,
			int n,
			int ndat,
			int lowfil,
			int lupfil,
			unsigned int* num_blocks)
{

  //number of lines treated by each block
  par->LinesPerBlock = num_lines;

  //number of total allowed elements of a input line
  unsigned int num_elem_tot=MAX_SHARED_SIZE/sizeof(float)/num_lines; //between 960 and 60
  
  //number of elements of the output
  unsigned int num_elem_max=min(num_elem_tot-lowfil-lupfil-1,n); //between 944 and 44 for 16 fil

  //number of pieces in which a line is divided
  unsigned int linecuts=(n-1)/num_elem_max+1;

  //number of half warps treating the single line cut
  unsigned int linehws=num_threads/HALF_WARP_SIZE/num_lines; //it is assumed they are multiples

  //total number of blocks we must have
  *num_blocks=((ndat-1)/num_lines + 1)*linecuts;

  //number of elements treated by each block 
  par->ElementsPerBlock = linehws*((n/linecuts-1)/linehws+1);

  //number of elements which should be calculated by each thread in a half-warp
  par->HalfWarpCalculatedElements = par -> ElementsPerBlock/linehws; //multiples

  //number of elements which should be copied by each thread in a half-warp
  par->HalfWarpCopiedElements = (par -> ElementsPerBlock+lowfil+lupfil)/linehws + 1;

  //the last half-warp copies less elements in general
  par->LastHalfWarpCopiedElements = par -> ElementsPerBlock+lowfil+lupfil+1
    - par -> HalfWarpCopiedElements;
  
  //lowfil and lupfil parameters
  par->lowfil = lowfil;
  par->lupfil = lupfil;

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


}

//1D convolution of multiple lines by the same block
__global__ void conv1d_stride(unsigned int n,unsigned int ndat,float *psi_in,float *psi_out)
{
  
  const unsigned int num_lines = par.LinesPerBlock;
  const unsigned int num_elem = par.ElementsPerBlock;
  const unsigned int hwcopy_gen = par.HalfWarpCopiedElements;
  const unsigned int hwcalc_gen = par.HalfWarpCalculatedElements;
  const unsigned int lowfil = par.lowfil;
  const unsigned int lupfil = par.lupfil;
  const unsigned int linecuts = (n-1)/num_elem + 1;
  //number of elements calculated by a half-warp
  const unsigned int hwelems = HALF_WARP_SIZE/num_lines; 

  const unsigned int tid = threadIdx.x; //thread id
  const unsigned int bid = blockIdx.x+gridDim.x*blockIdx.y; // 2D for more than 2^16 lines
  const unsigned int numhw = blockDim.x/HALF_WARP_SIZE;

  //line treated by the given block
  unsigned int lineOffset = min((bid/linecuts)*num_lines,ndat-num_lines);
  //starting element treated by the block
  unsigned int elemOffset = min((bid & linecuts)*num_elem,n-num_elem);
  //line treated by the given thread
  unsigned int thline = tid % num_lines; //always HWS % numlines =0
  //half-warp id
  const unsigned int hwid = tid / HALF_WARP_SIZE; // =tid >> 4 if HWS=16
  //tid within the HW
  const unsigned int tid_hw = tid & (HALF_WARP_SIZE -1); //HWS is a power of
						  //two
  const unsigned int hwcopy_eff = (hwid == numhw-1?par.LastHalfWarpCopiedElements:hwcopy_gen);

  __shared__ float psi_sh[MAX_SHARED_SIZE];  

  //write data in shared memory
  //element treated by the given thread in n-axis
  unsigned int thelem = tid_hw/num_lines + hwid * hwelems * hwcopy_gen;
  //base element for the given thread in shared memory
  unsigned int ShBaseElem = thline + num_lines*thelem;
  //alternative formulation
  //unsigned int ShBaseElem = thline + tid_hw + hwid *HALF_WARP_SIZE * hwcopy_gen;

  unsigned int InBaseElem = thline+lineOffset+(thelem + elemOffset)*ndat;

  int epsilon,npos;

  //NOTE: it is assumed that for non-first segments the starting
  //points is far enough for the filter to be contained
  //and the same for non-last segments.
  //in other terms: lenght of the line is always bigger than max(lowfil,lupfil)
  for(int i=0,ipos=elemOffset-lowfil+tid_hw+hwid*hwcopy_gen;i < hwcopy_eff; ++i)
    {
      epsilon=(ipos < 0 ? -1 : ipos/n);
      npos=ipos-epsilon*n;
      psi_sh[ShBaseElem+i*HALF_WARP_SIZE]=psi_in[InBaseElem+ndat*(hwelems*i+npos-elemOffset)];
      ipos += hwelems;
    }


  //end shared memory copy
  __syncthreads();

  //element treated by the given thread in n-axis
  thelem = tid_hw/num_lines + hwid *hwelems * hwcalc_gen;
  //base element for the given thread in shared memory
  ShBaseElem = thline + num_lines*thelem;

  unsigned int OutBaseElem = n*(thline+lineOffset)+ thelem + elemOffset;

  //perform convolution in shared memory 
  //each thread calculate a number of elements, identical for each half-warp
  for(int i=0;i < hwcalc_gen; ++i)
    {
      //values of the convolution
      register float conv = 0.f;
      #pragma unroll 20 //loop unrolling should be performed by hand
      for(int j=0;j < lowfil+lupfil+1;++j)
	{
	  conv += 
	    par.fil[j]*psi_sh[ShBaseElem + j*num_lines +i*HALF_WARP_SIZE];
	}
      psi_out[OutBaseElem+i*hwelems]=psi_sh[ShBaseElem+lowfil*num_lines]; //for testing only
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

  
  if(dogenconv(n1,
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

  //calculate the number of threads and blocks
  unsigned int num_threads = 256;
  unsigned int num_lines = min(16,ndat);
  unsigned int numBlocks;

  constantParameters(&parCPU,num_threads,num_lines,n,ndat,lowfil,lupfil,&numBlocks);

  unsigned int num_blocksx=numBlocks & 65535;
  unsigned int num_blocksy=numBlocks >> 16; // this should be
					    // corrected to allow
					    // prime values

  //send them to constant memory
  if(cudaMemcpyToSymbol(par,&parCPU, sizeof(par_t)) != 0)
    {
      printf("MemcpyToSymbol error\n");

      return 1;
    }
 
  //printf(" %i numthdandblock %i %i\n",num_threads,numBlockDim1,numBlockDim2); 

  //define the number of threads and blocks according to parameter definitions

  dim3  grid1(num_blocksx, num_blocksy, 1);  
  dim3  threads1(num_threads, 1, 1);
 
  //launch the kernel grid
  conv1d_stride <<< grid1, threads1 >>>(n,ndat,GPU_odata, GPU_idata);

  cudaThreadSynchronize();

  return 0;

}

/****/
