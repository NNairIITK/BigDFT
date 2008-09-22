#include <stdio.h>
#include <cutil.h>
#include <pthread.h>
#include <semaphore.h>
#include <sched.h>

#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)

typedef struct   _par
{
  unsigned int SIZE_SHARED_1 ;
  int SIZE_SHARED_2;
 

  fetchTabs_t fetchTabs;
  calcTabs_t calcTabs;

  
  unsigned  int sizeLineFirstBlocks;
  unsigned  int sizeLineLastBlock;


  int lowfil, lupfil; //structure of f_fct
  float f_fct[100];

  unsigned int lineLastBlock;


 

} par_t;

//create the parameters to be used for calculating the convolution
//with a given stride
void constantParameters(param_t* param,
		 unsigned int numThread,
		 unsigned int num_elements_dim1,
		 unsigned int num_elements_dim2,
		 unsigned int *numBlockDim1,
		 unsigned int *numBlockDim2,
		 const float *f_data,
		 int fsize,
		 int lowfil, 
		 int lupfil)
{

  //maximum size of the shared memory array
  //conceived for maximize occupancy on a hardware of compute
  //capability 1.2 and higher (1024 threads at same time on a given multiprocessor)
  const int MAX_SHARED_SIZE=3840; //16*240 =~ 3.9 kB

  const int HALF_WARP_SIZE=warpSize/2; //16 for all architectures

  //number of total elements of a input line
  unsigned int num_elem_tot=MAX_SHARED_SIZE/sizeof(float)/num_lines; //between 960 and 60
  
  //number of elements of the output
  unsigned int num_elem=min(num_elem_tot-lowfil-lupfil-1,n); //between 944 and 44

  //number of pieces in which a line is divided
  unsigned int linecuts=(n-1)/num_elem+1;
   
  //number of threads treating the single line
  unsigned int linethds=numthds/num_lines; //it is assumed they are multiples

  //number of elements calculated by the single thread (halfwarp-dependent value)
  unsigned int numcalc=num_elem/linethds;

  //filter values for this convolution
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



  // const int MAX_LINE_SIZE = 100;//255;
  param->lowfil = lowfil;
  param->lupfil = lupfil;


  param->SIZE_SHARED_2 = num_elements_dim2;
  param->SIZE_SHARED_1 = 16;//SIZE_SHARED_TOTAL/param->SIZE_SHARED_2 ;


 
  /* if(param->SIZE_SHARED_1 > 16)
     param->SIZE_SHARED_1 = 16;*/
 
  *numBlockDim1 = ceilf((float)num_elements_dim1 / param->SIZE_SHARED_1);



  // -------- line > 255 --------

  //we need 2 tabSeq & currPosTab. One for the firsts x block, and one
  //for the last
 
  int nbBlockDim2 = ceilf((float)num_elements_dim2/(MAX_LINE_SIZE - fsize));
  // param->size_shared_2_real = param->SIZE_SHARED_2 + fsize * nbBlockDim2;
  // printf("nbBlockDim2 %i\n",nbBlockDim2);

  unsigned int sizeLineFirstBlocks;
  unsigned int sizeLineLastBlock;

  if(nbBlockDim2 > 1) //a line is on more of one block
    {
      // sizeLineLastBlock = num_elements_dim2 - (MAX_LINE_SIZE - fsize)*(nbBlockDim2 - 1);
      //	    sizeLineFirstBlocks = MAX_LINE_SIZE - fsize;
      sizeLineFirstBlocks = ceilf((float)num_elements_dim2/nbBlockDim2);
      sizeLineLastBlock = num_elements_dim2 - (nbBlockDim2 - 1)*sizeLineFirstBlocks;
      
      
      //Creation of the first tab (for first blocks)
      //fetch tab (we fetch nf in more)
      createTabs(param->fetchTabs.tabSeq,
		 param->fetchTabs.currPosTab,
		 numThread,
		 sizeLineFirstBlocks + fsize,
		 0);
      
      createTabs(param->fetchTabs.tabSeq,
		 param->fetchTabs.currPosTab,
		 numThread,
		 sizeLineLastBlock + fsize,
		 1);
      
      //calc tabs
      createTabs(param->calcTabs.tabSeq,
		 param->calcTabs.currPosTab,
		 numThread,
		 sizeLineFirstBlocks ,
		 0);

      createTabs(param->calcTabs.tabSeq,
		 param->calcTabs.currPosTab,
		 numThread,
		 sizeLineLastBlock,
		 1);

      
      // printf("sizeLineLastBlock %i, sizeLineFirstBlocks %i\n",sizeLineLastBlock,sizeLineFirstBlocks);
    }
  else
    {
      createTabs(param->fetchTabs.tabSeq,
		 param->fetchTabs.currPosTab,
		 numThread,
		 num_elements_dim2 + fsize,
		 0);

      createTabs(param->calcTabs.tabSeq,
		 param->calcTabs.currPosTab,
		 numThread,
		 num_elements_dim2,
		 0);

      sizeLineFirstBlocks = num_elements_dim2;
      sizeLineLastBlock = sizeLineFirstBlocks;
      
    }

  *numBlockDim2 = nbBlockDim2;
  param->sizeLineFirstBlocks = sizeLineFirstBlocks;
  param->sizeLineLastBlock = sizeLineLastBlock;

  // -------- end line > 255 ----

  //if last block is not full (lineNumber%16 !=0)
  
  param->lineLastBlock = num_elements_dim1  % 16 ;
  if( param->lineLastBlock == 0)
    param->lineLastBlock = 16;
  else
    --(param->lineLastBlock);
  //---- end last block not full

  for(int i=0;i<fsize;++i)
    {
      param->f_fct[i] = f_data[i];
    }


  // const int tailleTab = (2*numThread)/32; //(2*numThread)/32;

}

//1D convolution of multiple lines by the same block
__global__ void conv1d_stride(unsigned int n,unsigned int ndat,float *psi_in,float *psi_out,
			int lowfil,int lupfil)
{
  
  const unsigned int num_lines = par.LinesPerBlock;
  const unsigned int num_elem = par.ElementsPerBlock;
  const unsigned int hwcopy_gen = par.HalfWarpCopiedElements;
  const unsigned int hwcalc_gen = par.HalfWarpCalculatedElements;
  const unsigned int lowfil = par.lowfil;
  const unsigned int linecuts = (n-1)/num_elem + 1;
  //number of elements calculated by a half-warp
  const unsigned int hwelems = HALF_WARP_SIZE/num_lines; 

  const unsigned int tid = threadIdx.x; //thread id
  const unsigned int bid = blockIdx.x;  //block id (can be 2-dim)

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

  __shared__ float psi_sh[SIZE_SHARED_MAX];  

  //write data in shared memory
  //element treated by the given thread in n-axis
  unsigned int thelem = tid_hw/num_lines + hwid * hwelems * hwcopy_gen;
  //base element for the given thread in shared memory
  unsigned int ShBaseElem = thline + num_lines*thelem;
  //alternative formulation
  //unsigned int ShBaseElem = thline + tid_hw + hwid *HALF_WARP_SIZE * hwcopy_gen;

  unsigned int InBaseElem = thline+lineOffset+(thelem + elemOffset)*ndat;

  //NOTE: it is assumed that for non-first segments the starting
  //points is far enough for the filter to be contained
  //and the same for non-last segments.
  //in other terms: lenght of the line is always bigger than max(lowfil,lupfil)
  for(int i=0,ipos=elemOffset-lowfil+tid_hw;i < hwcopy_eff; ++i)
    {
      epsilon=(ipos < 0 ? -1 : ipos/n);
      npos=ipos-epsilon*n;
      psi_sh[ShBaseElem+i*HALF_WARP_SIZE]=psi_in[InBaseElem+ndat*(hwelems*i+npos-elemOffset)];
      ipos += hwelems;
    }


  //end shared memory copy
  __syncthreads();

  //element treated by the given thread in n-axis
  thelem = tid_hw/num_lines + hwid *hwlines * hwcalc_gen;
  //base element for the given thread in shared memory
  ShBaseElem = thline + num_lines*thelem;

  unsigned int OutBaseElem = n*(thline+lineOffset)+ thelem + elemOffset;


  //perform convolution in shared memory 
  //each thread calculate a number of elements, identical for each half-warp
  for(i=0;i < hwcalc_eff; ++i)
    {
      //values of the convolution
      register float conv = 0.f;
      #pragma unroll 20 //loop unrolling should be performed by hand
      for(int j=0;j < lowfil+lupfil+1;++j)
	{
	  conv += 
	    par.fil[j]*psi_sh[ShBaseElem + j*num_lines +i*HALF_WARP_SIZE];
	}
      psi_out[OutBaseElem+i*hwlines]=psi_sh[ShBaseElem+lowfil*num_lines]; //for testing only
    }
}
