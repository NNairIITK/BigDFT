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
#define MAX_SHARED_SIZE 4096 //16*256 4 kB (should be =~ 3.9 kB, try also 3072)
#define HALF_WARP_SIZE 16 // for all architectures

typedef struct   _par
{
  unsigned int LinesPerBlock;
  unsigned int ElementsPerBlock;
  unsigned int ElementsPerHalfWarp;
  //unsigned int HalfWarpCalculatedElements;
  //unsigned int HalfWarpCopiedElements;
  //unsigned int LastHalfWarpCopiedElements;

  int lowfil, lupfil; //structure of f_fct
  int thline[HALF_WARP_SIZE]; //line considered within the half-warp
  int hwelem_calc[16]; //maximum number of half warps
  int hwelem_copy[16]; //maximum number of half-warps
  int hwoffset_calc[16]; //maximum number of half warps
  int hwoffset_copy[16]; //maximum number of half-warps
  
  float fil[16];

} par_t;

__constant__ par_t par;

int dogenconv(int ndat,
	      int n, 
	      float *GPU_idata,
	      float *GPU_odata,
	      int lowfil,
	      int lupfil);

void correctSequence(int thds,int elem,int * tab);

//create the parameters to be used for calculating the convolution
//with a given stride
void constantParameters(par_t* par,
			unsigned int num_halfwarps,
			unsigned int num_lines,
			int n,
			int ndat,
			int lowfil,
			int lupfil,
			unsigned int* linecuts,
			unsigned int* num_blocks)

{

  //number of lines treated by each block
  par->LinesPerBlock = num_lines;

  //number of total allowed elements of a input line
  unsigned int num_elem_tot=MAX_SHARED_SIZE/sizeof(float)/num_lines; //between 1024 and 64
  
  //number of elements of the output
  unsigned int num_elem_max=min(num_elem_tot-lowfil-lupfil-1,n); //between 1008 and 48 for 16-fil

  //number of pieces in which a line is divided
  *linecuts=(n-1)/num_elem_max+1;

  //total number of blocks we must have
  *num_blocks=((ndat-1)/num_lines + 1)* (*linecuts);

  
  //number of elements treated by the single half-warp
  par -> ElementsPerHalfWarp = HALF_WARP_SIZE/num_lines; //it is assumed they are multiples

  /*
  //minimum number of elements treated by a single block
  //this should not be bigger than num_elem_max
  unsigned int num_elem_min=(num_halfwarps*HALF_WARP_SIZE)/num_lines;
  */

  //printf("num_elem_tot %i,num_elem_max %i,linecuts %i,hwelems %i,num_blocks %i, num_elem_min %i \n",
  //num_elem_tot,num_elem_max,*linecuts,hwelems,*num_blocks,num_elem_min);

  //number of elements treated by each block 
  par->ElementsPerBlock = min(HALF_WARP_SIZE*(((n-1)/(*linecuts))/HALF_WARP_SIZE+1),n);

  for(int j=0;j < HALF_WARP_SIZE ; ++j)
    {
      par->thline[j]= j & (num_lines - 1); //num_lines always a power of two       
    }

  //define the sequences of the number of elements
  correctSequence(num_halfwarps,par->ElementsPerBlock,par->hwelem_calc);

  correctSequence(num_halfwarps,par->ElementsPerBlock+lowfil+lupfil+1,par->hwelem_copy);

  //define the offsets
  for(int j=0,pos_calc=0,pos_copy=0;j < num_halfwarps ; ++j)
    {
      par->hwoffset_calc[j]=pos_calc;
      par->hwoffset_copy[j]=pos_copy;
      pos_calc+=par->hwelem_calc[j];
      pos_copy+=par->hwelem_copy[j];
    }
 

  /*
  //number of elements which should be calculated by each thread in a half-warp
  par->HalfWarpCalculatedElements = par -> ElementsPerBlock/num_elem_min; //multiples
  
  
  //the last half-warp copies less elements in general
  par->LastHalfWarpCopiedElements = 
  */
  
  //lowfil and lupfil parameters
  par->lowfil = lowfil;
  par->lupfil = lupfil;

  //printf("ElementsPerBlock %i,HalfWarpCalculatedElements %i,HalfWarpCopiedElements %i,LastHalfWarpCalcElements %i, LastHalfWarpCopiedElements %i \n",
  //par->ElementsPerBlock,par->hwelem_calc[0],par->hwelem_copy[0],
  //par->hwelem_calc[num_halfwarps-1],par->hwelem_copy[num_halfwarps-1]);

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

//uniformise the tabular of the number of elements treated by each
//thread (counterpart of uniformiseTab)
//divide the operations in three steps: the general, particular, last
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



//1D convolution of multiple lines by the same block
__global__ void conv1d_stride(int n,int ndat,float *psi_in,float *psi_out)
{
  

  //lines treated by each block, constant values
  const unsigned int num_lines = par.LinesPerBlock;
  const unsigned int num_elem = par.ElementsPerBlock;
  //const unsigned int hwcopy_gen = par.hwelem_copy[0];//HalfWarpCopiedElements;
  //const unsigned int hwcalc_gen = par.hwelem_calc[0];//HalfWarpCalculatedElements;
  const unsigned int lowfil = par.lowfil;
  //number of output elements of a single half-warp
  const unsigned int hwelems = par.ElementsPerHalfWarp; 

  //line cut id
  const unsigned int linecutid = blockIdx.x;
  //block id (in ndat axis)
  const unsigned int bid = blockIdx.y;
  //line treated by the given block
  unsigned int lineOffset = min(bid*num_lines,ndat-num_lines);
  //starting element treated by the block
  unsigned int elemOffset = min(linecutid*num_elem,n-num_elem);

  //half-warp id
  const unsigned int hwid = threadIdx.y;
  //tid within the HW
  const unsigned int tid_hw = threadIdx.x;

  __shared__ float psi_sh[MAX_SHARED_SIZE/sizeof(float)];

  //line treated by the given thread
  unsigned int thline = par.thline[tid_hw];
  //write data in shared memory
  //element treated by the given thread in n-axis
  unsigned int thelem = tid_hw/num_lines + par.hwoffset_copy[hwid];//hwid * hwelems * hwcopy_gen;
  //base element for the given thread in shared memory
  //unsigned int ShBaseElem = thline + num_lines*thelem;
  //alternative formulation
  //unsigned int ShBaseElem = tid_hw + hwid *HALF_WARP_SIZE * hwcopy_gen;
  //alternative
  unsigned int ShBaseElem = tid_hw + num_lines*par.hwoffset_copy[hwid];

  unsigned int InBaseElem = thline+lineOffset;//+(thelem + elemOffset)*ndat;

  unsigned int hwcopy_eff = par.hwelem_copy[hwid];
  //(hwid == blockDim.y-1?par.LastHalfWarpCopiedElements:hwcopy_gen);
  
  int epsilon,npos;

  //NOTE: it is assumed that for non-first segments the starting
  //points is far enough for the filter to be contained
  //and the same for non-last segments.
  //in other terms: lenght of the line is always bigger than
  //max(lowfil,lupfil)

  for(int i=0,ipos=elemOffset-lowfil+thelem;i < hwcopy_eff ; ++i)
    {
      epsilon=(ipos < 0 ? -1 : ipos/n);
      npos=ipos-epsilon*n;
      psi_sh[ShBaseElem+i*HALF_WARP_SIZE]=psi_in[InBaseElem+ndat*npos];
      //psi_sh[thelem+i*hwelems][thline]=psi_in[InBaseElem+ndat*npos];
      ipos += hwelems;
    }


  //end shared memory copy
  __syncthreads();

  //element treated by the given thread in n-axis
  thelem = tid_hw/num_lines + par.hwoffset_calc[hwid];//hwid *hwelems * hwcalc_gen;
  //base element for the given thread in shared memory
  //ShBaseElem = thline + num_lines*thelem;
  ShBaseElem = tid_hw + num_lines*par.hwoffset_calc[hwid];//hwid *HALF_WARP_SIZE * hwcalc_gen;

  unsigned int hwcalc_eff = par.hwelem_calc[hwid];

  unsigned int OutBaseElem = n*(thline+lineOffset)+ thelem + elemOffset;

  //perform convolution in shared memory 
  //each thread calculate a number of elements, identical for each
  //half-warp
  //#pragma unroll 5 (to be tested if it is important)
  for(int i=0;i < hwcalc_eff; ++i)
    {
      //values of the convolution
      register float conv =
      //hand-unrolled loop (16 elements for this filter)
	
	par.fil[0]*psi_sh[ShBaseElem               +i*HALF_WARP_SIZE] +
	par.fil[1]*psi_sh[ShBaseElem +   num_lines +i*HALF_WARP_SIZE] +
	par.fil[2]*psi_sh[ShBaseElem + 2*num_lines +i*HALF_WARP_SIZE] +
	par.fil[3]*psi_sh[ShBaseElem + 3*num_lines +i*HALF_WARP_SIZE] +
	par.fil[4]*psi_sh[ShBaseElem + 4*num_lines +i*HALF_WARP_SIZE] +
	par.fil[5]*psi_sh[ShBaseElem + 5*num_lines +i*HALF_WARP_SIZE] +
	par.fil[6]*psi_sh[ShBaseElem + 6*num_lines +i*HALF_WARP_SIZE] +
	par.fil[7]*psi_sh[ShBaseElem + 7*num_lines +i*HALF_WARP_SIZE] +
	par.fil[8]*psi_sh[ShBaseElem + 8*num_lines +i*HALF_WARP_SIZE] +
	par.fil[9]*psi_sh[ShBaseElem + 9*num_lines +i*HALF_WARP_SIZE] +
	par.fil[10]*psi_sh[ShBaseElem + 10*num_lines +i*HALF_WARP_SIZE] +
	par.fil[11]*psi_sh[ShBaseElem + 11*num_lines +i*HALF_WARP_SIZE] +
	par.fil[12]*psi_sh[ShBaseElem + 12*num_lines +i*HALF_WARP_SIZE] +
	par.fil[13]*psi_sh[ShBaseElem + 13*num_lines +i*HALF_WARP_SIZE] +
	par.fil[14]*psi_sh[ShBaseElem + 14*num_lines +i*HALF_WARP_SIZE] +
	par.fil[15]*psi_sh[ShBaseElem + 15*num_lines +i*HALF_WARP_SIZE];

	/*
	par.fil[0]*psi_sh[thelem +    i*hwelems][thline]  +
	par.fil[1]*psi_sh[thelem + 1 +i*hwelems][thline]  +
	par.fil[2]*psi_sh[thelem + 2 +i*hwelems][thline]  +
	par.fil[3]*psi_sh[thelem + 3 +i*hwelems][thline]  +
	par.fil[4]*psi_sh[thelem + 4 +i*hwelems][thline]  +
	par.fil[5]*psi_sh[thelem + 5 +i*hwelems][thline]  +
	par.fil[6]*psi_sh[thelem + 6 +i*hwelems][thline]  +
	par.fil[7]*psi_sh[thelem + 7 +i*hwelems][thline]  +
	par.fil[8]*psi_sh[thelem + 8 +i*hwelems][thline]  +
	par.fil[9]*psi_sh[thelem + 9 +i*hwelems][thline]  +
	par.fil[10]*psi_sh[thelem + 10 +i*hwelems][thline] +
	par.fil[11]*psi_sh[thelem + 11 +i*hwelems][thline] +
	par.fil[12]*psi_sh[thelem + 12 +i*hwelems][thline] +
	par.fil[13]*psi_sh[thelem + 13 +i*hwelems][thline] +
	par.fil[14]*psi_sh[thelem + 14 +i*hwelems][thline] +
	par.fil[15]*psi_sh[thelem + 15 +i*hwelems][thline] ;
	*/
      
      psi_out[OutBaseElem+i*hwelems]=conv;
      //psi_sh[ShBaseElem+(lowfil)*num_lines+i*HALF_WARP_SIZE]; //for testing only
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
  unsigned int num_halfwarps = 16;
  unsigned int num_lines = min(16,ndat); //hard coded for the moment
  unsigned int numBlocks,linecuts;

  constantParameters(&parCPU,num_halfwarps,num_lines,n,ndat,lowfil,lupfil,
		     &linecuts,&numBlocks);

  //printf("num_blocksx %i, num_blocksy %i\n",linecuts,numBlocks);

  //send them to constant memory
  if(cudaMemcpyToSymbol(par,&parCPU, sizeof(par_t)) != 0)
    {
      printf("MemcpyToSymbol error\n");

      return 1;
    }
 
  //printf(" %i numthdandblock %i %i\n",num_threads,numBlockDim1,numBlockDim2); 

  //define the number of threads and blocks according to parameter definitions

  dim3  grid1(linecuts, numBlocks, 1);  
  dim3  threads1(HALF_WARP_SIZE, num_halfwarps , 1);
 
  //launch the kernel grid
  conv1d_stride <<< grid1, threads1 >>>(n,ndat,GPU_idata, GPU_odata);

  cudaThreadSynchronize();

  return 0;

}
