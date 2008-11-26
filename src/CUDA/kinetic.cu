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
#define LOWFIL 14
#define LUPFIL 14


//convolution filters
#define KFIL0   -3.5536922899131901941296809374f
#define KFIL1    2.2191465938911163898794546405f
#define KFIL2   -0.6156141465570069496314853949f
#define KFIL3    0.2371780582153805636239247476f
#define KFIL4   -0.0822663999742123340987663521f
#define KFIL5    0.02207029188482255523789911295638968409f
#define KFIL6   -0.409765689342633823899327051188315485e-2f
#define KFIL7    0.45167920287502235349480037639758496e-3f
#define KFIL8   -0.2398228524507599670405555359023135e-4f
#define KFIL9    2.0904234952920365957922889447361e-6f
#define KFIL10  -3.7230763047369275848791496973044e-7f
#define KFIL11  -1.05857055496741470373494132287e-8f
#define KFIL12  -5.813879830282540547959250667e-11f
#define KFIL13   2.70800493626319438269856689037647576e-13f
#define KFIL14  -6.924474940639200152025730585882e-18f


/*
#define KFIL0   0.e-3f 
#define KFIL1   1.e-3f 
#define KFIL2   2.e-3f 
#define KFIL3   3.e-3f 
#define KFIL4   4.e-3f 
#define KFIL5   5.e-3f 
#define KFIL6   6.e-3f 
#define KFIL7   7.e-3f 
#define KFIL8   8.e-3f 
#define KFIL9   9.e-3f 
#define KFIL10 10.e-3f
#define KFIL11 11.e-3f
#define KFIL12 12.e-3f
#define KFIL13 13.e-3f
#define KFIL14 14.e-3f
*/

typedef struct  _parK
{
  unsigned int ElementsPerBlock;

  int thline[HALF_WARP_SIZE]; //line considered by a thread within the half-warp
  int thelem[HALF_WARP_SIZE]; //elements considered by a thread within the half-warp
  int hwelem_calc[16]; //maximum number of half warps
  int hwelem_copy[16]; //maximum number of half-warps
  int hwoffset_calc[16]; //maximum number of half warps
  int hwoffset_copy[16]; //maximum number of half-warps
  float scale;

} parK_t;

__constant__ parK_t par[3];

//declare the texture for binding the input psi
//texture<float> psi_tex;

void correctSequence(int thds,int elem,int * tab);

int kineticfilter(int n1,int n2, int n3,
		  float h1,float h2,float h3,
		  float *psi,
		  float *work1,
		  float *work2,
		  float *epot);

// Magic Filter parameters to be used for calculating the convolution
void KParameters(parK_t* par,
		 unsigned int* num_halfwarps,
		 int n,
		 int ndat,
		 float hgrid,
		 unsigned int* linecuts,
		 unsigned int* num_blocks)

{

  //number of total allowed elements of a input line
  unsigned int num_elem_tot=MAX_SHARED_SIZE/sizeof(float)/NUM_LINES; //between1024and64
  
  //number of elements of the output
  unsigned int num_elem_max=min(num_elem_tot-LOWFIL-LUPFIL-1,n); //between 996 and 35 

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

  //define the scale factor to be applied to the convolution
  par->scale=0.5f/(hgrid*hgrid);
  
}

//1D convolution of multiple lines in the same block
__global__ void kinetic1d(int n,int ndat, 
			  float *psi_in, float *psi_out, float *psi_in_rot, 
			  int idim)
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

  for(int i=0,ipos=elemOffset-LOWFIL+thelem;i < par[idim].hwelem_copy[hwid] ; ++i)
    {
      epsilon=(ipos < 0 ? -1 : ipos/n);
      npos=ipos-epsilon*n;
      psi_sh[ShBaseElem]=psi_in[BaseElem+ndat*npos];

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
	//hand-unrolled loop 
	//order changed for increasing the precision
	KFIL14*(psi_sh[ShBaseElem               ] +
		psi_sh[ShBaseElem + 28*NUM_LINES]) +
	KFIL13*(psi_sh[ShBaseElem +    NUM_LINES]  +
		psi_sh[ShBaseElem + 27*NUM_LINES]) +
	KFIL12*(psi_sh[ShBaseElem +  2*NUM_LINES]  +
		psi_sh[ShBaseElem + 26*NUM_LINES]) +
	KFIL11*(psi_sh[ShBaseElem +  3*NUM_LINES]  +
		psi_sh[ShBaseElem + 25*NUM_LINES]) +
	KFIL10*(psi_sh[ShBaseElem +  4*NUM_LINES]  +
		psi_sh[ShBaseElem + 24*NUM_LINES]) +
	KFIL9 *(psi_sh[ShBaseElem +  5*NUM_LINES]  +
		psi_sh[ShBaseElem + 23*NUM_LINES]) +
	KFIL8 *(psi_sh[ShBaseElem +  6*NUM_LINES]  +
		psi_sh[ShBaseElem + 22*NUM_LINES]) +
	KFIL7 *(psi_sh[ShBaseElem +  7*NUM_LINES]  +
		psi_sh[ShBaseElem + 21*NUM_LINES]) +
	KFIL6 *(psi_sh[ShBaseElem +  8*NUM_LINES]  +
		psi_sh[ShBaseElem + 20*NUM_LINES]) +
	KFIL5 *(psi_sh[ShBaseElem +  9*NUM_LINES]  +
		psi_sh[ShBaseElem + 19*NUM_LINES]) +
	KFIL4 *(psi_sh[ShBaseElem + 10*NUM_LINES]  +
		psi_sh[ShBaseElem + 18*NUM_LINES]) +
	KFIL3 *(psi_sh[ShBaseElem + 11*NUM_LINES]  +
		psi_sh[ShBaseElem + 17*NUM_LINES]) +
	KFIL2 *(psi_sh[ShBaseElem + 12*NUM_LINES]  +
		psi_sh[ShBaseElem + 16*NUM_LINES]) +
	KFIL1 *(psi_sh[ShBaseElem + 13*NUM_LINES]  +
		psi_sh[ShBaseElem + 15*NUM_LINES]) +
	KFIL0 * psi_sh[ShBaseElem + 14*NUM_LINES];

      psi_out[BaseElem]=-par[idim].scale*conv;//another array should
					      //be added
      //psi_sh[ShBaseElem+LOWFIL*NUM_LINES]; //for testing only
      psi_in_rot[BaseElem]=psi_sh[ShBaseElem+LOWFIL*NUM_LINES]; 

      ShBaseElem += HALF_WARP_SIZE;
      BaseElem += HW_ELEM;
      

    }
}


extern "C" 
void kineticterm_(int *n1,int *n2,int *n3,
		  float *hx,float *hy,float *hz,
		  float **psi,float **work1,float **work2,
		  float *ekin) 

{

  
  if(kineticfilter(*n1+1,*n2+1,*n3+1,
		   *hx,*hy,*hz,
		   *psi,*work1,*work2,
		   ekin) != 0)
    {
      printf("ERROR: GPU kineticfilter\n ");
      return;
    } 
  return; 
}


int kineticfilter(int n1,int n2, int n3,
		   float h1,float h2,float h3,
		   float *psi,
		   float *work1,
		   float *work2,
		   float *epot)
{

  //create the parameters
  parK_t parCPU[3];

  //calculate the number of threads and blocks
  unsigned int numBlocks,linecuts,num_halfwarps;

  //calculate the parameters in constant memory for each of the 1D convolution
  //define the number of threads and blocks according to parameter definitions
  KParameters(&parCPU[2],&num_halfwarps,n3,n1*n2,h3,&linecuts,&numBlocks);
  dim3  grid3(linecuts,  numBlocks, 1);  
  dim3  threads3(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n3,n1*n2);

  KParameters(&parCPU[1],&num_halfwarps,n2,n1*n3,h2,&linecuts,&numBlocks);
  dim3  grid2(linecuts,  numBlocks, 1);  
  dim3  threads2(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n2,n1*n3);

  KParameters(&parCPU[0],&num_halfwarps,n1,n2*n3,h1,&linecuts,&numBlocks);
  dim3  grid1(linecuts,  numBlocks, 1);  
  dim3  threads1(HALF_WARP_SIZE, num_halfwarps , 1);

  //printf("num_blocksx %i, num_blocksy %i, halfwarps %i,n1,ndat, %i %i\n",
  //linecuts,numBlocks,num_halfwarps,n1,n3*n2);

  //send them to constant memory, once and for all
  if(cudaMemcpyToSymbol(*par,&parCPU, 3*sizeof(parK_t)) != 0)
    {
      printf("MemcpyToSymbol error\n");

      return 1;
    }

  kinetic1d <<< grid3, threads3 >>>(n3,n1*n2,psi,work1,work2,2);
  cudaThreadSynchronize();

  //kinetic1d <<< grid2, threads2 >>>(n2,n1*n3,work2,psi,1);
  //cudaThreadSynchronize();

  //kinetic1d <<< grid1, threads1 >>>(n1,n2*n3,psi,pot,work,0);
  //cudaThreadSynchronize();

  return 0;

}

/****/
