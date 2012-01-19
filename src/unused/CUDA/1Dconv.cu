/****u* CUDA/1Dconv.cu
**
** AUTHOR
**  Matthieu Ospici
**
** CHANGELOG
** Started on  Tue Mar 11 14:54:55 2008 Matthieu Ospici
** Last update Thu Sep 18 18:22:55 2008 Luigi Genovese
** Last update Thu Dec 25 2008 Thierry Deutsch: comment cutil.h to avoid the use of libcutil.a
**
** SOURCE
*/

#include <stdio.h>

//#include <cutil.h>
//#include <multithreading.h>
#include <pthread.h>
#include <semaphore.h>
#include <sched.h>

#include <time.h>
#include <sys/time.h>



#include "convolution.h"

//#include "newcase.cu"
//#include "conv_pot_shared_kernel.cu"
//#include "deffct.h" //convseria

#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)

#define SIZE_SHARED 288


__global__ void conv_lg(unsigned int n,unsigned int ndat,float *psi_in,float *psi_out,
			int lowfil,int lupfil)
{

  const unsigned int tid = threadIdx.x; //thread id
  const unsigned int bid = blockIdx.x;  //block id (can be 2-dim)

  const unsigned int numthreads = blockDim.x;


  //things to be done in this version of the kernel:
  //1)pass to constant memory the parameters, including the filters
  //2)implement texture fetching for non-coalesced addressed in read
  //memory
  //3)eliminate non-warp dependencies

  //all these elements can be decided with the block sizes
  //number of parts in which the line is separated
  // (n-1)/256+1
  int linecuts = ((n-1) >> 8)+1;

  //line cut in which the block is placed
  int icut = bid % linecuts; // to be controlled if modulo operator is
			     // not too slow

  //start point in n axis of the line treated by the block
  int istart = icut * (n/linecuts);

  //end point
  int iend = min((icut+1)*(n/linecuts),n);

  //number of elements
  int num_elem = iend-istart;

  //actual column treated by the block, in ndat axis
  int idat = bid / linecuts;

  //starting element in the global memory for the OUTPUT array
  int startelem = n*idat + istart;

  //distance from the halfwarp base element to guarantee coalesced
  //access in global memory
  int startd = startelem & 15;  // bitwise modulo 16
  
  //the first element which should be written in the output array
  //should be processed by thread startd

  __shared__ float psi_sh[SIZE_SHARED];  

  //in case the line cut is not the first, correct the starting point
  //such that the left elements are copied also

  //NOTE: it is assumed that for non-first segments the starting
  //points is far enough for the filter to be contained
  //and the same for non-last segments.
  //in other terms: lenght of the line is always bigger than max(lowfil,lupfil)

  int ileftcopystart,irightcopystart;
  int thcpyleft1,thcpyright1,thcpyright0,thcpyleft0;

  //start copy for the central points. the first thread copy the first
  //point since we are guaranteed there are enough threads
  //NOTE: thread 0 use the shared memory bank lowfil+startd & 15

  if (tid < num_elem)
    {
      psi_sh[tid+lowfil+startd]=psi_in[(ndat+1)*(istart+tid)+idat];
    }

  thcpyleft0 = 0; //temporary assignments, to be changed to avoid bank conflicts
  thcpyleft1 = lowfil;

  thcpyright0 = numthreads - lupfil;
  thcpyright1 = numthreads;

  //copy the rightmost points. if there are enough threads continue
  //with them, otherwise use threads of the same shared memory banks
  //than the initial ones
  irightcopystart = (icut < linecuts -1 ? istart+num_elem : 0) - thcpyright0;
  ileftcopystart = (icut > 0 ? istart-lowfil : n-lowfil) - thcpyleft0;
  
  if (tid >=thcpyleft0 && tid < thcpyleft1)
    {
      psi_sh[tid+startd-thcpyleft0]=
	psi_in[(ndat+1)*(ileftcopystart+tid)+idat];
    }

  if (tid >=thcpyright0 && tid < thcpyright1)
    {
      psi_sh[tid+startd+lowfil+num_elem-thcpyright0]=
	psi_in[(ndat+1)*(irightcopystart+tid)+idat];
    }


  //thread 0 starts writing at memory bank startd+1

  //there is now the copy of the central points; the first central
  //point will be copied by thread lowfil, then we continue until the 
  //minimum between numthreads and num_elem is reached
  //then recuperate the elements which are lacking
  //in case of need the starting thread of recuperation should access
  //the same bank as the first copy
  //bank associated to thcpycenter1 =(thcpycenter1 + startd) & 15: 
  //end shared memory copy
  __syncthreads();


  //copy the filters in shared memory for the moment
  //while they have to be put in constant memory space
 
  __shared__ float fil[16];

  fil[0] = 8.4334247333529341094733325815816e-7f;
  fil[1] =-0.1290557201342060969516786758559028e-4f;
  fil[2] = 0.8762984476210559564689161894116397e-4f;
  fil[3] =-0.30158038132690463167163703826169879e-3f;
  fil[4] = 0.174723713672993903449447812749852942e-2f;
  fil[5] =-0.942047030201080385922711540948195075e-2f;
  fil[6] = 0.2373821463724942397566389712597274535e-1f;
  fil[7] = 0.612625895831207982195380597e-1f;
  fil[8] = 0.9940415697834003993178616713f;
  fil[9] =-0.604895289196983516002834636e-1f;
  fil[10]=-0.2103025160930381434955489412839065067e-1f;
  fil[11]= 0.1337263414854794752733423467013220997e-1f;
  fil[12]=-0.344128144493493857280881509686821861e-2f;
  fil[13]= 0.49443227688689919192282259476750972e-3f;
  fil[14]=-0.5185986881173432922848639136911487e-4f;
  fil[15]= 2.72734492911979659657715313017228e-6f;


  //perform convolution in shared memory and write results in the
  //lowfil-scaled address
  register float conv = 0.f;
  #pragma unroll 20
  for(int j=0;j < lowfil+lupfil+1;++j)
    {
       conv += fil[j]*psi_sh[tid+j];
    }
  
  //write in global memory by taking care of the coalescing
  if (tid >= startd && tid < startd + num_elem)
    {
      psi_out[startelem-startd+tid]=conv;//psi_sh[tid+lowfil];
    }

}

__constant__ param_t param;

//__global__ void conv_shared(unsigned int n1,unsigned int n2,float *t_out,float *t_in,int nf);

__global__ void conv_shared2(unsigned int n1,unsigned int n2,float *t_out,float *t_in,int nf)
{

__shared__ float shared_temp[SIZE_SHARED_TOTAL];




  const unsigned int thid = threadIdx.x; //ID of current thread
  const unsigned int bidx = blockIdx.x;
  const unsigned int SIZE_SHARED_1 = param.SIZE_SHARED_1; 
  const unsigned int lineNumber = thid % 16;  


  int lineNumberFetch = lineNumber;

  //such conditionals are warp-independent and do not create problems
  if(bidx == gridDim.x - 1)
    {
      lineNumberFetch = (int)lineNumber - (int)param.lineLastBlock;

      if(lineNumberFetch <= 0)
	lineNumberFetch += (int)param.lineLastBlock;
      else
	lineNumberFetch = 0;
    }


  //copy into memory

  const unsigned int  step = (gridDim.y == 1) ? 0 : blockIdx.y/(gridDim.y - 1);
  const unsigned int TO_COPY =  param.fetchTabs.tabSeq[step][thid/16];    
  const int offset = param.fetchTabs.currPosTab[step][thid/16]; // columns offset (per thread)         
  const unsigned int TO_COPY_CALC =  param.calcTabs.tabSeq[step][thid/16];       
  const unsigned int offset_calc = param.calcTabs.currPosTab[step][thid/16]; // calc collums offset (per thread)


  int baseOffset;
      
  if(gridDim.y == 1)
    baseOffset = 0; //only one block
  else
    {
      baseOffset = param.sizeLineFirstBlocks * blockIdx.y;
    }
  
  
  
  for(int i=0,mod = baseOffset + offset + param.lowfil; i < TO_COPY; ++i)
    {
      if(mod >= param.SIZE_SHARED_2)
	{ 
	  mod = mod - param.SIZE_SHARED_2;
	}
      else if(mod < 0)
	{
	  mod = (param.SIZE_SHARED_2 + mod);
	}
      
      shared_temp[( i + offset)*16 + lineNumber] = t_in[(lineNumberFetch + SIZE_SHARED_1 *(bidx)) + mod*n1];
      
      
      ++mod;
    }
  
  
  __syncthreads();
  
  
  for(int i=0;i < TO_COPY_CALC;++i)
    {
      register float tmp = 0;
      
      #pragma unroll 20
      for(unsigned int j=0 ;j < nf;++j)
	{
	  
	  
	  tmp += shared_temp[lineNumber + (i + offset_calc + j)*16] * param.f_fct[j];

	}   
      t_out[(lineNumberFetch + SIZE_SHARED_1 *(bidx))*n2 +(i+offset_calc + baseOffset)] = tmp;
    }

}



__global__ void conv_shared_multi(unsigned int n1,unsigned int n2,float *t_out,float *t_in,int nf);

int dooldconv(int n1,
	       int n2, 
	       float *GPU_idata,
	       float *GPU_odata,
	       const float *filters,
	       int fsize,
	       int lowfil,
	       int lupfil,
	       unsigned int num_threads);

void donewconv(int n,
	       int ndat,
	       float *psi_in,
	       float *psi_out,
	       int lowfil,
	       int lupfil);



#include "convSeria.h"

//pthread_mutex_t mutexTime = PTHREAD_MUTEX_INITIALIZER; //for timers

//pthread_mutex_t mutexCond;
//pthread_cond_t pcond;

/*static sem_t sem;
static sem_t sem2;
static sem_t sem3;


static pthread_t *threads; //beark beark !!!*/


//********** INTERFACE BETWEEN FORTRAN & C ************


/*
n1,n2,n3 : fortran dimension
t_in : in array
t_out : out array
f_data : the f function
lowfil
lupfil
 */

extern "C" 
void previous1dconv_(int *n,
		     int *ndat,
		      float *t_in, 
		      float *t_out,
		      float *f_data,
		      int *lowfil, 
		      int *lupfil)
{

  const int fsize = *lupfil - *lowfil + 1;

  multiTab_t m_dataIn;
  multiTab_t m_dataOut;
  


  m_dataIn.n1 =  1;
  m_dataIn.n2 =  *ndat;
  m_dataIn.n3 = *n + 1;

  m_dataOut = m_dataIn;
  
  m_dataIn.tab = t_in;
  m_dataOut.tab = t_out;



  unsigned int mem_size1 = m_dataIn.n1 * m_dataIn.n2 * m_dataIn.n3 * sizeof(float);
  
  unsigned int mem_size2 = m_dataOut.n1 * m_dataOut.n2 * m_dataOut.n3 * sizeof(float);


  if(cudaMalloc( (void**) &(m_dataIn.GPU_data), mem_size1) != 0)
    {
      printf("erreur malloc GPU_odata\n");
      printf(" %i \n",mem_size1);
      return;
    }
 
  if( cudaMalloc( (void**) &(m_dataOut.GPU_data), mem_size2) != 0) //only the data
    {
      printf("erreur malloc GPU_idata\n");        
      return;
    } 



  evalPerfGPU_t evPerf;
  



  conv1dGPU(&m_dataIn,
	    &m_dataOut,
	    f_data,
	    fsize,
	    *lowfil,
	    *lupfil,
	    256,
	    &evPerf);
  
}


extern "C"
void n1dconv_(int *nx,
	      int *nconv,
	      float **data_in, 
	      float **data_out, 
	      float *filters, 
	      int *lowfil, 
	      int *lupfil)
{
  //define the number of threads and blocks according to parameter
  //definitions
  const int ndat = *nconv -1;
  const int n = *nx+1;

  donewconv(n,ndat,*data_in,*data_out,*lowfil,*lupfil);

  return;
}
     
void donewconv(int n,
	       int ndat,
	       float *psi_in,
	       float *psi_out,
	       int lowfil,
	       int lupfil)
{
  //declare the texture for binding the input psi
  //texture<float, 2> psi_tex;

  const int num_threads = min(64* ((n+15)/64 + 1),256);


  //printf(" %i numthdandblock %i %i\n",num_threads,ndat+1); 
  //printf(" %i numthds %i \n",num_threads,(n-1) >> 8);

  //cudaBindTexture(NULL,psi_tex,psi_in,n*ndat*sizeof(float));

  dim3 grid1((ndat+1)*(n/num_threads), 1, 1);  
  dim3 threads1(num_threads, 1, 1);

  //launch the kernel grid
  conv_lg <<< grid1, threads1 >>>(n,ndat, psi_in, psi_out,lowfil, lupfil);

  //cudaUnbindTexture(psi_tex);

  cudaThreadSynchronize();

  return;  

}

//interface, only the 1d convolution
extern "C" 
void m1dconv_(int *n, 
	      int *ndat, 
	      float **data_in, 
	      float **data_out, 
	      float *filters, 
	      int *lowfil, 
	      int *lupfil)
{
  const int fsize = *lupfil - *lowfil + 1;

  const int n1 = *ndat;
  const int n2 = *n+1;

  
  if(dooldconv(n1,
	       n2, 
	       *data_in,
	       *data_out,
	       filters,
	       fsize,
	       *lowfil,
	       *lupfil,
	       256) != 0)
    {
      return;
    } 
  return; 
}

int dooldconv(int n1,
	      int n2, 
	      float *GPU_idata,
	      float *GPU_odata,
	      const float *filters,
	      int fsize,
	      int lowfil,
	      int lupfil,
	      unsigned int num_threads)
{

  //create the parameters
  param_t paramToSend;
  unsigned int numBlockDim1,numBlockDim2;

  createParam(&paramToSend,num_threads,n1,n2,&numBlockDim1,&numBlockDim2,
	      filters,fsize,lowfil,lupfil);

  //printf("numBlocks1 %i, num_blocks2 %i\n",numBlockDim1,numBlockDim2);

  //send them to constant memory
  
  if(cudaMemcpyToSymbol(param,&paramToSend, sizeof(param_t)) != 0)
    {

      printf("MemcpyToSymbol error\n");

      return 1;
    }
 
  //printf(" %i numthdandblock %i %i\n",num_threads,numBlockDim1,numBlockDim2); 

  //define the number of threads and blocks according to parameter definitions

  dim3  grid1(numBlockDim1, numBlockDim2, 1);  
  dim3  threads1(num_threads, 1, 1);
 
  //launch the kernel grid
  conv_shared2 <<< grid1, threads1 >>>(n1,n2,GPU_odata, GPU_idata,fsize);

  cudaThreadSynchronize();

  return 0;

}


//********** END INTERFACE BETWEEN FORTRAN & C ************



int conv1dGPU(multiTab_t* m_dataIn,
	      multiTab_t* m_dataOut,
	      const float *f_data,
	      int fsize,
	      int lowfil,
	      int lupfil,
	      unsigned int num_threads,
	      evalPerfGPU_t *evPerf)

{

  //******* copy CPU -> GPU
 
 
  const unsigned int mem_size = m_dataIn->n1 *m_dataIn->n2 * m_dataIn->n3 * sizeof(float) ;

 

  // unsigned int timer;
  float* res_data = m_dataOut->tab;


  //evPerf->host_malloc = cutGetTimerValueM(timer);

  //allocation
  float* GPU_idata = m_dataIn->GPU_data;  
  // printf("GPU_idata %p\n",GPU_idata);
  float* GPU_odata = m_dataOut->GPU_data;  

 
  evPerf->cuda_malloc =0;

  //copy

  if(cudaMemcpy( GPU_idata, m_dataIn->tab, mem_size, cudaMemcpyHostToDevice)  != 0)
    {
      printf("erreur memcpy GPU_idata\n");
      cudaFree(GPU_odata);
      cudaFree(GPU_idata);
      return 1;
    }


  //evPerf->GPU_trsf_CPUGPU = cutGetTimerValueM(timer);


  param_t paramToSend;
  unsigned int numBlockDim1,numBlockDim2;
  
  //************ Code execution

 // ----------- STEP 1 --------



    unsigned int n1 = m_dataIn->n1*m_dataIn->n2;
    unsigned int n2 = m_dataIn->n3; 


  createParam(&paramToSend,num_threads,n1,n2,&numBlockDim1,&numBlockDim2,f_data,fsize,lowfil,lupfil);

 
  if( cudaMemcpyToSymbol(param,&paramToSend, sizeof(param_t)) != 0)
    {

      printf("Erreur copy to symbol param\n");

      cudaFree(GPU_odata);
      cudaFree(GPU_idata);
      return 1;
      }
  
  dim3  grid1(numBlockDim1, numBlockDim2, 1);  
  dim3  threads1(num_threads, 1, 1);


 
  //  cutStartTimer(timer);
 
   conv_shared_multi <<< grid1, threads1 >>>(n1,n2,GPU_odata, GPU_idata,fsize);
    cudaThreadSynchronize();
 
    //   cutStopTimer(timer);
     // ----------- END STEP 1 --------


    //evPerf->GPU_calc = cutGetTimerValueM(timer);

  //************* COPY back to CPU
 

  if( cudaMemcpy( res_data, 
		  GPU_odata,
		  mem_size,
		  cudaMemcpyDeviceToHost) != 0)
    {
      cudaFree(GPU_odata);
      cudaFree(GPU_idata);
      printf("Erreur memcpy devtohost\n");
      return 1;
    }


  //evPerf->GPU_trsf_GPUCPU = cutGetTimerValueM(timer);
  


  cudaError_t err = cudaGetLastError();   
  if( cudaSuccess != err)
    {                                                
      printf("Cuda error:  in file '%s' in line %i : %s.\n",
	      __FILE__, 
	      __LINE__, 
	      cudaGetErrorString( err) );

      cudaFree(GPU_odata);
      cudaFree(GPU_idata);
      return 1;              

    }


 
  cudaFree(GPU_odata);
  cudaFree(GPU_idata);
  return 0;
}


