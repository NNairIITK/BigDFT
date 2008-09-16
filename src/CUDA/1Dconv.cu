/*
** convolution.cu
** 
** Made by (Matthieu Ospici)
** Login   <ospici@badiane>
** 
** Started on  Tue Mar 11 14:54:55 2008 Matthieu Ospici
** Last update Fri Jun 27 17:12:28 2008 Luigi Genovese
*/

#include <stdio.h>
#include <cutil.h>
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

  //number of parts in which the line is separated
  int linecuts = 1;//((n-1)/numthreads + 1);

  //line cut in which the block is placed
  int icut = 0;//bid % linecuts; // to be controlled if modulo operator is
			     // not too slow

  //start point in n axis of the line treated by the block
  int istart = 0;//icut * (n/linecuts);

  //end point
  int iend = n;//min((icut+1)*(n/linecuts),n);

  //number of elements
  int num_elem = n;//iend-istart;

  //actual column treated by the block, in ndat axis
  int idat = bid;// / linecuts;

  //starting element in the global memory for the OUTPUT array
  int startelem = n*idat + istart;

  //distance from the halfwarp base element to guarantee coalesced
  //access in global memory
  int startd = 0;//startelem & 15;  // bitwise modulo 16
  
  //the first element which should be written in the output array
  //should be processed by thread startd

  __shared__ float psi_sh[SIZE_SHARED];  

  //in case the line cut is not the first, correct the starting point
  //such that the left elements are copied also
  //NOTE: it is assumed that for non-first segments the starting
  //points is far enough for the filter to be contained

  //int startdcorr=(icut > 0 ? lowfil : 0); 
  
  //copy of elements in shared memory
  if (tid < num_elem)
    {
      //bank conflicts for psi_sh: for each half-warp there is 
      //a linear addressing with stride 1, so ok.

      //coalesced accesses for psi_in: the access are completely
      //uncoalesced, must pass through texture fetches

      psi_sh[tid+lowfil+startd]=psi_in[(ndat+1)*(istart+tid)+idat];
    }


  //end shared memory copy
  __syncthreads();


  //perform convolution in shared memory and write results in the
  //lowfil-scaled address
  for(int j=0;j < lowfil+lupfil+1;++j)
    {
       register float tmp = 0;
       tmp = psi_sh[tid];
    }


  //write in global memory by taking care of the coalescing
  if (tid >= startd && tid < startd + num_elem)
    {
      psi_out[startelem-startd+tid]=psi_sh[tid+lowfil];
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

  //such conditionals should be avoided
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
  const int offset = param.fetchTabs.currPosTab[step][thid/16]; // collums offset (per thread)         
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
  const int num_threads = min(64* (n/64 + 1),256);


  printf(" %i numthds %i \n",num_threads,14 % 1);

  dim3 grid1(ndat+1, 1, 1);  
  dim3 threads1(num_threads, 1, 1);

  //launch the kernel grid
  conv_lg <<< grid1, threads1 >>>(n,ndat, psi_in, psi_out,lowfil, lupfil);

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

  //send them to constant memory
  
  if(cudaMemcpyToSymbol(param,&paramToSend, sizeof(param_t)) != 0)
    {

      printf("MemcpyToSymbol error\n");

      return 1;
    }
  
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


extern "C" 
void gpu_allocate__(int *nsize, //memory size
		    float **GPU_pointer, // pointer indicating the GPU address
		    int *ierr) // error code, 1 if failure
		    
{

  unsigned int mem_size = (*nsize)*sizeof(float);


  //allocate memory on GPU, return error code in case of problems
  *ierr=0;
  if(cudaMalloc( (void**) (GPU_pointer), mem_size) != 0)
    {
      printf("GPU allocation error \n");
      *ierr=1;
      return;
    }
}

extern "C" 
void gpu_deallocate__(float **GPU_pointer, // pointer indicating the GPU address
		      int *ierr) // error code, 1 if failure
{
  //deallocate memory on GPU, return error code in case of problems
  *ierr=0;
  if(cudaFree(*GPU_pointer) != 0)
    {
      printf("GPU deallocation error \n");
      *ierr=1;
      return;
    }
}


//Temporary send-receive operations, displacements to be added (other routines?)


extern "C"
void gpu_send__(int *nsize,
		float *CPU_pointer, 
		float **GPU_pointer,
		int *ierr)
{

  unsigned int mem_size = (*nsize)*sizeof(float);

  //copy V to GPU
  *ierr=0;
  if(cudaMemcpy(*GPU_pointer, CPU_pointer, mem_size, cudaMemcpyHostToDevice)  != 0)
    {
      printf("HostToDevice Memcpy error \n");
      *ierr=1;
      return;
    }

}

extern "C" 
void gpu_receive__(int *nsize,
		float *CPU_pointer, 
		float **GPU_pointer,
		int *ierr)
{

  unsigned int mem_size = (*nsize)*sizeof(float);

  //copy V to GPU
  *ierr=0;
  if(cudaMemcpy(CPU_pointer,*GPU_pointer, mem_size, cudaMemcpyDeviceToHost)  != 0)
    {
      printf("DeviceToHost Memcpy error \n");
      printf(" %i \n",mem_size);
      *ierr=1;
      return;
    }

}
