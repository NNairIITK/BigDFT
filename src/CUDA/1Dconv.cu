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

//#include "conv_shared_multi_kernel.cu"
//#include "deffct.h" //convseria

__constant__ param_t param;

__global__ void conv_shared_multi(unsigned int n1,unsigned int n2,float *t_out,float *t_in,int nf);



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
