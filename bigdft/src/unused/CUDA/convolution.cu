/****u* CUDA/convolution.cu
**
** 
** AUTHOR
**  Matthieu Ospici
**
** CHANGELOG
** Started on  Tue Mar 11 14:54:55 2008 Matthieu Ospici
**
** Last update Fri Jun 27 17:12:28 2008 Luigi Genovese
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

#include "conv_shared_multi_kernel.cu"
//#include "deffct.h" //convseria
#include "convSeria.h"

pthread_mutex_t mutexTime = PTHREAD_MUTEX_INITIALIZER; //for timers

//pthread_mutex_t mutexCond;
//pthread_cond_t pcond;

/*static sem_t sem;
static sem_t sem2;
static sem_t sem3;

static pthread_t *threads; 
//beark beark !!!
*/


// Functions in libcutil.a for timing
void cutCreateTimerM(unsigned int *timer)
{
  pthread_mutex_lock (&mutexTime);

  //cutCreateTimer(timer);
  
  pthread_mutex_unlock (&mutexTime);

}

void cutStopTimerM(unsigned int timer)
{
  pthread_mutex_lock (&mutexTime);

  //cutStopTimer(timer);

  pthread_mutex_unlock (&mutexTime);
}

void cutStartTimerM(unsigned int timer)
{
  pthread_mutex_lock (&mutexTime);

  //cutStartTimer(timer);

  pthread_mutex_unlock (&mutexTime);
}

void cutResetTimerM(unsigned int timer)
{
  pthread_mutex_lock (&mutexTime);

  //cutResetTimer(timer);

  pthread_mutex_unlock (&mutexTime);
}

float cutGetTimerValueM(unsigned int timer)
{
  pthread_mutex_lock (&mutexTime);

  //float tmp =  cutGetTimerValue(timer);
  float tmp = 0.0;
 
  pthread_mutex_unlock (&mutexTime);
  return tmp;

}


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
void intertamponcgpu_(int *n1,
		      int *n2, 
		      int *n3, 
		      float *t_in, 
		      float *t_out,
		      float *f_data,
		      int *lowfil, 
		      int *lupfil)
{

  const int fsize = *lupfil - *lowfil + 1;

  multiTab_t m_dataIn;
  multiTab_t m_dataOut;
  


  m_dataIn.n1 = *n1 + 1;
  m_dataIn.n2 = *n2 + 1;
  m_dataIn.n3 = *n3 + 1;

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
  



  conv3dGPU(&m_dataIn,
	    &m_dataOut,
	    f_data,
	    fsize,
	    *lowfil,
	    *lupfil,
	    256,
	    &evPerf);
  
}
//********** END INTERFACE BETWEEN FORTRAN & C ************





inline void printErrorPlus()
{
  cudaError_t err = cudaGetLastError();   
  if( cudaSuccess != err)
    {                                                
      printf("Cuda error:  in file '%s' in line %i : %s.\n",
	      __FILE__, 
	      __LINE__, 
	      cudaGetErrorString( err) );
           
    }
}


int conv3dGPU(multiTab_t* m_dataIn,
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

 unsigned int timer;
  cutCreateTimerM(&timer);
  cutResetTimerM(timer);

  cutStartTimerM(timer);
  float* res_data = m_dataOut->tab;



  cutStopTimerM(timer);


  evPerf->host_malloc = cutGetTimerValueM(timer);

  cutResetTimerM(timer);
 
  //  cutStartTimerM(timer);


  //allocation
  float* GPU_idata = m_dataIn->GPU_data;  
  // printf("GPU_idata %p\n",GPU_idata);
  float* GPU_odata = m_dataOut->GPU_data;  

 
  evPerf->cuda_malloc =0;

  cutResetTimerM(timer);
  //copy
  cutStartTimerM(timer);

  if(cudaMemcpy( GPU_idata, m_dataIn->tab, mem_size, cudaMemcpyHostToDevice)  != 0)
    {
      printf("erreur memcpy GPU_idata\n");
      cudaFree(GPU_odata);
      cudaFree(GPU_idata);
      return 1;
    }

  cutStopTimerM(timer);


  evPerf->GPU_trsf_CPUGPU = cutGetTimerValueM(timer);

  cutResetTimerM(timer);
  

  param_t paramToSend;
  unsigned int numBlockDim1,numBlockDim2;
  
  //************ Code execution

    cutStartTimerM(timer);

 // ----------- STEP 1 --------



    unsigned int n1 = m_dataIn->n1*m_dataIn->n2;
    unsigned int n2 = m_dataIn->n3; 


  createParam(&paramToSend,num_threads,n1,n2,&numBlockDim1,&numBlockDim2,f_data,fsize,lowfil,lupfil);

 
  if( cudaMemcpyToSymbol(param,&paramToSend, sizeof(param_t)) != 0)
    {

      printf("Erreur copy to symbol param\n");

      printErrorPlus();

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

  // ----------- STEP 2 --------

    n1 = m_dataIn->n3 * m_dataIn->n1;
    n2 = m_dataIn->n2;
 

   createParam(&paramToSend,num_threads,n1,n2,&numBlockDim1,&numBlockDim2,f_data,fsize,lowfil,lupfil);
 
 
  if( cudaMemcpyToSymbol(param,&paramToSend, sizeof(param_t)) != 0)
    {
      printf("Erreur copy to symbol param\n");
      printErrorPlus();
      cudaFree(GPU_odata);
      cudaFree(GPU_idata);
      return 1;
      }
  
  dim3  grid2(numBlockDim1, numBlockDim2, 1);  
  dim3  threads2(num_threads, 1, 1);
  


  conv_shared_multi <<< grid2, threads2 >>>(n1,n2,GPU_idata, GPU_odata,fsize);
  cudaThreadSynchronize();
  

  // ----------- END STEP 2 --------

  // ----------- STEP 3 --------

  n1 = m_dataIn->n3 * m_dataIn->n2;
  n2 = m_dataIn->n1 ;

  createParam(&paramToSend,num_threads,n1,n2,&numBlockDim1,&numBlockDim2,f_data,fsize,lowfil,lupfil);
 
 
  if( cudaMemcpyToSymbol(param,&paramToSend, sizeof(param_t)) != 0)
    {
      printf("Erreur copy to symbol param\n");
      printErrorPlus();
      cudaFree(GPU_odata);
      cudaFree(GPU_idata);
      return 1;
      }
  
  dim3  grid3(numBlockDim1, numBlockDim2, 1);  
  dim3  threads3(num_threads, 1, 1);

 
    conv_shared_multi <<< grid3, threads3 >>>(n1,n2,GPU_odata, GPU_idata,fsize);
   cudaThreadSynchronize();

  // ----------- END STEP 3 --------



    cutStopTimerM(timer);


  evPerf->GPU_calc = cutGetTimerValueM(timer);

  cutResetTimerM(timer);

  //************* COPY back to CPU
 

  cutStartTimerM(timer);
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
  cutStopTimerM(timer);


  evPerf->GPU_trsf_GPUCPU = cutGetTimerValueM(timer);
  
  cutResetTimerM(timer);


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
/****/
