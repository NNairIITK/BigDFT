/*
** convolution.cu
** 
** Made by (Matthieu Ospici)
** Login   <ospici@badiane>
** 
** Started on  Tue Mar 11 14:54:55 2008 Matthieu Ospici
** Last update Sun May 12 01:17:25 2002 Speed Blue
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

#include "conv_shared_multi_kernel.cu"
//#include "deffct.h" //convseria
#include "convSeria.h"

pthread_mutex_t mutexTime = PTHREAD_MUTEX_INITIALIZER; //for timers

//pthread_mutex_t mutexCond;
//pthread_cond_t pcond;

/*static sem_t sem;
static sem_t sem2;
static sem_t sem3;


static pthread_t *threads; //beark beark !!!*/


void cutCreateTimerM(unsigned int *timer)
{
   pthread_mutex_lock (&mutexTime);

  cutCreateTimer(timer);
  
  pthread_mutex_unlock (&mutexTime);

}

void cutStopTimerM(unsigned int timer)
{
   pthread_mutex_lock (&mutexTime);


  cutStopTimer(timer);
  pthread_mutex_unlock (&mutexTime);
}
void cutStartTimerM(unsigned int timer)
{
   pthread_mutex_lock (&mutexTime);

  cutStartTimer(timer);

  pthread_mutex_unlock (&mutexTime);
}
void cutResetTimerM(unsigned int timer)
{
   pthread_mutex_lock (&mutexTime);

  cutResetTimer(timer);

  pthread_mutex_unlock (&mutexTime);
}

float cutGetTimerValueM(unsigned int timer)
{
   pthread_mutex_lock (&mutexTime);

  float tmp =  cutGetTimerValue(timer);
 
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


  //alocation
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


void freeMultiTab(multiTab_t* m)
{
  free(m->tab);
}

//***************** PRIVATE MEMBER *****************
void uniformiseTab( int* tab,unsigned int *offsetTab,int sizeTab,int lineSize)
{

  unsigned int tmp = 0;


  if(tab[sizeTab - 1] < tab[0])
    {
      offsetTab[0] = 0;
      for(int i=1;i<sizeTab+1;++i)
	{
	  tmp += tab[i-1];
	  offsetTab[i] = tmp;
	}
      return; //the tab is OK
    }
  int i =0;
  while(tab[sizeTab - 1] > tab[0])
    {
      --tab[sizeTab - 1];
      ++tab[i];
      
      ++i;
      if(i == sizeTab - 1)
	i=0;
    }

  tmp = 0;
  offsetTab[0] = tmp;
  for(int i=1;i<sizeTab+1;++i)
    {
      tmp += tab[i-1];
      offsetTab[i] = tmp;
    }
}

void createTabs(int tabSeq[2][40],
		unsigned int currPosTab[2][40],
		unsigned int numThread,
		unsigned int sizeLine,
		short step) //0 or 1
{

  int seqNum = (2*numThread)/32;
  int toCopy = ceilf((float)sizeLine/seqNum);


  if((int)(sizeLine - toCopy*(numThread - 16)/16) < 0)
    --toCopy;


  for(int i=0;i< seqNum - 1;++i)
    {
      tabSeq[step][ i] = toCopy;
    }
  tabSeq[step][ seqNum - 1] = sizeLine - toCopy*(numThread - 16)/16;
 
 
  uniformiseTab(tabSeq[step],currPosTab[step],seqNum,sizeLine);
}


void createParam(param_t* param,
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
















