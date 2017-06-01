/****u* CUDA/convolutionPot.cu
**
** 
** AUTHOR
**  Matthieu Ospici
**
** CHANGELOG
** Started on  Tue Mar 11 14:54:55 2008 Matthieu Ospici
**
** Last update Fri Jun 27 17:09:09 2008 Luigi Genovese
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



#include "convolutionPot.h"

#include "conv_pot_shared_kernel.cu"
//#include  "conv_shared_multi_kernel.cu"
//#include "deffct.h" //convseria
#include "convSeria.h"

//pthread_mutex_t mutexTime = PTHREAD_MUTEX_INITIALIZER; //for timers

//pthread_mutex_t mutexCond;
//pthread_cond_t pcond;

/*static sem_t sem;
static sem_t sem2;
static sem_t sem3;


static pthread_t *threads; //beark beark !!!*/


/*void cutCreateTimerM(unsigned int *timer)
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

}*/

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
int cuda_alloc_mem__(int *nOrbp,
		    int *n1,
		    int *n2, 
		    int *n3, 
		    float *t_psiS, //size numIter * (dim1*dim2*dim3)
		    float *t_V,  //size dim1*dim2*dim3
		    
		    float **tp_psi_GPU, //pointer array : size numIter *
		    float **t_V_GPU,  //size dim1*dim2*dim3
		    float **t_Work_GPU
		   


) 
{

  int n1C = *n1 + 1;
  int n2C = *n2 + 1;
  int n3C = *n3 + 1;

  unsigned int mem_size = n1C * n2C * n3C *sizeof(float);


  //allocate all tp_psi_GPU
  for (int i=0;i < *nOrbp; ++i)
    {
     
      // float *tmpYo;
      
      if(cudaMalloc( (void**) &(tp_psi_GPU[i]), mem_size) != 0) //only the data
	{
	  printf("erreur malloc tp_psi_GPU :  %i\n",i);        
	  return 1;
	} 

      //  tp_psi_GPU[i] = tmpYo;
      
    }


  //alocate V
    if(cudaMalloc( (void**) (t_V_GPU), mem_size) != 0)
    {
      printf("erreur malloc p_V_GPU\n");
      return 1;
      }



 

//alocate Work

  if(cudaMalloc( (void**) (t_Work_GPU), mem_size) != 0)
    {
      printf("erreur malloc t_Work_GPU \n");
      return 1;
    }


  /// copy all psi to GPU memory
  for (int i=0;i < *nOrbp; ++i)
    {
      if(cudaMemcpy(tp_psi_GPU[i] , t_psiS + i*mem_size, mem_size, cudaMemcpyHostToDevice)  != 0)
	{
	  printf("erreur memcpy tp_psi_GPU[%i] \n",i);

	  return 1;
	}
    }
  

  //copy V to GPU

  if(cudaMemcpy(*t_V_GPU, t_V, mem_size, cudaMemcpyHostToDevice)  != 0)
    {
      printf("erreur memcpy t_V \n");

      return 1;
    }


  return 0;
}



extern "C" 
int cuda_psi_to_vpsi__(int *nOrbp,
		      int *n1,
		      int *n2, 
		      int *n3, 
		      
		      float **tp_psi_GPU, 
		      float **t_V_GPU,
		      float **t_Work_GPU,
		      
		      float *f_data1, //On the CPU memory
		      float *f_data2, //On the CPU memory
		      
		      
		      int *lowfil1,
		      int *lupfil1,
		      
		      
		      int *lowfil2,
		      int *lupfil2)		  		    
{


 


  const int n1C = *n1 + 1;
  const int n2C = *n2 + 1;
  const int n3C = *n3 + 1;
 

  const int fsize1 = *lupfil1 - *lowfil1 + 1;
  const int fsize2 = *lupfil2 - *lowfil2 + 1;



  for(int i=0 ; i<*nOrbp ; ++i)
    {
  
      if(conv3dGPU_pot(n1C,
		       n2C,
		       n3C,
		       
		       tp_psi_GPU[i],
		       *t_V_GPU,
		       *t_Work_GPU,
		       
		       
		       
		       f_data1,
		       f_data2,
		       
		       fsize1,
		       *lowfil1,
		       *lupfil1,
		       
		       fsize2,
		       *lowfil2,
		       *lupfil2,
		       
		       256,
		       NULL) != 0)
	{
	  return 1;
	}


    }

    return 0;

}


extern "C" 
int  cuda_fetch_vpsi__(int *nOrbp,
		      int *n1,
		      int *n2, 
		      int *n3, 
		      float **tp_psi_GPU,
		      float *t_psiS)
{

  const int mem_size = (*n1+1) * (*n2+1) * (*n3+1) * sizeof(float);

  for(int i=0 ; i<*nOrbp ; ++i)
    {
      if( cudaMemcpy( t_psiS + mem_size*i,
		      tp_psi_GPU[i],
		      mem_size,
		      cudaMemcpyDeviceToHost) != 0)
	{
	  
	  printf("Erreur memcpy devtohost\n");
	  return 1;
	}
   
    }

  return 0;

} 


extern "C" 
void cuda_deallocate_mem__(int *nOrbp,
			  float **tp_psi_GPU, 
			  float **t_V_GPU,
			  float **t_Work_GPU)
{
  for(int i=0 ; i<*nOrbp ; ++i)
    {
      cudaFree(tp_psi_GPU[i]);
    }

  cudaFree(*t_V_GPU);
  cudaFree(*t_Work_GPU);
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


int conv3dGPU_pot(int dn1,
		  int dn2,
		  int dn3,

		  float *t_psi_GPU,      		  
		  float *t_V_GPU,
		  float *t_Work_GPU,
		  
		  
		  const float *f_data1,
		  const float *f_data2,

		  int fsize1,
		  int lowfil1,
		  int lupfil1,

		  int fsize2,
		  int lowfil2,
		  int lupfil2,

		  unsigned int num_threads,
		  evalPerfGPU_t *evPerf)

{


  //******* copy CPU -> GPU
 
 
  //  const unsigned int mem_size = dn1 *dn2 * dn3 * sizeof(float) ;

  // float* res_data = m_dataOut->tab;

   /* unsigned int timer;
  cutCreateTimerM(&timer);
  cutResetTimerM(timer);

  cutStartTimerM(timer);




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

  cutResetTimerM(timer);*/
  

  param_t paramToSend;
  unsigned int numBlockDim1,numBlockDim2;
  
  //************ Code execution

  //  cutStartTimerM(timer);




 // ----------- STEP 1 --------
    unsigned int n1 = dn1*dn2;
    unsigned int n2 = dn3; 


  createParam(&paramToSend,num_threads,n1,n2,&numBlockDim1,&numBlockDim2,f_data1,fsize1,lowfil1,lupfil1);

  //  printf("%i,%i,%i\n",fsize1,lowfil1,lupfil1);
 
  if( cudaMemcpyToSymbol(param,&paramToSend, sizeof(param_t)) != 0)
    {

      printf("Erreur copy to symbol param\n");

      printErrorPlus();

 
      return 1;
      }
  
  dim3  grid1(numBlockDim1, numBlockDim2, 1);  
  dim3  threads1(num_threads, 1, 1);

 
 
  //  cutStartTimer(timer);
 
   conv_shared <<< grid1, threads1 >>>(n1,n2,t_Work_GPU, t_psi_GPU,fsize1);
     cudaThreadSynchronize();
     //  printf("END\n");
    //   cutStopTimer(timer);
  // ----------- END STEP 1 --------

  // ----------- STEP 2 --------

    n1 = dn3 * dn1;
    n2 = dn2;
 

   createParam(&paramToSend,num_threads,n1,n2,&numBlockDim1,&numBlockDim2,f_data1,fsize1,lowfil1,lupfil1);
 
 
  if( cudaMemcpyToSymbol(param,&paramToSend, sizeof(param_t)) != 0)
    {
      printf("Erreur copy to symbol param\n");
      printErrorPlus();

      return 1;
      }
  
  dim3  grid2(numBlockDim1, numBlockDim2, 1);  
  dim3  threads2(num_threads, 1, 1);
  


    conv_shared <<< grid2, threads2 >>>(n1,n2,t_psi_GPU, t_Work_GPU,fsize1);
    cudaThreadSynchronize();
  

  // ----------- END STEP 2 --------

  // ----------- STEP 3 CONV + POT --------

  n1 = dn3 * dn2;
  n2 = dn1 ;

  createParam(&paramToSend,num_threads,n1,n2,&numBlockDim1,&numBlockDim2,f_data1,fsize1,lowfil1,lupfil1);
 
 
  if( cudaMemcpyToSymbol(param,&paramToSend, sizeof(param_t)) != 0)
    {
      printf("Erreur copy to symbol param\n");
      printErrorPlus();

      return 1;
      }
  
  dim3  grid3(numBlockDim1, numBlockDim2, 1);  
  dim3  threads3(num_threads, 1, 1);

 
  conv_pot_shared <<< grid3, threads3 >>>(n1,n2,t_Work_GPU, t_psi_GPU, fsize1, t_V_GPU);
 cudaThreadSynchronize();

  // ----------- END STEP 3 --------




// ----------- STEP 4 --------
  n1 = dn1*dn2;
  n2 = dn3; 


  createParam(&paramToSend,num_threads,n1,n2,&numBlockDim1,&numBlockDim2,f_data2,fsize2,lowfil2,lupfil2);

 
  if( cudaMemcpyToSymbol(param,&paramToSend, sizeof(param_t)) != 0)
    {

      printf("Erreur copy to symbol param\n");

      printErrorPlus();

 
      return 1;
      }
  
  dim3  grid4(numBlockDim1, numBlockDim2, 1);  
  dim3  threads4(num_threads, 1, 1);


 
  //  cutStartTimer(timer);
 
    conv_shared <<< grid4, threads4 >>>(n1,n2,t_psi_GPU,t_Work_GPU,fsize2);
     cudaThreadSynchronize();


  // ----------- END STEP 4 --------




// ----------- STEP 5 --------


   n1 = dn3 * dn1;
    n2 = dn2;
 

   createParam(&paramToSend,num_threads,n1,n2,&numBlockDim1,&numBlockDim2,f_data2,fsize2,lowfil2,lupfil2);
 
 
  if( cudaMemcpyToSymbol(param,&paramToSend, sizeof(param_t)) != 0)
    {
      printf("Erreur copy to symbol param\n");
      printErrorPlus();

      return 1;
      }
  
  dim3  grid5(numBlockDim1, numBlockDim2, 1);  
  dim3  threads5(num_threads, 1, 1);
  


    conv_shared <<< grid5, threads5 >>>(n1,n2,t_Work_GPU,t_psi_GPU, fsize2);
    cudaThreadSynchronize();
  


  // ----------- END STEP 5 --------



// ----------- STEP 6 --------

  n1 = dn3 * dn2;
  n2 = dn1 ;

  createParam(&paramToSend,num_threads,n1,n2,&numBlockDim1,&numBlockDim2,f_data2,fsize2,lowfil2,lupfil2);
 
 
  if( cudaMemcpyToSymbol(param,&paramToSend, sizeof(param_t)) != 0)
    {
      printf("Erreur copy to symbol param\n");
      printErrorPlus();
 
      return 1;
      }
  
  dim3  grid6(numBlockDim1, numBlockDim2, 1);  
  dim3  threads6(num_threads, 1, 1);

 
    conv_shared <<< grid3, threads3 >>>(n1,n2,t_psi_GPU,t_Work_GPU, fsize2);
   cudaThreadSynchronize();


  // ----------- END STEP 6 --------



  /*   cutStopTimerM(timer);


  evPerf->GPU_calc = cutGetTimerValueM(timer);

  cutResetTimerM(timer);*/

  //************* COPY back to CPU
 

  /*  cutStartTimerM(timer);
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
  cudaFree(GPU_idata);*/
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
/****/
