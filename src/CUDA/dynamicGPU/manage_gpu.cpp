#include <iostream>

#include "manage_gpu.h"
#include "cudafct.h"

using std::cout;
using std::endl;

manage_gpu::manage_gpu( int currGPU)
{
  
  /* c_cudaGetDeviceCount(&numGPU);
  dim1Max = dim1Max_;
  dim2Max = dim2Max_;
  dim3Max = dim3Max_;*/

  numGPU = 1; //one thread only

  v_thread.resize(numGPU);

  t_sem1 = new sem_t[numGPU];
  t_sem2 = new sem_t[numGPU];
  t_threadP = new void*[numGPU];

  t_idata = new float*[numGPU];
  t_odata = new float*[numGPU];


  t_f_callThread = new fct_call*[numGPU];
 
  for(int i=0;i<numGPU;++i)
    {
      sem_init(&t_sem1[i], 0, 0);
      sem_init(&t_sem2[i], 0, 0);
    }

 
  //create threads
  // Initialize and set thread detached attribute
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  // pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);

  // int *tmpParam = new int[numGPU];
  
  s_p_InitThread  *tmpParam = new s_p_InitThread[numGPU];
  for(int i = 0; i < numGPU; i++)
    {
      tmpParam[i].tid = i;
      tmpParam[i].thi = this;
      pthread_create(&v_thread.at(i), &attr,threadRoutine, (void*)&tmpParam[i]);
     
    }


  //call thread for allocating memory
  //1. put fct_call_malloc in t_f_callThread

  for(int i=0;i<numGPU;++i)
    {
      t_f_callThread[i] = new fct_call_setDevice(currGPU);
    }

 
  //deblock all semaphores in order to execute fcts

  deblockAllThread();
  //wait for all semaphore reblock in order to be sure that all thread
  //has finish execution of doMallocThread

  waitForAllThread();


  //now all thread has finished, we can delete the two tabs
  // std::cout << "FINISH" << endl;
  //  delete[] tmpParamCudaMalloc;
  delete[] tmpParam;


  std::for_each(t_f_callThread,t_f_callThread+(numGPU-1),deleter());
  initGPU_control();
}

manage_gpu::~manage_gpu()
{

  std::for_each(v_gpu_control.begin(),v_gpu_control.end(),deleter());
  delete[] t_sem1;
  delete[] t_sem2;
  delete[] t_threadP;
  delete[] t_idata;
  delete[] t_odata;
  delete[] t_f_callThread;
}


void manage_gpu::initGPU_control()
{
 
  for(int i=0;i<numGPU;++i)
    {
        v_gpu_control.push_back(new gpu_control(*this,i));
    }
}


gpu_control& manage_gpu::getGPUcontrol(int gpuID)
{
  return *v_gpu_control.at(gpuID);
}


void manage_gpu::deblockOneThread(int tid)
{
  sem_post(&t_sem1[tid]);
}
void manage_gpu::waitForOneThread(int tid)
{
  sem_wait(&t_sem2[tid]);
}
void manage_gpu::deblockAllThread()
{
  for(int i=0;i<numGPU;++i)
    {
      sem_post(&t_sem1[i]);
    }
}
void manage_gpu::waitForAllThread()
{
  for( int i=0;i<numGPU;++i)
    {
      sem_wait(&t_sem2[i]);
    }
}

//**** BEGIN THREAD FUNCTIONs *****//


void manage_gpu::threadLoop(int tID)
{
  while(true)
    {
      sem_wait(&t_sem1[tID]);
      
      //    if(t_threadP[tID] != NULL)
      //p_fct(t_threadP[tID]); //call the function
      (*t_f_callThread[tID])(tID); //call the function, one objet per thread
      sem_post(&t_sem2[tID]);
    }
}



void* manage_gpu::threadRoutine(void* param)
{

  //reinterpret_cast<manage_gpu*>(param)->threadLoop();
  s_p_InitThread  *p = (reinterpret_cast<s_p_InitThread*>(param));
  

  //  std::cout << "thid real : " << pthread_self() << " logic : " << p->tid << std::endl;
  
  p->thi->threadLoop(p->tid); //the only way to call a non-static fct...
  
  
  return NULL;
}

//**** END THREAD FUNCTIONs *****//



//*** PARAM CLASS DEFS *****
void manage_gpu::fct_call_setDevice::operator()(int tID)
{
 
  /*  int CPUaff = tID;
  cpu_set_t cpus;
  CPU_ZERO(&cpus);
  
  CPU_SET(CPUaff, &cpus);
  
  std::cout << "on set proc - " << CPUaff << std::endl;

  if(sched_setaffinity(0, sizeof(cpu_set_t),&cpus) < 0)
    {
      printf("erreur setaffinity\n");
    }
  
  if(!CPU_ISSET(CPUaff, &cpus))
    {
      printf("*** CPU - %i - pas setté ***\n",CPUaff);
    }
  */
  
  // if(c_cudaSetDevice(1) != 0)
    if(c_cudaSetDevice(gpuID) != 0)
    cout << "Set device ERROR !" << std::endl;

  cout << pthread_self()<< "Set device : tid : " << gpuID << std::endl;

  //set proc
  int CPUaff2;

  if(gpuID != 0)
    CPUaff2 = 7;
  else
    CPUaff2 =1;

  cpu_set_t cpus;
  CPU_ZERO(&cpus);
  
  CPU_SET(CPUaff2, &cpus);
  
  std::cout << "on set proc - " << CPUaff2 << std::endl;

   if(sched_setaffinity(0, sizeof(cpu_set_t),&cpus) < 0)
    {
      printf("erreur setaffinity\n");
    }
  
  if(!CPU_ISSET(CPUaff2, &cpus))
    {
      printf("*** CPU - %i - pas setté ***\n",CPUaff2);
      }
  

  // cout << "cuda malloc " <<   memSize << "device " << gpuID<< endl;


  /*  if(c_cudaMalloc( (void**) &(base->t_idata[tID]), memSize) != 0)
    {
      printf("erreur malloc GPU_odata\n");
      return;
    }
 
  if( c_cudaMalloc( (void**) &(base->t_odata[tID]),  memSize) != 0)
    {
      printf("erreur malloc GPU_idata\n");        
      return;
      } */
}
