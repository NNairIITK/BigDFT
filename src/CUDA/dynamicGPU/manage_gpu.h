/*
** manage_gpu.h
** 
** Made by Matthieu Ospici
** Login   <mo219174@badiane>
** 
** Started on  Wed Apr 16 15:28:45 2008 Matthieu Ospici
** Last update Wed Apr 16 15:28:45 2008 Matthieu Ospici
*/

#ifndef   	MANAGE_GPU_H_
# define   	MANAGE_GPU_H_

#include <vector>
#include <algorithm>

#include <pthread.h>
#include <semaphore.h>



#include "gpu_control.h"
#include "fct_call.h"
#include "class_utils.h"






class manage_gpu
{
public:
  manage_gpu(int numGPU);

  ~manage_gpu();

  int getNumGPU() const {return numGPU;}

  gpu_control& getGPUcontrol(int gpuID);
 

  //private:


  int numGPU;
  unsigned int dim1Max,dim2Max,dim3Max;
  
 
  void deblockOneThread(int tid);
  void waitForOneThread(int tid);
  void deblockAllThread();
  void waitForAllThread();




  std::vector<gpu_control*> v_gpu_control;

  std::vector<pthread_t> v_thread;


  void initGPU_control();

  // STRUCTURES (parameter...)
  struct s_p_InitThread
  {
    manage_gpu *thi;
    int tid;
  };


  class fct_call_setDevice : public fct_call
  {
  public:
    fct_call_setDevice(int gpuID_):gpuID(gpuID_){}
    virtual void operator()(int); 
    virtual ~fct_call_setDevice(){};

  private:

    int gpuID;

  };

 
  //END STRUCTURES

  fct_call **t_f_callThread; //object used in thread loop in order to call
                           //differents functions, one object per
                           //thread, more safe



  //low level code (thread, synchro...)

  void threadLoop(int tID);
  static void* threadRoutine(void*);
  //  void doMallocThread(void *); //init memory, and wait



 

  sem_t *t_sem1; //not a vector because vector are not thread-safe
   sem_t *t_sem2; //not a vector because vector are not thread-safe

   void **t_threadP; //not a vector because vector are not thread-safe
  float **t_idata;
   float **t_odata;


  //  friend class gpu_control;
  //    friend class gpu_convolution ;
};


#endif 	    /* !MANAGE_GPU_H_ */ 
