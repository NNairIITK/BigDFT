#include <iostream>
#include <sched.h>
#include <pthread.h>
#include <semaphore.h>


#include "gpu_control.h"

#include "manage_gpu.h"



void gpu_control::fct_call_affinity::operator()(int)
{
 
  //set thread affinity
  /*  int CPUaff2 = 0;
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

  */
  
  //std::cout << "thid real : " << pthread_self() << std::endl;
}
void gpu_control::changeCPU_affinity(int CPU)
{
  mg.t_f_callThread[gpuID] = new fct_call_affinity(CPU);
  
  mg.deblockOneThread(gpuID);
  mg.waitForOneThread(gpuID);

}

void gpu_control::changeThreadOperation(fct_call *fctc)
{
  mg.t_f_callThread[gpuID] = fctc;
}

void gpu_control::deblockThread()
{
  mg.deblockOneThread(gpuID);
}

void gpu_control::waitForOneThread()
{
  mg.waitForOneThread(gpuID);
}

float* gpu_control::getIdata()
{
  return mg.t_idata[gpuID];
}

float* gpu_control::getOdata()
{
return mg.t_odata[gpuID];
}
