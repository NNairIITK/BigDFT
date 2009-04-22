#ifndef GPU_CONTROL_H
#define GPU_CONTROL_H


#include "fct_call.h"

class manage_gpu;




class gpu_control
{
public:
  gpu_control(manage_gpu& mg_,int gpuID_):gpuID(gpuID_),mg(mg_){}
  
  void changeCPU_affinity(int CPU);

  int getGPU_id() const {return gpuID;}


  void changeThreadOperation(fct_call*);
  void deblockThread();
  void waitForOneThread();

  float *getIdata();
  float *getOdata();
  
  // virtual void tmp(){};
  // virtual ~gpu_control(){};
protected:
  const int gpuID;
  manage_gpu& mg;


  class fct_call_affinity : public fct_call
  {
  public:
    fct_call_affinity(int CPUaff_):CPUaff(CPUaff_){}
    virtual void operator()(int tid) ; 
    virtual ~fct_call_affinity(){};

  private:
    int CPUaff;
 

  
  };

 
  //parameter class for dialling with thread
};

#endif
