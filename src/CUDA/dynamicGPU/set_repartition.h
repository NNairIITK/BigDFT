#ifndef __setrepartitionh__
#define __setrepartitionh__

#include "exceptions.h"
#include "manage_global_var.h"
class set_repartition
{

public:
  set_repartition(int numActors,int numGPU,int iproc_,global_gpu_attach* gga_)
    :NUM_ACTORS(numActors),NUM_GPU(numGPU),iproc(iproc_),gga(gga_){}

  virtual void do_repartition(int currNum) const = 0;

protected:

  const int NUM_ACTORS;
  const int NUM_GPU;
   const int iproc;
  global_gpu_attach* gga;
};


class set_repartition_shared : public set_repartition
{

public:
  set_repartition_shared(int numActors,int numGPU,int iproc,global_gpu_attach* gga_)
    :set_repartition(numActors,numGPU,iproc,gga_){}

  virtual void do_repartition(int currNum) const;
};

class set_repartition_static : public set_repartition
{

public:
  set_repartition_static(int numActors,int numGPU,int iproc,global_gpu_attach* gga_)
    :set_repartition(numActors,numGPU,iproc,gga_)
  {}

  virtual void do_repartition(int currNum) const;


};
#endif
