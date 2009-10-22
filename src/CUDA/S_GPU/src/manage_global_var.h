#ifndef __manageglobalvarh__
#define __manageglobalvarh__


//store if process/thread is attached to GPU
//this version is for process
class global_gpu_attach
{
public:
  global_gpu_attach()
    :isAttached(false){}
  void setIsAttached(bool v) {isAttached = v;}
  bool getIsAttached() const {return isAttached;}

private:
  bool isAttached;
};




#endif
