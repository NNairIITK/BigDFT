#ifndef FCT_CALL_IMPLEMENTATION_H
#define FCT_CALL_IMPLEMENTATION_H

#include "fct_call.h"
#include "init_network.h"
#include "sg_common_def.h"




class fct_call_trsf_CPU_GPU : public fct_call
{
public:
  fct_call_trsf_CPU_GPU(const void *src_, GPU_ptr dest_, size_t mem_size_)
 :src(src_), dest(dest_),mem_size(mem_size_){}
  virtual void operator()(int); 
  virtual ~fct_call_trsf_CPU_GPU(){};
  
private:
  
  const void *src;
  GPU_ptr dest;
  
  size_t mem_size;
  

  
};

class fct_call_trsf_GPU_CPU : public fct_call
{
public:
  fct_call_trsf_GPU_CPU(const GPU_ptr src_, void* dest_, size_t mem_size_)
 :src(src_), dest(dest_),mem_size(mem_size_){}
  virtual void operator()(int); 
  virtual ~fct_call_trsf_GPU_CPU(){};
  
private:
  
  const GPU_ptr src;
  void *dest;
  
  size_t mem_size;
  
  
 
};


class fct_call_memcpy : public fct_call
{
public:
  fct_call_memcpy(const void *src_, void* dest_, size_t mem_size_)
 :src(src_), dest(dest_),mem_size(mem_size_){}
  virtual void operator()(int); 
  virtual ~fct_call_memcpy(){};
  
private:
  
  const void *src;
  void *dest;  
  size_t mem_size;

};


class fct_call_calc_generic : public fct_call
{
public:
  fct_call_calc_generic(sg_callback_ptr f_call_, void *param_, size_t size_param_)
    :f_call(f_call_),size_param(size_param_)
  {malloc_and_copy(param_,size_param_);}
  virtual void operator()(int); 
  virtual ~fct_call_calc_generic();

  
private:
  void malloc_and_copy(void* src,size_t size);
  sg_callback_ptr f_call;
  void *local_param;
  size_t size_param;
};
#endif
