#ifndef CONVOLUTION_FCT_CALL_H
#define CONVOLUTION_FCT_CALL_H

#include "fct_call.h"
#include "init_network.h"
#include "manage_gpu.h"




class fct_call_trsf_CPU_GPU : public fct_call
{
public:
  fct_call_trsf_CPU_GPU(const void *src_, void* dest_, unsigned int memSize_,local_network& loc_net_)
 :src(src_), dest(dest_),memSize(memSize_),loc_net(loc_net_){}
  virtual void operator()(int); 
  virtual ~fct_call_trsf_CPU_GPU(){};
  
private:
  
  const void *src;
  void* dest;
  
  unsigned int memSize;
  
  local_network& loc_net;
  
};

class fct_call_trsf_GPU_CPU : public fct_call
{
public:
  fct_call_trsf_GPU_CPU(const void *src_, void* dest_, unsigned int memSize_,local_network& loc_net_)
 :src(src_), dest(dest_),memSize(memSize_),loc_net(loc_net_){}
  virtual void operator()(int); 
  virtual ~fct_call_trsf_GPU_CPU(){};
  
private:
  
  const void *src;
  void* dest;
  
  unsigned int memSize;
  
    local_network& loc_net;
 
};






class fct_call_trsf_memcpy_c_to_f : public fct_call
{
public:
  fct_call_trsf_memcpy_c_to_f(unsigned int mem_size_,
			      void *destFortran_,
			      void **src_)
    :mem_size(mem_size_),destFortran(destFortran_),src(src_){}

  virtual void operator()(int); 
  virtual ~fct_call_trsf_memcpy_c_to_f(){};
  
private:
  unsigned int mem_size;
  void *destFortran;
  void **src;
};



class fct_call_trsf_memcpy_f_to_c : public fct_call
{
public:
  fct_call_trsf_memcpy_f_to_c(unsigned int mem_size_,
			      void **dest_,
			      void *srcFortran_)
    :mem_size(mem_size_),dest(dest_),srcFortran(srcFortran_){}

  virtual void operator()(int); 
  virtual ~fct_call_trsf_memcpy_f_to_c(){};

private:
  unsigned int mem_size;
  void **dest;
  void *srcFortran;
};



template<typename T>
class fct_call_calc_precond : public fct_call
{
public:
  fct_call_calc_precond(int *n1_,int *n2_, int *n3_,int *npsi_,
			T *h1_,T *h2_,T *h3_,
			T **x_,int **keys_, 
			T **r_,T **b_,T **d_,
			T **work1_,T **work2_,T **work3_,
			T *c_,int *ncong_, T *gnrm_)
    :n1(*n1_),n2(*n2_),n3(*n3_),npsi(*npsi_),
     h1(*h1_),h2(*h2_),h3(*h3_),
     x(x_),keys(keys_),
     r(r_),b(b_),d(d_),
     work1(work1_),work2(work2_),work3(work3_),
     c(*c_),ncong(*ncong_),gnrm(gnrm_){}

  virtual void operator()(int); 
  virtual ~fct_call_calc_precond(){};

private:

  int n1,  n2, n3, npsi;
  T h1,h2,h3;
  T **x;
  int **keys;
  T **r, **b, **d;
  T **work1, **work2, **work3;
  T c;
  int ncong;
  T *gnrm;

};

template<typename T>
class fct_call_calc_locden : public fct_call
{
public:
  fct_call_calc_locden(int *n1_,int *n2_, int *n3_,int *norbp_,int *nspin_,
		       T *h1_,T *h2_,T *h3_,
		       T *occup_,T *spinsgn_,
		       T **psi_,int **keys_, 
		       T **work1_,T **work2_,
		       T **rho_)
    :n1(*n1_),n2(*n2_),n3(*n3_),norbp(*norbp_),nspin(*nspin_),
     h1(*h1_),h2(*h2_),h3(*h3_),
     occup(occup_),spinsgn(spinsgn_),
     psi(psi_),keys(keys_),
     work1(work1_),work2(work2_),
     rho(rho_){}


  virtual void operator()(int); 
  virtual ~fct_call_calc_locden(){};
private:

  int n1,n2, n3,norbp,nspin;
  T h1,h2,h3;
  T *occup, *spinsgn; //array
  T **psi;
  int **keys;
  T **work1;
  T **work2;
  T **rho;
};




template<typename T>
class fct_call_calc_hamiltonian: public fct_call
{
public:
  fct_call_calc_hamiltonian(int *n1_,int *n2_, int *n3_,
			    T *h1_,T *h2_,T *h3_,
			    T **psi_,T **pot_,int **keys_, 
			    T **work1_,T **work2_,T **work3_,
			    T *epot_sum_,T *ekin_sum_,
			    T *ocupGPU_)
    :n1(*n1_),n2(*n2_),n3(*n3_),
     h1(*h1_),h2(*h2_),h3(*h3_),
     psi(psi_),pot(pot_),keys(keys_),
     work1(work1_), work2(work2_),work3(work3_),
     epot_sum(epot_sum_),ekin_sum(ekin_sum_),
     ocupGPU(*ocupGPU_){}

  
  virtual void operator()(int); 
  virtual ~fct_call_calc_hamiltonian(){};
  
private:
  
  int n1, n2, n3;
  T h1, h2, h3;
  T **psi,  **pot;
  int **keys;
  T **work1, **work2, **work3;
  T *epot_sum, *ekin_sum;
  T ocupGPU; //in order to perform some calculations
  
 

  
};


class fct_call_pack_unpack : public fct_call
{
  
  virtual void operator()(int) {};
  virtual ~fct_call_pack_unpack(){};
};



/*class fct_call_alloc_pi : public fct_call
{
public:
  fct_call_alloc_pi(void **GPU_pointer_, unsigned int memSize_)
 :GPU_pointer(GPU_pointer_),memSize(memSize_){}
  virtual void operator()(int); 
  virtual ~fct_call_alloc_pi(){};
  
private:
  

  void **GPU_pointer;
  int memSize;
  
  
 
  
};

class fct_call_alloc : public fct_call
{
public:
  fct_call_alloc(void **GPU_pointer_, unsigned int memSize_)
 :GPU_pointer(GPU_pointer_),memSize(memSize_){}
  virtual void operator()(int); 
  virtual ~fct_call_alloc(){};
  
private:
  

  void **GPU_pointer;
  int memSize;
  
  
 
  
  };*/
#endif
