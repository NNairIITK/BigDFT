#include <iostream>
#include "convolution_fct_call.h"
 #include <unistd.h>
#include <pthread.h>

#include "trace_exec.h"





#include "cudafct.h"
#include "message.h"

extern trace_exec *tracer;


extern "C" 
void gpulocham_(int *n1,int *n2, int *n3,
		double *h1,double *h2,double *h3,
		double **psi,double **pot,int **keys, 
		double **work1,double **work2,double **work3,
		double *epot,double *ekin);

extern "C" 
void gpuprecond_(int *n1,int *n2, int *n3,int *npsi,
		 double *h1,double *h2,double *h3,
		 double **x,int **keys, 
		 double **r,double **b,double **d,
		 double **work1,double **work2,double **work3,
		 double *c,int *ncong, double *gnrm);


extern "C" 
void gpulocden_(int *n1,int *n2, int *n3,int *norbp,int *nspin,
		double *h1,double *h2,double *h3,
		double *occup,double *spinsgn,
		double **psi,int **keys, 
		double **work1,double **work2,
		double **rho);



void fct_call_trsf_CPU_GPU::operator()(int tID)
{
  // gettime();
  //  std::cout << pthread_self() << " BEG trsfING, tid : CPU_GPU" << tID << std::endl;


 
  // std::cout << "dest : " << dest <<",srdc : " << src << std::endl;
  write_trace wt = tracer->getNewWriteTrace(GPU_TRSF,"Trsf CPU->GPU");

  wt.begTimeCount();

  c_cuda_gpu_send_pi(dest,src,memSize);
  //gettime();
  //  std::cout << pthread_self() << " END trsfING, tid : CPU_GPU" << tID << std::endl;

  wt.endTimeCount();


  //send messages and debloc thread 

  /* message msg;
  msg.msgType = mTRANSF;
  msg.gpu_avail = loc_net.getCurrGPU();
  */

  //  loc_net.send_next(&msg);
  // loc_net.man_gpu->getGPUcontrol(0).waitForOneThread();
  
}

void fct_call_trsf_GPU_CPU::operator()(int tID)
{
  //gettime();
  //  std::cout << pthread_self() << " BEG trsfING, tid : GPU_CPU" << tID << std::endl;
  //  std::cout << "dest : " << dest <<",srdc : " << src << std::endl;
  // sleep(10);
  write_trace wt = tracer->getNewWriteTrace(GPU_TRSF,"Trsf GPU->CPU");

  wt.begTimeCount();


  c_cuda_gpu_recv_pi(dest,src,memSize);
  //gettime();
  //  std::cout << pthread_self() << "END  trsfING, tid : GPU_CPU" << tID << std::endl;

  wt.endTimeCount();


 message msg;
  msg.msgType = mTRANSF;
  msg.gpu_avail = loc_net.getCurrGPU();

  //   loc_net.send_next(&msg);
 // loc_net.man_gpu->getGPUcontrol(0).waitForOneThread();
}

void fct_call_trsf_memcpy_c_to_f::operator()(int tID)
{
  
 write_trace wt = tracer->getNewWriteTrace(CPU_MEMCPY,"Trsf c to f");

  wt.begTimeCount();
 
  // std::cout << "Begin mem C to F  cpy dest : " << destFortran << ", src : " << *src << " memsize " << mem_size << std::endl;

  memcpy(destFortran,*src, mem_size);

wt.endTimeCount();

  //  std::cout << "end memcpy C to F" << std::endl;
}

void fct_call_trsf_memcpy_f_to_c::operator()(int tID)
{
 
  write_trace wt = tracer->getNewWriteTrace(CPU_MEMCPY,"Trsf f to c");

  wt.begTimeCount();
  //  std::cout << "Begin mem F to C, cpy dest : " << *dest << ", src : " << srcFortran << " memsize " << mem_size << std::endl;
  
  memcpy(*dest,srcFortran, mem_size);
wt.endTimeCount();

  //std::cout << "end memcpy F to C" << std::endl;
}

template<typename T>
void fct_call_calc_hamiltonian<T>::operator()(int tID)
{
  std::cout << pthread_self() << " calcING, tid : " << tID << std::endl;

}

template<> 
void fct_call_calc_hamiltonian<double>::operator()(int tID)
{
 write_trace wt = tracer->getNewWriteTrace(GPU_CALC,"locham");

  wt.begTimeCount();

  //gettime();
  //  std::cout << pthread_self() << " BEG calcING, tid : " << tID << std::endl;
  double epotToAdd,ekinToAdd;
  gpulocham_(&n1,&n2,&n3,
	     &h1,&h2,&h3,
	     psi,pot,keys,			     
	     work1, work2, work3,
	     &epotToAdd,&ekinToAdd);

  // ekin_sum=ekin_sum+orbs%occup(iorb+orbs%isorb)*ekin
  //    epot_sum=epot_sum+orbs%occup(iorb+orbs%isorb)*epot

 
  *ekin_sum = *ekin_sum + ocupGPU*ekinToAdd;
 *epot_sum = *epot_sum + ocupGPU*epotToAdd;

wt.endTimeCount();

  //gettime();
  //  std::cout << pthread_self() << " END calcING, tid : " << tID << std::endl;


  //  message msg;
  // msg.msgType = mCALC;
  // msg.gpu_avail = loc_net.getCurrGPU();

  //  loc_net.send_next(&msg);
 //  loc_net.man_gpu->getGPUcontrol(0).waitForOneThread();

}


template<typename T>
void fct_call_calc_precond<T>::operator()(int tID)
{
  std::cout << pthread_self() << " YO YO YO ERROR : " << tID << std::endl;

}

template<> 
void fct_call_calc_precond<double>::operator()(int tID)
{
  //gettime();
// std::cout << pthread_self() << " BEG calcING PRECOND, tid : " << tID << std::endl;

// std::cout  << " in calc : h2 : " << h2 << ", npsi " << npsi << std::endl;
 write_trace wt = tracer->getNewWriteTrace(GPU_CALC,"precond");

  wt.begTimeCount();



  double gnrmToAdd;
  gpuprecond_(&n1,&n2, &n3,&npsi,
	      &h1,&h2,&h3,
	      x,keys, 
	      r,b,d,
	      work1,work2,work3,
	      &c,&ncong, &gnrmToAdd);

  *gnrm += gnrmToAdd;

wt.endTimeCount();

  //gettime();
  //  std::cout << pthread_self() << " END calcING PRECOND, tid : " << tID << std::endl;

}

template<typename T>
void fct_call_calc_locden<T>::operator()(int tID)
{
  std::cout << pthread_self() << " YO YO YO ERROR : " << tID << std::endl;

}

template<>
void fct_call_calc_locden<double>::operator()(int tID)
{

 write_trace wt = tracer->getNewWriteTrace(GPU_CALC,"locden");

  wt.begTimeCount();

  gpulocden_(&n1,&n2,&n3,&norbp,&nspin,
	     &h1,&h2,&h3,
	     occup,spinsgn,
	     psi,keys, 
	     work1,work2,
	     rho);
wt.endTimeCount();

}

/*void fct_call_alloc_pi::operator()(int tID)
{
  //std::cout << pthread_self() <<" PIptr : " <<  GPU_pointer << " mem si " << memSize << std::endl;
  c_cudaMallocHost(GPU_pointer,memSize);
  }

void fct_call_alloc::operator()(int tID)
{
  int ii;
  c_cuda_get_device(&ii);

  //  std::cout << pthread_self() <<" BEG  alloc, device ATTA : " <<ii << std::endl;
  // std::cout << pthread_self() <<" ptr : " <<  GPU_pointer << " mem si " << memSize << std::endl;
   c_cudaMalloc(GPU_pointer,memSize);

  
   //  std::cout << pthread_self() <<" END  alloc " << std::endl;
   }*/
