#include <iostream>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include <sstream>

#include "exceptions.h"

#include "init_network.h"
#include "message.h"
#include "manage_gpu.h"

#include "class_utils.h" //for deleter

#include "read_conf_file.h" 

#include "localqueu.h"
#include "convolution_fct_call.h"
#include "manage_cpu_affinity.h"

extern short gpu_precision;
unsigned int getPrecisionSize();


local_network *l = NULL;
localqueu *locq = NULL;

sem_unix *sem_gpu_CALC;
sem_unix *sem_gpu_TRSF;




void init_gpu_sharing(readConfFile &read_conf, int iproc, int *error)
{
  *error = 0;
  try
    {

      int mpi_tasks_per_node,num_GPU;
     
      //read file
      // readConfFile read_conf(NAME_FILE);
      read_conf.get("MPI_TASKS_PER_NODE",&mpi_tasks_per_node);
      read_conf.get("NUM_GPU",&num_GPU);
   
      manage_cpu_affinity mca(iproc);

      for(int i=0;i<num_GPU;++i)
	{
	  std::ostringstream iss;
	  std::string aff;
	  iss << "GPU_CPUS_AFF" << "_" << i;
	  read_conf.get(iss.str(),aff);
	  
	  mca.add_connexion(gpu_cpus_connexion(aff));

	}

      //init node
      l = new local_network(mpi_tasks_per_node,num_GPU,mca,iproc);

      locq = new localqueu();
      sem_gpu_CALC = l->getSemCalc();
      sem_gpu_TRSF = l->getSemTrsf();
    }

  catch(synchronization_error se)
    {
      std::cerr << "*** ERROR(s) DETECTED AT THE INITIALIZATION OF THE INTER-NODE COMMUNICATION SYSTEM ***" << std::endl;
      std::cerr << "ERROR MESSAGE : " << se.what() << std::endl;
      *error = 1;
    }

  catch(inter_node_communication_error ie)
    {
      std::cerr << "*** ERROR(s) DETECTED AT THE INITIALIZATION OF THE INTER-NODE COMMUNICATION SYSTEM ***" << std::endl;
      std::cerr << "ERROR MESSAGE : " << ie.what() << std::endl;
      *error = 1;
    }


  catch(read_not_found re)
    {
      std::cerr << "*** ERROR : INVALID CONFIG FILE. You have to set the number of mpi tasks per node and the number of GPU to use per node ***" << std::endl;
      std::cerr << "Missing information : " << re.what() << std::endl;
      *error = 1;
    }

  catch(file_not_found fe)
    {
      std::cerr << "*** ERROR : CONFIG FILE NOT FOUND" << std::endl;
      std::cerr << "File not found : " << fe.what() << std::endl;
      *error = 1;
    }



  catch(check_calc_error cce)
    {
      std::cerr << "*** ERROR : HARDWARE PROBLEME ON A CARD" << std::endl;
      std::cerr << "We have send calculations to a card and the result was bad. *** Hostname " << cce.what() << "***" << std::endl;
      *error = 1;
    }




  catch(std::exception e)
    {
      std::cerr << "*** ERROR(s) DETECTED AT THE INITIALIZATION OF THE INTER-NODE COMMUNICATION SYSTEM ***" << std::endl;
      std::cerr << "ERROR MESSAGE : " << e.what() << std::endl;
      *error = 1;
    }
 
}


extern "C" 
void stop_gpu_sharing__()
{
  //  std::for_each(locq.begin(),locq.end(),deleter());
  delete locq;
  delete l;
}


extern "C"
void create_stream__(gpu_stream **new_stream)
{
  *new_stream = new gpu_stream(l->getCurrGPU());
  locq->addStream(*new_stream);
}


extern "C"
void launch_all_streams__()
{


  
  l->messageLoopNetwork(*locq);


  //now we can empty stream
  locq->removeStreams();
}


extern "C"
void gpu_pack_unpack_stream__( int *ierr,
			       gpu_stream **stream)
{
  *ierr=0;

  fct_call_pack_unpack *pack_unpack= 
    new fct_call_pack_unpack();


  (*stream)->addOp(pack_unpack,PACK);
}


		 
extern "C"
void gpu_send_pi_stream__(int *nsize,
			  void **CPU_pointer, 
			  void **GPU_pointer,
			  int *ierr,
			  gpu_stream **stream)
{
  //  std::cout << "send_pi debut" << std::endl
  //	    << "dst : " << *GPU_pointer <<",src " << *CPU_pointer << std::endl;
  unsigned int mem_size = (*nsize)*sizeof(double);
  

  *ierr=0;

  fct_call_trsf_CPU_GPU *trsfCPU_GPU = 
    new fct_call_trsf_CPU_GPU(*CPU_pointer,*GPU_pointer,mem_size,*l);


  (*stream)->addOp(trsfCPU_GPU,TRANSF);
  //  std::cout << "send_pi END " << std::endl;
}

extern "C"
void gpu_receive_pi_stream__(int *nsize,
			     void **CPU_pointer, 
			     void **GPU_pointer,
			     int *ierr,
			     gpu_stream **stream)
{
  //  std::cout << "recv pi START " << std::endl;

  unsigned int mem_size = (*nsize)*sizeof(double);
 
  *ierr=0;

  fct_call_trsf_GPU_CPU *trsfGPU_CPU = 
    new fct_call_trsf_GPU_CPU(*GPU_pointer,*CPU_pointer,mem_size,*l);
 

  (*stream)->addOp(trsfGPU_CPU,TRANSF);
  // std::cout << "recv pi END " << std::endl;
}


extern "C"
void mem_copy_f_to_c_stream__(int *nsize, //memory size
			      void **dest,
			      void *srcFortran,
			      int *ierr,
			      gpu_stream **stream)
{
  *ierr=0;

  unsigned int mem_size = (*nsize)*getPrecisionSize();

  fct_call_trsf_memcpy_f_to_c *trsf_f_c =
    new fct_call_trsf_memcpy_f_to_c(mem_size,dest,srcFortran);


 (*stream)->addOp(trsf_f_c,TRANSF);


}

extern "C"
void mem_copy_c_to_f_stream__(int *nsize, //memory size
			      void *destFortran,
			      void **src,
			      int *ierr,
			      gpu_stream **stream)
{

  *ierr=0;

  unsigned int mem_size = (*nsize)*getPrecisionSize();
  fct_call_trsf_memcpy_c_to_f *trsf_c_f =
    new fct_call_trsf_memcpy_c_to_f(mem_size,destFortran,src);
 
  (*stream)->addOp(trsf_c_f,TRANSF);
}

//-------stream version of calculation-----

extern "C"
void gpulocden_stream__(int *n1,int *n2, int *n3,int *norbp,int *nspin,
			double *h1,double *h2,double *h3,
			double *occup,double *spinsgn,
			double **psi,int **keys, 
			double **work1,double **work2,
			double **rho,
			gpu_stream **stream)
{

  fct_call_calc_locden<double> *calclocden =
    new fct_call_calc_locden<double>(n1,n2,n3,norbp,nspin,
				     h1,h2,h3,
				     occup,spinsgn,
				     psi,keys, 
				     work1,work2,
				     rho);

 (*stream)->addOp(calclocden,CALC);

}

extern "C" 
void gpulocham_stream__(int *n1,int *n2, int *n3,
			double *h1,double *h2,double *h3,
			double **psi,double **pot,int **keys, 
			double **work1,double **work2,double **work3,
			double *epot_sum,double *ekin_sum,
			double *ocupGPU,
			gpu_stream **stream)
{
  
  fct_call_calc_hamiltonian<double> *calcham =
    new fct_call_calc_hamiltonian<double>(n1,n2,n3,
					  h1,h2,h3,
					  psi,pot,keys,
					  work1, work2, work3,
					  epot_sum,ekin_sum,ocupGPU);

 (*stream)->addOp(calcham,CALC);
}


extern "C" 
void gpuprecond_stream__(int *n1,int *n2, int *n3,int *npsi,
			 double *h1,double *h2,double *h3,
			 double **x,int **keys, 
			 double **r,double **b,double **d,
			 double **work1,double **work2,double **work3,
			 double *c,int *ncong, double *gnrm,
			 gpu_stream **stream)
{

  // std::cout  << " gpuprec_stream : h2 : " << *h2 << ", npsi " << *npsi << std::endl;
  fct_call_calc_precond<double> *calcprecond =
    new fct_call_calc_precond<double>(n1,n2, n3,npsi,
				      h1,h2,h3,
				      x,keys, 
				      r,b,d,
				      work1,work2,work3,
				      c,ncong, gnrm);

  (*stream)->addOp(calcprecond,CALC);
}




// =======================================================
extern "C" 
void mem_copy_f_to_c__(int *nsize, //memory size
		       void **dest,
		       void *srcFortran,
		       int *ierr) // error code, 1 if failure

		    
{

  unsigned int mem_size = (*nsize)*getPrecisionSize();
  *ierr=0;
  
  memcpy(*dest,srcFortran, mem_size);

  
}

extern "C" 
void mem_copy_c_to_f__(int *nsize, //memory size
		       void *destFortran,
		       void **src,
		       int *ierr) // error code, 1 if failure

		    
{

unsigned int mem_size = (*nsize)*getPrecisionSize();
  *ierr=0;

  memcpy(destFortran,*src, mem_size);

  
  }



/*extern "C" 
void gpu_allocate_new__(int *nsize, //memory size
			void **GPU_pointer, // pointer indicating the GPU address
			int *ierr) // error code, 1 if failure		    
{
  //  if(l->getCurr() == 0)
    {
      //      sem_gpu->P();
  unsigned int mem_size = (*nsize)*sizeof(double);

  std::cout <<  " Debut  gpu_allocate_new__" << std::endl;

  //allocate memory on GPU, return error code in case of problems
  *ierr=0;


  //  if(l->getCurr() == 0)
  // std::cout << "entré :" ;
  //  int age;
  // std::cin >> age;


      fct_call_alloc *calal =
    new fct_call_alloc(GPU_pointer,mem_size);
  l->man_gpu->getGPUcontrol(0).changeThreadOperation(calal); 
	      
	      
	      
  l->man_gpu->getGPUcontrol(0).deblockThread();
	      
	      
	      
  l->man_gpu->getGPUcontrol(0).waitForOneThread();
    
  std::cout <<  " FIN  gpu_allocate_new__" << std::endl;
  // sem_gpu->V();
    }
}


extern "C" 
void cpu_pinned_allocate_new__(int *nsize, //memory size
			       void **CPU_pointer,
			       int *ierr) // error code, 1 if failure

		    
{
  // if(l->getCurr() == 100)
    {
      //     sem_gpu->P();

  std::cout <<  " DEBUT  gpu_pinned_allocate_new__" << std::endl;

  unsigned int mem_size = (*nsize)*sizeof(double);
  *ierr=0;


  
  fct_call_alloc_pi *calalpi =
    new fct_call_alloc_pi(CPU_pointer,mem_size);
  l->man_gpu->getGPUcontrol(0).changeThreadOperation(calalpi); 
  
  
  
  l->man_gpu->getGPUcontrol(0).deblockThread();
  
	      
  
  l->man_gpu->getGPUcontrol(0).waitForOneThread();
 std::cout <<  " FIN  gpu_pinned_allocate_new__" << std::endl;
 // sem_gpu->V();
    }
    }*/
