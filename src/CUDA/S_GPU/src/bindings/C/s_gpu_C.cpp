#include <iostream>
#include <vector>
#include <sstream>

#include "exceptions.h"

#include "init_network.h"

#include "read_conf_file.h" 

#include "localqueu.h"
#include "manage_cpu_affinity.h"
#include "set_repartition.h"
#include "manage_global_var.h"
#include "fct_call_implementation.h"


#include <s_gpu.h>


global_gpu_attach *g_gpu_attach = NULL;

local_network *l = NULL;
localqueu *locq = NULL;

sem_unix *sem_gpu_CALC;
sem_unix *sem_gpu_TRSF;


extern "C"
int sg_init(bool *GPUshare, bool *useGPU,int iproc)
{
  *useGPU = false;

  try
    {
      g_gpu_attach = new global_gpu_attach();

      const char *NAME_FILE = "GPU.config";
      readConfFile read_conf(NAME_FILE);


      int mpi_tasks_per_node,num_GPU;
      int use_shared;
      //  int iconv_param,iblas_param;

      //read file

      read_conf.get("MPI_TASKS_PER_NODE",&mpi_tasks_per_node);
      read_conf.get("NUM_GPU",&num_GPU);
      read_conf.get("USE_SHARED",&use_shared);
      //  read_conf.get("USE_GPU_BLAS",&iblas_param);
      //  read_conf.get("USE_GPU_CONV",&iconv_param);


      manage_cpu_affinity mca(iproc);

      for(int i=0;i<num_GPU;++i)
	{
	  std::ostringstream iss;
	  std::string aff;
	  iss << "GPU_CPUS_AFF" << "_" << i;
	  read_conf.get(iss.str(),aff);
	  
	  mca.add_connexion(gpu_cpus_connexion(aff));

	}


      set_repartition *set_r;
      if(use_shared == 1)
	{
	  set_r = new set_repartition_shared(mpi_tasks_per_node,num_GPU,iproc,g_gpu_attach);
	  *GPUshare = true;
	}
      else
	{
	  set_r = new set_repartition_static(mpi_tasks_per_node,num_GPU,iproc,g_gpu_attach);
	  *GPUshare = false;
	}
      //init node
      l = new local_network(mpi_tasks_per_node,num_GPU,mca,set_r,iproc);



      if(iproc == 0)
	std::cout << "Check card on all nodes...." << std::endl;

      //disable GPU for tasks that not need it
      if(g_gpu_attach->getIsAttached())
	{
	  //check the card precision, in order to detect error
	  //call a fortran function in check_card/check_init.f90
	  // checker::runTestOne(); //check only if the card has one GPU...


	  *useGPU = true;

	}


      //print repartition affinity
      if(iproc == 0)      
	mca.print_affinity_matrix();


      delete set_r; //ugly, to change...
      locq = new localqueu();
      sem_gpu_CALC = l->getSemCalc();
      sem_gpu_TRSF = l->getSemTrsf();
    }

  catch(synchronization_error& se)
    {
      std::cerr << "*** ERROR(s) DETECTED AT THE INITIALIZATION OF THE INTER-NODE COMMUNICATION SYSTEM ***" << std::endl;
      std::cerr << "ERROR MESSAGE : " << se.what() << std::endl;
      return 1;
    }

  catch(inter_node_communication_error& ie)
    {
      std::cerr << "*** ERROR(s) DETECTED AT THE INITIALIZATION OF THE INTER-NODE COMMUNICATION SYSTEM ***" << std::endl;
      std::cerr << "ERROR MESSAGE : " << ie.what() << std::endl;
      return 1;
    }


  catch(read_not_found& re)
    {
      std::cerr << "*** ERROR : INVALID CONFIG FILE. You have to set the number of mpi tasks per node and the number of GPU to use per node ***" << std::endl;
      std::cerr << "Missing information : " << re.what() << std::endl;
      return 1;
    }

  catch(file_not_found& fe)
    {
      std::cerr << "*** ERROR : CONFIG FILE NOT FOUND" << std::endl;
      std::cerr << "File not found : " << fe.what() << std::endl;
      return 1;
    }



  catch(check_calc_error& cce)
    {
      std::cerr << "*** ERROR : HARDWARE PROBLEME ON A CARD" << std::endl;
      std::cerr << "We have send calculations to a card and the result was bad. *** Hostname " << cce.what() << "***" << std::endl;
      return 1;
    }




  catch(std::exception& e)
    {
      std::cerr << "*** ERROR(s) DETECTED AT THE INITIALIZATION OF THE INTER-NODE COMMUNICATION SYSTEM ***" << std::endl;
      std::cerr << "ERROR MESSAGE : " << e.what() << std::endl;
      return 1;
    }
 

  catch(...)
    {
      std::cerr<< "** Unexpected exception "<< std::endl;
      return 1;

    }

  return 0;

  // std::ostringstream ostr;
  // ostr << "trace_" << iproc;
  // tracer = new trace_exec(ostr.str(),false);
}



extern "C"
sg_stream_ptr create_stream()
{
  gpu_stream *new_stream = new gpu_stream(l->getCurrGPU());
  locq->addStream(new_stream);

  return (sg_stream_ptr)new_stream;
}


extern "C"
void sg_launch_all_streams()
{
  l->messageLoopNetwork(*locq);


  //now we can remove the empty streams
  locq->removeStreams();
}


extern "C"
int sg_gpu_send_arr(GPU_ptr dest, const void *src, size_t mem_size, sg_stream_ptr stream)
{
 
  

 
  fct_call_trsf_CPU_GPU *trsfCPU_GPU = 
    new fct_call_trsf_CPU_GPU(src,dest,mem_size);


  ((gpu_stream*)stream)->addOp(trsfCPU_GPU,TRANSF);

  return 0;
}


extern "C"
int sg_gpu_receiv_arr(void *dest, const GPU_ptr src, size_t mem_size, sg_stream_ptr stream)
{

  fct_call_trsf_GPU_CPU *trsfGPU_CPU = 
    new fct_call_trsf_GPU_CPU(src,dest,mem_size);
 

  ((gpu_stream*)stream)->addOp(trsfGPU_CPU,TRANSF);


  return 0;
}


extern "C"
int sg_mem_copy(void *dest, const void *src, size_t mem_size, sg_stream_ptr stream)
{
fct_call_memcpy *trsf_memcpy =
  new fct_call_memcpy(src,dest,mem_size);

 
 ((gpu_stream*)stream)->addOp(trsf_memcpy,TRANSF);

 return 0;
}


extern "C"
int sg_calc(sg_callback_ptr f_call, void *param, size_t param_size,sg_stream_ptr stream)
{
  
  fct_call_calc_generic *calc_generic =
    new fct_call_calc_generic(f_call,param,param_size);

  
  ((gpu_stream*)stream)->addOp(calc_generic,CALC);
  return 0;
}
