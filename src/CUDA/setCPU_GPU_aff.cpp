
#include <iostream>

#include <exception>

#include <sched.h>

#include "read_conf_file.h"

int c_cudaSetDevice(int device);




extern "C" void set_cpu_gpu_aff__(int *iproc, int *flag, int *flag_blas_conv)
{
  const std::string NameFile("cpu_gpu_aff.config");


  try
    {
      readConfFileGPU_CPU r(NameFile);

      *flag = 1;

      *flag_blas_conv = r.getFlag(*iproc);


      int aff_GPU_ID = r.getGPU(*iproc);
      int aff_CPU_ID = r.getCPU(*iproc);
   
      *flag = 0; //OK, we can use GPU for this MPI_ID because no execptions
      std::cout << "ID : " << *iproc << " CPU_ID : " << aff_CPU_ID << " GPU_ID : " << aff_GPU_ID << " FLAG BLAS = " << *flag_blas_conv << std::endl;

      //ok set affinity
      
      cpu_set_t cpus;
      CPU_ZERO(&cpus);
      
      CPU_SET(aff_CPU_ID, &cpus);
      
      //  std::cout << "on set proc - " << CPUaff << std::endl;
      
      if(sched_setaffinity(0, sizeof(cpu_set_t),&cpus) < 0)
	{
	  std::cout << "erreur setaffinity" << std::endl;
	}
      
      if(!CPU_ISSET(aff_CPU_ID, &cpus))
	{
	  std::cout << "*** CPU - " << aff_CPU_ID << " - pas setté ***" << std::endl;

	}
      

      //set GPU affinity
      
      if(c_cudaSetDevice(aff_GPU_ID) != 0)
	std::cout << "erreur affinity cuda" << std::endl;
    }
  catch(read_not_found_GPU e)
    {
      std::cout << "ID : " << *iproc << "UNDEF GPU, so we disable affinity" << std::endl;
    }
  catch(file_not_found e)
    {
      std::cerr<< "**ERROR File not found : " <<  NameFile << std::endl;
      //  return;
    }
  catch(read_not_found_CPU e)
    {
      std::cout << "ID : " << *iproc << "UNDEF CPU , so we disable affinity " << std::endl;
    }

  catch(read_not_found e)
    {
      std::cout << "STRAGE STRANGE " << std::endl;
    }
  

}
