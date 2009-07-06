#include <iostream>
#include <exception>
#include <string>
#include <sstream>


#include <sched.h>
#include "read_conf_file.h"

#include "trace_exec.h"

trace_exec *tracer;

int c_cudaSetDevice(int device);
void init_gpu_sharing(readConfFile&,int iproc,int *error);

const char *NAME_FILE = "GPU.config";

void set_cpu_gpu_aff(int iproc, int *flag, int *flag_blas_conv, int *error);


extern "C"
void init_lib__(int *iproc,int *error, int *iconv, int *iblas, bool * GPUshare)
{
  
 
  try
    {
      //Critical configuration information, if one is missing, exit
      readConfFile read_conf(NAME_FILE);
      int use_shared;
      
      read_conf.get("USE_SHARED",&use_shared);
      
 

      if(use_shared == 1) 
	{
	  if(*iproc == 0)
	    std::cout << "** GPU SHARING ENABLED" << std::endl;

	  init_gpu_sharing(read_conf,*iproc,error);
	  *GPUshare = true;
	  *iconv = 0;
	  *iblas = 0; //enable GPU convolution and GPU blas
	}
      else
	{
	  if(*iproc == 0)
	    std::cout << "** GPU SHARING *DISABLED*" << std::endl;

	  set_cpu_gpu_aff(*iproc, iconv, iblas, error);
	

	  *GPUshare = false;
	}
    }
  catch(read_not_found re)
    {
      std::cerr << "*** ERROR : INVALID CONFIG FILE. You have to set USE_SHARED to 1 or 0" << std::endl;
      std::cerr << "Missing information : " << re.what() << std::endl;
      *error = 1;
    }

  catch(file_not_found e)
    {
      std::cerr<< "**ERROR GPU configuration  file not found : " <<  e.what() << std::endl;
      *error = 1;

    }

  catch(...)
    {
      std::cerr<< "** Unexpected exception "<< std::endl;
      *error = 1;

    }




  std::ostringstream ostr;
  ostr << "trace_" << *iproc;
  tracer = new trace_exec(ostr.str(),false);


}


void set_cpu_gpu_aff(int iproc, int *flag, int *flag_blas_conv, int *error)
{
  const std::string NameFile(NAME_FILE);


  try
    {
      readConfFileGPU_CPU r(NameFile);

      *flag = 1;

      *flag_blas_conv = r.getFlag(iproc);


      int aff_GPU_ID = r.getGPU(iproc);
      int aff_CPU_ID = r.getCPU(iproc);
   
      *flag = 0; //OK, we can use GPU for this MPI_ID because no execptions
      std::cout << "ID : " << iproc << " CPU_ID : " << aff_CPU_ID << " GPU_ID : " << aff_GPU_ID << " FLAG BLAS = " << *flag_blas_conv << std::endl;

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
      std::cout << "ID : " << iproc << "UNDEF GPU, so we disable affinity" << std::endl;
    }
  catch(file_not_found e)
    {
      std::cerr<< "**ERROR File not found : " <<  NameFile << std::endl;
      //  return;
    }
  catch(read_not_found_CPU e)
    {
      std::cout << "ID : " << iproc << "UNDEF CPU , so we disable affinity " << std::endl;
    }

  catch(read_not_found e)
    {
      std::cout << "STRAGE STRANGE " << std::endl;
    }
  
  


}





