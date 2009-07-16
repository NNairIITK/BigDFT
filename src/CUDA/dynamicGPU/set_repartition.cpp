#include <vector>
#include <iostream>

#include "set_repartition.h"
#include "cudafct.h"

void set_repartition_shared::do_repartition(int currNum) const
{

  if(iproc == 0)
    std::cout << "shared repartition" << std::endl;
  int currGPU;

  //compute which card has each participant
  std::vector<int> tab(NUM_ACTORS,0);
  

  const int div = NUM_ACTORS / NUM_GPU;

  int numIt = div;
  int itGPU = 0;
  for(int i=0; i<NUM_ACTORS; ++i)
    {
      tab.at(i) = itGPU;
      
      if(iproc == 0)				
	std::cout << "Unix process (not MPI) " << i << " has GPU : " << itGPU << std::endl;
      
      if(i + 1 == numIt)
	{
	  numIt += div;
	  ++itGPU;
	}
      
    }
  
	

 currGPU = tab.at(currNum);
 c_cuda_setdevice(currGPU);


 gga->setIsAttached(true);
 //set GPU

}
void set_repartition_static::do_repartition(int currNum) const
{
  if(iproc == 0)
    std::cout << "static repartition" << std::endl;


  if(currNum < NUM_GPU)
    {

      std::cout << "Unix process (not MPI) " << currNum << " has GPU : " << currNum << std::endl;
     c_cuda_setdevice(currNum);

     gga->setIsAttached(true);
    }
}
