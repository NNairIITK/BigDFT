#include <string>
#include "checker.h"
#include "../class_utils.h"
#include "../exceptions.h"
//this function is on check_card/check_init.f90
extern "C"
void check_init__(int *error);


 void checker::runTestOne() throw (check_calc_error)
{
  int ret_error;
  check_init__(&ret_error);
  
  if(ret_error != 0)
    {
      
      std::string hostname;
      getHostName(hostname);
      throw check_calc_error(hostname);
      
    }
}
