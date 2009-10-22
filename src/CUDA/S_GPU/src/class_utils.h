#ifndef CLASS_UTILS_H
#define CLASS_UTILS_H

#include <sstream>

class deleter 
{ 
public:
   template <class T> void operator ()(T*& p) const 
   { 
      delete p;
      p = NULL;
   } 
}; 


inline void getHostName(std::string& h)
{
  const int HOST_NAME_SIZE = 300;
  char hostname[HOST_NAME_SIZE];
  gethostname(hostname,HOST_NAME_SIZE);

  h = hostname;
}



template<typename T>
T strTo(const std::string& str)
{
  T dest;
  std::istringstream iss( str );
  iss >> dest;
  return dest;
}




#endif
