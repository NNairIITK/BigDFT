#ifndef __cpputilsh__
#define __cpputilsh__

#include <sstream>
template<typename T>
T strTo(const std::string& str)
{
  T dest;
  std::istringstream iss( str );
  iss >> dest;
  return dest;
}

#endif
