#ifndef __checkcuda__
#define __checkcuda__

#include <string>
#include <sstream>
#include "class_utils.h"

/*
inline void getHostName(std::string& h)
{
  const int HOST_NAME_SIZE = 300;
  char hostname[HOST_NAME_SIZE];
  gethostname(hostname,HOST_NAME_SIZE);

  h = hostname;
}
*/

template<class E>
inline void check(bool i,const std::string& st,const char* file,int line)
{
  if(i)
    {
      //   std::string line(__LINE_);
      //    std::string file(__FILE__);
      std::ostringstream oss;

      std::string host;
      getHostName(host);

      oss << "ERROR : " << st << std::endl << "host : " << host << ", file : " << file << ", line : " << line << std::endl;
      throw E(oss.str());
    }
}

template<class E>
inline void check(bool i,const std::string& st)
{
  if(i)
    {
      //   std::string line(__LINE_);
      //    std::string file(__FILE__);
      std::ostringstream oss;

      std::string host;
      getHostName(host);

      oss << "ERROR : " << st << std::endl << "host : " << host  << std::endl;
      throw E(oss.str());
    }
}




template<class E>
inline void check_cuda_error()
{
  cudaError_t error = cudaGetLastError();


  check<E>(error != cudaSuccess,cudaGetErrorString(error));

}


class cuda_error: public std::exception
{
 public:
  cuda_error(std::string _msg) throw()
    :msg(_msg)
    {}

    virtual const char* what() const throw()
    {
      return msg.c_str();
    }



    virtual ~cuda_error() throw()
      {}

 private:

    std::string msg;       //Error description

};




#endif
