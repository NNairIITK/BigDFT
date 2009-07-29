#ifndef __exceptionh__
#define __exceptionh__

#include <string>
#include <exception>
#include <sstream>  
#include "class_utils.h"


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





class inter_node_communication_error: public std::exception
{
public:
    inter_node_communication_error(std::string _msg) throw()
         :msg(_msg)
    {}
 
     virtual const char* what() const throw()
     {
         return msg.c_str();
     }
     
 
    
    virtual ~inter_node_communication_error() throw()
    {}
 
private:
  
    std::string msg;       //Error description
  
};



class synchronization_error: public std::exception
{
public:
  synchronization_error(std::string _msg) throw()
    :msg(_msg)
  {}
 
  virtual const char* what() const throw()
  {
    return msg.c_str();
  }
     
 
    
    virtual ~synchronization_error() throw()
    {}
 
private:
  
    std::string msg;       //Error description
  
};



class check_calc_error: public std::exception
{
public:
  check_calc_error(std::string _msg) throw()
    :msg(_msg)
  {}
 
  virtual const char* what() const throw()
  {
    return msg.c_str();
  }
     
 
    
    virtual ~check_calc_error() throw()
    {}
 
private:
  
    std::string msg;       //Error description
  
};



class repartition_error: public std::exception
{
public:
  repartition_error(std::string _msg) throw()
    :msg(_msg)
  {}
 
  virtual const char* what() const throw()
  {
    return msg.c_str();
  }
     
 
    
    virtual ~repartition_error() throw()
    {}
 
private:
  
    std::string msg;       //Error description
  
};


class other_error: public std::exception
{
public:
  other_error(std::string _msg) throw()
    :msg(_msg)
  {}
 
  virtual const char* what() const throw()
  {
    return msg.c_str();
  }
     
 
    
    virtual ~other_error() throw()
    {}
 
private:
  
    std::string msg;       //Error description
  
};



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
