#ifndef __read_conf_exceptionh__
#define __read_conf_exceptionh__

#include <string>
#include <exception>
 
class file_not_found: public std::exception
{
public:
   file_not_found(std::string _msg) throw()
         :msg(_msg)
    {}
 
     virtual const char* what() const throw()
     {
         return msg.c_str();
     }
     
 
    
    virtual ~file_not_found() throw()
    {}
 
private:
  
    std::string msg;       //Error description
  
};



class read_not_found: public std::exception
{
public:
  read_not_found(std::string _msg) throw()
    :msg(_msg)
  {}
 
  virtual const char* what() const throw()
  {
    return msg.c_str();
  }
     
 
    
    virtual ~read_not_found() throw()
    {}
 
private:
  
    std::string msg;       //Error description
  
};



class read_not_found_GPU : public read_not_found
{
public:
  read_not_found_GPU(std::string _msg) throw()
    :read_not_found(_msg){}

  virtual ~read_not_found_GPU() throw()
  {}
};

class read_not_found_CPU : public read_not_found
{
public:
read_not_found_CPU(std::string _msg) throw()
    :read_not_found(_msg){}

  virtual ~read_not_found_CPU() throw()
  {}
};

#endif
