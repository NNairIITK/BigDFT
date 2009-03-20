#ifndef __exceptionh__
#define __exceptionh__

#include <string>
#include <exception>
 
class network_error: public std::exception
{
public:
    network_error(std::string _msg) throw()
         :msg(_msg)
    {}
 
     virtual const char* what() const throw()
     {
         return msg.c_str();
     }
     
 
    
    virtual ~network_error() throw()
    {}
 
private:
  
    std::string msg;       //Description de l'erreur
  
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
  
    std::string msg;       //Description de l'erreur
  
};
#endif
