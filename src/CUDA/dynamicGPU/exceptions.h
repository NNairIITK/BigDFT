#ifndef __exceptionh__
#define __exceptionh__

#include <string>
#include <exception>
 
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
#endif
