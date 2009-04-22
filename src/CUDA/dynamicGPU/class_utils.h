#ifndef CLASS_UTILS_H
#define CLASS_UTILS_H

class deleter 
{ 
public:
   template <class T> void operator ()(T*& p) const 
   { 
      delete p;
      p = NULL;
   } 
}; 


#endif
