/****h* fct_call.h/fct_call
*  DESCRIPTION
*   fct_call is a class designed to be extended in order to create new code for GPU
*   the virtual operator() must be extended
* AUTHOR
*   Matthieu Ospici
*  SOURCE
*/


#ifndef   	FCT_CALL_H_
# define   	FCT_CALL_H_


class fct_call
{
public:
  virtual void operator()(int)  = 0; 
  virtual ~fct_call(){};
};



#endif 	    /* !FCT_CALL_H_ */
