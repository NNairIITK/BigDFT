/*
** fct_call.h
** 
** Made by Matthieu Ospici
** Login   <mo219174@badiane>
** 
** Started on  Thu Apr 17 13:32:28 2008 Matthieu Ospici
** Last update Thu Apr 17 13:32:28 2008 Matthieu Ospici
*/

#ifndef   	FCT_CALL_H_
# define   	FCT_CALL_H_


class fct_call
{
public:
  virtual void operator()(int)  = 0; 
  virtual ~fct_call(){};
};



class fct_call_do_nothing : public fct_call
{
public:
  virtual void operator()(int) {}; 
  virtual ~fct_call_do_nothing(){};
};

#endif 	    /* !FCT_CALL_H_ */
