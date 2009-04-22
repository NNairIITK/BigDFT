#ifndef LOCALQUEU_H
#define LOCALQUEU_H

#include <vector>

#include "fct_call.h"

enum gpuOp {CALC, TRANSF, FINISH};

class gpu_stream
{
public:
  gpu_stream(int gpu_id_)
    :gpu_id(gpu_id_),is_gpu_id_set(true)//is setted...sni
  {}
  ~gpu_stream();

  void addOp(fct_call*, gpuOp op);
  fct_call* getNextFct_call();


  int getGpu_id() const {return gpu_id;}
  void setGpu_id(int id);

  bool getIs_gpu_id_set() const {return is_gpu_id_set;}

  gpuOp getOp_to_do() const;
 

private:
  struct struct_op_assoc
  {
    struct_op_assoc(fct_call* fct_c_,gpuOp op_)
      :fct_c(fct_c_),op(op_)
    {}
    fct_call* fct_c;
    gpuOp op;
  };

  int gpu_id; //computed on the gpu number gpu_id
  bool is_gpu_id_set;
  std::vector<struct_op_assoc*> v_stream;
 

 
  
  std::vector<struct_op_assoc*>::iterator next_step_it;
  
};

class localqueu
{
public:
  localqueu();

  bool isFinish() const;

  void addStream(gpu_stream* stream);
  void  removeStreams();
  fct_call* gpu_dispo(int gpu_dipo, gpuOp op);
private:

  
  std::vector<gpu_stream*> v_stream_queu;

  

};


#endif
