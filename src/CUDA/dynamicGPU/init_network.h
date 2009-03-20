#ifndef __initnet__
#define __initnet__

#include <sys/socket.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/types.h>    
#include <sys/socket.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <stdio.h>
#include <signal.h>
#include <netdb.h>
#include <errno.h>

#include <string>
#include <vector>

#include "message.h"
#include "exceptions.h"
#include "localqueu.h"
#include "manage_gpu.h"


class sem_unix
{
public:
  sem_unix(int semid,int currGPU); //constructor if we have aleready a sem_id
  sem_unix(const char *nomfic, int numGPU,int currGPU) throw (network_error);
  ~sem_unix();

  void P() throw (network_error);
  void V() throw (network_error);

  int getSemid() const {return semid;};
private:
  int semid;

  int currGPU;
  int numgpu;
  bool initHere;
};

class manage_end
{
public:
  manage_end(int);
  ~manage_end();
  void setPartFinish(int numPart);
  
  bool allFinish() const;
private:
  std::vector<bool> *part_matrix;  
};

//------------------------------------

class local_network
{
public:
  local_network(int num_local_mpi_node, int num_gpu) throw (network_error)
    :NUM_PARTICIPANTS(num_local_mpi_node),NUM_GPU(num_gpu)
  {init();}  

  ~local_network();

  void messageLoopNetwork(localqueu&) throw (network_error);



  inline int getVoisin() const {return voisin;}
  inline int getCurr() const {return currNum;}
  

  inline int getCurrGPU() const {return currGPU;}

  sem_unix *getSemCalc() {return sem_unix_gpu_CALC;}
  sem_unix *getSemTrsf() {return sem_unix_gpu_TRSF;}
  manage_gpu *man_gpu;

  int send_next(const message* msg) throw (network_error);
  int recv_prev(message* msg) throw (network_error);

private:
  void init() throw (network_error);



  void   doTRSF(fct_call *fct);
  void   doCALC(fct_call *fct);
  const int NUM_PARTICIPANTS;

  const int NUM_GPU;
  int sock;
  sockaddr_un toaddr;
  sockaddr_un servaddr;

  int currGPU;
  int currNum;
  int voisin;

  sem_unix *sem_unix_gpu_TRSF;
  sem_unix *sem_unix_gpu_CALC;
  enum op_thread {THREAD_CALC = 0, THREAD_TRSF};
 
};



#endif
