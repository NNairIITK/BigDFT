#ifndef __messageh__
#define __messageh__

enum msgType_t {mCALC, mTRANSF, mSTOP, mSEM, mTERMINATE };

struct message
{
  msgType_t msgType;
 
  int gpu_avail; //this token means gpu_avail is available
  
  // unsigned int memAvaillable;
  
  int node; //for mSTOP, node has finish, or mSEM, semID...
 

};


#endif
