#include <string>
#include <iostream>
#include <sstream>

#include <errno.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/socket.h>
#include <stdlib.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <stdio.h>
#include <signal.h>
#include <netdb.h>
#include <errno.h>
#include <sys/ipc.h>
#include <sys/sem.h>



#include "init_network.h"
#include "message.h"
#include "localqueu.h"
#include "manage_gpu.h"
#include "cudafct.h"
void local_network::init() throw (network_error)
{

  man_gpu = NULL;
  sem_unix_gpu_TRSF = NULL;
  sem_unix_gpu_CALC = NULL;


  const int SIZE_NAME = 100;
  const int SIZE_FICH_NAME = 40;
  
  char nomfic[SIZE_FICH_NAME];
  
  
  
  int participants[NUM_PARTICIPANTS];
  FILE * es = NULL;
  int lu = 0;
      
  int n,lectOK;
  
  
  char curr_sock_path[SIZE_NAME];
  struct  sockaddr_un servaddr; /* server adress (current address) */
  
  char to_sock_path[SIZE_NAME];
  //  struct  sockaddr_un toaddr; /* neigbour adress */
  
  bzero(&servaddr, sizeof(servaddr));
  bzero(&toaddr, sizeof(servaddr));
  
  
  
  try
    {    
      if ((sock = socket(AF_UNIX, SOCK_DGRAM, 0)) < 0) 
	{
	  sock = 0;
	  throw network_error("server: socket");
      
	}
      
      
      snprintf(curr_sock_path,SIZE_NAME,"/tmp/sock%i",getpid());
      unlink(curr_sock_path);
      
      servaddr.sun_family = AF_UNIX;
      strncpy(servaddr.sun_path, curr_sock_path,SIZE_NAME);
  
      if (bind(sock, (struct sockaddr*)&servaddr, sizeof(servaddr)) < 0) 
	{
	  throw network_error("server: socket");
	} 
      
      
      
      
      
      strncpy(nomfic, "/tmp/yo.yo", SIZE_FICH_NAME);
      es = fopen(nomfic, "a+");
      fprintf(es, "%i\n", getpid());
      //  printf( "%i\n", getpid());
      while(lu < NUM_PARTICIPANTS)
	{
	  lu = 0;
      
	  sleep(1);
	  rewind(es);
	  
	  do
	    {
	      lectOK = fscanf(es, "%i", &n);
	      if (lectOK == 1) 
		{
		  participants[lu] = n;
		  if(n == getpid())
		    {
		      currNum = lu;
		    }
		  ++lu;
		  // printf("%i-- %i\n", getpid(), n);
		}
	    }
	  while (lectOK == 1 && fgetc(es) != EOF);
	}
   
      
      
      voisin = (currNum + 1)%NUM_PARTICIPANTS;
      
      
      
      //ok we have the network topology !
      
      snprintf(to_sock_path,SIZE_NAME,"/tmp/sock%i",participants[voisin]);
      
      toaddr.sun_family = AF_UNIX;
      strncpy(toaddr.sun_path, to_sock_path,SIZE_NAME);
      

      //compute which card has each participant
      std::vector<int> tab(NUM_PARTICIPANTS,0);
      //    int tab[NUM_PARTICIPANTS];
      
      int div = NUM_PARTICIPANTS / NUM_GPU;
      int v=-1;
        for(int i=0;i<NUM_PARTICIPANTS;i+=div)
	{
	  ++v;
	    for(int j=0;j<div;++j)
	    {
	         tab.at(j+i) = v;
	      //    tab.at(0) = 0;
	      //   tab.at(1) = 0;
	      //   tab.at(2) = 0;
	        std::cout << "CURR GPU " << j+i << " " << v << std::endl;

	      // std::cout << "CURR GPU " << "O" << " " << tab.at(0) << std::endl;
	      // std::cout << "CURR GPU " << "1" << " " << tab.at(1) << std::endl;
	    }
	}
      
      currGPU = tab.at(getCurr());

      //	std::cout << "Cur GPu ù*** " << currGPU << std::endl;
      // currGPU = 1;//tab[getCurr()];



     

      // c_cuda_setdevice(currGPU);
      //create numGPU semaphore (only the process 0)

      if(currNum == 0)
	{
	  //create sems and send to all
	  sem_unix_gpu_TRSF = new sem_unix(nomfic,2*NUM_GPU,2*currGPU); //2sem per GPUs

	  sem_unix_gpu_CALC = new sem_unix(sem_unix_gpu_TRSF->getSemid(),2*currGPU + 1);

	  std::cout << "creat " << 2*NUM_GPU << "sem, TRSF : " <<2*currGPU<< " CALC : " << 2*currGPU + 1 << std::endl;

	  message msgSem;

	
	  msgSem.node = sem_unix_gpu_TRSF->getSemid();
	  std::cout << "On envoie !" << msgSem.node << std::endl;
	  send_next(&msgSem);
	
	  
	  recv_prev(&msgSem);
	  std::cout << "OK message bouclé" << msgSem.node << std::endl;
	}
      else
	{
	  message msgSemRcv;
	  recv_prev(&msgSemRcv);
	  sem_unix_gpu_TRSF = new sem_unix(msgSemRcv.node,2*currGPU);


	  sem_unix_gpu_CALC = new sem_unix(msgSemRcv.node,2*currGPU +1);

	  std::cout << "RECU,  sem, TRSF : " <<2*currGPU<< " CALC : " << 2*currGPU + 1 << std::endl;

	  std::cout << "On a recu , curGPU" << currGPU << "  curndode " << msgSemRcv.node << std::endl;

	  std::cout << "On re envoie !" << msgSemRcv.node << std::endl;
	  send_next(&msgSemRcv);
	  }
      fclose(es);
      remove(nomfic);

      c_cuda_setdevice(currGPU);
      //    man_gpu = new manage_gpu(currGPU); //snif...
      // man_gpu = new manage_gpu(1); //snif...
    }
  catch(...)
    {
      if(sock !=  0)
	close(sock);

      unlink(toaddr.sun_path);
      if(es != NULL)
	{
	  fclose(es);
	  remove(nomfic);
	}

      if(sem_unix_gpu_TRSF != NULL)
	delete sem_unix_gpu_TRSF;
      
      if(sem_unix_gpu_CALC != NULL)
	delete sem_unix_gpu_CALC;
      //  if(man_gpu != NULL)
      //	delete man_gpu;
      throw;
    }


}




 local_network::~local_network()
{
     close(sock);
     unlink(toaddr.sun_path);
     //  delete man_gpu;
     delete sem_unix_gpu_CALC;
     delete sem_unix_gpu_TRSF;
}
int local_network::send_next(const message *msg) throw (network_error)
{
  if(sendto(sock, msg, sizeof(message), 0, (struct sockaddr*)&toaddr, sizeof(toaddr)) < 0)
    {
      perror("sendto()");

      std::ostringstream oss;
      oss << "sendto(), msg : " << msg->msgType;


      throw network_error(oss.str());
      // return errno;
    }

  return 0;
}

int local_network::recv_prev(message *msg) throw (network_error)
{
  int n;
socklen_t tosize;
 if((n = recvfrom(sock, msg, sizeof(message), 0, NULL, &tosize)) < 0)
    {
      perror("recvfrom()");
      throw network_error("recvfrom()");
      
    }

  return n;
}



void local_network::messageLoopNetwork(localqueu& locq) throw (network_error)
{

  //  manage_end me(NUM_PARTICIPANTS);
 
  //  std::cout << "CUR gpu " <<currGPU<< std::endl;
  //  std::cout << "after me" << NUM_PARTICIPANTS  << std::endl;


  //initials messages (NUM_GPU for CALC and NUM_GPU for TRANSF
  //send only by 0

  /* if(getCurr() == 0)
    {
      message msgCalcInit;
      message msgTransfInit;
      msgCalcInit.msgType = mCALC;
      //    msgCalcInit.initiedBy = 0;
      
      msgTransfInit.msgType = mTRANSF;
      //  msgTransfInit.initiedBy = 0;

      for(int i=0 ; i<NUM_GPU ; ++i)
	{
	  msgCalcInit.gpu_avail = i;
	  msgTransfInit.gpu_avail = i;

	  //    msgCalcInit.gpu_avail = 1;
	  //  msgTransfInit.gpu_avail = 1;

	  // std::cout << "send msg : GPU :" << i << std::endl;
	  send_next(&msgCalcInit);
	  send_next(&msgTransfInit);
	}
      std::cout << "init msg sent : GPU :"  << std::endl;
      }*/


  // message msgRecv;
 

  //  bool mstopSend = false;
  //  bool mtermSend = false;
  fct_call *fct;
  //std::cout << "BEGIN OF DO  " << std::endl;
  do
    {
      //  std::cout << "Iteration !!  " << std::endl;
      if(locq.isFinish())
	break; //finish

      
      fct = locq.gpu_dispo(currGPU,TRANSF);
      if(fct != NULL)
	{
	  //we must do a TRSF
	  doTRSF(fct);
	}

      fct = locq.gpu_dispo(currGPU,CALC);
      if(fct != NULL)
	{
	  //we must do a CALC
	  doCALC(fct);
	}
      
      //  std::cout << "end of one iteration " << std::endl;
      /*  sem_unix_gpu_TRSF->P();

      std::cout << "on tente un TRANSF " << std::endl;
      fct = locq.gpu_dispo(currGPU,TRANSF);
      if(fct == NULL)
	{
	  sem_unix_gpu_TRSF->V();
	  goto bric;
	}
      std::cout << "on FAIT un TRANSF " << std::endl;
      (*fct)(666);
     

      sem_unix_gpu_TRSF->V();

    bric:
      
    if(locq.isFinish())
	break;

    sem_unix_gpu_CALC->P();
    
      std::cout << "on tente un CACL " << std::endl;
      fct = locq.gpu_dispo(currGPU,CALC);
      if(fct == NULL)
	{
	  sem_unix_gpu_CALC->V();
	  continue;
	}
      std::cout << "on FAIT un CALC " << std::endl;
     

      (*fct)(666);
 
      
      
      sem_unix_gpu_CALC->V();*/


      /*    deb:
      recv_prev(&msgRecv);// << std::endl;
      
      if(msgRecv.msgType == mTERMINATE)
	{
	  std::cout << "RECV mterminate " << std::endl;
	  if(getCurr() != 0)
	    send_next(&msgRecv); //send and imediatly quit
	  break;
	}

      
      if(getCurr() == 0 && me.allFinish() && mtermSend==false)
	{
	  //we must terminate all node
	  message msgTerm;
	  msgTerm.msgType = mTERMINATE;
	  std::cout << "send mterminate " << std::endl;
	  send_next(&msgTerm);
	  mtermSend = true;
	  continue;
	}

      if(getCurr() == 0 && me.allFinish())
	continue; //ignore all message because all node are finish
      
      //  std::cout << "recv : " << 
    

      if(!locq.isFinish())
	{
	  //msg recv, so we must do something !
	  fct_call *fct;
	  switch (msgRecv.msgType)
	    {
	    case mCALC:
	      //calc token
	      std::cout << "recv : mCALC " << std::endl;
	      fct = locq.gpu_dispo(msgRecv.gpu_avail,CALC);
	      if(fct == NULL)
		break;
	      
	      std::cout << "recv : mCALC not null " << std::endl;
	      //  (*fct)(getCurrGPU());
	         man_gpu->getGPUcontrol(0).changeThreadOperation(fct); 
	      
	      
	      
	         man_gpu->getGPUcontrol(0).deblockThread();
	      
		 //	 	 goto deb;
	      
			 	 	     man_gpu->getGPUcontrol(0).waitForOneThread();
	      break;
	    case mTRANSF:
	      std::cout << "recv : mTRANSF " << std::endl;
	      fct = locq.gpu_dispo(msgRecv.gpu_avail,TRANSF);
	      if(fct == NULL)
		break;
	      
	      std::cout << "recv : mTRANSF not null " << std::endl;
	      // (*fct)(getCurrGPU());
	        man_gpu->getGPUcontrol(0).changeThreadOperation(fct);
	      
	      
	      
	        man_gpu->getGPUcontrol(0).deblockThread();
	      
		//	goto deb;
		//thread deblocked on fct
			 man_gpu->getGPUcontrol(0).waitForOneThread();
	      break;
	    case mSTOP:
	   
std::cout << "recv1 : mSTOP , " << msgRecv.node<< std::endl;
	      me.setPartFinish(msgRecv.node);
	      if(msgRecv.node == getCurr())
		continue; //mSTOP must bedestroyed
	      break;
	    default:
	      throw network_error("Unknow message");
	    }
	  
	}
      else
	{
	  if(mstopSend == false)
	    {
	      std::cout << "Is finish, send mSTOP" << std::endl;
	      me.setPartFinish(getCurr());
	      
	      message msgStop;
	      msgStop.msgType = mSTOP;
	      msgStop.node = getCurr();
	      send_next(&msgStop);
	      mstopSend = true;
	      //is finish,so we must initiate end
	    }
	}

      if(msgRecv.msgType == mSTOP)
	{
	  std::cout << "recv2 : mSTOP , " << msgRecv.node<< std::endl;
	  me.setPartFinish(msgRecv.node);
	  if(msgRecv.node == getCurr())
	    continue; //mSTOP must bedestroyed
	}

     
      std::cout << "resend " << std::endl;
      send_next(&msgRecv);*/

      //  std::cout << "from : " << msg2.from << ", to " << msg2.to << std::endl;
    }
  while(true );

  //std::cout << "END OF DO  " << std::endl;
}

void   local_network::doTRSF(fct_call *fct)
{
  sem_unix_gpu_TRSF->P();
  
  //  std::cout << "on FAIT un TRANSF !!" << std::endl;
  (*fct)(666);
  

  sem_unix_gpu_TRSF->V();
}
void   local_network::doCALC(fct_call *fct)
{
  sem_unix_gpu_CALC->P();
  
  //  std::cout << "on FAIT un CALC !!" << std::endl;
  (*fct)(666);
  

  sem_unix_gpu_CALC->V();
}

//------------ manage end ----------------

manage_end::manage_end(int num)
{
  part_matrix = new std::vector<bool>(num,false);

}
manage_end::~manage_end()
{
  delete part_matrix;
}
void manage_end::setPartFinish(int numPart)
{
  part_matrix->at(numPart) = true;
}

bool manage_end::allFinish() const
{
  std::vector<bool>::const_iterator cit;
  
  for(cit = part_matrix->begin();cit < part_matrix->end() ; ++cit)
    {
      if(!(*cit))
	return false;
    }
  return true;
}


//----------------------- sem unix


sem_unix::sem_unix(const char *nomfic, int numGPU,int currGPU_) throw (network_error)
 {
   currGPU = currGPU_;

   const int semKEY = 33;
   union semun 
   {
     int              val;    /* Value for SETVAL */
     struct semid_ds *buf;    /* Buffer for IPC_STAT, IPC_SET */
     unsigned short  *array;  /* Array for GETALL, SETALL */
     //  struct seminfo  *__buf;  /* Buffer for IPC_INFO
     //			 (Linux-specific) */
   } arg_ctl;
   
   bzero(&arg_ctl,sizeof(semun));
   /*   union semun 
   {
     int val;
     struct semid_ds *buf;
     ushort *array;
     } arg_ctl;*/
   
	  
  
   semid = semget ( ftok (nomfic, semKEY), numGPU, IPC_CREAT | IPC_EXCL | 0666);
   
   if(semid == -1) 
     {
       perror("semid");
       std::cout << "Erreur semid " <<std::endl;
       semctl (semid , 0 , IPC_RMID , 0) ;
       throw network_error("Erreur semid");
     }

   arg_ctl.val = 1;

   for(int i=0;i<numGPU;++i)
     {
       if(semctl (semid, i, SETVAL, arg_ctl) == -1) 
	 {
	   std::cout << "Erreur semctl gpu :  " << i << std::endl;
	   semctl (semid , 0 , IPC_RMID , 0) ;
	   throw network_error("Erreur semctl");
	
	 }
     }
   numgpu = numGPU;
   initHere = true;
 }


sem_unix::sem_unix(int _semid,int _currGPU)
{
  semid = _semid;
  currGPU = _currGPU;
  initHere = false;
}
sem_unix::~sem_unix()
{
  if(initHere)
    {
      if(semctl (semid , 0 , IPC_RMID , 0) <0)
	{
	  std::cout << "erreur destructeur sem " << std::endl;
	}
      std::cout << "Remove sem : "  << std::endl; 
    }
}
 

void sem_unix::P() throw (network_error)
{
  // std::cout << "Do P() " << semid<< ", id " << currGPU <<  std::endl;
  struct sembuf sempar;
  bzero(&sempar,sizeof(sembuf));
  sempar.sem_num = currGPU;
  sempar.sem_op = -1;
  sempar.sem_flg = SEM_UNDO;
  if(semop (semid , &sempar , 1) < 0)
    {
      perror("P()");
      throw network_error("Erreur semop P()");
    }
  //  std::cout << "FINISH  P() " << semid<< ", id " << currGPU <<  std::endl;
}
/********************************************************/
void sem_unix::V() throw (network_error)
{
  //  std::cout << "Do V() " << semid<< ", id " << currGPU << std::endl;
  struct sembuf sempar;
  bzero(&sempar,sizeof(sembuf));
  sempar.sem_num = currGPU;
  sempar.sem_op = 1;
  sempar.sem_flg = SEM_UNDO;
  if(semop (semid , &sempar , 1) < 0)
    {
      perror("P()");
      throw network_error("Erreur semop V()");
    }
  //  std::cout << "FINISH  V() " << semid<< ", id " << currGPU <<  std::endl;
}
