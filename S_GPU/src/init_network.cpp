/**
 * @file init_network.cpp
 * @author Matthieu Ospici
 * 
 * @brief
 * C++ classes used for initializing the GPU sharing system
 * 
 * @section LICENSE
 * 
 * Copyright (C) 2010 BULL LIG CEA-INAC UJF
 *
 * This file is part of S_GPU library.
 * 
 * S_GPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * S_GPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with S_GPU.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG
#include <config.h>
#endif

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

void local_network::init ( manage_cpu_affinity& mca, const set_repartition* set_r, int _iproc ) throw ( inter_node_communication_error, check_calc_error )
{
	iproc = _iproc;
	// man_gpu = NULL;
	sem_unix_gpu_TRSF = NULL;
	sem_unix_gpu_CALC = NULL;
	const int SIZE_NAME = 100;
	const int SIZE_FICH_NAME = 40;
	char nomfic[SIZE_FICH_NAME];
	int participants[NUM_PARTICIPANTS];
	FILE * es = NULL;
	int lu = 0;
	int n, lectOK;
	char curr_sock_path[SIZE_NAME];
	struct  sockaddr_un servaddr; /* server adress (current address) */
	char to_sock_path[SIZE_NAME];
	//  struct  sockaddr_un toaddr; /* neigbour adress */

	bzero ( &servaddr, sizeof ( servaddr ) );
	bzero ( &toaddr, sizeof ( servaddr ) );
	try
	{
		if ( ( sock = socket ( AF_UNIX, SOCK_DGRAM, 0 ) ) < 0 )
		{
			sock = 0;
			throw inter_node_communication_error ( "server: socket" );
		}
		snprintf ( curr_sock_path, SIZE_NAME, "/tmp/sock%i", getpid() );
		unlink ( curr_sock_path );

		servaddr.sun_family = AF_UNIX;
		strncpy ( servaddr.sun_path, curr_sock_path, SIZE_NAME );

		if ( bind ( sock, ( struct sockaddr* ) &servaddr, sizeof ( servaddr ) ) < 0 )
		{
			throw inter_node_communication_error ( "server: bind" );
		}

		strncpy ( nomfic, "/tmp/yo.yo", SIZE_FICH_NAME );
		es = fopen ( nomfic, "a+" );

		//put lock on file

		if ( lockf ( fileno ( es ), F_LOCK, 0 ) < 0 )
			printf ( "** WARNING ** lockf take ERROR\n" );

		fprintf ( es, "%i\n", getpid() );

		lockf ( fileno ( es ), F_ULOCK, 0 );
		std::cout << "Waiting for registration of all process..." << std::endl;
		//  printf( "%i\n", getpid());
		while ( lu < NUM_PARTICIPANTS )
		{
			lu = 0;
			sleep ( 1 );
			rewind ( es );
			do
			{
				lectOK = fscanf ( es, "%i", &n );

				if ( lectOK == 1 )
				{
					if ( lu >= NUM_PARTICIPANTS )
					{
						throw inter_node_communication_error ( "Too much MPI process launched on node" );
					}
					participants[lu] = n;
					if ( n == getpid() )
					{
						currNum = lu;
					}
					++lu;
				}
			}
			while ( lectOK == 1 && fgetc ( es ) != EOF );
		}
		//  std::cout << "Process  registration done !" << std::endl;
		voisin = ( currNum + 1 ) % NUM_PARTICIPANTS;
		//ok we have the network topology !
		snprintf ( to_sock_path, SIZE_NAME, "/tmp/sock%i", participants[voisin] );

		toaddr.sun_family = AF_UNIX;
		strncpy ( toaddr.sun_path, to_sock_path, SIZE_NAME );

		//set repartition (static, shared...)
		set_r->do_repartition ( currNum, &currGPU );

		//create numGPU semaphore (only the process 0)

		if ( currNum == 0 )
		{
			//create sems and send to all
			sem_unix_gpu_TRSF = new sem_unix ( 2*NUM_GPU, 2*currGPU ); //2sem per GPUs
			sem_unix_gpu_CALC = new sem_unix();
			sem_unix_gpu_CALC->createFromExistingSem ( sem_unix_gpu_TRSF->getSemid(), 2*currGPU + 1 );
			message msgSem;
			msgSem.node = sem_unix_gpu_TRSF->getSemid();
			send_next ( &msgSem );
			recv_prev ( &msgSem );
			std::cout << "OK, all process has semaphores : " << msgSem.node << std::endl;
		}
		else
		{
			message msgSemRcv;
			recv_prev ( &msgSemRcv );
			sem_unix_gpu_TRSF = new sem_unix();
			sem_unix_gpu_TRSF->createFromExistingSem ( msgSemRcv.node, 2*currGPU );
			sem_unix_gpu_CALC = new sem_unix();
			sem_unix_gpu_CALC->createFromExistingSem ( msgSemRcv.node, 2*currGPU + 1 );
			send_next ( &msgSemRcv );
		}
		fclose ( es );
		remove ( nomfic );
		es = NULL;
		//set affinity
		//CPU
		mca.set_affinity ( currNum );
	}
	catch ( ... )
	{
		if ( sock !=  0 )
			close ( sock );
		unlink ( toaddr.sun_path );
		if ( es != NULL )
		{
			fclose ( es );
			remove ( nomfic );
		}
		if ( sem_unix_gpu_TRSF != NULL )
			delete sem_unix_gpu_TRSF;
		if ( sem_unix_gpu_CALC != NULL )
			delete sem_unix_gpu_CALC;
		//  if(man_gpu != NULL)
		//	delete man_gpu;
		throw;
	}
}

local_network::~local_network()
{
	close ( sock );
	unlink ( toaddr.sun_path );
	//  delete man_gpu;
	delete sem_unix_gpu_CALC;
	delete sem_unix_gpu_TRSF;
}

int local_network::send_next ( const message *msg ) throw ( inter_node_communication_error )
{
	if ( sendto ( sock, msg, sizeof ( message ), 0, ( struct sockaddr* ) &toaddr, sizeof ( toaddr ) ) < 0 )
	{
		std::ostringstream oss;
		oss << "sendto(), msg : " << msg->msgType;

		throw inter_node_communication_error ( oss.str() );
		// return errno;
	}
	return 0;
}

int local_network::recv_prev ( message *msg ) throw ( inter_node_communication_error )
{
	int n;
	socklen_t tosize;
	if ( ( n = recvfrom ( sock, msg, sizeof ( message ), 0, NULL, &tosize ) ) < 0 )
	{
		throw inter_node_communication_error ( "recvfrom()" );
	}
	return n;
}

void local_network::messageLoopNetwork ( localqueu& locq ) throw ( inter_node_communication_error )
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
	bool packed = false;
	do
	{
		//  std::cout << "Iteration !!  " << std::endl;
		if ( locq.isFinish() )
			break; //finish
		fct = locq.gpu_dispo ( currGPU, PACK );
		if ( fct != NULL )
		{
			//we have a packed operation !!
			// next operation must take sem, and the other don't have to take it while onother packed operation is seen
			if ( !packed )
			{
				//take sem
				sem_unix_gpu_TRSF->P();
				packed = true;
			}
			else
			{
				sem_unix_gpu_TRSF->V();
				packed = false;
			}
		}

		fct = locq.gpu_dispo ( currGPU, TRANSF );
		if ( fct != NULL )
		{
			//we must do a TRSF
			doTRSF ( fct, packed );
		}

		fct = locq.gpu_dispo ( currGPU, CALC );

		if ( fct != NULL )
		{
			//we must do a CALC
			doCALC ( fct, packed );
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
	while ( true );
	//std::cout << "END OF DO  " << std::endl;
}

void   local_network::doTRSF ( fct_call *fct, bool packed )
{
	if ( !packed )
		sem_unix_gpu_TRSF->P();

	//call the function (operator() overloaded)
	( *fct ) ( currGPU );

	if ( !packed )
		sem_unix_gpu_TRSF->V();
}

void   local_network::doCALC ( fct_call *fct, bool packed )
{
	if ( !packed )
		sem_unix_gpu_CALC->P();

	//call the function (operator() overloaded)
	( *fct ) ( currGPU );

	if ( !packed )
		sem_unix_gpu_CALC->V();
}

//------------ manage end ----------------

manage_end::manage_end ( int num )
{
	part_matrix = new std::vector<bool> ( num, false );
}

manage_end::~manage_end()
{
	delete part_matrix;
}

void manage_end::setPartFinish ( int numPart )
{
	part_matrix->at ( numPart ) = true;
}

bool manage_end::allFinish() const
{
	std::vector<bool>::const_iterator cit;
	for ( cit = part_matrix->begin(); cit < part_matrix->end() ; ++cit )
	{
		if ( ! ( *cit ) )
			return false;
	}
	return true;
}

//----------------------- sem unix

sem_unix::sem_unix ( int numGPU, int currGPU_ ) throw ( synchronization_error )
{
	currGPU = currGPU_;
	//  const int semKEY = 33;
	union semun
	{
		int              val;    /* Value for SETVAL */
		struct semid_ds *buf;    /* Buffer for IPC_STAT, IPC_SET */
		unsigned short  *array;  /* Array for GETALL, SETALL */
	} arg_ctl;

	bzero ( &arg_ctl, sizeof ( semun ) );

	//  semid = semget ( ftok (nomfic, semKEY), numGPU, IPC_CREAT | IPC_EXCL | 0666);
	semid = semget ( IPC_PRIVATE, numGPU, IPC_CREAT | IPC_EXCL | 0666 );

	if ( semid == -1 )
	{
		perror ( "Error init sem" );
		const int ERR_SIZE = 1024;
		char errorm[ERR_SIZE];
		strerror_r ( errno, errorm, ERR_SIZE );

		semctl ( semid , 0 , IPC_RMID , 0 ) ;
		throw synchronization_error ( errorm );
	}

	arg_ctl.val = 1;

	for ( int i = 0; i < numGPU; ++i )
	{
		if ( semctl ( semid, i, SETVAL, arg_ctl ) == -1 )
		{
			semctl ( semid , 0 , IPC_RMID , 0 ) ;
			throw synchronization_error ( "Semctl ERROR" );
		}
	}
	numgpu = numGPU;
	initHere = true;
}

sem_unix::sem_unix()
{
	initHere = false;
}

void sem_unix::createFromExistingSem ( int _semid, int _currGPU )
{
	semid = _semid;
	currGPU = _currGPU;
}

sem_unix::~sem_unix()
{
	if ( initHere )
	{
		if ( semctl ( semid , 0 , IPC_RMID , 0 ) < 0 )
		{
			std::cout << "Semaphore destruction error" << std::endl;
		}
		std::cout << "Remove sem : "  << semid << std::endl;
	}
}

void sem_unix::P() throw ( synchronization_error )
{
	struct sembuf sempar;
	bzero ( &sempar, sizeof ( sembuf ) );
	sempar.sem_num = currGPU;
	sempar.sem_op = -1;
	sempar.sem_flg = SEM_UNDO;

	if ( semop ( semid , &sempar , 1 ) < 0 )
	{
		perror ( "semop P()" );
		std::cerr << "sem-num " << currGPU << std::endl;
		throw synchronization_error ( "Semop P() error" );
	}
}

/********************************************************/

void sem_unix::V() throw ( synchronization_error )
{
	struct sembuf sempar;
	bzero ( &sempar, sizeof ( sembuf ) );
	sempar.sem_num = currGPU;
	sempar.sem_op = 1;
	sempar.sem_flg = SEM_UNDO;

	if ( semop ( semid , &sempar , 1 ) < 0 )
	{
		throw synchronization_error ( "Semop V() error" );
	}
}
