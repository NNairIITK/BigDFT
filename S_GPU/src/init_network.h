/**
 * @file init_network.h
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

#ifndef INIT_NETWORK_H_
#define INIT_NETWORK_H_

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <stdio.h>
#include <signal.h>
#include <netdb.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>

#include <string>
#include <vector>

#include "message.h"
#include "exceptions.h"
#include "localqueu.h"
#include "manage_cpu_affinity.h"
#include "set_repartition.h"

//encapsulation of an unix semaphore
class sem_unix
{
public:
	sem_unix(); //constructor if we have aleready a sem_id
	sem_unix ( int numGPU, int currGPU ) throw ( synchronization_error );
	~sem_unix();
	void P() throw ( synchronization_error );
	void V() throw ( synchronization_error );
	int getSemid() const
	{
		return semid;
	};
	void createFromExistingSem ( int semid, int currGPU );

private:
	int semid;
	int currGPU;
	int numgpu;
	bool initHere;
};

class manage_end
{
public:
	manage_end ( int );
	~manage_end();
	void setPartFinish ( int numPart );
	bool allFinish() const;

private:
	std::vector<bool> *part_matrix;
};

//------------------------------------

class local_network
{
public:
	local_network ( int num_local_mpi_node, int num_gpu, manage_cpu_affinity& mca, const set_repartition* set_r, int iproc ) throw ( inter_node_communication_error, check_calc_error )
			: NUM_PARTICIPANTS ( num_local_mpi_node ), NUM_GPU ( num_gpu )
	{
		init ( mca, set_r, iproc );
	}
	~local_network();
	void messageLoopNetwork ( localqueu& ) throw ( inter_node_communication_error );
	inline int getVoisin() const
	{
		return voisin;
	}
	inline int getCurr() const
	{
		return currNum;
	}
	inline int getCurrGPU() const
	{
		return currGPU;
	}
	sem_unix *getSemCalc()
	{
		return sem_unix_gpu_CALC;
	}
	sem_unix *getSemTrsf()
	{
		return sem_unix_gpu_TRSF;
	}
	// manage_gpu *man_gpu;
	int send_next ( const message* msg ) throw ( inter_node_communication_error );
	int recv_prev ( message* msg ) throw ( inter_node_communication_error );

private:
	void init ( manage_cpu_affinity&, const set_repartition*, int iproc ) throw ( inter_node_communication_error, check_calc_error );
	void   doTRSF ( fct_call *fct, bool );
	void   doCALC ( fct_call *fct, bool );
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
	int iproc;
};
#endif
