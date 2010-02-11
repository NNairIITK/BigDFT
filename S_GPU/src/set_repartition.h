/**
 * @file set_repartition.h
 * @author Matthieu Ospici
 * 
 * @brief
 * Handle repartition of the GPUs among the processes. Manage both shared and 
 * static repartition.
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

#ifndef SET_REPARTITION_H_
#define SET_REPARTITION_H_

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include "exceptions.h"
#include "manage_global_var.h"

class set_repartition
{
public:
	set_repartition ( int numActors, int numGPU, int iproc_, global_gpu_attach* gga_ )
			: NUM_ACTORS ( numActors ), NUM_GPU ( numGPU ), iproc ( iproc_ ), gga ( gga_ ) {}
	virtual void do_repartition ( int currNum, int *connectedGPU ) const = 0;
	virtual ~set_repartition() {}

protected:
	const int NUM_ACTORS;
	const int NUM_GPU;
	const int iproc;
	global_gpu_attach* gga;
};

class set_repartition_shared : public set_repartition
{
public:
	set_repartition_shared ( int numActors, int numGPU, int iproc, global_gpu_attach* gga_ )
			: set_repartition ( numActors, numGPU, iproc, gga_ ) {}
	virtual void do_repartition ( int currNum, int *connectedGPU ) const;
	virtual ~set_repartition_shared() {}
};

class set_repartition_static : public set_repartition
{
public:
	set_repartition_static ( int numActors, int numGPU, int iproc, global_gpu_attach* gga_ )
			: set_repartition ( numActors, numGPU, iproc, gga_ )
	{}
	virtual void do_repartition ( int currNum, int *connectedGPU ) const;
	virtual ~set_repartition_static() {}
};
#endif
