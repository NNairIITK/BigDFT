/**
 * @file set_repartition.cpp
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

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include <vector>
#include <iostream>

#include "set_repartition.h"
#include "CUDA/cudafct.h"

void set_repartition_shared::do_repartition ( int currNum, int *connectedGPU ) const
{
	if ( iproc == 0 )
		std::cout << "Shared repartition" << std::endl;
  //int currGPU;
	//compute which card has each participant
	std::vector<int> tab ( NUM_ACTORS, 0 );
	const int div = NUM_ACTORS / NUM_GPU;
	int numIt = div;
	int itGPU = 0;

	for ( int i = 0; i < NUM_ACTORS; ++i )
	{
		tab.at ( i ) = itGPU;
		if ( iproc == 0 )
			std::cout << "Unix process (not MPI) " << i << " has GPU : " << itGPU << std::endl;
		if ( i + 1 == numIt )
		{
			numIt += div;
			++itGPU;
		}
	}

	*connectedGPU = tab.at ( currNum );
	c_cuda_setdevice ( *connectedGPU );
	gga->setIsAttached ( true );
	//set GPU
}

void set_repartition_static::do_repartition ( int currNum, int *connectedGPU ) const
{
	if ( iproc == 0 )
		std::cout << "Static repartition" << std::endl;
	if ( currNum < NUM_GPU )
	{
		std::cout << "Unix process (not MPI) " << currNum << " has GPU : " << currNum << std::endl;
		c_cuda_setdevice ( currNum );
		gga->setIsAttached ( true );
		*connectedGPU = currNum;
	}
}
