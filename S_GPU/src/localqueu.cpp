/**
 * @file localqueu.cpp
 * @author Matthieu Ospici
 * 
 * @brief
 * Two classes are defined here:
 *   - gpu_stream: handle one stream for a GPU. As a reminder, a stream is a following of dependant operations
 *   which are executed on the GPU.
 *   - localqueu: the list containing all the streams on a node.
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

#include <iostream>
#include <algorithm>

#include "localqueu.h"
#include "exceptions.h"
#include "class_utils.h"

void gpu_stream::setGpu_id ( int id )
{
	if ( is_gpu_id_set == true )
		throw other_error ( "gpu_already_set" );
	else
	{
		gpu_id = id;
		is_gpu_id_set = true;
	}
}

fct_call* gpu_stream::getNextFct_call()
{
	if ( next_step_it != v_stream.end() )
	{
		return ( * ( next_step_it++ ) )->fct_c;
		//  ++next_step_it;
	}
	else return NULL;
}

gpuOp gpu_stream::getOp_to_do() const
{
	// std::cout << "get op nb elem " << v_stream.size() << std::endl;
	if ( next_step_it != v_stream.end() )
	{
		return ( *next_step_it )->op;
	}
	else
		return FINISH;
}


void gpu_stream::addOp ( fct_call* gs, gpuOp op )
{
	struct_op_assoc* t = new struct_op_assoc ( gs, op );
	v_stream.push_back ( t );
	next_step_it = v_stream.begin();
}
gpu_stream::~gpu_stream()
{
	std::for_each ( v_stream.begin(), v_stream.end(), deleter() );
}

//--------- local queue -------------

localqueu::localqueu()
{
}


void localqueu::addStream ( gpu_stream* stream )
{
	v_stream_queu.push_back ( stream );
}

void localqueu::removeStreams()
{
	std::for_each ( v_stream_queu.begin(), v_stream_queu.end(), deleter() );
	v_stream_queu.erase ( v_stream_queu.begin(), v_stream_queu.end() );
}

fct_call* localqueu::gpu_dispo ( int gpu_dipo, gpuOp op )
{
	// std::cout << "beg gpu_dispo " << v_stream_queu.size()  << std::endl;
	std::vector<gpu_stream*>::iterator it;
	for ( it = v_stream_queu.begin() ; it < v_stream_queu.end(); it++ )
	{
		//  std::cout << "it " << (*it)->getOp_to_do()  << std::endl;
		if ( ( *it )->getIs_gpu_id_set() == false && ( *it )->getOp_to_do() == op )
		{
			//ok, we can compute this stream on gpu and set a gpu to the stream
			//	  std::cout << "First case !!" << std::endl;
			( *it )->setGpu_id ( gpu_dipo );
			return  ( *it )->getNextFct_call();
		}
		else if ( ( *it )->getGpu_id() == gpu_dipo && ( *it )->getOp_to_do() == op  )
		{
			//	    std::cout << "Second case !!" << std::endl;
			//ok, we can compute this stream on gpu
			return ( *it )->getNextFct_call();
		}
	}
	return NULL; //not find...
}


bool localqueu::isFinish() const
{
	std::vector<gpu_stream*>::const_iterator cit;
	for ( cit = v_stream_queu.begin() ; cit < v_stream_queu.end(); cit++ )
	{
		if ( ( *cit )->getOp_to_do()  !=  FINISH )
			return false;
	}
	return true;
}
