/**
 * @file localqueu.h
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

#ifndef LOCALQUEU_H_
#define LOCALQUEU_H_

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include <vector>

#include "fct_call.h"

enum gpuOp {CALC, TRANSF, PACK, FINISH};

class gpu_stream
{
public:
	gpu_stream ( int gpu_id_ )
			: gpu_id ( gpu_id_ ), is_gpu_id_set ( true ) //is setted...sni
	{}
	~gpu_stream();
	void addOp ( fct_call*, gpuOp op );
	fct_call* getNextFct_call();
	int getGpu_id() const
	{
		return gpu_id;
	}
	void setGpu_id ( int id );

	bool getIs_gpu_id_set() const
	{
		return is_gpu_id_set;
	}
	gpuOp getOp_to_do() const;

private:
	struct struct_op_assoc
	{
		struct_op_assoc ( fct_call* fct_c_, gpuOp op_ )
				: fct_c ( fct_c_ ), op ( op_ )
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
	void addStream ( gpu_stream* stream );
	void  removeStreams();
	fct_call* gpu_dispo ( int gpu_dipo, gpuOp op );

private:
	std::vector<gpu_stream*> v_stream_queu;
};
#endif
