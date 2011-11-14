/**
 * @file fct_call.h
 * @author Matthieu Ospici
 * 
 * @brief
 * fct_call is a class designed to be extended in order to create new code for GPU
 * the virtual operator() must be extended
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
#ifndef   	FCT_CALL_H_
# define   	FCT_CALL_H_

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include "sg_common_def.h"

class fct_call
{
public:
#ifdef TRACE_S_GPU
	fct_call ( tracer_callbacks_t *_tcal ) : tcal ( _tcal ) {}
#endif
	virtual void operator() ( int )  = 0;
	virtual ~fct_call() {};

protected:
#ifdef TRACE_S_GPU
	tracer_callbacks_t *tcal;
#endif
};
#endif 	    /* !FCT_CALL_H_ */
