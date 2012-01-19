/**
 * @file manage_global_var.h
 * @author Matthieu Ospici
 * 
 * @brief
 * This file contains all the classes managing global variables.
 *   - global_gpu_attach: inform if process/thread is attached to GPU. This version is for process
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

#ifndef MANAGE_GLOBAL_VAR_H_
#define MANAGE_GLOBAL_VAR_H_

#ifdef HAVE_CONFIG
#include <config.h>
#endif

//store if process/thread is attached to GPU
//this version is for process
class global_gpu_attach
{
public:
	global_gpu_attach()
			: isAttached ( false ) {}
	void setIsAttached ( bool v )
	{
		isAttached = v;
	}
	bool getIsAttached() const
	{
		return isAttached;
	}

private:
	bool isAttached;
};
#endif
