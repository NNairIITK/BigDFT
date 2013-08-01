/**
 * @file message.h
 * @author Matthieu Ospici
 * 
 * @brief
 * Definition of the message C structure.
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

#ifndef MESSAGE_H_
#define MESSAGE_H_

#ifdef HAVE_CONFIG
#include <config.h>
#endif

enum msgType_t {mCALC, mTRANSF, mSTOP, mSEM, mTERMINATE };

struct message
{
	msgType_t msgType;
	int gpu_avail; //this token means gpu_avail is available
	// unsigned int memAvaillable;
	int node; //for mSTOP, node has finish, or mSEM, semID...
};
#endif
