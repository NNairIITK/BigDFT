/**
 * @file class_utils.h
 * @author Matthieu Ospici
 * 
 * @brief
 * Some classes usefull for S_GPU library.
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

#ifndef CLASS_UTILS_H_
#define CLASS_UTILS_H_

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include <sstream>
#include <unistd.h>

class deleter
{
public:
	template <class T> void operator () ( T*& p ) const
	{
		delete p;
		p = NULL;
	}
};
inline void getHostName ( std::string& h )
{
	const int HOST_NAME_SIZE = 300;
	char hostname[HOST_NAME_SIZE];
	gethostname ( hostname, HOST_NAME_SIZE );
	h = hostname;
}

template<typename T>
T strTo ( const std::string& str )
{
	T dest;
	std::istringstream iss ( str );
	iss >> dest;
	return dest;
}
#endif
