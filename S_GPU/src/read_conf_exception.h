/**
 * @file read_conf_exception.h
 * @author Matthieu Ospici
 * 
 * @brief
 * Handle exceptions for conf file reading
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

#ifndef READ_CONF_EXCEPTION_H_
#define READ_CONF_EXCEPTION_H_

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include <string>
#include <exception>

class file_not_found: public std::exception
{
public:
	file_not_found ( std::string _msg ) throw()
			: msg ( _msg )
	{}
	virtual const char* what() const throw()
	{
		return msg.c_str();
	}
	virtual ~file_not_found() throw()
	{}

private:
	std::string msg;       //Error description
};

class read_not_found: public std::exception
{
public:
	read_not_found ( std::string _msg ) throw()
			: msg ( _msg )
	{}
	virtual const char* what() const throw()
	{
		return msg.c_str();
	}
	virtual ~read_not_found() throw()
	{}

private:
	std::string msg;       //Error description
};

class read_not_found_GPU : public read_not_found
{
public:
	read_not_found_GPU ( std::string _msg ) throw()
			: read_not_found ( _msg ) {}
	virtual ~read_not_found_GPU() throw()
	{}
};

class read_not_found_CPU : public read_not_found
{
public:
	read_not_found_CPU ( std::string _msg ) throw()
			: read_not_found ( _msg ) {}
	virtual ~read_not_found_CPU() throw()
	{}
};
#endif
