/**
 * @file exceptions.h
 * @author Matthieu Ospici
 * 
 * @brief
 * Exception handler. They are principally used to check returns of CUDA functions.
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

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include <string>
#include <exception>
#include <sstream>

#include "class_utils.h"

template<class E>
inline void check ( bool i, const std::string& st )
{
	if ( i )
	{
		//   std::string line(__LINE_);
		//    std::string file(__FILE__);
		std::ostringstream oss;
		std::string host;
		getHostName ( host );
		oss << "ERROR : " << st << std::endl << "host : " << host  << std::endl;
		throw E ( oss.str() );
	}
}

template<class E>
inline void check ( bool i, const std::string& st, const char* file, int line )
{
	if ( i )
	{
		//   std::string line(__LINE_);
		//    std::string file(__FILE__);
		std::ostringstream oss;
		std::string host;
		getHostName ( host );
		oss << "ERROR : " << st << std::endl << "host : " << host << ", file : " << file << ", line : " << line << std::endl;
		throw E ( oss.str() );
	}
}

class inter_node_communication_error: public std::exception
{
public:
	inter_node_communication_error ( std::string _msg ) throw()
			: msg ( _msg )
	{}
	virtual const char* what() const throw()
	{
		return msg.c_str();
	}
	virtual ~inter_node_communication_error() throw()
	{}

private:
	std::string msg;       //Error description
};

class synchronization_error: public std::exception
{
public:
	synchronization_error ( std::string _msg ) throw()
			: msg ( _msg )
	{}
	virtual const char* what() const throw()
	{
		return msg.c_str();
	}
	virtual ~synchronization_error() throw()
	{}

private:
	std::string msg;       //Error description
};

class check_calc_error: public std::exception
{
public:
	check_calc_error ( std::string _msg ) throw()
			: msg ( _msg )
	{}
	virtual const char* what() const throw()
	{
		return msg.c_str();
	}
	virtual ~check_calc_error() throw()
	{}

private:
	std::string msg;       //Error description
};

class repartition_error: public std::exception
{
public:
	repartition_error ( std::string _msg ) throw()
			: msg ( _msg )
	{}
	virtual const char* what() const throw()
	{
		return msg.c_str();
	}
	virtual ~repartition_error() throw()
	{}

private:
	std::string msg;       //Error description
};

class other_error: public std::exception
{
public:
	other_error ( std::string _msg ) throw()
			: msg ( _msg )
	{}
	virtual const char* what() const throw()
	{
		return msg.c_str();
	}
	virtual ~other_error() throw()
	{}

private:
	std::string msg;       //Error description
};

class cuda_error: public std::exception
{
public:
	cuda_error ( std::string _msg ) throw()
			: msg ( _msg )
	{}
	virtual const char* what() const throw()
	{
		return msg.c_str();
	}
	virtual ~cuda_error() throw()
	{}

private:
	std::string msg;       //Error description
};
#endif
