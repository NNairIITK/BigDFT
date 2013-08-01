/**
 * @file read_conf_file.h
 * @author Matthieu Ospici
 * 
 * @brief
 * Configuration file management.
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

#ifndef   	READ_CONF_FILE_H_
#define   	READ_CONF_FILE_H_

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include <map>
#include <string>

#include "read_conf_exception.h"
#include "class_utils.h"

typedef std::map<std::string , std::string> mapFile_t;

class readConfFile
{
public:
	readConfFile ( const std::string& filename ) throw ( file_not_found );
	void get ( const std::string& key, std::string& value ) const throw ( read_not_found );
	void get ( const std::string& key, int *value ) const throw ( read_not_found );

private:
	mapFile_t mFile;
};

class readConfFileGPU_CPU : public readConfFile
{
public:
	readConfFileGPU_CPU ( const std::string& filename ) : readConfFile ( filename ) {};
	int getGPU ( int MPI_ID ) const throw ( read_not_found_GPU );
	int getCPU ( int MPI_ID ) const throw ( read_not_found_CPU );
	int getFlag ( int MPI_ID ) const throw(); //0 CUDA, 1 BLAS
};

#endif 	    /* !READ_CONF_FILE_H_ */
/****/
