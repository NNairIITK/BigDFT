/**
 * @file manage_cpu_affinity.h
 * @author Matthieu Ospici
 * 
 * @brief
 * Manage the affinity between the CPUs and the GPUs. At the initialization part of S_GPU,
 * using config information (in config file GPU.config), it is defined that porcesses located
 * on some CPUs will use specific one GPU
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

#ifndef MANAGE_CPU_AFFINITY_H_
#define MANAGE_CPU_AFFINITY_H_

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include <vector>
#include <list>
#include <string>

//define the ID of CPUs connected to one GPU
class gpu_cpus_connexion
{
public:
	//  gpu_cpus_connexion():num_connected(0){}
	gpu_cpus_connexion ( const std::string& st );
	void add_cpu_from_str ( const std::string& ); // string syntax : "cpu1,cpu2,....,cpun
	void add_cpu ( int cpuNum );
	int set_affinity ( int cpuID ) const; //set current tread/process affinity with an avaible cpu number
	int get_num_cpus() const
	{
		return cpus_connected.size();
	}

private:
	std::vector<int> cpus_connected;
};

//======================================

class manage_cpu_affinity
{
public:
	manage_cpu_affinity ( int _iproc ) : iproc ( _iproc ) {}
	void add_connexion ( const gpu_cpus_connexion& );
	//set the calling process affinity to a processor that we have to connect to gpu_to_atach GPU
	void set_affinity ( int cpu_aff ) ;
	void print_affinity_matrix() const;

private:
	std::list<gpu_cpus_connexion> connexions;
	std::vector<int> affinity_matrix; //store afinity
	int iproc;
};
#endif
