/**
 * @file manage_cpu_affinity.cpp
 * @author Matthieu Ospici
 * 
 * @brief
 * Manage the affinity between the CPUs and the GPUs. At the initialization part of S_GPU,
 * using config information (in config file GPU.config), it is defined that processes located
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

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include <iostream>
#include <sched.h>
#include <sstream>

#include "manage_cpu_affinity.h"
#include "class_utils.h"

void gpu_cpus_connexion::add_cpu ( int cpuNum )
{
	cpus_connected.push_back ( cpuNum );
}

gpu_cpus_connexion::gpu_cpus_connexion ( const std::string& st )
{
	add_cpu_from_str ( st );
}

void gpu_cpus_connexion::add_cpu_from_str ( const std::string& st )
{
	std::string token;
	std::istringstream iss ( st );
	while ( getline ( iss, token, ',' ) )
	{
		//      std::cout << "token : " << strTo<int>(token) << std::endl;
		cpus_connected.push_back ( strTo<int> ( token ) );
	}
}

int gpu_cpus_connexion::set_affinity ( int cpuID ) const
{
	int aff_CPU_ID = cpus_connected.at ( cpuID );
	//  std::cout << cpuID << "avec " << aff_CPU_ID << std::endl;
	cpu_set_t cpus;
	CPU_ZERO ( &cpus );
	CPU_SET ( aff_CPU_ID, &cpus );
	//  std::cout << "on set proc - " << CPUaff << std::endl;
	if ( sched_setaffinity ( 0, sizeof ( cpu_set_t ), &cpus ) < 0 )
	{
		std::cerr << "Error setaffinity" << std::endl;
	}

	if ( !CPU_ISSET ( aff_CPU_ID, &cpus ) )
	{
		std::cout << "*** CPU - " << aff_CPU_ID << " - is not set ***" << std::endl;
	}
	return aff_CPU_ID;
}

//=============================================================================

void manage_cpu_affinity::add_connexion ( const gpu_cpus_connexion& gcc )
{
	connexions.push_back ( gcc );
}

void manage_cpu_affinity::set_affinity ( int cpuID )
{
	std::list<gpu_cpus_connexion>::const_iterator cit = connexions.begin() ;
	int cur_cpu_id = 0;
	bool found = false;
	//  std::cout << "CPUID " << cpuID << std::endl;
	//  for (cit=connexions.begin(); cit != connexions.end(); ++cit)
	while ( !found && cit != connexions.end() )
	{
		//  cur_cpu_id  += cit->get_num_cpus();
		if ( cpuID < ( cur_cpu_id + cit->get_num_cpus() ) ) //this id is on this block
		{
			int test = cit->set_affinity ( cpuID - cur_cpu_id );
			affinity_matrix.push_back ( test );

			found = true;
		}
		cur_cpu_id  += cit->get_num_cpus();
		++cit;
	}
	if ( !found )
		affinity_matrix.push_back ( -1 );
	//connexions.at(gpu_to_atach).set_current_affinity();
}

void manage_cpu_affinity::print_affinity_matrix() const
{
	std::vector<int>::const_iterator cit = affinity_matrix.begin() ;
	int i;
	for ( cit = affinity_matrix.begin(), i = 0 ; cit != affinity_matrix.end() ; ++cit, ++i )
	{
		std::cout << "Logic CPU : " << i << " with real CPU " << *cit << std::endl;
	}
}
