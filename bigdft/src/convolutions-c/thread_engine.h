/*  Copyright (C) 2007 Erik Saule, Brice Videau

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/
#ifndef THREAD_ENGINE_H
#define THREAD_ENGINE_H 1
//#include "tool.h"

struct thread_engine_param {
  char * topology_file;
  char * placement_file;
  int cpu_number;
  int thread_number;
  int **thread_cpus;
  int *thread_cpus_number;
  int *main_cpus;
  int main_cpus_number;
  int memory_allocation;
  int real_time;
  int stack;
  int verbose;
  int use_spin;
};
extern struct thread_engine_param *engine_params;
#ifdef __cplusplus
 extern "C" {
#endif

void get_engine_opt(int argc, char *argv[]);
void init_thread_engine ( );
//void init_thread_engine_spin ( struct thread_engine_param * param );
void run_bench( void * (*main_program)(void * param), void * (*thread_program)(void * param), void ** params ) ;
//void run_bench_spin( void * (*main_program)(void * param), void * (*thread_program)(void * param), void ** params ) ;

#ifdef __cplusplus
 }
#endif
#endif
