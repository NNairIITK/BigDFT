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
#include <pthread.h>
#include <iostream>
#ifndef u32
typedef unsigned int u32;
#endif
//#include <asm/msr.h>
#include <sched.h>
#include <sys/types.h>
#include <unistd.h>
#include <semaphore.h>
#include <stdlib.h>
//#include <errno.h>
#include <libconfig.h++>
#include "thread_engine.h"

#define ASMPAUSE asm("" : : : "memory")


struct thread_data {
  void * ( * volatile thread_code)(void * param);
  void * volatile param;
  pthread_mutex_t mutex;
  pthread_cond_t cond;
  pthread_spinlock_t spin;
};

struct thread_param {
  int thread_number;
};

struct thread_engine_param *engine_params;
struct thread_param * thread_params;
struct thread_data * volatile * thread_datas;
pthread_t * threads;
sem_t semaphore;
pthread_spinlock_t global_spin;
volatile unsigned int thread_counter;


void get_topology( ) {
	if( engine_params->verbose ) std::cout << "reading topology" << std::endl;
        libconfig::Config topology;
        topology.readFile(engine_params->topology_file);
        engine_params->cpu_number = topology.lookup("topology.cpu_number");
}

void get_placement( ) {
	try {
	        libconfig::Config thread_placement;
		if( engine_params->verbose ) std::cout << "reading placement" << std::endl;
	        thread_placement.readFile(engine_params->placement_file);
        	libconfig::Setting& thread_cpus = thread_placement.lookup("placement.thread_cpus");
	        libconfig::Setting& main_cpus = thread_placement.lookup("placement.main_cpus");
	        engine_params->thread_number = thread_cpus.getLength();
	        engine_params->main_cpus = new int[main_cpus.getLength()];
	        engine_params->main_cpus_number = main_cpus.getLength();
	        for( int i=0; i<main_cpus.getLength(); i++ ) {
	                engine_params->main_cpus[i] = main_cpus[i];
	        }
	        engine_params->thread_cpus = new int*[thread_cpus.getLength()];
	        engine_params->thread_cpus_number = new int[thread_cpus.getLength()];
	        for( int i=0; i<thread_cpus.getLength(); i++ ) {
	                engine_params->thread_cpus[i] = new int[thread_cpus[i].getLength()];
	                engine_params->thread_cpus_number[i] = thread_cpus[i].getLength();
	                for( int j=0; j<thread_cpus[i].getLength(); j++) {
	                        engine_params->thread_cpus[i][j] = thread_cpus[i][j];
	                }
	        }
        	engine_params->memory_allocation = thread_placement.lookup("placement.memory_allocation");
                engine_params->use_spin = (bool)thread_placement.lookup("use_spin");
	} catch (libconfig::ParseException e) {
		std::cerr << "libconfig::ParseException" << std::endl;
		std::cerr <<e.getLine () <<" : "<<e.getError() << std::endl;
		exit(-1);
	}
}

void get_engine_opt(int argc, char *argv[]) {
	int ch;
	char *endptr;
        engine_params = new struct thread_engine_param;
        //errno = 0;

        engine_params->thread_number = 0;
	engine_params->cpu_number = 1;
	engine_params->real_time = false;
	engine_params->topology_file = NULL;
	engine_params->placement_file = NULL;
	engine_params->verbose = 0;

	while ((ch = getopt(argc, argv, "t:p:rsv:")) != -1) {
		switch (ch) {
                        case 't':
                                engine_params->topology_file = optarg;
                                break;
                        case 'p':
                                engine_params->placement_file = optarg;
                                break;
			case 'r':
				//real_time
				engine_params->real_time = true;
				break;
			case 'v':
				engine_params->verbose = strtol(optarg, &endptr, 10);
//				if( errno != 0)  {
//					perror("Invalid verbose level");
//				};
				break;
			default:
				printf ("Unknown parameter : %c\n",ch);
				break;
		}
	}
	if( engine_params->topology_file )
		get_topology();
	else {
		std::cerr << "no topology file given" << std::endl;
		exit(-1);
	}
	if( engine_params->placement_file )
		get_placement();
	else {
		std::cerr << "no placement file given" << std::endl;
		exit(-1);
	}
	argc -= optind;
	argv += optind;	
}

void * thread_master_function_spin( void *arg ) {
  struct thread_param *param;
  param = (struct thread_param *)arg;
  struct thread_data data;
  data.thread_code = NULL;
  data.param = NULL;
  pthread_spin_init(&(data.spin),PTHREAD_PROCESS_PRIVATE);
  thread_datas[param->thread_number] = &data;
  sem_post(&semaphore);
  while(1) {
    while (1) {
      while (data.thread_code == NULL){ASMPAUSE;}
      pthread_spin_lock(&(data.spin));//
      if (data.thread_code != NULL) break;
      pthread_spin_unlock(&(data.spin));
    }
    pthread_spin_unlock(&(data.spin));
    data.thread_code(data.param);
    data.param = NULL;
    data.thread_code = NULL;
    pthread_spin_lock(&global_spin);
    thread_counter--;
    pthread_spin_unlock(&global_spin);
  }
  return NULL;
}

void * thread_master_function( void *arg ) {
  struct thread_param *param;
  param = (struct thread_param *)arg;
  struct thread_data data;
  data.thread_code = NULL;
  data.param = NULL;
  pthread_mutex_init(&(data.mutex),NULL);
  pthread_cond_init(&(data.cond),NULL);
  thread_datas[param->thread_number] = &data;
  while ( 1 ) {
    sem_post(&semaphore);
    pthread_mutex_lock(&(data.mutex));
    while(data.thread_code == NULL)
      pthread_cond_wait(&(data.cond),&(data.mutex));
    pthread_mutex_unlock(&(data.mutex));

    data.thread_code(data.param);
    data.thread_code = NULL;
    data.param = NULL;

  }
  return NULL;
}

void set_cpu_sets(int *cpus, int cpus_number ) {
  cpu_set_t cpuset;
  cpu_set_t fullcpuset;
  CPU_ZERO(&fullcpuset);
  int cpu_number = engine_params->cpu_number;
  for( int i=0; i<cpu_number; i++) {
    CPU_SET(i,&fullcpuset);
  }
  if(cpu_number > 0) {
    if(engine_params->verbose) std::cout << "cpus :";
    CPU_ZERO(&cpuset);
    for( int j=0; j<cpus_number; j++ ){
      if(engine_params->verbose) std::cout << " "<< cpus[j];
      if(cpus[j]>=cpu_number || cpus[j]<0) { std::cerr << "Invalid CPU number : " << cpus[j] << std::endl; exit(1);}
      CPU_SET(cpus[j], &cpuset);
    }
    if(engine_params->verbose) std::cout << std::endl;
    if(pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset)<0)
      perror("pthread_setaffinity_np");
    sched_yield();
  } else {
    if(engine_params->verbose) std::cout << "all cpus" << std::endl;
    if(pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &fullcpuset)<0)
      perror("pthread_setaffinity_np");
    sched_yield();
  }
}

void init_thread_engine ( ) {
  if(engine_params->real_time) {       
    struct sched_param sched_p;
    sched_p.__sched_priority = sched_get_priority_max (SCHED_RR);
    if (sched_setscheduler (0, SCHED_RR, &sched_p) != 0)
      perror("sched_setscheduler");
  }
  thread_datas = new thread_data*[engine_params->thread_number];
  thread_params = new thread_param[engine_params->thread_number];
  threads = new pthread_t[engine_params->thread_number];
  if(engine_params->use_spin) {
    thread_counter=0;
    pthread_spin_init(&(global_spin),PTHREAD_PROCESS_PRIVATE);
  }
  sem_init(&semaphore,0,0);
  for( int i=0; i<engine_params->thread_number; i++ ) {
    if( engine_params->verbose ) std::cout << "binding thread : " << i+1 << std::endl; 
    set_cpu_sets(engine_params->thread_cpus[i], engine_params->thread_cpus_number[i]);
    thread_params[i].thread_number = i;
    if(engine_params->use_spin) {
      pthread_create( &(threads[i]), NULL, thread_master_function_spin, (void *) &(thread_params[i]) );
    } else {
      pthread_create( &(threads[i]), NULL, thread_master_function, (void *) &(thread_params[i]) );
    }
    sem_wait(&semaphore);
  }
  if( engine_params->verbose ) std::cout << "binding  main thread (0)" << std::endl; 
  set_cpu_sets(engine_params->main_cpus, engine_params->main_cpus_number);
}

void _run_bench_spin( void * (*main_program)(void * param), void * (*thread_program)(void * param), void ** params ) {
  for( int i=0; i<engine_params->thread_number; i++ ) {
    struct thread_data *d = thread_datas[i];
    pthread_spin_lock(&global_spin);
    thread_counter++;
    pthread_spin_unlock(&global_spin);
    pthread_spin_lock(&(d->spin));
    d->param = params[i+1];
    d->thread_code = thread_program;
    pthread_spin_unlock(&(d->spin));
  }
  main_program(params[0]);
  while (1) {
    while (thread_counter != 0){ASMPAUSE;}
    pthread_spin_lock(&global_spin);
    if (thread_counter == 0) break;
    pthread_spin_unlock(&global_spin);
  }
  pthread_spin_unlock(&global_spin);
}

void _run_bench( void * (*main_program)(void * param), void * (*thread_program)(void * param), void ** params ) {
  for( int i=0; i<engine_params->thread_number; i++ ) {
    struct thread_data *d = thread_datas[i];
    pthread_mutex_lock(&(d->mutex));
    d->param = params[i+1];
    d->thread_code = thread_program;
    pthread_cond_signal(&(d->cond));
    pthread_mutex_unlock(&(d->mutex));
  }
  main_program(params[0]);
  for( int i=0; i<engine_params->thread_number; i++ ) {
    sem_wait(&semaphore);
  }
}

void run_bench( void * (*main_program)(void * param), void * (*thread_program)(void * param), void ** params ) {
  if(engine_params->use_spin) {
    _run_bench_spin(main_program, thread_program, params);
  } else {
    _run_bench(main_program, thread_program, params);
  }
}
