//! @file
//!  Wrappers for OpenCL (??)
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include "OpenCL_wrappers.h"

size_t event_number;
size_t event_allocated;
event * event_list;


void FC_FUNC_(init_event_list,INIT_EVENT_LIST)(){
  event_list = NULL;
  event_allocated = 0;
  event_number = 0;
}

int addToEventList (event ev)
{
        if(event_number == event_allocated) // Are more refs required?
        {
                // Feel free to change the initial number of refs
                // and the rate at which refs are allocated.
                if (event_allocated == 0)
                        event_allocated = 100; // Start off with 3 refs
                else
                        event_allocated *= 2; // Double the number
                                                    // of refs allocated

                // Make the reallocation transactional
                // by using a temporary variable first
                void *_tmp = realloc(event_list, (event_allocated * sizeof(event)));

                // If the reallocation didn't go so well,
                // inform the user and bail out
                if (!_tmp)
                {
                        fprintf(stderr, "ERROR: Couldn't realloc memory!\n");
                        return(-1);
                }

                // Things are looking good so far
                event_list = (event*)_tmp;
        }

        event_list[event_number] = ev;
        event_number++;

        return event_number;
}

void FC_FUNC_(print_event_list,PRINT_EVENT_LIST)() {
#if PROFILING
  FILE * f;
  size_t i;
  cl_ulong  queued,submit,start,end;
  cl_int ciErrNum;
  event e;
  f = fopen("cl_profiling.yaml","w");
  fprintf(f,"---\n");
  for(i=0;i<event_number;i++){
    e = event_list[i];
    ciErrNum = clGetEventProfilingInfo(e.e, CL_PROFILING_COMMAND_QUEUED, sizeof(queued), &queued, NULL );
    oclErrorCheck(ciErrNum,"Failed to get profiling info!");    
    ciErrNum = clGetEventProfilingInfo(e.e, CL_PROFILING_COMMAND_SUBMIT, sizeof(submit), &submit, NULL );
    oclErrorCheck(ciErrNum,"Failed to get profiling info!");    
    ciErrNum = clGetEventProfilingInfo(e.e, CL_PROFILING_COMMAND_START, sizeof(start), &start, NULL );
    oclErrorCheck(ciErrNum,"Failed to get profiling info!");    
    ciErrNum = clGetEventProfilingInfo(e.e, CL_PROFILING_COMMAND_END, sizeof(end), &end, NULL );
    oclErrorCheck(ciErrNum,"Failed to get profiling info!");    
    fprintf(f,"-\n");
    fprintf(f,"  event_queued: %lu\n", queued);
    fprintf(f,"  event_submit: %lu\n", submit);
    fprintf(f,"  event_start: %lu\n", start);
    fprintf(f,"  event_end: %lu\n", end);
    fprintf(f,"  comment: %s\n", e.comment);
  }
  fprintf(f,"...\n");
  fclose(f);
#endif
}
