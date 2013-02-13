/*
!> @file
!!  Routines to access easily to the filesystem and to other C goodies.
!! @author
!!    Copyright (C) 2007-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
*/

#include <config.h>

#define _GNU_SOURCE

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_TIME_H
#include <time.h>
#endif

#ifndef HAVE_STRNDUP
char* strndup(const char *src, size_t len);
#endif

#ifndef HAVE_CLOCK_GETTIME
#define CLOCK_REALTIME 0
static int clock_gettime(int clk_id, struct timespec *tp)
{
  struct timeval now;
  int rv;

  tp->tv_sec = 0;
  tp->tv_nsec = 0;

  rv = gettimeofday(&now, NULL);
  if (rv != 0)
    return rv;
        
  tp->tv_sec = now.tv_sec;
  tp->tv_nsec = now.tv_usec * 1000;

  return 0;
}
#endif

void FC_FUNC(nanosec,NANOSEC)(unsigned long long int * t){
  struct timespec time;
  clock_gettime(CLOCK_REALTIME, &time);
  *t = time.tv_sec;
  *t *= 1000000000;
  *t += time.tv_nsec;
}

void FC_FUNC(getaddress, GETADDRESS)(void *ptr,char *address, int *lgaddress,
			     int* status)
{
  char buff[50]; //test buffer to check the length
  int lgt,lgCpy;

  memset(address,' ', sizeof(char) * (*lgaddress));

  lgt=sprintf(buff,"%p",(void*)ptr);
  if (lgt <= *lgaddress)
    {
      *status = 0;
      lgt=sprintf(address,"%p",(void*)ptr);
    }
  else
    *status = 1;
  return;

  //printf("\n test address = %p %d; \n", (void*)ptr,lgt);
  //return;
}

void FC_FUNC(getdir, GETDIR)(const char *dir, int *lgDir,
                             char *out, int *lgOut,
                             int *status)
{
  char *path;
  struct stat sb;
  int lgCpy;

  *status = 1;
  memset(out, ' ', sizeof(char) * (*lgOut));

  path = strndup(dir, (size_t)*lgDir);

  if (stat(path, &sb) == 0)
    {
      free(path);
      /* Dir exists already. */
      if (S_ISDIR(sb.st_mode))
        {
          *status = 0;
          lgCpy = ((*lgDir > *lgOut - 1)?*lgOut - 1:*lgDir);
          memcpy(out, dir, sizeof(char) * lgCpy);
          /* Add a '/' if not already present */
          if (out[lgCpy-1] != '/') { out[lgCpy] = '/'; };
        }
      else
        *status = 1;
      return;
    }

  /* Try to create it. */
#ifdef _WIN32
  if (mkdir(path) != 0)
#else
  if (mkdir(path, 0755) != 0)
#endif
    {
      free(path);
      *status = 2;
      return;
    }

  free(path);
  lgCpy = ((*lgDir > *lgOut - 1)?*lgOut - 1:*lgDir);
  memcpy(out, dir, sizeof(char) * lgCpy);
  out[lgCpy] = '/';
  *status = 0;
  return;
}

void FC_FUNC(delete, DELETE)(const char *f, int *lgF, int *status)
{
  char *path;

  path = strndup(f, (size_t)*lgF);
  *status = unlink(path);
  free(path);
}

void FC_FUNC(deldir, DELDIR)(const char *f, int *lgF, int *status)
{
  char *path;

  path = strndup(f, (size_t)*lgF);
  *status = rmdir(path);
  free(path);
}


void FC_FUNC(movefile, MOVEFILE)(const char *oldfile, int *lgoldfile, const char *newfile, int *lgnewfile, int *status)
{
  char *oldpath;
  char *newpath;

  oldpath = strndup(oldfile, (size_t)*lgoldfile);
  newpath = strndup(newfile, (size_t)*lgnewfile);
  *status = rename(oldpath,newpath);
  free(oldpath);
  free(newpath);
}
