/*
!> @file
!!  Routines to access easily to the filesystem (C part).
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
#include <stdlib.h>
#include <string.h>

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
          out[lgCpy] = '/';
        }
      else
        *status = 1;
      return;
    }

  /* Try to create it. */
  if (mkdir(path, 0755) != 0)
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
