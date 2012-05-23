/*
!> @file
!!  Routines to read YAML position files.
!! @author
!!    Copyright (C) 2007-2012 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
*/

#ifndef ATOMS_YAML_H
#define ATOMS_YAML_H

typedef struct _PosinpAtoms
{
  /* The cell. */
  unsigned int BC, Units;
  double acell[3], angdeg[3];

  /* The positions. */
  unsigned int nat, ntypes, units;
  double *rxyz;
  char **atomnames;
  unsigned int *iatype, *ifrztyp;
  int *igspin, *igchg;

  /* The forces. */
  /* To be done. */

  /* Additional data. */
  /* To be done. */
} PosinpAtoms;

typedef struct _PosinpList PosinpList;
struct _PosinpList
{
  PosinpList *next;
  PosinpAtoms *data;
};

PosinpList* posinp_yaml_parse(const char *filename);
void posinp_yaml_free_list(PosinpList *lst);

#endif
