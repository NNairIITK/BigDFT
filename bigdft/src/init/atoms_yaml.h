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

typedef enum
  {
    POSINP_BC_FREE,
    POSINP_BC_WIRE_Z,
    POSINP_BC_SURFACE_XZ,
    POSINP_BC_PERIODIC,
    POSINP_N_BC
  } PosinpBC;

typedef enum
  {
    POSINP_CELL_UNITS_BOHR,
    POSINP_CELL_UNITS_ANGSTROEM,
    POSINP_CELL_N_UNITS
  } PosinpCellUnits;
typedef enum
  {
    POSINP_COORD_UNITS_BOHR,
    POSINP_COORD_UNITS_ANGSTROEM,
    POSINP_COORD_UNITS_REDUCED,
    POSINP_COORD_N_UNITS
  } PosinpCoordUnits;
typedef enum
  {
    POSINP_FORCE_UNITS_HARTREE_PER_BOHR,
    POSINP_FORCE_UNITS_EV_PER_ANGSTROEM,
    POSINP_FORCE_N_UNITS
  } PosinpForceUnits;
typedef enum
  {
    POSINP_ENERG_UNITS_HARTREE,
    POSINP_ENERG_UNITS_EV,
    POSINP_ENERG_UNITS_RYDBERG,
    POSINP_ENERG_N_UNITS
  } PosinpEnergUnits;
typedef enum
  {
    POSINP_FROZEN_FREE,
    POSINP_FROZEN_FULL,
    POSINP_FROZEN_Y,
    POSINP_FROZEN_XZ,
    POSINP_N_FROZEN
  } PosinpFrozenType;

typedef struct _PosinpAtoms
{
  /* The cell. */
  PosinpBC BC;
  PosinpCellUnits Units;
  double acell[3], angdeg[3];

  /* The positions. */
  unsigned int nat, ntypes;
  PosinpCoordUnits units;
  double *rxyz;
  char **atomnames;
  unsigned int *iatype, *ifrztyp;
  int *igspin, *igchg;

  /* The forces. */
  PosinpForceUnits funits;
  double fnrm, maxval;
  double *fxyz;

  /* Additional data. */
  char *comment;
  PosinpEnergUnits eunits;
  double energy;
  unsigned int converged;
  double gnrm_wfn;
} PosinpAtoms;

typedef struct _PosinpList PosinpList;
struct _PosinpList
{
  PosinpList *next;
  PosinpAtoms *data;
};

PosinpList* posinp_yaml_parse(const char *filename, char **message);
void posinp_yaml_free_list(PosinpList *lst);

#endif
