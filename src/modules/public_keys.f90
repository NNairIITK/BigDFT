!> @file
!! Contains public parameters that have to be used here and there in the code
!! @author
!!    Copyright (C) 2007-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module public_keys
  implicit none

  public ! guess why?

  !> Parameters to avoid typos in dictionary keys
  character(len=*), parameter :: ASTRUCT_UNITS = 'units' 
  character(len=*), parameter :: ASTRUCT_CELL = 'cell' 
  character(len=*), parameter :: ASTRUCT_POSITIONS = 'positions' 
  character(len=*), parameter :: ASTRUCT_PROPERTIES = 'properties' 
  character(len=*), parameter :: GOUT_ENERGY = 'energy (Ha)' 
  character(len=*), parameter :: GOUT_FORCES = 'forces (Ha/Bohr)' 
  character(len=*), parameter :: FORMAT_KEY = 'format' 
  character(len=*), parameter :: OCCUPATION = 'occupation' 
  character(len=*), parameter :: FORMAT_YAML = 'yaml' 
  character(len=*), parameter :: RADII_KEY = 'Radii of active regions (AU)' 
  character(len=*), parameter :: LPSP_KEY = 'Local Pseudo Potential (HGH convention)' 
  character(len=*), parameter :: NLPSP_KEY = 'NonLocal PSP Parameters'
  character(len=*), parameter :: PSPXC_KEY = 'Pseudopotential XC'
  character(len=*), parameter :: PSP_TYPE = 'Pseudopotential type'
  character(len=*), parameter :: COARSE = 'Coarse'
  character(len=*), parameter :: COARSE_PSP = 'Coarse PSP'
  character(len=*), parameter :: FINE = 'Fine'
  character(len=*), parameter :: SOURCE_KEY = 'Source'
  character(len=*), parameter :: ATOMIC_NUMBER = 'Atomic number'
  character(len=*), parameter :: ELECTRON_NUMBER = 'No. of Electrons'


end module public_keys
