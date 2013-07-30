!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_xml_pawpseudo_types
!! NAME
!! m_xml_pseudo_types
!!
!! FUNCTION
!!  Module with pseudopotential structures for PAW XML format
!!
!! COPYRIGHT
!! Copyright (C) 2005-2013 ABINIT group (FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

#include "abi_common_for_bigdft.h"


module m_xml_pawpseudo_types

 use m_profiling

implicit none
!
! Data structures for a prototype pseudopotential
!
integer, parameter, private    :: dp = selected_real_kind(14)
!
!public  :: dump_pseudo
!
!-----------------------------------------------------------
type, public :: radial_grid_t
      logical           :: tread=.false.
      character(len=20)              :: eq
      character(len=4)               :: id
      real(kind=dp)                  :: aa
      real(kind=dp)                  :: bb
      real(kind=dp)                  :: dd
      integer                        :: nn
      integer                        :: istart
      integer                        :: iend
end type radial_grid_t

!-----------------------------------------------------------
! type, public :: ijfunc_t
!       real(kind=dp), dimension(:), pointer    :: data
! end type ijfunc_t
!-----------------------------------------------------------
type, public :: radialfunc_t
      logical           :: tread=.false.
      character(len=6)               ::grid
      character(len=6)               ::state
      real(kind=dp), dimension(:), pointer    :: data
end type radialfunc_t
!-----------------------------------------------------------
type, public :: shape_function_t
      logical           :: tread=.false.
      character(len=20)              :: gtype
      real(kind=dp)                  :: rc
      character(len=6)               :: grid
      integer                        :: lamb
      real(kind=dp), dimension(:,:), pointer    :: data
end type shape_function_t
!-----------------------------------------------------------

type, public :: state_t
      logical           :: tread=.false.
      character(len=6)               :: id
      real(kind=dp)                  :: ff
      real(kind=dp)                  :: rc
      real(kind=dp)                  :: ee
      integer                        :: nn
      integer                        :: ll
end type state_t
!-----------------------------------------------------------
type, public :: valence_states_t
      logical           :: tread=.false.
      integer                        :: nval
      type(state_t),dimension(:), pointer :: state
end type valence_states_t
!-----------------------------------------------------------

type, public :: generator_t
      logical           :: tread=.false.
      character(len=20)               :: gen
      character(len=20)               :: name
end type generator_t
!-----------------------------------------------------------
type, public :: xc_functional_t
      logical           :: tread=.false.
      character(len=12)               :: functionaltype
      character(len=12)               :: name
end type xc_functional_t
!-----------------------------------------------------------
type, public :: atom_t
      logical           :: tread=.false.
        character(len=2)        :: symbol
        real(kind=dp)           :: znucl
        real(kind=dp)           :: zion
        real(kind=dp)           :: zval
end type atom_t
!-----------------------------------------------------------

type, public :: paw_setup_t
      character(len=3)  :: version
      logical           :: tread=.false.
      integer                  :: ngrid
      character(len=4)         :: idgrid
      type(atom_t)             :: atom
      type(xc_functional_t)    :: xc_functional
      type(generator_t)        :: generator
      type(valence_states_t)    :: valence_states
      type(radial_grid_t) ,dimension(:), pointer     :: radial_grid
      type(shape_function_t)   :: shape_function
      type(radialfunc_t)   :: ae_core_density
      type(radialfunc_t)   :: pseudo_core_density
      type(radialfunc_t)   :: pseudo_valence_density
      type(radialfunc_t)   :: zero_potential
      type(radialfunc_t)   :: ae_core_kinetic_energy_density
      type(radialfunc_t)   :: pseudo_core_kinetic_energy_density
      type(radialfunc_t),dimension(:),pointer :: ae_partial_wave
      type(radialfunc_t),dimension(:),pointer :: pseudo_partial_wave
      type(radialfunc_t),dimension(:),pointer :: projector_function
      type(radialfunc_t)   :: kresse_joubert_local_ionic_potential
      type(radialfunc_t)   :: blochl_local_ionic_potential
      type(radialfunc_t)    :: kinetic_energy_differences
end type paw_setup_t



end module m_xml_pawpseudo_types
!!***
