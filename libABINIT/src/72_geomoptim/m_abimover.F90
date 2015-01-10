!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_abimover
!! NAME
!! m_abimover
!!
!! FUNCTION
!! This module contains type definitions for the molecular dynamics:
!! * mttk_type : data For Martyna et al. (TTK) MD integration scheme
!!
!! COPYRIGHT
!! Copyright (C) 2001-2014 ABINIT group (JJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_abimover

 use defs_basis

 implicit none
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/mttk_type
!! NAME
!! mttk_type
!!
!! FUNCTION
!! For Martyna et al. (TTK) reversible MD integration scheme and related data
!!
!! SOURCE

 type mttk_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

!Real (double precision) scalars

   real(dp) :: glogv
    !Logarithm of the volume

   real(dp) :: vlogv
    !Derivative of logv

!Real (double precision) arrays

  real(dp) :: gboxg(3,3)
   !Imbalance in pressure (see paper)

  real(dp) :: vboxg(3,3)
   !Velocity of log rprimd (see paper)

  real(dp), pointer :: glogs(:)
   ! glogs(nnos)
   ! Imbalance of kinetic energy

  real(dp), pointer :: vlogs(:)
   ! vlogs(nnos)
   ! Velocities of thermostat variables

  real(dp), pointer :: xlogs(:)
   ! xlogs(nnos)
   ! Positions of thermostat variables

 end type mttk_type
!!***

!----------------------------------------------------------------------

end module m_abimover
