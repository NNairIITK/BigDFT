!!****m* ABINIT/interfaces_14_hidewrite
!! NAME
!! interfaces_14_hidewrite
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/14_hidewrite
!!
!! COPYRIGHT
!! Copyright (C) 2010-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module interfaces_14_hidewrite

 implicit none

interface
 subroutine abi_wrtout(unit,msg,mode_paral)
  implicit none
  integer,intent(in) :: unit
  character(len=*),optional,intent(in) :: mode_paral
  character(len=*),intent(in) :: msg
 end subroutine abi_wrtout
end interface

interface
 subroutine abi_wrtout_myproc(unit,message,mpicomm)
  implicit none
  integer,intent(in),optional :: mpicomm
  integer,intent(in) :: unit
  character(len=*),intent(in) :: message
 end subroutine abi_wrtout_myproc
end interface

interface
 subroutine abi_write_lines(unit,message)
  implicit none
  integer,intent(in) :: unit
  character(len=*),intent(in) :: message
 end subroutine abi_write_lines
end interface

end module interfaces_14_hidewrite
!!***
