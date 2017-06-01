!!****m* ABINIT/abi_interfaces_lowlevel
!! NAME
!! abi_interfaces_lowlevel
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/abi_lowlevel
!!
!! COPYRIGHT
!! Copyright (C) 2010-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module abi_interfaces_lowlevel

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

interface
 subroutine abi_leave_new(mode_paral,exit_status,print_config)
  implicit none
  integer,intent(in),optional :: exit_status
  character(len=4),intent(in) :: mode_paral
  logical,intent(in),optional :: print_config
 end subroutine abi_leave_new
end interface

end module abi_interfaces_lowlevel
!!***
