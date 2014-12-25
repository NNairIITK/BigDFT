!!****m* ABINIT/interfaces_16_hideleave
!! NAME
!! interfaces_16_hideleave
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/16_hideleave
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

module interfaces_16_hideleave

 implicit none

interface
 subroutine abi_leave_new(mode_paral,exit_status,print_config)
  implicit none
  integer,intent(in),optional :: exit_status
  character(len=4),intent(in) :: mode_paral
  logical,intent(in),optional :: print_config
 end subroutine abi_leave_new
end interface

end module interfaces_16_hideleave
!!***
