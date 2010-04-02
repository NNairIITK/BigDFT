!!****m* ABINIT/interfaces_16_hideleave
!! NAME
!! interfaces_16_hideleave
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/16_hideleave
!!
!! COPYRIGHT
!! Copyright (C) 2010 ABINIT group
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

module interfaces_16_hideleave

 implicit none

interface
 subroutine leave_myproc
  implicit none
 end subroutine leave_myproc
end interface

interface
 subroutine leave_new(mode_paral)
  implicit none
  character(len=4),intent(in) :: mode_paral
 end subroutine leave_new
end interface

interface
 subroutine print_ierr(ierr,nrcalled,nrcalling)
  implicit none
  integer, intent(in) :: ierr
  character(len=*),intent(in) :: nrcalled
  character(len=*),intent(in) :: nrcalling
 end subroutine print_ierr
end interface

end module interfaces_16_hideleave
!!***
