!!****m* ABINIT/interfaces_18_timing
!! NAME
!! interfaces_18_timing
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/18_timing
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

module interfaces_18_timing

 implicit none

interface
 subroutine timab(nn,option,tottim)
  use defs_basis
  implicit none
  integer,intent(in) :: nn
  integer,intent(in) :: option
  real(dp),intent(out) :: tottim(2)
 end subroutine timab
end interface

interface
 subroutine time_accu(nn,return_ncount,tottim, totflops, totftimes)
  use defs_basis
  implicit none
  integer,intent(in) :: nn
  integer,intent(out) :: return_ncount
  real(dp),intent(out) :: totflops
  real(dp),intent(out) :: totftimes(2)
  real(dp),intent(out) :: tottim(2)
 end subroutine time_accu
end interface

interface
 subroutine timein(cpu,wall)
  use defs_basis
  implicit none
  real(dp),intent(out) :: cpu
  real(dp),intent(out) :: wall
 end subroutine timein
end interface

end module interfaces_18_timing
!!***
