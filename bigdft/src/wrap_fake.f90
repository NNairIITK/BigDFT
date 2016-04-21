!! @file 
!! Fake routines for ambertoold
!! @author 
!!    Copyright (C) 2010-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHOR
subroutine call_nab_gradient(rxyzint,fxyz,epot,icc)
implicit none
real(8), intent(in) :: rxyzint, fxyz, epot
integer, intent(in) :: icc
stop "Ambertools not present"
end subroutine

subroutine nab_init()
  use dictionaries
  call f_err_throw("Ambertools not present")
end subroutine
!!subroutine nab_init(nat,rxyz,fxyz,fnpdb,nfnpdb,l_sat,atomnamesdmy)
!!implicit none
!!integer, intent(in) :: nat, nfnpdb, l_sat
!!real(8), intent(in) :: rxyz(3,nat)
!!real(8), intent(in) :: fxyz(3,nat)
!!character(len=*), intent(in) :: fnpdb
!!character(len=*), intent(in) :: atomnamesdmy(nat)
!!stop "Ambertools not present"
!!end subroutine
