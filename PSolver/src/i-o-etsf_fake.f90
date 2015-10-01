!> @file
!!    File for handling the read-write of a given simulation in ETSF format, fake version
!!
!! @author
!!    G. Fisicaro, L. Genovese (September 2015)
!!    Copyright (C) 2002-2015 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!> Read a field in the ISF basis in the ETSF format
subroutine read_etsf(filename,geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho,&
     nat,rxyz)
  use PSbase
  implicit none
  character(len=*), intent(in) :: filename
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: nspin
  integer, intent(out) ::  n1i,n2i,n3i
  real(gp), intent(out) :: hxh,hyh,hzh
  real(dp), dimension(:,:), pointer :: rho
  real(gp), dimension(:,:), pointer,  optional :: rxyz
  integer, intent(out), optional ::  nat

  write(0, "(A)") "Illegal call to read_etsf(), not compiled with ETSF_IO support."
  stop

  !To avoid warnings from the compiler
  write(*,*) filename,geocode,nspin,n1i,n2i,n3i
  n1i=0
  n2i=0
  n3i=0
  hxh=0.0_gp
  hyh=0.0_gp
  hzh=0.0_gp
  rho(1,1)=0.0_gp
  if (present(nat)) then
      rxyz(1,1)=0.0_gp
     nat=0
  end if
END SUBROUTINE read_etsf

