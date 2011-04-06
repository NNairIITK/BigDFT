!> @file
!!  Fake routines for ETSF-IO
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!>   Write a field in the ISF basis in the ETSF format
subroutine write_etsf_density(filename,message,at,rxyz,n1,n2,n3,n1i,n2i,n3i,hxh,hyh,hzh,&
     x)
  use module_base
  use module_types

  implicit none
  character(len=*), intent(in) :: filename,message
  integer, intent(in) :: n1,n2,n3,n1i,n2i,n3i
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: at
  real(wp), dimension(n1i,n2i,n3i), intent(in) :: x
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  !local variables

  write(0, "(A)") "Illegal call to write_etsf_density(), not compiled with ETSF_IO support."
  stop
END SUBROUTINE write_etsf_density


!>   Read a field in the ISF basis in the ETSF format
!!
subroutine read_etsf(filename,geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,rho,&
     nat,rxyz)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: nspin
  integer, intent(out) ::  n1i,n2i,n3i
  real(gp), intent(out) :: hxh,hyh,hzh
  real(dp), dimension(:,:), pointer :: rho
  real(gp), dimension(:,:), pointer, optional :: rxyz
  integer, intent(out), optional ::  nat

  write(0, "(A)") "Illegal call to read_etsf(), not compiled with ETSF_IO support."
  stop
END SUBROUTINE read_etsf


subroutine read_waves_etsf(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz,  & 
     wfd,psi)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3
  real(gp), intent(in) :: hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(orbitals_data), intent(in) :: orbs
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), dimension(3,at%nat), intent(out) :: rxyz_old
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(out) :: psi
  character(len = *), intent(in) :: filename

  stop 'No ETSF support at compilation!'
END SUBROUTINE read_waves_etsf

subroutine write_waves_etsf(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,nat,rxyz,wfd,psi)
  use module_types
  use module_base
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3,nat
  real(gp), intent(in) :: hx,hy,hz
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd
  real(gp), dimension(3,nat), intent(in) :: rxyz
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(in) :: psi
  character(len = *), intent(in) :: filename

  stop 'No ETSF support at compilation!'
END SUBROUTINE write_waves_etsf
