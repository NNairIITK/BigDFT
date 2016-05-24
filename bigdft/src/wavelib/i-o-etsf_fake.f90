!> @file
!!  Fake routines for ETSF-IO
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 



!> Read wavefunctions using ETSF format
subroutine read_waves_etsf(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz,  & 
     wfd,psi)
  use module_defs, only: wp,gp,dp
  use module_types
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3
  real(gp), intent(in) :: hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(orbitals_data), intent(in) :: orbs
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(out) :: psi
  character(len = *), intent(in) :: filename

  stop 'No ETSF support at compilation!'

  !To avoid warnings from the compiler
  write(*,*) filename,iproc,n1,n2,n3,hx,hy,hz,rxyz(1,1)
  rxyz_old=0.0_gp
  psi(1,1)=0.0_dp
END SUBROUTINE read_waves_etsf


subroutine read_one_wave_etsf(iproc,filename,iorbp,isorb,nspinor,n1,n2,n3,&
     & hx,hy,hz,at,rxyz_old,rxyz,wfd,psi,eval)
  use PSbase
  use module_types
  implicit none
  integer, intent(in) :: iorbp,iproc,n1,n2,n3,nspinor,isorb
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(atoms_data), intent(in) :: at
  real(gp), intent(in) :: hx,hy,hz
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(dp), intent(out) :: eval
  real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
  real(dp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,nspinor), intent(out) :: psi
  character(len = *), intent(in) :: filename
  !To avoid warnings from the compiler
  rxyz_old = rxyz
  psi = 0.0_dp
  eval = 0.0_dp

  stop 'No ETSF support at compilation!'

  !To avoid warning from the compiler
  write (*,*) filename,wfd%nvctr_c,at%astruct%nat,hx,hy,hz
END SUBROUTINE read_one_wave_etsf

!> Write wavefunctions in ETSF format
subroutine write_waves_etsf(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,nat,rxyz,wfd,psi)
  use PSbase
  use module_types
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3,nat
  real(gp), intent(in) :: hx,hy,hz
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd
  real(gp), dimension(3,nat), intent(in) :: rxyz
  real(dp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(in) :: psi
  character(len = *), intent(in) :: filename

  stop 'No ETSF support at compilation!'

  !To avoid warnings from the compiler
  write(*,*) iproc,filename,n1,n2,n3,nat,hx,hy,hz,rxyz(1,1),psi(1,1)
END SUBROUTINE write_waves_etsf


subroutine read_psi_compress_etsf(ncid, iorbp, orbs, nvctr, wfd, psi, orblist)
  use PSbase
  use module_types
  implicit none
  integer, intent(in) :: iorbp, ncid
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(orbitals_data), intent(in) :: orbs
  integer, dimension(wfd%nvctr_c), intent(in) :: nvctr
  real(dp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(out) :: psi
  integer, dimension(orbs%norb), intent(in) :: orblist
  !To avoid warning from the compiler
  psi = 0.0_dp

  stop 'No ETSF support at compilation!'

  !To avoid warning from the compiler
  write(*,*) iorbp, ncid, wfd%nvctr_c, orbs%norb, nvctr, orblist
END SUBROUTINE read_psi_compress_etsf


subroutine read_psi_full_etsf(ncid, iorbp, orbs, n1, n2, n3, &
     & nvctr_c, nvctr, gcoord, psig, orblist)
  use PSbase
  use module_types
  implicit none
  integer, intent(in) :: iorbp, n1, n2, n3, nvctr_c, ncid
  type(orbitals_data), intent(in) :: orbs
  real(dp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(out) :: psig
  integer, dimension(3,nvctr_c), intent(in) :: gcoord
  integer, dimension(orbs%norb), intent(in) :: orblist
  integer, dimension(nvctr_c), intent(in) :: nvctr
  !To avoid warning from the compiler
  psig = 0.0_dp

  stop 'No ETSF support at compilation!'

  !To avoid warning from the compiler
  write(*,*) iorbp, n1, n2, n3, nvctr_c, ncid
  write(*,*) orbs%norb, gcoord, orblist, nvctr
END SUBROUTINE read_psi_full_etsf


subroutine write_psi_compress_etsf(ncid, iorbp, orbs, nvctr, wfd, psi)
  use PSbase
  use module_types
  implicit none
  integer, intent(in) :: iorbp, ncid
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(orbitals_data), intent(in) :: orbs
  integer, dimension(wfd%nvctr_c), intent(in) :: nvctr
  real(dp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(in) :: psi

  stop 'No ETSF support at compilation!'

  !To avoid warning from the compiler
  write(*,*) iorbp, ncid, wfd%nvctr_c,orbs%norb,nvctr, psi(1)
END SUBROUTINE write_psi_compress_etsf


subroutine readwavetoisf_etsf(lstat, filename, iorbp, hx, hy, hz, &
     & n1, n2, n3, nspinor, psiscf)
  use PSbase
  implicit none
  character(len = *), intent(in) :: filename
  integer, intent(in) :: iorbp
  integer, intent(out) :: n1, n2, n3, nspinor
  real(gp), intent(out) :: hx, hy, hz
  real(dp),  dimension(:,:,:,:), pointer, intent(out) :: psiscf
  logical, intent(out) :: lstat
  !To avoid warning from the compiler
  n1 = 0
  n2 = 0
  n3 = 0
  nspinor = 0
  hx = 0.0_gp
  hy = 0.0_gp
  hz = 0.0_gp
  lstat = .false.
  nullify(psiscf)

  stop 'No ETSF support at compilation!'

  !To avoid warning from the compiler
  write(*,*) filename, iorbp
END SUBROUTINE readwavetoisf_etsf


subroutine readwavedescr_etsf(lstat, filename, norbu, norbd, nkpt, nspinor)
  use PSbase
  implicit none
  character(len = *), intent(in) :: filename
  integer, intent(out) :: norbu, norbd, nkpt, nspinor
  logical, intent(out) :: lstat
  !To avoid warning from the compiler
  norbu = 0
  norbd = 0 
  nkpt = 0
  nspinor = 0
  lstat = .false.

  stop 'No ETSF support at compilation!'

  !To avoid warning from the compiler
  write(*,*) filename
END SUBROUTINE readwavedescr_etsf


subroutine read_pw_waves(filename, iproc, nproc, at, rxyz, Glr, orbs, psig, rhoij)
  use module_defs, only: gp, wp
  use module_atoms
  use locregs
  use module_types
  implicit none
  character(len = *), intent(in) :: filename
  integer, intent(in) :: iproc, nproc
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  type(locreg_descriptors), intent(in) :: Glr
  type(orbitals_data), intent(in) :: orbs
  real(wp), dimension(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, orbs%norbp), intent(out) :: psig
  real(wp), dimension(:,:,:), pointer, optional :: rhoij

  psig(1,1) = 0._wp

  stop 'No ETSF support at compilation!'

  !To avoid warning from the compiler
  write(*,*) filename, iproc, nproc, at%astruct%nat, rxyz(1,1), Glr%d%n1, orbs%norbp
  if (present(rhoij)) then
     write(*,*) associated(rhoij)
  end if
END SUBROUTINE read_pw_waves
