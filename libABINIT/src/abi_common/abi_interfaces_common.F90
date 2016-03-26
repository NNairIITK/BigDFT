!!****m* ABINIT/abi_interfaces_common
!! NAME
!! abi_interfaces_common
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/abi_common
!!
!! COPYRIGHT
!! Copyright (C) 2010 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!!
!! SOURCE

module abi_interfaces_common

 implicit none

interface
 subroutine abi_ewald(iproc,nproc,mpi_comm,eew,gmet,grewtn,natom,ntypat,rmet,typat,ucvol,xred,zion)
  use abi_defs_basis
  use abi_interfaces_lowlevel
  implicit none
  integer,intent(in) :: natom,iproc,nproc,mpi_comm
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: eew
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
  real(dp),intent(out) :: grewtn(3,natom)
 end subroutine abi_ewald
end interface

interface
 subroutine abi_ewald2(iproc,nproc,mpi_comm,gmet,natom,ntypat,rmet,rprimd,stress,&
&                  typat,ucvol,xred,zion)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: natom,iproc,nproc,mpi_comm
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: zion(ntypat)
  real(dp),intent(out) :: stress(6)
 end subroutine abi_ewald2
end interface

interface
 subroutine abi_fconv(fcart,iatfix,iexit,itime,natom,ntime,&
&                 optcell,strfact,strtarget,strten,tolmxf)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: itime
  integer,intent(in) :: natom
  integer,intent(in) :: ntime
  integer,intent(in) :: optcell
  integer,intent(inout) :: iexit
  real(dp),intent(in) :: strfact
  real(dp),intent(in) :: tolmxf
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: fcart(3,natom)
  real(dp),intent(in) :: strtarget(6)
  real(dp),intent(in) :: strten(6)
 end subroutine abi_fconv
end interface

interface
 subroutine abi_prtxvf(fcart,iatfix,iout,natom,prtvel,vel,xcart)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: natom
  integer,intent(in) :: prtvel
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: fcart(3,natom)
  real(dp),intent(in) :: vel(3,natom)
  real(dp),intent(in) :: xcart(3,natom)
 end subroutine abi_prtxvf
end interface

end module abi_interfaces_common
!!***
