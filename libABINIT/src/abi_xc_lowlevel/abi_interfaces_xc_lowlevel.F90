!!****m* ABINIT/abi_interfaces_xc_lowlevel
!! NAME
!! abi_interfaces_xc_lowlevel
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/abi_xc_lowlevel
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

module abi_interfaces_xc_lowlevel

 implicit none

interface
 subroutine abi_drivexc(exc,ixc,npts,nspden,order,rho_updn,vxc,ndvxc,ngr2,&
&                       nd2vxc,nvxcdgr,dvxc,d2vxc,grho2_updn,vxcgr,exexch,&
&                       lrho_updn,vxclrho,tau_updn,vxctau)
  use abi_defs_basis
  implicit none
  integer,intent(in),optional :: exexch
  integer,intent(in) :: ixc
  integer,intent(in) :: nd2vxc
  integer,intent(in) :: ndvxc
  integer,intent(in) :: ngr2
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  integer,intent(in) :: nvxcdgr
  integer,intent(in) :: order
  real(dp),intent(out),optional :: d2vxc(npts,nd2vxc)
  real(dp),intent(out),optional :: dvxc(npts,ndvxc)
  real(dp),intent(out) :: exc(npts)
  real(dp),intent(in),optional :: grho2_updn(npts,ngr2)
  real(dp),intent(in),optional :: lrho_updn(npts,nspden)
  real(dp),intent(in) :: rho_updn(npts,nspden)
  real(dp),intent(in),optional :: tau_updn(npts,nspden)
  real(dp),intent(out) :: vxc(npts,nspden)
  real(dp),intent(out),optional :: vxcgr(npts,nvxcdgr)
  real(dp),intent(out),optional :: vxclrho(npts,nspden)
  real(dp),intent(out),optional :: vxctau(npts,nspden)
 end subroutine abi_drivexc
end interface

interface
 subroutine abi_invcb(rhoarr,rspts,npts)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: npts
  real(dp),intent(in) :: rhoarr(npts)
  real(dp),intent(out) :: rspts(npts)
 end subroutine abi_invcb
end interface

interface
   subroutine abi_mkdenpos(iwarn,nfft,nspden,option,rhonow,xc_denpos)
     use abi_defs_basis
     implicit none
     integer,intent(in) :: nfft,nspden,option
     integer,intent(inout) :: iwarn
     real(dp),intent(in) :: xc_denpos
     !arrays
     real(dp),intent(inout) :: rhonow(nfft,nspden)
   end subroutine abi_mkdenpos
end interface

interface
 subroutine abi_size_dvxc(ixc,ndvxc,ngr2,nd2vxc,nspden,nvxcdgr,order)
  implicit none
  integer, intent(in) :: ixc
  integer, intent(out) :: nd2vxc
  integer, intent(out) :: ndvxc
  integer, intent(out) :: ngr2
  integer, intent(in) :: nspden
  integer, intent(out) :: nvxcdgr
  integer, intent(in) :: order
 end subroutine abi_size_dvxc
end interface

interface
 subroutine abi_xchcth(dvxcdgr,exci,grho2_updn,ixc,npts,nspden,order,rho_updn,vxci)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: ixc
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  integer,intent(in) :: order
  real(dp),intent(out) :: dvxcdgr(npts,2)
  real(dp),intent(out) :: exci(npts)
  real(dp),intent(in) :: grho2_updn(npts,2*nspden-1)
  real(dp),intent(in) :: rho_updn(npts,nspden)
  real(dp),intent(out) :: vxci(npts,nspden)
 end subroutine abi_xchcth
end interface

interface
 subroutine abi_xchelu(exc,npt,order,rspts,vxc,dvxc)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine abi_xchelu
end interface

interface
 subroutine abi_xclb(grho2_updn,npts,nspden,rho_updn,vxci)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  real(dp),intent(in) :: grho2_updn(npts,2*nspden-1)
  real(dp),intent(in) :: rho_updn(npts,nspden)
  real(dp),intent(inout) :: vxci(npts,nspden)
 end subroutine abi_xclb
end interface

interface
 subroutine abi_xcpbe(exci,npts,nspden,option,order,rho_updn,vxci,ndvxci,&
&                     ngr2,nd2vxci,d2vxci,dvxcdgr,dvxci,exexch,grho2_updn)
  use abi_defs_basis
  implicit none
  integer,intent(in),optional :: exexch
  integer,intent(in) :: nd2vxci
  integer,intent(in) :: ndvxci
  integer,intent(in) :: ngr2
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: order
  real(dp),intent(out),optional :: d2vxci(npts,nd2vxci)
  real(dp),intent(out),optional :: dvxcdgr(npts,3)
  real(dp),intent(out),optional :: dvxci(npts,ndvxci)
  real(dp),intent(out) :: exci(npts)
  real(dp),intent(in),optional :: grho2_updn(npts,ngr2)
  real(dp),intent(in) :: rho_updn(npts,nspden)
  real(dp),intent(out) :: vxci(npts,nspden)
 end subroutine abi_xcpbe
end interface

interface
 subroutine abi_xcpzca(exc,npt,order,rhor,rspts,vxc,dvxc)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rhor(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine abi_xcpzca
end interface

interface
 subroutine abi_xcspol(exc,npts,nspden,order,rspts,vxc,zeta,ndvxc,dvxc)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: ndvxc
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npts,ndvxc)
  real(dp),intent(out) :: exc(npts)
  real(dp),intent(in) :: rspts(npts)
  real(dp),intent(out) :: vxc(npts,nspden)
  real(dp),intent(in) :: zeta(npts)
 end subroutine abi_xcspol
end interface

interface
 subroutine abi_xctetr(exc,npt,order,rhor,rspts,vxc,d2vxc,dvxc)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: d2vxc(npt)
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rhor(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine abi_xctetr
end interface

interface
 subroutine abi_xcwign(exc,npt,order,rspts,vxc,dvxc)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine abi_xcwign
end interface

interface
 subroutine abi_xcxalp(exc,npt,order,rspts,vxc,dvxc)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine abi_xcxalp
end interface

end module abi_interfaces_xc_lowlevel
!!***
