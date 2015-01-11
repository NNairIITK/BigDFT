!!****m* ABINIT/interfaces_41_xc_lowlevel
!! NAME
!! interfaces_41_xc_lowlevel
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/41_xc_lowlevel
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

module interfaces_41_xc_lowlevel

 implicit none

interface
 subroutine drivexc(exc,ixc,npts,nspden,order,rho_updn,vxc,ndvxc,ngr2,nd2vxc,nvxcdgr,&  !Mandatory arguments
  &  dvxc,d2vxc,grho2_updn,vxcgr,exexch,lrho_updn,vxclrho,tau_updn,vxctau)    !Optional arguments
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
 end subroutine drivexc
end interface

interface
 subroutine invcb(rhoarr,rspts,npts)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: npts
  real(dp),intent(in) :: rhoarr(npts)
  real(dp),intent(out) :: rspts(npts)
 end subroutine invcb
end interface

interface
   subroutine mkdenpos(iwarn,nfft,nspden,option,rhonow,xc_denpos)
     use abi_defs_basis
     implicit none
     integer,intent(in) :: nfft,nspden,option
     integer,intent(inout) :: iwarn
     real(dp),intent(in) :: xc_denpos
     !arrays
     real(dp),intent(inout) :: rhonow(nfft,nspden)
   end subroutine mkdenpos
end interface

interface
 subroutine size_dvxc(ixc,ndvxc,ngr2,nd2vxc,nspden,nvxcdgr,order)
  implicit none
  integer, intent(in) :: ixc
  integer, intent(out) :: nd2vxc
  integer, intent(out) :: ndvxc
  integer, intent(out) :: ngr2
  integer, intent(in) :: nspden
  integer, intent(out) :: nvxcdgr
  integer, intent(in) :: order
 end subroutine size_dvxc
end interface

interface
 subroutine xchcth(dvxcdgr,exci,grho2_updn,ixc,npts,nspden,&  
  &  order,rho_updn,vxci)
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
 end subroutine xchcth
end interface

interface
 subroutine xchelu(exc,npt,order,rspts,vxc,dvxc)  ! dvxc is optional
  use abi_defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xchelu
end interface

interface
 subroutine xclb(grho2_updn,npts,nspden,rho_updn,vxci)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  real(dp),intent(in) :: grho2_updn(npts,2*nspden-1)
  real(dp),intent(in) :: rho_updn(npts,nspden)
  real(dp),intent(inout) :: vxci(npts,nspden)
 end subroutine xclb
end interface

interface
 subroutine xcpbe(exci,npts,nspden,option,order,rho_updn,vxci,ndvxci,ngr2,nd2vxci,&  !Mandatory Arguments
  &  d2vxci,dvxcdgr,dvxci,exexch,grho2_updn)                          !Optional Arguments
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
 end subroutine xcpbe
end interface

interface
 subroutine xcpzca(exc,npt,order,rhor,rspts,vxc,&  !Mandatory arguments
  &  dvxc)                            !Optional arguments
  use abi_defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rhor(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xcpzca
end interface

interface
 subroutine xcspol(exc,npts,nspden,order,rspts,vxc,zeta,ndvxc,&  !Mandatory arguments
  &  dvxc)                            !Optional arguments
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
 end subroutine xcspol
end interface

interface
 subroutine xctetr(exc,npt,order,rhor,rspts,vxc,&  !Mandatory arguments
  &  d2vxc,dvxc)                    !Optional arguments
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
 end subroutine xctetr
end interface

interface
 subroutine xcwign(exc,npt,order,rspts,vxc,&  !Mandatory arguments
  &  dvxc)                           !Optional arguments
  use abi_defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xcwign
end interface

interface
 subroutine xcxalp(exc,npt,order,rspts,vxc, dvxc)  ! dvxc is optional
  use abi_defs_basis
  implicit none
  integer,intent(in) :: npt
  integer,intent(in) :: order
  real(dp),intent(out),optional :: dvxc(npt)
  real(dp),intent(out) :: exc(npt)
  real(dp),intent(in) :: rspts(npt)
  real(dp),intent(out) :: vxc(npt)
 end subroutine xcxalp
end interface

end module interfaces_41_xc_lowlevel
!!***
