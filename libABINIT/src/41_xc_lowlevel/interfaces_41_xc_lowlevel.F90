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
!! To do that: config/scripts/abilint . .
!! 
!!
!! SOURCE

module interfaces_41_xc_lowlevel

 implicit none

interface
 subroutine drivexc(exc,ixc,npts,nspden,order,rho_updn,vxc,ndvxc,ngr2,nd2vxc,nvxcdgr,&  !Mandatory arguments
  &  dvxc,d2vxc,grho2_updn,vxcgr,exexch,lrho_updn,vxclrho,tau_updn,vxctau)    !Optional arguments
  use defs_basis
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
 subroutine echo_xc_name (ixc)
  implicit none
  integer, intent(in) :: ixc
 end subroutine echo_xc_name
end interface

interface
 subroutine gammapositron(gamma,grhocore2,grhoe2,igamma,ngr,npt,rhocore,rhoer,rhopr,usecore)
  use defs_basis
  implicit none
  integer,intent(in) :: igamma
  integer,intent(in) :: ngr
  integer,intent(in) :: npt
  integer,intent(in) :: usecore
  real(dp),intent(out) :: gamma(npt,2)
  real(dp),intent(in) :: grhocore2(ngr*usecore)
  real(dp),intent(in) :: grhoe2(ngr)
  real(dp),intent(in) :: rhocore(npt*usecore)
  real(dp),intent(in) :: rhoer(npt)
  real(dp),intent(in) :: rhopr(npt)
 end subroutine gammapositron
end interface

interface
 subroutine hartre(cplex,gmet,gsqcut,izero,mpi_enreg,nfft,ngfft,paral_kgb,qphon,rhog,vhartr)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: izero
  integer,intent(in) :: nfft
  integer,intent(in) :: paral_kgb
  real(dp),intent(in) :: gsqcut
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhog(2,nfft)
  real(dp),intent(out) :: vhartr(cplex*nfft)
 end subroutine hartre
end interface

interface
 subroutine invcb(rhoarr,rspts,npts)
  use defs_basis
  implicit none
  integer,intent(in) :: npts
  real(dp),intent(in) :: rhoarr(npts)
  real(dp),intent(out) :: rspts(npts)
 end subroutine invcb
end interface

interface
 subroutine mkcore(corstr,dyfrx2,grxc,mpi_enreg,natom,nfft,nspden,ntypat,n1,n1xccc,&  
  &  n2,n3,option,rprimd,typat,ucvol,vxc,xcccrc,xccc1d,xccc3d,xred)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: n1
  integer,intent(in) :: n1xccc
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: option
  type(mpi_type) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: corstr(6)
  real(dp),intent(out) :: dyfrx2(3,3,natom)
  real(dp),intent(out) :: grxc(3,natom)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: vxc(nfft,nspden)
  real(dp),intent(in) :: xccc1d(n1xccc,6,ntypat)
  real(dp),intent(inout) :: xccc3d(nfft)
  real(dp),intent(in) :: xcccrc(ntypat)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mkcore
end interface

interface
   subroutine mkdenpos(iwarn,nfft,nspden,option,rhonow,xc_denpos)
     use defs_basis
     implicit none
     integer,intent(in) :: nfft,nspden,option
     integer,intent(inout) :: iwarn
     real(dp),intent(in) :: xc_denpos
     !arrays
     real(dp),intent(inout) :: rhonow(nfft,nspden)
   end subroutine mkdenpos
end interface

interface
 subroutine mkvxc3(cplex,kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,n3xccc,option,&  
  &  paral_kgb,qphon,rhor1,rprimd,vxc1,xccc3d1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: option
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhor1(cplex*nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vxc1(cplex*nfft,nspden)
  real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
 end subroutine mkvxc3
end interface

interface
 subroutine mkvxcgga3(cplex,gprimd,kxc,mpi_enreg,nfft,ngfft,nkxc,&  
  &  nspden,paral_kgb,qphon,rhor1tmp,vxc1)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nkxc
  integer,intent(in) :: nspden
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhor1tmp(cplex*nfft,2)
  real(dp),intent(out) :: vxc1(cplex*nfft,nspden)
 end subroutine mkvxcgga3
end interface

interface
 subroutine phase(ngfft,ph)
  use defs_basis
  implicit none
  integer,intent(in) :: ngfft
  real(dp),intent(out) :: ph(2*ngfft)
 end subroutine phase
end interface

!!$interface
!!$ subroutine rhohxc(dtset,enxc,gsqcut,izero,kxc,mpi_enreg,nfft,ngfft,&  
!!$  &  nhat,nhatdim,nhatgr,nhatgrdim,nkxc,nk3xc,nspden,n3xccc,option,rhog,rhor,rprimd,&  
!!$  &  strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d,k3xc,&  
!!$  &  electronpositron,taug,taur) ! optional argument
!!$  use defs_basis
!!$  use defs_abitypes
!!$  use m_electronpositron
!!$  implicit none
!!$  integer,intent(in) :: izero
!!$  integer,intent(in) :: n3xccc
!!$  integer,intent(in) :: nfft
!!$  integer,intent(in) :: nhatdim
!!$  integer,intent(in) :: nhatgrdim
!!$  integer,intent(in) :: nk3xc
!!$  integer,intent(in) :: nkxc
!!$  integer,intent(in) :: nspden
!!$  integer,intent(in) :: option
!!$  integer,intent(in) :: usexcnhat
!!$  type(dataset_type),intent(in) :: dtset
!!$  type(electronpositron_type),pointer,optional :: electronpositron
!!$  real(dp),intent(out) :: enxc
!!$  real(dp),intent(in) :: gsqcut
!!$  type(mpi_type),intent(inout) :: mpi_enreg
!!$  real(dp),intent(out) :: vxcavg
!!$  integer,intent(in) :: ngfft(18)
!!$  real(dp),intent(out),optional :: k3xc(1:nfft,1:nk3xc)
!!$  real(dp),intent(out) :: kxc(nfft,nkxc)
!!$  real(dp),intent(in) :: nhat(nfft,nspden*nhatdim)
!!$  real(dp),intent(in) :: nhatgr(nfft,nspden,3*nhatgrdim)
!!$  real(dp),intent(in) :: rhog(2,nfft)
!!$  real(dp),intent(in) :: rhor(nfft,nspden)
!!$  real(dp),intent(in) :: rprimd(3,3)
!!$  real(dp),intent(out) :: strsxc(6)
!!$  real(dp),intent(in),optional :: taug(2,nfft*dtset%usekden)
!!$  real(dp),intent(in),optional :: taur(nfft,nspden*dtset%usekden)
!!$  real(dp),intent(out) :: vhartr(nfft)
!!$  real(dp),intent(out) :: vxc(nfft,nspden)
!!$  real(dp),intent(in) :: xccc3d(n3xccc)
!!$ end subroutine rhohxc
!!$end interface
!!$
!!$interface
!!$ subroutine rhohxcpositron(electronpositron,gprimd,kxcapn,mpi_enreg,nfft,ngfft,nkxc,nspden,n3xccc,&  
!!$  &  paral_kgb,rhor,strsxc,ucvol,vhartr,vxcapn,vxcavg,xccc3d)
!!$  use defs_basis
!!$  use defs_abitypes
!!$  use m_electronpositron
!!$  implicit none
!!$  integer,intent(in) :: n3xccc
!!$  integer,intent(in) :: nfft
!!$  integer,intent(in) :: nkxc
!!$  integer,intent(in) :: nspden
!!$  integer,intent(in) :: paral_kgb
!!$  type(electronpositron_type),pointer :: electronpositron
!!$  type(mpi_type),intent(inout) :: mpi_enreg
!!$  real(dp),intent(in) :: ucvol
!!$  real(dp),intent(out) :: vxcavg
!!$  integer,intent(in) :: ngfft(18)
!!$  real(dp),intent(in) :: gprimd(3,3)
!!$  real(dp),intent(out) :: kxcapn(nfft,nkxc)
!!$  real(dp),intent(in) :: rhor(nfft,nspden)
!!$  real(dp),intent(out) :: strsxc(6)
!!$  real(dp),intent(out) :: vhartr(nfft)
!!$  real(dp),intent(out) :: vxcapn(nfft,nspden)
!!$  real(dp),intent(in) :: xccc3d(n3xccc)
!!$ end subroutine rhohxcpositron
!!$end interface

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
 subroutine xc_kernel(Dtset,ixc,MPI_enreg,ngfft,nfft,nsppol,rhor,rprimd,npw,dim_kxcg,kxcg,gvec)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: dim_kxcg
  integer,intent(in) :: ixc
  integer,intent(in) :: nfft
  integer,intent(in) :: npw
  integer,intent(in) :: nsppol
  type(dataset_type),intent(in) :: Dtset
  type(mpi_type),intent(inout) :: MPI_enreg
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: gvec(3,npw)
  complex(gwpc),intent(out) :: kxcg(nfft,dim_kxcg)
  real(dp),intent(in) :: rhor(nfft,nsppol)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine xc_kernel
end interface

interface
 subroutine xcden (cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden,paral_kgb,qphon,rhor,rhonow,&  !Mandatory arguments
  lrhonow)              !Optional arguments
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: ishift
  integer,intent(in) :: nfft
  integer,intent(in) :: ngrad
  integer,intent(in) :: nspden
  integer,intent(in) :: paral_kgb
  type(mpi_type) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out),optional :: lrhonow(cplex*nfft,nspden)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(out) :: rhonow(cplex*nfft,nspden,ngrad*ngrad)
  real(dp),intent(in) :: rhor(cplex*nfft,nspden)
 end subroutine xcden
end interface

interface
 subroutine xchcth(dvxcdgr,exci,grho2_updn,ixc,npts,nspden,&  
  &  order,rho_updn,vxci)
  use defs_basis
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
  use defs_basis
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
  use defs_basis
  implicit none
  integer,intent(in) :: npts
  integer,intent(in) :: nspden
  real(dp),intent(in) :: grho2_updn(npts,2*nspden-1)
  real(dp),intent(in) :: rho_updn(npts,nspden)
  real(dp),intent(inout) :: vxci(npts,nspden)
 end subroutine xclb
end interface

interface
 subroutine xcmult (dnexcdn,nfft,ngrad,nspden,nspgrad,rhonow)
  use defs_basis
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: ngrad
  integer,intent(in) :: nspden
  integer,intent(in) :: nspgrad
  real(dp),intent(in) :: dnexcdn(nfft,nspgrad)
  real(dp),intent(inout) :: rhonow(nfft,nspden,ngrad*ngrad)
 end subroutine xcmult
end interface

interface
 subroutine xcpbe(exci,npts,nspden,option,order,rho_updn,vxci,ndvxci,ngr2,nd2vxci,&  !Mandatory Arguments
  &  d2vxci,dvxcdgr,dvxci,exexch,grho2_updn)                          !Optional Arguments
  use defs_basis
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
 subroutine xcpositron(fnxc,grhoe2,ixcpositron,ngr,npt,posdensity0_limit,rhoer,rhopr,vxce,vxcegr,vxcp,&  
  &  dvxce,dvxcp) ! optional arguments
  use defs_basis
  implicit none
  integer,intent(in) :: ixcpositron
  integer,intent(in) :: ngr
  integer,intent(in) :: npt
  logical,intent(in) :: posdensity0_limit
  real(dp),intent(out),optional :: dvxce(npt)
  real(dp),intent(out),optional :: dvxcp(npt)
  real(dp),intent(out) :: fnxc(npt)
  real(dp),intent(in) :: grhoe2(ngr)
  real(dp),intent(in) :: rhoer(npt)
  real(dp),intent(in) :: rhopr(npt)
  real(dp),intent(out) :: vxce(npt)
  real(dp),intent(out) :: vxcegr(ngr)
  real(dp),intent(out) :: vxcp(npt)
 end subroutine xcpositron
end interface

interface
 subroutine xcpot (cplex,dnexcdn,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden,&  
  &  nspgrad,paral_kgb,qphon,rhonow,vxc)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: ishift
  integer,intent(in) :: nfft
  integer,intent(in) :: ngrad
  integer,intent(in) :: nspden
  integer,intent(in) :: nspgrad
  integer,intent(in) :: paral_kgb
  type(mpi_type) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: dnexcdn(cplex*nfft,nspgrad)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rhonow(cplex*nfft,nspden,ngrad*ngrad)
  real(dp),intent(out) :: vxc(cplex*nfft,nspden)
 end subroutine xcpot
end interface

interface
 subroutine xcpzca(exc,npt,order,rhor,rspts,vxc,&  !Mandatory arguments
  &  dvxc)                            !Optional arguments
  use defs_basis
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
  use defs_basis
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
  use defs_basis
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
  use defs_basis
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
  use defs_basis
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
