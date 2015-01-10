!!****f* ABINIT/size_dvxc
!! NAME
!! size_dvxc
!!
!! FUNCTION
!! Give the size of the array dvxc(npts,ndvxc) and the second dimension of the d2vxc(npts,nd2vxc)
!! needed for the allocations depending on the routine which is called from the drivexc routine
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!! This routine has been written from rhohxc (DCA, XG, GMR, MF, GZ)
!!
!! INPUTS
!!  ixc= choice of exchange-correlation scheme
!!  order=gives the maximal derivative of Exc computed.
!!    1=usual value (return exc and vxc)
!!    2=also computes the kernel (return exc,vxc,kxc)
!!   -2=like 2, except (to be described)
!!    3=also computes the derivative of the kernel (return exc,vxc,kxc,k3xc)
!!
!! OUTPUT
!!  ndvxc size of the array dvxc(npts,ndvxc) for allocation
!!  ngr2 size of the array grho2_updn(npts,ngr2) for allocation
!!  nd2vxc size of the array d2vxc(npts,nd2vxc) for allocation
!!  nvxcdgr size of the array dvxcdgr(npts,nvxcdgr) for allocation
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine size_dvxc(ixc,ndvxc,ngr2,nd2vxc,nspden,nvxcdgr,order)

 use defs_basis

 implicit none

!Arguments----------------------
 integer, intent(in) :: ixc,nspden,order
 integer, intent(out) :: ndvxc,ngr2,nd2vxc,nvxcdgr

!Local variables----------------

! *************************************************************************

 nvxcdgr=0
 ndvxc=0
 nd2vxc=0
 ngr2=2*min(nspden,2)-1
 if (order**2 <= 1) then
   ndvxc=0
   nvxcdgr=0
   if (((ixc>=11 .and. ixc<=15) .or. (ixc==23)) .and. ixc/=13) nvxcdgr=3
   if (ixc==16.or.ixc==17.or.ixc==26.or.ixc==27) nvxcdgr=2
   if (ixc<0) nvxcdgr=3
 else
   if (ixc==1 .or. ixc==21 .or. ixc==22 .or. (ixc>=7 .and. ixc<=10) .or. ixc==13) then
!    new Teter fit (4/93) to Ceperley-Alder data, with spin-pol option    !routine xcspol
!    routine xcpbe, with different options (optpbe) and orders (order)
     ndvxc=min(nspden,2)+1
!    if (ixc>=7 .and. ixc<=10) nvxcdgr=3
   else if (ixc>=2 .and. ixc<=6) then
!    Perdew-Zunger fit to Ceperly-Alder data (no spin-pol)                !routine xcpzca
!    Teter fit (4/91) to Ceperley-Alder values (no spin-pol)              !routine xctetr
!    Wigner xc (no spin-pol)           !routine xcwign
!    Hedin-Lundqvist xc (no spin-pol)           !routine xchelu
!    X-alpha (no spin-pol)           !routine xcxalp
     ndvxc=1

   else if (ixc==12) then
!    routine xcpbe, with optpbe=-2 and different orders (order)
     ndvxc=8
     nvxcdgr=3
   else if ((ixc>=11 .and. ixc<=15 .and. ixc/=13) .or. (ixc==23)) then
!    routine xcpbe, with different options (optpbe) and orders (order)
     ndvxc=15
     nvxcdgr=3
     nd2vxc=1  !to be corrected when the calculation of d2vxcar will be implemented for these functionals
   else if(ixc==16 .or. ixc==17 .or. ixc==26 .or. ixc==27 .or. ixc<0 ) then
!    Should be 0
     ndvxc=0
     if (ixc==16 .or. ixc==17 .or. ixc==26 .or. ixc==27) nvxcdgr=2
   else if (ixc<0) then
     ndvxc=0
     nvxcdgr=3
   end if

!  definition of nd2vxc, second dimension of the array for third order derivatives
!  for the non spin polarized LDA case nd2vxc=1
   if (ixc==3) nd2vxc=1 
   if  ((ixc>=7 .and. ixc<=10) .or.(ixc==13)) then
     nd2vxc=3*nspden-2   !second dimension for the array of third order derivative
   end if

 end if

end subroutine size_dvxc
!!***
