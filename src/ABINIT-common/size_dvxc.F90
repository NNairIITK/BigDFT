!!****f* ABINIT/size_dvxc
!! NAME
!! size_dvxc
!!
!! FUNCTION
!! Give the size of the array dvxc(npts,ndvxc) 
!! needed for the allocations depending on the routine which is called from the drivexc routine
!!
!! COPYRIGHT
!! Copyright (C) 1998-2006 ABINIT group (TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.
!! This routine has been written from rhohxc_coll(DCA, XG, GMR, MF, GZ)
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
!!
!! PARENTS
!!      pawxc,pawxcm,rhohxc_coll
!!
!! CHILDREN
!!      
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine size_dvxc(ixc,ndvxc,nspden,nvxcdgr,order)

 implicit none

!Arguments----------------------
 integer, intent(in) :: ixc,nspden,order
 integer, intent(out) :: ndvxc,nvxcdgr

!Local variables----------------
 nvxcdgr=0
 ndvxc=0
 if (order**2 <= 1) then
    ndvxc=0
    nvxcdgr=0
    if (ixc>=11 .and. ixc<=15 .and. ixc/=13) nvxcdgr=3 
    if (ixc==16) nvxcdgr=2
 else
    if (ixc==1 .or. ixc==21 .or. ixc==22 .or. (ixc>=7 .and. ixc<=10) .or. ixc==13) then
       ! new Teter fit (4/93) to Ceperley-Alder data, with spin-pol option    !routine xcspol
       !routine xcpbe, with different options (optpbe) and orders (order)
       ndvxc=nspden+1
       !if (ixc>=7 .and. ixc<=10) nvxcdgr=3
    else if (ixc>=2 .and. ixc<=6) then
       ! Perdew-Zunger fit to Ceperly-Alder data (no spin-pol)                !routine xcpzca
       ! Teter fit (4/91) to Ceperley-Alder values (no spin-pol)              !routine xctetr
       ! Wigner xc (no spin-pol)           !routine xcwign
       ! Hedin-Lundqvist xc (no spin-pol)           !routine xchelu
       ! X-alpha (no spin-pol)           !routine xcxalp
       ndvxc=1
      
    else if (ixc==12) then
       !routine xcpbe, with optpbe=-2 and different orders (order)
       ndvxc=8
       nvxcdgr=3
    else if (ixc>=11 .and. ixc<=15 .and. ixc/=13) then
       !routine xcpbe, with different options (optpbe) and orders (order)
       ndvxc=15
       nvxcdgr=3
    else if(ixc==16 .or. ((ixc>=30).and.(ixc<=34)) ) then
       !Should be 0
       ndvxc=0
       if (ixc==16) nvxcdgr=2
    end if
 end if

end subroutine size_dvxc
!!***
