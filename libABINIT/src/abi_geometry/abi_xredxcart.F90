!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_xredxcart
!! NAME
!! abi_xredxcart
!!
!! FUNCTION
!! If option==1 :
!! Convert from dimensionless reduced coordinates xred(3,natom)
!! to cartesian coordinates xcart(3,natom) in bohr by using
!! xcart(mu,ia)=rprimd(mu,1)*xred(1,ia)
!!             +rprimd(mu,2)*xred(2,ia)
!!             +rprimd(mu,3)*xred(3,ia)
!!
!! If option==-1
!! Convert from cartesian coordinates xcart(3,natom) in bohr to
!! dimensionless reduced coordinates xred(3,natom) by using
!! xred(mu,ia)=(gprimd(1,mu)*xcart(1,ia)+gprimd(2,mu)*xcart(2,ia)+
!!       gprimd(3,mu)*xcart(3,ia))
!! where gprimd is the inverse of rprimd
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  natom=number of atoms in unit cell
!!  option=see above
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output (see above):
!!  xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!  xred(3,natom)=dimensionless reduced coordinates of atoms
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine abi_xredxcart(natom,option,rprimd,xcart,xred)

 use abi_defs_basis
 use abi_interfaces_lowlevel
 use abi_interfaces_numeric

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,option
!arrays
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: xcart(3,natom),xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,mu
 character(len=500) :: message
!arrays
 real(dp) :: gprimd(3,3)

! *************************************************************************

 if(option==1)then
   do iatom=1,natom
     do mu=1,3
       xcart(mu,iatom)=rprimd(mu,1)*xred(1,iatom)+rprimd(mu,2)*xred(2,iatom)+&
&       rprimd(mu,3)*xred(3,iatom)
     end do
   end do
 else if(option==-1)then
   call abi_matr3inv(rprimd,gprimd)
   do iatom=1,natom
     do mu=1,3
       xred(mu,iatom)= gprimd(1,mu)*xcart(1,iatom)+gprimd(2,mu)*xcart(2,iatom)+&
&       gprimd(3,mu)*xcart(3,iatom)
     end do
   end do
 else
   write(message, '(a,a,a,a,i4,a)' ) ch10,&
&   ' abi_xredxcart : BUG -',ch10,&
&   '  Option must be 1 or -1, while it is ',option,'.'
   call abi_wrtout(std_out,message,'COLL')
   call abi_leave_new('COLL')
 end if

end subroutine abi_xredxcart
!!***
