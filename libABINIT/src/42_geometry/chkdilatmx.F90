!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkdilatmx
!! NAME
!! chkdilatmx
!!
!! FUNCTION
!! Check whether the new rprimd does not give a too large number
!! of plane waves, compared to the one booked for rprimd, taking
!! into account the maximal dilatation dilatmx. Actually check whether
!! the new Fermi sphere is inside the old one, dilated.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dilatmx     = maximal dilatation factor (usually the input variable)
!!  rprimd      = new primitive vectors
!!  rprimd_orig = original primitive vectors (usually the input variable)
!!
!! OUTPUT
!!
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!!
!! PARENTS
!!      driver,scfcv
!!
!! CHILDREN
!!      leave_new,matr3eigval,matr3inv,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

subroutine chkdilatmx(dilatmx,rprimd,rprimd_orig)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: dilatmx
!arrays
 real(dp),intent(in) :: rprimd(3,3),rprimd_orig(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,mu
 real(dp) :: dilatmx_new
 character(len=500) :: message
!arrays
 real(dp) :: eigval(3),gprimd_orig(3,3),met(3,3),old_to_new(3,3)

! *************************************************************************

!DEBUG
!write(6,*)' chkdilatmx : enter '
!write(6,*)' rprimd_orig=',rprimd_orig
!write(6,*)' rprimd=',rprimd
!ENDDEBUG

!Generates gprimd
 call matr3inv(rprimd_orig,gprimd_orig)

!Find the matrix that transform an original xcart to xred, then
!to the new xcart
 do mu=1,3
   old_to_new(mu,:)=rprimd(mu,1)*gprimd_orig(:,1)+&
&   rprimd(mu,2)*gprimd_orig(:,2)+&
&   rprimd(mu,3)*gprimd_orig(:,3)
 end do

!The largest increase in length will be obtained thanks
!to the diagonalization of the corresponding metric matrix :
!it is the square root of its largest eigenvalue.
 do ii=1,3
   do jj=1,3
     met(ii,jj)=old_to_new(1,ii)*old_to_new(1,jj)+&
&     old_to_new(2,ii)*old_to_new(2,jj)+&
&     old_to_new(3,ii)*old_to_new(3,jj)
   end do
 end do
!DEBUG
!write(6,*)' met=',met
!ENDDEBUG
 call matr3eigval(eigval,met)

 dilatmx_new=sqrt(maxval(eigval(:)))

 if(dilatmx_new>dilatmx+tol6)then
   write(message,'(10a,es16.6,2a)') ch10,&
&   ' chkdilatmx: ERROR -',ch10,&
&   '  The new primitive vectors rprimd (an evolving quantity)',ch10,&
&   '  are too large with respect to the old rprimd and the accompanying dilatmx :',ch10,&
&   '  this large change of unit cell parameters is not allowed by the present value of dilatmx.',ch10,&
&   '  You need at least dilatmx=',dilatmx_new+tol6,ch10,&
&   '  Action : increase the input variable dilatmx.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!DEBUG
!write(6,*)' chkdilatmx : exit'
!write(6,*)' eigval=',eigval
!stop
!ENDDEBUG

end subroutine chkdilatmx
!!***
