!{\src2tex{textfont=tt}}
!!****f* ABINIT/symdet
!! NAME
!! symdet
!!
!! FUNCTION
!! Compute determinant of each input symmetry matrix sym(3,3,i)
!! and check that the determinant is always +/- 1.  Integer arithmetic.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nsym=number of symmetry operations
!! sym(3,3,nsym)=integer symmetry array
!!
!! OUTPUT
!! determinant(nsym)=determinant of each symmetry operation
!!
!! PARENTS
!!      remove_inversion,setsym,symanal,symspgr
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

subroutine symdet(determinant,nsym,sym)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
!arrays
 integer,intent(in) :: sym(3,3,nsym)
 integer,intent(out) :: determinant(nsym)

!Local variables-------------------------------
!scalars
 integer :: det,isym
 character(len=500) :: message

! *************************************************************************

 do isym=1,nsym
   det=sym(1,1,isym)*sym(2,2,isym)*sym(3,3,isym)+&
&   sym(2,1,isym)*sym(3,2,isym)*sym(1,3,isym)+&
&   sym(1,2,isym)*sym(2,3,isym)*sym(3,1,isym) - &
&   (sym(3,1,isym)*sym(2,2,isym)*sym(1,3,isym)+&
&   sym(2,1,isym)*sym(1,2,isym)*sym(3,3,isym)+&
&   sym(3,2,isym)*sym(2,3,isym)*sym(1,1,isym))
   if (abs(det)/=1) then
     write(message, '(a,a,a,a,i5,a,i10,a,a,a,a,a)' ) ch10,&
&     ' symdet: ERROR -',ch10,&
&     '  Abs(determinant) for symmetry number',isym,&
&     ' is',det,' .',ch10,&
&     '  For a legitimate symmetry, abs(determinant) must be 1.',ch10,&
&     '  Action : check your symmetry operations (symrel) in input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
   determinant(isym)=det
 end do

end subroutine symdet
!!***
