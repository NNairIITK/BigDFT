!{\src2tex{textfont=tt}}
!!****f* ABINIT/sym2cart
!! NAME
!! sym2cart
!!
!! FUNCTION
!! Routine called by the program optic
!! Convert to symmetry matrice in cartesian coordinates
!!
!! COPYRIGHT
!! Copyright (C) 2002-2010 ABINIT group (SSharma,MVer,VRecoules,TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!	gprimd(3,3)=dimensional primitive translations for reciprocal space
!!	nsym=number of symmetries in group
!!	rprimd(3,3)=dimensional real space primitive translations (bohr)
!!	symrel(3,3,nsym)=symmetry matrices in terms of real space
!!
!! OUTPUT
!!	symcart(3,3)=symmetry matrice in cartesian coordinates (reals)
!!
!! PARENTS
!!      optic
!!
!! CHILDREN
!!      dgemm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

subroutine sym2cart(gprimd,nsym,rprimd,symrel,symcart)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments -----------------------------------------------
! in
! out
!scalars
 integer,intent(in) :: nsym
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3)
 real(dp),intent(out) :: symcart(3,3,nsym)

!Local variables-------------------------------
!scalars
 integer :: isym
!arrays
 real(dp) :: rsym(3,3),rsymcart(3,3),tmp(3,3)

! *************************************************************************

 do isym=1,nsym
   rsym(:,:) = dble(symrel(:,:,isym))
!  write (*,*) 'rsym = ',rsym
   call dgemm('N','N',3,3,3,one,rprimd,3,rsym,  3,zero,tmp,     3)
   call dgemm('N','N',3,3,3,one,tmp,   3,gprimd,3,zero,rsymcart,3)
!  write (*,*) 'rsymcart = ',rsymcart
   symcart(:,:,isym) = rsymcart(:,:)
 end do

end subroutine sym2cart
!!***
