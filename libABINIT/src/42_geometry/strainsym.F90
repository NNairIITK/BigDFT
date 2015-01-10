!{\src2tex{textfont=tt}}
!!****f* ABINIT/strainsym
!! NAME
!! strainsym
!!
!! FUNCTION
!! For given order of point group, symmetrizes the strain tensor,
!! then produce primitive vectors based on the symmetrized strain.
!!
!! COPYRIGHT
!! Copyright (C) 2007-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nsym=order of group.
!! rprimd(3,3)= primitive vectors, to be symmetrized
!! rprimd0(3,3)= reference primitive vectors, already symmetrized
!! symrel(3,3,nsym)=symmetry operators in terms of action on primitive translations
!!
!! OUTPUT
!! rprimd_symm(3,3)= symmetrized primitive vectors
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine strainsym(nsym,rprimd0,rprimd,rprimd_symm,symrel)

 use defs_basis
 use abi_interfaces_linalg
 use interfaces_32_util

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: rprimd(3,3),rprimd0(3,3)
 real(dp),intent(out) :: rprimd_symm(3,3)

!Local variables-------------------------------
!scalars
 integer :: isym
!arrays
 integer :: symrel_it(3,3)
 real(dp) :: rprimd0_inv(3,3),strain(3,3),strain_symm(3,3),tmp_mat(3,3)

!**************************************************************************

!copy initial rprimd input and construct inverse
 rprimd0_inv = rprimd0
 call matrginv(rprimd0_inv,3,3)

!define strain as rprimd = strain * rprimd0 (in cartesian frame)
!so strain = rprimd * rprimd0^{-1}
!transform to triclinic frame with rprimd0^{-1} * strain * rprimd0
!giving strain as rprimd0^{-1} * rprimd
 call dgemm('N','N',3,3,3,one,rprimd0_inv,3,rprimd,3,zero,strain,3)

!loop over symmetry elements to obtain symmetrized strain matrix
 strain_symm = zero
 do isym = 1, nsym

!  this loop accumulates symrel^{-1}*strain*symrel into strain_symm

!  mati3inv gives the inverse transpose of symrel
   call mati3inv(symrel(:,:,isym),symrel_it)
   call dgemm('N','N',3,3,3,one,strain,3,dble(symrel(:,:,isym)),3,zero,tmp_mat,3)
   call dgemm('T','N',3,3,3,one,dble(symrel_it),3,tmp_mat,3,one,strain_symm,3)

 end do

!normalize by number of symmetry operations
 strain_symm = strain_symm/dble(nsym)

!this step is equivalent to r_new = r_old * strain * r_old^{-1} * r_old,
!that is, convert strain back to cartesian frame and then multipy by r_old,
!to get the r_new primitive vectors

 call dgemm('N','N',3,3,3,one,rprimd0,3,strain_symm,3,zero,rprimd_symm,3)

end subroutine strainsym
!!***
