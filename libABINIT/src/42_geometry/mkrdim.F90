!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkrdim
!! NAME
!! mkrdim
!!
!! FUNCTION
!!  Trivial subroutine to make dimensional real space
!!  primitive translations from length scales acell(3)
!!  and dimensionless translations rprim(3,3).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  acell(3)=unit cell length scales (bohr)
!!  rprim(3,3)=dimensionless real space primitive translations
!!
!! OUTPUT
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!              where: rprimd(i,j)=rprim(i,j)*acell(j)
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mkrdim(acell,rprim,rprimd)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3)
 real(dp),intent(out) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,jj

! *************************************************************************

 do ii=1,3
   do jj=1,3
     rprimd(ii,jj)=rprim(ii,jj)*acell(jj)
   end do
 end do

end subroutine mkrdim
!!***
