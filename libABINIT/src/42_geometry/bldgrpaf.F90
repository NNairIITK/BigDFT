!{\src2tex{textfont=tt}}
!!****f* ABINIT/bldgrpaf
!! NAME bldgrpaf
!! bldgrpaf
!!
!! FUNCTION
!! Yields all the magnetic symmetry operations starting from the symrel and
!!  the magnetic generators
!! It applies all the generators onto symrel, it obtains all the other operations.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2010 ABINIT group (RC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! msym = default number of symmetry operations
!! nogenaf = number of generators, number of operations to be applied onto themselves
!! nsym = number of symmetry operations
!! symafm = (anti)ferromagnetic part of symmetry operations
!! symrel = 3D matrix containg symmetry operations
!! symrel_magn = 3D matrix containg magnetic symmetry generators
!! tnons = 2D matrix containing translations associated
!!
!! OUTPUT
!! symrel = 3D matrix containg symmetry operations
!! tnons = 2D matrix containing translations associated
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!!
!! PARENTS
!! None
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

subroutine bldgrpaf(msym,nogenaf,nsym,symafm,symrel,symrel_magn,tnons)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msym,nogenaf,nsym
!arrays
 integer,intent(in) :: symrel_magn(3,3,3)
 integer,intent(inout) :: symafm(msym),symrel(3,3,msym)
 real(dp),intent(inout) :: tnons(3,msym)

!Local variables ------------------------------
!scalars
 integer :: ii,jj,kk
 real(dp) :: nastyzero
!arrays
 integer :: matrintoper(3,3)
 real(dp) :: matrinttransl(3)

! *************************************************************************

 nastyzero=0.1

!DEBUG
!write(6,*)' bldgrpaf : enter, builds the space group symmetry '
!write(6,*)' bldgrpaf : number of generators : ',nogenaf
!ENDDEBUG


!Find the magnetic generators within the nonmagnetic group operations
 do ii=1,nogenaf
   do jj=1,nsym
     if(sum(abs(symrel(:,:,jj)-symrel_magn(:,:,ii)))==0) then
       symafm(jj)=-1
       exit
     end if
   end do
 end do

 do ii=1,nsym
   do jj=1,nsym
     if (symafm(ii)*symafm(jj)==-1) then

!      Computes the new symmetry opreration, like:
!      !   $ { R1 | v1 }{ R2 | v2 } = { R1.R2 | v1+R1.v2 } $
       matrintoper(:,:) = matmul(symrel(:,:,ii),symrel(:,:,jj))
       matrinttransl(:) = tnons(:,ii)+matmul(symrel(:,:,ii),tnons(:,jj))

       do kk=1,nsym
         if( (sum(abs(symrel(:,:,kk)-matrintoper(:,:)))==0).and.&
&         (abs(tnons(1,kk)-matrinttransl(1))+&
&         abs(tnons(2,kk)-matrinttransl(2))+&
&         abs(tnons(3,kk)-matrinttransl(3)))<nastyzero ) then
           symafm(kk)=-1
           exit
         end if
       end do

     end if
   end do
 end do



end subroutine bldgrpaf
!!***
