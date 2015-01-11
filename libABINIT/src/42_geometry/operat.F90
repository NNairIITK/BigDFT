!{\src2tex{textfont=tt}}
!!****f* ABINIT/operat
!! NAME
!! operat
!!
!! FUNCTION
!! Computes the atomic position of all the atoms in the unit cell starting
!! with the symmetry operations and the atoms from the assymetric unit cell
!!
!! COPYRIGHT
!! Copyright (C) 1999-2009 ABINIT group (RC,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  natrd = number of atoms in the assymetric unit cell
!!  natom = total number of atoms (to be checked)
!!  nsym = number of symmetry operations
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,nsym)=symmetry operations in real space in terms
!!   of primitive translations
!!  tnons(3,nsym)=nonsymmorphic translations for symmetry operations
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  At input, for the assymetric unit cell
!!  spinat(3,1:natrd)=spin-magnetization of the atoms
!!  typat(1:natrd)=type integer for each atom in cell
!!  xred(3,1:natrd)=reduced dimensionless atomic coordinates
!!
!!  At output, for the complete unit cell
!!  spinat(3,1:natom)=spin-magnetization of the atoms
!!  typat(1:natom)=type integer for each atom in cell
!!  xred(3,1:natom)=reduced dimensionless atomic coordinates
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine operat(natom,natrd,nsym,spinat,symafm,symrel,tnons,typat,xred)

 use abi_defs_basis
 use abi_interfaces_lowlevel

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,natrd,nsym
!arrays
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)
 integer,intent(inout) :: typat(natom)
 real(dp),intent(in) :: tnons(3,nsym)
 real(dp),intent(inout) :: spinat(3,natom),xred(3,natom)

!Local variables ------------------------------
!scalars
 integer :: curat,flagch,flageq,ii,iij,jj,kk
 real(dp),parameter :: nastyzero=1.0d-8
 character(len=500) :: message
!arrays
 integer :: bcktypat(nsym*natrd)
 real(dp) :: bckat(3),bckspinat(3,nsym*natrd),bckxred(3,nsym*natrd)

! *************************************************************************

 curat=0

!Cycle over all the symmetry operations
 do ii=1,nsym

! Cycle over all the atoms in the assymetric unit cell
  do jj=1,natrd

!  Symmetry operation application
   bckat(:)=matmul(symrel(:,:,ii),xred(:,jj))+tnons(:,ii)

!  Normalization of the coordinates in [0,1)
   do iij=1,3
    do while (bckat(iij)<-nastyzero)
     bckat(iij)=bckat(iij)+1.0d0
    end do
    do while (bckat(iij)>=1.0d0-nastyzero)
     bckat(iij)=bckat(iij)-1.0d0
    end do
   end do

!  Check for duplicate atoms
   flagch=0
   do kk=1,curat
    flageq=0
    if ( abs(bckxred(1,kk)-bckat(1))<nastyzero  .and. &
&    abs(bckxred(2,kk)-bckat(2))<nastyzero  .and. &
&    abs(bckxred(3,kk)-bckat(3))<nastyzero       ) exit
    flagch=flagch+1
   end do

   if (flagch==curat) then
!   Add the obtained atom to the bckxred list
    curat=curat+1
    bckxred(:,curat)=bckat
    bcktypat(curat)=typat(jj)
    bckspinat(:,curat)=spinat(:,jj)*symafm(ii)
   end if

  end do

 end do

 if (curat>natom) then
  write(message, '(a,a,a,a,i3,a,a,i7,a,a,a,a)' ) ch10,&
&  ' operat : ERROR -',ch10,&
&  '  The number of atoms obtained from symmetries, ',curat,ch10,&
&  '  is greater than the input number of atoms, natom=',natom,ch10,&
&  '  This is not allowed.  ',ch10,&
&  '  Action : modify natom or the symmetry data in the input file.'
  call abi_wrtout(06,  message,'COLL')
  call abi_leave_new('COLL')
 end if

 if (curat<natom) then
  write(message, '(a,a,a,a,i3,a,a,i7,a,a,a,a)' ) ch10,&
&  ' abinit : ERROR -',ch10,&
&  '  operat : The number of atoms obtained from symmetries, ',curat,ch10,&
&  '  is lower than the input number of atoms, natom=',natom,ch10,&
&  '  This is not allowed.  ',ch10,&
&  '  Action : modify natom or the symmetry data in the input file.'
  call abi_wrtout(06,  message,'COLL')
  call abi_leave_new('COLL')
 end if

!Assignment of symmetry to xred
 xred(:,1:natom)=bckxred(:,1:natom)
 typat(1:natom)=bcktypat(1:natom)
 spinat(1:3,1:natom)=bckspinat(1:3,1:natom)

end subroutine operat
!!***
