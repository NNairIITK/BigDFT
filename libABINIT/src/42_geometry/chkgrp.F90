!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkgrp
!! NAME chkgrp
!! chkgrp
!!
!! FUNCTION
!! Checks that a set of input symmetries constitutes a group.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nsym = number of symmetry operations
!! symafm = (anti)ferromagnetic part of symmetry operations
!! symrel = 3D matrix containg symmetry operations
!!
!! OUTPUT
!!
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!! TODO
!! List input and output
!! SHOULD ALSO CHECK THE tnons !
!!
!! PARENTS
!!      chkinp,gensymspgr,m_bands_sym,m_bz_mesh,setsym
!!
!! CHILDREN
!!      abi_leave_new,abi_wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

subroutine chkgrp(nsym,symafm,symrel)

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
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)

!Local variables-------------------------------
!scalars
 integer :: ii,isym,jj,jsym,kk,ksym,symafmchk,testeq=1
 character(len=500) :: message
!arrays
 integer :: chk(3,3)

! *************************************************************************

!DEBUG
!write(6,*)' chkgrp : enter'
!write(6,*)'     isym         symrel            symafm '
!do isym=1,nsym
!write(6, '(i3,a,9i3,a,i3)' )isym,'   ',symrel(:,:,isym),'   ',symafm(isym)
!end do
!ENDDEBUG

 do isym=1,nsym
   do jsym=1,nsym

!    Compute the product of the two symmetries
     do ii=1,3
       do jj=1,3
         chk(ii,jj)=0
         do kk=1,3
           chk(ii,jj)=chk(ii,jj)+&
&           symrel(ii,kk,jsym)*symrel(kk,jj,isym)
         end do
       end do
     end do
     symafmchk=symafm(jsym)*symafm(isym)

!    Check that product array is one of original symmetries
     do ksym=1,nsym
       testeq=1
       do ii=1,3
         do jj=1,3
           if(chk(ii,jj)/=symrel(ii,jj,ksym))testeq=0
         end do
       end do
       if(symafmchk/=symafm(ksym))testeq=0
!      The test is positive
       if (testeq==1) exit
     end do

!    The test is positive
     if(testeq==1)exit

     write(message, '(a,a,a,a,2i3,a)' ) ch10,&
&     ' chkgrp : ERROR -',ch10,&
&     '  Error: product of symmetries',isym,jsym,' is not in group.'
     call abi_wrtout(std_out,message,'COLL')
     write(message, '(a,a,a,a,a)' ) &
&     '  This indicates that the input symmetry elements',ch10,&
&     '  do not possess closure under group composition.',ch10,&
&     '  Action : check symrel, symafm and fix them.'
     call abi_wrtout(std_out,message,'COLL')
     call abi_leave_new('COLL')

!    End loop on jsym. Note that an "exit" instruction is present inside the loop
   end do

!  End loop on isym
 end do

end subroutine chkgrp
!!***
