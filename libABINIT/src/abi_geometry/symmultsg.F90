!{\src2tex{textfont=tt}}
!!****f* ABINIT/symmultsg
!! NAME symmultsg
!! symmultsg
!!
!!
!! FUNCTION
!! Yields all the symmetry operations starting from the generators.
!! Applies all the generators onto themselves, and obtains all the other operations.
!! Iterates until it reaches nsym.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2010 ABINIT group (RC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nsym = number of symmetry operations
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym) = 3D matrix containg symmetry operations
!! tnons(3,nsym) = 2D matrix containing translations associated
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine symmultsg(nsym,symafm,symrel,tnons)

 use abi_defs_basis
 use abi_interfaces_lowlevel

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
!arrays
 integer :: symafm(nsym),symrel(3,3,nsym)
 real(dp) :: tnons(3,nsym)

!Local variables ------------------------------
!matrintoper(3,3) & matrinttransl(3) are intermediate arrays of the new
!      symmetry operations obtained, in order to check their uniqueness.
!flagop,flagtr = flags used during the checking of the similarity between
!      the obtained operation and the already existent ones
!ii,ijk,ijkl,jjj,kk = counters in the cycles
!scalars
 integer :: flagma,flagop,flagtr,ii,isym,jsym,matrintsymafm
 real(dp) :: nastyzero
 character(len=500) :: message
!arrays
 integer :: bcksymafm(nsym),bcksymrel(3,3,nsym),matrintoper(3,3)
 integer :: symequiv(nsym,nsym)
 real(dp) :: bcktnons(3,2*nsym),matrinttransl(3)

! *************************************************************************

 nastyzero=0.1

!Transfer the generators to bcksymrel
 do ii=1,nsym
   bcksymrel(:,:,ii)=symrel(:,:,ii)
   bcktnons(:,ii)=tnons(:,ii)
   bcksymafm(ii)=symafm(ii)
 end do

 symequiv(:,:)=0

!Simply iterate until the group is complete
 do isym=1,nsym            ! loop over symmetries
   do jsym=1,nsym           ! loop over symmetries

!    Computing block of the new symmetry operation according to:
!    !   $ { R1 | v1 }{ R2 | v2 } = { R1.R2 | v1+R1.v2 } $
     matrintoper(:,:) = matmul(bcksymrel(:,:,isym),bcksymrel(:,:,jsym))
     matrinttransl(:) = bcktnons(:,isym)+matmul(bcksymrel(:,:,isym),bcktnons(:,jsym))
     matrintsymafm    = bcksymafm(isym)*bcksymafm(jsym)

!    Rescaling translation between 0 and 1
     do ii=1,3
       if (matrinttransl(ii)>=0.99) then
         do while (matrinttransl(ii)>=0.99)
           matrinttransl(ii)=matrinttransl(ii)-1.0
         end do
       end if
       if (matrinttransl(ii)<0.0) then
         do while (matrinttransl(ii)<0.0)
           matrinttransl(ii)=matrinttransl(ii)+1.0
         end do
       end if
       if ( abs(matrinttransl(ii))<nastyzero) matrinttransl(ii)=0.0
       if ( abs(matrinttransl(ii)-1.0)<nastyzero) matrinttransl(ii)=0.0
     end do

!    Identify the resulting symmetry
     do ii=1,nsym

       flagop=0 ; flagtr=0 ; flagma=0

!      Check for rotation similarity
       if(sum((matrintoper-bcksymrel(:,:,ii))**2)==0)flagop=1

!      Check for translation similarity
       if(maxval((matrinttransl-bcktnons(:,ii))**2)<nastyzero**2)flagtr=1

!      Check for the ferromagnetic character
       if(matrintsymafm==symafm(ii)) flagma=1

       if(flagop+flagtr+flagma==3) then
         symequiv(isym,jsym)=ii
         exit
       end if

     end do

   end do
 end do

 write(6,*) ' Space group multiplication table'
 do isym=1,nsym
   write(6,*) ' Combined operations for the symmetry operation number: ',isym
   do ii=1,14
     if (nsym>ii*16) then
       write(message, '(1x,16i5)' )symequiv(isym,(ii-1)*16+1:ii*16)
       call abi_wrtout(std_out,message,'COLL')
     else
       if (nsym-(ii-1)*16 == 1)  write(message, '(1x,1i5)' )  symequiv(isym,(ii-1)*16+1)
       if (nsym-(ii-1)*16 == 2)  write(message, '(1x,2i5)' )  symequiv(isym,(ii-1)*16+1:(ii-1)*16+2)
       if (nsym-(ii-1)*16 == 3)  write(message, '(1x,3i5)' )  symequiv(isym,(ii-1)*16+1:(ii-1)*16+3)
       if (nsym-(ii-1)*16 == 4)  write(message, '(1x,4i5)' )  symequiv(isym,(ii-1)*16+1:(ii-1)*16+4)
       if (nsym-(ii-1)*16 == 6)  write(message, '(1x,6i5)' )  symequiv(isym,(ii-1)*16+1:(ii-1)*16+6)
       if (nsym-(ii-1)*16 == 8)  write(message, '(1x,8i5)' )  symequiv(isym,(ii-1)*16+1:(ii-1)*16+8)
       if (nsym-(ii-1)*16 == 12) write(message, '(1x,12i5)' ) symequiv(isym,(ii-1)*16+1:(ii-1)*16+12)
       call abi_wrtout(std_out,message,'COLL')
       exit
     end if
   end do
 end do

end subroutine symmultsg
!!***
