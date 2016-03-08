!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_holocell
!! NAME
!! abi_holocell
!!
!! FUNCTION
!! Examine whether the trial conventional cell described by cell_base
!! is coherent with the required holohedral group.
!! Note : for iholohedry=4, the tetragonal axis is not required to be
!! along the C axis.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2010 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cell_base(3,3)=basis vectors of the conventional cell
!!  iholohedry=required holohegral group
!!  iholohedry=1   triclinic      1bar
!!  iholohedry=2   monoclinic     2/m
!!  iholohedry=3   orthorhombic   mmm
!!  iholohedry=4   tetragonal     4/mmm
!!  iholohedry=5   trigonal       3bar m
!!  iholohedry=6   hexagonal      6/mmm
!!  iholohedry=7   cubic          m3bar m
!!
!! OUTPUT
!!  foundc=1 if the basis vectors supports the required holohedry ; =0 otherwise
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine abi_holocell(cell_base,foundc,iholohedry)

 use abi_defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iholohedry
 integer,intent(out) :: foundc
!arrays
 real(dp),intent(in) :: cell_base(3,3)

!Local variables ------------------------------
!scalars
 integer :: allequal,ii,orth
!arrays
 integer :: ang90(3),equal(3)
 real(dp) :: metric(3,3)

!**************************************************************************

 do ii=1,3
   metric(:,ii)=cell_base(1,:)*cell_base(1,ii)+&
&   cell_base(2,:)*cell_base(2,ii)+&
&   cell_base(3,:)*cell_base(3,ii)
 end do

!Examine the angles and vector lengths
 ang90(:)=0
 if(abs(metric(1,2))<tol8)ang90(3)=1
 if(abs(metric(1,3))<tol8)ang90(2)=1
 if(abs(metric(2,3))<tol8)ang90(1)=1
 orth=0
 if(ang90(1)==1 .and. ang90(2)==1 .and. ang90(3)==1) orth=1
 equal(:)=0
 if(abs(metric(1,1)-metric(2,2))<tol8)equal(3)=1
 if(abs(metric(1,1)-metric(3,3))<tol8)equal(2)=1
 if(abs(metric(2,2)-metric(3,3))<tol8)equal(1)=1
 allequal=0
 if(equal(1)==1 .and. equal(2)==1 .and. equal(3)==1) allequal=1

 foundc=0
 if(iholohedry==1)                                      foundc=1
 if(iholohedry==2 .and. ang90(1)+ang90(3)==2 )          foundc=1
 if(iholohedry==3 .and. orth==1)                        foundc=1
 if(iholohedry==4 .and. orth==1 .and.          &
& (equal(3)==1 .or. equal(2)==1 .or. equal(1)==1) ) foundc=1
 if(iholohedry==5 .and. allequal==1 .and. &
& (abs(metric(1,2)-metric(2,3))<tol8) .and. &
& (abs(metric(1,2)-metric(1,3))<tol8)         )      foundc=1
 if(iholohedry==6 .and. equal(3)==1 .and. &
& ang90(1)==1 .and. ang90(2)==1 .and. &
& (2*metric(1,2)-metric(1,1))<tol8           )      foundc=1
 if(iholohedry==7 .and. orth==1 .and. allequal==1)      foundc=1

end subroutine abi_holocell
!!***
