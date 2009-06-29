!{\src2tex{textfont=tt}}
!!****f* ABINIT/symlist_fcc
!! NAME
!! symlist_fcc
!!
!! FUNCTION
!! Determine the space group from the number and type of symmetry operations
!! FCC case
!!
!! COPYRIGHT
!! Copyright (C) 2000-2009 ABINIT group (RC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! additional_info=information that is needed beyond n_axes, in order
!!  to discriminate between specific space groups
!! brvltt=Bravais lattice type
!! nsym=actual number of symmetries
!! n_axes(31)=array containing the number of all the possible symmetry operations
!!
!! OUTPUT
!! spgroup=space group number ; returns 0 if not found
!!
!! NOTES
!!
!! The list of symmetry operations is for the conventional cell
!!
!! TODO
!! For the time being there are several groups where uncertainties still exist
!! This will be solved in the very next ABINIT version
!!
!! PARENTS
!!      symspgr
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine symlist_fcc(additional_info,brvltt,nsym,n_axes,spgroup)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: additional_info,brvltt,nsym
 integer,intent(out) :: spgroup
!arrays
 integer,intent(in) :: n_axes(31)

!Local variables-------------------------------
!character(len=500) :: message
!arrays
 integer :: n_axest(31)

!**************************************************************************

!DEBUG
!write(6,*) ' symlist_fcc : enter '
!write(6,*) ' nsym = ', nsym
!write(6,*) ' brvltt = ',brvltt
!write(6, '(a,10i3)' ) ' n_axes(1:10) =',n_axes(1:10)
!write(6, '(a,10i3)' ) ' n_axes(11:20)=',n_axes(11:20)
!write(6, '(a,11i3)' ) ' n_axes(21:31)=',n_axes(21:31)
!ENDDEBUG

 spgroup=0

 select case(nsym)

  case(16)

   n_axest=(/0,0,0,0,0,0,3,1,6,0,  0,0,0,0,0,0,0,0,0,6,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=22
   n_axest=(/0,0,0,0,0,0,3,1,2,0,  0,0,0,0,2,4,0,2,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=42
   n_axest=(/0,0,0,0,0,0,3,1,2,0,  0,0,0,0,0,0,8,0,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=43


  case(32)

   n_axest=(/0,0,0,0,4,0,3,1,6,0,  0,0,0,0,3,6,0,3,0,6,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=69
   n_axest=(/0,0,0,0,4,0,3,1,6,0,  0,0,0,0,0,0,12,0,0,6, 0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=70

  case(48)

   n_axest=(/0,0,0,0,0,0,3,1,6,32, 0,0,0,0,0,0,0,0,0,6, 0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=196

  case(96)

   n_axest=(/0,0,32,0,4,0,3,1,6,32,  0,0,0,0,3,6,0,3,0,6,   0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=202
   n_axest=(/0,0,32,0,4,0,3,1,6,32,  0,0,0,0,0,0,12,0,0,6,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=203
   n_axest=(/0,0,0,0,0,0,3,1,18,32,  0,12,0,0,0,0,0,0,0,18, 0,0,0,0,12,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=209
   n_axest=(/0,0,0,0,0,0,3,1,18,32,  0,0,0,0,0,0,0,0,0,18,  0,0,0,24,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=210
   n_axest=(/0,24,0,0,0,0,3,1,6,32,  0,0,0,0,9,0,0,0,15,6,   0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=216
   n_axest=(/0,24,0,0,0,0,3,1,6,32,  0,0,0,0,0,9,0,0,15,6,   0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=219
   n_axest=(/0,24,0,0,0,0,3,1,6,32,  0,0,0,0,0,9,0,3,12,6,   0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=219


  case(192)

!  Note that the identification of the mirror planes is still ambiguous for cF
   n_axest=(/0,24,32,0,4,0,3,1,18,32,  0,12,0,0,12,6,0,3,15,18, 0,0,0,0,12,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=225
   n_axest=(/0,24,32,0,4,0,3,1,18,32,  0,12,0,0,3,15,0,6,12,18,  0,0,0,0,12,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=226
   n_axest=(/0,24,32,0,4,0,3,1,18,32,  0,0,0,0,9,0,12,0,15,18,  0,0,0,24,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=227
   n_axest=(/0,24,32,0,4,0,3,1,18,32,  0,0,0,0,6,0,12,0,18,18,  0,0,0,24,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=227
   n_axest=(/0,24,32,0,4,0,3,1,18,32,  0,0,0,0,0,6,12,6,12,18,  0,0,0,24,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=228
   n_axest=(/0,24,32,0,4,0,3,1,18,32,  0,0,0,0,0,9,12,3,12,18,  0,0,0,24,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=228

 end select

end subroutine symlist_fcc
!!***
