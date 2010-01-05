!{\src2tex{textfont=tt}}
!!****f* ABINIT/symanal
!! NAME
!! symanal
!!
!! FUNCTION
!! Derive the name of the point group, from symrel.
!!
!! Also,indicate whether the input value
!! of the holohedry (iholohedry is contained in bravais)
!! is consistent with
!! the point group, as iholohedry was determined directly
!! from the lattice vectors, NOT taking into account
!! that the atomic positions might have broken the symmetry
!! of the lattice.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2009 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nsym=actual number of symmetries
!! symrel(3,3,nsym)= nsym symmetry operations in real space in terms
!! of primitive translations
!!
!! OUTPUT
!! problem= if 1, the holohedry here determined has a lower symmetry
!!  than the one expected from the Bravais lattice number. Otherwise 0.
!!
!! SIDE EFFECTS
!! Input/Output
!! bravais(11): bravais(1)=iholohedry
!!              bravais(2)=center
!!              bravais(3:11)=coordinates of rprimd in the axes
!!              of the conventional bravais lattice (*2 if center/=0)
!! character(len=5) ptgroup=symmetry point group
!!
!! PARENTS
!!      ingeo
!!
!! CHILDREN
!!      leave_new,symdet,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

subroutine symanal(bravais,nsym,problem,ptgroup,symrel)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
! use interfaces_14_hidewrite
! use interfaces_16_hideleave
! use interfaces_42_geometry, except_this_one => symanal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(out) :: problem
 character(len=5),intent(inout) :: ptgroup
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 integer,intent(inout) :: bravais(11)

!Local variables-------------------------------
!scalars
 integer :: iholohedry,inversion,iorder,isym
 character(len=500) :: message
!arrays
 integer :: identity(3,3),matrix(3,3),n_axes(-6:6),trial(3,3)
 integer,allocatable :: determinant(:),order(:),root_invers(:)
 character(len=2),allocatable :: ptsym(:)

!**************************************************************************

!DEBUG
!write(6,*)' symanal : enter'
!do isym=1,nsym
!write(6, '(i3,2x,9i3)' )isym,symrel(:,:,isym)
!end do
!ENDDEBUG

 identity(:,:)=0
 identity(1,1)=1 ; identity(2,2)=1 ; identity(3,3)=1
 n_axes(:)=0

 allocate(determinant(nsym),order(nsym),ptsym(nsym),root_invers(nsym))

!Get the determinant
 call symdet(determinant,nsym,symrel)

!Get the order of each the symmetry operation, as well as the maximal order
!Also, examine whether each symmetry operation is the inversion, or a root
!of the inversion (like -3)
!Finally, decide which kind of point symmetry operation it is
 do isym=1,nsym

  trial(:,:)=identity(:,:)
  matrix(:,:)=symrel(:,:,isym)
  order(isym)=0
  root_invers(isym)=0
  do iorder=1,6
   trial=matmul(matrix,trial)
   if(sum((trial-identity)**2)==0)then
    order(isym)=iorder
    exit
   end if
   if(sum((trial+identity)**2)==0)then
    root_invers(isym)=iorder
    if(iorder==1)inversion=isym
   end if
  end do
  if(order(isym)==0)then
   write(message, '(a,a,a,a,i4,a)' ) ch10,&
&   ' symanal : BUG -',ch10,&
&   '  The symmetry operation number',isym,' is not a root of unity'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
  end if

! determinant, order and root_invers are enough to determine the
! kind of symmetry operation
  ptsym(isym)='no'
  select case(order(isym))
   case(1)
    ptsym(isym)=' 1' ; n_axes(1)=n_axes(1)+1
   case(2)
    if(determinant(isym)== 1)then
     ptsym(isym)=' 2' ; n_axes(2)=n_axes(2)+1
    else if(determinant(isym)==-1 .and. root_invers(isym)==1)then
     ptsym(isym)='-1' ; n_axes(-1)=n_axes(-1)+1
    else if(determinant(isym)==-1 .and. root_invers(isym)==0)then
     ptsym(isym)='-2' ; n_axes(-2)=n_axes(-2)+1
    end if
   case(3)
    ptsym(isym)=' 3' ; n_axes(3)=n_axes(3)+1
   case(4)
    if(determinant(isym)== 1)then
     ptsym(isym)=' 4' ; n_axes(4)=n_axes(4)+1
    else if(determinant(isym)==-1)then
     ptsym(isym)='-4' ; n_axes(-4)=n_axes(-4)+1
    end if
   case(6)
    if(determinant(isym)== 1)then
     ptsym(isym)=' 6' ; n_axes(6)=n_axes(6)+1
    else if(determinant(isym)==-1 .and. root_invers(isym)==3)then
     ptsym(isym)='-3' ; n_axes(-3)=n_axes(-3)+1
    else if(determinant(isym)==-1 .and. root_invers(isym)==0)then
     ptsym(isym)='-6' ; n_axes(-6)=n_axes(-6)+1
    end if
  end select

  if(ptsym(isym)=='no')then
   write(message,'(a,a,a,a,i4,a,a,a,i4,a,a,i4,a,a,i4)' ) ch10,&
&   ' symanal : BUG -',ch10,&
&   '  The symmetry operation number',isym,' could not be identified',ch10,&
&   '  order(isym)      =',order(isym),ch10,&
&   '  determinant(isym)=',determinant(isym),ch10,&
&   '  root_invers(isym)=',root_invers(isym)

   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
  end if

 end do

 iholohedry=0
 if     (sum((n_axes-(/0,0,0,0,0,0, 0 ,1,0,0,0,0,0/))**2)==0)then
  ptgroup='    1' ; iholohedry=1
 else if(sum((n_axes-(/0,0,0,0,0,1, 0 ,1,0,0,0,0,0/))**2)==0)then
  ptgroup='   -1' ; iholohedry=1

 else if(sum((n_axes-(/0,0,0,0,0,0, 0 ,1,1,0,0,0,0/))**2)==0)then
  ptgroup='    2' ; iholohedry=2
 else if(sum((n_axes-(/0,0,0,0,1,0, 0 ,1,0,0,0,0,0/))**2)==0)then
  ptgroup='   -2' ; iholohedry=2
 else if(sum((n_axes-(/0,0,0,0,1,1, 0 ,1,1,0,0,0,0/))**2)==0)then
  ptgroup='  2/m' ; iholohedry=2

 else if(sum((n_axes-(/0,0,0,0,0,0, 0 ,1,3,0,0,0,0/))**2)==0)then
  ptgroup='  222' ; iholohedry=3
 else if(sum((n_axes-(/0,0,0,0,2,0, 0 ,1,1,0,0,0,0/))**2)==0)then
  ptgroup='  mm2' ; iholohedry=3
 else if(sum((n_axes-(/0,0,0,0,3,1, 0 ,1,3,0,0,0,0/))**2)==0)then
  ptgroup='  mmm' ; iholohedry=3

 else if(sum((n_axes-(/0,0,0,0,0,0, 0 ,1,1,0,2,0,0/))**2)==0)then
  ptgroup='    4' ; iholohedry=4
 else if(sum((n_axes-(/0,0,2,0,0,0, 0 ,1,1,0,0,0,0/))**2)==0)then
  ptgroup='   -4' ; iholohedry=4
 else if(sum((n_axes-(/0,0,2,0,1,1, 0 ,1,1,0,2,0,0/))**2)==0)then
  ptgroup='  4/m' ; iholohedry=4
 else if(sum((n_axes-(/0,0,0,0,0,0, 0 ,1,5,0,2,0,0/))**2)==0)then
  ptgroup='  422' ; iholohedry=4
 else if(sum((n_axes-(/0,0,0,0,4,0, 0 ,1,1,0,2,0,0/))**2)==0)then
  ptgroup='  4mm' ; iholohedry=4
 else if(sum((n_axes-(/0,0,2,0,2,0, 0 ,1,3,0,0,0,0/))**2)==0)then
  ptgroup=' -42m' ; iholohedry=4
 else if(sum((n_axes-(/0,0,2,0,5,1, 0 ,1,5,0,2,0,0/))**2)==0)then
  ptgroup='4/mmm' ; iholohedry=4

 else if(sum((n_axes-(/0,0,0,0,0,0, 0 ,1,0,2,0,0,0/))**2)==0)then
  ptgroup='    3' ; iholohedry=5
 else if(sum((n_axes-(/0,0,0,2,0,1, 0 ,1,0,2,0,0,0/))**2)==0)then
  ptgroup='   -3' ; iholohedry=5
 else if(sum((n_axes-(/0,0,0,0,0,0, 0 ,1,3,2,0,0,0/))**2)==0)then
  ptgroup='   32' ; iholohedry=5
 else if(sum((n_axes-(/0,0,0,0,3,0, 0 ,1,0,2,0,0,0/))**2)==0)then
  ptgroup='   3m' ; iholohedry=5
 else if(sum((n_axes-(/0,0,0,2,3,1, 0 ,1,3,2,0,0,0/))**2)==0)then
  ptgroup='  -3m' ; iholohedry=5

 else if(sum((n_axes-(/0,0,0,0,0,0, 0 ,1,1,2,0,0,2/))**2)==0)then
  ptgroup='    6' ; iholohedry=6
 else if(sum((n_axes-(/2,0,0,0,1,0, 0 ,1,0,2,0,0,0/))**2)==0)then
  ptgroup='   -6' ; iholohedry=6
 else if(sum((n_axes-(/2,0,0,2,1,1, 0 ,1,1,2,0,0,2/))**2)==0)then
  ptgroup='  6/m' ; iholohedry=6
 else if(sum((n_axes-(/0,0,0,0,0,0, 0 ,1,7,2,0,0,2/))**2)==0)then
  ptgroup='  622' ; iholohedry=6
 else if(sum((n_axes-(/0,0,0,0,6,0, 0 ,1,1,2,0,0,2/))**2)==0)then
  ptgroup='  6mm' ; iholohedry=6
 else if(sum((n_axes-(/2,0,0,0,4,0, 0 ,1,3,2,0,0,0/))**2)==0)then
  ptgroup=' -62m' ; iholohedry=6
 else if(sum((n_axes-(/2,0,0,2,7,1, 0 ,1,7,2,0,0,2/))**2)==0)then
  ptgroup='6/mmm' ; iholohedry=6

 else if(sum((n_axes-(/0,0,0,0,0,0, 0 ,1,3,8,0,0,0/))**2)==0)then
  ptgroup='   23' ; iholohedry=7
 else if(sum((n_axes-(/0,0,0,8,3,1, 0 ,1,3,8,0,0,0/))**2)==0)then
  ptgroup='  m-3' ; iholohedry=7
 else if(sum((n_axes-(/0,0,0,0,0,0, 0 ,1,9,8,6,0,0/))**2)==0)then
  ptgroup='  432' ; iholohedry=7
 else if(sum((n_axes-(/0,0,6,0,6,0, 0 ,1,3,8,0,0,0/))**2)==0)then
  ptgroup=' -43m' ; iholohedry=7
 else if(sum((n_axes-(/0,0,6,8,9,1, 0 ,1,9,8,6,0,0/))**2)==0)then
  ptgroup=' m-3m' ; iholohedry=7

 end if

 if(iholohedry==0)then
  write(message, '(a,a,a,a)' )ch10,&
&  ' symanal : BUG -',ch10,&
&  '  Could not find the point group'
  call wrtout(std_out,message,'COLL')
  call leave_new('COLL')
 end if

!DEBUG
!write(6, '(a,13i3)' )' symanal : n_axes(-6:6)=',n_axes(-6:6)
!write(6,*)' iholohedry, ptgroup=',iholohedry,',',ptgroup
!ENDDEBUG

!Examine the agreement with bravais(1)
!Warning : might change Bravais lattice hR to hP, if hexagonal axes
 problem=0
 select case (bravais(1))
  case(7)
   if(iholohedry<6)problem=1
   if(iholohedry==6)problem=2
  case(6)
   if(iholohedry<4)problem=1
   if(iholohedry==7 .or. iholohedry==4)problem=2
!  Here, change hR into hP
   if(iholohedry==5)iholohedry=6
  case(5)
   if(iholohedry<4)problem=1
   if(iholohedry==7 .or. iholohedry==6 .or. iholohedry==4)problem=2
  case(4)
   if(iholohedry<4)problem=1
   if(iholohedry>4)problem=2
  case(3)
   if(iholohedry<3)problem=1
   if(iholohedry>3)problem=2
  case(2)
   if(iholohedry<2)problem=1
   if(iholohedry>2)problem=2
  case(1)
   if(iholohedry>1)problem=2
 end select

 if(problem==1)then
  write(message, '(a,a,a,a,a,a,i3,a,a,a,i3,a,a,a)' )ch10,&
&  ' symanal : COMMENT -',ch10,&
&  '  The Bravais lattice determined only from the primitive',ch10,&
&  '  vectors, bravais(1)=',bravais(1),', is more symmetric',ch10,&
&  '  than the real one, iholohedry=',iholohedry,', obtained by taking into',ch10,&
&  '  account the atomic positions.'
  call wrtout(std_out,message,'COLL')
 end if

 if(problem==2)then
  write(message, '(6a,i3,3a,i3,7a)' )ch10,&
&  ' symanal : BUG -',ch10,&
&  '  The Bravais lattice determined only from the primitive',ch10,&
&  '  vectors (rprim or angdeg), bravais(1)=',bravais(1),', is not compatible',ch10,&
&  '  with the real one, iholohedry=',iholohedry,', obtained by taking into',ch10,&
&  '  account the atomic positions. This might be due to an insufficient',ch10,&
&  '  number of digits in the specification of rprim (at least 10),',ch10,&
&  '  or to an erroneous rprim or angdeg. If this is not the case, then ...'
  call wrtout(std_out,message,'COLL')
  call leave_new('COLL')
 end if

 bravais(1)=iholohedry

!DEBUG
!do isym=1,nsym
!write(6, '(a,3i5)' )&
!&  ' symanal : isym,determinant,order=',isym,determinant(isym),order(isym)
!end do
!ENDDEBUG

 deallocate(determinant,order,ptsym,root_invers)

end subroutine symanal
!!***
