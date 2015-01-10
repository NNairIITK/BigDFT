!{\src2tex{textfont=tt}}
!!****f* ABINIT/symspgr
!! NAME
!! symspgr
!!
!! FUNCTION
!! Using the type of each symmetry operation
!! (found in symplanes.f and symaxes.f):
!! proper symmetries 1,2,2_1,3,3_1,3_2,4,4_1,4_2,4_3,6,6_1,...6_5
!! improper symmetries -1,m,a,b,c,d,n,g,-3,-4,-6 ,
!! build an array with the number of such operations.
!! then, call symlist.f to identify the space group.
!! The identification is not unambiguous still ...
!!
!! COPYRIGHT
!! Copyright (C) 2000-2010 ABINIT group (RC, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! bravais(11): bravais(1)=iholohedry
!!              bravais(2)=center
!!              bravais(3:11)=coordinates of rprimd in the axes
!!              of the conventional bravais lattice (*2 if center/=0)
!! nsym=actual number of symmetries
!! symrel(3,3,nsym)= nsym symmetry operations in real space in terms
!!   of primitive translations
!! tnons(3,nsym)=nonsymmorphic translations for each symmetry (would
!!   be 0 0 0 each for a symmorphic space group)
!!
!! OUTPUT
!! spgroup=symmetry space group number
!!
!! NOTES
!! It is assumed that the symmetry operations will be entered in the
!! symrel tnons arrays, for the PRIMITIVE cell. The matrix of transformation
!! from the primitive cell to the conventional cell is described
!! in the array "bravais" (see symlatt.F90).
!! The present routine first make the transformation from the
!! primitive coordinates to the conventional ones, then eventually
!! generate additional symmetries, taking into account the
!! centering translations.
!! Then, the order and determinant of each symmetry operation
!! is determined.
!!
!! For proper symmetries (rotations), the
!! associated translation is also determined.
!! However, left or right handed screw rotations are
!! not (presently) distinguished, and will be attributed equally
!! to left or right.
!!
!! For the detailed description of the labelling of the axes,
!! see symaxes.f and symplanes.f
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine symspgr(bravais,nsym,spgroup,symrel,tnons,tolsym)

 use defs_basis
 use abi_interfaces_lowlevel
 use interfaces_42_geometry, except_this_one => symspgr

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(out) :: spgroup
 real(dp),intent(in) :: tolsym
!arrays
 integer,intent(in) :: bravais(11),symrel(3,3,nsym)
 real(dp),intent(inout) :: tnons(3,nsym)

!Local variables-------------------------------
!scalars
 integer :: additional_info,brvltt,center,direction=0,identified,iholohedry,ii
 integer :: iorder,ishift,isym,nshift,nsymconv,option,spgaxor,spgorig,sporder
 integer :: tnons_order,type_axis
 character(len=1) :: brvsb
 character(len=15) :: intsb,ptintsb,ptschsb,schsb
 character(len=35) :: intsbl
 character(len=500) :: message
!arrays
 integer :: identity(3,3),matrix(3,3),n_axes(31),n_axest(31),test_direction(3)
 integer :: trial(3,3)
 integer,allocatable :: determinant(:),order(:),symrelconv(:,:,:),t_axes(:)
 real(dp) :: axes(3,3),reduced(3),rprimdconv(3,3),trialt(3),vect(3,3)
 real(dp),allocatable :: shift(:,:),tnonsconv(:,:)

!**************************************************************************

!Initialize brvltt, from bravais(2) and bravais(1)
 center=bravais(2)
 iholohedry=bravais(1)
 brvltt=1
 if(center==-1)brvltt=2  ! Inner centering
 if(center==-3)brvltt=3  ! Face centering
 if(center==1)brvltt=5  ! A-Face centering
 if(center==2)brvltt=6  ! B-Face centering
 if(center==3)brvltt=4  ! C-Face centering
 if(iholohedry==5)brvltt=7  ! Rhombohedral

!Produce the symmetry operations, in the axis of the
!conventional cell
 nsymconv=nsym
 if(center/=0)nsymconv=2*nsymconv
 if(center==-3)nsymconv=4*nsym
 allocate(symrelconv(3,3,nsymconv),tnonsconv(3,nsymconv))

!Produce symrel and tnons in conventional axes,
!name them symrelconv and tnonsconv
 rprimdconv(:,1)=bravais(3:5)
 rprimdconv(:,2)=bravais(6:8)
 rprimdconv(:,3)=bravais(9:11)

 if(center/=0)rprimdconv(:,:)=rprimdconv(:,:)*half

 axes(:,:)=zero
 axes(1,1)=one ; axes(2,2)=one ; axes(3,3)=one
 symrelconv(:,:,1:nsym)=symrel(:,:,1:nsym)
!Note that the number of symmetry operations is still nsym
 call symrelrot(nsym,rprimdconv,axes,symrelconv,tolsym)
 option=1
 call xredxcart(nsym,option,rprimdconv,tnonsconv,tnons)
!Gives the associated translation, with components in the
!interval ]-0.5,0.5] .
 tnonsconv(:,1:nsym)=tnonsconv(:,1:nsym)-nint(tnonsconv(:,1:nsym)-tol6)

!If the Bravais lattice is centered, duplicate or quadruplicate
!the number of symmetry operations, using the Bravais
!lattice shifts
 nshift=1
 if(center/=0)nshift=2
 if(center==-3)nshift=4
 allocate(shift(3,nshift))
 shift(:,1)=zero
 if(center/=0 .and. center/=-3)then
   shift(:,2)=half
   if(center==1)shift(1,2)=zero
   if(center==2)shift(2,2)=zero
   if(center==3)shift(3,2)=zero
 else if(center==-3)then
   shift(:,2)=half ; shift(1,2)=zero
   shift(:,3)=half ; shift(2,3)=zero
   shift(:,4)=half ; shift(3,4)=zero
 end if ! center/=0 or -3
 if(nshift/=1)then
   do ishift=2,nshift
     symrelconv(:,:,(ishift-1)*nsym+1:ishift*nsym)=symrelconv(:,:,1:nsym)
     do isym=1,nsym
       tnonsconv(:,(ishift-1)*nsym+isym)=tnonsconv(:,isym)+shift(:,ishift)
     end do
   end do ! ishift
 end if ! nshift/=1

!At this stage, all the symmetry operations are available,
!expressed in the conventional axis, and also include
!the Bravais lattive translations, and associated operations...

 identity(:,:)=0
 identity(1,1)=1 ; identity(2,2)=1 ; identity(3,3)=1
 n_axes(:)=0

 allocate(determinant(nsymconv),order(nsymconv))

!Get the determinant
 call symdet(determinant,nsymconv,symrelconv)

!Get the order of each the symmetry operation, as well as the maximal order
!Also, examine whether each symmetry operation is the inversion, or a root
!of the inversion (like -3)
!Decide which kind of point symmetry operation it is
!Finally assign tnonsconv order and decide the space symmetry operation

 allocate(t_axes(nsymconv))

 do isym=1,nsymconv

   trial(:,:)=identity(:,:)
   matrix(:,:)=symrelconv(:,:,isym)

   order(isym)=0
   do iorder=1,6
     trial=matmul(matrix,trial)
     if(sum((trial-identity)**2)==0)then
       order(isym)=iorder
       exit
     end if
     if(sum((trial+identity)**2)==0)then
       order(isym)=iorder
       exit
     end if
   end do

   if(order(isym)==0)then
     write(message, '(a,a,a,a,i4,a)' ) ch10,&
&     ' symspgr : BUG -',ch10,&
&     '  The symmetry operation number',isym,' is not a root of unity'
     call abi_wrtout(std_out,message,'COLL')
     call abi_leave_new('COLL')
   end if

!  Determination of the characteristics of proper symmetries (rotations)
   if(determinant(isym)==1)then

!    Determine the translation vector associated to the rotations
!    and its order : apply the symmetry operation
!    then analyse the resulting vector.
     identified=0
     trialt(:)=zero
     do ii=1,order(isym)
       trialt(:)=matmul(symrelconv(:,:,isym),trialt(:)) &
&       +tnonsconv(:,isym)
     end do
!    Gives the associated translation, with components in the
!    interval [-0.5,0.5] .
     reduced(:)=trialt(:)-nint(trialt(:)-tol6)

     if(sum(abs(reduced(:)))<tol6)identified=1
     if( (center==1 .or. center==-3) .and. &
&     sum(abs(reduced(:)-(/zero,half,half/)))<tol6 )identified=2
     if( (center==2 .or. center==-3) .and. &
&     sum(abs(reduced(:)-(/half,zero,half/)))<tol6 )identified=3
     if( (center==3 .or. center==-3) .and. &
&     sum(abs(reduced(:)-(/half,half,zero/)))<tol6 )identified=4
     if(center==-1.and. sum(abs(reduced(:)-(/half,half,half/)))<tol6 )identified=5

!    If the symmetry operation has not been identified, there is a problem ...
     if(identified==0) then
       write(message, '(a,a,i3,a,3(a,3i4,a),a,3es22.12,a,a,3es22.12)' )&
&       ch10,' symspgr : problem with isym=',isym,ch10,&
&       '  symrelconv(:,1,isym)=',symrelconv(:,1,isym),ch10,&
&       '  symrelconv(:,2,isym)=',symrelconv(:,2,isym),ch10,&
&       '  symrelconv(:,3,isym)=',symrelconv(:,3,isym),ch10,&
&       '  tnonsconv(:,isym)=',tnonsconv(:,isym),ch10,&
&       '  trialt(:)=',trialt(:)
       call abi_wrtout(std_out,message,'COLL')
       write(message, '(4a,i4,2a)' )&
&       ch10,&
&       ' symspgr : BUG -',ch10,&
&       '  The space symmetry operation number',isym,ch10,&
&       '  is not a (translated) root of unity'
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if

!    Compute the translation vector associated with one rotation
     trialt(:)=trialt(:)/order(isym)
     trialt(:)=trialt(:)-nint(trialt(:)-tol6)

!    Analyse the resulting vector.
     identified=0
     do ii=1,order(isym)
       reduced(:)=ii*trialt(:)-nint(ii*trialt(:)-tol6)
       if(sum(abs(reduced(:)))<tol6)identified=1
       if( (center==1 .or. center==-3) .and. &
&       sum(abs(reduced(:)-(/zero,half,half/)))<tol6 )identified=2
       if( (center==2 .or. center==-3) .and. &
&       sum(abs(reduced(:)-(/half,zero,half/)))<tol6 )identified=3
       if( (center==3 .or. center==-3) .and. &
&       sum(abs(reduced(:)-(/half,half,zero/)))<tol6 )identified=4
       if(center==-1.and. sum(abs(reduced(:)-(/half,half,half/)))<tol6 )identified=5

       if(identified/=0)then
         tnons_order=ii
         exit
       end if
     end do ! ii

!    Determinant (here=+1, as we are dealing with proper symmetry operations),
!    order, tnons_order and identified are enough to
!    determine the kind of symmetry operation

     select case(order(isym))
       case(1)                       ! point symmetry 1
         if(identified==1) then
           type_axis=8                 ! 1
           write(message,'(a,i3,a)') &
&           ' symspgr : the symmetry operation no. ',isym,' is the identity'
           call abi_wrtout(6,message,'COLL')
         else
           type_axis=7                 ! t
           write(message,'(a,i3,a)') &
&           ' symspgr : the symmetry operation no. ',isym,' is a pure translation '
           call abi_wrtout(6,message,'COLL')
         end if
       case(2,3,4,6)                 ! point symmetry 2,3,4,6 - rotations
         call symaxes(center,iholohedry,&
&         isym,symrelconv(:,:,isym),order(isym),tnons_order,trialt,type_axis)
     end select

   else if (determinant(isym)==-1)then

!    Now, take care of the improper symmetry operations.
!    Their treatment is relatively easy, except for the mirror planes
     select case(order(isym))
       case(1)                       ! point symmetry 1
         type_axis=5                  ! -1
         write(message,'(a,i3,a)') &
&         ' symspgr : the symmetry operation no. ',isym,' is an inversion'
         call abi_wrtout(6,message,'COLL')
       case(2)                       ! point symmetry 2 - planes
         call symplanes(center,iholohedry,isym,&
&         symrelconv(:,:,isym),tnonsconv(:,isym),type_axis)
       case(3)                       ! point symmetry 3
         type_axis=3                  ! -3
         write(message,'(a,i3,a)') &
&         ' symspgr : the symmetry operation no. ',isym,' is a -3 axis '
         call abi_wrtout(6,message,'COLL')
       case(4)                       ! point symmetry 1
         type_axis=2                  ! -4
         write(message,'(a,i3,a)') &
&         ' symspgr : the symmetry operation no. ',isym,' is a -4 axis '
         call abi_wrtout(6,message,'COLL')
       case(6)                       ! point symmetry 1
         type_axis=1                  ! -6
         write(message,'(a,i3,a)') &
&         ' symspgr : the symmetry operation no. ',isym,' is a -6 axis '
         call abi_wrtout(6,message,'COLL')
     end select

   end if ! determinant==1 or -1

   t_axes(isym)=type_axis
   n_axes(type_axis)=n_axes(type_axis)+1

 end do ! isym=1,nsymconv

 if (sum(n_axes)-nsymconv/=0) then
   write(message, '(10a)' )ch10,&
   ' symspgr : WARNING -',ch10,&
   '  Not all the symmetries have been recognized. ',ch10,&
   '  This might be due either to an error in the input file',ch10,&
   '  or to a BUG in ABINIT',ch10,&
   '  Please contact the ABINIT group.'
   call abi_wrtout(6,message,'COLL')
 end if

!Treat cases in which the space group cannot be identified on the
!basis of n_axes one need additional informations
 if(brvltt==1)then
!  If the bravais lattice is primitive
   if(nsymconv==4)then
     n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,0,2,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0)then    ! Spgroup 27 (Pcc2) or 32 (Pba2)
       write(6,*)' symspgr: 27 or 32'
       additional_info=2
!      Select binary axis
       do isym=1,nsymconv
         if(t_axes(isym)==8)then
!          Find direction of binary axis
           if(symrelconv(1,1,isym)==1)direction=1
           if(symrelconv(2,2,isym)==1)direction=2
           if(symrelconv(3,3,isym)==1)direction=3
         end if
       end do
!      Examine the projection of the translation vector of the a, b or c mirror planes
!      onto the binary axis
       do isym=1,nsymconv
         if(t_axes(isym)==16)then
           if(abs(tnonsconv(direction,isym))>tol8)additional_info=1
         end if
       end do
     end if
   else if(nsymconv==8)then
     n_axest=(/0,0,0,0,1,0,0,1,1,0,  0,0,0,0,1,2,0,0,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0)then    ! Spgroup 55 (Pbam) or 57 (Pbcm)
       write(6,*)' symspgr: 55 or 57'
       additional_info=1
!      Select mirror plane m
       do isym=1,nsymconv
         if(t_axes(isym)==15)then
!          Find direction of mirror plane
           if(symrelconv(1,1,isym)==-1)direction=1
           if(symrelconv(2,2,isym)==-1)direction=2
           if(symrelconv(3,3,isym)==-1)direction=3
         end if
       end do
!      Examine the projection of the translation vector of the a, b, or c mirror planes
!      onto the binary axis
       do isym=1,nsymconv
         if(t_axes(isym)==16)then
           if(abs(tnonsconv(direction,isym))>tol8)additional_info=2
         end if
       end do
     end if
     n_axest=(/0,0,0,0,1,0,0,1,1,0,  0,0,0,0,0,2,0,1,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0)then    ! Spgroup 56 (Pccn) or 60 (Pbcn)
       write(6,*)' symspgr: 56 or 60'
       additional_info=1
!      Select mirror plane n
       do isym=1,nsymconv
         if(t_axes(isym)==18)then
!          Find direction of mirror plane
           if(symrelconv(1,1,isym)==-1)direction=1
           if(symrelconv(2,2,isym)==-1)direction=2
           if(symrelconv(3,3,isym)==-1)direction=3
         end if
       end do
!      Examine the projection of the translation vector of the a, b, or c mirror planes
!      onto the binary axis
       do isym=1,nsymconv
         if(t_axes(isym)==16)then
           if(abs(tnonsconv(direction,isym))<tol8)additional_info=2
         end if
       end do
     end if
   end if
 else if(brvltt==2)then
!  In the few next lines, use additional_info as a flag
   additional_info=0
!  If the bravais lattice is inner-centered
   if(nsymconv==8)then
!    Test spgroup 23 (I222) or 24 (I2_{1}2_{1}2_{1})
     n_axest=(/0,0,0,0,0,0,1,1,3,0,  0,0,0,0,0,0,0,0,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) additional_info=1
   else if(nsymconv==24)then
!    Test spgroup 197 (I23) or 199 (I2_{1}3)
     n_axest=(/0,0,0,0,0,0,1,1,3,16, 0,0,0,0,0,0,0,0,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) additional_info=1
   end if
   if(additional_info==1)then
     write(6,*)' symspgr: (23 or 24) or (197 or 199)'
!    Select the three binary axes (they might be 2 or 2_1 !)
     test_direction(:)=0
     do isym=1,nsymconv
       if(t_axes(isym)==20)then
!        Find direction of axis
         do direction=1,3
           if(symrelconv(direction,direction,isym)==1)then
             test_direction(direction)=1
             if(abs(tnonsconv(direction,isym))<tol8)then
               vect(:,direction)=tnonsconv(:,isym)
             else
               vect(:,direction)=tnonsconv(:,isym)+half
             end if
             vect(:,direction)=vect(:,direction)-nint(vect(:,direction)-tol8)
             vect(direction,direction)=zero
           end if
         end do ! direction=1,3
       end if ! if binary axis
     end do ! isym
     if(test_direction(1)/=1 .or. test_direction(2)/=1 .and. test_direction(3)/=1)then
       write(message, '(8a,3i4)' )ch10,&
&       ' symspgr : BUG -',ch10,&
&       '  For space groups 23, 24, 197 or 197, the three binary axes',ch10,&
&       '  are not equally partitioned along the x, y and z directions',ch10,&
&       '  test_direction(1:3)=',test_direction(:)
       call abi_wrtout(std_out,message,'COLL')
       call abi_leave_new('COLL')
     end if
     additional_info=1
     if(abs(vect(1,2)-vect(1,3))>tol8 .or. &
&     abs(vect(2,1)-vect(2,3))>tol8 .or. &
&     abs(vect(3,1)-vect(3,2))>tol8       ) additional_info=2
   end if ! additional informations are needed
 end if ! brvltt==1

 if(brvltt==0 .or. brvltt==1)then         ! Primitive
   call symlist_prim(additional_info,nsymconv,n_axes,spgroup)
 else if(brvltt==2)then
   call symlist_bcc(additional_info,nsymconv,n_axes,spgroup)
 else if(brvltt==3)then
   call symlist_fcc(nsymconv,n_axes,spgroup)
 else
   call symlist_others(brvltt,nsymconv,n_axes,spgroup)
 end if

 if(spgroup==0) then
   write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
   ' symspgr : WARNING -',ch10,&
   '  Could not find the space group.',ch10,&
   '  This often happens when the user selects a restricted set of symmetries ',ch10,&
   '  in the input file, instead of letting the code automatically find symmetries.'
   call abi_wrtout(6,message,'COLL')
 end if

 spgorig=1 ; spgaxor=1
 call spgdata(brvsb,intsb,intsbl,ptintsb,ptschsb,&
& schsb,spgaxor,spgroup,sporder,spgorig)

 if(spgroup/=0)then
   write(message, '(a,i4,2x,a,a,a,a,a)' ) ' symspgr : spgroup=',spgroup,&
&   trim(brvsb),trim(intsb),'   (=',trim(schsb),')'
   call abi_wrtout(6,message,'COLL')
 end if

 deallocate(determinant,order)
 deallocate(shift,symrelconv,tnonsconv,t_axes)

end subroutine symspgr
!!***
