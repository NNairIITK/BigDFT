!{\src2tex{textfont=tt}}
!!****f* ABINIT/xfpack
!! NAME
!! xfpack
!!
!! FUNCTION
!! If option=1, transfer xred, acell, and rprim to vin
!! If option=2, transfer vin  to xred, acell and rprim
!! If option=3, transfer fred and strten to vout
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell0(3)=reference length scales of primitive translations (bohr), needed
!!   for some values of optcell.
!! natom=number of atoms in cell
!! ndim=dimension of vin and vout arrays
!! nsym=order of group.
!! rprimd0(3,3)=reference real space primitive translations,
!!   needed for some values of optcell.
!! optcell=option for the optimisation of the unit cell. Described in abinit_help.
!!  Depending on its value, different part of strten, or acell and rprim
!!  are contained in vin and vout.
!! option= see above
!! strtarget(6)=target stresses ; they will be subtracted from strten when vout
!!  is computed.
!! symrel(3,3,nsym)=symmetry operators in terms of action on primitive translations
!! ucvol=unit cell volume (bohr^3), needed for some values of optcell.
!! ucvol0=reference unit cell volume (bohr^3), needed for some values of optcell.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output variables
!! acell(3)=length scales of primitive translations (bohr)
!! fred(3,natom)=grads of Etot wrt reduced coordinates (hartree)
!! rprim(3,3)=dimensionless real space primitive translations
!! strten(6)=components of the stress tensor (hartree/bohr^3)
!! vin(ndim)=vector that contains xred and some quantity derived
!!   from acell and rprim, depending on the value of optcell.
!! vout(ndim)=vector that contains fred and some quantity derived from
!!   strten, depending on the value of optcell, and taking care ot strtarget
!! xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! PARENTS
!!      brdmin,delocint,moldyn
!!
!! CHILDREN
!!      leave_new,matr3inv,metric,mkrdim,strainsym,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

subroutine xfpack(acell,acell0,fred,natom,ndim,nsym,optcell,option,rprim,rprimd0,&
& strtarget,strten,symrel,ucvol,ucvol0,vin,vout,xred)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ndim,nsym,optcell,option
 real(dp),intent(in) :: ucvol0
 real(dp),intent(inout) :: ucvol
!arrays
 integer,intent(in) :: symrel(3,3)
 real(dp),intent(in) :: acell0(3),rprimd0(3,3),strtarget(6)
 real(dp),intent(inout) :: acell(3),fred(3,natom),rprim(3,3),strten(6)
 real(dp),intent(inout) :: vin(ndim),vout(ndim),xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,kk
 real(dp) :: scale,strdiag
 character(len=500) :: message
!arrays
 real(dp) :: dstr(6),gmet(3,3),gprimd(3,3),gprimd0(3,3),rmet(3,3),rprimd(3,3)
 real(dp) :: rprimd_symm(3,3),scaling(3,3)

! *************************************************************************

 if(optcell==0 .and. ndim/=3*natom)then
   write(message,'(a,a,a,a,a,a,i4,a,i4,a)' )ch10,&
&   ' xfpack: BUG -',ch10,&
&   '  When optcell=0, ndim MUST be equal to 3*natom,',ch10,&
&   '  while ndim=',ndim,' and 3*natom=',3*natom,'.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if( (optcell==1 .or. optcell==4 .or. optcell==5 .or. optcell==6) &
& .and. ndim/=3*natom+1)then
   write(message,'(a,a,a,a,a,a,i4,a,i4,a)' )ch10,&
&   ' xfpack: BUG -',ch10,&
&   '  When optcell=1,4,5 or 6, ndim MUST be equal to 3*natom+1,',ch10,&
&   '  while ndim=',ndim,' and 3*natom+1=',3*natom+1,'.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if( (optcell==2 .or. optcell==3) &
& .and. ndim/=3*natom+6)then
   write(message,'(a,a,a,a,a,a,i4,a,i4,a)' )ch10,&
&   ' xfpack: BUG -',ch10,&
&   '  When optcell=2 or 3, ndim MUST be equal to 3*natom+6,',ch10,&
&   '  while ndim=',ndim,' and 3*natom+6=',3*natom+6,'.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if( optcell>=7 .and. ndim/=3*natom+3)then
   write(message,'(a,a,a,a,a,a,i4,a,i4,a)' )ch10,&
&   ' xfpack: BUG -',ch10,&
&   '  When optcell=7,8 or 9, ndim MUST be equal to 3*natom+3,',ch10,&
&   '  while ndim=',ndim,' and 3*natom+3=',3*natom+3,'.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if(option==1)then

!  Get vin from xred, acell, and rprim
   vin(1:3*natom)= reshape(xred(:,:), (/3*natom/) )
   if(optcell/=0)then
     call mkrdim(acell,rprim,rprimd)
     call strainsym(nsym,rprimd0,rprimd,rprimd_symm,symrel)
     call metric(gmet,gprimd,-1,rmet,rprimd_symm,ucvol)

     if(optcell==1)then

!      vin(3*natom+1)=ucvol**third
       vin(3*natom+1)=(ucvol/ucvol0)**third

     else if(optcell==2 .or. optcell==3 .or. optcell>=7)then

!      Generates gprimd0
       call matr3inv(rprimd0,gprimd0)
       do ii=1,3
         do jj=1,3
           scaling(ii,jj)=0.0_dp
           do kk=1,3
             scaling(ii,jj)=scaling(ii,jj)+rprimd_symm(ii,kk)*gprimd0(jj,kk)
           end do
         end do
       end do
!      Rescale if the volume must be preserved
       if(optcell==3)then
         scale=(ucvol0/ucvol)**third
         scaling(:,:)=scale*scaling(:,:)
       end if
       if(optcell==2 .or. optcell==3)then
         vin(3*natom+1)=scaling(1,1) ; vin(3*natom+4)=(scaling(2,3)+scaling(3,2))*0.5_dp
         vin(3*natom+2)=scaling(2,2) ; vin(3*natom+5)=(scaling(1,3)+scaling(3,1))*0.5_dp
         vin(3*natom+3)=scaling(3,3) ; vin(3*natom+6)=(scaling(1,2)+scaling(2,1))*0.5_dp
       else if(optcell>=7)then
         vin(3*natom+1)=scaling(1,1)
         vin(3*natom+2)=scaling(2,2)
         vin(3*natom+3)=scaling(3,3)
         if(optcell==7)vin(3*natom+1)=(scaling(2,3)+scaling(3,2))*0.5_dp
         if(optcell==8)vin(3*natom+2)=(scaling(1,3)+scaling(3,1))*0.5_dp
         if(optcell==9)vin(3*natom+3)=(scaling(1,2)+scaling(2,1))*0.5_dp
       end if

     else if(optcell==4 .or. optcell==5 .or. optcell==6)then

       vin(3*natom+1)=acell(optcell-3)/acell0(optcell-3)

     end if

   end if

 else if(option==2)then

!  Get xred, and eventually acell and rprim from vin
   xred(:,:)=reshape( vin(1:3*natom), (/3,natom/) )

   if(optcell==1)then

!    acell(:)=acell0(:)*vin(3*natom+1)/(ucvol0**third)
     acell(:)=acell0(:)*vin(3*natom+1)

   else if(optcell==2 .or. optcell==3 .or. optcell>=7 )then

     scaling(:,:)=0.0_dp
     scaling(1,1)=1.0_dp ; scaling(2,2)=1.0_dp ; scaling(3,3)=1.0_dp

     if(optcell==2 .or. optcell==3)then
       scaling(1,1)=vin(3*natom+1)
       scaling(2,2)=vin(3*natom+2)
       scaling(3,3)=vin(3*natom+3)
       scaling(2,3)=vin(3*natom+4) ; scaling(3,2)=vin(3*natom+4)
       scaling(1,3)=vin(3*natom+5) ; scaling(3,1)=vin(3*natom+5)
       scaling(1,2)=vin(3*natom+6) ; scaling(2,1)=vin(3*natom+6)
     else if(optcell==7)then
       scaling(2,2)=vin(3*natom+2) ; scaling(3,3)=vin(3*natom+3)
       scaling(2,3)=vin(3*natom+1) ; scaling(3,2)=vin(3*natom+1)
     else if(optcell==8)then
       scaling(1,1)=vin(3*natom+1) ; scaling(3,3)=vin(3*natom+3)
       scaling(1,3)=vin(3*natom+2) ; scaling(3,1)=vin(3*natom+2)
     else if(optcell==9)then
       scaling(1,1)=vin(3*natom+1) ; scaling(2,2)=vin(3*natom+2)
       scaling(1,2)=vin(3*natom+3) ; scaling(2,1)=vin(3*natom+3)
     end if
     do ii=1,3
       do jj=1,3
         rprimd(ii,jj)=0.0_dp
         do kk=1,3
           rprimd(ii,jj)=rprimd(ii,jj)+scaling(ii,kk)*rprimd0(kk,jj)
         end do
       end do
     end do
!    Rescale if the volume must be preserved
     if(optcell==3)then
       call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
       scale=(ucvol0/ucvol)**third
       rprimd(:,:)=scale*rprimd(:,:)
     end if
     call strainsym(nsym,rprimd0,rprimd,rprimd_symm,symrel)
!    Use a representation based on normalised rprim vectors
     do ii=1,3
       acell(ii)=sqrt(rprimd_symm(1,ii)**2+rprimd_symm(2,ii)**2+rprimd_symm(3,ii)**2)
       rprim(:,ii)=rprimd_symm(:,ii)/acell(ii)
     end do

   else if(optcell==4 .or. optcell==5 .or. optcell==6)then

     acell(:)=acell0(:) ; acell(optcell-3)=vin(3*natom+1)*acell0(optcell-3)

   end if

 else if(option==3)then

!  Get vout from fred and strten
   vout(1:3*natom)= reshape(fred(:,:), (/3*natom/) )
   dstr(:)=strten(:)-strtarget(:)

   if(optcell==1)then

     vout(3*natom+1)=( dstr(1)+dstr(2)+dstr(3))*ucvol

   else if(optcell==2 .or. optcell==3 .or. optcell>=7)then

!    Eventually take away the trace
     strdiag=0.0_dp
     if(optcell==3) strdiag=(dstr(1)+dstr(2)+dstr(3))/3.0_dp
     if(optcell==2 .or. optcell==3)then
       vout(3*natom+1:3*natom+3)=(dstr(1:3)-strdiag)*ucvol
!      For non-diagonal derivatives, must take into account
!      that eps(i,j) AND eps(j,i) are varied at the same time. Thus, derivative
!      is twice larger
       vout(3*natom+4:3*natom+6)=dstr(4:6)*ucvol*2.0_dp
     else if(optcell==7 .or. optcell==8 .or. optcell==9)then
!      Similar to case optcell==2 or optcell==3, but in 2 dimensions.
       vout(3*natom+1:3*natom+3)=dstr(1:3)*ucvol
       vout(3*natom+optcell-6)  =dstr(optcell-3)*ucvol*2.0_dp
     end if

   else if(optcell==4 .or. optcell==5 .or. optcell==6)then

     vout(3*natom+1)=dstr(optcell-3)*ucvol

   end if

 else

   write(message, '(a,a,a,a,a,a,i3,a)' )ch10,&
&   ' xfpack : BUG -',ch10,&
&   '  The only allowed values for option are 1, 2 and 3,',ch10,&
&   '  while it is found that option=',option,'.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')

 end if

end subroutine xfpack
!!***
