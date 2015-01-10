!{\src2tex{textfont=tt}}
!!****f* ABINIT/metric
!! NAME metric
!! metric
!!
!! FUNCTION
!! Compute first dimensional primitive translation vectors in reciprocal space
!!  gprimd from rprimd, and eventually writes out.
!! Then, computes metrics for real and recip space
!!  rmet and gmet using length
!!  dimensional primitive translation vectors
!!  in columns of rprimd(3,3) and gprimd(3,3).
!!  gprimd is the inverse transpose of rprimd.
!!  i.e. $ rmet_{i,j}= \sum_k ( rprimd_{k,i}*rprimd_{k,j} )  $
!!       $ gmet_{i,j}= \sum_k ( gprimd_{k,i}*gprimd_{k,j} )  $
!! Also computes unit cell volume ucvol in $\textrm{bohr}^3$
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  iout=unit number of output file.  If iout<0, do not write output.
!!
!! OUTPUT
!!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  rmet(3,3)=real space metric ($\textrm{bohr}^{2}$).
!!  ucvol=unit cell volume ($\textrm{bohr}^{3}$).
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine metric(gmet,gprimd,iout,rmet,rprimd,ucvol)

 use defs_basis
 use abi_interfaces_lowlevel
 use interfaces_32_util

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout
 real(dp),intent(out) :: ucvol
!arrays
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(out) :: gmet(3,3),gprimd(3,3),rmet(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,nu
 character(len=500) :: message
!arrays
 real(dp) :: angle(3)

! *************************************************************************

!Compute unit cell volume
 ucvol=rprimd(1,1)*(rprimd(2,2)*rprimd(3,3)-rprimd(3,2)*rprimd(2,3))+&
& rprimd(2,1)*(rprimd(3,2)*rprimd(1,3)-rprimd(1,2)*rprimd(3,3))+&
& rprimd(3,1)*(rprimd(1,2)*rprimd(2,3)-rprimd(2,2)*rprimd(1,3))

!Check that the input primitive translations are not
!linearly dependent (and none is zero); i.e. ucvol~=0
!Also ask that the mixed product is positive.
 if (abs(ucvol)<1.0d-12) then 
   write(message, '(8a)' ) ch10,&
&   ' metric : ERROR -',ch10,&
&   '  Input rprim and acell gives vanishing unit cell volume.',ch10,&
&   '  This indicates linear dependency between primitive lattice vectors',&
&   ch10,'  Action : correct either rprim or acell in input file.'
   call abi_wrtout(std_out,message,'COLL')
   call abi_leave_new('COLL')
 end if
 if (ucvol<0.0d0)then
   write(message,&
&   '(5a,3(a,3es16.6,a),7a)' ) ch10,&
&   ' metric : ERROR -',ch10,&
&   '  Current rprimd gives negative (R1xR2).R3 . ',ch10,&
&   '  Rprimd =',rprimd(:,1),ch10,&
&   '          ',rprimd(:,2),ch10,&
&   '          ',rprimd(:,3),ch10,&
&   '  Action : if the cell size and shape are fixed (optcell==0),',ch10,&
&   '   exchange two of the input rprim vectors;',ch10,&
&   '   if you are optimizing the cell size and shape (optcell/=0),',ch10,&
&   '   maybe the move was too large, and you might try to decrease strprecon.'
   call abi_wrtout(std_out,message,'COLL')
   call abi_leave_new('COLL')
 end if

!Generates gprimd
 call matr3inv(rprimd,gprimd)

!Write out rprimd, gprimd and ucvol
 if (iout>=0) then
   write(message, '(a,a)' )' Real(R)+Recip(G) ',&
&   'space primitive vectors, cartesian coordinates (Bohr,Bohr^-1):'
   call abi_wrtout(iout,message,'COLL')
   do nu=1,3
     write(message, '(1x,a,i1,a,3f11.7,2x,a,i1,a,3f11.7)' ) &
&     'R(',nu,')=',rprimd(:,nu)+tol10,&
&     'G(',nu,')=',gprimd(:,nu)+tol10
     call abi_wrtout(iout,message,'COLL')
   end do
   write(message, '(a,1p,e15.7,a)' ) &
&   ' Unit cell volume ucvol=',ucvol+tol10,' bohr^3'
   call abi_wrtout(iout,message,'COLL')
   call abi_wrtout(std_out,  message,'COLL')
 end if

!Compute real space metrics
 do ii=1,3
   rmet(ii,:)=rprimd(1,ii)*rprimd(1,:)+&
&   rprimd(2,ii)*rprimd(2,:)+&
&   rprimd(3,ii)*rprimd(3,:)
 end do

!Compute reciprocal space metrics
 do ii=1,3
   gmet(ii,:)=gprimd(1,ii)*gprimd(1,:)+&
&   gprimd(2,ii)*gprimd(2,:)+&
&   gprimd(3,ii)*gprimd(3,:)
 end do

!Write out the angles
 if (iout>=0) then
   angle(1)=acos(rmet(2,3)/sqrt(rmet(2,2)*rmet(3,3)))/two_pi*360.0d0
   angle(2)=acos(rmet(1,3)/sqrt(rmet(1,1)*rmet(3,3)))/two_pi*360.0d0
   angle(3)=acos(rmet(1,2)/sqrt(rmet(1,1)*rmet(2,2)))/two_pi*360.0d0
   write(message, '(a,3es16.8,a)' )&
&   ' Angles (23,13,12)=',angle(1:3),' degrees'
   call abi_wrtout(iout,message,'COLL')
   call abi_wrtout(std_out,message,'COLL')
 end if

end subroutine metric
!!***
