!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkprimit
!! NAME  chkprimit
!! chkprimit
!!
!!
!! FUNCTION
!! Check whether the cell is primitive or not.
!! If chkprim/=0 and the cell is non-primitive, stops.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! chkprim= if non-zero, check that the unit cell is primitive.
!! nsym=actual number of symmetries
!! symafm(nsym)= (anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym)= nsym symmetry operations in real space in terms
!!   of primitive translations
!!
!! OUTPUT
!!  multi=multiplicity of the unit cell
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine chkprimit(chkprim,multi,nsym,symafm,symrel)

 use abi_defs_basis
 use abi_interfaces_lowlevel

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chkprim,nsym
 integer,intent(out) :: multi
!arrays
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)

!Local variables-------------------------------
!scalars
 integer :: isym
 character(len=500) :: message

!**************************************************************************

!Loop over each symmetry operation of the Bravais lattice
!Find whether it is the identity, or a pure translation,
!without change of sign of the spin
 multi=0
 do isym=1,nsym
   if( abs(symrel(1,1,isym)-1)+&
&   abs(symrel(2,2,isym)-1)+&
&   abs(symrel(3,3,isym)-1)+&
&   abs(symrel(1,2,isym))+abs(symrel(2,1,isym))+&
&   abs(symrel(2,3,isym))+abs(symrel(3,2,isym))+&
&   abs(symrel(3,1,isym))+abs(symrel(1,3,isym))+&
&   abs(symafm(isym)-1)                           ==0 )then
     multi=multi+1
   end if
 end do

!Check whether the cell is primitive
 if(multi>1)then
   if(chkprim/=0)then
     write(message,'(a,a,a,a,a,a,i3,a,a,a,a,a,a,a,a,a)') ch10,&
&     ' chkprimit : ERROR -',ch10,&
&     '  According to the symmetry finder, the unit cell is',ch10,&
&     '  NOT primitive. The multiplicity is',multi,' .',ch10,&
&     '  The use of non-primitive unit cells is allowed',ch10,&
&     '  only when the input variable chkprim is 0.',ch10,&
&     '  Action : either change your unit cell (rprim or angdeg),',ch10,&
&     '  or set chkprim to 0.'
     call abi_wrtout(std_out,message,'COLL')
     call abi_leave_new('COLL')
   else
     write(message,'(a,a,a,a,a,a,i3,a,a,a)') ch10,&
&     ' chkprimit : COMMENT -',ch10,&
&     '  According to the symmetry finder, the unit cell is',ch10,&
&     '  not primitive, with multiplicity=',multi,'.',ch10,&
&     '  This is allowed, as the input variable chkprim is 0.'
     call abi_wrtout(std_out,message,'COLL')
   end if
 end if

end subroutine chkprimit
!!***
