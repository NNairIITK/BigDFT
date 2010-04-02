!{\src2tex{textfont=tt}}
!!****f* ABINIT/xredxcart
!! NAME
!! xredxcart
!!
!! FUNCTION
!! If option==1 :
!! Convert from dimensionless reduced coordinates xred(3,natom)
!! to cartesian coordinates xcart(3,natom) in bohr by using
!! xcart(mu,ia)=rprimd(mu,1)*xred(1,ia)
!!             +rprimd(mu,2)*xred(2,ia)
!!             +rprimd(mu,3)*xred(3,ia)
!!
!! If option==-1
!! Convert from cartesian coordinates xcart(3,natom) in bohr to
!! dimensionless reduced coordinates xred(3,natom) by using
!! xred(mu,ia)=(gprimd(1,mu)*xcart(1,ia)+gprimd(2,mu)*xcart(2,ia)+
!!       gprimd(3,mu)*xcart(3,ia))
!! where gprimd is the inverse of rprimd
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  natom=number of atoms in unit cell
!!  option=see above
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output (see above):
!!  xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!  xred(3,natom)=dimensionless reduced coordinates of atoms
!!
!! PARENTS
!!      afterscfloop,berryphase,berryphase_new,bonds_lgth_angles,brdmin,constrf
!!      delocint,diisrelax,driver,hirsh,ingeo,ionion_realspace,jvec_to_B
!!      localorb_S,m_crystal,make_efg_el,make_efg_ion,mklocl,mklocl_realspace
!!      moldyn,move,out1dm,out_geometry_xml,outqmc,outvars
!!      partial_dos_fractions,prcref,prcref_PMA,prtspgroup,rdddb9,relaxpol
!!      scphon,spin_current,symspgr,vtorho,wffile,wvl_init_type_proj
!!      wvl_init_type_wfs,wvl_memory,wvl_rwwf,wvl_setboxgeometry,wvl_vtorho
!!      wvl_wfsinp_reformat,wvl_wfsinp_scratch
!!
!! CHILDREN
!!      leave_new,matr3inv,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

subroutine xredxcart(natom,option,rprimd,xcart,xred)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,option
!arrays
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: xcart(3,natom),xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,mu
 character(len=500) :: message
!arrays
 real(dp) :: gprimd(3,3)

! *************************************************************************

 if(option==1)then
   do iatom=1,natom
     do mu=1,3
       xcart(mu,iatom)=rprimd(mu,1)*xred(1,iatom)+rprimd(mu,2)*xred(2,iatom)+&
&       rprimd(mu,3)*xred(3,iatom)
     end do
   end do
 else if(option==-1)then
   call matr3inv(rprimd,gprimd)
   do iatom=1,natom
     do mu=1,3
       xred(mu,iatom)= gprimd(1,mu)*xcart(1,iatom)+gprimd(2,mu)*xcart(2,iatom)+&
&       gprimd(3,mu)*xcart(3,iatom)
     end do
   end do
 else
   write(message, '(a,a,a,a,i4,a)' ) ch10,&
&   ' xredxcart : BUG -',ch10,&
&   '  Option must be 1 or -1, while it is ',option,'.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

end subroutine xredxcart
!!***
