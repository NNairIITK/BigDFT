!!****m* ABINIT/interfaces_56_recipspace
!! NAME
!! interfaces_56_recipspace
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/56_recipspace
!!
!! COPYRIGHT
!! Copyright (C) 2010 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!!
!! SOURCE

module interfaces_56_recipspace

 implicit none

interface
 subroutine getkgrid(iout,iscf,kpt,kptopt,kptrlatt,kptrlen,&  
  &  msym,nkpt,nkpt_computed,nshiftk,nsym,rprimd,shiftk,symafm,&  
  &  symrel,vacuum,wtk,kbz_p)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: iscf
  integer,intent(in) :: kptopt
  integer,intent(in) :: msym
  integer,intent(in) :: nkpt
  integer,intent(out) :: nkpt_computed
  integer,intent(inout) :: nshiftk
  integer,intent(in) :: nsym
  real(dp),intent(out) :: kptrlen
  integer,intent(inout) :: kptrlatt(3,3)
  integer,intent(in) :: vacuum(3)
  real(dp),optional,pointer :: kbz_p(:,:)
  real(dp),intent(out) :: kpt(3,nkpt)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: shiftk(3,8)
  integer,intent(in) :: symafm(msym)
  integer,intent(in) :: symrel(3,3,msym)
  real(dp),intent(out) :: wtk(nkpt)
 end subroutine getkgrid
end interface

interface
 subroutine irrzg(irrzon,nspden,nsppol,nsym,n1,n2,n3,phnons,&  
  &  symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(out) :: irrzon(n1*n2*n3,2,(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(out) :: phnons(2,n1*n2*n3,(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
 end subroutine irrzg
end interface

interface
 subroutine smpbz(brav,iout,kptrlatt,mkpt,nkpt,nshiftk,option,shiftk,spkpt)
  use defs_basis
  implicit none
  integer,intent(in) :: brav
  integer,intent(in) :: iout
  integer,intent(in) :: mkpt
  integer,intent(out) :: nkpt
  integer,intent(in) :: nshiftk
  integer,intent(in) :: option
  integer,intent(in) :: kptrlatt(3,3)
  real(dp),intent(in) :: shiftk(3,nshiftk)
  real(dp),intent(out) :: spkpt(3,mkpt)
 end subroutine smpbz
end interface

interface
 subroutine symkpt(gmet,indkpt1,kptns,nkpt,nkpt1,nsym1,option,&  
  &  symrc1,timrev,wtk,wtk_folded)
  use defs_basis
  implicit none
  integer,intent(in) :: nkpt
  integer,intent(out) :: nkpt1
  integer,intent(in) :: nsym1
  integer,intent(in) :: option
  integer,intent(in) :: timrev
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(out) :: indkpt1(nkpt)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: symrc1(3,3,nsym1)
  real(dp),intent(in) :: wtk(nkpt)
  real(dp),intent(out) :: wtk_folded(nkpt)
 end subroutine symkpt
end interface

interface
 subroutine testkgrid(bravais,iout,kptrlatt,kptrlen,&  
  &  msym,nshiftk,nsym,prtkpt,rprimd,shiftk,symafm,symrel,vacuum)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: msym
  integer,intent(out) :: nshiftk
  integer,intent(in) :: nsym
  integer,intent(in) :: prtkpt
  real(dp),intent(inout) :: kptrlen
  integer,intent(in) :: bravais(11)
  integer,intent(out) :: kptrlatt(3,3)
  integer,intent(in) :: vacuum(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: shiftk(3,8)
  integer,intent(in) :: symafm(msym)
  integer,intent(in) :: symrel(3,3,msym)
 end subroutine testkgrid
end interface

end module interfaces_56_recipspace
!!***
