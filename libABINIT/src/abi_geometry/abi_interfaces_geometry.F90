!!****m* ABINIT/abi_interfaces_geometry
!! NAME
!! abi_interfaces_geometry
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/42_geometry
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

module abi_interfaces_geometry

 implicit none

interface
 subroutine abi_bldgrp(msym,nogen,nsym,symafm,symrel,tnons)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(inout) :: nogen
  integer,intent(in) :: nsym
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine abi_bldgrp
end interface

interface
 subroutine abi_bonds_lgth_angles(coordn,fnameabo_app_geo,natom,ntypat,&  
  &  rprimd,typat,xred,znucl)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: coordn
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  character(len=fnlen),intent(in) :: fnameabo_app_geo
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine abi_bonds_lgth_angles
end interface

interface
 subroutine abi_chkgrp(nsym,symafm,symrel)
  implicit none
  integer,intent(in) :: nsym
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine abi_chkgrp
end interface

interface
 subroutine abi_chkprimit(chkprim,multi,nsym,symafm,symrel)
  implicit none
  integer,intent(in) :: chkprim
  integer,intent(out) :: multi
  integer,intent(in) :: nsym
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine abi_chkprimit
end interface

interface
 subroutine abi_getptgroupma(ptgroup,ptgroupha,ptgroupma)
  implicit none
  integer,intent(out) :: ptgroupma
  character(len=5),intent(in) :: ptgroup
  character(len=5),intent(in) :: ptgroupha
 end subroutine abi_getptgroupma
end interface

interface
 subroutine abi_holocell(cell_base,foundc,iholohedry)
  use abi_defs_basis
  implicit none
  integer,intent(out) :: foundc
  integer,intent(in) :: iholohedry
  real(dp),intent(in) :: cell_base(3,3)
 end subroutine abi_holocell
end interface

interface
 subroutine metric(gmet,gprimd,iout,rmet,rprimd,ucvol)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: iout
  real(dp),intent(out) :: ucvol
  real(dp),intent(out) :: gmet(3,3)
  real(dp),intent(out) :: gprimd(3,3)
  real(dp),intent(out) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine metric
end interface

interface
 subroutine abi_mkrdim(acell,rprim,rprimd)
  use abi_defs_basis
  implicit none
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(out) :: rprimd(3,3)
 end subroutine abi_mkrdim
end interface

interface
 subroutine abi_ptgmadata(ptgroupma,ptgrpmasb)
  implicit none
  integer,intent(in) :: ptgroupma
  character(len=10),intent(out) :: ptgrpmasb
 end subroutine abi_ptgmadata
end interface

interface
 subroutine abi_smallprim(metmin,minim,rprimd)
  use abi_defs_basis
  implicit none
  real(dp),intent(out) :: metmin(3,3)
  real(dp),intent(out) :: minim(3,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine abi_smallprim
end interface

interface
 subroutine abi_spgdata(brvsb,intsb,intsbl,ptintsb,ptschsb,&  
  &  schsb,spgaxor,spgroup,sporder,spgorig)
  implicit none
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(out) :: sporder
  character(len=1),intent(out) :: brvsb
  character(len=15),intent(out) :: intsb
  character(len=35),intent(out) :: intsbl
  character(len=15),intent(out) :: ptintsb
  character(len=15),intent(out) :: ptschsb
  character(len=15),intent(out) :: schsb
 end subroutine abi_spgdata
end interface

interface
 subroutine abi_strainsym(nsym,rprimd0,rprimd,rprimd_symm,symrel)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: nsym
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rprimd0(3,3)
  real(dp),intent(out) :: rprimd_symm(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine abi_strainsym
end interface

interface
 subroutine abi_strconv(frac,gprimd,cart)
  use abi_defs_basis
  implicit none
  real(dp),intent(out) :: cart(6)
  real(dp),intent(in) :: frac(6)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine abi_strconv
end interface

interface
 subroutine abi_symanal(bravais,chkprim,genafm,msym,nsym,ptgroupma,rprimd,spgroup,symafm,symrel,tnons,tolsym)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: chkprim
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  integer,intent(out) :: ptgroupma
  integer,intent(out) :: spgroup
  real(dp),intent(in) :: tolsym
  integer,intent(out) :: bravais(11)
  real(dp),intent(out) :: genafm(3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(msym)
  integer,intent(in) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine abi_symanal
end interface

interface
 subroutine abi_symatm(indsym,natom,nsym,symrec,tnons,tolsym,typat,xred)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  real(dp), intent(in) :: tolsym
  integer,intent(out) :: indsym(4,nsym,natom)
  integer,intent(in) :: symrec(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine abi_symatm
end interface

interface
 subroutine abi_symaxes(center,iholohedry,&  
  &  isym,isymrelconv,ordersym,tnons_order,trialt,type_axis)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: center
  integer,intent(in) :: iholohedry
  integer,intent(in) :: isym
  integer,intent(in) :: ordersym
  integer,intent(in) :: tnons_order
  integer,intent(out) :: type_axis
  integer,intent(in) :: isymrelconv(3,3)
  real(dp),intent(in) :: trialt(3)
 end subroutine abi_symaxes
end interface

interface
 subroutine abi_symbrav(bravais,msym,nsym,ptgroup,rprimd,symrel,tolsym)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  character(len=5),intent(out) :: ptgroup
  real(dp),intent(in) :: tolsym
  integer,intent(out) :: bravais(11)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,msym)
 end subroutine abi_symbrav
end interface

interface
 subroutine abi_symchk(difmin,eatom,natom,tratom,transl,trtypat,typat,xred)
  use abi_defs_basis
  implicit none
  integer,intent(out) :: eatom
  integer,intent(in) :: natom
  integer,intent(in) :: trtypat
  integer,intent(out) :: transl(3)
  real(dp),intent(out) :: difmin(3)
  real(dp),intent(in) :: tratom(3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine abi_symchk
end interface

interface
 subroutine abi_symdet(determinant,nsym,sym)
  implicit none
  integer,intent(in) :: nsym
  integer,intent(out) :: determinant(nsym)
  integer,intent(in) :: sym(3,3,nsym)
 end subroutine abi_symdet
end interface

interface
 subroutine abi_symfind(berryopt,efield,gprimd,jellslab,msym,natom,noncoll,nptsym,nsym,&  
  &  ptsymrel,spinat,symafm,symrel,tnons,tolsym,typat,use_inversion,xred)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: jellslab
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: noncoll
  integer,intent(in) :: nptsym
  integer,intent(out) :: nsym
  integer,intent(in) :: use_inversion
  real(dp),intent(in) :: tolsym
  real(dp),intent(in) :: efield(3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: ptsymrel(3,3,msym)
  real(dp),intent(in) :: spinat(3,natom)
  integer,intent(out) :: symafm(msym)
  integer,intent(out) :: symrel(3,3,msym)
  real(dp),intent(out) :: tnons(3,msym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine abi_symfind
end interface

interface
 subroutine abi_symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(out) :: nptsym
  real(dp),intent(in) :: tolsym
  integer,intent(out) :: bravais(11)
  integer,intent(out) :: ptsymrel(3,3,msym)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine abi_symlatt
end interface

interface
 subroutine abi_symlist_bcc(additional_info,nsym,n_axes,spgroup)
  implicit none
  integer,intent(in) :: additional_info
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  integer,intent(in) :: n_axes(31)
 end subroutine abi_symlist_bcc
end interface

interface
 subroutine abi_symlist_fcc(nsym,n_axes,spgroup)
  implicit none
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  integer,intent(in) :: n_axes(31)
 end subroutine abi_symlist_fcc
end interface

interface
 subroutine abi_symlist_others(brvltt,nsym,n_axes,spgroup)
  implicit none
  integer,intent(in) :: brvltt
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  integer,intent(in) :: n_axes(31)
 end subroutine abi_symlist_others
end interface

interface
 subroutine abi_symlist_prim(additional_info,nsym,n_axes,spgroup)
  implicit none
  integer,intent(in) :: additional_info
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  integer,intent(in) :: n_axes(31)
 end subroutine abi_symlist_prim
end interface

interface
 subroutine abi_symplanes(center,iholohedry,isym,isymrelconv,itnonsconv,type_axis)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: center
  integer,intent(in) :: iholohedry
  integer,intent(in) :: isym
  integer,intent(out) :: type_axis
  integer,intent(in) :: isymrelconv(3,3)
  real(dp),intent(in) :: itnonsconv(3)
 end subroutine abi_symplanes
end interface

interface
 subroutine abi_symptgroup(iholohedry,nsym,ptgroup,symrel)
  implicit none
  integer,intent(out) :: iholohedry
  integer,intent(in) :: nsym
  character(len=5),intent(out) :: ptgroup
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine abi_symptgroup
end interface

interface
 subroutine abi_symrelrot(nsym,rprimd,rprimd_new,symrel,tolsym)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: nsym
  real(dp),intent(in) :: tolsym
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rprimd_new(3,3)
  integer,intent(inout) :: symrel(3,3,nsym)
 end subroutine abi_symrelrot
end interface

interface
 subroutine abi_symsgcube(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
  &  spgroupma,symafm,symrel,tnons)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(out) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(out) :: symafm(msym)
  integer,intent(out) :: symrel(3,3,msym)
  real(dp),intent(out) :: tnons(3,msym)
 end subroutine abi_symsgcube
end interface

interface
 subroutine abi_symsghexa(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
  &  spgroupma,symafm,symrel,tnons)
  use abi_defs_basis
  implicit none
  integer,intent(out) :: brvltt
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(out) :: symafm(msym)
  integer,intent(out) :: symrel(3,3,msym)
  real(dp),intent(out) :: tnons(3,msym)
 end subroutine abi_symsghexa
end interface

interface
 subroutine abi_symsgmono(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
  &  spgroupma,symafm,symrel,tnons)
  use abi_defs_basis
  implicit none
  integer,intent(out) :: brvltt
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(out) :: symafm(msym)
  integer,intent(out) :: symrel(3,3,msym)
  real(dp),intent(out) :: tnons(3,msym)
 end subroutine abi_symsgmono
end interface

interface
 subroutine abi_symsgortho(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
  &  spgroupma,symafm,symrel,tnons)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(out) :: symafm(msym)
  integer,intent(out) :: symrel(3,3,msym)
  real(dp),intent(out) :: tnons(3,msym)
 end subroutine abi_symsgortho
end interface

interface
 subroutine abi_symsgtetra(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
  &  spgroupma,symafm,symrel,tnons)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(out) :: symafm(msym)
  integer,intent(out) :: symrel(3,3,msym)
  real(dp),intent(out) :: tnons(3,msym)
 end subroutine abi_symsgtetra
end interface

interface
 subroutine abi_symspgr(bravais,nsym,spgroup,symrel,tnons,tolsym)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  real(dp),intent(in) :: tolsym
  integer,intent(in) :: bravais(11)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(inout) :: tnons(3,nsym)
 end subroutine abi_symspgr
end interface

interface
 subroutine abi_xredxcart(natom,option,rprimd,xcart,xred)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: option
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xcart(3,natom)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine abi_xredxcart
end interface

end module abi_interfaces_geometry
!!***
