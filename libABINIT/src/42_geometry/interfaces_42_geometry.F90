!!****m* ABINIT/interfaces_42_geometry
!! NAME
!! interfaces_42_geometry
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

module interfaces_42_geometry

 implicit none

interface
 subroutine bldgrp(msym,nogen,nsym,symafm,symrel,tnons)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(inout) :: nogen
  integer,intent(in) :: nsym
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine bldgrp
end interface

interface
 subroutine bldgrpaf(msym,nogenaf,nsym,symafm,symrel,symrel_magn,tnons)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(in) :: nogenaf
  integer,intent(in) :: nsym
  integer,intent(in) :: symrel_magn(3,3,3)
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine bldgrpaf
end interface

interface
 subroutine bonds_lgth_angles(coordn,fnameabo_app_geo,natom,ntypat,&  
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
 end subroutine bonds_lgth_angles
end interface

interface
 subroutine chkdilatmx(dilatmx,rprimd,rprimd_orig)
  use abi_defs_basis
  implicit none
  real(dp),intent(in) :: dilatmx
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rprimd_orig(3,3)
 end subroutine chkdilatmx
end interface

interface
 subroutine chkgrp(nsym,symafm,symrel)
  implicit none
  integer,intent(in) :: nsym
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine chkgrp
end interface

interface
 subroutine chkorthsy(gprimd,iout,nsym,rmet,rprimd,symrel)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: nsym
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine chkorthsy
end interface

interface
 subroutine chkprimit(chkprim,multi,nsym,symafm,symrel)
  implicit none
  integer,intent(in) :: chkprim
  integer,intent(out) :: multi
  integer,intent(in) :: nsym
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine chkprimit
end interface

interface
 subroutine gensymshub(genafm,spgroup,spgroupma,shubnikov)
  use abi_defs_basis
  implicit none
  integer,intent(out) :: shubnikov
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  real(dp),intent(out) :: genafm(3)
 end subroutine gensymshub
end interface

interface
 subroutine gensymshub4(genafm,msym,nsym,symafm,symrel,tnons)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(inout) :: nsym
  real(dp),intent(in) :: genafm(3)
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine gensymshub4
end interface

interface
 subroutine gensymspgr(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,&  
  &  spgroup,spgroupma,symafm,symrel,tnons)
  use abi_defs_basis
  implicit none
  integer,intent(inout) :: brvltt
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
 end subroutine gensymspgr
end interface

interface
 subroutine getptgroupma(ptgroup,ptgroupha,ptgroupma)
  implicit none
  integer,intent(out) :: ptgroupma
  character(len=5),intent(in) :: ptgroup
  character(len=5),intent(in) :: ptgroupha
 end subroutine getptgroupma
end interface

interface
 subroutine getspinrot(rprimd,spinrot,symrel_conv)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: symrel_conv(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: spinrot(4)
 end subroutine getspinrot
end interface

interface
 subroutine gridgcart(gcart,gprimd,ngfft)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: gcart(ngfft(1),ngfft(2),ngfft(3),3)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine gridgcart
end interface

interface
 subroutine holocell(cell_base,foundc,iholohedry)
  use abi_defs_basis
  implicit none
  integer,intent(out) :: foundc
  integer,intent(in) :: iholohedry
  real(dp),intent(in) :: cell_base(3,3)
 end subroutine holocell
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
 subroutine mkrdim(acell,rprim,rprimd)
  use abi_defs_basis
  implicit none
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(out) :: rprimd(3,3)
 end subroutine mkrdim
end interface

interface
 subroutine operat(natom,natrd,nsym,spinat,symafm,symrel,tnons,typat,xred)
  use abi_defs_basis
  use abi_interfaces_lowlevel
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: natrd
  integer,intent(in) :: nsym
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  integer,intent(inout) :: typat(natom)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(inout) :: spinat(3,natom)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine operat
end interface

interface
 subroutine prtspgroup(bravais,genafm,iout,jdtset,ptgroupma,spgroup)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: jdtset
  integer,intent(in) :: ptgroupma
  integer,intent(in) :: spgroup
  integer,intent(in) :: bravais(11)
  real(dp),intent(inout) :: genafm(3)
 end subroutine prtspgroup
end interface

interface
 subroutine ptgmadata(ptgroupma,ptgrpmasb)
  implicit none
  integer,intent(in) :: ptgroupma
  character(len=10),intent(out) :: ptgrpmasb
 end subroutine ptgmadata
end interface

interface
 subroutine smallprim(metmin,minim,rprimd)
  use abi_defs_basis
  implicit none
  real(dp),intent(out) :: metmin(3,3)
  real(dp),intent(out) :: minim(3,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine smallprim
end interface

interface
 subroutine spgdata(brvsb,intsb,intsbl,ptintsb,ptschsb,&  
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
 end subroutine spgdata
end interface

interface
 subroutine strainsym(nsym,rprimd0,rprimd,rprimd_symm,symrel)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: nsym
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rprimd0(3,3)
  real(dp),intent(out) :: rprimd_symm(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine strainsym
end interface

interface
 subroutine strconv(frac,gprimd,cart)
  use abi_defs_basis
  implicit none
  real(dp),intent(out) :: cart(6)
  real(dp),intent(in) :: frac(6)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine strconv
end interface

interface
 subroutine stresssym(gprimd,nsym,stress,sym)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: nsym
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(inout) :: stress(6)
  integer,intent(in) :: sym(3,3,nsym)
 end subroutine stresssym
end interface

interface
 subroutine sym2cart(gprimd,nsym,rprimd,symrel,symcart)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: nsym
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: symcart(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine sym2cart
end interface

interface
 subroutine symanal(bravais,chkprim,genafm,msym,nsym,ptgroupma,rprimd,spgroup,symafm,symrel,tnons,tolsym)
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
 end subroutine symanal
end interface

interface
 subroutine symatm(indsym,natom,nsym,symrec,tnons,tolsym,typat,xred)
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
 end subroutine symatm
end interface

interface
 subroutine symaxes(center,iholohedry,&  
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
 end subroutine symaxes
end interface

interface
 subroutine symbrav(bravais,msym,nsym,ptgroup,rprimd,symrel,tolsym)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  character(len=5),intent(out) :: ptgroup
  real(dp),intent(in) :: tolsym
  integer,intent(out) :: bravais(11)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,msym)
 end subroutine symbrav
end interface

interface
 subroutine symchk(difmin,eatom,natom,tratom,transl,trtypat,typat,xred)
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
 end subroutine symchk
end interface

interface
 subroutine symdet(determinant,nsym,sym)
  implicit none
  integer,intent(in) :: nsym
  integer,intent(out) :: determinant(nsym)
  integer,intent(in) :: sym(3,3,nsym)
 end subroutine symdet
end interface

interface
 subroutine symfind(berryopt,efield,gprimd,jellslab,msym,natom,noncoll,nptsym,nsym,&  
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
 end subroutine symfind
end interface

interface
 subroutine symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(out) :: nptsym
  real(dp),intent(in) :: tolsym
  integer,intent(out) :: bravais(11)
  integer,intent(out) :: ptsymrel(3,3,msym)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine symlatt
end interface

interface
 subroutine symlist_bcc(additional_info,nsym,n_axes,spgroup)
  implicit none
  integer,intent(in) :: additional_info
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  integer,intent(in) :: n_axes(31)
 end subroutine symlist_bcc
end interface

interface
 subroutine symlist_fcc(nsym,n_axes,spgroup)
  implicit none
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  integer,intent(in) :: n_axes(31)
 end subroutine symlist_fcc
end interface

interface
 subroutine symlist_others(brvltt,nsym,n_axes,spgroup)
  implicit none
  integer,intent(in) :: brvltt
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  integer,intent(in) :: n_axes(31)
 end subroutine symlist_others
end interface

interface
 subroutine symlist_prim(additional_info,nsym,n_axes,spgroup)
  implicit none
  integer,intent(in) :: additional_info
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  integer,intent(in) :: n_axes(31)
 end subroutine symlist_prim
end interface

interface
 subroutine symmultsg(nsym,symafm,symrel,tnons)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: nsym
  integer :: symafm(nsym)
  integer :: symrel(3,3,nsym)
  real(dp) :: tnons(3,nsym)
 end subroutine symmultsg
end interface

interface
 subroutine symplanes(center,iholohedry,isym,isymrelconv,itnonsconv,type_axis)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: center
  integer,intent(in) :: iholohedry
  integer,intent(in) :: isym
  integer,intent(out) :: type_axis
  integer,intent(in) :: isymrelconv(3,3)
  real(dp),intent(in) :: itnonsconv(3)
 end subroutine symplanes
end interface

interface
 subroutine symptgroup(iholohedry,nsym,ptgroup,symrel)
  implicit none
  integer,intent(out) :: iholohedry
  integer,intent(in) :: nsym
  character(len=5),intent(out) :: ptgroup
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine symptgroup
end interface

interface
 subroutine symredcart(aprim,bprim,symcart,symred)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: symred(3,3)
  real(dp),intent(in) :: aprim(3,3)
  real(dp),intent(in) :: bprim(3,3)
  real(dp),intent(out) :: symcart(3,3)
 end subroutine symredcart
end interface

interface
 subroutine symrelrot(nsym,rprimd,rprimd_new,symrel,tolsym)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: nsym
  real(dp),intent(in) :: tolsym
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rprimd_new(3,3)
  integer,intent(inout) :: symrel(3,3,nsym)
 end subroutine symrelrot
end interface

interface
 subroutine symsgcube(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
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
 end subroutine symsgcube
end interface

interface
 subroutine symsghexa(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
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
 end subroutine symsghexa
end interface

interface
 subroutine symsgmono(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
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
 end subroutine symsgmono
end interface

interface
 subroutine symsgortho(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
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
 end subroutine symsgortho
end interface

interface
 subroutine symsgtetra(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
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
 end subroutine symsgtetra
end interface

interface
 subroutine symspgr(bravais,nsym,spgroup,symrel,tnons,tolsym)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  real(dp),intent(in) :: tolsym
  integer,intent(in) :: bravais(11)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(inout) :: tnons(3,nsym)
 end subroutine symspgr
end interface

interface
 subroutine symzat(indsym,natom,nsym,symrel,tnons,xred)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine symzat
end interface

interface
 subroutine xredxcart(natom,option,rprimd,xcart,xred)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: option
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xcart(3,natom)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine xredxcart
end interface

end module interfaces_42_geometry
!!***
