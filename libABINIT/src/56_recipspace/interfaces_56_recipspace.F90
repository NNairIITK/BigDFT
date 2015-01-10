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
 subroutine bound(dsqmax,dsqmin,gbound,gmet,kpt,ngfft,plane)
  use defs_basis
  implicit none
  integer,intent(out) :: plane
  real(dp),intent(out) :: dsqmax
  real(dp),intent(out) :: dsqmin
  integer,intent(out) :: gbound(3)
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: kpt(3)
 end subroutine bound
end interface

interface
 subroutine bound_new(dsqmax,dsqmin,gbound,gmet,kpt,ngfft,plane)
  use defs_basis
  implicit none
  integer,intent(out) :: plane
  real(dp),intent(out) :: dsqmax
  real(dp),intent(out) :: dsqmin
  integer,intent(out) :: gbound(3)
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: kpt(3)
 end subroutine bound_new
end interface

interface
 subroutine getfullg(nbase,nsym,pinv,sizepw,gbase,symrec,cnorm,maxpw,gbig,shlim,ierr)
  use defs_basis
  implicit none
  integer,intent(out) :: ierr
  integer,intent(out) :: maxpw
  integer,intent(in) :: nbase
  integer,intent(in) :: nsym
  integer,intent(in) :: pinv
  integer,intent(in) :: sizepw
  real(dp),intent(in) :: cnorm(nbase)
  integer,intent(in) :: gbase(3,nbase)
  integer,intent(out) :: gbig(3,sizepw)
  integer,intent(out) :: shlim(nbase)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine getfullg
end interface

interface
 subroutine get_full_kgrid(indkpt,klatt,kpt,kpt_fullbz,kptrlatt,nkpt,&  
  &  nkpt_fullbz,nshiftk,nsym,shiftk,symrel)
  use defs_basis
  implicit none
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkpt_fullbz
  integer,intent(in) :: nshiftk
  integer,intent(in) :: nsym
  integer,intent(in) :: kptrlatt(3,3)
  integer,intent(out) :: indkpt(nkpt_fullbz)
  real(dp),intent(in) :: klatt(3,3)
  real(dp),intent(in) :: kpt(3,nkpt)
  real(dp),intent(out) :: kpt_fullbz(3,nkpt_fullbz)
  real(dp),intent(in) :: shiftk(3,nshiftk)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine get_full_kgrid
end interface

interface
 subroutine get_irredg(npw_k,nsym,pinv,gprimd,symrec,gcurr,nbasek,gbasek,cnormk)
  use defs_basis
  implicit none
  integer,intent(out) :: nbasek
  integer,intent(in) :: npw_k
  integer,intent(in) :: nsym
  integer,intent(in) :: pinv
  real(dp),intent(out) :: cnormk(npw_k)
  integer,intent(out) :: gbasek(3,npw_k)
  integer,intent(in) :: gcurr(3,npw_k)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine get_irredg
end interface

interface
 subroutine get_tetra (indkpt,gprimd,klatt,kpt_fullbz,mtetra,nkpt_fullbz,&  
  &  ntetra,tetra_full,tetra_mult,tetra_wrap,vv)
  use defs_basis
  implicit none
  integer,intent(in) :: mtetra
  integer,intent(in) :: nkpt_fullbz
  integer,intent(out) :: ntetra
  real(dp),intent(out) :: vv
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indkpt(nkpt_fullbz)
  real(dp),intent(in) :: klatt(3,3)
  real(dp),intent(in) :: kpt_fullbz(3,nkpt_fullbz)
  integer,intent(out) :: tetra_full(4,2,mtetra)
  integer,intent(out) :: tetra_mult(mtetra)
  integer,intent(out) :: tetra_wrap(3,4,mtetra)
 end subroutine get_tetra
end interface

interface
 subroutine getcut(boxcut,ecut,gmet,gsqcut,iboxcut,iout,kpt,ngfft)
  use defs_basis
  implicit none
  integer,intent(in) :: iboxcut
  integer,intent(in) :: iout
  real(dp),intent(out) :: boxcut
  real(dp),intent(in) :: ecut
  real(dp),intent(out) :: gsqcut
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: kpt(3)
 end subroutine getcut
end interface

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
 subroutine getkpgnorm(gprimd,kpt,kg_k,kpgnorm,npw_k)
  use defs_basis
  implicit none
  integer,intent(in) :: npw_k
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg_k(3,npw_k)
  real(dp),intent(out) :: kpgnorm(npw_k)
  real(dp),intent(in) :: kpt(3)
 end subroutine getkpgnorm
end interface

interface
 subroutine getmpw(ecut,exchn2n3d,gmet,istwfk,kptns,mpi_enreg,mpw,nkpt)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: exchn2n3d
  integer,intent(out) :: mpw
  integer,intent(in) :: nkpt
  real(dp),intent(in) :: ecut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: istwfk(nkpt)
  real(dp),intent(in) :: kptns(3,nkpt)
 end subroutine getmpw
end interface

interface
 subroutine getng(boxcutmin,ecut,gmet,me_fft,mgfft,nfft,ngfft,nproc_fft,nsym,option_lob,paral_fft,symrel)
  use defs_basis
  implicit none
  integer,intent(in) :: me_fft
  integer,intent(out) :: mgfft
  integer,intent(out) :: nfft
  integer,intent(in) :: nproc_fft
  integer,intent(in) :: nsym
  integer,intent(in) :: option_lob
  integer,intent(in) :: paral_fft
  real(dp),intent(in) :: boxcutmin
  real(dp),intent(in) :: ecut
  integer,intent(inout) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine getng
end interface

interface
 subroutine getph(atindx,natom,n1,n2,n3,ph1d,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: natom
  integer,intent(in) :: atindx(natom)
  real(dp),intent(out) :: ph1d(2,(2*n1+1+2*n2+1+2*n3+1)*natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine getph
end interface

interface
 subroutine getwtk(kpt,nkpt,nsym,symrel,wtk)
  use defs_basis
  implicit none
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsym
  real(dp),intent(in) :: kpt(3,nkpt)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(out) :: wtk(nkpt)
 end subroutine getwtk
end interface

interface
 subroutine initylmg(gprimd,kg,kptns,mkmem,mpi_enreg,mpsang,mpw,nband,nkpt,&  
  &  npwarr,nsppol,optder,rprimd,unkg,unylm,ylm,ylm_gr)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: optder
  integer,intent(in) :: unkg
  integer,intent(in) :: unylm
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: ylm(mpw*mkmem,mpsang*mpsang)
  real(dp),intent(out) :: ylm_gr(mpw*mkmem,3+6*(optder/2),mpsang*mpsang)
 end subroutine initylmg
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
 subroutine kpgio(ecut,exchn2n3d,gmet,istwfk,kg,kgnam,kptns,mkmem,nband,nkpt,&  
  &  mode_paral,mpi_enreg,mpw,npwarr,npwtot,nsppol,unkg)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: exchn2n3d
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: unkg
  real(dp),intent(in) :: ecut
  character(len=fnlen),intent(in) :: kgnam
  character(len=4),intent(in) :: mode_paral
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: istwfk(nkpt)
  integer,intent(out) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(out) :: npwarr(nkpt)
  integer,intent(out) :: npwtot(nkpt)
 end subroutine kpgio
end interface

interface
 subroutine kpgsph(ecut,exchn2n3d,gmet,ikg,ikpt,istwf_k,kg,kpt,mkmem,mpi_enreg,mpw,npw)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: exchn2n3d
  integer,intent(in) :: ikg
  integer,intent(in) :: ikpt
  integer,intent(in) :: istwf_k
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(out) :: npw
  real(dp),intent(in) :: ecut
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(out) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kpt(3)
 end subroutine kpgsph
end interface

interface
 subroutine laplacian(gprimd,mpi_enreg,nfft,nfunc,ngfft,paral_kgb,rdfuncr,&  
  &  laplacerdfuncr,rdfuncg_out,laplacerdfuncg_out,g2cart_out,rdfuncg_in,g2cart_in)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: nfunc
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out),optional,target :: g2cart_in(nfft)
  real(dp),intent(out),optional,target :: g2cart_out(nfft)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out),optional,target :: laplacerdfuncg_out(2,nfft,nfunc)
  real(dp),intent(inout),optional :: laplacerdfuncr(nfft,nfunc)
  real(dp),intent(out),optional,target :: rdfuncg_in(2,nfft,nfunc)
  real(dp),intent(out),optional,target :: rdfuncg_out(2,nfft,nfunc)
  real(dp),intent(inout),optional,target :: rdfuncr(nfft,nfunc)
 end subroutine laplacian
end interface

interface
 subroutine merge_kgirr(nsym,pinv,nkpt,mpw,sizepw,symrec,nbasek,cnormk,gbasek,nbase,gbase,cnorm,ierr)
  use defs_basis
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: mpw
  integer,intent(out) :: nbase
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsym
  integer,intent(in) :: pinv
  integer,intent(in) :: sizepw
  real(dp),intent(out) :: cnorm(sizepw)
  real(dp),intent(in) :: cnormk(mpw,nkpt)
  integer,intent(out) :: gbase(3,sizepw)
  integer,intent(in) :: gbasek(3,mpw,nkpt)
  integer,intent(in) :: nbasek(nkpt)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine merge_kgirr
end interface

interface
 subroutine mkkin (ecut,ecutsm,effmass,gmet,kg,kinpw,kpt,npw)
  use defs_basis
  implicit none
  integer,intent(in) :: npw
  real(dp),intent(in) :: ecut
  real(dp),intent(in) :: ecutsm
  real(dp),intent(in) :: effmass
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg(3,npw)
  real(dp),intent(out) :: kinpw(npw)
  real(dp),intent(in) :: kpt(3)
 end subroutine mkkin
end interface

interface
 subroutine pmat2cart(eigen11,eigen12,eigen13,mband,nkpt,nsppol,pmat,rprimd)
  use defs_basis
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  real(dp),intent(in) :: eigen11(2,mband,mband,nkpt,nsppol)
  real(dp),intent(in) :: eigen12(2,mband,mband,nkpt,nsppol)
  real(dp),intent(in) :: eigen13(2,mband,mband,nkpt,nsppol)
  complex(dpc),intent(out) :: pmat(mband,mband,nkpt,3,nsppol)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine pmat2cart
end interface

interface
 subroutine setshells(ecut,npw,nsh,nsym,gmet,gprimd,symrel,tag,ucvol)
  use defs_basis
  implicit none
  integer,intent(inout) :: npw
  integer,intent(inout) :: nsh
  integer,intent(in) :: nsym
  real(dp),intent(inout) :: ecut
  character(len=*),intent(in) :: tag
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine setshells
end interface

interface
 subroutine setsym(indsym,irrzon,iscf,natom,&  
  &  nfft,ngfft,nspden,nsppol,nsym,phnons,&  
  &  symafm,symrec,symrel,tnons,typat,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: iscf
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ngfft(18)
  integer,intent(out) :: indsym(4,nsym,natom)
  integer,intent(out) :: irrzon(nfft,2,(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(out) :: phnons(2,nfft,(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: symafm(nsym)
  integer,intent(out) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine setsym
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
 subroutine symg(kg_diel,npwdiel,nsym,phdiel,sym_g,symrel,tmrev_g,tnons)
  use defs_basis
  implicit none
  integer,intent(in) :: npwdiel
  integer,intent(in) :: nsym
  integer,intent(in) :: kg_diel(3,npwdiel)
  real(dp),intent(out) :: phdiel(2,npwdiel,nsym)
  integer,intent(out) :: sym_g(npwdiel,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  integer,intent(out) :: tmrev_g(npwdiel)
  real(dp),intent(in) :: tnons(3,nsym)
 end subroutine symg
end interface

interface
 subroutine symkchk(kptns,nkpt,nsym,symrec,timrev)
  use defs_basis
  implicit none
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsym
  integer,intent(in) :: timrev
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine symkchk
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
 subroutine symq3(nsym,qpt,symq,symrec,timrev,prtvol)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  integer,intent(in),optional :: prtvol
  integer,intent(out) :: timrev
  real(dp),intent(in) :: qpt(3)
  integer,intent(out) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine symq3
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
