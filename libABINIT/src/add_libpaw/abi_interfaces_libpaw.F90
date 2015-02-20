!!****m* ABINIT/interfaces_42_libpaw
!! NAME
!! interfaces_42_libpaw
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/42_libpaw
!!
!! COPYRIGHT
!! Copyright (C) 2010-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!!
!!
!! SOURCE

#include "../libpaw/libpaw.h"

module abi_interfaces_libpaw

 implicit none

 interface
    subroutine pawinit(gnt_option,gsqcut_eff,lcutdens,lmix,mpsang,nphi,nsym,ntheta,&
         &                  pawang,pawrad,pawspnorb,pawtab,pawxcdev,xclevel,usepotzero)

      USE_DEFS
      use m_pawang, only : pawang_type
      use m_pawrad, only : pawrad_type
      use m_pawtab, only : pawtab_type

      implicit none

      !Arguments ---------------------------------------------
      !scalars
      integer,intent(in) :: gnt_option,lcutdens,lmix,mpsang,nphi,nsym,ntheta
      integer,intent(in) :: pawspnorb,pawxcdev,xclevel,usepotzero
      real(dp),intent(in) :: gsqcut_eff
      type(pawang_type),intent(inout) :: pawang
      !arrays
      type(pawrad_type),intent(in) :: pawrad(:)
      type(pawtab_type),target,intent(inout) :: pawtab(:)
    end subroutine pawinit
 end interface

 interface
    subroutine initrhoij(cplex,lexexch,lpawu,my_natom,natom,&
         &                    nspden,nspinor,nsppol,ntypat,pawrhoij,pawspnorb,pawtab,spinat,typat,&
         &                    ngrhoij,nlmnmix,use_rhoij_,use_rhoijres,& ! optional arguments
         &                    mpi_atmtab,mpi_comm_atom) ! optional arguments (parallelism)

      USE_DEFS
      use m_pawtab,      only : pawtab_type
      use m_pawrhoij,    only : pawrhoij_type

      implicit none

      !Arguments ---------------------------------------------
      !scalars
      integer,intent(in) :: cplex,my_natom,natom,nspden,nspinor,nsppol,ntypat,pawspnorb
      integer,intent(in),optional :: mpi_comm_atom,ngrhoij,nlmnmix,use_rhoij_,use_rhoijres
      character(len=500) :: message
      !arrays
      integer,intent(in) :: lexexch(ntypat),lpawu(ntypat)
      integer,intent(in) :: typat(natom)
      integer,optional,target,intent(in) :: mpi_atmtab(:)
      real(dp),intent(in) :: spinat(3,natom)
      type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
      type(pawtab_type),intent(in) :: pawtab(ntypat)
    end subroutine initrhoij
 end interface

 interface
    subroutine pawnhatfr(ider,idir,ipert,my_natom,natom,nspden,ntypat,&
         &                    pawang,pawfgrtab,pawrhoij,pawtab,rprimd, &
         &                    mpi_atmtab,mpi_comm_atom) ! optional arguments (parallelism)

      USE_DEFS

      use m_pawang,       only : pawang_type
      use m_pawtab,       only : pawtab_type
      use m_pawfgrtab,    only : pawfgrtab_type
      use m_pawrhoij,     only : pawrhoij_type
      
      implicit none

      !Arguments ------------------------------------
      !scalars
      integer,intent(in) :: ider,idir,ipert,my_natom,natom,nspden,ntypat
      integer,optional,intent(in) :: mpi_comm_atom
      type(pawang_type),intent(in) :: pawang
      !arrays
      integer,optional,target,intent(in) :: mpi_atmtab(:)
      real(dp),intent(in) :: rprimd(3,3)
      type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
      type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
      type(pawtab_type),intent(in) :: pawtab(ntypat)
    end subroutine pawnhatfr
 end interface

 interface
    subroutine mean_fftr(arraysp,meansp,nfft,nfftot,nspden,mpi_comm_sphgrid)

      USE_DEFS
      implicit none

      !Arguments ------------------------------------
      !scalars
      integer,intent(in) :: nfft,nfftot,nspden
      integer,intent(in),optional:: mpi_comm_sphgrid
      !arrays
      real(dp),intent(in) :: arraysp(nfft,nspden)
      real(dp),intent(out) :: meansp(nspden)
    end subroutine mean_fftr
 end interface

 interface
    subroutine pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,&
         &          my_natom,natom,nfft,ngfft,nhatgrdim,nspden,ntypat,pawang,pawfgrtab,&
         &          pawgrnhat,pawnhat,pawrhoij,pawrhoij0,pawtab,qphon,rprimd,ucvol,usewvl,xred,&
         &          mpi_atmtab,mpi_comm_atom,mpi_comm_fft,mpi_comm_wvl,me_g0,paral_kgb,distribfft) ! optional arguments

      USE_DEFS
      use m_distribfft,   only : distribfft_type

      use m_pawang,       only : pawang_type
      use m_pawtab,       only : pawtab_type
      use m_pawfgrtab,    only : pawfgrtab_type
      use m_pawrhoij,     only : pawrhoij_type

      implicit none

      !Arguments ---------------------------------------------
      !scalars
      integer,intent(in) :: cplex,ider,idir,ipert,izero,my_natom,natom,nfft
      integer,intent(in)  :: usewvl
      integer,intent(in) :: nhatgrdim,nspden,ntypat
      real(dp),intent(in) :: ucvol
      real(dp),intent(out) :: compch_fft
      type(pawang_type),intent(in) :: pawang

      integer,optional,intent(in) :: me_g0,mpi_comm_atom,mpi_comm_fft,mpi_comm_wvl,paral_kgb
      type(distribfft_type),optional,intent(in),target :: distribfft
      !arrays
      integer,intent(in) :: ngfft(18)
      real(dp),intent(in) :: gprimd(3,3),qphon(3),rprimd(3,3),xred(3,natom)
      real(dp),intent(out) :: pawgrnhat(cplex*nfft,nspden,3*nhatgrdim)
      real(dp),intent(inout) :: pawnhat(cplex*nfft,nspden) !vz_i
      type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
      type(pawrhoij_type),intent(in) :: pawrhoij(my_natom),pawrhoij0(my_natom)
      type(pawtab_type),intent(in) :: pawtab(ntypat)

      integer,optional,target,intent(in) :: mpi_atmtab(:)
    end subroutine pawmknhat
 end interface

 interface
    subroutine pawaccrhoij(atindx,cplex,cwaveprj,cwaveprj1,ipert,isppol,my_natom,natom,&
         &                       nspinor,occ_k,option,pawrhoij,usetimerev,wtk_k, &
         &                       occ_k_2, &
         &                       mpi_comm_atom,mpi_atmtab ) ! optional (parallelism)

      USE_DEFS

      use m_pawrhoij,   only : pawrhoij_type
      use m_pawcprj,    only : pawcprj_type

      implicit none

      !Arguments ---------------------------------------------
      !scalars
      integer,intent(in) :: cplex,ipert,isppol,my_natom,natom,nspinor,option
      integer,optional,intent(in) :: mpi_comm_atom
      logical,intent(in) :: usetimerev
      real(dp),intent(in) :: occ_k,wtk_k
      real(dp),optional,intent(in) :: occ_k_2
      !arrays
      integer,intent(in) :: atindx(natom)
      integer,optional,target,intent(in) :: mpi_atmtab(:)
      type(pawcprj_type),intent(in) :: cwaveprj(natom,nspinor),cwaveprj1(natom,nspinor)
      type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
    end subroutine pawaccrhoij
 end interface

 interface
    subroutine pawmkrho(compch_fft,cplex,gprimd,idir,indsym,ipert,&
         &          nfft, nfftc, ngfft, &
         &          my_natom,natom,nspden,nsym,ntypat,paral_kgb,pawang,pawfgrtab,pawprtvol,&
         &          pawrhoij,pawrhoij_unsym,&
         &          pawtab,qphon,rhopsg,rhopsr,rhor,rprimd,symafm,symrec,typat,ucvol,usewvl,xred,&
         &          comm_fft, comm_wvl, comm_atom, mpi_atmtab, pawang_sym,pawnhat,pawrhoij0) ! optional arguments

      USE_DEFS
      use m_pawang,   only : pawang_type
      use m_pawtab,   only : pawtab_type
      use m_pawfgrtab,only : pawfgrtab_type
      use m_pawrhoij, only : pawrhoij_type

      implicit none

      !Arguments ------------------------------------
      !scalars
      integer,intent(in) :: cplex,idir,ipert,my_natom,natom,nspden,nsym,ntypat,paral_kgb,pawprtvol
      integer, intent(in) :: nfft, nfftc, ngfft(18)
      integer,intent(in) :: usewvl
      real(dp),intent(in) :: ucvol
      real(dp),intent(out) :: compch_fft
      type(pawang_type),intent(in) :: pawang
      type(pawang_type),intent(in),optional :: pawang_sym
      integer,optional,intent(in) :: comm_atom, comm_fft, comm_wvl
      !arrays
      integer,intent(in) :: indsym(4,nsym,natom)
      integer,intent(in) :: symafm(nsym),symrec(3,3,nsym),typat(natom)
      real(dp),intent(in) :: gprimd(3,3),qphon(3),rprimd(3,3),xred(3,natom)
      real(dp),intent(inout),target,optional :: pawnhat(cplex*nfft,nspden) !vz_i
      real(dp),intent(inout) :: rhor(cplex*nfft,nspden)
      real(dp),intent(inout) :: rhopsg(2,nfftc),rhopsr(cplex*nfftc,nspden)
      type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
      type(pawrhoij_type),intent(inout),target :: pawrhoij(:)
      type(pawrhoij_type),intent(inout) :: pawrhoij_unsym(:)
      type(pawrhoij_type),intent(in),target,optional :: pawrhoij0(my_natom)
      type(pawtab_type),intent(in) :: pawtab(ntypat)
      integer,optional,target,intent(in) :: mpi_atmtab(:)
    end subroutine pawmkrho
 end interface

 interface
    subroutine pawdensities(compch_sph,cplex,iatom,lmselectin,lmselectout,lm_size,nhat1,nspden,nzlmopt,&
         &          opt_compch,opt_dens,opt_l,opt_print,pawang,pawprtvol,pawrad,pawrhoij,pawtab,rho1,trho1,&
         &          one_over_rad2) ! optional

      USE_DEFS

      use m_pawang, only : pawang_type
      use m_pawrad, only : pawrad_type
      use m_pawtab, only : pawtab_type
      use m_pawrhoij, only : pawrhoij_type

      implicit none

      !Arguments ---------------------------------------------
      !scalars
      integer,intent(in) :: cplex,iatom,lm_size,nspden,nzlmopt,opt_compch,opt_dens,opt_l,opt_print,pawprtvol
      ! jmb  real(dp),intent(out) :: compch_sph
      real(dp),intent(inout) :: compch_sph
      type(pawang_type),intent(in) :: pawang
      type(pawrad_type),intent(in) :: pawrad
      type(pawrhoij_type),intent(in) :: pawrhoij
      type(pawtab_type),intent(in) :: pawtab
      !arrays
      logical,intent(in) :: lmselectin(lm_size)
      logical,intent(inout) :: lmselectout(lm_size) !vz_i
      real(dp),intent(in),target,optional :: one_over_rad2(pawrad%mesh_size)
      real(dp),intent(out) :: nhat1(cplex*pawrad%mesh_size,lm_size,nspden*(1-((opt_dens+1)/2)))
      real(dp),intent(out) ::  rho1(cplex*pawrad%mesh_size,lm_size,nspden)
      real(dp),intent(out) :: trho1(cplex*pawrad%mesh_size,lm_size,nspden*(1-(opt_dens/2)))
    end subroutine pawdensities
 end interface

 interface
    subroutine pawdenpot(compch_sph,epaw,epawdc,ipert,ixc,&
         & my_natom,natom,nspden,ntypat,nzlmopt,option,paw_an,paw_an0,&
         & paw_ij,pawang,pawprtvol,pawrad,pawrhoij,pawspnorb,pawtab,pawxcdev,spnorbscl,xclevel,xc_denpos,ucvol,znucl,&
         & mpi_atmtab,mpi_comm_atom,vpotzero) ! optional arguments

      USE_DEFS
      use m_pawang,  only: pawang_type
      use m_pawrad,  only: pawrad_type
      use m_pawtab,  only: pawtab_type
      use m_paw_an, only : paw_an_type
      use m_paw_ij, only : paw_ij_type
      use m_pawrhoij,only: pawrhoij_type
      ! use m_electronpositron, only: electronpositron_type,electronpositron_calctype

      implicit none

      !Arguments ---------------------------------------------
      !scalars
      integer,intent(in) :: ipert,ixc,my_natom,natom,nspden,ntypat,nzlmopt,option,pawprtvol
      integer,intent(in) :: pawspnorb,pawxcdev,xclevel
      integer,optional,intent(in) :: mpi_comm_atom 
      real(dp), intent(in) :: spnorbscl,xc_denpos,ucvol
      real(dp),intent(out) :: compch_sph,epaw,epawdc
      ! xg  real(dp),intent(inout) :: compch_sph,epaw,epawdc
      real(dp),intent(out),optional :: vpotzero(2)
      ! type(electronpositron_type),pointer,optional :: electronpositron
      type(pawang_type),intent(in) :: pawang
      !arrays
      integer,optional,target,intent(in) :: mpi_atmtab(:)
      real(dp) :: znucl(ntypat)
      type(paw_an_type),intent(inout) :: paw_an(my_natom)
      type(paw_an_type), intent(in) :: paw_an0(my_natom)
      type(paw_ij_type),intent(inout) :: paw_ij(my_natom)
      type(pawrad_type),intent(in) :: pawrad(ntypat)
      type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
      type(pawtab_type),intent(in) :: pawtab(ntypat)
    end subroutine pawdenpot
 end interface

 interface
    subroutine wvl_nhatgrid(atindx1,geocode,h,i3s,natom,natom_tot,&
         & nattyp,ntypat,n1,n1i,n2,n2i,n3,n3pi,optcut,optgr0,optgr1,optgr2,optrad,&
         & pawfgrtab,pawtab,shift,rxyz)

      USE_DEFS

      use m_pawtab,       only : pawtab_type
      use m_pawfgrtab,    only : pawfgrtab_type

      implicit none

      !Arguments ---------------------------------------------
      !scalars
      integer,intent(in) :: i3s,natom,natom_tot,ntypat,optcut,optgr0,optgr1,optgr2,optrad
      integer,intent(in) :: n1,n2,n3,n1i,n2i,n3pi,shift
      real(dp),intent(in) :: h(3)
      character(1),intent(in) :: geocode
      !integer,intent(in),optional :: mpi_comm_wvl
      !arrays
      integer,intent(in) :: atindx1(natom),nattyp(ntypat)
      real(dp),intent(in) :: rxyz(3,natom)
      type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
      type(pawtab_type),intent(in) :: pawtab(ntypat)
    end subroutine wvl_nhatgrid
 end interface
end module abi_interfaces_libpaw
!!***
