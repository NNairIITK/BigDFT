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

end module abi_interfaces_libpaw
!!***
