!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paral_atom
!! NAME
!!  m_paral_atom
!!
!! FUNCTION
!!  This module provides routines and methods used to manage the parallelisation/distribution
!!  of PAW data over atomic sites
!!
!! COPYRIGHT
!! Copyright (C) 2012-2014 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paral_atom

 use defs_basis
 use m_profiling
 use m_errors
 use m_xmpi

 implicit none

 private

!public procedures.
 public :: get_my_natom
 public :: get_my_atmtab
 public :: free_my_atmtab
 public :: get_proc_atmtab
!!***

CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* m_paral_atom/get_my_natom
!! NAME
!! get_my_natom
!!
!! FUNCTION
!! Given the total number of atoms, return the number of atoms treated by current process
!!
!! COPYRIGHT
!! Copyright (C) 2012-2014 ABINIT group (MD,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  comm_atom=communicator over atoms
!!  natom=total number of atoms
!!
!! OUTPUT
!!  my_natom=number of atoms treated by current process
!!
!! PARENTS
!!      initmpi_atom,m_paw_an,m_paw_ij,m_pawfgrtab,m_pawrhoij
!!
!! CHILDREN
!!
!! SOURCE


subroutine get_my_natom(comm_atom,my_natom,natom)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_my_natom'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: comm_atom,natom
 integer,intent(out) :: my_natom
!arrays

!Local variables ---------------------------------------
!scalars
 integer ::  me,nproc
!arrays

! *************************************************************************

 my_natom=natom
 if (xmpi_paral==1) then
   if (comm_atom/=xmpi_self.and.comm_atom/=xmpi_comm_null)  then
     nproc=xcomm_size(comm_atom)
     me=xcomm_rank(comm_atom)
     my_natom=natom/nproc
     if (me<=(mod(natom,nproc)-1)) my_natom=natom/nproc + 1
   endif
 endif

end subroutine get_my_natom
!!***

!----------------------------------------------------------------------

!!****f* m_paral_atom/get_my_atmtab
!! NAME
!! get_my_atmtab
!!
!! FUNCTION
!! Given the total number of atoms and a MPI communicator return a table
!! containing the indexes of atoms treated by current processor.
!!
!! COPYRIGHT
!! Copyright (C) 2012-2014 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  comm_atom=communicator over atoms
!!  my_natom_ref= --optional-- a reference value for the local number of atoms
!!                (just for checking purposes)
!!  natom=total number of atoms
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  my_atmtab(:)=indexes of atoms treated by current process
!!               nothing is done if my_atmtab(:) is already associated to a target
!!  my_atmtab_allocated=true if my_atmtab is allocated
!!  paral_atom=flag controlling parallelism over atoms
!!
!! PARENTS
!!      accrho3,calc_efg,calc_fc,denfgr,eltxccore,ewald3,ewald4,expibi,expibr
!!      initmpi_atom,initrhoij,m_hamiltonian,m_paw_an,m_paw_ij,m_paw_pwaves_lmn
!!      m_pawdij,m_pawfgrtab,m_pawrhoij,make_efg_onsite,make_fc_paw,newfermie1
!!      nhatgrid,paw_mknewh0,pawaccrhoij,pawdenpot,pawenergy3,pawgrnl,pawmkaewf
!!      pawmknhat,pawmknhat_psipsi,pawnhatfr,pawprt,pawsushat,pawtwdij
!!      pawtwdij_1,pawtwdij_2a,pawtwdij_2b,pawtwdij_2c,pawtwdij_2d,pawtwdij_2e
!!      pawtwdij_2f,pawuj_red,qijb_kk,setnoccmmp,setrhoijpbe0,twexpibi
!!      twqijb_kk
!!
!! CHILDREN
!!
!! SOURCE


subroutine get_my_atmtab(comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,&
&                        my_natom_ref) ! optional argument


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_my_atmtab'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: comm_atom,natom
 integer,intent(in),optional :: my_natom_ref
 logical,intent(inout) :: my_atmtab_allocated,paral_atom
!arrays
 integer,pointer :: my_atmtab(:)
!Local variables ---------------------------------------
!scalars
 integer :: iatom,me,my_natom,natom_bef,nmod,nproc
 character(len=500) :: msg
!arrays

! *************************************************************************

 my_atmtab_allocated=.false.
 if (.not.paral_atom) return

 if (comm_atom==xmpi_self.or.comm_atom==xmpi_comm_null) paral_atom=.false.
 if (paral_atom)  then
   nproc=xcomm_size(comm_atom)
   paral_atom=(nproc>1)
   if (paral_atom) then
     if (.not.associated(my_atmtab)) then
!      Get local number of atoms
       me=xcomm_rank(comm_atom)
       my_natom=natom/nproc
       if (me<=(mod(natom,nproc)-1)) my_natom=natom/nproc + 1
!      Get table of atoms
       if (my_natom>0) then
         ABI_ALLOCATE(my_atmtab,(my_natom))
         my_atmtab_allocated=.true.
         if (my_natom==natom) then
           my_atmtab(1:my_natom)=(/(iatom,iatom=1,natom)/)
         else
!          The atoms are distributed contigously by egal part
!          (the rest is distributed on all the procs)
           nmod=mod(natom,nproc)
           if (me<=(nmod-1)) then
             natom_bef=me*(natom/nproc)+me
           else
             natom_bef=me*(natom/nproc)+nmod
           endif
           do iatom=1,my_natom
             my_atmtab(iatom)=iatom+natom_bef
           enddo
         end if
       end if
     else
       my_natom=size(my_atmtab)
     end if
     if (present(my_natom_ref).and.(my_natom>0)) then
       if (my_natom_ref/=size(my_atmtab)) then
         msg='my_atmtab should have a size equal to my_natom !'
         MSG_BUG(msg)
       end if
     end if
   end if
 endif

end subroutine get_my_atmtab
!!***

!----------------------------------------------------------------------

!!****f* m_paral_atom/free_my_atmtab
!! NAME
!! free_my_atmtab
!!
!! FUNCTION
!! Cleanly deallocate a table of atom indexes (my_atmtab)
!!
!! COPYRIGHT
!! Copyright (C) 2012-2014 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  my_atmtab_allocated=true if my_atmtab is allocated
!!  my_atmtab(:)=indexes of atoms treated by current process
!!               nothing is done if my_atmtab(:) is already associated to a target
!!
!! PARENTS
!!      calc_efg,calc_fc,denfgr,eltxccore,ewald3,ewald4,expibi,expibr,initrhoij
!!      m_hamiltonian,m_paw_an,m_paw_ij,m_paw_pwaves_lmn,m_pawdij,m_pawfgrtab
!!      m_pawrhoij,make_efg_onsite,make_fc_paw,newfermie1,nhatgrid,paw_mknewh0
!!      pawaccrhoij,pawdenpot,pawenergy3,pawgrnl,pawmkaewf,pawmknhat
!!      pawmknhat_psipsi,pawnhatfr,pawprt,pawsushat,pawtwdij,pawtwdij_1
!!      pawtwdij_2a,pawtwdij_2b,pawtwdij_2c,pawtwdij_2d,pawtwdij_2e,pawtwdij_2f
!!      pawuj_red,qijb_kk,setnoccmmp,setrhoijpbe0,twexpibi,twqijb_kk
!!
!! CHILDREN
!!
!! SOURCE

subroutine free_my_atmtab(my_atmtab,my_atmtab_allocated)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'free_my_atmtab'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 logical,intent(inout) :: my_atmtab_allocated
!arrays
 integer,pointer :: my_atmtab(:)
!Local variables ---------------------------------------
!scalars
!arrays

! *************************************************************************

 if (my_atmtab_allocated) then
   ABI_DEALLOCATE(my_atmtab)
   nullify(my_atmtab)
   my_atmtab_allocated=.false.
 end if

end subroutine free_my_atmtab
!!***

!----------------------------------------------------------------------

!!****f* m_paral_atom/get_proc_atmtab
!! NAME
!! get_proc_atmtab
!!
!! FUNCTION
!!  Given the total number of atoms and the size of a communicator,
!!  return a table containing the indexes of atoms treated by a processor (with index iproc)
!!
!! COPYRIGHT
!! Copyright (C) 2012-2013 ABINIT group (MD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  comm_atom_size= size of communicator (over atoms)
!!  iproc= rank of the processor
!!  natom= total number of atoms
!!
!! OUTPUT
!!  atmtab(natom_out)= indexes of atoms treated by process iproc
!!  natom_out= number of atoms treated by process iproc
!!
!! NOTES
!! In case of modification of the distribution of atom over proc,
!!   get_atmtab must be modified accordingly
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine get_proc_atmtab(iproc,atmtab,natom_out,natom,comm_atom_size)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_proc_atmtab'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer, intent(in) :: comm_atom_size,iproc,natom
!arrays
 integer,intent(out) :: natom_out
 integer,allocatable, intent(out):: atmtab(:)

!Local variables ---------------------------------------
!scalars
 integer :: iatom,natom_bef,nmod,nproc
!arrays

! *************************************************************************

 nproc=comm_atom_size

 natom_out=natom/nproc ; if (iproc<=(mod(natom,nproc)-1)) natom_out=natom/nproc+1

! Get table of atoms
 if (natom_out>0) then
   ABI_ALLOCATE(atmtab,(natom_out))
!  The atoms are distributed contigously by egal part
!  The rest is distributed on all the procs
!  (see get_my_atmtab)
   nmod=mod(natom,nproc)
   if (iproc<=(nmod-1)) then
     natom_bef=iproc*(natom/nproc)+iproc
   else
     natom_bef=iproc*(natom/nproc)+nmod
   end if
   do iatom=1,natom_out
     atmtab(iatom)=iatom+natom_bef
   end do

 else
   natom_out=0
   ABI_ALLOCATE(atmtab,(0))
 end if

end subroutine get_proc_atmtab
!!***

!----------------------------------------------------------------------

END MODULE m_paral_atom
!!***
