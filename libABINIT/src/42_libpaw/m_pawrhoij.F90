!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_pawrhoij
!! NAME
!!  m_pawrhoij
!!
!! FUNCTION
!!  This module contains the definition of the pawrhoij_type structured datatype,
!!  as well as related functions and methods.
!!  pawrhoij_type variables define rhoij occupancies matrixes used within PAW formalism.
!!
!! COPYRIGHT
!! Copyright (C) 2012-2013 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

#include "abi_common_for_bigdft.h"

MODULE m_pawrhoij

 use defs_basis
 use m_pawtab, only : pawtab_type
! use m_errors
use interfaces_12_hide_mpi
use interfaces_14_hidewrite
use interfaces_16_hideleave
 use m_profiling
 use m_xmpi
 use m_paral_atom

 implicit none

 private

!public procedures.
 public :: rhoij_alloc
 public :: rhoij_free
 public :: rhoij_nullify
 public :: rhoij_copy
 public :: rhoij_gather
 public :: rhoij_bcast
 public :: rhoij_io
 public :: rhoij_unpack
 public :: rhoij_init_unpacked
 public :: rhoij_destroy_unpacked
 public :: rhoij_mpi_sum
!!***

!!****t* m_pawrhoij/pawrhoij_type
!! NAME
!! pawrhoij_type
!!
!! FUNCTION
!! This structured datatype contains rhoij quantities (occucpancies)
!! and related data, used in PAW calculations.
!!
!! SOURCE

 type,public :: pawrhoij_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

!Integer scalars

  integer :: cplex
   ! cplex=1 if rhoij are real, 2 if rhoij are complex
   ! For GS calculations: rhoij are real provided that time-reversal can be used.

  integer :: itypat
   ! itypat=type of the atom

  integer :: lmn_size
   ! Number of (l,m,n) elements for the paw basis

  integer :: lmn2_size
   ! lmn2_size=lmn_size*(lmn_size+1)/2
   ! where lmn_size is the number of (l,m,n) elements for the paw basis

  integer :: lmnmix_sz
   ! lmnmix_sz=number of (lmn,lmn_prime) verifying l<=lmix and l_prime<=lmix
   !           i.e. number of rhoij elements being mixed during SCF cycle
   ! lmnmix_sz=0 if mixing data are not used

  integer :: ngrhoij
   ! First dimension of array grhoij

  integer :: nrhoijsel
   ! nrhoijsel
   ! Number of non-zero value of rhoij
   ! This is the size of rhoijp(:,:) (see below in this datastructure)

  integer :: nspden
   ! Number of spin-density components for rhoij (may be different from nspden for density)

  integer :: nspinor
   ! Number of spinorial components

  integer :: nsppol
   ! Number of independent spin-components

  integer :: use_rhoij_
   ! 1 if pawrhoij%rhoij_ is allocated

  integer :: use_rhoijp
   ! 1 if pawrhoij%rhoijp and pawrhoij%rhoijselect are allocated

  integer :: use_rhoijres
   ! 1 if pawrhoij%rhoijres is allocated

!Integer arrays

  integer, pointer :: kpawmix(:)
   ! kpawmix(lmnmix_sz)
   ! Indirect array selecting the elements of rhoij
   ! being mixed during SCF cycle

  integer, pointer :: rhoijselect(:)
   ! rhoijselect(lmn2_size)
   ! Indirect array selecting the non-zero elements of rhoij:
   ! rhoijselect(isel,ispden)=klmn if rhoij(klmn,ispden) is non-zero

!Real (real(dp)) arrays

  real(dp), pointer :: grhoij (:,:,:)
   ! grhoij(ngrhoij,cplex*lmn2_size,nspden)
   ! Gradients of Rho_ij wrt xred, strains, ... (non-packed storage)

  real(dp), pointer :: rhoij_ (:,:)
   ! rhoij_(cplex*lmn2_size,nspden)
   ! Array used to (temporary) store Rho_ij in a non-packed storage mode

  real(dp), pointer :: rhoijp (:,:)
   ! rhoijp(cplex*lmn2_size,nspden)
   ! Augmentation waves occupancies Rho_ij in PACKED STORAGE (only non-zero elements are stored)

  real(dp), pointer :: rhoijres (:,:)
   ! rhoijres(cplex*lmn2_size,nspden)
   ! Rho_ij residuals during SCF cycle (non-packed storage)

 end type pawrhoij_type
!!***

CONTAINS

!===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/rhoij_alloc
!! NAME
!! rhoij_alloc
!!
!! FUNCTION
!! Initialize and allocate a pawrhoij datastructure
!!
!! INPUTS
!! [comm_atom] = communicator over atoms  (OPTIONAL)
!! cplex=1 if rhoij is REAL,2 if COMPLEX
!! [my_atmtab(:)] = Index of atoms treated by current proc (OPTIONAL)
!! nspden=number of spin-components for rhoij
!! nsppol=number of spinorial components for rhoij
!! nsppol=number of independant spin-components for rhoij
!! typat(:)=types of atoms
!! [lmnsize(:)]=array of (l,m,n) sizes for rhoij for each type of atom (OPTIONAL)
!!              must be present if [pawtab] argument is not passed
!! [ngrhoij]=number of gradients to be allocated (OPTIONAL, default=0)
!! [nlmnmix]=number of rhoij elements to be mixed during SCF cycle (OPTIONAL, default=0)
!! [pawtab(:)] <type(pawtab_type)>=paw tabulated starting data (OPTIONAL)
!!             must be present if [lmnsize(:)] argument is not passed
!! [use_rhoij_]=1 if pawrhoij(:)%rhoij_ has to be allocated (OPTIONAL, default=0)
!! [use_rhoijp]=1 if pawrhoij(:)%rhoijp has to be allocated (OPTIONAL, default=1)
!!              (in that case, pawrhoij%rhoijselect is also allocated)
!! [use_rhoijres]=1 if pawrhoij(:)%rhoijres has to be allocated (OPTIONAL, default=0)
!!
!! SIDE EFFECTS
!! pawrhoij(:)<type(pawrhoij_type)>= rhoij datastructure
!!
!! NOTES
!! One of the two optional arguments lmnsize(:) or pawtab(:) must be present !
!! If both are present, only pawtab(:) is used.
!!
!! PARENTS
!!      bethe_salpeter,energy,extraprho,initrhoij,loper3,m_electronpositron
!!      m_header,m_pawrhoij,m_qparticles,nstpaw3,paw_qpscgw,respfn,scfcv3
!!      screening,setup_bse,setup_positron,setup_screening,setup_sigma,sigma
!!      vtorho
!!
!! CHILDREN
!!      xsum_mpi
!!
!! SOURCE

subroutine rhoij_alloc(pawrhoij,cplex,nspden,nspinor,nsppol,typat,&           ! Mandatory arguments
&                      lmnsize,ngrhoij,nlmnmix,pawtab,use_rhoij_,use_rhoijp,& ! Optional arguments
&                      use_rhoijres,mpi_comm_atom,mpi_atmtab)                 ! Optional arguments


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_alloc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nspden,nspinor,nsppol
 integer,optional,intent(in):: mpi_comm_atom,ngrhoij,nlmnmix,use_rhoij_,use_rhoijp,use_rhoijres
 integer,optional,target,intent(in)  :: mpi_atmtab(:)
!arrays
 integer,intent(in) :: typat(:)
 integer,optional,target,intent(in) :: lmnsize(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)
 type(pawtab_type),optional,intent(in)  :: pawtab(:)

!Local variables-------------------------------
!scalars
 integer :: irhoij,irhoij_,itypat,lmn2_size,nn1,natom,nrhoij
 logical :: has_rhoijp,my_atmtab_allocated,paral_atom
 character(len=500) :: msg
!array
 integer,pointer :: lmn_size(:),my_atmtab(:)

! *************************************************************************

 nrhoij=size(pawrhoij);natom=size(typat)
 if (nrhoij>natom) then
   msg=' wrong sizes (1) !'
   call wrtout(std_out,msg,'COLL')
   call leave_new('COLL')
 end if

!Select lmn_size for each atom type
 if (present(pawtab)) then
   nn1=size(pawtab)
   if (maxval(typat)>nn1) then
     msg=' wrong sizes (2) !'
     call wrtout(std_out,msg,'COLL')
     call leave_new('COLL')
   end if
   ABI_ALLOCATE(lmn_size,(nn1))
   do itypat=1,nn1
     lmn_size(itypat)=pawtab(itypat)%lmn_size
   end do
 else if (present(lmnsize)) then
   nn1=size(lmnsize)
   if (maxval(typat)>nn1) then
     msg=' wrong sizes (3) !'
     call wrtout(std_out,msg,'COLL')
     call leave_new('COLL')
   end if
   lmn_size => lmnsize
 else
   msg=' one of the 2 arguments pawtab or lmnsize must be present !'
   call wrtout(std_out,msg,'COLL')
   call leave_new('COLL')
 end if

!Set up parallelism over atoms
 paral_atom=(present(mpi_comm_atom).and.(nrhoij/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 call get_my_atmtab(mpi_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom)

 if (nrhoij>0) then
   do irhoij=1,nrhoij
     irhoij_=irhoij;if (paral_atom) irhoij_=my_atmtab(irhoij)
     itypat=typat(irhoij_)

     lmn2_size=lmn_size(itypat)*(lmn_size(itypat)+1)/2

!    Scalars initializations
     pawrhoij(irhoij)%cplex=cplex
     pawrhoij(irhoij)%itypat=itypat
     pawrhoij(irhoij)%lmn_size=lmn_size(itypat)
     pawrhoij(irhoij)%lmn2_size=lmn2_size
     pawrhoij(irhoij)%nspden=nspden
     pawrhoij(irhoij)%nspinor=nspinor
     pawrhoij(irhoij)%nsppol=nsppol
     pawrhoij(irhoij)%nrhoijsel=0
     pawrhoij(irhoij)%lmnmix_sz=0
     pawrhoij(irhoij)%ngrhoij=0
     pawrhoij(irhoij)%use_rhoij_=0
     pawrhoij(irhoij)%use_rhoijres=0
     nullify(pawrhoij(irhoij)%kpawmix)
     nullify(pawrhoij(irhoij)%grhoij)
     nullify(pawrhoij(irhoij)%rhoij_)
     nullify(pawrhoij(irhoij)%rhoijp)
     nullify(pawrhoij(irhoij)%rhoijres)
     nullify(pawrhoij(irhoij)%rhoijselect)

!    Arrays allocations
     has_rhoijp=.true.; if (present(use_rhoijp)) has_rhoijp=(use_rhoijp>0)
     if (has_rhoijp) then
       pawrhoij(irhoij)%use_rhoijp=1
       ABI_ALLOCATE(pawrhoij(irhoij)%rhoijselect,(lmn2_size))
       ABI_ALLOCATE(pawrhoij(irhoij)%rhoijp,(cplex*lmn2_size,nspden))
       pawrhoij(irhoij)%rhoijselect(:)=0
       pawrhoij(irhoij)%rhoijp(:,:)=zero
     end if

     if (present(ngrhoij)) then
       if (ngrhoij>0) then
         pawrhoij(irhoij)%ngrhoij=ngrhoij
         ABI_ALLOCATE(pawrhoij(irhoij)%grhoij,(ngrhoij,cplex*lmn2_size,nspden))
         pawrhoij(irhoij)%grhoij=zero
       end if
     end if
     if (present(nlmnmix)) then
       if (nlmnmix>0) then
         pawrhoij(irhoij)%lmnmix_sz=nlmnmix
         ABI_ALLOCATE(pawrhoij(irhoij)%kpawmix,(nlmnmix))
         pawrhoij(irhoij)%kpawmix=0
       end if
     end if
     if (present(use_rhoij_)) then
       if (use_rhoij_>0) then
         pawrhoij(irhoij)%use_rhoij_=use_rhoij_
         ABI_ALLOCATE(pawrhoij(irhoij)%rhoij_,(cplex*lmn2_size,nspden))
         pawrhoij(irhoij)%rhoij_=zero
       end if
     end if
     if (present(use_rhoijres)) then
       if (use_rhoijres>0) then
         pawrhoij(irhoij)%use_rhoijres=use_rhoijres
         ABI_ALLOCATE(pawrhoij(irhoij)%rhoijres,(cplex*lmn2_size,nspden))
         pawrhoij(irhoij)%rhoijres=zero
       end if
     end if

   end do
 end if

 if (present(pawtab)) then
   ABI_DEALLOCATE(lmn_size)
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine rhoij_alloc
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/rhoij_free
!! NAME
!! rhoij_free
!!
!! FUNCTION
!! Destroy a pawrhoij datastructure
!!
!! SIDE EFFECTS
!! pawrhoij(:)<type(pawrhoij_type)>= rhoij datastructure
!!
!! PARENTS
!!      bethe_salpeter,energy,gstate,gw_tools,loper3,m_electronpositron
!!      m_header,m_pawrhoij,m_scf_history,nstpaw3,pawgrnl,pawmkrhoij,pawprt
!!      respfn,scfcv3,screening,setup_bse,setup_positron,setup_screening
!!      setup_sigma,sigma,symrhoij,vtorho
!!
!! CHILDREN
!!      xsum_mpi
!!
!! SOURCE

subroutine rhoij_free(pawrhoij)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables-------------------------------
!scalars
 integer :: irhoij,nrhoij

! *************************************************************************

 nrhoij=size(pawrhoij)

 if (nrhoij>0) then
   do irhoij=1,nrhoij
     pawrhoij(irhoij)%nrhoijsel=0
     pawrhoij(irhoij)%ngrhoij=0
     pawrhoij(irhoij)%lmnmix_sz=0
     pawrhoij(irhoij)%use_rhoij_=0
     pawrhoij(irhoij)%use_rhoijp=0
     pawrhoij(irhoij)%use_rhoijres=0
     if (associated(pawrhoij(irhoij)%rhoijp))       then
       ABI_DEALLOCATE(pawrhoij(irhoij)%rhoijp)
       nullify(pawrhoij(irhoij)%rhoijp)
     end if
     if (associated(pawrhoij(irhoij)%rhoijselect))  then
       ABI_DEALLOCATE(pawrhoij(irhoij)%rhoijselect)
       nullify(pawrhoij(irhoij)%rhoijselect)
     end if
     if (associated(pawrhoij(irhoij)%grhoij))       then
       ABI_DEALLOCATE(pawrhoij(irhoij)%grhoij)
       nullify(pawrhoij(irhoij)%grhoij)
     end if
     if (associated(pawrhoij(irhoij)%kpawmix))      then
       ABI_DEALLOCATE(pawrhoij(irhoij)%kpawmix)
       nullify(pawrhoij(irhoij)%kpawmix)
     end if
     if (associated(pawrhoij(irhoij)%rhoij_))       then
       ABI_DEALLOCATE(pawrhoij(irhoij)%rhoij_)
       nullify(pawrhoij(irhoij)%rhoij_)
     end if
     if (associated(pawrhoij(irhoij)%rhoijres))     then
       ABI_DEALLOCATE(pawrhoij(irhoij)%rhoijres)
       nullify(pawrhoij(irhoij)%rhoijres)
     end if
   end do
 end if

end subroutine rhoij_free
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/rhoij_nullify
!! NAME
!! rhoij_nullify
!!
!! FUNCTION
!! Nullify (initialize to null) a pawrhoij datastructure
!!
!! SIDE EFFECTS
!! pawrhoij(:)<type(pawrhoij_type)>= rhoij datastructure
!!
!! PARENTS
!!      m_scf_history,pawgrnl,pawprt,symrhoij
!!
!! CHILDREN
!!      xsum_mpi
!!
!! SOURCE

subroutine rhoij_nullify(pawrhoij)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables-------------------------------
!scalars
 integer :: irhoij,nrhoij

! *************************************************************************

 nrhoij=size(pawrhoij)

 if (nrhoij>0) then
   do irhoij=1,nrhoij
     pawrhoij(irhoij)%nrhoijsel=0
     pawrhoij(irhoij)%ngrhoij=0
     pawrhoij(irhoij)%lmnmix_sz=0
     pawrhoij(irhoij)%use_rhoij_=0
     pawrhoij(irhoij)%use_rhoijp=0
     pawrhoij(irhoij)%use_rhoijres=0
     nullify(pawrhoij(irhoij)%kpawmix)
     nullify(pawrhoij(irhoij)%rhoijp)
     nullify(pawrhoij(irhoij)%rhoijselect)
     nullify(pawrhoij(irhoij)%grhoij)
     nullify(pawrhoij(irhoij)%rhoij_)
     nullify(pawrhoij(irhoij)%rhoijres)

   end do
 end if

end subroutine rhoij_nullify
!!***

!!****f* m_pawrhoij/rhoij_copy
!! NAME
!! rhoij_copy
!!
!! FUNCTION
!! Copy one pawrhoij datastructure into another
!! Can take into accound changes of dimensions
!! Can copy a shared pawrhoij into distributed ones (when parallelism is activated)
!!
!! INPUTS
!!  keep_cplex= optional argument (logical, default=.TRUE.)
!!              if .TRUE. pawrhoij_out(:)%cplex is NOT MODIFIED,
!!              even if different from pawrhoij_in(:)%cplex
!!  keep_itypat= optional argument (logical, default=.FALSE.)
!!               if .TRUE. pawrhoij_out(:)%ityp is NOT MODIFIED,
!!               even if different from pawrhoij_in(:)%ityp
!!  keep_nspden= optional argument (logical, default=.TRUE.)
!!               if .TRUE. pawrhoij_out(:)%nspden is NOT MODIFIED,
!!               even if different from pawrhoij_in(:)%nspden
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  mpi_comm_atom=--optional-- MPI communicator over atoms
!!  pawrhoij_in(:)<type(pawrhoij_type)>= input rhoij datastructure
!!
!! SIDE EFFECTS
!!  pawrhoij_out(:)<type(pawrhoij_type)>= output rhoij datastructure
!!
!! PARENTS
!!      bethe_salpeter,gstate,inwffil,ioarr,loper3,m_electronpositron,m_header
!!      m_pawrhoij,respfn,screening,setup_bse,setup_positron,setup_screening
!!      setup_sigma,sigma
!!
!! CHILDREN
!!      xsum_mpi
!!
!! SOURCE

subroutine rhoij_copy(pawrhoij_in,pawrhoij_copy, &
&          keep_cplex,keep_itypat,keep_nspden,& ! optional arguments
&          mpi_atmtab,mpi_comm_atom)            ! optional arguments (parallelism)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_copy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: mpi_comm_atom
 logical,intent(in),optional :: keep_cplex,keep_itypat,keep_nspden
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(pawrhoij_type),intent(in) :: pawrhoij_in(:)
 type(pawrhoij_type),intent(inout),target :: pawrhoij_copy(:)

!Local variables-------------------------------
!scalars
 integer :: cplex_in,cplex_out,dplex_in,dplex_out,i_in,i_out,ilmn
 integer :: irhoij,ispden,jrhoij,lmn2_size_out,lmnmix,my_nrhoij
 integer :: ngrhoij,nrhoij_in,nrhoij_max,nrhoij_out,nselect,nselect_out
 integer :: nspden_in,nspden_out,paral_case,use_rhoij_,use_rhoijp,use_rhoijres
 logical :: change_dim,keep_cplex_,keep_itypat_,keep_nspden_,my_atmtab_allocated,paral_atom
 character(len=500) :: msg
!arrays
 integer,pointer :: my_atmtab(:)
 integer,allocatable :: nlmn(:),typat(:)
 type(pawrhoij_type),pointer :: pawrhoij_out(:)

! *************************************************************************

!Retrieve sizes
 nrhoij_in=size(pawrhoij_in);nrhoij_out=size(pawrhoij_copy)

!Init flags
 keep_cplex_=.true.
 if (present(keep_cplex)) keep_cplex_=keep_cplex
 keep_itypat_=.false.
 if (present(keep_itypat)) keep_itypat_=keep_itypat
 keep_nspden_=.true.
 if (present(keep_nspden)) keep_nspden_=keep_nspden

!Set up parallelism over atoms
 paral_atom=(present(mpi_comm_atom));if (paral_atom) paral_atom=(xcomm_size(mpi_comm_atom)>1)
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_atmtab_allocated=.false.

!Determine in which case we are (parallelism, ...)
!No parallelism: a single copy operation
 paral_case=0;nrhoij_max=nrhoij_in
 pawrhoij_out => pawrhoij_copy
 if (paral_atom) then
   if (nrhoij_out<nrhoij_in) then ! Parallelism: the copy operation is a scatter
     call get_my_natom(mpi_comm_atom,my_nrhoij,nrhoij_in)
     if (my_nrhoij==nrhoij_out) then
       call get_my_atmtab(mpi_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,nrhoij_in)
       paral_case=1;nrhoij_max=nrhoij_out
       pawrhoij_out => pawrhoij_copy
     else
       msg=' nrhoij_out should be equal to my_natom !'
       call wrtout(std_out,msg,'COLL')
       call leave_new('COLL')
     end if
   else                            ! Parallelism: the copy operation is a gather
     call get_my_natom(mpi_comm_atom,my_nrhoij,nrhoij_out)
     if (my_nrhoij==nrhoij_in) then
       paral_case=2;nrhoij_max=nrhoij_in
       ABI_DATATYPE_ALLOCATE(pawrhoij_out,(nrhoij_in))
       if (nrhoij_in>0) then
         ABI_ALLOCATE(typat,(nrhoij_in))
         ABI_ALLOCATE(nlmn,(nrhoij_in))
         do irhoij=1,nrhoij_in
           typat(irhoij)=irhoij;nlmn(irhoij)=pawrhoij_in(irhoij)%lmn_size
         end do
         call rhoij_alloc(pawrhoij_out,pawrhoij_copy(1)%cplex,pawrhoij_copy(1)%nspden,&
&         pawrhoij_copy(1)%nspinor,pawrhoij_copy(1)%nsppol,typat,&
&         lmnsize=nlmn,ngrhoij=pawrhoij_copy(1)%ngrhoij,nlmnmix=pawrhoij_copy(1)%lmnmix_sz,&
&         use_rhoij_=pawrhoij_copy(1)%use_rhoij_,&
&         use_rhoijp=pawrhoij_copy(1)%use_rhoijp,&
&         use_rhoijres=pawrhoij_copy(1)%use_rhoijres)
         ABI_DEALLOCATE(typat)
         ABI_DEALLOCATE(nlmn)
       end if
     else
       msg=' nrhoij_in should be equal to my_natom !'
       call wrtout(std_out,msg,'COLL')
       call leave_new('COLL')
     end if
   end if
 end if

!Loop on rhoij components
 if (nrhoij_max>0) then
   do irhoij=1,nrhoij_max
     jrhoij=irhoij;if (paral_case==1) jrhoij=my_atmtab(irhoij)

     lmn2_size_out=pawrhoij_in(jrhoij)%lmn2_size
     cplex_in=pawrhoij_in(jrhoij)%cplex
     cplex_out=cplex_in;if(keep_cplex_)cplex_out=pawrhoij_out(irhoij)%cplex
     nspden_in=pawrhoij_in(jrhoij)%nspden
     nselect=pawrhoij_in(jrhoij)%nrhoijsel
     nselect_out=pawrhoij_out(irhoij)%nrhoijsel
     nspden_out=nspden_in;if(keep_nspden_)nspden_out=pawrhoij_out(irhoij)%nspden

     change_dim=(pawrhoij_out(irhoij)%cplex/=cplex_out.or. &
&     pawrhoij_out(irhoij)%lmn2_size/=lmn2_size_out.or. &
&     pawrhoij_out(irhoij)%nspden/=nspden_out.or. &
&     nselect/=nselect_out)
     dplex_in=cplex_in-1;dplex_out=cplex_out-1

!    Scalars
     pawrhoij_out(irhoij)%cplex=cplex_out+0
     pawrhoij_out(irhoij)%nspden=nspden_out+0
     pawrhoij_out(irhoij)%lmn2_size=lmn2_size_out+0
     pawrhoij_out(irhoij)%lmn_size=pawrhoij_in(jrhoij)%lmn_size+0
     if(.not.keep_itypat_) pawrhoij_out(irhoij)%itypat =pawrhoij_in(jrhoij)%itypat+0
     if(.not.keep_nspden_) pawrhoij_out(irhoij)%nsppol =pawrhoij_in(jrhoij)%nsppol+0
     if(.not.keep_nspden_) pawrhoij_out(irhoij)%nspinor=pawrhoij_in(jrhoij)%nspinor+0
     pawrhoij_out(irhoij)%nrhoijsel=nselect+0
!    if (pawrhoij_out(irhoij)%itypat/=pawrhoij_in(jrhoij)%itypat) then
!    write(unit=msg,fmt='(a,i3,a)') 'Type of atom ',jrhoij,' is different (dont copy it) !'
!    MSG_COMMENT(msg)
!    end if

!    Optional pointer: non-zero elements of rhoij
     use_rhoijp=pawrhoij_in(jrhoij)%use_rhoijp
     if (pawrhoij_out(irhoij)%use_rhoijp/=use_rhoijp) then
       if (pawrhoij_out(irhoij)%use_rhoijp>0)  then
         ABI_DEALLOCATE(pawrhoij_out(irhoij)%rhoijp)
         ABI_DEALLOCATE(pawrhoij_out(irhoij)%rhoijselect)
       end if
       if (use_rhoijp>0)  then
         ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoijp,(cplex_out*lmn2_size_out,nspden_out))
         ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoijselect,(lmn2_size_out))
       end if
       pawrhoij_out(irhoij)%use_rhoijp=use_rhoijp
     end if
     if (use_rhoijp>0) then
       if (change_dim) then
         if(associated(pawrhoij_out(irhoij)%rhoijp)) then
           ABI_DEALLOCATE(pawrhoij_out(irhoij)%rhoijp)
         end if
         ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoijp,(cplex_out*lmn2_size_out,nspden_out))
       end if
       if (cplex_out==cplex_in.and.nspden_out==nspden_in) then
         do ispden=1,nspden_out
           pawrhoij_out(irhoij)%rhoijp(1:cplex_out*nselect,ispden)=pawrhoij_in(jrhoij)%rhoijp(1:cplex_out*nselect,ispden)+zero
           if ( (nselect<lmn2_size_out) .and. &
&               (size(pawrhoij_out(irhoij)%rhoijp(:,ispden))==cplex_out*lmn2_size_out) ) then 
             pawrhoij_out(irhoij)%rhoijp(cplex_out*nselect+1:cplex_out*lmn2_size_out,ispden)=zero
           end if 
         end do
       else
         pawrhoij_out(irhoij)%rhoijp(:,:)=zero
         if (nspden_out==1) then
           if (nspden_in==2) then
             do ilmn=1,nselect
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijp(i_out,1)=pawrhoij_in(jrhoij)%rhoijp(i_in,1) &
&               +pawrhoij_in(jrhoij)%rhoijp(i_in,2)+zero
             end do
           else ! nspden_in==1 or nspden_in=4
             do ilmn=1,nselect
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijp(i_out,1)=pawrhoij_in(jrhoij)%rhoijp(i_in,1)+zero
             end do
           end if
         else if (nspden_out==2) then
           if (nspden_in==1) then
             do ilmn=1,nselect
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijp(i_out,1)=half*pawrhoij_in(jrhoij)%rhoijp(i_in,1)+zero
               pawrhoij_out(irhoij)%rhoijp(i_out,2)=pawrhoij_in(jrhoij)%rhoijp(i_out,1)
             end do
           else if (nspden_in==2) then
             do ilmn=1,nselect
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijp(i_out,1:2)=pawrhoij_in(jrhoij)%rhoijp(i_in,1:2)+zero
             end do
           else ! nspden_in==4
             do ilmn=1,nselect
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijp(i_out,1)=half*(pawrhoij_in(jrhoij)%rhoijp(i_in,1) &
&               +pawrhoij_in(jrhoij)%rhoijp(i_in,4))+zero
               pawrhoij_out(irhoij)%rhoijp(i_out,2)=half*(pawrhoij_in(jrhoij)%rhoijp(i_in,1) &
&               -pawrhoij_in(jrhoij)%rhoijp(i_in,4))+zero
             end do
           end if
         else if (nspden_out==4) then
           if (nspden_in==1) then
             do ilmn=1,nselect
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijp(i_out,1)=pawrhoij_in(jrhoij)%rhoijp(i_in,1)+zero
             end do
           else if (nspden_in==2) then
             do ilmn=1,nselect
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijp(i_out,1)=pawrhoij_in(jrhoij)%rhoijp(i_in,1) &
&               +pawrhoij_in(jrhoij)%rhoijp(i_in,2)+zero
               pawrhoij_out(irhoij)%rhoijp(i_out,4)=pawrhoij_in(jrhoij)%rhoijp(i_in,1) &
&               -pawrhoij_in(jrhoij)%rhoijp(i_in,2)+zero
             end do
           else ! nspden_in==4
             do ilmn=1,nselect
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijp(i_out,1:4)=pawrhoij_in(jrhoij)%rhoijp(i_in,1:4)+zero
             end do
           end if
         end if
       end if
     end if

!    Optional pointer: indexes for non-zero elements selection
     if (use_rhoijp>0) then
       if (change_dim) then
         if(associated(pawrhoij_out(irhoij)%rhoijselect)) then
           ABI_DEALLOCATE(pawrhoij_out(irhoij)%rhoijselect)
         end if
         ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoijselect,(lmn2_size_out))
       end if
       pawrhoij_out(irhoij)%rhoijselect(1:nselect)=pawrhoij_in(jrhoij)%rhoijselect(1:nselect)+0
       if ( (nselect<lmn2_size_out) .and. & 
&           (size(pawrhoij_out(irhoij)%rhoijselect(:))==lmn2_size_out) ) then 
         pawrhoij_out(irhoij)%rhoijselect(nselect+1:lmn2_size_out)=0
       end if 
     end if

!    Optional pointer: indexes of rhoij to be mixed
     lmnmix=pawrhoij_in(jrhoij)%lmnmix_sz
     if (pawrhoij_out(irhoij)%lmnmix_sz/=lmnmix) then
       if (pawrhoij_out(irhoij)%lmnmix_sz>0)  then
         ABI_DEALLOCATE(pawrhoij_out(irhoij)%kpawmix)
       end if
       if (lmnmix>0)  then
         ABI_ALLOCATE(pawrhoij_out(irhoij)%kpawmix,(lmnmix))
       end if
       pawrhoij_out(irhoij)%lmnmix_sz=lmnmix
     end if
     if (lmnmix>0) pawrhoij_out(irhoij)%kpawmix(1:lmnmix)=pawrhoij_in(jrhoij)%kpawmix(1:lmnmix)

!    Optional pointer: gradients of rhoij
     ngrhoij=pawrhoij_in(jrhoij)%ngrhoij
     if (pawrhoij_out(irhoij)%ngrhoij/=ngrhoij) then
       if (pawrhoij_out(irhoij)%ngrhoij>0)  then
         ABI_DEALLOCATE(pawrhoij_out(irhoij)%grhoij)
       end if
       if (ngrhoij>0)  then
         ABI_ALLOCATE(pawrhoij_out(irhoij)%grhoij,(ngrhoij,cplex_out*lmn2_size_out,nspden_out))
       end if
       pawrhoij_out(irhoij)%ngrhoij=ngrhoij
     end if
     if (ngrhoij>0) then
       if (change_dim) then
         ABI_DEALLOCATE(pawrhoij_out(irhoij)%grhoij)
         ABI_ALLOCATE(pawrhoij_out(irhoij)%grhoij,(ngrhoij,cplex_out*lmn2_size_out,nspden_out))
       end if
       if (cplex_out==cplex_in.and.nspden_out==nspden_in) then
         do ispden=1,nspden_out
           do ilmn=1,cplex_out*lmn2_size_out
             pawrhoij_out(irhoij)%grhoij(1:ngrhoij,ilmn,ispden)=pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,ilmn,ispden)
           end do
         end do
       else
         pawrhoij_out(irhoij)%grhoij(:,:,:)=zero
         if (nspden_out==1) then
           if (nspden_in==2) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,1)=pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,1) &
&               +pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,2)
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,1)=pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,1)
             end do
           end if
         else if (nspden_out==2) then
           if (nspden_in==1) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,1)=half*(pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,1) &
&               +pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,2))
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,2)=pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,1)
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,1)=half*(pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,1) &
&               +pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,4))
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,2)=half*(pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,1) &
&               -pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,4))
             end do
           end if
         else if (nspden_out==4) then
           if (nspden_in==1) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,1)=pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,1)
             end do
           else ! nspden_in==2
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,1)=half*(pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,1) &
&               +pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,2))
               pawrhoij_out(irhoij)%grhoij(1:ngrhoij,i_out,4)=half*(pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,1) &
&               -pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,i_in,2))
             end do
           end if
         end if
       end if
     end if

!    Optional pointer: residuals of rhoij
     use_rhoijres=pawrhoij_in(jrhoij)%use_rhoijres
     if (pawrhoij_out(irhoij)%use_rhoijres/=use_rhoijres) then
       if (pawrhoij_out(irhoij)%use_rhoijres>0)  then
         ABI_DEALLOCATE(pawrhoij_out(irhoij)%rhoijres)
       end if
       if (use_rhoijres>0)  then
         ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoijres,(cplex_out*lmn2_size_out,nspden_out))
       end if
       pawrhoij_out(irhoij)%use_rhoijres=use_rhoijres
     end if
     if (use_rhoijres>0) then
       if (change_dim) then
         ABI_DEALLOCATE(pawrhoij_out(irhoij)%rhoijres)
         ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoijres,(cplex_out*lmn2_size_out,nspden_out))
       end if
       if (cplex_out==cplex_in.and.nspden_out==nspden_in) then
         do ispden=1,nspden_out
           do ilmn=1,cplex_out*lmn2_size_out
             pawrhoij_out(irhoij)%rhoijres(ilmn,ispden)=pawrhoij_in(jrhoij)%rhoijres(ilmn,ispden)
           end do
         end do
       else
         pawrhoij_out(irhoij)%rhoijres(:,:)=zero
         if (nspden_out==1) then
           if (nspden_in==2) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijres(i_out,1)=pawrhoij_in(jrhoij)%rhoijres(i_in,1) &
&               +pawrhoij_in(jrhoij)%rhoijres(i_in,2)
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijres(i_out,1)=pawrhoij_in(jrhoij)%rhoijres(i_in,1)
             end do
           end if
         else if (nspden_out==2) then
           if (nspden_in==1) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijres(i_out,1)=half*(pawrhoij_in(jrhoij)%rhoijres(i_in,1) &
&               +pawrhoij_in(jrhoij)%rhoijres(i_in,2))
               pawrhoij_out(irhoij)%rhoijres(i_out,2)=pawrhoij_out(irhoij)%rhoijres(i_out,1)
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijres(i_out,1)=half*(pawrhoij_in(jrhoij)%rhoijres(i_in,1) &
&               +pawrhoij_in(jrhoij)%rhoijres(i_in,4))
               pawrhoij_out(irhoij)%rhoijres(i_out,2)=half*(pawrhoij_in(jrhoij)%rhoijres(i_in,1) &
&               -pawrhoij_in(jrhoij)%rhoijres(i_in,4))
             end do
           end if
         else if (nspden_out==4) then
           if (nspden_in==1) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijres(i_out,1)=pawrhoij_in(jrhoij)%rhoijres(i_in,1)
             end do
           else ! nspden_in==2
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoijres(i_out,1)=half*(pawrhoij_in(jrhoij)%rhoijres(i_in,1) &
&               +pawrhoij_in(jrhoij)%rhoijres(i_in,2))
               pawrhoij_out(irhoij)%rhoijres(i_out,4)=half*(pawrhoij_in(jrhoij)%rhoijres(i_in,1) &
&               -pawrhoij_in(jrhoij)%rhoijres(i_in,2))
             end do
           end if
         end if
       end if
     end if

!    Optional pointer: non-symmetrized rhoij
     use_rhoij_=pawrhoij_in(jrhoij)%use_rhoij_
     if (pawrhoij_out(irhoij)%use_rhoij_/=use_rhoij_) then
       if (pawrhoij_out(irhoij)%use_rhoij_>0)  then
         ABI_DEALLOCATE(pawrhoij_out(irhoij)%rhoij_)
       end if
       if (use_rhoij_>0)  then
         ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoij_,(cplex_out*lmn2_size_out,nspden_out))
       end if
       pawrhoij_out(irhoij)%use_rhoij_=use_rhoij_
     end if
     if (use_rhoij_>0) then
       if (change_dim) then
         if(associated(pawrhoij_out(irhoij)%rhoij_)) then
           ABI_DEALLOCATE(pawrhoij_out(irhoij)%rhoij_)
         end if
         ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoij_,(cplex_out*lmn2_size_out,nspden_out))
       end if
       if (cplex_out==cplex_in.and.nspden_out==nspden_in) then
         do ispden=1,nspden_out
           do ilmn=1,cplex_out*lmn2_size_out

             pawrhoij_out(irhoij)%rhoij_(ilmn,ispden)=pawrhoij_in(jrhoij)%rhoij_(ilmn,ispden)
           end do
         end do
       else
         pawrhoij_out(irhoij)%rhoij_(:,:)=zero
         if (nspden_out==1) then
           if (nspden_in==2) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoij_(i_out,1)=pawrhoij_in(jrhoij)%rhoij_(i_in,1) &
&               +pawrhoij_in(jrhoij)%rhoij_(i_in,2)
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoij_(i_out,1)=pawrhoij_in(jrhoij)%rhoij_(i_in,1)
             end do
           end if
         else if (nspden_out==2) then
           if (nspden_in==1) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoij_(i_out,1)=half*(pawrhoij_in(jrhoij)%rhoij_(i_in,1) &
&               +pawrhoij_in(jrhoij)%rhoij_(i_in,2))
               pawrhoij_out(irhoij)%rhoij_(i_out,2)=pawrhoij_out(irhoij)%rhoij_(i_out,1)
             end do
           else ! nspden_in==4
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoij_(i_out,1)=half*(pawrhoij_in(jrhoij)%rhoij_(i_in,1) &
&               +pawrhoij_in(jrhoij)%rhoij_(i_in,4))
               pawrhoij_out(irhoij)%rhoij_(i_out,2)=half*(pawrhoij_in(jrhoij)%rhoij_(i_in,1) &
&               -pawrhoij_in(jrhoij)%rhoij_(i_in,4))
             end do
           end if
         else if (nspden_out==4) then
           if (nspden_in==1) then
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoij_(i_out,1)=pawrhoij_in(jrhoij)%rhoij_(i_in,1)
             end do
           else ! nspden_in==2
             do ilmn=1,lmn2_size_out
               i_in=cplex_in*ilmn-dplex_in;i_out=cplex_out*ilmn-dplex_out
               pawrhoij_out(irhoij)%rhoij_(i_out,1)=half*(pawrhoij_in(jrhoij)%rhoij_(i_in,1) &
&               +pawrhoij_in(jrhoij)%rhoij_(i_in,2))
               pawrhoij_out(irhoij)%rhoij_(i_out,4)=half*(pawrhoij_in(jrhoij)%rhoij_(i_in,1) &
&               -pawrhoij_in(jrhoij)%rhoij_(i_in,2))
             end do
           end if
         end if
       end if
     end if

   end do ! irhoij
 end if

!Parallel case: do a gather if needed
 if (paral_case==2) then
   call rhoij_free(pawrhoij_copy)
   call rhoij_gather(pawrhoij_out,pawrhoij_copy,-1,mpi_comm_atom)
   call rhoij_free(pawrhoij_out)
   ABI_DATATYPE_DEALLOCATE(pawrhoij_out)

!  Sequential case: fill missing elements
 else if (paral_case==0) then
   if (nrhoij_in<nrhoij_out) then
     do irhoij=nrhoij_in+1,nrhoij_out
       pawrhoij_copy(irhoij)%nrhoijsel=0
       if (pawrhoij_copy(irhoij)%use_rhoijp>0) then
         pawrhoij_copy(irhoij)%rhoijselect=0
         pawrhoij_copy(irhoij)%rhoijp=zero
       end if
       if (pawrhoij_copy(irhoij)%lmnmix_sz>0) pawrhoij_copy(irhoij)%kpawmix=0
       if (pawrhoij_copy(irhoij)%ngrhoij>0) pawrhoij_copy(irhoij)%grhoij=zero
       if (pawrhoij_copy(irhoij)%use_rhoij_>0) pawrhoij_copy(irhoij)%rhoij_=zero
       if (pawrhoij_copy(irhoij)%use_rhoijres>0) pawrhoij_copy(irhoij)%rhoijres=zero
     end do
   end if
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine rhoij_copy
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/rhoij_gather
!! NAME
!! rhoij_gather
!!
!! FUNCTION
!!   (All)Gather pawrhoij datastructures
!!
!! INPUTS
!!  master=master communicator receiving data ; if -1 do a ALLGATHER
!!  mpi_comm_atom= communicator
!!  pawrhoij_in(:)<type(pawrhoij_type)>= input rhoij datastructures on every process
!!  with_grhoij   : optional argument (logical, default=.TRUE.)
!!                  TRUE if pawrhoij%grhoij field is included in the gather operation
!!  with_lmnmix   : optional argument (logical, default=.TRUE.)
!!                  TRUE if pawrhoij%lmnmix field is included in the gather operation
!!  with_rhoijp   : optional argument (logical, default=.TRUE.)
!!                  TRUE if pawrhoij%rhoijp and pawrhoij%rhoijselect fields
!!                       are included in the gather operation
!!  with_rhoijres : optional argument (logical, default=.TRUE.)
!!                 TRUE if pawrhoij%rhoijres field is included in the gather operation
!!  with_rhoij_   : optional argument (logical, default=.TRUE.)
!!                  TRUE if pawrhoij%rhoij_ field is included in the gather operation
!!
!! OUTPUT
!!  pawrhoij_gathered(:)<type(pawrhoij_type)>= output rhoij datastructure
!!
!! NOTES
!!  The gathered structure are ordered like in sequential mode.
!!
!! PARENTS
!!      m_pawrhoij,pawgrnl,pawprt,symrhoij
!!
!! CHILDREN
!!      xsum_mpi
!!
!! SOURCE


 subroutine rhoij_gather(pawrhoij_in,pawrhoij_gathered,master,mpi_comm_atom, &
 &          with_grhoij,with_lmnmix,with_rhoijp,with_rhoijres,with_rhoij_) ! optional arguments


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_gather'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: master,mpi_comm_atom
 logical,intent(in),optional :: with_grhoij,with_lmnmix,with_rhoijp,with_rhoijres,with_rhoij_
!arrays
 type(pawrhoij_type),intent(inout) :: pawrhoij_in(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij_gathered(:)
!Local variables-------------------------------
!scalars
 integer :: buf_dp_size,buf_dp_size_all,buf_int_size,buf_int_size_all
 integer :: cplex,ierr,ii,indx_dp,indx_int,irhoij,isp,jrhoij,lmn2_size,lmnmix,me_atom
 integer :: ngrhoij,nproc_atom,nrhoij_in,nrhoij_in_sum,nrhoij_out,nselect,nspden
 integer :: rhoij_size2,use_rhoijp,use_rhoijres,use_rhoij_
 logical :: my_atmtab_allocated,paral_atom
 logical :: with_grhoij_,with_lmnmix_,with_rhoijp_,with_rhoijres_,with_rhoij__
 character(len=500) :: msg
!arrays
 integer,allocatable :: buf_int(:),buf_int_all(:),count_dp(:),count_int(:)
 integer,allocatable :: disp_dp(:),disp_int(:)
 integer, pointer :: my_atmtab(:)
 real(dp),allocatable :: buf_dp(:),buf_dp_all(:)

! *************************************************************************

 nrhoij_in=size(pawrhoij_in);nrhoij_out=size(pawrhoij_gathered)

 nproc_atom=xcomm_size(mpi_comm_atom)
 me_atom=xcomm_rank(mpi_comm_atom)

 if (nproc_atom==1) then
   if (master==-1.or.me_atom==master) then
     call rhoij_copy(pawrhoij_in,pawrhoij_gathered,.false.,.false.,.false.)
   end if
   return
 end if

!Test on sizes
 nrhoij_in_sum=nrhoij_in
 call xsum_mpi(nrhoij_in_sum,mpi_comm_atom,ierr)
 if (master==-1) then
   if (nrhoij_out/=nrhoij_in_sum) then
     msg='Wrong sizes sum[nrhoij_ij]/=nrhoij_out !'
     call wrtout(std_out,msg,'COLL')
     call leave_new('COLL')
   end if
 else
   if (me_atom==master.and.nrhoij_out/=nrhoij_in_sum) then
     msg='(2) pawrhoij_gathered wrongly allocated !'
     call wrtout(std_out,msg,'COLL')
     call leave_new('COLL')
   end if
 end if

!Optional arguments
 with_grhoij_  =.true.;if (present(with_grhoij))  with_grhoij_  =with_grhoij
 with_lmnmix_  =.true.;if (present(with_lmnmix))  with_lmnmix_  =with_lmnmix
 with_rhoijp_  =.true.;if (present(with_rhoijp))  with_rhoijp_  =with_rhoijp
 with_rhoijres_=.true.;if (present(with_rhoijres))with_rhoijres_=with_rhoijres
 with_rhoij__  =.true.;if (present(with_rhoij_))  with_rhoij__  =with_rhoij_

!Retrieve table of atoms
 paral_atom=.true.;nullify(my_atmtab)
 call get_my_atmtab(mpi_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,nrhoij_in_sum, &
& my_natom_ref=nrhoij_in)

!Compute sizes of buffers
 buf_int_size=0;buf_dp_size=0
 nselect=0;lmnmix=0;ngrhoij=0;rhoij_size2=0
 use_rhoijp=0;use_rhoijres=0;use_rhoij_=0
 do irhoij=1,nrhoij_in
   cplex    =pawrhoij_in(irhoij)%cplex
   lmn2_size=pawrhoij_in(irhoij)%lmn2_size
   nspden   =pawrhoij_in(irhoij)%nspden
   if (with_lmnmix_) lmnmix=pawrhoij_in(irhoij)%lmnmix_sz
   if (with_grhoij_) ngrhoij=pawrhoij_in(irhoij)%ngrhoij
   if (with_rhoijp_) use_rhoijp=pawrhoij_in(irhoij)%use_rhoijp
   if (with_rhoijres_)use_rhoijres=pawrhoij_in(irhoij)%use_rhoijres
   if (with_rhoij__) use_rhoij_=pawrhoij_in(irhoij)%use_rhoij_
   buf_int_size=buf_int_size+15
   if (use_rhoijp>0) then
     nselect=pawrhoij_in(irhoij)%nrhoijsel
     buf_int_size=buf_int_size+nselect
     buf_dp_size=buf_dp_size + cplex*nselect*nspden
   end if
   if (lmnmix>0)       buf_int_size=buf_int_size+lmnmix
   if (ngrhoij>0)      buf_dp_size=buf_dp_size + cplex*lmn2_size*nspden*ngrhoij
   if (use_rhoijres>0) buf_dp_size=buf_dp_size + cplex*lmn2_size*nspden
   if (use_rhoij_>0) then
     rhoij_size2=size(pawrhoij_in(irhoij)%rhoij_,dim=2)
     buf_dp_size=buf_dp_size + cplex*lmn2_size*rhoij_size2
   end if
 end do

!Fill input buffers
 ABI_ALLOCATE(buf_int,(buf_int_size))
 ABI_ALLOCATE(buf_dp ,(buf_dp_size))
 indx_int=1;indx_dp =1
 lmnmix=0;ngrhoij=0;nselect=0;rhoij_size2=0
 use_rhoijp=0;use_rhoijres=0;use_rhoij_=0
 do irhoij=1,nrhoij_in
   cplex    =pawrhoij_in(irhoij)%cplex
   lmn2_size=pawrhoij_in(irhoij)%lmn2_size
   nspden   =pawrhoij_in(irhoij)%nspden
   if (with_lmnmix_) lmnmix=pawrhoij_in(irhoij)%lmnmix_sz
   if (with_grhoij_) ngrhoij=pawrhoij_in(irhoij)%ngrhoij
   if (with_rhoijp_) use_rhoijp=pawrhoij_in(irhoij)%use_rhoijp
   if (with_rhoijp_) nselect=pawrhoij_in(irhoij)%nrhoijsel
   if (with_rhoijres_)use_rhoijres=pawrhoij_in(irhoij)%use_rhoijres
   if (with_rhoij__) use_rhoij_  =pawrhoij_in(irhoij)%use_rhoij_
   if (with_rhoij__) rhoij_size2 =size(pawrhoij_in(irhoij)%rhoij_,dim=2)
   buf_int(indx_int)=my_atmtab(irhoij)               ;indx_int=indx_int+1
   buf_int(indx_int)=cplex                           ;indx_int=indx_int+1
   buf_int(indx_int)=lmn2_size                       ;indx_int=indx_int+1
   buf_int(indx_int)=nspden                          ;indx_int=indx_int+1
   buf_int(indx_int)=nselect                         ;indx_int=indx_int+1
   buf_int(indx_int)=lmnmix                          ;indx_int=indx_int+1
   buf_int(indx_int)=ngrhoij                         ;indx_int=indx_int+1
   buf_int(indx_int)=use_rhoijp                      ;indx_int=indx_int+1
   buf_int(indx_int)=use_rhoijres                    ;indx_int=indx_int+1
   buf_int(indx_int)=use_rhoij_                      ;indx_int=indx_int+1
   buf_int(indx_int)=rhoij_size2                     ;indx_int=indx_int+1
   buf_int(indx_int)=pawrhoij_in(irhoij)%itypat      ;indx_int=indx_int+1
   buf_int(indx_int)=pawrhoij_in(irhoij)%lmn_size    ;indx_int=indx_int+1
   buf_int(indx_int)=pawrhoij_in(irhoij)%nsppol      ;indx_int=indx_int+1
   buf_int(indx_int)=pawrhoij_in(irhoij)%nspinor     ;indx_int=indx_int+1
   if (use_rhoijp>0) then
     buf_int(indx_int:indx_int+nselect-1)=pawrhoij_in(irhoij)%rhoijselect(1:nselect)
     indx_int=indx_int+nselect
     do isp=1,nspden
       buf_dp(indx_dp:indx_dp+cplex*nselect-1)=pawrhoij_in(irhoij)%rhoijp(1:cplex*nselect,isp)
       indx_dp=indx_dp+cplex*nselect
     end do
   end if
   if (lmnmix>0) then
     buf_int(indx_int:indx_int+lmnmix-1)=pawrhoij_in(irhoij)%kpawmix(1:lmnmix)
     indx_int=indx_int+lmnmix
   end if
   if (ngrhoij>0) then
     do isp=1,nspden
       do ii=1,cplex*lmn2_size
         buf_dp(indx_dp:indx_dp+ngrhoij-1)=pawrhoij_in(irhoij)%grhoij(1:ngrhoij,ii,isp)
         indx_dp=indx_dp+ngrhoij
       end do
     end do
   end if
   if (use_rhoijres>0) then
     do isp=1,nspden
       buf_dp(indx_dp:indx_dp+cplex*lmn2_size-1)=pawrhoij_in(irhoij)%rhoijres(1:cplex*lmn2_size,isp)
       indx_dp=indx_dp+cplex*lmn2_size
     end do
   end if
   if (use_rhoij_>0) then
     do isp=1,rhoij_size2
       buf_dp(indx_dp:indx_dp+cplex*lmn2_size-1)=pawrhoij_in(irhoij)%rhoij_(1:cplex*lmn2_size,isp)
       indx_dp=indx_dp+cplex*lmn2_size
     end do
   end if
 end do

!Check
 if ((indx_int-1/=buf_int_size).or.(indx_dp-1/=buf_dp_size)) then
   write(msg,*) 'Wrong buffer sizes: buf_int_size=',buf_int_size,' buf_dp_size=',buf_dp_size
   call wrtout(std_out,msg,'COLL')
   call leave_new('COLL')
 end if

!Prepare communications
 ABI_ALLOCATE(count_int,(nproc_atom))
 ABI_ALLOCATE(disp_int,(nproc_atom))
 ABI_ALLOCATE(count_dp,(nproc_atom))
 ABI_ALLOCATE(disp_dp,(nproc_atom))
 call xallgather_mpi(buf_int_size,count_int,mpi_comm_atom,ierr)
 call xallgather_mpi(buf_dp_size ,count_dp ,mpi_comm_atom,ierr)
 disp_int(1)=0;disp_dp(1)=0
 do ii=2,nproc_atom
   disp_int(ii)=disp_int(ii-1)+count_int(ii-1)
   disp_dp (ii)=disp_dp (ii-1)+count_dp (ii-1)
 end do
 buf_int_size_all=sum(count_int)
 buf_dp_size_all =sum(count_dp)
 if (master==-1.or.me_atom==master) then
   ABI_ALLOCATE(buf_int_all,(buf_int_size_all))
   ABI_ALLOCATE(buf_dp_all ,(buf_dp_size_all))
 else
   ABI_ALLOCATE(buf_int_all,(1))
   ABI_ALLOCATE(buf_dp_all ,(1))
 end if

!Communicate
 if (master==-1) then
   call xallgatherv_mpi(buf_int,buf_int_size,buf_int_all,count_int,disp_int,mpi_comm_atom,ierr)
   call xallgatherv_mpi(buf_dp ,buf_dp_size ,buf_dp_all ,count_dp ,disp_dp ,mpi_comm_atom,ierr)
 else
   call xgatherv_mpi(buf_int,buf_int_size,buf_int_all,count_int,disp_int,master,mpi_comm_atom,ierr)
   call xgatherv_mpi(buf_dp ,buf_dp_size ,buf_dp_all ,count_dp ,disp_dp ,master,mpi_comm_atom,ierr)
 end if

!Retrieve data from output buffer
 if (master==-1.or.me_atom==master) then
   indx_int=1;indx_dp=1
   call rhoij_free(pawrhoij_gathered)
   do irhoij=1,nrhoij_out
     jrhoij      =buf_int_all(indx_int)    ;indx_int=indx_int+1
     cplex       =buf_int_all(indx_int)    ;indx_int=indx_int+1
     lmn2_size   =buf_int_all(indx_int)    ;indx_int=indx_int+1
     nspden      =buf_int_all(indx_int)    ;indx_int=indx_int+1
     nselect     =buf_int_all(indx_int)    ;indx_int=indx_int+1
     lmnmix      =buf_int_all(indx_int)    ;indx_int=indx_int+1
     ngrhoij     =buf_int_all(indx_int)    ;indx_int=indx_int+1
     use_rhoijp  =buf_int_all(indx_int)    ;indx_int=indx_int+1
     use_rhoijres=buf_int_all(indx_int)    ;indx_int=indx_int+1
     use_rhoij_  =buf_int_all(indx_int)    ;indx_int=indx_int+1
     rhoij_size2 =buf_int_all(indx_int)    ;indx_int=indx_int+1
     pawrhoij_gathered(jrhoij)%itypat=buf_int_all(indx_int)   ;indx_int=indx_int+1
     pawrhoij_gathered(jrhoij)%lmn_size=buf_int_all(indx_int) ;indx_int=indx_int+1
     pawrhoij_gathered(jrhoij)%nsppol=buf_int_all(indx_int)   ;indx_int=indx_int+1
     pawrhoij_gathered(jrhoij)%nspinor=buf_int_all(indx_int)  ;indx_int=indx_int+1
     pawrhoij_gathered(jrhoij)%cplex=cplex
     pawrhoij_gathered(jrhoij)%lmn2_size=lmn2_size
     pawrhoij_gathered(jrhoij)%nspden=nspden
     pawrhoij_gathered(jrhoij)%nrhoijsel=nselect
     pawrhoij_gathered(jrhoij)%lmnmix_sz=lmnmix
     pawrhoij_gathered(jrhoij)%ngrhoij=ngrhoij
     pawrhoij_gathered(jrhoij)%use_rhoijp=use_rhoijp
     pawrhoij_gathered(jrhoij)%use_rhoijres=use_rhoijres
     pawrhoij_gathered(jrhoij)%use_rhoij_=use_rhoij_
     if (use_rhoijp>0) then
       ABI_ALLOCATE(pawrhoij_gathered(jrhoij)%rhoijselect,(nselect))
       pawrhoij_gathered(jrhoij)%rhoijselect(1:nselect)=buf_int_all(indx_int:indx_int+nselect-1)
       indx_int=indx_int+nselect
       ABI_ALLOCATE(pawrhoij_gathered(jrhoij)%rhoijp,(cplex*nselect,nspden))
       do isp=1,nspden
         pawrhoij_gathered(jrhoij)%rhoijp(1:cplex*nselect,isp)=buf_dp_all(indx_dp:indx_dp+cplex*nselect-1)
         indx_dp=indx_dp+cplex*nselect
       end do
     end if
     if (lmnmix>0) then
       ABI_ALLOCATE(pawrhoij_gathered(jrhoij)%kpawmix,(lmnmix))
       pawrhoij_gathered(jrhoij)%kpawmix(1:lmnmix)=buf_int_all(indx_int:indx_int+lmnmix-1)
       indx_int=indx_int+lmnmix
     end if
     if (ngrhoij>0) then
       ABI_ALLOCATE(pawrhoij_gathered(jrhoij)%grhoij,(ngrhoij,cplex*lmn2_size,nspden))
       do isp=1,nspden
         do ii=1,cplex*lmn2_size
           pawrhoij_gathered(jrhoij)%grhoij(1:ngrhoij,ii,isp)=buf_dp_all(indx_dp:indx_dp+ngrhoij-1)
           indx_dp=indx_dp+ngrhoij
         end do
       end do
     end if
     if (use_rhoijres>0) then
       ABI_ALLOCATE(pawrhoij_gathered(jrhoij)%rhoijres,(cplex*lmn2_size,nspden))
       do isp=1,nspden
         pawrhoij_gathered(jrhoij)%rhoijres(1:cplex*lmn2_size,isp)=buf_dp_all(indx_dp:indx_dp+cplex*lmn2_size-1)
         indx_dp=indx_dp+cplex*lmn2_size
       end do
     end if
     if (use_rhoij_>0) then
       ABI_ALLOCATE(pawrhoij_gathered(jrhoij)%rhoij_,(cplex*lmn2_size,rhoij_size2))
       do isp=1,rhoij_size2
         pawrhoij_gathered(jrhoij)%rhoij_(1:cplex*lmn2_size,isp)=buf_dp_all(indx_dp:indx_dp+cplex*lmn2_size-1)
         indx_dp=indx_dp+cplex*lmn2_size
       end do
     end if
   end do
   if ((indx_int/=1+buf_int_size_all).or.(indx_dp/=1+buf_dp_size_all)) then
     write(msg,*) 'Wrong buffer sizes: buf_int_size_all=',buf_int_size_all,' buf_dp_size_all=',buf_dp_size_all
     call wrtout(std_out,msg,'COLL')
     call leave_new('COLL')
   end if
 end if

!Free memory
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)
 ABI_DEALLOCATE(buf_int)
 ABI_DEALLOCATE(count_int)
 ABI_DEALLOCATE(disp_int)
 ABI_DEALLOCATE(buf_dp)
 ABI_DEALLOCATE(count_dp)
 ABI_DEALLOCATE(disp_dp)
 ABI_DEALLOCATE(buf_int_all)
 ABI_DEALLOCATE(buf_dp_all)

end subroutine rhoij_gather
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/rhoij_bcast
!! NAME
!! rhoij_gather
!!
!! FUNCTION
!!   Broadcast pawrhoij datastructures
!!   Can take into account a distribution of data over a "atom" communicator
!!
!! INPUTS
!!  master=master communicator receiving data
!!  mpi_comm= MPI communicator
!!  mpi_comm_atom= --optional-- MPI communicator over atoms
!!  pawrhoij_in(:)<type(pawrhoij_type)>= input rhoij datastructures on master process
!!
!! OUTPUT
!!  pawrhoij_out(:)<type(pawrhoij_type)>= output rhoij datastructure on every process
!!    Eventually distributed according to mpi_comm_atom communicator
!!
!! PARENTS
!!      loper3,respfn
!!
!! CHILDREN
!!      xallgatherv_mpi,xgatherv_mpi,xscatter_mpi,xsum_mpi
!!
!! SOURCE

 subroutine rhoij_bcast(pawrhoij_in,pawrhoij_out,master,mpi_comm,mpi_comm_atom)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_bcast'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: master,mpi_comm
 integer,intent(in),optional :: mpi_comm_atom
!arrays
 type(pawrhoij_type),intent(inout) :: pawrhoij_in(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij_out(:)
!Local variables-------------------------------
!scalars
 integer :: buf_dp_size,buf_dp_size_all,buf_int_size,buf_int_size_all
 integer :: cplex,ierr,ii,indx_dp,indx_int,irhoij,isp,jrhoij,lmn2_size,lmnmix,me
 integer :: mpi_comm_atom_,ngrhoij,nproc,nproc_atom,nrhoij_in,nrhoij_out,nrhoij_out_all
 integer :: nselect,nspden,rhoij_size2,use_rhoijp,use_rhoijres,use_rhoij_
 logical :: my_atmtab_allocated,paral_atom
 character(len=500) :: msg
!arrays
 integer,allocatable :: atmtab(:),buf_int(:),buf_int_all(:)
 integer,allocatable :: count_dp(:),count_int(:),disp_dp(:),disp_int(:)
 integer,pointer :: my_atmtab(:)
 real(dp),allocatable :: buf_dp(:),buf_dp_all(:)

! *************************************************************************

!Load MPI "atom" distribution data
 mpi_comm_atom_=xmpi_self;nproc_atom=1
 if (present(mpi_comm_atom)) then
   mpi_comm_atom_=mpi_comm_atom
   nproc_atom=xcomm_size(mpi_comm_atom_)
   paral_atom=(nproc_atom>1)
   if (mpi_comm_atom_/=mpi_comm.and.nproc_atom/=1) then
     msg='wrong mpi_comm_atom communicator !'
     call wrtout(std_out,msg,'COLL')
     call leave_new('COLL')
   end if
 end if

!Load global MPI data
 me=xcomm_rank(mpi_comm)
 nproc=xcomm_size(mpi_comm)

!Just copy in case of a sequential run
 if (nproc==1.and.nproc_atom==1) then
   call rhoij_copy(pawrhoij_in,pawrhoij_out,.false.,.false.,.false.)
   return
 end if

!Retrieve and test pawrhoij sizes
 nrhoij_in=0;if (me==master) nrhoij_in=size(pawrhoij_in)
 nrhoij_out=size(pawrhoij_out);nrhoij_out_all=nrhoij_out
 if (paral_atom) then
   call xsum_mpi(nrhoij_out_all,mpi_comm_atom_,ierr)
 end if
 if (me==master.and.nrhoij_in/=nrhoij_out_all) then
   msg='pawrhoij_in or pawrhoij_out wrongly allocated !'
   call wrtout(std_out,msg,'COLL')
   call leave_new('COLL')
 end if

!Retrieve table(s) of atoms (if necessary)
 if (paral_atom) then
   nullify(my_atmtab)
   call get_my_atmtab(mpi_comm_atom_,my_atmtab,my_atmtab_allocated,paral_atom, &
&                     nrhoij_out_all,my_natom_ref=nrhoij_out)
   ABI_ALLOCATE(disp_int,(nproc_atom))
   ABI_ALLOCATE(count_int,(nproc_atom))
   call xallgather_mpi(nrhoij_out,count_int,mpi_comm_atom_,ierr)
   disp_int(1)=0
   do ii=2,nproc_atom
     disp_int(ii)=disp_int(ii-1)+count_int(ii-1)
   end do
   ABI_ALLOCATE(atmtab,(nrhoij_in))
   call xgatherv_mpi(my_atmtab,nrhoij_out,atmtab,count_int,disp_int,&
&                    master,mpi_comm_atom_,ierr)
   ABI_DEALLOCATE(disp_int)
   ABI_DEALLOCATE(count_int)
 end if

!Compute size of input buffers
 buf_int_size=0;buf_dp_size=0
 do irhoij=1,nrhoij_in
   jrhoij=irhoij;if (paral_atom) jrhoij=atmtab(irhoij)
   cplex       =pawrhoij_in(jrhoij)%cplex
   lmn2_size   =pawrhoij_in(jrhoij)%lmn2_size
   nspden      =pawrhoij_in(jrhoij)%nspden
   lmnmix      =pawrhoij_in(jrhoij)%lmnmix_sz
   ngrhoij     =pawrhoij_in(jrhoij)%ngrhoij
   use_rhoijp  =pawrhoij_in(jrhoij)%use_rhoijp
   use_rhoijres=pawrhoij_in(jrhoij)%use_rhoijres
   use_rhoij_  =pawrhoij_in(jrhoij)%use_rhoij_
   buf_int_size=buf_int_size+15
   if (ngrhoij>0) buf_dp_size=buf_dp_size+cplex*lmn2_size*nspden*ngrhoij
   if (use_rhoijres>0) buf_dp_size=buf_dp_size+cplex*lmn2_size*nspden
   if (use_rhoijp>0) then
     nselect=pawrhoij_in(jrhoij)%nrhoijsel
     buf_int_size=buf_int_size+nselect
     buf_dp_size=buf_dp_size+cplex*nselect*nspden
   end if
   if (use_rhoij_>0) then
     rhoij_size2=size(pawrhoij_in(jrhoij)%rhoij_,dim=2)
     buf_dp_size=buf_dp_size+cplex*lmn2_size*rhoij_size2
   end if
 end do

!Prepare buffers/tabs for communication
 if (paral_atom) then
   ABI_ALLOCATE(count_int,(nproc_atom))
   ABI_ALLOCATE(count_dp,(nproc_atom))
   ABI_ALLOCATE(disp_int,(nproc_atom))
   ABI_ALLOCATE(disp_dp,(nproc_atom))
   call xallgather_mpi(buf_int_size,count_int,mpi_comm_atom_,ierr)
   call xallgather_mpi(buf_dp_size ,count_dp ,mpi_comm_atom_,ierr)
   disp_int(1)=0;disp_dp(1)=0
   do ii=2,nproc_atom
     disp_int(ii)=disp_int(ii-1)+count_int(ii-1)
     disp_dp (ii)=disp_dp (ii-1)+count_dp (ii-1)
   end do
   buf_int_size_all=sum(count_int)
   buf_dp_size_all =sum(count_dp)
 else
   buf_int_size_all=buf_int_size
   buf_dp_size_all =buf_dp_size
 end if
 ABI_ALLOCATE(buf_int,(buf_int_size))
 ABI_ALLOCATE(buf_dp ,(buf_dp_size))
 if (me==master) then
   ABI_ALLOCATE(buf_int_all,(buf_int_size_all))
   ABI_ALLOCATE(buf_dp_all ,(buf_dp_size_all))
 else
   ABI_ALLOCATE(buf_int_all,(1))
   ABI_ALLOCATE(buf_dp_all ,(1))
 end if

!Fill input buffers
 if (me==master) then
   indx_int=1;indx_dp =1
   do irhoij=1,nrhoij_in
     jrhoij=irhoij;if (paral_atom) jrhoij=atmtab(irhoij)
     cplex       =pawrhoij_in(jrhoij)%cplex
     lmn2_size   =pawrhoij_in(jrhoij)%lmn2_size
     nspden      =pawrhoij_in(jrhoij)%nspden
     lmnmix      =pawrhoij_in(jrhoij)%lmnmix_sz
     ngrhoij     =pawrhoij_in(jrhoij)%ngrhoij
     use_rhoijp  =pawrhoij_in(jrhoij)%use_rhoijp
     nselect     =pawrhoij_in(jrhoij)%nrhoijsel
     use_rhoijres=pawrhoij_in(jrhoij)%use_rhoijres
     use_rhoij_  =pawrhoij_in(jrhoij)%use_rhoij_
     rhoij_size2 =size(pawrhoij_in(jrhoij)%rhoij_,dim=2)
     buf_int_all(indx_int)=jrhoij                      ;indx_int=indx_int+1 ! Not used !
     buf_int_all(indx_int)=cplex                       ;indx_int=indx_int+1
     buf_int_all(indx_int)=lmn2_size                   ;indx_int=indx_int+1
     buf_int_all(indx_int)=nspden                      ;indx_int=indx_int+1
     buf_int_all(indx_int)=nselect                     ;indx_int=indx_int+1
     buf_int_all(indx_int)=lmnmix                      ;indx_int=indx_int+1
     buf_int_all(indx_int)=ngrhoij                     ;indx_int=indx_int+1
     buf_int_all(indx_int)=use_rhoijp                  ;indx_int=indx_int+1
     buf_int_all(indx_int)=use_rhoijres                ;indx_int=indx_int+1
     buf_int_all(indx_int)=use_rhoij_                  ;indx_int=indx_int+1
     buf_int_all(indx_int)=rhoij_size2                 ;indx_int=indx_int+1
     buf_int_all(indx_int)=pawrhoij_in(jrhoij)%itypat  ;indx_int=indx_int+1
     buf_int_all(indx_int)=pawrhoij_in(jrhoij)%lmn_size;indx_int=indx_int+1
     buf_int_all(indx_int)=pawrhoij_in(jrhoij)%nsppol  ;indx_int=indx_int+1
     buf_int_all(indx_int)=pawrhoij_in(jrhoij)%nspinor ;indx_int=indx_int+1
     if (use_rhoijp>0) then
       buf_int_all(indx_int:indx_int+nselect-1)=pawrhoij_in(jrhoij)%rhoijselect(1:nselect)
       indx_int=indx_int+nselect
       do isp=1,nspden
         buf_dp_all(indx_dp:indx_dp+cplex*nselect-1)=pawrhoij_in(jrhoij)%rhoijp(1:cplex*nselect,isp)
         indx_dp=indx_dp+cplex*nselect
       end do
     end if
     if (lmnmix>0) then
       buf_int_all(indx_int:indx_int+lmnmix-1)=pawrhoij_in(jrhoij)%kpawmix(1:lmnmix)
       indx_int=indx_int+lmnmix
     end if
     if (ngrhoij>0) then
       do isp=1,nspden
         do ii=1,cplex*lmn2_size
           buf_dp_all(indx_dp:indx_dp+ngrhoij-1)=pawrhoij_in(jrhoij)%grhoij(1:ngrhoij,ii,isp)
           indx_dp=indx_dp+ngrhoij
         end do
       end do
     end if
     if (use_rhoijres>0) then
       do isp=1,nspden
         buf_dp_all(indx_dp:indx_dp+cplex*lmn2_size-1)=pawrhoij_in(jrhoij)%rhoijres(1:cplex*lmn2_size,isp)
         indx_dp=indx_dp+cplex*lmn2_size
       end do
     end if
     if (use_rhoij_>0) then
       do isp=1,rhoij_size2
         buf_dp_all(indx_dp:indx_dp+cplex*lmn2_size-1)=pawrhoij_in(jrhoij)%rhoij_(1:cplex*lmn2_size,isp)
         indx_dp=indx_dp+cplex*lmn2_size
       end do
     end if
   end do
! Check
  if ((indx_int-1/=buf_int_size_all).or.(indx_dp-1/=buf_dp_size_all)) then
    msg='(1) Wrong buffer sizes !'
    call wrtout(std_out,msg,'COLL')
    call leave_new('COLL')
  end if
 end if ! me=master

!Communicate
 if (paral_atom) then
   call xscatterv_mpi(buf_int_all,count_int,disp_int,buf_int,buf_int_size,master,mpi_comm,ierr)
   call xscatterv_mpi(buf_dp_all ,count_dp ,disp_dp ,buf_dp ,buf_dp_size ,master,mpi_comm,ierr)
 else
   buf_int=buf_int_all;buf_dp=buf_dp_all
   call xcast_mpi(buf_int,master,mpi_comm,ierr)
   call xcast_mpi(buf_dp ,master,mpi_comm,ierr)
 end if

!Retrieve data from output buffer
 indx_int=1;indx_dp=1
 call rhoij_free(pawrhoij_out)
 do irhoij=1,nrhoij_out
   jrhoij      =buf_int(indx_int);indx_int=indx_int+1 ! Not used !
   cplex       =buf_int(indx_int);indx_int=indx_int+1
   lmn2_size   =buf_int(indx_int);indx_int=indx_int+1
   nspden      =buf_int(indx_int);indx_int=indx_int+1
   nselect     =buf_int(indx_int);indx_int=indx_int+1
   lmnmix      =buf_int(indx_int);indx_int=indx_int+1
   ngrhoij     =buf_int(indx_int);indx_int=indx_int+1
   use_rhoijp  =buf_int(indx_int);indx_int=indx_int+1
   use_rhoijres=buf_int(indx_int);indx_int=indx_int+1
   use_rhoij_  =buf_int(indx_int);indx_int=indx_int+1
   rhoij_size2 =buf_int(indx_int);indx_int=indx_int+1
   pawrhoij_out(irhoij)%itypat=buf_int(indx_int)  ;indx_int=indx_int+1
   pawrhoij_out(irhoij)%lmn_size=buf_int(indx_int);indx_int=indx_int+1
   pawrhoij_out(irhoij)%nsppol=buf_int(indx_int)  ;indx_int=indx_int+1
   pawrhoij_out(irhoij)%nspinor=buf_int(indx_int) ;indx_int=indx_int+1
   pawrhoij_out(irhoij)%cplex=cplex
   pawrhoij_out(irhoij)%lmn2_size=lmn2_size
   pawrhoij_out(irhoij)%nspden=nspden
   pawrhoij_out(irhoij)%nrhoijsel=nselect
   pawrhoij_out(irhoij)%lmnmix_sz=lmnmix
   pawrhoij_out(irhoij)%ngrhoij=ngrhoij
   pawrhoij_out(irhoij)%use_rhoijp=use_rhoijp
   pawrhoij_out(irhoij)%use_rhoijres=use_rhoijres
   pawrhoij_out(irhoij)%use_rhoij_=use_rhoij_
   if (use_rhoijp>0) then
     ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoijselect,(nselect))
     pawrhoij_out(irhoij)%rhoijselect(1:nselect)=buf_int(indx_int:indx_int+nselect-1)
     indx_int=indx_int+nselect
     ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoijp,(cplex*nselect,nspden))
     do isp=1,nspden
       pawrhoij_out(irhoij)%rhoijp(1:cplex*nselect,isp)=buf_dp(indx_dp:indx_dp+cplex*nselect-1)
       indx_dp=indx_dp+cplex*nselect
     end do
   end if
   if (lmnmix>0) then
     ABI_ALLOCATE(pawrhoij_out(irhoij)%kpawmix,(lmnmix))
     pawrhoij_out(irhoij)%kpawmix(1:lmnmix)=buf_int(indx_int:indx_int+lmnmix-1)
     indx_int=indx_int+lmnmix
   end if
   if (ngrhoij>0) then
     ABI_ALLOCATE(pawrhoij_out(irhoij)%grhoij,(ngrhoij,cplex*lmn2_size,nspden))
     do isp=1,nspden
       do ii=1,cplex*lmn2_size
         pawrhoij_out(irhoij)%grhoij(1:ngrhoij,ii,isp)=buf_dp(indx_dp:indx_dp+ngrhoij-1)
         indx_dp=indx_dp+ngrhoij
       end do
     end do
   end if
   if (use_rhoijres>0) then
     ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoijres,(cplex*lmn2_size,nspden))
     do isp=1,nspden
       pawrhoij_out(irhoij)%rhoijres(1:cplex*lmn2_size,isp)=buf_dp(indx_dp:indx_dp+cplex*lmn2_size-1)
       indx_dp=indx_dp+cplex*lmn2_size
     end do
   end if
   if (use_rhoij_>0) then
     ABI_ALLOCATE(pawrhoij_out(irhoij)%rhoij_,(cplex*lmn2_size,rhoij_size2))
     do isp=1,rhoij_size2
       pawrhoij_out(irhoij)%rhoij_(1:cplex*lmn2_size,isp)=buf_dp(indx_dp:indx_dp+cplex*lmn2_size-1)
       indx_dp=indx_dp+cplex*lmn2_size
     end do
   end if
 end do
!Check
 if ((indx_int/=1+buf_int_size).or.(indx_dp/=1+buf_dp_size)) then
   msg='(2) Wrong buffer sizes !'
   call wrtout(std_out,msg,'COLL')
   call leave_new('COLL')
 end if

!Free memory
 ABI_DEALLOCATE(buf_int)
 ABI_DEALLOCATE(buf_dp)
 ABI_DEALLOCATE(buf_int_all)
 ABI_DEALLOCATE(buf_dp_all)
 if (paral_atom) then
   ABI_DEALLOCATE(count_int)
   ABI_DEALLOCATE(count_dp)
   ABI_DEALLOCATE(disp_int)
   ABI_DEALLOCATE(disp_dp)
   ABI_DEALLOCATE(atmtab)
   call free_my_atmtab(my_atmtab,my_atmtab_allocated)
 end if

end subroutine rhoij_bcast
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/rhoij_io
!! NAME
!! rhoij_io
!!
!! FUNCTION
!! IO method for pawrhoij datastructures.
!!
!! INPUTS
!!  unitfi=Unit number for IO file (already opened in the caller).
!!  nsppol_in=Number of independent spin polarizations. Only used for reading.
!!  nspinor_in=Number of spinorial components. Only used for reading.
!!  nspden_in=Number of spin-density components. only used for reading.
!!  nlmn_type(ntypat)= Number of (l,m,n) elements for the paw basis for each type of atom. Only used for reading.
!!  typat(natom) =Type of each atom.
!!  headform=Format of the abinit header (only used for reading as we need to know how to read
!!    the data. Writing is always done using the latest headform.
!!  rdwr_mode(len=*)=String defining the IO mode. Possible values (not case sensitive):
!!    "W"= For writing to unitfi
!!    "R"= For reading from unitfi
!!    "E"= For echoing.
!!    "D"= for debug
!!  [form(len=*)]= String defining the file format. Defaults to Fortran binary mode i.e., "unformatted"
!!  Other possible values are (case insensitive):
!!    "formatted"=For IO on a file open in formatted mode.
!!  [natinc]=Defines the increment in the loop over natom used for echoing the pawrhoij(natom) datastructures.
!!    If not specified, only the first and the last atom will be printed.
!!
!! SIDE EFFECTS
!!  pawrhoij(:)<type(pawrhoij_type)>= rhoij datastructure.
!!   if rdwr_mode="W", it will be written on unit unitfi using the file format given by form.
!!   if rdwr_mode="R", pawrhoij will be read and initialized from unit unitfi that has been
!!      opened with form=form.
!!   if rdwr_mode="E", the routines only echoes the content of the structure.
!!
!! PARENTS
!!      m_header,m_qparticles
!!
!! CHILDREN
!!      xsum_mpi
!!
!! SOURCE

subroutine rhoij_io(pawrhoij,unitfi,nsppol_in,nspinor_in,nspden_in,nlmn_type,typat,headform,rdwr_mode,form,natinc)

 use m_io_tools,   only : flush_unit
 use m_fstrings,   only : toupper
 use m_pawio, only: pawio_print_ij

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_io'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unitfi,headform,nspden_in,nspinor_in,nsppol_in
 integer,optional,intent(in) :: natinc
 character(len=*),intent(in) :: rdwr_mode
 character(len=*),optional,intent(in) :: form
!arrays
 integer,intent(in) :: typat(:),nlmn_type(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables-------------------------------
!scalars
 integer :: iatom,natom,ispden,bsize,ii,jj,nselect,my_cplex,my_natinc,my_nspden
 logical :: isbinary
!arrays
 integer,allocatable :: ibuffer(:),nsel44(:,:),nsel56(:)
 real(dp), allocatable :: buffer(:)

! *************************************************************************

 natom=SIZE(pawrhoij);if (natom==0) return

 isbinary=.TRUE.
 if (PRESENT(form)) then
   if (toupper(form)=="FORMATTED") isbinary=.FALSE.
 end if

 select case (rdwr_mode(1:1))

   case ("R","r") ! Reading the Rhoij tab.

     if ((headform>=44).and.(headform<56)) then

       ABI_ALLOCATE(nsel44,(nspden_in,natom))

       if (isbinary) then
         read(unitfi  ) ((nsel44(ispden,iatom),ispden=1,nspden_in),iatom=1,natom)
       else
         read(unitfi,*) ((nsel44(ispden,iatom),ispden=1,nspden_in),iatom=1,natom)
       end if

       call rhoij_alloc(pawrhoij,1,nspden_in,nspinor_in,nsppol_in,typat,lmnsize=nlmn_type)
       do iatom=1,natom
         pawrhoij(iatom)%nrhoijsel=nsel44(1,iatom)
       end do

       bsize=sum(nsel44)
       ABI_ALLOCATE(ibuffer,(bsize))
       ABI_ALLOCATE(buffer,(bsize))

       if (isbinary) then
         read(unitfi  ) ibuffer(:),buffer(:)
       else
         read(unitfi,*) ibuffer(:),buffer(:)
       end if

       ii=0
       do iatom=1,natom
         nselect=nsel44(1,iatom)
         pawrhoij(iatom)%rhoijselect(1:nselect)=ibuffer(ii+1:ii+nselect)
         do ispden=1,nspden_in
           pawrhoij(iatom)%rhoijp(1:nselect,ispden)=buffer(ii+1:ii+nselect)
           ii=ii+nselect
         end do
       end do
       ABI_DEALLOCATE(ibuffer)
       ABI_DEALLOCATE(buffer)
       ABI_DEALLOCATE(nsel44)

     else if (headform>=56) then
       ABI_ALLOCATE(nsel56,(natom))

       if (headform==56) then
         if (isbinary) then
           read(unitfi  ) (nsel56(iatom),iatom=1,natom),my_cplex
         else
           read(unitfi,*) (nsel56(iatom),iatom=1,natom),my_cplex
         end if
         my_nspden=nspden_in
       else
         if (isbinary) then
           read(unitfi  ) (nsel56(iatom),iatom=1,natom),my_cplex,my_nspden
         else
           read(unitfi,*) (nsel56(iatom),iatom=1,natom),my_cplex,my_nspden
         end if
       end if

       call rhoij_alloc(pawrhoij,my_cplex,my_nspden,nspinor_in,nsppol_in,typat,lmnsize=nlmn_type)
       do iatom=1,natom
         pawrhoij(iatom)%nrhoijsel=nsel56(iatom)
       end do
       bsize=sum(nsel56)
       ABI_ALLOCATE(ibuffer,(bsize))
       ABI_ALLOCATE(buffer,(bsize*nspden_in*my_cplex))

       if (isbinary) then
         read(unitfi  ) ibuffer(:),buffer(:)
       else
         read(unitfi,*) ibuffer(:),buffer(:)
       end if

       ii=0;jj=0
       do iatom=1,natom
         nselect=nsel56(iatom)
         pawrhoij(iatom)%rhoijselect(1:nselect)=ibuffer(ii+1:ii+nselect)
         ii=ii+nselect
         do ispden=1,nspden_in
           pawrhoij(iatom)%rhoijp(1:my_cplex*nselect,ispden)=buffer(jj+1:jj+my_cplex*nselect)
           jj=jj+my_cplex*nselect
         end do
       end do
       ABI_DEALLOCATE(ibuffer)
       ABI_DEALLOCATE(buffer)
       ABI_DEALLOCATE(nsel56)
     end if

   case ("W","w") ! Writing the Rhoij tab. Latest format is used.

     ABI_ALLOCATE(nsel56,(natom))
     my_cplex =pawrhoij(1)%cplex
     my_nspden=pawrhoij(1)%nspden
     do iatom=1,natom
       nsel56(iatom)=pawrhoij(iatom)%nrhoijsel
     end do

     if (isbinary) then
       write(unitfi  ) (nsel56(iatom),iatom=1,natom),my_cplex,my_nspden
     else
       write(unitfi,*) (nsel56(iatom),iatom=1,natom),my_cplex,my_nspden
     end if

     bsize=sum(nsel56)
     ABI_ALLOCATE(ibuffer,(bsize))
     ABI_ALLOCATE(buffer,(bsize*my_nspden*my_cplex))
     ii=0;jj=0
     do iatom=1,natom
       nselect=nsel56(iatom)
       ibuffer(ii+1:ii+nselect)=pawrhoij(iatom)%rhoijselect(1:nselect)
       ii=ii+nselect
       do ispden=1,my_nspden
         buffer(jj+1:jj+my_cplex*nselect)=pawrhoij(iatom)%rhoijp(1:my_cplex*nselect,ispden)
         jj=jj+my_cplex*nselect
       end do
     end do

     if (isbinary) then
       write(unitfi  ) ibuffer(:),buffer(:)
     else
       write(unitfi,*) ibuffer(:),buffer(:)
     end if

     ABI_DEALLOCATE(ibuffer)
     ABI_DEALLOCATE(buffer)
     ABI_DEALLOCATE(nsel56)

   case ("E","e") ! Echoing
     my_natinc=1; if(natom>1) my_natinc=natom-1
     if (PRESENT(natinc)) my_natinc = natinc ! user-defined increment.
     ABI_ALLOCATE(ibuffer,(0))
     do iatom=1,natom,my_natinc
       do ispden=1,pawrhoij(iatom)%nspden
         write(unitfi, '(a,i4,a,i1,a)' ) ' rhoij(',iatom,',',ispden,')=  (max 12 non-zero components will be written)'
         call pawio_print_ij(unitfi,pawrhoij(iatom)%rhoijp(:,ispden),&
&         pawrhoij(iatom)%nrhoijsel,&
&         pawrhoij(iatom)%cplex,&
&         pawrhoij(iatom)%lmn_size,-1,ibuffer,1,0,&
&         pawrhoij(iatom)%rhoijselect(:),-1.d0,1)
       end do
     end do
     ABI_DEALLOCATE(ibuffer)

   case default
     MSG_ERROR("Wrong rdwr_mode"//TRIM(rdwr_mode))
!     call wrtout(std_out,"Wrong rdwr_mode"//TRIM(rdwr_mode,'COLL')
!     call leave_new('COLL')

 end select

end subroutine rhoij_io
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/rhoij_unpack
!! NAME
!! rhoij_unpack
!!
!! FUNCTION
!!  Unpack the values store in rhoijp copying them to the rhoij_ array.
!!
!! SIDE EFFECTS
!!  rhoij(:) <pawrhoij_type)>= input/output datastructure
!!   * In output the rhoij_ array is filled with the values stored in the packed array rhoijp.
!!   * If use_rhoij_/=1, rhoij_ is allocated and the corresponding flag is set to 1.
!!
!! PARENTS
!!      paw_qpscgw
!!
!! CHILDREN
!!      xsum_mpi
!!
!! SOURCE

subroutine rhoij_unpack(rhoij)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_unpack'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(pawrhoij_type),intent(inout) :: rhoij(:)

!Local variables-------------------------------
 integer :: natom,iat,lmn2_size,isel,klmn,nspden,cplex

! *************************************************************************

 natom  = SIZE(rhoij) ; if (natom==0) return
 nspden = rhoij(1)%nspden    ! MT jan-2010: this should not be nspden but nsppol or 4 if nspden=4
 cplex  = rhoij(1)%cplex

 do iat=1,natom

   lmn2_size =rhoij(iat)%lmn2_size

   if (rhoij(iat)%use_rhoij_/=1) then ! Have to allocate rhoij
     ABI_ALLOCATE(rhoij(iat)%rhoij_,(cplex*lmn2_size,nspden))
     rhoij(iat)%use_rhoij_=1
   end if
   rhoij(iat)%rhoij_ = zero

   do isel=1,rhoij(iat)%nrhoijsel ! Looping over non-zero ij elements.
     klmn = rhoij(iat)%rhoijselect(isel)
     rhoij(iat)%rhoij_(klmn,:) = rhoij(iat)%rhoijp(isel,:)
   end do

 end do ! natom

end subroutine rhoij_unpack
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/rhoij_init_unpacked
!! NAME
!! rhoij_init_unpacked
!!
!! FUNCTION
!!  Initialize field of rhoij datastructure for unpacked values (pawrhoij%rhoij_ array)
!!
!! SIDE EFFECTS
!!  rhoij(:) <pawrhoij_type)>= input/output datastructure
!!   * In output the rhoij_ array is allocated
!!
!! PARENTS
!!      energy,pawmkrhoij,rhofermi3,vtorho3
!!
!! CHILDREN
!!      xsum_mpi
!!
!! SOURCE

subroutine rhoij_init_unpacked(rhoij)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_init_unpacked'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(pawrhoij_type),intent(inout) :: rhoij(:)

!Local variables-------------------------------
 integer :: iat,nrhoij,nsp2

! *************************************************************************

 nrhoij=SIZE(rhoij);if (nrhoij==0) return
 nsp2=rhoij(1)%nsppol;if (rhoij(1)%nspden==4) nsp2=4

 do iat=1,nrhoij

   if (associated(rhoij(iat)%rhoij_))  then
     ABI_DEALLOCATE(rhoij(iat)%rhoij_)
   end if
   ABI_ALLOCATE(rhoij(iat)%rhoij_,(rhoij(iat)%cplex*rhoij(iat)%lmn2_size,nsp2))
   rhoij(iat)%use_rhoij_=1
   rhoij(iat)%rhoij_=zero

 end do

end subroutine rhoij_init_unpacked
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/rhoij_destroy_unpacked
!! NAME
!! rhoij_destroy_unpacked
!!
!! FUNCTION
!!  Destroy field of rhoij datastructure for unpacked values (pawrhoij%rhoij_ array)
!!
!! SIDE EFFECTS
!!  rhoij(:) <pawrhoij_type)>= input/output datastructure
!!   * In output the rhoij_ array is deallocated
!!
!! PARENTS
!!      pawmkrho
!!
!! CHILDREN
!!      xsum_mpi
!!
!! SOURCE

subroutine rhoij_destroy_unpacked(rhoij)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_destroy_unpacked'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(pawrhoij_type),intent(inout) :: rhoij(:)

!Local variables-------------------------------
 integer :: iat,nrhoij

! *************************************************************************

 nrhoij=SIZE(rhoij);if (nrhoij==0) return

 do iat=1,nrhoij

   if (associated(rhoij(iat)%rhoij_))  then
     ABI_DEALLOCATE(rhoij(iat)%rhoij_)
   end if
   nullify(rhoij(iat)%rhoij_)
   rhoij(iat)%use_rhoij_=0

 end do

end subroutine rhoij_destroy_unpacked
!!***

!----------------------------------------------------------------------

!!****f* m_pawrhoij/rhoij_mpi_sum
!! NAME
!! rhoij_mpi_sum
!!
!! FUNCTION
!! Build the MPI sum of the unsymmetrized PAW rhoij_ (augmentation occupancies)
!! Remember:for each atom, rho_ij=Sum_{n,k} {occ(n,k)*<Cnk|p_i><p_j|Cnk>}
!!
!! INPUTS
!!  comm1=MPI communicator. Data will be MPI summed inside comm1
!!  [comm2]=second MPI communicator. If present, rhoij_ will be MPI summed inside comm2 after the collective sum in comm1.
!!
!! SIDE EFFECTS
!!  pawrhoij(:) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  Input: the data calculateed by this processor.
!!  Otput: the final MPI sum over comm1 and comm2.
!!
!! PARENTS
!!      energy,pawmkrhoij,wfd_pawrhoij
!!
!! CHILDREN
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

#include "abi_common_for_bigdft.h"

subroutine rhoij_mpi_sum(pawrhoij,comm1,comm2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhoij_mpi_sum'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: comm1
 integer,optional,intent(in) :: comm2
!arrays
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables ---------------------------------------
!scalars
 integer :: bufdim,iatom,ierr,isppol,jdim,nsp2,natom
 integer :: nproc1,nproc2
 !character(len=500) :: msg
!arrays
 integer,allocatable :: dimlmn(:)
 real(dp),allocatable :: buffer1(:),buffer2(:)

!************************************************************************

! DBG_ENTER("COLL")

 natom=SIZE(pawrhoij);if (natom==0) return

 nproc1 = xcomm_size(comm1)
 nproc2=1; if (PRESENT(comm2)) nproc2 = xcomm_size(comm2)
 if (nproc1==1.and.nproc2==1) RETURN

!Fill the MPI buffer from the local rhoij_
 ABI_ALLOCATE(dimlmn,(natom))
 dimlmn(1:natom)=pawrhoij(1:natom)%cplex*pawrhoij(1:natom)%lmn2_size
 nsp2=pawrhoij(1)%nsppol; if (pawrhoij(1)%nspden==4) nsp2=4
 bufdim=sum(dimlmn)*nsp2
 ABI_ALLOCATE(buffer1,(bufdim))
 ABI_ALLOCATE(buffer2,(bufdim))
 jdim=0
 do iatom=1,natom
   do isppol=1,nsp2
     buffer1(jdim+1:jdim+dimlmn(iatom))=pawrhoij(iatom)%rhoij_(:,isppol)
     jdim=jdim+dimlmn(iatom)
   end do
 end do
!
!Build sum of  pawrhoij%rhoij_
 call xsum_mpi(buffer1,buffer2,bufdim,comm1,ierr)      ! Sum over the first communicator.
 if (PRESENT(comm2)) call xsum_mpi(buffer2,comm2,ierr) ! Sum over the second communicator.
!
!Unpack the MPI packet filling rhoij_.
 jdim=0
 do iatom=1,natom
   do isppol=1,nsp2
     pawrhoij(iatom)%rhoij_(:,isppol)=buffer2(jdim+1:jdim+dimlmn(iatom))
     jdim=jdim+dimlmn(iatom)
   end do
 end do

 ABI_DEALLOCATE(buffer1)
 ABI_DEALLOCATE(buffer2)
 ABI_DEALLOCATE(dimlmn)

! DBG_EXIT("COLL")

end subroutine rhoij_mpi_sum
!!***

!----------------------------------------------------------------------

END MODULE m_pawrhoij
!!***
