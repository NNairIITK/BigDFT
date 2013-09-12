!> @file
!!  Module to handle the fragments of a system
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 
!> Module which defines some important structures and methods to manipulate embedding systems
module module_fragments
  use module_base, only: gp,wp
  use module_types
  use dynamic_memory
  implicit none

  private

  !interface operator(*)
  !   module procedure transform_fragment
  !end interface

   interface
      subroutine read_atomic_file(file,iproc,astruct,status,comment,energy,fxyz)
         !n(c) use module_base
         use module_types
         implicit none
         character(len=*), intent(in) :: file
         integer, intent(in) :: iproc
         type(atomic_structure), intent(inout) :: astruct
         integer, intent(out), optional :: status
         real(gp), intent(out), optional :: energy
         real(gp), dimension(:,:), pointer, optional :: fxyz
         character(len =*), intent(out), optional :: comment
      END SUBROUTINE read_atomic_file
   end interface


  !> information about the basis set metadata shared between cubic and ASF approaches
  type, public :: minimal_orbitals_data
     integer :: norb          !< Total number of orbitals per k point
     integer :: norbp         !< Total number of orbitals for the given processors
     integer :: isorb         !< Total number of orbitals for the given processors
     integer, dimension(:), pointer :: inwhichlocreg,onwhichatom !< associate the basis centers
     integer, dimension(:), pointer :: isorb_par,ispot
     integer, dimension(:,:), pointer :: norb_par
  end type minimal_orbitals_data

  type, public :: phi_array
     real(wp), dimension(:,:,:,:,:,:), pointer :: psig
  end type phi_array

  type, public :: fragment_basis
     integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
     integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
     type(local_zone_descriptors) :: Lzd
     type(minimal_orbitals_data) :: forbs
     real(wp), dimension(:), pointer :: psi_full
     type(phi_array), dimension(:), pointer :: phi
  end type fragment_basis

  type, public :: fragment_orbitals
     integer :: norb, nphi

  end type fragment_orbitals

  !> Defines the minimal information to identify a system building block
  type, public :: system_fragment
     integer :: nat_env !< environment atoms which complete fragment specifications
     real(gp), dimension(:,:), pointer :: rxyz_env !< position of atoms in environment (AU), external reference frame
     type(atomic_structure) :: astruct_frg !< Number of atoms, positions, atom type etc for fragment
     type(fragment_basis) :: fbasis !< fragment basis, associated only if coherent with positions, pointer - do we really want this to be a pointer?
     ! add coeffs and or kernel
     integer :: nksorb
     real(gp), dimension(:,:), pointer :: coeff    
  end type system_fragment

  !> Contains the rotation and translation (possibly deformation) which have to be applied to a given fragment
  type, public :: fragment_transformation
     !real(gp), dimension(3) ::  !< translation of fragment center
     real(gp), dimension(3) :: rot_center_new !< positions of the centers
     real(gp), dimension(3) :: rot_center !< center of rotation in original coordinates (might be fragment center or not)
     real(gp), dimension(3) :: rot_axis !< unit rotation axis (should have modulus one)
     real(gp) :: theta !< angle of rotation 
     character(len=2), dimension(:), pointer :: discrete_operations
  end type fragment_transformation

  !public operator(*)

  public :: fragment_null, fragment_free, init_fragments, minimal_orbitals_data_null, frag_center, find_frag_trans

contains

  ! nned to check somewhere if fragment linking is consistent - initialize linear from file?!
  ! initializes reference fragments (already nullified), if it isn't a fragment calculation sets to appropriate dummy values
  ! ignoring environment for now
  subroutine init_fragments(in,orbs,astruct,ref_frags)
    use module_types
    implicit none
    type(input_variables), intent(in) :: in
    type(orbitals_data), intent(in) :: orbs ! orbitals of full system, needed to set 'dummy' values
    type(atomic_structure), intent(in) :: astruct ! atomic structure of full system, needed to set 'dummy' values
    type(system_fragment), dimension(in%frag%nfrag_ref), intent(inout) :: ref_frags

    ! local variables
    integer :: ifrag

    if (in%lin%fragment_calculation) then
        ! read fragment posinps and initialize fragment, except for psi and lzds
        do ifrag=1,in%frag%nfrag_ref
           call init_fragment_from_file(ref_frags(ifrag),trim(in%dir_output)//trim(in%frag%label(ifrag)),in)
        end do

        ! check that fragments are sensible, i.e. correct number of atoms, atom types etc.
        call check_fragments(in,ref_frags,astruct)

     ! set appropriate default values as this is not a fragment calculation
     else
        ! nullify fragment
        ref_frags(1)=fragment_null()

        ! want all components of ref_frags(1)%fbasis%forbs to point towards appropriate component of orbs of full system
        call orbs_to_min_orbs_point(orbs,ref_frags(1)%fbasis%forbs)

        ! environment and coeffs
        call fragment_allocate(ref_frags(1))

        ! rest of fbasis isn't needed (I think!) so can remain nullified
  
        ! astruct - fill in other bits later
        ref_frags(1)%astruct_frg%nat=astruct%nat

     end if

  end subroutine init_fragments


  !> Initializes all of fragment except lzd using the fragment posinp and tmb files
  subroutine init_fragment_from_file(frag,frag_name,input) ! switch this to pure if possible
    use module_types
    !use module_interfaces
    implicit none
    type(system_fragment), intent(inout) :: frag
    character(len=*), intent(in) :: frag_name
    type(input_variables), intent(in) :: input

    ! local variables

    ! nullify fragment
    frag=fragment_null()

    ! read fragment positions
    call read_atomic_file(frag_name(1:len(frag_name)),bigdft_mpi%iproc,frag%astruct_frg)

    ! iproc, nproc, nspinor not needed yet, add in later
    call init_minimal_orbitals_data(bigdft_mpi%iproc, bigdft_mpi%nproc, 1, input, frag%astruct_frg, frag%fbasis%forbs)
    !call init_minimal_orbitals_data(iproc, nproc, nspinor, input, frag%astruct_frg, frag%fbasis%forbs)

    ! environment and coeffs
    call fragment_allocate(frag)

    ! allocate/initialize fragment basis...
    !currently orbitals are initialized via initAndUtils/init_orbitals_data_for_linear
    !which calls wavefunctions/orbitals_descriptors to assign parallel bits and superfluous stuff
    !and locreg_orbitals/assignToLocreg2 to give inwhichlocreg and onwhichatom - but we should really take onwhichatom from file?!

    !type, public :: fragment_basis
    !   integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
    !   integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
    !   type(local_zone_descriptors) :: Lzd
    !   type(minimal_orbitals_data) :: forbs

  end subroutine init_fragment_from_file


  ! sanity check on fragment definitions
  subroutine check_fragments(input,ref_frags,astruct)
    use module_types
    implicit none
    type(input_variables), intent(in) :: input
    type(atomic_structure), intent(in) :: astruct ! atomic structure of full system, needed to check fragments are sensible
    type(system_fragment), dimension(input%frag%nfrag_ref), intent(in) :: ref_frags

    integer :: ifrag, ifrag_ref, tot_frag_ats, isfat, iat
    logical :: fragments_ok

    fragments_ok=.true.

    tot_frag_ats=0
    do ifrag=1,input%frag%nfrag
       ifrag_ref = input%frag%frag_index(ifrag)
       tot_frag_ats = tot_frag_ats + ref_frags(ifrag_ref)%astruct_frg%nat
    end do

    if (tot_frag_ats /= astruct%nat) then
       fragments_ok=.false.
       write(*,*) 'Sum of atoms in fragments not equal to total atoms in system',tot_frag_ats,astruct%nat
    end if

    isfat=0
    do ifrag=1,input%frag%nfrag
       ifrag_ref=input%frag%frag_index(ifrag)

       do iat=1,ref_frags(ifrag_ref)%astruct_frg%nat
          if (ref_frags(ifrag_ref)%astruct_frg%atomnames(ref_frags(ifrag_ref)%astruct_frg%iatype(iat)) &
               /= astruct%atomnames(astruct%iatype(iat+isfat))) then
             fragments_ok=.false.
             write(*,*) 'Atom type for fragment ',ifrag,', reference fragment ',ifrag_ref,' atom number ',iat,&
                  ' does not match ',ref_frags(ifrag_ref)%astruct_frg%atomnames(ref_frags(ifrag_ref)%astruct_frg%iatype(iat)),&
                  astruct%atomnames(astruct%iatype(iat+isfat))
          end if
       end do

       isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat     
    end do

    if (.not. fragments_ok) stop 'Problem in fragment definitions'

  end subroutine check_fragments


    !type, public :: minimal_orbitals_data
    !   integer :: norb          !< Total number of orbitals per k point
    !   integer :: norbp         !< Total number of orbitals for the given processors
    !   integer :: isorb         !< Total number of orbitals for the given processors
    !   integer, dimension(:), pointer :: inwhichlocreg,onwhichatom !< associate the basis centers
    !   integer, dimension(:), pointer :: isorb_par,ispot
    !   integer, dimension(:,:), pointer :: norb_par
  !> just initializing norb for now, come back and do the rest later
  subroutine init_minimal_orbitals_data(iproc, nproc, nspinor, input, astruct, forbs)
    use module_base
    use module_types
    implicit none
  
    ! Calling arguments
    integer,intent(in) :: iproc, nproc, nspinor
    type(input_variables),intent(in) :: input
    type(atomic_structure),intent(in) :: astruct
    type(minimal_orbitals_data),intent(out) :: forbs
  
    ! Local variables
    integer :: norb, ityp, iat
    !integer :: norbu, norbd, nlr, ilr, iall, iorb, istat
    !integer,dimension(:),allocatable :: norbsPerLocreg, norbsPerAtom
    !real(kind=8),dimension(:,:),allocatable :: locregCenter
    !character(len=*),parameter :: subname='init_minimal_orbitals_data'

    call timing(iproc,'init_orbs_lin ','ON')
 
    ! Count the number of basis functions.
    !allocate(norbsPerAtom(astruct%nat), stat=istat)
    !call memocc(istat, norbsPerAtom, 'norbsPerAtom', subname)
    norb=0
    !nlr=0
    do iat=1,astruct%nat
       ityp=astruct%iatype(iat)
       !norbsPerAtom(iat)=input%lin%norbsPerType(ityp)
       norb=norb+input%lin%norbsPerType(ityp)
       !nlr=nlr+input%lin%norbsPerType(ityp)
    end do

    ! Distribute the basis functions among the processors.
    !norbu=norb
    !norbd=0
 
    forbs%norb=norb

    !call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, input%nspin, nspinor,&
    !     input%nkpt, input%kpt, input%wkpt, lorbs,.true.) !simple repartition
 

    !allocate(locregCenter(3,nlr), stat=istat)
    !call memocc(istat, locregCenter, 'locregCenter', subname)
  
    !ilr=0
    !do iat=1,astruct%nat
    !    ityp=astruct%iatype(iat)
    !    do iorb=1,input%lin%norbsPerType(ityp)
    !        ilr=ilr+1
    !        locregCenter(:,ilr)=astruct%rxyz(:,iat)
    !        ! DEBUGLR write(10,*) iorb,locregCenter(:,ilr)
    !    end do
    !end do
 
    !allocate(norbsPerLocreg(nlr), stat=istat)
    !call memocc(istat, norbsPerLocreg, 'norbsPerLocreg', subname)
    !norbsPerLocreg=1 !should be norbsPerLocreg
    
    !iall=-product(shape(forbs%inWhichLocreg))*kind(forbs%inWhichLocreg)
    !deallocate(forbs%inWhichLocreg, stat=istat)
    !call memocc(istat, iall, 'forbs%inWhichLocreg', subname)
    !call assignToLocreg2(iproc, nproc, forbs%norb, forbs%norb_par, astruct%nat, nlr, &
    !     input%nspin, norbsPerLocreg, locregCenter, forbs%inwhichlocreg)

    !iall=-product(shape(forbs%onwhichatom))*kind(forbs%onwhichatom)
    !deallocate(forbs%onwhichatom, stat=istat)
    !call memocc(istat, iall, 'forbs%onwhichatom', subname)
    !call assignToLocreg2(iproc, nproc, forbs%norb, forbs%norb_par, astruct%nat, astruct%nat, &
    !     input%nspin, norbsPerAtom, astruct%rxyz, forbs%onwhichatom)
  

    !iall=-product(shape(norbsPerLocreg))*kind(norbsPerLocreg)
    !deallocate(norbsPerLocreg, stat=istat)
    !call memocc(istat, iall, 'norbsPerLocreg', subname)
  
    !iall=-product(shape(locregCenter))*kind(locregCenter)
    !deallocate(locregCenter, stat=istat)
    !call memocc(istat, iall, 'locregCenter', subname)

    !iall=-product(shape(norbsPerAtom))*kind(norbsPerAtom)
    !deallocate(norbsPerAtom, stat=istat)
    !call memocc(istat, iall, 'norbsPerAtom', subname)


    call timing(iproc,'init_orbs_lin ','OF')

  end subroutine init_minimal_orbitals_data

  ! copy a given full orbs structure to minimal orbs structure
  subroutine orbs_to_min_orbs_copy(orbs,forbs)
    use module_types
    implicit none
    ! Calling arguments
    type(orbitals_data),intent(in):: orbs
    type(minimal_orbitals_data),intent(inout):: forbs

    ! Local variables
    integer:: iis1, iie1, iis2, iie2, i1, i2, istat, iall
    character(len=200) :: subname

    subname='orbs_to_min_orbs_copy'

    forbs%norb = orbs%norb
    forbs%norbp = orbs%norbp
    forbs%isorb = orbs%isorb

    if(associated(forbs%norb_par)) then
        iall=-product(shape(forbs%norb_par))*kind(forbs%norb_par)
        deallocate(forbs%norb_par, stat=istat)
        call memocc(istat, iall, 'forbs%norb_par', subname)
    end if
    if(associated(orbs%norb_par)) then
        iis1=lbound(orbs%norb_par,1)
        iie1=ubound(orbs%norb_par,1)
        iis2=lbound(orbs%norb_par,2)
        iie2=ubound(orbs%norb_par,2)
        allocate(forbs%norb_par(iis1:iie1,iis2:iie2), stat=istat)
        call memocc(istat, forbs%norb_par, 'forbs%norb_par', subname)
        do i1=iis1,iie1
           do i2 = iis2,iie2
            forbs%norb_par(i1,i2) = orbs%norb_par(i1,i2)
           end do
        end do
    end if

    if(associated(forbs%inwhichlocreg)) then
        iall=-product(shape(forbs%inwhichlocreg))*kind(forbs%inwhichlocreg)
        deallocate(forbs%inwhichlocreg, stat=istat)
        call memocc(istat, iall, 'forbs%inwhichlocreg', subname)
    end if
    if(associated(orbs%inwhichlocreg)) then
        iis1=lbound(orbs%inwhichlocreg,1)
        iie1=ubound(orbs%inwhichlocreg,1)
        allocate(forbs%inwhichlocreg(iis1:iie1), stat=istat)
        call memocc(istat, forbs%inwhichlocreg, 'forbs%inwhichlocreg', subname)
        do i1=iis1,iie1
            forbs%inwhichlocreg(i1) = orbs%inwhichlocreg(i1)
        end do
    end if

    if(associated(forbs%onwhichatom)) then
        iall=-product(shape(forbs%onwhichatom))*kind(forbs%onwhichatom)
        deallocate(forbs%onwhichatom, stat=istat)
        call memocc(istat, iall, 'forbs%onwhichatom', subname)
    end if
    if(associated(orbs%onwhichatom)) then
        iis1=lbound(orbs%onwhichatom,1)
        iie1=ubound(orbs%onwhichatom,1)
        allocate(forbs%onwhichatom(iis1:iie1), stat=istat)
        call memocc(istat, forbs%onwhichatom, 'forbs%onwhichatom', subname)
        do i1=iis1,iie1
            forbs%onwhichatom(i1) = orbs%onwhichatom(i1)
        end do
    end if

    if(associated(forbs%isorb_par)) then
        iall=-product(shape(forbs%isorb_par))*kind(forbs%isorb_par)
        deallocate(forbs%isorb_par, stat=istat)
        call memocc(istat, iall, 'forbs%isorb_par', subname)
    end if
    if(associated(orbs%isorb_par)) then
        iis1=lbound(orbs%isorb_par,1)
        iie1=ubound(orbs%isorb_par,1)
        allocate(forbs%isorb_par(iis1:iie1), stat=istat)
        call memocc(istat, forbs%isorb_par, 'forbs%isorb_par', subname)
        do i1=iis1,iie1
            forbs%isorb_par(i1) = orbs%isorb_par(i1)
        end do
    end if

    if(associated(forbs%ispot)) then
        iall=-product(shape(forbs%ispot))*kind(forbs%ispot)
        deallocate(forbs%ispot, stat=istat)
        call memocc(istat, iall, 'forbs%ispot', subname)
    end if
    if(associated(orbs%ispot)) then
        iis1=lbound(orbs%ispot,1)
        iie1=ubound(orbs%ispot,1)
        allocate(forbs%ispot(iis1:iie1), stat=istat)
        call memocc(istat, forbs%ispot, 'forbs%ispot', subname)
        do i1=iis1,iie1
            forbs%ispot(i1) = orbs%ispot(i1)
        end do
    end if

  end subroutine orbs_to_min_orbs_copy

  ! point minimal orbs structure to a given full orbs structure
  subroutine orbs_to_min_orbs_point(orbs,forbs)
    use module_types
    implicit none
    ! Calling arguments
    type(orbitals_data),intent(in):: orbs
    type(minimal_orbitals_data),intent(inout):: forbs

    forbs%norb = orbs%norb
    forbs%norbp = orbs%norbp
    forbs%isorb = orbs%isorb

    forbs%norb_par => orbs%norb_par
    forbs%inwhichlocreg => orbs%inwhichlocreg
    forbs%onwhichatom => orbs%onwhichatom
    forbs%isorb_par => orbs%isorb_par
    forbs%ispot => orbs%ispot
  end subroutine orbs_to_min_orbs_point


  pure function minimal_orbitals_data_null() result(forbs)
    implicit none
    type(minimal_orbitals_data) :: forbs
    call nullify_minimal_orbitals_data(forbs)
  end function minimal_orbitals_data_null

  pure subroutine nullify_minimal_orbitals_data(forbs)
    implicit none
    type(minimal_orbitals_data), intent(out) :: forbs

    forbs%norb=0
    forbs%norbp=0
    forbs%isorb=0
    nullify(forbs%inwhichlocreg)
    nullify(forbs%onwhichatom)
    nullify(forbs%isorb_par)
    nullify(forbs%ispot)
    nullify(forbs%norb_par)
  end subroutine nullify_minimal_orbitals_data


  pure function fragment_null() result(frag)
    implicit none
    type(system_fragment) :: frag

    frag%nat_env=0
    nullify(frag%rxyz_env)
    frag%nksorb=0
    nullify(frag%coeff)
    call nullify_atomic_structure(frag%astruct_frg)
    ! nullify fragment basis
    call nullify_fragment_basis(frag%fbasis)

  end function fragment_null

  pure function fragment_basis_null() result(basis)
    implicit none
    type(fragment_basis) :: basis
    call nullify_fragment_basis(basis)
  end function fragment_basis_null

  pure subroutine nullify_fragment_basis(basis)
    implicit none
    type(fragment_basis), intent(out) :: basis

    basis%npsidim_orbs=0
    basis%npsidim_comp=0
    call nullify_local_zone_descriptors(basis%lzd)
    call nullify_minimal_orbitals_data(basis%forbs)
    !basis%forbs=minimal_orbitals_data_null()
    nullify(basis%psi_full)
    nullify(basis%phi)
  end subroutine nullify_fragment_basis

  subroutine minimal_orbitals_data_free(forbs)
    implicit none
    type(minimal_orbitals_data), intent(inout) :: forbs

    if (associated(forbs%inwhichlocreg)) call f_free_ptr(forbs%inwhichlocreg)
    if (associated(forbs%onwhichatom)) call f_free_ptr(forbs%onwhichatom)
    if (associated(forbs%isorb_par)) call f_free_ptr(forbs%isorb_par)
    if (associated(forbs%ispot)) call f_free_ptr(forbs%ispot)
    if (associated(forbs%norb_par)) call f_free_ptr(forbs%norb_par)
    forbs=minimal_orbitals_data_null()
  end subroutine minimal_orbitals_data_free

  subroutine fragment_basis_free(basis)
    implicit none
    type(fragment_basis), intent(inout) :: basis
    character(len=200) :: subname

    subname='fragment_basis_free'
    call deallocate_local_zone_descriptors(basis%lzd,subname)
    call minimal_orbitals_data_free(basis%forbs)
    if (associated(basis%psi_full)) call f_free_ptr(basis%psi_full)
    basis=fragment_basis_null()
  end subroutine fragment_basis_free

  subroutine fragment_free(frag)
    implicit none
    type(system_fragment), intent(inout) :: frag
    character(len=200) :: subname
    integer :: i_all, i_stat

    subname='fragment_free'

    call deallocate_atomic_structure(frag%astruct_frg,subname)
    frag%astruct_frg=atomic_structure_null()
    call f_routine(id='fragment_free')
    if (associated(frag%rxyz_env)) call f_free_ptr(frag%rxyz_env)
    if (associated(frag%coeff)) call f_free_ptr(frag%coeff)
    call f_release_routine()
    call fragment_basis_free(frag%fbasis)
    frag=fragment_null()

  end subroutine fragment_free

  subroutine fragment_allocate(frag)
    implicit none
    type(system_fragment), intent(inout) :: frag

    call f_routine(id='fragment_allocate')

    frag%rxyz_env=f_malloc_ptr((/3,frag%nat_env/),id='frag%rxyz_env')
    frag%coeff=f_malloc_ptr((/frag%fbasis%forbs%norb,frag%fbasis%forbs%norb/),id='frag%coeff')

    call f_release_routine()

  end subroutine fragment_allocate


  pure function frg_center(frag)
    implicit none
    type(system_fragment), intent(in) :: frag
    real(gp), dimension(3) :: frg_center
    !local variables
    integer :: iat

    frg_center=0.0_gp
    do iat=1,frag%astruct_frg%nat
       frg_center=frg_center+frag%astruct_frg%rxyz(:,iat)
    end do
    frg_center=frg_center/real(frag%astruct_frg%nat,gp)

  end function frg_center

  !function transform_fragment(trans,frag) result(frag_new)
  !  implicit none
  !  type(fragment_transformation), intent(in) :: trans
  !  type(system_fragment), intent(in) :: frag
  !  type(system_fragment) :: frag_new
  !
  !  ! local variables
  !  integer :: iat
  !
  !  frag_new=fragment_null()
  !  frag_new%astruct_frg%nat=frag%astruct_frg%nat
  !  frag_new%nat_env=frag%nat_env
  !  frag_new%astruct_frg%ntypes=frag%astruct_frg%ntypes
  !
  !  ! allocate arrays here, leave fragment_basis nullified
  !  call fragment_allocate(frag_new)
  !
  !  ! do fragment first, then environment
  !  do iat=1,frag%astruct_frg%nat
  !     frag_new%astruct_frg%rxyz(:,iat)=rotate_vector(trans%rot_axis,&
  !          trans%theta,frag%astruct_frg%rxyz(:,iat)-trans%rot_center(:))
  !     frag_new%astruct_frg%rxyz(:,iat)=frag_new%astruct_frg%rxyz(:,iat)+trans%rot_center(:)+trans%dr(:)
  !  end do
  !
  !  do iat=1,frag%nat_env
  !     frag_new%rxyz_env(:,iat)=rotate_vector(trans%rot_axis,trans%theta,frag%rxyz_env(:,iat)-trans%rot_center(:))
  !     frag_new%rxyz_env(:,iat)=frag_new%rxyz_env(:,iat)+trans%rot_center(:)+trans%dr(:)
  !  end do
  !
  !  ! to complete should copy across iatype and atomnames but only using as a check to see how effective trans is
  !
  !end function transform_fragment



  !> Express the coordinates of a vector into a rotated reference frame
  pure function rotate_vector(newz,theta,vec) result(vecn)
     use module_base
     implicit none
     real(gp), intent(in) :: theta
     real(gp), dimension(3), intent(in) :: newz,vec
     real(gp), dimension(3) :: vecn
     !local variables
     real(gp) :: sint,cost,onemc,x,y,z

     !save recalculation
     sint=sin(theta)
     cost=cos(theta)
     onemc=1.0_gp-cost
     x=vec(1)
     y=vec(2)
     z=vec(3)

     vecn(1)=x*(cost + onemc*newz(1)**2) + y*(onemc*newz(1)*newz(2) - sint*newz(3)) &
          + z*(sint*newz(2) + onemc*newz(1)*newz(3))
     vecn(2)=y*(cost + onemc*newz(2)**2) + x*(onemc*newz(1)*newz(2) + sint*newz(3)) &
          + z*(-(sint*newz(1)) + onemc*newz(2)*newz(3))
     vecn(3)=z*(cost + onemc*newz(3)**2) + x*(onemc*newz(1)*newz(3) - sint*newz(2)) &
          + y*(sint*newz(1) + onemc*newz(2)*newz(3))

  end function rotate_vector

  !pure function transform_fragment_basis(trans,basis) result(basis_new)
  !  implicit none
  !  type(fragment_transformation), intent(in) :: trans
  !  type(fragment_basis), intent(in) :: basis
  !  type(fragment_basis) :: basis_new

  !  basis_new=fragment_basis_null()

  !  ! minimal_orbitals_data should remain the same
  !! want to have defined new lzd already?  in which case not a pure function...

  !!   integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
  !!   integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
  !!   type(local_zone_descriptors) :: Lzd
  !!   type(minimal_orbitals_data) :: forbs
  !!   real(wp), dimension(:), pointer :: psi

  !end function transform_fragment_basis

  subroutine find_frag_trans(nat,rxyz_ref,rxyz_new,frag_trans)
    use module_base
    use yaml_output
    implicit none
    integer, intent(in) :: nat !< fragment size
    real(gp), dimension(3,nat), intent(in) :: rxyz_ref,rxyz_new !<coordinates measured wrt rot_center
    type(fragment_transformation), intent(inout) :: frag_trans
    !local variables
    integer, parameter :: lwork=7*3
    integer :: info,iat,i_stat!,i
    real(gp) :: dets,J
    real(gp), dimension(3) :: SM_arr !< array of SVD and M array
    real(gp), dimension(lwork) :: work !< array of SVD and M array
    real(gp), dimension(3,nat) :: J_arr !< matrix for calculating Wahba's cost function
    real(gp), dimension(3,3) :: B_mat,R_mat,U_mat,VT_mat !<matrices of Wahba's problem
    character(len=100) :: subname

    subname='find_frag_trans'

    B_mat=0.0_gp
    R_mat=0.0_gp

    !all positions are of weight one for the moment
    call dgemm('N','T',3,3,nat,1.0_gp,rxyz_new,3,rxyz_ref,3,0.0_gp,B_mat,3)

    !find matrix of svd
    call dgesvd('A','A',3,3,B_mat,3,SM_arr,U_mat,3,VT_mat,3,work,lwork,info)
    if (f_err_raise(info/=0,'Problem in DGESVD')) return
    
    !multiply last line of VT_mat by det(U)*det(V)
    dets=det_33(U_mat)*det_33(VT_mat)
    VT_mat(3,:)=VT_mat(3,:)*dets

    !find rotation matrix
    call dgemm('N','N',3,3,3,1.0_gp,U_mat,3,VT_mat,3,0.0_gp,R_mat,3)

    !find the angle from R matrix
    frag_trans%theta=theta_from_r(R_mat)

    !find rot_axis
    frag_trans%rot_axis=axis_from_r(R_mat)

    !print*,'Rmat:',frag_trans%theta
    !do i=1,3
    !   write(*,'(3(F12.6,2x))') R_mat(i,:)
    !end do

    !print*,'Rcalc:',frag_trans%theta
    !write(*,'(3(F12.6,2x))') cos(frag_trans%theta)+frag_trans%rot_axis(1)**2*(1.0_gp-cos(frag_trans%theta)),&
    !     frag_trans%rot_axis(1)*frag_trans%rot_axis(2)*(1.0_gp-cos(frag_trans%theta))&
    !     -frag_trans%rot_axis(3)*sin(frag_trans%theta),&
    !     frag_trans%rot_axis(1)*frag_trans%rot_axis(3)*(1.0_gp-cos(frag_trans%theta))&
    !     +frag_trans%rot_axis(2)*sin(frag_trans%theta)
    !write(*,'(3(F12.6,2x))') frag_trans%rot_axis(2)*frag_trans%rot_axis(1)*(1.0_gp-cos(frag_trans%theta))&
    !     +frag_trans%rot_axis(3)*sin(frag_trans%theta),&
    !     cos(frag_trans%theta)+frag_trans%rot_axis(2)**2*(1.0_gp-cos(frag_trans%theta)),&
    !     frag_trans%rot_axis(2)*frag_trans%rot_axis(3)*(1.0_gp-cos(frag_trans%theta))&
    !     -frag_trans%rot_axis(1)*sin(frag_trans%theta)
    !write(*,'(3(F12.6,2x))') frag_trans%rot_axis(3)*frag_trans%rot_axis(1)*(1.0_gp-cos(frag_trans%theta))&
    !     -frag_trans%rot_axis(2)*sin(frag_trans%theta),&
    !     frag_trans%rot_axis(3)*frag_trans%rot_axis(2)*(1.0_gp-cos(frag_trans%theta))&
    !     +frag_trans%rot_axis(1)*sin(frag_trans%theta),&
    !     cos(frag_trans%theta)+frag_trans%rot_axis(3)**2*(1.0_gp-cos(frag_trans%theta))
    
    !to be verified that the cost function of Wahba's problem is little
    J_arr=rxyz_new
    call dgemm('N','N',3,nat,3,-1.0_gp,R_mat,3,rxyz_ref,3,1.0_gp,J_arr,3)
    J=0.0_gp
    do iat=1,nat
       J=J+J_arr(1,iat)**2+J_arr(2,iat)**2+J_arr(3,iat)**2
    end do

    if (J>1.0e-4) then
       print*,"Error, Wahba's cost function is too big",J,frag_trans%theta/(4.0_gp*atan(1.d0)/180.0_gp)
    end if

    !check the pertinence of the suggested rotation
    !if (abs(frag_trans%theta) > 60.d0*(4.0_gp*atan(1.d0)/180.0_gp)) print*,'frag_trans%theta=',frag_trans%theta/(4.0_gp*atan(1.d0)/180.0_gp)
	 !if  (f_err_raise(abs(frag_trans%theta) > 60.d0*(4.0_gp*atan(1.d0)/180.0_gp),'Angle frag_trans%theta not optimal (frag_trans%theta= '//&
	 !      yaml_toa(frag_trans%theta)//' )')) return

    ! reduce the angle if frag_trans%theta > 60 degrees
    !if (abs(frag_trans%theta) > 60.d0*(4.0_gp*atan(1.d0)/180.0_gp)) then
    !   write(*,*) 'before',frag_trans%rot_axis,frag_trans%theta/(4.0_gp*atan(1.d0)/180.0_gp)
    !   call find_discrete_operations(frag_trans,R_mat)
    !   write(*,*) 'after',frag_trans%rot_axis,frag_trans%theta/(4.0_gp*atan(1.d0)/180.0_gp)
    !else
       allocate(frag_trans%discrete_operations(0),stat=i_stat)
       call memocc(i_stat,frag_trans%discrete_operations,'frag_trans%discrete_operations',subname)
    !end if


  end subroutine find_frag_trans



  subroutine find_discrete_operations(frag_trans,R_mat)
    use module_base
    use yaml_output
    use dictionaries
    implicit none
    type(fragment_transformation), intent(inout) :: frag_trans
    real(gp), dimension(3,3), intent(inout) :: R_mat
    !local variables
    character(len=*), parameter :: subname='find_discrete_operations'
    integer :: ival, nval, i_stat, ix, iy, iz, list_len_max
    real(gp) :: min_theta, theta_orig, new_theta
    type(dictionary), pointer :: list, list_tmp
    integer, dimension(3,2) :: axis_type

    call dict_init(list)

    theta_orig=frag_trans%theta
    min_theta=theta_orig
    list_len_max=4

    do ix=3,0,-1
       do iy=0,3
          do iz=0,3

             axis_type(1,1)=ix  !type of operation
             axis_type(2,1)=iy  !type of operation
             axis_type(3,1)=iz  !type of operation
             axis_type(1,2)=1  !axis of operation
             axis_type(2,2)=2  !axis of operation
             axis_type(3,2)=3  !axis of operation

             call apply_discrete_operations_to_matrix(frag_trans%rot_axis,frag_trans%theta,new_theta,&
                  R_mat,axis_type,list,.false.)
!print*,'ix,iy,iz',ix,iy,iz,new_theta/(4.0_gp*atan(1.d0)/180.0_gp),min_theta/(4.0_gp*atan(1.d0)/180.0_gp)
             !write(*,'(F8.4,1x,L2,1x)',advance='no') new_theta/(4.0_gp*atan(1.d0)/180.0_gp)
             if (new_theta<min_theta) min_theta=new_theta
             if (new_theta/(4.0_gp*atan(1.d0)/180.0_gp)<=60.0_gp) then
                call dict_init(list_tmp)
                call apply_discrete_operations_to_matrix(frag_trans%rot_axis,frag_trans%theta,new_theta,&
                     R_mat,axis_type,list,.true.)
                exit
                !nval=dict_len(list_tmp)
!print*,'dict',nval,list_len_max
                !if (nval<list_len_max) then
                !    call dict_free(list) 
                !    call dict_init(list)
                !    list_len_max=nval
                !    do ival=0,nval-1
                !       call add(list,list_tmp//ival)
                !    end do
                !end if
                !call dict_free(list_tmp) 
                !if (list_len_max==1) exit
             end if

          end do
          if (min_theta/(4.0_gp*atan(1.d0)/180.0_gp)<=60.0_gp) exit
       end do
       if (min_theta/(4.0_gp*atan(1.d0)/180.0_gp)<=60.0_gp) exit
    end do

    ! copy from dict to frag_trans structure
    nval=dict_len(list)
    allocate(frag_trans%discrete_operations(nval),stat=i_stat)
    call memocc(i_stat,frag_trans%discrete_operations,'frag_trans%discrete_operations',subname)

    do ival=0,nval-1
       frag_trans%discrete_operations(ival+1)=list//ival
       print*,'filling discrete_ops',ival+1,frag_trans%discrete_operations(ival+1)
    end do

   call dict_free(list) 

   write(*,*) 'Old theta, new theta',theta_orig/(4.0_gp*atan(1.d0)/180.0_gp), min_theta/(4.0_gp*atan(1.d0)/180.0_gp),nval

   frag_trans%theta=min_theta
   frag_trans%rot_axis=axis_from_r(R_mat)

  end subroutine find_discrete_operations


  subroutine apply_discrete_operations_to_matrix(rot_axis,theta,new_theta,R_mat,axis_type,list,make_list)
    use module_base
    use dictionaries
    implicit none

    real(gp), intent(in) :: theta
    real(gp), dimension(3), intent(in) :: rot_axis
    real(gp), intent(out) :: new_theta
    real(gp), dimension(3,3), intent(inout) :: R_mat
    integer, dimension(3,2), intent(in) :: axis_type
    type(dictionary), pointer :: list ! list of operations performed
    logical, intent(in) :: make_list ! adds operations to a list
    !local variables
    real(gp), dimension(3,3) :: R_mat_new, R_mat2
    integer :: i

    R_mat2=R_mat

    do i=1,3,2
       if (axis_type(i,1)==0) then
          ! identity operation, no need to add to list
          R_mat_new=R_mat2
       else if (axis_type(i,1)==1) then
          call subtract_90_rotation(R_mat2,R_mat_new,axis_type(i,2))
          if (axis_type(i,2)==1.and.make_list) call add(list,'x1')
          if (axis_type(i,2)==2.and.make_list) call add(list,'y1')
          if (axis_type(i,2)==3.and.make_list) call add(list,'z1')
       else if (axis_type(i,1)==2) then
          call subtract_180_rotation(R_mat2,R_mat_new,axis_type(i,2))
          if (axis_type(i,2)==1.and.make_list) call add(list,'x2')
          if (axis_type(i,2)==2.and.make_list) call add(list,'y2')
          if (axis_type(i,2)==3.and.make_list) call add(list,'z2')
       else if (axis_type(i,1)==3) then
          call add_90_rotation(R_mat2,R_mat_new,axis_type(i,2))
          if (axis_type(i,2)==1.and.make_list) call add(list,'x3')
          if (axis_type(i,2)==2.and.make_list) call add(list,'y3')
          if (axis_type(i,2)==3.and.make_list) call add(list,'z3')
       end if

       new_theta=theta_from_r(R_mat_new)
       !newz_new=axis_from_r(R_mat_new)
        !if (make_list) print*,'discrete op',i,axis_type(i,2),new_theta!,axis_from_r(R_mat_new)

       if (i==3) exit

       if (axis_type(i+1,1)==0) then
          ! identity operation, no need to add to list
          R_mat2=R_mat_new
       else if (axis_type(i+1,1)==1) then
          call subtract_90_rotation(R_mat_new,R_mat2,axis_type(i+1,2))
          if (axis_type(i+1,2)==1.and.make_list) call add(list,'x1')
          if (axis_type(i+1,2)==2.and.make_list) call add(list,'y1')
          if (axis_type(i+1,2)==3.and.make_list) call add(list,'z1')
       else if (axis_type(i+1,1)==2) then
          call subtract_180_rotation(R_mat_new,R_mat2,axis_type(i+1,2))
          if (axis_type(i+1,2)==1.and.make_list) call add(list,'x2')
          if (axis_type(i+1,2)==2.and.make_list) call add(list,'y2')
          if (axis_type(i+1,2)==3.and.make_list) call add(list,'z2')
       else if (axis_type(i+1,1)==3) then
          call add_90_rotation(R_mat_new,R_mat2,axis_type(i+1,2))
          if (axis_type(i+1,2)==1.and.make_list) call add(list,'x3')
          if (axis_type(i+1,2)==2.and.make_list) call add(list,'y3')
          if (axis_type(i+1,2)==3.and.make_list) call add(list,'z3')
       end if

       new_theta=theta_from_r(R_mat2)
       !newz_new=axis_from_r(R_mat2)
       !if (make_list) print*,'discrete op',i+1,axis_type(i+1,2),new_theta!,axis_from_r(R_mat2)
    end do

    if (make_list) R_mat=R_mat_new

  end subroutine apply_discrete_operations_to_matrix


  subroutine subtract_180_rotation(R_mat,R_mat_new,xyz)
    implicit none

    real(gp), dimension(3,3), intent(in) :: R_mat
    real(gp), dimension(3,3), intent(out) :: R_mat_new
    integer, intent(in) :: xyz ! subtract 180 degree rotation around this axis

    real(gp), dimension(3,3) :: R_mat_180_xyz_inv

    R_mat_180_xyz_inv=0.0_gp

    if (xyz==1) then
       R_mat_180_xyz_inv(1,1)=1
       R_mat_180_xyz_inv(2,2)=-1
       R_mat_180_xyz_inv(3,3)=-1   
    else if (xyz==2) then
       R_mat_180_xyz_inv(1,1)=-1
       R_mat_180_xyz_inv(2,2)=1
       R_mat_180_xyz_inv(3,3)=-1   
    else if (xyz==3) then
       R_mat_180_xyz_inv(1,1)=-1
       R_mat_180_xyz_inv(2,2)=-1
       R_mat_180_xyz_inv(3,3)=1   
    else
       print*,'Error'
       stop
    end if

    call dgemm('N','N',3,3,3,1.0_gp,R_mat_180_xyz_inv,3,R_mat,3,0.0_gp,R_mat_new,3)

  end subroutine subtract_180_rotation


  subroutine subtract_90_rotation(R_mat,R_mat_new,xyz)
    implicit none

    real(gp), dimension(3,3), intent(in) :: R_mat
    real(gp), dimension(3,3), intent(out) :: R_mat_new
    integer, intent(in) :: xyz ! subtract 180 degree rotation around this axis

    real(gp), dimension(3,3) :: R_mat_90_xyz_inv

    R_mat_90_xyz_inv=0.0_gp

    if (xyz==1) then
       R_mat_90_xyz_inv(1,1)=1
       R_mat_90_xyz_inv(2,3)=1
       R_mat_90_xyz_inv(3,2)=-1   
    else if (xyz==2) then
       R_mat_90_xyz_inv(1,3)=-1
       R_mat_90_xyz_inv(2,2)=1
       R_mat_90_xyz_inv(3,1)=1   
    else if (xyz==3) then
       R_mat_90_xyz_inv(1,2)=1
       R_mat_90_xyz_inv(2,1)=-1
       R_mat_90_xyz_inv(3,3)=1   
    else
       print*,'Error'
       stop
    end if

    call dgemm('N','N',3,3,3,1.0_gp,R_mat_90_xyz_inv,3,R_mat,3,0.0_gp,R_mat_new,3)

  end subroutine subtract_90_rotation


  subroutine add_90_rotation(R_mat,R_mat_new,xyz)
    implicit none

    real(gp), dimension(3,3), intent(in) :: R_mat
    real(gp), dimension(3,3), intent(out) :: R_mat_new
    integer, intent(in) :: xyz ! subtract 180 degree rotation around this axis

    real(gp), dimension(3,3) :: R_mat_90_xyz_inv

    R_mat_90_xyz_inv=0.0_gp

    if (xyz==1) then
       R_mat_90_xyz_inv(1,1)=1
       R_mat_90_xyz_inv(2,3)=-1
       R_mat_90_xyz_inv(3,2)=1   
    else if (xyz==2) then
       R_mat_90_xyz_inv(1,3)=1
       R_mat_90_xyz_inv(2,2)=1
       R_mat_90_xyz_inv(3,1)=-1   
    else if (xyz==3) then
       R_mat_90_xyz_inv(1,2)=-1
       R_mat_90_xyz_inv(2,1)=1
       R_mat_90_xyz_inv(3,3)=1   
    else
       print*,'Error'
       stop
    end if

    call dgemm('N','N',3,3,3,1.0_gp,R_mat_90_xyz_inv,3,R_mat,3,0.0_gp,R_mat_new,3)

  end subroutine add_90_rotation

  pure function theta_from_r(R_mat) result(theta)
    implicit none
    real(gp), dimension(3,3), intent(in) :: R_mat

    real(gp) :: theta

    !local variables
    real(gp) :: tracem1

    tracem1=R_mat(1,1)+R_mat(2,2)+R_mat(3,3)-1.0_gp

    if (abs(tracem1) - 2.0_gp > -1.e-14_gp) then
      if (tracem1 > 0.0_gp) theta = 0.0_gp
      if (tracem1 <= 0.0_gp) theta = 180.0_gp*(4.0_gp*atan(1.d0)/180.0_gp)
    else
       theta=acos(0.5_gp*tracem1)
    end if

  end function theta_from_r

  function axis_from_r(R_mat) result(rot_axis)
    implicit none
    real(gp), dimension(3,3), intent(in) :: R_mat
    real(gp), dimension(3) :: rot_axis

    !local variables
    real(gp) :: dnrm2, norm

    rot_axis(1)=R_mat(3,2)-R_mat(2,3)    
    rot_axis(2)=R_mat(1,3)-R_mat(3,1)
    rot_axis(3)=R_mat(2,1)-R_mat(1,2)    

    !normalize it

    norm=dnrm2(3,rot_axis,1)
    !print*,'axis_from_r',norm,rot_axis
    if (norm>=1.e-5_gp) then
       call dscal(3,1.0_gp/norm,rot_axis,1)
       !print*,'axis_from_r2',norm,rot_axis
    else
       ! squares of rot_axis are diag 0.5(R+I), signs as before
       !this is not good as it gives a array of modulus bigger than one
       rot_axis(1:2)=0.0_gp
       rot_axis(3)=1.0_gp
       !rot_axis(1)=sign(dsqrt(0.5_gp*(R_mat(1,1)+1.0_gp)),rot_axis(1))
       !rot_axis(2)=sign(dsqrt(0.5_gp*(R_mat(2,2)+1.0_gp)),rot_axis(2))
       !rot_axis(3)=sign(dsqrt(0.5_gp*(R_mat(3,3)+1.0_gp)),rot_axis(3))
       !print*,'axis_from_r3',norm,rot_axis
    end if

  end function axis_from_r



  pure function frag_center(nat,rxyz) result(cen)
    implicit none
    integer, intent(in) :: nat 
    real(gp), dimension(3,nat), intent(in) :: rxyz
    real(gp), dimension(3) :: cen
    !local variables
    integer :: iat,i
    
    cen=0.0_gp
    if (nat > 0) then
       do iat=1,nat
          do i=1,3
             cen(i)=cen(i)+rxyz(i,iat)
          end do
       end do
       cen=cen/real(nat,gp)
    end if
    
  end function frag_center

  !>determinant of a 3x3 matrix
  pure function det_33(a) result(det)
    implicit none
    real(gp), dimension(3,3), intent(in) :: a
    real(gp) :: det

    det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) &
         + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3))  &
         + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
  end function det_33




  subroutine find_discrete_operations_random(rot_axis,theta,R_mat,discrete_ops)
    use module_base
    use yaml_output
    use dictionaries
    implicit none
    real(gp), dimension(3), intent(inout) :: rot_axis !< unit rotation axis (should have modulus one)
    real(gp), intent(inout) :: theta !< angle of rotation (ref becomes new)
    real(gp), dimension(3,3), intent(inout) :: R_mat
    character(len=2), dimension(:), pointer :: discrete_ops
    !local variables
    character(len=*), parameter :: subname='find_discrete_operations_random'
    integer :: j, ival, nval, i_stat
    real(gp) :: min_theta, theta_orig
    integer :: rand_size
    integer, allocatable, dimension(:) :: rand_seed
    real(kind=gp) :: rtime, rn
    character(len=10) :: sys_time 
    logical, dimension(3,3) :: stag
    type(dictionary), pointer :: list

    call dict_init(list)

    call random_seed(size=rand_size)
    allocate(rand_seed(1:rand_size))
    call date_and_time(time=sys_time)
    read(sys_time,*) rtime
    rand_seed=int(rtime*1000.0_dp)
    call random_seed(put=rand_seed)
    deallocate(rand_seed) 

    theta_orig=theta
    write(*,'(F8.4,2x)',advance='no') theta_orig/(4.0_gp*atan(1.d0)/180.0_gp)

    write(16,*) ''
    write(16,*) 'begin'
    write(16,*) theta_orig/(4.0_gp*atan(1.d0)/180.0_gp),rot_axis

    min_theta=theta_orig

    do j=1,100000
       stag=.false.
       call random_number(rn) 
       if (rn < min(j,10)*0.01_dp) then
          stag(1,1)=.true.
       else if (rn < min(j,10)*0.02_dp) then
          stag(1,2)=.true.
       else if (rn < min(j,10)*0.03_dp) then
          stag(1,3)=.true.
       else if (rn < min(j,10)*0.04_dp) then
          stag(2,1)=.true.
       else if (rn < min(j,10)*0.05_dp) then
          stag(2,2)=.true.
       else if (rn < min(j,10)*0.06_dp) then
          stag(2,3)=.true.
       else if (rn < min(j,10)*0.07_dp) then
          stag(3,1)=.true.
       else if (rn < min(j,10)*0.08_dp) then
          stag(3,2)=.true.
       else if (rn < min(j,10)*0.09_dp) then
          stag(3,3)=.true.
       end if
       call apply_discrete_operations_random(rot_axis,theta,min_theta,R_mat,stag,list)
       !write(*,'(F8.4,1x,L2,1x)',advance='no') theta/(4.0_gp*atan(1.d0)/180.0_gp),rn < min(j,10)*0.09_dp
       if (min_theta/(4.0_gp*atan(1.d0)/180.0_gp)<=45.0_gp) exit
    end do
    write(16,*) 'end'
    write(16,*) ''

    ! copy from dict to frag_trans structure
    nval=dict_len(list)
    allocate(discrete_ops(nval),stat=i_stat)
    call memocc(i_stat,discrete_ops,'discrete_ops',subname)

    do ival=0,nval-1
       discrete_ops(ival+1)=list//ival
    end do

   call dict_free(list) 

   write(*,*) 'Old theta, new theta',theta_orig/(4.0_gp*atan(1.d0)/180.0_gp), min_theta/(4.0_gp*atan(1.d0)/180.0_gp),nval

  end subroutine find_discrete_operations_random


  subroutine apply_discrete_operations_random(rot_axis,theta,mintheta,R_mat,stag,list)
    use module_base
    use dictionaries
    implicit none

    real(gp), intent(inout) :: theta
    real(gp), dimension(3), intent(inout) :: rot_axis
    real(gp), intent(inout) :: mintheta
    real(gp), dimension(3,3), intent(inout) :: R_mat
    logical, dimension(3,3), intent(in) :: stag
    type(dictionary), pointer :: list
    !local variables
    real(gp), dimension(3,3) :: R_mat_new
    integer :: i
    logical :: stag1,stag2,stag3

    do i=1,3
       call subtract_180_rotation(R_mat,R_mat_new,i)
       if (i==1) stag1=stag(1,1)
       if (i==2) stag1=stag(1,2)
       if (i==3) stag1=stag(1,3)
       if (theta_from_r(R_mat_new)<mintheta .or. stag1) then
          if (theta_from_r(R_mat_new)<mintheta) mintheta=theta_from_r(R_mat_new)
          theta=theta_from_r(R_mat_new)
          rot_axis=axis_from_r(R_mat_new)
          R_mat=R_mat_new
          if (i==1) call add(list,'x1')
          if (i==2) call add(list,'y1')
          if (i==3) call add(list,'z1')
          write(16,*) 'trans 180',i,stag1,theta/(4.0_gp*atan(1.d0)/180.0_gp),mintheta/(4.0_gp*atan(1.d0)/180.0_gp),rot_axis
          return
       end if

       call subtract_90_rotation(R_mat,R_mat_new,i)
       if (i==1) stag2=stag(2,1)
       if (i==2) stag2=stag(2,2)
       if (i==3) stag2=stag(2,3)
       if (theta_from_r(R_mat_new)<mintheta .or. stag2) then
          if (theta_from_r(R_mat_new)<mintheta) mintheta=theta_from_r(R_mat_new)
          theta=theta_from_r(R_mat_new)
          rot_axis=axis_from_r(R_mat_new)
          R_mat=R_mat_new
          if (i==1) call add(list,'x2')
          if (i==2) call add(list,'y2')
          if (i==3) call add(list,'z2')
          write(16,*) 'trans 90',i,stag2,theta/(4.0_gp*atan(1.d0)/180.0_gp),mintheta/(4.0_gp*atan(1.d0)/180.0_gp),rot_axis
          return
       end if

       call add_90_rotation(R_mat,R_mat_new,i)
       if (i==1) stag3=stag(3,1)
       if (i==2) stag3=stag(3,2)
       if (i==3) stag3=stag(3,3)
       if (theta_from_r(R_mat_new)<mintheta .or. stag3) then
          if (theta_from_r(R_mat_new)<mintheta) mintheta=theta_from_r(R_mat_new)
          theta=theta_from_r(R_mat_new)
          rot_axis=axis_from_r(R_mat_new)
          R_mat=R_mat_new
          if (i==1) call add(list,'x3')
          if (i==2) call add(list,'y3')
          if (i==3) call add(list,'z3')
          write(16,*) 'trans -90',i,stag3,theta/(4.0_gp*atan(1.d0)/180.0_gp),mintheta/(4.0_gp*atan(1.d0)/180.0_gp),rot_axis
          return
       end if
    end do

  end subroutine apply_discrete_operations_random


end module module_fragments












