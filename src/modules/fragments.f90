!> Define important structures and methods to manipulate embedding systems
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


  ! initializes all of fragment except lzd using the fragment posinp and tmb files
  subroutine init_fragment_from_file(frag,dir_name,input) ! switch this to pure if possible
    use module_types
    !use module_interfaces
    implicit none
    type(system_fragment), intent(inout) :: frag
    character(len=*), intent(in) :: dir_name
    type(input_variables), intent(in) :: input

    ! local variables
    integer :: iat

    ! nullify fragment
    frag=fragment_null()

    ! read fragment positions
    call read_atomic_file(dir_name(1:len(dir_name)-1),bigdft_mpi%iproc,frag%astruct_frg)


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

    !type, public :: minimal_orbitals_data
    !   integer :: norb          !< Total number of orbitals per k point
    !   integer :: norbp         !< Total number of orbitals for the given processors
    !   integer :: isorb         !< Total number of orbitals for the given processors
    !   integer, dimension(:), pointer :: inwhichlocreg,onwhichatom !< associate the basis centers
    !   integer, dimension(:), pointer :: isorb_par,ispot
    !   integer, dimension(:,:), pointer :: norb_par
  ! just initializing norb for now, come back and do the rest later
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
    integer :: norb, norbu, norbd, ityp, iat, ilr, istat, iall, iorb, nlr
    integer,dimension(:),allocatable :: norbsPerLocreg, norbsPerAtom
    real(kind=8),dimension(:,:),allocatable :: locregCenter
    character(len=*),parameter :: subname='init_minimal_orbitals_data'

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
    !frag%fbasis=fragment_basis_null()
    call nullify_fragment_basis(frag%fbasis)

  end function fragment_null

  pure function fragment_basis_null() result(basis)
    implicit none
    type(fragment_basis) :: basis

    basis%npsidim_orbs=0
    basis%npsidim_comp=0
    call nullify_local_zone_descriptors(basis%lzd)
    call nullify_minimal_orbitals_data(basis%forbs)
    !basis%forbs=minimal_orbitals_data_null()
    nullify(basis%psi_full)
    nullify(basis%phi)
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
    if (associated(frag%astruct_frg%rxyz)) then
       i_all=-product(shape(frag%astruct_frg%rxyz))*kind(frag%astruct_frg%rxyz)
       deallocate(frag%astruct_frg%rxyz,stat=i_stat)
       call memocc(i_stat,i_all,'frag%astruct_frg%rxyz',subname)
    end if
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

  subroutine find_frag_trans(nat,rxyz_ref,rxyz_new,rot_axis,theta)
    use module_base
    use yaml_output
    implicit none
    integer, intent(in) :: nat !< fragment size
    real(gp), dimension(3,nat), intent(in) :: rxyz_ref,rxyz_new !<coordinates measured wrt rot_center
    real(gp), dimension(3) :: rot_axis !< unit rotation axis (should have modulus one)
    real(gp) :: theta !< angle of rotation (ref becomes new)
    !local variables
    integer, parameter :: lwork=7*3
    integer :: info, i
    real(gp) :: dets,tracem1,dnrm2
    real(gp), dimension(3) :: SM_arr !< array of SVD and M array
    real(gp), dimension(lwork) :: work !< array of SVD and M array
    real(gp), dimension(3,3) :: B_mat,R_mat,U_mat,VT_mat !<matrices of Wahba's problem

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

    !to be verified that the cost function of Wahba's problem is little


    !find the angle from R matrix
    tracem1=R_mat(1,1)+R_mat(2,2)+R_mat(3,3)-1.0_gp

	 if (abs(tracem1) - 2.0_gp > -1.e-14_gp) then
	     if (tracem1 > 0.0_gp) theta = 0.0_gp
	     if (tracem1 <= 0.0_gp) theta = 180.0_gp*(4.0_gp*atan(1.d0)/180.0_gp)
    else
       theta=acos(0.5_gp*tracem1)
       !convert to degrees
       !theta=theta/(4.0_gp*atan(1.d0)/180.0_gp)
    end if
    !check the pertinence of the suggested rotation
	 if  (f_err_raise(abs(theta) > 60.d0*(4.0_gp*atan(1.d0)/180.0_gp),'Angle theta not optimal (theta= '//&
	        yaml_toa(theta)//' )')) return
    !find rot_axis
    rot_axis(1)=R_mat(3,2)-R_mat(2,3)    
    rot_axis(2)=R_mat(1,3)-R_mat(3,1)
    rot_axis(3)=R_mat(2,1)-R_mat(1,2)    
    !normalize it
    call dscal(3,1.0_gp/dnrm2(3,rot_axis,1),rot_axis,1)

    !print*,'Rmat:',theta
    !do i=1,3
    !   write(*,'(3(F12.6,2x))') R_mat(i,:)
    !end do

    !print*,'Rcalc:',theta
    !write(*,'(3(F12.6,2x))') cos(theta)+rot_axis(1)**2*(1.0_gp-cos(theta)),&
    !     rot_axis(1)*rot_axis(2)*(1.0_gp-cos(theta))-rot_axis(3)*sin(theta),&
    !     rot_axis(1)*rot_axis(3)*(1.0_gp-cos(theta))+rot_axis(2)*sin(theta)
    !write(*,'(3(F12.6,2x))') rot_axis(2)*rot_axis(1)*(1.0_gp-cos(theta))+rot_axis(3)*sin(theta),&
    !     cos(theta)+rot_axis(2)**2*(1.0_gp-cos(theta)),&
    !     rot_axis(2)*rot_axis(3)*(1.0_gp-cos(theta))-rot_axis(1)*sin(theta)
    !write(*,'(3(F12.6,2x))') rot_axis(3)*rot_axis(1)*(1.0_gp-cos(theta))-rot_axis(2)*sin(theta),&
    !     rot_axis(3)*rot_axis(2)*(1.0_gp-cos(theta))+rot_axis(1)*sin(theta),&
    !     cos(theta)+rot_axis(3)**2*(1.0_gp-cos(theta))
    


  end subroutine find_frag_trans

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

end module module_fragments

!newfrg=trans*frag


















