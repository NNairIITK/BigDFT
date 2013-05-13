!> Define important structures and methods to manipulate embedding systems
module fragments
  use module_base, only: gp,wp
  use module_types
  use dynamic_memory
  implicit none

  private

  interface operator(*)
     module procedure transform_fragment
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

  type, public :: fragment_basis
     integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
     integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
     type(local_zone_descriptors) :: Lzd
     type(minimal_orbitals_data) :: forbs
     real(wp), dimension(:), pointer :: psi
  end type fragment_basis

  !> Defines the minimal information to identify a system building block
  type, public :: system_fragment
     integer :: nat_frg !< atoms effectively belonging to the fragment
     integer :: nat_env !< environment atoms which complete fragment specifications
     integer :: ntypes  !< Number of atom types (including system and env)
     integer, dimension(:), pointer :: iatype  !< Type of the atoms (dimension nat_frag+nat_env)
     real(gp), dimension(:,:), pointer :: rxyz_frg !< position of atoms in fragment (AU), external reference frame
     real(gp), dimension(:,:), pointer :: rxyz_env !< position of atoms in environment (AU), external reference frame
     character(len=20), dimension(:), pointer :: atomnames !< Name of type of atoms
     type(fragment_basis), pointer :: frag_basis !< fragment basis, associated only if coherent with positions
  end type system_fragment

  !> Contains the rotation and translation (possibly deformation) which have to be applied to a given fragment
  type, public :: fragment_transformation
     real(gp), dimension(3) :: dr !< translation of fragment center
     real(gp), dimension(3) :: rot_center !< center of rotation in original coordinates (might be fragment center or not)
     real(gp), dimension(3) :: rot_axis !< unit rotation axis (should have modulus one)
     real(gp) :: theta !< angle of rotation 
  end type fragment_transformation

  public operator(*)

contains

  ! for reference fragments call init_fragment(frag,input%frag%frag_info(i,1),input%frag%frag_info(i,2),NTYPES,rxyz,&
  !                                            atomnames,iatype)
  ! rxyz, atomnames and iatype should come from the file fragmenti.xyz

  subroutine init_fragment(frag,nat_frg,nat_env,ntypes,rxyz,atomnames,iatype) ! switch this to pure if possible
    use module_types
    implicit none
    type(system_fragment), intent(inout) :: frag
    integer, intent(in) :: nat_frg, nat_env, ntypes
    real(gp), dimension(3,nat_frg+nat_env), intent(in) :: rxyz ! array containing positions of fragment atoms then environment atoms
    character(len=20), dimension(ntypes), intent(in) :: atomnames
    integer, dimension(nat_frg+nat_env), intent(in) :: iatype

    ! local variables
    integer :: iat

    ! nullify fragment
    frag=fragment_null()

    ! set basic integers
    frag%nat_frg=nat_frg
    frag%nat_env=nat_env
    frag%ntypes=ntypes

    ! allocate iatype, atomnames and rxyzs
    call fragment_allocate(frag)

    ! fill iatype, atomnames and rxyzs
    do iat=1,frag%nat_frg
       frag%rxyz_frg(:,iat)=rxyz(:,iat)
    end do
    
    do iat=frag%nat_frg+1,frag%nat_env
       frag%rxyz_env(:,iat)=rxyz(:,iat)
    end do

    frag%atomnames=atomnames
    frag%iatype=iatype
    

    ! allocate/initialize fragment basis...
    !currently orbitals are initialized via initAndUtils/init_orbitals_data_for_linear
    !which calls wavefunctions/orbitals_descriptors to assign parallel bits and superfluous stuff
    !and locreg_orbitals/assignToLocreg2 to give inwhichlocreg and onwhichatom - but we should really take onwhichatom from file?!

    !here we don't want to allocate psi yet, as for system fragments this won't be necessary, 
    !for reference fragments this can be done when it's filled (already nullified in fragment_basis_null)

    !for lzd do we want one for ref frags, one for whole system and one for system fragments or some kind of pointing?
    !need to think more about what should be replaced, not just what should be added in the way of initialization

    !type, public :: minimal_orbitals_data
    !   integer :: norb          !< Total number of orbitals per k point
    !   integer :: norbp         !< Total number of orbitals for the given processors
    !   integer :: isorb         !< Total number of orbitals for the given processors
    !   integer, dimension(:), pointer :: inwhichlocreg,onwhichatom !< associate the basis centers
    !   integer, dimension(:), pointer :: isorb_par,ispot
    !   integer, dimension(:,:), pointer :: norb_par

    !type, public :: fragment_basis
    !   integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
    !   integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
    !   type(local_zone_descriptors) :: Lzd
    !   type(minimal_orbitals_data) :: forbs


  end subroutine init_fragment


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

    frag%nat_frg=0
    frag%nat_env=0
    nullify(frag%rxyz_frg)
    nullify(frag%rxyz_env)

    ! nullify fragment basis
    frag%frag_basis=fragment_basis_null()

  end function fragment_null

  pure function fragment_basis_null() result(basis)
    implicit none
    type(fragment_basis) :: basis

    basis%npsidim_orbs=0
    basis%npsidim_comp=0
    call nullify_local_zone_descriptors(basis%lzd)
    call nullify_minimal_orbitals_data(basis%forbs)
    !basis%forbs=minimal_orbitals_data_null()
    nullify(basis%psi)

  end function fragment_basis_null

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
    character(len=256) :: subname  

    subname='fragment_basis_free'
    call deallocate_local_zone_descriptors(basis%lzd,subname)
    call minimal_orbitals_data_free(basis%forbs)
    if (associated(basis%psi)) call f_free_ptr(basis%psi)
    basis=fragment_basis_null()
  end subroutine fragment_basis_free

  subroutine fragment_free(frag)
    implicit none
    type(system_fragment), intent(inout) :: frag
  
    if (associated(frag%rxyz_frg)) call f_free_ptr(frag%rxyz_frg)
    if (associated(frag%rxyz_env)) call f_free_ptr(frag%rxyz_env)
    frag=fragment_null()
  end subroutine fragment_free

  subroutine fragment_allocate(frag)
    implicit none
    type(system_fragment), intent(inout) :: frag
  
    call f_routine(id='fragment_allocate')

    frag%iatype=f_malloc_ptr(frag%nat_frg+frag%nat_env,id='frag%iatype')
    !frag%atomnames=f_malloc_ptr(frag%ntypes,id='frag%atomnames') ! fmalloc objection to character here
    frag%rxyz_frg=f_malloc_ptr((/3,frag%nat_frg/),id='frag%rxyz_frg')
    frag%rxyz_env=f_malloc_ptr((/3,frag%nat_env/),id='frag%rxyz_env')

    call f_release_routine()

  end subroutine fragment_allocate


  pure function frg_center(frag)
    implicit none
    type(system_fragment), intent(in) :: frag
    real(gp), dimension(3) :: frg_center
    !local variables
    integer :: iat
    
    frg_center=0.0_gp
    do iat=1,frag%nat_frg
       frg_center=frg_center+frag%rxyz_frg(:,iat)
    end do
    frg_center=frg_center/real(frag%nat_frg,gp)

  end function frg_center

  function transform_fragment(trans,frag) result(frag_new)
    implicit none
    type(fragment_transformation), intent(in) :: trans
    type(system_fragment), intent(in) :: frag
    type(system_fragment) :: frag_new

    ! local variables
    integer :: iat

    frag_new=fragment_null()
    frag_new%nat_frg=frag%nat_frg
    frag_new%nat_env=frag%nat_env
    frag_new%ntypes=frag%ntypes

    ! allocate arrays here, leave fragment_basis nullified
    call fragment_allocate(frag_new)

    ! do fragment first, then environment
    do iat=1,frag%nat_frg
       frag_new%rxyz_frg(:,iat)=rotate_vector(trans%rot_axis,trans%theta,frag%rxyz_frg(:,iat)-trans%rot_center(:))
       frag_new%rxyz_frg(:,iat)=frag_new%rxyz_frg(:,iat)+trans%rot_center(:)+trans%dr(:)
    end do

    do iat=1,frag%nat_env
       frag_new%rxyz_env(:,iat)=rotate_vector(trans%rot_axis,trans%theta,frag%rxyz_env(:,iat)-trans%rot_center(:))
       frag_new%rxyz_env(:,iat)=frag_new%rxyz_env(:,iat)+trans%rot_center(:)+trans%dr(:)
    end do

    ! to complete should copy across iatype and atomnames but only using as a check to see how effective trans is

  end function transform_fragment



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


end module fragments

!newfrg=trans*frag





















