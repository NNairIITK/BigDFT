module psp_projectors_base
  use module_base
  use gaussians
  use locregs
  use locreg_operations
  implicit none

  private

  integer,parameter,public :: NCPLX_MAX = 2


  !> Non local pseudopotential descriptors
  type, public :: nonlocal_psp_descriptors
     integer :: mproj !< number of projectors for this descriptor
     real(gp) :: gau_cut !< cutting radius for the gaussian description of projectors.
     integer :: nlr !< total no. localization regions potentially interacting with the psp
     type(locreg_descriptors) :: plr !< localization region descriptor of a given projector (null if nlp=0)
     type(wfd_to_wfd), dimension(:), pointer :: tolr !<maskings for the locregs, dimension noverlap
     integer,dimension(:),pointer :: lut_tolr !< lookup table for tolr, dimension noverlap
     integer :: noverlap !< number of locregs which overlap with the projectors of the given atom
  end type nonlocal_psp_descriptors

  type,public :: workarrays_projectors
    real(wp),pointer,dimension(:,:,:,:) :: wprojx,wprojy,wprojz
    real(wp),pointer,dimension(:) :: wproj
    real(wp),pointer,dimension(:,:,:) :: work
  end type workarrays_projectors
  
  !> describe the information associated to the non-local part of Pseudopotentials
  type, public :: DFT_PSP_projectors
     logical :: on_the_fly             !< strategy for projector creation
     logical :: normalized             !< .true. if projectors are normalized to one.
     logical :: apply_gamma_target     !< apply the target identified by the gamma_mmp value
     integer :: nproj,nprojel,natoms   !< Number of projectors and number of elements
     real(gp) :: zerovol               !< Proportion of zero components.
     type(gaussian_basis_new) :: proj_G !< Store the projector representations in gaussians.
     real(wp), dimension(:), pointer :: proj !<storage space of the projectors in wavelet basis
     type(nonlocal_psp_descriptors), dimension(:), pointer :: pspd !<descriptor per projector, of size natom
     !> array to identify the order of the which are the atoms for which the density matrix is needed
     !! array of size natom,lmax
     integer, dimension(:,:), pointer :: iagamma
     !> density matrix for the required atoms, allocated from 1 to maxval(iagamma)
     real(wp), dimension(:,:,:,:,:), pointer :: gamma_mmp
     !>workspace for packing the wavefunctions in the case of multiple projectors
     real(wp), dimension(:), pointer :: wpack
     !> scalar product of the projectors and the wavefuntions, term by term (raw data)
     real(wp), dimension(:), pointer :: scpr
     !> full data of the scalar products
     real(wp), dimension(:), pointer :: cproj
     !> same quantity after application of the hamiltonian
     real(wp), dimension(:), pointer :: hcproj
     type(workarrays_projectors) :: wpr !< contains the workarrays for the projector creation end type DFT_PSP_projectors
  end type DFT_PSP_projectors

  public :: free_DFT_PSP_projectors
  public :: DFT_PSP_projectors_null
  public :: nonlocal_psp_descriptors_null!,free_pspd_ptr
  public :: workarrays_projectors_null, allocate_workarrays_projectors, deallocate_workarrays_projectors
contains

  pure function nonlocal_psp_descriptors_null() result(pspd)
    implicit none
    type(nonlocal_psp_descriptors) :: pspd
    call nullify_nonlocal_psp_descriptors(pspd)
  end function nonlocal_psp_descriptors_null

  pure subroutine nullify_nonlocal_psp_descriptors(pspd)
    use module_defs, only: UNINITIALIZED
    implicit none
    type(nonlocal_psp_descriptors), intent(out) :: pspd
    pspd%mproj=0
    pspd%gau_cut = UNINITIALIZED(pspd%gau_cut)
    pspd%nlr=0
    call nullify_locreg_descriptors(pspd%plr)
    nullify(pspd%tolr)
    nullify(pspd%lut_tolr)
    pspd%noverlap=0
  end subroutine nullify_nonlocal_psp_descriptors

  pure function DFT_PSP_projectors_null() result(nl)
    implicit none
    type(DFT_PSP_projectors) :: nl
    call nullify_DFT_PSP_projectors(nl)
  end function DFT_PSP_projectors_null

  pure subroutine nullify_DFT_PSP_projectors(nl)
    implicit none
    type(DFT_PSP_projectors), intent(out) :: nl
    nl%on_the_fly=.true.
    nl%normalized=.false.
    nl%apply_gamma_target=.false.
    nl%nproj=0
    nl%nprojel=0
    nl%natoms=0
    nl%zerovol=100.0_gp
    call nullify_gaussian_basis_new(nl%proj_G)
    call nullify_workarrays_projectors(nl%wpr)
    nullify(nl%iagamma)
    nullify(nl%gamma_mmp)
    nullify(nl%proj)
    nullify(nl%pspd)
    nullify(nl%wpack)
    nullify(nl%scpr)
    nullify(nl%cproj)
    nullify(nl%hcproj)
  end subroutine nullify_DFT_PSP_projectors

  !allocators


  !destructors
  subroutine deallocate_nonlocal_psp_descriptors(pspd)
    implicit none
    type(nonlocal_psp_descriptors), intent(inout) :: pspd
    !local variables
    call free_tolr_ptr(pspd%tolr)
    call f_free_ptr(pspd%lut_tolr)
    call deallocate_locreg_descriptors(pspd%plr)
  end subroutine deallocate_nonlocal_psp_descriptors

  subroutine free_pspd_ptr(pspd)
    implicit none
    type(nonlocal_psp_descriptors), dimension(:), pointer :: pspd
    !local variables
    integer :: iat

    if (.not. associated(pspd)) return
    do iat=lbound(pspd,1),ubound(pspd,1)
       call deallocate_nonlocal_psp_descriptors(pspd(iat))
    end do
    deallocate(pspd)
    nullify(pspd)

  end subroutine free_pspd_ptr

  subroutine deallocate_DFT_PSP_projectors(nl)
    implicit none
    type(DFT_PSP_projectors), intent(inout) :: nl

    call free_pspd_ptr(nl%pspd)
    nullify(nl%proj_G%rxyz)
    call gaussian_basis_free(nl%proj_G)
    call deallocate_workarrays_projectors(nl%wpr)
    call f_free_ptr(nl%iagamma)
    call f_free_ptr(nl%gamma_mmp)
    call f_free_ptr(nl%proj)
    call f_free_ptr(nl%wpack)
    call f_free_ptr(nl%scpr)
    call f_free_ptr(nl%cproj)
    call f_free_ptr(nl%hcproj)
  END SUBROUTINE deallocate_DFT_PSP_projectors


  subroutine free_DFT_PSP_projectors(nl)
    implicit none
    type(DFT_PSP_projectors), intent(inout) :: nl
    call deallocate_DFT_PSP_projectors(nl)
    call nullify_DFT_PSP_projectors(nl)
  end subroutine free_DFT_PSP_projectors


  pure function workarrays_projectors_null() result(wp)
    implicit none
    type(workarrays_projectors) :: wp
    call nullify_workarrays_projectors(wp)
  end function workarrays_projectors_null

  pure subroutine nullify_workarrays_projectors(wp)
    implicit none
    type(workarrays_projectors),intent(out) :: wp
    nullify(wp%wprojx)
    nullify(wp%wprojy)
    nullify(wp%wprojz)
    nullify(wp%wproj)
    nullify(wp%work)
  end subroutine nullify_workarrays_projectors

  subroutine allocate_workarrays_projectors(n1, n2, n3, wp)
    implicit none
    integer,intent(in) :: n1, n2, n3
    type(workarrays_projectors),intent(inout) :: wp
    integer,parameter :: nterm_max=20
    integer,parameter :: nw=65536
    wp%wprojx = f_malloc_ptr((/ 1.to.NCPLX_MAX, 0.to.n1, 1.to.2, 1.to.nterm_max /),id='wprojx')
    wp%wprojy = f_malloc_ptr((/ 1.to.NCPLX_MAX, 0.to.n2, 1.to.2, 1.to.nterm_max /),id='wprojy')
    wp%wprojz = f_malloc_ptr((/ 1.to.NCPLX_MAX, 0.to.n3, 1.to.2, 1.to.nterm_max /),id='wprojz')
    wp%wproj = f_malloc_ptr(NCPLX_MAX*(max(n1,n2,n3)+1)*2,id='wprojz')
    wp%work = f_malloc_ptr((/ 0.to.nw, 1.to.2, 1.to.2 /),id='work')
  end subroutine allocate_workarrays_projectors

  subroutine deallocate_workarrays_projectors(wp)
    implicit none
    type(workarrays_projectors),intent(inout) :: wp
    call f_free_ptr(wp%wprojx)
    call f_free_ptr(wp%wprojy)
    call f_free_ptr(wp%wprojz)
    call f_free_ptr(wp%wproj)
    call f_free_ptr(wp%work)
  end subroutine deallocate_workarrays_projectors

end module psp_projectors_base
