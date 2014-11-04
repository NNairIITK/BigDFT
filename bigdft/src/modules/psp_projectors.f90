!> @file
!! Datatypes and associated methods relativ s to the nonlocal projectors
!! @author
!!    Copyright (C) 2007-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module defining datatypes of the projectors as well as constructors and destructors
module psp_projectors
  use module_base
  use gaussians
  use locregs
  implicit none

  private

  !> Type of pseudopotential
  integer, parameter, public :: PSPCODE_UNINITIALIZED = 1
  integer, parameter, public :: PSPCODE_GTH = 2
  integer, parameter, public :: PSPCODE_HGH = 3
  integer, parameter, public :: PSPCODE_PAW = 7
  integer, parameter, public :: PSPCODE_HGH_K = 10
  integer, parameter, public :: PSPCODE_HGH_K_NLCC = 12


  !> Parameters identifying the different strategy for the application of a projector 
  !! in a localisation region
  integer, parameter, public :: PSP_APPLY_SKIP=0 !< The projector is not applied. This might happend when ilr and iat does not interact
  integer, parameter :: PSP_APPLY_MASK=1         !< Use mask arrays. The mask array has to be created before.
  integer, parameter :: PSP_APPLY_KEYS=2         !< Use keys. No mask nor packing. Equivalend to traditional application
  integer, parameter :: PSP_APPLY_MASK_PACK=3    !< Use masking and creates a pack arrays from them. 
                                                 !! Most likely this is the common usage for atoms
                                                 !! with lots of projectors and localization regions "close" to them
  integer, parameter :: PSP_APPLY_KEYS_PACK=4    !< Use keys and pack arrays. Useful especially when there is no memory to create a lot of packing arrays, 
                                                 !! for example when lots of lrs interacts with lots of atoms


  !> arrays defining how a given projector and a given wavefunction descriptor should interact
  type, public :: nlpsp_to_wfd
     integer :: strategy !< can be MASK,KEYS,MASK_PACK,KEYS_PACK,SKIP
     integer :: nmseg_c !< number of segments intersecting in the coarse region
     integer :: nmseg_f !< number of segments intersecting in the fine region
     integer, dimension(:,:), pointer :: mask !<mask array of dimesion 3,nmseg_c+nmseg_f for psp application
  end type nlpsp_to_wfd


  !> Non local pseudopotential descriptors
  type, public :: nonlocal_psp_descriptors
     integer :: mproj !< number of projectors for this descriptor
     real(gp) :: gau_cut !< cutting radius for the gaussian description of projectors.
     integer :: nlr !< total no. localization regions potentially interacting with the psp
     type(locreg_descriptors) :: plr !< localization region descriptor of a given projector (null if nlp=0)
     type(nlpsp_to_wfd), dimension(:), pointer :: tolr !<maskings for the locregs, dimension noverlap
     integer,dimension(:),pointer :: lut_tolr !< lookup table for tolr, dimension noverlap
     integer :: noverlap !< number of locregs which overlap with the projectors of the given atom
  end type nonlocal_psp_descriptors


  !> describe the information associated to the non-local part of Pseudopotentials
  type, public :: DFT_PSP_projectors 
     logical :: on_the_fly             !< strategy for projector creation
     logical :: normalized             !< .true. if projectors are normalized to one.
     integer :: nproj,nprojel,natoms   !< Number of projectors and number of elements
     real(gp) :: zerovol               !< Proportion of zero components.
     type(gaussian_basis_new) :: proj_G !< Store the projector representations is gaussians.
     real(wp), dimension(:), pointer :: proj !<storage space of the projectors in wavelet basis
     type(nonlocal_psp_descriptors), dimension(:), pointer :: pspd !<descriptor per projector, of size natom
     !>workspace for packing the wavefunctions in the case of multiple projectors
     real(wp), dimension(:), pointer :: wpack 
     !> scalar product of the projectors and the wavefuntions, term by term (raw data)
     real(wp), dimension(:), pointer :: scpr
     !> full data of the scalar products
     real(wp), dimension(:), pointer :: cproj
     !> same quantity after application of the hamiltonian
     real(wp), dimension(:), pointer :: hcproj
  end type DFT_PSP_projectors


  public :: free_DFT_PSP_projectors,update_nlpsp,hgh_psp_application,DFT_PSP_projectors_null
  public :: nonlocal_psp_descriptors_null,bounds_to_plr_limits,set_nlpsp_to_wfd,pregion_size
  public :: deallocate_nonlocal_psp_descriptors


contains


  !creators
  pure function nlpsp_to_wfd_null() result(tolr)
    implicit none
    type(nlpsp_to_wfd) :: tolr
    call nullify_nlpsp_to_wfd(tolr)
  end function nlpsp_to_wfd_null
  pure subroutine nullify_nlpsp_to_wfd(tolr)
    implicit none
    type(nlpsp_to_wfd), intent(out) :: tolr
    tolr%strategy=PSP_APPLY_SKIP
    tolr%nmseg_c=0
    tolr%nmseg_f=0
    nullify(tolr%mask)
  end subroutine nullify_nlpsp_to_wfd

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
    nl%nproj=0
    nl%nprojel=0
    nl%natoms=0
    nl%zerovol=100.0_gp
    call nullify_gaussian_basis_new(nl%proj_G)! = gaussian_basis_null()
    nullify(nl%proj)
    nullify(nl%pspd)
    nullify(nl%wpack)
    nullify(nl%scpr)
    nullify(nl%cproj)
    nullify(nl%hcproj)
  end subroutine nullify_DFT_PSP_projectors

  !allocators

  !destructors
  subroutine deallocate_nlpsp_to_wfd(tolr)
    implicit none
    type(nlpsp_to_wfd), intent(inout) :: tolr
    call f_free_ptr(tolr%mask)
  end subroutine deallocate_nlpsp_to_wfd


  subroutine deallocate_nonlocal_psp_descriptors(pspd)
    implicit none
    type(nonlocal_psp_descriptors), intent(inout) :: pspd
    !local variables
    integer :: ilr
    if (associated(pspd%tolr)) then
       do ilr=1,size(pspd%tolr)
          call deallocate_nlpsp_to_wfd(pspd%tolr(ilr))
          call nullify_nlpsp_to_wfd(pspd%tolr(ilr))
       end do
       deallocate(pspd%tolr)
       nullify(pspd%tolr)
    end if
    call deallocate_locreg_descriptors(pspd%plr)
    call f_free_ptr(pspd%lut_tolr)
  end subroutine deallocate_nonlocal_psp_descriptors


  subroutine deallocate_DFT_PSP_projectors(nl)
    implicit none
    type(DFT_PSP_projectors), intent(inout) :: nl
    !local variables
    integer :: iat

    if (associated(nl%pspd)) then
       do iat=1,nl%natoms
          call deallocate_nonlocal_psp_descriptors(nl%pspd(iat))
       end do
       deallocate(nl%pspd)
       nullify(nl%pspd)
    end if
    nullify(nl%proj_G%rxyz)
    call gaussian_basis_free(nl%proj_G)
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


  !then routines which are typical of the projector application or creation follow

  !> converts the bound of the lr descriptors in local bounds of the plr locregs
  pure subroutine bounds_to_plr_limits(thatway,icoarse,plr,nl1,nl2,nl3,nu1,nu2,nu3)
    implicit none
    logical, intent(in) :: thatway !< if .true., the plr descriptors has to be filled
    !! if .false., the nl bounds are filled from the plr
    integer, intent(in) :: icoarse !<controls whether to assign coarse or fine
    !!limits. Can be 1 or 2. 
    !!The 2 case cannot be doe before 
    !!the case with 1 has been filled
    type(locreg_descriptors), intent(inout) :: plr !<projectors locreg
    integer, intent(inout) :: nl1,nl2,nl3,nu1,nu2,nu3 !<lower and upper bounds of locregs

    if (thatway) then
       if (icoarse==1) then !coarse limits (to be done first)
          plr%ns1=nl1     
          plr%ns2=nl2       
          plr%ns3=nl3       

          plr%d%n1=nu1-nl1
          plr%d%n2=nu2-nl2
          plr%d%n3=nu3-nl3
       else if (icoarse == 2) then
          plr%d%nfl1=nl1-plr%ns1
          plr%d%nfl2=nl2-plr%ns2       
          plr%d%nfl3=nl3-plr%ns3

          plr%d%nfu1=nu1-plr%ns1
          plr%d%nfu2=nu2-plr%ns2
          plr%d%nfu3=nu3-plr%ns3
!!$         else
!!$            stop 'WRONG icoarse' 
       end if
    else
       if (icoarse==1) then !coarse limits
          nl1=plr%ns1
          nl2=plr%ns2
          nl3=plr%ns3

          nu1=plr%d%n1+plr%ns1
          nu2=plr%d%n2+plr%ns2
          nu3=plr%d%n3+plr%ns3
       else if (icoarse == 2) then
          nl1=plr%d%nfl1+plr%ns1
          nl2=plr%d%nfl2+plr%ns2
          nl3=plr%d%nfl3+plr%ns3

          nu1=plr%d%nfu1+plr%ns1
          nu2=plr%d%nfu2+plr%ns2
          nu3=plr%d%nfu3+plr%ns3
!!$         else
!!$            stop 'WRONG icoarse, false case' 
       end if
    end if
  end subroutine bounds_to_plr_limits


  !> Finds the size of the smallest subbox that contains a localization region made 
  !! out of atom centered spheres
  subroutine pregion_size(geocode,rxyz,radius,rmult,hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
    !use module_base, only: gp
    implicit none
    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
    integer, intent(in) :: n1,n2,n3
    real(gp), intent(in) :: hx,hy,hz,rmult,radius
    real(gp), dimension(3), intent(in) :: rxyz
    integer, intent(out) :: nl1,nu1,nl2,nu2,nl3,nu3
    !Local variables
    double precision, parameter :: eps_mach=1.d-12
    !n(c) real(kind=8) :: onem
    real(gp) :: cxmax,cymax,czmax,cxmin,cymin,czmin,rad

    rad=radius*rmult
    cxmax=rxyz(1)+rad ; cxmin=rxyz(1)-rad
    cymax=rxyz(2)+rad ; cymin=rxyz(2)-rad
    czmax=rxyz(3)+rad ; czmin=rxyz(3)-rad
    !n(c) onem=1.d0-eps_mach

    nl1=ceiling(real(cxmin/hx,kind=8) - eps_mach)   
    nl2=ceiling(real(cymin/hy,kind=8) - eps_mach)   
    nl3=ceiling(real(czmin/hz,kind=8) - eps_mach)   
    nu1=floor(real(cxmax/hx,kind=8) + eps_mach)  
    nu2=floor(real(cymax/hy,kind=8) + eps_mach)  
    nu3=floor(real(czmax/hz,kind=8) + eps_mach)  

    !for non-free BC the projectors are not allowed to be also outside the box
    if (geocode == 'F') then
       if (nl1 < 0)   stop 'nl1: projector region outside cell'
       if (nl2 < 0)   stop 'nl2: projector region outside cell'
       if (nl3 < 0)   stop 'nl3: projector region outside cell'
       if (nu1 > n1)   stop 'nu1: projector region outside cell'
       if (nu2 > n2)   stop 'nu2: projector region outside cell'
       if (nu3 > n3)   stop 'nu3: projector region outside cell'
    else if (geocode == 'S') then
       !correct the extremes if they run outside the box
       if (nl1 < 0 .or. nu1 > n1) then
          nl1=0
          nu1=n1
       end if
       if (nl2 < 0)   stop 'nl2: projector region outside cell'
       if (nu2 > n2)   stop 'nu2: projector region outside cell'
       if (nl3 < 0 .or. nu3 > n3) then
          nl3=0
          nu3=n3
       end if
    else if (geocode == 'P') then
       !correct the extremes if they run outside the box
       if (nl1 < 0 .or. nu1 > n1) then
          nl1=0
          nu1=n1
       end if
       if (nl2 < 0 .or. nu2 > n2) then
          nl2=0
          nu2=n2
       end if
       if (nl3 < 0 .or. nu3 > n3) then
          nl3=0
          nu3=n3
       end if
    end if

  END SUBROUTINE pregion_size

  !> routine to update the PSP descriptors as soon as the localization regions
  ! are modified
  subroutine update_nlpsp(nl,nlr,lrs,Glr,lr_mask)
    implicit none
    integer, intent(in) :: nlr
    type(locreg_descriptors), intent(in) :: Glr
    !>logical array of the localization regions active on site
    !it is true for all the elements corresponding to localisation 
    !! regions whose descriptors are calculated
    logical, dimension(nlr), intent(in) :: lr_mask
    type(locreg_descriptors), dimension(nlr), intent(in) :: lrs
    type(DFT_PSP_projectors), intent(inout) :: nl
    !local variables
    integer :: nbseg_dim,nkeyg_dim,iat,ilr
    integer, dimension(:), allocatable :: nbsegs_cf,keyg_lin


    call f_routine(id='update_nlpsp')

    !find allocating dimensions for work arrays
    nbseg_dim=0
    do iat=1,nl%natoms
       nbseg_dim=max(nbseg_dim,&
            nl%pspd(iat)%plr%wfd%nseg_c+nl%pspd(iat)%plr%wfd%nseg_f)
    end do
    nkeyg_dim=0
    do ilr=1,nlr
       nkeyg_dim=max(nkeyg_dim,lrs(ilr)%wfd%nseg_c+lrs(ilr)%wfd%nseg_f)
    end do

    !allocate the work arrays for building tolr array of structures
    nbsegs_cf=f_malloc(nbseg_dim,id='nbsegs_cf')
    keyg_lin=f_malloc(nkeyg_dim,id='keyg_lin')
    !reconstruct the projectors for any of the atoms
    do iat=1,nl%natoms
       !free the pre-existing array of structures
       if (associated(nl%pspd(iat)%tolr)) then
          do ilr=1,size(nl%pspd(iat)%tolr)
             call deallocate_nlpsp_to_wfd(nl%pspd(iat)%tolr(ilr))
             call nullify_nlpsp_to_wfd(nl%pspd(iat)%tolr(ilr))
          end do
          deallocate(nl%pspd(iat)%tolr)
          nullify(nl%pspd(iat)%tolr)
          call f_free_ptr(nl%pspd(iat)%lut_tolr)
       end if
       if (nl%pspd(iat)%mproj > 0) then 
          !then fill it again, if the locreg is demanded
          nl%pspd(iat)%nlr=nlr
          call set_nlpsp_to_wfd(Glr,nl%pspd(iat)%plr,&
               keyg_lin,nbsegs_cf,nl%pspd(iat)%noverlap,nl%pspd(iat)%lut_tolr,nl%pspd(iat)%tolr,lrs,lr_mask)
       end if
    end do

    call f_free(keyg_lin)
    call f_free(nbsegs_cf)

    call f_release_routine()

  end subroutine update_nlpsp

  !> initialize the information for matching the localisation region
  !! of each projector to all the localisation regions of the system
  subroutine set_nlpsp_to_wfd(Glr,plr,keyag_lin_cf,nbsegs_cf,noverlap,lut_tolr,tolr,lrs,lr_mask)
    implicit none
    type(locreg_descriptors), intent(in) :: Glr !<global simulation domain
    type(locreg_descriptors), intent(in) :: plr !<locreg of the projector
    !>work array: needed to build the mask array
    integer, dimension(plr%wfd%nseg_c+plr%wfd%nseg_f), intent(inout) :: nbsegs_cf
    !>work array: needed to put the unstrided keygloc of all locregs
    !! the dimension has to be maxval(lrs(:)%nseg_c+lrs(:)%nseg_f)
    integer, dimension(*), intent(inout) :: keyag_lin_cf
    !>structures which have to be filled to prepare projector applications
    integer,dimension(:),pointer,intent(out) :: lut_tolr !< lookup table
    integer,intent(out) :: noverlap !< dimension of the arrays lut_tolr and tolr
    type(nlpsp_to_wfd), dimension(:), pointer,intent(out) :: tolr 
    !> mask array which is associated to the localization regions of interest in this processor
    logical, dimension(:), optional, intent(in) :: lr_mask
    !> descriptors of all the localization regions of the simulation domain
    !! susceptible to interact with the projector
    type(locreg_descriptors), dimension(:), optional, intent(in) :: lrs
    !local variables
    logical :: overlap
    integer :: ilr,nlr, iilr, ioverlap

    call f_routine(id='set_nlpsp_to_wfd')

    ! Determine the size of tolr and initialize the corresponding lookup table
    nlr=1
    overlap=.true.
    if (present(lrs)) then
        nlr=size(lrs)
    end if
    if (nlr <=0) return

    ! Count how many overlaps exist
    noverlap=0
    do ilr=1,nlr
       !control if the projector overlaps with this locreg
       if (present(lrs)) then
          overlap=.true.
          if (present(lr_mask)) overlap=lr_mask(ilr)
          if (overlap) call check_overlap(lrs(ilr),plr,Glr,overlap)
          !if there is overlap, activate the strategy for the application
          if (overlap) then
              noverlap=noverlap+1
          end if
       else
          noverlap=noverlap+1
       end if
    end do
    lut_tolr = f_malloc_ptr(noverlap,id='lut_tolr')
    lut_tolr = PSP_APPLY_SKIP

    ! Now assign the values
    ioverlap=0
    do ilr=1,nlr
       !control if the projector overlaps with this locreg
       if (present(lrs)) then
          overlap=.true.
          if (present(lr_mask)) overlap=lr_mask(ilr)
          if (overlap) call check_overlap(lrs(ilr),plr,Glr,overlap)
          !if there is overlap, activate the strategy for the application
          if (overlap) then
              ioverlap=ioverlap+1
              lut_tolr(ioverlap)=ilr
          end if
       else
          ioverlap=ioverlap+1
          lut_tolr(ioverlap)=ilr
       end if
    end do
    if (ioverlap/=noverlap) stop 'ioverlap/=noverlap'



    if (present(lrs) .and. present(lr_mask)) then
       if (f_err_raise(nlr /= size(lr_mask),'The sizes of lr_mask and lrs should coincide',&
            err_name='BIGDFT_RUNTIME_ERROR')) return
    end if
    !allocate the pointer with the good size
    allocate(tolr(noverlap))
    !then for any of the localization regions check the strategy
    !for applying the projectors
    !ioverlap=0
    do ilr=1,noverlap
       iilr=lut_tolr(ilr)
       !if (iilr==PSP_APPLY_SKIP) cycle
       !ioverlap=ioverlap+1
       !this will set to PSP_APPLY_SKIP the projector application
       call nullify_nlpsp_to_wfd(tolr(ilr))
       !now control if the projector overlaps with this locreg
       if (present(lrs)) then
          overlap=.true.
          if (present(lr_mask)) overlap=lr_mask(iilr)
          if (overlap) call check_overlap(lrs(iilr),plr,Glr,overlap)
          !if there is overlap, activate the strategy for the application
          if (overlap) then
             call init_tolr(tolr(ilr),lrs(iilr)%wfd,plr%wfd,keyag_lin_cf,nbsegs_cf)
          end if
       else
          call init_tolr(tolr(ilr),Glr%wfd,plr%wfd,keyag_lin_cf,nbsegs_cf)
       end if
       !then the best strategy can be decided according to total number of 
       !common points
       !complete stategy, the packing array is created after first projector
       if (overlap) tolr(ilr)%strategy=PSP_APPLY_MASK_PACK
       !masking is used but packing is not created, 
       !useful when only one projector has to be applied
       !tolr(ilr)%strategy=PSP_APPLY_MASK
       !old scheme, even though mask arrays is created it is never used.
       !most likely this scheme is useful for debugging purposes
       !tolr(ilr)%strategy=PSP_APPLY_KEYS
    end do

    !!if (ioverlap/=noverlap) stop 'ioverlap/=noverlap'

    call f_release_routine()

  end subroutine set_nlpsp_to_wfd

  !> initialize the nlpsp_to_wfd descriptor starting from 
  !! the descriptors of the localization regions
  subroutine init_tolr(tolr,wfd_lr,wfd_p,keyag_lin_cf,nbsegs_cf)
    implicit none
    !>descriptors of the localization region to mask
    type(wavefunctions_descriptors), intent(in) :: wfd_lr
    !>descriptors of the projectors
    type(wavefunctions_descriptors), intent(in) :: wfd_p
    !> array of the unstrided keyglob starting points of wfd_w
    integer, dimension(wfd_lr%nseg_c+wfd_lr%nseg_f), intent(inout) :: keyag_lin_cf
    !> number of common segments of the wfd_w for each of the segment of wfd_p.
    integer, dimension(wfd_p%nseg_c+wfd_p%nseg_f), intent(inout) :: nbsegs_cf
    !> structure for apply the projector to the corresponding locreg
    type(nlpsp_to_wfd), intent(inout) :: tolr

    call f_routine(id='init_tolr')
    
    !calculate the size of the mask array
    call vcopy(wfd_lr%nseg_c+wfd_lr%nseg_f,&
         wfd_lr%keyglob(1,1),2,keyag_lin_cf(1),1)
    call to_zero(wfd_p%nseg_c+wfd_p%nseg_f,nbsegs_cf(1))
    call mask_sizes(wfd_lr,wfd_p,keyag_lin_cf,nbsegs_cf,&
         tolr%nmseg_c,tolr%nmseg_f)
    !then allocate and fill it
    tolr%mask=&
         f_malloc0_ptr((/3,tolr%nmseg_c+tolr%nmseg_f/),&
         id='mask')
    !and filled
    call init_mask(wfd_lr,wfd_p,keyag_lin_cf,nbsegs_cf,&
         tolr%nmseg_c,tolr%nmseg_f,tolr%mask)

    call f_release_routine()

  end subroutine init_tolr


  !>find the size of the mask array for a given couple plr - llr
  subroutine mask_sizes(wfd_w,wfd_p,keyag_lin_cf,nbsegs_cf,nmseg_c,nmseg_f)

    implicit none
    type(wavefunctions_descriptors), intent(in) :: wfd_w,wfd_p
    !> array of the unstrided keyglob starting points of wfd_w, pre-filled
    integer, dimension(wfd_w%nseg_c+wfd_w%nseg_f), intent(in) :: keyag_lin_cf
    !> number of common segments of the wfd_w for each of the segment of wfd_p.
    !! should be initialized to zero at input
    integer, dimension(wfd_p%nseg_c+wfd_p%nseg_f), intent(inout) :: nbsegs_cf
    integer, intent(out) :: nmseg_c,nmseg_f

    call count_wblas_segs(wfd_w%nseg_c,wfd_p%nseg_c,keyag_lin_cf(1),&
         wfd_w%keyglob(1,1),wfd_p%keyglob(1,1),nbsegs_cf(1))
    !  print *,'no of points',sum(nbsegs_cf),wfd_w%nseg_c,wfd_p%nseg_c
    call integrate_nseg(wfd_p%nseg_c,nbsegs_cf(1),nmseg_c)
    !  print *,'no of points',nmseg_c

    if (wfd_w%nseg_f >0 .and. wfd_p%nseg_f > 0 ) then
       call count_wblas_segs(wfd_w%nseg_f,wfd_p%nseg_f,keyag_lin_cf(wfd_w%nseg_c+1),&
            wfd_w%keyglob(1,wfd_w%nseg_c+1),wfd_p%keyglob(1,wfd_p%nseg_c+1),&
            nbsegs_cf(wfd_p%nseg_c+1))
       call integrate_nseg(wfd_p%nseg_f,nbsegs_cf(wfd_p%nseg_c+1),nmseg_f)
    else
       nmseg_f=0
    end if

  contains

    !> count the total number of segments and define the integral array of displacements
    pure subroutine integrate_nseg(mseg,msegs,nseg_tot)
      implicit none
      integer, intent(in) :: mseg
      integer, dimension(mseg), intent(inout) :: msegs
      integer, intent(out) :: nseg_tot
      !local variables
      integer :: iseg,jseg

      nseg_tot=0
      do iseg=1,mseg
         jseg=msegs(iseg)
         msegs(iseg)=nseg_tot
         nseg_tot=nseg_tot+jseg
      end do
    end subroutine integrate_nseg

  end subroutine mask_sizes

  !>fill the mask array which has been previoulsly allocated and cleaned
  subroutine init_mask(wfd_w,wfd_p,keyag_lin_cf,nbsegs_cf,nmseg_c,nmseg_f,mask)
    implicit none
    integer, intent(in) :: nmseg_c,nmseg_f
    type(wavefunctions_descriptors), intent(in) :: wfd_w,wfd_p
    !> array of the unstrided keyglob starting points of wfd_w, pre-filled
    integer, dimension(wfd_w%nseg_c+wfd_w%nseg_f), intent(in) :: keyag_lin_cf
    !> number of common segments of the wfd_w for each of the segment of wfd_p.
    !! should be created by mask_sizes routine
    integer, dimension(wfd_p%nseg_c+wfd_p%nseg_f), intent(in) :: nbsegs_cf
    !>masking array. On output, it indicates for any of the segments 
    !which are common between the wavefunction and the projector
    !the starting positions in the packed arrays of projectors and wavefunction
    !respectively
    integer, dimension(3,nmseg_c+nmseg_f), intent(inout) :: mask

    call fill_wblas_segs(wfd_w%nseg_c,wfd_p%nseg_c,nmseg_c,&
         nbsegs_cf(1),keyag_lin_cf(1),wfd_w%keyglob(1,1),wfd_p%keyglob(1,1),&
         wfd_w%keyvglob(1),wfd_p%keyvglob(1),mask(1,1))
    if (nmseg_f > 0) then
       call fill_wblas_segs(wfd_w%nseg_f,wfd_p%nseg_f,nmseg_f,&
            nbsegs_cf(wfd_p%nseg_c+1),keyag_lin_cf(wfd_w%nseg_c+1),&
            wfd_w%keyglob(1,wfd_w%nseg_c+1),wfd_p%keyglob(1,wfd_p%nseg_c+1),&
            wfd_w%keyvglob(wfd_w%nseg_c+1),wfd_p%keyvglob(wfd_p%nseg_c+1),&
            mask(1,nmseg_c+1))
    end if

  end subroutine init_mask

  !> applies a projector of HGH type, written in a localization region
  !! onto a wavefunction written in the same formalism
  !! uses the desctiptors for the application which have been defined previously
  ! replace the routine nl_HGH_application as it does not need allocating arrays anymore
  subroutine hgh_psp_application(hij,ncplx_p,n_p,wfd_p,proj,&
       ncplx_w,n_w,wfd_w,tolr,psi_pack,scpr,pdpsi,hpdpsi,psi,hpsi,eproj)
    implicit none
    integer, intent(in) :: ncplx_p !< number of complex components of the projector
    integer, intent(in) :: n_p !< number of elements of the projector
    integer, intent(in) :: ncplx_w !< number of complex components of the wavefunction
    integer, intent(in) :: n_w !< number of complex components of the wavefunction
    type(wavefunctions_descriptors), intent(in) :: wfd_p !< descriptors of projectors
    type(wavefunctions_descriptors), intent(in) :: wfd_w !< descriptors of wavefunction
    !> interaction between the wavefuntion and the psp projector
    type(nlpsp_to_wfd), intent(in) :: tolr
    !> matrix of nonlocal HGH psp
    real(gp), dimension(3,3,4), intent(in) :: hij
    !> components of the projectors, real and imaginary parts
    real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,ncplx_p,n_p), intent(in) :: proj 
    !> components of wavefunctions, real and imaginary parts
    real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,ncplx_w,n_w), intent(in) :: psi 
    !> components of wavefunctions after application, real and imaginary parts
    real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,ncplx_w,n_w), intent(inout) :: hpsi 
    !> workspaces for the packing array
    real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_w*ncplx_w), intent(inout) :: psi_pack
    !> array of the scalar product between the projectors and the wavefunctions
    real(wp), dimension(ncplx_w,n_w,ncplx_p,n_p), intent(inout) :: scpr
    !> array of the coefficients of the hgh projectors
    real(wp), dimension(max(ncplx_w,ncplx_p),n_w,n_p), intent(inout) :: pdpsi
    !> array of the coefficients of the hgh projectors multiplied by HGH matrix
    real(wp), dimension(max(ncplx_w,ncplx_p),n_w,n_p), intent(inout) :: hpdpsi
    real(gp), intent(out) :: eproj !< energy of the projectors

    eproj=0.0_gp
    !put to zero the array
    call to_zero((wfd_p%nvctr_c+7*wfd_p%nvctr_f)*n_w*ncplx_w,psi_pack(1,1))

    !here also the strategy can be considered
    call proj_dot_psi(n_p*ncplx_p,wfd_p,proj,n_w*ncplx_w,wfd_w,psi,&
         tolr%nmseg_c,tolr%nmseg_f,tolr%mask,psi_pack,scpr)

    !first create the coefficients for the application of the matrix
    !pdpsi = < p_i | psi >
    call full_coefficients('C',ncplx_p,n_p,'N',ncplx_w,n_w,scpr,'N',pdpsi)

    call apply_hij_coeff(hij,max(ncplx_w,ncplx_p)*n_w,n_p,pdpsi,hpdpsi)

    !then create the coefficients for the evaluation of the projector energy
    !pdpsi= < psi | p_i> = conj(< p_i | psi >)
    call full_coefficients('N',ncplx_p,n_p,'C',ncplx_w,n_w,scpr,'C',pdpsi)

    !the energy can be calculated here
    eproj=dot(max(ncplx_p,ncplx_w)*n_w*n_p,pdpsi(1,1,1),1,hpdpsi(1,1,1),1)

    !then the coefficients have to be transformed for the projectors
    call reverse_coefficients(ncplx_p,n_p,ncplx_w,n_w,hpdpsi,scpr)

    call scpr_proj_p_hpsi(n_p*ncplx_p,wfd_p,proj,n_w*ncplx_w,wfd_w,&
         tolr%nmseg_c,tolr%nmseg_f,tolr%mask,psi_pack,scpr,hpsi)

  end subroutine hgh_psp_application

  !> routine for applying the coefficients needed HGH-type PSP to the scalar product
  !! among wavefunctions and projectors. The coefficients are real therefore 
  !! there is no need to separate scpr in its real and imaginary part before
  pure subroutine apply_hij_coeff(hij,n_w,n_p,scpr,hscpr)
    use module_base, only: gp
    implicit none
    integer, intent(in) :: n_p,n_w
    real(gp), dimension(3,3,4), intent(in) :: hij
    real(gp), dimension(n_w,n_p), intent(in) :: scpr
    real(gp), dimension(n_w,n_p), intent(out) :: hscpr
    !local variables
    integer :: i,j,l,m,iproj,iw
    real(gp), dimension(7,3,4) :: cproj,dproj 
    logical, dimension(3,4) :: cont

!!$    !fill the hij matrix
!!$    call hgh_hij_matrix(npspcode,psppar,hij)

    !define the logical array to identify the point from which the block is finished
    do l=1,4
       do i=1,3
          cont(i,l)=(hij(i,i,l) /= 0.0_gp)
       end do
    end do
   

    reversed_loop: do iw=1,n_w
       dproj=0.0_gp

       iproj=1
       !fill the complete coefficients
       do l=1,4 !diagonal in l
          do i=1,3
             if (cont(i,l)) then !psppar(l,i) /= 0.0_gp) then
                do m=1,2*l-1
                   cproj(m,i,l)=scpr(iw,iproj)
                   iproj=iproj+1
                end do
             else
                do m=1,2*l-1
                   cproj(m,i,l)=0.0_gp
                end do
             end if
          end do
       end do

       !applies the hij matrix
       do l=1,4 !diagonal in l
          do i=1,3
             do j=1,3
                do m=1,2*l-1 !diagonal in m
                   dproj(m,i,l)=dproj(m,i,l)+&
                        hij(i,j,l)*cproj(m,j,l)
                end do
             end do
          end do
       end do

       !copy back the coefficient
       iproj=1
       !fill the complete coefficients
       do l=1,4 !diagonal in l
          do i=1,3
             if (cont(i,l)) then !psppar(l,i) /= 0.0_gp) then
                do m=1,2*l-1
                   hscpr(iw,iproj)=dproj(m,i,l)
                   iproj=iproj+1
                end do
             end if
          end do
       end do
    end do reversed_loop

  end subroutine apply_hij_coeff

  pure subroutine reverse_coefficients(ncplx_p,n_p,ncplx_w,n_w,pdpsi,scpr)
    implicit none
    integer, intent(in) :: ncplx_p,ncplx_w,n_p,n_w
    real(wp), dimension(max(ncplx_w,ncplx_p),n_w,n_p), intent(in) :: pdpsi
    real(wp), dimension(ncplx_w,n_w,ncplx_p,n_p), intent(out) :: scpr
    !local variables
    logical :: cplx_p,cplx_w,cplx_pw
    integer :: iw,ip,icplx

    cplx_p=ncplx_p==2
    cplx_w=ncplx_w==2
    cplx_pw=cplx_p .and. cplx_w

    if (cplx_pw) then
       do ip=1,n_p
          do iw=1,n_w
             scpr(1,iw,1,ip)=pdpsi(1,iw,ip)
             scpr(2,iw,1,ip)=pdpsi(2,iw,ip)
             scpr(1,iw,2,ip)=-pdpsi(2,iw,ip)
             scpr(2,iw,2,ip)=pdpsi(1,iw,ip)
          end do
       end do
       !copy the values, only one of the two might be 2
    else if (cplx_p) then
       do ip=1,n_p
          do icplx=1,ncplx_p
             do iw=1,n_w
                scpr(1,iw,icplx,ip)=pdpsi(icplx,iw,ip)
             end do
          end do
       end do
    else if (cplx_w) then
       do ip=1,n_p
          do iw=1,n_w
             do icplx=1,ncplx_w
                scpr(icplx,iw,1,ip)=pdpsi(icplx,iw,ip)
             end do
          end do
       end do
    else !real case
       do ip=1,n_p
          do iw=1,n_w
             scpr(1,iw,1,ip)=pdpsi(1,iw,ip)
          end do
       end do

    end if
  end subroutine reverse_coefficients

  !> Identify the coefficients
  pure subroutine full_coefficients(trans_p,ncplx_p,n_p,trans_w,ncplx_w,n_w,scpr,trans,pdpsi)
    implicit none
    integer, intent(in) :: ncplx_p,ncplx_w,n_p,n_w
    character(len=1), intent(in) :: trans_p,trans_w,trans
    real(wp), dimension(ncplx_w,n_w,ncplx_p,n_p), intent(in) :: scpr
    real(wp), dimension(max(ncplx_w,ncplx_p),n_w,n_p), intent(out) :: pdpsi
    !local variables
    logical :: cplx_p,cplx_w,cplx_pw
    integer :: iw,ip,ieps_p,ieps_w,ieps
    real(wp) :: prfr,prfi,pifr,pifi

    cplx_p=ncplx_p==2
    cplx_w=ncplx_w==2
    cplx_pw=cplx_p .and. cplx_w

    ieps_p=1
    if (trans_p=='C' .and. cplx_p) ieps_p=-1
    ieps_w=1
    if (trans_w=='C' .and. cplx_w) ieps_w=-1
    ieps=1
    if (trans=='C' .and. (cplx_p .or. cplx_w)) ieps=-1


    !the coefficients have to be transformed to the complex version
    if ((.not. cplx_p) .and. (.not.  cplx_w)) then
       !real case, simply copy the values
       do ip=1,n_p
          do iw=1,n_w
             pdpsi(1,iw,ip)=scpr(1,iw,1,ip)
          end do
       end do
    else
       !complex case, build real and imaginary part when applicable
       prfi=0.0_wp
       pifr=0.0_wp
       pifi=0.0_wp
       do ip=1,n_p
          do iw=1,n_w
             prfr=scpr(1,iw,1,ip)
             if (cplx_p) pifr=scpr(1,iw,2,ip)
             if (cplx_w) prfi=scpr(2,iw,1,ip)
             if (cplx_pw) pifi=scpr(2,iw,2,ip)   
             !real part
             pdpsi(1,iw,ip)=prfr-ieps_p*ieps_w*pifi
             !imaginary part
             pdpsi(2,iw,ip)=ieps*ieps_w*prfi+ieps*ieps_p*pifr
          end do
       end do
    end if

  end subroutine full_coefficients

end module psp_projectors


!> External routine as the psppar parameters are often passed by address
subroutine hgh_hij_matrix(npspcode,psppar,hij)
  use module_defs, only: gp
  use psp_projectors, only: PSPCODE_GTH, PSPCODE_HGH, PSPCODE_HGH_K, PSPCODE_HGH_K_NLCC, PSPCODE_PAW
  implicit none
  !Arguments
  integer, intent(in) :: npspcode
  real(gp), dimension(0:4,0:6), intent(in) :: psppar
  real(gp), dimension(3,3,4), intent(out) :: hij
  !Local variables
  integer :: l,i,j
  real(gp), dimension(2,2,3) :: offdiagarr

  !enter the coefficients for the off-diagonal terms (HGH case, npspcode=PSPCODE_HGH)
  offdiagarr(1,1,1)=-0.5_gp*sqrt(3._gp/5._gp)
  offdiagarr(2,1,1)=-0.5_gp*sqrt(100._gp/63._gp)
  offdiagarr(1,2,1)=0.5_gp*sqrt(5._gp/21._gp)
  offdiagarr(2,2,1)=0.0_gp !never used
  offdiagarr(1,1,2)=-0.5_gp*sqrt(5._gp/7._gp)  
  offdiagarr(2,1,2)=-7._gp/3._gp*sqrt(1._gp/11._gp)
  offdiagarr(1,2,2)=1._gp/6._gp*sqrt(35._gp/11._gp)
  offdiagarr(2,2,2)=0.0_gp !never used
  offdiagarr(1,1,3)=-0.5_gp*sqrt(7._gp/9._gp)
  offdiagarr(2,1,3)=-9._gp*sqrt(1._gp/143._gp)
  offdiagarr(1,2,3)=0.5_gp*sqrt(63._gp/143._gp)
  offdiagarr(2,2,3)=0.0_gp !never used

  !  call to_zero(3*3*4,hij(1,1,1))
  hij=0.0_gp

  do l=1,4
     !term for all npspcodes
     loop_diag: do i=1,3
        hij(i,i,l)=psppar(l,i) !diagonal term
        if ((npspcode == PSPCODE_HGH .and. l/=4 .and. i/=3) .or. &
             ((npspcode == PSPCODE_HGH_K .or. npspcode == PSPCODE_HGH_K_NLCC) .and. i/=3)) then !HGH(-K) case, offdiagonal terms
           loop_offdiag: do j=i+1,3
              if (psppar(l,j) == 0.0_gp) exit loop_offdiag
              !offdiagonal HGH term
              if (npspcode == PSPCODE_HGH) then !traditional HGH convention
                 hij(i,j,l)=offdiagarr(i,j-i,l)*psppar(l,j)
              else !HGH-K convention
                 hij(i,j,l)=psppar(l,i+j+1)
              end if
              hij(j,i,l)=hij(i,j,l) !symmetrization
           end do loop_offdiag
        end if
     end do loop_diag
  end do

end subroutine hgh_hij_matrix

!>routine to drive the application of the projector in HGH formalism
subroutine NL_HGH_application(hij,ncplx_p,n_p,wfd_p,proj,&
     ncplx_w,n_w,wfd_w,tolr,psi_pack,scpr,pdpsi,hpdpsi,psi,hpsi,eproj)
  use module_defs, only : gp,wp
  use psp_projectors, only: hgh_psp_application,nlpsp_to_wfd
  use locregs, only: wavefunctions_descriptors
  implicit none
  integer, intent(in) :: ncplx_p,n_p,ncplx_w,n_w
  !> interaction between the wavefuntion and the psp projector
  type(nlpsp_to_wfd), intent(in) :: tolr 
  type(wavefunctions_descriptors), intent(in) :: wfd_p !< descriptors of projectors
  type(wavefunctions_descriptors), intent(in) :: wfd_w !< descriptors of wavefunction
  !> matrix of nonlocal HGH psp
  real(gp), dimension(3,3,4), intent(in) :: hij
  real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,ncplx_p,n_p), intent(in) :: proj !< components of the projectors, real and imaginary parts
  real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,ncplx_w,n_w), intent(in) :: psi !< components of wavefunctions, real and imaginary parts

  real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,ncplx_w,n_w), intent(inout) :: hpsi !< components of wavefunctions, real and imaginary parts
  !> workspaces for the packing array
  real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_w*ncplx_w), intent(inout) :: psi_pack
  !> array of the scalar product between the projectors and the wavefunctions
  real(wp), dimension(ncplx_w,n_w,ncplx_p,n_p), intent(inout) :: scpr
  !> array of the coefficients of the hgh projectors
  real(wp), dimension(max(ncplx_w,ncplx_p),n_w,n_p), intent(inout) :: pdpsi
  !> array of the coefficients of the hgh projectors multiplied by HGH matrix
  real(wp), dimension(max(ncplx_w,ncplx_p),n_w,n_p), intent(inout) :: hpdpsi
  real(gp), intent(out) :: eproj

  call hgh_psp_application(hij,ncplx_p,n_p,wfd_p,proj,&
       ncplx_w,n_w,wfd_w,tolr,psi_pack,scpr,pdpsi,hpdpsi,psi,hpsi,eproj)

end subroutine NL_HGH_application
