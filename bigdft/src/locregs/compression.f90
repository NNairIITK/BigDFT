!> @file
!! Datatypes and associated methods relative to the descriptors for the compressions wavelet descriptions
!! @author
!!    Copyright (C) 2016-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!> Datatypes for wavefunction descriptors
module compression
  use module_defs, only: wp
  implicit none

  private 

  !> Parameters identifying the different strategy for the application of a projector
  !! in a localisation region
  integer, parameter :: STRATEGY_SKIP=0 !< The projector is not applied. This might happend when ilr and iat does not interact
  integer, parameter :: STRATEGY_MASK=1         !< Use mask arrays. The mask array has to be created before.
  integer, parameter :: STRATEGY_KEYS=2         !< Use keys. No mask nor packing. Equivalend to traditional application
  integer, parameter :: STRATEGY_MASK_PACK=3    !< Use masking and creates a pack arrays from them.
  !! Most likely this is the common usage for atoms
  !! with lots of projectors and localization regions "close" to them
  integer, parameter :: STRATEGY_KEYS_PACK=4    !< Use keys and pack arrays. Useful especially when there is no memory to create a lot of packing arrays,
  !! for example when lots of lrs interacts with lots of atoms


  !> Used for lookup table for compressed wavefunctions
  type, public :: wavefunctions_descriptors
     integer :: nvctr_c,nvctr_f,nseg_c,nseg_f
     integer, dimension(:,:), pointer :: keyglob
     integer, dimension(:,:), pointer :: keygloc
     integer, dimension(:), pointer :: keyvloc,keyvglob
  end type wavefunctions_descriptors

  !> arrays defining how a given projector and a given wavefunction descriptor should interact
  type, public :: wfd_to_wfd
     integer :: strategy !< can be STRATEGY_MASK,STRATEGY_KEYS,STRATEGY_MASK_PACK,STRATEGY_KEYS_PACK,STRATEGY_SKIP
     integer :: nmseg_c !< number of segments intersecting in the coarse region
     integer :: nmseg_f !< number of segments intersecting in the fine region
     integer, dimension(:,:), pointer :: mask !<mask array of dimesion 3,nmseg_c+nmseg_f for psp application
  end type wfd_to_wfd


  public :: allocate_wfd,deallocate_wfd,copy_wavefunctions_descriptors
  public :: deallocate_wfd_to_wfd,nullify_wfd
  public :: nullify_wfd_to_wfd,tolr_set_strategy
  public :: cproj_dot,cproj_pr_p_psi,pr_dot_psi
  public :: wfd_to_wfd_skip,free_tolr_ptr,init_tolr

contains

  pure function wfd_null() result(wfd)
    implicit none
    type(wavefunctions_descriptors) :: wfd
    call nullify_wfd(wfd)
  end function wfd_null

  pure subroutine nullify_wfd(wfd)
    implicit none
    type(wavefunctions_descriptors), intent(out) :: wfd
    wfd%nvctr_c=0
    wfd%nvctr_f=0
    wfd%nseg_c=0
    wfd%nseg_f=0
    nullify(wfd%keyglob)
    nullify(wfd%keygloc)
    nullify(wfd%keyvglob)
    nullify(wfd%keyvloc)
  end subroutine nullify_wfd

  !creators
  pure function wfd_to_wfd_null() result(tolr)
    implicit none
    type(wfd_to_wfd) :: tolr
    call nullify_wfd_to_wfd(tolr)
  end function wfd_to_wfd_null
  pure subroutine nullify_wfd_to_wfd(tolr)
    implicit none
    type(wfd_to_wfd), intent(out) :: tolr
    tolr%strategy=STRATEGY_SKIP
    tolr%nmseg_c=0
    tolr%nmseg_f=0
    nullify(tolr%mask)
  end subroutine nullify_wfd_to_wfd

  !destructors
  subroutine deallocate_wfd_to_wfd(tolr)
    use dynamic_memory
    implicit none
    type(wfd_to_wfd), intent(inout) :: tolr
    call f_free_ptr(tolr%mask)
  end subroutine deallocate_wfd_to_wfd

  subroutine free_tolr_ptr(tolr)
    implicit none
    type(wfd_to_wfd), dimension(:), pointer :: tolr
    !local variables
    integer :: ilr
    if (.not. associated(tolr)) return
    do ilr=lbound(tolr,1),ubound(tolr,1)
       call deallocate_wfd_to_wfd(tolr(ilr))
       call nullify_wfd_to_wfd(tolr(ilr))
    end do
    deallocate(tolr)
    nullify(tolr)
  end subroutine free_tolr_ptr

  pure function wfd_to_wfd_skip(tolr)
    implicit none
    type(wfd_to_wfd), intent(in) :: tolr
    logical :: wfd_to_wfd_skip

    wfd_to_wfd_skip = tolr%strategy == STRATEGY_SKIP
  end function wfd_to_wfd_skip


  !> here we should have already defined the number of segments
  subroutine allocate_wfd(wfd)
    use dynamic_memory
    implicit none
    type(wavefunctions_descriptors), intent(inout) :: wfd
    !local variables
    integer :: nsegs

    nsegs=max(1,wfd%nseg_c+wfd%nseg_f)
    wfd%keyvloc=f_malloc_ptr(nsegs,id='wfd%keyvloc')
    wfd%keyvglob=f_malloc_ptr(nsegs,id='wfd%keyvglob')
    wfd%keyglob=f_malloc_ptr((/2,nsegs/),id='wfd%keyglob')
    wfd%keygloc=f_malloc_ptr((/2,nsegs/),id='wfd%keygloc')
  END SUBROUTINE allocate_wfd

  !> De-Allocate wavefunctions_descriptors
  subroutine deallocate_wfd(wfd)
    use dynamic_memory
    implicit none
    type(wavefunctions_descriptors), intent(inout) :: wfd

    !in case the two objects points to the same target
    if (associated(wfd%keyglob, target = wfd%keygloc)) then
       !assuming that globals has been created afterwards
       nullify(wfd%keygloc)
       call f_free_ptr(wfd%keyglob)
    else
       call f_free_ptr(wfd%keygloc)
       call f_free_ptr(wfd%keyglob)
    end if
    if (associated(wfd%keyvloc, target= wfd%keyvglob)) then
       nullify(wfd%keyvloc)
       call f_free_ptr(wfd%keyvglob)
    else
       call f_free_ptr(wfd%keyvloc)
       call f_free_ptr(wfd%keyvglob)
    end if
  END SUBROUTINE deallocate_wfd


  subroutine copy_wavefunctions_descriptors(wfdin, wfdout)
    use dynamic_memory
    implicit none
    ! Calling arguments
    type(wavefunctions_descriptors), intent(in) :: wfdin
    type(wavefunctions_descriptors), intent(out) :: wfdout

    ! Local variables
    !integer:: istat,iis1, iie1, iis2, iie2,i1, i2, iall

    !nullify all pointers first
    call nullify_wfd(wfdout)

    wfdout%nvctr_c = wfdin%nvctr_c
    wfdout%nvctr_f = wfdin%nvctr_f
    wfdout%nseg_c = wfdin%nseg_c
    wfdout%nseg_f = wfdin%nseg_f

    !new method
    wfdout%keygloc=f_malloc_ptr(src_ptr=wfdin%keygloc,id='wfdout%keygloc')
    wfdout%keyglob=f_malloc_ptr(src_ptr=wfdin%keyglob,id='wfdout%keyglob')
    wfdout%keyvloc=f_malloc_ptr(src_ptr=wfdin%keyvloc,id='wfdout%keyvloc')
    wfdout%keyvglob=f_malloc_ptr(src_ptr=wfdin%keyvglob,id='wfdout%keyvglob')

!!$    !no need to insert lbounds as the allocation start from 1
!!$    if (associated(wfdin%keygloc)) wfdout%keygloc=f_malloc_ptr(src=wfdin%keygloc,id='wfdout%keygloc')
!!$    if (associated(wfdin%keyglob)) wfdout%keyglob=f_malloc_ptr(src=wfdin%keyglob,id='wfdout%keyglob')
!!$    if (associated(wfdin%keyvloc)) wfdout%keyvloc=f_malloc_ptr(src=wfdin%keyvloc,id='wfdout%keyvloc')
!!$    if (associated(wfdin%keyvglob))wfdout%keyvglob=f_malloc_ptr(src=wfdin%keyvglob,id='wfdout%keyvglob')


  end subroutine copy_wavefunctions_descriptors

  !> initialize the wfd_to_wfd descriptor starting from 
  !! the descriptors of the localization regions
  subroutine init_tolr(tolr,wfd_lr,wfd_p,keyag_lin_cf,nbsegs_cf)
    use dynamic_memory
    use wrapper_linalg
    use f_utils
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
    type(wfd_to_wfd), intent(inout) :: tolr

    call f_routine(id='init_tolr')

    !calculate the size of the mask array
    call vcopy(wfd_lr%nseg_c+wfd_lr%nseg_f,&
         wfd_lr%keyglob(1,1),2,keyag_lin_cf(1),1)
    call f_zero(nbsegs_cf)
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

  subroutine tolr_set_strategy(tolr,strategy)
    use dictionaries, only: f_err_throw
    implicit none
    character(len=*), intent(in) :: strategy
    type(wfd_to_wfd), intent(inout) :: tolr
    select case(trim(strategy))
    case('MASK_PACK','mask_pack')
       tolr%strategy=STRATEGY_MASK_PACK
    case('MASK','mask')
       tolr%strategy=STRATEGY_MASK
    case('KEYS','keys')
       tolr%strategy=STRATEGY_KEYS
    case('KEYS_PACK','keys_pack')
       tolr%strategy=STRATEGY_KEYS_PACK
    case default
       call f_err_throw('Unknown wfd_to_wfd strategy')
    end select
    
  end subroutine tolr_set_strategy

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

  subroutine pr_dot_psi(ncplx_p,n_p,wfd_p,pr,ncplx_w,n_w,wfd_w,psi,tolr,&
       wpack,scpr,cproj)
    use f_utils
    implicit none
    integer, intent(in) :: ncplx_p !< number of complex components of the projector
    integer, intent(in) :: n_p !< number of elements of the projector
    integer, intent(in) :: ncplx_w !< number of complex components of the wavefunction
    integer, intent(in) :: n_w !< number of complex components of the wavefunction
    type(wavefunctions_descriptors), intent(in) :: wfd_p !< descriptors of projectors
    type(wavefunctions_descriptors), intent(in) :: wfd_w !< descriptors of wavefunction
    !> interaction between the wavefuntion and the psp projector
    type(wfd_to_wfd), intent(in) :: tolr
    !> components of the projectors, real and imaginary parts
    real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,ncplx_p,n_p), intent(in) :: pr
    !> components of wavefunctions, real and imaginary parts
    real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,ncplx_w,n_w), intent(in) :: psi 
    !> workspaces for the packing array
    real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_w*ncplx_w), intent(out) :: wpack
    !> array of the scalar product between the projectors and the wavefunctions
    real(wp), dimension(ncplx_w,n_w,ncplx_p,n_p), intent(out) :: scpr
    !> array of the coefficients of the hgh projectors
    real(wp), dimension(max(ncplx_w,ncplx_p),n_w,n_p), intent(out) :: cproj

    call f_zero(wpack)
    !here also the strategy can be considered
    call proj_dot_psi(n_p*ncplx_p,wfd_p,pr,n_w*ncplx_w,wfd_w,psi,&
         tolr,wpack,scpr)
    !first create the coefficients for the application of the matrix
    !pdpsi = < p_i | psi >
    call full_coefficients('C',ncplx_p,n_p,'N',ncplx_w,n_w,scpr,'N',cproj)

  end subroutine pr_dot_psi

  !>perform the scalar product between two cproj arrays.
  !useful to calculate the nl projector energy
  subroutine cproj_dot(ncplx_p,n_p,ncplx_w,n_w,scpr,a,b,eproj)
    use wrapper_linalg, only: dot
    implicit none
    integer, intent(in) :: ncplx_p !< number of complex components of the projector
    integer, intent(in) :: n_p !< number of elements of the projector
    integer, intent(in) :: ncplx_w !< number of complex components of the wavefunction
    integer, intent(in) :: n_w !< number of complex components of the wavefunction
    !> array of the scalar product between the projectors and the wavefunctions
    real(wp), dimension(ncplx_w,n_w,ncplx_p,n_p), intent(in) :: scpr
    real(wp), dimension(max(ncplx_w,ncplx_p),n_w,n_p), intent(inout) :: a !<cproj, conjugated in output
    real(wp), dimension(max(ncplx_w,ncplx_p),n_w,n_p), intent(in) :: b !<cproj arrays
    real(wp), intent(out) :: eproj

    !then create the coefficients for the evaluation of the projector energy
    !pdpsi= < psi | p_i> = conj(< p_i | psi >)
    call full_coefficients('N',ncplx_p,n_p,'C',ncplx_w,n_w,scpr,'C',a)

    !the energy can be calculated here
    eproj=dot(max(ncplx_p,ncplx_w)*n_w*n_p,a(1,1,1),1,b(1,1,1),1)
  end subroutine cproj_dot

  !> update of the psi 
  subroutine cproj_pr_p_psi(cproj,ncplx_p,n_p,wfd_p,pr,ncplx_w,n_w,wfd_w,psi,tolr,&
       wpack,scpr)
    implicit none
    integer, intent(in) :: ncplx_p !< number of complex components of the projector
    integer, intent(in) :: n_p !< number of elements of the projector
    integer, intent(in) :: ncplx_w !< number of complex components of the wavefunction
    integer, intent(in) :: n_w !< number of complex components of the wavefunction
    type(wavefunctions_descriptors), intent(in) :: wfd_p !< descriptors of projectors
    type(wavefunctions_descriptors), intent(in) :: wfd_w !< descriptors of wavefunction
    !> interaction between the wavefuntion and the psp projector
    type(wfd_to_wfd), intent(in) :: tolr
    !> array of the coefficients of the hgh projectors
    real(wp), dimension(max(ncplx_w,ncplx_p),n_w,n_p), intent(in) :: cproj
    !> components of the projectors, real and imaginary parts
    real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,ncplx_p,n_p), intent(in) :: pr
    !> components of wavefunctions, real and imaginary parts
    real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,ncplx_w,n_w), intent(inout) :: psi 
    !> workspaces for the packing array
    real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_w*ncplx_w), intent(out) :: wpack
    !> array of the scalar product between the projectors and the wavefunctions
    real(wp), dimension(ncplx_w,n_w,ncplx_p,n_p), intent(out) :: scpr

    !then the coefficients have to be transformed for the projectors
    call reverse_coefficients(ncplx_p,n_p,ncplx_w,n_w,cproj,scpr)

    call scpr_proj_p_hpsi(n_p*ncplx_p,wfd_p,pr,n_w*ncplx_w,wfd_w,&
         tolr,wpack,scpr,psi)
  end subroutine cproj_pr_p_psi


  !> Performs the scalar product of a projector with a wavefunction each one writeen in Daubechies basis
  !! with its own descriptors.
  !! A masking array is then calculated to avoid the calculation of bitonic search for the scalar product
  !! If the number of projectors is bigger than 1 the wavefunction is also packed in the number of components
  !! of the projector to ease its successive application
  subroutine proj_dot_psi(n_p,wfd_p,proj,n_w,wfd_w,psi,tolr,psi_pack,scpr)
    use wrapper_linalg
    implicit none
    integer, intent(in) :: n_p !< number of projectors (real and imaginary part included)
    integer, intent(in) :: n_w !< number of wavefunctions (real and imaginary part included)
!!$      integer, intent(in) :: nmseg_c,nmseg_f !< segments of the masking array
    type(wavefunctions_descriptors), intent(in) :: wfd_p !< descriptors of projectors
    type(wavefunctions_descriptors), intent(in) :: wfd_w !< descriptors of wavefunction
    real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_p), intent(in) :: proj !< components of the projectors
    real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,n_w), intent(in) :: psi !< components of wavefunction
    type(wfd_to_wfd), intent(in) :: tolr !< datatype for strategy information
!!$      integer, dimension(3,nmseg_c+nmseg_f), intent(in) :: tolr%mask !<lookup array in the wfn segments
!!$      !indicating the points where data have to be taken for dot product
!!$      ! always produced. Has to be initialized to zero first
    real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_w), intent(inout) :: psi_pack !< packed array of psi in projector form
    !needed only when n_p is bigger than one 
    real(wp), dimension(n_w,n_p), intent(out) :: scpr !< array of the scalar product of all the components
    !local variables
    logical :: mask,pack!, parameter :: mask=.true.,pack=.true.
    integer :: is_w,is_sw,is_p,is_sp,iw,ip,is_sm
    !intensive routines
    external :: wpdot_keys_pack,wpdot_mask_pack

    if (tolr%strategy==STRATEGY_SKIP) return

    !calculate starting points of the fine regions
    !they have to be calculated considering that there could be no fine grid points
    !therefore the array values should not go out of bounds even though their value is actually not used
    is_w=wfd_w%nvctr_c+min(wfd_w%nvctr_f,1)
    is_sw=wfd_w%nseg_c+min(wfd_w%nseg_f,1)

    is_p=wfd_p%nvctr_c+min(wfd_p%nvctr_f,1)
    is_sp=wfd_p%nseg_c+min(wfd_p%nseg_f,1)

    is_sm=tolr%nmseg_c+min(tolr%nmseg_f,1)

    pack=(tolr%strategy==STRATEGY_MASK_PACK) .or. (tolr%strategy==STRATEGY_KEYS_PACK)
    mask=(tolr%strategy==STRATEGY_MASK_PACK) .or. (tolr%strategy==STRATEGY_MASK)

    if (pack) then
       if (.not. mask) then
          do iw=1,n_w
             call wpdot_keys_pack(wfd_w%nvctr_c,wfd_w%nvctr_f,wfd_w%nseg_c,wfd_w%nseg_f,&
                  wfd_w%keyvglob(1),wfd_w%keyvglob(is_sw),wfd_w%keyglob(1,1),wfd_w%keyglob(1,is_sw),&
                  psi(1,iw),psi(is_w,iw),&
                  wfd_p%nvctr_c,wfd_p%nvctr_f,wfd_p%nseg_c,wfd_p%nseg_f,&
                  wfd_p%keyvglob(1),wfd_p%keyvglob(is_sp),wfd_p%keyglob(1,1),wfd_p%keyglob(1,is_sp),&
                  proj(1,1),proj(is_p,1),&
                  psi_pack(1,iw),psi_pack(is_p,iw),scpr(iw,1))
          end do
       else 
          do iw=1,n_w
             call wpdot_mask_pack(wfd_w%nvctr_c,wfd_w%nvctr_f,tolr%nmseg_c,tolr%nmseg_f,&
                  tolr%mask(1,1),tolr%mask(1,is_sm),psi(1,iw),psi(is_w,iw),&
                  wfd_p%nvctr_c,wfd_p%nvctr_f,proj(1,1),proj(is_p,1),&
                  psi_pack(1,iw),psi_pack(is_p,iw),scpr(iw,1))
          end do
       end if

       !now that the packed array is constructed linear algebra routine can be used to calculate
       !use multithreaded dgemm or customized ones in the case of no OMP parallelized algebra
       !scpr(iw,ip) = < psi_iw| p_ip >
       if (n_p > 1) then
!!$            call f_gemm(trans_a='T',a=psi_pack,&
!!$                 shape_b=[wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_p-1],b=proj(1,2),&
!!$                 c=scpr(1,2),shape_c=[n_w,n_p-1])

          call gemm('T','N',n_w,n_p-1,wfd_p%nvctr_c+7*wfd_p%nvctr_f,1.0_wp,psi_pack(1,1),&
               wfd_p%nvctr_c+7*wfd_p%nvctr_f,proj(1,2),wfd_p%nvctr_c+7*wfd_p%nvctr_f,0.0_wp,&
               scpr(1,2),n_w)
       end if

    else
       do ip=1,n_p
          do iw=1,n_w
             call wpdot_keys(wfd_w%nvctr_c,wfd_w%nvctr_f,wfd_w%nseg_c,wfd_w%nseg_f,&
                  wfd_w%keyvglob(1),wfd_w%keyvglob(is_sw),wfd_w%keyglob(1,1),wfd_w%keyglob(1,is_sw),&
                  psi(1,iw),psi(is_w,iw),&
                  wfd_p%nvctr_c,wfd_p%nvctr_f,wfd_p%nseg_c,wfd_p%nseg_f,&
                  wfd_p%keyvglob(1),wfd_p%keyvglob(is_sp),wfd_p%keyglob(1,1),wfd_p%keyglob(1,is_sp),&
                  proj(1,ip),proj(is_p,ip),&
                  scpr(iw,ip))
          end do
       end do
    end if
  end subroutine proj_dot_psi


  !> Performs the update of a set of wavefunctions with a projector each one written in Daubechies basis
  !! with its own descriptors.
  !! A masking array is used calculated to avoid the calculation of bitonic search for the scalar product
  !! If the number of projectors is bigger than 1 the wavefunction is also given by packing in the number of components
  !! of the projector to ease its successive application
  subroutine scpr_proj_p_hpsi(n_p,wfd_p,proj,n_w,wfd_w,tolr,hpsi_pack,scpr,hpsi)
    use wrapper_linalg
    use f_utils
    implicit none
    integer, intent(in) :: n_p !< number of projectors (real and imaginary part included)
    integer, intent(in) :: n_w !< number of wavefunctions (real and imaginary part included)
!!$      integer, intent(in) :: nmseg_c,nmseg_f !< segments of the masking array
    type(wavefunctions_descriptors), intent(in) :: wfd_p !< descriptors of projectors
    type(wavefunctions_descriptors), intent(in) :: wfd_w !< descriptors of wavefunction
    real(wp), dimension(n_w,n_p), intent(in) :: scpr !< array of the scalar product of all the components
    real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_p), intent(in) :: proj !< components of the projectors
    type(wfd_to_wfd), intent(in) :: tolr
!!$      integer, dimension(3,nmseg_c+nmseg_f), intent(in) :: psi_mask !<lookup array in the wfn segments
    !indicating the points where data have to be taken for dot product
    ! always produced. Has to be initialized to zero first
    real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_w), intent(inout) :: hpsi_pack !< work array of hpsi in projector form
    !needed only when n_p is bigger than one 

    real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,n_w), intent(inout) :: hpsi !< wavefunction result
    !local variables
    logical :: mask,pack !parameter :: mask=.false.,pack=.true.
    external :: waxpy_mask_unpack
    integer :: is_w,is_sw,is_p,is_sp,iw,is_sm

    if (tolr%strategy==STRATEGY_SKIP) return


    is_w=wfd_w%nvctr_c+min(wfd_w%nvctr_f,1)
    is_sw=wfd_w%nseg_c+min(wfd_w%nseg_f,1)

    is_p=wfd_p%nvctr_c+min(wfd_p%nvctr_f,1)
    is_sp=wfd_p%nseg_c+min(wfd_p%nseg_f,1)

    is_sm=tolr%nmseg_c+min(tolr%nmseg_f,1)

    pack=(tolr%strategy==STRATEGY_MASK_PACK) .or. (tolr%strategy==STRATEGY_KEYS_PACK)
    mask=(tolr%strategy==STRATEGY_MASK_PACK) .or. (tolr%strategy==STRATEGY_MASK)


    if (pack) then
       !once the coefficients are determined fill the components of the wavefunction with the last projector
       !linear algebra up to the second last projector
       !|psi_iw>=O_iw,jp| p_jp>

       if (n_p > 1) then
!!$            call f_gemm(a=proj(1,1),shape_a=[wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_p-1],&
!!$                 b=scpr(1,1),shape_b=[n_w,n_p-1],trans_b='T',c=hpsi_pack)

          call gemm('N','T',wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_w,n_p-1,&
               1.0_wp,proj(1,1),wfd_p%nvctr_c+7*wfd_p%nvctr_f,&
               scpr(1,1),n_w,0.0_wp,&
               hpsi_pack(1,1),wfd_p%nvctr_c+7*wfd_p%nvctr_f)
       else
          call f_zero(hpsi_pack)
       end if

       !then last projector
       if (mask) then
          do iw=1,n_w
             call waxpy_mask_unpack(wfd_w%nvctr_c,wfd_w%nvctr_f,tolr%nmseg_c,tolr%nmseg_f,&
                  tolr%mask(1,1),tolr%mask(1,is_sm),hpsi_pack(1,iw),hpsi_pack(is_p,iw),&
                  hpsi(1,iw),hpsi(is_w,iw),&
                  wfd_p%nvctr_c,wfd_p%nvctr_f,proj(1,n_p),proj(is_p,n_p),&
                  scpr(iw,n_p))
          end do
       else
          do iw=1,n_w
             call waxpy_keys_unpack(wfd_w%nvctr_c,wfd_w%nvctr_f,wfd_w%nseg_c,wfd_w%nseg_f,&
                  wfd_w%keyvglob(1),wfd_w%keyvglob(is_sw),wfd_w%keyglob(1,1),wfd_w%keyglob(1,is_sw),&
                  hpsi(1,iw),hpsi(is_w,iw),&
                  wfd_p%nvctr_c,wfd_p%nvctr_f,wfd_p%nseg_c,wfd_p%nseg_f,&
                  wfd_p%keyvglob(1),wfd_p%keyvglob(is_sp),&
                  wfd_p%keyglob(1,1),wfd_p%keyglob(1,is_sp),&
                  proj(1,n_p),proj(is_p,n_p),&
                  hpsi_pack(1,iw),hpsi_pack(is_p,iw),scpr(iw,n_p))
          end do
       end if
    else

    end if

  end subroutine scpr_proj_p_hpsi

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


end module compression
