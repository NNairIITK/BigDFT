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
  use psp_projectors_base
  implicit none

  private

  public :: projector_has_overlap,get_proj_locreg
  public :: bounds_to_plr_limits,pregion_size,set_nlpsp_to_wfd
  public :: hgh_psp_application
  public :: update_nlpsp

  !routines which are typical of the projector application or creation follow
  contains

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
    use pseudopotentials, only: apply_hij_coeff
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
    !call to_zero((wfd_p%nvctr_c+7*wfd_p%nvctr_f)*n_w*ncplx_w,psi_pack(1,1))
    call f_zero(psi_pack)

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

  !> Performs the update of a set of wavefunctions with a projector each one written in Daubechies basis
  !! with its own descriptors.
  !! A masking array is used calculated to avoid the calculation of bitonic search for the scalar product
  !! If the number of projectors is bigger than 1 the wavefunction is also given by packing in the number of components
  !! of the projector to ease its successive application
  subroutine scpr_proj_p_hpsi(n_p,wfd_p,proj,n_w,wfd_w,&
       nmseg_c,nmseg_f,psi_mask,hpsi_pack,scpr,hpsi)
    implicit none
    integer, intent(in) :: n_p !< number of projectors (real and imaginary part included)
    integer, intent(in) :: n_w !< number of wavefunctions (real and imaginary part included)
    integer, intent(in) :: nmseg_c,nmseg_f !< segments of the masking array
    type(wavefunctions_descriptors), intent(in) :: wfd_p !< descriptors of projectors
    type(wavefunctions_descriptors), intent(in) :: wfd_w !< descriptors of wavefunction
    real(wp), dimension(n_w,n_p), intent(in) :: scpr !< array of the scalar product of all the components
    real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_p), intent(in) :: proj !< components of the projectors
    integer, dimension(3,nmseg_c+nmseg_f), intent(in) :: psi_mask !<lookup array in the wfn segments
    !indicating the points where data have to be taken for dot product
    ! always produced. Has to be initialized to zero first
    real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_w), intent(inout) :: hpsi_pack !< work array of hpsi in projector form
    !needed only when n_p is bigger than one 

    real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,n_w), intent(inout) :: hpsi !< wavefunction result
    !local variables
    logical, parameter :: mask=.false.,pack=.true.
    external :: waxpy_mask_unpack
    integer :: is_w,is_sw,is_p,is_sp,iw,is_sm

    is_w=wfd_w%nvctr_c+min(wfd_w%nvctr_f,1)
    is_sw=wfd_w%nseg_c+min(wfd_w%nseg_f,1)

    is_p=wfd_p%nvctr_c+min(wfd_p%nvctr_f,1)
    is_sp=wfd_p%nseg_c+min(wfd_p%nseg_f,1)

    is_sm=nmseg_c+min(nmseg_f,1)

    if (pack) then
       !once the coefficients are determined fill the components of the wavefunction with the last projector
       !linear algebra up to the second last projector
       !|psi_iw>=O_iw,jp| p_jp>

       if (n_p > 1) then
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
             call waxpy_mask_unpack(wfd_w%nvctr_c,wfd_w%nvctr_f,nmseg_c,nmseg_f,&
                  psi_mask(1,1),psi_mask(1,is_sm),hpsi_pack(1,iw),hpsi_pack(is_p,iw),&
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

    !call to_zero(max(ncplx_w,ncplx_p)*n_w*n_p,pdpsi(1,1,1))
    pdpsi=0.0_wp

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

  !> find the locreg that is associated to the given projector of atom iat
  !! for a locreg of label ilr. Shoudl the locreg not be found, the result is zero.
  function get_proj_locreg(nl,iat,ilr) result(iilr)
    implicit none
    integer, intent(in) :: iat,ilr
    type(DFT_PSP_projectors), intent(in) :: nl
    integer :: iilr
    !local variables
    integer :: jlr

    iilr=0
    do jlr=1,nl%pspd(iat)%noverlap
       if (nl%pspd(iat)%lut_tolr(jlr)==ilr) then
          iilr=jlr
          exit
       end if
    end do

  end function get_proj_locreg

  function projector_has_overlap(iat, ilr, llr, glr, nl) result(overlap)
    implicit none
    ! Calling arguments
    integer,intent(in) :: iat, ilr
    type(locreg_descriptors),intent(in) :: llr, glr
    type(DFT_PSP_projectors),intent(in) :: nl
    logical :: overlap
    ! Local variables
    logical :: goon
    integer :: mproj, jlr, iilr
  
    overlap = .false.
  
    ! Check whether the projectors of this atom have an overlap with locreg ilr
    iilr=get_proj_locreg(nl,iat,ilr)
    goon=iilr/=0
!!$    do jlr=1,nl%pspd(iat)%noverlap
!!$        if (nl%pspd(iat)%lut_tolr(jlr)==ilr) then
!!$            goon=.true.
!!$            iilr=jlr
!!$            exit
!!$        end if
!!$    end do
    if (.not.goon) return
  
    mproj=nl%pspd(iat)%mproj
    !no projector on this atom
    if(mproj == 0) return
    if(nl%pspd(iat)%tolr(iilr)%strategy == PSP_APPLY_SKIP) return
  
    call check_overlap(llr, nl%pspd(iat)%plr, glr, overlap)
  
    end function projector_has_overlap

    !> Performs the scalar product of a projector with a wavefunction each one writeen in Daubechies basis
    !! with its own descriptors.
    !! A masking array is then calculated to avoid the calculation of bitonic search for the scalar product
    !! If the number of projectors is bigger than 1 the wavefunction is also packed in the number of components
    !! of the projector to ease its successive application
    subroutine proj_dot_psi(n_p,wfd_p,proj,n_w,wfd_w,psi,nmseg_c,nmseg_f,psi_mask,psi_pack,scpr)
      implicit none
      integer, intent(in) :: n_p !< number of projectors (real and imaginary part included)
      integer, intent(in) :: n_w !< number of wavefunctions (real and imaginary part included)
      integer, intent(in) :: nmseg_c,nmseg_f !< segments of the masking array
      type(wavefunctions_descriptors), intent(in) :: wfd_p !< descriptors of projectors
      type(wavefunctions_descriptors), intent(in) :: wfd_w !< descriptors of wavefunction
      real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_p), intent(in) :: proj !< components of the projectors
      real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,n_w), intent(in) :: psi !< components of wavefunction
      integer, dimension(3,nmseg_c+nmseg_f), intent(in) :: psi_mask !<lookup array in the wfn segments
      !indicating the points where data have to be taken for dot product
      ! always produced. Has to be initialized to zero first
      real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_w), intent(inout) :: psi_pack !< packed array of psi in projector form
      !needed only when n_p is bigger than one 
      real(wp), dimension(n_w,n_p), intent(out) :: scpr !< array of the scalar product of all the components
      !local variables
      logical, parameter :: mask=.true.,pack=.true.
      integer :: is_w,is_sw,is_p,is_sp,iw,ip,is_sm
      !intensive routines
      external :: wpdot_keys_pack,wpdot_mask_pack

      !calculate starting points of the fine regions
      !they have to be calculated considering that there could be no fine grid points
      !therefore the array values should not go out of bounds even though their value is actually not used
      is_w=wfd_w%nvctr_c+min(wfd_w%nvctr_f,1)
      is_sw=wfd_w%nseg_c+min(wfd_w%nseg_f,1)

      is_p=wfd_p%nvctr_c+min(wfd_p%nvctr_f,1)
      is_sp=wfd_p%nseg_c+min(wfd_p%nseg_f,1)

      is_sm=nmseg_c+min(nmseg_f,1)

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
               call wpdot_mask_pack(wfd_w%nvctr_c,wfd_w%nvctr_f,nmseg_c,nmseg_f,&
                    psi_mask(1,1),psi_mask(1,is_sm),psi(1,iw),psi(is_w,iw),&
                    wfd_p%nvctr_c,wfd_p%nvctr_f,proj(1,1),proj(is_p,1),&
                    psi_pack(1,iw),psi_pack(is_p,iw),scpr(iw,1))
            end do
         end if

         !now that the packed array is constructed linear algebra routine can be used to calculate
         !use multithreaded dgemm or customized ones in the case of no OMP parallelized algebra
         !scpr(iw,ip) = < psi_iw| p_ip >
         if (n_p > 1) then
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
  
  end module psp_projectors

!>routine to drive the application of the projector in HGH formalism
subroutine NL_HGH_application(hij,ncplx_p,n_p,wfd_p,proj,&
     ncplx_w,n_w,wfd_w,tolr,psi_pack,scpr,pdpsi,hpdpsi,psi,hpsi,eproj)
  use module_defs, only : gp,wp
  use psp_projectors, only: hgh_psp_application
  use locregs, only: wavefunctions_descriptors
  use psp_projectors_base, only: nlpsp_to_wfd
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

!> Applies one real projector operator in the form |p> hp <p| onto a set of wavefunctions described by the same descriptors
!! accumulate the result on the array hpsi and calculate scpr @f$<p|psi_w>$@f such that energy can be expressed in the form @f$\sum_w <psi_w|p> hp <p|psi_w>@f$
subroutine apply_oneproj_operator(wfd_p,proj,hp,n_w,wfd_w,psi,hpsi,scpr)
  use module_base
  use module_types, only: wavefunctions_descriptors
  implicit none
  integer, intent(in) :: n_w !< complex components of the wavefunction
  real(wp), intent(in) :: hp !<coefficient of the projector operator
  type(wavefunctions_descriptors), intent(in) :: wfd_p !< descriptors of projectors
  type(wavefunctions_descriptors), intent(in) :: wfd_w !< descriptors of wavefunction
  !  real(gp), dimension(ncplx_o,ncomp_p,ncomp_p,ncomp_w), intent(in) :: hij !< matrix of operator in nonlocal projectors basis
  real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f), intent(in) :: proj !< components of the projector
  real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,n_w), intent(in) :: psi !< components of wavefunction
  real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,n_w), intent(inout) :: hpsi !<application of NL operator on psi
  real(wp), dimension(n_w), intent(out) :: scpr !<array of <p|psi_w>, to be used to evaluate energy terms
  !local variables
  character(len=*), parameter :: subname='apply_oneproj'
  integer :: is_w,is_sw,is_p,is_sp,iw
  integer, dimension(:,:), allocatable :: psi_mask
  !routines which are optimized in separate files
  external :: wpdot_keys,wpdot_mask,waxpy_mask

  call f_routine(id=subname)

  !calculate starting points of the fine regions
  !they have to be calculated considering that there could be no fine grid points
  !therefore the array values should not go out of bounds even though their value is actually not used
  is_w=wfd_w%nvctr_c+min(wfd_w%nvctr_f,1)
  is_sw=wfd_w%nseg_c+min(wfd_w%nseg_f,1)

  is_p=wfd_p%nvctr_c+min(wfd_p%nvctr_f,1)
  is_sp=wfd_p%nseg_c+min(wfd_p%nseg_f,1)

  !mask array to avoid multiple calls to bitonic search routines
  psi_mask=f_malloc0((/3,wfd_w%nseg_c+wfd_w%nseg_f/),id='psi_mask')
  call wpdot_keys(wfd_w%nvctr_c,wfd_w%nvctr_f,wfd_w%nseg_c,wfd_w%nseg_f,&
       wfd_w%keyvglob(1),wfd_w%keyvglob(is_sw),wfd_w%keyglob(1,1),wfd_w%keyglob(1,is_sw),&
       psi(1,1),psi(is_w,1),&
       wfd_p%nvctr_c,wfd_p%nvctr_f,wfd_p%nseg_c,wfd_p%nseg_f,&
       wfd_p%keyvglob(1),wfd_p%keyvglob(is_sp),wfd_p%keyglob(1,1),wfd_p%keyglob(1,is_sp),&
       proj(1),proj(is_p),&
       scpr(1))
  !use now mask arrays to calculate the rest of the scalar product
  do iw=2,n_w
  call wpdot_keys(wfd_w%nvctr_c,wfd_w%nvctr_f,wfd_w%nseg_c,wfd_w%nseg_f,&
       wfd_w%keyvglob(1),wfd_w%keyvglob(is_sw),&
       wfd_w%keyglob(1,1),wfd_w%keyglob(1,is_sw),&
       psi(1,iw),psi(is_w,iw),&
       wfd_p%nvctr_c,wfd_p%nvctr_f,wfd_p%nseg_c,wfd_p%nseg_f,&
       wfd_p%keyvglob(1),wfd_p%keyvglob(is_sp),&
       wfd_p%keyglob(1,1),wfd_p%keyglob(1,is_sp),&
       proj(1),proj(is_p),&
       scpr(iw))

!!$     call wpdot_mask(wfd_w%nvctr_c,wfd_w%nvctr_f,wfd_w%nseg_c,wfd_w%nseg_f,&
!!$          psi_mask(1,1),psi_mask(1,is_sw),psi(1,iw),psi(is_w,iw),&
!!$          wfd_p%nvctr_c,wfd_p%nvctr_f,proj(1),proj(is_p),scpr(iw))
  end do

  !then reduce the projector in the wavefunction
  do iw=1,n_w
     call waxpy(hp*scpr(iw),wfd_p%nvctr_c,wfd_p%nvctr_f,&
          wfd_p%nseg_c,wfd_p%nseg_f,&
          wfd_p%keyvglob(1),wfd_p%keyvglob(is_sp),&
          wfd_p%keyglob(1,1),wfd_p%keyglob(1,is_sp),proj(1),proj(is_p),& 
          wfd_w%nvctr_c,wfd_w%nvctr_f,wfd_w%nseg_c,wfd_w%nseg_f,&
          wfd_w%keyvglob(1),wfd_w%keyvglob(is_sw),&
          wfd_w%keyglob(1,1),wfd_w%keyglob(1,is_sw),&
          hpsi(1,iw),hpsi(is_w,iw))
!!$     call waxpy_mask(wfd_w%nvctr_c,wfd_w%nvctr_f,wfd_w%nseg_c,wfd_w%nseg_f,&
!!$          psi_mask(1,1),psi_mask(1,is_sw),hpsi(1,iw),hpsi(is_w,iw),&
!!$          wfd_p%nvctr_c,wfd_p%nvctr_f,proj(1),proj(is_p),&
!!$          hp*scpr(iw))
  end do

  call f_free(psi_mask)

  call f_release_routine()

end subroutine apply_oneproj_operator
