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
  public :: bounds_to_plr_limits,pregion_size
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
    use locreg_operations
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
       call free_tolr_ptr(nl%pspd(iat)%tolr)
       call f_free_ptr(nl%pspd(iat)%lut_tolr)
       if (nl%pspd(iat)%mproj > 0) then
          !then fill it again, if the locreg is demanded
          nl%pspd(iat)%nlr=nlr
          call set_wfd_to_wfd(Glr,nl%pspd(iat)%plr,&
               keyg_lin,nbsegs_cf,nl%pspd(iat)%noverlap,nl%pspd(iat)%lut_tolr,nl%pspd(iat)%tolr,lrs,lr_mask)
       end if
    end do

    call f_free(keyg_lin)
    call f_free(nbsegs_cf)

    call f_release_routine()

  end subroutine update_nlpsp

  !> applies a projector of HGH type, written in a localization region
  !! onto a wavefunction written in the same formalism
  !! uses the desctiptors for the application which have been defined previously
  ! replace the routine nl_HGH_application as it does not need allocating arrays anymore
  subroutine hgh_psp_application(hij,ncplx_p,n_p,wfd_p,proj,&
       ncplx_w,n_w,wfd_w,tolr,psi_pack,scpr,pdpsi,hpdpsi,psi,hpsi,eproj)
    use pseudopotentials, only: apply_hij_coeff,atomic_proj_coeff
    use locreg_operations
    implicit none
    integer, intent(in) :: ncplx_p !< number of complex components of the projector
    integer, intent(in) :: n_p !< number of elements of the projector
    integer, intent(in) :: ncplx_w !< number of complex components of the wavefunction
    integer, intent(in) :: n_w !< number of complex components of the wavefunction
    type(wavefunctions_descriptors), intent(in) :: wfd_p !< descriptors of projectors
    type(wavefunctions_descriptors), intent(in) :: wfd_w !< descriptors of wavefunction
    !> interaction between the wavefuntion and the psp projector
    type(wfd_to_wfd), intent(in) :: tolr
    !> matrix of nonlocal HGH psp
    !real(gp), dimension(3,3,4), intent(in) :: hij
    type(atomic_proj_coeff), dimension(3,3,4) :: hij
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

    call pr_dot_psi(ncplx_p,n_p,wfd_p,proj,ncplx_w,n_w,wfd_w,psi,tolr,&
         psi_pack,scpr,pdpsi)

!!$    !put to zero the array
!!$    !call to_zero((wfd_p%nvctr_c+7*wfd_p%nvctr_f)*n_w*ncplx_w,psi_pack(1,1))
!!$    call f_zero(psi_pack)
!!$
!!$    !here also the strategy can be considered
!!$    call proj_dot_psi(n_p*ncplx_p,wfd_p,proj,n_w*ncplx_w,wfd_w,psi,&
!!$         tolr,psi_pack,scpr)
!!$
!!$    !first create the coefficients for the application of the matrix
!!$    !pdpsi = < p_i | psi >
!!$    call full_coefficients('C',ncplx_p,n_p,'N',ncplx_w,n_w,scpr,'N',pdpsi)

    call apply_hij_coeff(hij,max(ncplx_w,ncplx_p)*n_w,n_p,pdpsi,hpdpsi)

    call cproj_dot(ncplx_p,n_p,ncplx_w,n_w,scpr,pdpsi,hpdpsi,eproj)

!!$    !then create the coefficients for the evaluation of the projector energy
!!$    !pdpsi= < psi | p_i> = conj(< p_i | psi >)
!!$    call full_coefficients('N',ncplx_p,n_p,'C',ncplx_w,n_w,scpr,'C',pdpsi)
!!$
!!$    !the energy can be calculated here
!!$    eproj=dot(max(ncplx_p,ncplx_w)*n_w*n_p,pdpsi(1,1,1),1,hpdpsi(1,1,1),1)

    call cproj_pr_p_psi(hpdpsi,ncplx_p,n_p,wfd_p,proj,ncplx_w,n_w,wfd_w,hpsi,tolr,&
         psi_pack,scpr)

!!$    !then the coefficients have to be transformed for the projectors
!!$    call reverse_coefficients(ncplx_p,n_p,ncplx_w,n_w,hpdpsi,scpr)
!!$
!!$    call scpr_proj_p_hpsi(n_p*ncplx_p,wfd_p,proj,n_w*ncplx_w,wfd_w,&
!!$         tolr,psi_pack,scpr,hpsi)

  end subroutine hgh_psp_application
 
!!$  !> Calculate the scalar product with the projectors of a given set of 
!!$  !! orbitals (or support functions) given in the same localization region
!!$  subroutine calculate_cproj(ncplx_p,n_p,wfd_p,proj,&
!!$       ncplx_w,n_w,wfd_w,tolr,psi_pack,scpr,psi,pdpsi,hpdpsi)
!!$    use pseudopotentials, only: apply_hij_coeff
!!$    implicit none
!!$    integer, intent(in) :: ncplx_p !< number of complex components of the projector
!!$    integer, intent(in) :: n_p !< number of elements of the projector
!!$    integer, intent(in) :: ncplx_w !< number of complex components of the wavefunction
!!$    integer, intent(in) :: n_w !< number of complex components of the wavefunction
!!$    type(wavefunctions_descriptors), intent(in) :: wfd_p !< descriptors of projectors
!!$    type(wavefunctions_descriptors), intent(in) :: wfd_w !< descriptors of wavefunction
!!$    !> interaction between the wavefuntion and the psp projector
!!$    type(nlpsp_to_wfd), intent(in) :: tolr
!!$    !> components of the projectors, real and imaginary parts
!!$    real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,ncplx_p,n_p), intent(in) :: proj
!!$    !> components of wavefunctions, real and imaginary parts
!!$    real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,ncplx_w,n_w), intent(in) :: psi
!!$    !> workspaces for the packing array
!!$    real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_w*ncplx_w), intent(inout) :: psi_pack
!!$    !> array of the scalar product between the projectors and the wavefunctions
!!$    real(wp), dimension(ncplx_w,n_w,ncplx_p,n_p), intent(inout) :: scpr
!!$    !> array of the coefficients of the hgh projectors
!!$    real(wp), dimension(max(ncplx_w,ncplx_p),n_w,n_p), intent(inout) :: pdpsi
!!$    !> array of the coefficients of the hgh projectors multiplied by HGH matrix
!!$    real(wp), dimension(max(ncplx_w,ncplx_p),n_w,n_p), intent(inout) :: hpdpsi
!!$
!!$    !put to zero the array
!!$
!!$    call f_zero(psi_pack)
!!$
!!$    !here also the PSP application strategy can be considered
!!$    call proj_dot_psi(n_p*ncplx_p,wfd_p,proj,n_w*ncplx_w,wfd_w,psi,&
!!$         tolr%nmseg_c,tolr%nmseg_f,tolr%mask,psi_pack,scpr)
!!$
!!$    !first create the coefficients for the application of the matrix
!!$    !pdpsi = < p_i | psi >
!!$    call full_coefficients('C',ncplx_p,n_p,'N',ncplx_w,n_w,scpr,'N',pdpsi)
!!$
!!$    !then create the coefficients for the evaluation of the projector energy
!!$    !pdpsi= < psi | p_i> = conj(< p_i | psi >)
!!$    call full_coefficients('N',ncplx_p,n_p,'C',ncplx_w,n_w,scpr,'C',pdpsi)
!!$
!!$  end subroutine calculate_cproj


!!$    !call to_zero(max(ncplx_w,ncplx_p)*n_w*n_p,pdpsi(1,1,1))
!!$    pdpsi=0.0_wp

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
    use locreg_operations, only: wfd_to_wfd_skip
    implicit none
    ! Calling arguments
    integer,intent(in) :: iat, ilr
    type(locreg_descriptors),intent(in) :: llr, glr
    type(DFT_PSP_projectors),intent(in) :: nl
    logical :: overlap
    ! Local variables
    logical :: goon
    integer :: mproj, iilr
!!$ integer :: jlr

    overlap = .false.

    ! Check whether the projectors of this atom have an overlap with locreg ilr
    iilr=get_proj_locreg(nl,iat,ilr)
    goon=iilr/=0
    if (.not.goon) return

    mproj=nl%pspd(iat)%mproj
    !no projector on this atom
    if(mproj == 0) return
    if(wfd_to_wfd_skip(nl%pspd(iat)%tolr(iilr))) return

    call check_overlap(llr, nl%pspd(iat)%plr, glr, overlap)

    end function projector_has_overlap
  
  end module psp_projectors

!>routine to drive the application of the projector in HGH formalism
subroutine NL_HGH_application(hij,ncplx_p,n_p,wfd_p,proj,&
     ncplx_w,n_w,wfd_w,tolr,psi_pack,scpr,pdpsi,hpdpsi,psi,hpsi,eproj)
  use module_defs, only : gp,wp
  use psp_projectors, only: hgh_psp_application
  use locregs, only: wavefunctions_descriptors
  use locreg_operations, only: wfd_to_wfd
  implicit none
  integer, intent(in) :: ncplx_p,n_p,ncplx_w,n_w
  !> interaction between the wavefuntion and the psp projector
  type(wfd_to_wfd), intent(in) :: tolr 
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

!!$  call hgh_psp_application(hij,ncplx_p,n_p,wfd_p,proj,&
!!$       ncplx_w,n_w,wfd_w,tolr,psi_pack,scpr,pdpsi,hpdpsi,psi,hpsi,eproj)

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
