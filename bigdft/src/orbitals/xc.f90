!> @file
!!  Wrapper around XC library routines (both ABINIT and LibXC).
!! @author
!!    Copyright (C) 2008-2010 ABINIT group (MOliveira)
!!    Copyright (C) 2008-2013 BigDFT group (DC)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module handling the eXchange-Correlation functionals
!! using ABINIT (old) and libxc library
module module_xc

  use module_base
  use xc_f90_types_m
  use libxc_funcs_m
  use xc_f90_lib_m
  use yaml_output
  use abi_interfaces_xc_lowlevel, only: abi_drivexc,abi_size_dvxc

  implicit none

  integer, public, parameter :: XC_ABINIT = 1  !< xc functionals from ABINIT
  integer, public, parameter :: XC_LIBXC  = 2  !< xc from libxc
  integer, public, parameter :: XC_MIXED  = 3  !> xc mixing the origin of the xc functionals

  integer, public, parameter :: XC_HARTREE = 0        !< IXC code for Hartree
  integer, public, parameter :: XC_HARTREE_FOCK = 100 !< IXC code for Hartree-Fock
  integer, public, parameter :: XC_NO_HARTREE = 1000  !< IXC code for no Hartree and no functional (ex. hydrogen atom)


  !> Structure containing the information to call the routines for the calculation of the xc functionals
  type, public :: libxc_functional
     !private
     type(xc_f90_pointer_t) :: conf !< the pointer used to call the library
     type(xc_f90_pointer_t) :: info !< Information about the functional
  end type libxc_functional

  !> Structure containing the information about the xc functionals
  type xc_info
     integer :: ixc        !< input XC code
     integer :: kind       !< ABINIT or LibXC
     integer :: family(2)  !< LDA, GGA, etc.
     integer :: id(2)      !< Identifier
     real(kind=8) :: default_alpha !< alpha set by input
     type(libxc_functional) :: funcs(2)
  end type xc_info

  !type(xc_info) :: xc      !< Global structure about the used xc functionals

  logical :: abinit_init = .false.                                   !< .True. if already ABINIT_XC_NAMES intialized
  !logical :: libxc_init = .false.                                   !< .True. if the libXC library has been initialized
  integer, parameter :: ABINIT_N_NAMES = 29                          !< Number of names for ABINIT_XC_NAMES
  character(len=500), dimension(0:ABINIT_N_NAMES) :: ABINIT_XC_NAMES !< Names of the xc functionals used by ABINIT

  private

  !> Public routines
  public :: xc_info, &
       &    xc_init, &
       &    xc_dump, &
       &    xc_init_rho, &
       &    xc_clean_rho, &
       &    xc_getvxc, &
       &    xc_isgga, &
       &    xc_exctXfac, &
       &    xc_end, &
       &    xc_get_name, &
       &    xc_get_id_from_name
  public :: xc_energy_new



contains

  subroutine obj_init_(xcObj, ixc, kind, nspden, default_alpha)
    implicit none

    !Arguments
    !scalars
    type(xc_info), intent(out) :: xcObj
    integer, intent(in) :: nspden
    integer, intent(in) :: kind
    integer, intent(in) :: ixc
    real(kind=8), intent(in), optional :: default_alpha

    !Local variables
    !scalars
    integer :: i, ierr

    ! **

    xcObj%ixc  = ixc
    xcObj%kind = kind
    if (kind == XC_LIBXC .or. kind == XC_MIXED) then
       xcObj%id(2) = abs(ixc) / 1000
       xcObj%id(1) = abs(ixc) - xcObj%id(2) * 1000
    else if (kind == XC_ABINIT) then
       ! ABINIT functional.
       xcObj%id(1) = abs(ixc)
       xcObj%id(2) = 0
    end if

    xcObj%family = 0
    if (xcObj%kind == XC_LIBXC .or. xcObj%kind == XC_MIXED) then
       ! LibXC case.

       do i = 1, 2
          if (xcObj%id(i) == 0) cycle
          ! Get XC functional family
          xcObj%family(i) = xc_f90_family_from_id(xcObj%id(i))
          select case (xcObj%family(i))
          case (XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
             call xc_f90_func_init(xcObj%funcs(i)%conf,xcObj%funcs(i)%info,xcObj%id(i),nspden)
          case default
             write(*,*) "Error: unsupported functional, change ixc."
             call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierr)
          end select
       end do

    else if (xcObj%kind == XC_ABINIT) then
       ! ABINIT case

       ! Get XC functional family
       if ((xcObj%id(1) >= 1 .and. xcObj%id(1) < 11) .or. xcObj%id(1) == 24) then
          xcObj%family(1) = XC_FAMILY_LDA
       else if (xcObj%id(1) >= 11 .and. xcObj%id(1) < 18) then
          xcObj%family(1) = XC_FAMILY_GGA
       else if (xcObj%id(1) >= 23 .and. xcObj%id(1) < 28) then
          xcObj%family(1) = XC_FAMILY_GGA
       else if (xcObj%id(1) >= 31 .and. xcObj%id(1) < 35) then
          xcObj%family(1) = XC_FAMILY_LDA
       else if (xcObj%id(1) == XC_HARTREE .or. xcObj%id(1) == XC_HARTREE_FOCK .or. xcObj%id(1) == XC_NO_HARTREE) then
          xcObj%family(1) = 0
       else
          write(*,*) "Error: unsupported functional, change ixc."
          call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierr)
       end if
    end if

    if(present(default_alpha)) then
      xcObj%default_alpha=default_alpha
    else
      xcObj%default_alpha=-1.0d0
    end if
    
  end subroutine obj_init_

  subroutine obj_free_(xcObj)
    implicit none

    type(xc_info), intent(inout) :: xcObj

    if (xcObj%kind == XC_LIBXC .or. xcObj%kind == XC_MIXED) then
       if (xcObj%id(1) > 0) call xc_f90_func_end(xcObj%funcs(1)%conf)
       if (xcObj%id(2) > 0) call xc_f90_func_end(xcObj%funcs(2)%conf)
    end if

  end subroutine obj_free_


  !> Give the name of the XC functional
  subroutine obj_get_name_(xcObj, name)
    implicit none
    !Arguments
    type(xc_info), intent(in) :: xcObj      !< XC objects
    character(len=500), intent(out) :: name !< XC functional name
    !Local variables
    integer :: i
    character(len = 500) :: messX, messC, messXC

    if (xcObj%kind == XC_LIBXC .or. xcObj%kind == XC_MIXED) then
       write(messX, "(A)") ""
       write(messC, "(A)") ""
       write(messXC, "(A)") ""
       do i = 1, 2
          if (xcObj%family(i) == 0) cycle
          select case(xc_f90_info_kind(xcObj%funcs(i)%info))
          case (XC_EXCHANGE)
             call xc_f90_info_name(xcObj%funcs(i)%info, messX)
          case (XC_CORRELATION)
             call xc_f90_info_name(xcObj%funcs(i)%info, messC)
          case (XC_EXCHANGE_CORRELATION)
             call xc_f90_info_name(xcObj%funcs(i)%info, messXC)
          end select
       end do
       if (len(trim(messXC)) > 0) then
          write(name, "(2A)")    "XC: ", trim(messXC)
       else
          if (len(trim(messC)) == 0) then
             write(name, "(2A)") "X-: ", trim(messX)
          else if (len(trim(messX)) == 0) then
             write(name, "(2A)") "-C: ", trim(messC)
          else
             if (trim(messX) == trim(messC)) then
                write(name, "(2A)") "XC: ", trim(messX)
             else
                write(name, "(4A)") "XC: ", trim(messX), ", ", trim(messC)
             end if
          end if
       end if
    else if (xcObj%kind == XC_ABINIT) then
       call obj_init_abinit_xc_names_()
       if (xcObj%id(1) == 1000) then
          !No Hartree and No functional
          write(name,"(A)") ABINIT_XC_NAMES(ABINIT_N_NAMES-1)
       else
          write(name, "(A)") ABINIT_XC_NAMES(min(xcObj%id(1),ABINIT_N_NAMES))
       end if
    end if
  end subroutine obj_get_name_


  !>  Dump XC info on screen.
  subroutine obj_dump_(xcObj)
    implicit none
    type(xc_info), intent(in) :: xcObj

    integer :: i, ii
    type(xc_f90_pointer_t) :: str1
    character(len=500) :: message

    ! Dump functional information
    call obj_get_name_(xcObj, message)
    !write(*,"(1x,a19,a65)") "Exchange-corr. ref.", "("//trim(message)//")"
    call yaml_map('Exchange-Correlation reference','"'//trim(message)//'"')
    if (xcObj%kind == XC_ABINIT) then
       call yaml_map('XC functional implementation','ABINIT')
       !write(*,"(1x,A84)") "XC functional provided by ABINIT."
    else
       call yaml_map('XC functional implementation','libXC')
       call yaml_sequence_open('Reference Papers')
       call yaml_sequence('"Comput. Phys. Commun. 183, 2272 (2012)"')
       do i = 1, 2
          if (xcObj%family(i) == 0) cycle
          
          ii = 0
          if (xcObj%id(i) > 0) then
             call xc_f90_info_refs(xcObj%funcs(i)%info,ii,str1,message)
             do while (ii >= 0)
                !write(*,"(1x,a1,1x,a82)") "|", trim(message)
                call yaml_sequence('"'//trim(message)//'"')
                call xc_f90_info_refs(xcObj%funcs(i)%info,ii,str1,message)
             end do
          end if
       end do
       call yaml_sequence_close()
    end if
  end subroutine obj_dump_

  !>  Get a name, corresponding to the given XC id.
  subroutine xc_get_name(name,ixc,kind)
    implicit none

    character(len=500), intent(out) :: name
    integer, intent(in) :: kind
    integer, intent(in) :: ixc

    type(xc_info) :: xcObj

    call obj_init_(xcObj, ixc, kind, 1)
    call obj_get_name_(xcObj, name)
    call obj_free_(xcObj)
  end subroutine xc_get_name

  !> Get the XC id from the name
  subroutine xc_get_id_from_name(id, name)
    implicit none
    
    integer, intent(out) :: id
    character(len = *), intent(in) :: name

    integer :: i

    i = index(name, "+")
    if (i > 0) then
       id = -xc_f90_functional_get_number(name(1:i - 1)) - 1000 * xc_f90_functional_get_number(name(i+1:len(name)))
    else
       id = -xc_f90_functional_get_number(name)
    end if
  end subroutine xc_get_id_from_name

  !>  Dump XC info on screen.
  subroutine xc_dump(ixc,kind,nspden)
    implicit none
    integer, intent(in) :: kind
    integer, intent(in) :: ixc
    integer, intent(in) :: nspden

    type(xc_info) :: xcObj

    call obj_init_(xcObj, ixc, kind, nspden)
    call obj_dump_(xcObj)
    call obj_free_(xcObj)
  end subroutine xc_dump

  !>  Initialize the desired XC functional, from LibXC.
  subroutine xc_init(xcObj,ixc,kind,nspden,default_alpha)
    implicit none

    !Arguments
    !scalars
    type(xc_info), intent(out) :: xcObj
    integer, intent(in) :: nspden
    integer, intent(in) :: kind
    integer, intent(in) :: ixc
    real(kind=8), intent(in), optional :: default_alpha
    !local variables
!!$    integer :: ixc_prev

    !check if we are trying to initialize the libXC to a different functional
!!$    if (libxc_init) then
!!$       ixc_prev=xc%id(1)+1000*xc%id(2)
!!$       if (f_err_raise(((xc%kind== XC_LIBXC .or. kind == XC_MIXED) .and. ixc/=-ixc_prev),&
!!$            'LibXC library has been already initialized with ixc='//trim(yaml_toa(ixc_prev))//&
!!$            ', finalize it first to reinitialize it with ixc='//trim(yaml_toa(ixc)),&
!!$            err_name='BIGDFT_RUNTIME_ERROR')) return
!!$    end if
!!$    libxc_init=.true.
    if(present(default_alpha)) then
      call obj_init_(xcObj, ixc, kind, nspden, default_alpha)
    else
      call obj_init_(xcObj, ixc, kind, nspden)
    end if
  end subroutine xc_init

  !> End usage of LibXC functional. Call LibXC end function,
  !! and deallocate module contents.
  subroutine xc_end(xcObj)
    implicit none
    type(xc_info), intent(inout) :: xcObj

    !here no exception is needed if finalization has already been done
    !if (libxc_init) then
    call obj_free_(xcObj)
    !   libxc_init=.false.
    !end if

  end subroutine xc_end

  !> Test function to identify whether the presently used functional
  !! is a GGA or not
  function xc_isgga(xcObj)
    implicit none
    type(xc_info), intent(in) :: xcObj

    logical :: xc_isgga

    if (any(xcObj%family == XC_FAMILY_GGA) .or. &
         & any(xcObj%family == XC_FAMILY_HYB_GGA)) then
       xc_isgga = .true.
    else
       xc_isgga = .false.
    end if
  end function xc_isgga

  !> Calculate the exchange-correlation factor (percentage) to add in the functional
  real(kind=8) function xc_exctXfac(xcObj)
    implicit none
    integer :: i
    real(kind=8) :: alpha,omega
    type(xc_info), intent(in) :: xcObj

    xc_exctXfac = 0.d0
    if (any(xcObj%family == XC_FAMILY_HYB_GGA)) then
      ! Was the factor overriden in the input file ?
      if(xcObj%default_alpha .ne. -1.0d0) then
         xc_exctXfac=xcObj%default_alpha
         !write(*,*) 'Alpha HF was overriden : ', xc_exctXfac
      else 
        do i=1,2
          !factors for the exact exchange contribution of different hybrid functionals
          if (xcObj%id(i) == 0) cycle
          select case (abs(xcObj%id(i)))
            case (XC_HYB_GGA_XC_HSE03, XC_HYB_GGA_XC_HSE06, XC_HYB_GGA_XC_HJS_PBE,&
                  XC_HYB_GGA_XC_HJS_PBE_SOL, XC_HYB_GGA_XC_HJS_B88, &
                  XC_HYB_GGA_XC_HJS_B97X, XC_HYB_GGA_XC_CAM_B3LYP, &
                  XC_HYB_GGA_XC_TUNED_CAM_B3LYP, XC_HYB_GGA_XC_CAMY_BLYP) !! SCREENED HYB
              call xc_f90_hyb_cam_coef(xcObj%funcs(i)%conf, omega, alpha, xc_exctXfac)
            case (XC_HARTREE_FOCK) !! PURE HF
              xc_exctXfac = 1.d0
            case default !! HYB
              call xc_f90_hyb_exx_coef(xcObj%funcs(i)%conf,xc_exctXfac)
              !!!write(6,*)"id,fact",xcObj%id(i),xc_exctXfac
           end select
         end do
      end if
    end if
  end function xc_exctXfac

  subroutine xc_init_rho(xcObj, n, rho, nproc)
    implicit none
    ! Arguments
    type(xc_info), intent(in) :: xcObj
    integer :: n,nproc
    real(dp) :: rho(n)

    if (xcObj%kind == XC_ABINIT) then
       call tenminustwenty(n,rho,nproc)
    else
       call f_zero(rho)
    end if
  end subroutine xc_init_rho

  subroutine xc_clean_rho(xcObj, n, rho, nproc)
    implicit none
    ! Arguments
    type(xc_info), intent(in) :: xcObj
    integer :: n,nproc
    real(dp) :: rho(n)

    integer :: i

    if (xcObj%kind == XC_ABINIT) then
       do i = 1, n, 1
          if (rho(i) < 1d-20) then
             rho(i) = 1d-20 / real(nproc,kind=dp)
          end if
       end do
    else
       do i = 1, n, 1
          if (rho(i) < 0.) then
             rho(i) = 0.0_dp
          end if
       end do
    end if
  end subroutine xc_clean_rho

  !> Return XC potential and energy, from input density (even gradient etc...)
  subroutine xc_getvxc(xcObj, npts,exc,nspden,rho,vxc,grho2,vxcgr,dvxci)
    implicit none

    !Arguments ------------------------------------
    type(xc_info), intent(in) :: xcObj
    integer, intent(in) :: npts,nspden
    real(dp),intent(out), dimension(npts) :: exc
    real(dp),intent(in), dimension(npts,nspden)  :: rho
    real(dp),intent(out), dimension(npts,nspden) :: vxc
    real(dp),intent(in) :: grho2(*)
    real(dp),intent(out) :: vxcgr(*)
    real(dp),intent(out), optional :: dvxci(npts,nspden + 1)

    !Local variables-------------------------------
    integer, parameter :: n_blocks = 128
    integer :: i, j, ipts, ixc, ndvxc, ngr2, nd2vxc, nvxcdgr, order
    integer :: ipte, nb
    integer :: ndvxc_tmp,ngr2_tmp,nd2vxc_tmp,nspden_tmp,nvxcdgr_tmp
    real(dp), dimension(n_blocks) :: exctmp(n_blocks)
    real(dp), dimension(nspden, n_blocks) :: rhotmp, vxctmp
    real(dp), allocatable :: rho_tmp(:,:), vxc_tmp(:,:), dvxci_tmp(:,:)
    real(dp), dimension(2*min(nspden,2)-1, n_blocks) :: sigma, vsigma
    real(dp), dimension(3, n_blocks) :: v2rho2
    real(dp), dimension(6, n_blocks) :: v2rhosigma, v2sigma2
    real(dp), allocatable :: rho_(:,:), exc_(:), vxc_(:,:)
    !n(c) character(len=*), parameter :: subname='xc_getvxc'

    call f_routine(id='xc_getvxc')

    if (xcObj%kind == XC_ABINIT) then
       ! ABINIT case, call drivexc
       ixc = xcObj%id(1)
       !Allocations of the exchange-correlation terms, depending on the ixc value
       order = 1
       if (present(dvxci) .and. nspden == 1) order = -2
       if (present(dvxci) .and. nspden == 2) order = +2

       call abi_size_dvxc(ixc,ndvxc,ngr2,nd2vxc,nspden,nvxcdgr,order)

       !let us apply ABINIT routines
       select case (xcObj%family(1))

       case (XC_FAMILY_LDA)
          if (order**2 <=1 .or. ixc >= 31 .and. ixc <= 34) then
             call abi_drivexc(exc,ixc,npts,nspden,order,rho,vxc,&
                  ndvxc,ngr2,nd2vxc,nvxcdgr)

          !Correction of the difference bw the nspin==1 and unpolarized nspin==2 case
          !For nspin==1, ABINIT routines do not respect f_xc = dvxc/drho (but they do for nspin==2).
          !We then force the code to do as if nspin==2.
          else if (nspden==1) then
             ! ABINIT case, call drivexc
             ixc = xcObj%id(1)
             !Allocations of the exchange-correlation terms, depending on the ixc value
             order = 1
             if ( present(dvxci) ) order = -2
             
             !temporary nspden==2, to force the code to do as if nspin==2
             nspden_tmp=2
             !define the size of the temporary variables
             call abi_size_dvxc(ixc,ndvxc_tmp,ngr2_tmp,nd2vxc_tmp,nspden_tmp,nvxcdgr_tmp,order)
             !allocate the temporary variables
             rho_tmp=f_malloc([npts,nspden_tmp],id='rho_tmp')
             vxc_tmp=f_malloc([npts,nspden_tmp],id='vxc_tmp')
             dvxci_tmp=f_malloc([npts,ndvxc_tmp],id='dvxci_tmp')

             !fill density
             call f_memcpy(n=npts,dest=rho_tmp(1,1),src=rho(1,1))
             call f_memcpy(n=npts,dest=rho_tmp(1,2),src=rho(1,1))
             !!fill potential
             !call f_memcpy(n=npts,dest=vxc_tmp(1,1),src=vxc(1,1))
             !call f_memcpy(n=npts,dest=vxc_tmp(1,2),src=vxc(1,1))

             !!$$the density of nspin==1 is too big by a factor 2 (contains the up and down contribution).
             !!$$rho_tmp=0.5*rho_tmp

             !fill vxc_tmp and dvxci_tmp
             call abi_drivexc(exc,ixc,npts,nspden_tmp,order,rho_tmp,vxc_tmp,&
                  ndvxc_tmp,ngr2_tmp,nd2vxc_tmp,nvxcdgr_tmp, dvxc=dvxci_tmp)

             !fill the true vxc
             call axpy(npts,1.d0,vxc_tmp(1,2),1,vxc_tmp(1,1),1) 
             call f_memcpy(n=npts,src=vxc_tmp(1,1),dest=vxc(1,1))
             !fill the true dvxci
             call f_memcpy(n=2*npts,src=dvxci_tmp(1,1),dest=dvxci(1,1))

             !free the temporary variables
             call f_free(rho_tmp)
             call f_free(vxc_tmp)
             call f_free(dvxci_tmp)

          else
             call abi_drivexc(exc,ixc,npts,nspden,order,rho,vxc,&
                  ndvxc,ngr2,nd2vxc,nvxcdgr,&
                  dvxc=dvxci)
          end if


       case (XC_FAMILY_GGA)
          !case with gradient, no big order
          if (ixc /= 13) then             
             call abi_drivexc(exc,ixc,npts,nspden,order,rho,&
                  vxc,ndvxc,ngr2,nd2vxc,nvxcdgr,&
                  grho2_updn=grho2,vxcgr=vxcgr) 
          else
             call abi_drivexc(exc,ixc,npts,nspden,order,rho,&
                  vxc,ndvxc,ngr2,nd2vxc,nvxcdgr,&
                  grho2_updn=grho2) 
          end if
       end select
    else if (xcObj%kind == XC_MIXED) then
       ! LibXC case with ABINIT rho distribution.
       ! Inititalize all relevant arrays to zero
       call f_zero(vxc) !=real(0.0,dp)
       call f_zero(exc) !=real(0.0,dp)
       if (xc_isgga(xcObj)) call f_zero(npts * 3, vxcgr(1))
       if (present(dvxci)) dvxci=real(0,dp)
       !Loop over points
       !$omp parallel do default(shared) &
       !!$omp shared(npts,rho,grho2,nspden) &
       !!$omp shared(xc,vxc,exc,vxcgr,dvxci) &
       !$omp private(ipts,ipte,nb,sigma,j,rhotmp,v2rho2,v2rhosigma,v2sigma2) &
       !$omp private(vsigma,i,exctmp,vxctmp)
       do ipts = 1, npts, n_blocks
          ipte = min(ipts + n_blocks - 1, npts)
          nb = ipte - ipts + 1

          ! Convert the quantities provided by ABINIT to the ones needed by libxc
          if (nspden == 1) then
             ! ABINIT passes rho_up in the spin-unpolarized case, while the libxc
             ! expects the total density
             rhotmp(1, 1:nb) = real(2,dp) * rho(ipts:ipte,1)
          else
             rhotmp(1, 1:nb) = rho(ipts:ipte,1)
             rhotmp(2, 1:nb) = rho(ipts:ipte,2)
          end if
          if (xc_isgga(xcObj)) then
             sigma=real(0,dp)
             if (nspden==1) then
                ! ABINIT passes |rho_up|^2 while Libxc needs |rho_tot|^2
                sigma(1, 1:nb) = real(4,dp) * grho2(ipts:ipte)
             else
                ! ABINIT passes |rho_up|^2, |rho_dn|^2, and |rho_tot|^2
                ! while Libxc needs |rho_up|^2, rho_up.rho_dn, and |rho_dn|^2
                do i = 1, nb
                   sigma(1, i) = grho2(ipts + i - 1)
                   sigma(2, i) = (grho2(ipts + i - 1 + 2 * npts) - &
                        & grho2(ipts + i - 1) - grho2(ipts + i - 1 + npts))/real(2,dp)
                   sigma(3, i) = grho2(ipts + i - 1 + npts)
                end do
             end if
          end if
          !Loop over functionals
          do i = 1,2
             if (xcObj%id(i) == 0) cycle

             !Get the potential (and possibly the energy)
             if (iand(xc_f90_info_flags(xcObj%funcs(i)%info), XC_FLAGS_HAVE_EXC) .ne. 0) then
                select case (xcObj%family(i))
                case (XC_FAMILY_LDA)
                   call xc_f90_lda_exc_vxc(xcObj%funcs(i)%conf,nb,rhotmp(1,1),exctmp(1),vxctmp(1,1))
                case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
                   call xc_f90_gga_exc_vxc(xcObj%funcs(i)%conf,nb,rhotmp(1,1),sigma(1,1),exctmp(1),vxctmp(1,1),vsigma(1,1))
                end select

             else
                exctmp=real(0,dp)
                select case (xcObj%family(i))
                case (XC_FAMILY_LDA)
                   call xc_f90_lda_vxc(xcObj%funcs(i)%conf,nb,rhotmp(1,1),vxctmp(1,1))
                case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
                   call xc_f90_gga_vxc(xcObj%funcs(i)%conf,nb,rhotmp(1,1),sigma(1,1),vxctmp(1,1),vsigma(1,1))
                end select
             end if

             if (present(dvxci)) then
                v2rho2     = real(0, dp)
                v2rhosigma = real(0, dp)
                v2sigma2   = real(0, dp)
                if (iand(xc_f90_info_flags(xcObj%funcs(i)%info), XC_FLAGS_HAVE_FXC) .ne. 0) then
                   select case (xcObj%family(i))
                   case (XC_FAMILY_LDA)
                      call xc_f90_lda_fxc(xcObj%funcs(i)%conf,nb,rhotmp(1,1),v2rho2(1,1))
                   case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
                      call xc_f90_gga_fxc(xcObj%funcs(i)%conf,nb,rhotmp(1,1),sigma(1,1), &
                           & v2rho2(1,1),v2rhosigma(1,1), v2sigma2(1,1))
                   end select
                end if
             end if

             exc(ipts:ipte) = exc(ipts:ipte) + exctmp(1:nb)
             do j = 1, nspden
                vxc(ipts:ipte,j) = vxc(ipts:ipte,j) + vxctmp(j, 1:nb)
             end do
             
             if (xc_isgga(xcObj)) then
                !Convert the quantities returned by Libxc to the ones needed by ABINIT
                if (nspden == 1) then
                      vxcgr(ipts + 2 * npts:ipte + 2 * npts) = &
                           & vxcgr(ipts + 2 * npts:ipte + 2 * npts) + &
                           & vsigma(1, 1:nb)*real(2,dp)
                else
                   !here a FPE has been found in XC test
                   vxcgr(ipts:ipte) = &
                        & vxcgr(ipts:ipte) + real(2,dp)*vsigma(1, 1:nb) - vsigma(2, 1:nb)
                   vxcgr(ipts + npts:ipte + npts) = &
                        & vxcgr(ipts + npts:ipte + npts) + &
                        & real(2,dp)*vsigma(3, 1:nb) - vsigma(2, 1:nb)
                   vxcgr(ipts + 2 * npts:ipte + 2 * npts) = &
                        & vxcgr(ipts + 2 * npts:ipte + 2 * npts) + vsigma(2, 1:nb)
                end if
             end if

             if (present(dvxci)) then
                !Convert the quantities returned by Libxc to the ones needed by ABINIT
                if (nspden == 1) then
                   dvxci(ipts:ipte,1) = dvxci(ipts:ipte,1) + v2rho2(1,1:nb)*real(2,dp)
                else
                   ! TODO
!!$                dvxci(ipts,1) = dvxci(ipts,1) + real(2,dp)*vsigma(1) - vsigma(2)
!!$                dvxci(ipts,2) = dvxci(ipts,2) + real(2,dp)*vsigma(3) - vsigma(2)
!!$                dvxci(ipts,3) = dvxci(ipts,3) + vsigma(2)
                end if
             end if

          end do
       end do
       !$omp end parallel do

    else if (xcObj%kind == XC_LIBXC) then
       ! Pure LibXC case.
       ! WARNING: LDA implementation only, first derivative, no fxc

       allocate(rho_(nspden, npts))
       allocate(exc_(npts))
       allocate(vxc_(nspden, npts))

       ! Convert the quantities provided by BigDFT to the ones needed by libxc
       if (nspden == 1) then
          rho_ = real(2,dp) * reshape(rho, (/ nspden, npts /))
       else
          call vcopy(npts * nspden, rho(1,1), 1, rho_(1,1), 1)
       end if

       exc(1) = UNINITIALIZED(1.d0)
       !Loop over functionals
       do i = 1,2
          if (xcObj%id(i) == 0) cycle

          !Get the potential (and possibly the energy)
          if (iand(xc_f90_info_flags(xcObj%funcs(i)%info), XC_FLAGS_HAVE_EXC) .ne. 0) then
             select case (xcObj%family(i))
             case (XC_FAMILY_LDA)
                call xc_f90_lda_exc_vxc(xcObj%funcs(i)%conf,npts,rho_(1,1),exc_(1),vxc_(1,1))
             end select
          else
             exc_=real(0,dp)
             select case (xcObj%family(i))
             case (XC_FAMILY_LDA)
                call xc_f90_lda_vxc(xcObj%funcs(i)%conf,npts,rho_(1,1),vxc_(1,1))
             end select
          end if

          if (exc(1) == UNINITIALIZED(1.d0)) then
             call vcopy(npts, exc_(1), 1, exc(1), 1)
             call vcopy(npts * nspden, vxc_(1,1), 1, vxc(1,1), 1)
          else
             call daxpy(npts, 1.d0, exc_, 1, exc, 1)
             call daxpy(npts * nspden, 1.d0, vxc_, 1, vxc, 1)
          end if
       end do
       deallocate(rho_)
       deallocate(exc_)
       deallocate(vxc_)
    else
       call f_err_throw("XC module not initialised, the object kind is '"//xcObj%kind//&
            "', which is neither XC_LIBXC ("+XC_LIBXC//") nor XC_ABINIT ("+XC_ABINIT//')',&
            err_name='BIGDFT_RUNTIME_ERROR')
    end if

    call f_release_routine()

  end subroutine xc_getvxc


  !> Initialize the names of the xc functionals used by ABINIT
  subroutine obj_init_abinit_xc_names_()
    if (abinit_init) return

    write(ABINIT_XC_NAMES( 0), "(A)") "XC: NO Semilocal XC (Hartree only)"
    write(ABINIT_XC_NAMES( 1), "(A)") "XC: Teter 93"
    write(ABINIT_XC_NAMES( 2), "(A)") "XC: Slater exchange, Perdew & Zunger"
    write(ABINIT_XC_NAMES( 3), "(A)") "XC: Teter 91"
    write(ABINIT_XC_NAMES( 4), "(A)") "XC: Slater exchange, Wigner"
    write(ABINIT_XC_NAMES( 5), "(A)") "XC: Slater exchange, Hedin & Lundqvist"
    write(ABINIT_XC_NAMES( 6), "(A)") "-C: Slater's Xalpha"
    write(ABINIT_XC_NAMES( 7), "(A)") "XC: Slater exchange, Perdew & Wang"
    write(ABINIT_XC_NAMES( 8), "(A)") "X-: Slater exchange"
    write(ABINIT_XC_NAMES( 9), "(A)") "XC: Slater exchange, Random Phase Approximation (RPA)"
    write(ABINIT_XC_NAMES(11), "(A)") "XC: Perdew, Burke & Ernzerhof"
    write(ABINIT_XC_NAMES(12), "(A)") "X-: Perdew, Burke & Ernzerhof"
    write(ABINIT_XC_NAMES(13), "(A)") "XC: van Leeuwen & Baerends"
    write(ABINIT_XC_NAMES(14), "(A)") "XC: Revised PBE from Zhang & Yang, Perdew, Burke & Ernzerhof"
    write(ABINIT_XC_NAMES(15), "(A)") "XC: Hammer, Hansen, and Nørskov, Perdew, Burke & Ernzerhof"
    write(ABINIT_XC_NAMES(16), "(A)") "XC: HCTH/93"
    write(ABINIT_XC_NAMES(17), "(A)") "XC: HCTH/120"
    write(ABINIT_XC_NAMES(23), "(A)") "X-: Wu & Cohen"
    write(ABINIT_XC_NAMES(26), "(A)") "XC: HCTH/147"
    write(ABINIT_XC_NAMES(27), "(A)") "XC: HCTH/407"
    write(ABINIT_XC_NAMES(28), "(A)") "No Hartree and XC terms"
    !Hartree-Fock (always the last - ixc=100)
    write(ABINIT_XC_NAMES(29), "(A)") "Hartree-Fock Exchange only"

    abinit_init = .true.
  end subroutine obj_init_abinit_xc_names_


  !> Calculate the XC terms from the given density in a distributed way.
  !! it assign also the proper part of the density to the zf array 
  !! which will be used for the core of the FFT procedure.
  !! Following the values of ixc and of sumpion, the array pot_ion is either summed or assigned
  !! to the XC potential, or even ignored.
  !!
  !! @warning
  !!    The dimensions of pot_ion must be compatible with geocode, datacode and ixc.
  !!    Since the arguments of these routines are indicated with the *,
  !!    it is IMPERATIVE to refer to PSolver routine for the correct allocation sizes.
  subroutine xc_energy_new(geocode,m1,m3,nxc,nwb,nxt,nwbl,nwbr,&
       nxcl,nxcr,xc,hx,hy,hz,rho,gradient,vxci,exc,vxc,order,ndvxc,dvxci,nspden,wbstr)

    implicit none
    !Arguments
    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
    integer, intent(in) :: m1,m3     !< Global dimensions in the three directions.
    integer, intent(in) :: nxc       !< Value of the effective distributed dimension in the third direction
    integer, intent(in) :: nwb       !< Enlarged dimension for calculating the WB correction
    integer, intent(in) :: nxt       !< Enlarged dimension for calculating the GGA case 
    !! (further enlarged for compatibility with WB correction if it is the case)
    integer, intent(in) :: nwbl,nwbr !< nwb=nxc+nxcl+nxcr-2, nwb+nwbl+nwbr=nxt.
    integer, intent(in) :: nxcl,nxcr !< Shifts in the three directions to be compatible with the relation
    !> eXchange-Correlation code. Indicates the XC functional to be used 
    !!   for calculating XC energies and potential. 
    !!   ixc=0 indicates that no XC terms are computed. 
    !!   The XC functional codes follow the ABINIT convention or if negative the libXC one.
    type(xc_info), intent(in) :: xc
    integer, intent(in) :: order,ndvxc,nspden
    real(gp), intent(in) :: hx,hy,hz                            !< Grid spacings. 
    real(dp), dimension(*), intent(in) :: gradient              !< of size 1 if not needed
    real(dp), dimension(m1,m3,nxt,nspden), intent(inout) :: rho !< Density in the distributed format, also in spin-polarised
    real(dp), dimension(m1,m3,nwb,nspden), intent(out) :: vxci
    real(dp), dimension(m1,m3,nwb,ndvxc), intent(out) :: dvxci
    real(dp), intent(out) :: exc,vxc                            !< XC energy and integral of @f$\rho V_{xc}@f$ respectively
    real(dp), dimension(6), intent(inout) :: wbstr

    !Local variables----------------
    character(len=*), parameter :: subname='xc_energy'
    real(dp), dimension(:,:,:), allocatable :: exci
    real(dp), dimension(:,:,:,:), allocatable :: dvxcdgr
    !real(dp), dimension(:,:,:,:,:), allocatable :: gradient
    real(dp) :: elocal,vlocal,rhov,sfactor
    integer :: npts,offset,ispden
    integer :: i1,i2,i3,j1,j2,j3,jp2,jppp2
    logical :: use_gradient

    call f_routine(id='xc_energy_new')

    !check for the dimensions
    if (nwb/=nxcl+nxc+nxcr-2 .or. nxt/=nwbr+nwb+nwbl) then
       call f_err_throw('The XC dimensions are not correct; see '//&
            'nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr'//yaml_toa([nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr]),&
            err_name='BIGDFT_RUNTIME_ERROR')
       return
    end if

    !starting point of the density array for the GGA cases in parallel
    offset=nwbl+1
    !divide by two the density to applicate it in the ABINIT xc routines
    use_gradient = xc_isgga(xc)

    if (use_gradient) then
       dvxcdgr = f_malloc((/ m1, m3, nwb, 3 /),id='dvxcdgr')
    else
       dvxcdgr = f_malloc((/ 1, 1, 1, 1 /),id='dvxcdgr')
    end if

    !Allocations
    exci = f_malloc((/ m1, m3, nwb /),id='exci')

    npts=m1*m3*nwb

    !do a separate calculation of the grid to allow for OMP parallelisation
    ! Do the calculation.
    if (abs(order) == 1) then
       call xc_getvxc(xc, npts,exci,nspden,rho(1,1,offset,1),vxci,gradient,dvxcdgr)
    else if (abs(order) == 2) then
       call xc_getvxc(xc, npts,exci,nspden,rho(1,1,offset,1),vxci,gradient,dvxcdgr,dvxci)
    end if
    wbstr(:)=0._dp
    if (use_gradient) then
       ! Do not calculate the White-Bird term in the Leeuwen Baerends XC case
       if (xc%ixc /= 13 .and. xc%ixc /= -160) then
          call vxcpostprocessing(geocode,m1,m3,nwb,nxc,nxcl,nxcr,nspden,3,gradient,&
               real(hx,dp),real(hy,dp),real(hz,dp),dvxcdgr,vxci,wbstr)
       end if

       !restore the density array in the good position if it was shifted for the parallel GGA
       !operation not necessarily needed, but related to the fact that the array has three
       !indices which make it difficult to treat
       !one should convert the operations with one indices arrays
       if (nspden==2 .and. nxt /= nwb) then
          j3=nwb+1
          do i3=nwb-nwbr,1,-1
             j3=j3-1
             do i2=1,m3
                do i1=1,m1
                   rho(i1,i2,nwbl+j3,2)=rho(i1,i2,i3,2)
                end do
             end do
          end do
          do i3=nxt,nwb+nwbl+1,-1 !we have nwbr points
             j3=j3-1
             do i2=1,m3
                do i1=1,m1
                   rho(i1,i2,nwbl+j3,2)=rho(i1,i2,i3,1)
                end do
             end do
          end do
       end if
    end if

    call f_free(dvxcdgr)

    !this part should be put out from this routine due to the Global distribution code
    exc=0.0_dp
    vxc=0.0_dp
    sfactor=1.0_dp
    if(nspden==1) sfactor=2.0_dp

    !compact the rho array into the total charge density
    !try to use dot and vcopy routines, more general
    ! e.g. exc=dot(m1*m3*nxc,exci(1,1,nxcl),1,rho(1,1,offset+nxcl-1,ispden),1)

    ispden=1
    do jp2=1,nxc
       j2=offset+jp2+nxcl-2
       jppp2=jp2+nxcl-1
       do j3=1,m3
          do j1=1,m1
             rhov=rho(j1,j3,j2,ispden)
             elocal=exci(j1,j3,jppp2)
             vlocal=vxci(j1,j3,jppp2,ispden)
             exc=exc+elocal*rhov
             vxc=vxc+vlocal*rhov
             rho(j1,j3,jp2,1)=sfactor*rhov!restore the original normalization
             !potxc(j1,j3,jp2,ispden)=real(vlocal,wp)
          end do
       end do
    end do
    !spin-polarised case
    if (nspden==2) then
       ispden=2
       do jp2=1,nxc
          j2=offset+jp2+nxcl-2
          jppp2=jp2+nxcl-1
          do j3=1,m3
             do j1=1,m1
                rhov=rho(j1,j3,j2,ispden)
                elocal=exci(j1,j3,jppp2)
                vlocal=vxci(j1,j3,jppp2,ispden)
                exc=exc+elocal*rhov
                vxc=vxc+vlocal*rhov
                rho(j1,j3,jp2,1)=rho(j1,j3,jp2,1)+sfactor*rhov
                !potxc(j1,j3,jp2,ispden)=real(vlocal,dp)
             end do
          end do
       end do
    end if

    !the two factor is due to the 
    !need of using the density of states in abinit routines
    exc=sfactor*real(hx*hy*hz,dp)*exc
    vxc=sfactor*real(hx*hy*hz,dp)*vxc

    !De-allocations
    call f_free(exci)

    call f_release_routine()

  END SUBROUTINE xc_energy_new

end module module_xc
