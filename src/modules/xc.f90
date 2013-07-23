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
  use interfaces_56_xc
  use yaml_output

  implicit none

  integer, public, parameter :: XC_ABINIT = 1  !< xc functionals from ABINIT
  integer, public, parameter :: XC_LIBXC  = 2  !< xc from libxc
  integer, public, parameter :: XC_MIXED  = 3  !> xc mixing the origin of the xc functionals

  !> Structure containing the information to call the routines for the calculation of the xc functionals
  type libxc_functional
     private
     type(xc_f90_pointer_t) :: conf !< the pointer used to call the library
     type(xc_f90_pointer_t) :: info !< Information about the functional
  end type libxc_functional

  !> Structure containing the information about the xc functionals
  type xc_info
     private
     integer :: kind       !< ABINIT or LibXC
     integer :: family(2)  !< LDA, GGA, etc.
     integer :: id(2)      !< Identifier

     type(libxc_functional) :: funcs(2)
  end type xc_info

  type(xc_info) :: xc      !< Global structure about the used xc functionals

  logical :: abinit_init = .false.            !< .True. if already ABINIT_XC_NAMES intialized
  character(len=500) :: ABINIT_XC_NAMES(0:28) !< Names of the xc functionals used by ABINIT

  private
  public :: xc_init, &
       &    xc_dump, &
       &    xc_init_rho, &
       &    xc_clean_rho, &
       &    xc_getvxc, &
       &    xc_isgga, &
       &    xc_exctXfac, &
       &    xc_end, &
       &    xc_get_name

contains

  subroutine obj_init_(xcObj, ixc, kind, nspden)
    implicit none

    !Arguments
    !scalars
    type(xc_info), intent(out) :: xcObj
    integer, intent(in) :: nspden
    integer, intent(in) :: kind
    integer, intent(in) :: ixc

    !Local variables
    !scalars
    integer :: i, ierr

    ! **

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
       else if (xcObj%id(1) == 0 .or. xcObj%id(1)==100) then
          xcObj%family(1) = 0
       else
          write(*,*) "Error: unsupported functional, change ixc."
          call MPI_ABORT(bigdft_mpi%mpi_comm,0,ierr)
       end if
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

  subroutine obj_get_name_(xcObj, name)
    implicit none

    character(len=500), intent(out) :: name
    type(xc_info), intent(in) :: xcObj
    
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
       write(name, "(A)") ABINIT_XC_NAMES(min(xcObj%id(1),28))
    end if
  end subroutine obj_get_name_

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

  !>  Initialize the desired XC functional, from LibXC.
  subroutine xc_init(ixc,kind,nspden)
    implicit none

    !Arguments
    !scalars
    integer, intent(in) :: nspden
    integer, intent(in) :: kind
    integer, intent(in) :: ixc

    call obj_init_(xc, ixc, kind, nspden)
  end subroutine xc_init

  !>  Dump XC info on screen.
  subroutine xc_dump()
    implicit none

    integer :: i, ii
    type(xc_f90_pointer_t) :: str
    character(len=500) :: message

    ! Dump functional information
    call obj_get_name_(xc, message)
    !write(*,"(1x,a19,a65)") "Exchange-corr. ref.", "("//trim(message)//")"
    call yaml_map('Exchange-Correlation reference','"'//trim(message)//'"')
    if (xc%kind == XC_ABINIT) then
       call yaml_map('XC functional implementation','ABINIT')
       !write(*,"(1x,A84)") "XC functional provided by ABINIT."
    else
       call yaml_map('XC functional implementation','libXC')
       call yaml_open_sequence('Reference Papers')
       call yaml_sequence('"Comput. Phys. Commun. 183, 2272 (2012)"')
       do i = 1, 2
          if (xc%family(i) == 0) cycle
          
          ii = 0
          if (xc%id(i) > 0) then
             call xc_f90_info_refs(xc%funcs(i)%info,ii,str,message)
             do while (ii >= 0)
                !write(*,"(1x,a1,1x,a82)") "|", trim(message)
                call yaml_sequence('"'//trim(message)//'"')
                call xc_f90_info_refs(xc%funcs(i)%info,ii,str,message)
             end do
          end if
       end do
       call yaml_close_sequence()
    end if
  end subroutine xc_dump

  !> End usage of LibXC functional. Call LibXC end function,
  !! and deallocate module contents.
  subroutine xc_end()
    implicit none

    call obj_free_(xc)
  end subroutine xc_end

  !> Test function to identify whether the presently used functional
  !! is a GGA or not
  function xc_isgga()
    implicit none

    logical :: xc_isgga

    if (any(xc%family == XC_FAMILY_GGA) .or. &
         & any(xc%family == XC_FAMILY_HYB_GGA)) then
       xc_isgga = .true.
    else
       xc_isgga = .false.
    end if
  end function xc_isgga

  real(kind=8) function xc_exctXfac()
    implicit none

    xc_exctXfac = 0.d0
    if (any(xc%family == XC_FAMILY_HYB_GGA)) then
       !factors for the exact exchange contribution of different hybrid functionals
       if (any(xc%id == XC_HYB_GGA_XC_PBEH)) then
          xc_exctXfac = 0.25d0 
       end if
    end if

    !hartree-fock value
    if (xc%id(1) == 100 .and. xc%id(2) == 0) then
       xc_exctXfac = 1.d0 
    end if

  end function xc_exctXfac

  subroutine xc_init_rho(n, rho, nproc)
    implicit none
    ! Arguments
    integer :: n,nproc
    real(dp) :: rho(n)

    if (xc%kind == XC_ABINIT) then
       call tenminustwenty(n,rho,nproc)
    else
       call razero(n,rho)
    end if
  end subroutine xc_init_rho

  subroutine xc_clean_rho(n, rho, nproc)
    implicit none
    ! Arguments
    integer :: n,nproc
    real(dp) :: rho(n)

    integer :: i

    if (xc%kind == XC_ABINIT) then
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
  subroutine xc_getvxc(npts,exc,nspden,rho,vxc,grho2,vxcgr,dvxci)
    implicit none

    !Arguments ------------------------------------
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
    real(dp), dimension(n_blocks) :: exctmp(n_blocks)
    real(dp), dimension(nspden, n_blocks) :: rhotmp, vxctmp
    real(dp), dimension(2*min(nspden,2)-1, n_blocks) :: sigma, vsigma
    real(dp), dimension(3, n_blocks) :: v2rho2
    real(dp), dimension(6, n_blocks) :: v2rhosigma, v2sigma2
    real(dp), allocatable :: rho_(:,:), exc_(:), vxc_(:,:)
    !n(c) character(len=*), parameter :: subname='xc_getvxc'

    if (xc%kind == XC_ABINIT) then
       ! ABINIT case, call drivexc
       ixc = xc%id(1)
       !Allocations of the exchange-correlation terms, depending on the ixc value
       order = 1
       if (present(dvxci) .and. nspden == 1) order = -2
       if (present(dvxci) .and. nspden == 2) order = +2

       call size_dvxc(ixc,ndvxc,ngr2,nd2vxc,nspden,nvxcdgr,order)

       !let us apply ABINIT routines
       select case (xc%family(1))
       case (XC_FAMILY_LDA)
          if (order**2 <=1 .or. ixc >= 31 .and. ixc <= 34) then
             call drivexc(exc,ixc,npts,nspden,order,rho,vxc,&
                  ndvxc,ngr2,nd2vxc,nvxcdgr)
          else
             call drivexc(exc,ixc,npts,nspden,order,rho,vxc,&
                  ndvxc,ngr2,nd2vxc,nvxcdgr,&
                  dvxc=dvxci)
          end if
       case (XC_FAMILY_GGA)
          !case with gradient, no big order
          if (ixc /= 13) then             
             call drivexc(exc,ixc,npts,nspden,order,rho,&
                  vxc,ndvxc,ngr2,nd2vxc,nvxcdgr,&
                  grho2_updn=grho2,vxcgr=vxcgr) 
          else
             call drivexc(exc,ixc,npts,nspden,order,rho,&
                  vxc,ndvxc,ngr2,nd2vxc,nvxcdgr,&
                  grho2_updn=grho2) 
          end if
       end select
    else if (xc%kind == XC_MIXED) then
       ! LibXC case with ABINIT rho distribution.
       ! Inititalize all relevant arrays to zero
       vxc=real(0,dp)
       exc=real(0,dp)
       if (xc_isgga()) call to_zero(npts * 3, vxcgr(1))
       if (present(dvxci)) dvxci=real(0,dp)

       !Loop over points
       !$omp parallel do default(private) &
       !$omp & shared(npts,rho,grho2,nspden) &
       !$omp & shared(xc,vxc,exc,vxcgr,dvxci) 
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
          if (xc_isgga()) then
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
             if (xc%id(i) == 0) cycle

             !Get the potential (and possibly the energy)
             if (iand(xc_f90_info_flags(xc%funcs(i)%info), XC_FLAGS_HAVE_EXC) .ne. 0) then
                select case (xc%family(i))
                case (XC_FAMILY_LDA)
                   call xc_f90_lda_exc_vxc(xc%funcs(i)%conf,nb,rhotmp(1,1),exctmp(1),vxctmp(1,1))
                case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
                   call xc_f90_gga_exc_vxc(xc%funcs(i)%conf,nb,rhotmp(1,1),sigma(1,1),exctmp(1),vxctmp(1,1),vsigma(1,1))
                end select

             else
                exctmp=real(0,dp)
                select case (xc%family(i))
                case (XC_FAMILY_LDA)
                   call xc_f90_lda_vxc(xc%funcs(i)%conf,nb,rhotmp(1,1),vxctmp(1,1))
                case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
                   call xc_f90_gga_vxc(xc%funcs(i)%conf,nb,rhotmp(1,1),sigma(1,1),vxctmp(1,1),vsigma(1,1))
                end select
             end if

             if (present(dvxci)) then
                v2rho2     = real(0, dp)
                v2rhosigma = real(0, dp)
                v2sigma2   = real(0, dp)
                if (iand(xc_f90_info_flags(xc%funcs(i)%info), XC_FLAGS_HAVE_FXC) .ne. 0) then
                   select case (xc%family(i))
                   case (XC_FAMILY_LDA)
                      call xc_f90_lda_fxc(xc%funcs(i)%conf,nb,rhotmp(1,1),v2rho2(1,1))
                   case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
                      call xc_f90_gga_fxc(xc%funcs(i)%conf,nb,rhotmp(1,1),sigma(1,1), &
                           & v2rho2(1,1),v2rhosigma(1,1), v2sigma2(1,1))
                   end select
                end if
             end if

             exc(ipts:ipte) = exc(ipts:ipte) + exctmp(1:nb)
             do j = 1, nspden
                vxc(ipts:ipte,j) = vxc(ipts:ipte,j) + vxctmp(j, 1:nb)
             end do
             
             if (xc_isgga()) then
                !Convert the quantities returned by Libxc to the ones needed by ABINIT
                if (nspden == 1) then
                      vxcgr(ipts + 2 * npts:ipte + 2 * npts) = &
                           & vxcgr(ipts + 2 * npts:ipte + 2 * npts) + &
                           & vsigma(1, 1:nb)*real(2,dp)
                else
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

    else if (xc%kind == XC_LIBXC) then
       ! Pure LibXC case.
       ! WARNING: LDA implementation only, first derivative, no fxc

       allocate(rho_(nspden, npts))
       allocate(exc_(npts))
       allocate(vxc_(nspden, npts))

       ! Convert the quantities provided by BigDFT to the ones needed by libxc
       if (nspden == 1) then
          rho_ = real(2,dp) * reshape(rho, (/ nspden, npts /))
       else
          call dcopy(npts * nspden, rho, 1, rho_, 1)
       end if

       exc(1) = UNINITIALIZED(1.d0)
       !Loop over functionals
       do i = 1,2
          if (xc%id(i) == 0) cycle

          !Get the potential (and possibly the energy)
          if (iand(xc_f90_info_flags(xc%funcs(i)%info), XC_FLAGS_HAVE_EXC) .ne. 0) then
             select case (xc%family(i))
             case (XC_FAMILY_LDA)
                call xc_f90_lda_exc_vxc(xc%funcs(i)%conf,npts,rho_(1,1),exc_(1),vxc_(1,1))
             end select
          else
             exc_=real(0,dp)
             select case (xc%family(i))
             case (XC_FAMILY_LDA)
                call xc_f90_lda_vxc(xc%funcs(i)%conf,npts,rho_(1,1),vxc_(1,1))
             end select
          end if

          if (exc(1) == UNINITIALIZED(1.d0)) then
             call dcopy(npts, exc_, 1, exc, 1)
             call dcopy(npts * nspden, vxc_, 1, vxc, 1)
          else
             call daxpy(npts, 1.d0, exc_, 1, exc, 1)
             call daxpy(npts * nspden, 1.d0, vxc_, 1, vxc, 1)
          end if
       end do
       deallocate(rho_)
       deallocate(exc_)
       deallocate(vxc_)
    else
       write(0,*) "ERROR: XC module not initialised."
    end if

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
    !Hatree-fock (always the last - ixc=100)
    write(ABINIT_XC_NAMES(28), "(A)") "Hartree-Fock Exchange only"

    abinit_init = .true.
  end subroutine obj_init_abinit_xc_names_

end module module_xc
