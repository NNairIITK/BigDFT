!> @file
!!  Wrapper around XC library routines (both BAINIT and LibXC).
!! @author
!!    Copyright (C) 2008-2010 ABINIT group (MOliveira)
!!    Copyright (C) 2008-2011 BigDFT group (DC)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

module module_xc

  use module_base
  use xc_f90_types_m
  use libxc_funcs_m
  use xc_f90_lib_m
  use interfaces_56_xc

  implicit none

  integer, public, parameter :: XC_ABINIT = 1
  integer, public, parameter :: XC_LIBXC  = 2
  integer, public, parameter :: XC_MIXED  = 3

  type libxc_functional
     private
     type(xc_f90_pointer_t) :: conf ! the pointer used to call the library
     type(xc_f90_pointer_t) :: info ! information about the functional
  end type libxc_functional

  type xc_info
     private
     integer :: kind       ! ABINIT or LibXC
     integer :: family(2)  ! LDA, GGA, etc.
     integer :: id(2)      ! identifier

     type(libxc_functional) :: funcs(2)
  end type xc_info

  type(xc_info) :: xc

  private
  public :: xc_init, &
       &    xc_dump, &
       &    xc_init_rho, &
       &    xc_clean_rho, &
       &    xc_getvxc, &
       &    xc_isgga, &
       &    xc_exctXfac, &
       &    xc_end

contains

  !>  Initialize the desired XC functional, from LibXC.
  subroutine xc_init(ixc,kind,nspden)

    implicit none

    !Arguments ------------------------------------
    !scalars
    integer, intent(in) :: nspden
    integer, intent(in) :: kind
    integer, intent(in) :: ixc

    !Local variables-------------------------------
    !scalars
    integer :: i, ierr

    ! *************************************************************************

    xc%kind = kind
    if (kind == XC_LIBXC .or. kind == XC_MIXED) then
       xc%id(2) = abs(ixc) / 1000
       xc%id(1) = abs(ixc) - xc%id(2) * 1000
    else if (kind == XC_ABINIT) then
       ! ABINIT functional.
       xc%id(1) = abs(ixc)
       xc%id(2) = 0
    end if

    xc%family = 0
    if (xc%kind == XC_LIBXC .or. xc%kind == XC_MIXED) then
       ! LibXC case.

       do i = 1, 2
          if (xc%id(i) == 0) cycle
          ! Get XC functional family
          xc%family(i) = xc_f90_family_from_id(xc%id(i))
          select case (xc%family(i))
          case (XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
             call xc_f90_func_init(xc%funcs(i)%conf,xc%funcs(i)%info,xc%id(i),nspden)
          case default
             write(*,*) "Error: unsupported functional, change ixc."
             call MPI_ABORT(MPI_COMM_WORLD,0,ierr)
          end select
       end do

    else if (xc%kind == XC_ABINIT) then
       ! ABINIT case

       ! Get XC functional family
       if ((xc%id(1) >= 1 .and. xc%id(1) < 11) .or. xc%id(1) == 24) then
          xc%family(1) = XC_FAMILY_LDA
       else if (xc%id(1) >= 11 .and. xc%id(1) < 18) then
          xc%family(1) = XC_FAMILY_GGA
       else if (xc%id(1) >= 23 .and. xc%id(1) < 28) then
          xc%family(1) = XC_FAMILY_GGA
       else if (xc%id(1) >= 31 .and. xc%id(1) < 35) then
          xc%family(1) = XC_FAMILY_LDA
       else if (xc%id(1) == 0) then
          xc%family(1) = 0
       else
          write(*,*) "Error: unsupported functional, change ixc."
          call MPI_ABORT(MPI_COMM_WORLD,0,ierr)
       end if
    end if
  end subroutine xc_init

  !>  Dump XC info on screen.
  subroutine xc_dump()
    implicit none

    integer :: ii, jj
    type(xc_f90_pointer_t) :: str
    character(len=500) :: message, message2

    ! Dump functional information
    if (xc%kind == XC_ABINIT) then
       write(*,"(1x,A41,1x,A1)") "XC functional provided by ABINIT.", "|"
    else
       if (xc%id(1) > 0) then
          call xc_f90_info_name(xc%funcs(1)%info,message)
       else
          write(message, "(A)") ""
       end if
       if (xc%id(2) > 0) then
          call xc_f90_info_name(xc%funcs(2)%info,message2)
       else
          write(message2, "(A)") ""
       end if
       write(*,"(1x,a41,1x,A1,1x,A40)") trim(message), "|", trim(message2)
       ii = 0
       jj = 0
       if (xc%id(1) > 0) then
          call xc_f90_info_refs(xc%funcs(1)%info,ii,str,message)
       else
          ii = -1
       end if
       if (xc%id(2) > 0) then
          call xc_f90_info_refs(xc%funcs(2)%info,jj,str,message2)
       else
          jj = -1
       end if
       do while (ii >= 0 .or. jj >= 0)
          if (ii >= 0) then
             write(*,"(1x,a41)", advance = "NO") trim(message)
             call xc_f90_info_refs(xc%funcs(1)%info,ii,str,message)
          else
             write(*,"(1x,a41)", advance = "NO") " "
          end if
          write(*, "(1x,A1,1x)", advance = "NO") "|"
          if (jj >= 0) then
             write(*,"(a40)") trim(message2)
             call xc_f90_info_refs(xc%funcs(2)%info,jj,str,message2)
          else
             write(*,"(a40)") " "
          end if
       end do
    end if
  end subroutine xc_dump

  !>  End usage of LibXC functional. Call LibXC end function,
  !!  and deallocate module contents.
  subroutine xc_end()
    implicit none

    if (xc%kind == XC_LIBXC .or. xc%kind == XC_MIXED) then
       if (xc%id(1) > 0) call xc_f90_func_end(xc%funcs(1)%conf)
       if (xc%id(2) > 0) call xc_f90_func_end(xc%funcs(2)%conf)
    end if
  end subroutine xc_end

  !>  Test function to identify whether the presently used functional
  !!  is a GGA or not
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
             rho(i) = 1d-20 / nproc
          end if
       end do
    else
       do i = 1, n, 1
          if (rho(i) < 0.) then
             rho(i) = 0.
          end if
       end do
    end if
  end subroutine xc_clean_rho

  !>  Return XC potential and energy, from input density (even gradient etc...)
  subroutine xc_getvxc(npts,exc,nspden,rho,vxc,grho2,vxcgr,dvxci)
    implicit none

    !Arguments ------------------------------------
    integer, intent(in) :: npts,nspden
    real(dp),intent(in)  :: rho(npts,nspden)
    real(dp),intent(out) :: vxc(npts,nspden), exc(npts)
    real(dp),intent(in) :: grho2(*)
    real(dp),intent(out) :: vxcgr(*)
    real(dp),intent(out), optional :: dvxci(npts,nspden + 1)

    !Local variables-------------------------------
    integer  :: i, j, ipts, ixc, ndvxc, ngr2, nd2vxc, nvxcdgr, order
    integer, parameter :: n_blocks = 128
    integer :: ipte, nb
    real(dp) :: rhotmp(nspden, n_blocks), exctmp(n_blocks), vxctmp(nspden, n_blocks)
    real(dp) :: sigma(2*min(nspden,2)-1, n_blocks), vsigma(2*min(nspden,2)-1, n_blocks)
    real(dp) :: v2rho2(3, n_blocks), v2rhosigma(6, n_blocks), v2sigma2(6, n_blocks)
    real(dp), allocatable :: rho_(:,:), exc_(:), vxc_(:,:)
    character(len=*), parameter :: subname='xc_getvxc'

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
    else if (xc%kind == XC_LIBXC) then
       ! Pure LibXC case.
       ! WARNING: LDA implementation only, first derivative.

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
    end if

  end subroutine xc_getvxc

end module module_xc
!!***
