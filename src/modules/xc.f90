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

  type libxc_functional
     private
     integer         :: family ! LDA, GGA, etc.
     integer         :: id     ! identifier

     type(xc_f90_pointer_t) :: conf ! the pointer used to call the library
     type(xc_f90_pointer_t) :: info ! information about the functional
  end type libxc_functional

  type(libxc_functional) :: funcs(2)

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
  subroutine xc_init(ixc,nspden)

    implicit none

    !Arguments ------------------------------------
    !scalars
    integer, intent(in) :: nspden
    integer, intent(in) :: ixc

    !Local variables-------------------------------
    !scalars
    integer :: i, ii, ierr, jj
    type(xc_f90_pointer_t) :: str
    character(len=500) :: message, message2

    ! *************************************************************************

    if (ixc < 0) then
       funcs(1)%id = -ixc/1000
       funcs(2)%id = -ixc - funcs(1)%id*1000
    else if (ixc > 0) then
       ! ABINIT functional.
       funcs(1)%id = -ixc
       funcs(2)%id = 0
    else
       funcs(1)%id = 0
       funcs(2)%id = 0
    end if

    do i = 1, 2
       funcs(i)%family = 0
       if (funcs(i)%id > 0) then
          ! LibXC case.

          ! Get XC functional family
          funcs(i)%family = xc_f90_family_from_id(funcs(i)%id)
          select case (funcs(i)%family)
          case (XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
             call xc_f90_func_init(funcs(i)%conf,funcs(i)%info,funcs(i)%id,nspden)
          case default
             write(*,*) "Error: unsupported functional, change ixc."
             call MPI_ABORT(MPI_COMM_WORLD,0,ierr)
          end select

       else if (funcs(i)%id < 0) then
          ! ABINIT case

          ! Get XC functional family
          if ((-funcs(i)%id >= 1 .and. -funcs(i)%id < 11) .or. -funcs(i)%id == 24) then
             funcs(i)%family = XC_FAMILY_LDA
          else if (-funcs(i)%id >= 11 .and. -funcs(i)%id < 18) then
             funcs(i)%family = XC_FAMILY_GGA
          else if (-funcs(i)%id >= 23 .and. -funcs(i)%id < 28) then
             funcs(i)%family = XC_FAMILY_GGA
          else if (-funcs(i)%id >= 31 .and. -funcs(i)%id < 35) then
             funcs(i)%family = XC_FAMILY_LDA
          else
             write(*,*) "Error: unsupported functional, change ixc."
             call MPI_ABORT(MPI_COMM_WORLD,0,ierr)
          end if
       end if
    end do
  end subroutine xc_init

  !>  Dump XC info on screen.
  subroutine xc_dump()
    implicit none

    integer :: ii, jj
    type(xc_f90_pointer_t) :: str
    character(len=500) :: message, message2

    ! Dump functional information
    if (funcs(1)%id < 0) then
       write(*,"(1x,A41,1x,A1)") "XC functional provided by ABINIT.", "|"
    else if (any(funcs(:)%id > 0)) then
       if (funcs(1)%id > 0) then
          call xc_f90_info_name(funcs(1)%info,message)
       else
          write(message, "(A)") ""
       end if
       if (funcs(2)%id > 0) then
          call xc_f90_info_name(funcs(2)%info,message2)
       else
          write(message2, "(A)") ""
       end if
       write(*,"(1x,a41,1x,A1,1x,A40)") trim(message), "|", trim(message2)
       ii = 0
       jj = 0
       if (funcs(1)%id > 0) then
          call xc_f90_info_refs(funcs(1)%info,ii,str,message)
       else
          ii = -1
       end if
       if (funcs(2)%id > 0) then
          call xc_f90_info_refs(funcs(2)%info,jj,str,message2)
       else
          jj = -1
       end if
       do while (ii >= 0 .or. jj >= 0)
          if (ii >= 0) then
             write(*,"(1x,a41)", advance = "NO") trim(message)
             call xc_f90_info_refs(funcs(1)%info,ii,str,message)
          else
             write(*,"(1x,a41)", advance = "NO") " "
          end if
          write(*, "(1x,A1,1x)", advance = "NO") "|"
          if (jj >= 0) then
             write(*,"(a40)") trim(message2)
             call xc_f90_info_refs(funcs(2)%info,jj,str,message2)
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

    integer :: i

    do i = 1, 2
       if (funcs(i)%id > 0) call xc_f90_func_end(funcs(i)%conf)
    end do
  end subroutine xc_end

  !>  Test function to identify whether the presently used functional
  !!  is a GGA or not
  function xc_isgga()
    implicit none

    logical :: xc_isgga

    if (any(funcs%family == XC_FAMILY_GGA) .or. &
         & any(funcs%family == XC_FAMILY_HYB_GGA)) then
       xc_isgga = .true.
    else
       xc_isgga = .false.
    end if
  end function xc_isgga

  real(kind=8) function xc_exctXfac()
    implicit none

    xc_exctXfac = 0.d0
    if (any(funcs%family == XC_FAMILY_HYB_GGA)) then
       !factors for the exact exchange contribution of different hybrid functionals
       if (any(funcs%id == XC_HYB_GGA_XC_PBEH)) then
          xc_exctXfac = 0.25d0 
       end if
    end if
  end function xc_exctXfac

  subroutine xc_init_rho(n, rho, nproc)
    implicit none
    ! Arguments
    integer :: n,nproc
    real(dp) :: rho(n)

    if (any(funcs%id < 0)) then
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

    if (any(funcs%id < 0)) then
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
    real(dp),intent(in), optional :: grho2(npts,2*min(nspden,2)-1)
    real(dp),intent(out), optional :: vxcgr(npts,3)
    real(dp),intent(out), optional :: dvxci(npts,nspden + 1)

    !Local variables-------------------------------
    integer  :: i, ipts, ixc, ndvxc, ngr2, nd2vxc, nvxcdgr, order
    real(dp) :: rhotmp(nspden), exctmp, sigma(3), vsigma(3), vxctmp(nspden)
    real(dp) :: v2rho2(3), v2rhosigma(6), v2sigma2(6)
    character(len=*), parameter :: subname='xc_getvxc'

    ! Inititalize all relevant arrays to zero
    vxc=real(0,dp)
    exc=real(0,dp)
    vxctmp=real(0,dp)
    exctmp=real(0,dp)
    if (present(vxcgr)) vxcgr=real(0,dp)
    if (present(dvxci)) dvxci=real(0,dp)

    if (any(funcs%id < 0)) then
       ! ABINIT case, call drivexc
       ixc = -funcs(1)%id
       !Allocations of the exchange-correlation terms, depending on the ixc value
       order = 1
       if (present(dvxci) .and. nspden == 1) order = -2
       if (present(dvxci) .and. nspden == 2) order = +2

       call size_dvxc(ixc,ndvxc,ngr2,nd2vxc,nspden,nvxcdgr,order)

       !let us apply ABINIT routines
       select case (funcs(1)%family)
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
    else
       ! LibXC case.

       !Loop over points
       do ipts = 1, npts

          ! Convert the quantities provided by ABINIT to the ones needed by libxc
          if (nspden == 1) then
             ! ABINIT passes rho_up in the spin-unpolarized case, while the libxc
             ! expects the total density
             rhotmp(1:nspden) = real(2,dp)*rho(ipts,1:nspden)
          else
             rhotmp(1:nspden) = rho(ipts,1:nspden)
          end if
          if (xc_isgga()) then
             sigma=real(0,dp)
             if (present(grho2)) then
                if (nspden==1) then
                   ! ABINIT passes |rho_up|^2 while Libxc needs |rho_tot|^2
                   sigma(1) = real(4,dp)*grho2(ipts,1)
                else
                   ! ABINIT passes |rho_up|^2, |rho_dn|^2, and |rho_tot|^2
                   ! while Libxc needs |rho_up|^2, rho_up.rho_dn, and |rho_dn|^2
                   sigma(1) = grho2(ipts,1)
                   sigma(2) = (grho2(ipts,3) - grho2(ipts,1) - grho2(ipts,2))/real(2,dp)
                   sigma(3) = grho2(ipts,2)
                end if
             end if
          end if

          !Loop over functionals
          do i = 1,2
             if (funcs(i)%id == 0) cycle

             !Get the potential (and possibly the energy)
             if (iand(xc_f90_info_flags(funcs(i)%info), XC_FLAGS_HAVE_EXC) .ne. 0) then
                select case (funcs(i)%family)
                case (XC_FAMILY_LDA)
                   call xc_f90_lda_exc_vxc(funcs(i)%conf,1,rhotmp(1),exctmp,vxctmp(1))
                case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
                   call xc_f90_gga_exc_vxc(funcs(i)%conf,1,rhotmp(1),sigma(1),exctmp,vxctmp(1),vsigma(1))
                end select

             else
                exctmp=real(0,dp)
                select case (funcs(i)%family)
                case (XC_FAMILY_LDA)
                   call xc_f90_lda_vxc(funcs(i)%conf,1,rhotmp(1),vxctmp(1))
                case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
                   call xc_f90_gga_vxc(funcs(i)%conf,1,rhotmp(1),sigma(1),vxctmp(1),vsigma(1))
                end select
             end if
             if (present(dvxci)) then
                v2rho2     = real(0, dp)
                v2rhosigma = real(0, dp)
                v2sigma2   = real(0, dp)
                if (iand(xc_f90_info_flags(funcs(i)%info), XC_FLAGS_HAVE_FXC) .ne. 0) then
                   select case (funcs(i)%family)
                   case (XC_FAMILY_LDA)
                      call xc_f90_lda_fxc(funcs(i)%conf,1,rhotmp(1),v2rho2(1))
                   case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
                      call xc_f90_gga_fxc(funcs(i)%conf,1,rhotmp(1),sigma(1), &
                           & v2rho2(1),v2rhosigma(1), v2sigma2(1))
                   end select
                end if
             end if

             exc(ipts) = exc(ipts) + exctmp
             vxc(ipts,1:nspden) = vxc(ipts,1:nspden) + vxctmp(1:nspden)

             if (xc_isgga() .and. present(vxcgr)) then
                !Convert the quantities returned by Libxc to the ones needed by ABINIT
                if (nspden == 1) then
                   vxcgr(ipts,3) = vxcgr(ipts,3) + vsigma(1)*real(2,dp)
                else
                   vxcgr(ipts,1) = vxcgr(ipts,1) + real(2,dp)*vsigma(1) - vsigma(2)
                   vxcgr(ipts,2) = vxcgr(ipts,2) + real(2,dp)*vsigma(3) - vsigma(2)
                   vxcgr(ipts,3) = vxcgr(ipts,3) + vsigma(2)
                end if
             end if

             if (present(dvxci)) then
                !Convert the quantities returned by Libxc to the ones needed by ABINIT
                if (nspden == 1) then
                   dvxci(ipts,1) = dvxci(ipts,1) + v2rho2(1)*real(2,dp)
                else
                   ! TODO
!!$                dvxci(ipts,1) = dvxci(ipts,1) + real(2,dp)*vsigma(1) - vsigma(2)
!!$                dvxci(ipts,2) = dvxci(ipts,2) + real(2,dp)*vsigma(3) - vsigma(2)
!!$                dvxci(ipts,3) = dvxci(ipts,3) + vsigma(2)
                end if
             end if

          end do
       end do
    end if

  end subroutine xc_getvxc

end module module_xc
!!***
