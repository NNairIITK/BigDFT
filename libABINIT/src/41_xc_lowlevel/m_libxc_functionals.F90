!{\src2tex{textfont=tt}}
!!****m* ABINIT/libxc_functionals
!! NAME
!!  libxc_functionals
!!
!! FUNCTION
!!  Module containing interfaces to the LibXC library, for exchange
!!  correlation potentials and energies. The interfacing between
!!  the ABINIT and LibXC formats and datastructures happens here.
!!  Also contains basic container datatype for LibXC interfacing.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2010 ABINIT group (MOliveira)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module libxc_functionals

 use defs_basis
 use abi_interfaces_lowlevel

#if defined HAVE_LIBXC
  use xc_f90_types_m
  use libxc_funcs_m
  use xc_f90_lib_m
#endif

  implicit none

  type libxc_functional
    private
    integer         :: family ! LDA, GGA, etc.
    integer         :: id     ! identifier

#if defined HAVE_LIBXC
    type(xc_f90_pointer_t) :: conf ! the pointer used to call the library
    type(xc_f90_pointer_t) :: info ! information about the functional
#endif
  end type libxc_functional

  type(libxc_functional) :: funcs(2)

  private
  public :: libxc_functionals_init, &
&      libxc_functionals_getvxc, &
&      libxc_functionals_isgga, &
&      libxc_functionals_ismgga, &
&      libxc_functionals_exctXfac, &
&      libxc_functionals_end

contains
!!*** 

!!****f* libxc_functionals/libxc_functionals_init
!! NAME
!!  libxc_functionals_init
!!
!! FUNCTION
!!  Initialize the desired XC functional, from LibXC.
!!  * Call the LibXC initializer
!!  * Fill preliminary fields in module structures.
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! SOURCE

  subroutine libxc_functionals_init(ixc,nspden)

    implicit none

!Arguments ------------------------------------
!scalars

    integer, intent(in) :: nspden
    integer, intent(in) :: ixc

!Local variables-------------------------------
!scalars

    integer :: i, ii
    character(len=500) :: message
#if defined HAVE_LIBXC
    type(xc_f90_pointer_t) :: str
#endif
! *************************************************************************

#if defined HAVE_LIBXC
    if (ixc < 0) then
       funcs(1)%id = -ixc/1000
       funcs(2)%id = -ixc - funcs(1)%id*1000
    else
       funcs(1)%id = 0
       funcs(2)%id = 0
    end if

    do i = 1, 2
      if (funcs(i)%id == 0) then
        funcs(i)%family = 0
        cycle
      end if

      ! Get XC functional family
      funcs(i)%family = xc_f90_family_from_id(funcs(i)%id)
      select case (funcs(i)%family)
      case (XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_HYB_GGA)
        call xc_f90_func_init(funcs(i)%conf,funcs(i)%info,funcs(i)%id,nspden)
      case default
        write(message, '(4a,i8,2a,i8,6a)' )ch10,&
             &    ' libxc_functionals_init : ERROR -',ch10,&
             &    '  Invalid IXC = ',ixc,ch10,&
             &    '  The LibXC functional family ',funcs(i)%family,&
             &    '  is currently unsupported by ABINIT',ch10,&
             &    '  (-1 means the family is unknown to the LibXC itself)',ch10,&
             &    '  Please consult the LibXC documentation',ch10
        call abi_wrtout(std_out,message,'COLL')
        call abi_leave_new('COLL')
      end select

      ! Dump functional information
      call xc_f90_info_name(funcs(i)%info,message)
      call abi_wrtout(std_out,message,'COLL')
      ii = 0
      call xc_f90_info_refs(funcs(i)%info,ii,str,message)
      do while (ii >= 0)
        call abi_wrtout(std_out,message,'COLL')
        call xc_f90_info_refs(funcs(i)%info,ii,str,message)
      end do
    end do
#else
    write(message, '(a,a,a,a)' ) ch10,&
         & ' wvl_init_type_wfs : LibXC library is not compiled.', ch10, &
         & '   Action, used the flag --enable-libxc when configuring.'
    call abi_wrtout(6,message,'COLL')
    call abi_leave_new('COLL')
#endif
  end subroutine libxc_functionals_init
!!***

!!****f* libxc_functionals/libxc_functionals_end
!! NAME
!!  libxc_functionals_end
!!
!! FUNCTION
!!  End usage of LibXC functional. Call LibXC end function,
!!  and deallocate module contents.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

  subroutine libxc_functionals_end()


    implicit none

    integer :: i
#ifndef HAVE_LIBXC
   character(len=500) :: message
#endif

#if defined HAVE_LIBXC
    do i = 1, 2
      if (funcs(i)%id == 0) cycle
      call xc_f90_func_end(funcs(i)%conf)
    end do
#else
    write(message, '(a,a,a,a)' ) ch10,&
         & ' wvl_init_type_wfs : LibXC library is not compiled.', ch10, &
         & '   Action, used the flag --enable-libxc when configuring.'
    call abi_wrtout(6,message,'COLL')
    call abi_leave_new('COLL')
#endif
  end subroutine libxc_functionals_end
!!*** 

!!****f* libxc_functionals/libxc_functionals_isgga
!! NAME
!!  libxc_functionals_isgga
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  is a GGA or not
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! SOURCE

  function libxc_functionals_isgga()

    implicit none

!Arguments ------------------------------------

!Local variables-------------------------------

    logical :: libxc_functionals_isgga

! *************************************************************************

#if defined HAVE_LIBXC
    if (any(funcs%family == XC_FAMILY_GGA) .or. any(funcs%family == XC_FAMILY_HYB_GGA)) then
      libxc_functionals_isgga = .true.
    else
      libxc_functionals_isgga = .false.
    end if
#else
    libxc_functionals_isgga = .false.
#endif
  end function libxc_functionals_isgga
!!*** 

!!****f* libxc_functionals/libxc_functionals_exctXfac
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! SOURCE

  real(kind=8) function libxc_functionals_exctXfac()


    implicit none

!Arguments ------------------------------------

!Local variables-------------------------------

! *************************************************************************

    libxc_functionals_exctXfac = 0.d0
#if defined HAVE_LIBXC
    if (any(funcs%family == XC_FAMILY_HYB_GGA)) then
       !factors for the exact exchange contribution of different hybrid functionals
       if (any(funcs%id == XC_HYB_GGA_XC_PBEH)) then
          libxc_functionals_exctXfac = 0.25d0 
       end if
    end if
#endif

  end function libxc_functionals_exctXfac
!!*** 

!!****f* libxc_functionals/libxc_functionals_ismgga
!! NAME
!!  libxc_functionals_ismgga
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  is a Meta-GGA or not
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

  function libxc_functionals_ismgga()

    implicit none

!Arguments ------------------------------------

!Local variables-------------------------------

    logical :: libxc_functionals_ismgga

! *************************************************************************
#if defined HAVE_LIBXC
    if (any(funcs%family == XC_FAMILY_MGGA)) then
      libxc_functionals_ismgga = .true.
    else
      libxc_functionals_ismgga = .false.
    end if
#else
    libxc_functionals_ismgga = .false.
#endif
  end function libxc_functionals_ismgga
!!*** 

!!****f* libxc_functionals/libxc_functionals_getvxc
!! NAME
!!  libxc_functionals_getvxc
!!
!! FUNCTION
!!  Return XC potential and energy, from input density (event gradient etc...)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

  subroutine libxc_functionals_getvxc(npts,exc,nspden,rho,vxc,grho2,vxcgr,lrho,vxclrho,tau,vxctau)

    implicit none

!Arguments ------------------------------------

    integer, intent(in) :: npts,nspden
    real(dp),intent(in)  :: rho(npts,nspden)
    real(dp),intent(out) :: vxc(npts,nspden), exc(npts)
    real(dp),intent(in),optional :: grho2(npts,2*min(nspden,2)-1)
    real(dp),intent(out),optional :: vxcgr(npts,3)
    real(dp),intent(in),optional :: lrho(npts,nspden)
    real(dp),intent(out),optional :: vxclrho(npts,nspden)
    real(dp),intent(in),optional :: tau(npts,nspden)
    real(dp),intent(out),optional :: vxctau(npts,nspden)

!Local variables-------------------------------

    integer  :: i, ipts
    real(dp) :: rhotmp(nspden), exctmp, sigma(3), vsigma(3), vxctmp(nspden)
    real(dp) :: lrhotmp(nspden), tautmp(nspden), vxclrhotmp(nspden), vxctautmp(nspden)
!  character(len=500) :: message

! *************************************************************************

#if defined HAVE_LIBXC
    ! Inititalize all relevant arrays to zero
    vxc=zero
    exc=zero
    vxctmp=zero
    exctmp=zero
    if (any(funcs%family == XC_FAMILY_GGA)) vxcgr=zero
    if (any(funcs%family == XC_FAMILY_MGGA)) then
      vxcgr=zero
      vxclrho=zero
      vxctau=zero
    end if

    !Loop over points
    do ipts = 1, npts

      ! Convert the quantities provided by ABINIT to the ones needed by libxc
      if (nspden == 1) then
        ! ABINIT passes rho_up in the spin-unpolarized case, while the libxc
        ! expects the total density
        rhotmp(1:nspden) = two*rho(ipts,1:nspden)
      else
        rhotmp(1:nspden) = rho(ipts,1:nspden)
      end if
      if (libxc_functionals_isgga()) then
        sigma=zero
        if (nspden==1) then
          ! ABINIT passes |rho_up|^2 while Libxc needs |rho_tot|^2
          sigma(1) = four*grho2(ipts,1)
        else
          ! ABINIT passes |rho_up|^2, |rho_dn|^2, and |rho_tot|^2
          ! while Libxc needs |rho_up|^2, rho_up.rho_dn, and |rho_dn|^2
          sigma(1) = grho2(ipts,1)
          sigma(2) = (grho2(ipts,3) - grho2(ipts,1) - grho2(ipts,2))/two
          sigma(3) = grho2(ipts,2)
        end if
      end if
      if (any(funcs%family == XC_FAMILY_MGGA)) then
        if (nspden==1) then
          lrhotmp(1:nspden) = two*lrho(ipts,1:nspden)
          tautmp(1:nspden) = four*tau(ipts,1:nspden)
        else

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
          case (XC_FAMILY_GGA)
            call xc_f90_gga_exc_vxc(funcs(i)%conf,1,rhotmp(1),sigma(1),exctmp,vxctmp(1),vsigma(1))
          case (XC_FAMILY_HYB_GGA)
            call xc_f90_gga_exc_vxc(funcs(i)%conf,1,rhotmp(1),sigma(1),exctmp,vxctmp(1),vsigma(1))
          case (XC_FAMILY_MGGA)
            call xc_f90_mgga_exc_vxc(funcs(i)%conf,1,rhotmp(1),sigma(1),lrhotmp(1),&
                        tautmp(1),exctmp,vxctmp(1),vsigma(1),vxclrhotmp(1),vxctautmp(1))
          end select

        else
          exctmp=zero
          select case (funcs(i)%family)
          case (XC_FAMILY_LDA)
             call xc_f90_lda_vxc(funcs(i)%conf,1,rhotmp(1),vxctmp(1))
          case (XC_FAMILY_GGA)
             call xc_f90_gga_vxc(funcs(i)%conf,1,rhotmp(1),sigma(1),vxctmp(1),vsigma(1))
          case (XC_FAMILY_HYB_GGA)
             call xc_f90_gga_vxc(funcs(i)%conf,1,rhotmp(1),sigma(1),vxctmp(1),vsigma(1))
          case (XC_FAMILY_MGGA)
             call xc_f90_mgga_vxc(funcs(i)%conf,1,rhotmp(1),sigma(1),lrhotmp(1),&
                  tautmp(1),vxctmp(1),vsigma(1),vxclrhotmp(1),vxctautmp(1))
          end select
        end if

        exc(ipts) = exc(ipts) + exctmp
        vxc(ipts,1:nspden) = vxc(ipts,1:nspden) + vxctmp(1:nspden)

        if (libxc_functionals_isgga() .or. any(funcs%family == XC_FAMILY_MGGA)) then
          !Convert the quantities returned by Libxc to the ones needed by ABINIT
          if (nspden == 1) then
            vxcgr(ipts,3) = vxcgr(ipts,3) + vsigma(1)*two
          else
            vxcgr(ipts,1) = vxcgr(ipts,1) + two*vsigma(1) - vsigma(2)
            vxcgr(ipts,2) = vxcgr(ipts,2) + two*vsigma(3) - vsigma(2)
            vxcgr(ipts,3) = vxcgr(ipts,3) + vsigma(2)
          end if
        end if
        if (any(funcs%family == XC_FAMILY_MGGA)) then
          vxclrho(ipts,1:nspden) = vxclrho(ipts,1:nspden) + vxclrhotmp(1:nspden)
          vxctau(ipts,1:nspden) = vxctau(ipts,1:nspden) + vxctautmp(1:nspden)
        end if

      end do

    end do

#endif
  end subroutine libxc_functionals_getvxc

end module 
!!***
