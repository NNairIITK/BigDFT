!{\src2tex{textfont=tt}}
!!****m* ABINIT/libxc_functionals
!! NAME
!!  libxc_functionals
!!
!! FUNCTION
!!  (to be provided)  
!!
!! COPYRIGHT
!! Copyright (C) 2008-2009 ABINIT group (MOliveira)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module libxc_functionals

  use defs_basis
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
    type(xc_f90_func_t) :: conf ! the pointer used to call the library
    type(xc_f90_info_t) :: info ! information about the functional
#endif

    integer         :: irel
    real(dp)        :: xalpha
  end type libxc_functional

  type(libxc_functional) :: funcs(2)

  private
  public :: libxc_functionals_init, &
&      libxc_functionals_getvxc, &
&      libxc_functionals_isgga, &
&      libxc_functionals_end

contains
!!*** 

!!****f* libxc_functionals/libxc_functionals_init
!! NAME
!!  libxc_functionals_init
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      xc_f90_gga_vxc,xc_f90_lda_vxc
!!
!! SOURCE

  subroutine libxc_functionals_init(ixc,nspden)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

    implicit none

!Arguments ------------------------------------
!scalars

    integer, intent(in) :: nspden
    integer, intent(in) :: ixc

!Local variables-------------------------------
!scalars

    integer :: i
    character(len=500) :: message

! *************************************************************************

#if defined HAVE_LIBXC
    funcs(1)%id = -ixc/1000
    funcs(2)%id = -ixc - funcs(1)%id*1000

    do i = 1, 2
      if (funcs(i)%id == 0) then
        funcs(i)%family = 0
        cycle
      end if

      ! Get XC functional family
      funcs(i)%family = xc_f90_family_from_id(funcs(i)%id)

      ! Extra variables
      if (funcs(i)%family == XC_FAMILY_LDA .and. funcs(i)%id == XC_LDA_X) then
        !Here we should decide if we use the relativistic corrections to the exchange part of the LDA. For now we will leave it as XC_NON_RELATIVISTIC
        funcs(i)%irel = XC_NON_RELATIVISTIC
      end if
      if (funcs(i)%family == XC_FAMILY_LDA .and. funcs(i)%id == XC_LDA_C_XALPHA) then
        !Here we should get the value of xalpha. For now we will set it to zero
        funcs(i)%xalpha=zero
      end if

      ! Init LibXC
      select case(funcs(i)%family)
      case(XC_FAMILY_LDA)
        if (funcs(i)%id == XC_LDA_X) then
          call xc_f90_lda_init(funcs(i)%conf,funcs(i)%info,funcs(i)%id,nspden,3,funcs(i)%irel)
        elseif (funcs(i)%id == XC_LDA_C_XALPHA) then
          call xc_f90_lda_init(funcs(i)%conf,funcs(i)%info,XC_LDA_C_XALPHA,nspden,3,funcs(i)%xalpha)
        else
          call xc_f90_lda_init(funcs(i)%conf,funcs(i)%info,funcs(i)%id,nspden)
        end if
      case(XC_FAMILY_GGA)
        call xc_f90_gga_init(funcs(i)%conf,funcs(i)%info,funcs(i)%id,nspden)
      case default
        write(message, '(4a,i8,2a,i8,6a)' )ch10,&
             &    ' libxc_functionals_init : ERROR -',ch10,&
             &    '  Invalid IXC = ',ixc,ch10,&
             &    '  The LibXC functional family ',funcs(i)%family,&
             &    '  is currently unsupported by ABINIT',ch10,&
             &    '  (-1 means the family is unknown to the LibXC itself)',ch10,&
             &    '  Please consult the LibXC documentation',ch10
        call wrtout(std_out,message,'COLL')
        call leave_new('COLL')
      end select

      ! Dump functional information
      call xc_f90_info_name(funcs(i)%info,message)
      call wrtout(std_out,message,'COLL')

    end do
#else
    write(message, '(a,a,a,a)' ) ch10,&
         & ' wvl_init_type_wfs : LibXC library is not compiled.', ch10, &
         & '   Action, used the flag --enable-libxc when configuring.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
#endif

  end subroutine libxc_functionals_init
!!***

!!****f* libxc_functionals/libxc_functionals_end
!! NAME
!!  libxc_functionals_end
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      xc_f90_gga_vxc,xc_f90_lda_vxc
!!
!! SOURCE
  subroutine libxc_functionals_end()


    implicit none

    integer :: i
    character(len=500) :: message

#if defined HAVE_LIBXC
    do i = 1, 2
      if (funcs(i)%id == 0) cycle
      select case (funcs(i)%family)
      case (XC_FAMILY_LDA)
        call xc_f90_lda_end(funcs(i)%conf)
      case (XC_FAMILY_GGA)
        call xc_f90_gga_end(funcs(i)%conf)
      end select
    end do
#else
    write(message, '(a,a,a,a)' ) ch10,&
         & ' wvl_init_type_wfs : LibXC library is not compiled.', ch10, &
         & '   Action, used the flag --enable-libxc when configuring.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
#endif

  end subroutine libxc_functionals_end
!!*** 

!!****f* libxc_functionals/libxc_functionals_isgga
!! NAME
!!  libxc_functionals_isgga
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

  function libxc_functionals_isgga()


    implicit none

!Arguments ------------------------------------

!Local variables-------------------------------

    logical :: libxc_functionals_isgga
    character(len=500) :: message

! *************************************************************************

#if defined HAVE_LIBXC
    if (any(funcs%family == XC_FAMILY_GGA)) then
      libxc_functionals_isgga = .true.
    else
      libxc_functionals_isgga = .false.
    end if
#else
    libxc_functionals_isgga = .false.
    write(message, '(a,a,a,a)' ) ch10,&
         & ' wvl_init_type_wfs : LibXC library is not compiled.', ch10, &
         & '   Action, used the flag --enable-libxc when configuring.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
#endif

  end function libxc_functionals_isgga
!!*** 

!!****f* libxc_functionals/libxc_functionals_getvxc
!! NAME
!!  libxc_functionals_getvxc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!      xc_f90_gga_vxc,xc_f90_lda_vxc
!!
!! SOURCE

  subroutine libxc_functionals_getvxc(exc,nspden,rho,vxc,grho2,vxcgr)


    implicit none

!Arguments ------------------------------------

    integer, intent(in) :: nspden
    real(dp),intent(in)  :: rho(nspden)
    real(dp),intent(out) :: vxc(nspden), exc
    real(dp),intent(in),optional :: grho2(2*min(nspden,2)-1)
    real(dp),intent(out),optional :: vxcgr(3)

!Local variables-------------------------------

    integer  :: i
    real(dp) :: rhotmp(nspden), exctmp, sigma(3), vsigma(3), vxctmp(nspden)
    character(len=500) :: message

! *************************************************************************

#if defined HAVE_LIBXC
    ! Inititalize all relevant arrays to zero
    vxc=zero
    exc=zero
    vxctmp=zero
    exctmp=zero
    if (any(funcs%family == XC_FAMILY_GGA)) vxcgr=zero

    ! ABINIT passes rho_up in the spin-unpolarized case, while the libxc expects the total density
    if (nspden == 1) then
      rhotmp(1:nspden) = two*rho(1:nspden)
    else
      rhotmp(1:nspden) = rho(1:nspden)
    end if

    do i = 1,2
      if (funcs(i)%id == 0) cycle

      select case (funcs(i)%family)
      case (XC_FAMILY_LDA)
        call xc_f90_lda_vxc(funcs(i)%conf,rhotmp(1),exctmp,vxctmp(1))

      case (XC_FAMILY_GGA)
        sigma=zero
        if (nspden==1) then
          ! ABINIT passes |rho_up|^2 while Libxc needs |rho_tot|^2
          sigma(1) = four*grho2(1)
        else
          ! ABINIT passes |rho_up|^2, |rho_dn|^2, and |rho_tot|^2
          ! while Libxc needs |rho_up|^2, rho_up.rho_dn, and |rho_dn|^2
          sigma(1) = grho2(1)
          sigma(2) = (grho2(3) - grho2(1) - grho2(2))/two
          sigma(3) = grho2(2)
        end if
        call xc_f90_gga_vxc(funcs(i)%conf,rhotmp(1),sigma(1),exctmp,vxctmp(1),vsigma(1))
        if (nspden == 1) then
          vxcgr(3) = vxcgr(3) + vsigma(1)*two
        else
          vxcgr(1) = vxcgr(1) + two*vsigma(1) - vsigma(2)
          vxcgr(2) = vxcgr(2) + two*vsigma(3) - vsigma(2)
          vxcgr(3) = vxcgr(3) + vsigma(2)
        end if

      end select
      exc=exc+exctmp
      vxc(1:nspden)=vxc(1:nspden)+vxctmp(1:nspden)

    end do
#else
    write(message, '(a,a,a,a)' ) ch10,&
         & ' wvl_init_type_wfs : LibXC library is not compiled.', ch10, &
         & '   Action, used the flag --enable-libxc when configuring.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
#endif

  end subroutine libxc_functionals_getvxc

end module 
!!***
