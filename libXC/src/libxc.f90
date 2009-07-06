# 1 "./libxc_master.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "./libxc_master.F90"
!! Copyright (C) 2003-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: libxc.f90 3550 2007-11-19 14:32:49Z marques $







module xc_f90_types_m



  integer, public, parameter :: xc_f90_kind = selected_real_kind(14)


  type xc_f90_func_t
    private
    integer, pointer :: func
  end type xc_f90_func_t

  type xc_f90_info_t
    private
    integer, pointer :: info
  end type xc_f90_info_t

end module xc_f90_types_m

module xc_f90_lib_m

  use xc_f90_types_m
  use libxc_funcs_m

  implicit none

  public
  public :: &
    xc_f90_func_t, &
    xc_f90_info_t, &
    xc_f90_info_number, &
    xc_f90_info_kind, &
    xc_f90_info_name, &
    xc_f90_info_family, &
    xc_f90_info_refs, &
    xc_f90_family_from_id, &
    xc_f90_lda_init, &
    xc_f90_lda, &
    xc_f90_lda_exc, &
    xc_f90_lda_vxc, &
    xc_f90_lda_fxc, &
    xc_f90_lda_kxc, &
    xc_f90_lda_end, &
    xc_f90_lda_c_1d_csc_set_params, &
    xc_f90_lda_c_2d_prm_set_params, &
    xc_f90_gga_init, &
    xc_f90_gga, &
    xc_f90_gga_exc, &
    xc_f90_gga_vxc, &
    xc_f90_gga_fxc, &
    xc_f90_gga_end, &
    xc_f90_gga_lb_set_params, &
    xc_f90_gga_lb_modified, &
    xc_f90_hyb_gga_init, &
    xc_f90_hyb_gga_end, &
    xc_f90_hyb_gga, &
    xc_f90_hyb_gga_exc, &
    xc_f90_hyb_gga_vxc, &
    xc_f90_hyb_gga_fxc, &
    xc_f90_hyb_gga_exx_coef, &
    xc_f90_lca_init, &
    xc_f90_lca_end, &
    xc_f90_lca, &
    xc_f90_mgga_init, &
    xc_f90_mgga, &
    xc_f90_mgga_exc, &
    xc_f90_mgga_vxc, &
    xc_f90_mgga_fxc, &
    xc_f90_mgga_end

  ! Families of xc functionals
  integer, public, parameter :: &
    XC_FAMILY_UNKNOWN = -1, &
    XC_FAMILY_NONE = 0, &
    XC_FAMILY_LDA = 1, &
    XC_FAMILY_GGA = 2, &
    XC_FAMILY_MGGA = 4, &
    XC_FAMILY_LCA = 8, &
    XC_FAMILY_OEP = 16, &
    XC_FAMILY_HYB_GGA = 32

  integer, public, parameter :: &
    XC_UNPOLARIZED = 1, & ! Spin unpolarized
    XC_POLARIZED = 2 ! Spin polarized

  integer, public, parameter :: &
    XC_NON_RELATIVISTIC = 0, & ! Functional includes or not relativistic
    XC_RELATIVISTIC = 1 ! corrections. Only available in some functionals.

  ! Kinds
  integer, public, parameter :: &
    XC_EXCHANGE = 0, &
    XC_CORRELATION = 1, &
    XC_EXCHANGE_CORRELATION = 2


  !----------------------------------------------------------------
  interface
    integer function xc_f90_info_number(info)
      use xc_f90_types_m
      type(xc_f90_info_t), intent(in) :: info
    end function xc_f90_info_number

    integer function xc_f90_info_kind(info)
      use xc_f90_types_m
      type(xc_f90_info_t), intent(in) :: info
    end function xc_f90_info_kind

    subroutine xc_f90_info_name(info, s)
      use xc_f90_types_m
      type(xc_f90_info_t), intent(in) :: info
      character(len=*), intent(out) :: s
    end subroutine xc_f90_info_name

    integer function xc_f90_info_family(info)
      use xc_f90_types_m
      type(xc_f90_info_t), intent(in) :: info
    end function xc_f90_info_family

    subroutine xc_f90_info_refs(info, n, s)
      use xc_f90_types_m
      type(xc_f90_info_t), intent(in) :: info
      integer, intent(inout) :: n
      character(len=*), intent(out) :: s
    end subroutine xc_f90_info_refs

    integer function xc_f90_family_from_id(id)
      use xc_f90_types_m
      integer, intent(in) :: id
    end function xc_f90_family_from_id
  end interface


  ! LDAs
  ! We will use the same public interface (xc_lda_init) for the 3 C procedures
  !----------------------------------------------------------------
  interface xc_f90_lda_init
    subroutine xc_f90_lda_init_(p, info, functional, nspin)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(out) :: p
      type(xc_f90_info_t), intent(out) :: info
      integer, intent(in) :: functional
      integer, intent(in) :: nspin
    end subroutine xc_f90_lda_init_

    subroutine xc_f90_lda_x_init(p, info, functional, nspin, dim, irel)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(out) :: p
      type(xc_f90_info_t), intent(out) :: info
      integer, intent(in) :: functional
      integer, intent(in) :: nspin ! XC_UNPOLARIZED or XC_POLARIZED
      integer, intent(in) :: dim ! 2 or 3 dimensions
      integer, intent(in) :: irel ! XC_NON_RELATIVISTIC or XC_RELATIVISTIC
    end subroutine xc_f90_lda_x_init

    subroutine xc_f90_lda_c_xalpha_init(p, info, functional, nspin, dim, alpha)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(out) :: p
      type(xc_f90_info_t), intent(out) :: info
      integer, intent(in) :: functional
      integer, intent(in) :: nspin ! XC_UNPOLARIZED or XC_POLARIZED
      integer, intent(in) :: dim ! 2 or 3 dimensions
      real(xc_f90_kind), intent(in) :: alpha ! Ec = alpha Ex
    end subroutine xc_f90_lda_c_xalpha_init
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_lda_end(p)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(inout) :: p
    end subroutine xc_f90_lda_end
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_lda(p, rho, zk, vrho, fxc, kxc)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(out) :: zk ! the energy per unit particle
      real(xc_f90_kind), intent(out) :: vrho ! v(nspin) the potential
      real(xc_f90_kind), intent(out) :: fxc ! v(nspin,nspin) the xc kernel
      real(xc_f90_kind), intent(out) :: kxc ! v(nspin,nspin,nspin) the derivative of xc kernel
    end subroutine xc_f90_lda

    subroutine xc_f90_lda_exc(p, rho, zk)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(out) :: zk ! the energy per unit particle
    end subroutine xc_f90_lda_exc

    subroutine xc_f90_lda_vxc(p, rho, e, v)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(out) :: e ! the energy per unit particle
      real(xc_f90_kind), intent(out) :: v ! v(nspin) the potential
    end subroutine xc_f90_lda_vxc

    subroutine xc_f90_lda_fxc(p, rho, fxc)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(out) :: fxc ! v(nspin,nspin) the xc kernel
    end subroutine xc_f90_lda_fxc

    subroutine xc_f90_lda_kxc(p, rho, kxc)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(out) :: kxc
    end subroutine xc_f90_lda_kxc
  end interface


  interface
    subroutine xc_f90_lda_c_1d_csc_set_params(p, bb)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(inout) :: p
      real(xc_f90_kind), intent(in) :: bb
    end subroutine xc_f90_lda_c_1d_csc_set_params
  end interface

  interface
    subroutine xc_f90_lda_c_2d_prm_set_params(p, N)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(inout) :: p
      real(xc_f90_kind), intent(in) :: N
    end subroutine xc_f90_lda_c_2d_prm_set_params
  end interface

  ! GGAs
  !----------------------------------------------------------------
  interface
    subroutine xc_f90_gga_init(p, info, functional, nspin)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(out) :: p
      type(xc_f90_info_t), intent(out) :: info
      integer, intent(in) :: functional
      integer, intent(in) :: nspin
    end subroutine xc_f90_gga_init
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_gga_end(p)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(inout) :: p
    end subroutine xc_f90_gga_end
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_gga(p, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(out) :: zk
      real(xc_f90_kind), intent(out) :: vrho
      real(xc_f90_kind), intent(out) :: vsigma
      real(xc_f90_kind), intent(out) :: v2rho2
      real(xc_f90_kind), intent(out) :: v2rhosigma
      real(xc_f90_kind), intent(out) :: v2sigma2
    end subroutine xc_f90_gga
  end interface


  interface
    subroutine xc_f90_gga_exc(p, rho, sigma, zk)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(out) :: zk
    end subroutine xc_f90_gga_exc
  end interface

  interface
    subroutine xc_f90_gga_vxc(p, rho, sigma, zk, vrho, vsigma)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(out) :: zk
      real(xc_f90_kind), intent(out) :: vrho
      real(xc_f90_kind), intent(out) :: vsigma
    end subroutine xc_f90_gga_vxc
  end interface

  interface
    subroutine xc_f90_gga_fxc(p, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(out) :: v2rho2
      real(xc_f90_kind), intent(out) :: v2rhosigma
      real(xc_f90_kind), intent(out) :: v2sigma2
    end subroutine xc_f90_gga_fxc
  end interface

  !----------------------------------------------------------------
  interface
    subroutine xc_f90_gga_lb_set_params(p, modified, threshold, ip, qtot)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      integer, intent(in) :: modified ! should we use the modified version
      real(xc_f90_kind), intent(in) :: threshold ! if so, the threshold to use the asymtotic version
      real(xc_f90_kind), intent(in) :: ip ! ionization potential
      real(xc_f90_kind), intent(in) :: qtot ! total charge
    end subroutine xc_f90_gga_lb_set_params
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_gga_lb_modified(p, rho, grho, r, dedd)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(in) :: grho ! grho(3,nspin) the gradient of the density
      real(xc_f90_kind), intent(in) :: r ! distance from center of finite system
      real(xc_f90_kind), intent(out) :: dedd
    end subroutine xc_f90_gga_lb_modified
  end interface


  ! Hybrids GGAs
  !----------------------------------------------------------------
  interface
    subroutine xc_f90_hyb_gga_init(p, info, functional, nspin)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(out) :: p
      type(xc_f90_info_t), intent(out) :: info
      integer, intent(in) :: functional
      integer, intent(in) :: nspin
    end subroutine xc_f90_hyb_gga_init
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_hyb_gga_end(p)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(inout) :: p
    end subroutine xc_f90_hyb_gga_end
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_hyb_gga(p, rho, grho, e, dedd, dedgd)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(in) :: grho ! grho(3,nspin) the gradient of the density
      real(xc_f90_kind), intent(out) :: e ! the energy per unit particle
      real(xc_f90_kind), intent(out) :: dedd ! dedd(nspin) the derivative of the energy
                                            ! in terms of the density
      real(xc_f90_kind), intent(out) :: dedgd ! and in terms of the gradient of the density
    end subroutine xc_f90_hyb_gga
  end interface

  interface
    subroutine xc_f90_hyb_gga_exc(p, rho, sigma, zk)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(out) :: zk
    end subroutine xc_f90_hyb_gga_exc
  end interface

  interface
    subroutine xc_f90_hyb_gga_vxc(p, rho, sigma, zk, vrho, vsigma)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(out) :: zk
      real(xc_f90_kind), intent(out) :: vrho
      real(xc_f90_kind), intent(out) :: vsigma
    end subroutine xc_f90_hyb_gga_vxc
  end interface

  interface
    subroutine xc_f90_hyb_gga_fxc(p, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(out) :: v2rho2
      real(xc_f90_kind), intent(out) :: v2rhosigma
      real(xc_f90_kind), intent(out) :: v2sigma2
    end subroutine xc_f90_hyb_gga_fxc
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_hyb_gga_exx_coef(p, coef)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(out) :: coef
    end subroutine xc_f90_hyb_gga_exx_coef
  end interface


  ! the meta-GGAs
  !----------------------------------------------------------------
  interface
    subroutine xc_f90_mgga_init(p, info, functional, nspin)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(out) :: p
      type(xc_f90_info_t), intent(out) :: info
      integer, intent(in) :: functional
      integer, intent(in) :: nspin
    end subroutine xc_f90_mgga_init
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_mgga_end(p)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(inout) :: p
    end subroutine xc_f90_mgga_end
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_mgga(p, rho, sigma, tau, zk, vrho, vsigma, vtau, &
      v2rho2, v2rhosigma, v2sigma2, v2rhotau, v2tausigma, v2tau2)

      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(in) :: tau
      real(xc_f90_kind), intent(out) :: zk
      real(xc_f90_kind), intent(out) :: vrho
      real(xc_f90_kind), intent(out) :: vsigma
      real(xc_f90_kind), intent(out) :: vtau
      real(xc_f90_kind), intent(out) :: v2rho2
      real(xc_f90_kind), intent(out) :: v2rhosigma
      real(xc_f90_kind), intent(out) :: v2sigma2
      real(xc_f90_kind), intent(out) :: v2rhotau
      real(xc_f90_kind), intent(out) :: v2tausigma
      real(xc_f90_kind), intent(out) :: v2tau2
    end subroutine xc_f90_mgga
  end interface


  interface
    subroutine xc_f90_mgga_exc(p, rho, sigma, tau, zk)

      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(in) :: tau
      real(xc_f90_kind), intent(out) :: zk
    end subroutine xc_f90_mgga_exc
  end interface


  interface
    subroutine xc_f90_mgga_vxc(p, rho, sigma, tau, zk, vrho, vsigma, vtau)

      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(in) :: tau
      real(xc_f90_kind), intent(out) :: zk
      real(xc_f90_kind), intent(out) :: vrho
      real(xc_f90_kind), intent(out) :: vsigma
      real(xc_f90_kind), intent(out) :: vtau
    end subroutine xc_f90_mgga_vxc
  end interface


  interface
    subroutine xc_f90_mgga_fxc(p, rho, sigma, tau, &
      v2rho2, v2rhosigma, v2sigma2, v2rhotau, v2tausigma, v2tau2)

      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho
      real(xc_f90_kind), intent(in) :: sigma
      real(xc_f90_kind), intent(in) :: tau
      real(xc_f90_kind), intent(out) :: v2rho2
      real(xc_f90_kind), intent(out) :: v2rhosigma
      real(xc_f90_kind), intent(out) :: v2sigma2
      real(xc_f90_kind), intent(out) :: v2rhotau
      real(xc_f90_kind), intent(out) :: v2tausigma
      real(xc_f90_kind), intent(out) :: v2tau2
    end subroutine xc_f90_mgga_fxc
  end interface

  ! the LCAs
  !----------------------------------------------------------------
  interface
    subroutine xc_f90_lca_init(p, info, functional, nspin)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(out) :: p
      type(xc_f90_info_t), intent(out) :: info
      integer, intent(in) :: functional
      integer, intent(in) :: nspin
    end subroutine xc_f90_lca_init

    subroutine xc_f90_lca_end(p)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(inout) :: p
    end subroutine xc_f90_lca_end
  end interface


  !----------------------------------------------------------------
  interface
    subroutine xc_f90_lca(p, rho, v, e, dedd, dedv)
      use xc_f90_types_m
      type(xc_f90_func_t), intent(in) :: p
      real(xc_f90_kind), intent(in) :: rho ! rho(nspin) the density
      real(xc_f90_kind), intent(in) :: v ! v(3,nspin) the vorticity
      real(xc_f90_kind), intent(out) :: e ! the energy per unit particle
      real(xc_f90_kind), intent(out) :: dedd ! dedd(nspin) the derivative of the energy
                                              ! in terms of the density
      real(xc_f90_kind), intent(out) :: dedv ! and in terms of the vorticity
    end subroutine xc_f90_lca
  end interface

end module xc_f90_lib_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
