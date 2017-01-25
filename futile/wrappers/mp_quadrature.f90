!> @file
!! Multipole preserving quadrature scheme
!! @author
!!    Copyright (C) 2016-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module multipole_preserving
  use dynamic_memory
  use f_precisions, only: gp => f_double
  implicit none

  private

  integer :: itype_scf=0                          !< Type of the interpolating SCF, 0= data unallocated
  integer :: n_scf=-1                             !< Number of points of the allocated data
  integer :: nrange_scf=0                         !< range of the integration
  real(gp), dimension(:), allocatable :: scf_data !< Values for the interpolating scaling functions points
  !> log of the minimum value of the scf data
  !! to avoid floating point exceptions while multiplying with it
  real(gp) :: mn_scf = 0.0_gp

  public :: initialize_real_space_conversion,finalize_real_space_conversion,scfdotf,mp_exp

  contains

    !> Prepare the array for the evaluation with the interpolating Scaling Functions
    !! one might add also the function to be converted and the 
    !! prescription for integrating knowing the scaling relation of the function
    subroutine initialize_real_space_conversion(npoints,isf_m,nmoms)
      implicit none
      integer, intent(in), optional :: npoints,isf_m,nmoms
      !local variables
      character(len=*), parameter :: subname='initialize_real_space_conversion'
      integer :: n_range,i,nmm
      real(gp) :: tt
      real(gp), dimension(:), allocatable :: x_scf !< to be removed in a future implementation

      if (present(isf_m)) then
         itype_scf=isf_m
      else
         itype_scf=16
      end if

      if (present(nmoms)) then
         nmm=nmoms
      else
         nmm=0
      end if

      if (present(npoints)) then
         n_scf=2*(itype_scf+nmm)*npoints
      else
         n_scf=2*(itype_scf+nmm)*(2**6)
      end if

      !allocations for scaling function data array
      x_scf = f_malloc(0.to.n_scf,id='x_scf')

      scf_data = f_malloc(0.to.n_scf,id='scf_data')

      !Build the scaling function external routine coming from Poisson Solver. To be customized accordingly
      !call scaling_function(itype_scf,n_scf,n_range,x_scf,scf_data)
      !call wavelet_function(itype_scf,n_scf,x_scf,scf_data)
      call ISF_family(itype_scf,nmm,n_scf,n_range,x_scf,scf_data)
      !stop 
      call f_free(x_scf)

      nrange_scf=n_range
      !define the log of the smallest nonzero value as the 
      !cutoff for multiplying with it
      !this means that the values which are 
      !lower than scf_data squared will be considered as zero
      mn_scf=epsilon(1.d0)**2 !just to put a "big" value
      do i=0,n_scf
         tt=scf_data(i)
         if (tt /= 0.0_gp .and. abs(tt) > sqrt(tiny(1.0))) then
            mn_scf=min(mn_scf,tt**2)
         end if
      end do

    end subroutine initialize_real_space_conversion


    !> Deallocate scf_data
    subroutine finalize_real_space_conversion()
      implicit none

      itype_scf=0
      n_scf=-1
      mn_scf=0.0_gp
      nrange_scf=0
      call f_free(scf_data)

    end subroutine finalize_real_space_conversion


    !> multipole-preserving gaussian function
    !! chooses between traditional exponential and scfdotf 
    !! according to the value of the exponent in units of the grid spacing
    !! the function is supposed to be x**pow*exp(-expo*x**2)
    !! where x=hgrid*j-x0
    !! @warning
    !! this function is also elemental to ease its evaluation, though 
    !! the usage for vector argument is discouraged: dedicated routines has to be 
    !! written to meet performance
    !! @todo 
    !!  Optimize it!
    elemental pure function mp_exp(hgrid,x0,expo,j,pow,modified)
      use numerics, only: safe_exp
      implicit none
      real(gp), intent(in) :: hgrid   !< Hgrid 
      real(gp), intent(in) :: x0      !< X value
      real(gp), intent(in) :: expo    !< Exponent of the gaussian
      logical, intent(in) :: modified !< Switch to scfdotf if true
      integer, intent(in) :: j        !< Location of the scf from x0
      integer, intent(in) :: pow      !< Exp(-expo*x**2)*(x**pow)
      real(gp) :: mp_exp
      !local variables
      real(gp) :: x

      !added failsafe to avoid segfaults
      if (modified .and. allocated(scf_data)) then
         mp_exp=scfdotf(j,hgrid,expo,x0,pow)
      else
         x=hgrid*j-x0
         mp_exp=safe_exp(-expo*x**2)
         if (pow /= 0) mp_exp=mp_exp*(x**pow)
      end if
    end function mp_exp


    !> This function calculates the scalar product between a ISF and a 
    !! input function, which is a gaussian times a power centered
    !! @f$g(x) = (x-x_0)^{pow} e^{-pgauss (x-x_0)}@f$
    !! here pure specifier is redundant
    !! we should add here the threshold from which the 
    !! normal function can be evaluated
    elemental pure function scfdotf(j,hgrid,pgauss,x0,pow) result(gint)
      use numerics, only: safe_exp
      implicit none
      !Arguments
      integer, intent(in) :: j !<value of the input result in the hgrid reference
      integer, intent(in) :: pow
      real(gp), intent(in) :: hgrid,pgauss,x0
      real(gp) :: gint
      !local variables
      integer :: i
      real(gp) :: x,absci,fabsci,dx
      gint=0.0_gp

      !Step grid for the integration
      !dx = real(2*itype_scf,gp)/real(n_scf,gp)
      dx = real(nrange_scf,gp)/real(n_scf,gp)
      !starting point for the x coordinate for integration
      !x  = real(j-itype_scf+1,gp)-dx
      x  = real(j-nrange_scf/2+1,gp)-dx

      !the loop can be unrolled to maximize performances
      if (pow /= 0) then
         do i=0,n_scf
            x=x+dx
            absci = x*hgrid - x0
            !here evaluate the function
            fabsci = absci**pow
            absci = -pgauss*absci*absci
            fabsci = fabsci*safe_exp(absci,underflow=mn_scf)
            !calculate the integral
            gint = gint + scf_data(i)*fabsci
            !       print *,'test',i,scf_data(i),fabsci,pgauss,pow,absci
         end do
      else
         do i=0,n_scf
            x=x+dx
            absci = x*hgrid - x0
            !          !here evaluate the function
            absci = -pgauss*absci*absci
            fabsci = safe_exp(absci,underflow=mn_scf)
            !calculate the integral
            !          fabsci= safe_gaussian(x0,x*hgrid,pgauss)
            !print *,'test',i,scf_data(i),fabsci,pgauss,absci,log(tiny(1.d0)),tiny(1.d0)
            gint = gint + scf_data(i)*fabsci

         end do
      end if
      gint = gint*dx

    end function scfdotf
    
end module multipole_preserving
