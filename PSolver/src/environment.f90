!> @file
!!    Modulefile for environment computation
!!
!! @author
!!    G. Fisicaro, L. Genovese (September 2015)
!!    Copyright (C) 2002-2015 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
module environment
  use f_enums, only: f_enumerator
  use f_precisions, only: f_double
  use numerics, only: safe_exp
  use dictionaries, only: f_err_throw
  implicit none

  private

  ! General precision, density and the potential types, to be moved in a low-levle module
  integer, parameter, public :: gp=f_double  !< general-type precision
  integer, parameter, public :: dp=f_double  !< density-type precision
  integer, parameter, public :: wp=f_double  !< potential-type precision
  real(gp), parameter :: pi=3.141592653589793238462643383279502884197_gp
  real(gp), parameter :: vacuum_eps=1.0_gp


  !> how to set the dielectric function
  integer, parameter, public :: PS_EPSILON_VACUUM = -1000
  integer, parameter, public :: PS_EPSILON_RIGID_CAVITY = 1001
  integer, parameter, public :: PS_EPSILON_SCCS = 1002

  integer, parameter :: PS_PCG = 1234
  integer, parameter :: PS_PI = 1432

  type(f_enumerator), public :: PS_NONE_ENUM=f_enumerator('vacuum',PS_EPSILON_VACUUM,null())
  type(f_enumerator), public :: PS_RIGID_ENUM=f_enumerator('rigid',PS_EPSILON_RIGID_CAVITY,null())
  type(f_enumerator), public :: PS_SCCS_ENUM=f_enumerator('sccs',PS_EPSILON_SCCS,null())

  type(f_enumerator), parameter, public :: PS_VAC_ENUM=f_enumerator('VAC',PS_EPSILON_VACUUM,null())
  type(f_enumerator), parameter, public :: PS_PI_ENUM=f_enumerator('PI',PS_PI,null())
  type(f_enumerator), parameter, public :: PS_PCG_ENUM=f_enumerator('PCG',PS_PCG,null())
  type(f_enumerator), parameter, public :: PS_PB_ENUM=f_enumerator('PB',PS_PCG,null()) !< poisson boltzmann

  
  !> define the cavity type
  type, public :: cavity_data
     real(gp) :: epsilon0 !< dielectriv constant of the medium
     real(gp) :: edensmax !<maximum value of the density for the cavity
     real(gp) :: edensmin !<minimum  value of the density for the cavity
     real(dp) :: gammaS !< surface tension of the solvent [dyn/cm]
     real(dp) :: alphaS !< proportionality factor for the repulsion free energy in term of the surface integral [dyn/cm]
     real(dp) :: betaV !<proportionality factor for the dispersion free energy in term of the volume integral [GPa]     
  end type cavity_data

  public :: cavity_init,eps,epsprime,epssecond,oneoeps,oneosqrteps,logepsprime,corr_term,nabla_u_square
  public :: nabla_u,div_u_i,nabla_u_and_square,cavity_default,nonvacuum_projection,update_rhopol

contains

  pure function cavity_default() result(c)
    implicit none
    type(cavity_data) :: c
    c%epsilon0= 78.36_gp !<water at ambient condition 
    c%edensmax = 0.005_gp !0.0050d0
    c%edensmin = 0.0001_gp
    c%gammaS = 72._gp ![dyn/cm]   
    c%alphaS = -22.0_gp ![dyn/cm]   end function cavity_default
    c%betaV = -0.35_gp ![GPa]     
  end function cavity_default

  !>initialize the cavity parameters
  function cavity_init(epsilon0,edensmax,edensmin,gammaS,alphaS,betaV) result(c)
    implicit none
    real(gp), intent(in), optional :: epsilon0,edensmax,edensmin,gammaS,alphaS,betaV
    type(cavity_data) :: c
    c=cavity_default()
    if (present(epsilon0)) c%epsilon0=epsilon0
    if (present(edensmax)) c%edensmax=edensmax
    if (present(edensmin)) c%edensmax=edensmin
    if (present(gammaS)) c%gammaS=gammaS
    if (present(alphaS)) c%alphaS=alphaS  
    if (present(betaV )) c%betaV =betaV 
  end function cavity_init

  pure function epsilon_transition(rho,pow,der,cavity) result(epsilon)
    character(len=*), intent(in) :: pow !<power to epsilon
    integer, intent(in) :: der !< derivative of epsolin
    real(dp), intent(in) :: rho
    type(cavity_data), intent(in) :: cavity
    real(gp) :: epsilon
    !local variables
    real(gp) :: fact1,fact2,fact0,r,zeta,w,dw,dzetao2pi,d2zetao2pi,dzeta2o2pi,d2w

    epsilon=-1.d0 !this value should always be overwritten
    fact0=cavity%edensmax/cavity%edensmin
    r=cavity%edensmax/rho
    fact1=1.0_gp/log(fact0)
    fact2=log(r)
    zeta=2.0_gp*pi*fact1*fact2
    w=(zeta-sin(zeta))/(2.d0*pi)
    dzetao2pi=-fact1/rho
    dw=(1.d0-cos(zeta))*dzetao2pi

    select case(pow)
    case('1')
       select case(der)
       case(0)
          epsilon=cavity%epsilon0*safe_exp(w)
       case(1)
          epsilon=cavity%epsilon0*safe_exp(w)*dw
       case(2)
          d2zetao2pi=fact1/rho**2
          dzeta2o2pi=2.0_gp*pi*dzetao2pi**2
          d2w=sin(zeta)*dzeta2o2pi+(1.d0-cos(zeta))*d2zetao2pi
          epsilon=cavity%epsilon0*safe_exp(w)*(dw*dw+d2w)
       end select
    case('-1/2')
       epsilon=safe_exp(-0.5_gp*w)/sqrt(cavity%epsilon0)
    case('-1')
       epsilon=safe_exp(-w)/cavity%epsilon0
    case('L')
       select case(der)
       case(1)
          epsilon=dw
       end select
    end select

  end function epsilon_transition

  pure function eps(rho,cavity)
    implicit none
    real(dp), intent(in) :: rho
    type(cavity_data), intent(in) :: cavity
    real(gp) :: eps

    !we are in a inner region
    if (rho > cavity%edensmax) then
       eps=vacuum_eps
    !we are in a outer region
    else if (rho < cavity%edensmin) then
       eps=cavity%epsilon0
    else
       eps=epsilon_transition(rho,'1',0,cavity)
    end if
  end function eps

  pure function oneosqrteps(rho,cavity) result(eps)
    implicit none
    real(dp), intent(in) :: rho
    type(cavity_data), intent(in) :: cavity
    real(gp) :: eps

    !we are in a inner region
    if (rho > cavity%edensmax) then
       eps=1.0_gp/sqrt(vacuum_eps)
       !we are in a outer region
    else if (rho < cavity%edensmin) then
       eps=1.0_gp/sqrt(cavity%epsilon0)
    else
       eps=epsilon_transition(rho,'-1/2',0,cavity)
    end if
  end function oneosqrteps

  pure function oneoeps(rho,cavity) result(eps)
    implicit none
    real(dp), intent(in) :: rho
    type(cavity_data), intent(in) :: cavity
    real(gp) :: eps

    !we are in a inner region
    if (rho > cavity%edensmax) then
       eps=1.0_gp/vacuum_eps
       !we are in a outer region
    else if (rho < cavity%edensmin) then
       eps=1.0_gp/cavity%epsilon0
    else
       eps=epsilon_transition(rho,'-1',0,cavity)
    end if
  end function oneoeps

  pure function epsprime(rho,cavity) result(eps)
    implicit none
    real(dp), intent(in) :: rho
    type(cavity_data), intent(in) :: cavity
    real(gp) :: eps

    !we are in a inner region
    if (rho > cavity%edensmax) then
       eps=0.0_gp
       !we are in a outer region
    else if (rho < cavity%edensmin) then
       eps=0.0_gp
    else
       eps=epsilon_transition(rho,'1',1,cavity)
    end if
  end function epsprime

  pure function epssecond(rho,cavity) result(eps)
    implicit none
    real(dp), intent(in) :: rho
    type(cavity_data), intent(in) :: cavity
    real(gp) :: eps

    !we are in a inner region
    if (rho > cavity%edensmax) then
       eps=0.0_gp
       !we are in a outer region
    else if (rho < cavity%edensmin) then
       eps=0.0_gp
    else
       eps=epsilon_transition(rho,'1',2,cavity)
    end if
  end function epssecond

  pure function logepsprime(rho,cavity) result(eps)
    implicit none
    real(dp), intent(in) :: rho
    type(cavity_data), intent(in) :: cavity
    real(gp) :: eps

    !we are in a inner region
    if (rho > cavity%edensmax) then
       eps=0.0_gp
       !we are in a outer region
    else if (rho < cavity%edensmin) then
       eps=0.0_gp
    else
       eps=epsilon_transition(rho,'L',1,cavity)
    end if
  end function logepsprime

  pure function corr_term(rho,nabla2rho,deltarho,cavity)
    implicit none
    real(gp), intent(in) :: rho !<density
    real(gp), intent(in) :: nabla2rho !<square of the density gradient
    real(gp), intent(in) :: deltarho !<square of the density gradient
    type(cavity_data), intent(in) :: cavity
    real(gp) :: corr_term
    !local variables
    real(gp) :: epspr

    epspr=epsprime(rho,cavity)
    corr_term=-0.125_gp/pi*(0.5_gp*nabla2rho*&
         (epspr*logepsprime(rho,cavity)-epssecond(rho,cavity))-epspr*deltarho)
  end function corr_term
  


  !>This routine computes 'nord' order accurate first derivatives 
  !! on a equally spaced grid with coefficients from 'Mathematica' program.
  !! 
  !! input:
  !! ngrid       = number of points in the grid, 
  !! u(ngrid)    = function values at the grid points
  !! 
  !! output:
  !! du(ngrid)   = first derivative values at the grid points
  subroutine nabla_u_and_square(geocode,n01,n02,n03,u,du,du2,nord,hgrids)
    implicit none


    !c..declare the pass
    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
    integer, intent(in) :: n01,n02,n03,nord
    real(kind=8), dimension(3), intent(in) :: hgrids
    real(kind=8), dimension(n01,n02,n03), intent(in) :: u
    real(kind=8), dimension(n01,n02,n03,3), intent(out) :: du !<nabla of u
    real(kind=8), dimension(n01,n02,n03) :: du2 !<square module of nabla u

    !c..local variables
    integer :: n,m,n_cell
    integer :: i,j,ib,i1,i2,i3,ii
    real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D, c1D_1, c1D_2, c1D_3
    real(kind=8) :: hx,hy,hz,d
    logical :: perx,pery,perz

    n = nord+1
    m = nord/2
    hx = hgrids(1)!acell/real(n01,kind=8)
    hy = hgrids(2)!acell/real(n02,kind=8)
    hz = hgrids(3)!acell/real(n03,kind=8)
    n_cell = max(n01,n02,n03)

    !buffers associated to the geocode
    !conditions for periodicity in the three directions
    perx=(geocode /= 'F')
    pery=(geocode == 'P')
    perz=(geocode /= 'F')

    ! Beware that n_cell has to be > than n.
    if (n_cell.lt.n) then
       call f_err_throw('Ngrid in has to be setted > than n=nord + 1')
       !stop
    end if

    ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
    !Only nord=2,4,6,8,16

    select case(nord)
    case(2,4,6,8,16)
       !O.K.
    case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
    end select

    do i=-m,m
       do j=-m,m
          c1D(i,j)=0.d0
       end do
    end do

    include 'FiniteDiffCorff.inc'
    c1D_1 = c1D/hx
    c1D_2 = c1D/hy
    c1D_3 = c1D/hz
!!!default(shared) private(i1,i2,i3,j,ii, d) 
    !$omp parallel do default(none) &
    !$omp private(i3,i2,i1,d,ii,j) &
    !$omp shared(du2,perx,m,n01,n02,n03,du,c1D_1,u)
    do i3=1,n03
       do i2=1,n02
          do i1=1,n01

             d = 0.0d0
             !du2(i1,i2,i3) = 0.0d0

             if (i1.le.m) then
                if (perx) then
                   do j=-m,m
                      ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                      d = d + c1D_1(j,0)*u(ii,i2,i3)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3)
                   end do
                end if
             else if (i1.gt.n01-m) then
                if (perx) then
                   do j=-m,m
                      ii=modulo(i1 + j - 1, n01 ) + 1
                      d = d + c1D_1(j,0)*u(ii,i2,i3)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3)
                   end do
                end if
             else
                do j=-m,m
                   d = d + c1D_1(j,0)*u(i1 + j,i2,i3)
                end do
             end if
             du(i1,i2,i3,1)= d
             du2(i1,i2,i3) = d*d

          end do
       end do
    end do
    !$omp end parallel do

    !$omp parallel do default(shared) private(i1,i2,i3,j,ii,d) 
    do i3=1,n03
       do i2=1,n02
          do i1=1,n01

             d = 0.0d0

             if (i2.le.m) then
                if (pery) then
                   do j=-m,m
                      ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                      d = d + c1D_2(j,0)*u(i1,ii,i3)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_2(j,i2-m-1)*u(i1,j+m+1,i3)
                   end do
                end if
             else if (i2.gt.n02-m) then
                if (pery) then
                   do j=-m,m
                      ii=modulo(i2 + j - 1, n02 ) + 1
                      d = d + c1D_2(j,0)*u(i1,ii,i3)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_2(j,i2-n02+m)*u(i1,n02 + j - m,i3)
                   end do
                end if
             else
                do j=-m,m
                   d = d + c1D_2(j,0)*u(i1,i2 + j,i3)
                end do
             end if
             du(i1,i2,i3,2)= d
             du2(i1,i2,i3) = du2(i1,i2,i3) + d*d
          end do
       end do
    end do
    !$omp end parallel do

    !$omp parallel do default(shared) private(i1,i2,i3,j,ii,d) 
    do i3=1,n03
       do i2=1,n02
          do i1=1,n01
             d = 0.0d0

             if (i3.le.m) then
                if (perz) then
                   do j=-m,m
                      ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                      d = d + c1D_3(j,0)*u(i1,i2,ii)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_3(j,i3-m-1)*u(i1,i2,j+m+1)
                   end do
                end if
             else if (i3.gt.n03-m) then
                if (perz) then
                   do j=-m,m
                      ii=modulo(i3 + j - 1, n03 ) + 1
                      d = d + c1D_3(j,0)*u(i1,i2,ii)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_3(j,i3-n03+m)*u(i1,i2,n03 + j - m)
                   end do
                end if
             else
                do j=-m,m
                   d = d + c1D_3(j,0)*u(i1,i2,i3 + j)
                end do
             end if
             du(i1,i2,i3,3)=d
             du2(i1,i2,i3) = du2(i1,i2,i3) + d*d

          end do
       end do
    end do
    !$omp end parallel do
  end subroutine nabla_u_and_square

  !>only calculate the square of the gradient (in orthorhombic cells)
  subroutine nabla_u_square(geocode,n01,n02,n03,u,du2,nord,hgrids)
    implicit none


    !c..declare the pass
    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
    integer, intent(in) :: n01,n02,n03,nord
    real(kind=8), dimension(3), intent(in) :: hgrids
    real(kind=8), dimension(n01,n02,n03), intent(in) :: u
    real(kind=8), dimension(n01,n02,n03), intent(out) :: du2 !<square module of nabla u

    !c..local variables
    integer :: n,m,n_cell
    integer :: i,j,ib,i1,i2,i3,ii
    real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D, c1D_1, c1D_2, c1D_3
    real(kind=8) :: hx,hy,hz,d
    logical :: perx,pery,perz

    n = nord+1
    m = nord/2
    hx = hgrids(1)!acell/real(n01,kind=8)
    hy = hgrids(2)!acell/real(n02,kind=8)
    hz = hgrids(3)!acell/real(n03,kind=8)
    n_cell = max(n01,n02,n03)

    !buffers associated to the geocode
    !conditions for periodicity in the three directions
    perx=(geocode /= 'F')
    pery=(geocode == 'P')
    perz=(geocode /= 'F')

    ! Beware that n_cell has to be > than n.
    if (n_cell.lt.n) then
       call f_err_throw('Ngrid in has to be setted > than n=nord + 1')
       !stop
    end if

    ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
    !Only nord=2,4,6,8,16

    select case(nord)
    case(2,4,6,8,16)
       !O.K.
    case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
    end select

    do i=-m,m
       do j=-m,m
          c1D(i,j)=0.d0
       end do
    end do

    include 'FiniteDiffCorff.inc'
    c1D_1 = c1D/hx
    c1D_2 = c1D/hy
    c1D_3 = c1D/hz
!!!default(shared) private(i1,i2,i3,j,ii, d) 
    !$omp parallel do default(none) &
    !$omp private(i3,i2,i1,d,ii,j) &
    !$omp shared(du2,perx,m,n01,n02,n03,c1D_1,u)
    do i3=1,n03
       do i2=1,n02
          do i1=1,n01

             d = 0.0d0
             !du2(i1,i2,i3) = 0.0d0

             if (i1.le.m) then
                if (perx) then
                   do j=-m,m
                      ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                      d = d + c1D_1(j,0)*u(ii,i2,i3)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3)
                   end do
                end if
             else if (i1.gt.n01-m) then
                if (perx) then
                   do j=-m,m
                      ii=modulo(i1 + j - 1, n01 ) + 1
                      d = d + c1D_1(j,0)*u(ii,i2,i3)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3)
                   end do
                end if
             else
                do j=-m,m
                   d = d + c1D_1(j,0)*u(i1 + j,i2,i3)
                end do
             end if
             !du(i1,i2,i3,1)= d
             du2(i1,i2,i3) = d*d

          end do
       end do
    end do
    !$omp end parallel do

    !$omp parallel do default(shared) private(i1,i2,i3,j,ii,d) 
    do i3=1,n03
       do i2=1,n02
          do i1=1,n01

             d = 0.0d0

             if (i2.le.m) then
                if (pery) then
                   do j=-m,m
                      ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                      d = d + c1D_2(j,0)*u(i1,ii,i3)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_2(j,i2-m-1)*u(i1,j+m+1,i3)
                   end do
                end if
             else if (i2.gt.n02-m) then
                if (pery) then
                   do j=-m,m
                      ii=modulo(i2 + j - 1, n02 ) + 1
                      d = d + c1D_2(j,0)*u(i1,ii,i3)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_2(j,i2-n02+m)*u(i1,n02 + j - m,i3)
                   end do
                end if
             else
                do j=-m,m
                   d = d + c1D_2(j,0)*u(i1,i2 + j,i3)
                end do
             end if
             !du(i1,i2,i3,2)= d
             du2(i1,i2,i3) = du2(i1,i2,i3) + d*d
          end do
       end do
    end do
    !$omp end parallel do

    !$omp parallel do default(shared) private(i1,i2,i3,j,ii,d) 
    do i3=1,n03
       do i2=1,n02
          do i1=1,n01
             d = 0.0d0

             if (i3.le.m) then
                if (perz) then
                   do j=-m,m
                      ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                      d = d + c1D_3(j,0)*u(i1,i2,ii)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_3(j,i3-m-1)*u(i1,i2,j+m+1)
                   end do
                end if
             else if (i3.gt.n03-m) then
                if (perz) then
                   do j=-m,m
                      ii=modulo(i3 + j - 1, n03 ) + 1
                      d = d + c1D_3(j,0)*u(i1,i2,ii)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_3(j,i3-n03+m)*u(i1,i2,n03 + j - m)
                   end do
                end if
             else
                do j=-m,m
                   d = d + c1D_3(j,0)*u(i1,i2,i3 + j)
                end do
             end if
             !du(i1,i2,i3,3)=d
             du2(i1,i2,i3) = du2(i1,i2,i3) + d*d

          end do
       end do
    end do
    !$omp end parallel do
  end subroutine nabla_u_square

  !>this routine computes 'nord' order accurate first derivatives 
  !!on a equally spaced grid with coefficients from 'Matematica' program.
  !!
  !!input:
  !!ngrid       = number of points in the grid, 
  !!u(ngrid)    = function values at the grid points
  !!
  !!output:
  !!du(ngrid)   = first derivative values at the grid points
  !!
  !!declare the pass
  subroutine div_u_i(geocode,n01,n02,n03,u,du,nord,hgrids,cc)
    implicit none

    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
    integer, intent(in) :: n01,n02,n03,nord
    real(kind=8), dimension(3), intent(in) :: hgrids
    real(kind=8), dimension(n01,n02,n03,3), intent(in) :: u !<vector field u_i
    real(kind=8), dimension(n01,n02,n03), intent(out) :: du !<divergence d_i u_i
    real(kind=8), dimension(n01,n02,n03), intent(out), optional :: cc !< u_i u_j d_i u_j , needed for the surface term where u_i=d_i rho

    !c..local variables
    integer :: n,m,n_cell
    integer :: i,j,ib,i1,i2,i3,ii
    real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D, c1D_1, c1D_2, c1D_3
    real(kind=8) :: hx,hy,hz,d1,uxy,uyz,uxz
    real(kind=8), parameter :: zero = 0.d0! 1.0d-11
    logical :: perx,pery,perz

    n = nord+1
    m = nord/2
    hx = hgrids(1)!acell/real(n01,kind=8)
    hy = hgrids(2)!acell/real(n02,kind=8)
    hz = hgrids(3)!acell/real(n03,kind=8)
    n_cell = max(n01,n02,n03)

    !buffers associated to the geocode
    !conditions for periodicity in the three directions
    perx=(geocode /= 'F')
    pery=(geocode == 'P')
    perz=(geocode /= 'F')

    ! Beware that n_cell has to be > than n.
    if (n_cell < n) then
       write(*,*)'ngrid in has to be set > than n=nord + 1'
       stop
    end if

    ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
    !Only nord=2,4,6,8,16

    select case(nord)
    case(2,4,6,8,16)
       !O.K.
    case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
    end select

    do i=-m,m
       do j=-m,m
          c1D(i,j)=0.d0
       end do
    end do

    include 'FiniteDiffCorff.inc'

    c1D_1 = c1D/hx
    c1D_2 = c1D/hy
    c1D_3 = c1D/hz

!!$  if (present(cc)) then
!!$     cc(i1,i2,i3) = (u(i1,i2,i3,1)**2)*d1+(u(i1,i2,i3,2)**2)*d2+(u(i1,i2,i3,3)**2)*d3+&
!!$          2.d0*u(i1,i2,i3,1)*u(i1,i2,i3,2)*uxy+2.d0*u(i1,i2,i3,2)*u(i1,i2,i3,3)*uyz+&
!!$          2.d0*u(i1,i2,i3,1)*u(i1,i2,i3,3)*uxz
!!$  end if


    if (present(cc)) then
       !$omp parallel do default(none) &
       !$omp private(i1,i2,i3,j,ii, d1,uxz,uxy) &
       !$omp shared(du,c1D_1,u,perx,m,hx,n01,n02,n03,cc)
       do i3=1,n03
          do i2=1,n02
             do i1=1,n01

                du(i1,i2,i3) = 0.0d0

                d1 = 0.d0
                uxy = 0.d0
                uxz = 0.d0
                if (i1.le.m) then
                   if (perx) then
                      do j=-m,m
                         ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                         d1 = d1 + c1D_1(j,0)*u(ii,i2,i3,1)
                         uxy = uxy + c1D_1(j,0)*u(ii,i2,i3,2)!/hx
                         uxz = uxz + c1D_1(j,0)*u(ii,i2,i3,3)!/hx
                      end do
                   else
                      do j=-m,m
                         d1 = d1 + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3,1)
                         uxy = uxy + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3,2)!/hx
                         uxz = uxz + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3,3)!/hx
                      end do
                   end if
                else if (i1.gt.n01-m) then
                   if (perx) then
                      do j=-m,m
                         ii=modulo(i1 + j - 1, n01 ) + 1
                         d1 = d1 + c1D_1(j,0)*u(ii,i2,i3,1)
                         uxy = uxy + c1D_1(j,0)*u(ii,i2,i3,2)!/hx
                         uxz = uxz + c1D_1(j,0)*u(ii,i2,i3,3)!/hx
                      end do
                   else
                      do j=-m,m
                         d1 = d1 + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3,1)
                         uxy = uxy + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3,2)!/hx
                         uxz = uxz + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3,3)!/hx
                      end do
                   end if
                else
                   do j=-m,m
                      d1 = d1 + c1D_1(j,0)*u(i1 + j,i2,i3,1)
                      uxy = uxy + c1D_1(j,0)*u(i1 + j,i2,i3,2)!/hx
                      uxz = uxz + c1D_1(j,0)*u(i1 + j,i2,i3,3)!/hx
                   end do
                end if
                !uxy=uxy/hx
                !uxz=uxz/hx

                du(i1,i2,i3) =d1
                cc(i1,i2,i3) = (u(i1,i2,i3,1)**2)*d1 + &
                     2.d0*u(i1,i2,i3,1)*u(i1,i2,i3,2)*uxy + &
                     2.d0*u(i1,i2,i3,1)*u(i1,i2,i3,3)*uxz

             end do
          end do
       end do
       !$omp end parallel do

       !shared) 
       !$omp parallel do default(none) &
       !$omp private(i1,i2,i3,j,ii,d1,uyz) &
       !$omp shared(n01,n02,n03,pery,m,c1D_2,u,du,cc)
       do i3=1,n03
          do i2=1,n02
             do i1=1,n01
                d1=0.d0
                uyz = 0.d0

                if (i2.le.m) then
                   if (pery) then
                      do j=-m,m
                         ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                         d1 = d1 + c1D_2(j,0)*u(i1,ii,i3,2)
                         uyz = uyz + c1D_2(j,0)*u(i1,ii,i3,3)!/hy
                      end do
                   else
                      do j=-m,m
                         d1 = d1 + c1D_2(j,i2-m-1)*u(i1,j+m+1,i3,2)
                         uyz = uyz + c1D_2(j,i2-m-1)*u(i1,j+m+1,i3,3)!/hy
                      end do
                   end if
                else if (i2.gt.n02-m) then
                   if (pery) then
                      do j=-m,m
                         ii=modulo(i2 + j - 1, n02 ) + 1
                         d1 = d1 + c1D_2(j,0)*u(i1,ii,i3,2)
                         uyz = uyz + c1D_2(j,0)*u(i1,ii,i3,3)!/hy
                      end do
                   else
                      do j=-m,m
                         d1 = d1 + c1D_2(j,i2-n02+m)*u(i1,n02 + j - m,i3,2)
                         uyz = uyz + c1D_2(j,i2-n02+m)*u(i1,n02 + j - m,i3,3)!/hy
                      end do
                   end if
                else
                   do j=-m,m
                      d1 = d1 + c1D_2(j,0)*u(i1,i2 + j,i3,2)
                      uyz = uyz + c1D_2(j,0)*u(i1,i2 + j,i3,3)!/hy
                   end do
                end if
                du(i1,i2,i3) = du(i1,i2,i3) + d1
                !uyz=uyz/hy
                cc(i1,i2,i3) = cc(i1,i2,i3) + (u(i1,i2,i3,2)**2)*d1+ &
                     2.d0*u(i1,i2,i3,2)*u(i1,i2,i3,3)*uyz
             end do
          end do
       end do
       !$omp end parallel do


       !(shared) private(i1,i2,i3,j,ii, d1) 
       !$omp parallel do default(none) &
       !$omp private(i1,i2,i3,j,ii,d1) &
       !$omp shared(n01,n02,n03,perz,m,c1D_3,u,du,cc)
       do i3=1,n03
          do i2=1,n02
             do i1=1,n01
                d1=0.d0
                if (i3.le.m) then
                   if (perz) then
                      do j=-m,m
                         ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                         d1 = d1 + c1D_3(j,0)*u(i1,i2,ii,3)
                      end do
                   else
                      do j=-m,m
                         d1 = d1 + c1D_3(j,i3-m-1)*u(i1,i2,j+m+1,3)
                      end do
                   end if
                else if (i3.gt.n03-m) then
                   if (perz) then
                      do j=-m,m
                         ii=modulo(i3 + j - 1, n03 ) + 1
                         d1 = d1 + c1D_3(j,0)*u(i1,i2,ii,3)
                      end do
                   else
                      do j=-m,m
                         d1 = d1 + c1D_3(j,i3-n03+m)*u(i1,i2,n03 + j - m,3)
                      end do
                   end if
                else
                   do j=-m,m
                      d1 = d1 + c1D_3(j,0)*u(i1,i2,i3 + j,3)
                   end do
                end if
                du(i1,i2,i3) = du(i1,i2,i3)+d1
                cc(i1,i2,i3) = cc(i1,i2,i3) + (u(i1,i2,i3,3)**2)*d1
             end do
          end do
       end do
       !$omp end parallel do
    else

       !$omp parallel do default(none) &
       !$omp private(i1,i2,i3,j,ii, d1) &
       !$omp shared(du,c1D_1,u,perx,m,hx,n01,n02,n03)
       do i3=1,n03
          do i2=1,n02
             do i1=1,n01

                du(i1,i2,i3) = 0.0d0

                d1 = 0.d0
                if (i1.le.m) then
                   if (perx) then
                      do j=-m,m
                         ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                         d1 = d1 + c1D_1(j,0)*u(ii,i2,i3,1)
                      end do
                   else
                      do j=-m,m
                         d1 = d1 + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3,1)
                      end do
                   end if
                else if (i1.gt.n01-m) then
                   if (perx) then
                      do j=-m,m
                         ii=modulo(i1 + j - 1, n01 ) + 1
                         d1 = d1 + c1D_1(j,0)*u(ii,i2,i3,1)
                      end do
                   else
                      do j=-m,m
                         d1 = d1 + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3,1)
                      end do
                   end if
                else
                   do j=-m,m
                      d1 = d1 + c1D_1(j,0)*u(i1 + j,i2,i3,1)
                   end do
                end if
                du(i1,i2,i3) =d1
             end do
          end do
       end do
       !$omp end parallel do

       !default(shared) private(i1,i2,i3,j,ii,d1,uyz) 
       !$omp parallel do default(none) &
       !$omp private(i1,i2,i3,j,ii,d1) &
       !$omp shared(c1D_2,u,n01,n02,n03,du,m,pery)
       do i3=1,n03
          do i2=1,n02
             do i1=1,n01
                d1=0.d0
                if (i2.le.m) then
                   if (pery) then
                      do j=-m,m
                         ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                         d1 = d1 + c1D_2(j,0)*u(i1,ii,i3,2)
                      end do
                   else
                      do j=-m,m
                         d1 = d1 + c1D_2(j,i2-m-1)*u(i1,j+m+1,i3,2)
                      end do
                   end if
                else if (i2.gt.n02-m) then
                   if (pery) then
                      do j=-m,m
                         ii=modulo(i2 + j - 1, n02 ) + 1
                         d1 = d1 + c1D_2(j,0)*u(i1,ii,i3,2)
                      end do
                   else
                      do j=-m,m
                         d1 = d1 + c1D_2(j,i2-n02+m)*u(i1,n02 + j - m,i3,2)
                      end do
                   end if
                else
                   do j=-m,m
                      d1 = d1 + c1D_2(j,0)*u(i1,i2 + j,i3,2)
                   end do
                end if
                du(i1,i2,i3) = du(i1,i2,i3) + d1
             end do
          end do
       end do
       !$omp end parallel do

       !$omp parallel do default(none) &
       !$omp private(i1,i2,i3,j,ii,d1) &
       !$omp shared(c1D_3,u,n01,n02,n03,du,m,perz)
       do i3=1,n03
          do i2=1,n02
             do i1=1,n01
                d1=0.d0
                if (i3.le.m) then
                   if (perz) then
                      do j=-m,m
                         ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                         d1 = d1 + c1D_3(j,0)*u(i1,i2,ii,3)
                      end do
                   else
                      do j=-m,m
                         d1 = d1 + c1D_3(j,i3-m-1)*u(i1,i2,j+m+1,3)
                      end do
                   end if
                else if (i3.gt.n03-m) then
                   if (perz) then
                      do j=-m,m
                         ii=modulo(i3 + j - 1, n03 ) + 1
                         d1 = d1 + c1D_3(j,0)*u(i1,i2,ii,3)
                      end do
                   else
                      do j=-m,m
                         d1 = d1 + c1D_3(j,i3-n03+m)*u(i1,i2,n03 + j - m,3)
                      end do
                   end if
                else
                   do j=-m,m
                      d1 = d1 + c1D_3(j,0)*u(i1,i2,i3 + j,3)
                   end do
                end if
                du(i1,i2,i3) = du(i1,i2,i3)+d1
             end do
          end do
       end do
       !$omp end parallel do
    end if

  end subroutine div_u_i

  subroutine nabla_u(geocode,n01,n02,n03,u,du,nord,hgrids)
    implicit none

    !c..this routine computes 'nord' order accurate first derivatives 
    !c..on a equally spaced grid with coefficients from 'Matematica' program.

    !c..input:
    !c..ngrid       = number of points in the grid, 
    !c..u(ngrid)    = function values at the grid points

    !c..output:
    !c..du(ngrid)   = first derivative values at the grid points

    !c..declare the pass
    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
    integer, intent(in) :: n01,n02,n03,nord
    real(kind=8), dimension(3), intent(in) :: hgrids
    real(kind=8), dimension(n01,n02,n03), intent(in) :: u !<scalar field
    real(kind=8), dimension(n01,n02,n03,3) :: du !< nabla d_i u

    !c..local variables
    integer :: n,m,n_cell,ii
    integer :: i,j,ib,i1,i2,i3
    real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D, c1D_1, c1D_2, c1D_3
    real(kind=8) :: hx,hy,hz, d
    logical :: perx,pery,perz

    n = nord+1
    m = nord/2
    hx = hgrids(1)!acell/real(n01,kind=8)
    hy = hgrids(2)!acell/real(n02,kind=8)
    hz = hgrids(3)!acell/real(n03,kind=8)
    n_cell = max(n01,n02,n03)

    !buffers associated to the geocode
    !conditions for periodicity in the three directions
    perx=(geocode /= 'F')
    pery=(geocode == 'P')
    perz=(geocode /= 'F')


    ! Beware that n_cell has to be > than n.
    if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
    end if

    ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
    !Only nord=2,4,6,8,16

    select case(nord)
    case(2,4,6,8,16)
       !O.K.
    case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
    end select

    do i=-m,m
       do j=-m,m
          c1D(i,j)=0.d0
       end do
    end do

    include 'FiniteDiffCorff.inc'

    c1D_1 = c1D/hx
    c1D_2 = c1D/hy
    c1D_3 = c1D/hz
    !$omp parallel do default(shared) private(i1,i2,i3,j,ii, d) 
    do i3=1,n03
       do i2=1,n02
          do i1=1,n01

             d= 0.0d0

             if (i1.le.m) then
                if (perx) then
                   do j=-m,m
                      ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                      d = d + c1D_1(j,0)*u(ii,i2,i3)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3)
                   end do
                end if
             else if (i1.gt.n01-m) then
                if (perx) then
                   do j=-m,m
                      ii=modulo(i1 + j - 1, n01 ) + 1
                      d = d + c1D_1(j,0)*u(ii,i2,i3)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3)
                   end do
                end if
             else
                do j=-m,m
                   d = d + c1D_1(j,0)*u(i1 + j,i2,i3)
                end do
             end if
             du(i1,i2,i3,1) = d
          end do
       end do
    end do
    !$omp end parallel do

    !$omp parallel do default(shared) private(i1,i2,i3,j,ii, d) 
    do i3=1,n03
       do i2=1,n02
          do i1=1,n01
             d = 0.0d0

             if (i2.le.m) then
                if (pery) then
                   do j=-m,m
                      ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                      d = d + c1D_2(j,0)*u(i1,ii,i3)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_2(j,i2-m-1)*u(i1,j+m+1,i3)
                   end do
                end if
             else if (i2.gt.n02-m) then
                if (pery) then
                   do j=-m,m
                      ii=modulo(i2 + j - 1, n02 ) + 1
                      d = d + c1D_2(j,0)*u(i1,ii,i3)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_2(j,i2-n02+m)*u(i1,n02 + j - m,i3)
                   end do
                end if
             else
                do j=-m,m
                   d = d + c1D_2(j,0)*u(i1,i2 + j,i3)
                end do
             end if
             du(i1,i2,i3,2)=d
          end do
       end do
    end do
    !$omp end parallel do

    !$omp parallel do default(shared) private(i1,i2,i3,j,ii, d) 
    do i3=1,n03
       do i2=1,n02
          do i1=1,n01
             d = 0.0d0
             if (i3.le.m) then
                if (perz) then
                   do j=-m,m
                      ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                      d = d + c1D_3(j,0)*u(i1,i2,ii)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_3(j,i3-m-1)*u(i1,i2,j+m+1)
                   end do
                end if
             else if (i3.gt.n03-m) then
                if (perz) then
                   do j=-m,m
                      ii=modulo(i3 + j - 1, n03 ) + 1
                      d = d + c1D_3(j,0)*u(i1,i2,ii)
                   end do
                else
                   do j=-m,m
                      d = d + c1D_3(j,i3-n03+m)*u(i1,i2,n03 + j - m)
                   end do
                end if
             else
                do j=-m,m
                   d = d + c1D_3(j,0)*u(i1,i2,i3 + j)
                end do
             end if
             du(i1,i2,i3,3)=d
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine nabla_u

  !> Like fssnord3DmatNabla but corrected such that the index goes at the beginning
  !! Multiplies also times (nabla epsilon)/(4pi*epsilon)= nabla (log(epsilon))/(4*pi)
  subroutine update_rhopol(geocode,n01,n02,n03,u,nord,hgrids,eta,dlogeps,rhopol,rhores2)
    !use module_defs, only: pi_param
    implicit none

    !c..this routine computes 'nord' order accurate first derivatives 
    !c..on a equally spaced grid with coefficients from 'Matematica' program.

    !c..input:
    !c..ngrid       = number of points in the grid, 
    !c..u(ngrid)    = function values at the grid points

    !c..output:
    !c..du(ngrid)   = first derivative values at the grid points

    !c..declare the pass

    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
    integer, intent(in) :: n01,n02,n03,nord
    real(kind=8), intent(in) :: eta
    real(kind=8), dimension(3), intent(in) :: hgrids
    real(kind=8), dimension(n01,n02,n03), intent(in) :: u
    real(kind=8), dimension(3,n01,n02,n03), intent(in) :: dlogeps
    real(kind=8), dimension(n01,n02,n03), intent(inout) :: rhopol
    real(kind=8), intent(out) :: rhores2

    !c..local variables
    integer :: n,m,n_cell
    integer :: i,j,ib,i1,i2,i3,isp,i1_max,i2_max,ii
    !real(kind=8), parameter :: oneo4pi=0.25d0/pi_param
    real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D,c1DF
    real(kind=8) :: hx,hy,hz,max_diff,fact,dx,dy,dz,res,rho
    real(kind=8) :: oneo4pi,rpoints
    logical :: perx,pery,perz

    oneo4pi=1.0d0/(16.d0*atan(1.d0))

    n = nord+1
    m = nord/2
    hx = hgrids(1)!acell/real(n01,kind=8)
    hy = hgrids(2)!acell/real(n02,kind=8)
    hz = hgrids(3)!acell/real(n03,kind=8)
    n_cell = max(n01,n02,n03)
    rpoints=product(real([n01,n02,n03],dp))

    !buffers associated to the geocode
    !conditions for periodicity in the three directions
    perx=(geocode /= 'F')
    pery=(geocode == 'P')
    perz=(geocode /= 'F')

    ! Beware that n_cell has to be > than n.
    if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
    end if

    ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
    !Only nord=2,4,6,8,16
    if (all(nord /=[2,4,6,8,16])) then
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
    end if

    do i=-m,m
       do j=-m,m
          c1D(i,j)=0.d0
          c1DF(i,j)=0.d0
       end do
    end do

    include 'FiniteDiffCorff.inc'

    rhores2=0.d0
    !$omp parallel do default(shared) &
    !$omp private(i1,i2,i3,j,ii, dx,dy,dz,res,rho)&
!!!!$omp shared(m,n01,n02,n03,perx,pery,perz,rhopol,u,hx,hy,hz,c1D,eta,oneo4pi,dlogeps) &
    !$omp reduction(+:rhores2)
    do i3=1,n03
       do i2=1,n02
          do i1=1,n01
             dx=0.d0
             if (i1.le.m) then
                if (perx) then
                   do j=-m,m
                      ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                      dx = dx + c1D(j,0)*u(ii,i2,i3)
                   end do
                else
                   do j=-m,m
                      dx = dx + c1D(j,i1-m-1)*u(j+m+1,i2,i3)
                   end do
                end if
             else if (i1.gt.n01-m) then
                if (perx) then
                   do j=-m,m
                      ii=modulo(i1 + j - 1, n01 ) + 1
                      dx = dx + c1D(j,0)*u(ii,i2,i3)
                   end do
                else
                   do j=-m,m
                      dx = dx + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3)
                   end do
                end if
             else
                do j=-m,m
                   dx = dx + c1D(j,0)*u(i1 + j,i2,i3)
                end do
             end if
             dx=dx/hx
!!$        end do
!!$     end do
!!$  end do
!!$  !$omp end parallel do
!!$
!!$  !$omp parallel do default(shared) private(i1,i2,i3,j,ii, dy) 
!!$  do i3=1,n03
!!$     do i2=1,n02
!!$        do i1=1,n01
             dy = 0.0d0
             if (i2.le.m) then
                if (pery) then
                   do j=-m,m
                      ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                      dy = dy + c1D(j,0)*u(i1,ii,i3)
                   end do
                else
                   do j=-m,m
                      dy = dy + c1D(j,i2-m-1)*u(i1,j+m+1,i3)
                   end do
                end if
             else if (i2.gt.n02-m) then
                if (pery) then
                   do j=-m,m
                      ii=modulo(i2 + j - 1, n02 ) + 1
                      dy = dy + c1D(j,0)*u(i1,ii,i3)
                   end do
                else
                   do j=-m,m
                      dy = dy + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3)
                   end do
                end if
             else
                do j=-m,m
                   dy = dy + c1D(j,0)*u(i1,i2 + j,i3)
                end do
             end if
             dy=dy/hy
!!$        end do
!!$     end do
!!$  end do
!!$  !$omp end parallel do
!!$
!!$  !$omp parallel do default(shared) private(i1,i2,i3,j,ii, dz) 
!!$  do i3=1,n03
!!$     do i2=1,n02
!!$        do i1=1,n01
             dz = 0.0d0
             if (i3.le.m) then
                if (perz) then
                   do j=-m,m
                      ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                      dz = dz + c1D(j,0)*u(i1,i2,ii)
                   end do
                else
                   do j=-m,m
                      dz = dz + c1D(j,i3-m-1)*u(i1,i2,j+m+1)
                   end do
                end if
             else if (i3.gt.n03-m) then
                if (perz) then
                   do j=-m,m
                      ii=modulo(i3 + j - 1, n03 ) + 1
                      dz = dz + c1D(j,0)*u(i1,i2,ii)
                   end do
                else
                   do j=-m,m
                      dz = dz + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m)
                   end do
                end if
             else
                do j=-m,m
                   dz = dz + c1D(j,0)*u(i1,i2,i3 + j)
                end do
             end if
             dz=dz/hz

             !retrieve the previous treatment
             res = dlogeps(1,i1,i2,i3)*dx + &
                  dlogeps(2,i1,i2,i3)*dy + dlogeps(3,i1,i2,i3)*dz
             res = res*oneo4pi
             rho=rhopol(i1,i2,i3)
             res=res-rho
             res=eta*res
             rhores2=rhores2+res*res
             rhopol(i1,i2,i3)=res+rho


          end do
       end do
    end do
    !$omp end parallel do

    !this part should now go inside the open loop

!!$  !$omp parallel do default(shared) private(i1,i2,i3,res,rho) &
!!$  !$omp reduction(+:rhores2) 
!!$  do i3=1,n03
!!$     do i2=1,n02
!!$        do i1=1,n01
!!$           !retrieve the previous treatment
!!$           res = dlogeps(1,i1,i2,i3)*dx + &
!!$                dlogeps(2,i1,i2,i3)*dy + dlogeps(3,i1,i2,i3)*dz
!!$           res = res*oneo4pi
!!$           rho=rhopol(i1,i2,i3)
!!$           res=res-rho
!!$           res=eta*res
!!$           rhores2=rhores2+res*res
!!$           rhopol(i1,i2,i3)=res+rho
!!$        end do
!!$     end do
!!$  end do
!!$  !$omp end parallel do

    !  rhores2=rhores2/rpoints

  end subroutine update_rhopol

  !> verify that the density is considerably zero in the region where epsilon is different from one
  subroutine nonvacuum_projection(n1,n23,rho,oneoeps,norm)
    implicit none
    integer, intent(in) :: n1,n23 !< parallelized box dimensions
    real(dp), dimension(n1,n23), intent(in) :: rho !<charge density
    !>inverse of epsilon (might also be inverse of sqrt(eps))
    real(dp), dimension(n1,n23), intent(in) :: oneoeps 
    real(dp), intent(out) :: norm !< \int of rho where epsilon /=1
    !local variables
    integer :: i1,i23
    real(dp), parameter :: tol= 5.d-1

    norm=0.0_dp
    !$omp parallel do default(shared) private(i1,i23)&
    !$omp reduction(+:norm)
    do i23=1,n23
       do i1=1,n1
          if (abs(oneoeps(i1,i23) - 1.0_dp) > tol) norm=norm+rho(i1,i23)
       end do
    end do
    !$omp end parallel do

  end subroutine nonvacuum_projection



end module environment

