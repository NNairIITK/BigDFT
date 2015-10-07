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
  use PSbase
  use numerics, only: safe_exp,twopi,oneotwopi,oneofourpi,Bohr_Ang,AU_GPa,&
       dyn_AU,oneoeightpi
  use dictionaries, only: f_err_throw
  implicit none

  private

  real(gp), parameter, public :: vacuum_eps=1.0_gp


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

  !conversion factors in AU

  !> dyn/cm into atomic units (5.291772109217d-9/8.238722514d-3)
  real(gp), parameter :: SurfAU=Bohr_Ang*1.e-8/dyn_AU
  

  
  !> define the cavity type
  type, public :: cavity_data
     real(gp) :: epsilon0 !< dielectriv constant of the medium
     real(gp) :: edensmax !<maximum value of the density for the cavity
     real(gp) :: edensmin !<minimum  value of the density for the cavity
     real(dp) :: gammaS !< surface tension of the solvent [dyn/cm]
     real(dp) :: alphaS !< proportionality factor for the repulsion free energy in term of the surface integral [dyn/cm]
     real(dp) :: betaV !<proportionality factor for the dispersion free energy in term of the volume integral [GPa]     
  end type cavity_data

  public :: cavity_init,eps,epsprime,epssecond,oneoeps,oneosqrteps,logepsprime,corr_term
  public :: cavity_default,surf_term,epsle0,epsl,d1eps,dlepsdrho_sccs

contains

  pure function cavity_default() result(c)
    implicit none
    type(cavity_data) :: c
    c%epsilon0= 78.36_gp !<water at ambient condition 
    c%edensmax = 0.005_gp !0.0050d0
    c%edensmin = 0.0001_gp
    c%gammaS = 72._gp*SurfAU ![dyn/cm]   
    c%alphaS = -22.0_gp*SurfAU ![dyn/cm]   end function cavity_default
    c%betaV = -0.35_gp/AU_GPa ![GPa]     
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
    if (present(gammaS)) c%gammaS=gammaS*SurfAU
    if (present(alphaS)) c%alphaS=alphaS*SurfAU  
    if (present(betaV )) c%betaV =betaV/AU_GPa
  end function cavity_init

  subroutine dump_cavity(cavity)
    use yaml_output
    implicit none
    type(cavity_data), intent(in) :: cavity

    call yaml_mapping_open('Cavity parameters')
    call yaml_map('Epsilon',cavity%epsilon0)
    call yaml_map('edensmax',cavity%edensmax)
    call yaml_map('edensmin',cavity%edensmin)
    call yaml_map('gammaS',cavity%gammaS)
    call yaml_map('alphaS',cavity%alphaS)
    call yaml_map('betaV',cavity%betaV)
    call yaml_mapping_close()
  end subroutine dump_cavity


  pure function epsilon_transition(rho,pow,der,cavity) result(eps)
    implicit none
    character(len=*), intent(in) :: pow !<power to epsilon
    integer, intent(in) :: der !< derivative of epsolin
    real(dp), intent(in) :: rho
    type(cavity_data), intent(in) :: cavity
    real(gp) :: eps
    !local variables
    real(gp) :: fact1,fact2,fact0,r,zeta,w,dw,dzetao2pi,d2zetao2pi,dzeta2o2pi,d2w,ep,s,l

    eps=-1.d0 !this value should always be overwritten
    fact0=cavity%edensmax/cavity%edensmin
    r=cavity%edensmax/rho
    fact1=1.0_gp/log(fact0)
    fact2=log(r)
    zeta=twopi*fact1*fact2
    s=sin(zeta)
    w=oneotwopi*(zeta-s)
    dzetao2pi=-fact1/rho
    dw=(1.d0-cos(zeta))*dzetao2pi
    ep=cavity%epsilon0**w
    l=log(cavity%epsilon0)

    select case(pow)
    case('1')
       select case(der)
       case(0)
          eps=ep
       case(1)
          eps=ep*dw*l
       case(2)
          d2zetao2pi=fact1/rho**2
          dzeta2o2pi=twopi*dzetao2pi**2
          d2w=sin(zeta)*dzeta2o2pi+(1.d0-cos(zeta))*d2zetao2pi
          !eps=cavity%epsilon0*safe_exp(w)*(dw*dw+d2w)
          eps=l*ep*(l*dw*dw+d2w)
       end select
    case('-1/2')
       eps=cavity%epsilon0**(-0.5_gp*w)
    case('-1')
       eps=1.0_gp/ep
    case('L')
       select case(der)
       case(1)
          eps=dw*l
       end select
    case('C')
       !calculate the term 1/2 epsprime*logepsprime -epssecond, needed for the correction term
       !tt=dw*l
       !eps=0.5_gp*ep*dw*l*dw*l-l*ep*(l*dw*dw+d2w)
       d2zetao2pi=fact1/rho**2
       dzeta2o2pi=twopi*dzetao2pi**2
       d2w=sin(zeta)*dzeta2o2pi+(1.d0-cos(zeta))*d2zetao2pi

       eps=-l*ep*(0.5_gp*l*dw*dw+d2w)
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
    real(gp) :: epspr,fact1,fact2,fact3,r,t,coeff,coeff1,dtx,ep,w,ct

    !we are in a inner region
    if (rho > cavity%edensmax) then
       corr_term=0.0_gp
       !we are in a outer region
    else if (rho < cavity%edensmin) then
       corr_term=0.0_gp
    else
       epspr=epsprime(rho,cavity)
       ct=epsilon_transition(rho,'C',0,cavity)
!!$       corr_term=-0.125_gp/pi*(nabla2rho*&
!!$            (0.5_gp*epspr*logepsprime(rho,cavity)-epssecond(rho,cavity))-epspr*deltarho)
       corr_term=-oneoeightpi*(nabla2rho*ct-epspr*deltarho)

!!$       !old definition of the correction term
!!$       fact1=2.d0*pi/log(cavity%edensmax/cavity%edensmin)
!!$       fact2=(log(cavity%epsilon0))/(2.d0*pi)
!!$       fact3=(log(cavity%epsilon0))/log(cavity%edensmax/cavity%edensmin)
!!$
!!$       r=fact1*(log(cavity%edensmax/rho))
!!$       t=fact2*(r-sin(r))
!!$       w=(r-sin(r))/(2.0_gp*pi)
!!$       ep=cavity%epsilon0**w
!!$       coeff=fact3*(1.d0-cos(r))
!!$       dtx=-coeff/rho  !first derivative of t wrt rho
!!$       coeff1=(0.5d0*(coeff**2)+fact3*fact1*sin(r)+coeff)/(rho**2)
!!$       !corr_term=(0.125d0/pi)*safe_exp(t)*(coeff1*nabla2rho+dtx*deltarho)
!!$       corr_term=(0.125d0/pi)*ep*(coeff1*nabla2rho+dtx*deltarho)
    end if
  end function corr_term
  
  !>surface term multiplied by epsilon m1
  function surf_term(rho,nabla2rho,deltarho,ccrho,cavity)
    implicit none
    real(gp), intent(in) :: rho !<density
    real(gp), intent(in) :: nabla2rho !<square of the density gradient
    real(gp), intent(in) :: deltarho !<square of the density gradient
    real(gp), intent(in) :: ccrho !< u_i u_j d_i u_j , needed for the surface term where u_i=d_i rho
    type(cavity_data), intent(in) :: cavity
    real(gp) :: surf_term
    !local variables
    real(gp) :: de,c1,d
    !we are in a inner region
    if (rho > cavity%edensmax) then
       surf_term=0.0_gp
       !we are in a outer region
    else if (rho < cavity%edensmin) then
       surf_term=0.0_gp
    else
       de=epsprime(rho,cavity)
       d=sqrt(nabla2rho)
       c1=(ccrho/nabla2rho-deltarho)/d
       surf_term=de*c1
    end if

  end function surf_term

  pure function depsoeps(r,rc,delta,epsilon0)
    implicit none
    real(kind=8), intent(in) :: r,rc,delta,epsilon0
    real(kind=8) :: depsoeps
    !local variables
    real(kind=8) :: d

    d=(r-rc)/delta
    depsoeps=(epsilon0-1.d0)/delta*exp(-d**2)/epsle0(r,rc,delta,epsilon0)
  end function depsoeps

  pure function d1eps(r,rc,delta)
    use numerics, only: safe_exp,pi
    implicit none
    real(kind=8), intent(in) :: r,rc,delta
    real(kind=8) :: d1eps
    !local variables
    real(kind=8) :: d

    d=(r-rc)/delta
    d1eps=(1.d0/(delta*sqrt(pi)))*max(safe_exp(-d**2),1.0d-24)
  end function d1eps

  pure function epsl(r,rc,delta)
    implicit none
    real(kind=8), intent(in) :: r,rc,delta
    real(kind=8) :: epsl
    !local variables
    real(kind=8) :: d

    d=(r-rc)/delta
    epsl=0.5d0*(erf(d)+1.d0)
  end function epsl

  pure function epsle0(r,rc,delta,epsilon0)
    implicit none
    real(kind=8), intent(in) :: r,rc,delta,epsilon0
    real(kind=8) :: epsle0
    !local variables
    real(kind=8) :: d

    d=(r-rc)/delta
    epsle0=0.5d0*(epsilon0-1.d0)*(erf(d)+1.d0)+1.d0
  end function epsle0

  !> calculate the Extra potential and add it to the Hartree one
  !!at the same time evaluate the energy of the extra term given the 
  !! electronic charge density, and add if needed the ionic potential
  subroutine add_Vextra(n1,n23,nabla2_pot,depsdrho,dsurfdrho,cavity,&
       only_es,rho,sumpion,pot_ion,pot,eVextra)
    implicit none
    !>if .true., the added potential only comes from the 
    !!electrostatic contribution
    logical, intent(in) :: only_es,sumpion
    integer, intent(in) :: n1,n23
    !> on input, square of the gradient of the potential.
    !! on output, extra term of the potential
    real(dp), dimension(n1,n23), intent(in) :: depsdrho,dsurfdrho,rho,pot_ion
    type(cavity_data), intent(in) :: cavity
    real(dp), dimension(n1,n23), intent(in) :: nabla2_pot
    real(dp), dimension(n1,n23), intent(out) :: pot
    real(dp), intent(out) :: eVextra
    !local variables
    integer :: i1,i23
    real(dp) :: ep,sp,rh,pt

    eVextra=0.0_dp
    if (only_es) then
       if (sumpion) then
          !$omp parallel do default(shared) private(i1,i23,ep,rh,pt) &
          !$omp reduction(+:eVextra)       
          do i23=1,n23
             do i1=1,n1
                ep=depsdrho(i1,i23)
                pt=-oneoeightpi*ep*nabla2_pot(i1,i23)
                pot(i1,i23)=pt+pot_ion(i1,i23)
                rh=rho(i1,i23)*pt
                eVextra=eVextra+rh
             end do
          end do
          !$omp end parallel do
       else
          !$omp parallel do default(shared) private(i1,i23,ep,rh,pt) &
          !$omp reduction(+:eVextra)       
          do i23=1,n23
             do i1=1,n1
                ep=depsdrho(i1,i23)
                pt=-oneoeightpi*ep*nabla2_pot(i1,i23)
                pot(i1,i23)=pt
                rh=rho(i1,i23)*pt
                eVextra=eVextra+rh
             end do
          end do
          !$omp end parallel do
       end if
    else
       if (sumpion) then
          !$omp parallel do default(shared) private(i1,i23,ep,sp,rh,pt)&
          !$omp reduction(+:eVextra)
          do i23=1,n23
             do i1=1,n1
                ep=depsdrho(i1,i23)
                sp=dsurfdrho(i1,i23)
                pt=-oneoeightpi*ep*nabla2_pot(i1,i23)+&
                     (cavity%alphaS+cavity%gammaS)*sp+&
                     cavity%betaV*ep/(1.d0-cavity%epsilon0)
                pot(i1,i23)=pt+pot_ion(i1,i23)
                rh=rho(i1,i23)*pt
                eVextra=eVextra+rh
             end do
          end do
          !$omp end parallel do
       else
          !$omp parallel do default(shared) private(i1,i23,ep,sp,rh,pt)&
          !$omp reduction(+:eVextra)
          do i23=1,n23
             do i1=1,n1
                ep=depsdrho(i1,i23)
                sp=dsurfdrho(i1,i23)
                pt=-oneoeightpi*ep*nabla2_pot(i1,i23)+&
                     (cavity%alphaS+cavity%gammaS)*sp+&
                     cavity%betaV*ep/(1.d0-cavity%epsilon0)
                pot(i1,i23)=pt
                rh=rho(i1,i23)*pt
                eVextra=eVextra+rh
             end do
          end do
          !$omp end parallel do
       end if
    end if
  end subroutine add_Vextra

  !> calculate dlogepsilon with respect to rho in the sccs case
  subroutine dlepsdrho_sccs(ndims,rho,nabla_rho,epsinner,dlogepsilon,cavity)
    implicit none
    integer, dimension(3), intent(in) :: ndims
    type(cavity_data), intent(in) :: cavity
    real(dp), dimension(ndims(1),ndims(2),ndims(3)), intent(in) :: epsinner,rho
    real(dp), dimension(ndims(1),ndims(2),ndims(3),3), intent(in) :: nabla_rho
    real(dp), dimension(3,ndims(1),ndims(2),ndims(3)), intent(out) :: dlogepsilon


    !local variables
    real(dp), parameter :: innervalue = 0.9d0 !to be defined differently
    integer :: i1,i2,i3,i,n01,n02,n03
    real(dp) :: logepspr

    !aliasing
    n01=ndims(1)
    n02=ndims(2)
    n03=ndims(3)

    do i3=1,n03
       do i2=1,n02
          do i1=1,n01
             if (epsinner(i1,i2,i3).gt.innervalue) then ! Check for inner sccs cavity value to fix as vacuum
                do i=1,3
                   dlogepsilon(i,i1,i2,i3)=0.d0 !dlogeps(i,i1,i2,i3)
                end do
             else
                logepspr=logepsprime(rho(i1,i2,i3),cavity)
                do i=1,3
                   dlogepsilon(i,i1,i2,i3)=nabla_rho(i1,i2,i3,i)*logepspr
                end do

             end if
          end do
       end do
    end do
  end subroutine dlepsdrho_sccs
    

end module environment

