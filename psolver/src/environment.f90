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
  integer, parameter :: PS_PB_NONE = 1433
  integer, parameter :: PS_PB_LINEAR = 1434
  integer, parameter :: PS_PB_STANDARD = 1435
  integer, parameter :: PS_PB_MODIFIED = 1436

  type(f_enumerator), public :: PS_NONE_ENUM=f_enumerator('vacuum',PS_EPSILON_VACUUM,null())
  type(f_enumerator), public :: PS_RIGID_ENUM=f_enumerator('rigid',PS_EPSILON_RIGID_CAVITY,null())
  type(f_enumerator), public :: PS_SCCS_ENUM=f_enumerator('sccs',PS_EPSILON_SCCS,null())
  type(f_enumerator), public :: PS_PB_NONE_ENUM=f_enumerator('PB_NONE',PS_PB_NONE,null()) !< poisson boltzmann
  type(f_enumerator), public :: PS_PB_LINEAR_ENUM=f_enumerator('PB_LINEAR',PS_PB_LINEAR,null())
  type(f_enumerator), public :: PS_PB_STANDARD_ENUM=f_enumerator('PB_STANDARD',PS_PB_STANDARD,null())
  type(f_enumerator), public :: PS_PB_MODIFIED_ENUM=f_enumerator('PB_MODIFIED',PS_PB_MODIFIED,null())


  type(f_enumerator), parameter, public :: PS_VAC_ENUM=f_enumerator('VAC',PS_EPSILON_VACUUM,null())
  type(f_enumerator), parameter, public :: PS_PI_ENUM=f_enumerator('PI',PS_PI,null())
  type(f_enumerator), parameter, public :: PS_PCG_ENUM=f_enumerator('PCG',PS_PCG,null())

  !>threshold for comparison with zero
  real(dp), parameter :: thr=1.d-15

  !conversion factors in AU

  !> dyn/cm into atomic units (5.291772109217d-9/8.238722514d-3)
  real(gp), parameter, public :: SurfAU=Bohr_Ang*1.e-8/dyn_AU
  
  !> define the cavity type
  type, public :: cavity_data
     real(gp) :: epsilon0 !< dielectriv constant of the medium
     real(gp) :: edensmax !<maximum value of the density for the cavity
     real(gp) :: edensmin !<minimum  value of the density for the cavity
     real(gp) :: delta !<parameter for the PCM cavity in the case of rigid
     real(dp) :: gammaS !< surface tension of the solvent [dyn/cm]
     real(dp) :: alphaS !< proportionality factor for the repulsion free energy in term of the surface integral [dyn/cm]
     real(dp) :: betaV !<proportionality factor for the dispersion free energy in term of the volume integral [GPa]     
  end type cavity_data

  public :: cavity_init,eps,epsprime,epssecond,oneoeps,oneosqrteps,logepsprime,corr_term
  public :: cavity_default,surf_term,epsle0,epsl,d1eps,dlepsdrho_sccs,add_Vextra,PB_iteration
  public :: rigid_cavity_arrays,rigid_cavity_forces,nabla2pot_epsm1

contains

  pure function cavity_default() result(c)
    implicit none
    type(cavity_data) :: c
    c%epsilon0= 78.36_gp !<water at ambient condition 
    c%edensmax = 0.005_gp !0.0050d0
    c%edensmin = 0.0001_gp
    c%delta = 2.0_gp
    c%gammaS = 72._gp*SurfAU ![dyn/cm]   
    c%alphaS = -22.0_gp*SurfAU ![dyn/cm]   end function cavity_default
    c%betaV = -0.35_gp/AU_GPa ![GPa]     
  end function cavity_default

  !>initialize the cavity parameters
  function cavity_init(epsilon0,edensmax,edensmin,delta,gammaS,alphaS,betaV) result(c)
    implicit none
    real(gp), intent(in), optional :: epsilon0,edensmax,edensmin,gammaS,alphaS,betaV,delta
    type(cavity_data) :: c
    c=cavity_default()
    if (present(delta) .and. (present(edensmax) .or. present(edensmin)))&
         call f_err_throw('Only delta or edensmax(min) should be present')
    if (present(epsilon0)) c%epsilon0=epsilon0
    if (present(edensmax)) c%edensmax=edensmax
    if (present(edensmin)) c%edensmax=edensmin
    if (present(delta)) c%delta=delta
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
    integer, intent(in) :: der !< derivative of epsilon
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

  !rigid cavity terms
  subroutine rigid_cavity_arrays(cavity,mesh,v,nat,rxyz,radii,eps,deps,dleps,corr)
    use box
    implicit none
    type(cavity_data), intent(in) :: cavity
    type(cell), intent(in) :: mesh
    integer, intent(in) :: nat !< number of centres defining the cavity
    real(gp), dimension(nat), intent(in) :: radii !< radii of each of the atoms
    !>array of the position in the reference frame of rxyz
    real(gp), dimension(3), intent(in) :: v
    !> position of all the atoms in the grid coordinates
    real(gp), dimension(3,nat), intent(in) :: rxyz
    real(gp), intent(out) :: eps,deps,corr
    real(gp), dimension(3), intent(out) :: dleps
    !local variables
    integer :: iat
    real(gp) :: ep,dcorrha,rad,eh,d1e,dlogh,d2e,d,d2ha,dd
    real(gp), dimension(3) :: dha

    ep=1.0_dp
    dha=0.0_dp
    dcorrha=0.0_dp
    loop_at: do iat=1,nat
       rad=radii(iat)
       d=minimum_distance(mesh,v,rxyz(1,iat))
       if (d.eq.0.d0) then
        d=1.0d-30
        eh=epsl(d,rad,cavity%delta)
        ep=ep*eh
        d1e=0.0_dp
        d2e=0.0_dp
       else
        eh=epsl(d,rad,cavity%delta)
        ep=ep*eh
        d1e=d1eps(d,rad,cavity%delta)
        d2e=d2eps(d,rad,cavity%delta)
       end if
       if (ep < thr) then
          ep=0.0_dp
          exit loop_at
       end if
       if (abs(d1e) < thr) then
          dlogh=0.0_gp
       else
          dlogh=d1e/epsl(d,rad,cavity%delta)
          dcorrha=dcorrha+d2e/eh-(d1e**2)/eh**2+2.0_gp*d1e/eh/d
          dha=dha+dlogh*closest_r(mesh,v,center=rxyz(:,iat))/d
       end if
    end do loop_at
    eps=(cavity%epsilon0-vacuum_eps)*ep+vacuum_eps
    dleps=(eps-vacuum_eps)/eps*dha
    d2ha=square(mesh,dha)
    deps=(eps-vacuum_eps)*sqrt(d2ha)
    !corr=0.5_gp*(eps-vacuum_eps)/eps*(0.5_gp*d2ha*(1+eps)/eps+dcorrha)
    dd=(eps-vacuum_eps)*(d2ha+dcorrha)
    corr=-oneoeightpi*(0.5_gp*(eps-vacuum_eps)**2/eps*d2ha-dd)
  end subroutine rigid_cavity_arrays

  subroutine rigid_cavity_forces(only_electrostatic,cavity,mesh,v,&
       nat,rxyz,radii,epsr,npot2,fxyz,deps)
    use box
    implicit none
    type(cavity_data), intent(in) :: cavity
    type(cell), intent(in) :: mesh
    logical, intent(in) :: only_electrostatic
    integer, intent(in) :: nat !< number of centres defining the cavity
    real(gp), intent(in) :: npot2,epsr
    real(gp), dimension(nat), intent(in) :: radii !< radii of each of the atoms
    !>array of the position in the reference frame of rxyz
    real(gp), dimension(3), intent(in) :: v
    !> position of all the atoms in the grid coordinates
    real(gp), dimension(3,nat), intent(in) :: rxyz
    real(gp), dimension(3,nat), intent(inout) :: fxyz !<forces array
    real(gp), dimension(3), intent(in) :: deps !<gradient of epsilon(r)
    !local variables
    integer :: iat,i,j
    real(gp) :: d,dlogh,rad,tt,ttV,ttS,hh,epsrm1,eps0m1,sqdeps,ep,mm
    real(gp), dimension(3) :: f_Vterm,f_Sterm,depsdRi,vr,ddloghdRi,vect

    eps0m1=cavity%epsilon0-vacuum_eps
    hh=mesh%volume_element
    do iat=1,nat
       d=minimum_distance(mesh,v,rxyz(1,iat))
       rad=radii(iat)
       if (d.eq.0.d0) then
        d=1.0d-30
        dlogh=0.0_dp
       else
        dlogh=d1eps(d,rad,cavity%delta)
       end if
       if (abs(dlogh) < thr) cycle
       ep=epsl(d,rad,cavity%delta)
       dlogh=dlogh/ep
       epsrm1=epsr-vacuum_eps 
       if (abs(epsrm1) < thr) cycle
       vr(:)=closest_r(mesh,v,center=rxyz(:,iat))
       depsdRi(:)=-epsrm1*dlogh/d*vr(:)
       tt=-oneoeightpi*npot2*hh
       !here the forces can be calculated
       fxyz(:,iat)=fxyz(:,iat)+tt*depsdRi(:) ! Electrostatic force
       if (.not. only_electrostatic) then 
        ttV=cavity%betaV/eps0m1*hh ! CAN BE DONE EXTERNAL TO THE LOOP (INTEGRAL)
        f_Vterm(:)=ttV*depsdRi(:) ! Force from the Volume term to the energy
        sqdeps=sqrt(square(mesh,deps))
        ttS=-(cavity%alphaS+cavity%gammaS)/eps0m1/sqdeps*hh ! CAN BE DONE EXTERNAL TO THE LOOP (INTEGRAL), no sqdeps
        ddloghdRi(:)=((dlogh**2)/d-d2eps(d,rad,cavity%delta)/ep)*vr(:)
        do i=1,3
         vect(i)=0.0_dp
         do j=1,3
          mm=depsdRi(i)*deps(j)/epsrm1+epsrm1*(ddloghdRi(i)*vr(j)/d)
          vect(i)=vect(i)+mm*deps(j)
         end do
        end do
        f_Sterm(:)=ttS*(vect(:)-epsrm1*dlogh/d*deps(:)) ! Force from the Surface term to the energy
        fxyz(:,iat)=fxyz(:,iat)+f_Vterm(:)+f_Sterm(:)
       end if
    end do
  end subroutine rigid_cavity_forces


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
    d1eps=(1.d0/(delta*sqrt(pi)))*safe_exp(-d**2)
  end function d1eps

  function d2eps(r,rc,delta)
    implicit none
    real(kind=8), intent(in) :: r,rc,delta
    real(kind=8) :: d2eps
    !local variables
    real(kind=8) :: d

    d=(r-rc)/delta
    d2eps=-2.0_gp*d/delta*d1eps(r,rc,delta)
    
  end function d2eps

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
       only_es,sumpion,pot_ion,pot)
    implicit none
    !>if .true., the added potential only comes from the 
    !!electrostatic contribution
    logical, intent(in) :: only_es,sumpion
    integer, intent(in) :: n1,n23
    !> on input, square of the gradient of the potential.
    !! on output, extra term of the potential
    real(dp), dimension(n1,n23), intent(in) :: depsdrho,dsurfdrho,pot_ion
    type(cavity_data), intent(in) :: cavity
    real(dp), dimension(n1,n23), intent(in) :: nabla2_pot
    real(dp), dimension(n1,n23), intent(out) :: pot
    !local variables
    integer :: i1,i23
    real(dp) :: ep,sp,pt

    if (only_es) then
       if (sumpion) then
          !$omp parallel do default(shared) private(i1,i23,ep,pt)
          do i23=1,n23
             do i1=1,n1
                ep=depsdrho(i1,i23)
                pt=-oneoeightpi*ep*nabla2_pot(i1,i23)
                pot(i1,i23)=pt+pot_ion(i1,i23)
             end do
          end do
          !$omp end parallel do
       else
          !$omp parallel do default(shared) private(i1,i23,ep,pt)
          do i23=1,n23
             do i1=1,n1
                ep=depsdrho(i1,i23)
                pt=-oneoeightpi*ep*nabla2_pot(i1,i23)
                pot(i1,i23)=pt
             end do
          end do
          !$omp end parallel do
       end if
    else
       if (sumpion) then
          !$omp parallel do default(shared) private(i1,i23,ep,sp,pt)
          do i23=1,n23
             do i1=1,n1
                ep=depsdrho(i1,i23)
                sp=dsurfdrho(i1,i23)
                pt=-oneoeightpi*ep*nabla2_pot(i1,i23)+&
                     (cavity%alphaS+cavity%gammaS)*sp+&
                     cavity%betaV*ep/(vacuum_eps-cavity%epsilon0)
                pot(i1,i23)=pt+pot_ion(i1,i23)
             end do
          end do
          !$omp end parallel do
       else
          !$omp parallel do default(shared) private(i1,i23,ep,sp,pt)
          do i23=1,n23
             do i1=1,n1
                ep=depsdrho(i1,i23)
                sp=dsurfdrho(i1,i23)
                pt=-oneoeightpi*ep*nabla2_pot(i1,i23)+&
                     (cavity%alphaS+cavity%gammaS)*sp+&
                     cavity%betaV*ep/(vacuum_eps-cavity%epsilon0)
                pot(i1,i23)=pt
             end do
          end do
          !$omp end parallel do
       end if
    end if
  end subroutine add_Vextra

  subroutine nabla2pot_epsm1(n1,n23,eps,nabla2_pot,np2em1)
    implicit none
    integer, intent(in) :: n1,n23
    real(dp), dimension(n1*n23), intent(in) :: eps,nabla2_pot
    real(dp), dimension(n1*n23), intent(out) :: np2em1
    !local variables
    integer :: i123

    
    !$omp parallel do default(shared) private(i123)
    do i123=1,n1*n23
       np2em1(i123)=(eps(i123)-vacuum_eps)*nabla2_pot(i123)
    end do
    !$omp end parallel do

  end subroutine nabla2pot_epsm1

  !> Calculation of the Poisson-Boltzmann function.
  !! following the definitions given in J. J. López-García, J. Horno, C. Grosse Langmuir 27, 13970-13974 (2011).
  pure function PB_charge(epsilon,cavity,pot) result(ions_conc)

    use numerics, only: safe_exp

    implicit none

    !> argument of the Poisson-Botzmann function
    real(dp), intent(in) :: pot,epsilon
    type(cavity_data), intent(in) :: cavity
    real(dp) :: ions_conc
    ! Values to be given in the cavity structure
    integer, parameter :: n_ions = 2 !< number of ionic species in the dielectric liquid system
    real(dp), dimension(n_ions) :: z_ions !< valence of ionic species
    real(dp), dimension(n_ions) :: c_ions !< bulk concentations of ionic species [mol/m^3]
    real(dp), dimension(n_ions) :: c_max  !< maximum local concentration that ionic species can attain [mol/m^3]
    real(dp), dimension(n_ions) :: r_ions !< effective ionic radius of ionic species [m]
    real(dp), parameter :: Temp = 300 ! Temperature of the liquid system [K]
    !> packing coefficient p = 1 for perfect packing, p = pi_greek/(3(2)^{1/2}) ≈ 0.74 for close packing,
    real(dp), parameter :: p = 0.74d0 
    !! p ≈ 0.64 for random close packing, and p = pi_greek/6 ≈ 0.52 for simple cubic packing.
    ! Nedeed constant
    real(dp), parameter :: n_avo = 6.0221412927d23 ! Avogadro's number [1/mol]
    real(dp), parameter :: k_b = 3.166811429d-6 ! Boltzmann constant in atomic unit [E_{H}/K]
    real(dp), parameter :: bohr = 5.291772109217d-11 ! m
    !local variables
    integer :: i,j
    real(dp) :: pi,fact,vol_bohr,K_bT,t,fact1,sumc,y,h,l
    real(dp), dimension(n_ions) :: c_ratio  !< c_ions/c_max
    integer, parameter :: PBeq=3 ! Set 1 for linear, 2 for standard, 3 for
    ! modified Poisson-Boltzmann equation.

    pi = 4.d0*datan(1.d0)
    k_bT = k_b*Temp
    vol_bohr=bohr*bohr*bohr
    fact=n_avo*vol_bohr
    fact1=(4.d0/3.d0)*pi*n_avo
    l=0.d0

    z_ions(1)=1.0d0
    z_ions(2)=-1.0d0
    c_ions(1)=100.0d0
    c_ions(2)=100.0d0
    r_ions(1)=3.0d-10
    r_ions(2)=3.0d-10
    sumc=0.d0
    do i=1,n_ions
       c_max(i)=p/(fact1*(r_ions(1)**3))
       c_ratio(i)=c_ions(i)/c_max(i)
       sumc=sumc+c_ratio(i)
    end do

    !--------------------------------------------------------
    if (PBeq.eq.1) then
       ! Linear Poisson-Boltzmann Equation.
       ions_conc = 0.d0
       do i=1,n_ions
          t = -z_ions(i)*pot/k_bT !*0.01d0
          ions_conc = ions_conc + z_ions(i)*c_ions(i)*t
       end do
       ions_conc = ions_conc*fact!*1.d3  
    else if (PBeq.eq.2) then
       ! Standard Poisson-Boltzmann Equation.
       ions_conc = 0.d0
       do i=1,n_ions
          t = -z_ions(i)*pot/k_bT!*0.01d0
          t=safe_exp(t) ! Comment this line for linear Poisson-Boltzmann Equation.
          ions_conc = ions_conc + z_ions(i)*c_ions(i)*t
       end do
       ions_conc = ions_conc*fact!*1.d3  
    else if (PBeq.eq.3) then  
       ! Modified Poisson-Boltzmann Equation.
       ions_conc = 0.d0
       do i=1,n_ions
          y=pot/k_bT ! correct one
          !y=pot/k_bT*0.05d0
          t = -z_ions(i)*y 
          h=0.d0
          do j=1,n_ions
             h=h+c_ratio(j)*safe_exp((z_ions(i)-z_ions(j))*y)
          end do
          l=safe_exp(z_ions(i)*y)*(1.d0-sumc)+h
          t=1.d0/l
          ions_conc = ions_conc + z_ions(i)*c_ions(i)*t 
       end do
       ions_conc = ions_conc*fact ! correct one
       !ions_conc = ions_conc*fact*5.0d2
    end if
    !--------------------------------------------------------
    !multiply the value for the cavity parameter
    ions_conc=((epsilon-1.0d0)/(cavity%epsilon0-1.0d0))*ions_conc

  end function PB_charge


  !> take the input density and the present potenital and provdes the new charge density for GPS operation
  subroutine PB_iteration(n1,n23,i23s,eta,cavity,density,pot,eps,rho_ions,rho_new,res_PB)
    implicit none
    integer, intent(in) :: n1,n23,i23s
    real(dp), intent(in) :: eta
    type(cavity_data), intent(in) :: cavity
    real(dp), dimension(n1,n23), intent(in) :: density,eps
    real(dp), dimension(n1,*), intent(in) :: pot
    real(dp), dimension(n1,n23), intent(inout) :: rho_ions
    real(dp), intent(out) :: res_PB
    real(dp), dimension(n1,n23), intent(out) :: rho_new
    !local variables
    integer :: i1,i23
    real(dp) :: res,rho

    res_PB=0.d0
    do i23=1,n23
       do i1=1,n1
          res=PB_charge(eps(i1,i23),cavity,pot(i1,i23s+i23)) ! Additional contribution to the Generalized Poisson operator
          ! for the Poisson-Boltzmann equation.
          rho=rho_ions(i1,i23)
          res=res-rho
          res=eta*res
          res_PB=res_PB+res*res
          rho_ions(i1,i23)=res+rho
          rho_new(i1,i23) = density(i1,i23) + rho_ions(i1,i23)
       end do
    end do

  end subroutine PB_iteration


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
