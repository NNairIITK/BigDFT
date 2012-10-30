!> @file
!!  Define the fortran types
!! @author
!!    Copyright (C) 2008-2011 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 
!> Modules which contains the handling of operations in Gaussian basis
!! Spherical harmonics are used in the cartesian form
module gaussians

  use module_base, only:gp

  private

  !> Structures of basis of gaussian functions
  type, public :: gaussian_basis
     integer :: nat  !< number of centers
     integer :: ncoeff !< number of total basis elements
     integer :: nshltot !< total number of shells (m quantum number ignored) 
     integer :: nexpo !< number of exponents (sum of the contractions)
     !storage units
     integer, dimension(:), pointer :: nshell !< number of shells for any of the centers
     integer, dimension(:), pointer :: ndoc,nam !< degree of contraction, angular momentum of any shell
     real(gp), dimension(:), pointer :: xp,psiat !<factors and values of the exponents (complex numbers are allowed)
     real(gp), dimension(:,:), pointer :: rxyz !<positions of the centers
  end type gaussian_basis

  public :: gaudim_check,normalize_shell,gaussian_overlap,kinetic_overlap

contains

  
  !>   Overlap matrix between two different basis structures
  subroutine gaussian_overlap(A,B,ovrlp)
    implicit none
    type(gaussian_basis), intent(in) :: A,B
    real(gp), dimension(A%ncoeff,B%ncoeff) :: ovrlp 
    !only lower triangular part for A%ncoeff=B%ncoeff
    !local variables
    integer, parameter :: niw=18,nrw=6
    integer :: ishell,iexpo,icoeff,iat,jat,isat,jsat,jshell
    integer :: iovrlp,jovrlp,jcoeff,jexpo
    integer :: ngA,ngB,lA,lB,mA,mB
    real(gp) :: dx,dy,dz
    integer, dimension(niw) :: iw
    real(gp), dimension(nrw) :: rw

    iovrlp=0
    ishell=0
    iexpo=1
    icoeff=1

    !loop on each shell (intensive calculation)
    do iat=1,A%nat
       do isat=1,A%nshell(iat)
          ishell=ishell+1
          ngA=A%ndoc(ishell)
          lA=A%nam(ishell)
          do mA=1,2*lA-1
             iovrlp=iovrlp+1

             jovrlp=0
             jshell=0
             jexpo=1
             jcoeff=1

             do jat=1,B%nat
                dx=B%rxyz(1,jat)-A%rxyz(1,iat)
                dy=B%rxyz(2,jat)-A%rxyz(2,iat)
                dz=B%rxyz(3,jat)-A%rxyz(3,iat)
                do jsat=1,B%nshell(jat)
                   jshell=jshell+1
                   ngB=B%ndoc(jshell)
                   lB=B%nam(jshell)
                   do mB=1,2*lB-1
                      jovrlp=jovrlp+1
                      !if ((jovrlp >= iovrlp .and. A%ncoeff == B%ncoeff) .or. &
                      !     A%ncoeff /= B%ncoeff ) then
                      call gbasovrlp(A%xp(iexpo),A%psiat(iexpo),&
                           B%xp(jexpo),B%psiat(jexpo),&
                           ngA,ngB,lA,mA,lB,mB,dx,dy,dz,&
                           niw,nrw,iw,rw,ovrlp(iovrlp,jovrlp))
                      !end if
                   end do
                   jexpo=jexpo+ngB
                   jcoeff=jcoeff+2*lB-1
                end do
             end do
          end do
          iexpo=iexpo+ngA
          icoeff=icoeff+2*lA-1
       end do
    end do

    call gaudim_check(iexpo,icoeff,ishell,A%nexpo,A%ncoeff,A%nshltot)
    call gaudim_check(jexpo,jcoeff,jshell,B%nexpo,B%ncoeff,B%nshltot)

  END SUBROUTINE gaussian_overlap

  !>   Overlap kinetic matrix between two different basis structures
  !!   the kinetic operator is applicated on the A basis structure
  subroutine kinetic_overlap(A,B,ovrlp)
    implicit none
    type(gaussian_basis), intent(in) :: A,B
    real(gp), dimension(A%ncoeff,B%ncoeff) :: ovrlp 
    !only lower triangular part for A%ncoeff=B%ncoeff
    !local variables
    integer, parameter :: niw=18,nrw=6
    integer :: ishell,iexpo,icoeff,iat,jat,isat,jsat,jshell
    integer :: iovrlp,jovrlp,jcoeff,jexpo
    integer :: ngA,ngB,lA,lB,mA,mB
    real(gp) :: dx,dy,dz
    integer, dimension(niw) :: iw
    real(gp), dimension(nrw) :: rw

    iovrlp=0
    ishell=0
    iexpo=1
    icoeff=1

    !loop on each shell (intensive calculation)
    do iat=1,A%nat
       do isat=1,A%nshell(iat)
          ishell=ishell+1
          ngA=A%ndoc(ishell)
          lA=A%nam(ishell)
          do mA=1,2*lA-1
             iovrlp=iovrlp+1

             jovrlp=0
             jshell=0
             jexpo=1
             jcoeff=1

             do jat=1,B%nat
                dx=B%rxyz(1,jat)-A%rxyz(1,iat)
                dy=B%rxyz(2,jat)-A%rxyz(2,iat)
                dz=B%rxyz(3,jat)-A%rxyz(3,iat)
                do jsat=1,B%nshell(jat)
                   jshell=jshell+1
                   ngB=B%ndoc(jshell)
                   lB=B%nam(jshell)
                   do mB=1,2*lB-1
                      jovrlp=jovrlp+1
                      if (jovrlp >= iovrlp .and. A%ncoeff == B%ncoeff .or. &
                           A%ncoeff /= B%ncoeff ) then
                         call kineticovrlp(A%xp(iexpo),A%psiat(iexpo),&
                              B%xp(jexpo),B%psiat(jexpo),&
                              ngA,ngB,lA,mA,lB,mB,dx,dy,dz,&
                              niw,nrw,iw,rw,ovrlp(iovrlp,jovrlp))
                      end if
                   end do
                   jexpo=jexpo+ngB
                   jcoeff=jcoeff+2*lB-1
                end do
             end do
          end do
          iexpo=iexpo+ngA
          icoeff=icoeff+2*lA-1
       end do
    end do

    call gaudim_check(iexpo,icoeff,ishell,A%nexpo,A%ncoeff,A%nshltot)
    call gaudim_check(jexpo,jcoeff,jshell,B%nexpo,B%ncoeff,B%nshltot)

  END SUBROUTINE kinetic_overlap

  !>   Calculates the scalar product between two shells
  !!   by considering only the nonzero coefficients
  !!   actual building block for calculating overlap matrix
  !!   inserted work arrays for calculation
  subroutine gbasovrlp(expo1,coeff1,expo2,coeff2,ng1,ng2,l1,m1,l2,m2,dx,dy,dz,&
       niw,nrw,iw,rw,ovrlp)
    implicit none
    integer, intent(in) :: ng1,ng2,l1,m1,l2,m2,niw,nrw
    real(gp), intent(in) :: dx,dy,dz
    integer, dimension(niw) :: iw
    real(gp), dimension(nrw) :: rw
    real(gp), dimension(ng1), intent(in) :: expo1,coeff1
    real(gp), dimension(ng2), intent(in) :: expo2,coeff2
    real(gp), intent(out) :: ovrlp
    !local variables
    integer :: i1,i2
    real(gp) :: a1,a2,c1,c2,govrlpr

    ovrlp=0.d0
    do i1=1,ng1
       a1=expo1(i1)
       a1=0.5_gp/a1**2
       c1=coeff1(i1)
       do i2=1,ng2
          a2=expo2(i2)
          a2=0.5_gp/a2**2
          c2=coeff2(i2)
          call gprod(a1,a2,dx,dy,dz,l1,m1,l2,m2,niw,nrw,iw,rw,govrlpr)
          govrlpr=c1*govrlpr*c2
          !print *,c1,c2,govrlpr
          ovrlp=ovrlp+govrlpr
       end do
    end do

  END SUBROUTINE gbasovrlp

  !>   Calculates the scalar product between two shells
  !!   by considering only the nonzero coefficients
  !!   actual building block for calculating overlap matrix
  !!   inserted work arrays for calculation
  subroutine kineticovrlp(expo1,coeff1,expo2,coeff2,ng1,ng2,l1,m1,l2,m2,dx,dy,dz,&
       niw,nrw,iw,rw,ovrlp)
    implicit none
    integer, intent(in) :: ng1,ng2,l1,m1,l2,m2,niw,nrw
    real(gp), intent(in) :: dx,dy,dz
    integer, dimension(niw) :: iw
    real(gp), dimension(nrw) :: rw
    real(gp), dimension(ng1), intent(in) :: expo1,coeff1
    real(gp), dimension(ng2), intent(in) :: expo2,coeff2
    real(gp), intent(out) :: ovrlp
    !local variables
    integer :: i1,i2
    real(gp) :: a1,a2,c1,c2,govrlpr

    ovrlp=0.d0
    do i1=1,ng1
       a1=expo1(i1)
       a1=0.5_gp/a1**2
       c1=coeff1(i1)
       do i2=1,ng2
          a2=expo2(i2)
          a2=0.5_gp/a2**2
          c2=coeff2(i2)
          call kinprod(a1,a2,dx,dy,dz,l1,m1,l2,m2,niw,nrw,iw,rw,govrlpr)
          govrlpr=c1*govrlpr*c2
          !print *,c1,c2,govrlpr
          ovrlp=ovrlp+govrlpr
       end do
    end do

  END SUBROUTINE kineticovrlp

  !>   Calculates a dot product between two differents gaussians times spherical harmonics
  !!   valid only for shell which belongs to different atoms, and with also dy/=0/=dx dz/=0
  !!   to be rearranged when only some of them is zero
  subroutine gprod(a1,a2,dx,dy,dz,l1,m1,l2,m2,niw,nrw,iw,rw,ovrlp)
    implicit none
    integer, intent(in) :: l1,l2,m1,m2,niw,nrw 
    real(gp), intent(in) :: a1,a2,dx,dy,dz
    integer, dimension(niw) :: iw !work array of the exponents of the two polynomials
    real(gp), dimension(nrw) :: rw !work array of the polynomials coefficients 
    real(gp), intent(out) :: ovrlp
    !local variables
    integer, parameter :: nx=3
    integer :: n1,n2,i1,i2,px,py,pz,qx,qy,qz
    real(gp) :: fx,fy,fz,fa,fb!,govrlp

    !calculates the number of different couples
    call calc_coeff_inguess(l1,m1,nx,n1,&
         iw(1),iw(nx+1),iw(2*nx+1),rw(1))
    call calc_coeff_inguess(l2,m2,nx,n2,&
         iw(3*nx+1),iw(4*nx+1),iw(5*nx+1),rw(n1+1))
    ovrlp=0.d0
    do i2=1,n2
       qx=iw(3*nx+i2)
       qy=iw(4*nx+i2)
       qz=iw(5*nx+i2)
       fb=rw(n1+i2)
       do i1=1,n1
          px=iw(i1)
          py=iw(nx+i1)
          pz=iw(2*nx+i1)
          fa=rw(i1)

          fx=govrlp(a1,a2,dx,px,qx)
          fy=govrlp(a1,a2,dy,py,qy)
          fz=govrlp(a1,a2,dz,pz,qz)

          ovrlp=ovrlp+fa*fb*fx*fy*fz
          !print *,i1,i2,fx,fy,fz,fa,fb
       end do
    end do

  END SUBROUTINE gprod

!!$
!!$  !> Calculates a dot product between two basis functions
!!$  !! Basis function is identified by its quantum numbers and sigmas
!!$  !! the contraction coefficients are also given
!!$  !! the difference in positions between the basis centers is given in 
!!$  !! orthogonal cartesian coordinates
!!$  subroutine gdot(ng1,d1,s1,n1,l1,m1,ng2,d2,s2,n2,l2,m2,dr,ovrlp)
!!$    implicit none
!!$    integer, intent(in) :: ng1,ng2,n1,n2,l1,l2,m1,m2
!!$    real(gp), dimension(3), intent(in) :: dr
!!$    real(gp), dimension(ng1), intent(in) :: s1,d1 !<they should be contiguous
!!$    real(gp), dimension(ng2), intent(in) :: s2,d2
!!$    real(gp), intent(out) :: ovrlp
!!$    !local variables
!!$    integer, parameter :: nx=3
!!$    integer :: n1,n2,i1,i2,px,py,pz,qx,qy,qz
!!$    real(gp) :: fx,fy,fz,fa,fb!,govrlp
!!$
!!$    !id of tensor product decompositions (derivatives to be added)
!!$    lm1=l1**2+m1
!!$    lm2=l2**2+m2
!!$    ovrlp=0.0_gp
!!$    do ig2=1,ng2
!!$       do ig1=1,ng1
!!$          integral=0.0_gp
!!$          do i2=1,ntpd(lm2,n2)
!!$             do i1=1,ntpd(lm1,n1)
!!$                f(1)=govrlp(s1(ig1),s2(ig2),dr(1),&
!!$                     pws(1,i1,lm1,n1),pws(1,i2,lm2,n2))
!!$                f(2)=govrlp(s1(ig1),s2(ig2),dr(2),&
!!$                     pws(2,i1,lm1,n1),pws(2,i2,lm2,n2))
!!$                f(3)=govrlp(s1(ig1),s2(ig2),dr(3),&
!!$                     pws(3,i1,lm1,n1),pws(3,i2,lm2,n2))
!!$                integral=integral+ftpd(i1,lm1,n1)*ftpd(i2,lm2,n2)*f(1)*f(2)*f(3)
!!$             end do
!!$          end do
!!$          ovrlp=ovrlp+d1(ig1)*d2(ig2)*integral
!!$       end do
!!$    end do
!!$
!!$  END SUBROUTINE gdot



  !>   Calculates @f$\int \exp^{-a1*x^2} x^l1 \exp^{-a2*(x-d)^2} (x-d)^l2 dx@f$
  !!   Uses gauint0 if d==0
  !!
  pure function govrlp(a1,a2,d,l1,l2)
    implicit none
    integer, intent(in) :: l1,l2
    real(gp), intent(in) :: a1,a2,d
    real(gp) :: govrlp
    !local variables
    integer :: p
    real(gp) :: prefac,stot,aeff,ceff,tt,fsum!,gauint,gauint0

    !quick check
    if (d==0.0_gp) then
       govrlp=gauint0(a1+a2,l1+l2)
       return
    end if

    !build the prefactor
    prefac=a1+a2
    prefac=a2/prefac
    prefac=a1*prefac
    prefac=-d**2*prefac
    prefac=dexp(prefac)

    !build the effective exponent and coefficients
    aeff=a1+a2
    ceff=a2*d
    ceff=ceff/aeff

    !build the first term in the sum
    stot=(-d)**l2
    stot=gauint(aeff,ceff,l1)*stot

    !perform the sum
    do p=1,l2/2
       tt=rfac(1,p)
       fsum=rfac(l2-p+1,l2)
       fsum=fsum/tt
       tt=(-d)**(l2-p)
       fsum=fsum*tt
       tt=gauint(aeff,ceff,l1+p)
       fsum=fsum*tt
       stot=stot+fsum
    end do
    do p=l2/2+1,l2
       tt=rfac(1,l2-p)
       fsum=rfac(p+1,l2)
       fsum=fsum/tt
       tt=(-d)**(l2-p)
       fsum=fsum*tt
       tt=gauint(aeff,ceff,l1+p)
       fsum=fsum*tt
       stot=stot+fsum
    end do

    !final result
    govrlp=prefac*stot
  END FUNCTION govrlp
  
  !>   Kinetic overlap between gaussians, based on cartesian coordinates
  !!   calculates a dot product between two differents gaussians times spherical harmonics
  !!   valid only for shell which belongs to different atoms, and with also dy/=0/=dx dz/=0
  !!   to be rearranged when only some of them is zero
  !!
  subroutine kinprod(a1,a2,dx,dy,dz,l1,m1,l2,m2,niw,nrw,iw,rw,ovrlp)
    implicit none
    integer, intent(in) :: l1,l2,m1,m2,niw,nrw 
    real(gp), intent(in) :: a1,a2,dx,dy,dz
    integer, dimension(niw) :: iw !work array of the exponents of the two polynomials
    real(gp), dimension(nrw) :: rw !work array of the polynomials coefficients 
    real(gp), intent(out) :: ovrlp
    !local variables
    integer, parameter :: nx=3
    integer :: n1,n2,i1,i2,px,py,pz,qx,qy,qz
    real(gp) :: fx,fy,fz,fa,fb,d2fx,d2fy,d2fz!,govrlp,kinovrlp

    !calculates the number of different couples
    call calc_coeff_inguess(l1,m1,nx,n1,&
         iw(1),iw(nx+1),iw(2*nx+1),rw(1))
    call calc_coeff_inguess(l2,m2,nx,n2,&
         iw(3*nx+1),iw(4*nx+1),iw(5*nx+1),rw(n1+1))
    ovrlp=0.d0
    do i2=1,n2
       qx=iw(3*nx+i2)
       qy=iw(4*nx+i2)
       qz=iw(5*nx+i2)
       fb=rw(n1+i2)
       do i1=1,n1
          px=iw(i1)
          py=iw(nx+i1)
          pz=iw(2*nx+i1)
          fa=rw(i1)

          fx=govrlp(a1,a2,dx,px,qx)
          fy=govrlp(a1,a2,dy,py,qy)
          fz=govrlp(a1,a2,dz,pz,qz)

          d2fx=kinovrlp(a1,a2,dx,px,qx)
          d2fy=kinovrlp(a1,a2,dy,py,qy)
          d2fz=kinovrlp(a1,a2,dz,pz,qz)

          ovrlp=ovrlp-0.5_gp*fa*fb*(d2fx*fy*fz+fx*d2fy*fz+fx*fy*d2fz)
          !print *,i1,i2,fx,fy,fz,fa,fb
       end do
    end do

  END SUBROUTINE kinprod

  !>   Calculates @f$\int d^2/dx^2(\exp^{-a1*x^2} x^l1) \exp^{-a2*(x-d)^2} (x-d)^l2 dx@f$
  !!   in terms of the govrlp function below
  pure function kinovrlp(a1,a2,d,l1,l2)
    implicit none
    integer, intent(in) :: l1,l2
    real(gp), intent(in) :: a1,a2,d
    real(gp) :: kinovrlp
    !local variables
    real(gp) :: fac,ovrlp!govrlp

    !case l1+2
    fac=4._gp*a1**2
    ovrlp=govrlp(a1,a2,d,l1+2,l2)
    kinovrlp=fac*ovrlp
    !case l1
    fac=2._gp*a1*real(2*l1+1,gp)
    ovrlp=govrlp(a1,a2,d,l1,l2)
    kinovrlp=kinovrlp-fac*ovrlp
    !case l1-2 (if applicable)
    if (l1 >=2) then
       fac=real(l1*(l1-1),gp)
       ovrlp=govrlp(a1,a2,d,l1-2,l2)
       kinovrlp=kinovrlp+fac*ovrlp
    end if
  END FUNCTION kinovrlp

  !>   Calculates @f$\int \exp^{-a*x^2} x^l dx@f$
  !!   this works for all l
  pure function gauint0(a,l)
    implicit none
    integer, intent(in) :: l
    real(gp), intent(in) :: a
    real(gp) :: gauint0
    !local variables
    real(gp), parameter :: gammaonehalf=1.772453850905516027298d0
    integer :: p
    real(gp) :: prefac,tt
    !build the prefactor
    prefac=sqrt(a)
    prefac=1.d0/prefac
    prefac=gammaonehalf*prefac**(l+1)

    p=l/2
    if (2*p < l) then
       gauint0=0.0_gp
       return
    end if
    tt=xfac(1,p,-0.5d0)
    !final result
    gauint0=prefac*tt

  END FUNCTION gauint0



  !>   Calculates @f$\int \exp^{-a*(x-c)^2} x^l dx@f$
  !!   this works ONLY when c/=0.d0
  !!
  !!
  pure function gauint(a,c,l)
    implicit none
    integer, intent(in) :: l
    real(gp), intent(in) :: a,c
    real(gp) :: gauint
    !local variables
    real(gp), parameter :: gammaonehalf=1.772453850905516027298d0
    integer :: p
    real(gp) :: prefac,stot,fsum,tt!,firstprod

    !quick check
    !if (c==0.0_gp) then
    !   stop 'gauint0 should be called'
    !end if

    !build the prefactor
    prefac=sqrt(a)
    prefac=1.d0/prefac
    prefac=gammaonehalf*prefac

    !the first term of the sum is one
    !but we have to multiply for the prefactor
    stot=c**l

    !if (c==0.0_gp .and. l==0) then
    !   do p=0,20
    !      print *,'stot,p',stot,a,p,gauint0(a,p)
    !   end do
    !end if

    !calculate the sum
    do p=1,l/4
       tt=rfac(p+1,2*p)
       fsum=rfac(l-2*p+1,l)
       fsum=fsum/tt
       tt=firstprod(p)
       fsum=fsum*tt
       tt=c**(l-2*p)
       tt=tt/a**p
       fsum=fsum*tt
       stot=stot+fsum
    end do
    do p=l/4+1,l/3
       tt=rfac(p+1,l-2*p)
       fsum=rfac(2*p+1,l)
       fsum=fsum/tt
       tt=firstprod(p)
       fsum=fsum*tt
       tt=c**(l-2*p)
       tt=tt/a**p
       fsum=fsum*tt
       stot=stot+fsum
    end do
    do p=l/3+1,l/2
       tt=rfac(l-2*p+1,p)
       fsum=rfac(2*p+1,l)
       fsum=fsum*tt
       tt=firstprod(p)
       fsum=fsum*tt
       tt=c**(l-2*p)
       tt=tt/a**p
       fsum=fsum*tt
       stot=stot+fsum
    end do

    !final result
    gauint=stot*prefac

  END FUNCTION gauint

  !>
  pure function firstprod(p)
    implicit none
    integer, intent(in) :: p
    real(gp) :: firstprod
    !local variables
    integer :: i
    real(gp) :: tt
    firstprod=1.0_gp
    do i=1,p
       tt=real(2*i,gp)
       tt=1.0_gp/tt
       tt=1.0_gp-tt
       firstprod=firstprod*tt
    end do
  END FUNCTION firstprod
  
  !>
  subroutine gaudim_check(iexpo,icoeff,ishell,nexpo,ncoeff,nshltot)
    implicit none
    integer, intent(in) :: iexpo,icoeff,ishell,nexpo,ncoeff,nshltot
    !check of the dimensions
    if (iexpo /= nexpo+1) then
       write(*,*)' ERROR: nexpo+1 <> iexpo',nexpo,iexpo
       stop
    else if (icoeff /= ncoeff+1) then
       write(*,*)' ERROR: ncoeff+1 <> icoeff',ncoeff,icoeff
       stop
    else if (ishell /= nshltot) then
       write(*,*)' ERROR: nshltot <> ishell',nshltot,ishell
       stop
    end if
  END SUBROUTINE gaudim_check

  !>   Normalize a given atomic shell following the angular momentum
  !!
  !!
  pure subroutine normalize_shell(ng,l,expo,coeff)
    implicit none
    integer, intent(in) :: ng,l
    real(gp), dimension(ng), intent(in) :: expo
    real(gp), dimension(ng), intent(inout) :: coeff
    !local variables
    integer :: i,j
    real(gp) :: norm,tt,e1,ex,c1,c2!,gauint0

    norm=0.d0
    do i=1,ng
       e1=expo(i)
       c1=coeff(i)
       do j=1,ng
          ex=expo(j)+e1
          c2=coeff(j)
          tt=gauint0(ex,2*l+2)
          norm=norm+c1*tt*c2
       end do
    end do
    norm=sqrt(0.5_gp*norm)
    norm=1.0_gp/norm
    do i=1,ng
       coeff(i)=coeff(i)*norm
    end do

    !print *,'l=',l,'norm=',norm

  END SUBROUTINE normalize_shell

  !> Factorial (float)
  pure function rfac(is,ie)
    implicit none
    integer, intent(in) :: is,ie
    real(gp) :: rfac
    !local variables
    integer :: i
    real(gp) :: tt
    rfac=1.d0
    do i=is,ie
       tt=real(i,gp)
       rfac=rfac*tt
    end do
  END FUNCTION rfac

  !> Partial factorial, with real shift
  !!With this function n!=xfac(1,n,0.d0)
  pure function xfac(is,ie,sh)
    implicit none
    integer, intent(in) :: is,ie
    real(gp), intent(in) :: sh
    real(gp) :: xfac
    !local variables
    integer :: i
    real(gp) :: tt
    xfac=1.d0
    do i=is,ie
       tt=real(i,gp)+sh
       xfac=xfac*tt
    end do
  END FUNCTION xfac

end module gaussians
