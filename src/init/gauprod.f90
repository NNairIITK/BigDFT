!!$!calculates the overlap matrix of the gaussian basis functions
!!$subroutine gaussian_overlap(ngau,nterm,expcent,coeff,lxyz,overlap)
!!$  implicit none
!!$  integer, intent(in) :: ngau,nterm
!!$  !exponents of the spherical harmonics nterm=primitve gaussians * spherical components
!!$  integer, dimension(3,nterm,ngau), intent(in) :: lxyz 
!!$  !coefficients of the spherical harmonics times primitive gaussian functions
!!$  real(gp), dimension(nterm,ngau), intent(in) :: coeff 
!!$  !exponent and center
!!$  real(gp), dimension(2,nterm,ngau), intent(in) :: expcent
!!$  
!!$end subroutine gaussian_overlap

!normalize a given atomic shell following the angular momentum
subroutine normalize_shell(ng,l,expo,coeff)
  use module_base
  implicit none
  integer, intent(in) :: ng,l
  real(gp), dimension(ng), intent(in) :: expo
  real(gp), dimension(ng), intent(inout) :: coeff
  !local variables
  integer :: i,j
  real(gp) :: norm,tt,e1,ex,c1,c2,gauint0

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

end subroutine normalize_shell

!calculates \int_0^\infty \exp^{-a*x^2} x^l dx
function gauinth(a,l)
  use module_base
  implicit none
  integer, intent(in) :: l
  real(gp), intent(in) :: a
  real(gp) :: gauinth
  !local variables
  real(gp), parameter :: gammaonehalf=1.772453850905516027298d0
  integer :: p
  real(gp) :: xfac,prefac,tt,firstprod,sh
  !build the prefactor
  prefac=sqrt(a)
  prefac=1.d0/prefac
  prefac=prefac**(l+1)
  p=1-l+2*(l/2)
  tt=0.5d0*gammaonehalf**p
  prefac=prefac*tt
  sh=-0.5d0*real(p,gp)
  p=l/2
  tt=xfac(1,p,sh)
  !final result
  gauinth=prefac*tt
  
end function gauinth

!calculates the scalar product between two shells
!by considering only the nonzero coefficients
!actual building block for calculating overlap matrix
!inserted work arrays for calculation
subroutine gbasovrlp(expo1,coeff1,expo2,coeff2,ng1,ng2,l1,m1,l2,m2,dx,dy,dz,&
     niw,nrw,iw,rw,ovrlp)
  use module_base
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
  
end subroutine gbasovrlp

!overlap matrix between two different basis structures
!works only if A%ncoeff==B%ncoeff
subroutine gaussian_overlap(A,B,ovrlp)
  use module_types
  type(gaussian_basis), intent(in) :: A,B
  real(gp), dimension(A%ncoeff,B%ncoeff) :: ovrlp !only lower triangular part for A%ncoeff=B%ncoeff

  !local variables
  integer, parameter :: niw=18,nrw=6
  integer :: ishell,iexpo,icoeff,ishellB,iexpoB,icoeffB,iat,jat,isat,isatB,jsat,jshell
  integer :: jstart,iovrlp,jovrlp
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
                    if (jovrlp >= iovrlp .and. A%ncoeff == B%ncoeff) then
                       call gbasovrlp(A%xp(iexpo),A%psiat(iexpo),B%xp(jexpo),B%psiat(jexpo),&
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
  
end subroutine gaussian_overlap


!calculates a dot product between two differents gaussians times spherical harmonics
!vaild only for shell which belongs to different atoms, and with also dy/=0/=dx dz/=0
!to be rearranged when only some of them is zero
subroutine gprod(a1,a2,dx,dy,dz,l1,m1,l2,m2,niw,nrw,iw,rw,ovrlp)
  use module_base
  implicit none
  integer, intent(in) :: l1,l2,m1,m2,niw,nrw 
  real(gp), intent(in) :: a1,a2,dx,dy,dz
  integer, dimension(niw) :: iw !work array of the exponents of the two polynomials
  real(gp), dimension(nrw) :: rw !work array of the polynomials coefficients 
  real(gp), intent(out) :: ovrlp
  !local variables
  integer, parameter :: nx=3
  integer :: n1,n2,i1,i2,px,py,pz,qx,qy,qz,i
  integer :: lx1,lx2,lx3,ly1,ly2,ly3,lz1,lz2,lz3
  integer :: mx1,mx2,mx3,my1,my2,my3,mz1,mz2,mz3
  real(gp) :: fx,fy,fz,fa,fb,govrlp,f1,f2,f3,g1,g2,g3

!!$  rw(1)=0.d0
!!$  rw(2)=0.d0
!!$  rw(3)=0.d0
!!$
!!$  do i=1,3*nx
!!$     iw(i)=0
!!$  end do
  
  
  !calculates the number of different couples
  call calc_coeff_inguess(l1,m1,nx,n1,&
       iw(1),iw(nx+1),iw(2*nx+1),rw(1))
!!$  lx1=iw(1)
!!$  lx2=iw(2)
!!$  lx3=iw(3)
!!$  ly1=iw(4)
!!$  ly2=iw(5)
!!$  ly3=iw(6)
!!$  lz1=iw(7)
!!$  lz2=iw(8)
!!$  lz3=iw(9)
!!$  
!!$  f1=rw(1)
!!$  f2=rw(2)
!!$  f3=rw(3)

  call calc_coeff_inguess(l2,m2,nx,n2,&
       iw(3*nx+1),iw(4*nx+1),iw(5*nx+1),rw(n1+1))

!!$  mx1=iw(1)
!!$  mx2=iw(2)
!!$  mx3=iw(3)
!!$  my1=iw(4)
!!$  my2=iw(5)
!!$  my3=iw(6)
!!$  mz1=iw(7)
!!$  mz2=iw(8)
!!$  mz3=iw(9)
!!$  
!!$  g1=rw(1)
!!$  g2=rw(2)
!!$  g3=rw(3)
!!$
!!$  !start unrolled loop
!!$  ovrlp=0.d0
!!$
!!$  fx=govrlp(a1,a2,dx,lx1,mx1)
!!$  fy=govrlp(a1,a2,dy,ly1,my1)
!!$  fz=govrlp(a1,a2,dz,lz1,mz1)
!!$
!!$  ovrlp=ovrlp+f1*g1*fx*fy*fz
!!$
!!$  fx=govrlp(a1,a2,dx,lx1,mx2)
!!$  fy=govrlp(a1,a2,dy,ly1,my2)
!!$  fz=govrlp(a1,a2,dz,lz1,mz2)
!!$
!!$  ovrlp=ovrlp+f1*g2*fx*fy*fz
!!$
!!$  fx=govrlp(a1,a2,dx,lx1,mx3)
!!$  fy=govrlp(a1,a2,dy,ly1,my3)
!!$  fz=govrlp(a1,a2,dz,lz1,mz3)
!!$
!!$  ovrlp=ovrlp+f1*g3*fx*fy*fz
!!$
!!$  fx=govrlp(a1,a2,dx,lx2,mx1)
!!$  fy=govrlp(a1,a2,dy,ly2,my1)
!!$  fz=govrlp(a1,a2,dz,lz2,mz1)
!!$
!!$  ovrlp=ovrlp+f2*g1*fx*fy*fz
!!$
!!$  fx=govrlp(a1,a2,dx,lx2,mx2)
!!$  fy=govrlp(a1,a2,dy,ly2,my2)
!!$  fz=govrlp(a1,a2,dz,lz2,mz2)
!!$
!!$  ovrlp=ovrlp+f2*g2*fx*fy*fz
!!$
!!$  fx=govrlp(a1,a2,dx,lx2,mx3)
!!$  fy=govrlp(a1,a2,dy,ly2,my3)
!!$  fz=govrlp(a1,a2,dz,lz2,mz3)
!!$
!!$  ovrlp=ovrlp+f2*g3*fx*fy*fz
!!$
!!$  fx=govrlp(a1,a2,dx,lx3,mx1)
!!$  fy=govrlp(a1,a2,dy,ly3,my1)
!!$  fz=govrlp(a1,a2,dz,lz3,mz1)
!!$
!!$  ovrlp=ovrlp+f3*g1*fx*fy*fz
!!$
!!$  fx=govrlp(a1,a2,dx,lx3,mx2)
!!$  fy=govrlp(a1,a2,dy,ly3,my2)
!!$  fz=govrlp(a1,a2,dz,lz3,mz2)
!!$
!!$  ovrlp=ovrlp+f3*g2*fx*fy*fz
!!$
!!$  fx=govrlp(a1,a2,dx,lx3,mx3)
!!$  fy=govrlp(a1,a2,dy,ly3,my3)
!!$  fz=govrlp(a1,a2,dz,lz3,mz3)
!!$
!!$  ovrlp=ovrlp+f3*g3*fx*fy*fz

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
 
end subroutine gprod

!calculates \int \exp^{-a1*x^2} x^l1 \exp^{-a2*(x-d)^2} (x-d)^l2 dx
function govrlp(a1,a2,d,l1,l2)
  use module_base
  implicit none
  integer, intent(in) :: l1,l2
  real(gp), intent(in) :: a1,a2,d
  real(gp) :: govrlp
  !local variables
  integer :: p
  real(gp) :: prefac,rfac,stot,aeff,ceff,tt,fsum,gauint

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
end function govrlp

!calculates \int \exp^{-a*(x-c)^2} x^l dx
!this works ALSO when c/=0.d0
function gauint(a,c,l)
  use module_base
  implicit none
  integer, intent(in) :: l
  real(gp), intent(in) :: a,c
  real(gp) :: gauint
  !local variables
  real(gp), parameter :: gammaonehalf=1.772453850905516027298d0
  integer :: p
  real(gp) :: rfac,prefac,xsum,stot,fsum,tt,firstprod
  !build the prefactor
  prefac=sqrt(a)
  prefac=1.d0/prefac
  prefac=gammaonehalf*prefac

  !the first term of the sum is one
  !but we have to multiply for the prefactor
  stot=c**l

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
  
end function gauint

!calculates \int \exp^{-a*x^2} x^l dx
!this works only when l is even (if not equal to zero)
function gauint0(a,l)
  use module_base
  implicit none
  integer, intent(in) :: l
  real(gp), intent(in) :: a
  real(gp) :: gauint0
  !local variables
  real(gp), parameter :: gammaonehalf=1.772453850905516027298d0
  integer :: p
  real(gp) :: xfac,prefac,tt,firstprod
  !build the prefactor
  prefac=sqrt(a)
  prefac=1.d0/prefac
  prefac=gammaonehalf*prefac**(l+1)

  p=l/2
  tt=xfac(1,p,-0.5d0)
  !final result
  gauint0=prefac*tt
  
end function gauint0

function firstprod(p)
  use module_base
  implicit none
  integer, intent(in) :: p
  real(gp) :: firstprod
  !local variables
  integer :: i
  real(gp) :: tt
  firstprod=1.d0
  do i=1,p
     tt=real(2*i,gp)
     tt=1.d0/tt
     tt=1.d0-tt
     firstprod=firstprod*tt
  end do
end function firstprod

function rfac(is,ie)
  use module_base
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
end function rfac

!with this function n!=xfac(1,n,0.d0)
function xfac(is,ie,sh)
  use module_base
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
end function xfac


!end of the interesting part


!the same function but with integer factorials (valid ONLY if l<=18)
!not a visible improvement in speed with respect to the analogous real
function gauinti(a,c,l)
  use module_base
  implicit none
  integer, intent(in) :: l
  real(gp), intent(in) :: a,c
  real(gp) :: gauinti
  !local variables
  real(gp), parameter :: gammaonehalf=1.772453850905516027298d0
  integer :: p,ifac
  real(gp) :: prefac,xsum,stot,fsum,tt,firstprod
  !build the prefactor
  prefac=sqrt(a)
  prefac=c**l/prefac
  prefac=gammaonehalf*prefac

  !object in the sum
  xsum=a*c**2
  xsum=1.d0/xsum

  !the first term of the sum is one
  stot=1.d0

  !calculate the sum
  do p=1,l/4
     tt=real(ifac(p+1,2*p),gp)
     fsum=real(ifac(l-2*p+1,l),gp)
     fsum=fsum/tt
     tt=firstprod(p)
     fsum=fsum*tt
     fsum=fsum*xsum**p
     stot=stot+fsum
  end do
  do p=l/4+1,l/3
     tt=real(ifac(p+1,l-2*p),gp)
     fsum=real(ifac(2*p+1,l),gp)
     fsum=fsum/tt
     tt=firstprod(p)
     fsum=fsum*tt
     fsum=fsum*xsum**p
     stot=stot+fsum
  end do
  do p=l/3+1,l/2
     tt=real(ifac(l-2*p+1,p),gp)
     fsum=real(ifac(2*p+1,l),gp)
     fsum=fsum*tt
     tt=firstprod(p)
     fsum=fsum*tt
     fsum=fsum*xsum**p
     stot=stot+fsum
  end do

  !final result
  gauinti=stot*prefac
  
end function gauinti


!valid if p<l/4 AND p/=0
function secondprod1(p,l)
  use module_base
  implicit none
  integer, intent(in) :: p,l
  real(gp) :: secondprod1
  !local variables
  integer :: i
  real(gp) :: tt,part1,rfac
  part1=rfac(p+1,2*p)
  !divide by the last value
  part1=real(l,gp)/part1
  tt=rfac(l-2*p+1,l-1)
!!$  part1=1.d0
!!$  do i=p+1,2*p !in the second case the bound must be changed here
!!$     tt=real(i,gp)
!!$     part1=part1*tt
!!$  end do
  secondprod1=tt*part1
end function secondprod1

!valid if p>=l/4 AND p<l/3
function secondprod2(p,l)
  use module_base
  implicit none
  integer, intent(in) :: p,l
  real(gp) :: secondprod2
  !local variables
  integer :: i
  real(gp) :: tt,part1,rfac
  part1=rfac(p+1,l-2*p)
  !divide by the last value
  part1=real(l,gp)/part1
  tt=rfac(2*p+1,l-1)
!!$  part1=1.d0
!!$  do i=p+1,2*p !in the second case the bound must be changed here
!!$     tt=real(i,gp)
!!$     part1=part1*tt
!!$  end do
  secondprod2=tt*part1
end function secondprod2


!integer version of factorial
function ifac(is,ie)
  implicit none
  integer, intent(in) :: is,ie
  integer :: ifac
  !local variables
  integer :: i
  ifac=1
  do i=is,ie
     ifac=ifac*i
  end do
end function ifac

!calculate the projection of norb wavefunctions on a gaussian basis set
subroutine wavelets_to_gaussians(geocode,norbp,n1,n2,n3,G,thetaphi,hx,hy,hz,wfd,psi,coeffs)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: norbp,n1,n2,n3
  real(gp), intent(in) :: hx,hy,hz
  type(gaussian_basis), intent(in) :: G
  type(wavefunctions_descriptors), intent(in) :: wfd
  real(gp), dimension(2,G%nat), intent(in) :: thetaphi
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
  real(wp), dimension(G%ncoeff,norbp), intent(out) :: coeffs
  !local variables
  integer :: iorb
  
  do iorb=1,norbp
     call orbital_projection(geocode,n1,n2,n3,G%nat,G%rxyz,thetaphi,&
          G%nshell,G%ndoc,G%nam,G%xp,G%psiat,G%nshltot,G%nexpo,G%ncoeff,&
          hx,hy,hz,wfd,psi(1,iorb),coeffs(1,iorb))
  end do
  
end subroutine wavelets_to_gaussians

!calculate the projection of a given orbital on a gaussian basis centered on 
!a set of points
subroutine orbital_projection(geocode,n1,n2,n3,nat,rxyz,thetaphi,nshell,ndoc,nam,xp,psiat,nshltot,nexpo,ncoeff,&
     hx,hy,hz,wfd,psi,coeffs)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n1,n2,n3,nat,nshltot,nexpo,ncoeff
  real(gp), intent(in) :: hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd
  integer, dimension(nat), intent(in) :: nshell
  integer, dimension(nshltot), intent(in) :: ndoc,nam
  real(gp), dimension(nexpo), intent(in) :: xp,psiat
  real(gp), dimension(2,nat), intent(in) :: thetaphi
  real(gp), dimension(3,nat), intent(in) :: rxyz
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(in) :: psi
  real(wp), dimension(ncoeff), intent(out) :: coeffs
  !local variables
  integer :: ishell,iexpo,icoeff,iat,isat,ng,l
  
  ishell=0
  iexpo=1
  icoeff=1
  do iat=1,nat
     do isat=1,nshell(iat)
        ishell=ishell+1
        ng=ndoc(ishell)
        l=nam(ishell)
        call lsh_projection(geocode,l,ng,xp(iexpo),psiat(iexpo),n1,n2,n3,rxyz(1,iat),&
             thetaphi(1,iat),hx,hy,hz,wfd,psi,coeffs(icoeff))
        iexpo=iexpo+ng
        icoeff=icoeff+2*l-1
     end do
  end do

  call gaudim_check(iexpo,icoeff,ishell,nexpo,ncoeff,nshltot)
  
end subroutine orbital_projection

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
end subroutine gaudim_check


!calculate the projection of a gaussian for a given eigenspace of spherical harmonics
!centered on a given point and rotated by theta(along y) and phi(along x)
subroutine lsh_projection(geocode,l,ng,xp,psiat,n1,n2,n3,rxyz,thetaphi,hx,hy,hz,&
     wfd,psi,coeffs)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: l,ng,n1,n2,n3
  real(gp), intent(in) :: hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd
  real(gp), dimension(2), intent(in) :: thetaphi
  real(gp), dimension(3), intent(in) :: rxyz
  real(gp), dimension(ng), intent(in) :: xp,psiat
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(in) :: psi
  real(wp), dimension(2*l-1), intent(out) :: coeffs
  !local variables
  integer, parameter :: nterm_max=3 !to be enlarged for generic theta and phi
  integer :: nterm,m
  integer, dimension(nterm_max) :: lx,ly,lz
  real(gp), dimension(nterm_max) :: fac_arr
  
  do m=1,2*l-1
     !for the moment theta and phi are ignored but a new routine should be made
     call calc_coeff_inguess(l,m,nterm_max,nterm,lx,ly,lz,fac_arr)

     call wavetogau(geocode,n1,n2,n3,ng,nterm,lx,ly,lz,fac_arr,xp,psiat,&
          rxyz(1),rxyz(2),rxyz(3),hx,hy,hz,wfd%nseg_c,wfd%nvctr_c,wfd%keyg(1,1),wfd%keyv(1),&
          wfd%nseg_f,wfd%nvctr_f,wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1),&
          psi(1),psi(wfd%nvctr_c+1),coeffs(m))

     !here we should modify the coefficients following the direction of the axis

  end do

end subroutine lsh_projection




! calculate the scalar product between a sum of gaussians times polynomials and a wavefunction
! \int dx dy dz 
!           \sum_i=1..ntp fac_arr(i) {
!                \sum_j=1..nterm psiat(j) [exp(-r^2/(2*(xp(j)^2)))] 
!                    *((x-rx)^lx(i) *(y-ry)^ly(i) * (z-rz)^lz(i) ))} psi(x,y,z)
! expressed in Daubechies Basis in the arrays psi_c, psi_f
subroutine wavetogau(geocode,n1,n2,n3,nterm,ntp,lx,ly,lz,fac_arr,xp,psiat,rx,ry,rz,hx,hy,hz, & 
     nseg_c,mvctr_c,keyg_c,keyv_c,nseg_f,mvctr_f,keyg_f,keyv_f,psi_c,psi_f,overlap)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n1,n2,n3,nterm,ntp,nseg_c,nseg_f,mvctr_c,mvctr_f
  real(gp), intent(in) :: rx,ry,rz,hx,hy,hz
  integer, dimension(ntp), intent(in) :: lx,ly,lz
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(gp), dimension(ntp), intent(in) :: fac_arr
  real(gp), dimension(nterm), intent(in) :: xp,psiat
  real(wp), dimension(mvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(in) :: psi_f
  real(wp), intent(out) :: overlap
  !local variables
  character(len=*), parameter :: subname='wavetogau'
  integer, parameter ::nw=16000
  logical :: perx,pery,perz
  integer:: iterm,itp,n_gau,ml1,mu1,ml2,mu2,ml3,mu3,i1,i2,i3,i_all,i_stat,iseg,ii,jj,j0,j1,i0,i
  real(wp) :: ovlp_c,ovlp_f1,ovlp_f2,ovlp_f3,ovlp_f4,ovlp_f5,ovlp_f6,ovlp_f7,ovlp
  real(gp) :: gau_a,te
  real(wp), dimension(:,:), allocatable :: work,wprojx,wprojy,wprojz

  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')


  allocate(wprojx(0:n1,2+ndebug),stat=i_stat)
  call memocc(i_stat,wprojx,'wprojx',subname)
  allocate(wprojy(0:n2,2+ndebug),stat=i_stat)
  call memocc(i_stat,wprojy,'wprojy',subname)
  allocate(wprojz(0:n3,2+ndebug),stat=i_stat)
  call memocc(i_stat,wprojz,'wprojz',subname)
  allocate(work(0:nw,2+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)

  overlap=0.0_wp

  do itp=1,ntp

     do iterm=1,nterm
        gau_a=xp(iterm)
        n_gau=lx(itp)
        call gauss_to_daub(hx,fac_arr(itp),rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw,&
             perx)
        n_gau=ly(itp)
        call gauss_to_daub(hy,1.0_gp,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw,pery)
        n_gau=lz(itp)
        call gauss_to_daub(hz,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw,&
             perz)

        !scalar product (in wavefunction precision)

        ! coarse part
        ovlp_c=0.0_wp
        do iseg=1,nseg_c
           jj=keyv_c(iseg)
           j0=keyg_c(1,iseg)
           j1=keyg_c(2,iseg)
           ii=j0-1
           i3=ii/((n1+1)*(n2+1))
           ii=ii-i3*(n1+1)*(n2+1)
           i2=ii/(n1+1)
           i0=ii-i2*(n1+1)
           i1=i0+j1-j0
           do i=i0,i1
              ovlp_c=ovlp_c+psi_c(i-i0+jj)*wprojx(i,1)*wprojy(i2,1)*wprojz(i3,1)
           enddo
        enddo

        ! fine part
        ovlp_f1=0.0_wp
        ovlp_f2=0.0_wp
        ovlp_f3=0.0_wp
        ovlp_f4=0.0_wp
        ovlp_f5=0.0_wp
        ovlp_f6=0.0_wp
        ovlp_f7=0.0_wp
        do iseg=1,nseg_f
           jj=keyv_f(iseg)
           j0=keyg_f(1,iseg)
           j1=keyg_f(2,iseg)
           ii=j0-1
           i3=ii/((n1+1)*(n2+1))
           ii=ii-i3*(n1+1)*(n2+1)
           i2=ii/(n1+1)
           i0=ii-i2*(n1+1)
           i1=i0+j1-j0
           do i=i0,i1
              !itp=itp+1
              ovlp_f1=ovlp_f1+psi_f(1,i-i0+jj)*wprojx(i,2)*wprojy(i2,1)*wprojz(i3,1)
              ovlp_f2=ovlp_f2+psi_f(2,i-i0+jj)*wprojx(i,1)*wprojy(i2,2)*wprojz(i3,1)
              ovlp_f3=ovlp_f3+psi_f(3,i-i0+jj)*wprojx(i,2)*wprojy(i2,2)*wprojz(i3,1)
              ovlp_f4=ovlp_f4+psi_f(4,i-i0+jj)*wprojx(i,1)*wprojy(i2,1)*wprojz(i3,2)
              ovlp_f5=ovlp_f5+psi_f(5,i-i0+jj)*wprojx(i,2)*wprojy(i2,1)*wprojz(i3,2)
              ovlp_f6=ovlp_f6+psi_f(6,i-i0+jj)*wprojx(i,1)*wprojy(i2,2)*wprojz(i3,2)
              ovlp_f7=ovlp_f7+psi_f(7,i-i0+jj)*wprojx(i,2)*wprojy(i2,2)*wprojz(i3,2)
           enddo
        enddo

        ovlp=ovlp_c+ovlp_f1+ovlp_f2+ovlp_f3+ovlp_f4+ovlp_f5+ovlp_f6+ovlp_f7

        overlap=overlap+ovlp

     end do


  end do

  i_all=-product(shape(wprojx))*kind(wprojx)
  deallocate(wprojx,stat=i_stat)
  call memocc(i_stat,i_all,'wprojx',subname)
  i_all=-product(shape(wprojy))*kind(wprojy)
  deallocate(wprojy,stat=i_stat)
  call memocc(i_stat,i_all,'wprojy',subname)
  i_all=-product(shape(wprojz))*kind(wprojz)
  deallocate(wprojz,stat=i_stat)
  call memocc(i_stat,i_all,'wprojz',subname)
  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)

end subroutine wavetogau
