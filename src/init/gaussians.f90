!overlap matrix between two different basis structures
subroutine gaussian_overlap(A,B,ovrlp)
  use module_base
  use module_types
  implicit none
  type(gaussian_basis), intent(in) :: A,B
  real(gp), dimension(A%ncoeff,B%ncoeff) :: ovrlp 
  !only lower triangular part for A%ncoeff=B%ncoeff
  !local variables
  integer, parameter :: niw=18,nrw=6
  integer :: ishell,iexpo,icoeff,ishellB,iexpoB,icoeffB,iat,jat,isat,isatB,jsat,jshell
  integer :: jstart,iovrlp,jovrlp,jcoeff,jexpo
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
                       call gbasovrlp(A%xp(iexpo),A%psiat(iexpo),&
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
  
end subroutine gaussian_overlap

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


!overlap kinetic matrix between two different basis structures
!the kinetic operator is applicated on the A basis structure
subroutine kinetic_overlap(A,B,ovrlp)
  use module_base
  use module_types
  implicit none
  type(gaussian_basis), intent(in) :: A,B
  real(gp), dimension(A%ncoeff,B%ncoeff) :: ovrlp 
  !only lower triangular part for A%ncoeff=B%ncoeff
  !local variables
  integer, parameter :: niw=18,nrw=6
  integer :: ishell,iexpo,icoeff,ishellB,iexpoB,icoeffB,iat,jat,isat,isatB,jsat,jshell
  integer :: jstart,iovrlp,jovrlp,jcoeff,jexpo
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
  
end subroutine kinetic_overlap

!calculates the scalar product between two shells
!by considering only the nonzero coefficients
!actual building block for calculating overlap matrix
!inserted work arrays for calculation
subroutine kineticovrlp(expo1,coeff1,expo2,coeff2,ng1,ng2,l1,m1,l2,m2,dx,dy,dz,&
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
        call kinprod(a1,a2,dx,dy,dz,l1,m1,l2,m2,niw,nrw,iw,rw,govrlpr)
        govrlpr=c1*govrlpr*c2
        !print *,c1,c2,govrlpr
        ovrlp=ovrlp+govrlpr
     end do
  end do
  
end subroutine kineticovrlp


!overlap kinetic matrix between two different basis structures
!the kinetic operator is applicated on the A basis structure
subroutine potential_overlap(A,B,pot,n1,n2,n3,hx,hy,hz,ovrlp)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(gp), intent(in) :: hx,hy,hz
  type(gaussian_basis), intent(in) :: A,B
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: pot
  real(gp), dimension(A%ncoeff,B%ncoeff), intent(out) :: ovrlp 
  !only lower triangular part for A%ncoeff=B%ncoeff
  !local variables 
 integer, parameter :: niw=18,nrw=6
  integer :: ishell,iexpo,icoeff,ishellB,iexpoB,icoeffB,iat,jat,isat,isatB,jsat,jshell
  integer :: jstart,iovrlp,jovrlp,jcoeff,jexpo
  integer :: ngA,ngB,lA,lB,mA,mB
  real(gp) :: rxa,rya,rza,rxb,ryb,rzb
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

           !here one may insert jat=iat if gaussians do not overlap
           do jat=1,B%nat
              rxa=A%rxyz(1,iat)
              rya=A%rxyz(2,iat)
              rza=A%rxyz(3,iat)

              rxb=B%rxyz(1,jat)
              ryb=B%rxyz(2,jat)
              rzb=B%rxyz(3,jat)

              do jsat=1,B%nshell(jat)
                 jshell=jshell+1
                 ngB=B%ndoc(jshell)
                 lB=B%nam(jshell)
                 do mB=1,2*lB-1
                    jovrlp=jovrlp+1
                    if (jovrlp >= iovrlp .and. A%ncoeff == B%ncoeff) then
                       
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
  
end subroutine potential_overlap

subroutine locpotovrlp(n1,n2,n3,pot,expo1,coeff1,expo2,coeff2,ng1,ng2,l1,m1,l2,m2,&
     rxa,rya,rza,rxb,ryb,rzb,niw,nrw,iw,rw,ovrlp)
  use module_base
  implicit none
  integer, intent(in) :: ng1,ng2,l1,m1,l2,m2,niw,nrw
  real(gp), intent(in) :: rxa,rya,rza,rxb,ryb,rzb
  integer, dimension(niw) :: iw
  real(gp), dimension(nrw) :: rw
  real(gp), dimension(ng1), intent(in) :: expo1,coeff1
  real(gp), dimension(ng2), intent(in) :: expo2,coeff2
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: pot
  real(gp), intent(out) :: ovrlp
  !local variables
  integer :: ii1,ii2,j1,j2,j3,i1,i2,i3,isx,isy,isz,iex,iey,iez
  real(gp) :: a1,a2,c1,c2,rmean,xmean,ymean,zmean,prodgau,cutoff,xa,ya,za,xb,yb,zb


  !calculates the polynomials which multiply each product of gaussians
  call calc_coeff_inguess(l1,m1,nx,n1,&
       iw(1),iw(nx+1),iw(2*nx+1),rw(1))
  call calc_coeff_inguess(l2,m2,nx,n2,&
       iw(3*nx+1),iw(4*nx+1),iw(5*nx+1),rw(n1+1))


  ovrlp=0.d0
  do ii1=1,ng1
     a1=expo1(ii1)
     a1=0.5_gp/a1**2
     c1=coeff1(ii1)
     do ii2=1,ng2
        a2=expo2(ii2)
        a2=0.5_gp/a2**2
        c2=coeff2(ii2)
        !calculate overall factor given by product of gaussian
        nexp=a1*a2/(a1+a2)
        expx=nexp*(rxa-rxb)**2
        expy=nexp*(rya-ryb)**2
        expz=nexp*(rza-rzb)**2
        factor=c1*c2*exp(-expx-expy-expz)
        ovrlp=0.0_gp
        if (factor > 1.e-8_gp) then

           xmean=rmean(a1,a2,rxa,rxb)
           ymean=rmean(a1,a2,rya,ryb)
           zmean=rmean(a1,a2,rza,rzb)
           cutoff=5._gp*sqrt(0.5_gp/(a1+a2))
           !limits for integration of the potential
           isx=floor((rx-cutoff)/hx)
           isy=floor((ry-cutoff)/hy)
           isz=floor((rz-cutoff)/hz)

           iex=ceiling((rx+cutoff)/hx)
           iey=ceiling((ry+cutoff)/hy)
           iez=ceiling((rz+cutoff)/hz)
          
           povrlp=0.0_gp
           do i3=isz,iez
              call ind_gauss(.false.,i3,0,n3,j3,goz) 
              if (goz) then
                 z=real(i3,gp)*hz-zmean
                 za=real(i3,gp)*hz-rza
                 zb=real(i3,gp)*hz-rzb
                 prodgaus=exp(-(a1+a2)*z**2)
                 do i2=isy,iey
                    call ind_gauss(.false.,i2,0,n2,j2,goy)
                    if (goy) then
                       y=real(i2,gp)*hy-ymean
                       ya=real(i2,gp)*hy-rya
                       yb=real(i2,gp)*hy-ryb
                       prodgaus=prodgaus*exp(-(a1+a2)*y**2)
                       do i1=isx,iex
                          call ind_gauss(.false.,i1,0,n1,j1,gox)
                          if (gox) then
                             x=real(i1,gp)*hx-zmean
                             xa=real(i1,gp)*hx-rxa
                             xb=real(i1,gp)*hx-rxb
                             prodgaus=prodgaus*exp(-(a1+a2)*x**2)
                             polb=0.0_gp
                             do i2=1,n2
                                qx=iw(3*nx+i2)
                                qy=iw(4*nx+i2)
                                qz=iw(5*nx+i2)
                                fb=rw(n1+i2)
                                polb=polb+fb*(xb**qx)*(yb**qy)*(zb**qz)
                             end do
                             pola=0.0_gp
                             do i1=1,n1
                                px=iw(i1)
                                py=iw(nx+i1)
                                pz=iw(2*nx+i1)
                                fa=rw(i1)
                                pola=pola+fa*(xa**px)*(ya**py)*(za**pz)
                             end do
                             prodgaus=prodgaus*polx*poly
                                end do
                             end do
                             povrlp=povrlp+pot(j1,j2,j3)*prodgaus
                          end if
                       enddo
                    end if
                 enddo
              end if
           end do
           ovrlp=factor*ovrlp
        end if
     end do
  end do
  
end subroutine locpotovrlp

function rmean(a1,a2,r1,r2)
  use module_base
  implicit none
  real(gp), intent(in) :: a1,a2,r1,r2
  real(gp) :: rmean
  
  rmean=a1*r1+a2*r2
  rmean=rmean/(a1+a2)

end function rmean

subroutine ind_gauss(periodic,i,is,n,j,go)
  implicit none
  logical, intent(in) :: periodic
  integer, intent(in) :: i,n,is
  logical, intent(out) :: go
  integer, intent(out) :: j

  if (periodic .and. is == 0) then
     go=.true.
     j=modulo(i,n+1)
  else
     j=i
     if (i >= is .and. i <= n1+is) then
        go=.true.
     else
        go=.false.
     end if
  end if

end subroutine ind_gauss

