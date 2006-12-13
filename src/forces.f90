subroutine local_forces(iproc,nproc,ntypes,nat,iatype,atomnames,rxyz,psppar,nelpsp,hgrid,&
     n1,n2,n3,rho,pot,floc)
! Calculates the local forces acting on the atoms belonging to iproc
  use libBigDFT
  
  implicit none
  !Arguments---------
  integer, intent(in) :: iproc,nproc,ntypes,nat,n1,n2,n3
  real(kind=8), intent(in) :: hgrid
  character*20, dimension(100), intent(in) :: atomnames
  real(kind=8), dimension(0:2,0:4,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(in) :: rho,pot
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(kind=8), dimension(3,nat), intent(out) :: floc
  !Local variables---------
  real(kind=8) :: hgridh,pi,prefactor,cutoff,rloc,Vel,rhoel
  real(kind=8) :: fxerf,fyerf,fzerf,fxion,fyion,fzion,fxgau,fygau,fzgau,forceleaked,forceloc
  real(kind=8) :: rx,ry,rz,x,y,z,arg,r2,xp,dist,tt
  integer :: ii,ix,iy,iz,i1,i2,i3,iat,jat,ityp,jtyp,nloc,iloc,istat
  !array of coefficients of the derivative
  real(kind=8), dimension(4) :: cprime 

  hgridh=hgrid*.5d0 
  pi=4.d0*atan(1.d0)



  do iat=1,nat
  if (mod(iat-1,nproc).eq.iproc) then
  write(*,*) iproc,' calculates local force on atom ',iat
     ityp=iatype(iat)
     !coordinates of the center
     rx=rxyz(1,iat) 
     ry=rxyz(2,iat) 
     rz=rxyz(3,iat)
     !nearest grid points to the center
     ix=nint(rx/hgridh)
     iy=nint(ry/hgridh)
     iz=nint(rz/hgridh)
     !inizialization of the forces
     !ion-ion term
     fxion=0.d0
     fyion=0.d0
     fzion=0.d0
     !ion-electron term, error function part
     fxerf=0.d0
     fyerf=0.d0
     fzerf=0.d0
     !ion-electron term, gaussian part
     fxgau=0.d0
     fygau=0.d0
     fzgau=0.d0

     !Derivative of the ion-ion energy
     do jat=1,iat-1
        dist=sqrt((rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2)
        jtyp=iatype(jat)
        !eion=eion+nelpsp(jtyp)*nelpsp(ityp)/dist
        fxion=fxion+nelpsp(jtyp)*(nelpsp(ityp)/(dist**3))*(rx-rxyz(1,jat))
        fyion=fyion+nelpsp(jtyp)*(nelpsp(ityp)/(dist**3))*(ry-rxyz(2,jat))
        fzion=fzion+nelpsp(jtyp)*(nelpsp(ityp)/(dist**3))*(rz-rxyz(3,jat))
     end do
     do jat=iat+1,nat
        dist=sqrt((rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2)
        jtyp=iatype(jat)
        fxion=fxion+nelpsp(jtyp)*(nelpsp(ityp)/(dist**3))*(rx-rxyz(1,jat))
        fyion=fyion+nelpsp(jtyp)*(nelpsp(ityp)/(dist**3))*(ry-rxyz(2,jat))
        fzion=fzion+nelpsp(jtyp)*(nelpsp(ityp)/(dist**3))*(rz-rxyz(3,jat))
     end do

     
     !building array of coefficients of the derivative of the gaussian part
     cprime(1)=2.d0*psppar(0,2,ityp)-psppar(0,1,ityp)
     cprime(2)=4.d0*psppar(0,3,ityp)-psppar(0,2,ityp)
     cprime(3)=6.d0*psppar(0,4,ityp)-psppar(0,3,ityp)
     cprime(4)=-psppar(0,4,ityp)

     ! determine number of local terms
     nloc=0
     do iloc=1,4
        if (psppar(0,iloc,ityp).ne.0.d0) nloc=iloc
     enddo

     !local part
     forceleaked=0.d0
     rloc=psppar(0,0,ityp)
     prefactor=nelpsp(ityp)/(2.d0*pi*sqrt(2.d0*pi)*rloc**5)
     !maximum extension of the gaussian
     cutoff=10.d0*rloc
     !nearest grid point to the cutoff
     ii=nint(cutoff/hgridh)
     !calculate the forces near the atom due to the error function part of the potential
     do i3=iz-ii,iz+ii
        do i2=iy-ii,iy+ii
           do i1=ix-ii,ix+ii
              x=i1*hgridh-rx
              y=i2*hgridh-ry
              z=i3*hgridh-rz
              r2=x**2+y**2+z**2
              arg=r2/rloc**2
              xp=exp(-.5d0*arg)
              if (i3.ge.-14 .and. i3.le.2*n3+16  .and.  & 
                   i2.ge.-14 .and. i2.le.2*n2+16  .and.  & 
                   i1.ge.-14 .and. i1.le.2*n1+16 ) then
                 !gaussian part
                 if (nloc /= 0) then
                    tt=cprime(nloc)
                    do iloc=nloc-1,1,-1
                       tt=arg*tt+cprime(iloc)
                    enddo
                    rhoel=rho(i1,i2,i3)
                    forceloc=xp*tt*rhoel
                    fxgau=fxgau+forceloc*x
                    fygau=fygau+forceloc*y
                    fzgau=fzgau+forceloc*z
                 end if
                 !error function part
                 Vel=pot(i1,i2,i3)
                 fxerf=fxerf+xp*Vel*x
                 fyerf=fyerf+xp*Vel*y
                 fzerf=fzerf+xp*Vel*z
              else
                 forceleaked=forceleaked+xp*(1+tt)
              endif
           end do
        end do
     end do

     !final result of the forces
     
    floc(1,iat)=fxion+(hgridh**3*prefactor)*fxerf+(hgridh**3/rloc**2)*fxgau
    floc(2,iat)=fyion+(hgridh**3*prefactor)*fyerf+(hgridh**3/rloc**2)*fygau
    floc(3,iat)=fzion+(hgridh**3*prefactor)*fzerf+(hgridh**3/rloc**2)*fzgau

  else

    floc(1,iat)=0.d0
    floc(2,iat)=0.d0
    floc(3,iat)=0.d0

  endif
  end do

  forceleaked=forceleaked*prefactor*hgridh**3
  if (iproc.eq.0) write(*,'(a,e21.14,1x,e10.3)') 'leaked force: ',forceleaked

end subroutine local_forces

subroutine nonlocal_forces(iproc,nproc,n1,n2,n3,nboxp_c,nboxp_f, & 
     ntypes,nat,norb,norb_p,istart,nprojel,nproj,&
     iatype,psppar,occup,nseg_c,nseg_f,nvctr_c,nvctr_f,nseg_p,nvctr_p,proj,  &
     keyg,keyv,keyg_p,keyv_p,psi,rxyz,radii_cf,cpmult,fpmult,hgrid,fsep)
!Calculates the nonlocal forces on all atoms arising from the wavefunctions belonging to iproc and ads them to the force array
  use libBigDFT
  
  implicit none
  !Arguments-------------
  integer, intent(in) :: iproc,nproc,ntypes,nat,norb,norb_p,istart,nprojel,nproj,nseg_c,nseg_f,nvctr_c,nvctr_f
  integer, intent(in) :: n1,n2,n3
  real(kind=8),intent(in) :: cpmult,fpmult,hgrid 
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(0:2*nat), intent(in) :: nseg_p,nvctr_p
  integer, dimension(2,3,nat), intent(in) :: nboxp_c,nboxp_f
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_p(2*nat)), intent(in) :: keyg_p
  integer, dimension(nseg_p(2*nat)), intent(in) :: keyv_p
  real(kind=8), dimension(0:2,0:4,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(ntypes,2), intent(in) :: radii_cf
  real(kind=8), dimension(norb), intent(in) :: occup
  real(kind=8), dimension(nprojel), intent(in) :: proj
  real(kind=8), dimension(nvctr_c+7*nvctr_f,norb_p), intent(in) :: psi
  real(kind=8), dimension(3,nat), intent(out) :: fsep
  !Local Variables--------------
  real(kind=8), parameter :: eps_mach=1.d-12
  integer, dimension(3) :: lx,ly,lz
  real(kind=8), dimension(:,:), allocatable :: derproj,fxyz_orb
  real(kind=8), dimension(:), allocatable :: auxproj_c,auxproj_f
  integer :: istart_c,istart_f,iproj,iat,ityp,i,j,l,nterm
  integer :: mvctr_c,mvctr_f,mbseg_c,mbseg_f,maseg_c,maseg_f,iseg_c,iseg_f,jseg_c,jseg_f
  integer :: nl1_c,nl2_c,nl3_c,nl1_f,nl2_f,nl3_f,nu1_c,nu2_c,nu3_c,nu1_f,nu2_f,nu3_f
  integer :: mavctr_c,mavctr_f,mbvctr_c,mbvctr_f,ipsi_c,ipsi_f,iorb,i_c,i_f
  real(kind=8) :: fpi,factor,gau_a,onem,scpr,scprp,tcprx,tcpry,tcprz,rx,ry,rz,eproj,fx,fy,fz

  allocate(derproj(nprojel,3))

  !create the derivative of the projectors
  istart_c=1
  iproj=0
  fpi=(4.d0*atan(1.d0))**(-.75d0)
  do iat=1,nat
     rx=rxyz(1,iat)
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)
     ityp=iatype(iat)
     
     mvctr_c=nvctr_p(2*iat-1)-nvctr_p(2*iat-2)
     mvctr_f=nvctr_p(2*iat  )-nvctr_p(2*iat-1)

     nl1_c=nboxp_c(1,1,iat) ; nu1_c=nboxp_c(2,1,iat)
     nl2_c=nboxp_c(1,2,iat) ; nu2_c=nboxp_c(2,2,iat)
     nl3_c=nboxp_c(1,3,iat) ; nu3_c=nboxp_c(2,3,iat)
     nl1_f=nboxp_f(1,1,iat) ; nu1_f=nboxp_f(2,1,iat)
     nl2_f=nboxp_f(1,2,iat) ; nu2_f=nboxp_f(2,2,iat)
     nl3_f=nboxp_f(1,3,iat) ; nu3_f=nboxp_f(2,3,iat)

     !allocation of the auxiliary arrays
     !in case of need they can be allocated only in the 
     !useful case (when we have the 2s or the 1p projectors)
     allocate(auxproj_c(mvctr_c),auxproj_f(7*mvctr_f))


     ! ONLY GTH PSP form (not HGH)
     do i=1,2
        do j=1,2
           if (psppar(i,j,ityp).ne.0.d0) then
              gau_a=psppar(i,0,ityp)
              do l=1,2*i-1

                 istart_f=istart_c+mvctr_c

                 if (i.eq.1 .and. j.eq.1) then    ! first s type projector
                    factor=-fpi/(sqrt(gau_a)**3)/gau_a**2

                    !derivative wrt x direction
                    nterm=1
                    lx(1)=1 ; ly(1)=0 ; lz(1)=0 
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,derproj(istart_c,1),derproj(istart_f,1))

                    !derivative wrt y direction
                    nterm=1
                    lx(1)=0 ; ly(1)=1 ; lz(1)=0 
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,derproj(istart_c,2),derproj(istart_f,2))

                    !derivative wrt z direction
                    nterm=1
                    lx(1)=0 ; ly(1)=0 ; lz(1)=1 
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,derproj(istart_c,3),derproj(istart_f,3))
                    
                 else if (i.eq.1 .and. j.eq.2) then   ! second s type projector
                    !derivative wrt x direction (first part)
                    nterm=3 
                    lx(1)=3 ; ly(1)=0 ; lz(1)=0 
                    lx(2)=1 ; ly(2)=2 ; lz(2)=0 
                    lx(3)=1 ; ly(3)=0 ; lz(3)=2 
                    factor=-sqrt(4.d0/15.d0)*fpi/(sqrt(gau_a)**7)/gau_a**2
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,derproj(istart_c,1),derproj(istart_f,1))
                    !derivative wrt x direction (second part)
                    nterm=1 
                    lx(1)=1 ; ly(1)=0 ; lz(1)=0 
                    factor=2.d0*sqrt(4.d0/15.d0)*fpi/(sqrt(gau_a)**7)
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,auxproj_c,auxproj_f)
                    !assemble the two parts
                    do i_c=1,mvctr_c
                       derproj(istart_c-1+i_c,1)=&
                            derproj(istart_c-1+i_c,1)+auxproj_c(i_c)
                    end do
                    do i_f=1,7*mvctr_f
                       derproj(istart_f-1+i_f,1)=&
                            derproj(istart_f-1+i_c,1)+auxproj_f(i_f)
                    end do

                    !derivative wrt y direction (first part)
                    nterm=3 
                    lx(1)=2 ; ly(1)=1 ; lz(1)=0 
                    lx(2)=0 ; ly(2)=3 ; lz(2)=0 
                    lx(3)=0 ; ly(3)=1 ; lz(3)=2 
                    factor=-sqrt(4.d0/15.d0)*fpi/(sqrt(gau_a)**7)/gau_a**2
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,derproj(istart_c,2),derproj(istart_f,2))
                    !derivative wrt y direction (second part)
                    nterm=1 
                    lx(1)=0 ; ly(1)=1 ; lz(1)=0 
                    factor=2.d0*sqrt(4.d0/15.d0)*fpi/(sqrt(gau_a)**7)
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,auxproj_c,auxproj_f)
                    !assemble the two parts
                    do i_c=1,mvctr_c
                       derproj(istart_c-1+i_c,2)=&
                            derproj(istart_c-1+i_c,2)+auxproj_c(i_c)
                    end do
                    do i_f=1,7*mvctr_f
                       derproj(istart_f-1+i_f,1)=&
                            derproj(istart_f-1+i_c,2)+auxproj_f(i_f)
                    end do

                    !derivative wrt z direction (first part)
                    nterm=3 
                    lx(1)=2 ; ly(1)=0 ; lz(1)=1 
                    lx(2)=0 ; ly(2)=2 ; lz(2)=1 
                    lx(3)=0 ; ly(3)=0 ; lz(3)=3 
                    factor=-sqrt(4.d0/15.d0)*fpi/(sqrt(gau_a)**7)/gau_a**2
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,derproj(istart_c,3),derproj(istart_f,3))
                    !derivative wrt z direction (second part)
                    nterm=1 
                    lx(1)=0 ; ly(1)=0 ; lz(1)=1 
                    factor=2.d0*sqrt(4.d0/15.d0)*fpi/(sqrt(gau_a)**7)
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,auxproj_c,auxproj_f)
                    !assemble the two parts
                    do i_c=1,mvctr_c
                       derproj(istart_c-1+i_c,3)=&
                            derproj(istart_c-1+i_c,3)+auxproj_c(i_c)
                    end do
                    do i_f=1,7*mvctr_f
                       derproj(istart_f-1+i_f,3)=&
                            derproj(istart_f-1+i_c,3)+auxproj_f(i_f)
                    end do


                 else if (i.eq.2 .and. l.eq.1) then  ! px type projector
                    !derivative wrt x direction (first part)
                    nterm=1 
                    lx(1)=2 ; ly(1)=0 ; lz(1)=0 
                    factor=-sqrt(2.d0)*fpi/(sqrt(gau_a)**5)/gau_a**2
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,derproj(istart_c,1),derproj(istart_f,1))
                    !derivative wrt x direction (second part)
                    nterm=1 
                    lx(1)=0 ; ly(1)=0 ; lz(1)=0 
                    factor=sqrt(2.d0)*fpi/(sqrt(gau_a)**5)
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,auxproj_c,auxproj_f)
                    !assemble the two parts
                    do i_c=1,mvctr_c
                       derproj(istart_c-1+i_c,1)=&
                            derproj(istart_c-1+i_c,1)+auxproj_c(i_c)
                    end do
                    do i_f=1,7*mvctr_f
                       derproj(istart_f-1+i_f,1)=&
                            derproj(istart_f-1+i_c,1)+auxproj_f(i_f)
                    end do

                    !derivative wrt y direction
                    nterm=1
                    lx(1)=1 ; ly(1)=1 ; lz(1)=0
                    factor=-sqrt(2.d0)*fpi/(sqrt(gau_a)**5)/gau_a**2
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,derproj(istart_c,2),derproj(istart_f,2))

                    !derivative wrt z direction
                    nterm=1
                    lx(1)=1 ; ly(1)=0 ; lz(1)=1
                    factor=-sqrt(2.d0)*fpi/(sqrt(gau_a)**5)/gau_a**2
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,derproj(istart_c,3),derproj(istart_f,3))


                 else if (i.eq.2 .and. l.eq.2) then  ! py type projector

                    !derivative wrt x direction
                    nterm=1
                    lx(1)=1 ; ly(1)=1 ; lz(1)=0
                    factor=-sqrt(2.d0)*fpi/(sqrt(gau_a)**5)/gau_a**2
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,derproj(istart_c,1),derproj(istart_f,1))


                    !derivative wrt y direction (first part)
                    nterm=1 
                    lx(1)=0 ; ly(1)=2 ; lz(1)=0 
                    factor=-sqrt(2.d0)*fpi/(sqrt(gau_a)**5)/gau_a**2
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,derproj(istart_c,2),derproj(istart_f,2))
                    !derivative wrt y direction (second part)
                    nterm=1 
                    lx(1)=0 ; ly(1)=0 ; lz(1)=0 
                    factor=sqrt(2.d0)*fpi/(sqrt(gau_a)**5)
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,auxproj_c,auxproj_f)
                    !assemble the two parts
                    do i_c=1,mvctr_c
                       derproj(istart_c-1+i_c,2)=&
                            derproj(istart_c-1+i_c,2)+auxproj_c(i_c)
                    end do
                    do i_f=1,7*mvctr_f
                       derproj(istart_f-1+i_f,2)=&
                            derproj(istart_f-1+i_c,2)+auxproj_f(i_f)
                    end do

                    !derivative wrt z direction
                    nterm=1
                    lx(1)=0 ; ly(1)=1 ; lz(1)=1
                    factor=-sqrt(2.d0)*fpi/(sqrt(gau_a)**5)/gau_a**2
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,derproj(istart_c,3),derproj(istart_f,3))

                 else if (i.eq.2 .and. l.eq.3) then  ! pz type projector

                    !derivative wrt x direction
                    nterm=1
                    lx(1)=1 ; ly(1)=0 ; lz(1)=1
                    factor=-sqrt(2.d0)*fpi/(sqrt(gau_a)**5)/gau_a**2
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,derproj(istart_c,1),derproj(istart_f,1))

                    !derivative wrt y direction
                    nterm=1
                    lx(1)=0 ; ly(1)=1 ; lz(1)=1
                    factor=-sqrt(2.d0)*fpi/(sqrt(gau_a)**5)/gau_a**2
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,derproj(istart_c,2),derproj(istart_f,2))


                    !derivative wrt z direction (first part)
                    nterm=1 
                    lx(1)=0 ; ly(1)=0 ; lz(1)=2 
                    factor=-sqrt(2.d0)*fpi/(sqrt(gau_a)**5)/gau_a**2
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,derproj(istart_c,3),derproj(istart_f,3))
                    !derivative wrt z direction (second part)
                    nterm=1 
                    lx(1)=0 ; ly(1)=0 ; lz(1)=0 
                    factor=sqrt(2.d0)*fpi/(sqrt(gau_a)**5)
                    call crtproj(nterm,n1,n2,n3, & 
                         nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
                         radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                         mvctr_c,mvctr_f,auxproj_c,auxproj_f)
                    !assemble the two parts
                    do i_c=1,mvctr_c
                       derproj(istart_c-1+i_c,3)=&
                            derproj(istart_c-1+i_c,3)+auxproj_c(i_c)
                    end do
                    do i_f=1,7*mvctr_f
                       derproj(istart_f-1+i_f,2)=&
                            derproj(istart_f-1+i_c,3)+auxproj_f(i_f)
                    end do

                 else
                    stop 'PSP format error'
                 end if

                 iproj=iproj+1

                 istart_c=istart_f+7*mvctr_f
                 if (istart_c.gt.istart) stop 'istart_c > istart'

              end do
           end if
        end do
     end do
     deallocate(auxproj_c,auxproj_f)
  end do
  if (iproj.ne.nproj) stop 'incorrect number of projectors created'
  ! projector part finished

  allocate(fxyz_orb(3,nat))
  print *,'end of the projector part'


!  fsep(:,:)=0.d0

  onem=1.d0-eps_mach
  ! loop over all my orbitals
  do iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)

     ! loop over all projectors
     fxyz_orb(:,:)=0.d0
     iproj=0
     eproj=0.d0
     istart_c=1
     do iat=1,nat
     fx=0.d0
     fy=0.d0
     fz=0.d0
        mbseg_c=nseg_p(2*iat-1)-nseg_p(2*iat-2)
        mbseg_f=nseg_p(2*iat  )-nseg_p(2*iat-1)
        jseg_c=nseg_p(2*iat-2)+1
        jseg_f=nseg_p(2*iat-1)+1
        mbvctr_c=nvctr_p(2*iat-1)-nvctr_p(2*iat-2)
        mbvctr_f=nvctr_p(2*iat  )-nvctr_p(2*iat-1)
        ityp=iatype(iat)
        ! ONLY GTH PSP form (not HGH)
        do i=1,2
           do j=1,2
              if (psppar(i,j,ityp).ne.0.d0) then
                 do l=1,2*i-1
                    iproj=iproj+1
                    istart_f=istart_c+mbvctr_c
                    call wdot(  &
                         nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                         keyg(1,1),keyg(1,nseg_c+1),psi(1,iorb-iproc*norb_p),psi(nvctr_c+1,iorb-iproc*norb_p),  &
                         mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                         keyg_p(1,jseg_c),keyg_p(1,jseg_f),proj(istart_c),proj(istart_f),scpr)

                    ! case with derivative in the x direction
                    call wdot(  &
                         nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                         keyg(1,1),keyg(1,nseg_c+1),psi(1,iorb-iproc*norb_p),psi(nvctr_c+1,iorb-iproc*norb_p), &
                         mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                         keyg_p(1,jseg_c),keyg_p(1,jseg_f),derproj(istart_c,1),derproj(istart_f,1),tcprx)

                    ! case with derivative in the y direction
                    call wdot(  &
                         nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                         keyg(1,1),keyg(1,nseg_c+1),psi(1,iorb-iproc*norb_p),psi(nvctr_c+1,iorb-iproc*norb_p), &
                         mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                         keyg_p(1,jseg_c),keyg_p(1,jseg_f),derproj(istart_c,2),derproj(istart_f,2),tcpry)

                    ! case with derivative in the z direction
                    call wdot(  &
                         nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                         keyg(1,1),keyg(1,nseg_c+1),psi(1,iorb-iproc*norb_p),psi(nvctr_c+1,iorb-iproc*norb_p), &
                         mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                         keyg_p(1,jseg_c),keyg_p(1,jseg_f),derproj(istart_c,3),derproj(istart_f,3),tcprz)


                    scprp=scpr*psppar(i,j,ityp)

                    fxyz_orb(1,iat)=fxyz_orb(1,iat)+scprp*tcprx
                    fxyz_orb(2,iat)=fxyz_orb(2,iat)+scprp*tcpry
                    fxyz_orb(3,iat)=fxyz_orb(3,iat)+scprp*tcprz

                    istart_c=istart_f+7*mbvctr_f
                 end do
              end if
           end do
        end do
     end do

     do iat=1,nat
        fsep(1,iat)=fsep(1,iat)+2*occup(iorb)*fxyz_orb(1,iat)
        fsep(2,iat)=fsep(2,iat)+2*occup(iorb)*fxyz_orb(2,iat)
        fsep(3,iat)=fsep(3,iat)+2*occup(iorb)*fxyz_orb(3,iat)
     end do

     if (iproj.ne.nproj) stop '1:applyprojectors'
     if (istart_c-1.ne.nprojel) stop '2:applyprojectors'
  end do

  deallocate(fxyz_orb,derproj)

end subroutine nonlocal_forces

