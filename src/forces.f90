subroutine local_forces(iproc,nproc,ntypes,nat,iatype,atomnames,rxyz,psppar,nelpsp,hgrid,&
     n1,n2,n3,rho,pot,floc)
! Calculates the local forces acting on the atoms belonging to iproc
  use libBigDFT
  
  implicit none
  !Arguments---------
  integer, intent(in) :: iproc,nproc,ntypes,nat,n1,n2,n3
  real(kind=8), intent(in) :: hgrid
  character*20, dimension(100), intent(in) :: atomnames
  real(kind=8), dimension(0:4,0:4,ntypes), intent(in) :: psppar
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(in) :: rho,pot
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(kind=8), dimension(3,nat), intent(out) :: floc
  !Local variables---------
  real(kind=8) :: hgridh,pi,prefactor,cutoff,rloc,Vel,rhoel
  real(kind=8) :: fxerf,fyerf,fzerf,fxion,fyion,fzion,fxgau,fygau,fzgau,forceleaked,forceloc
  real(kind=8) :: rx,ry,rz,x,y,z,arg,r2,xp,dist,tt
  integer :: ii,ix,iy,iz,i1,i2,i3,iat,jat,ityp,jtyp,nloc,iloc,i_all,i_stat
  !array of coefficients of the derivative
  real(kind=8), dimension(4) :: cprime 

  hgridh=hgrid*.5d0 
  pi=4.d0*atan(1.d0)

  do iat=1,nat
  if (mod(iat-1,nproc).eq.iproc) then
     write(*,'(1x,i0,a,i0)') iproc,' calculates local force on atom ',iat
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
                 forceleaked=forceleaked+xp*(1.d0+tt)
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
  if (iproc.eq.0) write(*,'(a,e21.14,1x,e10.3)') ' leaked force: ',forceleaked

end subroutine local_forces

subroutine nonlocal_forces(iproc,nproc,n1,n2,n3,nboxp_c,nboxp_f, & 
     ntypes,nat,norb,norb_p,istart,nprojel,nproj,&
     iatype,psppar,npspcode,occup,nseg_c,nseg_f,nvctr_c,nvctr_f,nseg_p,nvctr_p,proj,  &
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
  real(kind=8), dimension(0:4,0:4,ntypes), intent(in) :: psppar
  integer, dimension(ntypes), intent(in) :: npspcode
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(ntypes,2), intent(in) :: radii_cf
  real(kind=8), dimension(norb), intent(in) :: occup
  real(kind=8), dimension(nprojel), intent(in) :: proj
  real(kind=8), dimension(nvctr_c+7*nvctr_f,norb_p), intent(in) :: psi
  real(kind=8), dimension(3,nat), intent(out) :: fsep
  !Local Variables--------------
  real(kind=8), parameter :: eps_mach=1.d-12
  real(kind=8), dimension(:,:), allocatable :: derproj,fxyz_orb,fac_arr
  real(kind=8), dimension(:), allocatable :: auxproj_c,auxproj_f
  real(kind=8), dimension(:,:,:,:), allocatable :: scalprod
  integer, dimension(:,:,:), allocatable :: lxyz_arr
  integer, dimension(:), allocatable :: nterm_arr,lx,ly,lz
  integer :: istart_c,istart_f,iproj,iat,ityp,i,j,l,m,nterm
  integer :: istart_c_i,istart_f_i,istart_c_j,istart_f_j
  integer :: mvctr_c,mvctr_f,mbseg_c,mbseg_f,maseg_c,maseg_f,iseg_c,iseg_f,jseg_c,jseg_f
  integer :: nl1_c,nl2_c,nl3_c,nl1_f,nl2_f,nl3_f,nu1_c,nu2_c,nu3_c,nu1_f,nu2_f,nu3_f
  integer :: mavctr_c,mavctr_f,mbvctr_c,mbvctr_f,ipsi_c,ipsi_f,iorb,i_c,i_f
  real(kind=8) :: fpi,factor,gau_a,onem,scpr,scprp,tcprx,tcpry,tcprz,rx,ry,rz,eproj,fx,fy,fz
  real(kind=8) :: scpr_i,scpr_j,scprp_i,scprp_j,tcprx_i,tcprx_j,tcpry_i,tcpry_j,tcprz_i,tcprz_j
  real(kind=8) :: offdiagcoeff,hij
  integer :: idir,iadd,iterm,nterm_max,i_all,i_stat
  nterm_max=20 !if GTH nterm_max=4
  allocate(derproj(nprojel,3),stat=i_all)
  allocate(fac_arr(nterm_max,3),stat=i_stat)
  i_all=i_all+i_stat
  allocate(lxyz_arr(3,nterm_max,3),stat=i_stat)
  i_all=i_all+i_stat
  allocate(lx(nterm_max),stat=i_stat)
  i_all=i_all+i_stat
  allocate(ly(nterm_max),stat=i_stat)
  i_all=i_all+i_stat
  allocate(lz(nterm_max),stat=i_stat)
  i_all=i_all+i_stat
  allocate(nterm_arr(3),stat=i_stat)
  i_all=i_all+i_stat
  allocate(scalprod(0:3,4,3,7),stat=i_stat)
  if (i_all+i_stat /= 0) then
     write(*,*)' nonlocal_forces: problem of memory allocation'
     stop
  end if

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


     do l=1,4!compatibility between GTH and HGH form
        do i=1,3
           if (psppar(l,i,ityp).ne.0.d0) then
              gau_a=psppar(l,0,ityp)
              factor=sqrt(2.d0)*fpi/(sqrt(gau_a)**(2*(l-1)+4*i-1))
              do m=1,2*l-1

                 istart_f=istart_c+mvctr_c

                 call calc_coeff_derproj(l,i,m,nterm_max,gau_a,nterm_arr,lxyz_arr,fac_arr)

                 do idir=1,3
                    nterm=nterm_arr(idir)
                    fac_arr(1:nterm,idir)=factor*fac_arr(1:nterm,idir)
                    do iterm=1,nterm
                       lx(iterm)=lxyz_arr(1,iterm,idir)
                       ly(iterm)=lxyz_arr(2,iterm,idir)
                       lz(iterm)=lxyz_arr(3,iterm,idir)
                    end do

                    call crtproj(iproc,nterm,n1,n2,n3,nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,&
                         nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,radii_cf(iatype(iat),2),&
                         cpmult,fpmult,hgrid,gau_a,fac_arr(1,idir),rx,ry,rz,lx,ly,lz,&
                         mvctr_c,mvctr_f,derproj(istart_c,idir),derproj(istart_f,idir))

                 end do
                 
                 iproj=iproj+1

                 istart_c=istart_f+7*mvctr_f
                 if (istart_c.gt.istart) stop 'istart_c > istart'
                 
              end do
           end if
        end do
     end do

  end do
  if (iproj.ne.nproj) stop 'incorrect number of projectors created'
  ! projector part finished

  allocate(fxyz_orb(3,nat),stat=i_all)
  if (i_all /= 0) then
     write(*,*)' nonlocal_forces: problem of memory allocation'
     stop
  end if
!  print *,'end of the projector part'


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
     scalprod(:,:,:,:)=0.d0
     do l=1,4
        do i=1,3
           if (psppar(l,i,ityp).ne.0.d0) then
              do m=1,2*l-1
                 iproj=iproj+1
                 istart_f=istart_c+mbvctr_c
                 call wpdot(  &
                      nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                      keyg(1,1),keyg(1,nseg_c+1),&
                      psi(1,iorb-iproc*norb_p),psi(nvctr_c+1,iorb-iproc*norb_p),  &
                      mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                      keyg_p(1,jseg_c),keyg_p(1,jseg_f),&
                      proj(istart_c),proj(istart_f),scpr)

                 ! case with derivative in the x direction
                 call wpdot(  &
                      nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                      keyg(1,1),keyg(1,nseg_c+1),&
                      psi(1,iorb-iproc*norb_p),psi(nvctr_c+1,iorb-iproc*norb_p), &
                      mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                      keyg_p(1,jseg_c),keyg_p(1,jseg_f),&
                      derproj(istart_c,1),derproj(istart_f,1),tcprx)

                 ! case with derivative in the y direction
                 call wpdot(  &
                      nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                      keyg(1,1),keyg(1,nseg_c+1),&
                      psi(1,iorb-iproc*norb_p),psi(nvctr_c+1,iorb-iproc*norb_p), &
                      mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                      keyg_p(1,jseg_c),keyg_p(1,jseg_f),&
                      derproj(istart_c,2),derproj(istart_f,2),tcpry)

                 ! case with derivative in the z direction
                 call wpdot(  &
                      nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                      keyg(1,1),keyg(1,nseg_c+1),&
                      psi(1,iorb-iproc*norb_p),psi(nvctr_c+1,iorb-iproc*norb_p), &
                      mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                      keyg_p(1,jseg_c),keyg_p(1,jseg_f),&
                      derproj(istart_c,3),derproj(istart_f,3),tcprz)


                 scalprod(0,l,i,m)=scpr
                 scalprod(1,l,i,m)=tcprx
                 scalprod(2,l,i,m)=tcpry
                 scalprod(3,l,i,m)=tcprz

                 scprp=scpr*psppar(l,i,ityp)

                 fxyz_orb(1,iat)=fxyz_orb(1,iat)+scprp*tcprx
                 fxyz_orb(2,iat)=fxyz_orb(2,iat)+scprp*tcpry
                 fxyz_orb(3,iat)=fxyz_orb(3,iat)+scprp*tcprz

                 istart_c=istart_f+7*mbvctr_f
              end do
           end if
        end do
     end do
     !HGH case, offdiagonal terms
     if (npspcode(ityp) == 3) then
        do l=1,3
           do i=1,2
              if (psppar(l,i,ityp).ne.0.d0) then 
                 loop_j: do j=i+1,3
                    if (psppar(l,j,ityp) .eq. 0.d0) exit loop_j
                    !calculate the coefficients for the off-diagonal terms
                    if (l==1) then
                       if (i==1) then
                          if (j==2) offdiagcoeff=-0.5d0*sqrt(3.d0/5.d0)
                          if (j==3) offdiagcoeff=0.5d0*sqrt(5.d0/21.d0)
                       else
                          offdiagcoeff=-0.5d0*sqrt(100.d0/63.d0)
                       end if
                    else if (l==2) then
                       if (i==1) then
                          if (j==2) offdiagcoeff=-0.5d0*sqrt(5.d0/7.d0)
                          if (j==3) offdiagcoeff=1.d0/6.d0*sqrt(35.d0/11.d0)
                       else
                          offdiagcoeff=-7.d0/3.d0*sqrt(1.d0/11.d0)
                       end if
                    else if (l==3) then
                       if (i==1) then
                          if (j==2) offdiagcoeff=-0.5d0*sqrt(7.d0/9.d0)
                          if (j==3) offdiagcoeff=0.5d0*sqrt(63.d0/143.d0)
                       else
                          offdiagcoeff=-9.d0*sqrt(1.d0/143.d0)
                       end if
                    end if
                    hij=offdiagcoeff*psppar(l,j,ityp)
                    do m=1,2*l-1

                       !F_t= 2.d0*h_ij (<D_tp_i|psi><psi|p_j>+<p_i|psi><psi|D_tp_j>)
                       fxyz_orb(1,iat)=fxyz_orb(1,iat)+&
                            hij*(scalprod(0,l,i,m)*scalprod(1,l,j,m)+&
                            scalprod(1,l,i,m)*scalprod(0,l,j,m))
                       fxyz_orb(2,iat)=fxyz_orb(2,iat)+&
                            hij*(scalprod(0,l,i,m)*scalprod(2,l,j,m)+&
                            scalprod(2,l,i,m)*scalprod(0,l,j,m))
                       fxyz_orb(3,iat)=fxyz_orb(3,iat)+&
                            hij*(scalprod(0,l,i,m)*scalprod(3,l,j,m)+&
                            scalprod(3,l,i,m)*scalprod(0,l,j,m))
                    end do
                 end do loop_j
              end if
           end do
        end do
     end if
  end do

  do iat=1,nat
     fsep(1,iat)=fsep(1,iat)+occup(iorb)*2.d0*fxyz_orb(1,iat)
     fsep(2,iat)=fsep(2,iat)+occup(iorb)*2.d0*fxyz_orb(2,iat)
     fsep(3,iat)=fsep(3,iat)+occup(iorb)*2.d0*fxyz_orb(3,iat)
  end do

  if (iproj.ne.nproj) stop '1:applyprojectors'
  if (istart_c-1.ne.nprojel) stop '2:applyprojectors'
end do

  deallocate(fxyz_orb,stat=i_all)
  deallocate(derproj,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(lxyz_arr,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(nterm_arr,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(fac_arr,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(lx,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(ly,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(lz,stat=i_stat)
  i_all=i_all+i_stat
  deallocate(scalprod,stat=i_stat)
  if (i_all+i_stat /= 0) then
     write(*,*)' nonlocal_forces: problem of memory deallocation'
     stop
  end if
end subroutine nonlocal_forces

subroutine calc_coeff_derproj(l,i,m,nterm_max,rhol,nterm_arr,lxyz_arr,fac_arr)
  implicit none
  integer, intent(in) :: l,i,m,nterm_max
  integer, dimension(3), intent(out) :: nterm_arr
  real(kind=8), intent(in) :: rhol
  integer, dimension(3,nterm_max,3), intent(out) :: lxyz_arr
  real(kind=8), dimension(nterm_max,3), intent(out) :: fac_arr

if (l.eq.1 .and. i.eq.1 .and. m.eq.1) then
   nterm_arr(1)=1
   nterm_arr(2)=1
   nterm_arr(3)=1
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   fac_arr(1,1)=-0.7071067811865475244008444/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   fac_arr(1,2)=-0.7071067811865475244008444/rhol**2d0
   lxyz_arr(1,1,3)=0 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   fac_arr(1,3)=-0.7071067811865475244008444/rhol**2d0
else if (l.eq.1 .and. i.eq.2 .and. m.eq.1) then
   nterm_arr(1)=4
   nterm_arr(2)=4
   nterm_arr(3)=4
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=2 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=2
   fac_arr(1,1)=0.730296743340221484609293d0
   fac_arr(2,1)=-0.3651483716701107423046465/rhol**2d0
   fac_arr(3,1)=-0.3651483716701107423046465/rhol**2d0
   fac_arr(4,1)=-0.3651483716701107423046465/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=2
   fac_arr(1,2)=0.730296743340221484609293d0
   fac_arr(2,2)=-0.3651483716701107423046465/rhol**2d0
   fac_arr(3,2)=-0.3651483716701107423046465/rhol**2d0
   fac_arr(4,2)=-0.3651483716701107423046465/rhol**2d0
   lxyz_arr(1,1,3)=0 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=0 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=2 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=3
   fac_arr(1,3)=0.730296743340221484609293d0
   fac_arr(2,3)=-0.3651483716701107423046465/rhol**2d0
   fac_arr(3,3)=-0.3651483716701107423046465/rhol**2d0
   fac_arr(4,3)=-0.3651483716701107423046465/rhol**2d0
else if (l.eq.1 .and. i.eq.3 .and. m.eq.1) then
   nterm_arr(1)=9
   nterm_arr(2)=9
   nterm_arr(3)=9
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=2
   lxyz_arr(1,4,1)=5 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=4 ; lxyz_arr(3,6,1)=0
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=0 ; lxyz_arr(3,7,1)=2
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=2 ; lxyz_arr(3,8,1)=2
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=0 ; lxyz_arr(3,9,1)=4
   fac_arr(1,1)=0.3680349649825889161579343d0
   fac_arr(2,1)=0.3680349649825889161579343d0
   fac_arr(3,1)=0.3680349649825889161579343d0
   fac_arr(4,1)=-0.09200874124564722903948358/rhol**2d0
   fac_arr(5,1)=-0.1840174824912944580789672/rhol**2d0
   fac_arr(6,1)=-0.09200874124564722903948358/rhol**2d0
   fac_arr(7,1)=-0.1840174824912944580789672/rhol**2d0
   fac_arr(8,1)=-0.1840174824912944580789672/rhol**2d0
   fac_arr(9,1)=-0.09200874124564722903948358/rhol**2d0
   lxyz_arr(1,1,2)=2 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=2
   lxyz_arr(1,4,2)=4 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=2 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=5 ; lxyz_arr(3,6,2)=0
   lxyz_arr(1,7,2)=2 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=2
   lxyz_arr(1,8,2)=0 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=2
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=1 ; lxyz_arr(3,9,2)=4
   fac_arr(1,2)=0.3680349649825889161579343d0
   fac_arr(2,2)=0.3680349649825889161579343d0
   fac_arr(3,2)=0.3680349649825889161579343d0
   fac_arr(4,2)=-0.09200874124564722903948358/rhol**2d0
   fac_arr(5,2)=-0.1840174824912944580789672/rhol**2d0
   fac_arr(6,2)=-0.09200874124564722903948358/rhol**2d0
   fac_arr(7,2)=-0.1840174824912944580789672/rhol**2d0
   fac_arr(8,2)=-0.1840174824912944580789672/rhol**2d0
   fac_arr(9,2)=-0.09200874124564722903948358/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=3
   lxyz_arr(1,4,3)=4 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=2 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=4 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=2 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=3
   lxyz_arr(1,8,3)=0 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=3
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=0 ; lxyz_arr(3,9,3)=5
   fac_arr(1,3)=0.3680349649825889161579343d0
   fac_arr(2,3)=0.3680349649825889161579343d0
   fac_arr(3,3)=0.3680349649825889161579343d0
   fac_arr(4,3)=-0.09200874124564722903948358/rhol**2d0
   fac_arr(5,3)=-0.1840174824912944580789672/rhol**2d0
   fac_arr(6,3)=-0.09200874124564722903948358/rhol**2d0
   fac_arr(7,3)=-0.1840174824912944580789672/rhol**2d0
   fac_arr(8,3)=-0.1840174824912944580789672/rhol**2d0
   fac_arr(9,3)=-0.09200874124564722903948358/rhol**2d0
else if (l.eq.2 .and. i.eq.1 .and. m.eq.1) then
   nterm_arr(1)=2
   nterm_arr(2)=1
   nterm_arr(3)=1
   lxyz_arr(1,1,1)=0 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=0
   fac_arr(1,1)=1.d0
   fac_arr(2,1)=-1./rhol**2d0
   lxyz_arr(1,1,2)=1 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   fac_arr(1,2)=-1./rhol**2d0
   lxyz_arr(1,1,3)=1 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   fac_arr(1,3)=-1./rhol**2d0
else if (l.eq.2 .and. i.eq.1 .and. m.eq.2) then
   nterm_arr(1)=1
   nterm_arr(2)=2
   nterm_arr(3)=1
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   fac_arr(1,1)=-1./rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   fac_arr(1,2)=1.d0
   fac_arr(2,2)=-1./rhol**2d0
   lxyz_arr(1,1,3)=0 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   fac_arr(1,3)=-1./rhol**2d0
else if (l.eq.2 .and. i.eq.1 .and. m.eq.3) then
   nterm_arr(1)=1
   nterm_arr(2)=1
   nterm_arr(3)=2
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   fac_arr(1,1)=-1./rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   fac_arr(1,2)=-1./rhol**2d0
   lxyz_arr(1,1,3)=0 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=0 ; lxyz_arr(3,2,3)=2
   fac_arr(1,3)=1.d0
   fac_arr(2,3)=-1./rhol**2d0
else if (l.eq.2 .and. i.eq.2 .and. m.eq.1) then
   nterm_arr(1)=6
   nterm_arr(2)=4
   nterm_arr(3)=4
   lxyz_arr(1,1,1)=2 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=0 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=2
   lxyz_arr(1,4,1)=4 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=2 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=2 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=2
   fac_arr(1,1)=1.014185105674219893011542d0
   fac_arr(2,1)=0.3380617018914066310038473d0
   fac_arr(3,1)=0.3380617018914066310038473d0
   fac_arr(4,1)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(5,1)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(6,1)=-0.3380617018914066310038473/rhol**2d0
   lxyz_arr(1,1,2)=1 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=3 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=1 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=2
   fac_arr(1,2)=0.6761234037828132620076947d0
   fac_arr(2,2)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(3,2)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(4,2)=-0.3380617018914066310038473/rhol**2d0
   lxyz_arr(1,1,3)=1 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=3 ; lxyz_arr(2,2,3)=0 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=2 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=1 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=3
   fac_arr(1,3)=0.6761234037828132620076947d0
   fac_arr(2,3)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(3,3)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(4,3)=-0.3380617018914066310038473/rhol**2d0
else if (l.eq.2 .and. i.eq.2 .and. m.eq.2) then
   nterm_arr(1)=4
   nterm_arr(2)=6
   nterm_arr(3)=4
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=1 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=3 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=2
   fac_arr(1,1)=0.6761234037828132620076947d0
   fac_arr(2,1)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(3,1)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(4,1)=-0.3380617018914066310038473/rhol**2d0
   lxyz_arr(1,1,2)=2 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=0 ; lxyz_arr(3,3,2)=2
   lxyz_arr(1,4,2)=2 ; lxyz_arr(2,4,2)=2 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=4 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=2 ; lxyz_arr(3,6,2)=2
   fac_arr(1,2)=0.3380617018914066310038473d0
   fac_arr(2,2)=1.014185105674219893011542d0
   fac_arr(3,2)=0.3380617018914066310038473d0
   fac_arr(4,2)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(5,2)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(6,2)=-0.3380617018914066310038473/rhol**2d0
   lxyz_arr(1,1,3)=0 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=1 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=3 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=3
   fac_arr(1,3)=0.6761234037828132620076947d0
   fac_arr(2,3)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(3,3)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(4,3)=-0.3380617018914066310038473/rhol**2d0
else if (l.eq.2 .and. i.eq.2 .and. m.eq.3) then
   nterm_arr(1)=4
   nterm_arr(2)=4
   nterm_arr(3)=6
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=2 ; lxyz_arr(3,3,1)=1
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=3
   fac_arr(1,1)=0.6761234037828132620076947d0
   fac_arr(2,1)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(3,1)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(4,1)=-0.3380617018914066310038473/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=1
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=3
   fac_arr(1,2)=0.6761234037828132620076947d0
   fac_arr(2,2)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(3,2)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(4,2)=-0.3380617018914066310038473/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=2
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=4
   fac_arr(1,3)=0.3380617018914066310038473d0
   fac_arr(2,3)=0.3380617018914066310038473d0
   fac_arr(3,3)=1.014185105674219893011542d0
   fac_arr(4,3)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(5,3)=-0.3380617018914066310038473/rhol**2d0
   fac_arr(6,3)=-0.3380617018914066310038473/rhol**2d0
else if (l.eq.2 .and. i.eq.3 .and. m.eq.1) then
   nterm_arr(1)=12
   nterm_arr(2)=9
   nterm_arr(3)=9
   lxyz_arr(1,1,1)=4 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=2 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=2
   lxyz_arr(1,5,1)=0 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=0 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=4
   lxyz_arr(1,7,1)=6 ; lxyz_arr(2,7,1)=0 ; lxyz_arr(3,7,1)=0
   lxyz_arr(1,8,1)=4 ; lxyz_arr(2,8,1)=2 ; lxyz_arr(3,8,1)=0
   lxyz_arr(1,9,1)=2 ; lxyz_arr(2,9,1)=4 ; lxyz_arr(3,9,1)=0
   lxyz_arr(1,10,1)=4 ; lxyz_arr(2,10,1)=0 ; lxyz_arr(3,10,1)=2
   lxyz_arr(1,11,1)=2 ; lxyz_arr(2,11,1)=2 ; lxyz_arr(3,11,1)=2
   lxyz_arr(1,12,1)=2 ; lxyz_arr(2,12,1)=0 ; lxyz_arr(3,12,1)=4
   fac_arr(1,1)=0.3397647942917503630913594d0
   fac_arr(2,1)=0.4077177531501004357096312d0
   fac_arr(3,1)=0.06795295885835007261827187d0
   fac_arr(4,1)=0.4077177531501004357096312d0
   fac_arr(5,1)=0.1359059177167001452365437d0
   fac_arr(6,1)=0.06795295885835007261827187d0
   fac_arr(7,1)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(8,1)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(9,1)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(10,1)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(11,1)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(12,1)=-0.06795295885835007261827187/rhol**2d0
   lxyz_arr(1,1,2)=3 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=1 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=2
   lxyz_arr(1,4,2)=5 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=3 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=5 ; lxyz_arr(3,6,2)=0
   lxyz_arr(1,7,2)=3 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=2
   lxyz_arr(1,8,2)=1 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=2
   lxyz_arr(1,9,2)=1 ; lxyz_arr(2,9,2)=1 ; lxyz_arr(3,9,2)=4
   fac_arr(1,2)=0.2718118354334002904730875d0
   fac_arr(2,2)=0.2718118354334002904730875d0
   fac_arr(3,2)=0.2718118354334002904730875d0
   fac_arr(4,2)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(5,2)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(6,2)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(7,2)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(8,2)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(9,2)=-0.06795295885835007261827187/rhol**2d0
   lxyz_arr(1,1,3)=3 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=3
   lxyz_arr(1,4,3)=5 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=3 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=1 ; lxyz_arr(2,6,3)=4 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=3 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=3
   lxyz_arr(1,8,3)=1 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=3
   lxyz_arr(1,9,3)=1 ; lxyz_arr(2,9,3)=0 ; lxyz_arr(3,9,3)=5
   fac_arr(1,3)=0.2718118354334002904730875d0
   fac_arr(2,3)=0.2718118354334002904730875d0
   fac_arr(3,3)=0.2718118354334002904730875d0
   fac_arr(4,3)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(5,3)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(6,3)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(7,3)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(8,3)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(9,3)=-0.06795295885835007261827187/rhol**2d0
else if (l.eq.2 .and. i.eq.3 .and. m.eq.2) then
   nterm_arr(1)=9
   nterm_arr(2)=12
   nterm_arr(3)=9
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=1 ; lxyz_arr(3,3,1)=2
   lxyz_arr(1,4,1)=5 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=5 ; lxyz_arr(3,6,1)=0
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=1 ; lxyz_arr(3,7,1)=2
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=3 ; lxyz_arr(3,8,1)=2
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=1 ; lxyz_arr(3,9,1)=4
   fac_arr(1,1)=0.2718118354334002904730875d0
   fac_arr(2,1)=0.2718118354334002904730875d0
   fac_arr(3,1)=0.2718118354334002904730875d0
   fac_arr(4,1)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(5,1)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(6,1)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(7,1)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(8,1)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(9,1)=-0.06795295885835007261827187/rhol**2d0
   lxyz_arr(1,1,2)=4 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=4 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=2 ; lxyz_arr(2,4,2)=0 ; lxyz_arr(3,4,2)=2
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=2 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=0 ; lxyz_arr(3,6,2)=4
   lxyz_arr(1,7,2)=4 ; lxyz_arr(2,7,2)=2 ; lxyz_arr(3,7,2)=0
   lxyz_arr(1,8,2)=2 ; lxyz_arr(2,8,2)=4 ; lxyz_arr(3,8,2)=0
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=6 ; lxyz_arr(3,9,2)=0
   lxyz_arr(1,10,2)=2 ; lxyz_arr(2,10,2)=2 ; lxyz_arr(3,10,2)=2
   lxyz_arr(1,11,2)=0 ; lxyz_arr(2,11,2)=4 ; lxyz_arr(3,11,2)=2
   lxyz_arr(1,12,2)=0 ; lxyz_arr(2,12,2)=2 ; lxyz_arr(3,12,2)=4
   fac_arr(1,2)=0.06795295885835007261827187d0
   fac_arr(2,2)=0.4077177531501004357096312d0
   fac_arr(3,2)=0.3397647942917503630913594d0
   fac_arr(4,2)=0.1359059177167001452365437d0
   fac_arr(5,2)=0.4077177531501004357096312d0
   fac_arr(6,2)=0.06795295885835007261827187d0
   fac_arr(7,2)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(8,2)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(9,2)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(10,2)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(11,2)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(12,2)=-0.06795295885835007261827187/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=1 ; lxyz_arr(3,3,3)=3
   lxyz_arr(1,4,3)=4 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=2 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=5 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=2 ; lxyz_arr(2,7,3)=1 ; lxyz_arr(3,7,3)=3
   lxyz_arr(1,8,3)=0 ; lxyz_arr(2,8,3)=3 ; lxyz_arr(3,8,3)=3
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=1 ; lxyz_arr(3,9,3)=5
   fac_arr(1,3)=0.2718118354334002904730875d0
   fac_arr(2,3)=0.2718118354334002904730875d0
   fac_arr(3,3)=0.2718118354334002904730875d0
   fac_arr(4,3)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(5,3)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(6,3)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(7,3)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(8,3)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(9,3)=-0.06795295885835007261827187/rhol**2d0
else if (l.eq.2 .and. i.eq.3 .and. m.eq.3) then
   nterm_arr(1)=9
   nterm_arr(2)=9
   nterm_arr(3)=12
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=3
   lxyz_arr(1,4,1)=5 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=1
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=1
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=4 ; lxyz_arr(3,6,1)=1
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=0 ; lxyz_arr(3,7,1)=3
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=2 ; lxyz_arr(3,8,1)=3
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=0 ; lxyz_arr(3,9,1)=5
   fac_arr(1,1)=0.2718118354334002904730875d0
   fac_arr(2,1)=0.2718118354334002904730875d0
   fac_arr(3,1)=0.2718118354334002904730875d0
   fac_arr(4,1)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(5,1)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(6,1)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(7,1)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(8,1)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(9,1)=-0.06795295885835007261827187/rhol**2d0
   lxyz_arr(1,1,2)=2 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=3
   lxyz_arr(1,4,2)=4 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=1
   lxyz_arr(1,5,2)=2 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=1
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=5 ; lxyz_arr(3,6,2)=1
   lxyz_arr(1,7,2)=2 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=3
   lxyz_arr(1,8,2)=0 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=3
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=1 ; lxyz_arr(3,9,2)=5
   fac_arr(1,2)=0.2718118354334002904730875d0
   fac_arr(2,2)=0.2718118354334002904730875d0
   fac_arr(3,2)=0.2718118354334002904730875d0
   fac_arr(4,2)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(5,2)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(6,2)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(7,2)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(8,2)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(9,2)=-0.06795295885835007261827187/rhol**2d0
   lxyz_arr(1,1,3)=4 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=4 ; lxyz_arr(3,3,3)=0
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=4
   lxyz_arr(1,7,3)=4 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=2
   lxyz_arr(1,8,3)=2 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=2
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=4 ; lxyz_arr(3,9,3)=2
   lxyz_arr(1,10,3)=2 ; lxyz_arr(2,10,3)=0 ; lxyz_arr(3,10,3)=4
   lxyz_arr(1,11,3)=0 ; lxyz_arr(2,11,3)=2 ; lxyz_arr(3,11,3)=4
   lxyz_arr(1,12,3)=0 ; lxyz_arr(2,12,3)=0 ; lxyz_arr(3,12,3)=6
   fac_arr(1,3)=0.06795295885835007261827187d0
   fac_arr(2,3)=0.1359059177167001452365437d0
   fac_arr(3,3)=0.06795295885835007261827187d0
   fac_arr(4,3)=0.4077177531501004357096312d0
   fac_arr(5,3)=0.4077177531501004357096312d0
   fac_arr(6,3)=0.3397647942917503630913594d0
   fac_arr(7,3)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(8,3)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(9,3)=-0.06795295885835007261827187/rhol**2d0
   fac_arr(10,3)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(11,3)=-0.1359059177167001452365437/rhol**2d0
   fac_arr(12,3)=-0.06795295885835007261827187/rhol**2d0
else if (l.eq.3 .and. i.eq.1 .and. m.eq.1) then
   nterm_arr(1)=1
   nterm_arr(2)=2
   nterm_arr(3)=2
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=1
   fac_arr(1,1)=-1.414213562373095048801689/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=1
   fac_arr(1,2)=1.414213562373095048801689d0
   fac_arr(2,2)=-1.414213562373095048801689/rhol**2d0
   lxyz_arr(1,1,3)=0 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=1 ; lxyz_arr(3,2,3)=2
   fac_arr(1,3)=1.414213562373095048801689d0
   fac_arr(2,3)=-1.414213562373095048801689/rhol**2d0
else if (l.eq.3 .and. i.eq.1 .and. m.eq.2) then
   nterm_arr(1)=2
   nterm_arr(2)=1
   nterm_arr(3)=2
   lxyz_arr(1,1,1)=0 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=1
   fac_arr(1,1)=1.414213562373095048801689d0
   fac_arr(2,1)=-1.414213562373095048801689/rhol**2d0
   lxyz_arr(1,1,2)=1 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   fac_arr(1,2)=-1.414213562373095048801689/rhol**2d0
   lxyz_arr(1,1,3)=1 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=0 ; lxyz_arr(3,2,3)=2
   fac_arr(1,3)=1.414213562373095048801689d0
   fac_arr(2,3)=-1.414213562373095048801689/rhol**2d0
else if (l.eq.3 .and. i.eq.1 .and. m.eq.3) then
   nterm_arr(1)=2
   nterm_arr(2)=2
   nterm_arr(3)=1
   lxyz_arr(1,1,1)=0 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=1 ; lxyz_arr(3,2,1)=0
   fac_arr(1,1)=1.414213562373095048801689d0
   fac_arr(2,1)=-1.414213562373095048801689/rhol**2d0
   lxyz_arr(1,1,2)=1 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=1 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   fac_arr(1,2)=1.414213562373095048801689d0
   fac_arr(2,2)=-1.414213562373095048801689/rhol**2d0
   lxyz_arr(1,1,3)=1 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   fac_arr(1,3)=-1.414213562373095048801689/rhol**2d0
else if (l.eq.3 .and. i.eq.1 .and. m.eq.4) then
   nterm_arr(1)=3
   nterm_arr(2)=3
   nterm_arr(3)=2
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=2 ; lxyz_arr(3,3,1)=0
   fac_arr(1,1)=1.414213562373095048801689d0
   fac_arr(2,1)=-0.7071067811865475244008444/rhol**2d0
   fac_arr(3,1)=0.7071067811865475244008444/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=0
   fac_arr(1,2)=-1.414213562373095048801689d0
   fac_arr(2,2)=-0.7071067811865475244008444/rhol**2d0
   fac_arr(3,2)=0.7071067811865475244008444/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   fac_arr(1,3)=-0.7071067811865475244008444/rhol**2d0
   fac_arr(2,3)=0.7071067811865475244008444/rhol**2d0
else if (l.eq.3 .and. i.eq.1 .and. m.eq.5) then
   nterm_arr(1)=4
   nterm_arr(2)=4
   nterm_arr(3)=4
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=2 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=2
   fac_arr(1,1)=-0.816496580927726032732428d0
   fac_arr(2,1)=0.408248290463863016366214/rhol**2d0
   fac_arr(3,1)=0.408248290463863016366214/rhol**2d0
   fac_arr(4,1)=-0.816496580927726032732428/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=2
   fac_arr(1,2)=-0.816496580927726032732428d0
   fac_arr(2,2)=0.408248290463863016366214/rhol**2d0
   fac_arr(3,2)=0.408248290463863016366214/rhol**2d0
   fac_arr(4,2)=-0.816496580927726032732428/rhol**2d0
   lxyz_arr(1,1,3)=0 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=0 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=2 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=3
   fac_arr(1,3)=1.632993161855452065464856d0
   fac_arr(2,3)=0.408248290463863016366214/rhol**2d0
   fac_arr(3,3)=0.408248290463863016366214/rhol**2d0
   fac_arr(4,3)=-0.816496580927726032732428/rhol**2d0
else if (l.eq.3 .and. i.eq.2 .and. m.eq.1) then
   nterm_arr(1)=4
   nterm_arr(2)=6
   nterm_arr(3)=6
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=1 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=3 ; lxyz_arr(3,3,1)=1
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=3
   fac_arr(1,1)=0.7126966450997983591588093d0
   fac_arr(2,1)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(3,1)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(4,1)=-0.3563483225498991795794046/rhol**2d0
   lxyz_arr(1,1,2)=2 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=0 ; lxyz_arr(3,3,2)=3
   lxyz_arr(1,4,2)=2 ; lxyz_arr(2,4,2)=2 ; lxyz_arr(3,4,2)=1
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=4 ; lxyz_arr(3,5,2)=1
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=2 ; lxyz_arr(3,6,2)=3
   fac_arr(1,2)=0.3563483225498991795794046d0
   fac_arr(2,2)=1.069044967649697538738214d0
   fac_arr(3,2)=0.3563483225498991795794046d0
   fac_arr(4,2)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(5,2)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(6,2)=-0.3563483225498991795794046/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=1 ; lxyz_arr(3,3,3)=2
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=1 ; lxyz_arr(3,6,3)=4
   fac_arr(1,3)=0.3563483225498991795794046d0
   fac_arr(2,3)=0.3563483225498991795794046d0
   fac_arr(3,3)=1.069044967649697538738214d0
   fac_arr(4,3)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(5,3)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(6,3)=-0.3563483225498991795794046/rhol**2d0
else if (l.eq.3 .and. i.eq.2 .and. m.eq.2) then
   nterm_arr(1)=6
   nterm_arr(2)=4
   nterm_arr(3)=6
   lxyz_arr(1,1,1)=2 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=0 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=3
   lxyz_arr(1,4,1)=4 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=1
   lxyz_arr(1,5,1)=2 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=1
   lxyz_arr(1,6,1)=2 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=3
   fac_arr(1,1)=1.069044967649697538738214d0
   fac_arr(2,1)=0.3563483225498991795794046d0
   fac_arr(3,1)=0.3563483225498991795794046d0
   fac_arr(4,1)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(5,1)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(6,1)=-0.3563483225498991795794046/rhol**2d0
   lxyz_arr(1,1,2)=1 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=3 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=1
   lxyz_arr(1,4,2)=1 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=3
   fac_arr(1,2)=0.7126966450997983591588093d0
   fac_arr(2,2)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(3,2)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(4,2)=-0.3563483225498991795794046/rhol**2d0
   lxyz_arr(1,1,3)=3 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=2
   lxyz_arr(1,4,3)=3 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=1 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=1 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=4
   fac_arr(1,3)=0.3563483225498991795794046d0
   fac_arr(2,3)=0.3563483225498991795794046d0
   fac_arr(3,3)=1.069044967649697538738214d0
   fac_arr(4,3)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(5,3)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(6,3)=-0.3563483225498991795794046/rhol**2d0
else if (l.eq.3 .and. i.eq.2 .and. m.eq.3) then
   nterm_arr(1)=6
   nterm_arr(2)=6
   nterm_arr(3)=4
   lxyz_arr(1,1,1)=2 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=0 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=1 ; lxyz_arr(3,3,1)=2
   lxyz_arr(1,4,1)=4 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=2 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=2 ; lxyz_arr(2,6,1)=1 ; lxyz_arr(3,6,1)=2
   fac_arr(1,1)=1.069044967649697538738214d0
   fac_arr(2,1)=0.3563483225498991795794046d0
   fac_arr(3,1)=0.3563483225498991795794046d0
   fac_arr(4,1)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(5,1)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(6,1)=-0.3563483225498991795794046/rhol**2d0
   lxyz_arr(1,1,2)=3 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=1 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=0 ; lxyz_arr(3,3,2)=2
   lxyz_arr(1,4,2)=3 ; lxyz_arr(2,4,2)=2 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=1 ; lxyz_arr(2,5,2)=4 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=2 ; lxyz_arr(3,6,2)=2
   fac_arr(1,2)=0.3563483225498991795794046d0
   fac_arr(2,2)=1.069044967649697538738214d0
   fac_arr(3,2)=0.3563483225498991795794046d0
   fac_arr(4,2)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(5,2)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(6,2)=-0.3563483225498991795794046/rhol**2d0
   lxyz_arr(1,1,3)=1 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=3 ; lxyz_arr(2,2,3)=1 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=3 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=1 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=3
   fac_arr(1,3)=0.7126966450997983591588093d0
   fac_arr(2,3)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(3,3)=-0.3563483225498991795794046/rhol**2d0
   fac_arr(4,3)=-0.3563483225498991795794046/rhol**2d0
else if (l.eq.3 .and. i.eq.2 .and. m.eq.4) then
   nterm_arr(1)=6
   nterm_arr(2)=6
   nterm_arr(3)=6
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=2
   lxyz_arr(1,3,1)=5 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=4 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=0 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=2 ; lxyz_arr(3,6,1)=2
   fac_arr(1,1)=0.7126966450997983591588093d0
   fac_arr(2,1)=0.3563483225498991795794046d0
   fac_arr(3,1)=-0.1781741612749495897897023/rhol**2d0
   fac_arr(4,1)=0.1781741612749495897897023/rhol**2d0
   fac_arr(5,1)=-0.1781741612749495897897023/rhol**2d0
   fac_arr(6,1)=0.1781741612749495897897023/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=3 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=2
   lxyz_arr(1,3,2)=4 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=5 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=2 ; lxyz_arr(2,5,2)=1 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=3 ; lxyz_arr(3,6,2)=2
   fac_arr(1,2)=-0.7126966450997983591588093d0
   fac_arr(2,2)=-0.3563483225498991795794046d0
   fac_arr(3,2)=-0.1781741612749495897897023/rhol**2d0
   fac_arr(4,2)=0.1781741612749495897897023/rhol**2d0
   fac_arr(5,2)=-0.1781741612749495897897023/rhol**2d0
   fac_arr(6,2)=0.1781741612749495897897023/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=4 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=4 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=2 ; lxyz_arr(2,5,3)=0 ; lxyz_arr(3,5,3)=3
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=2 ; lxyz_arr(3,6,3)=3
   fac_arr(1,3)=0.3563483225498991795794046d0
   fac_arr(2,3)=-0.3563483225498991795794046d0
   fac_arr(3,3)=-0.1781741612749495897897023/rhol**2d0
   fac_arr(4,3)=0.1781741612749495897897023/rhol**2d0
   fac_arr(5,3)=-0.1781741612749495897897023/rhol**2d0
   fac_arr(6,3)=0.1781741612749495897897023/rhol**2d0
else if (l.eq.3 .and. i.eq.2 .and. m.eq.5) then
   nterm_arr(1)=9
   nterm_arr(2)=9
   nterm_arr(3)=9
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=2
   lxyz_arr(1,4,1)=5 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=4 ; lxyz_arr(3,6,1)=0
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=0 ; lxyz_arr(3,7,1)=2
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=2 ; lxyz_arr(3,8,1)=2
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=0 ; lxyz_arr(3,9,1)=4
   fac_arr(1,1)=-0.4114755998989117606962519d0
   fac_arr(2,1)=-0.4114755998989117606962519d0
   fac_arr(3,1)=0.205737799949455880348126d0
   fac_arr(4,1)=0.102868899974727940174063/rhol**2d0
   fac_arr(5,1)=0.205737799949455880348126/rhol**2d0
   fac_arr(6,1)=0.102868899974727940174063/rhol**2d0
   fac_arr(7,1)=-0.102868899974727940174063/rhol**2d0
   fac_arr(8,1)=-0.102868899974727940174063/rhol**2d0
   fac_arr(9,1)=-0.205737799949455880348126/rhol**2d0
   lxyz_arr(1,1,2)=2 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=2
   lxyz_arr(1,4,2)=4 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=2 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=5 ; lxyz_arr(3,6,2)=0
   lxyz_arr(1,7,2)=2 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=2
   lxyz_arr(1,8,2)=0 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=2
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=1 ; lxyz_arr(3,9,2)=4
   fac_arr(1,2)=-0.4114755998989117606962519d0
   fac_arr(2,2)=-0.4114755998989117606962519d0
   fac_arr(3,2)=0.205737799949455880348126d0
   fac_arr(4,2)=0.102868899974727940174063/rhol**2d0
   fac_arr(5,2)=0.205737799949455880348126/rhol**2d0
   fac_arr(6,2)=0.102868899974727940174063/rhol**2d0
   fac_arr(7,2)=-0.102868899974727940174063/rhol**2d0
   fac_arr(8,2)=-0.102868899974727940174063/rhol**2d0
   fac_arr(9,2)=-0.205737799949455880348126/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=3
   lxyz_arr(1,4,3)=4 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=2 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=4 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=2 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=3
   lxyz_arr(1,8,3)=0 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=3
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=0 ; lxyz_arr(3,9,3)=5
   fac_arr(1,3)=0.205737799949455880348126d0
   fac_arr(2,3)=0.205737799949455880348126d0
   fac_arr(3,3)=0.8229511997978235213925038d0
   fac_arr(4,3)=0.102868899974727940174063/rhol**2d0
   fac_arr(5,3)=0.205737799949455880348126/rhol**2d0
   fac_arr(6,3)=0.102868899974727940174063/rhol**2d0
   fac_arr(7,3)=-0.102868899974727940174063/rhol**2d0
   fac_arr(8,3)=-0.102868899974727940174063/rhol**2d0
   fac_arr(9,3)=-0.205737799949455880348126/rhol**2d0
else if (l.eq.3 .and. i.eq.3 .and. m.eq.1) then
   nterm_arr(1)=9
   nterm_arr(2)=12
   nterm_arr(3)=12
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=1 ; lxyz_arr(3,3,1)=3
   lxyz_arr(1,4,1)=5 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=1
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=1
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=5 ; lxyz_arr(3,6,1)=1
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=1 ; lxyz_arr(3,7,1)=3
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=3 ; lxyz_arr(3,8,1)=3
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=1 ; lxyz_arr(3,9,1)=5
   fac_arr(1,1)=0.2383947500094262395810797d0
   fac_arr(2,1)=0.2383947500094262395810797d0
   fac_arr(3,1)=0.2383947500094262395810797d0
   fac_arr(4,1)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(5,1)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(6,1)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(7,1)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(8,1)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(9,1)=-0.05959868750235655989526993/rhol**2d0
   lxyz_arr(1,1,2)=4 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=4 ; lxyz_arr(3,3,2)=1
   lxyz_arr(1,4,2)=2 ; lxyz_arr(2,4,2)=0 ; lxyz_arr(3,4,2)=3
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=2 ; lxyz_arr(3,5,2)=3
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=0 ; lxyz_arr(3,6,2)=5
   lxyz_arr(1,7,2)=4 ; lxyz_arr(2,7,2)=2 ; lxyz_arr(3,7,2)=1
   lxyz_arr(1,8,2)=2 ; lxyz_arr(2,8,2)=4 ; lxyz_arr(3,8,2)=1
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=6 ; lxyz_arr(3,9,2)=1
   lxyz_arr(1,10,2)=2 ; lxyz_arr(2,10,2)=2 ; lxyz_arr(3,10,2)=3
   lxyz_arr(1,11,2)=0 ; lxyz_arr(2,11,2)=4 ; lxyz_arr(3,11,2)=3
   lxyz_arr(1,12,2)=0 ; lxyz_arr(2,12,2)=2 ; lxyz_arr(3,12,2)=5
   fac_arr(1,2)=0.05959868750235655989526993d0
   fac_arr(2,2)=0.3575921250141393593716196d0
   fac_arr(3,2)=0.2979934375117827994763496d0
   fac_arr(4,2)=0.1191973750047131197905399d0
   fac_arr(5,2)=0.3575921250141393593716196d0
   fac_arr(6,2)=0.05959868750235655989526993d0
   fac_arr(7,2)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(8,2)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(9,2)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(10,2)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(11,2)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(12,2)=-0.05959868750235655989526993/rhol**2d0
   lxyz_arr(1,1,3)=4 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=5 ; lxyz_arr(3,3,3)=0
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=1 ; lxyz_arr(3,6,3)=4
   lxyz_arr(1,7,3)=4 ; lxyz_arr(2,7,3)=1 ; lxyz_arr(3,7,3)=2
   lxyz_arr(1,8,3)=2 ; lxyz_arr(2,8,3)=3 ; lxyz_arr(3,8,3)=2
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=5 ; lxyz_arr(3,9,3)=2
   lxyz_arr(1,10,3)=2 ; lxyz_arr(2,10,3)=1 ; lxyz_arr(3,10,3)=4
   lxyz_arr(1,11,3)=0 ; lxyz_arr(2,11,3)=3 ; lxyz_arr(3,11,3)=4
   lxyz_arr(1,12,3)=0 ; lxyz_arr(2,12,3)=1 ; lxyz_arr(3,12,3)=6
   fac_arr(1,3)=0.05959868750235655989526993d0
   fac_arr(2,3)=0.1191973750047131197905399d0
   fac_arr(3,3)=0.05959868750235655989526993d0
   fac_arr(4,3)=0.3575921250141393593716196d0
   fac_arr(5,3)=0.3575921250141393593716196d0
   fac_arr(6,3)=0.2979934375117827994763496d0
   fac_arr(7,3)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(8,3)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(9,3)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(10,3)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(11,3)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(12,3)=-0.05959868750235655989526993/rhol**2d0
else if (l.eq.3 .and. i.eq.3 .and. m.eq.2) then
   nterm_arr(1)=12
   nterm_arr(2)=9
   nterm_arr(3)=12
   lxyz_arr(1,1,1)=4 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=1
   lxyz_arr(1,4,1)=2 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=3
   lxyz_arr(1,5,1)=0 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=3
   lxyz_arr(1,6,1)=0 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=5
   lxyz_arr(1,7,1)=6 ; lxyz_arr(2,7,1)=0 ; lxyz_arr(3,7,1)=1
   lxyz_arr(1,8,1)=4 ; lxyz_arr(2,8,1)=2 ; lxyz_arr(3,8,1)=1
   lxyz_arr(1,9,1)=2 ; lxyz_arr(2,9,1)=4 ; lxyz_arr(3,9,1)=1
   lxyz_arr(1,10,1)=4 ; lxyz_arr(2,10,1)=0 ; lxyz_arr(3,10,1)=3
   lxyz_arr(1,11,1)=2 ; lxyz_arr(2,11,1)=2 ; lxyz_arr(3,11,1)=3
   lxyz_arr(1,12,1)=2 ; lxyz_arr(2,12,1)=0 ; lxyz_arr(3,12,1)=5
   fac_arr(1,1)=0.2979934375117827994763496d0
   fac_arr(2,1)=0.3575921250141393593716196d0
   fac_arr(3,1)=0.05959868750235655989526993d0
   fac_arr(4,1)=0.3575921250141393593716196d0
   fac_arr(5,1)=0.1191973750047131197905399d0
   fac_arr(6,1)=0.05959868750235655989526993d0
   fac_arr(7,1)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(8,1)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(9,1)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(10,1)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(11,1)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(12,1)=-0.05959868750235655989526993/rhol**2d0
   lxyz_arr(1,1,2)=3 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=1 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=3
   lxyz_arr(1,4,2)=5 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=1
   lxyz_arr(1,5,2)=3 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=1
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=5 ; lxyz_arr(3,6,2)=1
   lxyz_arr(1,7,2)=3 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=3
   lxyz_arr(1,8,2)=1 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=3
   lxyz_arr(1,9,2)=1 ; lxyz_arr(2,9,2)=1 ; lxyz_arr(3,9,2)=5
   fac_arr(1,2)=0.2383947500094262395810797d0
   fac_arr(2,2)=0.2383947500094262395810797d0
   fac_arr(3,2)=0.2383947500094262395810797d0
   fac_arr(4,2)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(5,2)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(6,2)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(7,2)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(8,2)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(9,2)=-0.05959868750235655989526993/rhol**2d0
   lxyz_arr(1,1,3)=5 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=3 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=4 ; lxyz_arr(3,3,3)=0
   lxyz_arr(1,4,3)=3 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=1 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=1 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=4
   lxyz_arr(1,7,3)=5 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=2
   lxyz_arr(1,8,3)=3 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=2
   lxyz_arr(1,9,3)=1 ; lxyz_arr(2,9,3)=4 ; lxyz_arr(3,9,3)=2
   lxyz_arr(1,10,3)=3 ; lxyz_arr(2,10,3)=0 ; lxyz_arr(3,10,3)=4
   lxyz_arr(1,11,3)=1 ; lxyz_arr(2,11,3)=2 ; lxyz_arr(3,11,3)=4
   lxyz_arr(1,12,3)=1 ; lxyz_arr(2,12,3)=0 ; lxyz_arr(3,12,3)=6
   fac_arr(1,3)=0.05959868750235655989526993d0
   fac_arr(2,3)=0.1191973750047131197905399d0
   fac_arr(3,3)=0.05959868750235655989526993d0
   fac_arr(4,3)=0.3575921250141393593716196d0
   fac_arr(5,3)=0.3575921250141393593716196d0
   fac_arr(6,3)=0.2979934375117827994763496d0
   fac_arr(7,3)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(8,3)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(9,3)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(10,3)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(11,3)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(12,3)=-0.05959868750235655989526993/rhol**2d0
else if (l.eq.3 .and. i.eq.3 .and. m.eq.3) then
   nterm_arr(1)=12
   nterm_arr(2)=12
   nterm_arr(3)=9
   lxyz_arr(1,1,1)=4 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=5 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=2 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=2
   lxyz_arr(1,5,1)=0 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=0 ; lxyz_arr(2,6,1)=1 ; lxyz_arr(3,6,1)=4
   lxyz_arr(1,7,1)=6 ; lxyz_arr(2,7,1)=1 ; lxyz_arr(3,7,1)=0
   lxyz_arr(1,8,1)=4 ; lxyz_arr(2,8,1)=3 ; lxyz_arr(3,8,1)=0
   lxyz_arr(1,9,1)=2 ; lxyz_arr(2,9,1)=5 ; lxyz_arr(3,9,1)=0
   lxyz_arr(1,10,1)=4 ; lxyz_arr(2,10,1)=1 ; lxyz_arr(3,10,1)=2
   lxyz_arr(1,11,1)=2 ; lxyz_arr(2,11,1)=3 ; lxyz_arr(3,11,1)=2
   lxyz_arr(1,12,1)=2 ; lxyz_arr(2,12,1)=1 ; lxyz_arr(3,12,1)=4
   fac_arr(1,1)=0.2979934375117827994763496d0
   fac_arr(2,1)=0.3575921250141393593716196d0
   fac_arr(3,1)=0.05959868750235655989526993d0
   fac_arr(4,1)=0.3575921250141393593716196d0
   fac_arr(5,1)=0.1191973750047131197905399d0
   fac_arr(6,1)=0.05959868750235655989526993d0
   fac_arr(7,1)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(8,1)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(9,1)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(10,1)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(11,1)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(12,1)=-0.05959868750235655989526993/rhol**2d0
   lxyz_arr(1,1,2)=5 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=3 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=4 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=3 ; lxyz_arr(2,4,2)=0 ; lxyz_arr(3,4,2)=2
   lxyz_arr(1,5,2)=1 ; lxyz_arr(2,5,2)=2 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=0 ; lxyz_arr(3,6,2)=4
   lxyz_arr(1,7,2)=5 ; lxyz_arr(2,7,2)=2 ; lxyz_arr(3,7,2)=0
   lxyz_arr(1,8,2)=3 ; lxyz_arr(2,8,2)=4 ; lxyz_arr(3,8,2)=0
   lxyz_arr(1,9,2)=1 ; lxyz_arr(2,9,2)=6 ; lxyz_arr(3,9,2)=0
   lxyz_arr(1,10,2)=3 ; lxyz_arr(2,10,2)=2 ; lxyz_arr(3,10,2)=2
   lxyz_arr(1,11,2)=1 ; lxyz_arr(2,11,2)=4 ; lxyz_arr(3,11,2)=2
   lxyz_arr(1,12,2)=1 ; lxyz_arr(2,12,2)=2 ; lxyz_arr(3,12,2)=4
   fac_arr(1,2)=0.05959868750235655989526993d0
   fac_arr(2,2)=0.3575921250141393593716196d0
   fac_arr(3,2)=0.2979934375117827994763496d0
   fac_arr(4,2)=0.1191973750047131197905399d0
   fac_arr(5,2)=0.3575921250141393593716196d0
   fac_arr(6,2)=0.05959868750235655989526993d0
   fac_arr(7,2)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(8,2)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(9,2)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(10,2)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(11,2)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(12,2)=-0.05959868750235655989526993/rhol**2d0
   lxyz_arr(1,1,3)=3 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=1 ; lxyz_arr(3,3,3)=3
   lxyz_arr(1,4,3)=5 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=3 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=1 ; lxyz_arr(2,6,3)=5 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=3 ; lxyz_arr(2,7,3)=1 ; lxyz_arr(3,7,3)=3
   lxyz_arr(1,8,3)=1 ; lxyz_arr(2,8,3)=3 ; lxyz_arr(3,8,3)=3
   lxyz_arr(1,9,3)=1 ; lxyz_arr(2,9,3)=1 ; lxyz_arr(3,9,3)=5
   fac_arr(1,3)=0.2383947500094262395810797d0
   fac_arr(2,3)=0.2383947500094262395810797d0
   fac_arr(3,3)=0.2383947500094262395810797d0
   fac_arr(4,3)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(5,3)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(6,3)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(7,3)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(8,3)=-0.1191973750047131197905399/rhol**2d0
   fac_arr(9,3)=-0.05959868750235655989526993/rhol**2d0
else if (l.eq.3 .and. i.eq.3 .and. m.eq.4) then
   nterm_arr(1)=13
   nterm_arr(2)=13
   nterm_arr(3)=12
   lxyz_arr(1,1,1)=5 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=3 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=2
   lxyz_arr(1,5,1)=1 ; lxyz_arr(2,5,1)=0 ; lxyz_arr(3,5,1)=4
   lxyz_arr(1,6,1)=7 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=0
   lxyz_arr(1,7,1)=5 ; lxyz_arr(2,7,1)=2 ; lxyz_arr(3,7,1)=0
   lxyz_arr(1,8,1)=3 ; lxyz_arr(2,8,1)=4 ; lxyz_arr(3,8,1)=0
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=6 ; lxyz_arr(3,9,1)=0
   lxyz_arr(1,10,1)=5 ; lxyz_arr(2,10,1)=0 ; lxyz_arr(3,10,1)=2
   lxyz_arr(1,11,1)=1 ; lxyz_arr(2,11,1)=4 ; lxyz_arr(3,11,1)=2
   lxyz_arr(1,12,1)=3 ; lxyz_arr(2,12,1)=0 ; lxyz_arr(3,12,1)=4
   lxyz_arr(1,13,1)=1 ; lxyz_arr(2,13,1)=2 ; lxyz_arr(3,13,1)=4
   fac_arr(1,1)=0.1787960625070696796858098d0
   fac_arr(2,1)=0.1191973750047131197905399d0
   fac_arr(3,1)=-0.05959868750235655989526993d0
   fac_arr(4,1)=0.2383947500094262395810797d0
   fac_arr(5,1)=0.05959868750235655989526993d0
   fac_arr(6,1)=-0.02979934375117827994763496/rhol**2d0
   fac_arr(7,1)=-0.02979934375117827994763496/rhol**2d0
   fac_arr(8,1)=0.02979934375117827994763496/rhol**2d0
   fac_arr(9,1)=0.02979934375117827994763496/rhol**2d0
   fac_arr(10,1)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(11,1)=0.05959868750235655989526993/rhol**2d0
   fac_arr(12,1)=-0.02979934375117827994763496/rhol**2d0
   fac_arr(13,1)=0.02979934375117827994763496/rhol**2d0
   lxyz_arr(1,1,2)=4 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=5 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=3 ; lxyz_arr(3,4,2)=2
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=1 ; lxyz_arr(3,5,2)=4
   lxyz_arr(1,6,2)=6 ; lxyz_arr(2,6,2)=1 ; lxyz_arr(3,6,2)=0
   lxyz_arr(1,7,2)=4 ; lxyz_arr(2,7,2)=3 ; lxyz_arr(3,7,2)=0
   lxyz_arr(1,8,2)=2 ; lxyz_arr(2,8,2)=5 ; lxyz_arr(3,8,2)=0
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=7 ; lxyz_arr(3,9,2)=0
   lxyz_arr(1,10,2)=4 ; lxyz_arr(2,10,2)=1 ; lxyz_arr(3,10,2)=2
   lxyz_arr(1,11,2)=0 ; lxyz_arr(2,11,2)=5 ; lxyz_arr(3,11,2)=2
   lxyz_arr(1,12,2)=2 ; lxyz_arr(2,12,2)=1 ; lxyz_arr(3,12,2)=4
   lxyz_arr(1,13,2)=0 ; lxyz_arr(2,13,2)=3 ; lxyz_arr(3,13,2)=4
   fac_arr(1,2)=0.05959868750235655989526993d0
   fac_arr(2,2)=-0.1191973750047131197905399d0
   fac_arr(3,2)=-0.1787960625070696796858098d0
   fac_arr(4,2)=-0.2383947500094262395810797d0
   fac_arr(5,2)=-0.05959868750235655989526993d0
   fac_arr(6,2)=-0.02979934375117827994763496/rhol**2d0
   fac_arr(7,2)=-0.02979934375117827994763496/rhol**2d0
   fac_arr(8,2)=0.02979934375117827994763496/rhol**2d0
   fac_arr(9,2)=0.02979934375117827994763496/rhol**2d0
   fac_arr(10,2)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(11,2)=0.05959868750235655989526993/rhol**2d0
   fac_arr(12,2)=-0.02979934375117827994763496/rhol**2d0
   fac_arr(13,2)=0.02979934375117827994763496/rhol**2d0
   lxyz_arr(1,1,3)=4 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=4 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=2 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=3
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=2 ; lxyz_arr(3,4,3)=3
   lxyz_arr(1,5,3)=6 ; lxyz_arr(2,5,3)=0 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=4 ; lxyz_arr(2,6,3)=2 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=2 ; lxyz_arr(2,7,3)=4 ; lxyz_arr(3,7,3)=1
   lxyz_arr(1,8,3)=0 ; lxyz_arr(2,8,3)=6 ; lxyz_arr(3,8,3)=1
   lxyz_arr(1,9,3)=4 ; lxyz_arr(2,9,3)=0 ; lxyz_arr(3,9,3)=3
   lxyz_arr(1,10,3)=0 ; lxyz_arr(2,10,3)=4 ; lxyz_arr(3,10,3)=3
   lxyz_arr(1,11,3)=2 ; lxyz_arr(2,11,3)=0 ; lxyz_arr(3,11,3)=5
   lxyz_arr(1,12,3)=0 ; lxyz_arr(2,12,3)=2 ; lxyz_arr(3,12,3)=5
   fac_arr(1,3)=0.1191973750047131197905399d0
   fac_arr(2,3)=-0.1191973750047131197905399d0
   fac_arr(3,3)=0.1191973750047131197905399d0
   fac_arr(4,3)=-0.1191973750047131197905399d0
   fac_arr(5,3)=-0.02979934375117827994763496/rhol**2d0
   fac_arr(6,3)=-0.02979934375117827994763496/rhol**2d0
   fac_arr(7,3)=0.02979934375117827994763496/rhol**2d0
   fac_arr(8,3)=0.02979934375117827994763496/rhol**2d0
   fac_arr(9,3)=-0.05959868750235655989526993/rhol**2d0
   fac_arr(10,3)=0.05959868750235655989526993/rhol**2d0
   fac_arr(11,3)=-0.02979934375117827994763496/rhol**2d0
   fac_arr(12,3)=0.02979934375117827994763496/rhol**2d0
else if (l.eq.3 .and. i.eq.3 .and. m.eq.5) then
   nterm_arr(1)=11
   nterm_arr(2)=11
   nterm_arr(3)=10
   lxyz_arr(1,1,1)=5 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=4
   lxyz_arr(1,5,1)=7 ; lxyz_arr(2,5,1)=0 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=5 ; lxyz_arr(2,6,1)=2 ; lxyz_arr(3,6,1)=0
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=4 ; lxyz_arr(3,7,1)=0
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=6 ; lxyz_arr(3,8,1)=0
   lxyz_arr(1,9,1)=3 ; lxyz_arr(2,9,1)=0 ; lxyz_arr(3,9,1)=4
   lxyz_arr(1,10,1)=1 ; lxyz_arr(2,10,1)=2 ; lxyz_arr(3,10,1)=4
   lxyz_arr(1,11,1)=1 ; lxyz_arr(2,11,1)=0 ; lxyz_arr(3,11,1)=6
   fac_arr(1,1)=-0.1032279548185018340124748d0
   fac_arr(2,1)=-0.2064559096370036680249495d0
   fac_arr(3,1)=-0.1032279548185018340124748d0
   fac_arr(4,1)=0.1032279548185018340124748d0
   fac_arr(5,1)=0.01720465913641697233541246/rhol**2d0
   fac_arr(6,1)=0.05161397740925091700623738/rhol**2d0
   fac_arr(7,1)=0.05161397740925091700623738/rhol**2d0
   fac_arr(8,1)=0.01720465913641697233541246/rhol**2d0
   fac_arr(9,1)=-0.05161397740925091700623738/rhol**2d0
   fac_arr(10,1)=-0.05161397740925091700623738/rhol**2d0
   fac_arr(11,1)=-0.03440931827283394467082492/rhol**2d0
   lxyz_arr(1,1,2)=4 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=5 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=4
   lxyz_arr(1,5,2)=6 ; lxyz_arr(2,5,2)=1 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=4 ; lxyz_arr(2,6,2)=3 ; lxyz_arr(3,6,2)=0
   lxyz_arr(1,7,2)=2 ; lxyz_arr(2,7,2)=5 ; lxyz_arr(3,7,2)=0
   lxyz_arr(1,8,2)=0 ; lxyz_arr(2,8,2)=7 ; lxyz_arr(3,8,2)=0
   lxyz_arr(1,9,2)=2 ; lxyz_arr(2,9,2)=1 ; lxyz_arr(3,9,2)=4
   lxyz_arr(1,10,2)=0 ; lxyz_arr(2,10,2)=3 ; lxyz_arr(3,10,2)=4
   lxyz_arr(1,11,2)=0 ; lxyz_arr(2,11,2)=1 ; lxyz_arr(3,11,2)=6
   fac_arr(1,2)=-0.1032279548185018340124748d0
   fac_arr(2,2)=-0.2064559096370036680249495d0
   fac_arr(3,2)=-0.1032279548185018340124748d0
   fac_arr(4,2)=0.1032279548185018340124748d0
   fac_arr(5,2)=0.01720465913641697233541246/rhol**2d0
   fac_arr(6,2)=0.05161397740925091700623738/rhol**2d0
   fac_arr(7,2)=0.05161397740925091700623738/rhol**2d0
   fac_arr(8,2)=0.01720465913641697233541246/rhol**2d0
   fac_arr(9,2)=-0.05161397740925091700623738/rhol**2d0
   fac_arr(10,2)=-0.05161397740925091700623738/rhol**2d0
   fac_arr(11,2)=-0.03440931827283394467082492/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=3
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=3
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=5
   lxyz_arr(1,4,3)=6 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=4 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=2 ; lxyz_arr(2,6,3)=4 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=0 ; lxyz_arr(2,7,3)=6 ; lxyz_arr(3,7,3)=1
   lxyz_arr(1,8,3)=2 ; lxyz_arr(2,8,3)=0 ; lxyz_arr(3,8,3)=5
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=2 ; lxyz_arr(3,9,3)=5
   lxyz_arr(1,10,3)=0 ; lxyz_arr(2,10,3)=0 ; lxyz_arr(3,10,3)=7
   fac_arr(1,3)=0.2064559096370036680249495d0
   fac_arr(2,3)=0.2064559096370036680249495d0
   fac_arr(3,3)=0.2064559096370036680249495d0
   fac_arr(4,3)=0.01720465913641697233541246/rhol**2d0
   fac_arr(5,3)=0.05161397740925091700623738/rhol**2d0
   fac_arr(6,3)=0.05161397740925091700623738/rhol**2d0
   fac_arr(7,3)=0.01720465913641697233541246/rhol**2d0
   fac_arr(8,3)=-0.05161397740925091700623738/rhol**2d0
   fac_arr(9,3)=-0.05161397740925091700623738/rhol**2d0
   fac_arr(10,3)=-0.03440931827283394467082492/rhol**2d0
else if (l.eq.4 .and. i.eq.1 .and. m.eq.1) then
   nterm_arr(1)=6
   nterm_arr(2)=4
   nterm_arr(3)=4
   lxyz_arr(1,1,1)=2 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=0 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=2
   lxyz_arr(1,4,1)=4 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=2 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=2 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=2
   fac_arr(1,1)=0.9486832980505137995996681d0
   fac_arr(2,1)=0.3162277660168379331998894d0
   fac_arr(3,1)=-1.264911064067351732799557d0
   fac_arr(4,1)=-0.3162277660168379331998894/rhol**2d0
   fac_arr(5,1)=-0.3162277660168379331998894/rhol**2d0
   fac_arr(6,1)=1.264911064067351732799557/rhol**2d0
   lxyz_arr(1,1,2)=1 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=3 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=1 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=2
   fac_arr(1,2)=0.6324555320336758663997787d0
   fac_arr(2,2)=-0.3162277660168379331998894/rhol**2d0
   fac_arr(3,2)=-0.3162277660168379331998894/rhol**2d0
   fac_arr(4,2)=1.264911064067351732799557/rhol**2d0
   lxyz_arr(1,1,3)=1 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=3 ; lxyz_arr(2,2,3)=0 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=2 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=1 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=3
   fac_arr(1,3)=-2.529822128134703465599115d0
   fac_arr(2,3)=-0.3162277660168379331998894/rhol**2d0
   fac_arr(3,3)=-0.3162277660168379331998894/rhol**2d0
   fac_arr(4,3)=1.264911064067351732799557/rhol**2d0
else if (l.eq.4 .and. i.eq.1 .and. m.eq.2) then
   nterm_arr(1)=4
   nterm_arr(2)=6
   nterm_arr(3)=4
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=1 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=3 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=2
   fac_arr(1,1)=0.6324555320336758663997787d0
   fac_arr(2,1)=-0.3162277660168379331998894/rhol**2d0
   fac_arr(3,1)=-0.3162277660168379331998894/rhol**2d0
   fac_arr(4,1)=1.264911064067351732799557/rhol**2d0
   lxyz_arr(1,1,2)=2 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=0 ; lxyz_arr(3,3,2)=2
   lxyz_arr(1,4,2)=2 ; lxyz_arr(2,4,2)=2 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=4 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=2 ; lxyz_arr(3,6,2)=2
   fac_arr(1,2)=0.3162277660168379331998894d0
   fac_arr(2,2)=0.9486832980505137995996681d0
   fac_arr(3,2)=-1.264911064067351732799557d0
   fac_arr(4,2)=-0.3162277660168379331998894/rhol**2d0
   fac_arr(5,2)=-0.3162277660168379331998894/rhol**2d0
   fac_arr(6,2)=1.264911064067351732799557/rhol**2d0
   lxyz_arr(1,1,3)=0 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=1 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=3 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=3
   fac_arr(1,3)=-2.529822128134703465599115d0
   fac_arr(2,3)=-0.3162277660168379331998894/rhol**2d0
   fac_arr(3,3)=-0.3162277660168379331998894/rhol**2d0
   fac_arr(4,3)=1.264911064067351732799557/rhol**2d0
else if (l.eq.4 .and. i.eq.1 .and. m.eq.3) then
   nterm_arr(1)=4
   nterm_arr(2)=4
   nterm_arr(3)=6
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=2 ; lxyz_arr(3,3,1)=1
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=3
   fac_arr(1,1)=1.549193338482966754071706d0
   fac_arr(2,1)=-0.7745966692414833770358531/rhol**2d0
   fac_arr(3,1)=-0.7745966692414833770358531/rhol**2d0
   fac_arr(4,1)=0.5163977794943222513572354/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=1
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=3
   fac_arr(1,2)=1.549193338482966754071706d0
   fac_arr(2,2)=-0.7745966692414833770358531/rhol**2d0
   fac_arr(3,2)=-0.7745966692414833770358531/rhol**2d0
   fac_arr(4,2)=0.5163977794943222513572354/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=2
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=4
   fac_arr(1,3)=0.7745966692414833770358531d0
   fac_arr(2,3)=0.7745966692414833770358531d0
   fac_arr(3,3)=-1.549193338482966754071706d0
   fac_arr(4,3)=-0.7745966692414833770358531/rhol**2d0
   fac_arr(5,3)=-0.7745966692414833770358531/rhol**2d0
   fac_arr(6,3)=0.5163977794943222513572354/rhol**2d0
else if (l.eq.4 .and. i.eq.1 .and. m.eq.4) then
   nterm_arr(1)=4
   nterm_arr(2)=3
   nterm_arr(3)=2
   lxyz_arr(1,1,1)=2 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=0 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=4 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=2 ; lxyz_arr(2,4,1)=2 ; lxyz_arr(3,4,1)=0
   fac_arr(1,1)=1.224744871391589049098642d0
   fac_arr(2,1)=-1.224744871391589049098642d0
   fac_arr(3,1)=-0.408248290463863016366214/rhol**2d0
   fac_arr(4,1)=1.224744871391589049098642/rhol**2d0
   lxyz_arr(1,1,2)=1 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=3 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=0
   fac_arr(1,2)=-2.449489742783178098197284d0
   fac_arr(2,2)=-0.408248290463863016366214/rhol**2d0
   fac_arr(3,2)=1.224744871391589049098642/rhol**2d0
   lxyz_arr(1,1,3)=3 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   fac_arr(1,3)=-0.408248290463863016366214/rhol**2d0
   fac_arr(2,3)=1.224744871391589049098642/rhol**2d0
else if (l.eq.4 .and. i.eq.1 .and. m.eq.5) then
   nterm_arr(1)=3
   nterm_arr(2)=4
   nterm_arr(3)=2
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=1 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=3 ; lxyz_arr(3,3,1)=0
   fac_arr(1,1)=-2.449489742783178098197284d0
   fac_arr(2,1)=1.224744871391589049098642/rhol**2d0
   fac_arr(3,1)=-0.408248290463863016366214/rhol**2d0
   lxyz_arr(1,1,2)=2 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=2 ; lxyz_arr(2,3,2)=2 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=4 ; lxyz_arr(3,4,2)=0
   fac_arr(1,2)=-1.224744871391589049098642d0
   fac_arr(2,2)=1.224744871391589049098642d0
   fac_arr(3,2)=1.224744871391589049098642/rhol**2d0
   fac_arr(4,2)=-0.408248290463863016366214/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=1
   fac_arr(1,3)=1.224744871391589049098642/rhol**2d0
   fac_arr(2,3)=-0.408248290463863016366214/rhol**2d0
else if (l.eq.4 .and. i.eq.1 .and. m.eq.6) then
   nterm_arr(1)=3
   nterm_arr(2)=3
   nterm_arr(3)=4
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=2 ; lxyz_arr(3,3,1)=1
   fac_arr(1,1)=2.d0
   fac_arr(2,1)=-1./rhol**2d0
   fac_arr(3,1)=rhol**(-2)
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=1
   fac_arr(1,2)=-2.d0
   fac_arr(2,2)=-1./rhol**2d0
   fac_arr(3,2)=rhol**(-2)
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=2 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=2
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=2 ; lxyz_arr(3,4,3)=2
   fac_arr(1,3)=1.d0
   fac_arr(2,3)=-1.d0
   fac_arr(3,3)=-1./rhol**2d0
   fac_arr(4,3)=rhol**(-2)
else if (l.eq.4 .and. i.eq.1 .and. m.eq.7) then
   nterm_arr(1)=2
   nterm_arr(2)=2
   nterm_arr(3)=2
   lxyz_arr(1,1,1)=0 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=1 ; lxyz_arr(3,2,1)=1
   fac_arr(1,1)=2.d0
   fac_arr(2,1)=-2./rhol**2d0
   lxyz_arr(1,1,2)=1 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=1 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=1
   fac_arr(1,2)=2.d0
   fac_arr(2,2)=-2./rhol**2d0
   lxyz_arr(1,1,3)=1 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=1 ; lxyz_arr(3,2,3)=2
   fac_arr(1,3)=2.d0
   fac_arr(2,3)=-2./rhol**2d0
else if (l.eq.4 .and. i.eq.2 .and. m.eq.1) then
   nterm_arr(1)=12
   nterm_arr(2)=9
   nterm_arr(3)=9
   lxyz_arr(1,1,1)=4 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=2 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=2
   lxyz_arr(1,5,1)=0 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=0 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=4
   lxyz_arr(1,7,1)=6 ; lxyz_arr(2,7,1)=0 ; lxyz_arr(3,7,1)=0
   lxyz_arr(1,8,1)=4 ; lxyz_arr(2,8,1)=2 ; lxyz_arr(3,8,1)=0
   lxyz_arr(1,9,1)=2 ; lxyz_arr(2,9,1)=4 ; lxyz_arr(3,9,1)=0
   lxyz_arr(1,10,1)=4 ; lxyz_arr(2,10,1)=0 ; lxyz_arr(3,10,1)=2
   lxyz_arr(1,11,1)=2 ; lxyz_arr(2,11,1)=2 ; lxyz_arr(3,11,1)=2
   lxyz_arr(1,12,1)=2 ; lxyz_arr(2,12,1)=0 ; lxyz_arr(3,12,1)=4
   fac_arr(1,1)=0.3178208630818641051489253d0
   fac_arr(2,1)=0.3813850356982369261787104d0
   fac_arr(3,1)=0.06356417261637282102978506d0
   fac_arr(4,1)=-0.5720775535473553892680656d0
   fac_arr(5,1)=-0.1906925178491184630893552d0
   fac_arr(6,1)=-0.2542566904654912841191402d0
   fac_arr(7,1)=-0.06356417261637282102978506/rhol**2d0
   fac_arr(8,1)=-0.1271283452327456420595701/rhol**2d0
   fac_arr(9,1)=-0.06356417261637282102978506/rhol**2d0
   fac_arr(10,1)=0.1906925178491184630893552/rhol**2d0
   fac_arr(11,1)=0.1906925178491184630893552/rhol**2d0
   fac_arr(12,1)=0.2542566904654912841191402/rhol**2d0
   lxyz_arr(1,1,2)=3 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=1 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=2
   lxyz_arr(1,4,2)=5 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=3 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=5 ; lxyz_arr(3,6,2)=0
   lxyz_arr(1,7,2)=3 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=2
   lxyz_arr(1,8,2)=1 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=2
   lxyz_arr(1,9,2)=1 ; lxyz_arr(2,9,2)=1 ; lxyz_arr(3,9,2)=4
   fac_arr(1,2)=0.2542566904654912841191402d0
   fac_arr(2,2)=0.2542566904654912841191402d0
   fac_arr(3,2)=-0.3813850356982369261787104d0
   fac_arr(4,2)=-0.06356417261637282102978506/rhol**2d0
   fac_arr(5,2)=-0.1271283452327456420595701/rhol**2d0
   fac_arr(6,2)=-0.06356417261637282102978506/rhol**2d0
   fac_arr(7,2)=0.1906925178491184630893552/rhol**2d0
   fac_arr(8,2)=0.1906925178491184630893552/rhol**2d0
   fac_arr(9,2)=0.2542566904654912841191402/rhol**2d0
   lxyz_arr(1,1,3)=3 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=3
   lxyz_arr(1,4,3)=5 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=3 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=1 ; lxyz_arr(2,6,3)=4 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=3 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=3
   lxyz_arr(1,8,3)=1 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=3
   lxyz_arr(1,9,3)=1 ; lxyz_arr(2,9,3)=0 ; lxyz_arr(3,9,3)=5
   fac_arr(1,3)=-0.3813850356982369261787104d0
   fac_arr(2,3)=-0.3813850356982369261787104d0
   fac_arr(3,3)=-1.017026761861965136476561d0
   fac_arr(4,3)=-0.06356417261637282102978506/rhol**2d0
   fac_arr(5,3)=-0.1271283452327456420595701/rhol**2d0
   fac_arr(6,3)=-0.06356417261637282102978506/rhol**2d0
   fac_arr(7,3)=0.1906925178491184630893552/rhol**2d0
   fac_arr(8,3)=0.1906925178491184630893552/rhol**2d0
   fac_arr(9,3)=0.2542566904654912841191402/rhol**2d0
else if (l.eq.4 .and. i.eq.2 .and. m.eq.2) then
   nterm_arr(1)=9
   nterm_arr(2)=12
   nterm_arr(3)=9
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=1 ; lxyz_arr(3,3,1)=2
   lxyz_arr(1,4,1)=5 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=5 ; lxyz_arr(3,6,1)=0
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=1 ; lxyz_arr(3,7,1)=2
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=3 ; lxyz_arr(3,8,1)=2
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=1 ; lxyz_arr(3,9,1)=4
   fac_arr(1,1)=0.2542566904654912841191402d0
   fac_arr(2,1)=0.2542566904654912841191402d0
   fac_arr(3,1)=-0.3813850356982369261787104d0
   fac_arr(4,1)=-0.06356417261637282102978506/rhol**2d0
   fac_arr(5,1)=-0.1271283452327456420595701/rhol**2d0
   fac_arr(6,1)=-0.06356417261637282102978506/rhol**2d0
   fac_arr(7,1)=0.1906925178491184630893552/rhol**2d0
   fac_arr(8,1)=0.1906925178491184630893552/rhol**2d0
   fac_arr(9,1)=0.2542566904654912841191402/rhol**2d0
   lxyz_arr(1,1,2)=4 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=4 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=2 ; lxyz_arr(2,4,2)=0 ; lxyz_arr(3,4,2)=2
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=2 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=0 ; lxyz_arr(3,6,2)=4
   lxyz_arr(1,7,2)=4 ; lxyz_arr(2,7,2)=2 ; lxyz_arr(3,7,2)=0
   lxyz_arr(1,8,2)=2 ; lxyz_arr(2,8,2)=4 ; lxyz_arr(3,8,2)=0
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=6 ; lxyz_arr(3,9,2)=0
   lxyz_arr(1,10,2)=2 ; lxyz_arr(2,10,2)=2 ; lxyz_arr(3,10,2)=2
   lxyz_arr(1,11,2)=0 ; lxyz_arr(2,11,2)=4 ; lxyz_arr(3,11,2)=2
   lxyz_arr(1,12,2)=0 ; lxyz_arr(2,12,2)=2 ; lxyz_arr(3,12,2)=4
   fac_arr(1,2)=0.06356417261637282102978506d0
   fac_arr(2,2)=0.3813850356982369261787104d0
   fac_arr(3,2)=0.3178208630818641051489253d0
   fac_arr(4,2)=-0.1906925178491184630893552d0
   fac_arr(5,2)=-0.5720775535473553892680656d0
   fac_arr(6,2)=-0.2542566904654912841191402d0
   fac_arr(7,2)=-0.06356417261637282102978506/rhol**2d0
   fac_arr(8,2)=-0.1271283452327456420595701/rhol**2d0
   fac_arr(9,2)=-0.06356417261637282102978506/rhol**2d0
   fac_arr(10,2)=0.1906925178491184630893552/rhol**2d0
   fac_arr(11,2)=0.1906925178491184630893552/rhol**2d0
   fac_arr(12,2)=0.2542566904654912841191402/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=1 ; lxyz_arr(3,3,3)=3
   lxyz_arr(1,4,3)=4 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=2 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=5 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=2 ; lxyz_arr(2,7,3)=1 ; lxyz_arr(3,7,3)=3
   lxyz_arr(1,8,3)=0 ; lxyz_arr(2,8,3)=3 ; lxyz_arr(3,8,3)=3
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=1 ; lxyz_arr(3,9,3)=5
   fac_arr(1,3)=-0.3813850356982369261787104d0
   fac_arr(2,3)=-0.3813850356982369261787104d0
   fac_arr(3,3)=-1.017026761861965136476561d0
   fac_arr(4,3)=-0.06356417261637282102978506/rhol**2d0
   fac_arr(5,3)=-0.1271283452327456420595701/rhol**2d0
   fac_arr(6,3)=-0.06356417261637282102978506/rhol**2d0
   fac_arr(7,3)=0.1906925178491184630893552/rhol**2d0
   fac_arr(8,3)=0.1906925178491184630893552/rhol**2d0
   fac_arr(9,3)=0.2542566904654912841191402/rhol**2d0
else if (l.eq.4 .and. i.eq.2 .and. m.eq.3) then
   nterm_arr(1)=9
   nterm_arr(2)=9
   nterm_arr(3)=12
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=3
   lxyz_arr(1,4,1)=5 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=1
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=1
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=4 ; lxyz_arr(3,6,1)=1
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=0 ; lxyz_arr(3,7,1)=3
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=2 ; lxyz_arr(3,8,1)=3
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=0 ; lxyz_arr(3,9,1)=5
   fac_arr(1,1)=0.6227991553292183767329405d0
   fac_arr(2,1)=0.6227991553292183767329405d0
   fac_arr(3,1)=0.1037998592215363961221568d0
   fac_arr(4,1)=-0.1556997888323045941832351/rhol**2d0
   fac_arr(5,1)=-0.3113995776646091883664703/rhol**2d0
   fac_arr(6,1)=-0.1556997888323045941832351/rhol**2d0
   fac_arr(7,1)=-0.05189992961076819806107838/rhol**2d0
   fac_arr(8,1)=-0.05189992961076819806107838/rhol**2d0
   fac_arr(9,1)=0.1037998592215363961221568/rhol**2d0
   lxyz_arr(1,1,2)=2 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=3
   lxyz_arr(1,4,2)=4 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=1
   lxyz_arr(1,5,2)=2 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=1
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=5 ; lxyz_arr(3,6,2)=1
   lxyz_arr(1,7,2)=2 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=3
   lxyz_arr(1,8,2)=0 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=3
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=1 ; lxyz_arr(3,9,2)=5
   fac_arr(1,2)=0.6227991553292183767329405d0
   fac_arr(2,2)=0.6227991553292183767329405d0
   fac_arr(3,2)=0.1037998592215363961221568d0
   fac_arr(4,2)=-0.1556997888323045941832351/rhol**2d0
   fac_arr(5,2)=-0.3113995776646091883664703/rhol**2d0
   fac_arr(6,2)=-0.1556997888323045941832351/rhol**2d0
   fac_arr(7,2)=-0.05189992961076819806107838/rhol**2d0
   fac_arr(8,2)=-0.05189992961076819806107838/rhol**2d0
   fac_arr(9,2)=0.1037998592215363961221568/rhol**2d0
   lxyz_arr(1,1,3)=4 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=4 ; lxyz_arr(3,3,3)=0
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=4
   lxyz_arr(1,7,3)=4 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=2
   lxyz_arr(1,8,3)=2 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=2
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=4 ; lxyz_arr(3,9,3)=2
   lxyz_arr(1,10,3)=2 ; lxyz_arr(2,10,3)=0 ; lxyz_arr(3,10,3)=4
   lxyz_arr(1,11,3)=0 ; lxyz_arr(2,11,3)=2 ; lxyz_arr(3,11,3)=4
   lxyz_arr(1,12,3)=0 ; lxyz_arr(2,12,3)=0 ; lxyz_arr(3,12,3)=6
   fac_arr(1,3)=0.1556997888323045941832351d0
   fac_arr(2,3)=0.3113995776646091883664703d0
   fac_arr(3,3)=0.1556997888323045941832351d0
   fac_arr(4,3)=0.1556997888323045941832351d0
   fac_arr(5,3)=0.1556997888323045941832351d0
   fac_arr(6,3)=-0.5189992961076819806107838d0
   fac_arr(7,3)=-0.1556997888323045941832351/rhol**2d0
   fac_arr(8,3)=-0.3113995776646091883664703/rhol**2d0
   fac_arr(9,3)=-0.1556997888323045941832351/rhol**2d0
   fac_arr(10,3)=-0.05189992961076819806107838/rhol**2d0
   fac_arr(11,3)=-0.05189992961076819806107838/rhol**2d0
   fac_arr(12,3)=0.1037998592215363961221568/rhol**2d0
else if (l.eq.4 .and. i.eq.2 .and. m.eq.4) then
   nterm_arr(1)=10
   nterm_arr(2)=8
   nterm_arr(3)=7
   lxyz_arr(1,1,1)=4 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=2 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=2
   lxyz_arr(1,5,1)=0 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=6 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=0
   lxyz_arr(1,7,1)=4 ; lxyz_arr(2,7,1)=2 ; lxyz_arr(3,7,1)=0
   lxyz_arr(1,8,1)=2 ; lxyz_arr(2,8,1)=4 ; lxyz_arr(3,8,1)=0
   lxyz_arr(1,9,1)=4 ; lxyz_arr(2,9,1)=0 ; lxyz_arr(3,9,1)=2
   lxyz_arr(1,10,1)=2 ; lxyz_arr(2,10,1)=2 ; lxyz_arr(3,10,1)=2
   fac_arr(1,1)=0.4103049699311091091141355d0
   fac_arr(2,1)=-0.4923659639173309309369626d0
   fac_arr(3,1)=-0.2461829819586654654684813d0
   fac_arr(4,1)=0.2461829819586654654684813d0
   fac_arr(5,1)=-0.2461829819586654654684813d0
   fac_arr(6,1)=-0.08206099398622182182282711/rhol**2d0
   fac_arr(7,1)=0.1641219879724436436456542/rhol**2d0
   fac_arr(8,1)=0.2461829819586654654684813/rhol**2d0
   fac_arr(9,1)=-0.08206099398622182182282711/rhol**2d0
   fac_arr(10,1)=0.2461829819586654654684813/rhol**2d0
   lxyz_arr(1,1,2)=3 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=1 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=2
   lxyz_arr(1,4,2)=5 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=3 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=5 ; lxyz_arr(3,6,2)=0
   lxyz_arr(1,7,2)=3 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=2
   lxyz_arr(1,8,2)=1 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=2
   fac_arr(1,2)=-0.3282439759448872872913084d0
   fac_arr(2,2)=-0.9847319278346618618739253d0
   fac_arr(3,2)=-0.4923659639173309309369626d0
   fac_arr(4,2)=-0.08206099398622182182282711/rhol**2d0
   fac_arr(5,2)=0.1641219879724436436456542/rhol**2d0
   fac_arr(6,2)=0.2461829819586654654684813/rhol**2d0
   fac_arr(7,2)=-0.08206099398622182182282711/rhol**2d0
   fac_arr(8,2)=0.2461829819586654654684813/rhol**2d0
   lxyz_arr(1,1,3)=3 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=5 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=3 ; lxyz_arr(2,4,3)=2 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=1 ; lxyz_arr(2,5,3)=4 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=3 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=3
   lxyz_arr(1,7,3)=1 ; lxyz_arr(2,7,3)=2 ; lxyz_arr(3,7,3)=3
   fac_arr(1,3)=0.1641219879724436436456542d0
   fac_arr(2,3)=-0.4923659639173309309369626d0
   fac_arr(3,3)=-0.08206099398622182182282711/rhol**2d0
   fac_arr(4,3)=0.1641219879724436436456542/rhol**2d0
   fac_arr(5,3)=0.2461829819586654654684813/rhol**2d0
   fac_arr(6,3)=-0.08206099398622182182282711/rhol**2d0
   fac_arr(7,3)=0.2461829819586654654684813/rhol**2d0
else if (l.eq.4 .and. i.eq.2 .and. m.eq.5) then
   nterm_arr(1)=8
   nterm_arr(2)=10
   nterm_arr(3)=7
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=1 ; lxyz_arr(3,3,1)=2
   lxyz_arr(1,4,1)=5 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=5 ; lxyz_arr(3,6,1)=0
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=1 ; lxyz_arr(3,7,1)=2
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=3 ; lxyz_arr(3,8,1)=2
   fac_arr(1,1)=-0.9847319278346618618739253d0
   fac_arr(2,1)=-0.3282439759448872872913084d0
   fac_arr(3,1)=-0.4923659639173309309369626d0
   fac_arr(4,1)=0.2461829819586654654684813/rhol**2d0
   fac_arr(5,1)=0.1641219879724436436456542/rhol**2d0
   fac_arr(6,1)=-0.08206099398622182182282711/rhol**2d0
   fac_arr(7,1)=0.2461829819586654654684813/rhol**2d0
   fac_arr(8,1)=-0.08206099398622182182282711/rhol**2d0
   lxyz_arr(1,1,2)=4 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=4 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=2 ; lxyz_arr(2,4,2)=0 ; lxyz_arr(3,4,2)=2
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=2 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=4 ; lxyz_arr(2,6,2)=2 ; lxyz_arr(3,6,2)=0
   lxyz_arr(1,7,2)=2 ; lxyz_arr(2,7,2)=4 ; lxyz_arr(3,7,2)=0
   lxyz_arr(1,8,2)=0 ; lxyz_arr(2,8,2)=6 ; lxyz_arr(3,8,2)=0
   lxyz_arr(1,9,2)=2 ; lxyz_arr(2,9,2)=2 ; lxyz_arr(3,9,2)=2
   lxyz_arr(1,10,2)=0 ; lxyz_arr(2,10,2)=4 ; lxyz_arr(3,10,2)=2
   fac_arr(1,2)=-0.2461829819586654654684813d0
   fac_arr(2,2)=-0.4923659639173309309369626d0
   fac_arr(3,2)=0.4103049699311091091141355d0
   fac_arr(4,2)=-0.2461829819586654654684813d0
   fac_arr(5,2)=0.2461829819586654654684813d0
   fac_arr(6,2)=0.2461829819586654654684813/rhol**2d0
   fac_arr(7,2)=0.1641219879724436436456542/rhol**2d0
   fac_arr(8,2)=-0.08206099398622182182282711/rhol**2d0
   fac_arr(9,2)=0.2461829819586654654684813/rhol**2d0
   fac_arr(10,2)=-0.08206099398622182182282711/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=4 ; lxyz_arr(2,3,3)=1 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=3 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=5 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=2 ; lxyz_arr(2,6,3)=1 ; lxyz_arr(3,6,3)=3
   lxyz_arr(1,7,3)=0 ; lxyz_arr(2,7,3)=3 ; lxyz_arr(3,7,3)=3
   fac_arr(1,3)=-0.4923659639173309309369626d0
   fac_arr(2,3)=0.1641219879724436436456542d0
   fac_arr(3,3)=0.2461829819586654654684813/rhol**2d0
   fac_arr(4,3)=0.1641219879724436436456542/rhol**2d0
   fac_arr(5,3)=-0.08206099398622182182282711/rhol**2d0
   fac_arr(6,3)=0.2461829819586654654684813/rhol**2d0
   fac_arr(7,3)=-0.08206099398622182182282711/rhol**2d0
else if (l.eq.4 .and. i.eq.2 .and. m.eq.6) then
   nterm_arr(1)=6
   nterm_arr(2)=6
   nterm_arr(3)=8
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=3
   lxyz_arr(1,3,1)=5 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=1
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=4 ; lxyz_arr(3,4,1)=1
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=0 ; lxyz_arr(3,5,1)=3
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=2 ; lxyz_arr(3,6,1)=3
   fac_arr(1,1)=0.8040302522073696603914988d0
   fac_arr(2,1)=0.4020151261036848301957494d0
   fac_arr(3,1)=-0.2010075630518424150978747/rhol**2d0
   fac_arr(4,1)=0.2010075630518424150978747/rhol**2d0
   fac_arr(5,1)=-0.2010075630518424150978747/rhol**2d0
   fac_arr(6,1)=0.2010075630518424150978747/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=3 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=3
   lxyz_arr(1,3,2)=4 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=1
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=5 ; lxyz_arr(3,4,2)=1
   lxyz_arr(1,5,2)=2 ; lxyz_arr(2,5,2)=1 ; lxyz_arr(3,5,2)=3
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=3 ; lxyz_arr(3,6,2)=3
   fac_arr(1,2)=-0.8040302522073696603914988d0
   fac_arr(2,2)=-0.4020151261036848301957494d0
   fac_arr(3,2)=-0.2010075630518424150978747/rhol**2d0
   fac_arr(4,2)=0.2010075630518424150978747/rhol**2d0
   fac_arr(5,2)=-0.2010075630518424150978747/rhol**2d0
   fac_arr(6,2)=0.2010075630518424150978747/rhol**2d0
   lxyz_arr(1,1,3)=4 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=4 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=2 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=2
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=2 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=4 ; lxyz_arr(2,5,3)=0 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=4 ; lxyz_arr(3,6,3)=2
   lxyz_arr(1,7,3)=2 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=4
   lxyz_arr(1,8,3)=0 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=4
   fac_arr(1,3)=0.2010075630518424150978747d0
   fac_arr(2,3)=-0.2010075630518424150978747d0
   fac_arr(3,3)=0.6030226891555272452936241d0
   fac_arr(4,3)=-0.6030226891555272452936241d0
   fac_arr(5,3)=-0.2010075630518424150978747/rhol**2d0
   fac_arr(6,3)=0.2010075630518424150978747/rhol**2d0
   fac_arr(7,3)=-0.2010075630518424150978747/rhol**2d0
   fac_arr(8,3)=0.2010075630518424150978747/rhol**2d0
else if (l.eq.4 .and. i.eq.2 .and. m.eq.7) then
   nterm_arr(1)=6
   nterm_arr(2)=6
   nterm_arr(3)=6
   lxyz_arr(1,1,1)=2 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=0 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=1 ; lxyz_arr(3,3,1)=3
   lxyz_arr(1,4,1)=4 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=1
   lxyz_arr(1,5,1)=2 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=1
   lxyz_arr(1,6,1)=2 ; lxyz_arr(2,6,1)=1 ; lxyz_arr(3,6,1)=3
   fac_arr(1,1)=1.206045378311054490587248d0
   fac_arr(2,1)=0.4020151261036848301957494d0
   fac_arr(3,1)=0.4020151261036848301957494d0
   fac_arr(4,1)=-0.4020151261036848301957494/rhol**2d0
   fac_arr(5,1)=-0.4020151261036848301957494/rhol**2d0
   fac_arr(6,1)=-0.4020151261036848301957494/rhol**2d0
   lxyz_arr(1,1,2)=3 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=1 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=0 ; lxyz_arr(3,3,2)=3
   lxyz_arr(1,4,2)=3 ; lxyz_arr(2,4,2)=2 ; lxyz_arr(3,4,2)=1
   lxyz_arr(1,5,2)=1 ; lxyz_arr(2,5,2)=4 ; lxyz_arr(3,5,2)=1
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=2 ; lxyz_arr(3,6,2)=3
   fac_arr(1,2)=0.4020151261036848301957494d0
   fac_arr(2,2)=1.206045378311054490587248d0
   fac_arr(3,2)=0.4020151261036848301957494d0
   fac_arr(4,2)=-0.4020151261036848301957494/rhol**2d0
   fac_arr(5,2)=-0.4020151261036848301957494/rhol**2d0
   fac_arr(6,2)=-0.4020151261036848301957494/rhol**2d0
   lxyz_arr(1,1,3)=3 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=1 ; lxyz_arr(3,3,3)=2
   lxyz_arr(1,4,3)=3 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=1 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=1 ; lxyz_arr(2,6,3)=1 ; lxyz_arr(3,6,3)=4
   fac_arr(1,3)=0.4020151261036848301957494d0
   fac_arr(2,3)=0.4020151261036848301957494d0
   fac_arr(3,3)=1.206045378311054490587248d0
   fac_arr(4,3)=-0.4020151261036848301957494/rhol**2d0
   fac_arr(5,3)=-0.4020151261036848301957494/rhol**2d0
   fac_arr(6,3)=-0.4020151261036848301957494/rhol**2d0
else if (l.eq.4 .and. i.eq.3 .and. m.eq.1) then
   nterm_arr(1)=20
   nterm_arr(2)=16
   nterm_arr(3)=16
   lxyz_arr(1,1,1)=6 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=4 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=2 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=0 ; lxyz_arr(2,4,1)=6 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=4 ; lxyz_arr(2,5,1)=0 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=2 ; lxyz_arr(2,6,1)=2 ; lxyz_arr(3,6,1)=2
   lxyz_arr(1,7,1)=0 ; lxyz_arr(2,7,1)=4 ; lxyz_arr(3,7,1)=2
   lxyz_arr(1,8,1)=2 ; lxyz_arr(2,8,1)=0 ; lxyz_arr(3,8,1)=4
   lxyz_arr(1,9,1)=0 ; lxyz_arr(2,9,1)=2 ; lxyz_arr(3,9,1)=4
   lxyz_arr(1,10,1)=0 ; lxyz_arr(2,10,1)=0 ; lxyz_arr(3,10,1)=6
   lxyz_arr(1,11,1)=8 ; lxyz_arr(2,11,1)=0 ; lxyz_arr(3,11,1)=0
   lxyz_arr(1,12,1)=6 ; lxyz_arr(2,12,1)=2 ; lxyz_arr(3,12,1)=0
   lxyz_arr(1,13,1)=4 ; lxyz_arr(2,13,1)=4 ; lxyz_arr(3,13,1)=0
   lxyz_arr(1,14,1)=2 ; lxyz_arr(2,14,1)=6 ; lxyz_arr(3,14,1)=0
   lxyz_arr(1,15,1)=6 ; lxyz_arr(2,15,1)=0 ; lxyz_arr(3,15,1)=2
   lxyz_arr(1,16,1)=4 ; lxyz_arr(2,16,1)=2 ; lxyz_arr(3,16,1)=2
   lxyz_arr(1,17,1)=2 ; lxyz_arr(2,17,1)=4 ; lxyz_arr(3,17,1)=2
   lxyz_arr(1,18,1)=4 ; lxyz_arr(2,18,1)=0 ; lxyz_arr(3,18,1)=4
   lxyz_arr(1,19,1)=2 ; lxyz_arr(2,19,1)=2 ; lxyz_arr(3,19,1)=4
   lxyz_arr(1,20,1)=2 ; lxyz_arr(2,20,1)=0 ; lxyz_arr(3,20,1)=6
   fac_arr(1,1)=0.06372694925323242808889581d0
   fac_arr(2,1)=0.1365577483997837744762053d0
   fac_arr(3,1)=0.08193464903987026468572318d0
   fac_arr(4,1)=0.009103849893318918298413687d0
   fac_arr(5,1)=-0.09103849893318918298413687d0
   fac_arr(6,1)=-0.1092461987198270195809642d0
   fac_arr(7,1)=-0.01820769978663783659682737d0
   fac_arr(8,1)=-0.1911808477596972842666874d0
   fac_arr(9,1)=-0.06372694925323242808889581d0
   fac_arr(10,1)=-0.03641539957327567319365475d0
   fac_arr(11,1)=-0.009103849893318918298413687/rhol**2d0
   fac_arr(12,1)=-0.02731154967995675489524106/rhol**2d0
   fac_arr(13,1)=-0.02731154967995675489524106/rhol**2d0
   fac_arr(14,1)=-0.009103849893318918298413687/rhol**2d0
   fac_arr(15,1)=0.01820769978663783659682737/rhol**2d0
   fac_arr(16,1)=0.03641539957327567319365475/rhol**2d0
   fac_arr(17,1)=0.01820769978663783659682737/rhol**2d0
   fac_arr(18,1)=0.06372694925323242808889581/rhol**2d0
   fac_arr(19,1)=0.06372694925323242808889581/rhol**2d0
   fac_arr(20,1)=0.03641539957327567319365475/rhol**2d0
   lxyz_arr(1,1,2)=5 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=3 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=5 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=3 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=2
   lxyz_arr(1,5,2)=1 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=1 ; lxyz_arr(3,6,2)=4
   lxyz_arr(1,7,2)=7 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=0
   lxyz_arr(1,8,2)=5 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=0
   lxyz_arr(1,9,2)=3 ; lxyz_arr(2,9,2)=5 ; lxyz_arr(3,9,2)=0
   lxyz_arr(1,10,2)=1 ; lxyz_arr(2,10,2)=7 ; lxyz_arr(3,10,2)=0
   lxyz_arr(1,11,2)=5 ; lxyz_arr(2,11,2)=1 ; lxyz_arr(3,11,2)=2
   lxyz_arr(1,12,2)=3 ; lxyz_arr(2,12,2)=3 ; lxyz_arr(3,12,2)=2
   lxyz_arr(1,13,2)=1 ; lxyz_arr(2,13,2)=5 ; lxyz_arr(3,13,2)=2
   lxyz_arr(1,14,2)=3 ; lxyz_arr(2,14,2)=1 ; lxyz_arr(3,14,2)=4
   lxyz_arr(1,15,2)=1 ; lxyz_arr(2,15,2)=3 ; lxyz_arr(3,15,2)=4
   lxyz_arr(1,16,2)=1 ; lxyz_arr(2,16,2)=1 ; lxyz_arr(3,16,2)=6
   fac_arr(1,2)=0.05462309935991350979048212d0
   fac_arr(2,2)=0.1092461987198270195809642d0
   fac_arr(3,2)=0.05462309935991350979048212d0
   fac_arr(4,2)=-0.0728307991465513463873095d0
   fac_arr(5,2)=-0.0728307991465513463873095d0
   fac_arr(6,2)=-0.1274538985064648561777916d0
   fac_arr(7,2)=-0.009103849893318918298413687/rhol**2d0
   fac_arr(8,2)=-0.02731154967995675489524106/rhol**2d0
   fac_arr(9,2)=-0.02731154967995675489524106/rhol**2d0
   fac_arr(10,2)=-0.009103849893318918298413687/rhol**2d0
   fac_arr(11,2)=0.01820769978663783659682737/rhol**2d0
   fac_arr(12,2)=0.03641539957327567319365475/rhol**2d0
   fac_arr(13,2)=0.01820769978663783659682737/rhol**2d0
   fac_arr(14,2)=0.06372694925323242808889581/rhol**2d0
   fac_arr(15,2)=0.06372694925323242808889581/rhol**2d0
   fac_arr(16,2)=0.03641539957327567319365475/rhol**2d0
   lxyz_arr(1,1,3)=5 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=3 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=4 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=3 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=3
   lxyz_arr(1,5,3)=1 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=3
   lxyz_arr(1,6,3)=1 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=5
   lxyz_arr(1,7,3)=7 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=1
   lxyz_arr(1,8,3)=5 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=1
   lxyz_arr(1,9,3)=3 ; lxyz_arr(2,9,3)=4 ; lxyz_arr(3,9,3)=1
   lxyz_arr(1,10,3)=1 ; lxyz_arr(2,10,3)=6 ; lxyz_arr(3,10,3)=1
   lxyz_arr(1,11,3)=5 ; lxyz_arr(2,11,3)=0 ; lxyz_arr(3,11,3)=3
   lxyz_arr(1,12,3)=3 ; lxyz_arr(2,12,3)=2 ; lxyz_arr(3,12,3)=3
   lxyz_arr(1,13,3)=1 ; lxyz_arr(2,13,3)=4 ; lxyz_arr(3,13,3)=3
   lxyz_arr(1,14,3)=3 ; lxyz_arr(2,14,3)=0 ; lxyz_arr(3,14,3)=5
   lxyz_arr(1,15,3)=1 ; lxyz_arr(2,15,3)=2 ; lxyz_arr(3,15,3)=5
   lxyz_arr(1,16,3)=1 ; lxyz_arr(2,16,3)=0 ; lxyz_arr(3,16,3)=7
   fac_arr(1,3)=-0.03641539957327567319365475d0
   fac_arr(2,3)=-0.0728307991465513463873095d0
   fac_arr(3,3)=-0.03641539957327567319365475d0
   fac_arr(4,3)=-0.2549077970129297123555832d0
   fac_arr(5,3)=-0.2549077970129297123555832d0
   fac_arr(6,3)=-0.2184923974396540391619285d0
   fac_arr(7,3)=-0.009103849893318918298413687/rhol**2d0
   fac_arr(8,3)=-0.02731154967995675489524106/rhol**2d0
   fac_arr(9,3)=-0.02731154967995675489524106/rhol**2d0
   fac_arr(10,3)=-0.009103849893318918298413687/rhol**2d0
   fac_arr(11,3)=0.01820769978663783659682737/rhol**2d0
   fac_arr(12,3)=0.03641539957327567319365475/rhol**2d0
   fac_arr(13,3)=0.01820769978663783659682737/rhol**2d0
   fac_arr(14,3)=0.06372694925323242808889581/rhol**2d0
   fac_arr(15,3)=0.06372694925323242808889581/rhol**2d0
   fac_arr(16,3)=0.03641539957327567319365475/rhol**2d0
else if (l.eq.4 .and. i.eq.3 .and. m.eq.2) then
   nterm_arr(1)=16
   nterm_arr(2)=20
   nterm_arr(3)=16
   lxyz_arr(1,1,1)=5 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=5 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=3 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=2
   lxyz_arr(1,5,1)=1 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=1 ; lxyz_arr(3,6,1)=4
   lxyz_arr(1,7,1)=7 ; lxyz_arr(2,7,1)=1 ; lxyz_arr(3,7,1)=0
   lxyz_arr(1,8,1)=5 ; lxyz_arr(2,8,1)=3 ; lxyz_arr(3,8,1)=0
   lxyz_arr(1,9,1)=3 ; lxyz_arr(2,9,1)=5 ; lxyz_arr(3,9,1)=0
   lxyz_arr(1,10,1)=1 ; lxyz_arr(2,10,1)=7 ; lxyz_arr(3,10,1)=0
   lxyz_arr(1,11,1)=5 ; lxyz_arr(2,11,1)=1 ; lxyz_arr(3,11,1)=2
   lxyz_arr(1,12,1)=3 ; lxyz_arr(2,12,1)=3 ; lxyz_arr(3,12,1)=2
   lxyz_arr(1,13,1)=1 ; lxyz_arr(2,13,1)=5 ; lxyz_arr(3,13,1)=2
   lxyz_arr(1,14,1)=3 ; lxyz_arr(2,14,1)=1 ; lxyz_arr(3,14,1)=4
   lxyz_arr(1,15,1)=1 ; lxyz_arr(2,15,1)=3 ; lxyz_arr(3,15,1)=4
   lxyz_arr(1,16,1)=1 ; lxyz_arr(2,16,1)=1 ; lxyz_arr(3,16,1)=6
   fac_arr(1,1)=0.05462309935991350979048212d0
   fac_arr(2,1)=0.1092461987198270195809642d0
   fac_arr(3,1)=0.05462309935991350979048212d0
   fac_arr(4,1)=-0.0728307991465513463873095d0
   fac_arr(5,1)=-0.0728307991465513463873095d0
   fac_arr(6,1)=-0.1274538985064648561777916d0
   fac_arr(7,1)=-0.009103849893318918298413687/rhol**2d0
   fac_arr(8,1)=-0.02731154967995675489524106/rhol**2d0
   fac_arr(9,1)=-0.02731154967995675489524106/rhol**2d0
   fac_arr(10,1)=-0.009103849893318918298413687/rhol**2d0
   fac_arr(11,1)=0.01820769978663783659682737/rhol**2d0
   fac_arr(12,1)=0.03641539957327567319365475/rhol**2d0
   fac_arr(13,1)=0.01820769978663783659682737/rhol**2d0
   fac_arr(14,1)=0.06372694925323242808889581/rhol**2d0
   fac_arr(15,1)=0.06372694925323242808889581/rhol**2d0
   fac_arr(16,1)=0.03641539957327567319365475/rhol**2d0
   lxyz_arr(1,1,2)=6 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=4 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=2 ; lxyz_arr(2,3,2)=4 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=6 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=4 ; lxyz_arr(2,5,2)=0 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=2 ; lxyz_arr(2,6,2)=2 ; lxyz_arr(3,6,2)=2
   lxyz_arr(1,7,2)=0 ; lxyz_arr(2,7,2)=4 ; lxyz_arr(3,7,2)=2
   lxyz_arr(1,8,2)=2 ; lxyz_arr(2,8,2)=0 ; lxyz_arr(3,8,2)=4
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=2 ; lxyz_arr(3,9,2)=4
   lxyz_arr(1,10,2)=0 ; lxyz_arr(2,10,2)=0 ; lxyz_arr(3,10,2)=6
   lxyz_arr(1,11,2)=6 ; lxyz_arr(2,11,2)=2 ; lxyz_arr(3,11,2)=0
   lxyz_arr(1,12,2)=4 ; lxyz_arr(2,12,2)=4 ; lxyz_arr(3,12,2)=0
   lxyz_arr(1,13,2)=2 ; lxyz_arr(2,13,2)=6 ; lxyz_arr(3,13,2)=0
   lxyz_arr(1,14,2)=0 ; lxyz_arr(2,14,2)=8 ; lxyz_arr(3,14,2)=0
   lxyz_arr(1,15,2)=4 ; lxyz_arr(2,15,2)=2 ; lxyz_arr(3,15,2)=2
   lxyz_arr(1,16,2)=2 ; lxyz_arr(2,16,2)=4 ; lxyz_arr(3,16,2)=2
   lxyz_arr(1,17,2)=0 ; lxyz_arr(2,17,2)=6 ; lxyz_arr(3,17,2)=2
   lxyz_arr(1,18,2)=2 ; lxyz_arr(2,18,2)=2 ; lxyz_arr(3,18,2)=4
   lxyz_arr(1,19,2)=0 ; lxyz_arr(2,19,2)=4 ; lxyz_arr(3,19,2)=4
   lxyz_arr(1,20,2)=0 ; lxyz_arr(2,20,2)=2 ; lxyz_arr(3,20,2)=6
   fac_arr(1,2)=0.009103849893318918298413687d0
   fac_arr(2,2)=0.08193464903987026468572318d0
   fac_arr(3,2)=0.1365577483997837744762053d0
   fac_arr(4,2)=0.06372694925323242808889581d0
   fac_arr(5,2)=-0.01820769978663783659682737d0
   fac_arr(6,2)=-0.1092461987198270195809642d0
   fac_arr(7,2)=-0.09103849893318918298413687d0
   fac_arr(8,2)=-0.06372694925323242808889581d0
   fac_arr(9,2)=-0.1911808477596972842666874d0
   fac_arr(10,2)=-0.03641539957327567319365475d0
   fac_arr(11,2)=-0.009103849893318918298413687/rhol**2d0
   fac_arr(12,2)=-0.02731154967995675489524106/rhol**2d0
   fac_arr(13,2)=-0.02731154967995675489524106/rhol**2d0
   fac_arr(14,2)=-0.009103849893318918298413687/rhol**2d0
   fac_arr(15,2)=0.01820769978663783659682737/rhol**2d0
   fac_arr(16,2)=0.03641539957327567319365475/rhol**2d0
   fac_arr(17,2)=0.01820769978663783659682737/rhol**2d0
   fac_arr(18,2)=0.06372694925323242808889581/rhol**2d0
   fac_arr(19,2)=0.06372694925323242808889581/rhol**2d0
   fac_arr(20,2)=0.03641539957327567319365475/rhol**2d0
   lxyz_arr(1,1,3)=4 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=5 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=3
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=3
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=1 ; lxyz_arr(3,6,3)=5
   lxyz_arr(1,7,3)=6 ; lxyz_arr(2,7,3)=1 ; lxyz_arr(3,7,3)=1
   lxyz_arr(1,8,3)=4 ; lxyz_arr(2,8,3)=3 ; lxyz_arr(3,8,3)=1
   lxyz_arr(1,9,3)=2 ; lxyz_arr(2,9,3)=5 ; lxyz_arr(3,9,3)=1
   lxyz_arr(1,10,3)=0 ; lxyz_arr(2,10,3)=7 ; lxyz_arr(3,10,3)=1
   lxyz_arr(1,11,3)=4 ; lxyz_arr(2,11,3)=1 ; lxyz_arr(3,11,3)=3
   lxyz_arr(1,12,3)=2 ; lxyz_arr(2,12,3)=3 ; lxyz_arr(3,12,3)=3
   lxyz_arr(1,13,3)=0 ; lxyz_arr(2,13,3)=5 ; lxyz_arr(3,13,3)=3
   lxyz_arr(1,14,3)=2 ; lxyz_arr(2,14,3)=1 ; lxyz_arr(3,14,3)=5
   lxyz_arr(1,15,3)=0 ; lxyz_arr(2,15,3)=3 ; lxyz_arr(3,15,3)=5
   lxyz_arr(1,16,3)=0 ; lxyz_arr(2,16,3)=1 ; lxyz_arr(3,16,3)=7
   fac_arr(1,3)=-0.03641539957327567319365475d0
   fac_arr(2,3)=-0.0728307991465513463873095d0
   fac_arr(3,3)=-0.03641539957327567319365475d0
   fac_arr(4,3)=-0.2549077970129297123555832d0
   fac_arr(5,3)=-0.2549077970129297123555832d0
   fac_arr(6,3)=-0.2184923974396540391619285d0
   fac_arr(7,3)=-0.009103849893318918298413687/rhol**2d0
   fac_arr(8,3)=-0.02731154967995675489524106/rhol**2d0
   fac_arr(9,3)=-0.02731154967995675489524106/rhol**2d0
   fac_arr(10,3)=-0.009103849893318918298413687/rhol**2d0
   fac_arr(11,3)=0.01820769978663783659682737/rhol**2d0
   fac_arr(12,3)=0.03641539957327567319365475/rhol**2d0
   fac_arr(13,3)=0.01820769978663783659682737/rhol**2d0
   fac_arr(14,3)=0.06372694925323242808889581/rhol**2d0
   fac_arr(15,3)=0.06372694925323242808889581/rhol**2d0
   fac_arr(16,3)=0.03641539957327567319365475/rhol**2d0
else if (l.eq.4 .and. i.eq.3 .and. m.eq.3) then
   nterm_arr(1)=16
   nterm_arr(2)=16
   nterm_arr(3)=20
   lxyz_arr(1,1,1)=5 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=1
   lxyz_arr(1,4,1)=3 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=3
   lxyz_arr(1,5,1)=1 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=3
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=5
   lxyz_arr(1,7,1)=7 ; lxyz_arr(2,7,1)=0 ; lxyz_arr(3,7,1)=1
   lxyz_arr(1,8,1)=5 ; lxyz_arr(2,8,1)=2 ; lxyz_arr(3,8,1)=1
   lxyz_arr(1,9,1)=3 ; lxyz_arr(2,9,1)=4 ; lxyz_arr(3,9,1)=1
   lxyz_arr(1,10,1)=1 ; lxyz_arr(2,10,1)=6 ; lxyz_arr(3,10,1)=1
   lxyz_arr(1,11,1)=5 ; lxyz_arr(2,11,1)=0 ; lxyz_arr(3,11,1)=3
   lxyz_arr(1,12,1)=3 ; lxyz_arr(2,12,1)=2 ; lxyz_arr(3,12,1)=3
   lxyz_arr(1,13,1)=1 ; lxyz_arr(2,13,1)=4 ; lxyz_arr(3,13,1)=3
   lxyz_arr(1,14,1)=3 ; lxyz_arr(2,14,1)=0 ; lxyz_arr(3,14,1)=5
   lxyz_arr(1,15,1)=1 ; lxyz_arr(2,15,1)=2 ; lxyz_arr(3,15,1)=5
   lxyz_arr(1,16,1)=1 ; lxyz_arr(2,16,1)=0 ; lxyz_arr(3,16,1)=7
   fac_arr(1,1)=0.1337987216011345233133409d0
   fac_arr(2,1)=0.2675974432022690466266818d0
   fac_arr(3,1)=0.1337987216011345233133409d0
   fac_arr(4,1)=0.1189321969787862429451919d0
   fac_arr(5,1)=0.1189321969787862429451919d0
   fac_arr(6,1)=-0.01486652462234828036814899d0
   fac_arr(7,1)=-0.02229978693352242055222348/rhol**2d0
   fac_arr(8,1)=-0.06689936080056726165667044/rhol**2d0
   fac_arr(9,1)=-0.06689936080056726165667044/rhol**2d0
   fac_arr(10,1)=-0.02229978693352242055222348/rhol**2d0
   fac_arr(11,1)=-0.02973304924469656073629797/rhol**2d0
   fac_arr(12,1)=-0.05946609848939312147259594/rhol**2d0
   fac_arr(13,1)=-0.02973304924469656073629797/rhol**2d0
   fac_arr(14,1)=0.007433262311174140184074493/rhol**2d0
   fac_arr(15,1)=0.007433262311174140184074493/rhol**2d0
   fac_arr(16,1)=0.01486652462234828036814899/rhol**2d0
   lxyz_arr(1,1,2)=4 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=5 ; lxyz_arr(3,3,2)=1
   lxyz_arr(1,4,2)=2 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=3
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=3
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=1 ; lxyz_arr(3,6,2)=5
   lxyz_arr(1,7,2)=6 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=1
   lxyz_arr(1,8,2)=4 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=1
   lxyz_arr(1,9,2)=2 ; lxyz_arr(2,9,2)=5 ; lxyz_arr(3,9,2)=1
   lxyz_arr(1,10,2)=0 ; lxyz_arr(2,10,2)=7 ; lxyz_arr(3,10,2)=1
   lxyz_arr(1,11,2)=4 ; lxyz_arr(2,11,2)=1 ; lxyz_arr(3,11,2)=3
   lxyz_arr(1,12,2)=2 ; lxyz_arr(2,12,2)=3 ; lxyz_arr(3,12,2)=3
   lxyz_arr(1,13,2)=0 ; lxyz_arr(2,13,2)=5 ; lxyz_arr(3,13,2)=3
   lxyz_arr(1,14,2)=2 ; lxyz_arr(2,14,2)=1 ; lxyz_arr(3,14,2)=5
   lxyz_arr(1,15,2)=0 ; lxyz_arr(2,15,2)=3 ; lxyz_arr(3,15,2)=5
   lxyz_arr(1,16,2)=0 ; lxyz_arr(2,16,2)=1 ; lxyz_arr(3,16,2)=7
   fac_arr(1,2)=0.1337987216011345233133409d0
   fac_arr(2,2)=0.2675974432022690466266818d0
   fac_arr(3,2)=0.1337987216011345233133409d0
   fac_arr(4,2)=0.1189321969787862429451919d0
   fac_arr(5,2)=0.1189321969787862429451919d0
   fac_arr(6,2)=-0.01486652462234828036814899d0
   fac_arr(7,2)=-0.02229978693352242055222348/rhol**2d0
   fac_arr(8,2)=-0.06689936080056726165667044/rhol**2d0
   fac_arr(9,2)=-0.06689936080056726165667044/rhol**2d0
   fac_arr(10,2)=-0.02229978693352242055222348/rhol**2d0
   fac_arr(11,2)=-0.02973304924469656073629797/rhol**2d0
   fac_arr(12,2)=-0.05946609848939312147259594/rhol**2d0
   fac_arr(13,2)=-0.02973304924469656073629797/rhol**2d0
   fac_arr(14,2)=0.007433262311174140184074493/rhol**2d0
   fac_arr(15,2)=0.007433262311174140184074493/rhol**2d0
   fac_arr(16,2)=0.01486652462234828036814899/rhol**2d0
   lxyz_arr(1,1,3)=6 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=4 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=2 ; lxyz_arr(2,3,3)=4 ; lxyz_arr(3,3,3)=0
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=6 ; lxyz_arr(3,4,3)=0
   lxyz_arr(1,5,3)=4 ; lxyz_arr(2,5,3)=0 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=2 ; lxyz_arr(2,6,3)=2 ; lxyz_arr(3,6,3)=2
   lxyz_arr(1,7,3)=0 ; lxyz_arr(2,7,3)=4 ; lxyz_arr(3,7,3)=2
   lxyz_arr(1,8,3)=2 ; lxyz_arr(2,8,3)=0 ; lxyz_arr(3,8,3)=4
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=2 ; lxyz_arr(3,9,3)=4
   lxyz_arr(1,10,3)=0 ; lxyz_arr(2,10,3)=0 ; lxyz_arr(3,10,3)=6
   lxyz_arr(1,11,3)=6 ; lxyz_arr(2,11,3)=0 ; lxyz_arr(3,11,3)=2
   lxyz_arr(1,12,3)=4 ; lxyz_arr(2,12,3)=2 ; lxyz_arr(3,12,3)=2
   lxyz_arr(1,13,3)=2 ; lxyz_arr(2,13,3)=4 ; lxyz_arr(3,13,3)=2
   lxyz_arr(1,14,3)=0 ; lxyz_arr(2,14,3)=6 ; lxyz_arr(3,14,3)=2
   lxyz_arr(1,15,3)=4 ; lxyz_arr(2,15,3)=0 ; lxyz_arr(3,15,3)=4
   lxyz_arr(1,16,3)=2 ; lxyz_arr(2,16,3)=2 ; lxyz_arr(3,16,3)=4
   lxyz_arr(1,17,3)=0 ; lxyz_arr(2,17,3)=4 ; lxyz_arr(3,17,3)=4
   lxyz_arr(1,18,3)=2 ; lxyz_arr(2,18,3)=0 ; lxyz_arr(3,18,3)=6
   lxyz_arr(1,19,3)=0 ; lxyz_arr(2,19,3)=2 ; lxyz_arr(3,19,3)=6
   lxyz_arr(1,20,3)=0 ; lxyz_arr(2,20,3)=0 ; lxyz_arr(3,20,3)=8
   fac_arr(1,3)=0.02229978693352242055222348d0
   fac_arr(2,3)=0.06689936080056726165667044d0
   fac_arr(3,3)=0.06689936080056726165667044d0
   fac_arr(4,3)=0.02229978693352242055222348d0
   fac_arr(5,3)=0.08919914773408968220889392d0
   fac_arr(6,3)=0.1783982954681793644177878d0
   fac_arr(7,3)=0.08919914773408968220889392d0
   fac_arr(8,3)=-0.03716631155587070092037247d0
   fac_arr(9,3)=-0.03716631155587070092037247d0
   fac_arr(10,3)=-0.1040656723564379625770429d0
   fac_arr(11,3)=-0.02229978693352242055222348/rhol**2d0
   fac_arr(12,3)=-0.06689936080056726165667044/rhol**2d0
   fac_arr(13,3)=-0.06689936080056726165667044/rhol**2d0
   fac_arr(14,3)=-0.02229978693352242055222348/rhol**2d0
   fac_arr(15,3)=-0.02973304924469656073629797/rhol**2d0
   fac_arr(16,3)=-0.05946609848939312147259594/rhol**2d0
   fac_arr(17,3)=-0.02973304924469656073629797/rhol**2d0
   fac_arr(18,3)=0.007433262311174140184074493/rhol**2d0
   fac_arr(19,3)=0.007433262311174140184074493/rhol**2d0
   fac_arr(20,3)=0.01486652462234828036814899/rhol**2d0
else if (l.eq.4 .and. i.eq.3 .and. m.eq.4) then
   nterm_arr(1)=18
   nterm_arr(2)=15
   nterm_arr(3)=14
   lxyz_arr(1,1,1)=6 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=4 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=2 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=0 ; lxyz_arr(2,4,1)=6 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=4 ; lxyz_arr(2,5,1)=0 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=2 ; lxyz_arr(2,6,1)=2 ; lxyz_arr(3,6,1)=2
   lxyz_arr(1,7,1)=0 ; lxyz_arr(2,7,1)=4 ; lxyz_arr(3,7,1)=2
   lxyz_arr(1,8,1)=2 ; lxyz_arr(2,8,1)=0 ; lxyz_arr(3,8,1)=4
   lxyz_arr(1,9,1)=0 ; lxyz_arr(2,9,1)=2 ; lxyz_arr(3,9,1)=4
   lxyz_arr(1,10,1)=8 ; lxyz_arr(2,10,1)=0 ; lxyz_arr(3,10,1)=0
   lxyz_arr(1,11,1)=6 ; lxyz_arr(2,11,1)=2 ; lxyz_arr(3,11,1)=0
   lxyz_arr(1,12,1)=4 ; lxyz_arr(2,12,1)=4 ; lxyz_arr(3,12,1)=0
   lxyz_arr(1,13,1)=2 ; lxyz_arr(2,13,1)=6 ; lxyz_arr(3,13,1)=0
   lxyz_arr(1,14,1)=6 ; lxyz_arr(2,14,1)=0 ; lxyz_arr(3,14,1)=2
   lxyz_arr(1,15,1)=4 ; lxyz_arr(2,15,1)=2 ; lxyz_arr(3,15,1)=2
   lxyz_arr(1,16,1)=2 ; lxyz_arr(2,16,1)=4 ; lxyz_arr(3,16,1)=2
   lxyz_arr(1,17,1)=4 ; lxyz_arr(2,17,1)=0 ; lxyz_arr(3,17,1)=4
   lxyz_arr(1,18,1)=2 ; lxyz_arr(2,18,1)=2 ; lxyz_arr(3,18,1)=4
   fac_arr(1,1)=0.08227113772079145865717289d0
   fac_arr(2,1)=-0.05876509837199389904083778d0
   fac_arr(3,1)=-0.1762952951159816971225133d0
   fac_arr(4,1)=-0.03525905902319633942450267d0
   fac_arr(5,1)=0.1175301967439877980816756d0
   fac_arr(6,1)=-0.1410362360927853576980107d0
   fac_arr(7,1)=-0.07051811804639267884900533d0
   fac_arr(8,1)=0.03525905902319633942450267d0
   fac_arr(9,1)=-0.03525905902319633942450267d0
   fac_arr(10,1)=-0.01175301967439877980816756/rhol**2d0
   fac_arr(11,1)=0.01175301967439877980816756/rhol**2d0
   fac_arr(12,1)=0.05876509837199389904083778/rhol**2d0
   fac_arr(13,1)=0.03525905902319633942450267/rhol**2d0
   fac_arr(14,1)=-0.02350603934879755961633511/rhol**2d0
   fac_arr(15,1)=0.04701207869759511923267022/rhol**2d0
   fac_arr(16,1)=0.07051811804639267884900533/rhol**2d0
   fac_arr(17,1)=-0.01175301967439877980816756/rhol**2d0
   fac_arr(18,1)=0.03525905902319633942450267/rhol**2d0
   lxyz_arr(1,1,2)=5 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=3 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=5 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=3 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=2
   lxyz_arr(1,5,2)=1 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=1 ; lxyz_arr(3,6,2)=4
   lxyz_arr(1,7,2)=7 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=0
   lxyz_arr(1,8,2)=5 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=0
   lxyz_arr(1,9,2)=3 ; lxyz_arr(2,9,2)=5 ; lxyz_arr(3,9,2)=0
   lxyz_arr(1,10,2)=1 ; lxyz_arr(2,10,2)=7 ; lxyz_arr(3,10,2)=0
   lxyz_arr(1,11,2)=5 ; lxyz_arr(2,11,2)=1 ; lxyz_arr(3,11,2)=2
   lxyz_arr(1,12,2)=3 ; lxyz_arr(2,12,2)=3 ; lxyz_arr(3,12,2)=2
   lxyz_arr(1,13,2)=1 ; lxyz_arr(2,13,2)=5 ; lxyz_arr(3,13,2)=2
   lxyz_arr(1,14,2)=3 ; lxyz_arr(2,14,2)=1 ; lxyz_arr(3,14,2)=4
   lxyz_arr(1,15,2)=1 ; lxyz_arr(2,15,2)=3 ; lxyz_arr(3,15,2)=4
   fac_arr(1,2)=-0.02350603934879755961633511d0
   fac_arr(2,2)=-0.2350603934879755961633511d0
   fac_arr(3,2)=-0.211554354139178036547016d0
   fac_arr(4,2)=-0.09402415739519023846534044d0
   fac_arr(5,2)=-0.2820724721855707153960213d0
   fac_arr(6,2)=-0.07051811804639267884900533d0
   fac_arr(7,2)=-0.01175301967439877980816756/rhol**2d0
   fac_arr(8,2)=0.01175301967439877980816756/rhol**2d0
   fac_arr(9,2)=0.05876509837199389904083778/rhol**2d0
   fac_arr(10,2)=0.03525905902319633942450267/rhol**2d0
   fac_arr(11,2)=-0.02350603934879755961633511/rhol**2d0
   fac_arr(12,2)=0.04701207869759511923267022/rhol**2d0
   fac_arr(13,2)=0.07051811804639267884900533/rhol**2d0
   fac_arr(14,2)=-0.01175301967439877980816756/rhol**2d0
   fac_arr(15,2)=0.03525905902319633942450267/rhol**2d0
   lxyz_arr(1,1,3)=5 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=3 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=4 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=3 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=3
   lxyz_arr(1,5,3)=1 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=3
   lxyz_arr(1,6,3)=7 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=5 ; lxyz_arr(2,7,3)=2 ; lxyz_arr(3,7,3)=1
   lxyz_arr(1,8,3)=3 ; lxyz_arr(2,8,3)=4 ; lxyz_arr(3,8,3)=1
   lxyz_arr(1,9,3)=1 ; lxyz_arr(2,9,3)=6 ; lxyz_arr(3,9,3)=1
   lxyz_arr(1,10,3)=5 ; lxyz_arr(2,10,3)=0 ; lxyz_arr(3,10,3)=3
   lxyz_arr(1,11,3)=3 ; lxyz_arr(2,11,3)=2 ; lxyz_arr(3,11,3)=3
   lxyz_arr(1,12,3)=1 ; lxyz_arr(2,12,3)=4 ; lxyz_arr(3,12,3)=3
   lxyz_arr(1,13,3)=3 ; lxyz_arr(2,13,3)=0 ; lxyz_arr(3,13,3)=5
   lxyz_arr(1,14,3)=1 ; lxyz_arr(2,14,3)=2 ; lxyz_arr(3,14,3)=5
   fac_arr(1,3)=0.04701207869759511923267022d0
   fac_arr(2,3)=-0.09402415739519023846534044d0
   fac_arr(3,3)=-0.1410362360927853576980107d0
   fac_arr(4,3)=0.04701207869759511923267022d0
   fac_arr(5,3)=-0.1410362360927853576980107d0
   fac_arr(6,3)=-0.01175301967439877980816756/rhol**2d0
   fac_arr(7,3)=0.01175301967439877980816756/rhol**2d0
   fac_arr(8,3)=0.05876509837199389904083778/rhol**2d0
   fac_arr(9,3)=0.03525905902319633942450267/rhol**2d0
   fac_arr(10,3)=-0.02350603934879755961633511/rhol**2d0
   fac_arr(11,3)=0.04701207869759511923267022/rhol**2d0
   fac_arr(12,3)=0.07051811804639267884900533/rhol**2d0
   fac_arr(13,3)=-0.01175301967439877980816756/rhol**2d0
   fac_arr(14,3)=0.03525905902319633942450267/rhol**2d0
else if (l.eq.4 .and. i.eq.3 .and. m.eq.5) then
   nterm_arr(1)=15
   nterm_arr(2)=18
   nterm_arr(3)=14
   lxyz_arr(1,1,1)=5 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=5 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=3 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=2
   lxyz_arr(1,5,1)=1 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=1 ; lxyz_arr(3,6,1)=4
   lxyz_arr(1,7,1)=7 ; lxyz_arr(2,7,1)=1 ; lxyz_arr(3,7,1)=0
   lxyz_arr(1,8,1)=5 ; lxyz_arr(2,8,1)=3 ; lxyz_arr(3,8,1)=0
   lxyz_arr(1,9,1)=3 ; lxyz_arr(2,9,1)=5 ; lxyz_arr(3,9,1)=0
   lxyz_arr(1,10,1)=1 ; lxyz_arr(2,10,1)=7 ; lxyz_arr(3,10,1)=0
   lxyz_arr(1,11,1)=5 ; lxyz_arr(2,11,1)=1 ; lxyz_arr(3,11,1)=2
   lxyz_arr(1,12,1)=3 ; lxyz_arr(2,12,1)=3 ; lxyz_arr(3,12,1)=2
   lxyz_arr(1,13,1)=1 ; lxyz_arr(2,13,1)=5 ; lxyz_arr(3,13,1)=2
   lxyz_arr(1,14,1)=3 ; lxyz_arr(2,14,1)=1 ; lxyz_arr(3,14,1)=4
   lxyz_arr(1,15,1)=1 ; lxyz_arr(2,15,1)=3 ; lxyz_arr(3,15,1)=4
   fac_arr(1,1)=-0.211554354139178036547016d0
   fac_arr(2,1)=-0.2350603934879755961633511d0
   fac_arr(3,1)=-0.02350603934879755961633511d0
   fac_arr(4,1)=-0.2820724721855707153960213d0
   fac_arr(5,1)=-0.09402415739519023846534044d0
   fac_arr(6,1)=-0.07051811804639267884900533d0
   fac_arr(7,1)=0.03525905902319633942450267/rhol**2d0
   fac_arr(8,1)=0.05876509837199389904083778/rhol**2d0
   fac_arr(9,1)=0.01175301967439877980816756/rhol**2d0
   fac_arr(10,1)=-0.01175301967439877980816756/rhol**2d0
   fac_arr(11,1)=0.07051811804639267884900533/rhol**2d0
   fac_arr(12,1)=0.04701207869759511923267022/rhol**2d0
   fac_arr(13,1)=-0.02350603934879755961633511/rhol**2d0
   fac_arr(14,1)=0.03525905902319633942450267/rhol**2d0
   fac_arr(15,1)=-0.01175301967439877980816756/rhol**2d0
   lxyz_arr(1,1,2)=6 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=4 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=2 ; lxyz_arr(2,3,2)=4 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=6 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=4 ; lxyz_arr(2,5,2)=0 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=2 ; lxyz_arr(2,6,2)=2 ; lxyz_arr(3,6,2)=2
   lxyz_arr(1,7,2)=0 ; lxyz_arr(2,7,2)=4 ; lxyz_arr(3,7,2)=2
   lxyz_arr(1,8,2)=2 ; lxyz_arr(2,8,2)=0 ; lxyz_arr(3,8,2)=4
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=2 ; lxyz_arr(3,9,2)=4
   lxyz_arr(1,10,2)=6 ; lxyz_arr(2,10,2)=2 ; lxyz_arr(3,10,2)=0
   lxyz_arr(1,11,2)=4 ; lxyz_arr(2,11,2)=4 ; lxyz_arr(3,11,2)=0
   lxyz_arr(1,12,2)=2 ; lxyz_arr(2,12,2)=6 ; lxyz_arr(3,12,2)=0
   lxyz_arr(1,13,2)=0 ; lxyz_arr(2,13,2)=8 ; lxyz_arr(3,13,2)=0
   lxyz_arr(1,14,2)=4 ; lxyz_arr(2,14,2)=2 ; lxyz_arr(3,14,2)=2
   lxyz_arr(1,15,2)=2 ; lxyz_arr(2,15,2)=4 ; lxyz_arr(3,15,2)=2
   lxyz_arr(1,16,2)=0 ; lxyz_arr(2,16,2)=6 ; lxyz_arr(3,16,2)=2
   lxyz_arr(1,17,2)=2 ; lxyz_arr(2,17,2)=2 ; lxyz_arr(3,17,2)=4
   lxyz_arr(1,18,2)=0 ; lxyz_arr(2,18,2)=4 ; lxyz_arr(3,18,2)=4
   fac_arr(1,2)=-0.03525905902319633942450267d0
   fac_arr(2,2)=-0.1762952951159816971225133d0
   fac_arr(3,2)=-0.05876509837199389904083778d0
   fac_arr(4,2)=0.08227113772079145865717289d0
   fac_arr(5,2)=-0.07051811804639267884900533d0
   fac_arr(6,2)=-0.1410362360927853576980107d0
   fac_arr(7,2)=0.1175301967439877980816756d0
   fac_arr(8,2)=-0.03525905902319633942450267d0
   fac_arr(9,2)=0.03525905902319633942450267d0
   fac_arr(10,2)=0.03525905902319633942450267/rhol**2d0
   fac_arr(11,2)=0.05876509837199389904083778/rhol**2d0
   fac_arr(12,2)=0.01175301967439877980816756/rhol**2d0
   fac_arr(13,2)=-0.01175301967439877980816756/rhol**2d0
   fac_arr(14,2)=0.07051811804639267884900533/rhol**2d0
   fac_arr(15,2)=0.04701207869759511923267022/rhol**2d0
   fac_arr(16,2)=-0.02350603934879755961633511/rhol**2d0
   fac_arr(17,2)=0.03525905902319633942450267/rhol**2d0
   fac_arr(18,2)=-0.01175301967439877980816756/rhol**2d0
   lxyz_arr(1,1,3)=4 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=5 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=3
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=3
   lxyz_arr(1,6,3)=6 ; lxyz_arr(2,6,3)=1 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=4 ; lxyz_arr(2,7,3)=3 ; lxyz_arr(3,7,3)=1
   lxyz_arr(1,8,3)=2 ; lxyz_arr(2,8,3)=5 ; lxyz_arr(3,8,3)=1
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=7 ; lxyz_arr(3,9,3)=1
   lxyz_arr(1,10,3)=4 ; lxyz_arr(2,10,3)=1 ; lxyz_arr(3,10,3)=3
   lxyz_arr(1,11,3)=2 ; lxyz_arr(2,11,3)=3 ; lxyz_arr(3,11,3)=3
   lxyz_arr(1,12,3)=0 ; lxyz_arr(2,12,3)=5 ; lxyz_arr(3,12,3)=3
   lxyz_arr(1,13,3)=2 ; lxyz_arr(2,13,3)=1 ; lxyz_arr(3,13,3)=5
   lxyz_arr(1,14,3)=0 ; lxyz_arr(2,14,3)=3 ; lxyz_arr(3,14,3)=5
   fac_arr(1,3)=-0.1410362360927853576980107d0
   fac_arr(2,3)=-0.09402415739519023846534044d0
   fac_arr(3,3)=0.04701207869759511923267022d0
   fac_arr(4,3)=-0.1410362360927853576980107d0
   fac_arr(5,3)=0.04701207869759511923267022d0
   fac_arr(6,3)=0.03525905902319633942450267/rhol**2d0
   fac_arr(7,3)=0.05876509837199389904083778/rhol**2d0
   fac_arr(8,3)=0.01175301967439877980816756/rhol**2d0
   fac_arr(9,3)=-0.01175301967439877980816756/rhol**2d0
   fac_arr(10,3)=0.07051811804639267884900533/rhol**2d0
   fac_arr(11,3)=0.04701207869759511923267022/rhol**2d0
   fac_arr(12,3)=-0.02350603934879755961633511/rhol**2d0
   fac_arr(13,3)=0.03525905902319633942450267/rhol**2d0
   fac_arr(14,3)=-0.01175301967439877980816756/rhol**2d0
else if (l.eq.4 .and. i.eq.3 .and. m.eq.6) then
   nterm_arr(1)=13
   nterm_arr(2)=13
   nterm_arr(3)=16
   lxyz_arr(1,1,1)=5 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=1
   lxyz_arr(1,4,1)=3 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=3
   lxyz_arr(1,5,1)=1 ; lxyz_arr(2,5,1)=0 ; lxyz_arr(3,5,1)=5
   lxyz_arr(1,6,1)=7 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=1
   lxyz_arr(1,7,1)=5 ; lxyz_arr(2,7,1)=2 ; lxyz_arr(3,7,1)=1
   lxyz_arr(1,8,1)=3 ; lxyz_arr(2,8,1)=4 ; lxyz_arr(3,8,1)=1
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=6 ; lxyz_arr(3,9,1)=1
   lxyz_arr(1,10,1)=5 ; lxyz_arr(2,10,1)=0 ; lxyz_arr(3,10,1)=3
   lxyz_arr(1,11,1)=1 ; lxyz_arr(2,11,1)=4 ; lxyz_arr(3,11,1)=3
   lxyz_arr(1,12,1)=3 ; lxyz_arr(2,12,1)=0 ; lxyz_arr(3,12,1)=5
   lxyz_arr(1,13,1)=1 ; lxyz_arr(2,13,1)=2 ; lxyz_arr(3,13,1)=5
   fac_arr(1,1)=0.1727334068350121925245643d0
   fac_arr(2,1)=0.1151556045566747950163762d0
   fac_arr(3,1)=-0.05757780227833739750818811d0
   fac_arr(4,1)=0.2303112091133495900327524d0
   fac_arr(5,1)=0.05757780227833739750818811d0
   fac_arr(6,1)=-0.02878890113916869875409405/rhol**2d0
   fac_arr(7,1)=-0.02878890113916869875409405/rhol**2d0
   fac_arr(8,1)=0.02878890113916869875409405/rhol**2d0
   fac_arr(9,1)=0.02878890113916869875409405/rhol**2d0
   fac_arr(10,1)=-0.05757780227833739750818811/rhol**2d0
   fac_arr(11,1)=0.05757780227833739750818811/rhol**2d0
   fac_arr(12,1)=-0.02878890113916869875409405/rhol**2d0
   fac_arr(13,1)=0.02878890113916869875409405/rhol**2d0
   lxyz_arr(1,1,2)=4 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=5 ; lxyz_arr(3,3,2)=1
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=3 ; lxyz_arr(3,4,2)=3
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=1 ; lxyz_arr(3,5,2)=5
   lxyz_arr(1,6,2)=6 ; lxyz_arr(2,6,2)=1 ; lxyz_arr(3,6,2)=1
   lxyz_arr(1,7,2)=4 ; lxyz_arr(2,7,2)=3 ; lxyz_arr(3,7,2)=1
   lxyz_arr(1,8,2)=2 ; lxyz_arr(2,8,2)=5 ; lxyz_arr(3,8,2)=1
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=7 ; lxyz_arr(3,9,2)=1
   lxyz_arr(1,10,2)=4 ; lxyz_arr(2,10,2)=1 ; lxyz_arr(3,10,2)=3
   lxyz_arr(1,11,2)=0 ; lxyz_arr(2,11,2)=5 ; lxyz_arr(3,11,2)=3
   lxyz_arr(1,12,2)=2 ; lxyz_arr(2,12,2)=1 ; lxyz_arr(3,12,2)=5
   lxyz_arr(1,13,2)=0 ; lxyz_arr(2,13,2)=3 ; lxyz_arr(3,13,2)=5
   fac_arr(1,2)=0.05757780227833739750818811d0
   fac_arr(2,2)=-0.1151556045566747950163762d0
   fac_arr(3,2)=-0.1727334068350121925245643d0
   fac_arr(4,2)=-0.2303112091133495900327524d0
   fac_arr(5,2)=-0.05757780227833739750818811d0
   fac_arr(6,2)=-0.02878890113916869875409405/rhol**2d0
   fac_arr(7,2)=-0.02878890113916869875409405/rhol**2d0
   fac_arr(8,2)=0.02878890113916869875409405/rhol**2d0
   fac_arr(9,2)=0.02878890113916869875409405/rhol**2d0
   fac_arr(10,2)=-0.05757780227833739750818811/rhol**2d0
   fac_arr(11,2)=0.05757780227833739750818811/rhol**2d0
   fac_arr(12,2)=-0.02878890113916869875409405/rhol**2d0
   fac_arr(13,2)=0.02878890113916869875409405/rhol**2d0
   lxyz_arr(1,1,3)=6 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=4 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=2 ; lxyz_arr(2,3,3)=4 ; lxyz_arr(3,3,3)=0
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=6 ; lxyz_arr(3,4,3)=0
   lxyz_arr(1,5,3)=4 ; lxyz_arr(2,5,3)=0 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=4 ; lxyz_arr(3,6,3)=2
   lxyz_arr(1,7,3)=2 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=4
   lxyz_arr(1,8,3)=0 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=4
   lxyz_arr(1,9,3)=6 ; lxyz_arr(2,9,3)=0 ; lxyz_arr(3,9,3)=2
   lxyz_arr(1,10,3)=4 ; lxyz_arr(2,10,3)=2 ; lxyz_arr(3,10,3)=2
   lxyz_arr(1,11,3)=2 ; lxyz_arr(2,11,3)=4 ; lxyz_arr(3,11,3)=2
   lxyz_arr(1,12,3)=0 ; lxyz_arr(2,12,3)=6 ; lxyz_arr(3,12,3)=2
   lxyz_arr(1,13,3)=4 ; lxyz_arr(2,13,3)=0 ; lxyz_arr(3,13,3)=4
   lxyz_arr(1,14,3)=0 ; lxyz_arr(2,14,3)=4 ; lxyz_arr(3,14,3)=4
   lxyz_arr(1,15,3)=2 ; lxyz_arr(2,15,3)=0 ; lxyz_arr(3,15,3)=6
   lxyz_arr(1,16,3)=0 ; lxyz_arr(2,16,3)=2 ; lxyz_arr(3,16,3)=6
   fac_arr(1,3)=0.02878890113916869875409405d0
   fac_arr(2,3)=0.02878890113916869875409405d0
   fac_arr(3,3)=-0.02878890113916869875409405d0
   fac_arr(4,3)=-0.02878890113916869875409405d0
   fac_arr(5,3)=0.1727334068350121925245643d0
   fac_arr(6,3)=-0.1727334068350121925245643d0
   fac_arr(7,3)=0.1439445056958434937704703d0
   fac_arr(8,3)=-0.1439445056958434937704703d0
   fac_arr(9,3)=-0.02878890113916869875409405/rhol**2d0
   fac_arr(10,3)=-0.02878890113916869875409405/rhol**2d0
   fac_arr(11,3)=0.02878890113916869875409405/rhol**2d0
   fac_arr(12,3)=0.02878890113916869875409405/rhol**2d0
   fac_arr(13,3)=-0.05757780227833739750818811/rhol**2d0
   fac_arr(14,3)=0.05757780227833739750818811/rhol**2d0
   fac_arr(15,3)=-0.02878890113916869875409405/rhol**2d0
   fac_arr(16,3)=0.02878890113916869875409405/rhol**2d0
else if (l.eq.4 .and. i.eq.3 .and. m.eq.7) then
   nterm_arr(1)=12
   nterm_arr(2)=12
   nterm_arr(3)=12
   lxyz_arr(1,1,1)=4 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=5 ; lxyz_arr(3,3,1)=1
   lxyz_arr(1,4,1)=2 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=3
   lxyz_arr(1,5,1)=0 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=3
   lxyz_arr(1,6,1)=0 ; lxyz_arr(2,6,1)=1 ; lxyz_arr(3,6,1)=5
   lxyz_arr(1,7,1)=6 ; lxyz_arr(2,7,1)=1 ; lxyz_arr(3,7,1)=1
   lxyz_arr(1,8,1)=4 ; lxyz_arr(2,8,1)=3 ; lxyz_arr(3,8,1)=1
   lxyz_arr(1,9,1)=2 ; lxyz_arr(2,9,1)=5 ; lxyz_arr(3,9,1)=1
   lxyz_arr(1,10,1)=4 ; lxyz_arr(2,10,1)=1 ; lxyz_arr(3,10,1)=3
   lxyz_arr(1,11,1)=2 ; lxyz_arr(2,11,1)=3 ; lxyz_arr(3,11,1)=3
   lxyz_arr(1,12,1)=2 ; lxyz_arr(2,12,1)=1 ; lxyz_arr(3,12,1)=5
   fac_arr(1,1)=0.2878890113916869875409405d0
   fac_arr(2,1)=0.3454668136700243850491286d0
   fac_arr(3,1)=0.05757780227833739750818811d0
   fac_arr(4,1)=0.3454668136700243850491286d0
   fac_arr(5,1)=0.1151556045566747950163762d0
   fac_arr(6,1)=0.05757780227833739750818811d0
   fac_arr(7,1)=-0.05757780227833739750818811/rhol**2d0
   fac_arr(8,1)=-0.1151556045566747950163762/rhol**2d0
   fac_arr(9,1)=-0.05757780227833739750818811/rhol**2d0
   fac_arr(10,1)=-0.1151556045566747950163762/rhol**2d0
   fac_arr(11,1)=-0.1151556045566747950163762/rhol**2d0
   fac_arr(12,1)=-0.05757780227833739750818811/rhol**2d0
   lxyz_arr(1,1,2)=5 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=3 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=4 ; lxyz_arr(3,3,2)=1
   lxyz_arr(1,4,2)=3 ; lxyz_arr(2,4,2)=0 ; lxyz_arr(3,4,2)=3
   lxyz_arr(1,5,2)=1 ; lxyz_arr(2,5,2)=2 ; lxyz_arr(3,5,2)=3
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=0 ; lxyz_arr(3,6,2)=5
   lxyz_arr(1,7,2)=5 ; lxyz_arr(2,7,2)=2 ; lxyz_arr(3,7,2)=1
   lxyz_arr(1,8,2)=3 ; lxyz_arr(2,8,2)=4 ; lxyz_arr(3,8,2)=1
   lxyz_arr(1,9,2)=1 ; lxyz_arr(2,9,2)=6 ; lxyz_arr(3,9,2)=1
   lxyz_arr(1,10,2)=3 ; lxyz_arr(2,10,2)=2 ; lxyz_arr(3,10,2)=3
   lxyz_arr(1,11,2)=1 ; lxyz_arr(2,11,2)=4 ; lxyz_arr(3,11,2)=3
   lxyz_arr(1,12,2)=1 ; lxyz_arr(2,12,2)=2 ; lxyz_arr(3,12,2)=5
   fac_arr(1,2)=0.05757780227833739750818811d0
   fac_arr(2,2)=0.3454668136700243850491286d0
   fac_arr(3,2)=0.2878890113916869875409405d0
   fac_arr(4,2)=0.1151556045566747950163762d0
   fac_arr(5,2)=0.3454668136700243850491286d0
   fac_arr(6,2)=0.05757780227833739750818811d0
   fac_arr(7,2)=-0.05757780227833739750818811/rhol**2d0
   fac_arr(8,2)=-0.1151556045566747950163762/rhol**2d0
   fac_arr(9,2)=-0.05757780227833739750818811/rhol**2d0
   fac_arr(10,2)=-0.1151556045566747950163762/rhol**2d0
   fac_arr(11,2)=-0.1151556045566747950163762/rhol**2d0
   fac_arr(12,2)=-0.05757780227833739750818811/rhol**2d0
   lxyz_arr(1,1,3)=5 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=3 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=5 ; lxyz_arr(3,3,3)=0
   lxyz_arr(1,4,3)=3 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=1 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=1 ; lxyz_arr(2,6,3)=1 ; lxyz_arr(3,6,3)=4
   lxyz_arr(1,7,3)=5 ; lxyz_arr(2,7,3)=1 ; lxyz_arr(3,7,3)=2
   lxyz_arr(1,8,3)=3 ; lxyz_arr(2,8,3)=3 ; lxyz_arr(3,8,3)=2
   lxyz_arr(1,9,3)=1 ; lxyz_arr(2,9,3)=5 ; lxyz_arr(3,9,3)=2
   lxyz_arr(1,10,3)=3 ; lxyz_arr(2,10,3)=1 ; lxyz_arr(3,10,3)=4
   lxyz_arr(1,11,3)=1 ; lxyz_arr(2,11,3)=3 ; lxyz_arr(3,11,3)=4
   lxyz_arr(1,12,3)=1 ; lxyz_arr(2,12,3)=1 ; lxyz_arr(3,12,3)=6
   fac_arr(1,3)=0.05757780227833739750818811d0
   fac_arr(2,3)=0.1151556045566747950163762d0
   fac_arr(3,3)=0.05757780227833739750818811d0
   fac_arr(4,3)=0.3454668136700243850491286d0
   fac_arr(5,3)=0.3454668136700243850491286d0
   fac_arr(6,3)=0.2878890113916869875409405d0
   fac_arr(7,3)=-0.05757780227833739750818811/rhol**2d0
   fac_arr(8,3)=-0.1151556045566747950163762/rhol**2d0
   fac_arr(9,3)=-0.05757780227833739750818811/rhol**2d0
   fac_arr(10,3)=-0.1151556045566747950163762/rhol**2d0
   fac_arr(11,3)=-0.1151556045566747950163762/rhol**2d0
   fac_arr(12,3)=-0.05757780227833739750818811/rhol**2d0
else
   stop 'PSP format error'
end if
end subroutine calc_coeff_derproj
