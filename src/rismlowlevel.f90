!!****f* BigDFT/rsimslowlevel
!! DESCRIPTION
!!   In this file, we have the analytic routines for the calculation of the overlap of short-range functions
!!***


!!****f* BigDFT/kinetic_overlap_h
!! FUNCTION
!!   Overlap kinetic matrix between two different basis structures
!!   the kinetic operator is applicated on the A basis structure
!!   The basis structure is supposed to be based on s-functions times Hermite polynomials
!! SOURCE
!!
subroutine kinetic_overlap_h(A,B,ovrlp)
  use module_base
  use module_types
  implicit none
  type(gaussian_basis), intent(in) :: A,B
  real(gp), dimension(A%ncoeff,B%ncoeff) :: ovrlp 
  !only lower triangular part for A%ncoeff=B%ncoeff
  !local variables
  integer, parameter :: niw=304,nrw=304
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
        if (lA /= 1 ) stop 'only s case is allowed for Hermite polynomials, A' 
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
                 if (lA /= 1 ) stop 'only s case is allowed for Hermite polynomials, B' 
                 do mB=1,2*lB-1
                    jovrlp=jovrlp+1
                    if (jovrlp >= iovrlp .and. A%ncoeff == B%ncoeff .or. &
                         A%ncoeff /= B%ncoeff ) then
                       call kineticovrlp_h(A%xp(iexpo),A%psiat(iexpo),&
                            B%xp(jexpo),B%psiat(jexpo),&
                            ngA,ngB,lA,isat,lB,jsat,dx,dy,dz,&
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
  
END SUBROUTINE kinetic_overlap_h
!!***


!!****f* BigDFT/kineticovrlp_h
!! FUNCTION
!!   Calculates the scalar product between two shells
!!   by considering only the nonzero coefficients
!!   actual building block for calculating overlap matrix
!!   inserted work arrays for calculation
!!   Only Hermite polynomials of r^2 are to be considered
!! SOURCE
!!
subroutine kineticovrlp_h(expo1,coeff1,expo2,coeff2,ng1,ng2,l1,ih1,l2,ih2,dx,dy,dz,&
     niw,nrw,iw,rw,ovrlp)
  use module_base
  implicit none
  integer, intent(in) :: ng1,ng2,l1,ih1,l2,ih2,niw,nrw
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
        call kinprod_h(a1,a2,dx,dy,dz,l1,ih1,l2,ih2,niw,nrw,iw,rw,govrlpr)
        govrlpr=c1*govrlpr*c2
        !print *,c1,c2,govrlpr
        ovrlp=ovrlp+govrlpr
     end do
  end do
  
END SUBROUTINE kineticovrlp_h
!!***


!!****f* BigDFT/kinprod_h
!! FUNCTION
!!  Kinetic overlap between gaussians, based on cartesian coordinates
!!  calculates a dot product between two differents gaussians times spherical harmonics
!!  only hermite polynomials
!! SOURCE
!!
subroutine kinprod_h(a1,a2,dx,dy,dz,l1,ih1,l2,ih2,niw,nrw,iw,rw,ovrlp)
  use module_base
  implicit none
  integer, intent(in) :: l1,l2,ih1,ih2,niw,nrw 
  real(gp), intent(in) :: a1,a2,dx,dy,dz
  integer, dimension(niw) :: iw !work array of the exponents of the two polynomials
  real(gp), dimension(nrw) :: rw !work array of the polynomials coefficients 
  real(gp), intent(out) :: ovrlp
  !local variables
  integer, parameter :: nx=48
  integer :: n1,n2,i1,i2,px,py,pz,qx,qy,qz
  real(gp) :: fx,fy,fz,fa,fb,govrlp,kinovrlp,d2fx,d2fy,d2fz

  !calculate the coefficients of the hermite polynomials
  !calculates the number of different couples
  call calc_coeff_hermite_r2(l1,ih1,nx,n1,&
       iw(1),iw(nx+1),iw(2*nx+1),rw(1))
  call calc_coeff_hermite_r2(l2,ih2,nx,n2,&
       iw(3*nx+1),iw(4*nx+1),iw(5*nx+1),rw(n1+1))

  ovrlp=0.0_gp
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
 
END SUBROUTINE kinprod_h
!!***


!!****f* BigDFT/calc_coeff_hermite_r2
!! FUNCTION
!!   
!!
!! SOURCE
!!
subroutine calc_coeff_hermite_r2(l,ih,nterm_max,nterm,lx,ly,lz,fac_arr)
  use module_base
  implicit none
  integer, intent(in) :: l,ih,nterm_max
  integer, intent(out) :: nterm
  integer, dimension(nterm_max), intent(out) :: lx,ly,lz
  real(gp), dimension(nterm_max), intent(out) :: fac_arr
  !local variables
  integer, parameter :: norder_max=16
  integer :: ntermr,iterm
  integer, dimension(norder_max) :: lr
  integer, dimension(norder_max) :: nfac_arr

  if (l /=1) stop 'error coeff_hermite'

  !take the coefficients for the herimt polynomials
  select case (ih)
  case (1)
     ntermr=1
     lr(1)=0
     nfac_arr(1)=1
  case (2)
     ntermr=1
     lr(1)=1
     nfac_arr(1)=1
  case (3)
     ntermr=2
     lr(1)=0
     lr(1)=2
     nfac_arr(1)=1
     nfac_arr(2)=-1
     !to be continued with all the coefficients
  case default
     stop 'order of hermite polynomial not provided'
  end select

  !and write them on the output arrays 
  nterm=1
  do iterm=1,ntermr
     if (lr(iterm) == 0) then
        lx(nterm)=0
        ly(nterm)=0
        lz(nterm)=0
        fac_arr(nterm)=real(nfac_arr(iterm),gp)
        nterm=nterm+1
     else
        lx(nterm)=2*lr(iterm) ;lx(nterm+1)=0          ;lx(nterm+2)=0
        ly(nterm)=0           ;ly(nterm+1)=2*lr(iterm);ly(nterm+2)=0          
        lz(nterm)=0           ;lz(nterm+1)=0          ;lz(nterm+2)=2*lr(iterm)          
        fac_arr(nterm)=real(nfac_arr(iterm),gp)
        fac_arr(nterm+1)=fac_arr(nterm)
        fac_arr(nterm+2)=fac_arr(nterm)
        nterm=nterm+3
     end if
  end do
  !last term
  nterm=nterm-1

end subroutine calc_coeff_hermite_r2
!!***

subroutine gaussians_to_wavelets_new_h(iproc,nproc,lr,orbs,hx,hy,hz,G,wfn_gau,psi)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbs
  type(gaussian_basis), intent(in) :: G
  real(wp), dimension(G%ncoeff,orbs%nspinor,orbs%norbp), intent(in) :: wfn_gau
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(out) :: psi
  !local variables
  integer :: iorb,ierr,ispinor,ncplx
  real(dp) :: normdev,tt,scpr,totnorm
  real(gp) :: kx,ky,kz

  if(iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no')&
       'Writing wavefunctions in wavelet form...'

  normdev=0.0_dp
  tt=0.0_dp
  do iorb=1,orbs%norbp
     !features of the k-point ikpt
     kx=orbs%kpts(1,orbs%iokpt(iorb))
     ky=orbs%kpts(2,orbs%iokpt(iorb))
     kz=orbs%kpts(3,orbs%iokpt(iorb))

     !evaluate the complexity of the k-point
     if (kx**2 + ky**2 + kz**2 == 0.0_gp) then
        ncplx=1
     else
        ncplx=2
     end if
     totnorm=0.0_dp
     do ispinor=1,orbs%nspinor,ncplx
        !print *,'start',ispinor,ncplx,iorb,orbs%nspinor
        !the Block wavefunctions are exp(-Ikr) psi(r) (with MINUS k)
        call gaussians_to_wavelets_orb(ncplx,lr,hx,hy,hz,kx,ky,kz,G,&
             wfn_gau(1,ispinor,iorb),psi(1,ispinor,iorb))
        !print *,'end',ispinor,ncplx,iorb,orbs%nspinor
        call wnrm_wrap(ncplx,lr%wfd%nvctr_c,lr%wfd%nvctr_f,psi(1,ispinor,iorb),scpr) 
        totnorm=totnorm+scpr
     end do
     !write(*,'(1x,a,i5,1pe14.7,i3)')'norm of orbital ',iorb,totnorm,ncplx
     do ispinor=1,orbs%nspinor
        call wscal_wrap(lr%wfd%nvctr_c,lr%wfd%nvctr_f,real(1.0_dp/sqrt(totnorm),wp),&
             psi(1,ispinor,iorb))
     end do
     tt=max(tt,abs(1.0_dp-totnorm))
     !print *,'iorb,norm',totnorm
  end do

  if (iproc ==0  .and. verbose > 1) write(*,'(1x,a)')'done.'
  !renormalize the orbitals
  !calculate the deviation from 1 of the orbital norm
  if (nproc > 1) then
     call MPI_REDUCE(tt,normdev,1,mpidtypd,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  else
     normdev=tt
  end if
  if (iproc ==0) write(*,'(1x,a,1pe12.2)')&
       'Deviation from normalization of the imported orbitals',normdev

END SUBROUTINE gaussians_to_wavelets_new_h


subroutine gaussians_to_wavelets_orb_h(ncplx,lr,hx,hy,hz,kx,ky,kz,G,wfn_gau,psi)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ncplx
  real(gp), intent(in) :: hx,hy,hz,kx,ky,kz
  type(locreg_descriptors), intent(in) :: lr
  type(gaussian_basis), intent(in) :: G
  real(wp), dimension(G%ncoeff), intent(in) :: wfn_gau
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*ncplx), intent(out) :: psi
  !local variables
  character(len=*), parameter :: subname='gaussians_to_wavelets_orb'
  integer, parameter :: nterm_max=48,maxsizeKB=2048,nw=65536
  logical :: perx,pery,perz
  integer :: i_stat,i_all,ishell,iexpo,icoeff,iat,isat,ng,l,m,i,nterm,ig
  integer :: nterms_max,nterms,iscoeff,iterm,n_gau,ml1,mu1,ml2,mu2,ml3,mu3
  real(gp) :: rx,ry,rz,gau_a
  integer, dimension(nterm_max) :: lx,ly,lz
  real(gp), dimension(nterm_max) :: fac_arr
  real(wp), allocatable, dimension(:,:,:) :: work
  real(wp), allocatable, dimension(:,:,:,:) :: wx,wy,wz

  !calculate nterms_max:
  !allows only maxsizeKB per one-dimensional array
  !(for a grid of dimension 100 nterms_max=655)
  !bu with at least ngx*nterm_max ~= 100 elements
  nterms_max=max(maxsizeKB*1024/(2*ncplx*max(lr%d%n1,lr%d%n2,lr%d%n3)),100)

  allocate(work(0:nw,2,2+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)

  allocate(wx(ncplx,0:lr%d%n1,2,nterms_max+ndebug),stat=i_stat)
  call memocc(i_stat,wx,'wx',subname)
  allocate(wy(ncplx,0:lr%d%n2,2,nterms_max+ndebug),stat=i_stat)
  call memocc(i_stat,wy,'wy',subname)
  allocate(wz(ncplx,0:lr%d%n3,2,nterms_max+ndebug),stat=i_stat)
  call memocc(i_stat,wz,'wz',subname)

  !conditions for periodicity in the three directions
  perx=(lr%geocode /= 'F')
  pery=(lr%geocode == 'P')
  perz=(lr%geocode /= 'F')

  !initialize the wavefunction
  call razero((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*ncplx,psi)

  !calculate the number of terms for this orbital
  nterms=0
  !loop over the atoms
  ishell=0
  iexpo=1
  icoeff=1
  iscoeff=1
  iterm=1
  do iat=1,G%nat
     rx=G%rxyz(1,iat)
     ry=G%rxyz(2,iat)
     rz=G%rxyz(3,iat)
     !loop over the number of shells of the atom type
     do isat=1,G%nshell(iat)
        ishell=ishell+1
        !the degree of contraction of the basis function
        !is the same as the ng value of the createAtomicOrbitals routine
        ng=G%ndoc(ishell)
        !angular momentum of the basis set(shifted for compatibility with BigDFT routines
        l=G%nam(ishell)
        if (l/=1) stop 'error l gautowav_h'
        !print *,iproc,iat,ishell,G%nam(ishell),G%nshell(iat)
        !multiply the values of the gaussian contraction times the orbital coefficient

        do m=1,2*l-1
           call calc_coeff_inguess(l,isat,nterm_max,nterm,lx,ly,lz,fac_arr)
           !control whether the basis element may be
           !contribute to some of the orbital of the processor
           if (wfn_gau(icoeff) /= 0.0_wp) then
              if (nterms + nterm*ng > nterms_max) then
                 !accumulate wavefuncton
                 call wfn_from_tensprod(lr,ncplx,nterms,wx,wy,wz,psi)
                 iterm=1
                 nterms=0
              end if
              !assign the arrays
              !make sure that the coefficients returned by 
              !gauss_to_daub are zero outside [ml:mr] 
              do ig=1,ng
                 do i=1,nterm
                    !print *,iat,ig,i,fac_arr(i),wfn_gau(icoeff),G%xp(iexpo+ig-1)
                    gau_a=G%xp(iexpo+ig-1)
                    n_gau=lx(i)
                    !print *,'x',gau_a,nterm,ncplx,kx,ky,kz,ml1,mu1,lr%d%n1
                    call gauss_to_daub_k(hx,kx*hx,ncplx,fac_arr(i),rx,gau_a,n_gau,&
                         lr%d%n1,ml1,mu1,&
                         wx(1,0,1,iterm),work,nw,perx) 
                    n_gau=ly(i)
                    !print *,'y',ml2,mu2,lr%d%n2
                    call gauss_to_daub_k(hy,ky*hy,ncplx,wfn_gau(icoeff),ry,gau_a,n_gau,&
                         lr%d%n2,ml2,mu2,&
                         wy(1,0,1,iterm),work,nw,pery) 
                    n_gau=lz(i) 
                    !print *,'z',ml3,mu3,lr%d%n3
                    call gauss_to_daub_k(hz,kz*hz,ncplx,G%psiat(iexpo+ig-1),rz,gau_a,n_gau,&
                         lr%d%n3,ml3,mu3,&
                         wz(1,0,1,iterm),work,nw,perz)
                    iterm=iterm+1
                 end do
              end do
              nterms=nterms+nterm*ng
           end if
           icoeff=icoeff+1
        end do
        iexpo=iexpo+ng
     end do
  end do

  call gaudim_check(iexpo,icoeff,ishell,G%nexpo,G%ncoeff,G%nshltot)

  !accumulate wavefuncton
  call wfn_from_tensprod(lr,ncplx,nterms,wx,wy,wz,psi)
!psi=1.d0
  i_all=-product(shape(wx))*kind(wx)
  deallocate(wx,stat=i_stat)
  call memocc(i_stat,i_all,'wx',subname)
  i_all=-product(shape(wy))*kind(wy)
  deallocate(wy,stat=i_stat)
  call memocc(i_stat,i_all,'wy',subname)
  i_all=-product(shape(wz))*kind(wz)
  deallocate(wz,stat=i_stat)
  call memocc(i_stat,i_all,'wz',subname)

  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)



END SUBROUTINE gaussians_to_wavelets_orb_h

!!****f* BigDFT/gaussian_overlap
!! FUNCTION
!!   Overlap matrix between two different basis structures
!!   The first one is a gaussian hermite basis
!! SOURCE
!!
subroutine gaussian_overlap_h(A,B,ovrlp)
  use module_base
  use module_types
  implicit none
  type(gaussian_basis), intent(in) :: A,B
  real(gp), dimension(A%ncoeff,B%ncoeff) :: ovrlp 
  !only lower triangular part for A%ncoeff=B%ncoeff
  !local variables
  integer, parameter :: niw=302,nrw=302
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
        if (lA /= 1) stop 'only s terms supported, gauovrlp'
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
                    if ((jovrlp >= iovrlp .and. A%ncoeff == B%ncoeff) .or. &
                         A%ncoeff /= B%ncoeff ) then
                       call gbasovrlp_h(A%xp(iexpo),A%psiat(iexpo),&
                            B%xp(jexpo),B%psiat(jexpo),&
                            ngA,ngB,lA,isat,lB,mB,dx,dy,dz,&
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
  
END SUBROUTINE gaussian_overlap_h
!!***


!!****f* BigDFT/gbasovrlp
!! FUNCTION
!!   Calculates the scalar product between two shells
!!   by considering only the nonzero coefficients
!!   actual building block for calculating overlap matrix
!!   inserted work arrays for calculation
!!   The first are Hermite polynomial basis
!! SOURCE
!!
subroutine gbasovrlp_h(expo1,coeff1,expo2,coeff2,ng1,ng2,l1,ih1,l2,m2,dx,dy,dz,&
     niw,nrw,iw,rw,ovrlp)
  use module_base
  implicit none
  integer, intent(in) :: ng1,ng2,l1,ih1,l2,m2,niw,nrw
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
        call gprod_h(a1,a2,dx,dy,dz,l1,ih1,l2,m2,niw,nrw,iw,rw,govrlpr)
        govrlpr=c1*govrlpr*c2
        !print *,c1,c2,govrlpr
        ovrlp=ovrlp+govrlpr
     end do
  end do
  
END SUBROUTINE gbasovrlp_h
!!***
!calculates a dot product between two differents gaussians times spherical harmonics
!the first one is supposed to be hermite polynomial matrix
subroutine gprod_h(a1,a2,dx,dy,dz,l1,ih1,l2,m2,niw,nrw,iw,rw,ovrlp)
  use module_base
  implicit none
  integer, intent(in) :: l1,l2,ih1,m2,niw,nrw 
  real(gp), intent(in) :: a1,a2,dx,dy,dz
  integer, dimension(niw) :: iw !work array of the exponents of the two polynomials
  real(gp), dimension(nrw) :: rw !work array of the polynomials coefficients 
  real(gp), intent(out) :: ovrlp
  !local variables
  integer, parameter :: nx=48
  integer :: n1,n2,i1,i2,px,py,pz,qx,qy,qz
  real(gp) :: fx,fy,fz,fa,fb,govrlp

  !calculates the number of different couples
  call calc_coeff_hermite_r2(l1,ih1,nx,n1,&
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
 
END SUBROUTINE gprod_h
