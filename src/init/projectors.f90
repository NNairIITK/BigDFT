!fill the proj array with the PSP projectors or their derivatives, following idir value
subroutine fill_projectors(geocode,iproc,n1,n2,n3,hx,hy,hz,cpmult,fpmult,at,rxyz,radii_cf,&
     nlpspd,proj,idir)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,n1,n2,n3,idir
  real(gp), intent(in) :: hx,hy,hz,cpmult,fpmult
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), dimension(at%ntypes,2), intent(in) :: radii_cf
  real(wp), dimension(nlpspd%nprojel), intent(out) :: proj
  !local variables
  integer, parameter :: nterm_max=20 !if GTH nterm_max=4
  integer :: istart_c,istart_f,mvctr_c,mvctr_f
  integer :: nl1_c,nl1_f,nl2_c,nl2_f,nl3_c,nl3_f,nu1_c,nu1_f,nu2_c,nu2_f,nu3_c,nu3_f
  integer :: iat,i,l,m,iproj,ityp,nterm,nwarnings,iterm
  real(gp) :: fpi,factor,gau_a,rx,ry,rz
  real(dp) :: scpr
  integer, dimension(3) :: nterm_arr
  integer, dimension(nterm_max) :: lx,ly,lz
  integer, dimension(3,nterm_max,3) :: lxyz_arr
  real(gp), dimension(nterm_max) :: factors
  real(gp), dimension(nterm_max,3) :: fac_arr

  if (iproc.eq.0 .and. nlpspd%nproj /=0 .and. idir==0) write(*,'(1x,a)',advance='no') &
       'Calculating wavelets expansion of projectors...'
  !warnings related to the projectors norm
  nwarnings=0
  !allocate these vectors up to the maximum size we can get
  istart_c=1
  iproj=0
  fpi=(4.0_gp*atan(1.0_gp))**(-.75_gp)
  do iat=1,at%nat
     rx=rxyz(1,iat) 
     ry=rxyz(2,iat) 
     rz=rxyz(3,iat)
     ityp=at%iatype(iat)

     mvctr_c=nlpspd%nvctr_p(2*iat-1)-nlpspd%nvctr_p(2*iat-2)
     mvctr_f=nlpspd%nvctr_p(2*iat  )-nlpspd%nvctr_p(2*iat-1)

     nl1_c=nlpspd%nboxp_c(1,1,iat)
     nl2_c=nlpspd%nboxp_c(1,2,iat)
     nl3_c=nlpspd%nboxp_c(1,3,iat)
     nl1_f=nlpspd%nboxp_f(1,1,iat)
     nl2_f=nlpspd%nboxp_f(1,2,iat)
     nl3_f=nlpspd%nboxp_f(1,3,iat)

     nu1_c=nlpspd%nboxp_c(2,1,iat)
     nu2_c=nlpspd%nboxp_c(2,2,iat)
     nu3_c=nlpspd%nboxp_c(2,3,iat)
     nu1_f=nlpspd%nboxp_f(2,1,iat)
     nu2_f=nlpspd%nboxp_f(2,2,iat)
     nu3_f=nlpspd%nboxp_f(2,3,iat)

     !decide the loop bounds
     do l=1,4 !generic case, also for HGHs (for GTH it will stop at l=2)
        do i=1,3 !generic case, also for HGHs (for GTH it will stop at i=2)
           if (at%psppar(l,i,ityp) /= 0.0_gp) then
              gau_a=at%psppar(l,0,ityp)
              factor=sqrt(2.0_gp)*fpi/(sqrt(gau_a)**(2*(l-1)+4*i-1))
              do m=1,2*l-1
                 istart_f=istart_c+mvctr_c

                 if (idir==0) then !normal projector calculation case
                    call calc_coeff_proj(l,i,m,nterm_max,nterm,lx,ly,lz,factors)
                    
                    factors(1:nterm)=factor*factors(1:nterm)
                 else !calculation of projector derivative
                    call calc_coeff_derproj(l,i,m,nterm_max,gau_a,nterm_arr,lxyz_arr,fac_arr)

                    nterm=nterm_arr(idir)
                    do iterm=1,nterm
                       factors(iterm)=factor*fac_arr(iterm,idir)
                       lx(iterm)=lxyz_arr(1,iterm,idir)
                       ly(iterm)=lxyz_arr(2,iterm,idir)
                       lz(iterm)=lxyz_arr(3,iterm,idir)
                    end do                    
                 end if

                 call crtproj(geocode,iproc,nterm,n1,n2,n3,nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c, &
                      & nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,radii_cf(at%iatype(iat),2), & 
                      & cpmult,fpmult,hx,hy,hz,gau_a,factors,rx,ry,rz,lx,ly,lz, & 
                      & mvctr_c,mvctr_f,proj(istart_c),proj(istart_f))

                 iproj=iproj+1

                 ! testing
                 if (idir == 0) then
                    call wnrm(mvctr_c,mvctr_f,proj(istart_c), &
                         & proj(istart_f),scpr)
                    if (abs(1.d0-scpr) > 1.d-2) then
                       if (abs(1.d0-scpr) > 1.d-1) then
                          if (iproc == 0) then
                             write(*,'(1x,a)')'error found!'
                             write(*,'(1x,a,i4,a,a6,a,i1,a,i1,a,f4.3)')&
                                  'The norm of the nonlocal PSP for atom n=',iat,&
                                  ' (',trim(at%atomnames(at%iatype(iat))),&
                                  ') labeled by l=',l,' m=',m,' is ',scpr
                             write(*,'(1x,a)')&
                                  'while it is supposed to be about 1.0. Control PSP data or reduce grid spacing.'
                          end if
                          stop
                       else
                          nwarnings=nwarnings+1
                       end if
                    end if
                    !do iterm=1,nterm
                    !   if (iproc.eq.0) write(*,'(1x,a,i0,1x,a,1pe10.3,3(1x,i0))') &
                    !        'projector: iat,atomname,gau_a,lx,ly,lz ', & 
                    !        iat,trim(at%atomnames(at%iatype(iat))),gau_a,lx(iterm),ly(iterm),lz(iterm)
                    !enddo
                 end if
                 !end testing

                 istart_c=istart_f+7*mvctr_f
                 if (istart_c > nlpspd%nprojel+1) stop 'istart_c > nprojel+1'
              enddo
           endif
        enddo
     enddo
  enddo
  if (iproj.ne.nlpspd%nproj) stop 'incorrect number of projectors created'
  ! projector part finished

  if (iproc == 0 .and. nlpspd%nproj /=0 .and. idir == 0) then
     if (nwarnings == 0) then
        write(*,'(1x,a)')'done.'
     else
        write(*,'(1x,a,i0,a)')'found ',nwarnings,' warnings.'
        write(*,'(1x,a)')'Some projectors may be too rough.'
        write(*,'(1x,a,f6.3)')&
             'Consider the possibility of reducing hgrid for having a more accurate run.'
     end if
  end if

end subroutine fill_projectors


subroutine numb_proj(ityp,ntypes,psppar,npspcode,mproj)
  ! Determines the number of projectors (valid for GTH and HGH pseudopotentials)
  implicit real(kind=8) (a-h,o-z)
  dimension psppar(0:4,0:6,ntypes),npspcode(ntypes)

  mproj=0
  if (npspcode(ityp) == 2) then !GTH
     do l=1,2 
        do i=1,2 
           if (psppar(l,i,ityp).ne.0.d0) mproj=mproj+2*l-1
        enddo
     enddo
  else if (npspcode(ityp) == 3 .or. npspcode(ityp) == 10) then !HGH and HGH-K
     do l=1,4 
        do i=1,3 
           if (psppar(l,i,ityp).ne.0.d0) mproj=mproj+2*l-1
        enddo
     enddo
  end if
END SUBROUTINE numb_proj

subroutine crtproj(geocode,iproc,nterm,n1,n2,n3, & 
     nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,  & 
     radius_f,cpmult,fpmult,hx,hy,hz,gau_a,fac_arr,rx,ry,rz,lx,ly,lz, & 
     mvctr_c,mvctr_f,proj_c,proj_f)
  ! returns the compressed form of a Gaussian projector 
  ! x^lx * y^ly * z^lz * exp (-1/(2*gau_a^2) *((x-cntrx)^2 + (y-cntry)^2 + (z-cntrz)^2 ))
  ! in the arrays proj_c, proj_f
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nterm,n1,n2,n3,mvctr_c,mvctr_f
  integer, intent(in) :: nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f
  real(kind=8), intent(in) :: radius_f,cpmult,fpmult,hx,hy,hz,gau_a,rx,ry,rz
  integer, dimension(nterm), intent(in) :: lx,ly,lz
  real(kind=8), dimension(nterm), intent(in) :: fac_arr
  real(kind=8), dimension(mvctr_c), intent(out) :: proj_c
  real(kind=8), dimension(7,mvctr_f), intent(out) :: proj_f
  !local variables
  character(len=*), parameter :: subname='crtproj'
  integer, parameter :: nw=16000
  logical :: perx,pery,perz !variables controlling the periodicity in x,y,z
  integer :: iterm,n_gau,ml1,ml2,ml3,mu1,mu2,mu3,i1,i2,i3,mvctr,i_all,i_stat,j1,j2,j3
  real(kind=8) :: rad_c,rad_f,factor,err_norm,dz2,dy2,dx2,te,d2,cpmult_max,fpmult_max
  real(kind=8), allocatable, dimension(:,:,:) :: wprojx,wprojy,wprojz
  real(kind=8), allocatable, dimension(:,:) :: work

  allocate(wprojx(0:n1,2,nterm+ndebug),stat=i_stat)
  call memocc(i_stat,wprojx,'wprojx',subname)
  allocate(wprojy(0:n2,2,nterm+ndebug),stat=i_stat)
  call memocc(i_stat,wprojy,'wprojy',subname)
  allocate(wprojz(0:n3,2,nterm+ndebug),stat=i_stat)
  call memocc(i_stat,wprojz,'wprojz',subname)
  allocate(work(0:nw,2+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)

  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')


  rad_c=radius_f*cpmult
  rad_f=radius_f*fpmult

  !print *,''
  !print *,'limits',nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f

  ! make sure that the coefficients returned by CALL GAUSS_TO_DAUB are zero outside [ml:mr] 
  err_norm=0.d0 
  do iterm=1,nterm
     factor=fac_arr(iterm)
     n_gau=lx(iterm) 
     CALL GAUSS_TO_DAUB(hx,factor,rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1,iterm),te,work,nw,perx) 
     err_norm=max(err_norm,te) 
     n_gau=ly(iterm) 
     CALL GAUSS_TO_DAUB(hy,1.d0,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1,iterm),te,work,nw,pery) 
     err_norm=max(err_norm,te) 
     n_gau=lz(iterm) 
     CALL GAUSS_TO_DAUB(hz,1.d0,rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1,iterm),te,work,nw,perz) 
     err_norm=max(err_norm,te) 
     if (iproc.eq.0 .and. geocode == 'F')  then
        if (ml1.gt.min(nl1_c,nl1_f)) write(*,*) 'Projector box larger than needed: ml1'
        if (ml2.gt.min(nl2_c,nl2_f)) write(*,*) 'Projector box larger than needed: ml2'
        if (ml3.gt.min(nl3_c,nl3_f)) write(*,*) 'Projector box larger than needed: ml3'
        if (mu1.lt.max(nu1_c,nu1_f)) write(*,*) 'Projector box larger than needed: mu1'
        if (mu2.lt.max(nu2_c,nu2_f)) write(*,*) 'Projector box larger than needed: mu2'
        if (mu3.lt.max(nu3_c,nu3_f)) write(*,*) 'Projector box larger than needed: mu3'
!!$        !approximate maximum value for cpmult,fpmult
!!$        cpmult_max=max((nu1_c-nl1_c)*hx*0.5d0/radius_f,(nu2_c-nl2_c)*hy*0.5d0/radius_f,(nu3_c-nl3_c)*hz*0.5d0/radius_f)
!!$        fpmult_max=max((nu1_f-nl1_f)*hx*0.5d0/radius_f,(nu2_f-nl2_f)*hy*0.5d0/radius_f,(nu3_f-nl3_f)*hz*0.5d0/radius_f)
!!$        print *,'cpmult_max,fpmult_max=',cpmult_max,fpmult_max
     endif

     !print *,iterm,ml1,mu1,ml2,mu2,ml3,mu3,n1,n2,n3,rad_c,rad_f

  end do
  !        if (iproc.eq.0) write(*,*) 'max err_norm ',err_norm

  ! First term: coarse projector components
  !write(*,*) 'rad_c=',rad_c
  
  mvctr=0

  do i3=nl3_c,nu3_c
     dz2=d2(i3,hz,n3,rz,perz)
     do i2=nl2_c,nu2_c
        dy2=d2(i2,hy,n2,ry,pery)
        do i1=nl1_c,nu1_c
           dx2=d2(i1,hx,n1,rx,perx)
           if (dx2+(dy2+dz2) <= rad_c**2) then
              mvctr=mvctr+1
              proj_c(mvctr)=wprojx(i1,1,1)*wprojy(i2,1,1)*wprojz(i3,1,1)
           endif
        enddo
     enddo
  enddo

  if (mvctr /=  mvctr_c) then
     write(*,'(1x,a,i0,1x,i0)')' ERROR: mvctr >< mvctr_c ',mvctr,mvctr_c
     stop
  end if

  ! First term: fine projector components
  mvctr=0
  do i3=nl3_f,nu3_f
     dz2=d2(i3,hz,n3,rz,perz)
     do i2=nl2_f,nu2_f
        dy2=d2(i2,hy,n2,ry,pery)
        do i1=nl1_f,nu1_f
           dx2=d2(i1,hx,n1,rx,perx)
           if (dx2+(dy2+dz2) <= rad_f**2) then
              mvctr=mvctr+1
              proj_f(1,mvctr)=wprojx(i1,2,1)*wprojy(i2,1,1)*wprojz(i3,1,1)
              proj_f(2,mvctr)=wprojx(i1,1,1)*wprojy(i2,2,1)*wprojz(i3,1,1)
              proj_f(3,mvctr)=wprojx(i1,2,1)*wprojy(i2,2,1)*wprojz(i3,1,1)
              proj_f(4,mvctr)=wprojx(i1,1,1)*wprojy(i2,1,1)*wprojz(i3,2,1)
              proj_f(5,mvctr)=wprojx(i1,2,1)*wprojy(i2,1,1)*wprojz(i3,2,1)
              proj_f(6,mvctr)=wprojx(i1,1,1)*wprojy(i2,2,1)*wprojz(i3,2,1)
              proj_f(7,mvctr)=wprojx(i1,2,1)*wprojy(i2,2,1)*wprojz(i3,2,1)
           endif
        enddo
     enddo
  enddo
  if (mvctr /= mvctr_f) then
     write(*,'(1x,a,i0,1x,i0)')' ERROR: mvctr >< mvctr_f ',mvctr,mvctr_f
     stop 
  end if

  do iterm=2,nterm

     ! Other terms: coarse projector components
     mvctr=0
     do i3=nl3_c,nu3_c
        dz2=d2(i3,hz,n3,rz,perz)
        do i2=nl2_c,nu2_c
           dy2=d2(i2,hy,n2,ry,pery)
           do i1=nl1_c,nu1_c
              dx2=d2(i1,hx,n1,rx,perx)
              if (dx2+(dy2+dz2).le.rad_c**2) then
              mvctr=mvctr+1
                 proj_c(mvctr)=proj_c(mvctr)+wprojx(i1,1,iterm)*wprojy(i2,1,iterm)*wprojz(i3,1,iterm)
              endif
           enddo
        enddo
     enddo

     ! Other terms: fine projector components
     mvctr=0
     do i3=nl3_f,nu3_f
        dz2=d2(i3,hz,n3,rz,perz)
        do i2=nl2_f,nu2_f
           dy2=d2(i2,hy,n2,ry,pery)
           do i1=nl1_f,nu1_f
              dx2=d2(i1,hx,n1,rx,perx)
              if (dx2+(dy2+dz2).le.rad_f**2) then
                 mvctr=mvctr+1
                 proj_f(1,mvctr)=&
                      proj_f(1,mvctr)+wprojx(i1,2,iterm)*wprojy(i2,1,iterm)*wprojz(i3,1,iterm)
                 proj_f(2,mvctr)=&
                      proj_f(2,mvctr)+wprojx(i1,1,iterm)*wprojy(i2,2,iterm)*wprojz(i3,1,iterm)
                 proj_f(3,mvctr)=&
                      proj_f(3,mvctr)+wprojx(i1,2,iterm)*wprojy(i2,2,iterm)*wprojz(i3,1,iterm)
                 proj_f(4,mvctr)=&
                      proj_f(4,mvctr)+wprojx(i1,1,iterm)*wprojy(i2,1,iterm)*wprojz(i3,2,iterm)
                 proj_f(5,mvctr)=&
                      proj_f(5,mvctr)+wprojx(i1,2,iterm)*wprojy(i2,1,iterm)*wprojz(i3,2,iterm)
                 proj_f(6,mvctr)=&
                      proj_f(6,mvctr)+wprojx(i1,1,iterm)*wprojy(i2,2,iterm)*wprojz(i3,2,iterm)
                 proj_f(7,mvctr)=&
                      proj_f(7,mvctr)+wprojx(i1,2,iterm)*wprojy(i2,2,iterm)*wprojz(i3,2,iterm)
              endif
           enddo
        enddo
     enddo

  enddo

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

end subroutine crtproj

function d2(i,h,n,r,periodic)
  implicit none
  logical, intent(in) :: periodic !determine whether the dimension is periodic or not
  integer, intent(in) :: i,n
  real(kind=8), intent(in) :: h,r
  real(kind=8) :: d2
  !local variables
  integer :: j

  if (periodic) then
     !take the point which has the minimum distance  comparing with the periodic image
     !which is i-+(n+1) if it is on the second or on the first half respectively
     j=i-(2*(i/(n/2+1))-1)*(n+1)
     d2=min((real(i,kind=8)*h-r)**2,(real(j,kind=8)*h-r)**2)
  else
     d2=(real(i,kind=8)*h-r)**2
  end if
  
end function d2

subroutine pregion_size(geocode,rxyz,radius,rmult,hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
  ! finds the size of the smallest subbox that contains a localization region made 
  ! out of atom centered spheres
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n1,n2,n3
  real(kind=8), intent(in) :: hx,hy,hz,rmult,radius
  real(kind=8), dimension(3), intent(in) :: rxyz
  integer, intent(out) :: nl1,nu1,nl2,nu2,nl3,nu3
  !local variables
  real(kind=8), parameter :: eps_mach=1.d-12
  real(kind=8) :: cxmax,cymax,czmax,cxmin,cymin,czmin,rad,onem

  rad=radius*rmult
  cxmax=rxyz(1)+rad ; cxmin=rxyz(1)-rad
  cymax=rxyz(2)+rad ; cymin=rxyz(2)-rad
  czmax=rxyz(3)+rad ; czmin=rxyz(3)-rad
  onem=1.d0-eps_mach

  nl1=ceiling(cxmin/hx - eps_mach)   
  nl2=ceiling(cymin/hy - eps_mach)   
  nl3=ceiling(czmin/hz - eps_mach)   
  nu1=floor(cxmax/hx + eps_mach)  
  nu2=floor(cymax/hy + eps_mach)  
  nu3=floor(czmax/hz + eps_mach)  

  !for non-free BC the projectors are allowed to be also outside the box
  if (geocode == 'F') then
     if (nl1 < 0)   stop 'nl1: projector region outside cell'
     if (nl2 < 0)   stop 'nl2: projector region outside cell'
     if (nl3 < 0)   stop 'nl3: projector region outside cell'
     if (nu1 > n1)   stop 'nu1: projector region outside cell'
     if (nu2 > n2)   stop 'nu2: projector region outside cell'
     if (nu3 > n3)   stop 'nu3: projector region outside cell'
  else
     !correct the extremes if they run outside the box
     if (nl1 < 0 .or. nu1 > n1) then
        nl1=0
        nu1=n1
     end if
     if (nl2 < 0 .or. nu2 > n2) then
        nl2=0
        nu2=n2
     end if
     if (nl3 < 0 .or. nu3 > n3) then
        nl3=0
        nu3=n3
     end if
  end if

END SUBROUTINE pregion_size

subroutine calc_coeff_proj(l,i,m,nterm_max,nterm,lx,ly,lz,fac_arr)
  
  implicit none
  integer, intent(in) :: l,i,m,nterm_max
  integer, intent(out) :: nterm
  integer, dimension(nterm_max), intent(out) :: lx,ly,lz
  real(kind=8), dimension(nterm_max), intent(out) :: fac_arr

  if (l.eq.1 .and. i.eq.1 .and. m.eq.1) then
     nterm=1
     lx(1)=0 ; ly(1)=0 ; lz(1)=0
     fac_arr(1)=0.7071067811865475244008444d0
  else if (l.eq.1 .and. i.eq.2 .and. m.eq.1) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=0.3651483716701107423046465d0
     fac_arr(2)=0.3651483716701107423046465d0
     fac_arr(3)=0.3651483716701107423046465d0
  else if (l.eq.1 .and. i.eq.3 .and. m.eq.1) then
     nterm=6
     lx(1)=4 ; ly(1)=0 ; lz(1)=0
     lx(2)=2 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=4 ; lz(3)=0
     lx(4)=2 ; ly(4)=0 ; lz(4)=2
     lx(5)=0 ; ly(5)=2 ; lz(5)=2
     lx(6)=0 ; ly(6)=0 ; lz(6)=4
     fac_arr(1)=0.09200874124564722903948358d0
     fac_arr(2)=0.1840174824912944580789672d0
     fac_arr(3)=0.09200874124564722903948358d0
     fac_arr(4)=0.1840174824912944580789672d0
     fac_arr(5)=0.1840174824912944580789672d0
     fac_arr(6)=0.09200874124564722903948358d0
  else if (l.eq.2 .and. i.eq.1 .and. m.eq.1) then
     nterm=1
     lx(1)=1 ; ly(1)=0 ; lz(1)=0
     fac_arr(1)=1.000000000000000000000000d0
  else if (l.eq.2 .and. i.eq.1 .and. m.eq.2) then
     nterm=1
     lx(1)=0 ; ly(1)=1 ; lz(1)=0
     fac_arr(1)=1.000000000000000000000000d0
  else if (l.eq.2 .and. i.eq.1 .and. m.eq.3) then
     nterm=1
     lx(1)=0 ; ly(1)=0 ; lz(1)=1
     fac_arr(1)=1.000000000000000000000000d0
  else if (l.eq.2 .and. i.eq.2 .and. m.eq.1) then
     nterm=3
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=0.3380617018914066310038473d0
     fac_arr(2)=0.3380617018914066310038473d0
     fac_arr(3)=0.3380617018914066310038473d0
  else if (l.eq.2 .and. i.eq.2 .and. m.eq.2) then
     nterm=3
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=1 ; lz(3)=2
     fac_arr(1)=0.3380617018914066310038473d0
     fac_arr(2)=0.3380617018914066310038473d0
     fac_arr(3)=0.3380617018914066310038473d0
  else if (l.eq.2 .and. i.eq.2 .and. m.eq.3) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=0 ; lz(3)=3
     fac_arr(1)=0.3380617018914066310038473d0
     fac_arr(2)=0.3380617018914066310038473d0
     fac_arr(3)=0.3380617018914066310038473d0
  else if (l.eq.2 .and. i.eq.3 .and. m.eq.1) then
     nterm=6
     lx(1)=5 ; ly(1)=0 ; lz(1)=0
     lx(2)=3 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=4 ; lz(3)=0
     lx(4)=3 ; ly(4)=0 ; lz(4)=2
     lx(5)=1 ; ly(5)=2 ; lz(5)=2
     lx(6)=1 ; ly(6)=0 ; lz(6)=4
     fac_arr(1)=0.06795295885835007261827187d0
     fac_arr(2)=0.1359059177167001452365437d0
     fac_arr(3)=0.06795295885835007261827187d0
     fac_arr(4)=0.1359059177167001452365437d0
     fac_arr(5)=0.1359059177167001452365437d0
     fac_arr(6)=0.06795295885835007261827187d0
  else if (l.eq.2 .and. i.eq.3 .and. m.eq.2) then
     nterm=6
     lx(1)=4 ; ly(1)=1 ; lz(1)=0
     lx(2)=2 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=5 ; lz(3)=0
     lx(4)=2 ; ly(4)=1 ; lz(4)=2
     lx(5)=0 ; ly(5)=3 ; lz(5)=2
     lx(6)=0 ; ly(6)=1 ; lz(6)=4
     fac_arr(1)=0.06795295885835007261827187d0
     fac_arr(2)=0.1359059177167001452365437d0
     fac_arr(3)=0.06795295885835007261827187d0
     fac_arr(4)=0.1359059177167001452365437d0
     fac_arr(5)=0.1359059177167001452365437d0
     fac_arr(6)=0.06795295885835007261827187d0
  else if (l.eq.2 .and. i.eq.3 .and. m.eq.3) then
     nterm=6
     lx(1)=4 ; ly(1)=0 ; lz(1)=1
     lx(2)=2 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=4 ; lz(3)=1
     lx(4)=2 ; ly(4)=0 ; lz(4)=3
     lx(5)=0 ; ly(5)=2 ; lz(5)=3
     lx(6)=0 ; ly(6)=0 ; lz(6)=5
     fac_arr(1)=0.06795295885835007261827187d0
     fac_arr(2)=0.1359059177167001452365437d0
     fac_arr(3)=0.06795295885835007261827187d0
     fac_arr(4)=0.1359059177167001452365437d0
     fac_arr(5)=0.1359059177167001452365437d0
     fac_arr(6)=0.06795295885835007261827187d0
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.1) then
     nterm=1
     lx(1)=0 ; ly(1)=1 ; lz(1)=1
     fac_arr(1)=1.414213562373095048801689d0
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.2) then
     nterm=1
     lx(1)=1 ; ly(1)=0 ; lz(1)=1
     fac_arr(1)=1.414213562373095048801689d0
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.3) then
     nterm=1
     lx(1)=1 ; ly(1)=1 ; lz(1)=0
     fac_arr(1)=1.414213562373095048801689d0
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.4) then
     nterm=2
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     fac_arr(1)=0.7071067811865475244008444d0
     fac_arr(2)=-0.7071067811865475244008444d0
  else if (l.eq.3 .and. i.eq.1 .and. m.eq.5) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=-0.4082482904638630163662140d0
     fac_arr(2)=-0.4082482904638630163662140d0
     fac_arr(3)=0.8164965809277260327324280d0
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.1) then
     nterm=3
     lx(1)=2 ; ly(1)=1 ; lz(1)=1
     lx(2)=0 ; ly(2)=3 ; lz(2)=1
     lx(3)=0 ; ly(3)=1 ; lz(3)=3
     fac_arr(1)=0.3563483225498991795794046d0
     fac_arr(2)=0.3563483225498991795794046d0
     fac_arr(3)=0.3563483225498991795794046d0
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.2) then
     nterm=3
     lx(1)=3 ; ly(1)=0 ; lz(1)=1
     lx(2)=1 ; ly(2)=2 ; lz(2)=1
     lx(3)=1 ; ly(3)=0 ; lz(3)=3
     fac_arr(1)=0.3563483225498991795794046d0
     fac_arr(2)=0.3563483225498991795794046d0
     fac_arr(3)=0.3563483225498991795794046d0
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.3) then
     nterm=3
     lx(1)=3 ; ly(1)=1 ; lz(1)=0
     lx(2)=1 ; ly(2)=3 ; lz(2)=0
     lx(3)=1 ; ly(3)=1 ; lz(3)=2
     fac_arr(1)=0.3563483225498991795794046d0
     fac_arr(2)=0.3563483225498991795794046d0
     fac_arr(3)=0.3563483225498991795794046d0
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.4) then
     nterm=4
     lx(1)=4 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=4 ; lz(2)=0
     lx(3)=2 ; ly(3)=0 ; lz(3)=2
     lx(4)=0 ; ly(4)=2 ; lz(4)=2
     fac_arr(1)=0.1781741612749495897897023d0
     fac_arr(2)=-0.1781741612749495897897023d0
     fac_arr(3)=0.1781741612749495897897023d0
     fac_arr(4)=-0.1781741612749495897897023d0
  else if (l.eq.3 .and. i.eq.2 .and. m.eq.5) then
     nterm=6
     lx(1)=4 ; ly(1)=0 ; lz(1)=0
     lx(2)=2 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=4 ; lz(3)=0
     lx(4)=2 ; ly(4)=0 ; lz(4)=2
     lx(5)=0 ; ly(5)=2 ; lz(5)=2
     lx(6)=0 ; ly(6)=0 ; lz(6)=4
     fac_arr(1)=-0.1028688999747279401740630d0
     fac_arr(2)=-0.2057377999494558803481260d0
     fac_arr(3)=-0.1028688999747279401740630d0
     fac_arr(4)=0.1028688999747279401740630d0
     fac_arr(5)=0.1028688999747279401740630d0
     fac_arr(6)=0.2057377999494558803481260d0
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.1) then
     nterm=6
     lx(1)=4 ; ly(1)=1 ; lz(1)=1
     lx(2)=2 ; ly(2)=3 ; lz(2)=1
     lx(3)=0 ; ly(3)=5 ; lz(3)=1
     lx(4)=2 ; ly(4)=1 ; lz(4)=3
     lx(5)=0 ; ly(5)=3 ; lz(5)=3
     lx(6)=0 ; ly(6)=1 ; lz(6)=5
     fac_arr(1)=0.05959868750235655989526993d0
     fac_arr(2)=0.1191973750047131197905399d0
     fac_arr(3)=0.05959868750235655989526993d0
     fac_arr(4)=0.1191973750047131197905399d0
     fac_arr(5)=0.1191973750047131197905399d0
     fac_arr(6)=0.05959868750235655989526993d0
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.2) then
     nterm=6
     lx(1)=5 ; ly(1)=0 ; lz(1)=1
     lx(2)=3 ; ly(2)=2 ; lz(2)=1
     lx(3)=1 ; ly(3)=4 ; lz(3)=1
     lx(4)=3 ; ly(4)=0 ; lz(4)=3
     lx(5)=1 ; ly(5)=2 ; lz(5)=3
     lx(6)=1 ; ly(6)=0 ; lz(6)=5
     fac_arr(1)=0.05959868750235655989526993d0
     fac_arr(2)=0.1191973750047131197905399d0
     fac_arr(3)=0.05959868750235655989526993d0
     fac_arr(4)=0.1191973750047131197905399d0
     fac_arr(5)=0.1191973750047131197905399d0
     fac_arr(6)=0.05959868750235655989526993d0
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.3) then
     nterm=6
     lx(1)=5 ; ly(1)=1 ; lz(1)=0
     lx(2)=3 ; ly(2)=3 ; lz(2)=0
     lx(3)=1 ; ly(3)=5 ; lz(3)=0
     lx(4)=3 ; ly(4)=1 ; lz(4)=2
     lx(5)=1 ; ly(5)=3 ; lz(5)=2
     lx(6)=1 ; ly(6)=1 ; lz(6)=4
     fac_arr(1)=0.05959868750235655989526993d0
     fac_arr(2)=0.1191973750047131197905399d0
     fac_arr(3)=0.05959868750235655989526993d0
     fac_arr(4)=0.1191973750047131197905399d0
     fac_arr(5)=0.1191973750047131197905399d0
     fac_arr(6)=0.05959868750235655989526993d0
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.4) then
     nterm=8
     lx(1)=6 ; ly(1)=0 ; lz(1)=0
     lx(2)=4 ; ly(2)=2 ; lz(2)=0
     lx(3)=2 ; ly(3)=4 ; lz(3)=0
     lx(4)=0 ; ly(4)=6 ; lz(4)=0
     lx(5)=4 ; ly(5)=0 ; lz(5)=2
     lx(6)=0 ; ly(6)=4 ; lz(6)=2
     lx(7)=2 ; ly(7)=0 ; lz(7)=4
     lx(8)=0 ; ly(8)=2 ; lz(8)=4
     fac_arr(1)=0.02979934375117827994763496d0
     fac_arr(2)=0.02979934375117827994763496d0
     fac_arr(3)=-0.02979934375117827994763496d0
     fac_arr(4)=-0.02979934375117827994763496d0
     fac_arr(5)=0.05959868750235655989526993d0
     fac_arr(6)=-0.05959868750235655989526993d0
     fac_arr(7)=0.02979934375117827994763496d0
     fac_arr(8)=-0.02979934375117827994763496d0
  else if (l.eq.3 .and. i.eq.3 .and. m.eq.5) then
     nterm=7
     lx(1)=6 ; ly(1)=0 ; lz(1)=0
     lx(2)=4 ; ly(2)=2 ; lz(2)=0
     lx(3)=2 ; ly(3)=4 ; lz(3)=0
     lx(4)=0 ; ly(4)=6 ; lz(4)=0
     lx(5)=2 ; ly(5)=0 ; lz(5)=4
     lx(6)=0 ; ly(6)=2 ; lz(6)=4
     lx(7)=0 ; ly(7)=0 ; lz(7)=6
     fac_arr(1)=-0.01720465913641697233541246d0
     fac_arr(2)=-0.05161397740925091700623738d0
     fac_arr(3)=-0.05161397740925091700623738d0
     fac_arr(4)=-0.01720465913641697233541246d0
     fac_arr(5)=0.05161397740925091700623738d0
     fac_arr(6)=0.05161397740925091700623738d0
     fac_arr(7)=0.03440931827283394467082492d0
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.1) then
     nterm=3
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=0.3162277660168379331998894d0
     fac_arr(2)=0.3162277660168379331998894d0
     fac_arr(3)=-1.264911064067351732799557d0
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.2) then
     nterm=3
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=1 ; lz(3)=2
     fac_arr(1)=0.3162277660168379331998894d0
     fac_arr(2)=0.3162277660168379331998894d0
     fac_arr(3)=-1.264911064067351732799557d0
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.3) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=0 ; lz(3)=3
     fac_arr(1)=0.7745966692414833770358531d0
     fac_arr(2)=0.7745966692414833770358531d0
     fac_arr(3)=-0.5163977794943222513572354d0
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.4) then
     nterm=2
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     fac_arr(1)=0.4082482904638630163662140d0
     fac_arr(2)=-1.224744871391589049098642d0
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.5) then
     nterm=2
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     fac_arr(1)=-1.224744871391589049098642d0
     fac_arr(2)=0.4082482904638630163662140d0
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.6) then
     nterm=2
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     fac_arr(1)=1.000000000000000000000000d0
     fac_arr(2)=-1.000000000000000000000000d0
  else if (l.eq.4 .and. i.eq.1 .and. m.eq.7) then
     nterm=1
     lx(1)=1 ; ly(1)=1 ; lz(1)=1
     fac_arr(1)=2.000000000000000000000000d0
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.1) then
     nterm=6
     lx(1)=5 ; ly(1)=0 ; lz(1)=0
     lx(2)=3 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=4 ; lz(3)=0
     lx(4)=3 ; ly(4)=0 ; lz(4)=2
     lx(5)=1 ; ly(5)=2 ; lz(5)=2
     lx(6)=1 ; ly(6)=0 ; lz(6)=4
     fac_arr(1)=0.06356417261637282102978506d0
     fac_arr(2)=0.1271283452327456420595701d0
     fac_arr(3)=0.06356417261637282102978506d0
     fac_arr(4)=-0.1906925178491184630893552d0
     fac_arr(5)=-0.1906925178491184630893552d0
     fac_arr(6)=-0.2542566904654912841191402d0
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.2) then
     nterm=6
     lx(1)=4 ; ly(1)=1 ; lz(1)=0
     lx(2)=2 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=5 ; lz(3)=0
     lx(4)=2 ; ly(4)=1 ; lz(4)=2
     lx(5)=0 ; ly(5)=3 ; lz(5)=2
     lx(6)=0 ; ly(6)=1 ; lz(6)=4
     fac_arr(1)=0.06356417261637282102978506d0
     fac_arr(2)=0.1271283452327456420595701d0
     fac_arr(3)=0.06356417261637282102978506d0
     fac_arr(4)=-0.1906925178491184630893552d0
     fac_arr(5)=-0.1906925178491184630893552d0
     fac_arr(6)=-0.2542566904654912841191402d0
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.3) then
     nterm=6
     lx(1)=4 ; ly(1)=0 ; lz(1)=1
     lx(2)=2 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=4 ; lz(3)=1
     lx(4)=2 ; ly(4)=0 ; lz(4)=3
     lx(5)=0 ; ly(5)=2 ; lz(5)=3
     lx(6)=0 ; ly(6)=0 ; lz(6)=5
     fac_arr(1)=0.1556997888323045941832351d0
     fac_arr(2)=0.3113995776646091883664703d0
     fac_arr(3)=0.1556997888323045941832351d0
     fac_arr(4)=0.05189992961076819806107838d0
     fac_arr(5)=0.05189992961076819806107838d0
     fac_arr(6)=-0.1037998592215363961221568d0
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.4) then
     nterm=5
     lx(1)=5 ; ly(1)=0 ; lz(1)=0
     lx(2)=3 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=4 ; lz(3)=0
     lx(4)=3 ; ly(4)=0 ; lz(4)=2
     lx(5)=1 ; ly(5)=2 ; lz(5)=2
     fac_arr(1)=0.08206099398622182182282711d0
     fac_arr(2)=-0.1641219879724436436456542d0
     fac_arr(3)=-0.2461829819586654654684813d0
     fac_arr(4)=0.08206099398622182182282711d0
     fac_arr(5)=-0.2461829819586654654684813d0
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.5) then
     nterm=5
     lx(1)=4 ; ly(1)=1 ; lz(1)=0
     lx(2)=2 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=5 ; lz(3)=0
     lx(4)=2 ; ly(4)=1 ; lz(4)=2
     lx(5)=0 ; ly(5)=3 ; lz(5)=2
     fac_arr(1)=-0.2461829819586654654684813d0
     fac_arr(2)=-0.1641219879724436436456542d0
     fac_arr(3)=0.08206099398622182182282711d0
     fac_arr(4)=-0.2461829819586654654684813d0
     fac_arr(5)=0.08206099398622182182282711d0
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.6) then
     nterm=4
     lx(1)=4 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=4 ; lz(2)=1
     lx(3)=2 ; ly(3)=0 ; lz(3)=3
     lx(4)=0 ; ly(4)=2 ; lz(4)=3
     fac_arr(1)=0.2010075630518424150978747d0
     fac_arr(2)=-0.2010075630518424150978747d0
     fac_arr(3)=0.2010075630518424150978747d0
     fac_arr(4)=-0.2010075630518424150978747d0
  else if (l.eq.4 .and. i.eq.2 .and. m.eq.7) then
     nterm=3
     lx(1)=3 ; ly(1)=1 ; lz(1)=1
     lx(2)=1 ; ly(2)=3 ; lz(2)=1
     lx(3)=1 ; ly(3)=1 ; lz(3)=3
     fac_arr(1)=0.4020151261036848301957494d0
     fac_arr(2)=0.4020151261036848301957494d0
     fac_arr(3)=0.4020151261036848301957494d0
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.1) then
     nterm=10
     lx(1)=7 ; ly(1)=0 ; lz(1)=0
     lx(2)=5 ; ly(2)=2 ; lz(2)=0
     lx(3)=3 ; ly(3)=4 ; lz(3)=0
     lx(4)=1 ; ly(4)=6 ; lz(4)=0
     lx(5)=5 ; ly(5)=0 ; lz(5)=2
     lx(6)=3 ; ly(6)=2 ; lz(6)=2
     lx(7)=1 ; ly(7)=4 ; lz(7)=2
     lx(8)=3 ; ly(8)=0 ; lz(8)=4
     lx(9)=1 ; ly(9)=2 ; lz(9)=4
     lx(10)=1 ; ly(10)=0 ; lz(10)=6
     fac_arr(1)=0.009103849893318918298413687d0
     fac_arr(2)=0.02731154967995675489524106d0
     fac_arr(3)=0.02731154967995675489524106d0
     fac_arr(4)=0.009103849893318918298413687d0
     fac_arr(5)=-0.01820769978663783659682737d0
     fac_arr(6)=-0.03641539957327567319365475d0
     fac_arr(7)=-0.01820769978663783659682737d0
     fac_arr(8)=-0.06372694925323242808889581d0
     fac_arr(9)=-0.06372694925323242808889581d0
     fac_arr(10)=-0.03641539957327567319365475d0
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.2) then
     nterm=10
     lx(1)=6 ; ly(1)=1 ; lz(1)=0
     lx(2)=4 ; ly(2)=3 ; lz(2)=0
     lx(3)=2 ; ly(3)=5 ; lz(3)=0
     lx(4)=0 ; ly(4)=7 ; lz(4)=0
     lx(5)=4 ; ly(5)=1 ; lz(5)=2
     lx(6)=2 ; ly(6)=3 ; lz(6)=2
     lx(7)=0 ; ly(7)=5 ; lz(7)=2
     lx(8)=2 ; ly(8)=1 ; lz(8)=4
     lx(9)=0 ; ly(9)=3 ; lz(9)=4
     lx(10)=0 ; ly(10)=1 ; lz(10)=6
     fac_arr(1)=0.009103849893318918298413687d0
     fac_arr(2)=0.02731154967995675489524106d0
     fac_arr(3)=0.02731154967995675489524106d0
     fac_arr(4)=0.009103849893318918298413687d0
     fac_arr(5)=-0.01820769978663783659682737d0
     fac_arr(6)=-0.03641539957327567319365475d0
     fac_arr(7)=-0.01820769978663783659682737d0
     fac_arr(8)=-0.06372694925323242808889581d0
     fac_arr(9)=-0.06372694925323242808889581d0
     fac_arr(10)=-0.03641539957327567319365475d0
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.3) then
     nterm=10
     lx(1)=6 ; ly(1)=0 ; lz(1)=1
     lx(2)=4 ; ly(2)=2 ; lz(2)=1
     lx(3)=2 ; ly(3)=4 ; lz(3)=1
     lx(4)=0 ; ly(4)=6 ; lz(4)=1
     lx(5)=4 ; ly(5)=0 ; lz(5)=3
     lx(6)=2 ; ly(6)=2 ; lz(6)=3
     lx(7)=0 ; ly(7)=4 ; lz(7)=3
     lx(8)=2 ; ly(8)=0 ; lz(8)=5
     lx(9)=0 ; ly(9)=2 ; lz(9)=5
     lx(10)=0 ; ly(10)=0 ; lz(10)=7
     fac_arr(1)=0.02229978693352242055222348d0
     fac_arr(2)=0.06689936080056726165667044d0
     fac_arr(3)=0.06689936080056726165667044d0
     fac_arr(4)=0.02229978693352242055222348d0
     fac_arr(5)=0.02973304924469656073629797d0
     fac_arr(6)=0.05946609848939312147259594d0
     fac_arr(7)=0.02973304924469656073629797d0
     fac_arr(8)=-0.007433262311174140184074493d0
     fac_arr(9)=-0.007433262311174140184074493d0
     fac_arr(10)=-0.01486652462234828036814899d0
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.4) then
     nterm=9
     lx(1)=7 ; ly(1)=0 ; lz(1)=0
     lx(2)=5 ; ly(2)=2 ; lz(2)=0
     lx(3)=3 ; ly(3)=4 ; lz(3)=0
     lx(4)=1 ; ly(4)=6 ; lz(4)=0
     lx(5)=5 ; ly(5)=0 ; lz(5)=2
     lx(6)=3 ; ly(6)=2 ; lz(6)=2
     lx(7)=1 ; ly(7)=4 ; lz(7)=2
     lx(8)=3 ; ly(8)=0 ; lz(8)=4
     lx(9)=1 ; ly(9)=2 ; lz(9)=4
     fac_arr(1)=0.01175301967439877980816756d0
     fac_arr(2)=-0.01175301967439877980816756d0
     fac_arr(3)=-0.05876509837199389904083778d0
     fac_arr(4)=-0.03525905902319633942450267d0
     fac_arr(5)=0.02350603934879755961633511d0
     fac_arr(6)=-0.04701207869759511923267022d0
     fac_arr(7)=-0.07051811804639267884900533d0
     fac_arr(8)=0.01175301967439877980816756d0
     fac_arr(9)=-0.03525905902319633942450267d0
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.5) then
     nterm=9
     lx(1)=6 ; ly(1)=1 ; lz(1)=0
     lx(2)=4 ; ly(2)=3 ; lz(2)=0
     lx(3)=2 ; ly(3)=5 ; lz(3)=0
     lx(4)=0 ; ly(4)=7 ; lz(4)=0
     lx(5)=4 ; ly(5)=1 ; lz(5)=2
     lx(6)=2 ; ly(6)=3 ; lz(6)=2
     lx(7)=0 ; ly(7)=5 ; lz(7)=2
     lx(8)=2 ; ly(8)=1 ; lz(8)=4
     lx(9)=0 ; ly(9)=3 ; lz(9)=4
     fac_arr(1)=-0.03525905902319633942450267d0
     fac_arr(2)=-0.05876509837199389904083778d0
     fac_arr(3)=-0.01175301967439877980816756d0
     fac_arr(4)=0.01175301967439877980816756d0
     fac_arr(5)=-0.07051811804639267884900533d0
     fac_arr(6)=-0.04701207869759511923267022d0
     fac_arr(7)=0.02350603934879755961633511d0
     fac_arr(8)=-0.03525905902319633942450267d0
     fac_arr(9)=0.01175301967439877980816756d0
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.6) then
     nterm=8
     lx(1)=6 ; ly(1)=0 ; lz(1)=1
     lx(2)=4 ; ly(2)=2 ; lz(2)=1
     lx(3)=2 ; ly(3)=4 ; lz(3)=1
     lx(4)=0 ; ly(4)=6 ; lz(4)=1
     lx(5)=4 ; ly(5)=0 ; lz(5)=3
     lx(6)=0 ; ly(6)=4 ; lz(6)=3
     lx(7)=2 ; ly(7)=0 ; lz(7)=5
     lx(8)=0 ; ly(8)=2 ; lz(8)=5
     fac_arr(1)=0.02878890113916869875409405d0
     fac_arr(2)=0.02878890113916869875409405d0
     fac_arr(3)=-0.02878890113916869875409405d0
     fac_arr(4)=-0.02878890113916869875409405d0
     fac_arr(5)=0.05757780227833739750818811d0
     fac_arr(6)=-0.05757780227833739750818811d0
     fac_arr(7)=0.02878890113916869875409405d0
     fac_arr(8)=-0.02878890113916869875409405d0
  else if (l.eq.4 .and. i.eq.3 .and. m.eq.7) then
     nterm=6
     lx(1)=5 ; ly(1)=1 ; lz(1)=1
     lx(2)=3 ; ly(2)=3 ; lz(2)=1
     lx(3)=1 ; ly(3)=5 ; lz(3)=1
     lx(4)=3 ; ly(4)=1 ; lz(4)=3
     lx(5)=1 ; ly(5)=3 ; lz(5)=3
     lx(6)=1 ; ly(6)=1 ; lz(6)=5
     fac_arr(1)=0.05757780227833739750818811d0
     fac_arr(2)=0.1151556045566747950163762d0
     fac_arr(3)=0.05757780227833739750818811d0
     fac_arr(4)=0.1151556045566747950163762d0
     fac_arr(5)=0.1151556045566747950163762d0
     fac_arr(6)=0.05757780227833739750818811d0

  else
     stop 'PSP format error'
  endif
  
END SUBROUTINE calc_coeff_proj
