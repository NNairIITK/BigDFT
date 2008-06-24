subroutine readAtomicOrbitals(iproc,ngx,xp,psiat,occupat,ng,nl,nzatom,nelpsp,&
     & psppar,npspcode,norbe,norbsc,atomnames,ntypes,iatype,natpol,iasctype,nat,natsc,nspin,&
     & scorb,norbsc_arr)
  use module_base
  implicit none
  integer, intent(in) :: ngx,iproc,ntypes,nspin,nat,natsc
  character(len=20), dimension(ntypes), intent(in) :: atomnames
  integer, dimension(ntypes), intent(in) :: nzatom,nelpsp,npspcode,iasctype
  integer, dimension(nat), intent(in) :: iatype,natpol
  real(gp), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  integer, intent(out) :: norbe,norbsc
  logical, dimension(4,2,natsc), intent(out) :: scorb
  integer, dimension(ntypes), intent(out) :: ng
  integer, dimension(4,ntypes), intent(out) :: nl
  integer, dimension(natsc+1,nspin), intent(out) :: norbsc_arr
  real(gp), dimension(ngx,ntypes), intent(out) :: xp
  real(gp), dimension(5,ntypes), intent(out) :: occupat
  real(gp), dimension(ngx,5,ntypes), intent(out) :: psiat
  !local variables
  character(len = 20) :: pspatomname
  logical :: exists,found
  integer :: ity,i,j,l,ifile,ng_fake,ierror,iatsc,iat,ipow,lsc,inorbsc,iorbsc_count,niasc,nlsc
  real(kind=8) :: sccode

  ! Read the data file.
  nl(1:4,1:ntypes) = 0
  ng(1:ntypes) = 0
  xp(1:ngx,1:ntypes) = 0.0_gp
  psiat(1:ngx,1:5,1:ntypes) = 0.0_gp
  occupat(1:5,1:ntypes)= 0.0_gp

  ! Test if the file 'inguess.dat exists (no more used)
  inquire(file='inguess.dat',exist=exists)
  !exists=.false.
  if (exists) then
     open(unit=24,file='inguess.dat',form='formatted',action='read',status='old')
  end if

  loop_assign: do ity=1,ntypes

     if (exists) then
        rewind(24)
     end if
     found = .false.

     loop_find: do
        if (.not.exists) then
           !          The file 'inguess.dat' does not exist: automatic generation
           exit loop_find
        end if
        read(24,'(a)',iostat=ierror) pspatomname
        if (ierror /= 0) then
           !          Read error or end of file
           exit loop_find
        end if

        if (pspatomname .eq. atomnames(ity)) then
           if (iproc.eq.0) then
              write(*,'(1x,a,a,a)') 'input wavefunction data for atom ',trim(atomnames(ity)),&
                   ' found'
           end if
           found = .true.
           read(24,*) nl(1,ity),(occupat(i,ity),i=1,nl(1,ity)),  &
              nl(2,ity),(occupat(i,ity),i=1+nl(1,ity),nl(2,ity)+nl(1,ity)) ,&
              nl(3,ity),(occupat(i,ity),i=1+nl(2,ity)+nl(1,ity),nl(3,ity)+nl(2,ity)+nl(1,ity)),&
              nl(4,ity),(occupat(i,ity),&
              i=1+nl(3,ity)+nl(2,ity)+nl(1,ity),nl(4,ity)+nl(3,ity)+nl(2,ity)+nl(1,ity))
           !print *,nl(:,ity),occupat(:,ity)
           if (nl(1,ity)+nl(2,ity)+nl(3,ity)+nl(4,ity).gt.5) then
              print *,'error: number of valence orbitals too big'
              print *,nl(:,ity),occupat(:,ity)
              stop
           end if
           read(24,*) ng(ity)
           !print *, pspatomnames(ity),(nl(l,ity),l=1,4),ng(ity),ngx,npsp
           if (ng(ity).gt.ngx) stop 'enlarge ngx'
           !read(24,'(30(e12.5))') (xp(i,ity)  ,i=1,ng(ity))
           read(24,*) (xp(i,ity)  ,i=1,ng(ity))
           do i=1,ng(ity) 
              read(24,*) (psiat(i,j,ity),j=1,nl(1,ity)+nl(2,ity)+nl(3,ity)+nl(4,ity))
           enddo

           exit loop_find

        else

           read(24,*)
           read(24,*)ng_fake
           read(24,*) 
           do i=1,ng_fake
              read(24,*)
           enddo

        end if

     enddo loop_find

     if (.not.found) then

        if (iproc.eq.0) then
           write(*,'(1x,a,a6,a)',advance='no')&
                'Input wavefunction data for atom ',trim(atomnames(ity)),&
                ' NOT found, automatic generation...'
        end if

        !the default value for the gaussians is chosen to be 21
        ng(ity)=21

        call iguess_generator(iproc,nzatom(ity),nelpsp(ity),psppar(0,0,ity),npspcode(ity),&
             ng(ity)-1,nl(1,ity),5,occupat(1:5,ity),xp(1:ng(ity),ity),psiat(1:ng(ity),1:5,ity))

!!$        !values obtained from the input guess generator in iguess.dat format
!!$        !write these values on a file in the case of the HGH-K pseudo, for check
!!$        if (iproc .eq. 0 .and. npspcode(ity)==10) then
!!$           open(unit=12,file='inguess.new',status='unknown')
!!$
!!$           !write(*,*)' --------COPY THESE VALUES INSIDE inguess.dat--------'
!!$           write(12,*)trim(atomnames(ity))//' (remove _lda)'
!!$           write(12,*)nl(1,ity),(occupat(i,ity),i=1,nl(1,ity)),&
!!$                nl(2,ity),(occupat(i+nl(1,ity),ity),i=1,nl(2,ity)),&
!!$                nl(3,ity),(occupat(i+nl(1,ity)+nl(2,ity),ity),i=1,nl(3,ity)),&
!!$                nl(4,ity),(occupat(i+nl(1,ity)+nl(2,ity)+nl(3,ity),ity),i=1,nl(4,ity))
!!$           write(12,*)ng(ity)
!!$           write(12,'(30(e12.5))')xp(1:ng(ity),ity)
!!$           do j=1,ng(ity)
!!$              write(12,*)(psiat(j,i,ity),i=1,nl(1,ity)+nl(2,ity)+nl(3,ity)+nl(4,ity))
!!$           end do
!!$           !print *,' --------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^--------'
!!$        end if
        if (iproc.eq.0) write(*,'(1x,a)')'done.'

     end if

  end do loop_assign

  close(unit=24)

  ! number of orbitals, total and semicore
  !HERE WE MUST INSERT THE SPIN DEPENDENCE
  norbe=0
  norbsc=0
  iatsc=0
  scorb(:,:,:)=.false.
  do iat=1,nat
     ity=iatype(iat)
     norbe=norbe+nl(1,ity)+3*nl(2,ity)+5*nl(3,ity)+7*nl(4,ity)
     if (iasctype(ity)/=0) then !the atom has some semicore orbitals
        iatsc=iatsc+1
        niasc=iasctype(ity)
        !count the semicore orbitals for this atom
        iorbsc_count=0
        do lsc=4,1,-1
           nlsc=niasc/4**(lsc-1)
           iorbsc_count=iorbsc_count+nlsc*(2*lsc-1)
           do i=1,nlsc
              scorb(lsc,i,iatsc)=.true.
           end do
           niasc=niasc-nlsc*4**(lsc-1)
        end do
!!$        !here we should assign the new way of determining semicore orbitals 
!!$        iatsc=iatsc+1
!!$        sccode=real(iasctype(ity),kind=8)
!!$        inorbsc=ceiling(dlog(sccode)/dlog(10.d0))
!!$        if (sccode==1.d0) inorbsc=1
!!$        ipow=inorbsc-1
!!$        scorb(:,:,iatsc)=.false.
!!$        !count the semicore orbitals for this atom
!!$        iorbsc_count=0
!!$        do i=1,inorbsc
!!$           lsc=floor(sccode/10.d0**ipow)
!!$           iorbsc_count=iorbsc_count+(2*lsc-1)
!!$           scorb(lsc,1,iatsc)=.true.
!!$           sccode=sccode-real(lsc,kind=8)*10.d0**ipow
!!$           ipow=ipow-1
!!$           !print *,iasctype(ity),inorbsc,lsc
!!$        end do

        norbsc_arr(iatsc,1)=iorbsc_count
        norbsc=norbsc+iorbsc_count
        !if (iproc == 0) print *,iasctype(ity),iorbsc_count,norbsc,scorb(:,:,iatsc)
     end if
  end do

  !orbitals which are non semicore
  norbsc_arr(natsc+1,1)=norbe-norbsc

  !duplicate the values in the case of spin-polarization
  if (nspin == 2) norbsc_arr(:,2)=norbsc_arr(:,1)

  if (iproc ==0) then
     write(*,'(1x,a,i0,a)')'Generating ',nspin*norbe,' Atomic Input Orbitals'
     if (norbsc /=0)   write(*,'(1x,a,i0,a)')'  of which ',nspin*norbsc,' are semicore orbitals'
  end if

END SUBROUTINE readAtomicOrbitals

subroutine createAtomicOrbitals(geocode,iproc,nproc,atomnames,&
     & nat,rxyz,norbe,norbep,norbsc,occupe,occupat,ngx,xp,psiat,ng,nl,&
     & nvctr_c,nvctr_f,n1,n2,n3,hx,hy,hz,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nseg_c,nseg_f,&
     & keyg,keyv,iatype,ntypes,iasctype,natsc,natpol,nspin,psi,eks,scorb)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: nat,norbe,norbep,ngx,iproc,nproc
  integer, intent(in) :: nvctr_c,nvctr_f,n1,n2,n3,nseg_c,nseg_f
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,ntypes
  integer, intent(in) :: norbsc,natsc,nspin
  real(gp), intent(in) :: hx,hy,hz
  character(len=20), dimension(ntypes), intent(in) :: atomnames
  logical, dimension(4,2,natsc), intent(in) :: scorb
  integer, dimension(ntypes), intent(in) :: iasctype
  integer, dimension(nat), intent(in) :: iatype,natpol
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  real(gp), dimension(3,nat), intent(in) :: rxyz
  integer, dimension(ntypes), intent(inout) :: ng
  integer, dimension(4,ntypes), intent(inout) :: nl
  real(gp), dimension(ngx,ntypes), intent(inout) :: xp
  real(gp), dimension(5,ntypes), intent(inout) :: occupat
  real(gp), dimension(ngx,5,ntypes), intent(inout) :: psiat
  real(gp), intent(out) :: eks
  real(gp), dimension(nspin*norbe), intent(out) :: occupe
  real(wp), dimension(nvctr_c+7*nvctr_f,norbep), intent(out) :: psi
  !local variables
  character(len=*), parameter :: subname= 'createAtomicOrbitals'
  integer, parameter :: nterm_max=3,noccmax=2,nlmax=4
  logical :: myorbital,polarised
  integer :: iatsc,i_all,i_stat,ispin,ipolres,ipolorb,ichg
  integer :: iorb,jorb,iat,ity,i,ictot,inl,l,m,nctot,nterm
  real(wp) :: scprw
  real(dp) :: scpr
  real(gp) :: rx,ry,rz,ek,occshell
  logical, dimension(nlmax,noccmax) :: semicore
  integer, dimension(2) :: iorbsc,iorbv
  integer, dimension(nterm_max) :: lx,ly,lz
  real(gp), dimension(nterm_max) :: fac_arr
  real(gp), dimension(:), allocatable :: psiatn

  allocate(psiatn(ngx+ndebug),stat=i_stat)
  call memocc(i_stat,psiatn,'psiatn',subname)

  eks=0.0_gp
  iorb=0
  iatsc=0

  !initialise the orbital counters
  iorbsc(1)=0
  iorbv(1)=norbsc
  !used in case of spin-polarisation, ignored otherwise
  iorbsc(2)=norbe
  iorbv(2)=norbsc+norbe


  if (iproc ==0) then
     write(*,'(1x,a)',advance='no')'Calculating AIO wavefunctions...'
  end if

  do iat=1,nat

     rx=rxyz(1,iat)
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)

     ity=iatype(iat)

     !here we can evaluate whether the atom has semicore orbitals and with
     !which value(s) of l
     if (iasctype(ity)/=0) then !the atom has some semicore orbitals
        iatsc=iatsc+1
        semicore(:,:)=scorb(:,:,iatsc)
     else
        semicore(:,:)=.false.
     end if

     !calculate the atomic input orbitals
     ictot=0
     nctot=nl(1,ity)+nl(2,ity)+nl(3,ity)+nl(4,ity)
     if (iorbsc(1)+nctot > norbe .and. iorbv(1)+nctot > norbe) then
        print *,'transgpw occupe',nl(:,ity),norbe
        stop
     end if
     polarised=.false.
     call charge_and_spol(natpol(iat),ichg,ipolres)
     do l=1,4
        do inl=1,nl(l,ity)
           ictot=ictot+1
           !contribution to the kinetic energy given by the electrons in this shell
           call atomkin(l-1,ng(ity),xp(1,ity),psiat(1,ictot,ity),psiatn,ek)
           eks=eks+ek*occupat(ictot,ity)
           if (nint(occupat(ictot,ity)) /=  2*(2*l-1) ) then
              !this is a polarisable orbital
!!$              if (polarised) then
!!$                 if (iproc == 0) then
!!$                    write(*,'(1x,a)')&
!!$                         'ERROR: only one polarisable orbital is allowed, check electronic configuration'
!!$                    stop
!!$                 end if
!!$              end if
              polarised=.true.
              !assuming that the control of the allowed polarisation is already done
              ipolorb=min(ipolres,int(occupat(ictot,ity)))
              ipolres=ipolres-ipolorb
              !this check can be inserted also elsewhere
              if (ipolres < 0) then
                 if(iproc==0) write(*,'(1x,4(a,i0))')&
                      'Too high polarisation for atom number= ',iat,&
                      ' Inserted=',modulo(natpol(iat),1000)-100,' Assigned=',ipolorb,&
                      ' the maximum is=',nint(occupat(ictot,ity))
                 stop
              end if

           else
              !check for odd values of the occupation number
              if (mod(nint(occupat(ictot,ity)),2) /= 0) then
                 if (iproc == 0) write(*,'(1x,a)')&
                      'The occupation number in the case of closed shells must be even'
                 stop
              end if
           end if
           do ispin=1,nspin
              !the order of the orbitals (iorb,jorb) must put in the beginning
              !the semicore orbitals
              if (semicore(l,inl)) then
                 !the orbital is semi-core
                 iorb=iorbsc(ispin)
                 !print *,'iproc, SEMICORE orbital, iat,l',iproc,iat,l
                 !the occupation number is divided by two in the case of spin polarisation
                 if (nspin==2) then
                    occshell=0.5_gp*occupat(ictot,ity)
                 else
                    occshell=occupat(ictot,ity)
                 end if
              else
                 !normal case, the orbital is a valence orbital
                 iorb=iorbv(ispin)
                 occshell=occupat(ictot,ity)                 
                 if (nspin==2) then
                    if (polarised) then
                       occshell=0.5_gp*(occshell+real(1-2*(ispin-1),gp)*ipolorb)
                    else
                       occshell=0.5_gp*occshell
                    end if
                 end if
              end if

              do m=1,2*l-1
                 iorb=iorb+1
                 jorb=iorb-iproc*norbep
                 occupe(iorb)=occshell/real(2*l-1,gp)
                 if (myorbital(iorb,nspin*norbe,iproc,nproc)) then
                    !this will calculate the proper spherical harmonics
                    call calc_coeff_inguess(l,m,nterm_max,nterm,lx,ly,lz,fac_arr)
                    call crtonewave(geocode,n1,n2,n3,ng(ity),nterm,lx,ly,lz,fac_arr,&
                         xp(1,ity),psiatn,rx,ry,rz,hx,hy,hz, & 
                         0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                         nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,&
                         keyg(1,nseg_c+1),keyv(nseg_c+1),&
                         psi(1,jorb),psi(nvctr_c+1,jorb))
                    !renormalise wavefunction in case of too crude box
                    call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
                    !in the periodic case the function is not always normalised
                    !write(*,'(1x,a24,a7,2(a3,i1),a16,i4,i4,1x,1pe14.7)')&
                    !     'ATOMIC INPUT ORBITAL for atom',trim(atomnames(ity)),&
                    !     'l=',l,'m=',m,'iorb,jorb,norm',iorb,jorb,scpr 
                    scprw=real(1.0_dp/sqrt(scpr),wp)
                    call wscal(nvctr_c,nvctr_f,scprw,psi(1,jorb),psi(nvctr_c+1,jorb))
                    !call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
                    !print *,'newnorm', scpr,occupe(iorb),occshell,ictot
                 endif
              end do
              if (semicore(l,inl)) then
                    !increase semicore orbitals
                    iorbsc(ispin)=iorb
              else
                 !increase valence orbitals
                 iorbv(ispin)=iorb
              end if
           end do
        end do
     end do
     
     if (ictot /= nctot) stop 'createAtomic orbitals: error (nctot)'

  end do
  if (iorbsc(1) /= norbsc) then
     stop 'createAtomic orbitals: error (iorbsc)'
  end if
  if (iorbv(1) /= norbe) stop 'createAtomic orbitals: error (iorbv)'
  if (iatsc /= natsc) stop 'createAtomic orbitals: error (iatsc)'

  if (nspin==2) then
     if (iorbsc(2) /= norbsc+norbe) stop 'createAtomic orbitals: error (iorbsc) nspin=2'
     if (iorbv(2) /= 2*norbe) stop 'createAtomic orbitals: error (iorbv) nspin=2'
  end if

  i_all=-product(shape(psiatn))*kind(psiatn)
  deallocate(psiatn,stat=i_stat)
  call memocc(i_stat,i_all,'psiatn',subname)

  if (iproc ==0) then
     write(*,'(1x,a)')'done.'
  end if

end subroutine createAtomicOrbitals


! calculates the kinetic energy of an atomic wavefunction expressed in Gaussians
! the output psiatn is a normalized version of psiat
subroutine atomkin(l,ng,xp,psiat,psiatn,ek)
  use module_base
  implicit none
  integer, intent(in) :: l,ng
  real(gp), dimension(ng), intent(in) :: xp,psiat
  real(gp), intent(out) :: ek
  real(gp), dimension(ng), intent(out) :: psiatn
  !local variables
  integer :: i,j
  real(gp) :: gml,tt,xpi,xpj,d,sxp,const,hij,sij

  !        gml=.5d0*gamma(.5d0+l)
  gml = 0.0_gp
  if (l.eq.0) then 
     gml=0.88622692545275801365_gp
  else if (l.eq.1) then 
     gml=0.44311346272637900682_gp
  else if (l.eq.2) then 
     gml=0.66467019408956851024_gp
  else if (l.eq.3) then 
     gml=1.6616754852239212756_gp
  else
     stop 'atomkin'
  endif

  ek=0.0_gp
  tt=0.0_gp
  do i=1,ng
     xpi=.5_gp/xp(i)**2
     do j=1,ng
        xpj=.5_gp/xp(j)**2
        d=xpi+xpj
        sxp=1.0_gp/d
        const=gml*sqrt(sxp)**(2*l+1)
        ! kinetic energy  matrix element hij
        hij=.5_gp*const*sxp**2* ( 3._gp*xpi*xpj +                  &
             real(l,gp)*(6._gp*xpi*xpj-xpi**2-xpj**2) -        &
             real(l**2,gp)*(xpi-xpj)**2  ) + .5_gp*real(l,gp)*(real(l,gp)+1._gp)*const
        sij=const*sxp*(real(l,gp)+.5_gp)
        ek=ek+hij*psiat(i)*psiat(j)
        tt=tt+sij*psiat(i)*psiat(j)
     enddo
  enddo

  if (abs(tt-1._gp).gt.1.e-2_gp) write(*,*) 'presumably wrong inguess data',l,tt
  ! energy expectation value
  ek=ek/tt
  !write(*,*) 'ek=',ek,tt,l,ng
  ! scale atomic wavefunction
  tt=sqrt(1._gp/tt)
!!$        if (l.eq.0) then  ! multiply with 1/sqrt(4*pi)
!!$        tt=tt*0.28209479177387814347_gp
!!$        else if (l.eq.1) then  ! multiply with sqrt(3/(4*pi))
!!$        tt=tt*0.48860251190291992159_gp
!!$        !decide the value of the normalization to be used
!!$        endif
  do i=1,ng
     psiatn(i)=psiat(i)*tt
  enddo

end subroutine atomkin


subroutine calc_coeff_inguess(l,m,nterm_max,nterm,lx,ly,lz,fac_arr)
  use module_base
  implicit none
  integer, intent(in) :: l,m,nterm_max
  integer, intent(out) :: nterm
  integer, dimension(nterm_max), intent(out) :: lx,ly,lz
  real(gp), dimension(nterm_max), intent(out) :: fac_arr

  if (l.eq.1 .and. m.eq.1) then
     nterm=1
     lx(1)=0 ; ly(1)=0 ; lz(1)=0
     fac_arr(1)=0.28209479177387814347_gp

  else if (l.eq.2  .and. m.eq.1) then
     nterm=1
     lx(1)=1 ; ly(1)=0 ; lz(1)=0
     fac_arr(1)=0.48860251190291992159_gp
  else if (l.eq.2  .and. m.eq.2) then
     nterm=1
     lx(1)=0 ; ly(1)=1 ; lz(1)=0
     fac_arr(1)=0.48860251190291992159_gp
  else if (l.eq.2  .and. m.eq.3) then
     nterm=1
     lx(1)=0 ; ly(1)=0 ; lz(1)=1
     fac_arr(1)=0.48860251190291992159_gp

  else if (l.eq.3  .and. m.eq.1) then
     nterm=1
     lx(1)=0 ; ly(1)=1 ; lz(1)=1
     fac_arr(1)=1.092548430592079_gp
  else if (l.eq.3  .and. m.eq.2) then
     nterm=1
     lx(1)=1 ; ly(1)=0 ; lz(1)=1
     fac_arr(1)=1.092548430592079_gp
  else if (l.eq.3  .and. m.eq.3) then
     nterm=1
     lx(1)=1 ; ly(1)=1 ; lz(1)=0
     fac_arr(1)=1.092548430592079_gp
  else if (l.eq.3  .and. m.eq.4) then
     nterm=2
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     fac_arr(1)=0.5462742152960396_gp
     fac_arr(2)=-0.5462742152960396_gp
  else if (l.eq.3  .and. m.eq.5) then 
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=-0.3153915652525201_gp
     fac_arr(2)=-0.3153915652525201_gp
     fac_arr(3)=2._gp*0.3153915652525201_gp

  else if (l.eq.4  .and. m.eq.1) then
     nterm=3
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=0.4570457994644658_gp
     fac_arr(2)=0.4570457994644658_gp
     fac_arr(3)=-4._gp*0.4570457994644658_gp
  else if (l.eq.4  .and. m.eq.2) then
     nterm=3
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=1 ; lz(3)=2
     fac_arr(1)=0.4570457994644658_gp
     fac_arr(2)=0.4570457994644658_gp
     fac_arr(3)=-4._gp*0.4570457994644658_gp
  else if (l.eq.4  .and. m.eq.3) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=0 ; lz(3)=3
     fac_arr(1)=3._gp*0.3731763325901154_gp
     fac_arr(2)=3._gp*0.3731763325901154_gp
     fac_arr(3)=-2._gp*0.3731763325901154_gp
  else if (l.eq.4  .and. m.eq.4) then
     nterm=2
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     fac_arr(1)=0.5900435899266436_gp
     fac_arr(2)=-3._gp*0.5900435899266436_gp
  else if (l.eq.4  .and. m.eq.5) then
     nterm=2
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     fac_arr(1)=-3._gp*0.5900435899266436_gp
     fac_arr(2)=0.5900435899266436_gp
  else if (l.eq.4  .and. m.eq.6) then
     nterm=2
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     fac_arr(1)=1.445305721320277_gp
     fac_arr(2)=-1.445305721320277_gp
  else if (l.eq.4  .and. m.eq.7) then
     nterm=1
     lx(1)=1 ; ly(1)=1 ; lz(1)=1
     fac_arr(1)=2.890611442640554_gp
  else
     print *,'l,m',l,m
     stop 'input guess format error'
  endif

END SUBROUTINE calc_coeff_inguess


subroutine iguess_generator(iproc,izatom,ielpsp,psppar,npspcode,ng,nl,nmax_occ,occupat,expo,psiat)
  use module_base
  implicit none
  integer, intent(in) :: iproc,izatom,ielpsp,ng,npspcode,nmax_occ
  real(gp), dimension(0:4,0:6), intent(in) :: psppar
  integer, dimension(4), intent(out) :: nl
  real(gp), dimension(ng+1), intent(out) :: expo
  real(gp), dimension(nmax_occ), intent(out) :: occupat
  real(gp), dimension(ng+1,nmax_occ), intent(out) :: psiat
  !local variables
  character(len=*), parameter :: subname='iguess_generator'
  character(len=27) :: string 
  character(len=2) :: symbol
  integer, parameter :: lmax=3,n_int=100,noccmax=2
  real(gp), parameter :: fact=4.0_gp
  integer, dimension(6,4) :: neleconf
  real(gp), dimension(3) :: gpot
  real(gp), dimension(6) :: ott
  real(gp), dimension(noccmax,lmax+1) :: occup,aeval,chrg,res
  real(gp), dimension(:), allocatable :: xp,alps
  real(gp), dimension(:,:), allocatable :: vh,hsep,ofdcoef
  real(gp), dimension(:,:,:), allocatable :: psi
  real(gp), dimension(:,:,:,:), allocatable :: rmt
  logical :: exists
  integer :: lpx,ncount,nsccode,mxpl,mxchg
  integer :: l,i,j,iocc,il,lwrite,i_all,i_stat
  real(gp) :: alpz,alpl,rcov,rprb,zion,rij,a,a0,a0in,tt,ehomo

  !filename = 'psppar.'//trim(atomname)

  lpx=0
  lpx_determination: do i=1,4
     if (psppar(i,0) == 0.0_gp) then
     exit lpx_determination
     else
        lpx=i-1
     end if
  end do lpx_determination

  allocate(alps(lpx+1+ndebug),stat=i_stat)
  call memocc(i_stat,alps,'alps',subname)
  allocate(hsep(6,lpx+1+ndebug),stat=i_stat)
  call memocc(i_stat,hsep,'hsep',subname)

  !assignation of radii and coefficients of the local part
  alpz=psppar(0,0)
  alpl=psppar(0,0)
  alps(1:lpx+1)=psppar(1:lpx+1,0)
  gpot(1:3)=psppar(0,1:3)

  !assignation of the coefficents for the nondiagonal terms
  if (npspcode == 2) then !GTH case
     do l=1,lpx+1
        hsep(1,l)=psppar(l,1)
        hsep(2,l)=0.0_gp
        hsep(3,l)=psppar(l,2)
        hsep(4,l)=0.0_gp
        hsep(5,l)=0.0_gp
        hsep(6,l)=psppar(l,3)
     end do
  else if (npspcode == 3) then !HGH case
     allocate(ofdcoef(3,4+ndebug),stat=i_stat)
     call memocc(i_stat,ofdcoef,'ofdcoef',subname)

     ofdcoef(1,1)=-0.5_gp*sqrt(3._gp/5._gp) !h2
     ofdcoef(2,1)=0.5_gp*sqrt(5._gp/21._gp) !h4
     ofdcoef(3,1)=-0.5_gp*sqrt(100.0_gp/63._gp) !h5

     ofdcoef(1,2)=-0.5_gp*sqrt(5._gp/7._gp) !h2
     ofdcoef(2,2)=1._gp/6._gp*sqrt(35._gp/11._gp) !h4
     ofdcoef(3,2)=-7._gp/3._gp*sqrt(1._gp/11._gp) !h5

     ofdcoef(1,3)=-0.5_gp*sqrt(7._gp/9._gp) !h2
     ofdcoef(2,3)=0.5_gp*sqrt(63._gp/143._gp) !h4
     ofdcoef(3,3)=-9._gp*sqrt(1._gp/143._gp) !h5

     ofdcoef(1,4)=0.0_gp !h2
     ofdcoef(2,4)=0.0_gp !h4
     ofdcoef(3,4)=0.0_gp !h5

     !define the values of hsep starting from the pseudopotential file
     do l=1,lpx+1
        hsep(1,l)=psppar(l,1)
        hsep(2,l)=psppar(l,2)*ofdcoef(1,l)
        hsep(3,l)=psppar(l,2)
        hsep(4,l)=psppar(l,3)*ofdcoef(2,l)
        hsep(5,l)=psppar(l,3)*ofdcoef(3,l)
        hsep(6,l)=psppar(l,3)
     end do
     i_all=-product(shape(ofdcoef))*kind(ofdcoef)
     deallocate(ofdcoef,stat=i_stat)
     call memocc(i_stat,i_all,'ofdcoef',subname)
  else if (npspcode == 10) then !HGH-K case
     do l=1,lpx+1
        hsep(1,l)=psppar(l,1) !h11
        hsep(2,l)=psppar(l,4) !h12
        hsep(3,l)=psppar(l,2) !h22
        hsep(4,l)=psppar(l,5) !h13
        hsep(5,l)=psppar(l,6) !h23
        hsep(6,l)=psppar(l,3) !h33
     end do
  end if

  !Now the treatment of the occupation number
  call eleconf(izatom,ielpsp,symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg)

  occup(:,:)=0.0_gp
   do l=0,lmax
     iocc=0
     do i=1,6
        ott(i)=real(neleconf(i,l+1),gp)
        if (ott(i) > 0.0_gp) then
           iocc=iocc+1
            if (iocc > noccmax) stop 'iguess_generator: noccmax too small'
           occup(iocc,l+1)=ott(i)
        endif
     end do
     nl(l+1)=iocc
  end do

  !allocate arrays for the gatom routine
  allocate(vh(4*(ng+1)**2,4*(ng+1)**2+ndebug),stat=i_stat)
  call memocc(i_stat,vh,'vh',subname)
  allocate(psi(0:ng,noccmax,lmax+1+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)
  allocate(xp(0:ng+ndebug),stat=i_stat)
  call memocc(i_stat,xp,'xp',subname)
  allocate(rmt(n_int,0:ng,0:ng,lmax+1+ndebug),stat=i_stat)
  call memocc(i_stat,rmt,'rmt',subname)

  zion=real(ielpsp,gp)

  !can be switched on for debugging
  !if (iproc.eq.0) write(*,'(1x,a,a7,a9,i3,i3,a9,i3,f5.2)')&
  !     'Input Guess Generation for atom',trim(atomname),&
  !     'Z,Zion=',izatom,ielpsp,'ng,rprb=',ng+1,rprb

  rij=3._gp
  ! exponents of gaussians
  a0in=alpz
  a0=a0in/rij
  !       tt=sqrt(sqrt(2._gp))
  tt=2._gp**.3_gp
  do i=0,ng
     a=a0*tt**i
     xp(i)=.5_gp/a**2
  end do

  ! initial guess
  do l=0,lmax
     do iocc=1,noccmax
        do i=0,ng
           psi(i,iocc,l+1)=0.0_gp
        end do
     end do
  end do

  call crtvh(ng,lmax,xp,vh,rprb,fact,n_int,rmt)

  call gatom(rcov,rprb,lmax,lpx,noccmax,occup,&
       zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
       aeval,ng,psi,res,chrg)

  !post-treatment of the inguess data
  do i=1,ng+1
     expo(i)=sqrt(0.5_gp/xp(i-1))
  end do

  i=0
  do l=1,4
     do iocc=1,nl(l)
        i=i+1
        occupat(i)=occup(iocc,l)
        do j=1,ng+1
           psiat(j,i)=psi(j-1,iocc,l)
        end do
     end do
  end do

  i_all=-product(shape(vh))*kind(vh)
  deallocate(vh,stat=i_stat)
  call memocc(i_stat,i_all,'vh',subname)
  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)
  i_all=-product(shape(xp))*kind(xp)
  deallocate(xp,stat=i_stat)
  call memocc(i_stat,i_all,'xp',subname)
  i_all=-product(shape(rmt))*kind(rmt)
  deallocate(rmt,stat=i_stat)
  call memocc(i_stat,i_all,'rmt',subname)
  i_all=-product(shape(hsep))*kind(hsep)
  deallocate(hsep,stat=i_stat)
  call memocc(i_stat,i_all,'hsep',subname)
  i_all=-product(shape(alps))*kind(alps)
  deallocate(alps,stat=i_stat)
  call memocc(i_stat,i_all,'alps',subname)

END SUBROUTINE iguess_generator


subroutine gatom(rcov,rprb,lmax,lpx,noccmax,occup,&
                 zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,nintp,&
                 aeval,ng,psi,res,chrg)
  use module_base, only: gp
  implicit real(gp) (a-h,o-z)
  logical :: noproj
  integer, parameter :: n_int=100
  dimension psi(0:ng,noccmax,lmax+1),aeval(noccmax,lmax+1),&
       hh(0:ng,0:ng),ss(0:ng,0:ng),eval(0:ng),evec(0:ng,0:ng),&
       aux(2*ng+2),&
       gpot(3),hsep(6,lpx+1),rmt(n_int,0:ng,0:ng,lmax+1),&
       pp1(0:ng,lpx+1),pp2(0:ng,lpx+1),pp3(0:ng,lpx+1),alps(lpx+1),&
       potgrd(n_int),&
       rho(0:ng,0:ng,lmax+1),rhoold(0:ng,0:ng,lmax+1),xcgrd(n_int),&
       occup(noccmax,lmax+1),chrg(noccmax,lmax+1),&
       vh(0:ng,0:ng,4,0:ng,0:ng,4),&
       res(noccmax,lmax+1),xp(0:ng)
  if (nintp.ne.n_int) stop 'n_int><nintp'

  do l=0,lmax
     if (occup(1,l+1).gt.0._gp) lcx=l
  end do
  !write(6,*) 'lcx',lcx
 
  noproj=.true.
  do l=1,lpx+1
     noproj = noproj .and. (alps(l) .eq. 0._gp)
  end do


! projectors, just in case
  if (.not. noproj) then
     do l=0,lpx
        gml1=sqrt( gamma(real(l,gp)+1.5_gp) / (2._gp*alps(l+1)**(2*l+3)) )
        gml2=sqrt( gamma(real(l,gp)+3.5_gp) / (2._gp*alps(l+1)**(2*l+7)) )&
            /(real(l,gp)+2.5_gp)
        gml3=sqrt( gamma(real(l,gp)+5.5_gp) / (2._gp*alps(l+1)**(2*l+11)) )&
            /((real(l,gp)+3.5_gp)*(real(l,gp)+4.5_gp))
        tt=1._gp/(2._gp*alps(l+1)**2)
        do i=0,ng
           ttt=1._gp/(xp(i)+tt)
           pp1(i,l+1)=gml1*(sqrt(ttt)**(2*l+3))
           pp2(i,l+1)=gml2*ttt*(sqrt(ttt)**(2*l+3))
           pp3(i,l+1)=gml3*ttt**2*(sqrt(ttt)**(2*l+3))
        end do
     end do
  else
     pp1(:,:)=0._gp
     pp2(:,:)=0._gp
     pp3(:,:)=0._gp
  end if

  do l=0,lmax
     do j=0,ng
        do i=0,ng
           rho(i,j,l+1)=0._gp
        end do
     end do
  end do

  evsum=1.d30
  big_loop: do it=1,50
     evsumold=evsum
     evsum=0._gp
     
! coefficients of charge density
     do l=0,lmax
        do j=0,ng
           do i=0,ng
              rhoold(i,j,l+1)=rho(i,j,l+1)
              rho(i,j,l+1)=0._gp        
           end do
        end do
     end do

     do l=0,lmax
        do iocc=1,noccmax
           if (occup(iocc,l+1).gt.0._gp) then
              do j=0,ng
                 do i=0,ng
                    rho(i,j,l+1)=rho(i,j,l+1) + &
                         psi(i,iocc,l+1)*psi(j,iocc,l+1)*occup(iocc,l+1)
                 end do
              end do
           end if
        end do
     end do


     rmix=.5_gp
     if (it.eq.1) rmix=1._gp
     do l=0,lmax
        do j=0,ng
           do i=0,ng
              tt=rmix*rho(i,j,l+1) + (1._gp-rmix)*rhoold(i,j,l+1)
              rho(i,j,l+1)=tt
           end do
        end do
     end do

! XC potential on grid
!        do k=1,n_int
!           xcgrd(k)=0._gp
!        end do
!        do l=0,lmax
!           do j=0,ng
!              do i=0,ng
!                 do k=1,n_int
!                    xcgrd(k)=xcgrd(k)+rmt(k,i,j,l+1)*rho(i,j,l+1)
!                 end do
!              end do
!           end do
!        end do
     call DGEMV('N',n_int,(lcx+1)*(ng+1)**2,1._gp,&
                rmt,n_int,rho,1,0._gp,xcgrd,1)

     dr=fact*rprb/real(n_int,gp)
     do k=1,n_int
        r=(real(k,gp)-.5_gp)*dr
! divide by 4 pi
        tt=xcgrd(k)*0.07957747154594768_gp
! multiply with r^2 to speed up calculation of matrix elements
        xcgrd(k)=emuxc(tt)*r**2
     end do

     loop_l: do l=0,lmax
        gml=.5_gp*gamma(.5_gp+real(l,gp))

!  lower triangles only
        loop_i: do i=0,ng
           loop_j: do j=0,i
              d=xp(i)+xp(j)
              sxp=1._gp/d
              const=gml*sqrt(sxp)**(2*l+1)
! overlap
              ss(i,j)=const*sxp*(real(l,gp)+.5_gp)
! kinetic energy
              hh(i,j)=.5_gp*const*sxp**2* ( 3._gp*xp(i)*xp(j) +&
                   real(l,gp)*(6._gp*xp(i)*xp(j)-xp(i)**2-xp(j)**2) -&
                   real(l,gp)**2*(xp(i)-xp(j))**2  ) + .5_gp*real(l,gp)*(real(l,gp)+1._gp)*const
! potential energy from parabolic potential
              hh(i,j)=hh(i,j) +&
                   .5_gp*const*sxp**2*(real(l,gp)+.5_gp)*(real(l,gp)+1.5_gp)/rprb**4 
! hartree potential from ionic core charge
              tt=sqrt(1._gp+2._gp*alpz**2*d)
              if (l.eq.0) then
                 hh(i,j)=hh(i,j) -zion/(2._gp*d*tt)
              else if (l.eq.1) then
                 hh(i,j)=hh(i,j) -zion* &
                      (1._gp + 3._gp*alpz**2*d)/(2._gp*d**2*tt**3)
              else if (l.eq.2) then
                 hh(i,j)=hh(i,j) -zion* &
                      (2._gp + 10._gp*alpz**2*d + 15._gp*alpz**4*d**2)/(2._gp*d**3*tt**5)
              else if (l.eq.3) then
                 hh(i,j)=hh(i,j) -zion*3._gp* &
                      (2._gp +14._gp*alpz**2*d +35._gp*alpz**4*d**2 +35._gp*alpz**6*d**3)/&
                      (2._gp*d**4*tt**7)
              else 
                 stop 'l too big'
              end if
! potential from repulsive gauss potential
              tt=alpl**2/(.5_gp+d*alpl**2)
              hh(i,j)=hh(i,j)+ gpot(1)*.5_gp*gamma(1.5_gp+real(l,gp))*tt**(1.5_gp+real(l,gp))&
                   + (gpot(2)/alpl**2)*.5_gp*gamma(2.5_gp+real(l,gp))*tt**(2.5_gp+real(l,gp))&
                   + (gpot(3)/alpl**4)*.5_gp*gamma(3.5_gp+real(l,gp))*tt**(3.5_gp+real(l,gp))
! separable terms
              if (l.le.lpx) then
                 hh(i,j)=hh(i,j) + pp1(i,l+1)*hsep(1,l+1)*pp1(j,l+1)&
                      + pp1(i,l+1)*hsep(2,l+1)*pp2(j,l+1)&
                      + pp2(i,l+1)*hsep(2,l+1)*pp1(j,l+1)&
                      + pp2(i,l+1)*hsep(3,l+1)*pp2(j,l+1)&
                      + pp1(i,l+1)*hsep(4,l+1)*pp3(j,l+1)&
                      + pp3(i,l+1)*hsep(4,l+1)*pp1(j,l+1)&
                      + pp2(i,l+1)*hsep(5,l+1)*pp3(j,l+1)&
                      + pp3(i,l+1)*hsep(5,l+1)*pp2(j,l+1)&
                      + pp3(i,l+1)*hsep(6,l+1)*pp3(j,l+1)
              end if
! hartree potential from valence charge distribution
!              tt=0._gp
!              do lp=0,lcx
!                 do jp=0,ng
!                    do ip=0,ng
!                       tt=tt + vh(ip,jp,lp+1,i,j,l+1)*rho(ip,jp,lp+1)
!                    end do
!                 end do
!              end do
              tt=DDOT((lcx+1)*(ng+1)**2,vh(0,0,1,i,j,l+1),1,rho(0,0,1),1)
              hh(i,j)=hh(i,j) + tt
! potential from XC potential
              dr=fact*rprb/real(n_int,gp)
!              tt=0._gp
!              do k=1,n_int
!                 tt=tt+xcgrd(k)*rmt(k,i,j,l+1)
!              end do
              tt=DDOT(n_int,rmt(1,i,j,l+1),1,xcgrd(1),1)
              hh(i,j)=hh(i,j)+tt*dr
           end do loop_j
        end do loop_i

! ESSL
!        call DSYGV(1,hh,ng+1,ss,ng+1,eval,evec,ng+1,ng+1,aux,2*ng+2)
! LAPACK
        call DSYGV(1,'V','L',ng+1,hh,ng+1,ss,ng+1,eval,evec,(ng+1)**2,info)
        if (info.ne.0) write(6,*) 'LAPACK',info
        do iocc=0,noccmax-1
           do i=0,ng
              evec(i,iocc)=hh(i,iocc)
           end do
        end do
! end LAPACK
        do iocc=1,noccmax
           evsum=evsum+eval(iocc-1)
           aeval(iocc,l+1)=eval(iocc-1)
           do i=0,ng
              psi(i,iocc,l+1)=evec(i,iocc-1)
           end do
        end do
!        write(6,*) 'eval',l
!55      format(5(e14.7))
!        write(6,55) eval 
!        write(6,*) 'diff eval'
!        write(6,55) (eval(i)-eval(i-1),i=1,ng)
!        write(6,*) 'evec',l
!33      format(10(e9.2))
!        do i=0,ng
!           write(6,33) (evec(i,iocc),iocc=0,noccmax-1)
!        end do

     end do loop_l

     tt=abs(evsum-evsumold)
!     write(6,*) 'evdiff',it,tt
     if (tt.lt.1.e-12_gp) then
         exit big_loop
     end if
  end do big_loop
! End of the big loop

  call resid(lmax,lpx,noccmax,rprb,xp,aeval,psi,rho,ng,res,&
             zion,alpz,alpl,gpot,pp1,pp2,pp3,alps,hsep,fact,n_int,&
             potgrd,xcgrd,noproj)

! charge up to radius rcov
  if (lmax.gt.3) stop 'cannot calculate chrg'
  do l=0,lmax
     do iocc=1,noccmax
        chrg(iocc,l+1)=0._gp
     end do
  end do

  do iocc=1,noccmax
     do j=0,ng
        do i=0,ng
           d=xp(i)+xp(j)
           sd=sqrt(d)
           terf=derf(sd*rcov) 
           texp=exp(-d*rcov**2)

           tt=0.4431134627263791_gp*terf/sd**3 - 0.5_gp*rcov*texp/d
           chrg(iocc,1)=chrg(iocc,1) + psi(i,iocc,1)*psi(j,iocc,1)*tt
           if (lmax.eq.0) then
              cycle
           end if
           tt=0.6646701940895686_gp*terf/sd**5 + &
              (-0.75_gp*rcov*texp - 0.5_gp*d*rcov**3*texp)/d**2
           chrg(iocc,2)=chrg(iocc,2) + psi(i,iocc,2)*psi(j,iocc,2)*tt
           if (lmax.eq.1) then
               cycle
           end if
           tt=1.661675485223921_gp*terf/sd**7 + &
              (-1.875_gp*rcov*texp-1.25_gp*d*rcov**3*texp-.5_gp*d**2*rcov**5*texp) &
              /d**3
           chrg(iocc,3)=chrg(iocc,3) + psi(i,iocc,3)*psi(j,iocc,3)*tt
           if (lmax.eq.2) then
              cycle
           end if
           tt=5.815864198283725_gp*terf/sd**9 + &
              (-6.5625_gp*rcov*texp - 4.375_gp*d*rcov**3*texp - &
              1.75_gp*d**2*rcov**5*texp - .5_gp*d**3*rcov**7*texp)/d**4
           chrg(iocc,4)=chrg(iocc,4) + psi(i,iocc,4)*psi(j,iocc,4)*tt
        end do
     end do
  end do

! writing lines suppressed
!!$        write(66,*)  lmax+1
!!$        write(66,*) ' #LINETYPE{1324}' 
!!$        write(66,*) ' $' 
!!$  do l=0,lmax
!!$           write(66,*) ' 161'
!!$     r=0._gp
!!$     do
!!$        tt= wave(ng,l,xp,psi(0,1,l+1),r)
!!$              write(66,*) r,tt
!!$        r=r+.025_gp
!!$        if(r > 4.00001_gp) exit
!!$     end do
!!$  end do
! writing lines suppressed
!!$        write(67,*) min(lmax+1,3)
!!$        write(67,*) ' #LINETYPE{132}'
!!$        write(67,*) ' #TITLE{FOURIER}' 
!!$        write(67,*) ' $'
  dr=6.28_gp/rprb/200._gp
!!$        write(67,*) ' 200'
  rk=0._gp
  loop_rk1: do 
     tt=0._gp
     do i=0,ng
        texp=exp(-.25_gp*rk**2/xp(i))
!        texp=exp(-.5_gp*energy/xp(i))
        sd=sqrt(xp(i))
        tt=tt+psi(i,1,1)*0.4431134627263791_gp*texp/sd**3
     end do
!!$           write(67,*) rk,tt
     rk=rk+dr
     if(rk > 6.28_gp/rprb-.5_gp*dr) exit loop_rk1
  end do loop_rk1
  if (lmax.ge.1) then
!!$           write(67,*) ' 200'
     rk=0._gp
     loop_rk2: do 
        tt=0._gp
        do i=0,ng
           texp=exp(-.25_gp*rk**2/xp(i))
           sd=sqrt(xp(i))
           tt=tt+psi(i,1,2)*0.2215567313631895_gp*rk*texp/sd**5
        end do
!!$              write(67,*) rk,tt
        rk=rk+dr
        if (rk > 6.28_gp/rprb-.5_gp*dr) exit loop_rk2
     end do loop_rk2
  end if
  if (lmax.ge.2) then
!!$           write(67,*) ' 200'
     rk=0._gp
     do 
        tt=0._gp
        do i=0,ng
           texp=exp(-.25_gp*rk**2/xp(i))
           sd=sqrt(xp(i))
           tt=tt+psi(i,1,3)*0.1107783656815948_gp*rk**2*texp/sd**7
        end do
!!$              write(67,*) rk,tt
        rk=rk+dr
        if (rk > 6.28_gp/rprb-.5_gp*dr) exit
     end do
  end if

END SUBROUTINE gatom



subroutine resid(lmax,lpx,noccmax,rprb,xp,aeval,psi,rho,&
                 ng,res,zion,alpz,alpl,gpot,pp1,pp2,pp3,alps,hsep,fact,n_int,&
                 potgrd,xcgrd,noproj)
  use module_base, only: gp
  implicit real(gp) (a-h,o-z)
  logical :: noproj
  dimension psi(0:ng,noccmax,lmax+1),rho(0:ng,0:ng,lmax+1),&
       gpot(3),pp1(0:ng,lmax+1),pp2(0:ng,lmax+1),pp3(0:ng,lmax+1),&
       alps(lmax+1),hsep(6,lmax+1),res(noccmax,lmax+1),xp(0:ng),&
       xcgrd(n_int),aeval(noccmax,lmax+1),potgrd(n_int)
  
! potential on grid 
  dr=fact*rprb/real(n_int,gp)
  do k=1,n_int
     r=(real(k,gp)-.5_gp)*dr
     potgrd(k)= .5_gp*(r/rprb**2)**2 - &
          zion*derf(r/(sqrt(2._gp)*alpz))/r &
          + exp(-.5_gp*(r/alpl)**2)*&
          ( gpot(1) + gpot(2)*(r/alpl)**2 + gpot(3)*(r/alpl)**4 )&
          + xcgrd(k)/r**2
     do j=0,ng
        do i=0,ng
           spi=1.772453850905516_gp
           d=xp(i)+xp(j)
           sd=sqrt(d)
           tx=exp(-d*r**2)
           tt=spi*derf(sd*r)
           u_gp=tt/(4._gp*sd**3*r)
           potgrd(k)=potgrd(k)+u_gp*rho(i,j,1)
           ud1=-tx/(4._gp*d**2) + 3._gp*tt/(8._gp*sd**5*r)
           if (lmax.ge.1) potgrd(k)=potgrd(k)+ud1*rho(i,j,2)
           ud2=-tx*(7._gp + 2._gp*d*r**2)/(8._gp*d**3) +&
               15._gp*tt/(16._gp*sd**7*r)
           if (lmax.ge.2) potgrd(k)=potgrd(k)+ud2*rho(i,j,3)
           ud3=-tx*(57._gp+22._gp*d*r**2+4._gp*d**2*r**4)/(16._gp*d**4) + &
               105._gp*tt/(32._gp*sd**9*r)
           if (lmax.ge.3) potgrd(k)=potgrd(k)+ud3*rho(i,j,4)
        end do
     end do
  end do

  loop_ll: do ll=0,lmax
     if (ll.le.lpx .and. .not. noproj) then
        rnrm1=1._gp/sqrt(.5_gp*gamma(real(ll,gp)+1.5_gp)*alps(ll+1)**(2*ll+3))
        rnrm2=1._gp/sqrt(.5_gp*gamma(real(ll,gp)+3.5_gp)*alps(ll+1)**(2*ll+7))
        rnrm3=1._gp/sqrt(.5_gp*gamma(real(ll,gp)+5.5_gp)*alps(ll+1)**(2*ll+11))
     end if
     loop_iocc: do iocc=1,noccmax
! separable part
        if (ll.le.lpx) then
           scpr1=DDOT(ng+1,psi(0,iocc,ll+1),1,pp1(0,ll+1),1)
           scpr2=DDOT(ng+1,psi(0,iocc,ll+1),1,pp2(0,ll+1),1)
           scpr3=DDOT(ng+1,psi(0,iocc,ll+1),1,pp3(0,ll+1),1)
        end if
        res(iocc,ll+1)=0._gp
        loop_j: do j=1,n_int
! wavefunction on grid
           r=(real(j,gp)-.5_gp)*dr
           psigrd = wave(ng,ll,xp,psi(0,iocc,ll+1),r)
! kinetic energy        
           rkin=0._gp
           do i=0,ng
              rkin=rkin + psi(i,iocc,ll+1) *  (&
                   xp(i)*(3._gp+2._gp*real(ll,gp)-2._gp*xp(i)*r**2)*exp(-xp(i)*r**2) )
           end do
           rkin=rkin*r**ll
! separable part
           if (ll.le.lpx .and. .not. noproj) then
              sep =& 
                   (scpr1*hsep(1,ll+1) + scpr2*hsep(2,ll+1) + scpr3*hsep(4,ll+1))&
                   *rnrm1*r**ll*exp(-.5_gp*(r/alps(ll+1))**2)   +&
                   (scpr1*hsep(2,ll+1) + scpr2*hsep(3,ll+1) + scpr3*hsep(5,ll+1))&
                   *rnrm2*r**(ll+2)*exp(-.5_gp*(r/alps(ll+1))**2)   +&
                   (scpr1*hsep(4,ll+1) + scpr2*hsep(5,ll+1) + scpr3*hsep(6,ll+1))&
                   *rnrm3*r**(ll+4)*exp(-.5_gp*(r/alps(ll+1))**2)
           else
              sep=0._gp
           end if
! residue
           tt=rkin+sep+(potgrd(j)-aeval(iocc,ll+1))*psigrd
!384        format(6(e12.5))
!12        format(i2,i2,e9.2,3(e12.5),e10.3)
           res(iocc,ll+1)=res(iocc,ll+1) + tt**2*dr
        end do loop_j
     end do loop_iocc
  end do loop_ll
!  do l=0,lmax
!     do iocc=1,noccmax
!        write(6,*) 'res',l,iocc,res(iocc,l+1)
!     end do
!  end do

END SUBROUTINE resid


subroutine crtvh(ng,lmax,xp,vh,rprb,fact,n_int,rmt)
  use module_base, only: gp
  implicit real(gp) (a-h,o-z)
  dimension vh(0:ng,0:ng,0:3,0:ng,0:ng,0:3),xp(0:ng),&
            rmt(n_int,0:ng,0:ng,lmax+1)
  if (lmax.gt.3) stop 'crtvh'

  dr=fact*rprb/real(n_int,gp)
  do l=0,lmax
     do k=1,n_int
        r=(real(k,gp)-.5_gp)*dr
        do j=0,ng
           do i=0,ng
              rmt(k,i,j,l+1)=(r**2)**l*exp(-(xp(i)+xp(j))*r**2)
           end do
        end do
     end do
  end do

  loop_j: do j=0,ng
     loop_i: do i=0,ng
        c=xp(i)+xp(j)
        loop_jp: do jp=0,ng
           loop_ip: do ip=0,ng
              d=xp(ip)+xp(jp)
              scpd=sqrt(c+d)
              vh(ip,jp,0,i,j,0)=0.2215567313631895_gp/(c*d*scpd)
              vh(ip,jp,1,i,j,0)=&
                   .1107783656815948_gp*(2._gp*c+3._gp*d)/(c*d**2*scpd**3)
              vh(ip,jp,2,i,j,0)=.05538918284079739_gp*&
                   (8._gp*c**2+20._gp*c*d+15._gp*d**2)/(c*d**3*scpd**5)
              vh(ip,jp,3,i,j,0)=.0830837742611961_gp*&
              (16._gp*c**3+56._gp*c**2*d+70._gp*c*d**2+35._gp*d**3)/&
                   (c*d**4*scpd**7)

              vh(ip,jp,0,i,j,1)=&
                   .1107783656815948_gp*(3._gp*c+2._gp*d)/(c**2*d*scpd**3)
              vh(ip,jp,1,i,j,1)=&
                   .05538918284079739_gp*(6._gp*c**2+15._gp*c*d+6._gp*d**2)/&
                   (c**2*d**2*scpd**5)
              vh(ip,jp,2,i,j,1)=.02769459142039869_gp*&
                   (24._gp*c**3+84._gp*c**2*d+105._gp*c*d**2+30._gp*d**3)/&
                   (c**2*d**3*scpd**7)
              vh(ip,jp,3,i,j,1)=0.04154188713059803_gp*&
                   (48._gp*c**4+216._gp*c**3*d+378._gp*c**2*d**2+&
                   315._gp*c*d**3+70._gp*d**4)/(c**2*d**4*scpd**9)

              vh(ip,jp,0,i,j,2)=&
                   .05538918284079739_gp*(15._gp*c**2+20._gp*c*d+8._gp*d**2)/&
                   (c**3*d*scpd**5)
              vh(ip,jp,1,i,j,2)=.02769459142039869_gp*&
                   (30._gp*c**3+105._gp*c**2*d+84._gp*c*d**2+24._gp*d**3)/&
                   (c**3*d**2*scpd**7)
              vh(ip,jp,2,i,j,2)=&
                   .2077094356529901_gp*(8._gp*c**4+36._gp*c**3*d+63._gp*c**2*d**2+&
                   36._gp*c*d**3+8._gp*d**4)/(c**3*d**3*scpd**9)
              vh(ip,jp,3,i,j,2)=&
                   .1038547178264951_gp*(48._gp*c**5+264._gp*c**4*d+594._gp*c**3*d**2+&
                   693._gp*c**2*d**3+308._gp*c*d**4+56._gp*d**5)/&
                   (c**3*d**4*scpd**11)

              vh(ip,jp,0,i,j,3)=.0830837742611961_gp*&
                   (35._gp*c**3+70._gp*c**2*d+56._gp*c*d**2+16._gp*d**3)/&
                   (c**4*d*scpd**7)
              vh(ip,jp,1,i,j,3)=&
                   .04154188713059803_gp*(70._gp*c**4+315._gp*c**3*d+378._gp*c**2*d**2+&
                   216._gp*c*d**3+48._gp*d**4)/(c**4*d**2*scpd**9)
              vh(ip,jp,2,i,j,3)=&
                   .1038547178264951_gp*(56._gp*c**5+308._gp*c**4*d+693._gp*c**3*d**2+&
                   594._gp*c**2*d**3+264._gp*c*d**4+48._gp*d**5)/&
                   (c**4*d**3*scpd**11)
              vh(ip,jp,3,i,j,3)=&
                   1.090474537178198_gp*(16._gp*c**6+104._gp*c**5*d+286._gp*c**4*d**2+&
                   429._gp*c**3*d**3+286._gp*c**2*d**4+104._gp*c*d**5+16._gp*d**6)/&
                   (c**4*d**4*scpd**13)
           end do loop_ip
        end do loop_jp
     end do loop_i
  end do loop_j

END SUBROUTINE crtvh

 function wave(ng,ll,xp,psi,r)
  use module_base, only: gp
  implicit real(gp) (a-h,o-z)
  dimension psi(0:ng),xp(0:ng)

  wave=0._gp
  do i=0,ng
     wave=wave + psi(i)*exp(-xp(i)*r**2)
  end do
  if(ll>0)then
     wave=wave*r**ll
  endif
end function wave


function emuxc(rho)
  use module_base, only: gp
  implicit real(gp) (a-h,o-z)
  parameter (a0p=.4581652932831429_gp,&
       a1p=2.217058676663745_gp,&
       a2p=0.7405551735357053_gp,&
       a3p=0.01968227878617998_gp)
  parameter (b1p=1.0_gp,&
       b2p=4.504130959426697_gp,&
       b3p=1.110667363742916_gp,&
       b4p=0.02359291751427506_gp)
  parameter (rsfac=.6203504908994000_gp,ot=1._gp/3._gp)
  parameter (c1=4._gp*a0p*b1p/3.0_gp,  c2=5.0_gp*a0p*b2p/3.0_gp+a1p*b1p,&
       c3=2.0_gp*a0p*b3p+4.0_gp*a1p*b2p/3.0_gp+2.0_gp*a2p*b1p/3.0_gp,&
       c4=7.0_gp*a0p*b4p/3.0_gp+5.0_gp*a1p*b3p/3.0_gp+a2p*b2p+a3p*b1p/3.0_gp,&
       c5=2.0_gp*a1p*b4p+4.0_gp*a2p*b3p/3.0_gp+2.0_gp*a3p*b2p/3.0_gp,&
       c6=5.0_gp*a2p*b4p/3.0_gp+a3p*b3p,c7=4.0_gp*a3p*b4p/3.0_gp)
  if(rho.lt.1.e-24_gp) then
    emuxc=0._gp
  else
    if(rho.lt.0._gp) write(6,*) ' rho less than zero',rho
    rs=rsfac*rho**(-ot)
    top=-rs*(c1+rs*(c2+rs*(c3+rs*(c4+rs*(c5+rs*(c6+rs*c7))))))
    bot=rs*(b1p+rs*(b2p+rs*(b3p+rs*b4p)))
    emuxc=top/(bot*bot)
  end if
end function emuxc


! restricted version of the Gamma function
 function gamma(x)
  use module_base, only: gp
  implicit real(gp) (a-h,o-z)

  if (x.le.0._gp) stop 'wrong argument for gamma'
  if (mod(x,1._gp).eq.0._gp) then
     ii=int(x)
     do i=2,ii
        gamma=gamma*real(i-1,gp)
     end do
  else if (mod(x,.5_gp).eq.0._gp) then
     ii=int(x-.5_gp)
!     gamma=sqrt(3.14159265358979_gp)
     gamma=1.772453850905516027_gp
     do i=1,ii
        gamma=gamma*(real(i,gp)-.5_gp)
     end do
  else
     stop 'wrong argument for gamma'
  end if
end function gamma

!  call psitospi(iproc,nproc,norbe,norbep,norbsc,nat,&
!       wfd%nvctr_c,wfd%nvctr_f,at%iatype,at%ntypes,&
!       at%iasctype,at%natsc,at%natpol,nspin,spinsgne,psi)
subroutine psitospi0(iproc,nproc,norbe,norbep,norbsc,nat,&
     & nvctr_c,nvctr_f,iatype,ntypes, &
     iasctype,natsc,natpol,nspin,spinsgne,psi)
  use module_base
  implicit none
  integer, intent(in) :: norbe,norbep,iproc,nproc,nat
  integer, intent(in) :: nvctr_c,nvctr_f
  integer, intent(in) :: ntypes
  integer, intent(in) :: norbsc,natsc,nspin
  integer, dimension(ntypes), intent(in) :: iasctype
  integer, dimension(nat), intent(in) :: iatype,natpol
  integer, dimension(norbe*nspin), intent(in) :: spinsgne
  real(kind=8), dimension(nvctr_c+7*nvctr_f,norbep*nspin), intent(out) :: psi
  !local variables
  character(len=*), parameter :: subname='psitospi0'
  logical :: myorbital,polarised
  integer :: iatsc,i_all,i_stat,ispin,ipolres,ipolorb,nvctr
  integer :: iorb,jorb,iat,ity,i
  real(kind=8) :: facu,facd
  real(kind=8), dimension(:,:), allocatable :: psi_o
  integer, dimension(2) :: iorbsc,iorbv

  !initialise the orbital counters
  iorbsc(1)=0
  iorbv(1)=norbsc
  !used in case of spin-polarisation, ignored otherwise
  iorbsc(2)=norbe
  iorbv(2)=norbsc+norbe


  if (iproc ==0) then
     write(*,'(1x,a)',advance='no')'Transforming AIO to spinors...'
  end if
  
  nvctr=nvctr_c+7*nvctr_f
  allocate(psi_o(nvctr,norbep+ndebug),stat=i_stat)
  call memocc(i_stat,psi_o,'psi_o',subname)

  do iorb=1,norbep
     do i=1,nvctr
        psi_o(i,iorb)=psi(i,iorb)
     end do
  end do

 
  call razero(nvctr*nspin*norbep,psi)
  

  do iorb=1,norbe
     jorb=iorb-iproc*norbep
     if (myorbital(iorb,nspin*norbe,iproc,nproc)) then
        if(spinsgne(jorb)>0.0d0) then
           facu=1.0d0
           facd=0.0d0
        else
           facu=0.0d0
           facd=1.0d0
        end if
        do i=1,nvctr
           psi(i,iorb*4-3) = facu*psi_o(i,iorb)
           psi(i,iorb*4-2) = .0d0*psi_o(i,iorb)
           psi(i,iorb*4-1) = facd*psi_o(i,iorb)
           psi(i,iorb*4)   = .0d0*psi_o(i,iorb)
        end do
     end if
  end do
  i_all=-product(shape(psi_o))*kind(psi_o)
  deallocate(psi_o,stat=i_stat)
  call memocc(i_stat,i_all,'psi_o',subname)

  if (iproc ==0) then
     write(*,'(1x,a)')'done.'
  end if

END SUBROUTINE psitospi0
