subroutine inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,Glr,nvirt,nspin,&
     orbs,orbse,orbsv,norbsc_arr,locrad,G,psigau,eks)
  use module_base
  use module_types
  use module_interfaces, except_this_one => inputguess_gaussian_orbitals
  implicit none
  integer, intent(in) :: iproc,nproc,nspin
  integer, intent(inout) :: nvirt
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: Glr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), intent(out) :: eks
  integer, dimension(at%natsc+1,nspin), intent(out) :: norbsc_arr
  real(gp), dimension(at%nat), intent(out) :: locrad
  type(orbitals_data), intent(out) :: orbse,orbsv
  type(gaussian_basis), intent(out) :: G
  real(wp), dimension(:,:,:), pointer :: psigau
  !local variables
  character(len=*), parameter :: subname='inputguess_gaussian_orbitals'
  integer, parameter :: ngx=31
  integer :: norbe,norbme,norbyou,i_stat,i_all,norbsc,nvirte
  integer :: ispin,jproc,ist,jpst,nspinorfororbse,noncoll
  type(communications_arrays) :: commsv
  logical, dimension(:,:,:), allocatable :: scorb
  integer, dimension(:), allocatable :: ng,iorbtolr
  integer, dimension(:,:), allocatable :: nl
  !real(gp), dimension(:), allocatable :: occupe,spinsgneovrlp
  real(gp), dimension(:,:), allocatable :: xp,occupat
  real(gp), dimension(:,:,:), allocatable :: psiat

  allocate(xp(ngx,at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,xp,'xp',subname)
  allocate(psiat(ngx,5,at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,psiat,'psiat',subname)
  allocate(occupat(5,at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,occupat,'occupat',subname)
  allocate(ng(at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,ng,'ng',subname)
  allocate(nl(4,at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,nl,'nl',subname)
  allocate(scorb(4,2,at%natsc+ndebug),stat=i_stat)
  call memocc(i_stat,scorb,'scorb',subname)

  !Generate the input guess via the inguess_generator
  !here we should allocate the gaussian basis descriptors 
  !the prescriptions can be found in the creation of psp basis
  call readAtomicOrbitals(iproc,ngx,xp,psiat,occupat,ng,nl,at,norbe,norbsc,nspin,&
       scorb,norbsc_arr,locrad)

  !in the non-collinear case the number of orbitals double
  if (orbs%nspinor == 4) then
     noncoll=2
  else
     noncoll=1
  end if

  if (iproc ==0) then
     write(*,'(1x,a,i0,a)')'Generating ',nspin*noncoll*norbe,' Atomic Input Orbitals'
     if (norbsc /=0)   write(*,'(1x,a,i0,a)')'  of which ',nspin*noncoll*norbsc,&
          ' are semicore orbitals'
  end if


  if (nvirt /= 0) then
     !Check for max number of virtual orbitals
     !the unoccupied orbitals available as a LCAO
     !this is well defined only for closed-shell systems
     nvirte=noncoll*norbe-max(orbs%norbu,orbs%norbd)
     if(nvirt == nvirte .and. nvirt/=0 .and. iproc==0) then
        write(*,'(1x,a)')&
             "WARNING: A smaller number of virtual orbitals may be needed for better convergence."
        write(*,'(1x,a,i0)')'         Put nvirte= ',nvirte
     end if
     if (nvirte < nvirt) then
        nvirt=nvirte
        if(iproc==0) write(*,'(1x,a,i3)')&
             "WARNING: Number of virtual orbitals is too large. New value: ",nvirt
     end if
  end if
  !no Davidson calculation if nvirt=0
  if (nvirt==0) nvirte=0

  !create the orbitals descriptors, for virtual and inputguess orbitals
  allocate(orbsv%norb_par(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,orbsv%norb_par,'orbsv%norb_par',subname)
  !davidson treatment for spin-pol case should be reworked
  call orbitals_descriptors(iproc,nproc,nvirte,nvirte,0,1,orbsv)
  !allocate the arrays and fill them properly
  allocate(orbsv%occup(orbsv%norb+ndebug),stat=i_stat)
  call memocc(i_stat,orbsv%occup,'orbsv%occup',subname)
  allocate(orbsv%spinsgn(orbsv%norb+ndebug),stat=i_stat)
  call memocc(i_stat,orbsv%spinsgn,'orbsv%spinsgn',subname)
  orbsv%occup(1:orbsv%norb)=1.0_gp
  orbsv%spinsgn(1:orbsv%norb)=1.0_gp

  !allocate communications arrays for virtual orbitals
  !warning: here the aim is just to calculate npsidim, should be fixed
  call allocate_comms(nproc,commsv,subname)
  call orbitals_communicators(iproc,nproc,Glr,orbsv,commsv)  
  call deallocate_comms(commsv,subname)

  !deallocation if no davidson calculation
  if (nvirt == 0) then
     i_all=-product(shape(orbsv%norb_par))*kind(orbsv%norb_par)
     deallocate(orbsv%norb_par,stat=i_stat)
     call memocc(i_stat,i_all,'orbsv%norb_par',subname)
  end if


  !!!orbitals descriptor for inguess orbitals
  nspinorfororbse=orbs%nspinor

  allocate(orbse%norb_par(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,orbse%norb_par,'orbse%norb_par',subname)
  !the number of orbitals to be considered is doubled 
  !in the case of a spin-polarised calculation
  !also for non-collinear case
  !nspin*noncoll is always <= 2
  call orbitals_descriptors(iproc,nproc,nspin*noncoll*norbe,noncoll*norbe,(nspin-1)*norbe,&
       nspinorfororbse,orbse)
  !allocate the arrays and fill them properly
  allocate(orbse%occup(orbse%norb+ndebug),stat=i_stat)
  call memocc(i_stat,orbse%occup,'orbse%occup',subname)
  allocate(orbse%spinsgn(orbse%norb+ndebug),stat=i_stat)
  call memocc(i_stat,orbse%spinsgn,'orbse%spinsgn',subname)
  ist=1
  do ispin=1,nspin
     orbse%spinsgn(ist:ist+norbe-1)=real(1-2*(ispin-1),gp)
     ist=norbe+1
  end do

  !this is the distribution procedure for cubic code
  !should be referred to another routine
  if (iproc == 0 .and. nproc > 1) then
     jpst=0
     do jproc=0,nproc-1
        norbme=orbse%norb_par(jproc)
        norbyou=orbse%norb_par(min(jproc+1,nproc-1))
        if (norbme /= norbyou .or. jproc == nproc-1) then
           !this is a screen output that must be modified
           write(*,'(3(a,i0),a)')&
                ' Processes from ',jpst,' to ',jproc,' treat ',norbme,' inguess orbitals '
           jpst=jproc+1
        end if
     end do
     !write(*,'(3(a,i0),a)')&
     !     ' Processes from ',jpst,' to ',nproc-1,' treat ',norbyou,' inguess orbitals '
  end if

  !allocate the gaussian coefficients for the number of orbitals which is needed
  allocate(psigau(norbe,orbse%nspinor,orbse%isorb+orbse%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,psigau,'psigau',subname)
  allocate(iorbtolr(orbse%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,iorbtolr,'iorbtolr',subname)

  !fill just the interesting part of the orbital
  call AtomicOrbitals(iproc,nproc,at,rxyz,norbe,orbse,norbsc,occupat,&
       ngx,xp,psiat,ng,nl,nspin,eks,scorb,G,&
       psigau(1,1,min(orbse%isorb+1,orbse%norb)),&
       iorbtolr)

  i_all=-product(shape(scorb))*kind(scorb)
  deallocate(scorb,stat=i_stat)
  call memocc(i_stat,i_all,'scorb',subname)
  i_all=-product(shape(xp))*kind(xp)
  deallocate(xp,stat=i_stat)
  call memocc(i_stat,i_all,'xp',subname)
  i_all=-product(shape(psiat))*kind(psiat)
  deallocate(psiat,stat=i_stat)
  call memocc(i_stat,i_all,'psiat',subname)
  i_all=-product(shape(occupat))*kind(occupat)
  deallocate(occupat,stat=i_stat)
  call memocc(i_stat,i_all,'occupat',subname)
  i_all=-product(shape(ng))*kind(ng)
  deallocate(ng,stat=i_stat)
  call memocc(i_stat,i_all,'ng',subname)
  i_all=-product(shape(nl))*kind(nl)
  deallocate(nl,stat=i_stat)
  call memocc(i_stat,i_all,'nl',subname)
  i_all=-product(shape(iorbtolr))*kind(iorbtolr)
  deallocate(iorbtolr,stat=i_stat)
  call memocc(i_stat,i_all,'iorbtolr',subname)


end subroutine inputguess_gaussian_orbitals



subroutine readAtomicOrbitals(iproc,ngx,xp,psiat,occupat,ng,nl,at,norbe,norbsc,nspin,&
     & scorb,norbsc_arr,locrad)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ngx,iproc,nspin
  integer, intent(out) :: norbe,norbsc
  type(atoms_data), intent(inout) :: at
  logical, dimension(4,2,at%natsc), intent(out) :: scorb
  integer, dimension(at%ntypes), intent(out) :: ng
  integer, dimension(4,at%ntypes), intent(out) :: nl
  integer, dimension(at%natsc+1,nspin), intent(out) :: norbsc_arr
  real(gp), dimension(ngx,at%ntypes), intent(out) :: xp
  real(gp), dimension(5,at%ntypes), intent(out) :: occupat
  real(gp), dimension(ngx,5,at%ntypes), intent(out) :: psiat
  real(gp), dimension(at%nat), intent(out) :: locrad
  !local variables
  character(len=*), parameter :: subname='readAtomicOrbitals'
  integer, parameter :: nmax=6,lmax=3
  character(len=2) :: symbol
  character(len=20) :: pspatomname
  integer :: ity,i,j,l,ifile,ng_fake,ierror,iatsc,iat,ipow,lsc,inorbsc
  integer :: ichg,nsccode,ispol,mxpl,mxchg,i_all,i_stat,ityp,ishell,ictotpsi
  integer :: norbat,iorbsc_count,niasc,nlsc,iexpo,ig,ishltmp
  real(gp) :: rcov,rprb,ehomo
  integer, dimension(nmax,0:lmax) :: neleconf

  ! Read the data file.
  nl(1:4,1:at%ntypes) = 0
  ng(1:at%ntypes) = 0
  xp(1:ngx,1:at%ntypes) = 0.0_gp
  psiat(1:ngx,1:5,1:at%ntypes) = 0.0_gp
  occupat(1:5,1:at%ntypes)= 0.0_gp

  loop_assign: do ity=1,at%ntypes

     if (iproc == 0 .and. verbose > 1) then
        write(*,'(1x,a,a6,a)',advance='no')&
             'Input wavefunction data for atom ',trim(at%atomnames(ity)),&
             ' NOT found, automatic generation...'
     end if
     
     !the default value for the gaussians is chosen to be 21
     ng(ity)=21
     
     call iguess_generator(iproc,at%nzatom(ity),at%nelpsp(ity),at%psppar(0,0,ity),&
          at%npspcode(ity),&
          ng(ity)-1,nl(1,ity),5,occupat(1:5,ity),xp(1:ng(ity),ity),&
          psiat(1:ng(ity),1:5,ity))
     
     if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)')'done.'
  
  end do loop_assign

  close(unit=24)

  ! number of orbitals, total and semicore
  norbe=0
  norbsc=0
  iatsc=0
  scorb(:,:,:)=.false.
  do iat=1,at%nat
     ity=at%iatype(iat)
     !experimental case: multiply by two for non-collinear case
     !norbat=nl(1,ity)+3*nl(2,ity)+5*nl(3,ity)+7*nl(4,ity)
     norbat=(nl(1,ity)+3*nl(2,ity)+5*nl(3,ity)+7*nl(4,ity))
     norbe=norbe+norbat
     nsccode=at%iasctype(ity)
     !calculate the localisation radius for the input orbitals 
     call eleconf(at%nzatom(ity),at%nelpsp(ity),symbol,rcov,rprb,ehomo,&
          neleconf,nsccode,mxpl,mxchg)
     locrad(iat)=5._gp/sqrt(abs(2._gp*ehomo))
     call charge_and_spol(at%natpol(iat),ichg,ispol)
     !correct in the case of input charge positioning
     if (ichg /=0) then
        call correct_semicore(at%atomnames(ity),6,3,ichg,neleconf,nsccode)
     end if
     if (nsccode/=0) then !the atom has some semicore orbitals
        iatsc=iatsc+1
        niasc=nsccode
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
        norbsc_arr(iatsc,1)=iorbsc_count
        norbsc=norbsc+iorbsc_count
        !if (iproc == 0) write(*,*) iat,nsccode,iorbsc_count,norbsc,scorb(:,:,iatsc)
     end if
  end do

  !orbitals which are non semicore
  norbsc_arr(at%natsc+1,1)=norbe-norbsc

  !duplicate the values in the case of spin-polarization
  if (nspin == 2) norbsc_arr(:,2)=norbsc_arr(:,1)

END SUBROUTINE readAtomicOrbitals

subroutine createAtomicOrbitals(iproc,nproc,at,&
     rxyz,norbe,norbep,norbsc,occupe,occupat,ngx,xp,psiat,ng,nl,&
     wfd,n1,n2,n3,hx,hy,hz,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nspin,psi,eks,scorb)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: norbe,norbep,ngx,iproc,nproc,n1,n2,n3
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, intent(in) :: norbsc,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(wavefunctions_descriptors), intent(in) :: wfd
  logical, dimension(4,2,at%natsc), intent(in) :: scorb
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  integer, dimension(at%ntypes), intent(inout) :: ng
  integer, dimension(4,at%ntypes), intent(inout) :: nl
  real(gp), dimension(ngx,at%ntypes), intent(inout) :: xp
  real(gp), dimension(5,at%ntypes), intent(inout) :: occupat
  real(gp), dimension(ngx,5,at%ntypes), intent(inout) :: psiat
  real(gp), intent(out) :: eks
  real(gp), dimension(nspin*norbe), intent(out) :: occupe
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbep), intent(out) :: psi
  !local variables
  character(len=*), parameter :: subname= 'createAtomicOrbitals'
  integer, parameter :: nterm_max=3,noccmax=2,nlmax=4,nlevmax=6
  character(len=2) :: symbol
  logical :: myorbital,polarised
  integer :: iatsc,i_all,i_stat,ispin,ipolres,ipolorb,ichg,nsccode,mxpl,mxchg
  integer :: iorb,jorb,iat,ity,i,ictot,inl,l,m,nctot,nterm,iocc,ictotpsi
  real(wp) :: scprw
  real(dp) :: scpr
  real(gp) :: rx,ry,rz,ek,occshell,rcov,rprb,ehomo,shelloccup
  logical, dimension(nlmax,noccmax) :: semicore
  integer, dimension(2) :: iorbsc,iorbv
  integer, dimension(nlevmax,0:(nlmax-1)) :: neleconf
  integer, dimension(nterm_max) :: lx,ly,lz
  integer, dimension(nlmax) :: nlat
  real(gp), dimension(5) :: occupatat
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

  if (iproc == 0 .and. verbose > 1) then
     write(*,'(1x,a)',advance='no')'Calculating AIO wavefunctions...'
  end if

  do iat=1,at%nat

     rx=rxyz(1,iat)
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)

     ity=at%iatype(iat)

     !copy the nls and the occupation numbers
     nlat(:)=nl(:,ity)
     occupatat(:)=occupat(:,ity)

     !calculate the charge to be placed on an atom then correct
     !the semicore value if it is the case
     call charge_and_spol(at%natpol(iat),ichg,ipolres)
     nsccode=at%iasctype(ity)
     if (ichg /=0) then
        call eleconf(at%nzatom(ity),at%nelpsp(ity),symbol,rcov,rprb,ehomo,&
             neleconf,nsccode,mxpl,mxchg)
        call correct_semicore(at%atomnames(ity),nlevmax,nlmax-1,ichg,neleconf,nsccode)
        !we should then correct the occupatat and the nlat arrays for this atom
        iocc=0
        do l=0,nlmax-1
           inl=0
           do i=1,nlevmax
              if (neleconf(i,l) > 0) then
                 iocc=iocc+1
                 inl=inl+1
                 if (inl > noccmax) stop 'iguess_generator: noccmax too small'
                 occupatat(iocc)=real(neleconf(i,l),gp)
              endif
           end do
           nlat(l+1)=inl
        end do
     end if

     !the scorb array was already corrected in readAtomicOrbitals routine
     if (nsccode/=0) then !the atom has some semicore orbitals
        iatsc=iatsc+1
        semicore(:,:)=scorb(:,:,iatsc)
     else
        semicore(:,:)=.false.
     end if

     !calculate the atomic input orbitals
     ictot=0
     ictotpsi=0
     nctot=nl(1,ity)+nl(2,ity)+nl(3,ity)+nl(4,ity)
     if (iorbsc(1)+nctot > norbe .and. iorbv(1)+nctot > norbe) then
        write(*,*) 'transgpw occupe',nlat(:),norbe
        stop
     end if
     polarised=.false.
     do l=1,4
        do inl=1,nl(l,ity)
           ictotpsi=ictotpsi+1
           ictot=ictot+1
           !this can happen if we charge an atom such that an occupied shell disappears
           shelloccup=occupatat(ictot)
           if (inl > nlat(l)) then
              shelloccup=0.0_gp
              ictot=ictot-1
           end if
           !contribution to the kinetic energy given by the electrons in this shell
           call atomkin(l-1,ng(ity),xp(1,ity),psiat(1,ictotpsi,ity),psiatn,ek)
           eks=eks+ek*shelloccup
           if (nint(shelloccup) /=  2*(2*l-1) ) then
              !this is a polarisable orbital
              polarised=.true.
              !assuming that the control of the allowed polarisation is already done
              ipolorb=min(ipolres,int(shelloccup))
              ipolres=ipolres-ipolorb
              !this check can be inserted also elsewhere
              if (ipolres < 0) then
                 if(iproc==0) write(*,'(1x,4(a,i0))')&
                      'Too high polarisation for atom number= ',iat,&
                      ' Inserted=',modulo(at%natpol(iat),1000)-100,' Assigned=',ipolorb,&
                      ' the maximum is=',nint(shelloccup)
                 stop
              end if

           else
              !check for odd values of the occupation number
              if (mod(nint(shelloccup),2) /= 0) then
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
                 !write(*,*) 'iproc, SEMICORE orbital, iat,l',iproc,iat,l
                 !the occupation number is divided by two in the case of spin polarisation
                 if (nspin==2) then
                    occshell=0.5_gp*shelloccup
                 else
                    occshell=shelloccup
                 end if
              else
                 !normal case, the orbital is a valence orbital
                 iorb=iorbv(ispin)
                 occshell=shelloccup                 
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
                    call crtonewave(at%geocode,n1,n2,n3,ng(ity),nterm,lx,ly,lz,fac_arr,&
                         xp(1,ity),psiatn,rx,ry,rz,hx,hy,hz, & 
                         0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                         wfd%nseg_c,wfd%nvctr_c,wfd%keyg(1,1),wfd%keyv(1),wfd%nseg_f,wfd%nvctr_f,&
                         wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1),&
                         psi(1,jorb),psi(wfd%nvctr_c+1,jorb))
                    !renormalise wavefunction in case of too crude box
                    call wnrm(wfd%nvctr_c,wfd%nvctr_f,psi(1,jorb),psi(wfd%nvctr_c+1,jorb),scpr) 
                    !in the periodic case the function is not always normalised
                    !write(*,'(1x,a24,a7,2(a3,i1),a16,i4,i4,1x,1pe14.7)')&
                    !     'ATOMIC INPUT ORBITAL for atom',trim(at%atomnames(ity)),&
                    !     'l=',l,'m=',m,'iorb,jorb,norm',iorb,jorb,scpr
                    scprw=real(1.0_dp/sqrt(scpr),wp)
                    call wscal(wfd%nvctr_c,wfd%nvctr_f,scprw,psi(1,jorb),psi(wfd%nvctr_c+1,jorb))
                    !call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
                    !write(*,*) 'newnorm', scpr,occupe(iorb),occshell,ictot
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
     
     if (ictotpsi /= nctot) stop 'createAtomic orbitals: error (nctot)'

  end do
  if (iorbsc(1) /= norbsc) then
     write(*,*) iorbsc(1),norbsc
     stop 'createAtomic orbitals: error (iorbsc)'
  end if
  if (iorbv(1)/= norbe) stop 'createAtomic orbitals: error (iorbv)'
  if (iatsc /= at%natsc) stop 'createAtomic orbitals: error (iatsc)'

  if (nspin==2) then
     if (iorbsc(2)/= norbsc+norbe) stop 'createAtomic orbitals: error (iorbsc) nspin=2'
     if (iorbv(2) /= 2*norbe) stop 'createAtomic orbitals: error (iorbv) nspin=2'
  end if

  i_all=-product(shape(psiatn))*kind(psiatn)
  deallocate(psiatn,stat=i_stat)
  call memocc(i_stat,i_all,'psiatn',subname)

  if (iproc == 0 .and. verbose > 1) then
     write(*,'(1x,a)')'done.'
  end if

end subroutine createAtomicOrbitals

subroutine AtomicOrbitals(iproc,nproc,at,rxyz,norbe,orbse,norbsc,occupat,&
     ngx,xp,psiat,ng,nl,nspin,eks,scorb,G,gaucoeff,iorbtolr)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: norbe,ngx,iproc,nproc
  integer, intent(in) :: norbsc,nspin
  type(atoms_data), intent(in) :: at
  logical, dimension(4,2,at%natsc), intent(in) :: scorb
  real(gp), dimension(3,at%nat), intent(in), target :: rxyz
  type(orbitals_data), intent(inout) :: orbse
  integer, dimension(at%ntypes), intent(inout) :: ng
  integer, dimension(4,at%ntypes), intent(inout) :: nl
  real(gp), dimension(ngx,at%ntypes), intent(inout) :: xp
  real(gp), dimension(5,at%ntypes), intent(inout) :: occupat
  real(gp), dimension(ngx,5,at%ntypes), intent(inout) :: psiat
  type(gaussian_basis), intent(out) :: G
  real(gp), intent(out) :: eks
  integer, dimension(orbse%norbp), intent(out) :: iorbtolr !assign the localisation region
  real(wp), dimension(norbe,orbse%nspinor,orbse%norbp), intent(out) :: gaucoeff !norbe=G%ncoeff
  !local variables
  character(len=*), parameter :: subname= 'AtomicOrbitals'
  integer, parameter :: nterm_max=3,noccmax=2,nlmax=4,nlevmax=6
  character(len=2) :: symbol
  logical :: myorbital,polarised,orbpol_nc
  integer :: iatsc,i_all,i_stat,ispin,ipolres,ipolorb,ichg,nsccode,mxpl,mxchg,iexpo,ishltmp
  integer :: iorb,jorb,iat,ity,i,ictot,inl,l,m,nctot,nterm,iocc,ictotpsi,ishell,icoeff
  integer :: noncoll,ig,ispinor,icoll,norbpol_nc
  real(wp) :: scprw
  real(dp) :: scpr
  real(gp) :: rx,ry,rz,ek,occshell,rcov,rprb,ehomo,shelloccup,mx,my,mz,ma,mb,mc,md
  real(gp) :: mnorm,fac,occres
  logical, dimension(nlmax,noccmax) :: semicore
  integer, dimension(2) :: iorbsc,iorbv
  integer, dimension(nlevmax,0:(nlmax-1)) :: neleconf
  integer, dimension(nterm_max) :: lx,ly,lz
  integer, dimension(nlmax) :: nlat
  real(gp), dimension(5) :: occupatat
  real(gp), dimension(nterm_max) :: fac_arr
  real(gp), dimension(:), allocatable :: psiatn
  real(gp), dimension(:,:), allocatable :: atmoments

  if (iproc == 0 .and. verbose > 1) then
     write(*,'(1x,a)',advance='no')'Calculating AIO wavefunctions...'
  end if

  !gaussian basis structure informations
  !insert these things in the loops above
  !here we can create a gaussian basis structure for the input guess functions
  !the number of gaussian centers are thus nat
  G%nat=at%nat
  G%rxyz => rxyz
  !copy the parsed values in the gaussian structure
  !count also the total number of shells
  allocate(G%nshell(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,G%nshell,'G%nshell',subname)
  
  !if non-collinear it is like nspin=1 but with the double of orbitals
  if (orbse%nspinor == 4) then
     noncoll=2
  else
     noncoll=1
  end if

  G%nshltot=0
  do iat=1,at%nat
     ity=at%iatype(iat)
     G%nshell(iat)=(nl(1,ity)+nl(2,ity)+nl(3,ity)+nl(4,ity))
     G%nshltot=G%nshltot+G%nshell(iat)
  end do

  allocate(G%ndoc(G%nshltot+ndebug),stat=i_stat)
  call memocc(i_stat,G%ndoc,'G%ndoc',subname)
  allocate(G%nam(G%nshltot+ndebug),stat=i_stat)
  call memocc(i_stat,G%nam,'G%nam',subname)

  !assign shell IDs and count the number of exponents and coefficients
  G%nexpo=0
  G%ncoeff=0
  ishell=0
  do iat=1,at%nat
     ity=at%iatype(iat)
     ishltmp=0
     do l=1,4
        do i=1,nl(l,ity)
           ishell=ishell+1
           ishltmp=ishltmp+1
           G%ndoc(ishell)=ng(ity)
           G%nam(ishell)=l
           G%nexpo=G%nexpo+ng(ity)
           G%ncoeff=G%ncoeff+2*l-1
        end do
     end do
     if (ishltmp /= G%nshell(iat)) then
        write(*,*)'ERROR: ishelltmp <> nshell',ishell,G%nshell(iat)
        stop 
     end if
  end do

  if (norbe /= G%ncoeff) then
     write(*,*)'ERROR: norbe /= G%ncoeff',norbe,G%ncoeff
     stop 
  end if
  
  do jorb=1,orbse%norbp
     do ispinor=1,orbse%nspinor
        do icoeff=1,G%ncoeff
           gaucoeff(icoeff,ispinor,jorb)=0.0_wp
        end do
     end do
  end do

  !allocate and assign the exponents and the coefficients
  allocate(G%psiat(G%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,G%psiat,'G%psiat',subname)

  allocate(G%xp(G%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,G%xp,'G%xp',subname)

  allocate(psiatn(ngx+ndebug),stat=i_stat)
  call memocc(i_stat,psiatn,'psiatn',subname)

  !read the atomic moments if non-collinear
  !WARNING: units are not good for the moment
  if (orbse%nspinor == 4) then
     allocate(atmoments(3,at%nat+ndebug),stat=i_stat)
     call memocc(i_stat,atmoments,'atmoments',subname)

     open(unit=22,file='moments',form='formatted',action='read',status='old')
     !this part can be transferred on the atomic orbitals section
     do iat=1,at%nat
        read(unit=22,fmt=*,iostat=i_stat) mx,my,mz
        if (i_stat /= 0) then
           write(unit=*,fmt='(a,i0,a,i0,a)') &
                'The file "moments" has the line ',iat,&
                ' which have not 3 numbers for the atom ',iat,'.'
           stop 'The file "moments" is not correct!'
        end if
        atmoments(1,iat)=mx
        atmoments(2,iat)=my
        atmoments(3,iat)=mz
     end do
  end if


  eks=0.0_gp
  iorb=0
  iatsc=0

  !initialise the orbital counters
  iorbsc(1)=0
  iorbv(1)=norbsc
  !used in case of spin-polarisation, ignored otherwise
  iorbsc(2)=norbe
  iorbv(2)=norbsc+norbe

  ishell=0
  iexpo=0
  icoeff=1
  do iat=1,at%nat

     ity=at%iatype(iat)

     !copy the nls and the occupation numbers
     nlat(:)=nl(:,ity)
     occupatat(:)=occupat(:,ity)

     !calculate the charge to be placed on an atom then correct
     !the semicore value if it is the case
     call charge_and_spol(at%natpol(iat),ichg,ipolres)
     nsccode=at%iasctype(ity)
     if (ichg /=0) then
        call eleconf(at%nzatom(ity),at%nelpsp(ity),symbol,rcov,rprb,ehomo,&
             neleconf,nsccode,mxpl,mxchg)
        call correct_semicore(at%atomnames(ity),nlevmax,nlmax-1,ichg,neleconf,nsccode)
        !we should then correct the occupatat and the nlat arrays for this atom
        iocc=0
        do l=0,nlmax-1
           inl=0
           do i=1,nlevmax
              if (neleconf(i,l) > 0) then
                 iocc=iocc+1
                 inl=inl+1
                 if (inl > noccmax) stop 'iguess_generator: noccmax too small'
                 occupatat(iocc)=real(neleconf(i,l),gp)
              endif
           end do
           nlat(l+1)=inl
        end do
     end if

     !the scorb array was already corrected in readAtomicOrbitals routine
     if (nsccode/=0) then !the atom has some semicore orbitals
        iatsc=iatsc+1
        semicore(:,:)=scorb(:,:,iatsc)
     else
        semicore(:,:)=.false.
     end if

     !calculate the atomic input orbitals
     ictot=0
     ictotpsi=0
     nctot=(nl(1,ity)+nl(2,ity)+nl(3,ity)+nl(4,ity))
     if (iorbsc(1)+nctot > norbe .and. iorbv(1)+nctot > norbe) then
        print *,'transgpw occupe',nlat(:),norbe
        stop
     end if
     polarised=.false.
     
     do l=1,4
        do inl=1,nl(l,ity)
           ictotpsi=ictotpsi+1
           ictot=ictot+1
           ishell=ishell+1
           shelloccup=occupatat(ictot)
           !this can happen if we charge an atom such that an occupied shell disappears
           if (inl > nlat(l)) then
              shelloccup=0.0_gp
              ictot=ictot-1
           end if
           !contribution to the kinetic energy given by the electrons in this shell
           call atomkin(l-1,ng(ity),xp(1,ity),psiat(1,ictotpsi,ity),psiatn,ek)
           eks=eks+ek*shelloccup
           do ig=1,G%ndoc(ishell)
              iexpo=iexpo+1
              G%psiat(iexpo)=psiatn(ig)
              G%xp(iexpo)=xp(ig,ity)
           end do
           !decide the polarisation of the orbital by changing the population
           if (nint(shelloccup) /=  2*(2*l-1) ) then
              !this is a polarisable orbital
              polarised=.true.
              !assuming that the control of the allowed polarisation is already done
              ipolorb=min(ipolres,int(shelloccup))
              ipolres=ipolres-ipolorb
              !this check can be inserted also elsewhere
              if (ipolres < 0) then
                 if(iproc==0) write(*,'(1x,4(a,i0))')&
                      'Too high polarisation for atom number= ',iat,&
                      ' Inserted=',modulo(at%natpol(iat),1000)-100,' Assigned=',ipolorb,&
                      ' the maximum is=',nint(shelloccup)
                 stop
              end if

           else
              !check for odd values of the occupation number
              if (mod(nint(shelloccup),2) /= 0) then
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
                 if (nspin==2 .or. orbse%nspinor==4) then
                    occshell=0.5_gp*shelloccup
                 else
                    occshell=shelloccup
                 end if
              else
                 !normal case, the orbital is a valence orbital
                 iorb=iorbv(ispin)
                 occshell=shelloccup                 
                 if (nspin==2 .or. orbse%nspinor==4) then
                    if (polarised) then
                       occshell=0.5_gp*(occshell+real(1-2*(ispin-1),gp)*ipolorb)
                    else
                       occshell=0.5_gp*occshell
                    end if
                 end if
              end if

              !residue for the occupation number, to be used for
              !non-collinear case 
              occres=occshell
              !number of orbitals which will be polarised in this shell
              norbpol_nc=2*l-1
              do m=1,2*l-1
                 !each orbital has two electrons in the case of the 
                 !non-collinear case
                 do icoll=1,noncoll !non-trivial only for nspinor=4
                    iorb=iorb+1
                    jorb=iorb-orbse%isorb
                    
                    !the occupation number rule changes for non-collinear
                    if (orbse%nspinor == 4) then
                    !for each orbital of the shell, use the Hund rule
                    !for determining the occupation
                       if (ceiling(occres) >= real(2*l-1,gp)) then
                          !the orbital is not polarised
                          orbpol_nc=.false.
                          !assign the occupation to one (Hund's rule)
                          orbse%occup(iorb)=1.0_gp
                          !once finished this orbital (icoll=2), lower the occres
                          !and the number of polarisable orbitals
                          if (icoll==2) then
                             occres=occres-1.0_gp
                             norbpol_nc=norbpol_nc-1
                          end if
                       else
                          !in that case the orbital is polarised
                          !polarise uniformly the orbital following the direction
                          !indicated by the atmoments array
                          orbpol_nc=.true.
                          !the occupation number is assigned uniformly on the 
                          !remaining orbitals
                          !occupate only one of the orbitals
                          if (icoll ==1) then
                             orbse%occup(iorb)=2.0_gp*occres/real(norbpol_nc,gp)
                          else
                             orbse%occup(iorb)=0.0_gp
                          end if
                       end if
                       !print *,iorb,l,m,occshell,ceiling(occshell),norbpol_nc,orbpol_nc,orbse%occup(iorb)
                    else
                       orbse%occup(iorb)=occshell/real(2*l-1,gp)
                    end if
                    if (orbse%isorb < iorb .and. iorb <= orbse%isorb+orbse%norbp) then
                       if (orbse%nspinor == 1) then
                          do ispinor=1,orbse%nspinor
                             !here we put only the case nspinor==1
                             !we can put a phase for check with the complex wavefunction
                             gaucoeff(icoeff,ispinor,jorb)=1.0_wp
                          end do
                       else if (orbse%nspinor == 2) then
                          do ispinor=1,orbse%nspinor
                             !we can put a phase for check with the complex wavefunction
                             gaucoeff(icoeff,ispinor,jorb)=1.0_wp/sqrt(2.0_wp)
                          end do
                       else if (orbse%nspinor == 4) then
                          !assign the input orbitals according to the atomic moments
                          fac=0.5_gp

                          !if the shell is not polarised
                          !put one electron up and the other down
                          !otherwise, put a small unbalancing on the orbitals
                          !such that the momentum of the
                          if (orbpol_nc) then
                             !in the case of unoccupied orbital, 
                             !choose the orthogonal direction
                             mx=atmoments(1,iat)
                             my=atmoments(2,iat)
                             mz=atmoments(3,iat)

                             if (orbse%occup(iorb) == 0.0_gp) then
                                mx=-mx
                                my=-my
                                mz=-mz
                             end if
                          else
                             mx=0.0_gp
                             my=0.0_gp
                             mz=1.0_gp-2.0_gp*real(icoll-1,gp)
                          end if

                          mnorm=sqrt(mx**2+my**2+mz**2)
                          if (mnorm /= 0.0_gp) then
                             mx=mx/mnorm
                             my=my/mnorm
                             mz=mz/mnorm
                          end if

                          ma=0.0_gp
                          mb=0.0_gp
                          mc=0.0_gp
                          md=0.0_gp

                          if(mz > 0.0_gp) then 
                             ma=ma+mz
                          else
                             mc=mc+abs(mz)
                          end if
                          if(mx > 0.0_gp) then 
                             ma=ma+fac*mx
                             mb=mb+fac*mx
                             mc=mc+fac*mx
                             md=md+fac*mx
                          else
                             ma=ma-fac*abs(mx)
                             mb=mb-fac*abs(mx)
                             mc=mc+fac*abs(mx)
                             md=md+fac*abs(mx)
                          end if
                          if(my > 0.0_gp) then 
                             ma=ma+fac*my
                             mb=mb-fac*my
                             mc=mc+fac*my
                             md=md+fac*my
                          else
                             ma=ma-fac*abs(my)
                             mb=mb+fac*abs(my)
                             mc=mc+fac*abs(my)
                             md=md+fac*abs(my)
                          end if
                          if(mx==0.0_gp .and. my==0.0_gp .and. mz==0.0_gp) then
                             ma=1.0_gp/sqrt(2.0_gp)
                             mb=0.0_gp
                             mc=1.0_gp/sqrt(2.0_gp)
                             md=0.0_gp
                          end if

                          !assign the gaussian coefficients for each
                          !spinorial component
                          gaucoeff(icoeff,1,jorb)=real(ma,wp)
                          gaucoeff(icoeff,2,jorb)=real(mb,wp)
                          gaucoeff(icoeff,3,jorb)=real(mc,wp)
                          gaucoeff(icoeff,4,jorb)=real(md,wp)
                          !print *,'here',iat,jorb,gaucoeff(icoeff,:,jorb)
                       end if
                       
                       !associate to each orbital the reference localisation region
                       !here identified by the atom
                       iorbtolr(jorb)=iat 
                    endif
                 end do
                 icoeff=icoeff+1
              end do
              if (semicore(l,inl)) then
                 !increase semicore orbitals
                 iorbsc(ispin)=iorb
              else
                 !increase valence orbitals
                 iorbv(ispin)=iorb
              end if
              icoeff=icoeff-(2*l-1)
           end do
           icoeff=icoeff+(2*l-1)
        end do
     end do
     if (ictotpsi /= nctot) stop 'Atomic orbitals: error (nctot)'
  end do
  if (iexpo /= G%nexpo) then
     write(*,*)'ERROR: iexpo <> nexpo',iexpo,G%nexpo
     stop 
  end if

  !print *,'icoeff,ncoeff',icoeff,G%ncoeff

  if (iorbsc(1) /= norbsc) then
     print *,iorbsc(1),norbsc
     stop 'Atomic orbitals: error (iorbsc)'
  end if
  if (iorbv(1)/= noncoll*norbe) stop 'Atomic orbitals: error (iorbv)'
  if (iatsc /= at%natsc) stop 'Atomic orbitals: error (iatsc)'

  if (nspin==2) then
     if (iorbsc(2)/= norbsc+norbe) stop 'createAtomic orbitals: error (iorbsc) nspin=2'
     if (iorbv(2) /= 2*norbe) stop 'createAtomic orbitals: error (iorbv) nspin=2'
  end if

  !create a gaussian basis descriptor with all the information about the 


  i_all=-product(shape(psiatn))*kind(psiatn)
  deallocate(psiatn,stat=i_stat)
  call memocc(i_stat,i_all,'psiatn',subname)


  if (orbse%nspinor == 4) then
     i_all=-product(shape(atmoments))*kind(atmoments)
     deallocate(atmoments,stat=i_stat)
     call memocc(i_stat,i_all,'atmoments',subname)
  end if


  if (iproc ==0 .and. verbose > 1) then
     write(*,'(1x,a)')'done.'
  end if

end subroutine AtomicOrbitals


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
     write(*,*) 'l,m',l,m
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



! AMmodif  modified iguess_generator for absorbing atom


! The following two routines are copyrighted. I found them on the web at http://jin.ece.uiuc.edu/routines/routines.html
! when lloking for incomplete gamma function gamma(a,x)

        SUBROUTINE INCOG(A,X,GIN,GIM,GIP)
!
!       ===================================================
!       Purpose: Compute the incomplete gamma function
!                r(a,x),  and P(a,x)
!       Input :  a   --- Parameter ( a 
!                x   --- Argument 
!       Output:  GIN --- r(a,x)
!                GIM ---  gamma(a,x)
!                GIP --- P(a,x)
!       Routine called: GAMMA for computing 
!       ===================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XAM=-X+A*DLOG(X)
        IF (XAM.GT.700.0.OR.A.GT.170.0) THEN
           WRITE(*,*)'a and/or x too large'
           STOP
        ENDIF
        IF (X.EQ.0.0) THEN
           GIN=0.0
           GA= GAMMA(A)
           GIM=GA
           GIP=0.0
        ELSE IF (X.LE.1.0+A) THEN
           S=1.0D0/A
           R=S
           DO 10 K=1,60
              R=R*X/(A+K)
              S=S+R
              IF (DABS(R/S).LT.1.0D-15) GO TO 15
10         CONTINUE
15         GIN=DEXP(XAM)*S
           GA= GAMMA(A)
           GIP=GIN/GA
           GIM=GA-GIN
        ELSE IF (X.GT.1.0+A) THEN
           T0=0.0D0
           DO 20 K=60,1,-1
              T0=(K-A)/(1.0D0+K/(X+T0))
20         CONTINUE
           GIM=DEXP(XAM)/(X+T0)
           GA=GAMMA(A)
           GIN=GA-GIM
           GIP=1.0D0-GIM/GA
        ENDIF
        END

!!$        SUBROUTINE GAMMA(X,GA)
!!$!
!!$!       ==================================================
!!$!       Purpose: Compute gamma function)
!!$!       Input :  x  --- Argument of
!!$!                       ( x is not equal to 0,-1,-2,)
!!$!       Output:  GA ---
!!$!       ==================================================
!!$!
!!$        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!!$        DIMENSION G(26)
!!$        PI=3.141592653589793D0
!!$        IF (X.EQ.INT(X)) THEN
!!$           IF (X.GT.0.0D0) THEN
!!$              GA=1.0D0
!!$              M1=X-1
!!$              DO 10 K=2,M1
!!$10               GA=GA*K
!!$           ELSE
!!$              GA=1.0D+300
!!$           ENDIF
!!$        ELSE
!!$           IF (DABS(X).GT.1.0D0) THEN
!!$              Z=DABS(X)
!!$              M=INT(Z)
!!$              R=1.0D0
!!$              DO 15 K=1,M
!!$15               R=R*(Z-K)
!!$              Z=Z-M
!!$           ELSE
!!$              Z=X
!!$           ENDIF
!!$           DATA G/1.0D0,0.5772156649015329D0,   -0.6558780715202538D0, -0.420026350340952D-1,  0.1665386113822915D0,-.421977345555443D-1,     -.96219715278770D-2, .72189432466630D-2, -.11651675918591D-2, -.2152416741149D-3,     .1280502823882D-3, -.201348547807D-4, -.12504934821D-5, .11330272320D-5,               -.2056338417D-6, .61160950D-8,   .50020075D-8, -.11812746D-8,               .1043427D-9, .77823D-11,               -.36968D-11, .51D-12,               -.206D-13, -.54D-14, .14D-14, .1D-15/
!!$           GR=G(26)
!!$           DO 20 K=25,1,-1
!!$20            GR=GR*Z+G(K)
!!$           GA=1.0D0/(GR*Z)
!!$           IF (DABS(X).GT.1.0D0) THEN
!!$              GA=GA*R
!!$              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
!!$           ENDIF
!!$        ENDIF
!!$        RETURN
!!$        END


! end of copyrighted routines



subroutine find_pfproj( Nsol,Ngrid,rgrid,  psi1s, psigrid, real_start, psigrid_pseudo)
  use module_base
  implicit none
 
  integer, intent(in) :: Nsol,Ngrid, real_start
  real(gp), intent(in) :: psi1s(Ngrid), rgrid(Ngrid)
  real(gp), intent(inout) :: psigrid(Ngrid,Nsol),psigrid_pseudo(Ngrid,Nsol)
  ! ------------------------------------------------------------------------------------
  real(gp) dumgrid(Ngrid), coeffs(Nsol), dumgrid2(Ngrid), mass, mass_pseudo
  integer segno(Nsol), segno_pseudo(Nsol)
  integer i,j,k
  print *, " coefficienti su autostati " 
  do i=1, Nsol
     do k=1, Ngrid
        dumgrid(k)=psigrid(k,i)*psi1s(k)
     enddo
     call integrate(dumgrid, dumgrid2, rgrid, Ngrid)
     coeffs(i)=dumgrid2(Ngrid)
     print *, i, coeffs(i)

     mass=0.0D0
     mass_pseudo=0.0D0

     do k=1, Ngrid
        if( abs(psigrid(k,i))>mass) mass= abs(psigrid(k,i))
        if( abs(psigrid_pseudo(k,i))>mass_pseudo) mass_pseudo= abs(psigrid_pseudo(k,i))
     enddo
     

     do k=Ngrid, 1,-1
        if( abs(psigrid(k,i))>mass*0.01) then
           if( psigrid(k,i).gt.0.0D0) then
              segno(i)=1.0
           else
              segno(i)=-1
           endif
           exit
        endif
     enddo


     do k=Ngrid, 1,-1
        if( abs(psigrid_pseudo(k,i))>mass_pseudo*0.01) then
           if( psigrid_pseudo(k,i).gt.0.0D0) then
              segno_pseudo(i)=1
           else
              segno_pseudo(i)=-1
           endif
           exit
        endif
     enddo

  enddo



 

 call  DGEMM('N','N',Ngrid ,1,   Nsol,1.0d0 ,psigrid , Ngrid ,coeffs ,Nsol, 0.0D0 , dumgrid , Ngrid)


  do i=1,Nsol
     coeffs(i)=coeffs(i)*segno(i)*segno_pseudo(i)
  enddo

  call  DGEMM('N','N',Ngrid ,1,   Nsol-real_start+1  ,1.0d0 ,psigrid_pseudo , Ngrid ,&
       coeffs(real_start) ,Nsol-real_start+1, 0.0D0 , dumgrid2 , Ngrid)
  print * , " OK " 

  
  psigrid(:,1)=dumgrid
  psigrid_pseudo(:,1)=dumgrid2



  return
end subroutine find_pfproj
  


subroutine find_Scoeffs(Norder,Scoeffs,ng,gcoeffs,expo,noccmax,lmax,psi,&
     iocc_for_j,abs_final_L,psi1s,R, ene_for_j)
  use module_base
  implicit none
 
  integer, intent(in) :: Norder,ng, noccmax,  lmax
  real(gp), intent(in) :: R
  real(gp), dimension(ng+1), intent(in) :: expo
  integer, dimension(Norder), intent(in) :: iocc_for_j

  real(gp), dimension(0:ng,noccmax,lmax+1+ndebug), intent(in) :: psi
  real(gp), dimension(0:ng), intent(in) :: psi1s

  real(gp), dimension(Norder), intent(in) :: ene_for_j
  real(gp), dimension(Norder), intent(out) :: Scoeffs
  real(gp), dimension(0:ng), intent(out) :: gcoeffs
  integer, intent(in):: abs_final_L

  ! -----------------------------------------------------------
  real(gp)  ::   Soverlap(Norder,Norder)
  real(gp)  ::   Score(Norder)
  integer :: i,j,k,l,n,m, iw, INFO, LWORK, nord
  real(gp) :: a1,b1,a2,b2, A,B, spi, pi
  real(gp) :: W(Norder), WORK(4000)
  real(gp)::  sum, totalpow, gin, gim, gip, ggg, gamma

!!$  print *, " Norder, " , Norder
!!$  print *, " ng " , ng
!!$  print *, " noccmax " , noccmax
!!$  print *, " -------------- " 
!!$
!!$
!!$  print *, "  abs_final_L+1 ",  abs_final_L+1," lmax ", lmax, " iocc ",  iocc_for_j

  lwork=4000

  spi=1.772453850905516_gp
  pi=spi*spi


  do nord=1, Norder
     print *,  nord, " 1 ene_for_j  ",  ene_for_j(nord) 
  enddo




  !! prova con np

  totalpow=2.0_gp + abs_final_L +1

  do n=1,Norder
!!$     print *, n
     Score(n)=0.0
     iw = iocc_for_j(n)
     do i=0,ng
        a1 = psi1s(i)

        b1 = 1.0/2.0/expo(i+1)/expo(i+1)
        do j=0,ng
!!$           print *, "j ", j, " iw ", iw,"  abs_final_L+1 ",  abs_final_L+1, " iocc ",   iocc_for_j
           
           a2 = psi(j, iw, abs_final_L+1)
!!$           print *, " OK " 
           b2 = 1.0/2.0/expo(j+1)/expo(j+1)
!!$           print *, " ? " 
           A=a1*a2
           B=b1+b2

           call INCOG(  (1.0_gp+totalpow)/2.0_gp  ,B*R*R  ,GIN,GIM,GIP)
           ggg= gamma( (1.0_gp+totalpow)/2.0_gp  )
           Score(n)=Score(n)+A*0.5_gp* ( R**(1.0+totalpow) *(B*R*R)**(-0.5-totalpow/2 )*( ggg-GIM ) )
        enddo
     enddo
  enddo

  do nord=1, Norder
     print *,  nord, " 2 ene_for_j  ",  ene_for_j(nord) 
  enddo


!!$  print *, " A " 
  
  totalpow=2.0_gp + 2*abs_final_L

  do m=1,Norder
     do n=1,Norder
        Soverlap(m,n)=0.0
        do i=0,ng
           iw = iocc_for_j(m)
           a1 = psi(i, iw, abs_final_L+1)
           b1 = 1.0/2.0/expo(i+1)/expo(i+1)
           do j=0,ng
              iw = iocc_for_j(n)
              a2 = psi(j, iw, abs_final_L+1)
              b2 = 1.0/2.0/expo(j+1)/expo(j+1)
              
              A=a1*a2
              B=b1+b2
              

              call INCOG(  (1.0_gp+totalpow)/2.0_gp  ,B*R*R  ,GIN,GIM,GIP)
              ggg= gamma( (1.0_gp+totalpow)/2.0_gp   )
              Soverlap(m,n)=Soverlap(m,n)+A*0.5_gp* ( R**(1.0+totalpow) *(B*R*R)**(-0.5-totalpow/2 )*( ggg-GIM ) )
              
           enddo
        enddo

     enddo
  enddo
  


  print *, Score

  do nord=1, Norder
     print *,  nord, " 3 ene_for_j  ",  ene_for_j(nord) 
  enddo


!!$  print *, " B " 


  call DSYEV( 'V', 'L', Norder, Soverlap, Norder, W, WORK, -1, INFO )

  

  call DSYEV( 'V', 'L', Norder, Soverlap, Norder, W, WORK, LWORK, INFO )
  




  do nord=1, Norder
     print *,  nord, " 4 ene_for_j  ",  ene_for_j(nord) 
  enddo





  print *, " DSYEV on Soverlap , eigenvalues : ", W
  print *, " Soverlap " 
  print *, Soverlap

  Scoeffs(:)=0.0
  do n=1, Norder
     if( W(n).gt. 0.0000001*W(Norder)    ) then
        sum=0.0
        do i=1, Norder
           sum =  sum+ Soverlap(i,n)*Score(i)
        enddo
        sum=sum/W(n)
        do i=1, Norder
           Scoeffs(i)=Scoeffs(i)+Soverlap(i,n)*sum
        enddo
     endif
  enddo



!!$
!!$  do nord=1, Norder
!!$     Scoeffs(nord)=Score(nord)
!!$  enddo
!!$


  do nord=1, Norder
     print *,  nord, " 5 ene_for_j  ",  ene_for_j(nord) 
  enddo






  
  gcoeffs(:)=0.0
  do n=1, Norder
     do i=0,ng
        iw = iocc_for_j(n)
        gcoeffs(i) = gcoeffs(i) + Scoeffs(n)* psi(i, iw, abs_final_L+1)
     enddo
  enddo

  return 
end subroutine find_scoeffs








subroutine     find_Scoeffs_grid( ng,  expo, Ngrid, rgrid, psi1s , gcoeffs , l )
  use module_base
  implicit none
 
  integer, intent(in) :: ng,  l, Ngrid
  real(gp), intent(in) :: rgrid(Ngrid), psi1s(Ngrid)
  real(gp), dimension(ng+1), intent(in) :: expo
  real(gp), dimension(0:ng), intent(out) :: gcoeffs

  ! -----------------------------------------------------------
  real(gp)  ::   Soverlap(0:ng,0:ng)
  real(gp)  ::   Score(0:ng), dumgrid(Ngrid), dumgrid2(Ngrid)
  integer :: i,j,k,n,m, iw, INFO, LWORK, nord
  real(gp) :: a1,b1,a2,b2, A,B, spi, pi
  real(gp) :: W(0:ng), WORK(3*(ng+1)*(ng+1))
  real(gp)::  sum, totalpow, gin, gim, gip, ggg, gamma



  lwork= 3*(ng+1)*(ng+1)

  spi=1.772453850905516_gp
  pi=spi*spi

  totalpow=2.0_gp + 2*l  
  do i=0,ng
     b1 = 1.0/2.0/expo(i+1)/expo(i+1)
     do j=0,ng
        Soverlap(i,j)=0.0
        b2 = 1.0/2.0/expo(j+1)/expo(j+1)
        
        B=b1+b2
        
        ggg= gamma( (1.0_gp+totalpow)/2.0_gp   )
        Soverlap(i,j)=0.5_gp* ggg * B**(-(1.0D0+totalpow)/2 )
     enddo
     do k=1, Ngrid
        dumgrid(k)= psi1s(k)*exp( -b1*rgrid(k)*rgrid(k) ) *rgrid(k)**( l+1 ) 
     enddo
     call integrate(dumgrid, dumgrid2, rgrid, Ngrid)
     Score(i) = dumgrid2(Ngrid)
  enddo
  

  call DSYEV( 'V', 'L', ng+1, Soverlap(0,0) , ng+1, W(0), WORK, LWORK, INFO )
  



  gcoeffs(:)=0.0
  do n=0, ng
     if( W(n).gt. 1.0e-13*W(ng)    ) then
        sum=0.0
        do i=0, ng
           sum =  sum+ Soverlap(i,n)*Score(i)
        enddo
        sum=sum/W(n)
        do i=0,ng
           gcoeffs(i)=gcoeffs(i)+Soverlap(i,n)*sum
        enddo
     endif
  enddo
  

  return 
end subroutine find_Scoeffs_grid




subroutine dump_1gauwf_on_radgrid(prefix, ng , expo,psi   ,lpow   )
  use module_base
  implicit none
 
  character(*) , intent(in) ::  prefix
  integer, intent(in) :: ng,lpow
  real(gp), dimension(ng+1), intent(in) :: expo
  real(gp), dimension(0:ng), intent(in) :: psi

  ! local
  integer, parameter :: n_int=100
  character(200)  filename
  integer l,i, iocc,ig
  real(gp) r,sum



  write(filename,'(a)') prefix


  open(unit=22,file=filename)
  do i=1,2000!0.01,20,0.01
     r=0.01_gp*real(i-1,gp)+0.01_gp
     sum=0.0_gp
     do ig = 0,ng
        sum=sum+psi(ig)*exp( -r*r/2.0/expo(ig+1)/expo(ig+1) )
     enddo
     write(22,*) r, sum*(r**lpow)
  enddo
  close(unit=22)

return 
end subroutine dump_1gauwf_on_radgrid

real*8  function  value_at_r(r, ng , expo,psi     )
  use module_base, only: gp

  implicit none
 
  real(gp) , intent(in) ::  r
  integer, intent(in) :: ng
  real(gp), dimension(ng+1), intent(in) :: expo
  real(gp), dimension(0:ng), intent(in) :: psi

  ! local
  integer, parameter :: n_int=100

  integer l,i, iocc,ig
  real(8) sum



  sum=0.0
  do ig = 0,ng
     sum=sum+psi(ig)*exp( -r*r/2.0/expo(ig+1)/expo(ig+1) )
  enddo
  value_at_r=sum

  return 
end function value_at_r


subroutine dump_gauwf_on_radgrid(prefix  ,ng ,noccmax , lmax , expo,psi,aeval, occup     )
  use module_base, only: gp
  implicit none
 
  character(*) , intent(in) ::  prefix
  integer, intent(in) :: ng,noccmax, lmax
  real(gp), dimension(ng+1), intent(in) :: expo
  real(gp), dimension(0:ng,noccmax,lmax+1), intent(in) :: psi
  real(gp), dimension(noccmax,lmax+1  ), intent(in) ::  aeval,occup

  ! local
  integer, parameter :: n_int=100
  character(200)  filename
  integer l,i,j,iocc,ig
  real(8) r,sum

  do i=1,noccmax
     do l=0,lmax
  
        write(filename,'(a,a1,i1,a1,i1)') prefix,'_',i,'_',l
        
        open(unit=22,file=filename)
        do j=1,2000!0.01,20,0.01
           r=0.01_gp*real(j-1,gp)+0.01_gp
           sum=0.0
           do ig = 0,ng
              sum=sum+psi(ig,i,l+1)*exp( -r*r/2.0/expo(ig+1)/expo(ig+1) )
           enddo
           write(22,*) r, sum
        enddo
        close(unit=22)

        write(filename,'(a,a1,i1,a1,i1,a)') prefix,'_',i,'_',l,'_coeff'
        open(unit=22,file=filename)
        do ig = 0,ng
           write(22,*) psi(ig,i,l+1)
        enddo
        

        close(unit=22)

        


     enddo
  enddo
return
end subroutine dump_gauwf_on_radgrid



subroutine abs_generator_modified(iproc,izatom,ielpsp,psppar,npspcode,ng, noccmax, lmax ,expo,psi, aeval, occup, psp_modifier, &
     Nsol, Labs, Ngrid,Egrid,  rgrid , psigrid  )
  use module_base, only: gp, memocc
  implicit none
  integer, intent(in) :: iproc,izatom,ielpsp,ng,npspcode,noccmax, lmax, Nsol, labs, Ngrid
  real(gp), dimension(0:4,0:6), intent(in) :: psppar
  integer, intent(in) :: psp_modifier
  

  real(gp), dimension(ng+1), intent(out) :: expo

  integer, parameter :: n_int=100

  real(gp), dimension(0:ng,noccmax,lmax+1), intent(out) :: psi, Egrid(Nsol), rgrid(Ngrid), psigrid(Ngrid,Nsol  )
  real(gp), dimension(noccmax,lmax+1  ), intent(out) ::  aeval,occup
  
  !local variables
  character(len=*), parameter :: subname='iguess_generator'
  character(len=27) :: string 
  character(len=2) :: symbol
  real(gp), parameter :: fact=4.0_gp
  integer, dimension(6,4) :: neleconf
  real(gp), dimension(3) :: gpot
  real(gp), dimension(6) :: ott
  real(gp), dimension(noccmax,lmax+1) ::chrg,res
  real(gp), dimension(:), allocatable :: xp,alps
  real(gp), dimension(:,:), allocatable :: vh,hsep,ofdcoef

  real(gp), dimension(:,:,:,:), allocatable :: rmt
  logical :: exists
  integer :: lpx,ncount,nsccode,mxpl,mxchg
  integer :: l,i,j,iocc,il,lwrite,i_all,i_stat
  real(gp) :: alpz,alpl,rcov,rprb,zion,rij,a,a0,a0in,tt,ehomo
  real*8 value_at_r
  integer :: igrid

  !filename = 'psppar.'//trim(atomname)

  lpx=0
  if (psp_modifier.ne.0) then
     lpx = 0
  else
     lpx_determination: do i=1,4
        if (psppar(i,0) == 0.0_gp) then
           exit lpx_determination
        else
           lpx=i-1
        end if
     end do lpx_determination
  endif
  allocate(alps(lpx+1),stat=i_stat)
  call memocc(i_stat,alps,'alps',subname)
  allocate(hsep(6,lpx+1),stat=i_stat)
  call memocc(i_stat,hsep,'hsep',subname)

  !assignation of radii and coefficients of the local part


  if (psp_modifier.ne.0) then
     alpz=0.001_gp
     alpl=alpz
     alps(1:lpx+1)=0.0_gp
     gpot(1:3)=0.0_gp
     zion= izatom 
  else
     alpz=psppar(0,0)
     alpl=psppar(0,0)
     alps(1:lpx+1)=psppar(1:lpx+1,0)
     gpot(1:3)=psppar(0,1:3)
     zion=real(ielpsp,gp)
  endif

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
     allocate(ofdcoef(3,4),stat=i_stat)
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
  
  if(psp_modifier.ne.0) then
     call modified_eleconf(izatom,ielpsp,symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg)
  else
     call eleconf(izatom,ielpsp,symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg)
  endif

  occup(:,:)=0.000000000_gp
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

  end do

  !allocate arrays for the gatom routine
  allocate(vh(4*(ng+1)**2,4*(ng+1)**2),stat=i_stat)
  call memocc(i_stat,vh,'vh',subname)

  allocate(xp(0:ng),stat=i_stat)
  call memocc(i_stat,xp,'xp',subname)
  allocate(rmt(n_int,0:ng,0:ng,lmax+1),stat=i_stat)
  call memocc(i_stat,rmt,'rmt',subname)


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


  print *, "alpz ", alpz
  print *, "alpl ", alpl
  print *,  " gpot " , gpot
  print *, "lpx " , lpx
  print *, " alps", alps
  print *, " zion", zion
  print *, " hsep  ", hsep
  print *, " occup " , occup


!!$  call gatom(rcov,rprb,lmax,lpx,noccmax,occup,&
!!$       zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
!!$       aeval,ng,psi,res,chrg)



  call gatom_modified(rcov,rprb,lmax,lpx,noccmax,occup,&
                 zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
                 aeval,ng,psi,res,chrg,&
                 Nsol, Labs, Ngrid,Egrid,  rgrid , psigrid )

  print *, " --------------- AEVAl "
  print *, aeval

  !post-treatment of the inguess data
  do i=1,ng+1
     expo(i)=sqrt(0.5_gp/xp(i-1))
  end do


  do l=0,lmax
     do iocc=1,noccmax
        if( value_at_r(rprb, ng , expo,psi(0,iocc,l+1)).lt.0.0     ) then
           do i=0,ng
              psi(i,iocc,l+1)=-psi(i,iocc,l+1)
           enddo
        endif
     enddo
  enddo
  

  i_all=-product(shape(vh))*kind(vh)
  deallocate(vh,stat=i_stat)
  call memocc(i_stat,i_all,'vh',subname)
  i_all=-product(shape(psi))*kind(psi)

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

END SUBROUTINE abs_generator_modified














subroutine abs_generator(iproc,izatom,ielpsp,psppar,npspcode,ng, noccmax, lmax ,expo,psi, aeval, occup, psp_modifier )
  use module_base, only: gp, memocc
  implicit none
  integer, intent(in) :: iproc,izatom,ielpsp,ng,npspcode,noccmax, lmax
  real(gp), dimension(0:4,0:6), intent(in) :: psppar
  integer, intent(in) :: psp_modifier
  

  real(gp), dimension(ng+1), intent(out) :: expo

  integer, parameter :: n_int=100

  real(gp), dimension(0:ng,noccmax,lmax+1), intent(out) :: psi
  real(gp), dimension(noccmax,lmax+1  ), intent(out) ::  aeval,occup
  
  !local variables
  character(len=*), parameter :: subname='iguess_generator'
  character(len=27) :: string 
  character(len=2) :: symbol
  real(gp), parameter :: fact=4.0_gp
  integer, dimension(6,4) :: neleconf
  real(gp), dimension(3) :: gpot
  real(gp), dimension(6) :: ott
  real(gp), dimension(noccmax,lmax+1) ::chrg,res
  real(gp), dimension(:), allocatable :: xp,alps
  real(gp), dimension(:,:), allocatable :: vh,hsep,ofdcoef

  real(gp), dimension(:,:,:,:), allocatable :: rmt
  logical :: exists
  integer :: lpx,ncount,nsccode,mxpl,mxchg
  integer :: l,i,j,iocc,il,lwrite,i_all,i_stat
  real(gp) :: alpz,alpl,rcov,rprb,zion,rij,a,a0,a0in,tt,ehomo
  real*8 value_at_r

  !filename = 'psppar.'//trim(atomname)

  lpx=0
  if (psp_modifier.ne.0) then
     lpx = 0
  else
     lpx_determination: do i=1,4
        if (psppar(i,0) == 0.0_gp) then
           exit lpx_determination
        else
           lpx=i-1
        end if
     end do lpx_determination
  endif
  allocate(alps(lpx+1),stat=i_stat)
  call memocc(i_stat,alps,'alps',subname)
  allocate(hsep(6,lpx+1),stat=i_stat)
  call memocc(i_stat,hsep,'hsep',subname)

  !assignation of radii and coefficients of the local part


  if (psp_modifier.ne.0) then
     alpz=0.001_gp
     alpl=alpz
     alps(1:lpx+1)=0.0_gp
     gpot(1:3)=0.0_gp
     zion= izatom 
  else
     alpz=psppar(0,0)
     alpl=psppar(0,0)
     alps(1:lpx+1)=psppar(1:lpx+1,0)
     gpot(1:3)=psppar(0,1:3)
     zion=real(ielpsp,gp)
  endif

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
     allocate(ofdcoef(3,4),stat=i_stat)
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
  
  if(psp_modifier.ne.0) then
     call modified_eleconf(izatom,ielpsp,symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg)
  else
     call eleconf(izatom,ielpsp,symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg)
  endif

  occup(:,:)=0.000000000_gp
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

  end do

  !allocate arrays for the gatom routine
  allocate(vh(4*(ng+1)**2,4*(ng+1)**2),stat=i_stat)
  call memocc(i_stat,vh,'vh',subname)

  allocate(xp(0:ng),stat=i_stat)
  call memocc(i_stat,xp,'xp',subname)
  allocate(rmt(n_int,0:ng,0:ng,lmax+1),stat=i_stat)
  call memocc(i_stat,rmt,'rmt',subname)


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


  print *, "alpz ", alpz
  print *, "alpl ", alpl
  print *,  " gpot " , gpot
  print *, "lpx " , lpx
  print *, " alps", alps
  print *, " zion", zion
  print *, " hsep  ", hsep
  print *, " occup " , occup


  call gatom(rcov,rprb,lmax,lpx,noccmax,occup,&
       zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
       aeval,ng,psi,res,chrg)



  print *, " --------------- AEVAl "
  print *, aeval

  !post-treatment of the inguess data
  do i=1,ng+1
     expo(i)=sqrt(0.5_gp/xp(i-1))
  end do


  do l=0,lmax
     do iocc=1,noccmax
        if( value_at_r(rprb, ng , expo,psi(0,iocc,l+1)).lt.0.0     ) then
           do i=0,ng
              psi(i,iocc,l+1)=-psi(i,iocc,l+1)
           enddo
        endif
     enddo
  enddo
  

  i_all=-product(shape(vh))*kind(vh)
  deallocate(vh,stat=i_stat)
  call memocc(i_stat,i_all,'vh',subname)
  i_all=-product(shape(psi))*kind(psi)

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

END SUBROUTINE abs_generator


! AMmodif end




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


! ------------------------------------------------
  


! -----------------------------------------------




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



subroutine integrate(f,fint,x,Nx)
  use module_base, only: gp
  implicit none
  real(gp), intent(in) :: f(0:Nx-1), x(0:Nx-1)
  real(gp), intent(out):: fint(0:Nx-1)
  integer, intent(in):: Nx
  ! -------------------------
  real(gp) sfin, sgro
  integer i
  real(gp) d,h
  real(gp) p,q,r
  real(gp) fp,fq,fr
  
  
  
  if(Nx.lt.3) then
     print *, " error, Nx<3 in routine integrate "
     stop
  endif
  
  
  sfin=0.0_gp
  sgro=0.0_gp
  fint(0)=0.0_gp
  
  do i=2, Nx-1,2
     sgro=sgro+ (f(i-2)+ f(i))*(x(i)-x(i-2))
     sfin=sfin+(f(i-2) + f(i-1))*(x(i-1)-x(i-2))
     sfin=sfin+(f(i-1)+f(i))*(x(i)-x(i-1))
     fint(i)=(4*sfin-sgro)/6.0_gp
  enddo
  
  
  
  do i=1, Nx-1,2
     
     if( i .lt. Nx-1) then
        d=x(i)-x(i-1)
        h=x(i+1)-x(i)
        
        p=d*(2*d+3*h)/6/(d+h)
        q=d*(d+3*h)/6./h
        r=-d*d*d/6./h/(d+h)
        
        fp=f(i-1)
        fq=f(i)
        fr=f(i+1)
        
        fint(i)=fint(i-1)+p*fp+q*fq+r*fr
     else
        h=x(i-1)-x(i-2)
        d=x(i)-x(i-1)
        p=-d*d*d/6./h/(d+h)
        q=d*(d+3*h)/6./h
        r=d*(2*d+3*h)/6/(d+h)
        
        fp=f(i-2)
        fq=f(i-1)
        fr=f(i)
        fint(i)=fint(i-1)+p*fp+q*fq+r*fr
     endif
     
  enddo
  return
end subroutine integrate


real*8 function pow(x,n)
  real*8 x
  integer n
  pow=x**n
  return
end function pow

real*8 function phase( E, N, rgrid, V, nonloc, y,l,normalize,Z, onlyout)
  use module_base, only: gp
  implicit none
  integer N, normalize, onlyout,l
  real*8 E,rgrid(N),V(N), nonloc(N), Z, y(N)
  ! ----------------------------------------------
  integer ii, i,j
  real*8  ypi
  
  integer yNcross, yflag
  real*8  dh,dl,dh2,dl2,dp,dl3,dhdl2,dh2dl,dh3,add1,add2,deno,num
  real*8  Gh,Gl,G0,Nlc,Nla,Nlb
  real*8  ca,cb,ep,em,y1,y2
  
  real*8 Ga,Gb,Gc
  real*8 ya,yb,yc
  real*8 r, PI
  
  real*8 fact,norm, func(N), funcdum(N), pow
  integer count
  
 



  PI=4*atan(1.0d0)
  
  phase=0
  yNcross=0
  yflag=0
  
  ! -- Calcul du point intermediaire où le V est minimum
  if( onlyout.eq.1 ) then
     ii=N-1  
  else
     ii=0
     do i=1 , N-1
        if(V(i+1).ge. V(i) )  then
           ii=i
           if(i.eq.1)  then
              ii=10
              if(ii.ge.N-3) then
                 print *, "  attenzione !!! if(I>=n-3) in schro.cc ";
                 stop
              endif
           endif
           exit
        endif
     enddo
     
     if(ii.eq.0) then 
        ii=N-10;
        print *, " attenzione !!!I=N-1 in schro.cc  "
        print *, " l est " ,  l
        stop
     endif
     
     
     !  if(I<100) I=100;
     
  endif
  


  ! print * , " propagation de 1 a " , ii
  
  ! //------------ Propagation de  1 a I  ----------------------
  
  do i=1,2
     if(rgrid(i).eq.0.0 .and. l.ne.-1) then
        y(i)=0
     else if (l.eq.-4) then
        y(i)=0.0_gp
        print *," metto tutto a zero "
     else 
        y(i) = exp((l+1)*log(rgrid(i)))
     endif

  enddo
  
  
  do i=2,ii
     dh = rgrid(i+1)-rgrid(i)
     dl = rgrid(i)-rgrid(i-1)
     G0 = 2*(V(i)  -E)
     Nlb=nonloc(i)
     if(dabs(G0*dh*dl).lt.1.0)  then
        dh2= dh*dh
        dl2= dl*dl
        dp = dh*dl
        
        Gh = 2*(V(i+1)-E)
        Gl = 2*(V(i-1)-E)
        
        Nlc=nonloc(i+1)
        Nla=nonloc(i-1)
        
        R  = dh2 + dl2 + dp 
        
        dh3=dh2*dh
        dhdl2=dl*dp
        dh2dl=dh*dp
        dl3=dl*dl2
        
        ! // ********************************************************************************
        ! // ** 17 luglio 1998. Aggiunta della parte per il termine non locale.
        ! // ** Vedere file numerov.mathematica
        ! // ** ( h^3*(-Nla  + Nlb ) + l^3*(Nlb  - Nlc) + 
        ! // **      h^2*l*(Nla+4*Nlbd+Nlc)+h*l^2*(Nla+4*Nlb+Nlc))/l*(12-(h^2+h*l-l^2)*Gc))
        
        deno=(dl*(1._gp-Gh*(R-2*dl2)/12._gp))
        
        add1=( y(i)*(dh+(dl*Gh*(R-2*dl2)+G0*(dh+dl)*(R+2*dh*dl))/12) &
             -y(i-1)*dh*(1.-Gl*(R-2*dh2)/12.) &
             )/deno;
        
        
        add2=( dh3*(-Nla  + Nlb ) + dl3*(Nlb  - Nlc) + &
             dh2dl*(Nla  + 4*Nlb + Nlc ) + dhdl2*(Nla + 4*Nlb + Nlc) )/deno/6
        
        y(i+1) = y(i)+add1+add2


     else
        print *, " occorre prendere asso piu piccolo "
        stop
     endif
     
     
     if(dabs(y(i+1)).gt.1.0D40 .and. l.ne.-4) then 
        fact= 1./dabs(y(i+1)) 
        count =0
        do j=i+1,1,-1
           if(y(j).eq.0) then
              count=count+1
           else 
              count =0
           endif
           if(count .gt. 2) exit
           y(j)=y(j)*fact 
        enddo
     endif
     if( (i+1) .le.ii ) then !  // LA PHASE EST CALCULE EN I
        if((y(i)*y(i+1)).lt.0) yNcross=yNcross+1
        if(yflag == 0 .and. (y(i)*y(i+1)) == 0)then
           yflag=1-yflag
           yNcross=yNcross+1
        endif
     endif
  enddo
  
  if( y(ii).eq.0.0) then 
     print *, " y[I] == 0.0 dans Schroedinger xs"
     stop
  endif
  
  !  // ypI=(y[I+1]*dl2-y[I-1]*dh2+y[I]*(dh2-dl2))/((dh+dl)*dl*dh); // vecchia maniera di calcolare yp
  !  // sostituita con la maniera seguente, presa dal metodo numerov (vedi numerov.nb)
  
  
  
  
  i=ii
  dh = rgrid(i+1)-rgrid(i)
  dl = rgrid(i  )-rgrid(i-1)
  
  Nlc=nonloc(i+1)
  Nla=nonloc(i-1)
  Nlb=nonloc(i)
  
  Gb = 2*(V(i)  -E);
  if(dabs(G0*dh*dl).gt.1.0) then
     print *, " problema fabs(G0*dh*dl)>1.0 in calcolo di yp in I "
     stop
  end if

  Gc = 2*(V(i+1)-E)
  Ga = 2*(V(i-1)-E)
  
  ya=y(ii-1)
  yb=y(ii)
  yc=y(ii+1)
  
  ypI=(dh*(Ga*pow(dl,2)*(-12*dh + Gc*pow(dh,3) - 6*dl + &
       Gc*pow(dh,2)*dl) - &
       6 *(-12*dh + Gc*pow(dh,3) - 12*dl + 2*Gc*pow(dh,2)*dl))*ya + &
       (dh + dl)*(Gb*pow(dl,2)* &
       (-24*dh + 2*Gc*pow(dh,3) - 6*dl + 3*Gc*pow(dh,2)*dl) + &
       6 *(-12*dh + Gc*pow(dh,3) + Gc*pow(dh,2)*dl - &
       Gc*dh*pow(dl,2) + Gc*pow(dl,3)))*yb)/&
       (6.*dh*dl*(dh + dl)*(-12 + Gc*(pow(dh,2) + dh*dl - pow(dl,2)))) 
  
  ! // *****************************************************************************
  ! // ** Termine aggiuntivo
  ! // **
  ! // 6*l^3*Nlc + 
  ! // Nla*(-12*h^2*l - 6*h*l^2 + h^3*l*(h + l)*Gc) + 
  ! // Nlb*(-24*h^2*l - 30*h*l^2 - 6*l^3 +      2*h^3*l*(h + l)*Gc + 3*h^2*l^2*(h + l)*Gc)
  ! // / al denominatore
  ! //
  ! //  6*h*(h + l)*(-12 + (h^2 + h*l - l^2)*Gc)
  
  deno= 6*dh*(dh + dl)*(-12 + (dh*dh + dh*dl - dl*dl)*Gc)
	    
  num=  6*dl*dl*dl*Nlc&
       +Nla*(-12*dh*dh*dl-6*dh*dl*dl+dh*dh*dh*dl*(dh + dl)*Gc)&
       +Nlb*(-24*dh*dh*dl - 30*dh*dl*dl - &
       6*dl*dl*dl+2*dh*dh*dh*dl*(dh + dl)*Gc +3*dh*dh*dl*dl*(dh+dl)*Gc)

  ypI=ypI+num/deno
	    
  
  if( onlyout.eq.1) then
     ypI = ypI+ dh*(  y(ii)*Gb  +y(ii+1)*Gc )/2 
     phase=ypI
     return
  endif


  ! print * , "  y(ii), ypI " , y(ii), ypI
	  
	  
  if(dabs(y(ii)) .gt. dabs( ypI ) ) then
     phase=-atan(ypI/y(ii))
  else
     r = y(ii)/ypI ;
     if(  r .gt.0.  ) then
        phase=-( PI/2.0 - atan(r) )
     else
        phase=- ( - PI/2. - atan( r ));
     endif
  endif


  ! print *, " phase " 

  
  if( dabs(y(ii) ) .gt. dabs(ypI)) then
     y(ii+1)=y(ii+1)/y(ii)
     do i=1, ii
        y(i)=y(i)/y(ii)
     enddo
  else
     do i=1, ii+1
        y(i)=y(i)/ypI
     enddo
  endif
  
  ! //------------ Propagation de   I   à  rinf --------------------------
  
  do i=N,N-1,-1 
     y(i)=N-i    ! //y[i]=exp(-sqrt(-2*E)*r[i]);
  enddo
  
  do i=N-1, ii,-1
     dh = rgrid(i)-rgrid(i-1)
     dl = rgrid(i+1)-rgrid(i)
     G0 = 2*(V(i)  -E)
     Nlb=nonloc(i)
     
     if(dabs(G0*dh*dl).lt.1.0) then
        dh2= dh*dh
        dl2=dl*dl
        dp = dh*dl
        
        Gh = 2*(V(i-1)-E)
        Gl = 2*(V(i+1)-E)
        
        Nlc=nonloc(i-1)
        Nla=nonloc(i+1)
		  
        R  = dh2 + dl2 + dp 
		  
        dh3=dh2*dh
        dhdl2=dl*dp
        dh2dl=dh*dp
        dl3=dl*dl2
        
        ! // ********************************************************************************
        ! // ** 17 luglio 1998. Aggiunta della parte per il termine non locale.
        ! // ** Vedere file numerov.mathematica
        ! // ** ( h^3*(-Nla  + Nlb ) + l^3*(Nlb  - Nlc) + 
        ! // **      h^2*l*(Nla+4*Nlbd+Nlc)+h*l^2*(Nla+4*Nlb+Nlc))/l*(12-(h^2+h*l-l^2)*Gc))
        
        deno=(dl*(1.-Gh*(R-2*dl2)/12.))
        
        add1=( y(i)*(dh+(dl*Gh*(R-2*dl2)+G0*(dh+dl)*(R+2*dh*dl))/12)&
             -y(i+1)*dh*(1.-Gl*(R-2*dh2)/12.) &
             )/deno
        
        add2=( dh3*(-Nla  + Nlb ) + dl3*(Nlb  - Nlc) + &
             dh2dl*(Nla  + 4*Nlb + Nlc ) + dhdl2*(Nla + 4*Nlb + Nlc) )/deno
		  
        y(i-1) = y(i)+add1+add2
        
        ! /*  
        ! y[i-1] = y[i]+( y[i]*(dh+(dl*Gh*(R-2*dl2)+G0*(dh+dl)*(R+2*dh*dl))/12)
        ! -y[i+1]*dh*(1.-Gl*(R-2*dh2)/12.)
        ! )/
        ! (dl*(1.-Gh*(R-2*dl2)/12.));
        ! */
     else
        print *, "usare un pass piu piccolo "
        stop
     endif
 
     

     if(dabs(y(i-1)).gt.1.0D8) then 
        fact= 1./dabs(y(i-1)) 
        count =0
        do j=i-1, N 
           if(y(j).eq.0) then
              count=count+1
           else 
              count =0
           endif
           if(count.ne.0) exit
           y(j)=y(j)*fact 
        enddo
     endif

     if( (i-1) .ge.ii ) then!  // LA PHASE EST CALCULE EN I
        if((y(i)*y(i-1) ).lt.0) yNcross=yNcross+1
        if(yflag.eq.1 .and. (y(i)*y(i-1)).eq.0) then
           yflag=1-yflag
           yNcross=yNcross+1
        endif
     endif
  enddo
  !  // ypI=(y[I+1]*dh2-y[I-1]*dl2+y[I]*(dl2-dh2))/((dh+dl)*dl*dh);
  !  // sostituita con la maniera seguente, presa dal metodo numerov (vedi numerov.nb)
  ! {

  i=ii
  dh = rgrid(i+1)-rgrid(i)
  dl = rgrid(i  )-rgrid(i-1)
  
  Nlc=nonloc(i+1)
  Nla=nonloc(i-1)
  Nlb=nonloc(i)
  
  
  Gb = 2*(V(i)  -E)
  if(dabs(G0*dh*dl).gt.1.0) then
     print *, " problema fabs(G0*dh*dl)>1.0 in calcolo di yp in I "
     stop
  endif
  Gc = 2*(V(i+1)-E)
  Ga = 2*(V(i-1)-E)

  ya=y(ii-1)
  yb=y(ii)
  yc=y(ii+1)
	    
  ypI=(dh*(Ga*pow(dl,2)*(-12*dh + Gc*pow(dh,3) - 6*dl + &
       Gc*pow(dh,2)*dl) - &
       6*(-12*dh + Gc*pow(dh,3) - 12*dl + 2*Gc*pow(dh,2)*dl))*ya + &
       (dh + dl)*(Gb*pow(dl,2)*&
       (-24*dh + 2*Gc*pow(dh,3) - 6*dl + 3*Gc*pow(dh,2)*dl) + &
       6*(-12*dh + Gc*pow(dh,3) + Gc*pow(dh,2)*dl - &
       Gc*dh*pow(dl,2) + Gc*pow(dl,3)))*yb)/&
       (6.*dh*dl*(dh + dl)*(-12 + Gc*(pow(dh,2) + dh*dl - pow(dl,2)))) 
	    
  ! // *****************************************************************************
  !	    // ** Termine aggiuntivo
  !	    // **
  !	    // 6*l^3*Nlc + 
  !	    // Nla*(-12*h^2*l - 6*h*l^2 + h^3*l*(h + l)*Gc) + 
  !	    // Nlb*(-24*h^2*l - 30*h*l^2 - 6*l^3 +      2*h^3*l*(h + l)*Gc + 3*h^2*l^2*(h + l)*Gc)
  !	    // / al denominatore
  !	    //
  !	    //  6*h*(h + l)*(-12 + (h^2 + h*l - l^2)*Gc)
	    
  deno= 6*dh*(dh + dl)*(-12 + (dh*dh + dh*dl - dl*dl)*Gc)
	    
  num=  6*dl*dl*dl*Nlc &
       +Nla*(-12*dh*dh*dl-6*dh*dl*dl+dh*dh*dh*dl*(dh + dl)*Gc)&
       +Nlb*(-24*dh*dh*dl - 30*dh*dl*dl - &
       6*dl*dl*dl+2*dh*dh*dh*dl*(dh + dl)*Gc +3*dh*dh*dl*dl*(dh+dl)*Gc&
       )
	    
  ypI=ypI+num/deno
	    
  !}
	
  phase=phase+PI*yNcross 
	

 
  if(dabs(y(ii)) .gt. dabs( ypI ) ) Then
     phase =phase+atan(ypI/y(ii))
  else
     r = y(ii)/ypI 
     if(  r .gt.0.  ) then
        phase =phase+( PI/2. - atan(r) )
     else
        phase = phase+ ( - PI/2. - atan( r ))
     endif
  endif
	  
	  
  if( dabs(y(ii)) .gt. dabs(ypI)) then
     y(ii-1)=y(ii-1)/y(ii)
     do i=N,ii,-1
        y(i)=y(i)/y(ii)
     enddo
  else
     do i=N,ii-1,-1
        y(i)=y(i)/ypI
     enddo
  endif
	  
	  
  ! //------------------------------------------------------------------------------
  if(normalize.eq.1) then
     do i=1,N
        func(i)=y(i)*y(i)
     enddo
     
     call integrate(func,funcdum,rgrid,N)
     norm=sqrt(funcdum(N) )
     print *, " norm " , norm
     do i=1,N
        y(i)=y(i)/norm
     enddo
  endif
  return 
end function phase
	


subroutine schro(E, r,  V,nonloc, y, NGRID, nsol, l,  Z)
  implicit none
  integer nsol,ngrid
  real*8 E,r(ngrid),v(ngrid), nonloc(ngrid)
  real*8 y(ngrid)
  real*8 l,Z
  ! -------------------------------
  real*8 Elow, Ehigh, Eguess
  real*8 pathh, pathl, fase
  integer i
  real*8 Phase, but
  real*8 PI



  PI=4.0*atan(1.0)
  Elow=0;

  do i=1, NGRID
      if(V(i) .lt. Elow) Elow=V(i)
   enddo
  if( Elow .lt. -Z*Z) Elow = -Z*Z
  

  Ehigh = 0.0;


  
  pathh = Phase(Ehigh,NGRID,r,v,nonloc,y,  int(l) ,0, Z,1);
  pathl = Phase(Elow ,NGRID, r,v,nonloc,y,  int(l) ,0,  Z,1);
 
  if( pathl.gt.0.0)  then
      print *, " pathl>0.0 " 
      stop
   endif
  

      
   but= PI*(nsol-l-1)
   if(but .gt. pathh) then
      Ehigh = (but+1)*(but+1)/r(NGRID-1)/r(NGRID-1)
      do while( but .gt. Phase(Ehigh ,NGRID, r,V,nonloc,y, int(l) , 0, Z,1 )	 )
         Ehigh =2*Ehigh;
      enddo
   endif

   do while(1.gt.0) 

      Eguess = (Elow+Ehigh)/2

      ! print *, " Eguess ", Eguess

      fase=  Phase(Eguess ,NGRID, r,V,nonloc,y,int(l) , 0, Z ,0)
      
      ! print *, " fase " , fase



      if( fase.gt.but) then
         Ehigh=Eguess
      else
         Elow=Eguess
      endif

      if(dabs(but-fase).lt.1.0e-8) exit
   enddo

   
   fase  = Phase(Eguess,NGRID, r,v,nonloc,y, int(l),1 ,  Z  ,0)
   E=Eguess;
   print *, E, fase
   print *, y(1), y(2), y(Ngrid)



   return
 end subroutine schro











subroutine gatom_modified(rcov,rprb,lmax,lpx,noccmax,occup,&
                 zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,nintp,&
                 aeval,ng,psi,res,chrg,&
                 Nsol, Labs, Ngrid,Egrid,  rgrid , psigrid )
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
       res(noccmax,lmax+1),xp(0:ng),& 
       rgrid(Ngrid), psigrid(Ngrid, Nsol),psigrid_naked(Ngrid,Nsol), projgrid(Ngrid,3), &
       rhogrid(Ngrid), potgrid(Ngrid), &
       vxcgrid(Ngrid), &
       dumgrid1(Ngrid),dumgrid2(Ngrid),  &
       Egrid(nsol), ppgrid(Nsol,3), work(nsol*nsol*2), &
       H(Nsol, Nsol)

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
        
        if(l.eq.Labs) then
           do iorder=1,3
              do igrid=1,Ngrid
                 r=rgrid(igrid)
                 projgrid(igrid,iorder) = (  exp( -r*r*tt)*(r**l *r) )* r**(2*(i-1))
              end do
              
              do igrid=1, Ngrid
                 dumgrid1(igrid) =   projgrid(igrid,iorder)*projgrid(igrid,iorder)
              enddo
              call integrate( dumgrid1, dumgrid2, rgrid, Ngrid)
              do igrid=1, Ngrid
                 projgrid(igrid,iorder) = projgrid(igrid,iorder) / sqrt( dumgrid2(Ngrid) ) 
              enddo
           enddo
        endif
     end do
  else
     pp1(:,:)=0._gp
     pp2(:,:)=0._gp
     pp3(:,:)=0._gp
     projgrid(:,:)=0.0_gp
  end if

  do l=0,lmax
     do j=0,ng
        do i=0,ng
           rho(i,j,l+1)=0._gp
        end do
     end do
  end do


  do igrid=1,Ngrid
     rhogrid(igrid)=0.0_gp
  enddo

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

     do igrid=1,Ngrid
        r=rgrid(igrid)
        rhogrid(igrid)=0.0_gp
        do l=0, lmax
           do iocc=1,noccmax
              if (occup(iocc,l+1).gt.0._gp) then
                 dum=0.0_gp
                 do j=0,ng
                    dum=dum + psi(j,iocc,l+1)*exp(-r*r *   xp(j) )* r**(l+1)
                 end do
                 rhogrid(igrid)=rhogrid(igrid)+dum*dum*occup(iocc,l+1)
              end if
           end do
        enddo
        dum = rhogrid(igrid)/r/r *0.07957747154594768_gp
        vxcgrid(igrid)=emuxc(dum)
     enddo
     
  
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

     do igrid=1, Ngrid
        r=rgrid(igrid)
        potgrid(igrid) =0.5_gp*r*r  /    rprb**4 
        potgrid(igrid) = potgrid(igrid) - zion/r * derf( r/alpz/sqrt(2.0)   )
        rr = r/alpz
        potgrid(igrid) = potgrid(igrid) + exp(-0.5_gp * rr**2 )*( gpot(1)+gpot(2)*rr**2 + gpot(3)*rr**4 )
     enddo
! poisson per potgrid
     call integrate( rhogrid, dumgrid1, rgrid, Ngrid)
     do igrid=1, Ngrid
        potgrid(igrid)=potgrid(igrid)+dumgrid1(igrid)/rgrid(igrid)
     enddo

     do igrid=1, Ngrid
        dumgrid1(igrid) = rhogrid(igrid)/rgrid(igrid)
     enddo
     call integrate( dumgrid1, dumgrid2, rgrid, Ngrid)
     do igrid=1, Ngrid
        potgrid(igrid)=potgrid(igrid)-dumgrid2(igrid) + dumgrid2(Ngrid)
     enddo
! ------------------------------------------------------------------------------------



 

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
              if (1.eq.1) then
              hh(i,j)=hh(i,j)+ gpot(1)*.5_gp*gamma(1.5_gp+real(l,gp))*tt**(1.5_gp+real(l,gp))&
                   + (gpot(2)/alpl**2)*.5_gp*gamma(2.5_gp+real(l,gp))*tt**(2.5_gp+real(l,gp))&
                   + (gpot(3)/alpl**4)*.5_gp*gamma(3.5_gp+real(l,gp))*tt**(3.5_gp+real(l,gp))
           endif
! separable terms
              if (1.eq.1 .and. l.le.lpx) then
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



  dumgrid1(:)=0.0_gp
  do igrid=1, ngrid
     r=rgrid(igrid)
     potgrid(igrid)=potgrid(igrid)+ 0.5_gp*labs*(labs+1.0_gp)/r/r
  enddo
  


  do isol=1,nsol
     print *, " chiamo schro per isol=", isol
     call schro(Egrid(isol) , rgrid ,  potgrid , dumgrid1, psigrid_naked(:,isol) , ngrid , isol ,real( labs,gp) ,  zion)
     print *, " Egrid(", isol,")=", Egrid(isol)
  enddo
  

  H(:,:)=0.0D0
  do i=1,Nsol
     H(i,i)=Egrid(i)
     do iproj=1,3
        do igrid=1,Ngrid
           dumgrid1(igrid)=psigrid_naked(igrid,i)*projgrid(igrid,iproj)
        enddo
        call integrate(dumgrid1,dumgrid2,rgrid,Ngrid)
        ppgrid(i,iproj)=dumgrid2(Ngrid)
     enddo
  enddo

  do i=1,Nsol
     do j=1, Nsol
        H(i,j)=H(i,j)+ ppgrid(i,1)*hsep(1,labs+1)*ppgrid(j,1)&
                      + ppgrid(i,1)*hsep(2,labs+1)*ppgrid(j,2)&
                      + ppgrid(i,2)*hsep(2,labs+1)*ppgrid(j,1)&
                      + ppgrid(i,2)*hsep(3,labs+1)*ppgrid(j,2)&
                      + ppgrid(i,1)*hsep(4,labs+1)*ppgrid(j,3)&
                      + ppgrid(i,3)*hsep(4,labs+1)*ppgrid(j,1)&
                      + ppgrid(i,2)*hsep(5,labs+1)*ppgrid(j,3)&
                      + ppgrid(i,3)*hsep(5,labs+1)*ppgrid(j,2)&
                      + ppgrid(i,3)*hsep(6,labs+1)*ppgrid(j,3)

        do igrid=1,Ngrid
           dumgrid1(igrid)=psigrid_naked(igrid,i)*psigrid_naked(igrid,j)*vxcgrid(igrid)
        enddo
        call integrate(dumgrid1,dumgrid2,rgrid,Ngrid)
        H(i,j)=H(i,j)+dumgrid2(Ngrid)
     enddo
  enddo
  


  call DSYEV('V','U', Nsol, H, Nsol,Egrid , WORK, Nsol*Nsol*2, INFO)
  
  call  DGEMM('N','N',Ngrid ,Nsol,   Nsol,1.0d0 ,psigrid_naked, Ngrid ,H,Nsol, 0.0D0 , psigrid , Ngrid)

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


! ------------------------------------------------
  


! -----------------------------------------------




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

END SUBROUTINE gatom_modified



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
     gamma=1.0_gp
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
  real(kind=8), dimension(nvctr_c+7*nvctr_f,norbep*nspin), intent(inout) :: psi
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
