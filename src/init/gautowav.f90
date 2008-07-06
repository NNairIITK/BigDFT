subroutine dual_gaussian_coefficients(norbp,G,coeffs)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: norbp
  type(gaussian_basis), intent(in) :: G
  real(gp), dimension(G%ncoeff,norbp), intent(inout) :: coeffs !warning: the precision here should be wp
  !local variables
  character(len=*), parameter :: subname='dual_gaussian_coefficients'
  integer :: nwork,info,i_stat,i_all
  integer, dimension(:), allocatable :: iwork
  real(gp), dimension(:), allocatable :: ovrlp,work
  
  allocate(iwork(6+ndebug),stat=i_stat)
  call memocc(i_stat,iwork,'iwork',subname)
  allocate(ovrlp(G%ncoeff*G%ncoeff+ndebug),stat=i_stat)
  call memocc(i_stat,ovrlp,'ovrlp',subname)

  !temporary allocation of the work array, workspace query in dsysv
  allocate(work(1+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)

  call dsysv('U',G%ncoeff,norbp,ovrlp(1),G%ncoeff,iwork(1),coeffs(1,1),G%ncoeff,&
       work(1),-1,info)
  nwork=work(1)

  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)
  allocate(work(nwork+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)

  
  call gaussian_overlap(G,G,ovrlp)
  call dsysv('U',G%ncoeff,norbp,ovrlp(1),G%ncoeff,iwork(1),coeffs(1,1),G%ncoeff,&
       work,nwork,info)

  i_all=-product(shape(iwork))*kind(iwork)
  deallocate(iwork,stat=i_stat)
  call memocc(i_stat,i_all,'iwork',subname)
  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)
  i_all=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ovrlp',subname)

end subroutine dual_gaussian_coefficients

!control the accuracy of the expansion in gaussian
subroutine check_gaussian_expansion(geocode,iproc,nproc,norb,norbp,&
     n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hx,hy,hz,wfd,psi,G,coeffs)
  use module_base
  use module_types
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(gp), intent(in) :: hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(gaussian_basis), intent(in) :: G
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(in) :: psi
  real(wp), dimension(G%ncoeff,norbp), intent(in) :: coeffs
  !local variables
  character(len=*), parameter :: subname='check_gaussian_expansion'
  integer :: iorb,i_stat,i_all,i,j,ierr
  real(wp) :: maxdiffp,maxdiff,orbdiff
  real(wp), dimension(:), allocatable :: workpsi

  allocate(workpsi((wfd%nvctr_c+7*wfd%nvctr_f)*norbp+ndebug),stat=i_stat)
  call memocc(i_stat,workpsi,'workpsi',subname)

  call gaussians_to_wavelets(geocode,iproc,nproc,norb,norbp,&
     n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hx,hy,hz,wfd,G,coeffs,workpsi)

  maxdiffp=0.0_wp
  do iorb=1,norbp
     orbdiff=0.0_wp
     if (iorb+iproc*norbp <= norb) then
        do i=1,wfd%nvctr_c+7*wfd%nvctr_f
           j=i+(iorb-1)*wfd%nvctr_c+7*wfd%nvctr_f
           orbdiff=max(orbdiff,(psi(i,iorb)-workpsi(j))**2)
        end do
     end if
     maxdiffp=max(maxdiffp,orbdiff)
     print *,'iproc,iorb,orbdiff',iorb,orbdiff
  end do

  call MPI_GATHER(maxdiffp,1,mpidtypw,maxdiff,1,mpidtypw,0,MPI_COMM_WORLD,ierr)

  if (iproc == 0) then
     write(*,*)' Mean L2 norm of gaussian-wavelet difference:',sqrt(maxdiff)/real(norb,wp)
  end if
  i_all=-product(shape(workpsi))*kind(workpsi)
  deallocate(workpsi,stat=i_stat)
  call memocc(i_stat,i_all,'workpsi',subname)


end subroutine check_gaussian_expansion

subroutine parse_cp2k_files(iproc,basisfile,orbitalfile,nat,ntypes,norb,iatype,rxyz,CP2K,wfn_cp2k)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: basisfile,orbitalfile
  integer, intent(in) :: norb,iproc,nat,ntypes
  integer, dimension(nat), intent(in) :: iatype
  real(gp), dimension(3,nat), target, intent(in) :: rxyz
  type(gaussian_basis), intent(out) :: CP2K
  real(wp), dimension(:,:), pointer :: wfn_cp2k
  !local variables
  character(len=*), parameter :: subname='parse_cp2k_files'
  character(len=6) :: string,symbol
  character(len=100) :: line
  integer, parameter :: nterm_max=3
  integer :: ngx,nbx,npgf,nst,nend,ng,lshell,num,mmax,myshift,icbas,isbas,nbas,nco,i,ipar,ipg,jat
  integer :: iorb,jorb,iat,ityp,l,m,nterm,i_all,i_stat,ibas,ig,iset,jbas,iterm,ishell,lmax,m1,m2
  integer :: ierr,isat,iexpo,icoeff,iam
  real(dp) :: tt,normdev
  real(gp) :: exponent,coefficient
  integer, dimension(:), allocatable :: nshell,iorbtmp,iw
  integer, dimension(:,:), allocatable :: nam,ndoc
  real(gp), dimension(:), allocatable :: ctmp
  real(gp), dimension(:,:,:), allocatable :: contcoeff,expo
  real(wp), dimension(:,:,:,:), allocatable :: cimu

  !parse the output of CP2K to read the basis set information

  if (iproc==0) write(*,'(1x,a)',advance='no')&
       'Reading Basis Set information and wavefunctions coefficients...'

  ngx=0
  nbx=0
  lmax=0

  allocate(nshell(ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,nshell,'nshell',subname)

  open(unit=35,file=trim(basisfile),action='read')

  !read the lines for analyzing the maximum number of primitive gaussian functions
  !and also for the maximum angular momentum
  ityp=0
  ngx=0
  nbx=0
  lmax=0
  ipg=0
  for_ngx: do
     if (ityp > ntypes) exit for_ngx
     read(35,'(a100)')line
     !analyzing the different possibilities
     read(line,*,iostat=i_stat)tt,string,symbol
     if (i_stat == 0 .and. string=='Atomic' .and. symbol=='kind:') then
        ityp=ityp+1
        if (ityp > 1) then
           nshell(ityp-1)=ishell
           nbx=max(nbx,ishell)
        end if
        ishell=0
        cycle for_ngx
     end if
     read(line,*,iostat=i_stat)iset,num,num,num,exponent,coefficient
     if (i_stat==0) then
        !print *,num,exponent,coefficient
        ishell=ishell+1
        lmax=max(lmax,num)
        ngx=max(ngx,ipg)
        !print *,ishell,ipg,lmax
        ipg=1
        cycle for_ngx
     end if
     read(line,*,iostat=i_stat)exponent,coefficient
     if (i_stat==0) then
        ipg=ipg+1
        cycle for_ngx
     end if
  end do for_ngx

  !now store the values
  rewind(35)

  !here allocate arrays
  allocate(nam(nbx,ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,nam,'nam',subname)
  allocate(ndoc(nbx,ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,ndoc,'ndoc',subname)
  allocate(contcoeff(ngx,nbx,ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,contcoeff,'contcoeff',subname)
  allocate(expo(ngx,nbx,ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,expo,'expo',subname)
  
  ityp=0
  ipg=0
  store_basis: do
     if (ityp > ntypes) exit store_basis
     read(35,'(a100)')line
     !analyzing the different possibilities
     read(line,*,iostat=i_stat)tt,string,symbol
     if (i_stat == 0 .and. string=='Atomic' .and. symbol=='kind:') then
        ityp=ityp+1
        if (ityp > 1) then 
           ndoc(ishell,ityp-1)=ipg
        end if
        ishell=0
        cycle store_basis
     end if
     read(line,*,iostat=i_stat)iset,num,num,num,exponent,coefficient
     if (i_stat==0) then
        !print *,num,exponent,coefficient
        ishell=ishell+1
        nam(ishell,ityp)=num
        lmax=max(lmax,num)
        if (ishell > 1) ndoc(ishell-1,ityp)=ipg
        expo(1,ishell,ityp)=exponent
        contcoeff(1,ishell,ityp)=coefficient
        ipg=1
        cycle store_basis
     end if
     read(line,*,iostat=i_stat)exponent,coefficient
     if (i_stat==0) then
        ipg=ipg+1
        expo(ipg,ishell,ityp)=exponent
        contcoeff(ipg,ishell,ityp)=coefficient
        cycle store_basis
     end if
  end do store_basis

  !close the file of the basis definition
  close(35)

  !renormalize the coefficients in each shell
  do ityp=1,ntypes
     do ishell=1,nshell(ityp)
        call normalize_shell(ndoc(ishell,ityp),nam(ishell,ityp),&
             expo(1,ishell,ityp),contcoeff(1,ishell,ityp))
     end do
  end do


  !the number of gaussian centers are thus nat
  CP2K%nat=nat
  CP2K%rxyz => rxyz
  !copy the parsed values in the gaussian structure
  !count also the total number of shells
  allocate(CP2K%nshell(nat+ndebug),stat=i_stat)
  call memocc(i_stat,CP2K%nshell,'CP2K%nshell',subname)
  
  CP2K%nshltot=0
  do iat=1,nat
     ityp=iatype(iat)
     CP2K%nshell(iat)=nshell(ityp)
     CP2K%nshltot=CP2K%nshltot+nshell(ityp)
  end do

  allocate(CP2K%ndoc(CP2K%nshltot+ndebug),stat=i_stat)
  call memocc(i_stat,CP2K%ndoc,'CP2K%ndoc',subname)
  allocate(CP2K%nam(CP2K%nshltot+ndebug),stat=i_stat)
  call memocc(i_stat,CP2K%nam,'CP2K%nam',subname)

  !assign shell IDs and count the number of exponents and coefficients
  CP2K%nexpo=0
  CP2K%ncoeff=0
  ishell=0
  do iat=1,nat
     ityp=iatype(iat)
     do isat=1,CP2K%nshell(iat)
        ishell=ishell+1
        CP2K%ndoc(ishell)=ndoc(isat,ityp)
        CP2K%nam(ishell)=nam(isat,ityp)+1
        CP2K%nexpo=CP2K%nexpo+ndoc(isat,ityp)
        CP2K%ncoeff=CP2K%ncoeff+2*nam(isat,ityp)+1
     end do
  end do

  !allocate and assign the exponents and the coefficients
  allocate(CP2K%xp(CP2K%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,CP2K%xp,'CP2K%xp',subname)
  allocate(CP2K%psiat(CP2K%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,CP2K%psiat,'CP2K%psiat',subname)

  ishell=0
  iexpo=0
  do iat=1,nat
     ityp=iatype(iat)
     do isat=1,CP2K%nshell(iat)
        ishell=ishell+1
        do ig=1,CP2K%ndoc(ishell)
           iexpo=iexpo+1
           CP2K%psiat(iexpo)=contcoeff(ig,isat,ityp)
           CP2K%xp(iexpo)=sqrt(0.5_gp/expo(ig,isat,ityp))
        end do
     end do
  end do


!!$  !print the found values
!!$  do ityp=1,ntypes
!!$     do ishell=1,nshell(ityp)
!!$        print *,'ityp=',ityp,'ishell=',ishell,'l=',nam(ishell,ityp),'ndoc=',ndoc(ishell,ityp)
!!$        do ipg=1,ndoc(ishell,ityp)
!!$           print *,'expo=',expo(ipg,ishell,ityp),'coeff=',contcoeff(ipg,ishell,ityp)
!!$        end do
!!$     end do
!!$  end do


  i_all=-product(shape(contcoeff))*kind(contcoeff)
  deallocate(contcoeff,stat=i_stat)
  call memocc(i_stat,i_all,'contcoeff',subname)
  i_all=-product(shape(expo))*kind(expo)
  deallocate(expo,stat=i_stat)
  call memocc(i_stat,i_all,'expo',subname)


  mmax=2*lmax+1
  !now read the coefficients of the gaussian converged orbitals
  open(unit=36,file=trim(orbitalfile),action='read')
  !here there is the orbital label, for the moment it is assumed to vary between 1 and 4
  allocate(ctmp(10+ndebug),stat=i_stat)
  call memocc(i_stat,ctmp,'ctmp',subname)
  allocate(iorbtmp(10+ndebug),stat=i_stat)
  call memocc(i_stat,iorbtmp,'iorbtmp',subname)
  allocate(cimu(mmax,nbx,nat,norb+ndebug),stat=i_stat)
  call memocc(i_stat,cimu,'cimu',subname)

  read(36,*)
  read_line1: do
     read(36,'(a100)')line
     !analyse how many orbitals are contained in a given line
     read_orbitals1: do ipar=10,1,-1
        read(line,*,iostat=i_stat)(iorbtmp(i),i=1,ipar)
        if (i_stat==0) then
           read(line,*)nst
           exit read_line1
        end if
     end do read_orbitals1
  end do read_line1
  nend=nst+ipar-1

!!$  nst=1
!!$  nend=4
  jat=0
  ishell=1
  jbas=0
  iat=nat
  !now read the data to assign the coefficients
  store_coeff: do
     read(36,'(a100)')line
     !choose between different cases
     read(line,*,iostat=i_stat)ibas,iat,symbol,string,(ctmp(iorb),iorb=1,nend-nst+1)
     if (i_stat==0) then
        !print *,line,nst,nend
        if (jat==iat) then
           jbas=jbas+1
           if (jbas > 2*nam(ishell,iatype(iat))+1) then
              jbas=1
              ishell=ishell+1
              if (ishell > nshell(iatype(iat))) then
                 if (iproc==0) write(*,'(1x,a,i0,a)')&
                      'Problem in the gaucoeff.dat file, the number of shells of atom ',iat ,&
                      ' is incoherent'
                 stop
              end if
           end if
        else
           jbas=1
           ishell=1
        end if
       symbol=trim(string)
       do iorb=nst,nend
          cimu(jbas+myshift(symbol),ishell,iat,iorb)=ctmp(iorb-nst+1)
       end do
       jat=iat
       if (jbas==2*nam(ishell,iatype(iat))+1 .and. ishell==nshell(iatype(iat))&
            .and. iat==nat .and. nend==norb) then
          exit store_coeff
       else
          cycle store_coeff
       end if
     end if

     read_orbitals: do ipar=10,1,-1
        read(line,*,iostat=i_stat)(iorbtmp(i),i=1,ipar)
        if (i_stat==0) then
           read(line,*)nst
           nend=nst+ipar-1
           if (jat/=nat) then
              if (iproc==0)write(*,'(1x,a,i0,a)')&
                   'Problem in the gaucoeff.dat file, only ',iat ,' atoms processed'
              stop
           else
              cycle store_coeff
           end if
        end if
     end do read_orbitals

  end do store_coeff
  close(36)

!!$  !print the found values
!!$  do iat=1,nat
!!$     ityp=iatype(iat)
!!$     do ishell=1,nshell(ityp)
!!$        do jbas=1,2*nam(ishell,ityp)+1
!!$           print *,iat,ishell,nam(ishell,ityp),jbas,(cimu(jbas,ishell,iat,iorb),iorb=1,norb)
!!$        end do
!!$     end do
!!$  end do

  !allocate and assign the coefficients of each orbital
  allocate(wfn_cp2k(CP2K%ncoeff,norb+ndebug),stat=i_stat)
  call memocc(i_stat,wfn_cp2k,'wfn_cp2k',subname)
  do iorb=1,norb
     icoeff=0
     ishell=0
     do iat=1,CP2K%nat
        do isat=1,CP2K%nshell(iat)
           ishell=ishell+1
           do iam=1,2*CP2K%nam(ishell)-1
              icoeff=icoeff+1
              wfn_cp2k(icoeff,iorb)=cimu(iam,isat,iat,iorb)
           end do
        end do
     end do
  end do

  call gaudim_check(1,icoeff+1,ishell,0,CP2K%ncoeff,CP2K%nshltot)

  i_all=-product(shape(ctmp))*kind(ctmp)
  deallocate(ctmp,stat=i_stat)
  call memocc(i_stat,i_all,'ctmp',subname)
  i_all=-product(shape(iorbtmp))*kind(iorbtmp)
  deallocate(iorbtmp,stat=i_stat)
  call memocc(i_stat,i_all,'iorbtmp',subname)
  i_all=-product(shape(nshell))*kind(nshell)
  deallocate(nshell,stat=i_stat)
  call memocc(i_stat,i_all,'nshell',subname)
  i_all=-product(shape(nam))*kind(nam)
  deallocate(nam,stat=i_stat)
  call memocc(i_stat,i_all,'nam',subname)
  i_all=-product(shape(ndoc))*kind(ndoc)
  deallocate(ndoc,stat=i_stat)
  call memocc(i_stat,i_all,'ndoc',subname)
  i_all=-product(shape(cimu))*kind(cimu)
  deallocate(cimu,stat=i_stat)
  call memocc(i_stat,i_all,'cimu',subname)

  if (iproc==0) then
     write(*,'(1x,a)')'done.'
  end if


end subroutine parse_cp2k_files

subroutine gaussians_to_wavelets(geocode,iproc,nproc,norb,norbp,&
     n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hx,hy,hz,wfd,G,wfn_gau,psi)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3 
  real(gp), intent(in) :: hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(gaussian_basis), intent(in) :: G
  real(wp), dimension(G%ncoeff,norb), intent(in) :: wfn_gau
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(out) :: psi

  !local variables
  character(len=*), parameter :: subname='gaussians_to_wavelets'
  integer, parameter :: nterm_max=3
  logical :: myorbital
  integer :: i_stat,i_all,ishell,iexpo,icoeff,iat,isat,ng,l,m,iorb,jorb,i,nterm,ierr,ig
  real(dp) :: normdev,tt,scpr
  real(gp) :: rx,ry,rz
  integer, dimension(nterm_max) :: lx,ly,lz
  real(gp), dimension(nterm_max) :: fac_arr
  real(wp), dimension(:), allocatable :: tpsi

  if(iproc == 0) write(*,'(1x,a)',advance='no')'Writing wavefunctions in wavelet form '

  allocate(tpsi(wfd%nvctr_c+7*wfd%nvctr_f+ndebug),stat=i_stat)
  call memocc(i_stat,tpsi,'tpsi',subname)

  !initialize the wavefunction
  call razero((wfd%nvctr_c+7*wfd%nvctr_f)*norbp,psi)
  !this can be changed to be passed only once to all the gaussian basis
  !eks=0.d0
  !loop over the atoms
  ishell=0
  iexpo=1
  icoeff=1
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
        !multiply the values of the gaussian contraction times the orbital coefficient
        do m=1,2*l-1
           call calc_coeff_inguess(l,m,nterm_max,nterm,lx,ly,lz,fac_arr)
!!$           !this kinetic energy is not reliable
!!$           eks=eks+ek*occup(iorb)*cimu(m,ishell,iat,iorb)
           call crtonewave(geocode,n1,n2,n3,ng,nterm,lx,ly,lz,fac_arr,G%xp(iexpo),G%psiat(iexpo),&
                rx,ry,rz,hx,hy,hz,0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                wfd%nseg_c,wfd%nvctr_c,wfd%keyg,wfd%keyv,wfd%nseg_f,wfd%nvctr_f,&
                wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1),&
                tpsi(1),tpsi(wfd%nvctr_c+1))
           !sum the result inside the orbital wavefunction
           !loop over the orbitals
           do iorb=1,norb
              if (myorbital(iorb,norb,iproc,nproc)) then
                 jorb=iorb-iproc*norbp
                 do i=1,wfd%nvctr_c+7*wfd%nvctr_f
                    !for this also daxpy BLAS can be used
                    psi(i,jorb)=psi(i,jorb)+wfn_gau(icoeff,iorb)*tpsi(i)
                 end do
              end if
           end do
           icoeff=icoeff+1
        end do
        iexpo=iexpo+ng
     end do
     if (iproc == 0) then
        write(*,'(a)',advance='no') &
             repeat('.',(iat*40)/G%nat-((iat-1)*40)/G%nat)
     end if
  end do

  call gaudim_check(iexpo,icoeff,ishell,G%nexpo,G%ncoeff,G%nshltot)

  if (iproc ==0 ) write(*,'(1x,a)')'done.'
  !renormalize the orbitals
  !calculate the deviation from 1 of the orbital norm
  normdev=0.0_dp
  tt=0.0_dp
  do iorb=1,norb
     if (myorbital(iorb,norb,iproc,nproc)) then
        jorb=iorb-iproc*norbp
        call wnrm(wfd%nvctr_c,wfd%nvctr_f,psi(1,jorb),psi(wfd%nvctr_c+1,jorb),scpr) 
        call wscal(wfd%nvctr_c,wfd%nvctr_f,real(1.0_dp/sqrt(scpr),wp),psi(1,jorb),psi(wfd%nvctr_c+1,jorb))
        print *,'norm of orbital ',iorb,scpr
        tt=max(tt,abs(1.0_dp-scpr))
     end if
  end do
  if (nproc > 1) then
     call MPI_REDUCE(tt,normdev,1,mpidtypd,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  else
     normdev=tt
  end if
  if (iproc ==0 ) write(*,'(1x,a,1pe12.2)')&
       'Deviation from normalization of the imported orbitals',normdev

  i_all=-product(shape(tpsi))*kind(tpsi)
  deallocate(tpsi,stat=i_stat)
  call memocc(i_stat,i_all,'tpsi',subname)


end subroutine gaussians_to_wavelets


subroutine gautowav(geocode,iproc,nproc,nat,ntypes,norb,norbp,n1,n2,n3,&
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     nvctr_c,nvctr_f,nseg_c,nseg_f,keyg,keyv,iatype,occup,rxyz,hx,hy,hz,psi,eks)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: norb,norbp,iproc,nproc,nat,ntypes
  integer, intent(in) :: nvctr_c,nvctr_f,n1,n2,n3,nseg_c,nseg_f
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nat), intent(in) :: iatype
  real(gp), intent(in) :: hx,hy,hz
  real(gp), dimension(3,nat), intent(in) :: rxyz
  real(gp), dimension(norb), intent(in) :: occup
  real(gp), intent(out) :: eks
  real(wp), dimension(nvctr_c+7*nvctr_f,norbp), intent(out) :: psi
  !local variables
  character(len=*), parameter :: subname='gautowav'
  logical :: myorbital
  character(len=6) :: string,symbol
  character(len=100) :: line
  integer, parameter :: nterm_max=3
  integer :: ngx,nbx,npgf,nst,nend,ng,lshell,num,mmax,myshift,icbas,isbas,nbas,nco,i,ipar,ipg,jat
  integer :: iorb,jorb,iat,ityp,l,m,nterm,i_all,i_stat,ibas,ig,iset,jbas,iterm,ishell,lmax,m1,m2
  integer :: ierr
  real(dp) :: tt,normdev
  real(gp) :: rx,ry,rz
  real(gp) :: exponent,coefficient,scpr,ek
  integer, dimension(nterm_max) :: lx,ly,lz
  real(gp), dimension(nterm_max) :: fac_arr
  integer, dimension(:), allocatable :: nshell,iorbtmp,iw
  integer, dimension(:,:), allocatable :: nam,ndoc
  real(wp), dimension(:), allocatable :: tpsi,ctmp
  real(gp), dimension(:), allocatable :: psiatn,xp
  real(gp), dimension(:,:,:), allocatable :: contcoeff,expo
  real(wp), dimension(:,:,:,:), allocatable :: cimu

  !parse the output of CP2K to read the basis set information

  if (iproc==0) write(*,'(1x,a)',advance='no')&
       'Reading Basis Set information and wavefunctions coefficients...'

  ngx=0
  nbx=0
  lmax=0

  allocate(nshell(ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,nshell,'nshell',subname)

  open(unit=35,file='gaubasis.dat',action='read')

  !read the lines for analyzing the maximum number of primitive gaussian functions
  !and als for the maximum angular momentum
  ityp=0
  ngx=0
  nbx=0
  lmax=0
  ipg=0
  for_ngx: do
     if (ityp > ntypes) exit for_ngx
     read(35,'(a100)')line
     !analyzing the different possibilities
     read(line,*,iostat=i_stat)tt,string,symbol
     if (i_stat == 0 .and. string=='Atomic' .and. symbol=='kind:') then
        ityp=ityp+1
        if (ityp > 1) then
           nshell(ityp-1)=ishell
           nbx=max(nbx,ishell)
        end if
        ishell=0
        cycle for_ngx
     end if
     read(line,*,iostat=i_stat)iset,num,num,num,exponent,coefficient
     if (i_stat==0) then
        !print *,num,exponent,coefficient
        ishell=ishell+1
        lmax=max(lmax,num)
        ngx=max(ngx,ipg)
        !print *,ishell,ipg,lmax
        ipg=1
        cycle for_ngx
     end if
     read(line,*,iostat=i_stat)exponent,coefficient
     if (i_stat==0) then
        ipg=ipg+1
        cycle for_ngx
     end if
  end do for_ngx

  !now store the values
  rewind(35)

  !here allocate arrays
  allocate(nam(nbx,ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,nam,'nam',subname)
  allocate(ndoc(nbx,ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,ndoc,'ndoc',subname)
  allocate(contcoeff(ngx,nbx,ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,contcoeff,'contcoeff',subname)
  allocate(expo(ngx,nbx,ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,expo,'expo',subname)
  
  ityp=0
  ipg=0
  store_basis: do
     if (ityp > ntypes) exit store_basis
     read(35,'(a100)')line
     !analyzing the different possibilities
     read(line,*,iostat=i_stat)tt,string,symbol
     if (i_stat == 0 .and. string=='Atomic' .and. symbol=='kind:') then
        ityp=ityp+1
        if (ityp > 1) then 
           ndoc(ishell,ityp-1)=ipg
        end if
        ishell=0
        cycle store_basis
     end if
     read(line,*,iostat=i_stat)iset,num,num,num,exponent,coefficient
     if (i_stat==0) then
        !print *,num,exponent,coefficient
        ishell=ishell+1
        nam(ishell,ityp)=num
        lmax=max(lmax,num)
        if (ishell > 1) ndoc(ishell-1,ityp)=ipg
        expo(1,ishell,ityp)=exponent
        contcoeff(1,ishell,ityp)=coefficient
        ipg=1
        cycle store_basis
     end if
     read(line,*,iostat=i_stat)exponent,coefficient
     if (i_stat==0) then
        ipg=ipg+1
        expo(ipg,ishell,ityp)=exponent
        contcoeff(ipg,ishell,ityp)=coefficient
        cycle store_basis
     end if
  end do store_basis

  !close the file of the basis definition
  close(35)

  !renormalize the coefficients in each shell
  do ityp=1,ntypes
     do ishell=1,nshell(ityp)
        call normalize_shell(ndoc(ishell,ityp),nam(ishell,ityp),&
             expo(1,ishell,ityp),contcoeff(1,ishell,ityp))
     end do
  end do


!!$  !print the found values
!!$  do ityp=1,ntypes
!!$     do ishell=1,nshell(ityp)
!!$        print *,'ityp=',ityp,'ishell=',ishell,'l=',nam(ishell,ityp),'ndoc=',ndoc(ishell,ityp)
!!$        do ipg=1,ndoc(ishell,ityp)
!!$           print *,'expo=',expo(ipg,ishell,ityp),'coeff=',contcoeff(ipg,ishell,ityp)
!!$        end do
!!$     end do
!!$  end do

!!$  !here we can start calculate the overlap matrix between the different basis functions
!!$  !as an example we can try to calculate the overlap in one shell
!!$  allocate(iw(18),stat=i_stat)
!!$  call memocc(i_stat,product(shape(iw))*kind(iw),'iw','gautowav')
!!$  allocate(rw(6),stat=i_stat)
!!$  call memocc(i_stat,product(shape(rw))*kind(rw),'rw','gautowav')
!!$
!!$  do ityp=1,ntypes
!!$     do ishell=1,nshell(ityp)
!!$        !perform the scalar product internally to the shell
!!$        do m1=1,2*nam(ishell,ityp)+1
!!$           do m2=1,2*nam(ishell,ityp)+1
!!$              call gbasovrlp(expo(1,ishell,ityp),contcoeff(1,ishell,ityp),&
!!$                   expo(1,ishell,ityp),contcoeff(1,ishell,ityp),&
!!$                   ndoc(ishell,ityp),ndoc(ishell,ityp),&
!!$                   nam(ishell,ityp)+1,m1,nam(ishell,ityp)+1,m2,&
!!$                   0.d0,0.d0,0.d0,&
!!$                   18,6,iw,rw,ovrlp)
!!$              if (iproc==0) then
!!$                 print *,ityp,ishell,nam(ishell,ityp),m1,m2,ovrlp!&
!!$                      !contcoeff(1:ndoc(ishell,ityp),ishell,ityp),ovrlp
!!$              end if
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  i_all=-product(shape(iw))*kind(iw)
!!$  deallocate(iw,stat=i_stat)
!!$  call memocc(i_stat,i_all,'iw','gautowav')
!!$  i_all=-product(shape(rw))*kind(rw)
!!$  deallocate(rw,stat=i_stat)
!!$  call memocc(i_stat,i_all,'rw','gautowav')

!!$subroutine basis_ovrlp(nat,norb,nbx,ngx,lmax,ntypes,nam,ndoc,contcoeff,expo,cimu)
!!$  
!!$  
!!$end subroutine basis_ovrlp



  mmax=2*lmax+1
  !now read the coefficients of the gaussian converged orbitals
  open(unit=36,file='gaucoeff.dat',action='read')
  !here there is the orbital label, for the moment it is assumed to vary between 1 and 4
  allocate(ctmp(10+ndebug),stat=i_stat)
  call memocc(i_stat,ctmp,'ctmp',subname)
  allocate(iorbtmp(10+ndebug),stat=i_stat)
  call memocc(i_stat,iorbtmp,'iorbtmp',subname)
  allocate(cimu(mmax,nbx,nat,norb+ndebug),stat=i_stat)
  call memocc(i_stat,cimu,'cimu',subname)

  read(36,*)
  read_line1: do
     read(36,'(a100)')line
     !analyse how many orbitals are contained in a given line
     read_orbitals1: do ipar=10,1,-1
        read(line,*,iostat=i_stat)(iorbtmp(i),i=1,ipar)
        if (i_stat==0) then
           read(line,*)nst
           exit read_line1
        end if
     end do read_orbitals1
  end do read_line1
  nend=nst+ipar-1

!!$  nst=1
!!$  nend=4
  jat=0
  ishell=1
  jbas=0
  iat=nat
  !now read the data to assign the coefficients
  store_coeff: do
     read(36,'(a100)')line
     !choose between different cases
     read(line,*,iostat=i_stat)ibas,iat,symbol,string,(ctmp(iorb),iorb=1,nend-nst+1)
     if (i_stat==0) then
        !print *,line,nst,nend
        if (jat==iat) then
           jbas=jbas+1
           if (jbas > 2*nam(ishell,iatype(iat))+1) then
              jbas=1
              ishell=ishell+1
              if (ishell > nshell(iatype(iat))) then
                 if (iproc==0) write(*,'(1x,a,i0,a)')&
                      'Problem in the gaucoeff.dat file, the number of shells of atom ',iat ,&
                      ' is incoherent'
                 stop
              end if
           end if
        else
           jbas=1
           ishell=1
        end if
       symbol=trim(string)
       do iorb=nst,nend
          cimu(jbas+myshift(symbol),ishell,iat,iorb)=ctmp(iorb-nst+1)
       end do
       jat=iat
       if (jbas==2*nam(ishell,iatype(iat))+1 .and. ishell==nshell(iatype(iat))&
            .and. iat==nat .and. nend==norb) then
          exit store_coeff
       else
          cycle store_coeff
       end if
     end if

     read_orbitals: do ipar=10,1,-1
        read(line,*,iostat=i_stat)(iorbtmp(i),i=1,ipar)
        if (i_stat==0) then
           read(line,*)nst
           nend=nst+ipar-1
           if (jat/=nat) then
              if (iproc==0)write(*,'(1x,a,i0,a)')&
                   'Problem in the gaucoeff.dat file, only ',iat ,' atoms processed'
              stop
           else
              cycle store_coeff
           end if
        end if
     end do read_orbitals

  end do store_coeff
  close(36)

!!$  !print the found values
!!$  do iat=1,nat
!!$     ityp=iatype(iat)
!!$     do ishell=1,nshell(ityp)
!!$        do jbas=1,2*nam(ishell,ityp)+1
!!$           print *,iat,ishell,nam(ishell,ityp),jbas,(cimu(jbas,ishell,iat,iorb),iorb=1,norb)
!!$        end do
!!$     end do
!!$  end do

  i_all=-product(shape(ctmp))*kind(ctmp)
  deallocate(ctmp,stat=i_stat)
  call memocc(i_stat,i_all,'ctmp',subname)
  i_all=-product(shape(iorbtmp))*kind(iorbtmp)
  deallocate(iorbtmp,stat=i_stat)
  call memocc(i_stat,i_all,'iorbtmp',subname)


  !now apply this basis set information to construct the wavelets wavefunctions

  if (iproc==0) then
     write(*,'(1x,a)')'done.'
     write(*,'(1x,a)',advance='no')'Writing wavefunctions in wavelet form '
  end if

  allocate(psiatn(ngx+ndebug),stat=i_stat)
  call memocc(i_stat,psiatn,'psiatn',subname)
  allocate(xp(ngx+ndebug),stat=i_stat)
  call memocc(i_stat,xp,'xp',subname)
  allocate(tpsi(nvctr_c+7*nvctr_f+ndebug),stat=i_stat)
  call memocc(i_stat,tpsi,'tpsi',subname)

  !initialize the wavefunction
  call razero((nvctr_c+7*nvctr_f)*norbp,psi)
  !this can be changed to be passed only once to all the gaussian basis
  !eks=0.d0
  !loop over the atoms
  do iat=1,nat
     ityp=iatype(iat)
     rx=rxyz(1,iat)
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)
     !loop over the number of shells of the atom type
     do ishell=1,nshell(ityp)
        !the degree of contraction of the basis function
        !is the same as the ng value of the createAtomicOrbitals routine
        ng=ndoc(ishell,ityp)
        !angular momentum of the basis set(shifted for compatibility with BigDFT routines
        l=nam(ishell,ityp)+1
        !amplitude coefficients (contraction coefficients of the basis times
        !the amplitude of this basis in that orbital)
        !exponents for the gaussian expansion adapted following the convention
        !of the routine gauss_to_daub
        do ig=1,ng
           psiatn(ig)=contcoeff(ig,ishell,ityp)
           xp(ig)=sqrt(0.5_gp/expo(ig,ishell,ityp))
        end do
        !multiply the values of the gaussian contraction times the orbital coefficient
        do m=1,2*l-1
           call calc_coeff_inguess(l,m,nterm_max,nterm,lx,ly,lz,fac_arr)
!!$           !this kinetic energy is not reliable
!!$           eks=eks+ek*occup(iorb)*cimu(m,ishell,iat,iorb)
           call crtonewave(geocode,n1,n2,n3,ng,nterm,lx,ly,lz,fac_arr,xp,psiatn,&
                rx,ry,rz,hx,hy,hz,0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,&
                keyg(1,nseg_c+1),keyv(nseg_c+1),&
                tpsi(1),tpsi(nvctr_c+1))
           !sum the result inside the orbital wavefunction
           !loop over the orbitals
           do iorb=1,norb
              if (myorbital(iorb,norb,iproc,nproc)) then
                 jorb=iorb-iproc*norbp
                 do i=1,nvctr_c+7*nvctr_f
                    !for this also daxpy BLAS can be used
                    psi(i,jorb)=psi(i,jorb)+cimu(m,ishell,iat,iorb)*tpsi(i)
                 end do
              end if
           end do
        end do
     end do
     if (iproc == 0) then
        write(*,'(a)',advance='no') &
             repeat('.',(iat*40)/nat-((iat-1)*40)/nat)
     end if
  end do
  if (iproc ==0 ) write(*,'(1x,a)')'done.'
  !renormalize the orbitals
  !calculate the deviation from 1 of the orbital norm
  normdev=0.0_dp
  tt=0.0_dp
  do iorb=1,norb
     if (myorbital(iorb,norb,iproc,nproc)) then
        jorb=iorb-iproc*norbp
        call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
        call wscal(nvctr_c,nvctr_f,real(1.0_dp/sqrt(scpr),wp),psi(1,jorb),psi(nvctr_c+1,jorb))
        !print *,'norm of orbital ',iorb,scpr
        tt=max(tt,abs(1.0_dp-scpr))
     end if
  end do
  if (nproc > 1) then
     call MPI_REDUCE(tt,normdev,1,mpidtypd,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  else
     normdev=tt
  end if
  if (iproc ==0 ) write(*,'(1x,a,1pe12.2)')&
       'Deviation from normalization of the imported orbitals',normdev

  !now we have to evaluate the eigenvalues of this hamiltonian

  i_all=-product(shape(tpsi))*kind(tpsi)
  deallocate(tpsi,stat=i_stat)
  call memocc(i_stat,i_all,'tpsi',subname)

  i_all=-product(shape(nshell))*kind(nshell)
  deallocate(nshell,stat=i_stat)
  call memocc(i_stat,i_all,'nshell',subname)
  i_all=-product(shape(nam))*kind(nam)
  deallocate(nam,stat=i_stat)
  call memocc(i_stat,i_all,'nam',subname)
  i_all=-product(shape(ndoc))*kind(ndoc)
  deallocate(ndoc,stat=i_stat)
  call memocc(i_stat,i_all,'ndoc',subname)
  i_all=-product(shape(contcoeff))*kind(contcoeff)
  deallocate(contcoeff,stat=i_stat)
  call memocc(i_stat,i_all,'contcoeff',subname)
  i_all=-product(shape(expo))*kind(expo)
  deallocate(expo,stat=i_stat)
  call memocc(i_stat,i_all,'expo',subname)
  i_all=-product(shape(cimu))*kind(cimu)
  deallocate(cimu,stat=i_stat)
  call memocc(i_stat,i_all,'cimu',subname)

  i_all=-product(shape(xp))*kind(xp)
  deallocate(xp,stat=i_stat)
  call memocc(i_stat,i_all,'xp',subname)
  i_all=-product(shape(psiatn))*kind(psiatn)
  deallocate(psiatn,stat=i_stat)
  call memocc(i_stat,i_all,'psiatn',subname)

end subroutine gautowav

!calculate the shift between the spherical harmonics of CP2K and the one of BigDFT
function myshift(symbol)
  implicit none
  character(len=5), intent(in) :: symbol
  integer :: myshift
  if (symbol(2:2)=='s') then
     myshift=0
  else if (symbol(2:2)=='p') then
     if( symbol(3:3)=='y') then
        myshift=1
     else if( symbol(3:3)=='z') then
        myshift=1
     else if( symbol(3:3)=='x') then
        myshift=-2
     end if
  else if ( symbol(2:2)=='d') then
     if( symbol(3:4)=='-2') then
        myshift=2
     else if( symbol(3:4)=='-1') then
        myshift=-1
     else if( symbol(3:3)=='0') then
        myshift=2
     else if( symbol(3:4)=='+1') then
        myshift=-2
     else if( symbol(3:4)=='+2') then
        myshift=-1
     end if
  else if ( symbol(2:2)=='f') then
     if( symbol(3:4)=='-3') then
        myshift=4
     else if( symbol(3:4)=='-2') then
        myshift=5
     else if( symbol(3:4)=='-1') then
        myshift=-1
     else if( symbol(3:3)=='0') then
        myshift=-1
     else if( symbol(3:4)=='+1') then
        myshift=-4
     else if( symbol(3:4)=='+2') then
        myshift=0
     else if( symbol(3:4)=='+3') then
        myshift=-3
     end if
  else if ( symbol(2:2)=='g') then
     write(*,'(1x,a)')'the orbitals of type g are not yet implemented in BigDFT'
     stop
  end if

end function myshift

logical function myorbital(iorb,norbe,iproc,nproc)
  implicit none
  integer, intent(in) :: iorb,norbe,iproc,nproc
  !local variables
  real(kind=8), parameter :: eps_mach=1.d-12
  integer :: norbep
  real(kind=8) :: tt

  tt=dble(norbe)/dble(nproc)
  norbep=int((1.d0-eps_mach*tt) + tt)
  if (iorb >= iproc*norbep+1 .and. iorb <= min((iproc+1)*norbep,norbe)) then
     myorbital=.true.
  else
     myorbital=.false.
  endif

end function myorbital

! returns an input guess orbital that is a Gaussian centered at a Wannier center
! exp (-1/(2*gau_a^2) *((x-cntrx)^2 + (y-cntry)^2 + (z-cntrz)^2 ))
! in the arrays psi_c, psi_f
subroutine crtonewave(geocode,n1,n2,n3,nterm,ntp,lx,ly,lz,fac_arr,xp,psiat,rx,ry,rz,hx,hy,hz, & 
     nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,  & 
     nseg_c,mvctr_c,keyg_c,keyv_c,nseg_f,mvctr_f,keyg_f,keyv_f,psi_c,psi_f)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n1,n2,n3,nterm,ntp,nseg_c,nseg_f,mvctr_c,mvctr_f
  integer, intent(in) :: nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f
  real(gp), intent(in) :: rx,ry,rz,hx,hy,hz
  integer, dimension(ntp), intent(in) :: lx,ly,lz
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(gp), dimension(ntp), intent(in) :: fac_arr
  real(gp), dimension(nterm), intent(in) :: xp,psiat
  real(wp), dimension(mvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(out) :: psi_f
  !local variables
  character(len=*), parameter :: subname='crtonewave'
  integer, parameter ::nw=16000
  logical :: perx,pery,perz
  integer:: iterm,itp,n_gau,ml1,mu1,ml2,mu2,ml3,mu3,i1,i2,i3,i_all,i_stat,iseg,ii,jj,j0,j1,i0,i
  real(gp) :: gau_a,te
  real(wp), dimension(:,:), allocatable :: work,wprojx,wprojy,wprojz
  real(wp), dimension(:,:,:), allocatable :: psig_c
  real(wp), dimension(:,:,:,:), allocatable :: psig_f

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
  allocate(psig_c(nl1_c:nu1_c,nl2_c:nu2_c,nl3_c:nu3_c+ndebug),stat=i_stat)
  call memocc(i_stat,psig_c,'psig_c',subname)
  allocate(psig_f(7,nl1_f:nu1_f,nl2_f:nu2_f,nl3_f:nu3_f+ndebug),stat=i_stat)
  call memocc(i_stat,psig_f,'psig_f',subname)

  !print *,'limits',nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f

  iterm=1
  itp=1
  gau_a=xp(iterm)
  n_gau=lx(itp)
  CALL GAUSS_TO_DAUB(hx,fac_arr(itp),rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw,perx)
  n_gau=ly(itp)
  CALL GAUSS_TO_DAUB(hy,1.0_gp,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw,pery)
  n_gau=lz(itp)
  CALL GAUSS_TO_DAUB(hz,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw,perz)

  ! First term: coarse projector components
  do i3=nl3_c,nu3_c
     do i2=nl2_c,nu2_c
        do i1=nl1_c,nu1_c
           psig_c(i1,i2,i3)=wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
        enddo
     enddo
  enddo

  ! First term: fine projector components
  do i3=nl3_f,nu3_f
     do i2=nl2_f,nu2_f
        do i1=nl1_f,nu1_f
           psig_f(1,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,1)
           psig_f(2,i1,i2,i3)=wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,1)
           psig_f(3,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,1)
           psig_f(4,i1,i2,i3)=wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,2)
           psig_f(5,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,2)
           psig_f(6,i1,i2,i3)=wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,2)
           psig_f(7,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,2)
        enddo
     enddo
  enddo

  do iterm=2,nterm
     gau_a=xp(iterm)
     n_gau=lx(itp)
     CALL GAUSS_TO_DAUB(hx,fac_arr(itp),rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw,perx)
     n_gau=ly(itp)
     CALL GAUSS_TO_DAUB(hy,1.0_gp,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw,pery)
     n_gau=lz(itp)
     CALL GAUSS_TO_DAUB(hz,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw,perz)

     ! First term: coarse projector components
     do i3=nl3_c,nu3_c
        do i2=nl2_c,nu2_c
           do i1=nl1_c,nu1_c
              psig_c(i1,i2,i3)=psig_c(i1,i2,i3)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
           enddo
        enddo
     enddo

     ! First term: fine projector components
     do i3=nl3_f,nu3_f
        do i2=nl2_f,nu2_f
           do i1=nl1_f,nu1_f
              psig_f(1,i1,i2,i3)=psig_f(1,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,1)
              psig_f(2,i1,i2,i3)=psig_f(2,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,1)
              psig_f(3,i1,i2,i3)=psig_f(3,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,1)
              psig_f(4,i1,i2,i3)=psig_f(4,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,2)
              psig_f(5,i1,i2,i3)=psig_f(5,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,2)
              psig_f(6,i1,i2,i3)=psig_f(6,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,2)
              psig_f(7,i1,i2,i3)=psig_f(7,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,2)
           enddo
        enddo
     enddo

  end do

  do itp=2,ntp

     do iterm=1,nterm
        gau_a=xp(iterm)
        n_gau=lx(itp)
        CALL GAUSS_TO_DAUB(hx,fac_arr(itp),rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw,&
             perx)
        n_gau=ly(itp)
        CALL GAUSS_TO_DAUB(hy,1.0_gp,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw,pery)
        n_gau=lz(itp)
        CALL GAUSS_TO_DAUB(hz,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw,&
             perz)

        ! First term: coarse projector components
        do i3=nl3_c,nu3_c
           do i2=nl2_c,nu2_c
              do i1=nl1_c,nu1_c
                 psig_c(i1,i2,i3)=psig_c(i1,i2,i3)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
              enddo
           enddo
        enddo

        ! First term: fine projector components
        do i3=nl3_f,nu3_f
           do i2=nl2_f,nu2_f
              do i1=nl1_f,nu1_f
                 psig_f(1,i1,i2,i3)=psig_f(1,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,1)
                 psig_f(2,i1,i2,i3)=psig_f(2,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,1)
                 psig_f(3,i1,i2,i3)=psig_f(3,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,1)
                 psig_f(4,i1,i2,i3)=psig_f(4,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,2)
                 psig_f(5,i1,i2,i3)=psig_f(5,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,2)
                 psig_f(6,i1,i2,i3)=psig_f(6,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,2)
                 psig_f(7,i1,i2,i3)=psig_f(7,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,2)
              enddo
           enddo
        enddo

     end do


  end do


  !wavefunction compression

  !itp=0
  ! coarse part
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
        !itp=itp+1
        psi_c(i-i0+jj)=psig_c(i,i2,i3)
     enddo
  enddo

  !print *,'nvctr_c',itp,mvctr_c

  !itp=0
  ! fine part
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
        psi_f(1,i-i0+jj)=psig_f(1,i,i2,i3)
        psi_f(2,i-i0+jj)=psig_f(2,i,i2,i3)
        psi_f(3,i-i0+jj)=psig_f(3,i,i2,i3)
        psi_f(4,i-i0+jj)=psig_f(4,i,i2,i3)
        psi_f(5,i-i0+jj)=psig_f(5,i,i2,i3)
        psi_f(6,i-i0+jj)=psig_f(6,i,i2,i3)
        psi_f(7,i-i0+jj)=psig_f(7,i,i2,i3)
     enddo
  enddo

  !print *,'nvctr_f',itp,mvctr_f

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
  i_all=-product(shape(psig_c))*kind(psig_c)
  deallocate(psig_c,stat=i_stat)
  call memocc(i_stat,i_all,'psig_c',subname)
  i_all=-product(shape(psig_f))*kind(psig_f)
  deallocate(psig_f,stat=i_stat)
  call memocc(i_stat,i_all,'psig_f',subname)

END SUBROUTINE crtonewave
