!> @file
!!  Routines to check the accuracy of the gaussian expansion
!! @author
!!    Copyright (C) 2007-2011 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>  Control the accuracy of the expansion in gaussian
subroutine check_gaussian_expansion(iproc,nproc,orbs,Lzd,psi,G,coeffs)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  type(gaussian_basis), intent(in) :: G
  real(wp), dimension(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs%norbp), intent(in) :: psi
  real(wp), dimension(G%ncoeff,orbs%norbp), intent(in) :: coeffs
  !local variables
  character(len=*), parameter :: subname='check_gaussian_expansion'
  integer :: iorb,i_stat,i_all,i,j,ierr
  real(wp) :: maxdiffp,maxdiff,orbdiff
  real(wp), dimension(:), allocatable :: workpsi

  allocate(workpsi((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,workpsi,'workpsi',subname)

  !call gaussians_to_wavelets(iproc,nproc,lr%geocode,orbs,lr%d,hx,hy,hz,&
  !     lr%wfd,G,coeffs,workpsi)

  call gaussians_to_wavelets_new(iproc,nproc,Lzd,orbs,G,coeffs,workpsi)

  maxdiffp=0.0_wp
  do iorb=1,orbs%norbp
     orbdiff=0.0_wp
     !if (iorb+iproc*norbp <= norb) then
        do i=1,Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f
           j=i+(iorb-1)*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)
           orbdiff=max(orbdiff,(psi(i,iorb)-workpsi(j))**2)
        end do
     !end if
     maxdiffp=max(maxdiffp,orbdiff)
     !print *,'iproc,iorb,orbdiff',iorb,orbdiff
  end do

  if (nproc > 1) then
     call MPI_REDUCE(maxdiffp,maxdiff,1,mpidtypw,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  else
     maxdiff=maxdiffp
  end if

  if (iproc == 0) then
     write(*,'(1x,a,1pe12.4)')'Mean L2 norm of gaussian-wavelet difference:',&
          sqrt(maxdiff/real(orbs%norb,wp))
  end if
  i_all=-product(shape(workpsi))*kind(workpsi)
  deallocate(workpsi,stat=i_stat)
  call memocc(i_stat,i_all,'workpsi',subname)

END SUBROUTINE check_gaussian_expansion


!> Parse the output of CP2K to read the basis set information
subroutine parse_cp2k_files(iproc,basisfile,orbitalfile,nat,ntypes,orbs,iatype,rxyz,&
     CP2K,wfn_cp2k)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: basisfile,orbitalfile
  integer, intent(in) :: iproc,nat,ntypes
  type(orbitals_data), intent(in) :: orbs
  integer, dimension(nat), intent(in) :: iatype
  real(gp), dimension(3,nat), target, intent(in) :: rxyz
  type(gaussian_basis), intent(out) :: CP2K
  real(wp), dimension(:,:), pointer :: wfn_cp2k
  !local variables
  character(len=*), parameter :: subname='parse_cp2k_files'
  character(len=6) :: string,symbol
  character(len=100) :: line
  !n(c) integer, parameter :: nterm_max=3
  integer :: ngx,nbx,nst,nend,num,mmax,myshift,i,ipar,ipg,jat
  integer :: iorb,jorb,iat,ityp,i_all,i_stat,ibas,ig,iset,jbas,ishell,lmax
  integer :: isat,iexpo,icoeff,iam
  real(dp) :: tt
  real(gp) :: exponent,coefficient
  integer, dimension(:), allocatable :: nshell,iorbtmp
  integer, dimension(:,:), allocatable :: nam,ndoc
  real(gp), dimension(:), allocatable :: ctmp
  real(gp), dimension(:,:,:), allocatable :: contcoeff,expo
  real(wp), dimension(:,:,:,:), allocatable :: cimu


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


!!!  !print the found values
!!!  do ityp=1,ntypes
!!!     do ishell=1,nshell(ityp)
!!!        print *,'ityp=',ityp,'ishell=',ishell,'l=',nam(ishell,ityp),'ndoc=',ndoc(ishell,ityp)
!!!        do ipg=1,ndoc(ishell,ityp)
!!!           print *,'expo=',expo(ipg,ishell,ityp),'coeff=',contcoeff(ipg,ishell,ityp)
!!!        end do
!!!     end do
!!!  end do


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
  allocate(cimu(mmax,nbx,nat,orbs%norb+ndebug),stat=i_stat)
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

!!!  nst=1
!!!  nend=4
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
                 !if (iproc==0) 
                    write(*,'(1x,a,i0,a)')&
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
            .and. iat==nat .and. nend==orbs%norb) then
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
              !if (iproc==0)
                    write(*,'(1x,a,i0,a)')&
                   'Problem in the gaucoeff.dat file, only ',iat ,' atoms processed'
              stop
           else
              cycle store_coeff
           end if
        end if
     end do read_orbitals

  end do store_coeff
  close(36)

!!!  !print the found values
!!!  do iat=1,nat
!!!     ityp=iatype(iat)
!!!     do ishell=1,nshell(ityp)
!!!        do jbas=1,2*nam(ishell,ityp)+1
!!!           print *,iat,ishell,nam(ishell,ityp),jbas,(cimu(jbas,ishell,iat,iorb),iorb=1,norb)
!!!        end do
!!!     end do
!!!  end do

  !allocate and assign the coefficients of each orbital
  allocate(wfn_cp2k(CP2K%ncoeff,orbs%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,wfn_cp2k,'wfn_cp2k',subname)
  do iorb=1,orbs%norbp
     jorb=iorb+orbs%isorb
     icoeff=0
     ishell=0
     do iat=1,CP2K%nat
        do isat=1,CP2K%nshell(iat)
           ishell=ishell+1
           do iam=1,2*CP2K%nam(ishell)-1
              icoeff=icoeff+1
              if (jorb <= orbs%norb) wfn_cp2k(icoeff,iorb)=cimu(iam,isat,iat,jorb)
           end do
        end do
     end do
     call gaudim_check(1,icoeff+1,ishell,0,CP2K%ncoeff,CP2K%nshltot)
  end do

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

END SUBROUTINE parse_cp2k_files


subroutine gaussians_to_wavelets(iproc,nproc,geocode,orbs,grid,hx,hy,hz,wfd,G,wfn_gau,psi)
  use module_base
  use module_types
  use yaml_output
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc
  real(gp), intent(in) :: hx,hy,hz
  type(grid_dimensions), intent(in) :: grid
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(orbitals_data), intent(in) :: orbs
  type(gaussian_basis), intent(in) :: G
  real(wp), dimension(G%ncoeff,orbs%nspinor,orbs%norbp), intent(in) :: wfn_gau
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(out) :: psi

  !local variables
  character(len=*), parameter :: subname='gaussians_to_wavelets'
  integer, parameter :: nterm_max=3
  logical :: maycalc
  integer :: i_stat,i_all,ishell,iexpo,icoeff,iat,isat,ng,l,m,iorb,jorb,nterm,ierr,ispinor
  real(dp) :: normdev,tt,scpr,totnorm
  real(gp) :: rx,ry,rz
  integer, dimension(nterm_max) :: lx,ly,lz
  real(gp), dimension(nterm_max) :: fac_arr
  real(wp), dimension(:), allocatable :: tpsi
  
  !if(iproc == 0 .and. verbose > 1) then
     !write(*,'(1x,a)',advance='no')'Writing wavefunctions in wavelet form '
  !end if

  allocate(tpsi(wfd%nvctr_c+7*wfd%nvctr_f+ndebug),stat=i_stat)
  call memocc(i_stat,tpsi,'tpsi',subname)

  !initialize the wavefunction
  call razero((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%norbp*orbs%nspinor,psi)
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
        !print *,iproc,iat,ishell,G%nam(ishell),G%nshell(iat)
        !multiply the values of the gaussian contraction times the orbital coefficient

        do m=1,2*l-1
           call calc_coeff_inguess(l,m,nterm_max,nterm,lx,ly,lz,fac_arr)
           !control whether the basis element may be
           !contribute to some of the orbital of the processor
           maycalc=.false.
           loop_calc: do iorb=1,orbs%norb*orbs%nkpts
              if (orbs%isorb < iorb .and. iorb <= orbs%isorb+orbs%norbp) then
                 jorb=iorb-orbs%isorb
                 do ispinor=1,orbs%nspinor
                    if (wfn_gau(icoeff,ispinor,jorb) /= 0.0_wp) then
                       maycalc=.true.
                       exit loop_calc
                    end if
                 end do
              end if
           end do loop_calc
           if (maycalc) then
              call crtonewave(geocode,grid%n1,grid%n2,grid%n3,ng,nterm,lx,ly,lz,fac_arr,&
                   G%xp(iexpo),G%psiat(iexpo),&
                   rx,ry,rz,hx,hy,hz,&
                   0,grid%n1,0,grid%n2,0,grid%n3,&
                   grid%nfl1,grid%nfu1,grid%nfl2,grid%nfu2,grid%nfl3,grid%nfu3,  & 
                   wfd%nseg_c,wfd%nvctr_c,wfd%keygloc,wfd%keyvloc,wfd%nseg_f,wfd%nvctr_f,&
                   wfd%keygloc(1,wfd%nseg_c+min(1,wfd%nseg_f)),&
                   wfd%keyvloc(wfd%nseg_c+min(1,wfd%nseg_f)),&
                   tpsi(1),tpsi(wfd%nvctr_c+min(1,wfd%nvctr_f)))
           end if
           !sum the result inside the orbital wavefunction
           !loop over the orbitals
           do iorb=1,orbs%norb*orbs%nkpts
              if (orbs%isorb < iorb .and. iorb <= orbs%isorb+orbs%norbp) then
                 jorb=iorb-orbs%isorb
                 do ispinor=1,orbs%nspinor
                    call axpy(wfd%nvctr_c+7*wfd%nvctr_f,wfn_gau(icoeff,ispinor,jorb),&
                         tpsi(1),1,psi(1,ispinor,jorb),1)
                 end do
              end if
           end do
           icoeff=icoeff+1
        end do
        iexpo=iexpo+ng
     end do
!     if (iproc == 0 .and. verbose > 1) then
!        write(*,'(a)',advance='no') &
!             repeat('.',(iat*40)/G%nat-((iat-1)*40)/G%nat)
!     end if
  end do

  call gaudim_check(iexpo,icoeff,ishell,G%nexpo,G%ncoeff,G%nshltot)

  if (iproc ==0  .and. verbose > 1) then
     call yaml_map('Wavelet conversion succeeded',.true.)
     !write(*,'(1x,a)')'done.'
  end if
  !renormalize the orbitals
  !calculate the deviation from 1 of the orbital norm
  normdev=0.0_dp
  tt=0.0_dp
  do iorb=1,orbs%norb*orbs%nkpts
     if (orbs%isorb < iorb .and. iorb <= orbs%isorb+orbs%norbp) then
        jorb=iorb-orbs%isorb
        totnorm=0.0_dp
       do ispinor=1,orbs%nspinor !to be verified in case of nspinor=4
          call wnrm_wrap(1,wfd%nvctr_c,wfd%nvctr_f,psi(1,ispinor,jorb),scpr) 
           totnorm=totnorm+scpr

           !print *,'AAA',iproc,iorb,ispinor,scpr,jorb
        end do
        do ispinor=1,orbs%nspinor !to be verified in case of nspinor=4
           call wscal_wrap(wfd%nvctr_c,wfd%nvctr_f,real(1.0_dp/sqrt(totnorm),wp),&
                psi(1,ispinor,jorb))

           !call wnrm_wrap(wfd%nvctr_c,wfd%nvctr_f,psi(1,ispinor,jorb),scpr) 
           !print *,'BBB',iproc,iorb,ispinor,scpr

        end do
        !write(*,'(1x,a,i5,1pe14.7)')'norm of orbital ',iorb,totnorm
        tt=max(tt,abs(1.0_dp-totnorm))
     end if
  end do
  if (nproc > 1) then
     call MPI_REDUCE(tt,normdev,1,mpidtypd,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  else
     normdev=tt
  end if
  if (iproc ==0) then
     !write(*,'(1x,a,1pe12.2)')&
     !  'Deviation from normalization of the imported orbitals',normdev
     call yaml_map('Deviation from normalization',normdev,fmt='(1pe12.2)')
  end if

  i_all=-product(shape(tpsi))*kind(tpsi)
  deallocate(tpsi,stat=i_stat)
  call memocc(i_stat,i_all,'tpsi',subname)

END SUBROUTINE gaussians_to_wavelets

subroutine gaussians_to_wavelets_new(iproc,nproc,Lzd,orbs,G,wfn_gau,psi)
  use module_base
  use module_types
  use yaml_output
  implicit none
  integer, intent(in) :: iproc,nproc
  type(local_zone_descriptors), intent(in) :: Lzd
  type(orbitals_data), intent(in) :: orbs
  type(gaussian_basis), intent(in) :: G
  real(wp), dimension(G%ncoeff,orbs%nspinor,orbs%norbp), intent(in) :: wfn_gau
  real(wp), dimension(orbs%npsidim_orbs), intent(out) :: psi

  !local variables
  integer :: iorb,ierr,ispinor,ncplx,ind,ind2,ilr
  real(dp) :: normdev,tt,scpr,totnorm
  real(gp) :: kx,ky,kz

  !if(iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no')&
  !     'Writing wavefunctions in wavelet form...'
  
  normdev=0.0_dp
  tt=0.0_dp  
  ind = 1
  ind2 = 1
  do iorb=1,orbs%norbp
     ilr = orbs%inWhichLocreg(iorb+orbs%isorb)
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
        !if (iproc == 0)print *,'start',ispinor,ncplx,iorb+orbs%isorb,orbs%nspinor
        !the Block wavefunctions are exp(-Ikr) psi(r) (with MINUS k)
        call gaussians_to_wavelets_orb(ncplx,Lzd%Llr(ilr),&
             Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),kx,ky,kz,G,&
             wfn_gau(1,ispinor,iorb),psi(ind))

        !if (iproc == 0)print *,'end',ispinor,ncplx,iorb+orbs%isorb,orbs%nspinor
        call wnrm_wrap(ncplx,Lzd%Llr(ilr)%wfd%nvctr_c,Lzd%Llr(ilr)%wfd%nvctr_f,psi(ind),scpr) 
        totnorm=totnorm+scpr
        ind = ind + (Lzd%Llr(ilr)%wfd%nvctr_c + 7*Lzd%Llr(ilr)%wfd%nvctr_f)*ncplx
     end do
     !write(*,'(1x,a,i5,1pe14.7,i3)')'norm of orbital ',iorb,totnorm,ncplx
     do ispinor=1,orbs%nspinor
        call wscal_wrap(Lzd%Llr(ilr)%wfd%nvctr_c,Lzd%Llr(ilr)%wfd%nvctr_f,real(1.0_dp/sqrt(totnorm),wp),&
             psi(ind2))
        ind2 = ind2 + (Lzd%Llr(ilr)%wfd%nvctr_c + 7*Lzd%Llr(ilr)%wfd%nvctr_f)
     end do
     tt=max(tt,abs(1.0_dp-totnorm))
     !print *,'iorb,norm',totnorm
  end do

!  if (iproc ==0  .and. verbose > 1) write(*,'(1x,a)')'done.'
  if (iproc ==0  .and. verbose > 1) then
     call yaml_map('Wavelet conversion succeeded',.true.)
     !write(*,'(1x,a)')'done.'
  end if

  !renormalize the orbitals
  !calculate the deviation from 1 of the orbital norm
  if (nproc > 1) then
     call MPI_REDUCE(tt,normdev,1,mpidtypd,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  else
     normdev=tt
  end if
  if (iproc ==0) then
     !write(*,'(1x,a,1pe12.2)')&
     !  'Deviation from normalization of the imported orbitals',normdev
     call yaml_map('Deviation from normalization',normdev,fmt='(1pe12.2)')
  end if

END SUBROUTINE gaussians_to_wavelets_new



subroutine gaussians_to_wavelets_orb(ncplx,lr,hx,hy,hz,kx,ky,kz,G,wfn_gau,psi)
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
  integer, parameter :: nterm_max=3,maxsizeKB=2048,nw=65536
  logical :: perx,pery,perz
  integer :: i_stat,i_all,ishell,iexpo,icoeff,iat,isat,ng,l,m,i,nterm,ig
  integer :: nterms_max,nterms,iterm,n_gau,ml1,mu1,ml2,mu2,ml3,mu3 !n(c) iscoeff
  real(gp) :: rx,ry,rz,gau_a
  integer, dimension(nterm_max) :: lx,ly,lz
  real(gp), dimension(nterm_max) :: fac_arr
  real(wp), allocatable, dimension(:,:,:) :: work
  real(wp), allocatable, dimension(:,:,:,:) :: wx,wy,wz

  !calculate nterms_max:
  !allows only maxsizeKB per one-dimensional array
  !(for a grid of dimension 100 nterms_max=655)
  !but with at least ngx*nterm_max ~= 100 elements
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
  call to_zero((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*ncplx,psi(1))

  !calculate the number of terms for this orbital
  nterms=0
  !loop over the atoms
  ishell=0
  iexpo=1
  icoeff=1
  !n(c) iscoeff=1
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
        !print *,iproc,iat,ishell,G%nam(ishell),G%nshell(iat)
        !multiply the values of the gaussian contraction times the orbital coefficient

        do m=1,2*l-1
           call calc_coeff_inguess(l,m,nterm_max,nterm,lx,ly,lz,fac_arr)
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
                         lr%ns1,lr%d%n1,ml1,mu1,&
                         wx(1,0,1,iterm),work,nw,perx) 
                    n_gau=ly(i)
                    !print *,'y',ml2,mu2,lr%d%n2
                    call gauss_to_daub_k(hy,ky*hy,ncplx,wfn_gau(icoeff),ry,gau_a,n_gau,&
                         lr%ns2,lr%d%n2,ml2,mu2,&
                         wy(1,0,1,iterm),work,nw,pery) 
                    n_gau=lz(i) 
                    !print *,'z',ml3,mu3,lr%d%n3
                    call gauss_to_daub_k(hz,kz*hz,ncplx,G%psiat(iexpo+ig-1),rz,gau_a,n_gau,&
                         lr%ns3,lr%d%n3,ml3,mu3,&
                         wz(1,0,1,iterm),work,nw,perz)
                    iterm=iterm+1
                 end do
              end do
              nterms=nterms+nterm*ng
              !print *,'nterms',nterms,nterms_max
           end if
           icoeff=icoeff+1
        end do
        iexpo=iexpo+ng
     end do
  end do

  call gaudim_check(iexpo,icoeff,ishell,G%nexpo,G%ncoeff,G%nshltot)

  !accumulate wavefunction
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

END SUBROUTINE gaussians_to_wavelets_orb




subroutine gaussians_c_to_wavelets_orb(ncplx,lr,hx,hy,hz,kx,ky,kz,G,wfn_gau,psi, cutoff)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ncplx
  real(gp), intent(in) :: hx,hy,hz,kx,ky,kz
  type(locreg_descriptors), intent(in) :: lr
  type(gaussian_basis_c), intent(in) :: G
  real(wp), dimension(G%ncoeff), intent(in) :: wfn_gau
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*ncplx), intent(out) :: psi
  real(gp) cutoff

  !local variables
  character(len=*), parameter :: subname='gaussians_to_wavelets_orb'
  integer, parameter :: nterm_max=3,maxsizeKB=2048,nw=65536
  logical :: perx,pery,perz
  integer :: i_stat,i_all,ishell,iexpo,icoeff,iat,isat,ng,l,m,i,nterm,ig
  integer :: nterms_max,nterms,iscoeff,iterm,n_gau,ml1,mu1,ml2,mu2,ml3,mu3
  real(gp) :: rx,ry,rz,gau_a, gau_bf
  integer, dimension(nterm_max) :: lx,ly,lz
  real(gp), dimension(nterm_max) :: fac_arr
  real(wp), allocatable, dimension(:,:,:, :) :: work
  real(wp), allocatable, dimension(:,  :,:,:,:) :: wx,wy,wz
  real(wp), allocatable, dimension(:,:) :: cossinfacts
  integer :: ncplxC

  ncplxC=2

  ! if ( ncplx.ne.1) then
  !    stop ' ncplx must be 1 in  actual version of  gaussians_c_to_wavelets_orb'
  ! end if


  !calculate nterms_max:
  !allows only maxsizeKB per one-dimensional array
  !(for a grid of dimension 100 nterms_max=655)
  !but with at least ngx*nterm_max ~= 100 elements
  nterms_max=max(maxsizeKB*1024/(2*ncplxC*max(lr%d%n1,lr%d%n2,lr%d%n3)),100)

  allocate(work(0:nw,2,2, ncplx+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)

  allocate(wx( ncplxC, ncplx,0:lr%d%n1,2,nterms_max+ndebug),stat=i_stat)
  call memocc(i_stat,wx,'wx',subname)
  allocate(wy( ncplxC, ncplx,0:lr%d%n2,2,nterms_max+ndebug),stat=i_stat)
  call memocc(i_stat,wy,'wy',subname)
  allocate(wz(ncplxC, ncplx,0:lr%d%n3,2,nterms_max+ndebug),stat=i_stat)
  call memocc(i_stat,wz,'wz',subname)

  allocate(   cossinfacts(1:2, 1:nterms_max+ndebug) ,stat=i_stat)
  call memocc(i_stat,cossinfacts,'cossinfacts',subname)


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
        l=abs(G%nam(ishell)  )
        !print *,iproc,iat,ishell,G%nam(ishell),G%nshell(iat)
        !multiply the values of the gaussian contraction times the orbital coefficient
    

        do m=1,2*l-1
           call calc_coeff_inguess(l,m,nterm_max,nterm,lx,ly,lz,fac_arr)
           !control whether the basis element may be
           !contribute to some of the orbital of the processor

           if (wfn_gau(icoeff) /= 0.0_wp) then
              if (nterms + nterm*ng > nterms_max) then
                 call wfn_from_tensprod_cossin(lr, ncplx,cossinfacts , nterms,wx,wy,wz,psi)
                 iterm=1
                 nterms=0
              end if
              !assign the arrays
              !make sure that the coefficients returned by 
              !gauss_to_daub are zero outside [ml:mr] 
              do ig=1,ng

                 gau_a= sqrt( 1.0_gp/REAL(2.0_gp*G%expof(iexpo+ig-1))   )
                 gau_bf = AIMAG ( G%expof(iexpo+ig-1) )

                 do i=1,nterm
                    
                    n_gau=lx(i)

                    call gauss_c_to_daub_k(hx,kx,ncplx,gau_bf ,ncplxC,fac_arr(i), &
                         rx,gau_a,  n_gau,&
                         lr%d%n1,ml1,mu1,&
                         wx(1,1,0,1,iterm),work,nw,perx, cutoff) 

                    n_gau=ly(i)
                    !print *,'y',ml2,mu2,lr%d%n2
                    call gauss_c_to_daub_k(hy,ky,ncplx,gau_bf,ncplxC,wfn_gau(icoeff), &
                         ry,gau_a,n_gau,&
                         lr%d%n2,ml2,mu2,&
                         wy(1,1,0,1,iterm),work,nw,pery, cutoff) 
                    n_gau=lz(i) 
                    !print *,'z',ml3,mu3,lr%d%n3
                    call gauss_c_to_daub_k(hz,kz,ncplx,gau_bf,ncplxC,  1.0_wp,  &
                         rz,gau_a,n_gau,&
                         lr%d%n3,ml3,mu3,&
                         wz(1,1,0,1,iterm),work,nw,perz, cutoff)

                    cossinfacts(1,iterm)= REAL( G%psiat(iexpo+ig-1))
                    cossinfacts(2,iterm)= AIMAG(G%psiat(iexpo+ig-1)) 

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


  call wfn_from_tensprod_cossin(lr, ncplx,  cossinfacts    ,nterms,wx,wy,wz,psi)





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

  i_all=-product(shape(cossinfacts))*kind(cossinfacts)
  deallocate(cossinfacts,stat=i_stat)
  call memocc(i_stat,i_all,'cossinfacts',subname)


  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)



END SUBROUTINE gaussians_c_to_wavelets_orb



!> Accumulate 3d wavefunction in complex form from a tensor produc decomposition
!! universal routine which should be used for all gautowav operations
subroutine wfn_from_tensprod(lr,ncplx,nterm,wx,wy,wz,psi)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ncplx,nterm
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(ncplx,0:lr%d%n1,2,nterm), intent(in) :: wx
  real(wp), dimension(ncplx,0:lr%d%n2,2,nterm), intent(in) :: wy
  real(wp), dimension(ncplx,0:lr%d%n3,2,nterm), intent(in) :: wz
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*ncplx), intent(inout) :: psi
  !local variables
  integer :: iseg,i,i0,i1,i2,i3,jj,ind_c,ind_f,iterm,nvctr
  real(wp) :: re_cmplx_prod,im_cmplx_prod
  !!$  integer :: ithread,nthread,omp_get_thread_num,omp_get_num_threads

  !the filling of the wavefunction should be different if ncplx==1 or 2
  !split such as to avoid intensive call to if statements

  !!$omp parallel default(private) shared(lr%nseg_c,lr%wfd%keyv,lr%wfd%keyg,lr%d) &
  !!$omp shared(psi,wx,wy,wz,lr%wfd%nvctr_c) &
  !!$omp shared(nterm,lr%wfd%nvctr_f,lr%wfd%nseg_f)

  !!$	ithread=omp_get_thread_num()
  !!$	nthread=omp_get_num_threads()
  if (ncplx == 1) then

     !!$  if(ithread .eq. 0) then
     ! Other terms: coarse projector components
     ! coarse part
     nvctr=0
     do iseg=1,lr%wfd%nseg_c
        call segments_to_grid(lr%wfd%keyvloc(iseg),lr%wfd%keygloc(1,iseg),lr%d,i0,i1,i2,i3,jj)
        do i=i0,i1
           ind_c=i-i0+jj
           do iterm=1,nterm
              psi(ind_c)=psi(ind_c)+&
                   wx(1,i,1,iterm)*wy(1,i2,1,iterm)*wz(1,i3,1,iterm)
           end do
           nvctr=nvctr+1
        end do
     end do

     if (nvctr /=  lr%wfd%nvctr_c) then
        write(*,'(1x,a,i0,1x,i0)')' ERROR: nvctr /= nvctr_c ',nvctr,lr%wfd%nvctr_c
        stop
     end if
     !!$  end if

     !!$  if(ithread .eq. 1 .or. nthread .eq. 1) then
     ! Other terms: fine projector components
     nvctr=0
     do iseg=lr%wfd%nseg_c+1,lr%wfd%nseg_c+lr%wfd%nseg_f
        call segments_to_grid(lr%wfd%keyvloc(iseg),lr%wfd%keygloc(1,iseg),lr%d,i0,i1,i2,i3,jj)
        do i=i0,i1
           ind_f=lr%wfd%nvctr_c+7*(i-i0+jj-1)
           do iterm=1,nterm
              psi(ind_f+1)=psi(ind_f+1)+&
                   wx(1,i,2,iterm)*wy(1,i2,1,iterm)*wz(1,i3,1,iterm)
              psi(ind_f+2)=psi(ind_f+2)+&
                   wx(1,i,1,iterm)*wy(1,i2,2,iterm)*wz(1,i3,1,iterm)
              psi(ind_f+3)=psi(ind_f+3)+&
                   wx(1,i,2,iterm)*wy(1,i2,2,iterm)*wz(1,i3,1,iterm)
              psi(ind_f+4)=psi(ind_f+4)+&
                   wx(1,i,1,iterm)*wy(1,i2,1,iterm)*wz(1,i3,2,iterm)
              psi(ind_f+5)=psi(ind_f+5)+&
                   wx(1,i,2,iterm)*wy(1,i2,1,iterm)*wz(1,i3,2,iterm)
              psi(ind_f+6)=psi(ind_f+6)+&
                   wx(1,i,1,iterm)*wy(1,i2,2,iterm)*wz(1,i3,2,iterm)
              psi(ind_f+7)=psi(ind_f+7)+&
                   wx(1,i,2,iterm)*wy(1,i2,2,iterm)*wz(1,i3,2,iterm)
           end do
           nvctr=nvctr+1
        end do
     end do
     if (nvctr /= lr%wfd%nvctr_f) then
        write(*,'(1x,a,i0,1x,i0)')' ERROR: nvctr /= nvctr_f ',nvctr,lr%wfd%nvctr_f
        stop 
     end if
     !!$  end if
  else if (ncplx ==2) then

     !part with real and imaginary part
     !modify the openMP statements such as to benefit from parallelisation

     !!$  if(ithread .eq. 0) then
     ! Other terms: coarse projector components
     ! coarse part
     nvctr=0
     do iseg=1,lr%wfd%nseg_c
        call segments_to_grid(lr%wfd%keyvloc(iseg),lr%wfd%keygloc(1,iseg),lr%d,i0,i1,i2,i3,jj)
        do i=i0,i1
           ind_c=i-i0+jj
           do iterm=1,nterm
              psi(ind_c)=psi(ind_c)+re_cmplx_prod(&
                   wx(1,i,1,iterm),wy(1,i2,1,iterm),wz(1,i3,1,iterm))
           end do
           nvctr=nvctr+1
        end do
     end do
     if (nvctr /=  lr%wfd%nvctr_c) then
        write(*,'(1x,a,i0,1x,i0)')' ERROR: nvctr /= nvctr_c ',nvctr,lr%wfd%nvctr_c
        stop
     end if
     !!$  end if

     !!$  if(ithread .eq. 1 .or. nthread .eq. 1) then
     ! Other terms: fine projector components
     nvctr=0
     do iseg=lr%wfd%nseg_c+1,lr%wfd%nseg_c+lr%wfd%nseg_f
        call segments_to_grid(lr%wfd%keyvloc(iseg),lr%wfd%keygloc(1,iseg),lr%d,i0,i1,i2,i3,jj)
        do i=i0,i1
           ind_f=lr%wfd%nvctr_c+7*(i-i0+jj-1)
           do iterm=1,nterm
              psi(ind_f+1)=psi(ind_f+1)+re_cmplx_prod(&
                   wx(1,i,2,iterm),wy(1,i2,1,iterm),wz(1,i3,1,iterm))
              psi(ind_f+2)=psi(ind_f+2)+re_cmplx_prod(&
                   wx(1,i,1,iterm),wy(1,i2,2,iterm),wz(1,i3,1,iterm))
              psi(ind_f+3)=psi(ind_f+3)+re_cmplx_prod(&
                   wx(1,i,2,iterm),wy(1,i2,2,iterm),wz(1,i3,1,iterm))
              psi(ind_f+4)=psi(ind_f+4)+re_cmplx_prod(&
                   wx(1,i,1,iterm),wy(1,i2,1,iterm),wz(1,i3,2,iterm))
              psi(ind_f+5)=psi(ind_f+5)+re_cmplx_prod(&
                   wx(1,i,2,iterm),wy(1,i2,1,iterm),wz(1,i3,2,iterm))
              psi(ind_f+6)=psi(ind_f+6)+re_cmplx_prod(&
                   wx(1,i,1,iterm),wy(1,i2,2,iterm),wz(1,i3,2,iterm))
              psi(ind_f+7)=psi(ind_f+7)+re_cmplx_prod(&
                   wx(1,i,2,iterm),wy(1,i2,2,iterm),wz(1,i3,2,iterm))
           end do
           nvctr=nvctr+1
        end do
     end do
     if (nvctr /= lr%wfd%nvctr_f) then
        write(*,'(1x,a,i0,1x,i0)')' ERROR: nvctr /= nvctr_f ',nvctr,lr%wfd%nvctr_f
        stop 
     end if
     !!$  end if
     
     !now the imaginary part
     
     !!$  if((ithread == 0 .and. nthread <= 2) .or. ithread == 2) then 
     ! Other terms: coarse projector components
     ! coarse part
     do iseg=1,lr%wfd%nseg_c
        call segments_to_grid(lr%wfd%keyvloc(iseg),lr%wfd%keygloc(1,iseg),lr%d,i0,i1,i2,i3,jj)
        do i=i0,i1
           ind_c=lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+i-i0+jj
           do iterm=1,nterm
              psi(ind_c)=psi(ind_c)+im_cmplx_prod(&
                   wx(1,i,1,iterm),wy(1,i2,1,iterm),wz(1,i3,1,iterm))
           end do
        end do
     end do

     !!$  end if

     !!$  if((ithread .eq. 1 .and. nthread <=3) .or. nthread .eq. 1 .or. ithread == 3) then
     ! Other terms: fine projector components
     do iseg=lr%wfd%nseg_c+1,lr%wfd%nseg_c+lr%wfd%nseg_f
        call segments_to_grid(lr%wfd%keyvloc(iseg),lr%wfd%keygloc(1,iseg),lr%d,i0,i1,i2,i3,jj)
        do i=i0,i1
           ind_f=lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+lr%wfd%nvctr_c+7*(i-i0+jj-1)
           do iterm=1,nterm
              psi(ind_f+1)=psi(ind_f+1)+im_cmplx_prod(&
                   wx(1,i,2,iterm),wy(1,i2,1,iterm),wz(1,i3,1,iterm))
              psi(ind_f+2)=psi(ind_f+2)+im_cmplx_prod(&
                   wx(1,i,1,iterm),wy(1,i2,2,iterm),wz(1,i3,1,iterm))
              psi(ind_f+3)=psi(ind_f+3)+im_cmplx_prod(&
                   wx(1,i,2,iterm),wy(1,i2,2,iterm),wz(1,i3,1,iterm))
              psi(ind_f+4)=psi(ind_f+4)+im_cmplx_prod(&
                   wx(1,i,1,iterm),wy(1,i2,1,iterm),wz(1,i3,2,iterm))
              psi(ind_f+5)=psi(ind_f+5)+im_cmplx_prod(&
                   wx(1,i,2,iterm),wy(1,i2,1,iterm),wz(1,i3,2,iterm))
              psi(ind_f+6)=psi(ind_f+6)+im_cmplx_prod(&
                   wx(1,i,1,iterm),wy(1,i2,2,iterm),wz(1,i3,2,iterm))
              psi(ind_f+7)=psi(ind_f+7)+im_cmplx_prod(&
                   wx(1,i,2,iterm),wy(1,i2,2,iterm),wz(1,i3,2,iterm))
           end do
        end do
     end do
     !!$  end if
  end if

  !!$omp end parallel

END SUBROUTINE wfn_from_tensprod


function re_re_cmplx_prod(a,b,c)
  use module_base
  implicit none
  real(wp), dimension(2,2), intent(in) :: a,b,c
  real(wp) :: re_re_cmplx_prod
  real(wp) :: re_cmplx_prod
  
  re_re_cmplx_prod=re_cmplx_prod( a(1,1),b(1,1),c(1,1)) &
       -re_cmplx_prod( a(1,1),b(1,2),c(1,2)) &
       -re_cmplx_prod( a(1,2),b(1,1),c(1,2)) &
       -re_cmplx_prod( a(1,2),b(1,2),c(1,1))
END FUNCTION re_re_cmplx_prod


function im_re_cmplx_prod(a,b,c)
  use module_base
  implicit none
  real(wp), dimension(2,2), intent(in) :: a,b,c
  real(wp) :: im_re_cmplx_prod
  real(wp) :: re_cmplx_prod
  
  im_re_cmplx_prod=-re_cmplx_prod(a(1,2),b(1,2),c(1,2)) &
                   +re_cmplx_prod(a(1,2),b(1,1),c(1,1)) &
                   +re_cmplx_prod(a(1,1),b(1,2),c(1,1)) &
                   +re_cmplx_prod(a(1,1),b(1,1),c(1,2))
  
END FUNCTION im_re_cmplx_prod


function re_im_cmplx_prod(a,b,c)
  use module_base
  implicit none
  real(wp), dimension(2,2), intent(in) :: a,b,c
  real(wp) :: re_im_cmplx_prod
  real(wp) :: im_cmplx_prod
  
  re_im_cmplx_prod=im_cmplx_prod( a(1,1),b(1,1),c(1,1)) &
       -im_cmplx_prod( a(1,1),b(1,2),c(1,2)) &
       -im_cmplx_prod( a(1,2),b(1,1),c(1,2)) &
       -im_cmplx_prod( a(1,2),b(1,2),c(1,1))
  
END FUNCTION re_im_cmplx_prod


function im_im_cmplx_prod(a,b,c)
  use module_base
  implicit none
  real(wp), dimension(2,2), intent(in) :: a,b,c
  real(wp) :: im_im_cmplx_prod
  real(wp) :: im_cmplx_prod
  
  im_im_cmplx_prod=-im_cmplx_prod(a(1,2),b(1,2),c(1,2)) &
                   +im_cmplx_prod(a(1,2),b(1,1),c(1,1)) &
                   +im_cmplx_prod(a(1,1),b(1,2),c(1,1)) &
                   +im_cmplx_prod(a(1,1),b(1,1),c(1,2))  
END FUNCTION im_im_cmplx_prod


!> Accumulate 3d projector in real form from a tensor produc decomposition
!! using complex gaussians
subroutine wfn_from_tensprod_cossin(lr,ncplx,  cossinfacts ,nterm,wx,wy,wz,psi)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nterm, ncplx
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(2,ncplx,0:lr%d%n1,2,nterm), intent(in) :: wx
  real(wp), dimension(2,ncplx,0:lr%d%n2,2,nterm), intent(in) :: wy
  real(wp), dimension(2,ncplx,0:lr%d%n3,2,nterm), intent(in) :: wz
  real(wp) :: cossinfacts(2,nterm)


  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*ncplx), intent(inout) :: psi
  !local variables
  integer :: iseg,i,i0,i1,i2,i3,jj,ind_c,ind_f,iterm,nvctr
  real(wp) :: re_cmplx_prod,im_cmplx_prod
  real(wp) :: re_re_cmplx_prod,re_im_cmplx_prod,im_re_cmplx_prod,im_im_cmplx_prod 

  !!$omp parallel default(private) shared(lr%nseg_c,lr%wfd%keyv,lr%wfd%keyg,lr%d) &
  !!$omp shared(psi,wx,wy,wz,lr%wfd%nvctr_c) &
  !!$omp shared(nterm,lr%wfd%nvctr_f,lr%wfd%nseg_f)
  if (ncplx == 1) then
     nvctr=0
     do iseg=1,lr%wfd%nseg_c
        call segments_to_grid(lr%wfd%keyvglob(iseg),lr%wfd%keyglob(1,iseg),lr%d,i0,i1,i2,i3,jj)
        do i=i0,i1
           ind_c=i-i0+jj
           do iterm=1,nterm
              psi(ind_c)=psi(ind_c)+re_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,1,iterm))*cossinfacts(1,iterm)
           end do
           nvctr=nvctr+1
        end do
     end do
     if (nvctr /=  lr%wfd%nvctr_c) then
        write(*,'(1x,a,i0,1x,i0)')' ERROR: nvctr >< nvctr_c ',nvctr,lr%wfd%nvctr_c
        stop
     end if
     !!$  end if

     !!$  if(ithread .eq. 1 .or. nthread .eq. 1) then
     ! Other terms: fine projector components
     nvctr=0
     do iseg=lr%wfd%nseg_c+1,lr%wfd%nseg_c+lr%wfd%nseg_f
        call segments_to_grid(lr%wfd%keyvglob(iseg),lr%wfd%keyglob(1,iseg),lr%d,i0,i1,i2,i3,jj)
        do i=i0,i1
           ind_f=lr%wfd%nvctr_c+7*(i-i0+jj-1)
           do iterm=1,nterm
              psi(ind_f+1)=psi(ind_f+1)+re_cmplx_prod(&
                   wx(1, 1,i,2,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,1,iterm))*cossinfacts(1,iterm)
              psi(ind_f+2)=psi(ind_f+2)+re_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,1,iterm))*cossinfacts(1,iterm)
              psi(ind_f+3)=psi(ind_f+3)+re_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,1,iterm))*cossinfacts(1,iterm)
              psi(ind_f+4)=psi(ind_f+4)+re_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,2,iterm))*cossinfacts(1,iterm)
              psi(ind_f+5)=psi(ind_f+5)+re_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,2,iterm))*cossinfacts(1,iterm)
              psi(ind_f+6)=psi(ind_f+6)+re_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,2,iterm))*cossinfacts(1,iterm)
              psi(ind_f+7)=psi(ind_f+7)+re_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,2,iterm))*cossinfacts(1,iterm)
           end do
           nvctr=nvctr+1
        end do
     end do
     if (nvctr /= lr%wfd%nvctr_f) then
        write(*,'(1x,a,i0,1x,i0)')' ERROR: nvctr >< nvctr_f ',nvctr,lr%wfd%nvctr_f
        stop 
     end if
     !!$  end if
     
     !now the imaginary part
     
     !!$  if((ithread == 0 .and. nthread <= 2) .or. ithread == 2) then 
     ! Other terms: coarse projector components
     ! coarse part
     do iseg=1,lr%wfd%nseg_c
        call segments_to_grid(lr%wfd%keyvglob(iseg),lr%wfd%keyglob(1,iseg),lr%d,i0,i1,i2,i3,jj)
        do i=i0,i1
           ind_c=i-i0+jj
           do iterm=1,nterm
              psi(ind_c)=psi(ind_c)+im_cmplx_prod(&
                   wx(1,1, i,1,iterm),wy(1,1, i2,1,iterm),wz(1,1, i3,1,iterm))*cossinfacts(2,iterm)
           end do
        end do
     end do

     !!$  end if

     !!$  if((ithread .eq. 1 .and. nthread <=3) .or. nthread .eq. 1 .or. ithread == 3) then
     ! Other terms: fine projector components
     do iseg=lr%wfd%nseg_c+1,lr%wfd%nseg_c+lr%wfd%nseg_f
        call segments_to_grid(lr%wfd%keyvglob(iseg),lr%wfd%keyglob(1,iseg),lr%d,i0,i1,i2,i3,jj)
        do i=i0,i1
           ind_f=lr%wfd%nvctr_c+7*(i-i0+jj-1)
           do iterm=1,nterm
              psi(ind_f+1)=psi(ind_f+1)+im_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,1,iterm))*cossinfacts(2,iterm)
              psi(ind_f+2)=psi(ind_f+2)+im_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,1,iterm))*cossinfacts(2,iterm)
              psi(ind_f+3)=psi(ind_f+3)+im_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,1,iterm))*cossinfacts(2,iterm)
              psi(ind_f+4)=psi(ind_f+4)+im_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,2,iterm))*cossinfacts(2,iterm)
              psi(ind_f+5)=psi(ind_f+5)+im_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,2,iterm))*cossinfacts(2,iterm)
              psi(ind_f+6)=psi(ind_f+6)+im_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,2,iterm))*cossinfacts(2,iterm)
              psi(ind_f+7)=psi(ind_f+7)+im_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,2,iterm))*cossinfacts(2,iterm)
           end do
        end do
     end do
  else if (ncplx ==2) then




     nvctr=0
     do iseg=1,lr%wfd%nseg_c
        call segments_to_grid(lr%wfd%keyvglob(iseg),lr%wfd%keyglob(1,iseg),lr%d,i0,i1,i2,i3,jj)

        do i=i0,i1
           ind_c=i-i0+jj
           do iterm=1,nterm
              psi(ind_c)=psi(ind_c)+re_re_cmplx_prod(&
                   wx(1,1, i,1,iterm),wy(1,1, i2,1,iterm),wz(1,1, i3,1,iterm))*cossinfacts(1,iterm)
           end do
           nvctr=nvctr+1
        end do

        do i=i0,i1
           ind_c= lr%wfd%nvctr_c + 7*lr%wfd%nvctr_f + i-i0+jj 
           do iterm=1,nterm
              psi(ind_c)=psi(ind_c)+im_re_cmplx_prod(&
                   wx(1,1, i,1,iterm),wy(1,1, i2,1,iterm),wz(1,1, i3,1,iterm))*cossinfacts(1,iterm)
           end do
        end do



     end do
     if (nvctr /=  lr%wfd%nvctr_c) then
        write(*,'(1x,a,i0,1x,i0)')' ERROR: nvctr >< nvctr_c ',nvctr,lr%wfd%nvctr_c
        stop
     end if
     !!$  end if

     !!$  if(ithread .eq. 1 .or. nthread .eq. 1) then
     ! Other terms: fine projector components
     nvctr=0
     do iseg=lr%wfd%nseg_c+1,lr%wfd%nseg_c+lr%wfd%nseg_f
        call segments_to_grid(lr%wfd%keyvglob(iseg),lr%wfd%keyglob(1,iseg),lr%d,i0,i1,i2,i3,jj)
        do i=i0,i1
           ind_f=lr%wfd%nvctr_c+7*(i-i0+jj-1)
           do iterm=1,nterm
              psi(ind_f+1)=psi(ind_f+1)+re_re_cmplx_prod(&
                   wx(1, 1,i,2,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,1,iterm))*cossinfacts(1,iterm)
              psi(ind_f+2)=psi(ind_f+2)+re_re_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,1,iterm))*cossinfacts(1,iterm)
              psi(ind_f+3)=psi(ind_f+3)+re_re_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,1,iterm))*cossinfacts(1,iterm)
              psi(ind_f+4)=psi(ind_f+4)+re_re_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,2,iterm))*cossinfacts(1,iterm)
              psi(ind_f+5)=psi(ind_f+5)+re_re_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,2,iterm))*cossinfacts(1,iterm)
              psi(ind_f+6)=psi(ind_f+6)+re_re_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,2,iterm))*cossinfacts(1,iterm)
              psi(ind_f+7)=psi(ind_f+7)+re_re_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,2,iterm))*cossinfacts(1,iterm)
           end do
           nvctr=nvctr+1
        end do

        do i=i0,i1
           ind_f=lr%wfd%nvctr_c + 7*lr%wfd%nvctr_f +  lr%wfd%nvctr_c+7*(i-i0+jj-1)
           do iterm=1,nterm
              psi(ind_f+1)=psi(ind_f+1)+im_re_cmplx_prod(&
                   wx(1, 1,i,2,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,1,iterm))*cossinfacts(1,iterm)
              psi(ind_f+2)=psi(ind_f+2)+im_re_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,1,iterm))*cossinfacts(1,iterm)
              psi(ind_f+3)=psi(ind_f+3)+im_re_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,1,iterm))*cossinfacts(1,iterm)
              psi(ind_f+4)=psi(ind_f+4)+im_re_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,2,iterm))*cossinfacts(1,iterm)
              psi(ind_f+5)=psi(ind_f+5)+im_re_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,2,iterm))*cossinfacts(1,iterm)
              psi(ind_f+6)=psi(ind_f+6)+im_re_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,2,iterm))*cossinfacts(1,iterm)
              psi(ind_f+7)=psi(ind_f+7)+im_re_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,2,iterm))*cossinfacts(1,iterm)
           end do
        end do


     end do
     if (nvctr /= lr%wfd%nvctr_f) then
        write(*,'(1x,a,i0,1x,i0)')' ERROR: nvctr >< nvctr_f ',nvctr,lr%wfd%nvctr_f
        stop 
     end if
     !!$  end if
     
     !now the imaginary part
     
     !!$  if((ithread == 0 .and. nthread <= 2) .or. ithread == 2) then 
     ! Other terms: coarse projector components
     ! coarse part
     do iseg=1,lr%wfd%nseg_c
        call segments_to_grid(lr%wfd%keyvglob(iseg),lr%wfd%keyglob(1,iseg),lr%d,i0,i1,i2,i3,jj)

        do i=i0,i1
           ind_c=i-i0+jj
           do iterm=1,nterm
              psi(ind_c)=psi(ind_c)+re_im_cmplx_prod(&
                   wx(1,1, i,1,iterm),wy(1,1, i2,1,iterm),wz(1,1, i3,1,iterm))*cossinfacts(2,iterm)
           end do
        end do

        do i=i0,i1
           ind_c=lr%wfd%nvctr_c + 7*lr%wfd%nvctr_f + i-i0+jj
           do iterm=1,nterm
              psi(ind_c)=psi(ind_c)+im_im_cmplx_prod(&
                   wx(1,1, i,1,iterm),wy(1,1, i2,1,iterm),wz(1,1, i3,1,iterm))*cossinfacts(2,iterm)
           end do
        end do




     end do

     !!$  end if

     !!$  if((ithread .eq. 1 .and. nthread <=3) .or. nthread .eq. 1 .or. ithread == 3) then
     ! Other terms: fine projector components
     do iseg=lr%wfd%nseg_c+1,lr%wfd%nseg_c+lr%wfd%nseg_f
        call segments_to_grid(lr%wfd%keyvglob(iseg),lr%wfd%keyglob(1,iseg),lr%d,i0,i1,i2,i3,jj)

        do i=i0,i1
           ind_f=lr%wfd%nvctr_c+7*(i-i0+jj-1)
           do iterm=1,nterm
              psi(ind_f+1)=psi(ind_f+1)+re_im_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,1,iterm))*cossinfacts(2,iterm)
              psi(ind_f+2)=psi(ind_f+2)+re_im_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,1,iterm))*cossinfacts(2,iterm)
              psi(ind_f+3)=psi(ind_f+3)+re_im_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,1,iterm))*cossinfacts(2,iterm)
              psi(ind_f+4)=psi(ind_f+4)+re_im_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,2,iterm))*cossinfacts(2,iterm)
              psi(ind_f+5)=psi(ind_f+5)+re_im_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,2,iterm))*cossinfacts(2,iterm)
              psi(ind_f+6)=psi(ind_f+6)+re_im_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,2,iterm))*cossinfacts(2,iterm)
              psi(ind_f+7)=psi(ind_f+7)+re_im_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,2,iterm))*cossinfacts(2,iterm)
           end do
        end do

        do i=i0,i1
           ind_f=lr%wfd%nvctr_c + 7*lr%wfd%nvctr_f +lr%wfd%nvctr_c+7*(i-i0+jj-1)
           do iterm=1,nterm
              psi(ind_f+1)=psi(ind_f+1)+im_im_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,1,iterm))*cossinfacts(2,iterm)
              psi(ind_f+2)=psi(ind_f+2)+im_im_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,1,iterm))*cossinfacts(2,iterm)
              psi(ind_f+3)=psi(ind_f+3)+im_im_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,1,iterm))*cossinfacts(2,iterm)
              psi(ind_f+4)=psi(ind_f+4)+im_im_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,2,iterm))*cossinfacts(2,iterm)
              psi(ind_f+5)=psi(ind_f+5)+im_im_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,1,iterm),wz(1,1,i3,2,iterm))*cossinfacts(2,iterm)
              psi(ind_f+6)=psi(ind_f+6)+im_im_cmplx_prod(&
                   wx(1,1,i,1,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,2,iterm))*cossinfacts(2,iterm)
              psi(ind_f+7)=psi(ind_f+7)+im_im_cmplx_prod(&
                   wx(1,1,i,2,iterm),wy(1,1,i2,2,iterm),wz(1,1,i3,2,iterm))*cossinfacts(2,iterm)
           end do
        end do

     end do

  end if

  !!$omp end parallel

END SUBROUTINE wfn_from_tensprod_cossin




subroutine segments_to_grid(keyv,keyg,grid,i0,i1,i2,i3,jj)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: keyv
  integer, dimension(2), intent(in) :: keyg
  type(grid_dimensions), intent(in) :: grid
  integer, intent(out) :: i0,i1,i2,i3,jj
  !local variables
  integer :: j0,j1,ii

  jj=keyv
  j0=keyg(1)
  j1=keyg(2)
  ii=j0-1
  i3=ii/((grid%n1+1)*(grid%n2+1))
  ii=ii-i3*(grid%n1+1)*(grid%n2+1)
  i2=ii/(grid%n1+1)
  i0=ii-i2*(grid%n1+1)
  i1=i0+j1-j0
END SUBROUTINE segments_to_grid


!temporary creation, better to put in standby
!!!subroutine sumrho_gaussians(geocode,iproc,nproc,norb,norbp,nspin,nspinor,&
!!!     n1i,n2i,n3i,n3d,i3s,hxh,hyh,hzh,G,gaucoeff,occup,spinsgn,rho)
!!!  use module_base
!!!  use module_types
!!!  implicit none
!!!  character(len=1), intent(in) :: geocode
!!!  integer, intent(in) :: iproc,nproc,norb,norbp,n1i,n2i,n3i 
!!!  real(gp), intent(in) :: hx,hy,hz
!!!  type(wavefunctions_descriptors), intent(in) :: wfd
!!!  type(gaussian_basis), intent(in) :: G
!!!  real(wp), dimension(G%ncoeff,norbp), intent(in) :: wfn_gau
!!!  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp), intent(out) :: psi
!!!
!!!  !local variables
!!!  character(len=*), parameter :: subname='gaussians_to_wavelets'
!!!  integer, parameter :: nterm_max=3
!!!  logical :: myorbital
!!!  integer :: i_stat,i_all,ishell,iexpo,icoeff,iat,isat,ng,l,m,iorb,jorb,i,nterm,ierr,ig
!!!  real(dp) :: normdev,tt,scpr
!!!  real(gp) :: rx,ry,rz
!!!  integer, dimension(nterm_max) :: lx,ly,lz
!!!  real(gp), dimension(nterm_max) :: fac_arr
!!!  real(wp), dimension(:), allocatable :: tpsi
!!!
!!!  if(iproc == 0) write(*,'(1x,a)',advance='no')'Writing wavefunctions in wavelet form '
!!!
!!!  allocate(tpsi(wfd%nvctr_c+7*wfd%nvctr_f+ndebug),stat=i_stat)
!!!  call memocc(i_stat,tpsi,'tpsi',subname)
!!!
!!!  !initialize the wavefunction
!!!  call razero((wfd%nvctr_c+7*wfd%nvctr_f)*norbp,psi)
!!!  !this can be changed to be passed only once to all the gaussian basis
!!!  !eks=0.d0
!!!  !loop over the atoms
!!!  ishell=0
!!!  iexpo=1
!!!  icoeff=1
!!!  do iat=1,G%nat
!!!     rx=G%rxyz(1,iat)
!!!     ry=G%rxyz(2,iat)
!!!     rz=G%rxyz(3,iat)
!!!     !loop over the number of shells of the atom type
!!!     do isat=1,G%nshell(iat)
!!!        ishell=ishell+1
!!!        !the degree of contraction of the basis function
!!!        !is the same as the ng value of the createAtomicOrbitals routine
!!!        ng=G%ndoc(ishell)
!!!        !angular momentum of the basis set(shifted for compatibility with BigDFT routines
!!!        l=G%nam(ishell)
!!!        !multiply the values of the gaussian contraction times the orbital coefficient
!!!        do m=1,2*l-1
!!!           call calc_coeff_inguess(l,m,nterm_max,nterm,lx,ly,lz,fac_arr)
           !this kinetic energy is not reliable
!!eks=eks+ek*occup(iorb)*cimu(m,ishell,iat,iorb)
!!!           call crtonewave(geocode,n1,n2,n3,ng,nterm,lx,ly,lz,fac_arr,G%xp(iexpo),G%psiat(iexpo),&
!!!                rx,ry,rz,hx,hy,hz,0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
!!!                wfd%nseg_c,wfd%nvctr_c,wfd%keyg,wfd%keyv,wfd%nseg_f,wfd%nvctr_f,&
!!!                wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1),&
!!!                tpsi(1),tpsi(wfd%nvctr_c+1))
!!!           !sum the result inside the orbital wavefunction
!!!           !loop over the orbitals
!!!           do iorb=1,norb
!!!              if (myorbital(iorb,norb,iproc,nproc)) then
!!!                 jorb=iorb-iproc*norbp
!!!                 do i=1,wfd%nvctr_c+7*wfd%nvctr_f
!!!                    !for this also daxpy BLAS can be used
!!!                    psi(i,jorb)=psi(i,jorb)+wfn_gau(icoeff,jorb)*tpsi(i)
!!!                 end do
!!!              end if
!!!           end do
!!!           icoeff=icoeff+1
!!!        end do
!!!        iexpo=iexpo+ng
!!!     end do
!!!     if (iproc == 0) then
!!!        write(*,'(a)',advance='no') &
!!!             repeat('.',(iat*40)/G%nat-((iat-1)*40)/G%nat)
!!!     end if
!!!  end do
!!!
!!!  call gaudim_check(iexpo,icoeff,ishell,G%nexpo,G%ncoeff,G%nshltot)
!!!
!!!  if (iproc ==0 ) write(*,'(1x,a)')'done.'
!!!  !renormalize the orbitals
!!!  !calculate the deviation from 1 of the orbital norm
!!!  normdev=0.0_dp
!!!  tt=0.0_dp
!!!  do iorb=1,norb
!!!     if (myorbital(iorb,norb,iproc,nproc)) then
!!!        jorb=iorb-iproc*norbp
!!!        call wnrm(wfd%nvctr_c,wfd%nvctr_f,psi(1,jorb),psi(wfd%nvctr_c+1,jorb),scpr) 
!!!        call wscal(wfd%nvctr_c,wfd%nvctr_f,real(1.0_dp/sqrt(scpr),wp),psi(1,jorb),psi(wfd%nvctr_c+1,jorb))
!!!        !print *,'norm of orbital ',iorb,scpr
!!!        tt=max(tt,abs(1.0_dp-scpr))
!!!     end if
!!!  end do
!!!  if (nproc > 1) then
!!!     call MPI_REDUCE(tt,normdev,1,mpidtypd,MPI_MAX,0,MPI_COMM_WORLD,ierr)
!!!  else
!!!     normdev=tt
!!!  end if
!!!  if (iproc ==0 ) write(*,'(1x,a,1pe12.2)')&
!!!       'Deviation from normalization of the imported orbitals',normdev
!!!
!!!  i_all=-product(shape(tpsi))*kind(tpsi)
!!!  deallocate(tpsi,stat=i_stat)
!!!  call memocc(i_stat,i_all,'tpsi',subname)
!!!
!!!
!!!  !Creates charge density arising from the ionic PSP cores
!!!  if (n3pi >0 ) then
!!!
!!!     !conditions for periodicity in the three directions
!!!     perx=(geocode /= 'F')
!!!     pery=(geocode == 'P')
!!!     perz=(geocode /= 'F')
!!!
!!!     call ext_buffers(perx,nbl1,nbr1)
!!!     call ext_buffers(pery,nbl2,nbr2)
!!!     call ext_buffers(perz,nbl3,nbr3)
!!!
!!!     call razero(n1i*n2i*n3pi,pot_ion)
!!!
!!!     do iat=1,nat
!!!        ityp=iatype(iat)
!!!        rx=rxyz(1,iat) 
!!!        ry=rxyz(2,iat)
!!!        rz=rxyz(3,iat)
!!!
!!!        rloc=psppar(0,0,ityp)
!!!        charge=real(nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
!!!        cutoff=10.d0*rloc
!!!
!!!        isx=floor((rx-cutoff)/hxh)
!!!        isy=floor((ry-cutoff)/hyh)
!!!        isz=floor((rz-cutoff)/hzh)
!!!
!!!        iex=ceiling((rx+cutoff)/hxh)
!!!        iey=ceiling((ry+cutoff)/hyh)
!!!        iez=ceiling((rz+cutoff)/hzh)
!!!
!!!        !these nested loops will be used also for the actual ionic forces, to be recalculated
!!!        do i3=isz,iez
!!!           z=real(i3,kind=8)*hzh-rz
!!!           call ind_positions(perz,i3,n3,j3,goz) 
!!!           j3=j3+nbl3+1
!!!           do i2=isy,iey
!!!              y=real(i2,kind=8)*hyh-ry
!!!              call ind_positions(pery,i2,n2,j2,goy)
!!!              do i1=isx,iex
!!!                 x=real(i1,kind=8)*hxh-rx
!!!                 call ind_positions(perx,i1,n1,j1,gox)
!!!                 r2=x**2+y**2+z**2
!!!                 arg=r2/rloc**2
!!!                 xp=exp(-.5d0*arg)
!!!                 if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
!!!                    ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
!!!                    pot_ion(ind)=pot_ion(ind)-xp*charge
!!!                 else if (.not. goz ) then
!!!                    rholeaked=rholeaked+xp*charge
!!!                 endif
!!!              enddo
!!!           enddo
!!!        enddo
!!!
!!!     enddo
!!!
!!!  end if
!!!
!!!  ! Check
!!!  tt=0.d0
!!!  do j3=1,n3pi
!!!     do i2= -nbl2,2*n2+1+nbr2
!!!        do i1= -nbl1,2*n1+1+nbr1
!!!           ind=i1+1+nbl1+(i2+nbl2)*n1i+(j3-1)*n1i*n2i
!!!           tt=tt+pot_ion(ind)
!!!        enddo
!!!     enddo
!!!  enddo
!!!
!!!  tt=tt*hxh*hyh*hzh
!!!  rholeaked=rholeaked*hxh*hyh*hzh
!!!
!!!  !print *,'test case input_rho_ion',iproc,i3start,i3end,n3pi,2*n3+16,tt
!!!
!!!  if (nproc > 1) then
!!!     charges_mpi(1)=tt
!!!     charges_mpi(2)=rholeaked
!!!     call MPI_ALLREDUCE(charges_mpi(1),charges_mpi(3),2,MPI_double_precision,  &
!!!          MPI_SUM,MPI_COMM_WORLD,ierr)
!!!     tt_tot=charges_mpi(3)
!!!     rholeaked_tot=charges_mpi(4)
!!!  else
!!!     tt_tot=tt
!!!     rholeaked_tot=rholeaked
!!!  end if
!!!
!!!  if (iproc.eq.0) write(*,'(1x,a,f26.12,2x,1pe10.3)') &
!!!       'total ionic charge, leaked charge ',tt_tot,rholeaked_tot
!!!
!!!
!!!END SUBROUTINE sumrho_gaussians


!> Parse the output of CP2K to read the basis set information
subroutine gautowav(geocode,iproc,nproc,nat,ntypes,norb,norbp,n1,n2,n3,&
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     nvctr_c,nvctr_f,nseg_c,nseg_f,keyg,keyv,iatype,rxyz,hx,hy,hz,psi) !n(c) occup (arg:l-5)
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
  !n(c) real(gp), dimension(norb), intent(in) :: occup
  !real(gp), intent(out) :: eks
  real(wp), dimension(nvctr_c+7*nvctr_f,norbp), intent(out) :: psi
  !local variables
  character(len=*), parameter :: subname='gautowav'
  logical :: myorbital
  character(len=6) :: string,symbol
  character(len=100) :: line
  integer, parameter :: nterm_max=3
  integer :: ngx,nbx,nst,nend,ng,num,mmax,myshift,i,ipar,ipg,jat
  integer :: iorb,jorb,iat,ityp,l,m,nterm,i_all,i_stat,ibas,ig,iset,jbas,ishell,lmax
  integer :: ierr
  real(dp) :: tt,normdev
  real(gp) :: rx,ry,rz
  real(gp) :: exponent,coefficient,scpr
  integer, dimension(nterm_max) :: lx,ly,lz
  real(gp), dimension(nterm_max) :: fac_arr
  integer, dimension(:), allocatable :: nshell,iorbtmp
  integer, dimension(:,:), allocatable :: nam,ndoc
  real(wp), dimension(:), allocatable :: tpsi,ctmp
  real(gp), dimension(:), allocatable :: psiatn,xp
  real(gp), dimension(:,:,:), allocatable :: contcoeff,expo
  real(wp), dimension(:,:,:,:), allocatable :: cimu


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
  ishell=0
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


!!!  !print the found values
!!!  do ityp=1,ntypes
!!!     do ishell=1,nshell(ityp)
!!!        print *,'ityp=',ityp,'ishell=',ishell,'l=',nam(ishell,ityp),'ndoc=',ndoc(ishell,ityp)
!!!        do ipg=1,ndoc(ishell,ityp)
!!!           print *,'expo=',expo(ipg,ishell,ityp),'coeff=',contcoeff(ipg,ishell,ityp)
!!!        end do
!!!     end do
!!!  end do

!!!  !here we can start calculate the overlap matrix between the different basis functions
!!!  !as an example we can try to calculate the overlap in one shell
!!!  allocate(iw(18),stat=i_stat)
!!!  call memocc(i_stat,product(shape(iw))*kind(iw),'iw','gautowav')
!!!  allocate(rw(6),stat=i_stat)
!!!  call memocc(i_stat,product(shape(rw))*kind(rw),'rw','gautowav')
!!!
!!!  do ityp=1,ntypes
!!!     do ishell=1,nshell(ityp)
!!!        !perform the scalar product internally to the shell
!!!        do m1=1,2*nam(ishell,ityp)+1
!!!           do m2=1,2*nam(ishell,ityp)+1
!!!              call gbasovrlp(expo(1,ishell,ityp),contcoeff(1,ishell,ityp),&
!!!                   expo(1,ishell,ityp),contcoeff(1,ishell,ityp),&
!!!                   ndoc(ishell,ityp),ndoc(ishell,ityp),&
!!!                   nam(ishell,ityp)+1,m1,nam(ishell,ityp)+1,m2,&
!!!                   0.d0,0.d0,0.d0,&
!!!                   18,6,iw,rw,ovrlp)
!!!              if (iproc==0) then
!!!                 print *,ityp,ishell,nam(ishell,ityp),m1,m2,ovrlp!&
!!!                      !contcoeff(1:ndoc(ishell,ityp),ishell,ityp),ovrlp
!!!              end if
!!!           end do
!!!        end do
!!!     end do
!!!  end do
!!!
!!!  i_all=-product(shape(iw))*kind(iw)
!!!  deallocate(iw,stat=i_stat)
!!!  call memocc(i_stat,i_all,'iw','gautowav')
!!!  i_all=-product(shape(rw))*kind(rw)
!!!  deallocate(rw,stat=i_stat)
!!!  call memocc(i_stat,i_all,'rw','gautowav')

!!!subroutine basis_ovrlp(nat,norb,nbx,ngx,lmax,ntypes,nam,ndoc,contcoeff,expo,cimu)
!!!  
!!!  
!!!END SUBROUTINE basis_ovrlp


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

!!!  nst=1
!!!  nend=4
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
                 !if (iproc==0) 
                  write(*,'(1x,a,i0,a)')&
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
              !if (iproc==0) 
                   write(*,'(1x,a,i0,a)')&
                   'Problem in the gaucoeff.dat file, only ',iat ,' atoms processed'
              stop
           else
              cycle store_coeff
           end if
        end if
     end do read_orbitals

  end do store_coeff
  close(36)

!!!  !print the found values
!!!  do iat=1,nat
!!!     ityp=iatype(iat)
!!!     do ishell=1,nshell(ityp)
!!!        do jbas=1,2*nam(ishell,ityp)+1
!!!           print *,iat,ishell,nam(ishell,ityp),jbas,(cimu(jbas,ishell,iat,iorb),iorb=1,norb)
!!!        end do
!!!     end do
!!!  end do

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
!!!           !this kinetic energy is not reliable
!!!           eks=eks+ek*occup(iorb)*cimu(m,ishell,iat,iorb)
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

END SUBROUTINE gautowav


!> Calculate the shift between the spherical harmonics of CP2K and the one of BigDFT
function myshift(symbol)
  implicit none
  character(len=5), intent(in) :: symbol
  integer :: myshift
  myshift=0
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


!> Returns an input guess orbital that is a Gaussian centered at a Wannier center
!! @f$ exp (-1/(2*gau_a^2) *((x-cntrx)^2 + (y-cntry)^2 + (z-cntrz)^2 )) @f$
!! in the arrays psi_c, psi_f
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
  integer, parameter ::nw=32000
  logical :: perx,pery,perz
  integer:: iterm,itp,n_gau,ml1,mu1,ml2,mu2,ml3,mu3,i1,i2,i3,i_all,i_stat,iseg,ii,jj,j0,j1,i0,i
  real(gp) :: gau_a,te
  real(wp), dimension(0:nw,2) :: work
  real(wp), dimension(:,:), allocatable :: wprojx,wprojy,wprojz
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
  allocate(psig_c(nl1_c:nu1_c,nl2_c:nu2_c,nl3_c:nu3_c+ndebug),stat=i_stat)
  call memocc(i_stat,psig_c,'psig_c',subname)
  allocate(psig_f(7,nl1_f:nu1_f,nl2_f:nu2_f,nl3_f:nu3_f+ndebug),stat=i_stat)
  call memocc(i_stat,psig_f,'psig_f',subname)

  !print *,'limits',nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f

  iterm=1
  itp=1
  gau_a=xp(iterm)
  n_gau=lx(itp)
  call gauss_to_daub(hx,fac_arr(itp),rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw,perx)
  n_gau=ly(itp)
  call gauss_to_daub(hy,1.0_gp,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw,pery)
  n_gau=lz(itp)
  call gauss_to_daub(hz,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw,perz)
!$omp parallel default(private) shared(nl3_c,nu3_c,nl2_c,nu2_c,nl1_c,nu1_c,wprojx,wprojy,wprojz) &
!$omp shared(nl3_f,nu3_f,nl2_f,nu2_f,nl1_f,nu1_f,psig_c,psig_f)
  ! First term: coarse projector components
!$omp do
  do i3=nl3_c,nu3_c
     do i2=nl2_c,nu2_c
        do i1=nl1_c,nu1_c
           psig_c(i1,i2,i3)=wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
        enddo
     enddo
  enddo
!$omp enddo
  ! First term: fine projector components
!$omp do
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
!$omp enddo
!$omp end parallel
  do iterm=2,nterm
     gau_a=xp(iterm)
     n_gau=lx(itp)
     call gauss_to_daub(hx,fac_arr(itp),rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw,perx)
     n_gau=ly(itp)
     call gauss_to_daub(hy,1.0_gp,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw,pery)
     n_gau=lz(itp)
     call gauss_to_daub(hz,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw,perz)

!$omp parallel default(private) shared(nl3_c,nu3_c,nl2_c,nu2_c,nl1_c,nu1_c,wprojx,wprojy,wprojz) &
!$omp shared(nl3_f,nu3_f,nl2_f,nu2_f,nl1_f,nu1_f,psig_c,psig_f)
  ! First term: coarse projector components
!$omp do
     do i3=nl3_c,nu3_c
        do i2=nl2_c,nu2_c
           do i1=nl1_c,nu1_c
              psig_c(i1,i2,i3)=psig_c(i1,i2,i3)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
           enddo
        enddo
     enddo
!$omp enddo

     ! First term: fine projector components
!$omp do
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
!$omp enddo
!$omp end parallel

  end do

  do itp=2,ntp

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

!$omp parallel default(private) shared(nl3_c,nu3_c,nl2_c,nu2_c,nl1_c,nu1_c,wprojx,wprojy,wprojz) &
!$omp shared(nl3_f,nu3_f,nl2_f,nu2_f,nl1_f,nu1_f,psig_c,psig_f)
  ! First term: coarse projector components
!$omp do
        do i3=nl3_c,nu3_c
           do i2=nl2_c,nu2_c
              do i1=nl1_c,nu1_c
                 psig_c(i1,i2,i3)=psig_c(i1,i2,i3)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
              enddo
           enddo
        enddo
!$omp enddo
        ! First term: fine projector components
!$omp do
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
!$omp enddo
!$omp end parallel
     end do


  end do


!$omp parallel default(private) shared(nseg_c,keyv_c,keyg_c,n1,n2,psi_c,nseg_f,keyv_f,keyg_f) &
!$omp shared(psi_f,psig_c,psig_f)

  !wavefunction compression

  !itp=0
  ! coarse part
!$omp do
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
!$omp enddo
  !print *,'nvctr_c',itp,mvctr_c

  !itp=0
  ! fine part
!$omp do
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
!$omp enddo
!$omp end parallel
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
  i_all=-product(shape(psig_c))*kind(psig_c)
  deallocate(psig_c,stat=i_stat)
  call memocc(i_stat,i_all,'psig_c',subname)
  i_all=-product(shape(psig_f))*kind(psig_f)
  deallocate(psig_f,stat=i_stat)
  call memocc(i_stat,i_all,'psig_f',subname)

END SUBROUTINE crtonewave
