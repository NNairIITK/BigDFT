subroutine initInputguessConfinement(iproc, nproc, at, Glr, input, lin, lig, rxyz, nscatterarr, tag)
  ! Input wavefunctions are found by a diagonalization in a minimal basis set
  ! Each processors write its initial wavefunctions into the wavefunction file
  ! The files are then read by readwave
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => initInputguessConfinement
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc,nproc
  type(atoms_data),intent(inout) :: at
  type(locreg_descriptors),intent(in) :: Glr
  type(input_variables)::input
  type(linearParameters),intent(inout):: lin
  type(linearInputGuess),intent(inout):: lig
  integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(gp),dimension(3,at%nat),intent(in):: rxyz
  integer,intent(inout):: tag

  ! Local variables
  character(len=*), parameter :: subname='initInputguessConfinement'
  real(gp), dimension(:),allocatable:: locrad
  integer,dimension(:),allocatable:: norbsPerAt
  integer, parameter :: nmax=6,lmax=3,noccmax=2,nelecmax=32
  integer:: ist, iadd, ii, jj, norbtot, istat, iall, iat, nspin_ig, norbat


  ! Nullify the local zone descriptors.
  call nullify_local_zone_descriptors(lig%lzdig)
  call nullify_local_zone_descriptors(lig%lzdGauss)

  ! Allocate some arrays we need for the input guess.
  allocate(locrad(at%nat+ndebug),stat=istat)
  call memocc(istat,locrad,'locrad',subname)
  allocate(norbsPerAt(at%nat), stat=istat)
  call memocc(istat, norbsPerAt, 'norbsPerAt', subname)

  ! Number of localization regions
  lig%lzdig%nlr=at%nat
  lig%lzdGauss%nlr=at%nat

  ! Spin for inputguess orbitals.
  if (input%nspin == 4) then
     nspin_ig=1
  else
     nspin_ig=input%nspin
  end if

  ! Determine how many atomic orbitals we have. Maybe we have to increase this number to more than
  ! its 'natural' value.
  norbat=0
  ist=0
  norbtot=0
  do iat=1,at%nat
      ii=lin%norbsPerType(at%iatype(iat))
      iadd=0
      do 
          ! Count the number of atomic orbitals and increase the number if necessary until we have more
          ! (or equal) atomic orbitals than basis functions per atom.
          jj=1*nint(at%aocc(1,iat))+3*nint(at%aocc(3,iat))+5*nint(at%aocc(7,iat))+7*nint(at%aocc(13,iat))
          if(jj>=ii) then
              ! we have enough atomic orbitals
              exit
          else
              ! add additional orbitals
              iadd=iadd+1
              select case(iadd)
                  case(1) 
                      at%aocc(1,iat)=1.d0
                  case(2) 
                      at%aocc(3,iat)=1.d0
                  case(3) 
                      at%aocc(7,iat)=1.d0
                  case(4) 
                      at%aocc(13,iat)=1.d0
                  case default 
                      write(*,'(1x,a)') 'ERROR: more than 16 basis functions per atom are not possible!'
                      stop
              end select
          end if
      end do
      norbsPerAt(iat)=jj
      norbat=norbat+norbsPerAt(iat)
      norbtot=norbtot+jj
  end do

  ! Nullify the orbitals_data type and then determine its values.
  call nullify_orbitals_data(lig%orbsig)
  !call orbitals_descriptors(iproc, nproc, norbtot, norbtot, 0, &
  !     input%nspin, lin%orbs%nspinor, lin%orbs%nkpts, lin%orbs%kpts, lin%orbs%kwgts, lig%orbsig)
  call orbitals_descriptors_forLinear(iproc, nproc, norbtot, norbtot, 0, &
       input%nspin, lin%orbs%nspinor, lin%orbs%nkpts, lin%orbs%kpts, lin%orbs%kwgts, lig%orbsig)
  call repartitionOrbitals(iproc, nproc, lig%orbsig%norb, lig%orbsig%norb_par, &
       lig%orbsig%norbp, lig%orbsig%isorb_par, lig%orbsig%isorb, lig%orbsig%onWhichMPI)


  ! lzdig%orbs%inWhichLocreg has been allocated in orbitals_descriptors. Since it will again be allcoated
  ! in assignToLocreg2, deallocate it first.
  iall=-product(shape(lig%orbsig%inWhichLocreg))*kind(lig%orbsig%inWhichLocreg)
  deallocate(lig%orbsig%inWhichLocreg,stat=istat)
  call memocc(istat,iall,'lig%orbsig%inWhichLocreg',subname)

  ! Assign the orbitals to the localization regions.
  call assignToLocreg2(iproc, at%nat, lig%lzdig%nlr, input%nspin, norbsPerAt, rxyz, lig%orbsig)

  ! Maybe this could be moved to another subroutine? Or be omitted at all?
  allocate(lig%orbsig%eval(lin%orbs%norb), stat=istat)
  call memocc(istat, lig%orbsig%eval, 'lig%orbsig%eval', subname)
  lig%orbsig%eval=-.5d0


  ! Nullify the locreg_descriptors and then copy Glr to it.
  call nullify_locreg_descriptors(lig%lzdGauss%Glr)
  call copy_locreg_descriptors(Glr, lig%lzdGauss%Glr, subname)

  ! Determine the localization regions.
  call initLocregs(iproc, nproc, at%nat, rxyz, lig%lzdig, lig%orbsig, input, Glr, lin%locrad, lin%locregShape)
  !call initLocregs(iproc, at%nat, rxyz, lin, input, Glr, phi, lphi)
  call copy_locreg_descriptors(Glr, lig%lzdig%Glr, subname)

  ! Determine the localization regions for the atomic orbitals, which have a different localization radius.
  locrad=max(12.d0,maxval(lin%locrad(:)))
  call nullify_orbitals_data(lig%orbsGauss)
  call copy_orbitals_data(lig%orbsig, lig%orbsGauss, subname)
  call initLocregs(iproc, nproc, at%nat, rxyz, lig%lzdGauss, lig%orbsGauss, input, Glr, locrad, lin%locregShape)

  ! Initialize the parameters needed for the orthonormalization of the atomic orbitals.
  !!! Attention: this is initialized for lzdGauss and not for lzdig!
  !!call initCommsOrtho(iproc, nproc, lig%lzdGauss, lig%orbsGauss, lig%orbsGauss%inWhichLocreg, &
  !!     input, lig%op, lig%comon, tag)
  call initCommsOrtho(iproc, nproc, lig%lzdig, lig%orbsig, lig%orbsig%inWhichLocreg, &
       input, lin%locregShape, lig%op, lig%comon, tag)

  ! Initialize the parameters needed for communicationg the potential.
  call copy_locreg_descriptors(Glr, lig%lzdig%Glr, subname)
  call initializeCommunicationPotential(iproc, nproc, nscatterarr, lig%orbsig, lig%lzdig, lig%comgp, &
       lig%orbsig%inWhichLocreg, tag)

  !!! Attention: this is initialized for lzdGauss and not for lzdig!
  !!call initMatrixCompression(iproc, nproc, lig%orbsig, lig%op, lig%mad)
  !!!call initCompressedMatmul2(lig%orbsig%norb, lig%mad%nseg, lig%mad%keyg, lig%mad%nsegmatmul, &
  !!!      lig%mad%keygmatmul, lig%mad%keyvmatmul)
  !!call initCompressedMatmul3(lig%orbsig%norb, lig%mad)

  !call initMatrixCompression(iproc, nproc, lig%orbsig, lig%op, lig%mad)
  call initMatrixCompression(iproc, nproc, lig%lzdig%nlr, lig%orbsig, &
       lig%op%noverlaps, lig%op%overlaps, lig%mad)
  call initCompressedMatmul3(lig%orbsig%norb, lig%mad)


  ! Deallocate the local arrays.
  iall=-product(shape(locrad))*kind(locrad)
  deallocate(locrad,stat=istat)
  call memocc(istat,iall,'locrad',subname)
  iall=-product(shape(norbsPerAt))*kind(norbsPerAt)
  deallocate(norbsPerAt,stat=istat)
  call memocc(istat,iall,'norbsPerAt',subname)

END SUBROUTINE initInputguessConfinement

!>   input guess wavefunction diagonalization
subroutine inputguessConfinement(iproc, nproc, at, &
     comms, Glr, input, lin, orbs, rxyz,denspot, rhopotold,&
     nlpspd, proj, GPU,  &
     tag, lphi, ehart, eexcu, vexcu)
  ! Input wavefunctions are found by a diagonalization in a minimal basis set
  ! Each processors write its initial wavefunctions into the wavefunction file
  ! The files are then read by readwave
  use module_base
  use module_interfaces, exceptThisOne => inputguessConfinement
  use module_types
  use Poisson_Solver
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc
  type(atoms_data), intent(inout) :: at
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(locreg_descriptors), intent(in) :: Glr
  type(communications_arrays), intent(in) :: comms
  type(GPU_pointers), intent(inout) :: GPU
  type(DFT_local_fields), intent(inout) :: denspot
  type(input_variables):: input
  type(linearParameters),intent(inout):: lin
  type(orbitals_data),intent(in):: orbs
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
  real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin),intent(inout) ::  rhopotold
  integer,intent(inout):: tag
  real(8),dimension(max(lin%orbs%npsidim_orbs,lin%orbs%npsidim_comp)),intent(out):: lphi
  real(8),intent(out):: ehart, eexcu, vexcu

  ! Local variables
  type(gaussian_basis):: G !basis for davidson IG
  character(len=*), parameter :: subname='inputguessConfinement'
  integer :: istat,iall,iat,nspin_ig,iorb,nvirt,norbat
  real(gp) :: hxh,hyh,hzh,eks,epot_sum,ekin_sum,eexctX,eproj_sum,eSIC_DC,t1,t2,time,tt,ddot,dsum
  integer, dimension(:,:), allocatable :: norbsc_arr
  !real(wp), dimension(:), allocatable :: potxc
  !real(dp), dimension(:,:), pointer :: rho_p
  real(gp), dimension(:), allocatable :: locrad
  real(wp), dimension(:,:,:), pointer :: psigau
  real(8),dimension(:),allocatable:: lchi, lchi2
  real(8),dimension(:,:),allocatable::  lhchi
  real(8), dimension(:,:,:),allocatable:: ham3
  integer, dimension(:),allocatable:: norbsPerAt, onWhichAtomTemp, mapping, inversemapping
  logical,dimension(:),allocatable:: covered
  !real(8),dimension(:),pointer:: lpot
  integer, parameter :: nmax=6,lmax=3,noccmax=2,nelecmax=32
  logical:: withConfinement, ovrlpx, ovrlpy, ovrlpz
  logical,dimension(:),allocatable:: doNotCalculate, skip
  integer :: ist,jst,jorb,iiAt,i,iadd,ii,jj,ilr,ind1,ind2
  integer :: ldim,gdim,ierr,jlr,kk,iiorb,ndim_lhchi,ii_orbs,ii_comp
  integer :: is1,ie1,is2,ie2,is3,ie3,js1,je1,js2,je2,js3,je3,nlocregPerMPI,jproc,jlrold
  integer:: norbTarget,norbpTemp,isorbTemp, nprocTemp, ncount
  integer,dimension(:),allocatable:: norb_parTemp, onWhichMPITemp
  type(confpot_data), dimension(:), allocatable :: confdatarr
  real(dp),dimension(6) :: xcstr

  !wvl+PAW objects
  integer::iatyp
  type(gaussian_basis),dimension(at%ntypes)::proj_G
  type(paw_objects)::paw

  !nullify paw objects:
  do iatyp=1,at%ntypes
  call nullify_gaussian_basis(proj_G(iatyp))
  end do
  paw%usepaw=0  !not using PAW
  call nullify_paw_objects(paw)

  if (iproc == 0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------- Input Wavefunctions Creation'
  end if

  ! Allocate some arrays we need for the input guess.
  allocate(norbsc_arr(at%natsc+1,input%nspin+ndebug),stat=istat)
  call memocc(istat,norbsc_arr,'norbsc_arr',subname)
  allocate(locrad(at%nat+ndebug),stat=istat)
  call memocc(istat,locrad,'locrad',subname)
  allocate(norbsPerAt(at%nat), stat=istat)
  call memocc(istat, norbsPerAt, 'norbsPerAt', subname)
  allocate(mapping(lin%lig%orbsig%norb), stat=istat)
  call memocc(istat, mapping, 'mapping', subname)
  allocate(covered(lin%lig%orbsig%norb), stat=istat)
  call memocc(istat, covered, 'covered', subname)
  allocate(inversemapping(lin%lig%orbsig%norb), stat=istat)
  call memocc(istat, inversemapping, 'inversemapping', subname)

  ! Number of localization regions.
  lin%lig%lzdig%nlr=at%nat
  lin%lig%lzdGauss%nlr=at%nat

  ! Spin for inputguess orbitals
  if (input%nspin == 4) then
     nspin_ig=1
  else
     nspin_ig=input%nspin
  end if

  ! Determine how many atomic orbitals we have. Maybe we have to increase this number to more than
  ! its 'natural' value.
  norbat=0
  ist=0
  do iat=1,at%nat
      ii=lin%norbsPerType(at%iatype(iat))
      iadd=0
      do 
          ! Count the number of atomic orbitals and increase the number if necessary until we have more
          ! (or equal) atomic orbitals than basis functions per atom.
          jj=1*nint(at%aocc(1,iat))+3*nint(at%aocc(3,iat))+&
               5*nint(at%aocc(7,iat))+7*nint(at%aocc(13,iat))
          if(jj>=ii) then
              ! we have enough atomic orbitals
              exit
          else
              ! add additional orbitals
              iadd=iadd+1
              select case(iadd)
                  case(1) 
                      at%aocc(1,iat)=1.d0
                  case(2) 
                      at%aocc(3,iat)=1.d0
                  case(3) 
                      at%aocc(7,iat)=1.d0
                  case(4) 
                      at%aocc(13,iat)=1.d0
                  case default 
                      write(*,'(1x,a)') 'ERROR: more than 16 basis functions per atom are not possible!'
                      stop
              end select
          end if
      end do
      norbsPerAt(iat)=jj
      norbat=norbat+norbsPerAt(iat)
  end do


  ! Create the atomic orbitals in a Gaussian basis. Deallocate lin%lig%orbsig, since it will be
  ! recalculated in inputguess_gaussian_orbitals.


  ! This array gives a mapping from the 'natural' orbital distribution (i.e. simply counting up the atoms) to
  ! our optimized orbital distribution (determined by in orbs%inwhichlocreg).
  iiorb=0
  covered=.false.
  do iat=1,at%nat
      do iorb=1,norbsPerAt(iat)
          iiorb=iiorb+1
          ! Search the corresponding entry in inwhichlocreg
          do jorb=1,lin%lig%orbsig%norb
              if(covered(jorb)) cycle
              jlr=lin%lig%orbsig%inwhichlocreg(jorb)
              if( lin%lzd%llr(jlr)%locregCenter(1)==rxyz(1,iat) .and. &
                  lin%lzd%llr(jlr)%locregCenter(2)==rxyz(2,iat) .and. &
                  lin%lzd%llr(jlr)%locregCenter(3)==rxyz(3,iat) ) then
                  covered(jorb)=.true.
                  mapping(iiorb)=jorb
                  !if(iproc==0) write(666,*) iiorb, mapping(iiorb)
                  exit
              end if
          end do
      end do
  end do

  ! Inverse mapping
  do iorb=1,lin%lig%orbsig%norb
      do jorb=1,lin%lig%orbsig%norb
          if(mapping(jorb)==iorb) then
              inversemapping(iorb)=jorb
              !if(iproc==0) write(888,*) iorb, inversemapping(iorb)
              exit
          end if
      end do
  end do



  nvirt=0
  call deallocate_orbitals_data(lin%lig%orbsig,subname)
  

  call inputguess_gaussian_orbitals_forLinear(iproc,nproc,lin%lig%orbsig%norb,at,rxyz,nvirt,nspin_ig,&
       lin%lig%lzdig%nlr, norbsPerAt, mapping, &
       orbs,lin%lig%orbsig,norbsc_arr,locrad,G,psigau,eks)
  ! Since inputguess_gaussian_orbitals overwrites lin%lig%orbsig,we again have to assign the correct value (neeed due to
  ! a different orbital distribution.
  call repartitionOrbitals(iproc,nproc,lin%lig%orbsig%norb,lin%lig%orbsig%norb_par,&
       lin%lig%orbsig%norbp,lin%lig%orbsig%isorb_par,lin%lig%orbsig%isorb,lin%lig%orbsig%onWhichMPI)

  ! Maybe this could be moved to another subroutine? Or be omitted at all?
  allocate(lin%lig%orbsig%eval(lin%orbs%norb), stat=istat)
  call memocc(istat, lin%lig%orbsig%eval, 'lin%lig%orbsig%eval', subname)
  lin%lig%orbsig%eval=-.5d0


  !dimension of the wavefunctions
  call wavefunction_dimension(lin%lig%lzdGauss,lin%lig%orbsig)

  ! Copy lin%lig%orbsig to lin%lig%orbsGauss, but keep the size of the orbitals in lin%lig%orbsGauss.
  call deallocate_orbitals_data(lin%lig%orbsGauss,subname)
  ii_orbs=lin%lig%orbsGauss%npsidim_orbs
  ii_comp=lin%lig%orbsGauss%npsidim_comp
  call nullify_orbitals_data(lin%lig%orbsGauss)
  call copy_orbitals_data(lin%lig%orbsig,lin%lig%orbsGauss,subname)
  lin%lig%orbsGauss%npsidim_orbs=ii_orbs
  lin%lig%orbsGauss%npsidim_comp=ii_comp

  ! Allcoate the array holding the orbitals. lchi2 are the atomic orbitals with the larger cutoff, whereas
  ! lchi are the atomic orbitals with the smaller cutoff.
  !print *,'here',lin%lig%orbsGauss%npsidim_orbs,lin%lig%orbsGauss%npsidim_comp
  allocate(lchi2(max(lin%lig%orbsGauss%npsidim_orbs,lin%lig%orbsGauss%npsidim_comp)),stat=istat)
  call memocc(istat,lchi2,'lchi2',subname)
  lchi2=0.d0

  ! Grid spacing on fine grid.
  hxh=.5_gp*input%hx
  hyh=.5_gp*input%hy
  hzh=.5_gp*input%hz

  ! Assign the size of the orbitals to the new variable lpsidimtot.
  !lin%lig%lzdig%lpsidimtot=lin%lig%orbsig%npsidim
  !lin%lig%lzdGauss%lpsidimtot=lin%lig%orbsGauss%npsidim

  ! Transform the atomic orbitals to the wavelet basis.
  lchi2=0.d0
  call gaussians_to_wavelets_new(iproc, nproc, lin%lig%lzdGauss, lin%lig%orbsGauss, input%hx, input%hy, input%hz, G, &
       psigau(1,1,min(lin%lig%orbsGauss%isorb+1, lin%lig%orbsGauss%norb)), lchi2)

  iall=-product(shape(psigau))*kind(psigau)
  deallocate(psigau,stat=istat)
  call memocc(istat,iall,'psigau',subname)

  call deallocate_gwf(G,subname)

  !restore wavefunction dimension
  call wavefunction_dimension(lin%lig%lzdig,lin%lig%orbsig)



  allocate(lchi(max(lin%lig%orbsig%npsidim_orbs,lin%lig%orbsig%npsidim_comp)+ndebug),stat=istat)
  call memocc(istat,lchi,'lchi',subname)
  lchi=0.d0
  !write(*,*) 'iproc, lin%lig%orbsig%npsidim+ndebug', iproc, lin%lig%orbsig%npsidim+ndebug
  !write(*,*) 'iproc, lin%lig%orbsGauss%npsidim', iproc, lin%lig%orbsGauss%npsidim

  ! Transform chi to the localization region. This requires that the localizatin region of lchi2 is larger than that
  ! of lchi.
  ind1=1
  ind2=1
  do iorb=1,lin%lig%orbsGauss%norbp
      !ilr = lin%lig%orbsig%inWhichLocregp(iorb)
      ilr = lin%lig%orbsig%inWhichLocreg(lin%lig%orbsig%isorb+iorb)
      ldim=lin%lig%lzdig%Llr(ilr)%wfd%nvctr_c+7*lin%lig%lzdig%Llr(ilr)%wfd%nvctr_f
      gdim=lin%lig%lzdGauss%llr(ilr)%wfd%nvctr_c+7*lin%lig%lzdGauss%llr(ilr)%wfd%nvctr_f
      call psi_to_locreg2(iproc, nproc, ldim, gdim, lin%lig%lzdig%llr(ilr), lin%lig%lzdGauss%llr(ilr), lchi2(ind1), lchi(ind2))
      ind1=ind1+lin%lig%lzdGauss%llr(ilr)%wfd%nvctr_c+7*lin%lig%lzdGauss%llr(ilr)%wfd%nvctr_f
      ind2=ind2+lin%lig%lzdig%Llr(ilr)%wfd%nvctr_c+7*lin%lig%lzdig%Llr(ilr)%wfd%nvctr_f
  end do
  if(ind1/=lin%lig%orbsGauss%npsidim_orbs+1) then
      write(*,'(2(a,i8),i8)') 'ERROR on process ',iproc,&
           ': ind1/=lin%lig%orbsGauss%npsidim+1',ind1,lin%lig%orbsGauss%npsidim_orbs+1
      stop
  end if
  if(ind2/=lin%lig%orbsig%npsidim_orbs+1) then
      write(*,'(2(a,i8),i8)') 'ERROR on process ',iproc,&
           ': ind2/=lin%lig%orbsig%npsidim+1',ind2,lin%lig%orbsig%npsidim_orbs+1
      stop
  end if

  ! Always use the exact Loewdin method.
  call orthonormalizeAtomicOrbitalsLocalized2(iproc, nproc, 0, lin%nItOrtho, lin%blocksize_pdsyev, &
       lin%blocksize_pdgemm, lin%lig%lzdig, lin%lig%orbsig, lin%lig%comon, &
       lin%lig%op, input, lin%lig%mad, lchi)

  ! Deallocate locrad, which is not used any longer.
  iall=-product(shape(locrad))*kind(locrad)
  deallocate(locrad,stat=istat)
  call memocc(istat,iall,'locrad',subname)

  !change again wavefunction dimension
  call wavefunction_dimension(lin%lig%lzdGauss,lin%lig%orbsig)


  call deallocate_orbitals_data(lin%lig%orbsGauss, subname)


  ! Create the potential. First calculate the charge density.
  if(iproc==0) write(*,'(1x,a)',advance='no') 'Calculating charge density...'

  call sumrho(iproc,nproc,lin%lig%orbsig,lin%lig%lzdGauss,&
       hxh,hyh,hzh,denspot%dpcom%nscatterarr,&
       GPU,at%sym,denspot%rhod,lchi2,denspot%rho_psi,inversemapping)
  call communicate_density(iproc,nproc,input%nspin,hxh,hyh,hzh,lin%lig%lzdGauss,&
       denspot%rhod,denspot%dpcom%nscatterarr,denspot%rho_psi,denspot%rhov)

  if(iproc==0) write(*,'(a)') 'done.'

  !restore wavefunction dimension
  call wavefunction_dimension(lin%lig%lzdig,lin%lig%orbsig)


  if(trim(lin%mixingmethod)=='dens') then
      call dcopy(max(glr%d%n1i*glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
  end if


  iall=-product(shape(lchi2))*kind(lchi2)
  deallocate(lchi2, stat=istat)
  call memocc(istat, iall, 'lchi2',subname)


  call deallocate_local_zone_descriptors(lin%lig%lzdGauss, subname)

  call updatePotential(iproc,nproc,at%geocode,input%ixc,input%nspin,hxh,hyh,hzh,Glr,denspot,ehart,eexcu,vexcu)

!!$  
!!$  if(orbs%nspinor==4) then
!!$     !this wrapper can be inserted inside the poisson solver 
!!$     call PSolverNC(at%geocode,'D',iproc,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,&
!!$          nscatterarr(iproc,1),& !this is n3d
!!$          input%ixc,hxh,hyh,hzh,&
!!$          rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,4)
!!$  else
!!$     !Allocate XC potential
!!$     if (nscatterarr(iproc,2) >0) then
!!$        allocate(potxc(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)*input%nspin+ndebug),stat=istat)
!!$        call memocc(istat,potxc,'potxc',subname)
!!$     else
!!$        allocate(potxc(1+ndebug),stat=istat)
!!$        call memocc(istat,potxc,'potxc',subname)
!!$     end if
!!$
!!$     call XC_potential(at%geocode,'D',iproc,nproc,&
!!$          Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,input%ixc,hxh,hyh,hzh,&
!!$          rhopot,eexcu,vexcu,input%nspin,rhocore,potxc,xcstr)
!!$
!!$
!!$     if( iand(potshortcut,4)==0) then
!!$        call H_potential(at%geocode,'D',iproc,nproc,&
!!$             Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
!!$             rhopot,pkernel,pot_ion,ehart,0.0_dp,.true.)
!!$     endif
!!$
!!$
!!$     !sum the two potentials in rhopot array
!!$     !fill the other part, for spin, polarised
!!$     if (input%nspin == 2) then
!!$        call dcopy(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),rhopot(1),1,&
!!$             rhopot(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)+1),1)
!!$     end if
!!$     !spin up and down together with the XC part
!!$     call axpy(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)*input%nspin,1.0_dp,potxc(1),1,&
!!$          rhopot(1),1)
!!$
!!$
!!$     iall=-product(shape(potxc))*kind(potxc)
!!$     deallocate(potxc,stat=istat)
!!$     call memocc(istat,iall,'potxc',subname)
!!$
!!$  end if

  if(trim(lin%mixingmethod)=='pot') then
      call dcopy(max(glr%d%n1i*glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
  end if



  !call dcopy(lin%lig%orbsig%npsidim,psi,1,hpsi,1)
  if (input%exctxpar == 'OP2P') eexctX = -99.0_gp

  
  ! Set localnorb, i.e. the number of orbitals a given process has in a specific loalization region.
  do ilr=1,lin%lig%lzdig%nlr
      lin%lig%lzdig%Llr(ilr)%localnorb=0
      do iorb=1,lin%lig%orbsig%norbp
          !if(lin%lig%orbsig%inWhichLocregp(iorb)==ilr) then
          if(lin%lig%orbsig%inWhichLocreg(lin%lig%orbsig%isorb+iorb)==ilr) then
              lin%lig%lzdig%Llr(ilr)%localnorb = lin%lig%lzdig%Llr(ilr)%localnorb+1
          end if
      end do
  end do


  ! Post the messages for the communication of the potential.
  !ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
  call allocateCommunicationsBuffersPotential(lin%lig%comgp, subname)
  call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, lin%lig%comgp)

  ! Gather the potential
  !!! IS THIS DONE BY full_local_potential???
  call gatherPotential(iproc, nproc, lin%lig%comgp)

  ! Apply the Hamiltonian for each atom.
  ! onWhichAtomTemp indicates that all orbitals feel the confining potential
  ! centered on atom iat.
  allocate(onWhichAtomTemp(lin%lig%orbsig%norb), stat=istat)
  call memocc(istat,onWhichAtomTemp,'onWhichAtomTemp',subname)
  allocate(doNotCalculate(lin%lig%lzdig%nlr), stat=istat)
  call memocc(istat, doNotCalculate, 'doNotCalculate', subname)
  allocate(skip(lin%lig%lzdig%nlr), stat=istat)
  call memocc(istat, skip, 'skip', subname)


  ! Determine for how many localization regions we need a Hamiltonian application.
  ndim_lhchi=0
  do iat=1,at%nat
      call getIndices(lin%lig%lzdig%llr(iat), is1, ie1, is2, ie2, is3, ie3)
      skip(iat)=.true.
      do jorb=1,lin%lig%orbsig%norbp
          !onWhichAtomTemp(jorb)=iat
          onWhichAtomTemp(lin%lig%orbsig%isorb+jorb)=iat
          !jlr=lin%lig%orbsig%inWhichLocregp(jorb)
          jlr=lin%lig%orbsig%inWhichLocreg(lin%lig%orbsig%isorb+jorb)
          if(lin%lig%orbsig%inWhichlocreg(jorb+lin%lig%orbsig%isorb)/=jlr) stop 'this should not happen'
          call getIndices(lin%lig%lzdig%llr(jlr), js1, je1, js2, je2, js3, je3)
          ovrlpx = ( is1<=je1 .and. ie1>=js1 )
          ovrlpy = ( is2<=je2 .and. ie2>=js2 )
          ovrlpz = ( is3<=je3 .and. ie3>=js3 )
          if(ovrlpx .and. ovrlpy .and. ovrlpz) then
              skip(iat)=.false.
          end if
      end do
      if(.not.skip(iat)) then
          ndim_lhchi=ndim_lhchi+1
      end if
  end do


  allocate(lhchi(max(lin%lig%orbsig%npsidim_orbs,lin%lig%orbsig%npsidim_comp),ndim_lhchi),stat=istat)
  call memocc(istat, lhchi, 'lhchi', subname)
  lhchi=0.d0


  if(iproc==0) write(*,'(1x,a)') 'Hamiltonian application for all atoms. This may take some time.'
  call mpi_barrier(mpi_comm_world,ierr)
  call cpu_time(t1)


  call local_potential_dimensions(lin%lig%lzdig,lin%lig%orbsig,denspot%dpcom%ngatherarr(0,1))

  call full_local_potential(iproc,nproc,lin%lig%orbsig,lin%lig%lzdig,2,&
       denspot%dpcom,denspot%rhov,denspot%pot_full,lin%lig%comgp)

!!$  call full_local_potential(iproc,nproc,&
!!$       lin%lig%lzdig%glr%d%n1i*lin%lig%lzdig%glr%d%n2i*nscatterarr(iproc,2),&
!!$       lin%lig%lzdig%glr%d%n1i*lin%lig%lzdig%glr%d%n2i*lin%lig%lzdig%glr%d%n3i,input%nspin,&
!!$       lin%lig%lzdig%glr%d%n1i*lin%lig%lzdig%glr%d%n2i*nscatterarr(iproc,1)*input%nspin,0,&
!!$       lin%lig%orbsig,lin%lig%lzdig,2,ngatherarr,rhopot,lpot,lin%lig%comgp)

  allocate(lin%lig%lzdig%doHamAppl(lin%lig%lzdig%nlr), stat=istat)
  call memocc(istat, lin%lig%lzdig%doHamAppl, 'lin%lig%lzdig%doHamAppl', subname)
  withConfinement=.true.
  ii=0
  do iat=1,at%nat
      doNotCalculate=.true.
      lin%lig%lzdig%doHamAppl=.false.
      !!call mpi_barrier(mpi_comm_world, ierr)
      call getIndices(lin%lig%lzdig%llr(iat), is1, ie1, is2, ie2, is3, ie3)
      skip(iat)=.true.
      do jorb=1,lin%lig%orbsig%norbp
          !onWhichAtomTemp(jorb)=iat
          onWhichAtomTemp(lin%lig%orbsig%isorb+jorb)=iat
          !jlr=onWhichAtomp(jorb)
          !jlr=lin%lig%orbsig%inWhichLocregp(jorb)
          jlr=lin%lig%orbsig%inWhichLocreg(lin%lig%orbsig%isorb+jorb)
          call getIndices(lin%lig%lzdig%llr(jlr), js1, je1, js2, je2, js3, je3)
          ovrlpx = ( is1<=je1 .and. ie1>=js1 )
          ovrlpy = ( is2<=je2 .and. ie2>=js2 )
          ovrlpz = ( is3<=je3 .and. ie3>=js3 )
          if(ovrlpx .and. ovrlpy .and. ovrlpz) then
              doNotCalculate(jlr)=.false.
              lin%lig%lzdig%doHamAppl(jlr)=.true.
              skip(iat)=.false.
          else
              doNotCalculate(jlr)=.true.
              lin%lig%lzdig%doHamAppl(jlr)=.false.
          end if
      end do
      if(iproc==0) write(*,'(3x,a,i0,a)', advance='no') 'atom ', iat, '... '

      if(.not.skip(iat)) then
          ii=ii+1
          if(lin%nItInguess>0) then
             allocate(confdatarr(lin%lig%orbsig%norbp))
             call define_confinement_data(confdatarr,lin%lig%orbsig,rxyz,at,&
                  input%hx,input%hy,input%hz,lin,lin%lig%lzdig,onWhichAtomTemp)
             call to_zero(lin%lig%orbsig%npsidim_orbs,lhchi(1,ii))
             call LocalHamiltonianApplication(iproc,nproc,at,lin%lig%orbsig,&
                  input%hx,input%hy,input%hz,&
                  lin%lig%lzdig,confdatarr,denspot%dpcom%ngatherarr,denspot%pot_full,lchi,lhchi(1,ii),&
                  ekin_sum,epot_sum,eexctX,eSIC_DC,input%SIC,GPU,&
                  pkernel=denspot%pkernelseq)
             call NonLocalHamiltonianApplication(iproc,at,lin%lig%orbsig,&
                  input%hx,input%hy,input%hz,rxyz,&
                  proj,lin%lig%lzdig,nlpspd,lchi,lhchi(1,ii),eproj_sum,proj_G,paw)
             deallocate(confdatarr)
          end if
      else
          if(iproc==0) write(*,'(3x,a)', advance='no') 'no Hamiltonian application required... '
      end if
      if(iproc==0) write(*,'(a)') 'done.'
  end do


  iall=-product(shape(denspot%pot_full))*kind(denspot%pot_full)
  deallocate(denspot%pot_full, stat=istat)
  call memocc(istat, iall, 'denspot%pot_full', subname)
   if(ii/=ndim_lhchi) then
      write(*,'(a,i0,a,2(a2,i0))') 'ERROR on process ',iproc,': ii/=ndim_lhchi',ii,ndim_lhchi
      stop
  end if
  call mpi_barrier(mpi_comm_world, ierr)
  call cpu_time(t2)
  time=t2-t1
  if (verbose > 2) then
     if(iproc==0) write(*,'(1x,a,es10.3)') 'time for applying potential:', time
  end if

  ! Deallocate the buffers needed for communication the potential.
  call deallocateCommunicationsBuffersPotential(lin%lig%comgp, subname)
  ! Deallocate the parameters needed for the communication of the potential.
  !call deallocate_p2pCommsGatherPot(lin%lig%comgp, subname)
  call deallocate_p2pComms(lin%lig%comgp, subname)



  ! The input guess is possibly performed only with a subset of all processes.
  if(lin%norbsPerProcIG>lin%orbs%norb) then
      norbTarget=lin%orbs%norb
  else
     norbTarget=lin%norbsperProcIG
  end if
  nprocTemp=ceiling(dble(lin%orbs%norb)/dble(norbTarget))
  nprocTemp=min(nprocTemp,nproc)
  if(iproc==0) write(*,'(a,i0,a)') 'The minimization is performed using ', nprocTemp, ' processes.'

  ! Create temporary norb_parTemp, onWhichMPITemp
  allocate(norb_parTemp(0:nprocTemp-1), stat=istat)
  call memocc(istat, norb_parTemp, 'norb_parTemp', subname)
  norb_parTemp=0
  tt=dble(lin%orbs%norb)/dble(nprocTemp)
  ii=floor(tt)
  ! ii is now the number of orbitals that every process has. Distribute the remaining ones.
  norb_parTemp(0:nprocTemp-1)=ii
  kk=lin%orbs%norb-nprocTemp*ii
  norb_parTemp(0:kk-1)=ii+1

  allocate(onWhichMPITemp(lin%orbs%norb), stat=istat)
  call memocc(istat, onWhichMPITemp, 'onWhichMPITemp', subname)
  iiorb=0
  do jproc=0,nprocTemp-1
      do iorb=1,norb_parTemp(jproc)
          iiorb=iiorb+1
          onWhichMPITemp(iiorb)=jproc
      end do
  end do

  ! Calculate the number of different matrices that have to be stored on a given MPI process.
  jlrold=0
  nlocregPerMPI=0
  do jorb=1,lin%orbs%norb
      jlr=lin%orbs%inWhichLocreg(jorb)
      !jproc=lin%orbs%onWhichMPI(jorb)
      jproc=onWhichMPITemp(jorb)
      !if(iproc==0) write(*,'(a,5i7)') 'jorb, jlr, jlrold, jproc, nlocregPerMPI', jorb, jlr, jlrold, jproc, nlocregPerMPI
      if(iproc==jproc) then
          if(jlr/=jlrold) then
              nlocregPerMPI=nlocregPerMPI+1
              jlrold=jlr
          end if
      end if
  end do





  ! Calculate the Hamiltonian matrix.
  call cpu_time(t1)
  allocate(ham3(lin%lig%orbsig%norb,lin%lig%orbsig%norb,nlocregPerMPI), stat=istat)
  call memocc(istat,ham3,'ham3',subname)

  if(lin%nItInguess>0) then
      call getHamiltonianMatrix6(iproc, nproc, nprocTemp, lin%lig%lzdig, lin%lig%orbsig, lin%orbs, &
           onWhichMPITemp, input, lin%lig%orbsig%inWhichLocreg, ndim_lhchi, &
           nlocregPerMPI, lchi, lhchi, skip, lin%lig%mad, lin%memoryForCommunOverlapIG, lin%locregShape, tag, ham3)
  end if


  iall=-product(shape(lhchi))*kind(lhchi)
  deallocate(lhchi, stat=istat)
  call memocc(istat, iall, 'lhchi',subname)


  ! Build the orbitals phi as linear combinations of the atomic orbitals.
  call buildLinearCombinationsLocalized3(iproc, nproc, lin%lig%orbsig, lin%orbs, lin%comms, at, Glr, input, lin%norbsPerType, &
       lin%lig%orbsig%inWhichLocreg, lchi, lphi, rxyz, lin%orbs%inWhichLocreg, lin, lin%lig%lzdig, nlocregPerMPI, tag, ham3, &
       lin%lig%comon, lin%lig%op, lin%lig%mad)
  !call cpu_time(t2)
  !!time=t2-t1
  !!call mpiallred(time, 1, mpi_sum, mpi_comm_world, ierr)
  !!if(iproc==0) write(*,'(1x,a,es10.3)') 'time for "buildLinearCombinations":', time/dble(nproc)

  if(iproc==0) write(*,'(1x,a)') '------------------------------------------------------------- Input guess generated.'


  ! Deallocate all local arrays.

  ! Deallocate all types that are not needed any longer.
  call deallocate_local_zone_descriptors(lin%lig%lzdig, subname)
  call deallocate_orbitals_data(lin%lig%orbsig, subname)
  call deallocate_matrixDescriptors(lin%lig%mad, subname)
  call deallocate_overlapParameters(lin%lig%op, subname)
  !call deallocate_p2pCommsOrthonormality(lin%lig%comon, subname)
  call deallocate_p2pComms(lin%lig%comon, subname)

  ! Deallocate all remaining local arrays.
  iall=-product(shape(norbsc_arr))*kind(norbsc_arr)
  deallocate(norbsc_arr,stat=istat)
  call memocc(istat,iall,'norbsc_arr',subname)

  iall=-product(shape(onWhichAtomTemp))*kind(onWhichAtomTemp)
  deallocate(onWhichAtomTemp, stat=istat)
  call memocc(istat, iall, 'onWhichAtomTemp',subname)
  
  iall=-product(shape(norbsPerAt))*kind(norbsPerAt)
  deallocate(norbsPerAt, stat=istat)
  call memocc(istat, iall, 'norbsPerAt',subname)

  iall=-product(shape(lchi))*kind(lchi)
  deallocate(lchi, stat=istat)
  call memocc(istat, iall, 'lchi',subname)

  iall=-product(shape(doNotCalculate))*kind(doNotCalculate)
  deallocate(doNotCalculate, stat=istat)
  call memocc(istat, iall, 'doNotCalculate',subname)

  iall=-product(shape(skip))*kind(skip)
  deallocate(skip, stat=istat)
  call memocc(istat, iall, 'skip',subname)

  iall=-product(shape(ham3))*kind(ham3)
  deallocate(ham3, stat=istat)
  call memocc(istat, iall, 'ham3',subname)

  iall=-product(shape(norb_parTemp))*kind(norb_parTemp)
  deallocate(norb_parTemp, stat=istat)
  call memocc(istat, iall, 'norb_parTemp',subname)

  iall=-product(shape(onWhichMPITemp))*kind(onWhichMPITemp)
  deallocate(onWhichMPITemp, stat=istat)
  call memocc(istat, iall, 'onWhichMPITemp',subname)

  iall=-product(shape(mapping))*kind(mapping)
  deallocate(mapping, stat=istat)
  call memocc(istat, iall, 'mapping',subname)

  iall=-product(shape(covered))*kind(covered)
  deallocate(covered, stat=istat)
  call memocc(istat, iall, 'covered',subname)

  iall=-product(shape(inversemapping))*kind(inversemapping)
  deallocate(inversemapping, stat=istat)
  call memocc(istat, iall, 'inversemapping',subname)

END SUBROUTINE inputguessConfinement





subroutine orthonormalizeAtomicOrbitalsLocalized2(iproc, nproc, methTransformOverlap, nItOrtho, blocksize_dsyev, &
           blocksize_pdgemm, lzd, orbs, comon, op, input, mad, lchi)

!
! Purpose:
! ========
!  Orthonormalizes the atomic orbitals chi using a Lowedin orthonormalization.
!
! Calling arguments:
!    

use module_base
use module_types
use module_interfaces, exceptThisOne => orthonormalizeAtomicOrbitalsLocalized2
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, methTransformOverlap, nItOrtho, blocksize_dsyev, blocksize_pdgemm
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: orbs
type(input_variables),intent(in):: input
type(p2pComms),intent(inout):: comon
type(overlapParameters),intent(inout):: op
type(matrixDescriptors),intent(in):: mad
real(8),dimension(orbs%npsidim_comp),intent(inout):: lchi

! Local variables
integer:: iorb, jorb, istat, iall, lwork, info, nvctrp, ierr, tag, ilr
real(8),dimension(:,:),allocatable:: ovrlp
character(len=*),parameter:: subname='orthonormalizeAtomicOrbitalsLocalized'


! Initialize the communication parameters.
!tag=5000
!call initCommsOrtho(iproc, nproc, lzd, orbs, orbs%inWhichLocreg, input, op, comon, tag)
allocate(ovrlp(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, ovrlp, 'ovrlp', subname)

call orthonormalizeLocalized(iproc, nproc, methTransformOverlap, nItOrtho, blocksize_dsyev, blocksize_pdgemm, &
     orbs, op, comon, lzd, mad, lchi, ovrlp)

iall=-product(shape(ovrlp))*kind(ovrlp)
deallocate(ovrlp, stat=istat)
call memocc(istat, iall, 'ovrlp', subname)

!call deallocate_overlapParameters(op, subname)
!call deallocate_p2pCommsOrthonormality(comon, subname)

end subroutine orthonormalizeAtomicOrbitalsLocalized2


!!subroutine orthonormalizeCoefficients(orbs, orbsig, coeff)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!type(orbitals_data),intent(in):: orbs, orbsig
!!real(8),dimension(orbsig%norb,orbs%norb),intent(inout):: coeff
!!
!!! Local variables
!!integer:: iorb, jorb, istat, iall, lwork, info
!!real(8),dimension(:),allocatable:: work, eval
!!real(8),dimension(:,:),allocatable:: ovrlp, coeffTemp
!!real(8),dimension(:,:,:),allocatable:: tempArr
!!character(len=*),parameter:: subname='orthonormalizeCoefficients'
!!real(8):: ddot
!!
!!        allocate(ovrlp(orbs%norb,orbs%norb), stat=istat)
!!        call memocc(istat, ovrlp, 'ovrlp', subname)
!!        allocate(eval(orbs%norb), stat=istat)
!!        call memocc(istat, eval, 'eval', subname)
!!        allocate(tempArr(orbs%norb,orbs%norb,2), stat=istat)
!!        call memocc(istat, tempArr, 'tempArr', subname)
!!        allocate(coeffTemp(orbsig%norb,orbs%norb), stat=istat)
!!        call memocc(istat, coeffTemp, 'coeffTemp', subname)
!!
!!
!!        !!! Orthonormalize the coefficient vectors (Gram-Schmidt).
!!        !!do iorb=1,orbs%norb
!!        !!    do jorb=1,iorb-1
!!        !!        tt=ddot(orbsig%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
!!        !!        call daxpy(orbsig%norb, -tt, coeff(1,jorb), 1, coeff(1,iorb), 1)
!!        !!    end do
!!        !!    tt=dnrm2(orbsig%norb, coeff(1,iorb), 1)
!!        !!    call dscal(orbsig%norb, 1/tt, coeff(1,iorb), 1)
!!        !!end do
!!
!!        !!! Orthonormalize the coefficient vectors (Loewdin).
!!        !!do iorb=1,orbs%norb
!!        !!    do jorb=1,orbs%norb
!!        !!        ovrlp(iorb,jorb)=ddot(orbsig%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
!!        !!    end do
!!        !!end do
!!
!!        allocate(work(1), stat=istat)
!!        call memocc(istat, work, 'work', subname)
!!        call dsyev('v', 'l', orbs%norb, ovrlp(1,1), orbs%norb, eval, work, -1, info)
!!        lwork=work(1)
!!        iall=-product(shape(work))*kind(work)
!!        deallocate(work, stat=istat)
!!        call memocc(istat, iall, 'work', subname)
!!        allocate(work(lwork), stat=istat)
!!        call memocc(istat, work, 'work', subname)
!!        call dsyev('v', 'l', orbs%norb, ovrlp(1,1), orbs%norb, eval, work, lwork, info)
!!        iall=-product(shape(work))*kind(work)
!!        deallocate(work, stat=istat)
!!        call memocc(istat, iall, 'work', subname)
!!
!!        ! Calculate S^{-1/2}. 
!!        ! First calulate ovrlp*diag(1/sqrt(evall)) (ovrlp is the diagonalized overlap
!!        ! matrix and diag(1/sqrt(evall)) the diagonal matrix consisting of the inverse square roots of the eigenvalues...
!!        do iorb=1,orbs%norb
!!            do jorb=1,orbs%norb
!!                tempArr(jorb,iorb,1)=ovrlp(jorb,iorb)*1.d0/sqrt(eval(iorb))
!!            end do
!!        end do
!!
!!        ! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
!!        ! This will give S^{-1/2}.
!!        call dgemm('n', 't', orbs%norb, orbs%norb, orbs%norb, 1.d0, ovrlp(1,1), &
!!        orbs%norb, tempArr(1,1,1), orbs%norb, 0.d0, &
!!        tempArr(1,1,2), orbs%norb)
!!
!!        ! Now calculate the orthonormal orbitals by applying S^{-1/2} to the orbitals.
!!        ! This requires the use of a temporary variable phidTemp.
!!        call dgemm('n', 'n', orbsig%norb, orbs%norb, orbs%norb, 1.d0, coeff(1,1), &
!!             orbsig%norb, tempArr(1,1,2), orbs%norb, 0.d0, &
!!             coeffTemp(1,1), orbsig%norb)
!!        
!!        ! Now copy the orbitals from the temporary variable to phid.
!!        call dcopy(orbs%norb*orbsig%norb, coeffTemp(1,1), 1, coeff(1,1), 1)
!!
!!        iall=-product(shape(ovrlp))*kind(ovrlp)
!!        deallocate(ovrlp, stat=istat)
!!        call memocc(istat, iall, 'ovrlp', subname)
!!
!!        iall=-product(shape(eval))*kind(eval)
!!        deallocate(eval, stat=istat)
!!        call memocc(istat, iall, 'eval', subname)
!!
!!        iall=-product(shape(tempArr))*kind(tempArr)
!!        deallocate(tempArr, stat=istat)
!!        call memocc(istat, iall, 'tempArr', subname)
!!
!!        iall=-product(shape(coeffTemp))*kind(coeffTemp)
!!        deallocate(coeffTemp, stat=istat)
!!        call memocc(istat, iall, 'coeffTemp', subname)
!!
!!end subroutine orthonormalizeCoefficients



subroutine initializeInguessParameters(iproc, orbs, orbsig, newComm, ip)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc
type(orbitals_data),intent(in):: orbs, orbsig
integer,intent(in):: newComm
type(inguessParameters),intent(inout):: ip

! Local variables
integer:: ii, kk, jproc, istat, ierr, norbTarget, iorb, iiorb
real(8):: tt
character(len=*),parameter:: subname='initializeInguessParameters'


  ip%norb=orbs%norb
  ip%norbtot=orbsig%norb

  ! In order to symplify the transposing/untransposing, the orbitals are padded with zeros such that 
  ! they can be distributed evenly over all processes when being transposed. The new length of the 
  ! orbitals after this padding is then given by ip%norbtotPad.
  ip%norbtotPad=ip%norbtot
  do
      if(mod(ip%norbtotPad, ip%nproc)==0) exit
      ip%norbtotPad=ip%norbtotPad+1
  end do

  ! Distribute the orbitals among the processes.
  allocate(ip%norb_par(0:ip%nproc-1), stat=istat)
  call memocc(istat, ip%norb_par, 'ip%norb_par', subname)
  ip%norb_par=0
  tt=dble(ip%norb)/dble(ip%nproc)
  ii=floor(tt)
  ! ii is now the number of orbitals that every process has. Distribute the remaining ones.
  ip%norb_par(0:ip%nproc-1)=ii
  kk=ip%norb-ip%nproc*ii
  ip%norb_par(0:kk-1)=ii+1

  ! ip%onWhichMPI indicates on which MPI process a given orbital is located.
  allocate(ip%onWhichMPI(ip%norb), stat=istat)
  call memocc(istat, ip%onWhichMPI, 'ip%onWhichMPI', subname)
  iiorb=0
  do jproc=0,ip%nproc-1
      do iorb=1,ip%norb_par(jproc)
          iiorb=iiorb+1
          ip%onWhichMPI(iiorb)=jproc
      end do
  end do


  ! Starting orbital for each process
  allocate(ip%isorb_par(0:ip%nproc-1), stat=istat)
  call memocc(istat, ip%isorb_par, 'ip%isorb_par', subname)
  ip%isorb_par=0
  ii=0
  do jproc=0,iproc-1
     ii=ii+ip%norb_par(jproc)
  end do
  !reference orbital for process
  ip%isorb=ii
  ip%isorb_par(iproc)=ip%isorb
  call mpiallred(ip%isorb_par(0), ip%nproc, mpi_sum, newComm, ierr)
  


  ! Calculate the number of elements that each process has when the vectors are transposed.
  ! nvctrp is the total number, nvctrp_nz is the nonzero numbers.
  allocate(ip%nvctrp_nz(0:ip%nproc-1), stat=istat)
  call memocc(istat, ip%nvctrp_nz, 'ip%nvctrp_nz', subname)
  tt=ip%norbtot/dble(ip%nproc)
  ii=floor(tt)
  ! ii is now the number of elements that every process has. Distribute the remaining ones.
  ip%nvctrp_nz=ii
  kk=ip%norbtot-ip%nproc*ii
  ip%nvctrp_nz(0:kk-1)=ii+1
  ! Check wheter this distribution is correct
  ii=0
  do jproc=0,ip%nproc-1
     ii=ii+ip%nvctrp_nz(jproc)
  end do
  if(ii/=ip%norbtot) then
     if(iproc==0) write(*,'(3x,a)') 'ERROR: wrong partition of ip%norbtot!'
     call mpi_barrier(newComm, ierr)
     stop
  end if

  ! With the padded zeros, the elements can be distributed evenly.
  ip%nvctrp=ip%norbtotPad/ip%nproc

  ! Define the values for the mpi_alltoallv.
  ! sendcounts: number of elements that a given  process sends to another process.
  ! senddispls: offset of the starting index on a given process for the send operation to another process.
  allocate(ip%sendcounts(0:ip%nproc-1), stat=istat)
  call memocc(istat, ip%sendcounts, 'ip%sendcounts', subname)
  allocate(ip%senddispls(0:ip%nproc-1), stat=istat)
  call memocc(istat, ip%senddispls, 'ip%senddispls', subname)
  ii=0
  do jproc=0,ip%nproc-1
      ip%sendcounts(jproc)=ip%nvctrp*ip%norb_par(iproc)
      ip%senddispls(jproc)=ii
      ii=ii+ip%sendcounts(jproc)
  end do
  ! recvcounts: number of elements that a given process receives from another process.
  ! recvdispls: offset of the starting index on a given process for the receive operation from another process.
  allocate(ip%recvcounts(0:ip%nproc-1), stat=istat)
  call memocc(istat, ip%recvcounts, 'ip%recvcounts', subname)
  allocate(ip%recvdispls(0:ip%nproc-1), stat=istat)
  call memocc(istat, ip%recvdispls, 'ip%recvdispls', subname)
  ii=0
  do jproc=0,ip%nproc-1
      ip%recvcounts(jproc)=ip%nvctrp*ip%norb_par(jproc)
      ip%recvdispls(jproc)=ii
      ii=ii+ip%recvcounts(jproc)
  end do

  ! Determine the size of the work array needed for the transposition.
  ip%sizeWork=max(ip%norbtotPad*ip%norb_par(iproc),sum(ip%recvcounts(:)))

end subroutine initializeInguessParameters




subroutine getHamiltonianMatrix6(iproc, nproc, nprocTemp, lzdig, orbsig, orbs, onWhichMPITemp, &
input, onWhichAtom, ndim_lhchi, nlocregPerMPI, lchi, lhchi, skip, mad, memoryForCommunOverlapIG, locregShape, &
tagout, ham)
use module_base
use module_types
use module_interfaces, exceptThisOne => getHamiltonianMatrix6
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nprocTemp, ndim_lhchi, nlocregPerMPI
type(local_zone_descriptors),intent(in):: lzdig
type(orbitals_data),intent(in):: orbsig, orbs
integer,dimension(orbs%norb),intent(in):: onWhichMPITemp
type(input_variables),intent(in):: input
integer,dimension(orbsig%norb),intent(in):: onWhichAtom
!real(8),dimension(lzdig%orbs%npsidim),intent(in):: chi
!real(8),dimension(lzdig%orbs%npsidim,nat),intent(in):: hchi
real(8),dimension(max(orbsig%npsidim_orbs,orbsig%npsidim_comp)),intent(in):: lchi
real(8),dimension(max(orbsig%npsidim_orbs,orbsig%npsidim_comp),ndim_lhchi),intent(in):: lhchi
logical,dimension(lzdig%nlr),intent(in):: skip
type(matrixDescriptors),intent(in):: mad
integer,intent(in):: memoryForCommunOverlapIG
character(len=1),intent(in):: locregShape
integer,intent(inout):: tagout
!logical,dimension(lin%lig%lzdig%nlr,0:nproc-1),intent(in):: skipGlobal
real(8),dimension(orbsig%norb,orbsig%norb,nlocregPerMPI),intent(out):: ham

! Local variables
integer:: sizeChi, istat, iorb, ilr, iall, ind1, ind2, ldim, gdim, iat, jproc, ilrold, iilr, iatold, iiorb, jlr, ii
integer:: jorb, ierr, noverlaps, iiat, iioverlap, ioverlap, tagx, availableMemory, jj, i, ist, jst, nshift
integer:: irecv, isend, nrecv, nsend, tag, tag0, jjproc, ind, imat, imatold, jjprocold
type(overlapParameters):: op
type(p2pComms):: comon
real(8),dimension(:,:),allocatable:: hamTemp
character(len=*),parameter:: subname='getHamiltonianMatrix6'
real(8),dimension(:,:),allocatable:: hamTempCompressed, hamTempCompressed2
integer,dimension(:),allocatable:: displs, sendcounts, sendrequests, recvrequests
real(8):: tt1, tt2, tt3


call nullify_p2pcomms(comon) 

allocate(sendcounts(0:nproc-1), stat=istat)
call memocc(istat, sendcounts, 'sendcounts', subname)
allocate(displs(0:nproc-1), stat=istat)
call memocc(istat, displs, 'displs', subname)

call getCommunArraysMatrixCompression(iproc, nproc, orbsig, mad, sendcounts, displs)
availableMemory=memoryForCommunOverlapIG*1048576
availableMemory=availableMemory/8 ! double precision
ii=maxval(sendcounts)
noverlaps=max(availableMemory/ii,1)
if(iproc==0) write(*,'(1x,a,i0,a)') 'the specified memory allows to overlap ', noverlaps,' iterations with communication'
noverlaps=min(noverlaps,lzdig%nlr)


allocate(hamTempCompressed(sendcounts(iproc),noverlaps), stat=istat)
call memocc(istat, hamTempCompressed, 'hamTempCompressed', subname)
allocate(hamTempCompressed2(mad%nvctr,nlocregPerMPI), stat=istat)
call memocc(istat, hamTempCompressed2, 'ovrlpCompressed2', subname)

allocate(hamTemp(orbsig%norb,orbsig%norb), stat=istat)
call memocc(istat, hamTemp, 'hamTemp', subname)

! Initialize the parameters for calculating the matrix.
call nullify_p2pComms(comon)
call initCommsOrtho(iproc, nproc, lzdig, orbsig, onWhichAtom, input, locregShape, op, comon, tagout)


call allocateCommuncationBuffersOrtho(comon, subname)

! Put lphi in the sendbuffer, i.e. lphi will be sent to other processes' receive buffer.
! Then post the messages and gather them.
!call extractOrbital2(iproc, nproc, orbsig, orbsig%npsidim, onWhichAtom, lzdig, op, lchi, comon)
call extractOrbital3(iproc, nproc, orbsig, orbsig%npsidim_orbs, onWhichAtom, lzdig, op, lchi, comon%nsendBuf, comon%sendBuf)
!!call postCommsOverlap(iproc, nproc, comon)
call postCommsOverlapNew(iproc, nproc, orbsig, op, lzdig, lchi, comon, tt1, tt2)
!call gatherOrbitals2(iproc, nproc, comon)
call collectnew(iproc, nproc, comon, mad, op, orbsig, lzdig, comon%nsendbuf, &
     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, tt1, tt2, tt3)



if(iproc==0) write(*,'(1x,a)') 'Calculating Hamiltonian matrix for all atoms. This may take some time.'
ilrold=0
iilr=0
iatold=0
ii=0
imatold=1
imat=0
do iat=1,lzdig%nlr

    if(iproc==0) write(*,'(3x,a,i0,a)', advance='no') 'Calculating matrix for atom ', iat, '... '

    ioverlap=mod(iat-1,noverlaps)+1


    ! Put lhphi to the sendbuffer, so we can the calculate <lphi|lhphi>
    if(.not.skip(iat)) then
        ii=ii+1
        call extractOrbital3(iproc, nproc, orbsig, orbsig%npsidim_orbs, onWhichAtom, lzdig, op, &
             lhchi(1,ii), comon%nsendBuf, comon%sendBuf)
        call calculateOverlapMatrix3Partial(iproc, nproc, orbsig, op, onWhichAtom, comon%nsendBuf, comon%sendBuf, &
             comon%nrecvBuf, comon%recvBuf, mad, hamTemp(1,1))
        !!do istat=1,orbsig%norb
        !!    do iall=1,orbsig%norb
        !!        write(40000+1000*iproc+iat,*) istat,iall,hamTemp(iall,istat)
        !!    end do
        !!end do
        call compressMatrixPerProcess(iproc, nproc, orbsig, mad, hamTemp, sendcounts(iproc), hamTempCompressed(1,ioverlap))

    else
        call razero(sendcounts(iproc), hamTempCompressed(1,ioverlap))
    end if
    if(iproc==0) write(*,'(a)') 'done.'

    
    if(ioverlap==noverlaps .or. iat==lzdig%nlr) then
        call timing(iproc,'ig_matric_comm','ON')
        
        ! Communicate the matrices calculated so far.
        if(iproc==0) write(*,'(1x,a)',advance='no') 'communicating matrices...'

        ! jj indicates how many matrices ar eto be communicated.
        jj=mod(iat-1,noverlaps)+1
        if(iproc==0) write(*,*) 'jj', jj

        ! nshift indicates how much the following loops do i=1,jj deviate from the outer loop on iat.
        nshift=iat-jj

        ! First determine the number of sends / receives for each process.
        nsend=0
        nrecv=0
        ilrold=-1
        jjprocold=-1
        do iorb=1,orbs%norb
            ilr=orbs%inWhichLocreg(iorb)
            jjproc=onWhichMPITemp(iorb)
            do iioverlap=1,jj
                iiat=iioverlap+nshift
                if(iproc<nprocTemp) then
                    if(ilr==ilrold .and. jjproc==jjprocold) cycle !Otherwise we would communicate the same again
                    if(ilr==iiat) then
                        ! Send this matrix to process jproc.
                        if(iproc==jjproc) then
                            do jproc=0,nproc-1
                                nrecv=nrecv+1
                            end do
                            nsend=nsend+1
                        else
                            nsend=nsend+1
                        end if
                    end if
                end if
            end do
            ilrold=ilr
            jjprocold=jjproc
        end do


        allocate(sendrequests(nsend), stat=istat)
        call memocc(istat, sendrequests, 'sendrequests', subname)
        allocate(recvrequests(nrecv), stat=istat)
        call memocc(istat, recvrequests, 'recvrequests', subname)

        ! Now communicate the matrices.
        tag0=1
        isend=0
        irecv=0
        ilrold=-1
        jjprocold=-1
        do iorb=1,orbs%norb
            ilr=orbs%inWhichLocreg(iorb)
            jjproc=onWhichMPITemp(iorb)
            do iioverlap=1,jj
                iiat=iioverlap+nshift
                ! Check whether this MPI needs this matrix. Since only nprocTemp processes will be involved
                ! in calculating the input guess, this check has to be done only for those processes.
                if(iproc<nprocTemp) then
                    if(ilr==ilrold .and. jjproc==jjprocold) cycle
                    if(ilr==iiat) then
                        ! Send to process jproc
                       if(iproc==jjproc .and. nproc > 1) then
                          imat=imat+1
                          do jproc=0,nproc-1
                             tag=tag0+jproc
                             irecv=irecv+1
                             !write(*,'(3(a,i0))') 'process ',iproc,' receives data from process ',jproc,' with tag ',tag
                             call mpi_irecv(hamTempCompressed2(displs(jproc)+1,imat), sendcounts(jproc), &
                                  mpi_double_precision, jproc, tag, mpi_comm_world, recvrequests(irecv), ierr)
                          end do
                          tag=tag0+iproc
                          isend=isend+1
                          !write(*,'(3(a,i0))') 'process ',iproc,' sends data to process ',jjproc,' with tag ',tag
                          call mpi_isend(hamTempCompressed(1,iioverlap), sendcounts(iproc), &
                               mpi_double_precision, jjproc, tag, mpi_comm_world, sendrequests(isend), ierr)
                       !else if (nproc ==1) then
                       !   call vcopy(sendcounts(iproc),hamTempCompressed(1,iioverlap),1,&
                        !       hamTempCompressed2(displs(jproc)+1,imat),1)
                       else if (nproc >1) then
                          tag=tag0+iproc
                          isend=isend+1
                          !write(*,'(3(a,i0))') 'Aprocess ',iproc,' sends data to process ',jjproc,' with tag ',tag
                          call mpi_isend(hamTempCompressed(1,iioverlap), sendcounts(iproc), &
                               mpi_double_precision, jjproc, tag, mpi_comm_world, sendrequests(isend), ierr)
                       else if (nproc == 1) then
                          imat=imat+1
                          call vcopy(sendcounts(iproc),hamTempCompressed(1,iioverlap),1,&
                               hamTempCompressed2(displs(iproc)+1,imat),1)
                       end if
                        tag0=tag0+1
                    end if
                end if
            end do
            ilrold=ilr
            jjprocold=jjproc
        end do
        
        ! Wait for the communication to complete
        if (nproc > 1) then
           isend=0
           waitForSend: do
              call mpi_waitany(nsend-isend, sendrequests(1), ind, mpi_status_ignore, ierr)
              isend=isend+1
              do i=ind,nsend-isend
                 sendrequests(i)=sendrequests(i+1)
              end do
              if(isend==nsend) exit waitForSend
           end do waitForSend

           irecv=0
           waitForRecv: do
              call mpi_waitany(nrecv-irecv, recvrequests(1), ind, mpi_status_ignore, ierr)
              irecv=irecv+1
              do i=ind,nrecv-irecv
                 recvrequests(i)=recvrequests(i+1)
              end do
              if(irecv==nrecv) exit waitForrecv
           end do waitForRecv
        end if
     ! Uncompress the matrices
     do i=imatold,imat
        !call uncompressMatrix(orbs%norb, mad, hamTempCompressed2(1,i), ham(1,1,i))
        call uncompressMatrix(orbsig%norb, mad, hamTempCompressed2(1,i), ham(1,1,i))
     end do
     imatold=imat+1

     iall=-product(shape(sendrequests))*kind(sendrequests)
     deallocate(sendrequests, stat=istat)
     call memocc(istat, iall, 'sendrequests', subname)
     iall=-product(shape(recvrequests))*kind(recvrequests)
     deallocate(recvrequests, stat=istat)
     call memocc(istat, iall, 'recvrequests', subname)

     if(iproc==0) write(*,'(a)') ' done.'

     call timing(iproc,'ig_matric_comm','OF')

  end if

end do



call mpi_barrier(mpi_comm_world, ierr)


!!if(imat/=ndim_lhchi) then
!!    write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': imat/=ndim_lhchi',imat,ndim_lhchi
!!    stop
!!end if

if(imat/=nlocregPerMPI .and. nproc >1) then
  write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': imat/=nlocregPerMPI',imat,nlocregPerMPI
  stop
end if
call deallocate_overlapParameters(op, subname)
call deallocate_p2pComms(comon, subname)
!call deallocateCommuncationBuffersOrtho(comon, subname)

iall=-product(shape(hamTempCompressed))*kind(hamTempCompressed)
deallocate(hamTempCompressed, stat=istat)
call memocc(istat, iall, 'hamTempCompressed', subname)
iall=-product(shape(hamTempCompressed2))*kind(hamTempCompressed2)
deallocate(hamTempCompressed2, stat=istat)
call memocc(istat, iall, 'hamTempCompressed2', subname)
iall=-product(shape(sendcounts))*kind(sendcounts)
deallocate(sendcounts, stat=istat)
call memocc(istat, iall, 'sendcounts', subname)
iall=-product(shape(displs))*kind(displs)
deallocate(displs, stat=istat)
call memocc(istat, iall, 'displs', subname)

iall=-product(shape(hamTemp))*kind(hamTemp)
deallocate(hamTemp, stat=istat)
call memocc(istat, iall, 'hamTemp', subname)

end subroutine getHamiltonianMatrix6





subroutine determineLocalizationRegions(iproc, nproc, nlr, norb, at, onWhichAtomAll, locrad, rxyz, lzd, hx, hy, hz, mlr)
use module_base
use module_types
use module_interfaces, exceptThisOne => determineLocalizationRegions
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nlr, norb
type(atoms_data),intent(in):: at
integer,dimension(norb),intent(in):: onWhichAtomAll
real(8),dimension(at%nat),intent(in):: locrad
real(8),dimension(3,at%nat),intent(in):: rxyz
type(local_zone_descriptors),intent(in):: lzd
real(8),intent(in):: hx, hy, hz
type(matrixLocalizationRegion),dimension(:),pointer,intent(out):: mlr

! Local variables
integer:: ilr, jlr, jorb, ii, istat, is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3
real(8):: cut, tt
logical:: ovrlpx, ovrlpy, ovrlpz
character(len=*),parameter:: subname='determineLocalizationRegions'


allocate(mlr(nlr), stat=istat)
do ilr=1,nlr
  call nullify_matrixLocalizationRegion(mlr(ilr))
end do

!!!! THIS WAS THE ORIGINAL
!!! Count for each localization region the number of matrix elements within the cutoff.
!!do ilr=1,nlr
!!    mlr(ilr)%norbinlr=0
!!    call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
!!    do jorb=1,norb
!!        jlr=onWhichAtomAll(jorb)
!!        call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
!!        ovrlpx = ( is1<=je1 .and. ie1>=js1 )
!!        ovrlpy = ( is2<=je2 .and. ie2>=js2 )
!!        ovrlpz = ( is3<=je3 .and. ie3>=js3 )
!!        if(ovrlpx .and. ovrlpy .and. ovrlpz) then
!!            mlr(ilr)%norbinlr=mlr(ilr)%norbinlr+1
!!        end if
!!    end do
!!    !if(iproc==0) write(*,'(a,2i8)') 'ilr, mlr(ilr)%norbinlr', ilr, mlr(ilr)%norbinlr
!!    allocate(mlr(ilr)%indexInGlobal(mlr(ilr)%norbinlr), stat=istat)
!!    call memocc(istat, mlr(ilr)%indexInGlobal, 'mlr(ilr)%indexInGlobal', subname)
!!    !if(iproc==0) write(*,'(a,i4,i7)') 'ilr, mlr(ilr)%norbinlr', ilr, mlr(ilr)%norbinlr
!!end do
!!
!!
!!! Now determine the indices of the elements with an overlap.
!!do ilr=1,nlr
!!    ii=0
!!    call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
!!    do jorb=1,norb
!!        jlr=onWhichAtomAll(jorb)
!!        call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
!!        ovrlpx = ( is1<=je1 .and. ie1>=js1 )
!!        ovrlpy = ( is2<=je2 .and. ie2>=js2 )
!!        ovrlpz = ( is3<=je3 .and. ie3>=js3 )
!!        if(ovrlpx .and. ovrlpy .and. ovrlpz) then
!!            ii=ii+1
!!            mlr(ilr)%indexInGlobal(ii)=jorb
!!            !if(iproc==0) write(*,'(a,3i8)') 'ilr, ii, mlr(ilr)%indexInGlobal(ii)', ilr, ii, mlr(ilr)%indexInGlobal(ii)
!!        end if
!!    end do
!!    if(ii/=mlr(ilr)%norbinlr) then
!!        write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ', iproc, ': ii/=mlr(ilr)%norbinlr', ii, mlr(ilr)%norbinlr
!!    end if
!!    !if(iproc==0) write(*,'(a,i6,200i5)') 'ilr, mlr(ilr)%indexInGlobal(ii)', ilr, mlr(ilr)%indexInGlobal(:)
!!end do

!! THIS IS NEW
! Count for each localization region the number of matrix elements within the cutoff.
do ilr=1,nlr
  mlr(ilr)%norbinlr=0
  call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
  do jorb=1,norb
     jlr=onWhichAtomAll(jorb)
     !call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
     js1=floor(rxyz(1,jlr)/hx)
     je1=ceiling(rxyz(1,jlr)/hx)
     js2=floor(rxyz(2,jlr)/hy)
     je2=ceiling(rxyz(2,jlr)/hy)
     js3=floor(rxyz(3,jlr)/hz)
     je3=ceiling(rxyz(3,jlr)/hz)
     ovrlpx = ( is1<=je1 .and. ie1>=js1 )
     ovrlpy = ( is2<=je2 .and. ie2>=js2 )
     ovrlpz = ( is3<=je3 .and. ie3>=js3 )
     if(ovrlpx .and. ovrlpy .and. ovrlpz) then
        mlr(ilr)%norbinlr=mlr(ilr)%norbinlr+1
     end if
  end do
  !if(iproc==0) write(*,'(a,2i8)') 'ilr, mlr(ilr)%norbinlr', ilr, mlr(ilr)%norbinlr
  allocate(mlr(ilr)%indexInGlobal(mlr(ilr)%norbinlr), stat=istat)
  call memocc(istat, mlr(ilr)%indexInGlobal, 'mlr(ilr)%indexInGlobal', subname)
  !if(iproc==0) write(*,'(a,i4,i7)') 'ilr, mlr(ilr)%norbinlr', ilr, mlr(ilr)%norbinlr
end do


! Now determine the indices of the elements with an overlap.
do ilr=1,nlr
  ii=0
  call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
  do jorb=1,norb
     jlr=onWhichAtomAll(jorb)
     !call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
     js1=floor(rxyz(1,jlr)/hx)
     je1=ceiling(rxyz(1,jlr)/hx)
     js2=floor(rxyz(2,jlr)/hy)
     je2=ceiling(rxyz(2,jlr)/hy)
     js3=floor(rxyz(3,jlr)/hz)
     je3=ceiling(rxyz(3,jlr)/hz)
     ovrlpx = ( is1<=je1 .and. ie1>=js1 )
     ovrlpy = ( is2<=je2 .and. ie2>=js2 )
     ovrlpz = ( is3<=je3 .and. ie3>=js3 )
     if(ovrlpx .and. ovrlpy .and. ovrlpz) then
        ii=ii+1
        mlr(ilr)%indexInGlobal(ii)=jorb
        !if(iproc==0) write(*,'(a,3i8)') 'ilr, ii, mlr(ilr)%indexInGlobal(ii)', ilr, ii, mlr(ilr)%indexInGlobal(ii)
     end if
  end do
  if(ii/=mlr(ilr)%norbinlr) then
     write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ', iproc, ': ii/=mlr(ilr)%norbinlr', ii, mlr(ilr)%norbinlr
  end if
  !if(iproc==0) write(*,'(a,i6,200i5)') 'ilr, mlr(ilr)%indexInGlobal(ii)', ilr, mlr(ilr)%indexInGlobal(:)
end do

end subroutine determineLocalizationRegions








subroutine extractMatrix3(iproc, nproc, norb, norbp, orbstot, onWhichAtomPhi, onWhichMPI, nmat, ham, matmin, hamextract)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nmat, norb, norbp
type(orbitals_data),intent(in):: orbstot
integer,dimension(norb),intent(in):: onWhichAtomPhi, onWhichMPI
real(8),dimension(orbstot%norb,orbstot%norb,nmat),intent(in):: ham
type(matrixMinimization),intent(inout):: matmin
real(8),dimension(:,:,:),pointer,intent(out):: hamextract

! Local variables
integer:: jorb, jlr, jproc, jlrold, jjorb, ind, indlarge, jnd, jndlarge, ii, istat, jjlr
character(len=*),parameter:: subname='extractMatrix'

allocate(matmin%inWhichLocregExtracted(norbp), stat=istat)
call memocc(istat, matmin%inWhichLocregExtracted, 'matmin%inWhichLocregExtracted', subname)
allocate(matmin%inWhichLocregOnMPI(norbp), stat=istat)
call memocc(istat, matmin%inWhichLocregOnMPI, 'matmin%inWhichLocregOnMPI', subname)

! Allocate the matrices holding the extracted quantities. In principle this matrix
! has a different size for each localization region. To simplify the program, allocate them
! with the same size for all localization regions on a given MPI process.
matmin%norbmax=0
jlrold=-1
matmin%nlrp=0 ! localization regions per process
jjorb=0
do jorb=1,norb
  jlr=onWhichAtomPhi(jorb)
  jproc=onWhichMPI(jorb)
  if(iproc==jproc) then
     jjorb=jjorb+1
     if(jlr/=jlrold) then
        matmin%nlrp=matmin%nlrp+1
     end if
     if(matmin%mlr(jlr)%norbinlr>matmin%norbmax) then
        matmin%norbmax=matmin%mlr(jlr)%norbinlr
     end if
     matmin%inWhichLocregExtracted(jjorb)=jlr
     matmin%inWhichLocregOnMPI(jjorb)=matmin%nlrp
     jlrold=jlr
  end if
end do

allocate(matmin%indexInLocreg(matmin%nlrp), stat=istat)
call memocc(istat, matmin%indexInLocreg, 'matmin%indexInLocreg', subname)

! Allocate the matrix
allocate(hamextract(matmin%norbmax,matmin%norbmax,matmin%nlrp), stat=istat)
call memocc(istat, hamextract, 'hamextract', subname)
hamextract=0.d0

! Exctract the data from the large Hamiltonian.
jlrold=-1
jjlr=0
do jorb=1,norb
  jlr=onWhichAtomPhi(jorb)
  jproc=onWhichMPI(jorb)
  if(iproc==jproc) then
     if(jlr/=jlrold) then
        jjlr=jjlr+1
        matmin%indexInLocreg(jjlr)=jlr
        ! To make it work for both input guess (where we have nmat>1 different matrices) and
        ! for the iterative diagonalization (where we have only nmat=1 matrix).
        ii=min(jlr,nmat)
        do ind=1,matmin%mlr(jlr)%norbinlr
           indlarge=matmin%mlr(jlr)%indexInGlobal(ind)
           do jnd=1,matmin%mlr(jlr)%norbinlr
              jndlarge=matmin%mlr(jlr)%indexInGlobal(jnd)
              hamextract(jnd,ind,jjlr)=ham(jndlarge,indlarge,jjlr)
              !!write(30000+iproc,'(6i8,es20.10)') jnd,ind,jjlr,jndlarge,indlarge,jjlr,hamextract(jnd,ind,jjlr)
           end do
        end do
        jlrold=jlr
     end if
  end if
end do

if(jjlr/=nmat) then
  write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': jjlr/nmat',jjlr,nmat
  stop
end if


end subroutine extractMatrix3





subroutine vectorGlobalToLocal(norbtot, mlr, vglobal, vlocal)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: norbtot
type(matrixLocalizationRegion),intent(in):: mlr
real(8),dimension(norbtot),intent(in):: vglobal
real(8),dimension(mlr%norbinlr),intent(out):: vlocal

! Local variables
integer:: ilocal, iglobal

do ilocal=1,mlr%norbinlr
  iglobal=mlr%indexInGlobal(ilocal)
  vlocal(ilocal)=vglobal(iglobal)
end do


end subroutine vectorGlobalToLocal




subroutine vectorLocalToGlobal(norbtot, mlr, vlocal, vglobal)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: norbtot
type(matrixLocalizationRegion),intent(in):: mlr
real(8),dimension(mlr%norbinlr),intent(in):: vlocal
real(8),dimension(norbtot),intent(out):: vglobal

! Local variables
integer:: ilocal, iglobal

vglobal=0.d0
do ilocal=1,mlr%norbinlr
  iglobal=mlr%indexInGlobal(ilocal)
  vglobal(iglobal)=vlocal(ilocal)
end do


end subroutine vectorLocalToGlobal





subroutine determineOverlapRegionMatrix(iproc, nproc, lzd, mlr, orbs, orbstot, onWhichAtom, onWhichAtomPhi, comom)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: orbs, orbstot
integer,dimension(orbstot%norb),intent(in):: onWhichAtom
integer,dimension(orbs%norb),intent(in):: onWhichAtomPhi
type(matrixLocalizationRegion),dimension(lzd%nlr),intent(in):: mlr
type(p2pCommsOrthonormalityMatrix),intent(out):: comom

! Local variables
integer:: ilr, jlr, klr, novrlp, korb, istat, jlrold, jjlr, jjorb, jorb, kkorb, lorb, iorb, jorbout, iiorb, iorbout
integer:: is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3, ks1, ke1, ks2, ke2, ks3, ke3, ilrold
logical:: ovrlpx_ki, ovrlpy_ki, ovrlpz_ki, ovrlpx_kj, ovrlpy_kj, ovrlpz_kj, ovrlpx, ovrlpy, ovrlpz
logical:: overlapFound
character(len=*),parameter:: subname='determineOverlapRegionMatrix'


!allocate(comom%noverlap(lzd%nlr), stat=istat)
allocate(comom%noverlap(orbs%norb), stat=istat)
call memocc(istat, comom%noverlap, 'comom%noverlap', subname)

! First count the number of overlapping localization regions for each localization region.
!do ilr=1,lzd%nlr
do iorbout=1,orbs%norb
  ilr=orbs%inwhichlocreg(iorbout)
  call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
  novrlp=0
  do jorbout=1,orbs%norb
     jlr=onWhichAtomPhi(jorbout)
!!!! THIS IS THE ORIGINAL
     !!call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
     !!ovrlpx = ( is1<=je1 .and. ie1>=js1 )
     !!ovrlpy = ( is2<=je2 .and. ie2>=js2 )
     !!ovrlpz = ( is3<=je3 .and. ie3>=js3 )
     !!if(ovrlpx .and. ovrlpy .and. ovrlpz) then
     !!    novrlp=novrlp+1
     !!end if
     !! THIS IS NEW ############################
     ! Check whether there is a common element.
     outloop1: do iorb=1,mlr(ilr)%norbinlr
        iiorb=mlr(ilr)%indexInGlobal(iorb)
        do jorb=1,mlr(jlr)%norbinlr
           jjorb=mlr(jlr)%indexInGlobal(jorb)
           if(iiorb==jjorb) then
              novrlp=novrlp+1
              exit outloop1
           end if
        end do
     end do outloop1
  end do
  !comom%noverlap(ilr)=novrlp
  comom%noverlap(iorbout)=novrlp
  !!if(iproc==0) write(*,*) 'ilr, comom%noverlap(ilr)', ilr, comom%noverlap(ilr) 
end do


!allocate(comom%overlaps(maxval(comom%noverlap(:)),lzd%nlr), stat=istat)
allocate(comom%overlaps(maxval(comom%noverlap(:)),orbs%norb), stat=istat)
call memocc(istat, comom%overlaps, 'comom%overlaps', subname)
!do ilr=1,lzd%nlr
do iorbout=1,orbs%norb
  ilr=orbs%inwhichlocreg(iorbout)
  !comom%overlaps(:,ilr)=0
  comom%overlaps(:,iorbout)=0
  call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
  novrlp=0
  do jorbout=1,orbs%norb
     jlr=onWhichAtomPhi(jorbout)
!!!! THIS IS THE ORIGINAL
     !!call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
     !!ovrlpx = ( is1<=je1 .and. ie1>=js1 )
     !!ovrlpy = ( is2<=je2 .and. ie2>=js2 )
     !!ovrlpz = ( is3<=je3 .and. ie3>=js3 )
     !!if(ovrlpx .and. ovrlpy .and. ovrlpz) then
     !!    novrlp=novrlp+1
     !!    comom%overlaps(novrlp,ilr)=jorbout
     !!    if(iproc==0) write(*,'(2(a,i0))') 'locreg ',ilr,' overlaps with orbital ',jorbout
     !!end if
     ! THIS IS NEW ############################
     ! Check whether there is a common element.
     outloop2: do iorb=1,mlr(ilr)%norbinlr
        iiorb=mlr(ilr)%indexInGlobal(iorb)
        do jorb=1,mlr(jlr)%norbinlr
           jjorb=mlr(jlr)%indexInGlobal(jorb)
           if(iiorb==jjorb) then
              novrlp=novrlp+1
              !comom%overlaps(novrlp,ilr)=jorbout
              comom%overlaps(novrlp,iorbout)=jorbout
              !if(iproc==0) write(*,'(2(a,i0))') 'locreg ',ilr,' overlaps with orbital ',jorbout
              exit outloop2
           end if
        end do
     end do outloop2
  end do
  !!if(iproc==0) write(*,'(a,i4,3x,100i5)') 'ilr, comom%overlaps(,ilr)', ilr, comom%overlaps(:,ilr) 
end do

allocate(comom%olr(maxval(comom%noverlap(:)),lzd%nlr), stat=istat)
do ilr=1,lzd%nlr
  do iorb=1,maxval(comom%noverlap(:))
     call nullify_matrixLocalizationRegion(comom%olr(iorb,ilr))
  end do
end do



! Now determine which orbitals (corresponding to basis functions) will be in the overlap localization region.
!do ilr=1,lzd%nlr
ilrold=-1
do iorbout=1,orbs%norb
  ilr=orbs%inwhichlocreg(iorbout)
  if(ilr==ilrold) cycle
  ilrold=ilr
  call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
  comom%olr(:,ilr)%norbinlr=0
  !do jorbout=1,comom%noverlap(ilr)
  do jorbout=1,comom%noverlap(iorbout)
     !jjorb=comom%overlaps(jorbout,ilr)
     jjorb=comom%overlaps(jorbout,iorbout)
     jlr=onWhichAtomPhi(jjorb)
!!!! THIS IS THE ORIGINAL
     !!call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
     !!do korb=1,mlr(jlr)%norbinlr
     !!    lorb=mlr(jlr)%indexInGlobal(korb)
     !!    klr=onWhichAtom(lorb)
     !!    call getIndices(lzd%llr(klr), ks1, ke1, ks2, ke2, ks3, ke3)
     !!    ovrlpx_ki = ( ks1<=ie1 .and. ke1>=is1 )
     !!    ovrlpy_ki = ( ks2<=ie2 .and. ke2>=is2 )
     !!    ovrlpz_ki = ( ks3<=ie3 .and. ke3>=is3 )
     !!    ovrlpx_kj = ( ks1<=je1 .and. ke1>=js1 )
     !!    ovrlpy_kj = ( ks2<=je2 .and. ke2>=js2 )
     !!    ovrlpz_kj = ( ks3<=je3 .and. ke3>=js3 )
     !!    ovrlpx = ( ovrlpx_ki .and. ovrlpx_kj )
     !!    ovrlpy = ( ovrlpy_ki .and. ovrlpy_kj )
     !!    ovrlpz = ( ovrlpz_ki .and. ovrlpz_kj )
     !!    if(ovrlpx .and. ovrlpy .and. ovrlpz) then
     !!        comom%olr(jorbout,ilr)%norbinlr=comom%olr(jorbout,ilr)%norbinlr+1
     !!    end if
     !!end do
     ! THIS IS NEW ############################
     ! Check whether there is a common element.
     do iorb=1,mlr(ilr)%norbinlr
        iiorb=mlr(ilr)%indexInGlobal(iorb)
        do jorb=1,mlr(jlr)%norbinlr
           jjorb=mlr(jlr)%indexInGlobal(jorb)
           if(iiorb==jjorb) then
              novrlp=novrlp+1
              comom%olr(jorbout,ilr)%norbinlr=comom%olr(jorbout,ilr)%norbinlr+1
              !exit
           end if
        end do
     end do
     allocate(comom%olr(jorbout,ilr)%indexInGlobal(comom%olr(jorbout,ilr)%norbinlr), stat=istat)
     call memocc(istat, comom%olr(jorbout,ilr)%indexInGlobal, 'comom%olr(jorbout,ilr)%indexInGlobal', subname)
  end do
end do

!!do ilr=1,lzd%nlr
!!    do jorb=1,comom%noverlap(ilr)
!!        if(iproc==0) write(*,'(a,2i5,2i8)') 'ilr, jjorb, comom%overlaps(jorb,ilr), comom%olr(jorb,ilr)%norbinlr', ilr, jorb, comom%overlaps(jorb,ilr), comom%olr(jorb,ilr)%norbinlr
!!    end do
!!end do


! Determine the indices to switch from global region to localization region.
!do ilr=1,lzd%nlr
ilrold=-1
do iorbout=1,orbs%norb
  ilr=orbs%inwhichlocreg(iorbout)
  if(ilr==ilrold) cycle
  ilrold=ilr
  call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
  !do jorbout=1,comom%noverlap(ilr)
  do jorbout=1,comom%noverlap(iorbout)
     !jjorb=comom%overlaps(jorbout,ilr)
     jjorb=comom%overlaps(jorbout,iorbout)
     jlr=onWhichAtomPhi(jjorb)
!!!! THIS IS THE ORIGINAL
     !!call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
     !!kkorb=0
     !!comom%olr(jorbout,ilr)%indexInGlobal(:)=0
     !!do korb=1,mlr(jlr)%norbinlr
     !!    lorb=mlr(jlr)%indexInGlobal(korb)
     !!    klr=onWhichAtom(lorb)
     !!    call getIndices(lzd%llr(klr), ks1, ke1, ks2, ke2, ks3, ke3)
     !!    ovrlpx_ki = ( ks1<=ie1 .and. ke1>=is1 )
     !!    ovrlpy_ki = ( ks2<=ie2 .and. ke2>=is2 )
     !!    ovrlpz_ki = ( ks3<=ie3 .and. ke3>=is3 )
     !!    ovrlpx_kj = ( ks1<=je1 .and. ke1>=js1 )
     !!    ovrlpy_kj = ( ks2<=je2 .and. ke2>=js2 )
     !!    ovrlpz_kj = ( ks3<=je3 .and. ke3>=js3 )
     !!    ovrlpx = ( ovrlpx_ki .and. ovrlpx_kj )
     !!    ovrlpy = ( ovrlpy_ki .and. ovrlpy_kj )
     !!    ovrlpz = ( ovrlpz_ki .and. ovrlpz_kj )
     !!    if(ovrlpx .and. ovrlpy .and. ovrlpz) then
     !!        kkorb=kkorb+1
     !!        comom%olr(jorbout,ilr)%indexInGlobal(kkorb)=korb
     !!        !if(iproc==0) write(*,'(a,4i9)') 'orbitals in overlap region: ilr, jlr, kkorb, comom%olr(jorbout,ilr)%indexInGlobal(kkorb)', ilr, jlr, kkorb, comom%olr(jorbout,ilr)%indexInGlobal(kkorb)
     !!    end if
     !!end do
     ! THIS IS NEW ############################
     ! Check whether there is a common element.
     kkorb=0
     do iorb=1,mlr(ilr)%norbinlr
        iiorb=mlr(ilr)%indexInGlobal(iorb)
        do jorb=1,mlr(jlr)%norbinlr
           jjorb=mlr(jlr)%indexInGlobal(jorb)
           if(iiorb==jjorb) then
              kkorb=kkorb+1
              !comom%olr(jorbout,ilr)%indexInGlobal(kkorb)=iiorb
              !comom%olr(jorbout,ilr)%indexInGlobal(kkorb)=iorb
              comom%olr(jorbout,ilr)%indexInGlobal(kkorb)=jorb
              !if(iproc==0) write(*,'(a,4i9)') 'orbitals in overlap region: ilr, jlr, kkorb, comom%olr(jorbout,ilr)%indexInGlobal(kkorb)', ilr, jlr, kkorb, comom%olr(jorbout,ilr)%indexInGlobal(kkorb)
              !exit
           end if
        end do
     end do
  end do
end do

! With these indices it is possible to extract data from the global region to the
! overlap region. For example: comom%olr(jorb,ilr) allows to extract data from orbital
! jorb (the jorb-th orbital overlapping with region ilr) to the overlap region. To expand
! this overlap region to the whole region ilr, we need comom%(iorb,jlr), where jlr is the 
! localization region of jorb and iorb the iorb-th overlapping orbital of region jlr.
! This information is stored in comom%olrForExpansion:
! comom%olrForExpansion(1,jorb,ilr)=jlr
! comom%olrForExpansion(2,jorb,ilr)=iorb
allocate(comom%olrForExpansion(2,maxval(comom%noverlap(:)),lzd%nlr), stat=istat)
call memocc(istat, comom%olrForExpansion, 'comom%olrForExpansion', subname)
comom%olrForExpansion=55555
!do ilr=1,lzd%nlr
ilrold=-1
do iorbout=1,orbs%norb
  ilr=orbs%inwhichlocreg(iorbout)
  if(ilr==ilrold) cycle
  ilrold=ilr
  !do iorb=1,comom%noverlap(ilr)
  do iorb=1,comom%noverlap(iorbout)
     !jorb=comom%overlaps(iorb,ilr)
     jorb=comom%overlaps(iorb,iorbout)
     !jlr=onWhichAtom(jorb)
     jlr=onWhichAtomPhi(jorb)
     comom%olrForExpansion(1,iorb,ilr)=jlr
     !do korb=1,comom%noverlap(jlr)
     do korb=1,comom%noverlap(jorb)
        !kkorb=comom%overlaps(korb,jlr)
        kkorb=comom%overlaps(korb,jorb)
        !klr=onWhichAtom(kkorb)
        klr=onWhichAtomPhi(kkorb)
        !if(iproc==0) write(*,'(a,5i9)') 'ilr, iorb, jlr, korb, klr', ilr, iorb, jlr, korb, klr
        if(klr==ilr) then
           comom%olrForExpansion(2,iorb,ilr)=korb
        end if
     end do
     !if(iproc==0) write(*,'(a,4i8)') 'ilr, iorb, comom%olrForExpansion(1,iorb,ilr), comom%olrForExpansion(2,iorb,ilr)', ilr, iorb, comom%olrForExpansion(1,iorb,ilr), comom%olrForExpansion(2,iorb,ilr)
  end do
end do




end subroutine determineOverlapRegionMatrix





subroutine initCommsMatrixOrtho(iproc, nproc, norb, norb_par, isorb_par, onWhichAtomPhi, onWhichMPI, tag, comom)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, norb
integer,dimension(norb),intent(in):: onWhichAtomPhi, onWhichMPI
integer,dimension(0:nproc-1),intent(in):: norb_par, isorb_par
integer,intent(inout):: tag
type(p2pCommsOrthonormalityMatrix),intent(inout):: comom

! Local variables
integer:: jlrold,jproc,jj,jorb,jjorb,jlr,jjmax,istat,jkorb,mpisource,mpidest,istsource,istdest,ncount,korb,iall,kkorb
integer:: iorb, irecv, isend
integer,dimension(:),allocatable:: istsourcearr, istdestarr
character(len=*),parameter:: subname='initCommsMatrixOrtho'


allocate(istsourcearr(0:nproc-1), stat=istat)
call memocc(istat, istsourcearr, 'istsourcearr', subname)
istsourcearr=1
allocate(istdestarr(0:nproc-1), stat=istat)
call memocc(istat, istdestarr, 'istdestarr', subname)
istdestarr=1

comom%nrecvbuf=0
allocate(comom%indexInRecvBuf(maxval(comom%noverlap(:)),0:nproc-1), stat=istat)
call memocc(istat, comom%indexInRecvBuf, 'comom%indexInRecvBuf', subname)

allocate(comom%noverlapProc(0:nproc-1), stat=istat)
call memocc(istat, comom%noverlapProc, 'comom%noverlapProc', subname)

! Count how many orbitals each process will receive
do jproc=0,nproc-1
  jlrold=0
  comom%noverlapProc(jproc)=0
  do jorb=1,norb_par(jproc)
     jjorb=isorb_par(jproc)+jorb
     jlr=onWhichAtomPhi(jjorb)
     if(jlr==jlrold) cycle
     !do korb=1,comom%noverlap(jlr)
     do korb=1,comom%noverlap(jjorb)
        comom%noverlapProc(jproc)=comom%noverlapProc(jproc)+1
     end do
     jlrold=jlr
  end do
end do

allocate(comom%comarr(8,maxval(comom%noverlapProc(:)),0:nproc-1), stat=istat)
call memocc(istat, comom%comarr, 'comom%comarr', subname)
allocate(comom%overlapsProc(maxval(comom%noverlapProc(:)),0:nproc-1), stat=istat)
call memocc(istat, comom%overlapsProc, 'comom%overlapsProc', subname)
allocate(comom%communComplete(maxval(comom%noverlapProc(:)),0:nproc-1), stat=istat)
call memocc(istat, comom%communComplete, 'comom%communComplete', subname)

comom%nsendBuf=0
comom%nrecvBuf=0
do jproc=0,nproc-1
  jkorb=0
  jlrold=0
  do jorb=1,norb_par(jproc)
     jjorb=isorb_par(jproc)+jorb
     jlr=onWhichAtomPhi(jjorb)
     if(jlr==jlrold) cycle
     !do korb=1,comom%noverlap(jlr)
     do korb=1,comom%noverlap(jjorb)
        jkorb=jkorb+1
        !kkorb=comom%overlaps(korb,jlr)
        kkorb=comom%overlaps(korb,jjorb)
        mpidest=jproc
        mpisource=onWhichMPI(kkorb)
        istsource=istsourcearr(mpisource)
        ncount=comom%olr(korb,jlr)%norbinlr 
        istdest=istdestarr(mpidest)
        tag=tag+1
        comom%overlapsProc(jkorb,jproc)=kkorb
        call setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comom%comarr(1,jkorb,jproc))
        if(iproc==mpisource) then
           comom%nsendBuf=comom%nsendBuf+ncount
           !write(*,'(3(a,i0),6x,a,4i8)') 'process ',iproc,' adds ',ncount,' elements to comom%nsendBuf. Value after that: ',comom%nsendBuf,'. jjorb, jlr, korb, kkorb', jjorb, jlr, korb, kkorb
        end if
        if(iproc==mpidest) then
           !write(*,'(2(a,i0),a)') 'process ',iproc,' adds ',ncount,' elements to comom%nrecvBuf.'
           comom%nrecvbuf=comom%nrecvbuf+ncount
           comom%indexInRecvBuf(korb,jproc)=istdest
        end if
        istdestarr(mpidest)=istdestarr(mpidest)+ncount
        istsourcearr(mpisource)=istsourcearr(mpisource)+ncount
     end do
     jlrold=jlr
  end do
end do

allocate(comom%recvBuf(comom%nrecvbuf), stat=istat)
call memocc(istat, comom%recvBuf, 'comom%recvBuf', subname)
allocate(comom%sendBuf(comom%nsendbuf), stat=istat)
call memocc(istat, comom%sendBuf, 'comom%sendBuf', subname)

iall=-product(shape(istsourcearr))*kind(istsourcearr)
deallocate(istsourcearr, stat=istat)
call memocc(istat, iall, 'istsourcearr', subname)
iall=-product(shape(istdestarr))*kind(istdestarr)
deallocate(istdestarr, stat=istat)
call memocc(istat, iall, 'istdestarr', subname)


irecv=0
do jproc=0,nproc-1
  do iorb=1,comom%noverlapProc(jproc)
     mpidest=comom%comarr(4,iorb,jproc)
     ! The orbitals are on different processes, so we need a point to point communication.
     if(iproc==mpidest) then
        irecv=irecv+1
     end if
  end do
end do
! Number of receives per process, will be used later
comom%nrecv=irecv

isend=0
do jproc=0,nproc-1
  do iorb=1,comom%noverlapProc(jproc)
     mpisource=comom%comarr(1,iorb,jproc)
     ! The orbitals are on different processes, so we need a point to point communication.
     if(iproc==mpisource) then
        isend=isend+1
     end if
  end do
end do
! Number of sends per process, will be used later
comom%nsend=isend

allocate(comom%requests(max(comom%nrecv,comom%nsend),2), stat=istat)
call memocc(istat, comom%requests, 'comom%requests', subname)

end subroutine initCommsMatrixOrtho




subroutine orthonormalizeVectors(iproc, nproc, comm, nItOrtho, methTransformOverlap, blocksize_dsyev, blocksize_pdgemm, &
  orbs, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, isorb, nlr, newComm, mad, mlr, vec, comom)
use module_base
use module_types
use module_interfaces, exceptThisOne => orthonormalizeVectors
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, comm, nItOrtho, methTransformOverlap, blocksize_dsyev, blocksize_pdgemm
integer,intent(in):: norbmax, norbp, isorb, nlr, newComm
type(orbitals_data),intent(in):: orbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom, onWhichMPI
integer,dimension(0:nproc-1),intent(in):: isorb_par
type(matrixDescriptors),intent(in):: mad
type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
real(8),dimension(norbmax,norbp),intent(inout):: vec
type(p2pCommsOrthonormalityMatrix),intent(inout):: comom

! Local variables
integer:: noverlaps, iorb, iiorb, ilr, istat, ilrold, jorb, iall, it, iorbmax, jorbmax, i
real(8):: tt, dnrm2, dev
real(8),dimension(:,:),allocatable:: vecOvrlp, ovrlp
character(len=*),parameter:: subname='orthonormalizeVectors'
real(8),dimension(orbs%norb):: vecglobal

noverlaps=0
ilrold=0
do iorb=1,norbp
  iiorb=isorb+iorb
  ilr=onWhichAtom(iiorb)
  if(ilr/=ilrold) then
     !noverlaps=noverlaps+comom%noverlap(ilr)
     noverlaps=noverlaps+comom%noverlap(iiorb)
  end if
  ilrold=ilr
end do
allocate(vecOvrlp(norbmax,noverlaps), stat=istat)
call memocc(istat, vecOvrlp, 'vecOvrlp', subname)
allocate(ovrlp(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, ovrlp, 'ovrlp', subname)

do it=1,nItOrtho

  call extractToOverlapregion(iproc, nproc, orbs%norb, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, vec, comom)
  !call postCommsVectorOrthonormalization(iproc, nproc, newComm, comom)
  !call gatherVectors(iproc, nproc, newComm, comom)
  call postCommsVectorOrthonormalizationNew(iproc, nproc, newComm, comom)
  call gatherVectorsNew(iproc, nproc, comom)

  call expandFromOverlapregion(iproc, nproc, isorb, norbp, orbs, onWhichAtom, comom, norbmax, noverlaps, vecOvrlp)

  ! Calculate the overlap matrix.
  call calculateOverlap(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, orbs%norb, comom, mlr, onWhichAtom, &
       vec, vecOvrlp, newComm, ovrlp)
  dev=0.d0
  iorbmax=0
  jorbmax=0
  do iorb=1,orbs%norb
     do jorb=1,orbs%norb
        !if(iproc==0) write(300,'(2i8,es15.7)') iorb, jorb, ovrlp(jorb,iorb)
        if(iorb==jorb) then
           tt=abs(1.d0-ovrlp(jorb,iorb))
        else
           tt=abs(ovrlp(jorb,iorb))
        end if
        if(tt>dev) then
           dev=tt
           iorbmax=iorb
           jorbmax=jorb
        end if
     end do
  end do
  if(iproc==0) then
     write(*,'(a,es14.6,2(2x,i0))') 'max deviation from unity, position:',dev, iorbmax, jorbmax
  end if
  call overlapPowerMinusOneHalf(iproc, nproc, comm, methTransformOverlap, blocksize_dsyev, blocksize_pdgemm, orbs%norb, mad, ovrlp)



  call orthonormalLinearCombinations(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, &
       orbs%norb, comom, mlr, onWhichAtom, vecOvrlp, ovrlp, vec)

  ! Normalize the vectors
  do iorb=1,norbp
     tt=dnrm2(norbmax, vec(1,iorb), 1)
     call dscal(norbmax, 1/tt, vec(1,iorb), 1)
  end do

end do

iall=-product(shape(ovrlp))*kind(ovrlp)
deallocate(ovrlp, stat=istat)
call memocc(istat, iall, 'ovrlp', subname)
iall=-product(shape(vecOvrlp))*kind(vecOvrlp)
deallocate(vecOvrlp, stat=istat)
call memocc(istat, iall, 'vecOvrlp', subname)

end subroutine orthonormalizeVectors




subroutine orthoconstraintVectors(iproc, nproc, methTransformOverlap, correctionOrthoconstraint, blocksize_pdgemm, orbs, &
  onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, isorb, nlr, newComm, mlr, mad, vec, grad, comom, trace)
use module_base
use module_types
use module_interfaces, exceptThisOne => orthoconstraintVectors
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, methTransformOverlap, correctionOrthoconstraint, blocksize_pdgemm, norbmax
integer,intent(in):: norbp, isorb, nlr, newComm
type(orbitals_data),intent(in):: orbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom, onWhichMPI
integer,dimension(0:nproc-1),intent(in):: isorb_par
type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
type(matrixDescriptors),intent(in):: mad
real(8),dimension(norbmax,norbp),intent(inout):: vec, grad
type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
real(8),intent(out):: trace

! Local variables
integer:: noverlaps, iorb, iiorb, ilr, istat, ilrold, jorb, iall
real(8),dimension(:,:),allocatable:: gradOvrlp, vecOvrlp, lagmat, ovrlp
character(len=*),parameter:: subname='orthoconstraintVectors'
real(8):: ddot

noverlaps=0
ilrold=0
do iorb=1,norbp
  iiorb=isorb+iorb
  ilr=onWhichAtom(iiorb)
  if(ilr/=ilrold) then
     !noverlaps=noverlaps+comom%noverlap(ilr)
     noverlaps=noverlaps+comom%noverlap(iiorb)
  end if
  ilrold=ilr
end do
allocate(gradOvrlp(norbmax,noverlaps), stat=istat)
call memocc(istat, gradOvrlp, 'gradOvrlp', subname)
allocate(vecOvrlp(norbmax,noverlaps), stat=istat)
call memocc(istat, vecOvrlp, 'vecOvrlp', subname)
allocate(lagmat(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, lagmat, 'lagmat', subname)
allocate(ovrlp(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, ovrlp, 'ovrlp', subname)

call extractToOverlapregion(iproc, nproc, orbs%norb, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, grad, comom)

!call postCommsVectorOrthonormalization(iproc, nproc, newComm, comom)
!call gatherVectors(iproc, nproc, newComm, comom)
call postCommsVectorOrthonormalizationNew(iproc, nproc, newComm, comom)
call gatherVectorsNew(iproc, nproc, comom)

call expandFromOverlapregion(iproc, nproc, isorb, norbp, orbs, onWhichAtom, comom, norbmax, noverlaps, gradOvrlp)

! Calculate the Lagrange multiplier matrix <vec|grad>.
call calculateOverlap(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, orbs%norb, comom, mlr, onWhichAtom, vec,&
    gradOvrlp, newComm, lagmat)
!subroutine calculateOverlap(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb, comom, mlr, onWhichAtom, vec,&
!           vecOvrlp, newComm, ovrlp)
trace=0.d0
do iorb=1,orbs%norb
  trace=trace+lagmat(iorb,iorb)
end do

! Now we also have to calculate the overlap matrix.
call extractToOverlapregion(iproc, nproc, orbs%norb, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, vec, comom)
!call postCommsVectorOrthonormalization(iproc, nproc, newComm, comom)
!call gatherVectors(iproc, nproc, newComm, comom)
call postCommsVectorOrthonormalizationNew(iproc, nproc, newComm, comom)
call gatherVectorsNew(iproc, nproc, comom)
call expandFromOverlapregion(iproc, nproc, isorb, norbp, orbs, onWhichAtom, comom, norbmax, noverlaps, vecOvrlp)
call calculateOverlap(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, orbs%norb, comom, mlr, onWhichAtom, vec,&
    vecOvrlp, newComm, ovrlp)

! Now apply the orthoconstraint.
call applyOrthoconstraintVectors(iproc, nproc, methTransformOverlap, correctionOrthoconstraint, &
    blocksize_pdgemm, newComm, orbs%norb, &
    norbmax, norbp, isorb, nlr, noverlaps, onWhichAtom, vecOvrlp, ovrlp, lagmat, comom, mlr, mad, orbs, grad)

!call transformOverlapMatrix(iproc, nproc, orbs%norb, ovrlp)
!call orthonormalLinearCombinations(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, orbs%norb, comom, mlr, onWhichAtom, vecOvrlp, ovrlp, vec)

iall=-product(shape(lagmat))*kind(lagmat)
deallocate(lagmat, stat=istat)
call memocc(istat, iall, 'lagmat', subname)

iall=-product(shape(gradOvrlp))*kind(gradOvrlp)
deallocate(gradOvrlp, stat=istat)
call memocc(istat, iall, 'gradOvrlp', subname)

iall=-product(shape(ovrlp))*kind(ovrlp)
deallocate(ovrlp, stat=istat)
call memocc(istat, iall, 'ovrlp', subname)

iall=-product(shape(vecOvrlp))*kind(vecOvrlp)
deallocate(vecOvrlp, stat=istat)
call memocc(istat, iall, 'vecOvrlp', subname)

end subroutine orthoconstraintVectors





subroutine postCommsVectorOrthonormalization(iproc, nproc, newComm, comom)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, newComm
type(p2pCommsOrthonormalityMatrix),intent(inout):: comom

! Local variables
integer:: nsends, nreceives, jproc, iorb, mpisource, istsource, ncount, mpidest, istdest, tag, ierr

nsends=0
nreceives=0
comom%communComplete=.false.
do jproc=0,nproc-1
  do iorb=1,comom%noverlapProc(jproc)
     mpisource=comom%comarr(1,iorb,jproc)
     istsource=comom%comarr(2,iorb,jproc)
     ncount=comom%comarr(3,iorb,jproc)
     mpidest=comom%comarr(4,iorb,jproc)
     istdest=comom%comarr(5,iorb,jproc)
     tag=comom%comarr(6,iorb,jproc)
     !if(iproc==0) write(*,'(a,4i9)') 'jproc, iorb, mpisource, mpidest', jproc, iorb, mpisource, mpidest
     if(mpisource/=mpidest) then
        ! The orbitals are on different processes, so we need a point to point communication.
        if(iproc==mpisource) then
           !write(*,'(6(a,i0))') 'process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
           call mpi_isend(comom%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, newComm,&
                comom%comarr(7,iorb,jproc), ierr)
           !call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, lin%comsr%comarr(8,iorb,jproc), ierr)
           comom%comarr(8,iorb,jproc)=mpi_request_null !is this correct?
           nsends=nsends+1
        else if(iproc==mpidest) then
           !write(*,'(6(a,i0))') 'Aprocess ', mpidest, ' receives ', ncount,&
           !     ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
           call mpi_irecv(comom%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, newComm,&
                comom%comarr(8,iorb,jproc), ierr)
           comom%comarr(7,iorb,jproc)=mpi_request_null !is this correct?
           nreceives=nreceives+1
        else
           comom%comarr(7,iorb,jproc)=mpi_request_null
           comom%comarr(8,iorb,jproc)=mpi_request_null
        end if
     else
        ! The orbitals are on the same process, so simply copy them.
        if(iproc==mpisource) then
           call dcopy(ncount, comom%sendBuf(istsource), 1, comom%recvBuf(istdest), 1)
           !write(*,'(6(a,i0))') 'process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag
           comom%comarr(7,iorb,jproc)=mpi_request_null
           comom%comarr(8,iorb,jproc)=mpi_request_null
           nsends=nsends+1
           nreceives=nreceives+1
           comom%communComplete(iorb,iproc)=.true.
        else
           comom%comarr(7,iorb,jproc)=mpi_request_null
           comom%comarr(8,iorb,jproc)=mpi_request_null
        end if

     end if
  end do
end do



end subroutine postCommsVectorOrthonormalization





subroutine postCommsVectorOrthonormalizationNew(iproc, nproc, newComm, comom)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, newComm
type(p2pCommsOrthonormalityMatrix),intent(inout):: comom

! Local variables
integer:: isend, irecv, jproc, iorb, mpisource, istsource, ncount, mpidest, istdest, tag, ierr

irecv=0
comom%communComplete=.false.
do jproc=0,nproc-1
  do iorb=1,comom%noverlapProc(jproc)
     mpisource=comom%comarr(1,iorb,jproc)
     istsource=comom%comarr(2,iorb,jproc)
     ncount=comom%comarr(3,iorb,jproc)
     mpidest=comom%comarr(4,iorb,jproc)
     istdest=comom%comarr(5,iorb,jproc)
     tag=comom%comarr(6,iorb,jproc)
     ! The orbitals are on different processes, so we need a point to point communication.
     if(iproc==mpidest .and. nproc > 1) then
        !write(*,'(6(a,i0))') 'Bprocess ', mpidest, ' receives ', ncount,&
        !     ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
        tag=iorb
        irecv=irecv+1
        call mpi_irecv(comom%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, newComm,&
             comom%requests(irecv,2), ierr)
     end if
  end do
end do

! Number of receives per process, will be used later
comom%nrecv=irecv

isend=0
do jproc=0,nproc-1
  do iorb=1,comom%noverlapProc(jproc)
     mpisource=comom%comarr(1,iorb,jproc)
     istsource=comom%comarr(2,iorb,jproc)
     ncount=comom%comarr(3,iorb,jproc)
     mpidest=comom%comarr(4,iorb,jproc)
     istdest=comom%comarr(5,iorb,jproc)
     tag=comom%comarr(6,iorb,jproc)
     ! The orbitals are on different processes, so we need a point to point communication.
     if(iproc==mpisource .and. nproc > 1) then
        !write(*,'(6(a,i0))') 'process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
        tag=iorb
        isend=isend+1
        call mpi_isend(comom%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, newComm,&
             comom%requests(isend,1), ierr)
     else if (nproc ==1) then
        call vcopy(ncount,comom%sendBuf(istsource),1,comom%recvBuf(istdest),1)
     end if
  end do
end do
! Number of sends per process, will be used later
comom%nsend=isend


end subroutine postCommsVectorOrthonormalizationNew






subroutine gatherVectors(iproc, nproc, newComm, comom)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, newComm
type(p2pCommsOrthonormalityMatrix),intent(inout):: comom

! Local variables
integer:: jorb, mpisource, mpidest, nfast, nslow, nsameproc, ierr, jproc
integer,dimension(mpi_status_size):: stat
logical:: sendComplete, receiveComplete

!!! Check whether the communications have completed.
!!nfast=0
!!nsameproc=0
!!testLoop: do
!!    do jproc=0,nproc-1
!!        do jorb=1,comom%noverlapProc(jproc)
!!            if(comom%communComplete(jorb,jproc)) cycle
!!            call mpi_test(comom%comarr(7,jorb,jproc), sendComplete, stat, ierr)
!!            call mpi_test(comom%comarr(8,jorb,jproc), receiveComplete, stat, ierr)
!!            if(sendComplete .and. receiveComplete) comom%communComplete(jorb,jproc)=.true.
!!            if(comom%communComplete(jorb,jproc)) then
!!                !write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', jorb
!!                mpisource=comom%comarr(1,jorb,jproc)
!!                mpidest=comom%comarr(4,jorb,jproc)
!!                if(mpisource/=mpidest) then
!!                    nfast=nfast+1
!!                else
!!                    nsameproc=nsameproc+1
!!                end if
!!            end if
!!        end do
!!    end do
!!    ! If we made it until here, either all all the communication is
!!    ! complete or we better wait for each single orbital.
!!    exit testLoop
!!end do testLoop


! Wait for the communications that have not completed yet
nslow=0
do jproc=0,nproc-1
  do jorb=1,comom%noverlapProc(jproc)
     if(comom%communComplete(jorb,jproc)) cycle
     !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
     nslow=nslow+1
     call mpi_wait(comom%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
     call mpi_wait(comom%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
     comom%communComplete(jorb,jproc)=.true.
  end do
end do

!call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
!call mpiallred(nfast, 1, mpi_sum, newComm, ierr)
!call mpiallred(nslow, 1, mpi_sum, newComm, ierr)
!call mpiallred(nsameproc, 1, mpi_sum, newComm, ierr)
!if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
!                       nfast, ' could be overlapped with computation.'
!if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'


end subroutine gatherVectors






subroutine gatherVectorsNew(iproc, nproc, comom)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(p2pCommsOrthonormalityMatrix),intent(inout):: comom

! Local variables
integer:: i, nsend, nrecv, ind, ierr

if (nproc >1) then
  nsend=0
  waitLoopSend: do
     call mpi_waitany(comom%nsend-nsend, comom%requests(1,1), ind, mpi_status_ignore, ierr)
     nsend=nsend+1
     do i=ind,comom%nsend-nsend
        comom%requests(i,1)=comom%requests(i+1,1)
     end do
     if(nsend==comom%nsend) exit waitLoopSend
  end do waitLoopSend


  nrecv=0
  waitLoopRecv: do
     call mpi_waitany(comom%nrecv-nrecv, comom%requests(1,2), ind, mpi_status_ignore, ierr)
     nrecv=nrecv+1
     do i=ind,comom%nrecv-nrecv
        comom%requests(i,2)=comom%requests(i+1,2)
     end do
     if(nrecv==comom%nrecv) exit waitLoopRecv
  end do waitLoopRecv
end if

end subroutine gatherVectorsNew








subroutine extractToOverlapregion(iproc, nproc, norb, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, vec, comom)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, norbmax, norb, norbp
integer,dimension(norb),intent(in):: onWhichAtom, onWhichMPI
integer,dimension(0:nproc-1),intent(in):: isorb_par
real(8),dimension(norbmax,norbp),intent(in):: vec
type(p2pCommsOrthonormalityMatrix),intent(inout):: comom

! Local variables
integer:: ilrold, iiprocold, ist, ilr, iiproc, jjorb, jorb, iorb, jjproc, korb, ind, i, jjlr

!write(*,*) 'iproc, norbmax, norbp', iproc, norbmax, norbp

ilrold=-1
iiprocold=-1
ist=0
do iorb=1,norb
    ilr=onWhichAtom(iorb)
    iiproc=onWhichMPI(iorb)
    if(ilr==ilrold .and.  iiproc==iiprocold) cycle ! otherwise we would extract the same again
    !do jorb=1,comom%noverlap(ilr)
    do jorb=1,comom%noverlap(iorb)
        !jjorb=comom%overlaps(jorb,ilr)
        jjorb=comom%overlaps(jorb,iorb)
        jjlr=onWhichAtom(jjorb)
        jjproc=onWhichMPI(jjorb)
        korb=jjorb-isorb_par(jjproc)
        !if(iproc==0) write(*,'(a,8i8)') 'iorb, ilr, iiproc, jorb, jjorb, jjlr, jjproc, korb', iorb, ilr, iiproc, jorb, jjorb, jjlr, jjproc, korb
        if(iproc==jjproc) then
            !write(*,'(a,6i9)') 'iproc, jorb, jjorb, ilr, comom%overlaps(jorb,ilr), comom%olr(jorb,ilr)%norbinlr', iproc, jorb, jjorb, ilr, comom%overlaps(jorb,ilr), comom%olr(jorb,ilr)%norbinlr
            !write(*,'(3(a,i0),6x,a,4i8)') 'process ',iproc,' adds ',comom%olr(jorb,ilr)%norbinlr,' elements to position ',ist,'. iorb, ilr, jorb, jjorb', iorb, ilr, jorb, jjorb
            do i=1,comom%olr(jorb,ilr)%norbinlr
                ind=comom%olr(jorb,ilr)%indexInGlobal(i)
                comom%sendBuf(ist+i)=vec(ind,korb)
            end do
            ist=ist+comom%olr(jorb,ilr)%norbinlr
        end if
    end do
    ilrold=ilr
    iiprocold=iiproc
end do

if(ist/=comom%nsendBuf) then
    write(*,'(1x,a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': ist/=comom%nsendBuf',ist,comom%nsendBuf
end if


end subroutine extractToOverlapregion





subroutine expandFromOverlapregion(iproc, nproc, isorb, norbp, orbs, onWhichAtom, comom, norbmax, noverlaps, vecOvrlp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, isorb, norbp, norbmax, noverlaps
type(orbitals_data),intent(in):: orbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(p2pCommsOrthonormalityMatrix),intent(in):: comom
real(8),dimension(norbmax,noverlaps),intent(out):: vecOvrlp

! Local variables
integer:: ilrold, ist, iorb, iiorb, ilr, jorb, klr, korb, i, ind, ijorb


vecOvrlp=0.d0
ilrold=0
ist=0
ijorb=0
do iorb=1,norbp
    iiorb=isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) cycle
    !do jorb=1,comom%noverlap(ilr)
    do jorb=1,comom%noverlap(iiorb)
        ijorb=ijorb+1
        klr=comom%olrForExpansion(1,jorb,ilr)
        korb=comom%olrForExpansion(2,jorb,ilr)
        !if(iproc==0) write(*,'(a,4i8)') 'iorb, jorb, comom%olrForExpansion(1,jorb,ilr), comom%olrForExpansion(2,jorb,ilr)', iorb, jorb, comom%olrForExpansion(1,jorb,ilr), comom%olrForExpansion(2,jorb,ilr)
        do i=1,comom%olr(korb,klr)%norbinlr
            ind=comom%olr(korb,klr)%indexInGlobal(i)
            vecOvrlp(ind,ijorb)=comom%recvBuf(ist+i)
            !if(iproc==4) write(*,'(a,9i9)') 'iorb, iiorb, ilr, jorb, klr, korb, ist, i, ind', iorb, iiorb, ilr, jorb, klr, korb, ist, i, ind
        end do
        ist=ist+comom%olr(korb,klr)%norbinlr
    end do
    ilrold=ilr
end do

if(ist/=comom%nrecvBuf) then
    write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': ist/=comom%nrecvBuf',ist,comom%nrecvBuf
end if


end subroutine expandFromOverlapregion



subroutine calculateOverlap(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb, comom, mlr, onWhichAtom, vec,&
           vecOvrlp, newComm, ovrlp)
use module_base
use module_types
implicit none

! Calling arguments 
integer,intent(in):: iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb, newComm
type(p2pCommsOrthonormalityMatrix),intent(in):: comom
type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
integer,dimension(norb),intent(in):: onWhichAtom
real(8),dimension(norbmax,norbp),intent(in):: vec
real(8),dimension(norbmax,noverlaps),intent(in):: vecOvrlp
real(8),dimension(norb,norb),intent(out):: ovrlp

! Local variables
integer:: ijorb, ilrold, ilr, iorb, iiorb, ncount, jjorb, jorb, ierr
real(8):: ddot



ovrlp=0.d0

ijorb=0
ilrold=0
do iorb=1,norbp
    iiorb=isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) then
         ! Put back index if we are in the same localization region, since then we can use the same vecOvrlp again.
         !ijorb=ijorb-comom%noverlap(ilr) 
         ijorb=ijorb-comom%noverlap(iiorb) 
    end if
    ncount=mlr(ilr)%norbinlr
    !do jorb=1,comom%noverlap(ilr)
    do jorb=1,comom%noverlap(iiorb)
        ijorb=ijorb+1
        !jjorb=comom%overlaps(jorb,ilr)
        jjorb=comom%overlaps(jorb,iiorb)
        ovrlp(iiorb,jjorb)=ddot(ncount, vec(1,iorb), 1, vecOvrlp(1,ijorb), 1)
    end do
    ilrold=ilr
end do

call mpiallred(ovrlp(1,1), norb**2, mpi_sum, newComm, ierr)


end subroutine calculateOverlap




subroutine orthonormalLinearCombinations(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb, comom, mlr, onWhichAtom,&
           vecOvrlp, ovrlp, vec)
use module_base
use module_types
implicit none

! Calling arguments 
integer,intent(in):: iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb
type(p2pCommsOrthonormalityMatrix),intent(in):: comom
type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
integer,dimension(norb),intent(in):: onWhichAtom
real(8),dimension(norbmax,noverlaps),intent(in):: vecOvrlp
real(8),dimension(norb,norb),intent(in):: ovrlp
real(8),dimension(norbmax,norbp),intent(inout):: vec

! Local variables
integer:: ijorb, ilrold, ilr, iorb, iiorb, ncount, jjorb, jorb, ierr, istat, iall
real(8):: ddot
real(8),dimension(:,:),allocatable:: vecTemp
character(len=*),parameter:: subname='orthonormalLinearCombinations'

allocate(vecTemp(norbmax,norbp), stat=istat)
call memocc(istat, vecTemp, 'vecTemp',subname)

call dcopy(norbmax*norbp, vec(1,1), 1, vecTemp(1,1), 1)

vec=0.d0

ijorb=0
ilrold=0
do iorb=1,norbp
    iiorb=isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) then
         ! Put back index if we are in the same localization region, since then we can use the same vecOvrlp again.
         !ijorb=ijorb-comom%noverlap(ilr) 
         ijorb=ijorb-comom%noverlap(iiorb) 
    end if
    ncount=mlr(ilr)%norbinlr
    !do jorb=1,comom%noverlap(ilr)
    do jorb=1,comom%noverlap(iiorb)
        ijorb=ijorb+1
        !jjorb=comom%overlaps(jorb,ilr)
        jjorb=comom%overlaps(jorb,iiorb)
        call daxpy(ncount, ovrlp(jjorb,iiorb), vecOvrlp(1,ijorb), 1, vec(1,iorb), 1)
    end do
    ilrold=ilr
end do

iall=-product(shape(vecTemp))*kind(vecTemp)
deallocate(vecTemp, stat=istat)
call memocc(istat, iall, 'vecTemp', subname)

end subroutine orthonormalLinearCombinations




subroutine applyOrthoconstraintVectors(iproc, nproc, methTransformOverlap, correctionOrthoconstraint, blocksize_pdgemm, &
           comm, norb, norbmax, norbp, isorb, nlr, noverlaps, onWhichAtom, vecOvrlp, ovrlp, &
           lagmat, comom, mlr, mad, orbs, grad)
use module_base
use module_types
use module_interfaces, exceptThisOne => applyOrthoconstraintVectors
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, methTransformOverlap, correctionOrthoconstraint, blocksize_pdgemm
integer,intent(in):: comm, norb, norbmax, norbp, isorb, nlr, noverlaps
integer,dimension(norb),intent(in):: onWhichAtom
real(8),dimension(norbmax,noverlaps),intent(in):: vecOvrlp
real(8),dimension(norb,norb),intent(in):: ovrlp
real(8),dimension(norb,norb),intent(inout):: lagmat
type(p2pCommsOrthonormalityMatrix),intent(in):: comom
type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
type(matrixDescriptors),intent(in):: mad
type(orbitals_data),intent(in):: orbs
real(8),dimension(norbmax,norbp),intent(inout):: grad

! Local variables
integer:: info, iorb, ilrold, iiorb, jjorb, ilr, ncount, jorb, ijorb, istat, iall
real(8),dimension(:,:),allocatable:: ovrlp2, ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans
real(8):: tt
character(len=*),parameter:: subname='applyOrthoconstraintVectors'


allocate(ovrlp_minus_one_lagmat(norb,norb), stat=istat)
call memocc(istat, ovrlp_minus_one_lagmat, 'ovrlp_minus_one_lagmat', subname)
allocate(ovrlp_minus_one_lagmat_trans(norb,norb), stat=istat)
call memocc(istat, ovrlp_minus_one_lagmat_trans, 'ovrlp_minus_one_lagmat_trans', subname)
allocate(ovrlp2(norb,norb), stat=istat)
call memocc(istat, ovrlp2, 'ovrlp2', subname)

call dcopy(norb**2, ovrlp(1,1), 1, ovrlp2(1,1), 1)
!!! Invert the overlap matrix
!!call dpotrf('l', norb, ovrlp2(1,1), norb, info)
!!if(info/=0) then
!!    write(*,'(1x,a,i0)') 'ERROR in dpotrf, info=',info
!!    stop
!!end if
!!call dpotri('l', norb, ovrlp2(1,1), norb, info)
!!if(info/=0) then
!!    write(*,'(1x,a,i0)') 'ERROR in dpotri, info=',info
!!    stop
!!end if

correctionIf: if(correctionOrthoconstraint==0) then
    ! Invert the overlap matrix
    call overlapPowerMinusOne(iproc, nproc, methTransformOverlap, norb, mad, orbs, ovrlp2)
    
    ! Multiply the Lagrange multiplier matrix with S^-1/2.
    ! First fill the upper triangle.
    do iorb=1,norb
        do jorb=1,iorb-1
            ovrlp2(jorb,iorb)=ovrlp2(iorb,jorb)
        end do
    end do
    if(blocksize_pdgemm<0) then
        !!call dgemm('n', 'n', norb, norb, norb, 1.d0, ovrlp2(1,1), norb, lagmat(1,1), norb, &
        !!     0.d0, ovrlp_minus_one_lagmat(1,1), norb)
        !!call dgemm('n', 't', norb, norb, norb, 1.d0, ovrlp2(1,1), norb, lagmat(1,1), norb, &
        !!     0.d0, ovrlp_minus_one_lagmat_trans(1,1), norb)
        !!call dsymm('l', 'l', norb, norb, 1.d0, ovrlp2(1,1), norb, lagmat(1,1), norb, &
        !!     0.d0, ovrlp_minus_one_lagmat(1,1), norb)
        ovrlp_minus_one_lagmat=0.d0
        !!call dgemm_compressed2(iproc, nproc, norb, mad%nsegline, mad%nseglinemax, mad%keygline, mad%nsegmatmul, mad%keygmatmul, &
        !!     ovrlp2, lagmat, ovrlp_minus_one_lagmat)
        call dgemm_compressed_parallel(iproc, nproc, norb, mad%nsegline, mad%nseglinemax,&
             mad%keygline, mad%nsegmatmul, mad%keygmatmul, &
             orbs%norb_par, orbs%isorb_par, orbs%norbp, ovrlp2, lagmat, ovrlp_minus_one_lagmat)
        !ovrlp_minus_one_lagmat=lagmat
        ! Transpose lagmat
        do iorb=1,norb
            do jorb=iorb+1,norb
                tt=lagmat(jorb,iorb)
                lagmat(jorb,iorb)=lagmat(iorb,jorb)
                lagmat(iorb,jorb)=tt
            end do
        end do
        !!call dsymm('l', 'l', norb, norb, 1.d0, ovrlp2(1,1), norb, lagmat(1,1), norb, &
        !!     0.d0, ovrlp_minus_one_lagmat_trans(1,1), norb)
        !ovrlp_minus_one_lagmat_trans=lagmat
        ovrlp_minus_one_lagmat_trans=0.d0
        !!call dgemm_compressed2(iproc, nproc, norb, mad%nsegline, mad%nseglinemax, mad%keygline, mad%nsegmatmul, mad%keygmatmul, &
        !!     ovrlp2, lagmat, ovrlp_minus_one_lagmat_trans)
        call dgemm_compressed_parallel(iproc, nproc, norb, mad%nsegline, mad%nseglinemax,&
             mad%keygline, mad%nsegmatmul, mad%keygmatmul, &
             orbs%norb_par, orbs%isorb_par, orbs%norbp, ovrlp2, lagmat, ovrlp_minus_one_lagmat_trans)
    
    else
        call dsymm_parallel(iproc, nproc, blocksize_pdgemm, comm, 'l', 'l', norb, norb, 1.d0, ovrlp2(1,1), &
             norb, lagmat(1,1), norb, &
             0.d0, ovrlp_minus_one_lagmat(1,1), norb)
        ! Transpose lagmat
        do iorb=1,norb
            do jorb=iorb+1,norb
                tt=lagmat(jorb,iorb)
                lagmat(jorb,iorb)=lagmat(iorb,jorb)
                lagmat(iorb,jorb)=tt
            end do
        end do
        call dsymm_parallel(iproc, nproc, blocksize_pdgemm, comm, 'l', 'l', norb, norb, 1.d0, ovrlp2(1,1), norb, &
            lagmat(1,1), norb, 0.d0, ovrlp_minus_one_lagmat_trans(1,1), norb)
    end if
else if(correctionOrthoconstraint==1) then correctionIf
    do iorb=1,norb
        do jorb=1,norb
            ovrlp_minus_one_lagmat(jorb,iorb)=lagmat(jorb,iorb)
            ovrlp_minus_one_lagmat_trans(jorb,iorb)=lagmat(iorb,jorb)
        end do
    end do
end if correctionIf


ilrold=-1
ijorb=0
do iorb=1,norbp
    iiorb=isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) then
        ! Set back the index of lphiovrlp, since we again need the same orbitals.
        !ijorb=ijorb-comom%noverlap(ilr)
        ijorb=ijorb-comom%noverlap(iiorb)
    end if
    ncount=mlr(ilr)%norbinlr
    !do jorb=1,comom%noverlap(ilr)
    do jorb=1,comom%noverlap(iiorb)
        ijorb=ijorb+1
        !jjorb=comom%overlaps(jorb,ilr)
        jjorb=comom%overlaps(jorb,iiorb)
        call daxpy(ncount, -.5d0*ovrlp_minus_one_lagmat(jjorb,iiorb), vecOvrlp(1,ijorb), 1, grad(1,iorb), 1)
        call daxpy(ncount, -.5d0*ovrlp_minus_one_lagmat_trans(jjorb,iiorb), vecOvrlp(1,ijorb), 1, grad(1,iorb), 1)
    end do
    ilrold=ilr
end do


iall=-product(shape(ovrlp_minus_one_lagmat))*kind(ovrlp_minus_one_lagmat)
deallocate(ovrlp_minus_one_lagmat, stat=istat)
call memocc(istat, iall, 'ovrlp_minus_one_lagmat', subname)
iall=-product(shape(ovrlp_minus_one_lagmat_trans))*kind(ovrlp_minus_one_lagmat_trans)
deallocate(ovrlp_minus_one_lagmat_trans, stat=istat)
call memocc(istat, iall, 'ovrlp_minus_one_lagmat_trans', subname)
iall=-product(shape(ovrlp2))*kind(ovrlp2)
deallocate(ovrlp2, stat=istat)
call memocc(istat, iall, 'ovrlp2', subname)


end subroutine applyOrthoconstraintVectors




subroutine buildLinearCombinations(iproc, nproc, lzdig, lzd, orbsig, orbs, input, coeff, lchi, locregShape, &
           tag, comonig, opig, madig, lphi)
use module_base
use module_types
use module_interfaces, exceptThisOne => buildLinearCombinations
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(local_zone_descriptors),intent(in):: lzdig, lzd
type(orbitals_data),intent(in):: orbsig, orbs
type(input_variables),intent(in):: input
real(8),dimension(orbsig%norb,orbs%norb),intent(in):: coeff
real(8),dimension(orbsig%npsidim_orbs),intent(in):: lchi
character(len=1),intent(in):: locregShape
integer,intent(inout):: tag
type(p2pComms):: comonig
type(overlapParameters):: opig
type(matrixDescriptors):: madig
real(8),dimension(orbs%npsidim_orbs),intent(out):: lphi

! Local variables
integer:: istat, iall, ist, jst, ilr, ilrold, iorb, iiorb, ncount, jorb, jjorb
!type(overlapParameters):: op
!type(p2pCommsOrthonormality):: comon
real(8),dimension(:),allocatable:: lchiovrlp
character(len=*),parameter:: subname='buildLinearCombinations'
!type(matrixDescriptors):: mad !just for calling collectnew, not really needed
real(8),dimension(:,:),allocatable:: ttmat
real(8):: tt1, tt2, tt3

!tag=10000
!call initCommsOrtho(iproc, nproc, lzdig, orbsig, orbsig%inWhichLocreg, input, locregShape, op, comon, tag)
allocate(lchiovrlp(opig%ndim_lphiovrlp), stat=istat)
call memocc(istat, lchiovrlp, 'lchiovrlp',subname)

call allocateCommuncationBuffersOrtho(comonig, subname)
!call extractOrbital2(iproc,nproc,orbsig,orbsig%npsidim,orbsig%inWhichLocreg,lzdig,op,lchi,comon)
call extractOrbital3(iproc,nproc,orbsig,orbsig%npsidim_orbs,orbsig%inWhichLocreg,&
     lzdig,opig,lchi,comonig%nsendBuf,comonig%sendBuf)
!call postCommsOverlap(iproc, nproc, comon)
call postCommsOverlapNew(iproc, nproc, orbsig, opig, lzdig, lchi, comonig, tt1, tt2)
!call gatherOrbitals2(iproc, nproc, comon)
!!allocate(ttmat(orbsig%norb,orbsig%norb))
call collectnew(iproc, nproc, comonig, madig, opig, orbsig, lzdig, comonig%nsendbuf, &
     comonig%sendbuf, comonig%nrecvbuf, comonig%recvbuf, tt1, tt2, tt3)
!!deallocate(ttmat)
call expandOrbital2(iproc, nproc, orbsig, input, orbsig%inWhichLocreg, lzdig, opig, comonig, lchiovrlp)
call deallocateCommuncationBuffersOrtho(comonig, subname)



lphi=0.d0

ist=1
jst=1
ilrold=-1
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=orbs%inWhichLocreg(iiorb)
    if(ilr==ilrold) then
        ! Set back the index of lphiovrlp, since we again need the same orbitals.
        jst=jst-opig%noverlaps(iiorb-1)*ncount
    end if
    !write(*,'(a,6i13)') 'iproc, iorb, iiorb, op%noverlaps(iiorb), ilr, jst', iproc, iorb, iiorb, op%noverlaps(iiorb), ilr, jst
    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    do jorb=1,opig%noverlaps(iiorb)
        jjorb=opig%overlaps(jorb,iiorb)
        !call daxpy(ncount, ovrlp(jjorb,iiorb), lphiovrlp(jst), 1, lphi(ist), 1)
        call daxpy(ncount, coeff(jjorb,iiorb), lchiovrlp(jst), 1, lphi(ist), 1)
        jst=jst+ncount
    end do

    ist=ist+ncount
    ilrold=ilr

end do

!!if(ist/=orbs%npsidim+1) then
!!    write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': ist/=orbs%npsidim+1',ist,orbs%npsidim+1
!!    stop
!!end if
if(ist>orbs%npsidim_orbs+1) then
    write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': ist/=orbs%npsidim_orbs+1',ist,orbs%npsidim_orbs+1
    stop
end if



!call deallocate_overlapParameters(op, subname)
!call deallocate_p2pCommsOrthonormality(comon, subname)


iall=-product(shape(lchiovrlp))*kind(lchiovrlp)
deallocate(lchiovrlp, stat=istat)
call memocc(istat, iall, 'lchiovrlp', subname)



end subroutine buildLinearCombinations


subroutine buildLinearCombinationsVariable(iproc, nproc, lzdig, lzd, orbsig, orbs, input, coeff, lchi, tag, lphi)
use module_base
use module_types
use module_interfaces, exceptThisOne => buildLinearCombinationsVariable
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(local_zone_descriptors),intent(in):: lzdig, lzd
type(orbitals_data),intent(in):: orbsig, orbs
type(input_variables),intent(in):: input
real(8),dimension(orbsig%norb,orbs%norb),intent(in):: coeff
real(8),dimension(orbsig%npsidim_orbs),intent(in):: lchi
integer,intent(inout):: tag
real(8),dimension(orbs%npsidim_orbs),intent(out):: lphi

! Local variables
integer:: istat, iall, ist, jst, ilr, ilrold, iorb, iiorb, ncount, jorb, jjorb, ii
type(overlapParameters):: op
type(p2pComms):: comon
real(8),dimension(:),allocatable:: lchiovrlp
character(len=*),parameter:: subname='buildLinearCombinationsVariable'

!tag=10000
call initCommsOrthoVariable(iproc, nproc, lzdig, orbs, orbsig, orbsig%inWhichLocreg, input, op, comon, tag)
allocate(lchiovrlp(op%ndim_lphiovrlp), stat=istat)
call memocc(istat, lchiovrlp, 'lchiovrlp',subname)

call allocateCommuncationBuffersOrtho(comon, subname)
call extractOrbital2Variable(iproc, nproc, orbs, orbsig, orbsig%npsidim_orbs, lzdig, op, lchi, comon)
call postCommsOverlap(iproc, nproc, comon)
call gatherOrbitals2(iproc, nproc, comon)
!call mpi_barrier(mpi_comm_world, ist)
!stop
call expandOrbital2Variable(iproc, nproc, orbs, input, lzdig, op, comon, lchiovrlp)

call deallocateCommuncationBuffersOrtho(comon, subname)



lphi=0.d0

ist=1
jst=1
ilrold=-1
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=orbs%inWhichLocreg(iiorb)
    if(ilr==ilrold) then
        ! Set back the index of lphiovrlp, since we again need the same orbitals.
        jst=jst-op%noverlaps(iiorb-1)*ncount
    end if
    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        !call daxpy(ncount, ovrlp(jjorb,iiorb), lphiovrlp(jst), 1, lphi(ist), 1)
        call daxpy(ncount, coeff(jjorb,iiorb), lchiovrlp(jst), 1, lphi(ist), 1)
        jst=jst+ncount
    end do

    ist=ist+ncount
    ilrold=ilr

end do

if(ist/=orbs%npsidim_orbs+1) then
    write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,&
         ': ist/=orbsig%npsidim+1',ist,orbsig%npsidim_orbs+1
    stop
end if



call deallocate_overlapParameters(op, subname)
call deallocate_p2pComms(comon, subname)


iall=-product(shape(lchiovrlp))*kind(lchiovrlp)
deallocate(lchiovrlp, stat=istat)
call memocc(istat, iall, 'lchiovrlp', subname)



end subroutine buildLinearCombinationsVariable




subroutine buildLinearCombinationsLocalized3(iproc, nproc, orbsig, orbs, comms, at, Glr, input, norbsPerType, &
           onWhichAtom, lchi, lphi, rxyz, onWhichAtomPhi, lin, lzdig, nlocregPerMPI, tag, ham3, comonig, opig, madig)
!
use module_base
use module_types
use module_interfaces, exceptThisOne => buildLinearCombinationsLocalized3
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nlocregPerMPI
type(orbitals_data),intent(in):: orbsig, orbs
type(communications_arrays),intent(in):: comms
type(atoms_data),intent(in):: at
type(locreg_descriptors),intent(in):: Glr
type(input_variables),intent(in):: input
type(linearParameters),intent(in):: lin
type(local_zone_descriptors),intent(inout):: lzdig
integer,dimension(at%ntypes):: norbsPerType
integer,dimension(orbsig%norb),intent(in):: onWhichAtom
real(8),dimension(max(orbsig%npsidim_orbs,orbsig%npsidim_comp)):: lchi
real(8),dimension(max(lin%orbs%npsidim_orbs,lin%orbs%npsidim_comp)):: lphi
real(8),dimension(3,at%nat):: rxyz
integer,dimension(orbs%norb):: onWhichAtomPhi
!!real(8),dimension(orbsig%norb,orbsig%norb,at%nat),intent(inout):: ham
integer,intent(inout):: tag
real(8),dimension(orbsig%norb,orbsig%norb,nlocregPerMPI),intent(inout):: ham3
type(p2pComms):: comonig
type(overlapParameters):: opig
type(matrixDescriptors):: madig

! Local variables
integer:: iorb, jorb, korb, iat, ist, jst, nvctrp, iall, istat, ierr, infoCoeff, k, l,it, iiAt, jjAt, methTransformOverlap, iiorb
real(8),dimension(:),allocatable:: alpha, coeffPad, coeff2, gradTemp, gradOld, fnrmArr, fnrmOvrlpArr, fnrmOldArr, grad
real(8),dimension(:,:),allocatable:: coeff, lagMat, lcoeff, lgrad, lgradold
integer,dimension(:),allocatable:: recvcounts, displs, norb_par
real(8):: ddot, cosangle, tt, dnrm2, fnrm, meanAlpha, cut, trace, traceOld, fnrmMax, valin, valout, tt2
logical:: converged
character(len=*),parameter:: subname='buildLinearCombinationsLocalized'
real(4):: ttreal, builtin_rand
integer:: wholeGroup, newGroup, newComm, norbtot, isx
integer,dimension(:),allocatable:: newID
integer:: ii, jproc, norbTarget, sendcount, ilr, iilr, ilrold, jlr
type(inguessParameters):: ip
real(8),dimension(:,:,:),pointer:: hamextract
type(p2pCommsOrthonormalityMatrix):: comom
type(matrixMinimization):: matmin
logical:: same
type(localizedDIISParameters):: ldiis
type(matrixDescriptors):: mad

  if(iproc==0) write(*,'(1x,a)') '------------------------------- Minimizing trace in the basis of the atomic orbitals'

  ! Allocate the local arrays that are hold by all processes.
  allocate(coeff(orbsig%norb,orbs%norb), stat=istat)
  call memocc(istat, coeff, 'coeff', subname)

  call nullify_matrixMinimization(matmin)


  ! Determine the number of processes we need, which will be stored in ip%nproc.
  ! This is the only variable in ip that all processes (also those which do not
  ! participate in the minimization of the trace) will know. The other variables
  ! are known only by the active processes and will be determined by initializeInguessParameters.
  if(lin%norbsPerProcIG>orbs%norb) then
      norbTarget=orbs%norb
  else
     norbTarget=lin%norbsperProcIG
  end if
  ip%nproc=ceiling(dble(orbs%norb)/dble(norbTarget))
  ip%nproc=min(ip%nproc,nproc)
  if(iproc==0) write(*,'(a,i0,a)') 'The minimization is performed using ', ip%nproc, ' processes.'


  ! Create the new communicator newComm.
  allocate(newID(0:ip%nproc-1), stat=istat)
  call memocc(istat, newID, 'newID', subname)
  do jproc=0,ip%nproc-1
     newID(jproc)=jproc
  end do
  call mpi_comm_group(mpi_comm_world, wholeGroup, ierr)
  call mpi_group_incl(wholeGroup, ip%nproc, newID, newGroup, ierr)
  call mpi_comm_create(mpi_comm_world, newGroup, newComm, ierr)

  ! Everything inside this if statements is only executed by the processes in newComm.
  processIf: if(iproc<ip%nproc) then

      ! Initialize the parameters for performing tha calculations in parallel.
      call initializeInguessParameters(iproc, orbs, orbsig, newComm, ip)
    
      ! Allocate the local arrays.
      call allocateArrays()

      call determineLocalizationRegions(iproc, ip%nproc, lzdig%nlr, orbsig%norb, at, onWhichAtom, &
           lin%locrad, rxyz, lin%lzd, input%hx, input%hy, input%hz, matmin%mlr)
      call extractMatrix3(iproc, ip%nproc, lin%orbs%norb, ip%norb_par(iproc), orbsig, onWhichAtomPhi, &
           ip%onWhichMPI, nlocregPerMPI, ham3, matmin, hamextract)

      call determineOverlapRegionMatrix(iproc, ip%nproc, lin%lzd, matmin%mlr, lin%orbs, orbsig, &
           onWhichAtom, onWhichAtomPhi, comom)
      !tag=1
      call initCommsMatrixOrtho(iproc, ip%nproc, lin%orbs%norb, ip%norb_par, ip%isorb_par, &
           onWhichAtomPhi, ip%onWhichMPI, tag, comom)
 
      call nullify_matrixDescriptors(mad)
      call initMatrixCompression(iproc, nproc, lzdig%nlr, orbsig, comom%noverlap, comom%overlaps, mad)
      call initCompressedMatmul3(orbsig%norb, mad)



      allocate(lcoeff(matmin%norbmax,ip%norb_par(iproc)), stat=istat)
      call memocc(istat, lcoeff, 'lcoeff', subname)
      allocate(lgrad(matmin%norbmax,ip%norb_par(iproc)), stat=istat)
      call memocc(istat, lgrad, 'lgrad', subname)
      allocate(lgradold(matmin%norbmax,ip%norb_par(iproc)), stat=istat)
      call memocc(istat, lgradold, 'lgradold', subname)
    
      ! Initialize the coefficient vectors. Put random number to places where it is
      ! reasonable (i.e. close to the atom where the basis function is centered).
      ! Make sure that the random initialization is done in the same way independent
      ! of the number of preocesses that are used.
      !call initRandomSeed(0, 1)

    
      coeffPad=0.d0
      ii=0
      do jproc=0,ip%nproc-1
          do iorb=1,ip%norb_par(jproc)
              iiAt=onWhichAtomPhi(ip%isorb_par(jproc)+iorb)
              ! Do not fill up to the boundary of the localization region, but only up to one fourth of it.
              cut=0.0625d0*lin%locrad(at%iatype(iiAt))**2
              do jorb=1,ip%norbtot
                  ii=ii+1
                  !call random_number(ttreal) ! Always call random_number to make it independent of the number of processes.
                  ttreal=builtin_rand(ii)
                  if(iproc==jproc) then
                      jjAt=onWhichAtom(jorb)
                      tt = (rxyz(1,iiat)-rxyz(1,jjAt))**2 + (rxyz(2,iiat)-rxyz(2,jjAt))**2 + (rxyz(3,iiat)-rxyz(3,jjAt))**2
                      if(tt>cut) then
                           coeffPad((iorb-1)*ip%norbtotPad+jorb)=0.d0
                      else
                          coeffPad((iorb-1)*ip%norbtotPad+jorb)=dble(ttreal)
                      end if
                  end if
              end do
          end do
      end do

    
      
      ! Initial step size for the optimization
      alpha=5.d-1
    
      ! Flag which checks convergence.
      converged=.false.
    
      if(iproc==0) write(*,'(1x,a)') '============================== optmizing coefficients =============================='
    
      ! The optimization loop.
    
      ! Transform to localization regions
      do iorb=1,ip%norb_par(iproc)
          ilr=matmin%inWhichLocregExtracted(iorb)
          if(ilr/=orbs%inWhichLocreg(iorb+orbs%isorb)) then
              write(*,'(a,2i6,3x,2i8)') &
                   'THIS IS STRANGE -- iproc, iorb, ilr, orbs%inWhichLocreg(iorb+orbs%isorb)',&
                   iproc, iorb, ilr, orbs%inWhichLocreg(iorb+orbs%isorb)
          end if
          call vectorGlobalToLocal(ip%norbtotPad, matmin%mlr(ilr), coeffPad((iorb-1)*ip%norbtotPad+1), lcoeff(1,iorb))
      end do



      if(lin%nItInguess==0) then
          ! Orthonormalize the coefficients.
          methTransformOverlap=0
          call orthonormalizeVectors(iproc, ip%nproc, newComm, lin%nItOrtho, methTransformOverlap, &
               lin%blocksize_pdsyev, lin%blocksize_pdgemm, &
               lin%orbs, onWhichAtomPhi, ip%onWhichMPI, ip%isorb_par, matmin%norbmax, ip%norb_par(iproc), ip%isorb_par(iproc), &
               lin%lzd%nlr, newComm, mad, matmin%mlr, lcoeff, comom)
      end if


      isx=1
      call initializeDIIS_inguess(isx, ip%norb_par(iproc), matmin, onWhichAtom(ip%isorb+1), ldiis)

      !!allocate(lagmatdiag(ip%norb_par(iproc)), stat=istat)
    
      iterLoop: do it=1,lin%nItInguess
    
          if (iproc==0 .and. mod(it,1)==0) then
              write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
          endif

          if(it<=2) then
              methTransformOverlap=0
          else
              methTransformOverlap=lin%methTransformOverlap
          end if
          !methTransformOverlap=0
          !methTransformOverlap=lin%methTransformOverlap
    
    
          ! Orthonormalize the coefficients.
          call orthonormalizeVectors(iproc, ip%nproc, newComm, lin%nItOrtho, methTransformOverlap, &
               lin%blocksize_pdsyev, lin%blocksize_pdgemm, &
               lin%orbs, onWhichAtomPhi, ip%onWhichMPI, ip%isorb_par, matmin%norbmax, ip%norb_par(iproc), ip%isorb_par(iproc), &
               lin%lzd%nlr, newComm, mad, matmin%mlr, lcoeff, comom)
          ilr=onWhichAtom(ip%isorb+1)

    
          ! Calculate the gradient grad.
          ilrold=0
          iilr=0
          do iorb=1,ip%norb_par(iproc)
              ilr=onWhichAtom(ip%isorb+iorb)
              iilr=matmin%inWhichLocregOnMPI(iorb)
              !!if(ilr>ilrold) then
              !!    iilr=iilr+1
              !!    ilrold=ilr
              !!end if
              call dgemv('n',matmin%mlr(ilr)%norbinlr,matmin%mlr(ilr)%norbinlr,1.d0,&
                   hamextract(1,1,iilr),matmin%norbmax,lcoeff(1,iorb),1,0.d0,lgrad(1,iorb),1)
          end do
          ilr=onWhichAtom(ip%isorb+1)

      
          if(it>1) then
              traceOld=trace
          else
              traceOld=1.d10
          end if
          ! Apply the orthoconstraint to the gradient. To do so first calculate the Lagrange
          ! multiplier matrix.

          !!do iorb=1,ip%norb_par(iproc)
          !!    lagmatdiag(iorb)=ddot(matmin%mlr(ilr)%norbinlr, lcoeff(1,iorb), 1, lgrad(1,iorb), 1)
          !!end do

          call orthoconstraintVectors(iproc, ip%nproc, methTransformOverlap, lin%correctionOrthoconstraint, &
               lin%blocksize_pdgemm, &
               lin%orbs, onWhichAtomPhi, ip%onWhichMPI, ip%isorb_par, &
               matmin%norbmax, ip%norb_par(iproc), ip%isorb_par(iproc), lin%lzd%nlr, newComm, &
               matmin%mlr, mad, lcoeff, lgrad, comom, trace)
          !!do jorb=1,matmin%norbmax
          !!    write(660+iproc,'(100f15.5)') (lgrad(jorb,iorb), iorb=1,ip%norb_par(iproc))
          !!end do
          ilr=onWhichAtom(ip%isorb+1)
          ! Calculate the gradient norm.
          fnrm=0.d0
          do iorb=1,ip%norb_par(iproc)
              ilr=onWhichAtom(ip%isorb+iorb)
              iilr=matmin%inWhichLocregOnMPI(iorb)
              fnrmArr(iorb)=ddot(matmin%mlr(ilr)%norbinlr, lgrad(1,iorb), 1, lgrad(1,iorb), 1)

              if(it>1) fnrmOvrlpArr(iorb)=ddot(matmin%mlr(ilr), lgrad(1,iorb), 1, lgradold(1,iorb), 1)
          end do
          call dcopy(ip%norb_par(iproc)*matmin%norbmax, lgrad(1,1), 1, lgradold(1,1), 1)
    
          ! Keep the gradient for the next iteration.
          if(it>1) then
              call dcopy(ip%norb_par(iproc), fnrmArr(1), 1, fnrmOldArr(1), 1)
          end if

          fnrmMax=0.d0
          meanAlpha=0.d0
          do iorb=1,ip%norb_par(iproc)

              fnrm=fnrm+fnrmArr(iorb)
              if(fnrmArr(iorb)>fnrmMax) fnrmMax=fnrmArr(iorb)
              if(it>1) then
              ! Adapt step size for the steepest descent minimization.
                  tt=fnrmOvrlpArr(iorb)/sqrt(fnrmArr(iorb)*fnrmOldArr(iorb))
                  if(tt>.9d0) then
                      alpha(iorb)=alpha(iorb)*1.1d0
                  else
                      alpha(iorb)=alpha(iorb)*.5d0
                  end if
              end if
              meanAlpha=meanAlpha+alpha(iorb)
          end do
         call mpiallred(fnrm, 1, mpi_sum, newComm, ierr)
         call mpiallred(fnrmMax, 1, mpi_max, newComm, ierr)
         fnrm=sqrt(fnrm)
         fnrmMax=sqrt(fnrmMax)
    
         ! Determine the mean step size for steepest descent iterations.
         call mpiallred(meanAlpha, 1, mpi_sum, newComm, ierr)
         meanAlpha=meanAlpha/dble(ip%norb)

          ! Precondition the gradient.
          !if(it>20) then
              do iorb=1,ip%norb_par(iproc)
                  ilr=onWhichAtom(ip%isorb+iorb)
                  iilr=matmin%inWhichLocregOnMPI(iorb)
                  !tt=ddot(matmin%mlr(ilr)%norbinlr, lcoeff(1,iorb), 1, lgrad(1,iorb), 1)
                  !tt=lagmatdiag(iorb)
                  !write(80+iproc,*) it, iorb, tt
                  !do istat=1,matmin%mlr(ilr)%norbinlr
                  !    write(90+iproc,'(3i8,2es16.8)') it, iorb, istat, lgrad(istat,iorb), lcoeff(istat,iorb)
                  !end do
                  !call preconditionGradient2(matmin%mlr(ilr)%norbinlr, matmin%norbmax, hamextract(1,1,iilr), tt, lgrad(1,iorb))
                  call preconditionGradient(matmin%mlr(ilr)%norbinlr, matmin%norbmax, hamextract(1,1,iilr), tt, lgrad(1,iorb))
              end do
          !!end if
      
    
          ! Write some informations to the screen, but only every 1000th iteration.
          if(iproc==0 .and. mod(it,1)==0) write(*,'(1x,a,es11.2,es22.13,es10.2)') 'fnrm, trace, mean alpha', &
              fnrm, trace, meanAlpha
          
          ! Check for convergence.
          if(fnrm<1.d-3) then
              if(iproc==0) write(*,'(1x,a,i0,a)') 'converged in ', it, ' iterations.'
              if(iproc==0) write(*,'(3x,a,2es14.5)') 'Final values for fnrm, trace:', fnrm, trace
              converged=.true.
              infoCoeff=it
              ! Transform back to global ragion.
              do iorb=1,ip%norb_par(iproc)
                  ilr=matmin%inWhichLocregExtracted(iorb)
                  call vectorLocalToGlobal(ip%norbtotPad, matmin%mlr(ilr), lcoeff(1,iorb), coeffPad((iorb-1)*ip%norbtotPad+1))
              end do
              exit
          end if
      
          ! Quit if the maximal number of iterations is reached.
          if(it==lin%nItInguess) then
              if(iproc==0) write(*,'(1x,a,i0,a)') 'WARNING: not converged within ', it, &
                  ' iterations! Exiting loop due to limitations of iterations.'
              if(iproc==0) write(*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, trace: ', fnrm, trace
              infoCoeff=-1
              ! Transform back to global region.
              do iorb=1,ip%norb_par(iproc)
                  ilr=matmin%inWhichLocregExtracted(iorb)
                  call vectorLocalToGlobal(ip%norbtotPad, matmin%mlr(ilr), lcoeff(1,iorb), coeffPad((iorb-1)*ip%norbtotPad+1))
              end do
              exit
          end if
    
          ! Improve the coefficients (by steepet descent).
          do iorb=1,ip%norb_par(iproc)
              ilr=onWhichAtom(ip%isorb+iorb)
              call daxpy(matmin%mlr(ilr)%norbinlr, -alpha(iorb), lgrad(1,iorb), 1, lcoeff(1,iorb), 1)
          end do
          !!if (ldiis%isx > 0) then
          !!    ldiis%mis=mod(ldiis%is,ldiis%isx)+1
          !!    ldiis%is=ldiis%is+1
          !!end if
          !!do iorb=1,ip%norb_par(iproc)
          !!    ilr=onWhichAtom(ip%isorb+iorb)
          !!    call dscal(matmin%mlr(ilr)%norbinlr, -alpha(iorb), lgrad(1,iorb), 1)
          !!end do
          !!call optimizeDIIS_inguess(iproc, ip%nproc, ip%norb_par(iproc), onWhichAtom(ip%isorb+1), matmin, lgrad, lcoeff, ldiis)

    
    
      end do iterLoop

      call deallocateDIIS(ldiis)
    
    
      if(iproc==0) write(*,'(1x,a)') '===================================================================================='
    
    
      ! Cut out the zeros
      allocate(coeff2(ip%norbtot*ip%norb_par(iproc)), stat=istat)
      call memocc(istat, coeff2, 'coeff2', subname)
      do iorb=1,ip%norb_par(iproc)
          call dcopy(ip%norbtot, coeffPad((iorb-1)*ip%norbtotPad+1), 1, coeff2((iorb-1)*ip%norbtot+1), 1)
      end do

      call deallocateArrays()

      iall=-product(shape(hamextract))*kind(hamextract)
      deallocate(hamextract, stat=istat)
      call memocc(istat, iall, 'hamextract', subname)

      call deallocate_matrixDescriptors(mad, subname)


   end if processIf
  
  call mpi_barrier(mpi_comm_world, ierr)
  
  !! Allocate coeff2 for those processes which did not allocate it
  !! during the previous if statement.
  if(iproc>=ip%nproc) then
      allocate(coeff2(1), stat=istat)
      call memocc(istat, coeff2, 'coeff2', subname)
  end if
  
  
  ! Now collect all coefficients on all processes.
  allocate(recvcounts(0:nproc-1), stat=istat)
  call memocc(istat, recvcounts, 'recvcounts', subname)
  allocate(displs(0:nproc-1), stat=istat)
  call memocc(istat, displs, 'displs', subname)
  
  !!! Send ip%norb_par and ip%norbtot to all processes.
  allocate(norb_par(0:ip%nproc-1), stat=istat)
  call memocc(istat, norb_par, 'norb_par', subname)
  if(iproc==0) then
      do jproc=0,ip%nproc-1
          norb_par(jproc)=ip%norb_par(jproc)
      end do
      norbtot=ip%norbtot
  end if
  call mpi_bcast(norb_par(0), ip%nproc, mpi_integer, 0, mpi_comm_world, ierr)
  call mpi_bcast(norbtot, 1, mpi_integer, 0, mpi_comm_world, ierr)
  
  ! Define the parameters, for the mpi_allgatherv.
  ii=0
  do jproc=0,ip%nproc-1
      recvcounts(jproc)=norbtot*norb_par(jproc)
      displs(jproc)=ii
      ii=ii+recvcounts(jproc)
  end do
  do jproc=ip%nproc,nproc-1
      recvcounts(jproc)=0
      displs(jproc)=ii
      ii=ii+recvcounts(jproc)
  end do
  if(iproc<ip%nproc) then
      sendcount=ip%norbtot*ip%norb_par(iproc)
  else
      sendcount=0
  end if

  ! Gather the coefficients.
  if (nproc > 1) then
     call mpi_allgatherv(coeff2(1), sendcount, mpi_double_precision, coeff(1,1), recvcounts, &
          displs, mpi_double_precision, mpi_comm_world, ierr)
  else
     call vcopy(sendcount,coeff2(1),1,coeff(1,1),1)
  end if

  ! Deallocate stuff which is not needed any more.
  if(iproc<ip%nproc) then
      call deallocate_inguessParameters(ip, subname)
      call deallocate_p2pCommsOrthonormalityMatrix(comom, subname)
      call deallocate_matrixMinimization(matmin,subname)

      iall=-product(shape(lcoeff))*kind(lcoeff)
      deallocate(lcoeff, stat=istat)
      call memocc(istat, iall, 'lcoeff', subname)

      iall=-product(shape(lgrad))*kind(lgrad)
      deallocate(lgrad, stat=istat)
      call memocc(istat, iall, 'lgrad', subname)

      iall=-product(shape(lgradold))*kind(lgradold)
      deallocate(lgradold, stat=istat)
      call memocc(istat, iall, 'lgradold', subname)
  end if

  iall=-product(shape(newID))*kind(newID)
  deallocate(newID, stat=istat)
  call memocc(istat, iall, 'newID', subname)

  iall=-product(shape(coeff2))*kind(coeff2)
  deallocate(coeff2, stat=istat)
  call memocc(istat, iall, 'coeff2', subname)

  iall=-product(shape(recvcounts))*kind(recvcounts)
  deallocate(recvcounts, stat=istat)
  call memocc(istat, iall, 'recvcounts', subname)

  iall=-product(shape(displs))*kind(displs)
  deallocate(displs, stat=istat)
  call memocc(istat, iall, 'displs', subname)

  iall=-product(shape(norb_par))*kind(norb_par)
  deallocate(norb_par, stat=istat)
  call memocc(istat, iall, 'norb_par', subname)

  ! Now every process has all coefficients, so we can build the linear combinations.
  ! Do this in a localized way -- TEST
  !call buildLinearCombinations(iproc, nproc, lzdig, lin%lzd, orbsig, lin%orbs, input, coeff, lchi, lphi)
  !call buildLinearCombinationsVariable(iproc, nproc, lzdig, lin%lzd, orbsig, lin%orbs, input, coeff, lchi, lphi)

  ! Now every process has all coefficients, so we can build the linear combinations.
  ! If the number of atomic orbitals is the same as the number of trace minimizing orbitals, we can use the
  ! first one (uses less memory for the overlap descriptors), otherwise the second one.
  ! So first check which version we have to use.
  same=.true.
  if(orbsig%norb==lin%orbs%norb) then
      do iorb=1,orbsig%norb
          ilr=orbsig%inWhichLocreg(iorb)
          jlr=lin%orbs%inWhichLocreg(iorb)
          if(ilr/=jlr) then
              same=.false.
              exit
          end if
      end do
  else
      same=.false.
  end if
  if(same) then
      call buildLinearCombinations(iproc, nproc, lzdig, lin%lzd, orbsig, lin%orbs, input, coeff, lchi, lin%locregShape, tag, &
           comonig, opig, madig, lphi)
  else
      !! THIS WAS THE ORIGINAL, BUT NOT WORKING.
      call buildLinearCombinationsVariable(iproc, nproc, lzdig, lin%lzd, orbsig, lin%orbs, input, coeff, lchi, tag, lphi)
  end if

  ! Deallocate the remaining local array.
  iall=-product(shape(coeff))*kind(coeff)
  deallocate(coeff, stat=istat)
  call memocc(istat, iall, 'coeff', subname)
  
  
  contains

    subroutine allocateArrays()
      allocate(coeffPad(max(ip%norbtotPad*ip%norb_par(iproc), ip%nvctrp*ip%norb)), stat=istat)
      call memocc(istat, coeffPad, 'coeffPad', subname)
      !!allocate(HamPad(ip%norbtotPad,ip%norbtotPad,at%nat), stat=istat)
      !!call memocc(istat, HamPad, 'HamPad', subname)
      allocate(grad(max(ip%norbtotPad*ip%norb_par(iproc), ip%nvctrp*ip%norb)), stat=istat)
      call memocc(istat, grad, 'grad', subname)
      allocate(gradOld(max(ip%norbtotPad*ip%norb_par(iproc), ip%nvctrp*ip%norb)), stat=istat)
      call memocc(istat, gradOld, 'gradOld', subname)
      allocate(fnrmArr(ip%norb), stat=istat)
      call memocc(istat, fnrmArr, 'fnrmArr', subname)
      allocate(fnrmOvrlpArr(ip%norb), stat=istat)
      call memocc(istat, fnrmOvrlpArr, 'fnrmOvrlpArr', subname)
      allocate(fnrmOldArr(ip%norb), stat=istat)
      call memocc(istat, fnrmOldArr, 'fnrmOldArr', subname)
      allocate(alpha(orbs%norb), stat=istat)
      call memocc(istat, alpha, 'alpha', subname)
      allocate(lagMat(orbs%norb,orbs%norb), stat=istat)
      call memocc(istat, lagMat, 'lagMat', subname)
    end subroutine allocateArrays


    subroutine deallocateArrays()
      iall=-product(shape(grad))*kind(grad)
      deallocate(grad, stat=istat)
      call memocc(istat, iall, 'grad', subname)

      iall=-product(shape(gradOld))*kind(gradOld)
      deallocate(gradOld, stat=istat)
      call memocc(istat, iall, 'gradOld', subname)

      iall=-product(shape(alpha))*kind(alpha)
      deallocate(alpha, stat=istat)
      call memocc(istat, iall, 'alpha', subname)

      iall=-product(shape(lagMat))*kind(lagMat)
      deallocate(lagMat, stat=istat)
      call memocc(istat, iall, 'lagMat', subname)
     
      iall=-product(shape(coeffPad))*kind(coeffPad)
      deallocate(coeffPad, stat=istat)
      call memocc(istat, iall, 'coeffPad', subname)

      !!iall=-product(shape(HamPad))*kind(HamPad)
      !!deallocate(HamPad, stat=istat)
      !!call memocc(istat, iall, 'HamPad', subname)

      iall=-product(shape(fnrmArr))*kind(fnrmArr)
      deallocate(fnrmArr, stat=istat)
      call memocc(istat, iall, 'fnrmArr', subname)

      iall=-product(shape(fnrmOvrlpArr))*kind(fnrmOvrlpArr)
      deallocate(fnrmOvrlpArr, stat=istat)
      call memocc(istat, iall, 'fnrmOvrlpArr', subname)

      iall=-product(shape(fnrmOldArr))*kind(fnrmOldArr)
      deallocate(fnrmOldArr, stat=istat)
      call memocc(istat, iall, 'fnrmOldArr', subname)


    end subroutine deallocateArrays


end subroutine buildLinearCombinationsLocalized3




subroutine preconditionGradient(nel, neltot, ham, cprec, grad)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: nel, neltot
  real(8),dimension(neltot,neltot),intent(in):: ham
  real(8),intent(in):: cprec
  real(8),dimension(nel),intent(inout):: grad
  
  ! Local variables
  integer:: iel, jel, info, istat, iall
  complex(8),dimension(:,:),allocatable:: mat
  complex(8),dimension(:),allocatable:: rhs
  integer,dimension(:),allocatable:: ipiv
  character(len=*),parameter:: subname='preconditionGradient'
  
  allocate(mat(nel,nel), stat=istat)
  !call memocc(istat, mat, 'mat', subname)
  allocate(rhs(nel), stat=istat)
  !call memocc(istat, mat, 'mat', subname)
  
  ! Build the matrix to be inverted (extract it from ham, which might have a larger dimension)
  do iel=1,nel
      do jel=1,nel
          mat(jel,iel) = cmplx(ham(jel,iel),0.d0,kind=8)
      end do
      mat(iel,iel)=mat(iel,iel)+cmplx(.5d0,-1.d-1,kind=8)
      !mat(iel,iel)=mat(iel,iel)-cprec
      rhs(iel)=grad(iel)
  end do
  
  
  allocate(ipiv(nel), stat=istat)
  call memocc(istat, ipiv, 'ipiv', subname)
  
  call zgesv(nel, 1, mat(1,1), nel, ipiv, rhs(1), nel, info)
  if(info/=0) then
      stop 'ERROR in dgesv'
  end if
  !call dcopy(nel, rhs(1), 1, grad(1), 1)
  do iel=1,nel
      grad(iel)=real(rhs(iel))
  end do
  
  iall=-product(shape(ipiv))*kind(ipiv)
  deallocate(ipiv, stat=istat)
  call memocc(istat, iall, 'ipiv', subname)
  
  iall=-product(shape(mat))*kind(mat)
  deallocate(mat, stat=istat)
  !call memocc(istat, iall, 'mat', subname)
  
  iall=-product(shape(rhs))*kind(rhs)
  deallocate(rhs, stat=istat)
  !call memocc(istat, iall, 'rhs', subname)

end subroutine preconditionGradient
