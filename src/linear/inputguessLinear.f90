subroutine initInputguessConfinement(iproc, nproc, at, lzd, orbs, collcom_reference, &
           Glr, input, hx, hy, hz, lin, tmbig, tmbgauss, rxyz, nscatterarr, tag)
  ! Input wavefunctions are found by a diagonalization in a minimal basis set
  ! Each processors write its initial wavefunctions into the wavefunction file
  ! The files are then read by readwave
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => initInputguessConfinement
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc,nproc
  real(gp), intent(in) :: hx, hy, hz
  type(atoms_data),intent(inout) :: at
  type(local_zone_descriptors),intent(in):: lzd
  type(orbitals_data),intent(in):: orbs
  type(collective_comms),intent(in):: collcom_reference
  type(locreg_descriptors),intent(in) :: Glr
  type(input_variables)::input
  type(linearInputParameters),intent(inout):: lin
  type(DFT_wavefunction),intent(out):: tmbig, tmbgauss
  integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(gp),dimension(3,at%nat),intent(in):: rxyz
  integer,intent(inout):: tag

  ! Local variables
  character(len=*), parameter :: subname='initInputguessConfinement'
  real(gp), dimension(:),allocatable:: locrad
  real(gp),dimension(:,:),allocatable:: locregCenter
  integer,dimension(:),allocatable:: norbsPerAt, norbsPerLocreg
  integer, parameter :: nmax=6,lmax=3,noccmax=2,nelecmax=32
  integer:: ist, iadd, ii, jj, norbtot, istat, iall, iat, nspin_ig, norbat, ityp, ilr, iorb
 


  ! Nullify the local zone descriptors.
  call nullify_local_zone_descriptors(tmbig%lzd)
  call nullify_local_zone_descriptors(tmbgauss%lzd)
  call nullify_orbitals_data(tmbig%orbs)
  call nullify_orbitals_data(tmbgauss%orbs)
  call nullify_matrixDescriptors(tmbig%mad)
  call nullify_overlapParameters(tmbig%op)
  call nullify_p2pComms(tmbig%comon)

  ! Allocate some arrays we need for the input guess.
  !!allocate(locrad(at%nat+ndebug),stat=istat)
  !!call memocc(istat,locrad,'locrad',subname)
  allocate(norbsPerAt(at%nat), stat=istat)
  call memocc(istat, norbsPerAt, 'norbsPerAt', subname)

  tmbig%lzd%hgrids(:)=lzd%hgrids(:)
  tmbgauss%lzd%hgrids(:)=lzd%hgrids(:)

  ! Number of localization regions
  !!tmbig%lzd%nlr=at%nat
  !!tmbgauss%lzd%nlr=at%nat

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

  tmbig%lzd%nlr=norbtot
  tmbgauss%lzd%nlr=at%nat

  allocate(norbsPerLocreg(tmbig%lzd%nlr), stat=istat)
  call memocc(istat, norbsPerLocreg, 'norbsPerLocreg', subname)

  allocate(locrad(tmbig%lzd%nlr),stat=istat)
  call memocc(istat,locrad,'locrad',subname)

  norbsPerLocreg=1

  ! Nullify the orbitals_data type and then determine its values.
  call nullify_orbitals_data(tmbig%orbs)
  !call orbitals_descriptors(iproc, nproc, norbtot, norbtot, 0, &
  !     input%nspin, orbs%nspinor, orbs%nkpts, orbs%kpts, orbs%kwgts, tmbig%orbs)
!!$  call orbitals_descriptors_forLinear(iproc, nproc, norbtot, norbtot, 0, &
!!$       input%nspin, orbs%nspinor, orbs%nkpts, orbs%kpts, orbs%kwgts, tmbig%orbs)
!!$  call repartitionOrbitals(iproc, nproc, tmbig%orbs%norb, tmbig%orbs%norb_par, &
!!$       tmbig%orbs%norbp, tmbig%orbs%isorb_par, tmbig%orbs%isorb, tmbig%orbs%onWhichMPI)

  call orbitals_descriptors(iproc, nproc, norbtot, norbtot, 0, &
       input%nspin, orbs%nspinor, orbs%nkpts, orbs%kpts, orbs%kwgts, tmbig%orbs,.true.) !simple repartition
  iall=-product(shape(tmbig%orbs%onwhichatom))*kind(tmbig%orbs%onwhichatom)
  deallocate(tmbig%orbs%onwhichatom,stat=istat)
  call memocc(istat,iall,'tmbig%orbs%onwhichatom',subname)

  call assignToLocreg2(iproc, nproc, tmbig%orbs%norb, tmbig%orbs%norb_par, at%nat, at%nat, &
       input%nspin, norbsPerAt, rxyz, tmbig%orbs%onwhichatom)
       !if(iproc==0) write(*,'(a,100i5)') 'tmbig%orbs%onwhichatom',tmbig%orbs%onwhichatom


  allocate(locregCenter(3,tmbig%orbs%norb), stat=istat)
  call memocc(istat, locregCenter, 'locregCenter', subname)

  ilr=0
  do iat=1,at%nat
      ityp=at%iatype(iat)
      !do iorb=1,lin%norbsPerType(ityp)
      do iorb=1,norbsPerAt(iat)
          ilr=ilr+1
          locregCenter(:,ilr)=rxyz(:,iat)
          locrad(ilr)=lin%locrad_type(ityp)
      end do
  end do
  if(ilr/=tmbig%lzd%nlr) then
      write(*,'(a,2(2x,i0))') 'ERROR: ilr/=tmbig%lzd%nlr',ilr,tmbig%lzd%nlr
      stop
  end if
  !if(iproc==0) write(*,'(a,100f)') 'lin%locrad_type',lin%locrad_type



  ! lzdig%orbs%inWhichLocreg has been allocated in orbitals_descriptors. Since it will again be allcoated
  ! in assignToLocreg2, deallocate it first.
  iall=-product(shape(tmbig%orbs%inWhichLocreg))*kind(tmbig%orbs%inWhichLocreg)
  deallocate(tmbig%orbs%inWhichLocreg,stat=istat)
  call memocc(istat,iall,'tmbig%orbs%inWhichLocreg',subname)

  ! Assign the orbitals to the localization regions.
  !!call assignToLocreg2(iproc, nproc, tmbig%orbs%norb, tmbig%orbs%norb_par, at%nat, tmbig%lzd%nlr, &
  !!     input%nspin, norbsPerAt, rxyz, tmbig%orbs%inwhichlocreg)
  call assignToLocreg2(iproc, nproc, tmbig%orbs%norb, tmbig%orbs%norb_par, at%nat, tmbig%lzd%nlr, &
       input%nspin, norbsPerLocreg, locregCenter, tmbig%orbs%inwhichlocreg)

  ! Maybe this could be moved to another subroutine? Or be omitted at all?
  allocate(tmbig%orbs%eval(orbs%norb), stat=istat)
  call memocc(istat, tmbig%orbs%eval, 'tmbig%orbs%eval', subname)
  tmbig%orbs%eval=-.5d0


  ! Nullify the locreg_descriptors and then copy Glr to it.
  call nullify_locreg_descriptors(tmbgauss%lzd%Glr)
  call copy_locreg_descriptors(Glr, tmbgauss%lzd%Glr, subname)

  ! Determine the localization regions.
  call initLocregs(iproc, nproc, tmbig%lzd%nlr, locregCenter, hx, hy, hz, tmbig%lzd, &
       tmbig%orbs, Glr, locrad, lin%locregShape)
  call copy_locreg_descriptors(Glr, tmbig%lzd%Glr, subname)

  ! Determine the localization regions for the atomic orbitals, which have a different localization radius.
  locrad=max(12.d0,maxval(lin%locrad(:)))
  call nullify_orbitals_data(tmbgauss%orbs)
  call copy_orbitals_data(tmbig%orbs, tmbgauss%orbs, subname)
  ! lzdig%orbs%inWhichLocreg has been allocated in orbitals_descriptors. Since it will again be allcoated
  ! in assignToLocreg2, deallocate it first.
  iall=-product(shape(tmbgauss%orbs%inWhichLocreg))*kind(tmbgauss%orbs%inWhichLocreg)
  deallocate(tmbgauss%orbs%inWhichLocreg,stat=istat)
  call memocc(istat,iall,'tmbgauss%orbs%inWhichLocreg',subname)
  ! Assign the orbitals to the localization regions.
  call assignToLocreg2(iproc, nproc, tmbgauss%orbs%norb, tmbgauss%orbs%norb_par, at%nat, tmbgauss%lzd%nlr, &
       input%nspin, norbsPerAt, rxyz, tmbgauss%orbs%inwhichlocreg)

  call initLocregs(iproc, nproc, tmbgauss%lzd%nlr, rxyz, input%hx, input%hy, input%hz, tmbgauss%lzd, &
       tmbgauss%orbs, Glr, locrad, lin%locregShape)

  ! Initialize the parameters needed for the orthonormalization of the atomic orbitals.
  !!! Attention: this is initialized for lzdGauss and not for lzdig!
  call initCommsOrtho(iproc, nproc, input%nspin, hx, hy, hz, tmbig%lzd, tmbig%lzd, tmbig%orbs, tmbig%orbs, &
       tmbig%orbs%inWhichLocreg, lin%locregShape, tmbig%op, tmbig%comon)

  ! Initialize the parameters needed for communicationg the potential.
  call copy_locreg_descriptors(Glr, tmbig%lzd%Glr, subname)
  call nullify_p2pComms(tmbig%comgp)
  call initialize_communication_potential(iproc, nproc, nscatterarr, tmbig%orbs, tmbig%lzd, tmbig%comgp)


  call initMatrixCompression(iproc, nproc, tmbig%lzd%nlr, tmbig%orbs, &
       tmbig%op%noverlaps, tmbig%op%overlaps, tmbig%mad)
  call initCompressedMatmul3(tmbig%orbs%norb, tmbig%mad)

  call nullify_collective_comms(tmbig%collcom)
  call init_collective_comms(iproc, nproc, tmbig%orbs, tmbig%lzd, tmbig%collcom, collcom_reference)

  ! Deallocate the local arrays.
  iall=-product(shape(locrad))*kind(locrad)
  deallocate(locrad,stat=istat)
  call memocc(istat,iall,'locrad',subname)
  iall=-product(shape(norbsPerAt))*kind(norbsPerAt)
  deallocate(norbsPerAt,stat=istat)
  call memocc(istat,iall,'norbsPerAt',subname)
  iall=-product(shape(norbsPerLocreg))*kind(norbsPerLocreg)
  deallocate(norbsPerLocreg,stat=istat)
  call memocc(istat,iall,'norbsPerLocreg',subname)
  iall=-product(shape(locregCenter))*kind(locregCenter)
  deallocate(locregCenter,stat=istat)
  call memocc(istat,iall,'locregCenter',subname)

END SUBROUTINE initInputguessConfinement





!>   input guess wavefunction diagonalization
subroutine inputguessConfinement(iproc, nproc, at, &
     input, hx, hy, hz, lzd, lorbs, rxyz, denspot, rhopotold,&
     nlpspd, proj, GPU, lphi,orbs,tmb)
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
  real(gp), intent(in) :: hx, hy, hz
  type(atoms_data), intent(inout) :: at
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(GPU_pointers), intent(inout) :: GPU
  type(DFT_local_fields), intent(inout) :: denspot
  type(input_variables),intent(inout):: input
  type(local_zone_descriptors),intent(inout):: lzd
  type(orbitals_data),intent(in):: lorbs
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
  real(dp),dimension(max(lzd%glr%d%n1i*lzd%glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin),intent(inout) ::  rhopotold
  real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(out):: lphi
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(inout):: tmb

  ! Local variables
  type(gaussian_basis):: G !basis for davidson IG
  character(len=*), parameter :: subname='inputguessConfinement'
  integer :: istat,iall,iat,nspin_ig,iorb,nvirt,norbat,ilrl,ilrg,tag
  real(gp) :: hxh,hyh,hzh,eks
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(gp), dimension(:), allocatable :: locrad
  real(wp), dimension(:,:,:), pointer :: psigau
  real(8),dimension(:),allocatable:: lchi, lchi2
  real(8),dimension(:,:),allocatable::  lhchi, locregCenter
  real(8), dimension(:,:,:),allocatable:: ham3
  integer, dimension(:),allocatable:: norbsPerAt, onWhichAtomTemp, mapping, inversemapping
  logical,dimension(:),allocatable:: covered
  !real(8),dimension(:),pointer:: lpot
  integer, parameter :: nmax=6,lmax=3,noccmax=2,nelecmax=32
  logical:: isoverlap
  logical,dimension(:),allocatable:: skip
  integer :: ist,jst,jorb,iiAt,i,iadd,ii,jj,ilr,ind1,ind2,ityp,owa,owa_old,ii_old
  integer :: ldim,gdim,jlr,iiorb,ndim_lhchi
  integer :: nlocregPerMPI,jproc,jlrold,infoCoeff
  !!integer,dimension(:),allocatable:: norb_parTemp, onWhichMPITemp
  type(confpot_data), dimension(:), allocatable :: confdatarr
  type(energy_terms) :: energs
  real(dp),dimension(6) :: xcstr
  type(DFT_wavefunction):: tmbig, tmbgauss

  if (iproc == 0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------- Input Wavefunctions Creation'
  end if

  ! Initialize evrything
  tag=1
  call initInputguessConfinement(iproc, nproc, at, lzd, lorbs, tmb%collcom, lzd%glr, input, hx, hy, hz, input%lin, &
       tmbig, tmbgauss, rxyz, denspot%dpcom%nscatterarr, tag)

  ! Allocate some arrays we need for the input guess.
  allocate(norbsc_arr(at%natsc+1,input%nspin+ndebug),stat=istat)
  call memocc(istat,norbsc_arr,'norbsc_arr',subname)
  allocate(locrad(at%nat+ndebug),stat=istat)
  call memocc(istat,locrad,'locrad',subname)
  allocate(norbsPerAt(at%nat), stat=istat)
  call memocc(istat, norbsPerAt, 'norbsPerAt', subname)
  allocate(mapping(tmbig%orbs%norb), stat=istat)
  call memocc(istat, mapping, 'mapping', subname)
  allocate(covered(tmbig%orbs%norb), stat=istat)
  call memocc(istat, covered, 'covered', subname)
  allocate(inversemapping(tmbig%orbs%norb), stat=istat)
  call memocc(istat, inversemapping, 'inversemapping', subname)


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
      ii=input%lin%norbsPerType(at%iatype(iat))
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

  ! Number of localization regions.
  tmbig%lzd%nlr=norbat
  tmbgauss%lzd%nlr=at%nat

  allocate(locregCenter(3,tmbig%lzd%nlr), stat=istat)
  call memocc(istat, locregCenter, 'locregCenter', subname)

  ilr=0
  do iat=1,at%nat
      ityp=at%iatype(iat)
      do iorb=1,norbsPerAt(iat)
          ilr=ilr+1
          locregCenter(:,ilr)=rxyz(:,iat)
      end do
  end do

  ! Create the atomic orbitals in a Gaussian basis. Deallocate tmbig%orbs, since it will be
  ! recalculated in inputguess_gaussian_orbitals.


  ! This array gives a mapping from the 'natural' orbital distribution (i.e. simply counting up the atoms) to
  ! our optimized orbital distribution (determined by in orbs%inwhichlocreg).
  iiorb=0
  covered=.false.
  do iat=1,at%nat
      do iorb=1,norbsPerAt(iat)
          iiorb=iiorb+1
          ! Search the corresponding entry in inwhichlocreg
          do jorb=1,tmbgauss%orbs%norb
              if(covered(jorb)) cycle
              jlr=tmbgauss%orbs%inwhichlocreg(jorb)
              if( tmbgauss%lzd%llr(jlr)%locregCenter(1)==rxyz(1,iat) .and. &
                  tmbgauss%lzd%llr(jlr)%locregCenter(2)==rxyz(2,iat) .and. &
                  tmbgauss%lzd%llr(jlr)%locregCenter(3)==rxyz(3,iat) ) then
                  covered(jorb)=.true.
                  mapping(iiorb)=jorb
                  exit
              end if
          end do
      end do
  end do

  ! Inverse mapping
  do iorb=1,tmbgauss%orbs%norb
      do jorb=1,tmbgauss%orbs%norb
          if(mapping(jorb)==iorb) then
              inversemapping(iorb)=jorb
              exit
          end if
      end do
  end do



  nvirt=0
  call deallocate_orbitals_data(tmbgauss%orbs,subname)

  call inputguess_gaussian_orbitals_forLinear(iproc,nproc,tmbgauss%orbs%norb,at,rxyz,nvirt,nspin_ig,&
       tmbgauss%lzd%nlr, norbsPerAt, mapping, &
       lorbs,tmbgauss%orbs,norbsc_arr,locrad,G,psigau,eks)
  ! Since inputguess_gaussian_orbitals overwrites tmbig%orbs,we again have to assign the correct value (neeed due to
  ! a different orbital distribution.
  !LG: It seems that this routine is already called in the previous routine. Commenting it out should leave things unchanged
  call repartitionOrbitals(iproc,nproc,tmbgauss%orbs%norb,tmbgauss%orbs%norb_par,&
       tmbgauss%orbs%norbp,tmbgauss%orbs%isorb_par,tmbgauss%orbs%isorb,tmbgauss%orbs%onWhichMPI)


  !dimension of the wavefunctions
  call wavefunction_dimension(tmbgauss%lzd,tmbgauss%orbs)
  call wavefunction_dimension(tmbig%lzd,tmbig%orbs)


  ! Allcoate the array holding the orbitals. lchi2 are the atomic orbitals with the larger cutoff, whereas
  ! lchi are the atomic orbitals with the smaller cutoff.
  allocate(lchi2(max(tmbgauss%orbs%npsidim_orbs,tmbgauss%orbs%npsidim_comp)),stat=istat)
  call memocc(istat,lchi2,'lchi2',subname)
  lchi2=0.d0

  ! Grid spacing on fine grid.
  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz

  ! Assign the size of the orbitals to the new variable lpsidimtot.
  tmbgauss%lzd%hgrids(1)=hx
  tmbgauss%lzd%hgrids(2)=hy
  tmbgauss%lzd%hgrids(3)=hz
  ! Transform the atomic orbitals to the wavelet basis.
  lchi2=0.d0
  call gaussians_to_wavelets_new(iproc,nproc,tmbgauss%lzd,tmbgauss%orbs,G,&
       psigau(1,1,min(tmbgauss%orbs%isorb+1,tmbgauss%orbs%norb)),lchi2)

  iall=-product(shape(psigau))*kind(psigau)
  deallocate(psigau,stat=istat)
  call memocc(istat,iall,'psigau',subname)

  call deallocate_gwf(G,subname)

  !restore wavefunction dimension
  call wavefunction_dimension(tmbig%lzd,tmbig%orbs)



  allocate(lchi(max(tmbig%orbs%npsidim_orbs,tmbig%orbs%npsidim_comp)+ndebug),stat=istat)
  call memocc(istat,lchi,'lchi',subname)
  lchi=0.d0

  ! Transform chi to the localization region. This requires that the localizatin region of lchi2 is larger than that
  ! of lchi.
  ind1=1
  ind2=1
  do iorb=1,tmbgauss%orbs%norbp
      ilrl = tmbig%orbs%inWhichLocreg(tmbig%orbs%isorb+iorb)
      ilrg = tmbgauss%orbs%inWhichLocreg(tmbgauss%orbs%isorb+iorb)
      ldim=tmbig%lzd%Llr(ilrl)%wfd%nvctr_c+7*tmbig%lzd%Llr(ilrl)%wfd%nvctr_f
      gdim=tmbgauss%lzd%llr(ilrg)%wfd%nvctr_c+7*tmbgauss%lzd%llr(ilrg)%wfd%nvctr_f
      call psi_to_locreg2(iproc, nproc, ldim, gdim, tmbig%lzd%llr(ilrl), tmbgauss%lzd%llr(ilrg), lchi2(ind1), lchi(ind2))
      ind1=ind1+tmbgauss%lzd%llr(ilrg)%wfd%nvctr_c+7*tmbgauss%lzd%llr(ilrg)%wfd%nvctr_f
      ind2=ind2+tmbig%lzd%Llr(ilrl)%wfd%nvctr_c+7*tmbig%lzd%Llr(ilrl)%wfd%nvctr_f
  end do
  if(tmbgauss%orbs%norbp>0 .and. ind1/=tmbgauss%orbs%npsidim_orbs+1) then
      write(*,'(2(a,i8),i8)') 'ERROR on process ',iproc,&
           ': ind1/=tmbgauss%orbs%npsidim+1',ind1,tmbgauss%orbs%npsidim_orbs+1
      stop
  end if
  if(tmbig%orbs%norbp>0 .and. ind2/=tmbig%orbs%npsidim_orbs+1) then
      write(*,'(2(a,i8),i8)') 'ERROR on process ',iproc,&
           ': ind2/=tmbig%orbs%npsidim+1',ind2,tmbig%orbs%npsidim_orbs+1
      stop
  end if

  ! Always use the exact Loewdin method.
  call orthonormalizeAtomicOrbitalsLocalized2(iproc, nproc, 0, input%lin%nItOrtho, &
       tmbig%lzd, tmbig%orbs, tmbig%comon, &
       tmbig%op, input, tmbig%mad, tmbig%collcom, tmb%orthpar, tmb%wfnmd%bpo, lchi)

  ! Deallocate locrad, which is not used any longer.
  iall=-product(shape(locrad))*kind(locrad)
  deallocate(locrad,stat=istat)
  call memocc(istat,iall,'locrad',subname)

  !change again wavefunction dimension
  call wavefunction_dimension(tmbig%lzd,tmbig%orbs)
  call wavefunction_dimension(tmbgauss%lzd,tmbgauss%orbs)




  ! Create the potential. First calculate the charge density.
  if(iproc==0) write(*,'(1x,a)',advance='no') 'Calculating charge density...'


  call sumrho(iproc,nproc,tmbgauss%orbs,tmbgauss%lzd,&
       hxh,hyh,hzh,denspot%dpcom%nscatterarr,&
       GPU,at%sym,denspot%rhod,lchi2,denspot%rho_psi,inversemapping)
  call communicate_density(iproc,nproc,input%nspin,hxh,hyh,hzh,tmbgauss%lzd,&
       denspot%rhod,denspot%dpcom%nscatterarr,denspot%rho_psi,denspot%rhov,.false.)

  if(iproc==0) write(*,'(a)') 'done.'

  !restore wavefunction dimension
  call wavefunction_dimension(tmbig%lzd,tmbig%orbs)


  if(trim(input%lin%mixingmethod)=='dens') then
      call dcopy(max(lzd%glr%d%n1i*lzd%glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
  end if


  iall=-product(shape(lchi2))*kind(lchi2)
  deallocate(lchi2, stat=istat)
  call memocc(istat, iall, 'lchi2',subname)


  call deallocate_local_zone_descriptors(tmbgauss%lzd, subname)

  call updatePotential(iproc,nproc,at%geocode,input%ixc,input%nspin,hxh,hyh,hzh,lzd%glr,denspot,&
       energs%eh,energs%exc,energs%evxc)

!!$  
!!$  if(orbs%nspinor==4) then
!!$     !this wrapper can be inserted inside the poisson solver 
!!$     call PSolverNC(at%geocode,'D',iproc,nproc,lzd%glr%d%n1i,lzd%glr%d%n2i,lzd%glr%d%n3i,&
!!$          nscatterarr(iproc,1),& !this is n3d
!!$          input%ixc,hxh,hyh,hzh,&
!!$          rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,4)
!!$  else
!!$     !Allocate XC potential
!!$     if (nscatterarr(iproc,2) >0) then
!!$        allocate(potxc(lzd%glr%d%n1i*lzd%glr%d%n2i*nscatterarr(iproc,2)*input%nspin+ndebug),stat=istat)
!!$        call memocc(istat,potxc,'potxc',subname)
!!$     else
!!$        allocate(potxc(1+ndebug),stat=istat)
!!$        call memocc(istat,potxc,'potxc',subname)
!!$     end if
!!$
!!$     call XC_potential(at%geocode,'D',iproc,nproc,&
!!$          lzd%glr%d%n1i,lzd%glr%d%n2i,lzd%glr%d%n3i,input%ixc,hxh,hyh,hzh,&
!!$          rhopot,eexcu,vexcu,input%nspin,rhocore,potxc,xcstr)
!!$
!!$
!!$     if( iand(potshortcut,4)==0) then
!!$        call H_potential(at%geocode,'D',iproc,nproc,&
!!$             lzd%glr%d%n1i,lzd%glr%d%n2i,lzd%glr%d%n3i,hxh,hyh,hzh,&
!!$             rhopot,pkernel,pot_ion,ehart,0.0_dp,.true.)
!!$     endif
!!$
!!$
!!$     !sum the two potentials in rhopot array
!!$     !fill the other part, for spin, polarised
!!$     if (input%nspin == 2) then
!!$        call dcopy(lzd%glr%d%n1i*lzd%glr%d%n2i*nscatterarr(iproc,2),rhopot(1),1,&
!!$             rhopot(lzd%glr%d%n1i*lzd%glr%d%n2i*nscatterarr(iproc,2)+1),1)
!!$     end if
!!$     !spin up and down together with the XC part
!!$     call axpy(lzd%glr%d%n1i*lzd%glr%d%n2i*nscatterarr(iproc,2)*input%nspin,1.0_dp,potxc(1),1,&
!!$          rhopot(1),1)
!!$
!!$
!!$     iall=-product(shape(potxc))*kind(potxc)
!!$     deallocate(potxc,stat=istat)
!!$     call memocc(istat,iall,'potxc',subname)
!!$
!!$  end if

  if(trim(input%lin%mixingmethod)=='pot') then
      call dcopy(max(lzd%glr%d%n1i*lzd%glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
  end if

  !call dcopy(tmbig%orbs%npsidim,psi,1,hpsi,1)
  if (input%exctxpar == 'OP2P') energs%eexctX = uninitialized(energs%eexctX)

  
  ! Set localnorb, i.e. the number of orbitals a given process has in a specific loalization region.
  do ilr=1,tmbig%lzd%nlr
      tmbig%lzd%Llr(ilr)%localnorb=0
      do iorb=1,tmbig%orbs%norbp
          !if(tmbig%orbs%inWhichLocregp(iorb)==ilr) then
          if(tmbig%orbs%inWhichLocreg(tmbig%orbs%isorb+iorb)==ilr) then
              tmbig%lzd%Llr(ilr)%localnorb = tmbig%lzd%Llr(ilr)%localnorb+1
          end if
      end do
  end do


  ! Post the messages for the communication of the potential.
  call allocateCommunicationsBuffersPotential(tmbig%comgp, subname)
  !!call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, tmbig%comgp)
  call post_p2p_communication(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, &
       tmbig%comgp%nrecvbuf, tmbig%comgp%recvbuf, tmbig%comgp)


  ! Apply the Hamiltonian for each atom.
  ! onWhichAtomTemp indicates that all orbitals feel the confining potential
  ! centered on atom iat.
  allocate(onWhichAtomTemp(tmbig%orbs%norb), stat=istat)
  call memocc(istat,onWhichAtomTemp,'onWhichAtomTemp',subname)
  allocate(skip(lzd%nlr), stat=istat)
  call memocc(istat, skip, 'skip', subname)


  ! Determine for how many localization regions we need a Hamiltonian application.
  ndim_lhchi=0
  do ilr=1,lzd%nlr
      skip(ilr)=.true.
      do jorb=1,tmbig%orbs%norbp
          onWhichAtomTemp(tmbig%orbs%isorb+jorb)=tmbig%orbs%inwhichlocreg(ilr)
          jlr=tmbig%orbs%inWhichLocreg(tmbig%orbs%isorb+jorb)
          if(tmbig%orbs%inWhichlocreg(jorb+tmbig%orbs%isorb)/=jlr) stop 'this should not happen'
          call check_overlap_cubic_periodic(tmb%lzd%Glr,tmb%lzd%llr(ilr),tmbig%lzd%llr(jlr),isoverlap)
           if(isoverlap) then
              skip(ilr)=.false.
          end if
      end do
      if(.not.skip(ilr)) then
          ndim_lhchi=ndim_lhchi+1
      end if
  end do


  allocate(lhchi(max(tmbig%orbs%npsidim_orbs,tmbig%orbs%npsidim_comp),ndim_lhchi),stat=istat)
  call memocc(istat, lhchi, 'lhchi', subname)
  lhchi=0.d0


  if(iproc==0) write(*,'(1x,a)') 'Hamiltonian application for all locregs. This may take some time.'


  call local_potential_dimensions(tmbig%lzd,tmbig%orbs,denspot%dpcom%ngatherarr(0,1))

  !!tmbig%comgp%communication_complete=.false.
  call full_local_potential(iproc,nproc,tmbig%orbs,tmbig%lzd,2,&
       denspot%dpcom,denspot%rhov,denspot%pot_work,tmbig%comgp)

  tmbig%lzd%hgrids(1)=hx
  tmbig%lzd%hgrids(2)=hy
  tmbig%lzd%hgrids(3)=hz


  allocate(tmbig%lzd%doHamAppl(tmbig%lzd%nlr), stat=istat)
  call memocc(istat, tmbig%lzd%doHamAppl, 'tmbig%lzd%doHamAppl', subname)
  ii=0
  owa_old=-1
  do ilr=1,lzd%nlr
      tmbig%lzd%doHamAppl=.false.
      skip(ilr)=.true.
      do jorb=1,tmbig%orbs%norbp
          onWhichAtomTemp(tmbig%orbs%isorb+jorb)=lorbs%onwhichatom(ilr)
          owa=lorbs%onwhichatom(ilr)
          jlr=tmbig%orbs%inWhichLocreg(tmbig%orbs%isorb+jorb)
          call check_overlap_cubic_periodic(tmb%lzd%Glr,tmb%lzd%llr(lorbs%inwhichlocreg(ilr)),&
               tmbig%lzd%llr(jlr),isoverlap)
          if(isoverlap) then
              tmbig%lzd%doHamAppl(jlr)=.true.
              skip(ilr)=.false.
          else
              tmbig%lzd%doHamAppl(jlr)=.false.
          end if
      end do
      if(iproc==0) write(*,'(3x,a,i0,a)', advance='no') 'locreg ', ilr, '... '

      if(.not.skip(ilr)) then
          ii=ii+1
          if(input%lin%nItInguess>0) then
             allocate(confdatarr(tmbig%orbs%norbp))
             call define_confinement_data(confdatarr,tmbig%orbs,rxyz,at,&
                  hx,hy,hz,input%lin%confpotorder,&
                  input%lin%potentialprefac_lowaccuracy,tmbig%lzd,onWhichAtomTemp)
             call to_zero(tmbig%orbs%npsidim_orbs,lhchi(1,ii))
             if(owa/=owa_old) then
                 call LocalHamiltonianApplication(iproc,nproc,at,tmbig%orbs,&
                      tmbig%lzd,confdatarr,denspot%dpcom%ngatherarr,denspot%pot_work,lchi,lhchi(1,ii),&
                      energs,input%SIC,GPU,.false.,&
                      pkernel=denspot%pkernelseq)
                 call NonLocalHamiltonianApplication(iproc,at,tmbig%orbs,&
                      rxyz,&
                      proj,tmbig%lzd,nlpspd,lchi,lhchi(1,ii),energs%eproj)
                 ii_old=ii
             else
                 call dcopy(tmbig%orbs%npsidim_orbs, lhchi(1,ii_old), 1, lhchi(1,ii), 1)
             end if
             deallocate(confdatarr)
          end if
      else
          if(iproc==0) write(*,'(3x,a)', advance='no') 'no Hamiltonian application required... '
      end if
      if(iproc==0) write(*,'(a)') 'done.'
      owa_old=owa
  end do


  ! Deallocate the buffers needed for communication the potential.
  call deallocateCommunicationsBuffersPotential(tmbig%comgp, subname)
  ! Deallocate the parameters needed for the communication of the potential.
  call deallocate_p2pComms(tmbig%comgp, subname)

  iall=-product(shape(denspot%pot_work))*kind(denspot%pot_work)
  deallocate(denspot%pot_work, stat=istat)
  call memocc(istat, iall, 'denspot%pot_work', subname)
   if(ii/=ndim_lhchi) then
      write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': ii/=ndim_lhchi',ii,ndim_lhchi
      stop
  end if



  ! Calculate the number of different matrices that have to be stored on a given MPI process.
  jlrold=0
  nlocregPerMPI=0
  do jorb=1,lorbs%norb
      jlr=lorbs%inWhichLocreg(jorb)
      jproc=lorbs%onWhichMPI(jorb)
      if(iproc==jproc) then
          if(jlr/=jlrold) then
              nlocregPerMPI=nlocregPerMPI+1
              jlrold=jlr
          end if
      end if
  end do





  ! Calculate the Hamiltonian matrix.
  allocate(ham3(tmbig%orbs%norb,tmbig%orbs%norb,nlocregPerMPI), stat=istat)
  call memocc(istat,ham3,'ham3',subname)

  if(input%lin%nItInguess>0) then
      call getHamiltonianMatrix6(iproc, nproc, lzd, tmbig%lzd, tmbig%orbs, lorbs, &
           input, hx, hy, hz, tmbig%orbs%inWhichLocreg, ndim_lhchi, &
           nlocregPerMPI, lchi, lhchi, skip, tmbig%mad, input%lin%memoryForCommunOverlapIG, input%lin%locregShape, tag, ham3)
  end if

  iall=-product(shape(lhchi))*kind(lhchi)
  deallocate(lhchi, stat=istat)
  call memocc(istat, iall, 'lhchi',subname)


  ! Build the orbitals phi as linear combinations of the atomic orbitals.
  call buildLinearCombinationsLocalized3(iproc, nproc, nlocregPerMPI, hx, hy, hz, &
           tmb, tmbig, at, input, lchi, locregCenter, rxyz, ham3, lphi)

  ! Calculate the coefficients
  call allocateCommunicationsBuffersPotential(tmb%comgp, subname)
  call post_p2p_communication(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, &
       tmb%comgp%nrecvbuf, tmb%comgp%recvbuf, tmb%comgp)
  call get_coeff(iproc,nproc,lzd,orbs,at,rxyz,denspot,GPU,infoCoeff,energs%ebs,nlpspd,proj,&
       tmb%wfnmd%bpo%blocksize_pdsyev,tmb%wfnmd%bpo%nproc_pdsyev,&
       hx,hy,hz,input%SIC,tmb,tmb)
  ! Deallocate the buffers needed for the communication of the potential.
  call deallocateCommunicationsBuffersPotential(tmb%comgp, subname)
  

  if(iproc==0) write(*,'(1x,a)') '------------------------------------------------------------- Input guess generated.'


  ! Deallocate all local arrays.

  ! Deallocate all types that are not needed any longer.
  call deallocate_local_zone_descriptors(tmbig%lzd, subname)
  call deallocate_orbitals_data(tmbig%orbs, subname)
  call deallocate_orbitals_data(tmbgauss%orbs, subname)
  call deallocate_matrixDescriptors(tmbig%mad, subname)
  call deallocate_overlapParameters(tmbig%op, subname)
  call deallocate_p2pComms(tmbig%comon, subname)
  call deallocate_collective_comms(tmbig%collcom, subname)

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

  iall=-product(shape(skip))*kind(skip)
  deallocate(skip, stat=istat)
  call memocc(istat, iall, 'skip',subname)

  iall=-product(shape(ham3))*kind(ham3)
  deallocate(ham3, stat=istat)
  call memocc(istat, iall, 'ham3',subname)

  iall=-product(shape(mapping))*kind(mapping)
  deallocate(mapping, stat=istat)
  call memocc(istat, iall, 'mapping',subname)

  iall=-product(shape(covered))*kind(covered)
  deallocate(covered, stat=istat)
  call memocc(istat, iall, 'covered',subname)

  iall=-product(shape(inversemapping))*kind(inversemapping)
  deallocate(inversemapping, stat=istat)
  call memocc(istat, iall, 'inversemapping',subname)

  iall=-product(shape(locregCenter))*kind(locregCenter)
  deallocate(locregCenter, stat=istat)
  call memocc(istat, iall, 'locregCenter',subname)

END SUBROUTINE inputguessConfinement





subroutine orthonormalizeAtomicOrbitalsLocalized2(iproc, nproc, methTransformOverlap, nItOrtho, &
           lzd, orbs, comon, op, input, mad, collcom, orthpar, bpo, lchi)

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
integer,intent(in):: iproc, nproc, methTransformOverlap, nItOrtho
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: orbs
type(input_variables),intent(in):: input
type(p2pComms),intent(inout):: comon
type(overlapParameters),intent(inout):: op
type(matrixDescriptors),intent(in):: mad
type(collective_comms),intent(in):: collcom
type(orthon_data),intent(in):: orthpar
type(basis_performance_options),intent(in):: bpo
real(8),dimension(orbs%npsidim_orbs),intent(inout):: lchi

! Local variables
integer:: iorb, jorb, istat, iall, lwork, info, nvctrp, ierr, tag, ilr
real(8),dimension(:,:),allocatable:: ovrlp
character(len=*),parameter:: subname='orthonormalizeAtomicOrbitalsLocalized2'


! Initialize the communication parameters.
allocate(ovrlp(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, ovrlp, 'ovrlp', subname)

call orthonormalizeLocalized(iproc, nproc, methTransformOverlap, nItOrtho, &
     orbs, op, comon, lzd, mad, collcom, orthpar, bpo, lchi, ovrlp)

iall=-product(shape(ovrlp))*kind(ovrlp)
deallocate(ovrlp, stat=istat)
call memocc(istat, iall, 'ovrlp', subname)


end subroutine orthonormalizeAtomicOrbitalsLocalized2

subroutine getHamiltonianMatrix6(iproc, nproc, lzd, lzdig, orbsig, orbs, &
input, hx, hy, hz, onWhichAtom, ndim_lhchi, nlocregPerMPI, lchi, lhchi, skip, mad, memoryForCommunOverlapIG, locregShape, &
tagout, ham)
use module_base
use module_types
use module_interfaces, exceptThisOne => getHamiltonianMatrix6
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, ndim_lhchi, nlocregPerMPI
real(gp),intent(in) :: hx, hy, hz
type(local_zone_descriptors),intent(in):: lzd, lzdig
type(orbitals_data),intent(in):: orbsig, orbs
type(input_variables),intent(in):: input
integer,dimension(orbsig%norb),intent(in):: onWhichAtom
!real(8),dimension(lzdig%orbs%npsidim),intent(in):: chi
!real(8),dimension(lzdig%orbs%npsidim,nat),intent(in):: hchi
real(8),dimension(orbsig%npsidim_orbs),intent(in):: lchi
real(8),dimension(orbsig%npsidim_orbs,ndim_lhchi),intent(in):: lhchi
logical,dimension(lzd%nlr),intent(in):: skip
type(matrixDescriptors),intent(in):: mad
integer,intent(in):: memoryForCommunOverlapIG
character(len=1),intent(in):: locregShape
integer,intent(inout):: tagout
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
noverlaps=min(noverlaps,lzd%nlr)


allocate(hamTempCompressed(max(sendcounts(iproc),1),noverlaps), stat=istat)
call memocc(istat, hamTempCompressed, 'hamTempCompressed', subname)
allocate(hamTempCompressed2(mad%nvctr,nlocregPerMPI), stat=istat)
call memocc(istat, hamTempCompressed2, 'ovrlpCompressed2', subname)

allocate(hamTemp(orbsig%norb,orbsig%norb), stat=istat)
call memocc(istat, hamTemp, 'hamTemp', subname)

! Initialize the parameters for calculating the matrix.
call nullify_p2pComms(comon)
call initCommsOrtho(iproc, nproc, input%nspin, hx, hy, hz, lzdig, lzdig, orbsig, orbsig, &
     onWhichAtom, locregShape, op, comon)


call allocateCommuncationBuffersOrtho(comon, subname)

! Put lphi in the sendbuffer, i.e. lphi will be sent to other processes' receive buffer.
! Then post the messages and gather them.
call extractOrbital3(iproc, nproc, orbsig, orbsig, orbsig%npsidim_orbs, onWhichAtom, lzdig, lzdig, op, op, &
     lchi, comon%nsendBuf, comon%sendBuf)
call post_p2p_communication(iproc, nproc, comon%nsendbuf, comon%sendbuf, &
     comon%nrecvbuf, comon%recvbuf, comon)
call wait_p2p_communication(iproc, nproc, comon)



if(iproc==0) write(*,'(1x,a)') 'Calculating Hamiltonian matrix for all atoms. This may take some time.'
ilrold=0
iilr=0
iatold=0
ii=0
imatold=1
imat=0
!do iat=1,lzdig%nlr
do iat=1,lzd%nlr

    if(iproc==0) write(*,'(3x,a,i0,a)', advance='no') 'Calculating matrix for locreg ', iat, '... '

    ioverlap=mod(iat-1,noverlaps)+1


    ! Put lhphi to the sendbuffer, so we can the calculate <lphi|lhphi>
    if(.not.skip(iat)) then
        ii=ii+1
        call extractOrbital3(iproc, nproc, orbsig, orbsig, orbsig%npsidim_orbs, onWhichAtom, lzdig, lzdig, op, op, &
             lhchi(1,ii), comon%nsendBuf, comon%sendBuf)
        call calculateOverlapMatrix3Partial(iproc, nproc, orbsig, op, onWhichAtom, comon%nsendBuf, comon%sendBuf, &
             comon%nrecvBuf, comon%recvBuf, mad, hamTemp(1,1))
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

        ! nshift indicates how much the following loops do i=1,jj deviate from the outer loop on iat.
        nshift=iat-jj

        ! First determine the number of sends / receives for each process.
        nsend=0
        nrecv=0
        ilrold=-1
        jjprocold=-1
        do iorb=1,orbs%norb
            ilr=orbs%inWhichLocreg(iorb)
            jjproc=orbs%onWhichMPI(iorb)
            do iioverlap=1,jj
                iiat=iioverlap+nshift
                if(iproc<nproc) then
                    if(ilr==ilrold .and. jjproc==jjprocold) cycle !Otherwise we would communicate the same again
                    if(ilr==iiat) then
                        ! Send this matrix to process jproc.
                        if(iproc==jjproc) then
                            do jproc=0,nproc-1
                                !nrecv=nrecv+1
                                if(orbs%norb_par(jproc,0)>0) nrecv=nrecv+1 !otherwise process jproc has no data and should not communicate...
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
        ! WARNING: Here we don't use the standard and unique tags available through p2p_tags. It should not matter
        ! since there is no other p2p communication going on at the moment, but still this could be improved...
        tag0=1
        isend=0
        irecv=0
        ilrold=-1
        jjprocold=-1
        do iorb=1,orbs%norb
            ilr=orbs%inWhichLocreg(iorb)
            jjproc=orbs%onWhichMPI(iorb)
            do iioverlap=1,jj
                iiat=iioverlap+nshift
                ! Check whether this MPI needs this matrix. Since only nprocTemp processes will be involved
                ! in calculating the input guess, this check has to be done only for those processes.
                if(ilr==ilrold .and. jjproc==jjprocold) cycle
                if(ilr==iiat) then
                   ! Send to process jproc
                   if(iproc==jjproc .and. nproc > 1) then
                      imat=imat+1
                      do jproc=0,nproc-1
                         if(orbs%norb_par(jproc,0)==0) cycle !process jproc has no data and should not communicate...
                         tag=tag0+jproc
                         irecv=irecv+1
                         call mpi_irecv(hamTempCompressed2(displs(jproc)+1,imat), sendcounts(jproc), &
                              mpi_double_precision, jproc, tag, mpi_comm_world, recvrequests(irecv), ierr)
                      end do
                      tag=tag0+iproc
                      isend=isend+1
                      call mpi_isend(hamTempCompressed(1,iorb), sendcounts(iproc), &
                           mpi_double_precision, jjproc, tag, mpi_comm_world, sendrequests(isend), ierr)
                   else if (nproc >1) then
                      tag=tag0+iproc
                      isend=isend+1
                      call mpi_isend(hamTempCompressed(1,iorb), sendcounts(iproc), &
                           mpi_double_precision, jjproc, tag, mpi_comm_world, sendrequests(isend), ierr)
                   else if (nproc == 1) then
                      imat=imat+1
                      call vcopy(sendcounts(iproc),hamTempCompressed(1,iorb),1,&
                           hamTempCompressed2(displs(iproc)+1,imat),1)
                   end if
                    tag0=tag0+1
                end if
            end do
            ilrold=ilr
            jjprocold=jjproc
        end do

        ! Wait for the communication to complete
        if (nproc > 1) then
          isend=0
          waitForSend: do
             if(isend==nsend) exit waitForSend
             call mpi_waitany(nsend-isend, sendrequests(1), ind, mpi_status_ignore, ierr)
             isend=isend+1
             do i=ind,nsend-isend
                sendrequests(i)=sendrequests(i+1)
             end do
          end do waitForSend

          irecv=0
          waitForRecv: do
             if(irecv==nrecv) exit waitForrecv
             call mpi_waitany(nrecv-irecv, recvrequests(1), ind, mpi_status_ignore, ierr)
             irecv=irecv+1
             do i=ind,nrecv-irecv
                recvrequests(i)=recvrequests(i+1)
             end do
          end do waitForRecv
        end if

     ! Uncompress the matrices
     do i=imatold,imat
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


if(imat/=nlocregPerMPI .and. nproc >1) then
  write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': imat/=nlocregPerMPI',imat,nlocregPerMPI
  stop
end if
call deallocate_overlapParameters(op, subname)
call deallocate_p2pComms(comon, subname)

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





subroutine determineLocalizationRegions(iproc, nproc, nlr, norb, at, onWhichAtomAll, locrad, rxyz, lzd, lzdig, hx, hy, hz, mlr)
use module_base
use module_types
use module_interfaces, exceptThisOne => determineLocalizationRegions
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nlr, norb
type(atoms_data),intent(in):: at
integer,dimension(norb),intent(in):: onWhichAtomAll
real(8),dimension(at%nat),intent(in):: locrad
real(8),dimension(3,nlr),intent(in):: rxyz
type(local_zone_descriptors),intent(in):: lzd, lzdig
real(8),intent(in):: hx, hy, hz
type(matrixLocalizationRegion),dimension(:),pointer,intent(out):: mlr

! Local variables
integer:: ilr, jlr, jorb, ii, istat
!integer::  is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3
real(8):: cut, tt
!logical:: ovrlpx, ovrlpy, ovrlpz
logical:: isoverlap
character(len=*),parameter:: subname='determineLocalizationRegions'


allocate(mlr(nlr), stat=istat)
do ilr=1,nlr
  call nullify_matrixLocalizationRegion(mlr(ilr))
end do
! Count for each localization region the number of matrix elements within the cutoff.
do ilr=1,nlr
  mlr(ilr)%norbinlr=0
  do jorb=1,norb
     jlr=onWhichAtomAll(jorb)
     call check_overlap_cubic_periodic(lzd%Glr,lzd%Llr(ilr),lzdig%Llr(jlr),isoverlap)     
     if(isoverlap) then
        mlr(ilr)%norbinlr=mlr(ilr)%norbinlr+1
     end if
  end do
  allocate(mlr(ilr)%indexInGlobal(mlr(ilr)%norbinlr), stat=istat)
  call memocc(istat, mlr(ilr)%indexInGlobal, 'mlr(ilr)%indexInGlobal', subname)
end do


! Now determine the indices of the elements with an overlap.
do ilr=1,nlr
  ii=0
  do jorb=1,norb
     jlr=onWhichAtomAll(jorb)
      call check_overlap_cubic_periodic(lzd%Glr,lzd%Llr(ilr),lzdig%Llr(jlr),isoverlap)
      if(isoverlap) then
        ii=ii+1
        mlr(ilr)%indexInGlobal(ii)=jorb
      end if
  end do
  if(ii/=mlr(ilr)%norbinlr) then
     write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ', iproc, ': ii/=mlr(ilr)%norbinlr', ii, mlr(ilr)%norbinlr
  end if
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
  novrlp=0
  do jorbout=1,orbs%norb
     jlr=onWhichAtomPhi(jorbout)
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
  comom%noverlap(iorbout)=novrlp
end do


allocate(comom%overlaps(maxval(comom%noverlap(:)),orbs%norb), stat=istat)
call memocc(istat, comom%overlaps, 'comom%overlaps', subname)
do iorbout=1,orbs%norb
  ilr=orbs%inwhichlocreg(iorbout)
  comom%overlaps(:,iorbout)=0
  novrlp=0
  do jorbout=1,orbs%norb
     jlr=onWhichAtomPhi(jorbout)
     ! Check whether there is a common element.
     outloop2: do iorb=1,mlr(ilr)%norbinlr
        iiorb=mlr(ilr)%indexInGlobal(iorb)
        do jorb=1,mlr(jlr)%norbinlr
           jjorb=mlr(jlr)%indexInGlobal(jorb)
           if(iiorb==jjorb) then
              novrlp=novrlp+1
              comom%overlaps(novrlp,iorbout)=jorbout
              exit outloop2
           end if
        end do
     end do outloop2
  end do
end do

allocate(comom%olr(maxval(comom%noverlap(:)),lzd%nlr), stat=istat)
do ilr=1,lzd%nlr
  do iorb=1,maxval(comom%noverlap(:))
     call nullify_matrixLocalizationRegion(comom%olr(iorb,ilr))
  end do
end do



! Now determine which orbitals (corresponding to basis functions) will be in the overlap localization region.
ilrold=-1
do iorbout=1,orbs%norb
  ilr=orbs%inwhichlocreg(iorbout)
  if(ilr==ilrold) cycle
  ilrold=ilr
  comom%olr(:,ilr)%norbinlr=0
  do jorbout=1,comom%noverlap(iorbout)
     jjorb=comom%overlaps(jorbout,iorbout)
     jlr=onWhichAtomPhi(jjorb)
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



! Determine the indices to switch from global region to localization region.
ilrold=-1
do iorbout=1,orbs%norb
  ilr=orbs%inwhichlocreg(iorbout)
  if(ilr==ilrold) cycle
  ilrold=ilr
  do jorbout=1,comom%noverlap(iorbout)
     jjorb=comom%overlaps(jorbout,iorbout)
     jlr=onWhichAtomPhi(jjorb)
     ! Check whether there is a common element.
     kkorb=0
     do iorb=1,mlr(ilr)%norbinlr
        iiorb=mlr(ilr)%indexInGlobal(iorb)
        do jorb=1,mlr(jlr)%norbinlr
           jjorb=mlr(jlr)%indexInGlobal(jorb)
           if(iiorb==jjorb) then
              kkorb=kkorb+1
              comom%olr(jorbout,ilr)%indexInGlobal(kkorb)=jorb
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
ilrold=-1
do iorbout=1,orbs%norb
  ilr=orbs%inwhichlocreg(iorbout)
  if(ilr==ilrold) cycle
  ilrold=ilr
  do iorb=1,comom%noverlap(iorbout)
     jorb=comom%overlaps(iorb,iorbout)
     jlr=onWhichAtomPhi(jorb)
     comom%olrForExpansion(1,iorb,ilr)=jlr
     do korb=1,comom%noverlap(jorb)
        kkorb=comom%overlaps(korb,jorb)
        klr=onWhichAtomPhi(kkorb)
        if(klr==ilr) then
           comom%olrForExpansion(2,iorb,ilr)=korb
        end if
     end do
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
     do korb=1,comom%noverlap(jjorb)
        comom%noverlapProc(jproc)=comom%noverlapProc(jproc)+1
     end do
     jlrold=jlr
  end do
end do

allocate(comom%comarr(6,maxval(comom%noverlapProc(:)),0:nproc-1), stat=istat)
call memocc(istat, comom%comarr, 'comom%comarr', subname)
allocate(comom%overlapsProc(maxval(comom%noverlapProc(:)),0:nproc-1), stat=istat)
call memocc(istat, comom%overlapsProc, 'comom%overlapsProc', subname)

comom%nsendBuf=0
comom%nrecvBuf=0
do jproc=0,nproc-1
  jkorb=0
  jlrold=0
  do jorb=1,norb_par(jproc)
     jjorb=isorb_par(jproc)+jorb
     jlr=onWhichAtomPhi(jjorb)
     if(jlr==jlrold) cycle
     do korb=1,comom%noverlap(jjorb)
        jkorb=jkorb+1
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
        end if
        if(iproc==mpidest) then
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



subroutine orthonormalizeVectors(iproc, nproc, comm, nItOrtho, methTransformOverlap, &
  orbs, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, isorb, nlr, newComm, mad, mlr, vec, comom, &
  collcom, orthpar, bpo)
use module_base
use module_types
use module_interfaces, exceptThisOne => orthonormalizeVectors
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, comm, nItOrtho, methTransformOverlap
integer,intent(in):: norbmax, norbp, isorb, nlr, newComm
type(orbitals_data),intent(in):: orbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom, onWhichMPI
integer,dimension(0:nproc-1),intent(in):: isorb_par
type(matrixDescriptors),intent(in):: mad
type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
real(8),dimension(norbmax,norbp),intent(inout):: vec
type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
type(collective_comms),intent(in):: collcom
type(orthon_data),intent(in):: orthpar
type(basis_performance_options),intent(in):: bpo

! Local variables
integer:: noverlaps, iorb, iiorb, ilr, istat, ilrold, jorb, iall, it, iorbmax, jorbmax, i, ist, ncnt
real(8):: tt, dnrm2, dev
real(8),dimension(:,:),allocatable:: vecOvrlp, ovrlp
character(len=*),parameter:: subname='orthonormalizeVectors'
real(8),dimension(orbs%norb):: vecglobal
real(8),dimension(:),allocatable:: psit_c, psit_f, psittemp_c, psittemp_f, vec_compr

allocate(ovrlp(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, ovrlp, 'ovrlp', subname)

do it=1,nItOrtho

    if(bpo%communication_strategy_overlap==COMMUNICATION_P2P) then
        
        noverlaps=0
        ilrold=0
        do iorb=1,norbp
          iiorb=isorb+iorb
          ilr=onWhichAtom(iiorb)
          if(ilr/=ilrold) then
             noverlaps=noverlaps+comom%noverlap(iiorb)
          end if
          ilrold=ilr
        end do
        allocate(vecOvrlp(norbmax,noverlaps), stat=istat)
        call memocc(istat, vecOvrlp, 'vecOvrlp', subname)
        
        
        call extractToOverlapregion(iproc, nproc, orbs%norb, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, vec, comom)
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
    
    else if(bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then
    
        allocate(vec_compr(collcom%ndimpsi_c), stat=istat)
        call memocc(istat, vec_compr, 'vec_compr', subname)
        ist=1
        do iorb=1,orbs%norbp
            iiorb=orbs%isorb+iorb
            ilr=orbs%inwhichlocreg(iiorb)
            call dcopy(mlr(ilr)%norbinlr, vec(1,iorb), 1, vec_compr(ist), 1)
            ist=ist+mlr(ilr)%norbinlr
        end do
    
        allocate(psit_c(sum(collcom%nrecvcounts_c)), stat=istat)
        call memocc(istat, psit_c, 'psit_c', subname)
        allocate(psit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
        call memocc(istat, psit_f, 'psit_f', subname)
        call transpose_localized(iproc, nproc, orbs, collcom, vec_compr, psit_c, psit_f)
        call calculate_overlap_transposed(iproc, nproc, orbs, mad, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp)
    
        dev=0.d0
        iorbmax=0
        jorbmax=0
        do iorb=1,orbs%norb
           do jorb=1,orbs%norb
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
    
    end if
    
    
    call overlapPowerMinusOneHalf(iproc, nproc, comm, methTransformOverlap, &
         orthpar%blocksize_pdsyev, orthpar%blocksize_pdgemm, orbs%norb, mad, ovrlp)
        
    if(bpo%communication_strategy_overlap==COMMUNICATION_P2P) then
        
        call orthonormalLinearCombinations(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, &
             orbs%norb, comom, mlr, onWhichAtom, vecOvrlp, ovrlp, vec)
        
        iall=-product(shape(vecOvrlp))*kind(vecOvrlp)
        deallocate(vecOvrlp, stat=istat)
        call memocc(istat, iall, 'vecOvrlp', subname)
    
    else if(bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then
    
    
        allocate(psittemp_c(sum(collcom%nrecvcounts_c)), stat=istat)
        call memocc(istat, psittemp_c, 'psittemp_c', subname)
        allocate(psittemp_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
        call memocc(istat, psittemp_f, 'psittemp_f', subname)
        call dcopy(sum(collcom%nrecvcounts_c), psit_c, 1, psittemp_c, 1)
        call dcopy(7*sum(collcom%nrecvcounts_f), psit_f, 1, psittemp_f, 1)
        call build_linear_combination_transposed(orbs%norb, ovrlp, collcom, psittemp_c, psittemp_f, .true., psit_c, psit_f)
        call untranspose_localized(iproc, nproc, orbs, collcom, psit_c, psit_f, vec_compr)
    
        ist=1
        do iorb=1,orbs%norbp
            iiorb=orbs%isorb+iorb
            ilr=orbs%inwhichlocreg(iiorb)
            call dcopy(mlr(ilr)%norbinlr, vec_compr(ist), 1, vec(1,iorb), 1)
            ist=ist+mlr(ilr)%norbinlr
        end do
    
        iall=-product(shape(vec_compr))*kind(vec_compr)
        deallocate(vec_compr, stat=istat)
        call memocc(istat, iall, 'vec_compr', subname)
    
        iall=-product(shape(psittemp_c))*kind(psittemp_c)
        deallocate(psittemp_c, stat=istat)
        call memocc(istat, iall, 'psittemp_c', subname)
        iall=-product(shape(psittemp_f))*kind(psittemp_f)
        deallocate(psittemp_f, stat=istat)
        call memocc(istat, iall, 'psittemp_f', subname)
        iall=-product(shape(psit_c))*kind(psit_c)
        deallocate(psit_c, stat=istat)
        call memocc(istat, iall, 'psit_c', subname)
        iall=-product(shape(psit_f))*kind(psit_f)
        deallocate(psit_f, stat=istat)
        call memocc(istat, iall, 'psit_f', subname)
    
    end if

    ! Normalize the vectors
    do iorb=1,orbs%norbp
       iiorb=orbs%isorb+iorb
       ilr=orbs%inwhichlocreg(iiorb)
       ncnt=mlr(ilr)%norbinlr
       tt=dnrm2(ncnt, vec(1,iorb), 1)
       call dscal(ncnt, 1/tt, vec(1,iorb), 1)
    end do

end do

iall=-product(shape(ovrlp))*kind(ovrlp)
deallocate(ovrlp, stat=istat)
call memocc(istat, iall, 'ovrlp', subname)

end subroutine orthonormalizeVectors




subroutine orthoconstraintVectors(iproc, nproc, methTransformOverlap, correctionOrthoconstraint, orbs, &
           onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, isorb, nlr, newComm, mlr, mad, vec, grad, comom, trace, &
           collcom, orthpar, bpo)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => orthoconstraintVectors
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, methTransformOverlap, correctionOrthoconstraint, norbmax
  integer,intent(in):: norbp, isorb, nlr, newComm
  type(orbitals_data),intent(in):: orbs
  integer,dimension(orbs%norb),intent(in):: onWhichAtom, onWhichMPI
  integer,dimension(0:nproc-1),intent(in):: isorb_par
  type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(norbmax,norbp),intent(inout):: vec, grad
  type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
  real(8),intent(out):: trace
  type(collective_comms),intent(in):: collcom
  type(orthon_data),intent(in):: orthpar
  type(basis_performance_options),intent(in):: bpo
  
  ! Local variables
  integer:: noverlaps, iorb, iiorb, ilr, istat, ilrold, jorb, iall, ijorb, ncount, jjorb, ist
  real(8),dimension(:,:),allocatable:: gradOvrlp, vecOvrlp, lagmat, ovrlp
  character(len=*),parameter:: subname='orthoconstraintVectors'
  real(8):: ddot
  real(8),dimension(:,:),allocatable:: ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans
  real(8),dimension(:),allocatable:: psit_c, psit_f, hpsit_c, hpsit_f, vec_compr, grad_compr
  
  allocate(ovrlp_minus_one_lagmat(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp_minus_one_lagmat, 'ovrlp_minus_one_lagmat', subname)
  allocate(ovrlp_minus_one_lagmat_trans(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp_minus_one_lagmat_trans, 'ovrlp_minus_one_lagmat_trans', subname)
  allocate(ovrlp(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp, 'ovrlp', subname)
  allocate(lagmat(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, lagmat, 'lagmat', subname)
  
  
  if(bpo%communication_strategy_overlap==COMMUNICATION_P2P) then
      
      noverlaps=0
      ilrold=0
      do iorb=1,norbp
        iiorb=isorb+iorb
        ilr=onWhichAtom(iiorb)
        if(ilr/=ilrold) then
           noverlaps=noverlaps+comom%noverlap(iiorb)
        end if
        ilrold=ilr
      end do
      allocate(gradOvrlp(norbmax,noverlaps), stat=istat)
      call memocc(istat, gradOvrlp, 'gradOvrlp', subname)
      allocate(vecOvrlp(norbmax,noverlaps), stat=istat)
      call memocc(istat, vecOvrlp, 'vecOvrlp', subname)
      
      call extractToOverlapregion(iproc, nproc, orbs%norb, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, grad, comom)
      
      call postCommsVectorOrthonormalizationNew(iproc, nproc, newComm, comom)
      call gatherVectorsNew(iproc, nproc, comom)
      
      call expandFromOverlapregion(iproc, nproc, isorb, norbp, orbs, onWhichAtom, comom, norbmax, noverlaps, gradOvrlp)
      
      ! Calculate the Lagrange multiplier matrix <vec|grad>.
      call calculateOverlap(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, orbs%norb, comom, mlr, onWhichAtom, vec,&
          gradOvrlp, newComm, lagmat)
      
      ! Now we also have to calculate the overlap matrix.
      call extractToOverlapregion(iproc, nproc, orbs%norb, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, vec, comom)
      call postCommsVectorOrthonormalizationNew(iproc, nproc, newComm, comom)
      call gatherVectorsNew(iproc, nproc, comom)
      call expandFromOverlapregion(iproc, nproc, isorb, norbp, orbs, onWhichAtom, comom, norbmax, noverlaps, vecOvrlp)
      call calculateOverlap(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, orbs%norb, comom, mlr, onWhichAtom, vec,&
          vecOvrlp, newComm, ovrlp)
  
  else if(bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then
  
      allocate(vec_compr(collcom%ndimpsi_c), stat=istat)
      call memocc(istat, vec_compr, 'vec_compr', subname)
      allocate(grad_compr(collcom%ndimpsi_c), stat=istat)
      call memocc(istat, grad_compr, 'grad_compr', subname)
      ist=1
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          call dcopy(mlr(ilr)%norbinlr, vec(1,iorb), 1, vec_compr(ist), 1)
          call dcopy(mlr(ilr)%norbinlr, grad(1,iorb), 1, grad_compr(ist), 1)
          ist=ist+mlr(ilr)%norbinlr
      end do
  
      allocate(psit_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, psit_c, 'psit_c', subname)
      allocate(psit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, psit_f, 'psit_f', subname)
      allocate(hpsit_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, hpsit_c, 'hpsit_c', subname)
      allocate(hpsit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, hpsit_f, 'hpsit_f', subname)
      call transpose_localized(iproc, nproc, orbs, collcom, vec_compr, psit_c, psit_f)
      call transpose_localized(iproc, nproc, orbs, collcom, grad_compr, hpsit_c, hpsit_f)
      call calculate_overlap_transposed(iproc, nproc, orbs, mad, collcom, psit_c, psit_c, psit_f, psit_f, ovrlp)
      call calculate_overlap_transposed(iproc, nproc, orbs, mad, collcom, psit_c, hpsit_c, psit_f, hpsit_f, lagmat)
  
  end if
      
  ! Now apply the orthoconstraint.
  call applyOrthoconstraintNonorthogonal2(iproc, nproc, methTransformOverlap, bpo%blocksize_pdgemm, &
           correctionOrthoconstraint, orbs, &
           lagmat, ovrlp, mad, ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans)
      
  if(bpo%communication_strategy_overlap==COMMUNICATION_P2P) then
  
      ilrold=-1
      ijorb=0
      do iorb=1,norbp
          iiorb=isorb+iorb
          ilr=onWhichAtom(iiorb)
          if(ilr==ilrold) then
              ! Set back the index of lphiovrlp, since we again need the same orbitals.
              ijorb=ijorb-comom%noverlap(iiorb)
          end if
          ncount=mlr(ilr)%norbinlr
          do jorb=1,comom%noverlap(iiorb)
              ijorb=ijorb+1
              jjorb=comom%overlaps(jorb,iiorb)
              call daxpy(ncount, -.5d0*ovrlp_minus_one_lagmat(jjorb,iiorb), vecOvrlp(1,ijorb), 1, grad(1,iorb), 1)
              call daxpy(ncount, -.5d0*ovrlp_minus_one_lagmat_trans(jjorb,iiorb), vecOvrlp(1,ijorb), 1, grad(1,iorb), 1)
          end do
          ilrold=ilr
      end do
  
      iall=-product(shape(gradOvrlp))*kind(gradOvrlp)
      deallocate(gradOvrlp, stat=istat)
      call memocc(istat, iall, 'gradOvrlp', subname)
      
      iall=-product(shape(vecOvrlp))*kind(vecOvrlp)
      deallocate(vecOvrlp, stat=istat)
      call memocc(istat, iall, 'vecOvrlp', subname)
  
  else if(bpo%communication_strategy_overlap==COMMUNICATION_COLLECTIVE) then
  
        do iorb=1,orbs%norb
            do jorb=1,orbs%norb
                ovrlp(jorb,iorb)=-.5d0*ovrlp_minus_one_lagmat(jorb,iorb)
            end do
        end do
        call build_linear_combination_transposed(orbs%norb, ovrlp, collcom, psit_c, psit_f, .false., hpsit_c, hpsit_f)
        
        do iorb=1,orbs%norb
            do jorb=1,orbs%norb
                ovrlp(jorb,iorb)=-.5d0*ovrlp_minus_one_lagmat_trans(jorb,iorb)
            end do
        end do
        call build_linear_combination_transposed(orbs%norb, ovrlp, collcom, psit_c, psit_f, .false., hpsit_c, hpsit_f)
  
      call untranspose_localized(iproc, nproc, orbs, collcom, psit_c, psit_f, vec_compr)
      call untranspose_localized(iproc, nproc, orbs, collcom, hpsit_c, hpsit_f, grad_compr)
  
  
      ist=1
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          call dcopy(mlr(ilr)%norbinlr, vec_compr(ist), 1, vec(1,iorb), 1)
          call dcopy(mlr(ilr)%norbinlr, grad_compr(ist), 1, grad(1,iorb), 1)
          ist=ist+mlr(ilr)%norbinlr
      end do
  
      iall=-product(shape(vec_compr))*kind(vec_compr)
      deallocate(vec_compr, stat=istat)
      call memocc(istat, iall, 'vec_compr', subname)
      iall=-product(shape(grad_compr))*kind(grad_compr)
      deallocate(grad_compr, stat=istat)
      call memocc(istat, iall, 'grad_compr', subname)
  
      iall=-product(shape(psit_c))*kind(psit_c)
      deallocate(psit_c, stat=istat)
      call memocc(istat, iall, 'psit_c', subname)
      iall=-product(shape(psit_f))*kind(psit_f)
      deallocate(psit_f, stat=istat)
      call memocc(istat, iall, 'psit_f', subname)
      iall=-product(shape(hpsit_c))*kind(hpsit_c)
      deallocate(hpsit_c, stat=istat)
      call memocc(istat, iall, 'hpsit_c', subname)
      iall=-product(shape(hpsit_f))*kind(hpsit_f)
      deallocate(hpsit_f, stat=istat)
      call memocc(istat, iall, 'hpsit_f', subname)
  
  end if
  
  trace=0.d0
  do iorb=1,orbs%norb
    trace=trace+lagmat(iorb,iorb)
  end do
  
  
  iall=-product(shape(ovrlp_minus_one_lagmat))*kind(ovrlp_minus_one_lagmat)
  deallocate(ovrlp_minus_one_lagmat, stat=istat)
  call memocc(istat, iall, 'ovrlp_minus_one_lagmat', subname)
  iall=-product(shape(ovrlp_minus_one_lagmat_trans))*kind(ovrlp_minus_one_lagmat_trans)
  deallocate(ovrlp_minus_one_lagmat_trans, stat=istat)
  call memocc(istat, iall, 'ovrlp_minus_one_lagmat_trans', subname)
  
  iall=-product(shape(lagmat))*kind(lagmat)
  deallocate(lagmat, stat=istat)
  call memocc(istat, iall, 'lagmat', subname)
  
  iall=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp, stat=istat)
  call memocc(istat, iall, 'ovrlp', subname)

end subroutine orthoconstraintVectors


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
!!comom%communComplete=.false.
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
        !write(*,'(6(a,i0))') 'process ', mpidest, ' receives ', ncount,&
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
        !write(*,'(6(a,i0))') 'process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, &
        !     ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
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
  if(comom%nsend>0) then
      waitLoopSend: do
         call mpi_waitany(comom%nsend-nsend, comom%requests(1,1), ind, mpi_status_ignore, ierr)
         nsend=nsend+1
         do i=ind,comom%nsend-nsend
            comom%requests(i,1)=comom%requests(i+1,1)
         end do
         if(nsend==comom%nsend) exit waitLoopSend
      end do waitLoopSend
  end if


  nrecv=0
  if(comom%nrecv>0) then
      waitLoopRecv: do
         call mpi_waitany(comom%nrecv-nrecv, comom%requests(1,2), ind, mpi_status_ignore, ierr)
         nrecv=nrecv+1
         do i=ind,comom%nrecv-nrecv
            comom%requests(i,2)=comom%requests(i+1,2)
         end do
         if(nrecv==comom%nrecv) exit waitLoopRecv
      end do waitLoopRecv
  end if
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

if(norbp>0) call dcopy(norbmax*norbp, vec(1,1), 1, vecTemp(1,1), 1)

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



subroutine buildLinearCombinations_new(iproc, nproc, lzdig, lzd, orbsig, orbs, coeff, lchi, &
           collcomig, collcom, lphi)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => buildLinearCombinations_new
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(local_zone_descriptors),intent(in):: lzdig, lzd
  type(orbitals_data),intent(in):: orbsig, orbs
  real(8),dimension(orbsig%norb,orbs%norb),intent(in):: coeff
  real(8),dimension(orbsig%npsidim_orbs),intent(in):: lchi
  type(collective_comms),intent(in):: collcomig, collcom
  real(8),dimension(orbs%npsidim_orbs),intent(out):: lphi
  
  ! Local variables
  integer:: istat, iall, ist, jst, ilr, ilrold, iorb, iiorb, ncount, jorb, jjorb, korb, kkorb, klr, iwa, kwa
  integer:: i0, ipt, ii, jj, i, j, j0
  real(8),dimension(:),allocatable:: lchiovrlp, lchiovrlp2, chit_c, chit_f, phit_c, phit_f
  character(len=*),parameter:: subname='buildLinearCombinations_new'
  real(8),dimension(:,:),allocatable:: ttmat
  real(8):: tt1, tt2, tt3
  type(p2pComms):: comon_tmb_ig, comon_ig_tmb
  type(overlapParameters):: op_tmb_ig, op_ig_tmb
  
  allocate(chit_c(sum(collcomig%nrecvcounts_c)), stat=istat)
  call memocc(istat, chit_c, 'chit_c', subname)
  allocate(chit_f(7*sum(collcomig%nrecvcounts_f)), stat=istat)
  call memocc(istat, chit_f, 'chit_f', subname)
  call transpose_localized(iproc, nproc, orbsig, collcomig, lchi, chit_c, chit_f, lzdig)
  
  allocate(phit_c(sum(collcom%nrecvcounts_c)), stat=istat)
  call memocc(istat, phit_c, 'phit_c', subname)
  allocate(phit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
  call memocc(istat, phit_f, 'phit_f', subname)
  
  phit_c=0.d0
  phit_f=0.d0
  
  i0=0
  j0=0
  do ipt=1,collcom%nptsp_c 
      ii=collcom%norb_per_gridpoint_c(ipt) 
      jj=collcomig%norb_per_gridpoint_c(ipt) 
      do i=1,ii
          iiorb=collcom%indexrecvorbital_c(i0+i)
          do j=1,jj
              jjorb=collcomig%indexrecvorbital_c(j0+j)
              phit_c(i0+i)=phit_c(i0+i)+coeff(jjorb,iiorb)*chit_c(j0+j)
          end do
      end do
      i0=i0+ii
      j0=j0+jj
  end do
  
  i0=0
  j0=0
  do ipt=1,collcom%nptsp_f
      ii=collcom%norb_per_gridpoint_f(ipt)
      jj=collcomig%norb_per_gridpoint_f(ipt)
      do i=1,ii
          iiorb=collcom%indexrecvorbital_f(i0+i)
          do j=1,jj
              jjorb=collcomig%indexrecvorbital_f(j0+j)
              phit_f(7*(i0+i)-6)=phit_f(7*(i0+i)-6)+coeff(jjorb,iiorb)*chit_f(7*(j0+j)-6)
              phit_f(7*(i0+i)-5)=phit_f(7*(i0+i)-5)+coeff(jjorb,iiorb)*chit_f(7*(j0+j)-5)
              phit_f(7*(i0+i)-4)=phit_f(7*(i0+i)-4)+coeff(jjorb,iiorb)*chit_f(7*(j0+j)-4)
              phit_f(7*(i0+i)-3)=phit_f(7*(i0+i)-3)+coeff(jjorb,iiorb)*chit_f(7*(j0+j)-3)
              phit_f(7*(i0+i)-2)=phit_f(7*(i0+i)-2)+coeff(jjorb,iiorb)*chit_f(7*(j0+j)-2)
              phit_f(7*(i0+i)-1)=phit_f(7*(i0+i)-1)+coeff(jjorb,iiorb)*chit_f(7*(j0+j)-1)
              phit_f(7*(i0+i)-0)=phit_f(7*(i0+i)-0)+coeff(jjorb,iiorb)*chit_f(7*(j0+j)-0)
          end do
      end do
      i0=i0+ii
      j0=j0+jj
  end do
  
  call untranspose_localized(iproc, nproc, orbs, collcom, phit_c, phit_f, lphi, lzd)
  
  iall=-product(shape(chit_c))*kind(chit_c)
  deallocate(chit_c, stat=istat)
  call memocc(istat, iall, 'chit_c', subname)
  iall=-product(shape(chit_f))*kind(chit_f)
  deallocate(chit_f, stat=istat)
  call memocc(istat, iall, 'chit_f', subname)
  
  iall=-product(shape(phit_c))*kind(phit_c)
  deallocate(phit_c, stat=istat)
  call memocc(istat, iall, 'phit_c', subname)
  iall=-product(shape(phit_f))*kind(phit_f)
  deallocate(phit_f, stat=istat)
  call memocc(istat, iall, 'phit_f', subname)
  
end subroutine buildLinearCombinations_new


subroutine buildLinearCombinationsLocalized3(iproc, nproc, nlocregPerMPI, hx, hy, hz, &
           tmb, tmbig, at, input, lchi, locregCenter, rxyz, ham, lphi)
use module_base
use module_types
use module_interfaces, exceptThisOne => buildLinearCombinationsLocalized3
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nlocregPerMPI
real(gp), intent(in) :: hx, hy, hz
type(DFT_wavefunction),intent(in):: tmb, tmbig
type(atoms_data),intent(in):: at
type(input_variables),intent(in):: input
real(8),dimension(tmbig%orbs%npsidim_orbs),intent(in):: lchi
real(8),dimension(3,tmbig%lzd%nlr),intent(in):: locregCenter
real(8),dimension(3,at%nat),intent(in):: rxyz
real(8),dimension(tmbig%orbs%norb,tmbig%orbs%norb,nlocregPerMPI),intent(in):: ham
real(8),dimension(tmb%orbs%npsidim_orbs),intent(out):: lphi

! Local variables
integer:: iorb, jorb, iall, istat, ierr, infoCoeff
integer:: it, iiAt, jjAt, methTransformOverlap, tag
real(8),dimension(:),allocatable:: alpha, coeffPad, fnrmArr, fnrmOvrlpArr, fnrmOldArr
real(8),dimension(:,:),allocatable:: coeff, lagMat, lcoeff, lgrad, lgradold
integer,dimension(:),allocatable:: recvcounts, displs, norb_par
real(8):: ddot, cosangle, tt, dnrm2, fnrm, meanAlpha, cut, trace, traceOld, fnrmMax
logical:: converged
character(len=*),parameter:: subname='buildLinearCombinationsLocalized3'
real(4):: ttreal, builtin_rand
integer:: norbtot, isx, iiiat
integer:: ii, jproc, sendcount, ilr, iilr, ilrold, jlr
real(8),dimension(:,:,:),pointer:: hamextract
type(p2pCommsOrthonormalityMatrix):: comom
type(matrixMinimization):: matmin
type(localizedDIISParameters):: ldiis
type(matrixDescriptors):: mad
type(collective_comms):: collcom_vectors



  if(iproc==0) write(*,'(1x,a)') '------------------------------- Minimizing trace in the basis of the atomic orbitals'

  call nullify_matrixMinimization(matmin)
  tag=1


  if(iproc==0) write(*,'(a,i0,a)') 'The minimization is performed using ', nproc, ' processes.'


  ! Allocate the local arrays.
  call allocateArrays()

  call determineLocalizationRegions(iproc, nproc, tmb%lzd%nlr, tmbig%orbs%norb, at, tmbig%orbs%inwhichlocreg, &
       input%lin%locrad, locregCenter, tmb%lzd, tmbig%lzd, hx, hy, hz, matmin%mlr)
  call extractMatrix3(iproc, nproc, tmb%orbs%norb, tmb%orbs%norbp, tmbig%orbs, tmb%orbs%inwhichlocreg, &
       tmb%orbs%onwhichmpi, nlocregPerMPI, ham, matmin, hamextract)

  call determineOverlapRegionMatrix(iproc, nproc, tmb%lzd, matmin%mlr, tmb%orbs, tmbig%orbs, &
       tmbig%orbs%inwhichlocreg, tmb%orbs%inwhichlocreg, comom)

  call initCommsMatrixOrtho(iproc, nproc, tmb%orbs%norb, tmb%orbs%norb_par, tmb%orbs%isorb_par, &
       tmb%orbs%inwhichlocreg, tmb%orbs%onwhichmpi, tag, comom)

  call nullify_matrixDescriptors(mad)
  call initMatrixCompression(iproc, nproc, tmbig%lzd%nlr, tmb%orbs, comom%noverlap, comom%overlaps, mad)
  call initCompressedMatmul3(tmb%orbs%norb, mad)

  call nullify_collective_comms(collcom_vectors)
  call init_collective_comms_vectors(iproc, nproc, tmb%lzd%nlr, tmb%orbs, tmbig%orbs, matmin%mlr, collcom_vectors)


  allocate(lcoeff(matmin%norbmax,tmb%orbs%norbp), stat=istat)
  call memocc(istat, lcoeff, 'lcoeff', subname)
  allocate(lgrad(matmin%norbmax,tmb%orbs%norbp), stat=istat)
  call memocc(istat, lgrad, 'lgrad', subname)
  allocate(lgradold(matmin%norbmax,tmb%orbs%norbp), stat=istat)
  call memocc(istat, lgradold, 'lgradold', subname)

  ! Initialize the coefficient vectors. Put random number to places where it is
  ! reasonable (i.e. close to the atom where the basis function is centered).
  ! Make sure that the random initialization is done in the same way independent
  ! of the number of preocesses that are used.
  !call initRandomSeed(0, 1)


  coeffPad=0.d0
  ii=0
  do jproc=0,nproc-1
      do iorb=1,tmb%orbs%norb_par(jproc,0)
          iiAt=tmb%orbs%inwhichlocreg(tmb%orbs%isorb_par(jproc)+iorb)
          iiiAt=tmb%orbs%onwhichatom(tmb%orbs%isorb_par(jproc)+iorb)
          ! Do not fill up to the boundary of the localization region, but only up to one fifth of it.
          !cut=0.0625d0*lin%locrad(at%iatype(iiAt))**2
          cut=0.04d0*input%lin%locrad(at%iatype(iiiAt))**2
          do jorb=1,tmbig%orbs%norb
              ii=ii+1
              ttreal=builtin_rand(ii)
              if(iproc==jproc) then
                  jjAt=tmbig%orbs%inwhichlocreg(jorb)
                  tt = (rxyz(1,iiiat)-locregCenter(1,jjAt))**2 + &
                       (rxyz(2,iiiat)-locregCenter(2,jjAt))**2 + &
                       (rxyz(3,iiiat)-locregCenter(3,jjAt))**2
                  if(tt>cut) then
                       coeffPad((iorb-1)*tmbig%orbs%norb+jorb)=0.d0
                  else
                      coeffPad((iorb-1)*tmbig%orbs%norb+jorb)=dble(ttreal)
                  end if
              end if
          end do
      end do
  end do


  
  ! Initial step size for the optimization
  alpha=5.d-1

  ! Flag which checks convergence.
  converged=.false.

  if(iproc==0) write(*,'(1x,a)') '============================== optimizing coefficients =============================='

  ! The optimization loop.

  ! Transform to localization regions
  do iorb=1,tmb%orbs%norbp
      ilr=matmin%inWhichLocregExtracted(iorb)
      if(ilr/=tmb%orbs%inWhichLocreg(iorb+tmb%orbs%isorb)) then
          write(*,'(a,2i6,3x,2i8)') &
               'THIS IS STRANGE -- iproc, iorb, ilr, tmb%orbs%inWhichLocreg(iorb+tmb%orbs%isorb)',&
               iproc, iorb, ilr, tmb%orbs%inWhichLocreg(iorb+tmb%orbs%isorb)
      end if
      call vectorGlobalToLocal(tmbig%orbs%norb, matmin%mlr(ilr), coeffPad((iorb-1)*tmbig%orbs%norb+1), lcoeff(1,iorb))
  end do



  if(input%lin%nItInguess==0) then
      ! Orthonormalize the coefficients.
      methTransformOverlap=0
      call orthonormalizeVectors(iproc, nproc, mpi_comm_world, input%lin%nItOrtho, methTransformOverlap, &
           tmb%orbs, tmb%orbs%inwhichlocreg, tmb%orbs%onwhichmpi, tmb%orbs%isorb_par, &
           matmin%norbmax, tmb%orbs%norbp, tmb%orbs%isorb_par(iproc), &
           tmb%lzd%nlr, mpi_comm_world, mad, matmin%mlr, lcoeff, comom, collcom_vectors, tmb%orthpar, tmb%wfnmd%bpo)
  end if


  !!isx=1
  !!if(tmb%orbs%norbp>0) then
  !!    ! otherwise it makes no sense...
  !!    call initializeDIIS_inguess(isx, tmb%orbs%norbp, matmin, tmb%orbs%inwhichlocreg(tmb%orbs%isorb+1), ldiis)
  !!end if



  iterLoop: do it=1,input%lin%nItInguess

      if (iproc==0 .and. mod(it,1)==0) then
          write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
      endif

      if(it<=2) then
          methTransformOverlap=0
      else
          methTransformOverlap=input%lin%methTransformOverlap
      end if


      ! Orthonormalize the coefficients.
      call orthonormalizeVectors(iproc, nproc, mpi_comm_world, input%lin%nItOrtho, methTransformOverlap, &
           tmb%orbs, tmb%orbs%inwhichlocreg, tmb%orbs%onwhichmpi, tmb%orbs%isorb_par, &
           matmin%norbmax, tmb%orbs%norbp, tmb%orbs%isorb_par(iproc), &
           tmb%lzd%nlr, mpi_comm_world, mad, matmin%mlr, lcoeff, comom, collcom_vectors, tmb%orthpar, tmb%wfnmd%bpo)


      ! Calculate the gradient grad.
      ilrold=0
      iilr=0
      do iorb=1,tmb%orbs%norbp
          ilr=tmb%orbs%inwhichlocreg(tmb%orbs%isorb+iorb)
          iilr=matmin%inWhichLocregOnMPI(iorb)
          call dgemv('n',matmin%mlr(ilr)%norbinlr,matmin%mlr(ilr)%norbinlr,1.d0,&
               hamextract(1,1,iilr),matmin%norbmax,lcoeff(1,iorb),1,0.d0,lgrad(1,iorb),1)
      end do

  
      if(it>1) then
          traceOld=trace
      else
          traceOld=1.d10
      end if
      ! Apply the orthoconstraint to the gradient. To do so first calculate the Lagrange
      ! multiplier matrix.
      call orthoconstraintVectors(iproc, nproc, methTransformOverlap, input%lin%correctionOrthoconstraint, &
           tmb%orbs, tmb%orbs%inwhichlocreg, tmb%orbs%onwhichmpi, tmb%orbs%isorb_par, &
           matmin%norbmax, tmb%orbs%norbp, tmb%orbs%isorb_par(iproc), tmb%lzd%nlr, mpi_comm_world, &
           matmin%mlr, mad, lcoeff, lgrad, comom, trace, collcom_vectors, tmb%orthpar, tmb%wfnmd%bpo)
      ! Calculate the gradient norm.
      fnrm=0.d0
      do iorb=1,tmb%orbs%norbp
          ilr=tmb%orbs%inwhichlocreg(tmb%orbs%isorb+iorb)
          iilr=matmin%inWhichLocregOnMPI(iorb)
          fnrmArr(iorb)=ddot(matmin%mlr(ilr)%norbinlr, lgrad(1,iorb), 1, lgrad(1,iorb), 1)

          if(it>1) fnrmOvrlpArr(iorb)=ddot(matmin%mlr(ilr)%norbinlr, lgrad(1,iorb), 1, lgradold(1,iorb), 1)
      end do
      if(tmb%orbs%norbp>0) call dcopy(tmb%orbs%norbp*matmin%norbmax, lgrad(1,1), 1, lgradold(1,1), 1)

      ! Keep the gradient for the next iteration.
      if(it>1 .and. (tmb%orbs%norbp>0)) then
          call dcopy(tmb%orbs%norbp, fnrmArr(1), 1, fnrmOldArr(1), 1)
      end if

      fnrmMax=0.d0
      meanAlpha=0.d0
      do iorb=1,tmb%orbs%norbp

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
     call mpiallred(fnrm, 1, mpi_sum, mpi_comm_world, ierr)
     call mpiallred(fnrmMax, 1, mpi_max, mpi_comm_world, ierr)
     fnrm=sqrt(fnrm)
     fnrmMax=sqrt(fnrmMax)

     ! Determine the mean step size for steepest descent iterations.
     call mpiallred(meanAlpha, 1, mpi_sum, mpi_comm_world, ierr)
     meanAlpha=meanAlpha/dble(tmb%orbs%norb)

      ! Precondition the gradient.
      do iorb=1,tmb%orbs%norbp
          ilr=tmb%orbs%inwhichlocreg(tmb%orbs%isorb+iorb)
          iilr=matmin%inWhichLocregOnMPI(iorb)
          call preconditionGradient(matmin%mlr(ilr)%norbinlr, matmin%norbmax, hamextract(1,1,iilr), tt, lgrad(1,iorb))
      end do
  

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
          do iorb=1,tmb%orbs%norbp
              ilr=matmin%inWhichLocregExtracted(iorb)
              call vectorLocalToGlobal(tmbig%orbs%norb, matmin%mlr(ilr), lcoeff(1,iorb), coeffPad((iorb-1)*tmbig%orbs%norb+1))
          end do
          exit
      end if
  
      ! Quit if the maximal number of iterations is reached.
      if(it==input%lin%nItInguess) then
          if(iproc==0) write(*,'(1x,a,i0,a)') 'WARNING: not converged within ', it, &
              ' iterations! Exiting loop due to limitations of iterations.'
          if(iproc==0) write(*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, trace: ', fnrm, trace
          infoCoeff=-1
          ! Transform back to global region.
          do iorb=1,tmb%orbs%norbp
              ilr=matmin%inWhichLocregExtracted(iorb)
              call vectorLocalToGlobal(tmbig%orbs%norb, matmin%mlr(ilr), lcoeff(1,iorb), coeffPad((iorb-1)*tmbig%orbs%norb+1))
          end do
          exit
      end if

      ! Improve the coefficients (by steepet descent).
      do iorb=1,tmb%orbs%norbp
          ilr=tmb%orbs%inwhichlocreg(tmb%orbs%isorb+iorb)
          call daxpy(matmin%mlr(ilr)%norbinlr,-alpha(iorb), lgrad(1,iorb), 1, lcoeff(1,iorb), 1)
      end do


  end do iterLoop

  if(iproc==0) write(*,'(1x,a)') '===================================================================================='



  iall=-product(shape(hamextract))*kind(hamextract)
  deallocate(hamextract, stat=istat)
  call memocc(istat, iall, 'hamextract', subname)

  call deallocate_matrixDescriptors(mad, subname)


  ! Now collect all coefficients on all processes.
  allocate(recvcounts(0:nproc-1), stat=istat)
  call memocc(istat, recvcounts, 'recvcounts', subname)
  allocate(displs(0:nproc-1), stat=istat)
  call memocc(istat, displs, 'displs', subname)
  
  ! Define the parameters, for the mpi_allgatherv.
  ii=0
  do jproc=0,nproc-1
      recvcounts(jproc)=tmbig%orbs%norb*tmb%orbs%norb_par(jproc,0)
      displs(jproc)=ii
      ii=ii+recvcounts(jproc)
  end do
  sendcount=tmbig%orbs%norb*tmb%orbs%norbp

  ! Allocate the local arrays that are hold by all processes.
  allocate(coeff(tmbig%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, coeff, 'coeff', subname)

  ! Gather the coefficients.
  if (nproc > 1) then
     call mpi_allgatherv(coeffPad(1), sendcount, mpi_double_precision, coeff(1,1), recvcounts, &
          displs, mpi_double_precision, mpi_comm_world, ierr)
  else
     call vcopy(sendcount,coeffPad(1),1,coeff(1,1),1)
  end if

  call deallocateArrays()

  ! Deallocate stuff which is not needed any more.
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

  iall=-product(shape(recvcounts))*kind(recvcounts)
  deallocate(recvcounts, stat=istat)
  call memocc(istat, iall, 'recvcounts', subname)

  iall=-product(shape(displs))*kind(displs)
  deallocate(displs, stat=istat)
  call memocc(istat, iall, 'displs', subname)


  ! Now every process has all coefficients, so we can build the linear combinations.
  call buildLinearCombinations_new(iproc, nproc, tmbig%lzd, tmb%lzd, tmbig%orbs, tmb%orbs, coeff, lchi, &
       tmbig%collcom, tmb%collcom, lphi)

  ! Deallocate the remaining local array.
  iall=-product(shape(coeff))*kind(coeff)
  deallocate(coeff, stat=istat)
  call memocc(istat, iall, 'coeff', subname)

  call deallocate_collective_comms(collcom_vectors, subname)
  
  
  contains

    subroutine allocateArrays()
      allocate(coeffPad(max(tmbig%orbs%norb*tmb%orbs%norbp,1)), stat=istat)
      call memocc(istat, coeffPad, 'coeffPad', subname)
      allocate(fnrmArr(tmb%orbs%norb), stat=istat)
      call memocc(istat, fnrmArr, 'fnrmArr', subname)
      allocate(fnrmOvrlpArr(tmb%orbs%norb), stat=istat)
      call memocc(istat, fnrmOvrlpArr, 'fnrmOvrlpArr', subname)
      allocate(fnrmOldArr(tmb%orbs%norb), stat=istat)
      call memocc(istat, fnrmOldArr, 'fnrmOldArr', subname)
      allocate(alpha(tmb%orbs%norb), stat=istat)
      call memocc(istat, alpha, 'alpha', subname)
      allocate(lagMat(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
      call memocc(istat, lagMat, 'lagMat', subname)
    end subroutine allocateArrays


    subroutine deallocateArrays()
      iall=-product(shape(alpha))*kind(alpha)
      deallocate(alpha, stat=istat)
      call memocc(istat, iall, 'alpha', subname)

      iall=-product(shape(lagMat))*kind(lagMat)
      deallocate(lagMat, stat=istat)
      call memocc(istat, iall, 'lagMat', subname)
     
      iall=-product(shape(coeffPad))*kind(coeffPad)
      deallocate(coeffPad, stat=istat)
      call memocc(istat, iall, 'coeffPad', subname)

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
