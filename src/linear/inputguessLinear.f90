!> @file
!! Input guess wavefunctions for linear version
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Input wavefunctions are found by a diagonalization in a minimal basis set
!! Each processors write its initial wavefunctions into the wavefunction file
!! The files are then read by readwave
subroutine initInputguessConfinement(iproc, nproc, at, lzd, orbs, collcom_reference, &
           Glr, input, hx, hy, hz, lin, tmb, tmbgauss, rxyz, nscatterarr)

  use module_base
  use module_types
  use module_interfaces, exceptThisOne => initInputguessConfinement
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc,nproc
  real(gp), intent(in) :: hx, hy, hz
  type(atoms_data),intent(inout) :: at
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  type(collective_comms),intent(in) :: collcom_reference
  type(locreg_descriptors),intent(in) :: Glr
  type(input_variables), intent(in) :: input
  type(linearInputParameters),intent(in) :: lin
  type(DFT_wavefunction),intent(in) :: tmb
  type(DFT_wavefunction),intent(out) :: tmbgauss
  integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(gp),dimension(3,at%nat),intent(in) :: rxyz

  ! Local variables
  character(len=*), parameter :: subname='initInputguessConfinement'
  real(gp), dimension(:),allocatable :: locrad
  real(gp),dimension(:,:),allocatable :: locregCenter
  integer,dimension(:),allocatable :: norbsPerAt, norbsPerLocreg
  integer, parameter :: nmax=6,lmax=3,noccmax=2,nelecmax=32
  integer :: ist, iadd, ii, jj, norbtot, istat, iall, iat, nspin_ig, norbat, ityp, ilr, iorb, ndim

  ! maybe use kswfn_init_comm for initialization?

  nullify(tmbgauss%psi)
  nullify(tmbgauss%hpsi)
  nullify(tmbgauss%psit)
  nullify(tmbgauss%psit_c)
  nullify(tmbgauss%psit_f)
  nullify(tmbgauss%spsi)
  nullify(tmbgauss%gaucoeffs)
  tmbgauss%can_use_transposed=.false.

  ! Nullify the local zone descriptors.
  call nullify_local_zone_descriptors(tmbgauss%lzd)
  call nullify_orbitals_data(tmbgauss%orbs)

  ! Allocate some arrays we need for the input guess.
  !!allocate(locrad(at%nat+ndebug),stat=istat)
  !!call memocc(istat,locrad,'locrad',subname)
  allocate(norbsPerAt(at%nat), stat=istat)
  call memocc(istat, norbsPerAt, 'norbsPerAt', subname)

  tmbgauss%lzd%hgrids(:)=lzd%hgrids(:)

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

  tmbgauss%lzd%nlr=at%nat

  allocate(norbsPerLocreg(norbtot), stat=istat)
  call memocc(istat, norbsPerLocreg, 'norbsPerLocreg', subname)

  allocate(locrad(norbtot),stat=istat)
  call memocc(istat,locrad,'locrad',subname)

  norbsPerLocreg=1

  ! Nullify the orbitals_data type and then determine its values.
  allocate(locregCenter(3,norbtot), stat=istat)
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

  ! Nullify the locreg_descriptors and then copy Glr to it.
  call nullify_locreg_descriptors(tmbgauss%lzd%Glr)
  call copy_locreg_descriptors(Glr, tmbgauss%lzd%Glr, subname)

  ! Determine the localization regions for the atomic orbitals, which have a different localization radius.
  !locrad=max(12.d0,maxval(lin%locrad(:)))
  locrad=lin%locrad
  !locrad=max(1.d0,maxval(lin%locrad(:)))
  !call nullify_orbitals_data(tmbgauss%orbs)
  call copy_orbitals_data(tmb%orbs, tmbgauss%orbs, subname)

  ! lzdig%orbs%inWhichLocreg has been allocated in orbitals_descriptors. Since it will again be allcoated
  ! in assignToLocreg2, deallocate it first.
  !!iall=-product(shape(tmbgauss%orbs%inWhichLocreg))*kind(tmbgauss%orbs%inWhichLocreg)
  !!deallocate(tmbgauss%orbs%inWhichLocreg,stat=istat)
  !!call memocc(istat,iall,'tmbgauss%orbs%inWhichLocreg',subname)
  !!! Assign the orbitals to the localization regions.
  !!call assignToLocreg2(iproc, nproc, tmbgauss%orbs%norb, tmbgauss%orbs%norb_par, at%nat, tmbgauss%lzd%nlr, &
  !!     input%nspin, norbsPerAt, rxyz, tmbgauss%orbs%inwhichlocreg)

  ! Take inwhichlocreg from tmb (otherwise there might be problems after the restart...
  do iorb=1,tmbgauss%orbs%norb
      tmbgauss%orbs%inwhichlocreg(iorb)=tmb%orbs%onwhichatom(iorb)
  end do
  !!tmbgauss%orbs%onwhichatom=tmb%orbs%onwhichatom
  

  call initLocregs(iproc, nproc, tmbgauss%lzd%nlr, rxyz, input%hx, input%hy, input%hz, tmbgauss%lzd, &
       tmbgauss%orbs, Glr, locrad, 's')

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
subroutine inputguessConfinement(iproc, nproc, inputpsi, at, &
     input, hx, hy, hz, lzd, lorbs, rxyz, denspot, rhopotold,&
     nlpspd, proj, GPU, lphi,orbs,tmb, tmblarge,energs,overlapmatrix)
  ! Input wavefunctions are found by a diagonalization in a minimal basis set
  ! Each processors write its initial wavefunctions into the wavefunction file
  ! The files are then read by readwave
  use module_base
  use module_interfaces, exceptThisOne => inputguessConfinement
  use module_types
  use Poisson_Solver
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,inputpsi
  real(gp), intent(in) :: hx, hy, hz
  type(atoms_data), intent(inout) :: at
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(GPU_pointers), intent(inout) :: GPU
  type(DFT_local_fields), intent(inout) :: denspot
  type(input_variables),intent(in) :: input
  type(local_zone_descriptors),intent(inout) :: lzd
  type(orbitals_data),intent(in) :: lorbs
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
  real(dp),dimension(max(lzd%glr%d%n1i*lzd%glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin),intent(inout) ::  rhopotold
  real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(out) :: lphi
  type(orbitals_data),intent(inout) :: orbs
  type(DFT_wavefunction),intent(inout) :: tmb
  type(DFT_wavefunction),intent(inout) :: tmblarge
  type(energy_terms),intent(inout) :: energs
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(out):: overlapmatrix

  ! Local variables
  type(gaussian_basis) :: G !basis for davidson IG
  character(len=*), parameter :: subname='inputguessConfinement'
  integer :: istat,iall,iat,nspin_ig,iorb,nvirt,norbat,ilrl,ilrg
  real(gp) :: hxh,hyh,hzh,eks,fnrm,V3prb, x0
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(gp), dimension(:), allocatable :: locrad
  real(wp), dimension(:,:,:), pointer :: psigau
  real(8),dimension(:),allocatable :: lchi, lchi2
  real(8),dimension(:,:),allocatable::  lhchi, locregCenter!, density_kernel!, ovrlp
  real(8), dimension(:,:,:),allocatable :: ham
  integer, dimension(:),allocatable :: norbsPerAt, onWhichAtomTemp, mapping, inversemapping
  logical,dimension(:),allocatable :: covered
  integer, parameter :: nmax=6,lmax=3,noccmax=2,nelecmax=32
  logical :: isoverlap
  logical,dimension(:),allocatable :: skip
  integer :: ist,jorb,iadd,ii,jj,ilr,ind1,ind2,ityp,owa,owa_old,ii_old
  integer :: ldim,gdim,jlr,iiorb,ndim_lhchi
  integer :: nlocregPerMPI,jproc,jlrold,infoCoeff
  !!integer,dimension(:),allocatable :: norb_parTemp, onWhichMPITemp
  type(confpot_data), dimension(:), allocatable :: confdatarr
!!  real(dp),dimension(6) :: xcstr
  type(DFT_wavefunction) :: tmbgauss
  type(GPU_pointers) :: GPUe
  character(len=2) :: symbol
  real(kind=8) :: rcov,rprb,ehomo,amu                                          
  real(kind=8) :: neleconf(nmax,0:lmax)                                        
  integer :: nsccode,mxpl,mxchg

  real(8),dimension(:),allocatable :: locrad_tmp
  !!real(8),dimension(:,:),allocatable:: locregCenter



  ! Initialize everything
  call timing(iproc,'init_inguess  ','ON')
  call initInputguessConfinement(iproc, nproc, at, lzd, lorbs, tmb%collcom, lzd%glr, input, hx, hy, hz, input%lin, &
       tmb, tmbgauss, rxyz, denspot%dpbox%nscatterarr)
  call timing(iproc,'init_inguess  ','OF')

  ! Allocate some arrays we need for the input guess.
  allocate(norbsc_arr(at%natsc+1,input%nspin+ndebug),stat=istat)
  call memocc(istat,norbsc_arr,'norbsc_arr',subname)
  allocate(locrad(at%nat+ndebug),stat=istat)
  call memocc(istat,locrad,'locrad',subname)
  allocate(norbsPerAt(at%nat), stat=istat)
  call memocc(istat, norbsPerAt, 'norbsPerAt', subname)
  allocate(mapping(tmb%orbs%norb), stat=istat)
  call memocc(istat, mapping, 'mapping', subname)
  allocate(covered(tmb%orbs%norb), stat=istat)
  call memocc(istat, covered, 'covered', subname)
  allocate(inversemapping(tmb%orbs%norb), stat=istat)
  call memocc(istat, inversemapping, 'inversemapping', subname)

  GPUe = GPU

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
  tmbgauss%lzd%nlr=at%nat

  allocate(locregCenter(3,norbat), stat=istat)
  call memocc(istat, locregCenter, 'locregCenter', subname)

  ilr=0
  do iat=1,at%nat
      ityp=at%iatype(iat)
      do iorb=1,norbsPerAt(iat)
          ilr=ilr+1
          locregCenter(:,ilr)=rxyz(:,iat)
      end do
  end do


  ! This array gives a mapping from the 'natural' orbital distribution (i.e. simply counting up the atoms) to
  ! our optimized orbital distribution (determined by in orbs%inwhichlocreg).
  iiorb=0
  covered=.false.
  do iat=1,at%nat
      do iorb=1,norbsPerAt(iat)
          iiorb=iiorb+1
          ! Search the corresponding entry in inwhichlocreg
          do jorb=1,tmb%orbs%norb
              if(covered(jorb)) cycle
              jlr=tmb%orbs%inwhichlocreg(jorb)
              if( tmb%lzd%llr(jlr)%locregCenter(1)==rxyz(1,iat) .and. &
                  tmb%lzd%llr(jlr)%locregCenter(2)==rxyz(2,iat) .and. &
                  tmb%lzd%llr(jlr)%locregCenter(3)==rxyz(3,iat) ) then
                  covered(jorb)=.true.
                  mapping(iiorb)=jorb
                  exit
              end if
          end do
      end do
  end do

  ! Inverse mapping
  do iorb=1,tmb%orbs%norb
      do jorb=1,tmb%orbs%norb
          if(mapping(jorb)==iorb) then
              inversemapping(iorb)=jorb
              exit
          end if
      end do
  end do



  nvirt=0
  call deallocate_orbitals_data(tmbgauss%orbs,subname)

  do ityp=1,at%ntypes
     call eleconf(at%nzatom(ityp),at%nelpsp(ityp),symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg,amu)
     if(4.d0*rprb>input%lin%locrad_type(ityp)) then
         if(iproc==0) write(*,'(3a,es10.2)') 'WARNING: locrad for atom type ',trim(symbol), &
                      ' is too small; minimal value is ',4.d0*rprb
     end if
     if(input%lin%potentialPrefac_lowaccuracy(ityp)>0.d0) then
         x0=(70.d0/input%lin%potentialPrefac_lowaccuracy(ityp))**.25d0
         if(iproc==0) write(*,'(a,a,2es11.2,es12.3)') 'type, 4.d0*rprb, x0, input%lin%locrad_type(ityp)', &
                      trim(symbol),4.d0*rprb, x0, input%lin%locrad_type(ityp)
         V3prb=input%lin%potentialPrefac_lowaccuracy(ityp)*(4.d0*rprb)**4
         if(iproc==0) write(*,'(a,es14.4)') 'V3prb',V3prb
     end if
  end do


  call inputguess_gaussian_orbitals_forLinear(iproc,nproc,tmb%orbs%norb,at,rxyz,nvirt,nspin_ig,&
       at%nat, norbsPerAt, mapping, &
       lorbs,tmbgauss%orbs,norbsc_arr,locrad,G,psigau,eks,input%lin%potentialPrefac_lowaccuracy)
  ! Take inwhichlocreg from tmb (otherwise there might be problems after the restart...
  do iorb=1,tmb%orbs%norb
      tmbgauss%orbs%inwhichlocreg(iorb)=tmb%orbs%onwhichatom(iorb)
  end do


  !dimension of the wavefunctions
  call wavefunction_dimension(tmbgauss%lzd,tmbgauss%orbs)


  ! Allcoate the array holding the orbitals. lchi2 are the atomic orbitals with the larger cutoff, whereas
  ! lchi are the atomic orbitals with the smaller cutoff.
  allocate(lchi2(max(tmbgauss%orbs%npsidim_orbs,tmbgauss%orbs%npsidim_comp)),stat=istat)
  call memocc(istat,lchi2,'lchi2',subname)
  !!lchi2=0.d0
  call to_zero(max(tmbgauss%orbs%npsidim_orbs,tmbgauss%orbs%npsidim_comp), lchi2(1))

  ! Grid spacing on fine grid.
  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz

  ! Assign the size of the orbitals to the new variable lpsidimtot.
  tmbgauss%lzd%hgrids(1)=hx
  tmbgauss%lzd%hgrids(2)=hy
  tmbgauss%lzd%hgrids(3)=hz
  ! Transform the atomic orbitals to the wavelet basis.
  call gaussians_to_wavelets_new(iproc,nproc,tmbgauss%lzd,tmbgauss%orbs,G,&
       psigau(1,1,min(tmb%orbs%isorb+1,tmb%orbs%norb)),lchi2)

  iall=-product(shape(psigau))*kind(psigau)
  deallocate(psigau,stat=istat)
  call memocc(istat,iall,'psigau',subname)

  call deallocate_gwf(G,subname)

  !restore wavefunction dimension



  call to_zero(max(lorbs%npsidim_orbs,lorbs%npsidim_comp), lphi(1))

  ! Transform chi to the localization region. This requires that the localizatin region of lchi2 is larger than that
  ! of lchi.
  ind1=1
  ind2=1
  do iorb=1,tmb%orbs%norbp
      ilrl = tmb%orbs%inWhichLocreg(tmb%orbs%isorb+iorb)
      !ilrg = tmbgauss%orbs%inWhichLocreg(tmbgauss%orbs%isorb+iorb)
      ilrg = tmb%orbs%onwhichatom(tmb%orbs%isorb+iorb)
      ldim=tmb%lzd%Llr(ilrl)%wfd%nvctr_c+7*tmb%lzd%Llr(ilrl)%wfd%nvctr_f
      gdim=tmbgauss%lzd%llr(ilrg)%wfd%nvctr_c+7*tmbgauss%lzd%llr(ilrg)%wfd%nvctr_f
      if (tmb%lzd%llr(ilrl)%locregCenter(1) /= tmbgauss%lzd%llr(ilrg)%locregCenter(1) .or. &
          tmb%lzd%llr(ilrl)%locregCenter(1) /= tmbgauss%lzd%llr(ilrg)%locregCenter(1) .or. &
          tmb%lzd%llr(ilrl)%locregCenter(1) /= tmbgauss%lzd%llr(ilrg)%locregCenter(1)) then
          stop 'ERROR: locreg centers are not identical!'
      end if
      call psi_to_locreg2(iproc, nproc, ldim, gdim, tmb%lzd%llr(ilrl), tmbgauss%lzd%llr(ilrg), lchi2(ind1), lphi(ind2))
      ind1=ind1+gdim
      ind2=ind2+ldim
  end do

  if(tmbgauss%orbs%norbp>0 .and. ind1/=tmbgauss%orbs%npsidim_orbs+1) then
      write(*,'(2(a,i8),i8)') 'ERROR on process ',iproc,&
           ': ind1/=tmbgauss%orbs%npsidim+1',ind1,tmbgauss%orbs%npsidim_orbs+1
      stop
  end if
  if(tmb%orbs%norbp>0 .and. ind2/=tmb%orbs%npsidim_orbs+1) then
      write(*,'(2(a,i8),i8)') 'ERROR on process ',iproc,&
           ': ind2/=tmb%orbs%npsidim+1',ind2,tmb%orbs%npsidim_orbs+1
      stop
  end if


  ! Deallocate locrad, which is not used any longer.
  iall=-product(shape(locrad))*kind(locrad)
  deallocate(locrad,stat=istat)
  call memocc(istat,iall,'locrad',subname)

  !change again wavefunction dimension
  call wavefunction_dimension(tmbgauss%lzd,tmbgauss%orbs)




  ! Create the potential. First calculate the charge density.
  !!call sumrho(denspot%dpbox,tmbgauss%orbs,tmbgauss%lzd,GPUe,at%sym,denspot%rhod,&
  !!     lchi2,denspot%rho_psi,inversemapping)
  do iorb=1,tmb%orbs%norb
      tmb%orbs%occup(iorb)=tmbgauss%orbs%occup(iorb)
  end do
  call sumrho(denspot%dpbox,tmb%orbs,tmb%lzd,GPUe,at%sym,denspot%rhod,&
       lphi,denspot%rho_psi,inversemapping)
  call communicate_density(denspot%dpbox,input%nspin,&!hxh,hyh,hzh,tmbgauss%lzd,&
       denspot%rhod,denspot%rho_psi,denspot%rhov,.false.)


  if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE) then
      call dcopy(max(lzd%glr%d%n1i*lzd%glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
  end if


  iall=-product(shape(lchi2))*kind(lchi2)
  deallocate(lchi2, stat=istat)
  call memocc(istat, iall, 'lchi2',subname)

  call deallocate_local_zone_descriptors(tmbgauss%lzd, subname)

  call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)

  if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
      call dcopy(max(lzd%glr%d%n1i*lzd%glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
  end if

  if (input%exctxpar == 'OP2P') energs%eexctX = uninitialized(energs%eexctX)


  call get_coeff(iproc,nproc,LINEAR_MIXDENS_SIMPLE,lzd,orbs,at,rxyz,denspot,GPU,infoCoeff,energs%ebs,nlpspd,proj,&
       input%SIC,tmb,fnrm,overlapmatrix,.true.,.false.,&
       tmblarge)

  ! Important: Don't use for the rest of the code
  tmblarge%can_use_transposed = .false.

  if(associated(tmblarge%psit_c)) then
      iall=-product(shape(tmblarge%psit_c))*kind(tmblarge%psit_c)
      deallocate(tmblarge%psit_c, stat=istat)
      call memocc(istat, iall, 'tmblarge%psit_c', subname)
  end if
  if(associated(tmblarge%psit_f)) then
      iall=-product(shape(tmblarge%psit_f))*kind(tmblarge%psit_f)
      deallocate(tmblarge%psit_f, stat=istat)
      call memocc(istat, iall, 'tmblarge%psit_f', subname)
  end if
  

  if(iproc==0) write(*,'(1x,a)') '------------------------------------------------------------- Input guess generated.'
  
  ! Deallocate all local arrays.

  ! Deallocate all types that are not needed any longer.
  call deallocate_orbitals_data(tmbgauss%orbs, subname)

  ! Deallocate all remaining local arrays.
  iall=-product(shape(norbsc_arr))*kind(norbsc_arr)
  deallocate(norbsc_arr,stat=istat)
  call memocc(istat,iall,'norbsc_arr',subname)

  iall=-product(shape(norbsPerAt))*kind(norbsPerAt)
  deallocate(norbsPerAt, stat=istat)
  call memocc(istat, iall, 'norbsPerAt',subname)

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





