!> @file 
!!   Input guess for the linear version
!! @author
!!   Copyright (C) 2011-2012 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Input guess wavefunction diagonalization
!! Input wavefunctions are found by a diagonalization in a minimal basis set
!! Each processors write its initial wavefunctions into the wavefunction file
!! The files are then read by readwave
subroutine inputguessConfinement(iproc, nproc, at, input, hx, hy, hz, &
     rxyz, nlpsp, GPU, orbs, kswfn, tmb, denspot, rhopotold, energs, &
     locregcenters)
  use module_base
  use module_interfaces, exceptThisOne => inputguessConfinement
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use yaml_output
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc
  real(gp), intent(in) :: hx, hy, hz
  type(atoms_data), intent(inout) :: at
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  type(GPU_pointers), intent(inout) :: GPU
  type(input_variables),intent(in) :: input
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  type(orbitals_data),intent(inout) :: orbs
  type(DFT_wavefunction),intent(inout) :: kswfn, tmb
  type(DFT_local_fields), intent(inout) :: denspot
  real(dp), dimension(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin),intent(inout) :: rhopotold
  type(energy_terms),intent(inout) :: energs
  real(kind=8),dimension(3,at%astruct%nat),intent(in),optional :: locregcenters

  ! Local variables
  type(gaussian_basis) :: G !basis for davidson IG
  character(len=*), parameter :: subname='inputguessConfinement'
  integer :: istat,iall,iat,nspin_ig,iorb,nvirt,norbat,matrixindex_in_compressed
  real(gp) :: hxh,hyh,hzh,eks,fnrm,V3prb,x0
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(gp), dimension(:), allocatable :: locrad
  real(wp), dimension(:,:,:), pointer :: psigau
  integer, dimension(:),allocatable :: norbsPerAt, mapping, inversemapping, minorbs_type, maxorbs_type
  logical,dimension(:),allocatable :: covered, type_covered
  real(kind=8),dimension(:,:),allocatable :: aocc
  integer :: ist,jorb,iadd,ii,jj,ityp,itype,iortho
  integer :: jlr,iiorb
  integer :: infoCoeff
  type(orbitals_data) :: orbs_gauss
  type(GPU_pointers) :: GPUe
  character(len=2) :: symbol
  real(kind=8) :: rcov,rprb,ehomo,amu,pnrm
  integer :: nsccode,mxpl,mxchg
  type(mixrhopotDIISParameters) :: mixdiis
  type(sparseMatrix) :: ham_small ! for FOE
  logical :: finished
  type(confpot_data),dimension(:),allocatable :: confdatarrtmp
  real(kind=8),dimension(:),allocatable :: philarge
  integer :: npsidim_large, sdim, ldim, ists, istl, ilr, nspin, info_basis_functions
  real(kind=8) :: ratio_deltas, trace, trace_old, fnrm_tmb
  logical :: ortho_on, reduce_conf
  type(localizedDIISParameters) :: ldiis
  real(wp), dimension(:,:,:), pointer :: mom_vec_fake

  call nullify_orbitals_data(orbs_gauss)
  nullify(mom_vec_fake)

  ! Allocate some arrays we need for the input guess.
  allocate(norbsc_arr(at%natsc+1,input%nspin+ndebug),stat=istat)
  call memocc(istat,norbsc_arr,'norbsc_arr',subname)
  allocate(locrad(at%astruct%nat+ndebug),stat=istat)
  call memocc(istat,locrad,'locrad',subname)
  allocate(norbsPerAt(at%astruct%nat), stat=istat)
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

  ! Keep the natural occupations
  allocate(aocc(32,at%astruct%nat),stat=istat)
  call memocc(istat,aocc,'aocc',subname)
  call vcopy(32*at%astruct%nat, at%aocc(1,1), 1, aocc(1,1), 1)

  ! Determine how many atomic orbitals we have. Maybe we have to increase this number to more than
  ! its 'natural' value.
  norbat=0
  ist=0
  do iat=1,at%astruct%nat
      ii=input%lin%norbsPerType(at%astruct%iatype(iat))
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

  ! This array gives a mapping from the 'natural' orbital distribution (i.e. simply counting up the atoms) to
  ! our optimized orbital distribution (determined by in orbs%inwhichlocreg).
  iiorb=0
  covered=.false.
  if (present(locregcenters)) then
      do iat=1,at%astruct%nat
          do iorb=1,norbsPerAt(iat)
              iiorb=iiorb+1
              ! Search the corresponding entry in inwhichlocreg
              do jorb=1,tmb%orbs%norb
                  if(covered(jorb)) cycle
                  jlr=tmb%orbs%inwhichlocreg(jorb)
                  if( tmb%lzd%llr(jlr)%locregCenter(1)==locregcenters(1,iat) .and. &
                      tmb%lzd%llr(jlr)%locregCenter(2)==locregcenters(2,iat) .and. &
                      tmb%lzd%llr(jlr)%locregCenter(3)==locregcenters(3,iat) ) then
                      covered(jorb)=.true.
                      mapping(iiorb)=jorb
                      exit
                  end if
              end do
          end do
      end do
  else
      do iat=1,at%astruct%nat
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
  end if

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

  !!do iorb=1,tmb%orbs%norb
  !!    ilr=tmb%orbs%inwhichlocreg(iorb)
  !!    write(500+10*iproc+0,*) tmb%lzd%llr(ilr)%locregcenter(1:3)
  !!    write(500+10*iproc+1,*) tmb%ham_descr%lzd%llr(ilr)%locregcenter(1:3)
  !!end do

! THIS OUTPUT SHOULD PROBABLY BE KEPT, BUT IS COMMENTED FOR THE MOMENT AS IT DOES NOT
! SEEM TO BE RELEVANT ANY MORE
  !!do ityp=1,at%astruct%ntypes
  !!   call eleconf(at%nzatom(ityp),at%nelpsp(ityp),symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg,amu)
  !!   if(4.d0*rprb>input%lin%locrad_type(ityp)) then
  !!       if(iproc==0) write(*,'(3a,es10.2)') 'WARNING: locrad for atom type ',trim(symbol), &
  !!                    ' is too small; minimal value is ',4.d0*rprb
  !!   end if
  !!   if(input%lin%potentialPrefac_ao(ityp)>0.d0) then
  !!       x0=(70.d0/input%lin%potentialPrefac_ao(ityp))**.25d0
  !!       if(iproc==0) write(*,'(a,a,2es11.2,es12.3)') 'type, 4.d0*rprb, x0, input%lin%locrad_type(ityp)', &
  !!                    trim(symbol),4.d0*rprb, x0, input%lin%locrad_type(ityp)
  !!       V3prb=input%lin%potentialPrefac_ao(ityp)*(4.d0*rprb)**4
  !!       if(iproc==0) write(*,'(a,es14.4)') 'V3prb',V3prb
  !!   end if
  !!end do

! THIS IS SOMETHING EXPERIMENTAL TO ESTIMATE THE CONVERGENCE THRESHOLD. NOT
! WORKING WELL, BUT STILL TO BE KEPT AS A TEMPLATE
!!  ! #######################################################################
!!  ! Estimate convergence criterion: kinetic energy for Gaussians and for
!!  ! wavelets (i.e. with cutoff)
!!  call inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,nvirt,nspin_ig,&
!!       tmb%orbs,orbs_gauss,norbsc_arr,locrad,G,psigau,eks,mapping)!,1.d-7*input%lin%potentialPrefac_ao)
!!  if (iproc==0) write(*,*) 'eks',eks
!!
!!  ! Create the potential. First calculate the charge density.
!!  do iorb=1,tmb%orbs%norb
!!      !if (iproc==0) write(*,*) 'WARNING: use mapping for occupation numbers!'
!!      !tmb%orbs%occup(iorb)=orbs_gauss%occup(iorb)
!!      tmb%orbs%occup(iorb)=orbs_gauss%occup(inversemapping(iorb))
!!  end do
!!
!!  ! Transform the atomic orbitals to the wavelet basis.
!!
!!  !!if (.false.) then
!!      ! linear version
!!
!!      if (orbs_gauss%norb/=tmb%orbs%norb) stop 'orbs%gauss%norb does not match tmbs%orbs%norb'
!!      orbs_gauss%inwhichlocreg=tmb%orbs%inwhichlocreg
!!      call wavefunction_dimension(tmb%lzd,orbs_gauss)
!!      call to_zero(max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%psi(1))
!!      call gaussians_to_wavelets_new(iproc,nproc,tmb%lzd,orbs_gauss,G,&
!!           psigau(1,1,min(tmb%orbs%isorb+1,tmb%orbs%norb)),tmb%psi)
!!
!!      ! Calculate kinetic energy
!!      allocate(confdatarrtmp(tmb%orbs%norbp))
!!      call default_confinement_data(confdatarrtmp,tmb%orbs%norbp)
!!
!!      call small_to_large_locreg(iproc, tmb%npsidim_orbs, &
!!           tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
!!           tmb%orbs, tmb%psi, tmb%ham_descr%psi)
!!      if (tmb%ham_descr%npsidim_orbs > 0) call to_zero(tmb%ham_descr%npsidim_orbs,tmb%hpsi(1))
!!
!!      call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
!!           tmb%ham_descr%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%ham_descr%psi,tmb%hpsi,&
!!           energs,input%SIC,GPU,3,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
!!      call SynchronizeHamiltonianApplication(nproc,tmb%ham_descr%npsidim_orbs,tmb%orbs,tmb%ham_descr%lzd,GPU,tmb%hpsi,&
!!           energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
!!
!!  !!else
!!  !!    ! cubic version
!!
!!  !!    if (orbs_gauss%norb/=tmb%orbs%norb) stop 'orbs%gauss%norb does not match tmbs%orbs%norb'
!!  !!    orbs_gauss%inwhichlocreg=tmb%orbs%inwhichlocreg
!!  !!    call wavefunction_dimension(tmb%lzd,orbs_gauss)
!!  !!    call to_zero(max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%psi(1))
!!  !!    call gaussians_to_wavelets_new(iproc,nproc,tmb%lzd,orbs_gauss,G,&
!!  !!         psigau(1,1,min(tmb%orbs%isorb+1,tmb%orbs%norb)),tmb%psi)
!!
!!  !!    ! Calculate kinetic energy
!!  !!    allocate(confdatarrtmp(tmb%orbs%norbp))
!!  !!    call default_confinement_data(confdatarrtmp,tmb%orbs%norbp)
!!
!!  !!    call small_to_large_locreg(iproc, tmb%npsidim_orbs, &
!!  !!         tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
!!  !!         tmb%orbs, tmb%psi, tmb%ham_descr%psi)
!!  !!    if (tmb%ham_descr%npsidim_orbs > 0) call to_zero(tmb%ham_descr%npsidim_orbs,tmb%hpsi(1))
!!
!!  !!    call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
!!  !!         tmb%ham_descr%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%ham_descr%psi,tmb%hpsi,&
!!  !!         energs,input%SIC,GPU,3,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
!!  !!    call SynchronizeHamiltonianApplication(nproc,tmb%ham_descr%npsidim_orbs,tmb%orbs,tmb%ham_descr%lzd,GPU,tmb%hpsi,&
!!  !!         energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
!!
!!  !!end if
!!
!!  if (iproc==0) write(*,*) 'eks, energs%ekin', eks, energs%ekin
!!  if (iproc==0) write(*,*) 'conv crit:', abs(eks-energs%ekin)/dble(tmb%orbs%norb)
!!  deallocate(confdatarrtmp)
!!  iall=-product(shape(psigau))*kind(psigau)
!!  deallocate(psigau,stat=istat)
!!  call memocc(istat,iall,'psigau',subname)
!!
!!  iall=-product(shape(orbs_gauss%onwhichatom))*kind(orbs_gauss%onwhichatom)
!!  deallocate(orbs_gauss%onwhichatom,stat=istat)
!!  call memocc(istat,iall,'orbs_gauss%onwhichatom',subname)
!!
!!  iall=-product(shape(orbs_gauss%norb_par))*kind(orbs_gauss%norb_par)
!!  deallocate(orbs_gauss%norb_par,stat=istat)
!!  call memocc(istat,iall,'orbs_gauss%norb_par',subname)
!!
!!  iall=-product(shape(orbs_gauss%kpts))*kind(orbs_gauss%kpts)
!!  deallocate(orbs_gauss%kpts,stat=istat)
!!  call memocc(istat,iall,'psigau',subname)
!!
!!  iall=-product(shape(orbs_gauss%spinsgn))*kind(orbs_gauss%spinsgn)
!!  deallocate(orbs_gauss%spinsgn,stat=istat)
!!  call memocc(istat,iall,'orbs_gauss%spinsgn',subname)
!!
!!  iall=-product(shape(orbs_gauss%ikptproc))*kind(orbs_gauss%ikptproc)
!!  deallocate(orbs_gauss%ikptproc,stat=istat)
!!  call memocc(istat,iall,'orbs_gauss%ikptproc',subname)
!!
!!  iall=-product(shape(orbs_gauss%kwgts))*kind(orbs_gauss%kwgts)
!!  deallocate(orbs_gauss%kwgts,stat=istat)
!!  call memocc(istat,iall,'orbs_gauss%kwgts',subname)
!!
!!  iall=-product(shape(orbs_gauss%occup))*kind(orbs_gauss%occup)
!!  deallocate(orbs_gauss%occup,stat=istat)
!!  call memocc(istat,iall,'orbs_gauss%occup',subname)
!!
!!  iall=-product(shape(orbs_gauss%inwhichlocreg))*kind(orbs_gauss%inwhichlocreg)
!!  deallocate(orbs_gauss%inwhichlocreg,stat=istat)
!!  call memocc(istat,iall,'orbs_gauss%inwhichlocreg',subname)
!!
!!  iall=-product(shape(orbs_gauss%iokpt))*kind(orbs_gauss%iokpt)
!!  deallocate(orbs_gauss%iokpt,stat=istat)
!!  call memocc(istat,iall,'orbs_gauss%iokpt',subname)
!!
!!  iall=-product(shape(orbs_gauss%ispot))*kind(orbs_gauss%ispot)
!!  deallocate(orbs_gauss%ispot,stat=istat)
!!  call memocc(istat,iall,'orbs_gauss%ispot',subname)
!!
!!  iall=-product(shape(orbs_gauss%isorb_par))*kind(orbs_gauss%isorb_par)
!!  deallocate(orbs_gauss%isorb_par,stat=istat)
!!  call memocc(istat,iall,'orbs_gauss%isorb_par',subname)
!!
!!  iall=-product(shape(G%ndoc))*kind(G%ndoc)
!!  deallocate(G%ndoc,stat=istat)
!!  call memocc(istat,iall,'G%ndoc',subname)
!!
!!  iall=-product(shape(G%nshell))*kind(G%nshell)
!!  deallocate(G%nshell,stat=istat)
!!  call memocc(istat,iall,'G%xp',subname)
!!
!!  iall=-product(shape(G%xp))*kind(G%xp)
!!  deallocate(G%xp,stat=istat)
!!  call memocc(istat,iall,'G%xp',subname)
!!
!!  iall=-product(shape(G%psiat))*kind(G%psiat)
!!  deallocate(G%psiat,stat=istat)
!!  call memocc(istat,iall,'G%psiat',subname)
!!
!!  iall=-product(shape(G%nam))*kind(G%nam)
!!  deallocate(G%nam,stat=istat)
!!  call memocc(istat,iall,'G%nam',subname)
!!
!!
!!
!!  !!call f_free(tmb%orbs%onwhichatom)
!!  !!call f_free(tmb%orbs%norb_par)
!!  !!call f_free(tmb%orbs%kpts)
!!  !!call f_free(tmb%orbs%spinsgn)
!!  !!call f_free(tmb%orbs%ikptproc)
!!  !!call f_free(tmb%orbs%kwgts)
!!  !!call f_free(tmb%orbs%occup)
!!  !!call f_free(tmb%orbs%inwhichlocreg)
!!  !!call f_free(tmb%orbs%iokpt)
!!  !!call f_free(tmb%orbs%ispot)
!!  !!call f_free(tmb%orbs%isorb_par)
!!  !!call f_free(G%ndoc)
!!  !!call f_free(G%nshell)
!!  !!call f_free(G%xp)
!!  !!call f_free(G%psiat)
!!
!!
!!
!!
!!  ! #######################################################################




  call inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,nvirt,nspin_ig,&
       tmb%orbs,orbs_gauss,norbsc_arr,locrad,G,psigau,eks,2,mapping,input%lin%potentialPrefac_ao)
  !!call inputguess_gaussian_orbitals_forLinear(iproc,nproc,tmb%orbs%norb,at,rxyz,nvirt,nspin_ig,&
  !!     tmb%lzd%nlr,norbsPerAt,mapping, &
  !!     tmb%orbs,orbs_gauss,norbsc_arr,locrad,G,psigau,eks,input%lin%potentialPrefac_ao)

  ! Take inwhichlocreg from tmb (otherwise there might be problems after the restart...
  !do iorb=1,tmb%orbs%norb
  !    orbs_gauss%inwhichlocreg(iorb)=tmb%orbs%onwhichatom(iorb)
  !end do


  ! Grid spacing on fine grid.
  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz

  ! Transform the atomic orbitals to the wavelet basis.
  if (orbs_gauss%norb/=tmb%orbs%norb) then
     print*,'orbs_gauss%norb does not match tmbs%orbs%norb',orbs_gauss%norb,tmb%orbs%norb
     stop 
  end if
  orbs_gauss%inwhichlocreg=tmb%orbs%inwhichlocreg
  call wavefunction_dimension(tmb%lzd,orbs_gauss)
  call to_zero(max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%psi(1))
  call gaussians_to_wavelets_new(iproc,nproc,tmb%lzd,orbs_gauss,G,&
       psigau(1,1,min(tmb%orbs%isorb+1,tmb%orbs%norb)),tmb%psi)

  iall=-product(shape(psigau))*kind(psigau)
  deallocate(psigau,stat=istat)
  call memocc(istat,iall,'psigau',subname)

  call deallocate_gwf(G,subname)

  ! Deallocate locrad, which is not used any longer.
  iall=-product(shape(locrad))*kind(locrad)
  deallocate(locrad,stat=istat)
  call memocc(istat,iall,'locrad',subname)

  ! Create the potential. First calculate the charge density.
  do iorb=1,tmb%orbs%norb
      !if (iproc==0 .and. iorb==1) write(*,*) 'WARNING: use mapping for occupation numbers!'
      !tmb%orbs%occup(iorb)=orbs_gauss%occup(iorb)
      tmb%orbs%occup(iorb)=orbs_gauss%occup(inversemapping(iorb))
  end do

  !!call sumrho(denspot%dpbox,tmb%orbs,tmb%lzd,GPUe,at%sym,denspot%rhod,&
  !!     tmb%psi,denspot%rho_psi,inversemapping)
  !!call communicate_density(denspot%dpbox,input%nspin,&!hxh,hyh,hzh,tmbgauss%lzd,&
  !!     denspot%rhod,denspot%rho_psi,denspot%rhov,.false.)

  !Put the Density kernel to identity for now
  call to_zero(tmb%linmat%denskern%nvctr, tmb%linmat%denskern%matrix_compr(1))
  do iorb=1,tmb%orbs%norb
     ii=matrixindex_in_compressed(tmb%linmat%denskern,iorb,iorb)
     !tmb%linmat%denskern%matrix_compr(ii)=1.d0*tmb%orbs%occup(inversemapping(iorb))
     tmb%linmat%denskern%matrix_compr(ii)=1.d0*tmb%orbs%occup(iorb)
  end do

  !Calculate the density in the new scheme
  call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, max(tmb%npsidim_orbs,tmb%npsidim_comp), &
       tmb%orbs, tmb%psi, tmb%collcom_sr)
  call sumrho_for_TMBs(iproc, nproc, tmb%Lzd%hgrids(1), tmb%Lzd%hgrids(2), tmb%Lzd%hgrids(3), &
       tmb%collcom_sr, tmb%linmat%denskern, tmb%Lzd%Glr%d%n1i*tmb%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)

  !!do istat=1,size(denspot%rhov)
  !!    write(300+iproc,*) istat, denspot%rhov(istat)
  !!end do 
  !!call mpi_finalize(istat)
  !!stop


  if (input%lin%mixing_after_inputguess) then
      if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE .or. input%lin%scf_mode==LINEAR_FOE &
           .or. input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
          call dcopy(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
      end if
  end if
  call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)

  if (input%lin%mixing_after_inputguess) then
      if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
          call dcopy(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
      end if
  end if
  if (input%exctxpar == 'OP2P') energs%eexctX = uninitialized(energs%eexctX)



  !!! PLOT ###########################################################################
  !!    hxh=0.5d0*tmb%lzd%hgrids(1)
  !!    hyh=0.5d0*tmb%lzd%hgrids(2)
  !!    hzh=0.5d0*tmb%lzd%hgrids(3)
  !!    npsidim_large=tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
  !!    allocate(philarge((tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f)*tmb%orbs%norbp))
  !!    philarge=0.d0
  !!    ists=1
  !!    istl=1
  !!    do iorb=1,tmb%orbs%norbp
  !!        ilr = tmb%orbs%inWhichLocreg(tmb%orbs%isorb+iorb)
  !!        sdim=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
  !!        ldim=tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
  !!        nspin=1 !this must be modified later
  !!        call Lpsi_to_global2(iproc, sdim, ldim, tmb%orbs%norb, tmb%orbs%nspinor, nspin, tmb%lzd%glr, &
  !!             tmb%lzd%llr(ilr), tmb%psi(ists), philarge(istl))
  !!        ists=ists+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
  !!        istl=istl+tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
  !!    end do
  !!    call plotOrbitals(iproc, tmb, philarge, at%astruct%nat, rxyz, hxh, hyh, hzh, 1, 'orbs')
  !!    deallocate(philarge)
  !!! END PLOT #######################################################################




  if (.not. input%lin%iterative_orthogonalization) then
      ! Standard orthonomalization
      !!if(iproc==0) write(*,*) 'calling orthonormalizeLocalized (exact)'
      if (iproc==0) call yaml_map('orthonormalization of input guess','standard')
      ! CHEATING here and passing tmb%linmat%denskern instead of tmb%linmat%inv_ovrlp
      !write(*,'(a,i4,4i8)') 'IG: iproc, lbound, ubound, minval, maxval',&
      !iproc, lbound(tmb%linmat%inv_ovrlp%matrixindex_in_compressed_fortransposed,2),&
      !ubound(tmb%linmat%inv_ovrlp%matrixindex_in_compressed_fortransposed,2),&
      !minval(tmb%collcom%indexrecvorbital_c),maxval(tmb%collcom%indexrecvorbital_c)
      !!if (iproc==0) write(*,*) 'WARNING: no ortho in inguess'
      call orthonormalizeLocalized(iproc, nproc, -1, tmb%npsidim_orbs, tmb%orbs, tmb%lzd, tmb%linmat%ovrlp, tmb%linmat%inv_ovrlp, &
           tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)
            
 else

     ! Iterative orthonomalization
     !!if(iproc==0) write(*,*) 'calling generalized orthonormalization'
     if (iproc==0) call yaml_map('orthonormalization of input guess','generalized')
     allocate(maxorbs_type(at%astruct%ntypes),stat=istat)
     call memocc(istat,maxorbs_type,'maxorbs_type',subname)
     allocate(minorbs_type(at%astruct%ntypes),stat=istat)
     call memocc(istat,minorbs_type,'minorbs_type',subname)
     allocate(type_covered(at%astruct%ntypes),stat=istat)
     call memocc(istat,type_covered,'type_covered',subname)
     minorbs_type(1:at%astruct%ntypes)=0
     iortho=0
     ortho_loop: do
         finished=.true.
         type_covered=.false.
         do iat=1,at%astruct%nat
             itype=at%astruct%iatype(iat)
             if (type_covered(itype)) cycle
             type_covered(itype)=.true.
             jj=1*ceiling(aocc(1,iat))+3*ceiling(aocc(3,iat))+&
                  5*ceiling(aocc(7,iat))+7*ceiling(aocc(13,iat))
             maxorbs_type(itype)=jj
             if (jj<input%lin%norbsPerType(at%astruct%iatype(iat))) then
                 finished=.false.
                 if (ceiling(aocc(1,iat))==0) then
                     aocc(1,iat)=1.d0
                 else if (ceiling(aocc(3,iat))==0) then
                     aocc(3,iat)=1.d0
                 else if (ceiling(aocc(7,iat))==0) then
                     aocc(7,iat)=1.d0
                 else if (ceiling(aocc(13,iat))==0) then
                     aocc(13,iat)=1.d0
                 end if
             end if
         end do
         if (iortho>0) then
             call gramschmidt_subset(iproc, nproc, -1, tmb%npsidim_orbs, &                                  
                  tmb%orbs, at, minorbs_type, maxorbs_type, tmb%lzd, tmb%linmat%ovrlp, &
                  tmb%linmat%inv_ovrlp, tmb%collcom, tmb%orthpar, &
                  tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)
         end if
         call orthonormalize_subset(iproc, nproc, -1, tmb%npsidim_orbs, &                                  
              tmb%orbs, at, minorbs_type, maxorbs_type, tmb%lzd, tmb%linmat%ovrlp, &
              tmb%linmat%inv_ovrlp, tmb%collcom, tmb%orthpar, &
              tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)
         if (finished) exit ortho_loop
         iortho=iortho+1
         minorbs_type(1:at%astruct%ntypes)=maxorbs_type(1:at%astruct%ntypes)+1
     end do ortho_loop
     iall=-product(shape(maxorbs_type))*kind(maxorbs_type)
     deallocate(maxorbs_type,stat=istat)
     call memocc(istat, iall,'maxorbs_type',subname)
     iall=-product(shape(minorbs_type))*kind(minorbs_type)
     deallocate(minorbs_type,stat=istat)
     call memocc(istat, iall,'minorbs_type',subname)
     iall=-product(shape(type_covered))*kind(type_covered)
     deallocate(type_covered,stat=istat)
     call memocc(istat, iall,'type_covered',subname)

 end if


 iall=-product(shape(aocc))*kind(aocc)
 deallocate(aocc,stat=istat)
 call memocc(istat, iall,'aocc',subname)

 !!call orthonormalizeLocalized(iproc, nproc, -1, tmb%npsidim_orbs, tmb%orbs, tmb%lzd, tmb%linmat%ovrlp, tmb%linmat%inv_ovrlp, &
 !!     tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)
 !!call mpi_finalize(istat)
 !!stop




 if (input%experimental_mode) then
     ! NEW: TRACE MINIMIZATION WITH ORTHONORMALIZATION ####################################
     ortho_on=.true.
     call initializeDIIS(input%lin%DIIS_hist_lowaccur, tmb%lzd, tmb%orbs, ldiis)
     ldiis%alphaSD=input%lin%alphaSD
     ldiis%alphaDIIS=input%lin%alphaDIIS
     energs%eexctX=0.d0 !temporary fix
     trace_old=0.d0 !initialization
     if (iproc==0) then
         !call yaml_close_map()
         call yaml_comment('Extended input guess for experimental mode',hfill='-')
         call yaml_open_map('Extended input guess')
         call yaml_open_sequence('support function optimization',label=&
                                           'it_supfun'//trim(adjustl(yaml_toa(0,fmt='(i3.3)'))))
     end if
     call getLocalizedBasis(iproc,nproc,at,orbs,rxyz,denspot,GPU,trace,trace_old,fnrm_tmb,&
         info_basis_functions,nlpsp,input%lin%scf_mode,ldiis,input%SIC,tmb,energs, &
         input%lin%nItPrecond,TARGET_FUNCTION_IS_TRACE,input%lin%correctionOrthoconstraint,&
         50,&
         ratio_deltas,ortho_on,input%lin%extra_states,0,1.d-3,input%experimental_mode,input%lin%early_stop)
     reduce_conf=.true.
     call yaml_close_sequence()
     call yaml_close_map()
     call deallocateDIIS(ldiis)
     !call yaml_open_map()
     ! END NEW ############################################################################
 end if


  call nullify_sparsematrix(ham_small) ! nullify anyway

  !!if (iproc==0) then
  !!    call yaml_close_map()
  !!end if

  if (iproc==0) then
      !call yaml_open_sequence('First kernel')
      !call yaml_open_sequence('kernel optimization',label=&
      !                          'it_kernel'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)'))))
      !call yaml_sequence(advance='no')
!      call yaml_open_map('Input Guess kernel ')
!      call yaml_map('Generation method',input%lin%scf_mode) 
      !call yaml_sequence(advance='no')
      !call yaml_open_map(flow=.false.)
      !call yaml_comment('kernel iter:'//yaml_toa(0,fmt='(i6)'),hfill='-')
  end if

  if (input%lin%scf_mode==LINEAR_FOE) then

      call sparse_copy_pattern(tmb%linmat%ovrlp,ham_small,iproc,subname)
      allocate(ham_small%matrix_compr(ham_small%nvctr), stat=istat)
      call memocc(istat, ham_small%matrix_compr, 'ham_small%matrix_compr', subname)

      call get_coeff(iproc,nproc,LINEAR_FOE,orbs,at,rxyz,denspot,GPU,infoCoeff,energs,nlpsp,&
           input%SIC,tmb,fnrm,.true.,.false.,.true.,ham_small,0,0,0,0,input%lin%order_taylor,input%calculate_KS_residue)

      if (input%lin%scf_mode==LINEAR_FOE) then ! deallocate ham_small
         call deallocate_sparsematrix(ham_small,subname)
      end if

  else
      call get_coeff(iproc,nproc,LINEAR_MIXDENS_SIMPLE,orbs,at,rxyz,denspot,GPU,infoCoeff,energs,nlpsp,&
           input%SIC,tmb,fnrm,.true.,.false.,.true.,ham_small,0,0,0,0,input%lin%order_taylor,input%calculate_KS_residue)

      call dcopy(kswfn%orbs%norb,tmb%orbs%eval(1),1,kswfn%orbs%eval(1),1)
      call evaltoocc(iproc,nproc,.false.,input%tel,kswfn%orbs,input%occopt)
      if (bigdft_mpi%iproc ==0) then
         call write_eigenvalues_data(0.1d0,kswfn%orbs,mom_vec_fake)
      end if

  end if



  call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, max(tmb%npsidim_orbs,tmb%npsidim_comp), &
       tmb%orbs, tmb%psi, tmb%collcom_sr)


  if (iproc==0) then
      call yaml_open_map('Hamiltonian update',flow=.true.)
     ! Use this subroutine to write the energies, with some
     ! fake number
     ! to prevent it from writing too much
    call write_energies(0,0,energs,0.d0,0.d0,'',.true.)
  end if


  call sumrho_for_TMBs(iproc, nproc, tmb%Lzd%hgrids(1), tmb%Lzd%hgrids(2), tmb%Lzd%hgrids(3), &
       tmb%collcom_sr, tmb%linmat%denskern, tmb%Lzd%Glr%d%n1i*tmb%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)

  !!!call plot_density(iproc,nproc,'initial',at,rxyz,denspot%dpbox,input%nspin,denspot%rhov)

  ! Mix the density.
  if (input%lin%mixing_after_inputguess .and. (input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE .or. input%lin%scf_mode==LINEAR_FOE)) then
     if (input%experimental_mode) then
         !if (iproc==0) write(*,*) 'WARNING: TAKE 1.d0 MIXING PARAMETER!'
         if (iproc==0) call yaml_map('INFO mixing parameter for this step',1.d0)
         call mix_main(iproc, nproc, 0, input, tmb%Lzd%Glr, 1.d0, &
              denspot, mixdiis, rhopotold, pnrm)
     else
         call mix_main(iproc, nproc, 0, input, tmb%Lzd%Glr, input%lin%alpha_mix_lowaccuracy, &
              denspot, mixdiis, rhopotold, pnrm)
     end if
  end if

  if(input%lin%scf_mode/=LINEAR_MIXPOT_SIMPLE) then
      call dcopy(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
  end if
  if (iproc==0) call yaml_newline()
  call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
  if(iproc==0) call yaml_close_map()
  ! Mix the potential.
  if (input%lin%mixing_after_inputguess .and. input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
     call mix_main(iproc, nproc, 0, input, tmb%Lzd%Glr, input%lin%alpha_mix_lowaccuracy, &
          denspot, mixdiis, rhopotold, pnrm)
  end if


  if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
      call dcopy(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
  end if


  ! Important: Don't use for the rest of the code
  tmb%ham_descr%can_use_transposed = .false.

  if(associated(tmb%ham_descr%psit_c)) then
      iall=-product(shape(tmb%ham_descr%psit_c))*kind(tmb%ham_descr%psit_c)
      deallocate(tmb%ham_descr%psit_c, stat=istat)
      call memocc(istat, iall, 'tmb%ham_descr%psit_c', subname)
  end if
  if(associated(tmb%ham_descr%psit_f)) then
      iall=-product(shape(tmb%ham_descr%psit_f))*kind(tmb%ham_descr%psit_f)
      deallocate(tmb%ham_descr%psit_f, stat=istat)
      call memocc(istat, iall, 'tmb%ham_descr%psit_f', subname)
  end if
  
  !if (iproc==0) then
  !    call yaml_close_map()
      !call yaml_close_sequence()
      !call yaml_close_sequence()
  !end if
  !!if(iproc==0) write(*,'(1x,a)') '------------------------------------------------------------- Input guess generated.'
  if (iproc==0) call yaml_comment('Input guess generated',hfill='=')
  
  ! Deallocate all local arrays.

  ! Deallocate all types that are not needed any longer.
  call deallocate_orbitals_data(orbs_gauss, subname)

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

END SUBROUTINE inputguessConfinement
