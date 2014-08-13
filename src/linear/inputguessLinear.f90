!> @file 
!!   Input guess for the linear version
!! @author
!!   Copyright (C) 2011-2013 BigDFT group 
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
  use gaussians, only: gaussian_basis, deallocate_gwf, nullify_gaussian_basis
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use yaml_output
  use sparsematrix_base, only: sparse_matrix, sparse_matrix_null, deallocate_sparse_matrix
  use sparsematrix_init, only: matrixindex_in_compressed
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
  real(dp), dimension(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin),intent(inout) :: rhopotold
  type(energy_terms),intent(inout) :: energs
  real(kind=8),dimension(3,at%astruct%nat),intent(in),optional :: locregcenters

  ! Local variables
  type(gaussian_basis) :: G !basis for davidson IG
  character(len=*), parameter :: subname='inputguessConfinement'
  integer :: istat,iall,iat,nspin_ig,iorb,nvirt,norbat,methTransformOverlap
  real(gp) :: hxh,hyh,hzh,eks,fnrm,V3prb,x0,tt
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(gp), dimension(:), allocatable :: locrad
  real(wp), dimension(:,:,:), pointer :: psigau
  integer, dimension(:), allocatable :: norbsPerAt, mapping, inversemapping, minorbs_type, maxorbs_type
  logical, dimension(:), allocatable :: covered, type_covered
  !real(kind=8), dimension(:,:), allocatable :: aocc
  integer, dimension(:,:), allocatable :: nl_copy 
  integer :: ist,jorb,iadd,ii,jj,ityp,itype,iortho
  integer :: jlr,iiorb,ispin
  integer :: infoCoeff, jproc
  type(orbitals_data) :: orbs_gauss
  type(GPU_pointers) :: GPUe
  character(len=2) :: symbol
  real(kind=8) :: rcov,rprb,pnrm
  !real(kind=8) :: ehomo,amu
  integer :: nsccode,mxpl,mxchg,inl
  type(mixrhopotDIISParameters) :: mixdiis
  logical :: finished, can_use_ham
  !type(confpot_data), dimension(:), allocatable :: confdatarrtmp
  integer :: info_basis_functions, order_taylor, i, ilr, iii
  real(kind=8) :: ratio_deltas, trace, trace_old, fnrm_tmb
  logical :: ortho_on, reduce_conf, rho_negative
  type(localizedDIISParameters) :: ldiis
  real(wp), dimension(:,:,:), pointer :: mom_vec_fake

  call f_routine(id=subname)

  call nullify_orbitals_data(orbs_gauss)
  call nullify_gaussian_basis(G)
  nullify(mom_vec_fake)

  ! Allocate some arrays we need for the input guess.
  norbsc_arr = f_malloc((/at%natsc+1,input%nspin/),id='norbsc_arr')
  locrad = f_malloc(at%astruct%nat,id='locrad')
  norbsPerAt = f_malloc(at%astruct%nat,id='norbsPerAt')
  mapping = f_malloc(tmb%orbs%norb,id='mapping')
  covered = f_malloc(tmb%orbs%norb,id='covered')
  inversemapping = f_malloc(tmb%orbs%norb,id='inversemapping')


  GPUe = GPU

  ! Spin for inputguess orbitals
  if (input%nspin == 4) then
     nspin_ig=1
  else
     nspin_ig=input%nspin
  end if

  ! Keep the natural occupations
  
  nl_copy=f_malloc((/0.to.3,1.to.at%astruct%nat/),id='nl_copy')
  do iat=1,at%astruct%nat
     nl_copy(:,iat)=at%aoig(iat)%nl
  end do

!!$  allocate(aocc(32,at%astruct%nat),stat=istat)
!!$  call memocc(istat,aocc,'aocc',subname)
!!$  call vcopy(32*at%astruct%nat, at%aocc(1,1), 1, aocc(1,1), 1)


  ! Determine how many atomic orbitals we have. Maybe we have to increase this number to more than
  ! its 'natural' value.
  norbat=0
  ist=0
  do iat=1,at%astruct%nat
      ii=input%lin%norbsPerType(at%astruct%iatype(iat))
      jj=at%aoig(iat)%nao
      if (jj < ii) then
         call f_err_throw('The number of basis functions asked per type'//&
              ' is exceeding the number of IG atomic orbitals'//&
              ', modify the electronic configuration of input atom '//&
              trim(at%astruct%atomnames(at%astruct%iatype(iat))),&
              err_name='BIGDFT_INPUT_VARIABLES_ERROR')
         call f_release_routine()
         return
      end if

!!$      iadd=0
!!$      do 
!!$          ! Count the number of atomic orbitals and increase the number if necessary until we have more
!!$          ! (or equal) atomic orbitals than basis functions per atom.
!!$         !jj=1*nint(at%aocc(1,iat))+3*nint(at%aocc(3,iat))+&
!!$         !      5*nint(at%aocc(7,iat))+7*nint(at%aocc(13,iat))
!!$         jj=sum(at%aoig(iat)%nl)
!!$
!!$          if(jj>=ii) then
!!$              ! we have enough atomic orbitals
!!$              exit
!!$          else
!!$             ! add additional orbitals
!!$             iadd=iadd+1
!!$             select case(iadd)
!!$             case(1) 
!!$                at%aocc(1,iat)=1.d0
!!$             case(2) 
!!$                at%aocc(3,iat)=1.d0
!!$             case(3) 
!!$                at%aocc(7,iat)=1.d0
!!$             case(4) 
!!$                at%aocc(13,iat)=1.d0
!!$             case default 
!!$                write(*,'(1x,a)') 'ERROR: more than 16 basis functions per atom are not possible!'
!!$                stop
!!$             end select
!!$          end if
!!$       end do

      norbsPerAt(iat)=jj
      norbat=norbat+norbsPerAt(iat)
  end do

  ! This array gives a mapping from the 'natural' orbital distribution (i.e. simply counting up the atoms) to
  ! our optimized orbital distribution (determined by in orbs%inwhichlocreg).
  iiorb=0
  covered=.false.
  if (present(locregcenters)) then
      do ispin=1,input%nspin
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
      end do
  else
      do ispin=1,input%nspin
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
  tmb%can_use_transposed=.false.

  ii=0
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      if (tmb%orbs%spinsgn(iiorb)>0.d0) then
          do i=1,tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
              ii=ii+1
              write(600,'(a,4i9,es16.7)') 'iorb, iiorb, ilr, i, val', iorb, iiorb, ilr, i, tmb%psi(ii)
          end do
      else
          do i=1,tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
              ii=ii+1
              write(610,'(a,4i9,es16.7)') 'iorb, iiorb, ilr, i, val', iorb, iiorb, ilr, i, tmb%psi(ii)
          end do
      end if
  end do

  call f_free_ptr(psigau)

  call deallocate_gwf(G)
  ! Deallocate locrad, which is not used any longer.
  call f_free(locrad)

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
  !call to_zero(tmb%linmat%denskern%nvctr, tmb%linmat%denskern%matrix_compr(1))
  call to_zero(tmb%linmat%l%nvctr*input%nspin, tmb%linmat%kernel_%matrix_compr(1))
  do iorb=1,tmb%orbs%norb
     !ii=matrixindex_in_compressed(tmb%linmat%denskern,iorb,iorb)
     ii=matrixindex_in_compressed(tmb%linmat%l,iorb,iorb)
     !tmb%linmat%denskern%matrix_compr(ii)=1.d0*tmb%orbs%occup(inversemapping(iorb))
     !tmb%linmat%denskern%matrix_compr(ii)=1.d0*tmb%orbs%occup(iorb)
     tmb%linmat%kernel_%matrix_compr(ii)=1.d0*tmb%orbs%occup(iorb)
  end do

  !Calculate the density in the new scheme
  call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, max(tmb%npsidim_orbs,tmb%npsidim_comp), &
       tmb%orbs, tmb%psi, tmb%collcom_sr)
  call sumrho_for_TMBs(iproc, nproc, tmb%Lzd%hgrids(1), tmb%Lzd%hgrids(2), tmb%Lzd%hgrids(3), &
       tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, denspot%dpbox%ndimrhopot, &
       denspot%rhov, rho_negative)
  if (rho_negative) then
      call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
      !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
      !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
      !!call clean_rho(iproc, nproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
  end if


  !!do istat=1,size(denspot%rhov)
  !!    write(300+iproc,*) istat, denspot%rhov(istat)
  !!end do 
  !!call mpi_finalize(istat)
  !!stop


  jj=0
  do ispin=1,input%nspin
      do ii=1,size(denspot%rhov)/input%nspin
          jj=jj+1
          write(9600+10*iproc+ispin,'(a,i9,es16.5)') 'ii, val', ii, denspot%rhov(jj)
      end do
  end do

  if (input%lin%mixing_after_inputguess) then
      if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE .or. input%lin%scf_mode==LINEAR_FOE &
           .or. input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
          call vcopy(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
          ! initial setting of the old charge density
          call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,0.0d0,denspot%mix,&
               denspot%rhov,1,denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
               at%astruct%cell_dim(1)*at%astruct%cell_dim(2)*at%astruct%cell_dim(3),&
               pnrm,denspot%dpbox%nscatterarr)
      end if
  end if
  !do ii=1,size(denspot%rhov)
  !    write(9600+iproc,'(a,2i9,es16.5)') 'ii, mod(ii-1,size(denspot%rhov)/2)+1, val', &
  !        ii, mod(ii-1,size(denspot%rhov)/2)+1, denspot%rhov(ii)
  !end do
  jj=0
  do ispin=1,input%nspin
      do ii=1,size(denspot%rhov)/input%nspin
          jj=jj+1
          write(9700+10*iproc+ispin,'(a,i9,es16.5)') 'ii, val', ii, denspot%rhov(jj)
          !if (ispin==2) then
          !    denspot%rhov(jj)=denspot%rhov(jj-size(denspot%rhov)/input%nspin)
          !    write(9800+10*iproc+ispin,'(a,i9,es16.5)') 'ii, val', ii, denspot%rhov(jj)
          !end if
      end do
  end do
  call updatePotential(input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
  write(*,*) 'after first updatePotential'

  !!write(*,'(a,4i8)') 'iproc, denspot%dpbox%n3d, denspot%dpbox%n3p, denspot%dpbox%nscatterarr(iproc,2)', &
  !!                    iproc, denspot%dpbox%n3d, denspot%dpbox%n3p, denspot%dpbox%nscatterarr(iproc,2)
  !!iall=0
  !!do jproc=0,nproc-1
  !!    do istat=1,tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%nscatterarr(jproc,2)
  !!        iall=iall+1
  !!        if (iproc==jproc) write(500+iproc,*) iall, denspot%rhov(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%i3xcsh+istat)
  !!    end do
  !!end do
  !!call mpi_finalize(istat)
  !!stop

  if (input%lin%mixing_after_inputguess) then
      if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
          call vcopy(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
          ! initial setting of the old charge density
          call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,0.0d0,denspot%mix,&
               denspot%rhov,1,denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
               at%astruct%cell_dim(1)*at%astruct%cell_dim(2)*at%astruct%cell_dim(3),&
               pnrm,denspot%dpbox%nscatterarr)
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
      if (iproc==0) call yaml_map('orthonormalization of input guess','standard')
      ! CHEATING here and passing tmb%linmat%denskern instead of tmb%linmat%inv_ovrlp
      !write(*,'(a,i4,4i8)') 'IG: iproc, lbound, ubound, minval, maxval',&
      !iproc, lbound(tmb%linmat%inv_ovrlp%matrixindex_in_compressed_fortransposed,2),&
      !ubound(tmb%linmat%inv_ovrlp%matrixindex_in_compressed_fortransposed,2),&
      !minval(tmb%collcom%indexrecvorbital_c),maxval(tmb%collcom%indexrecvorbital_c)
      !!if (iproc==0) write(*,*) 'WARNING: no ortho in inguess'
      methTransformOverlap=-1

      !iii=0
      !do iorb=1,tmb%orbs%norb
      !    ilr=tmb%orbs%inwhichlocreg(iorb)
      !    ii=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      !    if (tmb%orbs%spinsgn(iorb)>0.d0) then
      !        do i=1,ii
      !            iii=iii+1
      !            write(550,*) 'i, iorb, val', i, iorb, tmb%psi(iii)
      !        end do
      !    else
      !        do i=1,ii
      !            iii=iii+1
      !            write(551,*) 'i, iorb, val', i, iorb, tmb%psi(iii)
      !        end do
      !    end if
      !end do
      call orthonormalizeLocalized(iproc, nproc, methTransformOverlap, 1.d0, tmb%npsidim_orbs, tmb%orbs, tmb%lzd, &
           tmb%linmat%s, tmb%linmat%l, &
           tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed, &
           tmb%foe_obj)
            
 else
     ! Iterative orthonomalization
     !!if(iproc==0) write(*,*) 'calling generalized orthonormalization'
     if (iproc==0) call yaml_map('orthonormalization of input guess','generalized')
     maxorbs_type = f_malloc(at%astruct%ntypes,id='maxorbs_type')
     minorbs_type = f_malloc(at%astruct%ntypes,id='minorbs_type')
     type_covered = f_malloc(at%astruct%ntypes,id='type_covered')
     minorbs_type(1:at%astruct%ntypes)=0
     iortho=0
     ortho_loop: do
         finished=.true.
         type_covered=.false.
         do iat=1,at%astruct%nat
             itype=at%astruct%iatype(iat)
             if (type_covered(itype)) cycle
             type_covered(itype)=.true.
             !jj=1*ceiling(aocc(1,iat))+3*ceiling(aocc(3,iat))+&
             !     5*ceiling(aocc(7,iat))+7*ceiling(aocc(13,iat))
             jj=nl_copy(0,iat)+3*nl_copy(1,iat)+5*nl_copy(2,iat)+7*nl_copy(3,iat)
             maxorbs_type(itype)=jj
             !should not enter in the conditional below due to the raise of the exception above
             if (jj<input%lin%norbsPerType(at%astruct%iatype(iat))) then
                 finished=.false.
                 increase_count: do inl=1,4
                    if (nl_copy(inl,iat)==0) then
                       nl_copy(inl,iat)=1
                       call f_err_throw('InputguessLinear: Should not be here',&
                            err_name='BIGDFT_RUNTIME_ERROR')
                       exit increase_count
                    end if
                 end do increase_count
!!$                 if (ceiling(aocc(1,iat))==0) then
!!$                     aocc(1,iat)=1.d0
!!$                 else if (ceiling(aocc(3,iat))==0) then
!!$                     aocc(3,iat)=1.d0
!!$                 else if (ceiling(aocc(7,iat))==0) then
!!$                     aocc(7,iat)=1.d0
!!$                 else if (ceiling(aocc(13,iat))==0) then
!!$                     aocc(13,iat)=1.d0
!!$                 end if
             end if
         end do
         if (iortho>0) then
             call gramschmidt_subset(iproc, nproc, -1, tmb%npsidim_orbs, &                                  
                  tmb%orbs, at, minorbs_type, maxorbs_type, tmb%lzd, tmb%linmat%s, &
                  tmb%linmat%l, tmb%collcom, tmb%orthpar, &
                  tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)
         end if
         call orthonormalize_subset(iproc, nproc, -1, tmb%npsidim_orbs, &                                  
              tmb%orbs, at, minorbs_type, maxorbs_type, tmb%lzd, tmb%linmat%s, &
              tmb%linmat%l, tmb%collcom, tmb%orthpar, &
              tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)
         if (finished) exit ortho_loop
         iortho=iortho+1
         minorbs_type(1:at%astruct%ntypes)=maxorbs_type(1:at%astruct%ntypes)+1
     end do ortho_loop
     call f_free(maxorbs_type)
     call f_free(minorbs_type)
     call f_free(type_covered)

 end if

 call f_free(nl_copy)
 !!!!! adding some noise
 !!Write(*,*) 'warning: add some noise!'
 !!do istat=1,size(tmb%psi)
 !!    call random_number(tt)
 !!    tt=tt-0.5d0
 !!    tt=tt*0.6d0
 !!    tmb%psi(istat)=tmb%psi(istat)*(1.d0+tt)
 !!end do
 !!tmb%can_use_transposed=.false.


!!$ iall=-product(shape(aocc))*kind(aocc)
!!$ deallocate(aocc,stat=istat)
!!$ call memocc(istat, iall,'aocc',subname)

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
         !call yaml_mapping_close()
         call yaml_comment('Extended input guess for experimental mode',hfill='-')
         call yaml_mapping_open('Extended input guess')
         call yaml_sequence_open('support function optimization',label=&
                                           'it_supfun'//trim(adjustl(yaml_toa(0,fmt='(i3.3)'))))
     end if
     order_taylor=input%lin%order_taylor ! since this is intent(inout)
     call getLocalizedBasis(iproc,nproc,at,orbs,rxyz,denspot,GPU,trace,trace_old,fnrm_tmb,&
         info_basis_functions,nlpsp,input%lin%scf_mode,ldiis,input%SIC,tmb,energs, &
         input%lin%nItPrecond,TARGET_FUNCTION_IS_TRACE,input%lin%correctionOrthoconstraint,&
         50,&
         ratio_deltas,ortho_on,input%lin%extra_states,0,1.d-3,input%experimental_mode,input%lin%early_stop,&
         input%lin%gnrm_dynamic, input%lin%min_gnrm_for_dynamic, &
         can_use_ham, order_taylor, input%lin%max_inversion_error, input%kappa_conv, input%method_updatekernel,&
         input%purification_quickreturn, input%correction_co_contra)
     reduce_conf=.true.
     call yaml_sequence_close()
     call yaml_mapping_close()
     call deallocateDIIS(ldiis)
     !call yaml_mapping_open()
     ! END NEW ############################################################################
 end if



  !!if (iproc==0) then
  !!    call yaml_mapping_close()
  !!end if

  if (iproc==0) then
      !call yaml_sequence_open('First kernel')
      !call yaml_sequence_open('kernel optimization',label=&
      !                          'it_kernel'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)'))))
      !call yaml_sequence(advance='no')
!      call yaml_mapping_open('Input Guess kernel ')
!      call yaml_map('Generation method',input%lin%scf_mode) 
      !call yaml_sequence(advance='no')
      !call yaml_mapping_open(flow=.false.)
      !call yaml_comment('kernel iter:'//yaml_toa(0,fmt='(i6)'),hfill='-')
  end if

  ii=0
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      if (tmb%orbs%spinsgn(iiorb)>0.d0) then
          do i=1,tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
              ii=ii+1
              write(1600,'(a,4i9,es16.7)') 'iorb, iiorb, ilr, i, val', iorb, iiorb, ilr, i, tmb%psi(ii)
          end do
      else
          do i=1,tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
              ii=ii+1
              write(1610,'(a,4i9,es16.7)') 'iorb, iiorb, ilr, i, val', iorb, iiorb, ilr, i, tmb%psi(ii)
          end do
      end if
  end do

  order_taylor=input%lin%order_taylor ! since this is intent(inout)
  if (input%lin%scf_mode==LINEAR_FOE) then
      call get_coeff(iproc,nproc,LINEAR_FOE,orbs,at,rxyz,denspot,GPU,infoCoeff,energs,nlpsp,&
           input%SIC,tmb,fnrm,.true.,.false.,.true.,0,0,0,0,order_taylor,input%lin%max_inversion_error,&
           input%purification_quickreturn,&
           input%calculate_KS_residue,input%calculate_gap)
  else

      call get_coeff(iproc,nproc,LINEAR_MIXDENS_SIMPLE,orbs,at,rxyz,denspot,GPU,infoCoeff,energs,nlpsp,&
           input%SIC,tmb,fnrm,.true.,.false.,.true.,0,0,0,0,order_taylor,input%lin%max_inversion_error,&
           input%purification_quickreturn,&
           input%calculate_KS_residue,input%calculate_gap)

      !call vcopy(kswfn%orbs%norb,tmb%orbs%eval(1),1,kswfn%orbs%eval(1),1)
      ! Keep the ocupations for the moment.. maybe to be activated later (with a better if statement)
      if (input%Tel > 0.0_gp) then
          call evaltoocc(iproc,nproc,.false.,input%tel,kswfn%orbs,input%occopt)
      end if
      if (bigdft_mpi%iproc ==0) then
         call write_eigenvalues_data(0.1d0,kswfn%orbs,mom_vec_fake)
      end if

  end if


  call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, max(tmb%npsidim_orbs,tmb%npsidim_comp), &
       tmb%orbs, tmb%psi, tmb%collcom_sr)

  if (iproc==0) then
      call yaml_mapping_open('Hamiltonian update',flow=.true.)
     ! Use this subroutine to write the energies, with some
     ! fake number
     ! to prevent it from writing too much
    call write_energies(0,0,energs,0.d0,0.d0,'',.true.)
  end if

  call sumrho_for_TMBs(iproc, nproc, tmb%Lzd%hgrids(1), tmb%Lzd%hgrids(2), tmb%Lzd%hgrids(3), &
       tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, denspot%dpbox%ndimrhopot, &
       denspot%rhov, rho_negative)
  if (rho_negative) then
      call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
      !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
      !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
      !!call clean_rho(iproc, nproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
  end if

  !!!call plot_density(iproc,nproc,'initial',at,rxyz,denspot%dpbox,input%nspin,denspot%rhov)

  ! Mix the density.
  if (input%lin%mixing_after_inputguess .and. &
          (input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE .or. input%lin%scf_mode==LINEAR_FOE)) then
     if (input%experimental_mode) then
         !if (iproc==0) write(*,*) 'WARNING: TAKE 1.d0 MIXING PARAMETER!'
         if (iproc==0) call yaml_map('INFO mixing parameter for this step',1.d0)
         !!call mix_main(iproc, nproc, input%lin%scf_mode, 0, input, tmb%Lzd%Glr, 1.d0, &
         !!     denspot, mixdiis, rhopotold, pnrm)
         call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,0.d0,denspot%mix,&
              denspot%rhov,1,denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
              at%astruct%cell_dim(1)*at%astruct%cell_dim(2)*at%astruct%cell_dim(3),&
              pnrm,denspot%dpbox%nscatterarr)
     else
         !!call mix_main(iproc, nproc, input%lin%scf_mode, 0, input, tmb%Lzd%Glr, input%lin%alpha_mix_lowaccuracy, &
         !!     denspot, mixdiis, rhopotold, pnrm)
         call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,1.d0-input%lin%alpha_mix_lowaccuracy,denspot%mix,&
              denspot%rhov,2,denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
              at%astruct%cell_dim(1)*at%astruct%cell_dim(2)*at%astruct%cell_dim(3),&
              pnrm,denspot%dpbox%nscatterarr)
     end if
  end if

  if(input%lin%scf_mode/=LINEAR_MIXPOT_SIMPLE) then
      call vcopy(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
      ! initial setting of the old charge density
      call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,0.d0,denspot%mix,&
           denspot%rhov,1,denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
           at%astruct%cell_dim(1)*at%astruct%cell_dim(2)*at%astruct%cell_dim(3),&
           pnrm,denspot%dpbox%nscatterarr)
  end if
  if (iproc==0) call yaml_newline()
  call updatePotential(input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
  if(iproc==0) call yaml_mapping_close()
  ! Mix the potential.
  if (input%lin%mixing_after_inputguess .and. input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
     !!call mix_main(iproc, nproc, input%lin%scf_mode, 0, input, tmb%Lzd%Glr, input%lin%alpha_mix_lowaccuracy, &
     !!     denspot, mixdiis, rhopotold, pnrm)
     call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,1.d0-input%lin%alpha_mix_lowaccuracy,denspot%mix,&
          denspot%rhov,2,denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
          at%astruct%cell_dim(1)*at%astruct%cell_dim(2)*at%astruct%cell_dim(3),&
          pnrm,denspot%dpbox%nscatterarr)
  end if


  if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
      call vcopy(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, denspot%rhov(1), 1, rhopotold(1), 1)
      ! initial setting of the old potential
      call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,0.d0,denspot%mix,&
           denspot%rhov,1,denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
           at%astruct%cell_dim(1)*at%astruct%cell_dim(2)*at%astruct%cell_dim(3),&
           pnrm,denspot%dpbox%nscatterarr)
  end if


  ! Important: Don't use for the rest of the code
  tmb%ham_descr%can_use_transposed = .false.

  !if(associated(tmb%ham_descr%psit_c)) then
  !    call f_free_ptr(tmb%ham_descr%psit_c)
  !end if
  !if(associated(tmb%ham_descr%psit_f)) then
  !    call f_free_ptr(tmb%ham_descr%psit_f)
  !end if
  
  !if (iproc==0) then
  !    call yaml_mapping_close()
      !call yaml_sequence_close()
      !call yaml_sequence_close()
  !end if
  !!if(iproc==0) write(*,'(1x,a)') '------------------------------------------------------------- Input guess generated.'
  if (iproc==0) call yaml_comment('Input guess generated',hfill='=')
  
  ! Deallocate all local arrays.

  ! Deallocate all types that are not needed any longer.
  call deallocate_orbitals_data(orbs_gauss)

  ! Deallocate all remaining local arrays.
  call f_free(norbsc_arr)
  call f_free(norbsPerAt)
  call f_free(mapping)
  call f_free(covered)
  call f_free(inversemapping)

  call f_release_routine()

END SUBROUTINE inputguessConfinement
