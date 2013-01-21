!> @file
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

subroutine get_coeff(iproc,nproc,scf_mode,orbs,at,rxyz,denspot,GPU,infoCoeff,&
    ebs,nlpspd,proj,SIC,tmb,fnrm,calculate_overlap_matrix,communicate_phi_for_lsumrho,&
    tmblarge,calculate_ham,ldiis_coeff)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => get_coeff, exceptThisOneA => writeonewave
  use Poisson_Solver
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, scf_mode
   type(orbitals_data),intent(inout) :: orbs
  type(atoms_data),intent(in) :: at
  real(kind=8),dimension(3,at%nat),intent(in) :: rxyz
  type(DFT_local_fields), intent(inout) :: denspot
  type(GPU_pointers),intent(inout) :: GPU
  integer,intent(out) :: infoCoeff
  real(kind=8),intent(out) :: ebs
  real(kind=8),intent(inout) :: fnrm
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd
  real(wp),dimension(nlpspd%nprojel),intent(inout) :: proj
  type(SIC_data),intent(in) :: SIC
  type(DFT_wavefunction),intent(inout) :: tmb
  logical,intent(in):: calculate_overlap_matrix, communicate_phi_for_lsumrho
  type(DFT_wavefunction),intent(inout) :: tmblarge
  logical,intent(in) :: calculate_ham
  type(localizedDIISParameters),intent(inout),optional :: ldiis_coeff

  ! Local variables 
  integer :: istat, iall, iorb, jorb, korb, info, iiorb, ierr, ii, iseg
  integer :: isegsmall, iseglarge, iismall, iilarge, i, is, ie
  real(kind=8),dimension(:),allocatable :: hpsit_c, hpsit_f, ovrlp_compr_small, ham_compr_small
  real(kind=8),dimension(:,:),allocatable :: ham, overlapmatrix, density_kernel
  real(kind=8),dimension(:,:,:),allocatable :: matrixElements
  type(confpot_data),dimension(:),allocatable :: confdatarrtmp
  type(energy_terms) :: energs
  character(len=*),parameter :: subname='get_coeff'
  real(kind=8) :: tmprtr


  if(calculate_ham) then
      call local_potential_dimensions(tmblarge%lzd,tmblarge%orbs,denspot%dpbox%ngatherarr(0,1))
      call start_onesided_communication(iproc, nproc, max(denspot%dpbox%ndimpot,1), denspot%rhov, &
           tmblarge%comgp%nrecvbuf, tmblarge%comgp%recvbuf, tmblarge%comgp, tmblarge%lzd)
  end if

  ! Calculate the overlap matrix if required.
  if(calculate_overlap_matrix) then
      if(.not.tmb%can_use_transposed) then
          if(.not.associated(tmb%psit_c)) then
              allocate(tmb%psit_c(sum(tmb%collcom%nrecvcounts_c)), stat=istat)
              call memocc(istat, tmb%psit_c, 'tmb%psit_c', subname)
          end if
          if(.not.associated(tmb%psit_f)) then
              allocate(tmb%psit_f(7*sum(tmb%collcom%nrecvcounts_f)), stat=istat)
              call memocc(istat, tmb%psit_f, 'tmb%psit_f', subname)
          end if
          call transpose_localized(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
          tmb%can_use_transposed=.true.
      end if
      ! use tmblarge%sparsemat since the matrix compression must be done the large version
      call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmblarge%sparsemat, tmb%collcom, tmb%psit_c, &
           tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%ovrlp%matrix_compr)
  end if

  ! Post the p2p communications for the density. (must not be done in inputguess)
  if(communicate_phi_for_lsumrho) then
      call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, tmb%orbs, tmb%psi, tmb%collcom_sr)
  end if

  if(iproc==0) write(*,'(1x,a)') '----------------------------------- Determination of the orbitals in this new basis.'

  ! Calculate the Hamiltonian matrix if it is not already present.
  if(calculate_ham) then

      allocate(confdatarrtmp(tmb%orbs%norbp))
      call default_confinement_data(confdatarrtmp,tmb%orbs%norbp)

      call small_to_large_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmb%orbs, tmblarge%orbs, tmb%psi, tmblarge%psi)

      if (tmblarge%orbs%npsidim_orbs > 0) call to_zero(tmblarge%orbs%npsidim_orbs,tmblarge%hpsi(1))

      call NonLocalHamiltonianApplication(iproc,at,tmblarge%orbs,rxyz,&
           proj,tmblarge%lzd,nlpspd,tmblarge%psi,tmblarge%hpsi,energs%eproj)
      ! only kinetic as waiting for communications
      call LocalHamiltonianApplication(iproc,nproc,at,tmblarge%orbs,&
           tmblarge%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmblarge%psi,tmblarge%hpsi,&
           energs,SIC,GPU,3,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmblarge%comgp)
      call full_local_potential(iproc,nproc,tmblarge%orbs,tmblarge%lzd,2,denspot%dpbox,denspot%rhov,denspot%pot_work, &
           tmblarge%comgp)
      !call wait_p2p_communication(iproc, nproc, tmblarge%comgp)
      ! only potential
      call LocalHamiltonianApplication(iproc,nproc,at,tmblarge%orbs,&
           tmblarge%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmblarge%psi,tmblarge%hpsi,&
           energs,SIC,GPU,2,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmblarge%comgp)
      call timing(iproc,'glsynchham1','ON') !lr408t
      call SynchronizeHamiltonianApplication(nproc,tmblarge%orbs,tmblarge%lzd,GPU,tmblarge%hpsi,&
           energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
      call timing(iproc,'glsynchham1','OF') !lr408t
      deallocate(confdatarrtmp)

      !DEBUG
      !!if(iproc==0) then
      !! print *,'Ekin,Epot,Eproj,Eh,Exc,Evxc',energs%ekin,energs%epot,energs%eproj,energs%eh,energs%exc,energs%evxc
      !!end if
      !END DEBUG

      iall=-product(shape(denspot%pot_work))*kind(denspot%pot_work)
      deallocate(denspot%pot_work, stat=istat)
      call memocc(istat, iall, 'denspot%pot_work', subname)

      if(iproc==0) write(*,'(1x,a)') 'Hamiltonian application done.'


      ! Calculate the matrix elements <phi|H|phi>.
      if(.not.tmblarge%can_use_transposed) then
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

          allocate(tmblarge%psit_c(tmblarge%collcom%ndimind_c), stat=istat)
          call memocc(istat, tmblarge%psit_c, 'tmblarge%psit_c', subname)
          allocate(tmblarge%psit_f(7*tmblarge%collcom%ndimind_f), stat=istat)
          call memocc(istat, tmblarge%psit_f, 'tmblarge%psit_f', subname)
          call transpose_localized(iproc, nproc, tmblarge%orbs,  tmblarge%collcom, &
               tmblarge%psi, tmblarge%psit_c, tmblarge%psit_f, tmblarge%lzd)
          tmblarge%can_use_transposed=.true.
      end if

      allocate(hpsit_c(tmblarge%collcom%ndimind_c))
      call memocc(istat, hpsit_c, 'hpsit_c', subname)
      allocate(hpsit_f(7*tmblarge%collcom%ndimind_f))
      call memocc(istat, hpsit_f, 'hpsit_f', subname)
      call transpose_localized(iproc, nproc, tmblarge%orbs,  tmblarge%collcom, &
           tmblarge%hpsi, hpsit_c, hpsit_f, tmblarge%lzd)
      call calculate_overlap_transposed(iproc, nproc, tmblarge%orbs, tmblarge%sparsemat, tmblarge%collcom, &
           tmblarge%psit_c, hpsit_c, tmblarge%psit_f, hpsit_f, tmb%linmat%ham%matrix_compr)
      iall=-product(shape(hpsit_c))*kind(hpsit_c)
      deallocate(hpsit_c, stat=istat)
      call memocc(istat, iall, 'hpsit_c', subname)
      iall=-product(shape(hpsit_f))*kind(hpsit_f)
      deallocate(hpsit_f, stat=istat)
      call memocc(istat, iall, 'hpsit_f', subname)

  else
      if(iproc==0) write(*,*) 'No Hamiltonian application required.'
  end if


  if (scf_mode/=LINEAR_FOE) then
      allocate(ham(tmblarge%orbs%norb,tmblarge%orbs%norb), stat=istat)
      call memocc(istat, ham, 'ham', subname)
      call uncompressMatrix(tmblarge%orbs%norb, tmblarge%sparsemat, tmb%linmat%ham%matrix_compr, ham)
      allocate(overlapmatrix(tmblarge%orbs%norb,tmblarge%orbs%norb), stat=istat)
      call memocc(istat, overlapmatrix, 'overlapmatrix', subname)
      call uncompressMatrix(tmblarge%orbs%norb, tmblarge%sparsemat, tmb%linmat%ovrlp%matrix_compr, overlapmatrix)
  end if

  ! DEBUG LR
  !!if (iproc==0) then
  !!   open(11)
  !!   open(12)
  !!   do iorb=1,tmb%orbs%norb
  !!      do jorb=1,tmb%orbs%norb
  !!          write(11+iproc,*) iorb,jorb,ham(jorb,iorb)
  !!          write(12+iproc,*) iorb,jorb,overlapmatrix(jorb,iorb)
  !!      end do
  !!   end do
  !!   close(11)
  !!   close(12)
  !!end if
  ! END DEBUG LR

  ! Diagonalize the Hamiltonian.
  if(scf_mode==LINEAR_MIXPOT_SIMPLE .or. scf_mode==LINEAR_MIXDENS_SIMPLE) then
      ! Keep the Hamiltonian and the overlap since they will be overwritten by the diagonalization.
      allocate(matrixElements(tmb%orbs%norb,tmb%orbs%norb,2), stat=istat)
      call memocc(istat, matrixElements, 'matrixElements', subname)
      call dcopy(tmb%orbs%norb**2, ham(1,1), 1, matrixElements(1,1,1), 1)
      call dcopy(tmb%orbs%norb**2, overlapmatrix(1,1), 1, matrixElements(1,1,2), 1)
      if(tmb%orthpar%blocksize_pdsyev<0) then
          if(iproc==0) write(*,'(1x,a)',advance='no') 'Diagonalizing the Hamiltonian, sequential version... '
          call diagonalizeHamiltonian2(iproc, tmb%orbs, matrixElements(1,1,1), matrixElements(1,1,2), tmb%orbs%eval)
      else
          if(iproc==0) write(*,'(1x,a)',advance='no') 'Diagonalizing the Hamiltonian, parallel version... '
          call dsygv_parallel(iproc, nproc, tmb%orthpar%blocksize_pdsyev, tmb%orthpar%nproc_pdsyev, &
               bigdft_mpi%mpi_comm, 1, 'v', 'l',tmb%orbs%norb, &
               matrixElements(1,1,1), tmb%orbs%norb, matrixElements(1,1,2), tmb%orbs%norb, tmb%orbs%eval, info)
      end if
      if(iproc==0) write(*,'(a)') 'done.'

      do iorb=1,orbs%norb
          call dcopy(tmb%orbs%norb, matrixElements(1,iorb,1), 1, tmb%wfnmd%coeff(1,iorb), 1)
      end do
      infoCoeff=0

      ! Write some eigenvalues. Don't write all, but only a few around the last occupied orbital.
      if(iproc==0) then
          write(*,'(1x,a)') '-------------------------------------------------'
          write(*,'(1x,a)') 'some selected eigenvalues:'
          do iorb=max(orbs%norb-8,1),min(orbs%norb+8,tmb%orbs%norb)
              if(iorb==orbs%norb) then
                  write(*,'(3x,a,i0,a,es20.12,a)') 'eval(',iorb,')= ',tmb%orbs%eval(iorb),'  <-- last occupied orbital'
              else if(iorb==orbs%norb+1) then
                  write(*,'(3x,a,i0,a,es20.12,a)') 'eval(',iorb,')= ',tmb%orbs%eval(iorb),'  <-- first virtual orbital'
              else
                  write(*,'(3x,a,i0,a,es20.12)') 'eval(',iorb,')= ',tmb%orbs%eval(iorb)
              end if
          end do
          write(*,'(1x,a)') '-------------------------------------------------'
          write(*,'(1x,a,2es24.16)') 'lowest, highest ev:',tmb%orbs%eval(1),tmb%orbs%eval(tmb%orbs%norb)
      end if

      ! keep the eigenvalues for the preconditioning - instead should take h_alpha,alpha for both cases
      call vcopy(tmb%orbs%norb, tmb%orbs%eval(1), 1, tmblarge%orbs%eval(1), 1)
      ! instead just use -0.5 everywhere
      !tmb%orbs%eval(:) = -0.5_dp
      !tmblarge%orbs%eval(:) = -0.5_dp

      iall=-product(shape(matrixElements))*kind(matrixElements)
      deallocate(matrixElements, stat=istat)
      call memocc(istat, iall, 'matrixElements', subname)
  else if (scf_mode==LINEAR_DIRECT_MINIMIZATION) then
      if(.not.present(ldiis_coeff)) stop 'ldiis_coeff must be present for scf_mode==LINEAR_DIRECT_MINIMIZATION'
      call optimize_coeffs(iproc, nproc, orbs, ham, overlapmatrix, tmb, ldiis_coeff, fnrm)
  end if


  if (scf_mode/=LINEAR_FOE) then

      allocate(density_kernel(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
      call memocc(istat, density_kernel, 'density_kernel', subname)
      call calculate_density_kernel(iproc, nproc, .true., orbs, tmb%orbs, &
           tmb%wfnmd%coeff, density_kernel)
      call compress_matrix_for_allreduce(tmblarge%orbs%norb, tmblarge%sparsemat, density_kernel, tmb%linmat%denskern%matrix_compr)
      iall=-product(shape(density_kernel))*kind(density_kernel)
      deallocate(density_kernel, stat=istat)
      call memocc(istat, iall, 'density_kernel', subname)

      ebs=0.d0
      ii=0
      do iseg=1,tmblarge%sparsemat%nseg
          do jorb=tmblarge%sparsemat%keyg(1,iseg),tmblarge%sparsemat%keyg(2,iseg)
              ii=ii+1
              ebs = ebs + tmb%linmat%denskern%matrix_compr(ii)*tmb%linmat%ham%matrix_compr(ii)
          end do  
      end do

      ! Calculate the KS eigenvalues - needed for Pulay
      call to_zero(orbs%norb, orbs%eval(1))
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          do jorb=1,tmb%orbs%norb
              do korb=1,tmb%orbs%norb
                  orbs%eval(iiorb) = orbs%eval(iiorb) + &
                                     tmb%wfnmd%coeff(jorb,iiorb)*tmb%wfnmd%coeff(korb,iiorb)*ham(jorb,korb)
              end do
          end do
      end do
      call mpiallred(orbs%eval(1), orbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)

      iall=-product(shape(ham))*kind(ham)
      deallocate(ham, stat=istat)
      call memocc(istat, iall, 'ham', subname)
      iall=-product(shape(overlapmatrix))*kind(overlapmatrix)
      deallocate(overlapmatrix, stat=istat)
      call memocc(istat, iall, 'overlapmatrix', subname)

  else ! foe

      allocate(ovrlp_compr_small(tmb%sparsemat%nvctr), stat=istat)
      call memocc(istat, ovrlp_compr_small, 'ovrlp_compr_small', subname)
      allocate(ham_compr_small(tmb%sparsemat%nvctr), stat=istat)
      call memocc(istat, ham_compr_small, 'ham_compr_small', subname)

      iismall=0
      iseglarge=1
      do isegsmall=1,tmb%sparsemat%nseg
          do
              is=max(tmb%sparsemat%keyg(1,isegsmall),tmblarge%sparsemat%keyg(1,iseglarge))
              ie=min(tmb%sparsemat%keyg(2,isegsmall),tmblarge%sparsemat%keyg(2,iseglarge))
              iilarge=tmblarge%sparsemat%keyv(iseglarge)-tmblarge%sparsemat%keyg(1,iseglarge)
              do i=is,ie
                  iismall=iismall+1
                  ovrlp_compr_small(iismall)=tmb%linmat%ovrlp%matrix_compr(iilarge+i)
                  ham_compr_small(iismall)=tmb%linmat%ham%matrix_compr(iilarge+i)
              end do
              if (ie>=is) exit
              iseglarge=iseglarge+1
          end do
      end do

      tmprtr=0.d0
      call foe(iproc, nproc, tmb, tmblarge, orbs, tmb%wfnmd%evlow, tmb%wfnmd%evhigh, &
           tmb%wfnmd%fscale, tmb%wfnmd%ef, tmprtr, 2, &
           ham_compr_small, ovrlp_compr_small, tmb%wfnmd%bisection_shift, tmb%linmat%denskern%matrix_compr, ebs)
      ! Eigenvalues not available, therefore take -.5d0
      tmb%orbs%eval=-.5d0
      tmblarge%orbs%eval=-.5d0

      iall=-product(shape(ovrlp_compr_small))*kind(ovrlp_compr_small)
      deallocate(ovrlp_compr_small, stat=istat)
      call memocc(istat, iall, 'ovrlp_compr_small', subname)
      iall=-product(shape(ham_compr_small))*kind(ham_compr_small)
      deallocate(ham_compr_small, stat=istat)
      call memocc(istat, iall, 'ham_compr_small', subname)

  end if


end subroutine get_coeff



subroutine getLocalizedBasis(iproc,nproc,at,orbs,rxyz,denspot,GPU,trH,trH_old,&
    fnrm,infoBasisFunctions,nlpspd,scf_mode, proj,ldiis,SIC,tmb,tmblarge,energs_base,&
    reduce_conf,fix_supportfunctions,nit_precond,target_function,&
    correction_orthoconstraint,nit_basis,deltaenergy_multiplier_TMBexit,deltaenergy_multiplier_TMBfix)
  !
  ! Purpose:
  ! ========
  !   Calculates the localized basis functions phi. These basis functions are obtained by adding a
  !   quartic potential centered on the atoms to the ordinary Hamiltonian. The eigenfunctions are then
  !   determined by minimizing the trace until the gradient norm is below the convergence criterion.
  use module_base
  use module_types
  use module_interfaces, except_this_one => getLocalizedBasis, except_this_one_A => writeonewave
  !  use Poisson_Solver
  !use allocModule
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  integer,intent(out) :: infoBasisFunctions
  type(atoms_data), intent(in) :: at
  type(orbitals_data) :: orbs
  real(kind=8),dimension(3,at%nat) :: rxyz
  type(DFT_local_fields), intent(inout) :: denspot
  type(GPU_pointers), intent(inout) :: GPU
  real(kind=8),intent(out) :: trH, fnrm
  real(kind=8),intent(inout) :: trH_old
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd
  integer,intent(in) :: scf_mode
  real(wp),dimension(nlpspd%nprojel),intent(inout) :: proj
  type(localizedDIISParameters),intent(inout) :: ldiis
  type(DFT_wavefunction),target,intent(inout) :: tmb
  type(SIC_data) :: SIC !<parameters for the SIC methods
  type(DFT_wavefunction),target,intent(inout) :: tmblarge
  type(energy_terms),intent(in) :: energs_base
  logical,intent(out) :: reduce_conf, fix_supportfunctions
  integer, intent(in) :: nit_precond, target_function, correction_orthoconstraint, nit_basis
  real(kind=8),intent(in) :: deltaenergy_multiplier_TMBexit, deltaenergy_multiplier_TMBfix
 
  ! Local variables
  real(kind=8) :: fnrmMax, meanAlpha, ediff, noise, alpha_max, delta_energy, delta_energy_prev
  integer :: iorb, istat,ierr,it,iall,it_tot, ncount
  real(kind=8),dimension(:),allocatable :: alpha,fnrmOldArr,alphaDIIS, hpsit_c_tmp, hpsit_f_tmp, hpsi_noconf, psidiff
  real(kind=8),dimension(:),allocatable :: hpsi_noprecond
  real(kind=8),dimension(:,:),allocatable :: ovrlp, coeff_old, kernel
  logical :: energy_increased, overlap_calculated
  character(len=*),parameter :: subname='getLocalizedBasis'
  real(kind=8),dimension(:),pointer :: lhphiold, lphiold, hpsit_c, hpsit_f
  type(energy_terms) :: energs
  real(8),dimension(2):: reducearr
  real(gp) :: econf


  ! Allocate all local arrays.
  call allocateLocalArrays()

  ! setting lhphiold to zero for calculate_energy_and_gradient_linear - why is this needed?
  call to_zero(max(tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp),lhphiold(1))

  call timing(iproc,'getlocbasinit','ON')
  tmb%can_use_transposed=.false.
  if(iproc==0) write(*,'(1x,a)') '======================== Creation of the basis functions... ========================'

  alpha=ldiis%alphaSD
  alphaDIIS=ldiis%alphaDIIS
  ldiis%resetDIIS=.false.
  ldiis%immediateSwitchToSD=.false.
  noise=0.d0
 
  call timing(iproc,'getlocbasinit','OF')

  overlap_calculated=.false.
  it=0
  it_tot=0
  call local_potential_dimensions(tmblarge%lzd,tmblarge%orbs,denspot%dpbox%ngatherarr(0,1))
  call start_onesided_communication(iproc, nproc, max(denspot%dpbox%ndimpot,1), denspot%rhov, &
       tmblarge%comgp%nrecvbuf, tmblarge%comgp%recvbuf, tmblarge%comgp, tmblarge%lzd)

  reduce_conf=.false.
  fix_supportfunctions=.false.
  delta_energy_prev=1.d100

  iterLoop: do
      it=it+1
      it=max(it,1) !since it could become negative (2 is subtracted if the loop cycles)
      it_tot=it_tot+1

      fnrmMax=0.d0
      fnrm=0.d0
  
      if (iproc==0) then
          write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
      endif

      ! Calculate the unconstrained gradient by applying the Hamiltonian.
      if (tmblarge%orbs%npsidim_orbs > 0)  call to_zero(tmblarge%orbs%npsidim_orbs,tmblarge%hpsi(1))
      call small_to_large_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmb%orbs, tmblarge%orbs, &
           tmb%psi, tmblarge%psi)

      call NonLocalHamiltonianApplication(iproc,at,tmblarge%orbs,rxyz,&
           proj,tmblarge%lzd,nlpspd,tmblarge%psi,tmblarge%hpsi,energs%eproj)
      ! only kinetic because waiting for communications
      call LocalHamiltonianApplication(iproc,nproc,at,tmblarge%orbs,&
           tmblarge%lzd,tmblarge%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,tmblarge%psi,tmblarge%hpsi,&
           energs,SIC,GPU,3,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmblarge%comgp)
      call full_local_potential(iproc,nproc,tmblarge%orbs,tmblarge%lzd,2,denspot%dpbox,denspot%rhov,denspot%pot_work, &
           tmblarge%comgp)
      ! only potential
      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          call vcopy(tmblarge%orbs%npsidim_orbs, tmblarge%hpsi(1), 1, hpsi_noconf(1), 1)
          call LocalHamiltonianApplication(iproc,nproc,at,tmblarge%orbs,&
               tmblarge%lzd,tmblarge%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,tmblarge%psi,tmblarge%hpsi,&
               energs,SIC,GPU,2,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmblarge%comgp,&
               hpsi_noconf=hpsi_noconf,econf=econf)
          if (nproc>1) then
              call mpiallred(econf, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
          end if
      else
          call LocalHamiltonianApplication(iproc,nproc,at,tmblarge%orbs,&
               tmblarge%lzd,tmblarge%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,tmblarge%psi,tmblarge%hpsi,&
               energs,SIC,GPU,2,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmblarge%comgp)
      end if


      if (target_function==TARGET_FUNCTION_IS_HYBRID .and. iproc==0) then
          write(*,*) 'econf, econf/tmb%orbs%norb',econf, econf/tmb%orbs%norb
      end if

      call timing(iproc,'glsynchham2','ON')
      call SynchronizeHamiltonianApplication(nproc,tmblarge%orbs,tmblarge%lzd,GPU,tmblarge%hpsi,&
           energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
      call timing(iproc,'glsynchham2','OF')



      ! Apply the orthoconstraint to the gradient. This subroutine also calculates the trace trH.
      if(iproc==0) then
          write(*,'(a)', advance='no') ' Orthoconstraint... '
      end if

      call copy_orthon_data(tmb%orthpar, tmblarge%orthpar, subname)


      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          call transpose_localized(iproc, nproc, tmblarge%orbs, tmblarge%collcom, hpsi_noconf, hpsit_c, hpsit_f, tmblarge%lzd)
      else
          call transpose_localized(iproc, nproc, tmblarge%orbs, tmblarge%collcom, tmblarge%hpsi, hpsit_c, hpsit_f, tmblarge%lzd)
      end if

      ncount=sum(tmblarge%collcom%nrecvcounts_c)
      if(ncount>0) call dcopy(ncount, hpsit_c(1), 1, hpsit_c_tmp(1), 1)
      ncount=7*sum(tmblarge%collcom%nrecvcounts_f)
      if(ncount>0) call dcopy(ncount, hpsit_f(1), 1, hpsit_f_tmp(1), 1)

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
         call calculate_energy_and_gradient_linear(iproc, nproc, it, tmb%linmat%denskern%matrix_compr, &
              ldiis, fnrmOldArr, alpha, trH, trH_old, fnrm, fnrmMax, meanAlpha, alpha_max, &
              energy_increased, tmb, lhphiold, tmblarge, overlap_calculated, energs_base, hpsit_c, hpsit_f, &    
              nit_precond, target_function, correction_orthoconstraint, hpsi_noprecond)
      else
         call calculate_energy_and_gradient_linear(iproc, nproc, it, tmb%linmat%denskern%matrix_compr, &
              ldiis, fnrmOldArr, alpha, trH, trH_old, fnrm, fnrmMax, meanAlpha, alpha_max, &
              energy_increased, tmb, lhphiold, tmblarge, overlap_calculated, energs_base, hpsit_c, hpsit_f, &    
              nit_precond, target_function, correction_orthoconstraint)
      end if


      ediff=trH-trH_old

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          if (iproc==0) write(*,*) 'ediff, delta_energy_prev', ediff, delta_energy_prev
          if ((ediff>deltaenergy_multiplier_TMBexit*delta_energy_prev .or. energy_increased) .and. it>1) then
              if (iproc==0) write(*,*) 'reduce the confinement'
              reduce_conf=.true.
          end if
      end if


      if ((ediff>deltaenergy_multiplier_TMBfix*delta_energy_prev .and. .not.energy_increased) .and. it>1 .and. &
          target_function==TARGET_FUNCTION_IS_HYBRID) then
          if (iproc==0) write(*,*) 'Will fix the support functions'
          fix_supportfunctions=.true.
      end if

      !!delta_energy_prev=delta_energy

      if (energy_increased) then
          if (iproc==0) write(*,*) 'WARNING: ENERGY INCREASED'
          tmblarge%can_use_transposed=.false.
          call dcopy(tmb%orbs%npsidim_orbs, lphiold(1), 1, tmb%psi(1), 1)
          if (scf_mode/=LINEAR_FOE) then
              ! Recalculate the kernel with the old coefficients
              call dcopy(orbs%norb*tmb%orbs%norb, coeff_old(1,1), 1, tmb%wfnmd%coeff(1,1), 1)
              allocate(kernel(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
              call memocc(istat, kernel, 'kernel', subname)
              call calculate_density_kernel(iproc, nproc, .true., orbs, tmb%orbs, &
                   tmb%wfnmd%coeff, kernel)
              call compress_matrix_for_allreduce(tmblarge%orbs%norb, tmblarge%sparsemat, &
                   kernel, tmb%linmat%denskern%matrix_compr)
              iall=-product(shape(kernel))*kind(kernel)
              deallocate(kernel, stat=istat)
              call memocc(istat, iall, 'kernel', subname)
          end if
          trH_old=0.d0
          it=it-2 !go back one iteration (minus 2 since the counter was increased)
          if(associated(tmblarge%psit_c)) then
              iall=-product(shape(tmblarge%psit_c))*kind(tmblarge%psit_c)
              deallocate(tmblarge%psit_c, stat=istat)
              call memocc(istat, iall, 'tmblarge%psit_c', subname)
          end if
          if(associated(tmblarge%psit_f)) then
              iall=-product(shape(tmblarge%psit_f))*kind(tmblarge%psit_f)
              deallocate(tmblarge%psit_f, stat=istat)
              call memocc(istat, iall, 'tmblarge%psit_f', subname)
              tmblarge%can_use_transposed=.false.
          end if
          if(iproc==0) write(*,*) 'it_tot',it_tot
          overlap_calculated=.false.
          ! print info here anyway for debugging
          if (iproc==0) write(*,'(1x,a,i6,2es15.7,f17.10,2es13.4)') 'iter, fnrm, fnrmMax, ebs, diff, noise level', &
          it, fnrm, fnrmMax, trH, ediff,noise
          if(it_tot<3*nit_basis) cycle
      end if 


      ! Write some information to the screen.
      if(iproc==0 .and. target_function==TARGET_FUNCTION_IS_TRACE) &
          write(*,'(1x,a,i6,2es15.7,f17.10,es13.4)') 'iter, fnrm, fnrmMax, trace, diff', &
          it, fnrm, fnrmMax, trH, ediff
      if(iproc==0 .and. target_function==TARGET_FUNCTION_IS_ENERGY) &
          write(*,'(1x,a,i6,2es15.7,f17.10,es13.4)') 'iter, fnrm, fnrmMax, ebs, diff', &
          it, fnrm, fnrmMax, trH, ediff
      if(iproc==0 .and. target_function==TARGET_FUNCTION_IS_HYBRID) &
          write(*,'(1x,a,i6,2es15.7,f17.10,es13.4)') 'iter, fnrm, fnrmMax, hybrid, diff', &
          it, fnrm, fnrmMax, trH, ediff
      if(it>=nit_basis .or. it_tot>=3*nit_basis .or. reduce_conf) then
          if(it>=nit_basis .and. .not.energy_increased) then
              if(iproc==0) write(*,'(1x,a,i0,a)') 'WARNING: not converged within ', it, &
                  ' iterations! Exiting loop due to limitations of iterations.'
              if(iproc==0 .and. target_function==TARGET_FUNCTION_IS_TRACE) &
                  write(*,'(1x,a,2es15.7,f15.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
              if(iproc==0 .and. target_function==TARGET_FUNCTION_IS_ENERGY) &
                  write(*,'(1x,a,2es15.7,f15.7)') 'Final values for fnrm, fnrmMax, ebs: ', fnrm, fnrmMax, trH
              infoBasisFunctions=0
          else if(it_tot>=3*nit_basis) then
              if(iproc==0) write(*,'(1x,a,i0,a)') 'WARNING: there seem to be some problems, exiting now...'
              if(iproc==0 .and. target_function==TARGET_FUNCTION_IS_TRACE) &
                  write(*,'(1x,a,2es15.7,f15.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
              if(iproc==0 .and. target_function==TARGET_FUNCTION_IS_ENERGY) &
                  write(*,'(1x,a,2es15.7,f15.7)') 'Final values for fnrm, fnrmMax, ebs: ', fnrm, fnrmMax, trH
              infoBasisFunctions=-1
          else if (reduce_conf) then
              if (iproc==0) then
                  write(*,'(1x,a,2es15.7,f15.7)') 'Final values for fnrm, fnrmMax, hybrid: ', fnrm, fnrmMax, trH
              end if
              infoBasisFunctions=0
          end if
          if(iproc==0) write(*,'(1x,a)') '============================= Basis functions created. ============================='
          if (infoBasisFunctions>=0) then
              ! Calculate the Hamiltonian matrix, since we have all quantities ready. This matrix can then be used in the first
              ! iteration of get_coeff.
              call calculate_overlap_transposed(iproc, nproc, tmblarge%orbs, tmblarge%sparsemat, tmblarge%collcom, &
                   tmblarge%psit_c, hpsit_c_tmp, tmblarge%psit_f, hpsit_f_tmp, tmb%linmat%ham%matrix_compr)
          end if

          exit iterLoop
      end if
      trH_old=trH

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
         call hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, tmblarge, &
              lphiold, alpha, trH, meanAlpha, alpha_max, alphaDIIS, psidiff)
      else
         call hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, tmblarge, &
              lphiold, alpha, trH, meanAlpha, alpha_max, alphaDIIS)
      end if


      overlap_calculated=.false.
      ! It is now not possible to use the transposed quantities, since they have changed.
      if(tmblarge%can_use_transposed) then
          iall=-product(shape(tmblarge%psit_c))*kind(tmblarge%psit_c)
          deallocate(tmblarge%psit_c, stat=istat)
          call memocc(istat, iall, 'tmblarge%psit_c', subname)
          iall=-product(shape(tmblarge%psit_f))*kind(tmblarge%psit_f)
          deallocate(tmblarge%psit_f, stat=istat)
          call memocc(istat, iall, 'tmblarge%psit_f', subname)
          tmblarge%can_use_transposed=.false.
      end if


      ! Estimate the energy change, that is to be expected in the next optimization
      ! step, given by the product of the force and the "displacement" .
      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          call estimate_energy_change(tmb%orbs, tmb%lzd, psidiff, hpsi_noprecond, delta_energy)
          if (iproc==0) write(*,*) 'delta_energy', delta_energy
          delta_energy_prev=delta_energy
      end if



      ! Copy the coefficients to coeff_old. The coefficients will be modified in reconstruct_kernel.
      if (scf_mode/=LINEAR_FOE) then
          call dcopy(orbs%norb*tmb%orbs%norb, tmb%wfnmd%coeff(1,1), 1, coeff_old(1,1), 1)
      end if

      if(scf_mode/=LINEAR_FOE) then
          allocate(ovrlp(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
          call memocc(istat, ovrlp, 'ovrlp', subname)
          call reconstruct_kernel(iproc, nproc, 1, tmb%orthpar%blocksize_pdsyev, tmb%orthpar%blocksize_pdgemm, &
               orbs, tmb, tmblarge, ovrlp, overlap_calculated, tmb%linmat%denskern%matrix_compr)
          iall=-product(shape(ovrlp))*kind(ovrlp)
          deallocate(ovrlp, stat=istat)
          call memocc(istat, iall, 'ovrlp', subname)
      end if
      if(iproc==0) then
          write(*,'(a)') 'done.'
      end if


  end do iterLoop

  ! Deallocate potential
  iall=-product(shape(denspot%pot_work))*kind(denspot%pot_work)
  deallocate(denspot%pot_work, stat=istat)
  call memocc(istat, iall, 'denspot%pot_work', subname)


  ! Keep the values for the next iteration
  reducearr(1)=0.d0
  reducearr(2)=0.d0
  do iorb=1,tmb%orbs%norbp
      reducearr(1)=reducearr(1)+alpha(iorb)
      reducearr(2)=reducearr(2)+alphaDIIS(iorb)
  end do
  call mpiallred(reducearr(1), 2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  reducearr(1)=reducearr(1)/dble(tmb%orbs%norb)
  reducearr(2)=reducearr(2)/dble(tmb%orbs%norb)

  ldiis%alphaSD=reducearr(1)
  ldiis%alphaDIIS=reducearr(2)


  ! Deallocate all local arrays.
  call deallocateLocalArrays()

contains



    subroutine allocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine allocates all local arrays.
    !
      allocate(alpha(tmb%orbs%norbp), stat=istat)
      call memocc(istat, alpha, 'alpha', subname)

      allocate(alphaDIIS(tmb%orbs%norbp), stat=istat)
      call memocc(istat, alphaDIIS, 'alphaDIIS', subname)

      allocate(fnrmOldArr(tmb%orbs%norb), stat=istat)
      call memocc(istat, fnrmOldArr, 'fnrmOldArr', subname)

      allocate(tmb%hpsi(max(tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp)), stat=istat)
      call memocc(istat, tmb%hpsi, 'tmb%hpsi', subname)
    
      allocate(lhphiold(max(tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp)), stat=istat)
      call memocc(istat, lhphiold, 'lhphiold', subname)

      allocate(lphiold(size(tmb%psi)), stat=istat)
      call memocc(istat, lphiold, 'lphiold', subname)

      allocate(hpsit_c(sum(tmblarge%collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, hpsit_c, 'hpsit_c', subname)

      allocate(hpsit_f(7*sum(tmblarge%collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, hpsit_f, 'hpsit_f', subname)

      allocate(hpsit_c_tmp(sum(tmblarge%collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, hpsit_c_tmp, 'hpsit_c_tmp', subname)

      allocate(hpsit_f_tmp(7*sum(tmblarge%collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, hpsit_f_tmp, 'hpsit_f_tmp', subname)

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
         allocate(hpsi_noconf(tmblarge%orbs%npsidim_orbs), stat=istat)
         call memocc(istat, hpsi_noconf, 'hpsi_noconf', subname)

         allocate(psidiff(tmb%orbs%npsidim_orbs), stat=istat)
         call memocc(istat, psidiff, 'psidiff', subname)

         allocate(hpsi_noprecond(tmb%orbs%npsidim_orbs), stat=istat)
         call memocc(istat, hpsi_noprecond, 'hpsi_noprecond', subname)
      end if

      if (scf_mode/=LINEAR_FOE) then
          allocate(coeff_old(tmb%orbs%norb,orbs%norb), stat=istat)
          call memocc(istat, coeff_old, 'coeff_old', subname)
      end if


    end subroutine allocateLocalArrays


    subroutine deallocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine deallocates all local arrays.
    !
      iall=-product(shape(alpha))*kind(alpha)
      deallocate(alpha, stat=istat)
      call memocc(istat, iall, 'alpha', subname)

      iall=-product(shape(alphaDIIS))*kind(alphaDIIS)
      deallocate(alphaDIIS, stat=istat)
      call memocc(istat, iall, 'alphaDIIS', subname)

      iall=-product(shape(fnrmOldArr))*kind(fnrmOldArr)
      deallocate(fnrmOldArr, stat=istat)
      call memocc(istat, iall, 'fnrmOldArr', subname)

      iall=-product(shape(tmb%hpsi))*kind(tmb%hpsi)
      deallocate(tmb%hpsi, stat=istat)
      call memocc(istat, iall, 'tmb%hpsi', subname)

      iall=-product(shape(lhphiold))*kind(lhphiold)
      deallocate(lhphiold, stat=istat)
      call memocc(istat, iall, 'lhphiold', subname)

      iall=-product(shape(lphiold))*kind(lphiold)
      deallocate(lphiold, stat=istat)
      call memocc(istat, iall, 'lphiold', subname)

      iall=-product(shape(hpsit_c))*kind(hpsit_c)
      deallocate(hpsit_c, stat=istat)
      call memocc(istat, iall, 'hpsit_c', subname)

      iall=-product(shape(hpsit_f))*kind(hpsit_f)
      deallocate(hpsit_f, stat=istat)
      call memocc(istat, iall, 'hpsit_f', subname)

      iall=-product(shape(hpsit_c_tmp))*kind(hpsit_c_tmp)
      deallocate(hpsit_c_tmp, stat=istat)
      call memocc(istat, iall, 'hpsit_c_tmp', subname)

      iall=-product(shape(hpsit_f_tmp))*kind(hpsit_f_tmp)
      deallocate(hpsit_f_tmp, stat=istat)
      call memocc(istat, iall, 'hpsit_f_tmp', subname)

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
         iall=-product(shape(hpsi_noconf))*kind(hpsi_noconf)
         deallocate(hpsi_noconf, stat=istat)
         call memocc(istat, iall, 'hpsi_noconf', subname)

         iall=-product(shape(psidiff))*kind(psidiff)
         deallocate(psidiff, stat=istat)
         call memocc(istat, iall, 'psidiff', subname)

         iall=-product(shape(hpsi_noprecond))*kind(hpsi_noprecond)
         deallocate(hpsi_noprecond, stat=istat)
         call memocc(istat, iall, 'hpsi_noprecond', subname)
      end if

      if (scf_mode/=LINEAR_FOE) then
          iall=-product(shape(coeff_old))*kind(coeff_old)
          deallocate(coeff_old, stat=istat)
          call memocc(istat, iall, 'coeff_old', subname)
      end if

    end subroutine deallocateLocalArrays


end subroutine getLocalizedBasis



subroutine improveOrbitals(iproc, tmb, ldiis, alpha)
  use module_base
  use module_types
  use module_interfaces, except_this_one => improveOrbitals
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc
  type(DFT_wavefunction),intent(inout) :: tmb
  type(localizedDIISParameters),intent(inout) :: ldiis
  real(kind=8),dimension(tmb%orbs%norbp),intent(in) :: alpha
  
  ! Local variables
  integer :: istart, iorb, iiorb, ilr, ncount
  

  if(ldiis%isx==0) then ! steepest descents
      call timing(iproc,'optimize_SD   ','ON')
      istart=1
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          call daxpy(ncount, -alpha(iorb), tmb%hpsi(istart), 1, tmb%psi(istart), 1)
          istart=istart+ncount
      end do
      call timing(iproc,'optimize_SD   ','OF')
  else! DIIS
      ldiis%mis=mod(ldiis%is,ldiis%isx)+1
      ldiis%is=ldiis%is+1
      if(ldiis%alphaDIIS/=1.d0) then
          call dscal(max(tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp), ldiis%alphaDIIS, tmb%hpsi, 1)
      end if
      call optimizeDIIS(iproc, tmb%orbs, tmb%orbs, tmb%lzd, tmb%hpsi, tmb%psi, ldiis)
  end if

end subroutine improveOrbitals


subroutine my_geocode_buffers(geocode,nl1,nl2,nl3)
  implicit none
  integer, intent(out) :: nl1,nl2,nl3
  character(len=1), intent(in) :: geocode
  !local variables
  logical :: perx,pery,perz
  integer :: nr1,nr2,nr3

  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nl1,nr1)
  call ext_buffers(pery,nl2,nr2)
  call ext_buffers(perz,nl3,nr3)

end subroutine my_geocode_buffers



subroutine diagonalizeHamiltonian2(iproc, orbs, HamSmall, ovrlp, eval)
  !
  ! Purpose:
  ! ========
  !   Diagonalizes the Hamiltonian HamSmall and makes sure that all MPI processes give
  !   the same result. This is done by requiring that the first entry of each vector
  !   is positive.
  !
  ! Calling arguments:
  ! ==================
  !   Input arguments:
  !   ----------------
  !     iproc     process ID
  !     nproc     number of MPI processes
  !     orbs      type describing the physical orbitals psi
  !   Input / Putput arguments
  !     HamSmall  on input: the Hamiltonian
  !               on exit: the eigenvectors
  !   Output arguments
  !     eval      the associated eigenvalues 
  !
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc
  type(orbitals_data), intent(inout) :: orbs
  real(kind=8),dimension(orbs%norb, orbs%norb),intent(inout) :: HamSmall
  real(kind=8),dimension(orbs%norb, orbs%norb),intent(in) :: ovrlp
  real(kind=8),dimension(orbs%norb),intent(out) :: eval

  ! Local variables
  integer :: lwork, info, istat, iall
  real(kind=8),dimension(:),allocatable :: work
  character(len=*),parameter :: subname='diagonalizeHamiltonian'

  call timing(iproc,'diagonal_seq  ','ON')

  ! DEBUG: print hamiltonian and overlap matrices
  !if (iproc==0) then
  !   open(10)
  !   open(11)
  !   do iorb=1,orbs%norb
  !      do jorb=1,orbs%norb
  !         write(10,*) iorb,jorb,HamSmall(iorb,jorb)
  !         write(11,*) iorb,jorb,ovrlp(iorb,jorb)
  !      end do
  !      write(10,*) ''
  !      write(11,*) ''
  !   end do
  !   close(10)
  !   close(11)
  !end if
  ! DEBUG: print hamiltonian and overlap matrices

  ! Get the optimal work array size
  lwork=-1 
  allocate(work(100), stat=istat)
  call memocc(istat, work, 'work', subname)
  call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 
  lwork=int(work(1))

  ! Deallocate the work array and reallocate it with the optimal size
  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  call memocc(istat, iall, 'work', subname)
  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
  call memocc(istat, work, 'work', subname)

  ! Diagonalize the Hamiltonian
  call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 
 
  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  call memocc(istat, iall, 'work', subname)

  call timing(iproc,'diagonal_seq  ','OF')

end subroutine diagonalizeHamiltonian2

subroutine small_to_large_locreg(iproc, nproc, lzdsmall, lzdlarge, orbssmall, orbslarge, phismall, philarge)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzdsmall, lzdlarge
  type(orbitals_data),intent(in) :: orbssmall, orbslarge
  real(kind=8),dimension(orbssmall%npsidim_orbs),intent(in) :: phismall
  real(kind=8),dimension(orbslarge%npsidim_orbs),intent(out) :: philarge
  
  ! Local variables
  integer :: ists, istl, iorb, ilr, ilrlarge, sdim, ldim, nspin
       call timing(iproc,'small2large','ON') ! lr408t 
  ! No need to put arrays to zero, Lpsi_to_global2 will handle this.
  call to_zero(orbslarge%npsidim_orbs, philarge(1))
  ists=1
  istl=1
  do iorb=1,orbslarge%norbp
      ilr = orbssmall%inWhichLocreg(orbssmall%isorb+iorb)
      ilrlarge = orbslarge%inWhichLocreg(orbslarge%isorb+iorb)
      sdim=lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
      ldim=lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
      nspin=1 !this must be modified later
      call Lpsi_to_global2(iproc, sdim, ldim, orbssmall%norb, orbssmall%nspinor, nspin, lzdlarge%llr(ilrlarge), &
           lzdsmall%llr(ilr), phismall(ists), philarge(istl))
      ists=ists+lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
      istl=istl+lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
  end do
  if(orbssmall%norbp>0 .and. ists/=orbssmall%npsidim_orbs+1) then
      write(*,'(3(a,i0))') 'ERROR on process ',iproc,': ',ists,'=ists /= orbssmall%npsidim_orbs+1=',orbssmall%npsidim_orbs+1
      stop
  end if
  if(orbslarge%norbp>0 .and. istl/=orbslarge%npsidim_orbs+1) then
      write(*,'(3(a,i0))') 'ERROR on process ',iproc,': ',istl,'=istk /= orbslarge%npsidim_orbs+1=',orbslarge%npsidim_orbs+1
      stop
  end if
       call timing(iproc,'small2large','OF') ! lr408t 
end subroutine small_to_large_locreg


subroutine large_to_small_locreg(iproc, nproc, lzdsmall, lzdlarge, orbssmall, orbslarge, philarge, phismall)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzdsmall, lzdlarge
  type(orbitals_data),intent(in) :: orbssmall, orbslarge
  real(kind=8),dimension(orbslarge%npsidim_orbs),intent(in) :: philarge
  real(kind=8),dimension(orbssmall%npsidim_orbs),intent(out) :: phismall
  
  ! Local variables
  integer :: istl, ists, ilr, ilrlarge, ldim, gdim, iorb
       call timing(iproc,'large2small','ON') ! lr408t   
  ! Transform back to small locreg
  ! No need to this array to zero, since all values will be filled with a value during the copy.
  !!call to_zero(orbssmall%npsidim_orbs, phismall(1))
  ists=1
  istl=1
  do iorb=1,orbssmall%norbp
      ilr = orbssmall%inWhichLocreg(orbssmall%isorb+iorb)
      ilrlarge = orbslarge%inWhichLocreg(orbslarge%isorb+iorb)
      ldim=lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
      gdim=lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
      call psi_to_locreg2(iproc, ldim, gdim, lzdsmall%llr(ilr), lzdlarge%llr(ilrlarge), &
           philarge(istl:istl+gdim-1), phismall(ists:ists+ldim-1))
      ists=ists+lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
      istl=istl+lzdlarge%llr(ilrlarge)%wfd%nvctr_c+7*lzdlarge%llr(ilrlarge)%wfd%nvctr_f
  end do

  if(orbssmall%norbp>0 .and. ists/=orbssmall%npsidim_orbs+1) stop 'ists/=orbssmall%npsidim_orbs+1'
  if(orbslarge%norbp>0 .and. istl/=orbslarge%npsidim_orbs+1) stop 'istl/=orbslarge%npsidim_orbs+1'
       call timing(iproc,'large2small','OF') ! lr408t 
end subroutine large_to_small_locreg







subroutine communicate_basis_for_density_collective(iproc, nproc, lzd, orbs, lphi, collcom_sr)
  use module_base
  use module_types
  use module_interfaces, except_this_one => communicate_basis_for_density_collective
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  real(kind=8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in) :: lphi
  type(collective_comms),intent(inout) :: collcom_sr
  
  ! Local variables
  integer :: ist, istr, iorb, iiorb, ilr, istat, iall
  real(kind=8),dimension(:),allocatable :: psir, psirwork, psirtwork
  type(workarr_sumrho) :: w
  character(len=*),parameter :: subname='communicate_basis_for_density_collective'

  call timing(iproc,'commbasis4dens','ON')

  allocate(psir(collcom_sr%ndimpsi_c), stat=istat)
  call memocc(istat, psir, 'psir', subname)

  ! Allocate the communication buffers for the calculation of the charge density.
  !call allocateCommunicationbufferSumrho(iproc, comsr, subname)
  ! Transform all orbitals to real space.
  ist=1
  istr=1
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inWhichLocreg(iiorb)
      call initialize_work_arrays_sumrho(lzd%Llr(ilr), w)
      call daub_to_isf(lzd%Llr(ilr), w, lphi(ist), psir(istr))
      call deallocate_work_arrays_sumrho(w)
      ist = ist + lzd%Llr(ilr)%wfd%nvctr_c + 7*lzd%Llr(ilr)%wfd%nvctr_f
      istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
  end do
  if(istr/=collcom_sr%ndimpsi_c+1) then
      write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=collcom_sr%ndimpsi_c+1'
      stop
  end if

  allocate(psirwork(collcom_sr%ndimpsi_c), stat=istat)
  call memocc(istat, psirwork, 'psirwork', subname)

  call transpose_switch_psir(orbs, collcom_sr, psir, psirwork)

  iall=-product(shape(psir))*kind(psir)
  deallocate(psir, stat=istat)
  call memocc(istat, iall, 'psir', subname)

  allocate(psirtwork(collcom_sr%ndimind_c), stat=istat)
  call memocc(istat, psirtwork, 'psirtwork', subname)

  call transpose_communicate_psir(iproc, nproc, collcom_sr, psirwork, psirtwork)

  iall=-product(shape(psirwork))*kind(psirwork)
  deallocate(psirwork, stat=istat)
  call memocc(istat, iall, 'psirwork', subname)

  call transpose_unswitch_psirt(collcom_sr, psirtwork, collcom_sr%psit_c)

  iall=-product(shape(psirtwork))*kind(psirtwork)
  deallocate(psirtwork, stat=istat)
  call memocc(istat, iall, 'psirtwork', subname)

  call timing(iproc,'commbasis4dens','OF')

end subroutine communicate_basis_for_density_collective





subroutine DIISorSD(iproc, it, trH, tmbopt, ldiis, alpha, alphaDIIS, lphioldopt)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, it
  real(kind=8),intent(in) :: trH
  type(DFT_wavefunction),intent(inout) :: tmbopt
  type(localizedDIISParameters),intent(inout) :: ldiis
  real(kind=8),dimension(tmbopt%orbs%norbp),intent(inout) :: alpha, alphaDIIS
  real(kind=8),dimension(max(tmbopt%orbs%npsidim_comp,tmbopt%orbs%npsidim_orbs)),intent(out) :: lphioldopt
  
  ! Local variables
  integer :: idsx, ii, offset, istdest, iorb, iiorb, ilr, ncount, istsource
  

  ! Purpose:
  ! ========
  !   This subroutine decides whether one should use DIIS or variable step size
  !   steepest descent to improve the orbitals. In the beginning we start with DIIS
  !   with history length lin%DIISHistMax. If DIIS becomes unstable, we switch to
  !   steepest descent. If the steepest descent iterations are successful, we switch
  !   back to DIIS, but decrease the DIIS history length by one. However the DIIS
  !   history length is limited to be larger or equal than lin%DIISHistMin.


  ! If we swicthed to SD in the previous iteration, reset this flag.
  if(ldiis%switchSD) ldiis%switchSD=.false.
  !if(iproc==0) write(*,'(a,2es15.6,l5)') 'trH, ldiis%trmin, ldiis%resetDIIS', trH, ldiis%trmin, ldiis%resetDIIS

  ! Now come some checks whether the trace is descreasing or not. This further decides
  ! whether we should use DIIS or SD.

  ! Determine wheter the trace is decreasing (as it should) or increasing.
  ! This is done by comparing the current value with diisLIN%energy_min, which is
  ! the minimal value of the trace so far.
  !if(iproc==0) write(*,*) 'trH, ldiis%trmin', trH, ldiis%trmin
  if(trH<=ldiis%trmin .and. .not.ldiis%resetDIIS) then
      ! Everything ok
      ldiis%trmin=trH
      ldiis%switchSD=.false.
      ldiis%itBest=it
      ldiis%icountSDSatur=ldiis%icountSDSatur+1
      ldiis%icountDIISFailureCons=0
      !if(iproc==0) write(*,*) 'everything ok, copy last psi...'
      call dcopy(size(tmbopt%psi), tmbopt%psi(1), 1, lphioldopt(1), 1)

      ! If we are using SD (i.e. diisLIN%idsx==0) and the trace has been decreasing
      ! for at least 10 iterations, switch to DIIS. However the history length is decreased.
      if(ldiis%icountSDSatur>=10 .and. ldiis%isx==0 .or. ldiis%immediateSwitchToSD) then
          ldiis%icountSwitch=ldiis%icountSwitch+1
          idsx=max(ldiis%DIISHistMin,ldiis%DIISHistMax-ldiis%icountSwitch)
          if(idsx>0) then
              if(iproc==0) write(*,'(1x,a,i0)') 'switch to DIIS with new history length ', idsx
              ldiis%icountSDSatur=0
              ldiis%icountSwitch=0
              ldiis%icountDIISFailureTot=0
              ldiis%icountDIISFailureCons=0
              ldiis%is=0
              ldiis%switchSD=.false.
              ldiis%trmin=1.d100
              ldiis%trold=1.d100
              alpha=ldiis%alphaSD
              alphaDIIS=ldiis%alphaDIIS
              ldiis%icountDIISFailureTot=0
              ldiis%icountDIISFailureCons=0
              ldiis%immediateSwitchToSD=.false.
          end if
      end if
  else
      ! The trace is growing.
      ! Count how many times this occurs and (if we are using DIIS) switch to SD after 3 
      ! total failures or after 2 consecutive failures.
      ldiis%icountDIISFailureCons=ldiis%icountDIISFailureCons+1
      ldiis%icountDIISFailureTot=ldiis%icountDIISFailureTot+1
      ldiis%icountSDSatur=0
      if((ldiis%icountDIISFailureCons>=2 .or. ldiis%icountDIISFailureTot>=3 .or. ldiis%resetDIIS) .and. ldiis%isx>0) then
          ! Switch back to SD.
          alpha=ldiis%alphaSD
          if(iproc==0) then
              if(ldiis%icountDIISFailureCons>=2) write(*,'(1x,a,i0,a,es10.3)') 'DIIS failed ', &
                  ldiis%icountDIISFailureCons, ' times consecutively. Switch to SD with stepsize', alpha(1)
              if(ldiis%icountDIISFailureTot>=3) write(*,'(1x,a,i0,a,es10.3)') 'DIIS failed ', &
                  ldiis%icountDIISFailureTot, ' times in total. Switch to SD with stepsize', alpha(1)
              if(ldiis%resetDIIS) write(*,'(1x,a)') 'reset DIIS due to flag'
          end if
          if(ldiis%resetDIIS) then
              ldiis%resetDIIS=.false.
              ldiis%immediateSwitchToSD=.true.
              ldiis%trmin=1.d100
          end if
          ! Otherwise there could be problems due to the orthonormalization (which sligtly increases 
          ! value of the target function)
          ldiis%trmin=1.d100
          ! Try to get back the orbitals of the best iteration. This is possible if
          ! these orbitals are still present in the DIIS history.
          if(it-ldiis%itBest<ldiis%isx) then
             if(iproc==0) then
                 if(iproc==0) write(*,'(1x,a,i0,a)')  'Recover the orbitals from iteration ', &
                     ldiis%itBest, ' which are the best so far.'
             end if
             ii=modulo(ldiis%mis-(it-ldiis%itBest),ldiis%mis)
             offset=0
             istdest=1
             !if(iproc==0) write(*,*) 'copy DIIS history psi...'
             do iorb=1,tmbopt%orbs%norbp
                 iiorb=tmbopt%orbs%isorb+iorb
                 ilr=tmbopt%orbs%inWhichLocreg(iiorb)
                 ncount=tmbopt%lzd%llr(ilr)%wfd%nvctr_c+7*tmbopt%lzd%llr(ilr)%wfd%nvctr_f
                 istsource=offset+ii*ncount+1
                 call dcopy(ncount, ldiis%phiHist(istsource), 1, tmbopt%psi(istdest), 1)
                 call dcopy(ncount, ldiis%phiHist(istsource), 1, lphioldopt(istdest), 1)
                 offset=offset+ldiis%isx*ncount
                 istdest=istdest+ncount
             end do
         else
             !if(iproc==0) write(*,*) 'copy last psi...'
             call dcopy(size(tmbopt%psi), tmbopt%psi(1), 1, lphioldopt(1), 1)
         end if
         ldiis%isx=0
         ldiis%switchSD=.true.
      end if
      ! to indicate that no orthonormalization is required... (CHECK THIS!)
      if(ldiis%isx==0) ldiis%switchSD=.true. 
  end if

end subroutine DIISorSD




subroutine reconstruct_kernel(iproc, nproc, iorder, blocksize_dsyev, blocksize_pdgemm, orbs, tmb, &
           tmblarge, ovrlp_tmb, overlap_calculated, kernel_compr)
  use module_base
  use module_types
  use module_interfaces, except_this_one => reconstruct_kernel
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, iorder, blocksize_dsyev, blocksize_pdgemm
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(inout):: tmb, tmblarge
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(out):: ovrlp_tmb
  logical,intent(out):: overlap_calculated
  real(8),dimension(tmblarge%sparsemat%nvctr),intent(out):: kernel_compr

  ! Local variables
  integer:: istat, ierr, iall
  real(8),dimension(:,:),allocatable:: coeff_tmp, ovrlp_tmp, ovrlp_coeff, kernel
  real(8),dimension(:),allocatable :: ovrlp_compr
  character(len=*),parameter:: subname='reconstruct_kernel'
  integer,parameter :: ALLGATHERV=1, ALLREDUCE=2
  integer,parameter:: communication_strategy=ALLREDUCE
  !!integer :: iorb,jorb

  call timing(iproc,'renormCoefComp','ON')

  allocate(coeff_tmp(tmb%orbs%norb,max(orbs%norb,1)), stat=istat)
  call memocc(istat, coeff_tmp, 'coeff_tmp', subname)
  allocate(ovrlp_tmp(orbs%norb,max(orbs%norbp,1)), stat=istat)
  call memocc(istat, ovrlp_tmp, 'ovrlp_tmp', subname)
  allocate(ovrlp_coeff(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)

  if(iproc==0) then
      write(*,'(a)',advance='no') 'coeff renormalization...'
  end if

  ! Calculate the overlap matrix between the TMBs.
  if(.not.tmb%can_use_transposed) then
      if(associated(tmb%psit_c)) then
          iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
          deallocate(tmb%psit_c, stat=istat)
          call memocc(istat, iall, 'tmb%psit_c', subname)
      end if
      if(associated(tmb%psit_f)) then
          iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
          deallocate(tmb%psit_f, stat=istat)
          call memocc(istat, iall, 'tmb%psit_f', subname)
      end if
      allocate(tmb%psit_c(sum(tmb%collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, tmb%psit_c, 'tmb%psit_c', subname)
      allocate(tmb%psit_f(7*sum(tmb%collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, tmb%psit_f, 'tmb%psit_f', subname)
      call transpose_localized(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
      tmb%can_use_transposed=.true.
  end if
  call timing(iproc,'renormCoefComp','OF')
  allocate(ovrlp_compr(tmblarge%sparsemat%nvctr))
  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmblarge%sparsemat, tmb%collcom, &
       tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, ovrlp_compr)
  call uncompressMatrix(tmb%orbs%norb, tmblarge%sparsemat, ovrlp_compr, ovrlp_tmb)
  deallocate(ovrlp_compr)
  call timing(iproc,'renormCoefComp','ON')
  overlap_calculated=.true.

  !write(*,*) 'iproc, orbs%isorb', iproc, orbs%isorb
  ! Calculate the overlap matrix among the coefficients with resct to ovrlp_tmb.
  if (communication_strategy==ALLGATHERV) then
      if (orbs%norbp>0 )then
          call dgemm('n', 'n', tmb%orbs%norb, orbs%norbp, tmb%orbs%norb, 1.d0, ovrlp_tmb(1,1), tmb%orbs%norb, &
               tmb%wfnmd%coeff(1,orbs%isorb+1), tmb%orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
      end if
      if (orbs%norbp>0 )then
          call dgemm('t', 'n', orbs%norb, orbs%norbp, tmb%orbs%norb, 1.d0,  tmb%wfnmd%coeff(1,1), tmb%orbs%norb, &
               coeff_tmp(1,1), tmb%orbs%norb, 0.d0, ovrlp_tmp(1,1), orbs%norb)
      end if
      call timing(iproc,'renormCoefComp','OF')
      call timing(iproc,'renormCoefComm','ON')
      ! Gather together the complete matrix
      if (nproc>1) then
         call mpi_allgatherv(ovrlp_tmp(1,1), orbs%norb*orbs%norbp, mpi_double_precision, ovrlp_coeff(1,1), &
              orbs%norb*orbs%norb_par(:,0), orbs%norb*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      else
         call vcopy(orbs%norb*orbs%norb,ovrlp_tmp(1,1),1,ovrlp_coeff(1,1),1)
      end if
      call timing(iproc,'renormCoefComm','OF')
      call timing(iproc,'renormCoefComp','ON')
  end if

  if (communication_strategy==ALLREDUCE) then
      if (tmb%orbs%norbp>0) then
          call dgemm('n', 'n', tmb%orbs%norbp, orbs%norb, tmb%orbs%norb, 1.d0, ovrlp_tmb(tmb%orbs%isorb+1,1), &
               tmb%orbs%norb, tmb%wfnmd%coeff(1,1), tmb%orbs%norb, 0.d0, coeff_tmp, tmb%orbs%norbp)
          call dgemm('t', 'n', orbs%norb, orbs%norb, tmb%orbs%norbp, 1.d0, tmb%wfnmd%coeff(tmb%orbs%isorb+1,1), &
               tmb%orbs%norb, coeff_tmp, tmb%orbs%norbp, 0.d0, ovrlp_coeff, orbs%norb)
      else
          call to_zero(orbs%norb**2, ovrlp_coeff(1,1))
      end if
      if (nproc>1) then
          call timing(iproc,'renormCoefComp','OF')
          call timing(iproc,'renormCoefComm','ON')
          call mpiallred(ovrlp_coeff(1,1), orbs%norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
          call timing(iproc,'renormCoefComm','OF')
          call timing(iproc,'renormCoefComp','ON')
      end if
  end if



  ! Recalculate the kernel.
  call overlapPowerMinusOneHalf_old(iproc, nproc, bigdft_mpi%mpi_comm, iorder, &
       blocksize_dsyev, blocksize_pdgemm, orbs%norb, orbs%norbp, orbs%isorb, ovrlp_coeff)

  ! Build the new linear combinations
  if (communication_strategy==ALLGATHERV) then
      if (orbs%norbp>0 )then
          call dgemm('n', 'n', tmb%orbs%norb, orbs%norbp, orbs%norb, 1.d0, tmb%wfnmd%coeff(1,1), tmb%orbs%norb, &
               ovrlp_coeff(1,orbs%isorb+1), orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
      end if
      call timing(iproc,'renormCoefComp','OF')
      call timing(iproc,'renormCoefComm','ON')
      if (nproc>1) then
         call mpi_allgatherv(coeff_tmp(1,1), tmb%orbs%norb*orbs%norbp, mpi_double_precision, &
              tmb%wfnmd%coeff(1,1), tmb%orbs%norb*orbs%norb_par(:,0), tmb%orbs%norb*orbs%isorb_par, &
              mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      else
         call vcopy(tmb%orbs%norb*orbs%norb,coeff_tmp(1,1),1,tmb%wfnmd%coeff(1,1),1)
      end if
      call timing(iproc,'renormCoefComm','OF')
  end if

  if (communication_strategy==ALLREDUCE) then
     if (orbs%norbp>0) then
         call dgemm('n', 'n', tmb%orbs%norb, orbs%norb, orbs%norbp, 1.d0, tmb%wfnmd%coeff(1,orbs%isorb+1), &
              tmb%orbs%norb, ovrlp_coeff(orbs%isorb+1,1), orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
     else
         call to_zero(tmb%orbs%norb*orbs%norb, coeff_tmp(1,1))
     end if
     call timing(iproc,'renormCoefComp','OF')
     call timing(iproc,'renormCoefComm','ON')
     if (nproc>1) then
         call mpi_allreduce(coeff_tmp(1,1), tmb%wfnmd%coeff(1,1), tmb%orbs%norb*orbs%norb, mpi_double_precision, &
              mpi_sum, bigdft_mpi%mpi_comm, ierr)
     else
         call vcopy(tmb%orbs%norb*orbs%norb, coeff_tmp(1,1), 1, tmb%wfnmd%coeff(1,1), 1)
     end if
      call timing(iproc,'renormCoefComm','OF')
  end if


  !call dcopy(tmb%orbs%norb*orbs%norb, coeff_tmp(1,1), 1, tmb%wfnmd%coeff(1,1), 1)

  !!! Check normalization
  !!call dgemm('n', 'n', tmb%orbs%norb, orbs%norb, tmb%orbs%norb, 1.d0, ovrlp_tmb(1,1), tmb%orbs%norb, &
  !!     tmb%wfnmd%coeff(1,1), tmb%orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
  !!do iorb=1,orbs%norb
  !!    do jorb=1,orbs%norb
  !!        tt=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, coeff_tmp(1,jorb), 1)
  !!        tt2=ddot(tmb%orbs%norb, coeff_tmp(1,iorb), 1, tmb%wfnmd%coeff(1,jorb), 1)
  !!        tt3=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, tmb%wfnmd%coeff(1,jorb), 1)
  !!        if(iproc==0) write(200,'(2i6,3es15.5)') iorb, jorb, tt, tt2, tt3
  !!    end do
  !!end do


  ! Recalculate the kernel
  allocate(kernel(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, kernel, 'kernel', subname)
  call calculate_density_kernel(iproc, nproc, .true., orbs, tmb%orbs, &
       tmb%wfnmd%coeff, kernel)

  !DEBUG LR
  !!if (iproc==0) then
  !!   open(10)
  !!   do iorb=1,tmb%orbs%norb
  !!      do jorb=1,tmb%orbs%norb
  !!         write(10,*) iorb,jorb,kernel(iorb,jorb)
  !!      end do
  !!   end do
  !!   close(10)
  !!end if
  !END DEBUG LR

  call compress_matrix_for_allreduce(tmblarge%orbs%norb, tmblarge%sparsemat, kernel, kernel_compr)
  iall=-product(shape(kernel))*kind(kernel)
  deallocate(kernel,stat=istat)
  call memocc(istat,iall,'kernel',subname)

  iall=-product(shape(coeff_tmp))*kind(coeff_tmp)
  deallocate(coeff_tmp,stat=istat)
  call memocc(istat,iall,'coeff_tmp',subname)

  iall=-product(shape(ovrlp_tmp))*kind(ovrlp_tmp)
  deallocate(ovrlp_tmp,stat=istat)
  call memocc(istat,iall,'ovrlp_tmp',subname)

  iall=-product(shape(ovrlp_coeff))*kind(ovrlp_coeff)
  deallocate(ovrlp_coeff,stat=istat)
  call memocc(istat,iall,'ovrlp_coeff',subname)


end subroutine reconstruct_kernel


!> Estimate the energy change, given by the product of the force and the "displacement" .
subroutine estimate_energy_change(orbs, lzd, psidiff, hpsi_noprecond, delta_energy)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  real(kind=8),dimension(orbs%npsidim_orbs),intent(in) :: psidiff, hpsi_noprecond
  real(kind=8),intent(out) :: delta_energy

  ! Local variables
  integer :: ist, iorb, iiorb, ilr, ncount, ierr
  real(kind=8) :: tt, ddot

  ist=1
  delta_energy=0.d0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
      tt=ddot(ncount, psidiff(ist), 1, hpsi_noprecond(ist), 1)
      delta_energy=delta_energy+4.0d0*tt
      ist=ist+ncount
  end do
  call mpiallred(delta_energy, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)

end subroutine estimate_energy_change
