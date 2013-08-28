!> @file
!! Linear version: Handle local basis set
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

subroutine get_coeff(iproc,nproc,scf_mode,orbs,at,rxyz,denspot,GPU,infoCoeff,&
    ebs,nlpspd,proj,SIC,tmb,fnrm,calculate_overlap_matrix,communicate_phi_for_lsumrho,&
    calculate_ham,ham_small,extra_states,convcrit_dmin,nitdmin,curvefit_dmin,ldiis_coeff,cdft)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => get_coeff, exceptThisOneA => writeonewave
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use constrained_dft
  use diis_sd_optimization
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, scf_mode
  type(orbitals_data),intent(inout) :: orbs
  type(atoms_data),intent(in) :: at
  real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
  type(DFT_local_fields), intent(inout) :: denspot
  type(GPU_pointers),intent(inout) :: GPU
  integer,intent(out) :: infoCoeff
  real(kind=8),intent(inout) :: ebs
  real(kind=8),intent(inout) :: fnrm
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd
  real(wp),dimension(nlpspd%nprojel),intent(inout) :: proj
  type(SIC_data),intent(in) :: SIC
  type(DFT_wavefunction),intent(inout) :: tmb
  logical,intent(in):: calculate_overlap_matrix, communicate_phi_for_lsumrho
  logical,intent(in) :: calculate_ham
  type(sparseMatrix), intent(inout) :: ham_small ! for foe only
  type(DIIS_obj),intent(inout),optional :: ldiis_coeff ! for dmin only
  integer, intent(in), optional :: nitdmin ! for dmin only
  real(kind=gp), intent(in), optional :: convcrit_dmin ! for dmin only
  logical, intent(in), optional :: curvefit_dmin ! for dmin only
  type(cdft_data),intent(inout),optional :: cdft
  integer, intent(in) :: extra_states


  ! Local variables 
  integer :: istat, iall, iorb, jorb, info, ind_ham, ind_denskern
  integer :: isegsmall, iseglarge, iismall, iilarge, i, is, ie, matrixindex_in_compressed
  real(kind=8),dimension(:),allocatable :: hpsit_c, hpsit_f
  real(kind=8),dimension(:,:,:),allocatable :: matrixElements
  type(confpot_data),dimension(:),allocatable :: confdatarrtmp
  type(energy_terms) :: energs

  character(len=*),parameter :: subname='get_coeff'
  real(kind=gp) :: tmprtr
  real(kind=gp),dimension(:,:),allocatable :: coeff_orig
  real(kind=8) :: deviation
  integer :: iat, iiorb, jjorb


  if(calculate_ham) then
      call local_potential_dimensions(tmb%ham_descr%lzd,tmb%orbs,denspot%dpbox%ngatherarr(0,1))
      call start_onesided_communication(iproc, nproc, max(denspot%dpbox%ndimpot,1), denspot%rhov, &
           tmb%ham_descr%comgp%nrecvbuf, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, tmb%ham_descr%lzd)
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
          call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
               tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
          tmb%can_use_transposed=.true.
      end if

      call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
           tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%ovrlp)
    
  end if

  ! Temporaray: check deviation from unity
  allocate(tmb%linmat%ovrlp%matrix(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  call uncompressMatrix(iproc,tmb%linmat%ovrlp)
  call deviation_from_unity(iproc, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, deviation)
  if (iproc==0) then
      write(*,'(a,es16.6)') 'max dev from unity', deviation
  end if
  deallocate(tmb%linmat%ovrlp%matrix, stat=istat)

  ! Post the p2p communications for the density. (must not be done in inputguess)
  if(communicate_phi_for_lsumrho) then
      call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, max(tmb%npsidim_orbs,tmb%npsidim_comp), &
           tmb%orbs, tmb%psi, tmb%collcom_sr)
  end if

  if(iproc==0) write(*,'(1x,a)') '----------------------------------- Determination of the orbitals in this new basis.'

  ! Calculate the Hamiltonian matrix if it is not already present.
  if(calculate_ham) then
      allocate(confdatarrtmp(tmb%orbs%norbp))
      call default_confinement_data(confdatarrtmp,tmb%orbs%norbp)

      call small_to_large_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
           tmb%orbs, tmb%psi, tmb%ham_descr%psi)

      if (tmb%ham_descr%npsidim_orbs > 0) call to_zero(tmb%ham_descr%npsidim_orbs,tmb%hpsi(1))

      call NonLocalHamiltonianApplication(iproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,rxyz,&
           proj,tmb%ham_descr%lzd,nlpspd,tmb%ham_descr%psi,tmb%hpsi,energs%eproj)
      ! only kinetic as waiting for communications
      call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
           tmb%ham_descr%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%ham_descr%psi,tmb%hpsi,&
           energs,SIC,GPU,3,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
      call full_local_potential(iproc,nproc,tmb%orbs,tmb%ham_descr%lzd,2,denspot%dpbox,denspot%rhov,denspot%pot_work, &
           tmb%ham_descr%comgp)
      !call wait_p2p_communication(iproc, nproc, tmb%ham_descr%comgp)
      ! only potential
      call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
           tmb%ham_descr%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%ham_descr%psi,tmb%hpsi,&
           energs,SIC,GPU,2,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
      call timing(iproc,'glsynchham1','ON')
      call SynchronizeHamiltonianApplication(nproc,tmb%ham_descr%npsidim_orbs,tmb%orbs,tmb%ham_descr%lzd,GPU,tmb%hpsi,&
           energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
      call timing(iproc,'glsynchham1','OF')
      deallocate(confdatarrtmp)

      !DEBUG
      !if(iproc==0) then
      !  print *,'Ekin,Epot,Eproj,Eh,Exc,Evxc',energs%ekin,energs%epot,energs%eproj,energs%eh,energs%exc,energs%evxc
      !end if
      !END DEBUG

      iall=-product(shape(denspot%pot_work))*kind(denspot%pot_work)
      deallocate(denspot%pot_work, stat=istat)
      call memocc(istat, iall, 'denspot%pot_work', subname)

      if(iproc==0) write(*,'(1x,a)') 'Hamiltonian application done.'

      ! Calculate the matrix elements <phi|H|phi>.
      if(.not.tmb%ham_descr%can_use_transposed) then
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

          allocate(tmb%ham_descr%psit_c(tmb%ham_descr%collcom%ndimind_c), stat=istat)
          call memocc(istat, tmb%ham_descr%psit_c, 'tmb%ham_descr%psit_c', subname)
          allocate(tmb%ham_descr%psit_f(7*tmb%ham_descr%collcom%ndimind_f), stat=istat)
          call memocc(istat, tmb%ham_descr%psit_f, 'tmb%ham_descr%psit_f', subname)
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               tmb%ham_descr%psi, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, tmb%ham_descr%lzd)
          tmb%ham_descr%can_use_transposed=.true.
      end if

      allocate(hpsit_c(tmb%ham_descr%collcom%ndimind_c))
      call memocc(istat, hpsit_c, 'hpsit_c', subname)
      allocate(hpsit_f(7*tmb%ham_descr%collcom%ndimind_f))
      call memocc(istat, hpsit_f, 'hpsit_f', subname)
      call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
           tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
      call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
           tmb%ham_descr%psit_c, hpsit_c, tmb%ham_descr%psit_f, hpsit_f, tmb%linmat%ham)
      iall=-product(shape(hpsit_c))*kind(hpsit_c)
      deallocate(hpsit_c, stat=istat)
      call memocc(istat, iall, 'hpsit_c', subname)
      iall=-product(shape(hpsit_f))*kind(hpsit_f)
      deallocate(hpsit_f, stat=istat)
      call memocc(istat, iall, 'hpsit_f', subname)

      !! experimental by SM
      !if (iproc==0) write(*,*) 'deleting additional entries im ham.. SM'
      !allocate(tmb%linmat%ham%matrix(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
      !call memocc(istat, tmb%linmat%ham%matrix, 'tmb%linmat%ham%matrix', subname)
      !call uncompressMatrix(iproc,tmb%linmat%ham)
      !iorb=0
      !do iat=1,at%astruct%nat
      !    if (iproc==0) write(*,*) 'iat, at%astruct%iatype(iat)', iat, at%astruct%iatype(iat)
      !    if (at%astruct%iatype(iat)==1) then
      !        iiorb=4
      !        jjorb=9
      !    else if (at%astruct%iatype(iat)==2) then
      !        iiorb=1
      !        jjorb=1
      !    else
      !        stop 'wrong type'
      !    end if
      !    do i=1,jjorb
      !        iorb=iorb+1
      !        if (i>iiorb) then
      !            tmb%linmat%ham%matrix(:,iorb)=0.d0
      !            tmb%linmat%ham%matrix(iorb,:)=0.d0
      !        end if
      !    end do
      !end do
      !call compress_matrix_for_allreduce(iproc,tmb%linmat%ham)
      !iall=-product(shape(tmb%linmat%ham%matrix))*kind(tmb%linmat%ham%matrix)
      !deallocate(tmb%linmat%ham%matrix, stat=istat)
      !call memocc(istat, iall, 'tmb%linmat%ham%matrix', subname)


      if (scf_mode==LINEAR_FOE) then
         ! NOT ENTIRELY GENERAL HERE - assuming ovrlp is small and ham is large, converting ham to match ovrlp
         call timing(iproc,'FOE_init','ON') !lr408t
         iismall=0
         iseglarge=1
         do isegsmall=1,tmb%linmat%ovrlp%nseg
            do
               is=max(tmb%linmat%ovrlp%keyg(1,isegsmall),tmb%linmat%ham%keyg(1,iseglarge))
               ie=min(tmb%linmat%ovrlp%keyg(2,isegsmall),tmb%linmat%ham%keyg(2,iseglarge))
               iilarge=tmb%linmat%ham%keyv(iseglarge)-tmb%linmat%ham%keyg(1,iseglarge)
               do i=is,ie
                  iismall=iismall+1
                  ham_small%matrix_compr(iismall)=tmb%linmat%ham%matrix_compr(iilarge+i)
               end do
               if (ie>=is) exit
               iseglarge=iseglarge+1
            end do
         end do
         call timing(iproc,'FOE_init','OF') !lr408t
      end if

  else
      if(iproc==0) write(*,*) 'No Hamiltonian application required.'
  end if

  ! CDFT: add V*w_ab to Hamiltonian here - assuming ham and weight matrix have the same sparsity...
  if (present(cdft)) then
     call daxpy(tmb%linmat%ham%nvctr,cdft%lag_mult,cdft%weight_matrix%matrix_compr,1,tmb%linmat%ham%matrix_compr,1)   
  end if

  if (scf_mode/=LINEAR_FOE) then
      allocate(tmb%linmat%ham%matrix(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
      call memocc(istat, tmb%linmat%ham%matrix, 'tmb%linmat%ham%matrix', subname)
      call uncompressMatrix(iproc,tmb%linmat%ham)
      allocate(tmb%linmat%ovrlp%matrix(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
      call memocc(istat, tmb%linmat%ovrlp%matrix, 'tmb%linmat%ovrlp%matrix', subname)
      call uncompressMatrix(iproc,tmb%linmat%ovrlp)
  end if

  ! Diagonalize the Hamiltonian.
  if(scf_mode==LINEAR_MIXPOT_SIMPLE .or. scf_mode==LINEAR_MIXDENS_SIMPLE) then
      ! Keep the Hamiltonian and the overlap since they will be overwritten by the diagonalization.
      allocate(matrixElements(tmb%orbs%norb,tmb%orbs%norb,2), stat=istat)
      call memocc(istat, matrixElements, 'matrixElements', subname)
      call dcopy(tmb%orbs%norb**2, tmb%linmat%ham%matrix(1,1), 1, matrixElements(1,1,1), 1)
      call dcopy(tmb%orbs%norb**2, tmb%linmat%ovrlp%matrix(1,1), 1, matrixElements(1,1,2), 1)
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

      call dcopy(tmb%orbs%norb*tmb%orbs%norb, matrixElements(1,1,1), 1, tmb%coeff(1,1), 1)
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
      ! instead just use -0.5 everywhere
      !tmb%orbs%eval(:) = -0.5_dp

      iall=-product(shape(matrixElements))*kind(matrixElements)
      deallocate(matrixElements, stat=istat)
      call memocc(istat, iall, 'matrixElements', subname)
  else if (scf_mode==LINEAR_DIRECT_MINIMIZATION) then
     if(.not.present(ldiis_coeff)) stop 'ldiis_coeff must be present for scf_mode==LINEAR_DIRECT_MINIMIZATION'
     ! call routine which updates coeffs for tmb%orbs%norb or orbs%norb depending on if extra states are required
     if (extra_states>0) then
        call optimize_coeffs_extra(iproc, nproc, orbs, tmb, ldiis_coeff, fnrm, convcrit_dmin, nitdmin, ebs, &
             curvefit_dmin, extra_states)
     else
        call optimize_coeffs(iproc, nproc, orbs, tmb, ldiis_coeff, fnrm, convcrit_dmin, nitdmin, ebs, curvefit_dmin)
     end if
  end if

  ! CDFT: subtract V*w_ab from Hamiltonian so that we are calculating the correct energy
  if (present(cdft)) then
     call daxpy(tmb%linmat%ham%nvctr,-cdft%lag_mult,cdft%weight_matrix%matrix_compr,1,tmb%linmat%ham%matrix_compr,1)   
  end if

  if (scf_mode/=LINEAR_FOE) then
      call timing(iproc,'getlocbasinit','ON') !lr408t
      ! Calculate the band structure energy and update kernel
      if (scf_mode/=LINEAR_DIRECT_MINIMIZATION) then
         call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%denskern,tmb%linmat%ham,ebs,tmb%coeff,orbs,tmb%orbs,.true.)
      else if (present(cdft)) then
         ! for directmin we have the kernel already, but only the CDFT function not actual energy for CDFT
         call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%denskern,tmb%linmat%ham,ebs,tmb%coeff,orbs,tmb%orbs,.false.)
      end if
      call timing(iproc,'getlocbasinit','OF') !lr408t

      iall=-product(shape(tmb%linmat%ham%matrix))*kind(tmb%linmat%ham%matrix)
      deallocate(tmb%linmat%ham%matrix, stat=istat)
      call memocc(istat, iall, 'tmb%linmat%ham%matrix', subname)
      iall=-product(shape(tmb%linmat%ovrlp%matrix))*kind(tmb%linmat%ovrlp%matrix)
      deallocate(tmb%linmat%ovrlp%matrix, stat=istat)
      call memocc(istat, iall, 'tmb%linmat%ovrlp%matrix', subname)

  else ! foe

      tmprtr=0.d0
      call foe(iproc, nproc, tmb%orbs, tmb%foe_obj, &
           tmprtr, 2, ham_small, tmb%linmat%ovrlp, tmb%linmat%denskern, ebs)
      ! Eigenvalues not available, therefore take -.5d0
      tmb%orbs%eval=-.5d0

  end if


end subroutine get_coeff



subroutine getLocalizedBasis(iproc,nproc,at,orbs,rxyz,denspot,GPU,trH,trH_old,&
    fnrm,infoBasisFunctions,nlpspd,scf_mode, proj,ldiis,SIC,tmb,energs_base,&
    reduce_conf,fix_supportfunctions,nit_precond,target_function,&
    correction_orthoconstraint,nit_basis,deltaenergy_multiplier_TMBexit,deltaenergy_multiplier_TMBfix,&
    ratio_deltas,ortho_on,extra_states,itout)
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
  real(kind=8),dimension(3,at%astruct%nat) :: rxyz
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
  type(energy_terms),intent(in) :: energs_base
  logical,intent(out) :: reduce_conf, fix_supportfunctions
  integer, intent(in) :: nit_precond, target_function, correction_orthoconstraint, nit_basis
  real(kind=8),intent(in) :: deltaenergy_multiplier_TMBexit, deltaenergy_multiplier_TMBfix
  real(kind=8),intent(out) :: ratio_deltas
  logical, intent(inout) :: ortho_on
  integer, intent(in) :: extra_states
  integer,intent(in) :: itout
 
  ! Local variables
  real(kind=8) :: fnrmMax, meanAlpha, ediff, noise, alpha_max, delta_energy, delta_energy_prev
  integer :: iorb, istat, ierr, it, iall, it_tot, ncount, jorb
  real(kind=8),dimension(:),allocatable :: alpha,fnrmOldArr,alphaDIIS, hpsit_c_tmp, hpsit_f_tmp, hpsi_noconf, psidiff
  real(kind=8),dimension(:),allocatable :: hpsi_noprecond, occup_tmp, kernel_compr_tmp
  real(kind=8),dimension(:,:),allocatable :: coeff_old
  logical :: energy_increased, overlap_calculated
  character(len=*),parameter :: subname='getLocalizedBasis'
  real(kind=8),dimension(:),pointer :: lhphiold, lphiold, hpsit_c, hpsit_f, hpsi_small
  type(energy_terms) :: energs
  real(8),dimension(2):: reducearr
  real(gp) :: econf, ediff_sum, delta_energy_prev_sum
  real(kind=8),dimension(3,3) :: interpol_matrix, tmp_matrix
  real(kind=8),dimension(3) :: interpol_vector, interpol_solution
  integer :: i, ist, iiorb, ilr, ii, info
  real(kind=8) :: tt, ddot, d2e, ttt, energy_first
  integer,dimension(3) :: ipiv
  real(kind=8),dimension(:,:),allocatable :: psi_old
  real(kind=8),dimension(:),allocatable :: psi_tmp
  real(kind=8),dimension(3),save :: d2e_arr_out
  integer,save :: isatur_out
  integer :: isatur_in, correction_orthoconstraint_local
  logical :: stop_optimization, energy_increased_previous



  allocate(psi_old(size(tmb%psi),3))
  allocate(psi_tmp(size(tmb%psi)))
  psi_old(:,1)=tmb%psi
  if (itout==1) then
      isatur_out=0
      d2e_arr_out=0
  end if
  stop_optimization=.false.


  ! Allocate all local arrays.
  call allocateLocalArrays()

  !!!EXPERIMENTAL
  !!    call orthonormalizeLocalized(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%npsidim_orbs, tmb%orbs, tmb%lzd, &
  !!         tmb%linmat%ovrlp, tmb%linmat%inv_ovrlp, tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, &
  !!         tmb%can_use_transposed)
  !!    ortho_on=.false.
  !!!END EXPERIMENTAL

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
  !ortho=.true.
  call local_potential_dimensions(tmb%ham_descr%lzd,tmb%orbs,denspot%dpbox%ngatherarr(0,1))
  call start_onesided_communication(iproc, nproc, max(denspot%dpbox%ndimpot,1), denspot%rhov, &
       tmb%ham_descr%comgp%nrecvbuf, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, tmb%ham_descr%lzd)

  reduce_conf=.false.
  fix_supportfunctions=.false.
  delta_energy_prev=1.d100

  ediff_sum=0.d0
  delta_energy_prev_sum=0.d0
  energy_increased_previous=.false.

  isatur_in=0
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
      if (tmb%ham_descr%npsidim_orbs > 0)  call to_zero(tmb%ham_descr%npsidim_orbs,tmb%hpsi(1))
      call small_to_large_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
           tmb%orbs, tmb%psi, tmb%ham_descr%psi)

      call NonLocalHamiltonianApplication(iproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,rxyz,&
           proj,tmb%ham_descr%lzd,nlpspd,tmb%ham_descr%psi,tmb%hpsi,energs%eproj)
      ! only kinetic because waiting for communications
      call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
           tmb%ham_descr%lzd,tmb%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%ham_descr%psi,tmb%hpsi,&
           energs,SIC,GPU,3,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
      call full_local_potential(iproc,nproc,tmb%orbs,tmb%ham_descr%lzd,2,denspot%dpbox,denspot%rhov,denspot%pot_work, &
           tmb%ham_descr%comgp)
      ! only potential
      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          call vcopy(tmb%ham_descr%npsidim_orbs, tmb%hpsi(1), 1, hpsi_noconf(1), 1)
          call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
               tmb%ham_descr%lzd,tmb%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%ham_descr%psi,tmb%hpsi,&
               energs,SIC,GPU,2,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmb%ham_descr%comgp,&
               hpsi_noconf=hpsi_noconf,econf=econf)
          if (nproc>1) then
              call mpiallred(econf, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
          end if
      else
          call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
               tmb%ham_descr%lzd,tmb%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%ham_descr%psi,tmb%hpsi,&
               energs,SIC,GPU,2,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
      end if


      if (target_function==TARGET_FUNCTION_IS_HYBRID .and. iproc==0) then
          write(*,*) 'econf, econf/tmb%orbs%norb',econf, econf/tmb%orbs%norb
      end if

      call timing(iproc,'glsynchham2','ON')
      call SynchronizeHamiltonianApplication(nproc,tmb%ham_descr%npsidim_orbs,tmb%orbs,tmb%ham_descr%lzd,GPU,tmb%hpsi,&
           energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
      call timing(iproc,'glsynchham2','OF')

      ! Apply the orthoconstraint to the gradient. This subroutine also calculates the trace trH.
      if(iproc==0) then
          write(*,'(a)', advance='no') ' Orthoconstraint... '
      end if

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               hpsi_noconf, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
      else
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
      end if

      ncount=sum(tmb%ham_descr%collcom%nrecvcounts_c)
      if(ncount>0) call dcopy(ncount, hpsit_c(1), 1, hpsit_c_tmp(1), 1)
      ncount=7*sum(tmb%ham_descr%collcom%nrecvcounts_f)
      if(ncount>0) call dcopy(ncount, hpsit_f(1), 1, hpsit_f_tmp(1), 1)

      ! optimize the tmbs for a few extra states
      if (target_function==TARGET_FUNCTION_IS_ENERGY.and.extra_states>0) then
          allocate(kernel_compr_tmp(tmb%linmat%denskern%nvctr), stat=istat)
          call memocc(istat, kernel_compr_tmp, 'kernel_compr_tmp', subname)
          call vcopy(tmb%linmat%denskern%nvctr, tmb%linmat%denskern%matrix_compr(1), 1, kernel_compr_tmp(1), 1)
          !allocate(occup_tmp(tmb%orbs%norb), stat=istat)
          !call memocc(istat, occup_tmp, 'occup_tmp', subname)
          !call vcopy(tmb%orbs%norb, tmb%orbs%occup(1), 1, occup_tmp(1), 1)
          !call razero(tmb%orbs%norb,tmb%orbs%occup(1))
          !call vcopy(orbs%norb, orbs%occup(1), 1, tmb%orbs%occup(1), 1)
          !! occupy the next few states - don't need to preserve the charge as only using for support function optimization
          !do iorb=1,tmb%orbs%norb
          !   if (tmb%orbs%occup(iorb)==1.0_gp) then
          !      tmb%orbs%occup(iorb)=2.0_gp
          !   else if (tmb%orbs%occup(iorb)==0.0_gp) then
          !      do jorb=iorb,min(iorb+extra_states-1,tmb%orbs%norb)
          !         tmb%orbs%occup(jorb)=2.0_gp
          !      end do
          !      exit
          !   end if
          !end do
          call calculate_density_kernel(iproc, nproc, .true., tmb%orbs, tmb%orbs, tmb%coeff, tmb%linmat%denskern)
      end if

      correction_orthoconstraint_local=correction_orthoconstraint
      if(.not.ortho_on) then
          correction_orthoconstraint_local=2
      end if

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
         call calculate_energy_and_gradient_linear(iproc, nproc, it, ldiis, fnrmOldArr, alpha, trH, trH_old, fnrm, fnrmMax, &
              meanAlpha, alpha_max, energy_increased, tmb, lhphiold, overlap_calculated, energs_base, &
              hpsit_c, hpsit_f, nit_precond, target_function, correction_orthoconstraint, .false., hpsi_small, hpsi_noprecond)
      else
         call calculate_energy_and_gradient_linear(iproc, nproc, it, ldiis, fnrmOldArr, alpha, trH, trH_old, &
              fnrm, fnrmMax, meanAlpha, alpha_max, energy_increased, tmb, lhphiold, overlap_calculated, &
              energs_base, hpsit_c, hpsit_f, nit_precond, target_function, correction_orthoconstraint, .false., hpsi_small)
      end if

      if (it_tot==1) then
          energy_first=trH
      end if
      if (iproc==0) write(*,'(a,3es16.7)') 'trH, energy_first, (trH-energy_first)/energy_first', &
                                            trH, energy_first, (trH-energy_first)/energy_first
      if ((trH-energy_first)/energy_first>1.d-5) then
          stop_optimization=.true.
          if (iproc==0) write(*,'(a,3es16.7)') 'new stopping crit: trH, energy_first, (trH-energy_first)/energy_first', &
                                                trH, energy_first, (trH-energy_first)/energy_first
      end if


      ! NEW

      if (it_tot==1) then
          if (itout>3) then
              d2e_arr_out(1)=d2e_arr_out(2)
              d2e_arr_out(2)=d2e_arr_out(3)
          end if
          ii=max(min(itout,3),1)
          d2e_arr_out(ii)=trH
      end if


      if (it>3) then
          do i=1,3
              interpol_matrix(1,i)=interpol_matrix(2,i)
              interpol_matrix(2,i)=interpol_matrix(3,i)
              !interpol_matrix(3,i)=interpol_matrix(4,i)
              !interpol_matrix(4,i)=interpol_matrix(5,i)
          end do
          interpol_vector(1)=interpol_vector(2)
          interpol_vector(2)=interpol_vector(3)
          !interpol_vector(3)=interpol_vector(4)
          !interpol_vector(4)=interpol_vector(5)
          psi_old(:,1)=psi_old(:,2)
          psi_old(:,2)=psi_old(:,3)
          !psi_old(:,3)=psi_old(:,4)
          !psi_old(:,4)=psi_old(:,5)
      end if
      psi_tmp=tmb%psi-psi_old(:,1)
      ist=1
      tt=0.d0
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          tt=tt+ddot(ncount, psi_tmp(ist), 1, psi_tmp(ist), 1)
          ist=ist+ncount
      end do
      call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      tt=sqrt(tt/dble(tmb%orbs%norb))
      ii=max(min(it,3),1)
      !interpol_matrix(ii,1)=tt**4
      !interpol_matrix(ii,2)=tt**3
      !interpol_matrix(ii,3)=tt**2
      !interpol_matrix(ii,4)=tt
      !interpol_matrix(ii,5)=tt
      interpol_vector(ii)=trH
      psi_old(:,ii)=tmb%psi

      ! Solve the linear system interpol_matrix*interpol_solution=interpol_vector
      if (it>=3) then
          ttt=0.d0
          do i=1,3
              if (i>1) then
                  !psi_tmp=tmb%psi-psi_old(:,i)
                  !psi_tmp=psi_old(:,i)-psi_old(:,1)
                  psi_tmp=psi_old(:,i)-psi_old(:,i-1)
                  ist=1
                  tt=0.d0
                  do iorb=1,tmb%orbs%norbp
                      iiorb=tmb%orbs%isorb+iorb
                      ilr=tmb%orbs%inwhichlocreg(iiorb)
                      ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
                      tt=tt+ddot(ncount, psi_tmp(ist), 1, psi_tmp(ist), 1)
                      ist=ist+ncount
                  end do
                  call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
                  tt=sqrt(tt/dble(tmb%orbs%norb))
              else
                  tt=0.d0
              end if
              ttt=ttt+tt
              !interpol_matrix(i,1)=ttt**4
              !interpol_matrix(i,2)=ttt**3
              !interpol_matrix(i,3)=ttt**2
              !interpol_matrix(i,4)=ttt
              !interpol_matrix(i,5)=1
              interpol_matrix(i,1)=ttt**2
              interpol_matrix(i,2)=ttt
              interpol_matrix(i,3)=1.d0
          end do
          do i=1,ii
              interpol_solution(i)=interpol_vector(i)
              tmp_matrix(i,1)=interpol_matrix(i,1)
              tmp_matrix(i,2)=interpol_matrix(i,2)
              tmp_matrix(i,3)=interpol_matrix(i,3)
              !tmp_matrix(i,4)=interpol_matrix(i,4)
              !tmp_matrix(i,5)=interpol_matrix(i,5)
          end do
          if (iproc==0) then
              do i=1,3
                  write(*,'(3es14.6,5x,es14.6)') tmp_matrix(i,1:3), interpol_solution(i)
              end do
          end if
          call dgesv(ii, 1, tmp_matrix, 3, ipiv, interpol_solution, 3, info)
          if (info/=0) then
             if (iproc==0) write(*,'(1x,a,i0)') 'ERROR in dgesv (FOE), info=',info
          end if
          if (iproc==0) write(*,'(a,3es14.7)') 'interpol_solution(1:3)',interpol_solution(1:3)
      end if

      !d2e=6.d0*interpol_solution(1)*interpol_matrix(ii,3)+2.d0*interpol_solution(2)
      !d2e = 12.d0*interpol_solution(1)*interpol_matrix(ii,4)**2 + 6.d0*interpol_solution(2)*interpol_matrix(ii,4) + 2.d0*interpol_solution(3)
      d2e = 2.d0*interpol_solution(1)
      tt = abs(interpol_vector(1))*(interpol_vector(1) - 2.d0*interpol_vector(2) + interpol_vector(3))
      ttt = abs(d2e_arr_out(1))*(d2e_arr_out(1) - 2.d0*d2e_arr_out(2) + d2e_arr_out(3))
      !tt=tt/dble(tmb%orbs%norb)
      ttt=ttt/dble(tmb%orbs%norb)
      if (iproc==0) write(*,'(a,2es14.5)') 'tt, ttt', tt, ttt
      if (itout>=3 .and. it >=3) then
          !if (abs(tt)<1.d-4 .and. .not.energy_increased) then
          !if (abs(tt)<1.d-4 .and. abs(ttt)<1.d-3) then
          !if (tt>0.d0 .and. tt<1.d-1 .and. ttt>0.d0 .and. ttt<1.d1 .and. .false.) then
          !!if (tt>0.d0 .and. tt<1.d-1) then
          !!    if (iproc==0) write(*,*) 'SWITCH OFF ORTHO'
          !!    ortho_on=.false.
          !!end if
          !if (abs(tt)<1.d-5 .and. .not.energy_increased) then
          !if (abs(tt)<1.d-1 .and. .not.energy_increased) then
          if (abs(tt)<1.d-2 .and. .not.energy_increased) then
              isatur_in=isatur_in+1
          else
              isatur_in=0
          end if
          if (abs(ttt)<1.d-4) then
              isatur_out=isatur_out+1
          else
              isatur_out=0
          end if
          if (iproc==0) then
              write(*,'(a,3es12.4,2i4)') 'd2e',d2e,tt, ttt, isatur_in, isatur_out
          end if

          !if (isatur_in>=2) then
          if (isatur_in>=200) then
              stop_optimization=.true.
          end if

          !!if (isatur_in>=2 .and. isatur_out>=2) then
          !!    if (iproc==0) write(*,*) 'new fixing criterion'
          !!    fix_supportfunctions=.true.
          !!end if
      end if

      if (target_function==TARGET_FUNCTION_IS_ENERGY.and.extra_states>0) then
          !call vcopy(tmb%orbs%norb, occup_tmp(1), 1, tmb%orbs%occup(1), 1)
          !iall=-product(shape(occup_tmp))*kind(occup_tmp)
          !deallocate(occup_tmp, stat=istat)
          !call memocc(istat, iall, 'occup_tmp', subname)
          call vcopy(tmb%linmat%denskern%nvctr, kernel_compr_tmp(1), 1, tmb%linmat%denskern%matrix_compr(1), 1)
          iall=-product(shape(kernel_compr_tmp))*kind(kernel_compr_tmp)
          deallocate(kernel_compr_tmp, stat=istat)
          call memocc(istat, iall, 'kernel_compr_tmp', subname)
      end if

      ediff=trH-trH_old

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          if (.not.energy_increased .and. .not.energy_increased_previous) then
              ratio_deltas=ediff/delta_energy_prev
          end if
          if (.not.energy_increased_previous .and. it>1) then
              ediff_sum=ediff_sum+ediff
              delta_energy_prev_sum=delta_energy_prev_sum+delta_energy_prev
          end if
          if (ldiis%switchSD) then
              !!ratio_deltas=0.5d0
              !!if (iproc==0) write(*,*) 'WARNING: TEMPORARY FIX for ratio_deltas!'
          end if
          !if (iproc==0) write(*,*) 'WARNING: HACK FOR ratio_deltas, set to 1.d0!!'
          !if (iproc==0) write(*,*) 'WARNING: HACK FOR ratio_deltas, set to 0.5d0*(ratio_deltas+1.d0)!!'
          !ratio_deltas=0.5d0*(ratio_deltas+1.d0)
          !ratio_deltas=1.d0
          if (iproc==0) write(*,*) 'ediff, delta_energy_prev', ediff, delta_energy_prev
          if (iproc==0) write(*,*) 'ratio_deltas',ratio_deltas
          if ((ediff>deltaenergy_multiplier_TMBexit*delta_energy_prev .or. energy_increased) .and. it>1) then
          !if ((it>=nit_basis .or.  energy_increased) .and. it>1) then
              if (iproc==0) write(*,*) 'reduce the confinement'
              reduce_conf=.true.
          end if
      end if

      if (energy_increased) then
          energy_increased_previous=.true.
      else
          energy_increased_previous=.false.
      end if


      !!if ((ediff>deltaenergy_multiplier_TMBfix*delta_energy_prev .and. .not.energy_increased) .and. it>1 .and. &
      !!    target_function==TARGET_FUNCTION_IS_HYBRID) then
      !!    if (iproc==0) write(*,*) 'Will fix the support functions'
      !!    fix_supportfunctions=.true.
      !!end if

      !!delta_energy_prev=delta_energy

      if (energy_increased) then
          if (iproc==0) write(*,*) 'WARNING: ENERGY INCREASED'
          tmb%ham_descr%can_use_transposed=.false.
          call dcopy(tmb%npsidim_orbs, lphiold(1), 1, tmb%psi(1), 1)
          if (scf_mode/=LINEAR_FOE) then
              ! Recalculate the kernel with the old coefficients
              call dcopy(tmb%orbs%norb*tmb%orbs%norb, coeff_old(1,1), 1, tmb%coeff(1,1), 1)
              call calculate_density_kernel(iproc, nproc, .true., orbs, tmb%orbs, &
                   tmb%coeff, tmb%linmat%denskern)
          end if
          trH_old=0.d0
          it=it-2 !go back one iteration (minus 2 since the counter was increased)
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
          if(iproc==0) write(*,*) 'it_tot',it_tot
          overlap_calculated=.false.
          ! print info here anyway for debugging
          if (iproc==0) write(*,'(1x,a,i6,2es15.7,f17.10,2es13.4)') 'iter, fnrm, fnrmMax, ebs, diff, noise level', &
          it, fnrm, fnrmMax, trH, ediff,noise
          if (it_tot<2*nit_basis) then ! just in case the step size is the problem
             cycle
          else if(it_tot<3*nit_basis) then ! stop orthonormalizing the tmbs
             if (iproc==0) write(*,*) 'WARNING: SWITCHING OFF ORTHO COMMENTED'
             !if (iproc==0) write(*,'(a)') 'Energy increasing, switching off orthonormalization of tmbs'
             !ortho_on=.false.
             !alpha=alpha*5.0d0/3.0d0 ! increase alpha to make up for decrease from previous iteration
          end if
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
      !!if(it>=nit_basis .or. it_tot>=3*nit_basis .or. reduce_conf) then
      !!    if(it>=nit_basis .and. .not.energy_increased) then
      !if(it>=nit_basis .or. it_tot>=3*nit_basis) then
      if(it>=nit_basis .or. it_tot>=3*nit_basis .or. stop_optimization) then
          if(it>=nit_basis) then
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
          else if (stop_optimization) then
              if (iproc==0) then
                  write(*,'(1x,a,2es15.7,f15.7)') 'Final values for fnrm, fnrmMax, hybrid: ', fnrm, fnrmMax, trH
              end if
              infoBasisFunctions=0
          !!else if (reduce_conf) then
          !!    if (iproc==0) then
          !!        write(*,'(1x,a,2es15.7,f15.7)') 'Final values for fnrm, fnrmMax, hybrid: ', fnrm, fnrmMax, trH
          !!    end if
          !!    infoBasisFunctions=0
          end if
          if(iproc==0) write(*,'(1x,a)') '============================= Basis functions created. ============================='
          if (infoBasisFunctions>=0) then
              ! Calculate the Hamiltonian matrix, since we have all quantities ready. This matrix can then be used in the first
              ! iteration of get_coeff.
              call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
                   tmb%ham_descr%psit_c, hpsit_c_tmp, tmb%ham_descr%psit_f, hpsit_f_tmp, tmb%linmat%ham)
          end if

          exit iterLoop
      end if
      trH_old=trH

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
         call hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, &
              lphiold, alpha, trH, meanAlpha, alpha_max, alphaDIIS, hpsi_small, ortho_on, psidiff)
      else
         call hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, &
              lphiold, alpha, trH, meanAlpha, alpha_max, alphaDIIS, hpsi_small, ortho_on)
      end if


      overlap_calculated=.false.
      ! It is now not possible to use the transposed quantities, since they have changed.
      if(tmb%ham_descr%can_use_transposed) then
          iall=-product(shape(tmb%ham_descr%psit_c))*kind(tmb%ham_descr%psit_c)
          deallocate(tmb%ham_descr%psit_c, stat=istat)
          call memocc(istat, iall, 'tmb%ham_descr%psit_c', subname)
          iall=-product(shape(tmb%ham_descr%psit_f))*kind(tmb%ham_descr%psit_f)
          deallocate(tmb%ham_descr%psit_f, stat=istat)
          call memocc(istat, iall, 'tmb%ham_descr%psit_f', subname)
          tmb%ham_descr%can_use_transposed=.false.
      end if

      ! Estimate the energy change, that is to be expected in the next optimization
      ! step, given by the product of the force and the "displacement" .
      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          call estimate_energy_change(tmb%npsidim_orbs, tmb%orbs, tmb%lzd, psidiff, hpsi_noprecond, delta_energy)
          ! This is a hack...
          if (energy_increased) then
              delta_energy=1.d100
              !ratio_deltas=1.d100
          end if
          if (iproc==0) write(*,*) 'delta_energy', delta_energy
          delta_energy_prev=delta_energy
      end if

      ! Copy the coefficients to coeff_old. The coefficients will be modified in reconstruct_kernel.
      if (scf_mode/=LINEAR_FOE) then
          call dcopy(tmb%orbs%norb*tmb%orbs%norb, tmb%coeff(1,1), 1, coeff_old(1,1), 1)
      end if

      if(scf_mode/=LINEAR_FOE) then
          call reconstruct_kernel(iproc, nproc, 1, tmb%orthpar%blocksize_pdsyev, tmb%orthpar%blocksize_pdgemm, &
               orbs, tmb, overlap_calculated)
      else
          call purify_kernel(iproc, nproc, tmb, overlap_calculated)
      end if
      if(iproc==0) then
          write(*,'(a)') 'done.'
      end if

  end do iterLoop

  if (iproc==0) write(*,'(a,3es16.6)') 'NEW RATIO: ediff_sum, delta_energy_prev_sum, ratio_deltas', &
      ediff_sum, delta_energy_prev_sum, ratio_deltas
  ratio_deltas=ediff_sum/delta_energy_prev_sum

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

      allocate(hpsi_small(max(tmb%npsidim_orbs,tmb%npsidim_comp)), stat=istat)
      call memocc(istat, hpsi_small, 'hpsi_small', subname)
    
      allocate(lhphiold(max(tmb%npsidim_orbs,tmb%npsidim_comp)), stat=istat)
      call memocc(istat, lhphiold, 'lhphiold', subname)

      allocate(lphiold(size(tmb%psi)), stat=istat)
      call memocc(istat, lphiold, 'lphiold', subname)

      allocate(hpsit_c(sum(tmb%ham_descr%collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, hpsit_c, 'hpsit_c', subname)

      allocate(hpsit_f(7*sum(tmb%ham_descr%collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, hpsit_f, 'hpsit_f', subname)

      allocate(hpsit_c_tmp(sum(tmb%ham_descr%collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, hpsit_c_tmp, 'hpsit_c_tmp', subname)

      allocate(hpsit_f_tmp(7*sum(tmb%ham_descr%collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, hpsit_f_tmp, 'hpsit_f_tmp', subname)

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
         allocate(hpsi_noconf(tmb%ham_descr%npsidim_orbs), stat=istat)
         call memocc(istat, hpsi_noconf, 'hpsi_noconf', subname)

         allocate(psidiff(tmb%npsidim_orbs), stat=istat)
         call memocc(istat, psidiff, 'psidiff', subname)

         allocate(hpsi_noprecond(tmb%npsidim_orbs), stat=istat)
         call memocc(istat, hpsi_noprecond, 'hpsi_noprecond', subname)
      end if

      if (scf_mode/=LINEAR_FOE) then
          allocate(coeff_old(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
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

      iall=-product(shape(hpsi_small))*kind(hpsi_small)
      deallocate(hpsi_small, stat=istat)
      call memocc(istat, iall, 'hpsi_small', subname)

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



subroutine improveOrbitals(iproc, tmb, ldiis, alpha, gradient)
  use module_base
  use module_types
  use module_interfaces, except_this_one => improveOrbitals
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc
  type(DFT_wavefunction),intent(inout) :: tmb
  type(localizedDIISParameters),intent(inout) :: ldiis
  real(kind=8),dimension(tmb%orbs%norbp),intent(in) :: alpha
  real(kind=wp),dimension(max(tmb%npsidim_orbs,tmb%npsidim_comp)),intent(inout) :: gradient
  
  ! Local variables
  integer :: istart, iorb, iiorb, ilr, ncount

  if(ldiis%isx==0) then ! steepest descents
      call timing(iproc,'optimize_SD   ','ON')
      istart=1
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          call daxpy(ncount, -alpha(iorb), gradient(istart), 1, tmb%psi(istart), 1)
          istart=istart+ncount
      end do
      call timing(iproc,'optimize_SD   ','OF')
  else! DIIS
      ldiis%mis=mod(ldiis%is,ldiis%isx)+1
      ldiis%is=ldiis%is+1
      if(ldiis%alphaDIIS/=1.d0) then
          call dscal(max(tmb%npsidim_orbs,tmb%npsidim_comp), ldiis%alphaDIIS, gradient, 1)
      end if
      call optimizeDIIS(iproc, max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%lzd, gradient, tmb%psi, ldiis)
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
  if(info/=0)then
    write(*,*) 'ERROR: dsygv in diagonalizeHamiltonian2, info=',info,'N=',orbs%norb
  end if

  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  call memocc(istat, iall, 'work', subname)

  call timing(iproc,'diagonal_seq  ','OF')

end subroutine diagonalizeHamiltonian2

subroutine small_to_large_locreg(iproc, npsidim_orbs_small, npsidim_orbs_large, lzdsmall, lzdlarge, &
       orbs, phismall, philarge)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, npsidim_orbs_small, npsidim_orbs_large
  type(local_zone_descriptors),intent(in) :: lzdsmall, lzdlarge
  type(orbitals_data),intent(in) :: orbs
  real(kind=8),dimension(npsidim_orbs_small),intent(in) :: phismall
  real(kind=8),dimension(npsidim_orbs_large),intent(out) :: philarge
  
  ! Local variables
  integer :: ists, istl, iorb, ilr, sdim, ldim, nspin
       call timing(iproc,'small2large','ON') ! lr408t 
  ! No need to put arrays to zero, Lpsi_to_global2 will handle this.
  call to_zero(npsidim_orbs_large, philarge(1))
  ists=1
  istl=1
  do iorb=1,orbs%norbp
      ilr = orbs%inWhichLocreg(orbs%isorb+iorb)
      sdim=lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
      ldim=lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f
      nspin=1 !this must be modified later
      call Lpsi_to_global2(iproc, sdim, ldim, orbs%norb, orbs%nspinor, nspin, lzdlarge%llr(ilr), &
           lzdsmall%llr(ilr), phismall(ists), philarge(istl))
      ists=ists+lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
      istl=istl+lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f
  end do
  if(orbs%norbp>0 .and. ists/=npsidim_orbs_small+1) then
      write(*,'(3(a,i0))') 'ERROR on process ',iproc,': ',ists,'=ists /= npsidim_orbs_small+1=',npsidim_orbs_small+1
      stop
  end if
  if(orbs%norbp>0 .and. istl/=npsidim_orbs_large+1) then
      write(*,'(3(a,i0))') 'ERROR on process ',iproc,': ',istl,'=istk /= npsidim_orbs_large+1=',npsidim_orbs_large+1
      stop
  end if
       call timing(iproc,'small2large','OF') ! lr408t 
end subroutine small_to_large_locreg


subroutine large_to_small_locreg(iproc, npsidim_orbs_small, npsidim_orbs_large, lzdsmall, lzdlarge, &
       orbs, philarge, phismall)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, npsidim_orbs_small, npsidim_orbs_large
  type(local_zone_descriptors),intent(in) :: lzdsmall, lzdlarge
  type(orbitals_data),intent(in) :: orbs
  real(kind=8),dimension(npsidim_orbs_large),intent(in) :: philarge
  real(kind=8),dimension(npsidim_orbs_small),intent(out) :: phismall
  
  ! Local variables
  integer :: istl, ists, ilr, ldim, gdim, iorb
       call timing(iproc,'large2small','ON') ! lr408t   
  ! Transform back to small locreg
  ! No need to this array to zero, since all values will be filled with a value during the copy.
  !!call to_zero(npsidim_orbs_small, phismall(1))
  ists=1
  istl=1
  do iorb=1,orbs%norbp
      ilr = orbs%inWhichLocreg(orbs%isorb+iorb)
      ldim=lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
      gdim=lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f
      call psi_to_locreg2(iproc, ldim, gdim, lzdsmall%llr(ilr), lzdlarge%llr(ilr), &
           philarge(istl:istl+gdim-1), phismall(ists:ists+ldim-1))
      ists=ists+lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
      istl=istl+lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f
  end do

  if(orbs%norbp>0 .and. ists/=npsidim_orbs_small+1) stop 'ists/=npsidim_orbs_small+1'
  if(orbs%norbp>0 .and. istl/=npsidim_orbs_large+1) stop 'istl/=npsidim_orbs_large+1'
       call timing(iproc,'large2small','OF') ! lr408t 
end subroutine large_to_small_locreg







subroutine communicate_basis_for_density_collective(iproc, nproc, lzd, npsidim, orbs, lphi, collcom_sr)
  use module_base
  use module_types
  use module_interfaces, except_this_one => communicate_basis_for_density_collective
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, npsidim
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  real(kind=8),dimension(npsidim),intent(in) :: lphi
  type(collective_comms),intent(inout) :: collcom_sr
  
  ! Local variables
  integer :: ist, istr, iorb, iiorb, ilr, istat, iall
  real(kind=8),dimension(:),allocatable :: psir, psirwork, psirtwork
  type(workarr_sumrho) :: w
  character(len=*),parameter :: subname='comm_basis_for_dens_coll'

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

  call transpose_switch_psir(collcom_sr, psir, psirwork)

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
  real(kind=8),dimension(max(tmbopt%npsidim_orbs,tmbopt%npsidim_comp)),intent(out):: lphioldopt
  
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
  if(trH<=ldiis%trmin+1.d-12*abs(ldiis%trmin) .and. .not.ldiis%resetDIIS) then !1.d-12 is here to tolerate some noise...
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

subroutine reconstruct_kernel(iproc, nproc, inversion_method, blocksize_dsyev, blocksize_pdgemm, &
           orbs, tmb, overlap_calculated)
  use module_base
  use module_types
  use module_interfaces, except_this_one => reconstruct_kernel
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, inversion_method, blocksize_dsyev, blocksize_pdgemm
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(inout):: tmb
  logical,intent(inout):: overlap_calculated

  ! Local variables
  integer:: istat, iall
  character(len=*),parameter:: subname='reconstruct_kernel'

  !call timing(iproc,'renormCoefComp','ON')

  ! Calculate the overlap matrix between the TMBs.
  if(.not. overlap_calculated) then
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
         call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
              tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
         tmb%can_use_transposed=.true.
     end if
     !call timing(iproc,'renormCoefComp','OF')

     call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, &
          tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%ovrlp)

     !call timing(iproc,'renormCoefComp','ON')
     overlap_calculated=.true.
  end if

  allocate(tmb%linmat%ovrlp%matrix(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, tmb%linmat%ovrlp%matrix, 'tmb%linmat%ovrlp%matrix', subname)
  call uncompressMatrix(iproc,tmb%linmat%ovrlp)
  call reorthonormalize_coeff(iproc, nproc, orbs%norb, blocksize_dsyev, blocksize_pdgemm, inversion_method, &
       tmb%orbs, tmb%linmat%ovrlp, tmb%coeff, orbs)

  iall=-product(shape(tmb%linmat%ovrlp%matrix))*kind(tmb%linmat%ovrlp%matrix)
  deallocate(tmb%linmat%ovrlp%matrix, stat=istat)
  call memocc(istat, iall, 'tmb%linmat%ovrlp%matrix', subname)

  ! Recalculate the kernel
  call calculate_density_kernel(iproc, nproc, .true., orbs, tmb%orbs, tmb%coeff, tmb%linmat%denskern)

end subroutine reconstruct_kernel

!> Passing sparse ovrlp, but for now assuming ovrlp%matrix will be allocated and filled if using dense
subroutine reorthonormalize_coeff(iproc, nproc, norb, blocksize_dsyev, blocksize_pdgemm, inversion_method, basis_orbs, &
           basis_overlap, coeff, orbs)
  use module_base
  use module_types
  use module_interfaces, except_this_one => reorthonormalize_coeff
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc, nproc, norb
  integer, intent(in) :: blocksize_dsyev, blocksize_pdgemm, inversion_method
  type(orbitals_data), intent(in) :: basis_orbs   !number of basis functions
  type(orbitals_data), optional, intent(in) :: orbs   !Kohn-Sham orbitals that will be orthonormalized and their parallel distribution
  type(sparseMatrix),intent(in) :: basis_overlap
  real(kind=8),dimension(basis_orbs%norb,basis_orbs%norb),intent(inout) :: coeff
  ! Local variables
  integer :: ierr, istat, iall, ind, iorb, korb, llorb, jorb
  integer :: npts_per_proc, ind_start, ind_end, indc
  real(kind=8), dimension(:,:), allocatable :: coeff_tmp, ovrlp_coeff, ovrlp_coeff2, coefftrans
  character(len=*),parameter:: subname='reorthonormalize_coeff'
  !integer :: iorb, jorb !DEBUG
  real(kind=8) :: tt!, tt2, tt3, ddot   !DEBUG
  logical :: dense
  integer,parameter :: ALLGATHERV=1, ALLREDUCE=2
  integer, parameter :: communication_strategy=ALLGATHERV

  call mpi_barrier(bigdft_mpi%mpi_comm, ierr) ! to check timings
  call timing(iproc,'renormCoefCom1','ON')

  !if (present(orbs)) then
  !   communication_strategy=ALLREDUCE
  !else
  !   communication_strategy=ALLGATHERV
  !end if

  allocate(ovrlp_coeff(norb,norb), stat=istat)
  call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)

  allocate(coeff_tmp(basis_orbs%norbp,max(norb,1)), stat=istat)
  call memocc(istat, coeff_tmp, 'coeff_tmp', subname)

  if(iproc==0) then
      write(*,'(a)',advance='no') 'coeff renormalization...'
  end if

  dense=.true.

  if (dense) then
     ! Calculate the overlap matrix among the coefficients with respect to basis_overlap.
     if (basis_orbs%norbp>0) then
         call dgemm('n', 'n', basis_orbs%norbp, norb, basis_orbs%norb, 1.d0, basis_overlap%matrix(basis_orbs%isorb+1,1), &
              basis_orbs%norb, coeff(1,1), basis_orbs%norb, 0.d0, coeff_tmp, basis_orbs%norbp)
         call dgemm('t', 'n', norb, norb, basis_orbs%norbp, 1.d0, coeff(basis_orbs%isorb+1,1), &
              basis_orbs%norb, coeff_tmp, basis_orbs%norbp, 0.d0, ovrlp_coeff, norb)
      else
         call to_zero(norb**2, ovrlp_coeff(1,1))
      end if
  else ! sparse - still less efficient than dense, also needs moving to a subroutine

     call to_zero(norb**2, ovrlp_coeff(1,1))
     npts_per_proc = nint(real(basis_overlap%nvctr + basis_overlap%full_dim1,dp) / real(nproc*2,dp))
     ind_start = 1+iproc*npts_per_proc
     ind_end = (iproc+1)*npts_per_proc
     if (iproc==nproc-1) ind_end = basis_overlap%nvctr!ceiling(0.5d0*real(basis_overlap%nvctr + basis_overlap%full_dim1,dp))

     indc=0
     do ind = 1, basis_overlap%nvctr
        korb = basis_overlap%orb_from_index(1,ind)
        llorb = basis_overlap%orb_from_index(2,ind)
        if (korb<llorb) cycle ! so still only doing half
        indc = indc + 1
        if (indc < ind_start .or. indc > ind_end) cycle

        do iorb=1,norb
             if (llorb==korb) then
                tt=basis_overlap%matrix_compr(ind)*coeff(korb,iorb)
                do jorb=iorb,norb
                    ovrlp_coeff(jorb,iorb)=ovrlp_coeff(jorb,iorb) &
                         +coeff(llorb,jorb)*tt
                end do
             else
                do jorb=iorb,norb
                    ovrlp_coeff(jorb,iorb)=ovrlp_coeff(jorb,iorb) &
                         +(coeff(llorb,iorb)*coeff(korb,jorb)+coeff(llorb,jorb)*coeff(korb,iorb))&
                         *basis_overlap%matrix_compr(ind)
                end do
             end if
         end do
     end do

     ! use symmetry to calculate other half
     do iorb=1,norb
        do jorb=iorb+1,norb
           ovrlp_coeff(iorb,jorb) = ovrlp_coeff(jorb,iorb)
        end do
     end do

  end if !sparse/dense

  if (nproc>1) then
      call timing(iproc,'renormCoefCom1','OF')
      call timing(iproc,'renormCoefComm','ON')
      call mpiallred(ovrlp_coeff(1,1), norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      call timing(iproc,'renormCoefComm','OF')
      call timing(iproc,'renormCoefCom1','ON')
  end if

  ! Recalculate the kernel.
  allocate(ovrlp_coeff2(norb,norb), stat=istat)
  call memocc(istat, ovrlp_coeff2, 'ovrlp_coeff2', subname)

  call timing(iproc,'renormCoefCom1','OF')

  call overlapPowerPlusMinusOneHalf_old(iproc, nproc, bigdft_mpi%mpi_comm, inversion_method, &
       blocksize_dsyev, blocksize_pdgemm, norb, ovrlp_coeff, ovrlp_coeff2, .false., orbs)
  call timing(iproc,'renormCoefCom2','ON')

  iall=-product(shape(ovrlp_coeff))*kind(ovrlp_coeff)
  deallocate(ovrlp_coeff,stat=istat)
  call memocc(istat,iall,'ovrlp_coeff',subname)

  ! Build the new linear combinations
  !call dgemm('n', 'n', basis_orbs%norb, orbs%norb, orbs%norb, 1.d0, coeff(1,1), basis_orbs%norb, &
  !     ovrlp_coeff2(1,1), orbs%norb, 0.d0, coeff_tmp(1,1), basis_orbs%norb)
  !call dcopy(basis_orbs%norb*orbs%norb,coeff_tmp(1,1),1,coeff(1,1),1)

  ! Build the new linear combinations - all gather would be better, but allreduce easier for now

  iall=-product(shape(coeff_tmp))*kind(coeff_tmp)
  deallocate(coeff_tmp,stat=istat)
  call memocc(istat,iall,'coeff_tmp',subname)

  if (communication_strategy==ALLREDUCE) then
     allocate(coeff_tmp(basis_orbs%norb,orbs%norb), stat=istat)
     call memocc(istat, coeff_tmp, 'coeff_tmp', subname)

     if (orbs%norbp>0) then
        call dgemm('n', 't', basis_orbs%norb, orbs%norb, orbs%norbp, 1.d0, coeff(1,orbs%isorb+1), basis_orbs%norb, &
             ovrlp_coeff2(1,orbs%isorb+1), orbs%norb, 0.d0, coeff_tmp(1,1), basis_orbs%norb)
     else
        call to_zero(basis_orbs%norb*orbs%norb, coeff_tmp(1,1))
     end if

     if (nproc>1) then
         call mpiallred(coeff_tmp(1,1), basis_orbs%norb*orbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)
     end if
     call dcopy(basis_orbs%norb*orbs%norb,coeff_tmp(1,1),1,coeff(1,1),1)
  else
     allocate(coeff_tmp(norb,basis_orbs%norbp), stat=istat)
     call memocc(istat, coeff_tmp, 'coeff_tmp', subname)
     ! need to transpose so we can allgather - NOT VERY ELEGANT
     if (basis_orbs%norbp>0) then
         call dgemm('n', 't', norb, basis_orbs%norbp, norb, 1.d0, ovrlp_coeff2(1,1), norb, &
             coeff(1+basis_orbs%isorb,1), basis_orbs%norb, 0.d0, coeff_tmp(1,1), norb)
     end if

     allocate(coefftrans(norb,basis_orbs%norb), stat=istat)
     call memocc(istat, coefftrans, 'coefftrans', subname)

     ! gather together
     if(nproc > 1) then
        call mpi_allgatherv(coeff_tmp(1,1), basis_orbs%norbp*norb, mpi_double_precision, coefftrans(1,1), &
           norb*basis_orbs%norb_par(:,0), norb*basis_orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
     else
        call dcopy(basis_orbs%norbp*norb,coeff_tmp(1,1),1,coefftrans(1,1),1)
     end if

     ! untranspose coeff
     do iorb=1,norb
        do jorb=1,basis_orbs%norb
           coeff(jorb,iorb) = coefftrans(iorb,jorb)
        end do
     end do
     iall=-product(shape(coefftrans))*kind(coefftrans)
     deallocate(coefftrans,stat=istat)
     call memocc(istat,iall,'coefftrans',subname)
  end if

  call timing(iproc,'renormCoefCom2','OF')

  !!!DEBUG
  !!! Check normalization
  !!iall=-product(shape(coeff_tmp))*kind(coeff_tmp)
  !!deallocate(coeff_tmp,stat=istat)
  !!call memocc(istat,iall,'coeff_tmp',subname)
  !!allocate(coeff_tmp(basis_orbs%norb,basis_orbs%norb),stat=istat)
  !!call memocc(istat,coeff_tmp,'coeff_tmp',subname)
  !!call dgemm('n', 'n', basis_orbs%norb, basis_orbs%norb, basis_orbs%norb, 1.d0, basis_overlap(1,1), basis_orbs%norb, &
  !!     coeff(1,1), basis_orbs%norb, 0.d0, coeff_tmp(1,1), basis_orbs%norb)
  !!do iorb=1,basis_orbs%norb
  !!    do jorb=1,basis_orbs%norb
  !!        tt=ddot(basis_orbs%norb, coeff(1,iorb), 1, coeff_tmp(1,jorb), 1)
  !!        tt2=ddot(basis_orbs%norb, coeff_tmp(1,iorb), 1, coeff(1,jorb), 1)
  !!        tt3=ddot(basis_orbs%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
  !!        if(iproc==0) write(200,'(2i6,3es15.5)') iorb, jorb, tt, tt2, tt3
  !!    end do
  !!end do
  !!! END DEBUG



  iall=-product(shape(ovrlp_coeff2))*kind(ovrlp_coeff2)
  deallocate(ovrlp_coeff2,stat=istat)
  call memocc(istat,iall,'ovrlp_coeff2',subname)

  iall=-product(shape(coeff_tmp))*kind(coeff_tmp)
  deallocate(coeff_tmp,stat=istat)
  call memocc(istat,iall,'coeff_tmp',subname)

end subroutine reorthonormalize_coeff


!> Estimate the energy change, given by the product of the force and the "displacement" .
subroutine estimate_energy_change(npsidim_orbs, orbs, lzd, psidiff, hpsi_noprecond, delta_energy)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer, intent(in) :: npsidim_orbs
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  real(kind=8),dimension(npsidim_orbs),intent(in) :: psidiff, hpsi_noprecond
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



subroutine purify_kernel(iproc, nproc, tmb, overlap_calculated)
  use module_base
  use module_types
  use module_interfaces, except_this_one => purify_kernel
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(DFT_wavefunction),intent(inout):: tmb
  logical,intent(inout):: overlap_calculated

  ! Local variables
  integer :: istat, iall, it, lwork, info, iorb, jorb
  real(kind=8),dimension(:,:),allocatable :: ks, ksk, ksksk, kernel, overlap
  real(kind=8),dimension(:),allocatable :: eval, work
  character(len=*),parameter :: subname='purify_kernel'
  real(kind=8) :: dnrm2

  ! Calculate the overlap matrix between the TMBs.
  if(.not. overlap_calculated) then
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
         call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
              tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
         tmb%can_use_transposed=.true.
     end if
     !call timing(iproc,'renormCoefComp','OF')

     call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, &
          tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%ovrlp)

     !call timing(iproc,'renormCoefComp','ON')
     overlap_calculated=.true.
  end if

  allocate(tmb%linmat%ovrlp%matrix(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, tmb%linmat%ovrlp%matrix, 'tmb%linmat%ovrlp%matrix', subname)
  allocate(tmb%linmat%denskern%matrix(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, tmb%linmat%denskern%matrix, 'tmb%linmat%denskern%matrix', subname)
  call uncompressMatrix(iproc,tmb%linmat%ovrlp)
  call uncompressMatrix(iproc,tmb%linmat%denskern)

  allocate(ks(tmb%orbs%norb,tmb%orbs%norb),stat=istat)
  call memocc(istat, ks, 'ks', subname) 
  allocate(ksk(tmb%orbs%norb,tmb%orbs%norb),stat=istat)
  call memocc(istat, ksk, 'ksk', subname) 
  allocate(ksksk(tmb%orbs%norb,tmb%orbs%norb),stat=istat)
  call memocc(istat, ksksk, 'ksksk', subname) 

  allocate(kernel(tmb%orbs%norb,tmb%orbs%norb),stat=istat)
  !call memocc(istat, kernel, 'kernel', subname) 
  allocate(overlap(tmb%orbs%norb,tmb%orbs%norb),stat=istat)
  !call memocc(istat, overlap, 'overlap', subname) 
  allocate(eval(tmb%orbs%norb),stat=istat)
  !call memocc(istat, eval, 'eval', subname) 
  lwork=10*tmb%orbs%norb
  allocate(work(lwork),stat=istat)
  !call memocc(istat, work, 'work', subname) 

  tmb%linmat%denskern%matrix=0.5d0*tmb%linmat%denskern%matrix

  do it=1,3

      !!if (iproc==0) then
      !!    do iorb=1,tmb%orbs%norb
      !!        do jorb=1,tmb%orbs%norb
      !!            if (abs(tmb%linmat%ovrlp%matrix(iorb,jorb))-abs(tmb%linmat%ovrlp%matrix(jorb,iorb))>1.d-20) then
      !!                write(*,'(a,2es18.8)') 'NOT SYMM', tmb%linmat%ovrlp%matrix(iorb,jorb), tmb%linmat%ovrlp%matrix(jorb,iorb)
      !!            end if
      !!        end do
      !!    end do
      !!end if

      !!tmb%linmat%ovrlp%matrix=0.d0
      !!do istat=1,tmb%orbs%norb
      !!    tmb%linmat%ovrlp%matrix(istat,istat)=1.1d0
      !!end do

      kernel=tmb%linmat%denskern%matrix
      overlap=tmb%linmat%ovrlp%matrix

      !!call dsygv(1, 'n', 'l', tmb%orbs%norb, kernel, tmb%orbs%norb, overlap, tmb%orbs%norb, eval, work, lwork, info)
      !!if (info==0) then
      !!    if (iproc==0) then
      !!        write(*,*) 'eval min: ',eval(1)
      !!        write(*,*) 'eval max: ',eval(tmb%orbs%norb)
      !!    end if
      !!else
      !!    if (iproc==0) write(*,*) 'ERROR dsygv: info=',info
      !!    call mpi_finalize(info)
      !!    stop
      !!end if

      call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, tmb%linmat%denskern%matrix(1,1), tmb%orbs%norb, &
                 tmb%linmat%ovrlp%matrix(1,1), tmb%orbs%norb, 0.d0, ks(1,1), tmb%orbs%norb) 
      call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, ks(1,1), tmb%orbs%norb, &
                 tmb%linmat%denskern%matrix(1,1), tmb%orbs%norb, 0.d0, ksk(1,1), tmb%orbs%norb)
      call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, ks(1,1), tmb%orbs%norb, &
                 ksk(1,1), tmb%orbs%norb, 0.d0, ksksk(1,1), tmb%orbs%norb)
      if (iproc==0) write(*,*) 'PURIFYING THE KERNEL'
      !!if (iproc==0) then
      !!    do istat=1,tmb%orbs%norb
      !!        do iall=1,tmb%orbs%norb
      !!            write(200+iproc,*) istat, iall, tmb%linmat%denskern%matrix(iall,istat)
      !!            write(300+iproc,*) istat, iall, ks(iall,istat)
      !!            write(400+iproc,*) istat, iall, ksk(iall,istat)
      !!            write(500+iproc,*) istat, iall, ksksk(iall,istat)
      !!        end do
      !!    end do
      !!end if
      overlap=ksk-tmb%linmat%denskern%matrix
      if (iproc==0) write(*,*) 'diff from idempotency', dnrm2(tmb%orbs%norb**2, overlap, 1)
      tmb%linmat%denskern%matrix = 3*ksk-2*ksksk
      !tmb%linmat%denskern%matrix = tmb%linmat%denskern%matrix/1.1d0
      !!if (iproc==0) then
      !!    do istat=1,tmb%orbs%norb
      !!        do iall=1,tmb%orbs%norb
      !!            write(600+iproc,*) istat, iall, tmb%linmat%denskern%matrix(iall,istat)
      !!        end do
      !!    end do
      !!end if

      !call compress_matrix_for_allreduce(iproc,tmb%linmat%denskern)
      !call uncompressMatrix(iproc,tmb%linmat%denskern)

  end do
  tmb%linmat%denskern%matrix=2.0d0*tmb%linmat%denskern%matrix
  !!if (iproc==0) then
  !!    do istat=1,tmb%orbs%norb
  !!        do iall=1,tmb%orbs%norb
  !!            write(700+iproc,*) istat, iall, tmb%linmat%denskern%matrix(iall,istat)
  !!        end do
  !!    end do
  !!end if 
  iall = -product(shape(ks))*kind(ks)
  deallocate(ks,stat=istat)
  call memocc(istat, iall, 'ks', subname)
  iall = -product(shape(ksk))*kind(ksk)
  deallocate(ksk,stat=istat)
  call memocc(istat, iall, 'ksk', subname)
  iall = -product(shape(ksksk))*kind(ksksk)
  deallocate(ksksk,stat=istat)
  call memocc(istat, iall, 'ksksk', subname)

  call compress_matrix_for_allreduce(iproc,tmb%linmat%denskern)

  iall=-product(shape(tmb%linmat%ovrlp%matrix))*kind(tmb%linmat%ovrlp%matrix)
  deallocate(tmb%linmat%ovrlp%matrix, stat=istat)
  call memocc(istat, iall, 'tmb%linmat%ovrlp%matrix', subname)
  iall=-product(shape(tmb%linmat%denskern%matrix))*kind(tmb%linmat%denskern%matrix)
  deallocate(tmb%linmat%denskern%matrix, stat=istat)
  call memocc(istat, iall, 'tmb%linmat%denskern%matrix', subname)

end subroutine purify_kernel
