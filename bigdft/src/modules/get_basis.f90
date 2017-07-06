module get_basis
  implicit none
  private

  !> Public routines
  public :: getLocalizedBasis
  public :: allocate_precond_arrays
  public :: deallocate_precond_arrays

  contains


    subroutine getLocalizedBasis(iproc,nproc,at,orbs,rxyz,denspot,GPU,trH,trH_old,&
            fnrm_tmb,infoBasisFunctions,nlpsp,scf_mode,ldiis,SIC,tmb,energs_base,do_iterative_orthogonalization,sf_per_type,&
        nit_precond,target_function,&
        correction_orthoconstraint,nit_basis,&
        ratio_deltas,ortho_on,extra_states,itout,conv_crit,experimental_mode,early_stop,&
        gnrm_dynamic, min_gnrm_for_dynamic, can_use_ham, order_taylor, max_inversion_error, kappa_conv, &
        correction_co_contra, &
        precond_convol_workarrays, precond_workarrays, &
        wt_philarge,  wt_hphi, wt_phi, fnrm, energs_work, frag_calc, reset_DIIS_history, &
        cdft, input_frag, ref_frags, hphi_pspandkin, eproj, ekin)
      !
      ! Purpose:
      ! ========
      !   Calculates the localized basis functions phi. These basis functions are obtained by adding a
      !   quartic potential centered on the atoms to the ordinary Hamiltonian. The eigenfunctions are then
      !   determined by minimizing the trace until the gradient norm is below the convergence criterion.
      use module_base
      use module_types
      use yaml_output
      use get_kernel, only: renormalize_kernel
      use module_interfaces, only: LocalHamiltonianApplication, SynchronizeHamiltonianApplication
      use io, only: write_energies
      use communications_base, only: work_transpose, TRANSPOSE_FULL, TRANSPOSE_POST, TRANSPOSE_GATHER
      use communications, only: transpose_localized, untranspose_localized, start_onesided_communication, &
                                synchronize_onesided_communication
      use rhopotential, only: full_local_potential
      use sparsematrix_base, only: assignment(=), sparsematrix_malloc, sparsematrix_malloc_ptr, SPARSE_FULL, &
                                   SPARSE_TASKGROUP
      use constrained_dft, only: cdft_data
      use fragment_base, only: fragmentInputParameters
      use module_fragments, only: system_fragment
      use sparsematrix,only: gather_matrix_from_taskgroups_inplace, extract_taskgroup_inplace
      use transposed_operations, only: calculate_overlap_transposed
      use matrix_operations, only: overlapPowerGeneral, check_taylor_order
      use public_enums
      use locreg_operations
      use locregs_init, only: small_to_large_locreg
      use coeffs, only: calculate_density_kernel
      !  use Poisson_Solver
      !use allocModule
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      integer,intent(inout) :: order_taylor
      real(kind=8),intent(in) :: max_inversion_error
      integer,intent(out) :: infoBasisFunctions
      type(atoms_data), intent(in) :: at
      type(orbitals_data), intent(inout) :: orbs
      real(kind=8),dimension(3,at%astruct%nat) :: rxyz
      type(DFT_local_fields), intent(inout) :: denspot
      type(GPU_pointers), intent(inout) :: GPU
      real(kind=8),intent(out) :: trH, fnrm_tmb
      real(kind=8),intent(inout) :: trH_old
      type(DFT_PSP_projectors),intent(inout) :: nlpsp
      integer,intent(in) :: scf_mode
      type(localizedDIISParameters),intent(inout) :: ldiis
      type(DFT_wavefunction),target,intent(inout) :: tmb
      type(SIC_data) :: SIC !<parameters for the SIC methods
      type(energy_terms),intent(in) :: energs_base
      logical,intent(in) :: do_iterative_orthogonalization
      integer,dimension(at%astruct%ntypes),intent(in) :: sf_per_type
      integer, intent(in) :: nit_precond, target_function, correction_orthoconstraint, nit_basis
      real(kind=8),intent(out) :: ratio_deltas
      logical, intent(inout) :: ortho_on
      integer, intent(in) :: extra_states
      integer,intent(in) :: itout
      real(kind=8),intent(in) :: conv_crit, early_stop, gnrm_dynamic, min_gnrm_for_dynamic, kappa_conv
      logical,intent(in) :: experimental_mode
      logical,intent(out) :: can_use_ham
      logical,intent(in) :: correction_co_contra
      type(workarrays_quartic_convolutions),dimension(tmb%orbs%norbp),intent(inout) :: precond_convol_workarrays
      type(workarr_precond),dimension(tmb%orbs%norbp),intent(inout) :: precond_workarrays
      type(work_transpose),intent(inout) :: wt_philarge, wt_hphi, wt_phi
      type(work_mpiaccumulate),intent(inout) :: fnrm, energs_work
      logical, intent(in) :: frag_calc, reset_DIIS_history
      !these must all be present together
      type(cdft_data),intent(inout),optional :: cdft
      type(fragmentInputParameters),optional,intent(in) :: input_frag
      type(system_fragment), dimension(:), optional, intent(in) :: ref_frags
      real(kind=8),dimension(tmb%ham_descr%npsidim_orbs),intent(inout),optional :: hphi_pspandkin
      real(kind=8),intent(inout),optional :: eproj, ekin
     
      ! Local variables
      integer :: iorb, it, it_tot, ncount, ncharge, ii, kappa_satur, nit_exit, ispin, jproc
      !integer :: jorb, nspin
      !real(kind=8),dimension(:),allocatable :: occup_tmp
      real(kind=8) :: meanAlpha, ediff_best, alpha_max, delta_energy, delta_energy_prev, ediff
      real(kind=8),dimension(:),allocatable :: alpha,fnrmOldArr,alphaDIIS,occup_tmp
      real(kind=8),dimension(:),allocatable :: hpsit_c_tmp, hpsit_f_tmp, hpsi_tmp, psidiff, tmparr1, tmparr2
      real(kind=8),dimension(:),allocatable :: delta_energy_arr, hpsi_noprecond, kernel_compr_tmp, kernel_best, hphi_nococontra
      logical :: energy_increased, overlap_calculated, energy_diff, energy_increased_previous, complete_reset, even
      logical :: calculate_inverse, allow_increase, recovered_old_kernel
      real(kind=8),dimension(:),pointer :: lhphiold, lphiold, hpsit_c, hpsit_f, hpsi_small
      type(energy_terms) :: energs
      real(kind=8), dimension(2):: reducearr
      real(gp) :: econf, dynamic_convcrit, kappa_mean
      real(kind=8) :: energy_first, trH_ref, charge, fnrm_old
      real(kind=8),dimension(3),save :: kappa_history
      integer,save :: nkappa_history
      logical :: auxiliary_arguments_present
      logical,save :: has_already_converged
      logical,dimension(7) :: exit_loop
      type(matrices) :: ovrlp_old
      integer :: iiorb, ilr, i, ist
      real(kind=8) :: max_error, mean_error
      integer,dimension(:),allocatable :: n3p
      character(len=6) :: label
      integer,dimension(3) :: power
    
    
    
      call f_routine(id='getLocalizedBasis')
    
      !!fnrm = work_mpiaccumulate_null()
      !!fnrm%ncount = 1
      !!call allocate_work_mpiaccumulate(fnrm)
    
      !!energs_work = work_mpiaccumulate_null()
      !!energs_work%ncount = 4
      !!call allocate_work_mpiaccumulate(energs_work)
    
      energs = energy_terms_null()
      delta_energy_arr=f_malloc(nit_basis+6,id='delta_energy_arr')
      kernel_best=sparsematrix_malloc(tmb%linmat%smat(3),iaction=SPARSE_TASKGROUP,id='kernel_best')
      energy_diff=.false.
    
      ovrlp_old%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%smat(1), &
                               iaction=SPARSE_TASKGROUP, id='ovrlp_old%matrix_compr')
    
      ! Allocate all local arrays.
      call allocateLocalArrays()
    
    
      call timing(iproc,'getlocbasinit','ON')
      tmb%can_use_transposed=.false.
    
      alpha=ldiis%alphaSD
      alphaDIIS=ldiis%alphaDIIS
      ldiis%resetDIIS=.false.
      ldiis%immediateSwitchToSD=.false.
      allow_increase=.false.
    
     
      call timing(iproc,'getlocbasinit','OF')
    
      overlap_calculated=.false.
      it=0
      it_tot=0
      !ortho=.true.
      call local_potential_dimensions(iproc,tmb%ham_descr%lzd,tmb%orbs,denspot%xc,denspot%dpbox%ngatherarr(0,1))
      n3p = f_malloc(0.to.nproc-1,id='n3p')
      do jproc=0,nproc-1
          n3p(jproc) = max(denspot%dpbox%nscatterarr(jproc,2),1)
      end do
      call start_onesided_communication(iproc, nproc, denspot%dpbox%mesh%ndims(1), denspot%dpbox%mesh%ndims(2), &
           n3p, denspot%rhov, &
           tmb%ham_descr%comgp%nrecvbuf*tmb%ham_descr%comgp%nspin, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, &
           tmb%ham_descr%lzd)
      call f_free(n3p)
    
      delta_energy_prev=1.d100
    
      energy_increased_previous=.false.
      ratio_deltas=1.d0
      ediff_best=1.d0
      ediff=1.d0
      delta_energy_prev=1.d0
      delta_energy_arr=1.d0
      trH_ref=trH_old
      dynamic_convcrit=1.d-100
      kappa_satur=0
      fnrm_old = 0.d0
    
    
      ! Count whether there is an even or an odd number of electrons
      charge=0.d0
      do iorb=1,orbs%norb
          charge=charge+orbs%occup(iorb)
      end do
      ncharge=nint(charge)
      even=(mod(ncharge,2)==0)
    
      !!! Purify the initial kernel (only when necessary and if there is an even number of electrons)
      !!if (target_function/=TARGET_FUNCTION_IS_TRACE .and. even .and. &
      !!    (scf_mode==LINEAR_FOE .or. scf_mode==LINEAR_PEXSI)) then
      !!    if (iproc==0) then
      !!        call yaml_sequence(advance='no')
      !!        call yaml_mapping_open(flow=.true.)
      !!        call yaml_map('Initial kernel purification',.true.)
      !!    end if
      !!    overlap_calculated=.true.
      !!    do ispin=1,tmb%linmat%smat(3)%nspin
      !!        call purify_kernel(iproc, nproc, tmb, overlap_calculated, 1, 30, order_taylor, &
      !!             max_inversion_error, purification_quickreturn, ispin)
      !!    end do
      !!    if (iproc==0) call yaml_mapping_close()
      !!end if
    
      if (itout==0) then
          nkappa_history=0
          kappa_history=0.d0
          has_already_converged=.false.
      end if
    
      ! Check the optional arguments
      if (any((/present(hphi_pspandkin),present(eproj),present(ekin)/))) then
          if (all((/present(hphi_pspandkin),present(eproj),present(ekin)/))) then
              auxiliary_arguments_present = .true.
          else
              call f_err_throw('The arguments hphi_pspandkin, eproj and ekin miust be present at the same time',&
                   err_name='BIGDFT_RUNTIME_ERROR')
          end if
      else
          auxiliary_arguments_present = .false.
      end if
    
      energy_increased = .false.
      recovered_old_kernel = .false.
      complete_reset = .false.
    
      iterLoop: do
    
    
          it=it+1
          it=max(it,1) !since it could become negative (2 is subtracted if the loop cycles)
          it_tot=it_tot+1
    
          fnrm%sendbuf(1)=0.d0
          fnrm%receivebuf(1)=0.d0
      
          if (iproc==0) then
              call yaml_sequence(advance='no')
              call yaml_mapping_open(flow=.true.)
              call yaml_comment('iter:'//yaml_toa(it,fmt='(i6)'),hfill='-')
              !!if (target_function==TARGET_FUNCTION_IS_TRACE) then
              !!    call yaml_map('target function','TRACE')
              !!else if (target_function==TARGET_FUNCTION_IS_ENERGY) then
              !!    call yaml_map('target function','ENERGY')
              !!else if (target_function==TARGET_FUNCTION_IS_HYBRID) then
              !!    call yaml_map('target function','HYBRID')
              !!end if
              ! Reset the DIIS history
              if (it_tot==1) then
                  if (reset_DIIS_history) then
                      ldiis%is = 0
                      if (iproc==0) call yaml_map('reset DIIS history',.true.)
                  else
                      if (iproc==0) call yaml_map('reset DIIS history',.false.)
                  end if
              end if
          end if
    
          ! Synchronize the mpi_get before starting a new communication
          call synchronize_onesided_communication(iproc, nproc, tmb%ham_descr%comgp)
    
          ! Start the communication
          call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
               TRANSPOSE_POST, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd, wt_phi)
    
          ! Calculate the unconstrained gradient by applying the Hamiltonian.
          if (tmb%ham_descr%npsidim_orbs > 0)  call f_zero(tmb%ham_descr%npsidim_orbs,tmb%hpsi(1))
          call small_to_large_locreg(iproc, tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%isorb, tmb%orbs%inwhichlocreg, &
               tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
               tmb%psi, tmb%ham_descr%psi)
    
          ! Start the nonblocking transposition (the results will be gathered in
          ! orthoconstraintNonorthogonal)
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               TRANSPOSE_POST, tmb%ham_descr%psi, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, tmb%ham_descr%lzd, &
               wt_philarge)
    
          call NonLocalHamiltonianApplication(iproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
               tmb%ham_descr%lzd,nlpsp,tmb%ham_descr%psi,tmb%hpsi,energs%eproj,tmb%paw)
          ! only kinetic because waiting for communications
          call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
               tmb%ham_descr%lzd,tmb%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,&
               & tmb%ham_descr%psi,tmb%hpsi,energs,SIC,GPU,3,denspot%xc,&
               & pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
               & potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
          !!if (auxiliary_arguments_present .and. &
          !!    (target_function==TARGET_FUNCTION_IS_ENERGY .or. target_function==TARGET_FUNCTION_IS_HYBRID)) then
          if (auxiliary_arguments_present) then
              if (tmb%ham_descr%npsidim_orbs > 0) then
                  call f_memcpy(src=tmb%hpsi, dest=hphi_pspandkin)
                  eproj = energs%eproj
                  ekin = energs%ekin
              end if
          end if
          call full_local_potential(iproc,nproc,tmb%orbs,tmb%ham_descr%lzd,2,denspot%dpbox,&
               & denspot%xc,denspot%rhov,denspot%pot_work,tmb%ham_descr%comgp)
          ! only potential
          if (target_function==TARGET_FUNCTION_IS_HYBRID) then
              call vcopy(tmb%ham_descr%npsidim_orbs, tmb%hpsi(1), 1, hpsi_tmp(1), 1)
              call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
                   tmb%ham_descr%lzd,tmb%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,&
                   & tmb%ham_descr%psi,tmb%hpsi,energs,SIC,GPU,2,denspot%xc,&
                   & pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
                   & potential=denspot%rhov,comgp=tmb%ham_descr%comgp,&
                   hpsi_noconf=hpsi_tmp,econf=econf)
    
              !!if (nproc>1) then
              !!    call mpiallred(econf, 1, mpi_sum, bigdft_mpi%mpi_comm)
              !!end if
    
          else
              call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
                   tmb%ham_descr%lzd,tmb%confdatarr,denspot%dpbox%ngatherarr,&
                   & denspot%pot_work,tmb%ham_descr%psi,tmb%hpsi,energs,SIC,GPU,2,denspot%xc,&
                   & pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
                   & potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
          end if
    
    
          !!if (target_function==TARGET_FUNCTION_IS_HYBRID .and. iproc==0) then
          !!    write(*,*) 'econf, econf/tmb%orbs%norb',econf, econf/tmb%orbs%norb
          !!end if
    
          call timing(iproc,'glsynchham2','ON')
          call SynchronizeHamiltonianApplication(nproc,tmb%ham_descr%npsidim_orbs,&
               tmb%orbs,tmb%ham_descr%lzd,GPU,denspot%xc,tmb%hpsi,&
               energs,energs_work)
          call timing(iproc,'glsynchham2','OF')
    
          if (iproc==0) then
              call yaml_map('Hamiltonian Applied',.true.)
          end if
    
          !if (iproc==0) write(*,'(a,5es16.6)') 'ekin, eh, epot, eproj, eex', &
          !              energs%ekin, energs%eh, energs%epot, energs%eproj, energs%exc
    
    
          ! Start the communication
          if (target_function==TARGET_FUNCTION_IS_HYBRID) then
              call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
                   TRANSPOSE_POST, hpsi_tmp, hpsit_c, hpsit_f, tmb%ham_descr%lzd, wt_hphi)
          else
              call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
                   TRANSPOSE_POST, tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd, wt_hphi)
          end if
    
          ! Gather the data
          call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
               TRANSPOSE_GATHER, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd, wt_phi)
    
          if (nproc>1) then
              ! Wait for the communication of energs_work on root
              call mpi_fenceandfree(energs_work%window)
          end if
    
          ! Copy the value, only necessary on root
          if (iproc==0) then
              energs%ekin = energs_work%receivebuf(1)
              energs%epot = energs_work%receivebuf(2)
              energs%eproj = energs_work%receivebuf(3)
              energs%evsic = energs_work%receivebuf(4)
          end if
    
          ! Use this subroutine to write the energies, with some fake number
          ! to prevent it from writing too much
          if (iproc==0) then
              call write_energies(0,energs,0.d0,0.d0,'',only_energies=.true.,label='Components')
          end if
          if (iproc==0) then
              call yaml_map('Orthoconstraint',.true.)
          end if
    
          !if (target_function==TARGET_FUNCTION_IS_HYBRID .and. .not.energy_increased) then
          ! Only need to renormalize the kernel if it is actually used.
          if (target_function==TARGET_FUNCTION_IS_HYBRID .or. target_function==TARGET_FUNCTION_IS_ENERGY) then
              tmb%ham_descr%can_use_transposed=.false.
              call f_memcpy(src=tmb%linmat%ovrlp_%matrix_compr, dest=ovrlp_old%matrix_compr)
    
              call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, &
                   tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%smat(1), tmb%linmat%auxs, tmb%linmat%ovrlp_)
              !if (iproc==0) call yaml_newline()
              !if (iproc==0) call yaml_sequence_open('kernel update by renormalization')
              if (it==1 .or. energy_increased .or. .not.experimental_mode .or. complete_reset) then
                  ! Calculate S^1/2, as it can not be taken from memory
                  power(1)=2
                  !if (iproc==0) call yaml_warning('call overlapPowerGeneral')
                  !!if (iproc==0) write(*,*) 'call overlapPowerGeneral'
                  call overlapPowerGeneral(iproc, nproc,bigdft_mpi%mpi_comm,&
                       order_taylor, 1, power(1), -1, &
                       imode=1, ovrlp_smat=tmb%linmat%smat(1), inv_ovrlp_smat=tmb%linmat%smat(3), &
                       ovrlp_mat=ovrlp_old, inv_ovrlp_mat=tmb%linmat%ovrlppowers_(1), &
                       verbosity=0, &
                       check_accur=order_taylor<1000, max_error=max_error, mean_error=mean_error, &
                       ice_obj=tmb%ice_obj)
                  call check_taylor_order(iproc, mean_error, max_inversion_error, order_taylor)
              end if
              !if (.not.energy_increased) then
              if (.not.recovered_old_kernel .and. .not.complete_reset) then
                  !if (iproc==0) call yaml_warning('call renormalize_kernel')
                  !!if (iproc==0) write(*,*) 'sum(S), sum(Sold), sum(S-1)',&
                  !!sum(tmb%linmat%ovrlp_%matrix_compr), sum(ovrlp_old%matrix_compr), sum(tmb%linmat%ovrlppowers_(1)%matrix_compr)
                  !!if (iproc==0) write(*,*) 'before renorm: sum(K)', sum(tmb%linmat%kernel_%matrix_compr)
                  call renormalize_kernel(iproc, nproc, order_taylor, max_inversion_error, tmb, tmb%linmat%ovrlp_, ovrlp_old)
                  !!if (iproc==0) write(*,*) 'after renorm: sum(K)', sum(tmb%linmat%kernel_%matrix_compr)
              else
                  ! Calculate S^1/2 for the overlap matrix
                   power=(/2,-2,1/)
                   call overlapPowerGeneral(iproc, nproc, bigdft_mpi%mpi_comm, &
                        order_taylor, 3, power, -1, &
                        imode=1, ovrlp_smat=tmb%linmat%smat(1), inv_ovrlp_smat=tmb%linmat%smat(3), &
                        ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=tmb%linmat%ovrlppowers_, &
                        verbosity=0, &
                        check_accur=order_taylor<1000, max_error=max_error, mean_error=mean_error, &
                        ice_obj=tmb%ice_obj)
                   call check_taylor_order(iproc, mean_error, max_inversion_error, order_taylor)
              end if
          end if
    
          ! Gather the data
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               TRANSPOSE_GATHER, hpsi_tmp, hpsit_c, hpsit_f, tmb%ham_descr%lzd, wt_hphi)
          ncount=tmb%ham_descr%collcom%ndimind_c
          if(ncount>0) call vcopy(ncount, hpsit_c(1), 1, hpsit_c_tmp(1), 1)
          ncount=7*tmb%ham_descr%collcom%ndimind_f
          if(ncount>0) call vcopy(ncount, hpsit_f(1), 1, hpsit_f_tmp(1), 1)
    
          ! optimize the tmbs for a few extra states
          if (target_function==TARGET_FUNCTION_IS_ENERGY.and.extra_states>0) then
              kernel_compr_tmp = sparsematrix_malloc(tmb%linmat%smat(3), iaction=SPARSE_TASKGROUP, id='kernel_compr_tmp')
              !call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(3), tmb%linmat%kernel_)
              call vcopy(tmb%linmat%smat(3)%nvctrp_tg*tmb%linmat%smat(3)%nspin, &
                   tmb%linmat%kernel_%matrix_compr(1), 1, kernel_compr_tmp(1), 1)
              occup_tmp=f_malloc(orbs%norb,id='occup_tmp')
              ! make a copy of 'correct' occup
              call vcopy(orbs%norb, orbs%occup(1), 1, occup_tmp(1), 1)
              ! occupy the extra states - don't need to preserve the charge as only using for support function optimization
              do iorb=1,orbs%norb
                 orbs%occup(iorb)=2.0_gp
              end do
              call calculate_density_kernel(iproc, nproc, bigdft_mpi%mpi_comm, .true., &
                   orbs%norbp, orbs%isorb, orbs%norbu, orbs%norb, orbs%occup, &
                   tmb%coeff, tmb%linmat%smat(3), tmb%linmat%kernel_)
              call vcopy(orbs%norb, occup_tmp(1), 1, orbs%occup(1), 1)
              call f_free(occup_tmp)
       
              !call extract_taskgroup_inplace(tmb%linmat%smat(3), tmb%linmat%kernel_)
              !call transform_sparse_matrix(tmb%linmat%denskern, tmb%linmat%denskern_large, 'large_to_small')
          end if
    
          ! use hpsi_tmp as temporary array for hpsi_noprecond, even if it is allocated with a larger size
          !write(*,*) 'calling calc_energy_and.., correction_co_contra',correction_co_contra
          calculate_inverse = (target_function/=TARGET_FUNCTION_IS_HYBRID) .or. energy_increased
          !!call extract_taskgroup_inplace(tmb%linmat%smat(3), tmb%linmat%kernel_)
          call calculate_energy_and_gradient_linear(iproc, nproc, it, ldiis, fnrmOldArr, &
               fnrm_old, alpha, trH, trH_old, fnrm, &
               meanAlpha, alpha_max, energy_increased, tmb, lhphiold, overlap_calculated, energs_base, &
               hpsit_c, hpsit_f, nit_precond, target_function, correction_orthoconstraint, hpsi_small, &
               experimental_mode, calculate_inverse, &
               correction_co_contra, recovered_old_kernel, hpsi_noprecond=hpsi_tmp, norder_taylor=order_taylor, &
               max_inversion_error=max_inversion_error, &
               precond_convol_workarrays=precond_convol_workarrays, precond_workarrays=precond_workarrays, &
               wt_hphi=wt_hphi, wt_philarge=wt_philarge, &
               cdft=cdft, input_frag=input_frag, ref_frags=ref_frags)
          !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(3), tmb%linmat%kernel_)
          !fnrm_old=fnrm
    
    
          if (it_tot==1) then
              energy_first=trH
          end if
          if (experimental_mode) then
              if (iproc==0) call yaml_map('rel D',(trH-energy_first)/energy_first,fmt='(es9.2)')
              if ((trH-energy_first)<0.d0 .and. abs((trH-energy_first)/energy_first)>early_stop .and. itout>0) then
                  energy_diff=.true.
              end if
          end if
    
          if (target_function==TARGET_FUNCTION_IS_ENERGY.and.extra_states>0) then
              !if (tmb%linmat%smat(3)%nspin>1) stop 'THIS IS NOT TESTED FOR SPIN POLARIZED SYSTEMS!'
              call vcopy(tmb%linmat%smat(3)%nvctrp_tg*tmb%linmat%smat(3)%nspin, kernel_compr_tmp(1), 1, &
                   tmb%linmat%kernel_%matrix_compr(1), 1)
              call f_free(kernel_compr_tmp)
          end if
    
          ediff=trH-trH_old
          ediff_best=trH-trH_ref
    
          if (it>1 .and. (target_function==TARGET_FUNCTION_IS_HYBRID .or. experimental_mode)) then
              if (.not.energy_increased .and. .not.energy_increased_previous) then
                  if (.not.ldiis%switchSD) then
                      ratio_deltas=ediff_best/delta_energy_prev
                  else
                      ratio_deltas=ediff_best/delta_energy_arr(ldiis%itBest)
                  end if
              else
                  ! use a default value
                  if (iproc==0) then
                      call yaml_warning('use a fake value for kappa')
                      call yaml_newline()
                  end if
                  ratio_deltas=0.5d0
              end if
              if (ldiis%switchSD) then
                  !!ratio_deltas=0.5d0
                  !!if (iproc==0) write(*,*) 'WARNING: TEMPORARY FIX for ratio_deltas!'
              end if
              if (iproc==0) call yaml_map('kappa',ratio_deltas,fmt='(es10.3)')
              if (target_function==TARGET_FUNCTION_IS_HYBRID) then
                  !if (ratio_deltas>0.d0) then
                  !if (ratio_deltas>1.d-12) then
                  if (.not.energy_increased .and. .not.energy_increased_previous) then
                      if (iproc==0) call yaml_map('kappa to history',.true.)
                      nkappa_history=nkappa_history+1
                      ii=mod(nkappa_history-1,3)+1
                      kappa_history(ii)=ratio_deltas
                  else
                      if (iproc==0) call yaml_map('kappa to history',.false.)
                  end if
                  !!if (nkappa_history>=3) then
                  !!    kappa_mean=sum(kappa_history)/3.d0
                  !!    if (iproc==0) call yaml_map('mean kappa',kappa_mean,fmt='(es10.3)')
                  !!    dynamic_convcrit=conv_crit/kappa_mean
                  !!    if (iproc==0) call yaml_map('dynamic conv crit',dynamic_convcrit,fmt='(es9.2)')
                  !!end if
              end if
          end if
          if (target_function==TARGET_FUNCTION_IS_HYBRID) then
              if (nkappa_history>=3) then
                  kappa_mean=sum(kappa_history)/3.d0
                  if (iproc==0) call yaml_map('mean kappa',kappa_mean,fmt='(es10.3)')
                  !dynamic_convcrit=conv_crit/kappa_mean
                  dynamic_convcrit=gnrm_dynamic/kappa_mean
                  if (iproc==0) call yaml_map('dynamic conv crit',dynamic_convcrit,fmt='(es9.2)')
              end if
          end if
    
          if (energy_increased) then
              energy_increased_previous=.true.
          else
              energy_increased_previous=.false.
          end if
    
    
    
          !!delta_energy_prev=delta_energy
    
          ! Wait for the communication of fnrm on root
          if (nproc>1) then
              call mpi_fenceandfree(fnrm%window)
          end if
          fnrm%receivebuf(1)=sqrt(fnrm%receivebuf(1)/dble(tmb%orbs%norb))
    
          ! The other processes need to get fnrm as well. The fence will be later as only iproc=0 has to write.
          if (nproc>1) then
              if (iproc==0) fnrm%sendbuf(1) = fnrm%receivebuf(1)
              fnrm%window = mpiwindow(1, fnrm%sendbuf(1), bigdft_mpi%mpi_comm)
              if (iproc/=0) then
                  call mpiget(fnrm%receivebuf(1), 1, 0, int(0,kind=mpi_address_kind), fnrm%window)
              end if
          end if
    
          if (energy_increased .and. ldiis%isx==0 .and. (.not. allow_increase)) then
              !if (iproc==0) write(*,*) 'WARNING: ENERGY INCREASED'
              !if (iproc==0) call yaml_warning('The target function increased, D='&
              !              //trim(adjustl(yaml_toa(trH-ldiis%trmin,fmt='(es10.3)'))))
              if (nproc>1) then
                  call mpi_fenceandfree(fnrm%window)
              end if
              fnrm_old=fnrm%receivebuf(1)
              if (iproc==0) then
                  call yaml_newline()
                  call yaml_map('iter',it,fmt='(i5)')
                  call yaml_map('fnrm',fnrm%receivebuf(1),fmt='(es9.2)')
                  call yaml_mapping_open('Omega')
                  select case (target_function)
                  case (TARGET_FUNCTION_IS_TRACE)
                      label = 'TRACE'
                  case (TARGET_FUNCTION_IS_ENERGY)
                      label = 'ENERGY'
                  case (TARGET_FUNCTION_IS_HYBRID)
                      label = 'HYBRID'
                  case default
                      call f_err_throw('wrong target function', err_name='BIGDFT_RUNTIME_ERROR')
                  end select
                  call yaml_map(trim(label),trH,fmt='(es22.15)')
                  call yaml_mapping_close()
                  call yaml_map('D',ediff,fmt='(es9.2)')
                  call yaml_map('D best',ediff_best,fmt='(es9.2)')
              end if
              tmb%ham_descr%can_use_transposed=.false.
              call vcopy(tmb%npsidim_orbs, lphiold(1), 1, tmb%psi(1), 1)
              can_use_ham=.false.
              call vcopy(tmb%linmat%smat(3)%nvctrp_tg*tmb%linmat%smat(3)%nspin, kernel_best(1), 1, &
                   tmb%linmat%kernel_%matrix_compr(1), 1)
              if (iproc==0) then
                  call yaml_map('Recovering old support functions and kernel',.true.)
              end if
              recovered_old_kernel = .true.
              !ldiis%switchSD = .true.
              !write(*,*) 'cut alpha 0.6 main'
              !alpha(:) = alpha(:)*0.6d0
              !if (iproc==0) call yaml_warning('set recovered_old_kernel to true')
    
    
              ! Recalculate the matrix powers
              power=(/2,-2,1/)
              call overlapPowerGeneral(iproc, nproc, bigdft_mpi%mpi_comm, &
                   order_taylor, 3, power, -1, &
                   imode=1, ovrlp_smat=tmb%linmat%smat(1), inv_ovrlp_smat=tmb%linmat%smat(3), &
                   ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=tmb%linmat%ovrlppowers_, &
                   verbosity=0, &
                   check_accur=order_taylor<1000, max_error=max_error, mean_error=mean_error, &
                   ice_obj=tmb%ice_obj)
              call check_taylor_order(iproc, mean_error, max_inversion_error, order_taylor)
    
              trH_old=0.d0
              it=it-2 !go back one iteration (minus 2 since the counter was increased)
              overlap_calculated=.false.
              ! print info here anyway for debugging
              if (it_tot<2*nit_basis) then ! just in case the step size is the problem
                  call yaml_mapping_close()
                  call yaml_flush_document()
                  !call bigdft_utils_flush(unit=6)
                  ! This is to avoid memory leaks
                  call untranspose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
                       TRANSPOSE_GATHER, hpsit_c, hpsit_f, hpsi_tmp, tmb%ham_descr%lzd, wt_philarge)
    
                  ! for fragment calculations tmbs may be far from orthonormality so allow an increase in energy for a few iterations
                  if (frag_calc .and. itout<4) then
                     allow_increase=.true.
                  end if
                  cycle
              else if(it_tot<3*nit_basis .and. .not.experimental_mode) then ! stop orthonormalizing the tmbs
                 if (iproc==0) call yaml_newline()
                 if (iproc==0) call yaml_warning('Energy increasing, switching off orthonormalization of tmbs')
                 ortho_on=.false.
                 alpha=alpha*5.0d0/3.0d0 ! increase alpha to make up for decrease from previous iteration
              end if
          else
              can_use_ham=.true.
              recovered_old_kernel = .false.
              !if (iproc==0) call yaml_warning('set recovered_old_kernel to false')
          end if 
    
    
          ! information on the progress of the optimization
          if (iproc==0) then
              call yaml_newline()
              call yaml_map('iter',it,fmt='(i5)')
              call yaml_map('fnrm',fnrm%receivebuf(1),fmt='(es9.2)')
              !call yaml_map('Omega',trH,fmt='(es22.15)')
              select case (target_function)
              case (TARGET_FUNCTION_IS_TRACE)
                  label = 'TRACE'
              case (TARGET_FUNCTION_IS_ENERGY)
                  label = 'ENERGY'
              case (TARGET_FUNCTION_IS_HYBRID)
                  label = 'HYBRID'
              case default
                  call f_err_throw('wrong target function', err_name='BIGDFT_RUNTIME_ERROR')
              end select
              call yaml_mapping_open('Omega')
              call yaml_map(trim(label),trH,fmt='(es22.15)')
              call yaml_mapping_close()
              call yaml_map('D',ediff,fmt='(es9.2)')
              call yaml_map('D best',ediff_best,fmt='(es9.2)')
          end if
    
          ! Add some extra iterations if DIIS failed (max 6 failures are allowed before switching to SD)
          nit_exit=min(nit_basis+ldiis%icountDIISFailureTot,nit_basis+6)
    
          ! Normal case
          if (.not.energy_increased .or. ldiis%isx/=0 .or. allow_increase) then
              if (nproc>1) then
                  call mpi_fenceandfree(fnrm%window)
              end if
              fnrm_old=fnrm%receivebuf(1)
          end if
    
          ! Determine whether the loop should be exited
          exit_loop(1) = (it>=nit_exit)
          exit_loop(2) = (it_tot>=3*nit_basis)
          exit_loop(3) = energy_diff
          exit_loop(4) = (fnrm%receivebuf(1)<conv_crit .and. experimental_mode)
          exit_loop(5) = (experimental_mode .and. fnrm%receivebuf(1)<dynamic_convcrit &
                         .and. fnrm%receivebuf(1)<min_gnrm_for_dynamic &
                         .and. (it>1 .or. has_already_converged)) ! first overall convergence not allowed in a first iteration
          exit_loop(6) = (itout==0 .and. it>1 .and. ratio_deltas<kappa_conv .and.  ratio_deltas>0.d0)
          if (ratio_deltas>0.d0 .and. ratio_deltas<1.d-1) then
              kappa_satur=kappa_satur+1
          else
              kappa_satur=0
          end if
          exit_loop(7) = (.false. .and. itout>0 .and. kappa_satur>=2)
    
          if(any(exit_loop)) then
              if(exit_loop(1)) then
                  infoBasisFunctions=-1
                  if(iproc==0) call yaml_map('exit criterion','net number of iterations')
              end if
              if (exit_loop(2)) then
                  infoBasisFunctions=-2
                  if (iproc==0) call yaml_map('exit criterion','total number of iterations')
              end if
              if (exit_loop(3)) then
                  infoBasisFunctions=it
                  if (iproc==0) call yaml_map('exit criterion','energy difference')
              end if
              if (exit_loop(4)) then
                  if (iproc==0) call yaml_map('exit criterion','gradient')
                  infoBasisFunctions=it
              end if
              if (exit_loop(5)) then
                  if (iproc==0) call yaml_map('exit criterion','dynamic gradient')
                  infoBasisFunctions=it
                  has_already_converged=.true.
              end if
              if (exit_loop(6)) then
                  infoBasisFunctions=it
                  if (iproc==0) call yaml_map('exit criterion','extended input guess')
              end if
              if (exit_loop(7)) then
                  infoBasisFunctions=it
                  if (iproc==0) call yaml_map('exit criterion','kappa')
              end if
              if (can_use_ham) then
                  ! Calculate the Hamiltonian matrix, since we have all quantities ready. This matrix can then be used in the first
                  ! iteration of get_coeff.
                  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
                       tmb%ham_descr%psit_c, hpsit_c_tmp, tmb%ham_descr%psit_f, hpsit_f_tmp, &
                       tmb%linmat%smat(2), tmb%linmat%auxm, tmb%linmat%ham_)
                  !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(2), tmb%linmat%ham_)
              end if
    
              if (iproc==0) then
                  !yaml output
                  call yaml_mapping_close() !iteration
                  call yaml_flush_document()
                  !call bigdft_utils_flush(unit=6)
              end if
    
              ! This is to avoid memory leaks
              ! Gather together the data (was posted in orthoconstraintNonorthogonal)
              ! Give hpsit_c and hpsit_f, this should not matter if GATHER is specified.
              ! To be modified later
              call untranspose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
                   TRANSPOSE_GATHER, hpsit_c, hpsit_f, hpsi_tmp, tmb%ham_descr%lzd, wt_philarge)
    
              exit iterLoop
          end if
          trH_old=trH
    
          if (ldiis%isx>0) then
              ldiis%mis=mod(ldiis%is,ldiis%isx)+1 !to store the energy at the correct location in the history
          end if
          call hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, at, do_iterative_orthogonalization, sf_per_type, &
               lphiold, alpha, trH, meanAlpha, alpha_max, alphaDIIS, hpsi_small, ortho_on, psidiff, &
               experimental_mode, order_taylor, max_inversion_error, trH_ref, kernel_best, complete_reset)
    
    
          overlap_calculated=.false.
          ! It is now not possible to use the transposed quantities, since they have changed.
          if(tmb%ham_descr%can_use_transposed) then
              tmb%ham_descr%can_use_transposed=.false.
          end if
    
          ! Gather together the data (was posted in orthoconstraintNonorthogonal)
          ! Give hpsit_c and hpsit_f, this should not matter if GATHER is specified.
          ! To be modified later
          hphi_nococontra = f_malloc(tmb%ham_descr%npsidim_orbs,id='hphi_nococontra')
          call untranspose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               TRANSPOSE_GATHER, hpsit_c, hpsit_f, hphi_nococontra, tmb%ham_descr%lzd, wt_philarge)
          call large_to_small_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
               tmb%orbs, hphi_nococontra, hpsi_tmp)
          call f_free(hphi_nococontra)
    
    
          ! Estimate the energy change, that is to be expected in the next optimization
          ! step, given by the product of the force and the "displacement" .
          if (target_function==TARGET_FUNCTION_IS_HYBRID .or. experimental_mode) then
              call estimate_energy_change(tmb%npsidim_orbs, tmb%orbs, tmb%lzd, tmb%linmat%smat(3)%nspin, psidiff, &
                   hpsi_tmp,delta_energy)
              ! This is a hack...
              if (energy_increased) then
                  delta_energy=1.d100
                  !ratio_deltas=1.d100
              end if
              !if (iproc==0) write(*,*) 'delta_energy', delta_energy
              delta_energy_prev=delta_energy
              delta_energy_arr(max(it,1))=delta_energy !max since the counter was decreased if there are problems, might lead to wrong results otherwise
          end if
    
    
          ! Only need to reconstruct the kernel if it is actually used.
          !SM: Do we really need this here? We have already renormalized the kernel above in renormalize_kernel...
          !!if ((target_function/=TARGET_FUNCTION_IS_TRACE .or. scf_mode==LINEAR_DIRECT_MINIMIZATION) &
          !!     .and. .not.complete_reset ) then
          !!    if(scf_mode/=LINEAR_FOE .and. scf_mode/=LINEAR_PEXSI) then
          !!        call reconstruct_kernel(iproc, nproc, order_taylor, tmb%orthpar%blocksize_pdsyev, &
          !!             tmb%orthpar%blocksize_pdgemm, orbs, tmb, overlap_calculated)
          !!        if (iproc==0) call yaml_map('reconstruct kernel',.true.)
          !!    else if (experimental_mode .and. .not.complete_reset) then
          !!    end if
          !!end if
    
          if (iproc==0) then
              call yaml_mapping_close() !iteration
              call yaml_flush_document()
              !call bigdft_utils_flush(unit=6)
          end if
    
    
      end do iterLoop
    
      ! Write the final results
      if (iproc==0) then
          call yaml_sequence(label='final_supfun'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)'))),advance='no')
          call yaml_mapping_open(flow=.true.)
          call yaml_comment('iter:'//yaml_toa(it,fmt='(i6)'),hfill='-')
          !!if (target_function==TARGET_FUNCTION_IS_TRACE) then
          !!    call yaml_map('target function','TRACE')
          !!else if (target_function==TARGET_FUNCTION_IS_ENERGY) then
          !!    call yaml_map('target function','ENERGY')
          !!else if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          !!    call yaml_map('target function','HYBRID')
          !!end if
          call write_energies(0,energs,0.d0,0.d0,'',only_energies=.true.,label='Components')
          call yaml_newline()
          call yaml_map('nit',it,fmt='(i5)')
          call yaml_map('fnrm',fnrm%receivebuf(1),fmt='(es9.2)')
          !call yaml_map('Omega',trH,fmt='(es22.15)')
          select case (target_function)
          case (TARGET_FUNCTION_IS_TRACE)
              label = 'TRACE'
          case (TARGET_FUNCTION_IS_ENERGY)
              label = 'ENERGY'
          case (TARGET_FUNCTION_IS_HYBRID)
              label = 'HYBRID'
          case default
              call f_err_throw('wrong target function', err_name='BIGDFT_RUNTIME_ERROR')
          end select
          call yaml_mapping_open('Omega')
          call yaml_map(trim(label),trH,fmt='(es22.15)')
          call yaml_mapping_close()
          !call yaml_map('D',ediff,fmt='(es9.2)')
          !call yaml_map('D best',ediff_best,fmt='(es9.2)')
          call yaml_map('D total',trH-energy_first,fmt='(es9.2)')
          call yaml_mapping_close() !iteration
          call yaml_flush_document()
          !call bigdft_utils_flush(unit=6)
      end if
    
    
      !!if (iproc==0) then
      !!    call yaml_comment('Support functions created')
      !!end if
    
    
      ! Deallocate potential
      call f_free_ptr(denspot%pot_work)
    
    
      ! Keep the values for the next iteration
      reducearr(1)=0.d0
      reducearr(2)=0.d0
      do iorb=1,tmb%orbs%norbp
          reducearr(1)=reducearr(1)+alpha(iorb)
          reducearr(2)=reducearr(2)+alphaDIIS(iorb)
      end do
    
      if (nproc > 1) then
          call mpiallred(reducearr, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
    
      reducearr(1)=reducearr(1)/dble(tmb%orbs%norb)
      reducearr(2)=reducearr(2)/dble(tmb%orbs%norb)
    
      ldiis%alphaSD=reducearr(1)
      ldiis%alphaDIIS=reducearr(2)
    
    
      ! Deallocate all local arrays.
      call deallocateLocalArrays()
      call f_free(delta_energy_arr)
      call f_free(kernel_best)
      call f_free_ptr(ovrlp_old%matrix_compr)
    
      fnrm_tmb = fnrm%receivebuf(1)
      !!call deallocate_work_mpiaccumulate(fnrm)
      !!call deallocate_work_mpiaccumulate(energs_work)
    
      call f_release_routine()
    
    contains
    
    
        subroutine allocateLocalArrays()
        !
        ! Purpose:
        ! ========
        !   This subroutine allocates all local arrays.
        !
        logical :: with_confpot
        integer :: iiorb, ilr, ncplx
        real(gp) :: kx, ky, kz
    
          alpha = f_malloc(tmb%orbs%norbp,id='alpha')
          alphaDIIS = f_malloc(tmb%orbs%norbp,id='alphaDIIS')
          fnrmOldArr = f_malloc(tmb%orbs%norbp,id='fnrmOldArr')
          hpsi_small = f_malloc_ptr(max(tmb%npsidim_orbs, tmb%npsidim_comp),id='hpsi_small')
          lhphiold = f_malloc_ptr(max(tmb%npsidim_orbs, tmb%npsidim_comp),id='lhphiold')
          lphiold = f_malloc_ptr(size(tmb%psi),id='lphiold')
          hpsit_c = f_malloc_ptr(tmb%ham_descr%collcom%ndimind_c,id='hpsit_c')
          hpsit_f = f_malloc_ptr(7*tmb%ham_descr%collcom%ndimind_f,id='hpsit_f')
          hpsit_c_tmp = f_malloc(tmb%ham_descr%collcom%ndimind_c,id='hpsit_c_tmp')
          hpsit_f_tmp = f_malloc(7*tmb%ham_descr%collcom%ndimind_f,id='hpsit_f_tmp')
          hpsi_tmp = f_malloc(tmb%ham_descr%npsidim_orbs,id='hpsi_tmp')
          psidiff = f_malloc(tmb%npsidim_orbs,id='psidiff')
          !hpsi_noprecond = f_malloc(tmb%npsidim_orbs,id='hpsi_noprecond')
    
    
          !!allocate(precond_convol_workarrays(tmb%orbs%norbp))
          !!allocate(precond_workarrays(tmb%orbs%norbp))
          !!do iorb=1,tmb%orbs%norbp
          !!    iiorb=tmb%orbs%isorb+iorb
          !!    ilr=tmb%orbs%inwhichlocreg(iiorb)
          !!    with_confpot = (tmb%confdatarr(iorb)%prefac/=0.d0)
          !!    call init_local_work_arrays(tmb%lzd%llr(ilr)%d%n1, tmb%lzd%llr(ilr)%d%n2, tmb%lzd%llr(ilr)%d%n3, &
          !!         tmb%lzd%llr(ilr)%d%nfl1, tmb%lzd%llr(ilr)%d%nfu1, &
          !!         tmb%lzd%llr(ilr)%d%nfl2, tmb%lzd%llr(ilr)%d%nfu2, &
          !!         tmb%lzd%llr(ilr)%d%nfl3, tmb%lzd%llr(ilr)%d%nfu3, &
          !!         with_confpot, precond_convol_workarrays(iorb))
          !!    kx=tmb%orbs%kpts(1,tmb%orbs%iokpt(iorb))
          !!    ky=tmb%orbs%kpts(2,tmb%orbs%iokpt(iorb))
          !!    kz=tmb%orbs%kpts(3,tmb%orbs%iokpt(iorb))
          !!    if (kx**2+ky**2+kz**2 > 0.0_gp .or. tmb%orbs%nspinor==2 ) then
          !!       ncplx=2
          !!    else
          !!       ncplx=1
          !!    end if
          !!    call allocate_work_arrays(tmb%lzd%llr(ilr)%geocode, tmb%lzd%llr(ilr)%hybrid_on, &
          !!         ncplx, tmb%lzd%llr(ilr)%d, precond_workarrays(iorb))
          !!end do
    
    
        end subroutine allocateLocalArrays
    
    
        subroutine deallocateLocalArrays()
        !
        ! Purpose:
        ! ========
        !   This subroutine deallocates all local arrays.
        !
        integer :: iiorb, ilr, ncplx
        real(gp) :: kx, ky, kz
    
        call f_free(alpha)
        call f_free(alphaDIIS)
        call f_free(fnrmOldArr)
        call f_free_ptr(hpsi_small)
        call f_free_ptr(lhphiold)
        call f_free_ptr(lphiold)
        call f_free_ptr(hpsit_c)
        call f_free_ptr(hpsit_f)
        call f_free(hpsit_c_tmp)
        call f_free(hpsit_f_tmp)
        call f_free(hpsi_tmp)
        call f_free(psidiff)
        !call f_free(hpsi_noprecond)
        !!do iorb=1,tmb%orbs%norbp
        !!    iiorb=tmb%orbs%isorb+iorb
        !!    ilr=tmb%orbs%inwhichlocreg(iiorb)
        !!    call deallocate_workarrays_quartic_convolutions(precond_convol_workarrays(iorb))
        !!    kx=tmb%orbs%kpts(1,tmb%orbs%iokpt(iorb))
        !!    ky=tmb%orbs%kpts(2,tmb%orbs%iokpt(iorb))
        !!    kz=tmb%orbs%kpts(3,tmb%orbs%iokpt(iorb))
        !!    if (kx**2+ky**2+kz**2 > 0.0_gp .or. tmb%orbs%nspinor==2 ) then
        !!       ncplx=2
        !!    else
        !!       ncplx=1
        !!    end if
        !!    call deallocate_work_arrays(tmb%lzd%llr(ilr)%geocode, tmb%lzd%llr(ilr)%hybrid_on, &
        !!         ncplx, precond_workarrays(iorb))
        !!end do
        !!deallocate(precond_convol_workarrays)
        !!deallocate(precond_workarrays)
    
        end subroutine deallocateLocalArrays
    
    
    end subroutine getLocalizedBasis



    subroutine improveOrbitals(iproc, nproc, tmb, nspin, ldiis, alpha, gradient, experimental_mode)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nspin
      type(DFT_wavefunction),intent(inout) :: tmb
      type(localizedDIISParameters),intent(inout) :: ldiis
      real(kind=8),dimension(tmb%orbs%norbp),intent(in) :: alpha
      real(kind=wp),dimension(tmb%npsidim_orbs),intent(inout) :: gradient
      logical,intent(in) :: experimental_mode
      
      ! Local variables
      integer :: istart, iorb, iiorb, ilr, ncount
    
      call f_routine(id='improveOrbitals')
    
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
              if (tmb%orbs%norbp>0) call dscal(tmb%npsidim_orbs, ldiis%alphaDIIS, gradient, 1)
          end if
          call optimizeDIIS(iproc, nproc, max(tmb%npsidim_orbs,tmb%npsidim_comp), &
               tmb%orbs, nspin, tmb%lzd, gradient, tmb%psi, ldiis, &
               experimental_mode)
      end if
    
      call f_release_routine()
    
    end subroutine improveOrbitals



    subroutine DIISorSD(iproc, it, trH, tmbopt, ldiis, alpha, alphaDIIS, lphioldopt, trH_ref, kernel_best, complete_reset)
      use module_base
      use module_types
      use yaml_output
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, it
      real(kind=8),intent(in) :: trH
      type(DFT_wavefunction),intent(inout) :: tmbopt
      type(localizedDIISParameters),intent(inout) :: ldiis
      real(kind=8),dimension(tmbopt%orbs%norbp),intent(inout) :: alpha, alphaDIIS
      real(kind=8),dimension(max(tmbopt%npsidim_orbs,tmbopt%npsidim_comp)),intent(out):: lphioldopt
      real(kind=8),intent(out) :: trH_ref
      real(kind=8),dimension(tmbopt%linmat%smat(3)%nvctrp_tg*tmbopt%linmat%smat(3)%nspin),intent(inout) :: kernel_best
      logical,intent(out) :: complete_reset
      
      ! Local variables
      integer :: idsx, ii, offset, istdest, iorb, iiorb, ilr, ncount, istsource
      character(len=2) :: numfail_char
      character(len=10) :: stepsize_char
      
    
      ! Purpose:
      ! ========
      !   This subroutine decides whether one should use DIIS or variable step size
      !   steepest descent to improve the orbitals. In the beginning we start with DIIS
      !   with history length lin%DIISHistMax. If DIIS becomes unstable, we switch to
      !   steepest descent. If the steepest descent iterations are successful, we switch
      !   back to DIIS, but decrease the DIIS history length by one. However the DIIS
      !   history length is limited to be larger or equal than lin%DIISHistMin.
    
      ! indicates whether both the support functions and the kernel have been reset
      complete_reset=.false.
    
      ! history of the energy
      if (ldiis%isx>0) then
          ldiis%energy_hist(ldiis%mis)=trH
      end if
      !!write(*,'(a,10es14.6)') 'ldiis%energy_hist', ldiis%energy_hist
    
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
          trH_ref=trH
          call vcopy(tmbopt%linmat%smat(3)%nvctrp_tg*tmbopt%linmat%smat(3)%nspin, &
               tmbopt%linmat%kernel_%matrix_compr(1), 1, kernel_best(1), 1)
          !if(iproc==0) write(*,*) 'everything ok, copy last psi...'
          call vcopy(size(tmbopt%psi), tmbopt%psi(1), 1, lphioldopt(1), 1)
    
          ! If we are using SD (i.e. diisLIN%idsx==0) and the trace has been decreasing
          ! for at least 10 iterations, switch to DIIS. However the history length is decreased.
          if(ldiis%icountSDSatur>=10 .and. ldiis%isx==0 .or. ldiis%immediateSwitchToSD) then
              ldiis%icountSwitch=ldiis%icountSwitch+1
              idsx=max(ldiis%DIISHistMin,ldiis%DIISHistMax-ldiis%icountSwitch)
              if(idsx>0) then
                  if(iproc==0) call yaml_map('Switch to DIIS with new history length',idsx)
                  !write(*,'(1x,a,i0)') 'switch to DIIS with new history length ', idsx
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
          if (ldiis%isx>0) then
              ldiis%icountDIISFailureCons=ldiis%icountDIISFailureCons+1
              ldiis%icountDIISFailureTot=ldiis%icountDIISFailureTot+1
          end if
          ldiis%icountSDSatur=0
          if((ldiis%icountDIISFailureCons>=4 .or. ldiis%icountDIISFailureTot>=6 .or. ldiis%resetDIIS) .and. ldiis%isx>0) then
              ! Switch back to SD.
              alpha=ldiis%alphaSD
              if(iproc==0) then
                  !if(ldiis%icountDIISFailureCons>=4) write(*,'(1x,a,i0,a,es10.3)') 'DIIS failed ', &
                  !    ldiis%icountDIISFailureCons, ' times consecutively. Switch to SD with stepsize', alpha(1)
                  write(numfail_char,'(i2.2)') ldiis%icountDIISFailureCons
                  write(stepsize_char,'(es10.3)') alpha(1)
                  if(ldiis%icountDIISFailureCons>=4) then
                      call yaml_warning('DIIS failed '//numfail_char//' times consecutively. &
                           &Switch to SD with stepsize'//stepsize_char//'.')
                      call yaml_newline()
                      !!write(*,'(1x,a,i0,a,es10.3)') 'DIIS failed ', &
                      !!ldiis%icountDIISFailureCons, ' times consecutively. Switch to SD with stepsize', alpha(1)
                  end if
                  !!if(ldiis%icountDIISFailureTot>=6) write(*,'(1x,a,i0,a,es10.3)') 'DIIS failed ', &
                  !!    ldiis%icountDIISFailureTot, ' times in total. Switch to SD with stepsize', alpha(1)
                  if(ldiis%icountDIISFailureTot>=6) then
                      call yaml_warning('DIIS failed '//numfail_char//' times in total. &
                           &Switch to SD with stepsize'//stepsize_char//'.' )
                      call yaml_newline()
                  end if
                  if(ldiis%resetDIIS) then
                      call yaml_warning('reset DIIS due to flag')
                      call yaml_newline()
                      !write(*,'(1x,a)') 'reset DIIS due to flag'
                  end if
                  
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
                      !!if(iproc==0) write(*,'(1x,a,i0,a)')  'Recover the orbitals from iteration ', &
                      !!    ldiis%itBest, ' which are the best so far.'
                      if (iproc==0) then
                          call yaml_map('Take best TMBs from history',ldiis%itBest)
                      end if
                  end if
                  ii=modulo(ldiis%mis-(it-ldiis%itBest)-1,ldiis%isx)+1
                  !if (iproc==0) write(*,*) 'ii',ii
                  offset=0
                  istdest=1
                  !if(iproc==0) write(*,*) 'copy DIIS history psi...'
                  do iorb=1,tmbopt%orbs%norbp
                      iiorb=tmbopt%orbs%isorb+iorb
                      ilr=tmbopt%orbs%inWhichLocreg(iiorb)
                      ncount=tmbopt%lzd%llr(ilr)%wfd%nvctr_c+7*tmbopt%lzd%llr(ilr)%wfd%nvctr_f
                      istsource=offset+(ii-1)*ncount+1
                      call vcopy(ncount, ldiis%phiHist(istsource), 1, tmbopt%psi(istdest), 1)
                      call vcopy(ncount, ldiis%phiHist(istsource), 1, lphioldopt(istdest), 1)
                      !if (iproc==0 .and. iorb==1) write(*,*) 'istsource, istdest, val', istsource, istdest, tmbopt%psi(istdest)
                      offset=offset+ldiis%isx*ncount
                      istdest=istdest+ncount
                  end do
                  trH_ref=ldiis%energy_hist(ii)
                  !!if (iproc==0) write(*,*) 'take energy from entry',ii
                  call vcopy(tmbopt%linmat%smat(3)%nvctrp_tg*tmbopt%linmat%smat(3)%nspin, &
                       kernel_best(1), 1, tmbopt%linmat%kernel_%matrix_compr(1), 1)
                  !!call vcopy(tmbopt%linmat%smat(3)%nvctr, kernel_best(1), 1, tmbopt%linmat%denskern_large%matrix_compr(1), 1)
                  complete_reset=.true.
              else
                  !if(iproc==0) write(*,*) 'copy last psi...'
                  call vcopy(size(tmbopt%psi), tmbopt%psi(1), 1, lphioldopt(1), 1)
                  trH_ref=trH
              end if
              ldiis%isx=0
              ldiis%switchSD=.true.
          end if
          ! to indicate that no orthonormalization is required... (CHECK THIS!)
          !if(ldiis%isx==0) ldiis%switchSD=.true. 
      end if
    
    end subroutine DIISorSD



    !> Estimate the energy change, given by the product of the force and the "displacement" .
    subroutine estimate_energy_change(npsidim_orbs, orbs, lzd, nspin, psidiff, hpsi_noprecond, delta_energy)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer, intent(in) :: npsidim_orbs, nspin
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      real(kind=8),dimension(npsidim_orbs),intent(in) :: psidiff, hpsi_noprecond
      real(kind=8),intent(out) :: delta_energy
    
      ! Local variables
      integer :: ist, iorb, iiorb, ilr, ncount
      real(kind=8) :: tt, ddot
    
      call f_routine(id='estimate_energy_change')
    
      ist=1
      delta_energy=0.d0
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
          tt=ddot(ncount, psidiff(ist), 1, hpsi_noprecond(ist), 1)
          delta_energy=delta_energy+2.0d0*tt
          ist=ist+ncount
      end do
    
      if (nspin==1) then
          delta_energy = 2.d0*delta_energy
      end if
    
      if (bigdft_mpi%nproc > 1) then
          call mpiallred(delta_energy, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
    
      call f_release_routine()
    
    end subroutine estimate_energy_change



    subroutine allocate_precond_arrays(orbs, lzd, confdatarr, precond_convol_workarrays, precond_workarrays)
      use module_base, only: gp
      use module_types
      use locreg_operations
      implicit none
      ! Calling arguments
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      type(confpot_data),dimension(orbs%norbp),intent(in) ::  confdatarr
      type(workarrays_quartic_convolutions),dimension(:),pointer,intent(inout) :: precond_convol_workarrays
      type(workarr_precond),dimension(:),pointer,intent(inout) :: precond_workarrays
    
      ! Local variables
      integer :: iorb, iiorb, ilr, ncplx
      real(kind=8) :: kx, ky, kz
      logical :: with_confpot
    
      allocate(precond_convol_workarrays(orbs%norbp))
      allocate(precond_workarrays(orbs%norbp))
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          with_confpot = (confdatarr(iorb)%prefac/=0.d0)
          call init_local_work_arrays(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
               lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
               lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
               lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
               with_confpot, precond_convol_workarrays(iorb))
          kx=orbs%kpts(1,orbs%iokpt(iorb))
          ky=orbs%kpts(2,orbs%iokpt(iorb))
          kz=orbs%kpts(3,orbs%iokpt(iorb))
          if (kx**2+ky**2+kz**2 > 0.0_gp .or. orbs%nspinor==2 ) then
             ncplx=2
          else
             ncplx=1
          end if
          call allocate_work_arrays(lzd%llr(ilr)%geocode, lzd%llr(ilr)%hybrid_on, &
               ncplx, lzd%llr(ilr)%d, precond_workarrays(iorb))
      end do
    
    end subroutine allocate_precond_arrays
    
    
    subroutine deallocate_precond_arrays(orbs, lzd, precond_convol_workarrays, precond_workarrays)
      use module_base, only: gp
      use module_types
      use locreg_operations
      implicit none
      ! Calling arguments
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      type(workarrays_quartic_convolutions),dimension(:),pointer,intent(inout) :: precond_convol_workarrays
      type(workarr_precond),dimension(:),pointer,intent(inout) :: precond_workarrays
    
      ! Local variables
      integer :: iorb, iiorb, ilr, ncplx
      real(kind=8) :: kx, ky, kz
    
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          call deallocate_workarrays_quartic_convolutions(precond_convol_workarrays(iorb))
          kx=orbs%kpts(1,orbs%iokpt(iorb))
          ky=orbs%kpts(2,orbs%iokpt(iorb))
          kz=orbs%kpts(3,orbs%iokpt(iorb))
          if (kx**2+ky**2+kz**2 > 0.0_gp .or. orbs%nspinor==2 ) then
             ncplx=2
          else
             ncplx=1
          end if
          call deallocate_work_arrays(lzd%llr(ilr)%geocode, lzd%llr(ilr)%hybrid_on, &
               ncplx, precond_workarrays(iorb))
      end do
      deallocate(precond_convol_workarrays)
      deallocate(precond_workarrays)
    
    end subroutine deallocate_precond_arrays


    subroutine calculate_energy_and_gradient_linear(iproc, nproc, it, &
               ldiis, fnrmOldArr, fnrm_old, alpha, trH, trHold, fnrm, alpha_mean, alpha_max, &
               energy_increased, tmb, lhphiold, overlap_calculated, &
               energs, hpsit_c, hpsit_f, nit_precond, target_function, correction_orthoconstraint, &
               hpsi_small, experimental_mode, calculate_inverse, correction_co_contra, recovered_old_kernel, &
               hpsi_noprecond, norder_taylor, max_inversion_error, precond_convol_workarrays, precond_workarrays,&
               wt_hphi, wt_philarge, &
               cdft, input_frag, ref_frags)
      use module_base
      use module_types
      use yaml_output
      use module_interfaces, only: preconditionall2
      use communications_base, only: work_transpose, TRANSPOSE_FULL, TRANSPOSE_GATHER
      use communications, only: transpose_localized, untranspose_localized
      use sparsematrix_base, only: matrices, matrices_null, deallocate_matrices, &
                                   sparsematrix_malloc_ptr, assignment(=), &
                                   sparsematrix_malloc, SPARSE_TASKGROUP
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: transform_sparse_matrix
      use constrained_dft, only: cdft_data
      use module_fragments, only: system_fragment,fragmentInputParameters
      use transposed_operations, only: calculate_overlap_transposed, build_linear_combination_transposed
      use public_enums
      use orthonormalization, only: orthoconstraintNonorthogonal
      use locreg_operations
      use coeffs, only: calculate_kernel_and_energy
      implicit none
    
      ! Calling arguments
      integer, intent(in) :: iproc, nproc, it
      integer,intent(inout) :: norder_taylor
      real(kind=8),intent(in) :: max_inversion_error
      type(DFT_wavefunction), target, intent(inout):: tmb
      type(localizedDIISParameters), intent(inout) :: ldiis
      real(kind=8), dimension(tmb%orbs%norbp), intent(inout) :: fnrmOldArr
      real(kind=8),intent(inout) :: fnrm_old
      real(kind=8), dimension(tmb%orbs%norbp), intent(inout) :: alpha
      real(kind=8), intent(out):: trH, alpha_mean, alpha_max
      type(work_mpiaccumulate), intent(inout):: fnrm
      real(kind=8), intent(in):: trHold
      logical,intent(out) :: energy_increased
      real(kind=8), dimension(tmb%npsidim_orbs), intent(inout):: lhphiold
      logical, intent(inout):: overlap_calculated
      type(energy_terms), intent(in) :: energs
      real(kind=8),dimension(tmb%ham_descr%collcom%ndimind_c) :: hpsit_c
      real(kind=8),dimension(7*tmb%ham_descr%collcom%ndimind_f) :: hpsit_f
      integer, intent(in) :: nit_precond, target_function, correction_orthoconstraint
      logical, intent(in) :: experimental_mode, calculate_inverse, correction_co_contra, recovered_old_kernel
      real(kind=8), dimension(tmb%npsidim_orbs), intent(out) :: hpsi_small
      real(kind=8), dimension(tmb%npsidim_orbs),intent(out) :: hpsi_noprecond
      type(workarrays_quartic_convolutions),dimension(tmb%orbs%norbp),intent(inout) :: precond_convol_workarrays
      type(workarr_precond),dimension(tmb%orbs%norbp),intent(inout) :: precond_workarrays
      type(work_transpose),intent(inout) :: wt_hphi
      type(work_transpose),intent(inout) :: wt_philarge
      !!!these must all be present together
      type(cdft_data),intent(inout),optional :: cdft
      type(fragmentInputParameters),optional,intent(in) :: input_frag
      type(system_fragment), dimension(:), optional, intent(in) :: ref_frags
    
    
      ! Local variables
      integer :: iorb, iiorb, ilr, ncount, ierr, ist, ncnt, istat, iall, ii, jjorb, i
      integer :: lwork, info, ishift,ispin, iseg, request
      real(kind=8) :: ddot, tt, fnrmOvrlp_tot, fnrm_tot, fnrmold_tot, tt2, trkw, trH_sendbuf
      real(kind=8), dimension(:), pointer :: hpsittmp_c, hpsittmp_f
      real(kind=8), dimension(:), allocatable :: hpsi_conf, hpsit_c_orig, hpsit_f_orig
      real(kind=8), dimension(:), pointer :: kernel_compr_tmp
      real(kind=8), dimension(:), allocatable :: prefac, tmparr
      real(kind=8),dimension(2) :: reducearr
      integer,dimension(2) :: irowcol
      real(wp), dimension(2) :: garray
      real(dp) :: gnrm,gnrm_zero,gnrmMax,gnrm_old ! for preconditional2, replace with fnrm eventually, but keep separate for now
      type(matrices) :: matrixm
      real(kind=8),dimension(:),pointer :: cdft_gradt_c, cdft_gradt_f, cdft_grad, cdft_grad_small
      !type(work_transpose) :: wt_hphi
    
      call f_routine(id='calculate_energy_and_gradient_linear')
    
      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          hpsi_conf = f_malloc(tmb%npsidim_orbs,id='hpsi_conf')
          call large_to_small_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
               tmb%orbs, tmb%hpsi, hpsi_conf)
          call timing(iproc,'buildgrad_mcpy','ON')
          ist=1
          do iorb=1,tmb%orbs%norbp
              iiorb=tmb%orbs%isorb+iorb
              ilr=tmb%orbs%inwhichlocreg(iiorb)
              ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
              tt=ddot(ncount, hpsi_conf(ist), 1, tmb%psi(ist), 1)
              call daxpy(ncount, -tt, tmb%psi(ist), 1, hpsi_conf(ist), 1)
              ist=ist+ncount
          end do
          call timing(iproc,'buildgrad_mcpy','OF')
      end if
    
      ! by default no quick exit
      energy_increased=.false.
    
      if (target_function==TARGET_FUNCTION_IS_ENERGY .and. iproc==0) then
          !!do i=1,size(hpsit_c)
          !!    write(4231,'(a,i8,2es14.6)') 'i, tmb%ham_descr%psit_c(i), hpsit_c(i)', i, tmb%ham_descr%psit_c(i), hpsit_c(i)
          !!end do
          !!do i=1,size(hpsit_f)
          !!    write(4232,'(a,i8,2es14.6)') 'i, tmb%ham_descr%psit_f(i), hpsit_f(i)', i, tmb%ham_descr%psit_f(i), hpsit_f(i)
          !!end do
      end if
    
      hpsittmp_c = f_malloc_ptr(tmb%ham_descr%collcom%ndimind_c,id='hpsittmp_c')
      hpsittmp_f = f_malloc_ptr(7*tmb%ham_descr%collcom%ndimind_f,id='hpsittmp_f')
    
      if(target_function==TARGET_FUNCTION_IS_ENERGY .or. &
         target_function==TARGET_FUNCTION_IS_HYBRID) then
         !write(*,*) 'sum(K)',sum(tmb%linmat%kernel_%matrix_compr)
         call build_gradient(iproc, nproc, tmb, target_function, hpsit_c, hpsit_f, hpsittmp_c, hpsittmp_f)
      end if


      !@NEW Calculate Omega in a different way ####################################
      hpsit_c_orig = f_malloc(tmb%ham_descr%collcom%ndimind_c,id='hpsit_c_orig')
      hpsit_f_orig = f_malloc(7*tmb%ham_descr%collcom%ndimind_f,id='hpsit_f_orig')
      call f_memcpy(src=hpsit_c, dest=hpsit_c_orig)
      call f_memcpy(src=hpsit_f, dest=hpsit_f_orig)
      !############################################################################
    
     
      ! WARNING: TO BE CHECKED!!!!
      ! For the non polarized case, a factor of two is already included in the
      ! kernel. Therefore explicitely add this factor for the polarized case 
      ! in order to make the two cases analogous.
      if (tmb%linmat%smat(3)%nspin==2 .and. &
          (target_function==TARGET_FUNCTION_IS_ENERGY .or. target_function==TARGET_FUNCTION_IS_HYBRID)) then
          if (iproc==0) call yaml_warning('multiply the gradient by 2.0, check this!')
          !hpsit_c=2.d0*hpsit_c
          !hpsit_f=2.d0*hpsit_f
          call vscal(tmb%ham_descr%collcom%ndimind_c, 2.d0, hpsit_c(1), 1)
          call vscal(7*tmb%ham_descr%collcom%ndimind_f, 2.d0, hpsit_f(1), 1)
      end if
    
    
      if (correction_co_contra) then
          !@NEW correction for contra / covariant gradient
    
          if (target_function/=TARGET_FUNCTION_IS_HYBRID) then
              call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
                   TRANSPOSE_FULL, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
              tmb%can_use_transposed=.true.
    
              call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
                   tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%smat(1), tmb%linmat%auxs, tmb%linmat%ovrlp_)
              !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(1), tmb%linmat%ovrlp_)
          end if
          call vcopy(tmb%ham_descr%collcom%ndimind_c, hpsit_c(1), 1, hpsittmp_c(1), 1)
          call vcopy(7*tmb%ham_descr%collcom%ndimind_f, hpsit_f(1), 1, hpsittmp_f(1), 1)
    
          ! Transform to the larger sparse region in order to be compatible with tmb%ham_descr%collcom.
          ! To this end use ham_.
          !if (iproc==0) write(*,*) 'sum(S)',sum(tmb%linmat%ovrlp_%matrix_compr)
          call transform_sparse_matrix(iproc, tmb%linmat%smat(1), tmb%linmat%smat(2), SPARSE_TASKGROUP, 'small_to_large', &
               smat_in=tmb%linmat%ovrlp_%matrix_compr, lmat_out=tmb%linmat%ham_%matrix_compr)
    
          !tmparr = sparsematrix_malloc(tmb%linmat%smat(2),iaction=SPARSE_FULL,id='tmparr')
          !call vcopy(tmb%linmat%smat(2)%nvctr, tmb%linmat%ham_%matrix_compr(1), 1, tmparr(1), 1)
          !call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(2), tmb%linmat%ham_)
          call build_linear_combination_transposed(tmb%ham_descr%collcom, &
               tmb%linmat%smat(2), tmb%linmat%auxm, tmb%linmat%ham_, hpsittmp_c, hpsittmp_f, .true., hpsit_c, hpsit_f, iproc)
          !call vcopy(tmb%linmat%smat(2)%nvctr, tmparr(1), 1, tmb%linmat%ham_%matrix_compr(1), 1)
          !call f_free(tmparr)
    
    
          !@END NEW correction for contra / covariant gradient
      end if
    
      call f_free_ptr(hpsittmp_c)
      call f_free_ptr(hpsittmp_f)
    
    
      ! Calculate the overlap matrix if necessary
      if (.not.correction_co_contra) then
          call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
               TRANSPOSE_FULL, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
          tmb%can_use_transposed=.true.
          call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
               tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%smat(1), tmb%linmat%auxs, tmb%linmat%ovrlp_)
          !call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(1), tmb%linmat%ovrlp_)
      end if
    
    
      !!if (target_function==TARGET_FUNCTION_IS_ENERGY .and. iproc==0) then
      !!    ist=0
      !!    do iorb=1,tmb%orbs%norbp
      !!        iiorb=tmb%orbs%isorb+iorb
      !!        ilr=tmb%orbs%inwhichlocreg(iiorb)
      !!        ncount=tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_f
      !!        do i=1,ncount
      !!            ist=ist+1
      !!            if (tmb%orbs%spinsgn(iiorb)>0.d0) then
      !!                write(4211,'(a,2i10,f8.1,2es16.7)') 'iiorb, ist, spin, vals', iiorb, ist, tmb%orbs%spinsgn(iiorb), tmb%ham_descr%psi(ist), tmb%hpsi(ist)
      !!            else
      !!                write(4212,'(a,2i10,f8.1,2es16.7)') 'iiorb, ist, spin, val', iiorb, ist, tmb%orbs%spinsgn(iiorb), tmb%ham_descr%psi(ist), tmb%hpsi(ist)
      !!            end if
      !!        end do
      !!    end do
      !!    do i=1,size(hpsit_c)
      !!        write(4221,'(a,i8,2es14.6)') 'i, tmb%ham_descr%psit_c(i), hpsit_c(i)', i, tmb%ham_descr%psit_c(i), hpsit_c(i)
      !!    end do
      !!    do i=1,size(hpsit_f)
      !!        write(4222,'(a,i8,2es14.6)') 'i, tmb%ham_descr%psit_f(i), hpsit_f(i)', i, tmb%ham_descr%psit_f(i), hpsit_f(i)
      !!    end do
      !!end if

      !!if (iproc==0) write(*,*) 'before orthoconstr: sum(tmb%hpsi)',sum(tmb%hpsi)
    
      call orthoconstraintNonorthogonal(iproc, nproc, tmb%ham_descr%lzd, &
           tmb%ham_descr%npsidim_orbs, tmb%ham_descr%npsidim_comp, &
           tmb%orbs, tmb%ham_descr%collcom, tmb%orthpar, tmb%ice_obj, correction_orthoconstraint, &
           tmb%linmat, tmb%ham_descr%psi, tmb%hpsi, &
           tmb%linmat%smat(2), tmb%linmat%auxm, tmb%linmat%ham_, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, &
           hpsit_c, hpsit_f, hpsit_c_orig, hpsit_f_orig, &
           tmb%ham_descr%can_use_transposed, &
           overlap_calculated, calculate_inverse, norder_taylor, max_inversion_error, &
           tmb%npsidim_orbs, tmb%lzd, hpsi_noprecond, wt_philarge, wt_hphi)

      call f_free(hpsit_c_orig)
      call f_free(hpsit_f_orig)
    
    
      ! Calculate trace (or band structure energy, resp.)
      ! Here the calculation of the trace is only initiated, trH is not yet available
      call calculate_trace_start()
      call untranspose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
           TRANSPOSE_GATHER, hpsit_c, hpsit_f, tmb%hpsi, tmb%ham_descr%lzd, wt_hphi)

      !!if (iproc==0) write(*,*) 'after orthoconstr: sum(tmb%hpsi)',sum(tmb%hpsi)
    
    
      !EXPERIMENTAL and therefore deactivated
      ! SM: WARNING trH is only available after the call to calculate_trace_finish
      !add CDFT gradient, or at least an approximation thereof
      if (present(cdft).and..false.) then
         if (.not.present(input_frag).or..not.present(ref_frags)) stop 'input_frag, ref_frags and cdft must be present together'
         cdft_gradt_c = f_malloc_ptr(tmb%ham_descr%collcom%ndimind_c,id='cdft_gradt_c')
         cdft_gradt_f = f_malloc_ptr(7*tmb%ham_descr%collcom%ndimind_f,id='cdft_gradt_f')
         !calculate gradient (1st order taylor), assume S needs calculating though might not be needed
         !print*,size(cdft_gradt_c),size(cdft_gradt_f),tmb%ham_descr%collcom%ndimind_c,7*tmb%ham_descr%collcom%ndimind_f
    
    
         call calculate_weight_matrix_lowdin_gradient(cdft%weight_matrix,cdft%weight_matrix_,cdft%ifrag_charged,&
              tmb,input_frag,ref_frags,.true.,.true.,norder_taylor,cdft_gradt_c,cdft_gradt_f)
         !add gradient to hpsi_t
         !print*,'corr',cdft%lag_mult**2*ddot(tmb%ham_descr%collcom%ndimind_c, cdft_gradt_c(1), 1, cdft_gradt_c(1), 1),&
         !     cdft%lag_mult**2*ddot(7*tmb%ham_descr%collcom%ndimind_f, cdft_gradt_f(1), 1, cdft_gradt_f(1), 1)
         !print*,'orig',ddot(tmb%ham_descr%collcom%ndimind_c, hpsit_c(1), 1, hpsit_c(1), 1),&
         !          ddot(7*tmb%ham_descr%collcom%ndimind_f, hpsit_f(1), 1, hpsit_f(1), 1)
         !call daxpy(tmb%ham_descr%collcom%ndimind_c,cdft%lag_mult,cdft_gradt_c,1,hpsit_c,1)
         !call daxpy(7*tmb%ham_descr%collcom%ndimind_f,cdft%lag_mult,cdft_gradt_f,1,hpsit_f,1)
         !print*,'after',ddot(tmb%ham_descr%collcom%ndimind_c, hpsit_c(1), 1, hpsit_c(1), 1),&
         !        ddot(7*tmb%ham_descr%collcom%ndimind_f, hpsit_f(1), 1, hpsit_f(1), 1)
         !call dcopy(tmb%ham_descr%collcom%ndimind_c,cdft_gradt_c,1,hpsit_c,1)
         !call dcopy(7*tmb%ham_descr%collcom%ndimind_f,cdft_gradt_f,1,hpsit_f,1)
    
         cdft_grad=f_malloc_ptr(tmb%ham_descr%npsidim_orbs,id='cdft_grad')
         call untranspose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
              TRANSPOSE_FULL, cdft_gradt_c, cdft_gradt_f, cdft_grad, tmb%ham_descr%lzd)
         !print*,ddot(tmb%ham_descr%npsidim_orbs, cdft_grad(1), 1, cdft_grad(1), 1)
    
         if (.false.) then
            cdft_grad_small=f_malloc_ptr(tmb%npsidim_orbs,id='cdft_grad_small')
            !no point keeping in tmblarge for now as fd will only be in small
            call large_to_small_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
                 tmb%orbs, cdft_grad, cdft_grad_small)
            !print CDFT gradient
            !open(10+iproc)
            !do iorb=19000,21500 !1,tmb%npsidim_orbs
            !write(10+iproc,*) cdft_grad_small(iorb)
            !end do
            !close(10+iproc)
    
            call calculate_weight_matrix_lowdin_gradient_fd(cdft%weight_matrix,cdft%weight_matrix_,cdft%ifrag_charged,&
                 tmb,input_frag,ref_frags,.true.,.true.,norder_taylor,cdft_grad_small)
            !call untranspose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
            !     cdft_gradt_c, cdft_gradt_f, cdft_grad, tmb%ham_descr%lzd)
    
            !open(20+iproc)
            !do iorb=19000,21500 !1,tmb%npsidim_orbs
            !write(20+iproc,*) cdft_grad_small(iorb)
            !end do
            !close(20+iproc)
    
            !call f_free_ptr(cdft_grad_small)
    
            !call mpi_finalize(bigdft_mpi%mpi_comm)
            !stop
         end if
    
         call f_free_ptr(cdft_gradt_c)
         call f_free_ptr(cdft_gradt_f)   
         !call daxpy(tmb%ham_descr%npsidim_orbs,cdft%lag_mult,cdft_grad,1,tmb%hpsi,1)
         call dcopy(tmb%ham_descr%npsidim_orbs,cdft_grad,1,tmb%hpsi,1)
         !call dscal(tmb%ham_descr%npsidim_orbs,cdft%lag_mult,tmb%hpsi)
    
         call f_free_ptr(cdft_grad)
      end if
    
    
      !!if (target_function==TARGET_FUNCTION_IS_ENERGY .and. iproc==0) then
      !!    ist=0
      !!    do iorb=1,tmb%orbs%norbp
      !!        iiorb=tmb%orbs%isorb+iorb
      !!        ilr=tmb%orbs%inwhichlocreg(iiorb)
      !!        ncount=tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_f
      !!        do i=1,ncount
      !!            ist=ist+1
      !!            if (tmb%orbs%spinsgn(iiorb)>0.d0) then
      !!                write(4201,'(a,2i10,f8.1,es16.7)') 'iiorb, ist, spin, vals', iiorb, ist, tmb%orbs%spinsgn(iiorb), tmb%hpsi(ist)
      !!            else
      !!                write(4202,'(a,2i10,f8.1,es16.7)') 'iiorb, ist, spin, val', iiorb, ist, tmb%orbs%spinsgn(iiorb), tmb%hpsi(ist)
      !!            end if
      !!        end do
      !!    end do
      !!end if
    
    
    
      !temporary debug
      !call daxpy(tmb%npsidim_orbs,cdft%lag_mult,cdft_grad_small,1,hpsi_small,1)
      !call f_free_ptr(cdft_grad_small)
    
      !!if (target_function==TARGET_FUNCTION_IS_ENERGY .and. iproc==0) then
      !!    ist=0
      !!    do iorb=1,tmb%orbs%norbp
      !!        iiorb=tmb%orbs%isorb+iorb
      !!        ilr=tmb%orbs%inwhichlocreg(iiorb)
      !!        ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      !!        do i=1,ncount
      !!            ist=ist+1
      !!            if (tmb%orbs%spinsgn(iiorb)>0.d0) then
      !!                write(4301,'(a,2i10,f8.1,es16.7)') 'iiorb, ist, spin, vals', iiorb, ist, tmb%orbs%spinsgn(iiorb), hpsi_small(ist)
      !!            else
      !!                write(4302,'(a,2i10,f8.1,es16.7)') 'iiorb, ist, spin, val', iiorb, ist, tmb%orbs%spinsgn(iiorb), hpsi_small(ist)
      !!            end if
      !!        end do
      !!    end do
      !!end if
    
      !experimental
      if (present(cdft).and..false.) then
         !only correct energy not gradient for now
         !can give tmb%orbs twice as ksorbs is only used for recalculating the kernel
         stop 'MAKE SURE THAN calculate_kernel_and_energy IS CALLED APPRORIATELY:'
         call calculate_kernel_and_energy(iproc,nproc,bigdft_mpi%mpi_comm,tmb%linmat%smat(3),cdft%weight_matrix, &
              tmb%linmat%kernel_,cdft%weight_matrix_,trkw,tmb%coeff, &
              tmb%orbs%norbp, tmb%orbs%isorb, tmb%orbs%norbu, tmb%orbs%norb, tmb%orbs%occup, .false.)
         !cdft%charge is always constant (as is lagmult in this loop) so could in theory be ignored as in optimize_coeffs
         trH = trH + cdft%lag_mult*(trkw - cdft%charge)
         if (iproc==0) print*,'trH,trH+V(trkw-N),V(trkw-N)',trH-cdft%lag_mult*(trkw - cdft%charge),&
              trH,cdft%lag_mult*(trkw - cdft%charge)
       end if
    
    
    
    
      call large_to_small_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
           tmb%orbs, tmb%hpsi, hpsi_small)
    
      ! Copy the small gradient back to tmb%hpsi as a temporary array 
      call vcopy(tmb%npsidim_orbs, hpsi_small(1), 1, tmb%hpsi(1), 1)
    
      ! Do the preconditioning while trH is communicated
      if(target_function==TARGET_FUNCTION_IS_HYBRID) then
         ist=1
         do iorb=1,tmb%orbs%norbp
            iiorb=tmb%orbs%isorb+iorb
            ilr = tmb%orbs%inWhichLocreg(iiorb)
            ncnt=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
    
            tt=ddot(ncnt, hpsi_conf(ist), 1, hpsi_small(ist), 1)
            tt=tt/ddot(ncnt, hpsi_conf(ist), 1, hpsi_conf(ist), 1)
            call daxpy(ncnt, -tt, hpsi_conf(ist), 1, hpsi_small(ist), 1)
            call dscal(ncnt, tt, hpsi_conf(ist), 1)
    
            ist=ist+ncnt
         end do
    
         call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
              tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
              nit_precond,tmb%npsidim_orbs,hpsi_conf,tmb%confdatarr,gnrm,gnrm_zero, &
              precond_convol_workarrays, precond_workarrays)
    
         ! temporarily turn confining potential off...
         prefac = f_malloc(tmb%orbs%norbp,id='prefac')
         prefac(:)=tmb%confdatarr(:)%prefac
         tmb%confdatarr(:)%prefac=0.0d0
         call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
              tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
              nit_precond,tmb%npsidim_orbs,hpsi_small,tmb%confdatarr,gnrm,gnrm_zero, & ! prefac should be zero
              precond_convol_workarrays, precond_workarrays)
         call daxpy(tmb%npsidim_orbs, 1.d0, hpsi_conf(1), 1, hpsi_small(1), 1)
         ! ...revert back to correct value
         tmb%confdatarr(:)%prefac=prefac
    
         call f_free(prefac)
         call f_free(hpsi_conf)
      else
         call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
              tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
              nit_precond,tmb%npsidim_orbs,hpsi_small,tmb%confdatarr,gnrm,gnrm_zero,&
              precond_convol_workarrays, precond_workarrays)
      end if
    
    
      ! Here the calculation of the trace is completed, trH is available
      call calculate_trace_finish()
    
      ! WARNING: TO BE CHECKED!!!!
      ! For the polarized case, a factor of two was already multiplied to the
      ! gradient (see above). Therefore now undo this again for the band structure
      ! energy in order to get the analogous result to the non-polarized case.
      if (tmb%linmat%smat(3)%nspin==2 .and. &
          (target_function==TARGET_FUNCTION_IS_ENERGY .or. target_function==TARGET_FUNCTION_IS_HYBRID)) then
          if (iproc==0) call yaml_warning('divide the band stucture energy by 2.0, check this!')
          trH=0.5d0*trH
      end if
    
      ! trH is now the total energy (name is misleading, correct this)
      ! Multiply by 2 because when minimizing trace we don't have kernel
      if(tmb%orbs%nspin==1 .and. target_function==TARGET_FUNCTION_IS_TRACE) trH=2.d0*trH
      !if (iproc==0) call yaml_map('Omega old',trH)
      !!if (iproc==0) write(*,'(a,6es17.8)') 'eh, exc, evxc, eexctX, eion, edisp', &
      !!    energs%eh,energs%exc,energs%evxc,energs%eexctX,energs%eion,energs%edisp
      trH=trH-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
      !!if (iproc==0) write(*,*) 'trH',trH
    
    
      ! Determine whether the target function is increasing
      !if(.not. ldiis%switchSD .and. ldiis%isx==0) then
      if(.not. ldiis%switchSD .and. .not.recovered_old_kernel) then
          if(trH > ldiis%trmin+1.d-12*abs(ldiis%trmin)) then !1.d-12 is here to tolerate some noise...
              !!if(iproc==0) write(*,'(1x,a,es18.10,a,es18.10)') &
              !!    'WARNING: the target function is larger than its minimal value reached so far:',trH,' > ', ldiis%trmin
              if (iproc==0) then
                  call yaml_newline()
                  call yaml_warning('target function larger than its minimal value reached so far, &
                      &D='//trim(yaml_toa(trH-ldiis%trmin,fmt='(1es10.3)')))!//'. &
                      !&Decrease step size and restart with previous TMBs')
              end if
              !if(iproc==0) write(*,'(1x,a)') 'Decrease step size and restart with previous TMBs'
              energy_increased=.true.
          end if
      end if
    
    
      ! Determine the gradient norm and its maximal component. In addition, adapt the
      ! step size for the steepest descent minimization (depending on the angle 
      ! between the current gradient and the one from the previous iteration).
      ! This is of course only necessary if we are using steepest descent and not DIIS.
      ! if newgradient is true, the angle criterion cannot be used and the choice whether to
      ! decrease or increase the step size is only based on the fact whether the trace decreased or increased.
      fnrm%sendbuf(1)=0.d0
      fnrmOvrlp_tot=0.d0
      fnrmOld_tot=0.d0
      ist=1
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          tt = ddot(ncount, tmb%hpsi(ist), 1, tmb%hpsi(ist), 1)
          fnrm%sendbuf(1) = fnrm%sendbuf(1) + tt
          if(it>1) then
              tt2=ddot(ncount, tmb%hpsi(ist), 1, lhphiold(ist), 1)
              fnrmOvrlp_tot = fnrmOvrlp_tot + tt2
              if(ldiis%isx==0 .and. .not.ldiis%switchSD) then
                  ! Adapt step size for the steepest descent minimization.
                  if (.not.experimental_mode) then
                      tt2=tt2/sqrt(tt*fnrmOldArr(iorb))
                      ! apply thresholds so that alpha never goes below around 1.d-2 and above around 2
                      if(tt2>.6d0 .and. trH<trHold .and. alpha(iorb)<1.8d0) then
                          !write(*,*) 'incr alpha 1.1 sub'
                          alpha(iorb)=alpha(iorb)*1.1d0
                      else if (alpha(iorb)>1.7d-3) then
                          !write(*,*) 'cut alpha 0.6 sub'
                          alpha(iorb)=alpha(iorb)*.6d0
                      end if
                  end if
              end if
          end if
          fnrmOldArr(iorb)=tt
          ist=ist+ncount
      end do
      
    
      call communicate_fnrm()
    
    
      if (experimental_mode .and. it>1 .and. ldiis%isx==0 .and. .not.ldiis%switchSD) then
          !do iorb=1,tmb%orbs%norbp
          !    fnrmOld_tot=fnrmOld_tot+fnrmOldArr(iorb)
          !end do
          reducearr(1) = fnrmOvrlp_tot
          reducearr(2) = fnrm%sendbuf(1)
          if (nproc>1) then
              call mpiallred(reducearr, mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if
          !tt2=fnrmOvrlp_tot/sqrt(fnrm*fnrmOld_tot)
          tt2=reducearr(1)/sqrt(reducearr(2)*fnrm_old)
          ! apply thresholds so that alpha never goes below around 1.d-2 and above around 2
          if(tt2>.6d0 .and. trH<trHold .and. alpha(1)<1.8d0) then ! take alpha(1) since the value is the same for all
              alpha(:)=alpha(:)*1.1d0
              !write(*,*) 'incr alpha 1.1 sub'
          else if (alpha(1)>1.7d-3) then
              !write(*,*) 'cut alpha 0.6 sub'
              alpha(:)=alpha(:)*.6d0
          end if
      end if
    
      !fnrm_old=fnrm%receivebuf(1) ! This value will be used in th next call to this routine
    
      !!fnrm%sendbuf(1)=sqrt(fnrm%sendbuf(1)/dble(tmb%orbs%norb))
    
    
    
      ! Determine the mean step size for steepest descent iterations.
      call communicate_alpha()
      tt=sum(alpha)
      alpha_max=maxval(alpha)
      if (nproc > 1) then
         call mpiallred(tt, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
         call mpiallred(alpha_max, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
      end if
      alpha_mean=tt/dble(tmb%orbs%norb)
    
      ! Copy the gradient (will be used in the next iteration to adapt the step size).
      call vcopy(tmb%npsidim_orbs, tmb%hpsi(1), 1, lhphiold(1), 1)
      !!call timing(iproc,'buildgrad_mcpy','OF')
    
      ! if energy has increased or we only wanted to calculate the energy, not gradient, we can return here
      ! rather than calculating the preconditioning for nothing
      ! Do this only for steepest descent
      if ((energy_increased) .and. target_function/=TARGET_FUNCTION_IS_HYBRID .and. ldiis%isx==0) return
    
    
    
      !!@if(target_function==TARGET_FUNCTION_IS_HYBRID) then
      !!@   ist=1
      !!@   do iorb=1,tmb%orbs%norbp
      !!@      iiorb=tmb%orbs%isorb+iorb
      !!@      ilr = tmb%orbs%inWhichLocreg(iiorb)
      !!@      ncnt=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
    
      !!@      tt=ddot(ncnt, hpsi_conf(ist), 1, hpsi_small(ist), 1)
      !!@      tt=tt/ddot(ncnt, hpsi_conf(ist), 1, hpsi_conf(ist), 1)
      !!@      call daxpy(ncnt, -tt, hpsi_conf(ist), 1, hpsi_small(ist), 1)
      !!@      call dscal(ncnt, tt, hpsi_conf(ist), 1)
    
      !!@      ist=ist+ncnt
      !!@   end do
    
      !!@   call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
      !!@        tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
      !!@        nit_precond,tmb%npsidim_orbs,hpsi_conf,tmb%confdatarr,gnrm,gnrm_zero, &
      !!@        precond_convol_workarrays, precond_workarrays)
    
      !!@   ! temporarily turn confining potential off...
      !!@   prefac = f_malloc(tmb%orbs%norbp,id='prefac')
      !!@   prefac(:)=tmb%confdatarr(:)%prefac
      !!@   tmb%confdatarr(:)%prefac=0.0d0
      !!@   call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
      !!@        tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
      !!@        nit_precond,tmb%npsidim_orbs,hpsi_small,tmb%confdatarr,gnrm,gnrm_zero, & ! prefac should be zero
      !!@        precond_convol_workarrays, precond_workarrays)
      !!@   call daxpy(tmb%npsidim_orbs, 1.d0, hpsi_conf(1), 1, hpsi_small(1), 1)
      !!@   ! ...revert back to correct value
      !!@   tmb%confdatarr(:)%prefac=prefac
    
      !!@   call f_free(prefac)
      !!@   call f_free(hpsi_conf)
      !!@else
      !!@   call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
      !!@        tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
      !!@        nit_precond,tmb%npsidim_orbs,hpsi_small,tmb%confdatarr,gnrm,gnrm_zero,&
      !!@        precond_convol_workarrays, precond_workarrays)
      !!@end if
    
      !!if (target_function==TARGET_FUNCTION_IS_ENERGY .and. iproc==0) then
      !!    ist=0
      !!    do iorb=1,tmb%orbs%norbp
      !!        iiorb=tmb%orbs%isorb+iorb
      !!        ilr=tmb%orbs%inwhichlocreg(iiorb)
      !!        ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      !!        do i=1,ncount
      !!            ist=ist+1
      !!            if (tmb%orbs%spinsgn(iiorb)>0.d0) then
      !!                write(4401,'(a,2i10,f8.1,es16.7)') 'iiorb, ist, spin, vals', iiorb, ist, tmb%orbs%spinsgn(iiorb), hpsi_small(ist)
      !!            else
      !!                write(4402,'(a,2i10,f8.1,es16.7)') 'iiorb, ist, spin, val', iiorb, ist, tmb%orbs%spinsgn(iiorb), hpsi_small(ist)
      !!            end if
      !!        end do
      !!    end do
      !!end if
    
      if (iproc==0) then
          call yaml_map('Preconditioning',.true.)
      end if
    
      call f_release_routine()
    
    
      contains
    
        subroutine calculate_trace_start()
          implicit none
    
          ! Local variables
          integer :: iorb, iiorb
    
          call f_routine(id='calculate_trace_start')
          call timing(iproc,'calctrace_comp','ON')
    
          trH=0.d0
          do iorb=1,tmb%orbs%norbp
             iiorb=tmb%orbs%isorb+iorb
             ii=matrixindex_in_compressed(tmb%linmat%smat(2),iiorb,iiorb)
             trH = trH + tmb%linmat%ham_%matrix_compr(ii-tmb%linmat%smat(2)%isvctrp_tg)
             !write(*,*) 'iorb, tr', iorb, tmb%linmat%ham_%matrix_compr(ii-tmb%linmat%smat(2)%isvctrp_tg)
          end do
          call timing(iproc,'calctrace_comp','OF')
          call timing(iproc,'calctrace_comm','ON')
          if (nproc>1) then
              trH_sendbuf = trH
              call mpiiallred(trH_sendbuf, trH, 1, mpi_sum, bigdft_mpi%mpi_comm, request)
          end if
          call timing(iproc,'calctrace_comm','OF')
    
          call f_release_routine()
    
        end subroutine calculate_trace_start
    
    
        subroutine calculate_trace_finish()
          implicit none
    
          ! Local variables
          integer :: iorb, iiorb
    
          call f_routine(id='calculate_trace_finish')
          call timing(iproc,'calctrace_comm','ON')
          if (nproc>1) then
              call mpiwait(request)
          end if
          call timing(iproc,'calctrace_comm','OF')
    
          call f_release_routine()
    
        end subroutine calculate_trace_finish
    
    
      subroutine communicate_fnrm()
        implicit none
        ! Local variables
        integer :: jproc
        call f_routine(id='communicate_fnrm')
        !!if (nproc > 1) then
        !!   call mpiallred(fnrm, 1, mpi_sum, bigdft_mpi%mpi_comm)
        !!end if
    
        if (nproc>1) then
            fnrm%receivebuf = 0.d0
            fnrm%window = mpiwindow(1, fnrm%receivebuf(1), bigdft_mpi%mpi_comm)
            call mpiaccumulate(origin=fnrm%sendbuf(1), count=1, target_rank=0, &
                 target_disp=int(0,kind=mpi_address_kind),op=mpi_sum, window=fnrm%window)
        else
            fnrm%receivebuf(1) = fnrm%sendbuf(1)
        end if
    
        call f_release_routine()
      end subroutine communicate_fnrm
    
      subroutine communicate_alpha()
        implicit none
        call f_routine(id='communicate_alpha')
        tt=sum(alpha)
        alpha_max=maxval(alpha)
        if (nproc > 1) then
           call mpiallred(tt, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
           call mpiallred(alpha_max, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
        end if
        call f_release_routine()
      end subroutine communicate_alpha
    
    end subroutine calculate_energy_and_gradient_linear



    subroutine hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, at, do_iterative_orthonormalization, sf_per_type, &
               lphiold, alpha, trH, alpha_mean, alpha_max, alphaDIIS, hpsi_small, ortho, psidiff, &
               experimental_mode, order_taylor, max_inversion_error, trH_ref, kernel_best, complete_reset)
      use module_base
      use module_types
      use yaml_output
      use orthonormalization, only: orthonormalizeLocalized, iterative_orthonormalization
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, it
      integer,intent(inout) :: order_taylor
      real(kind=8),intent(in) :: max_inversion_error
      type(localizedDIISParameters), intent(inout) :: ldiis
      type(DFT_wavefunction), target,intent(inout) :: tmb
      type(atoms_data),intent(in) :: at
      logical,intent(in) :: do_iterative_orthonormalization
      integer,dimension(at%astruct%ntypes),intent(in) :: sf_per_type 
      real(kind=8), dimension(tmb%npsidim_orbs), intent(inout) :: lphiold
      real(kind=8), intent(in) :: trH, alpha_mean, alpha_max
      real(kind=8), dimension(tmb%orbs%norbp), intent(inout) :: alpha, alphaDIIS
      real(kind=8), dimension(tmb%npsidim_orbs), intent(inout) :: hpsi_small
      real(kind=8), dimension(tmb%npsidim_orbs), optional,intent(out) :: psidiff
      logical, intent(in) :: ortho, experimental_mode
      real(kind=8),intent(out) :: trH_ref
      real(kind=8),dimension(tmb%linmat%smat(3)%nvctrp_tg*tmb%linmat%smat(3)%nspin),intent(inout) :: kernel_best
      logical,intent(out) :: complete_reset
    
      ! Local variables
      integer :: i, iorb, ilr, ist, iiorb, ncount
      character(len=*), parameter :: subname='hpsitopsi_linear'
      real(kind=8), dimension(:), allocatable :: norm
      real(kind=8) :: ddot, dnrm2, tt
    
      call f_routine(id='hpsitopsi_linear')
    
      call DIISorSD(iproc, it, trH, tmb, ldiis, alpha, alphaDIIS, lphiold, trH_ref, kernel_best, complete_reset)
      if(iproc==0) then
          call yaml_newline()
          call yaml_mapping_open('Optimization',flow=.true.)
          if(ldiis%isx>0) then
              call yaml_map('algorithm','DIIS')
              call yaml_map('history length',ldiis%isx)
              call yaml_map('consecutive failures',ldiis%icountDIISFailureCons)
              call yaml_map('total failures',ldiis%icountDIISFailureTot)
          else
              call yaml_map('algorithm','SD')
              call yaml_map('mean alpha',alpha_mean,fmt='(es9.3)')
              call yaml_map('max alpha',alpha_max,fmt='(es9.3)')
              call yaml_map('consecutive successes',ldiis%icountSDSatur)
          end if
          call yaml_mapping_close()
          call yaml_newline()
      end if
    
      ! Improve the orbitals, depending on the choice made above.
      if (present(psidiff)) call vcopy(tmb%npsidim_orbs, tmb%psi(1), 1, psidiff(1), 1)
      if(.not.ldiis%switchSD) then
          call improveOrbitals(iproc, nproc, tmb, tmb%linmat%smat(1)%nspin, ldiis, alpha, hpsi_small, experimental_mode)
      else
          if (iproc==0) then
              call yaml_warning('no improvement of the orbitals, recalculate gradient')
              call yaml_newline()
          end if
      end if
    
      ! The transposed quantities can now not be used any more...
      if(tmb%can_use_transposed) then
          tmb%can_use_transposed=.false.
      end if
    
    
      if (.not.ortho .and. iproc==0) then
          call yaml_map('Orthogonalization',.false.)
      end if
    
      if(.not.ldiis%switchSD.and.ortho) then
          if (present(psidiff)) then
              do i=1,tmb%npsidim_orbs
                  psidiff(i)=tmb%psi(i)-psidiff(i)
              end do 
          end if
    
          if (do_iterative_orthonormalization) then
              call iterative_orthonormalization(iproc, nproc, 1, order_taylor, at, tmb%linmat%smat(1)%nspin, sf_per_type, tmb)
          else
              call orthonormalizeLocalized(iproc, nproc, order_taylor, max_inversion_error, tmb%npsidim_orbs, tmb%orbs, tmb%lzd, &
                   tmb%linmat%smat(1), tmb%linmat%auxs, tmb%linmat%smat(3), tmb%linmat%auxl, &
                   tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, &
                   tmb%can_use_transposed, ice_obj=tmb%ice_obj)
           end if
          if (iproc == 0) then
              call yaml_map('Orthogonalization',.true.)
          end if
      else if (experimental_mode .or. .not.ldiis%switchSD) then
          ist=1
          do iorb=1,tmb%orbs%norbp
              iiorb=tmb%orbs%isorb+iorb
              ilr=tmb%orbs%inwhichlocreg(iiorb)
              ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
              tt=dnrm2(ncount, tmb%psi(ist), 1)
              tt=1.d0/tt
              call dscal(ncount, tt, tmb%psi(ist), 1)
              ist=ist+ncount
          end do
          if (iproc == 0) then
              call yaml_map('Normalization',.true.)
              !call yaml_map('Normalization',.false.)
          end if
          if (present(psidiff)) then
              do i=1,tmb%npsidim_orbs
                  psidiff(i)=tmb%psi(i)-psidiff(i)
              end do 
          end if
      end if
    
      ! Emit that new wavefunctions are ready.
      if (tmb%c_obj /= 0) then
         call kswfn_emit_psi(tmb, it, 0, iproc, nproc)
      end if
    
      call f_release_routine
    
    end subroutine hpsitopsi_linear
    
    
    
    subroutine build_gradient(iproc, nproc, tmb, target_function, hpsit_c, hpsit_f, hpsittmp_c, hpsittmp_f)
      use module_base
      use module_types
      use sparsematrix_base, only: sparsematrix_malloc_ptr, assignment(=), &
                                   sparsematrix_malloc, SPARSE_TASKGROUP
      use sparsematrix, only: gather_matrix_from_taskgroups_inplace
      use communications_base, only: TRANSPOSE_FULL
      use communications, only: transpose_localized
      use transposed_operations, only: build_linear_combination_transposed
      use public_enums
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, target_function
      type(DFT_wavefunction),intent(inout) :: tmb
      real(kind=8),dimension(tmb%ham_descr%collcom%ndimind_c),intent(inout) :: hpsit_c
      real(kind=8),dimension(7*tmb%ham_descr%collcom%ndimind_f),intent(inout) :: hpsit_f
      real(kind=8),dimension(tmb%ham_descr%collcom%ndimind_c),intent(out) :: hpsittmp_c !<workarray
      real(kind=8),dimension(7*tmb%ham_descr%collcom%ndimind_f),intent(out) :: hpsittmp_f !<workarray
    
      ! Local variables
      integer :: ispin, ishift, iseg, ii, i, ist, ncount, iorb, ilr, isegstart, isegend, ierr, iiorb
      integer,dimension(2) :: irowcol
      real(kind=8),dimension(:),pointer :: kernel_compr_tmp
      real(kind=8),dimension(:),pointer :: matrix_local
      real(kind=8),dimension(:),allocatable :: tmparr
      integer,parameter :: ALLGATHERV=51, GET=52, GLOBAL_MATRIX=101, SUBMATRIX=102
      integer,parameter :: comm_strategy=GET
      integer,parameter :: data_strategy=SUBMATRIX!GLOBAL_MATRIX
    
          call f_routine(id='build_gradient')
          call timing(iproc,'buildgrad_mcpy','ON')
    
          if(tmb%ham_descr%collcom%ndimind_c>0) &
              call vcopy(tmb%ham_descr%collcom%ndimind_c, hpsit_c(1), 1, hpsittmp_c(1), 1)
          if(tmb%ham_descr%collcom%ndimind_f>0) &
              call vcopy(7*tmb%ham_descr%collcom%ndimind_f, hpsit_f(1), 1, hpsittmp_f(1), 1)
    
          if (target_function==TARGET_FUNCTION_IS_HYBRID) then
              kernel_compr_tmp = sparsematrix_malloc_ptr(tmb%linmat%smat(3),iaction=SPARSE_TASKGROUP,id='kernel_compr_tmp')
              do ispin=1,tmb%linmat%smat(3)%nspin
                  !call vcopy(tmb%linmat%smat(3)%nvctr*tmb%linmat%smat(3)%nspin, tmb%linmat%kernel_%matrix_compr(1), 1, kernel_compr_tmp(1), 1)
                  ist = (ispin-1)*tmb%linmat%smat(3)%nvctrp_tg + 1
                  call vcopy(tmb%linmat%smat(3)%nvctrp_tg, tmb%linmat%kernel_%matrix_compr(ist), 1, kernel_compr_tmp(ist), 1)
              end do
              if (data_strategy==GLOBAL_MATRIX) then
                  stop 'build_gradient: option GLOBAL_MATRIX deprecated'
                  !!isegstart = tmb%linmat%smat(3)%istsegline(tmb%linmat%smat(3)%isfvctr+1)
                  !!isegend = tmb%linmat%smat(3)%istsegline(tmb%linmat%smat(3)%isfvctr+tmb%linmat%smat(3)%nfvctrp) + &
                  !!          tmb%linmat%smat(3)%nsegline(tmb%linmat%smat(3)%isfvctr+tmb%linmat%smat(3)%nfvctrp)-1
                  !!matrix_local = f_malloc_ptr(tmb%linmat%smat(3)%nvctrp,id='matrix_local')
                  !!do ispin=1,tmb%linmat%smat(3)%nspin
                  !!    ishift=(ispin-1)*tmb%linmat%smat(3)%nvctr
                  !!    !$omp parallel default(none) &
                  !!    !$omp shared(isegstart,isegend,tmb,matrix_local,kernel_compr_tmp,ishift) &
                  !!    !$omp private(iseg,ii,i,irowcol)
                  !!    !$omp do
                  !!    do iseg=isegstart,isegend
                  !!        ii=tmb%linmat%smat(3)%keyv(iseg)
                  !!        ! A segment is always on one line, therefore no double loop
                  !!        do i=tmb%linmat%smat(3)%keyg(1,1,iseg),tmb%linmat%smat(3)%keyg(2,1,iseg)
                  !!            if(i==tmb%linmat%smat(3)%keyg(1,2,iseg)) then
                  !!                matrix_local(ii-tmb%linmat%smat(3)%isvctr)=0.d0
                  !!            else
                  !!                matrix_local(ii-tmb%linmat%smat(3)%isvctr)=kernel_compr_tmp(ii+ishift)
                  !!            end if
                  !!            ii=ii+1
                  !!        end do
                  !!    end do
                  !!    !$omp end do
                  !!    !$omp end parallel
                  !!    if (nproc>1) then
                  !!         call timing(iproc,'buildgrad_mcpy','OF')
                  !!         call timing(iproc,'buildgrad_comm','ON')
                  !!         !!call mpi_allgatherv(matrix_local(1), tmb%linmat%smat(3)%nvctrp, mpi_double_precision, &
                  !!         !!     tmb%linmat%kernel_%matrix_compr(ishift+1), tmb%linmat%smat(3)%nvctr_par, &
                  !!         !!     tmb%linmat%smat(3)%isvctr_par, mpi_double_precision, &
                  !!         !!     bigdft_mpi%mpi_comm, ierr)
                  !!         if (comm_strategy==ALLGATHERV) then
                  !!             call mpi_allgatherv(matrix_local(1), tmb%linmat%smat(3)%nvctrp, mpi_double_precision, &
                  !!                  tmb%linmat%kernel_%matrix_compr(ishift+1), tmb%linmat%smat(3)%nvctr_par, &
                  !!                  tmb%linmat%smat(3)%isvctr_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
                  !!             call f_free_ptr(matrix_local)
                  !!         else if (comm_strategy==GET) then
                  !!             !!call mpiget(iproc, nproc, bigdft_mpi%mpi_comm, tmb%linmat%smat(3)%nvctrp, matrix_local, &
                  !!             !!     tmb%linmat%smat(3)%nvctr_par, tmb%linmat%smat(3)%isvctr_par, &
                  !!             !!     tmb%linmat%smat(3)%nvctr, tmb%linmat%kernel_%matrix_compr(ishift+1:ishift+tmb%linmat%smat(3)%nvctr))
                  !!             call mpi_get_to_allgatherv(matrix_local(1), tmb%linmat%smat(3)%nvctrp, &
                  !!                  tmb%linmat%kernel_%matrix_compr(ishift+1), &
                  !!                  tmb%linmat%smat(3)%nvctr_par, tmb%linmat%smat(3)%isvctr_par, bigdft_mpi%mpi_comm)
                  !!         else
                  !!             stop 'build_gradient: wrong communication strategy'
                  !!         end if
                  !!         call timing(iproc,'buildgrad_comm','OF')
                  !!         call timing(iproc,'buildgrad_mcpy','ON')
                  !!         if (ispin==tmb%linmat%smat(3)%nspin) call f_free_ptr(matrix_local)
                  !!     else
                  !!         call vcopy(tmb%linmat%smat(3)%nvctr, matrix_local(1), 1, &
                  !!              tmb%linmat%kernel_%matrix_compr(ishift+1), 1)
                  !!     end if
                  !!end do
              else if (data_strategy==SUBMATRIX) then
                  do ispin=1,tmb%linmat%smat(3)%nspin
                      ishift=(ispin-1)*tmb%linmat%smat(3)%nvctrp_tg-tmb%linmat%smat(3)%isvctrp_tg
                      !$omp parallel default(none) &
                      !$omp shared(isegstart,isegend,tmb,matrix_local,kernel_compr_tmp,ishift) &
                      !$omp private(iseg,ii,i,irowcol)
                      !$omp do
                      do iseg=tmb%linmat%smat(3)%istartendseg_t(1),tmb%linmat%smat(3)%istartendseg_t(2)
                          ii=tmb%linmat%smat(3)%keyv(iseg)
                          ! A segment is always on one line, therefore no double loop
                          do i=tmb%linmat%smat(3)%keyg(1,1,iseg),tmb%linmat%smat(3)%keyg(2,1,iseg) !this is too much, but for the moment ok
                              if(i==tmb%linmat%smat(3)%keyg(1,2,iseg)) then
                                  tmb%linmat%kernel_%matrix_compr(ii+ishift)=0.d0
                              else
                                  tmb%linmat%kernel_%matrix_compr(ii+ishift)=kernel_compr_tmp(ii+ishift)
                              end if
                              ii=ii+1
                          end do
                      end do
                      !$omp end do
                      !$omp end parallel
                  end do
              else
                  stop 'build_gradient: wrong data strategy'
              end if
    
    
    
              ist=1
              do iorb=tmb%orbs%isorb+1,tmb%orbs%isorb+tmb%orbs%norbp
                  ilr=tmb%orbs%inwhichlocreg(iorb)
                  if (tmb%orbs%spinsgn(iorb)>0.d0) then
                      ispin=1
                  else
                      ispin=2
                  end if
                  iiorb = mod(iorb-1,tmb%linmat%smat(3)%nfvctr)+1 ! spin-independent index
                  ishift=(ispin-1)*tmb%linmat%smat(3)%nvctr-tmb%linmat%smat(3)%isvctrp_tg
                  isegstart = tmb%linmat%smat(3)%istsegline(iiorb)
                  isegend = tmb%linmat%smat(3)%istsegline(iiorb) + tmb%linmat%smat(3)%nsegline(iiorb) - 1
                  do iseg=isegstart,isegend
                  !do iseg=1, tmb%linmat%smat(3)%nseg
                      ii=tmb%linmat%smat(3)%keyv(iseg)
                      ! A segment is always on one line, therefore no double loop
                      do i=tmb%linmat%smat(3)%keyg(1,1,iseg),tmb%linmat%smat(3)%keyg(2,1,iseg)
                          if (tmb%linmat%smat(3)%keyg(1,2,iseg)/=iiorb) stop 'tmb%linmat%smat(3)%keyg(1,2,iseg)/=iiorb'
                          if(i==tmb%linmat%smat(3)%keyg(1,2,iseg)) then
                              ncount=tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_f
                              !write(*,*) 'iorb, ii, ishift, ist', iorb, ii, ishift, ist
                              call dscal(ncount, kernel_compr_tmp(ii+ishift), tmb%hpsi(ist), 1)
                              ist=ist+ncount
                          end if
                          ii=ii+1
                      end do
                  end do
              end do
              call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
                   TRANSPOSE_FULL, tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
    
              !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(2), tmb%linmat%ham_)
              call build_linear_combination_transposed(tmb%ham_descr%collcom, &
                   tmb%linmat%smat(3), tmb%linmat%auxl, tmb%linmat%kernel_, hpsittmp_c, &
                   hpsittmp_f, .false., hpsit_c, hpsit_f, iproc)
              ! copy correct kernel back
              do ispin=1,tmb%linmat%smat(3)%nspin
                  !call vcopy(tmb%linmat%smat(3)%nvctr*tmb%linmat%smat(3)%nspin, kernel_compr_tmp(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
                  ist = (ispin-1)*tmb%linmat%smat(3)%nvctr + 1
                  call vcopy(tmb%linmat%smat(3)%nvctrp_tg, kernel_compr_tmp(ist), 1, tmb%linmat%kernel_%matrix_compr(ist), 1)
              end do
              call f_free_ptr(kernel_compr_tmp)
          else
              !!tmparr = sparsematrix_malloc(tmb%linmat%smat(2),iaction=SPARSE_FULL,id='tmparr')
              !!call vcopy(tmb%linmat%smat(2)%nvctr, tmb%linmat%ham_%matrix_compr(1), 1, tmparr(1), 1)
              !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%smat(2), tmb%linmat%ham_)
              call build_linear_combination_transposed(tmb%ham_descr%collcom, &
                   tmb%linmat%smat(3), tmb%linmat%auxl, tmb%linmat%kernel_, hpsittmp_c, hpsittmp_f, .true., hpsit_c, hpsit_f, iproc)
              !!call vcopy(tmb%linmat%smat(2)%nvctr, tmparr(1), 1, tmb%linmat%ham_%matrix_compr(1), 1)
              !!call f_free(tmparr)
          end if
    
          call timing(iproc,'buildgrad_mcpy','OF')
          call f_release_routine()
    
    end subroutine build_gradient


    subroutine large_to_small_locreg(iproc, npsidim_orbs_small, npsidim_orbs_large, lzdsmall, lzdlarge, &
           orbs, philarge, phismall)
      use module_base
      use module_types
      use locreg_operations, only: psi_to_locreg2
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
      !!call f_zero(npsidim_orbs_small, phismall(1))
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



end module get_basis
