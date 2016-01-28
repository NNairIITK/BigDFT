!> @file
!!  Routines used by the linear scaling version
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine linearScaling(iproc,nproc,KSwfn,tmb,at,input,shift,rxyz,denspot,rhopotold,nlpsp,GPU,&
           energs,energy,fpulay,infocode,ref_frags,cdft, &
           fdisp, fion)
 
  use module_base
  use module_types
  use module_interfaces, only: allocate_precond_arrays, deallocate_precond_arrays, &
       & getLocalizedBasis, get_coeff, loewdin, sumrho, write_eigenvalues_data, &
       & write_energies, write_orbital_density,inputguessconfinement
  use yaml_output
  use module_fragments
  use constrained_dft
  use diis_sd_optimization
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use communications_base, only: work_transpose, allocate_p2pComms_buffer, &
                                 deallocate_p2pComms_buffer, &
                                 work_transpose_null, allocate_work_transpose, deallocate_work_transpose, &
                                 TRANSPOSE_POST, TRANSPOSE_GATHER
  use communications, only: synchronize_onesided_communication, transpose_localized, untranspose_localized
  !use communications_init, only: orbitals_communicators
  use sparsematrix_base, only: sparse_matrix, matrices, sparse_matrix_null, deallocate_sparse_matrix, &
                               matrices_null, allocate_matrices, deallocate_matrices, &
                               sparsematrix_malloc, sparsematrix_malloc_ptr, assignment(=), SPARSE_FULL, DENSE_FULL, &
                               SPARSE_TASKGROUP, DENSE_PARALLEL, sparsematrix_malloc0
  use matrix_operations, only: deviation_from_unity_parallel
  use sparsematrix, only: gather_matrix_from_taskgroups_inplace, extract_taskgroup_inplace, uncompress_matrix2, &
                          delete_coupling_terms, transform_sparse_matrix, matrix_matrix_mult_wrapper, &
                          uncompress_matrix_distributed2
  use communications, only: transpose_localized, start_onesided_communication
  use sparsematrix_init, only: matrixindex_in_compressed
  use io, only: writemywaves_linear, writemywaves_linear_fragments, write_linear_matrices, write_linear_coefficients
  use postprocessing_linear, only: loewdin_charge_analysis, &
                                   build_ks_orbitals
  use rhopotential, only: updatePotential, sumrho_for_TMBs, corrections_for_negative_charge
  use locreg_operations, only: workarrays_quartic_convolutions,workarr_precond
  use locregs_init, only: small_to_large_locreg
  use public_enums
  use multipole, only: multipole_analysis_driver, projector_for_charge_analysis, &
                       support_function_gross_multipoles, potential_from_charge_multipoles, &
                       calculate_rpowerx_matrices
  use transposed_operations, only: calculate_overlap_transposed
  use matrix_operations, only: overlapPowerGeneral
  use foe, only: fermi_operator_expansion
  use foe_base, only: foe_data_set_real
  use rhopotential, only: full_local_potential
  use transposed_operations, only: calculate_overlap_transposed
  use bounds, only: geocode_buffers
  use orthonormalization, only : orthonormalizeLocalized
  use multipole_base, only: lmax, external_potential_descriptors, deallocate_external_potential_descriptors
  use orbitalbasis
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(atoms_data),intent(inout) :: at
  type(input_variables),intent(in) :: input ! need to hack to be inout for geopt changes
  real(kind=8),dimension(3),intent(in) :: shift
  real(kind=8),dimension(3,at%astruct%nat),intent(inout) :: rxyz
  real(kind=8),dimension(3,at%astruct%nat),intent(out) :: fpulay
  type(DFT_local_fields), intent(inout) :: denspot
  real(gp), dimension(*), intent(inout) :: rhopotold
  type(DFT_PSP_projectors),intent(inout) :: nlpsp
  type(GPU_pointers),intent(inout) :: GPU
  type(energy_terms),intent(inout) :: energs
  real(gp), dimension(:), pointer :: rho,pot
  real(kind=8),intent(out) :: energy
  type(DFT_wavefunction),intent(inout),target :: tmb
  type(DFT_wavefunction),intent(inout),target :: KSwfn
  integer,intent(out) :: infocode
  type(system_fragment), dimension(:), pointer :: ref_frags ! for transfer integrals
  type(cdft_data), intent(inout) :: cdft
  real(kind=8),dimension(3,at%astruct%nat),intent(in) :: fdisp, fion

  !Local variables
  real(kind=8) :: pnrm,trace,trace_old,fnrm_tmb
  integer :: infoCoeff,istat,it_scc,itout,info_scf,i,iorb
  character(len=*), parameter :: subname='linearScaling'
  real(kind=8),dimension(:),allocatable :: rhopotold_out, eval
  real(kind=8) :: energyold, energyDiff, energyoldout, fnrm_pulay, convCritMix, convCritMix_init
  type(localizedDIISParameters) :: ldiis
  type(DIIS_obj) :: ldiis_coeff, vdiis
  logical :: can_use_ham, update_phi, locreg_increased, reduce_conf, orthonormalization_on, update_kernel
  logical :: fix_support_functions
  integer :: itype, istart, nit_lowaccuracy, nit_highaccuracy, k, l
  integer :: ldiis_coeff_hist, nitdmin, lwork
  logical :: ldiis_coeff_changed, overlap_calculated
  integer :: mix_hist, info_basis_functions, nit_scc, cur_it_highaccuracy, nit_scc_changed
  real(kind=8) :: pnrm_out, alpha_mix, ratio_deltas, convcrit_dmin, tt1, tt2, ehart_ps
  logical :: lowaccur_converged, exit_outer_loop, calculate_overlap, invert_overlap_matrix
  real(kind=8),dimension(:),allocatable :: locrad, kernel_orig, ovrlp_large, tmpmat, inv_full, inv_cropped, ovrlp_orig
  real(kind=8),dimension(:),allocatable :: kernel_cropped, ovrlp_cropped, ham_cropped, kernel_foe_cropped, ham_large
  real(kind=8),dimension(:),allocatable :: tmpmat_compr, hamtilde_compr, work
  integer:: target_function, nit_basis, ieval
  logical :: keep_value
  type(workarrays_quartic_convolutions),dimension(:),pointer :: precond_convol_workarrays
  type(workarr_precond),dimension(:),pointer :: precond_workarrays
  type(work_transpose) :: wt_philarge, wt_hpsinoprecond, wt_hphi, wt_phi
  integer,dimension(:,:),allocatable :: ioffset_isf
  integer :: is1, is2, is3, ie1, ie2, ie3, i1, i2, i3, ii, jj, info, ist
  real(kind=8),dimension(:),pointer :: hpsit_c, hpsit_f
  
  real(kind=gp) :: ebs, vgrad_old, vgrad, valpha, vold, vgrad2, vold_tmp, conv_crit_TMB, best_charge_diff, cdft_charge_thresh
  real(kind=gp), allocatable, dimension(:,:) :: coeff_tmp
  integer :: jorb, cdft_it, nelec, iat, ityp, norder_taylor, ispin, ishift
  integer :: dmin_diag_it, dmin_diag_freq, ioffset, nl1, nl2, nl3
  logical :: reorder, rho_negative
  logical :: write_fragments, write_full_system
  real(wp), dimension(:,:,:), pointer :: mom_vec_fake
  type(matrices) :: weight_matrix_
  real(kind=8) :: sign_of_energy_change
  integer :: nit_energyoscillation
  integer(kind=8) :: nsize
  type(work_mpiaccumulate) :: fnrm_work, energs_work
  integer :: ilr, iiorb, iiat
  real(kind=8),dimension(3) :: rr
  real(kind=8),dimension(1,1) :: K_H
  real(kind=8),dimension(4,4) :: K_O
  integer :: j, ind, n, ishifts, ishiftm, iq
  real(kind=8),dimension(:,:),allocatable :: tempmat, all_evals, ham_small, projector_small, ovrlp_small, theta, coeffs
  real(kind=8),dimension(:,:),pointer :: com
  real(kind=8),dimension(:,:),allocatable :: tmat, coeff, ovrlp_full
  real(kind=8),dimension(:),allocatable :: projector_compr
  real(kind=8),dimension(:,:,:),allocatable :: matrixElements, coeff_all, multipoles_out
  real(kind=8),dimension(:,:,:),pointer :: multipoles
  real(kind=8),dimension(:,:,:,:),allocatable :: test_pot
  type(external_potential_descriptors) :: ep
  type(orbital_basis) :: ob
  real(8),dimension(:),allocatable :: rho_tmp, tmparr
  real(8) :: tt, ddot, max_error, mean_error, r2, occ, tot_occ, ef, ef_low, ef_up, q, fac
  type(matrices),dimension(24) :: rpower_matrix
  character(len=20) :: method, do_ortho

  real(kind=8),dimension(:,:),allocatable :: ovrlp_fullp
  real(kind=8) :: max_deviation, mean_deviation, max_deviation_p, mean_deviation_p

  !integer :: ind, ilr, iorbp, iiorb, indg, npsidim_global
  !real(kind=gp), allocatable, dimension(:) :: gpsi, psit_large_c, psit_large_f, kpsit_c, kpsit_f, kpsi, kpsi_small
  !real(kind=gp), pointer, dimension(:) :: gpsi_all, gpsi_virt
  !type(orbitals_data) :: fakeorbs
  !type(comms_cubic) :: comms


  !!rho_tmp = f_malloc(size(denspot%rhov),id='rho_tmp')

  call timing(iproc,'linscalinit','ON')

  call f_routine(id='linear_scaling')

  call allocate_local_arrays()



  ! Allocate the communications buffers needed for the communications of the potential and
  ! post the messages. This will send to each process the part of the potential that this process
  ! needs for the application of the Hamlitonian to all orbitals on that process.
  call allocate_p2pComms_buffer(tmb%comgp)

  pnrm=1.d100
  pnrm_out=1.d100
  energyold=0.d0
  energyoldout=0.d0
  energy=0.d0
  energs%ebs=0.0d0
  target_function=TARGET_FUNCTION_IS_TRACE
  lowaccur_converged=.false.
  info_basis_functions=-1
  fix_support_functions=.false.
  cur_it_highaccuracy=0
  trace_old=0.0d0
  ldiis_coeff_hist=input%lin%dmin_hist_lowaccuracy
  reduce_conf=.false.
  ldiis_coeff_changed = .false.
  orthonormalization_on=.true.
  dmin_diag_it=0
  dmin_diag_freq=-1
  reorder=.false.
  nullify(mom_vec_fake)
  norder_taylor=input%lin%order_taylor
  sign_of_energy_change = -1.d0
  nit_energyoscillation = 0
  keep_value = .false.
  write_fragments = .false.
  write_full_system = .false.
  if (input%lin%output_fragments == OUTPUT_FRAGMENTS_AND_FULL .or. input%lin%output_fragments == OUTPUT_FRAGMENTS_ONLY) then
     write_fragments = .true.
  end if
  if (input%lin%output_fragments == OUTPUT_FRAGMENTS_AND_FULL .or. input%lin%output_fragments == OUTPUT_FULL_ONLY) then
     write_full_system = .true.
  end if


  ! Allocate the communication arrays for the calculation of the charge density.

  if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then  
     ldiis_coeff%alpha_coeff=input%lin%alphaSD_coeff
!!$     call initialize_DIIS_coeff(ldiis_coeff_hist, ldiis_coeff)
!!$     call allocate_DIIS_coeff(tmb, ldiis_coeff)
     if (input%lin%extra_states==0) then
        call DIIS_set(ldiis_coeff_hist,0.1_gp,tmb%orbs%norb*KSwfn%orbs%norbp,1,ldiis_coeff)
     else
        call DIIS_set(ldiis_coeff_hist,0.1_gp,tmb%orbs%norb*tmb%orbs%norbp,1,ldiis_coeff)
     end if
  end if

  tmb%can_use_transposed=.false.
  !nullify(tmb%psit_c)
  !nullify(tmb%psit_f)

  call timing(iproc,'linscalinit','OF')

  ! Check the quality of the input guess
  call check_inputguess()

  call timing(iproc,'linscalinit','ON')

  if(input%lin%scf_mode/=LINEAR_MIXPOT_SIMPLE) then
      call vcopy(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim),rhopotold(1),1,rhopotold_out(1),1)
      !!call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,alpha_mix,denspot%mix,&
      !!     denspot%rhov,1,denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
      !!     at%astruct%cell_dim(1)*at%astruct%cell_dim(2)*at%astruct%cell_dim(3),&
      !!     pnrm,denspot%dpbox%nscatterarr)
  else
      call vcopy(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim),denspot%rhov(1),1,rhopotold_out(1),1)
      call vcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,1) &
            *input%nspin, denspot%rhov(1), 1, rhopotOld(1), 1)
  end if


  !!!! These tests are already done in cluster, keep them for consitency with the test references
  !!!if (iproc==0) call yaml_mapping_open('Checking Communications of Minimal Basis')
  !!!call check_communications_locreg(iproc,nproc,tmb%orbs,input%nspin,tmb%lzd, &
  !!!     tmb%collcom,tmb%linmat%s,tmb%linmat%ovrlp_, &
  !!!     tmb%npsidim_orbs,tmb%npsidim_comp)
  !!!if (iproc==0) call yaml_mapping_close()
  !!!write(*,*) 'after 1st check, sums', sum(tmb%linmat%ovrlp_%matrix_compr), sum(tmb%linmat%kernel_%matrix_compr)

  !!!if (iproc==0) call yaml_mapping_open('Checking Communications of Enlarged Minimal Basis')
  !!!call check_communications_locreg(iproc,nproc,tmb%orbs,input%nspin,tmb%ham_descr%lzd, &
  !!!     tmb%ham_descr%collcom,tmb%linmat%m,tmb%linmat%ham_, &
  !!!     tmb%ham_descr%npsidim_orbs,tmb%ham_descr%npsidim_comp)
  !!!if (iproc ==0) call yaml_mapping_close()
  !!!write(*,*) 'after 2nd check, sums', sum(tmb%linmat%ovrlp_%matrix_compr), sum(tmb%linmat%kernel_%matrix_compr)


  ! CDFT: calculate w_ab here given w(r)
  ! CDFT: first check that we aren't updating the basis at any point and we don't have any low acc iterations
  if (input%lin%constrained_dft) then
     !if (nit_lowaccuracy>0 .or. input%lin%nItBasis_highaccuracy>1) then
     !   stop 'Basis cannot be updated for now in constrained DFT calculations and no low accuracy is allowed'
     !end if

     weight_matrix_ = matrices_null()
     call allocate_matrices(tmb%linmat%m, allocate_full=.false., matname='weight_matrix_', mat=weight_matrix_)
     weight_matrix_%matrix_compr=cdft%weight_matrix_%matrix_compr

     !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
     !!call extract_taskgroup_inplace(tmb%linmat%m, weight_matrix_)
     call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%l,tmb%linmat%m, &
          tmb%linmat%kernel_,weight_matrix_,&
          ebs,tmb%coeff,KSwfn%orbs,tmb%orbs,.false.)
     !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
     !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%m, weight_matrix_)

     !tmb%linmat%denskern_large%matrix_compr = tmb%linmat%kernel_%matrix_compr
     call deallocate_matrices(weight_matrix_)

     call timing(iproc,'linscalinit','OF')
     call timing(iproc,'constraineddft','ON')
     vgrad_old=ebs-cdft%charge

     !!if (iproc==0) write(*,'(a,4(ES16.6e3,2x))') 'N, Tr(KW), Tr(KW)-N, V*(Tr(KW)-N)',&
     !!     cdft%charge,ebs,vgrad_old,cdft%lag_mult*vgrad_old
     if (iproc==0) then
         call yaml_mapping_open('CDFT infos')
         call yaml_map('N',cdft%charge,fmt='(es16.6e3)')
         call yaml_map('Tr(KW)',ebs,fmt='(es16.6e3)')
         call yaml_map('Tr(KW)-N',vgrad_old,fmt='(es16.6e3)')
         call yaml_map('V*(Tr(KW)-N)',cdft%lag_mult*vgrad_old,fmt='(es16.6e3)')
         call yaml_mapping_close()
     end if
     vgrad_old=abs(vgrad_old)
     valpha=0.5_gp
     !best_charge_diff=vgrad_old
     coeff_tmp=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='coeff_tmp')
     cdft_charge_thresh=input%lin%cdft_conv_crit
     call timing(iproc,'constraineddft','OF')
     call timing(iproc,'linscalinit','ON')
  end if


  ! modify tmb%orbs%occup, as we normally use orbs%occup elsewhere
  if (input%lin%extra_states>0) then
     call f_zero(tmb%orbs%norb,tmb%orbs%occup(1))
     call vcopy(KSwfn%orbs%norb, KSwfn%orbs%occup(1), 1, tmb%orbs%occup(1), 1)
     ! occupy the next few states - don't need to preserve the charge as only using for support function optimization
     do iorb=1,tmb%orbs%norb
        if (tmb%orbs%occup(iorb)==1.0_gp) then
           tmb%orbs%occup(iorb)=2.0_gp
        else if (tmb%orbs%occup(iorb)==0.0_gp) then
           do jorb=iorb,min(iorb+input%lin%extra_states-1,tmb%orbs%norb)
             tmb%orbs%occup(jorb)=2.0_gp
           end do
           exit
        end if
     end do
  else
     ! only use tmb%orbs%occup for calculating energy components, otherwise using KSwfn%orbs%occup
     !SM: This is just to make sure that the value of a open shell calculation is equivalent to a closed shell calculations.
     ! Maybe one should change this to 2 and 1...
     if (input%nspin==1) then
         tmb%orbs%occup=1.0d0
     else
         tmb%orbs%occup=0.5d0
     end if
  end if

  ! if we want to ignore read in coeffs and diag at start - EXPERIMENTAL
  ! return to this point - don't need in all fragment cases, just those where we did an extra get_coeff in init
  if ((input%lin%diag_start .or. input%lin%fragment_calculation) .and. (input%inputPsiId .hasattr. 'FILE')) then !==INPUT_PSI_DISK_LINEAR) then
     ! Calculate the charge density.
     !!tmparr = sparsematrix_malloc(tmb%linmat%l,iaction=SPARSE_FULL,id='tmparr')
     !!call vcopy(tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(1), 1, tmparr(1), 1)
     !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
     call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
          tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, denspot%dpbox%ndimrhopot, &
          denspot%rhov, rho_negative)
     !!call vcopy(tmb%linmat%l%nvctr, tmparr(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
     !!call f_free(tmparr)

     if (rho_negative) then
         call corrections_for_negative_charge(iproc, nproc, at, denspot)
         !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
         !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
         !!call clean_rho(iproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
     end if
     ! Calculate the new potential.
     !if(iproc==0) write(*,'(1x,a)') '---------------------------------------------------------------- Updating potential.'
     !if (iproc==0) call yaml_map('update potential',.true.)
     if (iproc==0) call yaml_mapping_open('update pot',flow=.true.)
     call updatePotential(input%nspin,denspot,energs)!%eh,energs%exc,energs%evxc)
     if (iproc==0) call yaml_mapping_close()
  end if

  call timing(iproc,'linscalinit','OF')


  ! Add one iteration if no low accuracy is desired since we need then a first fake iteration, with istart=0
  istart = min(1,nit_lowaccuracy)
  infocode=0 !default value
  ! This is the main outer loop. Each iteration of this loop consists of a first loop in which the basis functions
  ! are optimized and a consecutive loop in which the density is mixed.
  if (iproc==0) then
      call yaml_comment('Self-Consistent Cycle',hfill='-')
      call yaml_sequence_open('Ground State Optimization')
  end if
  outerLoop: do itout=istart,nit_lowaccuracy+nit_highaccuracy

      call timing(iproc,'linscalinit','ON')
      if (input%lin%nlevel_accuracy==2) then
          ! Check whether the low accuracy part (i.e. with strong confining potential) has converged.
          call check_whether_lowaccuracy_converged(itout, nit_lowaccuracy, input%lin%lowaccuracy_conv_crit, &
               lowaccur_converged, pnrm_out)
          !!if (iproc==0) write(*,*) 'lowaccur_converged,associated(precond_workarrays)',lowaccur_converged,associated(precond_workarrays)
          !!if (lowaccur_converged) then
          !!    call deallocate_precond_arrays(tmb%orbs, tmb%lzd, precond_convol_workarrays, precond_workarrays)
          !!end if
          ! Set all remaining variables that we need for the optimizations of the basis functions and the mixing.
          call set_optimization_variables(input, at, tmb%orbs, tmb%lzd%nlr, tmb%orbs%onwhichatom, tmb%confdatarr, &
               convCritMix_init, lowaccur_converged, nit_scc, mix_hist, alpha_mix, locrad, target_function, nit_basis, &
               convcrit_dmin, nitdmin, conv_crit_TMB)
          if (itout==1) then
              call allocate_precond_arrays(tmb%orbs, tmb%lzd, tmb%confdatarr, precond_convol_workarrays, precond_workarrays)
              wt_philarge = work_transpose_null()
              wt_hpsinoprecond = work_transpose_null()
              wt_hphi = work_transpose_null()
              wt_phi = work_transpose_null()
              call allocate_work_transpose(nproc, tmb%ham_descr%collcom, wt_philarge)
              call allocate_work_transpose(nproc, tmb%ham_descr%collcom, wt_hpsinoprecond)
              call allocate_work_transpose(nproc, tmb%ham_descr%collcom, wt_hphi)
              call allocate_work_transpose(nproc, tmb%collcom, wt_phi)
          end if
          !!if (iproc==0) write(*,*) 'AFTER: lowaccur_converged,associated(precond_workarrays)',lowaccur_converged,associated(precond_workarrays)
          if (keep_value) then
              nit_scc = nit_scc_changed
          end if
      else if (input%lin%nlevel_accuracy==1 .and. itout==1) then
          call set_variables_for_hybrid(iproc, tmb%lzd%nlr, input, at, tmb%orbs, &
               lowaccur_converged, tmb%damping_factor_confinement, tmb%confdatarr, &
               target_function, nit_basis, nit_scc, mix_hist, locrad, alpha_mix, convCritMix_init, conv_crit_TMB)
               convcrit_dmin=input%lin%convCritDmin_highaccuracy
               nitdmin=input%lin%nItdmin_highaccuracy
          call allocate_precond_arrays(tmb%orbs, tmb%lzd, tmb%confdatarr, precond_convol_workarrays, precond_workarrays)
          wt_philarge = work_transpose_null()
          wt_hpsinoprecond = work_transpose_null()
          wt_hphi = work_transpose_null()
          wt_phi = work_transpose_null()
          call allocate_work_transpose(nproc, tmb%ham_descr%collcom, wt_philarge)
          call allocate_work_transpose(nproc, tmb%ham_descr%collcom, wt_hpsinoprecond)
          call allocate_work_transpose(nproc, tmb%ham_descr%collcom, wt_hphi)
          call allocate_work_transpose(nproc, tmb%collcom, wt_phi)
      end if
      call timing(iproc,'linscalinit','OF')

      ! Do one fake iteration if no low accuracy is desired.
      if(nit_lowaccuracy==0 .and. itout==0) then
          lowaccur_converged=.false.
          cur_it_highaccuracy=0
      end if

      if(lowaccur_converged) cur_it_highaccuracy=cur_it_highaccuracy+1

      if(cur_it_highaccuracy==1) then
          if (iproc==0) then
              call yaml_comment('Adjustments for high accuracy',hfill='=')
          end if
          ! nit_scc might have been overwritten, so recalculate it and set the ! flag to false
          if (input%lin%nlevel_accuracy==2) then
              ! Set all remaining variables that we need for the optimizations of the basis functions and the mixing.
              call set_optimization_variables(input, at, tmb%orbs, tmb%lzd%nlr, tmb%orbs%onwhichatom, tmb%confdatarr, &
                   convCritMix_init, lowaccur_converged, nit_scc, mix_hist, alpha_mix, locrad, target_function, nit_basis, &
                   convcrit_dmin, nitdmin, conv_crit_TMB)
              keep_value = .false.
          end if
          ! Adjust the confining potential if required.
          call adjust_locregs_and_confinement(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
               at, input, rxyz, KSwfn, tmb, denspot, nlpsp, ldiis, locreg_increased, lowaccur_converged, locrad)
          orthonormalization_on=.true.

          call deallocate_precond_arrays(tmb%orbs, tmb%lzd, precond_convol_workarrays, precond_workarrays)
          call deallocate_work_transpose(wt_philarge)
          call deallocate_work_transpose(wt_hpsinoprecond)
          call deallocate_work_transpose(wt_hphi)
          call deallocate_work_transpose(wt_phi)

          call allocate_precond_arrays(tmb%orbs, tmb%lzd, tmb%confdatarr, precond_convol_workarrays, precond_workarrays)
          wt_philarge = work_transpose_null()
          wt_hpsinoprecond = work_transpose_null()
          wt_hphi = work_transpose_null()
          wt_phi = work_transpose_null()
          call allocate_work_transpose(nproc, tmb%ham_descr%collcom, wt_philarge)
          call allocate_work_transpose(nproc, tmb%ham_descr%collcom, wt_hpsinoprecond)
          call allocate_work_transpose(nproc, tmb%ham_descr%collcom, wt_hphi)
          call allocate_work_transpose(nproc, tmb%collcom, wt_phi)

          ! is this really necessary if the locrads haven't changed?  we should check this!
          ! for now for CDFT don't do the extra get_coeffs, as don't want to add extra CDFT loop here
          if (target_function==TARGET_FUNCTION_IS_HYBRID) then
              if (iproc==0) write(*,*) 'WARNING: COMMENTED THESE LINES'
          else
             ! Calculate a new kernel since the old compressed one has changed its shape due to the locrads
             ! being different for low and high accuracy.
             update_phi=.true.
             tmb%can_use_transposed=.false.   !check if this is set properly!
             ! NB nothing is written to screen for this get_coeff
             if (.not. input%lin%constrained_dft) then
                ! Check whether the overlap matrix must be calculated and inverted (otherwise it has already been done)
                !!calculate_overlap = ((update_phi .and. .not.input%correction_co_contra) .or. cur_it_highaccuracy==1)
                !!invert_overlap_matrix = ((input%method_updatekernel/=UPDATE_BY_FOE .and. &
                !!                         input%method_updatekernel/=UPDATE_BY_RENORMALIZATION) .or. &
                !!                         cur_it_highaccuracy==1)
                calculate_overlap = ((update_phi .and. .not.input%correction_co_contra) .or. cur_it_highaccuracy==1)
                invert_overlap_matrix = calculate_overlap
                !invert_overlap_matrix = (.not.(target_function==TARGET_FUNCTION_IS_HYBRID .and. &
                !                          (input%method_updatekernel==UPDATE_BY_FOE .or. &
                !                         input%method_updatekernel==UPDATE_BY_RENORMALIZATION)))
                !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
                call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                     infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,calculate_overlap,invert_overlap_matrix,update_phi,&
                     .true.,input%lin%extra_states,itout,0,0,norder_taylor,input%lin%max_inversion_error,&
                     input%calculate_KS_residue,input%calculate_gap,energs_work,.false.,input%lin%coeff_factor,&
                     input%lin%pexsi_npoles, input%lin%pexsi_mumin,input%lin%pexsi_mumax,input%lin%pexsi_mu, &
                     input%lin%pexsi_temperature,input%lin%pexsi_tol_charge, &
                     convcrit_dmin,nitdmin,input%lin%curvefit_dmin,ldiis_coeff)
                !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
             end if
          end if

          ! Some special treatement if we are in the high accuracy part
          call adjust_DIIS_for_high_accuracy(input, denspot, lowaccur_converged, &
               ldiis_coeff_hist, ldiis_coeff_changed)
      end if


      if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then 
         !call initialize_DIIS_coeff(ldiis_coeff_hist, ldiis_coeff)
         call DIIS_free(ldiis_coeff)
         if (input%lin%extra_states==0) then
            call DIIS_set(ldiis_coeff_hist,0.1_gp,tmb%orbs%norb*KSwfn%orbs%norbp,1,ldiis_coeff)
         else
            call DIIS_set(ldiis_coeff_hist,0.1_gp,tmb%orbs%norb*tmb%orbs%norbp,1,ldiis_coeff)
         end if

         ! need to reallocate DIIS matrices to adjust for changing history length
!!$         if (ldiis_coeff_changed) then
!!$            call deallocateDIIS(ldiis_coeff)
!!$            call allocate_DIIS_coeff(tmb, ldiis_coeff)
!!$            ldiis_coeff_changed = .false.
!!$         end if
      end if

      if(itout>1 .or. (nit_lowaccuracy==0 .and. itout==1)) then
          call deallocateDIIS(ldiis)
      end if
      if (lowaccur_converged) then
          call initializeDIIS(input%lin%DIIS_hist_highaccur, tmb%lzd, tmb%orbs, ldiis)
      else
          call initializeDIIS(input%lin%DIIS_hist_lowaccur, tmb%lzd, tmb%orbs, ldiis)
      end if
      if(itout==1) then
          ldiis%alphaSD=input%lin%alphaSD
          ldiis%alphaDIIS=input%lin%alphaDIIS
      end if

      ! Do nothing if no low accuracy is desired.
      if (nit_lowaccuracy==0 .and. itout==0) then
          !!if (associated(tmb%psit_c)) then
          !!    call f_free_ptr(tmb%psit_c)
          !!end if
          !!if (associated(tmb%psit_f)) then
          !!    call f_free_ptr(tmb%psit_f)
          !!end if
          tmb%can_use_transposed=.false.
          if (iproc==0) then
              call yaml_sequence(advance='no')
              call yaml_sequence_open('fake iteration',label=&
                                 'it_fake'//trim(adjustl(yaml_toa(0,fmt='(i3.3)'))))
              call yaml_sequence(label='final_fake'//trim(adjustl(yaml_toa(0,fmt='(i3.3)'))),advance='no')
              call yaml_mapping_open(flow=.true.)
              call yaml_map('fake iteration','bridge low accuracy')
              call yaml_mapping_close
              call yaml_sequence_close()
          end if
          cycle outerLoop
      end if

      ! The basis functions shall be optimized except if it seems to saturate
      update_phi=.true.
      if(fix_support_functions) then
          update_phi=.false.
          tmb%can_use_transposed=.false.   !check if this is set properly!
      end if

      ! Improve the trace minimizing orbitals.
       if(update_phi) then
           if (target_function==TARGET_FUNCTION_IS_HYBRID .and. reduce_conf) then
               if (input%lin%reduce_confinement_factor>0.d0) then
                   !if (iproc==0) write(*,'(1x,a,es8.1)') 'Multiply the confinement prefactor by',input%lin%reduce_confinement_factor
                   if (iproc==0) call yaml_map('multiplicator for the confinement',input%lin%reduce_confinement_factor)
                   tmb%confdatarr(:)%prefac=input%lin%reduce_confinement_factor*tmb%confdatarr(:)%prefac
               else
                   if (ratio_deltas<=1.d0 .and. ratio_deltas>0.d0) then
                       !if (iproc==0) write(*,'(1x,a,es8.1)') 'Multiply the confinement prefactor by',ratio_deltas
                       !if (iproc==0) call yaml_map('multiplicator for the confinement',ratio_deltas)
                       if (iproc==0) call yaml_map('multiplicator for the confinement',ratio_deltas**2)
                       tmb%confdatarr(:)%prefac=(ratio_deltas**2)*tmb%confdatarr(:)%prefac
                   else if (ratio_deltas>1.d0) then
                       !if (iproc==0) write(*,*) 'WARNING: ratio_deltas>1!. Using 0.5 instead'
                       !if (iproc==0) write(*,'(1x,a,es8.1)') 'Multiply the confinement prefactor by',0.5d0
                       if (iproc==0) call yaml_warning('ratio_deltas>1, using 1.0 instead')
                       if (iproc==0) call yaml_newline()
                       if (iproc==0) call yaml_map('multiplicator for the confinement',1.0d0)
                       tmb%confdatarr(:)%prefac=1.0d0*tmb%confdatarr(:)%prefac
                   else if (ratio_deltas<=0.d0) then
                       !if (iproc==0) write(*,*) 'WARNING: ratio_deltas<=0.d0!. Using 0.5 instead'
                       !if (iproc==0) write(*,'(1x,a,es8.1)') 'Multiply the confinement prefactor by',0.5d0
                       if (iproc==0) call yaml_warning('ratio_deltas<=0.d0, using 0.5 instead')
                       if (iproc==0) call yaml_newline()
                       if (iproc==0) call yaml_map('multiplicator for the confinement',0.5d0)
                       tmb%confdatarr(:)%prefac=0.5d0*tmb%confdatarr(:)%prefac
                   end if
               end if
               !if (iproc==0) write(*,'(a,es18.8)') 'tmb%confdatarr(1)%prefac',tmb%confdatarr(1)%prefac
               !if (iproc==0) write(*,*) "WARNING: DON'T LET THE PREFACTOR GO BELOW 1.D-5"
               !tmb%confdatarr(:)%prefac=max(tmb%confdatarr(:)%prefac,2.4d-5)
           end if

           !!if (pnrm_out<5.d-9) then
           !!    if (iproc==0) write(*,*) 'outswitch off ortho'
           !!    orthonormalization_on=.false.
           !!end if
           !!if (sum(tmb%confdatarr(:)%prefac)==0.d0) then
           !!    if (iproc==0) write(*,*) 'WARNING: modifi nit_basis'
           !!    nit_basis=100
           !!end if
           !if (iproc==0) write(*,*) 'upper bound for prefac: 1.d-5'
           !tmb%confdatarr(:)%prefac=max(tmb%confdatarr(:)%prefac,1.d-5)

           !if (itout<=20) then
           !    if (iproc==0) write(*,*) 'set ldiis%isx=0)'
           !    ldiis%isx=0
           !end if
           if (input%experimental_mode) then
               if (iproc==0) call yaml_warning('No orthogonalizing of the support functions')
               !if (iproc==0) write(*,*) 'WARNING: set orthonormalization_on to false'
               orthonormalization_on=.false.
           end if
           if (iproc==0) then
               call yaml_comment('support function optimization',hfill='=')
           end if
           if (iproc==0) then
               call yaml_sequence(advance='no')
               call yaml_sequence_open('support function optimization',label=&
                              'it_supfun'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)'))))
           end if
           !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
           if (input%lin%constrained_dft) then
              call getLocalizedBasis(iproc,nproc,at,KSwfn%orbs,rxyz,denspot,GPU,trace,trace_old,fnrm_tmb,&
                  info_basis_functions,nlpsp,input%lin%scf_mode,ldiis,input%SIC,tmb,energs,&
                  input%lin%iterative_orthogonalization,input%lin%norbsPerType,&
                  input%lin%nItPrecond,target_function,input%lin%correctionOrthoconstraint,&
                  nit_basis,&
                  ratio_deltas,orthonormalization_on,input%lin%extra_states,itout,conv_crit_TMB,input%experimental_mode,&
                  input%lin%early_stop, input%lin%gnrm_dynamic, input%lin%min_gnrm_for_dynamic, &
                  can_use_ham, norder_taylor, input%lin%max_inversion_error, input%kappa_conv,&
                  input%correction_co_contra, &
                  precond_convol_workarrays, precond_workarrays, &
                  wt_philarge, wt_hpsinoprecond, wt_hphi, wt_phi, fnrm_work, energs_work, input%lin%fragment_calculation, &
                  cdft, input%frag, ref_frags)
           else
              call getLocalizedBasis(iproc,nproc,at,KSwfn%orbs,rxyz,denspot,GPU,trace,trace_old,fnrm_tmb,&
                  info_basis_functions,nlpsp,input%lin%scf_mode,ldiis,input%SIC,tmb,energs,&
                  input%lin%iterative_orthogonalization,input%lin%norbsPerType,&
                  input%lin%nItPrecond,target_function,input%lin%correctionOrthoconstraint,&
                  nit_basis,&
                  ratio_deltas,orthonormalization_on,input%lin%extra_states,itout,conv_crit_TMB,input%experimental_mode,&
                  input%lin%early_stop, input%lin%gnrm_dynamic, input%lin%min_gnrm_for_dynamic, &
                  can_use_ham, norder_taylor, input%lin%max_inversion_error, input%kappa_conv,&
                  input%correction_co_contra, precond_convol_workarrays, precond_workarrays, &
                  wt_philarge, wt_hpsinoprecond, wt_hphi, wt_phi, fnrm_work, energs_work, input%lin%fragment_calculation)
              !if (iproc==0) call yaml_scalar('call boundary analysis')
              call get_boundary_weight(iproc, nproc, tmb%orbs, tmb%lzd, at, &
                   input%crmult, tmb%npsidim_orbs, tmb%psi, 1.d-2)
           end if
           !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
           reduce_conf=.true.
           if (iproc==0) then
               call yaml_sequence_close()
           end if

    !       !update weight matrix following basis optimization (could add check to ensure this is really necessary)
    !       if (input%lin%constrained_dft) then
    !          if (trim(cdft%method)=='fragment_density') then ! fragment density approach
    !             if (input%lin%diag_start) stop 'fragment_density not allowed'
    !          else if (trim(cdft%method)=='lowdin') then ! direct weight matrix approach
    !             !should already have overlap matrix so no need to recalculate
    !             call calculate_weight_matrix_lowdin_wrapper(cdft,tmb,input,ref_frags,.false.,input%lin%order_taylor)
    !          else 
    !             stop 'Error invalid method for calculating CDFT weight matrix'
    !          end if
    !       end if

           tmb%can_use_transposed=.false. !since basis functions have changed...

           if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) ldiis_coeff%alpha_coeff=input%lin%alphaSD_coeff !reset to default value

           ! I think this is causing a memory leak somehow in certain cases (possibly only with fragment calculations?)
           if ((input%inputPsiId .hasattr. 'MEMORY') &!input%inputPsiId==101 
                .and. info_basis_functions<=-2 .and. itout==1 .and. (.not.input%lin%fragment_calculation)) then
               ! There seem to be some convergence problems after a restart. Better to quit
               ! and start with a new AO input guess.
               if (iproc==0) write(*,'(1x,a)') 'There are convergence problems after the restart. &
                                                &Start over again with an AO input guess.'
               !!if (associated(tmb%psit_c)) then
               !!    call f_free_ptr(tmb%psit_c)
               !!end if
               !!if (associated(tmb%psit_f)) then
               !!    call f_free_ptr(tmb%psit_f)
               !!end if
               infocode=2
               exit outerLoop
           end if
       end if

    ! Check whether we can use the Hamiltonian matrix from the TMB optimization
    ! for the first step of the coefficient optimization
    !can_use_ham=.true.
    ! can_use_ham was set in getLocalizedBasis
    if(target_function==TARGET_FUNCTION_IS_TRACE) then
        do itype=1,at%astruct%ntypes
            if(input%lin%potentialPrefac_lowaccuracy(itype)/=0.d0) then
                can_use_ham=.false.
                exit
            end if
        end do
    else if(target_function==TARGET_FUNCTION_IS_ENERGY) then
        do itype=1,at%astruct%ntypes
            if(input%lin%potentialPrefac_highaccuracy(itype)/=0.d0) then
                can_use_ham=.false.
                exit
            end if
        end do
    !!else if(target_function==TARGET_FUNCTION_IS_HYBRID) then
    !!    do itype=1,at%astruct%ntypes
    !!        if(input%lin%potentialPrefac_lowaccuracy(itype)/=0.d0) then
    !!            can_use_ham=.false.
    !!            exit
    !!        end if
    !!    end do
    end if

    if (iproc==0) then
        call yaml_comment('kernel optimization',hfill='=')
        !call yaml_sequence(advance='no')
        !if (input%lin%constrained_dft) then
        !    call yaml_mapping_open('kernel optimization',label=&
        !         'it_kernel'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//&
        !         '_'//trim(adjustl(yaml_toa(cdft_it,fmt='(i3.3)'))))
        !else
        !    call yaml_mapping_open('kernel optimization',label=&
        !         'it_kernel'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)'))))
        !end if
        call yaml_sequence(advance='no')
        if (input%lin%constrained_dft) then
            call yaml_sequence_open('kernel optimization',label=&
                 'it_kernel'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//&
                 '_'//trim(adjustl(yaml_toa(cdft_it,fmt='(i3.3)'))))
        else
            call yaml_sequence_open('kernel optimization',label=&
                 'it_kernel'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)'))))
        end if
    end if
    ! @NEW: adjust the convergence criterion for the kernel optimization
    ! The better the support functions are converged, the better the kernel sould be converged
    ! make this optional, otherwise sometimes the threshold becomes too long and it takes a really long time to converge
    if (input%adjust_kernel_threshold) then
       convCritMix = convCritMix_init*fnrm_tmb
    else
       convCritMix = convCritMix_init
    end if

    call scf_kernel(nit_scc, .false., update_phi)

    ! Write the final results
    if (iproc==0) then
        if (input%lin%constrained_dft) then
            call yaml_sequence(label='final_kernel'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//&
                 '_'//trim(adjustl(yaml_toa(cdft_it,fmt='(i3.3)'))),advance='no')
        else
            call yaml_sequence(label='final_kernel'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)'))),advance='no')
        end if
        call yaml_mapping_open(flow=.true.)
        call yaml_comment('iter:'//yaml_toa(it_scc,fmt='(i6)'),hfill='-')
        call printSummary()
        call yaml_mapping_close() !iteration
        call yaml_flush_document()
        !call bigdft_utils_flush(unit=6)
    end if

    ! Close sequence for the optimization steps
    if (iproc==0) then
        call yaml_sequence_close()
    end if

    if(tmb%can_use_transposed) then
        !!call f_free_ptr(tmb%psit_c)
        !!call f_free_ptr(tmb%psit_f)
        tmb%can_use_transposed=.false.
    end if

    call print_info(.false.)

    energyoldout=energy

    if (input%intermediate_forces) then
        call intermediate_forces()
    end if

    call check_for_exit()
    if(exit_outer_loop) exit outerLoop

    if(pnrm_out<input%lin%support_functions_converged.and.lowaccur_converged) then
        !if(iproc==0) write(*,*) 'fix the support functions from now on'
        if (iproc==0) call yaml_map('fix the support functions from now on',.true.)
        fix_support_functions=.true.
    end if

  end do outerLoop


  call deallocate_precond_arrays(tmb%orbs, tmb%lzd, precond_convol_workarrays, precond_workarrays)

  ! Moved down for tests
  call deallocate_work_transpose(wt_philarge)
  call deallocate_work_transpose(wt_hpsinoprecond)
  call deallocate_work_transpose(wt_hphi)
  call deallocate_work_transpose(wt_phi)


  if (input%wf_extent_analysis) then
      ioffset_isf = f_malloc((/3,tmb%orbs%norbp/),id='ioffset_isf')
      do iorb=1,tmb%orbs%norbp
          iiorb = tmb%orbs%isorb + iorb
          ilr = tmb%orbs%inwhichlocreg(iiorb)
          call geocode_buffers(tmb%lzd%Llr(ilr)%geocode, tmb%lzd%glr%geocode, nl1, nl2, nl3)
          ioffset_isf(1,iorb) = tmb%lzd%llr(ilr)%nsi1 - nl1 - 1
          ioffset_isf(2,iorb) = tmb%lzd%llr(ilr)%nsi2 - nl2 - 1
          ioffset_isf(3,iorb) = tmb%lzd%llr(ilr)%nsi3 - nl3 - 1
          !write(*,'(a,i8,2es16.8)') 'iorb, rxyzConf(3), locregcenter(3)', iorb, tmb%confdatarr(iorb)%rxyzConf(3), tmb%lzd%llr(ilr)%locregcenter(3)
      end do
      call analyze_wavefunctions('Support functions extent analysis', 'local', &
           tmb%lzd, tmb%orbs, tmb%npsidim_orbs, tmb%psi, ioffset_isf)
      call f_free(ioffset_isf)
  end if


  if (input%write_orbitals>0) then
      call build_ks_orbitals(iproc, nproc, tmb, KSwfn, at, rxyz, denspot, GPU, &
               energs, nlpsp, input, norder_taylor,&
               energy, energyDiff, energyold)
      !call write_orbital_density(iproc, .false., input%lin%plotBasisFunctions, 'KS', &
      !     KSwfn%orbs%npsidim_orbs, KSwfn%psi, KSwfn%orbs, KSwfn%lzd, at)

      !ioffset_isf = f_malloc((/3,orbs%norbp/),id='ioffset_isf')
      !do iorb=1,orbs%norbp
      !    !iiorb = tmb%orbs%isorb + iorb
      !    !ilr = tmb%orbs%inwhichlocreg(iiorb)
      !    !call geocode_buffers(tmb%lzd%Llr(ilr)%geocode, tmb%lzd%glr%geocode, nl1, nl2, nl3)
      !    ioffset_isf(1,iorb) = 0 !tmb%lzd%llr(ilr)%nsi1 - nl1 - 1
      !    ioffset_isf(2,iorb) = 0 !tmb%lzd%llr(ilr)%nsi2 - nl2 - 1
      !    ioffset_isf(3,iorb) = 0 !tmb%lzd%llr(ilr)%nsi3 - nl3 - 1
      !    !write(*,'(a,3es16.8)') 'iorb, rxyzConf(3), locregcenter(3)', iorb, tmb%confdatarr(iorb)%rxyzConf(3), tmb%lzd%llr(ilr)%locregcenter(3)
      !end do
      !call analyze_wavefunctions('global', tmb%lzd, orbs, KSwfn%orbs%npsidim_orbs, %psi, ioffset_isf)
      !call f_free(ioffset_isf)
  end if




  !TEMPORARY, to be cleaned/removed
  !!!missing occs but otherwise ok? (at least doesn't crash, add occs and recheck by comparing with lowdin/do linear to cubic
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! here as I guess we're deallocating something we need later
  !call nullify_orbitals_data(fakeorbs)
  !call copy_orbitals_data(tmb%orbs, fakeorbs, subname)
  !call orbitals_communicators(iproc, nproc, tmb%lzd%glr, fakeorbs, comms)
  !
  !npsidim_global=max(tmb%orbs%norbp*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f), &
  !                   tmb%orbs%norb*comms%nvctr_par(iproc,0)*fakeorbs%nspinor)
  !gpsi_all=f_malloc_ptr(npsidim_global,id='gpsi_all')
  !call build_ks_orbitals_laura_tmp(iproc, nproc, tmb, KSwfn, at, rxyz, denspot, GPU, &
  !         energs, nlpsp, input, norder_taylor, &
  !         energy, energyDiff, energyold, npsidim_global, gpsi_all)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   update_kernel=.false.

!switch this off for now, needs cleaning and stabilizing
if (.false.) then
   ! not sure if we always want to do this when writing to disk? or all fragment calculations?
   if (mod(input%lin%plotBasisFunctions,10) /= WF_FORMAT_NONE .and. input%lin%fragment_calculation) then

  ovrlp_fullp = sparsematrix_malloc(tmb%linmat%l,iaction=DENSE_PARALLEL,id='ovrlp_fullp')
  max_deviation=0.d0
  mean_deviation=0.d0
  do ispin=1,tmb%linmat%s%nspin
      ishift=(ispin-1)*tmb%linmat%s%nvctrp_tg
      call uncompress_matrix_distributed2(iproc, tmb%linmat%s, DENSE_PARALLEL, &
           tmb%linmat%ovrlp_%matrix_compr(ishift+1:), ovrlp_fullp)
      call deviation_from_unity_parallel(iproc, nproc, tmb%linmat%s%nfvctr, tmb%linmat%s%nfvctrp, &
           tmb%linmat%s%isfvctr, ovrlp_fullp, &
           tmb%linmat%s, max_deviation_p, mean_deviation_p)
      max_deviation = max_deviation + max_deviation_p/real(tmb%linmat%s%nspin,kind=8)
      mean_deviation = mean_deviation + mean_deviation_p/real(tmb%linmat%s%nspin,kind=8)
  end do
  call f_free(ovrlp_fullp)
  if (iproc==0) then
      call yaml_map('max dev from unity',max_deviation,fmt='(es9.2)')
      call yaml_map('mean dev from unity',mean_deviation,fmt='(es9.2)')
  end if

tmb%can_use_transposed=.false.
      !call orthonormalizeLocalized(iproc, nproc, norder_taylor, input%lin%max_inversion_error, tmb%npsidim_orbs, tmb%orbs, tmb%lzd, &
      !     tmb%linmat%s, tmb%linmat%l, tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)
      call orthonormalizeLocalized(iproc, nproc, norder_taylor, input%lin%max_inversion_error, tmb%npsidim_orbs, &
           tmb%orbs, tmb%lzd, tmb%linmat%s, tmb%linmat%l, tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, &
           tmb%can_use_transposed)

  call deallocate_matrices(tmb%linmat%ovrlp_)
  tmb%linmat%ovrlp_ = matrices_null()
  call allocate_matrices(tmb%linmat%s, allocate_full=.false., matname='tmb%linmat%ovrlp_', mat=tmb%linmat%ovrlp_)
  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, &
       tmb%linmat%s, tmb%linmat%ovrlp_)

  ovrlp_fullp = sparsematrix_malloc(tmb%linmat%l,iaction=DENSE_PARALLEL,id='ovrlp_fullp')
  max_deviation=0.d0
  mean_deviation=0.d0
  do ispin=1,tmb%linmat%s%nspin
      ishift=(ispin-1)*tmb%linmat%s%nvctrp_tg
      call uncompress_matrix_distributed2(iproc, tmb%linmat%s, DENSE_PARALLEL, &
           tmb%linmat%ovrlp_%matrix_compr(ishift+1:), ovrlp_fullp)
      call deviation_from_unity_parallel(iproc, nproc, tmb%linmat%s%nfvctr, tmb%linmat%s%nfvctrp, &
           tmb%linmat%s%isfvctr, ovrlp_fullp, &
           tmb%linmat%s, max_deviation_p, mean_deviation_p)
      max_deviation = max_deviation + max_deviation_p/real(tmb%linmat%s%nspin,kind=8)
      mean_deviation = mean_deviation + mean_deviation_p/real(tmb%linmat%s%nspin,kind=8)
  end do
  call f_free(ovrlp_fullp)
  if (iproc==0) then
      call yaml_map('max dev from unity',max_deviation,fmt='(es9.2)')
      call yaml_map('mean dev from unity',mean_deviation,fmt='(es9.2)')
  end if

tmb%can_use_transposed=.false.
      !call orthonormalizeLocalized(iproc, nproc, norder_taylor, input%lin%max_inversion_error, tmb%npsidim_orbs, tmb%orbs, tmb%lzd, &
      !     tmb%linmat%s, tmb%linmat%l, tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)
      call orthonormalizeLocalized(iproc, nproc, norder_taylor, input%lin%max_inversion_error, tmb%npsidim_orbs, &
           tmb%orbs, tmb%lzd, tmb%linmat%s, tmb%linmat%l, tmb%collcom, tmb%orthpar, tmb%psi, &
           tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)

  call deallocate_matrices(tmb%linmat%ovrlp_)
  tmb%linmat%ovrlp_ = matrices_null()
  call allocate_matrices(tmb%linmat%s, allocate_full=.false., matname='tmb%linmat%ovrlp_', mat=tmb%linmat%ovrlp_)
  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, &
       tmb%linmat%s, tmb%linmat%ovrlp_)

  ovrlp_fullp = sparsematrix_malloc(tmb%linmat%l,iaction=DENSE_PARALLEL,id='ovrlp_fullp')
  max_deviation=0.d0
  mean_deviation=0.d0
  do ispin=1,tmb%linmat%s%nspin
      ishift=(ispin-1)*tmb%linmat%s%nvctrp_tg
      call uncompress_matrix_distributed2(iproc, tmb%linmat%s, DENSE_PARALLEL, &
           tmb%linmat%ovrlp_%matrix_compr(ishift+1:), ovrlp_fullp)
      call deviation_from_unity_parallel(iproc, nproc, tmb%linmat%s%nfvctr, tmb%linmat%s%nfvctrp, &
           tmb%linmat%s%isfvctr, ovrlp_fullp, &
           tmb%linmat%s, max_deviation_p, mean_deviation_p)
      max_deviation = max_deviation + max_deviation_p/real(tmb%linmat%s%nspin,kind=8)
      mean_deviation = mean_deviation + mean_deviation_p/real(tmb%linmat%s%nspin,kind=8)
  end do
  call f_free(ovrlp_fullp)
  if (iproc==0) then
      call yaml_map('max dev from unity',max_deviation,fmt='(es9.2)')
      call yaml_map('mean dev from unity',mean_deviation,fmt='(es9.2)')
  end if


      update_phi=.true.
      !update_kernel=.true.
      !not sure if this is overkill...
      call scf_kernel(nit_scc, .false., update_phi)
   end if
end if

  ! Diagonalize the matrix for the FOE/direct min case to get the coefficients. Only necessary if
  ! the Pulay forces are to be calculated, or if we are printing eigenvalues for restart
  if ((input%lin%scf_mode==LINEAR_FOE.or.input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION)& 
       .and. (input%lin%pulay_correction.or.mod(input%lin%plotBasisFunctions,10) /= WF_FORMAT_NONE&
       .or. input%lin%diag_end .or. mod(input%lin%output_coeff_format,10) /= WF_FORMAT_NONE)) then

       !!if (input%lin%scf_mode==LINEAR_FOE) then
       !!    tmb%coeff=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%coeff')
       !!end if

       !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
       call get_coeff(iproc,nproc,LINEAR_MIXDENS_SIMPLE,KSwfn%orbs,at,rxyz,denspot,GPU,&
           infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,update_phi,.true.,.false.,&
           .true.,input%lin%extra_states,itout,0,0,norder_taylor,input%lin%max_inversion_error,&
           input%calculate_KS_residue,input%calculate_gap,energs_work,update_kernel,input%lin%coeff_factor, &
           input%lin%pexsi_npoles,input%lin%pexsi_mumin,input%lin%pexsi_mumax,input%lin%pexsi_mu, &
           input%lin%pexsi_temperature,input%lin%pexsi_tol_charge)
       !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)

       !!if (input%lin%scf_mode==LINEAR_FOE) then
       !!    call f_free_ptr(tmb%coeff)
       !!end if

       if (bigdft_mpi%iproc ==0) then 
          call write_eigenvalues_data(0.1d0,tmb%orbs,mom_vec_fake)
       end if
  end if

  !!if (input%kernel_analysis) then
  !!    call analyze_kernel(iproc, nproc, KSwfn, tmb)
  !!end if

  ! only print eigenvalues if they have meaning, i.e. diag or the case above
  if (input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE.or.input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE) then
     if (bigdft_mpi%iproc ==0) then 
        call write_eigenvalues_data(0.1d0,tmb%orbs,mom_vec_fake)
     end if
  end if


  !! TEST ##########################
  if (input%lin%new_pulay_correction) then
      !!if (input%lin%scf_mode==LINEAR_FOE .and. ) then
      !!    tmb%coeff=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%coeff')
      !!end if
      !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
      call get_coeff(iproc,nproc,LINEAR_MIXDENS_SIMPLE,KSwfn%orbs,at,rxyz,denspot,GPU,&
           infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,update_phi,.true.,.false.,&
           .true.,input%lin%extra_states,itout,0,0,norder_taylor,input%lin%max_inversion_error,&
           input%calculate_KS_residue,input%calculate_gap,energs_work,.false.,input%lin%coeff_factor, &
           input%lin%pexsi_npoles,input%lin%pexsi_mumin,input%lin%pexsi_mumax,input%lin%pexsi_mu,&
           input%lin%pexsi_temperature,input%lin%pexsi_tol_charge)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
      !!call scalprod_on_boundary(iproc, nproc, tmb, kswfn%orbs, at, fpulay)
      call pulay_correction_new(iproc, nproc, tmb, kswfn%orbs, at, fpulay)
      !!if (input%lin%scf_mode==LINEAR_FOE) then
      !!    call f_free_ptr(tmb%coeff)
      !!end if
  end if
  !! END TEST ######################


  ! only do if explicitly activated, but still check for fragment calculation
  if (input%coeff_weight_analysis .and. input%lin%fragment_calculation .and. input%frag%nfrag>1) then
     ! unless we already did a diagonalization, the coeffs will probably be nonsensical in this case, so print a warning
     ! maybe just don't do it in this case?  or do for the whole kernel and not just coeffs?
     if (input%lin%kernel_restart_mode==LIN_RESTART_KERNEL .or.  input%lin%kernel_restart_mode==LIN_RESTART_DIAG_KERNEL) then
        if (iproc==0) call yaml_warning('Output of coeff weight analysis might be nonsensical when restarting from kernel')
     end if
     call coeff_weight_analysis(iproc, nproc, input, KSwfn%orbs, tmb, ref_frags)
  end if


  if (input%lin%constrained_dft) then
     call cdft_data_free(cdft)
     call f_free(coeff_tmp)
  end if


  ! print the final summary
  call print_info(.true.)



  if (iproc==0) call yaml_sequence_close()

  if (input%foe_gap) then
      call calculate_gap_FOE(iproc, nproc, input, KSwfn%orbs, tmb)
  end if



!!  !!!# DIAGONALIZATION TEST ####################################################
!!
!!
!!
!!  !!if (nproc/=1) stop 'only implemented in serial'
!!
!!  matrixElements = f_malloc((/ tmb%linmat%m%nfvctr,tmb%linmat%m%nfvctr,2 /),id='matrixElements')
!!  eval = f_malloc(tmb%linmat%l%nfvctr,id='eval')
!!  all_evals = f_malloc((/at%astruct%nat*tmb%linmat%l%nfvctr,2/),id='eval')
!!
!!  theta = f_malloc0((/at%astruct%nat,tmb%orbs%norbp/),id='theta')
!!
!!  ovrlp_full = f_malloc((/tmb%linmat%m%nfvctr,tmb%linmat%m%nfvctr/),id='ovrlp_full')
!!
!!  ovrlp_large = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='ovrlp_large')
!!  ham_large = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='ham_large')
!!  tmpmat_compr = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='tmpmat_compr')
!!  hamtilde_compr = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='hamtilde_compr')
!!  call transform_sparse_matrix(tmb%linmat%s, tmb%linmat%l, tmb%linmat%ovrlp_%matrix_compr, ovrlp_large, 'small_to_large')
!!  call transform_sparse_matrix(tmb%linmat%m, tmb%linmat%l, tmb%linmat%ham_%matrix_compr, ham_large, 'small_to_large')
!!
!!  tmat = f_malloc0((/tmb%linmat%m%nfvctr,tmb%linmat%m%nfvctr/),id='tmat')
!!
!!  projector_compr = sparsematrix_malloc0(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='projector_compr')
!!
!!  com = f_malloc((/3,tmb%orbs%norbp/),id='com')
!!  call supportfunction_centers(at%astruct%nat, rxyz, size(tmb%psi), tmb%psi, tmb%collcom_sr%ndimpsi_c, &
!!       tmb%orbs, tmb%lzd, com)
!!
!!  ! Diagonalize the entire Hamiltonian of the system
!!  coeff = f_malloc((/ tmb%linmat%m%nfvctr,tmb%linmat%m%nfvctr /),id='coeff')
!!  coeff_all = f_malloc((/ tmb%linmat%m%nfvctr,tmb%linmat%m%nfvctr, at%astruct%nat /),id='coeff_all')
!!
!!  !!do ispin=1,tmb%linmat%m%nspin
!!  !!    if (ispin>1) stop 'not implemented for ispin>1'
!!  !!    ishifts = (ispin-1)*tmb%linmat%s%nvctrp_tg
!!  !!    ishiftm = (ispin-1)*tmb%linmat%m%nvctrp_tg
!!  !!    call f_zero(tmb%linmat%m%nfvctr**2, matrixElements(1,1,1))
!!  !!    tempmat = sparsematrix_malloc(tmb%linmat%m, iaction=DENSE_PARALLEL, id='tempmat')
!!  !!    call uncompress_matrix_distributed2(iproc, tmb%linmat%m, DENSE_PARALLEL, &
!!  !!         tmb%linmat%ham_%matrix_compr(ishiftm+1:ishiftm+tmb%linmat%m%nvctrp_tg), tempmat)
!!  !!    if (tmb%linmat%m%nfvctrp>0) then
!!  !!        call vcopy(tmb%linmat%m%nfvctr*tmb%linmat%m%nfvctrp, tempmat(1,1), 1, &
!!  !!             matrixElements(1,tmb%linmat%m%isfvctr+1,1), 1)
!!  !!    end if
!!  !!    call f_free(tempmat)
!!  !!    if (nproc>1) then
!!  !!        call mpiallred(matrixElements(1,1,1), tmb%linmat%m%nfvctr**2, &
!!  !!             mpi_sum, comm=bigdft_mpi%mpi_comm)
!!  !!    end if
!!
!!  !!    call f_zero(tmb%linmat%s%nfvctr**2, matrixElements(1,1,2))
!!  !!    tempmat = sparsematrix_malloc(tmb%linmat%s, iaction=DENSE_PARALLEL, id='tempmat')
!!  !!    call uncompress_matrix_distributed2(iproc, tmb%linmat%s, DENSE_PARALLEL, &
!!  !!         tmb%linmat%ovrlp_%matrix_compr(ishifts+1:), tempmat)
!!  !!    if (tmb%linmat%m%nfvctrp>0) then
!!  !!        call vcopy(tmb%linmat%s%nfvctr*tmb%linmat%s%nfvctrp, tempmat(1,1), 1, &
!!  !!             matrixElements(1,tmb%linmat%s%isfvctr+1,2), 1)
!!  !!    end if
!!  !!    call f_free(tempmat)
!!  !!    if (nproc>1) then
!!  !!        call mpiallred(matrixElements(1,1,2), tmb%linmat%s%nfvctr**2, &
!!  !!             mpi_sum, comm=bigdft_mpi%mpi_comm)
!!  !!    end if
!!  !!    call diagonalizeHamiltonian2(iproc, tmb%linmat%m%nfvctr, &
!!  !!         matrixElements(1,1,1), matrixElements(1,1,2), eval)
!!  !!    call vcopy(tmb%linmat%m%nfvctr*tmb%linmat%m%nfvctr, matrixElements(1,1,1), 1, coeff(1,1), 1)
!!  !!end do
!!
!!
!!
!!
!!  do iat=1,at%astruct%nat
!!
!!      itype = at%astruct%iatype(iat)
!!
!!      !!%%%! Modify the confinement
!!      !!%%%do iorb=1,tmb%orbs%norb
!!      !!%%%    iiat = tmb%orbs%onwhichatom(iorb)
!!      !!%%%    tmb%confdatarr(iorb)%rxyzConf(1:3) = rxyz(1:3,iat)
!!      !!%%%    tmb%confdatarr(iorb)%potorder = 0
!!      !!%%%    !if (tmb%orbs%onwhichatom(iorb)==iiat) then
!!      !!%%%    !    tmb%confdatarr(iorb)%prefac = 0.d-2
!!      !!%%%    !else
!!      !!%%%        r2 = (rxyz(1,iat)-rxyz(1,iiat))**2 + &
!!      !!%%%             (rxyz(2,iat)-rxyz(2,iiat))**2 + &
!!      !!%%%             (rxyz(3,iat)-rxyz(3,iiat))**2
!!      !!%%%        tmb%confdatarr(iorb)%prefac = 1.d3*r2**4
!!      !!%%%    !end if
!!      !!%%%end do
!!
!!      !!%%%! Start the communication
!!      !!%%%call local_potential_dimensions(iproc,tmb%ham_descr%lzd,tmb%orbs,denspot%xc,denspot%dpbox%ngatherarr(0,1))
!!      !!%%%call start_onesided_communication(iproc, nproc, denspot%dpbox%ndims(1), denspot%dpbox%ndims(2), &
!!      !!%%%     max(denspot%dpbox%nscatterarr(:,2),1), denspot%rhov, &
!!      !!%%%     tmb%ham_descr%comgp%nrecvbuf*tmb%ham_descr%comgp%nspin, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, &
!!      !!%%%     tmb%ham_descr%lzd)
!!
!!      !!%%%! Synchronize the mpi_get before starting a new communication
!!      !!%%%call synchronize_onesided_communication(iproc, nproc, tmb%ham_descr%comgp)
!!
!!      !!%%%!!! Start the communication
!!      !!%%%!!call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
!!      !!%%%!!     TRANSPOSE_POST, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd, wt_phi)
!!
!!      !!%%%! Calculate the unconstrained gradient by applying the Hamiltonian.
!!      !!%%%if (tmb%ham_descr%npsidim_orbs > 0)  call f_zero(tmb%ham_descr%npsidim_orbs,tmb%hpsi(1))
!!      !!%%%call small_to_large_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
!!      !!%%%     tmb%orbs, tmb%psi, tmb%ham_descr%psi)
!!
!!      !!%%%! Start the nonblocking transposition (the results will be gathered in
!!      !!%%%! orthoconstraintNonorthogonal)
!!      !!%%%call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
!!      !!%%%     TRANSPOSE_POST, tmb%ham_descr%psi, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, tmb%ham_descr%lzd, &
!!      !!%%%     wt_philarge)
!!
!!      !!%%%call NonLocalHamiltonianApplication(iproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
!!      !!%%%     tmb%ham_descr%lzd,nlpsp,tmb%ham_descr%psi,tmb%hpsi,energs%eproj,tmb%paw)
!!      !!%%%! only kinetic because waiting for communications
!!      !!%%%call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
!!      !!%%%     tmb%ham_descr%lzd,tmb%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,&
!!      !!%%%     & tmb%ham_descr%psi,tmb%hpsi,energs,input%SIC,GPU,3,denspot%xc,&
!!      !!%%%     & pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
!!      !!%%%     & potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
!!      !!%%%call full_local_potential(iproc,nproc,tmb%orbs,tmb%ham_descr%lzd,2,denspot%dpbox,&
!!      !!%%%     & denspot%xc,denspot%rhov,denspot%pot_work,tmb%ham_descr%comgp)
!!      !!%%%! only potential
!!      !!%%%!if (target_function==TARGET_FUNCTION_IS_HYBRID) then
!!      !!%%%!    call vcopy(tmb%ham_descr%npsidim_orbs, tmb%hpsi(1), 1, hpsi_tmp(1), 1)
!!      !!%%%!    call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
!!      !!%%%!         tmb%ham_descr%lzd,tmb%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,&
!!      !!%%%!         & tmb%ham_descr%psi,tmb%hpsi,energs,input%SIC,GPU,2,denspot%xc,&
!!      !!%%%!         & pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
!!      !!%%%!         & potential=denspot%rhov,comgp=tmb%ham_descr%comgp,&
!!      !!%%%!         hpsi_noconf=hpsi_tmp,econf=econf)
!!
!!      !!%%%!    !!if (nproc>1) then
!!      !!%%%!    !!    call mpiallred(econf, 1, mpi_sum, bigdft_mpi%mpi_comm)
!!      !!%%%!    !!end if
!!
!!      !!%%%!else
!!      !!%%%    call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
!!      !!%%%         tmb%ham_descr%lzd,tmb%confdatarr,denspot%dpbox%ngatherarr,&
!!      !!%%%         & denspot%pot_work,tmb%ham_descr%psi,tmb%hpsi,energs,input%SIC,GPU,2,denspot%xc,&
!!      !!%%%         & pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
!!      !!%%%         & potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
!!      !!%%%!end if
!!
!!
!!      !!%%%!!if (target_function==TARGET_FUNCTION_IS_HYBRID .and. iproc==0) then
!!      !!%%%!!    write(*,*) 'econf, econf/tmb%orbs%norb',econf, econf/tmb%orbs%norb
!!      !!%%%!!end if
!!
!!      !!%%%call timing(iproc,'glsynchham2','ON')
!!      !!%%%call SynchronizeHamiltonianApplication(nproc,tmb%ham_descr%npsidim_orbs,tmb%orbs,tmb%ham_descr%lzd,GPU,denspot%xc,tmb%hpsi,&
!!      !!%%%     energs,energs_work)
!!      !!%%%call timing(iproc,'glsynchham2','OF')
!!
!!      !!%%%if (iproc==0) then
!!      !!%%%    call yaml_map('Hamiltonian Applied',.true.)
!!      !!%%%end if
!!
!!      !!%%%!if (iproc==0) write(*,'(a,5es16.6)') 'ekin, eh, epot, eproj, eex', &
!!      !!%%%!              energs%ekin, energs%eh, energs%epot, energs%eproj, energs%exc
!!
!!      !!%%%hpsit_c = f_malloc_ptr(tmb%ham_descr%collcom%ndimind_c,id='hpsit_c')
!!      !!%%%hpsit_f = f_malloc_ptr(7*tmb%ham_descr%collcom%ndimind_f,id='hpsit_f')
!!
!!
!!      !!%%%! Start the communication
!!      !!%%%!!if (target_function==TARGET_FUNCTION_IS_HYBRID) then
!!      !!%%%!!    call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
!!      !!%%%!!         TRANSPOSE_POST, hpsi_tmp, hpsit_c, hpsit_f, tmb%ham_descr%lzd, wt_hphi)
!!      !!%%%!!else
!!      !!%%%    call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
!!      !!%%%         TRANSPOSE_POST, tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd, wt_hphi)
!!      !!%%%!!end if
!!
!!      !!%%%! Gather the data
!!      !!%%%call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
!!      !!%%%     TRANSPOSE_GATHER, tmb%ham_descr%psi, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, tmb%ham_descr%lzd, &
!!      !!%%%     wt_philarge)
!!      !!%%%call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
!!      !!%%%     TRANSPOSE_GATHER, tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd, wt_hphi)
!!
!!      !!%%%call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
!!      !!%%%     tmb%ham_descr%psit_c, hpsit_c, tmb%ham_descr%psit_f, hpsit_f, tmb%linmat%m, tmb%linmat%ham_)
!!
!!      do ispin=1,tmb%linmat%m%nspin
!!               ishifts = (ispin-1)*tmb%linmat%s%nvctrp_tg
!!               ishiftm = (ispin-1)*tmb%linmat%m%nvctrp_tg
!!               call f_zero(tmb%linmat%m%nfvctr**2, matrixElements(1,1,1))
!!               tempmat = sparsematrix_malloc(tmb%linmat%m, iaction=DENSE_PARALLEL, id='tempmat')
!!               call uncompress_matrix_distributed2(iproc, tmb%linmat%m, DENSE_PARALLEL, &
!!                    tmb%linmat%ham_%matrix_compr(ishiftm+1:ishiftm+tmb%linmat%m%nvctrp_tg), tempmat)
!!               if (tmb%linmat%m%nfvctrp>0) then
!!                   call vcopy(tmb%linmat%m%nfvctr*tmb%linmat%m%nfvctrp, tempmat(1,1), 1, &
!!                        matrixElements(1,tmb%linmat%m%isfvctr+1,1), 1)
!!               end if
!!               call f_free(tempmat)
!!               if (nproc>1) then
!!                   call mpiallred(matrixElements(1,1,1), tmb%linmat%m%nfvctr**2, &
!!                        mpi_sum, comm=bigdft_mpi%mpi_comm)
!!               end if
!!
!!               call f_zero(tmb%linmat%s%nfvctr**2, matrixElements(1,1,2))
!!               tempmat = sparsematrix_malloc(tmb%linmat%s, iaction=DENSE_PARALLEL, id='tempmat')
!!               call uncompress_matrix_distributed2(iproc, tmb%linmat%s, DENSE_PARALLEL, &
!!                    tmb%linmat%ovrlp_%matrix_compr(ishifts+1:), tempmat)
!!               if (tmb%linmat%m%nfvctrp>0) then
!!                   call vcopy(tmb%linmat%s%nfvctr*tmb%linmat%s%nfvctrp, tempmat(1,1), 1, &
!!                        matrixElements(1,tmb%linmat%s%isfvctr+1,2), 1)
!!               end if
!!               call f_free(tempmat)
!!               if (nproc>1) then
!!                   call mpiallred(matrixElements(1,1,2), tmb%linmat%s%nfvctr**2, &
!!                        mpi_sum, comm=bigdft_mpi%mpi_comm)
!!               end if
!!
!!               call vcopy(tmb%linmat%m%nfvctr*tmb%linmat%m%nfvctr, matrixElements(1,1,2), 1, ovrlp_full(1,1), 1)
!!               
!!               ! Add the specical confinement
!!               do iorb=1,tmb%linmat%m%nfvctr
!!                   do jorb=1,tmb%linmat%m%nfvctr
!!                        rr(1) = 0.5d0*(com(1,iorb)+com(1,jorb))
!!                        rr(2) = 0.5d0*(com(2,iorb)+com(2,jorb))
!!                        rr(3) = 0.5d0*(com(3,iorb)+com(3,jorb))
!!                        r2 = (rr(1)-rxyz(1,iat))**2 + (rr(2)-rxyz(2,iat))**2 + (rr(3)-rxyz(3,iat))**2
!!                        matrixElements(jorb,iorb,1) = matrixElements(jorb,iorb,1) + 1.d0*r2**2*matrixElements(jorb,iorb,2)
!!                   end do
!!               end do
!!
!!
!!               call diagonalizeHamiltonian2(iproc, tmb%linmat%m%nfvctr, &
!!                    matrixElements(1,1,1), matrixElements(1,1,2), eval)
!!               call vcopy(tmb%linmat%l%nfvctr, eval(1), 1, all_evals((iat-1)*tmb%linmat%l%nfvctr+1,1), 1)
!!               all_evals((iat-1)*tmb%linmat%l%nfvctr+1:iat*tmb%linmat%l%nfvctr,2) = real(iat,kind=8)
!!               call yaml_map('eval',eval)
!!               call vcopy(tmb%linmat%l%nfvctr**2, matrixElements(1,1,1), 1, coeff_all(1,1,iat), 1)
!!
!!               projector_small = f_malloc0((/tmb%linmat%m%nfvctr,tmb%linmat%m%nfvctr/),id='projector_small')
!!
!!               do i=1,tmb%linmat%m%nfvctr
!!                   do j=1,tmb%linmat%m%nfvctr
!!                       if (itype==2) then
!!                           n=1
!!                           do ieval=1,1
!!                               projector_small(j,i) = projector_small(j,i) + matrixElements(j,ieval,1)*matrixElements(i,ieval,1)
!!                           end do
!!                       else if (itype==1) then
!!                           n=4
!!                           do ieval=1,4
!!                               projector_small(j,i) = projector_small(j,i) + matrixElements(j,ieval,1)*matrixElements(i,ieval,1)
!!                           end do
!!                       end if
!!                   end do
!!               end do
!!               do i=1,tmb%linmat%m%nfvctr
!!                   do j=1,tmb%linmat%m%nfvctr
!!                       write(*,*) 'i, j, proj_small', i, j, projector_small(j,i)
!!                   end do
!!               end do
!!
!!               ! Calculate T_ij
!!               do i=1,tmb%linmat%m%nfvctr
!!                   do j=1,tmb%linmat%m%nfvctr
!!                       tt = 0.d0
!!                       do k=1,tmb%linmat%m%nfvctr
!!                           do l=1,tmb%linmat%m%nfvctr
!!                               tt = tt + coeff(k,i)*ovrlp_full(k,l)*matrixElements(l,j,1)
!!                           end do
!!                       end do
!!                       tmat(i,j) = tt
!!                   end do
!!               end do
!!
!!
!!
!!
!!               !!! Test: diagonalize the projector
!!               ham_small = f_malloc((/tmb%linmat%m%nfvctr,tmb%linmat%m%nfvctr/),id='ham_small')
!!               ham_small = projector_small
!!               lwork = 10*tmb%linmat%m%nfvctr
!!               work = f_malloc(lwork,id='work')
!!               call dsyev('v', 'l', tmb%linmat%m%nfvctr, ham_small, tmb%linmat%m%nfvctr, eval, work, lwork, info)
!!               if (info/=0) call f_err_throw('problem in dsyev', err_name='BIGDFT_LINALG_ERROR')
!!               call yaml_map('atom projector '//adjustl(trim(yaml_toa(iat))),eval,fmt='(es12.5)')
!!               call f_free(work)
!!               call f_free(ham_small)
!!
!!               do i=1,tmb%linmat%m%nfvctr
!!                   ii = i !ist + i
!!                   do j=1,tmb%linmat%m%nfvctr
!!                       jj = j !ist + j
!!                       ind=matrixindex_in_compressed(tmb%linmat%l, ii, jj)
!!                       projector_compr(ind) = projector_compr(ind) + projector_small(j,i)
!!                   end do
!!               end do
!!
!!               do i=1,tmb%linmat%m%nfvctr
!!                   write(*,'(a,12es10.2)') 'evec', matrixElements(i,1:tmb%linmat%l%nfvctr,1)
!!               end do
!!               if (iat==1) then
!!                   !O
!!                   do i=1,3
!!                       do j=1,tmb%linmat%m%nfvctr
!!                           theta(iat,j) = theta(iat,j) + abs(matrixElements(j,i,1))
!!                       end do
!!                   end do
!!               else
!!                   !H
!!                   do i=1,1
!!                       do j=1,tmb%linmat%m%nfvctr
!!                           theta(iat,j) = theta(iat,j) + abs(matrixElements(j,i,1))
!!                       end do
!!                   end do
!!               end if
!!
!!               call f_free(projector_small)
!!      end do
!!
!!  end do
!!
!!  write(*,*) 'TMAT'
!!  do i=1,tmb%linmat%m%nfvctr
!!      write(*,'(12es11.3)') tmat(i,1:tmb%linmat%m%nfvctr)
!!  end do
!!
!!  ! Calculate tmat^T * tmat, use matrixElements as workarray
!!  do i=1,tmb%linmat%m%nfvctr
!!      do j=1,tmb%linmat%m%nfvctr
!!          tt = 0.d0
!!          do k=1,tmb%linmat%m%nfvctr
!!              tt = tt + tmat(k,j)*tmat(k,i)
!!          end do
!!          matrixElements(j,i,1) = tt
!!      end do
!!  end do
!!  write(*,*) 'CHECK UNITARITY'
!!  do i=1,tmb%linmat%m%nfvctr
!!      write(*,'(12es11.3)') matrixElements(i,1:tmb%linmat%m%nfvctr,1)
!!  end do
!!
!!
!!  do i=1,size(projector_compr)
!!      write(*,*) 'i, projector_compr(i)', i, projector_compr(i)
!!  end do
!!
!!  !!! Test: diagonalize the projector
!!  ham_small = f_malloc((/tmb%linmat%m%nfvctr,tmb%linmat%m%nfvctr/),id='ham_small')
!!  call uncompress_matrix2(iproc, nproc, tmb%linmat%l, projector_compr, ham_small)
!!  lwork = 10*tmb%linmat%m%nfvctr
!!  work = f_malloc(lwork,id='work')
!!  call dsyev('v', 'l', tmb%linmat%m%nfvctr, ham_small, tmb%linmat%m%nfvctr, eval, work, lwork, info)
!!  if (info/=0) call f_err_throw('problem in dsyev', err_name='BIGDFT_LINALG_ERROR')
!!  call yaml_map('entire projector '//adjustl(trim(yaml_toa(iat))),eval,fmt='(es12.5)')
!!  call f_free(work)
!!  call f_free(ham_small)
!!
!!  ! Calculate K * S * P * S. Use  hamtilde_compr as workarray
!!  call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, &
!!       projector_compr, ovrlp_large, tmpmat_compr)
!!  hamtilde_compr = tmpmat_compr
!!  call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, &
!!       ovrlp_large, hamtilde_compr, tmpmat_compr)
!!  hamtilde_compr = tmpmat_compr
!!  call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, &
!!       tmb%linmat%kernel_%matrix_compr, hamtilde_compr, tmpmat_compr)
!!
!!  ! Calculate the partial traces
!!  ii = 1
!!  do 
!!      iat = tmb%linmat%l%on_which_atom(ii)
!!      itype = at%astruct%iatype(iat)
!!      n = input%lin%norbsPerType(itype)
!!      tt=0.d0
!!      do i=ii,ii+n-1
!!          ind=matrixindex_in_compressed(tmb%linmat%l, i, i)
!!          tt = tt + tmpmat_compr(ind)
!!      end do
!!      write(*,*) 'iat, itype, charge', iat, itype, tt
!!      ii = ii + n
!!      if (ii==tmb%linmat%l%nfvctr+1) exit
!!  end do
!!
!!  do j=1,tmb%linmat%m%nfvctr
!!      tt = 0.d0
!!      do iat=1,at%astruct%nat
!!          tt = tt + theta(iat,j)
!!      end do
!!      tt = 1.d0/tt
!!      call dscal(at%astruct%nat, tt, theta(1,j), 1)
!!      call yaml_map('theta '//adjustl(trim(yaml_toa(j))),theta(:,j))
!!  end do
!!
!!  do i=1,at%astruct%nat*tmb%linmat%l%nfvctr
!!     write(*,*) 'i, all_evals',i,all_evals(i,1),all_evals(i,2)
!!  end do
!!
!!
!!  do i=1,at%astruct%nat*tmb%linmat%l%nfvctr
!!      ! add i-1 since we are only searching in the subarray
!!      ii = minloc(all_evals(i:at%astruct%nat*tmb%linmat%l%nfvctr,1),1) + (i-1)
!!      write(*,*) 'i, ii',i, ii
!!      tt = all_evals(i,1)
!!      all_evals(i,1) = all_evals(ii,1)
!!      all_evals(ii,1) = tt
!!      tt = all_evals(i,2)
!!      all_evals(i,2) = all_evals(ii,2)
!!      all_evals(ii,2) = tt
!!  end do
!!
!!  ! Check this
!!  if (input%nspin==1) then
!!      fac = 0.5d0
!!  else if (input%nspin==2) then
!!      fac = 1.d0
!!  else 
!!      stop 'wrong nspin'
!!  end if
!!  q = 0.d0
!!  do iat=1,at%astruct%nat
!!      itype = at%astruct%iatype(iat)
!!      q = q + ceiling(fac*real(at%nelpsp(itype),kind=8))
!!  end do
!!  iq = nint(q)
!!
!!
!!  ! Determine how many states will be included
!!  !ef = all_evals(KSwfn%orbs%norb,1) !just a initial guess
!!  !ef_low = all_evals(1,1)
!!  !ef_up = all_evals(at%astruct%nat*tmb%linmat%l%nfvctr,1)
!!
!!  !!do
!!  !!    tot_occ = 0.d0
!!  !!    do i=1,at%astruct%nat*tmb%linmat%l%nfvctr
!!  !!        occ = 1.d0/(1.d0+safe_exp( (all_evals(i,1)-ef)*(1.d0/1.d-2) ) ) 
!!  !!        write(*,*) 'i, occ', i, occ
!!  !!        tot_occ = tot_occ + occ
!!  !!    end do
!!  !!    write(*,*) 'ef, tot_occ, low, up',ef, tot_occ,ef_low, ef_up
!!  !!    if (abs(tot_occ-q)<1.d-6) exit
!!  !!    if ((tot_occ-q)>0.d0) then
!!  !!        ef_up = ef
!!  !!        ef = 0.5d0*(ef+ef_low)
!!  !!    else
!!  !!        ef_low = ef
!!  !!        ef = 0.5d0*(ef+ef_up)
!!  !!    end if
!!  !!end do
!!  ef = all_evals(1,1)
!!  do 
!!      ef = ef + 1.d-3
!!      occ = 1.d0/(1.d0+safe_exp( (all_evals(iq,1)-ef)*(1.d0/1.d-2) ) ) 
!!      write(*,*) 'ef, occ', ef, occ
!!      if (abs(occ-1.d0)<1.d-6) exit
!!  end do
!!  call yaml_map('ef',ef)
!!
!!  do i=1,at%astruct%nat*tmb%linmat%l%nfvctr
!!      occ = 1.d0/(1.d0+safe_exp( (all_evals(i,1)-ef)*(1.d0/1.d-2) ) ) 
!!      !occ = 1.d0/(1.d0+safe_exp( (all_evals(i,1)-all_evals(iq,1))*(1.d0/1.d-2) ) ) 
!!      call yaml_map('ordered',(/real(i,kind=8),all_evals(i,1),all_evals(i,2),occ/))
!!  end do
!!
!!
!!  projector_compr = 0.d0
!!  do iat=1,at%astruct%nat
!!
!!      projector_small = f_malloc0((/tmb%linmat%m%nfvctr,tmb%linmat%m%nfvctr/),id='projector_small')
!!
!!      do i=1,tmb%linmat%m%nfvctr
!!          do j=1,tmb%linmat%m%nfvctr
!!               ii = 0
!!               do ieval=1,at%astruct%nat*tmb%linmat%m%nfvctr
!!                   if (nint(all_evals(ieval,2))/=iat) cycle
!!                   ii = ii + 1
!!                   occ = 1.d0/(1.d0+safe_exp( (all_evals(ieval,1)-ef)*(1.d0/1.d-2) ) )
!!                   write(*,*) 'iat, ieval, occ, ii', iat, ieval, occ, ii
!!                   projector_small(j,i) = projector_small(j,i) + occ*coeff_all(j,ii,iat)*coeff_all(i,ii,iat)
!!               end do
!!          end do
!!      end do
!!
!!      do i=1,tmb%linmat%m%nfvctr
!!          do j=1,tmb%linmat%m%nfvctr
!!              write(*,*) 'i, j, proj_small', i, j, projector_small(j,i)
!!          end do
!!      end do
!!
!!
!!      !!! Test: diagonalize the projector
!!      ham_small = f_malloc((/tmb%linmat%m%nfvctr,tmb%linmat%m%nfvctr/),id='ham_small')
!!      ham_small = projector_small
!!      lwork = 10*tmb%linmat%m%nfvctr
!!      work = f_malloc(lwork,id='work')
!!      call dsyev('v', 'l', tmb%linmat%m%nfvctr, ham_small, tmb%linmat%m%nfvctr, eval, work, lwork, info)
!!      if (info/=0) call f_err_throw('problem in dsyev', err_name='BIGDFT_LINALG_ERROR')
!!      call yaml_map('atom projector '//adjustl(trim(yaml_toa(iat))),eval,fmt='(es12.5)')
!!      call f_free(work)
!!      call f_free(ham_small)
!!
!!      do i=1,tmb%linmat%m%nfvctr
!!          ii = i !ist + i
!!          do j=1,tmb%linmat%m%nfvctr
!!              jj = j !ist + j
!!              ind=matrixindex_in_compressed(tmb%linmat%l, ii, jj)
!!              projector_compr(ind) = projector_compr(ind) + projector_small(j,i)
!!          end do
!!      end do
!!
!!      !!do i=1,tmb%linmat%m%nfvctr
!!      !!    write(*,'(a,12es10.2)') 'evec', matrixElements(i,1:tmb%linmat%l%nfvctr,1)
!!      !!end do
!!      !!if (iat==1) then
!!      !!    !O
!!      !!    do i=1,3
!!      !!        do j=1,tmb%linmat%m%nfvctr
!!      !!            theta(iat,j) = theta(iat,j) + abs(matrixElements(j,i,1))
!!      !!        end do
!!      !!    end do
!!      !!else
!!      !!    !H
!!      !!    do i=1,1
!!      !!        do j=1,tmb%linmat%m%nfvctr
!!      !!            theta(iat,j) = theta(iat,j) + abs(matrixElements(j,i,1))
!!      !!        end do
!!      !!    end do
!!      !!end if
!!
!!      call f_free(projector_small)
!!
!!  end do
!!
!!  do i=1,size(projector_compr)
!!      write(*,*) 'i, projector_compr(i)', i, projector_compr(i)
!!  end do
!!
!!
!!  !!! Test: diagonalize the projector
!!  ham_small = f_malloc((/tmb%linmat%m%nfvctr,tmb%linmat%m%nfvctr/),id='ham_small')
!!  call uncompress_matrix2(iproc, nproc, tmb%linmat%l, projector_compr, ham_small)
!!  lwork = 10*tmb%linmat%m%nfvctr
!!  work = f_malloc(lwork,id='work')
!!  call dsyev('v', 'l', tmb%linmat%m%nfvctr, ham_small, tmb%linmat%m%nfvctr, eval, work, lwork, info)
!!  if (info/=0) call f_err_throw('problem in dsyev', err_name='BIGDFT_LINALG_ERROR')
!!  call yaml_map('entire projector '//adjustl(trim(yaml_toa(iat))),eval,fmt='(es12.5)')
!!  call f_free(work)
!!  call f_free(ham_small)
!!
!!  ! Calculate K * S * P * S. Use  hamtilde_compr as workarray
!!  call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, &
!!       projector_compr, ovrlp_large, tmpmat_compr)
!!  hamtilde_compr = tmpmat_compr
!!  call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, &
!!       ovrlp_large, hamtilde_compr, tmpmat_compr)
!!  hamtilde_compr = tmpmat_compr
!!  call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, &
!!       tmb%linmat%kernel_%matrix_compr, hamtilde_compr, tmpmat_compr)
!!
!!  ! Calculate the partial traces
!!  ii = 1
!!  do 
!!      iat = tmb%linmat%l%on_which_atom(ii)
!!      itype = at%astruct%iatype(iat)
!!      n = input%lin%norbsPerType(itype)
!!      tt=0.d0
!!      do i=ii,ii+n-1
!!          ind=matrixindex_in_compressed(tmb%linmat%l, i, i)
!!          tt = tt + tmpmat_compr(ind)
!!      end do
!!      write(*,*) 'iat, itype, charge', iat, itype, tt
!!      ii = ii + n
!!      if (ii==tmb%linmat%l%nfvctr+1) exit
!!  end do
!!  
!!  !call deallocate_work_transpose(wt_philarge)
!!  !call deallocate_work_transpose(wt_hpsinoprecond)
!!  !call deallocate_work_transpose(wt_hphi)
!!  !call deallocate_work_transpose(wt_phi)
!!
!!  !# END DIAGONALIZATION TEST ################################################


  !!# @VERYNEW TEST ##############################################################
  !! Calculate S^1/1 and S^-1/2
  !call overlapPowerGeneral(iproc, nproc, norder_taylor, 2, (/2,-2/), -1, &
  !      imode=1, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
  !      ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=tmb%linmat%ovrlppowers_, &
  !      check_accur=.true., max_error=max_error, mean_error=mean_error)

  !! Calculate S^-1/2 * H S^-1/2
  !ham_large = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='ham_large')
  !tmpmat_compr = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='tmpmat_compr')
  !hamtilde_compr = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='hamtilde_compr')
  !call transform_sparse_matrix(tmb%linmat%m, tmb%linmat%l, tmb%linmat%ham_%matrix_compr, ham_large, 'small_to_large')
  !call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, ham_large, tmb%linmat%ovrlppowers_(2)%matrix_compr, tmpmat_compr)
  !call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, tmb%linmat%ovrlppowers_(2)%matrix_compr, tmpmat_compr, hamtilde_compr)

  !! Diagonalize the small matrices for each atom and calculate a projector on the occupied subspace
  !projector_compr = sparsematrix_malloc0(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='projector_compr')
  !do iat=1,at%astruct%nat
  !    itype = at%astruct%iatype(iat)
  !    n = input%lin%norbsPerType(itype)
  !    ham_small = f_malloc((/n,n/),id='ham_small')
  !    !ovrlp_small = f_malloc((/n,n/),id='ovrlp_small')
  !    projector_small = f_malloc0((/n,n/),id='projector_small')
  !    eval = f_malloc(n,id='eval')
  !    lwork = 10*n
  !    work = f_malloc(lwork,id='work')
  !    ! Determine the first support function belonging to this atom
  !    !write(*,*) 'tmb%linmat%l%on_which_atom',tmb%linmat%l%on_which_atom
  !    do i=1,tmb%linmat%l%nfvctr
  !        if (tmb%linmat%l%on_which_atom(i)==iat) then
  !            ist = i - 1
  !            exit
  !        end if
  !    end do
  !    do i=1,n
  !        ii = ist + i
  !        do j=1,n
  !            jj = ist + j
  !            ind=matrixindex_in_compressed(tmb%linmat%l, ii, jj)
  !            write(*,*) 'iat, i, j, ii, jj, ind, val', iat, i, j, ii, jj, ind, hamtilde_compr(ind)
  !            ham_small(j,i) = hamtilde_compr(ind)
  !            !ind=matrixindex_in_compressed(tmb%linmat%m, ii, jj)
  !            !ham_small(j,i) = tmb%linmat%ham_%matrix_compr(ind)
  !            !ind=matrixindex_in_compressed(tmb%linmat%s, ii, jj)
  !            !ovrlp_small(j,i) = tmb%linmat%ovrlp_%matrix_compr(ind)
  !        end do
  !    end do
  !    do i=1,n
  !        write(*,'(4es15.5)') ham_small(i,1:n)
  !    end do
  !    !call diagonalizeHamiltonian2(iproc, n, ham_small, ovrlp_small, eval)
  !    call dsyev('v', 'l', n, ham_small, n, eval, work, lwork, info)
  !    if (info/=0) call f_err_throw('problem in dsyev', err_name='BIGDFT_LINALG_ERROR')
  !    call yaml_map('atom '//adjustl(trim(yaml_toa(iat))),eval,fmt='(es12.5)')
  !    do i=1,n
  !        write(*,'(4es15.5)') ham_small(i,1:n)
  !    end do
  !    do i=1,n
  !        do j=1,n
  !            if (itype==2) then
  !                do ieval=1,1
  !                    projector_small(j,i) = projector_small(j,i) + ham_small(j,ieval)*ham_small(i,ieval)
  !                end do
  !            else if (itype==1) then
  !                do ieval=1,4
  !                    projector_small(j,i) = projector_small(j,i) + ham_small(j,ieval)*ham_small(i,ieval)
  !                end do
  !            end if
  !        end do
  !    end do

  !    ! Test: diagonalize the projector
  !    ham_small = projector_small
  !    call dsyev('v', 'l', n, ham_small, n, eval, work, lwork, info)
  !    if (info/=0) call f_err_throw('problem in dsyev', err_name='BIGDFT_LINALG_ERROR')
  !    call yaml_map('atom projector '//adjustl(trim(yaml_toa(iat))),eval,fmt='(es12.5)')

  !    do i=1,n
  !        ii = ist + i
  !        do j=1,n
  !            jj = ist + j
  !            ind=matrixindex_in_compressed(tmb%linmat%l, ii, jj)
  !            projector_compr(ind) = projector_small(j,i)
  !        end do
  !    end do
  !    call f_free(projector_small)
  !    call f_free(ham_small)
  !    call f_free(work)
  !    call f_free(eval)
  !end do
  !do i=1,size(projector_compr)
  !    write(*,*) 'i, projector_compr(i)', i, projector_compr(i)
  !end do

  !! Calculate K * S^1/2 * P * S^1/2. Use  hamtilde_compr as workarray
  !call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, &
  !     projector_compr, tmb%linmat%ovrlppowers_(1)%matrix_compr, tmpmat_compr)
  !hamtilde_compr = tmpmat_compr
  !call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, &
  !     tmb%linmat%ovrlppowers_(1)%matrix_compr, hamtilde_compr, tmpmat_compr)
  !hamtilde_compr = tmpmat_compr
  !call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, &
  !     tmb%linmat%kernel_%matrix_compr, hamtilde_compr, tmpmat_compr)

  !! Calculate the partial traces
  !ii = 1
  !do 
  !    iat = tmb%linmat%l%on_which_atom(ii)
  !    itype = at%astruct%iatype(iat)
  !    n = input%lin%norbsPerType(itype)
  !    tt=0.d0
  !    do i=ii,ii+n-1
  !        ind=matrixindex_in_compressed(tmb%linmat%l, i, i)
  !        tt = tt + tmpmat_compr(ind)
  !    end do
  !    write(*,*) 'iat, itype, charge', iat, itype, tt
  !    ii = ii + n
  !    if (ii==tmb%linmat%l%nfvctr+1) exit
  !end do
  !!# @END VERYNEW TEST ##########################################################


  !!# TEST ################################################################
  !update_phi = .false.
  !! Set to zero all term which couple different atoms
  !do i=1,tmb%linmat%s%nvctr
  !    write(*,*) 'i, val', i, tmb%linmat%ovrlp_%matrix_compr(i)
  !end do
  !ovrlp_orig = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%s, id='ovrlp_orig')
  !kernel_orig = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='kernel_orig')
  !inv_full = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='inv_full')
  !inv_cropped = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='inv_cropped')

  !kernel_cropped = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='kernel_cropped')
  !ovrlp_cropped = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%s, id='ovrlp_cropped')
  !ham_cropped = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%m, id='ham_cropped')
  !kernel_foe_cropped = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='kernel_foe_cropped')

  !call f_memcpy(src=tmb%linmat%kernel_%matrix_compr, dest=kernel_orig)
  !call f_memcpy(src=tmb%linmat%ovrlp_%matrix_compr, dest=ovrlp_orig)
  !! Calculate S^1/2 for the original overlap matrix
  !call overlapPowerGeneral(iproc, nproc, norder_taylor, 1, (/2/), -1, &
  !      imode=1, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
  !      ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=tmb%linmat%ovrlppowers_(2), &
  !      check_accur=.true., max_error=max_error, mean_error=mean_error)
  !inv_full = tmb%linmat%ovrlppowers_(2)%matrix_compr

  !call delete_coupling_terms(iproc, nproc, tmb%linmat%s, tmb%linmat%ovrlp_%matrix_compr)
  !call delete_coupling_terms(iproc, nproc, tmb%linmat%m, tmb%linmat%ham_%matrix_compr)
  !call delete_coupling_terms(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_%matrix_compr)

  !ovrlp_large = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='ovrlp_large')
  !tmpmat = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='tmpmat')
  !call transform_sparse_matrix(tmb%linmat%s, tmb%linmat%l, tmb%linmat%ovrlp_%matrix_compr, ovrlp_large, 'small_to_large')

  !ovrlp_cropped = tmb%linmat%ovrlp_%matrix_compr
  !ham_cropped = tmb%linmat%ham_%matrix_compr
  !kernel_cropped = tmb%linmat%kernel_%matrix_compr

  !! Calculate S^1/2 for the cropped overlap matrix
  !call overlapPowerGeneral(iproc, nproc, norder_taylor, 1, (/2/), -1, &
  !      imode=1, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
  !      ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=tmb%linmat%ovrlppowers_(1), &
  !      check_accur=.true., max_error=max_error, mean_error=mean_error)
  !inv_cropped = tmb%linmat%ovrlppowers_(1)%matrix_compr

  !kernel_foe_cropped = 0.d0
  !ii = 1
  !do 
  !    tmb%linmat%ovrlp_%matrix_compr = 0.d0
  !    tmb%linmat%ham_%matrix_compr = 0.d0
  !    tmb%linmat%kernel_%matrix_compr = 0.d0
  !    ! Overlap must have 1 one diagonal, for inversion
  !    do i=1,tmb%linmat%l%nfvctr
  !        ind=matrixindex_in_compressed(tmb%linmat%l, i, i)
  !        tmb%linmat%ovrlp_%matrix_compr(ind) = 1.d0
  !    end do
  !    iat = tmb%linmat%l%on_which_atom(ii)
  !    itype = at%astruct%iatype(iat)
  !    n = input%lin%norbsPerType(itype)
  !    do i=ii,ii+n-1
  !        do j=ii,ii+n-1
  !            ind=matrixindex_in_compressed(tmb%linmat%l, i, j)
  !            tmb%linmat%ovrlp_%matrix_compr(ind) = ovrlp_cropped(ind)
  !            tmb%linmat%ham_%matrix_compr(ind) = ham_cropped(ind)
  !            tmb%linmat%kernel_%matrix_compr(ind) = kernel_cropped(ind)
  !        end do
  !    end do
  !    if (n==1) then
  !        call foe_data_set_real(tmb%foe_obj,"charge",1.d0,1)
  !    else if (n==4) then
  !        call foe_data_set_real(tmb%foe_obj,"charge",6.d0,1)
  !    end if

  !    ! HERE FOE #########################################
  !    call fermi_operator_expansion(iproc, nproc, 0.d0, &
  !         energs%ebs, norder_taylor, 1.d-8, .true., &
  !         .true., 2, FOE_ACCURATE, &
  !         trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//'-'//trim(adjustl(yaml_toa(-999,fmt='(i3.3)'))), &
  !         tmb, tmb%linmat%ham_, tmb%linmat%ovrlp_, tmb%linmat%kernel_, tmb%foe_obj)
  !    write(*,*) 'tmb%linmat%kernel_%matrix_compr(1)',tmb%linmat%kernel_%matrix_compr(1)
  !    write(*,*) 'tmb%linmat%ham_%matrix_compr(1)',tmb%linmat%ham_%matrix_compr(1)
  !    write(*,*) 'tmb%linmat%ovrlp_%matrix_compr(1)',tmb%linmat%ovrlp_%matrix_compr(1)
  !    do i=ii,ii+n-1
  !        do j=ii,ii+n-1
  !            ind=matrixindex_in_compressed(tmb%linmat%l, i, j)
  !            if (n==1) then
  !                kernel_foe_cropped(ind) = tmb%linmat%kernel_%matrix_compr(ind)
  !            else if (n==4) then
  !                kernel_foe_cropped(ind) = 0.5d0*tmb%linmat%kernel_%matrix_compr(ind)
  !            end if
  !        end do
  !    end do
  !    ii = ii + n
  !    if (ii==tmb%linmat%l%nfvctr+1) exit
  !    !! HERE DIAG #########################################


  !    !! END OPTIONS #######################################
  !end do
  !tmb%linmat%kernel_%matrix_compr = kernel_foe_cropped



  !write(*,*) '==='
  !do i=1,tmb%linmat%s%nvctr
  !    write(*,*) 'i, val', i, tmb%linmat%ovrlp_%matrix_compr(i)
  !end do
  !!call scf_kernel(1, .true., update_phi)
  !do i=1,tmb%linmat%l%nvctr
  !    write(*,*) 'i, val ovrlp_large', i, ovrlp_large(i)
  !    write(*,*) 'i, val inv full', i, inv_full(i)
  !    write(*,*) 'i, val inv cropped', i, inv_cropped(i)
  !    write(*,*) 'i, val kernel', i, tmb%linmat%kernel_%matrix_compr(i)
  !    write(*,*) 'i, val kernel orig', i, kernel_orig(i)
  !end do
  !!tmb%linmat%kernel_%matrix_compr = 0.5d0*tmb%linmat%kernel_%matrix_compr
  !! use ovrlp_large as workarray
  !! NEW $$$$$$$$$$$$$$$$$
  !ovrlp_large = inv_cropped
  !call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_%matrix_compr, ovrlp_large, tmpmat)
  !ovrlp_large = tmpmat
  !call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, inv_cropped, ovrlp_large, tmpmat)
  !ovrlp_large = tmpmat
  !call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, inv_full, ovrlp_large, tmpmat)
  !ovrlp_large = tmpmat
  !call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, kernel_orig, ovrlp_large, tmpmat)
  !ovrlp_large = tmpmat
  !call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, inv_full, ovrlp_large, tmpmat)
  !!! OLD $$$$$$$$$$$$$$$$$
  !!ovrlp_large = inv_full
  !!call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, 0.5d0*kernel_orig, ovrlp_large, tmpmat)
  !!ovrlp_large = tmpmat
  !!call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, inv_full, ovrlp_large, tmpmat)
  !!ovrlp_large = tmpmat
  !!call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, inv_full, ovrlp_large, tmpmat)
  !!ovrlp_large = tmpmat
  !!call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, kernel_orig, ovrlp_large, tmpmat)
  !!ovrlp_large = tmpmat
  !!call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, inv_full, ovrlp_large, tmpmat)

  !!! ALTERNATIVE VERSION
  !!ovrlp_large = kernel_foe_cropped
  !!call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, inv_cropped, ovrlp_large, tmpmat)
  !!ovrlp_large = tmpmat
  !!call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, inv_full, ovrlp_large, tmpmat)
  !!ovrlp_large = tmpmat
  !!call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, kernel_orig, ovrlp_large, tmpmat)
  !!tmb%linmat%kernel_%matrix_compr = tmpmat
  !!overlap_calculated = .true.
  !!call foe_data_set_real(tmb%foe_obj,"charge",8.d0,1)
  !!tmb%linmat%ovrlp_%matrix_compr = ovrlp_orig
  !!call purify_kernel(iproc, nproc, tmb, overlap_calculated, 20, 20, norder_taylor, &
  !!     1.d-8, .false., 1)
  !!ovrlp_large = inv_full
  !!call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_%matrix_compr, ovrlp_large, tmpmat)
  !!ovrlp_large = tmpmat
  !!call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, inv_full, ovrlp_large, tmpmat)



  !ii = 1
  !do 
  !    iat = tmb%linmat%l%on_which_atom(ii)
  !    itype = at%astruct%iatype(iat)
  !    n = input%lin%norbsPerType(itype)
  !    tt=0.d0
  !    do i=ii,ii+n-1
  !        ind=matrixindex_in_compressed(tmb%linmat%l, i, i)
  !        tt = tt + tmpmat(ind)
  !    end do
  !    write(*,*) 'iat, itype, charge', iat, itype, tt
  !    ii = ii + n
  !    if (ii==tmb%linmat%l%nfvctr+1) exit
  !end do

  !!!call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, inv_ovrlp_mat=tmb%linmat%ovrlppowers_(2)%matrix_compr(i), kernel_orig, tmpmat)
  !!!call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, inv_ovrlp_mat=tmb%linmat%ovrlppowers_(2)%matrix_compr(i), kernel_orig, tmpmat)
  !!!call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, kernel_orig, tmpmat, tmb%linmat%kernel_%matrix_compr)
  !!# END TEST ############################################################


  !!# TEST 2 ################################################################
  !update_phi = .false.
  !! Set to zero all term which couple different atoms
  !do i=1,tmb%linmat%s%nvctr
  !    write(*,*) 'i, val', i, tmb%linmat%ovrlp_%matrix_compr(i)
  !end do
  !kernel_orig = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='kernel_orig')
  !call f_memcpy(src=tmb%linmat%kernel_%matrix_compr, dest=kernel_orig)
  !!call delete_coupling_terms(iproc, nproc, tmb%linmat%s, tmb%linmat%ovrlp_%matrix_compr)
  !!call delete_coupling_terms(iproc, nproc, tmb%linmat%m, tmb%linmat%ham_%matrix_compr)
  !!call delete_coupling_terms(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_%matrix_compr)
  !!write(*,*) '==='
  !!do i=1,tmb%linmat%s%nvctr
  !!    write(*,*) 'i, val', i, tmb%linmat%ovrlp_%matrix_compr(i)
  !!end do
  !!call scf_kernel(20, .true., update_phi)
  !K_H(1,1) =  9.9999999965251773E-01
  !K_O(1,1) =  2.0000051616997556E+00
  !K_O(2,1) = -1.7379613199227223E-03
  !K_O(3,1) = -1.7378989154477071E-03
  !K_O(4,1) = -1.7377247360473144E-03
  !K_O(1,2) = -1.7379613199235630E-03
  !K_O(2,2) =  1.3333348599588517E+00
  !K_O(3,2) = -6.1859793911177513E-05
  !K_O(4,2) = -6.1895423401417824E-05
  !K_O(1,3) = -1.7378989154472383E-03
  !K_O(2,3) = -6.1859793902466261E-05
  !K_O(3,3) =  1.3333348570644272E+00
  !K_O(4,3) = -6.1914916440734041E-05
  !K_O(1,4) = -1.7377247360467066E-03
  !K_O(2,4) = -6.1895423367217710E-05
  !K_O(3,4) = -6.1914916487668014E-05
  !K_O(4,4) =  1.3333348527961533E+00

  !!call f_zero(tmb%linmat%kernel_%matrix_compr)
  !!ii = 1
  !!do 
  !!    iat = tmb%linmat%l%on_which_atom(ii)
  !!    itype = at%astruct%iatype(iat)
  !!    n = input%lin%norbsPerType(itype)
  !!    do i=ii,ii+n-1
  !!        do j=ii,ii+n-1
  !!            ind=matrixindex_in_compressed(tmb%linmat%l, i, j)
  !!            if (n==1) then
  !!                tmb%linmat%kernel_%matrix_compr(ind) = K_H(i-ii+1,j-ii+1)
  !!            else if (n==4) then
  !!                tmb%linmat%kernel_%matrix_compr(ind) = K_O(i-ii+1,j-ii+1)
  !!            end if
  !!        end do 
  !!    end do
  !!    ii = ii + n
  !!    if (ii==tmb%linmat%l%nfvctr+1) exit
  !!end do
  !tmb%linmat%kernel_%matrix_compr = 0.5d0*tmb%linmat%kernel_%matrix_compr
  !ovrlp_large = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='ovrlp_large')
  !tmpmat = sparsematrix_malloc(iaction=SPARSE_TASKGROUP, smat=tmb%linmat%l, id='tmpmat')
  !call transform_sparse_matrix(tmb%linmat%s, tmb%linmat%l, tmb%linmat%ovrlp_%matrix_compr, ovrlp_large, 'small_to_large')
  !do i=1,tmb%linmat%l%nvctr
  !    write(*,*) 'i, val ovrlp_large', i, ovrlp_large(i)
  !    write(*,*) 'i, val kernel', i, tmb%linmat%kernel_%matrix_compr(i)
  !    write(*,*) 'i, val kernel orig', i, kernel_orig(i)
  !end do
  !call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, ovrlp_large, tmb%linmat%kernel_%matrix_compr, tmpmat)
  !call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, kernel_orig, tmpmat, tmb%linmat%kernel_%matrix_compr)

  !call matrix_matrix_mult_wrapper(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_%matrix_compr, ovrlp_large, tmpmat)
  !ii = 1
  !do 
  !    iat = tmb%linmat%l%on_which_atom(ii)
  !    itype = at%astruct%iatype(iat)
  !    n = input%lin%norbsPerType(itype)
  !    tt=0.d0
  !    do i=ii,ii+n-1
  !        ind=matrixindex_in_compressed(tmb%linmat%l, i, i)
  !        tt = tt + tmpmat(ind)
  !    end do
  !    write(*,*) 'iat, itype, charge', iat, itype, tt
  !    ii = ii + n
  !    if (ii==tmb%linmat%l%nfvctr+1) exit
  !end do

  !!# END TEST 2 ############################################################

  if (input%loewdin_charge_analysis) then
      if (iproc==0) then
          call yaml_mapping_open('Charge analysis, projector approach')
      end if
      !!call projector_for_charge_analysis(at, tmb%linmat%s, tmb%linmat%m, tmb%linmat%l, &
      !!     tmb%linmat%ovrlp_, tmb%linmat%ham_, tmb%linmat%kernel_, tmb%lzd, &
      !!     rxyz, input%lin%norbsPerType, centers_provided=.false., &
      !!     nphirdim=tmb%collcom_sr%ndimpsi_c, psi=tmb%psi, orbs=tmb%orbs)
      !com = f_malloc_ptr((/3,tmb%linmat%s%nfvctr/),id='com')
      !do i=1,tmb%linmat%s%nfvctr
      !    iat = tmb%linmat%s%on_which_atom(i)
      !    com(1:3,i) = rxyz(1:3,iat)
      !end do

      !!! TEST ##########################################################
      !!multipoles_out = f_malloc((/-lmax.to.lmax,0.to.lmax,1.to.at%astruct%nat/),id='multipoles_out')
      !!call multipoles_from_density(iproc, nproc, at, tmb%lzd, tmb%linmat%s, tmb%linmat%l, tmb%orbs, &
      !!     tmb%npsidim_orbs, tmb%psi, input%lin%norbsPerType, tmb%collcom, tmb%collcom_sr, tmb%orthpar, &
      !!     tmb%linmat%ovrlp_, tmb%linmat%kernel_, meth_overlap=norder_taylor, &
      !!     multipoles_out=multipoles_out)
      !!! END TEST ######################################################

      !! Calculate the support function multipoles
      !multipoles = f_malloc0((/-lmax.to.lmax,0.to.lmax,1.to.tmb%orbs%norb/),id='multipoles')
      !call support_function_gross_multipoles(iproc, tmb, at, denspot, multipoles)

      !write(*,*) 'call with multipoles'
      !call projector_for_charge_analysis(at, tmb%linmat%s, tmb%linmat%m, tmb%linmat%l, &
      !     tmb%linmat%ovrlp_, tmb%linmat%ham_, tmb%linmat%kernel_, &
      !     rxyz, calculate_centers=.false., multipoles=multipoles)



      ! @ NEW ##################################################################################################
      ! Calculate the matrices <phi|r**x|phi>
      do i=1,24
          rpower_matrix(i) = matrices_null()
          rpower_matrix(i)%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%s, SPARSE_FULL, id='rpower_matrix(i)%matrix_compr')
      end do
      call calculate_rpowerx_matrices(iproc, nproc, tmb%npsidim_orbs, tmb%collcom_sr%ndimpsi_c, tmb%lzd, &
           tmb%orbs, tmb%collcom, tmb%psi, tmb%linmat%s, rpower_matrix)
      ! @ END NEW ##############################################################################################
      call projector_for_charge_analysis(at, tmb%linmat%s, tmb%linmat%m, tmb%linmat%l, &
           tmb%linmat%ovrlp_, tmb%linmat%ham_, tmb%linmat%kernel_, &
           rxyz, calculate_centers=.false., write_output=.true., ortho='yes', &
           rpower_matrix=rpower_matrix, orbs=tmb%orbs)
      do i=1,24
          call deallocate_matrices(rpower_matrix(i))
      end do
      !call f_free(multipoles)
      !call f_free(multipoles_out)

      !call f_free_ptr(com)
      if (iproc==0) then
          call yaml_mapping_close()
      end if

      !call loewdin_charge_analysis(iproc, tmb, at, denspot, calculate_overlap_matrix=.true., &
      !     calculate_ovrlp_half=.true., meth_overlap=0)
      !theta = f_malloc((/at%astruct%nat,tmb%orbs%norbp/),id='theta')
      !!! This check is here to prevent inconsictencies between orbs distribution and matrix distribution, to be fixed
      !!if (tmb%orbs%norbp/=tmb%linmat%s%nfvctrp) then
      !!    call f_err_throw('tmb%orbs%norbp/=tmb%linmat%s%nfvctrp',err_name='BIGDFT_RUNTIME_ERROR')
      !!end if
      !call calculate_theta(at%astruct%nat, rxyz, size(tmb%psi), tmb%psi, tmb%collcom_sr%ndimpsi_c, &
      !     tmb%orbs, tmb%lzd, theta)
      !write(*,*) 'theta',theta
      if (iproc==0) then
          call yaml_mapping_open('Charge analysis, Loewdin approach')
      end if
      call loewdin_charge_analysis(iproc, tmb, at, denspot, calculate_overlap_matrix=.true., &
           calculate_ovrlp_half=.true., meth_overlap=norder_taylor, blocksize=tmb%orthpar%blocksize_pdsyev)!, &
           !ntheta=tmb%orbs%norbp, istheta=tmb%orbs%isorb, theta=theta)
      if (iproc==0) then
          call yaml_mapping_close()
      end if
      !!call support_function_multipoles(iproc, tmb, at, denspot)

      ! THIS IS COMMENTED FOR THE MOMENT #############################################################3

  end if

  if (input%support_function_multipoles) then
      call support_function_gross_multipoles(iproc, nproc, tmb, at, shift, denspot)
  end if

  if (input%lin%charge_multipoles>0) then
    !!$$ UNCOMMENT FOR TEST  !!write(200+iproc,*) tmb%linmat%ovrlp_%matrix_compr
    !!$$ UNCOMMENT FOR TEST  !!write(210+iproc,*) tmb%linmat%kernel_%matrix_compr

    !!$$ UNCOMMENT FOR TEST  ! TEST ################################################
    !!$$ UNCOMMENT FOR TEST  call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
    !!$$ UNCOMMENT FOR TEST       tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, &
    !!$$ UNCOMMENT FOR TEST       denspot%dpbox%ndimrhopot, &
    !!$$ UNCOMMENT FOR TEST       denspot%rhov, rho_negative)
    !!$$ UNCOMMENT FOR TEST  if (rho_negative) then
    !!$$ UNCOMMENT FOR TEST      call corrections_for_negative_charge(iproc, nproc, at, denspot)
    !!$$ UNCOMMENT FOR TEST  end if
    !!$$ UNCOMMENT FOR TEST  is3 = denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+1
    !!$$ UNCOMMENT FOR TEST  ie3 = denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,2)
    !!$$ UNCOMMENT FOR TEST  is2 = 1
    !!$$ UNCOMMENT FOR TEST  ie2 = denspot%dpbox%ndims(2)
    !!$$ UNCOMMENT FOR TEST  is1 = 1
    !!$$ UNCOMMENT FOR TEST  ie1 = denspot%dpbox%ndims(1)
    !!$$ UNCOMMENT FOR TEST  ii = 0
    !!$$ UNCOMMENT FOR TEST  do i3=is3,ie3
    !!$$ UNCOMMENT FOR TEST      do i2=is2,ie2
    !!$$ UNCOMMENT FOR TEST          do i1=is1,ie1
    !!$$ UNCOMMENT FOR TEST              ii = ii + 1
    !!$$ UNCOMMENT FOR TEST              write(190+iproc,'(3(a,i6),a,es18.8)') 'i1= ',i1,' i2= ',i2,' i3= ',i3,' val= ',denspot%rhov(ii)
    !!$$ UNCOMMENT FOR TEST          end do
    !!$$ UNCOMMENT FOR TEST      end do
    !!$$ UNCOMMENT FOR TEST  end do

    !!$$ UNCOMMENT FOR TEST  !write(*,*) 'BEFORE: sum(rhov)',sum(denspot%rhov)
    !!$$ UNCOMMENT FOR TEST  !write(*,*) 'BEFORE: sum(V_ext)',sum(denspot%V_ext)
    !!$$ UNCOMMENT FOR TEST  call H_potential('D',denspot%pkernel,denspot%rhov,denspot%V_ext,ehart_ps,0.0_dp,.true.,&
    !!$$ UNCOMMENT FOR TEST       quiet=denspot%PSquiet)!,rho_ion=denspot%rho_ion)
    !!$$ UNCOMMENT FOR TEST  !write(*,*) 'AFTER: sum(rhov)',sum(denspot%rhov)
    !!$$ UNCOMMENT FOR TEST  is3 = denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+1
    !!$$ UNCOMMENT FOR TEST  ie3 = denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,2)
    !!$$ UNCOMMENT FOR TEST  is2 = 1
    !!$$ UNCOMMENT FOR TEST  ie2 = denspot%dpbox%ndims(2)
    !!$$ UNCOMMENT FOR TEST  is1 = 1
    !!$$ UNCOMMENT FOR TEST  ie1 = denspot%dpbox%ndims(1)
    !!$$ UNCOMMENT FOR TEST  ii = 0
    !!$$ UNCOMMENT FOR TEST  do i3=is3,ie3
    !!$$ UNCOMMENT FOR TEST      do i2=is2,ie2
    !!$$ UNCOMMENT FOR TEST          do i1=is1,ie1
    !!$$ UNCOMMENT FOR TEST              ii = ii + 1
    !!$$ UNCOMMENT FOR TEST              write(200+iproc,'(3(a,i6),a,es18.8)') 'i1= ',i1,' i2= ',i2,' i3= ',i3,' val= ',denspot%rhov(ii)
    !!$$ UNCOMMENT FOR TEST          end do
    !!$$ UNCOMMENT FOR TEST      end do
    !!$$ UNCOMMENT FOR TEST  end do


      !if (input%lin%charge_multipoles==1) then
      !    call multipoles_from_density(iproc, nproc, at, tmb%lzd, tmb%linmat%s, tmb%linmat%l, tmb%orbs, &
      !         tmb%npsidim_orbs, tmb%psi, input%lin%norbsPerType, tmb%collcom, tmb%collcom_sr, tmb%orthpar, &
      !         tmb%linmat%ovrlp_, tmb%linmat%kernel_, meth_overlap=norder_taylor)
      !    !write(300+iproc,*) tmb%linmat%ovrlp_%matrix_compr
      !    !write(310+iproc,*) tmb%linmat%kernel_%matrix_compr
      multipoles = f_malloc_ptr((/-lmax.to.lmax,0.to.lmax,1.to.at%astruct%nat/),id='multipoles')
      if (input%lin%charge_multipoles/=0) then
          select case (input%lin%charge_multipoles)
          case (1,11)
              method='loewdin'
          case (2,12) 
              method='projector'
          case default
              call f_err_throw('wrong value of charge_multipoles')
          end select
          select case (input%lin%charge_multipoles)
          case (1,2)
              do_ortho='yes'
          case (11,12) 
              do_ortho='no'
          case default
              call f_err_throw('wrong value of charge_multipoles')
          end select
      end if
      !if (input%lin%charge_multipoles==1) then
          call multipole_analysis_driver(iproc, nproc, lmax, tmb%npsidim_orbs, tmb%psi, &
               max(tmb%collcom_sr%ndimpsi_c,1), at, tmb%lzd%hgrids, &
               tmb%orbs, tmb%linmat%s, tmb%linmat%m, tmb%linmat%l, tmb%collcom, tmb%collcom_sr, tmb%lzd, &
               denspot, tmb%orthpar, tmb%linmat%ovrlp_, tmb%linmat%ham_, tmb%linmat%kernel_, rxyz, &
               method=method, do_ortho=do_ortho, shift=shift, nsigma=input%nsigma, ixc=input%ixc, ep=ep )
      !!else if (input%lin%charge_multipoles==2) then
      !!    call multipole_analysis_driver(iproc, nproc, lmax, tmb%npsidim_orbs, tmb%psi, &
      !!         max(tmb%collcom_sr%ndimpsi_c,1), at, tmb%lzd%hgrids, &
      !!         tmb%orbs, tmb%linmat%s, tmb%linmat%m, tmb%linmat%l, tmb%collcom, tmb%collcom_sr, tmb%lzd, &
      !!         denspot, tmb%orthpar, tmb%linmat%ovrlp_, tmb%linmat%ham_, tmb%linmat%kernel_, rxyz, &
      !!         method='projector', do_ortho=do_ortho, shift=shift, nsigma=input%nsigma, ixc=input%ixc, ep=ep)
      !!end if
      !call get_optimal_sigmas(iproc, nproc, KSwfn, tmb, at, input, ep, shift, denspot)
      !!# TEST ######################################################################################################
      !test_pot = f_malloc0((/size(denspot%V_ext,1),size(denspot%V_ext,2),size(denspot%V_ext,3),2/),id='test_pot')
      !call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
      !     tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, &
      !     denspot%dpbox%ndimrhopot, &
      !     denspot%rhov, rho_negative)
      !if (rho_negative) then
      !    call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
      !end if
      !call H_potential('D',denspot%pkernel,denspot%rhov,denspot%V_ext,ehart_ps,0.0_dp,.true.,&
      !     quiet=denspot%PSquiet)!,rho_ion=denspot%rho_ion)
      !call dcopy(size(denspot%V_ext,1)*size(denspot%V_ext,2)*size(denspot%V_ext,3), &
      !     denspot%rhov(1), 1, test_pot(1,1,1,1), 1)
      !!call axpy(size(denspot%V_ext,1)*size(denspot%V_ext,2)*size(denspot%V_ext,3), &
      !!     1.d0, denspot%V_ext(1,1,1,1), 1, test_pot(1,1,1,1), 1)
      !call dcopy(size(denspot%V_ext,1)*size(denspot%V_ext,2)*size(denspot%V_ext,3), &
      !     denspot%rhov(1), 1, test_pot(1,1,1,2), 1)
      !call potential_from_charge_multipoles(iproc, nproc, at, denspot, ep, 1, &
      !     denspot%dpbox%ndims(1), 1, denspot%dpbox%ndims(2), &
      !     denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+1, &
      !     denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,2), &
      !     denspot%dpbox%hgrids(1),denspot%dpbox%hgrids(2),denspot%dpbox%hgrids(3), shift, test_pot(:,:,:,2))
      !do i3=1,size(denspot%V_ext,3)
      !    do i2=1,size(denspot%V_ext,2)
      !        do i1=1,size(denspot%V_ext,1)
      !            write(800,*) 'i1, i2, i3, vals', i1, i2, i3, test_pot(i1,i2,i3,1), test_pot(i1,i2,i3,2)
      !        end do
      !    end do
      !end do
      !call f_free(test_pot)
      !!# TEST ######################################################################################################
      call f_free_ptr(multipoles)
      call deallocate_external_potential_descriptors(ep)
  end if


  ! Deallocate everything that is not needed any more.
  if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) call DIIS_free(ldiis_coeff)!call deallocateDIIS(ldiis_coeff)
  call deallocateDIIS(ldiis)
  !!call wait_p2p_communication(iproc, nproc, tmb%comgp)
  call synchronize_onesided_communication(iproc, nproc, tmb%comgp)
  call deallocate_p2pComms_buffer(tmb%comgp)


  if (input%lin%pulay_correction .and. .not.input%lin%new_pulay_correction) then
      !if (iproc==0) write(*,'(1x,a)') 'WARNING: commented correction_locrad!'
      if (iproc==0) call yaml_warning('commented correction_locrad')
      !!! Testing energy corrections due to locrad
      !!call correction_locrad(iproc, nproc, tmblarge, KSwfn%orbs,tmb%coeff) 
      ! Calculate Pulay correction to the forces
      !!if (input%lin%scf_mode==LINEAR_FOE) then
      !!    tmb%coeff=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%coeff')
      !!end if
      call pulay_correction(iproc, nproc, KSwfn%orbs, at, rxyz, nlpsp, input%SIC, denspot, GPU, tmb, fpulay)
      !!if (input%lin%scf_mode==LINEAR_FOE) then
      !!    call f_free_ptr(tmb%coeff)
      !!end if
  else if (.not.input%lin%new_pulay_correction) then
      call f_zero(fpulay)
  end if

  if(tmb%ham_descr%can_use_transposed) then
      call f_free_ptr(tmb%ham_descr%psit_c)
      call f_free_ptr(tmb%ham_descr%psit_f)
      tmb%ham_descr%can_use_transposed=.false.
  end if
  ! here or cluster, not sure which is best
  deallocate(tmb%confdatarr, stat=istat)


  !Write the linear wavefunctions to file if asked
  if (mod(input%lin%plotBasisFunctions,10) /= WF_FORMAT_NONE) then
     nelec=0
     do iat=1,at%astruct%nat
        ityp=at%astruct%iatype(iat)
        nelec=nelec+at%nelpsp(ityp)
     enddo
     if (write_full_system) then
        call writemywaves_linear(iproc,trim(input%dir_output) // 'minBasis',mod(input%lin%plotBasisFunctions,10),&
             max(tmb%npsidim_orbs,tmb%npsidim_comp),tmb%Lzd,tmb%orbs,nelec,at,rxyz,tmb%psi,tmb%linmat%l%nfvctr,tmb%coeff)
        !!call write_linear_matrices(iproc,nproc,input%imethod_overlap,trim(input%dir_output),&
        !!     mod(input%lin%plotBasisFunctions,10),tmb,at,rxyz,input%lin%calculate_onsite_overlap)
        !!call write_linear_coefficients(0, trim(input%dir_output)//'KS_coeffs.bin', at, rxyz, &
        !!     tmb%linmat%l%nfvctr, tmb%orbs%norb, tmb%linmat%l%nspin, tmb%coeff, tmb%orbs%eval)
        if (input%lin%plotBasisFunctions>20) then
            call write_orbital_density(iproc, .true., mod(input%lin%plotBasisFunctions,10), 'SupFun', &
                 tmb%npsidim_orbs, tmb%psi, input, tmb%orbs, KSwfn%lzd, at, rxyz, .false., tmb%lzd)
        else if (input%lin%plotBasisFunctions>10) then
            call write_orbital_density(iproc, .true., mod(input%lin%plotBasisFunctions,10), 'SupFunDens', &
                 tmb%npsidim_orbs, tmb%psi, input, tmb%orbs, KSwfn%lzd, at, rxyz, .true., tmb%lzd)

        end if
     end if
     !write as fragments - for now don't write matrices, think later if this is useful/worth the effort
     !(now kernel is done internally in writemywaves)
     if (write_fragments .and. input%lin%fragment_calculation) then
        call writemywaves_linear_fragments(iproc,'minBasis',mod(input%lin%plotBasisFunctions,10),&
             max(tmb%npsidim_orbs,tmb%npsidim_comp),tmb%Lzd,tmb%orbs,nelec,at,rxyz,tmb%psi,tmb%coeff, &
             trim(input%dir_output),input%frag,ref_frags,tmb%linmat,norder_taylor,input%lin%max_inversion_error,&
             tmb%orthpar,input%lin%frag_num_neighbours,input%lin%frag_neighbour_cutoff)

!      call orthonormalizeLocalized(iproc, nproc, norder_taylor, input%lin%max_inversion_error, tmb%npsidim_orbs, tmb%orbs, tmb%lzd, &
!           tmb%linmat%s, tmb%linmat%l, tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)
     end if
  end if
  ! Write the sparse matrices
  if (mod(input%lin%output_mat_format,10) /= WF_FORMAT_NONE) then
      call timing(iproc,'write_matrices','ON')
      call write_linear_matrices(iproc,nproc,input%imethod_overlap,trim(input%dir_output),&
           input%lin%output_mat_format,tmb,at,rxyz,norder_taylor, &
           input%lin%calculate_onsite_overlap, write_SminusonehalfH=.true.)

      !temporary at the moment - to eventually be moved to more appropriate location
      !tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
      !call tmb_overlap_onsite(iproc, nproc, input%imethod_overlap, at, tmb, rxyz)
      !call tmb_overlap_onsite_rotate(iproc, nproc, input, at, tmb, rxyz, ref_frags)
      !call f_free_ptr(tmb%linmat%ovrlp_%matrix)
      call timing(iproc,'write_matrices','OF')
  end if

  ! Write the KS coefficients
  if (mod(input%lin%output_coeff_format,10) /= WF_FORMAT_NONE) then
      call write_linear_coefficients(0, trim(input%dir_output)//'KS_coeffs.bin', at, rxyz, &
           tmb%linmat%l%nfvctr, tmb%orbs%norb, tmb%linmat%l%nspin, tmb%coeff, tmb%orbs%eval)
  end if


       ! debug
       !tmb%linmat%kernel_%matrix = sparsematrix_malloc_ptr(tmb%linmat%l, DENSE_FULL, id='tmb%linmat%kernel__%matrix')
       !!call uncompress_matrix(bigdft_mpi%iproc,tmb%linmat%kernel_)
       !call uncompress_matrix2(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_%matrix_compr, tmb%linmat%kernel_%matrix)
       !if (iproc==0) then
       !   do iorb=1,tmb%orbs%norb
       !      do jorb=1,tmb%orbs%norb
       !         write(33,*) iorb,jorb,tmb%coeff(iorb,jorb),tmb%linmat%kernel_%matrix(iorb,jorb,1)
       !      end do
       !   end do
       !   write(33,*) ''
       !end if 
       !call f_free_ptr(tmb%linmat%kernel_%matrix) 
       !! end debug

  ! TEMPORARY DEBUG - plot in global box - CHECK WITH REFORMAT ETC IN LRs
  !nullify(gpsi_virt)
  !if (.false.) then
  !   gpsi_all=f_malloc_ptr(tmb%orbs%norbp*(tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f),id='gpsi_all')
  !   call to_zero(tmb%orbs%norbp*(tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f),gpsi_all)
  !   psit_large_c=f_malloc(tmb%ham_descr%collcom%ndimind_c,id='psit_large_c')
  !   psit_large_f=f_malloc(7*tmb%ham_descr%collcom%ndimind_f,id='psit_large_f')
  !   kpsit_c=f_malloc(tmb%ham_descr%collcom%ndimind_c,id='psit_large_c')
  !   kpsit_f=f_malloc(7*tmb%ham_descr%collcom%ndimind_f,id='psit_large_f')
  !   kpsi=f_malloc(tmb%ham_descr%npsidim_orbs,id='kpsi')
  !   kpsi_small=f_malloc(tmb%npsidim_orbs,id='kpsi_small')
  !   !NB might need kernel without occs to avoid including occs twice, unless tmb%occs are all 1...
  !   call small_to_large_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
  !        tmb%orbs, tmb%psi, tmb%ham_descr%psi)
  !   call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
  !        tmb%ham_descr%psi, psit_large_c, psit_large_f, tmb%ham_descr%lzd)
  !   call build_linear_combination_transposed(tmb%ham_descr%collcom, &
  !        tmb%linmat%l, tmb%linmat%kernel_, psit_large_c, psit_large_f, .true., kpsit_c, kpsit_f, iproc)
  !   call untranspose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
  !        kpsit_c, kpsit_f, kpsi, tmb%ham_descr%lzd)
  !   call large_to_small_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
  !        tmb%orbs, kpsi, kpsi_small)
  !   call f_free(psit_large_c)
  !   call f_free(psit_large_f)
  !   call f_free(kpsit_c)
  !   call f_free(kpsit_f)
  !   call f_free(kpsi)
  !end if

  !ind=1
  !!indg=1
  !gpsi=f_malloc((tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f),id='gpsi')
  !do iorbp=1,tmb%orbs%norbp
  !   iiorb=iorbp+tmb%orbs%isorb
  !   ilr = tmb%orbs%inwhichlocreg(iiorb)
  !
  !   call to_zero(tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f,gpsi)
  !   !call to_zero(tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f,gpsi_all(indg))
  !
  !!call vcopy(tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f,gpsi(1),1,gpsi_all(indg),1)
  !
  !   !call Lpsi_to_global2(iproc, tmb%Lzd%Llr(ilr)%wfd%nvctr_c+7*tmb%Lzd%Llr(ilr)%wfd%nvctr_f, &
  !   !     tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f, &
  !   !     1, 1, 1, tmb%Lzd%glr, tmb%Lzd%Llr(ilr), kpsi_small(ind), gpsi_all(indg))
  !
  !   call Lpsi_to_global2(iproc, tmb%Lzd%Llr(ilr)%wfd%nvctr_c+7*tmb%Lzd%Llr(ilr)%wfd%nvctr_f, &
  !        tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f, &
  !        1, 1, 1, tmb%Lzd%glr, tmb%Lzd%Llr(ilr), tmb%psi(ind), gpsi)
  ! 
  !   call plot_wf(trim(input%dir_output)//trim(adjustl(yaml_toa(iiorb))),1,at,1.0_dp,tmb%Lzd%glr,&
  !        tmb%Lzd%hgrids(1),tmb%Lzd%hgrids(2),tmb%Lzd%hgrids(3),rxyz,gpsi)
  !   !call plot_wf(trim(adjustl(orbname)),1,at,1.0_dp,tmb%Lzd%Llr(ilr),&
  !   !     tmb%Lzd%hgrids(1),tmb%Lzd%hgrids(2),tmb%Lzd%hgrids(3),rxyz,tmb%psi)
  ! 
  !   ind = ind + tmb%Lzd%Llr(ilr)%wfd%nvctr_c+7*tmb%Lzd%Llr(ilr)%wfd%nvctr_f
  !   !indg = indg + tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f
  !end do
  !call f_free(gpsi)
  !!call f_free(kpsi_small)

  !fakeorbs%norb=0
  !call local_analysis(iproc,nproc,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),at,rxyz,tmb%lzd%glr,&
  !     tmb%orbs,fakeorbs,gpsi_all,gpsi_virt)
  !call f_free_ptr(gpsi_all)
  ! END DEBUG 


  ! check why this is here!
  !!tmparr = sparsematrix_malloc(tmb%linmat%l,iaction=SPARSE_FULL,id='tmparr')
  !!call vcopy(tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(1), 1, tmparr(1), 1)
  !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
  call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
       tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, denspot%dpbox%ndimrhopot, &
       denspot%rhov, rho_negative)
  !!call vcopy(tmb%linmat%l%nvctr, tmparr(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
  !!call f_free(tmparr)
  if (rho_negative) then
      call corrections_for_negative_charge(iproc, nproc, at, denspot)
      !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
      !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
      !!call clean_rho(iproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
  end if


  ! Otherwise there are some problems... Check later.
  KSwfn%psi = f_malloc_ptr(1,id='KSwfn%psi')
  nullify(KSwfn%psit)

  nullify(rho,pot)

  call deallocate_local_arrays()
  call f_release_routine()

  call timing(bigdft_mpi%mpi_comm,'WFN_OPT','PR')

  contains

    !> This loop is simply copied down here such that it can again be called
    !! in a post-processing way.
    subroutine scf_kernel(nit_scc, remove_coupling_terms, update_phi)
       implicit none

       ! Calling arguments
       integer,intent(in) :: nit_scc
       logical,intent(in) :: remove_coupling_terms !<set the matrix elements coupling different atoms to zero
       logical,intent(inout) :: update_phi

       ! Local variables
       logical :: calculate_overlap, invert_overlap_matrix

       kernel_loop : do it_scc=1,nit_scc
           dmin_diag_it=dmin_diag_it+1
           ! If the hamiltonian is available do not recalculate it
           ! also using update_phi for calculate_overlap_matrix and communicate_phi_for_lsumrho
           ! since this is only required if basis changed
           if (iproc==0) then
              !if (it_scc==nit_scc) then
              !   call yaml_sequence(label='final_kernel'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)'))),advance='no')
              !else
                 call yaml_sequence(advance='no')
              !end if
              !call yaml_mapping_open(flow=.false.)
              call yaml_comment('kernel iter:'//yaml_toa(it_scc,fmt='(i6)'),hfill='-')
           end if
           ! Check whether the overlap matrix must be calculated and inverted (otherwise it has already been done)
           calculate_overlap = ((update_phi .and. .not.input%correction_co_contra))! .or. cur_it_highaccuracy==1)
           invert_overlap_matrix = (.not.target_function==TARGET_FUNCTION_IS_HYBRID .and. it_scc==1)
                                    !cur_it_highaccuracy==1)
           !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
           if(update_phi .and. can_use_ham) then! .and. info_basis_functions>=0) then
              if (input%lin%constrained_dft) then
                 !Allocate weight matrix which is used in the CDFT loop
                 weight_matrix_ = matrices_null()
                 call allocate_matrices(tmb%linmat%m, allocate_full=.false., matname='weight_matrix_', mat=weight_matrix_)
                 weight_matrix_%matrix_compr=cdft%weight_matrix_%matrix_compr
                 call extract_taskgroup_inplace(tmb%linmat%m, weight_matrix_)
                 !PB: This resets the DIIS history, effectively breaking DIIS.
                 !PB: What should be done is storing of both constrained and unconstrained matrices,
                 !PB: so we determine the extrapolation coefficients from the constrained and apply them
                 !PB: to the unconstrained to get next matrix.
                 if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then 
                    call DIIS_free(ldiis_coeff)
                    if (input%lin%extra_states==0) then
                       call DIIS_set(ldiis_coeff_hist,0.1_gp,tmb%orbs%norb*KSwfn%orbs%norbp,1,ldiis_coeff)
                    else
                       call DIIS_set(ldiis_coeff_hist,0.1_gp,tmb%orbs%norb*tmb%orbs%norbp,1,ldiis_coeff)
                    end if
                 end if
                 ! Set the DIIS for the constrained
                 call DIIS_set(30,valpha,1,1,vdiis)
                 call vcopy(tmb%orbs%norb**2,tmb%coeff(1,1),1,coeff_tmp(1,1),1)
                 vold=cdft%lag_mult

                 ! The self consistency cycle. Here we try to get a self consistent density/potential with the fixed basis.
                 cdft_loop : do cdft_it=1,100
                    call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                         infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,calculate_overlap,invert_overlap_matrix,update_phi,&
                         .false.,input%lin%extra_states,itout,it_scc,cdft_it,norder_taylor,input%lin%max_inversion_error,&
                         input%calculate_KS_residue,input%calculate_gap,energs_work,remove_coupling_terms,input%lin%coeff_factor,&
                      input%lin%pexsi_npoles,input%lin%pexsi_mumin,input%lin%pexsi_mumax,input%lin%pexsi_mu,& 
                      input%lin%pexsi_temperature,input%lin%pexsi_tol_charge, &
                         convcrit_dmin,nitdmin,input%lin%curvefit_dmin,ldiis_coeff,reorder,cdft)
                    call get_lagrange_mult(cdft_it,vgrad) 
                    ! CDFT: exit when W is converged wrt both V and rho
                    if (abs(vgrad) < cdft_charge_thresh .or. target_function==TARGET_FUNCTION_IS_TRACE) then
                       exit
                    end if
                 end do cdft_loop
                 call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%m, weight_matrix_)
                 call deallocate_matrices(weight_matrix_)
                 call DIIS_free(vdiis)
              else
                 call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                      infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,calculate_overlap,invert_overlap_matrix,update_phi,&
                      .false.,input%lin%extra_states,itout,it_scc,cdft_it,norder_taylor,input%lin%max_inversion_error,&
                      input%calculate_KS_residue,input%calculate_gap,energs_work,remove_coupling_terms,input%lin%coeff_factor,&
                      input%lin%pexsi_npoles,input%lin%pexsi_mumin,input%lin%pexsi_mumax,input%lin%pexsi_mu,& 
                      input%lin%pexsi_temperature,input%lin%pexsi_tol_charge, &
                      convcrit_dmin,nitdmin,input%lin%curvefit_dmin,ldiis_coeff,reorder)
              end if
           else
              if (input%lin%constrained_dft) then
                 !Allocate weight matrix which is used in the CDFT loop
                 weight_matrix_ = matrices_null()
                 call allocate_matrices(tmb%linmat%m, allocate_full=.false., matname='weight_matrix_', mat=weight_matrix_)
                 weight_matrix_%matrix_compr=cdft%weight_matrix_%matrix_compr
                 call extract_taskgroup_inplace(tmb%linmat%m, weight_matrix_)
                 !PB: This resets the DIIS history, effectively breaking DIIS.
                 !PB: What should be done is storing of both constrained and unconstrained matrices,
                 !PB: so we determine the extrapolation coefficients from the constrained and apply them
                 !PB: to the unconstrained to get next matrix.
                 if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then 
                    call DIIS_free(ldiis_coeff)
                    if (input%lin%extra_states==0) then
                       call DIIS_set(ldiis_coeff_hist,0.1_gp,tmb%orbs%norb*KSwfn%orbs%norbp,1,ldiis_coeff)
                    else
                       call DIIS_set(ldiis_coeff_hist,0.1_gp,tmb%orbs%norb*tmb%orbs%norbp,1,ldiis_coeff)
                    end if
                 end if
                 ! Set the DIIS for the constrained
                 call DIIS_set(30,valpha,1,1,vdiis)
                 call vcopy(tmb%orbs%norb**2,tmb%coeff(1,1),1,coeff_tmp(1,1),1)
                 vold=cdft%lag_mult

                 ! The self consistency cycle. Here we try to get a self consistent density/potential with the fixed basis.
                 cdft_loop1 : do cdft_it=1,100
                    call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                         infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,calculate_overlap,invert_overlap_matrix,update_phi,&
                         .true.,input%lin%extra_states,itout,it_scc,cdft_it,norder_taylor,input%lin%max_inversion_error,&
                         input%calculate_KS_residue,input%calculate_gap,energs_work,remove_coupling_terms,input%lin%coeff_factor,&
                      input%lin%pexsi_npoles,input%lin%pexsi_mumin,input%lin%pexsi_mumax,input%lin%pexsi_mu,& 
                      input%lin%pexsi_temperature,input%lin%pexsi_tol_charge, &
                         convcrit_dmin,nitdmin,input%lin%curvefit_dmin,ldiis_coeff,reorder,cdft)
                    call get_lagrange_mult(cdft_it,vgrad)
                    ! CDFT: exit when W is converged wrt both V and rho
                    if (abs(vgrad) < cdft_charge_thresh .or. target_function==TARGET_FUNCTION_IS_TRACE) then
                       exit
                    end if
                 end do cdft_loop1
                 call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%m, weight_matrix_)
                 call deallocate_matrices(weight_matrix_)
                 call DIIS_free(vdiis)
              else
                 call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                      infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,calculate_overlap,invert_overlap_matrix,update_phi,&
                      .true.,input%lin%extra_states,itout,it_scc,cdft_it,norder_taylor,input%lin%max_inversion_error,&
                      input%calculate_KS_residue,input%calculate_gap,energs_work,remove_coupling_terms,input%lin%coeff_factor,&
                      input%lin%pexsi_npoles,input%lin%pexsi_mumin,input%lin%pexsi_mumax,input%lin%pexsi_mu,& 
                      input%lin%pexsi_temperature,input%lin%pexsi_tol_charge, &
                      convcrit_dmin,nitdmin,input%lin%curvefit_dmin,ldiis_coeff,reorder)
              end if
           end if
           !do i=1,tmb%linmat%l%nvctr
           !    write(*,*) 'i, lernel', i, tmb%linmat%kernel_%matrix_compr(i)
           !end do
           !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)



           ! Since we do not update the basis functions anymore in this loop
           update_phi = .false.

           !EXPERIMENTAL (currently switched off)
           ! every so often during direct min want to diagonalize - figure out a good way to specify how often...
           if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION.and.iproc==0.and.dmin_diag_freq>=0) &
                print*,'COUNTDOWN',dmin_diag_freq-dmin_diag_it
           if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION&!.and.(it_scc==nit_scc.or.pnrm<convCritMix)&
                .and.dmin_diag_it>=dmin_diag_freq.and.dmin_diag_freq/=-1) then
              reorder=.true.
              !call get_coeff(iproc,nproc,LINEAR_MIXDENS_SIMPLE,KSwfn%orbs,at,rxyz,denspot,GPU,&
              !     infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,update_phi,update_phi,&
              !     .true.,ham_small,input%lin%extra_states)
              ! just diagonalize with optimized states?
              dmin_diag_it=0
           !else if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION.and.it_scc==nit_scc.and.dmin_diag_it>=dmin_diag_freq) then
           !   if (iproc==0) print*,'NOTCOUNTDOWN',pnrm,convcrit_dmin*10.0d0
           else
              reorder=.false.
           end if
           !END EXPERIMENTAL

           ! CDFT: this is the real energy here as we subtracted the constraint term from the Hamiltonian before calculating ebs
           ! Calculate the total energy.
           !if(iproc==0) write(*,'(a,9es14.6)') 'energs', &
           !    energs%ebs,energs%ekin, energs%epot, energs%eh,energs%exc,energs%evxc,energs%eexctX,energs%eion,energs%edisp
           energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
           energyDiff=energy-energyold
           energyold=energy

           ! update alpha_coeff for direct minimization steepest descents
           if(input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION .and. it_scc>1 .and.&
                ldiis_coeff%idsx == 0 .and. (.not. input%lin%curvefit_dmin)) then
              ! apply a cap so that alpha_coeff never goes below around 1.d-2 or above 2
              !SM <=1.d-10 due to the tests...
              if (energyDiff<=1.d-10 .and. ldiis_coeff%alpha_coeff < 1.8d0) then
                 ldiis_coeff%alpha_coeff=1.1d0*ldiis_coeff%alpha_coeff
              else if (ldiis_coeff%alpha_coeff > 1.7d-3) then
                 ldiis_coeff%alpha_coeff=0.5d0*ldiis_coeff%alpha_coeff
              end if
              !!if(iproc==0) write(*,*) ''
              !!if(iproc==0) write(*,*) 'alpha, energydiff',ldiis_coeff%alpha_coeff,energydiff
              if (iproc==0) then
                  call yaml_map('alpha',ldiis_coeff%alpha_coeff)
                  call yaml_map('energydiff',energydiff)
              end if
           end if

           ! Calculate the charge density.
           if (iproc==0) then
               call yaml_mapping_open('Hamiltonian update',flow=.true.)
               ! Use this subroutine to write the energies, with some
               ! fake number
               ! to prevent it from writing too much
               call write_energies(0,0,energs,0.d0,0.d0,'',.true.)
           end if
           !!tmparr = sparsematrix_malloc(tmb%linmat%l,iaction=SPARSE_FULL,id='tmparr')
           !!call vcopy(tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(1), 1, tmparr(1), 1)
           !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
           call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
                tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, &
                denspot%dpbox%ndimrhopot, &
                denspot%rhov, rho_negative)
           !!call vcopy(tmb%linmat%l%nvctr, tmparr(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
           !!call f_free(tmparr)
           if (rho_negative) then
               call corrections_for_negative_charge(iproc, nproc, at, denspot)
               !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
               !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
               !!call clean_rho(iproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
           end if


           ! Mix the density.
           if (input%lin%scf_mode/=LINEAR_MIXPOT_SIMPLE) then
              ! use it_scc+1 since we already have the density from the input guess as iteration 1
              !!rho_tmp=denspot%rhov
              call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,1.d0-alpha_mix,denspot%mix,&
                   denspot%rhov,it_scc+1,denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
                   at%astruct%cell_dim(1)*at%astruct%cell_dim(2)*at%astruct%cell_dim(3),&
                   pnrm,denspot%dpbox%nscatterarr)
               !!rho_tmp=rho_tmp-denspot%rhov
               !!tt=ddot(size(rho_tmp),rho_tmp,1,rho_tmp,1)
               !!call mpiallred(tt,1,mpi_sum,bigdft_mpi%mpi_comm)
               !!tt=tt/dble(denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%ndims(3))
               !!if (iproc==0) write(*,*) 'delta rho',tt
                   !!write(*,*) 'after mix_rhopot 1.1: pnrm', pnrm
              !SM: to make sure that the result is analogous for polarized and non-polarized calculations, to be checked...
              !write(*,*) 'old pnrm',pnrm
              !!tt1=sum(denspot%dpbox%nscatterarr(:,1))
              !!tt2=sum(denspot%dpbox%nscatterarr(:,2))
              !!pnrm = pnrm*sqrt(tt2/tt1)
              pnrm=pnrm*sqrt(real(denspot%mix%nspden,kind=8))
                   !!write(*,*) 'after mix_rhopot 1.2: pnrm', pnrm
              !write(*,*) 'new pnrm',pnrm
              call check_negative_rho(input%nspin, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, &
                   denspot%rhov, rho_negative)
              if (rho_negative) then
                  call corrections_for_negative_charge(iproc, nproc, at, denspot)
              end if

              if ((pnrm<convCritMix .or. it_scc==nit_scc)) then
                 ! calculate difference in density for convergence criterion of outer loop
                 ! ioffset is the buffer which is present for GGA calculations
                 ioffset=KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%i3xcsh
                 pnrm_out=0.d0
                 do ispin=1,input%nspin
                     ! ishift gives the start of the spin down component
                     ishift=(ispin-1)*KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d 
                     do i=1,KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p
                        pnrm_out=pnrm_out+(denspot%rhov(ishift+ioffset+i)-rhopotOld_out(ishift+ioffset+i))**2
                     end do
                 end do
                 ! To make the residue for the polarized and non-polarized case analogous
                 if (input%nspin==2) then
                     pnrm_out = pnrm_out*2.d0
                 end if

                 if (nproc > 1) then
                    call mpiallred(pnrm_out, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
                 end if
                 !pnrm_out = pnrm_out/dble(input%nspin)

                 nsize = int(KSwfn%Lzd%Glr%d%n1i,kind=8)*int(KSwfn%Lzd%Glr%d%n2i,kind=8)*int(KSwfn%Lzd%Glr%d%n3i,kind=8)
                 pnrm_out=sqrt(pnrm_out)/real(nsize,kind=8)
                 call vcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, &
                      denspot%rhov(1), 1, rhopotOld_out(1), 1)
              end if

           end if

           ! Calculate the new potential.
           !!if(iproc==0) write(*,'(1x,a)') '---------------------------------------------------------------- Updating potential.'
           if (iproc==0) then
!               if (iproc==0) call yaml_mapping_open('pot',flow=.true.)
               !call yaml_map('update potential',.true.)
           end if
           if (iproc==0) call yaml_newline()
           

           call updatePotential(input%nspin,denspot,energs)!%eh,energs%exc,energs%evxc)
           if (iproc==0) call yaml_mapping_close()

           !ii = 0
           !do i3=1,denspot%dpbox%ndims(3)
           !    do i2=1,denspot%dpbox%ndims(2)
           !        do i1=1,denspot%dpbox%ndims(1)
           !            ii = ii + 1
           !            write(200,*) 'vals', i1, i2, i3, denspot%rhov(ii)
           !        end do
           !    end do
           !end do
           !close(200)

           !!call mpi_finalize(ii)
           !!stop


           ! update occupations wrt eigenvalues (NB for directmin these aren't guaranteed to be true eigenvalues)
           ! switch off for FOE at the moment
           ! switch off for directmin too, unless we decide to reactivate calculating the expectation values the output is meaningless
           if (input%lin%scf_mode/=LINEAR_FOE .and. input%lin%scf_mode/=LINEAR_PEXSI .and. &
               input%lin%scf_mode/=LINEAR_DIRECT_MINIMIZATION) then
               !call vcopy(kswfn%orbs%norb,tmb%orbs%eval(1),1,kswfn%orbs%eval(1),1)
               ! Copy the spin up eigenvalues (or all in the case of a non-polarized calculation)
               call vcopy(kswfn%orbs%norbu,tmb%orbs%eval(1),1,kswfn%orbs%eval(1),1)
               if (input%nspin==2) then
                   ! Copy the spin down eigenvalues
                   call vcopy(kswfn%orbs%norbd,tmb%orbs%eval(tmb%linmat%l%nfvctr+1),1,kswfn%orbs%eval(kswfn%orbs%norbu+1),1)
               end if
               ! Keep the ocupations for the moment.. maybe to be activated later (with a better if statement)
               if (input%Tel > 0.0_gp) then
                   call evaltoocc(iproc,nproc,.false.,input%tel,kswfn%orbs,input%occopt)
               end if
              if (bigdft_mpi%iproc ==0) then 
                 call write_eigenvalues_data(0.1d0,kswfn%orbs,mom_vec_fake)
              end if
           end if

           ! Mix the potential
           if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
              call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,1.d0-alpha_mix,denspot%mix,&
                   denspot%rhov,it_scc+1,denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
                   at%astruct%cell_dim(1)*at%astruct%cell_dim(2)*at%astruct%cell_dim(3),&
                   pnrm,denspot%dpbox%nscatterarr)
                  !write(*,*) 'after mix_rhopot 1.1: pnrm', pnrm

              !SM: to make sure that the result is analogous for polarized and non-polarized calculations, to be checked...
              pnrm=pnrm*sqrt(real(denspot%mix%nspden,kind=8))
                  !write(*,*) 'after mix_rhopot 1.2: pnrm', pnrm
              if (pnrm<convCritMix .or. it_scc==nit_scc) then
                 ! calculate difference in density for convergence criterion of outer loop
                 ! There is no ioffset (unlike to the case of density mixing)
                 ! since also for GGA calculations there is no offset for the potential
                 pnrm_out=0.d0
                 do ispin=1,input%nspin
                     ! ishift gives the start of the spin down component
                     ishift=(ispin-1)*KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p 
                     do i=1,KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p
                        pnrm_out=pnrm_out+(denspot%rhov(ishift+i)-rhopotOld_out(ishift+i))**2
                     end do
                 end do
                 ! To make the residue for the polarized and non-polarized case analogous
                 if (input%nspin==2) then
                     pnrm_out = pnrm_out*2.d0
                 end if

                 if (nproc > 1) then
                    call mpiallred(pnrm_out, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
                 end if

                 pnrm_out=sqrt(pnrm_out)/(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*KSwfn%Lzd%Glr%d%n3i)!)*input%nspin)
                 call vcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, &
                      denspot%rhov(1), 1, rhopotOld_out(1), 1) 
              end if
           end if
 
           ! CDFT: need to pass V*w_ab to get_coeff so that it can be added to H_ab and the correct KS eqn can therefore be solved
           ! CDFT: for the first iteration this will be some initial guess for V (or from the previous outer loop)
           ! CDFT: all this will be in some extra CDFT loop
           if(input%lin%constrained_dft) then

           end if
           ! CDFT: end of CDFT loop to find V which correctly imposes constraint and corresponding density

           ! Keep the support functions fixed if they converged and the density
           ! change is below the tolerance already in the very first iteration
           if(it_scc==1 .and. pnrm<convCritMix .and.  info_basis_functions>0) then
              fix_support_functions=.true.
           end if

           ! Write some informations.
           call printSummary()

           if (pnrm<convCritMix.and.input%lin%scf_mode/=LINEAR_DIRECT_MINIMIZATION) then
               info_scf=it_scc
               if (iproc==0) then
                   !yaml output
                   !call yaml_mapping_close() !iteration
                  call yaml_flush_document()
                  !call bigdft_utils_flush(unit=6)
               end if
               exit
           else if (pnrm<convCritMix.and.input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
               if (iproc==0) then
                   !yaml output
                   !call yaml_mapping_close() !iteration
                  call yaml_flush_document()
                  !call bigdft_utils_flush(unit=6)
               end if
              exit
           !else if (pnrm<convCritMix.and.reorder) then
           !    exit
           !else if (pnrm<convCritMix) then
           !    reorder=.true.
           !    dmin_diag_it=0
           else
               info_scf=-1
           end if

           if (iproc==0) then
               !yaml output
               !call yaml_mapping_close() !iteration
              call yaml_flush_document()
              ! call bigdft_utils_flush(unit=6)
           end if

       end do kernel_loop

    end subroutine scf_kernel



    subroutine allocate_local_arrays()

      locrad = f_malloc(tmb%lzd%nlr,id='locrad')
      ! Allocate the old charge density (used to calculate the variation in the charge density)
      rhopotold_out = f_malloc(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim),id='rhopotold_out')

      fnrm_work = work_mpiaccumulate_null()
      fnrm_work%ncount = 1
      call allocate_work_mpiaccumulate(fnrm_work)
      energs_work = work_mpiaccumulate_null()
      energs_work%ncount = 4
      call allocate_work_mpiaccumulate(energs_work)

    end subroutine allocate_local_arrays


    subroutine deallocate_local_arrays()

      call f_free(locrad)
      call f_free(rhopotold_out)
      call deallocate_work_mpiaccumulate(energs_work)
      call deallocate_work_mpiaccumulate(fnrm_work)

    end subroutine deallocate_local_arrays


    subroutine check_inputguess()
      real(kind=8) :: dnrm2
      if (input%inputPsiId .hasattr. 'MEMORY') then !==101) then           !should we put 102 also?

          if (input%lin%pulay_correction) then
             ! Check the input guess by calculation the Pulay forces.

             call f_zero(tmb%ham_descr%npsidim_orbs,tmb%ham_descr%psi(1))
             call small_to_large_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
                  tmb%orbs, tmb%psi, tmb%ham_descr%psi)

             ! add get_coeff here
             ! - need some restructuring/reordering though, or addition of lots of extra initializations?!

             ! Calculate Pulay correction to the forces
             call pulay_correction(iproc, nproc, KSwfn%orbs, at, rxyz, nlpsp, input%SIC, denspot, GPU, tmb, fpulay)
             fnrm_pulay=dnrm2(3*at%astruct%nat, fpulay, 1)/sqrt(dble(at%astruct%nat))

             !if (iproc==0) write(*,*) 'fnrm_pulay',fnrm_pulay
             if (iproc==0) call yaml_map('fnrm Pulay',fnrm_pulay)

             if (fnrm_pulay>1.d-1) then !1.d-10
                if (iproc==0) then
                    call yaml_warning('The pulay force is too large after the restart. &
                         &Start over again with an AO input guess.')
                end if
                !!if (associated(tmb%psit_c)) then
                !!    call f_free_ptr(tmb%psit_c)
                !!end if
                !!if (associated(tmb%psit_f)) then
                !!    call f_free_ptr(tmb%psit_f)
                !!end if
                tmb%can_use_transposed=.false.
                nit_lowaccuracy=input%lin%nit_lowaccuracy
                nit_highaccuracy=input%lin%nit_highaccuracy
                !!input%lin%highaccuracy_conv_crit=1.d-8
                call inputguessConfinement(iproc, nproc, at, input, &
                     KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
                     rxyz, nlpsp, GPU, KSwfn%orbs, KSwfn, tmb, denspot, rhopotold, energs)
                     energs%eexctX=0.0_gp

                !already done in inputguess
                      ! CHEATING here and passing tmb%linmat%denskern instead of tmb%linmat%inv_ovrlp
                !call orthonormalizeLocalized(iproc, nproc, 0, tmb%npsidim_orbs, tmb%orbs, tmb%lzd, tmb%linmat%ovrlp, &
                !     tmb%linmat%denskern, tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)
             else if (fnrm_pulay>1.d-2) then ! 1.d2 1.d-2
                !!if (iproc==0) write(*,'(1x,a)') 'The pulay forces are rather large, so start with low accuracy.'
                if (iproc==0) then
                    call yaml_warning('The pulay forces are rather large, so start with low accuracy.')
                end if
                nit_lowaccuracy=input%lin%nit_lowaccuracy
                nit_highaccuracy=input%lin%nit_highaccuracy
             else if (fnrm_pulay>1.d-10) then !1d-10
                if (iproc==0) then
                    !write(*,'(1x,a)') 'The pulay forces are fairly large, so reoptimising basis with high accuracy only.'
                    call yaml_warning('The pulay forces are fairly large, so reoptimising basis with high accuracy only.')
                end if
                nit_lowaccuracy=0
                nit_highaccuracy=input%lin%nit_highaccuracy
             else
                if (iproc==0) write(*,'(1x,a)') &
                     'The pulay forces are fairly small, so not reoptimising basis.'
                    nit_lowaccuracy=0
                    nit_highaccuracy=0
             end if
          else
              ! Calculation of Pulay forces not possible, so always start with low accuracy
              call f_zero(fpulay)
              if (input%lin%scf_mode==LINEAR_FOE) then
                 nit_lowaccuracy=input%lin%nit_lowaccuracy
              else
                 ! double check that nit_high /= 0 and we're not in hybrid mode (otherwise we won't be doiung any iterations...)
                 if (input%lin%nit_highaccuracy /=0 .and. input%lin%nlevel_accuracy==2) then
                    nit_lowaccuracy=0
                 else
                    nit_lowaccuracy=input%lin%nit_lowaccuracy
                 end if
              end if
              nit_highaccuracy=input%lin%nit_highaccuracy
          end if
          if (input%lin%scf_mode==LINEAR_FOE .and. nit_lowaccuracy==0) then
              stop 'for FOE and restart, nit_lowaccuracy must not be zero!'
          end if
      else   !if not using restart just do the lowaccuracy like normal
          nit_lowaccuracy=input%lin%nit_lowaccuracy
          nit_highaccuracy=input%lin%nit_highaccuracy
      end if
    end subroutine check_inputguess

    subroutine check_for_exit()
      implicit none

      exit_outer_loop=.false.

      if (input%lin%nlevel_accuracy==2) then
          if(lowaccur_converged) then
              !cur_it_highaccuracy=cur_it_highaccuracy+1
              if(cur_it_highaccuracy==nit_highaccuracy) then
                  exit_outer_loop=.true.
              else if (pnrm_out<input%lin%highaccuracy_conv_crit) then
                  exit_outer_loop=.true.
              end if
          end if
      else if (input%lin%nlevel_accuracy==1) then
          if (itout==nit_lowaccuracy .or. pnrm_out<input%lin%highaccuracy_conv_crit) then
              exit_outer_loop=.true.
          end if
      end if

    end subroutine check_for_exit

    !> Print a short summary of some values calculated during the last iteration in the self
    !! consistency cycle.
    subroutine printSummary()
      implicit none

      if(iproc==0) then
          call yaml_sequence_open('summary',flow=.true.)
          call yaml_mapping_open()
          if(input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
              call yaml_map('kernel method','DMIN')
          else if (input%lin%scf_mode==LINEAR_FOE) then
              call yaml_map('kernel method','FOE')
          else if (input%lin%scf_mode==LINEAR_PEXSI) then
              call yaml_map('kernel method','PEXSI')
          else
              call yaml_map('kernel method','DIAG')
          end if

          if (input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE .or.  input%lin%scf_mode==LINEAR_FOE &
              .or. input%lin%scf_mode==LINEAR_PEXSI .or. input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
              call yaml_map('mix entity','DENS')
          else if (input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
              call yaml_map('mix entity','POT')
          end if
          call yaml_map('mix hist',mix_hist)
          call yaml_map('conv crit',convCritMix,fmt='(es8.2)')

          call yaml_newline()
          call yaml_map('iter',it_scc,fmt='(i6)')
          call yaml_map('delta',pnrm,fmt='(es9.2)')
          call yaml_map('energy',energy,fmt='(es24.17)')
          call yaml_map('D',energyDiff,fmt='(es10.3)')

          if (input%lin%constrained_dft) then
              call yaml_map('Tr(KW)',ebs,fmt='(es14.4)')
          end if     

          call yaml_mapping_close()
          call yaml_sequence_close()
      end if

    end subroutine printSummary

    !> Print a short summary of some values calculated during the last iteration in the self
    !! consistency cycle.
    subroutine print_info(final)
      implicit none

      real(kind=8) :: energyDiff, mean_conf
      logical, intent(in) :: final
      integer :: ii

      energyDiff = energy - energyoldout

      mean_conf=0.d0
      do iorb=1,tmb%orbs%norbp
          mean_conf=mean_conf+tmb%confdatarr(iorb)%prefac
      end do

      if (nproc > 1) then
         call mpiallred(mean_conf, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      mean_conf=mean_conf/dble(tmb%orbs%norb)

      ! Print out values related to two iterations of the outer loop.
      if(iproc==0.and.(.not.final)) then

          call yaml_comment('Summary of both steps',hfill='=')
          call yaml_sequence_open('self consistency summary',label=&
              'it_sc'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)'))))
          call yaml_sequence(advance='no')
          call yaml_mapping_open(flow=.true.)
          call yaml_map('iter',itout)
          if(target_function==TARGET_FUNCTION_IS_TRACE) then
              call yaml_map('Omega','TRACE')
          else if(target_function==TARGET_FUNCTION_IS_ENERGY) then
              call yaml_map('Omega','ENERGY')
          else if(target_function==TARGET_FUNCTION_IS_HYBRID) then
              call yaml_map('Omega','HYBRID')
          end if
          if (target_function==TARGET_FUNCTION_IS_HYBRID) then
              call yaml_map('mean conf prefac',mean_conf,fmt='(es8.2)')
              call yaml_map('damping',tmb%damping_factor_confinement,fmt='(es8.2)')
          end if
          if(info_basis_functions<=0) then
              call yaml_warning('support function optimization not converged')
              call yaml_newline()
          else
              call yaml_map('iterations to converge support functions',info_basis_functions)
          end if
          if(input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
              call yaml_map('kernel optimization','DIRMIN')
          else if (input%lin%scf_mode==LINEAR_FOE) then
              call yaml_map('kernel optimization','FOE')
          else
              call yaml_map('kernel optimization','DIAG')
          end if
          if(info_scf<0) then
              call yaml_warning('density optimization not converged')
              call yaml_newline()
          else
              call yaml_map('iterations to converge kernel optimization',info_scf)
              call yaml_newline()
          end if
          if(input%lin%scf_mode/=LINEAR_MIXPOT_SIMPLE) then
             if (.not. lowaccur_converged) then
                 call yaml_map('iter low',itout)
             else
                 call yaml_map('iter high',itout)
             end if
                 call yaml_map('delta out',pnrm_out,fmt='(es10.3)')
                 call yaml_map('energy',energy,fmt='(es27.17)')
                 call yaml_map('D',energyDiff,fmt='(es10.3)')
          else if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
             if (.not. lowaccur_converged) then
                 call yaml_map('iter low',itout)
             else
                 call yaml_map('iter high',itout)
             end if
             call yaml_map('delta out',pnrm_out,fmt='(es10.3)')
             call yaml_map('energy',energy,fmt='(es27.17)')
             call yaml_map('D',energyDiff,fmt='(es10.3)')
          end if
          call yaml_mapping_close()

          !when convergence is reached, use this block
      else if (iproc==0.and.final) then
          call yaml_comment('final results',hfill='=')
          call yaml_sequence_open('self consistency summary')
          call yaml_sequence(advance='no')
          call yaml_mapping_open(flow=.true.)
          call yaml_map('iter',itout)
          call write_energies(0,0,energs,0.d0,0.d0,'',.true.)
          if (input%lin%scf_mode/=LINEAR_MIXPOT_SIMPLE) then
             if (.not. lowaccur_converged) then
                 call yaml_map('iter low',itout)
                 call yaml_map('delta out',pnrm_out,fmt='(es10.3)')
                 call yaml_map('energy',energy,fmt='(es27.17)')
                 call yaml_map('D',energyDiff,fmt='(es10.3)')
                 call yaml_comment('FINAL')
                 call yaml_mapping_close()
             else
                 call yaml_map('iter high',itout)
                 call yaml_map('delta out',pnrm_out,fmt='(es10.3)')
                 call yaml_map('energy',energy,fmt='(es27.17)')
                 call yaml_map('D',energyDiff,fmt='(es10.3)')
                 call yaml_comment('FINAL')
                 call yaml_mapping_close()
             end if
          else if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
             if (.not. lowaccur_converged) then
                 call yaml_map('iter low',itout)
                 call yaml_map('delta out',pnrm_out,fmt='(es10.3)')
                 call yaml_map('energy',energy,fmt='(es27.17)')
                 call yaml_map('D',energyDiff,fmt='(es10.3)')
                 call yaml_comment('FINAL')
                 call yaml_mapping_close()
             else
                 call yaml_map('iter high',itout)
                 call yaml_map('delta out',pnrm_out,fmt='(es10.3)')
                 call yaml_map('energy',energy,fmt='(es27.17)')
                 call yaml_map('D',energyDiff,fmt='(es10.3)')
                 call yaml_comment('FINAL')
                 call yaml_mapping_close()
             end if
          end if
       end if

       call yaml_flush_document()
       !call bigdft_utils_flush(unit=6)
    call yaml_sequence_close()

    if (.not.final .and. input%adjust_kernel_iterations) then
        ! Determine whether the sign of the energy change is the same as in the previous iteration
        ! (i.e. whether the energy continues to increase or decrease)
        tt = sign(energyDiff,sign_of_energy_change)
        if (energyDiff/=0.d0) then
            if (tt/energyDiff>0.d0) then
                ! same sign, everything ok
            else if (abs(energyDiff/energy)>1.d-7) then
                nit_energyoscillation = nit_energyoscillation + 1
                if (iproc==0) then
                    call yaml_warning('oscillation of the energy, increase counter')
                    call yaml_map('energy_scillation_counter',nit_energyoscillation)
                end if
            end if
        end if
        if (nit_energyoscillation>1) then
            nit_scc = nit_scc + 1
            ii = 3*max(input%lin%nitSCCWhenFixed_lowaccuracy,input%lin%nitSCCWhenFixed_highaccuracy)
            if (nit_scc>ii) then
                if (iproc==0) call yaml_map('nit_scc reached maximum, reset to',ii)
                nit_scc = ii
            end if
            nit_energyoscillation = 0
            if (iproc==0) call yaml_map('new nit_scc',nit_scc)
            ! Needed for low/high accuracy... maybe to be cleaned
            nit_scc_changed = nit_scc
            keep_value = .true.
        end if
    end if


    ! Determine the sign of the energy change
    sign_of_energy_change = sign(1.d0,energyDiff)


    end subroutine print_info


    subroutine intermediate_forces()
      use module_forces, only: clean_forces
      implicit none
      ! Local variables
      real(kind=8) :: eh_tmp, exc_tmp, evxc_tmp, eexctX_tmp
      real(kind=8) :: fnoise, pressure, ehart_fake
      real(kind=8),dimension(6) :: ewaldstr, hstrten, xcstr, strten
      real(kind=8),dimension(:),allocatable :: rhopot_work
          real(kind=8),dimension(:,:),allocatable :: fxyz

      ! TEST: calculate forces here ####################################################
      fxyz=f_malloc((/3,at%astruct%nat/),id='fxyz')
      ewaldstr=1.d100
      hstrten=1.d100
      xcstr=1.d100
      eh_tmp=energs%eh
      exc_tmp=energs%exc
      evxc_tmp=energs%evxc
      eexctX_tmp=energs%eexctX

      ioffset=kswfn%Lzd%Glr%d%n1i*kswfn%Lzd%Glr%d%n2i*denspot%dpbox%i3xcsh

      if (denspot%dpbox%ndimpot>0) then
          denspot%pot_work=f_malloc_ptr(denspot%dpbox%ndimpot,id='denspot%dpbox%ndimpot')

      else
          denspot%pot_work=f_malloc_ptr(1,id='denspot%dpbox%ndimpot')
      end if
      if (denspot%dpbox%ndimrhopot>0) then
          rhopot_work=f_malloc(denspot%dpbox%ndimrhopot,id='rhopot_work')
      else
          rhopot_work=f_malloc(1,id='rhopot_work')
      end if
      call vcopy(denspot%dpbox%ndimrhopot,denspot%rhov(1),1,rhopot_work(1),1)


      !!tmparr = sparsematrix_malloc(tmb%linmat%l,iaction=SPARSE_FULL,id='tmparr')
      !!call vcopy(tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(1), 1, tmparr(1), 1)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
      call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
           tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, denspot%dpbox%ndimrhopot, &
           denspot%rhov, rho_negative)
      !!call vcopy(tmb%linmat%l%nvctr, tmparr(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
      !!call f_free(tmparr)
      if (rho_negative) then
          call corrections_for_negative_charge(iproc, nproc, at, denspot)
          !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
          !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
          !!call clean_rho(iproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      end if
      !denspot%rho_work=f_malloc_ptr(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim),id='denspot%rho_work')
      denspot%rho_work=f_malloc_ptr(denspot%dpbox%ndimpot,id='denspot%rho_work')
      !call vcopy(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim),&
      !     denspot%rhov(1),1,denspot%rho_work(1),1)
      call vcopy(denspot%dpbox%ndimpot,denspot%rhov(ioffset+1),1,denspot%rho_work(1),1)
      if (denspot%dpbox%nrhodim==2) then
          call axpy(denspot%dpbox%ndimpot,1.d0,denspot%rhov(ioffset+denspot%dpbox%ndimpot+1),1,denspot%rho_work(1),1)
      end if

      call updatePotential(input%nspin,denspot,energs)!%eh,energs%exc,energs%evxc)

      ! Density already present in denspot%rho_work
      call vcopy(denspot%dpbox%ndimpot,denspot%rho_work(1),1,denspot%pot_work(1),1)
      call H_potential('D',denspot%pkernel,denspot%pot_work,denspot%pot_work,ehart_fake,&
           0.0_dp,.false.,stress_tensor=hstrten)
      
      KSwfn%psi=f_malloc_ptr(1,id='KSwfn%psi')
      fpulay=0.d0
      !this is associated but not used in the routine for linear scaling
      call orbital_basis_associate(ob,orbs=KSwfn%orbs,Lzd=KSwfn%Lzd)
      call calculate_forces(iproc,nproc,denspot%pkernel%mpi_env%nproc,KSwfn%Lzd%Glr,at,ob,nlpsp,rxyz,& 
           KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
           denspot%dpbox,&
           denspot%dpbox%i3s+denspot%dpbox%i3xcsh,denspot%dpbox%n3p,denspot%dpbox%nrhodim,&
           .false.,denspot%dpbox%ngatherarr,denspot%rho_work,&
           denspot%pot_work,denspot%V_XC,size(KSwfn%psi),KSwfn%psi,fion,fdisp,fxyz,&
           input%calculate_strten,ewaldstr,hstrten,xcstr,strten,pressure,denspot%psoffset,1,tmb,fpulay)
      call orbital_basis_release(ob)
      call clean_forces(iproc,at%astruct,rxyz,fxyz,fnoise)
      if (iproc == 0) call write_forces(at%astruct,fxyz)

      call f_free(fxyz)
      call f_free_ptr(KSwfn%psi)

      call vcopy(denspot%dpbox%ndimrhopot,rhopot_work(1),1,denspot%rhov(1),1)
      energs%eh=eh_tmp
      energs%exc=exc_tmp
      energs%evxc=evxc_tmp
      energs%eexctX=eexctX_tmp

      call f_free(rhopot_work)
      call f_free_ptr(denspot%rho_work)
      call f_free_ptr(denspot%pot_work)


    end subroutine intermediate_forces


    
    subroutine get_lagrange_mult(cdft_it, vgrad)
       implicit none

       ! Calling arguments
       integer, intent(in) :: cdft_it
       real(kind=gp), intent(out) :: vgrad

       
         ! CDFT: Calculate gradient of V=Tr[Kw]-Nc 
         call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%l,tmb%linmat%m, &
              tmb%linmat%kernel_,weight_matrix_,&
              ebs,tmb%coeff,KSwfn%orbs,tmb%orbs,.false.)
         vgrad=ebs-cdft%charge

         ! CDFT: update V (maximizing E wrt V)
         ! CDFT: we updated the kernel in get_coeff so 1st deriv of W wrt V becomes Tr[Kw]-Nc as in CONQUEST
         ! CDFT: 2nd deriv more problematic?
         ! CDFT: use simplest possible scheme for now

         !if (iproc==0) write(*,*) ''
         !if (iproc==0) write(*,'(a,I4,2x,6(ES12.4e2,2x),2(ES16.6e2,2x))') &
         !     'itc, N, Tr(KW), Tr(KW)-N, V*(Tr(KW)-N), V, Vold, EBS, energy',&
         !     cdft_it,cdft%charge,ebs,vgrad,cdft%lag_mult*vgrad,cdft%lag_mult,vold,energs%ebs,energy
         if (iproc==0) then
            call yaml_sequence_open('CDFT',flow=.true.)
            call yaml_map('itc',cdft_it)
            call yaml_map('N',cdft%charge,fmt='(es12.2)')
            call yaml_map('Tr(KW)',ebs,fmt='(es14.4)')
            !call yaml_map('Tr(KW)-N',vgrad)
            call yaml_map('Vc',cdft%lag_mult,fmt='(es12.2)')
            call yaml_map('energy',energy,fmt='(es14.4)')
            call yaml_sequence_close()
         end if

         ! CDFT: exit when W is converged wrt both V and rho
         if (abs(vgrad) < cdft_charge_thresh .or. target_function==TARGET_FUNCTION_IS_TRACE) then
            return
         end if

         call timing(iproc,'constraineddft','ON')

         ! assuming density mixing/no mixing not potential mixing
         call vcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, &
              rhopotOld_out(1), 1, rhopotOld(1), 1) 

         call vcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, &
              rhopotOld_out(1), 1, denspot%rhov(1), 1)
         call timing(iproc,'constraineddft','OF')

         call updatePotential(input%nspin,denspot,energs)!%eh,energs%exc,energs%evxc)

         call timing(iproc,'constraineddft','ON')
         ! reset coeffs as well
         call vcopy(tmb%orbs%norb**2,coeff_tmp(1,1),1,tmb%coeff(1,1),1)

         if (.false.) then ! diis
            vdiis%mids=mod(vdiis%ids,vdiis%idsx)+1
            vdiis%ids=vdiis%ids+1
            vold=cdft%lag_mult
            call diis_opt(0,1,1,0,1,(/0/),(/1/),1,&
               cdft%lag_mult,-vgrad,vdiis) 
            !call diis_opt(iproc,nproc,1,0,1,(/iproc/),(/1/),1,&
            !   cdft%lag_mult,-vgrad,vdiis) 
         else if (.false.) then !sd
            if (abs(vgrad)<abs(vgrad_old)) then
               valpha=valpha*1.1d0
            else
               valpha=valpha*0.6d0
            end if
            vold=cdft%lag_mult
            cdft%lag_mult=cdft%lag_mult+valpha*vgrad
         else if (cdft_it==1) then !first step newton
            vold=cdft%lag_mult
            !debug:
            !if (iproc==0) write(*,'(a,I4,2x,6(ES16.6e3,2x))') 'itn, V, Vg',&
            !     cdft_it,cdft%lag_mult,vgrad
            cdft%lag_mult=cdft%lag_mult*2.0_gp
         else ! newton
            vgrad2=(vgrad-vgrad_old)/(cdft%lag_mult-vold)
            !debug:
            !if (iproc==0) write(*,'(a,I4,2x,6(ES16.6e3,2x))') 'itn, V, Vold, Vg, Vgold, Vg2, Vg/Vg2',&
            !     cdft_it,cdft%lag_mult,vold,vgrad,vgrad_old,vgrad2,vgrad/vgrad2
            vold_tmp=cdft%lag_mult
            cdft%lag_mult=vold-vgrad_old/vgrad2
            vold=vold_tmp
         end if

         vgrad_old=vgrad
         call timing(iproc,'constraineddft','OF')


    end subroutine get_lagrange_mult

end subroutine linearScaling



!note that this isn't 'correct' for fake fragment case
subroutine output_fragment_rotations(iproc,nat,rxyz,iformat,filename,input_frag,ref_frags)
  use module_base
  use module_types
  use yaml_output
  use module_fragments
  !use internal_io
  use public_enums
  implicit none

  integer, intent(in) :: iproc, iformat, nat
  character(len=*), intent(in) :: filename
  real(gp), dimension(3,nat), intent(in) :: rxyz
  type(fragmentInputParameters), intent(in) :: input_frag
  type(system_fragment), dimension(input_frag%nfrag_ref), intent(in) :: ref_frags
  !Local variables
  real(gp), parameter :: W_tol=1.e-3_gp
  integer :: ifrag, jfrag, ifrag_ref, jfrag_ref, iat, isfat, jsfat,itoo_big
  real(kind=gp), dimension(:,:), allocatable :: rxyz_ref, rxyz_new
  real(kind=gp) :: null_axe, error
  type(fragment_transformation) :: frag_trans
  character(len=*), parameter :: subname='output_fragment_rotations'

  if (iproc==0) then

     null_axe=1.0d0/dsqrt(3.0d0)

     if(iformat == WF_FORMAT_PLAIN) then
        open(99, file=filename//'rotations.bin', status='unknown',form='formatted')
        write(99,'(a)') '#Label, name, label, name, angle, axis'
     else
        open(99, file=filename//'rotations.bin', status='unknown',form='unformatted')
        write(99) '#Label, name, label, name, angle, axis'
     end if

     jsfat=0
     itoo_big=0
     do jfrag=1,input_frag%nfrag
        jfrag_ref=input_frag%frag_index(jfrag)
        isfat=0
        do ifrag=1,input_frag%nfrag
           ifrag_ref=input_frag%frag_index(ifrag)

           ! only calculate rotations if same type of reference fragments
           if (jfrag_ref/=ifrag_ref) then
              if (iformat==WF_FORMAT_PLAIN) then
                 write(99,'(2(a,1x,I5,1x),F12.6,2x,3(F12.6,1x),6(1x,F18.6))')  trim(input_frag%label(ifrag_ref)),ifrag,&
                      trim(input_frag%label(jfrag_ref)),jfrag,-1.0d0,null_axe,null_axe,null_axe,&
                      -1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0
              else
                 write(99) trim(input_frag%label(ifrag_ref)),ifrag,&
                      trim(input_frag%label(jfrag_ref)),jfrag,-1.0d0,null_axe,null_axe,null_axe,&
                      -1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0
              end if
              isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat
              cycle
           else if (ifrag==jfrag) then
              if (iformat==WF_FORMAT_PLAIN) then
                 write(99,'(2(a,1x,I5,1x),F12.6,2x,3(F12.6,1x),6(1x,F18.6))')  trim(input_frag%label(ifrag_ref)),ifrag,&
                      trim(input_frag%label(jfrag_ref)),jfrag,0.0d0,null_axe,null_axe,null_axe,&
                      0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0
              else
                 write(99)  trim(input_frag%label(ifrag_ref)),ifrag,&
                      trim(input_frag%label(jfrag_ref)),jfrag,0.0d0,null_axe,null_axe,null_axe,&
                      0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0
              end if
              isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat
              cycle
           end if

           rxyz_ref = f_malloc((/3,ref_frags(ifrag_ref)%astruct_frg%nat/),id='rxyz_ref')
           rxyz_new = f_malloc((/3,ref_frags(ifrag_ref)%astruct_frg%nat/),id='rxyz_new')

           do iat=1,ref_frags(ifrag_ref)%astruct_frg%nat
              rxyz_new(:,iat)=rxyz(:,isfat+iat)
              rxyz_ref(:,iat)=rxyz(:,jsfat+iat)
           end do

           ! use center of fragment for now, could later change to center of symmetry
           frag_trans%rot_center=frag_center(ref_frags(jfrag_ref)%astruct_frg%nat,rxyz_ref)
           frag_trans%rot_center_new=frag_center(ref_frags(ifrag_ref)%astruct_frg%nat,rxyz_new)

           ! shift rxyz wrt center of rotation
           do iat=1,ref_frags(ifrag_ref)%astruct_frg%nat
              rxyz_ref(:,iat)=rxyz_ref(:,iat)-frag_trans%rot_center
              rxyz_new(:,iat)=rxyz_new(:,iat)-frag_trans%rot_center_new
           end do

           call find_frag_trans(ref_frags(ifrag_ref)%astruct_frg%nat,rxyz_ref,rxyz_new,frag_trans,error)
           if (error > W_tol) call f_increment(itoo_big)
           call f_free(rxyz_ref)
           call f_free(rxyz_new)

           if (iformat==WF_FORMAT_PLAIN) then
              write(99,'(2(a,1x,I5,1x),F12.6,2x,3(F12.6,1x),6(1x,F18.6))') trim(input_frag%label(ifrag_ref)),ifrag,&
                   trim(input_frag%label(jfrag_ref)),jfrag,frag_trans%theta/(4.0_gp*atan(1.d0)/180.0_gp),frag_trans%rot_axis,&
                   frag_trans%rot_center,frag_trans%rot_center_new
           else
              write(99) trim(input_frag%label(ifrag_ref)),ifrag,&
                   trim(input_frag%label(jfrag_ref)),jfrag,frag_trans%theta/(4.0_gp*atan(1.d0)/180.0_gp),frag_trans%rot_axis,&
                   frag_trans%rot_center,frag_trans%rot_center_new
           end if
           isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat
        end do
        jsfat=jsfat+ref_frags(jfrag_ref)%astruct_frg%nat
     end do

     close(99)
     if (itoo_big > 0) call yaml_warning('Found (again) '//itoo_big//' warning of high Wahba cost functions')
  end if

end subroutine output_fragment_rotations

!> Check the relative weight which the support functions have at the
!! boundaries of the localization regions.
subroutine get_boundary_weight(iproc, nproc, orbs, lzd, atoms, crmult, nsize_psi, psi, crit)
  use module_base
  use module_types, only: orbitals_data, local_zone_descriptors
  use module_atoms, only: atoms_data
  use yaml_output
  use dynamic_memory
  use locreg_operations, only: boundary_weight
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  type(atoms_data),intent(in) :: atoms
  real(kind=8),intent(in) :: crmult
  integer,intent(in) :: nsize_psi
  real(kind=8),dimension(nsize_psi),intent(in) :: psi
  real(kind=8),intent(in) :: crit

  ! Local variables
  integer :: iorb, iiorb, ilr, ind, iat, iatype, nwarnings
  real(kind=8) :: atomrad, rad, boundary, weight_normalized, maxweight, meanweight
  real(kind=8),dimension(:),allocatable :: maxweight_types, meanweight_types
  integer,dimension(:),allocatable :: nwarnings_types, nsf_per_type

  call f_routine(id='get_boundary_weight')

  maxweight_types = f_malloc0(atoms%astruct%ntypes,id='maxweight_types')
  meanweight_types = f_malloc0(atoms%astruct%ntypes,id='maxweight_types')
  nwarnings_types = f_malloc0(atoms%astruct%ntypes,id='nwarnings_types')
  nsf_per_type = f_malloc0(atoms%astruct%ntypes,id='nsf_per_type')

  if (iproc==0) then
     call yaml_sequence(advance='no')
  end if


  nwarnings = 0
  maxweight = 0.d0
  meanweight = 0.d0
  if (orbs%norbp>0) then
     ind = 1
     do iorb=1,orbs%norbp
        iiorb = orbs%isorb + iorb
        ilr = orbs%inwhichlocreg(iiorb)

        iat = orbs%onwhichatom(iiorb)
        iatype = atoms%astruct%iatype(iat)
        atomrad = atoms%radii_cf(iatype,1)*crmult
        rad = atoms%radii_cf(atoms%astruct%iatype(iat),1)*crmult

        nsf_per_type(iatype) = nsf_per_type(iatype ) + 1

        weight_normalized = boundary_weight(lzd%hgrids,lzd%glr,lzd%llr(ilr),rad,psi(ind))

        ind = ind + lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f

        meanweight = meanweight + weight_normalized
        maxweight = max(maxweight,weight_normalized)
        meanweight_types(iatype) = meanweight_types(iatype) + weight_normalized
        maxweight_types(iatype) = max(maxweight_types(iatype),weight_normalized)
        if (weight_normalized>crit) then
           nwarnings = nwarnings + 1
           nwarnings_types(iatype) = nwarnings_types(iatype) + 1
        end if
     end do
     if (ind/=nsize_psi+1) then
        call f_err_throw('ind/=nsize_psi ('//trim(yaml_toa(ind))//'/='//trim(yaml_toa(nsize_psi))//')', &
             err_name='BIGDFT_RUNTIME_ERROR')
     end if
  end if

  ! Sum up among all tasks... could use workarrays
  if (nproc>1) then
     call mpiallred(nwarnings, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
     call mpiallred(meanweight, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
     call mpiallred(maxweight, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
     call mpiallred(nwarnings_types, mpi_sum, comm=bigdft_mpi%mpi_comm)
     call mpiallred(meanweight_types, mpi_sum, comm=bigdft_mpi%mpi_comm)
     call mpiallred(maxweight_types, mpi_max, comm=bigdft_mpi%mpi_comm)
     call mpiallred(nsf_per_type, mpi_sum, comm=bigdft_mpi%mpi_comm)
  end if
  meanweight = meanweight/real(orbs%norb,kind=8)
  do iatype=1,atoms%astruct%ntypes
     meanweight_types(iatype) = meanweight_types(iatype)/real(nsf_per_type(iatype),kind=8)
  end do
  if (iproc==0) then
     call yaml_sequence_open('Check boundary values')
     call yaml_sequence(advance='no')
     call yaml_mapping_open(flow=.true.)
     call yaml_map('type','overall')
     call yaml_map('mean / max value',(/meanweight,maxweight/),fmt='(2es9.2)')
     call yaml_map('warnings',nwarnings)
     call yaml_mapping_close()
     do iatype=1,atoms%astruct%ntypes
        call yaml_sequence(advance='no')
        call yaml_mapping_open(flow=.true.)
        call yaml_map('type',trim(atoms%astruct%atomnames(iatype)))
        call yaml_map('mean / max value',(/meanweight_types(iatype),maxweight_types(iatype)/),fmt='(2es9.2)')
        call yaml_map('warnings',nwarnings_types(iatype))
        call yaml_mapping_close()
     end do
     call yaml_sequence_close()
  end if

  ! Print the warnings
  if (nwarnings>0) then
     if (iproc==0) then
        call yaml_warning('The support function localization radii might be too small, got'&
             &//trim(yaml_toa(nwarnings))//' warnings')
     end if
  end if

  call f_free(maxweight_types)
  call f_free(meanweight_types)
  call f_free(nwarnings_types)
  call f_free(nsf_per_type)

  call f_release_routine()

end subroutine get_boundary_weight
