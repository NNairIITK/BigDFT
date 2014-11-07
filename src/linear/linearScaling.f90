!> @file
!!  Routines used by the linear scaling version
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine linearScaling(iproc,nproc,KSwfn,tmb,at,input,rxyz,denspot,rhopotold,nlpsp,GPU,&
           energs,energy,fpulay,infocode,ref_frags,cdft, &
           fdisp, fion)
 
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => linearScaling
  use yaml_output
  use module_fragments
  use constrained_dft
  use diis_sd_optimization
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use communications_base, only: allocate_p2pComms_buffer, &
                                 deallocate_p2pComms_buffer
  use communications, only: synchronize_onesided_communication
  use sparsematrix_base, only: sparse_matrix, sparse_matrix_null, deallocate_sparse_matrix, &
                               matrices_null, allocate_matrices, deallocate_matrices, &
                               sparsematrix_malloc, sparsematrix_malloc_ptr, assignment(=), SPARSE_FULL
  use sparsematrix, only: gather_matrix_from_taskgroups_inplace, extract_taskgroup_inplace
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(atoms_data),intent(inout) :: at
  type(input_variables),intent(in) :: input ! need to hack to be inout for geopt changes
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
  real(kind=8),dimension(:),allocatable :: rhopotold_out
  real(kind=8) :: energyold, energyDiff, energyoldout, fnrm_pulay, convCritMix, convCritMix_init
  type(localizedDIISParameters) :: ldiis
  type(DIIS_obj) :: ldiis_coeff, vdiis
  logical :: can_use_ham, update_phi, locreg_increased, reduce_conf, orthonormalization_on
  logical :: fix_support_functions
  integer :: itype, istart, nit_lowaccuracy, nit_highaccuracy
  integer :: ldiis_coeff_hist, nitdmin
  logical :: ldiis_coeff_changed
  integer :: mix_hist, info_basis_functions, nit_scc, cur_it_highaccuracy
  real(kind=8) :: pnrm_out, alpha_mix, ratio_deltas, convcrit_dmin, tt1, tt2
  logical :: lowaccur_converged, exit_outer_loop, calculate_overlap, invert_overlap_matrix
  real(kind=8),dimension(:),allocatable :: locrad
  integer:: target_function, nit_basis
  
  real(kind=gp) :: ebs, vgrad_old, vgrad, valpha, vold, vgrad2, vold_tmp, conv_crit_TMB, best_charge_diff, cdft_charge_thresh
  real(kind=gp), allocatable, dimension(:,:) :: coeff_tmp
  integer :: jorb, cdft_it, nelec, iat, ityp, norder_taylor, ispin, ishift
  integer :: dmin_diag_it, dmin_diag_freq, ioffset
  logical :: reorder, rho_negative
  real(wp), dimension(:,:,:), pointer :: mom_vec_fake
  type(matrices) :: weight_matrix_
  real(kind=8) :: sign_of_energy_change
  integer :: nit_energyoscillation

  real(8),dimension(:),allocatable :: rho_tmp, tmparr
  real(8) :: tt, ddot

  !!rho_tmp = f_malloc(size(denspot%rhov),id='rho_tmp')

  call timing(iproc,'linscalinit','ON') !lr408t

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

  call timing(iproc,'linscalinit','OF') !lr408t

  ! Check the quality of the input guess
  call check_inputguess()

  call timing(iproc,'linscalinit','ON') !lr408t

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
     cdft_charge_thresh=1.e-2
     call timing(iproc,'constraineddft','OF')
  end if


  ! modify tmb%orbs%occup, as we normally use orbs%occup elsewhere
  if (input%lin%extra_states>0) then
     call to_zero(tmb%orbs%norb,tmb%orbs%occup(1))
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
     tmb%orbs%occup=1.0d0
  end if

  ! if we want to ignore read in coeffs and diag at start - EXPERIMENTAL
  if (input%lin%diag_start .and. input%inputPsiId==INPUT_PSI_DISK_LINEAR) then
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
         call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
         !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
         !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
         !!call clean_rho(iproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
     end if
     ! Calculate the new potential.
     !if(iproc==0) write(*,'(1x,a)') '---------------------------------------------------------------- Updating potential.'
     !if (iproc==0) call yaml_map('update potential',.true.)
     if (iproc==0) call yaml_mapping_open('update pot',flow=.true.)
     call updatePotential(input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
  end if

  call timing(iproc,'linscalinit','OF') !lr408t


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

      if (input%lin%nlevel_accuracy==2) then
          ! Check whether the low accuracy part (i.e. with strong confining potential) has converged.
          call check_whether_lowaccuracy_converged(itout, nit_lowaccuracy, input%lin%lowaccuracy_conv_crit, &
               lowaccur_converged, pnrm_out)
          ! Set all remaining variables that we need for the optimizations of the basis functions and the mixing.
          call set_optimization_variables(input, at, tmb%orbs, tmb%lzd%nlr, tmb%orbs%onwhichatom, tmb%confdatarr, &
               convCritMix_init, lowaccur_converged, nit_scc, mix_hist, alpha_mix, locrad, target_function, nit_basis, &
               convcrit_dmin, nitdmin, conv_crit_TMB)
      else if (input%lin%nlevel_accuracy==1 .and. itout==1) then
          call set_variables_for_hybrid(iproc, tmb%lzd%nlr, input, at, tmb%orbs, &
               lowaccur_converged, tmb%damping_factor_confinement, tmb%confdatarr, &
               target_function, nit_basis, nit_scc, mix_hist, locrad, alpha_mix, convCritMix_init, conv_crit_TMB)
               convcrit_dmin=input%lin%convCritDmin_highaccuracy
               nitdmin=input%lin%nItdmin_highaccuracy
      end if

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
          ! Adjust the confining potential if required.
          call adjust_locregs_and_confinement(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
               at, input, rxyz, KSwfn, tmb, denspot, nlpsp, ldiis, locreg_increased, lowaccur_converged, locrad)
          orthonormalization_on=.true.


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
                     input%purification_quickreturn,&
                     input%calculate_KS_residue,input%calculate_gap,&
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
                  info_basis_functions,nlpsp,input%lin%scf_mode,ldiis,input%SIC,tmb,energs, &
                  input%lin%nItPrecond,target_function,input%lin%correctionOrthoconstraint,&
                  nit_basis,&
                  ratio_deltas,orthonormalization_on,input%lin%extra_states,itout,conv_crit_TMB,input%experimental_mode,&
                  input%lin%early_stop, input%lin%gnrm_dynamic, input%lin%min_gnrm_for_dynamic, &
                  can_use_ham, norder_taylor, input%lin%max_inversion_error, input%kappa_conv,&
                  input%method_updatekernel,input%purification_quickreturn, &
                  input%correction_co_contra, cdft, input%frag, ref_frags)
           else
              call getLocalizedBasis(iproc,nproc,at,KSwfn%orbs,rxyz,denspot,GPU,trace,trace_old,fnrm_tmb,&
                  info_basis_functions,nlpsp,input%lin%scf_mode,ldiis,input%SIC,tmb,energs, &
                  input%lin%nItPrecond,target_function,input%lin%correctionOrthoconstraint,&
                  nit_basis,&
                  ratio_deltas,orthonormalization_on,input%lin%extra_states,itout,conv_crit_TMB,input%experimental_mode,&
                  input%lin%early_stop, input%lin%gnrm_dynamic, input%lin%min_gnrm_for_dynamic, &
                  can_use_ham, norder_taylor, input%lin%max_inversion_error, input%kappa_conv,&
                  input%method_updatekernel,input%purification_quickreturn, &
                  input%correction_co_contra)
           end if
           !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
           reduce_conf=.true.
           if (iproc==0) then
               call yaml_sequence_close()
           end if

           !update weight matrix following basis optimization (could add check to ensure this is really necessary)
           if (input%lin%constrained_dft) then
              if (trim(cdft%method)=='fragment_density') then ! fragment density approach
                 if (input%lin%diag_start) stop 'fragment_density not allowed'
              else if (trim(cdft%method)=='lowdin') then ! direct weight matrix approach
                 !should already have overlap matrix so no need to recalculate
                 call calculate_weight_matrix_lowdin_wrapper(cdft,tmb,input,ref_frags,.false.,input%lin%order_taylor)
              else 
                 stop 'Error invalid method for calculating CDFT weight matrix'
              end if
           end if

           tmb%can_use_transposed=.false. !since basis functions have changed...

           if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) ldiis_coeff%alpha_coeff=input%lin%alphaSD_coeff !reset to default value

           ! I think this is causing a memory leak somehow in certain cases (possibly only with fragment calculations?)
           if (input%inputPsiId==101 .and. info_basis_functions<=-2 .and. itout==1 .and. (.not.input%lin%fragment_calculation)) then
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



      if (input%lin%constrained_dft) then
         call DIIS_set(30,valpha,1,1,vdiis)
         call vcopy(tmb%orbs%norb**2,tmb%coeff(1,1),1,coeff_tmp(1,1),1)
         vold=cdft%lag_mult
      end if
      ! CDFT: need to pass V*w_ab to get_coeff so that it can be added to H_ab and the correct KS eqn can therefore be solved
      ! CDFT: for the first iteration this will be some initial guess for V (or from the previous outer loop)
      ! CDFT: all this will be in some extra CDFT loop
      cdft_loop : do cdft_it=1,100
         if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION .and. input%lin%constrained_dft) then 
            call DIIS_free(ldiis_coeff)
            if (input%lin%extra_states==0) then
               call DIIS_set(ldiis_coeff_hist,0.1_gp,tmb%orbs%norb*KSwfn%orbs%norbp,1,ldiis_coeff)
            else
               call DIIS_set(ldiis_coeff_hist,0.1_gp,tmb%orbs%norb*tmb%orbs%norbp,1,ldiis_coeff)
            end if
         end if
         ! The self consistency cycle. Here we try to get a self consistent density/potential with the fixed basis.
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
         convCritMix = convCritMix_init*fnrm_tmb
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
             invert_overlap_matrix = (.not.(target_function==TARGET_FUNCTION_IS_HYBRID .and. &
                                       (input%method_updatekernel==UPDATE_BY_FOE .or. &
                                      input%method_updatekernel==UPDATE_BY_RENORMALIZATION)) .and. &
                                      it_scc==1)
                                      !cur_it_highaccuracy==1)
             !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
             if(update_phi .and. can_use_ham) then! .and. info_basis_functions>=0) then
                if (input%lin%constrained_dft) then
                   call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                        infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,calculate_overlap,invert_overlap_matrix,update_phi,&
                        .false.,input%lin%extra_states,itout,it_scc,cdft_it,norder_taylor,input%lin%max_inversion_error,&
                        input%purification_quickreturn,&
                        input%calculate_KS_residue,input%calculate_gap,&
                        convcrit_dmin,nitdmin,input%lin%curvefit_dmin,ldiis_coeff,reorder,cdft)
                else
                   call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                        infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,calculate_overlap,invert_overlap_matrix,update_phi,&
                        .false.,input%lin%extra_states,itout,it_scc,cdft_it,norder_taylor,input%lin%max_inversion_error,&
                        input%purification_quickreturn,&
                        input%calculate_KS_residue,input%calculate_gap,&
                        convcrit_dmin,nitdmin,input%lin%curvefit_dmin,ldiis_coeff,reorder)
                end if
             else
                if (input%lin%constrained_dft) then
                   call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                        infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,calculate_overlap,invert_overlap_matrix,update_phi,&
                        .true.,input%lin%extra_states,itout,it_scc,cdft_it,norder_taylor,input%lin%max_inversion_error,&
                        input%purification_quickreturn,&
                        input%calculate_KS_residue,input%calculate_gap,&
                        convcrit_dmin,nitdmin,input%lin%curvefit_dmin,ldiis_coeff,reorder,cdft)
                else
                   call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                        infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,calculate_overlap,invert_overlap_matrix,update_phi,&
                        .true.,input%lin%extra_states,itout,it_scc,cdft_it,norder_taylor,input%lin%max_inversion_error,&
                        input%purification_quickreturn,&
                        input%calculate_KS_residue,input%calculate_gap,&
                        convcrit_dmin,nitdmin,input%lin%curvefit_dmin,ldiis_coeff,reorder)
                end if
             end if
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
             !if(iproc==0) write(*,'(a,7es14.6)') 'energs', &
             !    energs%ebs,energs%eh,energs%exc,energs%evxc,energs%eexctX,energs%eion,energs%edisp
             energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
             energyDiff=energy-energyold
             energyold=energy

             ! update alpha_coeff for direct minimization steepest descents
             if(input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION .and. it_scc>1 .and.&
                  ldiis_coeff%idsx == 0 .and. (.not. input%lin%curvefit_dmin)) then
                ! apply a cap so that alpha_coeff never goes below around 1.d-2 or above 2
                if (energyDiff<0.d0 .and. ldiis_coeff%alpha_coeff < 1.8d0) then
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
                 call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
                 !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
                 !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
                 !!call clean_rho(iproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
             end if

             if (input%lin%constrained_dft) then
                !call timing(iproc,'constraineddft','ON')
                ! CDFT: see how satisfaction of constraint varies as kernel is updated
                ! CDFT: calculate Tr[Kw]-Nc
                weight_matrix_ = matrices_null()
                call allocate_matrices(tmb%linmat%m, allocate_full=.false., matname='weight_matrix_', mat=weight_matrix_)
                weight_matrix_%matrix_compr=cdft%weight_matrix_%matrix_compr

                !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
                call extract_taskgroup_inplace(tmb%linmat%m, weight_matrix_)
                call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%l,tmb%linmat%m, &
                     tmb%linmat%kernel_,weight_matrix_,&
                     ebs,tmb%coeff,KSwfn%orbs,tmb%orbs,.false.)
                !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
                call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%m, weight_matrix_)

                !tmb%linmat%denskern_large%matrix_compr = tmb%linmat%kernel_%matrix_compr
                call deallocate_matrices(weight_matrix_)
                !call timing(iproc,'constraineddft','OF')
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
                call check_negative_rho(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, &
                     denspot%rhov, rho_negative)
                if (rho_negative) then
                    call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
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
                      call mpiallred(pnrm_out, 1, mpi_sum, bigdft_mpi%mpi_comm)
                   end if
                   !pnrm_out = pnrm_out/dble(input%nspin)

                   pnrm_out=sqrt(pnrm_out)/(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*KSwfn%Lzd%Glr%d%n3i)
                   !only want to copy across when CDFT loop has also converged
                   if (.not. input%lin%constrained_dft .or. (ebs-cdft%charge < cdft_charge_thresh)) then
                      call vcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, &
                        denspot%rhov(1), 1, rhopotOld_out(1), 1)
                   end if
                end if

             end if

             ! Calculate the new potential.
             !!if(iproc==0) write(*,'(1x,a)') '---------------------------------------------------------------- Updating potential.'
             if (iproc==0) then
!                 if (iproc==0) call yaml_mapping_open('pot',flow=.true.)
                 !call yaml_map('update potential',.true.)
             end if
             if (iproc==0) call yaml_newline()
             

             call updatePotential(input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
             if (iproc==0) call yaml_mapping_close()


             ! update occupations wrt eigenvalues (NB for directmin these aren't guaranteed to be true eigenvalues)
             ! switch off for FOE at the moment
             ! switch off for directmin too, unless we decide to reactivate calculating the expectation values the output is meaningless
             if (input%lin%scf_mode/=LINEAR_FOE .and. input%lin%scf_mode/=LINEAR_DIRECT_MINIMIZATION) then
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
                      call mpiallred(pnrm_out, 1, mpi_sum, bigdft_mpi%mpi_comm)
                   end if

                   pnrm_out=sqrt(pnrm_out)/(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*KSwfn%Lzd%Glr%d%n3i)!)*input%nspin)
                   !only want to copy across when CDFT loop has also converged
                   if (.not. input%lin%constrained_dft .or. (ebs-cdft%charge < cdft_charge_thresh)) then
                      call vcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, &
                           denspot%rhov(1), 1, rhopotOld_out(1), 1) 
                   end if
                end if
             end if

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

         if (input%lin%constrained_dft) then

            !! CDFT: see how satisfaction of constraint varies as kernel is updated
            !! CDFT: calculate Tr[Kw]-Nc
            !call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%denskern,cdft%weight_matrix,&
            !     ebs,tmb%coeff,KSwfn%orbs,tmb%orbs,.false.)

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
            if (abs(vgrad) < cdft_charge_thresh) then
               !call vcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, &
               !     denspot%rhov(1), 1, rhopotOld_out(1), 1) 
               exit
            end if
            !reset to best previous coeffs, not necessarily original coeffs
            !if (abs(vgrad) < best_charge_diff) call vcopy(tmb%orbs%norb**2,tmb%coeff(1,1),1,coeff_tmp(1,1),1)
            !best_charge_diff=min(best_charge_diff,abs(vgrad))

            call timing(iproc,'constraineddft','ON')
            !! CHECK HERE WHETHER n3d is correct!!
            ! reset rhopotold (to zero) to ensure we don't exit immediately if V only changes a little
            !call to_zero(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, rhopotOld(1)) 

            ! assuming density mixing/no mixing not potential mixing
            call vcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, &
                 rhopotOld_out(1), 1, rhopotOld(1), 1) 

            call vcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, &
                 rhopotOld_out(1), 1, denspot%rhov(1), 1)
            call timing(iproc,'constraineddft','OF')

            call updatePotential(input%nspin,denspot,energs%eh,energs%exc,energs%evxc)

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

         ! if not constrained DFT exit straight away
         else
            exit
         end if
      end do cdft_loop
      if (input%lin%constrained_dft) call DIIS_free(vdiis)
      ! CDFT: end of CDFT loop to find V which correctly imposes constraint and corresponding density

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



  if (input%write_orbitals) then
      call build_ks_orbitals(iproc, nproc, tmb, KSwfn, at, rxyz, denspot, GPU, &
               energs, nlpsp, input, norder_taylor,&
               energy, energyDiff, energyold)
  end if


  ! Diagonalize the matrix for the FOE/direct min case to get the coefficients. Only necessary if
  ! the Pulay forces are to be calculated, or if we are printing eigenvalues for restart
  if ((input%lin%scf_mode==LINEAR_FOE.or.input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION)& 
       .and. (input%lin%pulay_correction.or.input%lin%plotBasisFunctions /= WF_FORMAT_NONE&
       .or. input%lin%diag_end)) then

       !!if (input%lin%scf_mode==LINEAR_FOE) then
       !!    tmb%coeff=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%coeff')
       !!end if

       !!call extract_taskgroup_inplace(tmb%linmat%l, tmb%linmat%kernel_)
       call get_coeff(iproc,nproc,LINEAR_MIXDENS_SIMPLE,KSwfn%orbs,at,rxyz,denspot,GPU,&
           infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,update_phi,.true.,.false.,&
           .true.,input%lin%extra_states,itout,0,0,norder_taylor,input%lin%max_inversion_error,&
           input%purification_quickreturn,&
           input%calculate_KS_residue,input%calculate_gap)
       !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)

       !!if (input%lin%scf_mode==LINEAR_FOE) then
       !!    call f_free_ptr(tmb%coeff)
       !!end if

       if (bigdft_mpi%iproc ==0) then 
          call write_eigenvalues_data(0.1d0,tmb%orbs,mom_vec_fake)
       end if
  end if

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
           input%purification_quickreturn,&
           input%calculate_KS_residue,input%calculate_gap)
      !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
      !!call scalprod_on_boundary(iproc, nproc, tmb, kswfn%orbs, at, fpulay)
      call pulay_correction_new(iproc, nproc, tmb, kswfn%orbs, at, fpulay)
      !!if (input%lin%scf_mode==LINEAR_FOE) then
      !!    call f_free_ptr(tmb%coeff)
      !!end if
  end if
  !! END TEST ######################


  if (input%lin%fragment_calculation .and. input%frag%nfrag>1) then
     call coeff_weight_analysis(iproc, nproc, input, KSwfn%orbs, tmb, ref_frags)
  end if


  if (input%lin%constrained_dft) then
     call cdft_data_free(cdft)
     call f_free(coeff_tmp)
  end if


  ! print the final summary
  call print_info(.true.)



  if (iproc==0) call yaml_sequence_close()

  if (input%loewdin_charge_analysis) then
      call loewdin_charge_analysis(iproc, tmb, at, denspot, calculate_overlap_matrix=.true., &
           calculate_ovrlp_half=.true., meth_overlap=0)
      call support_function_multipoles(iproc, tmb, at, denspot)
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
      call to_zero(3*at%astruct%nat, fpulay(1,1))
  end if

  if(tmb%ham_descr%can_use_transposed) then
      call f_free_ptr(tmb%ham_descr%psit_c)
      call f_free_ptr(tmb%ham_descr%psit_f)
      tmb%ham_descr%can_use_transposed=.false.
  end if
  ! here or cluster, not sure which is best
  deallocate(tmb%confdatarr, stat=istat)


  !Write the linear wavefunctions to file if asked, also write Hamiltonian and overlap matrices
  if (input%lin%plotBasisFunctions /= WF_FORMAT_NONE) then
     nelec=0
     do iat=1,at%astruct%nat
        ityp=at%astruct%iatype(iat)
        nelec=nelec+at%nelpsp(ityp)
     enddo
     call writemywaves_linear(iproc,trim(input%dir_output) // 'minBasis',input%lin%plotBasisFunctions,&
          max(tmb%npsidim_orbs,tmb%npsidim_comp),tmb%Lzd,tmb%orbs,nelec,at,rxyz,tmb%psi,tmb%coeff)
     call write_linear_matrices(iproc,nproc,input%imethod_overlap,trim(input%dir_output),&
          input%lin%plotBasisFunctions,tmb,at,rxyz)
  end if


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
      call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
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

    subroutine allocate_local_arrays()

      locrad = f_malloc(tmb%lzd%nlr,id='locrad')
      ! Allocate the old charge density (used to calculate the variation in the charge density)
      rhopotold_out = f_malloc(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim),id='rhopotold_out')

    end subroutine allocate_local_arrays


    subroutine deallocate_local_arrays()

      call f_free(locrad)
      call f_free(rhopotold_out)

    end subroutine deallocate_local_arrays


    subroutine check_inputguess()
      real(kind=8) :: dnrm2
      if (input%inputPsiId==101) then           !should we put 102 also?

          if (input%lin%pulay_correction) then
             ! Check the input guess by calculation the Pulay forces.

             call to_zero(tmb%ham_descr%npsidim_orbs,tmb%ham_descr%psi(1))
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
              call to_zero(3*at%astruct%nat, fpulay(1,1))
              if (input%lin%scf_mode==LINEAR_FOE) then
                 nit_lowaccuracy=input%lin%nit_lowaccuracy
              else
                 nit_lowaccuracy=0
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
          else
              call yaml_map('kernel method','DIAG')
          end if

          if (input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE .or.  input%lin%scf_mode==LINEAR_FOE &
              .or. input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
              call yaml_map('mix entity','DENS')
          else if (input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
              call yaml_map('mix entity','POT')
          end if
          call yaml_map('mix hist',mix_hist)
          call yaml_map('conv crit',convCritMix,fmt='(es8.2)')


          if (input%lin%constrained_dft) then
             if (iproc==0) then
                 call yaml_newline()
                 call yaml_map('iter',it_scc,fmt='(i6)')
                 call yaml_map('delta',pnrm,fmt='(es9.2)')
                 call yaml_map('energy',energy,fmt='(es24.17)')
                 call yaml_map('D',energyDiff,fmt='(es10.3)')
                 call yaml_map('Tr(KW)',ebs,fmt='(es14.4)')
                 call yaml_mapping_close()
             end if
          else
             if (iproc==0) then
                 call yaml_newline()
                 call yaml_map('iter',it_scc,fmt='(i6)')
                 call yaml_map('delta',pnrm,fmt='(es9.2)')
                 call yaml_map('energy',energy,fmt='(es24.17)')
                 call yaml_map('D',energyDiff,fmt='(es10.3)')
                 call yaml_mapping_close()
             end if
          end if     
          call yaml_sequence_close()
      end if




    end subroutine printSummary

    !> Print a short summary of some values calculated during the last iteration in the self
    !! consistency cycle.
    subroutine print_info(final)
      implicit none

      real(kind=8) :: energyDiff, mean_conf
      logical, intent(in) :: final

      energyDiff = energy - energyoldout

      mean_conf=0.d0
      do iorb=1,tmb%orbs%norbp
          mean_conf=mean_conf+tmb%confdatarr(iorb)%prefac
      end do

      if (nproc > 1) then
         call mpiallred(mean_conf, 1, mpi_sum, bigdft_mpi%mpi_comm)
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

    ! Determine whether the sign of the energy change is teh same as in the previous iteration
    ! (i.e. whether the energy continues to increase or decrease)
    tt = sign(energyDiff,sign_of_energy_change)
    if (tt/energyDiff>0.d0) then
        ! same sign, everything ok
    else if (abs(energyDiff/energy)>1.d-7) then
        nit_energyoscillation = nit_energyoscillation + 1
        if (iproc==0) then
            call yaml_warning('oscillation of the energy, increase counter')
            call yaml_map('energy_scillation_counter',nit_energyoscillation)
        end if
    end if
    if (nit_energyoscillation>1) then
        nit_scc = nit_scc + 1
        nit_energyoscillation = 0
        if (iproc==0) call yaml_map('new nit_scc',nit_scc)
    end if

    ! Determine the sign of the energy change
    sign_of_energy_change = sign(1.d0,energyDiff)


    end subroutine print_info


    subroutine intermediate_forces()

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
          call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
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

      call updatePotential(input%nspin,denspot,energs%eh,energs%exc,energs%evxc)

      ! Density already present in denspot%rho_work
      call vcopy(denspot%dpbox%ndimpot,denspot%rho_work(1),1,denspot%pot_work(1),1)
      call H_potential('D',denspot%pkernel,denspot%pot_work,denspot%pot_work,ehart_fake,&
           0.0_dp,.false.,stress_tensor=hstrten)

      
      KSwfn%psi=f_malloc_ptr(1,id='KSwfn%psi')


      fpulay=0.d0
      call calculate_forces(iproc,nproc,denspot%pkernel%mpi_env%nproc,KSwfn%Lzd%Glr,at,KSwfn%orbs,nlpsp,rxyz,& 
           KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
           denspot%dpbox%i3s+denspot%dpbox%i3xcsh,denspot%dpbox%n3p,&
           denspot%dpbox%nrhodim,.false.,denspot%dpbox%ngatherarr,denspot%rho_work,&
           denspot%pot_work,denspot%V_XC,size(KSwfn%psi),KSwfn%psi,fion,fdisp,fxyz,&
           ewaldstr,hstrten,xcstr,strten,fnoise,pressure,denspot%psoffset,1,tmb,fpulay)
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

end subroutine linearScaling




subroutine output_fragment_rotations(iproc,nat,rxyz,iformat,filename,input_frag,ref_frags)
  use module_base
  use module_types
  use yaml_output
  use module_fragments
  use internal_io
  use module_interfaces
  implicit none

  integer, intent(in) :: iproc, iformat, nat
  character(len=*), intent(in) :: filename
  real(gp), dimension(3,nat), intent(in) :: rxyz
  type(fragmentInputParameters), intent(in) :: input_frag
  type(system_fragment), dimension(input_frag%nfrag_ref), intent(in) :: ref_frags
  !Local variables
  integer :: ifrag, jfrag, ifrag_ref, jfrag_ref, iat, isfat, jsfat
  real(kind=gp), dimension(:,:), allocatable :: rxyz_ref, rxyz_new
  real(kind=gp) :: null_axe
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

           call find_frag_trans(ref_frags(ifrag_ref)%astruct_frg%nat,rxyz_ref,rxyz_new,frag_trans)

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

   end if

end subroutine output_fragment_rotations
