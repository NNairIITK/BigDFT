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
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(atoms_data),intent(inout) :: at
  type(input_variables),intent(in) :: input ! need to hack to be inout for geopt changes
  real(8),dimension(3,at%astruct%nat),intent(inout) :: rxyz
  real(8),dimension(3,at%astruct%nat),intent(out) :: fpulay
  type(DFT_local_fields), intent(inout) :: denspot
  real(gp), dimension(*), intent(inout) :: rhopotold
  type(DFT_PSP_projectors),intent(inout) :: nlpsp
  type(GPU_pointers),intent(inout) :: GPU
  type(energy_terms),intent(inout) :: energs
  real(gp), dimension(:), pointer :: rho,pot
  real(8),intent(out) :: energy
  type(DFT_wavefunction),intent(inout),target :: tmb
  type(DFT_wavefunction),intent(inout),target :: KSwfn
  integer,intent(out) :: infocode
  type(system_fragment), dimension(:), pointer :: ref_frags ! for transfer integrals
  type(cdft_data), intent(inout) :: cdft
  real(kind=8),dimension(3,at%astruct%nat),intent(in) :: fdisp, fion
  
  real(8) :: pnrm,trace,trace_old,fnrm_tmb
  integer :: infoCoeff,istat,iall,it_scc,itout,info_scf,i,ierr,iorb
  character(len=*), parameter :: subname='linearScaling'
  real(8),dimension(:),allocatable :: rhopotold_out
  real(8) :: energyold, energyDiff, energyoldout, fnrm_pulay, convCritMix
  type(mixrhopotDIISParameters) :: mixdiis
  type(localizedDIISParameters) :: ldiis!, ldiis_coeff
  type(DIIS_obj) :: ldiis_coeff, vdiis
  logical :: can_use_ham, update_phi, locreg_increased, reduce_conf, orthonormalization_on
  logical :: fix_support_functions, check_initialguess
  integer :: itype, istart, nit_lowaccuracy, nit_highaccuracy
  real(8),dimension(:),allocatable :: locrad_tmp
  integer :: ldiis_coeff_hist, nitdmin
  logical :: ldiis_coeff_changed
  integer :: mix_hist, info_basis_functions, nit_scc, cur_it_highaccuracy
  real(8) :: pnrm_out, alpha_mix, ratio_deltas, convcrit_dmin
  logical :: lowaccur_converged, exit_outer_loop
  real(8),dimension(:),allocatable :: locrad
  integer:: target_function, nit_basis
  type(sparse_matrix) :: ham_small
  integer :: isegsmall, iseglarge, iismall, iilarge, is, ie
  integer :: matrixindex_in_compressed
  
  real(kind=gp) :: ebs, vgrad_old, vgrad, valpha, vold, vgrad2, vold_tmp, conv_crit_TMB
  real(kind=gp), allocatable, dimension(:,:) :: coeff_tmp
  integer :: ind_denskern, ind_ham, jorb, cdft_it, nelec, iat, ityp, ifrag, ifrag_charged, ifrag_ref, isforb, itmb
  integer :: dmin_diag_it, dmin_diag_freq, ioffset
  logical :: reorder, rho_negative
  real(wp), dimension(:,:,:), pointer :: mom_vec_fake

  !!! EXPERIMENTAL ############################################
  type(sparse_matrix) :: denskern_init
  real(8),dimension(:),allocatable :: rho_init, rho_init_old, philarge
  real(8) :: tt, ddot, tt_old, meanconf_der, weight_boundary, weight_tot
  integer :: idens_cons, ii, sdim, ldim, npsidim_large, ists, istl, nspin, unitname, ilr
  real(8),dimension(10000) :: meanconf_array
  character(len=5) :: num
  character(len=50) :: filename
  real(kind=8),dimension(:,:),allocatable :: phi_delta
  !!! #########################################################

  ! DEBUG - for calculating centres
  type(workarr_sumrho) :: w
  real(gp), allocatable, dimension(:,:,:,:) :: psir
  integer :: ind, i_all, i_stat, nspinor, ix, iy, iz, iix, iiy, iiz
  real(gp) :: psix, psiy, psiz, xcent, ycent, zcent

  character(len=12) :: orbname
  real(gp), allocatable, dimension(:) :: psi2, gpsi, gpsi2
  real(gp), allocatable, dimension(:,:,:,:) :: psir2
  real(gp) :: tmb_diff, max_tmb_diff, cut
  integer :: j, k, n1i, n2i, n3i, i1, i2, i3, num_points, num_points_tot

  integer :: ist, iiorb, ncount
  real(kind=8) :: fnoise, pressure, ehart_fake, dnrm2
  real(kind=8),dimension(:,:),allocatable :: fxyz

  type(matrices) :: weight_matrix_

  call timing(iproc,'linscalinit','ON') !lr408t

  call f_routine(id='linear_scaling')

  call allocate_local_arrays()

  !!if(iproc==0) then
  !!    write(*,'(1x,a)') repeat('*',84)
  !!    write(*,'(1x,a)') '****************************** LINEAR SCALING VERSION ******************************'
  !!end if

  ! Allocate the communications buffers needed for the communications of the potential and
  ! post the messages. This will send to each process the part of the potential that this process
  ! needs for the application of the Hamlitonian to all orbitals on that process.
  call allocate_p2pComms_buffer(tmb%comgp)

  ! Initialize the DIIS mixing of the potential if required.
  if(input%lin%mixHist_lowaccuracy>0) then
      call initializeMixrhopotDIIS(input%lin%mixHist_lowaccuracy, denspot%dpbox%ndimrhopot, mixdiis)
  end if

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
  check_initialguess=.true.
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

  cut=maxval(tmb%lzd%llr(:)%locrad)

  !call nullify_sparse_matrix(ham_small) ! nullify anyway
  ham_small=sparse_matrix_null()

  if (input%lin%scf_mode==LINEAR_FOE) then ! allocate ham_small
     call sparse_copy_pattern(tmb%linmat%s,ham_small,iproc,subname)
     ham_small%matrix_compr = sparsematrix_malloc_ptr(ham_small,iaction=SPARSE_FULL,id='ham_small%matrix_compr')
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
  nullify(tmb%psit_c)
  nullify(tmb%psit_f)

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

  if (iproc==0) call yaml_open_map('Checking Communications of Minimal Basis')
  call check_communications_locreg(iproc,nproc,tmb%orbs,tmb%Lzd,tmb%collcom, &
       tmb%npsidim_orbs,tmb%npsidim_comp)
  if (iproc==0) call yaml_close_map()

  if (iproc==0) call yaml_open_map('Checking Communications of Enlarged Minimal Basis')
  call check_communications_locreg(iproc,nproc,tmb%orbs,tmb%ham_descr%lzd,tmb%ham_descr%collcom, &
       tmb%ham_descr%npsidim_orbs,tmb%ham_descr%npsidim_comp)
  if (iproc ==0) call yaml_close_map()


  ! CDFT: calculate w_ab here given w(r)
  ! CDFT: first check that we aren't updating the basis at any point and we don't have any low acc iterations
  if (input%lin%constrained_dft) then
     if (nit_lowaccuracy>0 .or. input%lin%nItBasis_highaccuracy>1) then
        stop 'Basis cannot be updated for now in constrained DFT calculations and no low accuracy is allowed'
     end if

     weight_matrix_ = matrices_null()
     call allocate_matrices(tmb%linmat%m, allocate_full=.false., matname='weight_matrix_', mat=weight_matrix_)
     weight_matrix_%matrix_compr=cdft%weight_matrix%matrix_compr
     call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%l,tmb%linmat%m, &
          tmb%linmat%kernel_,weight_matrix_,&
          ebs,tmb%coeff,KSwfn%orbs,tmb%orbs,.false.)
     !tmb%linmat%denskern_large%matrix_compr = tmb%linmat%kernel_%matrix_compr
     call deallocate_matrices(weight_matrix_)

     call timing(iproc,'constraineddft','ON')
     vgrad_old=ebs-cdft%charge

     !!if (iproc==0) write(*,'(a,4(ES16.6e3,2x))') 'N, Tr(KW), Tr(KW)-N, V*(Tr(KW)-N)',&
     !!     cdft%charge,ebs,vgrad_old,cdft%lag_mult*vgrad_old
     if (iproc==0) then
         call yaml_open_map('CDFT infos')
         call yaml_map('N',cdft%charge,fmt='(es16.6e3)')
         call yaml_map('Tr(KW)',ebs,fmt='(es16.6e3)')
         call yaml_map('Tr(KW)-N',vgrad_old,fmt='(es16.6e3)')
         call yaml_map('V*(Tr(KW)-N)',cdft%lag_mult*vgrad_old,fmt='(es16.6e3)')
         call yaml_close_map()
     end if
     vgrad_old=abs(vgrad_old)
     valpha=0.5_gp

     coeff_tmp=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='coeff_tmp')
     call timing(iproc,'constraineddft','OF')
  end if

  !!! EXPERIMENTAL #######################
  !!denskern_init=tmb%linmat%denskern
  !!nullify(denskern_init%matrix_compr)
  !!!nullify(denskern_init%matrix)
  !!allocate(denskern_init%matrix_compr(size(tmb%linmat%denskern%matrix_compr)))
  !!!allocate(denskern_init%matrix(size(tmb%linmat%denskern%matrix)))
  !!call vcopy(size(tmb%linmat%denskern%matrix_compr), tmb%linmat%denskern%matrix_compr, 1, denskern_init%matrix_compr, 1)
  !!!call vcopy(size(tmb%linmat%denskern%matrix), tmb%linmat%denskern%matrix, 1, denskern_init%matrix, 1)
  !!allocate(rho_init(size(denspot%rhov)))
  !!allocate(rho_init_old(size(denspot%rhov)))
  !!tt_old=1.d100
  !!rho_init=0.d0
  !!rho_init_old=0.d0
  !!idens_cons=0
  !!! ####################################

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
     call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
          tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, &
          denspot%rhov, rho_negative)
     if (rho_negative) then
         call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
         !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
         !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
         !!call clean_rho(iproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
     end if
     ! Calculate the new potential.
     !if(iproc==0) write(*,'(1x,a)') '---------------------------------------------------------------- Updating potential.'
     !if (iproc==0) call yaml_map('update potential',.true.)
     if (iproc==0) call yaml_open_map('update pot',flow=.true.)
     call updatePotential(input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
  end if

  call timing(iproc,'linscalinit','OF') !lr408t

  !! DEBUG - check centres
  !ind=1
  !do iorb=1,tmb%orbs%norbp
  !   iat=tmb%orbs%onwhichatom(iorb+tmb%orbs%isorb)
  !   ilr=tmb%orbs%inwhichlocreg(iorb+tmb%orbs%isorb)
  !
  !   allocate(psir(tmb%lzd%llr(ilr)%d%n1i, tmb%lzd%llr(ilr)%d%n2i, tmb%lzd%llr(ilr)%d%n3i, 1+ndebug),stat=i_stat)
  !   call memocc(i_stat,psir,'psir',subname)
  !   call initialize_work_arrays_sumrho(tmb%lzd%llr(ilr),w)
  !
  !   call daub_to_isf(tmb%lzd%llr(ilr),w,tmb%psi(ind),psir)
  !
  !   xcent=0.0d0
  !   ycent=0.0d0
  !   zcent=0.0d0
  !   do iz=1,tmb%lzd%llr(ilr)%d%n3i
  !      iiz=iz-15+tmb%lzd%llr(ilr)%nsi3
  !      do iy=1,tmb%lzd%llr(ilr)%d%n2i
  !         iiy=iy-15+tmb%lzd%llr(ilr)%nsi2
  !         do ix=1,tmb%lzd%llr(ilr)%d%n1i
  !            iix=ix-15+tmb%lzd%llr(ilr)%nsi1
  !            psix=psir(ix,iy,iz,1)*(iix*tmb%lzd%hgrids(1)*0.5d0)
  !            psiy=psir(ix,iy,iz,1)*(iiy*tmb%lzd%hgrids(2)*0.5d0)
  !            psiz=psir(ix,iy,iz,1)*(iiz*tmb%lzd%hgrids(3)*0.5d0)
  !            xcent=xcent+psir(ix,iy,iz,1)*psix
  !            ycent=ycent+psir(ix,iy,iz,1)*psiy
  !            zcent=zcent+psir(ix,iy,iz,1)*psiz
  !         end do
  !      end do
  !   end do
  !
  !   write(*,'(a,4I4,3(F12.8,x),3(F8.4,x))') 'iproc,iorb,ilr,iat,(xcent,ycent,zcent)-locregcenter,xcent,ycent,zcent',&
  !        iproc,iorb+tmb%orbs%isorb,ilr,iat,xcent-tmb%lzd%llr(ilr)%locregcenter(1),&
  !        ycent-tmb%lzd%llr(ilr)%locregcenter(2),zcent-tmb%lzd%llr(ilr)%locregcenter(3),&
  !        xcent,ycent,zcent
  !
  !   ind=ind+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
  !   call deallocate_work_arrays_sumrho(w)
  !   i_all=-product(shape(psir))*kind(psir)
  !   deallocate(psir,stat=i_stat)
  !   call memocc(i_stat,i_all,'psir',subname)
  !end do
  !! END DEBUG - check centres

  ! Add one iteration if no low accuracy is desired since we need then a first fake iteration, with istart=0
  istart = min(1,nit_lowaccuracy)
  infocode=0 !default value
  ! This is the main outer loop. Each iteration of this loop consists of a first loop in which the basis functions
  ! are optimized and a consecutive loop in which the density is mixed.
  if (iproc==0) then
      call yaml_comment('Self-Consistent Cycle',hfill='-')
      call yaml_open_sequence('Ground State Optimization')
  end if
  outerLoop: do itout=istart,nit_lowaccuracy+nit_highaccuracy

      if (input%lin%nlevel_accuracy==2) then
          ! Check whether the low accuracy part (i.e. with strong confining potential) has converged.
          call check_whether_lowaccuracy_converged(itout, nit_lowaccuracy, input%lin%lowaccuracy_conv_crit, &
               lowaccur_converged, pnrm_out)
          ! Set all remaining variables that we need for the optimizations of the basis functions and the mixing.
          call set_optimization_variables(input, at, tmb%orbs, tmb%lzd%nlr, tmb%orbs%onwhichatom, tmb%confdatarr, &
               convCritMix, lowaccur_converged, nit_scc, mix_hist, alpha_mix, locrad, target_function, nit_basis, &
               convcrit_dmin, nitdmin, conv_crit_TMB)
      else if (input%lin%nlevel_accuracy==1 .and. itout==1) then
          call set_variables_for_hybrid(tmb%lzd%nlr, input, at, tmb%orbs, lowaccur_converged, tmb%confdatarr, &
               target_function, nit_basis, nit_scc, mix_hist, locrad, alpha_mix, convCritMix, conv_crit_TMB)
               convcrit_dmin=input%lin%convCritDmin_highaccuracy
               nitdmin=input%lin%nItdmin_highaccuracy

         !! lowaccur_converged=.false.
         !! do iorb=1,tmb%orbs%norbp
         !!     ilr=tmb%orbs%inwhichlocreg(tmb%orbs%isorb+iorb)
         !!     iiat=tmb%orbs%onwhichatom(tmb%orbs%isorb+iorb)
         !!     tmb%confdatarr(iorb)%prefac=input%lin%potentialPrefac_lowaccuracy(at%astruct%iatype(iiat))
         !! end do
         !! target_function=TARGET_FUNCTION_IS_HYBRID
         !! nit_basis=input%lin%nItBasis_lowaccuracy
         !! nit_scc=input%lin%nitSCCWhenFixed_lowaccuracy
         !! mix_hist=input%lin%mixHist_lowaccuracy
         !! do ilr=1,tmb%lzd%nlr
         !!     locrad(ilr)=input%lin%locrad_lowaccuracy(ilr)
         !! end do
         !! alpha_mix=input%lin%alpha_mix_lowaccuracy
         !! convCritMix=input%lin%convCritMix_lowaccuracy
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

          if (locreg_increased .and. input%lin%scf_mode==LINEAR_FOE) then ! deallocate ham_small
             call deallocate_sparse_matrix(ham_small,subname)
             !call nullify_sparse_matrix(ham_small)
             ham_small=sparse_matrix_null()
             call sparse_copy_pattern(tmb%linmat%s,ham_small,iproc,subname)
             ham_small%matrix_compr = sparsematrix_malloc_ptr(ham_small,iaction=SPARSE_FULL,id='ham_small%matrix_compr')
          end if

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
                call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                     infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,update_phi,update_phi,&
                     .true.,ham_small,input%lin%extra_states,itout,0,0,input%lin%order_taylor,&
                     input%purification_quickreturn,input%adjust_FOE_temperature,&
                     input%calculate_KS_residue,input%calculate_gap,&
                     convcrit_dmin,nitdmin,input%lin%curvefit_dmin,ldiis_coeff)
             end if
          end if

          ! Some special treatement if we are in the high accuracy part
          call adjust_DIIS_for_high_accuracy(input, denspot, mixdiis, lowaccur_converged, &
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
          if (associated(tmb%psit_c)) then
              call f_free_ptr(tmb%psit_c)
          end if
          if (associated(tmb%psit_f)) then
              call f_free_ptr(tmb%psit_f)
          end if
          tmb%can_use_transposed=.false.
          if (iproc==0) then
              call yaml_sequence(advance='no')
              call yaml_open_sequence('fake iteration',label=&
                                 'it_fake'//trim(adjustl(yaml_toa(0,fmt='(i3.3)'))))
              call yaml_sequence(label='final_fake'//trim(adjustl(yaml_toa(0,fmt='(i3.3)'))),advance='no')
              call yaml_open_map(flow=.true.)
              call yaml_map('fake iteration','bridge low accuracy')
              call yaml_close_map
              call yaml_close_sequence()
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
               call yaml_open_sequence('support function optimization',label=&
                              'it_supfun'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)'))))
           end if
           !!if (itout<=2) then
               call getLocalizedBasis(iproc,nproc,at,KSwfn%orbs,rxyz,denspot,GPU,trace,trace_old,fnrm_tmb,&
                   info_basis_functions,nlpsp,input%lin%scf_mode,ldiis,input%SIC,tmb,energs, &
                   input%lin%nItPrecond,target_function,input%lin%correctionOrthoconstraint,&
                   nit_basis,&
                   ratio_deltas,orthonormalization_on,input%lin%extra_states,itout,conv_crit_TMB,input%experimental_mode,&
                   input%lin%early_stop, input%lin%gnrm_dynamic, input%lin%min_gnrm_for_dynamic, &
                   can_use_ham, input%lin%order_taylor, input%kappa_conv,&
                   input%method_updatekernel,input%purification_quickreturn, input%adjust_FOE_temperature, &
                   input%correction_co_contra)
               reduce_conf=.true.
           !!else
           !!    cut=cut-0.5d0
           !!    if (iproc==0) write(*,'(a,f7.2)') 'new cutoff:', cut
           !!    call cut_at_boundaries(cut, tmb)
           !!    ist=1
           !!    do iorb=1,tmb%orbs%norbp
           !!        iiorb=tmb%orbs%isorb+iorb
           !!        ilr=tmb%orbs%inwhichlocreg(iiorb)
           !!        ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
           !!        tt=dnrm2(ncount, tmb%psi(ist), 1, tmb)
           !!        tt=1/tt
           !!        !call dscal(ncount, tt, tmb%psi(ist), 1)
           !!        tt=dnrm2(ncount, tmb%psi(ist), 1, tmb)
           !!        write(*,*) 'iiorb, tt', iiorb, tt
           !!        ist=ist+ncount
           !!    end do
           !!end if
           if (iproc==0) then
               call yaml_close_sequence()
           end if

           !!! WRITE SUPPORT FUNCTIONS TO DISK ############################################
           !!npsidim_large=tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f                                                 
           !!allocate(philarge((tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f)*tmb%orbs%norbp))                          
           !!philarge=0.d0
           !!ists=1                                                                                                          
           !!istl=1
           !!do iorb=1,tmb%orbs%norbp
           !!    ilr = tmb%orbs%inWhichLocreg(tmb%orbs%isorb+iorb)                                                           
           !!    sdim=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f                                            
           !!    ldim=tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f                                                      
           !!    nspin=1 !this must be modified later
           !!    call Lpsi_to_global2(iproc, sdim, ldim, tmb%orbs%norb, tmb%orbs%nspinor, nspin, tmb%lzd%glr, &              
           !!         tmb%lzd%llr(ilr), tmb%psi(ists), philarge(istl))                                                       
           !!    write(num,'(i5.5)') tmb%orbs%isorb+iorb
           !!    filename='supfun_'//num
           !!    unitname=100*iproc+5
           !!    open(unit=unitname,file=trim(filename))
           !!    do i=1,tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
           !!        write(unitname,'(es25.17)') philarge(istl+i-1)
           !!    end do
           !!    close(unit=unitname)
           !!    ists=ists+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f                                       
           !!    istl=istl+tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f                                                 
           !!end do
           !!deallocate(philarge)
           !!! ############################################################################

           tmb%can_use_transposed=.false. !since basis functions have changed...

           if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) ldiis_coeff%alpha_coeff=input%lin%alphaSD_coeff !reset to default value

           if (input%inputPsiId==101 .and. info_basis_functions<=-2 .and. itout==1) then
               ! There seem to be some convergence problems after a restart. Better to quit
               ! and start with a new AO input guess.
               if (iproc==0) write(*,'(1x,a)') 'There are convergence problems after the restart. &
                                                &Start over again with an AO input guess.'
               if (associated(tmb%psit_c)) then
                   call f_free_ptr(tmb%psit_c)
               end if
               if (associated(tmb%psit_f)) then
                   call f_free_ptr(tmb%psit_f)
               end if
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

      if (can_use_ham .and. input%lin%scf_mode==LINEAR_FOE) then ! copy ham to ham_small here already as it won't be changing
        ! NOT ENTIRELY GENERAL HERE - assuming ovrlp is small and ham is large, converting ham to match ovrlp

         call timing(iproc,'FOE_init','ON') !lr408t

         iismall=0
         iseglarge=1
         do isegsmall=1,tmb%linmat%s%nseg
            do
               is=max(tmb%linmat%s%keyg(1,isegsmall),tmb%linmat%m%keyg(1,iseglarge))
               ie=min(tmb%linmat%s%keyg(2,isegsmall),tmb%linmat%m%keyg(2,iseglarge))
               iilarge=tmb%linmat%m%keyv(iseglarge)-tmb%linmat%m%keyg(1,iseglarge)
               do i=is,ie
                  iismall=iismall+1
                  ham_small%matrix_compr(iismall)=tmb%linmat%ham_%matrix_compr(iilarge+i)
               end do
               if (ie>=is) exit
               iseglarge=iseglarge+1
            end do
         end do

         call timing(iproc,'FOE_init','OF') !lr408t

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
             !    call yaml_open_map('kernel optimization',label=&
             !         'it_kernel'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//&
             !         '_'//trim(adjustl(yaml_toa(cdft_it,fmt='(i3.3)'))))
             !else
             !    call yaml_open_map('kernel optimization',label=&
             !         'it_kernel'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)'))))
             !end if
             call yaml_sequence(advance='no')
             if (input%lin%constrained_dft) then
                 call yaml_open_sequence('kernel optimization',label=&
                      'it_kernel'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//&
                      '_'//trim(adjustl(yaml_toa(cdft_it,fmt='(i3.3)'))))
             else
                 call yaml_open_sequence('kernel optimization',label=&
                      'it_kernel'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)'))))
             end if
         end if
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
                !call yaml_open_map(flow=.false.)
                call yaml_comment('kernel iter:'//yaml_toa(it_scc,fmt='(i6)'),hfill='-')
             end if
             if(update_phi .and. can_use_ham) then! .and. info_basis_functions>=0) then
                !!! TEST ###############################################################
                !!phi_delta=f_malloc0((/tmb%npsidim_orbs,3/),id='phi_delta')
                !!! Get the values of the support functions on the boundary of the localization region
                !!call extract_boundary(tmb, phi_delta, num_points, num_points_tot)
                !!weight_boundary=ddot(3*tmb%npsidim_orbs, phi_delta(1,1), 1, phi_delta(1,1), 1)
                !!call mpiallred(weight_boundary, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
                !!weight_boundary=sqrt(weight_boundary/tmb%orbs%norb)
                !!weight_tot=ddot(tmb%npsidim_orbs, tmb%psi(1), 1, tmb%psi(1), 1)
                !!call mpiallred(weight_tot, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
                !!weight_tot=sqrt(weight_tot/tmb%orbs%norb)
                !!call mpiallred(num_points, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
                !!call mpiallred(num_points_tot, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
                !!if (iproc==0) write(*,'(a,3es12.4,2I10)') 'weight boundary, weight tot, ratio, num points', &
                !!    weight_boundary, weight_tot, weight_boundary/weight_tot, num_points, num_points_tot
                !!call f_free(phi_delta)
                !!! END TEST ###########################################################
                if (input%lin%constrained_dft) then
                   call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                        infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,update_phi,update_phi,&
                        .false.,ham_small,input%lin%extra_states,itout,it_scc,cdft_it,input%lin%order_taylor,&
                        input%purification_quickreturn,input%adjust_FOE_temperature,&
                        input%calculate_KS_residue,input%calculate_gap,&
                        convcrit_dmin,nitdmin,input%lin%curvefit_dmin,ldiis_coeff,reorder,cdft)
                else
                   call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                        infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,update_phi,update_phi,&
                        .false.,ham_small,input%lin%extra_states,itout,it_scc,cdft_it,input%lin%order_taylor,&
                        input%purification_quickreturn,input%adjust_FOE_temperature,&
                        input%calculate_KS_residue,input%calculate_gap,&
                        convcrit_dmin,nitdmin,input%lin%curvefit_dmin,ldiis_coeff,reorder)
                end if
             else
                if (input%lin%constrained_dft) then
                   call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                        infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,update_phi,update_phi,&
                        .true.,ham_small,input%lin%extra_states,itout,it_scc,cdft_it,input%lin%order_taylor,&
                        input%purification_quickreturn,input%adjust_FOE_temperature,&
                        input%calculate_KS_residue,input%calculate_gap,&
                        convcrit_dmin,nitdmin,input%lin%curvefit_dmin,ldiis_coeff,reorder,cdft)
                else
                   call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                        infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,update_phi,update_phi,&
                        .true.,ham_small,input%lin%extra_states,itout,it_scc,cdft_it,input%lin%order_taylor,&
                        input%purification_quickreturn,input%adjust_FOE_temperature,&
                        input%calculate_KS_residue,input%calculate_gap,&
                        convcrit_dmin,nitdmin,input%lin%curvefit_dmin,ldiis_coeff,reorder)
                end if
             end if


             !!! TEMPORARY ##########################################################################
             !!do ii=1,tmb%linmat%denskern%nvctr
             !!     iorb = tmb%linmat%denskern%orb_from_index(1,ii)
             !!     jorb = tmb%linmat%denskern%orb_from_index(2,ii)
             !!     if (iproc==0) write(*,*) 'iorb, jorb, denskern', iorb, jorb, tmb%linmat%denskern%matrix_compr(ii)
             !!  end do
             !!! END TEMPORARY ######################################################################


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
             !if(iproc==0) print *,'energs',energs%ebs,energs%eh,energs%exc,energs%evxc,energs%eexctX,energs%eion,energs%edisp
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
                 call yaml_open_map('Hamiltonian update',flow=.true.)
                 ! Use this subroutine to write the energies, with some
                 ! fake number
                 ! to prevent it from writing too much
                 call write_energies(0,0,energs,0.d0,0.d0,'',.true.)
             end if
             call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
                  tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, &
                  KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, &
                  denspot%rhov, rho_negative)
             if (rho_negative) then
                 call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
                 !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
                 !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
                 !!call clean_rho(iproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
             end if


             ! Mix the density.
             if (input%lin%scf_mode/=LINEAR_MIXPOT_SIMPLE) then
                !if (iproc==0) then
                !    call yaml_map('density mixing; history',mix_hist)
                !end if
                !!write(*,'(a,2es16.9)') 'before mix: sum(denspot%rhov), sum(f_fftgr)', &
                !!                                    sum(denspot%rhov), sum(denspot%mix%f_fftgr(:,:,denspot%mix%i_vrespc(1)))
                !!write(*,'(a,2es16.9)') 'before mix: sum(denspot%rhov), sum(rhopotold)', &
                !!                                    sum(denspot%rhov), sum(rhopotold(1:max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim)))
                !!if (it_scc==1) then
                !!    call mix_main(iproc, nproc, input%lin%scf_mode, mix_hist, input, KSwfn%Lzd%Glr, alpha_mix, &
                !!         denspot, mixdiis, rhopotold, pnrm)
                !!else
                    ! use it_scc+1 since we already have the density from the input guess as iteration 1
                    call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,1.d0-alpha_mix,denspot%mix,&
                         denspot%rhov,it_scc+1,denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
                         at%astruct%cell_dim(1)*at%astruct%cell_dim(2)*at%astruct%cell_dim(3),&
                         pnrm,denspot%dpbox%nscatterarr)
                    call check_negative_rho(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, &
                         denspot%rhov, rho_negative)
                    if (rho_negative) then
                        call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
                    end if
                !!end if
                !!write(*,'(a,2es16.9)') 'after mix: sum(denspot%rhov), sum(f_fftgr)', &
                !!                                   sum(denspot%rhov), sum(denspot%mix%f_fftgr(:,:,denspot%mix%i_vrespc(1)))
                !!write(*,'(a,2es16.9)') 'after mix: sum(denspot%rhov), sum(rhopotold)', &
                !!                                    sum(denspot%rhov), sum(rhopotold(1:max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim)))

                if ((pnrm<convCritMix .or. it_scc==nit_scc) .and. (.not. input%lin%constrained_dft)) then
                   ! calculate difference in density for convergence criterion of outer loop
                   ! ioffset is the buffer which is present for GGA calculations
                   ioffset=KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%i3xcsh
                   pnrm_out=0.d0
                   do i=1,KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p
                      pnrm_out=pnrm_out+(denspot%rhov(ioffset+i)-rhopotOld_out(ioffset+i))**2
                   end do

                   if (nproc > 1) then
                      call mpiallred(pnrm_out, 1, mpi_sum, bigdft_mpi%mpi_comm)
                   end if

                   pnrm_out=sqrt(pnrm_out)/(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*KSwfn%Lzd%Glr%d%n3i*input%nspin)
                   call vcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, &
                     denspot%rhov(1), 1, rhopotOld_out(1), 1)
                end if

             end if

             ! Calculate the new potential.
             !!if(iproc==0) write(*,'(1x,a)') '---------------------------------------------------------------- Updating potential.'
             if (iproc==0) then
!                 if (iproc==0) call yaml_open_map('pot',flow=.true.)
                 !call yaml_map('update potential',.true.)
             end if
             if (iproc==0) call yaml_newline()
             

             call updatePotential(input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
             if (iproc==0) call yaml_close_map()


             ! update occupations wrt eigenvalues (NB for directmin these aren't guaranteed to be true eigenvalues)
             ! switch off for FOE at the moment
             if (input%lin%scf_mode/=LINEAR_FOE) then
                call vcopy(kswfn%orbs%norb,tmb%orbs%eval(1),1,kswfn%orbs%eval(1),1)
                call evaltoocc(iproc,nproc,.false.,input%tel,kswfn%orbs,input%occopt)
                if (bigdft_mpi%iproc ==0) then 
                   call write_eigenvalues_data(0.1d0,kswfn%orbs,mom_vec_fake)
                end if
             end if

             ! Mix the potential
             if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
                !if (iproc==0) then
                !    call yaml_map('potential mixing; history',mix_hist)
                !end if
                !!call mix_main(iproc, nproc, input%lin%scf_mode, mix_hist, input, KSwfn%Lzd%Glr, alpha_mix, &
                !!     denspot, mixdiis, rhopotold, pnrm)
                call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,1.d0-alpha_mix,denspot%mix,&
                     denspot%rhov,it_scc+1,denspot%dpbox%ndims(1),denspot%dpbox%ndims(2),denspot%dpbox%ndims(3),&
                     at%astruct%cell_dim(1)*at%astruct%cell_dim(2)*at%astruct%cell_dim(3),&
                     pnrm,denspot%dpbox%nscatterarr)
                if (pnrm<convCritMix .or. it_scc==nit_scc .and. (.not. input%lin%constrained_dft)) then
                   ! calculate difference in density for convergence criterion of outer loop
                   pnrm_out=0.d0
                   ! for the potential no buffers are present
                   !ioffset=KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%i3xcsh
                   ioffset=0
                   do i=1,KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p
                      pnrm_out=pnrm_out+(denspot%rhov(i+ioffset)-rhopotOld_out(i+ioffset))**2
                   end do

                   if (nproc > 1) then
                      call mpiallred(pnrm_out, 1, mpi_sum, bigdft_mpi%mpi_comm)
                   end if

                   pnrm_out=sqrt(pnrm_out)/(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*KSwfn%Lzd%Glr%d%n3i*input%nspin)
                   call vcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, &
                        denspot%rhov(1), 1, rhopotOld_out(1), 1) 
                end if
             end if

             ! Keep the support functions fixed if they converged and the density
             ! change is below the tolerance already in the very first iteration
             if(it_scc==1 .and. pnrm<convCritMix .and.  info_basis_functions>0) then
                fix_support_functions=.true.
             end if

             if (input%lin%constrained_dft) then
                !call timing(iproc,'constraineddft','ON')
                ! CDFT: see how satisfaction of constraint varies as kernel is updated
                ! CDFT: calculate Tr[Kw]-Nc
                weight_matrix_ = matrices_null()
                call allocate_matrices(tmb%linmat%m, allocate_full=.false., matname='weight_matrix_', mat=weight_matrix_)
                weight_matrix_%matrix_compr=cdft%weight_matrix%matrix_compr
                call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%l,tmb%linmat%m, &
                     tmb%linmat%kernel_,weight_matrix_,&
                     ebs,tmb%coeff,KSwfn%orbs,tmb%orbs,.false.)
                !tmb%linmat%denskern_large%matrix_compr = tmb%linmat%kernel_%matrix_compr
                call deallocate_matrices(weight_matrix_)
                !call timing(iproc,'constraineddft','OF')
             end if

             ! Write some informations.
             call printSummary()

             if (pnrm<convCritMix.and.input%lin%scf_mode/=LINEAR_DIRECT_MINIMIZATION) then
                 info_scf=it_scc
                 if (iproc==0) then
                     !yaml output
                     !call yaml_close_map() !iteration
                     call bigdft_utils_flush(unit=6)
                 end if
                 exit
             else if (pnrm<convCritMix.and.input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
                 if (iproc==0) then
                     !yaml output
                     !call yaml_close_map() !iteration
                     call bigdft_utils_flush(unit=6)
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
                 !call yaml_close_map() !iteration
                 call bigdft_utils_flush(unit=6)
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
             call yaml_open_map(flow=.true.)
             call yaml_comment('iter:'//yaml_toa(it_scc,fmt='(i6)'),hfill='-')
             call printSummary()
             call yaml_close_map() !iteration
             call bigdft_utils_flush(unit=6)
         end if

          ! Close sequence for the optimization steps
          if (iproc==0) then
              call yaml_close_sequence()
          end if

         if (input%lin%constrained_dft) then
            call timing(iproc,'constraineddft','ON')
            !! CDFT: see how satisfaction of constraint varies as kernel is updated
            !! CDFT: calculate Tr[Kw]-Nc
            !call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%denskern,cdft%weight_matrix,&
            !     ebs,tmb%coeff,KSwfn%orbs,tmb%orbs,.false.)

            !! CHECK HERE WHETHER n3d is correct!!
            ! reset rhopotold (to zero) to ensure we don't exit immediately if V only changes a little
            !call to_zero(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, rhopotOld(1)) 
            call vcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, &
                 rhopotOld_out(1), 1, rhopotOld(1), 1) 

            call vcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, &
                 rhopotOld_out(1), 1, denspot%rhov(1), 1)
            call timing(iproc,'constraineddft','OF')

            call updatePotential(input%nspin,denspot,energs%eh,energs%exc,energs%evxc)

            call timing(iproc,'constraineddft','ON')
            ! reset coeffs as well
            call vcopy(tmb%orbs%norb**2,coeff_tmp(1,1),1,tmb%coeff(1,1),1)

            vgrad=ebs-cdft%charge

            ! CDFT: update V (maximizing E wrt V)
            ! CDFT: we updated the kernel in get_coeff so 1st deriv of W wrt V becomes Tr[Kw]-Nc as in CONQUEST
            ! CDFT: 2nd deriv more problematic?
            ! CDFT: use simplest possible scheme for now

            if (iproc==0) write(*,*) ''
            if (iproc==0) write(*,'(a,I4,2x,6(ES12.4e2,2x),2(ES16.6e2,2x))') &
                 'itc, N, Tr(KW), Tr(KW)-N, V*(Tr(KW)-N), V, Vold, EBS, energy',&
                 cdft_it,cdft%charge,ebs,vgrad,cdft%lag_mult*vgrad,cdft%lag_mult,vold,energs%ebs,energy

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
               if (iproc==0) write(*,'(a,I4,2x,6(ES16.6e3,2x))') 'itn, V, Vg',&
                    cdft_it,cdft%lag_mult,vgrad
               cdft%lag_mult=cdft%lag_mult*2.0_gp
            else ! newton
               vgrad2=(vgrad-vgrad_old)/(cdft%lag_mult-vold)
               if (iproc==0) write(*,'(a,I4,2x,6(ES16.6e3,2x))') 'itn, V, Vold, Vg, Vgold, Vg2, Vg/Vg2',&
                    cdft_it,cdft%lag_mult,vold,vgrad,vgrad_old,vgrad2,vgrad/vgrad2
               vold_tmp=cdft%lag_mult
               cdft%lag_mult=vold-vgrad_old/vgrad2
               vold=vold_tmp
            end if

            vgrad_old=vgrad

            call timing(iproc,'constraineddft','OF')

            ! CDFT: exit when W is converged wrt both V and rho
            if (abs(ebs-cdft%charge) < 1.0e-2) exit

         ! if not constrained DFT exit straight away
         else
            exit
         end if
      end do cdft_loop
      if (input%lin%constrained_dft) call DIIS_free(vdiis)
      ! CDFT: end of CDFT loop to find V which correctly imposes constraint and corresponding density

      if(tmb%can_use_transposed) then
          call f_free_ptr(tmb%psit_c)
          call f_free_ptr(tmb%psit_f)
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
               energs, nlpsp, input, &
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

       call get_coeff(iproc,nproc,LINEAR_MIXDENS_SIMPLE,KSwfn%orbs,at,rxyz,denspot,GPU,&
           infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,update_phi,.false.,&
           .true.,ham_small,input%lin%extra_states,itout,0,0,input%lin%order_taylor,&
           input%purification_quickreturn,input%adjust_FOE_temperature,&
           input%calculate_KS_residue,input%calculate_gap)

       !!if (input%lin%scf_mode==LINEAR_FOE) then
       !!    call f_free_ptr(tmb%coeff)
       !!end if
  end if

       if (bigdft_mpi%iproc ==0) then 
          call write_eigenvalues_data(0.1d0,tmb%orbs,mom_vec_fake)
       end if


  !! TEST ##########################
  if (input%lin%new_pulay_correction) then
      !!if (input%lin%scf_mode==LINEAR_FOE .and. ) then
      !!    tmb%coeff=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%coeff')
      !!end if
      call get_coeff(iproc,nproc,LINEAR_MIXDENS_SIMPLE,KSwfn%orbs,at,rxyz,denspot,GPU,&
          infoCoeff,energs,nlpsp,input%SIC,tmb,pnrm,update_phi,.false.,&
          .true.,ham_small,input%lin%extra_states,itout,0,0,input%lin%order_taylor,&
          input%purification_quickreturn,input%adjust_FOE_temperature,&
          input%calculate_KS_residue,input%calculate_gap)
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

  if (input%lin%scf_mode==LINEAR_FOE) then ! deallocate ham_small
     call deallocate_sparse_matrix(ham_small,subname)
  end if

  if (input%lin%constrained_dft) then
     call cdft_data_free(cdft)
     call f_free(coeff_tmp)
  end if


  ! print the final summary
  call print_info(.true.)



  if (iproc==0) call yaml_close_sequence()

  if (input%loewdin_charge_analysis) then
      call loewdin_charge_analysis(iproc, tmb, at, denspot, calculate_overlap_matrix=.true., &
           calculate_ovrlp_half=.true., meth_overlap=0)
      call support_function_multipoles(iproc, tmb, at, denspot)
  end if


  ! Deallocate everything that is not needed any more.
  if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) call DIIS_free(ldiis_coeff)!call deallocateDIIS(ldiis_coeff)
  call deallocateDIIS(ldiis)
  if(input%lin%mixHist_highaccuracy>0) then
      call deallocateMixrhopotDIIS(mixdiis)
  end if
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
     call write_linear_matrices(iproc,nproc,trim(input%dir_output),input%lin%plotBasisFunctions,tmb,at,rxyz)
  end if

  ! not necessarily the best place for it
  !if (input%lin%fragment_calculation) then
  !   !input%lin%plotBasisFunctions
  !   call output_fragment_rotations(iproc,at%astruct%nat,rxyz,1,trim(input%dir_output),input%frag,ref_frags)
  !end if 

  !DEBUG
  !ind=1
  !do iorb=1,tmb%orbs%norbp
  !   write(orbname,*) iorb
  !   ilr=tmb%orbs%inwhichlocreg(iorb+tmb%orbs%isorb)
  !   call plot_wf(trim(adjustl(orbname)),1,at,1.0_dp,tmb%lzd%llr(ilr),KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),&
  !        KSwfn%Lzd%hgrids(3),rxyz,tmb%psi(ind:ind+tmb%Lzd%Llr(ilr)%wfd%nvctr_c+7*tmb%Lzd%Llr(ilr)%wfd%nvctr_f))
  !   ind=ind+tmb%Lzd%Llr(ilr)%wfd%nvctr_c+7*tmb%Lzd%Llr(ilr)%wfd%nvctr_f
  !end do
  ! END DEBUG

  !!! write tmbs in isf format as well
  !!if (input%lin%plotBasisFunctions /= WF_FORMAT_NONE) then
  !!   ! DEBUG - daub_to_isf, write_cube, read_cube, isf_to_daub check the same as starting psi
  !!   ind=1
  !!   !allocate(psi2(tmb%npsidim_orbs),stat=i_stat)
  !!   !call memocc(i_stat,psi2,'psi2',subname)
  !!   do iorb=1,tmb%orbs%norbp
  !!      iat=tmb%orbs%onwhichatom(iorb+tmb%orbs%isorb)
  !!      ilr=tmb%orbs%inwhichlocreg(iorb+tmb%orbs%isorb)
  !!   
  !!      allocate(psir(tmb%lzd%llr(ilr)%d%n1i, tmb%lzd%llr(ilr)%d%n2i, tmb%lzd%llr(ilr)%d%n3i, 1+ndebug),stat=i_stat)
  !!      call memocc(i_stat,psir,'psir',subname)
  !!      !allocate(psir2(tmb%lzd%llr(ilr)%d%n1i, tmb%lzd%llr(ilr)%d%n2i, tmb%lzd%llr(ilr)%d%n3i, 1+ndebug),stat=i_stat)
  !!      !call memocc(i_stat,psir,'psir2',subname)
  !!      call initialize_work_arrays_sumrho(tmb%lzd%llr(ilr),w)
  !!   
  !!      call daub_to_isf(tmb%lzd%llr(ilr),w,tmb%psi(ind),psir)
  !!   
  !!      write(orbname,*) iorb+tmb%orbs%isorb
  !!      !call write_cube_fields('tmbisf'//trim(adjustl(orbname)),'tmb in isf',at,1.0d0,rxyz,&
  !!      !     tmb%lzd%llr(ilr)%d%n1i,tmb%lzd%llr(ilr)%d%n2i,tmb%lzd%llr(ilr)%d%n3i,&
  !!      !     tmb%lzd%llr(ilr)%nsi1,tmb%lzd%llr(ilr)%nsi2,tmb%lzd%llr(ilr)%nsi3,&
  !!      !     tmb%Lzd%hgrids(1)*0.5d0,tmb%Lzd%hgrids(2)*0.5d0,tmb%Lzd%hgrids(3)*0.5d0,&
  !!      !     1.0_gp,psir,1,0.0_gp,psir)

  !!      open(99,file=trim(input%dir_output)//'tmbisf'//trim(adjustl(orbname))//'.dat',&
  !!                form="unformatted",status='unknown')
  !!      write(99) 'Tmb in isf format, to be used in conjunction with minbasis files'
  !!      write(99) tmb%lzd%llr(ilr)%d%n1i,tmb%lzd%llr(ilr)%d%n2i,tmb%lzd%llr(ilr)%d%n3i
  !!      write(99) tmb%lzd%llr(ilr)%nsi1,tmb%lzd%llr(ilr)%nsi2,tmb%lzd%llr(ilr)%nsi3
  !!      do k=1,tmb%lzd%llr(ilr)%d%n3i
  !!         do j=1,tmb%lzd%llr(ilr)%d%n2i
  !!            do i=1,tmb%lzd%llr(ilr)%d%n1i
  !!                 write(99) psir(i,j,k,1)
  !!            end do
  !!         end do
  !!      end do
  !!      close(99)

  !!      !!call read_cube_field('tmbisf'//trim(adjustl(orbname)),tmb%lzd%llr(ilr)%geocode,&
  !!      !!     tmb%lzd%llr(ilr)%d%n1i,tmb%lzd%llr(ilr)%d%n2i,tmb%lzd%llr(ilr)%d%n3i,psir2)

  !!      !open(370,file='tmbisf'//trim(adjustl(orbname))//'.dat')
  !!      !do i=1,tmb%lzd%llr(ilr)%d%n1i
  !!      !do j=1,tmb%lzd%llr(ilr)%d%n2i
  !!      !do k=1,tmb%lzd%llr(ilr)%d%n3i
  !!      !   read(370,*) psir2(i,j,k,1)
  !!      !end do
  !!      !end do
  !!      !end do
  !!      !close(370)

  !!      !call to_zero(tmb%npsidim_orbs,psi2)
  !!      !call isf_to_daub(tmb%lzd%llr(ilr),w,psir2,psi2(ind))
  !!   
  !!      !!tmb_diff=0.0d0
  !!      !!max_tmb_diff=0.0d0
  !!      !!do i=1,tmb%lzd%llr(ilr)%d%n1i
  !!      !!do j=1,tmb%lzd%llr(ilr)%d%n2i
  !!      !!do k=1,tmb%lzd%llr(ilr)%d%n3i
  !!      !!   tmb_diff=tmb_diff+dabs(psir(i,j,k,1)-psir2(i,j,k,1))
  !!      !!   max_tmb_diff=max(max_tmb_diff,dabs(psir(i,j,k,1)-psir2(i,j,k,1)))
  !!      !!!   write(370+iorb+tmb%orbs%isorb,*) psir(i,j,k,1),psir2(i,j,k,1),dabs(psir(i,j,k,1)-psir2(i,j,k,1))
  !!      !!end do
  !!      !!end do
  !!      !!end do
  !!      !!print*,'tmbr diff',iorb+tmb%orbs%isorb,tmb_diff,max_tmb_diff

  !!      !tmb_diff=0.0d0
  !!      !max_tmb_diff=0.0d0
  !!      !n1i=tmb%lzd%llr(ilr)%d%n1i
  !!      !n2i=tmb%lzd%llr(ilr)%d%n2i
  !!      !n3i=tmb%lzd%llr(ilr)%d%n3i
  !!      !do i=0,tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f-1
  !!      !   i3=(i/(n1i*n2i))+1
  !!      !   i2=(i-(i3-1)*n1i*n2i)/n1i+1
  !!      !   i1=mod(i,n1i)+1
  !!      !   tmb_diff=tmb_diff+dabs(tmb%psi(ind+i)-psi2(ind+i))
  !!      !   max_tmb_diff=max(max_tmb_diff,dabs(tmb%psi(ind+i)-psi2(ind+i)))
  !!      !   !write(270+iorb+tmb%orbs%isorb,*) tmb%psi(ind+i),psi2(ind+i),dabs(tmb%psi(ind+i)-psi2(ind+i))
  !!      !   !if (dabs(tmb%psi(ind+i)-psi2(ind+i))>1.0d-5) print*,'large error',iorb+tmb%orbs%isorb,&
  !!      !   !     tmb%psi(ind+i),psi2(ind+i),dabs(tmb%psi(ind+i)-psi2(ind+i)),i1,i2,i3,n1i,n2i,n3i
  !!      !end do
  !!      !print*,'tmb diff',iorb+tmb%orbs%isorb,tmb_diff/(tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f),max_tmb_diff

  !!      ind=ind+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
  !!      call deallocate_work_arrays_sumrho(w)
  !!      i_all=-product(shape(psir))*kind(psir)
  !!      deallocate(psir,stat=i_stat)
  !!      call memocc(i_stat,i_all,'psir',subname)
  !!   end do
  !!end if


  ! check why this is here!
  call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
       tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, &
       denspot%rhov, rho_negative)
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
      locrad_tmp = f_malloc(tmb%lzd%nlr,id='locrad_tmp')

    end subroutine allocate_local_arrays


    subroutine deallocate_local_arrays()

      call f_free(locrad)
      call f_free(locrad_tmp)
      call f_free(rhopotold_out)

    end subroutine deallocate_local_arrays


    subroutine check_inputguess()
      real(8) :: dnrm2
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
                if (associated(tmb%psit_c)) then
                    call f_free_ptr(tmb%psit_c)
                end if
                if (associated(tmb%psit_f)) then
                    call f_free_ptr(tmb%psit_f)
                end if
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
          call yaml_open_sequence('summary',flow=.true.)
          call yaml_open_map()
          if(input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
              call yaml_map('kernel optimization','DMIN')
          else if (input%lin%scf_mode==LINEAR_FOE) then
              call yaml_map('kernel optimization','FOE')
          else
              call yaml_map('kernel optimization','DIAG')
          end if

          if (input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE .or.  input%lin%scf_mode==LINEAR_FOE &
              .or. input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
              call yaml_map('mixing quantity','DENS')
          else if (input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
              call yaml_map('mixing quantity','POT')
          end if
          call yaml_map('mix hist',mix_hist)


          if (input%lin%constrained_dft) then
             if (iproc==0) then
                 call yaml_newline()
                 call yaml_map('iter',it_scc,fmt='(i6)')
                 call yaml_map('delta',pnrm,fmt='(es9.2)')
                 call yaml_map('energy',energy,fmt='(es24.17)')
                 call yaml_map('D',energyDiff,fmt='(es10.3)')
                 call yaml_map('Tr(KW)',ebs,fmt='(es14.4)')
                 call yaml_close_map()
             end if
          else
             if (iproc==0) then
                 call yaml_newline()
                 call yaml_map('iter',it_scc,fmt='(i6)')
                 call yaml_map('delta',pnrm,fmt='(es9.2)')
                 call yaml_map('energy',energy,fmt='(es24.17)')
                 call yaml_map('D',energyDiff,fmt='(es10.3)')
                 call yaml_close_map()
             end if
          end if     
          call yaml_close_sequence()
      end if

    end subroutine printSummary

    !> Print a short summary of some values calculated during the last iteration in the self
    !! consistency cycle.
    subroutine print_info(final)
      implicit none

      real(8) :: energyDiff, mean_conf
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
          call yaml_open_sequence('self consistency summary',label=&
              'it_sc'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)'))))
          call yaml_sequence(advance='no')
          call yaml_open_map(flow=.true.)
          call yaml_map('iter',itout)
          if(target_function==TARGET_FUNCTION_IS_TRACE) then
              call yaml_map('target function','TRACE')
          else if(target_function==TARGET_FUNCTION_IS_ENERGY) then
              call yaml_map('target function','ENERGY')
          else if(target_function==TARGET_FUNCTION_IS_HYBRID) then
              call yaml_map('target function','HYBRID')
          end if
          if (target_function==TARGET_FUNCTION_IS_HYBRID) then
              call yaml_map('mean conf prefac',mean_conf,fmt='(es9.2)')
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
          call yaml_close_map()

          !when convergence is reached, use this block
      else if (iproc==0.and.final) then
          call yaml_comment('final results',hfill='=')
          call yaml_open_sequence('self consistency summary')
          call yaml_sequence(advance='no')
          call yaml_open_map(flow=.true.)
          call yaml_map('iter',itout)
          call write_energies(0,0,energs,0.d0,0.d0,'',.true.)
          if (input%lin%scf_mode/=LINEAR_MIXPOT_SIMPLE) then
             if (.not. lowaccur_converged) then
                 call yaml_map('iter low',itout)
                 call yaml_map('delta out',pnrm_out,fmt='(es10.3)')
                 call yaml_map('energy',energy,fmt='(es27.17)')
                 call yaml_map('D',energyDiff,fmt='(es10.3)')
                 call yaml_comment('FINAL')
                 call yaml_close_map()
             else
                 call yaml_map('iter high',itout)
                 call yaml_map('delta out',pnrm_out,fmt='(es10.3)')
                 call yaml_map('energy',energy,fmt='(es27.17)')
                 call yaml_map('D',energyDiff,fmt='(es10.3)')
                 call yaml_comment('FINAL')
                 call yaml_close_map()
             end if
          else if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
             if (.not. lowaccur_converged) then
                 call yaml_map('iter low',itout)
                 call yaml_map('delta out',pnrm_out,fmt='(es10.3)')
                 call yaml_map('energy',energy,fmt='(es27.17)')
                 call yaml_map('D',energyDiff,fmt='(es10.3)')
                 call yaml_comment('FINAL')
                 call yaml_close_map()
             else
                 call yaml_map('iter high',itout)
                 call yaml_map('delta out',pnrm_out,fmt='(es10.3)')
                 call yaml_map('energy',energy,fmt='(es27.17)')
                 call yaml_map('D',energyDiff,fmt='(es10.3)')
                 call yaml_comment('FINAL')
                 call yaml_close_map()
             end if
          end if
       end if

    call bigdft_utils_flush(unit=6)
    call yaml_close_sequence()


    end subroutine print_info


    subroutine intermediate_forces()

      ! Local variables
      real(kind=8) :: eh_tmp, exc_tmp, evxc_tmp, eexctX_tmp
      real(kind=8),dimension(6) :: ewaldstr, hstrten, xcstr, strten
      real(kind=8),dimension(:),allocatable :: rhopot_work

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
          denspot%pot_work=f_malloc_ptr(denspot%dpbox%ndimpot+ndebug,id='denspot%dpbox%ndimpot+ndebug')
      else
          denspot%pot_work=f_malloc_ptr(1+ndebug,id='denspot%dpbox%ndimpot+ndebug')
      end if
      if (denspot%dpbox%ndimrhopot>0) then
          rhopot_work=f_malloc(denspot%dpbox%ndimrhopot+ndebug,id='rhopot_work')
      else
          rhopot_work=f_malloc(1+ndebug,id='rhopot_work')
      end if
      call vcopy(denspot%dpbox%ndimrhopot,denspot%rhov(1),1,rhopot_work(1),1)


      call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
           tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, &
           denspot%rhov, rho_negative)
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
  integer :: i_stat, i_all, ifrag, jfrag, ifrag_ref, jfrag_ref, iat, isfat, jsfat
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
