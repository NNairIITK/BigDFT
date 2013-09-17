!> @file
!!  Routines used by the linear scaling version
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine linearScaling(iproc,nproc,KSwfn,tmb,at,input,rxyz,denspot,rhopotold,nlpspd,proj,GPU,&
           energs,energy,fpulay,infocode,ref_frags,cdft)
 
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => linearScaling
  use yaml_output
  use module_fragments
  use constrained_dft
  use diis_sd_optimization
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(atoms_data),intent(inout) :: at
  type(input_variables),intent(in) :: input ! need to hack to be inout for geopt changes
  real(8),dimension(3,at%astruct%nat),intent(inout) :: rxyz
  real(8),dimension(3,at%astruct%nat),intent(out) :: fpulay
  type(DFT_local_fields), intent(inout) :: denspot
  real(gp), dimension(*), intent(inout) :: rhopotold
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd
  real(wp),dimension(nlpspd%nprojel),intent(inout) :: proj
  type(GPU_pointers),intent(inout) :: GPU
  type(energy_terms),intent(inout) :: energs
  real(gp), dimension(:), pointer :: rho,pot
  real(8),intent(out) :: energy
  type(DFT_wavefunction),intent(inout),target :: tmb
  type(DFT_wavefunction),intent(inout),target :: KSwfn
  integer,intent(out) :: infocode
  type(system_fragment), dimension(:), pointer :: ref_frags ! for transfer integrals
  type(cdft_data), intent(inout) :: cdft
  
  real(8) :: pnrm,trace,trace_old,fnrm_tmb
  integer :: infoCoeff,istat,iall,it_scc,itout,info_scf,i,ierr,iorb
  character(len=*), parameter :: subname='linearScaling'
  real(8),dimension(:),allocatable :: rhopotold_out
  real(8) :: energyold, energyDiff, energyoldout, fnrm_pulay, convCritMix
  type(mixrhopotDIISParameters) :: mixdiis
  type(localizedDIISParameters) :: ldiis!, ldiis_coeff
  type(DIIS_obj) :: ldiis_coeff, vdiis
  logical :: can_use_ham, update_phi, locreg_increased, reduce_conf, orthonormalization_on
  logical :: fix_support_functions, check_initialguess, fix_supportfunctions
  integer :: itype, istart, nit_lowaccuracy, nit_highaccuracy
  real(8),dimension(:),allocatable :: locrad_tmp
  integer :: ldiis_coeff_hist, nitdmin
  logical :: ldiis_coeff_changed
  integer :: mix_hist, info_basis_functions, nit_scc, cur_it_highaccuracy
  real(8) :: pnrm_out, alpha_mix, ratio_deltas, convcrit_dmin
  logical :: lowaccur_converged, exit_outer_loop
  real(8),dimension(:),allocatable :: locrad
  integer:: target_function, nit_basis
  type(sparseMatrix) :: ham_small
  integer :: isegsmall, iseglarge, iismall, iilarge, is, ie
  integer :: matrixindex_in_compressed
  
  real(kind=gp) :: ebs, vgrad_old, vgrad, valpha, vold, vgrad2, vold_tmp
  real(kind=gp), allocatable, dimension(:,:) :: coeff_tmp
  integer :: ind_denskern, ind_ham, jorb, cdft_it, nelec, iat, ityp, ifrag, ifrag_charged, ifrag_ref, isforb, itmb

  !!! EXPERIMENTAL ############################################
  type(sparseMatrix) :: denskern_init
  real(8),dimension(:),allocatable :: rho_init, rho_init_old, philarge
  real(8) :: tt, ddot, tt_old, meanconf_der
  integer :: idens_cons, ii, sdim, ldim, npsidim_large, ists, istl, nspin, unitname, ilr
  real(8),dimension(10000) :: meanconf_array
  character(len=5) :: num
  character(len=50) :: filename
  !!! #########################################################

  call timing(iproc,'linscalinit','ON') !lr408t

  call f_routine(id='linear_scaling')

  call allocate_local_arrays()

  if(iproc==0) then
      write(*,'(1x,a)') repeat('*',84)
      write(*,'(1x,a)') '****************************** LINEAR SCALING VERSION ******************************'
  end if

  ! Allocate the communications buffers needed for the communications of the potential and
  ! post the messages. This will send to each process the part of the potential that this process
  ! needs for the application of the Hamlitonian to all orbitals on that process.
  call allocateCommunicationsBuffersPotential(tmb%comgp, subname)

  ! Initialize the DIIS mixing of the potential if required.
  if(input%lin%mixHist_lowaccuracy>0) then
      call initializeMixrhopotDIIS(input%lin%mixHist_lowaccuracy, denspot%dpbox%ndimpot, mixdiis)
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

  call nullify_sparsematrix(ham_small) ! nullify anyway

  if (input%lin%scf_mode==LINEAR_FOE) then ! allocate ham_small
     call sparse_copy_pattern(tmb%linmat%ovrlp,ham_small,iproc,subname)
     allocate(ham_small%matrix_compr(ham_small%nvctr), stat=istat)
     call memocc(istat, ham_small%matrix_compr, 'ham_small%matrix_compr', subname)
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

  ! we already have psit in the other case, and in fact the overlap, so eventually could reuse that as well
  if (.not. (input%lin%constrained_dft .and. trim(cdft%method)=='lowdin')) then
     tmb%can_use_transposed=.false.
     nullify(tmb%psit_c)
     nullify(tmb%psit_f)
  end if

  call timing(iproc,'linscalinit','OF') !lr408t

  ! Check the quality of the input guess
  call check_inputguess()

  call timing(iproc,'linscalinit','ON') !lr408t

  if(input%lin%scf_mode/=LINEAR_MIXPOT_SIMPLE) then
      call dcopy(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim),rhopotold(1),1,rhopotold_out(1),1)
  else
      call dcopy(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim),denspot%rhov(1),1,rhopotold_out(1),1)
      call dcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p,1) &
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
     call timing(iproc,'constraineddft','ON')
     if (nit_lowaccuracy>0 .or. input%lin%nItBasis_highaccuracy>1) then
        stop 'Basis cannot be updated for now in constrained DFT calculations and no low accuracy is allowed'
     end if

     call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%denskern,cdft%weight_matrix,&
          ebs,tmb%coeff,KSwfn%orbs,tmb%orbs,.false.)
     vgrad_old=ebs-cdft%charge

     if (iproc==0) print*,'Tr(KW), Tr(KW)-N, V*(Tr(KW)-N)',ebs,vgrad_old,cdft%lag_mult*(ebs-cdft%charge)
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
  !!call dcopy(size(tmb%linmat%denskern%matrix_compr), tmb%linmat%denskern%matrix_compr, 1, denskern_init%matrix_compr, 1)
  !!!call dcopy(size(tmb%linmat%denskern%matrix), tmb%linmat%denskern%matrix, 1, denskern_init%matrix, 1)
  !!allocate(rho_init(size(denspot%rhov)))
  !!allocate(rho_init_old(size(denspot%rhov)))
  !!tt_old=1.d100
  !!rho_init=0.d0
  !!rho_init_old=0.d0
  !!idens_cons=0
  !!! ####################################

  ! modify tmb%orbs%occup, as we normally use orbs%occup elsewhere
  if (input%lin%extra_states>0) then
     call razero(tmb%orbs%norb,tmb%orbs%occup(1))
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
  end if

  call timing(iproc,'linscalinit','OF') !lr408t

  ! Add one iteration if no low accuracy is desired since we need then a first fake iteration, with istart=0
  istart = min(1,nit_lowaccuracy)
  infocode=0 !default value
  ! This is the main outer loop. Each iteration of this loop consists of a first loop in which the basis functions
  ! are optimized and a consecutive loop in which the density is mixed.
  outerLoop: do itout=istart,nit_lowaccuracy+nit_highaccuracy

      if (input%lin%nlevel_accuracy==2) then
          ! Check whether the low accuracy part (i.e. with strong confining potential) has converged.
          call check_whether_lowaccuracy_converged(itout, nit_lowaccuracy, input%lin%lowaccuracy_conv_crit, &
               lowaccur_converged, pnrm_out)
          ! Set all remaining variables that we need for the optimizations of the basis functions and the mixing.
          call set_optimization_variables(input, at, tmb%orbs, tmb%lzd%nlr, tmb%orbs%onwhichatom, tmb%confdatarr, &
               convCritMix, lowaccur_converged, nit_scc, mix_hist, alpha_mix, locrad, target_function, nit_basis, &
               convcrit_dmin, nitdmin)
      else if (input%lin%nlevel_accuracy==1 .and. itout==1) then
          call set_variables_for_hybrid(tmb%lzd%nlr, input, at, tmb%orbs, lowaccur_converged, tmb%confdatarr, &
               target_function, nit_basis, nit_scc, mix_hist, locrad, alpha_mix, convCritMix)
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
          ! Adjust the confining potential if required.
          call adjust_locregs_and_confinement(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
               at, input, rxyz, KSwfn, tmb, denspot, ldiis, locreg_increased, lowaccur_converged, locrad)
          orthonormalization_on=.true.

          if (locreg_increased .and. input%lin%scf_mode==LINEAR_FOE) then ! deallocate ham_small
             call deallocate_sparsematrix(ham_small,subname)
             call nullify_sparsematrix(ham_small)
             call sparse_copy_pattern(tmb%linmat%ovrlp,ham_small,iproc,subname)
             allocate(ham_small%matrix_compr(ham_small%nvctr), stat=istat)
             call memocc(istat, ham_small%matrix_compr, 'ham_small%matrix_compr', subname)
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
                     infoCoeff,energs%ebs,nlpspd,proj,input%SIC,tmb,pnrm,update_phi,update_phi,&
                     .true.,ham_small,input%lin%extra_states,convcrit_dmin,nitdmin,input%lin%curvefit_dmin,ldiis_coeff)
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
              iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
              deallocate(tmb%psit_c, stat=istat)
              call memocc(istat, iall, 'tmb%psit_c', subname)
          end if
          if (associated(tmb%psit_f)) then
              iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
              deallocate(tmb%psit_f, stat=istat)
              call memocc(istat, iall, 'tmb%psit_f', subname)
          end if
          tmb%can_use_transposed=.false.
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
                   if (iproc==0) write(*,'(1x,a,es8.1)') 'Multiply the confinement prefactor by',input%lin%reduce_confinement_factor
                   tmb%confdatarr(:)%prefac=input%lin%reduce_confinement_factor*tmb%confdatarr(:)%prefac
               else
                   if (ratio_deltas<=1.d0 .and. ratio_deltas>0.d0) then
                       if (iproc==0) write(*,'(1x,a,es8.1)') 'Multiply the confinement prefactor by',ratio_deltas
                       tmb%confdatarr(:)%prefac=ratio_deltas*tmb%confdatarr(:)%prefac
                   else if (ratio_deltas>1.d0) then
                       if (iproc==0) write(*,*) 'WARNING: ratio_deltas>1!. Using 0.5 instead'
                       if (iproc==0) write(*,'(1x,a,es8.1)') 'Multiply the confinement prefactor by',0.5d0
                       tmb%confdatarr(:)%prefac=0.5d0*tmb%confdatarr(:)%prefac
                   else if (ratio_deltas<=0.d0) then
                       if (iproc==0) write(*,*) 'WARNING: ratio_deltas<=0.d0!. Using 0.5 instead'
                       if (iproc==0) write(*,'(1x,a,es8.1)') 'Multiply the confinement prefactor by',0.5d0
                       tmb%confdatarr(:)%prefac=0.5d0*tmb%confdatarr(:)%prefac
                   end if
               end if
               !if (iproc==0) write(*,'(a,es18.8)') 'tmb%confdatarr(1)%prefac',tmb%confdatarr(1)%prefac
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
           !if (iproc==0) write(*,*) 'WARNING: set orthonormalization_on to false'
           !orthonormalization_on=.false.
           call getLocalizedBasis(iproc,nproc,at,KSwfn%orbs,rxyz,denspot,GPU,trace,trace_old,fnrm_tmb,&
               info_basis_functions,nlpspd,input%lin%scf_mode,proj,ldiis,input%SIC,tmb,energs, &
               reduce_conf,fix_supportfunctions,input%lin%nItPrecond,target_function,input%lin%correctionOrthoconstraint,&
               nit_basis,input%lin%deltaenergy_multiplier_TMBexit,input%lin%deltaenergy_multiplier_TMBfix,&
               ratio_deltas,orthonormalization_on,input%lin%extra_states,itout)

           ! WRITE SUPPORT FUNCTIONS TO DISK ############################################
           npsidim_large=tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f                                                 
           allocate(philarge((tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f)*tmb%orbs%norbp))                          
           philarge=0.d0
           ists=1                                                                                                          
           istl=1
           do iorb=1,tmb%orbs%norbp
               ilr = tmb%orbs%inWhichLocreg(tmb%orbs%isorb+iorb)                                                           
               sdim=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f                                            
               ldim=tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f                                                      
               nspin=1 !this must be modified later
               call Lpsi_to_global2(iproc, sdim, ldim, tmb%orbs%norb, tmb%orbs%nspinor, nspin, tmb%lzd%glr, &              
                    tmb%lzd%llr(ilr), tmb%psi(ists), philarge(istl))                                                       
               write(num,'(i5.5)') tmb%orbs%isorb+iorb
               filename='supfun_'//num
               unitname=100*iproc+5
               open(unit=unitname,file=trim(filename))
               do i=1,tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
                   write(unitname,'(es25.17)') philarge(istl+i-1)
               end do
               close(unit=unitname)
               ists=ists+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f                                       
               istl=istl+tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f                                                 
           end do
           deallocate(philarge)
           ! ############################################################################

           tmb%can_use_transposed=.false. !since basis functions have changed...

           if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) ldiis_coeff%alpha_coeff=input%lin%alphaSD_coeff !reset to default value

           if (input%inputPsiId==101 .and. info_basis_functions<0 .and. itout==1) then
               ! There seem to be some convergence problems after a restart. Better to quit
               ! and start with a new AO input guess.
               if (iproc==0) write(*,'(1x,a)') 'There are convergence problems after the restart. &
                                                &Start over again with an AO input guess.'
               if (associated(tmb%psit_c)) then
                   iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
                   deallocate(tmb%psit_c, stat=istat)
                   call memocc(istat, iall, 'tmb%psit_c', subname)
               end if
               if (associated(tmb%psit_f)) then
                   iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
                   deallocate(tmb%psit_f, stat=istat)
                   call memocc(istat, iall, 'tmb%psit_f', subname)
               end if
               infocode=2
               exit outerLoop
           end if
       end if

       ! Check whether we can use the Hamiltonian matrix from the TMB optimization
       ! for the first step of the coefficient optimization
       can_use_ham=.true.
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
       end if

      if (can_use_ham .and. input%lin%scf_mode==LINEAR_FOE) then ! copy ham to ham_small here already as it won't be changing
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


      if (input%lin%constrained_dft) then
         call DIIS_set(30,valpha,1,1,vdiis)
         call dcopy(tmb%orbs%norb**2,tmb%coeff(1,1),1,coeff_tmp,1)
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
         kernel_loop : do it_scc=1,nit_scc
             ! If the hamiltonian is available do not recalculate it
             ! also using update_phi for calculate_overlap_matrix and communicate_phi_for_lsumrho
             ! since this is only required if basis changed
             if(update_phi .and. can_use_ham .and. info_basis_functions>=0) then
                if (input%lin%constrained_dft) then
                   call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                        infoCoeff,energs%ebs,nlpspd,proj,input%SIC,tmb,pnrm,update_phi,update_phi,&
                        .false.,ham_small,input%lin%extra_states,convcrit_dmin,nitdmin,&
                        input%lin%curvefit_dmin,ldiis_coeff,cdft)
                else
                   call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                        infoCoeff,energs%ebs,nlpspd,proj,input%SIC,tmb,pnrm,update_phi,update_phi,&
                        .false.,ham_small,input%lin%extra_states,convcrit_dmin,nitdmin,&
                        input%lin%curvefit_dmin,ldiis_coeff)
                end if
             else
                if (input%lin%constrained_dft) then
                   call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                        infoCoeff,energs%ebs,nlpspd,proj,input%SIC,tmb,pnrm,update_phi,update_phi,&
                        .true.,ham_small,input%lin%extra_states,convcrit_dmin,nitdmin,&
                        input%lin%curvefit_dmin,ldiis_coeff,cdft)
                else
                   call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                        infoCoeff,energs%ebs,nlpspd,proj,input%SIC,tmb,pnrm,update_phi,update_phi,&
                        .true.,ham_small,input%lin%extra_states,convcrit_dmin,nitdmin,&
                        input%lin%curvefit_dmin,ldiis_coeff)
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
                if(iproc==0) write(*,*) ''
                if(iproc==0) write(*,*) 'alpha, energydiff',ldiis_coeff%alpha_coeff,energydiff
             end if

             ! Calculate the charge density.
             call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
                  tmb%collcom_sr, tmb%linmat%denskern, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)

             !!! EXPERIMENTAL ################################
             !!if (iproc==0) then
             !!     write(*,*) 'sum1',sum(tmb%linmat%denskern%matrix_compr)
             !!     write(*,*) 'sum2',sum(denskern_init%matrix_compr)
             !!end if
             !!call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
             !!     tmb%collcom_sr, denskern_init, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, rho_init)
             !!     !if (iproc==0) then
             !!     !     write(200,*) denspot%rhov
             !!     !     write(201,*) rho_init
             !!     !end if
             !!!call daxpy(size(denspot%rhov), -1.d0, denspot%rhov(1), 1, rho_init(1), 1)
             !!!tt=ddot(size(denspot%rhov), rho_init(1), 1, rho_init(1), 1)
             !!!call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
             !!!tt=sqrt(tt)/(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*KSwfn%Lzd%Glr%d%n3i*input%nspin)
             !!!if (it_scc==1) then
             !!!    if (iproc==0) write(*,*) 'derivative tt', tt-tt_old
             !!!    tt_old=tt
             !!!end if
             !!!if (iproc==0) write(*,*) 'iproc, tt', iproc, tt
             !!tt=0.d0
             !!do i=1,KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p
             !!    tt=tt+(rho_init(i)-rho_init_old(i))**2
             !!end do
             !!call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
             !!tt=sqrt(tt)/(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*KSwfn%Lzd%Glr%d%n3i*input%nspin)
             !!if (it_scc==1) then
             !!    if (iproc==0) write(*,*) 'charge diff:',tt
             !!    rho_init_old=rho_init
             !!    if (tt<1.d-8) then
             !!        idens_cons=idens_cons+1
             !!    else
             !!        idens_cons=0
             !!    end if
             !!    if (idens_cons>=3) then
             !!        fix_supportfunctions=.true.
             !!        if (iproc==0) write(*,*) 'new convergence criterion: will fix support functions'
             !!    end if
             !!end if

             !!! #############################################

             ! Mix the density.
             if (input%lin%scf_mode/=LINEAR_MIXPOT_SIMPLE) then
                call mix_main(iproc, nproc, mix_hist, input, KSwfn%Lzd%Glr, alpha_mix, &
                     denspot, mixdiis, rhopotold, pnrm)
                if ((pnrm<convCritMix .or. it_scc==nit_scc) .and. (.not. input%lin%constrained_dft)) then
                   ! calculate difference in density for convergence criterion of outer loop
                   pnrm_out=0.d0
                   do i=1,KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p
                      pnrm_out=pnrm_out+(denspot%rhov(i)-rhopotOld_out(i))**2
                   end do
                   call mpiallred(pnrm_out, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
                   pnrm_out=sqrt(pnrm_out)/(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*KSwfn%Lzd%Glr%d%n3i*input%nspin)
                      call dcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, &
                     denspot%rhov(1), 1, rhopotOld_out(1), 1)
                end if
             end if

             ! Calculate the new potential.
             if(iproc==0) write(*,'(1x,a)') '---------------------------------------------------------------- Updating potential.'
             call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)

             ! Mix the potential
             if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
                call mix_main(iproc, nproc, mix_hist, input, KSwfn%Lzd%Glr, alpha_mix, &
                     denspot, mixdiis, rhopotold, pnrm)
                if (pnrm<convCritMix .or. it_scc==nit_scc .and. (.not. input%lin%constrained_dft)) then
                   ! calculate difference in density for convergence criterion of outer loop
                   pnrm_out=0.d0
                   do i=1,KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p
                      pnrm_out=pnrm_out+(denspot%rhov(i)-rhopotOld_out(i))**2
                   end do
                   call mpiallred(pnrm_out, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
                   pnrm_out=sqrt(pnrm_out)/(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*KSwfn%Lzd%Glr%d%n3i*input%nspin)
                   call dcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, &
                        denspot%rhov(1), 1, rhopotOld_out(1), 1) 
                end if
             end if

             ! Keep the support functions fixed if they converged and the density
             ! change is below the tolerance already in the very first iteration
             if(it_scc==1 .and. pnrm<convCritMix .and.  info_basis_functions>0) then
                fix_support_functions=.true.
             end if

             ! Write some informations.
             call printSummary()

             if(pnrm<convCritMix) then
                 info_scf=it_scc
                 exit
             else
                 info_scf=-1
             end if

         end do kernel_loop

         if (input%lin%constrained_dft) then
            call timing(iproc,'constraineddft','ON')
            ! CDFT: see how satisfaction of constraint varies as kernel is updated
            ! CDFT: calculate Tr[Kw]-Nc
            call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%denskern,cdft%weight_matrix,&
                 ebs,tmb%coeff,KSwfn%orbs,tmb%orbs,.false.)
            ! reset rhopotold (to zero) to ensure we don't exit immediately if V only changes a little
            !call razero(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, rhopotOld(1)) 
            call dcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, &
                 rhopotOld_out(1), 1, rhopotOld(1), 1) 

            call dcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, &
                 rhopotOld_out(1), 1, denspot%rhov(1), 1) 
            call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)

            vgrad=ebs-cdft%charge

            ! CDFT: update V (maximizing E wrt V)
            ! CDFT: we updated the kernel in get_coeff so 1st deriv of W wrt V becomes Tr[Kw]-Nc as in CONQUEST
            ! CDFT: 2nd deriv more problematic?
            ! CDFT: use simplest possible scheme for now

            if (iproc==0) write(*,*) ''
            if (iproc==0) write(*,'(a,I4,2x,6(ES16.6e3,2x))') 'itc, Tr(KW), Tr(KW)-N, V*(Tr(KW)-N), V, Vold, EBS',&
                 cdft_it,ebs,vgrad,cdft%lag_mult*vgrad,cdft%lag_mult,vold,energs%ebs

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

            call dcopy(tmb%orbs%norb**2,coeff_tmp(1,1),1,tmb%coeff(1,1),1)
            !if (abs(abs(vgrad)-abs(vgrad_old))>0.1d0) call dcopy(tmb%orbs%norb**2,coeff_tmp(1,1),1,tmb%coeff(1,1),1)
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
          iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
          deallocate(tmb%psit_c, stat=istat)
          call memocc(istat, iall, 'tmb%psit_c', subname)
          iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
          deallocate(tmb%psit_f, stat=istat)
          call memocc(istat, iall, 'tmb%psit_f', subname)
      end if

      call print_info(.false.)

      energyoldout=energy

      call check_for_exit()
      if(exit_outer_loop) exit outerLoop

      if(pnrm_out<input%lin%support_functions_converged.and.lowaccur_converged .or. &
         fix_supportfunctions) then
          if(iproc==0) write(*,*) 'fix the support functions from now on'
          fix_support_functions=.true.
      end if

  end do outerLoop


  ! Diagonalize the matrix for the FOE/direct min case to get the coefficients. Only necessary if
  ! the Pulay forces are to be calculated, or if we are printing eigenvalues for restart
  if ((input%lin%scf_mode==LINEAR_FOE.or.input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION)& 
       .and. (input%lin%pulay_correction.or.input%lin%plotBasisFunctions /= WF_FORMAT_NONE)) then

       call get_coeff(iproc,nproc,LINEAR_MIXDENS_SIMPLE,KSwfn%orbs,at,rxyz,denspot,GPU,&
           infoCoeff,energs%ebs,nlpspd,proj,input%SIC,tmb,pnrm,update_phi,.false.,&
           .true.,ham_small,input%lin%extra_states,convcrit_dmin,nitdmin,input%lin%curvefit_dmin,ldiis_coeff)
  end if

  if (input%lin%scf_mode==LINEAR_FOE) then ! deallocate ham_small
     call deallocate_sparsematrix(ham_small,subname)
  end if

  if (input%lin%constrained_dft) then
     call cdft_data_free(cdft)
     call f_free(coeff_tmp)
  end if


  ! print the final summary
  call print_info(.true.)

  ! Deallocate everything that is not needed any more.
  if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) call DIIS_free(ldiis_coeff)!call deallocateDIIS(ldiis_coeff)
  call deallocateDIIS(ldiis)
  if(input%lin%mixHist_highaccuracy>0) then
      call deallocateMixrhopotDIIS(mixdiis)
  end if
  !!call wait_p2p_communication(iproc, nproc, tmb%comgp)
  call synchronize_onesided_communication(iproc, nproc, tmb%comgp)
  call deallocateCommunicationsBuffersPotential(tmb%comgp, subname)


  if (input%lin%pulay_correction) then
      if (iproc==0) write(*,'(1x,a)') 'WARNING: commented correction_locrad!'
      !!! Testing energy corrections due to locrad
      !!call correction_locrad(iproc, nproc, tmblarge, KSwfn%orbs,tmb%coeff) 
      ! Calculate Pulay correction to the forces
      call pulay_correction(iproc, nproc, KSwfn%orbs, at, rxyz, nlpspd, proj, input%SIC, denspot, GPU, tmb, fpulay)
  else
      call to_zero(3*at%astruct%nat, fpulay(1,1))
  end if

  if(tmb%ham_descr%can_use_transposed) then
      iall=-product(shape(tmb%ham_descr%psit_c))*kind(tmb%ham_descr%psit_c)
      deallocate(tmb%ham_descr%psit_c, stat=istat)
      call memocc(istat, iall, 'tmb%ham_descr%psit_c', subname)
      iall=-product(shape(tmb%ham_descr%psit_f))*kind(tmb%ham_descr%psit_f)
      deallocate(tmb%ham_descr%psit_f, stat=istat)
      call memocc(istat, iall, 'tmb%ham_descr%psit_f', subname)
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
  !   call output_fragment_rotations(iproc,nproc,at%astruct%nat,rxyz,1,trim(input%dir_output),input%frag,ref_frags)
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

  ! check why this is here!
  call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
       tmb%collcom_sr, tmb%linmat%denskern, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)


  ! Otherwise there are some problems... Check later.
  allocate(KSwfn%psi(1),stat=istat)
  call memocc(istat,KSwfn%psi,'KSwfn%psi',subname)
  nullify(KSwfn%psit)

  nullify(rho,pot)

  call deallocate_local_arrays()
  call f_release_routine()

  call timing(bigdft_mpi%mpi_comm,'WFN_OPT','PR')

  contains

    subroutine allocate_local_arrays()

      allocate(locrad(tmb%lzd%nlr), stat=istat)
      call memocc(istat, locrad, 'locrad', subname)

      ! Allocate the old charge density (used to calculate the variation in the charge density)
      allocate(rhopotold_out(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim)),stat=istat)
      call memocc(istat, rhopotold_out, 'rhopotold_out', subname)

      allocate(locrad_tmp(tmb%lzd%nlr), stat=istat)
      call memocc(istat, locrad_tmp, 'locrad_tmp', subname)

    end subroutine allocate_local_arrays


    subroutine deallocate_local_arrays()

      iall=-product(shape(locrad))*kind(locrad)
      deallocate(locrad, stat=istat)
      call memocc(istat, iall, 'locrad', subname)

      iall=-product(shape(locrad_tmp))*kind(locrad_tmp)
      deallocate(locrad_tmp, stat=istat)
      call memocc(istat, iall, 'locrad_tmp', subname)

      iall=-product(shape(rhopotold_out))*kind(rhopotold_out)
      deallocate(rhopotold_out, stat=istat)
      call memocc(istat, iall, 'rhopotold_out', subname)

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
             call pulay_correction(iproc, nproc, KSwfn%orbs, at, rxyz, nlpspd, proj, input%SIC, denspot, GPU, tmb, fpulay)
             fnrm_pulay=dnrm2(3*at%astruct%nat, fpulay, 1)/sqrt(dble(at%astruct%nat))

             if (iproc==0) write(*,*) 'fnrm_pulay',fnrm_pulay

             if (fnrm_pulay>1.d-1) then !1.d-10
                if (iproc==0) write(*,'(1x,a)') 'The pulay force is too large after the restart. &
                                                   &Start over again with an AO input guess.'
                if (associated(tmb%psit_c)) then
                    iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
                    deallocate(tmb%psit_c, stat=istat)
                    call memocc(istat, iall, 'tmb%psit_c', subname)
                end if
                if (associated(tmb%psit_f)) then
                    iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
                    deallocate(tmb%psit_f, stat=istat)
                    call memocc(istat, iall, 'tmb%psit_f', subname)
                end if
                tmb%can_use_transposed=.false.
                nit_lowaccuracy=input%lin%nit_lowaccuracy
                nit_highaccuracy=input%lin%nit_highaccuracy
                !!input%lin%highaccuracy_conv_crit=1.d-8
                call inputguessConfinement(iproc, nproc, at, input, &
                     KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
                     rxyz, nlpspd, proj, GPU, KSwfn%orbs, tmb, denspot, rhopotold, energs)
                     energs%eexctX=0.0_gp

                !already done in inputguess
                      ! CHEATING here and passing tmb%linmat%denskern instead of tmb%linmat%inv_ovrlp
                !call orthonormalizeLocalized(iproc, nproc, 0, tmb%npsidim_orbs, tmb%orbs, tmb%lzd, tmb%linmat%ovrlp, &
                !     tmb%linmat%denskern, tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)
             else if (fnrm_pulay>1.d-2) then ! 1.d2 1.d-2
                if (iproc==0) write(*,'(1x,a)') 'The pulay forces are rather large, so start with low accuracy.'
                nit_lowaccuracy=input%lin%nit_lowaccuracy
                nit_highaccuracy=input%lin%nit_highaccuracy
             else if (fnrm_pulay>1.d-10) then !1d-10
                if (iproc==0) write(*,'(1x,a)') &
                     'The pulay forces are fairly large, so reoptimising basis with high accuracy only.'
                nit_lowaccuracy=0
                nit_highaccuracy=input%lin%nit_highaccuracy
             else
                if (iproc==0) write(*,'(1x,a)') &
                     'The pulay forces are fairly small, so not reoptimising basis.'
                    nit_lowaccuracy=0
                    nit_highaccuracy=0
                    !!nit_highaccuracy=2
                    !!input%lin%nItBasis_highaccuracy=2
                    !!input%lin%nitSCCWhenFixed_highaccuracy=100
                    !!input%lin%highaccuracy_conv_crit=1.d-8
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
          write(*,'(1x,a)') repeat('+',92 + int(log(real(it_SCC))/log(10.)))
          write(*,'(1x,a,i0,a)') 'at iteration ', it_SCC, ' of the density optimization:'
          !!if(infoCoeff<0) then
          !!    write(*,'(3x,a)') '- WARNING: coefficients not converged!'
          !!else if(infoCoeff>0) then
          !!    write(*,'(3x,a,i0,a)') '- coefficients converged in ', infoCoeff, ' iterations.'
          if(input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
              write(*,'(3x,a)') 'coefficients / kernel obtained by direct minimization.'
          else if (input%lin%scf_mode==LINEAR_FOE) then
              write(*,'(3x,a)') 'kernel obtained by Fermi Operator Expansion'
          else
              write(*,'(3x,a)') 'coefficients / kernel obtained by diagonalization.'
          end if
          if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE .or. input%lin%scf_mode==LINEAR_FOE) then
              write(*,'(3x,a,3x,i0,2x,es13.7,es27.17,es14.4)') 'it, Delta DENS, energy, energyDiff', &
                   it_SCC, pnrm, energy, energyDiff
          else if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
              write(*,'(3x,a,3x,i0,2x,es13.7,es27.17,es14.4)') 'it, Delta POT, energy, energyDiff', &
                   it_SCC, pnrm, energy, energyDiff
          else if(input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
              write(*,'(3x,a,3x,i0,2x,es13.7,es27.17,es14.4)') 'it, Delta DENS, energy, energyDiff', &
                   it_SCC, pnrm, energy, energyDiff
          end if
          write(*,'(1x,a)') repeat('+',92 + int(log(real(it_SCC))/log(10.)))
      end if

    end subroutine printSummary

    !> Print a short summary of some values calculated during the last iteration in the self
    !! consistency cycle.
    subroutine print_info(final)
      implicit none

      real(8) :: energyDiff, mean_conf
      logical, intent(in) :: final

      energyDiff = energy - energyoldout

      if(target_function==TARGET_FUNCTION_IS_HYBRID) then
          mean_conf=0.d0
          do iorb=1,tmb%orbs%norbp
              mean_conf=mean_conf+tmb%confdatarr(iorb)%prefac
          end do
          call mpiallred(mean_conf, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
          mean_conf=mean_conf/dble(tmb%orbs%norb)
      end if

      ! Print out values related to two iterations of the outer loop.
      if(iproc==0.and.(.not.final)) then

          !Before convergence
          write(*,'(3x,a,7es20.12)') 'ebs, ehart, eexcu, vexcu, eexctX, eion, edisp', &
              energs%ebs, energs%eh, energs%exc, energs%evxc, energs%eexctX, energs%eion, energs%edisp
          if(input%lin%scf_mode/=LINEAR_MIXPOT_SIMPLE) then
             if (.not. lowaccur_converged) then
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                      'itoutL, Delta DENSOUT, energy, energyDiff', itout, pnrm_out, energy, &
                      energyDiff
             else
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                      'itoutH, Delta DENSOUT, energy, energyDiff', itout, pnrm_out, energy, &
                      energyDiff
             end if
          else if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
             if (.not. lowaccur_converged) then
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                      'itoutL, Delta POTOUT, energy energyDiff', itout, pnrm_out, energy, energyDiff
             else
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                      'itoutH, Delta POTOUT, energy energyDiff', itout, pnrm_out, energy, energyDiff
             end if
          end if

          !when convergence is reached, use this block
          write(*,'(1x,a)') repeat('#',92 + int(log(real(itout))/log(10.)))
          write(*,'(1x,a,i0,a)') 'at iteration ', itout, ' of the outer loop:'
          write(*,'(3x,a)') '> basis functions optimization:'
          if(target_function==TARGET_FUNCTION_IS_TRACE) then
              write(*,'(5x,a)') '- target function is trace'
          else if(target_function==TARGET_FUNCTION_IS_ENERGY) then
              write(*,'(5x,a)') '- target function is energy'
          else if(target_function==TARGET_FUNCTION_IS_HYBRID) then
              write(*,'(5x,a,es8.2)') '- target function is hybrid; mean confinement prefactor = ',mean_conf
          end if
          if(info_basis_functions<=0) then
              write(*,'(5x,a)') '- WARNING: basis functions not converged!'
          else
              write(*,'(5x,a,i0,a)') '- basis functions converged in ', info_basis_functions, ' iterations.'
          end if
          write(*,'(5x,a,es15.6,2x,es10.2)') 'Final values: target function, fnrm', trace, fnrm_tmb
          write(*,'(3x,a)') '> density optimization:'
          if(input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
              write(*,'(5x,a)') '- using direct minimization.'
          else if (input%lin%scf_mode==LINEAR_FOE) then
              write(*,'(5x,a)') '- using Fermi Operator Expansion / mixing.'
          else
              write(*,'(5x,a)') '- using diagonalization / mixing.'
          end if
          if(info_scf<0) then
              write(*,'(5x,a)') '- WARNING: density optimization not converged!'
          else
              write(*,'(5x,a,i0,a)') '- density optimization converged in ', info_scf, ' iterations.'
          end if
          if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE .or. input%lin%scf_mode==LINEAR_FOE) then
              write(*,'(5x,a,3x,i0,es12.2,es27.17)') 'FINAL values: it, Delta DENS, energy', itout, pnrm, energy
          else if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
              write(*,'(5x,a,3x,i0,es12.2,es27.17)') 'FINAL values: it, Delta POT, energy', itout, pnrm, energy
          else if(input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
              write(*,'(5x,a,3x,i0,es12.2,es27.17)') 'FINAL values: it, Delta DENS, energy', itout, pnrm, energy
          end if
          write(*,'(3x,a,es14.6)') '> energy difference to last iteration:', energyDiff
          write(*,'(1x,a)') repeat('#',92 + int(log(real(itout))/log(10.)))
      else if (iproc==0.and.final) then
          !Before convergence
          write(*,'(3x,a,7es20.12)') 'ebs, ehart, eexcu, vexcu, eexctX, eion, edisp', &
              energs%ebs, energs%eh, energs%exc, energs%evxc, energs%eexctX, energs%eion, energs%edisp
          if(input%lin%scf_mode/=LINEAR_MIXPOT_SIMPLE) then
             if (.not. lowaccur_converged) then
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4,3x,a)')&
                      'itoutL, Delta DENSOUT, energy, energyDiff', itout, pnrm_out, energy, &
                      energyDiff,'FINAL'
             else
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4,3x,a)')&
                      'itoutH, Delta DENSOUT, energy, energyDiff', itout, pnrm_out, energy, &
                      energyDiff,'FINAL'
             end if
          else if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
             if (.not. lowaccur_converged) then
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4,3x,a)')&
                      'itoutL, Delta POTOUT, energy energyDiff', itout, pnrm_out, energy, energyDiff,'FINAL'
             else
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4,3x,a)')&
                      'itoutH, Delta POTOUT, energy energyDiff', itout, pnrm_out, energy, energyDiff,'FINAL'
             end if
          end if
       end if

    ! WARNING HACK S.M.
    meanconf_array(itout)=mean_conf
    if (itout>=4) then
        meanconf_der = 11.d0/6.d0*meanconf_array(itout-0) &
                      -      3.d0*meanconf_array(itout-1) &
                      + 3.d0/2.d0*meanconf_array(itout-2) &
                      - 1.d0/3.d0*meanconf_array(itout-3)
        if (iproc==0) write(*,'(a,es16.5)') 'meanconf_der',meanconf_der
        if (iproc==0) write(*,'(a,es16.5)') 'abs(meanconf_der)/mean_conf',abs(meanconf_der)/mean_conf
    end if
    !if (mean_conf<1.d-15 .and. .false.) then
    !if (mean_conf<1.d-15) then
    !if (mean_conf<1.d-10 .and. abs(meanconf_der)<1.d-15) then
    if (mean_conf<1.d-10 .and. abs(meanconf_der)/mean_conf>1.d0 .and. .false.) then
    !if (itout>=38) then
        !if (iproc==0) write(*,*) 'WARNING MODIFY CONF'
        !tmb%confdatarr(:)%prefac=0.d0
        if (iproc==0) write(*,*) 'WARNING MODIFY nit_basis'
        nit_basis=0
    end if
    if (mean_conf<1.d-5) then
    !if (itout>=13) then
        if (iproc==0) write(*,*) 'outswitch off ortho'
        orthonormalization_on=.false.
    end if


    end subroutine print_info


end subroutine linearScaling



subroutine set_optimization_variables(input, at, lorbs, nlr, onwhichatom, confdatarr, &
     convCritMix, lowaccur_converged, nit_scc, mix_hist, alpha_mix, locrad, target_function, nit_basis, &
     convcrit_dmin, nitdmin)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: nlr
  type(orbitals_data),intent(in) :: lorbs
  type(input_variables),intent(in) :: input
  type(atoms_data),intent(in) :: at
  integer,dimension(lorbs%norb),intent(in) :: onwhichatom
  type(confpot_data),dimension(lorbs%norbp),intent(inout) :: confdatarr
  real(kind=8), intent(out) :: convCritMix, alpha_mix, convcrit_dmin
  logical, intent(in) :: lowaccur_converged
  integer, intent(out) :: nit_scc, mix_hist, nitdmin
  real(kind=8), dimension(nlr), intent(out) :: locrad
  integer, intent(out) :: target_function, nit_basis

  ! Local variables
  integer :: iorb, ilr, iiat

  if(lowaccur_converged) then
      do iorb=1,lorbs%norbp
          iiat=onwhichatom(lorbs%isorb+iorb)
          confdatarr(iorb)%prefac=input%lin%potentialPrefac_highaccuracy(at%astruct%iatype(iiat))
      end do
      target_function=TARGET_FUNCTION_IS_ENERGY
      nit_basis=input%lin%nItBasis_highaccuracy
      nit_scc=input%lin%nitSCCWhenFixed_highaccuracy
      mix_hist=input%lin%mixHist_highaccuracy
      do ilr=1,nlr
          locrad(ilr)=input%lin%locrad_highaccuracy(ilr)
      end do
      alpha_mix=input%lin%alpha_mix_highaccuracy
      convCritMix=input%lin%convCritMix_highaccuracy
      convcrit_dmin=input%lin%convCritDmin_highaccuracy
      nitdmin=input%lin%nItdmin_highaccuracy
  else
      do iorb=1,lorbs%norbp
          iiat=onwhichatom(lorbs%isorb+iorb)
          confdatarr(iorb)%prefac=input%lin%potentialPrefac_lowaccuracy(at%astruct%iatype(iiat))
      end do
      target_function=TARGET_FUNCTION_IS_TRACE
      nit_basis=input%lin%nItBasis_lowaccuracy
      nit_scc=input%lin%nitSCCWhenFixed_lowaccuracy
      mix_hist=input%lin%mixHist_lowaccuracy
      do ilr=1,nlr
          locrad(ilr)=input%lin%locrad_lowaccuracy(ilr)
      end do
      alpha_mix=input%lin%alpha_mix_lowaccuracy
      convCritMix=input%lin%convCritMix_lowaccuracy
      convcrit_dmin=input%lin%convCritDmin_lowaccuracy
      nitdmin=input%lin%nItdmin_lowaccuracy
  end if

  !!! new hybrid version... not the best place here
  !!if (input%lin%nit_highaccuracy==-1) then
  !!    do iorb=1,lorbs%norbp
  !!        ilr=lorbs%inwhichlocreg(lorbs%isorb+iorb)
  !!        iiat=onwhichatom(lorbs%isorb+iorb)
  !!        confdatarr(iorb)%prefac=input%lin%potentialPrefac_lowaccuracy(at%astruct%iatype(iiat))
  !!    end do
  !!    wfnmd%bs%target_function=TARGET_FUNCTION_IS_HYBRID
  !!    wfnmd%bs%nit_basis_optimization=input%lin%nItBasis_lowaccuracy
  !!    wfnmd%bs%conv_crit=input%lin%convCrit_lowaccuracy
  !!    nit_scc=input%lin%nitSCCWhenFixed_lowaccuracy
  !!    mix_hist=input%lin%mixHist_lowaccuracy
  !!    do ilr=1,nlr
  !!        locrad(ilr)=input%lin%locrad_lowaccuracy(ilr)
  !!    end do
  !!    alpha_mix=input%lin%alpha_mix_lowaccuracy
  !!end if

end subroutine set_optimization_variables



subroutine adjust_locregs_and_confinement(iproc, nproc, hx, hy, hz, at, input, &
           rxyz, KSwfn, tmb, denspot, ldiis, locreg_increased, lowaccur_converged, locrad)
  use module_base
  use module_types
  use module_interfaces, except_this_one => adjust_locregs_and_confinement
  implicit none
  
  ! Calling argument
  integer,intent(in) :: iproc, nproc
  real(8),intent(in) :: hx, hy, hz
  type(atoms_data),intent(in) :: at
  type(input_variables),intent(in) :: input
  real(8),dimension(3,at%astruct%nat),intent(in):: rxyz
  type(DFT_wavefunction),intent(inout) :: KSwfn, tmb
  type(DFT_local_fields),intent(inout) :: denspot
  type(localizedDIISParameters),intent(inout) :: ldiis
  logical, intent(out) :: locreg_increased
  logical, intent(in) :: lowaccur_converged
  real(8), dimension(tmb%lzd%nlr), intent(inout) :: locrad

  ! Local variables
  integer :: iall, istat, ilr, npsidim_orbs_tmp, npsidim_comp_tmp
  real(kind=8),dimension(:,:),allocatable :: locregCenter
  real(kind=8),dimension(:),allocatable :: lphilarge
  type(local_zone_descriptors) :: lzd_tmp
  character(len=*), parameter :: subname='adjust_locregs_and_confinement'

  locreg_increased=.false.
  if(lowaccur_converged ) then
      do ilr = 1, tmb%lzd%nlr
         if(input%lin%locrad_highaccuracy(ilr) /= input%lin%locrad_lowaccuracy(ilr)) then
             if(iproc==0) write(*,'(1x,a)') 'Increasing the localization radius for the high accuracy part.'
             locreg_increased=.true.
             exit
         end if
      end do
  end if

  if(locreg_increased) then
     !tag=1
     !call wait_p2p_communication(iproc, nproc, tmb%comgp)
     call synchronize_onesided_communication(iproc, nproc, tmb%comgp)
     call deallocate_p2pComms(tmb%comgp, subname)

     call deallocate_collective_comms(tmb%collcom, subname)
     call deallocate_collective_comms(tmb%collcom_sr, subname)

     call nullify_local_zone_descriptors(lzd_tmp)
     call copy_local_zone_descriptors(tmb%lzd, lzd_tmp, subname)
     call deallocate_local_zone_descriptors(tmb%lzd, subname)

     npsidim_orbs_tmp = tmb%npsidim_orbs
     npsidim_comp_tmp = tmb%npsidim_comp

     call deallocate_foe(tmb%foe_obj, subname)

     call deallocate_sparseMatrix(tmb%linmat%denskern, subname)
     call deallocate_sparseMatrix(tmb%linmat%inv_ovrlp, subname)
     call deallocate_sparseMatrix(tmb%linmat%ovrlp, subname)
     call deallocate_sparseMatrix(tmb%linmat%ham, subname)

     allocate(locregCenter(3,lzd_tmp%nlr), stat=istat)
     call memocc(istat, locregCenter, 'locregCenter', subname)
     do ilr=1,lzd_tmp%nlr
        locregCenter(:,ilr)=lzd_tmp%llr(ilr)%locregCenter
     end do

     !temporary,  moved from update_locreg
     tmb%orbs%eval=-0.5_gp
     call update_locreg(iproc, nproc, lzd_tmp%nlr, locrad, locregCenter, lzd_tmp%glr, .false., &
          denspot%dpbox%nscatterarr, hx, hy, hz, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%lzd, &
          tmb%npsidim_orbs, tmb%npsidim_comp, tmb%comgp, tmb%collcom, tmb%foe_obj, tmb%collcom_sr)

     iall=-product(shape(locregCenter))*kind(locregCenter)
     deallocate(locregCenter, stat=istat)
     call memocc(istat, iall, 'locregCenter', subname)

     ! calculate psi in new locreg
     allocate(lphilarge(tmb%npsidim_orbs), stat=istat)
     call memocc(istat, lphilarge, 'lphilarge', subname)
     call to_zero(tmb%npsidim_orbs, lphilarge(1))
     call small_to_large_locreg(iproc, npsidim_orbs_tmp, tmb%npsidim_orbs, lzd_tmp, tmb%lzd, &
          tmb%orbs, tmb%psi, lphilarge)

     call deallocate_local_zone_descriptors(lzd_tmp, subname)
     iall=-product(shape(tmb%psi))*kind(tmb%psi)
     deallocate(tmb%psi, stat=istat)
     call memocc(istat, iall, 'tmb%psi', subname)
     allocate(tmb%psi(tmb%npsidim_orbs), stat=istat)
     call memocc(istat, tmb%psi, 'tmb%psi', subname)
     call dcopy(tmb%npsidim_orbs, lphilarge(1), 1, tmb%psi(1), 1)
     iall=-product(shape(lphilarge))*kind(lphilarge)
     deallocate(lphilarge, stat=istat)
     call memocc(istat, iall, 'lphilarge', subname) 
     
     call update_ldiis_arrays(tmb, subname, ldiis)

     ! Emit that lzd has been changed.
     if (tmb%c_obj /= 0) then
        call kswfn_emit_lzd(tmb, iproc, nproc)
     end if

     ! Now update hamiltonian descriptors
     !call destroy_new_locregs(iproc, nproc, tmblarge)

     ! to eventually be better sorted - replace with e.g. destroy_hamiltonian_descriptors
     call synchronize_onesided_communication(iproc, nproc, tmb%ham_descr%comgp)
     call deallocate_p2pComms(tmb%ham_descr%comgp, subname)
     call deallocate_local_zone_descriptors(tmb%ham_descr%lzd, subname)
     call deallocate_collective_comms(tmb%ham_descr%collcom, subname)

     call deallocate_auxiliary_basis_function(subname, tmb%ham_descr%psi, tmb%hpsi)
     if(tmb%ham_descr%can_use_transposed) then
        iall=-product(shape(tmb%ham_descr%psit_c))*kind(tmb%ham_descr%psit_c)
        deallocate(tmb%ham_descr%psit_c, stat=istat)
        call memocc(istat, iall, 'tmb%ham_descr%psit_c', subname)
        iall=-product(shape(tmb%ham_descr%psit_f))*kind(tmb%ham_descr%psit_f)
        deallocate(tmb%ham_descr%psit_f, stat=istat)
        call memocc(istat, iall, 'tmb%ham_descr%psit_f', subname)
        tmb%ham_descr%can_use_transposed=.false.
     end if
     
     deallocate(tmb%confdatarr, stat=istat)

     call create_large_tmbs(iproc, nproc, KSwfn, tmb, denspot, input, at, rxyz, lowaccur_converged)

     ! Update sparse matrices
     call initSparseMatrix(iproc, nproc, tmb%ham_descr%lzd, tmb%orbs, input, tmb%linmat%ham)
     call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
          tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%ham)
     call initSparseMatrix(iproc, nproc, tmb%lzd, tmb%orbs, input, tmb%linmat%ovrlp)
     call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
          tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%ovrlp)
     !call initSparseMatrix(iproc, nproc, tmb%ham_descr%lzd, tmb%orbs, tmb%linmat%inv_ovrlp)
     call initSparseMatrix(iproc, nproc, tmb%ham_descr%lzd, tmb%orbs, input, tmb%linmat%denskern)
     call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
          tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%denskern)
     call nullify_sparsematrix(tmb%linmat%inv_ovrlp)
     call sparse_copy_pattern(tmb%linmat%denskern,tmb%linmat%inv_ovrlp,iproc,subname) ! save recalculating
     !call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%inv_ovrlp)

     allocate(tmb%linmat%denskern%matrix_compr(tmb%linmat%denskern%nvctr), stat=istat)
     call memocc(istat, tmb%linmat%denskern%matrix_compr, 'tmb%linmat%denskern%matrix_compr', subname)
     allocate(tmb%linmat%ham%matrix_compr(tmb%linmat%ham%nvctr), stat=istat)
     call memocc(istat, tmb%linmat%ham%matrix_compr, 'tmb%linmat%ham%matrix_compr', subname)
     allocate(tmb%linmat%ovrlp%matrix_compr(tmb%linmat%ovrlp%nvctr), stat=istat)
     call memocc(istat, tmb%linmat%ovrlp%matrix_compr, 'tmb%linmat%ovrlp%matrix_compr', subname)
     !allocate(tmb%linmat%inv_ovrlp%matrix_compr(tmb%linmat%inv_ovrlp%nvctr), stat=istat)
     !call memocc(istat, tmb%linmat%inv_ovrlp%matrix_compr, 'tmb%linmat%inv_ovrlp%matrix_compr', subname)

  else ! no change in locrad, just confining potential that needs updating

     call define_confinement_data(tmb%confdatarr,tmb%orbs,rxyz,at,&
          tmb%ham_descr%lzd%hgrids(1),tmb%ham_descr%lzd%hgrids(2),tmb%ham_descr%lzd%hgrids(3),&
          4,input%lin%potentialPrefac_highaccuracy,tmb%ham_descr%lzd,tmb%orbs%onwhichatom)

  end if

end subroutine adjust_locregs_and_confinement



subroutine adjust_DIIS_for_high_accuracy(input, denspot, mixdiis, lowaccur_converged, ldiis_coeff_hist, ldiis_coeff_changed)
  use module_base
  use module_types
  use module_interfaces, except_this_one => adjust_DIIS_for_high_accuracy
  implicit none
  
  ! Calling arguments
  type(input_variables),intent(in) :: input
  type(DFT_local_fields),intent(inout) :: denspot
  type(mixrhopotDIISParameters),intent(inout) :: mixdiis
  logical, intent(in) :: lowaccur_converged
  integer, intent(inout) :: ldiis_coeff_hist
  logical, intent(out) :: ldiis_coeff_changed  

  if(lowaccur_converged) then
     if(input%lin%mixHist_lowaccuracy==0 .and. input%lin%mixHist_highaccuracy>0) then
        call initializeMixrhopotDIIS(input%lin%mixHist_highaccuracy, denspot%dpbox%ndimpot, mixdiis)
     else if(input%lin%mixHist_lowaccuracy>0 .and. input%lin%mixHist_highaccuracy==0) then
        call deallocateMixrhopotDIIS(mixdiis)
     end if
     if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
        ! check whether ldiis_coeff_hist arrays will need reallocating due to change in history length
        if (ldiis_coeff_hist /= input%lin%dmin_hist_highaccuracy) then
           ldiis_coeff_changed=.true.
        else
           ldiis_coeff_changed=.false.
        end if
        ldiis_coeff_hist=input%lin%dmin_hist_highaccuracy
     end if
  else
     if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
        ldiis_coeff_changed=.false.
     end if
  end if
  
end subroutine adjust_DIIS_for_high_accuracy


subroutine check_whether_lowaccuracy_converged(itout, nit_lowaccuracy, lowaccuracy_convcrit, &
     lowaccur_converged, pnrm_out)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: itout
  integer,intent(in) :: nit_lowaccuracy
  real(8),intent(in) :: lowaccuracy_convcrit
  logical, intent(inout) :: lowaccur_converged
  real(kind=8), intent(in) :: pnrm_out
  
  if(.not.lowaccur_converged .and. &
       (itout>=nit_lowaccuracy+1 .or. pnrm_out<lowaccuracy_convcrit)) then
     lowaccur_converged=.true.
     !cur_it_highaccuracy=0
  end if 

end subroutine check_whether_lowaccuracy_converged

subroutine pulay_correction(iproc, nproc, orbs, at, rxyz, nlpspd, proj, SIC, denspot, GPU, tmb, fpulay)
  use module_base
  use module_types
  use module_interfaces, except_this_one => pulay_correction
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(atoms_data),intent(in) :: at
  real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd
  real(wp),dimension(nlpspd%nprojel),intent(inout) :: proj
  type(SIC_data),intent(in) :: SIC
  type(DFT_local_fields), intent(inout) :: denspot
  type(GPU_pointers),intent(inout) :: GPU
  type(DFT_wavefunction),intent(inout) :: tmb
  real(kind=8),dimension(3,at%astruct%nat),intent(out) :: fpulay

  ! Local variables
  integer:: istat, iall, ierr, iialpha, jorb
  integer:: iorb, ii, iseg, isegstart, isegend
  integer:: jat, jdir, ibeta
  !!integer :: ialpha, iat, iiorb
  real(kind=8) :: kernel, ekernel
  real(kind=8),dimension(:),allocatable :: lhphilarge, psit_c, psit_f, hpsit_c, hpsit_f, lpsit_c, lpsit_f
  type(sparseMatrix) :: dovrlp(3), dham(3)
  type(energy_terms) :: energs
  type(confpot_data),dimension(:),allocatable :: confdatarrtmp
  character(len=*),parameter :: subname='pulay_correction'

  ! Begin by updating the Hpsi
  call local_potential_dimensions(tmb%ham_descr%lzd,tmb%orbs,denspot%dpbox%ngatherarr(0,1))

  allocate(lhphilarge(tmb%ham_descr%npsidim_orbs), stat=istat)
  call memocc(istat, lhphilarge, 'lhphilarge', subname)
  call to_zero(tmb%ham_descr%npsidim_orbs,lhphilarge(1))

  !!call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
  !!     tmb%ham_descr%comgp%nrecvbuf, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, tmb%ham_descr%lzd)
  call start_onesided_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
       tmb%ham_descr%comgp%nrecvbuf, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, tmb%ham_descr%lzd)

  allocate(confdatarrtmp(tmb%orbs%norbp))
  call default_confinement_data(confdatarrtmp,tmb%orbs%norbp)


  call NonLocalHamiltonianApplication(iproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,rxyz,&
       proj,tmb%ham_descr%lzd,nlpspd,tmb%ham_descr%psi,lhphilarge,energs%eproj)

  ! only kinetic because waiting for communications
  call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
       tmb%ham_descr%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%ham_descr%psi,lhphilarge,&
       energs,SIC,GPU,3,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
  call full_local_potential(iproc,nproc,tmb%orbs,tmb%ham_descr%lzd,2,denspot%dpbox,denspot%rhov,denspot%pot_work, &
       tmb%ham_descr%comgp)
  ! only potential
  call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
       tmb%ham_descr%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%ham_descr%psi,lhphilarge,&
       energs,SIC,GPU,2,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmb%ham_descr%comgp)

  call timing(iproc,'glsynchham1','ON') !lr408t
  call SynchronizeHamiltonianApplication(nproc,tmb%ham_descr%npsidim_orbs,tmb%orbs,tmb%ham_descr%lzd,GPU,lhphilarge,&
       energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
  call timing(iproc,'glsynchham1','OF') !lr408t
  deallocate(confdatarrtmp)
  

  ! Now transpose the psi and hpsi
  allocate(lpsit_c(tmb%ham_descr%collcom%ndimind_c))
  call memocc(istat, lpsit_c, 'lpsit_c', subname)
  allocate(lpsit_f(7*tmb%ham_descr%collcom%ndimind_f))
  call memocc(istat, lpsit_f, 'lpsit_f', subname)
  allocate(hpsit_c(tmb%ham_descr%collcom%ndimind_c))
  call memocc(istat, hpsit_c, 'hpsit_c', subname)
  allocate(hpsit_f(7*tmb%ham_descr%collcom%ndimind_f))
  call memocc(istat, hpsit_f, 'hpsit_f', subname)
  allocate(psit_c(tmb%ham_descr%collcom%ndimind_c))
  call memocc(istat, psit_c, 'psit_c', subname)
  allocate(psit_f(7*tmb%ham_descr%collcom%ndimind_f))
  call memocc(istat, psit_f, 'psit_f', subname)

  call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
       tmb%ham_descr%psi, lpsit_c, lpsit_f, tmb%ham_descr%lzd)

  call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
       lhphilarge, hpsit_c, hpsit_f, tmb%ham_descr%lzd)

  !now build the derivative and related matrices <dPhi_a | H | Phi_b> and <dPhi_a | Phi_b>

  ! DOVRLP AND DHAM SHOULD HAVE DIFFERENT SPARSITIES, BUT TO MAKE LIFE EASIER KEEPING THEM THE SAME FOR NOW
  ! also array of structure a bit inelegant at the moment
  do jdir = 1, 3
    call nullify_sparsematrix(dovrlp(jdir))
    call nullify_sparsematrix(dham(jdir))
    call sparse_copy_pattern(tmb%linmat%ham,dovrlp(jdir),iproc,subname) 
    call sparse_copy_pattern(tmb%linmat%ham,dham(jdir),iproc,subname)
    allocate(dham(jdir)%matrix_compr(dham(jdir)%nvctr), stat=istat)
    call memocc(istat, dham(jdir)%matrix_compr, 'dham%matrix_compr', subname)
    allocate(dovrlp(jdir)%matrix_compr(dovrlp(jdir)%nvctr), stat=istat)
    call memocc(istat, dovrlp(jdir)%matrix_compr, 'dovrlp%matrix_compr', subname)

    call get_derivative(jdir, tmb%ham_descr%npsidim_orbs, tmb%ham_descr%lzd%hgrids(1), tmb%orbs, &
         tmb%ham_descr%lzd, tmb%ham_descr%psi, lhphilarge)

    call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
         lhphilarge, psit_c, psit_f, tmb%ham_descr%lzd)

    call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom,&
         psit_c, lpsit_c, psit_f, lpsit_f, dovrlp(jdir))

    call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom,&
         psit_c, hpsit_c, psit_f, hpsit_f, dham(jdir))
  end do


  !DEBUG
  !!print *,'iproc,tmb%orbs%norbp',iproc,tmb%orbs%norbp
  !!if(iproc==0)then
  !!do iorb = 1, tmb%orbs%norb
  !!   do iiorb=1,tmb%orbs%norb
  !!      !print *,'Hamiltonian of derivative: ',iorb, iiorb, (matrix(iorb,iiorb,jdir),jdir=1,3)
  !!      print *,'Overlap of derivative: ',iorb, iiorb, (dovrlp(iorb,iiorb,jdir),jdir=1,3)
  !!   end do
  !!end do
  !!end if
  !!!Check if derivatives are orthogonal to functions
  !!if(iproc==0)then
  !!  do iorb = 1, tmbder%orbs%norb
  !!     !print *,'overlap of derivative: ',iorb, (dovrlp(iorb,iiorb),iiorb=1,tmb%orbs%norb)
  !!     do iiorb=1,tmbder%orbs%norb
  !!         write(*,*) iorb, iiorb, dovrlp(iorb,iiorb)
  !!     end do
  !!  end do
  !!end if
  !END DEBUG

   ! needs generalizing if dovrlp and dham are to have different structures
   call to_zero(3*at%astruct%nat, fpulay(1,1))
   do jdir=1,3
     !do ialpha=1,tmb%orbs%norb
     if (tmb%orbs%norbp>0) then
         isegstart=dham(jdir)%istsegline(tmb%orbs%isorb_par(iproc)+1)
         if (tmb%orbs%isorb+tmb%orbs%norbp<tmb%orbs%norb) then
             isegend=dham(jdir)%istsegline(tmb%orbs%isorb_par(iproc+1)+1)-1
         else
             isegend=dham(jdir)%nseg
         end if
         do iseg=isegstart,isegend
              ii=dham(jdir)%keyv(iseg)-1
              do jorb=dham(jdir)%keyg(1,iseg),dham(jdir)%keyg(2,iseg)
                  ii=ii+1
                  iialpha = (jorb-1)/tmb%orbs%norb + 1
                  ibeta = jorb - (iialpha-1)*tmb%orbs%norb
                  jat=tmb%orbs%onwhichatom(iialpha)
                  kernel = 0.d0
                  ekernel= 0.d0
                  do iorb=1,orbs%norb
                      kernel  = kernel+orbs%occup(iorb)*tmb%coeff(iialpha,iorb)*tmb%coeff(ibeta,iorb)
                      ekernel = ekernel+tmb%orbs%eval(iorb)*orbs%occup(iorb) &
                           *tmb%coeff(iialpha,iorb)*tmb%coeff(ibeta,iorb) 
                  end do
                  fpulay(jdir,jat)=fpulay(jdir,jat)+&
                         2.0_gp*(kernel*dham(jdir)%matrix_compr(ii)-ekernel*dovrlp(jdir)%matrix_compr(ii))
              end do
         end do
     end if
   end do 

   call mpiallred(fpulay(1,1), 3*at%astruct%nat, mpi_sum, bigdft_mpi%mpi_comm, ierr)

  if(iproc==0) then
       do jat=1,at%astruct%nat
           write(*,'(a,i5,3es16.6)') 'iat, fpulay', jat, fpulay(1:3,jat)
       end do
  end if

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
  iall=-product(shape(lpsit_c))*kind(lpsit_c)
  deallocate(lpsit_c, stat=istat)
  call memocc(istat, iall, 'lpsit_c', subname)
  iall=-product(shape(lpsit_f))*kind(lpsit_f)
  deallocate(lpsit_f, stat=istat)
  call memocc(istat, iall, 'lpsit_f', subname)

  iall=-product(shape(lhphilarge))*kind(lhphilarge)
  deallocate(lhphilarge, stat=istat)
  call memocc(istat, iall, 'lhphilarge', subname)

  iall=-product(shape(denspot%pot_work))*kind(denspot%pot_work)
  deallocate(denspot%pot_work, stat=istat)
  call memocc(istat, iall, 'denspot%pot_work', subname)

  do jdir=1,3
     call deallocate_sparseMatrix(dovrlp(jdir),subname)
     call deallocate_sparseMatrix(dham(jdir),subname)
  end do

  if(iproc==0) write(*,'(1x,a)') 'done.'

end subroutine pulay_correction



!!subroutine derivative_coeffs_from_standard_coeffs(orbs, tmb, tmbder)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  type(orbitals_data),intent(in) :: orbs
!!  type(DFT_wavefunction),intent(in) :: tmb
!!  type(DFT_wavefunction),intent(out) :: tmbder
!!
!!  ! Local variables
!!  integer :: iorb, jorb, jjorb
!!
!!  call to_zero(tmbder%orbs%norb*orbs%norb, tmbder%coeff(1,1))
!!  do iorb=1,orbs%norb
!!      jjorb=0
!!      do jorb=1,tmbder%orbs%norb,4
!!          jjorb=jjorb+1
!!          tmbder%coeff(jorb,iorb)=tmb%coeff(jjorb,iorb)
!!      end do
!!  end do
!!
!!end subroutine derivative_coeffs_from_standard_coeffs



subroutine set_variables_for_hybrid(nlr, input, at, orbs, lowaccur_converged, confdatarr, &
           target_function, nit_basis, nit_scc, mix_hist, locrad, alpha_mix, convCritMix)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: nlr
  type(input_variables),intent(in) :: input
  type(atoms_data),intent(in) :: at
  type(orbitals_data),intent(in) :: orbs
  logical,intent(out) :: lowaccur_converged
  type(confpot_data),dimension(orbs%norbp),intent(inout) :: confdatarr
  integer,intent(out) :: target_function, nit_basis, nit_scc, mix_hist
  real(kind=8),dimension(nlr),intent(out) :: locrad
  real(kind=8),intent(out) :: alpha_mix, convCritMix

  ! Local variables
  integer :: iorb, ilr, iiat

  lowaccur_converged=.false.
  do iorb=1,orbs%norbp
      ilr=orbs%inwhichlocreg(orbs%isorb+iorb)
      iiat=orbs%onwhichatom(orbs%isorb+iorb)
      confdatarr(iorb)%prefac=input%lin%potentialPrefac_lowaccuracy(at%astruct%iatype(iiat))
  end do
  target_function=TARGET_FUNCTION_IS_HYBRID
  nit_basis=input%lin%nItBasis_lowaccuracy
  nit_scc=input%lin%nitSCCWhenFixed_lowaccuracy
  mix_hist=input%lin%mixHist_lowaccuracy
  do ilr=1,nlr
      locrad(ilr)=input%lin%locrad_lowaccuracy(ilr)
  end do
  alpha_mix=input%lin%alpha_mix_lowaccuracy
  convCritMix=input%lin%convCritMix_lowaccuracy

end subroutine set_variables_for_hybrid



! calculation of cSc and cHc using original coeffs (HOMO and LUMO only) and new Hamiltonian and overlap matrices
subroutine calc_transfer_integrals_old(iproc,nproc,input_frag,ref_frags,orbs,ham,ovrlp)
  use module_base
  use module_types
  use yaml_output
  use module_fragments
  use internal_io
  use module_interfaces
  implicit none

  integer, intent(in) :: iproc, nproc
  type(fragmentInputParameters), intent(in) :: input_frag
  type(system_fragment), dimension(input_frag%nfrag_ref), intent(in) :: ref_frags
  type(orbitals_data), intent(in) :: orbs
  type(sparseMatrix), intent(inout) :: ham, ovrlp
  !Local variables
  character(len=*), parameter :: subname='calc_transfer_integrals'
  integer :: i_stat, i_all, ifrag, jfrag, ntmb_tot, ind, itmb, ifrag_ref, ierr, ih, jh
  !integer :: jfrag_ref, jtmb
  integer, allocatable, dimension(:) :: homo
  real(gp), allocatable, dimension(:,:) :: homo_coeffs

  allocate(homo(input_frag%nfrag), stat=i_stat)
  call memocc(i_stat, homo, 'homo', subname)

  do ifrag=1,input_frag%nfrag
     ifrag_ref=input_frag%frag_index(ifrag)
     homo(ifrag)=ceiling(ref_frags(ifrag_ref)%nelec/2.0_gp)
  end do

  ntmb_tot=ham%full_dim1!=orbs%norb
  allocate(homo_coeffs(ntmb_tot,input_frag%nfrag), stat=i_stat)
  call memocc(i_stat, homo_coeffs, 'homo_coeffs', subname)

  if (input_frag%nfrag/=2) stop 'Error, only 2 fragments may currently be considered for transfer integral calculation'
  ! activate site energies only in case of more fragments

  if (iproc==0) write(*,*) 'HOMO and LUMO are defined as those of the neutral fragment'

  ! combine individual homo coeffs into a big ntmb_tot x input_frag%nfrag array

  ind=ref_frags(input_frag%frag_index(1))%fbasis%forbs%norb
  do ih=-1,2
     if (homo(input_frag%frag_index(1))+ih>ref_frags(input_frag%frag_index(1))%fbasis%forbs%norb) cycle

     call to_zero(ntmb_tot, homo_coeffs(1,1))

     do itmb=1,ref_frags(input_frag%frag_index(1))%fbasis%forbs%norb
        homo_coeffs(itmb,1)=ref_frags(input_frag%frag_index(1))%coeff(itmb,homo(input_frag%frag_index(1))+ih)
     end do

     do jh=-1,2
        if (homo(input_frag%frag_index(2))+jh>ref_frags(input_frag%frag_index(2))%fbasis%forbs%norb) cycle     

        call to_zero(ntmb_tot, homo_coeffs(1,2))

        do itmb=1,ref_frags(input_frag%frag_index(2))%fbasis%forbs%norb
           homo_coeffs(ind+itmb,2)=ref_frags(input_frag%frag_index(2))%coeff(itmb,homo(input_frag%frag_index(2))+jh)
        end do

        if (iproc==0) then
           if (ih<0) then
              write(*,'(a,I2)',advance='NO') 'Fragment 1 HOMO-',abs(ih)
           else if (ih==0) then
              write(*,'(a)',advance='NO') 'Fragment 1 HOMO'
           else if (ih==1) then
              write(*,'(a)',advance='NO') 'Fragment 1 LUMO'
           else
              write(*,'(a,I2)',advance='NO') 'Fragment 1 LUMO+',ih-1
           end if
        end if

        if (iproc==0) then
           if (jh<0) then
              write(*,'(a,I2,a)') ', fragment 2 HOMO-',abs(jh),'.  '
           else if (jh==0) then
              write(*,'(a)') ', fragment 2 HOMO.  '
           else if (jh==1) then
              write(*,'(a)') ', fragment 2 LUMO.  '
           else
              write(*,'(a,I2)') ', fragment 2 LUMO+',jh-1,'.  '
           end if
        end if

        call calc_transfer_integral_old(iproc,nproc,input_frag,orbs,ham,ovrlp,homo_coeffs)

     end do
  end do

  i_all = -product(shape(homo_coeffs))*kind(homo_coeffs)
  deallocate(homo_coeffs,stat=i_stat)
  call memocc(i_stat,i_all,'homo_coeffs',subname)

  i_all = -product(shape(homo))*kind(homo)
  deallocate(homo,stat=i_stat)
  call memocc(i_stat,i_all,'homo',subname)

end subroutine calc_transfer_integrals_old





! calculation of cSc and cHc using original coeffs (HOMO and LUMO only) and new Hamiltonian and overlap matrices
subroutine calc_transfer_integral_old(iproc,nproc,input_frag,orbs,ham,ovrlp,homo_coeffs)
  use module_base
  use module_types
  use yaml_output
  use module_fragments
  use internal_io
  use module_interfaces
  implicit none

  integer, intent(in) :: iproc, nproc
  type(fragmentInputParameters), intent(in) :: input_frag
  type(orbitals_data), intent(in) :: orbs
  type(sparseMatrix), intent(inout) :: ham, ovrlp
  real(kind=gp), dimension(ovrlp%full_dim1,input_frag%nfrag), intent(in) :: homo_coeffs
  !Local variables
  character(len=*), parameter :: subname='calc_transfer_integral'
  integer :: i_stat, i_all, ifrag, jfrag, ntmb_tot, ind, itmb, ierr, i, j
  !integer :: jfrag_ref, jtmb
  real(gp), allocatable, dimension(:,:) :: homo_ham, homo_ovrlp, coeff_tmp
  real(gp) :: orthog_energy


  ! make the coeff copies more efficient?

  allocate(coeff_tmp(orbs%norbp,input_frag%nfrag), stat=i_stat)
  call memocc(i_stat, coeff_tmp, 'coeff_tmp', subname)

  allocate(homo_ham(input_frag%nfrag,input_frag%nfrag), stat=i_stat)
  call memocc(i_stat, homo_ham, 'homo_ham', subname)
  allocate(homo_ovrlp(input_frag%nfrag,input_frag%nfrag), stat=i_stat)
  call memocc(i_stat, homo_ovrlp, 'homo_ovrlp', subname)

  allocate(ham%matrix(ham%full_dim1,ham%full_dim1), stat=i_stat)
  call memocc(i_stat, ham%matrix, 'ham%matrix', subname)
  call uncompressMatrix(iproc,ham)

  !DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  !rows op(a) and c, cols op(b) and c, cols op(a) and rows op(b)
  call to_zero(input_frag%nfrag**2, homo_ham(1,1))
  if (orbs%norbp>0) then
     call dgemm('n', 'n', orbs%norbp, input_frag%nfrag, orbs%norb, 1.d0, &
          ham%matrix(orbs%isorb+1,1),orbs%norb, &
          homo_coeffs(1,1), orbs%norb, 0.d0, &
          coeff_tmp, orbs%norbp)
     call dgemm('t', 'n', input_frag%nfrag, input_frag%nfrag, orbs%norbp, 1.d0, homo_coeffs(orbs%isorb+1,1), &
          orbs%norb, coeff_tmp, orbs%norbp, 0.d0, homo_ham, input_frag%nfrag)
  end if


  if (nproc>1) then
      call mpiallred(homo_ham(1,1), input_frag%nfrag**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if

  i_all=-product(shape(ham%matrix))*kind(ham%matrix)
  deallocate(ham%matrix, stat=i_stat)
  call memocc(i_stat, i_all, 'ham%matrix', subname)

  allocate(ovrlp%matrix(ovrlp%full_dim1,ovrlp%full_dim1), stat=i_stat)
  call memocc(i_stat, ovrlp%matrix, 'ovrlp%matrix', subname)
  call uncompressMatrix(iproc,ovrlp)

  call to_zero(input_frag%nfrag**2, homo_ovrlp(1,1))
  if (orbs%norbp>0) then
     call dgemm('n', 'n', orbs%norbp, input_frag%nfrag, orbs%norb, 1.d0, ovrlp%matrix(orbs%isorb+1,1), &
          orbs%norb, homo_coeffs(1,1), orbs%norb, 0.d0, coeff_tmp, orbs%norbp)
     call dgemm('t', 'n', input_frag%nfrag, input_frag%nfrag, orbs%norbp, 1.d0, homo_coeffs(orbs%isorb+1,1), &
          orbs%norb, coeff_tmp, orbs%norbp, 0.d0, homo_ovrlp, input_frag%nfrag)
  end if

  if (nproc>1) then
      call mpiallred(homo_ovrlp(1,1), input_frag%nfrag**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if

  i_all=-product(shape(ovrlp%matrix))*kind(ovrlp%matrix)
  deallocate(ovrlp%matrix, stat=i_stat)
  call memocc(i_stat, i_all, 'ovrlp%matrix', subname)

  i_all = -product(shape(coeff_tmp))*kind(coeff_tmp)
  deallocate(coeff_tmp,stat=i_stat)
  call memocc(i_stat,i_all,'coeff_tmp',subname)

  ! output results
  if (iproc==0) write(*,'(a)') '-----------------------------------------------------------------------------------------'
  if (input_frag%nfrag/=2) then
     if (iproc==0) write(*,*) 'Transfer integrals and site energies:'
     if (iproc==0) write(*,*) 'frag i, frag j, energy, overlap'
     do jfrag=1,input_frag%nfrag
        do ifrag=1,input_frag%nfrag
           if (iproc==0) write(*,'(2(I5,1x),1x,2(F16.12,1x))') jfrag, ifrag, homo_ham(jfrag,ifrag), homo_ovrlp(jfrag,ifrag)
        end do
     end do
  else ! include orthogonalized results as well
     !if (iproc==0) write(*,*) 'Site energies:'
     !if (iproc==0) write(*,*) 'frag i, energy, overlap, orthog energy'
     !i=1
     !j=2
     !orthog_energy= (0.5_gp/(1.0_gp-homo_ovrlp(i,j)**2)) &
     !             * ( (homo_ham(i,i)+homo_ham(j,j)) - 2.0_gp*homo_ham(i,j)*homo_ovrlp(i,j) &
     !             + (homo_ham(i,i)-homo_ham(j,j))*dsqrt(1.0_gp-homo_ovrlp(i,j)**2) )
     !if (iproc==0) write(*,'((I5,1x),1x,3(F16.12,1x))') 1, homo_ham(1,1), homo_ovrlp(1,1), orthog_energy
     !orthog_energy= (0.5_gp/(1.0_gp-homo_ovrlp(i,j)**2)) &
     !             * ( (homo_ham(i,i)+homo_ham(j,j)) - 2.0_gp*homo_ham(i,j)*homo_ovrlp(i,j) &
     !             - (homo_ham(i,i)-homo_ham(j,j))*dsqrt(1.0_gp-homo_ovrlp(i,j)**2) )
     !if (iproc==0) write(*,'((I5,1x),1x,3(F16.12,1x))') 2, homo_ham(2,2), homo_ovrlp(2,2), orthog_energy

     if (iproc==0) write(*,*) 'Transfer integrals:'
     if (iproc==0) write(*,*) 'frag i, frag j, energy, overlap, orthog energy'
     i=1
     j=2
     orthog_energy=(homo_ham(i,j)-0.5_gp*(homo_ham(i,i)+homo_ham(j,j))*homo_ovrlp(i,j))/(1.0_gp-homo_ovrlp(i,j)**2)
     if (iproc==0) write(*,'(2(I5,1x),1x,3(F16.12,1x))') 1, 2, homo_ham(1,2), homo_ovrlp(1,2),orthog_energy
     i=2
     j=1
     orthog_energy=(homo_ham(i,j)-0.5_gp*(homo_ham(i,i)+homo_ham(j,j))*homo_ovrlp(i,j))/(1.0_gp-homo_ovrlp(i,j)**2)
     if (iproc==0) write(*,'(2(I5,1x),1x,3(F16.12,1x))') 1, 2, homo_ham(2,1), homo_ovrlp(2,1),orthog_energy

  end if
  if (iproc==0) write(*,'(a)') '-----------------------------------------------------------------------------------------'

  i_all = -product(shape(homo_ham))*kind(homo_ham)
  deallocate(homo_ham,stat=i_stat)
  call memocc(i_stat,i_all,'homo_ham',subname)
  i_all = -product(shape(homo_ovrlp))*kind(homo_ovrlp)
  deallocate(homo_ovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'homo_ovrlp',subname)

end subroutine calc_transfer_integral_old


! calculation of cSc and cHc using original coeffs and new Hamiltonian and overlap matrices
! parallelization to be improved
! also have already uncompressed and recompressed ovrlp, so could change this
subroutine calc_transfer_integral(iproc,nproc,nstates,orbs,ham,ovrlp,homo_coeffs1,homo_coeffs2,homo_ham,homo_ovrlp)
  use module_base
  use module_types
  use yaml_output
  use module_fragments
  use internal_io
  use module_interfaces
  implicit none

  integer, intent(in) :: iproc, nproc, nstates
  type(orbitals_data), intent(in) :: orbs
  type(sparseMatrix), intent(inout) :: ham, ovrlp
  real(kind=gp), dimension(ovrlp%full_dim1,nstates), intent(in) :: homo_coeffs1, homo_coeffs2
  real(kind=gp), dimension(nstates), intent(inout) :: homo_ham, homo_ovrlp

  !Local variables
  integer :: i_stat, i_all, ifrag, jfrag, ntmb_tot, ind, itmb, ierr, i, j, istate
  real(gp), allocatable, dimension(:,:) :: coeff_tmp
  real(gp) :: orthog_energy

  coeff_tmp=f_malloc((/orbs%norbp,nstates/), id='coeff_tmp')

  !DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  !rows op(a) and c, cols op(b) and c, cols op(a) and rows op(b)
  !ham%matrix=f_malloc_ptr((/ham%full_dim1,ham%full_dim1/), id='ham%matrix')
  !call uncompressMatrix(iproc,ham)
  if (orbs%norbp>0) then
     do istate=1,nstates
        call dgemm('n', 'n', orbs%norbp, 1, orbs%norb, 1.d0, &
             ham%matrix(orbs%isorb+1,1),orbs%norb, &
             homo_coeffs1(1,istate), orbs%norb, 0.d0, &
             coeff_tmp(1,istate), orbs%norbp)
        call dgemm('t', 'n', 1, 1, orbs%norbp, 1.d0, homo_coeffs2(orbs%isorb+1,istate), &
             orbs%norb, coeff_tmp(1,istate), orbs%norbp, 0.d0, homo_ham(istate), 1)
     end do
  else
     call razero(nstates,homo_ham(1))
  end if

  if (nproc>1) then
      call mpiallred(homo_ham(1), nstates, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if

  !call f_free_ptr(ham%matrix)

  !ovrlp%matrix=f_malloc_ptr((/ovrlp%full_dim1,ovrlp%full_dim1/), id='ovrlp%matrix')
  !call uncompressMatrix(iproc,ovrlp)

  if (orbs%norbp>0) then
     do istate=1,nstates
        call dgemm('n', 'n', orbs%norbp, 1, orbs%norb, 1.d0, ovrlp%matrix(orbs%isorb+1,1), &
             orbs%norb, homo_coeffs1(1,istate), orbs%norb, 0.d0, coeff_tmp(1,istate), orbs%norbp)
        call dgemm('t', 'n', 1, 1, orbs%norbp, 1.d0, homo_coeffs2(orbs%isorb+1,istate), &
             orbs%norb, coeff_tmp(1,istate), orbs%norbp, 0.d0, homo_ovrlp(istate), 1)
     end do
  else
     call razero(nstates,homo_ovrlp(1))
  end if

  if (nproc>1) then
      call mpiallred(homo_ovrlp(1), nstates, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if

  !call f_free_ptr(ovrlp%matrix)
  call f_free(coeff_tmp)

end subroutine calc_transfer_integral


! calculation of cSc and cHc using original coeffs and new Hamiltonian and overlap matrices
! parallelization to be improved
! only calculates transfer integrals if we have two fragments
! occs are for neutral reference fragments...
subroutine calc_site_energies_transfer_integrals(iproc,nproc,input_frag,ref_frags,orbs,ham,ovrlp)
  use module_base
  use module_types
  use yaml_output
  use module_fragments
  use internal_io
  use module_interfaces
  implicit none

  integer, intent(in) :: iproc, nproc
  type(fragmentInputParameters), intent(in) :: input_frag
  type(orbitals_data), intent(in) :: orbs
  type(sparseMatrix), intent(inout) :: ham, ovrlp
  type(system_fragment), dimension(input_frag%nfrag_ref), intent(in) :: ref_frags
  !Local variables
  integer :: i_stat, i_all, ifrag, jfrag, ntmb_tot, ind, itmb, ierr, i, j, nstates, istate, ih, ifrag_ref
  integer :: ifrag_ref1, ifrag_ref2, homo1, homo2, jh, above_lumo, iind, jind, norb_tmp
  !integer :: jfrag_ref, jtmb
  real(gp), allocatable, dimension(:,:) :: coeffs_tmp, homo_coeffs, coeffs_orthog
  real(gp), allocatable, dimension(:) :: frag_sum, homo_ham, homo_ovrlp
  real(gp), allocatable, dimension(:) :: frag_sum_orthog, homo_ham_orthog, homo_ovrlp_orthog, eval_sum
  real(gp) :: frag_sum_tot, frag_sum_tot_orthog, eval_sum_tot, orthog_energy
  real(gp), dimension(1) :: trans_int_energy, trans_int_energy_orthog, trans_int_ovrlp, trans_int_ovrlp_orthog
  character(len=8) :: str

  call timing(iproc,'transfer_int','ON')
  call f_routine(id='calc_site_energies_transfer_integrals')

  nstates=0
  above_lumo=3
  do ifrag=1,input_frag%nfrag
     ifrag_ref= input_frag%frag_index(ifrag)
     nstates=nstates+min(ceiling((ref_frags(ifrag_ref)%nelec+1)/2.0_gp)+above_lumo,ref_frags(ifrag_ref)%fbasis%forbs%norb)
  end do

  homo_ham=f_malloc(nstates,id='homo_ham')
  homo_ovrlp=f_malloc(nstates,id='homo_ovrlp')
  homo_coeffs=f_malloc0((/ovrlp%full_dim1,nstates/), id='homo_coeffs')
  !coeffs_tmp=f_malloc((/ovrlp%full_dim1,ovrlp%full_dim1/), id='coeffs_tmp')
  ovrlp%matrix=f_malloc_ptr((/ovrlp%full_dim1,ovrlp%full_dim1/), id='ovrlp%matrix')
  call uncompressMatrix(iproc,ovrlp)

  istate=1
  ind=1
  do ifrag=1,input_frag%nfrag
     ifrag_ref=input_frag%frag_index(ifrag)

     norb_tmp=min(ceiling((ref_frags(ifrag_ref)%nelec+1)/2.0_gp)+above_lumo,ref_frags(ifrag_ref)%fbasis%forbs%norb)

     !call razero(ovrlp%full_dim1**2,coeffs_tmp(1,1),1)
     do ih=1,norb_tmp
        call dcopy(ref_frags(ifrag_ref)%fbasis%forbs%norb,ref_frags(ifrag_ref)%coeff(1,ih),1,homo_coeffs(ind,istate+ih-1),1)
        !call dcopy(ref_frags(ifrag_ref)%fbasis%forbs%norb,ref_frags(ifrag_ref)%coeff(1,ih),1,coeffs_tmp(ind,ih),1)
     end do

     !call reorthonormalize_coeff(iproc, nproc, norb_tmp, -8, -8, 0, orbs, ovrlp, coeffs_tmp(1,1))
     !call dcopy(orbs%norb*norb_tmp,coeffs_tmp(1,1),1,homo_coeffs(1,istate),1)

     istate=istate+norb_tmp
     ind=ind+ref_frags(ifrag_ref)%fbasis%forbs%norb
  end do
  !call f_free(coeffs_tmp)

  ham%matrix=f_malloc_ptr((/ham%full_dim1,ham%full_dim1/), id='ham%matrix')
  call uncompressMatrix(iproc,ham)
  call calc_transfer_integral(iproc,nproc,nstates,orbs,ham,ovrlp,homo_coeffs,homo_coeffs,homo_ham,homo_ovrlp)

  ! orthogonalize
  coeffs_tmp=f_malloc0((/orbs%norb,orbs%norb/), id='coeffs_orthog')
  call dcopy(orbs%norb*nstates,homo_coeffs(1,1),1,coeffs_tmp(1,1),1)
  call reorthonormalize_coeff(iproc, nproc, nstates, -8, -8, 0, orbs, ovrlp, coeffs_tmp(1,1))
  coeffs_orthog=f_malloc((/orbs%norb,nstates/), id='coeffs_orthog')
  call dcopy(orbs%norb*nstates,coeffs_tmp(1,1),1,coeffs_orthog(1,1),1)
  call f_free(coeffs_tmp)

  homo_ham_orthog=f_malloc(nstates, id='homo_ham_orthog')
  homo_ovrlp_orthog=f_malloc(nstates, id='homo_ovrlp_orthog')

  call calc_transfer_integral(iproc,nproc,nstates,orbs,ham,ovrlp,coeffs_orthog,coeffs_orthog,&
       homo_ham_orthog,homo_ovrlp_orthog)

  frag_sum=f_malloc0(nstates, id='frag_sum')
  frag_sum_orthog=f_malloc0(nstates, id='frag_sum_orthog')
  eval_sum=f_malloc0(nstates, id='eval_sum')

  if (iproc==0) write(*,'(a)') '-------------------------------------------------------------------------------------------------'
  if (iproc==0) write(*,*) 'Site energies:'

  if (iproc==0) write(*,*) 'state, energy, orthog energy, frag eval, overlap, orthog overlap, occ'
  istate=1
  frag_sum_tot=0
  frag_sum_tot_orthog=0
  eval_sum_tot=0
  do ifrag=1,input_frag%nfrag
     ifrag_ref=input_frag%frag_index(ifrag)
     if (iproc==0) write(*,'(a,i3)') trim(input_frag%label(ifrag_ref)),ifrag
     do ih=1,min(ceiling((ref_frags(ifrag_ref)%nelec+1)/2.0_gp)+above_lumo,ref_frags(ifrag_ref)%fbasis%forbs%norb)
        if (ih<ceiling(ref_frags(ifrag_ref)%nelec/2.0_gp)) then
           write(str,'(I2)') abs(ih-ceiling(ref_frags(ifrag_ref)%nelec/2.0_gp))
           if (iproc==0) write(*,'(a8)',advance='NO') ' HOMO-'//trim(adjustl(str))
        else if (ih==ceiling(ref_frags(ifrag_ref)%nelec/2.0_gp)) then
           if (iproc==0) write(*,'(a8)',advance='NO') ' HOMO'
        else if (ih==ceiling(ref_frags(ifrag_ref)%nelec/2.0_gp)+1) then
           if (iproc==0) write(*,'(a8)',advance='NO') ' LUMO'
        else
           write(str,'(I2)') ih-1-ceiling(ref_frags(ifrag_ref)%nelec/2.0_gp)
           if (iproc==0) write(*,'(a8)',advance='NO') ' LUMO+'//trim(adjustl(str))
        end if
        if (iproc==0) write(*,'(1x,5(F20.12,1x))',advance='NO') homo_ham(istate), homo_ham_orthog(istate), &
             ref_frags(ifrag_ref)%eval(ih), homo_ovrlp(istate), homo_ovrlp_orthog(istate)
        if (ih<ceiling(ref_frags(ifrag_ref)%nelec/2.0_gp)) then
           frag_sum(ifrag)=frag_sum(ifrag)+homo_ham(istate)
           frag_sum_orthog(ifrag)=frag_sum_orthog(ifrag)+homo_ham_orthog(istate)
           eval_sum(ifrag)=eval_sum(ifrag)+ref_frags(ifrag_ref)%eval(ih)
           if (iproc==0) write(*,'(1x,F4.2)') 2.0_gp
        else if (ih==ceiling(ref_frags(ifrag_ref)%nelec/2.0_gp)) then
           if (mod(real(ref_frags(ifrag_ref)%nelec,gp),2.0_gp)/=0.0_gp) then
              frag_sum(ifrag)=frag_sum(ifrag)+0.5_gp*homo_ham(istate)
              frag_sum_orthog(ifrag)=frag_sum_orthog(ifrag)+0.5_gp*homo_ham_orthog(istate)
              eval_sum(ifrag)=eval_sum(ifrag)+0.5_gp*ref_frags(ifrag_ref)%eval(ih)
              if (iproc==0) write(*,'(1x,F4.2)') 1.0_gp
           else
              frag_sum(ifrag)=frag_sum(ifrag)+homo_ham(istate)
              frag_sum_orthog(ifrag)=frag_sum_orthog(ifrag)+homo_ham_orthog(istate)
              eval_sum(ifrag)=eval_sum(ifrag)+ref_frags(ifrag_ref)%eval(ih)
              if (iproc==0) write(*,'(1x,F4.2)') 2.0_gp
           end if
        else
           if (iproc==0) write(*,'(1x,F4.2)') 0.0_gp
        end if
        istate=istate+1
     end do
     if (iproc==0) write(*,'(9x,3(F20.12,1x))') 2.0_gp*frag_sum(ifrag),&
          2.0_gp*frag_sum_orthog(ifrag),2.0_gp*eval_sum(ifrag)
       if (iproc==0) write(*,'(a)') '------------------------------------------------------------------------'//&
            '-------------------------'
     frag_sum_tot=frag_sum_tot+frag_sum(ifrag)
     frag_sum_tot_orthog=frag_sum_tot_orthog+frag_sum_orthog(ifrag)
     eval_sum_tot=eval_sum_tot+eval_sum(ifrag)
  end do

  if (iproc==0) write(*,'(9x,3(F20.12,1x))') 2.0_gp*frag_sum_tot, 2.0_gp*frag_sum_tot_orthog,2.0_gp*eval_sum_tot
  if (iproc==0) write(*,'(a)') '-------------------------------------------------------------------------------------------------'

  call f_free(eval_sum)
  call f_free(frag_sum)
  call f_free(frag_sum_orthog)
  call f_free(homo_ham_orthog)
  call f_free(homo_ovrlp_orthog)

  if (input_frag%nfrag>=2) then
     if (iproc==0) write(*,*) 'Transfer integrals (HOMO and LUMO are defined as those of the neutral fragment):'
     if (iproc==0) write(*,*) 'state1, state2, energy, orthog energy, orthog energy2, overlap, orthog overlap, occ1, occ2'
     iind=0
     do ifrag=1,input_frag%nfrag
        ifrag_ref1=input_frag%frag_index(ifrag)
        homo1=ceiling((ref_frags(ifrag_ref1)%nelec)/2.0_gp)

        jind=0
        do jfrag=1,ifrag
           ifrag_ref2=input_frag%frag_index(jfrag)
           homo2=ceiling((ref_frags(ifrag_ref2)%nelec)/2.0_gp)

           do jh=-above_lumo,1+above_lumo
              if (homo2+jh>ref_frags(ifrag_ref2)%fbasis%forbs%norb) cycle  
              if (homo2+jh<1) cycle  
              do ih=-above_lumo,1+above_lumo
                 if (homo1+ih>ref_frags(ifrag_ref1)%fbasis%forbs%norb) cycle
                 if (homo1+ih<1) cycle  

                 i=homo1+ih+iind
                 j=homo2+jh+jind
                      
                 if (iproc==0) then
                    if (ih<0) then
                       write(str,'(I2)') abs(ih)
                       write(*,'(a,I3,a8)',advance='NO') trim(input_frag%label(ifrag_ref1)),ifrag,' HOMO-'//trim(adjustl(str))
                    else if (ih==0) then
                       write(*,'(a,I3,a8)',advance='NO') trim(input_frag%label(ifrag_ref1)),ifrag,' HOMO  '
                    else if (ih==1) then
                       write(*,'(a,I3,a8)',advance='NO') trim(input_frag%label(ifrag_ref1)),ifrag,' LUMO  '
                    else
                       write(str,'(I2)') ih-1
                       write(*,'(a,I3,a8)',advance='NO') trim(input_frag%label(ifrag_ref1)),ifrag,' LUMO+'//trim(adjustl(str))
                    end if
                 end if

                 if (iproc==0) then
                    if (jh<0) then
                       write(str,'(I2)') abs(jh)
                       write(*,'(3x,a,I3,a8)',advance='NO') trim(input_frag%label(ifrag_ref2)),jfrag,&
                            ' HOMO-'//trim(adjustl(str))
                    else if (jh==0) then
                       write(*,'(3x,a,I3,a8)',advance='NO') trim(input_frag%label(ifrag_ref2)),jfrag,' HOMO  '
                    else if (jh==1) then
                       write(*,'(3x,a,I3,a8)',advance='NO') trim(input_frag%label(ifrag_ref2)),jfrag,' LUMO  '
                    else
                       write(str,'(I2)') jh-1
                       write(*,'(3x,a,I3,a8)',advance='NO') trim(input_frag%label(ifrag_ref2)),jfrag,&
                            ' LUMO+'//trim(adjustl(str))
                    end if
                 end if

                 call calc_transfer_integral(iproc,nproc,1,orbs,ham,ovrlp,homo_coeffs(1,i),homo_coeffs(1,j),&
                      trans_int_energy(1),trans_int_ovrlp(1))
                 call calc_transfer_integral(iproc,nproc,1,orbs,ham,ovrlp,coeffs_orthog(1,i),coeffs_orthog(1,j),&
                      trans_int_energy_orthog(1),trans_int_ovrlp_orthog(1))

                 orthog_energy=(trans_int_energy(1)-0.5_gp*(homo_ham(i)+homo_ham(j))*trans_int_ovrlp(1))&
                      /(1.0_gp-trans_int_ovrlp(1)**2)
      
                 if (iproc==0) write(*,'(2x,5(F16.12,1x))',advance='NO') trans_int_energy(1), &
                      trans_int_energy_orthog(1), orthog_energy, trans_int_ovrlp(1), trans_int_ovrlp_orthog(1)

                 if (homo1+ih<ceiling(ref_frags(ifrag_ref1)%nelec/2.0_gp)) then
                    if (iproc==0) write(*,'(1x,F4.2)',advance='NO') 2.0_gp
                 else if (homo1+ih==ceiling(ref_frags(ifrag_ref1)%nelec/2.0_gp)) then
                    if (mod(real(ref_frags(ifrag_ref1)%nelec,gp),2.0_gp)/=0.0_gp) then
                       if (iproc==0) write(*,'(1x,F4.2)',advance='NO') 1.0_gp
                    else
                       if (iproc==0) write(*,'(1x,F4.2)',advance='NO') 2.0_gp
                    end if
                 else
                    if (iproc==0) write(*,'(1x,F4.2)',advance='NO') 0.0_gp
                 end if

                 if (homo2+jh<ceiling(ref_frags(ifrag_ref2)%nelec/2.0_gp)) then
                    if (iproc==0) write(*,'(1x,F4.2)') 2.0_gp
                 else if (homo2+jh==ceiling(ref_frags(ifrag_ref2)%nelec/2.0_gp)) then
                    if (mod(real(ref_frags(ifrag_ref2)%nelec,gp),2.0_gp)/=0.0_gp) then
                       if (iproc==0) write(*,'(1x,F4.2)') 1.0_gp
                    else
                       if (iproc==0) write(*,'(1x,F4.2)') 2.0_gp
                    end if
                 else
                    if (iproc==0) write(*,'(1x,F4.2)') 0.0_gp
                 end if

              end do
           end do
           if (iproc==0) write(*,'(a)') '------------------------------------------------------------------------'//&
               '-------------------------'
           jind=jind+min(ceiling((ref_frags(ifrag_ref2)%nelec+1)/2.0_gp)+above_lumo,ref_frags(ifrag_ref2)%fbasis%forbs%norb)
        end do
        iind=iind+min(ceiling((ref_frags(ifrag_ref1)%nelec+1)/2.0_gp)+above_lumo,ref_frags(ifrag_ref1)%fbasis%forbs%norb)
     end do
  end if

  call f_free_ptr(ham%matrix)
  call f_free_ptr(ovrlp%matrix)

  call f_free(homo_ham)
  call f_free(homo_ovrlp)
  call f_free(homo_coeffs)
  call f_free(coeffs_orthog)

  call f_release_routine()
  call timing(iproc,'transfer_int','OF')

end subroutine calc_site_energies_transfer_integrals



subroutine output_fragment_rotations(iproc,nproc,nat,rxyz,iformat,filename,input_frag,ref_frags)
  use module_base
  use module_types
  use yaml_output
  use module_fragments
  use internal_io
  use module_interfaces
  implicit none

  integer, intent(in) :: iproc, nproc, iformat, nat
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

           allocate(rxyz_ref(3,ref_frags(ifrag_ref)%astruct_frg%nat), stat=i_stat)
           call memocc(i_stat, rxyz_ref, 'rxyz_ref', subname)
           allocate(rxyz_new(3,ref_frags(ifrag_ref)%astruct_frg%nat), stat=i_stat)
           call memocc(i_stat, rxyz_new, 'rxyz_ref', subname)

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

           i_all = -product(shape(rxyz_ref))*kind(rxyz_ref)
           deallocate(rxyz_ref,stat=i_stat)
           call memocc(i_stat,i_all,'rxyz_ref',subname)
           i_all = -product(shape(rxyz_new))*kind(rxyz_new)
           deallocate(rxyz_new,stat=i_stat)
           call memocc(i_stat,i_all,'rxyz_new',subname)

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
