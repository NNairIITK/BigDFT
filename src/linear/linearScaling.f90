subroutine linearScaling(iproc,nproc,KSwfn,tmb,tmblarge,at,input,&
           rxyz,denspot,rhopotold,nlpspd,proj,GPU,&
           energs,energy,fpulay,infocode)
 
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => linearScaling
  use yaml_output
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(atoms_data),intent(inout) :: at
  type(input_variables),intent(in) :: input
  real(8),dimension(3,at%nat),intent(inout) :: rxyz
  real(8),dimension(3,at%nat),intent(out) :: fpulay
  type(DFT_local_fields), intent(inout) :: denspot
  real(gp), dimension(:), intent(inout) :: rhopotold
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd
  real(wp),dimension(nlpspd%nprojel),intent(inout) :: proj
  type(GPU_pointers),intent(in out) :: GPU
  type(energy_terms),intent(inout) :: energs
  real(gp), dimension(:), pointer :: rho,pot
  real(8),intent(out) :: energy
  type(DFT_wavefunction),intent(inout),target :: tmb, tmblarge
  type(DFT_wavefunction),intent(inout),target :: KSwfn
  integer,intent(out) :: infocode
  
  real(8) :: pnrm,trace,trace_old,fnrm_tmb
  integer :: infoCoeff,istat,iall,it_scc,itout,info_scf,nsatur,i,ierr,it_coeff_opt
  character(len=*),parameter :: subname='linearScaling'
  real(8),dimension(:),allocatable :: rhopotold_out
  real(8) :: energyold, energyDiff, energyoldout, fnrm_pulay, convCritMix
  type(mixrhopotDIISParameters) :: mixdiis
  type(localizedDIISParameters) :: ldiis, ldiis_coeff
  logical :: can_use_ham, update_phi, locreg_increased
  logical :: fix_support_functions, check_initialguess
  integer :: itype, istart, nit_lowaccuracy, nit_highaccuracy
  real(8),dimension(:),allocatable :: locrad_tmp, ham_compr, overlapmatrix_compr
  !!real(8),dimension(:),allocatable :: rhotest
  !!real(kind=8),dimension(:,:),allocatable :: density_kernel
  integer :: ldiis_coeff_hist
  logical :: ldiis_coeff_changed
  integer :: mix_hist, info_basis_functions, nit_scc, cur_it_highaccuracy
  real(8) :: pnrm_out, alpha_mix
  logical :: lowaccur_converged, exit_outer_loop
  real(8),dimension(:),allocatable :: locrad
  integer:: target_function, nit_basis

  call timing(iproc,'linscalinit','ON') !lr408t

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
  if(input%lin%mixHist_lowaccuracy>0 .and. input%lin%scf_mode/=LINEAR_DIRECT_MINIMIZATION) then
      call initializeMixrhopotDIIS(input%lin%mixHist_lowaccuracy, denspot%dpbox%ndimpot, mixdiis)
  end if

  !!!! TEST #######################################
  !!    call init_collective_comms_sumro(iproc, nproc, tmb%lzd, tmb%orbs, collcom_sr)
  !!!! END TEST ###################################

  pnrm=1.d100
  pnrm_out=1.d100
  energyold=0.d0
  energyoldout=0.d0
  target_function=TARGET_FUNCTION_IS_TRACE
  lowaccur_converged=.false.
  info_basis_functions=-1
  nsatur=0
  fix_support_functions=.false.
  check_initialguess=.true.
  cur_it_highaccuracy=0
  trace_old=0.0d0
  ldiis_coeff_hist=input%lin%mixHist_lowaccuracy
  it_coeff_opt=0

  ! Allocate the communication arrays for the calculation of the charge density.
  !!call allocateCommunicationbufferSumrho(iproc, tmb%comsr, subname)

  ! take the eigenvalues from the input guess for the preconditioning 
  call vcopy(tmb%orbs%norb, tmb%orbs%eval(1), 1, tmblarge%orbs%eval(1), 1)

  call timing(iproc,'linscalinit','OF') !lr408t

  if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then  
     call initialize_DIIS_coeff(ldiis_coeff_hist, ldiis_coeff)
     call allocate_DIIS_coeff(tmb, KSwfn%orbs, ldiis_coeff)
  end if

  ! Should be removed by passing tmblarge to restart
  !!if(input%inputPsiId  == INPUT_PSI_MEMORY_LINEAR .or. input%inputPsiId  == INPUT_PSI_DISK_LINEAR) then
  !!   call create_large_tmbs(iproc, nproc, tmb, denspot, input, at, rxyz, lowaccur_converged, &
  !!        tmblarge)

  !!   ! Set to zero the large wavefunction. Later only the inner part will be filled. It must be made sure
  !!   ! that the outer part is not modified!
  !!   if (tmblarge%orbs%npsidim_orbs > 0) call to_zero(tmblarge%orbs%npsidim_orbs,tmblarge%psi(1))
  !!end if

  ! Orthogonalize the input guess minimal basis functions using exact calculation of S^-1/2
  tmb%can_use_transposed=.false.
  nullify(tmb%psit_c)
  nullify(tmb%psit_f)
  if(iproc==0) write(*,*) 'calling orthonormalizeLocalized (exact)'

  ! Give tmblarge%mad since this is the correct matrix description
  call orthonormalizeLocalized(iproc, nproc, -1, tmb%orbs, tmb%lzd, tmblarge%mad, tmb%collcom, &
       tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)

  ! Check the quality of the input guess
  call check_inputguess()

  if(input%lin%scf_mode/=LINEAR_MIXPOT_SIMPLE) then
      call dcopy(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim),rhopotold(1),1,rhopotold_out(1),1)
  else
      call dcopy(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim),denspot%rhov(1),1,rhopotold_out(1),1)
      call dcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p,1) &
            *input%nspin, denspot%rhov(1), 1, rhopotOld(1), 1)
  end if

  if (iproc==0) call yaml_open_map('Checking Communications of Minimal Basis')
  call check_communications_locreg(iproc,nproc,tmb%orbs,tmb%Lzd,tmb%collcom)
  if (iproc==0) call yaml_close_map()

  if (iproc==0) call yaml_open_map('Checking Communications of Enlarged Minimal Basis')
  call check_communications_locreg(iproc,nproc,tmblarge%orbs,&
       tmblarge%Lzd,tmblarge%collcom)
  if (iproc ==0) call yaml_close_map()


  ! Add one iteration if no low accuracy is desired since we need then a first fake iteration, with istart=0
  istart = min(1,nit_lowaccuracy)
  infocode=0 !default value
  ! This is the main outer loop. Each iteration of this loop consists of a first loop in which the basis functions
  ! are optimized and a consecutive loop in which the density is mixed.
  outerLoop: do itout=istart,nit_lowaccuracy+nit_highaccuracy

      ! Check whether the low accuracy part (i.e. with strong confining potential) has converged.
      call check_whether_lowaccuracy_converged(itout, nit_lowaccuracy, input%lin%lowaccuracy_conv_crit, &
           lowaccur_converged, pnrm_out)
      ! Set all remaining variables that we need for the optimizations of the basis functions and the mixing.
      call set_optimization_variables(input, at, tmb%orbs, tmb%lzd%nlr, tmb%orbs%onwhichatom, tmb%confdatarr, &
           convCritMix, lowaccur_converged, nit_scc, mix_hist, alpha_mix, locrad, target_function, nit_basis)

      ! Do one fake iteration if no low accuracy is desired.
      if(nit_lowaccuracy==0 .and. itout==0) then
          lowaccur_converged=.false.
          cur_it_highaccuracy=0
      end if

      if(lowaccur_converged) cur_it_highaccuracy=cur_it_highaccuracy+1

      if(cur_it_highaccuracy==1) then
          ! Adjust the confining potential if required.
          call adjust_locregs_and_confinement(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
               at, input, tmb, denspot, ldiis, locreg_increased, lowaccur_converged, locrad)
          ! Reajust tmblarge also
          if(locreg_increased) then
             call destroy_new_locregs(iproc, nproc, tmblarge)
             call deallocate_auxiliary_basis_function(subname, tmblarge%psi, tmblarge%hpsi)
             if(tmblarge%can_use_transposed) then
                 iall=-product(shape(tmblarge%psit_c))*kind(tmblarge%psit_c)
                 deallocate(tmblarge%psit_c, stat=istat)
                 call memocc(istat, iall, 'tmblarge%psit_c', subname)
                 iall=-product(shape(tmblarge%psit_f))*kind(tmblarge%psit_f)
                 deallocate(tmblarge%psit_f, stat=istat)
                 call memocc(istat, iall, 'tmblarge%psit_f', subname)
             end if
             deallocate(tmblarge%confdatarr, stat=istat)

             iall=-product(shape(ham_compr))*kind(ham_compr)
             deallocate(ham_compr, stat=istat)
             call memocc(istat, iall, 'ham_compr', subname)
             iall=-product(shape(overlapmatrix_compr))*kind(overlapmatrix_compr)
             deallocate(overlapmatrix_compr, stat=istat)
             call memocc(istat, iall, 'overlapmatrix_compr', subname)

             call create_large_tmbs(iproc, nproc, tmb, denspot, input, at, rxyz, lowaccur_converged, &
                  tmblarge)
             call init_collective_comms(iproc, nproc, tmb%orbs, tmb%lzd, tmblarge%mad, tmb%collcom)
             call init_collective_comms(iproc, nproc, tmblarge%orbs, tmblarge%lzd, tmblarge%mad, tmblarge%collcom)
             call init_collective_comms_sumro(iproc, nproc, tmb%lzd, tmb%orbs, tmblarge%mad, &
                  denspot%dpbox%nscatterarr, tmb%collcom_sr)

             allocate(ham_compr(tmblarge%mad%nvctr), stat=istat)
             call memocc(istat, ham_compr, 'ham_compr', subname)
             allocate(overlapmatrix_compr(tmblarge%mad%nvctr), stat=istat)
             call memocc(istat, overlapmatrix_compr, 'overlapmatrix_compr', subname)

          else
             call define_confinement_data(tmblarge%confdatarr,tmblarge%orbs,rxyz,at,&
                   tmblarge%lzd%hgrids(1),tmblarge%lzd%hgrids(2),tmblarge%lzd%hgrids(3),&
                   4,input%lin%potentialPrefac_highaccuracy,tmblarge%lzd,tmblarge%orbs%onwhichatom)
          end if
          ! Calculate a new kernel since the old compressed one has changed its shape due to the locrads
          ! being different for low and high accuracy.
          update_phi=.true.
          tmb%can_use_transposed=.false.   !check if this is set properly!
          call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
               infoCoeff,energs%ebs,nlpspd,proj,input%SIC,tmb,pnrm,update_phi,update_phi,&
               tmblarge,ham_compr,overlapmatrix_compr,.true.,it_coeff_opt,ldiis_coeff=ldiis_coeff)
      end if

      ! Some special treatement if we are in the high accuracy part
      call adjust_DIIS_for_high_accuracy(input, denspot, mixdiis, lowaccur_converged, &
           ldiis_coeff_hist, ldiis_coeff_changed)

      if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then 
         call initialize_DIIS_coeff(ldiis_coeff_hist, ldiis_coeff)
         ! need to reallocate DIIS matrices to adjust for changing history length
         if (ldiis_coeff_changed) then
            call deallocateDIIS(ldiis_coeff)
            call allocate_DIIS_coeff(tmb, KSwfn%orbs, ldiis_coeff)
         end if
      end if

      if(itout>1 .or. (nit_lowaccuracy==0 .and. itout==1)) then
          call deallocateDIIS(ldiis)
      end if
      if (lowaccur_converged) then
          call initializeDIIS(input%lin%DIIS_hist_highaccur, tmb%lzd, tmb%orbs, ldiis)
      else
          call initializeDIIS(input%lin%DIIS_hist_lowaccur, tmb%lzd, tmb%orbs, ldiis)
      end if
      ldiis%DIISHistMin=0
      if (lowaccur_converged) then
          ldiis%DIISHistMax=input%lin%DIIS_hist_highaccur
      else
          ldiis%DIISHistMax=input%lin%DIIS_hist_lowaccur
      end if
      ldiis%switchSD=.false.
      ldiis%trmin=1.d100
      ldiis%trold=1.d100
      ldiis%icountSDSatur=0
      ldiis%icountSwitch=0
      ldiis%icountDIISFailureTot=0
      ldiis%icountDIISFailureCons=0
      ldiis%is=0
      ldiis%switchSD=.false.
      ldiis%trmin=1.d100
      ldiis%trold=1.d100
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
           !!call uncompressMatrix(tmb%orbs%norb, tmblarge%mad, tmb%wfnmd%density_kernel_compr, tmb%wfnmd%density_kernel)
           !!do istat=1,tmb%orbs%norb
           !!    do iall=1,tmb%orbs%norb
           !!        !write(200+iproc,*) istat, iall, tmb%wfnmd%density_kernel(iall,istat)
           !!        !read(200+iproc,*) iorb, iiorb, tmb%wfnmd%density_kernel(iall,istat)
           !!    end do
           !!end do 
           call getLocalizedBasis(iproc,nproc,at,KSwfn%orbs,rxyz,denspot,GPU,trace,trace_old,fnrm_tmb,&
               info_basis_functions,nlpspd,input%lin%scf_mode,proj,ldiis,input%SIC,tmb,tmblarge,energs, &
               ham_compr,input%lin%nItPrecond,target_function,input%lin%correctionOrthoconstraint,nit_basis)
           if(info_basis_functions>0) then
               nsatur=nsatur+1
           end if
           tmb%can_use_transposed=.false. !since basis functions have changed...

           it_coeff_opt=0
           if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) ldiis_coeff%alpha_coeff=0.2d0 !reset to default value

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
           do itype=1,at%ntypes
               if(input%lin%potentialPrefac_lowaccuracy(itype)/=0.d0) then
                   can_use_ham=.false.
                   exit
               end if
           end do
       else if(target_function==TARGET_FUNCTION_IS_ENERGY) then
           do itype=1,at%ntypes
               if(input%lin%potentialPrefac_highaccuracy(itype)/=0.d0) then
                   can_use_ham=.false.
                   exit
               end if
           end do
       end if
  
       !ADDITIONAL MIXING BEFORE GET_COEFF
       !if (itout>1.and..false.) then
       !   !call reconstruct_kernel(iproc, nproc, 1, tmb%orthpar%blocksize_pdsyev, tmb%orthpar%blocksize_pdgemm, &
       !   !     KSwfn%orbs, tmb, overlapmatrix, overlap_calculated, tmb%wfnmd%density_kernel) 
       !
       !   call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
       !        tmb%orbs, tmb%collcom_sr, tmb%wfnmd%density_kernel, &
       !        KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
       !!if (.false.) then
       !! Mix the density.
       !if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE .or. input%lin%scf_mode==LINEAR_FOE) then
       !   lscv%compare_outer_loop = pnrm<input%lin%convCritMix .or. it_scc==nit_scc
       !   call mix_main(iproc, nproc, mix_hist, lscv%compare_outer_loop, input, KSwfn%Lzd%Glr, alpha_mix, &
       !        denspot, mixdiis, rhopotold, rhopotold_out, pnrm, pnrm_out)
       !end if
       !!end if
       !  
       !! Calculate the new potential.
       !if(iproc==0) write(*,'(1x,a)') '---------------------------------------------------------------- Updating potential.'
       !call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
       !if (.false.) then
       !! Mix the potential
       !if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
       !   lscv%compare_outer_loop = pnrm<input%lin%convCritMix .or. it_scc==nit_scc
       !   call mix_main(iproc, nproc, mix_hist, lscv%compare_outer_loop, input, KSwfn%Lzd%Glr, alpha_mix, &
       !        denspot, mixdiis, rhopotold, rhopotold_out, pnrm, pnrm_out)
       !end if 
       !end if
       !end if

      ! The self consistency cycle. Here we try to get a self consistent density/potential with the fixed basis.
      !call init_collective_comms(iproc, nproc, tmb%orbs, tmb%lzd, tmblarge%mad, tmb%collcom)
      kernel_loop : do it_scc=1,nit_scc

          ! If the hamiltonian is available do not recalculate it
          ! also using update_phi for calculate_overlap_matrix and communicate_phi_for_lsumrho
          ! since this is only required if basis changed
          if(update_phi .and. can_use_ham .and. info_basis_functions>=0) then
              call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                   infoCoeff,energs%ebs,nlpspd,proj,input%SIC,tmb,pnrm,update_phi,update_phi,&
                   tmblarge,ham_compr,overlapmatrix_compr,.false.,it_coeff_opt,ldiis_coeff=ldiis_coeff)
          else
              call get_coeff(iproc,nproc,input%lin%scf_mode,KSwfn%orbs,at,rxyz,denspot,GPU,&
                   infoCoeff,energs%ebs,nlpspd,proj,input%SIC,tmb,pnrm,update_phi,update_phi,&
                   tmblarge,ham_compr,overlapmatrix_compr,.true.,it_coeff_opt,ldiis_coeff=ldiis_coeff)
          end if

          ! Since we do not update the basis functions anymore in this loop
          update_phi = .false.

          ! Calculate the total energy.
          !if(iproc==0) print *,'energs',energs%ebs,energs%eh,energs%exc,energs%evxc,energs%eexctX,energs%eion,energs%edisp
          energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
          energyDiff=energy-energyold
          energyold=energy

          ! Calculate the charge density.
          !!call sumrhoForLocalizedBasis2(iproc, nproc, &
          !!     tmb%lzd, tmb%orbs, tmb%comsr, &
          !!     tmb%wfnmd%density_kernel, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, &
          !!     denspot%rhov, at, denspot%dpbox%nscatterarr)
          !!allocate(rhotest(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d), stat=istat)
          !!allocate(density_kernel(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
          !!call memocc(istat, density_kernel, 'density_kernel', subname)
          !!call uncompressMatrix(tmb%orbs%norb, tmblarge%mad, tmb%wfnmd%density_kernel_compr, density_kernel)
          call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
               tmb%orbs, tmblarge%mad, tmb%collcom_sr, tmb%wfnmd%density_kernel_compr, &
               KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
          !!iall=-product(shape(density_kernel))*kind(density_kernel)
          !!deallocate(density_kernel, stat=istat)
          !!call memocc(istat, iall, 'density_kernel', subname)

          !!do istat=1,KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d
          !!    write(1000+iproc,*) istat, denspot%rhov(istat)
          !!    write(2000+iproc,*) istat, rhotest(istat)
          !!end do
          !!deallocate(rhotest, stat=istat)

          ! Mix the density.
          if (input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE .or. input%lin%scf_mode==LINEAR_FOE) then
             call mix_main(iproc, nproc, mix_hist, input, KSwfn%Lzd%Glr, alpha_mix, &
                  denspot, mixdiis, rhopotold, pnrm)
          end if
 
          if (input%lin%scf_mode/=LINEAR_MIXPOT_SIMPLE .and.(pnrm<convCritMix .or. it_scc==nit_scc)) then
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

          ! Calculate the new potential.
          if(iproc==0) write(*,'(1x,a)') '---------------------------------------------------------------- Updating potential.'
          call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)

          ! Mix the potential
          if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
             call mix_main(iproc, nproc, mix_hist, input, KSwfn%Lzd%Glr, alpha_mix, &
                  denspot, mixdiis, rhopotold, pnrm)
             if (pnrm<convCritMix .or. it_scc==nit_scc) then
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

      if(tmb%can_use_transposed) then
          iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
          deallocate(tmb%psit_c, stat=istat)
          call memocc(istat, iall, 'tmb%psit_c', subname)
          iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
          deallocate(tmb%psit_f, stat=istat)
          call memocc(istat, iall, 'tmb%psit_f', subname)
      end if

      call print_info()

      energyoldout=energy

      call check_for_exit()
      if(exit_outer_loop) exit outerLoop

      if(pnrm_out<input%lin%support_functions_converged.and.lowaccur_converged) then
          if(iproc==0) write(*,*) 'fix the support functions from now on'
          fix_support_functions=.true.
      end if

  end do outerLoop

  ! Diagonalize the matrix for the FOE case to get the coefficients. Only necessary if
  ! the Pulay forces are to be calculated.
  if (input%lin%scf_mode==LINEAR_FOE .and. input%lin%pulay_correction) then
      call get_coeff(iproc,nproc,LINEAR_MIXDENS_SIMPLE,KSwfn%orbs,at,rxyz,denspot,GPU,&
           infoCoeff,energs%ebs,nlpspd,proj,input%SIC,tmb,pnrm,update_phi,.false.,&
           tmblarge,ham_compr,overlapmatrix_compr,.true.,it_coeff_opt,ldiis_coeff=ldiis_coeff)
  end if

  ! Deallocate everything that is not needed any more.
  if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) call deallocateDIIS(ldiis_coeff)
  call deallocateDIIS(ldiis)
  if(input%lin%mixHist_highaccuracy>0 .and. input%lin%scf_mode/=LINEAR_DIRECT_MINIMIZATION) then
      call deallocateMixrhopotDIIS(mixdiis)
  end if
  !!call wait_p2p_communication(iproc, nproc, tmb%comgp)
  call synchronize_onesided_communication(iproc, nproc, tmb%comgp)
  call deallocateCommunicationsBuffersPotential(tmb%comgp, subname)


  if (input%lin%pulay_correction) then
      if (iproc==0) write(*,'(1x,a)') 'WARNING: commented correction_locrad!'
      !!! Testing energy corrections due to locrad
      !!call correction_locrad(iproc, nproc, tmblarge, KSwfn%orbs,tmb%wfnmd%coeff) 
      ! Calculate Pulay correction to the forces
      call pulay_correction(iproc, nproc, KSwfn%orbs, at, rxyz, nlpspd, proj, input%SIC, denspot, GPU, tmb, &
           tmblarge, fpulay)
  else
      call to_zero(3*at%nat, fpulay(1,1))
  end if

  if(tmblarge%can_use_transposed) then
      iall=-product(shape(tmblarge%psit_c))*kind(tmblarge%psit_c)
      deallocate(tmblarge%psit_c, stat=istat)
      call memocc(istat, iall, 'tmblarge%psit_c', subname)
      iall=-product(shape(tmblarge%psit_f))*kind(tmblarge%psit_f)
      deallocate(tmblarge%psit_f, stat=istat)
      call memocc(istat, iall, 'tmblarge%psit_f', subname)
  end if
  deallocate(tmblarge%confdatarr, stat=istat)


  !Write the linear wavefunctions to file if asked
  if(input%lin%plotBasisFunctions /= WF_FORMAT_NONE) then
    call writemywaves_linear(iproc,trim(input%dir_output) // 'minBasis',input%lin%plotBasisFunctions,tmb%Lzd,&
       tmb%orbs,KSwfn%orbs%norb,at,rxyz,tmb%psi,tmb%wfnmd%coeff,KSwfn%orbs%eval)
  end if

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

  !!!!!call communicate_basis_for_density(iproc, nproc, tmb%lzd, tmb%orbs, tmb%psi, tmb%comsr)
  !!!call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, tmb%orbs, tmb%psi, tmb%collcom_sr)
  !!!call calculate_density_kernel(iproc, nproc, .true., tmb%orbs%norb, KSwfn%orbs, tmb%orbs, &
  !!!     tmb%wfnmd%coeff, tmb%wfnmd%density_kernel)
  !!call sumrhoForLocalizedBasis2(iproc, nproc, tmb%lzd, &
  !!     tmb%orbs, tmb%comsr, tmb%wfnmd%density_kernel, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, &
  !!     denspot%rhov, at,denspot%dpbox%nscatterarr)
  !!allocate(density_kernel(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  !!call memocc(istat, density_kernel, 'density_kernel', subname)
  !!call uncompressMatrix(tmb%orbs%norb, tmblarge%mad, tmb%wfnmd%density_kernel_compr, density_kernel)
  call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
       tmb%orbs, tmblarge%mad, tmb%collcom_sr, tmb%wfnmd%density_kernel_compr, &
       KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
  !!iall=-product(shape(density_kernel))*kind(density_kernel)
  !!deallocate(density_kernel, stat=istat)
  !!call memocc(istat, iall, 'density_kernel', subname)

  !!if (iproc==0) then
  !!    do istat=1,size(tmb%wfnmd%density_kernel_compr)
  !!        write(*,'(a,i6,es20.10)') 'istat, tmb%wfnmd%density_kernel_compr(istat)', istat, tmb%wfnmd%density_kernel_compr(istat)
  !!    end do
  !!end if

  !call destroy_new_locregs(iproc, nproc, tmblarge)
  !call deallocate_auxiliary_basis_function(subname, tmblarge%psi, tmblarge%hpsi)

  !!call deallocateCommunicationbufferSumrho(tmb%comsr, subname)

  !!! allocating here instead of input_wf to save memory
  !!allocate(KSwfn%psi(max(KSwfn%orbs%npsidim_comp,KSwfn%orbs%npsidim_orbs)+ndebug),stat=istat)
  !!call memocc(istat,KSwfn%psi,'KSwfn%psi',subname)


  !!! Build global orbitals psi (the physical ones).
  !!if(nproc>1) then
  !!   allocate(KSwfn%psit(max(KSwfn%orbs%npsidim_orbs,KSwfn%orbs%npsidim_comp)), stat=istat)
  !!   call memocc(istat, KSwfn%psit, 'KSwfn%psit', subname)
  !!else
  !!   KSwfn%psit => KSwfn%psi
  !!end if
  !!call transformToGlobal(iproc, nproc, tmb%lzd, tmb%orbs, KSwfn%orbs, KSwfn%comms, input, &
  !!     tmb%wfnmd%coeff, tmb%psi, KSwfn%psi, KSwfn%psit)
  !!if(nproc>1) then
  !!   iall=-product(shape(KSwfn%psit))*kind(KSwfn%psit)
  !!   deallocate(KSwfn%psit, stat=istat)
  !!   call memocc(istat, iall, 'KSwfn%psit', subname)
  !!else
  !!   nullify(KSwfn%psit)
  !!end if


  ! Otherwise there are some problems... Check later.
  allocate(KSwfn%psi(1),stat=istat)
  call memocc(istat,KSwfn%psi,'KSwfn%psi',subname)
  nullify(KSwfn%psit)

  nullify(rho,pot)

  call deallocate_local_arrays()

  call timing(iproc,'WFN_OPT','PR')

  contains

    subroutine allocate_local_arrays()

      allocate(locrad(tmb%lzd%nlr), stat=istat)
      call memocc(istat, locrad, 'locrad', subname)

      ! Allocate the old charge density (used to calculate the variation in the charge density)
      allocate(rhopotold_out(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim)),stat=istat)
      call memocc(istat, rhopotold_out, 'rhopotold_out', subname)

      allocate(locrad_tmp(tmb%lzd%nlr), stat=istat)
      call memocc(istat, locrad_tmp, 'locrad_tmp', subname)

      allocate(ham_compr(tmblarge%mad%nvctr), stat=istat)
      call memocc(istat, ham_compr, 'ham_compr', subname)

      allocate(overlapmatrix_compr(tmblarge%mad%nvctr), stat=istat)
      call memocc(istat, overlapmatrix_compr, 'overlapmatrix_compr', subname)


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

      iall=-product(shape(ham_compr))*kind(ham_compr)
      deallocate(ham_compr, stat=istat)
      call memocc(istat, iall, 'ham_compr', subname)

      iall=-product(shape(overlapmatrix_compr))*kind(overlapmatrix_compr)
      deallocate(overlapmatrix_compr, stat=istat)
      call memocc(istat, iall, 'overlapmatrix_compr', subname)

    end subroutine deallocate_local_arrays


    subroutine check_inputguess()
      real(8) :: dnrm2
      if (input%inputPsiId==101) then           !should we put 102 also?

          if (input%lin%pulay_correction) then
             ! Check the input guess by calculation the Pulay forces.

             call to_zero(tmblarge%orbs%npsidim_orbs,tmblarge%psi(1))
             call small_to_large_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmb%orbs, tmblarge%orbs, tmb%psi, tmblarge%psi)

             ! add get_coeff here
             ! - need some restructuring/reordering though, or addition of lots of extra initializations?!

             ! Calculate Pulay correction to the forces
             call pulay_correction(iproc, nproc, KSwfn%orbs, at, rxyz, nlpspd, proj, input%SIC, denspot, GPU, tmb, &
                  tmblarge, fpulay)
             fnrm_pulay=dnrm2(3*at%nat, fpulay, 1)/sqrt(dble(at%nat))

             if (iproc==0) write(*,*) 'fnrm_pulay',fnrm_pulay

             if (fnrm_pulay>1.d-1) then !1.d3 1.d-1
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
                call inputguessConfinement(iproc, nproc, at, input, &
                     KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
                     rxyz, nlpspd, proj, GPU, KSwfn%orbs, tmb, tmblarge, denspot, rhopotold, energs)
                     energs%eexctX=0.0_gp
                ! Give tmblarge%mad since this is the correct matrix description
                call orthonormalizeLocalized(iproc, nproc, 0, tmb%orbs, tmb%lzd, tmblarge%mad, tmb%collcom, &
                     tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)
             else if (fnrm_pulay>1.d-2) then ! 1.d2 1.d-2
                if (iproc==0) write(*,'(1x,a)') 'The pulay forces are rather large, so start with low accuracy.'
                nit_lowaccuracy=input%lin%nit_lowaccuracy
                nit_highaccuracy=input%lin%nit_highaccuracy
             else if (fnrm_pulay>1.d-10) then
                if (iproc==0) write(*,'(1x,a)') &
                     'The pulay forces are fairly large, so reoptimising basis with high accuracy only.'
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
              call to_zero(3*at%nat, fpulay(1,1))
              nit_lowaccuracy=input%lin%nit_lowaccuracy!0
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
  
      if(lowaccur_converged) then
          !cur_it_highaccuracy=cur_it_highaccuracy+1
          if(cur_it_highaccuracy==nit_highaccuracy) then
              exit_outer_loop=.true.
          else if (pnrm_out<input%lin%highaccuracy_conv_crit) then
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
              write(*,'(3x,a)') 'coefficients obtained by direct minimization.'
          else
              write(*,'(3x,a)') 'coefficients obtained by diagonalization.'
          end if
          if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE .or. input%lin%scf_mode==LINEAR_FOE) then
              write(*,'(3x,a,3x,i0,2x,es13.7,es27.17,es14.4)') 'it, Delta DENS, energy, energyDiff', &
                   it_SCC, pnrm, energy, energyDiff
          else if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
              write(*,'(3x,a,3x,i0,2x,es13.7,es27.17,es14.4)') 'it, Delta POT, energy, energyDiff', &
                   it_SCC, pnrm, energy, energyDiff
          else if(input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
              write(*,'(3x,a,3x,i0,2x,es13.7,es27.17,es14.4)') 'it, fnrm coeff, energy, energyDiff', &
                   it_SCC, pnrm, energy, energyDiff
          end if
          write(*,'(1x,a)') repeat('+',92 + int(log(real(it_SCC))/log(10.)))
      end if

    end subroutine printSummary

    !> Print a short summary of some values calculated during the last iteration in the self
    !! consistency cycle.
    subroutine print_info()
      implicit none

      real(8) :: energyDiff

      energyDiff = energy - energyoldout

      ! Print out values related to two iterations of the outer loop.
      if(iproc==0) then

          !Before convergence
          write(*,'(3x,a,7es18.10)') 'ebs, ehart, eexcu, vexcu, eexctX, eion, edisp', &
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
              write(*,'(5x,a,3x,i0,es12.2,es27.17)') 'FINAL values: it, fnrm coeff, energy', itout, pnrm, energy
          end if
          write(*,'(3x,a,es14.6)') '> energy difference to last iteration:', energyDiff
          write(*,'(1x,a)') repeat('#',92 + int(log(real(itout))/log(10.)))
      end if

    end subroutine print_info

end subroutine linearScaling



subroutine transformToGlobal(iproc,nproc,lzd,lorbs,orbs,comms,input,coeff,lphi,psi,psit)
use module_base
use module_types
use module_interfaces, exceptThisOne => transformToGlobal
implicit none

! Calling arguments
integer,intent(in) :: iproc, nproc
type(local_zone_descriptors),intent(in) :: lzd
type(orbitals_data),intent(in) :: lorbs, orbs
type(communications_arrays) :: comms
type(input_variables),intent(in) :: input
real(8),dimension(lorbs%norb,orbs%norb),intent(in) :: coeff
real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(inout) :: lphi
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),target,intent(out) :: psi
real(8),dimension(:),pointer,intent(inout) :: psit

! Local variables
integer :: ind1, ind2, istat, iall, iorb, ilr, ldim, gdim, nvctrp
real(8),dimension(:),pointer :: phiWork
real(8),dimension(:),allocatable :: phi
character(len=*),parameter :: subname='transformToGlobal'
type(orbitals_data) :: gorbs
type(communications_arrays) :: gcomms

  call nullify_orbitals_data(gorbs)
  call copy_orbitals_data(lorbs, gorbs, subname)
  call orbitals_communicators(iproc,nproc,lzd%glr,gorbs,gcomms)

  allocate(phi(max(gorbs%npsidim_orbs,gorbs%npsidim_comp)+ndebug), stat=istat)
  call memocc(istat, phi, 'phi', subname)
  allocate(phiWork(max(size(phi),size(psi))), stat=istat)
  call memocc(istat, phiWork, 'phiWork', subname)

  ind1=1
  ind2=1
  if (max(gorbs%npsidim_orbs,gorbs%npsidim_comp) > 0) &
       call to_zero(max(gorbs%npsidim_orbs,gorbs%npsidim_comp),phi(1))

  do iorb=1,lorbs%norbp
      ilr = lorbs%inWhichLocreg(lorbs%isorb+iorb)
      ldim=lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
      gdim=lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f

      call Lpsi_to_global2(iproc,nproc,ldim,gdim,lorbs%norb,lorbs%nspinor,input%nspin,lzd%Glr,&
           lzd%Llr(ilr),lphi(ind2),phi(ind1))
      ind1=ind1+lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f
      ind2=ind2+lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f

  end do
  call transpose_v(iproc, nproc, gorbs, lzd%Glr%wfd, gcomms, phi, work=phiWork)

  if(iproc==0) then
      write(*,'(1x,a)', advance='no') '------------------------------------- Building linear combinations... '
  end if
  ! Build the extended orbital psi as a linear combination of localized basis functions phi. for real O(N)
  ! this has to replaced, but at the moment it is still needed.
  !call buildWavefunctionModified(iproc, nproc, orbs, gorbs, comms, gcomms, phi, psi, coeff)

  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  call dgemm('n', 'n', nvctrp, orbs%norb, lorbs%norb, 1.d0, phi(1), nvctrp, coeff(1,1), &
       lorbs%norb, 0.d0, psi(1), nvctrp)

  ! not used in linearscaling
  !if(nproc>1) then
  !    call dcopy(orbs%npsidim_comp, psi, 1, psit, 1)
  !else
  !    psit => psi
  !end if

  call untranspose_v(iproc, nproc, orbs, lzd%Glr%wfd, comms, psi, work=phiWork)

  if(iproc==0) write(*,'(a)') 'done.'


  iall=-product(shape(phi))*kind(phi)
  deallocate(phi, stat=istat)
  call memocc(istat, iall, 'phi', subname)
  iall=-product(shape(phiWork))*kind(phiWork)
  deallocate(phiWork, stat=istat)
  call memocc(istat, iall, 'phiWork', subname)

  call deallocate_orbitals_data(gorbs, subname)
  call deallocate_communications_arrays(gcomms, subname)


end subroutine transformToGlobal



subroutine set_optimization_variables(input, at, lorbs, nlr, onwhichatom, confdatarr, &
     convCritMix, lowaccur_converged, nit_scc, mix_hist, alpha_mix, locrad, target_function, nit_basis)
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
  real(kind=8), intent(out) :: convCritMix, alpha_mix
  logical, intent(in) :: lowaccur_converged
  integer, intent(out) :: nit_scc, mix_hist
  real(kind=8), dimension(nlr), intent(out) :: locrad
  integer, intent(out) :: target_function, nit_basis

  ! Local variables
  integer :: iorb, ilr, iiat

  if(lowaccur_converged) then
      do iorb=1,lorbs%norbp
          iiat=onwhichatom(lorbs%isorb+iorb)
          confdatarr(iorb)%prefac=input%lin%potentialPrefac_highaccuracy(at%iatype(iiat))
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
  else
      do iorb=1,lorbs%norbp
          iiat=onwhichatom(lorbs%isorb+iorb)
          confdatarr(iorb)%prefac=input%lin%potentialPrefac_lowaccuracy(at%iatype(iiat))
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
  end if

end subroutine set_optimization_variables



subroutine adjust_locregs_and_confinement(iproc, nproc, hx, hy, hz, &
           at, input, tmb, denspot, ldiis, locreg_increased, lowaccur_converged, locrad)
  use module_base
  use module_types
  use module_interfaces, except_this_one => adjust_locregs_and_confinement
  implicit none
  
  ! Calling argument
  integer,intent(in) :: iproc, nproc
  real(8),intent(in) :: hx, hy, hz
  type(atoms_data),intent(in) :: at
  type(input_variables),intent(in) :: input
  type(DFT_wavefunction),intent(inout) :: tmb
  type(DFT_local_fields),intent(inout) :: denspot
  type(localizedDIISParameters),intent(inout) :: ldiis
  logical, intent(out) :: locreg_increased
  logical, intent(in) :: lowaccur_converged
  real(8), dimension(tmb%lzd%nlr), intent(inout) :: locrad

  ! Local variables
  integer :: ilr

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
      call redefine_locregs_quantities(iproc, nproc, hx, hy, hz, at, input, locrad, .true., tmb%lzd, tmb, denspot, ldiis)
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
     if (input%lin%scf_mode/=LINEAR_DIRECT_MINIMIZATION) then
        if(input%lin%mixHist_lowaccuracy==0 .and. input%lin%mixHist_highaccuracy>0) then
           call initializeMixrhopotDIIS(input%lin%mixHist_highaccuracy, denspot%dpbox%ndimpot, mixdiis)
        else if(input%lin%mixHist_lowaccuracy>0 .and. input%lin%mixHist_highaccuracy==0) then
           call deallocateMixrhopotDIIS(mixdiis)
        end if
     else
        ! check whether ldiis_coeff_hist arrays will need reallocating due to change in history length
        if (ldiis_coeff_hist /= input%lin%mixHist_highaccuracy) then
           ldiis_coeff_changed=.true.
        else
           ldiis_coeff_changed=.false.
        end if
        ldiis_coeff_hist=input%lin%mixHist_highaccuracy
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

subroutine pulay_correction(iproc, nproc, orbs, at, rxyz, nlpspd, proj, SIC, denspot, GPU, tmb, &
           tmblarge, fpulay)
  use module_base
  use module_types
  use module_interfaces, except_this_one => pulay_correction
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(atoms_data),intent(in) :: at
  real(kind=8),dimension(3,at%nat),intent(in) :: rxyz
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd
  real(wp),dimension(nlpspd%nprojel),intent(inout) :: proj
  type(SIC_data),intent(in) :: SIC
  type(DFT_local_fields), intent(inout) :: denspot
  type(GPU_pointers),intent(inout) :: GPU
  type(DFT_wavefunction),intent(in) :: tmb
  type(DFT_wavefunction),intent(inout) :: tmblarge
  real(kind=8),dimension(3,at%nat),intent(out) :: fpulay

  ! Local variables
  integer:: istat, iall, ierr, iialpha, jorb
  integer:: iorb, ii, iseg, isegstart, isegend
  integer:: jat, jdir, ibeta
  !!integer :: ialpha, iat, iiorb
  real(kind=8) :: kernel, ekernel
  real(kind=8),dimension(:),allocatable :: lhphilarge, psit_c, psit_f, hpsit_c, hpsit_f, lpsit_c, lpsit_f
  real(kind=8),dimension(:,:),allocatable :: matrix_compr, dovrlp_compr
  type(energy_terms) :: energs
  type(confpot_data),dimension(:),allocatable :: confdatarrtmp
  character(len=*),parameter :: subname='pulay_correction'

  ! Begin by updating the Hpsi
  call local_potential_dimensions(tmblarge%lzd,tmblarge%orbs,denspot%dpbox%ngatherarr(0,1))

  allocate(lhphilarge(tmblarge%orbs%npsidim_orbs), stat=istat)
  call memocc(istat, lhphilarge, 'lhphilarge', subname)
  call to_zero(tmblarge%orbs%npsidim_orbs,lhphilarge(1))

  !!call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
  !!     tmblarge%comgp%nrecvbuf, tmblarge%comgp%recvbuf, tmblarge%comgp, tmblarge%lzd)
  call start_onesided_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
       tmblarge%comgp%nrecvbuf, tmblarge%comgp%recvbuf, tmblarge%comgp, tmblarge%lzd)

  allocate(confdatarrtmp(tmblarge%orbs%norbp))
  call default_confinement_data(confdatarrtmp,tmblarge%orbs%norbp)

  allocate(tmblarge%lzd%doHamAppl(tmblarge%lzd%nlr), stat=istat)
  call memocc(istat, tmblarge%lzd%doHamAppl, 'tmblarge%lzd%doHamAppl', subname)
  tmblarge%lzd%doHamAppl=.true.

  call NonLocalHamiltonianApplication(iproc,at,tmblarge%orbs,rxyz,&
       proj,tmblarge%lzd,nlpspd,tmblarge%psi,lhphilarge,energs%eproj)

  ! only kinetic because waiting for communications
  call LocalHamiltonianApplication(iproc,nproc,at,tmblarge%orbs,&
       tmblarge%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmblarge%psi,lhphilarge,&
       energs,SIC,GPU,3,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmblarge%comgp)
  call full_local_potential(iproc,nproc,tmblarge%orbs,tmblarge%lzd,2,denspot%dpbox,denspot%rhov,denspot%pot_work, &
       tmblarge%comgp)
  ! only potential
  call LocalHamiltonianApplication(iproc,nproc,at,tmblarge%orbs,&
       tmblarge%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmblarge%psi,lhphilarge,&
       energs,SIC,GPU,2,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmblarge%comgp)

  call timing(iproc,'glsynchham1','ON') !lr408t
  call SynchronizeHamiltonianApplication(nproc,tmblarge%orbs,tmblarge%lzd,GPU,lhphilarge,&
       energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
  call timing(iproc,'glsynchham1','OF') !lr408t
  deallocate(confdatarrtmp)
  

  ! Now transpose the psi and hpsi
  allocate(lpsit_c(tmblarge%collcom%ndimind_c))
  call memocc(istat, lpsit_c, 'lpsit_c', subname)
  allocate(lpsit_f(7*tmblarge%collcom%ndimind_f))
  call memocc(istat, lpsit_f, 'lpsit_f', subname)
  allocate(hpsit_c(tmblarge%collcom%ndimind_c))
  call memocc(istat, hpsit_c, 'hpsit_c', subname)
  allocate(hpsit_f(7*tmblarge%collcom%ndimind_f))
  call memocc(istat, hpsit_f, 'hpsit_f', subname)
  allocate(psit_c(tmblarge%collcom%ndimind_c))
  call memocc(istat, psit_c, 'psit_c', subname)
  allocate(psit_f(7*tmblarge%collcom%ndimind_f))
  call memocc(istat, psit_f, 'psit_f', subname)

  call transpose_localized(iproc, nproc, tmblarge%orbs,  tmblarge%collcom, &
       tmblarge%psi, lpsit_c, lpsit_f, tmblarge%lzd)

  call transpose_localized(iproc, nproc, tmblarge%orbs,  tmblarge%collcom, &
       lhphilarge, hpsit_c, hpsit_f, tmblarge%lzd)

  !now build the derivative and related matrices <dPhi_a | H | Phi_b> and <dPhi_a | Phi_b>
  allocate(matrix_compr(tmblarge%mad%nvctr,3), stat=istat)
  call memocc(istat, matrix_compr, 'matrix_compr', subname)
  allocate(dovrlp_compr(tmblarge%mad%nvctr,3), stat=istat)
  call memocc(istat, dovrlp_compr, 'dovrlp_compr', subname)
  jdir=1
  do jdir = 1, 3
     call get_derivative(jdir, tmblarge%orbs%npsidim_orbs, tmblarge%lzd%hgrids(1), tmblarge%orbs, &
          tmblarge%lzd, tmblarge%psi, lhphilarge)

     call transpose_localized(iproc, nproc, tmblarge%orbs,  tmblarge%collcom, &
          lhphilarge, psit_c, psit_f, tmblarge%lzd)

     call calculate_overlap_transposed(iproc, nproc, tmblarge%orbs, tmblarge%mad, tmblarge%collcom,&
          psit_c, lpsit_c, psit_f, lpsit_f, dovrlp_compr(1,jdir))

     call calculate_overlap_transposed(iproc, nproc, tmblarge%orbs, tmblarge%mad, tmblarge%collcom,&
          psit_c, hpsit_c, psit_f, hpsit_f, matrix_compr(1,jdir))
  end do


  !DEBUG
  !!print *,'iproc,tmblarge%orbs%norbp',iproc,tmblarge%orbs%norbp
  !!if(iproc==0)then
  !!do iorb = 1, tmblarge%orbs%norb
  !!   do iiorb=1,tmblarge%orbs%norb
  !!      !print *,'Hamiltonian of derivative: ',iorb, iiorb, (matrix(iorb,iiorb,jdir),jdir=1,3)
  !!      print *,'Overlap of derivative: ',iorb, iiorb, (dovrlp(iorb,iiorb,jdir),jdir=1,3)
  !!   end do
  !!end do
  !!end if
  !!!Check if derivatives are orthogonal to functions
  !!if(iproc==0)then
  !!  do iorb = 1, tmbder%orbs%norb
  !!     !print *,'overlap of derivative: ',iorb, (dovrlp(iorb,iiorb),iiorb=1,tmblarge%orbs%norb)
  !!     do iiorb=1,tmbder%orbs%norb
  !!         write(*,*) iorb, iiorb, dovrlp(iorb,iiorb)
  !!     end do
  !!  end do
  !!end if
  !END DEBUG

   call to_zero(3*at%nat, fpulay(1,1))
   do jdir=1,3
     !do ialpha=1,tmblarge%orbs%norb
     if (tmblarge%orbs%norbp>0) then
         isegstart=tmblarge%mad%istsegline(tmblarge%orbs%isorb_par(iproc)+1)
         if (tmblarge%orbs%isorb+tmblarge%orbs%norbp<tmblarge%orbs%norb) then
             isegend=tmblarge%mad%istsegline(tmblarge%orbs%isorb_par(iproc+1)+1)-1
         else
             isegend=tmblarge%mad%nseg
         end if
         do iseg=isegstart,isegend
              ii=tmblarge%mad%keyv(iseg)-1
              do jorb=tmblarge%mad%keyg(1,iseg),tmblarge%mad%keyg(2,iseg)
                  ii=ii+1
                  iialpha = (jorb-1)/tmblarge%orbs%norb + 1
                  ibeta = jorb - (iialpha-1)*tmblarge%orbs%norb
                  jat=tmblarge%orbs%onwhichatom(iialpha)
                  kernel = 0.d0
                  ekernel= 0.d0
                  do iorb=1,orbs%norb
                      kernel  = kernel+orbs%occup(iorb)*tmb%wfnmd%coeff(iialpha,iorb)*tmb%wfnmd%coeff(ibeta,iorb)
                      ekernel = ekernel+orbs%eval(iorb)*orbs%occup(iorb)*tmb%wfnmd%coeff(iialpha,iorb)*tmb%wfnmd%coeff(ibeta,iorb) 
                  end do
                  fpulay(jdir,jat)=fpulay(jdir,jat)+&
                         2.0_gp*(kernel*matrix_compr(ii,jdir)-ekernel*dovrlp_compr(ii,jdir))
              end do
         end do
     end if
     !!do ialpha=1,tmblarge%orbs%norbp
     !!  iialpha=tmblarge%orbs%isorb+ialpha
     !!  jat=tmblarge%orbs%onwhichatom(iialpha)
     !!  do ibeta=1,tmblarge%orbs%norb
     !!     kernel = 0.d0
     !!     ekernel= 0.d0
     !!     do iorb=1,orbs%norb
     !!       kernel  = kernel+orbs%occup(iorb)*tmb%wfnmd%coeff(iialpha,iorb)*tmb%wfnmd%coeff(ibeta,iorb)
     !!       ekernel = ekernel+orbs%eval(iorb)*orbs%occup(iorb)*tmb%wfnmd%coeff(iialpha,iorb)*tmb%wfnmd%coeff(ibeta,iorb) 
     !!     end do
     !!     !do iat=1,at%nat
     !!     !if(jat == iat ) then
     !!     !!fpulay(jdir,jat)=fpulay(jdir,jat)+&
     !!     !!       2.0_gp*(kernel*matrix(ibeta,iialpha,jdir)-ekernel*dovrlp(ibeta,iialpha,jdir))
     !!     fpulay(jdir,jat)=fpulay(jdir,jat)+&
     !!            2.0_gp*(kernel*matrix_compr(ind,jdir)-ekernel*dovrlp_compr(ind,jdir))
     !!     !else
     !!     !fpulay(jdir,iat)=fpulay(jdir,iat)-&
     !!     !       2.0_gp/at%nat*(kernel*matrix(ibeta,ialpha,jdir)-ekernel*dovrlp(ibeta,ialpha,jdir))
     !!     !end if
     !!     !end do
     !!  end do
     !!end do
   end do 

   call mpiallred(fpulay(1,1), 3*at%nat, mpi_sum, bigdft_mpi%mpi_comm, ierr)

  if(iproc==0) then
       do jat=1,at%nat
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

  iall=-product(shape(tmblarge%lzd%doHamAppl))*kind(tmblarge%lzd%doHamAppl)
  deallocate(tmblarge%lzd%doHamAppl, stat=istat)
  call memocc(istat, iall, 'tmblarge%lzd%doHamAppl', subname)

  iall=-product(shape(denspot%pot_work))*kind(denspot%pot_work)
  deallocate(denspot%pot_work, stat=istat)
  call memocc(istat, iall, 'denspot%pot_work', subname)

  iall=-product(shape(matrix_compr))*kind(matrix_compr)
  deallocate(matrix_compr, stat=istat)
  call memocc(istat, iall, 'matrix_compr', subname)

  iall=-product(shape(dovrlp_compr))*kind(dovrlp_compr)
  deallocate(dovrlp_compr, stat=istat)
  call memocc(istat, iall, 'dovrlp_compr', subname)

  if(iproc==0) write(*,'(1x,a)') 'done.'

end subroutine pulay_correction



subroutine derivative_coeffs_from_standard_coeffs(orbs, tmb, tmbder)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(orbitals_data),intent(in) :: orbs
  type(DFT_wavefunction),intent(in) :: tmb
  type(DFT_wavefunction),intent(out) :: tmbder

  ! Local variables
  integer :: iorb, jorb, jjorb

  call to_zero(tmbder%orbs%norb*orbs%norb, tmbder%wfnmd%coeff(1,1))
  do iorb=1,orbs%norb
      jjorb=0
      do jorb=1,tmbder%orbs%norb,4
          jjorb=jjorb+1
          tmbder%wfnmd%coeff(jorb,iorb)=tmb%wfnmd%coeff(jjorb,iorb)
      end do
  end do

end subroutine derivative_coeffs_from_standard_coeffs



subroutine derivatives_with_orthoconstraint(iproc, nproc, tmb, tmbder)
  use module_base
  use module_types
  use module_interfaces, except_this_one => derivatives_with_orthoconstraint
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(DFT_wavefunction),intent(in) :: tmb
  type(DFT_wavefunction),intent(inout) :: tmbder

  ! Local variables
  integer :: i0, j0, ii, jj, ipt, i, iiorb, jjorb, istat, iall, j
  real(8),dimension(:),allocatable :: psit_c, psit_f, psidert_c, psidert_f
  real(8),dimension(:,:),allocatable :: matrix
  character(len=*),parameter :: subname='derivatives_with_orthoconstraint'


  write(*,*) 'WARNING: in derivatives_with_orthoconstraint'

  allocate(psit_c(tmb%collcom%ndimind_c), stat=istat)
  call memocc(istat, psit_c, 'psit_c', subname)
  allocate(psit_f(7*tmb%collcom%ndimind_f), stat=istat)
  call memocc(istat, psit_f, 'psit_f', subname)

  allocate(psidert_c(tmbder%collcom%ndimind_c), stat=istat)
  call memocc(istat, psidert_c, 'psidert_c', subname)
  allocate(psidert_f(7*tmbder%collcom%ndimind_f), stat=istat)
  call memocc(istat, psidert_f, 'psidert_f', subname)

  do istat=1,size(tmbder%psi)
      write(200,*) istat, tmbder%psi(istat)
  end do

  ! Transpose the support functions
  call transpose_localized(iproc, nproc, tmb%orbs,  tmb%collcom, &
       tmb%psi, psit_c, psit_f, tmb%lzd)
  do istat=1,size(psit_c)
      write(201,*) psit_c(istat)
  end do
  do istat=1,size(psit_f)
      write(201,*) psit_f(istat)
  end do

  ! Transpose the derivatives
  call transpose_localized(iproc, nproc, tmbder%orbs,  tmbder%collcom, &
       tmbder%psi, psidert_c, psidert_f, tmb%lzd)


  allocate(matrix(tmbder%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, matrix, 'matrix', subname)

  ! Calculate the matrix <dphi_i|phi_j>
  call calculate_pulay_overlap(iproc, nproc, tmbder%orbs, tmb%orbs, tmbder%collcom, &
       tmb%collcom, psidert_c, psit_c, psidert_f, psit_f, matrix)
  do i=1,tmb%orbs%norb
      do j=1,tmbder%orbs%norb
          if(iproc==0) write(400,*) i,j,matrix(j,i)
      end do
  end do


  i0=0
  j0=0
  do ipt=1,tmb%collcom%nptsp_c 
      ii=tmb%collcom%norb_per_gridpoint_c(ipt) 
      jj=tmbder%collcom%norb_per_gridpoint_c(ipt) 
      do i=1,jj
          jjorb=tmbder%collcom%indexrecvorbital_c(j0+i)
          do j=1,ii
              iiorb=tmb%collcom%indexrecvorbital_c(i0+j)
              write(333,'(3i8,3es15.5)') jjorb, iiorb, ceiling(dble(jjorb)/3.d0), &
                                         5d0*matrix(jjorb,iiorb), psidert_c(j0+i), psit_c(i0+j)
              psidert_c(j0+i)=psidert_c(j0+i)-.5d0*matrix(jjorb,iiorb)*psit_c(i0+j)
              if (iiorb==ceiling(dble(jjorb)/3.d0)) then
                  psidert_c(j0+i)=psidert_c(j0+i)-.5d0*matrix(iiorb,jjorb)*psit_c(i0+j)
              end if
          end do
      end do
      i0=i0+ii
      j0=j0+jj
  end do

  i0=0
  j0=0
  do ipt=1,tmb%collcom%nptsp_f 
      ii=tmb%collcom%norb_per_gridpoint_f(ipt) 
      jj=tmbder%collcom%norb_per_gridpoint_f(ipt) 
      do i=1,jj
          jjorb=tmbder%collcom%indexrecvorbital_f(j0+i)
          do j=1,ii
              iiorb=tmb%collcom%indexrecvorbital_f(i0+j)
              psidert_f(7*(j0+i)-6)=psidert_f(7*(j0+i)-6)-.5d0*matrix(jjorb,iiorb)*psit_f(7*(i0+j)-6)
              psidert_f(7*(j0+i)-5)=psidert_f(7*(j0+i)-5)-.5d0*matrix(jjorb,iiorb)*psit_f(7*(i0+j)-5)
              psidert_f(7*(j0+i)-4)=psidert_f(7*(j0+i)-4)-.5d0*matrix(jjorb,iiorb)*psit_f(7*(i0+j)-4)
              psidert_f(7*(j0+i)-3)=psidert_f(7*(j0+i)-3)-.5d0*matrix(jjorb,iiorb)*psit_f(7*(i0+j)-3)
              psidert_f(7*(j0+i)-2)=psidert_f(7*(j0+i)-2)-.5d0*matrix(jjorb,iiorb)*psit_f(7*(i0+j)-2)
              psidert_f(7*(j0+i)-1)=psidert_f(7*(j0+i)-1)-.5d0*matrix(jjorb,iiorb)*psit_f(7*(i0+j)-1)
              psidert_f(7*(j0+i)-0)=psidert_f(7*(j0+i)-0)-.5d0*matrix(jjorb,iiorb)*psit_f(7*(i0+j)-0)
              if (iiorb==ceiling(dble(jjorb)/3.d0)) then
                  psidert_f(7*(j0+i)-6)=psidert_f(7*(j0+i)-6)-.5d0*matrix(iiorb,jjorb)*psit_f(7*(i0+j)-6)
                  psidert_f(7*(j0+i)-5)=psidert_f(7*(j0+i)-5)-.5d0*matrix(iiorb,jjorb)*psit_f(7*(i0+j)-5)
                  psidert_f(7*(j0+i)-4)=psidert_f(7*(j0+i)-4)-.5d0*matrix(iiorb,jjorb)*psit_f(7*(i0+j)-4)
                  psidert_f(7*(j0+i)-3)=psidert_f(7*(j0+i)-3)-.5d0*matrix(iiorb,jjorb)*psit_f(7*(i0+j)-3)
                  psidert_f(7*(j0+i)-2)=psidert_f(7*(j0+i)-2)-.5d0*matrix(iiorb,jjorb)*psit_f(7*(i0+j)-2)
                  psidert_f(7*(j0+i)-1)=psidert_f(7*(j0+i)-1)-.5d0*matrix(iiorb,jjorb)*psit_f(7*(i0+j)-1)
                  psidert_f(7*(j0+i)-0)=psidert_f(7*(j0+i)-0)-.5d0*matrix(iiorb,jjorb)*psit_f(7*(i0+j)-0)
              end if
          end do
      end do
      i0=i0+ii
      j0=j0+jj
  end do

  !! TEST ONLY
  call calculate_pulay_overlap(iproc, nproc, tmbder%orbs, tmb%orbs, tmbder%collcom, &
       tmb%collcom, psidert_c, psit_c, psidert_f, psit_f, matrix)
  do i=1,tmb%orbs%norb
      do j=1,tmbder%orbs%norb
          if(iproc==0) write(450,*) i,j,matrix(j,i)
      end do
  end do

  !!do istat=1,size(tmbder%psi)
  !!    write(200+iproc,*) istat, tmbder%psi(istat)
  !!end do

  ! Untranpose the derivatives
  call untranspose_localized(iproc, nproc, tmbder%orbs, tmbder%collcom, psidert_c, psidert_f, tmbder%psi, tmb%lzd)
  !!do istat=1,size(tmbder%psi)
  !!    write(300+iproc,*) istat, tmbder%psi(istat)
  !!end do
  do istat=1,size(tmbder%psi)
      write(250,*) istat, tmbder%psi(istat)
  end do

  iall=-product(shape(matrix))*kind(matrix)
  deallocate(matrix,stat=istat)
  call memocc(istat,iall,'matrix',subname)

  iall=-product(shape(psit_c))*kind(psit_c)
  deallocate(psit_c,stat=istat)
  call memocc(istat,iall,'psit_c',subname)

  iall=-product(shape(psit_f))*kind(psit_f)
  deallocate(psit_f,stat=istat)
  call memocc(istat,iall,'psit_f',subname)

  iall=-product(shape(psidert_c))*kind(psidert_c)
  deallocate(psidert_c,stat=istat)
  call memocc(istat,iall,'psidert_c',subname)

  iall=-product(shape(psidert_f))*kind(psidert_f)
  deallocate(psidert_f,stat=istat)
  call memocc(istat,iall,'psidert_f',subname)

end subroutine derivatives_with_orthoconstraint
