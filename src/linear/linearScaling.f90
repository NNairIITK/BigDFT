subroutine linearScaling(iproc,nproc,KSwfn,tmb,at,input,&
           rxyz,fion,fdisp,denspot,rhopotold,nlpspd,proj,GPU,&
           energs,scpot,energy,fpulay,infocode)

use module_base
use module_types
use module_interfaces, exceptThisOne => linearScaling
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(atoms_data),intent(inout):: at
type(input_variables),intent(in):: input
real(8),dimension(3,at%nat),intent(inout):: rxyz
real(8),dimension(3,at%nat),intent(in):: fion, fdisp
real(8),dimension(3,at%nat),intent(out):: fpulay
type(DFT_local_fields), intent(inout) :: denspot
real(gp), dimension(:), intent(inout) :: rhopotold
type(nonlocal_psp_descriptors),intent(in):: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
type(GPU_pointers),intent(in out):: GPU
type(energy_terms),intent(inout) :: energs
logical,intent(in):: scpot
real(gp), dimension(:), pointer :: rho,pot
real(8),intent(out):: energy
type(DFT_wavefunction),intent(inout),target:: tmb
type(DFT_wavefunction),intent(inout),target:: KSwfn
integer,intent(out):: infocode

type(linear_scaling_control_variables):: lscv
real(8):: pnrm,trace,fnrm_tmb
integer:: infoCoeff,istat,iall,it_scc,ilr,itout,scf_mode,info_scf,nsatur
character(len=*),parameter:: subname='linearScaling'
real(8),dimension(:),allocatable:: rhopotold_out
real(8):: energyold, energyDiff, energyoldout, dnrm2, fnrm_pulay
type(mixrhopotDIISParameters):: mixdiis
type(localizedDIISParameters):: ldiis, ldiis_coeff
logical:: calculate_overlap_matrix, can_use
logical:: fix_support_functions, check_initialguess
integer:: nit_highaccur, itype, istart, nit_lowaccuracy, iorb, iiorb
real(8),dimension(:,:),allocatable:: overlapmatrix, ham
real(8),dimension(:),allocatable :: locrad_tmp, eval
type(DFT_wavefunction):: tmblarge
real(8),dimension(:),pointer:: lhphilarge, lhphilargeold, lphilargeold


  call timing(iproc,'linscalinit','ON') !lr408t

  call allocate_local_arrays()

  if(iproc==0) then
      write(*,'(1x,a)') repeat('*',84)
      write(*,'(1x,a)') '****************************** LINEAR SCALING VERSION ******************************'
  end if


  if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE) then
      call dcopy(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim),rhopotold(1),1,rhopotold_out(1),1)
  end if

  if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
      call dcopy(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim),denspot%rhov(1),1,rhopotold_out(1),1)
  end if

  ! Copy the current potential
  if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
       call dcopy(max(KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p,1) &
            *input%nspin, denspot%rhov(1), 1, rhopotOld(1), 1)
  end if

  ! Allocate the communications buffers needed for the communications of teh potential and
  ! post the messages. This will send to each process the part of the potential that this process
  ! needs for the application of the Hamlitonian to all orbitals on that process.
  call allocateCommunicationsBuffersPotential(tmb%comgp, subname)

  ! Initialize the DIIS mixing of the potential if required.
  if(input%lin%mixHist_lowaccuracy>0) then
      call initializeMixrhopotDIIS(input%lin%mixHist_lowaccuracy, denspot%dpbox%ndimpot, mixdiis)
  end if

  ! Flag that indicates that the basis functions shall be improved in the following.
  tmb%wfnmd%bs%update_phi=.true.
  pnrm=1.d100
  lscv%pnrm_out=1.d100
  energyold=0.d0
  energyoldout=0.d0
  tmb%wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE
  lscv%lowaccur_converged=.false.
  lscv%info_basis_functions=-1
  lscv%enlarge_locreg=.false.
  nit_highaccur=0
  nsatur=0
  fix_support_functions=.false.
  check_initialguess=.true.

  !!do iorb=1,tmb%orbs%norb
  !!    iiorb=tmb%orbs%isorb+iorb
  !!    ilr=tmb%orbs%inwhichlocreg(iiorb)
  !!    call plotGrid(iproc, nproc, tmb%orbs%norb, 1, 1, iiorb, tmb%lzd%llr(ilr), tmb%lzd%glr, at, rxyz, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3))
  !!end do


  ! Allocate the communication arrays for the calculation of the charge density.
  call allocateCommunicationbufferSumrho(iproc, tmb%comsr, subname)

  call vcopy(tmb%orbs%norb, tmb%orbs%eval(1), 1, eval(1), 1)
  call timing(iproc,'linscalinit','OF') !lr408t


  ! This is the main outer loop. Each iteration of this loop consists of a first loop in which the basis functions
  ! are optimized and a consecutive loop in which the density is mixed.
  call initialize_DIIS_coeff(3, ldiis_coeff)
  call allocate_DIIS_coeff(tmb, KSwfn%orbs, ldiis_coeff)

  if (tmb%restart_method == LINEAR_HIGHACCURACY) then
      nit_lowaccuracy=0
  else
      nit_lowaccuracy=input%lin%nit_lowaccuracy
  end if

  ! Add one iteration if no low accuracy is desired since we need then a first fake iteration.
  if (nit_lowaccuracy>0) then
      istart=1
  else if (nit_lowaccuracy==0) then
      istart=0
  end if


  !!call plot_density(iproc,nproc,'potential-start',at,rxyz,denspot%dpbox,1,denspot%rhov)

  infocode=0 !default value

  outerLoop: do itout=istart,nit_lowaccuracy+input%lin%nit_highaccuracy

      ! First to some initialization and determine the value of some control parameters.
      ! The basis functions shall be optimized
      tmb%wfnmd%bs%update_phi=.true.
      ! Convergence criterion for the self consistency loop
      lscv%self_consistent=input%lin%convCritMix
      ! Check whether the low accuracy part (i.e. with strong confining potential) has converged.
      call check_whether_lowaccuracy_converged(itout, nit_lowaccuracy, input%lin%lowaccuray_converged, lscv)
      ! Set all remaining variables that we need for the optimizations of the basis functions and the mixing.
      call set_optimization_variables(input, at, tmb%orbs, tmb%lzd%nlr, tmb%orbs%onwhichatom, &
           tmb%confdatarr, tmb%wfnmd, lscv)

      ! Do one fake iteration if no low accuracy is desired.
      if(nit_lowaccuracy==0 .and. itout==0) then
          lscv%lowaccur_converged=.false.
          lscv%nit_highaccuracy=0
      end if

      if(lscv%lowaccur_converged) nit_highaccur=nit_highaccur+1
      if(nit_highaccur==1) lscv%enlarge_locreg=.true.

      ! Adjust the confining potential if required.
      call adjust_locregs_and_confinement(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
           input, tmb, denspot, ldiis, lscv)

      ! Some special treatement if we are in the high accuracy part
      call adjust_DIIS_for_high_accuracy(input, tmb, denspot, ldiis, mixdiis, lscv)

      call initialize_DIIS_coeff(3, ldiis_coeff)


      ! Now all initializations are done...
      if(nit_highaccur==1) then
          !!call plot_density(iproc,nproc,'potential-afterlowaccur',at,rxyz,denspot%dpbox,1,denspot%rhov)
          call destroy_new_locregs(iproc, nproc, tmblarge)
          call deallocate_auxiliary_basis_function(subname, tmblarge%psi, lhphilarge, lhphilargeold, lphilargeold)
          if(tmblarge%can_use_transposed) then
              iall=-product(shape(tmblarge%psit_c))*kind(tmblarge%psit_c)
              deallocate(tmblarge%psit_c, stat=istat)
              call memocc(istat, iall, 'tmblarge%psit_c', subname)
              iall=-product(shape(tmblarge%psit_f))*kind(tmblarge%psit_f)
              deallocate(tmblarge%psit_f, stat=istat)
              call memocc(istat, iall, 'tmblarge%psit_f', subname)
          end if
          deallocate(tmblarge%confdatarr, stat=istat)
      end if


      ! 0 is the fake iteration for no low accuracy.
      if(itout==0 .or. itout==1 .or. nit_highaccur==1) then
          call create_large_tmbs(iproc, nproc, tmb, eval, denspot, input, at, rxyz, lscv%lowaccur_converged, &
               tmblarge, lhphilarge, lhphilargeold, lphilargeold)
               ! Set to zero the large wavefunction. Later only the inner part will be filled. It must be made sure
               ! that the outer part is not modified!
               if (tmblarge%orbs%npsidim_orbs > 0) call to_zero(tmblarge%orbs%npsidim_orbs,tmblarge%psi(1))
      end if


      if(itout==1 .or. itout==0) then
          ! Orthonormalize the TMBs
          ! just to be sure...
          tmb%can_use_transposed=.false.
          nullify(tmb%psit_c)
          nullify(tmb%psit_f)
          if(iproc==0) write(*,*) 'calling orthonormalizeLocalized (exact)'
          call orthonormalizeLocalized(iproc, nproc, 0, tmb%orthpar%nItOrtho, &
               tmb%orbs, tmb%op, tmb%comon, tmb%lzd, &
               tmb%mad, tmb%collcom, tmb%orthpar, tmb%wfnmd%bpo, tmb%psi, tmb%psit_c, tmb%psit_f, &
               tmb%can_use_transposed)
      end if


      if(itout>1 .or. (nit_lowaccuracy==0 .and. itout==1)) then
          call deallocateDIIS(ldiis)
      end if
      if (lscv%lowaccur_converged) then
          call initializeDIIS(input%lin%DIIS_hist_highaccur, tmb%lzd, tmb%orbs, tmb%orbs%norb, ldiis)
      else
          call initializeDIIS(input%lin%DIIS_hist_lowaccur, tmb%lzd, tmb%orbs, tmb%orbs%norb, ldiis)
      end if
      ldiis%DIISHistMin=0
      if (lscv%lowaccur_converged) then
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


      if (input%inputPsiId==101 .and. tmb%restart_method == LINEAR_HIGHACCURACY .and. check_initialguess) then
          ! Calculate Pulay correction to the forces
          call pulay_correction(iproc, nproc, input, KSwfn%orbs, at, rxyz, nlpspd, proj, input%SIC, denspot, GPU, tmb, &
               tmblarge, fpulay)
          fnrm_pulay=dnrm2(3*at%nat, fpulay, 1)/sqrt(dble(at%nat))
          if (iproc==0) write(*,*) 'fnrm_pulay',fnrm_pulay
          check_initialguess=.false.
          if (fnrm_pulay>1.d-1) then
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
              infocode=2
              exit outerLoop
          end if
      end if


      ! The self consistency cycle. Here we try to get a self consistent density/potential.
      ! In the first lscv%nit_scc_when_optimizing iteration, the basis functions are optimized, whereas in the remaining
      ! iteration the basis functions are fixed.
      do it_scc=1,lscv%nit_scc

         ! Do nothing if no low accuracy is desired.
         if (nit_lowaccuracy==0 .and. itout==0) then
             iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
             deallocate(tmb%psit_c, stat=istat)
             call memocc(istat, iall, 'tmb%psit_c', subname)
             iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
             deallocate(tmb%psit_f, stat=istat)
             call memocc(istat, iall, 'tmb%psit_f', subname)
             tmb%can_use_transposed=.false.
             !call deallocateDIIS(ldiis)
             cycle outerLoop
         end if

          scf_mode=input%lin%scf_mode

          ! Do not update the TMB if it_scc>lscv%nit_scc_when_optimizing
          if(it_scc>1) tmb%wfnmd%bs%update_phi=.false.

          ! Stop the optimization if it seems to saturate
          if(nsatur>=tmb%wfnmd%bs%nsatur_outer .or. fix_support_functions) then
              tmb%wfnmd%bs%update_phi=.false.
              if(it_scc==1) then
                  tmb%can_use_transposed=.false.
              end if
          end if

         ! Improve the trace minimizing orbitals.
          if(tmb%wfnmd%bs%update_phi) then

              call getLocalizedBasis(iproc,nproc,at,KSwfn%orbs,rxyz,denspot,GPU,trace,fnrm_tmb,lscv%info_basis_functions,&
                  nlpspd,proj,ldiis,input%SIC,tmb, tmblarge, lhphilarge, energs, ham)
              if(lscv%info_basis_functions>0) then
                  nsatur=nsatur+1
              end if
              tmb%can_use_transposed=.false. !since basis functions have changed...

              tmb%wfnmd%nphi=tmb%orbs%npsidim_orbs
              tmb%wfnmd%it_coeff_opt=0
              tmb%wfnmd%alpha_coeff=.2d0 !reset to default value

              if (input%inputPsiId==101 .and. lscv%info_basis_functions<0 .and. itout==1) then
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
                  !!if (associated(tmblarge%psit_c)) then
                  !!    iall=-product(shape(tmblarge%psit_c))*kind(tmblarge%psit_c)
                  !!    deallocate(tmblarge%psit_c, stat=istat)
                  !!    call memocc(istat, iall, 'tmblarge%psit_c', subname)
                  !!end if
                  !!if (associated(tmblarge%psit_f)) then
                  !!    iall=-product(shape(tmblarge%psit_f))*kind(tmblarge%psit_f)
                  !!    deallocate(tmblarge%psit_f, stat=istat)
                  !!    call memocc(istat, iall, 'tmblarge%psit_f', subname)
                  !!end if
                  infocode=2
                  exit outerLoop
              end if

          end if


          ! Only communicate the TMB for sumrho if required (i.e. only if the TMB were optimized).
          if(it_scc<=1) then
              tmb%wfnmd%bs%communicate_phi_for_lsumrho=.true.
              calculate_overlap_matrix=.true.
          else
              tmb%wfnmd%bs%communicate_phi_for_lsumrho=.false.
              calculate_overlap_matrix=.false.
          end if

          scf_mode=input%lin%scf_mode

          ! Calculate the coefficients
          ! Check whether we can use the Hamiltonian matrix from the TMB optimization
          can_use=.true.
          if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) then
              do itype=1,at%ntypes
                  if(input%lin%potentialPrefac_lowaccuracy(itype)/=0.d0) then
                      can_use=.false.
                      exit
                  end if
              end do
          else if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
              do itype=1,at%ntypes
                  if(input%lin%potentialPrefac_highaccuracy(itype)/=0.d0) then
                      can_use=.false.
                      exit
                  end if
              end do
          end if
          if(tmb%wfnmd%bs%update_phi .and. can_use .and. lscv%info_basis_functions>=0) then
              call get_coeff(iproc,nproc,scf_mode,tmb%lzd,KSwfn%orbs,at,rxyz,denspot,GPU,infoCoeff,energs%ebs,nlpspd,proj,&
                   input%SIC,tmb,pnrm,overlapmatrix,calculate_overlap_matrix,&
                   tmblarge, lhphilarge, ham=ham, ldiis_coeff=ldiis_coeff)
          else
              call get_coeff(iproc,nproc,scf_mode,tmb%lzd,KSwfn%orbs,at,rxyz,denspot,GPU,infoCoeff,energs%ebs,nlpspd,proj,&
                   input%SIC,tmb,pnrm,overlapmatrix,calculate_overlap_matrix,&
                   tmblarge, lhphilarge, ldiis_coeff=ldiis_coeff)
          end if


          ! Calculate the total energy.
          !!write(*,'(a,7es14.5)') 'energs%ebs,energs%eh,energs%exc,energs%evxc,energs%eexctX,energs%eion,energs%edisp',&
          !!            energs%ebs,energs%eh,energs%exc,energs%evxc,energs%eexctX,energs%eion,energs%edisp
          energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
          energyDiff=energy-energyold
          energyold=energy


          ! Calculate the charge density.
          call sumrhoForLocalizedBasis2(iproc, nproc, &
               tmb%lzd, input, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), tmb%orbs, tmb%comsr, &
               tmb%wfnmd%density_kernel, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, &
               denspot%rhov, at, denspot%dpbox%nscatterarr)
          !!if(.not.lscv%lowaccur_converged) then
          !!    call plot_density(iproc,nproc,'density-afterlowaccur',at,rxyz,denspot%dpbox,1,denspot%rhov)
          !!end if


          ! Mix the density.
          if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE) then
           lscv%compare_outer_loop = pnrm<lscv%self_consistent .or. it_scc==lscv%nit_scc
           call mix_main(iproc, nproc, lscv%mix_hist, lscv%compare_outer_loop, input, KSwfn%Lzd%Glr, lscv%alpha_mix, &
                denspot, mixdiis, rhopotold, rhopotold_out, pnrm, lscv%pnrm_out)
          end if

          ! Calculate the new potential.
          if(iproc==0) write(*,'(1x,a)') '---------------------------------------------------------------- Updating potential.'
          call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)

          ! Mix the potential
          if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
           lscv%compare_outer_loop = pnrm<lscv%self_consistent .or. it_scc==lscv%nit_scc
           call mix_main(iproc, nproc, lscv%mix_hist, lscv%compare_outer_loop, input, KSwfn%Lzd%Glr, lscv%alpha_mix, &
                denspot, mixdiis, rhopotold, rhopotold_out, pnrm, lscv%pnrm_out)
          end if

          ! Keep the support functions fixed if they converged and the density
          ! change is below the tolerance already in the very first iteration
          if(it_scc==1 .and. pnrm<lscv%self_consistent .and.  lscv%info_basis_functions>0) then
              fix_support_functions=.true.
          end if

          ! Write some informations.
          call printSummary(iproc, it_scc, lscv%info_basis_functions, &
               infoCoeff, pnrm, energy, energyDiff, input%lin%scf_mode)
          if(pnrm<lscv%self_consistent) then
              info_scf=it_scc
              exit
          else
              info_scf=-1
          end if

          if(nsatur<tmb%wfnmd%bs%nsatur_outer .and. it_scc<1) then
              ! Deallocate the transposed TMBs
              if(tmb%can_use_transposed) then
                  iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
                  deallocate(tmb%psit_c, stat=istat)
                  call memocc(istat, iall, 'tmb%psit_c', subname)
                  iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
                  deallocate(tmb%psit_f, stat=istat)
                  call memocc(istat, iall, 'tmb%psit_f', subname)
              end if
          end if

      end do


      if(tmb%can_use_transposed) then
          iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
          deallocate(tmb%psit_c, stat=istat)
          call memocc(istat, iall, 'tmb%psit_c', subname)
          iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
          deallocate(tmb%psit_f, stat=istat)
          call memocc(istat, iall, 'tmb%psit_f', subname)
      end if

      ! Print out values related to two iterations of the outer loop.
      if(iproc==0) then
          write(*,'(3x,a,7es18.10)') 'ebs, ehart, eexcu, vexcu, eexctX, eion, edisp', &
              energs%ebs, energs%eh, energs%exc, energs%evxc, energs%eexctX, energs%eion, energs%edisp
          if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE) then
             if (.not. lscv%lowaccur_converged) then
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                      'itoutL, Delta DENSOUT, energy, energyDiff', itout, lscv%pnrm_out, energy, &
                      energy-energyoldout
             else
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                      'itoutH, Delta DENSOUT, energy, energyDiff', itout, lscv%pnrm_out, energy, &
                      energy-energyoldout
             end if
          else if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
             if (.not. lscv%lowaccur_converged) then
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                      'itoutL, Delta POTOUT, energy energyDiff', itout, lscv%pnrm_out, energy, energy-energyoldout
             else
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                      'itoutH, Delta POTOUT, energy energyDiff', itout, lscv%pnrm_out, energy, energy-energyoldout
             end if
          else if(input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
             if (.not. lscv%lowaccur_converged) then
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                      'itoutL, fnrm coeff, energy energyDiff', itout, lscv%pnrm_out, energy, energy-energyoldout
             else
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                      'itoutH, fnrm coeff, energy energyDiff', itout, lscv%pnrm_out, energy, energy-energyoldout
             end if
          end if
      end if
      call print_info(iproc, itout, lscv%info_basis_functions, info_scf, input%lin%scf_mode, tmb%wfnmd%bs%target_function, &
           fnrm_tmb, pnrm, trace, energy, energy-energyoldout)

      energyoldout=energy


      call check_for_exit(input, lscv)
      if(lscv%exit_outer_loop) exit outerLoop

      if(lscv%pnrm_out<input%lin%support_functions_converged) then
          if(iproc==0) write(*,*) 'fix the support functions from now on'
          fix_support_functions=.true.
      end if


  end do outerLoop



  ! Deallocate eberything that is not needed any more.
  call deallocateDIIS(ldiis_coeff)
  call deallocateDIIS(ldiis)
  if(input%lin%mixHist_highaccuracy>0) then
      call deallocateMixrhopotDIIS(mixdiis)
  end if
  call wait_p2p_communication(iproc, nproc, tmb%comgp)
  call deallocateCommunicationsBuffersPotential(tmb%comgp, subname)

  ! Calculate Pulay correction to the forces
  call pulay_correction(iproc, nproc, input, KSwfn%orbs, at, rxyz, nlpspd, proj, input%SIC, denspot, GPU, tmb, &
       tmblarge, fpulay)


  call destroy_new_locregs(iproc, nproc, tmblarge)
  call deallocate_auxiliary_basis_function(subname, tmblarge%psi, lhphilarge, lhphilargeold, lphilargeold)
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
       tmb%orbs,KSwfn%orbs%norb,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),at,rxyz,&
       tmb%psi,tmb%wfnmd%coeff)
   end if

  !!open(unit=230+iproc)
  !!    do istat=1,KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d
  !!        write(230+iproc,*) istat,denspot%rhov(istat)
  !!    end do
  !!close(unit=230+iproc)
 

  call communicate_basis_for_density(iproc, nproc, tmb%lzd, tmb%orbs, tmb%psi, tmb%comsr)
  call calculate_density_kernel(iproc, nproc, .true., tmb%wfnmd%ld_coeff, KSwfn%orbs, tmb%orbs, &
       tmb%wfnmd%coeff, tmb%wfnmd%density_kernel)
  call sumrhoForLocalizedBasis2(iproc, nproc, tmb%lzd, input, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
       tmb%orbs, tmb%comsr, tmb%wfnmd%density_kernel, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, &
       denspot%rhov, at,denspot%dpbox%nscatterarr)
  !!call plot_density(iproc,nproc,'density-end',at,rxyz,denspot%dpbox,1,denspot%rhov)
  !!open(unit=210+iproc)
  !!    do istat=1,KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d
  !!        write(210+iproc,*) istat,denspot%rhov(istat)
  !!    end do
  !!close(unit=210+iproc)

  call deallocateCommunicationbufferSumrho(tmb%comsr, subname)

  ! allocating here instead of input_wf to save memory
  allocate(KSwfn%psi(max(KSwfn%orbs%npsidim_comp,KSwfn%orbs%npsidim_orbs)+ndebug),stat=istat)
  call memocc(istat,KSwfn%psi,'KSwfn%psi',subname)

  ! Build global orbitals psi (the physical ones).
  if(nproc>1) then
     allocate(KSwfn%psit(max(KSwfn%orbs%npsidim_orbs,KSwfn%orbs%npsidim_comp)), stat=istat)
     call memocc(istat, KSwfn%psit, 'KSwfn%psit', subname)
  else
     KSwfn%psit => KSwfn%psi
  end if
  call transformToGlobal(iproc, nproc, tmb%lzd, tmb%orbs, KSwfn%orbs, KSwfn%comms, input, tmb%wfnmd%ld_coeff, &
       tmb%wfnmd%coeff, tmb%psi, KSwfn%psi, KSwfn%psit)
  if(nproc>1) then
     iall=-product(shape(KSwfn%psit))*kind(KSwfn%psit)
     deallocate(KSwfn%psit, stat=istat)
     call memocc(istat, iall, 'KSwfn%psit', subname)
  else
     nullify(KSwfn%psit)
  end if

  nullify(rho,pot)

  call deallocate_local_arrays()

  call timing(iproc,'WFN_OPT','PR')

  !!open(unit=200+iproc)
  !!    do istat=1,tmb%orbs%npsidim_orbs
  !!        write(200+iproc,*) istat,tmb%psi(istat)
  !!    end do
  !!close(unit=200+iproc)


  !!open(unit=220+iproc)
  !!    do istat=1,KSwfn%orbs%norb
  !!        do iall=1,tmb%orbs%norb
  !!            write(220+iproc,*) istat,iall,tmb%wfnmd%coeff(iall,istat)
  !!        end do
  !!    end do
  !!close(unit=220+iproc)


  contains

    subroutine allocate_local_arrays()

      allocate(lscv%locrad(tmb%lzd%nlr), stat=istat)
      call memocc(istat, lscv%locrad, 'lscv%locrad', subname)

      allocate(ham(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
      call memocc(istat, ham, 'ham', subname)

      allocate(eval(tmb%orbs%norb), stat=istat)
      call memocc(istat, eval, 'eval', subname)

      ! Allocate the old charge density (used to calculate the variation in the charge density)
      allocate(rhopotold_out(max(denspot%dpbox%ndimrhopot,denspot%dpbox%nrhodim)),stat=istat)
      call memocc(istat, rhopotold_out, 'rhopotold_out', subname)

      allocate(locrad_tmp(tmb%lzd%nlr), stat=istat)
      call memocc(istat, locrad_tmp, 'locrad_tmp', subname)

      allocate(overlapmatrix(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
      call memocc(istat, overlapmatrix, 'overlapmatrix', subname)


    end subroutine allocate_local_arrays


    subroutine deallocate_local_arrays()

      iall=-product(shape(lscv%locrad))*kind(lscv%locrad)
      deallocate(lscv%locrad, stat=istat)
      call memocc(istat, iall, 'lscv%locrad', subname)

      iall=-product(shape(eval))*kind(eval)
      deallocate(eval, stat=istat)
      call memocc(istat, iall, 'eval', subname)

      iall=-product(shape(ham))*kind(ham)
      deallocate(ham, stat=istat)
      call memocc(istat, iall, 'ham', subname)

      iall=-product(shape(overlapmatrix))*kind(overlapmatrix)
      deallocate(overlapmatrix, stat=istat)
      call memocc(istat, iall, 'overlapmatrix', subname)

      iall=-product(shape(locrad_tmp))*kind(locrad_tmp)
      deallocate(locrad_tmp, stat=istat)
      call memocc(istat, iall, 'locrad_tmp', subname)

      iall=-product(shape(rhopotold_out))*kind(rhopotold_out)
      deallocate(rhopotold_out, stat=istat)
      call memocc(istat, iall, 'rhopotold_out', subname)

    end subroutine deallocate_local_arrays

end subroutine linearScaling





subroutine printSummary(iproc, itSCC, infoBasisFunctions, infoCoeff, pnrm, energy, energyDiff, scf_mode)
use module_base
use module_types
!
! Purpose:
! ========
!   Print a short summary of some values calculated during the last iteration in the self
!   consistency cycle.
! 
! Calling arguments:
! ==================
!   Input arguments
!   ---------------
!
implicit none

! Calling arguments
integer,intent(in):: iproc, itSCC, infoBasisFunctions, infoCoeff, scf_mode
real(8),intent(in):: pnrm, energy, energyDiff

  if(iproc==0) then
      write(*,'(1x,a)') repeat('+',92 + int(log(real(itSCC))/log(10.)))
      write(*,'(1x,a,i0,a)') 'at iteration ', itSCC, ' of the density optimization:'
      !!if(infoCoeff<0) then
      !!    write(*,'(3x,a)') '- WARNING: coefficients not converged!'
      !!else if(infoCoeff>0) then
      !!    write(*,'(3x,a,i0,a)') '- coefficients converged in ', infoCoeff, ' iterations.'
      if(scf_mode==LINEAR_DIRECT_MINIMIZATION) then
          write(*,'(3x,a)') 'coefficients obtained by direct minimization.'
      else
          write(*,'(3x,a)') 'coefficients obtained by diagonalization.'
      end if
      !!end if
      !if(mixingMethod=='dens') then
      if(scf_mode==LINEAR_MIXDENS_SIMPLE) then
          write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)') 'it, Delta DENS, energy, energyDiff', itSCC, pnrm, energy, energyDiff
      !else if(mixingMethod=='pot') then
      else if(scf_mode==LINEAR_MIXPOT_SIMPLE) then
          write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)') 'it, Delta POT, energy, energyDiff', itSCC, pnrm, energy, energyDiff
      else if(scf_mode==LINEAR_DIRECT_MINIMIZATION) then
          write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)') 'it, fnrm coeff, energy, energyDiff', itSCC, pnrm, energy, energyDiff
      end if
      write(*,'(1x,a)') repeat('+',92 + int(log(real(itSCC))/log(10.)))
  end if

end subroutine printSummary



subroutine print_info(iproc, itout, info_tmb, info_coeff, scf_mode, target_function, &
           fnrm_tmb, pnrm, value_tmb, energy, energyDiff)
use module_base
use module_types
!
! Purpose:
! ========
!   Print a short summary of some values calculated during the last iteration in the self
!   consistency cycle.
! 
! Calling arguments:
! ==================
!   Input arguments
!   ---------------
!
implicit none

! Calling arguments
integer,intent(in):: iproc, itout, info_tmb, info_coeff, scf_mode, target_function
real(8),intent(in):: fnrm_tmb, pnrm, value_tmb, energy, energyDiff

  if(iproc==0) then
      write(*,'(1x,a)') repeat('#',92 + int(log(real(itout))/log(10.)))
      write(*,'(1x,a,i0,a)') 'at iteration ', itout, ' of the outer loop:'
      write(*,'(3x,a)') '> basis functions optimization:'
      if(target_function==TARGET_FUNCTION_IS_TRACE) then
          write(*,'(5x,a)') '- target function is trace'
      else if(target_function==TARGET_FUNCTION_IS_ENERGY) then
          write(*,'(5x,a)') '- target function is energy'
      end if
      if(info_tmb<=0) then
          write(*,'(5x,a)') '- WARNING: basis functions not converged!'
      else
          write(*,'(5x,a,i0,a)') '- basis functions converged in ', info_tmb, ' iterations.'
      end if
      write(*,'(5x,a,es15.6,2x,es10.2)') 'Final values: target function, fnrm', value_tmb, fnrm_tmb
      write(*,'(3x,a)') '> density optimization:'
      if(scf_mode==LINEAR_DIRECT_MINIMIZATION) then
          write(*,'(5x,a)') '- using direct minimization.'
      else
          write(*,'(5x,a)') '- using diagonalization / mixing.'
      end if
      if(info_coeff<0) then
          write(*,'(5x,a)') '- WARNING: density optimization not converged!'
      else
          write(*,'(5x,a,i0,a)') '- density optimization converged in ', info_coeff, ' iterations.'
      end if
      if(scf_mode==LINEAR_MIXDENS_SIMPLE) then
          write(*,'(5x,a,3x,i0,es12.2,es27.17)') 'FINAL values: it, Delta DENS, energy', itout, pnrm, energy
      else if(scf_mode==LINEAR_MIXPOT_SIMPLE) then
          write(*,'(5x,a,3x,i0,es12.2,es27.17)') 'FINAL values: it, Delta POT, energy', itout, pnrm, energy
      else if(scf_mode==LINEAR_DIRECT_MINIMIZATION) then
          write(*,'(5x,a,3x,i0,es12.2,es27.17)') 'FINAL values: it, fnrm coeff, energy', itout, pnrm, energy
      end if
      write(*,'(3x,a,es14.6)') '> energy difference to last iteration:', energyDiff
      write(*,'(1x,a)') repeat('#',92 + int(log(real(itout))/log(10.)))
  end if

end subroutine print_info



subroutine transformToGlobal(iproc,nproc,lzd,lorbs,orbs,comms,input,ld_coeff,coeff,lphi,psi,psit)
use module_base
use module_types
use module_interfaces, exceptThisOne => transformToGlobal
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, ld_coeff
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: lorbs, orbs
type(communications_arrays):: comms
type(input_variables),intent(in):: input
real(8),dimension(ld_coeff,orbs%norb),intent(in):: coeff
real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(inout):: lphi
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),target,intent(out):: psi
real(8),dimension(:),pointer,intent(inout):: psit

! Local variables
integer:: ind1, ind2, istat, iall, iorb, ilr, ldim, gdim, nvctrp
real(8),dimension(:),pointer:: phiWork
real(8),dimension(:),allocatable:: phi
character(len=*),parameter:: subname='transformToGlobal'
type(orbitals_data):: gorbs
type(communications_arrays):: gcomms

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



subroutine set_optimization_variables(input, at, lorbs, nlr, onwhichatom, confdatarr, wfnmd, lscv)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: nlr
  type(orbitals_data),intent(in):: lorbs
  type(input_variables),intent(in):: input
  type(atoms_data),intent(in):: at
  integer,dimension(lorbs%norb),intent(in):: onwhichatom
  type(confpot_data),dimension(lorbs%norbp),intent(inout):: confdatarr
  type(wfn_metadata),intent(inout):: wfnmd
  type(linear_scaling_control_variables),intent(inout):: lscv

  ! Local variables
  integer:: iorb, ilr, iiat

  if(lscv%lowaccur_converged) then

      do iorb=1,lorbs%norbp
          ilr=lorbs%inwhichlocreg(lorbs%isorb+iorb)
          iiat=onwhichatom(lorbs%isorb+iorb)
          confdatarr(iorb)%prefac=input%lin%potentialPrefac_highaccuracy(at%iatype(iiat))
      end do
      wfnmd%bs%target_function=TARGET_FUNCTION_IS_ENERGY
      wfnmd%bs%nit_basis_optimization=input%lin%nItBasis_highaccuracy
      wfnmd%bs%conv_crit=input%lin%convCrit_highaccuracy
      lscv%nit_scc=input%lin%nitSCCWhenFixed_highaccuracy
      !!lscv%nit_scc_when_optimizing=input%lin%nitSCCWhenOptimizing_highaccuracy
      lscv%mix_hist=input%lin%mixHist_highaccuracy
      do ilr=1,nlr
          lscv%locrad(ilr)=input%lin%locrad_highaccuracy(ilr)
      end do
      lscv%alpha_mix=input%lin%alpha_mix_highaccuracy
      !!if(wfnmd%bs%update_phi) then
      !!    lscv%alpha_mix=input%lin%alphaMixWhenOptimizing_highaccuracy
      !!else
      !!    lscv%alpha_mix=input%lin%alphaMixWhenFixed_highaccuracy
      !!end if
      !!if(.not.lscv%enlarge_locreg) lscv%enlarge_locreg=.true.


  else

      do iorb=1,lorbs%norbp
          ilr=lorbs%inwhichlocreg(lorbs%isorb+iorb)
          iiat=onwhichatom(lorbs%isorb+iorb)
          confdatarr(iorb)%prefac=input%lin%potentialPrefac_lowaccuracy(at%iatype(iiat))
      end do
      wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE
      wfnmd%bs%nit_basis_optimization=input%lin%nItBasis_lowaccuracy
      wfnmd%bs%conv_crit=input%lin%convCrit_lowaccuracy
      lscv%nit_scc=input%lin%nitSCCWhenFixed_lowaccuracy
      !!lscv%nit_scc_when_optimizing=input%lin%nitSCCWhenOptimizing_lowaccuracy
      lscv%mix_hist=input%lin%mixHist_lowaccuracy
      do ilr=1,nlr
          lscv%locrad(ilr)=input%lin%locrad_lowaccuracy(ilr)
      end do
      lscv%alpha_mix=input%lin%alpha_mix_lowaccuracy
      !!if(wfnmd%bs%update_phi) then
      !!    lscv%alpha_mix=input%lin%alphaMixWhenOptimizing_lowaccuracy
      !!else
      !!    lscv%alpha_mix=input%lin%alphaMixWhenFixed_lowaccuracy
      !!end if

  end if

end subroutine set_optimization_variables



subroutine adjust_locregs_and_confinement(iproc, nproc, hx, hy, hz, &
           input, tmb, denspot, ldiis, lscv)
  use module_base
  use module_types
  use module_interfaces, except_this_one => adjust_locregs_and_confinement
  implicit none
  
  ! Calling argument
  integer,intent(in):: iproc, nproc
  real(8),intent(in):: hx, hy, hz
  type(input_variables),intent(in):: input
  type(DFT_wavefunction),intent(inout):: tmb
  type(DFT_local_fields),intent(inout) :: denspot
  type(localizedDIISParameters),intent(inout):: ldiis
  type(linear_scaling_control_variables),intent(inout):: lscv

  ! Local variables
  integer :: ilr
  logical:: change, locreg_increased


  locreg_increased=.false.
  if(lscv%lowaccur_converged .and. lscv%enlarge_locreg) then
      change = .false.
      do ilr = 1, tmb%lzd%nlr
         if(input%lin%locrad_highaccuracy(ilr) /= input%lin%locrad_lowaccuracy(ilr)) then
             change = .true.
             exit
         end if
      end do
      if(change) then
         if(iproc==0) then
             write(*,'(1x,a)') 'Increasing the localization radius for the high accuracy part.'
         end if
         locreg_increased=.true.
      end if
      lscv%enlarge_locreg=.false. !flag to indicate that the locregs should not be increased any more in the following iterations
  end if
  if(locreg_increased) then
      call redefine_locregs_quantities(iproc, nproc, hx, hy, hz, lscv%locrad, .true., tmb%lzd, tmb, denspot, ldiis)
  end if

end subroutine adjust_locregs_and_confinement



subroutine adjust_DIIS_for_high_accuracy(input, tmb, denspot, ldiis, mixdiis, lscv)
  use module_base
  use module_types
  use module_interfaces, except_this_one => adjust_DIIS_for_high_accuracy
  implicit none
  
  ! Calling arguments
  type(input_variables),intent(in):: input
  type(DFT_wavefunction),intent(in):: tmb
  type(DFT_local_fields),intent(inout) :: denspot
  type(localizedDIISParameters),intent(inout):: ldiis
  type(mixrhopotDIISParameters),intent(inout):: mixdiis
  type(linear_scaling_control_variables),intent(inout):: lscv
  
  !!lscv%exit_outer_loop=.false.
  
  if(lscv%lowaccur_converged) then
      if(input%lin%mixHist_lowaccuracy==0 .and. input%lin%mixHist_highaccuracy>0) then
          call initializeMixrhopotDIIS(input%lin%mixHist_highaccuracy, denspot%dpbox%ndimpot, mixdiis)
      else if(input%lin%mixHist_lowaccuracy>0 .and. input%lin%mixHist_highaccuracy==0) then
          call deallocateMixrhopotDIIS(mixdiis)
      end if
  end if
  
end subroutine adjust_DIIS_for_high_accuracy


subroutine check_for_exit(input, lscv)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(input_variables),intent(in):: input
  type(linear_scaling_control_variables),intent(inout):: lscv

  lscv%exit_outer_loop=.false.
  
  if(lscv%lowaccur_converged) then
      lscv%nit_highaccuracy=lscv%nit_highaccuracy+1
      if(lscv%nit_highaccuracy==input%lin%nit_highaccuracy) then
          lscv%exit_outer_loop=.true.
      else if (lscv%pnrm_out<input%lin%highaccuracy_converged) then !lr408
          lscv%exit_outer_loop=.true.
      end if
  end if

end subroutine check_for_exit


!!function check_whether_derivatives_to_be_used(input, itout, lscv)
!!  use module_base
!!  use module_types
!!  implicit none
!!  
!!  ! Calling arguments
!!  type(input_variables),intent(in):: input
!!  integer,intent(in):: itout
!!  type(linear_scaling_control_variables),intent(in):: lscv
!!  logical:: check_whether_derivatives_to_be_used
!!
!!  ! Local variables
!!  logical:: withder
!!
!!  if(input%lin%mixedmode) then
!!      if( (.not.lscv%lowaccur_converged .and. &
!!           (itout==input%lin%nit_lowaccuracy+1 .or. lscv%pnrm_out<input%lin%lowaccuray_converged) ) &
!!          .or. lscv%lowaccur_converged ) then
!!          withder=.true.
!!      else
!!          withder=.false.
!!      end if
!!  else
!!      if(input%lin%useDerivativeBasisFunctions) then
!!          withder=.true.
!!      else
!!          withder=.false.
!!      end if
!!  end if
!!  check_whether_derivatives_to_be_used=withder
!!
!!end function check_whether_derivatives_to_be_used



subroutine check_whether_lowaccuracy_converged(itout, nit_lowaccuracy, lowaccuray_converged, lscv)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: itout
  integer,intent(in):: nit_lowaccuracy
  real(8),intent(in):: lowaccuray_converged
  type(linear_scaling_control_variables),intent(inout):: lscv
  
  if(.not.lscv%lowaccur_converged .and. &
     (itout>=nit_lowaccuracy+1 .or. lscv%pnrm_out<lowaccuray_converged)) then
      lscv%lowaccur_converged=.true.
      lscv%nit_highaccuracy=0
  end if 

end subroutine check_whether_lowaccuracy_converged

subroutine pulay_correction(iproc, nproc, input, orbs, at, rxyz, nlpspd, proj, SIC, denspot, GPU, tmb, &
           tmblarge, fpulay)
  use module_base
  use module_types
  use module_interfaces, except_this_one => pulay_correction
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(input_variables),intent(in):: input
  type(orbitals_data),intent(in):: orbs
  type(atoms_data),intent(in):: at
  real(kind=8),dimension(3,at%nat),intent(in):: rxyz
  type(nonlocal_psp_descriptors),intent(in):: nlpspd
  real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
  type(SIC_data),intent(in):: SIC
  type(DFT_local_fields), intent(inout) :: denspot
  type(GPU_pointers),intent(inout):: GPU
  type(DFT_wavefunction),intent(in):: tmb
  type(DFT_wavefunction),intent(inout):: tmblarge
  real(kind=8),dimension(3,at%nat),intent(out):: fpulay

  ! Local variables
  integer:: norb, norbu, norbd, istat, iall, nspin, tag, ierr
  integer:: iorb, iiorb, ilr, nlr, ityp
  integer:: jjorb, jat, jorbsmall, kkorb,  jdir, iat
  integer:: lorbsmall, ldir, lat, llorb
  type(DFT_wavefunction):: tmbder
  real(kind=8),dimension(:),allocatable:: lhphilarge, psit_c, psit_f, hpsit_c, hpsit_f, lpsit_c, lpsit_f!, phiLoc
  real(kind=8),dimension(:,:),allocatable:: matrix, locregCenter, dovrlp
  type(energy_terms) :: energs
  type(confpot_data),dimension(:),allocatable :: confdatarrtmp
  integer,dimension(:),allocatable:: norbsPerAtom, norbsPerLocreg
  character(len=*),parameter:: subname='pulay_correction'


  call nullify_orbitals_data(tmbder%orbs)
  call nullify_p2pComms(tmbder%comrp)
  call nullify_collective_comms(tmbder%collcom)
  call nullify_matrixDescriptors(tmbder%mad)
  call nullify_overlapParameters(tmbder%op)
  call nullify_p2pComms(tmbder%comon)

  norbu=3*tmb%orbs%norb
  norb=norbu
  norbd=0
  nspin=1
  tag=0

  allocate(norbsPerAtom(at%nat), stat=istat)
  call memocc(istat, norbsPerAtom, 'norbsPerAtom', subname)
  norb=0
  nlr=0
  do iat=1,at%nat
      ityp=at%iatype(iat)
      norbsPerAtom(iat)=3*input%lin%norbsPerType(ityp)
      norb=norb+3*input%lin%norbsPerType(ityp)
      nlr=nlr+input%lin%norbsPerType(ityp)
  end do


  allocate(norbsPerLocreg(nlr), stat=istat) 
  call memocc(istat, norbsPerLocreg, 'norbsPerLocreg', subname) 
  norbsPerLocreg=3 !should be norbsPerLocreg 
     


  call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, nspin, tmblarge%orbs%nspinor,&
       tmblarge%orbs%nkpts, tmblarge%orbs%kpts, tmblarge%orbs%kwgts, tmbder%orbs,.true.) !simple repartition

  iall=-product(shape(tmbder%orbs%onwhichatom))*kind(tmbder%orbs%inWhichLocreg) 
  deallocate(tmbder%orbs%onwhichatom, stat=istat) 
  call memocc(istat, iall, 'lorbs%onwhichatom', subname)
 
  call assignToLocreg2(iproc, nproc, tmbder%orbs%norb, tmbder%orbs%norb_par, at%nat, at%nat, &
       nspin, norbsPerAtom, rxyz, tmbder%orbs%onwhichatom)

  allocate(locregCenter(3,nlr), stat=istat)
  call memocc(istat, locregCenter, 'locregCenter', subname)
  ilr=0
  do iat=1,at%nat
      ityp=at%iatype(iat)
      do iorb=1,input%lin%norbsPerType(ityp)
          ilr=ilr+1
          locregCenter(:,ilr)=rxyz(:,iat)
      end do
  end do

  iall=-product(shape(tmbder%orbs%inWhichLocreg))*kind(tmbder%orbs%inWhichLocreg) 
  deallocate(tmbder%orbs%inWhichLocreg, stat=istat) 
  call memocc(istat, iall, 'lorbs%inWhichLocreg', subname) 
  call assignToLocreg2(iproc, nproc, tmbder%orbs%norb, tmbder%orbs%norb_par, at%nat, nlr, &
       nspin, norbsPerLocreg, locregCenter, tmbder%orbs%inwhichlocreg)

  iall=-product(shape(locregCenter))*kind(locregCenter)
  deallocate(locregCenter, stat=istat)
  call memocc(istat, iall, 'locregCenter', subname)


  !Calculate the dimension of psi
  call update_wavefunctions_size(tmblarge%lzd,tmbder%orbs,iproc,nproc)  !using tmblarge%lzd because the derivatives are also enlarged

  call initializeRepartitionOrbitals(iproc, nproc, tag, tmblarge%orbs, tmbder%orbs, tmblarge%lzd, tmbder%comrp)
  call init_collective_comms(iproc, nproc, tmbder%orbs, tmblarge%lzd, tmbder%collcom)

  call initCommsOrtho(iproc, nproc, nspin, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
       tmblarge%lzd, tmblarge%lzd, tmbder%orbs, 's', tmb%wfnmd%bpo, tmbder%op, tmbder%comon)

  allocate(tmbder%psi(max(tmbder%orbs%npsidim_orbs,tmbder%orbs%npsidim_comp)), stat=istat)
  call memocc(istat, tmbder%psi, 'tmbder%psi', subname)
  !if (tmblarge%orbs%npsidim_orbs> 0) call to_zero(tmblarge%orbs%npsidim_orbs,tmblarge%psi(1))
  call to_zero(tmblarge%orbs%npsidim_orbs,tmblarge%psi(1))

  call small_to_large_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmb%orbs, tmblarge%orbs, tmb%psi, tmblarge%psi)

  call getDerivativeBasisFunctions(iproc,nproc,tmblarge%lzd%hgrids(1),tmblarge%lzd,tmblarge%orbs,tmbder%orbs,tmbder%comrp,&
       max(tmblarge%orbs%npsidim_orbs,tmblarge%orbs%npsidim_comp),tmblarge%psi,tmbder%psi)

  ! modify the derivatives
  !!call derivatives_with_orthoconstraint(iproc, nproc, tmblarge, tmbder)

  ! Apply Hamiltonian to tmb%psi
  call local_potential_dimensions(tmblarge%lzd,tmblarge%orbs,denspot%dpbox%ngatherarr(0,1))

  allocate(lhphilarge(tmblarge%orbs%npsidim_orbs), stat=istat)
  call memocc(istat, lhphilarge, 'lhphilarge', subname)
  !if (tmblarge%orbs%npsidim_orbs > 0) call to_zero(tmblarge%orbs%npsidim_orbs,lhphilarge(1))
  call to_zero(tmblarge%orbs%npsidim_orbs,lhphilarge(1))


  call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
       tmblarge%comgp%nrecvbuf, tmblarge%comgp%recvbuf, tmblarge%comgp)

  allocate(confdatarrtmp(tmbder%orbs%norbp))
  call default_confinement_data(confdatarrtmp,tmbder%orbs%norbp)

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
  !call wait_p2p_communication(iproc, nproc, tmblarge%comgp)
  ! only potential
  call LocalHamiltonianApplication(iproc,nproc,at,tmblarge%orbs,&
       tmblarge%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmblarge%psi,lhphilarge,&
       energs,SIC,GPU,2,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmblarge%comgp)

  call timing(iproc,'glsynchham1','ON') !lr408t
  call SynchronizeHamiltonianApplication(nproc,tmblarge%orbs,tmblarge%lzd,GPU,lhphilarge,&
       energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
  call timing(iproc,'glsynchham1','OF') !lr408t
  deallocate(confdatarrtmp)

  ! Calculate matrix
  allocate(matrix(tmbder%orbs%norb,tmblarge%orbs%norb), stat=istat) 
  call memocc(istat, matrix, 'matrix', subname)
  call to_zero(tmbder%orbs%norb*tmblarge%orbs%norb,matrix(1,1))

  allocate(hpsit_c(tmblarge%collcom%ndimind_c))
  call memocc(istat, hpsit_c, 'hpsit_c', subname)
  allocate(hpsit_f(7*tmblarge%collcom%ndimind_f))
  call memocc(istat, hpsit_f, 'hpsit_f', subname)
  allocate(psit_c(tmbder%collcom%ndimind_c))
  call memocc(istat, psit_c, 'psit_c', subname)
  allocate(psit_f(7*tmbder%collcom%ndimind_f))
  call memocc(istat, psit_f, 'psit_f', subname)

  call transpose_localized(iproc, nproc, tmbder%orbs,  tmbder%collcom, &
       tmbder%psi, psit_c, psit_f, tmblarge%lzd)
  call transpose_localized(iproc, nproc, tmblarge%orbs,  tmblarge%collcom, &
       lhphilarge, hpsit_c, hpsit_f, tmblarge%lzd)

  call calculate_pulay_overlap(iproc, nproc, tmbder%orbs, tmblarge%orbs, tmbder%collcom, &
       tmblarge%collcom, psit_c, hpsit_c, psit_f, hpsit_f, matrix)
 
 
  !DEBUG
  !!if(iproc==0)then
  !!do iorb = 1, tmbder%orbs%norb
  !!   print *,'Hamiltonian of derivative: ',iorb, (matrix(iorb,iiorb),iiorb=1,tmblarge%orbs%norb)
  !!end do
  !!end if
  !END DEBUG

  allocate(dovrlp(tmbder%orbs%norb,tmblarge%orbs%norb), stat=istat) 
  call memocc(istat, dovrlp, 'dovrlp', subname)
  call to_zero(tmbder%orbs%norb*tmblarge%orbs%norb,dovrlp(1,1))
  allocate(lpsit_c(tmblarge%collcom%ndimind_c))
  call memocc(istat, lpsit_c, 'lpsit_c', subname)
  allocate(lpsit_f(7*tmblarge%collcom%ndimind_f))
  call memocc(istat, lpsit_f, 'lpsit_f', subname)

  call transpose_localized(iproc, nproc, tmblarge%orbs,  tmblarge%collcom, &
       tmblarge%psi, lpsit_c, lpsit_f, tmblarge%lzd)
  
  call calculate_pulay_overlap(iproc, nproc, tmbder%orbs, tmblarge%orbs, tmbder%collcom, &
       tmblarge%collcom, psit_c, lpsit_c, psit_f, lpsit_f, dovrlp)

  !!!DEBUG
  !!!Check if derivatives are orthogonal to functions
  !!if(iproc==0)then
  !!  do iorb = 1, tmbder%orbs%norb
  !!     !print *,'overlap of derivative: ',iorb, (dovrlp(iorb,iiorb),iiorb=1,tmblarge%orbs%norb)
  !!     do iiorb=1,tmbder%orbs%norb
  !!         write(*,*) iorb, iiorb, dovrlp(iorb,iiorb)
  !!     end do
  !!  end do
  !!end if
  !!!END DEBUG

  iall=-product(shape(hpsit_c))*kind(hpsit_c)
  deallocate(hpsit_c, stat=istat)
  call memocc(istat, iall, 'hpsit_c', subname)
  iall=-product(shape(hpsit_f))*kind(hpsit_f)
  deallocate(hpsit_f, stat=istat)
  call memocc(istat, iall, 'hpsit_f', subname)
  iall=-product(shape(psit_c))*kind(psit_c)
  deallocate(psit_c, stat=istat)
  call memocc(istat, iall, 'psit_c', subname)
  iall=-product(shape(psit_f))*kind(psit_f)
  deallocate(psit_f, stat=istat)
  call memocc(istat, iall, 'psit_f', subname)

  ! Calculate Pulay correction
  ! note since the basis functions are real, only multiply by two instead of taking the real conjugate
  call to_zero(3*at%nat, fpulay(1,1))
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      do jjorb=1,tmbder%orbs%norb
          jat=tmbder%orbs%onwhichatom(jjorb)
          jdir=mod(jjorb-1,3) + 1 ! get direction: x=1, y=2 or z=3 
          jorbsmall=ceiling(dble(jjorb)/3.d0)
          do kkorb=1,tmblarge%orbs%norb
             fpulay(jdir,jat) = fpulay(jdir,jat) - &
              4*tmb%wfnmd%coeff(jorbsmall,iiorb)*tmb%wfnmd%coeff(kkorb,iiorb)* &
              (matrix(jjorb,kkorb) - tmblarge%orbs%eval(iiorb)*dovrlp(jjorb,kkorb))
              !!do llorb=1,tmbder%orbs%norb
              !!    lat=tmbder%orbs%onwhichatom(llorb)
              !!    ldir=mod(llorb-1,3) + 1 ! get direction: x=1, y=2 or z=3 
              !!    lorbsmall=ceiling(dble(llorb)/3.d0)
              !!    if (lat/=jat) fpulay(ldir,lat) = fpulay(ldir,lat) + &
              !!    2*tmb%wfnmd%coeff(jorbsmall,iiorb)*tmb%wfnmd%coeff(lorbsmall,iiorb)*dovrlp(llorb,lorbsmall)* &
              !!     (matrix(jjorb,lorbsmall) - tmblarge%orbs%eval(iiorb)*dovrlp(jjorb,lorbsmall))
              !!end do
          end do
      end do
  end do
  call mpiallred(fpulay(1,1), 3*at%nat, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  if(iproc==0) then
       do iat=1,at%nat
           write(*,'(a,i5,3es16.6)') 'iat, fpulay', iat, fpulay(1:3,iat)
       end do
  end if

  iall=-product(shape(lpsit_c))*kind(lpsit_c)
  deallocate(lpsit_c, stat=istat)
  call memocc(istat, iall, 'lpsit_c', subname)
  iall=-product(shape(lpsit_f))*kind(lpsit_f)
  deallocate(lpsit_f, stat=istat)
  call memocc(istat, iall, 'lpsit_f', subname)
  iall=-product(shape(dovrlp))*kind(dovrlp)
  deallocate(dovrlp, stat=istat)
  call memocc(istat, iall, 'dovrlp', subname)

  iall=-product(shape(matrix))*kind(matrix)
  deallocate(matrix, stat=istat)
  call memocc(istat, iall, 'matrix', subname)

  iall=-product(shape(lhphilarge))*kind(lhphilarge)
  deallocate(lhphilarge, stat=istat)
  call memocc(istat, iall, 'lhphilarge', subname)

  iall=-product(shape(tmblarge%lzd%doHamAppl))*kind(tmblarge%lzd%doHamAppl)
  deallocate(tmblarge%lzd%doHamAppl, stat=istat)
  call memocc(istat, iall, 'tmblarge%lzd%doHamAppl', subname)

  iall=-product(shape(denspot%pot_work))*kind(denspot%pot_work)
  deallocate(denspot%pot_work, stat=istat)
  call memocc(istat, iall, 'denspot%pot_work', subname)

  if(iproc==0) write(*,'(1x,a)') 'done.'

  iall=-product(shape(norbsPerAtom))*kind(norbsPerAtom)
  deallocate(norbsPerAtom, stat=istat)
  call memocc(istat, iall, 'norbsPerAtom', subname)

  iall=-product(shape(norbsPerLocreg))*kind(norbsPerLocreg)
  deallocate(norbsPerLocreg, stat=istat)
  call memocc(istat, iall, 'norbsPerLocreg', subname)

  iall=-product(shape(tmbder%psi))*kind(tmbder%psi)
  deallocate(tmbder%psi, stat=istat)
  call memocc(istat, iall, 'tmbder%psi', subname)

  iall=-product(shape(tmbder%comon%comarr))*kind(tmbder%comon%comarr)
  deallocate(tmbder%comon%comarr, stat=istat)
  call memocc(istat, iall, 'tmbder%comon%comarr', subname)

  iall=-product(shape(tmbder%comon%noverlaps))*kind(tmbder%comon%noverlaps)
  deallocate(tmbder%comon%noverlaps, stat=istat)
  call memocc(istat, iall, 'tmbder%comon%noverlaps', subname)

  call deallocate_orbitals_data(tmbder%orbs, subname)
  call deallocate_p2pComms(tmbder%comrp, subname)
  call deallocate_collective_comms(tmbder%collcom, subname)
  call deallocate_overlapParameters(tmbder%op, subname)

end subroutine pulay_correction



subroutine derivative_coeffs_from_standard_coeffs(orbs, tmb, tmbder)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(in):: tmb
  type(DFT_wavefunction),intent(out):: tmbder

  ! Local variables
  integer:: iorb, jorb, jjorb

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
  integer,intent(in):: iproc, nproc
  type(DFT_wavefunction),intent(in):: tmb
  type(DFT_wavefunction),intent(inout):: tmbder

  ! Local variables
  integer:: i0, j0, ii, jj, ipt, i, iiorb, jjorb, istat, iall, j
  real(8),dimension(:),allocatable:: psit_c, psit_f, psidert_c, psidert_f
  real(8),dimension(:,:),allocatable:: matrix
  character(len=*),parameter:: subname='derivatives_with_orthoconstraint'


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
