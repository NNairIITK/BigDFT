subroutine linearScaling(iproc,nproc,Glr,orbs,comms,tmb,tmbder,at,input,hx,hy,hz,&
           rxyz,fion,fdisp,denspot,rhopotold,nlpspd,proj,GPU,&
           energs,scpot,psi,energy)
use module_base
use module_types
use module_interfaces, exceptThisOne => linearScaling
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(locreg_descriptors),intent(in) :: Glr
type(orbitals_data),intent(inout):: orbs
type(communications_arrays),intent(in) :: comms
type(atoms_data),intent(inout):: at
type(input_variables),intent(in):: input
real(8),dimension(3,at%nat),intent(inout):: rxyz
real(8),dimension(3,at%nat),intent(in):: fion, fdisp
type(DFT_local_fields), intent(inout) :: denspot
real(gp), dimension(:), intent(inout) :: rhopotold
type(nonlocal_psp_descriptors),intent(in):: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
type(GPU_pointers),intent(in out):: GPU
type(energy_terms),intent(inout) :: energs
real(gp),intent(in):: hx,hy,hz
logical,intent(in):: scpot
real(8),dimension(:),pointer,intent(out):: psi
real(gp), dimension(:), pointer :: rho,pot
real(8),intent(out):: energy
type(DFT_wavefunction),intent(inout),target:: tmb
type(DFT_wavefunction),intent(inout),target:: tmbder

type(linear_scaling_control_variables):: lscv
real(8):: pnrm,trace,fnrm_tmb
integer:: infoCoeff,istat,iall,it_scc,ilr,tag,itout,iorb,scf_mode,info_scf,nsatur
character(len=*),parameter:: subname='linearScaling'
real(8),dimension(:),pointer :: psit
real(8),dimension(:),allocatable:: rhopotold_out
real(8):: energyold, energyDiff, energyoldout
type(mixrhopotDIISParameters):: mixdiis
type(localizedDIISParameters):: ldiis, ldiis_coeff
type(DFT_wavefunction),pointer:: tmbmix
logical:: check_whether_derivatives_to_be_used,coeffs_copied, first_time_with_der,calculate_overlap_matrix
integer:: jorb, jjorb, iiat,nit_highaccur
real(8),dimension(:,:),allocatable:: overlapmatrix
real(8),dimension(:),allocatable :: locrad_tmp, eval
type(DFT_wavefunction):: tmblarge, tmblargeder, tmblarge2
real(8),dimension(:,:),allocatable:: locregCenter, locregCenterTemp, kernel, Umat
real(8),dimension(:),pointer:: lhphilarge, lhphilargeold, lphilargeold, lhphilargeder, lhphilargeoldder, lphilargeoldder
real(8),dimension(3,at%nat):: fpulay


  call timing(iproc,'linscalinit','ON') !lr408t

  if(iproc==0) then
      write(*,'(1x,a)') repeat('*',84)
      write(*,'(1x,a)') '****************************** LINEAR SCALING VERSION ******************************'
  end if



  allocate(lscv%locrad(tmb%lzd%nlr), stat=istat)
  call memocc(istat, lscv%locrad, 'lscv%locrad', subname)


  ! Allocate the old charge density (used to calculate the variation in the charge density)
  allocate(rhopotold_out(max(glr%d%n1i*glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin), stat=istat)
  call memocc(istat, rhopotold_out, 'rhopotold_out', subname)

  if(input%lin%nItInguess>0) then
      tmb%wfnmd%bs%communicate_phi_for_lsumrho=.true.
      tmb%wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE

      do ilr=1,tmb%lzd%nlr
          lscv%locrad(ilr)=max(input%lin%locrad_lowaccuracy(ilr),tmb%lzd%llr(ilr)%locrad)
      end do

      if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE) then
          rhopotold_out=rhopotold
      end if

      if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
          rhopotold_out=denspot%rhov
      end if

      ! Copy the current potential
      if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
           call dcopy(max(Glr%d%n1i*Glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotOld(1), 1)
      end if
  end if



  ! Allocate the communications buffers needed for the communications of teh potential and
  ! post the messages. This will send to each process the part of the potential that this process
  ! needs for the application of the Hamlitonian to all orbitals on that process.
  call allocateCommunicationsBuffersPotential(tmb%comgp, subname)
  call allocateCommunicationsBuffersPotential(tmbder%comgp, subname)

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
  lscv%reduce_convergence_tolerance=.false.
  tmb%wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE
  lscv%lowaccur_converged=.false.
  lscv%info_basis_functions=-1
  lscv%idecrease=0
  lscv%increase_locreg=0.d0
  lscv%decrease_factor_total=1.d10 !initialize to some large value
  lscv%ifail=0
  lscv%enlarge_locreg=.false.

  ! tmbmix is the types we use for the mixing. It will point to either tmb if we don't use the derivatives
  ! or to tmbder if we use the derivatives.
  if(input%lin%useDerivativeBasisFunctions) then
      tmbmix => tmbder
  else
      tmbmix => tmb
  end if

  !!! Check whether it is possible to have variable localization regions or not.
  !!if(tmb%wfnmd%bs%nit_unitary_loop==-1 .and. tmb%wfnmd%bs%locreg_enlargement==1.d0) then
      lscv%variable_locregs=.false.
  !!else
  !!    lscv%variable_locregs=.true.
  !!end if


  !!! Allocate the communication arrays for the calculation of the charge density.
  call allocateCommunicationbufferSumrho(iproc, tmb%comsr, subname)
  call allocateCommunicationbufferSumrho(iproc, tmbder%comsr, subname)

  !!! Initialize DIIS...
  !!!!if(.not.lscv%lowaccur_converged) then
  !!    call initializeDIIS(input%lin%DIIS_hist_lowaccur, tmb%lzd, tmb%orbs, tmb%orbs%norb, ldiis)
  !!    ldiis%DIISHistMin=0
  !!    ldiis%DIISHistMax=input%lin%DIIS_hist_lowaccur
  !!    ldiis%alphaSD=input%lin%alphaSD
  !!    ldiis%alphaDIIS=input%lin%alphaDIIS
  !!    ldiis%icountSDSatur=0
  !!    ldiis%icountSwitch=0
  !!    ldiis%icountDIISFailureTot=0
  !!    ldiis%icountDIISFailureCons=0
  !!    ldiis%is=0
  !!    ldiis%switchSD=.false.
  !!    ldiis%trmin=1.d100
  !!    ldiis%trold=1.d100
  !!!!end if


  ! just to be sure...
  nullify(tmb%psit_c)
  nullify(tmb%psit_f)
  nullify(tmbder%psit_c)
  nullify(tmbder%psit_f)


  allocate(eval(tmb%orbs%norb), stat=istat)
  call memocc(istat, eval, 'eval', subname)
  call vcopy(tmb%orbs%norb, tmb%orbs%eval(1), 1, eval(1), 1)
  call timing(iproc,'linscalinit','OF') !lr408t
  ! This is the main outer loop. Each iteration of this loop consists of a first loop in which the basis functions
  ! are optimized and a consecutive loop in which the density is mixed.
  coeffs_copied=.false.
  first_time_with_der=.false.
  nit_highaccur=0

  nsatur=0

  ! Set to zero the large wavefunction. Later only the inner part will be filled. It must be made sure
  ! that the outer part is not modified!
  if (tmblarge%orbs%npsidim_orbs > 0) call to_zero(tmblarge%orbs%npsidim_orbs,tmblarge%psi(1))

  outerLoop: do itout=1,input%lin%nit_lowaccuracy+input%lin%nit_highaccuracy

      ! First to some initialization and determine the value of some control parameters.

      ! The basis functions shall be optimized
      tmb%wfnmd%bs%update_phi=.true.

      ! Convergence criterion for the self consistency loop
      lscv%self_consistent=input%lin%convCritMix

      ! Check whether the low accuracy part (i.e. with strong confining potential) has converged.
      call check_whether_lowaccuracy_converged(itout, input, lscv)

      ! Check whether the derivatives shall be used or not.
      lscv%withder=check_whether_derivatives_to_be_used(input, itout, lscv)

      if(lscv%withder .and. lscv%lowaccur_converged .and. .not.coeffs_copied) then
          call to_zero(tmbder%orbs%norb*orbs%norb, tmbder%wfnmd%coeff(1,1))
          do iorb=1,orbs%norb
              jjorb=0
              do jorb=1,tmbder%orbs%norb,4
                  jjorb=jjorb+1
                  tmbder%wfnmd%coeff(jorb,iorb)=tmb%wfnmd%coeff(jjorb,iorb)
              end do
          end do
          coeffs_copied=.true.
      end if

      ! Set all remaining variables that we need for the optimizations of the basis functions and the mixing.
      call set_optimization_variables(input, at, tmb%orbs, tmb%lzd%nlr, tmb%orbs%onwhichatom, &
           tmb%confdatarr, tmb%wfnmd, lscv)
      call set_optimization_variables(input, at, tmbder%orbs, tmb%lzd%nlr, tmbder%orbs%onwhichatom, &
           tmbder%confdatarr, tmbder%wfnmd, lscv)

      if(lscv%lowaccur_converged) nit_highaccur=nit_highaccur+1
      if(nit_highaccur==1) lscv%enlarge_locreg=.true.

      ! Adjust the confining potential if required.
      call adjust_locregs_and_confinement(iproc, nproc, hx, hy, hz, &
           input, tmb, tmbder, denspot, ldiis, lscv)

      ! Some special treatement if we are in the high accuracy part
      call adjust_DIIS_for_high_accuracy(input, tmb, denspot, ldiis, mixdiis, lscv)

      if(lscv%withder) then
          call initialize_DIIS_coeff(3, tmbder, orbs, ldiis_coeff)
      else
          call initialize_DIIS_coeff(3, tmb, orbs, ldiis_coeff)
      end if

      ! Now all initializations are done...


      allocate(locregCenter(3,tmb%lzd%nlr), stat=istat)
      call memocc(istat, locregCenter, 'locregCenter', subname)
      allocate(locrad_tmp(tmb%lzd%nlr), stat=istat)
      call memocc(istat, locrad_tmp, 'locrad_tmp', subname)


      if(nit_highaccur==1) then
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


      if(itout==1 .or. nit_highaccur==1) then
          call create_large_tmbs(iproc, nproc, tmb, eval, denspot, input, at, rxyz, lscv%lowaccur_converged, &
               tmblarge, lhphilarge, lhphilargeold, lphilargeold)
      end if



      if(lscv%withder) then
          call update_locreg(iproc, nproc, tmb%lzd%nlr, locrad_tmp, tmbder%orbs%inwhichlocreg, locregCenter, tmb%lzd%glr, &
               tmb%wfnmd%bpo, .false., denspot%dpbox%nscatterarr, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
               tmbder%orbs, tmblargeder%lzd, tmblargeder%orbs, tmblargeder%op, tmblargeder%comon, &
               tmblargeder%comgp, tmblargeder%comsr, tmblargeder%mad, tmblargeder%collcom)
          call allocate_auxiliary_basis_function(max(tmblargeder%orbs%npsidim_comp,tmblargeder%orbs%npsidim_orbs), subname, &
               tmblargeder%psi, lhphilargeder, lhphilargeoldder, lphilargeoldder)
          call copy_basis_performance_options(tmbder%wfnmd%bpo, tmblargeder%wfnmd%bpo, subname)
          call copy_orthon_data(tmbder%orthpar, tmblargeder%orthpar, subname)
          tmblargeder%wfnmd%nphi=tmblargeder%orbs%npsidim_orbs
          tmblargeder%can_use_transposed=.false.
          nullify(tmblargeder%psit_c)
          nullify(tmblargeder%psit_f)
          allocate(tmblargeder%confdatarr(tmblargeder%orbs%norbp), stat=istat)
          call vcopy(tmb%orbs%norb, tmbder%orbs%onwhichatom(1), 1, tmblargeder%orbs%onwhichatom(1), 1)
          if(.not.lscv%lowaccur_converged) then
              call define_confinement_data(tmblargeder%confdatarr,tmblargeder%orbs,rxyz,at,&
                   tmblargeder%lzd%hgrids(1),tmblargeder%lzd%hgrids(2),tmblargeder%lzd%hgrids(3),&
                   input%lin%ConfPotOrder,input%lin%potentialPrefac_lowaccuracy,tmblargeder%lzd,tmblargeder%orbs%onwhichatom)
          else
              call define_confinement_data(tmblargeder%confdatarr,tmblargeder%orbs,rxyz,at,&
                   tmblargeder%lzd%hgrids(1),tmblargeder%lzd%hgrids(2),tmblargeder%lzd%hgrids(3),&
                   input%lin%ConfPotOrder,input%lin%potentialPrefac_highaccuracy,tmblargeder%lzd,tmblargeder%orbs%onwhichatom)
          end if
      end if

      if(itout==1) then
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

      !!if (tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_TRACE) then
      !!    tmb%wfnmd%bs%maxdev_ortho=1.d-20
      !!else if (tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
      !!    tmb%wfnmd%bs%maxdev_ortho=1.d0
      !!end if

      !!if(iproc==0) write(*,*) 'WARNING: input has wrong intent!'
      !!if(lscv%lowaccur_converged) then
      !!    if(iproc==0) write(*,*) 'WARNING: set DIIS history to 0!'
      !!    input%lin%DIISHistMax=0
      !!end if

      if(itout>1) call deallocateDIIS(ldiis)
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

      ! The self consistency cycle. Here we try to get a self consistent density/potential.
      ! In the first lscv%nit_scc_when_optimizing iteration, the basis functions are optimized, whereas in the remaining
      ! iteration the basis functions are fixed.
      do it_scc=1,lscv%nit_scc

          scf_mode=input%lin%scf_mode

          ! Do not update the TMB if it_scc>lscv%nit_scc_when_optimizing
          if(it_scc>lscv%nit_scc_when_optimizing) tmb%wfnmd%bs%update_phi=.false.

          ! Stop the optimization if it seems to saturate
          if(nsatur>=tmb%wfnmd%bs%nsatur_outer) then
              tmb%wfnmd%bs%update_phi=.false.
              if(it_scc==1) then
                  tmb%can_use_transposed=.false.
                  allocate(overlapmatrix(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
                  call memocc(istat, overlapmatrix, 'overlapmatrix', subname)
              end if
          end if

          !!call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
          !!     tmb%comgp%nrecvbuf, tmb%comgp%recvbuf, tmb%comgp)
          if(lscv%withder) then
              call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
                   tmbder%comgp%nrecvbuf, tmbder%comgp%recvbuf, tmbder%comgp)
          end if

         ! Improve the trace minimizing orbitals.
          if(tmb%wfnmd%bs%update_phi) then
              if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
                  do iorb=1,orbs%norb
                      call dcopy(tmb%orbs%norb, tmb%wfnmd%coeff_proj(1,iorb), 1, tmb%wfnmd%coeff(1,iorb), 1)
                  end do
              end if

              call getLocalizedBasis(iproc,nproc,at,orbs,rxyz,denspot,GPU,trace,fnrm_tmb,lscv%info_basis_functions,&
                  nlpspd,proj,ldiis,input%SIC,tmb, tmblarge, lhphilarge)
              if(lscv%info_basis_functions>0) then
                  nsatur=nsatur+1
              else
                  !!nsatur=0
              end if
              tmb%can_use_transposed=.false. !since basis functions have changed...
              tmbder%can_use_transposed=.false. !since basis functions have changed...

              tmb%wfnmd%nphi=tmb%orbs%npsidim_orbs
              !reset counter for optimization of coefficients (otherwise step size will be decreases...)
              tmb%wfnmd%it_coeff_opt=0
              tmbder%wfnmd%it_coeff_opt=0
              tmb%wfnmd%alpha_coeff=.2d0 !reset to default value
              tmbder%wfnmd%alpha_coeff=.2d0 !reset to default value

              !!! Reset DIIS if we are at the first iteration of the high accuracy regime
              !!! or if DIIS became unstable in the previous optimization of the TMBs.
              !!!!$if(nit_highaccur<=1 .or. ldiis%isx<input%lin%DIISHistMax) then
              !!    call deallocateDIIS(ldiis)
              !!    call initializeDIIS(input%lin%DIISHistMax, tmb%lzd, tmb%orbs, tmb%orbs%norb, ldiis)
              !!    ldiis%DIISHistMin=input%lin%DIISHistMin
              !!    ldiis%DIISHistMax=input%lin%DIISHistMax
              !!    ldiis%alphaSD=input%lin%alphaSD
              !!    ldiis%alphaDIIS=input%lin%alphaDIIS
              !!    ldiis%icountSDSatur=0
              !!    ldiis%icountSwitch=0
              !!    ldiis%icountDIISFailureTot=0
              !!    ldiis%icountDIISFailureCons=0
              !!    ldiis%is=0
              !!    ldiis%switchSD=.false.
              !!    ldiis%trmin=1.d100
              !!    ldiis%trold=1.d100
              !!!!$else
              !!!!$    ! Keep the history in the high accuracy case.
              !!!!$    ! Since the potential changes, the values of ldiis%trmin must be reset.
              !!!!$    ldiis%switchSD=.false.
              !!!!$    ldiis%trmin=1.d100
              !!!!$    ldiis%trold=1.d100
              !!!!$    ldiis%icountSDSatur=0
              !!!!$    ldiis%icountSwitch=0
              !!!!$    ldiis%icountDIISFailureTot=0
              !!!!$    ldiis%icountDIISFailureCons=0
              !!!!$end if
          end if

          if((lscv%locreg_increased .or. (lscv%variable_locregs .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY)) &
              .and. tmb%wfnmd%bs%update_phi) then
              ! Redefine some quantities if the localization region has changed.
              if(lscv%withder) then
                  call redefine_locregs_quantities(iproc, nproc, hx, hy, hz, tmb%lzd%llr(:)%locrad, &
                       .false., tmb%lzd, tmb, tmbder, denspot)
                  call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
                       tmbder%comgp%nrecvbuf, tmbder%comgp%recvbuf, tmbder%comgp)
              end if
          end if

          ! Decide whether we have to use the derivatives or not.
          if(lscv%withder) then
              tmbmix => tmbder
          else
              tmbmix => tmb
          end if

          ! Build the derivatives if required.
          if(tmb%wfnmd%bs%update_phi .or. it_scc==0) then
              if(tmbmix%wfnmd%bs%use_derivative_basis) then
                  if((lscv%locreg_increased .or. &
                      (lscv%variable_locregs .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY)) &
                      .and. tmb%wfnmd%bs%update_phi) then
                      call deallocate_p2pComms(tmbder%comrp, subname)
                      call initializeRepartitionOrbitals(iproc, nproc, tag, tmb%orbs, tmbder%orbs, tmb%lzd, tmbder%comrp)
                      tmbmix => tmbder
                  end if
                  if(iproc==0) write(*,'(1x,a)',advance='no') 'calculating derivative basis functions...'
                  call getDerivativeBasisFunctions(iproc,nproc,hx,tmb%lzd,tmb%orbs,tmbmix%orbs,tmbmix%comrp,&
                       max(tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp),tmb%psi,tmbmix%psi)
                  if(iproc==0) write(*,'(a)') 'done.'
                  !! TEST ###############################################################################################
                  !write(*,*) 'test: orthonormalize derivatives'
                  !!call orthonormalizeLocalized(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%orthpar%nItOrtho, &
                  !!     tmbder%orbs, tmbder%op, tmbder%comon, tmb%lzd, &
                  !!     tmbder%mad, tmbder%collcom, tmbder%orthpar, tmbder%wfnmd%bpo, tmbder%psi, tmbder%psit_c, tmbder%psit_f, &
                  !!     tmbder%can_use_transposed)
                  !!if(tmbder%can_use_transposed) then
                  !!    ! This is not optimal, these quantities will be recalculated...
                  !!    iall = -product(shape(tmbder%psit_c))*kind(tmbder%psit_c)
                  !!    deallocate(tmbder%psit_c,stat=istat)
                  !!    call memocc(istat,iall,'tmbder%psit_c',subname)
                  !!    iall = -product(shape(tmbder%psit_f))*kind(tmbder%psit_f)
                  !!    deallocate(tmbder%psit_f,stat=istat)
                  !!    call memocc(istat,iall,'tmbder%psit_f',subname)
                  !!end if
                  !! END TEST ###########################################################################################
              else
                  call dcopy(tmb%wfnmd%nphi, tmb%psi(1), 1, tmbmix%psi(1), 1)
              end if

              allocate(overlapmatrix(tmbmix%orbs%norb,tmbmix%orbs%norb), stat=istat)
              call memocc(istat, overlapmatrix, 'overlapmatrix', subname)


          end if

          ! Only communicate the TMB for sumrho if required (i.e. only if the TMB were optimized).
          if(it_scc<=lscv%nit_scc_when_optimizing) then
              tmbmix%wfnmd%bs%communicate_phi_for_lsumrho=.true.
              calculate_overlap_matrix=.true.
          else
              tmbmix%wfnmd%bs%communicate_phi_for_lsumrho=.false.
              calculate_overlap_matrix=.false.
          end if

          if(lscv%withder .and. .not.first_time_with_der) then
              first_time_with_der=.true.
              !scf_mode=LINEAR_MIXDENS_SIMPLE
              scf_mode=input%lin%scf_mode
              call transform_coeffs_to_derivatives(iproc, nproc, orbs, tmb%lzd, tmb, tmbder)
              tmb%wfnmd%alpha_coeff=1.d-2
              tmbder%wfnmd%alpha_coeff=1.d-2
          else
              scf_mode=input%lin%scf_mode
          end if

          ! Calculate the coefficients
          if(.not.lscv%withder) then
              call get_coeff(iproc,nproc,scf_mode,tmb%lzd,orbs,at,rxyz,denspot,GPU,infoCoeff,energs%ebs,nlpspd,proj,&
                   input%SIC,tmbmix,tmb,pnrm,overlapmatrix,calculate_overlap_matrix,&
                   tmblarge, lhphilarge, ldiis_coeff)
          else
              call get_coeff(iproc,nproc,scf_mode,tmb%lzd,orbs,at,rxyz,denspot,GPU,infoCoeff,energs%ebs,nlpspd,proj,&
                   input%SIC,tmbmix,tmb,pnrm,overlapmatrix,calculate_overlap_matrix,&
                   tmblargeder, lhphilargeder, ldiis_coeff)
          end if

          !!if(itout==1 .or. itout==15) then
          !!    if(Iproc==0) WRITE(*,*) 'WRITE KERNEL TO FILE'
          !!    do iorb=1,tmb%orbs%norb
          !!        do jorb=1,tmb%orbs%norb
          !!        if(itout==1) write(200,*) iorb,jorb,tmbmix%wfnmd%density_kernel(jorb,iorb)
          !!        if(itout==15) write(300,*) iorb,jorb,tmbmix%wfnmd%density_kernel(jorb,iorb)
          !!        end do
          !!    end do
          !!end if
              


          ! Calculate the total energy.
          energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
          energyDiff=energy-energyold
          energyold=energy


          ! Calculate the charge density.
          call sumrhoForLocalizedBasis2(iproc, nproc, &
               tmb%lzd, input, hx, hy ,hz, tmbmix%orbs, tmbmix%comsr, &
               tmbmix%wfnmd%density_kernel, Glr%d%n1i*Glr%d%n2i*denspot%dpbox%n3d, &
               denspot%rhov, at, denspot%dpbox%nscatterarr)

          ! Mix the density.
          if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE) then
           lscv%compare_outer_loop = pnrm<lscv%self_consistent .or. it_scc==lscv%nit_scc
           call mix_main(iproc, nproc, lscv%mix_hist, lscv%compare_outer_loop, input, glr, lscv%alpha_mix, &
                denspot, mixdiis, rhopotold, rhopotold_out, pnrm, lscv%pnrm_out)
          end if

          ! Calculate the new potential.
          if(iproc==0) write(*,'(1x,a)') '---------------------------------------------------------------- Updating potential.'
          call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)

          ! Mix the potential
          if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
           lscv%compare_outer_loop = pnrm<lscv%self_consistent .or. it_scc==lscv%nit_scc
           call mix_main(iproc, nproc, lscv%mix_hist, lscv%compare_outer_loop, input, glr, lscv%alpha_mix, &
                denspot, mixdiis, rhopotold, rhopotold_out, pnrm, lscv%pnrm_out)
          end if

          ! Make sure that the previous communication is complete (only do that if this check
          ! for completeness has not been done in get_coeff)
          if(tmbmix%wfnmd%bs%use_derivative_basis .and. .not.tmb%wfnmd%bs%update_phi) then
              call wait_p2p_communication(iproc, nproc, tmb%comgp)
          end if
          if(lscv%withder) then
              call wait_p2p_communication(iproc, nproc, tmbder%comgp)
          end if

          ! Write some informations.
          call printSummary(iproc, it_scc, lscv%info_basis_functions, &
               infoCoeff, pnrm, energy, energyDiff, input%lin%scf_mode)
          if(pnrm<lscv%self_consistent) then
              info_scf=it_scc
              lscv%reduce_convergence_tolerance=.true.
              exit
          else
              info_scf=-1
              lscv%reduce_convergence_tolerance=.false.
          end if

          if(nsatur<tmb%wfnmd%bs%nsatur_outer .and. it_scc<lscv%nit_scc_when_optimizing) then
              ! Deallocate the transposed TMBs
              if(tmbmix%can_use_transposed) then
                  iall=-product(shape(tmbmix%psit_c))*kind(tmbmix%psit_c)
                  deallocate(tmbmix%psit_c, stat=istat)
                  call memocc(istat, iall, 'tmbmix%psit_c', subname)
                  iall=-product(shape(tmbmix%psit_f))*kind(tmbmix%psit_f)
                  deallocate(tmbmix%psit_f, stat=istat)
                  call memocc(istat, iall, 'tmbmix%psit_f', subname)
              end if
              write(*,*) 'deallocating overlapmatrix'
              iall=-product(shape(overlapmatrix))*kind(overlapmatrix)
              deallocate(overlapmatrix, stat=istat)
              call memocc(istat, iall, 'overlapmatrix', subname)
          end if

      end do


    if(lscv%withder) then
        call destroy_new_locregs(iproc, nproc, tmblargeder)
        call deallocate_auxiliary_basis_function(subname, tmblargeder%psi, lhphilargeder, lhphilargeoldder, lphilargeoldder)
        if(tmblargeder%can_use_transposed) then
            iall=-product(shape(tmblargeder%psit_c))*kind(tmblargeder%psit_c)
            deallocate(tmblargeder%psit_c, stat=istat)
            call memocc(istat, iall, 'tmblargeder%psit_c', subname)
            iall=-product(shape(tmblargeder%psit_f))*kind(tmblargeder%psit_f)
            deallocate(tmblargeder%psit_f, stat=istat)
            call memocc(istat, iall, 'tmblargeder%psit_f', subname)
        end if
        deallocate(tmblargeder%confdatarr, stat=istat)
    end if

    

      iall=-product(shape(locregCenter))*kind(locregCenter)
      deallocate(locregCenter, stat=istat)
      call memocc(istat, iall, 'locregCenter', subname)
      iall=-product(shape(locrad_tmp))*kind(locrad_tmp)
      deallocate(locrad_tmp, stat=istat)
      call memocc(istat, iall, 'locrad_tmp', subname)



      call deallocateDIIS(ldiis_coeff)

      iall=-product(shape(overlapmatrix))*kind(overlapmatrix)
      deallocate(overlapmatrix, stat=istat)
      call memocc(istat, iall, 'overlapmatrix', subname)


      if(tmbmix%can_use_transposed) then
          iall=-product(shape(tmbmix%psit_c))*kind(tmbmix%psit_c)
          deallocate(tmbmix%psit_c, stat=istat)
          call memocc(istat, iall, 'tmbmix%psit_c', subname)
          iall=-product(shape(tmbmix%psit_f))*kind(tmbmix%psit_f)
          deallocate(tmbmix%psit_f, stat=istat)
          call memocc(istat, iall, 'tmbmix%psit_f', subname)
      end if

      ! Print out values related to two iterations of the outer loop.
      if(iproc==0) then
          write(*,'(3x,a,7es18.10)') 'ebs, ehart, eexcu, vexcu, eexctX, eion, edisp', &
              energs%ebs, energs%eh, energs%exc, energs%evxc, energs%eexctX, energs%eion, energs%edisp
          !if(trim(input%lin%mixingMethod)=='dens') then
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

  end do outerLoop

  !!! Calculate Pulay correction to the forces
  !!call pulay_correction(iproc, nproc, input, orbs, at, rxyz, nlpspd, proj, input%SIC, denspot, GPU, tmb, &
  !!         tmblarge, fpulay)


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
  ! Deallocate DIIS structures.
  call deallocateDIIS(ldiis)

  call deallocateCommunicationbufferSumrho(tmb%comsr, subname)
  call deallocateCommunicationbufferSumrho(tmbder%comsr, subname)

  call wait_p2p_communication(iproc, nproc, tmb%comgp)
  call deallocateCommunicationsBuffersPotential(tmb%comgp, subname)
  if(tmbder%wfnmd%bs%use_derivative_basis) then
     call wait_p2p_communication(iproc, nproc, tmbder%comgp)
     call deallocateCommunicationsBuffersPotential(tmbder%comgp, subname)
  end if
  iall=-product(shape(rhopotold_out))*kind(rhopotold_out)
  deallocate(rhopotold_out, stat=istat)
  call memocc(istat, iall, 'rhopotold_out', subname)

  if(input%lin%mixHist_highaccuracy>0) then
      call deallocateMixrhopotDIIS(mixdiis)
  end if

  !Write the linear wavefunctions to file if asked
  if(input%lin%plotBasisFunctions /= WF_FORMAT_NONE) then
    call writemywaves_linear(iproc,trim(input%dir_output) // 'minBasis',input%lin%plotBasisFunctions,tmb%Lzd,&
       tmbmix%orbs,orbs%norb,hx,hy,hz,at,rxyz,tmbmix%psi,tmbmix%wfnmd%coeff)
   end if
  ! Allocate the communication buffers for the calculation of the charge density.

  call allocateCommunicationbufferSumrho(iproc, tmbmix%comsr, subname)

  call communicate_basis_for_density(iproc, nproc, tmb%lzd, tmbmix%orbs, tmbmix%psi, tmbmix%comsr)
  call calculate_density_kernel(iproc, nproc, tmbmix%wfnmd%ld_coeff, orbs, tmbmix%orbs, &
       tmbmix%wfnmd%coeff, tmbmix%wfnmd%density_kernel)
  call sumrhoForLocalizedBasis2(iproc, nproc, tmb%lzd, input, hx, hy, hz, &
       tmbmix%orbs, tmbmix%comsr, tmbmix%wfnmd%density_kernel, Glr%d%n1i*Glr%d%n2i*denspot%dpbox%n3d, &
       denspot%rhov, at,denspot%dpbox%nscatterarr)

  call deallocateCommunicationbufferSumrho(tmbmix%comsr, subname)

  ! Build global orbitals psi (the physical ones).
  if(nproc>1) then
     allocate(psit(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
     call memocc(istat, psit, 'psit', subname)
  else
     psit => psi
  end if
  call transformToGlobal(iproc, nproc, tmb%lzd, tmbmix%orbs, orbs, comms, input, tmbmix%wfnmd%ld_coeff, &
       tmbmix%wfnmd%coeff, tmbmix%psi, psi, psit)
  if(nproc>1) then
     iall=-product(shape(psit))*kind(psit)
     deallocate(psit, stat=istat)
     call memocc(istat, iall, 'psit', subname)
  else
     nullify(psit)
  end if

  nullify(rho,pot)

  iall=-product(shape(lscv%locrad))*kind(lscv%locrad)
  deallocate(lscv%locrad, stat=istat)
  call memocc(istat, iall, 'lscv%locrad', subname)

  iall=-product(shape(eval))*kind(eval)
  deallocate(eval, stat=istat)
  call memocc(istat, iall, 'eval', subname)

  call timing(iproc,'WFN_OPT','PR')

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
      if(info_tmb<0) then
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
real(8),dimension(lorbs%npsidim_orbs),intent(inout):: lphi
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),target,intent(out):: psi
real(8),dimension(:),pointer,intent(out):: psit

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
  call transpose_v(iproc, nproc, lorbs, lzd%Glr%wfd, gcomms, phi, work=phiWork)


  if(iproc==0) then
      write(*,'(1x,a)', advance='no') '------------------------------------- Building linear combinations... '
  end if
  ! Build the extended orbital psi as a linear combination of localized basis functions phi. for real O(N)
  ! this has to replaced, but at the moment it is still needed.
  !call buildWavefunctionModified(iproc, nproc, orbs, gorbs, comms, gcomms, phi, psi, coeff)
  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  call dgemm('n', 'n', nvctrp, orbs%norb, lorbs%norb, 1.d0, phi(1), nvctrp, coeff(1,1), &
       lorbs%norb, 0.d0, psi(1), nvctrp)


  if(nproc>1) then
      call dcopy(orbs%npsidim_comp, psi, 1, psit, 1)
  else
      psit => psi
  end if

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
      lscv%nit_scc=input%lin%nitSCCWhenOptimizing_highaccuracy+input%lin%nitSCCWhenFixed_highaccuracy
      lscv%nit_scc_when_optimizing=input%lin%nitSCCWhenOptimizing_highaccuracy
      lscv%mix_hist=input%lin%mixHist_highaccuracy
      do ilr=1,nlr
          lscv%locrad(ilr)=input%lin%locrad_highaccuracy(ilr)
      end do
      if(wfnmd%bs%update_phi) then
          lscv%alpha_mix=input%lin%alphaMixWhenOptimizing_highaccuracy
      else
          lscv%alpha_mix=input%lin%alphaMixWhenFixed_highaccuracy
      end if
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
      lscv%nit_scc=input%lin%nitSCCWhenOptimizing_lowaccuracy+input%lin%nitSCCWhenFixed_lowaccuracy
      lscv%nit_scc_when_optimizing=input%lin%nitSCCWhenOptimizing_lowaccuracy
      lscv%mix_hist=input%lin%mixHist_lowaccuracy
      do ilr=1,nlr
          lscv%locrad(ilr)=input%lin%locrad_lowaccuracy(ilr)
      end do
      if(wfnmd%bs%update_phi) then
          lscv%alpha_mix=input%lin%alphaMixWhenOptimizing_lowaccuracy
      else
          lscv%alpha_mix=input%lin%alphaMixWhenFixed_lowaccuracy
      end if

  end if

end subroutine set_optimization_variables



subroutine adjust_locregs_and_confinement(iproc, nproc, hx, hy, hz, &
           input, tmb, tmbder, denspot, ldiis, lscv)
  use module_base
  use module_types
  use module_interfaces, except_this_one => adjust_locregs_and_confinement
  implicit none
  
  ! Calling argument
  integer,intent(in):: iproc, nproc
  real(8),intent(in):: hx, hy, hz
  type(input_variables),intent(in):: input
  type(DFT_wavefunction),intent(inout):: tmb, tmbder
  type(DFT_local_fields),intent(inout) :: denspot
  type(localizedDIISParameters),intent(inout):: ldiis
  type(linear_scaling_control_variables),intent(inout):: lscv

  ! Local variables
  integer :: ilr
  logical:: redefine_derivatives, change
  character(len=*),parameter:: subname='adjust_locregs_and_confinement'

  if(tmb%wfnmd%bs%confinement_decrease_mode==DECREASE_ABRUPT) then
      lscv%decrease_factor_total=1.d0
  else if(tmb%wfnmd%bs%confinement_decrease_mode==DECREASE_LINEAR) then
      if(lscv%info_basis_functions>0) then
          lscv%idecrease=lscv%idecrease+1
          lscv%ifail=0
      else
          lscv%ifail=lscv%ifail+1
      end if
      lscv%decrease_factor_total=1.d0-dble(lscv%idecrease)*input%lin%decrease_step
  end if
  if(tmbder%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) lscv%decrease_factor_total=1.d0
  if(iproc==0) write(*,'(1x,a,f6.2,a)') 'Changing the confining potential to ', &
      100.d0*lscv%decrease_factor_total,'% of its initial value.'
  tmb%confdatarr(:)%prefac=lscv%decrease_factor_total*tmb%confdatarr(:)%prefac


  lscv%locreg_increased=.false.
  redefine_derivatives=.false.
  !!if(lscv%ifail>=input%lin%increase_locrad_after .and. .not.lscv%lowaccur_converged .and. &
  !!   input%lin%locrad_increase_amount>0.d0) then
  !!    lscv%increase_locreg=lscv%increase_locreg+input%lin%locrad_increase_amount
  !!    !lscv%increase_locreg=lscv%increase_locreg+0.d0
  !!    if(iproc==0) then
  !!        write(*,'(1x,a)') 'It seems that the convergence criterion can not be reached with this localization radius.'
  !!        write(*,'(1x,a,f6.2)') 'The localization radius is increased by totally',lscv%increase_locreg
  !!    end if
  !!    lscv%ifail=0
  !!    lscv%locrad=lscv%locrad+lscv%increase_locreg
  !!    if(lscv%withder) then
  !!        redefine_derivatives=.true.
  !!    end if
  !!    lscv%locreg_increased=.true.
  !!end if

  !redefine_derivatives=.false.
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
         lscv%locreg_increased=.true.
      end if
      lscv%enlarge_locreg=.false. !flag to indicate that the locregs should not be increased any more in the following iterations
  end if
  if(lscv%locreg_increased) then
      call redefine_locregs_quantities(iproc, nproc, hx, hy, hz, lscv%locrad, .true., tmb%lzd, tmb, tmb, denspot, ldiis)
  end if
  if(redefine_derivatives) then
      call redefine_locregs_quantities(iproc, nproc, hx, hy, hz, tmb%lzd%llr(:)%locrad, .false., tmb%lzd, tmb, tmbder, denspot)
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
      !!lscv%nit_highaccuracy=lscv%nit_highaccuracy+1
      if(lscv%nit_highaccuracy==input%lin%nit_highaccuracy+1) then
          ! Deallocate DIIS structures.
          !!call deallocateDIIS(ldiis)
          !!lscv%exit_outer_loop=.true.
      end if
      !!! only use steepest descent if the localization regions may change
      !!if(input%lin%nItInnerLoop/=-1 .or. tmb%wfnmd%bs%locreg_enlargement/=1.d0) then
      !!    ldiis%isx=0
      !!end if
  
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


function check_whether_derivatives_to_be_used(input, itout, lscv)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(input_variables),intent(in):: input
  integer,intent(in):: itout
  type(linear_scaling_control_variables),intent(in):: lscv
  logical:: check_whether_derivatives_to_be_used

  ! Local variables
  logical:: withder

  if(input%lin%mixedmode) then
      if( (.not.lscv%lowaccur_converged .and. &
           (itout==input%lin%nit_lowaccuracy+1 .or. lscv%pnrm_out<input%lin%lowaccuray_converged) ) &
          .or. lscv%lowaccur_converged ) then
          withder=.true.
      else
          withder=.false.
      end if
  else
      if(input%lin%useDerivativeBasisFunctions) then
          withder=.true.
      else
          withder=.false.
      end if
  end if
  check_whether_derivatives_to_be_used=withder

end function check_whether_derivatives_to_be_used



subroutine check_whether_lowaccuracy_converged(itout, input, lscv)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: itout
  type(input_variables),intent(in):: input
  type(linear_scaling_control_variables),intent(inout):: lscv
  
  if(.not.lscv%lowaccur_converged .and. &
     (itout==input%lin%nit_lowaccuracy+1 .or. lscv%pnrm_out<input%lin%lowaccuray_converged .or. &
      lscv%decrease_factor_total<1.d0-input%lin%decrease_amount)) then
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
  real(8),dimension(3,at%nat),intent(in):: rxyz
  type(nonlocal_psp_descriptors),intent(in):: nlpspd
  real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
  type(SIC_data),intent(in):: SIC
  type(DFT_local_fields), intent(inout) :: denspot
  type(GPU_pointers),intent(inout):: GPU
  type(DFT_wavefunction),intent(in):: tmb
  type(DFT_wavefunction),intent(inout):: tmblarge
  real(8),dimension(3,at%nat),intent(out):: fpulay

  ! Local variables
  integer:: norb, norbu, norbd, npsin, istat, iall, nspin, tag, ierr
  integer:: i, ist_orig, ist_dest, iorb, iiorb, ilr, ncount, nlr, ii, ityp
  integer:: jjorb, jat, jorbsmall, kkorb, kat, korbsmall, jdir, kdir, iat, npsidim, ndim
  type(DFT_wavefunction):: tmbder
  real(8),dimension(:),allocatable:: phiLoc, lhphi, lhphilarge, psit_c, psit_f, hpsit_c, hpsit_f
  real(8),dimension(:,:),allocatable:: matrix, locregCenter
  type(energy_terms) :: energs
  type(confpot_data),dimension(:),allocatable :: confdatarrtmp
  integer,dimension(:),allocatable:: norbsPerAtom, norbsPerLocreg
  logical:: use_derivative_basis
  character(len=*),parameter:: subname='pulay_correction'


  call nullify_orbitals_data(tmbder%orbs)
  call nullify_p2pComms(tmbder%comrp)
  call nullify_collective_comms(tmbder%collcom)
  call nullify_matrixDescriptors(tmbder%mad)
  call nullify_overlapParameters(tmbder%op)
  call nullify_p2pComms(tmbder%comon)

  norbu=4*tmb%orbs%norb
  norb=norbu
  norbd=0
  nspin=1
  tag=0

  allocate(norbsPerAtom(at%nat), stat=istat)
  call memocc(istat, norbsPerAtom, 'norbsPerAtom', subname)
  norb=0
  nlr=0
  use_derivative_basis=.true.
  if(use_derivative_basis) then
      ii=4
  else
      ii=1
  end if
  do iat=1,at%nat
      ityp=at%iatype(iat)
      norbsPerAtom(iat)=input%lin%norbsPerType(ityp)
      norb=norb+ii*input%lin%norbsPerType(ityp)
      nlr=nlr+input%lin%norbsPerType(ityp)
  end do


  allocate(norbsPerLocreg(nlr), stat=istat) 
  call memocc(istat, norbsPerLocreg, 'norbsPerLocreg', subname) 
  norbsPerLocreg=ii !should be norbsPerLocreg 
     


  call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, nspin, tmblarge%orbs%nspinor,&
       tmblarge%orbs%nkpts, tmblarge%orbs%kpts, tmblarge%orbs%kwgts, tmbder%orbs,.true.) !simple repartition

  iall=-product(shape(tmbder%orbs%onwhichatom))*kind(tmbder%orbs%inWhichLocreg) 
  deallocate(tmbder%orbs%onwhichatom, stat=istat) 
  call memocc(istat, iall, 'lorbs%onwhichatom', subname) 
  call assignToLocreg2(iproc, nproc, tmbder%orbs%norb, tmbder%orbs%norb_par, at%nat, at%nat, &
       nspin, norbsPerAtom, rxyz, tmbder%orbs%onwhichatom)

  iall=-product(shape(tmbder%orbs%inWhichLocreg))*kind(tmbder%orbs%inWhichLocreg) 
  deallocate(tmbder%orbs%inWhichLocreg, stat=istat) 
  call memocc(istat, iall, 'lorbs%inWhichLocreg', subname) 


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


  call assignToLocreg2(iproc, nproc, tmbder%orbs%norb, tmbder%orbs%norb_par, at%nat, nlr, &
       nspin, norbsPerLocreg, locregCenter, tmbder%orbs%inwhichlocreg)

  npsidim = 0
  do iorb=1,tmbder%orbs%norbp
      iiorb=tmbder%orbs%isorb+iorb
      ilr=tmbder%orbs%inwhichlocreg(iiorb)
      npsidim = npsidim + tmblarge%lzd%llr(ilr)%wfd%nvctr_c+7*tmblarge%lzd%llr(ilr)%wfd%nvctr_f
  end do
  tmbder%orbs%npsidim_orbs=npsidim

  !!write(*,'(a,i5,3i10)') 'iproc, tmbder%orbs%norb, tmbder%orbs%norbp, tmbder%orbs%npsidim_orbs',&
  !!    iproc, tmbder%orbs%norb, tmbder%orbs%norbp, tmbder%orbs%npsidim_orbs
  !!write(*,'(a,i6,100i5)') 'iproc, iwl', iproc, tmbder%orbs%inwhichlocreg

  call initializeRepartitionOrbitals(iproc, nproc, tag, tmblarge%orbs, tmbder%orbs, tmblarge%lzd, tmbder%comrp)
  call init_collective_comms(iproc, nproc, tmbder%orbs, tmblarge%lzd, tmbder%collcom)

  call initCommsOrtho(iproc, nproc, nspin, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
       tmblarge%lzd, tmblarge%lzd, tmbder%orbs, 's', tmb%wfnmd%bpo, tmbder%op, tmbder%comon)

  ndim = maxval(tmbder%op%noverlaps)
  call initMatrixCompression(iproc, nproc, tmblarge%lzd%nlr, ndim, tmbder%orbs, &
       tmbder%op%noverlaps, tmbder%op%overlaps, tmbder%mad)
  call initCompressedMatmul3(iproc, tmbder%orbs%norb, tmbder%mad)


  allocate(tmbder%psi(tmbder%orbs%npsidim_orbs), stat=istat)
  call memocc(istat, tmbder%psi, 'tmbder%psi', subname)

  if (tmblarge%orbs%npsidim_orbs > 0) call to_zero(tmblarge%orbs%npsidim_orbs,tmblarge%psi(1))
  call small_to_large_locreg(iproc, nproc, tmb%lzd, tmblarge%lzd, tmb%orbs, tmblarge%orbs, tmb%psi, tmblarge%psi)

  call getDerivativeBasisFunctions(iproc,nproc,tmblarge%lzd%hgrids(1),tmblarge%lzd,tmblarge%orbs,tmbder%orbs,tmbder%comrp,&
       max(tmblarge%orbs%npsidim_orbs,tmblarge%orbs%npsidim_comp),tmblarge%psi,tmbder%psi)


  allocate(phiLoc(4*max(tmblarge%orbs%npsidim_orbs,tmblarge%orbs%npsidim_comp)), stat=istat)
  call memocc(istat, phiLoc, 'phiLoc', subname)


  ! Apply Hamiltonian to tmb%psi
  allocate(lhphi(tmb%orbs%npsidim_orbs), stat=istat)
  call memocc(istat, lhphi, 'lhphi', subname)

  if (tmb%orbs%npsidim_orbs > 0) call to_zero(tmb%orbs%npsidim_orbs,lhphi(1))
  call local_potential_dimensions(tmblarge%lzd,tmblarge%orbs,denspot%dpbox%ngatherarr(0,1))

  allocate(lhphilarge(tmblarge%orbs%npsidim_orbs), stat=istat)
  call memocc(istat, lhphilarge, 'lhphilarge', subname)

  if (tmblarge%orbs%npsidim_orbs > 0) call to_zero(tmblarge%orbs%npsidim_orbs,lhphilarge(1))

  call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
       tmblarge%comgp%nrecvbuf, tmblarge%comgp%recvbuf, tmblarge%comgp)

  allocate(confdatarrtmp(tmbder%orbs%norbp))
  call default_confinement_data(confdatarrtmp,tmbder%orbs%norbp)

  allocate(tmblarge%lzd%doHamAppl(tmblarge%lzd%nlr), stat=istat)
  call memocc(istat, tmblarge%lzd%doHamAppl, 'tmblarge%lzd%doHamAppl', subname)
  tmblarge%lzd%doHamAppl=.true.

  call NonLocalHamiltonianApplication(iproc,at,tmblarge%orbs,rxyz,&
       proj,tmblarge%lzd,nlpspd,tmblarge%psi,lhphilarge,energs%eproj)

  call LocalHamiltonianApplication(iproc,nproc,at,tmblarge%orbs,&
       tmblarge%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmblarge%psi,lhphilarge,&
       energs,SIC,GPU,.false.,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmblarge%comgp)
  call timing(iproc,'glsynchham1','ON') !lr408t
  call SynchronizeHamiltonianApplication(nproc,tmblarge%orbs,tmblarge%lzd,GPU,lhphilarge,&
       energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
  call timing(iproc,'glsynchham1','OF') !lr408t
  deallocate(confdatarrtmp)


  ! Copy tmblarge%psi to phiLoc
  ist_orig=1
  ist_dest=1
  do iorb=1,tmblarge%orbs%norbp
      iiorb=tmblarge%orbs%isorb+iorb
      ilr=tmblarge%orbs%inwhichlocreg(iiorb)
      ncount=tmblarge%lzd%llr(ilr)%wfd%nvctr_c+7*tmblarge%lzd%llr(ilr)%wfd%nvctr_f
      do i=1,4
          call dcopy(ncount, tmblarge%psi(ist_orig), 1, phiLoc(ist_dest), 1)
          ist_dest=ist_dest+ncount
      end do
      ist_orig=ist_orig+ncount
  end do

  ! Repartition phiLoc to tmbder%psi
  call post_p2p_communication(iproc, nproc, size(phiLoc), phiLoc, size(tmbder%psi), tmbder%psi, tmbder%comrp)
  call wait_p2p_communication(iproc, nproc, tmbder%comrp)


  ! Calculate matrix

  allocate(matrix(tmbder%orbs%norb,tmbder%orbs%norb), stat=istat)
  call memocc(istat, matrix, 'matrix', subname)

  allocate(hpsit_c(tmbder%collcom%ndimind_c))
  call memocc(istat, hpsit_c, 'hpsit_c', subname)
  allocate(hpsit_f(7*tmbder%collcom%ndimind_f))
  call memocc(istat, hpsit_f, 'hpsit_f', subname)
  allocate(psit_c(tmbder%collcom%ndimind_c))
  call memocc(istat, psit_c, 'psit_c', subname)
  allocate(psit_f(7*tmbder%collcom%ndimind_f))
  call memocc(istat, psit_f, 'psit_f', subname)

  call transpose_localized(iproc, nproc, tmbder%orbs,  tmbder%collcom, &
       tmbder%psi, psit_c, psit_f, tmblarge%lzd)
  call transpose_localized(iproc, nproc, tmbder%orbs,  tmbder%collcom, &
       lhphilarge, hpsit_c, hpsit_f, tmblarge%lzd)
  call calculate_overlap_transposed(iproc, nproc, tmbder%orbs, tmbder%mad, tmbder%collcom, &
       psit_c, hpsit_c, psit_f, hpsit_f, matrix)

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
  call to_zero(3*at%nat, fpulay(1,1))
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      do jjorb=1,tmbder%orbs%norb
          jat=tmbder%orbs%onwhichatom(jjorb)
          jdir=mod(jjorb-1,4) ! get direction x, y or z (0=no direction)
          jorbsmall=ceiling(dble(jjorb)/4.d0)
          do kkorb=1,tmbder%orbs%norb
              kat=tmbder%orbs%onwhichatom(kkorb)
              kdir=mod(kkorb-1,4) ! get direction x, y or z (0=no direction)
              korbsmall=ceiling(dble(kkorb)/4.d0)
              if(jdir>0) then
                  fpulay(jdir,jat) = fpulay(jdir,jat) + &
                      tmb%wfnmd%coeff(jorbsmall,iiorb)*tmb%wfnmd%coeff(korbsmall,iiorb)*matrix(jjorb,kkorb)
              end if
              if(kdir>0) then
                  fpulay(kdir,kat) = fpulay(kdir,kat) + &
                      tmb%wfnmd%coeff(jorbsmall,iiorb)*tmb%wfnmd%coeff(korbsmall,iiorb)*matrix(kkorb,jjorb)
              end if
          end do
      end do
  end do

  call mpiallred(fpulay(1,1), 3*at%nat, mpi_sum, mpi_comm_world, ierr)
  if(iproc==0) then
       do iat=1,at%nat
           write(*,'(a,i5,3es16.6)') 'iat, fpulay', iat, fpulay(1:3,iat)
       end do
  end if


  iall=-product(shape(matrix))*kind(matrix)
  deallocate(matrix, stat=istat)
  call memocc(istat, iall, 'matrix', subname)

  iall=-product(shape(lhphilarge))*kind(lhphilarge)
  deallocate(lhphilarge, stat=istat)
  call memocc(istat, iall, 'lhphilarge', subname)

  iall=-product(shape(lhphi))*kind(lhphi)
  deallocate(lhphi, stat=istat)
  call memocc(istat, iall, 'lhphi', subname)

  iall=-product(shape(tmblarge%lzd%doHamAppl))*kind(tmblarge%lzd%doHamAppl)
  deallocate(tmblarge%lzd%doHamAppl, stat=istat)
  call memocc(istat, iall, 'tmblarge%lzd%doHamAppl', subname)

  iall=-product(shape(denspot%pot_work))*kind(denspot%pot_work)
  deallocate(denspot%pot_work, stat=istat)
  call memocc(istat, iall, 'denspot%pot_work', subname)

  if(iproc==0) write(*,'(1x,a)') 'done.'




  iall=-product(shape(locregCenter))*kind(locregCenter)
  deallocate(locregCenter, stat=istat)
  call memocc(istat, iall, 'locregCenter', subname)

  iall=-product(shape(norbsPerAtom))*kind(norbsPerAtom)
  deallocate(norbsPerAtom, stat=istat)
  call memocc(istat, iall, 'norbsPerAtom', subname)

  iall=-product(shape(norbsPerLocreg))*kind(norbsPerLocreg)
  deallocate(norbsPerLocreg, stat=istat)
  call memocc(istat, iall, 'norbsPerLocreg', subname)

  iall=-product(shape(phiLoc))*kind(phiLoc)
  deallocate(phiLoc, stat=istat)
  call memocc(istat, iall, 'phiLoc', subname)

  iall=-product(shape(tmbder%psi))*kind(tmbder%psi)
  deallocate(tmbder%psi, stat=istat)
  call memocc(istat, iall, 'tmbder%psi', subname)

  call deallocate_orbitals_data(tmbder%orbs, subname)
  call deallocate_p2pComms(tmbder%comrp, subname)

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
