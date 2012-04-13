subroutine linearScaling(iproc,nproc,Glr,orbs,comms,at,input,hx,hy,hz,&
     rxyz,fion,fdisp,denspot,nlpspd,proj,GPU,&
     eion,edisp,eexctX,scpot,psi,psit,energy)
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
type(nonlocal_psp_descriptors),intent(in):: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
type(GPU_pointers),intent(in out):: GPU
real(gp),intent(in):: eion, edisp, eexctX,hx,hy,hz
logical,intent(in):: scpot
real(8),dimension(:),pointer,intent(out):: psi, psit
real(gp), dimension(:), pointer :: rho,pot
real(8),intent(out):: energy

type(linear_scaling_control_variables):: lscv
integer:: infoCoeff,istat,iall,it_scc,ilr,tag,itout,iorb,ist,iiorb,ncnt
real(8):: ebs,pnrm,ehart,eexcu,vexcu,trace
character(len=*),parameter:: subname='linearScaling'
real(8),dimension(:),allocatable:: rhopotOld, rhopotold_out
real(8):: energyold, energyDiff, energyoldout
type(mixrhopotDIISParameters):: mixdiis
type(localizedDIISParameters):: ldiis
type(DFT_wavefunction),target:: tmb
type(DFT_wavefunction),target:: tmbder
type(DFT_wavefunction),pointer:: tmbmix
logical:: check_whether_derivatives_to_be_used
real(8),dimension(:),allocatable:: psit_c, psit_f, philarge, lphiovrlp, psittemp_c, psittemp_f
real(8),dimension(:,:),allocatable:: ovrlp, philarge_root
integer:: jorb, ldim, sdim, ists, istl, nspin, ierr
real(8):: ddot, tt1, tt2, tt3


  if(iproc==0) then
      write(*,'(1x,a)') repeat('*',84)
      write(*,'(1x,a)') '****************************** LINEAR SCALING VERSION ******************************'
  end if



  ! Initialize everything related to the linear scaling version ###########################################################
  call lin_input_variables_new(iproc,trim(input%file_lin),input,at)

  tmbder%wfnmd%bs%use_derivative_basis=input%lin%useDerivativeBasisFunctions
  tmb%wfnmd%bs%use_derivative_basis=.false.

  call init_orbitals_data_for_linear(iproc, nproc, orbs%nspinor, input, at, glr, tmb%wfnmd%bs%use_derivative_basis, rxyz, &
       tmb%orbs)
  call orbitals_communicators(iproc, nproc, glr, tmb%orbs, tmb%comms)

  call init_orbitals_data_for_linear(iproc, nproc, orbs%nspinor, input, at, glr, tmbder%wfnmd%bs%use_derivative_basis, rxyz, &
       tmbder%orbs)
  call orbitals_communicators(iproc, nproc, glr, tmbder%orbs, tmbder%comms)

  if(iproc==0) call print_orbital_distribution(iproc, nproc, tmb%orbs, tmbder%orbs)

  call init_local_zone_descriptors(iproc, nproc, input, hx, hy, hz, glr, at, rxyz, tmb%orbs, tmbder%orbs, tmb%lzd)

  call update_wavefunctions_size(tmb%lzd,tmb%orbs)
  call update_wavefunctions_size(tmb%lzd,tmbder%orbs)

  call create_wfn_metadata('l', max(tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp), tmb%orbs%norb, &
       tmb%orbs%norb, orbs%norb, input, tmb%wfnmd)
  allocate(tmb%psi(tmb%wfnmd%nphi), stat=istat)
  call memocc(istat, tmb%psi, 'tmb%psi', subname)

  call create_wfn_metadata('l', max(tmbder%orbs%npsidim_orbs,tmbder%orbs%npsidim_comp), tmbder%orbs%norb, &
       tmbder%orbs%norb, orbs%norb, input, tmbder%wfnmd)
  allocate(tmbder%psi(tmbder%wfnmd%nphi), stat=istat)
  call memocc(istat, tmbder%psi, 'tmbder%psi', subname)

  tmbder%wfnmd%bs%use_derivative_basis=input%lin%useDerivativeBasisFunctions
  tmb%wfnmd%bs%use_derivative_basis=.false.

  tag=0
  call initCommsOrtho(iproc, nproc, input%nspin, hx, hy, hz, tmb%lzd, tmb%lzd, &
       tmb%orbs,  tmb%orbs, tmb%orbs%inWhichLocreg,&
       input%lin%locregShape, tmb%op, tmb%comon, tag)
  call initCommsOrtho(iproc, nproc, input%nspin, hx, hy, hz, tmb%lzd, tmb%lzd, &
       tmbder%orbs, tmbder%orbs, tmbder%orbs%inWhichLocreg, &
       input%lin%locregShape, tmbder%op, tmbder%comon, tag)
  
  call initializeCommunicationPotential(iproc, nproc, denspot%dpcom%nscatterarr, &
       tmb%orbs, tmb%lzd, tmb%comgp, tmb%orbs%inWhichLocreg, tag)
  call initializeCommunicationPotential(iproc, nproc, denspot%dpcom%nscatterarr, &
       tmbder%orbs, tmb%lzd, tmbder%comgp, tmbder%orbs%inWhichLocreg, tag)

  if(input%lin%useDerivativeBasisFunctions) then
      call initializeRepartitionOrbitals(iproc, nproc, tag, tmb%orbs, tmbder%orbs, tmb%lzd, tmbder%comrp)
      call initializeRepartitionOrbitals(iproc, nproc, tag, tmb%orbs, tmbder%orbs, tmb%lzd, tmb%comrp)
  else
      call nullify_p2pComms(tmbder%comrp)
      call nullify_p2pComms(tmb%comrp)
  end if


  call nullify_p2pcomms(tmb%comsr)
  call initializeCommsSumrho(iproc, nproc, denspot%dpcom%nscatterarr, tmb%lzd, tmb%orbs, tag, tmb%comsr)
  call nullify_p2pcomms(tmbder%comsr)
  call initializeCommsSumrho(iproc, nproc, denspot%dpcom%nscatterarr, tmb%lzd, tmbder%orbs, tag, tmbder%comsr)

  call initMatrixCompression(iproc, nproc, tmb%lzd%nlr, tmb%orbs, tmb%op%noverlaps, tmb%op%overlaps, tmb%mad)
  call initCompressedMatmul3(tmb%orbs%norb, tmb%mad)
  call initMatrixCompression(iproc, nproc, tmb%lzd%nlr, tmbder%orbs, &
       tmbder%op%noverlaps, tmbder%op%overlaps, tmbder%mad)
  call initCompressedMatmul3(tmbder%orbs%norb, tmbder%mad)

  allocate(tmb%confdatarr(tmb%orbs%norbp))
  call define_confinement_data(tmb%confdatarr,tmb%orbs,rxyz,at,&
       input%hx,input%hy,input%hz,input%lin%confpotorder,input%lin%potentialprefac_lowaccuracy,tmb%lzd,tmb%orbs%onwhichatom)

  allocate(tmbder%confdatarr(tmbder%orbs%norbp))
  call define_confinement_data(tmbder%confdatarr,tmbder%orbs,rxyz,at,&
       input%hx,input%hy,input%hz,input%lin%confpotorder,&
       input%lin%potentialprefac_lowaccuracy,tmb%lzd,tmbder%orbs%onwhichatom)

  call nullify_collective_comms(tmb%collcom)
  call nullify_collective_comms(tmbder%collcom)
  call init_collective_comms(iproc, nproc, tmb%orbs, tmb%lzd, tmb%collcom)
  call init_collective_comms(iproc, nproc, tmbder%orbs, tmb%lzd, tmbder%collcom)

  ! Now all initializations are done ######################################################################################


  ! Assign some values to orthpar
  tmb%orthpar%methTransformOverlap = tmb%wfnmd%bs%meth_transform_overlap
  tmb%orthpar%nItOrtho = input%lin%nItOrtho
  tmb%orthpar%blocksize_pdsyev = tmb%wfnmd%bpo%blocksize_pdsyev
  tmb%orthpar%blocksize_pdgemm = tmb%wfnmd%bpo%blocksize_pdgemm

  tmbder%orthpar%methTransformOverlap = tmb%wfnmd%bs%meth_transform_overlap
  tmbder%orthpar%nItOrtho = input%lin%nItOrtho
  tmbder%orthpar%blocksize_pdsyev = tmb%wfnmd%bpo%blocksize_pdsyev
  tmbder%orthpar%blocksize_pdgemm = tmb%wfnmd%bpo%blocksize_pdgemm


  ! Allocate the global orbitals psi and psit
  if(.not.input%lin%transformToGlobal) then
      ! psi and psit will not be calculated, so only allocate them with size 1
      orbs%npsidim_orbs=1
      orbs%npsidim_comp=1
  end if
  allocate(psi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
  call memocc(istat, psi, 'psi', subname)
  if(nproc>1) then
      allocate(psit(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
      call memocc(istat, psit, 'psit', subname)
  else
      psit => psi
  end if

  ! Allocate the old charge density (used to calculate the variation in the charge density)
  allocate(rhopotold(max(glr%d%n1i*glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin), stat=istat)
  call memocc(istat, rhopotold, 'rhopotold', subname)
  allocate(rhopotold_out(max(glr%d%n1i*glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin), stat=istat)
  call memocc(istat, rhopotold_out, 'rhopotold_out', subname)




  ! Generate the input guess for the TMB
  tmb%wfnmd%bs%update_phi=.false.
  call inputguessConfinement(iproc, nproc, at, input, hx, hy, hz, tmb%lzd, tmb%orbs, rxyz, denspot ,rhopotold, &
       nlpspd, proj, GPU,  tmb%psi, orbs, tmb)

  !! Now one could calculate the charge density like this. It is not done since we would in this way overwrite
  !! the potential from the input guess.     
  !call allocateCommunicationbufferSumrho(iproc, with_auxarray, tmb%comsr, subname)
  !call communicate_basis_for_density(iproc, nproc, tmb%lzd, tmb%orbs, tmb%psi, tmb%comsr)
  !call sumrhoForLocalizedBasis2(iproc, nproc, orbs%norb, tmb%lzd, input, hx, hy, hz, tmb%orbs, tmb%comsr, &
  !     tmb%wfnmd%ld_coeff, tmb%wfnmd%coeff, Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3d, denspot%rhov, at,denspot%dpcom%nscatterarr)
  !call deallocateCommunicationbufferSumrho(tmb%comsr, subname)


  ! Initialize the DIIS mixing of the potential if required.
  if(input%lin%mixHist_lowaccuracy>0) then
      call initializeMixrhopotDIIS(input%lin%mixHist_lowaccuracy, denspot%dpcom%ndimpot, mixdiis)
  end if

  !end of the initialization part, will later be moved to cluster
  call timing(iproc,'INIT','PR')

  allocate(lscv%locrad(tmb%lzd%nlr), stat=istat)
  call memocc(istat, lscv%locrad, 'lscv%locrad', subname)


  if(input%lin%nItInguess>0) then
      tmb%wfnmd%bs%communicate_phi_for_lsumrho=.true.
      tmb%wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE

      do ilr=1,tmb%lzd%nlr
          lscv%locrad(ilr)=input%lin%locrad_lowaccuracy(ilr)
      end do

      if(trim(input%lin%mixingMethod)=='dens') then
          rhopotold_out=rhopotold
      end if

      if(trim(input%lin%mixingMethod)=='pot') then
          rhopotold_out=denspot%rhov
      end if

      ! Copy the current potential
      if(trim(input%lin%mixingMethod)=='pot') then
           call dcopy(max(Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotOld(1), 1)
      end if
  end if


  ! Allocate the communications buffers needed for the communications of teh potential and
  ! post the messages. This will send to each process the part of the potential that this process
  ! needs for the application of the Hamlitonian to all orbitals on that process.
  call allocateCommunicationsBuffersPotential(tmb%comgp, subname)
  !!call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, tmb%comgp)
  call post_p2p_communication(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, &
       tmb%comgp%nrecvbuf, tmb%comgp%recvbuf, tmb%comgp)
  ! If we also use the derivative of the basis functions, also send the potential in this case. This is
  ! needed since the orbitals may be partitioned in a different way when the derivatives are used.
  if(tmbder%wfnmd%bs%use_derivative_basis .and. .not.input%lin%mixedMode) then
      call allocateCommunicationsBuffersPotential(tmbder%comgp, subname)
      !!call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, tmbder%comgp)
      call post_p2p_communication(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, &
           tmbder%comgp%nrecvbuf, tmbder%comgp%recvbuf, tmbder%comgp)
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

  ! tmbmix is the types we use for the mixing. It will point to either tmb if we don't use the derivatives
  ! or to tmbder if we use the derivatives.
  if(input%lin%useDerivativeBasisFunctions) then
      tmbmix => tmbder
  else
      tmbmix => tmb
  end if

  ! Check whether it is possible to have variable localization regions or not.
  if(tmb%wfnmd%bs%nit_unitary_loop==-1 .and. tmb%wfnmd%bs%locreg_enlargement==1.d0) then
      lscv%variable_locregs=.false.
  else
      lscv%variable_locregs=.true.
  end if

  ! This is the main outer loop. Each iteration of this loop consists of a first loop in which the basis functions
  ! are optimized and a consecutive loop in which the density is mixed.
  outerLoop: do itout=1,input%lin%nit_lowaccuracy+input%lin%nit_highaccuracy

      ! First to some initialization and determine the value of some control parameters.

      ! Initialize DIIS...
      call initializeDIIS(input%lin%DIISHistMax, tmb%lzd, tmb%orbs, tmb%orbs%norb, ldiis)
      ldiis%DIISHistMin=input%lin%DIISHistMin
      ldiis%DIISHistMax=input%lin%DIISHistMax
      ldiis%alphaSD=input%lin%alphaSD
      ldiis%alphaDIIS=input%lin%alphaDIIS

      ! The basis functions shall be optimized
      tmb%wfnmd%bs%update_phi=.true.

      ! Convergence criterion for the self consistency loop
      lscv%self_consistent=input%lin%convCritMix


      ! Check whether the low accuracy part (i.e. with strong confining potential) has converged.
      call check_whether_lowaccuracy_converged(itout, input, lscv)

      ! Check whether the derivatives shall be used or not.
      lscv%withder=check_whether_derivatives_to_be_used(input, itout, lscv)


      ! Set all remaining variables that we need for the optimizations of the basis functions and the mixing.
      call set_optimization_variables(input, at, tmb%orbs, tmb%lzd%nlr, tmb%orbs%onwhichatom, &
           tmb%confdatarr, tmb%wfnmd, lscv)
      call set_optimization_variables(input, at, tmbder%orbs, tmb%lzd%nlr, tmbder%orbs%onwhichatom, &
           tmbder%confdatarr, tmbder%wfnmd, lscv)


      ! Adjust the confining potential if required.
      call adjust_locregs_and_confinement(iproc, nproc, hx, hy, hz, &
           input, tmb, tmbder, denspot, ldiis, lscv)

      ! Somce special treatement if we are in the high accuracy part
      call adjust_DIIS_for_high_accuracy(input, tmb, denspot, ldiis, mixdiis, lscv)
      if(lscv%exit_outer_loop) exit outerLoop

      !!! Allocate the communication arrays for the calculation of the charge density.
      !!if(.not. lscv%locreg_increased) call allocateCommunicationbufferSumrho(iproc, tmb%comsr, subname)
      !!call allocateCommunicationbufferSumrho(iproc, tmbder%comsr, subname)
      if(lscv%locreg_increased) call deallocateCommunicationbufferSumrho(tmb%comsr, subname)

      ! Now all initializations are done...


      ! The self consistency cycle. Here we try to get a self consistent density/potential.
      ! In the first lscv%nit_scc_when_optimizing iteration, the basis functions are optimized, whereas in the remaining
      ! iteration the basis functions are fixed.
      do it_scc=1,lscv%nit_scc
          ! Do not update the TMB if it_scc>1
          if(it_scc>lscv%nit_scc_when_optimizing) tmb%wfnmd%bs%update_phi=.false.

          !!if(lscv%locreg_increased) call deallocateCommunicationbufferSumrho(tmb%comsr, subname)

          ! Update the basis functions (if wfnmd%bs%update_phi is true), calculate the Hamiltonian in this basis, and diagonalize it.
          ! This is a flag whether the basis functions shall be updated.
          if(tmb%wfnmd%bs%update_phi) then
              ! Improve the trace minimizing orbitals.
              if(itout>1 .and. tmbmix%wfnmd%bs%use_derivative_basis) then
                  do iorb=1,orbs%norb
                      call dcopy(tmb%orbs%norb, tmbmix%wfnmd%coeff_proj(1,iorb), 1, tmb%wfnmd%coeff(1,iorb), 1)
                  end do
              end if
              call getLocalizedBasis(iproc,nproc,at,orbs,rxyz,denspot,GPU,trace, lscv%info_basis_functions,&
                  nlpspd,proj,ldiis,input%SIC,lscv%locrad,tmb)
              tmb%wfnmd%nphi=tmb%orbs%npsidim_orbs
          end if

          ! Decide whether we have to use the derivatives or not.
          if(input%lin%mixedmode) then
              if(.not.lscv%withder) then
                  tmbmix => tmb
              else
                  ! We have to communicate the potential in the first iteration
                  if(it_scc==1) then
                      call allocateCommunicationsBuffersPotential(tmbder%comgp, subname)
                      !!call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, tmbder%comgp)
                      call post_p2p_communication(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, &
                           tmbder%comgp%nrecvbuf, tmbder%comgp%recvbuf, tmbder%comgp)
                  end if
                  tmbmix => tmbder
              end if
          end if
          if(tmbmix%wfnmd%bs%use_derivative_basis) then
              ! Cancel the communication of the potential for the TMB, since we need in the following
              ! only the potential for the TMB including the derivatives.
              !call cancelCommunicationPotential(iproc, nproc, tmb%comgp)
              !!if(.not.tmb%wfnmd%bs%update_phi) call gatherPotential(iproc, nproc, tmb%comgp)
              if(.not.tmb%wfnmd%bs%update_phi) call wait_p2p_communication(iproc, nproc, tmb%comgp)
              call deallocateCommunicationsBuffersPotential(tmb%comgp, subname)
          end if
          if((lscv%locreg_increased .or. (lscv%variable_locregs .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY)) &
              .and. tmb%wfnmd%bs%update_phi) then
              ! Redefine some quantities if the localization region has changed.
              if(lscv%withder) then
                  call redefine_locregs_quantities(iproc, nproc, hx, hy, hz, tmb%lzd, tmb, tmbder, denspot)
                  call deallocateCommunicationbufferSumrho(tmbder%comsr, subname)
                  tmbmix => tmbder
              else
                  call redefine_locregs_quantities(iproc, nproc, hx, hy, hz, tmb%lzd, tmb, tmb, denspot)
                  call deallocateCommunicationbufferSumrho(tmb%comsr, subname)
                  tmbmix => tmb
              end if
          end if


          ! Build the derivatives if required.
          if(tmb%wfnmd%bs%update_phi .or. it_scc==0) then
              if(tmbmix%wfnmd%bs%use_derivative_basis) then
                  if((lscv%locreg_increased .or. &
                      (lscv%variable_locregs .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY)) &
                      .and. tmb%wfnmd%bs%update_phi) then
                      call deallocate_p2pComms(tmbder%comrp, subname)
                      call nullify_p2pComms(tmbder%comrp)
                      call initializeRepartitionOrbitals(iproc, nproc, tag, tmb%orbs, tmbder%orbs, tmb%lzd, tmbder%comrp)
                      tmbmix => tmbder
                  end if
                  if(iproc==0) write(*,'(1x,a)',advance='no') 'calculating derivative basis functions...'
                  call getDerivativeBasisFunctions(iproc,nproc,hx,tmb%lzd,tmb%orbs,tmbmix%orbs,tmbmix%comrp,&
                       max(tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp),tmb%psi,tmbmix%psi)
                  if(iproc==0) write(*,'(a)') 'done.'
              else
                  call dcopy(tmb%wfnmd%nphi, tmb%psi(1), 1, tmbmix%psi(1), 1)
              end if
          end if

          ! Only communicate the TMB for sumrho in the first iteration.
          if(it_scc<=lscv%nit_scc_when_optimizing) then
              tmbmix%wfnmd%bs%communicate_phi_for_lsumrho=.true.
              call allocateCommunicationbufferSumrho(iproc, tmb%comsr, subname)
              call allocateCommunicationbufferSumrho(iproc, tmbder%comsr, subname)
          !else if(lscv%locreg_increased) then
          !    call allocateCommunicationbufferSumrho(iproc, tmb%comsr, subname)
          !    tmbmix%wfnmd%bs%communicate_phi_for_lsumrho=.false.
          else
              tmbmix%wfnmd%bs%communicate_phi_for_lsumrho=.false.
          end if
          !!if(lscv%locreg_increased) then
          !!    call allocateCommunicationbufferSumrho(iproc, tmb%comsr, subname)
          !!end if

          ! Calculate the coefficients
          call mpi_barrier(mpi_comm_world, istat)
          call get_coeff(iproc,nproc,tmb%lzd,orbs,at,rxyz,denspot,GPU,infoCoeff,ebs,nlpspd,proj,&
               tmbmix%wfnmd%bpo%blocksize_pdsyev,tmbder%wfnmd%bpo%nproc_pdsyev,&
               hx,hy,hz,input%SIC,tmbmix)
          ! Deallocate the buffers needed for the communication of the potential.
          if(lscv%withder) then
              call deallocateCommunicationsBuffersPotential(tmbder%comgp, subname)
          else
              call deallocateCommunicationsBuffersPotential(tmb%comgp, subname)
          end if


          ! Calculate the charge density.
          call sumrhoForLocalizedBasis2(iproc, nproc, orbs%norb,&
               tmb%lzd, input, hx, hy ,hz, tmbmix%orbs, tmbmix%comsr, &
               tmbmix%wfnmd%ld_coeff, tmbmix%wfnmd%coeff, Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3d, &
               denspot%rhov, at, denspot%dpcom%nscatterarr)

          ! Mix the density.
          if(trim(input%lin%mixingMethod)=='dens') then
           lscv%compare_outer_loop = pnrm<lscv%self_consistent .or. it_scc==lscv%nit_scc
           call mix_main(iproc, nproc, lscv%mix_hist, lscv%compare_outer_loop, input, glr, lscv%alpha_mix, &
                denspot, mixdiis, rhopotold, rhopotold_out, pnrm, lscv%pnrm_out)
          end if


          ! Calculate the new potential.
          if(iproc==0) write(*,'(1x,a)') '---------------------------------------------------------------- Updating potential.'
          call updatePotential(iproc,nproc,at%geocode,input%ixc,input%nspin,&
               0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,Glr,denspot,ehart,eexcu,vexcu)

          ! Calculate the total energy.
          energy=ebs-ehart+eexcu-vexcu-eexctX+eion+edisp
          energyDiff=energy-energyold
          energyold=energy


          ! Mix the potential
          if(trim(input%lin%mixingMethod)=='pot') then
           lscv%compare_outer_loop = pnrm<lscv%self_consistent .or. it_scc==lscv%nit_scc
           call mix_main(iproc, nproc, lscv%mix_hist, lscv%compare_outer_loop, input, glr, lscv%alpha_mix, &
                denspot, mixdiis, rhopotold, rhopotold_out, pnrm, lscv%pnrm_out)
          end if


          ! Post communications for gathering the potential
          call allocateCommunicationsBuffersPotential(tmb%comgp, subname)
          !!call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, tmb%comgp)
          call post_p2p_communication(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, &
               tmb%comgp%nrecvbuf, tmb%comgp%recvbuf, tmb%comgp)
          if(tmbmix%wfnmd%bs%use_derivative_basis) then
              call allocateCommunicationsBuffersPotential(tmbmix%comgp, subname)
              !!call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, tmbmix%comgp)
              call post_p2p_communication(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, &
                   tmbmix%comgp%nrecvbuf, tmbmix%comgp%recvbuf, tmbmix%comgp)
          end if

          ! Write some informations.
          call printSummary(iproc, it_scc, lscv%info_basis_functions, &
               infoCoeff, pnrm, energy, energyDiff, input%lin%mixingMethod)
          if(pnrm<lscv%self_consistent) then
              lscv%reduce_convergence_tolerance=.true.
              exit
          else
              lscv%reduce_convergence_tolerance=.false.
          end if

          if(it_scc<lscv%nit_scc_when_optimizing .or. it_scc==lscv%nit_scc) then
              call deallocateCommunicationbufferSumrho(tmb%comsr, subname)
              call deallocateCommunicationbufferSumrho(tmbder%comsr, subname)
          !else if(lscv%locreg_increased) then
          !    call deallocateCommunicationbufferSumrho(tmb%comsr, subname)
          end if
      end do

      !!call deallocateCommunicationbufferSumrho(tmb%comsr, subname)
      !!call deallocateCommunicationbufferSumrho(tmbder%comsr, subname)

      ! Print out values related to two iterations of the outer loop.
      if(iproc==0) then
          write(*,'(3x,a,7es18.10)') 'ebs, ehart, eexcu, vexcu, eexctX, eion, edisp', &
              ebs, ehart, eexcu, vexcu, eexctX, eion, edisp
          if(trim(input%lin%mixingMethod)=='dens') then
              write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                   'itout, Delta DENSOUT, energy, energyDiff', itout, lscv%pnrm_out, energy, energy-energyoldout
          else if(trim(input%lin%mixingMethod)=='pot') then
              write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                   'itout, Delta POTOUT, energy energyDiff', itout, lscv%pnrm_out, energy, energy-energyoldout
          end if
      end if
      energyoldout=energy

      ! Deallocate DIIS structures.
      call deallocateDIIS(ldiis)

  end do outerLoop


  !!call cancelCommunicationPotential(iproc, nproc, tmb%comgp)
  !!call gatherPotential(iproc, nproc, tmb%comgp)
  call wait_p2p_communication(iproc, nproc, tmb%comgp)
  call deallocateCommunicationsBuffersPotential(tmb%comgp, subname)
  if(tmbder%wfnmd%bs%use_derivative_basis) then
      !!call cancelCommunicationPotential(iproc, nproc, tmbder%comgp)
      !!call gatherPotential(iproc, nproc, tmbder%comgp)
      call wait_p2p_communication(iproc, nproc, tmbder%comgp)
      call deallocateCommunicationsBuffersPotential(tmbder%comgp, subname)
  end if

  iall=-product(shape(rhopotOld))*kind(rhopotOld)
  deallocate(rhopotOld, stat=istat)
  call memocc(istat, iall, 'rhopotold', subname)
  iall=-product(shape(rhopotold_out))*kind(rhopotold_out)
  deallocate(rhopotold_out, stat=istat)
  call memocc(istat, iall, 'rhopotold_out', subname)

  if(input%lin%mixHist_highaccuracy>0) then
      call deallocateMixrhopotDIIS(mixdiis)
  end if



  ! Allocate the communication buffers for the calculation of the charge density.
  call allocateCommunicationbufferSumrho(iproc, tmbmix%comsr, subname)
  call communicate_basis_for_density(iproc, nproc, tmb%lzd, tmbmix%orbs, tmbmix%psi, tmbmix%comsr)
  call sumrhoForLocalizedBasis2(iproc, nproc, orbs%norb, tmb%lzd, input, hx, hy, hz, tmbmix%orbs, tmbmix%comsr, &
       tmbmix%wfnmd%ld_coeff, tmbmix%wfnmd%coeff, Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3d, denspot%rhov, at,denspot%dpcom%nscatterarr)

  call deallocateCommunicationbufferSumrho(tmbmix%comsr, subname)

  ! Build global orbitals psi (the physical ones).
  if(input%lin%transformToGlobal) then
      call transformToGlobal(iproc, nproc, tmb%lzd, tmbmix%orbs, orbs, comms, input, tmbmix%wfnmd%ld_coeff, &
           tmbmix%wfnmd%coeff, tmbmix%psi, psi, psit)
  end if


  nullify(rho,pot)
  call destroy_DFT_wavefunction(tmb)
  call destroy_DFT_wavefunction(tmbder)
  call deallocate_local_zone_descriptors(tmb%lzd, subname)
  call deallocateBasicArraysInput(input%lin)
  deallocate(tmb%confdatarr)
  deallocate(tmbder%confdatarr)

  iall=-product(shape(lscv%locrad))*kind(lscv%locrad)
  deallocate(lscv%locrad, stat=istat)
  call memocc(istat, iall, 'lscv%locrad', subname)

  call timing(iproc,'WFN_OPT','PR')



end subroutine linearScaling





subroutine printSummary(iproc, itSCC, infoBasisFunctions, infoCoeff, pnrm, energy, energyDiff, mixingMethod)
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
integer,intent(in):: iproc, itSCC, infoBasisFunctions, infoCoeff
real(8),intent(in):: pnrm, energy, energyDiff
character(len=4),intent(in):: mixingMethod

  if(iproc==0) then
      write(*,'(1x,a)') repeat('#',92 + int(log(real(itSCC))/log(10.)))
      write(*,'(1x,a,i0,a)') 'at iteration ', itSCC, ' of the self consistency cycle:'
      if(infoBasisFunctions<0) then
          write(*,'(3x,a)') '- WARNING: basis functions not converged!'
      else
          write(*,'(3x,a,i0,a)') '- basis functions converged in ', infoBasisFunctions, ' iterations.'
      end if
      if(infoCoeff<0) then
          write(*,'(3x,a)') '- WARNING: coefficients not converged!'
      else if(infoCoeff>0) then
          write(*,'(3x,a,i0,a)') '- coefficients converged in ', infoCoeff, ' iterations.'
      else
          write(*,'(3x,a)') '- coefficients obtained by diagonalization.'
      end if
      if(mixingMethod=='dens') then
          write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)') 'it, Delta DENS, energy, energyDiff', itSCC, pnrm, energy, energyDiff
      else if(mixingMethod=='pot') then
          write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)') 'it, Delta POT, energy energyDiff', itSCC, pnrm, energy, energyDiff
      end if
      write(*,'(1x,a)') repeat('#',92 + int(log(real(itSCC))/log(10.)))
  end if

end subroutine printSummary



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

  else

      do iorb=1,lorbs%norbp
          ilr=lorbs%inwhichlocreg(lorbs%isorb+iorb)
          iiat=onwhichatom(lorbs%isorb+iorb)
          confdatarr(iorb)%prefac=input%lin%potentialPrefac_lowaccuracy(at%iatype(iiat))
      end do
      wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE
      wfnmd%bs%nit_basis_optimization=input%lin%nItBasis_lowaccuracy
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
  integer:: istat, iall
  logical:: redefine_derivatives, redefine_standard
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
  if(iproc==0) write(*,'(1x,a,f6.2,a)') 'Reduce the confining potential to ', &
      100.d0*lscv%decrease_factor_total,'% of its initial value.'
  tmb%confdatarr(:)%prefac=lscv%decrease_factor_total*tmb%confdatarr(:)%prefac


  lscv%locreg_increased=.false.
  redefine_derivatives=.false.
  redefine_standard=.false.
  if(lscv%ifail>=3 .and. .not.lscv%lowaccur_converged) then
      lscv%increase_locreg=lscv%increase_locreg+1.d0
      if(iproc==0) then
          write(*,'(1x,a)') 'It seems that the convergence criterion can not be reached with this localization radius.'
          write(*,'(1x,a,f6.2)') 'The localization radius is increased by totally',lscv%increase_locreg
      end if
      lscv%ifail=0
      lscv%locrad=lscv%locrad+lscv%increase_locreg
      if(lscv%withder) then
          redefine_derivatives=.true.
      else
          redefine_standard=.true.
      end if
      lscv%locreg_increased=.true.
  end if

  redefine_derivatives=.false.
  if(lscv%lowaccur_converged) then
      if(iproc==0) then
          write(*,'(1x,a)') 'Increasing the localization radius for the high accuracy part.'
      end if
      lscv%locreg_increased=.true.
      redefine_standard=.true.
  end if
  if(lscv%locreg_increased) then
      call enlarge_locreg(iproc, nproc, hx, hy, hz, .false., tmb%lzd, lscv%locrad, &
           ldiis, denspot, tmb%wfnmd%nphi, tmb%psi, tmb)
  end if
  if(redefine_standard) then
      ! Fake allocation
      allocate(tmb%comsr%sendbuf(1), stat=istat)
      call memocc(istat, tmb%comsr%sendbuf, 'tmb%comsr%sendbuf', subname)
      allocate(tmb%comsr%recvbuf(1), stat=istat)
      call memocc(istat, tmb%comsr%recvbuf, 'tmb%comsr%recvbuf', subname)
      call redefine_locregs_quantities(iproc, nproc, hx, hy, hz, tmb%lzd, tmb, tmb, denspot)
  end if
  if(redefine_derivatives) then
      ! Fake allocation
      allocate(tmbder%comsr%sendbuf(1), stat=istat)
      call memocc(istat, tmbder%comsr%sendbuf, 'tmbder%comsr%sendbuf', subname)
      allocate(tmbder%comsr%recvbuf(1), stat=istat)
      call memocc(istat, tmbder%comsr%recvbuf, 'tmbder%comsr%recvbuf', subname)
      call redefine_locregs_quantities(iproc, nproc, hx, hy, hz, tmb%lzd, tmb, tmbder, denspot)
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
  
  lscv%exit_outer_loop=.false.
  
  if(lscv%lowaccur_converged) then
      lscv%nit_highaccuracy=lscv%nit_highaccuracy+1
      if(lscv%nit_highaccuracy==input%lin%nit_highaccuracy+1) then
          ! Deallocate DIIS structures.
          call deallocateDIIS(ldiis)
          lscv%exit_outer_loop=.true.
      end if
      ! only use steepest descent if the localization regions may change
      if(input%lin%nItInnerLoop/=-1 .or. tmb%wfnmd%bs%locreg_enlargement/=1.d0) then
          ldiis%isx=0
      end if
  
      if(input%lin%mixHist_lowaccuracy==0 .and. input%lin%mixHist_highaccuracy>0) then
          call initializeMixrhopotDIIS(input%lin%mixHist_highaccuracy, denspot%dpcom%ndimpot, mixdiis)
      else if(input%lin%mixHist_lowaccuracy>0 .and. input%lin%mixHist_highaccuracy==0) then
          call deallocateMixrhopotDIIS(mixdiis)
      end if
  end if
  
end subroutine adjust_DIIS_for_high_accuracy



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
