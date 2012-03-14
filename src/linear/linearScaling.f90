subroutine linearScaling(iproc,nproc,Glr,orbs,comms,at,input,hx,hy,hz,&
     rxyz,fion,fdisp,denspot,nlpspd,proj,GPU,&
     eion,edisp,eexctX,scpot,psi,psit,energy,fxyz)
!
! Purpose:
! ========
!   Top-level subroutine for the linear scaling version.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc       process ID
!     nproc       total number of processes
!     n3d         ??
!     n3p         ??
!     i3s         ??
!     i3xcsh      ??
!     Glr         type describing the localization region
!     orbs        type describing the physical orbitals psi
!     comms       type containing the communication parameters for the physical orbitals psi
!     at          type containing the parameters for the atoms
!     input       type  containing some very general parameters
!     lin         type containing parameters for the linear version
!     rxyz        atomic positions
!     nscatterarr ??
!     ngatherarr  ??
!     nlpspd      ??
!     proj        ??
!     pkernelseq  ??
!     radii_cf    coarse and fine radii around the atoms
!     irrzon      ??
!     phnons      ??
!     pkernel     ??
!     pot_ion     the ionic potential
!     rhocore     ??
!     potxc       ??
!     PSquiet     flag to control the output from the Poisson solver
!     eion        ionic energy
!     edisp       dispersion energy
!     eexctX      ??
!     scpot       flag indicating whether we have a self consistent calculation
!     fion        ionic forces
!     fdisp       dispersion forces
!   Input / Output arguments
!   ------------------------
!     GPU         parameters for GPUs?
!     rhopot      the charge density
!   Output arguments:
!   -----------------
!     psi         the physical orbitals
!     psit        psi transposed
!     fxyz        the forces acting on the atom
!     energy
!     

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
!real(8),dimension(orbs),intent(out):: psi
real(8),dimension(:),pointer,intent(out):: psi, psit
real(gp), dimension(:), pointer :: rho,pot
real(8),intent(out):: energy
real(8),dimension(3,at%nat),intent(out):: fxyz
!real(8),intent(out):: fnoise

! Local variables
integer:: infoBasisFunctions,infoCoeff,istat,iall,itSCC,nitSCC,i,ierr,potshortcut,ist,istr,ilr,tag,itout
integer :: jproc,iat,j, nit_highaccuracy, mixHist, nitSCCWhenOptimizing, nit, npsidim,ityp
real(8):: ebs, ebsMod, pnrm, tt, ehart, eexcu, vexcu, alphaMix, trace
character(len=*),parameter:: subname='linearScaling'
real(8),dimension(:),allocatable:: rhopotOld, rhopotold_out, locrad
logical:: reduceConvergenceTolerance, communicate_lphi, with_auxarray, lowaccur_converged, withder, variable_locregs
real(8):: t1, t2, time, t1tot, t2tot, timetot, t1ig, t2ig, timeig, t1init, t2init, timeinit, ddot, dnrm2, pnrm_out
real(8):: t1scc, t2scc, timescc, t1force, t2force, timeforce, energyold, energyDiff, energyoldout, selfConsistent
integer:: iorb, ndimtot, iiat
type(mixrhopotDIISParameters):: mixdiis
type(workarr_sumrho):: w
!real(8),dimension(:,:),allocatable:: coeff_proj
type(localizedDIISParameters):: ldiis
type(confpot_data), dimension(:),pointer :: confdatarr, confdatarrder
real(8):: fnoise,pressure
real(gp), dimension(6) :: ewaldstr,strten,hstrten,xcstr
type(orthon_data):: orthpar
integer,dimension(:),pointer:: onwhichatom
integer,dimension(:),allocatable:: norbsPerAtom
!type(wfn_metadata):: wfnmd
type(DFT_wavefunction),target:: tmb
type(DFT_wavefunction),target:: tmbder
type(DFT_wavefunction),pointer:: tmbmix, tmbopt
type(local_zone_descriptors):: lzd
type(orbitals_data):: orbs_tmp


  if(iproc==0) then
      write(*,'(1x,a)') repeat('*',84)
      write(*,'(1x,a)') '****************************** LINEAR SCALING VERSION ******************************'
  end if
  call mpi_barrier(mpi_comm_world, ierr)
  t1tot=mpi_wtime()

  ! Initialize the parameters for the linear scaling version and allocate all arrays.
  tag=0
  call mpi_barrier(mpi_comm_world, ierr)
  t1init=mpi_wtime()


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

  call init_local_zone_descriptors(iproc, nproc, input, glr, at, rxyz, tmb%orbs, tmbder%orbs, lzd)

  call update_wavefunctions_size(lzd,tmb%orbs)
  call update_wavefunctions_size(lzd,tmbder%orbs)

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

  call initCommsOrtho(iproc, nproc, input%nspin, hx, hy, hz, lzd, tmb%orbs, tmb%orbs%inWhichLocreg,&
       input%lin%locregShape, tmb%op, tmb%comon, tag)
  call initCommsOrtho(iproc, nproc, input%nspin, hx, hy, hz, lzd, tmbder%orbs, tmbder%orbs%inWhichLocreg, &
       input%lin%locregShape, tmbder%op, tmbder%comon, tag)
  
  call initializeCommunicationPotential(iproc, nproc, denspot%dpcom%nscatterarr, &
       tmb%orbs, lzd, tmb%comgp, tmb%orbs%inWhichLocreg, tag)
  call initializeCommunicationPotential(iproc, nproc, denspot%dpcom%nscatterarr, &
       tmbder%orbs, lzd, tmbder%comgp, tmbder%orbs%inWhichLocreg, tag)

  if(input%lin%useDerivativeBasisFunctions) then
      call initializeRepartitionOrbitals(iproc, nproc, tag, tmb%orbs, tmbder%orbs, lzd, tmbder%comrp)
      call initializeRepartitionOrbitals(iproc, nproc, tag, tmb%orbs, tmbder%orbs, lzd, tmb%comrp)
  else
      call nullify_p2pComms(tmbder%comrp)
      call nullify_p2pComms(tmb%comrp)
  end if


  call nullify_p2pcomms(tmb%comsr)
  call initializeCommsSumrho(iproc, nproc, denspot%dpcom%nscatterarr, lzd, tmb%orbs, tag, tmb%comsr)
  call nullify_p2pcomms(tmbder%comsr)
  call initializeCommsSumrho(iproc, nproc, denspot%dpcom%nscatterarr, lzd, tmbder%orbs, tag, tmbder%comsr)

  call initMatrixCompression(iproc, nproc, lzd%nlr, tmb%orbs, tmb%op%noverlaps, tmb%op%overlaps, tmb%mad)
  call initCompressedMatmul3(tmb%orbs%norb, tmb%mad)
  call initMatrixCompression(iproc, nproc, lzd%nlr, tmbder%orbs, &
       tmbder%op%noverlaps, tmbder%op%overlaps, tmbder%mad)
  call initCompressedMatmul3(tmbder%orbs%norb, tmbder%mad)

  allocate(confdatarr(tmb%orbs%norbp))
  call define_confinement_data(confdatarr,tmb%orbs,rxyz,at,&
       input%hx,input%hy,input%hz,input%lin%confpotorder,input%lin%potentialprefac_lowaccuracy,lzd,tmb%orbs%onwhichatom)

  allocate(confdatarrder(tmbder%orbs%norbp))
  call define_confinement_data(confdatarrder,tmbder%orbs,rxyz,at,&
       input%hx,input%hy,input%hz,input%lin%confpotorder,input%lin%potentialprefac_lowaccuracy,lzd,tmbder%orbs%onwhichatom)
  ! Now all initializations are done ######################################################################################


  ! Assign some values to orthpar
  orthpar%methTransformOverlap = tmb%wfnmd%bs%meth_transform_overlap
  orthpar%nItOrtho = input%lin%nItOrtho
  orthpar%blocksize_pdsyev = tmb%wfnmd%bpo%blocksize_pdsyev
  orthpar%blocksize_pdgemm = tmb%wfnmd%bpo%blocksize_pdgemm


  call mpi_barrier(mpi_comm_world, ierr)
  t2init=mpi_wtime()
  timeinit=t2init-t1init

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
  potshortcut=0 ! What is this?
  call mpi_barrier(mpi_comm_world, ierr)
  t1ig=mpi_wtime()
  call inputguessConfinement(iproc, nproc, at, &
       input, hx, hy, hz, lzd, tmb%orbs, rxyz, denspot ,rhopotold, &
       nlpspd, proj, GPU, &
       tmb%psi)
  call mpi_barrier(mpi_comm_world, ierr)
  t2ig=mpi_wtime()
  timeig=t2ig-t1ig
  t1scc=mpi_wtime()


  ! Initialize the DIIS mixing of the potential if required.
  if(input%lin%mixHist_lowaccuracy>0) then
      call initializeMixrhopotDIIS(input%lin%mixHist_lowaccuracy, denspot%dpcom%ndimpot, mixdiis)
  end if

  !end of the initialization part, will later be moved to cluster
  call timing(iproc,'INIT','PR')

  allocate(locrad(lzd%nlr), stat=istat)
  call memocc(istat, locrad, 'locrad', subname)


  if(input%lin%nItInguess>0) then
      ! Post communications for gathering the potential.
      if(input%lin%mixedmode) tmb%wfnmd%bs%use_derivative_basis=.false.
      if(input%lin%mixedmode) tmbder%wfnmd%bs%use_derivative_basis=.false.
      !!call allocateCommunicationsBuffersPotential(tmb%comgp, subname)
      !!call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, tmb%comgp)
      !!if(tmbder%wfnmd%bs%use_derivative_basis) then
      !!    call allocateCommunicationsBuffersPotential(tmbder%comgp, subname)
      !!    call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, tmbder%comgp)
      !!end if

      ! Calculate the Hamiltonian in the basis of the trace minimizing orbitals. Do not improve
      ! the basis functions (therefore update is set to false).
      ! This subroutine will also post the point to point messages needed for the calculation
      ! of the charge density.
      tmb%wfnmd%bs%communicate_phi_for_lsumrho=.true.
      with_auxarray=.false.
      tmb%wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE

      do ilr=1,lzd%nlr
          locrad(ilr)=input%lin%locrad_lowaccuracy(ilr)
      end do

      !!if(input%lin%mixedmode) then
      !!    call allocateCommunicationbufferSumrho(iproc, with_auxarray, tmb%comsr, subname)
      !!    tmbder%wfnmd%bs%use_derivative_basis=.false.
      !!    call getLinearPsi(iproc, nproc, lzd, orbs, tmb%orbs, tmb%orbs, tmb%comsr, &
      !!        tmb%mad, tmb%mad, tmb%op, tmb%op, tmb%comon, tmb%comon, &
      !!        tmb%comgp, tmb%comgp, at, rxyz, &
      !!        denspot, GPU, &
      !!        infoBasisFunctions, infoCoeff, 0, ebs, nlpspd, proj, &
      !!        ldiis, &
      !!        orthpar, confdatarr, tmbder%wfnmd%bpo%blocksize_pdgemm, &
      !!        tmbder%comrp, tmbder%wfnmd%bpo%blocksize_pdsyev, tmbder%wfnmd%bpo%nproc_pdsyev, &
      !!        hx, hy, hz, input%SIC, locrad, tmb, tmbder)
      !!else
      !!    call allocateCommunicationbufferSumrho(iproc,with_auxarray,tmbder%comsr,subname)
      !!    call getLinearPsi(iproc,nproc,lzd,orbs,tmb%orbs,tmbder%orbs,tmbder%comsr,&
      !!        tmb%mad,tmbder%mad,tmb%op,tmbder%op,tmb%comon,&
      !!        tmbder%comon,tmb%comgp,tmbder%comgp,at,rxyz,&
      !!        denspot,GPU,&
      !!        infoBasisFunctions,infoCoeff,0, ebs,nlpspd,proj,&
      !!        ldiis,orthpar,confdatarr,& 
      !!        tmbder%wfnmd%bpo%blocksize_pdgemm,&
      !!        tmbder%comrp,tmbder%wfnmd%bpo%blocksize_pdsyev,tmbder%wfnmd%bpo%nproc_pdsyev,&
      !!        hx,hy,hz,input%SIC, locrad, tmb, tmbder)
      !!end if

      !!! Calculate the charge density.
      !!if(input%lin%mixedmode) then
      !!    call deallocateCommunicationbufferSumrho(tmb%comsr, subname)
      !!else
      !!    call deallocateCommunicationbufferSumrho(tmbder%comsr, subname)
      !!end if

      if(trim(input%lin%mixingMethod)=='dens') then
          rhopotold_out=rhopotold
      end if


      if(trim(input%lin%mixingMethod)=='pot') then
          !!if(input%lin%mixHist_lowaccuracy==0) then
          !!    call mixPotential(iproc, denspot%dpcom%n3p, Glr, input, &
          !!         input%lin%alphaMixWhenFixed_lowaccuracy, rhopotOld, denspot%rhov, pnrm)
          !!else 
          !!    ndimtot=lzd%Glr%d%n1i*lzd%Glr%d%n2i*lzd%Glr%d%n3i
          !!    mixdiis%mis=mod(mixdiis%is,mixdiis%isx)+1
          !!    mixdiis%is=mixdiis%is+1
          !!    call mixrhopotDIIS(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, rhopotold, mixdiis, ndimtot, &
          !!         input%lin%alphaMixWhenFixed_lowaccuracy, 2, pnrm)
          !!end if
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
  call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, tmb%comgp)
  ! If we also use the derivative of the basis functions, also send the potential in this case. This is
  ! needed since the orbitals may be partitioned in a different way when the derivatives are used.
  if(tmbder%wfnmd%bs%use_derivative_basis) then
      call allocateCommunicationsBuffersPotential(tmbder%comgp, subname)
      call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, tmbder%comgp)
  end if




  !if(nproc==1) allocate(psit(size(psi)))
  !nitSCC=input%lin%nitSCCWhenOptimizing+input%lin%nitSCCWhenFixed
  ! Flag that indicates that the basis functions shall be improved in the following.
  tmb%wfnmd%bs%update_phi=.true.
  pnrm=1.d100
  pnrm_out=1.d100
  energyold=0.d0
  energyoldout=0.d0
  reduceConvergenceTolerance=.false.
  tmb%wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE
  lowaccur_converged=.false.

  if(input%lin%useDerivativeBasisFunctions) then
      tmbmix => tmbder
  else
      tmbmix => tmb
  end if

  ! Check whether it is possible to have variable localization regions or not.
  if(tmb%wfnmd%bs%nit_unitary_loop==-1 .and. tmb%wfnmd%bs%locreg_enlargement==1.d0) then
      variable_locregs=.false.
  else
      variable_locregs=.true.
  end if

  outerLoop: do itout=1,input%lin%nit_lowaccuracy+input%lin%nit_highaccuracy

      ! First to some initialization and determine the value of some control parameters.

      ! Initialize DIIS...
      call initializeDIIS(input%lin%DIISHistMax, lzd, tmb%orbs, tmb%orbs%norb, ldiis)
      ldiis%DIISHistMin=input%lin%DIISHistMin
      ldiis%DIISHistMax=input%lin%DIISHistMax
      ldiis%alphaSD=input%lin%alphaSD
      ldiis%alphaDIIS=input%lin%alphaDIIS

      ! The basis functions shall be optimized
      tmb%wfnmd%bs%update_phi=.true.

      ! Convergence criterion for the self consistency looo
      selfConsistent=input%lin%convCritMix

      ! Check whether the derivatives shall be used or not.
      if(input%lin%mixedmode) then
          if( (.not.lowaccur_converged .and. &
               (itout==input%lin%nit_lowaccuracy+1 .or. pnrm_out<input%lin%lowaccuray_converged) ) &
              .or. lowaccur_converged ) then
              withder=.true.
          else
              withder=.false.
          end if
      end if

      ! Check whether the low accuracy part (i.e. with strong confining potential) has converged.
      if(.not.lowaccur_converged .and. (itout==input%lin%nit_lowaccuracy+1 .or. pnrm_out<input%lin%lowaccuray_converged)) then
          lowaccur_converged=.true.
          nit_highaccuracy=0
      end if 

      ! Set all remaining variables that we need for the optimizations of the basis functions and the mixing.
      call set_optimization_variables(lowaccur_converged, input, at, tmb%orbs, lzd%nlr, tmb%orbs%onwhichatom, &
           confdatarr, tmb%wfnmd, locrad, nitSCC, nitSCCWhenOptimizing, mixHist, alphaMix)
      call set_optimization_variables(lowaccur_converged, input, at, tmbder%orbs, lzd%nlr, tmbder%orbs%onwhichatom, &
           confdatarrder, tmbder%wfnmd, locrad, nitSCC, nitSCCWhenOptimizing, mixHist, alphaMix)

      if(tmb%wfnmd%bs%confinement_decrease_mode==DECREASE_ABRUPT) then
          tt=1.d0
      else if(tmb%wfnmd%bs%confinement_decrease_mode==DECREASE_LINEAR) then
          tt=1.d0-(dble(itout-1))/dble(input%lin%nit_lowaccuracy)
          if(iproc==0) write(*,'(1x,a,f6.2,a)') 'Reduce the confining potential to ',100.d0*tt,'% of its initial value.'
      end if
      confdatarr(:)%prefac=tt*confdatarr(:)%prefac
      if(iproc==0) write(*,*) 'confdatarr(1)%prefac',confdatarr(1)%prefac

      ! Somce special treatement if we are in the high accuracy part
      if(lowaccur_converged) then
          nit_highaccuracy=nit_highaccuracy+1
          if(nit_highaccuracy==input%lin%nit_highaccuracy+1) then
              ! Deallocate DIIS structures.
              call deallocateDIIS(ldiis)
              exit outerLoop
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

      ! Allocate the communication arrays for the calculation of the charge density.
      with_auxarray=.false.
      call allocateCommunicationbufferSumrho(iproc, with_auxarray, tmb%comsr, subname)
      call allocateCommunicationbufferSumrho(iproc, with_auxarray, tmbder%comsr, subname)

      ! Now all initializations are done...

      ! The self consistency cycle. Here we try to get a self consistent density/potential.
      ! In the first nitSCCWhenOptimizing iteration, the basis functions are optimized, whereas in the remaining
      ! iteration the basis functions are fixed.
      do itSCC=1,nitSCC
          if(itSCC>nitSCCWhenOptimizing) tmb%wfnmd%bs%update_phi=.false.
          if(itSCC==1) then
              tmb%wfnmd%bs%communicate_phi_for_lsumrho=.true.
              tmbder%wfnmd%bs%communicate_phi_for_lsumrho=.true.
          else
              tmb%wfnmd%bs%communicate_phi_for_lsumrho=.false.
              tmbder%wfnmd%bs%communicate_phi_for_lsumrho=.false.
          end if

          ! Update the basis functions (if wfnmd%bs%update_phi is true), calculate the Hamiltonian in this basis, and diagonalize it.
          ! This is a flag whether the basis functions shall be updated.
          if(tmb%wfnmd%bs%update_phi) then
              ! Improve the trace minimizing orbitals.
              if(itout>1 .and. tmbmix%wfnmd%bs%use_derivative_basis) then
                  do iorb=1,orbs%norb
                      call dcopy(tmb%orbs%norb, tmbmix%wfnmd%coeff_proj(1,iorb), 1, tmb%wfnmd%coeff(1,iorb), 1)
                  end do
              end if
              call getLocalizedBasis(iproc,nproc,at,lzd,tmb%orbs,orbs,tmb%comon,tmb%op,tmb%comgp,tmb%mad,rxyz,&
                  denspot,GPU,trace,&
                  infoBasisFunctions,nlpspd,proj,ldiis,&
                  orthpar,confdatarr,tmb%wfnmd%bpo%blocksize_pdgemm,&
                  hx,hy,hz,input%SIC,locrad,tmb)
              tmb%wfnmd%nphi=tmb%orbs%npsidim_orbs
          end if
          !!if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY .and. tmb%wfnmd%bs%update_phi) then
          !!    call update_locreg(iproc, nproc, tmbder%wfnmd%bs%use_derivative_basis, denspot, hx, hy, hz, &
          !!         tmb%orbs, lzd, tmbder%orbs, tmbder%op, tmbder%comon, tmb%comgp, tmbder%comgp, tmb%comsr, tmbder%mad)
          !!end if
          if(input%lin%mixedmode) then
              if(.not.withder) then
                  tmbder%wfnmd%bs%use_derivative_basis=.false.
                  tmbmix => tmb
              else
                  tmbder%wfnmd%bs%use_derivative_basis=.true.
                  ! We have to communicate the potential in the first iteration
                  if(itSCC==1) then
                      call allocateCommunicationsBuffersPotential(tmbder%comgp, subname)
                      call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, tmbder%comgp)
                  end if
                  tmbmix => tmbder
              end if
          end if
          if(tmbmix%wfnmd%bs%use_derivative_basis) then
              ! Cancel the communication of the potential for the TMB, since we need in the following
              ! only the potential for the TMB including the derivatives.
              call cancelCommunicationPotential(iproc, nproc, tmb%comgp)
              call deallocateCommunicationsBuffersPotential(tmb%comgp, subname)
          end if
          !!!if(iproc==0) then
          !!    do ilr=1,lzd%nlr
          !!        write(*,'(a,2i6,l5)') 'before: iproc, ilr, associated(lzd%llr(ilr)%bounds%kb%ibyz_c)', iproc, ilr, associated(lzd%llr(ilr)%bounds%kb%ibyz_c)
          !!    end do
          !!!end if
          if(variable_locregs .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY &
          !if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY &
              .and. tmb%wfnmd%bs%update_phi) then
              if(tmbmix%wfnmd%bs%use_derivative_basis) then
                  call nullify_orbitals_data(orbs_tmp)
                  call copy_orbitals_data(tmb%orbs, orbs_tmp, subname)
                  call update_locreg(iproc, nproc, tmbmix%wfnmd%bs%use_derivative_basis, denspot, hx, hy, hz, &
                       orbs_tmp, lzd, tmbmix%orbs, tmbmix%op, tmbmix%comon, tmb%comgp, tmbmix%comgp, tmbmix%comsr, tmbmix%mad)
                  call deallocate_orbitals_data(orbs_tmp, subname)

                  tmbmix%wfnmd%nphi=tmbmix%orbs%npsidim_orbs
                  tmb%wfnmd%basis_is=BASIS_IS_ENHANCED


                  ! Reallocate tmbmix%psi, since it might have a new shape
                  iall=-product(shape(tmbmix%psi))*kind(tmbmix%psi)
                  deallocate(tmbmix%psi, stat=istat)
                  call memocc(istat, iall, 'tmbmix%psi', subname)

                  allocate(tmbmix%psi(tmbmix%orbs%npsidim_orbs), stat=istat)
                  call memocc(istat, tmbmix%psi, 'tmbmix%psi', subname)

                  if(.not.tmbmix%wfnmd%bs%use_derivative_basis) call dcopy(tmb%wfnmd%nphi, tmb%psi(1), 1, tmbmix%psi(1), 1)
                  call cancelCommunicationPotential(iproc, nproc, tmbmix%comgp)
                  call deallocateCommunicationsBuffersPotential(tmbmix%comgp, subname)
              else
                  tag=1
                  call deallocateCommunicationbufferSumrho(tmbmix%comsr, subname)
                  call deallocate_p2pComms(tmbmix%comsr, subname)
                  call nullify_p2pComms(tmbmix%comsr)
                  call initializeCommsSumrho(iproc, nproc, denspot%dpcom%nscatterarr, lzd, tmbmix%orbs, tag, tmbmix%comsr)
                  call allocateCommunicationbufferSumrho(iproc, .false., tmbmix%comsr, subname)
              end if
              !!call cancelCommunicationPotential(iproc, nproc, tmb%comgp)
              !!call deallocateCommunicationsBuffersPotential(tmb%comgp, subname)
              call deallocate_p2pComms(tmbmix%comgp, subname)
              call nullify_p2pComms(tmbmix%comgp)
              call initializeCommunicationPotential(iproc, nproc, denspot%dpcom%nscatterarr, tmbmix%orbs, &
                   lzd, tmbmix%comgp, tmbmix%orbs%inWhichLocreg, tag)
              !!if(tmbder%wfnmd%bs%use_derivative_basis) then
              !!    call allocateCommunicationsBuffersPotential(tmbmix%comgp, subname)
              !!    call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, tmbmix%comgp)
              !!end if
              call allocateCommunicationsBuffersPotential(tmbmix%comgp, subname)
              call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, tmbmix%comgp)
          end if
          !!!if(iproc==0) then
          !!    do ilr=1,lzd%nlr
          !!        write(*,'(a,2i6,l5)') 'after: iproc, ilr, associated(lzd%llr(ilr)%bounds%kb%ibyz_c)', iproc, ilr, associated(lzd%llr(ilr)%bounds%kb%ibyz_c)
          !!    end do
          !!!end if


          if(tmb%wfnmd%bs%update_phi .or. itSCC==0) then
              if(tmbmix%wfnmd%bs%use_derivative_basis) then
                  if(variable_locregs .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY &
                      .and. tmb%wfnmd%bs%update_phi) then
                      call deallocate_p2pComms(tmbmix%comrp, subname)
                      call nullify_p2pComms(tmbmix%comrp)
                      call initializeRepartitionOrbitals(iproc, nproc, tag, tmb%orbs, tmbmix%orbs, lzd, tmbmix%comrp)
                  end if
                  if(iproc==0) write(*,'(1x,a)',advance='no') 'calculating derivative basis functions...'
                  call getDerivativeBasisFunctions(iproc,nproc,hx,lzd,tmb%orbs,tmbmix%orbs,tmbmix%comrp,&
                       max(tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp),tmb%psi,tmbmix%psi)
                  if(iproc==0) write(*,'(a)') 'done.'
              else
                  call dcopy(tmb%wfnmd%nphi, tmb%psi(1), 1, tmbmix%psi(1), 1)
              end if
          end if


          ! Calculate the coefficients
          call get_coeff(iproc,nproc,lzd,orbs,at,rxyz,denspot,GPU,infoCoeff,ebs,nlpspd,proj,&
               tmbmix%wfnmd%bpo%blocksize_pdsyev,tmbder%wfnmd%bpo%nproc_pdsyev,&
               hx,hy,hz,input%SIC,tmbmix)


          ! Calculate the charge density.
          call sumrhoForLocalizedBasis2(iproc, nproc, orbs%norb,&
               lzd, input, hx, hy ,hz, tmbmix%orbs, tmbmix%comsr, &
               tmbmix%wfnmd%ld_coeff, tmbmix%wfnmd%coeff, Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3d, &
               denspot%rhov, at, denspot%dpcom%nscatterarr)

          ! Mix the density.
          if(trim(input%lin%mixingMethod)=='dens') then
              if(mixHist==0) then
                  call mixPotential(iproc, denspot%dpcom%n3p, Glr, input, alphaMix, rhopotOld, denspot%rhov, pnrm)
              else 
                  ndimtot=lzd%Glr%d%n1i*lzd%Glr%d%n2i*lzd%Glr%d%n3i
                  mixdiis%mis=mod(mixdiis%is,mixdiis%isx)+1
                  mixdiis%is=mixdiis%is+1
                  call mixrhopotDIIS(iproc, nproc, denspot%dpcom%ndimpot,&
                       denspot%rhov, rhopotold, mixdiis, ndimtot, alphaMix, 1, pnrm)
              end if
              ! Determine the change in the density between this iteration and the last iteration in the outer loop.
              if(pnrm<selfConsistent .or. itSCC==nitSCC) then
                  pnrm_out=0.d0
                  do i=1,Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3p
                      pnrm_out=pnrm_out+(denspot%rhov(i)-rhopotOld_out(i))**2
                  end do
                  call mpiallred(pnrm_out, 1, mpi_sum, mpi_comm_world, ierr)
                  pnrm_out=sqrt(pnrm_out)/(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*input%nspin)
                  call dcopy(max(Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotOld_out(1), 1)
              end if
          end if

          ! Copy the current charge density.
          if(trim(input%lin%mixingMethod)=='dens') then
              call dcopy(max(Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotOld(1), 1)
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
              if(mixHist==0) then
                  call mixPotential(iproc, denspot%dpcom%n3p, Glr, input, alphaMix, rhopotOld, denspot%rhov, pnrm)
              else 
                  ndimtot=lzd%Glr%d%n1i*lzd%Glr%d%n2i*lzd%Glr%d%n3i
                  mixdiis%mis=mod(mixdiis%is,mixdiis%isx)+1
                  mixdiis%is=mixdiis%is+1
                  call mixrhopotDIIS(iproc, nproc, denspot%dpcom%ndimpot,&
                       denspot%rhov, rhopotold, mixdiis, ndimtot, alphaMix, 2, pnrm)
              end if

              ! Determine the change in the density between this iteration and the last iteration 
              ! of the previous iteration in the outer loop.
              if(pnrm<selfConsistent .or. itSCC==nitSCC) then
                  pnrm_out=0.d0
                  do i=1,Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3p
                      pnrm_out=pnrm_out+(denspot%rhov(i)-rhopotOld_out(i))**2
                  end do
                  call mpiallred(pnrm_out, 1, mpi_sum, mpi_comm_world, ierr)
                  pnrm_out=sqrt(pnrm_out)/(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*input%nspin)
                  call dcopy(max(Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotOld_out(1), 1)
              end if
          end if

          ! Copy the current potential
          if(trim(input%lin%mixingMethod)=='pot') then
               call dcopy(max(Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotOld(1), 1)
          end if

          ! Post communications for gathering the potential
          call allocateCommunicationsBuffersPotential(tmb%comgp, subname)
          call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, tmb%comgp)
          if(tmbmix%wfnmd%bs%use_derivative_basis) then
              call allocateCommunicationsBuffersPotential(tmbmix%comgp, subname)
              call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, tmbmix%comgp)
          end if

          ! Write some informations.
          call printSummary(iproc, itSCC, infoBasisFunctions, &
               infoCoeff, pnrm, energy, energyDiff, input%lin%mixingMethod)
          if(pnrm<selfConsistent) then
              reduceConvergenceTolerance=.true.
              exit
          else
              reduceConvergenceTolerance=.false.
          end if
      end do

      call deallocateCommunicationbufferSumrho(tmb%comsr, subname)
      call deallocateCommunicationbufferSumrho(tmbder%comsr, subname)

      ! Print out values related to two iterations of the outer loop.
      if(iproc==0) then
          write(*,'(3x,a,7es18.10)') 'ebs, ehart, eexcu, vexcu, eexctX, eion, edisp', &
              ebs, ehart, eexcu, vexcu, eexctX, eion, edisp
          if(trim(input%lin%mixingMethod)=='dens') then
              write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                   'itout, Delta DENSOUT, energy, energyDiff', itout, pnrm_out, energy, energy-energyoldout
          else if(trim(input%lin%mixingMethod)=='pot') then
              write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                   'itout, Delta POTOUT, energy energyDiff', itout, pnrm_out, energy, energy-energyoldout
          end if
      end if
      !!if(abs(pnrm_out)<lin%convCritMixOut) exit
      energyoldout=energy

      ! Deallocate DIIS structures.
      call deallocateDIIS(ldiis)

  end do outerLoop


  call cancelCommunicationPotential(iproc, nproc, tmb%comgp)
  call deallocateCommunicationsBuffersPotential(tmb%comgp, subname)
  !!if(wfnmd%bs%use_derivative_basis) then
  if(tmbder%wfnmd%bs%use_derivative_basis) then
      call cancelCommunicationPotential(iproc, nproc, tmbder%comgp)
      call deallocateCommunicationsBuffersPotential(tmbder%comgp, subname)
  end if

  iall=-product(shape(rhopotOld))*kind(rhopotOld)
  deallocate(rhopotOld, stat=istat)
  call memocc(istat, iall, 'rhopotold', subname)
  iall=-product(shape(rhopotold_out))*kind(rhopotold_out)
  deallocate(rhopotold_out, stat=istat)
  call memocc(istat, iall, 'rhopotold_out', subname)
  !!iall=-product(shape(onwhichatom))*kind(onwhichatom)
  !!deallocate(onwhichatom, stat=istat)
  !!call memocc(istat, iall, 'onwhichatom', subname)
  !!iall=-product(shape(norbsPerAtom))*kind(norbsPerAtom)
  !!deallocate(norbsPerAtom, stat=istat)
  !!call memocc(istat, iall, 'norbsPerAtom', subname)

  if(input%lin%mixHist_highaccuracy>0) then
      call deallocateMixrhopotDIIS(mixdiis)
  end if

  call mpi_barrier(mpi_comm_world, ierr)
  t2scc=mpi_wtime()
  timescc=t2scc-t1scc


  ! Allocate the communication buffers for the calculation of the charge density.
  with_auxarray=.false.
  call allocateCommunicationbufferSumrho(iproc, with_auxarray, tmbmix%comsr, subname)
  call communicate_basis_for_density(iproc, nproc, lzd, tmbmix%orbs, tmbmix%psi, tmbmix%comsr)
  call sumrhoForLocalizedBasis2(iproc, nproc, orbs%norb, lzd, input, hx, hy, hz, tmbmix%orbs, tmbmix%comsr, &
       tmbmix%wfnmd%ld_coeff, tmbmix%wfnmd%coeff, Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3d, denspot%rhov, at,denspot%dpcom%nscatterarr)

  call deallocateCommunicationbufferSumrho(tmbmix%comsr, subname)

  call mpi_barrier(mpi_comm_world, ierr)
  t1force=mpi_wtime()
  ! Build global orbitals psi (the physical ones).
  if(input%lin%transformToGlobal) then
      call transformToGlobal(iproc, nproc, lzd, tmbmix%orbs, orbs, comms, input, tmbmix%wfnmd%ld_coeff, &
           tmbmix%wfnmd%coeff, tmbmix%psi, psi, psit)
  end if



  ! Put the timings here since there is a crash in the forces.
  call mpi_barrier(mpi_comm_world, ierr)
  t2tot=mpi_wtime()
  timetot=t2tot-t1tot
  timeforce=0.d0
  if(iproc==0) write(*,'(1x,a)') '================================================'
  if(iproc==0) write(*,'(1x,a,es10.3,a)') 'total time for linear scaling version:',timetot,'s'
  if(iproc==0) write(*,'(3x,a)') 'of which:'
  if(iproc==0) write(*,'(13x,a,es10.3,a,f4.1,a)') '- initialization:',timeinit,'s (',timeinit/timetot*100.d0,'%)'
  if(iproc==0) write(*,'(13x,a,es10.3,a,f4.1,a)') '- input guess:',timeig,'s (',timeig/timetot*100.d0,'%)'
  if(iproc==0) write(*,'(13x,a,es10.3,a,f4.1,a)') '- self consistency cycle:',timescc,'s (',timescc/timetot*100.d0,'%)'
  if(iproc==0) write(*,'(13x,a,es10.3,a,f4.1,a)') '- forces:',timeforce,'s (',timeforce/timetot*100.d0,'%)'
  if(iproc==0) write(*,'(1x,a)') '================================================'


  ! Calculate the forces we get with psi.
  !!call calculateForcesSub(iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh, Glr, orbs, at, input, hx, hy, hz, &
  !! comms, lin, nlpspd, proj, ngatherarr, nscatterarr, GPU, irrzon, phnons, pkernel, rxyz, fion, fdisp,&
  !! lphi, coeff, rhopot, fxyz, fnoise,radii_cf)

  !!!!associate the density
  !!!rho => rhopot

  !!!!add an if statement which says whether the charge density has already been calculated
  !!!call density_and_hpot(iproc,nproc,at%geocode,at%sym,orbs,lzd,&
  !!!     0.5_gp*input%hx,0.5_gp*input%hy,0.5_gp*input%hz,nscatterarr,&
  !!!     pkernel,rhodsc,GPU,psi,rho,pot,hstrten)

  !!!!fake ewald stress tensor
  !!!ewaldstr=0.0_gp
  !!!xcstr=0.0_gp
  !!!call calculate_forces(iproc,nproc,Glr,at,orbs,nlpspd,rxyz,&
  !!!     input%hx,input%hy,input%hz,proj,i3s+i3xcsh,n3p,&
  !!!     input%nspin,.false.,ngatherarr,rho,pot,potxc,psi,fion,fdisp,fxyz,&
  !!!     ewaldstr,hstrten,xcstr,strten,fnoise,pressure,0.0_dp)

!  iall=-product(shape(pot))*kind(pot)
!  deallocate(pot,stat=istat)
!  call memocc(istat,iall,'pot',subname)
  !no need of deallocating rho
  nullify(rho,pot)

  !!if(iproc==0) then
  !!   write(*,'(1x,a)') 'Force values for all atoms in x, y, z direction.'
  !!   do iat=1,at%nat
  !!      write(*,'(3x,i0,1x,a6,1x,3(1x,es17.10))') &
  !!           iat,trim(at%atomnames(at%iatype(iat))),(fxyz(j,iat),j=1,3)
  !!   end do
  !!end if


!!$  call calculateForcesLinear(iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh, Glr, orbs, at, input, hx, hy, hz,&
!!$   comms, lin, nlpspd, proj, ngatherarr, nscatterarr, GPU, irrzon, phnons, pkernel, rxyz, fion, fdisp,&
!!$   rhopot, psi, fxyz, fnoise)
  !!call mpi_barrier(mpi_comm_world, ierr)
  t2force=mpi_wtime()
  timeforce=t2force-t1force



  !!! Deallocate all arrays related to the linear scaling version.
  !!!!call deallocateLinear(iproc, lin, lphi, coeff)
  !!call deallocate_linearParameters(lin, subname)

  !!call destroy_wfn_metadata(wfnmd)
  call destroy_DFT_wavefunction(tmb)
  call destroy_DFT_wavefunction(tmbder)
  call deallocate_local_zone_descriptors(lzd, subname)

  !call deallocateBasicArraysInput(at, input%lin)
  call deallocateBasicArraysInput(input%lin)

  deallocate(confdatarr)
  deallocate(confdatarrder)
  !!call deallocateBasicArrays(lin)

  iall=-product(shape(locrad))*kind(locrad)
  deallocate(locrad, stat=istat)
  call memocc(istat, iall, 'locrad', subname)

  !!iall=-product(shape(coeff_proj))*kind(coeff_proj)
  !!deallocate(coeff_proj, stat=istat)
  !!call memocc(istat, iall, 'coeff_proj', subname)

  ! End of linear scaling part, except of the forces.
  call timing(iproc,'WFN_OPT','PR')


  !!!!call mpi_barrier(mpi_comm_world, ierr)
  !!t2tot=mpi_wtime()
  !!timetot=t2tot-t1tot
  !!if(iproc==0) write(*,'(1x,a)') '================================================'
  !!if(iproc==0) write(*,'(1x,a,es10.3,a)') 'total time for linear scaling version:',timetot,'s'
  !!if(iproc==0) write(*,'(3x,a)') 'of which:'
  !!if(iproc==0) write(*,'(13x,a,es10.3,a,f4.1,a)') '- initialization:',timeinit,'s (',timeinit/timetot*100.d0,'%)'
  !!if(iproc==0) write(*,'(13x,a,es10.3,a,f4.1,a)') '- input guess:',timeig,'s (',timeig/timetot*100.d0,'%)'
  !!if(iproc==0) write(*,'(13x,a,es10.3,a,f4.1,a)') '- self consistency cycle:',timescc,'s (',timescc/timetot*100.d0,'%)'
  !!if(iproc==0) write(*,'(13x,a,es10.3,a,f4.1,a)') '- forces:',timeforce,'s (',timeforce/timetot*100.d0,'%)'
  !!if(iproc==0) write(*,'(1x,a)') '================================================'

end subroutine linearScaling




subroutine mixPotential(iproc, n3p, Glr, input, alphaMix, rhopotOld, rhopot, pnrm)
!
! Purpose:
! ========
!   Mixes the potential in order to get a self consistent potential.
!
! Calling arguments:
! ==================
!   Input arguments
!   ---------------
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, n3p
type(locreg_descriptors),intent(in) :: Glr
type(input_variables),intent(in):: input
real(8),intent(in):: alphaMix
real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin),intent(in):: rhopotOld
real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin),intent(in out):: rhopot
real(8),intent(out):: pnrm

! Local variables
integer:: i, ierr
real(8):: tt

  call timing(iproc,'mix_linear    ','ON')

  pnrm=0.d0
  tt=1.d0-alphaMix
  !do i=1,max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin
  !do i=1,max(Glr%d%n1i*Glr%d%n2i*n3p,1)
  do i=1,Glr%d%n1i*Glr%d%n2i*n3p
      pnrm=pnrm+(rhopot(i)-rhopotOld(i))**2
      rhopot(i)=tt*rhopotOld(i)+alphaMix*rhopot(i)
  end do
  call mpiallred(pnrm, 1, mpi_sum, mpi_comm_world, ierr)
  pnrm=sqrt(pnrm)/(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*input%nspin)

  call timing(iproc,'mix_linear    ','OF')

end subroutine mixPotential




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


subroutine cancelCommunicationPotential(iproc, nproc, comgp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
!type(p2pCommsGatherPot),intent(inout):: comgp
type(p2pComms),intent(inout):: comgp

! Local variables
integer:: jproc, kproc, ierr
integer,dimension(mpi_status_size):: stat
logical:: sendComplete, receiveComplete

! Cancel all communications. 
! It gives errors, therefore simply wait for the communications to complete.
do jproc=0,nproc-1
    do kproc=1,comgp%noverlaps(jproc)
        !call mpi_test(comgp%comarr(7,kproc,jproc), sendComplete, stat, ierr)
        !call mpi_test(comgp%comarr(8,kproc,jproc), receiveComplete, stat, ierr)
        !if(sendComplete .and. receiveComplete) cycle
        !call mpi_cancel(comgp%comarr(7,kproc,jproc), ierr)
        !call mpi_cancel(comgp%comarr(8,kproc,jproc), ierr)
        call mpi_wait(comgp%comarr(7,kproc,jproc), stat, ierr)
        call mpi_wait(comgp%comarr(8,kproc,jproc), stat, ierr)
    end do
end do

end subroutine cancelCommunicationPotential


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

  !do iall=0,nproc-1
  !    write(*,'(a,i5,4i12)') 'START transformToGlobal: iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)', iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)  
  !end do

  allocate(phi(max(gorbs%npsidim_orbs,gorbs%npsidim_comp)+ndebug), stat=istat)
  call memocc(istat, phi, 'phi', subname)
  allocate(phiWork(max(size(phi),size(psi))), stat=istat)
  call memocc(istat, phiWork, 'phiWork', subname)

  ind1=1
  ind2=1
!  phi=0.d0
  if (max(gorbs%npsidim_orbs,gorbs%npsidim_comp) > 0) &
       call to_zero(max(gorbs%npsidim_orbs,gorbs%npsidim_comp),phi(1))

  do iorb=1,lorbs%norbp
      !ilr = lorbs%inWhichLocregp(iorb)
      ilr = lorbs%inWhichLocreg(lorbs%isorb+iorb)
      ldim=lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
      gdim=lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f
      call Lpsi_to_global2(iproc,nproc,ldim,gdim,lorbs%norb,lorbs%nspinor,input%nspin,lzd%Glr,&
           lzd%Llr(ilr),lphi(ind2),phi(ind1))
      ind1=ind1+lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f
      ind2=ind2+lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
  end do
  !if(ind1/=lin%gorbs%npsidim+1) then
  !    write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': ind1/=lin%gorbs%npsidim',ind1,lin%gorbs%npsidim
  !end if
  !if(ind2/=lorbs%npsidim+1) then
  !    write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': ind1/=lorbs%npsidim',ind1,lorbs%npsidim
  !end if
  !do iall=0,nproc-1
  !    write(*,'(a,i5,4i12)') 'after loop: iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)', iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)  
  !end do
  call transpose_v(iproc, nproc, lorbs, lzd%Glr%wfd, gcomms, phi, work=phiWork)
  !do iall=0,nproc-1
  !    write(*,'(a,i5,4i12)') 'after transpose phi: iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)', iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)  
  !end do


  if(iproc==0) then
      write(*,'(1x,a)', advance='no') '------------------------------------- Building linear combinations... '
  end if
  ! Build the extended orbital psi as a linear combination of localized basis functions phi. for real O(N)
  ! this has to replaced, but at the moment it is still needed.
  !call buildWavefunctionModified(iproc, nproc, orbs, gorbs, comms, gcomms, phi, psi, coeff)
  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  !write(*,*) 'iproc, nvctrp', iproc, nvctrp
  !write(*,*) 'iproc, orbs%npsidim', iproc, orbs%npsidim
  call dgemm('n', 'n', nvctrp, orbs%norb, lorbs%norb, 1.d0, phi(1), nvctrp, coeff(1,1), &
       lorbs%norb, 0.d0, psi(1), nvctrp)

  !do iall=0,nproc-1
  !    write(*,'(a,i5,4i12)') 'after buildWavefunctionModified: iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)', iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)  
  !end do


  if(nproc>1) then
      call dcopy(orbs%npsidim_comp, psi, 1, psit, 1)
  else
      psit => psi
  end if

  !call untranspose_v(iproc, nproc, lorbs, lzd%Glr%wfd, gcomms, phi, work=phiWork)
!  do iall=0,nproc-1
!      write(*,'(a,i5,4i12)') 'after untranspose phi: iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)', iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)  
!  end do
!
  !do iall=0,nproc-1
  !    write(*,'(a,i5,4i12)') 'iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)', iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)  
  !end do
  !call mpi_barrier(mpi_comm_world, iall)
  !flush(6)
  !stop
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


subroutine create_wfn_metadata(mode, nphi, lnorb, llbnorb, norb, input, wfnmd)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  character(len=1),intent(in):: mode
  integer,intent(in):: nphi, lnorb, llbnorb, norb
  type(input_variables),intent(in):: input
  type(wfn_metadata),intent(out):: wfnmd

  ! Local variables
  integer:: istat
  character(len=*),parameter:: subname='create_wfn_metadata'

  ! Determine which variables we need, depending on the mode we are in.
  if(mode=='l') then
      ! linear scaling mode

      wfnmd%nphi=nphi
      !!wfnmd%nlbphi=nlbphi
      wfnmd%basis_is=BASIS_IS_ENHANCED !since always it is allocated with wfnmd%nlbphi
      wfnmd%ld_coeff=llbnorb !leading dimension of the coeff array

      !!allocate(wfnmd%phi(wfnmd%nlbphi), stat=istat)
      !!call memocc(istat, wfnmd%phi, 'wfnmd%phi', subname)

      !!allocate(wfnmd%phiRestart(wfnmd%nphi), stat=istat)
      !!call memocc(istat, wfnmd%phiRestart, 'wfnmd%phiRestart', subname)

      allocate(wfnmd%coeff(llbnorb,norb), stat=istat)
      call memocc(istat, wfnmd%coeff, 'wfnmd%coeff', subname)

      allocate(wfnmd%coeff_proj(lnorb,norb), stat=istat)
      call memocc(istat, wfnmd%coeff_proj, 'wfnmd%coeff_proj', subname)

      call init_basis_specifications(input, wfnmd%bs)
      call init_basis_performance_options(input, wfnmd%bpo)

  else if(mode=='c') then
      ! cubic scaling mode

      !!nullify(wfnmd%phi)
      !!nullify(wfnmd%phiRestart)
      nullify(wfnmd%coeff)
      !!nullify(wfnmd%coeff_proj)
  else
      stop 'wrong mode'
  end if

end subroutine create_wfn_metadata


subroutine destroy_wfn_metadata(wfnmd)
  use module_base
  use module_types
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(wfn_metadata),intent(inout):: wfnmd

  ! Local variables
  integer:: istat, iall
  character(len=*),parameter:: subname='destroy_wfn_metadata'

  !!call checkAndDeallocatePointer(wfnmd%phi)
  !!call checkAndDeallocatePointer(wfnmd%phiRestart)
  !!call checkAndDeallocatePointer(wfnmd%coeff)
  !!call checkAndDeallocatePointer(wfnmd%coeff_proj)

  !!iall=-product(shape(wfnmd%phi))*kind(wfnmd%phi)
  !!deallocate(wfnmd%phi, stat=istat)
  !!call memocc(istat, iall, 'wfnmd%phi', subname)

  !!iall=-product(shape(wfnmd%phiRestart))*kind(wfnmd%phiRestart)
  !!deallocate(wfnmd%phiRestart, stat=istat)
  !!call memocc(istat, iall, 'wfnmd%phiRestart', subname)

  iall=-product(shape(wfnmd%coeff))*kind(wfnmd%coeff)
  deallocate(wfnmd%coeff, stat=istat)
  call memocc(istat, iall, 'wfnmd%coeff', subname)

  iall=-product(shape(wfnmd%coeff_proj))*kind(wfnmd%coeff_proj)
  deallocate(wfnmd%coeff_proj, stat=istat)
  call memocc(istat, iall, 'wfnmd%coeff_proj', subname)

end subroutine destroy_wfn_metadata


subroutine init_basis_specifications(input, bs)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(input_variables),intent(in):: input
  type(basis_specifications),intent(out):: bs
  
  bs%update_phi=.false.
  bs%communicate_phi_for_lsumrho=.false.
  bs%use_derivative_basis=input%lin%useDerivativeBasisFunctions
  bs%conv_crit=input%lin%convCrit
  bs%target_function=TARGET_FUNCTION_IS_TRACE
  bs%meth_transform_overlap=input%lin%methTransformOverlap
  bs%nit_precond=input%lin%nitPrecond
  bs%locreg_enlargement=input%lin%factor_enlarge
  bs%nit_basis_optimization=input%lin%nItBasis_lowaccuracy
  bs%nit_unitary_loop=input%lin%nItInnerLoop
  bs%confinement_decrease_mode=input%lin%confinement_decrease_mode

end subroutine init_basis_specifications


subroutine init_basis_performance_options(input, bpo)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(input_variables),intent(in):: input
  type(basis_performance_options),intent(out):: bpo
  
  bpo%blocksize_pdgemm=input%lin%blocksize_pdgemm
  bpo%blocksize_pdsyev=input%lin%blocksize_pdsyev
  bpo%nproc_pdsyev=input%lin%nproc_pdsyev

end subroutine init_basis_performance_options



subroutine set_optimization_variables(lowaccur_converged, input, at, lorbs, nlr, onwhichatom, confdatarr, wfnmd, &
           locrad, nitSCC, nitSCCWhenOptimizing, mixHist, alphaMix)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  logical,intent(in):: lowaccur_converged
  integer,intent(in):: nlr
  type(orbitals_data),intent(in):: lorbs
  type(input_variables),intent(in):: input
  type(atoms_data),intent(in):: at
  integer,dimension(lorbs%norb),intent(in):: onwhichatom
  type(confpot_data),dimension(lorbs%norbp),intent(inout):: confdatarr
  type(wfn_metadata),intent(inout):: wfnmd
  real(8),dimension(nlr),intent(out):: locrad
  integer,intent(out):: nitSCC, nitSCCWhenOptimizing, mixHist
  real(8),intent(out):: alphaMix

  ! Local variables
  integer:: iorb, ilr, iiat

  if(lowaccur_converged) then

      do iorb=1,lorbs%norbp
          ilr=lorbs%inwhichlocreg(lorbs%isorb+iorb)
          iiat=onwhichatom(lorbs%isorb+iorb)
          confdatarr(iorb)%prefac=input%lin%potentialPrefac_highaccuracy(at%iatype(iiat))
      end do
      wfnmd%bs%target_function=TARGET_FUNCTION_IS_ENERGY
      wfnmd%bs%nit_basis_optimization=input%lin%nItBasis_highaccuracy
      nitSCC=input%lin%nitSCCWhenOptimizing_highaccuracy+input%lin%nitSCCWhenFixed_highaccuracy
      nitSCCWhenOptimizing=input%lin%nitSCCWhenOptimizing_highaccuracy
      mixHist=input%lin%mixHist_highaccuracy
      do ilr=1,nlr
          locrad(ilr)=input%lin%locrad_highaccuracy(ilr)
      end do
      if(wfnmd%bs%update_phi) then
          alphaMix=input%lin%alphaMixWhenOptimizing_highaccuracy
      else
          alphaMix=input%lin%alphaMixWhenFixed_highaccuracy
      end if

  else

      do iorb=1,lorbs%norbp
          ilr=lorbs%inwhichlocreg(lorbs%isorb+iorb)
          iiat=onwhichatom(lorbs%isorb+iorb)
          confdatarr(iorb)%prefac=input%lin%potentialPrefac_lowaccuracy(at%iatype(iiat))
      end do
      wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE
      wfnmd%bs%nit_basis_optimization=input%lin%nItBasis_lowaccuracy
      nitSCC=input%lin%nitSCCWhenOptimizing_lowaccuracy+input%lin%nitSCCWhenFixed_lowaccuracy
      nitSCCWhenOptimizing=input%lin%nitSCCWhenOptimizing_lowaccuracy
      mixHist=input%lin%mixHist_lowaccuracy
      do ilr=1,nlr
          locrad(ilr)=input%lin%locrad_lowaccuracy(ilr)
      end do
      if(wfnmd%bs%update_phi) then
          alphaMix=input%lin%alphaMixWhenOptimizing_lowaccuracy
      else
          alphaMix=input%lin%alphaMixWhenFixed_lowaccuracy
      end if

  end if

end subroutine set_optimization_variables



subroutine create_DFT_wavefunction(mode, nphi, lnorb, norb, input, wfn)
  use module_base
  use module_types
  use module_interfaces, except_this_one => create_DFT_wavefunction
  implicit none
  
  ! Calling arguments
  character(len=1),intent(in):: mode
  integer,intent(in):: nphi, lnorb, norb
  type(input_variables),intent(in):: input
  type(DFT_wavefunction),intent(out):: wfn

  ! Local variables
  integer:: istat
  character(len=*),parameter:: subname='create_DFT_wavefunction'

  call create_wfn_metadata(mode, nphi, lnorb, lnorb, norb, input, wfn%wfnmd)

  allocate(wfn%psi(wfn%wfnmd%nphi), stat=istat)
  call memocc(istat, wfn%psi, 'wfn%psi', subname)

end subroutine create_DFT_wavefunction



subroutine destroy_DFT_wavefunction(wfn)
  use module_base
  use module_types
  use module_interfaces, except_this_one => destroy_DFT_wavefunction
  use deallocatePointers
  implicit none
  
  ! Calling arguments
  type(DFT_wavefunction),intent(inout):: wfn

  ! Local variables
  integer:: istat, iall
  character(len=*),parameter:: subname='destroy_DFT_wavefunction'

  iall=-product(shape(wfn%psi))*kind(wfn%psi)
  deallocate(wfn%psi, stat=istat)
  call memocc(istat, iall, 'wfn%psi', subname)

  call deallocate_overlapParameters(wfn%op, subname)
  call deallocate_p2pComms(wfn%comon, subname)
  call deallocate_p2pComms(wfn%comgp, subname)
  call deallocate_p2pComms(wfn%comrp, subname)
  call deallocate_p2pComms(wfn%comsr, subname)
  call deallocate_matrixDescriptors(wfn%mad, subname)
  call deallocate_orbitals_data(wfn%orbs, subname)
  call deallocate_communications_arrays(wfn%comms, subname)
  call destroy_wfn_metadata(wfn%wfnmd)

end subroutine destroy_DFT_wavefunction


subroutine update_wavefunctions_size(lzd,orbs)
use module_base
use module_types
implicit none

! Calling arguments
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(inout):: orbs

! Local variables
integer:: npsidim, ilr, iorb

  npsidim = 0
  do iorb=1,orbs%norbp
   ilr=orbs%inwhichlocreg(iorb+orbs%isorb)
   npsidim = npsidim + lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
  end do
  orbs%npsidim_orbs=max(npsidim,1)

end subroutine update_wavefunctions_size
