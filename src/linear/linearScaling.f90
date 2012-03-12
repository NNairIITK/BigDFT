subroutine linearScaling(iproc,nproc,Glr,orbs,comms,at,input,hx,hy,hz,&
     lin,rxyz,fion,fdisp,denspot,nlpspd,proj,GPU,&
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
type(linearParameters),intent(inout):: lin
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
integer :: jproc,iat,j, nit_highaccuracy, mixHist, nitSCCWhenOptimizing, nit
real(8):: ebs, ebsMod, pnrm, tt, ehart, eexcu, vexcu, alphaMix
character(len=*),parameter:: subname='linearScaling'
real(8),dimension(:),allocatable:: rhopotOld, rhopotold_out, locrad
logical:: reduceConvergenceTolerance, communicate_lphi, with_auxarray, lowaccur_converged, withder
real(8):: t1, t2, time, t1tot, t2tot, timetot, t1ig, t2ig, timeig, t1init, t2init, timeinit, ddot, dnrm2, pnrm_out
real(8):: t1scc, t2scc, timescc, t1force, t2force, timeforce, energyold, energyDiff, energyoldout, selfConsistent
integer:: iorb, ndimtot, iiat
type(mixrhopotDIISParameters):: mixdiis
type(workarr_sumrho):: w
!real(8),dimension(:,:),allocatable:: coeff_proj
type(localizedDIISParameters):: ldiis
type(confpot_data), dimension(:),pointer :: confdatarr
real(8):: fnoise,pressure
real(gp), dimension(6) :: ewaldstr,strten,hstrten,xcstr
type(orthon_data):: orthpar
integer,dimension(:),pointer:: onwhichatom
type(wfn_metadata):: wfnmd
type(DFT_wavefunction):: tmb
type(DFT_wavefunction):: tmbder


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
  call allocateAndInitializeLinear(iproc, nproc, Glr, orbs, at, nlpspd, lin, &
       input, hx, hy, hz, rxyz, denspot%dpcom%nscatterarr, tag, confdatarr, onwhichatom)


  call create_wfn_metadata('l', max(lin%orbs%npsidim_orbs,lin%orbs%npsidim_comp), &
       max(lin%lb%orbs%npsidim_orbs,lin%lb%orbs%npsidim_comp), &
       lin%orbs%norb, lin%lb%orbs%norb, orbs%norb, input, wfnmd)

  call create_DFT_wavefunction('l', max(lin%orbs%npsidim_orbs,lin%orbs%npsidim_comp), &
       lin%orbs%norb, orbs%norb, input, tmb)

  call create_DFT_wavefunction('l', max(lin%lb%orbs%npsidim_orbs,lin%lb%orbs%npsidim_comp), &
       lin%lb%orbs%norb, orbs%norb, input, tmbder)


  !!lin%potentialPrefac=lin%potentialPrefac_lowaccuracy
  !!allocate(confdatarr(lin%orbs%norbp))
  !!!use a temporary array onwhichatom instead of inwhichlocreg
  !!
  !!call define_confinement_data(confdatarr,lin%orbs,rxyz,at,&
  !!     hx,hy,hz,lin,lin%lzd,lin%orbs%inWhichLocreg)


  !!orthpar%methTransformOverlap = wfnmd%bs%meth_transform_overlap
  orthpar%methTransformOverlap = tmb%wfnmd%bs%meth_transform_overlap
  orthpar%nItOrtho = lin%nItOrtho
  !!orthpar%blocksize_pdsyev = wfnmd%bpo%blocksize_pdsyev
  orthpar%blocksize_pdsyev = tmb%wfnmd%bpo%blocksize_pdsyev
  !!orthpar%blocksize_pdgemm = wfnmd%bpo%blocksize_pdgemm
  orthpar%blocksize_pdgemm = tmb%wfnmd%bpo%blocksize_pdgemm


  call mpi_barrier(mpi_comm_world, ierr)
  t2init=mpi_wtime()
  timeinit=t2init-t1init

  if(.not.lin%transformToGlobal) then
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
  allocate(rhopotold(max(glr%d%n1i*glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin), stat=istat)
  call memocc(istat, rhopotold, 'rhopotold', subname)
  allocate(rhopotold_out(max(glr%d%n1i*glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin), stat=istat)
  call memocc(istat, rhopotold_out, 'rhopotold_out', subname)
  !rhopotold_out=1.d100

  !!allocate(coeff_proj(lin%orbs%norb,orbs%norb), stat=istat)
  !!call memocc(istat, coeff_proj, 'coeff_proj', subname)

  !write(*,'(a,100i6)') 'lin%orbs%inwhichlocreg', lin%orbs%inwhichlocreg


  potshortcut=0 ! What is this?
  call mpi_barrier(mpi_comm_world, ierr)
  t1ig=mpi_wtime()
  call inputguessConfinement(iproc, nproc, at, &
       input, hx, hy, hz, lin%lzd, lin%orbs, rxyz, denspot ,rhopotold, &
       nlpspd, proj, GPU, &
       tmb%psi)
  call mpi_barrier(mpi_comm_world, ierr)
  t2ig=mpi_wtime()
  timeig=t2ig-t1ig
  t1scc=mpi_wtime()
  !lphi=wfnmd%phi
  !call dcopy(lin%orbs%npsidim_orbs, wfnmd%phi(1), 1, lphi(1), 1)
  call dcopy(tmb%wfnmd%nphi, tmb%psi(1), 1, wfnmd%phi(1), 1)


  ! Initialize the DIIS mixing of the potential if required.
  if(lin%mixHist_lowaccuracy>0) then
     !ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*denspot%dpcom%nscatterarr(iproc,2)
      call initializeMixrhopotDIIS(lin%mixHist_lowaccuracy, denspot%dpcom%ndimpot, mixdiis)
  end if

  !end of the initialization part, will later be moved to cluster
  call timing(iproc,'INIT','PR')

  allocate(locrad(lin%lzd%nlr), stat=istat)
  call memocc(istat, locrad, 'locrad', subname)


  if(lin%nItInguess>0) then
      ! Post communications for gathering the potential.
     !ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*denspot%dpcom%nscatterarr(iproc,2)
      call allocateCommunicationsBuffersPotential(lin%comgp, subname)
      call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, lin%comgp)
      !!if(wfnmd%bs%use_derivative_basis) then
      if(tmb%wfnmd%bs%use_derivative_basis) then
          call allocateCommunicationsBuffersPotential(lin%lb%comgp, subname)
          call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, lin%lb%comgp)
      end if

      ! Calculate the Hamiltonian in the basis of the trace minimizing orbitals. Do not improve
      ! the basis functions (therefore update is set to false).
      ! This subroutine will also post the point to point messages needed for the calculation
      ! of the charge density.
      !wfnmd%bs%communicate_phi_for_lsumrho=.true.
      wfnmd%bs%communicate_phi_for_lsumrho=.true.
      tmb%wfnmd%bs%communicate_phi_for_lsumrho=.true.
      with_auxarray=.false.
      lin%newgradient=.false.
      wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE
      tmb%wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE

      if(lin%newgradient) then
          do ilr=1,lin%lzd%nlr
              !locrad(ilr)=lin%locrad_lowaccuracy(ilr)
              locrad(ilr)=lin%locrad_highaccuracy(ilr)
          end do
      else
          do ilr=1,lin%lzd%nlr
              !locrad(ilr)=lin%locrad_highaccuracy(ilr)
              locrad(ilr)=lin%locrad_lowaccuracy(ilr)
          end do
      end if

      if(lin%mixedmode) then
          call allocateCommunicationbufferSumrho(iproc, with_auxarray, lin%comsr, subname)
          wfnmd%bs%use_derivative_basis=.false.
          tmb%wfnmd%bs%use_derivative_basis=.false.
          call getLinearPsi(iproc, nproc, lin%lzd, orbs, lin%orbs, lin%orbs, lin%comsr, &
              lin%mad, lin%mad, lin%op, lin%op, lin%comon, lin%comon, &
              lin%comgp, lin%comgp, at, rxyz, &
              denspot, GPU, &
              infoBasisFunctions, infoCoeff, 0, ebs, nlpspd, proj, &
              ldiis, &
              orthpar, confdatarr, wfnmd%bpo%blocksize_pdgemm, &
              lin%lb%comrp, wfnmd%bpo%blocksize_pdsyev, wfnmd%bpo%nproc_pdsyev, &
              hx, hy, hz, input%SIC, locrad, wfnmd, tmb, tmbder)
      else
          call allocateCommunicationbufferSumrho(iproc,with_auxarray,lin%lb%comsr,subname)
          call getLinearPsi(iproc,nproc,lin%lzd,orbs,lin%orbs,lin%lb%orbs,lin%lb%comsr,&
              lin%mad,lin%lb%mad,lin%op,lin%lb%op,lin%comon,&
              lin%lb%comon,lin%comgp,lin%lb%comgp,at,rxyz,&
              denspot,GPU,&
              infoBasisFunctions,infoCoeff,0, ebs,nlpspd,proj,&
              ldiis,orthpar,confdatarr,& 
              wfnmd%bpo%blocksize_pdgemm,&
              lin%lb%comrp,wfnmd%bpo%blocksize_pdsyev,wfnmd%bpo%nproc_pdsyev,&
              hx,hy,hz,input%SIC, locrad, wfnmd, tmb, tmbder)
      end if
      !!call getLinearPsi(iproc, nproc, input%nspin, lin%lzd, orbs, lin%orbs, lin%lb%orbs, lin%lb%comsr, &
      !!    lin%op, lin%lb%op, lin%comon, lin%lb%comon, comms, at, lin, rxyz, rxyz, &
      !!    nscatterarr, ngatherarr, rhopot, GPU, input, pkernelseq, phi, updatePhi, &
      !!    infoBasisFunctions, infoCoeff, 0, n3p, n3pi, n3d, pkernel, &
      !!    i3s, i3xcsh, ebs, coeff, lphi, radii_cf, nlpspd, proj, wfnmd%bs%communicate_phi_for_lsumrho, coeff_proj)

      ! Calculate the charge density.
      !!call cpu_time(t1)
      if(lin%mixedmode) then
          call deallocateCommunicationbufferSumrho(lin%comsr, subname)
      else
          call deallocateCommunicationbufferSumrho(lin%lb%comsr, subname)
      end if
      !!call cpu_time(t2)
      !!time=t2-t1
      !!call mpiallred(time, 1, mpi_sum, mpi_comm_world, ierr)
      !!time=time/dble(nproc)
      !!if(iproc==0) write(*,'(1x,a,es12.4)') 'time for sumrho:',time

      if(trim(lin%mixingMethod)=='dens') then
          !if(lin%mixHist==0) then
          !    !if(n3p>0) call mixPotential(iproc, n3p, Glr, input, lin, rhopotOld, rhopot, pnrm)
          !    call mixPotential(iproc, n3p, Glr, input, lin%alphaMixWhenFixed, rhopotOld, rhopot, pnrm)
          !else 
          !    ndimpot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
          !    ndimtot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*lin%lzd%Glr%d%n3i
          !    mixdiis%mis=mod(mixdiis%is,mixdiis%isx)+1
          !    mixdiis%is=mixdiis%is+1
          !    call mixrhopotDIIS(iproc, nproc, ndimpot, rhopot, rhopotold, mixdiis, ndimtot, lin%alphaMixWhenFixed, 1, pnrm)
          !end if
          !rhopotold_out=rhopot
          rhopotold_out=rhopotold
      end if


      if(trim(lin%mixingMethod)=='pot') then
          if(lin%mixHist_lowaccuracy==0) then
              call mixPotential(iproc, denspot%dpcom%n3p, Glr, input, &
                   lin%alphaMixWhenFixed_lowaccuracy, rhopotOld, denspot%rhov, pnrm)
          else 
              !ndimpot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*denspot%dpcom%nscatterarr(iproc,2)
              ndimtot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*lin%lzd%Glr%d%n3i
              mixdiis%mis=mod(mixdiis%is,mixdiis%isx)+1
              mixdiis%is=mixdiis%is+1
              call mixrhopotDIIS(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, rhopotold, mixdiis, ndimtot, &
                   lin%alphaMixWhenFixed_lowaccuracy, 2, pnrm)
          end if
          rhopotold_out=denspot%rhov
      end if

      ! Copy the current potential
      if(trim(lin%mixingMethod)=='pot') then
           call dcopy(max(Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotOld(1), 1)
      end if
  end if

  !!! If we mix the potential, copy the potential.
  !!if(trim(lin%mixingMethod)=='pot') then
  !!    call dcopy(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin, rhopot(1), 1, rhopotOld(1), 1)
  !!end if

  ! Allocate the communications buffers needed for the communications of teh potential and
  ! post the messages. This will send to each process the part of the potential that this process
  ! needs for the application of the Hamlitonian to all orbitals on that process.
  !ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
  call allocateCommunicationsBuffersPotential(lin%comgp, subname)
  call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, lin%comgp)
  ! If we also use the derivative of the basis functions, also send the potential in this case. This is
  ! needed since the orbitals may be partitioned in a different way when the derivatives are used.
  if(wfnmd%bs%use_derivative_basis) then
  !!if(tmb%wfnmd%bs%use_derivative_basis) then
      call allocateCommunicationsBuffersPotential(lin%lb%comgp, subname)
      call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, lin%lb%comgp)
  end if




  !if(nproc==1) allocate(psit(size(psi)))
  nitSCC=lin%nitSCCWhenOptimizing+lin%nitSCCWhenFixed
  ! Flag that indicates that the basis functions shall be improved in the following.
  wfnmd%bs%update_phi=.true.
  tmb%wfnmd%bs%update_phi=.true.
  pnrm=1.d100
  pnrm_out=1.d100
  energyold=0.d0
  energyoldout=0.d0
  reduceConvergenceTolerance=.false.
  lin%newgradient=.false.
  wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE
  tmb%wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE
  lowaccur_converged=.false.

  outerLoop: do itout=1,lin%nit_lowaccuracy+lin%nit_highaccuracy

      ! First to some initialization and determine the value of some control parameters.

      ! Initialize DIIS...
      call initializeDIIS(lin%DIISHistMax, lin%lzd, lin%orbs, lin%orbs%norb, ldiis)
      ldiis%DIISHistMin=lin%DIISHistMin
      ldiis%DIISHistMax=lin%DIISHistMax
      ldiis%alphaSD=lin%alphaSD
      ldiis%alphaDIIS=lin%alphaDIIS

      ! The basis functions shall be optimized
      wfnmd%bs%update_phi=.true.
      tmb%wfnmd%bs%update_phi=.true.

      ! Convergence criterion for the self consistency looo
      selfConsistent=lin%convCritMix

      ! Check whether the derivatives shall be used or not.
      if(lin%mixedmode) then
          if( (.not.lowaccur_converged .and. (itout==lin%nit_lowaccuracy+1 .or. pnrm_out<lin%lowaccuray_converged) ) &
              .or. lowaccur_converged ) then
              withder=.true.
          else
              withder=.false.
          end if
      end if

      ! Check whether the low accuracy part (i.e. with strong confining potential) has converged.
      if(.not.lowaccur_converged .and. (itout==lin%nit_lowaccuracy+1 .or. pnrm_out<lin%lowaccuray_converged)) then
          lowaccur_converged=.true.
          nit_highaccuracy=0
      end if 

      ! Set all remaining variables that we need for the optimizations of the basis functions and the mixing.
      call set_optimization_variables(lowaccur_converged, input, at, lin%orbs, lin%lzd%nlr, onwhichatom, confdatarr, wfnmd, &
           locrad, nitSCC, nitSCCWhenOptimizing, mixHist, alphaMix)
      call set_optimization_variables(lowaccur_converged, input, at, lin%orbs, lin%lzd%nlr, onwhichatom, confdatarr, tmb%wfnmd, &
           locrad, nitSCC, nitSCCWhenOptimizing, mixHist, alphaMix)

      !!if(wfnmd%bs%confinement_decrease_mode==DECREASE_ABRUPT) then
      if(tmb%wfnmd%bs%confinement_decrease_mode==DECREASE_ABRUPT) then
          tt=1.d0
      !!else if(wfnmd%bs%confinement_decrease_mode==DECREASE_LINEAR) then
      else if(tmb%wfnmd%bs%confinement_decrease_mode==DECREASE_LINEAR) then
          tt=1.d0-(dble(itout-1))/dble(lin%nit_lowaccuracy)
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
          !!if(input%lin%nItInnerLoop/=-1 .or. wfnmd%bs%locreg_enlargement/=1.d0) then
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
      call allocateCommunicationbufferSumrho(iproc, with_auxarray, lin%comsr, subname)
      call allocateCommunicationbufferSumrho(iproc, with_auxarray, lin%lb%comsr, subname)

      ! Now all initializations are done...

      ! The self consistency cycle. Here we try to get a self consistent density/potential.
      ! In the first nitSCCWhenOptimizing iteration, the basis functions are optimized, whereas in the remaining
      ! iteration the basis functions are fixed.
      do itSCC=1,nitSCC
          if(itSCC>nitSCCWhenOptimizing) wfnmd%bs%update_phi=.false.
          if(itSCC>nitSCCWhenOptimizing) tmb%wfnmd%bs%update_phi=.false.
          if(itSCC==1) then
              wfnmd%bs%communicate_phi_for_lsumrho=.true.
              tmb%wfnmd%bs%communicate_phi_for_lsumrho=.true.
          else
              wfnmd%bs%communicate_phi_for_lsumrho=.false.
              tmb%wfnmd%bs%communicate_phi_for_lsumrho=.false.
          end if
          !!write(*,*) 'ATTENTION DEBUG'
          !!wfnmd%bs%communicate_phi_for_lsumrho=.true.

          ! Update the basis functions (if wfnmd%bs%update_phi is true), calculate the Hamiltonian in this basis, and diagonalize it.
          if(lin%mixedmode) then
              if(.not.withder) then
                  wfnmd%bs%use_derivative_basis=.false.
                  tmb%wfnmd%bs%use_derivative_basis=.false.
                  call getLinearPsi(iproc,nproc,lin%lzd,orbs,lin%orbs,lin%orbs,lin%comsr,&
                      lin%mad,lin%mad,lin%op,lin%op,lin%comon,&
                      lin%comon,lin%comgp,lin%comgp,at,rxyz,&
                      denspot,GPU,&
                      infoBasisFunctions,infoCoeff,itScc,ebs,nlpspd,proj,&
                      ldiis,orthpar,confdatarr,&
                      wfnmd%bpo%blocksize_pdgemm,&
                      lin%lb%comrp,wfnmd%bpo%blocksize_pdsyev,wfnmd%bpo%nproc_pdsyev,&
                      hx,hy,hz,input%SIC, locrad, wfnmd, tmb, tmbder)
              else
                  wfnmd%bs%use_derivative_basis=.true.
                  tmb%wfnmd%bs%use_derivative_basis=.true.
                  call getLinearPsi(iproc,nproc,lin%lzd,orbs,lin%orbs,lin%lb%orbs,lin%lb%comsr,&
                      lin%mad,lin%lb%mad,lin%op,lin%lb%op,&
                      lin%comon,lin%lb%comon,lin%comgp,lin%lb%comgp,at,rxyz,&
                      denspot,GPU,&
                      infoBasisFunctions,infoCoeff,itScc,ebs,nlpspd,proj,&
                      ldiis,orthpar,confdatarr,&
                      wfnmd%bpo%blocksize_pdgemm,&
                      lin%lb%comrp,wfnmd%bpo%blocksize_pdsyev,wfnmd%bpo%nproc_pdsyev,&
                      hx,hy,hz,input%SIC, locrad, wfnmd, tmb, tmbder)
              end if
          else
              call getLinearPsi(iproc,nproc,lin%lzd,orbs,lin%orbs,lin%lb%orbs,lin%lb%comsr,&
                  lin%mad,lin%lb%mad,lin%op,lin%lb%op,lin%comon,&
                  lin%lb%comon,lin%comgp,lin%lb%comgp,at,rxyz,&
                  denspot,GPU,&
                  infoBasisFunctions,infoCoeff,itScc,ebs,nlpspd,proj,&
                  ldiis,orthpar,confdatarr,&
                  wfnmd%bpo%blocksize_pdgemm,&
                  lin%lb%comrp,wfnmd%bpo%blocksize_pdsyev,wfnmd%bpo%nproc_pdsyev,&
                  hx,hy,hz,input%SIC, locrad, wfnmd, tmb, tmbder)
          end if


          ! Calculate the charge density.
          if(lin%mixedmode) then
              if(.not.withder) then
                  call sumrhoForLocalizedBasis2(iproc, nproc, orbs%norb, &
                       lin%lzd, input, hx, hy, hz, lin%orbs, lin%comsr, &
                       wfnmd%ld_coeff, tmbder%wfnmd%coeff, Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3d, &
                       denspot%rhov, at, denspot%dpcom%nscatterarr)
               else
                  call sumrhoForLocalizedBasis2(iproc, nproc, orbs%norb,&
                       lin%lzd, input, hx, hy, hz, lin%lb%orbs, lin%lb%comsr, &
                       wfnmd%ld_coeff, tmbder%wfnmd%coeff, Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3d,&
                       denspot%rhov, at, denspot%dpcom%nscatterarr)
               end if
          else
              call sumrhoForLocalizedBasis2(iproc, nproc, orbs%norb,&
                   lin%lzd, input, hx, hy ,hz, lin%lb%orbs, lin%lb%comsr, &
                   wfnmd%ld_coeff, tmbder%wfnmd%coeff, Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3d, &
                   denspot%rhov, at, denspot%dpcom%nscatterarr)
          end if

          ! Mix the density.
          if(trim(lin%mixingMethod)=='dens') then
              if(mixHist==0) then
                  call mixPotential(iproc, denspot%dpcom%n3p, Glr, input, alphaMix, rhopotOld, denspot%rhov, pnrm)
              else 
                 !ndimpot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
                  ndimtot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*lin%lzd%Glr%d%n3i
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
          if(trim(lin%mixingMethod)=='dens') then
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
          if(trim(lin%mixingMethod)=='pot') then
              if(mixHist==0) then
                  call mixPotential(iproc, denspot%dpcom%n3p, Glr, input, alphaMix, rhopotOld, denspot%rhov, pnrm)
              else 
                 !ndimpot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
                  ndimtot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*lin%lzd%Glr%d%n3i
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
          if(trim(lin%mixingMethod)=='pot') then
               call dcopy(max(Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotOld(1), 1)
          end if

          ! Post communications for gathering the potential
          call allocateCommunicationsBuffersPotential(lin%comgp, subname)
          call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, lin%comgp)
          !if(wfnmd%bs%use_derivative_basis) then
          if(tmb%wfnmd%bs%use_derivative_basis) then
              call allocateCommunicationsBuffersPotential(lin%lb%comgp, subname)
              call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, lin%lb%comgp)
          end if

          ! Write some informations.
          call printSummary(iproc, itSCC, infoBasisFunctions, &
               infoCoeff, pnrm, energy, energyDiff, lin%mixingMethod)
          if(pnrm<selfConsistent) then
              reduceConvergenceTolerance=.true.
              exit
          else
              reduceConvergenceTolerance=.false.
          end if
      end do

      call deallocateCommunicationbufferSumrho(lin%comsr, subname)
      call deallocateCommunicationbufferSumrho(lin%lb%comsr, subname)

      ! Print out values related to two iterations of the outer loop.
      if(iproc==0) then
          write(*,'(3x,a,7es18.10)') 'ebs, ehart, eexcu, vexcu, eexctX, eion, edisp', &
              ebs, ehart, eexcu, vexcu, eexctX, eion, edisp
          if(trim(lin%mixingMethod)=='dens') then
              write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                   'itout, Delta DENSOUT, energy, energyDiff', itout, pnrm_out, energy, energy-energyoldout
          else if(trim(lin%mixingMethod)=='pot') then
              write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                   'itout, Delta POTOUT, energy energyDiff', itout, pnrm_out, energy, energy-energyoldout
          end if
      end if
      !!if(abs(pnrm_out)<lin%convCritMixOut) exit
      energyoldout=energy

      ! Deallocate DIIS structures.
      call deallocateDIIS(ldiis)

  end do outerLoop


  call cancelCommunicationPotential(iproc, nproc, lin%comgp)
  call deallocateCommunicationsBuffersPotential(lin%comgp, subname)
  if(wfnmd%bs%use_derivative_basis) then
  !!if(tmb%wfnmd%bs%use_derivative_basis) then
      call cancelCommunicationPotential(iproc, nproc, lin%lb%comgp)
      call deallocateCommunicationsBuffersPotential(lin%lb%comgp, subname)
  end if

  iall=-product(shape(rhopotOld))*kind(rhopotOld)
  deallocate(rhopotOld, stat=istat)
  call memocc(istat, iall, 'rhopotold', subname)
  iall=-product(shape(rhopotold_out))*kind(rhopotold_out)
  deallocate(rhopotold_out, stat=istat)
  call memocc(istat, iall, 'rhopotold_out', subname)
  iall=-product(shape(onwhichatom))*kind(onwhichatom)
  deallocate(onwhichatom, stat=istat)
  call memocc(istat, iall, 'onwhichatom', subname)

  if(lin%mixHist_highaccuracy>0) then
      call deallocateMixrhopotDIIS(mixdiis)
  end if

  call mpi_barrier(mpi_comm_world, ierr)
  t2scc=mpi_wtime()
  timescc=t2scc-t1scc


  ! Allocate the communication buffers for the calculation of the charge density.
  with_auxarray=.false.
  call allocateCommunicationbufferSumrho(iproc, with_auxarray, lin%lb%comsr, subname)
  call communicate_basis_for_density(iproc, nproc, lin%lzd, lin%lb%orbs, wfnmd%phi, lin%lb%comsr)
  call sumrhoForLocalizedBasis2(iproc, nproc, orbs%norb, lin%lzd, input, hx, hy, hz, lin%lb%orbs, lin%lb%comsr, &
       wfnmd%ld_coeff, tmbder%wfnmd%coeff, Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3d, denspot%rhov, at,denspot%dpcom%nscatterarr)

  call deallocateCommunicationbufferSumrho(lin%lb%comsr, subname)

  call mpi_barrier(mpi_comm_world, ierr)
  t1force=mpi_wtime()
  ! Build global orbitals psi (the physical ones).
  if(lin%transformToGlobal) then
      call transformToGlobal(iproc, nproc, lin, orbs, comms, input, tmbder%wfnmd%coeff, wfnmd%phi, psi, psit)
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
  !!!call density_and_hpot(iproc,nproc,at%geocode,at%sym,orbs,lin%Lzd,&
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



  ! Deallocate all arrays related to the linear scaling version.
  !!call deallocateLinear(iproc, lin, lphi, coeff)
  call deallocate_linearParameters(lin, subname)

  call destroy_wfn_metadata(wfnmd)
  call destroy_DFT_wavefunction(tmb)
  call destroy_DFT_wavefunction(tmbder)

  !call deallocateBasicArraysInput(at, input%lin)
  call deallocateBasicArraysInput(input%lin)

  deallocate(confdatarr)
  call deallocateBasicArrays(lin)

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


subroutine transformToGlobal(iproc,nproc,lin,orbs,comms,input,coeff,lphi,psi,psit)
use module_base
use module_types
use module_interfaces, exceptThisOne => transformToGlobal
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linearParameters),intent(in):: lin
type(orbitals_data),intent(in):: orbs
type(communications_arrays):: comms
type(input_variables),intent(in):: input
real(8),dimension(lin%lb%orbs%norb,orbs%norb),intent(in):: coeff
!real(8),dimension(max(lin%orbs%npsidim_orbs,lin%orbs%npsidim_comp)),intent(inout):: lphi
real(8),dimension(*),intent(inout):: lphi
!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: psi, psit
real(8),dimension(:),pointer,intent(out):: psi, psit

! Local variables
integer:: ind1, ind2, istat, iall, iorb, ilr, ldim, gdim, nvctrp
real(8),dimension(:),pointer:: phiWork
real(8),dimension(:),allocatable:: phi
character(len=*),parameter:: subname='transformToGlobal'

  !do iall=0,nproc-1
  !    write(*,'(a,i5,4i12)') 'START transformToGlobal: iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)', iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)  
  !end do

  allocate(phi(max(lin%lb%gorbs%npsidim_orbs,lin%lb%gorbs%npsidim_comp)+ndebug), stat=istat)
  call memocc(istat, phi, 'phi', subname)
  allocate(phiWork(max(size(phi),size(psi))), stat=istat)
  call memocc(istat, phiWork, 'phiWork', subname)

  ind1=1
  ind2=1
!  phi=0.d0
  if (max(lin%lb%gorbs%npsidim_orbs,lin%lb%gorbs%npsidim_comp) > 0) &
       call to_zero(max(lin%lb%gorbs%npsidim_orbs,lin%lb%gorbs%npsidim_comp),phi(1))

  do iorb=1,lin%lb%orbs%norbp
      !ilr = lin%lb%orbs%inWhichLocregp(iorb)
      ilr = lin%lb%orbs%inWhichLocreg(lin%lb%orbs%isorb+iorb)
      ldim=lin%lzd%Llr(ilr)%wfd%nvctr_c+7*lin%lzd%Llr(ilr)%wfd%nvctr_f
      gdim=lin%lzd%Glr%wfd%nvctr_c+7*lin%lzd%Glr%wfd%nvctr_f
      call Lpsi_to_global2(iproc,nproc,ldim,gdim,lin%lb%orbs%norb,lin%lb%orbs%nspinor,input%nspin,lin%lzd%Glr,&
           lin%lzd%Llr(ilr),lphi(ind2),phi(ind1))
      ind1=ind1+lin%lzd%Glr%wfd%nvctr_c+7*lin%lzd%Glr%wfd%nvctr_f
      ind2=ind2+lin%lzd%Llr(ilr)%wfd%nvctr_c+7*lin%lzd%Llr(ilr)%wfd%nvctr_f
  end do
  !if(ind1/=lin%gorbs%npsidim+1) then
  !    write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': ind1/=lin%gorbs%npsidim',ind1,lin%gorbs%npsidim
  !end if
  !if(ind2/=lin%lb%orbs%npsidim+1) then
  !    write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': ind1/=lin%lb%orbs%npsidim',ind1,lin%lb%orbs%npsidim
  !end if
  !do iall=0,nproc-1
  !    write(*,'(a,i5,4i12)') 'after loop: iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)', iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)  
  !end do
  call transpose_v(iproc, nproc, lin%lb%orbs, lin%lzd%Glr%wfd, lin%lb%comms, phi, work=phiWork)
  !do iall=0,nproc-1
  !    write(*,'(a,i5,4i12)') 'after transpose phi: iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)', iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)  
  !end do


  if(iproc==0) then
      write(*,'(1x,a)', advance='no') '------------------------------------- Building linear combinations... '
  end if
  ! Build the extended orbital psi as a linear combination of localized basis functions phi. for real O(N)
  ! this has to replaced, but at the moment it is still needed.
  !call buildWavefunctionModified(iproc, nproc, orbs, lin%lb%gorbs, comms, lin%lb%gcomms, phi, psi, coeff)
  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  !write(*,*) 'iproc, nvctrp', iproc, nvctrp
  !write(*,*) 'iproc, orbs%npsidim', iproc, orbs%npsidim
  call dgemm('n', 'n', nvctrp, orbs%norb, lin%lb%orbs%norb, 1.d0, phi(1), nvctrp, coeff(1,1), &
       lin%lb%orbs%norb, 0.d0, psi(1), nvctrp)

  !do iall=0,nproc-1
  !    write(*,'(a,i5,4i12)') 'after buildWavefunctionModified: iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)', iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)  
  !end do


  if(nproc>1) then
      call dcopy(orbs%npsidim_comp, psi, 1, psit, 1)
  else
      psit => psi
  end if

  call untranspose_v(iproc, nproc, lin%lb%orbs, lin%lzd%Glr%wfd, lin%lb%comms, phi, work=phiWork)
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
  call untranspose_v(iproc, nproc, orbs, lin%lzd%Glr%wfd, comms, psi, work=phiWork)

  if(iproc==0) write(*,'(a)') 'done.'


  iall=-product(shape(phi))*kind(phi)
  deallocate(phi, stat=istat)
  call memocc(istat, iall, 'phi', subname)
  iall=-product(shape(phiWork))*kind(phiWork)
  deallocate(phiWork, stat=istat)
  call memocc(istat, iall, 'phiWork', subname)

end subroutine transformToGlobal


subroutine create_wfn_metadata(mode, nphi, nlbphi, lnorb, llbnorb, norb, input, wfnmd)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  character(len=1),intent(in):: mode
  integer,intent(in):: nphi, nlbphi, lnorb, llbnorb, norb
  type(input_variables),intent(in):: input
  type(wfn_metadata),intent(out):: wfnmd

  ! Local variables
  integer:: istat
  character(len=*),parameter:: subname='create_wfn_metadata'

  ! Determine which variables we need, depending on the mode we are in.
  if(mode=='l') then
      ! linear scaling mode

      wfnmd%nphi=nphi
      wfnmd%nlbphi=nlbphi
      wfnmd%basis_is=BASIS_IS_ENHANCED !since always it is allocated with wfnmd%nlbphi
      wfnmd%ld_coeff=llbnorb !leading dimension of the coeff array

      allocate(wfnmd%phi(wfnmd%nlbphi), stat=istat)
      call memocc(istat, wfnmd%phi, 'wfnmd%phi', subname)

      allocate(wfnmd%phiRestart(wfnmd%nphi), stat=istat)
      call memocc(istat, wfnmd%phiRestart, 'wfnmd%phiRestart', subname)

      allocate(wfnmd%coeff(llbnorb,norb), stat=istat)
      call memocc(istat, wfnmd%coeff, 'wfnmd%coeff', subname)

      !!allocate(wfnmd%coeff_proj(lnorb,norb), stat=istat)
      !!call memocc(istat, wfnmd%coeff_proj, 'wfnmd%coeff_proj', subname)

      call init_basis_specifications(input, wfnmd%bs)
      call init_basis_performance_options(input, wfnmd%bpo)

  else if(mode=='c') then
      ! cubic scaling mode

      nullify(wfnmd%phi)
      nullify(wfnmd%phiRestart)
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

  iall=-product(shape(wfnmd%phi))*kind(wfnmd%phi)
  deallocate(wfnmd%phi, stat=istat)
  call memocc(istat, iall, 'wfnmd%phi', subname)

  iall=-product(shape(wfnmd%phiRestart))*kind(wfnmd%phiRestart)
  deallocate(wfnmd%phiRestart, stat=istat)
  call memocc(istat, iall, 'wfnmd%phiRestart', subname)

  iall=-product(shape(wfnmd%coeff))*kind(wfnmd%coeff)
  deallocate(wfnmd%coeff, stat=istat)
  call memocc(istat, iall, 'wfnmd%coeff', subname)

  !!iall=-product(shape(wfnmd%coeff_proj))*kind(wfnmd%coeff_proj)
  !!deallocate(wfnmd%coeff_proj, stat=istat)
  !!call memocc(istat, iall, 'wfnmd%coeff_proj', subname)

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

  call create_wfn_metadata(mode, nphi, nphi, lnorb, lnorb, norb, input, wfn%wfnmd)

  allocate(wfn%psi(wfn%wfnmd%nphi), stat=istat)
  call memocc(istat, wfn%psi, 'wfn%psi', subname)

end subroutine create_DFT_wavefunction



subroutine destroy_DFT_wavefunction(wfn)
  use module_base
  use module_types
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

  call destroy_wfn_metadata(wfn%wfnmd)

end subroutine destroy_DFT_wavefunction

