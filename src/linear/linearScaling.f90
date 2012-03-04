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
real(8),dimension(:,:),pointer:: coeff
real(8):: ebs, ebsMod, pnrm, tt, ehart, eexcu, vexcu, alphaMix
character(len=*),parameter:: subname='linearScaling'
real(8),dimension(:),allocatable:: rhopotOld, rhopotold_out
logical:: updatePhi, reduceConvergenceTolerance, communicate_lphi, with_auxarray, lowaccur_converged, withder
real(8),dimension(:),pointer:: lphi
real(8):: t1, t2, time, t1tot, t2tot, timetot, t1ig, t2ig, timeig, t1init, t2init, timeinit, ddot, dnrm2, pnrm_out
real(8):: t1scc, t2scc, timescc, t1force, t2force, timeforce, energyold, energyDiff, energyoldout, selfConsistent
integer:: iorb, ndimtot, iiat
type(mixrhopotDIISParameters):: mixdiis
type(workarr_sumrho):: w
real(8),dimension(:,:),allocatable:: coeff_proj
type(localizedDIISParameters):: ldiis
type(confpot_data), dimension(:),pointer :: confdatarr
real(8):: fnoise,pressure
real(gp), dimension(6) :: ewaldstr,strten,hstrten,xcstr
type(orthon_data):: orthpar
integer,dimension(:),pointer:: onwhichatom

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
       input, hx, hy, hz, rxyz, denspot%dpcom%nscatterarr, tag, coeff, lphi, confdatarr, onwhichatom)

  !!lin%potentialPrefac=lin%potentialPrefac_lowaccuracy
  !!allocate(confdatarr(lin%orbs%norbp))
  !!!use a temporary array onwhichatom instead of inwhichlocreg
  !!
  !!call define_confinement_data(confdatarr,lin%orbs,rxyz,at,&
  !!     hx,hy,hz,lin,lin%lzd,lin%orbs%inWhichLocreg)


  orthpar%methTransformOverlap = lin%methTransformOverlap
  orthpar%nItOrtho = lin%nItOrtho
  orthpar%blocksize_pdsyev = lin%blocksize_pdsyev
  orthpar%blocksize_pdgemm = lin%blocksize_pdgemm


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

  allocate(coeff_proj(lin%orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, coeff_proj, 'coeff_proj', subname)

  !write(*,'(a,100i6)') 'lin%orbs%inwhichlocreg', lin%orbs%inwhichlocreg


  potshortcut=0 ! What is this?
  call mpi_barrier(mpi_comm_world, ierr)
  t1ig=mpi_wtime()
  call inputguessConfinement(iproc, nproc, at, &
       input, hx, hy, hz, lin%lzd, lin%orbs, rxyz, denspot ,rhopotold, &
       nlpspd, proj, GPU, &
       lphi)
  call mpi_barrier(mpi_comm_world, ierr)
  t2ig=mpi_wtime()
  timeig=t2ig-t1ig
  t1scc=mpi_wtime()

  call deallocateBasicArraysInput(at, input%lin)

  ! Initialize the DIIS mixing of the potential if required.
  if(lin%mixHist_lowaccuracy>0) then
     !ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*denspot%dpcom%nscatterarr(iproc,2)
      call initializeMixrhopotDIIS(lin%mixHist_lowaccuracy, denspot%dpcom%ndimpot, mixdiis)
  end if

  !end of the initialization part, will later be moved to cluster
  call timing(iproc,'INIT','PR')

  if(lin%nItInguess>0) then
      ! Post communications for gathering the potential.
     !ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*denspot%dpcom%nscatterarr(iproc,2)
      call allocateCommunicationsBuffersPotential(lin%comgp, subname)
      call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, lin%comgp)
      if(lin%useDerivativeBasisFunctions) then
          call allocateCommunicationsBuffersPotential(lin%lb%comgp, subname)
          call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, lin%lb%comgp)
      end if

      ! Calculate the Hamiltonian in the basis of the trace minimizing orbitals. Do not improve
      ! the basis functions (therefore update is set to false).
      ! This subroutine will also post the point to point messages needed for the calculation
      ! of the charge density.
      updatePhi=.false.
      !communicate_lphi=.true.
      communicate_lphi=.true.
      with_auxarray=.false.
      lin%newgradient=.false.
      if(lin%mixedmode) then
          call allocateCommunicationbufferSumrho(iproc, with_auxarray, lin%comsr, subname)
          lin%useDerivativeBasisFunctions=.false.
          call getLinearPsi(iproc, nproc, lin%lzd, orbs, lin%orbs, lin%orbs, lin%comsr, &
              lin%mad, lin%mad, lin%op, lin%op, lin%comon, lin%comon, &
              lin%comgp, lin%comgp, at, rxyz, &
              denspot, GPU, updatePhi, &
              infoBasisFunctions, infoCoeff, 0, ebs, coeff, lphi, nlpspd, proj, &
              communicate_lphi, coeff_proj, ldiis, nit, lin%nItInnerLoop, &
              lin%newgradient, orthpar, confdatarr, lin%methTransformOverlap, lin%blocksize_pdgemm, &
              lin%convCrit, lin%nItPrecond, lin%useDerivativeBasisFunctions, lin%lphiRestart, &
              lin%lb%comrp, lin%blocksize_pdsyev, lin%nproc_pdsyev, &
              hx, hy, hz, input%SIC, input%lin%factor_enlarge)
      else
          call allocateCommunicationbufferSumrho(iproc,with_auxarray,lin%lb%comsr,subname)
          call getLinearPsi(iproc,nproc,lin%lzd,orbs,lin%orbs,lin%lb%orbs,lin%lb%comsr,&
              lin%mad,lin%lb%mad,lin%op,lin%lb%op,lin%comon,&
              lin%lb%comon,lin%comgp,lin%lb%comgp,at,rxyz,&
              denspot,GPU,updatePhi,&
              infoBasisFunctions,infoCoeff,0, ebs,coeff,lphi,nlpspd,proj,communicate_lphi,&
              coeff_proj,ldiis,nit,lin%nItInnerLoop,lin%newgradient,orthpar,confdatarr,& 
              lin%methTransformOverlap,lin%blocksize_pdgemm,lin%convCrit,lin%nItPrecond,&
              lin%useDerivativeBasisFunctions,lin%lphiRestart,lin%lb%comrp,lin%blocksize_pdsyev,lin%nproc_pdsyev,&
              hx,hy,hz,input%SIC, input%lin%factor_enlarge)
      end if
      !!call getLinearPsi(iproc, nproc, input%nspin, lin%lzd, orbs, lin%orbs, lin%lb%orbs, lin%lb%comsr, &
      !!    lin%op, lin%lb%op, lin%comon, lin%lb%comon, comms, at, lin, rxyz, rxyz, &
      !!    nscatterarr, ngatherarr, rhopot, GPU, input, pkernelseq, phi, updatePhi, &
      !!    infoBasisFunctions, infoCoeff, 0, n3p, n3pi, n3d, pkernel, &
      !!    i3s, i3xcsh, ebs, coeff, lphi, radii_cf, nlpspd, proj, communicate_lphi, coeff_proj)

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
  if(lin%useDerivativeBasisFunctions) then
      call allocateCommunicationsBuffersPotential(lin%lb%comgp, subname)
      call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, lin%lb%comgp)
  end if




  !if(nproc==1) allocate(psit(size(psi)))
  nitSCC=lin%nitSCCWhenOptimizing+lin%nitSCCWhenFixed
  ! Flag that indicates that the basis functions shall be improved in the following.
  updatePhi=.true.
  pnrm=1.d100
  pnrm_out=1.d100
  energyold=0.d0
  energyoldout=0.d0
  reduceConvergenceTolerance=.false.
  lin%newgradient=.false.
  lowaccur_converged=.false.

  outerLoop: do itout=1,lin%nit_lowaccuracy+lin%nit_highaccuracy


      call initializeDIIS(lin%DIISHistMax, lin%lzd, lin%orbs, lin%orbs%norb, ldiis)
      !!call initializeDIIS(lin%DIISHistMax, lin%lzdlarge, lin%orbslarge, lin%orbslarge%norb, ldiis)
      ldiis%DIISHistMin=lin%DIISHistMin
      ldiis%DIISHistMax=lin%DIISHistMax
      ldiis%alphaSD=lin%alphaSD
      ldiis%alphaDIIS=lin%alphaDIIS

      updatePhi=.true.
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

      ! Choose the correct confining potential and gradient method, depending on whether we are in the low accuracy
      ! or high accuracy part.
      if(lowaccur_converged) then
          !!lin%potentialPrefac = lin%potentialPrefac_highaccuracy
          do iorb=1,lin%orbs%norbp
              ilr=lin%orbs%inwhichlocreg(lin%orbs%isorb+iorb)
              iiat=onwhichatom(lin%orbs%isorb+iorb)
              !confdatarr(iorb)%prefac=lin%potentialPrefac_highaccuracy(at%iatype(ilr))
              confdatarr(iorb)%prefac=lin%potentialPrefac_highaccuracy(at%iatype(iiat))
          end do
          lin%newgradient=.true.
          nit_highaccuracy=nit_highaccuracy+1
          nit=lin%nItBasis_highaccuracy
          if(nit_highaccuracy==lin%nit_highaccuracy+1) then
              ! Deallocate DIIS structures.
              call deallocateDIIS(ldiis)
              exit outerLoop
          end if
          ! only use steepest descent if the localization regions may change
          if(lin%nItInnerLoop/=-1 .or. input%lin%factor_enlarge/=1.d0) then
              ldiis%isx=0
          end if

      else
          !!lin%potentialPrefac = lin%potentialPrefac_lowaccuracy
          do iorb=1,lin%orbs%norbp
              ilr=lin%orbs%inwhichlocreg(lin%orbs%isorb+iorb)
              iiat=onwhichatom(lin%orbs%isorb+iorb)
              !confdatarr(iorb)%prefac=lin%potentialPrefac_lowaccuracy(at%iatype(ilr))
              confdatarr(iorb)%prefac=lin%potentialPrefac_lowaccuracy(at%iatype(iiat))
          end do
          lin%newgradient=.false.
          nit=lin%nItBasis_lowaccuracy
      end if

      ! Allocate the communication arrays for the calculation of the charge density.
      with_auxarray=.false.
      call allocateCommunicationbufferSumrho(iproc, with_auxarray, lin%comsr, subname)
      call allocateCommunicationbufferSumrho(iproc, with_auxarray, lin%lb%comsr, subname)

      ! Optimize the basis functions and them mix the density / potential to reach self consistency.
      if(lowaccur_converged) then
          nitSCC=lin%nitSCCWhenOptimizing_lowaccuracy+lin%nitSCCWhenFixed_lowaccuracy
          nitSCCWhenOptimizing=lin%nitSCCWhenOptimizing_lowaccuracy
          mixHist=lin%mixHist_lowaccuracy
          if(lin%mixHist_lowaccuracy==0 .and. lin%mixHist_highaccuracy>0) then
             !ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
              call initializeMixrhopotDIIS(lin%mixHist_highaccuracy, denspot%dpcom%ndimpot, mixdiis)
          else if(lin%mixHist_lowaccuracy>0 .and. lin%mixHist_highaccuracy==0) then
              call deallocateMixrhopotDIIS(mixdiis)
          end if
      else
          nitSCC=lin%nitSCCWhenOptimizing_highaccuracy+lin%nitSCCWhenFixed_highaccuracy
          nitSCCWhenOptimizing=lin%nitSCCWhenOptimizing_highaccuracy
          mixHist=lin%mixHist_highaccuracy
      end if

      ! The self consistency cycle. Here we try to get a self consistent density/potential.
      ! In the first nitSCCWhenOptimizing iteration, the basis functions are optimized, whereas in the remaining
      ! iteration the basis functions are fixed.
      do itSCC=1,nitSCC
          if(itSCC>nitSCCWhenOptimizing) updatePhi=.false.
          if(itSCC==1) then
              communicate_lphi=.true.
          else
              communicate_lphi=.false.
          end if

          ! Update the basis functions (if updatePhi is true), diagonalize the Hamiltonian in this basis, and diagonalize it.
          if(lin%mixedmode) then
              if(.not.withder) then
                  lin%useDerivativeBasisFunctions=.false.
                  call getLinearPsi(iproc,nproc,lin%lzd,orbs,lin%orbs,lin%orbs,lin%comsr,&
                      lin%mad,lin%mad,lin%op,lin%op,lin%comon,&
                      lin%comon,lin%comgp,lin%comgp,at,rxyz,&
                      denspot,GPU,updatePhi,&
                      infoBasisFunctions,infoCoeff,itScc,ebs,coeff,lphi,nlpspd,proj,communicate_lphi,&
                      coeff_proj,ldiis,nit,lin%nItInnerLoop,lin%newgradient,orthpar,confdatarr,&
                      lin%methTransformOverlap,lin%blocksize_pdgemm,lin%convCrit,lin%nItPrecond,&
                      lin%useDerivativeBasisFunctions,lin%lphiRestart,lin%lb%comrp,lin%blocksize_pdsyev,lin%nproc_pdsyev,&
                      hx,hy,hz,input%SIC, input%lin%factor_enlarge)
              else
                  lin%useDerivativeBasisFunctions=.true.
                  call getLinearPsi(iproc,nproc,lin%lzd,orbs,lin%orbs,lin%lb%orbs,lin%lb%comsr,&
                      lin%mad,lin%lb%mad,lin%op,lin%lb%op,&
                      lin%comon,lin%lb%comon,lin%comgp,lin%lb%comgp,at,rxyz,&
                      denspot,GPU,updatePhi,&
                      infoBasisFunctions,infoCoeff,itScc,ebs,coeff,lphi,nlpspd,proj,communicate_lphi,&
                      coeff_proj,ldiis,nit,lin%nItInnerLoop,lin%newgradient,orthpar,confdatarr,&
                      lin%methTransformOverlap,lin%blocksize_pdgemm,lin%convCrit,lin%nItPrecond,&
                      lin%useDerivativeBasisFunctions,lin%lphiRestart,lin%lb%comrp,lin%blocksize_pdsyev,lin%nproc_pdsyev,&
                      hx,hy,hz,input%SIC, input%lin%factor_enlarge)
              end if
          else
              call getLinearPsi(iproc,nproc,lin%lzd,orbs,lin%orbs,lin%lb%orbs,lin%lb%comsr,&
                  lin%mad,lin%lb%mad,lin%op,lin%lb%op,lin%comon,&
                  lin%lb%comon,lin%comgp,lin%lb%comgp,at,rxyz,&
                  denspot,GPU,updatePhi,&
                  infoBasisFunctions,infoCoeff,itScc,ebs,coeff,lphi,nlpspd,proj,communicate_lphi,&
                  coeff_proj,ldiis,nit,lin%nItInnerLoop,lin%newgradient,orthpar,confdatarr,&
                  lin%methTransformOverlap,lin%blocksize_pdgemm,lin%convCrit,lin%nItPrecond,&
                  lin%useDerivativeBasisFunctions,lin%lphiRestart,lin%lb%comrp,lin%blocksize_pdsyev,lin%nproc_pdsyev,&
                  hx,hy,hz,input%SIC, input%lin%factor_enlarge)
          end if


          ! Potential from electronic charge density
          !!call mpi_barrier(mpi_comm_world, ierr)
          !!call cpu_time(t1)
          if(lin%mixedmode) then
              if(.not.withder) then
                  call sumrhoForLocalizedBasis2(iproc, nproc, orbs%norb, &
                       lin%lzd, input, hx, hy, hz, lin%orbs, lin%comsr, &
                       coeff, Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3d, &
                       denspot%rhov, at, denspot%dpcom%nscatterarr)
               else
                  call sumrhoForLocalizedBasis2(iproc, nproc, orbs%norb,&
                       lin%lzd, input, hx, hy, hz, lin%lb%orbs, lin%lb%comsr, &
                       coeff, Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3d,&
                       denspot%rhov, at, denspot%dpcom%nscatterarr)
               end if
          else
              call sumrhoForLocalizedBasis2(iproc, nproc, orbs%norb,&
                   lin%lzd, input, hx, hy ,hz, lin%lb%orbs, lin%lb%comsr, &
                   coeff, Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3d, &
                   denspot%rhov, at, denspot%dpcom%nscatterarr)
          end if
          !!call mpi_barrier(mpi_comm_world, ierr)
          !!call cpu_time(t2)
          !!time=t2-t1
          !!if(iproc==0) write(*,'(1x,a,es10.3)') 'time for sumrho:', time

          ! Mix the density.
          if(trim(lin%mixingMethod)=='dens') then
              if(updatePhi) then
                  if(lowaccur_converged) then
                      alphaMix=lin%alphaMixWhenOptimizing_highaccuracy
                  else
                      alphaMix=lin%alphaMixWhenOptimizing_lowaccuracy
                  end if
              else
                  if(lowaccur_converged) then
                      alphaMix=lin%alphaMixWhenFixed_highaccuracy
                  else
                      alphaMix=lin%alphaMixWhenFixed_lowaccuracy
                  end if
              end if
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
!!$          call updatePotential(iproc, nproc, denspot%dpcom%n3d, denspot%dpcom%n3p, Glr, orbs, at, input, lin, &
!!$              denspot%rhov, nscatterarr, pkernel, pot_ion, rhocore, potxc, PSquiet, &
!!$              coeff, ehart, eexcu, vexcu)

          ! Calculate the total energy.
          energy=ebs-ehart+eexcu-vexcu-eexctX+eion+edisp
          energyDiff=energy-energyold
          energyold=energy


          ! Mix the potential
          if(trim(lin%mixingMethod)=='pot') then
              if(updatePhi) then
                  if(lowaccur_converged) then
                      alphaMix=lin%alphaMixWhenOptimizing_highaccuracy
                  else
                      alphaMix=lin%alphaMixWhenOptimizing_lowaccuracy
                  end if
              else
                  if(lowaccur_converged) then
                      alphaMix=lin%alphaMixWhenFixed_highaccuracy
                  else
                      alphaMix=lin%alphaMixWhenFixed_lowaccuracy
                  end if
              end if
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
          !ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
          call allocateCommunicationsBuffersPotential(lin%comgp, subname)
          call postCommunicationsPotential(iproc, nproc, denspot%dpcom%ndimpot, denspot%rhov, lin%comgp)
          if(lin%useDerivativeBasisFunctions) then
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
  if(lin%useDerivativeBasisFunctions) then
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
      call communicate_basis_for_density(iproc, nproc, lin%lzd, lin%lb%orbs, lphi, lin%lb%comsr)
  !!!! Transform all orbitals to real space.
  !!!ist=1
  !!!istr=1
  !!!do iorb=1,lin%lb%orbs%norbp
  !!!    ilr=lin%lb%orbs%inWhichLocreg(lin%lb%orbs%isorb+iorb)
  !!!    call initialize_work_arrays_sumrho(lin%lzd%Llr(ilr), w)
  !!!    call daub_to_isf(lin%lzd%Llr(ilr), w, lphi(ist), lin%lb%comsr%sendBuf(istr))
  !!!    call deallocate_work_arrays_sumrho(w)
  !!!    ist = ist + lin%lzd%Llr(ilr)%wfd%nvctr_c + 7*lin%lzd%Llr(ilr)%wfd%nvctr_f
  !!!    istr = istr + lin%lzd%Llr(ilr)%d%n1i*lin%lzd%Llr(ilr)%d%n2i*lin%lzd%Llr(ilr)%d%n3i
  !!!end do
  !!!if(istr/=lin%lb%comsr%nsendBuf+1) then
  !!!    write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=lin%comsr%nsendBuf+1'
  !!!    stop
  !!!end if

  !!!! Post the MPI messages for the communication of sumrho. Since we use non blocking point
  !!!! to point communication, the program will continue immediately. The messages will be gathered
  !!!! in the subroutine sumrhoForLocalizedBasis2.
  !!!call postCommunicationSumrho2(iproc, nproc, lin%lb%comsr, lin%lb%comsr%sendBuf, lin%lb%comsr%recvBuf)
  call sumrhoForLocalizedBasis2(iproc, nproc, orbs%norb, lin%lzd, input, hx, hy, hz, lin%lb%orbs, lin%lb%comsr, &
       coeff, Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3d, denspot%rhov, at,denspot%dpcom%nscatterarr)

  call deallocateCommunicationbufferSumrho(lin%lb%comsr, subname)

  call mpi_barrier(mpi_comm_world, ierr)
  t1force=mpi_wtime()
  ! Build global orbitals psi (the physical ones).
  if(lin%transformToGlobal) then
      call transformToGlobal(iproc, nproc, lin, orbs, comms, input, coeff, lphi, psi, psit)
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
  call deallocateLinear(iproc, lin, lphi, coeff)
  deallocate(confdatarr)
  call deallocateBasicArrays(at,lin)

  iall=-product(shape(coeff_proj))*kind(coeff_proj)
  deallocate(coeff_proj, stat=istat)
  call memocc(istat, iall, 'coeff_proj', subname)

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
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: psi, psit

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


  call dcopy(orbs%npsidim_comp, psi, 1, psit, 1)

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
