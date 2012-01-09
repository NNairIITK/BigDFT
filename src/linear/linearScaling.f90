subroutine linearScaling(iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh, Glr, orbs, comms, at, input, rhodsc, lin, rxyz, fion, fdisp,&
    radii_cf, nscatterarr, ngatherarr, nlpspd, proj, rhopot, GPU, pkernelseq, irrzon, phnons, pkernel, pot_ion, rhocore, potxc,&
    PSquiet, eion, edisp, eexctX, scpot, psi, psit, energy, fxyz)
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
integer,intent(in):: iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh
type(locreg_descriptors),intent(in) :: Glr
type(orbitals_data),intent(inout):: orbs
type(communications_arrays),intent(in) :: comms
type(atoms_data),intent(inout):: at
type(linearParameters),intent(in out):: lin
type(input_variables),intent(in):: input
type(rho_descriptors),intent(in) :: rhodsc
real(8),dimension(3,at%nat),intent(inout):: rxyz
real(8),dimension(3,at%nat),intent(in):: fion, fdisp
real(8),dimension(at%ntypes,3),intent(in):: radii_cf
integer,dimension(0:nproc-1,4),intent(inout):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!integer,dimension(0:nproc-1,2),intent(in):: ngatherarr
integer,dimension(0:nproc-1,2),intent(inout):: ngatherarr
type(nonlocal_psp_descriptors),intent(in):: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin),intent(in out):: rhopot
type(GPU_pointers),intent(in out):: GPU
real(dp),dimension(:),pointer,intent(in):: pkernelseq
integer, dimension(lin%as%size_irrzon(1),lin%as%size_irrzon(2),lin%as%size_irrzon(3)),intent(in) :: irrzon 
real(dp), dimension(lin%as%size_phnons(1),lin%as%size_phnons(2),lin%as%size_phnons(3)),intent(in) :: phnons 
real(dp), dimension(lin%as%size_pkernel),intent(in):: pkernel
real(wp), dimension(lin%as%size_pot_ion),intent(inout):: pot_ion
!real(wp), dimension(lin%as%size_rhocore):: rhocore 
real(wp), dimension(:),pointer,intent(in):: rhocore                  
real(wp), dimension(lin%as%size_potxc(1),lin%as%size_potxc(2),lin%as%size_potxc(3),lin%as%size_potxc(4)),intent(inout):: potxc
character(len=3),intent(in):: PSquiet
real(gp),intent(in):: eion, edisp, eexctX
logical,intent(in):: scpot
!real(8),dimension(orbs),intent(out):: psi
real(8),dimension(:),pointer,intent(out):: psi, psit
real(8),intent(out):: energy
real(8),dimension(3,at%nat),intent(out):: fxyz
!real(8),intent(out):: fnoise
real(8):: fnoise



! Local variables
integer:: infoBasisFunctions, infoCoeff, istat, iall, itSCC, nitSCC, i, ierr, potshortcut, ndimpot, ist, istr, ilr, tag, itout
real(8),dimension(:),pointer:: phi, phid
real(8),dimension(:,:),pointer:: coeff, coeffd
real(8):: ebs, ebsMod, pnrm, tt, ehart, eexcu, vexcu, alphaMix, dampingForMixing
character(len=*),parameter:: subname='linearScaling'
real(8),dimension(:),allocatable:: rhopotOld, rhopotold_out
type(linearParameters):: lind
logical:: updatePhi, reduceConvergenceTolerance, communicate_lphi, with_auxarray
real(8),dimension(:),pointer:: lphi, lphir, phibuffr

integer,dimension(:,:),allocatable:: nscatterarrTemp !n3d,n3p,i3s+i3xcsh-1,i3xcsh
real(8),dimension(:),allocatable:: phiTemp, lphiold
real(wp),dimension(:),allocatable:: projTemp
real(8):: t1, t2, time, t1tot, t2tot, timetot, t1ig, t2ig, timeig, t1init, t2init, timeinit, ddot, dnrm2, pnrm_out
real(8):: t1scc, t2scc, timescc, t1force, t2force, timeforce, energyold, energyDiff, energyoldout, selfConsistent
character(len=11):: orbName
character(len=10):: comment, procName, orbNumber
integer:: iorb, istart, sizeLphir, sizePhibuffr, ndimtot, iiorb, ncount, ncnt
type(mixrhopotDIISParameters):: mixdiis
type(workarr_sumrho):: w
real(8),dimension(:,:),allocatable:: ovrlp, coeff_proj








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
  call allocateAndInitializeLinear(iproc, nproc, Glr, orbs, at, nlpspd, lin, phi, &
       input, rxyz, nscatterarr, tag, coeff, lphi)

!!! Determine lin%cutoffweight
!!ist=1
!!lphi=1.d0
!!do iorb=1,lin%orbs%norbp
!!    iiorb=lin%orbs%isorb+iorb
!!    ilr=lin%orbs%inwhichlocreg(iiorb)
!!    ncnt = lin%lzd%llr(ilr)%wfd%nvctr_c + 7*lin%lzd%llr(ilr)%wfd%nvctr_f
!!    tt=sqrt(dble(ncnt))
!!    call dscal(ncnt, 1/tt, lphi(ist), 1)
!!    ist = ist + ncnt
!!end do
!!allocate(lin%lzd%cutoffweight(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, lin%lzd%cutoffweight, 'lin%lzd%cutoffweight', subname)
!!call allocateCommuncationBuffersOrtho(lin%comon, subname)
!!call getMatrixElements2(iproc, nproc, lin%lzd, lin%orbs, lin%op, lin%comon, lphi, lphi, lin%mad, lin%lzd%cutoffweight)
!!call deallocateCommuncationBuffersOrtho(lin%comon, subname)
!!!!call getOverlapMatrix2(iproc, nproc, lin%lzd, lin%orbs, lin%comon, lin%op, lphi, lin%mad, lin%lzd%cutoffweight)
!!do iorb=1,lin%orbs%norb
!!    do iiorb=1,lin%orbs%norb
!!        write(90+iproc,*) iorb, iiorb, lin%lzd%cutoffweight(iiorb,iorb)
!!    end do
!!end do

  call mpi_barrier(mpi_comm_world, ierr)
  t2init=mpi_wtime()
  timeinit=t2init-t1init

  if(.not.lin%transformToGlobal) then
      ! psi and psit will not be calculated, so only allocate them with size 1
      orbs%npsidim=1
  end if
  allocate(psi(orbs%npsidim), stat=istat)
  call memocc(istat, psi, 'psi', subname)
  allocate(psit(orbs%npsidim), stat=istat)
  call memocc(istat, psit, 'psit', subname)
  allocate(rhopotold(max(glr%d%n1i*glr%d%n2i*n3p,1)*input%nspin), stat=istat)
  call memocc(istat, rhopotold, 'rhopotold', subname)
  allocate(rhopotold_out(max(glr%d%n1i*glr%d%n2i*n3p,1)*input%nspin), stat=istat)
  call memocc(istat, rhopotold_out, 'rhopotold_out', subname)
  !rhopotold_out=1.d100

  allocate(lin%coeffall(lin%lb%orbs%norb,orbs%norb+lin%norbvirt), stat=istat)
  call memocc(istat, lin%coeffall, 'lin%coeffall', subname)

  allocate(coeff_proj(lin%orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, coeff_proj, 'coeff_proj', subname)


  call prepare_lnlpspd(iproc, at, input, lin%orbs, rxyz, radii_cf, lin%locregShape, lin%lzd)

  potshortcut=0 ! What is this?
  call mpi_barrier(mpi_comm_world, ierr)
  t1ig=mpi_wtime()
  call inputguessConfinement(iproc, nproc, at, &
       comms, Glr, input, rhodsc, lin, orbs, rxyz, n3p, rhopot, rhopotold, rhocore, pot_ion,&
       nlpspd, proj, pkernel, pkernelseq, &
       nscatterarr, ngatherarr, potshortcut, irrzon, phnons, GPU, radii_cf, &
       tag, lphi, ehart, eexcu, vexcu)
  call mpi_barrier(mpi_comm_world, ierr)
  t2ig=mpi_wtime()
  timeig=t2ig-t1ig
  t1scc=mpi_wtime()

  !!! Copy lphi to lin%lpsi, don't know whether this is a good choice
  !!lin%lpsi=lphi


  ! Initialize the DIIS mixing of the potential if required.
  if(lin%mixHist>0) then
      ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
      call initializeMixrhopotDIIS(lin%mixHist, ndimpot, mixdiis)
  end if

  !if(lin%nItInguess>0 .and. .false.) then
  if(lin%nItInguess>0 .and. .true.) then
      ! Post communications for gathering the potential.
      ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
      call allocateCommunicationsBuffersPotential(lin%comgp, subname)
      call postCommunicationsPotential(iproc, nproc, ndimpot, rhopot, lin%comgp)
      if(lin%useDerivativeBasisFunctions) then
          call allocateCommunicationsBuffersPotential(lin%lb%comgp, subname)
          call postCommunicationsPotential(iproc, nproc, ndimpot, rhopot, lin%lb%comgp)
      end if

      ! Calculate the Hamiltonian in the basis of the trace minimizing orbitals. Do not improve
      ! the basis functions (therefore update is set to false).
      ! This subroutine will also post the point to point messages needed for the calculation
      ! of the charge density.
      updatePhi=.false.
      !communicate_lphi=.true.
      communicate_lphi=.true.
      with_auxarray=.false.
      call allocateCommunicationbufferSumrho(iproc, with_auxarray, lin%comsr, subname)
      call getLinearPsi(iproc, nproc, input%nspin, Glr, orbs, comms, at, lin, rxyz, rxyz, &
          nscatterarr, ngatherarr, rhopot, GPU, input, pkernelseq, phi, psi, psit, updatePhi, &
          infoBasisFunctions, infoCoeff, 0, n3p, n3pi, n3d, pkernel, &
          i3s, i3xcsh, ebs, coeff, lphi, radii_cf, nlpspd, proj, communicate_lphi, coeff_proj)

      ! Calculate the charge density.
      call cpu_time(t1)
      !call sumrhoForLocalizedBasis2(iproc, nproc, orbs, Glr, input, lin, coeff, lphi, Glr%d%n1i*Glr%d%n2i*n3d, &
      !     rhopot, at, nscatterarr)
      !do istat=1,Glr%d%n1i*Glr%d%n2i*n3d
      !    write(1200+iproc,*) istat, rhopot(istat)
      !end do

      !call sumrholinear_auxiliary(iproc, nproc, orbs, Glr, input, lin, coeff, phi, at, nscatterarr)
      !call sumrholinear_withauxiliary(iproc, nproc, orbs, Glr, input, lin, coeff, Glr%d%n1i*Glr%d%n2i*n3d, &
      !     rhopot, at, nscatterarr)
      !do istat=1,Glr%d%n1i*Glr%d%n2i*n3d
      !    write(1100+iproc,*) istat, rhopot(istat)
      !end do
      call deallocateCommunicationbufferSumrho(lin%comsr, subname)
      call cpu_time(t2)
      time=t2-t1
      call mpiallred(time, 1, mpi_sum, mpi_comm_world, ierr)
      time=time/dble(nproc)
      if(iproc==0) write(*,'(1x,a,es12.4)') 'time for sumrho:',time

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

      ! if we mix the density, copy the current charge density.
      !allocate(rhopotold(max(glr%d%n1i*glr%d%n2i*n3p,1)*input%nspin), stat=istat)
      !call memocc(istat, rhopotold, 'rhopotold', subname)
      !if(trim(lin%mixingmethod)=='dens') then
      !    call dcopy(max(glr%d%n1i*glr%d%n2i*n3p,1)*input%nspin, rhopot(1), 1, rhopotold(1), 1)
      !end if

      !! Calculate the potential we get with the current chareg density.
      !call updatePotential(iproc, nproc, n3d, n3p, Glr, orbs, at, input, lin, phi,  &
      !    rhopot, nscatterarr, pkernel, pot_ion, rhocore, potxc, PSquiet, &
      !    coeff, ehart, eexcu, vexcu)

      if(trim(lin%mixingMethod)=='pot') then
          if(lin%mixHist==0) then
              call mixPotential(iproc, n3p, Glr, input, lin%alphaMixWhenFixed, rhopotOld, rhopot, pnrm)
          else 
              ndimpot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
              ndimtot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*lin%lzd%Glr%d%n3i
              mixdiis%mis=mod(mixdiis%is,mixdiis%isx)+1
              mixdiis%is=mixdiis%is+1
              call mixrhopotDIIS(iproc, nproc, ndimpot, rhopot, rhopotold, mixdiis, ndimtot, lin%alphaMixWhenFixed, 2, pnrm)
          end if
          rhopotold_out=rhopot
      end if

      ! Copy the current potential
      if(trim(lin%mixingMethod)=='pot') then
           call dcopy(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin, rhopot(1), 1, rhopotOld(1), 1)
      end if
  end if

  !!! If we mix the potential, copy the potential.
  !!if(trim(lin%mixingMethod)=='pot') then
  !!    call dcopy(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin, rhopot(1), 1, rhopotOld(1), 1)
  !!end if

  ! Allocate the communications buffers needed for the communications of teh potential and
  ! post the messages. This will send to each process the part of the potential that this process
  ! needs for the application of the Hamlitonian to all orbitals on that process.
  ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
  call allocateCommunicationsBuffersPotential(lin%comgp, subname)
  call postCommunicationsPotential(iproc, nproc, ndimpot, rhopot, lin%comgp)

  ! If we also use the derivative of the basis functions, also send the potential in this case. This is
  ! needed since the orbitals may be partitioned in a different way when the derivatives are used.
  if(lin%useDerivativeBasisFunctions) then
      call allocateCommunicationsBuffersPotential(lin%lb%comgp, subname)
      call postCommunicationsPotential(iproc, nproc, ndimpot, rhopot, lin%lb%comgp)
  end if


  !! For mixing phi together with mixing the density -- experimental
  !lphiold=lphi


  if(nproc==1) allocate(psit(size(psi)))
  nitSCC=lin%nitSCCWhenOptimizing+lin%nitSCCWhenFixed
  ! Flag that indicates that the basis functions shall be improved in the following.
  updatePhi=.true.
  pnrm=1.d100
  energyold=0.d0
  energyoldout=0.d0
  !lin%getCoeff='new'
  reduceConvergenceTolerance=.false.
  !if(iproc==0) write(*,'(a,es9.2)') 'dampingForMixing',dampingForMixing


  lin%newgradient=.false.

  do itout=1,lin%nItOuterSCC

      updatePhi=.true.

      if(reduceConvergenceTolerance) lin%fixBasis=max(lin%fixBasis*lin%factorFixBasis,lin%minimalFixBasis)
      selfConsistent=max(lin%convCritMix,5.d-3*lin%fixBasis)

      if(iproc==0) write(*,'(a,es12.4,3x,es12.4)') &
           'DELTA DENS for fixing basis functions, reaching self consistency:',lin%fixBasis, selfConsistent

      if(lin%sumrho_fast) then
          with_auxarray=.true.
      else
          with_auxarray=.false.
      end if
      call allocateCommunicationbufferSumrho(iproc, with_auxarray, lin%comsr, subname)

      if(itout==lin%nit_lowaccuracy) then
           lin%potentialPrefac = 5.d-3*lin%potentialPrefac
           lin%nItBasisFirst = 5*lin%nItBasisFirst
           lin%nItBasis = 30*lin%nItBasis
           lin%newgradient=.true.
      end if

      do itSCC=1,nitSCC
          if(itSCC>1 .and. pnrm<lin%fixBasis .or. itSCC==lin%nitSCCWhenOptimizing) updatePhi=.false.
          if(itSCC==1) then
              communicate_lphi=.true.
          else
              communicate_lphi=.false.
          end if
          ! This subroutine gives back the new psi and psit, which are a linear combination of localized basis functions.
          call getLinearPsi(iproc, nproc, input%nspin, Glr, orbs, comms, at, lin, rxyz, rxyz, &
              nscatterarr, ngatherarr, rhopot, GPU, input, pkernelseq, phi, psi, psit, updatePhi, &
              infoBasisFunctions, infoCoeff, itScc, n3p, n3pi, n3d, pkernel, &
              i3s, i3xcsh, ebs, coeff, lphi, radii_cf, nlpspd, proj, communicate_lphi, coeff_proj)


          ! Potential from electronic charge density
          call mpi_barrier(mpi_comm_world, ierr)
          call cpu_time(t1)
          if(.not. lin%sumrho_fast) then
              call sumrhoForLocalizedBasis2(iproc, nproc, orbs, Glr, input, lin, coeff, phi, Glr%d%n1i*Glr%d%n2i*n3d, &
                   rhopot, at, nscatterarr)
          else
              if(itSCC==1) then
                  call sumrholinear_auxiliary(iproc, nproc, orbs, Glr, input, lin, coeff, phi, at, nscatterarr)
                  call sumrholinear_withauxiliary(iproc, nproc, orbs, Glr, input, lin, coeff, Glr%d%n1i*Glr%d%n2i*n3d, &
                       rhopot, at, nscatterarr)
              else
                  call sumrholinear_withauxiliary(iproc, nproc, orbs, Glr, input, lin, coeff, Glr%d%n1i*Glr%d%n2i*n3d, &
                       rhopot, at, nscatterarr)
              end if
          end if
          call mpi_barrier(mpi_comm_world, ierr)
          call cpu_time(t2)
          time=t2-t1
          !!call mpiallred(time, 1, mpi_sum, mpi_comm_world, ierr)
          if(iproc==0) write(*,'(1x,a,es10.3)') 'time for sumrho:', time

          ! Mix the density.
          if(trim(lin%mixingMethod)=='dens') then
              if(updatePhi) then
                  alphaMix=lin%alphaMixWhenOptimizing
              else
                  alphaMix=lin%alphaMixWhenFixed
              end if
              if(lin%mixHist==0) then
                  call mixPotential(iproc, n3p, Glr, input, alphaMix, rhopotOld, rhopot, pnrm)
              else 
                  ndimpot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
                  ndimtot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*lin%lzd%Glr%d%n3i
                  mixdiis%mis=mod(mixdiis%is,mixdiis%isx)+1
                  mixdiis%is=mixdiis%is+1
                  call mixrhopotDIIS(iproc, nproc, ndimpot, rhopot, rhopotold, mixdiis, ndimtot, alphaMix, 1, pnrm)
              end if
              if(pnrm<selfConsistent .or. itSCC==nitSCC) then
                  pnrm_out=0.d0
                  do i=1,Glr%d%n1i*Glr%d%n2i*n3p
                      pnrm_out=pnrm_out+(rhopot(i)-rhopotOld_out(i))**2
                  end do
                  call mpiallred(pnrm_out, 1, mpi_sum, mpi_comm_world, ierr)
                  pnrm_out=sqrt(pnrm_out)/(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*input%nspin)
                  call dcopy(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin, rhopot(1), 1, rhopotOld_out(1), 1)
              end if
          end if

          ! Copy the current charge density.
          if(trim(lin%mixingMethod)=='dens') then
              call dcopy(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin, rhopot(1), 1, rhopotOld(1), 1)
          end if

          ! Calculate the new potential.
          call updatePotential(iproc, nproc, n3d, n3p, Glr, orbs, at, input, lin, phi,  &
              rhopot, nscatterarr, pkernel, pot_ion, rhocore, potxc, PSquiet, &
              coeff, ehart, eexcu, vexcu)
          ! Calculate the total energy.
          energy=ebs-ehart+eexcu-vexcu-eexctX+eion+edisp
          energyDiff=energy-energyold
          energyold=energy


          ! Post communications for gathering the potential
          ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)

          ! Mix the potential
          if(trim(lin%mixingMethod)=='pot') then
              if(updatePhi) then
                  alphaMix=lin%alphaMixWhenOptimizing
              else
                  alphaMix=lin%alphaMixWhenFixed
              end if
              if(lin%mixHist==0) then
                  call mixPotential(iproc, n3p, Glr, input, alphaMix, rhopotOld, rhopot, pnrm)
              else 
                  ndimpot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
                  ndimtot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*lin%lzd%Glr%d%n3i
                  mixdiis%mis=mod(mixdiis%is,mixdiis%isx)+1
                  mixdiis%is=mixdiis%is+1
                  call mixrhopotDIIS(iproc, nproc, ndimpot, rhopot, rhopotold, mixdiis, ndimtot, alphaMix, 2, pnrm)
              end if
              if(pnrm<selfConsistent .or. itSCC==nitSCC) then
                  pnrm_out=0.d0
                  do i=1,Glr%d%n1i*Glr%d%n2i*n3p
                      pnrm_out=pnrm_out+(rhopot(i)-rhopotOld_out(i))**2
                  end do
                  call mpiallred(pnrm_out, 1, mpi_sum, mpi_comm_world, ierr)
                  pnrm_out=sqrt(pnrm_out)/(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*input%nspin)
                  call dcopy(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin, rhopot(1), 1, rhopotOld_out(1), 1)
              end if
          end if

          ! Copy the current potential
          if(trim(lin%mixingMethod)=='pot') then
               call dcopy(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin, rhopot(1), 1, rhopotOld(1), 1)
          end if

          ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
          call allocateCommunicationsBuffersPotential(lin%comgp, subname)
          call postCommunicationsPotential(iproc, nproc, ndimpot, rhopot, lin%comgp)
          if(lin%useDerivativeBasisFunctions) then
              call allocateCommunicationsBuffersPotential(lin%lb%comgp, subname)
              call postCommunicationsPotential(iproc, nproc, ndimpot, rhopot, lin%lb%comgp)
          end if

          ! Write some informations
          call printSummary(iproc, itSCC, infoBasisFunctions, infoCoeff, pnrm, energy, energyDiff, lin%mixingMethod)
          if(pnrm<selfConsistent) then
              reduceConvergenceTolerance=.true.
              exit
          else
              reduceConvergenceTolerance=.false.
          end if
      end do

      call deallocateCommunicationbufferSumrho(lin%comsr, subname)

      if(iproc==0) then
          if(trim(lin%mixingMethod)=='dens') then
              write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                   'itout, Delta DENSOUT, energy, energyDiff', itout, pnrm_out, energy, energy-energyoldout
          else if(trim(lin%mixingMethod)=='pot') then
              write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                   'itout, Delta POTOUT, energy energyDiff', itout, pnrm_out, energy, energy-energyoldout
          end if
      end if
      if(abs(pnrm_out)<lin%convCritMixOut) exit
      energyoldout=energy
  end do

  call cancelCommunicationPotential(iproc, nproc, lin%comgp)
  call deallocateCommunicationsBuffersPotential(lin%comgp, subname)
  if(lin%useDerivativeBasisFunctions) then
      call cancelCommunicationPotential(iproc, nproc, lin%lb%comgp)
      call deallocateCommunicationsBuffersPotential(lin%lb%comgp, subname)
  end if

  iall=-product(shape(rhopotOld))*kind(rhopotOld)
  deallocate(rhopotOld, stat=istat)
  call memocc(istat, iall, 'rhopotOld', subname)
  iall=-product(shape(rhopotold_out))*kind(rhopotold_out)
  deallocate(rhopotold_out, stat=istat)
  call memocc(istat, iall, 'rhopotold_out', subname)

  if(lin%mixHist>0) then
      call deallocateMixrhopotDIIS(mixdiis)
  end if

  call mpi_barrier(mpi_comm_world, ierr)
  t2scc=mpi_wtime()
  timescc=t2scc-t1scc


  ! Allocate the communication buffers for the calculation of the charge density.
  with_auxarray=.false.
  call allocateCommunicationbufferSumrho(iproc, with_auxarray, lin%comsr, subname)
  ! Transform all orbitals to real space.
  ist=1
  istr=1
  do iorb=1,lin%lb%orbs%norbp
      ilr=lin%lb%orbs%inWhichLocregp(iorb)
      call initialize_work_arrays_sumrho(lin%lzd%Llr(ilr), w)
      !call daub_to_isf(lin%lzd%llr(ilr), w, lphi(ist), lphir(istr))
      call daub_to_isf(lin%lzd%Llr(ilr), w, lphi(ist), lin%comsr%sendBuf(istr))
      call deallocate_work_arrays_sumrho(w)
      ist = ist + lin%lzd%Llr(ilr)%wfd%nvctr_c + 7*lin%lzd%Llr(ilr)%wfd%nvctr_f
      istr = istr + lin%lzd%Llr(ilr)%d%n1i*lin%lzd%Llr(ilr)%d%n2i*lin%lzd%Llr(ilr)%d%n3i
  end do
  if(istr/=lin%comsr%nsendBuf+1) then
      write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=lin%comsr%nsendBuf+1'
      stop
  end if

  ! Post the MPI messages for the communication of sumrho. Since we use non blocking point
  ! to point communication, the program will continue immediately. The messages will be gathered
  ! in the subroutine sumrhoForLocalizedBasis2.
  call postCommunicationSumrho2(iproc, nproc, lin, lin%comsr%sendBuf, lin%comsr%recvBuf)
  call sumrhoForLocalizedBasis2(iproc, nproc, orbs, Glr, input, lin, coeff, phi, Glr%d%n1i*Glr%d%n2i*n3d, &
       rhopot, at, nscatterarr)

  !call sumrholinear_auxiliary(iproc, nproc, orbs, Glr, input, lin, coeff, phi, at, nscatterarr)
  !call sumrholinear_withauxiliary(iproc, nproc, orbs, Glr, input, lin, coeff, Glr%d%n1i*Glr%d%n2i*n3d, &
  !     rhopot, at, nscatterarr)
  call deallocateCommunicationbufferSumrho(lin%comsr, subname)

  call mpi_barrier(mpi_comm_world, ierr)
  t1force=mpi_wtime()
  ! Build global orbitals psi (the physical ones).
  if(lin%transformToGlobal) then
      call transformToGlobal(iproc, nproc, lin, orbs, comms, input, coeff, lphi, psi, psit)
  end if


  ! Calculate the forces we get with psi.
  !!call calculateForcesSub(iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh, Glr, orbs, at, input, comms, lin, nlpspd, &
  !!    proj, ngatherarr, nscatterarr, GPU, irrzon, phnons, pkernel, rxyz, fion, fdisp, lphi, coeff, rhopot, &
  !!    fxyz, fnoise,radii_cf)

  call calculateForcesLinear(iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh, Glr, orbs, at, input, comms, lin, nlpspd, &
       proj, ngatherarr, nscatterarr, GPU, irrzon, phnons, pkernel, rxyz, fion, fdisp, rhopot, psi, fxyz, fnoise)
  call mpi_barrier(mpi_comm_world, ierr)
  t2force=mpi_wtime()
  timeforce=t2force-t1force


  call free_lnlpspd(lin%orbs, lin%lzd)

  ! Deallocate all arrays related to the linear scaling version.
  call deallocateLinear(iproc, lin, phi, lphi, coeff)

  iall=-product(shape(lin%coeffall))*kind(lin%coeffall)
  deallocate(lin%coeffall, stat=istat)
  call memocc(istat, iall, 'lin%coeffall', subname)

  iall=-product(shape(coeff_proj))*kind(coeff_proj)
  deallocate(coeff_proj, stat=istat)
  call memocc(istat, iall, 'coeff_proj', subname)



  call mpi_barrier(mpi_comm_world, ierr)
  t2tot=mpi_wtime()
  timetot=t2tot-t1tot
  if(iproc==0) write(*,'(1x,a)') '================================================'
  if(iproc==0) write(*,'(1x,a,es10.3,a)') 'total time for linear scaling version:',timetot,'s'
  if(iproc==0) write(*,'(3x,a)') 'of which:'
  if(iproc==0) write(*,'(13x,a,es10.3,a,f4.1,a)') '- initialization:',timeinit,'s (',timeinit/timetot*100.d0,'%)'
  if(iproc==0) write(*,'(13x,a,es10.3,a,f4.1,a)') '- input guess:',timeig,'s (',timeig/timetot*100.d0,'%)'
  if(iproc==0) write(*,'(13x,a,es10.3,a,f4.1,a)') '- self consistency cycle:',timescc,'s (',timescc/timetot*100.d0,'%)'
  if(iproc==0) write(*,'(13x,a,es10.3,a,f4.1,a)') '- forces:',timeforce,'s (',timeforce/timetot*100.d0,'%)'
  if(iproc==0) write(*,'(1x,a)') '================================================'

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
type(p2pCommsGatherPot),intent(inout):: comgp

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




subroutine transformToGlobal(iproc, nproc, lin, orbs, comms, input, coeff, lphi, psi, psit)
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
real(8),dimension(lin%lb%orbs%npsidim),intent(inout):: lphi
real(8),dimension(orbs%npsidim),intent(out):: psi, psit

! Local variables
integer:: ind1, ind2, istat, iall, iorb, ilr, ldim, gdim, nvctrp
real(8),dimension(:),pointer:: phiWork
real(8),dimension(:),allocatable:: phi
character(len=*),parameter:: subname='transformToGlobal'
  !do iall=0,nproc-1
  !    write(*,'(a,i5,4i12)') 'START transformToGlobal: iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)', iproc, comms%ncntt(iall), comms%ndsplt(iall), comms%ncntd(iall), comms%ndspld(iall)  
  !end do

  allocate(phi(lin%lb%gorbs%npsidim), stat=istat)
  call memocc(istat, phi, 'phi', subname)
  allocate(phiWork(max(size(phi),size(psi))), stat=istat)
  call memocc(istat, phiWork, 'phiWork', subname)

  ind1=1
  ind2=1
  phi=0.d0
  do iorb=1,lin%lb%orbs%norbp
      ilr = lin%lb%orbs%inWhichLocregp(iorb)
      ldim=lin%lzd%Llr(ilr)%wfd%nvctr_c+7*lin%lzd%Llr(ilr)%wfd%nvctr_f
      gdim=lin%lzd%Glr%wfd%nvctr_c+7*lin%lzd%Glr%wfd%nvctr_f
      call Lpsi_to_global2(iproc, nproc, ldim, gdim, lin%lb%orbs%norb, lin%lb%orbs%nspinor, input%nspin, lin%lzd%Glr,&
           lin%lzd%Llr(ilr), lphi(ind2), phi(ind1))
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


  call dcopy(orbs%npsidim, psi, 1, psit, 1)

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
