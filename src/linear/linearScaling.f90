subroutine linearScaling(iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh, Glr, orbs, comms, at, input, lin, rxyz, fion, fdisp, radii_cf, &
    nscatterarr, ngatherarr, nlpspd, proj, rhopot, GPU, pkernelseq, irrzon, phnons, pkernel, pot_ion, rhocore, potxc, PSquiet, &
    eion, edisp, eexctX, scpot, psi, psit, energy, fxyz)
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
type(orbitals_data),intent(in):: orbs
type(communications_arrays),intent(in) :: comms
type(atoms_data),intent(inout):: at
type(linearParameters),intent(in out):: lin
type(input_variables),intent(in):: input
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
real(8),dimension(orbs%npsidim),intent(out):: psi
real(8),dimension(:),pointer,intent(out):: psit
real(8),intent(out):: energy
real(8),dimension(3,at%nat),intent(out):: fxyz
!real(8),intent(out):: fnoise
real(8):: fnoise



! Local variables
integer:: infoBasisFunctions, infoCoeff, istat, iall, itSCC, nitSCC, i, ierr, potshortcut, ndimpot, ist, istr, ilr
real(8),dimension(:),pointer:: phi, phid
real(8),dimension(:,:),pointer:: coeff, coeffd
real(8):: ebs, ebsMod, pnrm, tt, ehart, eexcu, vexcu
character(len=*),parameter:: subname='linearScaling'
real(8),dimension(:),allocatable:: rhopotOld
type(linearParameters):: lind
logical:: updatePhi
real(8),dimension(:),pointer:: lphi, lphir, phibuffr

integer,dimension(:,:),allocatable:: nscatterarrTemp !n3d,n3p,i3s+i3xcsh-1,i3xcsh
real(8),dimension(:),allocatable:: phiTemp
real(wp),dimension(:),allocatable:: projTemp
real(8):: t1, t2, time

character(len=11):: orbName
character(len=10):: comment, procName, orbNumber
integer:: iorb, istart, sizeLphir, sizePhibuffr, ndimtot
type(mixrhopotDIISParameters):: mixdiis
type(workarr_sumrho):: w




  if(iproc==0) then
      write(*,'(x,a)') repeat('*',84)
      write(*,'(x,a)') '****************************** LINEAR SCALING VERSION ******************************'
  end if

  ! Initialize the parameters for the linear scaling version and allocate all arrays.
  call allocateAndInitializeLinear(iproc, nproc, Glr, orbs, at, nlpspd, lin, phi, &
       input, rxyz, nscatterarr, coeff, lphi)

  potshortcut=0 ! What is this?
  call inputguessConfinement(iproc, nproc, at, &
       comms, Glr, input, lin, orbs, rxyz, n3p, rhopot, rhocore, pot_ion,&
       nlpspd, proj, pkernel, pkernelseq, &
       nscatterarr, ngatherarr, potshortcut, irrzon, phnons, GPU, radii_cf, &
       lphi, ehart, eexcu, vexcu)
  !!do iall=1,size(rhopot)
  !!    read(10000+iproc,*) rhopot(iall)
  !!end do
  !!do iall=1,size(lphi)
  !!    !read(11000+iproc,*) lphi(iall)
  !!    write(500+iproc,*) lphi(iall)
  !!end do

  ! Post communications for gathering the potential
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
  call getLinearPsi(iproc, nproc, input%nspin, Glr, orbs, comms, at, lin, rxyz, rxyz, &
      nscatterarr, ngatherarr, rhopot, GPU, input, pkernelseq, phi, psi, psit, updatePhi, &
      infoBasisFunctions, infoCoeff, itScc, n3p, n3pi, n3d, pkernel, &
      i3s, i3xcsh, ebs, coeff, lphi, radii_cf, nlpspd, proj)

  ! Calculate the charge density.
  call cpu_time(t1)
  call sumrhoForLocalizedBasis2(iproc, nproc, orbs, Glr, input, lin, coeff, phi, Glr%d%n1i*Glr%d%n2i*n3d, &
       rhopot, at, nscatterarr)
  call deallocateCommunicationbufferSumrho(lin%comsr, subname)
  call cpu_time(t2)
  time=t2-t1
  call mpiallred(time, 1, mpi_sum, mpi_comm_world, ierr)
  time=time/dble(nproc)
  if(iproc==0) write(*,'(x,a,es12.4)') 'time for sumrho:',time

  ! If we mix the density, copy the current charge density.
  allocate(rhopotOld(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin), stat=istat)
  call memocc(istat, rhopotOld, 'rhopotOld', subname)
  if(trim(lin%mixingMethod)=='dens') then
      call dcopy(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin, rhopot(1), 1, rhopotOld(1), 1)
  end if

  ! Calculate the potential we get with the current chareg density.
  call updatePotential(iproc, nproc, n3d, n3p, Glr, orbs, at, input, lin, phi,  &
      rhopot, nscatterarr, pkernel, pot_ion, rhocore, potxc, PSquiet, &
      coeff, ehart, eexcu, vexcu)

  ! If we mix the potential, copy the potential.
  if(trim(lin%mixingMethod)=='pot') then
      call dcopy(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin, rhopot(1), 1, rhopotOld(1), 1)
  end if

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



  ! Initialize the DIIS mixing of the potential if required.
  if(lin%mixHist>0) then
      call initializeMixrhopotDIIS(lin%mixHist, ndimpot, mixdiis)
  end if

  if(nproc==1) allocate(psit(size(psi)))
  nitSCC=lin%nitSCC
  ! Flag that indicates that the basis functions shall be improved in the following.
  updatePhi=.true.
  do itSCC=1,nitSCC
      !if(itSCC==10) updatePhi=.false.
      if(itSCC>1 .and. pnrm<1.d-8) updatePhi=.false.
      !!if(itSCC>1 .and. pnrm<7.3d-9) lin%nItBasis=1
      ! This subroutine gives back the new psi and psit, which are a linear combination of localized basis functions.
      call getLinearPsi(iproc, nproc, input%nspin, Glr, orbs, comms, at, lin, rxyz, rxyz, &
          nscatterarr, ngatherarr, rhopot, GPU, input, pkernelseq, phi, psi, psit, updatePhi, &
          infoBasisFunctions, infoCoeff, itScc, n3p, n3pi, n3d, pkernel, &
          i3s, i3xcsh, ebs, coeff, lphi, radii_cf, nlpspd, proj)


      ! Copy the current potential
      if(trim(lin%mixingMethod)=='pot') then
           call dcopy(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin, rhopot(1), 1, rhopotOld(1), 1)
      end if


      ! Potential from electronic charge density
      !!call cpu_time(t1)
      call sumrhoForLocalizedBasis2(iproc, nproc, orbs, Glr, input, lin, coeff, phi, Glr%d%n1i*Glr%d%n2i*n3d, &
           rhopot, at, nscatterarr)
      call deallocateCommunicationbufferSumrho(lin%comsr, subname)
      !!call cpu_time(t2)
      !!time=t2-t1
      !!call mpiallred(time, 1, mpi_sum, mpi_comm_world, ierr)
      !!if(iproc==0) write(*,'(x,a,es10.3)') 'time for sumrho:', time/dble(nproc)

      ! Mix the density.
      if(trim(lin%mixingMethod)=='dens') then
          if(lin%mixHist==0) then
              call mixPotential(iproc, n3p, Glr, input, lin, rhopotOld, rhopot, pnrm)
          else 
              ndimpot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
              ndimtot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*lin%lzd%Glr%d%n3i
              mixdiis%mis=mod(mixdiis%is,mixdiis%isx)+1
              mixdiis%is=mixdiis%is+1
              call mixrhopotDIIS(iproc, nproc, ndimpot, rhopot, rhopotold, mixdiis, ndimtot, lin%alphaMix, 1, pnrm)
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


      ! Post communications for gathering the potential
      ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)

      !!!!! TEST  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!  ! Calculate the forces we get with psi.
      !!  allocate(nscatterarrTemp(0:nproc-1,4), stat=istat)
      !!  call memocc(istat, nscatterarrTemp, 'nscatterarrTemp', subname)
      !!  allocate(phiTemp(size(phi)), stat=istat)
      !!  call memocc(istat, phiTemp, 'phiTemp', subname)
      !!  allocate(projTemp(nlpspd%nprojel), stat=istat)
      !!  call memocc(istat, projTemp, 'projTemp', subname)
      !!  projTemp=proj
      !!  nscatterarrTemp=nscatterarr
      !!  phiTemp=phi
      !!  call calculateForcesSub(iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh, Glr, orbs, at, input, comms, lin, nlpspd, &
      !!      proj, ngatherarr, nscatterarr, GPU, irrzon, phnons, pkernel, rxyz, fion, fdisp, psi, phi, coeff, fxyz, fnoise)
      !!  proj=projTemp
      !!  nscatterarr=nscatterarrTemp
      !!  phi=phiTemp
      !!  iall=-product(shape(nscatterarrTemp))*kind(nscatterarrTemp)
      !!  deallocate(nscatterarrTemp, stat=istat)
      !!  call memocc(istat, iall, 'nscatterarrTemp', subname)
      !!  iall=-product(shape(phiTemp))*kind(phiTemp)
      !!  deallocate(phiTemp, stat=istat)
      !!  call memocc(istat, iall, 'phiTemp', subname)
      !!  iall=-product(shape(projTemp))*kind(projTemp)
      !!  deallocate(projTemp, stat=istat)
      !!  call memocc(istat, iall, 'projTemp', subname)
      !!!!!  TEST  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Mix the potential
      if(trim(lin%mixingMethod)=='pot') then
          if(lin%mixHist==0) then
           call mixPotential(iproc, n3p, Glr, input, lin, rhopotOld, rhopot, pnrm)
          else 
              ndimpot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
              ndimtot=lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*lin%lzd%Glr%d%n3i
              mixdiis%mis=mod(mixdiis%is,mixdiis%isx)+1
              mixdiis%is=mixdiis%is+1
              call mixrhopotDIIS(iproc, nproc, ndimpot, rhopot, rhopotold, mixdiis, ndimtot, lin%alphaMix, 2, pnrm)
          end if
      end if
      ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
      call allocateCommunicationsBuffersPotential(lin%comgp, subname)
      call postCommunicationsPotential(iproc, nproc, ndimpot, rhopot, lin%comgp)
      if(lin%useDerivativeBasisFunctions) then
          call allocateCommunicationsBuffersPotential(lin%lb%comgp, subname)
          call postCommunicationsPotential(iproc, nproc, ndimpot, rhopot, lin%lb%comgp)
      end if

      ! Write some informations
      call printSummary(iproc, itSCC, infoBasisFunctions, infoCoeff, pnrm, energy, lin%mixingMethod)
      if(pnrm<lin%convCritMix) exit
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

  if(lin%mixHist>0) then
      call deallocateMixrhopotDIIS(mixdiis)
  end if


  ! Allocate the communication buffers for the calculation of the charge density.
  call allocateCommunicationbufferSumrho(lin%comsr, subname)
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
  call deallocateCommunicationbufferSumrho(lin%comsr, subname)


  ! Calculate the forces we get with psi.
  call calculateForcesSub(iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh, Glr, orbs, at, input, comms, lin, nlpspd, &
      proj, ngatherarr, nscatterarr, GPU, irrzon, phnons, pkernel, rxyz, fion, fdisp, psi, phi, coeff, rhopot, fxyz, fnoise)


  ! Deallocate all arrays related to the linear scaling version.
  call deallocateLinear(iproc, lin, phi, lphi, coeff)



end subroutine linearScaling




subroutine mixPotential(iproc, n3p, Glr, input, lin, rhopotOld, rhopot, pnrm)
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
type(linearParameters),intent(in):: lin
real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin),intent(in):: rhopotOld
real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin),intent(in out):: rhopot
real(8),intent(out):: pnrm

! Local variables
integer:: i, ierr
real(8):: tt


  pnrm=0.d0
  tt=1.d0-lin%alphaMix
  !do i=1,max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin
  do i=1,max(Glr%d%n1i*Glr%d%n2i*n3p,1)
      pnrm=pnrm+(rhopot(i)-rhopotOld(i))**2
      rhopot(i)=tt*rhopotOld(i)+lin%alphaMix*rhopot(i)
  end do
  call mpiallred(pnrm, 1, mpi_sum, mpi_comm_world, ierr)
  pnrm=sqrt(pnrm)/(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*input%nspin)

end subroutine mixPotential




subroutine printSummary(iproc, itSCC, infoBasisFunctions, infoCoeff, pnrm, energy, mixingMethod)
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
real(8),intent(in):: pnrm, energy
character(len=4),intent(in):: mixingMethod

  if(iproc==0) then
      write(*,'(x,a)') repeat('#',66 + int(log(real(itSCC))/log(10.)))
      write(*,'(x,a,i0,a)') 'at iteration ', itSCC, ' of the self consistency cycle:'
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
          write(*,'(3x,a,3x,i0,es11.2,es27.17)') 'it, Delta DENS, energy ', itSCC, pnrm, energy
      else if(mixingMethod=='pot') then
          write(*,'(3x,a,3x,i0,es11.2,es27.17)') 'it, Delta POT, energy ', itSCC, pnrm, energy
      end if
      write(*,'(x,a)') repeat('#',66 + int(log(real(itSCC))/log(10.)))
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
