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
type(atoms_data),intent(in):: at
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

! new variables
type(locreg_descriptors),dimension(:),pointer:: Llr
integer:: nlr
integer,dimension(:,:),pointer :: outofzone


! Local variables
integer:: infoBasisFunctions, infoCoeff, istat, iall, itSCC, nitSCC, i, ierr
real(8),dimension(:),allocatable:: phi, phid
real(8),dimension(:,:),allocatable:: occupForInguess, coeff, coeffd
real(8):: ebsMod, alpha, pnrm, tt
character(len=*),parameter:: subname='linearScaling'
real(8),dimension(:),allocatable:: rhopotOld
type(linearParameters):: lind

integer,dimension(:,:),allocatable:: nscatterarrTemp !n3d,n3p,i3s+i3xcsh-1,i3xcsh
real(8),dimension(:),allocatable:: phiTemp
real(wp),dimension(:),allocatable:: projTemp




  if(iproc==0) then
      write(*,'(x,a)') repeat('*',84)
      write(*,'(x,a)') '****************************** LINEAR SCALING VERSION ******************************'
      write(*,'(x,a)') '********* Use the selfconsistent potential for the linear scaling version. *********'
  end if
  allocate(occupForInguess(32,at%nat))

  ! Initialize the parameters for the linear scaling version. This will not affect the parameters for the cubic version.
  call allocateAndInitializeLinear(iproc, nproc, Glr, orbs, at, lin, lind, phi, phid, &
      input, rxyz, occupForInguess, coeff, coeffd, nlr, Llr, outofzone)

  ! The next subroutine will create the variable wave function descriptors.
  ! It is not used at the moment.
  !!$if(iproc==0) write(*,'(x,a)') ' ~~~~~~~ Creating the variable wave function descriptors and testing them... ~~~~~~~'
  !!$call initializeLocRegLIN(iproc, nproc, Glr, lin, at, input, rxyz, radii_cf)
  !!$if(iproc==0) write(*,'(x,a)') '~~~~~~~~~~~~~~~~~~~~~~~ Descriptors created and test passed. ~~~~~~~~~~~~~~~~~~~~~~~'


  if(nproc==1) allocate(psit(size(psi)))
  nitSCC=lin%nitSCC
  alpha=.1d0
  allocate(rhopotOld(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin), stat=istat)
  call memocc(istat, rhopotOld, 'rhopotOld', subname)
  do itSCC=1,nitSCC
      ! This subroutine gives back the new psi and psit, which are a linear combination of localized basis functions.
      call getLinearPsi(iproc, nproc, input%nspin, nlr, LLr, Glr, orbs, comms, at, lin, lind, rxyz, rxyz, &
          nscatterarr, ngatherarr, nlpspd, proj, rhopot, GPU, input, pkernelseq, phi, phid, psi, psit, &
          infoBasisFunctions, infoCoeff, itScc, n3p, n3pi, n3d, irrzon, phnons, pkernel, pot_ion, rhocore, potxc, PSquiet, &
          i3s, i3xcsh, fion, fdisp, fxyz, eion, edisp, fnoise, ebsMod, coeff, coeffd)

      ! Calculate the energy that we get with psi.
      call dcopy(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin, rhopot(1), 1, rhopotOld(1), 1)
      call potentialAndEnergySub(iproc, nproc, n3d, n3p, nlr, Llr, Glr, orbs, at, input, lin, phi, psi, rxyz, rxyz, &
          rhopot, nscatterarr, ngatherarr, GPU, irrzon, phnons, pkernel, pot_ion, rhocore, potxc, PSquiet, &
          proj, nlpspd, pkernelseq, eion, edisp, eexctX, scpot, coeff, ebsMod, energy)

      !!! TEST  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Calculate the forces we get with psi.
        allocate(nscatterarrTemp(0:nproc-1,4), stat=istat)
        call memocc(istat, nscatterarrTemp, 'nscatterarrTemp', subname)
        allocate(phiTemp(size(phi)), stat=istat)
        call memocc(istat, phiTemp, 'phiTemp', subname)
        allocate(projTemp(nlpspd%nprojel), stat=istat)
        call memocc(istat, projTemp, 'projTemp', subname)
        projTemp=proj
        nscatterarrTemp=nscatterarr
        phiTemp=phi
        call calculateForcesSub(iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh, Glr, orbs, at, input, comms, lin, nlpspd, &
            proj, ngatherarr, nscatterarr, GPU, irrzon, phnons, pkernel, rxyz, fion, fdisp, psi, phi, coeff, fxyz, fnoise)
        proj=projTemp
        nscatterarr=nscatterarrTemp
        phi=phiTemp
        iall=-product(shape(nscatterarrTemp))*kind(nscatterarrTemp)
        deallocate(nscatterarrTemp, stat=istat)
        call memocc(istat, iall, 'nscatterarrTemp', subname)
        iall=-product(shape(phiTemp))*kind(phiTemp)
        deallocate(phiTemp, stat=istat)
        call memocc(istat, iall, 'phiTemp', subname)
        iall=-product(shape(projTemp))*kind(projTemp)
        deallocate(projTemp, stat=istat)
        call memocc(istat, iall, 'projTemp', subname)
      !!!  TEST  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Mix the potential
      call mixPotential(iproc, n3p, Glr, input, lin, rhopotOld, rhopot, pnrm)

      ! Write some informations
      call printSummary(iproc, itSCC, infoBasisFunctions, infoCoeff, pnrm, energy)
  end do
  iall=-product(shape(rhopotOld))*kind(rhopotOld)
  deallocate(rhopotOld, stat=istat)
  call memocc(istat, iall, 'rhopotOld', subname)


  ! Calculate the forces we get with psi.
  call calculateForcesSub(iproc, nproc, n3d, n3p, n3pi, i3s, i3xcsh, Glr, orbs, at, input, comms, lin, nlpspd, &
      proj, ngatherarr, nscatterarr, GPU, irrzon, phnons, pkernel, rxyz, fion, fdisp, psi, phi, coeff, fxyz, fnoise)

  ! Deallocate all arrays related to the linear scaling version.
  call deallocateLinear(iproc, lin, lind, phi, coeff, phid, coeffd)


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




subroutine printSummary(iproc, itSCC, infoBasisFunctions, infoCoeff, pnrm, energy)
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
      else
          write(*,'(3x,a,i0,a)') '- coefficients converged in ', infoCoeff, ' iterations.'
      end if
      write(*,'(3x,a,3x,i0,es11.2,es27.17)') 'it, Delta POT, energy ', itSCC, pnrm, energy
      write(*,'(x,a)') repeat('#',66 + int(log(real(itSCC))/log(10.)))
  end if

end subroutine printSummary
