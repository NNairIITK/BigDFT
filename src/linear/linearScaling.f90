subroutine linearScaling(iproc, nproc, n3d, n3p, i3s, i3xcsh, Glr, orbs, comms, at, input, lin, rxyz, fion, fdisp, radii_cf, &
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
integer,intent(in):: iproc, nproc, n3d, n3p, i3s, i3xcsh
type(locreg_descriptors),intent(in) :: Glr
type(orbitals_data),intent(in):: orbs
type(communications_arrays),intent(in) :: comms
type(atoms_data),intent(in):: at
type(linearParameters),intent(in out):: lin
type(input_variables),intent(in):: input
real(8),dimension(3,at%nat),intent(in):: rxyz, fion, fdisp
real(8),dimension(at%ntypes,3),intent(in):: radii_cf
integer,dimension(0:nproc-1,4),intent(inout):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
integer,dimension(0:nproc-1,2),intent(in):: ngatherarr
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
integer:: infoBasisFunctions
real(8),dimension(:),allocatable:: phi
real(8),dimension(:,:),allocatable:: occupForInguess
real(8):: ebsMod


  if(iproc==0) then
      write(*,'(x,a)') repeat('*',84)
      write(*,'(x,a)') '****************************** LINEAR SCALING VERSION ******************************'
      write(*,'(x,a)') '********* Use the selfconsistent potential for the linear scaling version. *********'
  end if
  allocate(occupForInguess(32,at%nat))

  ! Initialize the parameters for the linear scaling version. This will not affect the parameters for the cubic version.
  call allocateAndInitializeLinear(iproc, nproc, Glr, orbs, at, lin, phi, input, rxyz, occupForInguess)

  ! The next subroutine will create the variable wave function descriptors.
  ! It is not used at the moment.
  !!$if(iproc==0) write(*,'(x,a)') ' ~~~~~~~ Creating the variable wave function descriptors and testing them... ~~~~~~~'
  !!$call initializeLocRegLIN(iproc, nproc, Glr, lin, at, input, rxyz, radii_cf)
  !!$if(iproc==0) write(*,'(x,a)') '~~~~~~~~~~~~~~~~~~~~~~~ Descriptors created and test passed. ~~~~~~~~~~~~~~~~~~~~~~~'


  if(nproc==1) allocate(psit(size(psi)))
  ! This subroutine gives back the new psi and psit, which are a linear combination of localized basis functions.
  call getLinearPsi(iproc, nproc, input%nspin, Glr, orbs, comms, at, lin, rxyz, rxyz, &
      nscatterarr, ngatherarr, nlpspd, proj, rhopot, GPU, input, pkernelseq, phi, psi, psit, &
      infoBasisFunctions, n3p, n3d, irrzon, phnons, pkernel, pot_ion, rhocore, potxc, PSquiet, ebsMod)

  ! Calculate the energy that we get with psi.
  call potentialAndEnergySub(iproc, nproc, n3d, n3p, Glr, orbs, at, input, lin, psi, rxyz, &
      rhopot, nscatterarr, ngatherarr, GPU, irrzon, phnons, pkernel, pot_ion, rhocore, potxc, PSquiet, &
      proj, nlpspd, pkernelseq, eion, edisp, eexctX, scpot, ebsMod, energy)


  ! Calculate the forces we get with psi.
  call calculateForcesSub(iproc, nproc, n3p, i3s, i3xcsh, Glr, orbs, at, input, lin, nlpspd, proj, ngatherarr, nscatterarr, GPU, &
      irrzon, phnons, pkernel, rxyz, fion, fdisp, psi, fxyz, fnoise)

  ! Deallocate all arrays related to the linear scaling version.
  call deallocateLinear(iproc, lin, phi)


end subroutine linearScaling
