subroutine linearScaling(iproc, nproc, Glr, orbs, comms, at, input, lin, rxyz, nscatterarr, ngatherarr, &
    nlpspd, proj, rhopot, GPU, pkernelseq, psi, psit, radii_cf, n3d, n3p, i3s, i3xcsh, irrzon, phnons, pkernel, pot_ion, &
    rhocore, potxc, PSquiet, eion, edisp, eexctX, scpot, fxyz, fion, fdisp)
!
! Purpose:
! ========
!   Top-level subroutine for the linear scaling version.
!
! Calling arguments:
! ==================

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
type(linearParameters):: lin
type(input_variables),intent(in):: input
real(8),dimension(3,at%nat):: rxyz, fxyz, fion, fdisp
integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
integer,dimension(0:nproc-1,2),intent(in) :: ngatherarr
type(nonlocal_psp_descriptors),intent(in) :: nlpspd
real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
real(dp), dimension(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin) :: rhopot
type(GPU_pointers):: GPU
real(dp),dimension(:),pointer,intent(in) :: pkernelseq
real(8),dimension(orbs%npsidim),intent(out):: psi
real(8),dimension(:),pointer:: psit
real(8),dimension(at%ntypes,3),intent(in):: radii_cf
integer, dimension(lin%as%size_irrzon(1),lin%as%size_irrzon(2),lin%as%size_irrzon(3)) :: irrzon 
real(dp), dimension(lin%as%size_phnons(1),lin%as%size_phnons(2),lin%as%size_phnons(3)) :: phnons 
real(dp), dimension(lin%as%size_pkernel):: pkernel
real(wp), dimension(lin%as%size_pot_ion):: pot_ion
!real(wp), dimension(lin%as%size_rhocore):: rhocore 
real(wp), dimension(:),pointer:: rhocore                  
real(wp), dimension(lin%as%size_potxc(1),lin%as%size_potxc(2),lin%as%size_potxc(3),lin%as%size_potxc(4)):: potxc
character(len=3):: PSquiet
real(gp):: eion, edisp, eexctX
logical:: scpot



! Local variables
integer:: infoBasisFunctions
real(8),dimension(:),allocatable:: phi
real(8),dimension(:,:),allocatable:: rxyzParab
real(8),dimension(:,:),allocatable:: occupForInguess


  if(iproc==0) then
      write(*,'(x,a)') repeat('*',84)
      write(*,'(x,a)') '****************************** LINEAR SCALING VERSION ******************************'
      write(*,'(x,a)') '********* Use the selfconsistent potential for the linear scaling version. *********'
  end if
  ! Initialize the parameters for the linear scaling version. This will not affect the parameters for the cubic version.
  allocate(occupForInguess(32,at%nat))
  call allocateAndInitializeLinear(iproc, nproc, Glr, orbs, at, lin, phi, input, rxyz, occupForInguess)
  if(iproc==0) write(*,'(x,a)') ' ~~~~~~~ Creating the variable wave function descriptors and testing them... ~~~~~~~'
  call initializeLocRegLIN(iproc, nproc, Glr, lin, at, input, rxyz, radii_cf)
  if(iproc==0) write(*,'(x,a)') '~~~~~~~~~~~~~~~~~~~~~~~ Descriptors created and test passed. ~~~~~~~~~~~~~~~~~~~~~~~'
  !if(iproc==0) write(*,'(a)') 'trying to reproduce the result with the linear scaling version...'
  if(nproc==1) allocate(psit(size(psi)))
  if(.not.allocated(rxyzParab)) allocate(rxyzParab(3,at%nat))
  rxyzParab=rxyz
  ! This subroutine gives back the new psi and psit, which are a linear combination of localized basis functions.
  call getLinearPsi(iproc, nproc, input%nspin, Glr, orbs, comms, at, lin, rxyz, rxyzParab, &
      nscatterarr, ngatherarr, nlpspd, proj, rhopot, GPU, input, pkernelseq, phi, psi, psit, &
      infoBasisFunctions, n3p)

call potentialAndEnergySub(iproc, nproc, Glr, orbs, at, input, lin, psi, rhopot, &
    nscatterarr, ngatherarr, GPU, irrzon, phnons, pkernel, pot_ion, rhocore, potxc, PSquiet, &
    proj, nlpspd, pkernelseq, rxyz, eion, edisp, eexctX, scpot, n3d, n3p)

  ! Calculate the potential arising from the new psi and calculate the energy.
  !!$call potentialAndEnergy()

  !nscatterarrCorrect=nscatterarr
  !ngatherarrCorrect=ngatherarr
  !projCorrect=proj
  !fxyzOld=fxyz

call calculateForcesSub(iproc, nproc, Glr, orbs, at, input, lin, nlpspd, proj, ngatherarr, nscatterarr, GPU, &
    irrzon, phnons, pkernel, rxyz, fxyz, fion, fdisp, n3p, i3s, i3xcsh, psi)

  ! Calculate the forces arising from the new psi.
  !!$call calculateForces()

  !!$call deallocateLinear(iproc, lin, phi)





end subroutine linearScaling
