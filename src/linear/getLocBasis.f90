subroutine getLinearPsi(iproc, nproc, nspin, Glr, orbs, orbsLIN, comms, commsLIN, at, rxyz, rxyzParab, &
    nscatterarr, ngatherarr, nlpspd, proj, rhopot, GPU, input, pkernelseq, phi, psi, psit, &
    infoBasisFunctions, n3p)
!
! Purpose:
! ========
!   This subroutine creates the orbitals psi out of a linear combination of localized basis functions
!   phi. To do so, it proceeds as follows:
!    1. Create the basis functions (with subroutine 'getLocalizedBasis')
!    2. Write the Hamiltonian in this new basis.
!    3. Diagonalize this Hamiltonian matrix.
!    4. Build the new linear combinations. 
!   The basis functions are localized by adding a confining quartic potential to the ordinary DFT 
!   Hamiltonian. There is no self consistency cycle for the potential, i.e. the basis functionsi
!   are optimized with a fixed potential.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc           process ID
!     nproc           total number of processes
!     nspin           npsin==1 -> closed shell; npsin==2 -> spin polarized
!     Glr             type describing the localization region
!     orbs            type describing the physical orbitals psi
!     orbsLIN         type describing the basic functions phi
!     comms           type containing the communication parameters for the physical orbitals psi
!     commsLIN        type containing the communication parameters for the basis functions psi
!     at              type containing the paraneters for the atoms
!     rxyz            the atomic positions
!     rxyzParab       the center of the confinement potential (at the moment identical rxyz)
!     nscatterarr     ???
!     ngatherarr      ???
!     nlpsp           ???
!     proj            ???
!     rhopot          the charge density
!     GPU             parameters for GPUs
!     input           type containing some very general parameters
!     pkernelseq      ???
!     n3p             ???
!  Input/Output arguments
!  ---------------------
!     phi             the localized basis functions. It is assumed that they have been initialized
!                     somewhere else
!   Output arguments
!   ----------------
!     psi             the physical orbitals, which will be a linear combinations of the localized
!                     basis functions phi
!     psit            psi transposed
!     infoBasisFunctions  indicated wheter the basis functions converged to the specified limit (value is 0)
!                         or whether the iteration stopped due to the iteration limit (value is -1). This info
!                         is returned by 'getLocalizedBasis'
!
use module_base
use module_types
use module_interfaces, exceptThisOne => getLinearPsi

! Calling arguments
integer,intent(in):: iproc, nproc, nspin, n3p
type(locreg_descriptors), intent(in) :: Glr
type(orbitals_data), intent(in) :: orbs, orbsLIN
type(communications_arrays), intent(in) :: comms
type(communications_arrays), intent(in) :: commsLIN
type(atoms_data), intent(in) :: at
type(input_variables),intent(in):: input
real(8),dimension(3,at%nat),intent(in):: rxyz, rxyzParab
integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
integer,dimension(0:nproc-1,2),intent(in) :: ngatherarr
type(nonlocal_psp_descriptors),intent(in) :: nlpspd
real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
!real(dp), dimension(*), intent(inout) :: rhopot
!real(dp), dimension(sizeRhopot), intent(inout) :: rhopot
real(dp), dimension(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin) :: rhopot
type(GPU_pointers):: GPU
real(dp),dimension(:),pointer,intent(in) :: pkernelseq
real(8),dimension(orbsLIN%npsidim),intent(in out):: phi
real(8),dimension(orbs%npsidim),intent(out):: psi, psit
integer,intent(out):: infoBasisFunctions

! Local variables 
integer:: istat 
real(8),dimension(:),allocatable:: hphi, eval 
real(8),dimension(:,:),allocatable:: HamSmall 
real(8),dimension(:),pointer:: phiWorkPointer 
real(8):: epot_sum,ekin_sum,eexctX,eproj_sum, ddot, trace 
real(wp), dimension(:), pointer :: potential 
character(len=*),parameter:: subname='getLinearPsi' 
 


allocate(hphi(orbsLIN%npsidim), stat=istat) 
call memocc(istat, hphi, 'hphi', subname)
allocate(phiWorkPointer(max(size(phi),size(psi))), stat=istat)
call memocc(istat, phiWorkPointer, 'phiWorkPointer', subname)

call getLocalizedBasis(iproc, nproc, at, orbs, Glr, input, orbsLIN, commsLIN, rxyz, nspin, nlpspd, proj, &
    nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, hphi, trace, rxyzParab, &
    orbsLIN%DIISHistMin, orbsLIN%DIISHistMax, infoBasisFunctions)

!allocate the potential in the full box
call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*n3p,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,input%nspin,&
     orbs%norb,orbs%norbp,ngatherarr,rhopot,potential)

call HamiltonianApplication(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
     nlpspd,proj,Glr,ngatherarr,potential,&
     phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel=pkernelseq)

!deallocate potential
call free_full_potential(nproc,potential,subname)


call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, hphi, work=phiWorkPointer)

allocate(HamSmall(orbsLIN%norb,orbsLIN%norb), stat=istat)
call transformHam(iproc, nproc, orbsLIN, commsLIN, phi, hphi, HamSmall)

if(iproc==0) write(*,'(a)', advance='no') 'Linear Algebra... '
allocate(eval(orbsLIN%norb), stat=istat)
call diagonalizeHamiltonian(iproc, nproc, orbsLIN, HamSmall, eval)

call buildWavefunction(iproc, nproc, orbs, orbsLIN, comms, commsLIN, phi, psi, HamSmall)

call dcopy(orbs%npsidim, psi, 1, psit, 1)
if(iproc==0) write(*,'(a)') 'done.'


call untranspose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, psi, work=phiWorkPointer)


iall=-product(shape(phiWorkPointer))*kind(phiWorkPointer)
deallocate(phiWorkPointer, stat=istat)
call memocc(istat, iall, 'phiWorkPointer', subname)

iall=-product(shape(hphi))*kind(hphi)
deallocate(hphi, stat=istat)
call memocc(istat, iall, 'hphi', subname)



end subroutine getLinearPsi
