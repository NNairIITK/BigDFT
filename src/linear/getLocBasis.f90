subroutine getLinearPsi(iproc, nproc, nspin, Glr, orbs, comms, at, lin, rxyz, rxyzParab, &
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
!     comms           type containing the communication parameters for the physical orbitals psi
!     at              type containing the paraneters for the atoms
!     lin             type containing parameters for the linear version
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
type(orbitals_data), intent(in) :: orbs
type(communications_arrays), intent(in) :: comms
type(atoms_data), intent(in) :: at
type(linearParameters):: lin
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
real(8),dimension(lin%orbs%npsidim),intent(in out):: phi
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
 


allocate(hphi(lin%orbs%npsidim), stat=istat) 
call memocc(istat, hphi, 'hphi', subname)
allocate(phiWorkPointer(max(size(phi),size(psi))), stat=istat)
call memocc(istat, phiWorkPointer, 'phiWorkPointer', subname)



call getLocalizedBasis(iproc, nproc, at, orbs, Glr, input, lin, rxyz, nspin, nlpspd, proj, &
    nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, hphi, trace, rxyzParab, &
    infoBasisFunctions)

!allocate the potential in the full box
call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*n3p,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,input%nspin,&
     lin%orbs%norb,lin%orbs%norbp,ngatherarr,rhopot,potential)

call HamiltonianApplication(iproc,nproc,at,lin%orbs,input%hx,input%hy,input%hz,rxyz,&
     nlpspd,proj,Glr,ngatherarr,potential,&
     phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel=pkernelseq)

!deallocate potential
call free_full_potential(nproc,potential,subname)


call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWorkPointer)
call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, hphi, work=phiWorkPointer)

allocate(HamSmall(lin%orbs%norb,lin%orbs%norb), stat=istat)
call transformHam(iproc, nproc, lin%orbs, lin%comms, phi, hphi, HamSmall)

if(iproc==0) write(*,'(a)', advance='no') 'Linear Algebra... '
allocate(eval(lin%orbs%norb), stat=istat)
call diagonalizeHamiltonian(iproc, nproc, lin%orbs, HamSmall, eval)

call buildWavefunction(iproc, nproc, orbs, lin%orbs, comms, lin%comms, phi, psi, HamSmall)

call dcopy(orbs%npsidim, psi, 1, psit, 1)
if(iproc==0) write(*,'(a)') 'done.'


call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWorkPointer)
call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, psi, work=phiWorkPointer)


iall=-product(shape(phiWorkPointer))*kind(phiWorkPointer)
deallocate(phiWorkPointer, stat=istat)
call memocc(istat, iall, 'phiWorkPointer', subname)

iall=-product(shape(hphi))*kind(hphi)
deallocate(hphi, stat=istat)
call memocc(istat, iall, 'hphi', subname)



end subroutine getLinearPsi




subroutine initializeParameters(iproc, nproc, Glr, orbs, at, lin, phi, input, rxyz, occupForInguess)
!
! Purpose:
! ========
!   This subroutine initializes all parameters needed for the linear scaling version.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc           process ID
!     nproc           total number of processes
!     Glr             type describing the localization region
!     orbs            type describing the physical orbitals psi
!     at              type containing the paraneters for the atoms
!     lin             type containing parameters for the linear version
!     input           type containing some very general parameters
!     rxyz            the atomic positions
!     occuprForINguess  delete maybe
!  Output arguments
!  ---------------------
!     phi             the localized basis functions. They are only initialized here, but
!                       not normalized.
!
use module_base
use module_types
use module_interfaces, except_this_one => initializeParameters
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(locreg_descriptors), intent(in) :: Glr
type(orbitals_data), intent(inout) :: orbs
type(atoms_data), intent(in) :: at
type(linearParameters):: lin
type(input_variables), intent(in) :: input
real(8),dimension(3,at%nat),intent(in):: rxyz
real(8),dimension(32,at%nat):: occupForInguess
real(8),dimension(:),allocatable,intent(out):: phi

! Local variables
integer:: jproc, istat, iorb, jorb, ierr, iat, ityp, iall, norb_tot
integer:: norb, norbu, norbd
integer,dimension(:),allocatable:: norbsPerType, norbsPerAtom
character(len=*),parameter:: subname='initializeParameters'
character(len=20):: atomname
logical:: written


! Allocate all local arrays.
allocate(norbsPerType(at%ntypes), stat=istat)
call memocc(istat, at%ntypes, 'norbsPerType', subname)
allocate(norbsPerAtom(at%nat), stat=istat)
call memocc(istat, at%nat, 'norbsPerAtom', subname)

! Read in all parameters related to the linear scaling version and print them.
allocate(lin%potentialPrefac(at%ntypes), stat=istat)
call memocc(istat, at%ntypes, 'lin%potentialPrefac', subname)
open(unit=99, file='input.lin')
read(99,*) lin%nItMax
read(99,*) lin%convCrit
read(99,*) lin%DIISHistMin, lin%DIISHistMax
if(lin%DIISHistMin>lin%DIISHistMax) then
    if(iproc==0) write(*,'(a,i0,a,i0,a)') 'ERROR: DIISHistMin must not be larger than &
    & DIISHistMax, but you chose ', lin%DIISHistMin, ' and ', lin%DIISHistMax, '!'
    stop
end if
if(iproc==0) write(*,'(x,a)') '================== Input parameters. =================='
if(iproc==0) write(*,'(x,a,9x,a,3x,a,3x,a,4x,a,4x,a)') '| ', ' | ', 'number of', ' | ', 'prefactor for', ' |'
if(iproc==0) write(*,'(x,a,a,a,a,a,a,a)') '| ', 'atom type', ' | ', 'basis functions', ' | ', 'confinement potential', ' |'
do iat=1,at%ntypes
    read(99,*) atomname, norbsPerType(iat), lin%potentialPrefac(iat)
    if(iproc==0) write(*,'(x,a,4x,a,a,a,a,i0,7x,a,7x,es9.3,6x,a)') '| ', trim(atomname), repeat(' ', 6-len_trim(atomname)), '|',  &
        repeat(' ', 10-ceiling(log10(dble(norbsPerType(iat)+1)))), norbsPerType(iat), '|', lin%potentialPrefac(iat), ' |'
end do
if(iproc==0) write(*,'(x,a)') '--------------------------------------------------------'
if(iproc==0) write(*,'(x,a,a,a,a,a,a,a)') '| ', 'maximal number', ' | ', 'convergence', ' | ', 'DIIS history', ' |'
if(iproc==0) write(*,'(x,a,x,a,a,x,a,x,a,x,a,2x,a)') '| ', 'of iterations', ' | ', 'criterion', ' | ', 'min   max', ' |'
if(iproc==0) write(*,'(x,a,a,i0,5x,a,x,es9.3,x,a,a,i0,3x,a,i0,2x,a)') '| ', repeat(' ', 9-ceiling(log10(dble(lin%nItMax+1)))), lin%nItMax, &
    ' | ', lin%convCrit, ' | ', repeat(' ', 4-ceiling(log10(dble(lin%DIISHistMin+1)))), lin%DIISHistMin, &
    repeat(' ', 4-ceiling(log10(dble(lin%nItMax+1)))), lin%DIISHistMax, ' |'
if(iproc==0) write(*,'(x,a)') '======================================================='
close(unit=99)


! Assign to each atom its number of basis functions and count how many basis functions 
! we have in total.
do iat=1,at%nat
    ityp=at%iatype(iat)
    norbsPerAtom(iat)=norbsPerType(ityp)
    lin%orbs%norb=lin%orbs%norb+norbsPerAtom(iat)
end do


! Distribute the orbitals among the processors.
norb=lin%orbs%norb
norbu=norb
norbd=0
call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, orbs%nspinor, input%nkpt, input%kpt, input%wkpt, lin%orbs)
written=.false.
do jproc=1,nproc-1
    if(lin%orbs%norb_par(jproc)<lin%orbs%norb_par(jproc-1)) then
        if(iproc==0) write(*,'(x,a,5(i0,a))') 'Processes from 0 to ',jproc-1,' treat ',lin%orbs%norb_par(jproc-1), &
            ' orbitals, processes from ',jproc,' to ',nproc-1,' treat ',lin%orbs%norb_par(jproc),' orbitals.'
        written=.true.
        exit
    end if
end do
if(.not.written) then
    if(iproc==0) write(*,'(x,a,2(i0,a))') 'Processes from 0 to ',nproc-1,' treat ',lin%orbs%norbp,' orbitals.'
end if


! Decide which orbital is centered in which atom.
allocate(lin%onWhichAtom(lin%orbs%norbp), stat=istat)
call memocc(istat, nproc, 'lin%onWhichAtomnorb_par', subname)
call assignOrbitalsToAtoms(iproc, at%nat, lin, norbsPerAtom)

!write(*,'(a,i2,3x,200i4)') 'iproc, lin%onWhichAtom', iproc, lin%onWhichAtom


! orbsLIN%isorb is the 'first' orbital for a given MPI process.
norb_tot=0
do jproc=0,iproc-1
   norb_tot=norb_tot+lin%orbs%norb_par(jproc)
end do
!reference orbital for process
lin%orbs%isorb=norb_tot


if(associated(lin%orbs%eval)) then
    nullify(lin%orbs%eval)
end if
allocate(lin%orbs%eval(lin%orbs%norb), stat=istat)
lin%orbs%eval=-.5d0


! Assign the parameters needed for the communication to lin%comms
call orbitals_communicators(iproc,nproc,Glr,lin%orbs,lin%comms)


! Allocate phi and initialize it at random
allocate(phi(lin%orbs%npsidim), stat=istat)
call memocc(istat, lin%orbs%npsidim, 'phi', subname)
call initRandomSeed(iproc, 1)
call random_number(phi)

!write(*,*) 'calling createInputGuess'
!call createInputGuess(iproc, orbsLIN, Glr, input, at, rxyz, phi)


! Deallocate all local arrays
iall=-product(shape(norbsPerType))*kind(norbsPerType)
deallocate(norbsPerType, stat=istat)
call memocc(istat, iall, 'norbsPerType', subname)
iall=-product(shape(norbsPerAtom))*kind(norbsPerAtom)
deallocate(norbsPerAtom, stat=istat)
call memocc(istat, iall, 'norbsPerAtom', subname)


end subroutine initializeParameters



subroutine assignOrbitalsToAtoms(iproc, nat, lin, norbsPerAt)
!
! Purpose:
! ========
!   Assigns the orbitals to the atoms, using the array lin%onWhichAtom.
!   If orbital i is centered on atom j, we have lin%onWhichAtom(i)=j.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc         process ID
!     nat           number of atoms
!     norbsPerAt    indicates how many orbitals are centered on each atom.
!  Input / Output arguments
!  ---------------------
!     lin           type containing parameters for the linear version
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nat
type(linearParameters):: lin
integer,dimension(nat):: norbsPerAt

! Local variables
integer:: jproc, iiOrb, iorb, jorb, jat


! There are four counters:
!   jproc: indicates which MPI process is handling the basis function which is being treated
!   jat: counts the atom numbers
!   jorb: counts the orbitals handled by a given process
!   iiOrb: counts the number of orbitals for a given atoms thas has already been assigned
jproc=0
jat=1
jorb=0
iiOrb=0

do iorb=1,lin%orbs%norb

    ! Switch to the next MPI process if the numbers of orbitals for a given
    ! MPI process is reached.
    if(jorb==lin%orbs%norb_par(jproc)) then
        jproc=jproc+1
        jorb=0
    end if
    
    ! Switch to the next atom if the number of basis functions for this atom is reached.
    if(iiOrb==norbsPerAt(jat)) then
        jat=jat+1
        iiOrb=0
    end if
    jorb=jorb+1
    iiOrb=iiOrb+1
    if(iproc==jproc) lin%onWhichAtom(jorb)=jat
end do    

end subroutine assignOrbitalsToAtoms





subroutine getLocalizedBasis(iproc, nproc, at, orbs, Glr, input, lin, rxyz, nspin, nlpspd, &
    proj, nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, hphi, trH, rxyzParabola, &
    infoBasisFunctions)
!
! Purpose:
! ========
!   Calculates the localized basis functions phi. These basis functions are obtained by adding a
!   quartic potential centered on the atoms to the ordinary Hamiltonian. The eigenfunctions are then
!   determined by minimizing the trace until the gradient norm is below the convergence criterion.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc           process ID
!     nproc           total number of processes
!     at              type containing the paraneters for the atoms
!     orbs            type describing the physical orbitals psi
!     Glr             type describing the localization region
!     input           type containing some very general parameters
!     lin             type containing parameters for the linear version
!     rxyz            the atomic positions
!     nspin           npsin==1 -> closed shell; npsin==2 -> spin polarized
!     nlpsp           ???
!     proj            ???
!     nscatterarr     ???
!     ngatherarr      ???
!     rhopot          the charge density
!     GPU             parameters for GPUs
!     pkernelseq      ???
!     rxyzParab       the center of the confinement potential (at the moment identical rxyz)
!     n3p             ???
!  Input/Output arguments
!  ---------------------
!     phi             the localized basis functions. It is assumed that they have been initialized
!                     somewhere else
!   Output arguments
!   ----------------
!     hphi            the modified Hamiltonian applied to phi
!     trH             the trace of the Hamiltonian
!     infoBasisFunctions  indicated wheter the basis functions converged to the specified limit (value is 0)
!                         or whether the iteration stopped due to the iteration limit (value is -1). This info
!                         is returned by 'getLocalizedBasis'
!
! Calling arguments:
!   Input arguments
!   Output arguments
!    phi   the localized basis functions
!
use module_base
use module_types
use module_interfaces, except_this_one => getLocalizedBasis
  use Poisson_Solver
!use allocModule
implicit none

! Calling arguments
integer:: iproc, nproc, infoBasisFunctions
type(atoms_data), intent(in) :: at
type(orbitals_data):: orbs
type(locreg_descriptors), intent(in) :: Glr
type(input_variables):: input
type(linearParameters):: lin
real(8),dimension(3,at%nat):: rxyz, rxyzParabola
integer:: nspin
type(nonlocal_psp_descriptors), intent(in) :: nlpspd
real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
real(dp), dimension(*), intent(inout) :: rhopot
type(GPU_pointers), intent(inout) :: GPU
real(dp), dimension(:), pointer :: pkernelseq
real(8),dimension(lin%orbs%npsidim):: phi, hphi
real(8):: trH

! Local variables
real(8) ::epot_sum, ekin_sum, eexctX, eproj_sum
real(8):: tt, ddot, fnrm, fnrmMax, meanAlpha, gnrm, gnrm_zero, gnrmMax
integer:: iorb, icountSDSatur, icountSwitch, idsx, icountDIISFailureTot, icountDIISFailureCons, itBest
integer:: istat, istart, ierr, ii, it, nbasisPerAtForDebug, ncong, iall, nvctrp
real(8),dimension(:),allocatable:: hphiold, alpha, diag, fnrmOldArr
real(8),dimension(:,:),allocatable:: HamSmall, fnrmArr, fnrmOvrlpArr
real(8),dimension(:),pointer:: phiWorkPointer
logical:: quiet, allowDIIS
character(len=*),parameter:: subname='getLocalizedBasis'
type(diis_objects):: diisLIN



call allocateLocalArrays()



icountSDSatur=0
icountSwitch=0
icountDIISFailureTot=0
icountDIISFailureCons=0


! No DIIS in the beginning
call initializeDIISParameters(lin%DIISHistMax)
allowDIIS=.true.
if(allowDIIS) then
else
    diisLIN%idsx=0
    call deallocate_diis_objects(diisLIN,subname)
end if



if(iproc==0) write(*,'(a)') '============================ basis functions creation... ============================'
alpha=1.d-2
iterLoop: do it=1,lin%nItMax
    fnrmMax=0.d0
    fnrm=0.d0

    if (iproc==0) then
        write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
    endif


    ! Orthonormalize the orbitals.
    if(iproc==0) then
        write(*,'(x,a)', advance='no') 'Orthonormalization... '
    end if
    call orthogonalize(iproc, nproc, lin%orbs, lin%comms, Glr%wfd, phi, input)

    ! Untranspose phi
    call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWorkPointer)


    ! Calculate the unconstrained gradient.
    if(iproc==0) then
        write(*,'(a)', advance='no') 'Hamiltonian application... '
    end if
    call HamiltonianApplicationParabola(iproc,nproc,at,lin%orbs,lin,input%hx,input%hy,input%hz,rxyz,&
         nlpspd,proj,Glr,ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
         rhopot(1),&
         phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, rxyzParabola, pkernel=pkernelseq)


    ! Apply the orthoconstraint to the gradient. This subroutine also calculates the trace trH.
    if(iproc==0) then
        write(*,'(a)', advance='no') 'orthoconstraint... '
    end if
    call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, hphi, work=phiWorkPointer)
    call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWorkPointer)
    call orthoconstraintNotSymmetric(iproc, nproc, lin%orbs, lin%comms, Glr%wfd, phi, hphi, trH, diag)


    ! Calculate the norm of the gradient (fnrmArr) and determine the angle between the current gradient and that
    ! of the previous iteration (fnrmOvrlpArr).
    nvctrp=lin%comms%nvctr_par(iproc,1) ! 1 for k-point
    istart=1
    do iorb=1,lin%orbs%norb
        if(it>1) fnrmOvrlpArr(iorb,2)=ddot(nvctrp*orbs%nspinor, hphi(istart), 1, hphiold(istart), 1)
        fnrmArr(iorb,2)=ddot(nvctrp*orbs%nspinor, hphi(istart), 1, hphi(istart), 1)
        istart=istart+nvctrp*orbs%nspinor
    end do
    call mpi_allreduce(fnrmArr(1,2), fnrmArr(1,1), lin%orbs%norb, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
    call mpi_allreduce(fnrmOvrlpArr(1,2), fnrmOvrlpArr(1,1), lin%orbs%norb, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)

    ! Keep the gradient for the next iteration.
    if(it>1) then
        call dcopy(lin%orbs%norb, fnrmArr(1,1), 1, fnrmOldArr(1), 1)
    end if

    ! Determine the gradient norm and its maximal component. In addition, adapt the
    !  step size for the steepest descent minimization (depending on the angle 
    ! between the current gradient and the one from the previous iteration).
    ! This is of course only necessary if we are using steepest descent and not DIIS.
    do iorb=1,lin%orbs%norb
        fnrm=fnrm+fnrmArr(iorb,1)
        if(fnrmArr(iorb,1)>fnrmMax) fnrmMax=fnrmArr(iorb,1)
        if(it>1 .and. diisLIN%idsx==0 .and. .not.diisLIN%switchSD) then
        ! Adapt step size for the steepest descent minimization.
            tt=fnrmOvrlpArr(iorb,1)/sqrt(fnrmArr(iorb,1)*fnrmOldArr(iorb))
            if(tt>.7d0) then
                alpha(iorb)=alpha(iorb)*1.05d0
            else
                alpha(iorb)=alpha(iorb)*.5d0
            end if
        end if
    end do
    fnrm=sqrt(fnrm)
    fnrmMax=sqrt(fnrmMax)
    ! Copy the gradient (will be used in the next iteration to adapt the step size).
    call dcopy(lin%orbs%norb*nvctrp*orbs%nspinor, hphi(1), 1, hphiold(1), 1)

    ! Untranspose hphi.
    call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, hphi, work=phiWorkPointer)


    ! Precondition the gradient
    if(iproc==0) then
        write(*,'(a)') 'preconditioning. '
    end if
    ncong=10
    gnrm=1.d3 ; gnrm_zero=1.d3
    call preconditionallLIN(iproc, nproc, lin%orbs, lin, Glr, input%hx, input%hy, input%hz, &
        ncong, hphi, gnrm, gnrm_zero, gnrmMax,  at%nat, rxyz, at, it)

    tt=gnrm
    call mpi_allreduce(tt, gnrm, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
    gnrm=sqrt(gnrm)
    tt=gnrmMax
    call mpi_allreduce(tt, gnrmMax, 1, mpi_double_precision, mpi_max, mpi_comm_world, ierr)
    gnrmMax=sqrt(gnrmMax)


    tt=sum(alpha)
    meanAlpha=tt/dble(lin%orbs%norb)



    if(iproc==0) write(*,'(x,a,i6,2es15.7,f14.7)') 'iter, fnrm, fnrmMax, trace', it, fnrm, fnrmMax, trH
    if(iproc==0) write(1000,'(i6,2es15.7,f15.7,es12.4)') it, fnrm, fnrmMax, trH, meanAlpha
    if(fnrmMax<lin%convCrit .or. it>=lin%nItMax) then
        if(it>=lin%nItMax) then
            if(iproc==0) write(*,'(a,i0,a)') 'WARNING: not converged within ', it, &
                ' iterations! Exiting loop due to limitations of iterations.'
            if(iproc==0) write(*,'(a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
            infoBasisFunctions=1
        else
            if(iproc==0) then
                write(*,'(a,i0,a,2es15.7,f12.7)') 'converged in ', it, ' iterations.'
                write (*,'(a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
            end if
            infoBasisFunctions=0
        end if
        if(iproc==0) write(*,'(a)') '============================== basis functions created =============================='
        call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWorkPointer)
        !call plotOrbitals(iproc, lin%orbs, Glr, phi, at%nat, rxyz, lin%onWhichAtom, .5d0*input%hx, &
        !    .5d0*input%hy, .5d0*input%hz, 1)
        exit iterLoop
    end if


    call DIISorSD()
    if(iproc==0) then
        if(diisLIN%idsx>0) then
            write(*,'(x,3(a,i0))') 'DIIS informations: history length=',diisLIN%idsx, ', consecutive failures=', &
                icountDIISFailureCons, ', total failures=', icountDIISFailureTot
        else
            write(*,'(x,a,es9.3,a,i0)') 'steepest descent informations: mean alpha=', meanAlpha, &
            ', consecutive successes=', icountSDSatur
        end if
    end if
    if(.not. diisLIN%switchSD) call improve()



 

end do iterLoop








contains

    subroutine initializeDIISParameters(idsxHere)
    ! Purpose:
    ! ========
    !   Initializes all parameters needed for the DIIS procedure.
    !
    ! Calling arguments
    !   idsx    DIIS history length
    !
    implicit none
    
    ! Calling arguments
    integer:: idsxHere

    diisLIN%switchSD=.false.
    diisLIN%idiistol=0
    diisLIN%mids=1
    diisLIN%ids=0
    diisLIN%idsx=idsxHere
    diisLIN%energy_min=1.d10
    diisLIN%energy_old=1.d10
    diisLIN%energy=1.d10
    diisLIN%alpha=2.d0
    call allocate_diis_objects(diisLIN%idsx, lin%orbs%npsidim, 1, diisLIN, subname) ! 1 for k-points
    end subroutine initializeDIISParameters


    subroutine DIISorSD()
    if(diisLIN%switchSD) diisLIN%switchSD=.false.
    ! Determine wheter to use DIIS or SD
    if(trH<=diisLIN%energy_min) then
        ! everything ok
        diisLIN%energy_min=trH
        diisLIN%switchSD=.false.
        itBest=it
        icountSDSatur=icountSDSatur+1
        icountDIISFailureCons=0
        if(icountSDSatur>=10 .and. diisLIN%idsx==0 .and. allowDIIS) then
            ! switch back to DIIS 
            icountSwitch=icountSwitch+1
            idsx=max(lin%DIISHistMin,lin%DIISHistMax-icountSwitch)
            if(iproc==0) write(*,'(a,i0)') 'switch to DIIS with new history length ', idsx
            call initializeDIISParameters(idsx)
            icountDIISFailureTot=0
            icountDIISFailureCons=0
        end if
    else
        ! the trace is growing.
        ! Count how many times this occurs and switch to SD after 3 times.
        icountDIISFailureCons=icountDIISFailureCons+1
        icountDIISFailureTot=icountDIISFailureTot+1
        icountSDSatur=0
        if((icountDIISFailureCons>=2 .or. icountDIISFailureTot>=3) .and. diisLIN%idsx>0) then
            alpha=1.d0
            if(iproc==0) then
                if(icountDIISFailureCons>=2) write(*,'(x,a,i0,a,es10.3)') 'DIIS failed ', &
                    icountDIISFailureCons, ' times consecutievly. Switch to SD with stepsize', alpha(1)
                if(icountDIISFailureTot>=3) write(*,'(x,a,i0,a,es10.3)') 'DIIS failed ', &
                    icountDIISFailureTot, ' times in total. Switch to SD with stepsize', alpha(1)
            end if
            ! Try to get back the orbitals of the best iteration. This is possible if
            ! these orbitals are still present in the DIIS history.
            if(it-itBest<diisLIN%idsx) then
               if(iproc==0) then
                   if(iproc==0) write(*,'(x,a,i0,a)')  'Recover the orbitals from iteration ', &
                       itBest, ' which are the best so far.'
               end if
               ii=modulo(diisLIN%mids-(it-itBest),diisLIN%mids)
               nvctrp=lin%comms%nvctr_par(iproc,1) ! 1 for k-point
               call dcopy(lin%orbs%norb*nvctrp, diisLIN%psidst(ii*nvctrp*lin%orbs%norb+1), 1, phi(1), 1)
            end if
            call deallocate_diis_objects(diisLIN, subname)
            diisLIN%idsx=0
            diisLIN%switchSD=.true.
        end if
    end if
    end subroutine DIISorSD


    subroutine improve()
    ! For DIIS 
    if (diisLIN%idsx > 0) then
       diisLIN%mids=mod(diisLIN%ids,diisLIN%idsx)+1
       diisLIN%ids=diisLIN%ids+1
    end if

    ! Follow the gradient using steepest descent.
    ! The same, but transposed
    call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, hphi, work=phiWorkPointer)
    
    ! steepest descent
    if(diisLIN%idsx==0) then
        istart=1
        nvctrp=lin%comms%nvctr_par(iproc,1) ! 1 for k-point
        do iorb=1,lin%orbs%norb
            call daxpy(nvctrp*orbs%nspinor, -alpha(iorb), hphi(istart), 1, phi(istart), 1)
            istart=istart+nvctrp*orbs%nspinor
        end do
    else
        ! DIIS
        quiet=.true. ! less output
        call psimix(iproc, nproc, lin%orbs, lin%comms, diisLIN, hphi, phi, quiet)
    end if
    end subroutine improve



    subroutine allocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine allocates all local arrays.
    !

    allocate(hphiold(lin%orbs%npsidim), stat=istat)
    call memocc(istat, lin%orbs%npsidim, 'hphiold', subname)

    allocate(alpha(lin%orbs%norb), stat=istat)
    call memocc(istat, lin%orbs%norbp, 'alpha', subname)

    allocate(fnrmArr(lin%orbs%norb,2), stat=istat)
    call memocc(istat, lin%orbs%norb*2, 'fnrmArr', subname)

    allocate(fnrmOldArr(lin%orbs%norb), stat=istat)
    call memocc(istat, lin%orbs%norb, 'fnrmOldArr', subname)

    allocate(fnrmOvrlpArr(lin%orbs%norb,2), stat=istat)
    call memocc(istat, lin%orbs%norb*2, 'fnrmOvrlpArr', subname)

    allocate(phiWorkPointer(size(phi)), stat=istat)
    call memocc(istat, size(phi), 'phiWorkPointer', subname)
    
    allocate(diag(lin%orbs%norb), stat=istat)
    

    end subroutine allocateLocalArrays


    subroutine deallocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine deallocates all local arrays.
    !

    iall=-product(shape(hphiold))*kind(hphiold)
    deallocate(hphiold, stat=istat)
    call memocc(istat, iall, 'hphiold', subname)
    
    iall=-product(shape(alpha))*kind(alpha)
    deallocate(alpha, stat=istat)
    call memocc(istat, iall, 'alpha', subname)

    iall=-product(shape(phiWorkPointer))*kind(phiWorkPointer)
    deallocate(phiWorkPointer, stat=istat)
    call memocc(istat, iall, 'phiWorkPointer', subname)
    
    ! if diisLIN%idsx==0, these arrays have already been deallocated
    if(diisLIN%idsx>0 .and. lin%DIISHistMax>0) call deallocate_diis_objects(diisLIN,subname)

    end subroutine deallocateLocalArrays


end subroutine getLocalizedBasis






subroutine transformHam(iproc, nproc, orbs, comms, phi, hphi, HamSmall)
!
! Purpose:
! =======
!   Builds the Hamiltonian in the basis of the localized basis functions phi. To do so, it gets all basis
!   functions |phi_i> and H|phi_i> and then calculates H_{ij}=<phi_i|H|phi_j>. The basis functions phi are
!   provided in the transposed form.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc      process ID
!     nproc      total number of processes
!     orbs       type describing the basis functions psi
!     comms      type containing the communication parameters for the physical orbitals phi
!     phi        basis functions 
!     hphi       the Hamiltonian applied to the basis functions 
!   Output arguments:
!   -----------------
!     HamSmall   Hamiltonian in small basis
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data), intent(in) :: orbs
type(communications_arrays), intent(in) :: comms
real(8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor,orbs%norb), intent(in) :: phi, hphi
real(8),dimension(orbs%norb,orbs%norb),intent(out):: HamSmall

! Local variables
integer:: istat, ierr, nvctrp, iall
real(8),dimension(:,:),allocatable:: HamTemp
character(len=*),parameter:: subname='transformHam'



! Allocate a temporary array if there are several MPI processes
if(nproc>1) then
    allocate(HamTemp(orbs%norb,orbs%norb), stat=istat)
    call memocc(istat, orbs%norb*orbs%norb, 'HamTemp', subname)
end if

! nvctrp is the amount of each phi hold by the current process
nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor

! Build the Hamiltonian. In the parallel case, each process writes its Hamiltonian in HamTemp
! and a mpi_allreduce sums up the contribution from all processes.
if(nproc==1) then
    call dgemm('t', 'n', orbs%norb, orbs%norb, nvctrp, 1.d0, phi(1,1), nvctrp, &
               hphi(1,1), nvctrp, 0.d0, HamSmall(1,1), orbs%norb)
else
    call dgemm('t', 'n', orbs%norb, orbs%norb, nvctrp, 1.d0, phi(1,1), nvctrp, &
               hphi(1,1), nvctrp, 0.d0, HamTemp(1,1), orbs%norb)
end if
if(nproc>1) then
    call mpi_allreduce(HamTemp(1,1), HamSmall(1,1), orbs%norb**2, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
end if

if(nproc>1) then
   iall=-product(shape(HamTemp))*kind(HamTemp)
   deallocate(HamTemp,stat=istat)
   call memocc(istat, iall, 'HamTemp', subname)
end if

end subroutine transformHam




subroutine diagonalizeHamiltonian(iproc, nproc, orbs, HamSmall, eval)
!
! Purpose:
! ========
!   Diagonalizes the Hamiltonian HamSmall and makes sure that all MPI processes give
!   the same result. This is done by requiring that the first entry of each vector
!   is positive.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc     process ID
!     nproc     number of MPI processes
!     orbs      type describing the physical orbitals psi
!   Input / Putput arguments
!     HamSmall  on input: the Hamiltonian
!               on exit: the eigenvectors
!   Output arguments
!     eval      the associated eigenvalues 
!
use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc, nproc
type(orbitals_data), intent(inout) :: orbs
real(8),dimension(orbs%norb, orbs%norb):: HamSmall
real(8),dimension(orbs%norb):: eval

! Local variables
integer:: lwork, info, istat, i, iorb, jorb
real(8),dimension(:),allocatable:: work


! Diagonalize the Hamiltonian 
lwork=-1 
allocate(work(1), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
call dsyev('v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, eval(1), work(1), lwork, info) 
lwork=work(1) 
deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
call dsyev('v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, eval(1), work(1), lwork, info) 

! Make sure that the eigenvectors are the same for all MPI processes. To do so, require that 
! the first entry of each vector is positive.
do iorb=1,orbs%norb
    if(HamSmall(1,iorb)<0.d0) then
        do jorb=1,orbs%norb
            HamSmall(jorb,iorb)=-HamSmall(jorb,iorb)
        end do
    end if
end do



end subroutine diagonalizeHamiltonian





subroutine buildWavefunction(iproc, nproc, orbs, orbsLIN, comms, commsLIN, phi, psi, HamSmall)
!
! Purpose:
! =======
!   Builds the physical orbitals psi as a linear combination of the basis functions phi. The coefficients
!   for this linear combination are obtained by diagonalizing the Hamiltonian matrix HamSmall.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc      process ID
!     nproc      total number of processes
!     orbs       type describing the physical orbitals psi
!     orbsLIN    type describing the basis functions phi
!     comms      type containing the communication parameters for the physical orbitals psi
!     commsLIN   type containing the communication parameters for the basis functions phi
!     phi        the basis functions 
!     HamSmall   the  Hamiltonian matrix
!   Output arguments:
!   -----------------
!     psi        the physical orbitals 
!

use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc, nproc
type(orbitals_data), intent(in) :: orbs
type(orbitals_data), intent(in) :: orbsLIN
type(communications_arrays), intent(in) :: comms
type(communications_arrays), intent(in) :: commsLIN
real(8),dimension(sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor,orbsLIN%norb) :: phi
real(8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor,orbs%norb) :: psi
real(8),dimension(orbsLIN%norb,orbsLIN%norb):: HamSmall

! Local variables
integer:: nvctrp


nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
call dgemm('n', 'n', nvctrp, orbs%norb, orbsLIN%norb, 1.d0, phi(1,1), nvctrp, HamSmall(1,1), &
           orbsLIN%norb, 0.d0, psi(1,1), nvctrp)



end subroutine buildWavefunction
