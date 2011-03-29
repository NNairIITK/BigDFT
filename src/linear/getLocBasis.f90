subroutine getLinearPsi(iproc, nproc, nspin, Glr, orbs, comms, at, lin, rxyz, rxyzParab, &
    nscatterarr, ngatherarr, nlpspd, proj, rhopot, GPU, input, pkernelseq, phi, psi, psit, &
    infoBasisFunctions, n3p, ebsMod)
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
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nspin, n3p
type(locreg_descriptors),intent(in):: Glr
type(orbitals_data),intent(in) :: orbs
type(communications_arrays),intent(in) :: comms
type(atoms_data),intent(in):: at
type(linearParameters),intent(in):: lin
type(input_variables),intent(in):: input
real(8),dimension(3,at%nat),intent(in):: rxyz, rxyzParab
integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
integer,dimension(0:nproc-1,2),intent(in):: ngatherarr
type(nonlocal_psp_descriptors),intent(in):: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(in):: proj
real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin),intent(inout) :: rhopot
type(GPU_pointers),intent(inout):: GPU
real(dp),dimension(:),pointer,intent(in):: pkernelseq
real(8),dimension(lin%orbs%npsidim),intent(inout):: phi
real(8),dimension(orbs%npsidim),intent(out):: psi, psit
integer,intent(out):: infoBasisFunctions
real(8),intent(out):: ebsMod

! Local variables 
integer:: istat, iall 
real(8),dimension(:),allocatable:: hphi, eval 
real(8),dimension(:,:),allocatable:: HamSmall 
real(8),dimension(:,:,:),allocatable:: matrixElements
real(8),dimension(:),pointer:: phiWork 
real(8)::epot_sum,ekin_sum,eexctX,eproj_sum, ddot, trace 
real(wp),dimension(:),pointer:: potential 
character(len=*),parameter:: subname='getLinearPsi' 

  
  allocate(hphi(lin%orbs%npsidim), stat=istat) 
  call memocc(istat, hphi, 'hphi', subname)
  allocate(phiWork(max(size(phi),size(psi))), stat=istat)
  call memocc(istat, phiWork, 'phiWork', subname)
  allocate(matrixElements(lin%orbs%norb,lin%orbs%norb,2), stat=istat)
  call memocc(istat, matrixElements, 'matrixElements', subname)
  
  
  
  call getLocalizedBasis(iproc, nproc, at, orbs, Glr, input, lin, rxyz, nspin, nlpspd, proj, &
      nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, hphi, trace, rxyzParab, &
      infoBasisFunctions)


  call HamiltonianApplicationConfinement(iproc,nproc,at,lin%orbs,lin,input%hx,input%hy,input%hz,rxyz,&
       nlpspd,proj,Glr,ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
       rhopot(1),&
       phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, rxyzParab, pkernel=pkernelseq)
  call getMatrixElements(iproc, nproc, Glr, lin, phi, hphi, matrixElements)


  if(iproc==0) write(*,'(x,a)') '----------------------------------- Determination of the orbitals in this new basis.'
  
  !allocate the potential in the full box
  call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*n3p,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,input%nspin,&
       lin%orbs%norb,lin%orbs%norbp,ngatherarr,rhopot,potential)
  
  call HamiltonianApplication(iproc,nproc,at,lin%orbs,input%hx,input%hy,input%hz,rxyz,&
       nlpspd,proj,Glr,ngatherarr,potential,&
       phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel=pkernelseq)
  if(iproc==0) write(*,'(x,a)', advance='no') 'done.'
  
  !deallocate potential
  call free_full_potential(nproc,potential,subname)
  
  
  call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
  call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, hphi, work=phiWork)
  
  allocate(HamSmall(lin%orbs%norb,lin%orbs%norb), stat=istat)
  call memocc(istat, HamSmall, 'HamSmall', subname)
  call transformHam(iproc, nproc, lin%orbs, lin%comms, phi, hphi, HamSmall)
  
  if(iproc==0) write(*,'(a)', advance='no') ' Diagonalization... '
  allocate(eval(lin%orbs%norb), stat=istat)
  call memocc(istat, eval, 'eval', subname)
  call diagonalizeHamiltonian(iproc, nproc, lin%orbs, HamSmall, eval)
  if(iproc==0) write(*,'(a)', advance='no') 'done.'

  call modifiedBSEnergy(input%nspin, orbs, lin, HamSmall(1,1), matrixElements(1,1,1), ebsMod)
  
  if(iproc==0) write(*,'(a)', advance='no') ' Linear combinations... '
  call buildWavefunction(iproc, nproc, orbs, lin%orbs, comms, lin%comms, phi, psi, HamSmall)
  
  call dcopy(orbs%npsidim, psi, 1, psit, 1)
  if(iproc==0) write(*,'(a)') 'done.'
  
  
  call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
  call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, psi, work=phiWork)
  
  
  iall=-product(shape(phiWork))*kind(phiWork)
  deallocate(phiWork, stat=istat)
  call memocc(istat, iall, 'phiWork', subname)
  
  iall=-product(shape(hphi))*kind(hphi)
  deallocate(hphi, stat=istat)
  call memocc(istat, iall, 'hphi', subname)
  
  iall=-product(shape(HamSmall))*kind(HamSmall)
  deallocate(HamSmall, stat=istat)
  call memocc(istat, iall, 'HamSmall', subname)
  
  iall=-product(shape(eval))*kind(eval)
  deallocate(eval, stat=istat)
  call memocc(istat, iall, 'eval', subname)

  iall=-product(shape(matrixElements))*kind(matrixElements)
  deallocate(matrixElements, stat=istat)
  call memocc(istat, iall, 'matrixElements', subname)

end subroutine getLinearPsi




subroutine allocateAndInitializeLinear(iproc, nproc, Glr, orbs, at, lin, phi, input, rxyz, occupForInguess)
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
use module_interfaces, exceptThisOne => allocateAndInitializeLinear
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(locreg_descriptors), intent(in):: Glr
type(orbitals_data), intent(in):: orbs
type(atoms_data), intent(in):: at
type(linearParameters):: lin
type(input_variables), intent(in):: input
real(8),dimension(3,at%nat),intent(in):: rxyz
real(8),dimension(32,at%nat):: occupForInguess
real(8),dimension(:),allocatable,intent(out):: phi

! Local variables
integer:: jproc, istat, iorb, jorb, ierr, iat, ityp, iall, norb_tot
integer:: norb, norbu, norbd
integer,dimension(:),allocatable:: norbsPerType, norbsPerAtom
character(len=*),parameter:: subname='allocateAndInitializeLinear'
character(len=20):: atomname
logical:: written, fileExists


! Allocate all local arrays.
allocate(norbsPerType(at%ntypes), stat=istat)
call memocc(istat, norbsPerType, 'norbsPerType', subname)
allocate(norbsPerAtom(at%nat), stat=istat)
call memocc(istat, norbsPerAtom, 'norbsPerAtom', subname)

! Read in all parameters related to the linear scaling version and print them.
inquire(file='input.lin', exist=fileExists)
if(.not. fileExists) then
    if(iproc==0) write(*,'(x,a)') "ERROR: the file 'input.lin' must be present for the linear &
        & scaling version!"
    call mpi_barrier(mpi_comm_world, ierr)
    stop
end if
allocate(lin%potentialPrefac(at%ntypes), stat=istat)
call memocc(istat, lin%potentialPrefac, 'lin%potentialPrefac', subname)
open(unit=99, file='input.lin')
read(99,*) lin%nItMax
read(99,*) lin%convCrit
read(99,*) lin%DIISHistMin, lin%DIISHistMax, lin%alphaSD
read(99,*) lin%startWithSD, lin%startDIIS
read(99,*) lin%nItPrecond
read(99,*) lin%plotBasisFunctions
if(lin%DIISHistMin>lin%DIISHistMax) then
    if(iproc==0) write(*,'(a,i0,a,i0,a)') 'ERROR: DIISHistMin must not be larger than &
    & DIISHistMax, but you chose ', lin%DIISHistMin, ' and ', lin%DIISHistMax, '!'
    stop
end if
if(iproc==0) write(*,'(x,a)') '################################ Input parameters. ################################'
if(iproc==0) write(*,'(x,a,9x,a,3x,a,3x,a,4x,a,4x,a)') '| ', ' | ', 'number of', ' | ', 'prefactor for', ' |'
if(iproc==0) write(*,'(x,a,a,a,a,a,a,a)') '| ', 'atom type', ' | ', 'basis functions', ' | ', &
    'confinement potential', ' |'
do iat=1,at%ntypes
    read(99,*) atomname, norbsPerType(iat), lin%potentialPrefac(iat)
    if(iproc==0) write(*,'(x,a,4x,a,a,a,a,i0,7x,a,7x,es9.3,6x,a)') '| ', trim(atomname), &
        repeat(' ', 6-len_trim(atomname)), '|', repeat(' ', 10-ceiling(log10(dble(norbsPerType(iat)+1)))), &
         norbsPerType(iat), '|', lin%potentialPrefac(iat), ' |'
end do
close(unit=99)
if(iproc==0) write(*,'(x,a)') '---------------------------------------------------------'
if(iproc==0) write(*,'(x,a)') '| maximal number | convergence | iterations in  | plot  |'
if(iproc==0) write(*,'(x,a)') '|  of iterations |  criterion  | preconditioner | basis |'
if(iproc==0) write(*,'(x,a,a,i0,5x,a,x,es9.3,x,a,a,i0,a,l,a)') '| ', &
    repeat(' ', 9-ceiling(log10(dble(lin%nItMax+1)))), lin%nItMax, ' | ', lin%convCrit, ' | ', &
      repeat(' ', 8-ceiling(log10(dble(lin%nItPrecond+1)))), lin%nItPrecond, '       |  ', &
      lin%plotBasisFunctions, '   |'
if(iproc==0) write(*,'(x,a)') '---------------------------------------------------------'
if(iproc==0) write(*,'(x,a)') '| DIIS history | alpha SD |  start  | allow DIIS |'
if(iproc==0) write(*,'(x,a)') '|  min   max   |          | with SD |            |'
if(iproc==0) write(*,'(x,a,a,i0,3x,a,i0,2x,a,x,es8.2,x,a,l,a,x,es10.3,a)') '| ', &
    repeat(' ', 4-ceiling(log10(dble(lin%DIISHistMin+1)))), lin%DIISHistMin, &
    repeat(' ', 4-ceiling(log10(dble(lin%nItMax+1)))), lin%DIISHistMax, ' |', &
    lin%alphaSD, '|   ', lin%startWithSD, '    |', lin%startDIIS, ' |'
if(iproc==0) write(*,'(x,a)') '--------------------------------------------------'


! Assign to each atom its number of basis functions and count how many basis functions 
! we have in total.
lin%orbs%norb=0
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
        !if(iproc==0) write(*,'(x,a,5(i0,a))') '#| Processes from 0 to ',jproc-1,' treat ',lin%orbs%norb_par(jproc-1), &
        !    ' orbitals, processes from ',jproc,' to ',nproc-1,' treat ',lin%orbs%norb_par(jproc),' orbitals.'
        if(iproc==0) write(*,'(x,a,2(i0,a),a,a)') '| Processes from 0 to ',jproc-1,' treat ',&
            lin%orbs%norb_par(jproc-1), ' orbitals,', &
            repeat(' ', 14-ceiling(log10(dble(jproc)))-ceiling(log10(dble(lin%orbs%norb_par(jproc-1)+1)))), '|'
        if(iproc==0) write(*,'(x,a,3(i0,a),a,a)')  '| processes from ',jproc,' to ',nproc-1,' treat ', &
            lin%orbs%norb_par(jproc),' orbitals.', &
            repeat(' ', 16-ceiling(log10(dble(jproc+1)))-ceiling(log10(dble(nproc)))-&
            ceiling(log10(dble(lin%orbs%norb_par(jproc)+1)))), '|'
        written=.true.
        exit
    end if
end do
if(.not.written) then
    if(iproc==0) write(*,'(x,a,2(i0,a),a,a)') '| Processes from 0 to ',nproc-1, &
        ' treat ',lin%orbs%norbp,' orbitals. |'!, &
        !repeat(' ', 15-ceiling(log10(dble(nproc)))-ceiling(log10(dble(lin%orbs%norbp+1)))), '|'
end if
if(iproc==0) write(*,'(x,a)') '###################################################################################'


! Decide which orbital is centered in which atom.
allocate(lin%onWhichAtom(lin%orbs%norbp), stat=istat)
call memocc(istat, lin%onWhichAtom, 'lin%onWhichAtom', subname)
call assignOrbitalsToAtoms(iproc, at%nat, lin, norbsPerAtom)

!write(*,'(a,i2,3x,200i4)') 'iproc, lin%onWhichAtom', iproc, lin%onWhichAtom


! lin%orbs%isorb is the 'first' orbital for a given MPI process.
norb_tot=0
do jproc=0,iproc-1
   norb_tot=norb_tot+lin%orbs%norb_par(jproc)
end do
!reference orbital for process
lin%orbs%isorb=norb_tot


allocate(lin%orbs%eval(lin%orbs%norb), stat=istat)
call memocc(istat, lin%orbs%eval, 'lin%orbs%eval', subname)
lin%orbs%eval=-.5d0


! Assign the parameters needed for the communication to lin%comms
call orbitals_communicators(iproc,nproc,Glr,lin%orbs,lin%comms)


! Allocate phi and initialize it at random
allocate(phi(lin%orbs%npsidim), stat=istat)
call memocc(istat, phi, 'phi', subname)
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


end subroutine allocateAndInitializeLinear



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
!     infoBasisFunctions  indicates wheter the basis functions converged to the specified limit (value is 0)
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
!  use Poisson_Solver
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
real(8),dimension(:),allocatable:: hphiold, alpha, fnrmOldArr
real(8),dimension(:,:),allocatable:: HamSmall, fnrmArr, fnrmOvrlpArr
real(8),dimension(:),pointer:: phiWork
logical:: quiet, allowDIIS, startWithSD
character(len=*),parameter:: subname='getLocalizedBasis'
character(len=1):: message
type(diis_objects):: diisLIN


  ! Allocate all local arrays
  call allocateLocalArrays()
  
  ! Initialize the DIIS parameters 
  icountSDSatur=0
  icountSwitch=0
  icountDIISFailureTot=0
  icountDIISFailureCons=0
  call initializeDIISParameters(lin%DIISHistMax)
  if(lin%startWithSD) then
      allowDIIS=.false.
      diisLIN%switchSD=.false.
      startWithSD=.true.
  else
      allowDIIS=.true.
      startWithSD=.false.
  end if
  
  
  if(iproc==0) write(*,'(x,a)') '======================== Creation of the basis functions... ========================'

  ! Assign the step size for SD iterations.
  alpha=lin%alphaSD
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
      call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
  
  
      ! Calculate the unconstrained gradient.
      if(iproc==0) then
          write(*,'(a)', advance='no') 'Hamiltonian application... '
      end if
      call HamiltonianApplicationConfinement(iproc,nproc,at,lin%orbs,lin,input%hx,input%hy,input%hz,rxyz,&
           nlpspd,proj,Glr,ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
           rhopot(1),&
           phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, rxyzParabola, pkernel=pkernelseq)
  
  
      ! Apply the orthoconstraint to the gradient. This subroutine also calculates the trace trH.
      if(iproc==0) then
          write(*,'(a)', advance='no') 'orthoconstraint... '
      end if
      call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, hphi, work=phiWork)
      call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
      call orthoconstraintNotSymmetric(iproc, nproc, lin%orbs, lin%comms, Glr%wfd, phi, hphi, trH)
  
  
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
      ! step size for the steepest descent minimization (depending on the angle 
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
      call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, hphi, work=phiWork)
  
  
      ! Precondition the gradient
      if(iproc==0) then
          write(*,'(a)') 'preconditioning. '
      end if
      gnrm=1.d3 ; gnrm_zero=1.d3
      call choosePreconditioner(iproc, nproc, lin%orbs, lin, Glr, input%hx, input%hy, input%hz, &
          lin%nItPrecond, hphi, at%nat, rxyz, at, it)
  
      ! Determine the mean step size for steepest descent iterations.
      tt=sum(alpha)
      meanAlpha=tt/dble(lin%orbs%norb)
  
      ! Write some informations to the screen.
      if(iproc==0) write(*,'(x,a,i6,2es15.7,f14.7)') 'iter, fnrm, fnrmMax, trace', it, fnrm, fnrmMax, trH
      if(iproc==0) write(1000,'(i6,2es15.7,f15.7,es12.4)') it, fnrm, fnrmMax, trH, meanAlpha
      if(fnrmMax<lin%convCrit .or. it>=lin%nItMax) then
          if(it>=lin%nItMax) then
              if(iproc==0) write(*,'(x,a,i0,a)') 'WARNING: not converged within ', it, &
                  ' iterations! Exiting loop due to limitations of iterations.'
              if(iproc==0) write(*,'(x,a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
              infoBasisFunctions=1
          else
              if(iproc==0) then
                  write(*,'(x,a,i0,a,2es15.7,f12.7)') 'converged in ', it, ' iterations.'
                  write (*,'(x,a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
              end if
              infoBasisFunctions=0
          end if
          if(iproc==0) write(*,'(x,a)') '============================= Basis functions created. ============================='
          call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
          if(lin%plotBasisFunctions) then
              call plotOrbitals(iproc, lin%orbs, Glr, phi, at%nat, rxyz, lin%onWhichAtom, .5d0*input%hx, &
                  .5d0*input%hy, .5d0*input%hz, 1)
          end if
          exit iterLoop
      end if
  
  
      call DIISorSD()
      if(iproc==0) then
          if(diisLIN%idsx>0) then
              write(*,'(x,3(a,i0))') 'DIIS informations: history length=',diisLIN%idsx, ', consecutive failures=', &
                  icountDIISFailureCons, ', total failures=', icountDIISFailureTot
          else
              if(allowDIIS) then
                  message='y'
              else
                  message='n'
              end if
              write(*,'(x,a,es9.3,a,i0,a,a)') 'steepest descent informations: mean alpha=', meanAlpha, &
              ', consecutive successes=', icountSDSatur, ', DIIS=', message
          end if
      end if
      if(.not. diisLIN%switchSD) call improveOrbitals()
  
  
  end do iterLoop


  call deallocateLocalArrays()

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
    !
    ! Purpose:
    ! ========
    !   This subroutine decides whether one should use DIIS or variable step size
    !   steepest descent to improve the orbitals. In the beginning we start with DIIS
    !   with history length lin%DIISHistMax. If DIIS becomes unstable, we switch to
    !   steepest descent. If the steepest descent iterations are successful, we switch
    !   back to DIIS, but decrease the DIIS history length by one. However the DIIS
    !   history length is limited to be larger or equal than lin%DIISHistMin.
    !


      ! Do SD anyway if the force is to large
      if(fnrmMax<lin%startDIIS .and. .not.allowDIIS) then
          allowDIIS=.true.
          if(iproc==0) write(*,'(x,a)') 'The force is small enough to allow DIIS.'
          ! This is to get the correct DIIS history 
          ! (it is chosen as max(lin%DIISHistMin,lin%DIISHistMax-icountSwitch).
          icountSwitch=-1
      else if(fnrmMax>lin%startDIIS .and. diisLIN%switchSD) then
          allowDIIS=.false.
          if(iproc==0) write(*,'(x,a)') 'The force is to large to allow DIIS.'
      end if    
      !if(.not.allowDIIS) then
          if(startWithSD .and. diisLIN%idsx>0) then
              call deallocate_diis_objects(diisLIN, subname)
              diisLIN%idsx=0
              diisLIN%switchSD=.false.
              startWithSD=.false.
          end if
          if(.not.startWithSD .and. .not.allowDIIS .and. diisLIN%idsx>0) then
              if(iproc==0) write(*,'(x,a,es10.3)') 'the force is too large, switch to SD with stepsize', alpha(1)
              call deallocate_diis_objects(diisLIN, subname)
              diisLIN%idsx=0
              diisLIN%switchSD=.false.
          end if
      !else
          ! If we swicthed to SD in the previous iteration, reset this flag.
          if(diisLIN%switchSD) diisLIN%switchSD=.false.

          ! Determine wheter the trace is decreasing (as it should) or increasing.
          ! This is done by comparing the current value with diisLIN%energy_min, which is
          ! the minimal value of the trace so far.
          if(trH<=diisLIN%energy_min) then
              ! Everything ok
              diisLIN%energy_min=trH
              diisLIN%switchSD=.false.
              itBest=it
              icountSDSatur=icountSDSatur+1
              icountDIISFailureCons=0

              ! If we are using SD (i.e. diisLIN%idsx==0) and the trace has been decreasing
              ! for at least 10 iterations, switch to DIIS. However the history length is decreased.
              if(icountSDSatur>=10 .and. diisLIN%idsx==0 .and. allowDIIS) then
                  icountSwitch=icountSwitch+1
                  idsx=max(lin%DIISHistMin,lin%DIISHistMax-icountSwitch)
                  if(iproc==0) write(*,'(x,a,i0)') 'switch to DIIS with new history length ', idsx
                  call initializeDIISParameters(idsx)
                  icountDIISFailureTot=0
                  icountDIISFailureCons=0
              end if
          else
              ! The trace is growing.
              ! Count how many times this occurs and (if we are using DIIS) switch to SD after 3 
              ! total failures or after 2 consecutive failures.
              icountDIISFailureCons=icountDIISFailureCons+1
              icountDIISFailureTot=icountDIISFailureTot+1
              icountSDSatur=0
              if((icountDIISFailureCons>=2 .or. icountDIISFailureTot>=3) .and. diisLIN%idsx>0) then
                  ! Switch back to SD. The initial step size is 1.d0.
                  alpha=lin%alphaSD
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
      !end if

    end subroutine DIISorSD


    subroutine improveOrbitals()
    !
    ! Purpose:
    ! ========
    !   This subroutine improves the basis functions by following the gradient 
    ! For DIIS 
    if (diisLIN%idsx > 0) then
       diisLIN%mids=mod(diisLIN%ids,diisLIN%idsx)+1
       diisLIN%ids=diisLIN%ids+1
    end if

    ! Follow the gradient using steepest descent.
    ! The same, but transposed
    call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, hphi, work=phiWork)
    
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
    end subroutine improveOrbitals



    subroutine allocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine allocates all local arrays.
    !

      allocate(hphiold(lin%orbs%npsidim), stat=istat)
      call memocc(istat, hphiold, 'hphiold', subname)

      allocate(alpha(lin%orbs%norb), stat=istat)
      call memocc(istat, alpha, 'alpha', subname)

      allocate(fnrmArr(lin%orbs%norb,2), stat=istat)
      call memocc(istat, fnrmArr, 'fnrmArr', subname)

      allocate(fnrmOldArr(lin%orbs%norb), stat=istat)
      call memocc(istat, fnrmOldArr, 'fnrmOldArr', subname)

      allocate(fnrmOvrlpArr(lin%orbs%norb,2), stat=istat)
      call memocc(istat, fnrmOvrlpArr, 'fnrmOvrlpArr', subname)

      allocate(phiWork(size(phi)), stat=istat)
      call memocc(istat, phiWork, 'phiWork', subname)
      
    

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

      iall=-product(shape(fnrmArr))*kind(fnrmArr)
      deallocate(fnrmArr, stat=istat)
      call memocc(istat, iall, 'fnrmArr', subname)

      iall=-product(shape(fnrmOldArr))*kind(fnrmOldArr)
      deallocate(fnrmOldArr, stat=istat)
      call memocc(istat, iall, 'fnrmOldArr', subname)

      iall=-product(shape(fnrmOvrlpArr))*kind(fnrmOvrlpArr)
      deallocate(fnrmOvrlpArr, stat=istat)
      call memocc(istat, iall, 'fnrmOvrlpArr', subname)

      iall=-product(shape(phiWork))*kind(phiWork)
      deallocate(phiWork, stat=istat)
      call memocc(istat, iall, 'phiWork', subname)
      
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
      call memocc(istat, HamTemp, 'HamTemp', subname)
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
integer:: lwork, info, istat, iall, i, iorb, jorb
real(8),dimension(:),allocatable:: work
character(len=*),parameter:: subname='diagonalizeHamiltonian'

  ! Get the optimal work array size
  lwork=-1 
  allocate(work(1), stat=istat)
  call memocc(istat, work, 'work', subname)
  call dsyev('v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, eval(1), work(1), lwork, info) 
  lwork=work(1) 

  ! Deallocate the work array ane reallocate it with the optimal size
  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  call memocc(istat, iall, 'work', subname)
  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
  call memocc(istat, work, 'work', subname)

  ! Diagonalize the Hamiltonian
  call dsyev('v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, eval(1), work(1), lwork, info) 

  ! Deallocate the work array.
  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  call memocc(istat, iall, 'work', subname)
  
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




subroutine getMatrixElements(iproc, nproc, Glr, lin, phi, hphi, matrixElements)
!
! Purpose:
! ========
!
! Calling arguments:
! ==================
!
use module_base
use module_types
use module_interfaces
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(locreg_descriptors),intent(in):: Glr
type(linearParameters),intent(in):: lin
real(8),dimension(lin%orbs%npsidim),intent(inout):: phi, hphi
real(8),dimension(lin%orbs%norb,lin%orbs%norb,2),intent(out):: matrixElements

! Local variables
integer:: istart, jstart, nvctrp, iorb, jorb, istat, iall, ierr
real(8):: ddot
real(8),dimension(:),pointer:: phiWork
character(len=*),parameter:: subname='getMatrixELements'


  allocate(phiWork(lin%orbs%npsidim), stat=istat)
  call memocc(istat, phiWork, 'phiWork', subname)


  call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
  call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, hphi, work=phiWork)

  ! Calculate <phi_i|H_j|phi_j>
  nvctrp=sum(lin%comms%nvctr_par(iproc,1:lin%orbs%nkptsp))*lin%orbs%nspinor
  jstart=1
  do jorb=1,lin%orbs%norb
      istart=1
      do iorb=1,lin%orbs%norb
          matrixElements(iorb,jorb,2)=ddot(nvctrp, phi(istart), 1, hphi(jstart), 1)
          istart=istart+nvctrp
      end do
      jstart=jstart+nvctrp
  end do
  call mpi_allreduce(matrixElements(1,1,2), matrixElements(1,1,1), lin%orbs%norb**2, &
      mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!!if(iproc==0) then
!!    write(*,*) 'matrix Elements'
!!    do iorb=1,lin%orbs%norb
!!        write(*,'(80es9.2)') (matrixElements(iorb,jorb,1), jorb=1,lin%orbs%norb)
!!    end do
!!end if


  call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
  call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, hphi, work=phiWork)

  iall=-product(shape(phiWork))*kind(phiWork)
  deallocate(phiWork)
  call memocc(istat, iall, 'phiWork', subname)


  !!! Calculate the modified band structure energy
  !!tt=0.d0
  !!do iorb=1,orbs%norb
  !!    do jorb=1,orbsLIN%norb
  !!        do korb=1,orbsLIN%norb
  !!            tt=tt+HamSmall(korb,iorb)*HamSmall(jorb,iorb)*matrixElements(korb,jorb,1)
  !!        end do
  !!    end do
  !!end do
  !!if(present(ebs_mod)) then
  !!    if(nspin==1) ebs_mod=2.d0*tt ! 2 for closed shell
  !!end if



end subroutine getMatrixElements





subroutine modifiedBSEnergy(nspin, orbs, lin, HamSmall, matrixElements, ebsMod)
!
! Purpose:
! ========
!
! Calling arguments:
! ==================
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: nspin
type(orbitals_data),intent(in) :: orbs
type(linearParameters),intent(in):: lin
real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(in):: HamSmall, matrixElements
real(8),intent(out):: ebsMod

! Local variables
integer:: iorb, jorb, korb
real(8):: tt

  ! Calculate the modified band structure energy
  tt=0.d0
  do iorb=1,orbs%norb
      do jorb=1,lin%orbs%norb
          do korb=1,lin%orbs%norb
              tt=tt+HamSmall(korb,iorb)*HamSmall(jorb,iorb)*matrixElements(korb,jorb)
          end do
      end do
  end do
  if(nspin==1) then
      ebsMod=2.d0*tt ! 2 for closed shell
  else
      ebsMod=tt
  end if



end subroutine modifiedBSEnergy





subroutine deallocateLinear(iproc, lin, phi)
!
! Purpose:
! ========
!   Deallocates all array related to the linear scaling version which have not been 
!   deallocated so far.
!
! Calling arguments:
! ==================
!
use module_base
use module_types
use module_interfaces, exceptThisOne => deallocateLinear
implicit none

! Calling arguments
integer,intent(in):: iproc
type(linearParameters):: lin
real(8),dimension(:),allocatable:: phi

! Local variables
integer:: istat, iall, iorb
character(len=*),parameter:: subname='deallocateLinear'


  iall=-product(shape(lin%potentialPrefac))*kind(lin%potentialPrefac)
  deallocate(lin%potentialPrefac, stat=istat)
  call memocc(istat, iall, 'lin%potentialPrefac', subname)
  
  iall=-product(shape(lin%onWhichAtom))*kind(lin%onWhichAtom)
  deallocate(lin%onWhichAtom, stat=istat)
  call memocc(istat, iall, 'lin%onWhichAtom', subname)
  
  call deallocate_orbs(lin%orbs,subname)

  call deallocate_comms(lin%comms,subname)
  
  iall=-product(shape(phi))*kind(phi)
  deallocate(phi, stat=istat)
  call memocc(istat, iall, 'phi', subname)

  iall=-product(shape(lin%orbs%eval))*kind(lin%orbs%eval)
  deallocate(lin%orbs%eval, stat=istat)
  call memocc(istat, iall, 'lin%orbs%eval', subname)

  if(associated(lin%wfds)) then
      !iall=-product(shape(lin%wfds))*kind(lin%wfds)
      !deallocate(lin%wfds, stat=istat)
      !call memocc(istat, iall, 'lin%wfds', subname)
      do iorb=1,lin%orbs%norbp
          call deallocate_wfd(lin%wfds(iorb,iproc), subname)
      end do
      deallocate(lin%wfds, stat=istat)
  end if

  if(associated(lin%comms%nvctr_parLIN)) then
      iall=-product(shape(lin%comms%nvctr_parLIN))*kind(lin%comms%nvctr_parLIN)
      deallocate(lin%comms%nvctr_parLIN, stat=istat)
      call memocc(istat, iall, 'lin%comms%nvctr_parLIN', subname)
  end if

  if(associated(lin%comms%ncntdLIN)) then
      iall=-product(shape(lin%comms%ncntdLIN))*kind(lin%comms%ncntdLIN)
      deallocate(lin%comms%ncntdLIN, stat=istat)
      call memocc(istat, iall, 'lin%comms%ncntdLIN', subname)
  end if

  if(associated(lin%comms%ndspldLIN)) then
      iall=-product(shape(lin%comms%ndspldLIN))*kind(lin%comms%ndspldLIN)
      deallocate(lin%comms%ndspldLIN, stat=istat)
      call memocc(istat, iall, 'lin%comms%ndspldLIN', subname)
  end if

  if(associated(lin%comms%ncnttLIN)) then
      iall=-product(shape(lin%comms%ncnttLIN))*kind(lin%comms%ncnttLIN)
      deallocate(lin%comms%ncnttLIN, stat=istat)
      call memocc(istat, iall, 'lin%comms%ncnttLIN', subname)
  end if

  if(associated(lin%comms%ndspltLIN)) then
      iall=-product(shape(lin%comms%ndspltLIN))*kind(lin%comms%ndspltLIN)
      deallocate(lin%comms%ndspltLIN, stat=istat)
      call memocc(istat, iall, 'lin%comms%ndspltLIN', subname)
  end if

  if(associated(lin%MPIComms)) then
      iall=-product(shape(lin%MPIComms))*kind(lin%MPIComms)
      deallocate(lin%MPIComms, stat=istat)
      call memocc(istat, iall, 'lin%MPIComms', subname)
  end if

  if(associated(lin%procsInComm)) then
      iall=-product(shape(lin%procsInComm))*kind(lin%procsInComm)
      deallocate(lin%procsInComm, stat=istat)
      call memocc(istat, iall, 'lin%procsInComm', subname)
  end if

  if(associated(lin%norbPerComm)) then
      iall=-product(shape(lin%norbPerComm))*kind(lin%norbPerComm)
      deallocate(lin%norbPerComm, stat=istat)
      call memocc(istat, iall, 'lin%norbPerComm', subname)
  end if

end subroutine deallocateLinear
