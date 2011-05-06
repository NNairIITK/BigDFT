subroutine allocateAndInitializeLinear(iproc, nproc, Glr, orbs, at, lin, phi, input, rxyz, occupForInguess, coeff)
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
type(locreg_descriptors),intent(in):: Glr
type(orbitals_data),intent(in):: orbs
type(atoms_data),intent(in):: at
type(linearParameters),intent(inout):: lin
type(input_variables),intent(in):: input
real(8),dimension(3,at%nat),intent(in):: rxyz
real(8),dimension(32,at%nat):: occupForInguess
real(8),dimension(:),allocatable,intent(out):: phi
real(8),dimension(:,:),allocatable,intent(out):: coeff

! Local variables
integer:: jproc, istat, iorb, jorb, ierr, iat, ityp, iall, norb_tot
integer:: norb, norbu, norbd
integer,dimension(:),allocatable:: norbsPerType, norbsPerAtom
character(len=*),parameter:: subname='allocateAndInitializeLinear'
character(len=20):: atomname
logical:: written, fileExists
real(8),dimension(:),pointer:: phiWork
real :: ttreal


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
read(99,*) lin%nItBasisFirst, lin%nItBasis
read(99,*) lin%convCrit
read(99,*) lin%DIISHistMin, lin%DIISHistMax, lin%alphaSD
read(99,*) lin%startWithSD, lin%startDIIS
read(99,*) lin%nItPrecond
read(99,*) lin%getCoeff
read(99,*) lin%nItCoeff, lin%convCritCoeff
read(99,*) lin%nItSCC, lin%alphaMix
read(99,*) lin%plotBasisFunctions
call checkLinearParameters(iproc, lin)
if(iproc==0) write(*,'(x,a)') '################################# Input parameters #################################'
if(iproc==0) write(*,'(x,a)') '>>>> General parameters.'
if(iproc==0) write(*,'(4x,a,9x,a,3x,a,3x,a,4x,a,4x,a)') '| ', ' | ', 'number of', ' | ', 'prefactor for', ' |'
if(iproc==0) write(*,'(4x,a,a,a,a,a,a,a)') '| ', 'atom type', ' | ', 'basis functions', ' | ', &
    'confinement potential', ' |'
do iat=1,at%ntypes
    read(99,*) atomname, norbsPerType(iat), lin%potentialPrefac(iat)
    if(iproc==0) write(*,'(4x,a,4x,a,a,a,a,i0,7x,a,7x,es9.3,6x,a)') '| ', trim(atomname), &
        repeat(' ', 6-len_trim(atomname)), '|', repeat(' ', 10-ceiling(log10(dble(norbsPerType(iat)+1)+1.d-10))), &
         norbsPerType(iat), '|', lin%potentialPrefac(iat), ' |'
end do
close(unit=99)
if(iproc==0) write(*,'(4x,a)') '-------------------------------------------------------'
if(iproc==0) write(*,'(4x,a)') '|  number of iterations in the   | alpha mix |'
if(iproc==0) write(*,'(4x,a)') '|     selfconsistency cycle      |           |'
if(iproc==0) write(*,'(4x,a,a,i0,16x,a,x,es9.3,x,a)') '|', repeat(' ', 16-ceiling(log10(dble(lin%nItSCC+1)+1.d-10))), &
     lin%nItSCC, '|', lin%alphaMix, '|'
if(iproc==0) write(*,'(4x,a)') '----------------------------------------------'
if(iproc==0) write(*,'(x,a)') '>>>> Parameters for the optimization of the basis functions.'
if(iproc==0) write(*,'(4x,a)') '| maximal number | convergence | iterations in  | get coef- | plot  |'
if(iproc==0) write(*,'(4x,a)') '|  of iterations |  criterion  | preconditioner | ficients  | basis |'
if(iproc==0) write(*,'(4x,a)') '|  first   else  |             |                |           |       |'
if(iproc==0) write(*,'(4x,a,a,i0,3x,a,i0,2x,a,x,es9.3,x,a,a,i0,a,a,a,l,a)') '| ', &
    repeat(' ', 5-ceiling(log10(dble(lin%nItBasisFirst+1)+1.d-10))), lin%nItBasisFirst, &
    repeat(' ', 5-ceiling(log10(dble(lin%nItBasis+1)+1.d-10))), lin%nItBasis, &
      '| ', lin%convCrit, ' | ', &
      repeat(' ', 8-ceiling(log10(dble(lin%nItPrecond+1)+1.d-10))), lin%nItPrecond, '       |   ', &
      lin%getCoeff, '    |  ', &
      lin%plotBasisFunctions, '   |'
if(iproc==0) write(*,'(4x,a)') '---------------------------------------------------------------------'
if(iproc==0) write(*,'(4x,a)') '| DIIS history | alpha SD |  start  | allow DIIS |'
if(iproc==0) write(*,'(4x,a)') '|  min   max   |          | with SD |            |'
if(iproc==0) write(*,'(4x,a,a,i0,3x,a,i0,3x,a,x,es8.2,x,a,l,a,x,es10.3,a)') '|', &
    repeat(' ', 4-ceiling(log10(dble(lin%DIISHistMin+1)+1.d-10))), lin%DIISHistMin, &
    repeat(' ', 3-ceiling(log10(dble(lin%DIISHistMax+1)+1.d-10))), lin%DIISHistMax, ' |', &
    lin%alphaSD, '|   ', lin%startWithSD, '    |', lin%startDIIS, ' |'
if(iproc==0) write(*,'(4x,a)') '--------------------------------------------------'
if(iproc==0) write(*,'(x,a)') '>>>> Parameters for the optimization of the coefficients.'
if(iproc==0) write(*,'(4x,a)') '| maximal number | convergence |'
if(iproc==0) write(*,'(4x,a)') '|  of iterations |  criterion  |'
if(iproc==0) write(*,'(4x,a,a,i0,5x,a,x,es9.3,x,a)') '| ', &
    repeat(' ', 9-ceiling(log10(dble(lin%nItCoeff+1)+1.d-10))), lin%nItCoeff, ' | ', lin%convCritCoeff, ' | '
if(iproc==0) write(*,'(4x,a)') '--------------------------------'


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
call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, input%nspin, orbs%nspinor, input%nkpt, input%kpt, input%wkpt, lin%orbs)
written=.false.
if(iproc==0) write(*,'(x,a)') '>>>> Partition of the basis functions among the processes.'
do jproc=1,nproc-1
    if(lin%orbs%norb_par(jproc)<lin%orbs%norb_par(jproc-1)) then
        !if(iproc==0) write(*,'(x,a,5(i0,a))') '#| Processes from 0 to ',jproc-1,' treat ',lin%orbs%norb_par(jproc-1), &
        !    ' orbitals, processes from ',jproc,' to ',nproc-1,' treat ',lin%orbs%norb_par(jproc),' orbitals.'
        if(iproc==0) write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',jproc-1,' treat ',&
            lin%orbs%norb_par(jproc-1), ' orbitals,', &
            repeat(' ', 3-ceiling(log10(dble(jproc)))-ceiling(log10(dble(lin%orbs%norb_par(jproc-1)+1)))), '|'
        if(iproc==0) write(*,'(4x,a,3(i0,a),a,a)')  '| processes from ',jproc,' to ',nproc-1,' treat ', &
            lin%orbs%norb_par(jproc),' orbitals.', &
            repeat(' ', 4-ceiling(log10(dble(jproc+1)))-ceiling(log10(dble(nproc)))-&
            ceiling(log10(dble(lin%orbs%norb_par(jproc)+1)))), '|'
        written=.true.
        exit
    end if
end do
if(.not.written) then
    if(iproc==0) write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',nproc-1, &
        ' treat ',lin%orbs%norbp,' orbitals. |'!, &
        !repeat(' ', 15-ceiling(log10(dble(nproc)))-ceiling(log10(dble(lin%orbs%norbp+1)))), '|'
end if
if(iproc==0) write(*,'(x,a)') '####################################################################################'


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
!call initRandomSeed(iproc, 1)
call initRandomSeed(0, 1)
!call random_number(phi)
call randomWithinCutoff(iproc, lin%orbs, Glr, at, lin, input, rxyz, phi)
!call plotOrbitals(iproc, lin%orbs, Glr, phi, at%nat, rxyz, lin%onWhichAtom, .5d0*input%hx, &
!    .5d0*input%hy, .5d0*input%hz, 1)
!!allocate(phiWork(lin%orbs%npsidim), stat=istat)
!!call memocc(istat, phiWork, 'phiWork', subname)
!!call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
!!iall=-product(shape(phiWork))*kind(phiWork)
!!deallocate(phiWork, stat=istat)
!!call memocc(istat, iall, 'phiWork', subname)


!write(*,*) 'calling createInputGuess'
!call createInputGuess(iproc, orbsLIN, Glr, input, at, rxyz, phi)

! Allocate the coefficients for the linear combinations of the  orbitals and initialize
! them at random.
! Do this only on the root, since the calculations to determine coeff are not yet parallelized.
allocate(coeff(lin%orbs%norb,orbs%norb), stat=istat)
call memocc(istat, coeff, 'coeff', subname)
call initRandomSeed(0, 1)
if(iproc==0) then
    do iorb=1,orbs%norb
       do jorb=1,lin%orbs%norb
          call random_number(ttreal)
          coeff(jorb,iorb)=real(ttreal,kind=8)
       end do
    end do
end if


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





subroutine checkLinearParameters(iproc, lin)
!
! Purpose:
! ========
!  Checks some values contained in the variable lin on errors.
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc
type(linearParameters),intent(in):: lin

! Local variables
integer:: ierr


  if(lin%DIISHistMin>lin%DIISHistMax) then
      if(iproc==0) write(*,'(a,i0,a,i0,a)') 'ERROR: DIISHistMin must not be larger than &
      & DIISHistMax, but you chose ', lin%DIISHistMin, ' and ', lin%DIISHistMax, '!'
      call mpi_barrier(mpi_comm_world, ierr)
      stop
  end if

  if(trim(lin%getCoeff)/='min' .and. trim(lin%getCoeff)/='diag') then
      if(iproc==0) write(*,'(a,a,a)') "ERROR: lin%getCoeff can have the values 'diag' or 'min', &
          & but we found '", trim(lin%getCoeff), "'!"
      call mpi_barrier(mpi_comm_world, ierr)
      stop
end if


end subroutine checkLinearParameters






subroutine deallocateLinear(iproc, lin, phi, coeff)
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
type(linearParameters),intent(inout):: lin
real(8),dimension(:),allocatable,intent(inout):: phi
real(8),dimension(:,:),allocatable,intent(inout):: coeff

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

  iall=-product(shape(coeff))*kind(coeff)
  deallocate(coeff, stat=istat)
  call memocc(istat, iall, 'coeff', subname)

end subroutine deallocateLinear





subroutine randomWithinCutoff(iproc, orbs, Glr, at, lin, input, rxyz, phi)
!
! Purpose:
! ========
!   Initializes the basis functions phi to random number within a cutoff range around
!   the atom at which they are centered.
!   The cutoff radius id given by 
!     cut=1.d0/lin%potentialPrefac(at%iatype(iiAt))
!     cut=cut**.25d0
!
use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc
type(orbitals_data), intent(inout) :: orbs
type(locreg_descriptors), intent(in) :: Glr
type(atoms_data),intent(in):: at
type(linearParameters),intent(in):: lin
type(input_variables), intent(in):: input
real(8),dimension(3,at%nat):: rxyz
real(8),dimension((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp):: phi

integer:: ix, iy, iz, ix0, iy0, iz0, iiAt, jj, iorb, i1, i2, i3, istart, ii, istat
real(8),dimension(:),allocatable:: phir
real(8):: hx, hy, hz, hxh, hyh, hzh, kx, ky, kz, tt, tt2, cut
real :: ttreal
type(workarr_sumrho) :: w
type(workarr_locham):: w_lh


allocate(phir(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i), stat=istat)
phi=0.d0

call initialize_work_arrays_sumrho(Glr,w)
call initialize_work_arrays_locham(Glr, orbs%nspinor,w_lh)

hx=input%hx
hy=input%hy
hz=input%hz
hxh=.5d0*hx
hyh=.5d0*hy
hzh=.5d0*hz

! Initialize phi to zero.
phi=0.d0

istart=0

!!call the random number as many times as the number of orbitals before
!!so that to associate unambiguously a random number to a  component-orbital pair
do ii=1,orbs%isorb*orbs%nspinor*(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i)
   call random_number(ttreal)
end do

!write(*,*) 'write, orbs%nbasisp', orbs%norbp
    orbLoop: do iorb=1,orbs%norbp
        call daub_to_isf(Glr,w,phi(istart+1),phir(1))
        iiAt=lin%onWhichAtom(iorb)
        ix0=nint(rxyz(1,iiAt)/hxh)
        iy0=nint(rxyz(2,iiAt)/hyh)
        iz0=nint(rxyz(3,iiAt)/hzh)
        cut=1.d0/lin%potentialPrefac(at%iatype(iiAt))
        cut=cut**.25d0
!cut=80000.d0

        jj=0
        !!open(unit=(iproc+1)*1000000+it*1000+iorb*10+7)
        !!open(unit=(iproc+1)*1000000+it*1000+iorb*10+8)
        !!open(unit=(iproc+1)*1000000+it*1000+iorb*10+9)
        !do i3=1,Glr%d%n3i
        do i3=-14,Glr%d%n3i-15
            !do i2=1,Glr%d%n2i
            do i2=-14,Glr%d%n2i-15
                !do i1=1,Glr%d%n1i
                do i1=-14,Glr%d%n1i-15
                  jj=jj+1
                  ! ! z component of point jj
                  ! iz=jj/(Glr%d%n2i*Glr%d%n1i)
                  ! ! Subtract the 'lower' xy layers
                  ! ii=jj-iz*(Glr%d%n2i*Glr%d%n1i)
                  ! y component of point jj
                  ! iy=ii/Glr%d%n1i
                  ! ! Subtract the 'lower' y rows
                  ! ii=ii-iy*Glr%d%n1i
                  ! ! x component
                  ! ix=ii

                   tt=hxh**2*(i1-ix0)**2 + hyh**2*(i2-iy0)**2 + hzh**2*(i3-iz0)**2
                   tt=sqrt(tt)
                   if(tt<cut) then
                      call random_number(ttreal)
                      
                      !phir(jj)=1.0d0!real(ttreal,kind=8)
                      phir(jj)=real(ttreal,kind=8)
                   else
                      !write(100+iproc,*) tt, cut
                      call random_number(ttreal)
                      !phir(jj)=tt2*.002d0*exp(-(4.d0*lin%potentialPrefac(at%iatype(iiAt))*tt))
                      phir(jj)=0.d0
                   end if
                       
                   !!!if(iy==ix0 .and. iz==iz0) write((iproc+1)*1000000+it*1000+iorb*10+7,*) ix, phir(jj)
                   !!!! Write along y-axis
                   !!!if(ix==ix0 .and. iz==iz0) write((iproc+1)*1000000+it*1000+iorb*10+8,*) iy, phir(jj)
                   !!!! Write along z-axis
                   !!!if(ix==ix0 .and. iy==iy0) write((iproc+1)*1000000+it*1000+iorb*10+9,*) iz, phir(jj)


                end do
            end do
        end do
        !!close(unit=(iproc+1)*1000000+it*1000+iorb*10+7)
        !!close(unit=(iproc+1)*1000000+it*1000+iorb*10+8)
        !!close(unit=(iproc+1)*1000000+it*1000+iorb*10+9)

        kx=orbs%kpts(1,orbs%iokpt(iorb))
        ky=orbs%kpts(2,orbs%iokpt(iorb))
        kz=orbs%kpts(3,orbs%iokpt(iorb))
        call isf_to_daub(Glr, w, phir(1), phi(istart+1))

        istart=istart+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor


    end do orbLoop

call deallocate_work_arrays_sumrho(w)
call deallocate_work_arrays_locham(Glr, w_lh)
deallocate(phir, stat=istat)


end subroutine randomWithinCutoff











subroutine plotOrbitals(iproc, orbs, Glr, phi, nat, rxyz, onWhichAtom, hxh, hyh, hzh, it)
!
! Plots the orbitals
!
use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc
type(orbitals_data), intent(inout) :: orbs
type(locreg_descriptors), intent(in) :: Glr
real(8),dimension((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp):: phi
integer:: nat
real(8),dimension(3,nat):: rxyz
integer,dimension(orbs%norbp):: onWhichAtom
real(8):: hxh, hyh, hzh
integer:: it

integer:: ix, iy, iz, ix0, iy0, iz0, iiAt, jj, iorb, i1, i2, i3, istart, ii, istat
integer:: unit1, unit2, unit3
real(8),dimension(:),allocatable:: phir
type(workarr_sumrho) :: w
character(len=10):: c1, c2, c3
character(len=50):: file1, file2, file3

allocate(phir(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i), stat=istat)

call initialize_work_arrays_sumrho(Glr,w)

istart=0

unit1=10*iproc+7
unit2=10*iproc+8
unit3=10*iproc+9

!write(*,*) 'write, orbs%nbasisp', orbs%norbp
    orbLoop: do iorb=1,orbs%norbp
        phir=0.d0
        call daub_to_isf(Glr,w,phi(istart+1),phir(1))
        iiAt=onWhichAtom(iorb)
        ix0=nint(rxyz(1,iiAt)/hxh)
        iy0=nint(rxyz(2,iiAt)/hyh)
        iz0=nint(rxyz(3,iiAt)/hzh)

        jj=0
        write(c1,'(i0)') iproc
        write(c2,'(i0)') iorb
        write(c3,'(i0)') it
        file1='orbs_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_x'
        file2='orbs_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_y'
        file3='orbs_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_z'
        open(unit=unit1, file=trim(file1))
        open(unit=unit2, file=trim(file2))
        open(unit=unit3, file=trim(file3))
        do i3=1,Glr%d%n3i
            do i2=1,Glr%d%n2i
                do i1=1,Glr%d%n1i
                   jj=jj+1
                   ! z component of point jj
                   iz=jj/(Glr%d%n2i*Glr%d%n1i)
                   ! Subtract the 'lower' xy layers
                   ii=jj-iz*(Glr%d%n2i*Glr%d%n1i)
                   ! y component of point jj
                   iy=ii/Glr%d%n1i
                   ! Subtract the 'lower' y rows
                   ii=ii-iy*Glr%d%n1i
                   ! x component
                   ix=ii
!if(phir(jj)>1.d0) write(*,'(a,3i7,es15.6)') 'WARNING: ix, iy, iz, phir(jj)', ix, iy, iz, phir(jj)
                   if(iy==ix0 .and. iz==iz0) write(unit1,*) ix, phir(jj)
                   ! Write along y-axis
                   if(ix==ix0 .and. iz==iz0) write(unit2,*) iy, phir(jj)
                   ! Write along z-axis
                   if(ix==ix0 .and. iy==iy0) write(unit3,*) iz, phir(jj)


                end do
            end do
        end do
        close(unit=unit1)
        close(unit=unit2)
        close(unit=unit3)

        istart=istart+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor

    end do orbLoop

call deallocate_work_arrays_sumrho(w)
deallocate(phir, stat=istat)


end subroutine plotOrbitals
