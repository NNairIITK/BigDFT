!> This subroutine initializes all parameters needed for the linear scaling version
!! and allocate all arrays.
subroutine allocateAndInitializeLinear(iproc, nproc, Glr, orbs, at, nlpspd, lin, phi, &
    input, rxyz, nscatterarr, coeff, lphi)
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
type(nonlocal_psp_descriptors),intent(in):: nlpspd
type(linearParameters),intent(inout):: lin
type(input_variables),intent(in):: input
real(8),dimension(3,at%nat),intent(in):: rxyz
integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
real(8),dimension(:),pointer,intent(out):: phi
real(8),dimension(:,:),pointer,intent(out):: coeff
real(8),dimension(:),pointer,intent(out):: lphi

! Local variables
integer:: norb, norbu, norbd, istat, iat, ityp, iall, ilr, iorb, tag
integer,dimension(:),allocatable:: norbsPerAtom
character(len=*),parameter:: subname='allocateAndInitializeLinear'
character(len=20),dimension(:),allocatable:: atomNames


tag=0

! Allocate all local arrays.
allocate(atomNames(at%ntypes), stat=istat)
call memocc(istat, atomNames, 'atomNames', subname)
allocate(norbsPerAtom(at%nat), stat=istat)
call memocc(istat, norbsPerAtom, 'norbsPerAtom', subname)

! Number of localization regions.
lin%nlr=at%nat
lin%lzd%nlr=at%nat
lin%lb%lzd%nlr=at%nat

! Allocate the basic arrays that are needed for reading the input parameters.
call allocateBasicArrays(at, lin)

! Read in all parameters related to the linear scaling version.
call readLinearParameters(iproc, lin, at, atomNames)

! Count the number of basis functions.
norb=0
do iat=1,at%nat
    ityp=at%iatype(iat)
    norbsPerAtom(iat)=lin%norbsPerType(ityp)
    norb=norb+norbsPerAtom(iat)
end do

! Distribute the basis functions among the processors.
norbu=norb
norbd=0
call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, input%nspin, orbs%nspinor,&
     input%nkpt, input%kpt, input%wkpt, lin%orbs)
call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, input%nspin, orbs%nspinor,&
     input%nkpt, input%kpt, input%wkpt, lin%Lorbs)
call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, input%nspin, orbs%nspinor,&
     input%nkpt, input%kpt, input%wkpt, lin%lzd%orbs)


! Do the same again, but take into acount that we may also use the derivatives of the basis functions with
! respect to x,y,z. These informations will be stored in lin%lb%orbs. If we don't use the derivtaive, then
! lin%lb%orbs will be identical to lin%orbs.
if(.not. lin%useDerivativeBasisFunctions) then
    norb=lin%orbs%norb
    norbu=norb
    norbd=0
else
    norb=4*lin%orbs%norb
    norbu=norb
    norbd=0
end if
call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, input%nspin, orbs%nspinor, input%nkpt, input%kpt, input%wkpt, lin%lb%orbs)
call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, input%nspin, orbs%nspinor, input%nkpt, input%kpt, input%wkpt, lin%lb%Lorbs)
call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, input%nspin, orbs%nspinor, input%nkpt, input%kpt, input%wkpt, lin%lb%lzd%orbs)


! Assign the parameters needed for the communication to lin%comms. Again distinguish
! between the 'normal' basis and the 'large' basis inlcuding the derivtaives.
call orbitals_communicators(iproc,nproc,Glr,lin%orbs,lin%comms)
call orbitals_communicators(iproc,nproc,Glr,lin%lb%orbs,lin%lb%comms)
call orbitals_communicators(iproc,nproc,Glr,lin%lzd%orbs,lin%lzd%comms)
call orbitals_communicators(iproc,nproc,Glr,lin%lb%lzd%orbs,lin%lb%lzd%comms)


! Write all parameters related to the linear scaling version to the screen.
call writeLinearParameters(iproc, nproc, at, lin, atomNames, lin%norbsPerType)

! Allocate (almost) all remaining arrays.
call allocateLinArrays(lin)

! Decide which orbital is centered on which atom, again for the 'normal' and
! the 'large' basis.
call assignOrbitalsToAtoms(iproc, lin%orbs, at%nat, norbsPerAtom, lin%onWhichAtom, lin%onWhichAtomAll)
if(lin%useDerivativeBasisFunctions) norbsPerAtom=4*norbsPerAtom
call assignOrbitalsToAtoms(iproc, lin%lb%orbs, at%nat, norbsPerAtom, lin%lb%onWhichAtom, lin%lb%onWhichAtomAll)
call assignOrbitalsToAtoms(iproc, lin%lb%lzd%orbs, at%nat, norbsPerAtom, lin%lb%onWhichAtom, lin%lb%onWhichAtomAll)
if(lin%useDerivativeBasisFunctions) norbsPerAtom=norbsPerAtom/4

! This is the same as above, but with orbs%inWhichLocreg instead of lin%onWhichAtom
call assignToLocreg2(iproc, at%nat, lin%lzd%nlr, input%nspin, norbsPerAtom, lin%lzd%orbs)
if(lin%useDerivativeBasisFunctions) norbsPerAtom=4*norbsPerAtom
call assignToLocreg2(iproc, at%nat, lin%lb%lzd%nlr, input%nspin, norbsPerAtom, lin%lb%lzd%orbs)
if(lin%useDerivativeBasisFunctions) norbsPerAtom=norbsPerAtom/4

! Initialize the localization regions.
call initLocregs(iproc, at%nat, rxyz, lin, input, Glr, phi, lphi)

! Maybe this could be moved to another subroutine? Or be omitted at all?
allocate(lin%orbs%eval(lin%orbs%norb), stat=istat)
call memocc(istat, lin%orbs%eval, 'lin%orbs%eval', subname)
lin%orbs%eval=-.5d0
allocate(lin%Lorbs%eval(lin%Lorbs%norb), stat=istat)
call memocc(istat, lin%Lorbs%eval, 'lin%Lorbs%eval', subname)
lin%Lorbs%eval=-.5d0
allocate(lin%lb%orbs%eval(lin%lb%orbs%norb), stat=istat)
call memocc(istat, lin%lb%orbs%eval, 'lin%lb%orbs%eval', subname)
lin%lb%orbs%eval=-.5d0
allocate(lin%lb%Lorbs%eval(lin%lb%Lorbs%norb), stat=istat)
call memocc(istat, lin%lb%Lorbs%eval, 'lin%lb%Lorbs%eval', subname)
lin%lb%Lorbs%eval=-.5d0

! Initialize the coefficients.
call initCoefficients(iproc, orbs, lin, coeff)

! Initialize the parameters for the point to point communication for the
! calculation of the charge density.
call initializeCommsSumrho2(iproc, nproc, nscatterarr, lin, tag)

! Copy Glr to lin%lzd
lin%lzd%Glr = Glr
lin%lb%lzd%Glr = Glr

! Copy nlpspd to lin%lzd
lin%lzd%Gnlpspd = nlpspd
lin%lb%lzd%Gnlpspd = nlpspd

! Set localnorb
do ilr=1,lin%lzd%nlr
    lin%lzd%Llr(ilr)%localnorb=0
    do iorb=1,lin%lzd%orbs%norbp
        if(lin%onWhichAtom(iorb)==ilr) then
            lin%lzd%Llr(ilr)%localnorb = lin%lzd%Llr(ilr)%localnorb+1
        end if
    end do
end do
! The same for the derivatives
do ilr=1,lin%lzd%nlr
    lin%lb%lzd%Llr(ilr)%localnorb=0
    do iorb=1,lin%lb%lzd%orbs%norbp
        if(lin%lb%onWhichAtom(iorb)==ilr) then
            lin%lb%lzd%Llr(ilr)%localnorb = lin%lb%lzd%Llr(ilr)%localnorb+1
        end if
    end do
    !write(*,'(a,2i4,3x,i8)') 'iproc, ilr, lin%lb%lzd%Llr(ilr)%localnorb', iproc, ilr, lin%lb%lzd%Llr(ilr)%localnorb
end do
!write(*,'(a,i4,4x,100i6)') 'iproc, lin%lb%lzd%orbs%inwhichlocreg(:)', iproc, lin%lb%lzd%orbs%inwhichlocreg(:)

! Initialize the parameters for the communication for the
! potential.
!call initializeCommunicationPotential(iproc, nproc, nscatterarr, lin)
call initializeCommunicationPotential(iproc, nproc, nscatterarr, lin%orbs, lin%lzd, lin%comgp, lin%onWhichAtomAll, tag)
call initializeCommunicationPotential(iproc, nproc, nscatterarr, lin%lb%orbs, lin%lb%lzd, lin%comgp_lb, lin%lb%onWhichAtomAll, tag)

! Initialize the parameters for the communication for the orthonormalization.
!!call initCommsOrtho(iproc, nproc, lin)
call initCommsOrtho(iproc, nproc, lin%lzd, lin%onWhichAtomAll, input, lin%op, lin%comon, tag)
call initCommsOrtho(iproc, nproc, lin%lb%lzd, lin%lb%onWhichAtomAll, input, lin%op_lb, lin%comon_lb, tag)

! Restart array for the basis functions (only needed if we use the derivative basis functions).
allocate(lin%lphiRestart(lin%lzd%orbs%npsidim), stat=istat)
call memocc(istat, lin%lphiRestart, 'lin%lphiRestart', subname)

! Stores the Hamiltonian in the basis of the localized orbitals
allocate(lin%hamold(lin%lb%orbs%norb,lin%lb%orbs%norb), stat=istat)
call memocc(istat, lin%hamold, 'lin%hamold', subname)

! Deallocate all local arrays.
iall=-product(shape(atomNames))*kind(atomNames)
deallocate(atomNames, stat=istat)
call memocc(istat, iall, 'atomNames', subname)

iall=-product(shape(norbsPerAtom))*kind(norbsPerAtom)
deallocate(norbsPerAtom, stat=istat)
call memocc(istat, iall, 'norbsPerAtom', subname)


end subroutine allocateAndInitializeLinear




subroutine readLinearParameters(iproc, lin, at, atomNames)
  use module_base
  use module_types
  implicit none
  
  integer,intent(in):: iproc
  type(linearParameters):: lin
  type(atoms_data),intent(in):: at
  character(len=20),dimension(at%ntypes):: atomNames
  !integer,dimension(at%ntypes):: norbsPerType
  
  ! Local variables
  integer:: istat, itype, ierr, iall, iat
  logical:: fileExists
  character(len=*),parameter:: subname='readLinearParameters'
  real(8),dimension(:),allocatable:: locradType
  
  allocate(locradType(at%ntypes), stat=istat)
  call memocc(istat, locradType, 'locradType', subname)
    
  ! Open the input file and read in the parameters.
  inquire(file='input.lin', exist=fileExists)
  if(.not. fileExists) then
      if(iproc==0) write(*,'(x,a)') "ERROR: the file 'input.lin' must be present for the linear &
          & scaling version!"
      call mpi_barrier(mpi_comm_world, ierr)
      stop
  end if
  open(unit=99, file='input.lin')
  read(99,*) lin%nItBasisFirst, lin%nItBasis
  read(99,*) lin%convCrit
  read(99,*) lin%DIISHistMin, lin%DIISHistMax, lin%alphaDIIS, lin%alphaSD
  read(99,*) lin%startWithSD, lin%startDIIS
  read(99,*) lin%nItPrecond
  read(99,*) lin%getCoeff
  read(99,*) lin%nItOrtho, lin%convCritOrtho
  read(99,*) lin%nItCoeff, lin%convCritCoeff
  read(99,*) lin%nItSCC, lin%alphaMix, lin%convCritMix
  read(99,*) lin%useDerivativeBasisFunctions, lin%ConfPotOrder
  read(99,*) lin%nItInguess
  read(99,*) lin%plotBasisFunctions
  read(99,*) lin%norbsPerProcIG
  call checkLinearParameters(iproc, lin)
  do itype=1,at%ntypes
      read(99,*) atomNames(itype), lin%norbsPerType(itype), lin%potentialPrefac(itype), locradType(itype)
  end do
  close(unit=99)
  
  ! Assign the localization radius to each atom.
  do iat=1,at%nat
      itype=at%iatype(iat)
      lin%locrad(iat)=locradType(itype)
  end do
  
  
  iall=-product(shape(locradType))*kind(locradType)
  deallocate(locradType, stat=istat)
  call memocc(istat, iall, 'locradType', subname)

end subroutine readLinearParameters




subroutine writeLinearParameters(iproc, nproc, at, lin, atomNames, norbsPerType)
!
! Purpose:
! ========
!   Write all parameters concerning the linear scaling version to the screen.
!
! Calling arguments:
! ==================
use module_base
use module_types
implicit none

integer,intent(in):: iproc, nproc
type(atoms_data),intent(in):: at
type(linearParameters),intent(in):: lin
character(len=20),dimension(at%ntypes),intent(in):: atomNames
integer,dimension(at%ntypes),intent(in):: norbsPerType

! Local variables
integer:: itype, jproc, len1, len2, space1, space2
logical:: written


if(iproc==0) write(*,'(x,a)') '################################# Input parameters #################################'
if(iproc==0) write(*,'(x,a)') '>>>> General parameters.'
if(iproc==0) write(*,'(4x,a)') '|           |    number of    |     prefactor for     | localization |'
if(iproc==0) write(*,'(4x,a)') '| atom type | basis functions | confinement potential |    radius    |'
do itype=1,at%ntypes
    if(iproc==0) write(*,'(4x,a,4x,a,a,a,a,i0,7x,a,7x,es9.3,6x,a,3x,f8.4,3x,a)') '| ', trim(atomNames(itype)), &
        repeat(' ', 6-len_trim(atomNames(itype))), '|', repeat(' ', 10-ceiling(log10(dble(norbsPerType(itype)+1)+1.d-10))), &
         norbsPerType(itype), '|', lin%potentialPrefac(itype), ' |', lin%locrad(itype), '|'
end do
close(unit=99)
if(iproc==0) write(*,'(4x,a)') '----------------------------------------------------------------------'
if(iproc==0) write(*,'(4x,a)') '| iterations in | alpha mix | convergence crit. | use the derivative | order of conf. |'
if(iproc==0) write(*,'(4x,a)') '|  in SC cycle  |           |    for mixing     |  basis functions   |   potential    |'
if(iproc==0) write(*,'(4x,a,a,i0,5x,a,x,es9.3,x,a,5x,es9.3,5x,a,8x,l,10x,a,7x,i1,8x,a)') '|', &
     repeat(' ', 10-ceiling(log10(dble(lin%nItSCC+1)+1.d-10))), &
     lin%nItSCC, '|', lin%alphaMix, '|', lin%convCritMix, '|', lin%useDerivativeBasisFunctions, '|', lin%confPotOrder, '|'
if(iproc==0) write(*,'(4x,a)') '---------------------------------------------------------------------------------------'
if(iproc==0) write(*,'(4x,a)') '| iterations in | orbitals per |'
if(iproc==0) write(*,'(4x,a)') '|  input guess  |   process    |'
if(iproc==0) write(*,'(4x,a,a,i0,5x,a,a,i0,6x,a)') '|', repeat(' ', 10-ceiling(log10(dble(lin%nItInguess+1)+1.d-10))), &
     lin%nItInguess, '|', repeat(' ', 8-ceiling(log10(dble(lin%norbsPerProcIG+1)+1.d-10))), lin%norbsPerProcIG, '|'
if(iproc==0) write(*,'(4x,a)') '--------------------------------'
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
if(iproc==0) write(*,'(4x,a)') '| DIIS history | alpha DIIS | alpha SD |  start  | allow DIIS | orthonormalization: |'
if(iproc==0) write(*,'(4x,a)') '|  min   max   |            |          | with SD |            | nit max   conv crit |'
if(iproc==0) write(*,'(4x,a,a,i0,3x,a,i0,3x,a,2x,es8.2,2x,a,x,es8.2,x,a,l,a,x,es10.3,a,a,i0,7x,es7.1,2x,a)') '|', &
    repeat(' ', 4-ceiling(log10(dble(lin%DIISHistMin+1)+1.d-10))), lin%DIISHistMin, &
    repeat(' ', 3-ceiling(log10(dble(lin%DIISHistMax+1)+1.d-10))), lin%DIISHistMax, ' |', &
    lin%alphaDIIS, '|', lin%alphaSD, '|   ', lin%startWithSD, '    |', lin%startDIIS, ' |', &
    repeat(' ', 5-ceiling(log10(dble(lin%nItOrtho+1)+1.d-10))), lin%nItOrtho, lin%convCritOrtho, '|'
if(iproc==0) write(*,'(4x,a)') '-------------------------------------------------------------------------------------'
if(iproc==0) write(*,'(x,a)') '>>>> Parameters for the optimization of the coefficients.'
if(iproc==0) write(*,'(4x,a)') '| maximal number | convergence |'
if(iproc==0) write(*,'(4x,a)') '|  of iterations |  criterion  |'
if(iproc==0) write(*,'(4x,a,a,i0,5x,a,x,es9.3,x,a)') '| ', &
    repeat(' ', 9-ceiling(log10(dble(lin%nItCoeff+1)+1.d-10))), lin%nItCoeff, ' | ', lin%convCritCoeff, ' | '
if(iproc==0) write(*,'(4x,a)') '--------------------------------'



written=.false.
if(iproc==0) write(*,'(x,a)') '>>>> Partition of the basis functions among the processes.'
do jproc=1,nproc-1
    if(lin%orbs%norb_par(jproc)<lin%orbs%norb_par(jproc-1)) then
        len1=1+ceiling(log10(dble(jproc-1)+1.d-5))+ceiling(log10(dble(lin%orbs%norb_par(jproc-1)+1.d-5)))
        len2=ceiling(log10(dble(jproc)+1.d-5))+ceiling(log10(dble(nproc-1)+1.d-5))+&
             ceiling(log10(dble(lin%orbs%norb_par(jproc)+1.d-5)))
        if(len1>=len2) then
            space1=1
            space2=1+len1-len2
        else
            space1=1+len2-len1
            space2=1
        end if
        if(iproc==0) write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',jproc-1,' treat ',&
            lin%orbs%norb_par(jproc-1), ' orbitals,', repeat(' ', space1), '|'
        if(iproc==0) write(*,'(4x,a,3(i0,a),a,a)')  '| processes from ',jproc,' to ',nproc-1,' treat ', &
            lin%orbs%norb_par(jproc),' orbitals.', repeat(' ', space2), '|'
        written=.true.
        exit
    end if
end do
if(.not.written) then
    if(iproc==0) write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',nproc-1, &
        ' treat ',lin%orbs%norbp,' orbitals. |'!, &
end if
if(iproc==0) write(*,'(x,a)') '-----------------------------------------------'


written=.false.
if(iproc==0) write(*,'(x,a)') '>>>> Partition of the basis functions including the derivatives among the processes.'
do jproc=1,nproc-1
    if(lin%lb%orbs%norb_par(jproc)<lin%lb%orbs%norb_par(jproc-1)) then
        len1=1+ceiling(log10(dble(jproc-1)+1.d-5))+ceiling(log10(dble(lin%lb%orbs%norb_par(jproc-1)+1.d-5)))
        len2=ceiling(log10(dble(jproc)+1.d-5))+ceiling(log10(dble(nproc-1)+1.d-5))+&
             ceiling(log10(dble(lin%lb%orbs%norb_par(jproc)+1.d-5)))
        if(len1>=len2) then
            space1=1
            space2=1+len1-len2
        else
            space1=1+len2-len1
            space2=1
        end if
        if(iproc==0) write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',jproc-1,' treat ',&
            lin%lb%orbs%norb_par(jproc-1), ' orbitals,', repeat(' ', space1), '|'
        if(iproc==0) write(*,'(4x,a,3(i0,a),a,a)')  '| processes from ',jproc,' to ',nproc-1,' treat ', &
            lin%lb%orbs%norb_par(jproc),' orbitals.', repeat(' ', space2), '|'
        written=.true.
        exit
    end if
end do
if(.not.written) then
    if(iproc==0) write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',nproc-1, &
        ' treat ',lin%lb%orbs%norbp,' orbitals. |'!, &
end if
if(iproc==0) write(*,'(x,a)') '####################################################################################'


end subroutine writeLinearParameters




subroutine assignOrbitalsToAtoms(iproc, orbs, nat, norbsPerAt, onWhichAtom, onWhichAtomAll)
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
type(orbitals_data):: orbs
integer,dimension(nat):: norbsPerAt
integer,dimension(orbs%norbp):: onWhichAtom
integer,dimension(orbs%norb):: onWhichAtomAll

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
  
  do iorb=1,orbs%norb
  
      ! Switch to the next MPI process if the numbers of orbitals for a given
      ! MPI process is reached.
      if(jorb==orbs%norb_par(jproc)) then
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
      if(iproc==jproc) onWhichAtom(jorb)=jat

      ! Global assignment, i.e. without taking in account
      ! the various MPI processes.
      onWhichAtomAll(iorb)=jat
  end do    

end subroutine assignOrbitalsToAtoms



!!subroutine assignOrbitalsToAtoms2(iproc, orbs, nat, norbsPerAt, onWhichAtom, onWhichAtomAll)
!!!
!!! Purpose:
!!! ========
!!!   Assigns the orbitals to the atoms, using the array lin%onWhichAtom.
!!!   If orbital i is centered on atom j, we have lin%onWhichAtom(i)=j.
!!!
!!! Calling arguments:
!!! ==================
!!!   Input arguments:
!!!   ----------------
!!!     iproc         process ID
!!!     nat           number of atoms
!!!     norbsPerAt    indicates how many orbitals are centered on each atom.
!!!  Input / Output arguments
!!!  ---------------------
!!!     lin           type containing parameters for the linear version
!!!
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nat
!!type(orbitals_data):: orbs
!!integer,dimension(nat):: norbsPerAt
!!integer,dimension(orbs%norbp):: onWhichAtom
!!integer,dimension(orbs%norb):: onWhichAtomAll
!!
!!! Local variables
!!integer:: jproc, iiOrb, iorb, jorb, jat
!!
!!
!!  ! Count the total 'weighted orbitals', i.e. the sum of 1/(number of orbitals for this atom).
!!  weightLocreg=0.d0
!!  do iat=1,nat
!!      weightLocreg=weightLocreg+1.d0/dble(norbsPerAt(iat))
!!  end do
!!
!!
!!  ! There are four counters:
!!  !   jproc: indicates which MPI process is handling the basis function which is being treated
!!  !   jat: counts the atom numbers
!!  !   jorb: counts the orbitals handled by a given process
!!  !   iiOrb: counts the number of orbitals for a given atoms thas has already been assigned
!!  jproc=0
!!  jat=1
!!  jorb=0
!!  iiOrb=0
!!  
!!  do iorb=1,orbs%norb
!!  
!!      ! Switch to the next MPI process if the numbers of orbitals for a given
!!      ! MPI process is reached.
!!      if(jorb==orbs%norb_par(jproc)) then
!!          jproc=jproc+1
!!          jorb=0
!!      end if
!!      
!!      ! Get the atom number. First calculate the 'weighted localization region, i.e. 1/(number of basis functions for this atom).
!!      ! Switch to the next atom if the number of basis functions for this atom is reached.
!!      if(iiOrb==norbsPerAt(jat)) then
!!          jat=jat+1
!!          iiOrb=0
!!      end if
!!      jorb=jorb+1
!!      iiOrb=iiOrb+1
!!      if(iproc==jproc) onWhichAtom(jorb)=jat
!!
!!      ! Global assignment, i.e. without taking in account
!!      ! the various MPI processes.
!!      onWhichAtomAll(iorb)=jat
!!  end do    
!!
!!end subroutine assignOrbitalsToAtoms2





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
      if(iproc==0) write(*,'(x,a,i0,a,i0,a)') 'ERROR: DIISHistMin must not be larger than &
      & DIISHistMax, but you chose ', lin%DIISHistMin, ' and ', lin%DIISHistMax, '!'
      call mpi_barrier(mpi_comm_world, ierr)
      stop
  end if

  if(trim(lin%getCoeff)/='min' .and. trim(lin%getCoeff)/='diag') then
      if(iproc==0) write(*,'(x,a,a,a)') "ERROR: lin%getCoeff can have the values 'diag' or 'min', &
          & but we found '", trim(lin%getCoeff), "'!"
      call mpi_barrier(mpi_comm_world, ierr)
      stop
  end if

  if(lin%confPotOrder/=4 .and. lin%confPotOrder/=6) then
      if(iproc==0) write(*,'(x,a,i0,a)') 'ERROR: lin%confPotOrder can have the values 4 or 6, &
          & but we found ', lin%confPotOrder, '!'
      call mpi_barrier(mpi_comm_world, ierr)
      stop
  end if

end subroutine checkLinearParameters






subroutine deallocateLinear(iproc, lin, phi, lphi, coeff)
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
real(8),dimension(:),pointer,intent(inout):: phi, lphi
real(8),dimension(:,:),pointer,intent(inout):: coeff

! Local variables
integer:: istat, iall, iorb
character(len=*),parameter:: subname='deallocateLinear'


  iall=-product(shape(lin%potentialPrefac))*kind(lin%potentialPrefac)
  deallocate(lin%potentialPrefac, stat=istat)
  call memocc(istat, iall, 'lin%potentialPrefac', subname)
  
  iall=-product(shape(lin%onWhichAtom))*kind(lin%onWhichAtom)
  deallocate(lin%onWhichAtom, stat=istat)
  call memocc(istat, iall, 'lin%onWhichAtom', subname)

  iall=-product(shape(lin%norbsPerType))*kind(lin%norbsPerType)
  deallocate(lin%norbsPerType, stat=istat)
  call memocc(istat, iall, 'lin%norbsPerType', subname)

  
  call deallocate_orbs(lin%orbs,subname)

  call deallocate_comms(lin%comms,subname)

  call deallocate_orbs(lin%lb%orbs,subname)

  call deallocate_comms(lin%lb%comms,subname)
  
  iall=-product(shape(phi))*kind(phi)
  deallocate(phi, stat=istat)
  call memocc(istat, iall, 'phi', subname)

  iall=-product(shape(lphi))*kind(lphi)
  deallocate(lphi, stat=istat)
  call memocc(istat, iall, 'lphi', subname)

  iall=-product(shape(lin%orbs%eval))*kind(lin%orbs%eval)
  deallocate(lin%orbs%eval, stat=istat)
  call memocc(istat, iall, 'lin%orbs%eval', subname)

  iall=-product(shape(lin%lb%orbs%eval))*kind(lin%lb%orbs%eval)
  deallocate(lin%lb%orbs%eval, stat=istat)
  call memocc(istat, iall, 'lin%lb%orbs%eval', subname)

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

  if(associated(lin%outofzone)) then
      iall=-product(shape(lin%outofzone))*kind(lin%outofzone)
      deallocate(lin%outofzone, stat=istat)
      call memocc(istat, iall, 'lin%outofzone', subname)
  end if

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

integer:: ix, iy, iz, ix0, iy0, iz0, iiAt, jj, iorb, i1, i2, i3, istart, ii, istat, iall
real(8),dimension(:),allocatable:: phir
real(8):: hx, hy, hz, hxh, hyh, hzh, kx, ky, kz, tt, tt2, cut
real :: ttreal
type(workarr_sumrho) :: w
type(workarr_locham):: w_lh
character(len=*),parameter:: subname='randomWithinCutoff'


allocate(phir(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i), stat=istat)
call memocc(istat, phir, 'phir', subname)
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

    orbLoop: do iorb=1,orbs%norbp
        call daub_to_isf(Glr,w,phi(istart+1),phir(1))
        iiAt=lin%onWhichAtom(iorb)
        ix0=nint(rxyz(1,iiAt)/hxh)
        iy0=nint(rxyz(2,iiAt)/hyh)
        iz0=nint(rxyz(3,iiAt)/hzh)
        cut=1.d0/lin%potentialPrefac(at%iatype(iiAt))
        cut=cut**(1.d0/dble(lin%confPotOrder))
        !cut=cut**.166666d0
!cut=80000.d0

        jj=0
        do i3=-14,Glr%d%n3i-15
            do i2=-14,Glr%d%n2i-15
                do i1=-14,Glr%d%n1i-15
                  jj=jj+1

                   tt=hxh**2*(i1-ix0)**2 + hyh**2*(i2-iy0)**2 + hzh**2*(i3-iz0)**2
                   tt=sqrt(tt)
                   if(tt<cut) then
                      call random_number(ttreal)
                      phir(jj)=real(ttreal,kind=8)
                   else
                      !call random_number(ttreal)
                      phir(jj)=0.d0
                   end if
                end do
            end do
        end do

        kx=orbs%kpts(1,orbs%iokpt(iorb))
        ky=orbs%kpts(2,orbs%iokpt(iorb))
        kz=orbs%kpts(3,orbs%iokpt(iorb))
        call isf_to_daub(Glr, w, phir(1), phi(istart+1))

        istart=istart+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor


    end do orbLoop

! Deallocate everything.
call deallocate_work_arrays_sumrho(w)
call deallocate_work_arrays_locham(Glr, w_lh)
iall=-product(shape(phir))*kind(phir)
deallocate(phir, stat=istat)
call memocc(istat, iall, 'phir', subname)


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




subroutine cutoffOutsideLocreg(iproc, nproc, Glr, at, input, lin, rxyz, phi)
! Cut off everything outside the localization region by setting it to zero.
! Then do a orthonormalization.
use module_base
use module_types
use module_interfaces
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(locreg_descriptors),intent(in) :: Glr
type(atoms_data),intent(in):: at
type(input_variables),intent(in):: input
type(linearParameters),intent(in):: lin
real(8),dimension(3,at%nat),intent(in):: rxyz
real(8),dimension(lin%orbs%npsidim),intent(inout):: phi

! Local variables
integer:: iorb, ist, i1, i2, i3, jj, iiAt, istat, iall, ierr
real(8):: tt, cut, hxh, hyh, hzh, ttIn, ttOut, ttIntot, ttOuttot
type(workarr_sumrho) :: w
real(8),dimension(:),allocatable:: phir
real(8),dimension(:),pointer:: phiWork
character(len=*),parameter:: subname='cutoffOutsideLocreg'

!write(*,*) 'in cutoffOutsideLocreg'

call initialize_work_arrays_sumrho(Glr, w)
hxh=input%hx*.5d0
hyh=input%hy*.5d0
hzh=input%hz*.5d0


allocate(phir(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i), stat=istat)
call memocc(istat, phir, 'phir', subname)

ist=1
ttIntot=0.d0
ttOuttot=0.d0
do iorb=1,lin%orbs%norbp
    ! Transform the orbitals to real space.
    phir=0.d0
    call daub_to_isf(Glr, w, phi(ist), phir(1))
    
    iiAt=lin%onWhichAtom(iorb)
    cut=lin%locrad(iiAt)
    
    jj=0
    ttIn=0.d0
    ttOut=0.d0
    do i3=-14,Glr%d%n3i-15
        do i2=-14,Glr%d%n2i-15
            do i1=-14,Glr%d%n1i-15
               jj=jj+1
               tt = (hxh*i1-rxyz(1,iiAt))**2 + (hyh*i2-rxyz(2,iiAt))**2 + (hzh*i3-rxyz(3,iiAt))**2
               tt=sqrt(tt)
               if(tt>cut) then
                  !write(*,'(a,4i7,3es20.12)') 'iorb, i1, i2, i3, tt, cut, phir(jj)', iorb, i1, i2, i3, tt, cut, phir(jj)
                  ttOut=ttOut+phir(jj)**2
                  phir(jj)=0.d0
               else
                  ttIn=ttIn+phir(jj)**2
               end if
            end do
        end do
    end do
    
    call isf_to_daub(Glr, w, phir(1), phi(ist))

    !write(*,'(a,i7,2es20.12)') 'before: iorb, ttIn, ttOut', iorb, ttIn, ttOut
    ist=ist+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)

    ttIntot = ttIntot + ttIn
    ttOuttot = ttOuttot + ttOut

end do

call mpiallred(ttIntot, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(ttOuttot, 1, mpi_sum, mpi_comm_world, ierr)
if(iproc==0) write(*,'(x,a)') 'cutting of outside localization region:'
if(iproc==0) write(*,'(3x,a,2es17.8)') 'before cut; average weights in / out:', ttIntot/dble(lin%orbs%norb), ttOuttot/dble(lin%orbs%norb)


call mpi_barrier(mpi_comm_world, ierr)
allocate(phiWork(lin%orbs%npsidim), stat=istat)
call memocc(istat, phiWork, 'phiWork', subname)
call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
call orthogonalize(iproc, nproc, lin%orbs, lin%comms, Glr%wfd, phi, input)
call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
iall=-product(shape(phiWork))*kind(phiWork)
deallocate(phiWork, stat=istat)
call memocc(istat, iall, 'phiWork', subname)

! Check
ist=1
ttIntot=0.d0
ttOuttot=0.d0
do iorb=1,lin%orbs%norbp
    ! Transform the orbitals to real space.
    phir=0.d0
    call daub_to_isf(Glr, w, phi(ist), phir(1))
    
    iiAt=lin%onWhichAtom(iorb)
    cut=lin%locrad(iiAt)
    !write(*,'(a,2i8,es10.3)') 'iorb, iiAt, cut', iorb, iiAt, cut
    
    jj=0
    ttIn=0.d0
    ttOut=0.d0
    do i3=-14,Glr%d%n3i-15
        do i2=-14,Glr%d%n2i-15
            do i1=-14,Glr%d%n1i-15
               jj=jj+1
               tt = (hxh*i1-rxyz(1,iiAt))**2 + (hyh*i2-rxyz(2,iiAt))**2 + (hzh*i3-rxyz(3,iiAt))**2
               tt=sqrt(tt)
               if(tt>cut) then
                  ttOut = ttOut + phir(jj)**2
               else
                  ttIn = ttIn + phir(jj)**2
               end if
            end do
        end do
    end do
    
    !write(*,'(a,i7,2es20.12)') 'after: iorb, ttIn, ttOut', iorb, ttIn, ttOut
    ist=ist+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)

    ttIntot = ttIntot + ttIn
    ttOuttot = ttOuttot + ttOut

end do

call mpiallred(ttIntot, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(ttOuttot, 1, mpi_sum, mpi_comm_world, ierr)
if(iproc==0) write(*,'(3x,a,2es17.8)') 'after cut; average weights in / out:', ttIntot/dble(lin%orbs%norb), ttOuttot/dble(lin%orbs%norb)

iall=-product(shape(phir))*kind(phir)
deallocate(phir, stat=istat)
call memocc(istat, iall, 'phir', subname)

call deallocate_work_arrays_sumrho(w)

end subroutine cutoffOutsideLocreg






subroutine initializeCommsSumrho2(iproc, nproc, nscatterarr, lin, tag)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
type(linearParameters),intent(inout):: lin
integer,intent(inout):: tag

! Local variables
integer:: istat, jproc, is, ie, ioverlap, i3s, i3e, ilr, iorb, is3ovrlp, n3ovrlp
character(len=*),parameter:: subname='initializeCommsSumrho'


! First count the number of overlapping orbitals for each slice.
allocate(lin%comsr%noverlaps(0:nproc-1), stat=istat)
call memocc(istat, lin%comsr%noverlaps, 'lin%comsr%noverlaps', subname)
do jproc=0,nproc-1
    is=nscatterarr(jproc,3)-14
    ie=is+nscatterarr(jproc,1)-1
    !if(iproc==0) write(*,'(a,3i8)') 'jproc, is, ie', jproc, is, ie
    ioverlap=0
    do iorb=1,lin%lb%orbs%norb
        ilr=lin%lb%onWhichAtomAll(iorb)
        i3s=2*lin%Llr(ilr)%ns3-14
        i3e=i3s+lin%Llr(ilr)%d%n3i-1
        if(i3s<=ie .and. i3e>=is) then
            ioverlap=ioverlap+1
        end if
    end do
    lin%comsr%noverlaps(jproc)=ioverlap
    !if(iproc==0) write(*,'(a,2i8)') 'jproc, lin%comsr%noverlaps(jproc)', jproc, lin%comsr%noverlaps(jproc)
end do
! Do the initialization concerning the calculation of the charge density.
allocate(lin%comsr%istarr(0:nproc-1), stat=istat)
call memocc(istat, lin%comsr%istarr, 'lin%comsr%istarr', subname)
!allocate(lin%comsr%istrarr(lin%comsr%noverlaps(iproc)), stat=istat)
allocate(lin%comsr%istrarr(0:nproc-1), stat=istat)
call memocc(istat, lin%comsr%istrarr, 'lin%comsr%istrarr', subname)
allocate(lin%comsr%overlaps(lin%comsr%noverlaps(iproc)), stat=istat)
call memocc(istat, lin%comsr%overlaps, 'lin%comsr%overlaps', subname)

allocate(lin%comsr%comarr(9,maxval(lin%comsr%noverlaps),0:nproc-1), stat=istat)
call memocc(istat, lin%comsr%comarr, 'lin%coms%commsSumrho', subname)


lin%comsr%istarr=1
lin%comsr%istrarr=1
lin%comsr%nrecvBuf=0
do jproc=0,nproc-1
    is=nscatterarr(jproc,3)-14
    ie=is+nscatterarr(jproc,1)-1
    ioverlap=0
    do iorb=1,lin%lb%orbs%norb
        ilr=lin%lb%onWhichAtomAll(iorb)
        i3s=2*lin%Llr(ilr)%ns3-14
        i3e=i3s+lin%Llr(ilr)%d%n3i-1
        if(i3s<=ie .and. i3e>=is) then
            ioverlap=ioverlap+1
            tag=tag+1
            is3ovrlp=max(is,i3s) !start of overlapping zone in z direction
            n3ovrlp=min(ie,i3e)-max(is,i3s)+1  !extent of overlapping zone in z direction
            is3ovrlp=is3ovrlp-2*lin%Llr(ilr)%ns3+15
            !call setCommunicationInformation2(jproc, iorb, is3ovrlp, n3ovrlp, lin%comsr%istrarr(jproc), tag, lin, lin%comsr%comarr(1,ioverlap,jproc))
            call setCommunicationInformation2(jproc, iorb, is3ovrlp, n3ovrlp, lin%comsr%istrarr(jproc), tag, lin%nlr, lin%Llr, &
                 lin%lb%onWhichAtomAll, lin%lb%orbs, lin%comsr%comarr(1,ioverlap,jproc))
            if(iproc==jproc) then
                !lin%comsr%sizePhibuffr = lin%comsr%sizePhibuffr + lin%Llr(ilr)%d%n1i*lin%Llr(ilr)%d%n2i*n3ovrlp
                lin%comsr%nrecvBuf = lin%comsr%nrecvBuf + lin%Llr(ilr)%d%n1i*lin%Llr(ilr)%d%n2i*n3ovrlp
                lin%comsr%overlaps(ioverlap)=iorb
                                                        !lin%Llr(ilr)%d%n1i*lin%Llr(ilr)%d%n2i*lin%Llr(ilr)%d%n3i
            end if
            lin%comsr%istrarr(jproc) = lin%comsr%istrarr(jproc) + lin%Llr(ilr)%d%n1i*lin%Llr(ilr)%d%n2i*n3ovrlp
        end if
    end do
end do


allocate(lin%comsr%communComplete(maxval(lin%comsr%noverlaps(:)),0:nproc-1), stat=istat)
call memocc(istat, lin%comsr%communComplete, 'lin%comsr%communComplete', subname)
allocate(lin%comsr%computComplete(maxval(lin%comsr%noverlaps(:)),0:nproc-1), stat=istat)
call memocc(istat, lin%comsr%computComplete, 'lin%comsr%computComplete', subname)


! Calculate the dimension of the wave function for each process.
! Do it for both the compressed ('npsidim') and for the uncompressed real space
! ('npsidimr') case.
lin%comsr%nsendBuf=0
do iorb=1,lin%lb%orbs%norbp
    ilr=lin%lb%onWhichAtom(iorb)
    lin%comsr%nsendBuf = lin%comsr%nsendBuf + lin%Llr(ilr)%d%n1i*lin%Llr(ilr)%d%n2i*lin%Llr(ilr)%d%n3i*lin%lb%orbs%nspinor
end do

allocate(lin%comsr%sendBuf(lin%comsr%nsendBuf), stat=istat)
call memocc(istat, lin%comsr%sendBuf, 'lin%comsr%sendBuf', subname)
call razero(lin%comsr%nSendBuf, lin%comsr%sendBuf)

allocate(lin%comsr%recvBuf(lin%comsr%nrecvBuf), stat=istat)
call memocc(istat, lin%comsr%recvBuf, 'lin%comsr%recvBuf', subname)
call razero(lin%comsr%nrecvBuf, lin%comsr%recvBuf)

end subroutine initializeCommsSumrho2


subroutine allocateLinArrays(lin)
use module_base
use module_types
implicit none

! Calling arguments
type(linearParameters),intent(inout):: lin

! Local variables
integer:: istat
character(len=*),parameter:: subname='allocateLinArrays'


allocate(lin%onWhichAtom(lin%orbs%norbp), stat=istat)
call memocc(istat, lin%onWhichAtom, 'lin%onWhichAtom', subname)

allocate(lin%onWhichAtomAll(lin%orbs%norb), stat=istat)
call memocc(istat, lin%onWhichAtom, 'lin%onWhichAtomAll', subname)

allocate(lin%lb%onWhichAtom(lin%lb%orbs%norbp), stat=istat)
call memocc(istat, lin%lb%onWhichAtom, 'lin%lb%onWhichAtom', subname)

allocate(lin%lb%onWhichAtomAll(lin%lb%orbs%norb), stat=istat)
call memocc(istat, lin%lb%onWhichAtom, 'lin%lb%onWhichAtomAll', subname)

allocate(lin%phiRestart(lin%orbs%npsidim), stat=istat)
call memocc(istat, lin%phiRestart, 'lin%phiRestart', subname)


end subroutine allocateLinArrays




subroutine allocateBasicArrays(at, lin)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(atoms_data),intent(in):: at
  type(linearParameters),intent(inout):: lin
  
  ! Local variables
  integer:: istat
  character(len=*),parameter:: subname='allocateBasicArrays'
  
  allocate(lin%norbsPerType(at%ntypes), stat=istat)
  call memocc(istat, lin%norbsPerType, 'lin%norbsPerType', subname)
  
  allocate(lin%potentialPrefac(at%ntypes), stat=istat)
  call memocc(istat, lin%potentialPrefac, 'lin%potentialPrefac', subname)

  allocate(lin%locrad(lin%nlr),stat=istat)
  call memocc(istat,lin%locrad,'lin%locrad',subname)
  
end subroutine allocateBasicArrays



subroutine initLocregs(iproc, nat, rxyz, lin, input, Glr, phi, lphi)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nat
real(8),dimension(3,nat),intent(in):: rxyz
type(linearParameters),intent(inout):: lin
type(input_variables),intent(in):: input
type(locreg_descriptors),intent(in):: Glr
real(8),dimension(:),pointer:: phi, lphi

! Local variables
integer:: istat, npsidim, npsidimr, iorb, ilr
character(len=*),parameter:: subname='initLocregs'

! Allocate the array of localisation regions
allocate(lin%Llr(lin%nlr),stat=istat)
allocate(lin%lzd%Llr(lin%lzd%nlr),stat=istat)
allocate(lin%lb%lzd%Llr(lin%lzd%nlr),stat=istat)
allocate(lin%outofzone(3,lin%nlr),stat=istat)
call memocc(istat,lin%outofzone,'lin%outofzone',subname)


!! Write some physical information on the Glr
!if(iproc==0) then
!    write(*,'(x,a)') '>>>>> Global localization region:'
!    write(*,'(3x,a,3i6)')'Global region n1,n2,n3: ',Glr%d%n1,Glr%d%n2,Glr%d%n3
!    write(*,'(3x,a,3i6)')'Global region n1i,n2i,n3i: ',Glr%d%n1i,Glr%d%n2i,Glr%d%n3i
!    write(*,'(3x,a,3i6)')'Global region nfl1,nfl2,nfl3: ',Glr%d%nfl1,Glr%d%nfl2,Glr%d%nfl3
!    write(*,'(3x,a,3i6)')'Global region nfu1,nfu2,nfu3: ',Glr%d%nfu1,Glr%d%nfu2,Glr%d%nfu3
!    write(*,'(3x,a,f6.2,f6.2,f6.2)')'Global dimension (x,y,z):',Glr%d%n1*input%hx,Glr%d%n2*input%hy,Glr%d%n3*input%hz
!    write(*,'(3x,a,f10.2)')'Global volume: ',Glr%d%n1*input%hx*Glr%d%n2*input%hy*Glr%d%n3*input%hz
!    write(*,'(3x,a,4i10)')'Global statistics: nseg_c, nseg_f, nvctr_c, nvctr_f',Glr%wfd%nseg_c,Glr%wfd%nseg_f,Glr%wfd%nvctr_c,Glr%wfd%nvctr_f
!    write(*,'(x,a)') '----------------------------------------------------------------------------------------------'
!end if

 call determine_locreg_periodic(iproc, lin%nlr, rxyz, lin%locrad, input%hx, input%hy, input%hz, Glr, lin%Llr)
 call determine_locreg_periodic(iproc, lin%lzd%nlr, rxyz, lin%locrad, input%hx, input%hy, input%hz, Glr, lin%lzd%Llr)
 call determine_locreg_periodic(iproc, lin%lb%lzd%nlr, rxyz, lin%locrad, input%hx, input%hy, input%hz, Glr, lin%lb%lzd%Llr)

do ilr=1,lin%nlr
    if(iproc==0) write(*,'(x,a,i0)') '>>>>>>> zone ', ilr
    if(iproc==0) write(*,'(3x,a,4i10)') 'nseg_c, nseg_f, nvctr_c, nvctr_f', lin%Llr(ilr)%wfd%nseg_c, lin%Llr(ilr)%wfd%nseg_f, lin%Llr(ilr)%wfd%nvctr_c, lin%Llr(ilr)%wfd%nvctr_f
    if(iproc==0) write(*,'(3x,a,3i8)') 'lin%Llr(ilr)%d%n1i, lin%Llr(ilr)%d%n2i, lin%Llr(ilr)%d%n3i', lin%Llr(ilr)%d%n1i, lin%Llr(ilr)%d%n2i, lin%Llr(ilr)%d%n3i
    if(iproc==0) write(*,'(a,6i8)') 'lin%Llr(ilr)%d%nfl1,lin%Llr(ilr)%d%nfu1,lin%Llr(ilr)%d%nfl2,lin%Llr(ilr)%d%nfu2,lin%Llr(ilr)%d%nfl3,lin%Llr(ilr)%d%nfu3',&
    lin%Llr(ilr)%d%nfl1,lin%Llr(ilr)%d%nfu1,lin%Llr(ilr)%d%nfl2,lin%Llr(ilr)%d%nfu2,lin%Llr(ilr)%d%nfl3,lin%Llr(ilr)%d%nfu3
    if(iproc==0) write(*,*) '---------------------------------'
    if(iproc==0) write(*,'(3x,a,4i10)') 'nseg_c, nseg_f, nvctr_c, nvctr_f', lin%lb%lzd%Llr(ilr)%wfd%nseg_c, lin%lb%lzd%Llr(ilr)%wfd%nseg_f, lin%lb%lzd%Llr(ilr)%wfd%nvctr_c, lin%lb%lzd%Llr(ilr)%wfd%nvctr_f
    if(iproc==0) write(*,'(3x,a,3i8)') 'lin%lb%lzd%Llr(ilr)%d%n1i, lin%lb%lzd%Llr(ilr)%d%n2i, lin%lb%lzd%Llr(ilr)%d%n3i', lin%lb%lzd%Llr(ilr)%d%n1i, lin%lb%lzd%Llr(ilr)%d%n2i, lin%lb%lzd%Llr(ilr)%d%n3i
    if(iproc==0) write(*,'(a,6i8)') 'lin%lb%lzd%Llr(ilr)%d%nfl1,lin%lb%lzd%Llr(ilr)%d%nfu1,lin%lb%lzd%Llr(ilr)%d%nfl2,lin%lb%lzd%Llr(ilr)%d%nfu2,lin%lb%lzd%Llr(ilr)%d%nfl3,lin%lb%lzd%Llr(ilr)%d%nfu3',&
    lin%lb%lzd%Llr(ilr)%d%nfl1,lin%lb%lzd%Llr(ilr)%d%nfu1,lin%lb%lzd%Llr(ilr)%d%nfl2,lin%lb%lzd%Llr(ilr)%d%nfu2,lin%lb%lzd%Llr(ilr)%d%nfl3,lin%lb%lzd%Llr(ilr)%d%nfu3
end do

! Calculate the dimension of the wave function for each process.
npsidim=0
do iorb=1,lin%orbs%norbp
    ilr=lin%onWhichAtom(iorb)
    npsidim = npsidim + (lin%Llr(ilr)%wfd%nvctr_c+7*lin%Llr(ilr)%wfd%nvctr_f)*lin%orbs%nspinor
end do
lin%Lorbs%npsidim=npsidim
lin%lzd%orbs%npsidim=npsidim

if(.not. lin%useDerivativeBasisFunctions) then
    lin%lb%Lorbs%npsidim=npsidim
    lin%lb%lzd%orbs%npsidim=npsidim
else
    npsidim=0
    do iorb=1,lin%lb%orbs%norbp
        ilr=lin%lb%onWhichAtom(iorb)
        npsidim = npsidim + (lin%Llr(ilr)%wfd%nvctr_c+7*lin%Llr(ilr)%wfd%nvctr_f)*lin%lb%orbs%nspinor
        !npsidimr = npsidimr + lin%Llr(ilr)%d%n1i*lin%Llr(ilr)%d%n2i*lin%Llr(ilr)%d%n3i*lin%lb%orbs%nspinor
    end do
    lin%lb%Lorbs%npsidim=npsidim
    lin%lb%lzd%orbs%npsidim=npsidim
end if



allocate(phi(lin%lb%orbs%npsidim), stat=istat)
call memocc(istat, phi, 'phi', subname)

allocate(lphi(lin%lb%Lorbs%npsidim), stat=istat)
call memocc(istat, lphi, 'lphi', subname)

allocate(lin%lphiold(lin%lb%Lorbs%npsidim), stat=istat)
call memocc(istat, lin%lphiold, 'lin%lphiold', subname)

allocate(lin%lhphiold(lin%lb%Lorbs%npsidim), stat=istat)
call memocc(istat, lin%lhphiold, 'lin%lhphiold', subname)

end subroutine initLocregs



!> Does the same as initLocregs, but has as argumenst lzd instead of lin, i.e. all quantities are
!! are assigned to lzd%Llr etc. instead of lin%Llr. Can probably completely replace initLocregs.
!subroutine initLocregs2(iproc, nat, rxyz, lzd, input, Glr, locrad, phi, lphi)
subroutine initLocregs2(iproc, nat, rxyz, lzd, input, Glr, locrad)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nat
real(8),dimension(3,nat),intent(in):: rxyz
type(linear_zone_descriptors),intent(inout):: lzd
type(input_variables),intent(in):: input
type(locreg_descriptors),intent(in):: Glr
real(8),dimension(lzd%nlr),intent(in):: locrad
!real(8),dimension(:),pointer:: phi, lphi

! Local variables
integer:: istat, npsidim, npsidimr, iorb, ilr
character(len=*),parameter:: subname='initLocregs'

! Allocate the array of localisation regions
allocate(lzd%Llr(lzd%nlr),stat=istat)
!! ATTENATION: WHAT ABAOUT OUTOFZONE??
!allocate(lin%outofzone(3,lin%nlr),stat=istat)
!call memocc(istat,lin%outofzone,'lin%outofzone',subname)


!! Write some physical information on the Glr
!if(iproc==0) then
!    write(*,'(x,a)') '>>>>> Global localization region:'
!    write(*,'(3x,a,3i6)')'Global region n1,n2,n3: ',Glr%d%n1,Glr%d%n2,Glr%d%n3
!    write(*,'(3x,a,3i6)')'Global region n1i,n2i,n3i: ',Glr%d%n1i,Glr%d%n2i,Glr%d%n3i
!    write(*,'(3x,a,3i6)')'Global region nfl1,nfl2,nfl3: ',Glr%d%nfl1,Glr%d%nfl2,Glr%d%nfl3
!    write(*,'(3x,a,3i6)')'Global region nfu1,nfu2,nfu3: ',Glr%d%nfu1,Glr%d%nfu2,Glr%d%nfu3
!    write(*,'(3x,a,f6.2,f6.2,f6.2)')'Global dimension (x,y,z):',Glr%d%n1*input%hx,Glr%d%n2*input%hy,Glr%d%n3*input%hz
!    write(*,'(3x,a,f10.2)')'Global volume: ',Glr%d%n1*input%hx*Glr%d%n2*input%hy*Glr%d%n3*input%hz
!    write(*,'(3x,a,4i10)')'Global statistics: nseg_c, nseg_f, nvctr_c, nvctr_f',Glr%wfd%nseg_c,Glr%wfd%nseg_f,Glr%wfd%nvctr_c,Glr%wfd%nvctr_f
!    write(*,'(x,a)') '----------------------------------------------------------------------------------------------'
!end if

 call determine_locreg_periodic(iproc, lzd%nlr, rxyz, locrad, input%hx, input%hy, input%hz, Glr, lzd%Llr)

!do ilr=1,lin%nlr
!    if(iproc==0) write(*,'(x,a,i0)') '>>>>>>> zone ', ilr
!    if(iproc==0) write(*,'(3x,a,4i10)') 'nseg_c, nseg_f, nvctr_c, nvctr_f', lin%Llr(ilr)%wfd%nseg_c, lin%Llr(ilr)%wfd%nseg_f, lin%Llr(ilr)%wfd%nvctr_c, lin%Llr(ilr)%wfd%nvctr_f
!    if(iproc==0) write(*,'(3x,a,3i8)') 'lin%Llr(ilr)%d%n1i, lin%Llr(ilr)%d%n2i, lin%Llr(ilr)%d%n3i', lin%Llr(ilr)%d%n1i, lin%Llr(ilr)%d%n2i, lin%Llr(ilr)%d%n3i
!    if(iproc==0) write(*,'(a,6i8)') 'lin%Llr(ilr)%d%nfl1,lin%Llr(ilr)%d%nfu1,lin%Llr(ilr)%d%nfl2,lin%Llr(ilr)%d%nfu2,lin%Llr(ilr)%d%nfl3,lin%Llr(ilr)%d%nfu3',&
!    lin%Llr(ilr)%d%nfl1,lin%Llr(ilr)%d%nfu1,lin%Llr(ilr)%d%nfl2,lin%Llr(ilr)%d%nfu2,lin%Llr(ilr)%d%nfl3,lin%Llr(ilr)%d%nfu3
!end do

! Calculate the dimension of the wave function for each process.
! Do it for both the compressed ('npsidim') and for the uncompressed real space
! ('npsidimr') case.
npsidim=0
do iorb=1,lzd%orbs%norbp
    !ilr=lin%onWhichAtom(iorb)
    ilr=lzd%orbs%inWhichLocreg(iorb)
    npsidim = npsidim + (lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f)*lzd%orbs%nspinor
end do
!! WARNING: CHECHK THIS
!lin%Lorbs%npsidim=npsidim
lzd%orbs%npsidim=npsidim

!! WARNING: CHECK THIS
!if(.not. lin%useDerivativeBasisFunctions) then
!    lin%lb%Lorbs%npsidim=npsidim
!else
!    npsidim=0
!    do iorb=1,lin%lb%orbs%norbp
!        ilr=lin%lb%onWhichAtom(iorb)
!        npsidim = npsidim + (lin%Llr(ilr)%wfd%nvctr_c+7*lin%Llr(ilr)%wfd%nvctr_f)*lin%lb%orbs%nspinor
!        npsidimr = npsidimr + lin%Llr(ilr)%d%n1i*lin%Llr(ilr)%d%n2i*lin%Llr(ilr)%d%n3i*lin%lb%orbs%nspinor
!    end do
!    lin%lb%Lorbs%npsidim=npsidim
!end if


 !! WARNING: CHECKTHIS
!allocate(phi(lin%lb%orbs%npsidim), stat=istat)
!call memocc(istat, phi, 'phi', subname)
!
!allocate(lphi(lin%lb%Lorbs%npsidim), stat=istat)
!call memocc(istat, lphi, 'lphi', subname)

end subroutine initLocregs2


!> Allocate the coefficients for the linear combinations of the  orbitals and initialize
!! them at random.
!! Do this only on the root, since the calculations to determine coeff are not yet parallelized.
subroutine initCoefficients(iproc, orbs, lin, coeff)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc
  type(orbitals_data),intent(in):: orbs
  type(linearParameters),intent(in):: lin
  real(8),dimension(:,:),pointer,intent(out):: coeff
  
  ! Local variables
  integer:: iorb, jorb, istat
  real:: ttreal
  character(len=*),parameter:: subname='initCoefficients'
  
  
  allocate(coeff(lin%lb%orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, coeff, 'coeff', subname)
  
  call initRandomSeed(0, 1)
  if(iproc==0) then
      do iorb=1,orbs%norb
         do jorb=1,lin%lb%orbs%norb
            call random_number(ttreal)
            coeff(jorb,iorb)=real(ttreal,kind=8)
         end do
      end do
  end if

end subroutine initCoefficients




!!!!!> Copies lrin to lrout.
!!!! Can be used to copy Glr to lin%lzr%Glr, can probabaly be deleted
!!!! as soon as this initialization is done somewhere else.
!!subroutine allocateAndCopyLocreg(lrin, lrout)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!type(locreg_descriptors),intent(in):: lrin
!!type(locreg_descriptors),intent(out):: lrout
!!
!!
!!lrout%geocode = lrin%geocode
!!
!!lrout%hybrid_on = lrin%hyrbid_on
!!
!!lrout%ns1 = lrin%ns1 ; lrout%ns2 = lrin%ns2 ; lrout%ns3 = lrin%ns3
!!
!!lrout%nsi1 = lrin%nsi1 ; lrout%nsi2 = lrin%nsi2 ; lrout%nsi3 = lrin%nsi3
!!
!!lrout%Localnorb = lrin%Localnorb
!!
!!call dcopy(3, lrin%outofzone(1), 1, lrout%outofzone(1), 1)
!!
!!ii=size(lrin%projflg)
!!allocate(lrout%projflg(ii))
!!call dcopy(ii, lrin%projflg(1), 1, lrout%projflg(1), 1)
!!
!!lrout%Glr = lrin%Glr
!!
!!lrout%Gnlpspd = lrin%nlpspd
!!
!!lrout%orbs = lrin%orbs
!!
!!lrout%comms = lrin%comms
!!
!!
!!
!!!> Contains the information needed for describing completely a
!!!! wavefunction localisation region
!!  type, public :: locreg_descriptors
!!     character(len=1) :: geocode
!!     logical :: hybrid_on               !< interesting for global, periodic, localisation regions
!!     integer :: ns1,ns2,ns3             !< starting point of the localisation region in global coordinates
!!     integer :: nsi1,nsi2,nsi3          !< starting point of locreg for interpolating grid
!!     integer :: Localnorb               !< number of orbitals contained in locreg
!!     integer,dimension(3) :: outofzone  !< vector of points outside of the zone outside Glr for periodic systems
!!     integer,dimension(:),pointer :: projflg    !< atoms contributing nlpsp projectors to locreg
!!     type(grid_dimensions) :: d
!!     type(wavefunctions_descriptors) :: wfd
!!     type(convolutions_bounds) :: bounds
!!  end type locreg_descriptors
!!
!!
!!
!!end subroutine allocateAndCopyLocreg

