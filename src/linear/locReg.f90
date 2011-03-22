subroutine initializeLocRegLIN(iproc, nproc, lr, orbsLIN, lin, at, input, rxyz, radii_cf)
!
! Purpose:
! ========
!   Determines the localization regions for each orbital.
!
! Calling arguments:
! ==================
!
use module_base
use module_types
use module_interfaces, exceptThisOne => initializeLocRegLIN
implicit none 

! Calling arguments
integer:: iproc, nproc
type(locreg_descriptors):: lr
type(orbitals_data),intent(in):: orbsLIN
type(linearParameters):: lin
type(atoms_data),intent(in):: at
type(input_variables),intent(in):: input
real(8),dimension(3,at%nat):: rxyz
real(8),dimension(at%ntypes,3):: radii_cf
type(communications_arrays):: commsLIN

! Local variables
integer:: iorb, iiAt, iitype, istat, iall, ierr, ii, ngroups, norbPerGroup, jprocStart, jprocEnd, i
integer:: norbtot1, norbtot2, igroup, jproc, norbpMax, lproc, uproc, wholeGroup
real(8):: radius, radiusCut, idealSplit, npsidimTemp
logical:: passedTest
logical,dimension(:,:,:,:),allocatable:: logridCut_c
logical,dimension(:,:,:,:),allocatable:: logridCut_f
integer,dimension(:),allocatable:: norbPerGroupArr, newID, newGroup, newComm
integer,dimension(:,:),allocatable:: procsInGroup
integer,dimension(:,:,:),allocatable:: tempArr
real(8),dimension(:),pointer:: psiInit, psi, psiWork
character(len=*),parameter:: subname='initializeLocRegLIN'


!! WARNING: during this subroutine lin%orbs%npsidim may be modified. Therefore copy it here
! and assign it back at the end of the subroutine.
npsidimTemp=lin%orbs%npsidim

! First check wheter we have free boundary conditions.
if(lr%geocode/='F' .and. lr%geocode/='f') then
    if(iproc==0) write(*,'(a)') 'initializeLocRegLIN only implemented for free boundary conditions!'
    call mpi_barrier(mpi_comm_world, ierr)
    stop
end if


! Copy the whole type describing the localization region. This will also automatically allocate all pointers.
! Change later?
lin%lr=lr



! logridCut is a logical array that is true for a given grid point if this point is within the
! localization radius and false otherwise. It is, of course, different for each orbital.
! Maybe later this array can be changed such that it does not cover the whole simulation box,
! but only a part of it.
allocate(logridCut_c(0:lr%d%n1,0:lr%d%n2,0:lr%d%n3,lin%orbs%norbp), stat=istat)
call memocc(istat, logridCut_c, 'logridCut_c', subname)
allocate(logridCut_f(0:lr%d%n1,0:lr%d%n2,0:lr%d%n3,lin%orbs%norbp), stat=istat)
call memocc(istat, logridCut_f, 'logridCut_f', subname)

norbpMax=maxval(lin%orbs%norb_par)
allocate(lin%wfds(norbpMax,0:nproc-1), stat=istat)
!call memocc(istat, lin%wfds, 'lin%wfds', subname)
do jproc=0,nproc-1
    do iorb=1,norbpMax
        lin%wfds(iorb,jproc)%nseg_c=0
        lin%wfds(iorb,jproc)%nseg_f=0
        lin%wfds(iorb,jproc)%nvctr_c=0
        lin%wfds(iorb,jproc)%nvctr_f=0
    end do
end do

! Now comes the loop which determines the localization region for each orbital.
! radiusCut gives the cutoff radius. More precisely: the cutoff radius is given by
! radius*radiusCut, where radius is the coarse/fine radius of the given atom type.
! IMPORTANT: The subroutines used by this part are identical to the cubic part! 
radiusCut=4.d0
do iorb=1,lin%orbs%norbp

    iiAt=lin%onWhichAtom(iorb)
    iitype=at%iatype(iiAt)
    radius=radii_cf(1,iitype)
    ! Fill logridCut. The cutoff for the localization region is given by radiusCut*radius
    call fill_logridCut(lin%lr%geocode, lin%lr%d%n1, lin%lr%d%n2, lin%lr%d%n3, 0, lin%lr%d%n1, 0, lin%lr%d%n2, 0, lin%lr%d%n3, 0, 1, &
         1, 1, rxyz(1,iiAt), radius, radiusCut, input%hx, input%hy, input%hz, logridCut_c(0,0,0,iorb))

    ! Calculate the number of segments and the number of grid points for each orbital.
    call num_segkeys(lin%lr%d%n1, lin%lr%d%n2, lin%lr%d%n3, 0, lin%lr%d%n1, 0, lin%lr%d%n2, 0, lin%lr%d%n3, logridCut_c(0,0,0,iorb), &
        lin%wfds(iorb,iproc)%nseg_c, lin%wfds(iorb,iproc)%nvctr_c)

    ! Now the same procedure for the fine radius.
    radius=radii_cf(2,iitype)
    call fill_logridCut(lin%lr%geocode, lin%lr%d%n1, lin%lr%d%n2, lin%lr%d%n3, 0, lin%lr%d%n1, 0, lin%lr%d%n2, 0, lin%lr%d%n3, 0, 1, &
         1, 1, rxyz(1,iiAt), radius, radiusCut, input%hx, input%hy, input%hz, logridCut_f(0,0,0,iorb))

    ! Calculate the number of segments and the number of grid points for each orbital.
    call num_segkeys(lin%lr%d%n1, lin%lr%d%n2, lin%lr%d%n3, 0, lin%lr%d%n1, 0, lin%lr%d%n2, 0, lin%lr%d%n3, logridCut_f(0,0,0,iorb), &
        lin%wfds(iorb,iproc)%nseg_f, lin%wfds(iorb,iproc)%nvctr_f)


    ! Now fill the descriptors.
    call allocate_wfd(lin%wfds(iorb,iproc), subname)
    ! First the coarse part
    call segkeys(lin%lr%d%n1, lin%lr%d%n2, lin%lr%d%n3, 0, lin%lr%d%n1, 0, lin%lr%d%n2, 0, lin%lr%d%n3, logridCut_c(0,0,0,iorb), &
        lin%wfds(iorb,iproc)%nseg_c, lin%wfds(iorb,iproc)%keyg(1,1), lin%wfds(iorb,iproc)%keyv(1))
    ! And then the fine part
    ii=lin%wfds(iorb,iproc)%nseg_c+1
    call segkeys(lin%lr%d%n1, lin%lr%d%n2, lin%lr%d%n3, 0, lin%lr%d%n1, 0, lin%lr%d%n2, 0, lin%lr%d%n3, logridCut_f(0,0,0,iorb), &
        lin%wfds(iorb,iproc)%nseg_f, lin%wfds(iorb,iproc)%keyg(1,ii), lin%wfds(iorb,iproc)%keyv(ii))

end do

! Now each orbital knows only its own localization region. However we want each orbital to
! know the number of grid points of all the other orbitals. This is done in the following.
allocate(tempArr(1:norbpMax,0:nproc-1,2), stat=istat)
call memocc(istat, tempArr, 'tempArr', subname)
tempArr=0
! First copy the number of coarse grid points (of all orbital on iproc) to a auxiliary array and sum it up
! among all processes. Since each process filled only its own part, these arrays are then identical
! on all processes.
do iorb=1,lin%orbs%norbp
    tempArr(iorb,iproc,2)=lin%wfds(iorb,iproc)%nvctr_c
end do
call mpi_allreduce(tempArr(1,0,2), tempArr(1,0,1), norbpMax*nproc, mpi_integer, mpi_sum, mpi_comm_world, ierr)
! Now assign the number of coarse grid points (of all orbitals on all processes) to lin%wfds(iorb,jproc)%nvctr_c.
do jproc=0,nproc-1
    do iorb=1,lin%orbs%norb_par(jproc)
        lin%wfds(iorb,jproc)%nvctr_c=tempArr(iorb,jproc,1)
        tempArr(iorb,jproc,1)=0
    end do
end do
! Now do the same with the fine quantities.
do iorb=1,lin%orbs%norbp
    tempArr(iorb,iproc,2)=lin%wfds(iorb,iproc)%nvctr_f
end do
call mpi_allreduce(tempArr(1,0,2), tempArr(1,0,1), norbpMax*nproc, mpi_integer, mpi_sum, mpi_comm_world, ierr)
do jproc=0,nproc-1
    do iorb=1,lin%orbs%norb_par(jproc)
        lin%wfds(iorb,jproc)%nvctr_f=tempArr(iorb,jproc,1)
    end do
end do


! Now divide the system in parts, i.e. create new MPI communicators that include only
! a part of the orbitals. For instance, group 1 contains the MPI processes 0 to 10 and
! group 2 contains the MPI processes 11 to 20. We specify how many orbitals a given group
! should contain and then assign the best matching number of MPI processes to this group.

! norbPerGroup gives the ideal number of orbitals per group.
! If you don't want to deal with these groups and have only one group (as usual), simply put
! norbPerGroup equal to lin%orbs%norb
!norbPerGroup=30
norbPerGroup=lin%orbs%norb

! ngroups is the number of groups that we will have.
ngroups=nint(dble(lin%orbs%norb)/dble(norbPerGroup))
if(ngroups>nproc) then
    if(iproc==0) write(*,'(a,i0,a,i0,a)') 'WARNING: change ngroups from ', ngroups, ' to ', nproc,'!'
    ngroups=nproc
end if

! idealSplit gives the number of orbitals that would be assigned to each group
! in the ideal case.
idealSplit=dble(lin%orbs%norb)/dble(ngroups)

! Now distribute the orbitals to the groups. Do not split MPI processes, i.e.
! all orbitals for one process will remain with this proccess.
! The procedure is as follows: We weant to assign idealSplit orbitals to a group.
! To do so, we iterate through the MPI processes and sum up the number of orbitals.
! If we are at process k of this iteration, then norbtot1 gives the sum of the orbitals
! up to process k, and norbtot2 the orbitals up to process k+1. If norbtot2 is closer
! to idealSplit than norbtot1, we continue the iteration, otherwise we split the groups
! at process k.
!allocate(norbPerGroupArr(ngroups), stat=istat)
!allocate(procsInGroup(2,ngroups))
lin%ncomms=ngroups
allocate(lin%procsInComm(2,lin%ncomms), stat=istat)
call memocc(istat, lin%procsInComm, 'lin%procsInComm', subname)
allocate(lin%norbPerComm(lin%ncomms), stat=istat)
call memocc(istat, lin%norbPerComm, 'lin%norbPerComm', subname)
norbtot1=0
norbtot2=0
igroup=1
jprocStart=0
jprocEnd=0
do jproc=0,nproc-1
    if(igroup==ngroups) then
        ! This last group has to take all the rest
        do ii=jproc,nproc-1
            norbtot1=norbtot1+lin%orbs%norb_par(jproc)
            jprocEnd=jprocEnd+1
        end do
    else
        norbtot1=norbtot1+orbsLIN%norb_par(jproc)
        if(jproc<nproc-1) norbtot2=norbtot1+lin%orbs%norb_par(jproc+1)
        jprocEnd=jprocEnd+1
    end if
    if(abs(dble(norbtot1)-idealSplit)<abs(dble(norbtot2-idealSplit)) .or. igroup==ngroups) then
        ! Here is the split between two groups
        lin%norbPerComm(igroup)=norbtot1
        lin%procsInComm(1,igroup)=jprocStart
        lin%procsInComm(2,igroup)=jprocEnd-1
        norbtot1=0
        norbtot2=0
        jprocStart=jproc+1
        if(igroup==ngroups) exit
        igroup=igroup+1
    end if
end do



! Now create the new MPI communicators.
! These communicators will be contained in the array newComm. If you want to
! use MPI processes only for the processes in group igroup, you can use the
! ordinary MPI routines just with newComm(igroup) instead of mpi_comm_world.
allocate(newID(0:nproc), stat=istat)
call memocc(istat, newID, 'newID', subname)
allocate(newGroup(1:ngroups))
call memocc(istat, newGroup, 'newGroup', subname)
allocate(lin%MPIcomms(1:ngroups), stat=istat)
call memocc(istat, lin%MPIComms, 'lin%MPIComms', subname)
do jproc=0,nproc-1
    newID(jproc)=jproc
end do
call mpi_comm_group(mpi_comm_world, wholeGroup, ierr)
do igroup=1,ngroups
    call mpi_group_incl(wholeGroup, newID(lin%procsInComm(2,igroup))-newID(lin%procsInComm(1,igroup))+1,&
        newID(lin%procsInComm(1,igroup)), newGroup(igroup), ierr)
    call mpi_comm_create(mpi_comm_world, newGroup(igroup), lin%MPIComms(igroup), ierr)
end do



! Now create the parameters for the transposition.
! lproc and uproc give the first and last process ID of the processes
! in the communicator igroup.
do igroup=1,ngroups
    lproc=lin%procsInComm(1,igroup)
    uproc=lin%procsInComm(2,igroup)
    if(iproc>=lproc .and. iproc<=uproc) then
        call orbitalsCommunicatorsWithGroups(iproc, lproc, uproc, lin, lr, lin%orbs, lin%comms, lin%MPIComms(igroup), lin%norbPerComm(igroup))
    end if
end do


!!! Write out the parameters for the transposition.
!!do igroup=1,ngroups
!!    lproc=lin%procsInComm(1,igroup)
!!    uproc=lin%procsInComm(2,igroup)
!!    do jproc=lproc,uproc
!!        if(iproc>=lproc .and. iproc<=uproc) write(*,'(a,3i5,4i12)') 'iproc, igroup, jproc, lin%comms%ncntdLIN(jproc), &
!!            & lin%comms%ndspldLIN(jproc), lin%comms%ncnttLIN(jproc), lin%comms%ndspltLIN(jproc)', &
!!            iproc, igroup, jproc, lin%comms%ncntdLIN(jproc), lin%comms%ndspldLIN(jproc), &
!!            lin%comms%ncnttLIN(jproc), lin%comms%ndspltLIN(jproc)
!!    end do
!!end do


! Test the transposition. To do so, transpose and retranspose an array. If the transposition is
! correct, this should give the same result
allocate(psi(lin%orbs%npsidim), stat=istat)
allocate(psiInit(lin%orbs%npsidim), stat=istat)
allocate(psiWork(lin%orbs%npsidim), stat=istat)
call random_number(psi)
psiInit=psi
!ii=0
!do iorb=1,lin%orbs%norbp
!    do i=1,lin%wfds(iorb,iproc)%nvctr_c+7*lin%wfds(iorb,iproc)%nvctr_f
!        ii=ii+1
!        write(100+iproc,*) ii,psi(ii)
!    end do
!end do
do igroup=1,ngroups
    lproc=lin%procsInComm(1,igroup)
    uproc=lin%procsInComm(2,igroup)
    !write(*,'(a,4i8)') 'igroup, lproc, uproc, lin%norbPerComm(igroup)', igroup, lproc, uproc, lin%norbPerComm(igroup)
    if(iproc>=lproc .and. iproc<=uproc) then
        call transpose_vLIN(iproc, lproc, uproc, lin%norbPerComm(igroup), lin%orbs, lin%comms, psi, lr, lin%MPIComms(igroup), work=psiWork)
        call untranspose_vLIN(iproc, lproc, uproc, lin%norbPerComm(igroup), lin%orbs, lin%comms, psi, lr, lin%MPIComms(igroup), work=psiWork)
    end if
end do

passedTest=.true.
ii=0
do iorb=1,lin%orbs%norbp
    do i=1,lin%wfds(iorb,iproc)%nvctr_c+7*lin%wfds(iorb,iproc)%nvctr_f
        ii=ii+1
        !write(300+iproc,*) ii,psi(ii)
        if(psi(ii)/=psiInit(ii)) then
            passedTest=.false.
        end if
    end do
end do
if(passedTest) then
    write(*,'(x,a,i0,a)') 'transposition test passed on process ', iproc, '.'
else
    write(*,'(x,a,i0,a)') 'transposition test failes on process ', iproc, '.'
    write(*,'(x,a)') 'The program will stop now!'
    stop
end if
nullify(psi)
nullify(psiWork)



iall=-product(shape(logridCut_c))*kind(logridCut_c)
deallocate(logridCut_c, stat=istat)
call memocc(istat, iall, 'logridCut_c', subname)
iall=-product(shape(logridCut_f))*kind(logridCut_f)
deallocate(logridCut_f, stat=istat)
call memocc(istat, iall, 'logridCut_f', subname)
iall=-product(shape(tempArr))*kind(tempArr)
deallocate(tempArr, stat=istat)
call memocc(istat, iall, 'tempArr', subname)
iall=-product(shape(newID))*kind(newID)
deallocate(newID, stat=istat)
call memocc(istat, iall, 'newID', subname)
iall=-product(shape(newGroup))*kind(newGroup)
deallocate(newGroup)
call memocc(istat, iall, 'newGroup', subname)


!! WARNING: assign back the original value of lin%orbs%npsidim
lin%orbs%npsidim=npsidimTemp

end subroutine initializeLocRegLIN




!!****f* BigDFT/fill_logrid
!! FUNCTION
!!   set up an array logrid(i1,i2,i3) that specifies whether the grid point
!!   i1,i2,i3 is the center of a scaling function/wavelet
!!
!! SOURCE
!!
subroutine fill_logridCut(geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,  &
     ntypes,iatype,rxyz,radii,rmult,hx,hy,hz,logrid)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,ntypes
  real(gp), intent(in) :: rmult,hx,hy,hz
  integer, dimension(nat), intent(in) :: iatype
  real(gp), dimension(ntypes), intent(in) :: radii
  real(gp), dimension(3,nat), intent(in) :: rxyz
  logical, dimension(0:n1,0:n2,0:n3), intent(out) :: logrid
  !local variables
  real(kind=8), parameter :: eps_mach=1.d-12
  integer :: i1,i2,i3,iat,ml1,ml2,ml3,mu1,mu2,mu3,j1,j2,j3
  real(gp) :: dx,dy2,dz2,rad


  !some checks
  if (geocode /='F') then
     !the nbuf value makes sense only in the case of free BC
     if (nbuf /=0) then
        write(*,'(1x,a)')'ERROR: a nonzero value of nbuf is allowed only for Free BC (tails)'
        stop
     end if
     !the grid spacings must be the same
     if (hx/= hy .or. hy /=hz .or. hx/=hz) then
!        write(*,'(1x,a)')'ERROR: For Free BC the grid spacings must be the same'
     end if
  end if

  if (geocode == 'F') then
     do i3=nl3,nu3
        do i2=nl2,nu2
           do i1=nl1,nu1
              logrid(i1,i2,i3)=.false.
           enddo
        enddo
     enddo
  else !
     do i3=0,n3
        do i2=0,n2
           do i1=0,n1
              logrid(i1,i2,i3)=.false.
           enddo
        enddo
     enddo
  end if

  do iat=1,nat
     rad=radii(iatype(iat))*rmult+real(nbuf,gp)*hx
     if (rad /= 0.0_gp) then
        ml1=max(ceiling((rxyz(1,iat)-rad)/hx - eps_mach), nl1)
        ml2=max(ceiling((rxyz(2,iat)-rad)/hy - eps_mach), nl2)
        ml3=max(ceiling((rxyz(3,iat)-rad)/hz - eps_mach), nl3)
        mu1=min(floor((rxyz(1,iat)+rad)/hx + eps_mach), nu1)
        mu2=min(floor((rxyz(2,iat)+rad)/hy + eps_mach), nu2)
        mu3=min(floor((rxyz(3,iat)+rad)/hz + eps_mach), nu3)
        !for Free BC, there must be no incoherences with the previously calculated delimiters
        if (geocode == 'F') then
           if (ml1 < nl1) stop 'ml1 < nl1'
           if (ml2 < nl2) stop 'ml2 < nl2'
           if (ml3 < nl3) stop 'ml3 < nl3'

           if (mu1 > nu1) stop 'mu1 > nu1'
           if (mu2 > nu2) stop 'mu2 > nu2'
           if (mu3 > nu3) stop 'mu3 > nu3'
        end if
        !what follows works always provided the check before
!$omp parallel default(shared) private(i3,dz2,j3,i2,dy2,j2,i1,j1,dx)
!$omp do
        do i3=max(ml3,-n3/2-1),min(mu3,n3+n3/2+1)
           dz2=(real(i3,gp)*hz-rxyz(3,iat))**2
           j3=modulo(i3,n3+1)
           do i2=max(ml2,-n2/2-1),min(mu2,n2+n2/2+1)
              dy2=(real(i2,gp)*hy-rxyz(2,iat))**2
              j2=modulo(i2,n2+1)
              do i1=max(ml1,-n1/2-1),min(mu1,n1+n1/2+1)
                 j1=modulo(i1,n1+1)
                 dx=real(i1,gp)*hx-rxyz(1,iat)
                 if (dx**2+(dy2+dz2) <= rad**2) then
                    logrid(j1,j2,j3)=.true.
                 endif
              enddo
           enddo
        enddo
!$omp enddo
!$omp end parallel
  end if
  enddo

END SUBROUTINE fill_logridCut
!!***









!!****f* BigDFT/orbitals_communicatorsLIN_group
!! FUNCTION
!!   Partition the orbitals between processors to ensure load balancing
!!   the criterion will depend on GPU computation
!!   and/or on the sizes of the different localisation region
!! DESCRIPTION
!!   Calculate the number of elements to be sent to each process
!!   and the array of displacements
!!   Cubic strategy: 
!!      - the components are equally distributed among the wavefunctions
!!      - each processor has all the orbitals in transposed form
!!      - each wavefunction is equally distributed in its transposed form
!!      - this holds for each k-point, which regroups different processors
!!
!! SOURCE
!!
!subroutine orbitals_communicatorsLIN_group(iproc, lproc, uproc, lin, lr, orbs, comms, newComm, norbPerGroup)
subroutine orbitalsCommunicatorsWithGroups(iproc, lproc, uproc, lin, lr, orbs, comms, newComm, norbPerGroup)
  !calculate the number of elements to be sent to each process
  !and the array of displacements
  !cubic strategy: -the components are equally distributed among the wavefunctions
  !                -each processor has all the orbitals in transposed form
  !                -each wavefunction is equally distributed in its transposed form
  !                -this holds for each k-point, which regroups different processors
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc, lproc, uproc
  type(linearParameters):: lin
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(inout) :: orbs
  type(communications_arrays), intent(out) :: comms
  integer:: newComm
  integer:: norbPerGroup
  !local variables
  logical:: yesorb, yescomp
  integer:: jproc, nvctr_tot, ikpts, iorb, iorbp, jorb, norb_tot, ikpt, istat, iall, ii, outproc, nproc
  integer:: ncomp_res, nkptsp, ierr
  integer, dimension(:), allocatable:: mykpts
  logical, dimension(:), allocatable:: GPU_for_comp
  integer, dimension(:,:),allocatable:: nvctr_par, norb_par !for all the components and orbitals (with k-pts)
  integer, dimension(:,:,:,:),allocatable:: nvctr_parLIN !for all the components and orbitals (with k-pts)
  character(len=*),parameter:: subname='orbitalsCommunicatorsWithGroups'
  
    nproc=uproc-lproc+1



  !check of allocation of important arrays
  if (.not. associated(orbs%norb_par)) then
     write(*,*)'ERROR: norb_par array not allocated'
     stop
  end if
  

  ! Allocate the local arrays.
  allocate(nvctr_par(lproc:uproc,0:orbs%nkpts+ndebug),stat=istat)
  call memocc(istat,nvctr_par,'nvctr_par',subname)
  allocate(norb_par(lproc:uproc,0:orbs%nkpts+ndebug),stat=istat)
  call memocc(istat,norb_par,'norb_par',subname)
  allocate(mykpts(orbs%nkpts+ndebug),stat=istat)
  call memocc(istat,mykpts,'mykpts',subname)

  ! Initialise the arrays
  do ikpts=0,orbs%nkpts
     do jproc=lproc,uproc
        nvctr_par(jproc,ikpts)=0 
        norb_par(jproc,ikpts)=0 
     end do
  end do


  ! Distribute the orbitals among the processes, taking into acount k-points.
  ! norb_par(jproc,ikpts)=ii means that process jproc holds ii orbitals
  ! of k-point ikpts
  jorb=1
  ikpts=1
  do jproc=lproc,uproc
     do iorbp=1,orbs%norb_par(jproc)
        norb_par(jproc,ikpts)=norb_par(jproc,ikpts)+1
        if (mod(jorb,orbs%norb)==0) then
           ikpts=ikpts+1
        end if
        jorb=jorb+1
     end do
  end do

  allocate(nvctr_parLIN(maxval(norb_par(:,1)),lproc:uproc,lproc:uproc,0:orbs%nkpts+ndebug),stat=istat)
  call memocc(istat,nvctr_parLIN,'nvctr_parLIN',subname)
  nvctr_parLIN=0


  !create an array which indicate which processor has a GPU associated 
  !from the viewpoint of the BLAS routines
  allocate(GPU_for_comp(0:nproc-1+ndebug),stat=istat)
  call memocc(istat,GPU_for_comp,'GPU_for_comp',subname)

  if (nproc > 1 .and. .not. GPUshare) then
     call MPI_ALLGATHER(GPUblas,1,MPI_LOGICAL,GPU_for_comp(0),1,MPI_LOGICAL,&
          MPI_COMM_WORLD,ierr)
  else
     GPU_for_comp(0)=GPUblas
  end if

  iall=-product(shape(GPU_for_comp))*kind(GPU_for_comp)
  deallocate(GPU_for_comp,stat=istat)
  call memocc(istat,iall,'GPU_for_comp',subname)


  ! Distribute the orbitals belonging to outproc among the nproc processes in the communicator (done ny calling 
  ! 'parallel_repartition_with_kpoints'). The informations are stored in the array nvctr_parLIN. The meaning is
  ! the following:
  ! nvctr_parLIN(iorb,outproc,jproc,0)=ii means that orbital iorb of process outproc passes ii entries
  ! to process jproc when it is transposed.
  do outproc=lproc,uproc
     do iorb=1,norb_par(outproc,1) ! 1 for k-points
         call parallel_repartition_with_kpoints(nproc,orbs%nkpts,&
             (lin%wfds(iorb,outproc)%nvctr_c+7*lin%wfds(iorb,outproc)%nvctr_f),nvctr_par)
         do jproc=lproc,uproc
             nvctr_parLIN(iorb,outproc,jproc,0)=nvctr_par(jproc,0)
         end do
     end do
  end do


  ! Redistribute the orbitals among the processes considering k-points.
  ! If we have no k-points, this part does not change nvctr_parLIN.
  ! (If this is really the case, we could avoid this part using an if statement?)
  do outproc=lproc,uproc
    do iorb=1,norb_par(outproc,1) ! index 1 for k-points
        ikpts=1
        ncomp_res=(lin%wfds(iorb,outproc)%nvctr_c+7*lin%wfds(iorb,outproc)%nvctr_f)
        do jproc=lproc,uproc
           loop_comps: do
              if (nvctr_parLIN(iorb,outproc,jproc,0) >= ncomp_res) then
                 nvctr_parLIN(iorb,outproc,jproc,ikpts)= ncomp_res
                 ikpts=ikpts+1
                 nvctr_parLIN(iorb,outproc,jproc,0)=nvctr_parLIN(iorb,outproc,jproc,0)-ncomp_res
                 ncomp_res=(lin%wfds(iorb,outproc)%nvctr_c+7*lin%wfds(iorb,outproc)%nvctr_f)
              else
                 nvctr_parLIN(iorb,outproc,jproc,ikpts)= nvctr_parLIN(iorb,outproc,jproc,0)
                 if(nvctr_parLIN(iorb,outproc,jproc,ikpts)==0) write(*,'(a,3i6)') 'ATTENTION: iorb, outproc, jproc', iorb, &
                         outproc, jproc
                 ncomp_res=ncomp_res-nvctr_parLIN(iorb,outproc,jproc,0)
                 nvctr_parLIN(iorb,outproc,jproc,0)=0
                 exit loop_comps
              end if
              if (nvctr_parLIN(iorb,outproc,jproc,0) == 0 ) then
                 ncomp_res=(lin%wfds(iorb,outproc)%nvctr_c+7*lin%wfds(iorb,outproc)%nvctr_f)
                 exit loop_comps
              end if
           end do loop_comps
      end do
   end do
 end do
 

  ! Now we do some checks to make sure that the above distribution is correct.
  ! First we check whether the orbitals are correctly distributed among the MPI communicator
  ! when they are transposed.
  do ikpts=1,orbs%nkpts
    do outproc=lproc,uproc
      do iorb=1,orbs%norbp
         nvctr_tot=0
         do jproc=lproc,uproc
             nvctr_tot=nvctr_tot+nvctr_parLIN(iorb,outproc,jproc,ikpts)
         end do
         if(nvctr_tot /= lin%wfds(iorb,outproc)%nvctr_c+7*lin%wfds(iorb,outproc)%nvctr_f) then
            write(*,*)'ERROR: partition of components incorrect, iorb, kpoint:',iorb, ikpts
            stop
         end if
      end do
    end do
  end do
 
  ! Now we check whether the number of orbitals are correctly distributed.
  if (orbs%norb /= 0) then
     do ikpts=1,orbs%nkpts
        norb_tot=0
        do jproc=lproc,uproc
           norb_tot=norb_tot+norb_par(jproc,ikpts)
        end do
        if(norb_tot /= norbPerGroup) then
           write(*,*)'ERROR: partition of orbitals incorrect, kpoint:',ikpts
           write(*,'(a,3i8)') 'norbPerGroup, norb_tot, orbs%norb', norbPerGroup, norb_tot, orbs%norb
           stop
        end if
     end do
  end if
  
  
  ! WARNING: Make sure that this does not interfere with the 'normal' subroutine, since
  !          it might change the value of orbs%ikptproc(ikpts).
  !this function which associates a given k-point to a processor in the component distribution
  !the association is chosen such that each k-point is associated to only
  !one processor
  !if two processors treat the same k-point the processor which highest rank is chosen
  do ikpts=1,orbs%nkpts
     loop_jproc: do jproc=uproc,lproc,-1
        if (nvctr_par(jproc,ikpts) /= 0) then
           orbs%ikptproc(ikpts)=jproc
           exit loop_jproc
        end if
     end do loop_jproc
  end do
  
 
  ! WARNING: Make sure that this does not interfere with the 'normal' subroutine, since
  !          it might change the value of orbs%iskpts.
  !calculate the number of k-points treated by each processor in both
  ! the component distribution and the orbital distribution.
  nkptsp=0
  orbs%iskpts=-1
  do ikpts=1,orbs%nkpts
     if (nvctr_par(iproc,ikpts) /= 0 .or. norb_par(iproc,ikpts) /= 0) then
        if (orbs%iskpts == -1) orbs%iskpts=ikpts-1
        nkptsp=nkptsp+1
        mykpts(nkptsp) = ikpts
     end if
  end do
  orbs%nkptsp=nkptsp
 
 
  !print the distribution scheme ussed for this set of orbital
  !in the case of multiple k-points
  if (iproc == 0 .and. verbose > 1 .and. orbs%nkpts > 1) then
     call print_distribution_schemes(nproc,orbs%nkpts,norb_par(0,1),nvctr_par(0,1))
  end if
 
  !before printing the distribution schemes, check that the two distributions contain
  !the same k-points
  yesorb=.false.
  kpt_components: do ikpts=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikpts
     do jorb=1,orbs%norbp
        if (orbs%iokpt(jorb) == ikpt) yesorb=.true.
     end do
     if (.not. yesorb .and. orbs%norbp /= 0) then
        write(*,*)' ERROR: processor ', iproc,' kpt ',ikpt,&
             ' not found in the orbital distribution'
        stop
     end if
  end do kpt_components
 
  yescomp=.false.
  kpt_orbitals: do jorb=1,orbs%norbp
     ikpt=orbs%iokpt(jorb)   
     do ikpts=1,orbs%nkptsp
        if (orbs%iskpts+ikpts == ikpt) yescomp=.true.
     end do
     if (.not. yescomp) then
        write(*,*)' ERROR: processor ', iproc,' kpt,',ikpt,&
             'not found in the component distribution'
        stop
     end if
  end do kpt_orbitals


 
  ! Now comes the determination of the arrays needed for the communication.
  
  ! First copy the content of nvctr_parLIN to comms%nvctr_parLIN.
  allocate(comms%nvctr_parLIN(1:orbs%norb,lproc:uproc,lproc:uproc,orbs%nkptsp+ndebug),stat=istat)
  call memocc(istat,comms%nvctr_parLIN,'nvctr_parLIN',subname)
  !assign the partition of the k-points to the communication array
  do ikpts=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikpts!orbs%ikptsp(ikpts)
     do jproc=lproc,uproc
        do outproc=lproc,uproc
           do iorb=1,norb_par(outproc,ikpt)
              comms%nvctr_parLIN(iorb,outproc,jproc,ikpt)=nvctr_parLIN(iorb,outproc,jproc,ikpt) 
           end do
        end do
     end do
  end do
 
  ! Now come the send counts for the transposition (i.e. mpialltoallv).
  ! comms%ncntdLIN(jproc)=ii means that the current process (i.e. process iproc) passes
  ! totally ii elements to process jproc.
  allocate(comms%ncntdLIN(lproc:uproc), stat=istat)
  call memocc(istat, comms%ncntdLIN, 'comms%ncntdLIN', subname)
  comms%ncntdLIN=0
  do jproc=lproc,uproc
        do ikpts=1,orbs%nkpts
           ii=0
           do jorb=1,norb_par(iproc,ikpts)
               ii=ii+nvctr_parLIN(jorb,iproc,jproc,ikpts)*orbs%nspinor
           end do
           comms%ncntdLIN(jproc)=comms%ncntdLIN(jproc)+ii
        end do
  end do
 
  ! Now come the send displacements for the mpialltoallv.
  ! comms%ndspldLIN(jproc)=ii means that data sent from the current process (i.e. process iproc)
  ! to process jproc starts at location ii in the array hold by iproc.
  allocate(comms%ndspldLIN(lproc:uproc), stat=istat)
  call memocc(istat, comms%ndspldLIN, 'comms%ndspldLIN', subname)
  comms%ndspldLIN=0
  do jproc=lproc+1,uproc
     comms%ndspldLIN(jproc)=comms%ndspldLIN(jproc-1)+comms%ncntdLIN(jproc-1)
  end do
 
 
  ! Now come the receive counts for mpialltoallv.
  ! comms%ncnttLIN(jproc)=ii means that the current process (i.e. process iproc) receives
  ! ii elements from process jproc.
  allocate(comms%ncnttLIN(lproc:uproc), stat=istat)
  call memocc(istat, comms%ncnttLIN, 'comms%ncnttLIN', subname)
  comms%ncnttLIN=0
  do jproc=lproc,uproc
      do ikpts=1,orbs%nkpts
          ii=0
          do jorb=1,norb_par(jproc,ikpts)
              ii=ii+nvctr_parLIN(jorb,jproc,iproc,ikpts)*orbs%nspinor
          end do
          comms%ncnttLIN(jproc)=comms%ncnttLIN(jproc)+ii
      end do
  end do
  

  ! Now come the receive displacements for mpialltoallv.
  ! comms%ndspltLIN(jproc)=ii means that the data sent from process jproc to the current
  ! process (i.e. process iproc) start at the location ii in the array hold by iproc.
  allocate(comms%ndspltLIN(lproc:uproc), stat=istat)
  call memocc(istat, comms%ndspltLIN, 'comms%ndspltLIN', subname)
  comms%ndspltLIN=0
  do jproc=lproc+1,uproc
      comms%ndspltLIN(jproc)=comms%ndspltLIN(jproc-1)+comms%ncnttLIN(jproc-1)
  end do
 
 
  ! Deallocate the local arrays
  iall=-product(shape(nvctr_par))*kind(nvctr_par)
  deallocate(nvctr_par,stat=istat)
  call memocc(istat,iall,'nvctr_par',subname)
  iall=-product(shape(nvctr_parLIN))*kind(nvctr_parLIN)
  deallocate(nvctr_parLIN,stat=istat)
  call memocc(istat,iall,'nvctr_parLIN',subname)
  iall=-product(shape(norb_par))*kind(norb_par)
  deallocate(norb_par,stat=istat)
  call memocc(istat,iall,'norb_par',subname)
  iall=-product(shape(mykpts))*kind(mykpts)
  deallocate(mykpts,stat=istat)
  call memocc(istat,iall,'mykpts',subname)
 

  ! Calculate the dimension of the wave function for the given process.
  ! Take into account max one k-point per processor??
  ! WARNING: This changes the value of orbs%npsidim, which is then used in the context
  !          of the usual linear scaling version. Therefore this value will be changed back
  !          again at the end of the subroutine initializeLocRegLIN.
  ! Calculate the dimension of psi if it has to hodl its wave functions in the direct (i.e. not
  ! transposed) way.
  orbs%npsidim=0
  do iorb=1,orbs%norbp
      orbs%npsidim=orbs%npsidim+(lin%wfds(iorb,iproc)%nvctr_c+7*lin%wfds(iorb,iproc)%nvctr_f)*orbs%nspinor
  end do
  ! Eventually the dimension must be larger to hold all wavefunctions in the transposed way.
  ! Choose the maximum of these two numbers.
  orbs%npsidim=max(orbs%npsidim,sum(comms%ncnttLIN(lproc:uproc)))
 
  write(*,'(1x,a,i0,4x,i5)') &
       'LIN: Wavefunctions memory occupation for root processor (Bytes), iproc ',&
       orbs%npsidim*8, iproc
 

!END SUBROUTINE orbitals_communicatorsLIN_group
END SUBROUTINE orbitalsCommunicatorsWithGroups
!!***
