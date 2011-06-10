subroutine initializeLocRegLIN(iproc, nproc, lr, lin, at, input, rxyz, radii_cf)
!
! Purpose:
! ========
!   Determines the localization region for each orbital. In contrast to the usual
!   version, the localization regions can be different for each orbital.
!   In addition this subroutine creates the parameters needed for the communication
!   of the variable-size orbitals. At the end it performs a test of the transposition
!   to check the above initializations.
!   In addition this subroutine is able to split up the MPI process into different
!   communicators and perform the communication operations separatly for each
!   communicator.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc      process ID
!     nproc      number of processes
!     lr         type describing the localization region of the usual cubic version.
!                  This type is used to initialize some parameters of the new variable
!                  localization region descriptors. Maybe this can be done otherwise later?
!     at         type containing the paraneters for the atoms
!     input      type containing some very general parameters
!     rxyz       atomic positions
!     radii_cf   coase and fine radii around the atoms
!   Input / Output arguments
!   ------------------------
!     lin        type containing all parameters concerning the linear scaling version
!
use module_base
use module_types
use module_interfaces, exceptThisOne => initializeLocRegLIN
implicit none 

! Calling arguments
integer,intent(in):: iproc, nproc
type(locreg_descriptors),intent(in):: lr
type(linearParameters),intent(in out):: lin
type(atoms_data),intent(in):: at
type(input_variables),intent(in):: input
real(8),dimension(3,at%nat),intent(in):: rxyz
real(8),dimension(at%ntypes,3),intent(in):: radii_cf

! Local variables
integer:: iorb, iiAt, iitype, istat, iall, ierr, ii, ngroups, norbPerGroup, jprocStart, jprocEnd, i
integer:: norbtot1, norbtot2, igroup, jproc, norbpMax, lproc, uproc, wholeGroup
real(8):: radius, radiusCut, idealSplit, npsidimTemp
logical:: passedTest, passedTestAll
logical,dimension(:,:,:,:),allocatable:: logridCut_c
logical,dimension(:,:,:,:),allocatable:: logridCut_f
integer,dimension(:),allocatable:: norbPerGroupArr, newGroup, newComm, npsidimArr
integer,dimension(:,:),allocatable:: procsInGroup, newID
integer,dimension(:,:,:),allocatable:: tempArr
real(8),dimension(:),pointer:: psiInit, psi, psiWork
character(len=*),parameter:: subname='initializeLocRegLIN'


  ! WARNING: during this subroutine lin%orbs%npsidim may be modified. Therefore copy it here
  !          and assign it back at the end of the subroutine.
  !          Otherwise the linear scaling version which is executed after the call to this
  !          subroutine will not work any more!
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
  
  ! Allocate the wave function descriptors which will describe the orbitals.
  ! Since each orbital has its own localization region, we need an array.
  ! Since we allocate an array of types, there are problems with memocc...
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
      call fill_logridCut(lin%lr%geocode, lin%lr%d%n1, lin%lr%d%n2, lin%lr%d%n3, 0, lin%lr%d%n1, &
           0, lin%lr%d%n2, 0, lin%lr%d%n3, 0, 1, 1, 1, rxyz(1,iiAt), radius, radiusCut, &
           input%hx, input%hy, input%hz, logridCut_c(0,0,0,iorb))
  
      ! Calculate the number of segments and the number of grid points for each orbital.
      call num_segkeys(lin%lr%d%n1, lin%lr%d%n2, lin%lr%d%n3, 0, lin%lr%d%n1, 0, &
          lin%lr%d%n2, 0, lin%lr%d%n3, logridCut_c(0,0,0,iorb), lin%wfds(iorb,iproc)%nseg_c, &
          lin%wfds(iorb,iproc)%nvctr_c)
  
      ! Now the same procedure for the fine radius.
      radius=radii_cf(2,iitype)
      call fill_logridCut(lin%lr%geocode, lin%lr%d%n1, lin%lr%d%n2, lin%lr%d%n3, 0, lin%lr%d%n1, &
           0, lin%lr%d%n2, 0, lin%lr%d%n3, 0, 1, 1, 1, rxyz(1,iiAt), radius, radiusCut, &
           input%hx, input%hy, input%hz, logridCut_f(0,0,0,iorb))
  
      ! Calculate the number of segments and the number of grid points for each orbital.
      call num_segkeys(lin%lr%d%n1, lin%lr%d%n2, lin%lr%d%n3, 0, lin%lr%d%n1, 0, &
          lin%lr%d%n2, 0, lin%lr%d%n3, logridCut_f(0,0,0,iorb), lin%wfds(iorb,iproc)%nseg_f, &
          lin%wfds(iorb,iproc)%nvctr_f)
  
  
      ! Now fill the descriptors.
      call allocate_wfd(lin%wfds(iorb,iproc), subname)
      ! First the coarse part
      call segkeys(lin%lr%d%n1, lin%lr%d%n2, lin%lr%d%n3, 0, lin%lr%d%n1, 0, lin%lr%d%n2, &
          0, lin%lr%d%n3, logridCut_c(0,0,0,iorb), lin%wfds(iorb,iproc)%nseg_c, &
          lin%wfds(iorb,iproc)%keyg(1,1), lin%wfds(iorb,iproc)%keyv(1))
      ! And then the fine part
      ii=lin%wfds(iorb,iproc)%nseg_c+1
      call segkeys(lin%lr%d%n1, lin%lr%d%n2, lin%lr%d%n3, 0, lin%lr%d%n1, 0, lin%lr%d%n2, &
          0, lin%lr%d%n3, logridCut_f(0,0,0,iorb), lin%wfds(iorb,iproc)%nseg_f, &
          lin%wfds(iorb,iproc)%keyg(1,ii), lin%wfds(iorb,iproc)%keyv(ii))
  
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
  norbPerGroup=20
  !norbPerGroup=lin%orbs%norb
  
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
          norbtot1=norbtot1+lin%orbs%norb_par(jproc)
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
  
  if(iproc==0) write(*,'(x,a,i0,a)') 'Create ', lin%ncomms, ' new MPI communicators:'
  do igroup=1,lin%ncomms
      if(iproc==0) write(*,'(3x,a,i0,a,i0,a,i0,a)') '- Communicator ', igroup, ' includes the processes ', &
                       lin%procsInComm(1,igroup), ' to ', lin%procsInComm(2,igroup), '.'
  end do
  
  
  
  ! Now create the new MPI communicators.
  ! These communicators will be contained in the array newComm. If you want to
  ! use MPI processes only for the processes in group igroup, you can use the
  ! ordinary MPI routines just with newComm(igroup) instead of mpi_comm_world.
  allocate(newID(0:nproc,lin%ncomms), stat=istat)
  call memocc(istat, newID, 'newID', subname)
  allocate(newGroup(1:ngroups))
  call memocc(istat, newGroup, 'newGroup', subname)
  allocate(lin%MPIcomms(1:ngroups), stat=istat)
  call memocc(istat, lin%MPIComms, 'lin%MPIComms', subname)
  call mpi_comm_group(mpi_comm_world, wholeGroup, ierr)
  do igroup=1,ngroups
      do jproc=0,lin%procsInComm(2,igroup)-lin%procsInComm(1,igroup)+1
          newID(jproc,igroup)=lin%procsInComm(1,igroup)+jproc
      end do
      !call mpi_group_incl(wholeGroup, newID(lin%procsInComm(2,igroup))-newID(lin%procsInComm(1,igroup))+1,&
      !    newID(lin%procsInComm(1,igroup)), newGroup(igroup), ierr)
      call mpi_group_incl(wholeGroup, lin%procsInComm(2,igroup)-lin%procsInComm(1,igroup)+1,&
          newID(0,igroup), newGroup(igroup), ierr)
      call mpi_comm_create(mpi_comm_world, newGroup(igroup), lin%MPIComms(igroup), ierr)
  end do
  
  
  
  ! Now create the parameters for the transposition.
  ! lproc and uproc give the first and last process ID of the processes
  ! in the communicator igroup.
  if(iproc==0) write(*,'(x,a)') 'Wavefunction memory occupation:'
  do igroup=1,ngroups
      lproc=lin%procsInComm(1,igroup)
      uproc=lin%procsInComm(2,igroup)
      if(iproc>=lproc .and. iproc<=uproc) then
          call orbitalsCommunicatorsWithGroups(iproc, lproc, uproc, lin, lin%MPIComms(igroup), lin%norbPerComm(igroup))
      end if
  end do

  ! Write the memory occupation for the wave functions on each process.
  allocate(npsidimArr(0:nproc-1), stat=istat)
  call memocc(istat, npsidimArr, 'npsidimArr', subname)
  call mpi_gather(lin%orbs%npsidim, 1, mpi_integer, npsidimArr(0), 1, mpi_integer, 0, mpi_comm_world, ierr)
  do jproc=0,nproc-1
      if(iproc==0) write(*,'(3x,a,i0,a,i0,a)') '- process ', jproc, ': ', npsidimArr(jproc)*8, ' Bytes'
  end do
  iall=-product(shape(npsidimArr))*kind(npsidimArr)
  deallocate(npsidimArr, stat=istat)
  call memocc(istat, iall, 'npsidimArr', subname)
  
  
  
  ! Test the transposition. To do so, transpose and retranspose an array. If the transposition is
  ! correct, this should give the same result.
  if(iproc==0) write(*,'(x,a)') 'Testing whether the transposition works...'
  allocate(psi(lin%orbs%npsidim), stat=istat)
  call memocc(istat, psi, 'psi', subname)
  allocate(psiInit(lin%orbs%npsidim), stat=istat)
  call memocc(istat, psiInit, 'psiInit', subname)
  allocate(psiWork(lin%orbs%npsidim), stat=istat)
  call memocc(istat, psiWork, 'psiWork', subname)
  call random_number(psi)
  call dcopy(lin%orbs%npsidim, psi(1), 1, psiInit(1), 1)
  do igroup=1,ngroups
      lproc=lin%procsInComm(1,igroup)
      uproc=lin%procsInComm(2,igroup)
      if(iproc>=lproc .and. iproc<=uproc) then
          call transpose_vLIN(iproc, lproc, uproc, lin%orbs, lin%comms, psi, lin%MPIComms(igroup), work=psiWork)
          call untranspose_vLIN(iproc, lproc, uproc, lin%orbs, lin%comms, psi, lin%MPIComms(igroup), work=psiWork)
      end if
  end do
  
  ! Check whether the transposition was successful. First check on each process.
  passedTest=.true.
  ii=0
  do iorb=1,lin%orbs%norbp
      do i=1,lin%wfds(iorb,iproc)%nvctr_c+7*lin%wfds(iorb,iproc)%nvctr_f
          ii=ii+1
          if(psi(ii)/=psiInit(ii)) then
              passedTest=.false.
          end if
      end do
  end do
  
  ! Now check whether the test failes on any process.
  call mpi_allreduce(passedTest, passedTestAll, 1, mpi_logical, mpi_land, mpi_comm_world, ierr)
  if(passedTestAll) then
      if(iproc==0) write(*,'(3x,a)') 'Transposition test passed on all processes.'
  else
      do jproc=0,nproc-1
          if(iproc==jproc .and. .not.passedTest) write(*,'(3x,a,i0,a)') 'ERROR: transposition test failed on process ', iproc, '!'
      end do
      stop
  end if
  
  ! Deallocate all local arrays
  iall=-product(shape(psi))*kind(psi)
  deallocate(psi, stat=istat)
  call memocc(istat, iall, 'psi', subname)
  iall=-product(shape(psiInit))*kind(psiInit)
  deallocate(psiInit, stat=istat)
  call memocc(istat, iall, 'psiInit', subname)
  iall=-product(shape(psiWork))*kind(psiWork)
  deallocate(psiWork, stat=istat)
  call memocc(istat, iall, 'psiWork', subname)
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
  
  
  ! WARNING: assign back the original value of lin%orbs%npsidim due to the reasons
  !          explained in the beginning.
  lin%orbs%npsidim=npsidimTemp

end subroutine initializeLocRegLIN




subroutine fill_logridCut(geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,  &
     ntypes,iatype,rxyz,radii,rmult,hx,hy,hz,logrid)
!
! Purpose:
! ========
!   This subroutine creates a logical array logrid(i1,i2,i3) which specifies whether the grid point
!   i1,i2,i3 is the center of a scaling function/wavelet for the current orbital.
!   The criterion is based on the distance of the point i1,i2,i3 to the atom on which the current 
!   orbital is centered: If the distance is smaller than radii*rmult, then this point is within
!   the localization region (and thus carries a scaling function/wavelet), otherwise not.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     geocode    boundary conditions
!     n1         grid dimension in x direction is 0:n1
!     n2         grid dimension in y direction is 0:n2
!     n3         grid dimension in z direction is 0:n3
!     nl1        lower bound of the grid (in x direction) for which we determine the localization region
!     nl2        lower bound of the grid (in y direction) for which we determine the localization region
!     nl3        lower bound of the grid (in z direction) for which we determine the localization region
!     nu1        upper bound of the grid (in x direction) for which we determine the localization region
!     nu2        upper bound of the grid (in y direction) for which we determine the localization region
!     nu3        upper bound of the grid (in z direction) for which we determine the localization region
!     nbuf       ??
!     nat        number of atoms for which we should check wheter the grid point i1,i2,i3 is within 
!                  the localization region of these atoms. For the linear scaling version we should only
!                  take one atom, since the orbitals are centered on one atom.
!     ntypes     number of types of atom. For the linear scaling version, this should be equal to one for
!                the same reasons as nat.
!     iatype     array indicating to which atom type a given atom belongs. For the linear scaling version
!                this is one for again the same reasons.
!     rxyz       atomic positions
!     radii      part 1 of the cutoff radius (calculated by BigDFT)
!     rmult      part 2 of the cutoff radius (user specified)
!                -> the cutoff radius is then given by radii*rmult
!     hx         hgrid in x dimension
!     hy         hgrid in y dimension
!     hz         hgrid in z dimension
!   Output arguments:
!   -----------------
!     logrid     array indicating wheter a given grid point i1,i2,i3 is within the localization region or not.
!                If logrid(i1,i2,i3) is true, then the point i1,i2,i3 is within the region, otherwise not.
!
use module_base
implicit none

! Calling arguments
character(len=1),intent(in):: geocode
integer,intent(in):: n1, n2, n3, nl1, nu1, nl2, nu2, nl3, nu3, nbuf, nat, ntypes
real(gp),intent(in):: rmult, hx, hy, hz
integer,dimension(nat),intent(in):: iatype
real(gp),dimension(ntypes),intent(in):: radii
real(gp),dimension(3,nat),intent(in):: rxyz
logical,dimension(0:n1,0:n2,0:n3),intent(out):: logrid

! Local variables
real(kind=8),parameter:: eps_mach=1.d-12
integer:: i1, i2, i3, iat, ml1, ml2, ml3, mu1, mu2, mu3, j1, j2, j3
real(gp):: dx, dy2, dz2, rad


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





subroutine orbitalsCommunicatorsWithGroups(iproc, lproc, uproc, lin, newComm, norbPerComm)
!
! Purpose:
! ========
!   Creates the parameters needed for the communications, i.e. all orbitals belonging to the current communicator
!   are distributed among all processes in this communicator. The strategy is the same as for the normal
!   cubic approach:
!    -each processor has all the orbitals in transposed form
!    -each wavefunction is equally distributed in its transposed form
!   WARNING: This is subroutine is an adapted version of the original one, which includes also k-points.
!   Some fragments of this k-point structure are still present, but I think that there is still some
!   work to do if this subroutine should work for k-points!
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc         process ID
!     lproc         lowest process ID in the current communicator
!     uproc         highest process ID in the current communicator
!     newComm       the current MPI communicator
!     norbPerComm  number of orbitals in the current communicator
!   Input / Output arguments:
!   -------------------------
!     lin           type containing all parameters concerning the linear scaling version
!
use module_base
use module_types
implicit none

! Calling arguments
integer, intent(in) :: iproc, lproc, uproc, newComm, norbPerComm
type(linearParameters),intent(in out):: lin

! Local variables
logical:: yesorb, yescomp
integer:: jproc, nvctr_tot, ikpts, iorb, iorbp, jorb, norb_tot, ikpt, istat, iall, ii, outproc, nproc
integer:: ncomp_res, nkptsp, ierr
integer, dimension(:), allocatable:: mykpts
logical, dimension(:), allocatable:: GPU_for_comp
integer, dimension(:,:),allocatable:: nvctr_par, norb_par !for all the components and orbitals (with k-pts)
integer, dimension(:,:,:,:),allocatable:: nvctr_parLIN !for all the components and orbitals (with k-pts)
character(len=*),parameter:: subname='orbitalsCommunicatorsWithGroups'
  

  ! Number of processes in the current communicator
  nproc=uproc-lproc+1

  ! Check of allocation of important arrays
  if (.not. associated(lin%orbs%norb_par)) then
     write(*,*)'ERROR: norb_par array not allocated'
     stop
  end if
  

  ! Allocate the local arrays.
  allocate(nvctr_par(lproc:uproc,0:lin%orbs%nkpts+ndebug),stat=istat)
  call memocc(istat,nvctr_par,'nvctr_par',subname)
  allocate(norb_par(lproc:uproc,0:lin%orbs%nkpts+ndebug),stat=istat)
  call memocc(istat,norb_par,'norb_par',subname)
  allocate(mykpts(lin%orbs%nkpts+ndebug),stat=istat)
  call memocc(istat,mykpts,'mykpts',subname)

  ! Initialise the arrays
  do ikpts=0,lin%orbs%nkpts
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
     do iorbp=1,lin%orbs%norb_par(jproc)
        norb_par(jproc,ikpts)=norb_par(jproc,ikpts)+1
        if (mod(jorb,lin%orbs%norb)==0) then
           ikpts=ikpts+1
        end if
        jorb=jorb+1
     end do
  end do

  allocate(nvctr_parLIN(maxval(norb_par(:,1)),lproc:uproc,lproc:uproc,0:lin%orbs%nkpts+ndebug),stat=istat)
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
         call parallel_repartition_with_kpoints(nproc,lin%orbs%nkpts,&
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
  do ikpts=1,lin%orbs%nkpts
    do outproc=lproc,uproc
      do iorb=1,lin%orbs%norbp
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
  if (lin%orbs%norb /= 0) then
     do ikpts=1,lin%orbs%nkpts
        norb_tot=0
        do jproc=lproc,uproc
           norb_tot=norb_tot+norb_par(jproc,ikpts)
        end do
        if(norb_tot /= norbPerComm) then
           write(*,*)'ERROR: partition of orbitals incorrect, kpoint:',ikpts
           write(*,'(a,3i8)') 'norbPerComm, norb_tot, lin%orbs%norb', norbPerComm, norb_tot, lin%orbs%norb
           stop
        end if
     end do
  end if
  
  
  ! WARNING: Make sure that this does not interfere with the 'normal' subroutine, since
  !          it might change the value of lin%orbs%ikptproc(ikpts).
  !this function which associates a given k-point to a processor in the component distribution
  !the association is chosen such that each k-point is associated to only
  !one processor
  !if two processors treat the same k-point the processor which highest rank is chosen
  do ikpts=1,lin%orbs%nkpts
     loop_jproc: do jproc=uproc,lproc,-1
        if (nvctr_par(jproc,ikpts) /= 0) then
           lin%orbs%ikptproc(ikpts)=jproc
           exit loop_jproc
        end if
     end do loop_jproc
  end do
  
 
  ! WARNING: Make sure that this does not interfere with the 'normal' subroutine, since
  !          it might change the value of lin%orbs%iskpts.
  !calculate the number of k-points treated by each processor in both
  ! the component distribution and the orbital distribution.
  nkptsp=0
  lin%orbs%iskpts=-1
  do ikpts=1,lin%orbs%nkpts
     if (nvctr_par(iproc,ikpts) /= 0 .or. norb_par(iproc,ikpts) /= 0) then
        if (lin%orbs%iskpts == -1) lin%orbs%iskpts=ikpts-1
        nkptsp=nkptsp+1
        mykpts(nkptsp) = ikpts
     end if
  end do
  lin%orbs%nkptsp=nkptsp
 
 
  !print the distribution scheme ussed for this set of orbital
  !in the case of multiple k-points
  if (iproc == 0 .and. verbose > 1 .and. lin%orbs%nkpts > 1) then
     call print_distribution_schemes(nproc,lin%orbs%nkpts,norb_par(0,1),nvctr_par(0,1))
  end if
 
  !before printing the distribution schemes, check that the two distributions contain
  !the same k-points
  yesorb=.false.
  kpt_components: do ikpts=1,lin%orbs%nkptsp
     ikpt=lin%orbs%iskpts+ikpts
     do jorb=1,lin%orbs%norbp
        if (lin%orbs%iokpt(jorb) == ikpt) yesorb=.true.
     end do
     if (.not. yesorb .and. lin%orbs%norbp /= 0) then
        write(*,*)' ERROR: processor ', iproc,' kpt ',ikpt,&
             ' not found in the orbital distribution'
        stop
     end if
  end do kpt_components
 
  yescomp=.false.
  kpt_orbitals: do jorb=1,lin%orbs%norbp
     ikpt=lin%orbs%iokpt(jorb)   
     do ikpts=1,lin%orbs%nkptsp
        if (lin%orbs%iskpts+ikpts == ikpt) yescomp=.true.
     end do
     if (.not. yescomp) then
        write(*,*)' ERROR: processor ', iproc,' kpt,',ikpt,&
             'not found in the component distribution'
        stop
     end if
  end do kpt_orbitals


 
  ! Now comes the determination of the arrays needed for the communication.
  
  ! First copy the content of nvctr_parLIN to lin%comms%nvctr_parLIN.
  allocate(lin%comms%nvctr_parLIN(1:lin%orbs%norb,lproc:uproc,lproc:uproc,lin%orbs%nkptsp+ndebug),stat=istat)
  call memocc(istat,lin%comms%nvctr_parLIN,'nvctr_parLIN',subname)
  !assign the partition of the k-points to the communication array
  do ikpts=1,lin%orbs%nkptsp
     ikpt=lin%orbs%iskpts+ikpts!lin%orbs%ikptsp(ikpts)
     do jproc=lproc,uproc
        do outproc=lproc,uproc
           do iorb=1,norb_par(outproc,ikpt)
              lin%comms%nvctr_parLIN(iorb,outproc,jproc,ikpt)=nvctr_parLIN(iorb,outproc,jproc,ikpt) 
           end do
        end do
     end do
  end do
 
  ! Now come the send counts for the transposition (i.e. mpialltoallv).
  ! lin%comms%ncntdLIN(jproc)=ii means that the current process (i.e. process iproc) passes
  ! totally ii elements to process jproc.
  allocate(lin%comms%ncntdLIN(lproc:uproc), stat=istat)
  call memocc(istat, lin%comms%ncntdLIN, 'lin%comms%ncntdLIN', subname)
  lin%comms%ncntdLIN=0
  do jproc=lproc,uproc
        do ikpts=1,lin%orbs%nkpts
           ii=0
           do jorb=1,norb_par(iproc,ikpts)
               ii=ii+nvctr_parLIN(jorb,iproc,jproc,ikpts)*lin%orbs%nspinor
           end do
           lin%comms%ncntdLIN(jproc)=lin%comms%ncntdLIN(jproc)+ii
        end do
  end do
 
  ! Now come the send displacements for the mpialltoallv.
  ! lin%comms%ndspldLIN(jproc)=ii means that data sent from the current process (i.e. process iproc)
  ! to process jproc starts at location ii in the array hold by iproc.
  allocate(lin%comms%ndspldLIN(lproc:uproc), stat=istat)
  call memocc(istat, lin%comms%ndspldLIN, 'lin%comms%ndspldLIN', subname)
  lin%comms%ndspldLIN=0
  do jproc=lproc+1,uproc
     lin%comms%ndspldLIN(jproc)=lin%comms%ndspldLIN(jproc-1)+lin%comms%ncntdLIN(jproc-1)
  end do
 
 
  ! Now come the receive counts for mpialltoallv.
  ! lin%comms%ncnttLIN(jproc)=ii means that the current process (i.e. process iproc) receives
  ! ii elements from process jproc.
  allocate(lin%comms%ncnttLIN(lproc:uproc), stat=istat)
  call memocc(istat, lin%comms%ncnttLIN, 'lin%comms%ncnttLIN', subname)
  lin%comms%ncnttLIN=0
  do jproc=lproc,uproc
      do ikpts=1,lin%orbs%nkpts
          ii=0
          do jorb=1,norb_par(jproc,ikpts)
              ii=ii+nvctr_parLIN(jorb,jproc,iproc,ikpts)*lin%orbs%nspinor
          end do
          lin%comms%ncnttLIN(jproc)=lin%comms%ncnttLIN(jproc)+ii
      end do
  end do
  

  ! Now come the receive displacements for mpialltoallv.
  ! lin%comms%ndspltLIN(jproc)=ii means that the data sent from process jproc to the current
  ! process (i.e. process iproc) start at the location ii in the array hold by iproc.
  allocate(lin%comms%ndspltLIN(lproc:uproc), stat=istat)
  call memocc(istat, lin%comms%ndspltLIN, 'lin%comms%ndspltLIN', subname)
  lin%comms%ndspltLIN=0
  do jproc=lproc+1,uproc
      lin%comms%ndspltLIN(jproc)=lin%comms%ndspltLIN(jproc-1)+lin%comms%ncnttLIN(jproc-1)
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
  ! WARNING: This changes the value of lin%orbs%npsidim, which is then used in the context
  !          of the usual linear scaling version. Therefore this value will be changed back
  !          again at the end of the subroutine initializeLocRegLIN.
  ! Calculate the dimension of psi if it has to hodl its wave functions in the direct (i.e. not
  ! transposed) way.
  lin%orbs%npsidim=0
  do iorb=1,lin%orbs%norbp
      lin%orbs%npsidim=lin%orbs%npsidim+(lin%wfds(iorb,iproc)%nvctr_c+7*lin%wfds(iorb,iproc)%nvctr_f)*lin%orbs%nspinor
  end do
  ! Eventually the dimension must be larger to hold all wavefunctions in the transposed way.
  ! Choose the maximum of these two numbers.
  lin%orbs%npsidim=max(lin%orbs%npsidim,sum(lin%comms%ncnttLIN(lproc:uproc)))
 


END SUBROUTINE orbitalsCommunicatorsWithGroups


!determine a set of localisation regions from the centers and the radii.
!cut in cubes the global reference system
subroutine determine_locreg_periodic(iproc,nlr,cxyz,locrad,hx,hy,hz,Glr,Llr)!,outofzone)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc
  integer, intent(in) :: nlr
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: Glr
  real(gp), dimension(nlr), intent(in) :: locrad
  real(gp), dimension(3,nlr), intent(in) :: cxyz
  type(locreg_descriptors), dimension(nlr), intent(out) :: Llr
!  integer, dimension(3,nlr),intent(out) :: outofzone
  !local variables
  character(len=*), parameter :: subname='determine_locreg'
  logical :: Gperx,Gpery,Gperz,Lperx,Lpery,Lperz
  logical :: warningx,warningy,warningz
  integer :: Gnbl1,Gnbl2,Gnbl3,Gnbr1,Gnbr2,Gnbr3
  integer :: Lnbl1,Lnbl2,Lnbl3,Lnbr1,Lnbr2,Lnbr3
  integer :: ilr,isx,isy,isz,iex,iey,iez
  integer :: ln1,ln2,ln3
  integer :: ii !tests
  integer,dimension(3) :: outofzone
  real(gp) :: rx,ry,rz,cutoff  

  if (iproc == 0) then
     write(*,*)'Inside determine_locreg_periodic:'
  end if

  !initialize out of zone and logicals
  outofzone (:) = 0     
  warningx = .false.
  warningy = .false.
  warningz = .false.  

  !determine the limits of the different localisation regions
  do ilr=1,nlr

     rx=cxyz(1,ilr)
     ry=cxyz(2,ilr)
     rz=cxyz(3,ilr)

     cutoff=locrad(ilr)

     isx=floor((rx-cutoff)/hx)
     isy=floor((ry-cutoff)/hy)
     isz=floor((rz-cutoff)/hz)

     iex=ceiling((rx+cutoff)/hx)
     iey=ceiling((ry+cutoff)/hy)
     iez=ceiling((rz+cutoff)/hz)

     ln1 = iex-isx
     ln2 = iey-isy
     ln3 = iez-isz

     ! First check if localization region fits inside box
!!!     if (iproc == 0 .and. verbose > 1) then
!!!        if ((iex - isx >= Glr%d%n1 - 14) .and. (warningx .eqv. .false.)) then
!!!           write(*,*)'Width of direction x :',(iex - isx)*hx,' of localization region:',ilr
!!!           write(*,*)'is close or exceeds to the width of the simulation box:',Glr%d%n1*hx
!!!           write(*,*)'Increasing the simulation box is recommended. The code will use the '
!!!           write(*,*)'simulation box width. This is the only warning for x direction.'
!!!           warningx = .true.
!!!        end if
!!!        if ((iey - isy >= Glr%d%n2 - 14) .and. (warningy .eqv. .false.)) then
!!!           write(*,*)'Width of direction y :',(iey - isy)*hy,' of localization region:',ilr
!!!           write(*,*)'is close or exceeds to the width of the simulation box:',Glr%d%n2*hy,'.'
!!!           write(*,*)'Increasing the simulation box is recommended. The code will use the width'
!!!           write(*,*)'of the simulation box. This is the only warning for y direction.'
!!!           warningy = .true.
!!!        end if
!!!        if ((iez - isz >= Glr%d%n3 - 14) .and. (warningz .eqv. .false.)) then
!!!           write(*,*)'Width of direction z :',(iez - isz)*hz,' of localization region:',ilr
!!!           write(*,*)'is close or exceeds to the width of the simulation box:',Glr%d%n3*hz,'.'
!!!           write(*,*)'Increasing the simulation box is recommended. The code will use the width'
!!!           write(*,*)'of the simulation box. This is the only warning for z direction.'
!!!           warningz = .true.
!!!        end if 
!!!     end if

     ! Localization regions should always have free boundary conditions
     Llr(ilr)%geocode='F'

     !assign the starting/ending points and outofzone for the different
     ! geometries
     select case(Glr%geocode)
     case('F')
        isx=max(isx,Glr%ns1)
        isy=max(isy,Glr%ns2)
        isz=max(isz,Glr%ns3)

        iex=min(iex,Glr%ns1+Glr%d%n1)
        iey=min(iey,Glr%ns2+Glr%d%n2)
        iez=min(iez,Glr%ns3+Glr%d%n3)

     case('S')
        ! Get starting and ending for x direction     
        if (iex - isx >= Glr%d%n1) then       
           isx=Glr%ns1
           iex=Glr%ns1 + Glr%d%n1
        else
           isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
           iex= ln1 + isx
           if (iex > Glr%ns1+Glr%d%n1) then
              outofzone(1)=modulo(iex,Glr%d%n1+1)
           end if           
        end if
        
        ! Get starting and ending for y direction (perpendicular to surface)
        isy=max(isy,Glr%ns2)
        iey=min(iey,Glr%ns2 + Glr%d%n2)
        outofzone(2) = 0

        !Get starting and ending for z direction
        if (iez - isz >= Glr%d%n3) then
           isz=Glr%ns3 
           iez=Glr%ns3 + Glr%d%n3
        else
           isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
           iez= ln3 + isz
           if (iez > Glr%ns3+Glr%d%n3) then
              outofzone(3)=modulo(iez,Glr%d%n3+1)
           end if 
        end if

     case('P')
         ! Get starting and ending for x direction     
        if (iex - isx >= Glr%d%n1) then       
           isx=Glr%ns1
           iex=Glr%ns1 + Glr%d%n1
        else
           isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
           iex= ln1 + isx
           if (iex > Glr%ns1+Glr%d%n1) then
              outofzone(1)=modulo(iex,Glr%d%n1+1)
           end if           
        end if
        
        ! Get starting and ending for y direction (perpendicular to surface)
        if (iey - isy >= Glr%d%n2) then       
           isy=Glr%ns2
           iey=Glr%ns2 + Glr%d%n2
         else
           isy=modulo(isy,Glr%d%n2+1) + Glr%ns2
           iey= ln2 + isy
           if (iey > Glr%ns2+Glr%d%n2) then
              outofzone(2)=modulo(iey,Glr%d%n2+1)
           end if           
        end if

        !Get starting and ending for z direction
        if (iez - isz >= Glr%d%n3) then
           isz=Glr%ns3 
           iez=Glr%ns3 + Glr%d%n3
        else
           isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
           iez= ln3 + isz
           if (iez > Glr%ns3+Glr%d%n3) then
              outofzone(3)=modulo(iez,Glr%d%n3+1)
           end if 
        end if
     end select

     !values for the starting point of the cube for wavelet grid
     Llr(ilr)%ns1=isx
     Llr(ilr)%ns2=isy
     Llr(ilr)%ns3=isz

     !dimensions of the localisation region
     Llr(ilr)%d%n1=iex-isx
     Llr(ilr)%d%n2=iey-isy
     Llr(ilr)%d%n3=iez-isz

     !assign outofzone
     Llr(ilr)%outofzone(:) = outofzone(:)

     ! Set the conditions for ext_buffers (conditions for buffer size)
     Gperx=(Glr%geocode /= 'F')
     Gpery=(Glr%geocode == 'P')
     Gperz=(Glr%geocode /= 'F')
     Lperx=(Llr(ilr)%geocode /= 'F')
     Lpery=(Llr(ilr)%geocode == 'P')
     Lperz=(Llr(ilr)%geocode /= 'F')

     !calculate the size of the buffers of interpolating function grid
     call ext_buffers(Gperx,Gnbl1,Gnbr1)
     call ext_buffers(Gpery,Gnbl2,Gnbr2)
     call ext_buffers(Gperz,Gnbl3,Gnbr3)
     call ext_buffers(Lperx,Lnbl1,Lnbr1)
     call ext_buffers(Lpery,Lnbl2,Lnbr2)
     call ext_buffers(Lperz,Lnbl3,Lnbr3)

     !starting point of the region for interpolating functions grid
     Llr(ilr)%nsi1= 2 * Llr(ilr)%ns1 - (Lnbl1 - Gnbl1)
     Llr(ilr)%nsi2= 2 * Llr(ilr)%ns2 - (Lnbl2 - Gnbl2)
     Llr(ilr)%nsi3= 2 * Llr(ilr)%ns3 - (Lnbl3 - Gnbl3)

     !dimensions of the fine grid inside the localisation region
     Llr(ilr)%d%nfl1=max(isx,Glr%d%nfl1)-isx ! should we really substract isx (probably because the routines are coded with 0 as origin)?
     Llr(ilr)%d%nfl2=max(isy,Glr%d%nfl2)-isy
     Llr(ilr)%d%nfl3=max(isz,Glr%d%nfl3)-isz
     
     !NOTE: This will not work with symmetries (must change it)
     Llr(ilr)%d%nfu1=min(iex,Glr%d%nfu1)-isx
     Llr(ilr)%d%nfu2=min(iey,Glr%d%nfu2)-isy
     Llr(ilr)%d%nfu3=min(iez,Glr%d%nfu3)-isz

     !dimensions of the interpolating scaling functions grid (reduce to +2?, check with Luigi)
     Llr(ilr)%d%n1i=2*Llr(ilr)%d%n1+31
     Llr(ilr)%d%n2i=2*Llr(ilr)%d%n2+31
     Llr(ilr)%d%n3i=2*Llr(ilr)%d%n3+31

!DEBUG
!     write(*,*)'Description of zone:',ilr
!     write(*,*)'ns:',Llr(ilr)%ns1,Llr(ilr)%ns2,Llr(ilr)%ns3
!     write(*,*)'ne:',Llr(ilr)%ns1+Llr(ilr)%d%n1,Llr(ilr)%ns2+Llr(ilr)%d%n2,Llr(ilr)%ns3+Llr(ilr)%d%n3
!     write(*,*)'n:',Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3
!     write(*,*)'nfl:',Llr(ilr)%d%nfl1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfl3
!     write(*,*)'nfu:',Llr(ilr)%d%nfu1,Llr(ilr)%d%nfu2,Llr(ilr)%d%nfu3
!     write(*,*)'ni:',Llr(ilr)%d%n1i,Llr(ilr)%d%n2i,Llr(ilr)%d%n3i
!     write(*,*)'outofzone',ilr,':',outofzone(:)
!DEBUG

    ! construct the wavefunction descriptors (wfd)
     call determine_wfd_periodicity(ilr,nlr,Glr,Llr,outofzone)
     
!     print *,'Before locreg_bounds'
!     print *,'n:',Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3
!     print *,'nl,nu:',Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3
!     print *,'wfd(nseg):',Llr(ilr)%wfd%nseg_c,Llr(ilr)%wfd%nseg_f
!     print *,'wfd(nvctr):',Llr(ilr)%wfd%nvctr_c,Llr(ilr)%wfd%nvctr_f

     ! Sould check if nfu works properly... also relative to locreg!!
     !if the localisation region is isolated build also the bounds
     if (Llr(ilr)%geocode=='F') then
        call locreg_bounds(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,&
             Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,&
             Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3,Llr(ilr)%wfd,Llr(ilr)%bounds)
     end if
  end do !on ilr

  !after all localisation regions are determined draw them
  !call draw_locregs(nlr,hx,hy,hz,Llr)

END SUBROUTINE determine_locreg_periodic


!#############################################################################################################################################
!!****f* BigDFT/overlap_region
!#############################################################################################################################################
!! FUNCTION: Determines the the wavefunction descriptors,wfd, and fine grid upper limit of locreg taking into account the pediodicity
!!          
!! WARNING: We assign Llr%nfl and llr%nfu with respect to the origin of the local zone, like in determine_locreg. 
!!         
!! SOURCE:
!!
subroutine determine_wfd_periodicity(ilr,nlr,Glr,Llr)!,outofzone)

  use module_base
  use module_types
 
  implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: ilr,nlr
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),dimension(nlr),intent(inout) :: Llr  ! Localization grid descriptors 
  !########################################
  !Subroutine Array Arguments
  !########################################
!  integer,dimension(3,nlr),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr
  !#############################################
  !local variables
  !############################################
  integer :: ii
  integer,dimension(3) :: Gife,Gifs,iedir,isdir,Lifs,Life,period
  integer :: nseg_c,nseg_f,nvctr_c,nvctr_f      ! total number of sgements and elements
  character(len=*), parameter :: subname='determine_wfd_periodicity'

   !starting point of locreg (always inside locreg)
   isdir(1) = Llr(ilr)%ns1
   isdir(2) = Llr(ilr)%ns2
   isdir(3) = Llr(ilr)%ns3
   !ending point of locreg (can be outside the simulation box)
   iedir(1) = Llr(ilr)%ns1 + Llr(ilr)%d%n1
   iedir(2) = Llr(ilr)%ns2 + Llr(ilr)%d%n2
   iedir(3) = Llr(ilr)%ns3 + Llr(ilr)%d%n3
   ! starting and ending point of fine grid in Global region
   Gifs(1) = Glr%d%nfl1 + Glr%ns1
   Gifs(2) = Glr%d%nfl2 + Glr%ns2
   Gifs(3) = Glr%d%nfl3 + Glr%ns3
   Gife(1) = Glr%d%nfu1 + Glr%ns1
   Gife(2) = Glr%d%nfu2 + Glr%ns2
   Gife(3) = Glr%d%nfu3 + Glr%ns3
   ! periodicity
   period(1) = Glr%d%n1
   period(2) = Glr%d%n2
   period(3) = Glr%d%n3

   ! Determine starting point of the fine grid in locreg
   do ii=1,3
      if (Llr(ilr)%outofzone(ii) > 0) then
         ! When periodicity, we must check for 2 different situations:
         ! (1) : starting of locreg before or in fine grid zone
         if (isdir(ii) < Gife(ii)) Lifs(ii) = max(isdir(ii),Gifs(ii))-isdir(ii)
         ! (2) : starting point after fine grid
         if (isdir(ii) > Gife(ii)) Lifs(ii) = max(isdir(ii),Gifs(ii)+period(ii))-isdir(ii)
      else
          Lifs(ii) = max(isdir(ii),Gifs(ii))-isdir(ii)
      end if 
   end do

   ! Determine ending point of the fine grid in locreg
   do ii=1,3
      if(Llr(ilr)%outofzone(ii) > 0) then
         !When periodicity, we must check for three different situations:
         ! (1) : ending of locreg before fine grid zone
         if(iedir(ii) < (Gifs(ii) + period(ii))) Life(ii) = Gife(ii)-isdir(ii)
         ! (2) : ending of locreg in fine grid zone
         if(iedir(ii) > (Gifs(ii) + period(ii)) .and. iedir(ii) < (Gife(ii) + period(ii))) then
           Life(ii) = iedir(ii)-isdir(ii)
         end if
         ! (3) : ending of locreg after ending of fine grid zone
         if(iedir(ii) > (Gife(ii)+period(ii))) Life(ii) = Gife(ii) + period(ii)-isdir(ii)
      else
         Life(ii) = min(iedir(ii),Gife(ii))-isdir(ii)
      end if
   end do

   ! Assign values to Llr
   Llr(ilr)%d%nfl1 = Lifs(1)
   Llr(ilr)%d%nfl2 = Lifs(2)
   Llr(ilr)%d%nfl3 = Lifs(3)
   Llr(ilr)%d%nfu1 = Life(1)
   Llr(ilr)%d%nfu2 = Life(2)
   Llr(ilr)%d%nfu3 = Life(3)

   ! define the wavefunction descriptors inside the localisation region
   !coarse part
   call num_segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),isdir(2),&
          iedir(2),isdir(3),iedir(3),Glr%wfd%nseg_c,Glr%wfd%nvctr_c,&
          Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),nseg_c,nvctr_c,Llr(ilr)%outofzone(:))
   !fine part
   call num_segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),isdir(2),&
          iedir(2),isdir(3),iedir(3),Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
          Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),nseg_f,nvctr_f,Llr(ilr)%outofzone(:))

   ! Assign the values to Llr
   Llr(ilr)%wfd%nseg_c = nseg_c
   Llr(ilr)%wfd%nseg_f = nseg_f
   Llr(ilr)%wfd%nvctr_c= nvctr_c
   Llr(ilr)%wfd%nvctr_f= nvctr_f

   !allocate the wavefunction descriptors following the needs
   call allocate_wfd(Llr(ilr)%wfd,subname)

   !Now, fill the descriptors:
   !coarse part
   call segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),&
        isdir(2),iedir(2),isdir(3),iedir(3),&
        Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
        Llr(ilr)%wfd%nseg_c,Llr(ilr)%wfd%nvctr_c,&
        Llr(ilr)%wfd%keyg(1,1),Llr(ilr)%wfd%keyv(1),Llr(ilr)%outofzone(:))

   !fine part
   call segkeys_periodic(Glr%d%n1,Glr%d%n2,Glr%d%n3,isdir(1),iedir(1),&
        isdir(2),iedir(2),isdir(3),iedir(3),Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
        Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
        Llr(ilr)%wfd%nseg_f,Llr(ilr)%wfd%nvctr_f,&
        Llr(ilr)%wfd%keyg(1,Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),&
        Llr(ilr)%wfd%keyv(Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),Llr(ilr)%outofzone(:))

END SUBROUTINE determine_wfd_periodicity


!#############################################################################################################################################
!!****f* BigDFT/num_segkeys_periodic
!#############################################################################################################################################
!! FUNCTION: Calculates the number of segments and elements in localisation region
!!          
!! WARNING:   
!!         
!! SOURCE:
!!
subroutine num_segkeys_periodic(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc,outofzone)
  implicit none
  integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr
  integer, dimension(nseg), intent(in) :: keyv
  integer, dimension(2,nseg), intent(in) :: keyg
  integer, intent(out) :: nseg_loc,nvctr_loc
  integer, dimension(3),intent(in) :: outofzone
  !local variables
  logical :: lseg,go1,go2,go3
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,nsrt,nend,nvctr_check

  nvctr_loc=0
  !control variable
  nvctr_check=0
  !start and end points
  nsrt=0
  nend=0

  do iseg=1,nseg
     jj=keyv(iseg)
     j0=keyg(1,iseg)
     j1=keyg(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     lseg=.false.
     ! overlap conditions if zone completely inside simulation box
     go2 = (i2sc <= i2 .and. i2 <= i2ec)
     go3 = (i3sc <= i3 .and. i3 <= i3ec)
     ! overlap conditions if zone as components in other periodic cells
     if(outofzone(2) > 0) go2 = (i2 <= outofzone(2) .or. i2 >= i2sc)
     if(outofzone(3) > 0) go3 = (i3 <= outofzone(3) .or. i3 >= i3sc)

     do i=i0,i1
        nvctr_check=nvctr_check+1
        go1 = (i1sc <= i .and. i <= i1ec)
        if(outofzone(1) > 0) go1 = (i <= outofzone(1) .or. i >= i1sc)

        if (go1 .and. go2 .and. go3 ) then
           nvctr_loc=nvctr_loc+1
           if (.not. lseg) then
              nsrt=nsrt+1
           end if
           lseg=.true.
        else
           if (lseg) then
              nend=nend+1
              lseg=.false.
           end if
        end if
     end do
     if (lseg) then
        nend=nend+1
     end if
  end do
  nseg_loc=nend

  !check
  if (nend /= nsrt) then
     write(*,*) 'nend , nsrt',nend,nsrt
     stop 'nend <> nsrt'
  endif

  if (nvctr_check /= nvctr) then
     write(*,'(1x,a,2(i6))')&
          'ERROR: incorrect number of coarse points examined for reducing the localisation region',&
          nvctr_check,nvctr
     stop
  end if

END SUBROUTINE num_segkeys_periodic

subroutine segkeys_periodic(n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,keyg,keyv,&
     nseg_loc,nvctr_loc,keyg_loc,keyv_loc,outofzone)
  implicit none
  integer, intent(in) :: n1,n2,n3,i1sc,i1ec,i2sc,i2ec,i3sc,i3ec,nseg,nvctr,nseg_loc,nvctr_loc
  integer, dimension(nseg), intent(in) :: keyv
  integer, dimension(2,nseg), intent(in) :: keyg
  integer, dimension(3), intent(in) :: outofzone
  integer, dimension(nseg_loc), intent(out) :: keyv_loc
  integer, dimension(2,nseg_loc), intent(out) :: keyg_loc
  !local variables
  logical :: go1,go2,go3,lseg
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i,ind,nsrt,nend,nvctr_check,n1l,n2l,n3l,i1l,i2l,i3l
  integer :: ngridp

  !dimensions of the localisation region (O:nIl)
  ! must be smaller or equal to simulation box dimensions
  n1l=i1ec-i1sc
  n2l=i2ec-i2sc
  n3l=i3ec-i3sc

  !control variable
  nvctr_check=0
  !start and end points
  nsrt=0
  nend=0
  do iseg=1,nseg
     jj=keyv(iseg)
     j0=keyg(1,iseg)
     j1=keyg(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     lseg=.false.

     ! intersection condition if zone inside simulation box
     go2 = (i2sc <= i2 .and. i2 <= i2ec)
     go3 = (i3sc <= i3 .and. i3 <= i3ec) 
     ! intersection condition if zone has components outside simulation box (periodic)
     if(outofzone(2) > 0) go2 = (i2 <= outofzone(2) .or. i2 >= i2sc)
     if(outofzone(3) > 0) go3 = (i3 <= outofzone(3) .or. i3 >= i3sc)

     do i=i0,i1
        go1 = (i1sc <= i .and. i <= i1ec)
        if(outofzone(1) > 0) go1 = (i <= outofzone(1) .or. i >= i1sc)
        if (go1 .and. go2 .and. go3) then
          !index of the compressed function
          i1l=i-i1sc
          if(outofzone(1) > 0 .and. i <= outofzone(1))i1l = i - i1sc + n1 + 1  
          i2l=i2-i2sc
          if(outofzone(2) > 0 .and. i2 <= outofzone(2))i2l = i2 - i2sc + n2 + 1
          i3l=i3-i3sc
          if(outofzone(3) > 0 .and. i3 <= outofzone(3))i3l = i3 - i3sc + n3 + 1
          ngridp=i3l*((n1l+1)*(n2l+1)) + i2l*(n1l+1) + i1l+1

          nvctr_check=nvctr_check+1
          if (.not. lseg) then
!             print *,'         check:',i,i2,i3,i1l,i2l,i3l,ngridp
             nsrt=nsrt+1
             keyg_loc(1,nsrt)=ngridp
             keyv_loc(nsrt)=nvctr_check
          end if
          lseg=.true.
        else 
           if (lseg) then
!              print *,'in        else:',i,i2,i3,i1l,i2l,i3l,ngridp
              nend=nend+1
              keyg_loc(2,nend)=ngridp
              lseg=.false.
           end if
        end if
     end do
     if (lseg) then
!        print *,'in second else:',i,i2,i3,i1l,i2l,i3l,ngridp
        nend=nend+1
        keyg_loc(2,nend)=ngridp
     end if
  end do

  !check
  if (nvctr_check /= nvctr_loc .or. nend /= nsrt .or. nend /= nseg_loc) then
     print *,'global region statistics:',nseg,nvctr
     write(*,*)&
          'ERROR: problem in segkeys_periodic  ',&
          'nvctr_check:',nvctr_check,'nvctr_loc:',nvctr_loc,&
          'nend:',nend,'nsrt:',nsrt,'nseg_loc:',nseg_loc
     stop
  end if

END SUBROUTINE segkeys_periodic

!#############################################################################################################################################
!!****f* BigDFT/overlap_region
!#############################################################################################################################################
!! FUNCTION: Determines the number of intersection regions between locregs, taking into account the periodicity of the system.
!!          
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine get_number_of_overlap_region(alr,blr,Glr,isovrlp,Llr,nlr)!,outofzone)

  use module_base
  use module_types
 
  implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer, intent(in) :: alr,blr              ! index of the two localization regions
  integer, intent(in) :: nlr                  ! number of localization regions
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  integer, intent(out) :: isovrlp             ! Integer giving the number of overlaps (max 8 with periodicity)
  !########################################
  !Subroutine Array Arguments
  !########################################
!  integer,dimension(3,nlr),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr
  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr  ! Localization grid descriptors 
  !#############################################
  !local variables
  !############################################
  integer :: ii,azones,bzones,i_stat,i_all
  integer :: izones,jzones
  integer,allocatable :: astart(:,:),aend(:,:),bstart(:,:),bend(:,:)
  character(len=*), parameter :: subname='get_number_of_overlap_region'
  logical :: go1,go2,go3

  azones = 1
  bzones = 1
! Calculate the number of regions to cut alr and blr
  do ii=1,3
     if(Llr(alr)%outofzone(ii) > 0) azones = azones * 2
     if(Llr(blr)%outofzone(ii) > 0) bzones = bzones * 2
  end do

!write(*,*)'azones,bzones',azones,bzones
!write(*,*)'outofzone',alr,':',outofzone(:,alr)
!write(*,*)'outofzone',blr,':',outofzone(:,blr)

!allocate astart and aend
  allocate(astart(3,azones),stat=i_stat)
  call memocc(i_stat,astart,'astart',subname)
  allocate(aend(3,azones),stat=i_stat)
  call memocc(i_stat,aend,'aend',subname)

!FRACTURE THE FIRST LOCALIZATION REGION
  call fracture_periodic_zone(azones,Glr,Llr(alr),Llr(alr)%outofzone(:),astart,aend)

!allocate bstart and bend
  allocate(bstart(3,bzones),stat=i_stat)
  call memocc(i_stat,bstart,'bstart',subname)
  allocate(bend(3,bzones),stat=i_stat)
  call memocc(i_stat,bend,'bend',subname)

!FRACTURE SECOND LOCREG
  call fracture_periodic_zone(bzones,Glr,Llr(blr),Llr(blr)%outofzone(:),bstart,bend)

! Now check the number of overlapping zones
 isovrlp = 0
 do izones=1,azones
   do jzones=1,bzones
      go1 = (bstart(1,jzones) .le. aend(1,izones) .and. bend(1,jzones) .ge. astart(1,izones)) 
      go2 = (bstart(2,jzones) .le. aend(2,izones) .and. bend(2,jzones) .ge. astart(2,izones)) 
      go3 = (bstart(3,jzones) .le. aend(3,izones) .and. bend(3,jzones) .ge. astart(3,izones)) 
      if(go1 .and. go2 .and. go3) then
         isovrlp = isovrlp + 1
      end if
   end do
 end do

! Deallocation block
  i_all = -product(shape(astart))*kind(astart)
  deallocate(astart,stat=i_stat)
  call memocc(i_stat,i_all,'astart',subname)
  i_all = -product(shape(aend))*kind(aend)
  deallocate(aend,stat=i_stat)
  call memocc(i_stat,i_all,'aend',subname)
  i_all = -product(shape(bstart))*kind(bstart)
  deallocate(bstart,stat=i_stat)
  call memocc(i_stat,i_all,'bstart',subname)
  i_all = -product(shape(bend))*kind(bend)
  deallocate(bend,stat=i_stat)
  call memocc(i_stat,i_all,'bend',subname)
  
END SUBROUTINE get_number_of_overlap_region

!#############################################################################################################################################
!!****f* BigDFT/fracture_periodic_zone
!#############################################################################################################################################
!! FUNCTION: Divides the locreg into zones contained inside the simulation box, by applying the primitive vectors
!!           It returns: astart(3,nzones) which is the starting points of the different zones (max. 8)
!!                       aend(3,nzones) which is the ending points of the different zones (max. 8)
!!          
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine fracture_periodic_zone(nzones,Glr,Llr,outofzone,astart,aend)

  use module_base
  use module_types
 
  implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: nzones
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),intent(in) :: Llr  ! Localization grid descriptors 
  !########################################
  !Subroutine Array Arguments
  !########################################
  integer,dimension(3),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr
  integer,dimension(3,nzones),intent(out) :: astart !
  integer,dimension(3,nzones),intent(out) :: aend !
  !#############################################
  !local variables
  !############################################
  integer :: ii,index,jj
  integer,dimension(3) :: alrs,alre,Gend,Gstart,period
  character(len=*), parameter :: subname='fracture_periodic_zone'
  
! Start and end of Global region
  Gstart(1) = Glr%ns1 
  Gstart(2) = Glr%ns2
  Gstart(3) = Glr%ns3  
  Gend(1) = Glr%ns1 + Glr%d%n1
  Gend(2) = Glr%ns2 + Glr%d%n2
  Gend(3) = Glr%ns3 + Glr%d%n3

! Periodicity of the system
  period(1) = Glr%d%n1 + 1
  period(2) = Glr%d%n2 + 1
  period(3) = Glr%d%n3 + 1

! Start and end of local region
  alrs(1) = Llr%ns1
  alrs(2) = Llr%ns2
  alrs(3) = Llr%ns3
  alre(1) = Llr%ns1 + Llr%d%n1
  alre(2) = Llr%ns2 + Llr%d%n2
  alre(3) = Llr%ns3 + Llr%d%n3

!assign the first zone (necessarily without shift) and initialize the rest
  do ii=1,3
     astart(ii,:) = alrs(ii)
     aend(ii,:) = min(Gend(ii),alre(ii))
  end do

!assign the other zones
  index = 2
  do ii=1,3
     if(outofzone(ii) > 0) then    !Translation: X,Y,Z
        astart(ii,index) =  Gstart(ii)
        aend(ii,index) = modulo(alre(ii),period(ii))
        index = index + 1
     end if 
     do jj=ii+1,3
        if(outofzone(ii) > 0 .and. outofzone(jj) > 0) then  !Translation: X+Y,X+Z,Y+Z
           astart(ii,index) = Gstart(ii)
           astart(jj,index) = Gstart(jj)
           aend(ii,index) = modulo(alre(ii),period(ii))
           aend(jj,index) = modulo(alre(jj),period(jj))
           index = index + 1
        end if
     end do
  end do

  if(outofzone(1) > 0 .and. outofzone(2) > 0 .and. outofzone(3) > 0 ) then ! Translation: X+Y+Z
     astart(1,index) = Gstart(1)
     astart(2,index) = Gstart(2)
     astart(3,index) = Gstart(3)
     aend(1,index) = modulo(alre(1),period(1))
     aend(2,index) = modulo(alre(2),period(2))
     aend(3,index) = modulo(alre(3),period(3))
  end if

END SUBROUTINE fracture_periodic_zone


!##############################################################################################################################################
!!****f* BigDFT/get_overlap_region
!##############################################################################################################################################
!! FUNCTION Given two localization regions, A and B, this routine returns a localization region corresponding to the intersection of A & B. 
!!
!! SOURCE
!!
subroutine get_overlap_region_periodic(alr,blr,Glr,isovrlp,Llr,nlr,Olr)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer, intent(in) :: alr,blr              ! index of the two localization regions
  integer, intent(in) :: nlr                  ! number of localization regions
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  integer, intent(in) :: isovrlp              ! Number of overlap regions
  !########################################
  !Subroutine Array Arguments
  !########################################
  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr  ! Localization grid descriptors 
  type(locreg_descriptors),dimension(isovrlp),intent(out) :: Olr ! Overlap localization regions
  !#############################################
  !local variables
  !############################################
  integer :: axmin,axmax,aymin,aymax,azmin,azmax ! bounds of localization region A
  integer :: bxmin,bxmax,bymin,bymax,bzmin,bzmax ! bounds of localization region B
  integer :: isx,isy,isz,iex,iey,iez             ! bounds of the overlap region
  character(len=*), parameter :: subname='get_overlap_region_periodic'
  !# NEW
  integer :: ii,azones,bzones,i_stat,i_all,index
  integer :: izones,jzones
  integer,allocatable :: astart(:,:),aend(:,:),bstart(:,:),bend(:,:)
  logical :: go1,go2,go3

  azones = 1
  bzones = 1
! Calculate the number of regions to cut alr and blr
  do ii=1,3
     if(Llr(alr)%outofzone(ii) > 0) azones = azones * 2
     if(Llr(blr)%outofzone(ii) > 0) bzones = bzones * 2
  end do

!allocate astart and aend
  allocate(astart(3,azones),stat=i_stat)
  call memocc(i_stat,astart,'astart',subname)
  allocate(aend(3,azones),stat=i_stat)
  call memocc(i_stat,aend,'aend',subname)

!FRACTURE THE FIRST LOCALIZATION REGION
  call fracture_periodic_zone(azones,Glr,Llr(alr),Llr(alr)%outofzone(:),astart,aend)

!allocate bstart and bend
  allocate(bstart(3,bzones),stat=i_stat)
  call memocc(i_stat,bstart,'bstart',subname)
  allocate(bend(3,bzones),stat=i_stat)
  call memocc(i_stat,bend,'bend',subname)

!FRACTURE SECOND LOCREG
  call fracture_periodic_zone(bzones,Glr,Llr(blr),Llr(blr)%outofzone(:),bstart,bend)

! Now check the number of overlapping zones
  index = 0
  do izones=1,azones
    do jzones=1,bzones
      go1 = (bstart(1,jzones) .le. aend(1,izones) .and. bend(1,jzones) .ge. astart(1,izones)) 
      go2 = (bstart(2,jzones) .le. aend(2,izones) .and. bend(2,jzones) .ge. astart(2,izones)) 
      go3 = (bstart(3,jzones) .le. aend(3,izones) .and. bend(3,jzones) .ge. astart(3,izones)) 
      if(go1 .and. go2 .and. go3) then
        index = index + 1

! Now construct the Overlap localization region descriptor
! only if there is an overlap. The following only works
! when the previous test is successful. Note also that
! isx, isy and isz are necessarily in the Glr by construction
! of the Llrs, so don't need to test them.
         
        ! Determine the limits of the overlap region
        isx = max(astart(1,izones),bstart(1,jzones))
        isy = max(astart(2,izones),bstart(2,jzones))
        isz = max(astart(3,izones),bstart(3,jzones))

        iex = min(aend(1,izones),bend(1,jzones))
        iey = min(aend(2,izones),bend(2,jzones))
        iez = min(aend(3,izones),bend(3,jzones))

!       Checks to assign the geometric code of the overlap region (TO DO,could be interesting for Pascal?)
!       This could change the values of the bounds, so do it here
!       for now, in sandbox,put free boundary to all zones
        Olr(index)%geocode = 'F'  

!       Values for the starting point of the cube
        Olr(index)%ns1 = isx
        Olr(index)%ns2 = isy
        Olr(index)%ns3 = isz

!       Dimensions of the overlap region
        Olr(index)%d%n1 = iex - isx 
        Olr(index)%d%n2 = iey - isy 
        Olr(index)%d%n3 = iez - isz 
    
!       Dimensions of the fine grid inside the overlap region
        if (isx < iex) then
           Olr(index)%d%nfl1=max(isx,Glr%d%nfl1)-isx
           Olr(index)%d%nfu1=min(iex,Glr%d%nfu1)-isx
        else
           write(*,*)'Yet to be implemented (little effort?)'
           stop
        end if

        if (isy < iey) then
           Olr(index)%d%nfl2=max(isy,Glr%d%nfl2)-isy
           Olr(index)%d%nfu2=min(iey,Glr%d%nfu2)-isy
        else
           write(*,*)'Yet to be implemented (little effort?)'
           stop
        end if

        if (isz < iez) then
           Olr(index)%d%nfl3=max(isz,Glr%d%nfl3)-isz
           Olr(index)%d%nfu3=min(iez,Glr%d%nfu3)-isz
        else
           write(*,*)'Yet to be implemented (little effort?)'
           stop
        end if

!       Dimensions of the interpolating scaling function grid 
!       (geocode already taken into acount because it is simple)
        select case(Olr(index)%geocode)
        case('F')
          Olr(index)%d%n1i=2*Olr(index)%d%n1+31
          Olr(index)%d%n2i=2*Olr(index)%d%n2+31
          Olr(index)%d%n3i=2*Olr(index)%d%n3+31
        case('S')
          Olr(index)%d%n1i=2*Olr(index)%d%n1+2
          Olr(index)%d%n2i=2*Olr(index)%d%n2+31
          Olr(index)%d%n3i=2*Olr(index)%d%n3+2
        case('P')
          Olr(index)%d%n1i=2*Olr(index)%d%n1+2
          Olr(index)%d%n2i=2*Olr(index)%d%n2+2
          Olr(index)%d%n3i=2*Olr(index)%d%n3+2
        end select
 
!       Now define the wavefunction descriptors inside the overlap region
!       First calculate the number of points and segments for the region
!       Coarse part:
        call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
         Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
         Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c)
!       Fine part:
        call num_segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
         Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
         Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
         Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
         Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f)

!       Now allocate the wavefunction descriptors (keyg,keyv) following the needs
        call allocate_wfd(Olr(index)%wfd,subname)

!       At last, fill the wavefunction descriptors
!       Coarse part
        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
          Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
          Olr(index)%wfd%nseg_c,Olr(index)%wfd%nvctr_c,&
          Olr(index)%wfd%keyg(1,1),Olr(index)%wfd%keyv(1))
!       Fine part
        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
          Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
          Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
          Olr(index)%wfd%nseg_f,Olr(index)%wfd%nvctr_f,&
          Olr(index)%wfd%keyg(1,Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)),&
          Olr(index)%wfd%keyv(Olr(index)%wfd%nseg_c+min(1,Olr(index)%wfd%nseg_f)))

!       If the localisation region is isolated build also the bounds
        if (Olr(index)%geocode=='F') then
           call locreg_bounds(Olr(index)%d%n1,Olr(index)%d%n2,Olr(index)%d%n3,&
             Olr(index)%d%nfl1,Olr(index)%d%nfu1,Olr(index)%d%nfl2,Olr(index)%d%nfu2,&
             Olr(index)%d%nfl3,Olr(index)%d%nfu3,Olr(index)%wfd,Olr(index)%bounds)
     
        end if
     end if ! go1 .and. go2 .and. go3
   end do !jzones
 end do !izones

! Deallocation block
  i_all = -product(shape(astart))*kind(astart)
  deallocate(astart,stat=i_stat)
  call memocc(i_stat,i_all,'astart',subname)
  i_all = -product(shape(aend))*kind(aend)
  deallocate(aend,stat=i_stat)
  call memocc(i_stat,i_all,'aend',subname)
  i_all = -product(shape(bstart))*kind(bstart)
  deallocate(bstart,stat=i_stat)
  call memocc(i_stat,i_all,'bstart',subname)
  i_all = -product(shape(bend))*kind(bend)
  deallocate(bend,stat=i_stat)
  call memocc(i_stat,i_all,'bend',subname)

! Check on the number of zones
  if (index /= isovrlp) then
      write(*,*)&
          'ERROR: problem in get_overlap_region_periodic ',&
          'index:',index,'not equal to isovrlp:',isovrlp,&
          'The number of overlap descriptors constructed does not',&
          'correspond to the number of overlap regions.'
     stop
  end if

END SUBROUTINE get_overlap_region_periodic
!%***

subroutine assignToLocreg(iproc, natom, nlr, nspin, Localnorb, orbse)
  use module_base
  use module_types
  implicit none
  
  integer,intent(in):: nlr,iproc,nspin,natom
  integer,dimension(nlr),intent(in):: Localnorb
  type(orbitals_data),intent(inout):: orbse
  
  ! Local variables
  integer:: jproc, iiOrb, iorb, jorb, jat,i_stat
  character(len=*), parameter :: subname='assignToLocreg'

  allocate(orbse%inWhichLocreg(orbse%norbp),stat=i_stat)
  call memocc(i_stat,orbse%inWhichLocreg,'orbse%inWhichLocreg',subname)

  ! There are four counters:
  !   jproc: indicates which MPI process is handling the basis function which is being treated
  !   jat: counts the atom numbers
  !   jorb: counts the orbitals handled by a given process
  !   iiOrb: counts the number of orbitals for a given atom thas has already been assigned
  jproc=0
  jat=1
  jorb=0
  iiOrb=0

  do iorb=1,orbse%norb

      ! Switch to the next MPI process if the numbers of orbitals for a given
      ! MPI process is reached.
      if(jorb==orbse%norb_par(jproc)) then
          jproc=jproc+1
          jorb=0
      end if

      ! Switch to the next atom if the number of basis functions for this atom is reached.
      if(iiOrb==Localnorb(jat)) then
          jat=jat+1
          iiOrb=0
      end if
      if(jat > natom) then
        jat = 1
      end if
      jorb=jorb+1
      iiOrb=iiOrb+1
      if(iproc==jproc) orbse%inWhichLocreg(jorb)=jat
  end do

end subroutine assignToLocreg
