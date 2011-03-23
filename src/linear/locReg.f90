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
