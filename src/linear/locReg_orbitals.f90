
!============================================================================
!WARNING: assignToLocreg does not take into account the Kpts yet !!
!============================================================================
subroutine assignToLocreg(iproc,nproc,nspinor,nspin,atoms,orbs,Lzd)
  use module_base
  use module_types
  implicit none

  integer,intent(in):: iproc,nproc,nspin,nspinor
  type(atoms_data),intent(in) :: atoms 
  type(orbitals_data),intent(inout):: orbs
  type(local_zone_descriptors) :: Lzd
  ! Local variables
  integer :: jproc,iiOrb,iorb,jorb,jat,i_stat,orbsctot,orbsc,ispin
  integer :: iat,ind,i_all,noncoll,Lnorb,dimtot,ilr,npsidim,ierr
  character(len=*), parameter :: subname='assignToLocreg'
  integer, dimension(:), allocatable :: Localnorb
  integer, parameter :: nmax=6,lmax=3,noccmax=2,nelecmax=32
  integer, dimension(lmax+1) :: nmoments
  real(gp), dimension(noccmax,lmax) :: occup              !dummy variable


! in the non-collinear case the number of orbitals double
! nspin_ig*noncoll is always <= 2
  if (nspinor == 4) then
     noncoll=2
  else
     noncoll=1
  end if

  allocate(Localnorb(Lzd%nlr+ndebug),stat=i_stat)
  call memocc(i_stat,Localnorb,'Localnorb',subname)

! NOTES: WORKS ONLY BECAUSE Llr coincides with the atoms !!
! NOTES: K-Points??
  nmoments = 0
  do ilr = 1, Lzd%nlr
     call count_atomic_shells(lmax,noccmax,nelecmax,nspin,nspinor,atoms%aocc(1,ilr),occup,nmoments)
     Lnorb=(nmoments(1)+3*nmoments(2)+5*nmoments(3)+7*nmoments(4))
     Localnorb(ilr) = Lnorb
  end do

!!  already associated = 1 by default
!  allocate(orbs%inWhichLocregP(max(1,orbs%norb_par(iproc,0))),stat=i_stat)

! initialize inwhichlocreg
  orbs%inWhichLocreg = 0
  !orbs%inWhichLocregP = 0

  ind = 0
  jproc=0
  jat=1
  jorb=ind
  iiOrb=0
  do iorb=ind,orbs%norb

      ! Switch to the next MPI process if the numbers of orbitals for a given
      ! MPI process is reached.
      if(jorb==orbs%norb_par(jproc,0)) then
          jproc=jproc+1
          jorb=ind
          if (jproc==nproc) exit
      end if
      orbsc = 0
      ! Switch to the next atom if the number of basis functions for this atom is reached.
      if(iiOrb==(Localnorb(jat)-orbsc)*noncoll) then
          jat=jat+1
          iiOrb=0
      end if
      if(jat > atoms%nat) then
        jat = 1
      end if
      jorb=jorb+1
      iiOrb=iiOrb+1
      if(iproc==jproc .and. orbs%norb_par(jproc,0)> 0) then
         !orbs%inWhichLocregP(jorb)=jat
         orbs%inWhichLocreg(jorb+orbs%isorb)=jat
      end if
  end do
  call mpiallred(orbs%inWhichLocreg(1),orbs%norb,MPI_SUM,MPI_COMM_WORLD,ierr)


! Calculate the dimension of the total wavefunction
!!  dimtot = 0
!!  if(orbs%norbp > 0) then
!!     do iorb = 1,orbs%norbp
!!        ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
!!        npsidim = (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*nspinor
!!        dimtot = dimtot + npsidim
!!     end do
!!  else if (orbs%norbp == 0) then
!!       dimtot = orbs%npsidim
!!  end if
!!  Lzd%Lpsidimtot = dimtot

  i_all=-product(shape(Localnorb))*kind(Localnorb)
  deallocate(Localnorb,stat=i_stat)
  call memocc(i_stat,i_all,'Localnorb',subname)

end subroutine assignToLocreg

subroutine wavefunction_dimension(Lzd,orbs)
  use module_types
  implicit none
  type(local_zone_descriptors),intent(in) :: Lzd
  type(orbitals_data),intent(inout) :: orbs
  !local variables
  integer :: dimtot,iorb,ilr,npsidim

  dimtot = 0
  if(orbs%norbp > 0) then
     do iorb = 1,orbs%norbp
        ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
        npsidim = (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
        dimtot = dimtot + npsidim
     end do
  else if (orbs%norbp == 0) then
       dimtot = orbs%npsidim_orbs
  end if
  !Lzd%Lpsidimtot = dimtot
  orbs%npsidim_orbs=dimtot
end subroutine wavefunction_dimension


subroutine assignToLocreg2(iproc, nproc, norb, norb_par, natom, nlr, nspin, Localnorb, rxyz, inwhichlocreg)
  use module_base
  use module_types
  implicit none

  integer,intent(in):: nlr,iproc,nproc,nspin,natom,norb
  integer,dimension(nlr),intent(in):: Localnorb
  integer,dimension(0:nproc-1),intent(in):: norb_par
  !real(8),dimension(3,natom),intent(in):: rxyz
  real(8),dimension(3,nlr),intent(in):: rxyz
  integer,dimension(:),pointer, intent(out):: inwhichlocreg

  ! Local variables
  integer:: iat, jproc, iiOrb, iorb, jorb, jat, iiat, i_stat, i_all
  character(len=*), parameter :: subname='assignToLocreg'
  logical,dimension(:),allocatable:: covered
  real(8):: tt, dmin, minvalue, xmin, xmax, ymin, ymax, zmin, zmax
  integer:: iatxmin, iatxmax, iatymin, iatymax, iatzmin, iatzmax, idir
  real(8),dimension(3):: diff

!!!! NEW VERSION #################################################################
  !allocate(orbse%inWhichLocreg(orbse%norbp),stat=i_stat)
  allocate(inWhichLocreg(norb),stat=i_stat)
  call memocc(i_stat,inWhichLocreg,'inWhichLocreg',subname)
  inWhichLocreg=-1
  !allocate(orbse%inWhichLocregp(orbse%norbp),stat=i_stat)
  !call memocc(i_stat,orbse%inWhichLocregp,'orbse%inWhichLocregp',subname)
  allocate(covered(nlr), stat=i_stat)
  call memocc(i_stat, covered, 'covered', subname)


  ! Determine in which direction the system has its largest extent
  xmin=1.d100
  ymin=1.d100
  zmin=1.d100
  xmax=-1.d100
  ymax=-1.d100
  zmax=-1.d100
  do iat=1,nlr
  !write(*,'(a,2i8,3es16.7)') 'iproc, iat, rxyz(1,iat), rxyz(2,iat), rxyz(3,iat)', iproc, iat, rxyz(1,iat), rxyz(2,iat), rxyz(3,iat)
      if(rxyz(1,iat)<xmin) then
          xmin=rxyz(1,iat)
          iatxmin=iat
      end if
      if(rxyz(1,iat)>xmax) then
          xmax=rxyz(1,iat)
          iatxmax=iat
      end if
      if(rxyz(2,iat)<ymin) then
          ymin=rxyz(2,iat)
          iatymin=iat
      end if
      if(rxyz(2,iat)>ymax) then
          ymax=rxyz(2,iat)
          iatymax=iat
      end if
      if(rxyz(3,iat)<zmin) then
          zmin=rxyz(3,iat)
          iatzmin=iat
      end if
      if(rxyz(3,iat)>zmax) then
          zmax=rxyz(3,iat)
          iatzmax=iat
      end if
  end do

  diff(1)=xmax-xmin
  diff(2)=ymax-ymin
  diff(3)=zmax-zmin
  if(maxloc(diff,1)==1) then
      idir=1
      iiat=iatxmin
  else if(maxloc(diff,1)==2) then
      idir=2
      iiat=iatymin
  else if(maxloc(diff,1)==3) then
      idir=3
      iiat=iatzmin
  else
      stop 'ERROR: not possible to determine the maximal extent'
  end if
 

  !! Determine the atom with lowest z coordinate
  !zmin=1.d100
  !    do iat=1,natom
  !    if(rxyz(3,iat)<zmin) then
  !        zmin=rxyz(3,iat)
  !        iiat=iat
  !    end if
  !end do

  ! There are four counters:
  !   jproc: indicates which MPI process is handling the basis function which is being treated
  !   jat: counts the atom numbers
  !   jorb: counts the orbitals handled by a given process
  !   iiOrb: counts the number of orbitals for a given atom thas has already been assigned
  jproc=0
  !jat=iiat
  jat=1
  jorb=0
  iiOrb=0

  covered=.false.
  covered(iiat)=.true.
  inWhichLocreg(1)=iiat
  iiorb=1

  do iorb=2,norb

      ! Switch to the next MPI process if the numbers of orbitals for a given
      ! MPI process is reached.
      !if(jorb==norb_par(jproc,0)) then
      if(jorb==norb_par(jproc)) then
          jproc=jproc+1
          jorb=0
      end if

      ! Switch to the next atom if the number of basis functions for this atom is reached.
      !if(iiOrb==Localnorb(jat)) then
      !if(iproc==0) write(*,*) 'localnorb(iiat)',localnorb(iiat)
      if(iiOrb==Localnorb(iiat)) then
          iiOrb=0
          !jat=jat+1
          ! Determine the nearest atom which has not been covered yet.
          !covered(jat)=.true.
          dmin=1.d100
          minvalue=1.d100
          do iat=1,nlr
              !write(*,'(a,i8,a,l3)') 'iproc, iorb, minvalue, iiat, covered', iproc, ' covered(iat) ', covered(iat)
              !if(iproc==0 .and. nlr>12) write(*,'(a,2i6,l5,i7)') 'iorb, iat, covered(13), iiat', iorb, iat, covered(13), iiat
              if(covered(iat)) then
                  !!write(*,'(a,i8,a,i4)') 'iproc, iorb, minvalue, iiat, covered', iproc, 'cycles for iat=',iat
                  cycle
              end if
              tt = (rxyz(1,iat)-rxyz(1,jat))**2 + (rxyz(2,iat)-rxyz(2,jat))**2 + (rxyz(3,iat)-rxyz(3,jat))**2
              !if(tt<dmin) then
              if(rxyz(idir,iat)<minvalue) then
                  iiat=iat
                  dmin=tt
                  minvalue=rxyz(idir,iat)
              end if
          end do
          !jat=iiat
          jat=jat+1
          covered(iiat)=.true.
      end if
      if(jat > nlr) then
        jat = 1
      end if
      jorb=jorb+1
      iiOrb=iiOrb+1
      !if(iproc==jproc) orbse%inWhichLocregp(jorb)=jat
      !orbse%inWhichLocreg(iorb)=jat
      !if(iproc==0) write(*,'(a,2i8,es16.8,i8,20l3)') 'iproc, iorb, minvalue, iiat, covered', iproc, iorb, minvalue, iiat, covered
      inWhichLocreg(iorb)=iiat
  end do

  i_all=-product(shape(covered))*kind(covered)
  deallocate(covered,stat=i_stat)
  call memocc(i_stat,i_all,'covered',subname)

  !write(*,'(a,i3,3x,100i4)') 'iproc, orbse%inWhichLocreg', iproc, orbse%inWhichLocreg
  !write(*,'(a,i3,3x,100i4)') 'iproc, orbse%inWhichLocregp', iproc, orbse%inWhichLocregp



!!!! OLD VERSION #################################################################
!!  !allocate(orbse%inWhichLocreg(orbse%norbp),stat=i_stat)
!!  allocate(orbse%inWhichLocreg(orbse%norb),stat=i_stat)
!!  call memocc(i_stat,orbse%inWhichLocreg,'orbse%inWhichLocreg',subname)
!!  !allocate(orbse%inWhichLocregp(orbse%norbp),stat=i_stat)
!!  !call memocc(i_stat,orbse%inWhichLocregp,'orbse%inWhichLocregp',subname)
!!  allocate(covered(natom), stat=i_stat)
!!  call memocc(i_stat, covered, 'covered', subname)
!! 
!!
!!  ! Determine the atom with lowest z coordinate
!!  zmin=1.d100
!!      do iat=1,natom
!!      if(rxyz(3,iat)<zmin) then
!!          zmin=rxyz(3,iat)
!!          iiat=iat
!!      end if
!!  end do
!!
!!  ! There are four counters:
!!  !   jproc: indicates which MPI process is handling the basis function which is being treated
!!  !   jat: counts the atom numbers
!!  !   jorb: counts the orbitals handled by a given process
!!  !   iiOrb: counts the number of orbitals for a given atom thas has already been assigned
!!  jproc=0
!!  !jat=iiat
!!  jat=1
!!  jorb=0
!!  iiOrb=0
!!
!!  covered=.false.
!!
!!  do iorb=1,orbse%norb
!!
!!      ! Switch to the next MPI process if the numbers of orbitals for a given
!!      ! MPI process is reached.
!!      if(jorb==orbse%norb_par(jproc,0)) then
!!          jproc=jproc+1
!!          jorb=0
!!      end if
!!
!!      ! Switch to the next atom if the number of basis functions for this atom is reached.
!!      if(iiOrb==Localnorb(jat)) then
!!          iiOrb=0
!!          !jat=jat+1
!!          ! Determine the nearest atom which has not been covered yet.
!!          covered(jat)=.true.
!!          dmin=1.d100
!!          do iat=1,natom
!!              if(covered(iat)) cycle
!!              tt = (rxyz(1,iat)-rxyz(1,jat))**2 + (rxyz(2,iat)-rxyz(2,jat))**2 + (rxyz(3,iat)-rxyz(3,jat))**2
!!              if(tt<dmin) then
!!                  iiat=iat
!!                  dmin=tt
!!              end if
!!          end do
!!          !jat=iiat
!!          jat=jat+1
!!      end if
!!      if(jat > natom) then
!!        jat = 1
!!      end if
!!      jorb=jorb+1
!!      iiOrb=iiOrb+1
!!      !if(iproc==jproc) orbse%inWhichLocregp(jorb)=jat
!!      orbse%inWhichLocreg(iorb)=jat
!!  end do
!!
!!  i_all=-product(shape(covered))*kind(covered)
!!  deallocate(covered,stat=i_stat)
!!  call memocc(i_stat,i_all,'covered',subname)
!!
!!  !write(*,'(a,i3,3x,100i4)') 'iproc, orbse%inWhichLocreg', iproc, orbse%inWhichLocreg
!!  !write(*,'(a,i3,3x,100i4)') 'iproc, orbse%inWhichLocregp', iproc, orbse%inWhichLocregp

end subroutine assignToLocreg2

subroutine orbital_ordering_by_locreg(ilr,iproc,isorb_par,nproc,norb,norbp,nlr,norb_par,inwhichlocreg,orbs)
  use module_base
  use module_types

  implicit none
  integer, intent(in) :: ilr                                     ! number of this locreg
  integer, intent(in) :: iproc                                   ! number of this process
  integer, intent(in) :: nproc                                   ! number of processes
  integer, intent(in) :: norbp                                   ! number of orbitals on this processor
  integer, intent(in) :: nlr                                     ! number of localization regions
  integer, intent(in) :: norb                                    ! total number of orbitals
  integer, dimension(nproc), intent(in) :: isorb_par             ! reference orbital for the processes
  integer, dimension(nproc), intent(in) :: norb_par              ! number of orbitals on this processor
  integer, dimension(norbp,nproc),intent(in) :: inwhichlocreg    ! mapping of the orbitals on this processor to the locregs
  type(orbitals_data), intent(inout) :: orbs                     ! Local orbitals_data types
  ! local variables
  character(len=*), parameter :: subname='orbital_ordering_by_locreg'
  integer :: jproc,ii
  integer :: iorb,i_stat

  ! Allocate inWhichLocreg
  allocate(orbs%inWhichLocreg(norb),stat=i_stat)
  call memocc(i_stat,orbs%inWhichLocreg,'orbs%inWhichLocreg',subname)

  iorb = 0
     do jproc=1,nproc
        do ii=1,norb_par(jproc)
           if (inwhichlocreg(ii,jproc) .ne. ilr) cycle
           iorb = iorb + 1
           orbs%inWhichLocreg(iorb) = ii+isorb_par(jproc)
        end do
     end do

  if(iorb .ne. norb) then
    write(*,'(a,i4,a,i4)') 'Error in orbital_ordering_by_locreg: iorb=',iorb,'is not equal to norb=',norb
  end if

end subroutine orbital_ordering_by_locreg


