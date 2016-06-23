!> @file
!! Localization Region to orbitals
!! @author
!!    Copyright (C) 2011-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> assignToLocreg does not take into account the Kpts yet !!
!! @warning assignToLocreg does not take into account the Kpts yet !!
subroutine assignToLocreg(iproc,nproc,nspinor,nspin,atoms,orbs,Lzd)
  use module_base
  use module_types
  !use ao_inguess, only: count_atomic_shells, ao_nspin_ig
  implicit none

  integer,intent(in):: iproc,nproc,nspin,nspinor
  type(atoms_data),intent(in) :: atoms
  type(orbitals_data),intent(inout):: orbs
  type(local_zone_descriptors) :: Lzd
  ! Local variables
  integer :: jproc,iiOrb,iorb,jorb,jat,orbsc!,ispin
  integer :: ind,noncoll,ilr!,ierr,dimtot,iat,npsidim,Lnorb
  character(len=*), parameter :: subname='assignToLocreg'
  integer, dimension(:), allocatable :: Localnorb
  !integer, parameter :: lmax=3,noccmax=2,nelecmax=32
  !integer, dimension(lmax+1) :: nmoments
  !real(gp), dimension(noccmax,lmax+1) :: occup              !dummy variable


! in the non-collinear case the number of orbitals double
! nspin_ig*noncoll is always <= 2
  if (nspinor == 4) then
     noncoll=2
  else
     noncoll=1
  end if

  Localnorb = f_malloc(Lzd%nlr,id='Localnorb')

! NOTES: WORKS ONLY BECAUSE Llr coincides with the atoms !!
! NOTES: K-Points??
  !nmoments = 0
  do ilr = 1, Lzd%nlr
     !call count_atomic_shells(ao_nspin_ig(nspin,nspinor=nspinor),&
     !     atoms%aoig(ilr)%aocc,occup,nmoments)
     !Lnorb=(nmoments(1)+3*nmoments(2)+5*nmoments(3)+7*nmoments(4))
     Localnorb(ilr) =atoms%aoig(ilr)%nao! Lnorb
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
      if(jat > atoms%astruct%nat) then
        jat = 1
      end if
      jorb=jorb+1
      iiOrb=iiOrb+1
      if(iproc==jproc .and. orbs%norb_par(jproc,0)> 0) then
         !orbs%inWhichLocregP(jorb)=jat
         orbs%inWhichLocreg(jorb+orbs%isorb)=jat
      end if
  end do

  if (nproc > 1) then
     call mpiallred(orbs%inWhichLocreg(1),orbs%norb,MPI_SUM,comm=bigdft_mpi%mpi_comm)
  end if

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

  call f_free(Localnorb)

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


!!subroutine assignToLocreg2(iproc, nproc, norb, norbu, norb_par, natom, nlr, nspin, Localnorb, spinsgn, rxyz, inwhichlocreg)
!!  use module_base
!!  use sort, only: QsortC
!!  implicit none
!!
!!  integer,intent(in):: nlr,iproc,nproc,nspin,natom,norb,norbu
!!  integer,dimension(nlr),intent(in):: Localnorb
!!  real(kind=8),dimension(norb),intent(in):: spinsgn
!!  integer,dimension(0:nproc-1),intent(in):: norb_par
!!  real(8),dimension(3,nlr),intent(in):: rxyz
!!  integer,dimension(:),pointer, intent(out):: inwhichlocreg
!!
!!  ! Local variables
!!  integer:: iat, jproc, iiOrb, iorb, jorb, jat, iiat, i_stat, i_all, ispin, iispin, istart, iend, ilr, jjat
!!  character(len=*), parameter :: subname='assignToLocreg'
!!  logical,dimension(:),allocatable:: covered
!!  integer,dimension(:),allocatable :: irxyz_ordered
!!  real(kind=8), parameter :: tol=1.0d-6 
!!  real(8):: tt, dmin, minvalue, xmin, xmax, ymin, ymax, zmin, zmax
!!  integer:: iatxmin, iatxmax, iatymin, iatymax, iatzmin, iatzmax, idir, nthread, ithread
!!  real(8),dimension(3):: diff
!!  real(kind=8),dimension(:),allocatable :: rxyz_dir, minvalue_thread
!!  integer,dimension(:),allocatable :: iiat_thread
!!  !$ integer :: omp_get_thread_num, omp_get_max_threads
!!
!!  call f_routine(id='assignToLocreg2')
!!
!!!!!! NEW VERSION #################################################################
!!  !allocate(orbse%inWhichLocreg(orbse%norbp),stat=i_stat)
!!  inWhichLocreg = f_malloc_ptr(norb,id='inWhichLocreg')
!!  do iorb=1,norb
!!      inWhichLocreg(iorb) = -1
!!  end do
!!  !allocate(orbse%inWhichLocregp(orbse%norbp),stat=i_stat)
!!  !call memocc(i_stat,orbse%inWhichLocregp,'orbse%inWhichLocregp',subname)
!!  covered = f_malloc(nlr,id='covered')
!!
!!  nthread = 1
!!  !$ nthread = omp_get_max_threads()
!!  iiat_thread = f_malloc(0.to.nthread-1,id='iiat_thread')
!!  minvalue_thread = f_malloc(0.to.nthread-1,id='minvalue_thread')
!!
!!  spin_loop: do ispin=1,nspin
!!
!!      if (nlr==natom) then
!!          ! case for onwhichatom, i.e. distributing the orbitals to the atoms
!!          istart=1
!!          iend=nlr
!!      else
!!          if (ispin==1) then
!!              ! case for inwhichlocreg, i.e. distributing the orbitals to the locregs
!!              istart=1
!!              iend=norbu
!!          else
!!              istart=norbu+1
!!              iend=norb
!!          end if
!!      end if
!!
!!      ! Determine in which direction the system has its largest extent
!!      xmin=1.d100
!!      ymin=1.d100
!!      zmin=1.d100
!!      xmax=-1.d100
!!      ymax=-1.d100
!!      zmax=-1.d100
!!      !do iat=1,nlr
!!      do iat=istart,iend
!!          !!!SM: mixing nlr with orbs.. not ideal
!!          !!if (spinsgn(iat)>0.d0) then
!!          !!    iispin=1
!!          !!else
!!          !!    iispin=2
!!          !!end if
!!          !!if (ispin/=iispin) cycle
!!      !write(*,'(a,2i8,3es16.7)') 'iproc, iat, rxyz(1,iat), rxyz(2,iat), rxyz(3,iat)', iproc, iat, rxyz(1,iat), rxyz(2,iat), rxyz(3,iat)
!!          if(rxyz(1,iat)<xmin) then
!!              xmin=rxyz(1,iat)
!!              iatxmin=iat
!!          end if
!!          if(rxyz(1,iat)>xmax) then
!!              xmax=rxyz(1,iat)
!!              iatxmax=iat
!!          end if
!!          if(rxyz(2,iat)<ymin) then
!!              ymin=rxyz(2,iat)
!!              iatymin=iat
!!          end if
!!          if(rxyz(2,iat)>ymax) then
!!              ymax=rxyz(2,iat)
!!              iatymax=iat
!!          end if
!!          if(rxyz(3,iat)<zmin) then
!!              zmin=rxyz(3,iat)
!!              iatzmin=iat
!!          end if
!!          if(rxyz(3,iat)>zmax) then
!!              zmax=rxyz(3,iat)
!!              iatzmax=iat
!!          end if
!!      end do
!!
!!      diff(1)=xmax-xmin
!!      diff(2)=ymax-ymin
!!      diff(3)=zmax-zmin
!!      !First 4 ifs control if directions the same length to disambiguate (was random before)
!!      !else, just choose the biggest
!!      if(abs(diff(1)-diff(2)) < tol .and. diff(1) > diff(3)) then
!!        idir=1
!!        iiat=iatxmin
!!      else if(abs(diff(1)-diff(3)) < tol .and. diff(1) > diff(2)) then
!!        idir=1
!!        iiat=iatxmin
!!      else if(abs(diff(2)-diff(3)) < tol .and. diff(2) > diff(1)) then
!!        idir=2
!!        iiat=iatymin
!!      else if(abs(diff(1)-diff(3)) < tol .and. abs(diff(2)-diff(3)) < tol) then
!!        idir=1
!!        iiat=iatxmin
!!      else
!!         if(maxloc(diff,1)==1) then
!!             idir=1
!!             iiat=iatxmin
!!         else if(maxloc(diff,1)==2) then
!!             idir=2
!!             iiat=iatymin
!!         else if(maxloc(diff,1)==3) then
!!             idir=3
!!             iiat=iatzmin
!!         else
!!             stop 'ERROR: not possible to determine the maximal extent'
!!         end if
!!      end if
!!
!!      !! Determine the atom with lowest z coordinate
!!      !zmin=1.d100
!!      !    do iat=1,natom
!!      !    if(rxyz(3,iat)<zmin) then
!!      !        zmin=rxyz(3,iat)
!!      !        iiat=iat
!!      !    end if
!!      !end do
!!
!!      ! There are four counters:
!!      !   jproc: indicates which MPI process is handling the basis function which is being treated
!!      !   jat: counts the atom numbers
!!      !   jorb: counts the orbitals handled by a given process
!!      !   iiOrb: counts the number of orbitals for a given atom thas has already been assigned
!!      jproc=0
!!      !jat=iiat
!!      jat=1
!!      jorb=0
!!      iiOrb=0
!!
!!      !covered=.false.
!!      do ilr=1,nlr
!!          covered(ilr) = .false.
!!      end do
!!      covered(iiat)=.true.
!!      if (ispin==1) then
!!          inWhichLocreg(1)=iiat
!!      else
!!          inWhichLocreg(norbu+1)=iiat
!!      end if
!!      iiorb=1
!!
!!      !!! Lookup array which orders the atoms with respect to the direction idir
!!      !!rxyz_dir = f_malloc(nlr,id='rxyz_dir')
!!      !!irxyz_ordered = f_malloc(nlr,id='irxyz_ordered')
!!      !!do ilr=1,nlr
!!      !!    rxyz_dir(ilr) = rxyz(idir,ilr)
!!      !!    irxyz_ordered(ilr) = ilr
!!      !!end do
!!      !!do ilr=1,nlr
!!      !!    write(*,*) 'B: rxyz_dir(ilr), irxyz_ordered(ilr)', rxyz_dir(ilr), irxyz_ordered(ilr)
!!      !!end do
!!      !!call QsortC(rxyz_dir, irxyz_ordered)
!!      !!do ilr=1,nlr
!!      !!    write(*,*) 'A: rxyz_dir(ilr), irxyz_ordered(ilr)', rxyz_dir(ilr), irxyz_ordered(ilr)
!!      !!end do
!!      !!write(*,*) 'inWhichLocreg(1), rxyz(idir,inWhichLocreg(1)', inWhichLocreg(1), rxyz(idir,inWhichLocreg(1))
!!
!!
!!      ! SM: This loop has a quadratic scaling and probably becomes expensive for large systems
!!      jjat = 1
!!      do iorb=1,norb
!!
!!          if (iorb==1 .or. iorb==norbu+1) cycle !these values have already been assigned
!!
!!          !SM: mixing nlr with orbs.. not ideal
!!          if (spinsgn(iorb)>0.d0) then
!!              iispin=1
!!          else
!!              iispin=2
!!          end if
!!          if (ispin/=iispin) cycle
!!
!!          if(iiOrb==Localnorb(iiat)) then
!!              jjat = jjat + 1
!!              iiOrb=0
!!              minvalue_thread(:)=huge(minvalue_thread(:))
!!              ithread = 0
!!              !$omp parallel default(none) shared(istart,iend,covered,rxyz,idir,minvalue_thread,iiat_thread) &
!!              !$omp private(iat) firstprivate(ithread)
!!              !$ ithread = omp_get_thread_num()
!!              !$omp do 
!!              do iat=istart,iend
!!                  if(.not.covered(iat)) then
!!                      if(rxyz(idir,iat)<minvalue_thread(ithread)) then
!!                          !iiat=iat
!!                          !minvalue=rxyz(idir,iat)
!!                          iiat_thread(ithread)=iat
!!                          minvalue_thread(ithread)=rxyz(idir,iat)
!!                      end if
!!                  end if
!!$           call f_err_throw('ERROR: not possible to determine the maximal extent',&
!!$                & err_name='BIGDFT_RUN_TIME')
!!              end do
!!              !$omp end do
!!              !$omp end parallel
!!              minvalue = huge(minvalue)
!!              do ithread=0,nthread-1
!!                  if (minvalue_thread(ithread)<minvalue) then
!!                      iiat = iiat_thread(ithread)
!!                      minvalue = minvalue_thread(ithread)
!!                  end if
!!              end do
!!              covered(iiat)=.true.
!!              !!iiat = irxyz_ordered(jjat)
!!              !write(*,*) 'iorb, iiorb, iiat', iorb, iiorb, iiat, jjat
!!          end if
!!          iiOrb=iiOrb+1
!!          inWhichLocreg(iorb)=iiat
!!          !!write(*,*) 'iorb, inWhichLocreg(iorb), rxyz(idir,iiat)', iorb, inWhichLocreg(iorb), rxyz(idir,iiat)
!!      end do
!!
!!  end do spin_loop
!!
!!  call f_free(iiat_thread)
!!  call f_free(minvalue_thread)
!!
!!  call f_free(covered)
!!
!!  call f_release_routine()
!!
!!end subroutine assignToLocreg2
