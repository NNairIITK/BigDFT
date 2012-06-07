!###########################################################################
! NOT USED ANYMORE
! DOES NOT WORK WITH MODERN MODIFICATIONS
!##########################################################################
!!!subroutine determine_Lorbs(iproc,nproc,at,Lzd,norbsc_arr,nspin)
!!!  use module_base
!!!  use module_types
!!!  implicit none
!!!  integer, intent(in) :: iproc,nproc
!!!  integer, intent(in) :: nspin
!!!  type(atoms_data), intent(in) :: at
!!!  type(local_zone_descriptors), intent(inout) :: Lzd
!!!  integer, dimension(at%natsc+1,nspin), intent(in) :: norbsc_arr
!!!  ! local variables
!!!  integer :: ilr
!!!  integer :: Lnorb
!!!  integer :: nspin_ig
!!!  integer :: noncoll
!!!  integer :: npsidim
!!!  integer :: dimtot
!!!  integer :: i_stat,i_all,ierr
!!!  integer, parameter :: nmax=6,lmax=3,noccmax=2,nelecmax=32
!!!  integer, dimension(lmax+1) :: nmoments
!!!  integer, dimension(Lzd%nlr) :: norbsc
!!!  integer, dimension(:), allocatable :: Localnorb
!!!  integer, dimension(:,:),allocatable :: iwl
!!!  real(gp), dimension(noccmax,lmax) :: occup              !dummy variable
!!!  character(len=*), parameter :: subname='determine_Lorbs'
!!!
!!!!spin for orbitals
!!!  if (nspin == 4) then
!!!     nspin_ig=1
!!!  else
!!!     nspin_ig=nspin
!!!  end if
!!!
!!!! in the non-collinear case the number of orbitals double
!!!! nspin_ig*noncoll is always <= 2
!!!  if (Lzd%orbs%nspinor == 4) then
!!!     noncoll=2
!!!  else
!!!     noncoll=1
!!!  end if
!!!
!!!  ! allocate statements
!!!
!!!  allocate(Localnorb(Lzd%nlr+ndebug),stat=i_stat)
!!!  call memocc(i_stat,Localnorb,'Localnorb',subname)
!!!  allocate(Lzd%Lorbs(Lzd%nlr+ndebug),stat=i_stat) !cannot call memocc because they are types
!!!                                                  ! it seems to associate all the pointers.
!!!
!!!! Calculate the dimension of the total wavefunction
!!!! NOTES: WORKS ONLY BECAUSE Llr coincides with the atoms !!
!!!! NOTES: K-Points??
!!!  dimtot = 0
!!!  do ilr = 1, Lzd%nlr
!!!     call count_atomic_shells(lmax,noccmax,nelecmax,nspin_ig,Lzd%orbs%nspinor,at%aocc(1,ilr),occup,nmoments)
!!!     Lnorb=(nmoments(1)+3*nmoments(2)+5*nmoments(3)+7*nmoments(4))
!!!     Localnorb(ilr) = Lnorb
!!!     npsidim = (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*Lnorb*Lzd%orbs%nspinor*nspin_ig*noncoll
!!!     dimtot = dimtot + npsidim
!!!     Lzd%Lorbs(ilr)%npsidim = npsidim
!!!  end do
!!!
!!!! Right total dimension of the wavefunctions
!!!  Lzd%Lpsidimtot = dimtot
!!!
!!!  allocate(iwl(Lzd%orbs%norbp,nproc),stat=i_stat)
!!!  call memocc(i_stat,iwl,'iwl',subname)
!!!  call razero(Lzd%orbs%norbp*nproc,iwl)
!!!
!!!! Determine inwhichlocreg (removed semicore treatement... to DEBUG)
!!!!  call assignToLocreg(iproc,nproc,at%nat,at%natsc,Lzd%nlr,noncoll,nspin_ig,Localnorb,Lzd%orbs,norbsc_arr,at%iasctype,norbsc,iwl)
!!!
!!!  call mpiallred(Lzd%orbs%inWhichLocreg(1), Lzd%orbs%norb, mpi_sum, mpi_comm_world, ierr)
!!!
!!!! DEBUG for inWhichLocreg(ilr)
!!!!  print *,'at%iasctype:',at%iasctype,Lzd%orbs%norb
!!!!  do ilr=1,Lzd%nlr
!!!!    print *,'ilr,localnorb:',ilr,Localnorb(ilr)
!!!!  end do
!!!!  do ilr=1,Lzd%orbs%norbp
!!!!    write(*,*) 'iorb, iwl', ilr, Lzd%orbs%inWhichLocreg(ilr),Lzd%orbs%occup(ilr)
!!!!  end do
!!!!  print *,'Global region:',iproc
!!!!  print *,'norb,norbp,norbu,norbd,nspin,nspinor,isorb,npsidim,nkpts,nkptsp,iskpts,norbsc,efermi',&
!!!!          Lzd%orbs%norb,Lzd%orbs%norbp,Lzd%orbs%norbu,Lzd%orbs%norbd,Lzd%orbs%nspin,Lzd%orbs%nspinor,&
!!!!          Lzd%orbs%isorb,Lzd%orbs%npsidim,Lzd%orbs%nkpts,Lzd%orbs%nkptsp,Lzd%orbs%iskpts,Lzd%orbs%norbsc,Lzd%orbs%efermi
!!!!  print *,'Lzd%orbs%norb_par',Lzd%orbs%norb_par
!!!!  print *,'Lzd%orbs%iokpt',Lzd%orbs%iokpt
!!!!  print *,'Lzd%orbs%ikptproc',Lzd%orbs%ikptproc
!!!!  print *,'inwhichlocreg',Lzd%orbs%inwhichlocreg
!!!!  print *,'inwhichlocregP',Lzd%orbs%inwhichlocregP
!!!!  print *,'onwhichMPI',Lzd%orbs%onwhichMPI
!!!!  print *,'Lzd%orbs%isorb_par',Lzd%orbs%isorb_par
!!!!  print *,'Lzd%orbs%occup',Lzd%orbs%occup
!!!!  print *,'Lzd%orbs%spinsgn',Lzd%orbs%spinsgn
!!!!  print *,'Lzd%orbs%kwgts',Lzd%orbs%kwgts
!!!!  print *,'Lzd%orbs%kpts',Lzd%orbs%kpts
!!!! END DEBUG
!!!
!!!  do ilr=1,Lzd%nlr
!!!     ! Set our new variable
!!!     Lzd%Lorbs(ilr)%norbsc = norbsc(ilr)
!!!     ! sets almost all of Lorbs(ilr) components
!!!     call linear_orbitals_descriptors(ilr,Lzd%nlr,iwl,iproc,nproc,Lzd%orbs,nspin_ig*noncoll*Localnorb(ilr),noncoll*Localnorb(ilr),&
!!!        (nspin_ig-1)*noncoll*Localnorb(ilr),nspin_ig,Lzd%orbs%nspinor,Lzd%orbs%nkpts,Lzd%orbs%kpts,&
!!!        Lzd%orbs%kwgts,Lzd%Lorbs(ilr))
!!!
!!!!DEBUG
!!!!    print *,'ilr:',ilr,iproc
!!!!    print *,'norb,norbp,norbu,norbd,nspin,nspinor:',Lzd%Lorbs(ilr)%norb,Lzd%Lorbs(ilr)%norbp,Lzd%Lorbs(ilr)%norbu,&
!!!!             Lzd%Lorbs(ilr)%norbd,Lzd%Lorbs(ilr)%nspin,Lzd%Lorbs(ilr)%nspinor
!!!!    print *,'isorb,nkpts,nkptsp,iskpts:',Lzd%Lorbs(ilr)%isorb,Lzd%Lorbs(ilr)%nkpts,Lzd%Lorbs(ilr)%nkptsp,Lzd%Lorbs(ilr)%iskpts
!!!!    print *,'norb_par',Lzd%Lorbs(ilr)%norb_par
!!!!    print *,'iokpt',Lzd%Lorbs(ilr)%iokpt
!!!!    print *,'ikptproc',Lzd%Lorbs(ilr)%ikptproc
!!!!    print *,'isorb_par',Lzd%Lorbs(ilr)%isorb_par
!!!!    print *,'occup',Lzd%Lorbs(ilr)%occup
!!!!    print *,'spinsgn',Lzd%Lorbs(ilr)%spinsgn
!!!!    print *,'kwgts',Lzd%Lorbs(ilr)%kwgts
!!!!    print *,'kpts',Lzd%Lorbs(ilr)%kpts
!!!!     if (iproc ==0) then
!!!!        write(*,'(a,i0)')'For locreg:',ilr
!!!!        write(*,'(1x,a,i0,a)')'Distributing ',nspin_ig*noncoll*Localnorb(ilr),' Atomic Input Orbitals'
!!!!        if (norbsc(ilr) /=0)   write(*,'(1x,a,i0,a)')'  of which ',norbsc(ilr),&
!!!!             ' are semicore orbitals'
!!!!     end if
!!!!END DEBUG
!!!
!!!  end do
!!!
!!!  !Deallocations
!!!!  i_all=-product(shape(Localnorb))*kind(Localnorb)
!!!!  deallocate(Localnorb,stat=i_stat)
!!!!  call memocc(i_stat,i_all,'Localnorb',subname)
!!!!  i_all=-product(shape(iwl))*kind(iwl)
!!!!  deallocate(iwl,stat=i_stat)
!!!!  call memocc(i_stat,i_all,'iwl',subname)
!!!
!!!
!!!end subroutine determine_Lorbs

!> Define the descriptors of the orbitals from a given norb
!! It uses the cubic strategy for partitioning the orbitals
subroutine linear_orbitals_descriptors(ilr,nlr,iwl,iproc,nproc,Gorbs,norb,norbu,norbd,nspin,nspinor,nkpt,kpt,wkpt,orbs)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nlr
  integer, intent(in) :: iproc,nproc,norb,norbu,norbd,nkpt,nspin
  integer, intent(in) :: nspinor,ilr
  type(orbitals_data), intent(in) :: Gorbs
  type(orbitals_data), intent(out) :: orbs
  real(gp), dimension(nkpt), intent(in) :: wkpt
  real(gp), dimension(3,nkpt), intent(in) :: kpt
  integer, dimension(Gorbs%norbp,nproc), intent(in) :: iwl
  !local variables
  character(len=*), parameter :: subname='linear_orbitals_descriptors'
  integer :: iorb,jproc,norb_tot,ikpt,i_stat,jorb,ierr,i_all,iiorb,ii
  logical, dimension(:), allocatable :: GPU_for_orbs
  integer, dimension(:), allocatable :: mykpts
  integer, dimension(:,:), allocatable :: norb_par !(with k-pts)

  allocate(orbs%norb_par(0:nproc-1+ndebug,0:orbs%nkpts),stat=i_stat)
  call memocc(i_stat,orbs%norb_par,'orbs%norb_par',subname)

  !assign the value of the k-points
  orbs%nkpts=nkpt
  !allocate vectors related to k-points
  allocate(orbs%kpts(3,orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%kpts,'orbs%kpts',subname)
  allocate(orbs%kwgts(orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%kwgts,'orbs%kwgts',subname)
  orbs%kpts(:,1:nkpt) = kpt(:,:)
  orbs%kwgts(1:nkpt) = wkpt(:)

  ! Change the wavefunctions to complex if k-points are used (except gamma).
  orbs%nspinor=nspinor
  if (nspinor == 1) then
     if (maxval(abs(orbs%kpts)) > 0._gp) orbs%nspinor=2
     !nspinor=2 !fake, used for testing with gamma
  end if
  orbs%nspin = nspin

  !initialise the array
  do jproc=0,nproc-1
     do ikpt = 0, orbs%nkpts
        orbs%norb_par(jproc,ikpt)=0 !size 0 nproc-1
     end do
  end do

  !create an array which indicate which processor has a GPU associated 
  !from the viewpoint of the BLAS routines
  if (.not. GPUshare) then
     allocate(GPU_for_orbs(0:nproc-1+ndebug),stat=i_stat)
     call memocc(i_stat,GPU_for_orbs,'GPU_for_orbs',subname)

     if (nproc > 1) then
        call MPI_ALLGATHER(GPUconv,1,MPI_LOGICAL,GPU_for_orbs(0),1,MPI_LOGICAL,&
             MPI_COMM_WORLD,ierr)
     else
        GPU_for_orbs(0)=GPUconv
     end if

     i_all=-product(shape(GPU_for_orbs))*kind(GPU_for_orbs)
     deallocate(GPU_for_orbs,stat=i_stat)
     call memocc(i_stat,i_all,'GPU_for_orbs',subname)
  end if

  call parallel_repartition_with_kpoints(nproc,Gorbs%nkpts,Gorbs%norb,orbs%norb_par)

  !check the distribution
  norb_tot=0
  do jproc=0,iproc-1
     norb_tot=norb_tot+Gorbs%norb_par(jproc,0)
  end do
  !reference orbital for process
  orbs%isorb=norb_tot
  do jproc=iproc,nproc-1
     norb_tot=norb_tot+Gorbs%norb_par(jproc,0)
  end do

  if(norb_tot /= Gorbs%norb*Gorbs%nkpts) then
     write(*,*)'ERROR: partition of orbitals incorrect'
     write(*,*)orbs%norb_par(:,:),Gorbs%norb*Gorbs%nkpts
     stop
  end if


  !calculate the k-points related quantities
  allocate(norb_par(0:nproc-1,orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,norb_par,'norb_par',subname)
  allocate(mykpts(orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,mykpts,'mykpts',subname)

  call parallel_repartition_per_kpoints(iproc,nproc,Gorbs%nkpts,Gorbs%norb,orbs%norb_par,&
       orbs%nkptsp,mykpts,norb_par)
  if (orbs%norb_par(iproc,0) >0) then
     orbs%iskpts=mykpts(1)-1
  else
     orbs%iskpts=0
  end if


  !allocate(orbs%ikptsp(orbs%nkptsp+ndebug),stat=i_stat)
  !call memocc(i_stat,orbs%ikptsp,'orbs%ikptsp',subname)
  !orbs%ikptsp(1:orbs%nkptsp)=mykpts(1:orbs%nkptsp)


  i_all=-product(shape(norb_par))*kind(norb_par)
  deallocate(norb_par,stat=i_stat)
  call memocc(i_stat,i_all,'norb_par',subname)
  i_all=-product(shape(mykpts))*kind(mykpts)
  deallocate(mykpts,stat=i_stat)
  call memocc(i_stat,i_all,'mykpts',subname)

  !assign the values of the orbitals data
  orbs%norb=norb
  orbs%norbp=orbs%norb_par(iproc,0)
  orbs%norbu=norbu
  orbs%norbd=norbd

  allocate(orbs%iokpt(orbs%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%iokpt,'orbs%iokpt',subname)

  !assign the k-point to the given orbital, counting one orbital after each other
  jorb=0
  do ikpt=1,Gorbs%nkpts
     do iorb=1,Gorbs%norb
        jorb=jorb+1 !this runs over norb*nkpts values
        if (jorb > orbs%isorb .and. jorb <= orbs%isorb+orbs%norbp) then
           orbs%iokpt(jorb-orbs%isorb)=ikpt
        end if
     end do
  end do

  ! allocate inwhichlocreg 
  ! WARNING: assigntoLocreg does not work with Kpts yet!!
  allocate(orbs%inwhichlocreg(orbs%norbp*orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%inwhichlocreg,'orbs%inwhichlocreg',subname)


  ! Define inwhichlocreg for the locregs
  ! this gives the mapping of the local orbitals to the Global one,
  ! because orbitals are reordered by locreg
  call orbital_ordering_by_locreg(ilr,iproc,Gorbs%isorb_par,nproc,orbs%norb,Gorbs%norbp,nlr,&
       Gorbs%norb_par,iwl,orbs)
  !allocate occupation number and spinsign
  !fill them in normal way
  allocate(orbs%occup(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%occup,'orbs%occup',subname)
  allocate(orbs%spinsgn(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%spinsgn,'orbs%spinsgn',subname)
  do iorb=1,norb*orbs%nkpts
     orbs%occup(iorb)= Gorbs%occup(orbs%inWhichLocreg(iorb))
     orbs%spinsgn(iorb) = Gorbs%spinsgn(orbs%inWhichLocreg(iorb))
  end do

  !put a default value for the fermi energy
  orbs%efermi = UNINITIALIZED(orbs%efermi)

  !allocate the array which assign the k-point to processor in transposed version
  allocate(orbs%ikptproc(orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%ikptproc,'orbs%ikptproc',subname)

  ! For right now, the K-points are the same
  do ii=1,orbs%nkpts
     orbs%ikptproc(ii) = Gorbs%ikptproc(ii)
  end do

  ! Define to new arrays:
  ! - orbs%isorb_par is the same as orbs%isorb, but every process also knows
  !   the reference orbital of each other process.
  ! - orbs%onWhichMPI indicates on which MPI process a given orbital
  !   is located.
  allocate(orbs%isorb_par(0:nproc-1), stat=i_stat)
  call memocc(i_stat, orbs%isorb_par, 'orbs%isorb_par', subname)
  allocate(orbs%onWhichMPI(Gorbs%norb), stat=i_stat)
  call memocc(i_stat, orbs%onWhichMPI, 'orbs%onWhichMPI', subname)
  iiorb=0
  orbs%isorb_par=0
  do jproc=0,nproc-1
      do iorb=1,orbs%norb_par(jproc,0)
          iiorb=iiorb+1
          orbs%onWhichMPI(iiorb)=jproc
      end do
      if(iproc==jproc) then
          orbs%isorb_par(jproc)=orbs%isorb
      end if
  end do
  call mpiallred(orbs%isorb_par(0), nproc, mpi_sum, mpi_comm_world, ierr)

END SUBROUTINE linear_orbitals_descriptors

