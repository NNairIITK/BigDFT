!>   Generate the input guess via the inguess_generator
subroutine inputguess_gaussian_orbitals_withOnWhichAtom(iproc,nproc,at,rxyz,Glr,nvirt,nspin,&
     orbs,orbse,norbsc_arr,locrad,G,psigau,eks,onWhichAtom)
  use module_base
  use module_types
  use module_interfaces, except_this_one => inputguess_gaussian_orbitals_withOnWhichAtom
  implicit none
  integer, intent(in) :: iproc,nproc,nspin
  integer, intent(inout) :: nvirt
  type(atoms_data), intent(inout) :: at
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: Glr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), intent(out) :: eks
  integer, dimension(at%natsc+1,nspin), intent(out) :: norbsc_arr
  real(gp), dimension(at%nat), intent(out) :: locrad
  type(orbitals_data), intent(inout) :: orbse
  type(gaussian_basis), intent(out) :: G
  real(wp), dimension(:,:,:), pointer :: psigau
  integer,dimension(orbse%norb),intent(out):: onWhichAtom
  !local variables
  character(len=*), parameter :: subname='inputguess_gaussian_orbitals'
  integer, parameter :: ngx=31
  integer :: norbe,norbme,norbyou,i_stat,i_all,norbsc,nvirte,ikpt
  integer :: ispin,jproc,ist,jpst,nspinorfororbse,noncoll
  logical, dimension(:,:,:), allocatable :: scorb
  integer, dimension(:), allocatable :: iorbtolr


  allocate(scorb(4,2,at%natsc+ndebug),stat=i_stat)
  call memocc(i_stat,scorb,'scorb',subname)

  !Generate the input guess via the inguess_generator
  !here we should allocate the gaussian basis descriptors 
  !the prescriptions can be found in the creation of psp basis
  call readAtomicOrbitals_withOnWhichAtom(at,orbse,norbe,norbsc,nspin,orbs%nspinor,&
       scorb,norbsc_arr,locrad,onWhichAtom)

  !in the non-collinear case the number of orbitals double
  if (orbs%nspinor == 4) then
     noncoll=2
  else
     noncoll=1
  end if

  if (iproc ==0) then
     write(*,'(1x,a,i0,a)')'Generating ',nspin*noncoll*norbe,' Atomic Input Orbitals'
     if (norbsc /=0)   write(*,'(1x,a,i0,a)')'  of which ',nspin*noncoll*norbsc,&
          ' are semicore orbitals'
  end if

  if (nvirt /= 0) then
     do ispin=1,nspin
        !Check for max number of virtual orbitals
        !the unoccupied orbitals available as a LCAO
        !this is well defined only for closed-shell systems
        !if (ispin == 1) nvirte=min(noncoll*norbe-orbs%norbu,nvirt)
        !if (ispin == 2) nvirte=min(noncoll*norbe-orbs%norbd,nvirt)
        !alternative test, put always the limit to the number of elements of the input guess
        nvirte=noncoll*norbe
        if(nvirt < nvirte .and. iproc==0) then
           write(*,'(1x,a)')&
                "WARNING: A bigger number of virtual orbitals may be needed for better convergence."
           write(*,'(1x,a,i0)')'         Put nvirt= ',nvirte
        end if
        !if (nvirte < nvirt) then
        !   nvirt=nvirte
        !   if(iproc==0) write(*,'(1x,a,i3)')&
        !        "WARNING: Number of virtual orbitals is too large. New value: ",nvirt
        !end if
        !nvirt=min(nvirt,nvirte)
     end do
  end if

  !allocate communications arrays for virtual orbitals
  !warning: here the aim is just to calculate npsidim, should be fixed
  !call allocate_comms(nproc,orbsv,commsv,subname)
!!$  call orbitals_communicators(iproc,nproc,Glr,orbsv,commsv)  
!!$  call deallocate_comms(commsv,subname)

  !!!orbitals descriptor for inguess orbitals
  nspinorfororbse=orbs%nspinor

  !the number of orbitals to be considered is doubled 
  !in the case of a spin-polarised calculation
  !also for non-collinear case
  !nspin*noncoll is always <= 2
  call orbitals_descriptors(iproc,nproc,nspin*noncoll*norbe,noncoll*norbe,(nspin-1)*norbe, &
       nspin,nspinorfororbse,orbs%nkpts,orbs%kpts,orbs%kwgts,orbse,.false.)
  do ikpt = 1, orbse%nkpts
     ist=1 + (ikpt - 1 ) * nspin*noncoll*norbe
     do ispin=1,nspin
        orbse%spinsgn(ist:ist+norbe-1)=real(1-2*(ispin-1),gp)
        ist=ist+norbe
     end do
  end do

  !this is the distribution procedure for cubic code
  !should be referred to another routine
  if (iproc == 0 .and. nproc > 1) then
     jpst=0
     do jproc=0,nproc-1
        norbme=orbse%norb_par(jproc,0)
        norbyou=orbse%norb_par(min(jproc+1,nproc-1),0)
        if (norbme /= norbyou .or. jproc == nproc-1) then
           !this is a screen output that must be modified
           write(*,'(3(a,i0),a)')&
                ' Processes from ',jpst,' to ',jproc,' treat ',norbme,' inguess orbitals '
           jpst=jproc+1
        end if
     end do
     !write(*,'(3(a,i0),a)')&
     !     ' Processes from ',jpst,' to ',nproc-1,' treat ',norbyou,' inguess orbitals '
  end if

  !write(*,'(a,3i6)') 'iproc, orbse%isorb, orbse%norbp', iproc, orbse%isorb,orbse%norbp
  !write(*,'(a,3i6)') 'norbe, orbse%nspinor, orbse%isorb+orbse%norbp+ndebug', norbe, orbse%nspinor, orbse%isorb+orbse%norbp+ndebug
  !allocate the gaussian coefficients for the number of orbitals which is needed
  allocate(psigau(norbe,orbse%nspinor,orbse%isorb+orbse%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,psigau,'psigau',subname)
  allocate(iorbtolr(orbse%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,iorbtolr,'iorbtolr',subname)

  !fill just the interesting part of the orbital
  call AtomicOrbitals(iproc,at,rxyz,norbe,orbse,norbsc,nspin,eks,scorb,G,&
       psigau(1,1,min(orbse%isorb+1,orbse%norb)),&
       iorbtolr)

  i_all=-product(shape(scorb))*kind(scorb)
  deallocate(scorb,stat=i_stat)
  call memocc(i_stat,i_all,'scorb',subname)

  i_all=-product(shape(iorbtolr))*kind(iorbtolr)
  deallocate(iorbtolr,stat=i_stat)
  call memocc(i_stat,i_all,'iorbtolr',subname)

END SUBROUTINE inputguess_gaussian_orbitals_withOnWhichAtom

