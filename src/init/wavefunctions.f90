!> @file
!!  Routines related to the definition of the wavefunctions
!! @author
!!    Copyright (C) 2010-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Define the descriptors of the orbitals from a given norb
!! It uses the cubic strategy for partitioning the orbitals
!! @param basedist   optional argument indicating the base orbitals distribution to start from
subroutine orbitals_descriptors(iproc,nproc,norb,norbu,norbd,nspin,nspinor,nkpt,kpt,wkpt,&
     orbs,simple,basedist)
  use module_base
  use module_types
  implicit none
  logical, intent(in) :: simple !< simple calculation of the repartition
  integer, intent(in) :: iproc,nproc,norb,norbu,norbd,nkpt,nspin
  integer, intent(in) :: nspinor
  type(orbitals_data), intent(inout) :: orbs
  real(gp), dimension(nkpt), intent(in) :: wkpt
  real(gp), dimension(3,nkpt), intent(in) :: kpt
  integer, dimension(0:nproc-1,nkpt), intent(in), optional :: basedist !> optional argument indicating the base orbitals distribution to start from
  !local variables
  character(len=*), parameter :: subname='orbitals_descriptors'
  integer :: iorb,jproc,norb_tot,ikpt,i_stat,jorb,ierr,i_all,norb_base,iiorb,mpiflag
  logical, dimension(:), allocatable :: GPU_for_orbs
  integer, dimension(:,:), allocatable :: norb_par !(with k-pts)

  !eTS value, updated in evaltocc
  orbs%eTS=0.0_gp

  orbs%norb_par = f_malloc_ptr((/ 0.to.nproc-1 , 0.to.nkpt+ndebug /),id='orbs%norb_par')

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
  call to_zero(nproc*(nkpt+1),orbs%norb_par(0,0))

  !create an array which indicate which processor has a GPU associated 
  !from the viewpoint of the BLAS routines (deprecated, not used anymore)
  if (.not. GPUshare) then
     allocate(GPU_for_orbs(0:nproc-1+ndebug),stat=i_stat)
     call memocc(i_stat,GPU_for_orbs,'GPU_for_orbs',subname)
     
     if (nproc > 1) then
        call MPI_ALLGATHER(GPUconv,1,MPI_LOGICAL,GPU_for_orbs(0),1,MPI_LOGICAL,&
             bigdft_mpi%mpi_comm,ierr)
     else
        GPU_for_orbs(0)=GPUconv
     end if
     
     i_all=-product(shape(GPU_for_orbs))*kind(GPU_for_orbs)
     deallocate(GPU_for_orbs,stat=i_stat)
     call memocc(i_stat,i_all,'GPU_for_orbs',subname)
  end if

  allocate(norb_par(0:nproc-1,orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,norb_par,'norb_par',subname)

  !old system for calculating k-point repartition
!!$  call parallel_repartition_with_kpoints(nproc,orbs%nkpts,norb,orbs%norb_par)
!!$
!!$  !check the distribution
!!$  norb_tot=0
!!$  do jproc=0,iproc-1
!!$     norb_tot=norb_tot+orbs%norb_par(jproc)
!!$  end do
!!$  !reference orbital for process
!!$  orbs%isorb=norb_tot
!!$  do jproc=iproc,nproc-1
!!$     norb_tot=norb_tot+orbs%norb_par(jproc)
!!$  end do
!!$
!!$  if(norb_tot /= norb*orbs%nkpts) then
!!$     write(*,*)'ERROR: partition of orbitals incorrect, report bug.'
!!$     write(*,*)orbs%norb_par(:),norb*orbs%nkpts
!!$     stop
!!$  end if
!!$
!!$  !calculate the k-points related quantities
!!$  allocate(mykpts(orbs%nkpts+ndebug),stat=i_stat)
!!$  call memocc(i_stat,mykpts,'mykpts',subname)
!!$
!!$  call parallel_repartition_per_kpoints(iproc,nproc,orbs%nkpts,norb,orbs%norb_par,&
!!$       orbs%nkptsp,mykpts,norb_par)
!!$  if (orbs%norb_par(iproc) >0) then
!!$     orbs%iskpts=mykpts(1)-1
!!$  else
!!$     orbs%iskpts=0
!!$  end if
!!$  i_all=-product(shape(mykpts))*kind(mykpts)
!!$  deallocate(mykpts,stat=i_stat)
!!$  call memocc(i_stat,i_all,'mykpts',subname)

  !new system for k-point repartition
  norb_base=0
  if (present(basedist)) then
     !the first k-point takes the number of orbitals
     do jproc=0,nproc-1
        norb_base=norb_base+basedist(jproc,1)
     end do
     call components_kpt_distribution(nproc,orbs%nkpts,norb_base,norb,basedist,norb_par)
  else
     call kpts_to_procs_via_obj(nproc,orbs%nkpts,norb,norb_par)
  end if
  !assign the values for norb_par and check the distribution
  norb_tot=0
  do jproc=0,nproc-1
     if (jproc==iproc) orbs%isorb=norb_tot
     do ikpt=1,orbs%nkpts
        orbs%norb_par(jproc,0)=orbs%norb_par(jproc,0)+norb_par(jproc,ikpt)
        orbs%norb_par(jproc,ikpt)=norb_par(jproc,ikpt)
     end do
     norb_tot=norb_tot+orbs%norb_par(jproc,0)
  end do

  if(norb_tot /= norb*orbs%nkpts) then
     write(*,*)'ERROR: partition of orbitals incorrect, report bug.'
     write(*,*)orbs%norb_par(:,0),norb*orbs%nkpts
     stop
  end if

  !allocate(orbs%ikptsp(orbs%nkptsp+ndebug),stat=i_stat)
  !call memocc(i_stat,orbs%ikptsp,'orbs%ikptsp',subname)
  !orbs%ikptsp(1:orbs%nkptsp)=mykpts(1:orbs%nkptsp)

  !this array will be reconstructed in the orbitals_communicators routine
  i_all=-product(shape(norb_par))*kind(norb_par)
  deallocate(norb_par,stat=i_stat)
  call memocc(i_stat,i_all,'norb_par',subname)

  !assign the values of the orbitals data
  orbs%norb=norb
  orbs%norbp=orbs%norb_par(iproc,0)
  orbs%norbu=norbu
  orbs%norbd=norbd


 ! Modify these values
  if (simple) then
     call repartitionOrbitals2(iproc,nproc,orbs%norb,orbs%norb_par,&
          orbs%norbp,orbs%isorb)
  end if

  orbs%iokpt = f_malloc_ptr(orbs%norbp+ndebug,id='orbs%iokpt')

  !assign the k-point to the given orbital, counting one orbital after each other
  jorb=0
  do ikpt=1,orbs%nkpts
     do iorb=1,orbs%norb
        jorb=jorb+1 !this runs over norb*nkpts values
        if (jorb > orbs%isorb .and. jorb <= orbs%isorb+orbs%norbp) then
           orbs%iokpt(jorb-orbs%isorb)=ikpt
        end if
     end do
  end do

  !allocate occupation number and spinsign
  !fill them in normal way
  orbs%occup = f_malloc_ptr(orbs%norb*orbs%nkpts+ndebug,id='orbs%occup')
  call memocc(i_stat,orbs%spinsgn,'orbs%spinsgn',subname)
  orbs%occup(1:orbs%norb*orbs%nkpts)=1.0_gp 
  do ikpt=1,orbs%nkpts
     do iorb=1,orbs%norbu
        orbs%spinsgn(iorb+(ikpt-1)*orbs%norb)=1.0_gp
     end do
     do iorb=1,orbs%norbd
        orbs%spinsgn(iorb+orbs%norbu+(ikpt-1)*orbs%norb)=-1.0_gp
     end do
  end do

  !put a default value for the fermi energy
  orbs%efermi = UNINITIALIZED(orbs%efermi)
  !and also for the gap
  orbs%HLgap = UNINITIALIZED(orbs%HLgap)

  ! allocate inwhichlocreg
  orbs%inwhichlocreg = f_malloc_ptr(orbs%norb*orbs%nkpts,id='orbs%inwhichlocreg')
  ! default for inwhichlocreg (all orbitals are situated in the same locreg)
  orbs%inwhichlocreg = 1

  ! allocate onwhichatom
  orbs%onwhichatom = f_malloc_ptr(orbs%norb*orbs%nkpts,id='orbs%onwhichatom')
  ! default for onwhichatom (all orbitals are situated in the same locreg)
  orbs%onwhichatom = 1

  !initialize the starting point of the potential for each orbital (to be removed?)
  allocate(orbs%ispot(orbs%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%ispot,'orbs%ispot',subname)


  !allocate the array which assign the k-point to processor in transposed version
  orbs%ikptproc = f_malloc_ptr(orbs%nkpts+ndebug,id='orbs%ikptproc')

  ! Define two new arrays:
  ! - orbs%isorb_par is the same as orbs%isorb, but every process also knows
  !   the reference orbital of each other process.
  ! - orbs%onWhichMPI indicates on which MPI process a given orbital
  !   is located.
  orbs%isorb_par = f_malloc_ptr(0.to.nproc-1,id='orbs%isorb_par')
  iiorb=0
  orbs%isorb_par=0
  do jproc=0,nproc-1
      if(iproc==jproc) then
          orbs%isorb_par(jproc)=orbs%isorb
      end if
  end do

  !this mpiflag is added to make memguess working
  call MPI_Initialized(mpiflag,ierr)
  if(nproc >1 .and. mpiflag /= 0) &
       call mpiallred(orbs%isorb_par(0),nproc,mpi_sum,bigdft_mpi%mpi_comm)

END SUBROUTINE orbitals_descriptors


subroutine repartitionOrbitals(iproc,nproc,norb,norb_par,norbp,isorb_par,isorb,onWhichMPI)
  use module_types
  use module_base
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, norb
  integer,dimension(0:nproc-1),intent(out):: norb_par, isorb_par
  integer,dimension(norb),intent(out):: onWhichMPI
  integer,intent(out):: norbp, isorb

  ! Local variables
  integer:: ii, kk, iiorb, iorb, ierr, jproc
  real(8):: tt

  ! Determine norb_par
  norb_par=0
  tt=dble(norb)/dble(nproc)
  ii=floor(tt)
  ! ii is now the number of orbitals that every process has. Distribute the remaining ones.
  norb_par(0:nproc-1)=ii
  kk=norb-nproc*ii
  norb_par(0:kk-1)=ii+1

  ! Determine norbp
  norbp=norb_par(iproc)

  ! Determine isorb
  isorb=0
  do jproc=0,iproc-1
      isorb=isorb+norb_par(jproc)
  end do

  ! Determine onWhichMPI and isorb_par
  iiorb=0
  isorb_par=0
  do jproc=0,nproc-1
      do iorb=1,norb_par(jproc)
          iiorb=iiorb+1
          onWhichMPI(iiorb)=jproc
      end do
      if(iproc==jproc) then
          isorb_par(jproc)=isorb
      end if
  end do
  !call MPI_Initialized(mpiflag,ierr)
  if(nproc >1) &!mpiflag /= 0) 
       call mpiallred(isorb_par(0), nproc, mpi_sum, bigdft_mpi%mpi_comm)

end subroutine repartitionOrbitals


subroutine repartitionOrbitals2(iproc, nproc, norb, norb_par, norbp, isorb)
  use module_base
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, norb
  integer,dimension(0:nproc-1),intent(out):: norb_par
  integer,intent(out):: norbp, isorb

  ! Local variables
  integer:: ii, kk, jproc
  real(8):: tt

  ! Determine norb_par
  norb_par=0
  tt=dble(norb)/dble(nproc)
  ii=floor(tt)
  ! ii is now the number of orbitals that every process has. Distribute the remaining ones.
  norb_par(0:nproc-1)=ii
  kk=norb-nproc*ii
  norb_par(0:kk-1)=ii+1

  ! Determine norbp
  norbp=norb_par(iproc)

  ! Determine isorb
  isorb=0
  do jproc=0,iproc-1
      isorb=isorb+norb_par(jproc)
  end do


end subroutine repartitionOrbitals2

subroutine lzd_set_hgrids(Lzd, hgrids)
  use module_base
  use module_types
  implicit none
  type(local_zone_descriptors), intent(inout) :: Lzd
  real(gp), intent(in) :: hgrids(3)
  !initial values
  Lzd%hgrids = hgrids
END SUBROUTINE lzd_set_hgrids

!> Fill the arrays occup and spinsgn
!! if iunit /=0 this means that the file 'input.occ' does exist and it opens
subroutine occupation_input_variables(verb,iunit,nelec,norb,norbu,norbuempty,norbdempty,nspin,occup,spinsgn)
  use module_base
  use module_input
  use yaml_output
  implicit none
  ! Arguments
  logical, intent(in) :: verb
  integer, intent(in) :: nelec,nspin,norb,norbu,iunit,norbuempty,norbdempty
  real(gp), dimension(norb), intent(out) :: occup,spinsgn
  ! Local variables
  integer :: iorb,nt,ne,it,ierror,iorb1,i
  real(gp) :: rocc
  character(len=20) :: string
  character(len=100) :: line

  do iorb=1,norb
     spinsgn(iorb)=1.0_gp
  end do
  if (nspin/=1) then
     do iorb=1,norbu
        spinsgn(iorb)=1.0_gp
     end do
     do iorb=norbu+1,norb
        spinsgn(iorb)=-1.0_gp
     end do
  end if
  ! write(*,'(1x,a,5i4,30f6.2)')'Spins: ',norb,norbu,norbd,norbup,norbdp,(spinsgn(iorb),iorb=1,norb)

  ! First fill the occupation numbers by default
  nt=0
  if (nspin==1) then
     ne=(nelec+1)/2
     do iorb=1,ne
        it=min(2,nelec-nt)
        occup(iorb)=real(it,gp)
        nt=nt+it
     enddo
     do iorb=ne+1,norb
        occup(iorb)=0._gp
     end do
  else
     if (norbuempty+norbdempty == 0) then
        if (norb > nelec) then
           do iorb=1,min(norbu,norb/2+1)
              it=min(1,nelec-nt)
              occup(iorb)=real(it,gp)
              nt=nt+it
           enddo
           do iorb=min(norbu,norb/2+1)+1,norbu
              occup(iorb)=0.0_gp
           end do
           do iorb=norbu+1,norbu+min(norb-norbu,norb/2+1)
              it=min(1,nelec-nt)
              occup(iorb)=real(it,gp)
              nt=nt+it
           enddo
           do iorb=norbu+min(norb-norbu,norb/2+1)+1,norb
              occup(iorb)=0.0_gp
           end do
        else
           do iorb=1,norb
              occup(iorb)=1.0_gp
           end do
        end if
     else
        do iorb=1,norbu-norbuempty
           occup(iorb)=1.0_gp
        end do
        do iorb=norbu-norbuempty+1,norbu
           occup(iorb)=0.0_gp
        end do
        do iorb=1,norb-norbu-norbdempty
           occup(norbu+iorb)=1.0_gp
        end do
        do iorb=norb-norbu-norbdempty+1,norb-norbu
           occup(norbu+iorb)=0.0_gp
        end do
     end if
  end if
  ! Then read the file "input.occ" if does exist
  if (iunit /= 0) then
     nt=0
     do
        read(unit=iunit,fmt='(a100)',iostat=ierror) line
        if (ierror /= 0) then
           exit
        end if
        !Transform the line in case there are slashes (to ease the parsing)
        do i=1,len(line)
           if (line(i:i) == '/') then
              line(i:i) = ':'
           end if
        end do
        read(line,*,iostat=ierror) iorb,string
        call read_fraction_string(string,rocc,ierror) 
        if (ierror /= 0) then
           exit
        end if

        if (ierror/=0) then
           exit
        else
           nt=nt+1
           if (iorb<0 .or. iorb>norb) then
              !if (iproc==0) then
              write(*,'(1x,a,i0,a)') 'ERROR in line ',nt+1,' of the file "[name].occ"'
              write(*,'(10x,a,i0,a)') 'The orbital index ',iorb,' is incorrect'
              !end if
              stop
           elseif (rocc<0._gp .or. rocc>2._gp) then
              !if (iproc==0) then
              write(*,'(1x,a,i0,a)') 'ERROR in line ',nt+1,' of the file "[name].occ"'
              write(*,'(10x,a,f5.2,a)') 'The occupation number ',rocc,' is not between 0. and 2.'
              !end if
              stop
           else
              occup(iorb)=rocc
           end if
        end if
     end do
     if (verb) then
        call yaml_comment('('//adjustl(trim(yaml_toa(nt)))//'lines read)')
        !write(*,'(1x,a,i0,a)') &
        !     'The occupation numbers are read from the file "[name].occ" (',nt,' lines read)'
     end if
     close(unit=iunit)

     if (nspin/=1) then
!!!        !Check if the polarisation is respected (mpol)
!!!        rup=sum(occup(1:norbu))
!!!        rdown=sum(occup(norbu+1:norb))
!!!        if (abs(rup-rdown-real(norbu-norbd,gp))>1.e-6_gp) then
!!!           if (iproc==0) then
!!!              write(*,'(1x,a,f13.6,a,i0)') 'From the file "input.occ", the polarization ',rup-rdown,&
!!!                             ' is not equal to ',norbu-norbd
!!!           end if
!!!           stop
!!!        end if
        !Fill spinsgn
        do iorb=1,norbu
           spinsgn(iorb)=1.0_gp
        end do
        do iorb=norbu+1,norb
           spinsgn(iorb)=-1.0_gp
        end do
     end if
  end if
  if (verb) then 
     call yaml_sequence(advance='no')
     call yaml_open_map('Occupation Numbers',flow=.true.)
     !write(*,'(1x,a,t28,i8)') 'Total Number of Orbitals',norb
     iorb1=1
     rocc=occup(1)
     do iorb=1,norb
        if (occup(iorb) /= rocc) then
           if (iorb1 == iorb-1) then
              call yaml_map('Orbital No.'//trim(yaml_toa(iorb1)),rocc,fmt='(f6.4)')
              !write(*,'(1x,a,i0,a,f6.4)') 'occup(',iorb1,')= ',rocc
           else
           call yaml_map('Orbitals No.'//trim(yaml_toa(iorb1))//'-'//&
                adjustl(trim(yaml_toa(iorb-1))),rocc,fmt='(f6.4)')
           !write(*,'(1x,a,i0,a,i0,a,f6.4)') 'occup(',iorb1,':',iorb-1,')= ',rocc
           end if
           rocc=occup(iorb)
           iorb1=iorb
        end if
     enddo
     if (iorb1 == norb) then
        call yaml_map('Orbital No.'//trim(yaml_toa(norb)),occup(norb),fmt='(f6.4)')
        !write(*,'(1x,a,i0,a,f6.4)') 'occup(',norb,')= ',occup(norb)
     else
        call yaml_map('Orbitals No.'//trim(yaml_toa(iorb1))//'-'//&
             adjustl(trim(yaml_toa(norb))),occup(norb),fmt='(f6.4)')
        !write(*,'(1x,a,i0,a,i0,a,f6.4)') 'occup(',iorb1,':',norb,')= ',occup(norb)
     end if
     call yaml_close_map()
  endif

  !Check if sum(occup)=nelec
  rocc=sum(occup)
  if (abs(rocc-real(nelec,gp))>1.e-6_gp) then
     call yaml_warning('ERROR in determining the occupation numbers: the total number of electrons ' &
        & // trim(yaml_toa(rocc,fmt='(f13.6)')) // ' is not equal to' // trim(yaml_toa(nelec)))
     !if (iproc==0) then
     !write(*,'(1x,a,f13.6,a,i0)') 'ERROR in determining the occupation numbers: the total number of electrons ',rocc,&
     !     ' is not equal to ',nelec
     !end if
     stop
  end if

END SUBROUTINE occupation_input_variables



