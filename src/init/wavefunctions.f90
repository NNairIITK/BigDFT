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

  allocate(orbs%norb_par(0:nproc-1,0:nkpt+ndebug),stat=i_stat)
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

  allocate(orbs%iokpt(orbs%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%iokpt,'orbs%iokpt',subname)

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
  allocate(orbs%occup(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%occup,'orbs%occup',subname)
  allocate(orbs%spinsgn(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
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
  allocate(orbs%inwhichlocreg(orbs%norb*orbs%nkpts),stat=i_stat)
  call memocc(i_stat,orbs%inwhichlocreg,'orbs%inwhichlocreg',subname)
  ! default for inwhichlocreg (all orbitals are situated in the same locreg)
  orbs%inwhichlocreg = 1

  ! allocate onwhichatom
  allocate(orbs%onwhichatom(orbs%norb*orbs%nkpts),stat=i_stat)
  call memocc(i_stat,orbs%onwhichatom,'orbs%onwhichatom',subname)
  ! default for onwhichatom (all orbitals are situated in the same locreg)
  orbs%onwhichatom = 1

  !initialize the starting point of the potential for each orbital (to be removed?)
  allocate(orbs%ispot(orbs%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%ispot,'orbs%ispot',subname)


  !allocate the array which assign the k-point to processor in transposed version
  allocate(orbs%ikptproc(orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%ikptproc,'orbs%ikptproc',subname)

  ! Define two new arrays:
  ! - orbs%isorb_par is the same as orbs%isorb, but every process also knows
  !   the reference orbital of each other process.
  ! - orbs%onWhichMPI indicates on which MPI process a given orbital
  !   is located.
  allocate(orbs%isorb_par(0:nproc-1), stat=i_stat)
  call memocc(i_stat, orbs%isorb_par, 'orbs%isorb_par', subname)
  allocate(orbs%onWhichMPI(sum(orbs%norb_par(:,0))), stat=i_stat)
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
  !this mpiflag is added to make memguess working
  call MPI_Initialized(mpiflag,ierr)
  if(nproc >1 .and. mpiflag /= 0) &
       call mpiallred(orbs%isorb_par(0),nproc,mpi_sum,bigdft_mpi%mpi_comm,ierr)

END SUBROUTINE orbitals_descriptors


!> Partition the orbitals between processors to ensure load balancing
!! the criterion will depend on GPU computation
!! and/or on the sizes of the different localisation region.
!!
!! Calculate the number of elements to be sent to each process
!! and the array of displacements.
!! Cubic strategy: 
!!    - the components are equally distributed among the wavefunctions
!!    - each processor has all the orbitals in transposed form
!!    - each wavefunction is equally distributed in its transposed form
!!    - this holds for each k-point, which regroups different processors
subroutine orbitals_communicators(iproc,nproc,lr,orbs,comms,basedist)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(inout) :: orbs
  type(communications_arrays), intent(out) :: comms
  integer, dimension(0:nproc-1,orbs%nkpts), intent(in), optional :: basedist
  !local variables
  character(len=*), parameter :: subname='orbitals_communicators'
  logical :: yesorb,yescomp
  integer :: jproc,nvctr_tot,ikpts,iorbp,jorb,norb_tot,ikpt,i_stat,i_all
  integer :: nkptsp,ierr,kproc,jkpts,jkpte,jsorb,lubo,lubc,info,jkpt
  integer, dimension(:), allocatable :: mykpts
  logical, dimension(:), allocatable :: GPU_for_comp
  integer, dimension(:,:), allocatable :: nvctr_par,norb_par !<for all the components and orbitals (with k-pts)
  
  !check of allocation of important arrays
  if (.not. associated(orbs%norb_par)) then
     write(*,*)'ERROR: norb_par array not allocated'
     stop
  end if

  allocate(nvctr_par(0:nproc-1,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,nvctr_par,'nvctr_par',subname)
  allocate(norb_par(0:nproc-1,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,norb_par,'norb_par',subname)
  allocate(mykpts(orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,mykpts,'mykpts',subname)

  !initialise the arrays
  do ikpts=0,orbs%nkpts
     do jproc=0,nproc-1
        nvctr_par(jproc,ikpts)=0 
        norb_par(jproc,ikpts)=0 
     end do
  end do

  !calculate the same k-point distribution for the orbitals
  !assign the k-point to the given orbital, counting one orbital after each other
  jorb=1
  ikpts=1
  do jproc=0,nproc-1
     do iorbp=1,orbs%norb_par(jproc,0)
        norb_par(jproc,ikpts)=norb_par(jproc,ikpts)+1
        if (mod(jorb,orbs%norb)==0) then
           ikpts=ikpts+1
        end if
        jorb=jorb+1
     end do
  end do
  !some checks
  if (orbs%norb /= 0) then
     !check the distribution
     do ikpts=1,orbs%nkpts
        !print *,'partition',ikpts,orbs%nkpts,'ikpts',norb_par(:,ikpts)
        norb_tot=0
        do jproc=0,nproc-1
           norb_tot=norb_tot+norb_par(jproc,ikpts)
        end do
        if(norb_tot /= orbs%norb) then
           write(*,*)'ERROR: partition of orbitals incorrect, kpoint:',ikpts
           stop
        end if
     end do
  end if


  !balance the components between processors
  !in the most symmetric way
  !here the components are taken into account for all the k-points

  !create an array which indicate which processor has a GPU associated 
  !from the viewpoint of the BLAS routines (deprecated, not used anymore)
  allocate(GPU_for_comp(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,GPU_for_comp,'GPU_for_comp',subname)

  if (nproc > 1 .and. .not. GPUshare) then
     call MPI_ALLGATHER(GPUblas,1,MPI_LOGICAL,GPU_for_comp(0),1,MPI_LOGICAL,&
          bigdft_mpi%mpi_comm,ierr)
  else
     GPU_for_comp(0)=GPUblas
  end if

  i_all=-product(shape(GPU_for_comp))*kind(GPU_for_comp)
  deallocate(GPU_for_comp,stat=i_stat)
  call memocc(i_stat,i_all,'GPU_for_comp',subname)

  !old k-point repartition
!!$  !decide the repartition for the components in the same way as the orbitals
!!$  call parallel_repartition_with_kpoints(nproc,orbs%nkpts,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),nvctr_par)

!!$  ikpts=1
!!$  ncomp_res=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
!!$  do jproc=0,nproc-1
!!$     loop_comps: do
!!$        if (nvctr_par(jproc,0) >= ncomp_res) then
!!$           nvctr_par(jproc,ikpts)= ncomp_res
!!$           ikpts=ikpts+1
!!$           nvctr_par(jproc,0)=nvctr_par(jproc,0)-ncomp_res
!!$           ncomp_res=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
!!$        else
!!$           nvctr_par(jproc,ikpts)= nvctr_par(jproc,0)
!!$           ncomp_res=ncomp_res-nvctr_par(jproc,0)
!!$           nvctr_par(jproc,0)=0
!!$           exit loop_comps
!!$        end if
!!$        if (nvctr_par(jproc,0) == 0 ) then
!!$           ncomp_res=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
!!$           exit loop_comps
!!$        end if
!!$
!!$     end do loop_comps
!!$  end do

  !new k-point repartition
  if (present(basedist)) then
     do jkpt=1,orbs%nkpts
        do jproc=0,nproc-1
           nvctr_par(jproc,jkpt)=basedist(jproc,jkpt)
        end do
     end do
  else
     !first try the naive repartition
     call kpts_to_procs_via_obj(nproc,orbs%nkpts,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),nvctr_par(0,1))
  end if
  !then silently check whether the distribution agree
  info=-1
  call check_kpt_distributions(nproc,orbs%nkpts,orbs%norb,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),&
       norb_par(0,1),nvctr_par(0,1),info,lubo,lubc)
  if (info/=0 .and. .not. present(basedist)) then !redo the distribution based on the orbitals scheme
     info=-1
     call components_kpt_distribution(nproc,orbs%nkpts,orbs%norb,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),norb_par(0,1),nvctr_par(0,1))
     call check_kpt_distributions(nproc,orbs%nkpts,orbs%norb,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),&
          norb_par(0,1),nvctr_par(0,1),info,lubo,lubc)
  end if
  if (info /=0) then
     if (iproc==0) then
        write(*,*)'ERROR for nproc,nkpts,norb,nvctr',nproc,orbs%nkpts,orbs%norb,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
        call print_distribution_schemes(6,nproc,orbs%nkpts,norb_par(0,1),nvctr_par(0,1))
     end if
     call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
     stop
  end if

!write(*,'(a,i2,3x,8i7,i10)') 'iproc, nvctr_par(jproc), sum', iproc, (nvctr_par(jproc,1), jproc=0,nproc-1), sum(nvctr_par(:,1))
!write(*,*) 'iproc, (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%norbp', iproc, (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%norbp
  !some checks
  !check the distribution
  do ikpts=1,orbs%nkpts
     !print *,'iproc,cpts:',lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,nvctr_par(:,ikpts)
     nvctr_tot=0
     do jproc=0,nproc-1
        nvctr_tot=nvctr_tot+nvctr_par(jproc,ikpts)
     end do
     if(nvctr_tot /= lr%wfd%nvctr_c+7*lr%wfd%nvctr_f) then
        write(*,*)'ERROR: partition of components incorrect, kpoint:',ikpts
        stop
     end if
  end do

  !this function which associates a given k-point to a processor in the component distribution
  !the association is chosen such that each k-point is associated to only
  !one processor
  !if two processors treat the same k-point the processor which highest rank is chosen
  do ikpts=1,orbs%nkpts
     loop_jproc: do jproc=nproc-1,0,-1
        if (nvctr_par(jproc,ikpts) /= 0) then
           orbs%ikptproc(ikpts)=jproc
           exit loop_jproc
        end if
     end do loop_jproc
  end do
  
  !print*,'check',orbs%ikptproc(:)

!write(*,*) 'orbs%norb_par',orbs%norb_par

  !calculate the number of k-points treated by each processor in both
  ! the component distribution and the orbital distribution.
  !to have a correct distribution, a k-point should be divided between the same processors
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

!!$  allocate(orbs%ikptsp(orbs%nkptsp+ndebug),stat=i_stat)
!!$  call memocc(i_stat,orbs%ikptsp,'orbs%ikptsp',subname)
!!$  orbs%ikptsp(1:orbs%nkptsp)=mykpts(1:orbs%nkptsp)

  !print the distribution scheme used for this set of orbital
  !in the case of multiple k-points
  if (iproc == 0 .and. verbose > 1 .and. orbs%nkpts > 1) then
     call print_distribution_schemes(6,nproc,orbs%nkpts,norb_par(0,1),nvctr_par(0,1))
  end if

  !print *,iproc,orbs%nkptsp,orbs%norbp,orbs%norb,orbs%nkpts
  !call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
  !call MPI_FINALIZE(ierr)
  !stop
  !check that for any processor the orbital k-point repartition is contained into the components
  do jproc=0,nproc-1
     jsorb=0
     do kproc=0,jproc-1
        jsorb=jsorb+orbs%norb_par(kproc,0)
     end do
     jkpts=min(jsorb/orbs%norb+1,orbs%nkpts)
     if (nvctr_par(jproc,jkpts) == 0 .and. orbs%norb_par(jproc,0) /=0 ) then
        if (iproc ==0) write(*,*)'ERROR, jproc: ',jproc,' the orbital k-points distribution starts before the components one'
        !print *,jsorb,jkpts,jproc,orbs%iskpts,nvctr_par(jproc,jkpts)
        stop
     end if
     jkpte=min((jsorb+orbs%norb_par(jproc,0)-1)/orbs%norb+1,orbs%nkpts)
     if (nvctr_par(jproc,jkpte) == 0 .and. orbs%norb_par(jproc,0) /=0) then
        if (iproc ==0) write(*,*)'ERROR, jproc: ',jproc,&
             ' the orbital k-points distribution ends after the components one'
        print *,jsorb,jkpte,jproc,orbs%iskpts,orbs%nkptsp,nvctr_par(jproc,jkpte)
        stop
     end if
  end do

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
        call MPI_ABORT(bigdft_mpi%mpi_comm, ierr)
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
        call MPI_ABORT(bigdft_mpi%mpi_comm, ierr)
     end if
  end do kpt_orbitals

  !print *,'AAAAiproc',iproc,orbs%iskpts,orbs%iskpts+orbs%nkptsp

  !allocate communication arrays
  allocate(comms%nvctr_par(0:nproc-1,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,comms%nvctr_par,'nvctr_par',subname)

  allocate(comms%ncntd(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,comms%ncntd,'ncntd',subname)

  allocate(comms%ncntt(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,comms%ncntt,'ncntt',subname)
  allocate(comms%ndspld(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,comms%ndspld,'ndspld',subname)
  allocate(comms%ndsplt(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,comms%ndsplt,'ndsplt',subname)

  !assign the partition of the k-points to the communication array
  !calculate the number of componenets associated to the k-point
  do jproc=0,nproc-1
     comms%nvctr_par(jproc,0)=0
     do ikpt=1,orbs%nkpts
        comms%nvctr_par(jproc,0)=comms%nvctr_par(jproc,0)+&
             nvctr_par(jproc,ikpt) 
        comms%nvctr_par(jproc,ikpt)=nvctr_par(jproc,ikpt)
     end do
  end do
!!$  do ikpts=1,orbs%nkptsp
!!$     ikpt=orbs%iskpts+ikpts!orbs%ikptsp(ikpts)
!!$     do jproc=0,nproc-1
!!$        comms%nvctr_par(jproc,ikpts)=nvctr_par(jproc,ikpt) 
!!$     end do
!!$  end do

  !with this distribution the orbitals and the components are ordered following k-points
  !there must be no overlap for the components
  !here we will print out the k-points components distribution, in the transposed and in the direct way

  do jproc=0,nproc-1
     comms%ncntd(jproc)=0
     do ikpts=1,orbs%nkpts
        comms%ncntd(jproc)=comms%ncntd(jproc)+&
             nvctr_par(jproc,ikpts)*norb_par(iproc,ikpts)*orbs%nspinor
     end do
  end do
  comms%ndspld(0)=0
  do jproc=1,nproc-1
     comms%ndspld(jproc)=comms%ndspld(jproc-1)+comms%ncntd(jproc-1)
  end do
  !receive buffer
  do jproc=0,nproc-1
     comms%ncntt(jproc)=0
     do ikpts=1,orbs%nkpts
        comms%ncntt(jproc)=comms%ncntt(jproc)+&
             nvctr_par(iproc,ikpts)*norb_par(jproc,ikpts)*orbs%nspinor
     end do
  end do
  comms%ndsplt(0)=0
  do jproc=1,nproc-1
     comms%ndsplt(jproc)=comms%ndsplt(jproc-1)+comms%ncntt(jproc-1)
  end do

  !print *,'iproc,comms',iproc,comms%ncntd,comms%ndspld,comms%ncntt,comms%ndsplt

  i_all=-product(shape(nvctr_par))*kind(nvctr_par)
  deallocate(nvctr_par,stat=i_stat)
  call memocc(i_stat,i_all,'nvctr_par',subname)
  i_all=-product(shape(norb_par))*kind(norb_par)
  deallocate(norb_par,stat=i_stat)
  call memocc(i_stat,i_all,'norb_par',subname)
  i_all=-product(shape(mykpts))*kind(mykpts)
  deallocate(mykpts,stat=i_stat)
  call memocc(i_stat,i_all,'mykpts',subname)

  !calculate the dimension of the wavefunction
  !for the given processor (this is only the cubic strategy)
  orbs%npsidim_orbs=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%norb_par(iproc,0)*orbs%nspinor
  orbs%npsidim_comp=sum(comms%ncntt(0:nproc-1))
    
!!$  orbs%npsidim=max((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%norb_par(iproc,0)*orbs%nspinor,&
!!$       sum(comms%ncntt(0:nproc-1)))

END SUBROUTINE orbitals_communicators


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
       call mpiallred(isorb_par(0), nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)

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


subroutine inputs_parse_params(in, iproc, dump)
  use module_types
  use module_xc
  implicit none
  type(input_variables), intent(inout) :: in
  integer, intent(in) :: iproc
  logical, intent(in) :: dump

  ! Parse all values independant from atoms.
  call perf_input_variables(iproc,dump,trim(in%file_perf),in)
  call dft_input_variables_new(iproc,dump,trim(in%file_dft),in)
  call mix_input_variables_new(iproc,dump,trim(in%file_mix),in)
  call geopt_input_variables_new(iproc,dump,trim(in%file_geopt),in)
  call tddft_input_variables_new(iproc,dump,trim(in%file_tddft),in)
  call sic_input_variables_new(iproc,dump,trim(in%file_sic),in)

  ! Initialise XC calculation
  if (in%ixc < 0) then
     call xc_init(in%ixc, XC_MIXED, in%nspin)
  else
     call xc_init(in%ixc, XC_ABINIT, in%nspin)
  end if
end subroutine inputs_parse_params


subroutine inputs_parse_add(in, atoms, iproc, dump)
  use module_types
  implicit none
  type(input_variables), intent(inout) :: in
  type(atoms_data), intent(inout) :: atoms
  integer, intent(in) :: iproc
  logical, intent(in) :: dump

  ! Read k-points input variables (if given)
  call kpt_input_variables_new(iproc,dump,trim(in%file_kpt),in,atoms%sym,atoms%geocode, &
       & (/ atoms%alat1, atoms%alat2, atoms%alat3 /))

  ! Linear scaling (if given)
  in%lin%fragment_calculation=.false. ! to make sure that if we're not doing a linear calculation we don't read fragment information
  call lin_input_variables_new(iproc,dump .and. (in%inputPsiId == INPUT_PSI_LINEAR_AO .or. &
       & in%inputPsiId == INPUT_PSI_DISK_LINEAR), trim(in%file_lin),in,atoms)

  ! Fragment information (if given)
  call fragment_input_variables(iproc,dump .and. (in%inputPsiId == INPUT_PSI_LINEAR_AO .or. &
       & in%inputPsiId == INPUT_PSI_DISK_LINEAR), trim(in%file_frag),in,atoms)

end subroutine inputs_parse_add


