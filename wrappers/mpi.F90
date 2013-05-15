#if defined HAVE_CONFIG_H
#include <config.inc>
#endif

module wrapper_MPI
  ! TO BE REMOVED with f_malloc
  use m_profiling
  ! TO BE REMOVED with f_malloc

  implicit none

  ! MPI handling
#ifdef HAVE_MPI2
  logical, parameter :: have_mpi2 = .true.  !< Flag to use in the code to switch between MPI1 and MPI2
#else
  integer :: MPI_IN_PLACE = 0               !< Fake MPI_IN_PLACE variable to allow compilation in sumrho.
  logical, parameter :: have_mpi2 = .false. !< Flag to use in the code to switch between MPI1 and MPI2
#endif
  include 'mpif.h'      !< MPI definitions and datatypes for density and wavefunctions

  logical :: mpi_thread_funneled_is_supported=.false. !< control the OMP_NESTED based overlap, checked by bigdft_mpi_init below

  !> interface for MPI_ALLREDUCE routine
  interface mpiallred
     module procedure mpiallred_int,mpiallred_real, &
          & mpiallred_double,mpiallred_double_1,mpiallred_double_2,&
          & mpiallred_log
  end interface

  !> global MPI communicator
  type, public :: mpi_environment
     integer :: mpi_comm
     integer :: iproc,nproc
     integer :: igroup,ngroup
  end type mpi_environment

  public :: mpi_environment_null
  public :: mpi_environment_free
  public :: mpi_environment_set
  
contains

  function mpi_environment_null() result(mpi)
    implicit none
    type(mpi_environment) :: mpi
    mpi%mpi_comm=MPI_COMM_WORLD
    mpi%igroup=0
    mpi%ngroup=1
    mpi%iproc=0
    mpi%nproc=1
  end function mpi_environment_null

  subroutine mpi_environment_free(mpi_env)
    implicit none
    type(mpi_environment), intent(inout) :: mpi_env
    !local variables
    integer :: ierr

    if (mpi_env%ngroup > 1) call MPI_COMM_FREE(mpi_env%mpi_comm,ierr)
  end subroutine mpi_environment_free

  !> Set the MPI environment (i.e. taskgroup or MPI communicator)
  !! @param mpi_env   MPI environment (out)
  !! @param iproc     proc id
  !! @param nproc     total number of MPI processes
  !! @param mpi_comm  global MPI_communicator
  !! @param groupsize Number of MPI processes by (task)group
  !!                  if 0 one taskgroup (MPI_COMM_WORLD)
  subroutine mpi_environment_set(mpi_env,iproc,nproc,mpi_comm,groupsize)
    use yaml_output
    implicit none
    integer, intent(in) :: iproc,nproc,mpi_comm,groupsize
    type(mpi_environment), intent(out) :: mpi_env

    mpi_env=mpi_environment_null()

    mpi_env%igroup=0
    mpi_env%ngroup=1
    mpi_env%iproc=iproc
    mpi_env%nproc=nproc
    mpi_env%mpi_comm=mpi_comm

    if (nproc >1 .and. groupsize > 0) then
       if (nproc >1 .and. groupsize < nproc .and. mod(nproc,groupsize)==0) then
          mpi_env%igroup=iproc/groupsize
          mpi_env%ngroup=nproc/groupsize
          mpi_env%iproc=mod(iproc,groupsize)
          mpi_env%nproc=groupsize
          call create_group_comm(mpi_comm,nproc,mpi_env%igroup,mpi_env%nproc,mpi_env%mpi_comm)
          if (iproc == 0) then
             call yaml_map('Total No. of Taskgroups created',nproc/mpi_env%nproc)
          end if
       end if
    end if
  end subroutine mpi_environment_set

  !> create communicators associated to the groups of size group_size
  subroutine create_group_comm(base_comm,nproc_base,group_id,group_size,group_comm)
    use yaml_output
    implicit none
    integer, intent(in) :: base_comm,group_size,nproc_base,group_id
    integer, intent(out) :: group_comm
    !local variables
    character(len=*), parameter :: subname='create_group_comm'
    integer :: grp,ierr,i,j,base_grp,temp_comm,i_stat,i_all
    integer, dimension(:), allocatable :: group_list

    allocate(group_list(group_size+ndebug),stat=i_stat)
    call memocc(i_stat,group_list,'group_list',subname)

    !take the base group
    call MPI_COMM_GROUP(base_comm,base_grp,ierr)
    if (ierr /=0) then
       call yaml_warning('Problem in group creation, ierr:'//yaml_toa(ierr))
       call MPI_ABORT(base_comm,1,ierr)
    end if
    do i=0,nproc_base/group_size-1
       !define the new groups and thread_id
       do j=0,group_size-1
          group_list(j+1)=i*group_size+j
       enddo
       call MPI_GROUP_INCL(base_grp,group_size,group_list,grp,ierr)
       if (ierr /=0) then
          call yaml_warning('Problem in group inclusion, ierr:'//yaml_toa(ierr))
          call MPI_ABORT(base_comm,1,ierr)
       end if
       call MPI_COMM_CREATE(base_comm,grp,temp_comm,ierr)
       if (ierr /=0) then
          call yaml_warning('Problem in communicator creator, ierr:'//yaml_toa(ierr))
          call MPI_ABORT(base_comm,1,ierr)
       end if
       !print *,'i,group_id,temp_comm',i,group_id,temp_comm
       if (i.eq. group_id) group_comm=temp_comm
    enddo

    i_all=-product(shape(group_list ))*kind(group_list )
    deallocate(group_list,stat=i_stat)
    call memocc(i_stat,i_all,'group_list',subname)
  end subroutine create_group_comm

  !> Create a communicator between proc of same rank between the taskgroups.
  subroutine create_rank_comm(group_comm, rank_comm)
    use yaml_output
    implicit none
    integer, intent(in) :: group_comm
    integer, intent(out) :: rank_comm
    !local variables
    character(len=*), parameter :: subname='create_group_master'
    integer :: iproc_group, nproc, nproc_group, ngroups
    integer :: grp, new_grp, ierr, i_stat, i_all, i, j
    integer, dimension(:), allocatable :: lrank, ids

    call mpi_comm_rank(group_comm, iproc_group, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
    call mpi_comm_size(group_comm, nproc_group, ierr)
    ngroups = nproc / nproc_group

    ! Put in lrank the group rank of each process, indexed by global iproc.
    allocate(lrank(nproc+ndebug), stat = i_stat)
    call memocc(i_stat, lrank, 'lrank', subname)
    call mpi_allgather(iproc_group, 1, MPI_INTEGER, lrank, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    ! Put in ids, the global iproc of each process that share the same group iproc.
    allocate(ids(ngroups+ndebug), stat = i_stat)
    call memocc(i_stat, ids, 'ids', subname)
    j = 1
    do i = 1, nproc
       if (lrank(i) == iproc_group) then
          ids(j) = i - 1
          j = j + 1
       end if
    end do
    i_all=-product(shape(lrank ))*kind(lrank )
    deallocate(lrank,stat=i_stat)
    call memocc(i_stat,i_all,'lrank',subname)

!!$    call mpi_comm_rank(MPI_COMM_WORLD, iproc_group, ierr)
!!$    write(*,*) iproc_group, "->", ids
    
    ! Create a new comminucator for the list of ids.
    call mpi_comm_group(MPI_COMM_WORLD, grp, ierr)
    call mpi_group_incl(grp, ngroups, ids, new_grp, ierr)
    call mpi_group_free(grp, ierr)
    i_all=-product(shape(ids ))*kind(ids )
    deallocate(ids,stat=i_stat)
    call memocc(i_stat,i_all,'ids',subname)

    call mpi_comm_create(MPI_COMM_WORLD, new_grp, rank_comm, ierr)
    call mpi_group_free(new_grp, ierr)    
  END SUBROUTINE create_rank_comm

  !interface for MPI_ALLREDUCE operations
  subroutine mpiallred_int(buffer,ntot,mpi_op,mpi_comm,ierr)
    implicit none
    integer, intent(in) :: ntot,mpi_op,mpi_comm
    integer, intent(inout) :: buffer
    integer, intent(out) :: ierr
#ifdef HAVE_MPI2
    !case with MPI_IN_PLACE
    call MPI_ALLREDUCE(MPI_IN_PLACE,buffer,ntot,&
         MPI_INTEGER,mpi_op,mpi_comm,ierr)
#else
    !local variables
    character(len=*), parameter :: subname='mpi_allred'
    integer :: i_all,i_stat
    integer, dimension(:), allocatable :: copybuf

    !case without mpi_in_place
    allocate(copybuf(ntot+ndebug),stat=i_stat)
    call memocc(i_stat,copybuf,'copybuf',subname)

    !not appropriate for integers, to be seen if it works
    call scopy(ntot,buffer,1,copybuf,1) 
    ierr=0 !put just for MPIfake compatibility
    call MPI_ALLREDUCE(copybuf,buffer,ntot,&
         MPI_INTEGER,mpi_op,mpi_comm,ierr)

    i_all=-product(shape(copybuf))*kind(copybuf)
    deallocate(copybuf,stat=i_stat)
    call memocc(i_stat,i_all,'copybuf',subname)
#endif
    if (ierr /=0) stop 'MPIALLRED_INT'

  end subroutine mpiallred_int

  subroutine mpiallred_real(buffer,ntot,mpi_op,mpi_comm,ierr)
    implicit none
    integer, intent(in) :: ntot,mpi_op,mpi_comm
    real(kind=4), intent(inout) :: buffer
    integer, intent(out) :: ierr
#ifdef HAVE_MPI2
    !case with MPI_IN_PLACE
    call MPI_ALLREDUCE(MPI_IN_PLACE,buffer,ntot,&
         MPI_REAL,mpi_op,mpi_comm,ierr)
#else
    !local variables
    character(len=*), parameter :: subname='mpi_allred'
    integer :: i_all,i_stat
    real(kind=4), dimension(:), allocatable :: copybuf

    !case without mpi_in_place
    allocate(copybuf(ntot+ndebug),stat=i_stat)
    call memocc(i_stat,copybuf,'copybuf',subname)

    call scopy(ntot,buffer,1,copybuf,1) 
    ierr=0 !put just for MPIfake compatibility
    call MPI_ALLREDUCE(copybuf,buffer,ntot,&
         MPI_REAL,mpi_op,mpi_comm,ierr)

    i_all=-product(shape(copybuf))*kind(copybuf)
    deallocate(copybuf,stat=i_stat)
    call memocc(i_stat,i_all,'copybuf',subname)
#endif
    if (ierr /=0) stop 'MPIALLRED_REAL'

  end subroutine mpiallred_real

  subroutine mpiallred_double(buffer,ntot,mpi_op,mpi_comm,ierr)
    implicit none
    integer, intent(in) :: ntot,mpi_op,mpi_comm
    real(kind=8), intent(inout) :: buffer
    integer, intent(out) :: ierr
#ifdef HAVE_MPI2
    !case with MPI_IN_PLACE
    call MPI_ALLREDUCE(MPI_IN_PLACE,buffer,ntot,&
         MPI_DOUBLE_PRECISION,mpi_op,mpi_comm,ierr)
#else
    !local variables
    character(len=*), parameter :: subname='mpi_allred'
    integer :: i_all,i_stat
    real(kind=8), dimension(:), allocatable :: copybuf

    !case without mpi_in_place
    allocate(copybuf(ntot+ndebug),stat=i_stat)
    call memocc(i_stat,copybuf,'copybuf',subname)

    call dcopy(ntot,buffer,1,copybuf,1) 
    ierr=0 !put just for MPIfake compatibility
    call MPI_ALLREDUCE(copybuf,buffer,ntot,&
         MPI_DOUBLE_PRECISION,mpi_op,mpi_comm,ierr)

    i_all=-product(shape(copybuf))*kind(copybuf)
    deallocate(copybuf,stat=i_stat)
    call memocc(i_stat,i_all,'copybuf',subname)
#endif
    if (ierr /=0) stop 'MPIALLRED_DBL'
  end subroutine mpiallred_double

  subroutine mpiallred_double_1(buffer,ntot,mpi_op,mpi_comm,ierr)
    implicit none
    integer, intent(in) :: ntot,mpi_op,mpi_comm
    real(kind=8), dimension(:), intent(inout) :: buffer
    integer, intent(out) :: ierr
#ifdef HAVE_MPI2
    !case with MPI_IN_PLACE
    call MPI_ALLREDUCE(MPI_IN_PLACE,buffer,ntot,&
         MPI_DOUBLE_PRECISION,mpi_op,mpi_comm,ierr)
#else
    !local variables
    character(len=*), parameter :: subname='mpi_allred'
    integer :: i_all,i_stat
    real(kind=8), dimension(:), allocatable :: copybuf

    !case without mpi_in_place
    allocate(copybuf(ntot+ndebug),stat=i_stat)
    call memocc(i_stat,copybuf,'copybuf',subname)

    call dcopy(ntot,buffer,1,copybuf,1) 
    ierr=0 !put just for MPIfake compatibility
    call MPI_ALLREDUCE(copybuf,buffer,ntot,&
         MPI_DOUBLE_PRECISION,mpi_op,mpi_comm,ierr)

    i_all=-product(shape(copybuf))*kind(copybuf)
    deallocate(copybuf,stat=i_stat)
    call memocc(i_stat,i_all,'copybuf',subname)
#endif
    if (ierr /=0) stop 'MPIALLRED_DBL'
  end subroutine mpiallred_double_1

  subroutine mpiallred_double_2(buffer,ntot,mpi_op,mpi_comm,ierr)
    implicit none
    integer, intent(in) :: ntot,mpi_op,mpi_comm
    real(kind=8), dimension(:,:), intent(inout) :: buffer
    integer, intent(out) :: ierr
#ifdef HAVE_MPI2
    !case with MPI_IN_PLACE
    call MPI_ALLREDUCE(MPI_IN_PLACE,buffer,ntot,&
         MPI_DOUBLE_PRECISION,mpi_op,mpi_comm,ierr)
#else
    !local variables
    character(len=*), parameter :: subname='mpi_allred'
    integer :: i_all,i_stat
    real(kind=8), dimension(:), allocatable :: copybuf

    !case without mpi_in_place
    allocate(copybuf(ntot+ndebug),stat=i_stat)
    call memocc(i_stat,copybuf,'copybuf',subname)

    call dcopy(ntot,buffer,1,copybuf,1) 
    ierr=0 !put just for MPIfake compatibility
    call MPI_ALLREDUCE(copybuf,buffer,ntot,&
         MPI_DOUBLE_PRECISION,mpi_op,mpi_comm,ierr)

    i_all=-product(shape(copybuf))*kind(copybuf)
    deallocate(copybuf,stat=i_stat)
    call memocc(i_stat,i_all,'copybuf',subname)
#endif
    if (ierr /=0) stop 'MPIALLRED_DBL'
  end subroutine mpiallred_double_2

  subroutine mpiallred_log(buffer,ntot,mpi_op,mpi_comm,ierr)
    implicit none
    integer, intent(in) :: ntot,mpi_op,mpi_comm
    logical, intent(inout) :: buffer
    integer, intent(out) :: ierr
#ifdef HAVE_MPI2
    !case with MPI_IN_PLACE
    call MPI_ALLREDUCE(MPI_IN_PLACE,buffer,ntot,&
         MPI_LOGICAL,mpi_op,mpi_comm,ierr)
#else
    !local variables
    character(len=*), parameter :: subname='mpi_allred'
    integer :: i_all,i_stat
    logical, dimension(:), allocatable :: copybuf

    !case without mpi_in_place
    allocate(copybuf(ntot+ndebug),stat=i_stat)
    call memocc(i_stat,copybuf,'copybuf',subname)

    !not appropriate for logical, to be seen if it works
    call scopy(ntot,buffer,1,copybuf,1) 
    ierr=0 !put just for MPIfake compatibility
    call MPI_ALLREDUCE(copybuf,buffer,ntot,&
         MPI_LOGICAL,mpi_op,mpi_comm,ierr)

    i_all=-product(shape(copybuf))*kind(copybuf)
    deallocate(copybuf,stat=i_stat)
    call memocc(i_stat,i_all,'copybuf',subname)
#endif

    !inform and stop if an error occurs
    if (ierr /=0) stop 'MPIALLRED_LOG'

  end subroutine mpiallred_log

end module wrapper_MPI

subroutine bigdft_mpi_init(ierr)
  use wrapper_mpi
  implicit none
  integer, intent(out) :: ierr
#ifdef HAVE_MPI_INIT_THREAD
  integer :: provided
  call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
  if (ierr /= MPI_SUCCESS) then
     write(*,*)'BigDFT_mpi_INIT: Error in MPI_INIT_THREAD',ierr
  else if (provided < MPI_THREAD_FUNNELED) then
     !write(*,*)'WARNING: MPI_THREAD_FUNNELED not supported!',provided,ierr
     !call MPI_INIT(ierr)
  else
     mpi_thread_funneled_is_supported=.true.
  endif
#else
  call MPI_INIT(ierr)      
  if (ierr /= MPI_SUCCESS) then
     write(*,*)'BigDFT_mpi_INIT: Error in MPI_INIT_THREAD',ierr
  end if
#endif
end subroutine bigdft_mpi_init

!> Activates the nesting for UNBLOCK_COMMS performance case
subroutine bigdft_open_nesting(num_threads)
  use wrapper_mpi
  implicit none
  integer, intent(in) :: num_threads
#ifdef HAVE_MPI_INIT_THREAD
  !$ call OMP_SET_NESTED(.true.) 
  !$ call OMP_SET_MAX_ACTIVE_LEVELS(2)
  !$ call OMP_SET_NUM_THREADS(num_threads)
#else
  integer :: ierr,idummy
  write(*,*)'BigDFT_open_nesting is not active!'
  call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)
  idummy=num_threads
#endif
end subroutine bigdft_open_nesting

!> Activates the nesting for UNBLOCK_COMMS performance case
subroutine bigdft_close_nesting(num_threads)
  use wrapper_mpi
  implicit none
  integer, intent(in) :: num_threads
#ifdef HAVE_MPI_INIT_THREAD
  !$ call OMP_SET_MAX_ACTIVE_LEVELS(1) !redundant
  !$ call OMP_SET_NESTED(.false.) 
  !$ call OMP_SET_NUM_THREADS(num_threads)
#else 
  integer :: ierr,idummy
  write(*,*)'BigDFT_close_nesting is not active!'
  call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)
  idummy=num_threads
#endif
end subroutine bigdft_close_nesting
