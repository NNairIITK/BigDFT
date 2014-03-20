!> @file
!! Wrapper for the MPI call (this file is preprocessed.)
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


#if defined HAVE_CONFIG_H
#include <config.inc>
#endif

!> Module defining the routines which wrap the MPI calls
module wrapper_MPI
  ! TO BE REMOVED with f_malloc
  !use memory_profiling!, only: ndebug
  use dynamic_memory
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

  logical :: mpi_thread_funneled_is_supported=.false. !< Control the OMP_NESTED based overlap, checked by bigdft_mpi_init below

  !> Interface for MPI_ALLREDUCE routine
  interface mpiallred
     module procedure mpiallred_int,mpiallred_real, &
          & mpiallred_double,mpiallred_double_1,mpiallred_double_2,&
          & mpiallred_log
  end interface mpiallred

  !> Interface for MPI_ALLGATHERV routine
  interface mpiallgatherv
     module procedure mpiallgatherv_double
  end interface mpiallgatherv

  !> Global MPI communicator which contains all information related to the MPI process
  type, public :: mpi_environment
     integer :: mpi_comm !< MPI communicator
     integer :: iproc    !< Process Id
                         !! @ingroup RESERVED
     integer :: nproc    !< Number of MPI processes (in the given communicator)
                         !! @ingroup RESERVED
     integer :: igroup   !< MPI Group Id
     integer :: ngroup   !< Number of MPI groups
  end type mpi_environment

  public :: mpi_environment_null
  public :: mpi_environment_free
  public :: mpi_environment_set
  public :: mpi_environment_set1 !to be removed
  
contains

  pure function mpi_environment_null() result(mpi)
    implicit none
    type(mpi_environment) :: mpi
    mpi%mpi_comm=MPI_COMM_NULL !better to put an invalid comm?
    mpi%igroup=-1
    mpi%ngroup=-1
    mpi%iproc=-1
    mpi%nproc=-1
  end function mpi_environment_null

  subroutine mpi_environment_free(mpi_env)
    implicit none
    type(mpi_environment), intent(inout) :: mpi_env
    !local variables
    integer :: ierr

    if (mpi_env%mpi_comm /= MPI_COMM_WORLD .and. &
         mpi_env%mpi_comm /= MPI_COMM_NULL) &
         call MPI_COMM_FREE(mpi_env%mpi_comm,ierr)
    mpi_env=mpi_environment_null()
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
    !local variables
    integer :: j
    integer, dimension(:), allocatable :: group_list

    call f_routine(id='mpi_environment_set')
    mpi_env=mpi_environment_null()

!!$    mpi_env%igroup=0
!!$    mpi_env%ngroup=1
!!$    mpi_env%iproc=iproc
!!$    mpi_env%nproc=nproc
    mpi_env%mpi_comm=mpi_comm
    mpi_env%igroup=iproc/groupsize
    mpi_env%ngroup=nproc/groupsize
    mpi_env%iproc=mod(iproc,groupsize)
    mpi_env%nproc=groupsize
    if (groupsize /= nproc) then
       !define the strategy for the taskgroups
       group_list=f_malloc(groupsize,id='group_list')
       !iproc in the same group are close to each other
       do j=0,groupsize-1
          group_list(j+1)=mpi_env%igroup*groupsize+j
       enddo
       call create_group_comm(mpi_comm,groupsize,group_list,mpi_env%mpi_comm)
       if (iproc == 0) then
          call yaml_map('Total No. of Taskgroups created',nproc/mpi_env%nproc)
       end if
       call f_free(group_list)
    end if

    call f_release_routine()
  end subroutine mpi_environment_set

!!! PSolver n1-n2 plane mpi partitioning !!! 
  !> This is exactly like mpi_environment_set but it always creates groups
  !! the routine above should be modified accordingly
!!$  subroutine mpi_environment_set2(mpi_env,iproc,nproc,mpi_comm,groupsize)
!!$    use yaml_output
!!$    implicit none
!!$    integer, intent(in) :: iproc,nproc,mpi_comm,groupsize
!!$    type(mpi_environment), intent(out) :: mpi_env
!!$    !local variables
!!$    integer :: j
!!$    integer, dimension(:), allocatable :: group_list
!!$
!!$    call f_routine(id='mpi_environment_set2')
!!$    mpi_env=mpi_environment_null()
!!$
!!$    mpi_env%mpi_comm=mpi_comm
!!$
!!$    mpi_env%igroup=iproc/groupsize
!!$    mpi_env%ngroup=nproc/groupsize
!!$    mpi_env%iproc=mod(iproc,groupsize)
!!$    mpi_env%nproc=groupsize
!!$
!!$    !define the strategy for the taskgroups
!!$    group_list=f_malloc(groupsize,id='group_list')
!!$    !iproc in the same group are close to each other
!!$    do j=0,groupsize-1
!!$       group_list(j+1)=mpi_env%igroup*groupsize+j
!!$    enddo
!!$
!!$    call create_group_comm(mpi_comm,nproc,mpi_env%igroup,mpi_env%nproc,group_list,mpi_env%mpi_comm)
!!$!    if (iproc == 0) then
!!$!       call yaml_map('Total No. of Taskgroups created',nproc/mpi_env%nproc)
!!$!    end if
!!$    call f_free(group_list)
!!$    call f_release_routine()
!!$  end subroutine mpi_environment_set2

  !this is a different procedure to assign the iproc according to the groups.
  subroutine mpi_environment_set1(mpi_env,iproc,nproc,mpi_comm,groupsize,ngroup)
    use yaml_output
    implicit none
    integer, intent(in) :: iproc,nproc,mpi_comm,groupsize,ngroup
    type(mpi_environment), intent(out) :: mpi_env
    !local variables
    integer :: j
    integer, dimension(:), allocatable :: group_list

    call f_routine(id='mpi_environment_set1')

    mpi_env=mpi_environment_null()

    mpi_env%igroup=-1

    mpi_env%ngroup=ngroup
    if (iproc < groupsize*ngroup) mpi_env%igroup=mod(iproc,ngroup)
    mpi_env%iproc=iproc/ngroup
    mpi_env%nproc=groupsize
    mpi_env%mpi_comm=mpi_comm

    !define the strategy for the taskgroups
    group_list=f_malloc(groupsize,id='group_list')
    !round-robin strategy
    if (mpi_env%igroup >0) then
       do j=0,groupsize-1
          group_list(j+1)=mpi_env%igroup+j*mpi_env%ngroup
       enddo
    else
       !these processes have MPI_COMM_NULL
       group_list=-1
       mpi_env%mpi_comm=MPI_COMM_NULL
    end if

    !call create_group_comm1(mpi_comm,nproc,mpi_env%igroup,ngroup,mpi_env%nproc,mpi_env%mpi_comm)
    call create_group_comm(mpi_comm,mpi_env%nproc,group_list,mpi_env%mpi_comm)
!    if (iproc == 0) then
!       call yaml_map('Total No. of Taskgroups created',ngroup)
!    end if
    call f_free(group_list)
    call f_release_routine()
  end subroutine mpi_environment_set1

  !> create communicators associated to the groups of size group_size
  subroutine create_group_comm(base_comm,group_size,group_list,group_comm)
    use yaml_output
    use dictionaries

    implicit none
    integer, intent(in) :: base_comm,group_size
    integer, dimension(group_size), intent(in) :: group_list !< list of id of the group identified by group_id in units of base_comm
    integer, intent(out) :: group_comm
    !local variables
    integer :: grp,ierr,base_grp

    !take the base group
    call MPI_COMM_GROUP(base_comm,base_grp,ierr)
    if (f_err_raise(ierr/=0,'Problem in group creation, ierr:'//yaml_toa(ierr),&
         err_name='BIGDFT_MPI_ERROR')) return
    !create the groups with the list
    call MPI_GROUP_INCL(base_grp,group_size,group_list,grp,ierr)
    if (f_err_raise(ierr/=0,'Problem in group inclusion, ierr:'//yaml_toa(ierr),&
         err_name='BIGDFT_MPI_ERROR')) return
    !free base group
    call MPI_GROUP_FREE(base_grp,ierr)
    if (f_err_raise(ierr/=0,'Problem in base_group free, ierr:'//yaml_toa(ierr),&
         err_name='BIGDFT_MPI_ERROR')) return
    !create the communicator (the communicator can be also null)
    call MPI_COMM_CREATE(base_comm,grp,group_comm,ierr)
    if (f_err_raise(ierr/=0,'Problem in communicator creator, ierr:'//yaml_toa(ierr),&
         err_name='BIGDFT_MPI_ERROR')) return
    !free temporary group
    call MPI_GROUP_FREE(grp,ierr)
    if (f_err_raise(ierr/=0,'Problem in new_group free, ierr:'//yaml_toa(ierr),&
              err_name='BIGDFT_MPI_ERROR')) return

  end subroutine create_group_comm

!!! PSolver n1-n2 plane mpi partitioning !!! 
!this routine is like create_group_comm with a different group_list
subroutine create_group_comm1(base_comm,nproc_base,group_id,ngroup,group_size,group_comm)
  use yaml_output
  implicit none
  integer, intent(in) :: base_comm,group_size,nproc_base,group_id,ngroup
  integer, intent(out) :: group_comm
  !local variables
  character(len=*), parameter :: subname='create_group_comm'
  integer :: grp,ierr,i,j,base_grp,temp_comm!,i_stat,i_all
  integer, dimension(:), allocatable :: group_list

! allocate(group_list(group_size+ndebug),stat=i_stat)
  group_list = f_malloc(group_size,id='group_list')

  !take the base group
  call MPI_COMM_GROUP(base_comm,base_grp,ierr)
  if (ierr /=0) then
     call yaml_warning('Problem in group creation, ierr:'//yaml_toa(ierr))
     call MPI_ABORT(base_comm,1,ierr)
  end if
  do i=0,ngroup-1
     !define the new groups and thread_id
     do j=0,group_size-1
        group_list(j+1)=i+j*ngroup
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

!i_all=-product(shape(group_list ))*kind(group_list )
! deallocate(group_list,stat=i_stat)
  call f_free(group_list)
end subroutine create_group_comm1

  !> Create a communicator between proc of same rank between the taskgroups.
  subroutine create_rank_comm(group_comm, rank_comm)
    use yaml_output
    implicit none
    integer, intent(in) :: group_comm
    integer, intent(out) :: rank_comm
    !local variables
    character(len=*), parameter :: subname='create_group_master'
    integer :: iproc_group, nproc, nproc_group, ngroups
    integer :: ierr, i, j
    integer, dimension(:), allocatable :: lrank, ids

    call mpi_comm_rank(group_comm, iproc_group, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
    call mpi_comm_size(group_comm, nproc_group, ierr)
    ngroups = nproc / nproc_group

    ! Put in lrank the group rank of each process, indexed by global iproc.
!   allocate(lrank(nproc+ndebug), stat = i_stat)
    lrank = f_malloc(nproc,id='lrank')
    call mpi_allgather(iproc_group, 1, MPI_INTEGER, lrank, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    ! Put in ids, the global iproc of each process that share the same group iproc.
!   allocate(ids(ngroups+ndebug), stat = i_stat)
    ids = f_malloc(ngroups,id='ids')
    j = 1
    do i = 1, nproc
       if (lrank(i) == iproc_group) then
          ids(j) = i - 1
          j = j + 1
       end if
    end do
!  i_all=-product(shape(lrank ))*kind(lrank )
!   deallocate(lrank,stat=i_stat)
    call f_free(lrank)

!!$    call mpi_comm_rank(MPI_COMM_WORLD, iproc_group, ierr)
!!$    write(*,*) iproc_group, "->", ids
    
    ! Create a new comminucator for the list of ids.
    call create_group_comm(MPI_COMM_WORLD, ngroups, ids, rank_comm)
!  i_all=-product(shape(ids ))*kind(ids )
!   deallocate(ids,stat=i_stat)
    call f_free(ids)
  END SUBROUTINE create_rank_comm


  !interface for MPI_ALLGATHERV operations
  subroutine mpiallgatherv_double(buffer,counts,displs,me,mpi_comm,ierr)
    implicit none
    integer, dimension(:), intent(in) :: counts, displs
    integer, intent(in) :: mpi_comm, me
    real(kind=8), intent(inout) :: buffer
    integer, intent(out) :: ierr
#ifdef HAVE_MPI2
    !case with MPI_IN_PLACE
    call MPI_ALLGATHERV(MPI_IN_PLACE,counts(me),MPI_DOUBLE_PRECISION,&
         buffer,counts,displs,MPI_DOUBLE_PRECISION,mpi_comm,ierr)
#else
    !local variables
    real(kind=8), dimension(:), allocatable :: copybuf

    !Here we have a performence penalty by copying all buffer, instead of
    !just the send part, but I don't see how to get buffer(displs(me))
    copybuf = f_malloc(sum(counts),id='copybuf')

    call dcopy(sum(counts),buffer,1,copybuf,1) 
    ierr=0 !put just for MPIfake compatibility
    call MPI_ALLGATHERV(copybuf(1+displs(me+1)),counts(me+1),MPI_DOUBLE_PRECISION,&
         buffer,counts,displs,MPI_DOUBLE_PRECISION,mpi_comm,ierr)

    call f_free(copybuf)
#endif
    if (ierr /=0) stop 'MPIALLGATHERV_DBL'
  end subroutine mpiallgatherv_double

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
    !integer :: i_all,i_stat
    integer, dimension(:), allocatable :: copybuf

    !case without mpi_in_place
!   allocate(copybuf(ntot+ndebug),stat=i_stat)
    copybuf = f_malloc(ntot,id='copybuf')

    !not appropriate for integers, to be seen if it works
    call scopy(ntot,buffer,1,copybuf,1) 
    ierr=0 !put just for MPIfake compatibility
    call MPI_ALLREDUCE(copybuf,buffer,ntot,&
         MPI_INTEGER,mpi_op,mpi_comm,ierr)

!  i_all=-product(shape(copybuf))*kind(copybuf)
!   deallocate(copybuf,stat=i_stat)
    call f_free(copybuf)
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
    !integer :: i_all,i_stat
    real(kind=4), dimension(:), allocatable :: copybuf

    !case without mpi_in_place
!   allocate(copybuf(ntot+ndebug),stat=i_stat)
    copybuf = f_malloc(ntot,id='copybuf')

    call scopy(ntot,buffer,1,copybuf,1) 
    ierr=0 !put just for MPIfake compatibility
    call MPI_ALLREDUCE(copybuf,buffer,ntot,&
         MPI_REAL,mpi_op,mpi_comm,ierr)

!  i_all=-product(shape(copybuf))*kind(copybuf)
!   deallocate(copybuf,stat=i_stat)
    call f_free(copybuf)
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
    !integer :: i_all,i_stat
    real(kind=8), dimension(:), allocatable :: copybuf

    !case without mpi_in_place
!   allocate(copybuf(ntot+ndebug),stat=i_stat)
    copybuf = f_malloc(ntot,id='copybuf')

    call dcopy(ntot,buffer,1,copybuf,1) 
    ierr=0 !put just for MPIfake compatibility
    call MPI_ALLREDUCE(copybuf,buffer,ntot,&
         MPI_DOUBLE_PRECISION,mpi_op,mpi_comm,ierr)

!  i_all=-product(shape(copybuf))*kind(copybuf)
!   deallocate(copybuf,stat=i_stat)
    call f_free(copybuf)
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
    !integer :: i_all,i_stat
    real(kind=8), dimension(:), allocatable :: copybuf

    !case without mpi_in_place
!   allocate(copybuf(ntot+ndebug),stat=i_stat)
    copybuf = f_malloc(ntot,id='copybuf')

    call dcopy(ntot,buffer,1,copybuf,1) 
    ierr=0 !put just for MPIfake compatibility
    call MPI_ALLREDUCE(copybuf,buffer,ntot,&
         MPI_DOUBLE_PRECISION,mpi_op,mpi_comm,ierr)

!  i_all=-product(shape(copybuf))*kind(copybuf)
!   deallocate(copybuf,stat=i_stat)
    call f_free(copybuf)
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
    !integer :: i_all,i_stat
    real(kind=8), dimension(:), allocatable :: copybuf

    !case without mpi_in_place
!   allocate(copybuf(ntot+ndebug),stat=i_stat)
    copybuf = f_malloc(ntot,id='copybuf')

    call dcopy(ntot,buffer,1,copybuf,1) 
    ierr=0 !put just for MPIfake compatibility
    call MPI_ALLREDUCE(copybuf,buffer,ntot,&
         MPI_DOUBLE_PRECISION,mpi_op,mpi_comm,ierr)

!  i_all=-product(shape(copybuf))*kind(copybuf)
!   deallocate(copybuf,stat=i_stat)
    call f_free(copybuf)
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
    !integer :: i_all,i_stat
    logical, dimension(:), allocatable :: copybuf

    !case without mpi_in_place
!   allocate(copybuf(ntot+ndebug),stat=i_stat)
    copybuf = f_malloc(ntot,id='copybuf')

    !not appropriate for logical, to be seen if it works
    call scopy(ntot,buffer,1,copybuf,1) 
    ierr=0 !put just for MPIfake compatibility
    call MPI_ALLREDUCE(copybuf,buffer,ntot,&
         MPI_LOGICAL,mpi_op,mpi_comm,ierr)

!  i_all=-product(shape(copybuf))*kind(copybuf)
!   deallocate(copybuf,stat=i_stat)
    call f_free(copybuf)
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
  integer :: idummy
  write(*,*)'BigDFT_open_nesting is not active!'
  !call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)
  stop
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
  integer :: idummy
  write(*,*)'BigDFT_close_nesting is not active!'
  stop
  !call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)
  idummy=num_threads
#endif
end subroutine bigdft_close_nesting
