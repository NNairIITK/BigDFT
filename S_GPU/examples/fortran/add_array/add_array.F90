program mull_array
implicit none
include 'mpif.h'
!include 's_gpu.h'

integer, parameter :: MASTER_TASK = 0
integer, parameter :: MATRIX_W = 10
! The height of the matrix must be a multiple of "number of MPI tasks - 1"
integer, parameter :: MATRIX_H  = 10
integer, parameter :: MAX_TASKS = 128

! mpi variables
integer :: rank, nb_tasks
integer :: stats ( MPI_STATUS_SIZE, MAX_TASKS ), ireq ( MAX_TASKS ) 
integer :: stat ( MPI_STATUS_SIZE ), ierr

! s_gpu variables
integer :: usegpu
logical :: gpushare, initerror
real ( kind = 8 ) :: stream_ptr

! matrices variables
!   Master side
integer :: nb_rows_by_task, remainder , i, j
double precision, dimension ( MATRIX_H * MATRIX_W ) :: matrix_a, matrix_b, matrix_result
!   Client side
! host pointers
double precision, dimension ( : ), allocatable :: a_h, b_h, c_h
! pinned pointers
real ( kind = 8 ) :: a_p, b_p, c_p
! device pointers
real ( kind = 8 ) :: a_d, b_d, c_d

! MPI initialization
call MPI_INIT ( ierr )
call MPI_COMM_RANK ( MPI_COMM_WORLD, rank, ierr )
call MPI_COMM_SIZE ( MPI_COMM_WORLD, nb_tasks, ierr )

gpushare = .true.
! S_GPU initialization
! This initialization function must be called by all the MPI task
! It must be called befaore any other s_gpu function
call sg_init ( gpushare, usegpu, rank, initerror )

if ( nb_tasks > MAX_TASKS ) then
  print *, "Too many MPI tasks: ", nb_tasks
  stop
endif
if ( modulo(MATRIX_H, nb_tasks - 1 ) /= 0 ) then
  print *, "The height of the matrix must be a multiple of number of MPI tasks - 1"
  print *, "matrix height = ", MATRIX_H
  print *, "number of MPI tasks = ", nb_tasks
  stop
endif

! Splitting of the matrixes
nb_rows_by_task = MATRIX_H / ( nb_tasks - 1 )

! Master node computation
if ( rank == MASTER_TASK ) then

! Initialization of the matrices
  do i = 1, MATRIX_W
    do j = 1, MATRIX_H
      matrix_a ( ( i - 1 ) * MATRIX_W + j ) = i
      matrix_b ( ( i - 1 ) * MATRIX_W + j ) = j
    end do
  end do

! Print matrix
!  print *, "====="
!  print *, "Node ", rank
!  do i = 1, MATRIX_H
!    do j = 1, MATRIX_W
!      write ( *, '(F3.0,A1)', advance = 'no' ) matrix_a ( ( i - 1 ) * MATRIX_W + j ), " "
!    end do
!    print *, " "
!  end do
!  print *, "====="

! Send intervals of arrays A and B to worker processes
  do i = 1, nb_tasks - 1
    call MPI_ISEND ( matrix_a ( ( i - 1 ) * nb_rows_by_task * MATRIX_W + 1 ), &
    nb_rows_by_task * MATRIX_W, &
    MPI_DOUBLE_PRECISION, &
    i, &
    i, &
    MPI_COMM_WORLD, &
    ireq ( i ), &
    ierr )
  end do
  do i = 1, nb_tasks - 1
    call MPI_ISEND ( matrix_b ( ( i - 1 ) * nb_rows_by_task * MATRIX_W + 1 ), &
    nb_rows_by_task * MATRIX_W, &
    MPI_DOUBLE_PRECISION, &
    i, &
    i, &
    MPI_COMM_WORLD, &
    ireq ( nb_tasks - 1 + i ), &
    ierr )
  end do
  call MPI_WAITALL ( nb_tasks * 2 - 2, ireq, stats, ierr )

! Receive all the result intervals
  do i = 1, nb_tasks - 1
    call MPI_IRECV ( matrix_result ( ( i - 1 ) * nb_rows_by_task * MATRIX_W + 1 ), &
    nb_rows_by_task * MATRIX_W, &
    MPI_DOUBLE_PRECISION, &
    i, &
    i, &
    MPI_COMM_WORLD, &
    ireq ( i ), &
    ierr )
  end do
  call MPI_WAITALL ( nb_tasks - 1 , ireq, stats, ierr )

! Print results
!  print *, "====="
!  print *, "Results: "
!  do i = 1, MATRIX_W
!  do j = 1, MATRIX_H
!  write ( *, * ) matrix_result ( ( i - 1 ) * MATRIX_W + j ), " "
!  end do
!  print *, " "
!  end do
!  print *, "====="
else
  ! Client node computation

  ! Host pointers for arrays A, B and result
  ! Allocate memory to receive the interval on this MPI node
  allocate ( a_h ( nb_rows_by_task * MATRIX_W ) )
  allocate ( b_h ( nb_rows_by_task * MATRIX_W ) )
  allocate ( c_h ( nb_rows_by_task * MATRIX_W ) )

! Receive the interval of the matrices A and B
  call MPI_RECV ( a_h, nb_rows_by_task * MATRIX_W, &
  MPI_DOUBLE_PRECISION, &
  MASTER_TASK, &
  rank, &
  MPI_COMM_WORLD, &
  stat, &
  ierr )
  call MPI_RECV ( b_h, nb_rows_by_task * MATRIX_W, &
  MPI_DOUBLE_PRECISION, &
  MASTER_TASK, &
  rank, &
  MPI_COMM_WORLD, &
  stat, &
  ierr )

! Print interval
!  print *, "====="
!  print *, "Node ", rank
!  do i = 1, nb_rows_by_task
!    do j = 1, MATRIX_W
!      write ( *, * ) a_h ( ( i - 1 ) * MATRIX_W + j ), " "
!    end do
!    print *, " "
!  end do
!  print *, "====="

! In order to benefit the asynchroneous mem copy CPU -> GPU (and GPU -> CPU), the host memory must be a special memory
! called pinned memory. Before any transfer to the GPU, the host memory is first copied to a pinned memory area.
! Allocate pinned memory
! This operation is immediate and blocking
  call sg_cpu_malloc_pinned ( a_p, nb_rows_by_task * MATRIX_W, 8, ierr )
  call sg_cpu_malloc_pinned ( b_p, nb_rows_by_task * MATRIX_W, 8, ierr )
  call sg_cpu_malloc_pinned ( c_p, nb_rows_by_task * MATRIX_W, 8, ierr )

  if ( ierr == 1 ) then
    print *, "Problem when trying to allocate pinned host memory"
    stop
  end if

! Allocate memory on the GPU device
! This operation is immediate and blocking
  call sg_gpu_malloc ( a_d, nb_rows_by_task * MATRIX_W, 8, ierr )
  call sg_gpu_malloc ( b_d, nb_rows_by_task * MATRIX_W, 8, ierr )
  call sg_gpu_malloc ( c_d, nb_rows_by_task * MATRIX_W, 8, ierr )

  if ( ierr == 1 ) then
    print *, "Problem when trying to allocate GPU memory"
    stop
  end if

! Create the stream
! The GPU operations will be sent on this stream
  call sg_create_stream ( stream_ptr )

! Send the second half of the part of the matrices A and B on the GPU for computation
! This operation is added in the stream. The execution is delayed.
  call sg_gpu_send_mem ( a_d, &
  a_h ( nb_rows_by_task * MATRIX_W / 2 ), &
  a_p, &
  nb_rows_by_task * MATRIX_W / 2 + 1, &
  8, &
  stream_ptr, &
  ierr )
  call sg_gpu_send_mem ( b_d, &
  b_h ( nb_rows_by_task * MATRIX_W / 2 ), &
  b_p, nb_rows_by_task * MATRIX_W / 2 + 1, &
  8, &
  stream_ptr, &
  ierr )

! This function adds the callback function in the stream
! This indirection is needed because the callback input is a C struct
  call binding_callback_mull ( a_d, b_d, c_d, &
  nb_rows_by_task * MATRIX_W / 2 + 1, stream_ptr )

! Retrieve results from device and store it in host array
! This operation is added in the stream. The execution is delayed.
  call sg_gpu_recv_mem ( c_h ( nb_rows_by_task * MATRIX_W / 2 ), &
  c_d, &
  c_p, &
  nb_rows_by_task * MATRIX_W / 2 + 1, &
  8, &
  stream_ptr, &
  ierr )

! The operations added to the stream previously are now executed asynchroneously
! If several streams have been created, they would be executed concurently
  call sg_exec_all_streams()

! In the same time, the MPI node compute the other half of the matrix addition
  do i = 1, nb_rows_by_task / 2
    do j = 1, MATRIX_W
      c_h ( ( i - 1 ) * MATRIX_W + j ) = &
        a_h ( ( i - 1 ) * MATRIX_W + j ) * b_h ( ( i - 1 ) * MATRIX_W + j )
    end do
  end do

! The MPI node send the result of its interval to the master node
  call MPI_SEND ( c_h, nb_rows_by_task * MATRIX_W, MPI_DOUBLE_PRECISION, &
  MASTER_TASK, rank, MPI_COMM_WORLD, ierr )

endif

call MPI_FINALIZE ( ierr )


end program mull_array
