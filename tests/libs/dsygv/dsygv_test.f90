!> @file
!!  Program to test the dsygv routine
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Test the dsygv routine
program dsygv_test
  !use dsygv_interfaces
  implicit none
  include 'mpif.h'

  ! Variables
  integer :: iproc, nproc, ierr, n, istat, blocksize, nthread, ithread
  integer :: omp_get_num_threads
  real(kind=8) :: t1, t2
  real(kind=8),dimension(:),allocatable :: eval
  real(kind=8),dimension(:,:),allocatable :: A_init, S_init, A, S
  integer,parameter :: itype=1
  character(len=1),parameter :: jobz='v', uplo='l'

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, iproc, ierr)
  call mpi_comm_size(mpi_comm_world, nproc, ierr)

  nthread=1
  !$omp parallel shared(nthread)
  !$omp master
  !$nthread=omp_get_num_threads()
  !$omp end master
  !$omp end parallel

  if (iproc==0) write(*,'(1x,a,i0)') 'number of MPI processes: ',nproc
  if (iproc==0) write(*,'(1x,a,i0)') 'number of OpenMP threads: ',nthread

  if (iproc==0) write(*,'(1x,a)') '>>> reading parameters from file input.dat <<<'
  open(unit=1,file='input.dat')
  read(1,*) n
  read(1,*) blocksize
  close(unit=1)

  if (iproc==0) then
      write(*,'(1x,a,i0)') 'matrix size: ',n
      write(*,'(1x,a,i0)') 'scalapack block size: ',blocksize
  end if

  allocate(eval(n), stat=istat) ; if (istat/=0) stop 'ERROR in allocating'
  allocate(S_init(n,n), stat=istat) ; if (istat/=0) stop 'ERROR in allocating'
  allocate(A_init(n,n), stat=istat) ; if (istat/=0) stop 'ERROR in allocating'
  allocate(S(n,n), stat=istat) ; if (istat/=0) stop 'ERROR in allocating'
  allocate(A(n,n), stat=istat) ; if (istat/=0) stop 'ERROR in allocating'
  
  call init_matrices(n, A_init, S_init)

  call dcopy(n**2, A_init(1,1), 1, A(1,1), 1)
  call dcopy(n**2, S_init(1,1), 1, S(1,1), 1)

  call mpi_barrier(mpi_comm_world, ierr)
  t1=mpi_wtime()
  call dsygv_wrapper(itype, jobz, uplo, n, A, n, S, n, eval)
  call mpi_barrier(mpi_comm_world, ierr)
  t2=mpi_wtime()
  if(iproc==0) write(*,'(1x,a)') '============================================================='
  if(iproc==0) write(*,'(1x,a,es17.9,a,es10.4)') 'LAPACK: sum of eigenvalues',sum(eval),',  time=',t2-t1
  if(iproc==0) write(*,'(1x,a)') '============================================================='


  call dcopy(n**2, A_init(1,1), 1, A(1,1), 1)
  call dcopy(n**2, S_init(1,1), 1, S(1,1), 1)
  call mpi_barrier(mpi_comm_world, ierr)
  t1=mpi_wtime()
  call pdsygvx_wrapper(iproc, nproc, blocksize, mpi_comm_world, itype, jobz, uplo, n, A, n, S, n, eval)
  call mpi_barrier(mpi_comm_world, ierr)
  t2=mpi_wtime()
  if(iproc==0) write(*,'(1x,a)') '================================================================'
  if(iproc==0) write(*,'(1x,a,es17.9,a,es10.4)') 'SCALAPACK: sum of eigenvalues',sum(eval),',  time=',t2-t1
  if(iproc==0) write(*,'(1x,a)') '================================================================'

  deallocate(eval, stat=istat) ; if (istat/=0) stop 'ERROR in allocating'
  deallocate(S_init, stat=istat) ; if (istat/=0) stop 'ERROR in allocating'
  deallocate(A_init, stat=istat) ; if (istat/=0) stop 'ERROR in allocating'
  deallocate(S, stat=istat) ; if (istat/=0) stop 'ERROR in allocating'
  deallocate(A, stat=istat) ; if (istat/=0) stop 'ERROR in allocating'
  
  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_finalize(ierr)

end program dsygv_test
