!> @file
!!  Routines to wrap the dsygv call for the test of dsygv (Lapack or Scalapak).
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Wrap the dsygv call (dsygv test).
subroutine dsygv_wrapper(itype, jobz, uplo, n, A, ldA, S, ldS, eval)
  implicit none

  ! Calling arguments
  integer,intent(in) :: itype, n, ldA, ldS
  character(len=1),intent(in) :: jobz, uplo
  real(kind=8),dimension(ldA,n),intent(inout) :: A
  real(kind=8),dimension(ldS,n),intent(inout) :: S
  real(kind=8),dimension(n),intent(inout) :: eval

  ! Local variables
  integer :: lwork, info, istat
  real(kind=8),dimension(:),allocatable :: work

  lwork=100*n
  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating'
  call dsygv(itype, jobz, uplo, n, A(1,1), n, S(1,1), n, eval(1), work, lwork, info)
  if(info/=0) stop 'ERROR in dsygv'
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating'

end subroutine dsygv_wrapper


!> Routine wrapping the parallel dsygv (scalapack).
subroutine pdsygvx_wrapper(iproc, nproc, blocksize, comm, itype, jobz, uplo, n, a, lda, b, ldb, w)
  implicit none
  include 'mpif.h'
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, blocksize, comm, itype, n, lda, ldb
  character(len=1),intent(in) :: jobz, uplo
  real(kind=8),dimension(lda,n),intent(inout) :: a
  real(kind=8),dimension(ldb,n),intent(inout) :: b
  real(kind=8),dimension(n),intent(out) :: w
  
  ! Local variables
  integer :: ierr, mbrow, mbcol, i, j, istat, lwork, ii1, ii2, nproc_scalapack, iall, info
  integer :: nprocrow, nproccol, context, irow, icol, lnrow, lncol, numroc, liwork, nw_found, nw_computed
  real(kind=8) :: tt1, tt2
  real(kind=8),dimension(:,:),allocatable :: la, lb, lz
  real(kind=8),dimension(:),allocatable :: work, gap
  integer,dimension(9) :: desc_lz, desc_la, desc_lb
  integer,dimension(:),allocatable :: iwork, ifail, icluster
  
  
  ! Block size for scalapack
  mbrow=blocksize
  mbcol=blocksize
  
  ! Number of processes that will be involved in the calculation
  tt1=dble(n)/dble(mbrow)
  tt2=dble(n)/dble(mbcol)
  ii1=ceiling(tt1)
  ii2=ceiling(tt2)
  nproc_scalapack = min(ii1*ii2,nproc)
  if(iproc==0) write(*,'(a,i0,a)') 'scalapack will use ',nproc_scalapack,' processes.'
  
  ! process grid: number of processes per row and column
  tt1=sqrt(dble(nproc_scalapack))
  ii1=ceiling(tt1)
  do i=ii1,nproc_scalapack
      if(mod(nproc_scalapack,i)==0) then
          nprocrow=i
          exit
      end if
  end do
  nproccol=nproc_scalapack/nprocrow
  if(iproc==0) write(*,'(a,i0,a,i0,a)') 'calculation is done on process grid with dimension ',nprocrow,' x ',nproccol,'.'
  
  
  ! Initialize blacs context
  call blacs_get(-1, 0, context)
  call blacs_gridinit(context, 'r', nprocrow, nproccol )
  call blacs_gridinfo(context,nprocrow, nproccol, irow, icol)
  !write(*,*) 'iproc, irow, icol', iproc, irow, icol
  
  ! Initialize the matrix mat to zero for processes that don't do the calculation.
  ! For processes participating in the diagonalization, 
  ! it will be partially (only at the position that process was working on) overwritten with the result. 
  ! At the end we can the make an allreduce to get the correct result on all processes.
  if(irow==-1) call dscal(lda*n, 0.d0, a(1,1), 1)
  
  ! Everything that follows is only done if the current process is part of the grid.
  processIf: if(irow/=-1) then
      ! Determine the size of the matrix (lnrow x lncol):
      lnrow = max(numroc(n, mbrow, irow, 0, nprocrow),1)
      lncol = max(numroc(n, mbcol, icol, 0, nproccol),1)
      !write(*,'(a,i0,a,i0,a,i0)') 'iproc ',iproc,' will have a local matrix of size ',lnrow,' x ',lncol
  
      ! Initialize descriptor arrays.
      call descinit(desc_la, n, n, mbrow, mbcol, 0, 0, context, lnrow, info)
      call descinit(desc_lb, n, n, mbrow, mbcol, 0, 0, context, lnrow, info)
      call descinit(desc_lz, n, n, mbrow, mbcol, 0, 0, context, lnrow, info)
  
      ! Allocate the local array la
      allocate(la(lnrow,lncol), stat=istat)
      allocate(lb(lnrow,lncol), stat=istat)
  
      ! Copy the global array mat to the local array la.
      ! The same for lb and b, respectively.
      !call dcopy(n**2, a(1,1), 1, mat(1,1), 1)
      !call dcopy(n**2, b(1,1), 1, b(1,1), 1)
      do i=1,n
          do j=1,n
              call pdelset(la(1,1), j, i, desc_la, a(j,i))
              call pdelset(lb(1,1), j, i, desc_lb, b(j,i))
          end do
      end do
  
  
      ! Solve the generalized eigenvalue problem.
      allocate(lz(lnrow,lncol), stat=istat)
      allocate(ifail(n), stat=istat)
      allocate(icluster(2*nprocrow*nproccol), stat=istat)
      allocate(gap(nprocrow*nproccol), stat=istat)
  
      ! workspace query
      lwork=-1
      liwork=-1
      allocate(work(100), stat=istat)
      allocate(iwork(100), stat=istat)
      call pdsygvx(itype, jobz, 'a', uplo, n, la(1,1), 1, 1, desc_la, lb(1,1), 1, 1, &
                   desc_lb, 0.d0, 1.d0, 0, 1, -1.d0, nw_found, nw_computed, w(1), &
                   -1.d0, lz(1,1), 1, 1, desc_lz, work, lwork, iwork, liwork, &
                   ifail, icluster, gap, info)
      lwork=ceiling(work(1))
      lwork=lwork+n**2 !to be sure to have enough workspace, to be optimized later.
      liwork=iwork(1)
      liwork=liwork+n**2 !to be sure to have enough workspace, to be optimized later.
      !write(*,*) 'iproc, lwork, liwork', iproc, lwork, liwork
      deallocate(work, stat=istat)
      deallocate(iwork, stat=istat)
  
      allocate(work(lwork), stat=istat)
      allocate(iwork(liwork), stat=istat)
  
      call pdsygvx(1, 'v', 'a', 'l', n, la(1,1), 1, 1, desc_la, lb(1,1), 1, 1, &
                   desc_lb, 0.d0, 1.d0, 0, 1, -1.d0, nw_found, nw_computed, w(1), &
                   -1.d0, lz(1,1), 1, 1, desc_lz, work, lwork, iwork, liwork, &
                   ifail, icluster, gap, info)
      if(info/=0) then
          write(*,'(2(a,i0))') 'ERROR in pdsygvx on process ',iproc,', info=',info
      end if
  
      ! Gather together the eigenvectors from all processes and store them in mat.
      do i=1,n
          do j=1,n
              call pdelset2(a(j,i), lz(1,1), j, i, desc_la, 0.d0)
          end do
      end do
  
  
      deallocate(la, stat=istat)
  
      deallocate(lz, stat=istat)
  
      deallocate(lb, stat=istat)
  
      deallocate(work, stat=istat)
  
      deallocate(iwork, stat=istat)
  
      deallocate(ifail, stat=istat)
  
      deallocate(icluster, stat=istat)
  
      deallocate(gap, stat=istat)

      call blacs_gridexit(context)
  
  end if processIF
  
  ! Gather the eigenvectors on all processes. Use b as temprorary array.
  call dcopy(n**2, a(1,1), 1, b(1,1), 1)
  call mpi_allreduce(b(1,1), a(1,1), n**2, mpi_double_precision, mpi_sum, comm, ierr)
  
  ! Broadcast the eigenvalues if required. If nproc_scalapack==nproc, then all processes
  ! diagonalized the matrix and therefore have the eigenvalues.
  if(nproc_scalapack/=nproc) then
      call mpi_bcast(w(1), n, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_bcast(info, 1, mpi_integer, 0, mpi_comm_world, ierr)
  end if

end subroutine pdsygvx_wrapper
