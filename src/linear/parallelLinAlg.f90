!! ATTENTION: This works only if the matrices have the same sizes for all processes!!
subroutine dgemm_parallel(iproc, nproc, blocksize, comm, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
use module_base
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, blocksize, comm, m, n, k, lda, ldb, ldc
character(len=1),intent(in):: transa, transb
real(8),intent(in):: alpha, beta
real(8),dimension(lda,k),intent(in):: a
real(8),dimension(ldb,n),intent(in):: b
real(8),dimension(ldc,n),intent(out):: c

! Local variables
integer:: ierr, i, j, istat, iall, ii, ii1, ii2, mbrow, mbcol, nproc_scalapack, nprocrow, nproccol
integer:: context, irow, icol, numroc, lnrow, lncol, jproc, info
integer:: lnrow_a, lncol_a, lnrow_b, lncol_b, lnrow_c, lncol_c
real(8):: tt1, tt2
real(8),dimension(:,:),allocatable:: la, lb, lc
integer,dimension(9):: desc_lc, desc_la, desc_lb
character(len=*),parameter:: subname='dgemm_parallel'

  ! Block size for scalapack
  mbrow=blocksize
  mbcol=blocksize
  
  ! Number of processes that will be involved in the calculation
  tt1=dble(m)/dble(mbrow)
  tt2=dble(n)/dble(mbcol)
  ii1=ceiling(tt1)
  ii2=ceiling(tt2)
  nproc_scalapack = min(ii1*ii2,nproc)
  !if(iproc==0) write(*,'(a,i0,a)') 'scalapack will use ',nproc_scalapack,' processes.'
  
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
  
  
  ! Initialize blacs context,
  call blacs_get(-1, 0, context)
  call blacs_gridinit(context, 'r', nprocrow, nproccol )
  call blacs_gridinfo(context,nprocrow, nproccol, irow, icol)
  
  ! Initialize the result c to zero. For processes participating in the calculation, 
  ! c will be partially (only at the position that process was working on) overwritten with the result. 
  ! At the end we can the make an allreduce to get the correct result on all processes.
  if(irow==-1) c=0.d0
  
  ! Only execute this part if this process has a part of the matrix to work on. 
  processIf: if(irow/=-1) then

      ! Determine the size of the local matrix la (lnrow_a x lncol_a):
      lnrow_a = numroc(m, mbrow, irow, 0, nprocrow)
      lncol_a = numroc(k, mbcol, icol, 0, nproccol)
      !write(*,'(a,i0,a,i0,a,i0)') 'iproc ',iproc,' will have a local a of size ',lnrow_a,' x ',lncol_a

      ! Determine the size of the local matrix lb (lnrow_b x lncol_b):
      lnrow_b = numroc(k, mbrow, irow, 0, nprocrow)
      lncol_b = numroc(n, mbcol, icol, 0, nproccol)
      !write(*,'(a,i0,a,i0,a,i0)') 'iproc ',iproc,' will have a local b of size ',lnrow_b,' x ',lncol_b

      ! Determine the size of the local matrix lc (lnrow_c x lncol_c):
      lnrow_c = numroc(m, mbrow, irow, 0, nprocrow)
      lncol_c = numroc(n, mbcol, icol, 0, nproccol)
      !write(*,'(a,i0,a,i0,a,i0)') 'iproc ',iproc,' will have a local c of size ',lnrow_c,' x ',lncol_c
  
      ! Initialize descriptor arrays.
      call descinit(desc_la, m, k, mbrow, mbcol, 0, 0, context, lnrow_a, info)
      call descinit(desc_lb, k, n, mbrow, mbcol, 0, 0, context, lnrow_b, info)
      call descinit(desc_lc, m, n, mbrow, mbcol, 0, 0, context, lnrow_c, info)
  
      ! Allocate the local arrays
      allocate(la(lnrow_a,lncol_a), stat=istat)
      call memocc(istat, la, 'la', subname)
      allocate(lb(lnrow_b,lncol_b), stat=istat)
      call memocc(istat, lb, 'lb', subname)
      allocate(lc(lnrow_c,lncol_c), stat=istat)
      call memocc(istat, lc, 'lc', subname)
  
      ! Copy the global array a to the local array la.
      ! The same for b and lb, cpectively.
      do i=1,k
          do j=1,m
              call pdelset(la(1,1), j, i, desc_la, a(j,i))
          end do
      end do
      do i=1,n
          do j=1,k
              call pdelset(lb(1,1), j, i, desc_lb, b(j,i))
          end do
      end do
  
  
      ! Do the matrix matrix multiplication.
      call pdgemm(transa, transb, m, n, k, 1.d0, la(1,1), 1, 1, desc_la, lb(1,1), 1, 1, &
                  desc_lb, 0.d0, lc(1,1), 1, 1, desc_lc)
  
  
      ! Put the local result lc to the global result c.
      do i=1,n
          do j=1,m
              call pdelset2(c(j,i), lc(1,1), j, i, desc_lc, 0.d0)
          end do
      end do
  
      ! Deallocate the local arrays.
      iall=-product(shape(la))*kind(la)
      deallocate(la, stat=istat)
      call memocc(istat, iall, 'la', subname)
  
      iall=-product(shape(lb))*kind(lb)
      deallocate(lb, stat=istat)
      call memocc(istat, iall, 'lb', subname)
  
      iall=-product(shape(lc))*kind(lc)
      deallocate(lc, stat=istat)
      call memocc(istat, iall, 'lc', subname)
  
  end if processIf
  
  
  ! Gather the result on all processes.
  call mpiallred(c(1,1), m*n, mpi_sum, comm, ierr)


end subroutine dgemm_parallel





!! ATTENTION: This works only if the matrices have the same sizes for all processes!!
subroutine dsymm_parallel(iproc, nproc, blocksize, comm, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
use module_base
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, blocksize, comm, m, n, lda, ldb, ldc
character(len=1),intent(in):: side, uplo
real(8),intent(in):: alpha, beta
real(8),dimension(lda,m),intent(in):: a
real(8),dimension(ldb,n),intent(in):: b
real(8),dimension(ldc,n),intent(out):: c

! Local variables
integer:: ierr, i, j, istat, iall, ii, ii1, ii2, mbrow, mbcol, nproc_scalapack, nprocrow, nproccol
integer:: context, irow, icol, numroc, lnrow, lncol, jproc, info
integer:: lnrow_a, lncol_a, lnrow_b, lncol_b, lnrow_c, lncol_c
real(8):: tt1, tt2
real(8),dimension(:,:),allocatable:: la, lb, lc
integer,dimension(9):: desc_lc, desc_la, desc_lb
character(len=*),parameter:: subname='dgemm_parallel'

  ! Block size for scalapack
  mbrow=blocksize
  mbcol=blocksize
  
  ! Number of processes that will be involved in the calculation
  tt1=dble(m)/dble(mbrow)
  tt2=dble(n)/dble(mbcol)
  ii1=ceiling(tt1)
  ii2=ceiling(tt2)
  nproc_scalapack = min(ii1*ii2,nproc)
  !if(iproc==0) write(*,'(a,i0,a)') 'scalapack will use ',nproc_scalapack,' processes.'
  
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
  
  
  ! Initialize blacs context,
  call blacs_get(-1, 0, context)
  call blacs_gridinit(context, 'r', nprocrow, nproccol )
  call blacs_gridinfo(context,nprocrow, nproccol, irow, icol)
  
  ! Initialize the result c to zero. For processes participating in the calculation, 
  ! c will be partially (only at the position that process was working on) overwritten with the result. 
  ! At the end we can the make an allreduce to get the correct result on all processes.
  if(irow==-1) c=0.d0
  
  ! Only execute this part if this process has a part of the matrix to work on. 
  processIf: if(irow/=-1) then

      ! Determine the size of the local matrix la (lnrow_a x lncol_a):
      lnrow_a = numroc(m, mbrow, irow, 0, nprocrow)
      lncol_a = numroc(m, mbcol, icol, 0, nproccol)
      !write(*,'(a,i0,a,i0,a,i0)') 'iproc ',iproc,' will have a local a of size ',lnrow_a,' x ',lncol_a

      ! Determine the size of the local matrix lb (lnrow_b x lncol_b):
      lnrow_b = numroc(m, mbrow, irow, 0, nprocrow)
      lncol_b = numroc(n, mbcol, icol, 0, nproccol)
      !write(*,'(a,i0,a,i0,a,i0)') 'iproc ',iproc,' will have a local b of size ',lnrow_b,' x ',lncol_b

      ! Determine the size of the local matrix lc (lnrow_c x lncol_c):
      lnrow_c = numroc(m, mbrow, irow, 0, nprocrow)
      lncol_c = numroc(n, mbcol, icol, 0, nproccol)
      !write(*,'(a,i0,a,i0,a,i0)') 'iproc ',iproc,' will have a local c of size ',lnrow_c,' x ',lncol_c
  
      ! Initialize descriptor arrays.
      call descinit(desc_la, m, m, mbrow, mbcol, 0, 0, context, lnrow_a, info)
      call descinit(desc_lb, m, n, mbrow, mbcol, 0, 0, context, lnrow_b, info)
      call descinit(desc_lc, m, n, mbrow, mbcol, 0, 0, context, lnrow_c, info)
  
      ! Allocate the local arrays
      allocate(la(lnrow_a,lncol_a), stat=istat)
      call memocc(istat, la, 'la', subname)
      allocate(lb(lnrow_b,lncol_b), stat=istat)
      call memocc(istat, lb, 'lb', subname)
      allocate(lc(lnrow_c,lncol_c), stat=istat)
      call memocc(istat, lc, 'lc', subname)
  
      ! Copy the global array a to the local array la.
      ! The same for b and lb, cpectively.
      do i=1,m
          do j=1,m
              call pdelset(la(1,1), j, i, desc_la, a(j,i))
          end do
      end do
      do i=1,n
          do j=1,m
              call pdelset(lb(1,1), j, i, desc_lb, b(j,i))
          end do
      end do
  
  
      ! Do the matrix matrix multiplication.
      call pdsymm(side, uplo, m, n, alpha, la(1,1), 1, 1, desc_la, lb(1,1), 1, 1, &
                  desc_lb, beta, lc(1,1), 1, 1, desc_lc)
  
  
      ! Put the local result lc to the global result c.
      do i=1,n
          do j=1,m
              call pdelset2(c(j,i), lc(1,1), j, i, desc_lc, 0.d0)
          end do
      end do
  
      ! Deallocate the local arrays.
      iall=-product(shape(la))*kind(la)
      deallocate(la, stat=istat)
      call memocc(istat, iall, 'la', subname)
  
      iall=-product(shape(lb))*kind(lb)
      deallocate(lb, stat=istat)
      call memocc(istat, iall, 'lb', subname)
  
      iall=-product(shape(lc))*kind(lc)
      deallocate(lc, stat=istat)
      call memocc(istat, iall, 'lc', subname)
  
  end if processIf
  
  
  ! Gather the result on all processes.
  call mpiallred(c(1,1), m*n, mpi_sum, comm, ierr)


end subroutine dsymm_parallel


subroutine diagonalizeHamiltonianParallel(iproc, nproc, norb, ham, ovrlp, eval)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, norb
  real(8),dimension(norb,norb),intent(inout):: ham, ovrlp
  real(8),dimension(norb),intent(out):: eval
  
  ! Local variables
  integer:: ierr, mbrow, mbcol, i, j, istat, lwork, info, ii1, ii2, nproc_scalapack, iall
  integer:: nprocrow, nproccol, context, irow, icol, lnrow, lncol, numroc, jproc, liwork, neval_found, neval_computed
  real(8):: tt1, tt2
  real(8),dimension(:,:),allocatable:: lmat, loverlap, levec
  real(8),dimension(:),allocatable:: work, gap
  integer,dimension(9):: desc_levec, desc_lmat, desc_loverlap
  integer,dimension(:),allocatable:: iwork, ifail, icluster
  character(len=*),parameter:: subname='diagonalizeHamiltonianParallel'
  
  
  
  
  ! Block size for scalapack
  mbrow=64
  mbcol=64
  
  ! Number of processes that will be involved in the calculation
  tt1=dble(norb)/dble(mbrow)
  tt2=dble(norb)/dble(mbcol)
  ii1=ceiling(tt1)
  ii2=ceiling(tt2)
  nproc_scalapack = min(ii1*ii2,nproc)
  !nproc_scalapack = nproc
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
  if(irow==-1) ham=0.d0
  
  ! Everything that follows is only done if the current process is part of the grid.
  processIf: if(irow/=-1) then
      ! Determine the size of the matrix (lnrow x lncol):
      lnrow = numroc(norb, mbrow, irow, 0, nprocrow)
      lncol = numroc(norb, mbcol, icol, 0, nproccol)
      write(*,'(a,i0,a,i0,a,i0)') 'iproc ',iproc,' will have a local matrix of size ',lnrow,' x ',lncol
  
      ! Initialize descriptor arrays.
      call descinit(desc_lmat, norb, norb, mbrow, mbcol, 0, 0, context, lnrow, info)
      call descinit(desc_loverlap, norb, norb, mbrow, mbcol, 0, 0, context, lnrow, info)
      call descinit(desc_levec, norb, norb, mbrow, mbcol, 0, 0, context, lnrow, info)
  
      ! Allocate the local array lmat
      allocate(lmat(lnrow,lncol), stat=istat)
      call memocc(istat, lmat, 'lmat', subname)
      allocate(loverlap(lnrow,lncol), stat=istat)
      call memocc(istat, loverlap, 'loverlap', subname)
  
      ! Copy the global array mat to the local array lmat.
      ! The same for loverlap and overlap, respectively.
      !call dcopy(norb**2, ham(1,1), 1, mat(1,1), 1)
      !call dcopy(norb**2, ovrlp(1,1), 1, overlap(1,1), 1)
      do i=1,norb
          do j=1,norb
              call pdelset(lmat(1,1), j, i, desc_lmat, ham(j,i))
              call pdelset(loverlap(1,1), j, i, desc_loverlap, ovrlp(j,i))
          end do
      end do
  
  
      ! Solve the generalized eigenvalue problem.
      allocate(levec(lnrow,lncol), stat=istat)
      call memocc(istat, levec, 'levec', subname)
      allocate(ifail(norb), stat=istat)
      call memocc(istat, ifail, 'ifail', subname)
      allocate(icluster(2*nprocrow*nproccol), stat=istat)
      call memocc(istat, icluster, 'icluster', subname)
      allocate(gap(nprocrow*nproccol), stat=istat)
      call memocc(istat, gap, 'gap', subname)
  
      ! workspace query
      lwork=-1
      liwork=-1
      allocate(work(1), stat=istat)
      call memocc(istat, work, 'work', subname)
      allocate(iwork(1), stat=istat) ; if(istat/=0) stop 'ERROR in allocating'
      call memocc(istat, iwork, 'iwork', subname)
      call pdsygvx(1, 'v', 'a', 'l', norb, lmat(1,1), 1, 1, desc_lmat, loverlap(1,1), 1, 1, &
                   desc_loverlap, 0.d0, 1.d0, 0, 1, -1.d0, neval_found, neval_computed, eval(1), &
                   -1.d0, levec(1,1), 1, 1, desc_levec, work, lwork, iwork, liwork, &
                   ifail, icluster, gap, info)
      lwork=ceiling(work(1))
      liwork=iwork(1)
      !write(*,*) 'iproc, lwork, liwork', iproc, lwork, liwork
      iall=-product(shape(work))*kind(work)
      deallocate(work, stat=istat)
      call memocc(istat, iall, 'work', subname)
      iall=-product(shape(iwork))*kind(iwork)
      deallocate(iwork, stat=istat)
      call memocc(istat, iall, 'iwork', subname)
  
      allocate(work(lwork), stat=istat)
      call memocc(istat, work, 'work', subname)
      allocate(iwork(liwork), stat=istat)
      call memocc(istat, iwork, 'iwork', subname)
  
      call pdsygvx(1, 'v', 'a', 'l', norb, lmat(1,1), 1, 1, desc_lmat, loverlap(1,1), 1, 1, &
                   desc_loverlap, 0.d0, 1.d0, 0, 1, -1.d0, neval_found, neval_computed, eval(1), &
                   -1.d0, levec(1,1), 1, 1, desc_levec, work, lwork, iwork, liwork, &
                   ifail, icluster, gap, info)
  
      ! Gather together the eigenvectors from all processes and store them in mat.
      do i=1,norb
          do j=1,norb
              call pdelset2(ham(j,i), levec(1,1), j, i, desc_lmat, 0.d0)
          end do
      end do
  
  
      iall=-product(shape(lmat))*kind(lmat)
      deallocate(lmat, stat=istat)
      call memocc(istat, iall, 'lmat', subname)
  
      iall=-product(shape(levec))*kind(levec)
      deallocate(levec, stat=istat)
      call memocc(istat, iall, 'levec', subname)
  
      iall=-product(shape(loverlap))*kind(loverlap)
      deallocate(loverlap, stat=istat)
      call memocc(istat, iall, 'loverlap', subname)
  
      iall=-product(shape(work))*kind(work)
      deallocate(work, stat=istat)
      call memocc(istat, iall, 'work', subname)
  
      iall=-product(shape(iwork))*kind(iwork)
      deallocate(iwork, stat=istat)
      call memocc(istat, iall, 'iwork', subname)
  
      iall=-product(shape(ifail))*kind(ifail)
      deallocate(ifail, stat=istat)
      call memocc(istat, iall, 'ifail', subname)
  
      iall=-product(shape(icluster))*kind(icluster)
      deallocate(icluster, stat=istat)
      call memocc(istat, iall, 'icluster', subname)
  
      iall=-product(shape(gap))*kind(gap)
      deallocate(gap, stat=istat)
      call memocc(istat, iall, 'gap', subname)
  
  end if processIF
  
  ! Gather the eigenvectors on all processes.
  call mpiallred(ham(1,1), norb**2, mpi_sum, mpi_comm_world, ierr)
  
  ! Broadcast the eigenvalues if required. If nproc_scalapack==nproc, then all processes
  ! diagonalized the matrix and therefore have the eigenvalues.
  if(nproc_scalapack/=nproc) then
      call mpi_bcast(eval(1), norb, mpi_double_precision, 0, mpi_comm_world, ierr)
  end if


end subroutine diagonalizeHamiltonianParallel
