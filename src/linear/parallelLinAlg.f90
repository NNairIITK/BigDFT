subroutine dgemm_parallel(iproc, nproc, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
use module_base
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, m, n, k, lda, ldb, ldc
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
  mbrow=4
  mbcol=4
  
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
  call mpiallred(c(1,1), m*n, mpi_sum, mpi_comm_world, ierr)


end subroutine dgemm_parallel
