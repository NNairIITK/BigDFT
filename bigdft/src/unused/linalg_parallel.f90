    !> ATTENTION: This works only if the matrices have the same sizes for all processes!!
    subroutine dsymm_parallel(iproc, nproc, blocksize, comm, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    use module_base
    implicit none
    
    ! Calling arguments
    integer,intent(in) :: iproc, nproc, blocksize, comm, m, n, lda, ldb, ldc
    character(len=1),intent(in) :: side, uplo
    real(kind=8),intent(in) :: alpha, beta
    real(kind=8),dimension(lda,m),intent(in) :: a
    real(kind=8),dimension(ldb,n),intent(in) :: b
    real(kind=8),dimension(ldc,n),intent(out) :: c
    
    ! Local variables
    integer :: ierr, i, j, istat, iall, ii1, ii2, mbrow, mbcol, nproc_scalapack, nprocrow, nproccol
    integer :: context, irow, icol, numroc, info
    integer :: lnrow_a, lncol_a, lnrow_b, lncol_b, lnrow_c, lncol_c
    real(kind=8) :: tt1, tt2
    real(kind=8),dimension(:,:),allocatable :: la, lb, lc
    integer,dimension(9) :: desc_lc, desc_la, desc_lb
    character(len=*),parameter :: subname='dgemm_parallel'
    
    
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
      !if(irow==-1) call to_zero(ldc*n, c(1,1))
      if(irow==-1) call vscal(ldc*n,0.0_wp, c(1,1),1)
      
      ! Only execute this part if this process has a part of the matrix to work on. 
      processIf: if(irow/=-1) then
    
          ! Determine the size of the local matrix la (lnrow_a x lncol_a):
          lnrow_a = max(numroc(m, mbrow, irow, 0, nprocrow),1)
          lncol_a = max(numroc(m, mbcol, icol, 0, nproccol),1)
          !write(*,'(a,i0,a,i0,a,i0)') 'iproc ',iproc,' will have a local a of size ',lnrow_a,' x ',lncol_a
    
          ! Determine the size of the local matrix lb (lnrow_b x lncol_b):
          lnrow_b = max(numroc(m, mbrow, irow, 0, nprocrow),1)
          lncol_b = max(numroc(n, mbcol, icol, 0, nproccol),1)
          !write(*,'(a,i0,a,i0,a,i0)') 'iproc ',iproc,' will have a local b of size ',lnrow_b,' x ',lncol_b
    
          ! Determine the size of the local matrix lc (lnrow_c x lncol_c):
          lnrow_c = max(numroc(m, mbrow, irow, 0, nprocrow),1)
          lncol_c = max(numroc(n, mbcol, icol, 0, nproccol),1)
          !write(*,'(a,i0,a,i0,a,i0)') 'iproc ',iproc,' will have a local c of size ',lnrow_c,' x ',lncol_c
      
          ! Initialize descriptor arrays.
          call descinit(desc_la, m, m, mbrow, mbcol, 0, 0, context, lnrow_a, info)
          call descinit(desc_lb, m, n, mbrow, mbcol, 0, 0, context, lnrow_b, info)
          call descinit(desc_lc, m, n, mbrow, mbcol, 0, 0, context, lnrow_c, info)
      
          ! Allocate the local arrays
          la = f_malloc((/ lnrow_a, lncol_a /),id='la')
          lb = f_malloc((/ lnrow_b, lncol_b /),id='lb')
          lc = f_malloc((/ lnrow_c, lncol_c /),id='lc')
      
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
          call f_free(la)
          call f_free(lb)
          call f_free(lc)
    
          call blacs_gridexit(context)
      
      end if processIf
      
      
      ! Gather the result on all processes.
      if (nproc > 1) then
         call mpiallred(c(1,1), m*n, mpi_sum, comm)
      end if
    
      !call blacs_exit(0)
    
    end subroutine dsymm_parallel
