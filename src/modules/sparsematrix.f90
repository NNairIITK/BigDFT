module sparsematrix
  use sparsematrix_base
  implicit none

  private

  public :: compress_matrix_for_allreduce
  public :: uncompressMatrix

  contains

    subroutine compress_matrix_for_allreduce(iproc,sparsemat)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix
      implicit none
      
      ! Calling arguments
      integer, intent(in) :: iproc
      type(sparse_matrix),intent(inout) :: sparsemat
    
      ! Local variables
      integer :: jj, irow, jcol, jjj, ierr
    
      call timing(iproc,'compress_uncom','ON')
    
      if (sparsemat%parallel_compression==0.or.bigdft_mpi%nproc==1) then
         !$omp parallel do default(private) shared(sparsemat)
         do jj=1,sparsemat%nvctr
            irow = sparsemat%orb_from_index(1,jj)
            jcol = sparsemat%orb_from_index(2,jj)
            sparsemat%matrix_compr(jj)=sparsemat%matrix(irow,jcol)
         end do
         !$omp end parallel do
      else if (sparsemat%parallel_compression==1) then
         call to_zero(sparsemat%nvctr, sparsemat%matrix_compr(1))
         !$omp parallel do default(private) shared(sparsemat)
         do jj=1,sparsemat%nvctrp
            jjj=jj+sparsemat%isvctr
            irow = sparsemat%orb_from_index(1,jjj)
            jcol = sparsemat%orb_from_index(2,jjj)
            sparsemat%matrix_compr(jjj)=sparsemat%matrix(irow,jcol)
         end do
         !$omp end parallel do
         call mpiallred(sparsemat%matrix_compr(1), sparsemat%nvctr, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      else
         sparsemat%matrix_comprp=f_malloc_ptr((sparsemat%nvctrp),id='sparsemat%matrix_comprp')
         !$omp parallel do default(private) shared(sparsemat)
         do jj=1,sparsemat%nvctrp
            jjj=jj+sparsemat%isvctr
            irow = sparsemat%orb_from_index(1,jjj)
            jcol = sparsemat%orb_from_index(2,jjj)
            sparsemat%matrix_comprp(jj)=sparsemat%matrix(irow,jcol)
         end do
         !$omp end parallel do
         call mpi_allgatherv(sparsemat%matrix_comprp, sparsemat%nvctrp, mpi_double_precision, sparsemat%matrix_compr, &
              sparsemat%nvctr_par(:), sparsemat%isvctr_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
         call f_free_ptr(sparsemat%matrix_comprp)
      end if
    
      call timing(iproc,'compress_uncom','OF')
    
    end subroutine compress_matrix_for_allreduce



    subroutine uncompressMatrix(iproc,sparsemat)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix
      implicit none
      
      ! Calling arguments
      integer, intent(in) :: iproc
      type(sparse_matrix), intent(inout) :: sparsemat
      
      ! Local variables
      integer :: ii, irow, jcol, iii, ierr
    
      call timing(iproc,'compress_uncom','ON')
    
      if (sparsemat%parallel_compression==0.or.bigdft_mpi%nproc==1) then
         call to_zero(sparsemat%nfvctr**2, sparsemat%matrix(1,1))
         !$omp parallel do default(private) shared(sparsemat)
         do ii=1,sparsemat%nvctr
            irow = sparsemat%orb_from_index(1,ii)
            jcol = sparsemat%orb_from_index(2,ii)
            sparsemat%matrix(irow,jcol)=sparsemat%matrix_compr(ii)
         end do
         !$omp end parallel do
      else if (sparsemat%parallel_compression==1) then
         call to_zero(sparsemat%nfvctr**2, sparsemat%matrix(1,1))
         !$omp parallel do default(private) shared(sparsemat)
         do ii=1,sparsemat%nvctrp
            iii=ii+sparsemat%isvctr
            irow = sparsemat%orb_from_index(1,iii)
            jcol = sparsemat%orb_from_index(2,iii)
            sparsemat%matrix(irow,jcol)=sparsemat%matrix_compr(iii)
         end do
         !$omp end parallel do
         call mpiallred(sparsemat%matrix(1,1), sparsemat%nfvctr**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      else
         sparsemat%matrixp=f_malloc_ptr((/sparsemat%nfvctr,sparsemat%nfvctrp/),id='sparsemat%matrixp')
         call to_zero(sparsemat%nfvctr*sparsemat%nfvctrp, sparsemat%matrixp(1,1))
         !$omp parallel do default(private) shared(sparsemat)
         do ii=1,sparsemat%nvctrp
            iii=ii+sparsemat%isvctr
            irow = sparsemat%orb_from_index(1,iii)
            jcol = sparsemat%orb_from_index(2,iii) - sparsemat%isfvctr
            sparsemat%matrixp(irow,jcol)=sparsemat%matrix_compr(iii)
         end do
         !$omp end parallel do
         call mpi_allgatherv(sparsemat%matrixp, sparsemat%nfvctr*sparsemat%nfvctrp, mpi_double_precision, sparsemat%matrix, &
              sparsemat%nfvctr*sparsemat%nfvctr_par(:), sparsemat%nfvctr*sparsemat%isfvctr_par, mpi_double_precision, &
              bigdft_mpi%mpi_comm, ierr)
         call f_free_ptr(sparsemat%matrixp)
      end if
      sparsemat%can_use_dense=.true.  
    
      call timing(iproc,'compress_uncom','OF')
    
    end subroutine uncompressMatrix


end module sparsematrix
