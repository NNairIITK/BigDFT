module sparsematrix
  use sparsematrix_base
  implicit none

  private

  public :: compress_matrix_for_allreduce
  public :: uncompressMatrix
  public :: check_matrix_compression

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



    subroutine check_matrix_compression(iproc,sparsemat)
      use module_base
      use module_types
      use module_interfaces
      use yaml_output
      use sparsematrix_base, only: sparse_matrix
      use sparsematrix, only: compress_matrix_for_allreduce, uncompressMatrix
      implicit none
      integer,intent(in) :: iproc
      type(sparse_matrix),intent(inout) :: sparsemat
      !Local variables
      integer :: i_stat, i_all, jorb, irow, icol, iseg, ii
      character(len=*),parameter :: subname='check_matrix_compression'
      real(kind=8) :: maxdiff
      real(kind=8), parameter :: tol=1.e-10
    
    
      !!allocate(sparsemat%matrix(sparsemat%nfvctr,sparsemat%nfvctr),stat=i_stat)
      !!call memocc(i_stat,sparsemat%matrix,'sparsemat%matrix',subname)
    
      !!allocate(sparsemat%matrix_compr(sparsemat%nvctr),stat=i_stat)
      !!call memocc(i_stat,sparsemat%matrix_compr,'sparsemat%matrix_compr',subname)
    
      sparsemat%matrix=f_malloc_ptr((/sparsemat%nfvctr,sparsemat%nfvctr/),id='sparsemat%matrix')
      sparsemat%matrix_compr=f_malloc_ptr(sparsemat%nvctr,id='sparsemat%matrix_compr')
    
      call to_zero(sparsemat%nfvctr**2,sparsemat%matrix(1,1))
      do iseg = 1, sparsemat%nseg
         do jorb = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
            call get_indecies(jorb,irow,icol)
            !print *,'irow,icol',irow, icol,test_value_matrix(sparsemat%nfvctr, irow, icol)
            sparsemat%matrix(irow,icol) = test_value_matrix(sparsemat%nfvctr, irow, icol)
         end do
      end do
      
      call compress_matrix_for_allreduce(iproc,sparsemat)
    
      maxdiff = 0.d0
      do iseg = 1, sparsemat%nseg
         ii=0
         do jorb = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
            call get_indecies(jorb,irow,icol)
            maxdiff = max(abs(sparsemat%matrix_compr(sparsemat%keyv(iseg)+ii)&
                 -test_value_matrix(sparsemat%nfvctr, irow, icol)),maxdiff)
            ii=ii+1
         end do
      end do
    
      if (iproc==0) call yaml_map('Tolerances for this check',tol,fmt='(1pe25.17)')
    
      if(iproc==0) then
        if (maxdiff > tol) then
           call yaml_warning('COMPRESSION ERROR : difference of '//trim(yaml_toa(maxdiff,fmt='(1pe12.5)')))
        else
           call yaml_map('Maxdiff for compress', maxdiff,fmt='(1pe25.17)')
        end if
      end if
    
      call uncompressMatrix(iproc,sparsemat)
    
      maxdiff = 0.d0
      do iseg = 1, sparsemat%nseg
         do jorb = sparsemat%keyg(1,iseg), sparsemat%keyg(2,iseg)
            call get_indecies(jorb,irow,icol)
            maxdiff = max(abs(sparsemat%matrix(irow,icol)-test_value_matrix(sparsemat%nfvctr, irow, icol)),maxdiff) 
         end do
      end do
    
      if(iproc==0) then
        if (maxdiff > tol) then
           call yaml_warning('UNCOMPRESSION ERROR : difference of '//trim(yaml_toa(maxdiff,fmt='(1pe12.5)')))
        else
           call yaml_map('Maxdiff for uncompress', maxdiff,fmt='(1pe25.17)')
        end if
      end if
    
      !!i_all = -product(shape(sparsemat%matrix))*kind(sparsemat%matrix)
      !!deallocate(sparsemat%matrix,stat=i_stat)
      !!call memocc(i_stat,i_all,'sparsemat%matrix',subname)
      !!i_all = -product(shape(sparsemat%matrix_compr))*kind(sparsemat%matrix_compr)
      !!deallocate(sparsemat%matrix_compr,stat=i_stat)
      !!call memocc(i_stat,i_all,'sparsemat%matrix_compr',subname)
      call f_free_ptr(sparsemat%matrix)
      call f_free_ptr(sparsemat%matrix_compr)
    
    contains
       !> define a value for the wavefunction which is dependent of the indices
       function test_value_matrix(norb,iorb,jorb)
          use module_base
          implicit none
          integer, intent(in) :: norb,iorb,jorb
          real(kind=8) :: test_value_matrix
    
          test_value_matrix = norb*(iorb-1)+jorb
          !print *,iorb,jorb,test_value_matrix
       END FUNCTION test_value_matrix
    
       subroutine get_indecies(ind,irow,icol)
         implicit none
         integer, intent(in) :: ind
         integer, intent(out) :: irow, icol
    
         icol = (ind - 1) / sparsemat%nfvctr + 1
         irow = ind - (icol-1)*sparsemat%nfvctr
         !print *,'irow,icol',irow,icol
       END SUBROUTINE get_indecies 
    end subroutine check_matrix_compression


end module sparsematrix
