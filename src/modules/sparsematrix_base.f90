module sparsematrix_base
  use module_base
  implicit none

  private

  type,public :: sparseMatrix
      integer :: nvctr, nseg, nvctrp, isvctr, parallel_compression, nfvctr, nfvctrp, isfvctr
      integer,dimension(:),pointer :: keyv, nsegline, istsegline, isvctr_par, nvctr_par, isfvctr_par, nfvctr_par
      integer,dimension(:,:),pointer :: keyg
      real(kind=8),dimension(:),pointer :: matrix_compr,matrix_comprp
      real(kind=8),dimension(:,:),pointer :: matrix,matrixp
      integer,dimension(:,:),pointer :: matrixindex_in_compressed_arr, orb_from_index
      integer,dimension(:,:),pointer :: matrixindex_in_compressed_fortransposed
      logical :: store_index, can_use_dense
  end type sparseMatrix


  !> Public routines
  public :: sparsematrix_null


  contains

    !> Creators and destructors

    pure function sparsematrix_null() result(sparsemat)
      implicit none
      type(sparseMatrix) :: sparsemat
      call nullify_sparsematrix(sparsemat)
    end function sparsematrix_null


    pure subroutine nullify_sparsematrix(sparsemat)
      implicit none
      type(sparseMatrix),intent(out):: sparsemat
      nullify(sparsemat%keyv)
      nullify(sparsemat%nsegline)
      nullify(sparsemat%keyg)
      nullify(sparsemat%istsegline)
      nullify(sparsemat%matrix)
      nullify(sparsemat%matrix_compr)
      nullify(sparsemat%matrixp)
      nullify(sparsemat%matrix_comprp)
      nullify(sparsemat%matrixindex_in_compressed_arr)
      nullify(sparsemat%orb_from_index)
      nullify(sparsemat%matrixindex_in_compressed_fortransposed)
      nullify(sparsemat%nvctr_par)
      nullify(sparsemat%isvctr_par)
      nullify(sparsemat%nfvctr_par)
      nullify(sparsemat%isfvctr_par)
    
    end subroutine nullify_sparsematrix


    !!subroutine deallocate_sparseMatrix(sparsemat, subname)
    !!  use module_base 
    !!  implicit none
    !!  ! Calling arguments
    !!  type(sparseMatrix),intent(inout):: sparsemat
    !!  character(len=*),intent(in):: subname
    !!
    !!  call checkAndDeallocatePointer(sparseMat%keyg, 'sparseMat%keyg', subname)
    !!  call checkAndDeallocatePointer(sparseMat%keyv, 'sparseMat%keyv', subname)
    !!  call checkAndDeallocatePointer(sparseMat%nsegline, 'sparseMat%nsegline', subname)
    !!  call checkAndDeallocatePointer(sparseMat%istsegline, 'sparseMat%istsegline', subname)
    !!  call checkAndDeallocatePointer(sparseMat%matrix_compr, 'sparseMat%matrix_compr', subname)
    !!  call checkAndDeallocatePointer(sparseMat%matrix, 'sparseMat%matrix', subname)
    !!  !call checkAndDeallocatePointer(sparseMat%matrixindex_in_compressed, 'sparseMat%matrixindex_in_compressed', subname)
    !!  call checkAndDeallocatePointer(sparseMat%matrixindex_in_compressed_arr, 'sparseMat%matrixindex_in_compressed_arr', subname)
    !!  call checkAndDeallocatePointer(sparseMat%orb_from_index, 'sparseMat%orb_from_index', subname)
    !!  call checkAndDeallocatePointer(sparseMat%matrixindex_in_compressed_fortransposed, &
    !!       'sparseMat%matrixindex_in_compressed_fortransposed', subname)
    !!  call f_free_ptr(sparseMat%isvctr_par)
    !!  call f_free_ptr(sparseMat%nvctr_par)
    !!  call f_free_ptr(sparseMat%isfvctr_par)
    !!  call f_free_ptr(sparseMat%nfvctr_par)
    !!  call f_free_ptr(sparseMat%matrixp)
    !!  call f_free_ptr(sparseMat%matrix_comprp)
    !!
    !!end subroutine deallocate_sparseMatrix


end module sparsematrix_base
