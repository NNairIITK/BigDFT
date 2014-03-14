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
  public :: deallocate_sparseMatrix
  public :: sparsematrix_null
  public :: allocate_sparsematrix_keys
  public :: allocate_sparsematrix_basic


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


    subroutine allocate_sparsematrix_basic(store_index, norb, nproc, sparsemat)
      implicit none
      logical,intent(in) :: store_index
      integer,intent(in) :: norb, nproc
      type(sparseMatrix),intent(inout) :: sparsemat
      integer :: istat
      character(len=*),parameter :: subname='allocate_sparsematrix_basic'
      sparsemat%nsegline=f_malloc_ptr(norb,id='sparsemat%nsegline')
      sparsemat%istsegline=f_malloc_ptr(norb,id='sparsemat%istsegline')
      if (store_index) then
          sparsemat%matrixindex_in_compressed_arr=f_malloc_ptr((/norb,norb/),id='sparsemat%matrixindex_in_compressed_arr')
      end if
      sparsemat%nvctr_par=f_malloc_ptr((/0.to.nproc-1/),id='sparsemat%nvctr_par')
      sparsemat%isvctr_par=f_malloc_ptr((/0.to.nproc-1/),id='sparsemat%isvctr_par')
    end subroutine allocate_sparsematrix_basic

    subroutine allocate_sparsematrix_keys(sparsemat)
      implicit none
      type(sparseMatrix),intent(inout) :: sparsemat
      integer :: istat
      character(len=*),parameter :: subname='allocate_sparsematrix_keys'
      sparsemat%keyv=f_malloc_ptr(sparsemat%nseg,id='sparsemat%keyv')
      sparsemat%keyg=f_malloc_ptr((/2,sparsemat%nseg/),id='sparsemat%keyg')
      sparsemat%orb_from_index=f_malloc_ptr((/2,sparsemat%nvctr/),id='sparsemat%orb_from_index')
    end subroutine allocate_sparsematrix_keys


    subroutine deallocate_sparseMatrix(sparsemat, subname)
      use module_base 
      implicit none
      ! Calling arguments
      type(sparseMatrix),intent(inout):: sparsemat
      character(len=*),intent(in):: subname
    
      if (associated(sparseMat%keyg)) call f_free_ptr(sparseMat%keyg)
      if (associated(sparseMat%keyv)) call f_free_ptr(sparseMat%keyv)
      if (associated(sparseMat%nsegline)) call f_free_ptr(sparseMat%nsegline)
      if (associated(sparseMat%istsegline)) call f_free_ptr(sparseMat%istsegline)
      if (associated(sparseMat%matrixindex_in_compressed_fortransposed)) &
          call f_free_ptr(sparseMat%matrixindex_in_compressed_fortransposed)
      if (associated(sparseMat%matrixindex_in_compressed_arr)) &
          call f_free_ptr(sparseMat%matrixindex_in_compressed_arr)
      if (associated(sparseMat%matrix_compr)) call f_free_ptr(sparseMat%matrix_compr)
      if (associated(sparseMat%matrix)) call f_free_ptr(sparseMat%matrix)
      if (associated(sparseMat%isvctr_par)) call f_free_ptr(sparseMat%isvctr_par)
      if (associated(sparseMat%nvctr_par)) call f_free_ptr(sparseMat%nvctr_par)
      if (associated(sparseMat%isfvctr_par)) call f_free_ptr(sparseMat%isfvctr_par)
      if (associated(sparseMat%nfvctr_par)) call f_free_ptr(sparseMat%nfvctr_par)
      if (associated(sparseMat%matrixp)) call f_free_ptr(sparseMat%matrixp)
      if (associated(sparseMat%matrix_comprp)) call f_free_ptr(sparseMat%matrix_comprp)
      if (associated(sparseMat%orb_from_index)) call f_free_ptr(sparseMat%orb_from_index)
    
    end subroutine deallocate_sparseMatrix


end module sparsematrix_base
