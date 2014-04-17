module sparsematrix_base
  use module_base
  implicit none

  private

  !> Contains the parameters needed for the sparse matrix matrix multiplication
  type,public :: sparse_matrix_matrix_multiplication
      integer :: nout, nseq, nmaxsegk, nmaxvalk, nseg
      integer,dimension(:),pointer :: ivectorindex, nsegline, istsegline, indices_extract_sequential
      integer,dimension(:,:),pointer :: onedimindices, keyg
  end type sparse_matrix_matrix_multiplication

  type,public :: sparse_matrix
      integer :: nvctr, nseg, nvctrp, isvctr, parallel_compression, nfvctr, nfvctrp, isfvctr
      integer,dimension(:),pointer :: keyv, nsegline, istsegline, isvctr_par, nvctr_par, isfvctr_par, nfvctr_par
      integer,dimension(:,:),pointer :: keyg
      real(kind=8),dimension(:),pointer :: matrix_compr,matrix_comprp
      real(kind=8),dimension(:,:),pointer :: matrix,matrixp
      integer,dimension(:,:),pointer :: matrixindex_in_compressed_arr, orb_from_index
      integer,dimension(:,:),pointer :: matrixindex_in_compressed_fortransposed
      logical :: store_index, can_use_dense
      type(sparse_matrix_matrix_multiplication) :: smmm
  end type sparse_matrix


  !> Public routines
  public :: deallocate_sparse_matrix
  public :: sparse_matrix_null
  public :: allocate_sparse_matrix_keys
  public :: allocate_sparse_matrix_basic
  public :: allocate_sparse_matrix_matrices
  public :: allocate_sparse_matrix_matrix_multiplication

  !> Public constants
  integer,parameter,public :: SPARSE_FULL     = 51
  integer,parameter,public :: SPARSE_PARALLEL = 52
  integer,parameter,public :: DENSE_FULL      = 53
  integer,parameter,public :: DENSE_PARALLEL  = 54
  integer,parameter,public :: SPARSEMM_SEQ    = 55



  contains

    !> Creators and destructors

    pure function sparse_matrix_null() result(sparsemat)
      implicit none
      type(sparse_matrix) :: sparsemat
      call nullify_sparse_matrix(sparsemat)
    end function sparse_matrix_null


    pure subroutine nullify_sparse_matrix(sparsemat)
      implicit none
      type(sparse_matrix),intent(out):: sparsemat
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
      call nullify_sparse_matrix_matrix_multiplication(sparsemat%smmm) 
    end subroutine nullify_sparse_matrix

    pure subroutine nullify_sparse_matrix_matrix_multiplication(smmm)
      implicit none
      type(sparse_matrix_matrix_multiplication),intent(out):: smmm
      nullify(smmm%ivectorindex)
      nullify(smmm%onedimindices)
      nullify(smmm%nsegline)
      nullify(smmm%istsegline)
      nullify(smmm%keyg)
      nullify(smmm%indices_extract_sequential)
    end subroutine nullify_sparse_matrix_matrix_multiplication


    subroutine allocate_sparse_matrix_basic(store_index, norb, nproc, sparsemat)
      implicit none
      logical,intent(in) :: store_index
      integer,intent(in) :: norb, nproc
      type(sparse_matrix),intent(inout) :: sparsemat
      integer :: istat
      sparsemat%nsegline=f_malloc_ptr(norb,id='sparsemat%nsegline')
      sparsemat%istsegline=f_malloc_ptr(norb,id='sparsemat%istsegline')
      if (store_index) then
          sparsemat%matrixindex_in_compressed_arr=f_malloc_ptr((/norb,norb/),id='sparsemat%matrixindex_in_compressed_arr')
      end if
      sparsemat%nvctr_par=f_malloc_ptr((/0.to.nproc-1/),id='sparsemat%nvctr_par')
      sparsemat%isvctr_par=f_malloc_ptr((/0.to.nproc-1/),id='sparsemat%isvctr_par')
    end subroutine allocate_sparse_matrix_basic


    subroutine allocate_sparse_matrix_keys(sparsemat)
      implicit none
      type(sparse_matrix),intent(inout) :: sparsemat
      integer :: istat
      sparsemat%keyv=f_malloc_ptr(sparsemat%nseg,id='sparsemat%keyv')
      sparsemat%keyg=f_malloc_ptr((/2,sparsemat%nseg/),id='sparsemat%keyg')
      sparsemat%orb_from_index=f_malloc_ptr((/2,sparsemat%nvctr/),id='sparsemat%orb_from_index')
    end subroutine allocate_sparse_matrix_keys


    subroutine allocate_sparse_matrix_matrices(sparsemat,allocate_full)
      implicit none
      type(sparse_matrix),intent(inout) :: sparsemat
      logical,intent(in) :: allocate_full
      integer :: istat
      sparsemat%matrix_compr = f_malloc_ptr(sparsemat%nvctr,id='sparsemat%matrix_compr')
      sparsemat%matrix_comprp = f_malloc_ptr(sparsemat%nvctrp,id='sparsemat%matrix_comprp')
      if (allocate_full) sparsemat%matrix = f_malloc_ptr((/sparsemat%nfvctr,sparsemat%nfvctr/),id='sparsemat%matrix')
      sparsemat%matrixp = f_malloc_ptr((/sparsemat%nfvctr,sparsemat%nfvctrp/),id='sparsemat%matrixp')
    end subroutine allocate_sparse_matrix_matrices


    subroutine allocate_sparse_matrix_matrix_multiplication(norb, nseg, nsegline, istsegline, keyg, smmm)
      implicit none
      integer,intent(in) :: norb, nseg
      integer,dimension(norb),intent(in) :: nsegline, istsegline
      integer,dimension(2,nseg),intent(in) :: keyg
      type(sparse_matrix_matrix_multiplication),intent(inout):: smmm
      smmm%ivectorindex=f_malloc_ptr(smmm%nseq,id='smmm%ivectorindex')
      smmm%onedimindices=f_malloc_ptr((/4,smmm%nout/),id='smmm%onedimindices')
      smmm%nsegline=f_malloc_ptr(norb,id='smmm%nsegline')
      smmm%istsegline=f_malloc_ptr(norb,id='smmm%istsegline')
      smmm%keyg=f_malloc_ptr((/2,nseg/),id='smmm%istsegline')
      smmm%indices_extract_sequential=f_malloc_ptr((/smmm%nseq/),id='smmm%indices_extract_sequential')
    end subroutine allocate_sparse_matrix_matrix_multiplication


    subroutine deallocate_sparse_matrix(sparsemat, subname)
      use module_base 
      implicit none
      ! Calling arguments
      type(sparse_matrix),intent(inout):: sparsemat
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
      call deallocate_sparse_matrix_matrix_multiplication(sparsemat%smmm)
    end subroutine deallocate_sparse_matrix

    subroutine deallocate_sparse_matrix_matrix_multiplication(smmm)
      implicit none
      type(sparse_matrix_matrix_multiplication),intent(out):: smmm
      call f_free_ptr(smmm%ivectorindex)
      call f_free_ptr(smmm%onedimindices)
      call f_free_ptr(smmm%nsegline)
      call f_free_ptr(smmm%istsegline)
      call f_free_ptr(smmm%keyg)
      call f_free_ptr(smmm%indices_extract_sequential)
    end subroutine deallocate_sparse_matrix_matrix_multiplication


end module sparsematrix_base
