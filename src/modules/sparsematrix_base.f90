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


    subroutine allocate_sparsematrix_basic(store_index, norb, sparsemat)
      implicit none
      logical,intent(in) :: store_index
      integer,intent(in) :: norb
      type(sparseMatrix),intent(inout) :: sparsemat
      integer :: istat
      character(len=*),parameter :: subname='allocate_sparsematrix_basic'
      allocate(sparsemat%nsegline(norb), stat=istat)
      call memocc(istat, sparsemat%nsegline, 'sparsemat%nsegline', subname)
      allocate(sparsemat%istsegline(norb), stat=istat)
      call memocc(istat, sparsemat%istsegline, 'sparsemat%istsegline', subname)
      if (store_index) then
          allocate(sparsemat%matrixindex_in_compressed_arr(norb,norb), stat=istat)
          call memocc(istat, sparsemat%matrixindex_in_compressed_arr, 'sparsemat%matrixindex_in_compressed_arr', subname)
      end if
      sparsemat%nvctr_par=f_malloc_ptr((/0.to.nproc-1/),id='sparsemat%nvctr_par')
      sparsemat%isvctr_par=f_malloc_ptr((/0.to.nproc-1/),id='sparsemat%isvctr_par')
    end subroutine allocate_sparsematrix_basic

    subroutine allocate_sparsematrix_keys(sparsemat)
      implicit none
      type(sparseMatrix),intent(inout) :: sparsemat
      integer :: istat
      character(len=*),parameter :: subname='allocate_sparsematrix_keys'
      allocate(sparsemat%keyv(sparsemat%nseg), stat=istat)
      call memocc(istat, sparsemat%keyv, 'sparsemat%keyv', subname)
      allocate(sparsemat%keyg(2,sparsemat%nseg), stat=istat)
      call memocc(istat, sparsemat%keyg, 'sparsemat%keyg', subname)
    end subroutine allocate_sparsematrix_keys


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
