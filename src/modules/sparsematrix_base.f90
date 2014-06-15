!> @file
!!  File defining the structures to deal with the sparse matrices
!! @author
!!    Copyright (C) 2014-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module defining the basic operations with sparse matrices (creation and destruction)
module sparsematrix_base
  use module_base
  implicit none

  private

  !> Contains the matrices
  type,public :: matrices
      real(kind=8),dimension(:),pointer :: matrix_compr,matrix_comprp
      real(kind=8),dimension(:,:),pointer :: matrix,matrixp
  end type matrices

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


  type, public :: sparse_matrix_info_ptr
     character(len=32) :: id !< name of the sparse matrix array
     integer :: iaction !< action for allocation
     type(sparse_matrix), pointer :: smat !<information of the sparse matrix
  end type sparse_matrix_info_ptr

  type, public :: sparse_matrix_info
     character(len=32) :: id !< name of the sparse matrix array
     integer :: iaction !< action for allocation
     type(sparse_matrix), pointer :: smat !<information of the sparse matrix
  end type sparse_matrix_info

  type, public :: sparse_matrix_info0_ptr
     character(len=32) :: id !< name of the sparse matrix array
     integer :: iaction !< action for allocation
     type(sparse_matrix), pointer :: smat !<information of the sparse matrix
  end type sparse_matrix_info0_ptr

  type, public :: sparse_matrix_info0
     character(len=32) :: id !< name of the sparse matrix array
     integer :: iaction !< action for allocation
     type(sparse_matrix), pointer :: smat !<information of the sparse matrix
  end type sparse_matrix_info0

!!  interface sparsematrix_allocate
!!    module procedure sparsematrix_allocate_1D, sparsematrix_allocate_2D
!!  end interface sparsematrix_allocate

  interface assignment(=)
     module procedure allocate_smat_d1_ptr,allocate_smat_d2_ptr, &
                      allocate_smat_d1,allocate_smat_d2, &
                      allocate0_smat_d1_ptr,allocate0_smat_d2_ptr, &
                      allocate0_smat_d1,allocate0_smat_d2
  end interface


  !> Public routines
  public :: deallocate_sparse_matrix
  public :: sparse_matrix_null
  public :: allocate_sparse_matrix_keys
  public :: allocate_sparse_matrix_basic
  public :: allocate_sparse_matrix_matrices
  public :: allocate_sparse_matrix_matrix_multiplication
  public :: sparsematrix_malloc_ptr
  public :: sparsematrix_malloc
  public :: sparsematrix_malloc0_ptr
  public :: sparsematrix_malloc0
  public :: assignment(=)
  public :: allocate_matrices
  public :: deallocate_matrices
  public :: matrices_null

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
      sparsemat%keyv=f_malloc_ptr(sparsemat%nseg,id='sparsemat%keyv')
      sparsemat%keyg=f_malloc_ptr((/2,sparsemat%nseg/),id='sparsemat%keyg')
      sparsemat%orb_from_index=f_malloc_ptr((/2,sparsemat%nvctr/),id='sparsemat%orb_from_index')
    end subroutine allocate_sparse_matrix_keys


    subroutine allocate_sparse_matrix_matrices(sparsemat,allocate_full)
      implicit none
      type(sparse_matrix),intent(inout) :: sparsemat
      logical,intent(in) :: allocate_full
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


    pure function matrices_null() result(mat)
      implicit none
      type(matrices) :: mat
      call nullify_matrices(mat)
    end function matrices_null


    pure subroutine nullify_matrices(mat)
      implicit none
      type(matrices),intent(out):: mat
      nullify(mat%matrix)
      nullify(mat%matrix_compr)
      nullify(mat%matrixp)
      nullify(mat%matrix_comprp)
    end subroutine nullify_matrices

    subroutine allocate_matrices(sparsemat, allocate_full, matname, mat)
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in) :: sparsemat
      logical,intent(in) :: allocate_full
      character(len=*),intent(in) :: matname
      type(matrices),intent(out) :: mat

      mat%matrix_compr = sparsematrix_malloc_ptr(sparsemat, iaction=SPARSE_FULL, id=trim(matname)//'%matrix_compr')
      mat%matrix_comprp = sparsematrix_malloc_ptr(sparsemat, iaction=SPARSE_PARALLEL, id=trim(matname)//'%matrix_comprp')
      if (allocate_full) mat%matrix = sparsematrix_malloc_ptr(sparsemat, iaction=DENSE_FULL, id=trim(matname)//'%matrix')
      mat%matrixp = sparsematrix_malloc_ptr(sparsemat, iaction=DENSE_PARALLEL, id=trim(matname)//'%matrixp')
    end subroutine allocate_matrices


    subroutine deallocate_matrices(mat)
      implicit none

      ! Calling arguments
      type(matrices),intent(inout) :: mat

      call f_free_ptr(mat%matrix_compr)
      call f_free_ptr(mat%matrix_comprp)
      call f_free_ptr(mat%matrix)
      call f_free_ptr(mat%matrixp)
    end subroutine deallocate_matrices


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

    subroutine allocate_smat_d1_ptr(smat_ptr,smat_info_ptr)
      implicit none
      double precision,dimension(:),pointer,intent(inout) :: smat_ptr
      type(sparse_matrix_info_ptr), intent(in) :: smat_info_ptr

      select case (smat_info_ptr%iaction)
      case (SPARSE_FULL)
          smat_ptr = f_malloc_ptr(smat_info_ptr%smat%nvctr,id=smat_info_ptr%id)
      case (SPARSE_PARALLEL)
          smat_ptr = f_malloc_ptr(smat_info_ptr%smat%nvctrp,id=smat_info_ptr%id)
      case (SPARSEMM_SEQ)
          smat_ptr = f_malloc_ptr(smat_info_ptr%smat%smmm%nseq,id=smat_info_ptr%id)
      case default
          call f_err_throw('The action specified for the 1d matrix allocation is invalid',&
               err_name='BIGDFT_RUNTIME_ERROR')
      end select
    end subroutine allocate_smat_d1_ptr


    subroutine allocate_smat_d2_ptr(smat_ptr,smat_info_ptr)
      implicit none
      double precision,dimension(:,:),pointer,intent(inout) :: smat_ptr
      type(sparse_matrix_info_ptr), intent(in) :: smat_info_ptr

      select case (smat_info_ptr%iaction)
      case (DENSE_FULL)
          smat_ptr = f_malloc_ptr((/smat_info_ptr%smat%nfvctr,smat_info_ptr%smat%nfvctr/),id=smat_info_ptr%id)
      case (DENSE_PARALLEL)
          smat_ptr = f_malloc_ptr((/smat_info_ptr%smat%nfvctr,smat_info_ptr%smat%nvctrp/),id=smat_info_ptr%id)
      case default
         call f_err_throw('The action specified for the 2d matrix allocation is invalid',&
              err_name='BIGDFT_RUNTIME_ERROR')
      end select
    end subroutine allocate_smat_d2_ptr


    subroutine allocate_smat_d1(smat,smat_info)
      implicit none
      double precision,dimension(:),allocatable,intent(inout) :: smat
      type(sparse_matrix_info), intent(in) :: smat_info

      select case (smat_info%iaction)
      case (SPARSE_FULL)
          smat = f_malloc(smat_info%smat%nvctr,id=smat_info%id)
      case (SPARSE_PARALLEL)
          smat = f_malloc(smat_info%smat%nvctrp,id=smat_info%id)
      case (SPARSEMM_SEQ)
          smat = f_malloc(smat_info%smat%smmm%nseq,id=smat_info%id)
      case default
          call f_err_throw('The action specified for the 1d matrix allocation is invalid',&
               err_name='BIGDFT_RUNTIME_ERROR')
      end select
    end subroutine allocate_smat_d1


    subroutine allocate_smat_d2(smat,smat_info)
      implicit none
      double precision,dimension(:,:),allocatable,intent(inout) :: smat
      type(sparse_matrix_info), intent(in) :: smat_info

      select case (smat_info%iaction)
      case (DENSE_FULL)
          smat = f_malloc((/smat_info%smat%nfvctr,smat_info%smat%nfvctr/),id=smat_info%id)
      case (DENSE_PARALLEL)
          smat = f_malloc((/smat_info%smat%nfvctr,smat_info%smat%nvctrp/),id=smat_info%id)
      case default
         call f_err_throw('The action specified for the 2d matrix allocation is invalid',&
              err_name='BIGDFT_RUNTIME_ERROR')
      end select
    end subroutine allocate_smat_d2


    subroutine allocate0_smat_d1_ptr(smat_ptr,smat_info0_ptr)
      implicit none
      double precision,dimension(:),pointer,intent(inout) :: smat_ptr
      type(sparse_matrix_info0_ptr), intent(in) :: smat_info0_ptr

      select case (smat_info0_ptr%iaction)
      case (SPARSE_FULL)
          smat_ptr = f_malloc0_ptr(smat_info0_ptr%smat%nvctr,id=smat_info0_ptr%id)
          !call to_zero(smat_info0_ptr%smat%nvctr,smat_ptr(1))
      case (SPARSE_PARALLEL)
          smat_ptr = f_malloc0_ptr(smat_info0_ptr%smat%nvctrp,id=smat_info0_ptr%id)
          !call to_zero(smat_info0_ptr%smat%nvctrp,smat_ptr(1))
      case (SPARSEMM_SEQ)
          smat_ptr = f_malloc0_ptr(smat_info0_ptr%smat%smmm%nseq,id=smat_info0_ptr%id)
          !call to_zero(smat_info0_ptr%smat%smmm%nseq,smat_ptr(1))
      case default
          call f_err_throw('The action specified for the 1d matrix allocation is invalid',&
               err_name='BIGDFT_RUNTIME_ERROR')
      end select
    end subroutine allocate0_smat_d1_ptr


    subroutine allocate0_smat_d2_ptr(smat_ptr,smat_info0_ptr)
      implicit none
      double precision,dimension(:,:),pointer,intent(inout) :: smat_ptr
      type(sparse_matrix_info0_ptr), intent(in) :: smat_info0_ptr

      select case (smat_info0_ptr%iaction)
      case (DENSE_FULL)
          smat_ptr = f_malloc0_ptr((/smat_info0_ptr%smat%nfvctr,smat_info0_ptr%smat%nfvctr/),id=smat_info0_ptr%id)
          !call to_zero(smat_info0_ptr%smat%nfvctr*smat_info0_ptr%smat%nfvctr,smat_ptr(1,1))
      case (DENSE_PARALLEL)
          smat_ptr = f_malloc0_ptr((/smat_info0_ptr%smat%nfvctr,smat_info0_ptr%smat%nvctrp/),id=smat_info0_ptr%id)
          !call to_zero(smat_info0_ptr%smat%nfvctr*smat_info0_ptr%smat%nvctrp,smat_ptr(1,1))
      case default
         call f_err_throw('The action specified for the 2d matrix allocation is invalid',&
              err_name='BIGDFT_RUNTIME_ERROR')
      end select
    end subroutine allocate0_smat_d2_ptr


    subroutine allocate0_smat_d1(smat,smat_info0)
      implicit none
      double precision,dimension(:),allocatable,intent(inout) :: smat
      type(sparse_matrix_info0), intent(in) :: smat_info0

      select case (smat_info0%iaction)
      case (SPARSE_FULL)
          smat = f_malloc0(smat_info0%smat%nvctr,id=smat_info0%id)
          !call to_zero(smat_info0%smat%nvctr,smat(1))
      case (SPARSE_PARALLEL)
          smat = f_malloc0(smat_info0%smat%nvctrp,id=smat_info0%id)
          !call to_zero(smat_info0%smat%nvctrp,smat(1))
      case (SPARSEMM_SEQ)
          smat = f_malloc0(smat_info0%smat%smmm%nseq,id=smat_info0%id)
          !call to_zero(smat_info0%smat%smmm%nseq,smat(1))
      case default
          call f_err_throw('The action specified for the 1d matrix allocation is invalid',&
               err_name='BIGDFT_RUNTIME_ERROR')
      end select
    end subroutine allocate0_smat_d1


    subroutine allocate0_smat_d2(smat,smat_info0)
      implicit none
      double precision,dimension(:,:),allocatable,intent(inout) :: smat
      type(sparse_matrix_info0), intent(in) :: smat_info0

      select case (smat_info0%iaction)
      case (DENSE_FULL)
          smat = f_malloc0((/smat_info0%smat%nfvctr,smat_info0%smat%nfvctr/),id=smat_info0%id)
          !call to_zero(smat_info0%smat%nfvctr*smat_info0%smat%nfvctr,smat(1,1))
      case (DENSE_PARALLEL)
          smat = f_malloc0((/smat_info0%smat%nfvctr,smat_info0%smat%nvctrp/),id=smat_info0%id)
          !call to_zero(smat_info0%smat%nfvctr*smat_info0%smat%nvctrp,smat(1,1))
      case default
          call f_err_throw('The action specified for the 2d matrix allocation is invalid',&
               err_name='BIGDFT_RUNTIME_ERROR')
      end select
    end subroutine allocate0_smat_d2



    function sparsematrix_malloc_ptr(smat, iaction, id) result(smat_info_ptr)
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in), target :: smat
      integer :: iaction
      character(len=*),intent(in) :: id
      type(sparse_matrix_info_ptr) :: smat_info_ptr

      smat_info_ptr%id(1:len(smat_info_ptr%id))=id
      smat_info_ptr%iaction=iaction
      smat_info_ptr%smat=>smat
    end function sparsematrix_malloc_ptr


    function sparsematrix_malloc(smat, iaction, id) result(smat_info)
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in), target :: smat
      integer :: iaction
      character(len=*),intent(in) :: id
      type(sparse_matrix_info) :: smat_info

      smat_info%id(1:len(smat_info%id))=id
      smat_info%iaction=iaction
      smat_info%smat=>smat
    end function sparsematrix_malloc


    function sparsematrix_malloc0_ptr(smat, iaction, id) result(smat_info0_ptr)
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in), target :: smat
      integer :: iaction
      character(len=*),intent(in) :: id
      type(sparse_matrix_info0_ptr) :: smat_info0_ptr

      smat_info0_ptr%id(1:len(smat_info0_ptr%id))=id
      smat_info0_ptr%iaction=iaction
      smat_info0_ptr%smat=>smat
    end function sparsematrix_malloc0_ptr


    function sparsematrix_malloc0(smat, iaction, id) result(smat_info0)
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in), target :: smat
      integer :: iaction
      character(len=*),intent(in) :: id
      type(sparse_matrix_info0) :: smat_info0

      smat_info0%id(1:len(smat_info0%id))=id
      smat_info0%iaction=iaction
      smat_info0%smat=>smat
    end function sparsematrix_malloc0
end module sparsematrix_base
