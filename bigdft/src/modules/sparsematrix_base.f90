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
      real(kind=8),dimension(:,:,:),pointer :: matrix,matrixp
      real(kind=8) :: power !< power of the matrix; eg.: 0 => original matrix, -1 => inverse, 0.5 => ^1/2, etc.
  end type matrices

  !> Contains the parameters needed for the sparse matrix matrix multiplication
  type,public :: sparse_matrix_matrix_multiplication
      integer,dimension(:),pointer :: keyv
      integer,dimension(:,:,:),pointer :: keyg
      integer :: nout, nseq, nseg
      integer :: nfvctrp !< modified number of matrix columns per MPI task for an optimized load balancing during matmul
      integer :: isfvctr !< modified starting column of the matrix for an optimized load balancing during matmul
      integer :: nvctrp_mm !< modified number of compressed matrix elements per MPI task
      integer :: isvctr_mm !< modified starting entry of the compressed matrix elements
      integer :: isseg !< segment containing the first entry (i.e. isvctr+1)
      integer :: ieseg !< segment containing the last entry (i.e. isvctr+nvctrp)
      integer,dimension(:),pointer :: isvctr_mm_par, nvctr_mm_par !<array that contains the values of nvctrp_mm and isvctr_mm of all MPI tasks
      integer,dimension(:),pointer :: ivectorindex, ivectorindex_new, nsegline, istsegline, indices_extract_sequential
      integer,dimension(:,:),pointer :: onedimindices, onedimindices_new, line_and_column_mm, line_and_column
      !!integer,dimension(:,:,:),pointer :: keyg
      integer,dimension(2) :: istartendseg_mm !< starting and ending segments of the matrix subpart which is actually used for the multiplication
                                              !! WARNING: the essential bounds are given by istartend_mm, the segments are used to speed up the code
      integer,dimension(2) :: istartend_mm !< starting and ending indices of the matrix subpart which is actually used for the multiplication
      integer,dimension(2) :: istartend_mm_dj !< starting and ending indices of the submatrices (partitioned in disjoint chunks)
      !!integer :: ncl_smmm !< number of elements for the compress local after a sparse matrix matrix multiplication
      integer :: nccomm_smmm !<number of communications required for the compress distributed after a sparse matrix matrix multiplication
      integer,dimension(:,:),pointer :: luccomm_smmm !<lookup array for the communications required for the compress distributed after a sparse matrix matrix multiplication
      integer :: nvctrp !< number of compressed matrix elements per MPI task
      integer :: isvctr !< starting entry of the compressed matrix elements
      integer,dimension(:),pointer :: isvctr_par, nvctr_par !<array that contains the values of nvctrp and isvctr of all MPI tasks
      integer :: nconsecutive_max !< max number of blocks (i.e. consecutive entries) for the sparse matmul
      integer,dimension(:,:,:),pointer :: consecutive_lookup !< lookup arrays for these blocks
  end type sparse_matrix_matrix_multiplication

  type,public :: sparse_matrix
      integer :: nvctr, nseg, nvctrp, isvctr, parallel_compression, nfvctr, nfvctrp, isfvctr, nspin
      integer :: offset_matrixindex_in_compressed_fortransposed
      integer,dimension(:),pointer :: keyv, nsegline, istsegline, isvctr_par, nvctr_par, isfvctr_par, nfvctr_par
      integer,dimension(:,:,:),pointer :: keyg
      integer,dimension(:,:),pointer :: matrixindex_in_compressed_arr!, orb_from_index
      integer,dimension(:,:),pointer :: matrixindex_in_compressed_fortransposed
      logical :: store_index, can_use_dense
      type(sparse_matrix_matrix_multiplication) :: smmm
      integer :: ntaskgroup !< total number of MPI taskgroups to which this task belongs
      integer :: ntaskgroupp !< number of MPI taskgroups to which this task belongs
      !> (1:2,1,:) gives the start and end of the taskgroups (in terms of total matrix indices), 
      !! (1:2,2,:) gives the start and end of the disjoint submatrix handled by the taskgroup
      !! The last rank is of dimension ntaskgroup
      integer,dimension(:,:,:),pointer :: taskgroup_startend
      integer,dimension(:),pointer :: taskgroupid !< dimension ntaskgroupp, gives the ID of the taskgroups to which a task belongs
      integer,dimension(:,:),pointer :: inwhichtaskgroup !< dimension (2,0:nproc-1), tells in which taskgroup a given task is
      type(mpi_environment),dimension(:),pointer :: mpi_groups
      integer,dimension(2) :: istartendseg_t !< starting and ending indices of the matrix subpart which is actually used i
                                             !! for the transposed operation (overlap calculation / orthocontraint)
                                             !! WARNING: the essential bounds are given by istartend_t, the segments are you used to speed up the code
      integer,dimension(2) :: istartend_t !< starting and ending indices of the matrix subpart which is actually used i
                                          !! for the transposed operation (overlap calculation / orthocontraint
      integer :: nvctrp_tg !< size of the taskgroup matrix (per processor)
      integer :: isvctrp_tg !< offset of the taskgroup matrix, given with respect to the non-parallelized matrix
      integer,dimension(2) :: iseseg_tg !< first and last segment of the taskgroup sparse matrix
      integer,dimension(2) :: istartend_local !< first and last element of the sparse matrix which is actually used by a given MPI task
      integer,dimension(2) :: istartendseg_local !< first and last segment of the sparse matrix which is actually used by a given MPI task
      integer,dimension(:,:),pointer :: tgranks !< global task IDs of the tasks in each taskgroup
      integer,dimension(:),pointer :: nranks !< number of task on each taskgroup
      integer :: nccomm !<number of communications required for the compress distributed in the dense parallel format
      integer,dimension(:,:),pointer :: luccomm !<lookup array for the communications required for the compress distributed in the dense parallel format
      integer,dimension(:),pointer :: on_which_atom !<dimension ntmb, indicates to which atoms a row/column of the matrix belongs
      character(len=1) :: geocode !< boundary conditions F(ree), W(ire), S(urface), P(eriodic)
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


  interface assignment(=)
     module procedure allocate_smat_d1_ptr,allocate_smat_d2_ptr,allocate_smat_d3_ptr, &
                      allocate_smat_d1,allocate_smat_d2,allocate_smat_d3, &
                      allocate0_smat_d1_ptr,allocate0_smat_d2_ptr,allocate0_smat_d3_ptr,&
                      allocate0_smat_d1,allocate0_smat_d2,allocate0_smat_d3
  end interface


  !> Public routines
  public :: deallocate_sparse_matrix
  public :: sparse_matrix_null
  public :: allocate_sparse_matrix_keys
  public :: allocate_sparse_matrix_basic
  public :: allocate_sparse_matrix_matrix_multiplication
  public :: sparsematrix_malloc_ptr
  public :: sparsematrix_malloc
  public :: sparsematrix_malloc0_ptr
  public :: sparsematrix_malloc0
  public :: assignment(=)
  public :: allocate_matrices
  public :: deallocate_matrices
  public :: matrices_null
  public :: copy_sparse_matrix,copy_matrices

  !> Public constants
  integer,parameter,public :: SPARSE_TASKGROUP    = 50
  integer,parameter,public :: SPARSE_FULL         = 51
  integer,parameter,public :: SPARSE_PARALLEL     = 52
  integer,parameter,public :: SPARSE_MATMUL_SMALL = 53
  integer,parameter,public :: SPARSE_MATMUL_LARGE = 54
  integer,parameter,public :: DENSE_FULL          = 55
  integer,parameter,public :: DENSE_PARALLEL      = 56
  integer,parameter,public :: DENSE_MATMUL        = 57
  integer,parameter,public :: SPARSEMM_SEQ        = 58



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
      nullify(sparsemat%matrixindex_in_compressed_arr)
      !nullify(sparsemat%orb_from_index)
      nullify(sparsemat%matrixindex_in_compressed_fortransposed)
      nullify(sparsemat%nvctr_par)
      nullify(sparsemat%isvctr_par)
      nullify(sparsemat%nfvctr_par)
      nullify(sparsemat%isfvctr_par)
      nullify(sparsemat%taskgroup_startend)
      nullify(sparsemat%taskgroupid)
      nullify(sparsemat%mpi_groups)
      nullify(sparsemat%inwhichtaskgroup)
      nullify(sparsemat%tgranks)
      nullify(sparsemat%nranks)
      nullify(sparsemat%luccomm)
      nullify(sparsemat%on_which_atom)
      call nullify_sparse_matrix_matrix_multiplication(sparsemat%smmm) 
    end subroutine nullify_sparse_matrix

    pure subroutine nullify_sparse_matrix_matrix_multiplication(smmm)
      implicit none
      type(sparse_matrix_matrix_multiplication),intent(out):: smmm
      nullify(smmm%ivectorindex)
      nullify(smmm%ivectorindex_new)
      nullify(smmm%onedimindices)
      nullify(smmm%onedimindices_new)
      nullify(smmm%line_and_column_mm)
      nullify(smmm%line_and_column)
      nullify(smmm%nsegline)
      nullify(smmm%istsegline)
      !!nullify(smmm%keyg)
      nullify(smmm%indices_extract_sequential)
      nullify(smmm%nvctr_par)
      nullify(smmm%isvctr_par)
      nullify(smmm%nvctr_mm_par)
      nullify(smmm%isvctr_mm_par)
      nullify(smmm%luccomm_smmm)
      nullify(smmm%keyv)
      nullify(smmm%keyg)
      nullify(smmm%consecutive_lookup)
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
      sparsemat%on_which_atom=f_malloc_ptr(norb,id='sparsemat%on_which_atom')
    end subroutine allocate_sparse_matrix_basic


    subroutine allocate_sparse_matrix_keys(store_index, sparsemat)
      implicit none
      logical,intent(in) :: store_index
      type(sparse_matrix),intent(inout) :: sparsemat
      sparsemat%keyv=f_malloc_ptr(sparsemat%nseg,id='sparsemat%keyv')
      sparsemat%keyg=f_malloc_ptr((/2,2,sparsemat%nseg/),id='sparsemat%keyg')
      !!if (store_index) then
      !!    sparsemat%orb_from_index=f_malloc_ptr((/2,sparsemat%nvctr/),id='sparsemat%orb_from_index')
      !!end if
    end subroutine allocate_sparse_matrix_keys



    subroutine allocate_sparse_matrix_matrix_multiplication(nproc, norb, nseg, smmm)
      implicit none
      integer,intent(in) :: nproc, norb, nseg
      type(sparse_matrix_matrix_multiplication),intent(inout):: smmm
      smmm%ivectorindex=f_malloc_ptr(smmm%nseq,id='smmm%ivectorindex')
      smmm%ivectorindex_new=f_malloc_ptr(smmm%nseq,id='smmm%ivectorindex_new')
      smmm%onedimindices=f_malloc_ptr((/4,smmm%nout/),id='smmm%onedimindices')
      smmm%onedimindices_new=f_malloc_ptr((/4,smmm%nout/),id='smmm%onedimindices_new')
      !smmm%line_and_column_mm=f_malloc_ptr((/2,smmm%nvctrp_mm/),id='smmm%line_and_column_mm')
      smmm%nsegline=f_malloc_ptr(norb,id='smmm%nsegline')
      smmm%istsegline=f_malloc_ptr(norb,id='smmm%istsegline')
      smmm%indices_extract_sequential=f_malloc_ptr(smmm%nseq,id='smmm%indices_extract_sequential')
      smmm%nvctr_par=f_malloc_ptr(0.to.nproc-1,id='smmm%nvctr_par')
      smmm%isvctr_par=f_malloc_ptr(0.to.nproc-1,id='smmm%isvctr_par')
      smmm%nvctr_mm_par=f_malloc_ptr(0.to.nproc-1,id='smmm%nvctr_mm_par')
      smmm%isvctr_mm_par=f_malloc_ptr(0.to.nproc-1,id='smmm%isvctr_mm_par')
      smmm%keyv=f_malloc_ptr(nseg,id='smmm%keyv')
      smmm%keyg=f_malloc_ptr((/2,2,nseg/),id='smmm%keyg')
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
      ! smmm%nfvctrp is the number of columns per MPI task for an optimized load
      ! balancing during the sparse matrix matrix multiplications.
      if (sparsemat%smmm%nfvctrp>sparsemat%nfvctrp) then
          mat%matrixp = sparsematrix_malloc_ptr(sparsemat, iaction=DENSE_MATMUL, id=trim(matname)//'%matrixp')
      else
          mat%matrixp = sparsematrix_malloc_ptr(sparsemat, iaction=DENSE_PARALLEL, id=trim(matname)//'%matrixp')
      end if

      mat%power = 0.d0
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


    subroutine deallocate_sparse_matrix(sparsemat)
      use module_base 
      implicit none
      ! Calling arguments
      type(sparse_matrix),intent(inout):: sparsemat
      ! Local variables
      integer :: i, is, ie
      call f_free_ptr(sparseMat%keyg)
      call f_free_ptr(sparseMat%keyv)
      call f_free_ptr(sparseMat%nsegline)
      call f_free_ptr(sparseMat%istsegline)
      call f_free_ptr(sparseMat%matrixindex_in_compressed_fortransposed)
      call f_free_ptr(sparseMat%matrixindex_in_compressed_arr)
      call f_free_ptr(sparseMat%isvctr_par)
      call f_free_ptr(sparseMat%nvctr_par)
      call f_free_ptr(sparseMat%isfvctr_par)
      call f_free_ptr(sparseMat%nfvctr_par)
      !!if (associated(sparseMat%orb_from_index)) call f_free_ptr(sparseMat%orb_from_index)
      call f_free_ptr(sparseMat%taskgroup_startend)
      call f_free_ptr(sparseMat%taskgroupid)
      call deallocate_sparse_matrix_matrix_multiplication(sparsemat%smmm)
      call smat_release_mpi_groups(sparseMat%mpi_groups)
      call f_free_ptr(sparseMat%inwhichtaskgroup)
      call f_free_ptr(sparseMat%tgranks)
      call f_free_ptr(sparseMat%nranks)
      call f_free_ptr(sparseMat%luccomm)
      call f_free_ptr(sparseMat%on_which_atom)
    end subroutine deallocate_sparse_matrix

    subroutine smat_release_mpi_groups(mpi_groups)
      implicit none
      type(mpi_environment), dimension(:), pointer :: mpi_groups
      !local variables
      integer :: is,ie,i
      if (associated(mpi_groups)) then
         is=lbound(mpi_groups,1)
         ie=ubound(mpi_groups,1)
         do i=is,ie
            call release_mpi_environment(mpi_groups(i))
         end do
         deallocate(mpi_groups)
         nullify(mpi_groups)
      end if
    end subroutine smat_release_mpi_groups

    subroutine copy_sparse_matrix(smat_in, smat_out)
      use wrapper_MPI, only: mpi_environment_null,deepcopy_mpi_environment
      use dynamic_memory
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat_in
      type(sparse_matrix),intent(out) :: smat_out

      ! Local variables
      integer :: i, is, ie

      smat_out%nvctr = smat_in%nvctr
      smat_out%nseg = smat_in%nseg
      smat_out%nvctrp = smat_in%nvctrp
      smat_out%isvctr = smat_in%isvctr
      smat_out%parallel_compression = smat_in%parallel_compression
      smat_out%nfvctr = smat_in%nfvctr
      smat_out%nfvctrp = smat_in%nfvctrp
      smat_out%isfvctr = smat_in%isfvctr
      smat_out%nspin = smat_in%nspin
      smat_out%offset_matrixindex_in_compressed_fortransposed = smat_in%offset_matrixindex_in_compressed_fortransposed
      smat_out%store_index = smat_in%store_index
      smat_out%can_use_dense = smat_in%can_use_dense
      smat_out%ntaskgroup = smat_in%ntaskgroup
      smat_out%ntaskgroupp = smat_in%ntaskgroupp
      smat_out%istartendseg_t(1:2) = smat_in%istartendseg_t(1:2)
      smat_out%istartend_t(1:2) = smat_in%istartend_t(1:2)
      smat_out%nvctrp_tg = smat_in%nvctrp_tg
      smat_out%isvctrp_tg = smat_in%isvctrp_tg
      smat_out%iseseg_tg(1:2) = smat_in%iseseg_tg(1:2)
      smat_out%istartend_local(1:2) = smat_in%istartend_local(1:2)
      smat_out%istartendseg_local(1:2) = smat_in%istartendseg_local(1:2)

      smat_out%keyv=f_malloc_ptr(src_ptr=smat_in%keyv, id='smat_out%')
      smat_out%nsegline=f_malloc_ptr(src_ptr=smat_in%nsegline, id='smat_out%nsegline')
      smat_out%istsegline=f_malloc_ptr(src_ptr=smat_in%istsegline, id='smat_out%istsegline')
      smat_out%isvctr_par=f_malloc_ptr(src_ptr=smat_in%isvctr_par, id='smat_out%isvctr_par')
      smat_out%nvctr_par=f_malloc_ptr(src_ptr=smat_in%nvctr_par, id='smat_out%nvctr_par')
      smat_out%isfvctr_par=f_malloc_ptr(src_ptr=smat_in%isfvctr_par, id='smat_out%isfvctr_par')
      smat_out%nfvctr_par=f_malloc_ptr(src_ptr=smat_in%nfvctr_par, id='smat_out%nfvctr_par')
      smat_out%keyg=f_malloc_ptr(src_ptr=smat_in%keyg, id='smat_out%keyg')
      smat_out%matrixindex_in_compressed_arr=f_malloc_ptr(src_ptr=smat_in%matrixindex_in_compressed_arr, &
           id='smat_out%matrixindex_in_compressed_arr')
      smat_out%matrixindex_in_compressed_fortransposed=f_malloc_ptr(src_ptr=smat_in%matrixindex_in_compressed_fortransposed, &
           id='smat_out%matrixindex_in_compressed_fortransposed')
      smat_out%taskgroup_startend=f_malloc_ptr(src_ptr=smat_in%taskgroup_startend, id='smat_out%taskgroup_startend')
      smat_out%taskgroupid=f_malloc_ptr(src_ptr=smat_in%taskgroupid, id='smat_out%taskgroupid')
      smat_out%inwhichtaskgroup=f_malloc_ptr(src_ptr=smat_in%inwhichtaskgroup, id='smat_out%inwhichtaskgroup')
      smat_out%tgranks=f_malloc_ptr(src_ptr=smat_in%tgranks, id='smat_out%tgranks')
      smat_out%nranks=f_malloc_ptr(src_ptr=smat_in%nranks, id='smat_out%nranks')
      smat_out%luccomm=f_malloc_ptr(src_ptr=smat_in%luccomm, id='smat_out%luccomm')
      smat_out%on_which_atom=f_malloc_ptr(src_ptr=smat_in%on_which_atom, id='smat_out%on_which_atom')

      call copy_sparse_matrix_matrix_multiplication(smat_in%smmm, smat_out%smmm)
      call smat_release_mpi_groups(smat_out%mpi_groups)
      if (associated(smat_in%mpi_groups)) then
         is = lbound(smat_in%mpi_groups,1)
         ie = ubound(smat_in%mpi_groups,1)
         allocate(smat_out%mpi_groups(is:ie))
         do i=is,ie
            smat_out%mpi_groups(i) = mpi_environment_null()
            call deepcopy_mpi_environment(src=smat_in%mpi_groups(i), dest=smat_out%mpi_groups(i))
         end do
      else
         nullify(smat_out%mpi_groups)
      end if

    end subroutine copy_sparse_matrix


    subroutine copy_sparse_matrix_matrix_multiplication(smmm_in, smmm_out)
      !use copy_utils, only: allocate_and_copy
      use dynamic_memory
      implicit none

      ! Calling arguments
      type(sparse_matrix_matrix_multiplication),intent(in) :: smmm_in
      type(sparse_matrix_matrix_multiplication),intent(out) :: smmm_out
      smmm_out%nout = smmm_in%nout
      smmm_out%nseq = smmm_in%nseq
      smmm_out%nseg = smmm_in%nseg
      smmm_out%nfvctrp = smmm_in%nfvctrp
      smmm_out%isfvctr = smmm_in%isfvctr
      smmm_out%nvctrp_mm = smmm_in%nvctrp_mm
      smmm_out%isvctr_mm = smmm_in%isvctr_mm
      smmm_out%isseg = smmm_in%isseg
      smmm_out%ieseg = smmm_in%ieseg
      smmm_out%nccomm_smmm = smmm_in%nccomm_smmm
      smmm_out%nvctrp = smmm_in%nvctrp
      smmm_out%isvctr = smmm_in%isvctr
      smmm_out%nconsecutive_max = smmm_in%nconsecutive_max
      smmm_out%istartendseg_mm(1:2) = smmm_in%istartendseg_mm(1:2)
      smmm_out%istartend_mm(1:2) = smmm_in%istartend_mm(1:2)
      smmm_out%istartend_mm_dj(1:2) = smmm_in%istartend_mm_dj(1:2)

      !!converted  call allocate_and_copy(smmm_in%keyv, smmm_out%keyv, id='smmm_out%keyv')
      smmm_out%keyv=f_malloc_ptr(src_ptr=smmm_in%keyv, id='smmm_out%keyv')
      !!converted  call allocate_and_copy(smmm_in%keyg, smmm_out%keyg, id='smmm_out%keyg')
      smmm_out%keyg=f_malloc_ptr(src_ptr=smmm_in%keyg, id='smmm_out%keyg')
      !!converted  call allocate_and_copy(smmm_in%isvctr_mm_par, smmm_out%isvctr_mm_par, id='smmm_out%isvctr_mm_par')
      smmm_out%isvctr_mm_par=f_malloc_ptr(src_ptr=smmm_in%isvctr_mm_par, id='smmm_out%isvctr_mm_par')
      !!converted  call allocate_and_copy(smmm_in%nvctr_mm_par, smmm_out%nvctr_mm_par, id='smmm_out%nvctr_mm_par')
      smmm_out%nvctr_mm_par=f_malloc_ptr(src_ptr=smmm_in%nvctr_mm_par, id='smmm_out%nvctr_mm_par')
      !!converted  call allocate_and_copy(smmm_in%ivectorindex, smmm_out%ivectorindex, id='smmm_out%ivectorindex')
      smmm_out%ivectorindex=f_malloc_ptr(src_ptr=smmm_in%ivectorindex, id='smmm_out%ivectorindex')
      !!converted  call allocate_and_copy(smmm_in%ivectorindex_new, smmm_out%ivectorindex_new, id='smmm_out%ivectorindex_new')
      smmm_out%ivectorindex_new=f_malloc_ptr(src_ptr=smmm_in%ivectorindex_new, id='smmm_out%ivectorindex_new')
      !!converted  call allocate_and_copy(smmm_in%nsegline, smmm_out%nsegline, id='smmm_out%segline')
      smmm_out%nsegline=f_malloc_ptr(src_ptr=smmm_in%nsegline, id='smmm_out%segline')
      !!converted  call allocate_and_copy(smmm_in%istsegline, smmm_out%istsegline, id='smmm_out%istsegline')
      smmm_out%istsegline=f_malloc_ptr(src_ptr=smmm_in%istsegline, id='smmm_out%istsegline')
      !!converted  call allocate_and_copy(smmm_in%indices_extract_sequential, smmm_out%indices_extract_sequential, &
      smmm_out%indices_extract_sequential=f_malloc_ptr(src_ptr=smmm_in%indices_extract_sequential, &
           id='smmm_out%ndices_extract_sequential')
      !!converted  call allocate_and_copy(smmm_in%onedimindices, smmm_out%onedimindices, id='smmm_out%onedimindices')
      smmm_out%onedimindices=f_malloc_ptr(src_ptr=smmm_in%onedimindices, id='smmm_out%onedimindices')
      !!converted  call allocate_and_copy(smmm_in%onedimindices_new, smmm_out%onedimindices_new, id='smmm_out%onedimindices_new')
      smmm_out%onedimindices_new=f_malloc_ptr(src_ptr=smmm_in%onedimindices_new, id='smmm_out%onedimindices_new')
      !!converted  call allocate_and_copy(smmm_in%line_and_column_mm, smmm_out%line_and_column_mm, id='smmm_out%line_and_column_mm')
      smmm_out%line_and_column_mm=f_malloc_ptr(src_ptr=smmm_in%line_and_column_mm, id='smmm_out%line_and_column_mm')
      !!converted  call allocate_and_copy(smmm_in%line_and_column, smmm_out%line_and_column, id='smmm_out%line_and_column')
      smmm_out%line_and_column=f_malloc_ptr(src_ptr=smmm_in%line_and_column, id='smmm_out%line_and_column')
      !!converted  call allocate_and_copy(smmm_in%luccomm_smmm, smmm_out%luccomm_smmm, id='smmm_out%luccomm_smmm')
      smmm_out%luccomm_smmm=f_malloc_ptr(src_ptr=smmm_in%luccomm_smmm, id='smmm_out%luccomm_smmm')
      !!converted  call allocate_and_copy(smmm_in%isvctr_par, smmm_out%isvctr_par, id='smmm_out%isvctr_par')
      smmm_out%isvctr_par=f_malloc_ptr(src_ptr=smmm_in%isvctr_par, id='smmm_out%isvctr_par')
      !!converted  call allocate_and_copy(smmm_in%nvctr_par, smmm_out%nvctr_par, id='smmm_out%nvctr_par')
      smmm_out%nvctr_par=f_malloc_ptr(src_ptr=smmm_in%nvctr_par, id='smmm_out%nvctr_par')
      !!converted  call allocate_and_copy(smmm_in%consecutive_lookup, smmm_out%consecutive_lookup, id='smmm_out%consecutive_lookup')
      smmm_out%consecutive_lookup=f_malloc_ptr(src_ptr=smmm_in%consecutive_lookup, id='smmm_out%consecutive_lookup')
      !!call allocate_and_copy(smmm_in%keyg, smmm_out%keyg, id='smmm_out%keyg')
    end subroutine copy_sparse_matrix_matrix_multiplication

    subroutine copy_matrices(mat_in, mat_out)
      !use copy_utils, only: allocate_and_copy
      use dynamic_memory
      implicit none

      ! Calling arguments
      type(matrices),intent(in) :: mat_in
      type(matrices),intent(out) :: mat_out


      mat_out%matrix_compr=f_malloc_ptr(src_ptr=mat_in%matrix_compr, id='mat_out%matrix_compr')
      mat_out%matrix_comprp=f_malloc_ptr(src_ptr=mat_in%matrix_comprp, id='mat_out%matrix_comprp')
      mat_out%matrix=f_malloc_ptr(src_ptr=mat_in%matrix, id='mat_out%matrix')
      mat_out%matrixp=f_malloc_ptr(src_ptr=mat_in%matrixp, id='mat_out%matrixp')

    end subroutine copy_matrices


    subroutine deallocate_sparse_matrix_matrix_multiplication(smmm)
      implicit none
      type(sparse_matrix_matrix_multiplication),intent(out):: smmm
      call f_free_ptr(smmm%ivectorindex)
      call f_free_ptr(smmm%ivectorindex_new)
      call f_free_ptr(smmm%onedimindices)
      call f_free_ptr(smmm%onedimindices_new)
      call f_free_ptr(smmm%line_and_column_mm)
      call f_free_ptr(smmm%line_and_column)
      call f_free_ptr(smmm%nsegline)
      call f_free_ptr(smmm%istsegline)
      !!call f_free_ptr(smmm%keyg)
      call f_free_ptr(smmm%indices_extract_sequential)
      call f_free_ptr(smmm%nvctr_par)
      call f_free_ptr(smmm%isvctr_par)
      call f_free_ptr(smmm%nvctr_mm_par)
      call f_free_ptr(smmm%isvctr_mm_par)
      call f_free_ptr(smmm%luccomm_smmm)
      call f_free_ptr(smmm%keyv)
      call f_free_ptr(smmm%keyg)
      call f_free_ptr(smmm%consecutive_lookup)
    end subroutine deallocate_sparse_matrix_matrix_multiplication

    subroutine allocate_smat_d1_ptr(smat_ptr,smat_info_ptr)
      implicit none
      double precision,dimension(:),pointer,intent(inout) :: smat_ptr
      type(sparse_matrix_info_ptr), intent(in) :: smat_info_ptr

      select case (smat_info_ptr%iaction)
      case (SPARSE_TASKGROUP)
          smat_ptr = f_malloc_ptr(max(1,smat_info_ptr%smat%nvctrp_tg*smat_info_ptr%smat%nspin),id=smat_info_ptr%id)
      case (SPARSE_FULL)
          smat_ptr = f_malloc_ptr(smat_info_ptr%smat%nvctr*smat_info_ptr%smat%nspin,id=smat_info_ptr%id)
      case (SPARSE_PARALLEL)
          smat_ptr = f_malloc_ptr(max(1,smat_info_ptr%smat%nvctrp*smat_info_ptr%smat%nspin),id=smat_info_ptr%id)
      case (SPARSEMM_SEQ)
          smat_ptr = f_malloc_ptr(max(1,smat_info_ptr%smat%smmm%nseq*smat_info_ptr%smat%nspin),id=smat_info_ptr%id)
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
          smat_ptr = f_malloc_ptr((/smat_info_ptr%smat%nfvctr,smat_info_ptr%smat%nfvctrp/),id=smat_info_ptr%id)
      case (DENSE_MATMUL)
          smat_ptr = f_malloc_ptr((/smat_info_ptr%smat%nfvctr,smat_info_ptr%smat%smmm%nfvctrp/),id=smat_info_ptr%id)
      case default
         call f_err_throw('The action specified for the 2d matrix allocation is invalid',&
              err_name='BIGDFT_RUNTIME_ERROR')
      end select
    end subroutine allocate_smat_d2_ptr
    

    subroutine allocate_smat_d3_ptr(smat_ptr,smat_info_ptr)
      implicit none
      double precision,dimension(:,:,:),pointer,intent(inout) :: smat_ptr
      type(sparse_matrix_info_ptr), intent(in) :: smat_info_ptr

      select case (smat_info_ptr%iaction)
      case (DENSE_FULL)
          smat_ptr = f_malloc_ptr((/smat_info_ptr%smat%nfvctr,smat_info_ptr%smat%nfvctr,smat_info_ptr%smat%nspin/), &
                                  id=smat_info_ptr%id)
      case (DENSE_PARALLEL)
          smat_ptr = f_malloc_ptr((/smat_info_ptr%smat%nfvctr,smat_info_ptr%smat%nfvctrp,smat_info_ptr%smat%nspin/),&
                                  id=smat_info_ptr%id)
      case (DENSE_MATMUL)
          smat_ptr = f_malloc_ptr((/smat_info_ptr%smat%nfvctr,smat_info_ptr%smat%smmm%nfvctrp,smat_info_ptr%smat%nspin/),&
                                  id=smat_info_ptr%id)
      case default
         call f_err_throw('The action specified for the 3d matrix allocation is invalid',&
              err_name='BIGDFT_RUNTIME_ERROR')
      end select
    end subroutine allocate_smat_d3_ptr


    subroutine allocate_smat_d1(smat,smat_info)
      implicit none
      double precision,dimension(:),allocatable,intent(inout) :: smat
      type(sparse_matrix_info), intent(in) :: smat_info

      select case (smat_info%iaction)
      case (SPARSE_TASKGROUP)
          smat = f_malloc(max(1,smat_info%smat%nvctrp_tg*smat_info%smat%nspin),id=smat_info%id)
      case (SPARSE_FULL)
          smat = f_malloc(smat_info%smat%nvctr*smat_info%smat%nspin,id=smat_info%id)
      case (SPARSE_PARALLEL)
          smat = f_malloc(max(1,smat_info%smat%nvctrp*smat_info%smat%nspin),id=smat_info%id)
      case (SPARSEMM_SEQ)
          smat = f_malloc(max(1,smat_info%smat%smmm%nseq*smat_info%smat%nspin),id=smat_info%id)
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
          smat = f_malloc((/smat_info%smat%nfvctr,smat_info%smat%nfvctrp/),id=smat_info%id)
      case (DENSE_MATMUL)
          smat = f_malloc((/smat_info%smat%nfvctr,smat_info%smat%smmm%nfvctrp/),id=smat_info%id)
      case default
         call f_err_throw('The action specified for the 2d matrix allocation is invalid',&
              err_name='BIGDFT_RUNTIME_ERROR')
      end select
    end subroutine allocate_smat_d2


    subroutine allocate_smat_d3(smat,smat_info)
      implicit none
      double precision,dimension(:,:,:),allocatable,intent(inout) :: smat
      type(sparse_matrix_info), intent(in) :: smat_info

      select case (smat_info%iaction)
      case (DENSE_FULL)
          smat = f_malloc((/smat_info%smat%nfvctr,smat_info%smat%nfvctr,smat_info%smat%nspin/),id=smat_info%id)
      case (DENSE_PARALLEL)
          smat = f_malloc((/smat_info%smat%nfvctr,smat_info%smat%nfvctrp,smat_info%smat%nspin/),id=smat_info%id)
      case (DENSE_MATMUL)
          smat = f_malloc((/smat_info%smat%nfvctr,smat_info%smat%smmm%nfvctrp,smat_info%smat%nspin/),id=smat_info%id)
      case default
         call f_err_throw('The action specified for the 3d matrix allocation is invalid',&
              err_name='BIGDFT_RUNTIME_ERROR')
      end select
    end subroutine allocate_smat_d3


    subroutine allocate0_smat_d1_ptr(smat_ptr,smat_info0_ptr)
      implicit none
      double precision,dimension(:),pointer,intent(inout) :: smat_ptr
      type(sparse_matrix_info0_ptr), intent(in) :: smat_info0_ptr

      select case (smat_info0_ptr%iaction)
      case (SPARSE_TASKGROUP)
          smat_ptr = f_malloc0_ptr(max(1,smat_info0_ptr%smat%nvctrp_tg*smat_info0_ptr%smat%nspin),id=smat_info0_ptr%id)
      case (SPARSE_FULL)
          smat_ptr = f_malloc0_ptr(smat_info0_ptr%smat%nvctr*smat_info0_ptr%smat%nspin,id=smat_info0_ptr%id)
      case (SPARSE_PARALLEL)
          smat_ptr = f_malloc0_ptr(max(1,smat_info0_ptr%smat%nvctrp*smat_info0_ptr%smat%nspin),id=smat_info0_ptr%id)
      case (SPARSEMM_SEQ)
          smat_ptr = f_malloc0_ptr(max(1,smat_info0_ptr%smat%smmm%nseq*smat_info0_ptr%smat%nspin),id=smat_info0_ptr%id)
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
      case (DENSE_PARALLEL)
          smat_ptr = f_malloc0_ptr((/smat_info0_ptr%smat%nfvctr,smat_info0_ptr%smat%nfvctrp/),id=smat_info0_ptr%id)
      case (DENSE_MATMUL)
          smat_ptr = f_malloc0_ptr((/smat_info0_ptr%smat%nfvctr,smat_info0_ptr%smat%smmm%nfvctrp/),id=smat_info0_ptr%id)
      case default
         call f_err_throw('The action specified for the 2d matrix allocation is invalid',&
              err_name='BIGDFT_RUNTIME_ERROR')
      end select
    end subroutine allocate0_smat_d2_ptr
    

    subroutine allocate0_smat_d3_ptr(smat_ptr,smat_info0_ptr)
      implicit none
      double precision,dimension(:,:,:),pointer,intent(inout) :: smat_ptr
      type(sparse_matrix_info0_ptr), intent(in) :: smat_info0_ptr

      select case (smat_info0_ptr%iaction)
      case (DENSE_FULL)
          smat_ptr = f_malloc0_ptr((/smat_info0_ptr%smat%nfvctr,smat_info0_ptr%smat%nfvctr,smat_info0_ptr%smat%nspin/), &
                                   id=smat_info0_ptr%id)
      case (DENSE_PARALLEL)
          smat_ptr = f_malloc0_ptr((/smat_info0_ptr%smat%nfvctr,smat_info0_ptr%smat%nfvctrp,smat_info0_ptr%smat%nspin/), &
                                   id=smat_info0_ptr%id)
      case (DENSE_MATMUL)
          smat_ptr = f_malloc0_ptr((/smat_info0_ptr%smat%nfvctr,smat_info0_ptr%smat%smmm%nfvctrp,smat_info0_ptr%smat%nspin/), &
                                   id=smat_info0_ptr%id)
      case default
         call f_err_throw('The action specified for the 2d matrix allocation is invalid',&
              err_name='BIGDFT_RUNTIME_ERROR')
      end select
    end subroutine allocate0_smat_d3_ptr


    subroutine allocate0_smat_d1(smat,smat_info0)
      implicit none
      double precision,dimension(:),allocatable,intent(inout) :: smat
      type(sparse_matrix_info0), intent(in) :: smat_info0

      select case (smat_info0%iaction)
      case (SPARSE_TASKGROUP)
          smat = f_malloc0(max(1,smat_info0%smat%nvctrp_tg*smat_info0%smat%nspin),id=smat_info0%id)
      case (SPARSE_FULL)
          smat = f_malloc0(smat_info0%smat%nvctr*smat_info0%smat%nspin,id=smat_info0%id)
      case (SPARSE_PARALLEL)
          smat = f_malloc0(max(1,smat_info0%smat%nvctrp*smat_info0%smat%nspin),id=smat_info0%id)
      case (SPARSEMM_SEQ)
          smat = f_malloc0(max(1,smat_info0%smat%smmm%nseq*smat_info0%smat%nspin),id=smat_info0%id)
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
      case (DENSE_PARALLEL)
          smat = f_malloc0((/smat_info0%smat%nfvctr,smat_info0%smat%nfvctrp/),id=smat_info0%id)
      case (DENSE_MATMUL)
          smat = f_malloc0((/smat_info0%smat%nfvctr,smat_info0%smat%smmm%nfvctrp/),id=smat_info0%id)
      case default
          call f_err_throw('The action specified for the 2d matrix allocation is invalid',&
               err_name='BIGDFT_RUNTIME_ERROR')
      end select
    end subroutine allocate0_smat_d2


    subroutine allocate0_smat_d3(smat,smat_info0)
      implicit none
      double precision,dimension(:,:,:),allocatable,intent(inout) :: smat
      type(sparse_matrix_info0), intent(in) :: smat_info0

      select case (smat_info0%iaction)
      case (DENSE_FULL)
          smat = f_malloc0((/smat_info0%smat%nfvctr,smat_info0%smat%nfvctr,smat_info0%smat%nspin/),id=smat_info0%id)
      case (DENSE_PARALLEL)
          smat = f_malloc0((/smat_info0%smat%nfvctr,smat_info0%smat%nfvctrp,smat_info0%smat%nspin/),id=smat_info0%id)
      case (DENSE_MATMUL)
          smat = f_malloc0((/smat_info0%smat%nfvctr,smat_info0%smat%smmm%nfvctrp,smat_info0%smat%nspin/),id=smat_info0%id)
      case default
          call f_err_throw('The action specified for the 2d matrix allocation is invalid',&
               err_name='BIGDFT_RUNTIME_ERROR')
      end select
    end subroutine allocate0_smat_d3


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
