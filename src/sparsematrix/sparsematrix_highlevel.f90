module sparsematrix_highlevel

  private

  public :: sparse_matrix_and_matrices_init_from_file_ccs
  public :: sparse_matrix_init_from_file_ccs
  public :: sparse_matrix_init_from_data
  public :: matrices_init
  public :: matrices_set

  contains

    subroutine sparse_matrix_and_matrices_init_from_file_ccs(filename, iproc, nproc, smat, mat)
      use module_base
      use sparsematrix_base, only: sparse_matrix, matrices
      use sparsematrix_init, only: read_ccs_format
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      character(len=*),intent(in) :: filename
      type(sparse_matrix),intent(out) :: smat
      type(matrices),intent(out) :: mat
    
      ! Local variables
      integer :: nfvctr, nvctr
      integer,dimension(:),pointer :: col_ptr, row_ind
      real(kind=8),dimension(:),pointer :: val
    
      call f_routine(id='sparse_matrix_and_matrices_init_from_file_ccs')
    
      ! Read in the matrix
      call read_ccs_format(filename, nfvctr, nvctr, col_ptr, row_ind, val)
    
      ! Generate the sparse_matrix type
      call sparse_matrix_init_from_data(iproc, nproc, nfvctr, nvctr, row_ind, col_ptr, smat)
    
      ! Generate the matrices type
      call matrices_init_from_data(smat, val, mat)
    
      ! Deallocate the pointers
      call f_free_ptr(col_ptr)
      call f_free_ptr(row_ind)
      call f_free_ptr(val)
    
      call f_release_routine()
    
    end subroutine sparse_matrix_and_matrices_init_from_file_ccs
    
    
    subroutine sparse_matrix_init_from_file_ccs(filename, iproc, nproc, smat)
      use module_base
      use sparsematrix_base, only: sparse_matrix
      use sparsematrix_init, only: read_ccs_format
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      character(len=*),intent(in) :: filename
      type(sparse_matrix),intent(out) :: smat
    
      ! Local variables
      integer :: nfvctr, nvctr
      integer,dimension(:),pointer :: col_ptr, row_ind
      real(kind=8),dimension(:),pointer :: val
    
      call f_routine(id='sparse_matrix_and_matrices_init_from_file_ccs')
    
      ! Read in the matrix
      call read_ccs_format(filename, nfvctr, nvctr, col_ptr, row_ind, val)
    
      ! Generate the sparse_matrix type
      call sparse_matrix_init_from_data(iproc, nproc, nfvctr, nvctr, row_ind, col_ptr, smat)
    
      ! Deallocate the pointers
      call f_free_ptr(col_ptr)
      call f_free_ptr(row_ind)
      call f_free_ptr(val)
    
      call f_release_routine()
    
    end subroutine sparse_matrix_init_from_file_ccs
    
    
    subroutine sparse_matrix_init_from_data(iproc, nproc, nfvctr, nvctr, row_ind, col_ptr, smat)
      use module_base
      use sparsematrix_base, only: sparse_matrix
      use sparsematrix_init, only: ccs_to_sparsebigdft_short, &
                                   bigdft_to_sparsebigdft, init_matrix_taskgroups
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nfvctr, nvctr
      integer,dimension(nvctr),intent(in) :: row_ind
      integer,dimension(nfvctr),intent(in) :: col_ptr
      type(sparse_matrix),intent(out) :: smat
    
      ! Local variables
      integer :: nseg
      integer,dimension(:),pointer :: keyv
      integer,dimension(:,:,:),pointer :: keyg
    
      call f_routine(id='sparse_matrix_init_from_data')
    
      ! Convert the sparsity pattern to the BigDFT format
      call ccs_to_sparsebigdft_short(nfvctr, nvctr, row_ind, col_ptr, nseg, keyv, keyg)
    
      ! Create the sparse_matrix structure
      call bigdft_to_sparsebigdft(iproc, nproc, nfvctr, nvctr, nseg, keyg, smat)
    
      ! Deallocate the pointers
      call f_free_ptr(keyv)
      call f_free_ptr(keyg)
    
      call f_release_routine()
    
    end subroutine sparse_matrix_init_from_data
    
    
    subroutine matrices_init_from_data(smat, val, mat)
      use module_base
      use sparsematrix_base, only: sparse_matrix, matrices, matrices_null, &
                                   assignment(=), sparsematrix_malloc_ptr, SPARSE_FULL
      implicit none
    
      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      real(kind=8),dimension(smat%nvctr),intent(in) :: val
      type(matrices),intent(out) :: mat
    
      call f_routine(id='matrices_init_from_data')
    
      ! Create the matrices structure
      mat = matrices_null()
      mat%matrix_compr = sparsematrix_malloc_ptr(smat, iaction=SPARSE_FULL, id='mat%matrix_compr')
    
      ! Copy the content
      call f_memcpy(src=val, dest=mat%matrix_compr)
    
      call f_release_routine()
    
    end subroutine matrices_init_from_data
    
    
    subroutine matrices_init(smat, mat)
      use module_base
      use sparsematrix_base,only: sparse_matrix, matrices, &
                                  matrices_null, assignment(=), sparsematrix_malloc_ptr, SPARSE_FULL
      implicit none
    
      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(out) :: mat
    
      mat = matrices_null()
      mat%matrix_compr = sparsematrix_malloc_ptr(smat, iaction=SPARSE_FULL, id='mat%matrix_compr')
    
    end subroutine matrices_init
    
    
    subroutine matrices_set(smat, val, mat)
      use module_base
      use sparsematrix_base,only: sparse_matrix, matrices, &
                                  matrices_null, assignment(=), sparsematrix_malloc_ptr, SPARSE_FULL
      implicit none
    
      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      real(kind=8),dimension(:) :: val
      type(matrices),intent(inout) :: mat
    
      call f_routine(id='matrices_set')
    
      if (size(val)/=smat%nvctr) then
          call f_err_throw('The size of the array used to set the matrix contents is wrong: '&
               &//trim(yaml_toa(size(val)))//' instead of '//trim(yaml_toa(smat%nvctr)))
      end if
      call f_memcpy(src=val, dest=mat%matrix_compr)
    
      call f_release_routine()
    
    end subroutine matrices_set

end module sparsematrix_highlevel
