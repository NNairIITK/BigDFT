module sparsematrix_highlevel

  private

  public :: sparse_matrix_and_matrices_init_from_file_ccs
  public :: sparse_matrix_init_from_file_ccs
  public :: sparse_matrix_init_from_data
  public :: matrices_init
  public :: matrices_set
  public :: ccs_data_from_sparse_matrix
  public :: ccs_matrix_write
  public :: matrix_matrix_multiplication
  public :: matrix_chebyshev_expansion

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


    subroutine ccs_data_from_sparse_matrix(smat, row_ind, col_ptr)
      use module_base
      use sparsematrix_base, only: sparse_matrix
      use sparsematrix_init, only: sparsebigdft_to_ccs
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      integer,dimension(:),pointer,intent(inout) :: row_ind, col_ptr

      if (size(row_ind)/=smat%nvctr) then
          call f_err_throw('The size of the array row_ind is wrong: '&
               &//trim(yaml_toa(size(row_ind)))//' instead of '//trim(yaml_toa(smat%nvctr)))
      end if
      if (size(col_ptr)/=smat%nfvctr) then
          call f_err_throw('The size of the array col_ptr is wrong: '&
               &//trim(yaml_toa(size(col_ptr)))//' instead of '//trim(yaml_toa(smat%nfvctr)))
      end if
      call sparsebigdft_to_ccs(smat%nfvctr, smat%nvctr, smat%nseg, smat%keyg, row_ind, col_ptr)
    end subroutine ccs_data_from_sparse_matrix


    subroutine ccs_matrix_write(filename, smat, row_ind, col_ptr, mat)
      use module_base
      use sparsematrix_base, only: sparse_matrix, matrices
      use io, only: write_ccs_matrix
      implicit none

      ! Calling arguments
      character(len=*),intent(in) :: filename
      type(sparse_matrix),intent(in) :: smat
      integer,dimension(:),pointer,intent(inout) :: row_ind, col_ptr
      type(matrices),intent(in) :: mat

      if (size(row_ind)/=smat%nvctr) then
          call f_err_throw('The size of the array row_ind is wrong: '&
               &//trim(yaml_toa(size(row_ind)))//' instead of '//trim(yaml_toa(smat%nvctr)))
      end if
      if (size(col_ptr)/=smat%nfvctr) then
          call f_err_throw('The size of the array col_ptr is wrong: '&
               &//trim(yaml_toa(size(col_ptr)))//' instead of '//trim(yaml_toa(smat%nfvctr)))
      end if

      call write_ccs_matrix(filename, smat%nfvctr, smat%nvctr, row_ind, col_ptr, mat%matrix_compr)

    end subroutine ccs_matrix_write


    subroutine matrix_matrix_multiplication(iproc, nproc, smat, a, b, c)
      use module_base
      use sparsematrix_base, only: sparse_matrix, matrices
      use sparsematrix, only: matrix_matrix_mult_wrapper
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(inout) :: a
      type(matrices),intent(inout) :: b, c !b actually also in... 

      call f_routine(id='matrix_matrix_multiplication')

      if (size(a%matrix_compr)/=smat%nvctrp_tg) then
          call f_err_throw('The size of the array a%matrix_compr is wrong: '&
               &//trim(yaml_toa(size(a%matrix_compr)))//' instead of '//trim(yaml_toa(smat%nvctrp_tg)))
      end if
      if (size(b%matrix_compr)/=smat%nvctrp_tg) then
          call f_err_throw('The size of the array b%matrix_compr is wrong: '&
               &//trim(yaml_toa(size(b%matrix_compr)))//' instead of '//trim(yaml_toa(smat%nvctrp_tg)))
      end if
      if (size(c%matrix_compr)/=smat%nvctrp_tg) then
          call f_err_throw('The size of the array c%matrix_compr is wrong: '&
               &//trim(yaml_toa(size(c%matrix_compr)))//' instead of '//trim(yaml_toa(smat%nvctrp_tg)))
      end if
      call matrix_matrix_mult_wrapper(iproc, nproc, smat, a%matrix_compr, b%matrix_compr, c%matrix_compr)

      call f_release_routine()

    end subroutine matrix_matrix_multiplication


    subroutine matrix_chebyshev_expansion(iproc, nproc, ndeg, ncalc, ex, &
               smat_in, smat_out, mat_in, mat_out)
      use module_base
      use sparsematrix_base, only: sparse_matrix, matrices
      use ice, only: inverse_chebyshev_expansion
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, ndeg, ncalc
      type(sparse_matrix),intent(in) ::smat_in, smat_out
      real(kind=8),dimension(ncalc),intent(in) :: ex
      type(matrices),intent(in) :: mat_in
      type(matrices),dimension(ncalc),intent(inout) :: mat_out

      ! Local variables
      integer :: i

      call f_routine(id='matrix_chebyshev_expansion')

      ! Check the dimensions of the internal arrays
      if (size(mat_in%matrix_compr)/=smat_in%nvctrp_tg) then
          call f_err_throw('The size of the array mat_in%matrix_compr is wrong: '&
               &//trim(yaml_toa(size(mat_in%matrix_compr)))//' instead of '//trim(yaml_toa(smat_in%nvctrp_tg)))
      end if
      do i=1,ncalc
          if (size(mat_out(i)%matrix_compr)/=smat_out%nvctrp_tg) then
              call f_err_throw('The size of the array mat_out%matrix_compr is wrong: '&
                   &//trim(yaml_toa(size(mat_out(i)%matrix_compr)))//' instead of '//trim(yaml_toa(smat_out%nvctrp_tg)))
          end if
      end do

      ! Check that number of non-zero elements of smat_in is not larger than that of smat_out. This is just a minimal
      ! check to ensure that the sparsity pattern of smat_in is contained within the one of smat_out. However this check
      ! is not sufficient to make sure that this condition is fulfilled.
      if (smat_in%nvctr>smat_out%nvctr) then
          call f_err_throw('The number of non-zero elements of smat_in ('//&
               trim(yaml_toa(smat_in%nvctr))//') is larger than the one of smat_out ('//&
               trim(yaml_toa(smat_out%nvctr))//')')
      end if

      call inverse_chebyshev_expansion(iproc, nproc, ndeg, &
           smat_in, smat_out, ncalc, ex, mat_in, mat_out)

      call f_release_routine()

    end subroutine matrix_chebyshev_expansion

end module sparsematrix_highlevel
