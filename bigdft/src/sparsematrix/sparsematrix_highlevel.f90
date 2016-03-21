module sparsematrix_highlevel

  private

  public :: sparse_matrix_and_matrices_init_from_file_ccs
  public :: sparse_matrix_init_from_file_ccs
  public :: sparse_matrix_init_from_data_ccs
  public :: sparse_matrix_and_matrices_init_from_file_bigdft
  public :: sparse_matrix_init_from_file_bigdft
  public :: sparse_matrix_metadata_init_from_file
  public :: read_sparse_matrix_metadata
  public :: sparse_matrix_init_from_data_bigdft
  public :: matrices_init
  public :: matrices_set_values
  public :: matrices_get_values
  public :: matrices_get_size
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
      call sparse_matrix_init_from_data_ccs(iproc, nproc, nfvctr, nvctr, row_ind, col_ptr, smat)

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
      call sparse_matrix_init_from_data_ccs(iproc, nproc, nfvctr, nvctr, row_ind, col_ptr, smat)

      ! Deallocate the pointers
      call f_free_ptr(col_ptr)
      call f_free_ptr(row_ind)
      call f_free_ptr(val)

      call f_release_routine()

    end subroutine sparse_matrix_init_from_file_ccs


    subroutine sparse_matrix_init_from_data_ccs(iproc, nproc, nfvctr, nvctr, row_ind, col_ptr, smat)
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

      call f_routine(id='sparse_matrix_init_from_data_ccs')

      ! Convert the sparsity pattern to the BigDFT format
      call ccs_to_sparsebigdft_short(nfvctr, nvctr, row_ind, col_ptr, nseg, keyv, keyg)

      ! Create the sparse_matrix structure
      call bigdft_to_sparsebigdft(iproc, nproc, nfvctr, nvctr, nseg, keyg, smat)

      ! Deallocate the pointers
      call f_free_ptr(keyv)
      call f_free_ptr(keyg)

      call f_release_routine()

    end subroutine sparse_matrix_init_from_data_ccs


    subroutine sparse_matrix_and_matrices_init_from_file_bigdft(filename, iproc, nproc, smat, mat, &
               init_matmul)
      use module_base
      use sparsematrix_base, only: sparse_matrix, matrices
      use sparsematrix_init, only: bigdft_to_sparsebigdft
      use sparsematrix_io, only: read_sparse_matrix
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      character(len=*),intent(in) :: filename
      type(sparse_matrix),intent(out) :: smat
      type(matrices),intent(out) :: mat
      logical,intent(in),optional :: init_matmul
      ! Optional variables that are contained within the sparse matrix format
      !!integer,intent(out),optional :: nat, ntypes
      !!integer,dimension(:),pointer,intent(inout),optional :: nzatom, nelpsp, iatype
      !!character(len=20),dimension(:),pointer,intent(inout),optional :: atomnames
      !!real(kind=8),dimension(:,:),pointer,intent(inout),optional :: rxyz
      !!integer,dimension(:),pointer,intent(inout),optional :: on_which_atom

      ! Local variables
      integer :: nspin, nfvctr, nseg, nvctr
!      character(len=1) :: geocode
      integer,dimension(:),pointer :: keyv
      integer,dimension(:,:,:),pointer :: keyg
      real(kind=8),dimension(:),pointer :: val
      logical :: init_matmul_
!      integer,dimension(:),pointer :: on_which_atom_
!      real(kind=8),dimension(3) :: cell_dim

      call f_routine(id='sparse_matrix_and_matrices_init_from_file_bigdft')

      ! Read in the matrix
      call read_sparse_matrix(filename, nspin, nfvctr, nseg, nvctr, keyv, keyg, val)


      if (present(init_matmul)) then
          init_matmul_ = init_matmul
      else
          init_matmul_ = .true.
      end if

      ! Create the sparse_matrix structure
      call bigdft_to_sparsebigdft(iproc, nproc, nfvctr, nvctr, nseg, keyg, smat, &
           init_matmul=init_matmul_)!, nspin=nspin, geocode=geocode, cell_dim=cell_dim, on_which_atom=on_which_atom_)

      ! Generate the matrices type
      call matrices_init_from_data(smat, val, mat)

      ! Deallocate the pointers
      call f_free_ptr(keyv)
      call f_free_ptr(keyg)
      call f_free_ptr(val)


      call f_release_routine()

    end subroutine sparse_matrix_and_matrices_init_from_file_bigdft


    subroutine sparse_matrix_init_from_file_bigdft(filename, iproc, nproc, smat)
      use module_base
      use sparsematrix_base, only: sparse_matrix
      use sparsematrix_init, only: bigdft_to_sparsebigdft
      use sparsematrix_io, only: read_sparse_matrix
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      character(len=*),intent(in) :: filename
      type(sparse_matrix),intent(out) :: smat

      ! Local variables
      integer :: nspin, nfvctr, nseg, nvctr
      character(len=1) :: geocode
      integer,dimension(:),pointer :: keyv
      integer,dimension(:,:,:),pointer :: keyg
      real(kind=8),dimension(:),pointer :: val
      real(kind=8),dimension(3) :: cell_dim
      integer,dimension(:),pointer :: on_which_atom

      call f_routine(id='sparse_matrix_and_matrices_init_from_file_ccs')

      ! Read in the matrix
      call read_sparse_matrix(filename, nspin, nfvctr, nseg, nvctr, keyv, keyg, val)

      ! Create the sparse_matrix structure
      call bigdft_to_sparsebigdft(iproc, nproc, nfvctr, nvctr, nseg, keyg, smat, &
           nspin=nspin, geocode=geocode, cell_dim=cell_dim, on_which_atom=on_which_atom)

      ! Deallocate the pointers
      call f_free_ptr(keyv)
      call f_free_ptr(keyg)
      call f_free_ptr(val)

      call f_release_routine()

    end subroutine sparse_matrix_init_from_file_bigdft


    subroutine sparse_matrix_metadata_init_from_file(filename, smmd)
      use module_base
      use sparsematrix_base, only: sparse_matrix_metadata, sparse_matrix_metadata_null
      use sparsematrix_io, only: read_sparse_matrix_metadata
      implicit none

      ! Calling arguments
      character(len=*),intent(in) :: filename
      type(sparse_matrix_metadata),intent(out) :: smmd

      smmd = sparse_matrix_metadata_null()
      call read_sparse_matrix_metadata(filename, smmd%nfvctr, smmd%nat, smmd%ntypes, &
           smmd%units, smmd%geocode, smmd%cell_dim, smmd%nzatom, smmd%nelpsp, &
           smmd%atomnames, smmd%iatype, smmd%rxyz, smmd%on_which_atom)

    end subroutine sparse_matrix_metadata_init_from_file


    subroutine sparse_matrix_init_from_data_bigdft(iproc, nproc, nfvctr, nvctr, nseg, keyg, smat)
      use module_base
      use sparsematrix_base, only: sparse_matrix
      use sparsematrix_init, only: ccs_to_sparsebigdft_short, &
                                   bigdft_to_sparsebigdft, init_matrix_taskgroups
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nfvctr, nvctr, nseg
      integer,dimension(2,2,nseg),intent(in) :: keyg
      type(sparse_matrix),intent(out) :: smat

      call f_routine(id='sparse_matrix_init_from_data_ccs')

      ! Create the sparse_matrix structure
      call bigdft_to_sparsebigdft(iproc, nproc, nfvctr, nvctr, nseg, keyg, smat)

      call f_release_routine()

    end subroutine sparse_matrix_init_from_data_bigdft


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


    subroutine matrices_set_values(smat, val, mat)
      use module_base
      use sparsematrix_base,only: sparse_matrix, matrices, &
                                  matrices_null, assignment(=), sparsematrix_malloc_ptr, SPARSE_FULL
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      real(kind=8),dimension(:),intent(in) :: val
      type(matrices),intent(inout) :: mat

      call f_routine(id='matrices_set_values')

      if (size(val)/=smat%nvctr) then
          call f_err_throw('The size of the array used to set the matrix contents is wrong: '&
               &//trim(yaml_toa(size(val)))//' instead of '//trim(yaml_toa(smat%nvctr)))
      end if
      if (size(mat%matrix_compr)/=smat%nvctr) then
          call f_err_throw('The size of the matrix array which should be set is wrong: '&
               &//trim(yaml_toa(size(mat%matrix_compr)))//' instead of '//trim(yaml_toa(smat%nvctr)))
      end if
      call f_memcpy(src=val, dest=mat%matrix_compr)

      call f_release_routine()

    end subroutine matrices_set_values


    subroutine matrices_get_values(smat, mat, val)
      use module_base
      use sparsematrix_base,only: sparse_matrix, matrices, &
                                  matrices_null, assignment(=), sparsematrix_malloc_ptr, SPARSE_FULL
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat
      real(kind=8),dimension(:),intent(inout) :: val

      call f_routine(id='matrices_get_values')

      if (size(mat%matrix_compr)/=smat%nvctr) then
          call f_err_throw('The size of the matrix array which should be set is wrong: '&
               &//trim(yaml_toa(size(mat%matrix_compr)))//' instead of '//trim(yaml_toa(smat%nvctr)))
      end if
      if (size(val)/=smat%nvctr) then
          call f_err_throw('The size of the array used to set the matrix contents is wrong: '&
               &//trim(yaml_toa(size(val)))//' instead of '//trim(yaml_toa(smat%nvctr)))
      end if
      call f_memcpy(src=mat%matrix_compr, dest=val)

      call f_release_routine()

    end subroutine matrices_get_values


    function matrices_get_size(smat) result(s)
      use module_base
      use sparsematrix_base,only: sparse_matrix
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      integer :: s

      call f_routine(id='matrices_get_sizes')

      s = smat%nvctr

      call f_release_routine()

    end function matrices_get_size


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
      use sparsematrix_io, only: write_ccs_matrix
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
               smat_in, smat_out, mat_in, mat_out, npl_auto)
      use module_base
      use sparsematrix_base, only: sparse_matrix, matrices
      use ice, only: inverse_chebyshev_expansion_new
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, ndeg, ncalc
      type(sparse_matrix),intent(in) ::smat_in, smat_out
      real(kind=8),dimension(ncalc),intent(in) :: ex
      type(matrices),intent(in) :: mat_in
      type(matrices),dimension(ncalc),intent(inout) :: mat_out
      logical,intent(in),optional :: npl_auto

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

      if (present(npl_auto)) then
          !!call inverse_chebyshev_expansion(iproc, nproc, ndeg, &
          !!     smat_in, smat_out, ncalc, ex, mat_in, mat_out, npl_auto)
          call inverse_chebyshev_expansion_new(iproc, nproc, ndeg, &
               smat_in, smat_out, ncalc, ex, mat_in, mat_out, npl_auto=npl_auto)
      else
          call inverse_chebyshev_expansion_new(iproc, nproc, ndeg, &
               smat_in, smat_out, ncalc, ex, mat_in, mat_out)
      end if

      call f_release_routine()

    end subroutine matrix_chebyshev_expansion


    subroutine matrix_fermi_operator_expansion(iproc, nproc, ndeg, foe_obj, smat_s, smat_h, smat_k, &
               overlap, ham, overlap_minus_one_half, kernel, ebs, &
               calculate_minusonehalf, max_inversion_error, foe_verbosity, label)
      use module_base
      use sparsematrix_base, only: sparse_matrix, matrices
      use foe_base, only: foe_data
      use foe, only: fermi_operator_expansion_new
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      integer,intent(inout) :: ndeg
      type(foe_data),intent(inout) :: foe_obj
      type(sparse_matrix),intent(in) ::smat_s, smat_h, smat_k
      type(matrices),intent(in) :: overlap, ham
      type(matrices),dimension(1),intent(inout) :: overlap_minus_one_half
      type(matrices),intent(inout) :: kernel
      real(kind=8),intent(out) :: ebs
      logical,intent(in),optional :: calculate_minusonehalf
      real(kind=8),intent(in),optional :: max_inversion_error
      integer,intent(in),optional :: foe_verbosity
      character(len=*),intent(in),optional :: label

      ! Local variables
      real(kind=8) :: max_inversion_error_
      logical :: calculate_minusonehalf_
      integer :: foe_verbosity_
      character(len=100) :: label_

      call f_routine(id='matrix_fermi_operator_expansion')

      max_inversion_error_ = 1.d-8
      if (present(max_inversion_error)) max_inversion_error_ = max_inversion_error
      calculate_minusonehalf_ = .true.
      if (present(calculate_minusonehalf)) calculate_minusonehalf_ = calculate_minusonehalf
      foe_verbosity_ = 1
      if (present(foe_verbosity)) foe_verbosity_ = foe_verbosity
      label_ = '001'
      if (present(label)) label_ = trim(label)

      ! Check the dimensions of the internal arrays
      if (size(overlap%matrix_compr)/=smat_s%nvctrp_tg) then
          call f_err_throw('The size of the array overlap%matrix_compr is wrong: '&
               &//trim(yaml_toa(size(overlap%matrix_compr)))//' instead of '//trim(yaml_toa(smat_s%nvctrp_tg)))
      end if
      if (size(ham%matrix_compr)/=smat_h%nvctrp_tg) then
          call f_err_throw('The size of the array mat_out%matrix_compr is wrong: '&
               &//trim(yaml_toa(size(ham%matrix_compr)))//' instead of '//trim(yaml_toa(smat_h%nvctrp_tg)))
      end if
      if (size(overlap_minus_one_half(1)%matrix_compr)/=smat_k%nvctrp_tg) then
          call f_err_throw('The size of the array mat_out%matrix_compr is wrong: '&
               &//trim(yaml_toa(size(overlap_minus_one_half(1)%matrix_compr)))//' instead of '//trim(yaml_toa(smat_k%nvctrp_tg)))
      end if
      if (size(kernel%matrix_compr)/=smat_k%nvctrp_tg) then
          call f_err_throw('The size of the array mat_out%matrix_compr is wrong: '&
               &//trim(yaml_toa(size(kernel%matrix_compr)))//' instead of '//trim(yaml_toa(smat_k%nvctrp_tg)))
      end if

      ! Check that number of non-zero elements of smat_s is not larger than that of smat_h, and that the one of
      ! smat_h is not larger than that of smat_k. This is just a minimal
      ! check to ensure that the sparsity pattern of smat_s is contained within the one of smat_h and the one
      ! of smat_h within smat_k. However this check is not sufficient to make sure that this condition is fulfilled.
      if (smat_h%nvctr>smat_s%nvctr) then
          call f_err_throw('The number of non-zero elements of smat_h ('//&
               trim(yaml_toa(smat_h%nvctr))//') is larger than the one of smat_s ('//&
               trim(yaml_toa(smat_s%nvctr))//')')
      end if
      if (smat_k%nvctr>smat_h%nvctr) then
          call f_err_throw('The number of non-zero elements of smat_k ('//&
               trim(yaml_toa(smat_k%nvctr))//') is larger than the one of smat_h ('//&
               trim(yaml_toa(smat_h%nvctr))//')')
      end if

      call fermi_operator_expansion_new(iproc, nproc, &
           ebs, ndeg, max_inversion_error_, &
           calculate_minusonehalf_, foe_verbosity_, &
           label_, smat_s, smat_h, smat_k, ham, overlap, overlap_minus_one_half, kernel, foe_obj)

      call f_release_routine()

    end subroutine matrix_fermi_operator_expansion

end module sparsematrix_highlevel
