module sparsematrix_highlevel
  use sparsematrix_base
  private

  !> Initialization routines
  public :: sparse_matrix_and_matrices_init_from_file_ccs
  public :: sparse_matrix_init_from_file_ccs
  public :: sparse_matrix_init_from_data_ccs
  public :: sparse_matrix_and_matrices_init_from_file_bigdft
  public :: sparse_matrix_init_from_file_bigdft
  public :: sparse_matrix_metadata_init_from_file
  public :: matrices_init_from_file_bigdft
  public :: read_sparse_matrix_metadata
  public :: sparse_matrix_init_from_data_bigdft
  public :: matrices_init
  !> Data manipulation routines
  public :: matrices_set_values
  public :: matrices_get_values
  public :: matrices_get_size
  public :: ccs_data_from_sparse_matrix
  public :: ccs_matrix_write
  !> Wrapper for complete operations
  public :: matrix_matrix_multiplication
  public :: matrix_chebyshev_expansion
  public :: matrix_fermi_operator_expansion
  public :: get_selected_eigenvalues_from_FOE
  public :: trace_AB
  public :: trace_A

  contains

    subroutine sparse_matrix_and_matrices_init_from_file_ccs(filename, iproc, nproc, comm, smat, mat, &
               init_matmul)
      use sparsematrix_init, only: read_ccs_format
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      character(len=*),intent(in) :: filename
      type(sparse_matrix),intent(out) :: smat
      type(matrices),intent(out) :: mat
      logical,intent(in),optional :: init_matmul

      ! Local variables
      integer :: nfvctr, nvctr
      integer,dimension(:),pointer :: col_ptr, row_ind
      real(kind=mp),dimension(:),pointer :: val
      logical :: init_matmul_

      call f_routine(id='sparse_matrix_and_matrices_init_from_file_ccs')

      ! Read in the matrix
      call read_ccs_format(filename, nfvctr, nvctr, col_ptr, row_ind, val)

      if (present(init_matmul)) then
          init_matmul_ = init_matmul
      else
          init_matmul_ = .false.
      end if
    
      ! Generate the sparse_matrix type
      call sparse_matrix_init_from_data_ccs(iproc, nproc, comm, nfvctr, nvctr, row_ind, col_ptr, smat, init_matmul_)

      ! Generate the matrices type
      call matrices_init_from_data(smat, val, mat)

      ! Deallocate the pointers
      call f_free_ptr(col_ptr)
      call f_free_ptr(row_ind)
      call f_free_ptr(val)

      call f_release_routine()

    end subroutine sparse_matrix_and_matrices_init_from_file_ccs


    subroutine sparse_matrix_init_from_file_ccs(filename, iproc, nproc, comm, smat)
      use sparsematrix_init, only: read_ccs_format
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      character(len=*),intent(in) :: filename
      type(sparse_matrix),intent(out) :: smat

      ! Local variables
      integer :: nfvctr, nvctr
      integer,dimension(:),pointer :: col_ptr, row_ind
      real(kind=mp),dimension(:),pointer :: val

      call f_routine(id='sparse_matrix_and_matrices_init_from_file_ccs')

      ! Read in the matrix
      call read_ccs_format(filename, nfvctr, nvctr, col_ptr, row_ind, val)

      ! Generate the sparse_matrix type
      call sparse_matrix_init_from_data_ccs(iproc, nproc, comm, nfvctr, nvctr, row_ind, col_ptr, smat)

      ! Deallocate the pointers
      call f_free_ptr(col_ptr)
      call f_free_ptr(row_ind)
      call f_free_ptr(val)

      call f_release_routine()

    end subroutine sparse_matrix_init_from_file_ccs


    subroutine sparse_matrix_init_from_data_ccs(iproc, nproc, comm, nfvctr, nvctr, row_ind, col_ptr, smat, &
               init_matmul, nvctr_mult, row_ind_mult, col_ptr_mult)
      use sparsematrix_init, only: ccs_to_sparsebigdft_short, &
                                   bigdft_to_sparsebigdft, init_matrix_taskgroups
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, nfvctr, nvctr
      integer,dimension(nvctr),intent(in) :: row_ind
      integer,dimension(nfvctr),intent(in) :: col_ptr
      type(sparse_matrix),intent(out) :: smat
      logical,intent(in),optional :: init_matmul
      integer,intent(in),optional :: nvctr_mult
      integer,dimension(:),intent(in),optional :: row_ind_mult
      integer,dimension(:),intent(in),optional :: col_ptr_mult

      ! Local variables
      integer :: nseg, nseg_mult
      integer,dimension(:),pointer :: keyv, keyv_mult
      integer,dimension(:,:,:),pointer :: keyg, keyg_mult
      logical :: init_matmul_

      call f_routine(id='sparse_matrix_init_from_data_ccs')

      if (present(init_matmul)) then
          init_matmul_ = init_matmul
      else
          init_matmul_ = .false.
      end if

      if (init_matmul_) then
          if (.not.present(nvctr_mult)) then
              call f_err_throw("'nvctr_mult' not present",err_name='SPARSEMATRIX_INITIALIZATION_ERROR')
          end if
          if (.not.present(row_ind_mult)) then
              call f_err_throw("'row_ind_mult' not present",err_name='SPARSEMATRIX_INITIALIZATION_ERROR')
          end if
          if (.not.present(col_ptr_mult)) then
              call f_err_throw("'col_ptr_mult' not present",err_name='SPARSEMATRIX_INITIALIZATION_ERROR')
          end if
          if (size(row_ind_mult)/=nvctr_mult) then
              call f_err_throw("'col_ptr_mult' has wrong size",err_name='SPARSEMATRIX_INITIALIZATION_ERROR')
          end if
          if (size(col_ptr_mult)/=nfvctr) then
              call f_err_throw("'col_ptr_mult' has wrong size",err_name='SPARSEMATRIX_INITIALIZATION_ERROR')
          end if
      end if
    
      ! Convert the sparsity pattern to the BigDFT format
      call ccs_to_sparsebigdft_short(nfvctr, nvctr, row_ind, col_ptr, nseg, keyv, keyg)

      if (init_matmul_) then
          call ccs_to_sparsebigdft_short(nfvctr, nvctr_mult, row_ind_mult, col_ptr_mult, nseg_mult, keyv_mult, keyg_mult)
      end if

      ! Create the sparse_matrix structure
      if (init_matmul_) then
          call bigdft_to_sparsebigdft(iproc, nproc, comm, nfvctr, nvctr, nseg, keyg, smat, &
               init_matmul=init_matmul_, nseg_mult=nseg_mult, nvctr_mult=nvctr_mult, keyg_mult=keyg_mult)
      else
          call bigdft_to_sparsebigdft(iproc, nproc, comm, nfvctr, nvctr, nseg, keyg, smat, &
               init_matmul=init_matmul_)
      end if

      ! Deallocate the pointers
      call f_free_ptr(keyv)
      call f_free_ptr(keyg)
      if (init_matmul_) then
          call f_free_ptr(keyv_mult)
          call f_free_ptr(keyg_mult)
      end if

      call f_release_routine()

    end subroutine sparse_matrix_init_from_data_ccs


    subroutine sparse_matrix_and_matrices_init_from_file_bigdft(filename, iproc, nproc, comm, smat, mat, &
               init_matmul, filename_mult)
      use sparsematrix_init, only: bigdft_to_sparsebigdft
      use sparsematrix_io, only: read_sparse_matrix
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      character(len=*),intent(in) :: filename
      type(sparse_matrix),intent(out) :: smat
      type(matrices),intent(out) :: mat
      logical,intent(in),optional :: init_matmul
      character(len=*),intent(in),optional :: filename_mult
      ! Optional variables that are contained within the sparse matrix format
      !!integer,intent(out),optional :: nat, ntypes
      !!integer,dimension(:),pointer,intent(inout),optional :: nzatom, nelpsp, iatype
      !!character(len=20),dimension(:),pointer,intent(inout),optional :: atomnames
      !!real(kind=mp),dimension(:,:),pointer,intent(inout),optional :: rxyz
      !!integer,dimension(:),pointer,intent(inout),optional :: on_which_atom

      ! Local variables
      integer :: nspin, nfvctr, nseg, nvctr
!      character(len=1) :: geocode
      integer,dimension(:),pointer :: keyv
      integer,dimension(:,:,:),pointer :: keyg
      real(kind=mp),dimension(:),pointer :: val
      logical :: init_matmul_
!      integer,dimension(:),pointer :: on_which_atom_
!      real(kind=8),dimension(3) :: cell_dim
      real(kind=mp),dimension(:,:),pointer :: rxyz_
      real(kind=mp),dimension(3) :: cell_dim
      type(sparse_matrix) :: smat_mult

      call f_routine(id='sparse_matrix_and_matrices_init_from_file_bigdft')

      ! Read in the matrix
      call read_sparse_matrix(filename, nspin, nfvctr, nseg, nvctr, keyv, keyg, val)

      if (present(init_matmul)) then
          init_matmul_ = init_matmul
      else
          init_matmul_ = .false.
      end if

      if (init_matmul_) then
          if (.not.present(filename_mult)) then
              call f_err_throw("'filename_mult' not present",err_name='SPARSEMATRIX_INITIALIZATION_ERROR')
          end if
          !call sparse_matrix_init_from_file_bigdft(filename_mult, iproc, nproc, comm, smat_mult, init_matmul=.false.)
      end if


      !!! Create the sparse_matrix structure
      !!call bigdft_to_sparsebigdft(iproc, nproc, comm, nfvctr, nvctr, nseg, keyg, smat, &
      !!     init_matmul=init_matmul_)!, nspin=nspin, geocode=geocode, cell_dim=cell_dim, on_which_atom=on_which_atom_)

      if (init_matmul_) then
          call sparse_matrix_init_from_file_bigdft(filename, iproc, nproc, comm, smat, init_matmul_, filename_mult)
      else
          call sparse_matrix_init_from_file_bigdft(filename, iproc, nproc, comm, smat, init_matmul_)
      end if

      ! Generate the matrices type
      call matrices_init_from_data(smat, val, mat)

      ! Deallocate the pointers
      call f_free_ptr(keyv)
      call f_free_ptr(keyg)
      call f_free_ptr(val)


      call f_release_routine()

    end subroutine sparse_matrix_and_matrices_init_from_file_bigdft


    recursive subroutine sparse_matrix_init_from_file_bigdft(filename, iproc, nproc, comm, smat, init_matmul, filename_mult)
      use sparsematrix_init, only: bigdft_to_sparsebigdft
      use sparsematrix_io, only: read_sparse_matrix
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      character(len=*),intent(in) :: filename
      type(sparse_matrix),intent(out) :: smat
      logical,intent(in),optional :: init_matmul
      character(len=*),intent(in),optional :: filename_mult

      ! Local variables
      integer :: nspin, nfvctr, nseg, nvctr
      character(len=1) :: geocode
      integer,dimension(:),pointer :: keyv
      integer,dimension(:,:,:),pointer :: keyg
      real(kind=mp),dimension(:),pointer :: val
      real(kind=mp),dimension(3) :: cell_dim
      integer,dimension(:),pointer :: on_which_atom
      logical :: init_matmul_
      type(sparse_matrix) :: smat_mult

      call f_routine(id='sparse_matrix_and_matrices_init_from_file_bigdft')

      if (present(init_matmul)) then
          init_matmul_ = init_matmul
      else
          init_matmul_ = .false.
      end if

      ! Read in the matrix
      call read_sparse_matrix(filename, nspin, nfvctr, nseg, nvctr, keyv, keyg, val)

      if (init_matmul_) then
          if (.not.present(filename_mult)) then
              call f_err_throw("'filename_mult' not present",err_name='SPARSEMATRIX_INITIALIZATION_ERROR')
          end if
          call sparse_matrix_init_from_file_bigdft(filename_mult, iproc, nproc, comm, smat_mult, init_matmul=.false.)
      end if

      ! Create the sparse_matrix structure
      !!call bigdft_to_sparsebigdft(iproc, nproc, comm, nfvctr, nvctr, nseg, keyg, smat, &
      !!     nspin=nspin, geocode=geocode, cell_dim=cell_dim, on_which_atom=on_which_atom)
      if (init_matmul_) then
          call sparse_matrix_init_from_data_bigdft(iproc, nproc, comm, nfvctr, nvctr, nseg, keyg, smat, init_matmul_, &
               nseg_mult=smat_mult%nseg, nvctr_mult=smat_mult%nvctr, keyg_mult=smat_mult%keyg)
          call deallocate_sparse_matrix(smat_mult)
      else
          call sparse_matrix_init_from_data_bigdft(iproc, nproc, comm, nfvctr, nvctr, nseg, keyg, smat, init_matmul_)
      end if

      ! Deallocate the pointers
      call f_free_ptr(keyv)
      call f_free_ptr(keyg)
      call f_free_ptr(val)

      call f_release_routine()

    end subroutine sparse_matrix_init_from_file_bigdft


    subroutine matrices_init_from_file_bigdft(filename, iproc, nproc, comm, smat, mat)
      use sparsematrix_io, only: read_sparse_matrix
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      character(len=*),intent(in) :: filename
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(out) :: mat
    
      ! Local variables
      integer :: nspin, nfvctr, nseg, nvctr, iseg
      integer,dimension(:),pointer :: keyv
      integer,dimension(:,:,:),pointer :: keyg
      real(kind=mp),dimension(:),pointer :: val
    
      call f_routine(id='matrices_init_from_file_bigdft')
    
      ! Read in the matrix
      call read_sparse_matrix(filename, nspin, nfvctr, nseg, nvctr, keyv, keyg, val)

      ! Check that it is consistent with the provided sparse matrix type
      if (nspin/=smat%nspin) then
          call f_err_throw('nspin/=smat%nspin')
      end if
      if (nfvctr/=smat%nfvctr) then
          call f_err_throw('nfvctr/=smat%nfvctr')
      end if
      if (nseg/=smat%nseg) then
          call f_err_throw('nseg/=smat%nseg')
      end if
      if (nvctr/=smat%nvctr) then
          call f_err_throw('nvctr/=smat%nvctr')
      end if
      do iseg=1,nseg
          if (keyv(iseg)/=smat%keyv(iseg)) then
              call f_err_throw('keyv(iseg)/=smat%keyv(iseg)')
          end if
          if (keyg(1,1,iseg)/=smat%keyg(1,1,iseg)) then
              call f_err_throw('keyg(1,1,iseg)/=smat%keyg(1,1,iseg)')
          end if
          if (keyg(2,1,iseg)/=smat%keyg(2,1,iseg)) then
              call f_err_throw('keyg(2,1,iseg)/=smat%keyg(2,1,iseg)')
          end if
          if (keyg(1,2,iseg)/=smat%keyg(1,2,iseg)) then
              call f_err_throw('keyg(1,2,iseg)/=smat%keyg(1,2,iseg)')
          end if
          if (keyg(2,2,iseg)/=smat%keyg(2,2,iseg)) then
              call f_err_throw('keyg(2,2,iseg)/=smat%keyg(2,2,iseg)')
          end if
      end do

    
      ! Generate the matrices type
      call matrices_init_from_data(smat, val, mat)
    
      ! Deallocate the pointers
      call f_free_ptr(keyv)
      call f_free_ptr(keyg)
      call f_free_ptr(val)

    
      call f_release_routine()
    
    end subroutine matrices_init_from_file_bigdft


    subroutine sparse_matrix_metadata_init_from_file(filename, smmd)
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


    subroutine sparse_matrix_init_from_data_bigdft(iproc, nproc, comm, nfvctr, nvctr, nseg, keyg, smat, &
               init_matmul, nseg_mult, nvctr_mult, keyg_mult)
      use sparsematrix_init, only: ccs_to_sparsebigdft_short, &
                                   bigdft_to_sparsebigdft, init_matrix_taskgroups
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, nfvctr, nvctr, nseg
      integer,dimension(2,2,nseg),intent(in) :: keyg
      type(sparse_matrix),intent(out) :: smat
      logical,intent(in) :: init_matmul
      integer,intent(in),optional :: nseg_mult, nvctr_mult
      integer,dimension(:,:,:),intent(in),optional :: keyg_mult

      call f_routine(id='sparse_matrix_init_from_data_bigdft')

      if (init_matmul) then
          if (.not.present(nseg_mult)) then
              call f_err_throw("'nseg_mult' not present",err_name='SPARSEMATRIX_INITIALIZATION_ERROR')
          end if
          if (.not.present(nvctr_mult)) then
              call f_err_throw("'nvctr_mult' not present",err_name='SPARSEMATRIX_INITIALIZATION_ERROR')
          end if
          if (.not.present(nseg_mult)) then
              call f_err_throw("'nvctr_mult' not present",err_name='SPARSEMATRIX_INITIALIZATION_ERROR')
          end if
      end if

      ! Create the sparse_matrix structure
      if (init_matmul) then
          call bigdft_to_sparsebigdft(iproc, nproc, comm, nfvctr, nvctr, nseg, keyg, smat, &
               init_matmul=init_matmul, nseg_mult=nseg_mult, nvctr_mult=nvctr_mult, keyg_mult=keyg_mult)
      else
          call bigdft_to_sparsebigdft(iproc, nproc, comm, nfvctr, nvctr, nseg, keyg, smat)
      end if

      call f_release_routine()

    end subroutine sparse_matrix_init_from_data_bigdft


    subroutine matrices_init_from_data(smat, val, mat)
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      real(kind=mp),dimension(smat%nvctr),intent(in) :: val
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
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(out) :: mat

      mat = matrices_null()
      mat%matrix_compr = sparsematrix_malloc_ptr(smat, iaction=SPARSE_FULL, id='mat%matrix_compr')

    end subroutine matrices_init


    subroutine matrices_set_values(smat, val, mat)
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      real(kind=mp),dimension(:),intent(in) :: val
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
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat
      real(kind=mp),dimension(:),intent(inout) :: val

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
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      integer :: s

      call f_routine(id='matrices_get_sizes')

      s = smat%nvctr

      call f_release_routine()

    end function matrices_get_size


    subroutine ccs_data_from_sparse_matrix(smat, row_ind, col_ptr)
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


    subroutine matrix_chebyshev_expansion(iproc, nproc, comm, ncalc, ex, &
               smat_in, smat_out, mat_in, mat_out, npl_auto, ice_obj)
      use ice, only: inverse_chebyshev_expansion_new
      use foe_base, only: foe_data
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, ncalc
      type(sparse_matrix),intent(in) ::smat_in, smat_out
      real(kind=mp),dimension(ncalc),intent(in) :: ex
      type(matrices),intent(in) :: mat_in
      type(matrices),dimension(ncalc),intent(inout) :: mat_out
      logical,intent(in),optional :: npl_auto
      type(foe_data),intent(inout),optional,target :: ice_obj

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
          if (present(ice_obj)) then
              call inverse_chebyshev_expansion_new(iproc, nproc, comm, &
                   smat_in, smat_out, ncalc, ex, mat_in, mat_out, npl_auto=npl_auto, ice_objx=ice_obj)
          else
              call inverse_chebyshev_expansion_new(iproc, nproc, comm, &
                   smat_in, smat_out, ncalc, ex, mat_in, mat_out, npl_auto=npl_auto)
          end if
      else
          if (present(ice_obj)) then
              call inverse_chebyshev_expansion_new(iproc, nproc, comm, &
                   smat_in, smat_out, ncalc, ex, mat_in, mat_out, ice_objx=ice_obj)
          else
              call inverse_chebyshev_expansion_new(iproc, nproc, comm, &
                   smat_in, smat_out, ncalc, ex, mat_in, mat_out)
          end if
      end if

      call f_release_routine()

    end subroutine matrix_chebyshev_expansion


    subroutine matrix_fermi_operator_expansion(iproc, nproc, comm, foe_obj, ice_obj, smat_s, smat_h, smat_k, &
               overlap, ham, overlap_minus_one_half, kernel, ebs, &
               calculate_minusonehalf, foe_verbosity, symmetrize_kernel, calculate_energy_density_kernel, energy_kernel)
      use foe_base, only: foe_data
      use foe, only: fermi_operator_expansion_new
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      type(foe_data),intent(inout) :: foe_obj, ice_obj
      type(sparse_matrix),intent(in) ::smat_s, smat_h, smat_k
      type(matrices),intent(in) :: overlap, ham
      type(matrices),dimension(1),intent(inout) :: overlap_minus_one_half
      type(matrices),intent(inout) :: kernel
      real(kind=mp),intent(out) :: ebs
      logical,intent(in),optional :: calculate_minusonehalf, symmetrize_kernel, calculate_energy_density_kernel
      integer,intent(in),optional :: foe_verbosity
      type(matrices),intent(inout),optional :: energy_kernel

      ! Local variables
      logical :: calculate_minusonehalf_, symmetrize_kernel_, calculate_energy_density_kernel_
      integer :: foe_verbosity_

      call f_routine(id='matrix_fermi_operator_expansion')

      calculate_minusonehalf_ = .true.
      if (present(calculate_minusonehalf)) calculate_minusonehalf_ = calculate_minusonehalf
      foe_verbosity_ = 1
      if (present(foe_verbosity)) foe_verbosity_ = foe_verbosity
      symmetrize_kernel_ = .false.
      if (present(symmetrize_kernel)) symmetrize_kernel_ = symmetrize_kernel
      calculate_energy_density_kernel_ = .false.
      if (present(calculate_energy_density_kernel)) calculate_energy_density_kernel_ = calculate_energy_density_kernel

      ! Check the optional arguments
      if (calculate_energy_density_kernel_) then
          if (.not.present(energy_kernel)) then
              call f_err_throw('energy_kernel not present',err_name='SPARSEMATRIX_RUNTIME_ERROR')
          end if
      end if

      ! Check the dimensions of the internal arrays
      if (size(overlap%matrix_compr)/=smat_s%nvctrp_tg*smat_s%nspin) then
          call f_err_throw('The size of the array overlap%matrix_compr is wrong: '&
               &//trim(yaml_toa(size(overlap%matrix_compr)))//' instead of '//&
               &trim(yaml_toa(smat_s%nvctrp_tg*smat_s%nspin)))
      end if
      if (size(ham%matrix_compr)/=smat_h%nvctrp_tg*smat_h%nspin) then
          call f_err_throw('The size of the array mat_out%matrix_compr is wrong: '&
               &//trim(yaml_toa(size(ham%matrix_compr)))//' instead of '//&
               &trim(yaml_toa(smat_h%nvctrp_tg*smat_h%nspin)))
      end if
      if (size(overlap_minus_one_half(1)%matrix_compr)/=smat_k%nvctrp_tg*smat_k%nspin) then
          call f_err_throw('The size of the array mat_out%matrix_compr is wrong: '&
               &//trim(yaml_toa(size(overlap_minus_one_half(1)%matrix_compr)))//' instead of '//&
               &trim(yaml_toa(smat_k%nvctrp_tg*smat_k%nspin)))
      end if
      if (size(kernel%matrix_compr)/=smat_k%nvctrp_tg*smat_k%nspin) then
          call f_err_throw('The size of the array mat_out%matrix_compr is wrong: '&
               &//trim(yaml_toa(size(kernel%matrix_compr)))//' instead of '//&
               &trim(yaml_toa(smat_k%nvctrp_tg*smat_k%nspin)))
      end if

      ! Check that number of non-zero elements of smat_s is not larger than that of smat_h, and that the one of
      ! smat_h is not larger than that of smat_k. This is just a minimal
      ! check to ensure that the sparsity pattern of smat_s is contained within the one of smat_h and the one
      ! of smat_h within smat_k. However this check is not sufficient to make sure that this condition is fulfilled.
      if (smat_s%nvctr>smat_h%nvctr) then
          call f_err_throw('The number of non-zero elements of smat_s ('//&
               trim(yaml_toa(smat_s%nvctr))//') is larger than the one of smat_h ('//&
               trim(yaml_toa(smat_h%nvctr))//')')
      end if
      if (smat_h%nvctr>smat_k%nvctr) then
          call f_err_throw('The number of non-zero elements of smat_h ('//&
               trim(yaml_toa(smat_h%nvctr))//') is larger than the one of smat_k ('//&
               trim(yaml_toa(smat_k%nvctr))//')')
      end if

      if (calculate_energy_density_kernel_) then
          call fermi_operator_expansion_new(iproc, nproc, comm, &
               ebs, &
               calculate_minusonehalf_, foe_verbosity_, &
               smat_s, smat_h, smat_k, ham, overlap, overlap_minus_one_half, kernel, foe_obj, ice_obj, &
               symmetrize_kernel_, calculate_energy_density_kernel_, energy_kernel_=energy_kernel)
      else
          call fermi_operator_expansion_new(iproc, nproc, comm, &
               ebs, &
               calculate_minusonehalf_, foe_verbosity_, &
               smat_s, smat_h, smat_k, ham, overlap, overlap_minus_one_half, kernel, foe_obj, ice_obj, &
               symmetrize_kernel_, calculate_energy_density_kernel_)
      end if

      call f_release_routine()

    end subroutine matrix_fermi_operator_expansion



    subroutine get_selected_eigenvalues_from_FOE(iproc, nproc, comm, iev_min, iev_max, &
               smat_s, smat_h, smat_k, overlap, ham, overlap_minus_one_half, evals, &
               fscale, calculate_minusonehalf, foe_verbosity)
      use foe_base, only: foe_data
      use foe, only: get_selected_eigenvalues
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, iev_min, iev_max
      type(sparse_matrix),intent(in) ::smat_s, smat_h, smat_k
      type(matrices),intent(in) :: overlap, ham
      type(matrices),dimension(1),intent(inout) :: overlap_minus_one_half
      real(kind=mp),dimension(iev_min:iev_max),intent(out) :: evals
      real(mp),intent(in),optional :: fscale
      logical,intent(in),optional :: calculate_minusonehalf
      integer,intent(in),optional :: foe_verbosity

      ! Local variables
      logical :: calculate_minusonehalf_, symmetrize_kernel_
      integer :: foe_verbosity_
      real(mp) :: fscale_

      call f_routine(id='matrix_fermi_operator_expansion')

      calculate_minusonehalf_ = .true.
      if (present(calculate_minusonehalf)) calculate_minusonehalf_ = calculate_minusonehalf
      foe_verbosity_ = 1
      if (present(foe_verbosity)) foe_verbosity_ = foe_verbosity
      symmetrize_kernel_ = .false.
      fscale_ = 5.e-3_mp
      if (present(fscale)) fscale_ = fscale

      ! Check the dimensions of the internal arrays
      if (size(overlap%matrix_compr)/=smat_s%nvctrp_tg*smat_s%nspin) then
          call f_err_throw('The size of the array overlap%matrix_compr is wrong: '&
               &//trim(yaml_toa(size(overlap%matrix_compr)))//' instead of '//&
               &trim(yaml_toa(smat_s%nvctrp_tg*smat_s%nspin)))
      end if
      if (size(ham%matrix_compr)/=smat_h%nvctrp_tg*smat_h%nspin) then
          call f_err_throw('The size of the array mat_out%matrix_compr is wrong: '&
               &//trim(yaml_toa(size(ham%matrix_compr)))//' instead of '//&
               &trim(yaml_toa(smat_h%nvctrp_tg*smat_h%nspin)))
      end if
      if (size(overlap_minus_one_half(1)%matrix_compr)/=smat_k%nvctrp_tg*smat_k%nspin) then
          call f_err_throw('The size of the array mat_out%matrix_compr is wrong: '&
               &//trim(yaml_toa(size(overlap_minus_one_half(1)%matrix_compr)))//' instead of '//&
               &trim(yaml_toa(smat_k%nvctrp_tg*smat_k%nspin)))
      end if

      ! Check that number of non-zero elements of smat_s is not larger than that of smat_h, and that the one of 
      ! smat_h is not larger than that of smat_k. This is just a minimal
      ! check to ensure that the sparsity pattern of smat_s is contained within the one of smat_h and the one 
      ! of smat_h within smat_k. However this check is not sufficient to make sure that this condition is fulfilled.
      if (smat_s%nvctr>smat_h%nvctr) then
          call f_err_throw('The number of non-zero elements of smat_s ('//&
               trim(yaml_toa(smat_s%nvctr))//') is larger than the one of smat_h ('//&
               trim(yaml_toa(smat_h%nvctr))//')')
      end if
      if (smat_h%nvctr>smat_k%nvctr) then
          call f_err_throw('The number of non-zero elements of smat_h ('//&
               trim(yaml_toa(smat_h%nvctr))//') is larger than the one of smat_k ('//&
               trim(yaml_toa(smat_k%nvctr))//')')
      end if

      ! Check whether the eigenvalues to be calculated are within the possible range
      if (iev_min<1 .or. iev_max>smat_s%nfvctr) then
          if (iproc==0) then
              call f_err_throw('The required eigenvalues are outside of the possible range')
          end if
      end if


      call get_selected_eigenvalues(iproc, nproc, comm, calculate_minusonehalf_, foe_verbosity_, &
           iev_min, iev_max, fscale, &
           smat_s, smat_h, smat_k, ham, overlap, overlap_minus_one_half, evals)


      call f_release_routine()

    end subroutine get_selected_eigenvalues_from_FOE


    !> Calculates the trace of the spin component ispin of the matrix product amat*bmat.
    !! WARNING: It is mandatory that the sparsity pattern of amat be contained
    !! within the sparsity pattern of bmat!
    function trace_AB(iproc, nproc, comm, asmat, bsmat, amat, bmat, ispin)
      use sparsematrix, only: trace_sparse_matrix_product
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, ispin
      type(sparse_matrix),intent(in) :: asmat, bsmat
      type(matrices),intent(in) :: amat
      type(matrices),intent(in) :: bmat
      real(kind=mp) :: trace_AB

      ! Local variables
      integer :: iashift, ibshift

      iashift=(ispin-1)*asmat%nvctrp_tg
      ibshift=(ispin-1)*bsmat%nvctrp_tg

      trace_AB = trace_sparse_matrix_product(iproc, nproc, comm, asmat, bsmat, &
                 amat%matrix_compr(iashift+1:), &
                 bmat%matrix_compr(ibshift+1:))

    end function trace_AB


    !> Calculates the trace of the spin component ispin of the matrix amat
    function trace_A(iproc, nproc, comm, asmat, amat, ispin)
      use sparsematrix, only: trace_sparse_matrix
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, ispin
      type(sparse_matrix),intent(in) :: asmat
      type(matrices),intent(in) :: amat
      real(kind=mp) :: trace_A

      ! Local variables
      integer :: iashift

      iashift=(ispin-1)*asmat%nvctrp_tg

      trace_A = trace_sparse_matrix(iproc, nproc, comm, asmat, amat%matrix_compr(iashift+1:))

    end function trace_A

end module sparsematrix_highlevel
