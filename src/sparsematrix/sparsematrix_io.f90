module sparsematrix_io
  implicit none

  private

  public :: write_ccs_matrix
  public :: read_sparse_matrix
  public :: write_sparse_matrix

  contains

    subroutine write_ccs_matrix(filename, nfvctr, nvctr, row_ind, col_ptr, mat_compr)
      use module_base
      implicit none
      !Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(in) :: nfvctr !number of rows/columns
      integer,intent(in) :: nvctr !number of non-zero elements
      integer,dimension(nvctr),intent(in) :: row_ind
      integer,dimension(nfvctr),intent(in) :: col_ptr
      real(kind=8),dimension(nvctr),intent(in) :: mat_compr
      ! Local variables
      integer :: i, iunit

      iunit = 99
      call f_open_file(iunit, file=trim(filename), binary=.false.)

      write(iunit,'(4(i0,1x))') nfvctr, nfvctr, nvctr, 0
      write(iunit,'(100000(i0,1x))') (col_ptr(i),i=1,nfvctr),nvctr+1
      write(iunit,'(100000(i0,1x))') (row_ind(i),i=1,nvctr)
      do i=1,nvctr
          write(iunit,'(es24.15)') mat_compr(i)
      end do

      call f_close(iunit)

    end subroutine write_ccs_matrix


    subroutine read_sparse_matrix(filename, nspin, geocode, cell_dim, nfvctr, nseg, nvctr, keyv, keyg, mat_compr, &
               nat, ntypes, nzatom, nelpsp, atomnames, iatype, rxyz, on_which_atom)
      use module_base
      implicit none
      
      ! Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(out) :: nspin, nfvctr, nseg, nvctr
      character(len=1),intent(out) :: geocode
      real(kind=8),dimension(3),intent(out) :: cell_dim
      integer,dimension(:),pointer,intent(out) :: keyv
      integer,dimension(:,:,:),pointer,intent(out) :: keyg
      real(kind=8),dimension(:),pointer,intent(out) :: mat_compr
      integer,intent(out),optional :: nat, ntypes
      integer,dimension(:),pointer,intent(inout),optional :: nzatom, nelpsp, iatype
      character(len=20),dimension(:),pointer,intent(inout),optional :: atomnames
      real(kind=8),dimension(:,:),pointer,intent(inout),optional :: rxyz
      integer,dimension(:),pointer,intent(inout),optional :: on_which_atom

      ! Local variables
      integer :: iunit, dummy_int, iseg, icol, irow, jorb, ind, ispin, iat, ntypes_, nat_, itype
      real(kind=8) :: dummy_double
      character(len=20) :: dummy_char
      logical :: read_rxyz, read_on_which_atom

      call f_routine(id='read_sparse_matrix')

      if (present(nat) .and. present(ntypes) .and. present(nzatom) .and.  &
          present(nelpsp) .and. present(atomnames) .and. present(iatype) .and. present(rxyz)) then
          read_rxyz = .true.
      else if (present(nat) .or. present(ntypes) .or. present(nzatom) .or.  &
          present(nelpsp) .or. present(atomnames) .or. present(iatype) .or. present(rxyz)) then
          call f_err_throw("not all optional arguments were given", &
               err_name='BIGDFT_RUNTIME_ERROR')
      else
          read_rxyz = .false.
      end if
      
      if (present(on_which_atom)) then
          read_on_which_atom = .true.
      else
          read_on_which_atom = .false.
      end if

      iunit = 99
      call f_open_file(iunit, file=trim(filename), binary=.false.)

      if (read_rxyz) then
          read(iunit,*) nat, ntypes, nspin, geocode, cell_dim
          nzatom = f_malloc_ptr(ntypes,id='nzatom')
          nelpsp = f_malloc_ptr(ntypes,id='nelpsp')
          atomnames = f_malloc0_str_ptr(len(atomnames),ntypes,id='atomnames')

          do itype=1,ntypes
              read(iunit,*) nzatom(itype), nelpsp(itype), atomnames(itype)
          end do
          rxyz = f_malloc_ptr((/3,nat/),id='rxyz')
          iatype = f_malloc_ptr(nat,id='iatype')
          do iat=1,nat
              read(iunit,*) iatype(iat), rxyz(1,iat), rxyz(2,iat), rxyz(3,iat)
          end do
      else
          read(iunit,*) nat_, ntypes_, nspin, geocode
          do itype=1,ntypes_
              read(iunit,*) dummy_int, dummy_int, dummy_char
          end do
          do iat=1,nat_
              read(iunit,*) dummy_int, dummy_double, dummy_double, dummy_double
          end do
      end if
      read(iunit,*) nfvctr, nseg, nvctr
      keyv = f_malloc_ptr(nseg,id='keyv')
      keyg = f_malloc_ptr((/2,2,nseg/),id='keyg')
      do iseg=1,nseg
          read(iunit,*) keyv(iseg), keyg(1,1,iseg), keyg(2,1,iseg), keyg(1,2,iseg), keyg(2,2,iseg)
      end do
      mat_compr = f_malloc_ptr(nvctr*nspin,id='mat_compr')
      if (read_on_which_atom) then
          nullify(on_which_atom)
          on_which_atom = f_malloc_ptr(nfvctr,id='on_which_atom')
          ind = 0
          do ispin=1,nspin
              do iseg=1,nseg
                  icol = keyg(1,2,iseg)
                  do jorb=keyg(1,1,iseg),keyg(2,1,iseg)
                      irow = jorb
                      ind = ind + 1
                      read(iunit,*) mat_compr(ind), on_which_atom(irow), on_which_atom(icol)
                  end do
              end do
          end do
      else
          ind = 0
          do ispin=1,nspin
              do iseg=1,nseg
                  icol = keyg(1,2,iseg)
                  do jorb=keyg(1,1,iseg),keyg(2,1,iseg)
                      irow = jorb
                      ind = ind + 1
                      read(iunit,*) mat_compr(ind), dummy_int, dummy_int
                  end do
              end do
          end do
      end if

      call f_close(iunit)

      call f_release_routine()

    end subroutine read_sparse_matrix


    !> Write a sparse matrix to disk.
    !! ATTENTION: This routine must be called by all MPI tasks due to the fact that the matrix 
    !! in distributed among the matrix taksgroups
    subroutine write_sparse_matrix(nat, ntypes, iatype, rxyz, nzatom, nelpsp, atomnames, smat, mat, filename)
      use module_base
      use sparsematrix_base, only: sparse_matrix, matrices, SPARSE_FULL, &
                                   assignment(=), sparsematrix_malloc
      use sparsematrix, only: gather_matrix_from_taskgroups
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: nat, ntypes
      integer,dimension(nat),intent(in) :: iatype
      real(kind=8),dimension(3,nat),intent(in) :: rxyz
      integer,dimension(ntypes),intent(in) :: nzatom, nelpsp
      character(len=*),dimension(ntypes),intent(in) :: atomnames
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat
      character(len=*),intent(in) :: filename

      ! Local variables
      integer :: iunit, iseg, icol, irow, jorb, iat, jat, ind, ispin, itype
      real(kind=8),dimension(:),allocatable :: matrix_compr

      call f_routine(id='write_sparse_matrix')

      matrix_compr = sparsematrix_malloc(smat,iaction=SPARSE_FULL,id='matrix_compr')
      call gather_matrix_from_taskgroups(bigdft_mpi%iproc, bigdft_mpi%nproc, &
           smat, mat%matrix_compr, matrix_compr)

      if (bigdft_mpi%iproc==0) then

          iunit = 99
          call f_open_file(iunit, file=trim(filename), binary=.false.)

          write(iunit,'(i10,2i6,3x,a,3es24.16,a)') nat, ntypes, smat%nspin, smat%geocode, smat%cell_dim, &
              '   # number of atoms, number of atom types, nspin, geocode, cell_dim'
          do itype=1,ntypes
              write(iunit,'(2i8,3x,a,a)') nzatom(itype), nelpsp(itype), trim(atomnames(itype)), &
                  '   # nz, nelpsp, name'
          end do
          do iat=1,nat
              write(iunit,'(i5, 3es24.16,a,i0)') iatype(iat), rxyz(1:3,iat), '   # atom no. ',iat
          end do
          write(iunit,'(3i12,a)') smat%nfvctr, smat%nseg, smat%nvctr, '   # nfvctr, nseg, nvctr'
          do iseg=1,smat%nseg
              write(iunit,'(5i12,a)') smat%keyv(iseg), smat%keyg(1,1,iseg), smat%keyg(2,1,iseg), &
                  smat%keyg(1,2,iseg), smat%keyg(2,2,iseg), '   # keyv, keyg(1,1), keyg(2,1), keyg(1,2), keyg(2,2)'
          end do
          ind = 0
          do ispin=1,smat%nspin
              do iseg=1,smat%nseg
                  icol = smat%keyg(1,2,iseg)
                  iat = smat%on_which_atom(icol)
                  do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                      irow = jorb
                      jat = smat%on_which_atom(irow)
                      ind = ind + 1
                      write(iunit,'(es24.16,2i12,a)') matrix_compr(ind), jat, iat, '   # matrix, jat, iat'
                  end do
              end do
          end do

          call f_close(iunit)

          call f_free(matrix_compr)

      end if

      call f_release_routine()

    end subroutine write_sparse_matrix

end module sparsematrix_io
