module sparsematrix_io
  implicit none

  private

  public :: write_ccs_matrix
  public :: read_sparse_matrix

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


    subroutine read_sparse_matrix(filename, nspin, geocode, nfvctr, nseg, nvctr, keyv, keyg, mat_compr, &
               nat, ntypes, nzatom, nelpsp, atomnames, iatype, rxyz, on_which_atom)
      use module_base
      implicit none
      
      ! Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(out) :: nspin, nfvctr, nseg, nvctr
      character(len=1),intent(out) :: geocode
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
          read(iunit,*) nat, ntypes, nspin, geocode
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

end module sparsematrix_io
