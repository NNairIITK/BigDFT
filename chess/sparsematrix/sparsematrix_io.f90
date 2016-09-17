module sparsematrix_io
  use sparsematrix_base
  implicit none

  private

  public :: write_ccs_matrix
  public :: read_sparse_matrix
  public :: read_sparse_matrix_metadata
  public :: write_sparse_matrix
  public :: write_sparse_matrix_metadata
  public :: write_linear_coefficients
  public :: read_linear_coefficients
  public :: write_dense_matrix

  contains

    subroutine write_ccs_matrix(filename, nfvctr, nvctr, row_ind, col_ptr, mat_compr)
      use dynamic_memory
      use f_utils
      implicit none
      !Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(in) :: nfvctr !number of rows/columns
      integer,intent(in) :: nvctr !number of non-zero elements
      integer,dimension(nvctr),intent(in) :: row_ind
      integer,dimension(nfvctr),intent(in) :: col_ptr
      real(kind=mp),dimension(nvctr),intent(in) :: mat_compr
      ! Local variables
      integer :: i, iunit

      call f_routine(id='write_ccs_matrix')

      iunit = 99
      call f_open_file(iunit, file=trim(filename), binary=.false.)

      write(iunit,'(4(i0,1x))') nfvctr, nfvctr, nvctr, 0
      write(iunit,'(10000(i0,1x))') (col_ptr(i),i=1,nfvctr),nvctr+1
      write(iunit,'(10000(i0,1x))') (row_ind(i),i=1,nvctr)
      do i=1,nvctr
          write(iunit,'(es24.15)') mat_compr(i)
      end do

      call f_close(iunit)

      call f_release_routine()

    end subroutine write_ccs_matrix


    subroutine read_sparse_matrix(filename, nspin, nfvctr, nseg, nvctr, keyv, keyg, mat_compr)
      use dynamic_memory
      use f_utils
      implicit none
      
      ! Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(out) :: nspin, nfvctr, nseg, nvctr
      integer,dimension(:),pointer,intent(out) :: keyv
      integer,dimension(:,:,:),pointer,intent(out) :: keyg
      real(kind=mp),dimension(:),pointer,intent(out) :: mat_compr

      ! Local variables
      integer :: iunit, dummy_int, iseg, icol, irow, jorb, ind, ispin, iat, ntypes_, nat_, itype
      real(kind=mp) :: dummy_double
      character(len=20) :: dummy_char
      logical :: read_rxyz, read_on_which_atom

      call f_routine(id='read_sparse_matrix')


      iunit = 99
      call f_open_file(iunit, file=trim(filename), binary=.false.)

      read(iunit,*) nspin, nfvctr, nseg, nvctr
      keyv = f_malloc_ptr(nseg,id='keyv')
      keyg = f_malloc_ptr((/2,2,nseg/),id='keyg')

      do iseg=1,nseg
          read(iunit,*) keyv(iseg), keyg(1,1,iseg), keyg(2,1,iseg), keyg(1,2,iseg), keyg(2,2,iseg)
      end do

      mat_compr = f_malloc_ptr(nvctr*nspin,id='mat_compr')
      ind = 0
      do ispin=1,nspin
          do iseg=1,nseg
              icol = keyg(1,2,iseg)
              do jorb=keyg(1,1,iseg),keyg(2,1,iseg)
                  irow = jorb
                  ind = ind + 1
                  read(iunit,*) mat_compr(ind)
              end do
          end do
      end do


      call f_close(iunit)

      call f_release_routine()

    end subroutine read_sparse_matrix


    subroutine read_sparse_matrix_metadata(filename, nfvctr, nat, ntypes, units, geocode, cell_dim, &
               nzatom, nelpsp, atomnames, iatype, rxyz, on_which_atom)
      use dynamic_memory
      use f_utils
      implicit none
      
      ! Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(out) :: nfvctr
      character(len=20),intent(out) :: units
      character(len=1),intent(out) :: geocode
      real(kind=mp),dimension(3),intent(out) :: cell_dim
      integer,intent(out) :: nat, ntypes
      integer,dimension(:),pointer,intent(inout) :: nzatom, nelpsp, iatype
      character(len=20),dimension(:),pointer,intent(inout) :: atomnames
      real(kind=mp),dimension(:,:),pointer,intent(inout) :: rxyz
      integer,dimension(:),pointer,intent(inout) :: on_which_atom

      ! Local variables
      integer :: iunit, dummy_int, iseg, icol, irow, jorb, ind, ispin, iat, ntypes_, nat_, itype, i
      real(kind=mp) :: dummy_double
      character(len=20) :: dummy_char
      logical :: read_rxyz, read_on_which_atom

      call f_routine(id='read_sparse_matrix_metadata')

      iunit = 99
      call f_open_file(iunit, file=trim(filename), binary=.false.)

      read(iunit,*) nfvctr, nat, ntypes
      read(iunit,*) units
      read(iunit,*) geocode, cell_dim
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
      on_which_atom = f_malloc_ptr(nfvctr,id='on_which_atom')
      do i=1,nfvctr
          read(iunit,*) on_which_atom(i)
      end do

      call f_close(iunit)

      call f_release_routine()

    end subroutine read_sparse_matrix_metadata


    !> Write a sparse matrix to disk.
    !! ATTENTION: This routine must be called by all MPI tasks due to the fact that the matrix 
    !! is distributed among the matrix taksgroups
    subroutine write_sparse_matrix(iproc, nproc, comm, smat, mat, filename)
      use dynamic_memory
      use f_utils
      use sparsematrix, only: gather_matrix_from_taskgroups
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat
      character(len=*),intent(in) :: filename

      ! Local variables
      integer :: iunit, iseg, icol, irow, jorb, iat, jat, ind, ispin, itype
      real(kind=mp),dimension(:),allocatable :: matrix_compr

      call f_routine(id='write_sparse_matrix')

      matrix_compr = sparsematrix_malloc(smat,iaction=SPARSE_FULL,id='matrix_compr')
      call gather_matrix_from_taskgroups(iproc, nproc, comm, &
           smat, mat%matrix_compr, matrix_compr)

      if (iproc==0) then

          iunit = 99
          call f_open_file(iunit, file=trim(filename), binary=.false.)

          !!write(iunit,'(i10,2i6,3x,a,3es24.16,a)') nat, ntypes, smat%nspin, smat%geocode, smat%cell_dim, &
          !!    '   # number of atoms, number of atom types, nspin, geocode, cell_dim'
          !!do itype=1,ntypes
          !!    write(iunit,'(2i8,3x,a,a)') nzatom(itype), nelpsp(itype), trim(atomnames(itype)), &
          !!        '   # nz, nelpsp, name'
          !!end do
          !!do iat=1,nat
          !!    write(iunit,'(i5, 3es24.16,a,i0)') iatype(iat), rxyz(1:3,iat), '   # atom no. ',iat
          !!end do
          write(iunit,'(4i12,a)') smat%nspin, smat%nfvctr, smat%nseg, smat%nvctr, '   # nspin, nfvctr, nseg, nvctr'
          do iseg=1,smat%nseg
              write(iunit,'(5i12,a)') smat%keyv(iseg), smat%keyg(1,1,iseg), smat%keyg(2,1,iseg), &
                  smat%keyg(1,2,iseg), smat%keyg(2,2,iseg), '   # keyv, keyg(1,1), keyg(2,1), keyg(1,2), keyg(2,2)'
          end do
          ind = 0
          do ispin=1,smat%nspin
              do iseg=1,smat%nseg
                  icol = smat%keyg(1,2,iseg)
                  !iat = smat%on_which_atom(icol)
                  do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                      irow = jorb
                      !jat = smat%on_which_atom(irow)
                      ind = ind + 1
                      !write(iunit,'(es24.16,2i12,a)') matrix_compr(ind), jat, iat, '   # matrix, jat, iat'
                      write(iunit,'(es24.16,a,i0,a)') matrix_compr(ind),'   # matrix_compr(',ind,')'
                  end do
              end do
          end do

          call f_close(iunit)

          call f_free(matrix_compr)

          if (smat%smatmul_initialized) then
              iunit = 99
              call f_open_file(iunit, file=trim(filename)//'_matmul', binary=.false.)

              write(iunit,'(4i12,a)') smat%nspin, smat%nfvctr, smat%smmm%nseg, sum(smat%smmm%nvctr_par), &
                  '   # nspin, nfvctr, nseg, nvctr'
              do iseg=1,smat%smmm%nseg
                  write(iunit,'(5i12,a)') smat%smmm%keyv(iseg), smat%smmm%keyg(1,1,iseg), smat%smmm%keyg(2,1,iseg), &
                      smat%smmm%keyg(1,2,iseg), smat%smmm%keyg(2,2,iseg), &
                      '   # keyv, keyg(1,1), keyg(2,1), keyg(1,2), keyg(2,2)'
              end do
              ind = 0
              do ispin=1,smat%nspin
                  do iseg=1,smat%smmm%nseg
                      icol = smat%smmm%keyg(1,2,iseg)
                      !iat = smat%on_which_atom(icol)
                      do jorb=smat%smmm%keyg(1,1,iseg),smat%smmm%keyg(2,1,iseg)
                          irow = jorb
                          !jat = smat%on_which_atom(irow)
                          ind = ind + 1
                          !write(iunit,'(es24.16,2i12,a)') matrix_compr(ind), jat, iat, '   # matrix, jat, iat'
                          write(iunit,'(es24.16,a,i0,a)') 0123456789._mp,'   # matrix_compr(',ind,')'
                      end do
                  end do
              end do

              call f_close(iunit)
          end if

      end if

      call f_release_routine()

    end subroutine write_sparse_matrix


    subroutine write_sparse_matrix_metadata(iproc, nfvctr, nat, ntypes, units, geocode, cell_dim, iatype, &
               rxyz, nzatom, nelpsp, atomnames, on_which_atom, filename)
      use dynamic_memory
      use f_utils
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nfvctr, nat, ntypes
      character(len=*),intent(in) :: units
      character(len=1),intent(in) :: geocode
      real(kind=mp),dimension(3),intent(in) :: cell_dim
      integer,dimension(nat),intent(in) :: iatype
      real(kind=mp),dimension(3,nat),intent(in) :: rxyz
      integer,dimension(ntypes),intent(in) :: nzatom, nelpsp
      character(len=*),dimension(ntypes),intent(in) :: atomnames
      integer,dimension(nfvctr),intent(in) :: on_which_atom
      character(len=*),intent(in) :: filename

      ! Local variables
      integer :: iunit, iseg, icol, irow, jorb, iat, jat, ind, ispin, itype, i
      real(kind=mp),dimension(:),allocatable :: matrix_compr

      call f_routine(id='write_sparse_matrix_metadata')

      if (iproc==0) then

          iunit = 99
          call f_open_file(iunit, file=trim(filename), binary=.false.)

          write(iunit,'(i10,i8,i6,a)') nfvctr, nat, ntypes, &
              '   # matrix dimension, number of atoms, number of atom types'
          write(iunit,'(a,a)') trim(units), '  # units'
          write(iunit,'(a,3es24.16,a)') geocode, cell_dim, &
              '   # geocode, cell_dim'
          do itype=1,ntypes
              write(iunit,'(2i8,3x,a,a)') nzatom(itype), nelpsp(itype), trim(atomnames(itype)), &
                  '   # nz, nelpsp, name'
          end do
          do iat=1,nat
              write(iunit,'(i5,3es24.16,a,i0)') iatype(iat), rxyz(1:3,iat), '   # atom no. ',iat
          end do
          do i=1,nfvctr
              write(iunit,'(i8,a,i0,a)') on_which_atom(i), '   # on_which_atom(',i,')'
          end do

          call f_close(iunit)

      end if

      call f_release_routine()

    end subroutine write_sparse_matrix_metadata


    subroutine write_linear_coefficients(iproc, iroot, filename, verbosity, nat, rxyz, iatype, ntypes, nzatom, &
               nelpsp, atomnames, nfvctr, ntmb, nspin, coeff, eval)
      use yaml_output
      implicit none
      ! Calling arguments
      character(len=*),intent(in) :: filename
      !type(atoms_data),intent(in) :: at
      integer,intent(in) :: iproc, iroot, nat, verbosity, ntypes, nfvctr, ntmb, nspin
      real(mp), dimension(3,nat), intent(in) :: rxyz
      integer,dimension(nat),intent(in) :: iatype
      integer,dimension(ntypes),intent(in) :: nzatom, nelpsp
      character(len=20),dimension(ntypes),intent(in) :: atomnames
      real(mp), dimension(nfvctr,ntmb), intent(in) :: coeff
      real(mp), dimension(ntmb), intent(in) :: eval
      ! Local variables
      integer :: iunit, itype, iat, i, j
      logical :: scaled

      call f_routine(id='write_linear_coefficients')


      if (iproc==iroot) then

          iunit = 99
          call f_open_file(iunit, file=trim(filename), binary=.false.)
    
          ! Write the Header
          !write(iunit,'(i10,2i6,a)') at%astruct%nat, at%astruct%ntypes, nspin, &
          write(iunit,'(i10,2i6,a)') nat, ntypes, nspin, &
              '   # number of atoms, number of atom types, nspin'
          !do itype=1,at%astruct%ntypes
          do itype=1,ntypes
              !write(iunit,'(2i8,3x,a,a)') at%nzatom(itype), at%nelpsp(itype), trim(at%astruct%atomnames(itype)), &
              write(iunit,'(2i8,3x,a,a)') nzatom(itype), nelpsp(itype), trim(atomnames(itype)), &
                  '   # nz, nelpsp, name'
          end do
          !do iat=1,at%astruct%nat
          do iat=1,nat
              !write(iunit,'(i5, 3es24.16,a,i0)') at%astruct%iatype(iat), rxyz(1:3,iat), '   # atom no. ',iat
              write(iunit,'(i5, 3es24.16,a,i0)') iatype(iat), rxyz(1:3,iat), '   # atom no. ',iat
          end do
          write(iunit,'(2i12,a)') nfvctr, ntmb, '   # nfvctr, ntmb'
          do i=1,ntmb
              write(iunit,'(es24.16,a,i0)') eval(i), '   # eval no. ', i
          enddo
    
          ! Now write the coefficients
          do i=1,ntmb
             ! First element always positive, for consistency when using for transfer integrals;
             ! unless 1st element below some threshold, in which case first significant element.
             scaled = .false.
             do j=1,nfvctr
                if (abs(coeff(j,i))>1.0d-3) then
                   if (coeff(j,i)<0.0_mp) call dscal(ntmb,-1.0_mp,coeff(1,i),1)
                   scaled = .true.
                   exit
                end if
             end do
             if (.not.scaled) then
                 call yaml_warning('Consistency between the written coefficients not guaranteed')
             end if
    
             do j = 1,nfvctr
                 write(iunit,'(es24.16,2i9,a)') coeff(j,i), j, i, '   # coeff, j, i'
             end do
          end do  
          if (verbosity >= 2 .and. iproc==0) call yaml_map('Wavefunction coefficients written',.true.)

          call f_close(iunit)

      end if

      call f_release_routine()
    
    end subroutine write_linear_coefficients


    subroutine read_linear_coefficients(filename, nspin, nfvctr, ntmb, coeff, &
               nat, ntypes, nzatom, nelpsp, iatype, atomnames, rxyz, eval)
      use yaml_output
      implicit none
      ! Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(out) :: nspin, nfvctr, ntmb
      real(kind=8),dimension(:,:),pointer,intent(inout) :: coeff
      integer,intent(out),optional :: nat, ntypes
      integer,dimension(:),pointer,intent(inout),optional :: nzatom, nelpsp, iatype
      character(len=20),dimension(:),pointer,intent(inout),optional :: atomnames
      real(kind=8),dimension(:,:),pointer,intent(inout),optional :: rxyz
      real(kind=8),dimension(:),pointer,intent(inout),optional :: eval
      ! Local variables
      real(kind=8) :: dummy_double
      character(len=20) :: dummy_char
      integer :: iunit, itype, iat, i, j, dummy_int, ntypes_, nat_
      logical :: scaled, read_rxyz, read_eval

      call f_routine(id='read_linear_coefficients')

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

      if (present(eval)) then
          read_eval = .true.
      else
          read_eval = .false.
      end if

      iunit = 99
      call f_open_file(iunit, file=trim(filename), binary=.false.)
    
      ! Read the Header
      if (read_rxyz) then
          read(iunit,*) nat, ntypes, nspin
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
          read(iunit,*) nat_, ntypes_, nspin
          do itype=1,ntypes_
              read(iunit,*) dummy_int, dummy_int, dummy_char
          end do
          do iat=1,nat_
              read(iunit,*) dummy_int, dummy_double, dummy_double, dummy_double
          end do
      end if

      read(iunit,*) nfvctr, ntmb

      if (read_eval) then
          eval = f_malloc_ptr(ntmb,id='eval')
          do i=1,ntmb
              read(iunit,*) eval(i)
          end do
      else
          do i=1,ntmb
              read(iunit,*) dummy_double
          end do
      end if
    
      ! Now read the coefficients
      coeff = f_malloc_ptr((/nfvctr,ntmb/),id='coeff')

      do i=1,ntmb

         do j = 1,nfvctr
             read(iunit,*) coeff(j,i)
         end do

         ! First element always positive, for consistency when using for transfer integrals;
         ! unless 1st element below some threshold, in which case first significant element.
         scaled = .false.
         do j=1,nfvctr
            if (abs(coeff(j,i))>1.0d-3) then
               if (coeff(j,i)<0.0_mp) call dscal(ntmb,-1.0_mp,coeff(1,i),1)
               scaled = .true.
               exit
            end if
         end do
         if (.not.scaled) then
             call yaml_warning('Consistency between the written coefficients not guaranteed')
         end if
    
      end do  

      call f_close(iunit)


      call f_release_routine()
    
    end subroutine read_linear_coefficients


    !> Write a sparse matrix to disk, but in dense format
    !! ATTENTION: This routine must be called by all MPI tasks due to the fact that the matrix 
    !! in distributed among the matrix taksgroups
    subroutine write_dense_matrix(iproc, nproc, comm, smat, mat, filename, binary)
      use sparsematrix_base
      use sparsematrix, only: uncompress_matrix2
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(inout) :: mat
      character(len=*),intent(in) :: filename
      logical, intent(in) :: binary

      ! Local variables
      integer :: iunit, iseg, icol, irow, jorb, iat, jat, ind, ispin, itype, iorb
      real(kind=8),dimension(:),allocatable :: matrix_compr

      call f_routine(id='write_dense_matrix')


      mat%matrix = sparsematrix_malloc_ptr(smat, iaction=DENSE_FULL, id='mat%matrix')
      call uncompress_matrix2(iproc, nproc, comm, &
           smat, mat%matrix_compr, mat%matrix)

      if (iproc==0) then

          iunit = 99
          call f_open_file(iunit, file=trim(filename), binary=binary)

          !write(iunit,'(i10,2i6,a)') nat, ntypes, smat%nspin, &
          !    '   # number of atoms, number of atom types, nspin'
          !do itype=1,ntypes
          !    write(iunit,'(2i8,3x,a,a)') nzatom(itype), nelpsp(itype), trim(atomnames(itype)), &
          !        '   # nz, nelpsp, name'
          !end do
          !do iat=1,nat
          !    write(iunit,'(i5, 3es24.16,a,i0)') iatype(iat), rxyz(1:3,iat), '   # atom no. ',iat
          !end do
          !write(iunit,'(3i12,a)') smat%nfvctr, smat%nseg, smat%nvctr, '   # nfvctr, nseg, nvctr'
          !do iseg=1,smat%nseg
          !    write(iunit,'(5i12,a)') smat%keyv(iseg), smat%keyg(1,1,iseg), smat%keyg(2,1,iseg), &
          !        smat%keyg(1,2,iseg), smat%keyg(2,2,iseg), '   # keyv, keyg(1,1), keyg(2,1), keyg(1,2), keyg(2,2)'
          !end do
          !ind = 0
          !do ispin=1,smat%nspin
          !    do iseg=1,smat%nseg
          !        icol = smat%keyg(1,2,iseg)
          !        iat = smat%on_which_atom(icol)
          !        do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
          !            irow = jorb
          !            jat = smat%on_which_atom(irow)
          !            ind = ind + 1
          !            write(iunit,'(es24.16,2i12,a)') matrix_compr(ind), jat, iat, '   # matrix, jat, iat'
          !        end do
          !    end do
          !end do

          !unify the structure with write_sparse?
          !!!if (.not. binary) then
          !!!    write(iunit,'(a,3i10,a)') '#  ',smat%nfvctr, nat, smat%nspin, &
          !!!        '    number of basis functions, number of atoms, number of spins'
          !!!else
          !!!    write(iunit) '#  ',smat%nfvctr, nat, smat%nspin, &
          !!!        '    number of basis functions, number of atoms, number of spins'
          !!!end if
          !!!do iat=1,nat
          !!!    if (.not. binary) then
          !!!        write(iunit,'(a,3es24.16,a,i4.4)') '#  ',rxyz(1:3,iat), '   # position of atom no. ',iat
          !!!    else
          !!!        write(iunit) '#  ',rxyz(1:3,iat)
          !!!    end if
          !!!end do
          if (.not. binary) then
              write(iunit,'(2i10,a)') smat%nspin, smat%nfvctr, ' # nspin, nfvctr'
          else
              write(iunit) smat%nspin, smat%nfvctr, ' # nspin, nfvctr'
          end if
    
          do ispin=1,smat%nspin
             !!do iorb=1,smat%nfvctr
             !!   iat=onwhichatom(iorb)
             !!   do jorb=1,smat%nfvctr
             !!      jat=onwhichatom(jorb)
             !!      if (.not. binary) then
             !!         write(iunit,'(2(i6,1x),es19.12,2(1x,i6),a)') iorb,jorb,mat%matrix(iorb,jorb,ispin),iat,jat, &
             !!             '   # i, j, mat(i,j), iat, jat'
             !!      else
             !!         write(iunit) iorb,jorb,mat%matrix(iorb,jorb,ispin),iat,jat
             !!      end if
             !!   end do
             !!end do
             do iorb=1,smat%nfvctr
                do jorb=1,smat%nfvctr
                   if (.not. binary) then
                      write(iunit,'(2(i6,1x),es19.12,a)') iorb,jorb,mat%matrix(iorb,jorb,ispin), &
                          '   # i, j, mat(i,j)'
                   else
                      write(iunit) iorb,jorb,mat%matrix(iorb,jorb,ispin)
                   end if
                end do
             end do
          end do

          call f_close(iunit)

          call f_free_ptr(mat%matrix)
      end if

      call f_release_routine()

    end subroutine write_dense_matrix

end module sparsematrix_io
