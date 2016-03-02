!> @file
!! Test of the sparsematrix library
!! @author
!!    Copyright (C) 2015-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


program driver
  use module_base
  use sparsematrix_init, only: sparsebigdft_to_ccs, read_ccs_format
  use sparsematrix_io, only: read_sparse_matrix, write_ccs_matrix
  use module_atoms, only: atoms_data, atoms_data_null
  use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc_ptr, &
                               assignment(=), SPARSE_FULL
  use bigdft_run, only: bigdft_init
  use sparsematrix_init, only: bigdft_to_sparsebigdft, distribute_columns_on_processes_simple
  use sparsematrix, only: transform_sparsity_pattern
  implicit none

  ! Variables
  character(len=*), parameter :: filename='matrix.dat'
  !integer :: nspin, nfvctr, nvctr, nseg, i, iseg
  !character(len=1) :: geocode
  !integer, dimension(:), pointer :: keyv
  !integer, dimension(:,:,:), pointer :: keyg
  !real(kind=8), dimension(:), pointer :: mat_compr
  integer, dimension(:), pointer :: row_ind, col_ptr

  integer :: nspin_h, nfvctr_h, nseg_h, nvctr_h!, nfvctrp_h, isfvctr_h
  integer :: nspin_s, nfvctr_s, nseg_s, nvctr_s!, nfvctrp_s, isfvctr_s
  integer, dimension(:), pointer :: keyv_h, keyv_s, on_which_atom_h, on_which_atom_s
  integer, dimension(:,:,:), pointer :: keyg_h, keyg_s
  real(kind=8), dimension(:), pointer :: matrix_compr_h, matrix_compr_s, matrix_compr_sl
  type(atoms_data) :: at
  character(len=1) :: geocode_h, geocode_s
  type(sparse_matrix) :: smat_h, smat_s

  call f_lib_initialize()

  call bigdft_init()

  !call read_sparse_matrix('hamiltonian_sparse.bin', nspin_h, geocode_h, nfvctr_h, nseg_h, nvctr_h, keyv_h, keyg_h, mat_compr_h)
  !call read_sparse_matrix('overlap_sparse.bin', nspin_s, geocode_s, nfvctr_s, nseg_s, nvctr_s, keyv_s, keyg_s, mat_compr_s)

  at = atoms_data_null()

  call read_sparse_matrix('hamiltonian_sparse.bin', nspin_h, geocode_h, nfvctr_h, nseg_h, nvctr_h, keyv_h, keyg_h, &
       matrix_compr_h, at%astruct%nat, at%astruct%ntypes, at%nzatom, at%nelpsp, &
       at%astruct%atomnames, at%astruct%iatype, at%astruct%rxyz, on_which_atom=on_which_atom_h)
  at%refcnt=f_ref_new('atoms')
  !call distribute_columns_on_processes_simple(bigdft_mpi%iproc, bigdft_mpi%nproc, nfvctr_h, nfvctrp_h, isfvctr_h)
  call bigdft_to_sparsebigdft(bigdft_mpi%iproc, bigdft_mpi%nproc, nfvctr_h, nvctr_h, nseg_h, keyg_h, smat_h, &
       nspin_h, geocode_h, on_which_atom_h)
  call f_free_ptr(keyv_h)
  call f_free_ptr(keyg_h)

  call read_sparse_matrix('overlap_sparse.bin', nspin_s, geocode_s, nfvctr_s, nseg_s, nvctr_s, keyv_s, keyg_s, &
       matrix_compr_s, on_which_atom=on_which_atom_s)
  !call distribute_columns_on_processes_simple(bigdft_mpi%iproc, bigdft_mpi%nproc, nfvctr_s, nfvctrp_s, isfvctr_s)
  call bigdft_to_sparsebigdft(bigdft_mpi%iproc, bigdft_mpi%nproc, nfvctr_s, nvctr_s, nseg_s, keyg_s, smat_s, &
       nspin_s, geocode_s, on_which_atom_s)
  call f_free_ptr(keyv_s)
  call f_free_ptr(keyg_s)


  matrix_compr_sl = sparsematrix_malloc_ptr(smat_h, iaction=SPARSE_FULL, id='matrix_compr_sl')
  call transform_sparsity_pattern(smat_h%nfvctr, smat_s%smmm%nvctrp_mm, smat_s%smmm%isvctr_mm, &
       smat_s%nseg, smat_s%keyv, smat_s%keyg, smat_s%smmm%line_and_column_mm, &
       smat_h%smmm%nvctrp, smat_h%smmm%isvctr, smat_h%smmm%nseg, smat_h%smmm%keyv, smat_h%smmm%keyg, &
       smat_h%smmm%istsegline, 'small_to_large', matrix_compr_s, matrix_compr_sl)

  row_ind = f_malloc_ptr(smat_h%nvctr,id='row_ind')
  col_ptr = f_malloc_ptr(smat_h%nfvctr,id='col_ptr')
  call sparsebigdft_to_ccs(smat_h%nfvctr, smat_h%nvctr, smat_h%nseg, smat_h%keyg, row_ind, col_ptr)

  call write_ccs_matrix('overlap_sparse_PEXSI.bin', smat_h%nfvctr, smat_h%nvctr, row_ind, col_ptr, matrix_compr_sl)
  call write_ccs_matrix('hamiltonian_sparse_PEXSI.bin', smat_h%nfvctr, smat_h%nvctr, row_ind, col_ptr, matrix_compr_h)


  !!!OLD   call read_sparse_matrix(filename, nspin, geocode, nfvctr, nseg, nvctr, keyv, keyg, mat_compr)

  !!!OLD   do iseg=1,nseg
  !!!OLD       write(*,*) 'iseg, keyv, keyg', iseg, keyv(iseg), keyg(1:2,1:2,iseg)
  !!!OLD   end do

  !!!OLD   row_ind = f_malloc_ptr(nvctr,id='row_ind')
  !!!OLD   col_ptr = f_malloc_ptr(nfvctr,id='col_ptr')
  !!!OLD   call sparsebigdft_to_ccs(nfvctr, nvctr, nseg, keyg, row_ind, col_ptr)

  !!!OLD   call write_ccs_matrix('matrix_ccs.dat', nfvctr, nvctr, row_ind, col_ptr, mat_compr)
  !!!OLD
  !!!OLD   !do i=1,nvctr
  !!!OLD   !    write(*,*) 'i, row_ind', i, row_ind(i)
  !!!OLD   !end do
  !!!OLD   !do i=1,nfvctr
  !!!OLD   !    write(*,*) 'i, col_ptr', i, col_ptr(i)
  !!!OLD   !end do

  !!!OLD   call f_free_ptr(keyv)
  !!!OLD   call f_free_ptr(keyg)
  !!!OLD   call f_free_ptr(mat_compr)

  !!!OLD   call ccs_to_sparsebigdft_short(nfvctr, nvctr, row_ind, col_ptr, nseg, keyv, keyg)

  !!!OLD   write(*,*) 'nseg',nseg
  !!!OLD   do iseg=1,nseg
  !!!OLD       write(*,*) 'iseg, keyv, keyg', iseg, keyv(iseg), keyg(1:2,1:2,iseg)
  !!!OLD   end do

  !!!OLD   call f_free_ptr(row_ind)
  !!!OLD   call f_free_ptr(col_ptr)
  !!!OLD   call f_free_ptr(keyv)
  !!!OLD   call f_free_ptr(keyg)


  !!!OLD   !! NEW TEST
  !!!OLD   !call read_ccs_format('lap2dr.matrix', nfvctr, nvctr, col_ptr, row_ind, mat_compr)
  !!!OLD   !call write_ccs_matrix('matrix_ccs.dat', nfvctr, nvctr, row_ind, col_ptr, mat_compr)

  call f_lib_finalize()


end program driver
