!> @file
!!  
!! @author
!!   Copyright (C) 2016 CheSS developers
!!
!!   This file is part of CheSS.
!!   
!!   CheSS is free software: you can redistribute it and/or modify
!!   it under the terms of the GNU Lesser General Public License as published by
!!   the Free Software Foundation, either version 3 of the License, or
!!   (at your option) any later version.
!!   
!!   CheSS is distributed in the hope that it will be useful,
!!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU Lesser General Public License for more details.
!!   
!!   You should have received a copy of the GNU Lesser General Public License
!!   along with CheSS.  If not, see <http://www.gnu.org/licenses/>.


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
  use io, only: read_sparse_matrix, write_ccs_matrix
  use module_atoms,only: atoms_data, atoms_data_null, deallocate_atoms_data
  use sparsematrix_base, only: sparse_matrix, matrices, sparsematrix_malloc_ptr, &
                               assignment(=), SPARSE_FULL, SPARSE_TASKGROUP, matrices_null, &
                               deallocate_sparse_matrix, deallocate_matrices
  use bigdft_run, only: bigdft_init
  use sparsematrix_init, only: bigdft_to_sparsebigdft, distribute_columns_on_processes_simple
  use sparsematrix, only: transform_sparsity_pattern
  use foe_base, only: foe_data, foe_data_deallocate
  use foe_common, only: init_foe
  use foe, only: fermi_operator_expansion
  use ice, only: inverse_chebyshev_expansion_new
  implicit none

  ! Variables
  character(len=*),parameter :: filename='matrix.dat'
  integer :: nspin, nfvctr, nvctr, nseg, i, iseg
  character(len=1) :: geocode
  integer,dimension(:),pointer :: keyv
  integer,dimension(:,:,:),pointer :: keyg
  real(kind=8),dimension(:),pointer :: mat_compr
  integer,dimension(:),pointer :: row_ind, col_ptr

  integer :: nspin_h, nfvctr_h, nseg_h, nvctr_h, nfvctrp_h, isfvctr_h
  integer :: nspin_s, nfvctr_s, nseg_s, nvctr_s, nfvctrp_s, isfvctr_s
  integer :: nspin_l, nfvctr_l, nseg_l, nvctr_l, nfvctrp_l, isfvctr_l
  integer :: norder_polynomial
  integer,dimension(:),pointer :: keyv_h, keyv_s, keyv_l, on_which_atom_h, on_which_atom_s, on_which_atom_l
  integer,dimension(:,:,:),pointer :: keyg_h, keyg_s, keyg_l
  real(kind=8),dimension(:),pointer :: matrix_compr_h, matrix_compr_s, matrix_compr_l, matrix_compr_sl
  type(atoms_data) :: at
  character(len=1) :: geocode_h, geocode_s, geocode_l
  type(sparse_matrix) :: smat_h, smat_s, smat_l
  type(matrices) :: ham, overlap, kernel
  type(matrices),dimension(1) :: inv_overlap
  real(kind=8),dimension(2) :: charge
  real(kind=8) :: tmprtr, evlow, evhigh, fscale, ebs
  real(kind=8) :: ef_interpol_det, ef_interpol_chargediff, fscale_lowerbound, fscale_upperbound
  real(kind=8) :: max_inversion_error
  integer :: evbounds_nsatur, evboundsshrink_nsatur, foe_verbosity, order_taylor
  type(foe_data) :: foe_obj
  logical :: calculate_minusonehalf
  character(len=4) :: label

  call f_lib_initialize()

  call bigdft_init()

  call sparse_matrix_and_matrices_init_from_file_ccs('overlap_ccs.dat', bigdft_mpi%iproc, bigdft_mpi%nproc, smat_s, overlap)
  call sparse_matrix_and_matrices_init_from_file_ccs('hamiltonian_ccs.dat', bigdft_mpi%iproc, bigdft_mpi%nproc, smat_h, ham)
  call sparse_matrix_and_matrices_init_from_file_ccs('kernel_ccs.dat', bigdft_mpi%iproc, bigdft_mpi%nproc, smat_l, kernel)
  inv_overlap(1) = matrices_null()
  inv_overlap(1)%matrix_compr = sparsematrix_malloc_ptr(smat_l, iaction=SPARSE_FULL, id='inv_overlap%matrix_compr')

  order_taylor = 1020
  call inverse_chebyshev_expansion_new(bigdft_mpi%iproc, bigdft_mpi%nproc, norder_polynomial, &
       smat_s, smat_l, 1, (/-1.d0/), overlap, inv_overlap)

!!!FOE  charge = 10.d0
!!!FOE  tmprtr = 0.d0
!!!FOE  evbounds_nsatur = 100
!!!FOE  evboundsshrink_nsatur = 4
!!!FOE  evlow = -1.0d0
!!!FOE  evhigh = 1.0d0
!!!FOE  fscale = 1.d-2
!!!FOE  ef_interpol_det = 1.d-12
!!!FOE  ef_interpol_chargediff = 10.d0
!!!FOE  fscale_lowerbound = 5.d-3
!!!FOE  fscale_upperbound = 5.d-2
!!!FOE  nspin = 1
!!!FOE  call init_foe(bigdft_mpi%iproc, bigdft_mpi%nproc, nspin, charge, tmprtr, evbounds_nsatur, evboundsshrink_nsatur, &
!!!FOE       evlow, evhigh, fscale, ef_interpol_det, ef_interpol_chargediff, &
!!!FOE       fscale_lowerbound, fscale_upperbound, foe_obj)
!!!FOE
!!!FOE  max_inversion_error = 1.d-8
!!!FOE  calculate_minusonehalf = .true.
!!!FOE  foe_verbosity = 1
!!!FOE  label = 'test'
!!!FOE  call fermi_operator_expansion(bigdft_mpi%iproc, bigdft_mpi%nproc, &
!!!FOE       ebs, order_taylor, max_inversion_error, &
!!!FOE       calculate_minusonehalf, foe_verbosity, &
!!!FOE       label, smat_s, smat_h, smat_l, ham, overlap, inv_overlap, kernel, foe_obj)
!!!FOE
!!!FOE  call foe_data_deallocate(foe_obj)

  call deallocate_sparse_matrix(smat_s)
  call deallocate_sparse_matrix(smat_h)
  call deallocate_sparse_matrix(smat_l)
  call deallocate_matrices(ham)
  call deallocate_matrices(overlap)
  call deallocate_matrices(kernel)
  call deallocate_matrices(inv_overlap(1))

  call f_lib_finalize()


end program driver


subroutine sparse_matrix_and_matrices_init_from_file_ccs(filename, iproc, nproc, smat, mat)
  use module_base
  use sparsematrix_base, only: sparse_matrix, matrices
  use sparsematrix_init, only: read_ccs_format, ccs_to_sparsebigdft_short, &
                               bigdft_to_sparsebigdft
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
