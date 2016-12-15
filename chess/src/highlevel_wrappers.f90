!> @file
!!   File containing high level wrappers to test CheSS
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


module highlevel_wrappers
  use futile
  use wrapper_mpi
  use sparsematrix_base
  private

  !> Public routines
  public :: calculate_eigenvalues
  public :: solve_eigensystem_lapack

  contains

    subroutine calculate_eigenvalues(iproc, nproc, matrix_format, metadata_file, &
               overlap_file, hamiltonian_file, kernel_file, kernel_matmul_file, &
               iev_minx, iev_maxx, fscale)
      use sparsematrix_highlevel, only: sparse_matrix_metadata_init_from_file, &
                                        sparse_matrix_and_matrices_init_from_file_bigdft, &
                                        matrices_init, &
                                        get_selected_eigenvalues_from_FOE
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      character(len=*),intent(in) :: matrix_format, metadata_file
      character(len=*),intent(in) :: overlap_file, hamiltonian_file, kernel_file, kernel_matmul_file
      integer,intent(in) :: iev_minx, iev_maxx
      real(mp),intent(in) :: fscale

      ! Local variables
      integer :: iev, iev_min, iev_max
      type(sparse_matrix_metadata) :: smmd
      type(sparse_matrix) :: smat_s, smat_m, smat_l
      type(matrices),dimension(1) :: ovrlp_minus_one_half
      type(matrices) :: ovrlp_mat, hamiltonian_mat, kernel_mat
      real(kind=8),dimension(:),allocatable :: eval
      external :: gather_timings


      call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
      call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(overlap_file), &
           iproc, nproc, mpiworld(), smat_s, ovrlp_mat, &
           init_matmul=.false.)
      call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(hamiltonian_file), &
           iproc, nproc, mpiworld(), smat_m, hamiltonian_mat, &
           init_matmul=.false.)
      call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(kernel_file), &
           iproc, nproc, mpiworld(), smat_l, kernel_mat, &
           init_matmul=.true., filename_mult=trim(kernel_matmul_file))
      call matrices_init(smat_l, ovrlp_minus_one_half(1))

      !call timing(mpiworld(),'INIT','PR')
      call f_timing_checkpoint(ctr_name='INIT',mpi_comm=mpiworld(),nproc=mpisize(),&
                   gather_routine=gather_timings)
      iev_min = iev_minx
      iev_max = iev_maxx
      if (iev_min<1 .or. iev_min>smat_s%nfvctr .or. iev_max>smat_s%nfvctr .or. iev_max<1) then
          if (iproc==0) then
              call yaml_warning('The required eigenvalues are outside of the possible range, automatic ajustment')
          end if
      end if
      iev_min = max(iev_min,1)
      iev_min = min(iev_min,smat_s%nfvctr)
      iev_max = min(iev_max,smat_s%nfvctr)
      iev_max = max(iev_max,1)
      eval = f_malloc(iev_min.to.iev_max,id='eval')
      if (iproc==0) then
          call yaml_mapping_open('Calculating eigenvalues using FOE')
      end if
      call get_selected_eigenvalues_from_FOE(iproc, nproc, mpiworld(), &
           iev_min, iev_max, smat_s, smat_m, smat_l, ovrlp_mat, hamiltonian_mat, &
           ovrlp_minus_one_half, eval, fscale)

      !!call timing(mpiworld(),'CALC','PR')
      call f_timing_checkpoint(ctr_name='CALC',mpi_comm=mpiworld(),nproc=mpisize(),&
                   gather_routine=gather_timings)

      if (iproc==0) then
          call yaml_sequence_open('values')
          do iev=iev_min,iev_max
              call yaml_sequence(advance='no')
              call yaml_mapping_open(flow=.true.)
              call yaml_map('ID',iev,fmt='(i6.6)')
              call yaml_map('eval',eval(iev),fmt='(es12.5)')
              call yaml_mapping_close()
          end do
          call yaml_sequence_close()
          call yaml_mapping_close()
      end if


      call deallocate_sparse_matrix(smat_s)
      call deallocate_sparse_matrix(smat_m)
      call deallocate_sparse_matrix(smat_l)
      call deallocate_matrices(ovrlp_mat)
      call deallocate_matrices(hamiltonian_mat)
      call deallocate_matrices(kernel_mat)
      call deallocate_matrices(ovrlp_minus_one_half(1))
      call deallocate_sparse_matrix_metadata(smmd)

      call f_free(eval)

      call f_timing_checkpoint(ctr_name='LAST',mpi_comm=mpiworld(),nproc=mpisize(),&
           gather_routine=gather_timings)

    end subroutine calculate_eigenvalues


    subroutine solve_eigensystem_lapack(iproc, nproc, matrix_format, metadata_file, &
               overlap_file, hamiltonian_file, scalapack_blocksize, write_output, &
               coeff_file, evals_out, coeffs_out)
      use sparsematrix, only: uncompress_matrix, &
                              diagonalizehamiltonian2
      use sparsematrix_io, only: write_linear_coefficients
      use sparsematrix_highlevel, only: sparse_matrix_metadata_init_from_file, &
                                        sparse_matrix_and_matrices_init_from_file_bigdft
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, scalapack_blocksize
      character(len=*),intent(in) :: matrix_format, metadata_file
      character(len=*),intent(in) :: overlap_file, hamiltonian_file
      logical :: write_output
      character(len=*),intent(in),optional :: coeff_file
      real(mp),dimension(:),pointer,optional :: evals_out
      real(mp),dimension(:,:,:),pointer,optional :: coeffs_out

      ! Local variables
      integer :: iunit
      type(sparse_matrix_metadata) :: smmd
      type(sparse_matrix) :: smat_s, smat_m
      type(matrices) :: ovrlp_mat, hamiltonian_mat
      real(kind=8),dimension(:),allocatable :: eval
      external :: gather_timings

      call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(overlap_file), &
           iproc, nproc, mpiworld(), smat_s, ovrlp_mat, &
           init_matmul=.false.)
      call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
      call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, trim(hamiltonian_file), &
           iproc, nproc, mpiworld(), smat_m, hamiltonian_mat, &
           init_matmul=.false.)

      ovrlp_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_FULL, id='ovrlp_mat%matrix')
      call uncompress_matrix(iproc, nproc, &
           smat_s, inmat=ovrlp_mat%matrix_compr, outmat=ovrlp_mat%matrix)
      hamiltonian_mat%matrix = sparsematrix_malloc_ptr(smat_s, iaction=DENSE_FULL, id='hamiltonian_mat%matrix')
      call uncompress_matrix(iproc, nproc, &
           smat_m, inmat=hamiltonian_mat%matrix_compr, outmat=hamiltonian_mat%matrix)
      eval = f_malloc(smat_s%nfvctr,id='eval')

      if (iproc==0) then
          call yaml_comment('Diagonalizing the matrix',hfill='~')
      end if
      call diagonalizeHamiltonian2(iproc, nproc, mpiworld(), scalapack_blocksize, &
           smat_s%nfvctr, hamiltonian_mat%matrix, ovrlp_mat%matrix, eval)
      if (iproc==0) then
          call yaml_comment('Matrix successfully diagonalized',hfill='~')
      end if

      if (write_output) then
          if (.not.present(coeff_file)) then
              call f_err_throw("'coeff_file' is not present")
          end if
          iunit=99
          call f_open_file(iunit, file=trim(coeff_file), binary=.false.)
          call write_linear_coefficients(iproc, 0, trim(coeff_file), 2, smmd%nat, smmd%rxyz, &
               smmd%iatype, smmd%ntypes, smmd%nzatom, &
               smmd%nelpsp, smmd%atomnames, smat_s%nfvctr, &
               smat_s%nfvctr, smat_s%nspin, hamiltonian_mat%matrix, eval)
          call f_close(iunit)
      end if

      if (present(evals_out)) then
          call f_memcpy(src=eval, dest=evals_out)
      end if
      if (present(coeffs_out)) then
          call f_memcpy(src=hamiltonian_mat%matrix, dest=coeffs_out)
      end if

      call f_free(eval)
      call deallocate_matrices(ovrlp_mat)
      call deallocate_matrices(hamiltonian_mat)
      call deallocate_sparse_matrix(smat_s)
      call deallocate_sparse_matrix(smat_m)
      call deallocate_sparse_matrix_metadata(smmd)
      !call f_free_ptr(rxyz)
      !call f_free_ptr(iatype)
      !call f_free_ptr(nzatom)
      !call f_free_ptr(nelpsp)
      !call f_free_str_ptr(len(atomnames),atomnames)

    end subroutine solve_eigensystem_lapack

end module highlevel_wrappers
