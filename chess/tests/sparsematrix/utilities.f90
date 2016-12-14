!> @file
!!   File containing helper routines for the test drivers
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


module utilities
  private

  public :: get_ccs_data_from_file
  public :: calculate_error

  contains


    subroutine get_ccs_data_from_file(filename, nfvctr, nvctr, row_ind, col_ptr)
      use sparsematrix_init, only: read_ccs_format
      use dynamic_memory, only: f_free_ptr
      implicit none

      ! Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(out) :: nfvctr, nvctr
      integer,dimension(:),pointer,intent(out) :: row_ind
      integer,dimension(:),pointer,intent(out) :: col_ptr

      ! Local variables
      real(kind=8),dimension(:),pointer :: val

      call read_ccs_format(filename, nfvctr, nvctr, col_ptr, row_ind, val)

      call f_free_ptr(val)

    end subroutine get_ccs_data_from_file


    subroutine calculate_error(iproc, smat, mat, mat_ref, nthreshold, threshold, check_full_matrix, header)
      use futile
      use sparsematrix_base
      use sparsematrix_types, only: sparse_matrix, matrices
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nthreshold
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat, mat_ref
      real(mp),dimension(nthreshold),intent(in) :: threshold
      logical,intent(in) :: check_full_matrix
      character(len=*),intent(in) :: header
    
      ! Local variables
      real(mp) :: max_error, mean_error, max_error_rel, mean_error_rel, tt, tt_rel, tt_ref
      real(mp),dimension(nthreshold) :: max_error_rel_threshold, mean_error_rel_threshold
      integer,dimension(nthreshold) :: nrel_threshold
      integer :: i, j, ithreshold
    
    
      ! Sparse matrices
      call set_to_zero()
      do i=1,smat%nvctr
          tt = abs(mat%matrix_compr(i)-mat_ref%matrix_compr(i))
          tt_rel = tt/abs(mat_ref%matrix_compr(i))
          tt_ref = mat_ref%matrix_compr(i)
          call calculate_errors()
      end do
      mean_error = mean_error/real(smat%nvctr,kind=8)
      mean_error_rel = mean_error_rel/real(smat%nvctr,kind=8)
      do ithreshold=1,nthreshold
          mean_error_rel_threshold(ithreshold) = mean_error_rel_threshold(ithreshold)/real(nrel_threshold(ithreshold),kind=8)
      end do
      if (iproc==0) then
          call print_errors(header=trim(header)//' (only within the sparsity pattern)')
      end if
    
      ! Full matrices
      if (check_full_matrix) then
          call set_to_zero()
          do i=1,smat%nfvctr
              do j=1,smat%nfvctr
                  tt = abs(mat%matrix(j,i,1)-mat_ref%matrix(j,i,1))
                  tt_rel = tt/abs(mat_ref%matrix(j,i,1))
                  tt_ref = mat_ref%matrix(j,i,1)
                  call calculate_errors()
              end do
          end do
          mean_error = mean_error/real(smat%nvctr,kind=8)
          mean_error_rel = mean_error_rel/real(smat%nvctr,kind=8)
          do ithreshold=1,nthreshold
              mean_error_rel_threshold(ithreshold) = mean_error_rel_threshold(ithreshold)/real(nrel_threshold(ithreshold),kind=8)
          end do
          if (iproc==0) then
              call print_errors(header=trim(header)//' (for the entire matrix)')
          end if
      end if
    
    
      contains
    
        subroutine set_to_zero()
          implicit none
          max_error = 0.0_mp
          mean_error = 0.0_mp
          max_error_rel = 0.0_mp
          mean_error_rel = 0.0_mp
          max_error_rel_threshold(:) = 0.0_mp
          mean_error_rel_threshold(:) = 0.0_mp
          nrel_threshold(:) = 0
        end subroutine set_to_zero
    
        subroutine calculate_errors()
          implicit none
          mean_error = mean_error + tt
          max_error = max(max_error,tt)
          mean_error_rel = mean_error_rel + tt_rel
          max_error_rel = max(max_error_rel,tt_rel)
          do ithreshold=1,nthreshold
              if (abs(tt_ref)>threshold(ithreshold)) then
                  nrel_threshold(ithreshold) = nrel_threshold(ithreshold) + 1
                  mean_error_rel_threshold(ithreshold) = mean_error_rel_threshold(ithreshold) + tt_rel
                  max_error_rel_threshold(ithreshold) = max(max_error_rel_threshold(ithreshold),tt_rel)
              end if
          end do
        end subroutine calculate_errors
    
        subroutine print_errors(header)
          implicit none
          character(len=*),intent(in) :: header
          call yaml_mapping_open(trim(header))
          call yaml_mapping_open('absolute error')
          call yaml_map('max error',max_error,fmt='(es10.3)')
          call yaml_map('mean error',mean_error,fmt='(es10.3)')
          call yaml_mapping_close()
          call yaml_mapping_open('relative error')
          call yaml_map('max error relative',max_error_rel,fmt='(es10.3)')
          call yaml_map('mean error relative',mean_error_rel,fmt='(es10.3)')
          call yaml_mapping_close()
          !call yaml_mapping_open('relative error with threshold')
          call yaml_sequence_open('relative error with threshold')
          do ithreshold=1,nthreshold
              call yaml_sequence(advance='no')
              call yaml_mapping_open(flow=.true.)
              call yaml_map('threshold value',threshold(ithreshold),fmt='(es8.1)')
              call yaml_map('max error relative',max_error_rel_threshold(ithreshold),fmt='(es10.3)')
              call yaml_map('mean error relative',mean_error_rel_threshold(ithreshold),fmt='(es10.3)')
              call yaml_mapping_close()
          end do
          call yaml_sequence_close()
          call yaml_mapping_close()
        end subroutine print_errors
    end subroutine calculate_error

end module utilities
