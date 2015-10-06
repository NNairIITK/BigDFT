!> @file
!! Constrained DFT (based on linear version)
!! @author
!!    Copyright (C) 2007-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module to perform constrained DFT calculations
module constrained_dft
  use module_base
  use sparsematrix_base, only: sparse_matrix, matrices
  use dynamic_memory
  implicit none

  private


  type, public :: cdft_data
     real(wp), dimension(:), pointer :: weight_function ! the weight function defining the constraint
     type(sparse_matrix) :: weight_matrix ! matrix elements of the weight function between tmbs
     type(matrices) :: weight_matrix_ ! matrix elements of the weight function between tmbs
     integer :: ndim_dens ! the dimension of the weight function
     real(gp) :: charge ! defines the value of the charge which is to be constrained
     real(gp) :: lag_mult ! the Lagrange multiplier used to enforce the constraint
     character(len=100) :: method
     integer, dimension(2) :: ifrag_charged ! make it allocatable eventually to allow for more charged fragments
     integer :: nfrag_charged
  end type cdft_data

  public :: nullify_cdft_data, cdft_data_allocate, cdft_data_free, cdft_data_init

contains


  subroutine nullify_cdft_data(cdft)
    use sparsematrix_base, only: sparse_matrix_null, matrices_null
    implicit none
    type(cdft_data), intent(out) :: cdft
    cdft%charge=0
    cdft%lag_mult=0.0_gp
    cdft%ndim_dens=0
    nullify(cdft%weight_function)
    !call nullify_sparse_matrix(cdft%weight_matrix)
    cdft%weight_matrix = sparse_matrix_null()
    cdft%weight_matrix_ = matrices_null()
  end subroutine nullify_cdft_data

  subroutine cdft_data_free(cdft)
    use sparsematrix_base, only: deallocate_sparse_matrix, deallocate_matrices
    implicit none
    type(cdft_data), intent(inout) :: cdft

    character(len=200), parameter :: subname='cdft_data_free'

    !if (associated(cdft%weight_function)) call f_free_ptr(cdft%weight_function)
    call deallocate_sparse_matrix(cdft%weight_matrix)
    call deallocate_matrices(cdft%weight_matrix_)
    call nullify_cdft_data(cdft)
  end subroutine cdft_data_free

  subroutine cdft_data_allocate(cdft,ham)
    use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc_ptr, &
         SPARSE_FULL, assignment(=),copy_sparse_matrix
    implicit none
    type(cdft_data), intent(inout) :: cdft
    type(sparse_matrix), intent(in) :: ham

    character(len=200), parameter :: subname='cdft_data_allocate'
    !!integer :: istat

    call f_routine(id='cdft_data_allocate')
    !call sparse_copy_pattern(ham, cdft%weight_matrix, bigdft_mpi%iproc, subname)
    call copy_sparse_matrix(ham, cdft%weight_matrix)
    !cdft%weight_matrix_%matrix_compr=f_malloc_ptr(cdft%weight_matrix%nvctr,id='cdft%weight_matrix%matrix_compr')
    cdft%weight_matrix_%matrix_compr=sparsematrix_malloc_ptr(cdft%weight_matrix,iaction=SPARSE_FULL, &
                                                             id='cdft%weight_matrix%matrix_compr')
    call f_release_routine()

  end subroutine cdft_data_allocate

  subroutine cdft_data_init(cdft,input_frag,ndimrho,transfer_int)
   use fragment_base, only: fragmentInputParameters
    implicit none
    type(cdft_data), intent(inout) :: cdft
    type(fragmentInputParameters), intent(in) :: input_frag
    integer, intent(in) :: ndimrho
    logical, intent(in) :: transfer_int

    integer :: ifrag, icharged

    ! For non-transfer integral calculation only one fragment should be charged
    ! For transfer integral calculation two should have charge
    ! the value is interpreted as the charge difference and so both should have the same charge
    ! we therefore do a calculation with a +ve difference followed by a -ve difference
    cdft%nfrag_charged=0
    do ifrag=1,input_frag%nfrag
       if (input_frag%charge(ifrag)/=0) cdft%nfrag_charged=cdft%nfrag_charged+1
    end do

    if (transfer_int) then
       if (cdft%nfrag_charged/=2) then
          call f_err_throw(&
               'Error in constrained DFT, two fragments must be charged for transfer integral calculation')
          return
       end if
    else ! could generalize this later (by summing charges and fragments), but for now keep as simplest scenario
       if (cdft%nfrag_charged/=1) then
          call f_err_throw(&
               'Error in constrained DFT, exactly one fragment must have a non-zero charge value'//&
               ' unless this is a transfer integral calculation')
          return
       end if
    end if

    icharged=1
    cdft%ifrag_charged=0
    do ifrag=1,input_frag%nfrag
       if (input_frag%charge(ifrag)/=0) then
           cdft%ifrag_charged(icharged)=ifrag
           icharged=icharged+1
       end if
    end do

    if (cdft%nfrag_charged==2) then
       if (input_frag%charge(cdft%ifrag_charged(1))/=input_frag%charge(cdft%ifrag_charged(2))) then
          call f_err_throw('Error in constrained DFT, both fragments should have the same charge, '//& 
               'which is interpreted as the charge difference between the two')
          return
       end if
    end if

    cdft%charge=input_frag%charge(cdft%ifrag_charged(1))

    cdft%ndim_dens=ndimrho ! either size of fragment psi (add to fragment structure?), or size of entire simulation cell

    if (cdft%charge<0) then
       cdft%lag_mult=-0.05_gp ! pick some sensible initial value here
    else
       cdft%lag_mult=0.05_gp ! pick some sensible initial value here
    end if

    cdft%method='lowdin'
    !cdft%method='fragment_density'

  end subroutine cdft_data_init




end module constrained_dft
