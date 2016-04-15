!> @file
!!  Linear version: Define Chebyshev polynomials
!! @author
!!    Copyright (C) 2012-2015 BigDFT group <br>
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


module chebyshev
  use sparsematrix_base
  implicit none

  private

  public :: chebyshev_clean
  public :: chebyshev_fast


  contains
 
    !> Again assuming all matrices have same sparsity, still some tidying to be done
    subroutine chebyshev_clean(iproc, nproc, npl, cc, kernel, ham_compr, &
               calculate_SHS, workarr_compr, nsize_polynomial, ncalc, &
               fermi_new, penalty_ev_new, chebyshev_polynomials, emergency_stop, &
               invovrlp_compr)
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: sequential_acces_matrix_fast, sequential_acces_matrix_fast2, &
                              compress_matrix_distributed_wrapper, sparsemm_new
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npl, nsize_polynomial, ncalc
      real(8),dimension(npl,3,ncalc),intent(in) :: cc
      type(sparse_matrix), intent(in) :: kernel
      real(kind=mp),dimension(kernel%nvctrp_tg),intent(in) :: ham_compr
      logical,intent(in) :: calculate_SHS
      real(kind=mp),dimension(kernel%nvctrp_tg),intent(inout) :: workarr_compr
      real(kind=mp),dimension(kernel%smmm%nvctrp,ncalc),intent(out) :: fermi_new
      real(kind=mp),dimension(kernel%smmm%nvctrp,2),intent(out) :: penalty_ev_new
      real(kind=mp),dimension(nsize_polynomial,npl),intent(out) :: chebyshev_polynomials
      logical,dimension(2),intent(out) :: emergency_stop
      real(kind=mp),dimension(kernel%nvctrp_tg),intent(in),optional :: invovrlp_compr
      ! Local variables
      character(len=*),parameter :: subname='chebyshev_clean'
      integer :: iorb,iiorb, jorb, ipl, i, iline, icolumn, jj, j
      integer :: isegstart, isegend, iseg, ii, jjorb, icalc
      real(kind=mp),dimension(:,:),allocatable :: vectors_new
      real(kind=mp),dimension(:),allocatable :: mat_seq
      !!real(kind=mp),dimension(:,:),allocatable :: matrix!, fermi_new, penalty_ev_new
      real(kind=mp),dimension(:),allocatable :: matrix_new
      real(kind=mp) :: tt, ddot
      integer :: jproc
    
      !call timing(iproc, 'chebyshev_comp', 'ON')
      call f_timing(TCAT_CME_POLYNOMIALS,'ON')
      call f_routine(id='chebyshev_clean')

      if (.not.kernel%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if


      !!do j=1,npl
      !!    write(*,*) 'in cheby: j, cc(j,2,1), cc(j,3,1)', j, cc(j,2,1), cc(j,3,1)
      !!end do
    
    
      !mat_compr = f_malloc(kernel%nvctrp_tg,id='mat_compr')
    

      if (kernel%nfvctrp>0) then
          mat_seq = sparsematrix_malloc(kernel, iaction=SPARSEMM_SEQ, id='mat_seq')
          vectors_new = f_malloc0((/kernel%smmm%nvctrp,4/),id='vectors_new')
      end if
        
      
      if (calculate_SHS) then
          if (.not.present(invovrlp_compr)) then
              call f_err_throw('invovrlp_compr not present')
          end if
          matrix_new = f_malloc0(kernel%smmm%nvctrp,id='matrix')
          if (kernel%smmm%nvctrp>0) then
              call prepare_matrix(kernel, invovrlp_compr, matrix_new)
              call sequential_acces_matrix_fast2(kernel, ham_compr, mat_seq)
              call sparsemm_new(iproc, kernel, mat_seq, matrix_new(1), vectors_new(1,1))
              call f_zero(matrix_new)
              call sequential_acces_matrix_fast2(kernel, invovrlp_compr, mat_seq)
              call sparsemm_new(iproc, kernel, mat_seq, vectors_new(1,1), matrix_new(1))
          end if
          call compress_matrix_distributed_wrapper(iproc, nproc, kernel, SPARSE_MATMUL_LARGE, &
               matrix_new, workarr_compr)
      else
          call vcopy(kernel%nvctrp_tg, ham_compr(1), 1, workarr_compr(1), 1)
      end if
      
      if (kernel%smmm%nvctrp>0) then
          call sequential_acces_matrix_fast2(kernel, workarr_compr, mat_seq)
      end if

      !call f_free(mat_compr)
      
        
      if (kernel%smmm%nfvctrp>0) then
        
          ! No need to set to zero the 3rd and 4th entry since they will be overwritten
          ! by copies of the 1st entry.
          if (kernel%smmm%nfvctrp>0) then
              call f_zero(2*kernel%smmm%nvctrp, vectors_new(1,1))
          end if
          do i=1,kernel%smmm%nvctrp
              ii = kernel%smmm%isvctr + i
              iline = kernel%smmm%line_and_column(1,i)
              icolumn = kernel%smmm%line_and_column(2,i)
              if (iline==icolumn) vectors_new(i,1) = 1.d0
          end do
        
          if (kernel%smmm%nvctrp>0) then

              call f_zero(fermi_new)
              call f_zero(penalty_ev_new)
        
              call vcopy(kernel%smmm%nvctrp, vectors_new(1,1), 1, vectors_new(1,3), 1)
              call vcopy(kernel%smmm%nvctrp, vectors_new(1,1), 1, vectors_new(1,4), 1)

              call compress_polynomial_vector_new(iproc, nproc, nsize_polynomial, &
                   kernel%nfvctr, kernel%smmm%nfvctrp, kernel, &
                   vectors_new(1,4), chebyshev_polynomials(1,1))

              do icalc=1,ncalc
                  call axpy(kernel%smmm%nvctrp, 0.5d0*cc(1,1,icalc), vectors_new(1,4), 1, fermi_new(1,icalc), 1)
              end do
              !write(* *) ' before loop: sum(penalty_ev_new)', sum(penalty_ev_new(:,1)), sum(penalty_ev_new(:,2))
              !write(*,*) 'cc(1,2,1), cc(1,3,1)', cc(1,2,1), cc(1,3,1)
              call axpy(kernel%smmm%nvctrp, 0.5d0*cc(1,2,1), vectors_new(1,4), 1, penalty_ev_new(1,1), 1)
              call axpy(kernel%smmm%nvctrp, 0.5d0*cc(1,3,1), vectors_new(1,4), 1, penalty_ev_new(1,2), 1)
              !write(*,*) ' before loop: sum(penalty_ev_new)', sum(penalty_ev_new(:,1)), sum(penalty_ev_new(:,2))
            
              call sparsemm_new(iproc, kernel, mat_seq, vectors_new(1,3), vectors_new(1,1))
              call vcopy(kernel%smmm%nvctrp, vectors_new(1,1), 1, vectors_new(1,2), 1)
    

              call compress_polynomial_vector_new(iproc, nproc, nsize_polynomial, &
                   kernel%nfvctr, kernel%smmm%nfvctrp, kernel, &
                   vectors_new(1,2), chebyshev_polynomials(1,2))
              do icalc=1,ncalc
                  call axpy(kernel%smmm%nvctrp, cc(2,1,icalc), vectors_new(1,2), 1, fermi_new(1,icalc), 1)
              end do
              !write(*,*) ' before loop: sum(penalty_ev_new)', sum(penalty_ev_new(:,1)), sum(penalty_ev_new(:,2))
              call axpy(kernel%smmm%nvctrp, cc(2,2,1), vectors_new(1,2), 1, penalty_ev_new(1,1), 1)
              call axpy(kernel%smmm%nvctrp, cc(2,3,1), vectors_new(1,2), 1, penalty_ev_new(1,2), 1)
            
              !write(*,*) ' before loop: sum(penalty_ev_new)', sum(penalty_ev_new(:,1)), sum(penalty_ev_new(:,2))
            
              emergency_stop=.false.
              main_loop: do ipl=3,npl
                  call sparsemm_new(iproc, kernel, mat_seq, vectors_new(1,1), vectors_new(1,2))
                  call axbyz_kernel_vectors_new(kernel, 2.d0, vectors_new(1,2), -1.d0, vectors_new(1,4), vectors_new(1,3))
                  call compress_polynomial_vector_new(iproc, nproc, nsize_polynomial, &
                       kernel%nfvctr, kernel%smmm%nfvctrp, kernel, &
                       vectors_new(1,3), chebyshev_polynomials(1,ipl))
                  do icalc=1,ncalc
                      call axpy(kernel%smmm%nvctrp, cc(ipl,1,icalc), vectors_new(1,3), 1, fermi_new(1,icalc), 1)
                  end do
                  call axpy(kernel%smmm%nvctrp, cc(ipl,2,1), vectors_new(1,3), 1, penalty_ev_new(1,1), 1)
                  call axpy(kernel%smmm%nvctrp, cc(ipl,3,1), vectors_new(1,3), 1, penalty_ev_new(1,2), 1)

                  !write(*,*) 'in loop: sum(penalty_ev_new)', &
                  !    ipl, sum(penalty_ev_new(:,1)), sum(penalty_ev_new(:,2)), &
                  !    sum(vectors_new(:,3)), sum(fermi_new(:,1)), cc(ipl,2,1)
             
                  !if (mod(ipl,2)==1) then
                  !    tt=cc(ipl,3,1)
                  !else
                  !    tt=-cc(ipl,3,1)
                  !end if
                  !call axpy(kernel%smmm%nvctrp, tt, vectors_new(1,3), 1, penalty_ev_new(1,2), 1)
             
                  call vcopy(kernel%smmm%nvctrp, vectors_new(1,1), 1, vectors_new(1,4), 1)
                  call vcopy(kernel%smmm%nvctrp, vectors_new(1,3), 1, vectors_new(1,1), 1)
    
                  ! Check the norm of the columns of the kernel and set a flag if it explodes, which might
                  ! be a consequence of the eigenvalue bounds being to small. Only
                  ! check the first matrix to be calculated.
                  !! New: Do this check on the penalty matrix
                  !!emergency_stop(1) = check_emergency_stop(kernel%smmm%nvctrp, ncalc, penalty_ev_new(1,1))
                  !!emergency_stop(2) = check_emergency_stop(kernel%smmm%nvctrp, ncalc, penalty_ev_new(1,2))
                  ! New: Do this check on the Chebyshev polynomials
                  emergency_stop(1) = check_emergency_stop(kernel%smmm%nvctrp, ncalc, vectors_new(1,3))
                  if (any(emergency_stop)) then
                      exit main_loop
                  end if
                  !!do iorb=1,kernel%smmm%nfvctrp
                  !!    !!tt=ddot(kernel%nfvctr, fermi(1,iorb,1), 1, fermi(1,iorb,1), 1)
                  !!    tt=ddot(kernel%smmm%nvctrp, fermi_new(1,1), 1, fermi_new(1,1), 1)
                  !!    if (abs(tt)>1000.d0*kernel%smmm%nvctrp) then
                  !!        emergency_stop=.true.
                  !!        exit main_loop
                  !!    end if
                  !!end do
              end do main_loop
              !write(*,*) 'emergency_stop',emergency_stop
              !write(*,*) 'sum(penalty_ev_new)', sum(penalty_ev_new(:,1)), sum(penalty_ev_new(:,2))
        
          end if
    
        
          if (calculate_SHS .and. kernel%smmm%nfvctrp>0) then
              !!call f_free(matrix)
              call f_free(matrix_new)
          end if
          if (kernel%smmm%nfvctrp>0) then
              call f_free(mat_seq)
              call f_free(vectors_new)
          end if
    
      end if
    
      !call timing(iproc, 'chebyshev_comp', 'OF')
      call f_timing(TCAT_CME_POLYNOMIALS,'OF')
      call f_release_routine()
    
    end subroutine chebyshev_clean
    

    subroutine prepare_matrix(smat, invovrlp_compr, matrix)
      use sparsematrix_init, only: matrixindex_in_compressed
      use dynamic_memory
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      real(kind=mp),dimension(smat%nvctrp_tg),intent(in) :: invovrlp_compr
      real(kind=mp),dimension(smat%smmm%nvctrp),intent(inout) :: matrix

      ! Local variables
      integer :: i, ii, iline, icolumn, jj

      call f_routine(id='prepare_matrix')

      if (.not.smat%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if

      !$omp parallel &
      !$omp default(none) &
      !$omp shared(smat, matrix, invovrlp_compr) &
      !$omp private(i, ii, iline, icolumn, jj)
      !$omp do schedule(guided)
      do i=1,smat%smmm%nvctrp
          ii = smat%smmm%isvctr + i
          iline = smat%smmm%line_and_column(1,i)
          icolumn = smat%smmm%line_and_column(2,i)
          jj=matrixindex_in_compressed(smat, icolumn, iline)
          if (jj>0) then
              matrix(i) = invovrlp_compr(jj-smat%isvctrp_tg)
          else
              matrix(i) = 0.d0
          end if
      end do
      !$omp end do
      !$omp end parallel

      call f_release_routine()

    end subroutine prepare_matrix
    

    ! Performs z = a*x + b*y
    subroutine axbyz_kernel_vectors_new(smat, a, x_compr, b, y_compr, z_compr)
      implicit none
    
      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      real(8),intent(in) :: a, b
      real(kind=mp),dimension(smat%smmm%nvctrp),intent(in) :: x_compr, y_compr
      real(kind=mp),dimension(smat%smmm%nvctrp),intent(out) :: z_compr
    
      ! Local variables
      integer :: i, jorb, iorb, ii, iline, icolumn

      call f_routine(id='axbyz_kernel_vectors_new')

      if (.not.smat%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if
    
      !$omp parallel default(shared) private(i)
      !$omp do schedule(static)
      do i=1,smat%smmm%nvctrp
          z_compr(i) = a*x_compr(i)+b*y_compr(i)
      end do
      !$omp end do
      !$omp end parallel

      call f_release_routine()
    
    end subroutine axbyz_kernel_vectors_new
    
    
    
    
    
    
    subroutine chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
               norb, norbp, fermi, chebyshev_polynomials, ncalc, cc, kernel_compressed)
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nsize_polynomial, npl, norb, norbp, ncalc
      type(sparse_matrix),intent(in) :: fermi
      real(kind=mp),dimension(nsize_polynomial,npl),intent(in) :: chebyshev_polynomials
      real(kind=mp),dimension(npl,ncalc),intent(in) :: cc
      real(kind=mp),dimension(nsize_polynomial,ncalc),intent(out) :: kernel_compressed
    
      ! Local variables
      integer :: ipl, icalc
    
      call f_routine(id='chebyshev_fast')
      !call timing(iproc, 'chebyshev_comp', 'ON')    
      call f_timing(TCAT_CME_POLYNOMIALS,'ON')
    
      if (nsize_polynomial>0) then
          call f_zero(kernel_compressed)
    
          do icalc=1,ncalc
              !write(*,*) 'icalc, ipl, kernel_compressed(1,icalc)', icalc, 0, kernel_compressed(1,icalc)
              call axpy(nsize_polynomial, 0.5d0*cc(1,icalc), chebyshev_polynomials(1,1), 1, kernel_compressed(1,icalc), 1)
              !write(*,*) 'icalc, ipl, kernel_compressed(1,icalc)', icalc, 1, kernel_compressed(1,icalc), cc(1,icalc)
              do ipl=2,npl
                  !write(*,*) 'icalc, ipl, kernel_compressed(1,icalc)', icalc, ipl, kernel_compressed(1,icalc), cc(ipl,icalc)
                  call axpy(nsize_polynomial, cc(ipl,icalc), chebyshev_polynomials(1,ipl), 1, kernel_compressed(1,icalc), 1)
              end do
          end do
    
      end if
    
      !call timing(iproc, 'chebyshev_comp', 'OF')    
      call f_timing(TCAT_CME_POLYNOMIALS,'OF')
      call f_release_routine()

    end subroutine chebyshev_fast


    subroutine compress_polynomial_vector_new(iproc, nproc, nsize_polynomial, norb, norbp, &
               fermi, vector_compr, vector_compressed)
      use sparsematrix, only: transform_sparsity_pattern
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nsize_polynomial, norb, norbp
      type(sparse_matrix),intent(in) :: fermi
      real(kind=mp),dimension(fermi%smmm%nvctrp),intent(inout) :: vector_compr
      real(kind=mp),dimension(nsize_polynomial),intent(out) :: vector_compressed
    
      ! Local variables
      integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb, iel, i, iline, icolumn
      real(kind=mp),dimension(:,:),allocatable :: vector
    
      call f_routine(id='compress_polynomial_vector_new')

      if (.not.fermi%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if
    
      call transform_sparsity_pattern(iproc, fermi%nfvctr, fermi%smmm%nvctrp_mm, fermi%smmm%isvctr_mm, &
           fermi%nseg, fermi%keyv, fermi%keyg, fermi%smmm%line_and_column_mm, &
           fermi%smmm%nvctrp, fermi%smmm%isvctr, fermi%smmm%nseg, fermi%smmm%keyv, fermi%smmm%keyg, &
           fermi%smmm%istsegline, 'large_to_small', &
           matrix_s_out=vector_compressed, matrix_l_in=vector_compr)
    
      call f_release_routine()
    
    end subroutine compress_polynomial_vector_new


    function check_emergency_stop(nvctrp, ncalc, column) result(ces)
      implicit none

      ! Calling arguments
      integer,intent(in) :: nvctrp, ncalc
      real(kind=mp),dimension(nvctrp),intent(in) :: column
      logical :: ces

      ! Local variables
      integer :: i
      real(kind=mp) :: tt

      call f_routine(id='check_emergency_stop')

      ces = .false.
      do i=1,nvctrp
          !if (abs(column(i))>1.d4) then
          if (abs(column(i))>1.d8) then
              ces = .true.
          end if
          !!write(*,*) 'sum(column(:,icalc))',sum(column(:,icalc))
          !!tt = dot(nvctrp, column(1,icalc), 1, column(1,icalc), 1)
          !!write(*,*) 'tt',tt
          !!if (abs(tt)>100000.d0*real(nvctrp,kind=mp)) then
          !!    ces = .true.
          !!end if
      end do

      call f_release_routine()

    end function check_emergency_stop


end module chebyshev
