!> @file
!!  Linear version: Define Chebyshev polynomials
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


module chebyshev

  implicit none

  private

  public :: chebyshev_clean
  public :: chebyshev_fast


  contains
 
    !> Again assuming all matrices have same sparsity, still some tidying to be done
    subroutine chebyshev_clean(iproc, nproc, npl, cc, kernel, ham_compr, &
               invovrlp_compr, calculate_SHS, nsize_polynomial, ncalc, fermi, penalty_ev, chebyshev_polynomials, &
               emergency_stop)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc, assignment(=), &
                                   DENSE_MATMUL, SPARSEMM_SEQ,sparsematrix_malloc0
      use sparsematrix, only: sequential_acces_matrix_fast, sequential_acces_matrix_fast2, &
                              sparsemm, compress_matrix_distributed
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npl, nsize_polynomial, ncalc
      real(8),dimension(npl,3,ncalc),intent(in) :: cc
      type(sparse_matrix), intent(in) :: kernel
      real(kind=8),dimension(kernel%nvctrp_tg),intent(in) :: ham_compr
      real(kind=8),dimension(kernel%nvctrp_tg),intent(in) :: invovrlp_compr
      logical,intent(in) :: calculate_SHS
      real(kind=8),dimension(kernel%nfvctr,kernel%smmm%nfvctrp,ncalc),intent(out) :: fermi
      real(kind=8),dimension(kernel%nfvctr,kernel%smmm%nfvctrp,2),intent(out) :: penalty_ev
      real(kind=8),dimension(nsize_polynomial,npl),intent(out) :: chebyshev_polynomials
      logical,intent(out) :: emergency_stop
      ! Local variables
      integer :: iorb,iiorb, jorb, ipl, ierr, nseq, nmaxvalk
      integer :: isegstart, isegend, iseg, ii, jjorb, icalc
      character(len=*),parameter :: subname='chebyshev_clean'
      real(8), dimension(:,:,:), allocatable :: vectors
      real(kind=8),dimension(:),allocatable :: mat_seq, mat_compr
      real(kind=8),dimension(:,:),allocatable :: matrix
      real(kind=8) :: tt, ddot
    
      call timing(iproc, 'chebyshev_comp', 'ON')
      call f_routine(id='chebyshev_clean')
    
      !!kernel%nfvctr = kernel%nfvctr
      !!kernel%nfvctrp = kernel%nfvctrp
      !!kernel%isfvctr = kernel%isfvctr
    
      mat_compr = f_malloc(kernel%nvctrp_tg,id='mat_compr')
    
      if (calculate_SHS) then
          matrix = sparsematrix_malloc0(kernel, iaction=DENSE_MATMUL, id='matrix')
      end if
      if (kernel%nfvctrp>0) then
    
        
        
          mat_seq = sparsematrix_malloc(kernel, iaction=SPARSEMM_SEQ, id='mat_seq')
        
          if (calculate_SHS) then
              !!matrix = sparsematrix_malloc(kernel, iaction=DENSE_MATMUL, id='matrix')
        
              !if (kernel%smmm%nfvctrp>0) then
              !    call f_zero(kernel%nfvctr*kernel%smmm%nfvctrp, matrix(1,1))
              !end if
              !write(*,*) 'WARNING CHEBYSHEV: MODIFYING MATRIX MULTIPLICATION'
              if (kernel%smmm%nfvctrp>0) then
    
                  !do iseg=isegstart,isegend
                  do iseg=kernel%smmm%isseg,kernel%smmm%ieseg
                      ii=kernel%keyv(iseg)-1
                      ! A segment is always on one line, therefore no double loop
                      do jorb=kernel%keyg(1,1,iseg),kernel%keyg(2,1,iseg)
                          ii=ii+1
                          if (ii<kernel%smmm%isvctr+1) cycle
                          if (ii>kernel%smmm%isvctr+kernel%smmm%nvctrp) exit
                          iiorb = kernel%keyg(1,2,iseg)
                          jjorb = jorb
                          matrix(jjorb,iiorb-kernel%smmm%isfvctr)=invovrlp_compr(ii-kernel%isvctrp_tg)
                          !if (jjorb==iiorb) then
                          !    matrix(jjorb,iiorb-kernel%isfvctr)=1.d0
                          !else
                          !    matrix(jjorb,iiorb-kernel%isfvctr)=0.d0
                          !end if
                      end do
                  end do
              end if
          end if
        
        
        
    
    
        
          vectors = f_malloc0((/ kernel%nfvctr, kernel%smmm%nfvctrp, 4 /),id='vectors')
          !if (kernel%smmm%nfvctrp>0) then
          !    call f_zero(kernel%nfvctr*kernel%smmm%nfvctrp, vectors(1,1,1))
          !end if
    
        
      end if
        
      
      if (calculate_SHS) then
      
          if (kernel%smmm%nfvctrp>0) then
              call sequential_acces_matrix_fast2(kernel, ham_compr, mat_seq)
              call sparsemm(kernel, mat_seq, matrix(1,1), vectors(1,1,1))
              call f_zero(matrix)
              ! use mat_seq as workarray
              call sequential_acces_matrix_fast2(kernel, invovrlp_compr, mat_seq)
              call sparsemm(kernel, mat_seq, vectors(1,1,1), matrix(1,1))
              !call f_zero(kernel%nvctr, SHS(1))
          end if
          !call f_zero(kernel%nvctrp_tg, SHS(1))
          
          if (kernel%smmm%nfvctrp>0) then
              do iseg=kernel%smmm%isseg,kernel%smmm%ieseg
                  ii=kernel%keyv(iseg)-1
                  ! A segment is always on one line, therefore no double loop
                  do jorb=kernel%keyg(1,1,iseg),kernel%keyg(2,1,iseg)
                      ii=ii+1
                      if (ii<kernel%smmm%isvctr+1) cycle
                      if (ii>kernel%smmm%isvctr+kernel%smmm%nvctrp) exit
                      iiorb = kernel%keyg(1,2,iseg)
                      jjorb = jorb
                      mat_compr(ii-kernel%isvctrp_tg)=matrix(jjorb,iiorb-kernel%smmm%isfvctr)
                  end do
              end do
          end if
      
          if (nproc > 1) then
             !call mpiallred(SHS(1), kernel%nvctr, mpi_sum, bigdft_mpi%mpi_comm)
             call compress_matrix_distributed(iproc, nproc, kernel, DENSE_MATMUL, &
                  matrix, mat_compr)
          end if
          
      else
          call vcopy(kernel%nvctrp_tg, ham_compr(1), 1, mat_compr(1), 1)
      
      end if
      
      if (kernel%smmm%nfvctrp>0) then
          call sequential_acces_matrix_fast2(kernel, mat_compr, mat_seq)
      end if
      
    
      !!if (iproc==0) then
      !!    do istat=1,kernel%nvctr
      !!        write(300,*) ham_compr(istat), SHS(istat)
      !!    end do
      !!end if
        
      if (kernel%smmm%nfvctrp>0) then
        
          ! No need to set to zero the 3rd and 4th entry since they will be overwritten
          ! by copies of the 1st entry.
          if (kernel%smmm%nfvctrp>0) then
              call f_zero(2*kernel%nfvctr*kernel%smmm%nfvctrp, vectors(1,1,1))
          end if
          do iorb=1,kernel%smmm%nfvctrp
              iiorb=kernel%smmm%isfvctr+iorb
              vectors(iiorb,iorb,1)=1.d0
          end do
        
          if (kernel%smmm%nfvctrp>0) then
        
              call vcopy(kernel%nfvctr*kernel%smmm%nfvctrp, vectors(1,1,1), 1, vectors(1,1,3), 1)
              call vcopy(kernel%nfvctr*kernel%smmm%nfvctrp, vectors(1,1,1), 1, vectors(1,1,4), 1)
            
              ! apply(3/2 - 1/2 S) H (3/2 - 1/2 S)
              !!if (number_of_matmuls==three) then
              !!    call sparsemm(kernel, invovrlp_compr_seq, vectors(1,1,3), vectors(1,1,1))
              !!    call sparsemm(kernel, ham_compr_seq, vectors(1,1,1), vectors(1,1,3))
              !!    call sparsemm(kernel, invovrlp_compr_seq, vectors(1,1,3), vectors(1,1,1))
              !!else if (number_of_matmuls==one) then
                  call sparsemm(kernel, mat_seq, vectors(1,1,3), vectors(1,1,1))
              !!end if
            
            
              call vcopy(kernel%nfvctr*kernel%smmm%nfvctrp, vectors(1,1,1), 1, vectors(1,1,2), 1)
    
              !initialize fermi
              call f_zero(fermi)
              call f_zero(penalty_ev)
              call compress_polynomial_vector(iproc, nproc, nsize_polynomial, &
                   kernel%nfvctr, kernel%smmm%nfvctrp, kernel%smmm%isfvctr, kernel, &
                   vectors(1,1,4), chebyshev_polynomials(1,1))
              do icalc=1,ncalc
                  call axpy_kernel_vectors(kernel%smmm%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                       0.5d0*cc(1,1,icalc), vectors(1,1,4), fermi(:,1,icalc))
              end do
              call axpy_kernel_vectors(kernel%smmm%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                   0.5d0*cc(1,3,1), vectors(1,1,4), penalty_ev(:,1,1))
              call axpy_kernel_vectors(kernel%smmm%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                   0.5d0*cc(1,3,1), vectors(1,1,4), penalty_ev(:,1,2))
              call compress_polynomial_vector(iproc, nproc, nsize_polynomial, &
                   kernel%nfvctr, kernel%smmm%nfvctrp, kernel%smmm%isfvctr, kernel, &
                   vectors(1,1,2), chebyshev_polynomials(1,2))
              do icalc=1,ncalc
                  call axpy_kernel_vectors(kernel%smmm%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                       cc(2,1,icalc), vectors(1,1,2), fermi(:,1,icalc))
              end do
              call axpy_kernel_vectors(kernel%smmm%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                   cc(2,3,1), vectors(1,1,2), penalty_ev(:,1,1))
              call axpy_kernel_vectors(kernel%smmm%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                   -cc(2,3,1), vectors(1,1,2), penalty_ev(:,1,2))
            
            
              emergency_stop=.false.
              main_loop: do ipl=3,npl
                  ! apply (3/2 - 1/2 S) H (3/2 - 1/2 S)
                  !!if (number_of_matmuls==three) then
                  !!    call sparsemm(kernel, invovrlp_compr_seq, vectors(1,1,1), vectors(1,1,2))
                  !!    call sparsemm(kernel, ham_compr_seq, vectors(1,1,2), vectors(1,1,3))
                  !!    call sparsemm(kernel, invovrlp_compr_seq, vectors(1,1,3), vectors(1,1,2))
                  !!else if (number_of_matmuls==one) then
                      call sparsemm(kernel, mat_seq, vectors(1,1,1), vectors(1,1,2))
                  !!end if
                  call axbyz_kernel_vectors(kernel%smmm%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                       2.d0, vectors(1,1,2), -1.d0, vectors(1,1,4), vectors(1,1,3))
                  call compress_polynomial_vector(iproc, nproc, nsize_polynomial, &
                       kernel%nfvctr, kernel%smmm%nfvctrp, kernel%smmm%isfvctr, kernel, vectors(1,1,3), &
                       chebyshev_polynomials(1,ipl))
                  do icalc=1,ncalc
                      call axpy_kernel_vectors(kernel%smmm%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                           cc(ipl,1,icalc), vectors(1,1,3), fermi(:,1,icalc))
                  end do
                  call axpy_kernel_vectors(kernel%smmm%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                       cc(ipl,3,1), vectors(1,1,3), penalty_ev(:,1,1))
             
                  if (mod(ipl,2)==1) then
                      tt=cc(ipl,3,1)
                  else
                      tt=-cc(ipl,3,1)
                  end if
                  call axpy_kernel_vectors(kernel%smmm%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                       tt, vectors(1,1,3), penalty_ev(:,1,2))
             
                  call copy_kernel_vectors(kernel%smmm%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                       vectors(1,1,1), vectors(1,1,4))
                  call copy_kernel_vectors(kernel%smmm%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                       vectors(1,1,3), vectors(1,1,1))
    
                  ! Check the norm of the columns of the kernel and set a flag if it explodes, which might
                  ! be a consequence of the eigenvalue bounds being to small. Only
                  ! check the first matrix to be calculated.
                  do iorb=1,kernel%smmm%nfvctrp
                      tt=ddot(kernel%nfvctr, fermi(1,iorb,1), 1, fermi(1,iorb,1), 1)
                      if (abs(tt)>1.d3) then
                          emergency_stop=.true.
                          exit main_loop
                      end if
                  end do
              end do main_loop
        
          end if
    
     
        
          if (calculate_SHS .and. kernel%smmm%nfvctrp>0) then
              call f_free(matrix)
          end if
          if (kernel%smmm%nfvctrp>0) then
              call f_free(mat_seq)
              call f_free(vectors)
          end if
          call f_free(mat_compr)
    
      end if
    
      call timing(iproc, 'chebyshev_comp', 'OF')
      call f_release_routine()
    
    end subroutine chebyshev_clean
    
    
    
    ! Performs z = a*x + b*y
    subroutine axbyz_kernel_vectors(norbp, norb, nout, onedimindices, a, x, b, y, z)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: norbp, norb, nout
      integer,dimension(4,nout),intent(in) :: onedimindices
      real(8),intent(in) :: a, b
      real(kind=8),dimension(norb,norbp),intent(in) :: x, y
      real(kind=8),dimension(norb,norbp),intent(out) :: z
    
      ! Local variables
      integer :: i, jorb, iorb

      call f_routine(id='axbyz_kernel_vectors')
    
      !$omp parallel default(private) shared(nout, onedimindices,a, b, x, y, z)
      !$omp do
      do i=1,nout
          iorb=onedimindices(1,i)
          jorb=onedimindices(2,i)
          z(jorb,iorb)=a*x(jorb,iorb)+b*y(jorb,iorb)
      end do
      !$omp end do
      !$omp end parallel

      call f_release_routine()
    
    end subroutine axbyz_kernel_vectors
    
    
    
    
    subroutine copy_kernel_vectors(norbp, norb, nout, onedimindices, a, b)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: norbp, norb, nout
      integer,dimension(4,nout),intent(in) :: onedimindices
      real(kind=8),dimension(norb,norbp),intent(in) :: a
      real(kind=8),dimension(norb,norbp),intent(out) :: b
    
      ! Local variables
      integer :: i, jorb, iorb
    
      call f_routine(id='copy_kernel_vectors')
    
      !$omp parallel default(private) shared(nout, onedimindices,a, b)
      !$omp do
      do i=1,nout
          iorb=onedimindices(1,i)
          jorb=onedimindices(2,i)
          b(jorb,iorb)=a(jorb,iorb)
      end do
      !$omp end do
      !$omp end parallel
    
      call f_release_routine()
    
    end subroutine copy_kernel_vectors
    
    
    
    
    subroutine axpy_kernel_vectors(norbp, norb, nout, onedimindices, a, x, y)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: norbp, norb, nout
      integer,dimension(4,nout),intent(in) :: onedimindices
      real(kind=8),intent(in) :: a
      real(kind=8),dimension(norb,norbp),intent(in) :: x
      real(kind=8),dimension(norb,norbp),intent(inout) :: y
    
      ! Local variables
      integer :: i, jorb, iorb

      call f_routine(id='axpy_kernel_vectors')
    
      !$omp parallel default(private) shared(nout, onedimindices, y, x, a)
      !$omp do
      do i=1,nout
          iorb=onedimindices(1,i)
          jorb=onedimindices(2,i)
          y(jorb,iorb)=y(jorb,iorb)+a*x(jorb,iorb)
      end do
      !$omp end do
      !$omp end parallel
    
      call f_release_routine()
    
    end subroutine axpy_kernel_vectors
    
    
    
    subroutine chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
               norb, norbp, fermi, chebyshev_polynomials, ncalc, cc, kernelp)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc, assignment(=),&
           SPARSE_FULL,sparsematrix_malloc0
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nsize_polynomial, npl, norb, norbp, ncalc
      type(sparse_matrix),intent(in) :: fermi
      real(kind=8),dimension(nsize_polynomial,npl),intent(in) :: chebyshev_polynomials
      real(kind=8),dimension(npl,ncalc),intent(in) :: cc
      real(kind=8),dimension(norb,norbp,ncalc),intent(out) :: kernelp
    
      ! Local variables
      integer :: ipl, icalc
      real(kind=8),dimension(:),allocatable :: kernel_compressed
    
      call f_routine(id='chebyshev_fast')
    
      if (nsize_polynomial>0) then
          kernel_compressed = sparsematrix_malloc0(fermi, iaction=SPARSE_FULL, id='kernel_compressed')
    
          !call f_zero(nsize_polynomial,kernel_compressed(1))
          do icalc=1,ncalc
              call daxpy(nsize_polynomial, 0.5d0*cc(1,icalc), chebyshev_polynomials(1,1), 1, kernel_compressed(1), 1)
              do ipl=2,npl
                  call daxpy(nsize_polynomial, cc(ipl,icalc), chebyshev_polynomials(1,ipl), 1, kernel_compressed(1), 1)
              end do
              call uncompress_polynomial_vector(iproc, nproc, nsize_polynomial, &
                   fermi, kernel_compressed, kernelp(1,1,icalc))
          end do
    
    
          call f_free(kernel_compressed)
      end if
    
      call f_release_routine()
    
    end subroutine chebyshev_fast

end module chebyshev
