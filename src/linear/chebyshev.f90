!> @file
!!  Linear version: Define Chebyshev polynomials
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 
!> Again assuming all matrices have same sparsity, still some tidying to be done
subroutine chebyshev_clean(iproc, nproc, npl, cc, norb, norbp, isorb, isorb_par, foe_obj, kernel, ham_compr, &
           ovrlp_compr, calculate_SHS, nsize_polynomial, SHS, fermi, penalty_ev, chebyshev_polynomials, &
           emergency_stop)
  use module_base
  use module_types
  use module_interfaces, except_this_one => chebyshev_clean
  use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc, assignment(=), &
                               DENSE_PARALLEL, SPARSEMM_SEQ
  use sparsematrix, only: sequential_acces_matrix_fast, sparsemm
  use foe_base, only: foe_data
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, npl, nsize_polynomial, norb, norbp, isorb
  integer,dimension(0:nproc-1),intent(in) :: isorb_par
  real(8),dimension(npl,3),intent(in) :: cc
  type(foe_data),intent(in) :: foe_obj
  type(sparse_matrix), intent(in) :: kernel
  real(kind=8),dimension(kernel%nvctr),intent(in) :: ham_compr, ovrlp_compr
  logical,intent(in) :: calculate_SHS
  real(kind=8),dimension(kernel%nvctr),intent(inout) :: SHS
  real(kind=8),dimension(kernel%nfvctr,kernel%nfvctrp),intent(out) :: fermi
  real(kind=8),dimension(kernel%nfvctr,kernel%nfvctrp,2),intent(out) :: penalty_ev
  real(kind=8),dimension(nsize_polynomial,npl),intent(out) :: chebyshev_polynomials
  logical,intent(out) :: emergency_stop
  ! Local variables
  integer :: iorb,iiorb, jorb, ipl, ierr, nseq, nmaxsegk, nmaxvalk
  integer :: isegstart, isegend, iseg, ii, jjorb, nout
  character(len=*),parameter :: subname='chebyshev_clean'
  real(8), dimension(:,:,:), allocatable :: vectors
  real(kind=8),dimension(:),allocatable :: ham_compr_seq, ovrlp_compr_seq, SHS_seq
  real(kind=8),dimension(:,:),allocatable :: matrix
  real(kind=8) :: tt, ddot
  integer,dimension(:,:,:),allocatable :: istindexarr
  integer,dimension(:),allocatable :: ivectorindex
  integer,parameter :: one=1, three=3
  integer,parameter :: number_of_matmuls=one
  integer,dimension(:,:),pointer :: onedimindices

  call timing(iproc, 'chebyshev_comp', 'ON')
  call f_routine(id='chebyshev_clean')

  !!kernel%nfvctr = kernel%nfvctr
  !!kernel%nfvctrp = kernel%nfvctrp
  !!kernel%isfvctr = kernel%isfvctr

  if (kernel%nfvctrp>0) then

    
      ham_compr_seq = sparsematrix_malloc(kernel, iaction=SPARSEMM_SEQ, id='ham_compr_seq')
      ovrlp_compr_seq = sparsematrix_malloc(kernel, iaction=SPARSEMM_SEQ, id='ovrlp_compr_seq')
    
    
      if (number_of_matmuls==one) then
          matrix = sparsematrix_malloc(kernel, iaction=DENSE_PARALLEL, id='matrix')
          SHS_seq = sparsematrix_malloc(kernel, iaction=SPARSEMM_SEQ, id='SHS_seq')
    
          if (kernel%nfvctrp>0) then
              call to_zero(kernel%nfvctr*kernel%nfvctrp, matrix(1,1))
          end if
          !write(*,*) 'WARNING CHEBYSHEV: MODIFYING MATRIX MULTIPLICATION'
          if (kernel%nfvctrp>0) then
              isegstart=kernel%istsegline(kernel%isfvctr+1)
              if (kernel%isfvctr+kernel%nfvctrp<kernel%nfvctr) then
                  isegend=kernel%istsegline(kernel%isfvctr_par(iproc+1)+1)-1
              else
                  isegend=kernel%nseg
              end if
              do iseg=isegstart,isegend
                  ii=kernel%keyv(iseg)-1
                  do jorb=kernel%keyg(1,iseg),kernel%keyg(2,iseg)
                      ii=ii+1
                      iiorb = (jorb-1)/kernel%nfvctr + 1
                      jjorb = jorb - (iiorb-1)*kernel%nfvctr
                      matrix(jjorb,iiorb-kernel%isfvctr)=ovrlp_compr(ii)
                      !if (jjorb==iiorb) then
                      !    matrix(jjorb,iiorb-kernel%isfvctr)=1.d0
                      !else
                      !    matrix(jjorb,iiorb-kernel%isfvctr)=0.d0
                      !end if
                  end do
              end do
          end if
      end if
    
      !!call sequential_acces_matrix(kernel%nfvctr, kernel%nfvctrp, kernel%isfvctr, kernel%smmm%nseg, &
      !!     kernel%smmm%nsegline, kernel%smmm%istsegline, kernel%smmm%keyg, &
      !!     kernel, ham_compr, kernel%smmm%nseq, kernel%smmm%nmaxsegk, kernel%smmm%nmaxvalk, &
      !!     ham_compr_seq)
      call sequential_acces_matrix_fast(kernel, ham_compr, ham_compr_seq)
    
    
      !!call sequential_acces_matrix(kernel%nfvctr, kernel%nfvctrp, kernel%isfvctr, kernel%smmm%nseg, &
      !!     kernel%smmm%nsegline, kernel%smmm%istsegline, kernel%smmm%keyg, &
      !!     kernel, ovrlp_compr, kernel%smmm%nseq, kernel%smmm%nmaxsegk, kernel%smmm%nmaxvalk, &
      !!     ovrlp_compr_seq)
      call sequential_acces_matrix_fast(kernel, ovrlp_compr, ovrlp_compr_seq)


    
      vectors = f_malloc((/ kernel%nfvctr, kernel%nfvctrp, 4 /),id='vectors')
      if (kernel%nfvctrp>0) then
          call to_zero(kernel%nfvctr*kernel%nfvctrp, vectors(1,1,1))
      end if
    
  end if
    
  if (number_of_matmuls==one) then
  
      if (calculate_SHS) then
  
          if (kernel%nfvctrp>0) then
              call sparsemm(kernel, ham_compr_seq, matrix(1,1), vectors(1,1,1))
              call to_zero(kernel%nfvctrp*kernel%nfvctr, matrix(1,1))
              call sparsemm(kernel, ovrlp_compr_seq, vectors(1,1,1), matrix(1,1))
              !call to_zero(kernel%nvctr, SHS(1))
          end if
          call to_zero(kernel%nvctr, SHS(1))
          
          if (kernel%nfvctrp>0) then
              isegstart=kernel%istsegline(kernel%isfvctr+1)
              if (kernel%isfvctr+kernel%nfvctrp<kernel%nfvctr) then
                  isegend=kernel%istsegline(kernel%isfvctr_par(iproc+1)+1)-1
              else
                  isegend=kernel%nseg
              end if
              do iseg=isegstart,isegend
                  ii=kernel%keyv(iseg)-1
                  do jorb=kernel%keyg(1,iseg),kernel%keyg(2,iseg)
                      ii=ii+1
                      iiorb = (jorb-1)/kernel%nfvctr + 1
                      jjorb = jorb - (iiorb-1)*kernel%nfvctr
                      SHS(ii)=matrix(jjorb,iiorb-kernel%isfvctr)
                  end do
              end do
          end if
  
          if (nproc > 1) then
             call mpiallred(SHS(1), kernel%nvctr, mpi_sum, bigdft_mpi%mpi_comm)
          end if

      else
          ! This is quick and dirty...
          SHS = ham_compr
  
      end if
  
      if (kernel%nfvctrp>0) then
          !!call sequential_acces_matrix(kernel%nfvctr, kernel%nfvctrp, kernel%isfvctr, kernel%smmm%nseg, &
          !!     kernel%smmm%nsegline, kernel%smmm%istsegline, kernel%smmm%keyg, &
          !!     kernel, SHS, kernel%smmm%nseq, kernel%smmm%nmaxsegk, &
          !!     kernel%smmm%nmaxvalk, SHS_seq)
          call sequential_acces_matrix_fast(kernel, SHS, SHS_seq)
      end if
  
  end if

  !!if (iproc==0) then
  !!    do istat=1,kernel%nvctr
  !!        write(300,*) ham_compr(istat), SHS(istat)
  !!    end do
  !!end if
    
  if (kernel%nfvctrp>0) then
    
      ! No need to set to zero the 3rd and 4th entry since they will be overwritten
      ! by copies of the 1st entry.
      if (kernel%nfvctrp>0) then
          call to_zero(2*kernel%nfvctr*kernel%nfvctrp, vectors(1,1,1))
      end if
      do iorb=1,kernel%nfvctrp
          iiorb=kernel%isfvctr+iorb
          vectors(iiorb,iorb,1)=1.d0
      end do
    
      if (kernel%nfvctrp>0) then
    
          call vcopy(kernel%nfvctr*kernel%nfvctrp, vectors(1,1,1), 1, vectors(1,1,3), 1)
          call vcopy(kernel%nfvctr*kernel%nfvctrp, vectors(1,1,1), 1, vectors(1,1,4), 1)
        
          ! apply(3/2 - 1/2 S) H (3/2 - 1/2 S)
          if (number_of_matmuls==three) then
              call sparsemm(kernel, ovrlp_compr_seq, vectors(1,1,3), vectors(1,1,1))
              call sparsemm(kernel, ham_compr_seq, vectors(1,1,1), vectors(1,1,3))
              call sparsemm(kernel, ovrlp_compr_seq, vectors(1,1,3), vectors(1,1,1))
          else if (number_of_matmuls==one) then
              call sparsemm(kernel, SHS_seq, vectors(1,1,3), vectors(1,1,1))
          end if
        
        
          call vcopy(kernel%nfvctr*kernel%nfvctrp, vectors(1,1,1), 1, vectors(1,1,2), 1)
        
          !initialize fermi
          call to_zero(kernel%nfvctrp*kernel%nfvctr, fermi(1,1))
          call to_zero(2*kernel%nfvctr*kernel%nfvctrp, penalty_ev(1,1,1))
          call compress_polynomial_vector(iproc, nproc, nsize_polynomial, &
               kernel%nfvctr, kernel%nfvctrp, kernel%isfvctr, isorb_par, kernel, &
               vectors(1,1,4), chebyshev_polynomials(1,1))
          call axpy_kernel_vectors(kernel%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
               0.5d0*cc(1,1), vectors(1,1,4), fermi(:,1))
          call axpy_kernel_vectors(kernel%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
               0.5d0*cc(1,3), vectors(1,1,4), penalty_ev(:,1,1))
          call axpy_kernel_vectors(kernel%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
               0.5d0*cc(1,3), vectors(1,1,4), penalty_ev(:,1,2))
          call compress_polynomial_vector(iproc, nproc, nsize_polynomial, &
               kernel%nfvctr, kernel%nfvctrp, kernel%isfvctr, isorb_par, kernel, vectors(1,1,2), chebyshev_polynomials(1,2))
          call axpy_kernel_vectors(kernel%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
               cc(2,1), vectors(1,1,2), fermi(:,1))
          call axpy_kernel_vectors(kernel%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
               cc(2,3), vectors(1,1,2), penalty_ev(:,1,1))
          call axpy_kernel_vectors(kernel%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
               -cc(2,3), vectors(1,1,2), penalty_ev(:,1,2))
        
        
          emergency_stop=.false.
          main_loop: do ipl=3,npl
              ! apply (3/2 - 1/2 S) H (3/2 - 1/2 S)
              if (number_of_matmuls==three) then
                  call sparsemm(kernel, ovrlp_compr_seq, vectors(1,1,1), vectors(1,1,2))
                  call sparsemm(kernel, ham_compr_seq, vectors(1,1,2), vectors(1,1,3))
                  call sparsemm(kernel, ovrlp_compr_seq, vectors(1,1,3), vectors(1,1,2))
              else if (number_of_matmuls==one) then
                  call sparsemm(kernel, SHS_seq, vectors(1,1,1), vectors(1,1,2))
              end if
              call axbyz_kernel_vectors(kernel%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                   2.d0, vectors(1,1,2), -1.d0, vectors(1,1,4), vectors(1,1,3))
              call compress_polynomial_vector(iproc, nproc, nsize_polynomial, &
                   kernel%nfvctr, kernel%nfvctrp, kernel%isfvctr, isorb_par, kernel, vectors(1,1,3), &
                   chebyshev_polynomials(1,ipl))
              call axpy_kernel_vectors(kernel%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                   cc(ipl,1), vectors(1,1,3), fermi(:,1))
              call axpy_kernel_vectors(kernel%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                   cc(ipl,3), vectors(1,1,3), penalty_ev(:,1,1))
         
              if (mod(ipl,2)==1) then
                  tt=cc(ipl,3)
              else
                  tt=-cc(ipl,3)
              end if
              call axpy_kernel_vectors(kernel%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                   tt, vectors(1,1,3), penalty_ev(:,1,2))
         
              call copy_kernel_vectors(kernel%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                   vectors(1,1,1), vectors(1,1,4))
              call copy_kernel_vectors(kernel%nfvctrp, kernel%nfvctr, kernel%smmm%nout, kernel%smmm%onedimindices, &
                   vectors(1,1,3), vectors(1,1,1))

              ! Check the norm of the columns of the kernel and set a flag if it explodes, which might
              ! be a consequence of the eigenvalue bounds being to small.
              do iorb=1,kernel%nfvctrp
                  tt=ddot(kernel%nfvctr, fermi(1,iorb), 1, fermi(1,iorb), 1)
                  if (abs(tt)>1.d3) then
                      emergency_stop=.true.
                      exit main_loop
                  end if
              end do
          end do main_loop
    
      end if

 
    
      call f_free(vectors)
      call f_free(ham_compr_seq)
      call f_free(ovrlp_compr_seq)
    
      if (number_of_matmuls==one) then
          call f_free(matrix)
          call f_free(SHS_seq)
      end if

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

  !$omp parallel default(private) shared(nout, onedimindices,a, b, x, y, z)
  !$omp do
  do i=1,nout
      iorb=onedimindices(1,i)
      jorb=onedimindices(2,i)
      z(jorb,iorb)=a*x(jorb,iorb)+b*y(jorb,iorb)
  end do
  !$omp end do
  !$omp end parallel

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


  !$omp parallel default(private) shared(nout, onedimindices,a, b)
  !$omp do
  do i=1,nout
      iorb=onedimindices(1,i)
      jorb=onedimindices(2,i)
      b(jorb,iorb)=a(jorb,iorb)
  end do
  !$omp end do
  !$omp end parallel


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

  !$omp parallel default(private) shared(nout, onedimindices, y, x, a)
  !$omp do
  do i=1,nout
      iorb=onedimindices(1,i)
      jorb=onedimindices(2,i)
      y(jorb,iorb)=y(jorb,iorb)+a*x(jorb,iorb)
  end do
  !$omp end do
  !$omp end parallel


end subroutine axpy_kernel_vectors



subroutine chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
           norb, norbp, isorb, isorb_par, fermi, chebyshev_polynomials, cc, kernelp)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc, assignment(=), SPARSE_FULL
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nsize_polynomial, npl, norb, norbp, isorb
  integer,dimension(0:nproc-1),intent(in) :: isorb_par
  type(sparse_matrix),intent(in) :: fermi
  real(kind=8),dimension(nsize_polynomial,npl),intent(in) :: chebyshev_polynomials
  real(kind=8),dimension(npl),intent(in) :: cc
  real(kind=8),dimension(norb,norbp),intent(out) :: kernelp

  ! Local variables
  integer :: ipl, iall
  real(kind=8),dimension(:),allocatable :: kernel_compressed

  call f_routine(id='chebyshev_fast')

  if (nsize_polynomial>0) then
      kernel_compressed = sparsematrix_malloc(fermi, iaction=SPARSE_FULL, id='kernel_compressed')

      call to_zero(nsize_polynomial,kernel_compressed(1))
      !write(*,*) 'ipl, first element', 1, chebyshev_polynomials(1,1)
      call daxpy(nsize_polynomial, 0.5d0*cc(1), chebyshev_polynomials(1,1), 1, kernel_compressed(1), 1)
      do ipl=2,npl
      !write(*,*) 'ipl, first element', ipl, chebyshev_polynomials(1,ipl)
          call daxpy(nsize_polynomial, cc(ipl), chebyshev_polynomials(1,ipl), 1, kernel_compressed(1), 1)
      end do

      call uncompress_polynomial_vector(iproc, nproc, nsize_polynomial, &
           norb, norbp, isorb, isorb_par, fermi, kernel_compressed, kernelp)

      call f_free(kernel_compressed)
  end if

  call f_release_routine()

end subroutine chebyshev_fast
