!> @file
!!  Define a module to wrap the linear algebra routines
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Modules which defines wrappers for the linear alegra.
module wrapper_linalg
  !!@todo MOVE ME TO F_MALLOC
  ! MOVE ME TO F_MALLOC
  !> initialize to zero an array
  interface to_zero
     module procedure put_to_zero_simple, &
          & put_to_zero_double, put_to_zero_double_1, put_to_zero_double_2, &
          & put_to_zero_integer
  end interface
  ! MOVE ME TO F_MALLOC

  !> Flag for GPU computing, if CUDA libraries are present
  !! in that case if a GPU is present a given MPI processor may or not perform a GPU calculation
  !! this value can be changed in the read_input_variables routine
  logical :: GPUblas=.false.

  !> interfaces for LAPACK routines
  interface potrf
     module procedure potrf_simple,potrf_double
  end interface
  interface c_potrf
     module procedure c_potrf_simple,c_potrf_double
  end interface
  interface trtri
     module procedure trtri_simple,trtri_double
  end interface 
  interface c_trtri
     module procedure c_trtri_simple,c_trtri_double
  end interface
  interface syev
     module procedure syev_simple,syev_double
  end interface
  interface heev
     module procedure heev_simple,heev_double
  end interface
  interface sygv
     module procedure sygv_simple,sygv_double
  end interface
  interface hegv
     module procedure hegv_simple,hegv_double
  end interface
  interface gesv
     module procedure gesv_simple,gesv_double
  end interface
  interface c_gesv
     module procedure c_gesv_simple,c_gesv_double
  end interface


  !> interfaces for BLAS routines
  interface gemm
     module procedure gemm_simple,gemm_double
  end interface
  interface gemmsy
     module procedure gemm_simple,gemmsy_double_wrap
  end interface
  interface c_gemm
     module procedure c_gemm_simple,c_gemm_double
  end interface
  interface dot
     module procedure dot_simple,dot_double
  end interface
  interface dotc
     module procedure dotc_simple,dotc_double
  end interface
  interface nrm2
     module procedure nrm2_simple,nrm2_double
  end interface
  interface vscal
     module procedure scal_simple,scal_double
  end interface
  interface vcopy
     module procedure copy_integer,copy_simple,copy_double,copy_double_to_simple,&
          copy_complex_real_simple,copy_complex_real_double
  end interface vcopy
  interface c_vscal
     module procedure c_scal_simple,c_scal_double
  end interface c_vscal
  interface syrk
     module procedure syrk_simple,syrk_double
  end interface syrk
  interface herk
     module procedure herk_simple,herk_double
  end interface herk
  interface trmm
     module procedure trmm_simple,trmm_double
  end interface trmm
  interface c_trmm
     module procedure c_trmm_simple,c_trmm_double
  end interface c_trmm
  interface axpy
     module procedure axpy_simple,axpy_double,axpy_simple_to_double
  end interface axpy
  interface c_axpy
     module procedure c_axpy_simple,c_axpy_double
  end interface c_axpy

contains

  !> Interfaces for LAPACK routines
  !! @warning
  !!   In these interfaces the input arrays are declared as scalars,
  !!   so the passage of the arguments by addresses is compulsory when calling
  !!   these routines
  !> Cholesky factorization of a positive definite matrix
  subroutine potrf_simple(uplo,n,a,lda,info)
    implicit none
    character(len=1), intent(in) :: uplo
    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a
    !call to LAPACK routine
    call spotrf(uplo,n,a,lda,info)
  end subroutine potrf_simple

  subroutine potrf_double(uplo,n,a,lda,info)
    implicit none
    character(len=1), intent(in) :: uplo
    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a
    !call to LAPACK routine
    call dpotrf(uplo,n,a,lda,info)
  end subroutine potrf_double

  subroutine c_potrf_simple(uplo,n,a,lda,info)
    implicit none
    character(len=1), intent(in) :: uplo
    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a
    !call to LAPACK routine
    call cpotrf(uplo,n,a,lda,info)
  end subroutine c_potrf_simple

  subroutine c_potrf_double(uplo,n,a,lda,info)
    implicit none
    character(len=1), intent(in) :: uplo
    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a
    !call to LAPACK routine
    call zpotrf(uplo,n,a,lda,info)
  end subroutine c_potrf_double

  !TRiangular matrix Inverse
  subroutine trtri_simple(uplo,diag,n,a,lda,info)
    implicit none
    character(len=1), intent(in) :: uplo,diag
    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a
    !call to LAPACK routine
    call strtri(uplo,diag,n,a,lda,info)
  end subroutine trtri_simple

  subroutine trtri_double(uplo,diag,n,a,lda,info)
    implicit none
    character(len=1), intent(in) :: uplo,diag
    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a
    !call to LAPACK routine
    call dtrtri(uplo,diag,n,a,lda,info)
  end subroutine trtri_double

  subroutine c_trtri_simple(uplo,diag,n,a,lda,info)
    implicit none
    character(len=1), intent(in) :: uplo,diag
    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a
    !call to LAPACK routine
    call ctrtri(uplo,diag,n,a,lda,info)
  end subroutine c_trtri_simple

  subroutine c_trtri_double(uplo,diag,n,a,lda,info)
    implicit none
    character(len=1), intent(in) :: uplo,diag
    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a
    !call to LAPACK routine
    call ztrtri(uplo,diag,n,a,lda,info)
  end subroutine c_trtri_double

  subroutine syev_simple(jobz,uplo,n,a,lda,w,work,lwork,info)
    implicit none
    character(len=1), intent(in) :: jobz,uplo
    integer, intent(in) :: lda,lwork,n
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a,work
    real(kind=4), intent(out) :: w
    !call to LAPACK routine
    call ssyev(jobz,uplo,n,a,lda,w,work,lwork,info)
  end subroutine syev_simple

  subroutine syev_double(jobz,uplo,n,a,lda,w,work,lwork,info)
    implicit none
    character(len=1), intent(in) :: jobz,uplo
    integer, intent(in) :: lda,lwork,n
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a,work
    real(kind=8), intent(out) :: w
    !call to LAPACK routine
    call dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
  end subroutine syev_double

  subroutine heev_simple(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
    implicit none
    character(len=1), intent(in) :: jobz,uplo
    integer, intent(in) :: lda,lwork,n
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a,work,rwork
    real(kind=4), intent(out) :: w
    !call to LAPACK routine
    call cheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
  end subroutine heev_simple

  subroutine heev_double(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
    implicit none
    character(len=1), intent(in) :: jobz,uplo
    integer, intent(in) :: lda,lwork,n
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a,work,rwork
    real(kind=8), intent(out) :: w
    !call to LAPACK routine
    call zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
  end subroutine heev_double

  subroutine sygv_simple(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
    implicit none
    character(len=1), intent(in) :: jobz,uplo
    integer, intent(in) :: itype,lda,ldb,lwork,n
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a,b,work
    real(kind=4), intent(out) :: w
    !call to LAPACK routine
    call ssygv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
  end subroutine sygv_simple

  subroutine sygv_double(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
    implicit none
    character(len=1), intent(in) :: jobz,uplo
    integer, intent(in) :: itype,lda,ldb,lwork,n
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a,b,work
    real(kind=8), intent(out) :: w
    !call to LAPACK routine
    call dsygv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
  end subroutine sygv_double

  subroutine hegv_simple(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
    implicit none
    character(len=1), intent(in) :: jobz,uplo
    integer, intent(in) :: lda,lwork,n,itype,ldb
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a,work,rwork,b
    real(kind=4), intent(out) :: w
    !call to LAPACK routine
    call chegv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
  end subroutine hegv_simple

  subroutine gesv_double(n,nrhs,a,lda,ipiv,b,ldb,info)
    implicit none
    integer, intent(in) :: n,lda,nrhs,ldb
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a,b
    integer, intent(out) :: ipiv
    !call to LAPACK routine
    call dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
  end subroutine gesv_double

  subroutine gesv_simple(n,nrhs,a,lda,ipiv,b,ldb,info)
    implicit none
    integer, intent(in) :: n,lda,nrhs,ldb
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a,b
    integer, intent(out) :: ipiv
    !call to LAPACK routine
    call sgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
  end subroutine gesv_simple

  subroutine hegv_double(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
    implicit none
    character(len=1), intent(in) :: jobz,uplo
    integer, intent(in) :: lda,lwork,n,itype,ldb
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a,work,rwork,b
    real(kind=8), intent(out) :: w
    !call to LAPACK routine
    call zhegv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
  end subroutine hegv_double

  subroutine c_gesv_double(n,nrhs,a,lda,ipiv,b,ldb,info)
    implicit none
    integer, intent(in) :: n,lda,nrhs,ldb
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a,b
    integer, intent(out) :: ipiv
    !call to LAPACK routine
    call zgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
  end subroutine c_gesv_double

  subroutine c_gesv_simple(n,nrhs,a,lda,ipiv,b,ldb,info)
    implicit none
    integer, intent(in) :: n,lda,nrhs,ldb
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a,b
    integer, intent(out) :: ipiv
    !call to LAPACK routine
    call cgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
  end subroutine c_gesv_simple


  !> Interfaces for BLAS routines
  !! @warning
  !!         In these interfaces the input arrays are declared as scalars,
  !!         so the passage of the arguments by addresses is compulsory when calling
  !!         these routines

  !SCALe a vector by a constant
  subroutine scal_simple(n,da,dx,incx)
    implicit none
    integer, intent(in) :: incx,n
    real(kind=4), intent(in) :: da
    real(kind=4), intent(inout) :: dx
    !call to BLAS routine
    call SSCAL(n,da,dx,incx)
  end subroutine scal_simple

  subroutine scal_double(n,da,dx,incx)
    implicit none
    integer, intent(in) :: incx,n
    real(kind=8), intent(in) :: da
    real(kind=8), intent(inout) :: dx
    !call to BLAS routine
    call DSCAL(n,da,dx,incx)
  end subroutine scal_double

  subroutine c_scal_simple(n,da,dx,incx)
    implicit none
    integer, intent(in) :: incx,n
    real(kind=4), intent(in) :: da
    real(kind=4), intent(out) :: dx
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_CSCAL(n,da,dx,incx)
    else
       !call to BLAS routine
       call CSCAL(n,da,dx,incx)
    end if
  end subroutine c_scal_simple

  subroutine c_scal_double(n,da,dx,incx)
    implicit none
    integer, intent(in) :: incx,n
    real(kind=8), intent(in) :: da
    real(kind=8), intent(out) :: dx
    !call to BLAS routine
    call ZSCAL(n,da,dx,incx)
  end subroutine c_scal_double

  !copy the vector
  subroutine copy_complex_real_simple(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    complex(kind=4), intent(in) :: dx
    real(kind=4), intent(out) :: dy
    !call to BLAS routine
    call SCOPY(n,dx,incx,dy,incy)
  end subroutine copy_complex_real_simple

  subroutine copy_complex_real_double(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    complex(kind=8), intent(in) :: dx
    real(kind=8), intent(out) :: dy
    !call to BLAS routine
    call DCOPY(n,dx,incx,dy,incy)
  end subroutine copy_complex_real_double

  subroutine copy_integer(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    integer, intent(in) :: dx
    integer, intent(out) :: dy
    !custom blas routine
    call icopy(n,dx,incx,dy,incy)
  end subroutine copy_integer

  subroutine copy_simple(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    real(kind=4), intent(in) :: dx
    real(kind=4), intent(out) :: dy
    !call to BLAS routine
    call SCOPY(n,dx,incx,dy,incy)
  end subroutine copy_simple

  subroutine copy_double(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    real(kind=8), intent(in) :: dx
    real(kind=8), intent(out) :: dy
    !call to BLAS routine
    call DCOPY(n,dx,incx,dy,incy)
  end subroutine copy_double

  subroutine copy_double_to_simple(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    real(kind=8), intent(in) :: dx
    real(kind=4), intent(out) :: dy
    !call to custom routine
    call dscopy(n,dx,incx,dy,incy)
  end subroutine copy_double_to_simple

  subroutine trmm_simple(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    implicit none
    character(len=1), intent(in) :: side,uplo,transa,diag
    integer, intent(in) :: lda,ldb,m,n
    real(kind=4), intent(in) :: alpha
    real(kind=4), intent(in) :: a
    real(kind=4), intent(inout) :: b
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_STRMM(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    else
       !call to BLAS routine
       call STRMM(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    end if
  end subroutine trmm_simple

  subroutine trmm_double(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    implicit none
    character(len=1), intent(in) :: side,uplo,transa,diag
    integer, intent(in) :: lda,ldb,m,n
    real(kind=8), intent(in) :: alpha
    real(kind=8), intent(in) :: a
    real(kind=8), intent(inout) :: b
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_DTRMM(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    else
       !call to BLAS routine
       call DTRMM(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    end if
  end subroutine trmm_double

  subroutine c_trmm_simple(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    implicit none
    character(len=1), intent(in) :: side,uplo,transa,diag
    integer, intent(in) :: lda,ldb,m,n
    complex(kind=4), intent(in) :: alpha
    real(kind=4), intent(in) :: a
    real(kind=4), intent(inout) :: b
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_CTRMM(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    else
       !call to BLAS routine
       call CTRMM(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    end if
  end subroutine c_trmm_simple

  subroutine c_trmm_double(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    implicit none
    character(len=1), intent(in) :: side,uplo,transa,diag
    integer, intent(in) :: lda,ldb,m,n
    complex(kind=8), intent(in) :: alpha
    real(kind=8), intent(in) :: a
    real(kind=8), intent(inout) :: b
    !call to BLAS routine
    call ZTRMM(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
  end subroutine c_trmm_double

  subroutine axpy_simple(n,da,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    real(kind=4), intent(in) :: da
    real(kind=4), intent(in) :: dx
    real(kind=4), intent(inout) :: dy
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_SAXPY(n,da,dx,incx,dy,incy)
    else
       !call to BLAS routine
       call SAXPY(n,da,dx,incx,dy,incy)
    end if
  end subroutine axpy_simple

  subroutine axpy_double(n,da,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    real(kind=8), intent(in) :: da
    real(kind=8), intent(in) :: dx
    real(kind=8), intent(inout) :: dy
    !call to BLAS routine
    call DAXPY(n,da,dx,incx,dy,incy)
  end subroutine axpy_double

  subroutine axpy_simple_to_double(n,da,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    real(kind=8), intent(in) :: da
    real(kind=4), intent(in) :: dx
    real(kind=8), intent(inout) :: dy
    !call to custom routine, for mixed precision sum
    call dasxpdy(n,da,dx,incx,dy,incy)
  end subroutine axpy_simple_to_double

  subroutine c_axpy_simple(n,da,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    real(kind=4), intent(in) :: da
    real(kind=4), intent(in) :: dx
    real(kind=4), intent(inout) :: dy
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_CAXPY(n,da,dx,incx,dy,incy)
    else
       !call to BLAS routine
       call CAXPY(n,da,dx,incx,dy,incy)
    end if
  end subroutine c_axpy_simple

  subroutine c_axpy_double(n,da,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    real(kind=8), intent(in) :: da
    real(kind=8), intent(in) :: dx
    real(kind=8), intent(inout) :: dy
    !call to BLAS routine
    call ZAXPY(n,da,dx,incx,dy,incy)
  end subroutine c_axpy_double

  !euclidean dot product
  function dot_simple(n,sx,incx,sy,incy)
    implicit none
    integer, intent(in) :: n,incx,incy
    real(kind=4), intent(in) :: sx,sy
    real(kind=4) :: dot_simple
    !local variables
    real(kind=4) :: cublas_sdot,sdot
    if (GPUblas) then
       !call to CUBLAS function
       dot_simple=cublas_sdot(n,sx,incx,sy,incy)
    else
       !call to BLAS function
       dot_simple=sdot(n,sx,incx,sy,incy)
    end if
  end function dot_simple

  !euclidean dot product
  function dotc_simple(n,sx,incx,sy,incy)
    implicit none
    integer, intent(in) :: n,incx,incy
    complex(kind=4), intent(in) :: sx,sy
    complex(kind=4) :: dotc_simple
    !local variables
    complex(kind=4) :: cdotc
    !call to BLAS function
    dotc_simple=cdotc(n,sx,incx,sy,incy)
  end function dotc_simple

  function dot_double(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: n,incx,incy
    real(kind=8), intent(in) :: dx,dy
    real(kind=8) :: dot_double
    !local variables
    real(kind=8) :: cublas_ddot,ddot
    if (GPUblas) then
       !call to CUBLAS function
       dot_double=cublas_ddot(n,dx,incx,dy,incy)
    else
       !call to BLAS function
       dot_double=ddot(n,dx,incx,dy,incy)
    end if
  end function dot_double

  function dotc_double(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: n,incx,incy
    complex(kind=8), intent(in) :: dx,dy
    complex(kind=8) :: dotc_double
    !local variables
    complex(kind=8) :: zdotc
    !call to BLAS function
    dotc_double=zdotc(n,dx,incx,dy,incy)
  end function dotc_double

  !euclidean NoRM of a vector
  function nrm2_simple(n,x,incx)
    implicit none
    integer, intent(in) :: n,incx
    real(kind=4), intent(in) :: x
    real(kind=4) :: nrm2_simple
    !local variables
    real(kind=4) :: cublas_snrm2,snrm2
    if (GPUblas .and. n>10000) then
       !call to CUBLAS function
       nrm2_simple=cublas_snrm2(n,x,incx)
    else
       !call to BLAS function
       nrm2_simple=snrm2(n,x,incx)
    end if
  end function nrm2_simple

  function nrm2_double(n,x,incx)
    implicit none
    integer, intent(in) :: n,incx
    real(kind=8), intent(in) :: x
    real(kind=8) :: nrm2_double
    !local variables
    real(kind=8) :: cublas_dnrm2,dnrm2
    if (GPUblas .and. n>10000) then
       !call to CUBLAS function
       nrm2_double=cublas_dnrm2(n,x,incx)
    else
       !call to BLAS routine
       nrm2_double=dnrm2(n,x,incx)
    end if
  end function nrm2_double

  !GEneral Matrix-Matrix multiplication routines
  subroutine gemm_simple(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: transa,transb
    integer, intent(in) :: k,lda,ldb,ldc,m,n
    real(kind=4), intent(in) :: alpha,beta
    real(kind=4), intent(in) :: a
    real(kind=4), intent(in) :: b
    real(kind=4), intent(inout) :: c
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_SGEMM(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    else
       !call to BLAS routine
       call SGEMM(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    end if
  end subroutine gemm_simple

  subroutine gemm_double(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: transa,transb
    integer, intent(in) :: k,lda,ldb,ldc,m,n
    real(kind=8), intent(in) :: alpha,beta
    real(kind=8), intent(in) :: a
    real(kind=8), intent(in) :: b
    real(kind=8), intent(inout) :: c
    !call to BLAS routine
    if (GPUblas) then
       call cublas_DGEMM(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    else
       call DGEMM(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    end if
  end subroutine gemm_double

  subroutine gemmsy_double_wrap(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: transa,transb
    integer, intent(in) :: k,lda,ldb,ldc,m,n
    real(kind=8), intent(in) :: alpha,beta
    real(kind=8), intent(in) :: a
    real(kind=8), intent(in) :: b
    real(kind=8), intent(inout) :: c
    !call to BLAS routine
    call gemmsy_double(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  end subroutine gemmsy_double_wrap

  subroutine c_gemm_simple(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: transa,transb
    integer, intent(in) :: k,lda,ldb,ldc,m,n
    complex(kind=4), intent(in) :: alpha,beta
    real(kind=4), intent(in) :: a
    real(kind=4), intent(in) :: b
    real(kind=4), intent(inout) :: c
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_CGEMM(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    else
       !call to BLAS routine
       call CGEMM(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    end if
  end subroutine c_gemm_simple

  subroutine c_gemm_double(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: transa,transb
    integer, intent(in) :: k,lda,ldb,ldc,m,n
    complex(kind=8), intent(in) :: alpha,beta
    real(kind=8), intent(in) :: a
    real(kind=8), intent(in) :: b
    real(kind=8), intent(inout) :: c
    !call to BLAS routine
    call ZGEMM(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  end subroutine c_gemm_double

  !SYmmetric Rank K operation
  subroutine syrk_simple(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: trans,uplo
    integer, intent(in) :: k,lda,ldc,n
    real(kind=4), intent(in) :: alpha,beta
    real(kind=4), intent(in) :: a
    real(kind=4), intent(out) :: c 
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_SSYRK(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    else
       !call to BLAS routine
       call SSYRK(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    end if
  end subroutine syrk_simple

  subroutine syrk_double(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: trans,uplo
    integer, intent(in) :: k,lda,ldc,n
    real(kind=8), intent(in) :: alpha,beta
    real(kind=8), intent(in) :: a
    real(kind=8), intent(out) :: c 
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_DSYRK(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    else
       !call to BLAS routine
       call DSYRK(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    end if
  end subroutine syrk_double

  !HErmitian Rank K operation
  subroutine herk_simple(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: trans,uplo
    integer, intent(in) :: k,lda,ldc,n
    real(kind=4), intent(in) :: alpha,beta
    real(kind=4), intent(in) :: a
    real(kind=4), intent(out) :: c 
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_CHERK(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    else
       !call to BLAS routine
       call CHERK(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    end if
  end subroutine herk_simple

  subroutine herk_double(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: trans,uplo
    integer, intent(in) :: k,lda,ldc,n
    real(kind=8), intent(in) :: alpha,beta
    real(kind=8), intent(in) :: a
    real(kind=8), intent(out) :: c 
    !call to BLAS routine
    call ZHERK(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
  end subroutine herk_double

  subroutine put_to_zero_simple(n,da)
    implicit none
    integer, intent(in) :: n
    real(kind=4), intent(out) :: da
    logical :: within_openmp
    !$ logical :: omp_in_parallel, omp_get_nested
    within_openmp=.false.
    !$    within_openmp=omp_in_parallel() .or. omp_get_nested()

    !call to custom routine
    if (.not. within_openmp) call timing(0,'Init to Zero  ','IR') 
    call razero_simple(n,da)
    if (.not. within_openmp) call timing(0,'Init to Zero  ','RS') 
  end subroutine put_to_zero_simple

  !!@todo To remove this routine which is not conformed to the Fortran standard (TD)
  subroutine put_to_zero_double(n,da)
    implicit none
    integer, intent(in) :: n
    real(kind=8), intent(out) :: da
    logical :: within_openmp
    !$ logical :: omp_in_parallel, omp_get_nested
    within_openmp=.false.
    !$    within_openmp=omp_in_parallel() .or. omp_get_nested()

    !call to custom routine
    !if (.not. within_openmp) call timing(0,'Init to Zero  ','IR') 
    call razero(n,da)
    !if (.not. within_openmp) call timing(0,'Init to Zero  ','RS') 
  end subroutine put_to_zero_double

  subroutine put_to_zero_double_1(n,da)
    implicit none
    integer, intent(in) :: n
    real(kind=8), dimension(:), intent(out) :: da
    logical :: within_openmp
    !$ logical :: omp_in_parallel,omp_get_nested
    within_openmp=.false.
    !$    within_openmp=omp_in_parallel() .or. omp_get_nested()

    !call to custom routine
    if (.not. within_openmp) call timing(0,'Init to Zero  ','IR') 
    call razero(n,da)
    if (.not. within_openmp) call timing(0,'Init to Zero  ','RS') 
  end subroutine put_to_zero_double_1

  subroutine put_to_zero_double_2(n,da)
    implicit none
    integer, intent(in) :: n
    real(kind=8), dimension(:,:), intent(out) :: da
    logical :: within_openmp
    !$ logical :: omp_in_parallel,omp_get_nested
    within_openmp=.false.
    !$    within_openmp=omp_in_parallel() .or. omp_get_nested()

    !call to custom routine
    if (.not. within_openmp) call timing(0,'Init to Zero  ','IR') 
    call razero(n,da)
    if (.not. within_openmp) call timing(0,'Init to Zero  ','RS') 
  end subroutine put_to_zero_double_2

  subroutine put_to_zero_integer(n,da)
    implicit none
    integer, intent(in) :: n
    integer, intent(out) :: da
    logical :: within_openmp
    !$ logical :: omp_in_parallel, omp_get_nested
    within_openmp=.false.
    !$    within_openmp=omp_in_parallel() .or. omp_get_nested()

    !call to custom routine
    if (.not. within_openmp) call timing(0,'Init to Zero  ','IR') 
    call razero_integer(n,da)
    if (.not. within_openmp) call timing(0,'Init to Zero  ','RS') 
  end subroutine put_to_zero_integer

end module wrapper_linalg
