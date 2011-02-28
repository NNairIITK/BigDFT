!!****m* BigDFT/module_base
!!
!! DESCRIPTION
!!  Modules which contains the low level definitions, as well as some profiling procedures
!!
!! AUTHOR
!!    Luigi Genovese
!!
!! COPYRIGHT
!!    Copyright (C) 2008-2010 CEA, ESRF
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!! 
#if defined HAVE_CONFIG_H
#include <config.inc>
#endif

module module_defs

  use m_profiling

  implicit none  

  ! Include variables set from configure.
  include 'configure.inc'

  ! Verbosity of the output, control the level of writing (minimal by default)
  integer :: verbose=2

  ! General precision, density and the wavefunctions types
  integer, parameter :: gp=kind(1.0d0)  !general-type precision
  integer, parameter :: dp=kind(1.0d0)  !density-type precision
  integer, parameter :: wp=kind(1.0d0)  !wavefunction-type precision
  integer, parameter :: tp=kind(1.0d0)  !diis precision (single in this context, if double is only for non-regression)

  ! MPI definitions and datatypes for density and wavefunctions
  include 'mpif.h'
  integer, parameter :: mpidtypw=MPI_DOUBLE_PRECISION,mpidtypd=MPI_DOUBLE_PRECISION
  integer, parameter :: mpidtypg=MPI_DOUBLE_PRECISION
#ifdef HAVE_MPI2
  ! Flag to use in the code to switch between MPI1 and MPI2
  logical, parameter :: have_mpi2 = .true.
#else
  ! Fake MPI_IN_PLACE variable to allow compilation in sumrho.
  integer :: MPI_IN_PLACE = 0
  ! Flag to use in the code to switch between MPI1 and MPI2
  logical, parameter :: have_mpi2 = .false.
#endif
  !integer, parameter :: mpidtypw=MPI_REAL,mpidtypd=MPI_REAL !in case of single precision

  !Flag for GPU computing, if CUDA libraries are present
  !in that case if a GPU is present a given MPI processor may or not perform a GPU calculation
  !this value can be changed in the read_input_variables routine
  logical :: GPUconv=.false.,GPUblas=.false.,GPUshare=.true.

  !Flag for GPU computing, if OpenCL libraries are present
  !in that case if a GPU is present a given MPI processor may or not perform a GPU calculation
  !this value can be changed in the read_input_variables routine
  logical :: OCLconv=.false.
  logical :: ASYNCconv=.true.

  !Logical parameter for the projectors application strategy (true for distributed way)
  !if the projector allocation passes the memorylimit this is switched to true
  !inside localize_projectors routines
  logical :: DistProjApply=.true.

  ! Physical constants.
  real(gp), parameter :: bohr2ang = 0.5291772108_gp                     ! 1 AU in angstroem
  real(gp), parameter :: ha2ev = 27.21138386_gp                         ! 1 Ha in eV
  real(gp), parameter :: Ha_cmm1=219474.6313705_gp                      ! 1 Hartree, in cm^-1 (from abinit 5.7.x)
  real(gp), parameter :: amu_emass=1.660538782e-27_gp/9.10938215e-31_gp ! 1 atomic mass unit, in electronic mass

  !interface for MPI_ALLREDUCE routine
  interface mpiallred
     module procedure mpiallred_int,mpiallred_real,mpiallred_double
  end interface


  !interfaces for LAPACK routines
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


  !interfaces for BLAS routines
  interface gemm
     module procedure gemm_simple,gemm_double
  end interface
  interface c_gemm
     module procedure c_gemm_simple,c_gemm_double
  end interface
  interface dot
     module procedure dot_simple,dot_double
  end interface
  interface nrm2
     module procedure nrm2_simple,nrm2_double
  end interface
  interface vscal
     module procedure scal_simple,scal_double
  end interface
  interface c_vscal
     module procedure c_scal_simple,c_scal_double
  end interface
  interface syrk
     module procedure syrk_simple,syrk_double
  end interface
  interface herk
     module procedure herk_simple,herk_double
  end interface
  interface trmm
     module procedure trmm_simple,trmm_double
  end interface
  interface c_trmm
     module procedure c_trmm_simple,c_trmm_double
  end interface
  interface axpy
     module procedure axpy_simple,axpy_double
  end interface
  interface c_axpy
     module procedure c_axpy_simple,c_axpy_double
  end interface

  contains
    
    !interface for MPI_ALLREDUCE operations
    subroutine mpiallred_int(buffer,ntot,mpi_op,mpi_comm,ierr)
      implicit none
      integer, intent(in) :: ntot,mpi_op,mpi_comm
      integer, intent(in) :: buffer
      integer, intent(out) :: ierr
#ifdef HAVE_MPI2
      !case with MPI_IN_PLACE
      call MPI_ALLREDUCE(MPI_IN_PLACE,buffer,ntot,&
           MPI_INTEGER,mpi_op,mpi_comm,ierr)
#else
      !local variables
      character(len=*), parameter :: subname='mpi_allred'
      integer :: i_all,i_stat
      integer, dimension(:), allocatable :: copybuf

      !case without mpi_in_place
      allocate(copybuf(ntot+ndebug),stat=i_stat)
      call memocc(i_stat,copybuf,'copybuf',subname)

      !not appropriate for integers, to be seen if it works
      call scopy(ntot,buffer,1,copybuf,1) 

      call MPI_ALLREDUCE(copybuf,buffer,ntot,&
           MPI_INTEGER,mpi_op,mpi_comm,ierr)
      
      i_all=-product(shape(copybuf))*kind(copybuf)
      deallocate(copybuf,stat=i_stat)
      call memocc(i_stat,i_all,'copybuf',subname)
#endif
    end subroutine mpiallred_int

    subroutine mpiallred_real(buffer,ntot,mpi_op,mpi_comm,ierr)
      implicit none
      integer, intent(in) :: ntot,mpi_op,mpi_comm
      real(kind=4), intent(in) :: buffer
      integer, intent(out) :: ierr
#ifdef HAVE_MPI2
      !case with MPI_IN_PLACE
      call MPI_ALLREDUCE(MPI_IN_PLACE,buffer,ntot,&
           MPI_REAL,mpi_op,mpi_comm,ierr)
#else
      !local variables
      character(len=*), parameter :: subname='mpi_allred'
      integer :: i_all,i_stat
      real(kind=4), dimension(:), allocatable :: copybuf

      !case without mpi_in_place
      allocate(copybuf(ntot+ndebug),stat=i_stat)
      call memocc(i_stat,copybuf,'copybuf',subname)
      
      call scopy(ntot,buffer,1,copybuf,1) 

      call MPI_ALLREDUCE(copybuf,buffer,ntot,&
           MPI_REAL,mpi_op,mpi_comm,ierr)
      
      i_all=-product(shape(copybuf))*kind(copybuf)
      deallocate(copybuf,stat=i_stat)
      call memocc(i_stat,i_all,'copybuf',subname)
#endif
    end subroutine mpiallred_real

    subroutine mpiallred_double(buffer,ntot,mpi_op,mpi_comm,ierr)
      implicit none
      integer, intent(in) :: ntot,mpi_op,mpi_comm
      real(kind=8), intent(in) :: buffer
      integer, intent(out) :: ierr
#ifdef HAVE_MPI2
      !case with MPI_IN_PLACE
      call MPI_ALLREDUCE(MPI_IN_PLACE,buffer,ntot,&
           MPI_DOUBLE_PRECISION,mpi_op,mpi_comm,ierr)
#else
      !local variables
      character(len=*), parameter :: subname='mpi_allred'
      integer :: i_all,i_stat
      real(kind=8), dimension(:), allocatable :: copybuf

      !case without mpi_in_place
      allocate(copybuf(ntot+ndebug),stat=i_stat)
      call memocc(i_stat,copybuf,'copybuf',subname)
      
      call dcopy(ntot,buffer,1,copybuf,1) 

      call MPI_ALLREDUCE(copybuf,buffer,ntot,&
           MPI_DOUBLE_PRECISION,mpi_op,mpi_comm,ierr)
      
      i_all=-product(shape(copybuf))*kind(copybuf)
      deallocate(copybuf,stat=i_stat)
      call memocc(i_stat,i_all,'copybuf',subname)
#endif
    end subroutine mpiallred_double
    


    !interfaces for LAPACK routines
    !WARNING: in these interfaces the input arrays are declared as scalars,
    !         so the passage of the arguments by addresses is compulsory when calling
    !         these routines
    
    !Cholesky factorization of a positive definite matrix
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


    !interfaces for BLAS routines
    !WARNING: in these interfaces the input arrays are declared as scalars,
    !         so the passage of the arguments by addresses is compulsory when calling
    !         these routines

    !SCALe a vector by a constant
    subroutine scal_simple(n,da,dx,incx)
      implicit none
      integer, intent(in) :: incx,n
      real(kind=4), intent(in) :: da
      real(kind=4), intent(out) :: dx
      if (GPUblas) then
         !call to CUBLAS routine
         call cublas_SSCAL(n,da,dx,incx)
      else
         !call to BLAS routine
         call SSCAL(n,da,dx,incx)
      end if
    end subroutine scal_simple

    subroutine scal_double(n,da,dx,incx)
      implicit none
      integer, intent(in) :: incx,n
      real(kind=8), intent(in) :: da
      real(kind=8), intent(out) :: dx
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

    !euclidean NoRM of a vector
    function nrm2_simple(n,x,incx)
      implicit none
      integer, intent(in) :: n,incx
      real(kind=4), intent(in) :: x
      real(kind=4) :: nrm2_simple
      !local variables
      real(kind=4) :: cublas_snrm2,snrm2
      if (GPUblas) then
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
      if (GPUblas) then
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

end module module_defs
!!***
