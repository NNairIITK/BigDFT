!> @file
!!  File defining parameters for BigDFT package (processed by the build system)
!! @author
!!    Copyright (C) 2008-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!! @warning
!!   THIS FILE IS PROCESSED BY THE BUILD SYSTEM.

#if defined HAVE_CONFIG_H
#include <config.inc>
#endif

!> Modules which contains the low level definitions, as well as some profiling procedures
module module_defs

  use m_profiling

  implicit none  

  include 'configure.inc' !< Include variables set from configure.

  integer :: verbose=2    !< Verbosity of the output, control the level of writing (minimal by default)

  ! General precision, density and the wavefunctions types
  integer, parameter :: gp=kind(1.0d0)  !< general-type precision
  integer, parameter :: dp=kind(1.0d0)  !< density-type precision
  integer, parameter :: wp=kind(1.0d0)  !< wavefunction-type precision
  integer, parameter :: tp=kind(1.0d0)  !< diis precision (single in this context, if double is only for non-regression)

  include 'mpif.h'      !< MPI definitions and datatypes for density and wavefunctions

  integer, parameter :: mpidtypw=MPI_DOUBLE_PRECISION
  integer, parameter :: mpidtypd=MPI_DOUBLE_PRECISION
  integer, parameter :: mpidtypg=MPI_DOUBLE_PRECISION
  !integer, parameter :: mpidtypw=MPI_REAL,mpidtypd=MPI_REAL !in case of single precision

#ifdef HAVE_MPI2
  logical, parameter :: have_mpi2 = .true.  !< Flag to use in the code to switch between MPI1 and MPI2
#else
  integer :: MPI_IN_PLACE = 0               !< Fake MPI_IN_PLACE variable to allow compilation in sumrho.
  logical, parameter :: have_mpi2 = .false. !< Flag to use in the code to switch between MPI1 and MPI2
#endif

  logical :: mpi_thread_funneled_is_supported=.false. !< control the OMP_NESTED based overlap, checked by bigdft_mpi_init below

  !> Flag for GPU computing, if CUDA libraries are present
  !! in that case if a GPU is present a given MPI processor may or not perform a GPU calculation
  !! this value can be changed in the read_input_variables routine
  logical :: GPUconv=.false.,GPUblas=.false.,GPUshare=.true.

  !> Flag for GPU computing, if OpenCL libraries are present
  !! in that case if a GPU is present a given MPI processor may or not perform a GPU calculation
  !! this value can be changed in the read_input_variables routine
  logical :: OCLconv=.false.
  logical :: ASYNCconv=.true.

  !> Logical parameter for the projectors application strategy (true for distributed way)
  !! if the projector allocation passes the memorylimit this is switched to true
  !! inside localize_projectors routines
  logical :: DistProjApply=.true.

  !> experimental variables to test the add of new functionalities
  logical :: experimental_modulebase_var_onlyfion=.false.


  !> Physical constants.
  real(gp), parameter :: bohr2ang = 0.5291772108_gp                     !> 1 AU in angstroem
  real(gp), parameter :: ha2ev = 27.21138386_gp                         !> 1 Ha in eV
  real(gp), parameter :: Ha_cmm1=219474.6313705_gp                      !> 1 Hartree, in cm^-1 (from abinit 5.7.x)
  real(gp), parameter :: Ha_eV=27.21138386_gp                           !> 1 Hartree, in eV
  real(gp), parameter :: Ha_K=315774.65_gp                              !> 1Hartree, in Kelvin
  real(gp), parameter :: Ha_THz=6579.683920722_gp                       !> 1 Hartree, in THz
  real(gp), parameter :: Ha_J=4.35974394d-18                            !> 1 Hartree, in J
  real(gp), parameter :: e_Cb=1.602176487d-19                           !> minus the electron charge, in Coulomb
  real(gp), parameter :: kb_HaK=8.617343d-5/Ha_eV                       !> Boltzmann constant in Ha/K
  real(gp), parameter :: amu_emass=1.660538782e-27_gp/9.10938215e-31_gp !> 1 atomic mass unit, in electronic mass
  real(gp), parameter :: GPaoAU=29421.010901602753                       !> 1Ha/Bohr^3 in GPa

  !> Evergreens
  real(dp), parameter :: pi_param=3.141592653589793238462643383279502884197_dp

  !> Code constants.
  !real(gp), parameter :: UNINITIALISED = -123456789._gp

  !> interface for MPI_ALLREDUCE routine
  interface mpiallred
     module procedure mpiallred_int,mpiallred_real,mpiallred_double,mpiallred_log
  end interface

  !interface for uninitialized variable
  interface UNINITIALIZED
     module procedure uninitialized_dbl,uninitialized_int,uninitialized_real,uninitialized_long
  end interface

  !initialize to zero an array
  interface to_zero
     module procedure put_to_zero_simple, &
          & put_to_zero_double, put_to_zero_double_1, put_to_zero_double_2, &
          & put_to_zero_integer
  end interface


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
     module procedure axpy_simple,axpy_double,axpy_simple_to_double
  end interface
  interface c_axpy
     module procedure c_axpy_simple,c_axpy_double
  end interface

  interface
     subroutine bigdft_utils_flush(unit)
       integer, intent(in) :: unit
     end subroutine bigdft_utils_flush
  end interface

  contains

    subroutine bigdft_mpi_init(ierr)
      implicit none
      integer, intent(out) :: ierr
#ifdef HAVE_MPI_INIT_THREAD
      integer :: provided
      call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
      if (ierr /= MPI_SUCCESS) then
         write(*,*)'BigDFT_mpi_INIT: Error in MPI_INIT_THREAD',ierr
      else if (provided < MPI_THREAD_FUNNELED) then
         !write(*,*)'WARNING: MPI_THREAD_FUNNELED not supported!',provided,ierr
         !call MPI_INIT(ierr)
      else
          mpi_thread_funneled_is_supported=.true.
      endif
#else
      call MPI_INIT(ierr)      
      if (ierr /= MPI_SUCCESS) then
         write(*,*)'BigDFT_mpi_INIT: Error in MPI_INIT_THREAD',ierr
      end if
#endif
    end subroutine bigdft_mpi_init

    !> Activates the nesting for UNBLOCK_COMMS performance case
    subroutine bigdft_open_nesting(num_threads)
      implicit none
      integer, intent(in) :: num_threads
#ifdef HAVE_MPI_INIT_THREAD
      !$ call OMP_SET_NESTED(.true.) 
      !$ call OMP_SET_MAX_ACTIVE_LEVELS(2)
      !$ call OMP_SET_NUM_THREADS(num_threads)
#else
      integer :: ierr,idummy
      write(*,*)'BigDFT_open_nesting is not active!'
      call MPI_ABORT(MPI_COMM_WORLD,ierr)
      idummy=num_threads
#endif
    end subroutine bigdft_open_nesting

    !> Activates the nesting for UNBLOCK_COMMS performance case
    subroutine bigdft_close_nesting(num_threads)
      implicit none
      integer, intent(in) :: num_threads
#ifdef HAVE_MPI_INIT_THREAD
      !$ call OMP_SET_NESTED(.false.) 
      !$ call OMP_SET_NUM_THREADS(num_threads)
#else 
      integer :: ierr,idummy
      write(*,*)'BigDFT_close_nesting is not active!'
      call MPI_ABORT(MPI_COMM_WORLD,ierr)
      idummy=num_threads
#endif
    end subroutine bigdft_close_nesting

    
    !interface for MPI_ALLREDUCE operations
    subroutine mpiallred_int(buffer,ntot,mpi_op,mpi_comm,ierr)
      implicit none
      integer, intent(in) :: ntot,mpi_op,mpi_comm
      integer, intent(inout) :: buffer
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
      ierr=0 !put just for MPIfake compatibility
      call MPI_ALLREDUCE(copybuf,buffer,ntot,&
           MPI_INTEGER,mpi_op,mpi_comm,ierr)
      
      i_all=-product(shape(copybuf))*kind(copybuf)
      deallocate(copybuf,stat=i_stat)
      call memocc(i_stat,i_all,'copybuf',subname)
#endif
      if (ierr /=0) stop 'MPIALLRED_INT'

    end subroutine mpiallred_int

    subroutine mpiallred_real(buffer,ntot,mpi_op,mpi_comm,ierr)
      implicit none
      integer, intent(in) :: ntot,mpi_op,mpi_comm
      real(kind=4), intent(inout) :: buffer
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
      ierr=0 !put just for MPIfake compatibility
      call MPI_ALLREDUCE(copybuf,buffer,ntot,&
           MPI_REAL,mpi_op,mpi_comm,ierr)
      
      i_all=-product(shape(copybuf))*kind(copybuf)
      deallocate(copybuf,stat=i_stat)
      call memocc(i_stat,i_all,'copybuf',subname)
#endif
      if (ierr /=0) stop 'MPIALLRED_REAL'

    end subroutine mpiallred_real

    subroutine mpiallred_double(buffer,ntot,mpi_op,mpi_comm,ierr)
      implicit none
      integer, intent(in) :: ntot,mpi_op,mpi_comm
      real(kind=8), intent(inout) :: buffer
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
      ierr=0 !put just for MPIfake compatibility
      call MPI_ALLREDUCE(copybuf,buffer,ntot,&
           MPI_DOUBLE_PRECISION,mpi_op,mpi_comm,ierr)
      
      i_all=-product(shape(copybuf))*kind(copybuf)
      deallocate(copybuf,stat=i_stat)
      call memocc(i_stat,i_all,'copybuf',subname)
#endif
      if (ierr /=0) stop 'MPIALLRED_DBL'
    end subroutine mpiallred_double

    !interface for MPI_ALLREDUCE operations
    subroutine mpiallred_log(buffer,ntot,mpi_op,mpi_comm,ierr)
      implicit none
      integer, intent(in) :: ntot,mpi_op,mpi_comm
      logical, intent(inout) :: buffer
      integer, intent(out) :: ierr
#ifdef HAVE_MPI2
      !case with MPI_IN_PLACE
      call MPI_ALLREDUCE(MPI_IN_PLACE,buffer,ntot,&
           MPI_LOGICAL,mpi_op,mpi_comm,ierr)
#else
      !local variables
      character(len=*), parameter :: subname='mpi_allred'
      integer :: i_all,i_stat
      logical, dimension(:), allocatable :: copybuf

      !case without mpi_in_place
      allocate(copybuf(ntot+ndebug),stat=i_stat)
      call memocc(i_stat,copybuf,'copybuf',subname)

      !not appropriate for logical, to be seen if it works
      call scopy(ntot,buffer,1,copybuf,1) 
      ierr=0 !put just for MPIfake compatibility
      call MPI_ALLREDUCE(copybuf,buffer,ntot,&
           MPI_LOGICAL,mpi_op,mpi_comm,ierr)
      
      i_all=-product(shape(copybuf))*kind(copybuf)
      deallocate(copybuf,stat=i_stat)
      call memocc(i_stat,i_all,'copybuf',subname)
#endif

      !inform and stop if an error occurs
      if (ierr /=0) stop 'MPIALLRED_LOG'

    end subroutine mpiallred_log

    function uninitialized_int(one) 
      implicit none
      integer(kind = 4), intent(in) :: one
      integer(kind = 4) :: uninitialized_int
      integer :: foo
      foo = kind(one)
      uninitialized_int=-123456789
    end function uninitialized_int

    function uninitialized_long(one) 
      implicit none
      integer(kind = 8), intent(in) :: one
      integer(kind = 8) :: uninitialized_long
      integer :: foo
      foo = kind(one)
      uninitialized_long=-123456789
    end function uninitialized_long

    function uninitialized_real(one) 
      implicit none
      real(kind=4), intent(in) :: one
      real(kind=4) :: uninitialized_real
      integer :: foo
      foo = kind(one)
      uninitialized_real=-123456789.e0
    end function uninitialized_real

    function uninitialized_dbl(one) 
      implicit none
      real(kind=8), intent(in) :: one
      real(kind=8) :: uninitialized_dbl
      integer :: foo
      foo = kind(one)
      uninitialized_dbl=-123456789.d0
    end function uninitialized_dbl

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

    subroutine put_to_zero_simple(n,da)
      implicit none
      integer, intent(in) :: n
      real(kind=4), intent(out) :: da
      logical :: within_openmp
      !$ logical :: omp_in_parallel
      within_openmp=.false.
      !$    within_openmp=omp_in_parallel()

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
      !$ logical :: omp_in_parallel
      within_openmp=.false.
      !$    within_openmp=omp_in_parallel()

      !call to custom routine
      if (.not. within_openmp) call timing(0,'Init to Zero  ','IR') 
      call razero(n,da)
      if (.not. within_openmp) call timing(0,'Init to Zero  ','RS') 
    end subroutine put_to_zero_double

    subroutine put_to_zero_double_1(n,da)
      implicit none
      integer, intent(in) :: n
      real(kind=8), dimension(n), intent(out) :: da
      logical :: within_openmp
      !$ logical :: omp_in_parallel
      within_openmp=.false.
      !$    within_openmp=omp_in_parallel()

      !call to custom routine
      if (.not. within_openmp) call timing(0,'Init to Zero  ','IR') 
      call razero(n,da)
      if (.not. within_openmp) call timing(0,'Init to Zero  ','RS') 
    end subroutine put_to_zero_double_1

    subroutine put_to_zero_double_2(n,da)
      implicit none
      integer, intent(in) :: n
      real(kind=8), dimension(n,*), intent(out) :: da
      logical :: within_openmp
      !$ logical :: omp_in_parallel
      within_openmp=.false.
      !$    within_openmp=omp_in_parallel()

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
      !$ logical :: omp_in_parallel
      within_openmp=.false.
      !$    within_openmp=omp_in_parallel()

      !call to custom routine
      if (.not. within_openmp) call timing(0,'Init to Zero  ','IR') 
      call razero_integer(n,da)
      if (.not. within_openmp) call timing(0,'Init to Zero  ','RS') 
    end subroutine put_to_zero_integer

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

    function fnrm_denpot(x,cplex,nfft,nspden,opt_denpot,user_data)
      use m_ab6_mixing
      implicit none
      integer, intent(in) :: cplex,nfft,nspden,opt_denpot
      double precision, intent(in) :: x(*)
      integer, intent(in) :: user_data(:)

      integer :: ierr, ie, iproc, npoints, ishift
      double precision :: fnrm_denpot, ar, nrm_local, dnrm2

      ! In case of density, we use nscatterarr.
      if (opt_denpot == AB6_MIXING_DENSITY) then
         call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
         if (ierr /= 0) then
            call MPI_ABORT(MPI_COMM_WORLD, ierr, ie)
         end if
         npoints = cplex * user_data(2 * iproc + 1)
         ishift  =         user_data(2 * iproc + 2)
      else
         npoints = cplex * nfft
         ishift  = 0
      end if

      ! Norm on spin up and spin down
      nrm_local = dnrm2(npoints * min(nspden,2), x(1 + ishift), 1)
      nrm_local = nrm_local ** 2

      if (nspden==4) then
         ! Add the magnetisation
         ar = dnrm2(npoints * 2, x(1 + cplex * nfft * 2 + ishift), 1)
         ar = ar ** 2
         if (opt_denpot == 0) then
            if (cplex == 1) then
               nrm_local = nrm_local + 2.d0 * ar
            else
               nrm_local = nrm_local + ar
            end if
         else
            nrm_local = 0.5d0 * (nrm_local + ar)
         end if
      end if
      
      ! Summarize on processors
      fnrm_denpot = nrm_local
      call MPI_ALLREDUCE(nrm_local, fnrm_denpot, 1, &
           & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      if (ierr /= 0) then
         call MPI_ABORT(MPI_COMM_WORLD, ierr, ie)
      end if
    end function fnrm_denpot

    function fdot_denpot(x,y,cplex,nfft,nspden,opt_denpot,user_data)
      use m_ab6_mixing
      implicit none
      integer, intent(in) :: cplex,nfft,nspden,opt_denpot
      double precision, intent(in) :: x(*), y(*)
      integer, intent(in) :: user_data(:)

      integer :: ierr, ie, iproc, npoints, ishift
      double precision :: fdot_denpot, ar, dot_local, ddot

      ! In case of density, we use nscatterarr.
      if (opt_denpot == AB6_MIXING_DENSITY) then
         call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
         if (ierr /= 0) then
            call MPI_ABORT(MPI_COMM_WORLD, ierr, ie)
         end if
         npoints = cplex * user_data(2 * iproc + 1)
         ishift  =         user_data(2 * iproc + 2)
      else
         npoints = cplex * nfft
         ishift  = 0
      end if

      if (opt_denpot == 0 .or. opt_denpot == 1) then
         ! Norm on spin up and spin down
         dot_local = ddot(npoints * min(nspden,2), x(1 + ishift), 1, y(1 + ishift), 1)

         if (nspden==4) then
            ! Add the magnetisation
            ar = ddot(npoints * 2, x(1 + cplex * nfft * 2 + ishift), 1, &
                 & y(1 + cplex * nfft * 2 + ishift), 1)
            if (opt_denpot == 0) then
               if (cplex == 1) then
                  dot_local = dot_local + 2.d0 * ar
               else
                  dot_local = dot_local + ar
               end if
            else
               dot_local = 0.5d0 * (dot_local + ar)
            end if
         end if
      else
         if(nspden==1)then
            dot_local = ddot(npoints, x(1 + ishift), 1, y(1 + ishift), 1)
         else if(nspden==2)then
            ! This is the spin up contribution
            dot_local = ddot(npoints, x(1 + ishift + nfft), 1, y(1 + ishift), 1)
            ! This is the spin down contribution
            dot_local = dot_local + ddot(npoints, x(1 + ishift ), 1, y(1 + ishift+ nfft), 1)
         else if(nspden==4)then
            !  \rho{\alpha,\beta} V^{\alpha,\beta} =
            !  rho*(V^{11}+V^{22})/2$
            !  + m_x Re(V^{12})- m_y Im{V^{12}}+ m_z(V^{11}-V^{22})/2
            dot_local = 0.5d0 * (ddot(npoints, x(1 + ishift), 1, y(1 + ishift), 1) + &
                 & dot(npoints, x(1 + ishift), 1, y(1 + ishift + nfft), 1))
            dot_local = dot_local + 0.5d0 * ( &
                 & ddot(npoints, x(1 + ishift + 3 * nfft), 1, y(1 + ishift), 1) - &
                 & ddot(npoints, x(1 + ishift + 3 * nfft), 1, y(1 + ishift + nfft), 1))
            dot_local = dot_local + &
                 & ddot(npoints, x(1 + ishift + nfft), 1, y(1 + ishift + 2 * nfft), 1)
            dot_local = dot_local - &
                 & ddot(npoints, x(1 + ishift + 2 * nfft), 1, y(1 + ishift + 3 * nfft), 1)
         end if
      end if
      
      ! Summarize on processors
      fdot_denpot = dot_local
      call MPI_ALLREDUCE(dot_local, fdot_denpot, 1, &
           & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      if (ierr /= 0) then
         call MPI_ABORT(MPI_COMM_WORLD, ierr, ie)
      end if
    end function fdot_denpot

#ifndef HAVE_FC_FLUSH
    integer function flush(unit)
      integer, intent(in) :: unit

      flush = unit
      return
    end function flush
#endif

end module module_defs
