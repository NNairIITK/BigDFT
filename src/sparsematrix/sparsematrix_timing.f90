module sparsematrix_timing
  use time_profiling
  implicit none

  ! This module is public, such that all other modules using this one inherit all modules used in here

  ! Timings categories
  integer,public,save :: TCAT_SMAT_COMPRESSION = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_COMPRESSION_COMMUNICATION = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_TRANSFORMATION = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_MULTIPLICATION = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_INITIALIZATION = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_CME_AUXILIARY = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_CME_POLYNOMIALS = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_CME_COEFFICIENTS = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_HL_MATRIX_OPERATIONS = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_HL_MATRIX_COMMUNICATIONS = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_HL_MATRIX_CHECKS = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_HL_DGEMM = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_HL_DSYEV = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_HL_DSYGV = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_HL_DGESV = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_HL_DPOTRF = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_HL_DPOTRI = TIMING_UNINITIALIZED

  contains

    !> Switch on the timing categories for the sparse matrices
    subroutine sparsematrix_initialize_timing_categories()
      use time_profiling, only: f_timing_category_group,f_timing_category
      implicit none
      character(len=*), parameter :: smat_manip = 'sparsematrix manipulation'
      character(len=*), parameter :: smat_init = 'sparsematrix initialization'
      character(len=*), parameter :: smat_comm = 'sparsematrix communications'
      character(len=*), parameter :: cme = 'chebyshev matrix expansion'
      character(len=*), parameter :: smat_blaslapack = 'BLAS / LAPACK'
  
  
      ! Definition of the groups
      call f_timing_category_group(smat_manip, 'sparse matrix operations')
      call f_timing_category_group(smat_init, 'sparse matrix initialization')
      call f_timing_category_group(smat_comm, 'sparse matrix communications')
      call f_timing_category_group(cme, 'chebyshev matrix expansion')
      call f_timing_category_group(smat_blaslapack, 'BLAS/LAPACK')
  
      ! Define the timing categories
  
      ! Initialization timing
      call f_timing_category('sparse matrix initialization', smat_init, &
           'sparse matrix initialization', TCAT_SMAT_INITIALIZATION)
  
      ! Low level timings
      call f_timing_category('Sparse matrix compression', smat_manip, &
           '(un)compression of sparse matrices', TCAT_SMAT_COMPRESSION)
      call f_timing_category('Sparse matrix compression communication', smat_comm, &
           '(un)compression communication of sparse matrices', TCAT_SMAT_COMPRESSION_COMMUNICATION)
      call f_timing_category('Sparse matrix transformation', smat_manip, &
           'sparsity pattern transformation of sparse matrices', TCAT_SMAT_TRANSFORMATION)
  
      call f_timing_category('Sparse matrix multiplication', smat_manip, &
           'sparse matrix matrix multiplication', TCAT_SMAT_MULTIPLICATION)
  
      ! Chebyshev Matrix Expansion timing
      call f_timing_category('CME auxiliary', cme, &
           'Chebyshev matrix expansion auxiliary', TCAT_CME_AUXILIARY)
      call f_timing_category('CME polynomials', cme, &
           'Chebyshev matrix expansion polynomials', TCAT_CME_POLYNOMIALS)
      call f_timing_category('CME coefficients', cme, &
           'Chebyshev matrix expansion coefficients', TCAT_CME_COEFFICIENTS)
  
      ! High level timings
      call f_timing_category('highlevel matrix operations', smat_manip, &
           'highlevel matrix operations', TCAT_HL_MATRIX_OPERATIONS)
      call f_timing_category('highlevel matrix communications', smat_comm, &
           'highlevel matrix communications', TCAT_HL_MATRIX_COMMUNICATIONS)
      call f_timing_category('highlevel matrix checks', smat_manip, &
           'highlevel matrix checks', TCAT_HL_MATRIX_CHECKS)
  
      ! BLAS / LAPACK timing
      call f_timing_category('DGEMM', smat_blaslapack, &
           '(Sca)LAPACK DGEMM', TCAT_HL_DGEMM)
      call f_timing_category('DSYEV', smat_blaslapack, &
           '(Sca)LAPACK DSYEV', TCAT_SMAT_HL_DSYEV)
      call f_timing_category('DSYGV', smat_blaslapack, &
           '(Sca)LAPACK DSYGV', TCAT_SMAT_HL_DSYGV)
      call f_timing_category('DGESV', smat_blaslapack, &
           '(Sca)LAPACK DGESV', TCAT_SMAT_HL_DGESV)
      call f_timing_category('DPOTRF', smat_blaslapack, &
           '(Sca)LAPACK DPOTRF', TCAT_SMAT_HL_DPOTRF)
      call f_timing_category('DPOTRI', smat_blaslapack, &
           '(Sca)LAPACK DPOTRI', TCAT_SMAT_HL_DPOTRI)
  
    end subroutine sparsematrix_initialize_timing_categories

end module sparsematrix_timing
