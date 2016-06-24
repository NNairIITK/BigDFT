module sparsematrix_errorhandling
  implicit none

  private

  !> Errorcodes
  integer,save,public :: SPARSEMATRIX_ALLOCATION_ERROR
  integer,save,public :: SPARSEMATRIX_MANIPULATION_ERROR
  integer,save,public :: SPARSEMATRIX_RUNTIME_ERROR
  integer,save,public :: SPARSEMATRIX_INITIALIZATION_ERROR

  !> Public routine
  public :: sparsematrix_init_errors

  contains

    !> Define the sparsematrix errors
    subroutine sparsematrix_init_errors()
      use dictionaries
      implicit none

      call f_err_define('SPARSEMATRIX_ALLOCATION_ERROR',&
           'a problem occured during the allocation of a sparse matrix',&
           SPARSEMATRIX_ALLOCATION_ERROR,&
           err_action='Check the calling arguments of the allocation routine')

      call f_err_define('SPARSEMATRIX_MANIPULATION_ERROR',&
           'a problem occured during the manipulation ((un)compression,sparsity pattern transformation) of a sparse matrix',&
           SPARSEMATRIX_MANIPULATION_ERROR,&
           err_action='Check the calling arguments of the manipulation routine and the array sizes')

      call f_err_define('SPARSEMATRIX_RUNTIME_ERROR',&
           'a general problem related to sparse matrices occured during runtime',&
           SPARSEMATRIX_MANIPULATION_ERROR,&
           err_action='Check the dedicated error message')

      call f_err_define('SPARSEMATRIX_INITIALIZATION_ERROR',&
           'a problem related to the initialization of a sparse matrix occured',&
           SPARSEMATRIX_MANIPULATION_ERROR,&
           err_action='Check the calling arguments and the dedicated error message')
  
      end subroutine sparsematrix_init_errors

end module sparsematrix_errorhandling
