module interfaces
  implicit none

  interface

    subroutine init_matrices(n, A, S)
      implicit none
      integer,intent(in) :: n
      real(kind=8),dimension(n,n),intent(out) :: A, S
    end subroutine init_matrices

    subroutine dsygv_wrapper(itype, jobz, uplo, n, A, ldA, S, ldS, eval)
      implicit none
      integer,intent(in) :: itype, n, ldA, ldS
      character(len=1),intent(in) :: jobz, uplo
      real(kind=8),dimension(ldA,n),intent(inout) :: A
      real(kind=8),dimension(ldS,n),intent(inout) :: S
      real(kind=8),dimension(n),intent(inout) :: eval
    end subroutine dsygv_wrapper

    subroutine pdsygvx_wrapper(iproc, nproc, blocksize, comm, itype, jobz, uplo, n, a, lda, b, ldb, w)
      implicit none
      integer,intent(in) :: iproc, nproc, blocksize, comm, itype, n, lda, ldb
      character(len=1),intent(in) :: jobz, uplo
      real(kind=8),dimension(lda,n),intent(inout) :: a
      real(kind=8),dimension(ldb,n),intent(inout) :: b
      real(kind=8),dimension(n),intent(out) :: w
    end subroutine pdsygvx_wrapper

  end interface

end module interfaces
