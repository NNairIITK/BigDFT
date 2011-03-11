subroutine gemmsy_double(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,y,ldy)
  use module_base
  implicit none
  character(len=1), intent (in) :: transa, transb
  integer, intent(in) :: m, n, k, lda, ldb, ldy
  real(kind=8), intent(in) :: alpha,beta
  real(kind=8), intent(in) :: a
  real(kind=8), intent(in) :: b
  real(kind=8), intent(inout) :: y
  call gemm_double(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,y,ldy)

END SUBROUTINE gemmsy_double
