subroutine gemmsy_double(transa,transb,n,p,alpha,a,lda,b,ldb,beta,y,ldy)
  implicit none
  character(len=1), intent (in) :: transa, transb
  integer, intent(in) :: n, p, lda, ldb, ldy
  real(kind=8), intent(in) :: alpha,beta
  real(kind=8), intent(in) :: a
  real(kind=8), intent(in) :: b
  real(kind=8), intent(inout) :: y
  call dgemmsy(transa,transb,n,p,alpha,a,lda,b,ldb,beta,y,ldy)

end subroutine gemmsy_double
