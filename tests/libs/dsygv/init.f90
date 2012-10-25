subroutine init_matrices(n, A, S)
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: n
  real(kind=8),dimension(n,n),intent(out) :: A, S

  ! Local variables
  integer :: i, j, ii
  real(4) :: builtin_rand

  ii=0
  do i=1,n
      do j=i,n
          !A(j,i)=dble(j)
          A(j,i)=dble(builtin_rand(ii))
          A(i,j)=A(j,i)
          if(i==j) then
              S(i,j)=1.d0
          else if(abs(i-j)==1) then
              S(i,j)=1.d-1
          else if(abs(i-j)==2) then
              S(i,j)=1.d-2
          else
              S(i,j)=0.d0
          end if
          S(j,i)=S(i,j)
      end do
  end do

end subroutine init_matrices
