!!$subroutine get_i1(array,iadd)
!!$  implicit none
!!$  integer(kind=8), intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$
!!$  iadd=array
!!$end subroutine get_i1
!!$subroutine get_dp1(array,iadd)
!!$  implicit none
!!$  integer(kind=8), intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$
!!$  iadd=array
!!$end subroutine get_dp1
!!$subroutine get_dp2(array,iadd)
!!$  implicit none
!!$  integer(kind=8), intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$
!!$  iadd=array
!!$end subroutine get_dp2
!!$subroutine get_dp3(array,iadd)
!!$  implicit none
!!$  integer(kind=8), intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$
!!$  iadd=array
!!$end subroutine get_dp3
!!$subroutine get_dp4(array,iadd)
!!$  implicit none
!!$  integer(kind=8), intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$
!!$  iadd=array
!!$end subroutine get_dp4
!!$
!!$subroutine get_dp1_ptr(array,iadd)
!!$  implicit none
!!$  integer(kind=8), intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$
!!$  iadd=array
!!$end subroutine get_dp1_ptr
!!$
!!$subroutine get_dp2_ptr(array,iadd)
!!$  implicit none
!!$  integer(kind=8), intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$
!!$  iadd=array
!!$end subroutine get_dp2_ptr
!!$
!!$subroutine get_dp3_ptr(array,iadd)
!!$  implicit none
!!$  integer(kind=8), intent(in) :: array
!!$  integer(kind=8), intent(out) :: iadd
!!$
!!$  iadd=array
!!$end subroutine get_dp3_ptr

subroutine call_external(routine,args)
  implicit none
  external :: routine
  integer(kind=8), intent(in) :: args


  print *,'calling external, args',args
  
  if (args==0) then
     call routine()
  else
     call routine(args)
  end if
end subroutine call_external

subroutine call_external_f(routine)!,args)
  implicit none
  external :: routine
!  integer(kind=8), intent(in) :: args


!  print *,'calling external, args',args
  
!  if (args==0) then
     call routine()
!  else
!     call routine(args)
!  end if
end subroutine call_external_f

!> function which identify the address of the scalar object
!! associated to a unknown quantity
function f_loc(routine)
  implicit none
  external :: routine
  integer(kind=8) :: f_loc

  call getlongaddress(routine,f_loc)

end function f_loc

!> Callback routine for severe errors
subroutine f_err_severe_internal()
  implicit none
  stop 'Severe error, cannot proceed'
end subroutine f_err_severe_internal

