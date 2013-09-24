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

