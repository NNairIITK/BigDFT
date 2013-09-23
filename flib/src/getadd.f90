!> @file
!! Routines to get the address of objects
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Routine to call an external routine with an integer argument
subroutine call_external(routine,args)
  implicit none
  external :: routine                  !< Routine to be called
  integer(kind=8), intent(in) :: args  !< Argument of the called routine


  print *,'calling external, args',args
  
  if (args==0) then
     call routine()
  else
     call routine(args)
  end if
end subroutine call_external


!> Call the external routine with no argument
subroutine call_external_f(routine)!,args)
  implicit none
  external :: routine                  !< Routine to be called
!  integer(kind=8), intent(in) :: args


!  print *,'calling external, args',args
  
!  if (args==0) then

     call routine()
!  else
!     call routine(args)
!  end if

end subroutine call_external_f


!> Function which identify the address of the scalar object
!! associated to a unknown quantity
function f_loc(routine)
  implicit none
  external :: routine       !< Object
  integer(kind=8) :: f_loc  !< Address of the object routine

  call getlongaddress(routine,f_loc)

end function f_loc

