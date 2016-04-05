module f_python
  implicit none

  interface
     subroutine f_python_initialize(iproc, nproc, igroup, ngroup)
       use dictionaries
       implicit none
       integer, intent(in) :: iproc, nproc, igroup, ngroup
     end subroutine f_python_initialize
  end interface

  interface
     subroutine f_python_finalize()
       use dictionaries
       implicit none
     end subroutine f_python_finalize
  end interface

  interface
     subroutine f_python_execute_dict(dict, status)
       use dictionaries
       implicit none
       type(dictionary), pointer :: dict
       integer, intent(out) :: status
     end subroutine f_python_execute_dict
  end interface

  interface
     subroutine f_python_execute(script, status)
       use dictionaries
       implicit none
       character(len = *), intent(in) :: script
       integer, intent(out) :: status
     end subroutine f_python_execute
  end interface
end module f_python
