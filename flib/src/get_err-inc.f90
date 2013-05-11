!> piece of routine to identify the errors, consider two optional arguments err_id and err_name
!! included in error_handling.f90
  !local variables
  integer :: get_error 
  integer :: nerr,ierr,jerr
  character(len=dict_msg_len) :: name

  get_error=-1 !no error specified
  nerr=dict_len(dict_present_error)
  if (present(err_name)) then
     get_error=0
     do ierr=0,nerr-1
        !this one can be substituted by the values of the dictionary
        jerr=dict_present_error//ierr
        name=dict_key(dict_errors//jerr)
        if (trim(name)==trim(err_name)) then
           get_error=1 !name
           exit
        end if
     end do
  else if (present(err_id)) then
     get_error=0
     do ierr=0,nerr-1
        jerr=dict_present_error//ierr
        if (jerr==err_id) then
           get_error=2
           exit
        end if
     end do
  end if
