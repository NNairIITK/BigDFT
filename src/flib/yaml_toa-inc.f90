  !body of the yaml_to template. Used in yaml_strings.f90
  character(len=max_value_length) :: str
  character(len=*), optional, intent(in) :: fmt

  !print *,'here',data,fmt
  str=repeat(' ',max_value_length)
  if (present(fmt)) then
     write(str,fmt) data
     !if the format has failed the result is full of stars
     if (trim(str) == repeat('*',len_trim(str))) write(str,cnv_fmt(data)) data
  else
     write(str,cnv_fmt(data)) data
  end if
  !otherwise write it in free format
  if (trim(str) == repeat('*',len_trim(str))) write(str,*) data
  !print *,'hereagain',str,data,fmt
  str=yaml_adjust(str)
