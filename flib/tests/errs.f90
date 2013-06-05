subroutine test_error_handling()
  use yaml_output
  use dictionaries!error_handling
  implicit none
  !local variables
  integer :: ival,err1,ERR_TOTO,ERR_TITI,ERR_GRAVE
  external :: abort_toto,abort_titi,abort1,abort2

  call yaml_comment('Error Handling Module Test',hfill='~')
   
!!$  print *,'address',f_loc(abort1)
!!$  print *,'address',f_loc(ival)

  call f_err_initialize()
  
!!$!  call f_err_set_callback(abort1)
  call f_err_severe_override(abort2)
  
  call f_err_define(err_name='ERR_TOTO',&
       err_msg='This is the error message for the error of kind 1 and it is written extensively'//&
       ' on purpose to see whether yaml module prints it',&
       err_action='For this error, contact the routine developer at mail at univ dot gov',&
       err_id=ERR_TOTO,callback=abort_toto)

  call f_err_define(err_name='ERR_TITI',err_msg='test2',err_id=ERR_TITI,&
       callback=abort_titi,callback_data=f_loc(ival))

  call f_err_define(err_name='ERR_GRAVE',err_msg='test2',err_id=ERR_GRAVE,&
       callback=f_err_severe)
  call yaml_map("Raising the TOTO error, errcode",ERR_TOTO) 

  if (f_err_raise(.true.,'Extra message added',err_id=ERR_TOTO)) continue ! return
    call yaml_map("Raising the TOTO error, by name, without condition",'ERR_TOTO') 
  if (f_err_raise(err_msg='Extra message added again',err_name='ERR_TOTO')) continue ! return
  
  call yaml_map("Callback done, errcode",ERR_TOTO)

!  call f_err_severe_restore()
  if (f_err_raise(.true.,'Generic error raised, some message here')) continue ! return

  call f_err_clean()
  call f_err_set_callback(abort2)

  call yaml_map("Callback done",f_err_raise(.true.,'Now TITI error has been raised',err_id=ERR_TITI))
  call yaml_map("Error check value",f_err_check())
  call yaml_map("Error check code",f_err_check(err_id=ERR_TOTO))
  call yaml_map("Error check code2",f_err_check(err_id=ERR_TITI))

end subroutine test_error_handling

subroutine abort1()
  use yaml_output
  implicit none
  call yaml_comment('Ouille',hfill='!')
end subroutine abort1

subroutine abort2()
  use yaml_output
  implicit none
  call yaml_comment('Aie',hfill='!')
end subroutine abort2

subroutine abort_toto()
  use yaml_output
  implicit none
  call yaml_comment('TOTO',hfill='!')
end subroutine abort_toto

subroutine abort_titi()
  use yaml_output
  implicit none
  call yaml_comment('TITI',hfill='!')
end subroutine abort_titi
