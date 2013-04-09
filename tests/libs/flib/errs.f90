subroutine test_error_handling()
  use yaml_output
  use error_handling
  implicit none
  !local variables
  integer :: ival,err1,ERR_TOTO,ERR_TITI
  external :: abort_toto,abort_titi

  call yaml_comment('Error Handling Module Test',hfill='~')
   
  call f_err_initialize()
  
  call f_err_set_callback(abort1)
  
  call f_err_define(err_name='ERR_TOTO',err_msg='test1',err_id=ERR_TOTO,callback=abort_toto)
  call f_err_define(err_name='ERR_TITI',err_msg='test2',err_id=ERR_TITI,callback=abort_titi)
  
  if (f_err_raise(condition=.true.,err_id=ERR_TOTO)) continue ! return
  
  call yaml_map("Callback done, errcode",ERR_TOTO)
  !call yaml_comment("HERE1")
  call f_err_set_callback(abort2)
  !call yaml_comment("HERE2")
  call yaml_map("Callback done",f_err_raise(.true.,err_id=ERR_TITI))
  !call yaml_comment("HERE3")
  call yaml_map("Error check value",f_err_check())
  !call yaml_comment("HERE4")
  call yaml_map("Error check code",f_err_check(err_id=ERR_TOTO))
  call yaml_map("Error check code2",f_err_check(err_id=ERR_TITI))
  call f_err_finalize()

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
