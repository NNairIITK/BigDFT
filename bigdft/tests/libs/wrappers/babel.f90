program babel
  use futile
  use module_atoms
  character(len=*), parameter :: input1=&
       "  {name: input, shortname: i, default: None,"//&
       "  help_string: Input file,"//&
       "  help_dict: {Allowed values: filepath}}"
  character(len=*), parameter :: input2=&
       "  {name: output,shortname: o, default: outfile.xyz,"//&
       "  help_string: Output file,"//&
       "  help_dict: {Allowed values: filepath}}"

  character(len=*), parameter :: inputs=&
       '-'//input1//f_cr//&
       '-'//input2

  character(len=64) :: fin, fout
  type(dictionary), pointer :: dict,types,options,iter
  type(atomic_structure) :: astruct

  interface
     subroutine openbabel_load(d, f,ln)
       use dictionaries
       implicit none
       type(dictionary), pointer :: d
       character(len = *), intent(in) :: f
       integer, intent(in) :: ln
     end subroutine openbabel_load
     subroutine openbabel_dump(d, t, f,ln)
       use dictionaries
       implicit none
       type(dictionary), pointer :: d, t
       character(len = *), intent(in) :: f
       integer, intent(in) :: ln
     end subroutine openbabel_dump
  end interface

  call f_lib_initialize()

  call yaml_argparse(options,inputs)

  call yaml_comment('Welcome to the F90 openbabel wrapper',hfill='-')

  fin=options//'input'
  fout=options//'output'
  
  call yaml_mapping_open('Reading positions')

  !dict=>dict_new()
  call dict_init(dict)
  call openbabel_load(dict,fin,len_trim(fin))

  call yaml_map(fin,dict)
  call yaml_mapping_close()

  call astruct_dict_get_types(dict, types)
  nullify(iter)
  do while (iterating(iter, on = types))
     call set(iter, dict_key(iter))
  end do
  
  call openbabel_dump(dict,types, fout,len_trim(fout))
  call yaml_map('Positions dumped into file',fout)

  call dict_free(options,dict, types)

  ! Reload the generated file.
  astruct = atomic_structure_null()
  call set_astruct_from_file("outfile.xyz", 0, astruct)

  call dict_init(dict)
  call astruct_merge_to_dict(dict, astruct, astruct%rxyz)
  call deallocate_atomic_structure(astruct)

  call yaml_map("outfile.xyz", dict)

  call dict_free(dict)

  call f_lib_finalize()
end program babel
