program mhgps
    use module_base
    use module_types
    use module_interfaces
    use module_input_dicts
    use yaml_output
    use module_atoms, only: deallocate_atoms_data,&
                            atomic_structure,&
                            read_atomic_file=>set_astruct_from_file
    use module_global_variables
    use module_init
    use module_saddle
    implicit none
    integer :: bigdft_get_number_of_atoms,bigdft_get_number_of_orbitals
    character(len=*), parameter :: subname='mhgps'
    type(globals) :: glob
    integer :: iproc=1
    type(atomic_structure) :: astruct1, astruct2

    write(*,*)'test'
    call f_lib_initialize()

    call init_global_variables(glob)

!    call read_atomic_file('posinp',iproc,astruct)











    call f_lib_finalize()

end program
