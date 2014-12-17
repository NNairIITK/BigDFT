program mhgpstool
    use module_base
    use module_types
    use module_interfaces
    use yaml_output
    use module_atoms, only: set_astruct_from_file,astruct_dump_to_file
    use module_mhgpstool
    use module_userinput
    implicit none
    type(atoms_data) :: atoms
    integer :: nsadmax, nfolder
    character(len=500), allocatable :: folders(:)
    type(userinput) :: mhgps_uinp
    real(gp) :: en_delta_min, fp_delta_min
    real(gp) :: en_delta_sad, fp_delta_sad


!    call f_lib_initialize()
!
!    call read_folders(nfolder,folders)
!    call count_saddle_points(nfolder,folders,nsadmax)
!    call read_input(mhgps_uinp)
!    write(*,*) 'hello world'
!!    call set_astruct_from_file(trim(fileFrom),0,at%astruct,fcomment,energy,fxyz)
!
!    call f_free_str(500,folders)
!    call f_lib_finalize()
end program mhgpstool
