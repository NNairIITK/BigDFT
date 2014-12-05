program mhgpstool
    use module_base
    use module_types
    use module_interfaces
    use yaml_output
    use module_atoms, only: set_astruct_from_file,astruct_dump_to_file
    implicit none
    type(atoms_data) :: atoms
    integer :: nsadmax, nfolder
    character(len=500), allocatable :: folders(:)
    interface
       subroutine read_folders(nfolder, folders)
         use module_base
         implicit none
         !parameters
         integer, intent(out) :: nfolder
         character(len=500), allocatable :: folders(:)
       end subroutine read_folders
    end interface


    call f_lib_initialize()

    call read_folders(nfolder,folders)
    call count_saddle_points(nsadmax)
    write(*,*) 'hello world'
!    call set_astruct_from_file(trim(fileFrom),0,at%astruct,fcomment,energy,fxyz)

    call f_lib_finalize()
end program mhgpstool
subroutine read_folders(nfolder,folders)
    use module_base
    implicit none
    !parameters
    integer, intent(out) :: nfolder
    character(len=500), allocatable :: folders(:)
    !internal
    integer :: u, istat
    character(len=600) :: line
real(8), allocatable :: dmy(:)
    u=f_get_free_unit()
    open(u,file='mhgpstool.inp')
    nfolder=0
    do
        read(u,*,iostat=istat)line
        if(istat/=0)exit
        nfolder=nfolder+1
    enddo
    folders = f_malloc_str(500,(/1.to.nfolder/),id='folders')
    close(u)
end subroutine read_folders
subroutine count_saddle_points(nsad)
    implicit none
    !parameters
    integer, intent(in) :: nsad
    !internal
end subroutine count_saddle_points
