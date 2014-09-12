module module_io
    implicit none

    private

    public :: read_mode
    public :: write_mode

contains

subroutine read_mode(nat,filename,minmode)
    use module_base
    use module_types
    use module_atoms, only: atomic_structure,&
                            deallocate_atomic_structure,&
                            read_atomic_file=>set_astruct_from_file
    use module_global_variables, only: iproc
    implicit none
    !parameters
    integer, intent(in) :: nat
    character(len=*), intent(in) :: filename
    real(gp), intent(inout) :: minmode(3,nat)
    !internal
    type(atomic_structure):: astruct !< Contains all info


    call read_atomic_file(filename,iproc,astruct)
    if(nat/=astruct%nat)stop '(MHGPS) severe error in read_mode: '//&
                             'nat/=astruct%nat'
    if (trim(astruct%units) /= 'atomic'&
       .and. trim(astruct%units) /= 'atomicd0'&
       .and. trim(astruct%units) /= 'bohr'&
       .and. trim(astruct%units) /= 'bohrd0') then
        stop '(MHGPS) severe error units of mode-file must be atomic'
    endif

    call vcopy(3 * astruct%nat,astruct%rxyz(1,1),1,minmode(1,1), 1)
    call deallocate_atomic_structure(astruct)
end subroutine
subroutine write_mode(nat,filename,minmode,rotforce)
    use module_base, only: gp
    use module_types
    use module_interfaces
    use module_atoms, only: read_atomic_file=>set_astruct_from_file
    use module_global_variables, only: iproc, atoms, ixyz_int
    implicit none
    !parameters
    integer, intent(in) :: nat
    character(len=*), intent(in) :: filename
    real(gp), intent(in) :: minmode(3,nat)
    real(gp), intent(in), optional :: rotforce(3,nat)
    !internal
    character(len=7) :: comment='minmode'
    character(len=11) :: units

    units=atoms%astruct%units
    atoms%astruct%units='atomicd0'
    if(present(rotforce))then
        call write_atomic_file(filename,&
              0.0_gp,minmode(1,1),ixyz_int,&
              atoms,trim(comment),forces=rotforce(1,1))
    else
        call write_atomic_file(filename,&
              0.0_gp,minmode(1,1),ixyz_int,&
              atoms,trim(comment))
    endif
    atoms%astruct%units=units
end subroutine

end module
