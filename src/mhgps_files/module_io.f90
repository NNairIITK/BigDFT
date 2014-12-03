module module_io
    implicit none

    private

    public :: read_mode
    public :: write_mode
    public :: read_jobs
    public :: get_first_struct_file

contains

subroutine get_first_struct_file(filename)
    use module_base
    use module_global_variables, only: operation_mode,&
                                       currdir
    implicit none
    !parameter
    character(len=200), intent(out) :: filename
    !internal
    integer :: istat,u
    character(len=33) :: line
    logical :: exists


    inquire(file=trim(adjustl(currdir))//'/job_list',exist=exists)
    if(exists)then
        u=f_get_free_unit()
        open(unit=u,file=trim(adjustl(currdir))//'/job_list')
        read(u,'(a)',iostat=istat)line
        if(istat/=0)then
            call f_err_throw(trim(adjustl(currdir))//&
                                '/job_file empty')
        endif
        read(line,*)filename
        filename=trim(adjustl(currdir))//'/'//filename
        call check_struct_file_exists(filename)
        close(u)
    else
        write(filename,'(a,i3.3)')trim(adjustl(currdir))//'/pos',1
        call check_struct_file_exists(filename)
    endif
end subroutine
subroutine read_jobs(njobs,joblist)
    !reads jobs from file or from available xyz/ascii files.
    !job file is limited to 999 lines
    use module_base
    use module_global_variables, only: operation_mode,&
                                       currdir
    implicit none
    !parameters
    integer, intent(out) :: njobs
    character(len=100) :: joblist(2,999)
    !internal
    logical :: exists,fexists
    character(len=300) :: line
    character(len=6) :: filename
    integer :: ifile,iline,u,istat
    inquire(file=trim(adjustl(currdir))//'/job_list',exist=exists)
    njobs=0
    if(exists)then
        u=f_get_free_unit()
        open(unit=u,file=trim(adjustl(currdir))//'/job_list')
        if(trim(adjustl(operation_mode))=='connect'.or.&
                       trim(adjustl(operation_mode))=='guessonly')then
            do iline=1,999
                read(u,'(a)',iostat=istat)line
                if(istat/=0)exit
                read(line,*,iostat=istat)joblist(1,iline),joblist(2,iline)
                if(istat/=0)then
                    call f_err_throw(trim(adjustl(currdir))//&
                                 '/job_file wrong format for connect run')
                endif
                joblist(1,iline)=trim(adjustl(currdir))//'/'//joblist(1,iline)
                joblist(2,iline)=trim(adjustl(currdir))//'/'//joblist(2,iline)
                call check_struct_file_exists(trim(adjustl(joblist(1,iline))))
                call check_struct_file_exists(trim(adjustl(joblist(2,iline))))
                njobs=njobs+1
            enddo
            if(iline==1)then
                call f_err_throw(trim(adjustl(currdir))//&
                                '/job_file empty')
            endif
        else if(trim(adjustl(operation_mode))=='simple'.or.&
                trim(adjustl(operation_mode))=='hessian'.or.&
                trim(adjustl(operation_mode))=='minimize')then
            do iline=1,999
                read(u,'(a)',iostat=istat)line
                if(istat/=0)exit
                read(line,*)joblist(1,iline)
                joblist(1,iline)=trim(adjustl(currdir))//'/'//joblist(1,iline)
                call check_struct_file_exists(trim(adjustl(joblist(1,iline))))
                njobs=njobs+1
            enddo
        endif
        close(u)
    else
        if(trim(adjustl(operation_mode))=='connect'.or.&
                       trim(adjustl(operation_mode))=='guessonly')then
            do ifile=1,998
                write(filename,'(a,i3.3)')'pos',ifile
                joblist(1,ifile)=trim(adjustl(currdir))//'/'//filename
                call check_struct_file_exists(trim(adjustl(joblist(1,ifile))),fexists)
                if(.not.fexists)exit
                write(filename,'(a,i3.3)')'pos',ifile+1
                joblist(2,ifile)=trim(adjustl(currdir))//'/'//filename
                call check_struct_file_exists(trim(adjustl(joblist(2,ifile))),fexists)
                if(.not.fexists)exit
                njobs=njobs+1
            enddo
        else if(trim(adjustl(operation_mode))=='simple'.or.&
                trim(adjustl(operation_mode))=='hessian'.or.&
                trim(adjustl(operation_mode))=='minimize')then
            do ifile=1,999
                write(filename,'(a,i3.3)')'pos',ifile
                joblist(1,ifile)=trim(adjustl(currdir))//'/'//filename
                call check_struct_file_exists(trim(adjustl(joblist(1,ifile))),fexists)
                if(.not.fexists)exit
                njobs=njobs+1
            enddo
        endif
    endif

end subroutine

subroutine check_struct_file_exists(filename,exists)
    use module_base
    implicit none
    !parameter
    character(len=*), intent(in) :: filename
    logical, optional, intent(out) :: exists
    !internal
    logical :: xyzexists=.false.,asciiexists=.false.
    integer :: indx,inda

    if(present(exists))then
        exists=.true.
    endif

    indx=index(filename,'.xyz')
    inda=index(filename,'.ascii')
    if(indx==0)then
        inquire(file=trim(adjustl(filename))//'.xyz',exist=xyzexists)
    else
        inquire(file=trim(adjustl(filename)),exist=xyzexists)
    endif
    if(inda==0)then
        inquire(file=trim(adjustl(filename))//'.ascii',exist=asciiexists)
    else
        inquire(file=trim(adjustl(filename)),exist=asciiexists)
    endif
    if(.not. (xyzexists .or. asciiexists))then
        if(present(exists))then
            exists=.false.
        else
        call f_err_throw('File '//trim(adjustl(filename))//&
                        ' does not exist.')
        endif
    endif
end subroutine

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
    if(nat/=astruct%nat) &
         call f_err_throw('(MHGPS) severe error in read_mode: '//&
         'nat/=astruct%nat')
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
    use module_atoms, only: astruct_dump_to_file
    use module_global_variables, only: iproc, astruct=>astruct_ptr
    implicit none
    !parameters
    integer, intent(in) :: nat
    character(len=*), intent(in) :: filename
    real(gp), intent(in) :: minmode(3,nat)
    real(gp), intent(in), optional :: rotforce(3,nat)
    !internal
    character(len=7) :: comment='minmode'
    character(len=11) :: units

    units=astruct%units
    astruct%units='atomicd0'
    if(present(rotforce))then
       call astruct_dump_to_file(astruct,filename,trim(comment),&
            rxyz=minmode,forces=rotforce)
!!$
!!$        call write_atomic_file(filename,&
!!$              0.0_gp,minmode(1,1),ixyz_int,&
!!$              atoms,trim(comment),forces=rotforce(1,1))
    else
       call astruct_dump_to_file(astruct,filename,trim(comment),&
            rxyz=minmode)
!!$
!!$        call write_atomic_file(filename,&
!!$              0.0_gp,minmode(1,1),ixyz_int,&
!!$              atoms,trim(comment))
    endif
    astruct%units=units
end subroutine

end module
