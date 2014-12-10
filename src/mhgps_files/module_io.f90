!> @file
!!  Input/Output for minima hopping
!!
!! @author
!!    Copyright (C) 2013-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module handling input/output for minima hopping guided
!! path sampling
module module_io

    implicit none

    private

    public :: read_mode
    public :: write_mode
    public :: read_jobs
    public :: write_jobs
    public :: read_restart
    public :: write_restart
    public :: get_first_struct_file
    public :: print_logo_mhgps

contains
!=====================================================================
subroutine get_first_struct_file(filename)
    use module_base
    use module_global_variables, only: currdir
    implicit none
    !parameter
    character(len=200), intent(out) :: filename
    !local
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
!=====================================================================
subroutine write_jobs(njobs,joblist,idone)
    use module_base
    use module_global_variables, only: currdir
    implicit none
    !parameters
    integer, intent(in) :: njobs
    character(len=100)   :: joblist(2,999)
    integer, intent(in)  :: idone
    !local
    integer :: ijob, u

    u=f_get_free_unit()
    open(unit=u,file=trim(adjustl(currdir))//'/job_list_restart')
    do ijob=idone+1,njobs
        write(u,'(a,1x,a)')trim(adjustl(joblist(1,ijob)(10:))),&
                           trim(adjustl(joblist(2,ijob)(10:)))
    enddo
    close(u)
     
end subroutine write_jobs
!=====================================================================
subroutine read_restart(isad,isadprob,ntodo)
    use module_base
    implicit none
    !parameters
    integer, intent(out) :: isad, isadprob, ntodo
    !local
    integer :: u
    logical :: exists
    u=f_get_free_unit()
    inquire(file='restart',exist=exists)
    if(exists)then
        open(unit=u,file='restart')
        read(u,*)isad,isadprob
        read(u,*)ntodo
        close(u)
    else
        isad=0
        ntodo=0
    endif
end subroutine read_restart
!=====================================================================
subroutine write_restart(isad,isadprob,ntodo)
    use module_base
    implicit none
    integer, intent(in) :: isad, isadprob, ntodo
    !local
    integer :: u
    u=f_get_free_unit()
    open(unit=u,file='restart')
    write(u,*)isad,isadprob
    write(u,*)ntodo
    close(u)
end subroutine write_restart
!=====================================================================
subroutine read_jobs(uinp,njobs,joblist)
    !reads jobs from file or from available xyz/ascii files.
    !job file is limited to 999 lines
    use module_base
    use module_userinput, only: userinput
    use module_global_variables, only: currdir
    implicit none
    !parameters
    type(userinput), intent(in) :: uinp
    integer, intent(out) :: njobs
    character(len=100) :: joblist(2,999)
    !local
    logical :: exists,fexists,exists_restart
    character(len=300) :: line
    character(len=6) :: filename
    integer :: ifile,iline,u,istat
    character(len=50) :: jobfile
    inquire(file=trim(adjustl(currdir))//'/job_list',exist=exists)
    if(exists)jobfile=trim(adjustl(currdir))//'/'//'job_list'
    inquire(file=trim(adjustl(currdir))//'/job_list_restart',&
                 exist=exists_restart)
    if(exists_restart)jobfile=trim(adjustl(currdir))//'/'//&
                      'job_list_restart'
    njobs=0
    if(exists .or. exists_restart)then
        u=f_get_free_unit()
        open(unit=u,file=trim(adjustl(jobfile)))
        if(trim(adjustl(uinp%operation_mode))=='connect'.or.&
                       trim(adjustl(uinp%operation_mode))=='guessonly')then
            do iline=1,999
                read(u,'(a)',iostat=istat)line
                if(istat/=0)exit
                read(line,*,iostat=istat)joblist(1,iline),&
                                         joblist(2,iline)
                if(istat/=0)then
                    call f_err_throw(trim(adjustl(jobfile))//&
                         ' wrong format for connect run')
                endif
                joblist(1,iline)=trim(adjustl(currdir))//'/'//&
                                 joblist(1,iline)
                joblist(2,iline)=trim(adjustl(currdir))//'/'//&
                                 joblist(2,iline)
                call check_struct_file_exists(&
                     trim(adjustl(joblist(1,iline))))
                call check_struct_file_exists(&
                     trim(adjustl(joblist(2,iline))))
                njobs=njobs+1
            enddo
            if(iline==1)then
                call f_err_throw(trim(adjustl(jobfile))//' empty')
            endif
        else if(trim(adjustl(uinp%operation_mode))=='simple'.or.&
                trim(adjustl(uinp%operation_mode))=='hessian'.or.&
                trim(adjustl(uinp%operation_mode))=='minimize')then
            do iline=1,999
                read(u,'(a)',iostat=istat)line
                if(istat/=0)exit
                read(line,*)joblist(1,iline)
                joblist(1,iline)=trim(adjustl(currdir))//'/'//&
                                 joblist(1,iline)
                joblist(2,iline)=''
                call check_struct_file_exists(&
                     trim(adjustl(joblist(1,iline))))
                njobs=njobs+1
            enddo
        endif
        close(u)
    else
        if(trim(adjustl(uinp%operation_mode))=='connect'.or.&
                       trim(adjustl(uinp%operation_mode))=='guessonly')then
            do ifile=1,998
                write(filename,'(a,i3.3)')'pos',ifile
                joblist(1,ifile)=trim(adjustl(currdir))//'/'//filename
                call check_struct_file_exists(&
                     trim(adjustl(joblist(1,ifile))),fexists)
                if(.not.fexists)exit
                write(filename,'(a,i3.3)')'pos',ifile+1
                joblist(2,ifile)=trim(adjustl(currdir))//'/'//filename
                call check_struct_file_exists(&
                     trim(adjustl(joblist(2,ifile))),fexists)
                if(.not.fexists)exit
                njobs=njobs+1
            enddo
        else if(trim(adjustl(uinp%operation_mode))=='simple'.or.&
                trim(adjustl(uinp%operation_mode))=='hessian'.or.&
                trim(adjustl(uinp%operation_mode))=='minimize')then
            do ifile=1,999
                write(filename,'(a,i3.3)')'pos',ifile
                joblist(1,ifile)=trim(adjustl(currdir))//'/'//filename
                joblist(2,iline)=''
                call check_struct_file_exists(&
                     trim(adjustl(joblist(1,ifile))),fexists)
                if(.not.fexists)exit
                njobs=njobs+1
            enddo
        endif
    endif

end subroutine
!=====================================================================
subroutine check_struct_file_exists(filename,exists)
    use module_base
    implicit none
    !parameter
    character(len=*), intent(in) :: filename
    logical, optional, intent(out) :: exists
    !local
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
        inquire(file=trim(adjustl(filename))//'.ascii',&
                exist=asciiexists)
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
!=====================================================================
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
    !local
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
!=====================================================================
subroutine write_mode(nat,runObj,outs,filename,minmode,rotforce)
    use module_base, only: gp
    use module_types
    use module_interfaces
    use module_atoms, only: astruct_dump_to_file
    use module_global_variables, only: iproc
    use bigdft_run
    implicit none
    !parameters
    integer, intent(in) :: nat
    type(run_objects), intent(inout) :: runObj
    type(state_properties), intent(inout) :: outs
    character(len=*), intent(in) :: filename
    real(gp), intent(in) :: minmode(3,nat)
    real(gp), intent(in), optional :: rotforce(3,nat)
    !local
    character(len=7) :: comment='minmode'
    character(len=20) :: units

    units=bigdft_get_units(runObj)
    call bigdft_set_units(runObj,'atomicd0')
    if(present(rotforce))then
       call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
            filename,trim(comment),rxyz=minmode,forces=rotforce)
    else
       call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
            filename,trim(comment),rxyz=minmode)
    endif
    call bigdft_set_units(runObj,units)
end subroutine
!=====================================================================
subroutine print_logo_mhgps()
    use module_global_variables, only: mhgps_version
    use yaml_output
    implicit none

    call yaml_comment('(MHGPS) Minima Hopping Guided Path Sampling',hfill='=')
    
    call yaml_mapping_open('(MHGPS) logo')
    call yaml_scalar('(MHGPS)      ___           ___           ___           ___           ___      ') 
    call yaml_scalar('(MHGPS)     /\__\         /\__\         /\  \         /\  \         /\  \     ')
    call yaml_scalar('(MHGPS)    /::|  |       /:/  /        /::\  \       /::\  \       /::\  \    ')
    call yaml_scalar('(MHGPS)   /:|:|  |      /:/__/        /:/\:\  \     /:/\:\  \     /:/\ \  \   ')
    call yaml_scalar('(MHGPS)  /:/|:|__|__   /::\  \ ___   /:/  \:\  \   /::\~\:\  \   _\:\~\ \  \  ')
    call yaml_scalar('(MHGPS) /:/ |::::\__\ /:/\:\  /\__\ /:/__/_\:\__\ /:/\:\ \:\__\ /\ \:\ \ \__\ ')
    call yaml_scalar('(MHGPS) \/__/~~/:/  / \/__\:\/:/  / \:\  /\ \/__/ \/__\:\/:/  / \:\ \:\ \/__/ ')
    call yaml_scalar('(MHGPS)       /:/  /       \::/  /   \:\ \:\__\        \::/  /   \:\ \:\__\   ')
    call yaml_scalar('(MHGPS)      /:/  /        /:/  /     \:\/:/  /         \/__/     \:\/:/  /   ')
    call yaml_scalar('(MHGPS)     /:/  /        /:/  /       \::/  /                     \::/  /    ')
    call yaml_scalar('(MHGPS)     \/__/         \/__/         \/__/                       \/__/     ')
    call yaml_scalar('(MHGPS)                                                   as post-processing  ')
    call yaml_scalar('(MHGPS)')
    call yaml_scalar('(MHGPS)')
    !call print_logo()
    call yaml_mapping_close()
    call yaml_mapping_open('(MHGPS) Reference Paper')
    call yaml_scalar('(MHGPS) The Journal of Chemical Physics 140 (21):214102 (2014)')
    call yaml_mapping_close()
    call yaml_map('(MHGPS) Version Number',mhgps_version)
    call yaml_map('(MHGPS) Timestamp of this run',yaml_date_and_time_toa())
end subroutine print_logo_mhgps


end module
