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
    public :: check_struct_file_exists

contains
!=====================================================================
subroutine get_first_struct_file(mhgpsst,filename)
    use module_base
    use module_mhgps_state
    implicit none
    !parameter
    type(mhgps_state), intent(in) :: mhgpsst
    character(len=200), intent(out) :: filename
    !local
    integer :: istat,u
    character(len=33) :: line
    logical :: exists


    inquire(file=trim(adjustl(mhgpsst%currDir))//'/job_list',exist=exists)
    if(exists)then
        u=f_get_free_unit()
        open(unit=u,file=trim(adjustl(mhgpsst%currDir))//'/job_list')
        read(u,'(a)',iostat=istat)line
        if(istat/=0)then
            call f_err_throw(trim(adjustl(mhgpsst%currDir))//&
                                '/job_file empty')
        endif
        read(line,*)filename
        filename=trim(adjustl(mhgpsst%currDir))//'/'//filename
        call check_struct_file_exists(filename)
        close(u)
    else
        write(filename,'(a,i3.3)')trim(adjustl(mhgpsst%currDir))//'/pos',1
        call check_struct_file_exists(filename)
    endif
end subroutine
!=====================================================================
subroutine read_restart(mhgpsst,runObj)
    use module_base
    use bigdft_run
    use module_mhgps_state
    implicit none
    !parameters
    type(mhgps_state), intent(inout) :: mhgpsst
    type(run_objects), intent(in)     :: runObj
    !local
    integer :: u
    logical :: exists
    integer :: iatt, iat
    u=f_get_free_unit()
    inquire(file='restart',exist=exists)
    if(exists)then
        open(unit=u,file='restart')
        read(u,*)mhgpsst%ifolder
        read(u,*)mhgpsst%isad,mhgpsst%isadprob
        read(u,*)mhgpsst%nsad
        read(u,*)mhgpsst%ntodo
        read(u,*)mhgpsst%nrestart
        mhgpsst%nrestart=mhgpsst%nrestart+1
        read(u,*)mhgpsst%nattempted,mhgpsst%nattemptedmax
        mhgpsst%attempted_connections =&
             f_malloc([3,runObj%atoms%astruct%nat,2,mhgpsst%nattemptedmax]&
             ,id='mhgpsst%attempted_connections')
        do iatt=1,mhgpsst%nattempted
            do iat=1,runObj%atoms%astruct%nat
            read(u,*)mhgpsst%attempted_connections(1,iat,1,iatt),&
                     mhgpsst%attempted_connections(2,iat,1,iatt),&
                     mhgpsst%attempted_connections(3,iat,1,iatt)
            enddo
            do iat=1,runObj%atoms%astruct%nat
            read(u,*)mhgpsst%attempted_connections(1,iat,2,iatt),&
                     mhgpsst%attempted_connections(2,iat,2,iatt),&
                     mhgpsst%attempted_connections(3,iat,2,iatt)
            enddo
        enddo
        close(u)
    else
        mhgpsst%ifolder=1
        mhgpsst%isad=0
        mhgpsst%isadprob=0
        mhgpsst%ntodo=0
        mhgpsst%nsad=0
        mhgpsst%nrestart=0
        mhgpsst%nattempted=0
        mhgpsst%nattemptedmax=1000
        mhgpsst%attempted_connections =&
             f_malloc([3,runObj%atoms%astruct%nat,2,mhgpsst%nattemptedmax]&
             ,id='mhgpsst%attempted_connections')
    endif
end subroutine read_restart
!=====================================================================
subroutine write_restart(mhgpsst,runObj,cobj,writeJobList)
    use module_base
    use bigdft_run
    use module_mhgps_state
    use module_connect_object
    implicit none
    type(mhgps_state), intent(inout)              :: mhgpsst
    type(run_objects), intent(in)     :: runObj
    type(connect_object), optional, intent(in) :: cobj
    logical, optional, intent(in) :: writeJobList
    !local
    integer :: u
    integer :: iatt, iat
    logical :: wJl
    wJl=.true.
    u=f_get_free_unit()
    open(unit=u,file='restart')
    write(u,*)mhgpsst%ifolder
    write(u,*)mhgpsst%isad,mhgpsst%isadprob
    write(u,*)mhgpsst%nsad
    write(u,*)mhgpsst%ntodo
    write(u,*)mhgpsst%nrestart
    write(u,*)mhgpsst%nattempted, mhgpsst%nattemptedmax
    do iatt=1,mhgpsst%nattempted
        do iat=1,runObj%atoms%astruct%nat
        write(u,'(3(1x,es24.17))')mhgpsst%attempted_connections(1,iat,1,iatt),&
                                mhgpsst%attempted_connections(2,iat,1,iatt),&
                                mhgpsst%attempted_connections(3,iat,1,iatt)
        enddo
        do iat=1,runObj%atoms%astruct%nat
        write(u,'(3(1x,es24.17))')mhgpsst%attempted_connections(1,iat,2,iatt),&
                                mhgpsst%attempted_connections(2,iat,2,iatt),&
                                mhgpsst%attempted_connections(3,iat,2,iatt)
        enddo
    enddo
    close(u)
    
    if(present(writeJobList)) wJl=writeJobList
    if(wJl)then
        call write_jobs(mhgpsst,runObj,cobj)
    endif
end subroutine write_restart
!=====================================================================
subroutine write_jobs(mhgpsst,runObj,cobj)
    use module_base
    use bigdft_run
    use module_atoms, only: astruct_dump_to_file
    use module_mhgps_state
    use module_connect_object
    implicit none
    !parameters
    type(mhgps_state), intent(inout) :: mhgpsst
    type(run_objects), intent(in)     :: runObj
    type(connect_object), optional, intent(in) :: cobj
    !local
    integer :: ijob, u
    character(len=1) :: comment
    character(len=21)  :: filenameR, filenameL
    logical :: lw
    comment = ' ' 
    
    lw=.false.

    if(mhgpsst%ijob+1<=mhgpsst%njobs) lw=.true.
    if(present(cobj))then
        if(cobj%ntodo>=1) lw=.true.
    endif
    if(lw)then
        u=f_get_free_unit()
        open(unit=u,file=trim(adjustl(mhgpsst%currdir))//'/job_list_restart')
            if(present(cobj))then
                do ijob=cobj%ntodo,1,-1
                    !write files
                    write(filenameR,'(a,i5.5,a,i5.5,a)')'restart_',&
                          mhgpsst%nrestart,'_',ijob,'_R'
                    call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
                         trim(adjustl(mhgpsst%currDir))//'/'//trim(filenameR),&
                         comment,rxyz=cobj%todorxyz(1,1,2,ijob))
                    write(filenameL,'(a,i5.5,a,i5.5,a)')'restart_',&
                          mhgpsst%nrestart,'_',ijob,'_L'
                    call astruct_dump_to_file(bigdft_get_astruct_ptr(runObj),&
                         trim(adjustl(mhgpsst%currDir))//'/'//trim(filenameL),&
                         comment,rxyz=cobj%todorxyz(1,1,1,ijob))
                    !write job file
                    write(u,'(a,1x,a)')trim(adjustl(filenameL)),trim(adjustl(filenameR))
                enddo
                do ijob=mhgpsst%ijob+1,mhgpsst%njobs
                    write(u,'(a,1x,a)')trim(adjustl(mhgpsst%joblist(1,ijob)(10:))),&
                                       trim(adjustl(mhgpsst%joblist(2,ijob)(10:)))
                enddo
            else
                do ijob=mhgpsst%ijob+1,mhgpsst%njobs
                    if(trim(adjustl(mhgpsst%joblist(1,ijob)(10:16)))/='restart')then
                    write(u,'(a,1x,a)')trim(adjustl(mhgpsst%joblist(1,ijob)(10:))),&
                                       trim(adjustl(mhgpsst%joblist(2,ijob)(10:)))
                    endif
                enddo
            endif
        close(u)
    endif
     
end subroutine write_jobs
!=====================================================================
subroutine read_jobs(uinp,mhgpsst)
    !reads jobs from file or from available xyz/ascii files.
    !job file is limited to 999 lines
    use module_base
    use module_userinput, only: userinput
    use module_mhgps_state
    implicit none
    !parameters
    type(userinput), intent(in) :: uinp
    type(mhgps_state), intent(inout) :: mhgpsst
    !local
    logical :: exists,fexists,exists_restart
    character(len=300) :: line
    character(len=6) :: filename
    integer :: ifile,iline,u,istat
    character(len=50) :: jobfile
    !put a barrier for all the processes
    !otherwise, the joblist file may be updated, before
    !all procsesses have read the identical file
    call mpibarrier(bigdft_mpi%mpi_comm)

    inquire(file=trim(adjustl(mhgpsst%currDir))//'/job_list',exist=exists)
    if(exists)jobfile=trim(adjustl(mhgpsst%currDir))//'/'//'job_list'
    inquire(file=trim(adjustl(mhgpsst%currDir))//'/job_list_restart',&
                 exist=exists_restart)
    if(exists_restart)jobfile=trim(adjustl(mhgpsst%currDir))//'/'//&
                      'job_list_restart'
    mhgpsst%njobs=0
    iline=1
    if(exists .or. exists_restart)then
        u=f_get_free_unit()
        open(unit=u,file=trim(adjustl(jobfile)))
        if(trim(adjustl(uinp%operation_mode))=='connect'.or.&
                       trim(adjustl(uinp%operation_mode))=='guessonly')then
            do iline=1,999
                read(u,'(a)',iostat=istat)line
                if(istat/=0)exit
                read(line,*,iostat=istat)mhgpsst%joblist(1,iline),&
                                         mhgpsst%joblist(2,iline)
                if(istat/=0)then
                    call f_err_throw(trim(adjustl(jobfile))//&
                         ' wrong format for connect run')
                endif
                mhgpsst%joblist(1,iline)=trim(adjustl(mhgpsst%currDir))//'/'//&
                                 mhgpsst%joblist(1,iline)
                mhgpsst%joblist(2,iline)=trim(adjustl(mhgpsst%currDir))//'/'//&
                                 mhgpsst%joblist(2,iline)
                call check_struct_file_exists(&
                     trim(adjustl(mhgpsst%joblist(1,iline))))
                call check_struct_file_exists(&
                     trim(adjustl(mhgpsst%joblist(2,iline))))
                mhgpsst%njobs=mhgpsst%njobs+1
            enddo
            if(iline==1)then
                call f_err_throw(trim(adjustl(jobfile))//' empty')
            endif
        else if(trim(adjustl(uinp%operation_mode))=='simple'.or.&
                trim(adjustl(uinp%operation_mode))=='simpleandminimize'.or.&
                trim(adjustl(uinp%operation_mode))=='hessian'.or.&
                trim(adjustl(uinp%operation_mode))=='minimize')then
            do iline=1,999
                read(u,'(a)',iostat=istat)line
                if(istat/=0)exit
                read(line,*)mhgpsst%joblist(1,iline)
                mhgpsst%joblist(1,iline)=trim(adjustl(mhgpsst%currDir))//'/'//&
                                 mhgpsst%joblist(1,iline)
                mhgpsst%joblist(2,iline)=''
                call check_struct_file_exists(&
                     trim(adjustl(mhgpsst%joblist(1,iline))))
                mhgpsst%njobs=mhgpsst%njobs+1
            enddo
        else
            call f_err_throw('Operationmode '//&
                 trim(adjustl(uinp%operation_mode))//' unknown.')
        endif
        close(u)
    else
        if(trim(adjustl(uinp%operation_mode))=='connect'.or.&
                       trim(adjustl(uinp%operation_mode))=='guessonly')then
            do ifile=1,998
                write(filename,'(a,i3.3)')'pos',ifile
                mhgpsst%joblist(1,ifile)=trim(adjustl(mhgpsst%currDir))//'/'//filename
                call check_struct_file_exists(&
                     trim(adjustl(mhgpsst%joblist(1,ifile))),fexists)
                if(.not.fexists)exit
                write(filename,'(a,i3.3)')'pos',ifile+1
                mhgpsst%joblist(2,ifile)=trim(adjustl(mhgpsst%currDir))//'/'//filename
                call check_struct_file_exists(&
                     trim(adjustl(mhgpsst%joblist(2,ifile))),fexists)
                if(.not.fexists)exit
                mhgpsst%njobs=mhgpsst%njobs+1
            enddo
        else if(trim(adjustl(uinp%operation_mode))=='simple'.or.&
                trim(adjustl(uinp%operation_mode))=='simpleandminimize'.or.&
                trim(adjustl(uinp%operation_mode))=='hessian'.or.&
                trim(adjustl(uinp%operation_mode))=='minimize')then
            do ifile=1,999
                write(filename,'(a,i3.3)')'pos',ifile
                mhgpsst%joblist(1,ifile)=trim(adjustl(mhgpsst%currDir))//'/'//filename
                mhgpsst%joblist(2,iline)=''
                call check_struct_file_exists(&
                     trim(adjustl(mhgpsst%joblist(1,ifile))),fexists)
                if(.not.fexists)exit
                mhgpsst%njobs=mhgpsst%njobs+1
            enddo
        else
            call f_err_throw('Operationmode '//&
                 trim(adjustl(uinp%operation_mode))//' unknown.')
        endif
    endif

    !put a barrier for all the processes
    !otherwise, the joblist file may be updated, before
    !all procsesses have read the identical file
    call mpibarrier(bigdft_mpi%mpi_comm)
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
subroutine read_mode(mhgpsst,nat,filename,minmode)
    use module_base
    use module_types
    use module_atoms, only: atomic_structure,&
                            deallocate_atomic_structure,&
                            read_atomic_file=>set_astruct_from_file
    use module_mhgps_state
    implicit none
    !parameters
    type(mhgps_state), intent(in) :: mhgpsst
    integer, intent(in) :: nat
    character(len=*), intent(in) :: filename
    real(gp), intent(inout) :: minmode(3,nat)
    !local
    type(atomic_structure):: astruct !< Contains all info


    call read_atomic_file(filename,mhgpsst%iproc,astruct,&
         disableTrans=.true.)
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
subroutine write_mode(runObj,outs,filename,minmode,rotforce)
    use module_base, only: gp
    use module_types
    use module_interfaces
    use module_atoms, only: astruct_dump_to_file
    use bigdft_run
    implicit none
    !parameters
    type(run_objects), intent(inout) :: runObj
    type(state_properties), intent(inout) :: outs
    character(len=*), intent(in) :: filename
    real(gp), intent(in) :: minmode(3,runObj%atoms%astruct%nat)
    real(gp), intent(in), optional :: rotforce(3,runObj%atoms%astruct%nat)
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
subroutine print_logo_mhgps(mhgpsst)
    use yaml_output
    use module_mhgps_state
    implicit none
    !parameters
    type(mhgps_state), intent(in) :: mhgpsst

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
    call yaml_mapping_open('(MHGPS) Reference Papers')
    call yaml_scalar('(MHGPS) The Journal of Chemical Physics 140, 214102 (2014)')
    call yaml_scalar('(MHGPS) The Journal of Chemical Physics 142, 034112 (2015)')
    call yaml_mapping_close()
    call yaml_map('(MHGPS) Version Number',trim(adjustl(mhgpsst%mhgps_version)))
    call yaml_map('(MHGPS) Timestamp of this run',yaml_date_and_time_toa())
end subroutine print_logo_mhgps


end module
