!> @file
!!  Routines to initialize some important structures for run_object
!! @author
!!    Copyright (C) 2007-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> Read the options in the command line using get_command statement
subroutine command_line_information(mpi_groupsize,posinp_file,run_id,ierr)
  use module_types
  implicit none
  integer, intent(out) :: mpi_groupsize
  character(len=*), intent(out) :: posinp_file !< file for list of radicals
  character(len=*), intent(out) :: run_id !< file for radical name
  integer, intent(out) :: ierr !< error code
  !local variables
  integer :: ncommands,icommands
  character(len=256) :: command

  ierr=BIGDFT_SUCCESS
  posinp_file=repeat(' ',len(posinp_file))
  run_id=repeat(' ',len(run_id))
  !traditional scheme
  !if (ncommands == 0) then
     run_id='input'
  !end if

  mpi_groupsize=0
  
  !first see how many arguments are present
  ncommands=COMMAND_ARGUMENT_COUNT()

  do icommands=1,ncommands
     command=repeat(' ',len(command))
     call get_command_argument(icommands,value=command,status=ierr)
     if (ierr /= 0) return
     call find_command()
     if (ierr /= 0) return
  end do

contains

  subroutine find_command()
    implicit none
    integer :: ipos
    integer, external :: bigdft_error_ret

    if (index(command,'--taskgroup-size=') > 0) then
       if (mpi_groupsize /= 0) then
          ierr=bigdft_error_ret(BIGDFT_INVALID,'taskgroup size specified twice')
       end if
       ipos=index(command,'=')
       read(command(ipos+1:len(command)),*)mpi_groupsize
    else if (index(command,'--run-id=') > 0) then
       if (trim(run_id) /= 'input') then
          ierr=bigdft_error_ret(BIGDFT_INVALID,'run_id specified twice')
       end if
       ipos=index(command,'=')
       read(command(ipos+1:len(command)),*)run_id
    else if (index(command,'--runs-file=') > 0) then
       if (len_trim(posinp_file) > 0 .or. trim(run_id) /= 'input') then
          ierr=bigdft_error_ret(BIGDFT_INVALID,'posinp_file specified twice or run_id already known')
       end if
       ipos=index(command,'=')
       read(command(ipos+1:len(command)),*)posinp_file
    else if (index(command,'--') > 0 .and. icommands==1) then
       !help screen
       call help_screen()
       stop
    else if (icommands==1) then
       read(command,*,iostat=ierr)run_id
    else
       call help_screen()
       stop
    end if
  end subroutine find_command

  subroutine help_screen()
    write(*,*)' Usage of the command line instruction'
    write(*,*)' --taskgroup-size=<mpi_groupsize>'
    write(*,*)' --runs-file=<list_posinp filename>'
    write(*,*)' --run-id=<name of the run>: it can be also specified as unique argument'
    write(*,*)' --help : prints this help screen'
  end subroutine help_screen

END SUBROUTINE command_line_information

!> Initialization of acceleration (OpenCL)
subroutine init_material_acceleration(iproc,matacc,GPU)
  use module_base
  use module_types
  use yaml_output
  implicit none
  integer, intent(in):: iproc
  type(material_acceleration), intent(in) :: matacc
  type(GPU_pointers), intent(out) :: GPU
  !local variables
  integer :: iconv,iblas,initerror,ierror,useGPU,mproc,ierr,nproc_node
  logical :: noaccel
  integer(kind=8) :: context_address

  noaccel = .true.
  if (matacc%iacceleration == 1) then
     call MPI_COMM_SIZE(bigdft_mpi%mpi_comm,mproc,ierr)
     !initialize the id_proc per node
     call processor_id_per_node(iproc,mproc,GPU%id_proc,nproc_node)
     call sg_init(GPUshare,useGPU,iproc,nproc_node,initerror)
     if (useGPU == 1) then
        iconv = 1
        iblas = 1
     else
        iconv = 0
        iblas = 0
     end if
     if (initerror == 1) then
        call yaml_warning('(iproc=' // trim(yaml_toa(iproc,fmt='(i0)')) // &
        &    ') S_GPU library init failed, aborting...')
        !write(*,'(1x,a)')'**** ERROR: S_GPU library init failed, aborting...'
        call MPI_ABORT(bigdft_mpi%mpi_comm,initerror,ierror)
     end if

     if (iconv == 1) then
        !change the value of the GPU convolution flag defined in the module_base
        GPUconv=.true.
     end if
     if (iblas == 1) then
        !change the value of the GPU convolution flag defined in the module_base
        GPUblas=.true.
     end if

     if (iproc == 0) then
        call yaml_map('Material acceleration','CUDA',advance='no')
        call yaml_comment('iproc=0')
       ! write(*,'(1x,a)') 'CUDA support activated (iproc=0)'
    end if

    noaccel = .false.
  else if (matacc%iacceleration >= 2) then
     ! OpenCL convolutions are activated
     ! use CUBLAS for the linear algebra for the moment
     call MPI_COMM_SIZE(bigdft_mpi%mpi_comm,mproc,ierr)
     !initialize the id_proc per node
     call processor_id_per_node(iproc,mproc,GPU%id_proc,nproc_node)
     !initialize the opencl context for any process in the node
     !call MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)
     !do jproc=0,mproc-1
     !   call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
     !   if (iproc == jproc) then
     !      print '(a,a,i4,i4)','Initializing for node: ',trim(nodename_local),iproc,GPU%id_proc
     call init_acceleration_OCL(matacc,GPU)
     !   end if
     !end do
     !to avoid a representation of the address which is lower than tiny(1.d0)
     context_address=transfer(GPU%context,context_address)
     !if (GPU%context /= 0.) then
     if (context_address /= int(0,kind=8)) then
        GPU%ndevices=min(GPU%ndevices,nproc_node)
        if (iproc == 0) then
           call yaml_map('Material acceleration','OpenCL',advance='no')
           call yaml_comment('iproc=0')
           call yaml_mapping_open('Number of OpenCL devices per node',flow=.true.)
           call yaml_map('used',trim(yaml_toa(min(GPU%ndevices,nproc_node),fmt='(i0)')))
           call yaml_map('available',trim(yaml_toa(GPU%ndevices,fmt='(i0)')))
           !write(*,'(1x,a,i5,i5)') 'OpenCL support activated, No. devices per node (used, available):',&
           !     min(GPU%ndevices,nproc_node),GPU%ndevices
           call yaml_mapping_close()
        end if
        !the number of devices is the min between the number of processes per node
        GPU%ndevices=min(GPU%ndevices,nproc_node)
        GPU%OCLconv=.true.
        noaccel = .false.
     end if
  end if

  if (noaccel .and. iproc == 0) then
     call yaml_map('Material acceleration',.false.,advance='no')
     call yaml_comment('iproc=0')
     ! write(*,'(1x,a)') 'No material acceleration (iproc=0)'
  end if

END SUBROUTINE init_material_acceleration


subroutine release_material_acceleration(GPU)
  use module_base
  use module_types
  implicit none
  type(GPU_pointers), intent(inout) :: GPU
  
  if (GPUconv) then
     call sg_end()
  end if

  if (GPU%OCLconv) then
     call release_acceleration_OCL(GPU)
     GPU%OCLconv=.false.
  end if

END SUBROUTINE release_material_acceleration


!> Give the number of MPI processes per node (nproc_node) and before iproc (iproc_node)
subroutine processor_id_per_node(iproc,nproc,iproc_node,nproc_node)
  use module_base
  use module_types
  use dynamic_memory
  implicit none
  integer, intent(in) :: iproc, nproc
  integer, intent(out) :: iproc_node, nproc_node
  !local variables
  character(len=*), parameter :: subname='processor_id_per_node'
  integer :: ierr,namelen,jproc
  character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
  character(len=MPI_MAX_PROCESSOR_NAME), dimension(:), allocatable :: nodename

  call f_routine(id=subname)

  if (nproc == 1) then
     iproc_node=0
     nproc_node=1
  else
     nodename=f_malloc_str(MPI_MAX_PROCESSOR_NAME,0 .to. nproc-1,id='nodename')
     !allocate(nodename(0:nproc-1+ndebug),stat=i_stat)
     !call memocc(i_stat,nodename,'nodename',subname)
     
     !initalise nodenames
     do jproc=0,nproc-1
        nodename(jproc)=repeat(' ',MPI_MAX_PROCESSOR_NAME)
     end do

     call MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)

     !gather the result between all the process
     call MPI_ALLGATHER(nodename_local,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
          nodename(0),MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
          bigdft_mpi%mpi_comm,ierr)

     !found the processors which belong to the same node
     !before the processor iproc
     iproc_node=0
     do jproc=0,iproc-1
        if (trim(nodename(jproc)) == trim(nodename(iproc))) then
           iproc_node=iproc_node+1
        end if
     end do
     nproc_node=iproc_node
     do jproc=iproc,nproc-1
        if (trim(nodename(jproc)) == trim(nodename(iproc))) then
           nproc_node=nproc_node+1
        end if
     end do
     
     call f_free_str(MPI_MAX_PROCESSOR_NAME,nodename)
     !i_all=-product(shape(nodename))*kind(nodename)
     !deallocate(nodename,stat=i_stat)
     !call memocc(i_stat,i_all,'nodename',subname)
  end if
  call f_release_routine()
END SUBROUTINE processor_id_per_node

subroutine ensure_log_file(writing_directory, logfile, ierr)
  use yaml_output
  use yaml_strings
  implicit none
  character(len = *), intent(in) :: writing_directory, logfile
  integer, intent(out) :: ierr

  logical :: exists
  integer :: lgt
  character(len = 500) :: logfile_old, logfile_dir, filepath

  ierr = 0
  filepath = writing_directory//logfile
  !inquire for the existence of a logfile
  inquire(file=trim(filepath),exist=exists)
  if (exists) then
     logfile_old=writing_directory//'logfiles'
     call getdir(logfile_old,&
          len_trim(logfile_old),logfile_dir,len(logfile_dir),ierr)
     if (ierr /= 0) then
        write(*,*) "ERROR: cannot create writing directory '" //trim(logfile_dir) // "'."
        return
     end if
     logfile_old=trim(logfile_dir)//logfile
     !change the name of the existing logfile
     lgt=index(logfile_old,'.yaml')
     call buffer_string(logfile_old,len(logfile_old),&
          trim(adjustl(yaml_time_toa()))//'.yaml',lgt)
     call movefile(trim(filepath),len_trim(filepath),trim(logfile_old),len_trim(logfile_old),ierr)
     if (ierr /= 0) then
        write(*,*) "ERROR: cannot move logfile '"//trim(logfile)
        write(*,*) '                      into '//trim(logfile_old)// "'."
        return
     end if
     call yaml_map('<BigDFT> Logfile existing, renamed into',&
          trim(logfile_old),unit=6)
  end if

end subroutine ensure_log_file

subroutine create_log_file(dict)

  use module_base
  use module_types
  use module_input
  use yaml_strings
  use yaml_output
  use dictionaries
  use bigdft_run, only: bigdft_get_run_properties, bigdft_set_run_properties
  implicit none
  type(dictionary), pointer :: dict
  !local variables
  integer :: ierr,ierror,lgt,unit_log
  character(len = max_field_length) :: writing_directory, run_name
  character(len=500) :: logfilename,path
  integer :: iproc_node, nproc_node
  logical :: log_to_disk

  ! Get user input writing_directory.
  writing_directory = "."
  call bigdft_get_run_properties(dict, outdir_id = writing_directory)
  ! Create writing_directory and parents if needed and broadcast everything.
  if (trim(writing_directory) /= '.') then
     path=repeat(' ',len(path))
     !add the output directory in the directory name
     if (bigdft_mpi%iproc == 0 .and. trim(writing_directory) /= '.') then
        call getdir(writing_directory,&
             len_trim(writing_directory),path,len(path),ierr)
        if (ierr /= 0) then
           write(*,*) "ERROR: cannot create writing directory '"&
                //trim(writing_directory) // "'."
           call MPI_ABORT(bigdft_mpi%mpi_comm,ierror,ierr)
        end if
     end if
     call MPI_BCAST(path,len(path),MPI_CHARACTER,0,bigdft_mpi%mpi_comm,ierr)
     lgt=min(len(writing_directory),len(path))
     writing_directory(1:lgt)=path(1:lgt)
  end if
  ! Add trailing slash if missing.
  lgt = len_trim(writing_directory)
  if (writing_directory(lgt:lgt) /= "/") &
       & writing_directory(min(lgt+1, len(writing_directory)):min(lgt+1, len(writing_directory))) = "/"

  ! Test if logging on disk is required.
  log_to_disk = (bigdft_mpi%ngroup > 1)
  call bigdft_get_run_properties(dict, log_to_disk = log_to_disk) !< May overwrite with user choice

  ! Save modified infos in dict.
  call bigdft_set_run_properties(dict, outdir_id = writing_directory, log_to_disk = log_to_disk)

  ! Now, create the logfile if needed.
  if (bigdft_mpi%iproc == 0) then
     if (log_to_disk) then
        ! Get Create log file name.
        call bigdft_get_run_properties(dict, naming_id = run_name)
        logfilename = "log" // trim(run_name) // ".yaml"
        path = trim(writing_directory)//trim(logfilename)
        call yaml_map('<BigDFT> log of the run will be written in logfile',path,unit=6)
        ! Check if logfile is already connected.
        call yaml_stream_connected(trim(path), unit_log, ierr)
        if (ierr /= 0) then
           ! Move possible existing log file.
           call ensure_log_file(trim(writing_directory), trim(logfilename), ierr)
           if (ierr /= 0) call MPI_ABORT(bigdft_mpi%mpi_comm,ierror,ierr)
           ! Close active stream and logfile if any. (TO BE MOVED IN RUN_UPDATE TO AVOID CLOSURE OF UPLEVEL INSTANCE)
           call yaml_get_default_stream(unit_log)
           if (unit_log /= 6) call yaml_close_stream(unit_log, ierr)
           !Create stream and logfile
           call yaml_set_stream(filename=trim(path),record_length=92,istat=ierr)
           !create that only if the stream is not already present, otherwise print a warning
           if (ierr == 0) then
              call yaml_get_default_stream(unit_log)
              call input_set_stdout(unit=unit_log)
              !call memocc_set_stdout(unit=70)
           else
              call yaml_warning('Logfile '//trim(path)//' cannot be created, stream already present. Ignoring...')
           end if
        else
           call yaml_release_document(unit_log)
           call yaml_set_default_stream(unit_log, ierr)
        end if ! Logfile already connected
     else
        !use stdout, do not crash if unit is present
        call yaml_set_stream(record_length=92,istat=ierr)
     end if ! Need to create a named logfile.

     !start writing on logfile
     call yaml_new_document()
     !welcome screen
     call print_logo()
  end if ! Logfile is created by master proc only
    
  if (bigdft_mpi%nproc >1) call processor_id_per_node(bigdft_mpi%iproc,bigdft_mpi%nproc,iproc_node,nproc_node)

  if (bigdft_mpi%iproc==0) then
     if (bigdft_mpi%nproc >1) call yaml_map('MPI tasks of root process node',nproc_node)
     call print_configure_options()
  end if
END SUBROUTINE create_log_file
