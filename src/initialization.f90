!> @file
!!  Routines to initialize some important structures for run_object
!! @author
!!    Copyright (C) 2007-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> Do all initialisation for all different files of BigDFT. 
!! Set default values if not any.
!! Initialize memocc
!! @todo
!!   Should be better for debug purpose to read input.perf before
subroutine bigdft_set_input(radical,posinp,in,atoms)
  use module_base
  use module_types
  use module_interfaces, except_this_one => bigdft_set_input
  use module_input_keys
  use yaml_output
  use dictionaries, only: dictionary
  implicit none

  !Arguments
  character(len=*),intent(in) :: posinp
  character(len=*),intent(in) :: radical
  type(input_variables), intent(inout) :: in
  type(atoms_data), intent(out) :: atoms

  character(len=*), parameter :: subname='bigdft_set_input'
  type(dictionary), pointer :: dict
!!$  logical :: exist_list
  call f_routine(id=subname)
  
  nullify(dict)
  atoms=atoms_null()
  ! Read atomic file
  call read_atomic_file(trim(posinp),bigdft_mpi%iproc,atoms%astruct)

  !read the input file(s) and transform them in the 
  call read_input_dict_from_files(trim(radical), bigdft_mpi,dict)

  call standard_inputfile_names(in,trim(radical))
  call inputs_from_dict(in, atoms, dict, .true.)

  call dict_free(dict)

  ! Generate the description of input variables.
  !if (bigdft_mpi%iproc == 0) then
  !   call input_keys_dump_def(trim(in%writing_directory) // "/input_help.yaml")
  !end if

  ! Read associated pseudo files.
  call init_atomic_values((bigdft_mpi%iproc == 0), atoms, in%ixc)
  call read_atomic_variables(atoms, trim(in%file_igpop),in%nspin)

  ! Start the signaling loop in a thread if necessary.
  if (in%signaling .and. bigdft_mpi%iproc == 0) then
     call bigdft_signals_init(in%gmainloop, 2, in%domain, len(trim(in%domain)))
     call bigdft_signals_start(in%gmainloop, in%signalTimeout)
  end if

  if (bigdft_mpi%iproc == 0) then
     call print_general_parameters(in,atoms)
     !call write_input_parameters(inputs,atoms)
  end if

  !if other steps are supposed to be done leave the last_run to minus one
  !otherwise put it to one
  if (in%last_run == -1 .and. in%ncount_cluster_x <=1 .or. in%ncount_cluster_x <= 1) then
     in%last_run = 1
  end if

  call f_release_routine()
END SUBROUTINE bigdft_set_input

!> De-allocate the variable of type input_variables
subroutine bigdft_free_input(in)
  use module_base
  use module_types
  use yaml_output
  type(input_variables), intent(inout) :: in
  
  call free_input_variables(in)
  call f_lib_finalize()
  !free all yaml_streams active
  call yaml_close_all_streams()

end subroutine bigdft_free_input

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
     !print *,'test',ncommands,icommands,command
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
       if (len_trim(run_id) > 0) then
          ierr=bigdft_error_ret(BIGDFT_INVALID,'run_id specified twice')
       end if
       ipos=index(command,'=')
       read(command(ipos+1:len(command)),*)run_id
    else if (index(command,'--runs-file=') > 0) then
       if (len_trim(posinp_file) > 0 .or. len_trim(run_id) >0) then
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

end subroutine command_line_information

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

  else if (matacc%iacceleration >= 2) then
     ! OpenCL convolutions are activated
     ! use CUBLAS for the linear algebra for the moment
     if (.not. OCLconv) then
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
        GPU%ndevices=min(GPU%ndevices,nproc_node)
        if (iproc == 0) then
           call yaml_map('Material acceleration','OpenCL',advance='no')
           call yaml_comment('iproc=0')
           call yaml_open_map('Number of OpenCL devices per node',flow=.true.)
           call yaml_map('used',trim(yaml_toa(min(GPU%ndevices,nproc_node),fmt='(i0)')))
           call yaml_map('available',trim(yaml_toa(GPU%ndevices,fmt='(i0)')))
           !write(*,'(1x,a,i5,i5)') 'OpenCL support activated, No. devices per node (used, available):',&
           !     min(GPU%ndevices,nproc_node),GPU%ndevices
           call yaml_close_map()
        end if
        !the number of devices is the min between the number of processes per node
        GPU%ndevices=min(GPU%ndevices,nproc_node)
        OCLconv=.true.
     end if

  else
     if (iproc == 0) then
        call yaml_map('Material acceleration',.false.,advance='no')
        call yaml_comment('iproc=0')
        ! write(*,'(1x,a)') 'No material acceleration (iproc=0)'
     end if
  end if

END SUBROUTINE init_material_acceleration


subroutine release_material_acceleration(GPU)
  use module_base
  use module_types
  implicit none
  type(GPU_pointers), intent(out) :: GPU
  
  if (GPUconv) then
     call sg_end()
  end if

  if (OCLconv) then
     call release_acceleration_OCL(GPU)
     OCLconv=.false.
  end if

END SUBROUTINE release_material_acceleration


!> Give the number of MPI processes per node (nproc_node) and before iproc (iproc_node)
subroutine processor_id_per_node(iproc,nproc,iproc_node,nproc_node)
  use module_base
  use module_types
  use dynamic_memory
  implicit none
  integer, intent(in) :: iproc,nproc
  integer, intent(out) :: iproc_node,nproc_node
  !local variables
  character(len=*), parameter :: subname='processor_id_per_node'
  integer :: ierr,namelen,i_stat,i_all,jproc
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
