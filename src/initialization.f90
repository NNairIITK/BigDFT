!> @file
!!  Routines to initialize some important structures for run_object
!! @author
!!    Copyright (C) 2007-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> Routines to handle the argument objects of call_bigdft().
subroutine run_objects_nullify(runObj)
  use module_types
  implicit none
  type(run_objects), intent(out) :: runObj

  nullify(runObj%user_inputs)
  nullify(runObj%inputs)
  nullify(runObj%atoms)
  nullify(runObj%rst)
  nullify(runObj%radii_cf)
END SUBROUTINE run_objects_nullify

subroutine run_objects_free(runObj, subname)
  use module_types
  use module_base
  use yaml_output
  use dictionaries
  implicit none
  type(run_objects), intent(inout) :: runObj
  character(len = *), intent(in) :: subname

  integer :: i_all, i_stat

  if (associated(runObj%user_inputs)) then
     call dict_free(runObj%user_inputs)
  end if
  if (associated(runObj%rst)) then
     call free_restart_objects(runObj%rst,subname)
     deallocate(runObj%rst)
  end if
  if (associated(runObj%atoms)) then
     call deallocate_atoms(runObj%atoms,subname) 
     deallocate(runObj%atoms)
  end if
  if (associated(runObj%inputs)) then
     call free_input_variables(runObj%inputs)
     deallocate(runObj%inputs)
  end if
  if (associated(runObj%radii_cf)) then
     i_all=-product(shape(runObj%radii_cf))*kind(runObj%radii_cf)
     deallocate(runObj%radii_cf,stat=i_stat)
     call memocc(i_stat,i_all,'radii_cf',subname)
  end if
  ! to be inserted again soon call f_lib_finalize()
  !call yaml_close_all_streams()
END SUBROUTINE run_objects_free

subroutine run_objects_free_container(runObj)
  use module_types
  use module_base
  use yaml_output
  implicit none
  type(run_objects), intent(inout) :: runObj

  integer :: i_all, i_stat

  ! User inputs are always owned by run objects.
  if (associated(runObj%user_inputs)) then
     call dict_free(runObj%user_inputs)
  end if
  ! Radii_cf are always owned by run objects.
  if (associated(runObj%radii_cf)) then
     i_all=-product(shape(runObj%radii_cf))*kind(runObj%radii_cf)
     deallocate(runObj%radii_cf,stat=i_stat)
     call memocc(i_stat,i_all,'radii_cf',"run_objects_free_container")
  end if
  ! Currently do nothing except nullifying everything.
  call run_objects_nullify(runObj)
END SUBROUTINE run_objects_free_container

subroutine run_objects_init_from_files(runObj, radical, posinp)
  use module_types
  use module_interfaces, only: atoms_new, rst_new, inputs_new, inputs_from_dict
  use module_input_dicts, only: user_dict_from_files
  implicit none
  type(run_objects), intent(out) :: runObj
  character(len = *), intent(in) :: radical, posinp

  integer :: i_stat
  integer(kind = 8) :: dummy

  call run_objects_nullify(runObj)

  call user_dict_from_files(runObj%user_inputs, radical, posinp, bigdft_mpi)

  call atoms_new(runObj%atoms)
  call inputs_new(runObj%inputs)
  call inputs_from_dict(runObj%inputs, runObj%atoms, runObj%user_inputs, .true.)

  ! Generate the description of input variables.
  !if (bigdft_mpi%iproc == 0) then
  !   call input_keys_dump_def(trim(in%writing_directory) // "/input_help.yaml")
  !end if

  if (bigdft_mpi%iproc == 0) then
     call print_general_parameters(runObj%inputs,runObj%atoms)
  end if

  allocate(runObj%radii_cf(runObj%atoms%astruct%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,runObj%radii_cf,'radii_cf',"run_objects_init_from_files")
  call read_radii_variables(runObj%atoms, runObj%radii_cf, &
       & runObj%inputs%crmult, runObj%inputs%frmult, runObj%inputs%projrad)

  call rst_new(dummy, runObj%rst)
  call restart_objects_new(runObj%rst)
  call restart_objects_set_mode(runObj%rst, runObj%inputs%inputpsiid)
  call restart_objects_set_nat(runObj%rst, runObj%atoms%astruct%nat, "run_objects_init_from_files")
  call restart_objects_set_mat_acc(runObj%rst, bigdft_mpi%iproc, runObj%inputs%matacc)

  ! Start the signaling loop in a thread if necessary.
  if (runObj%inputs%signaling .and. bigdft_mpi%iproc == 0) then
     call bigdft_signals_init(runObj%inputs%gmainloop, 2, &
          & runObj%inputs%domain, len_trim(runObj%inputs%domain))
     call bigdft_signals_start(runObj%inputs%gmainloop, runObj%inputs%signalTimeout)
  end if
END SUBROUTINE run_objects_init_from_files

!!$subroutine run_objects_nullify_from_objects(runObj, atoms, inputs, rst)
!!$  use module_types
!!$  implicit none
!!$  type(run_objects), intent(out) :: runObj
!!$  type(atoms_data), intent(in), target :: atoms
!!$  type(input_variables), intent(in), target, optional :: inputs
!!$  type(restart_objects), intent(in), target, optional :: rst
!!$
!!$  call run_objects_nullify(runObj)
!!$  runObj%atoms  => atoms
!!$  if (present(inputs)) then
!!$     runObj%inputs => inputs
!!$  else
!!$  end if
!!$  if (present(rst)) then
!!$     runObj%rst    => rst
!!$  else
!!$  end if
!!$END SUBROUTINE run_objects_nullify_from_objects

subroutine run_objects_associate(runObj, inputs, atoms, rst, rxyz0)
  use module_types
  implicit none
  type(run_objects), intent(out) :: runObj
  type(input_variables), intent(in), target :: inputs
  type(atoms_data), intent(in), target :: atoms
  type(restart_objects), intent(in), target :: rst
  real(gp), intent(in), optional :: rxyz0

  integer :: i_stat

  call run_objects_free_container(runObj)
  runObj%atoms  => atoms
  runObj%inputs => inputs
  runObj%rst    => rst
  if (present(rxyz0)) then
     call vcopy(3 * atoms%astruct%nat, rxyz0, 1, runObj%atoms%astruct%rxyz(1,1), 1)
  end if

  allocate(runObj%radii_cf(runObj%atoms%astruct%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,runObj%radii_cf,'radii_cf',"run_objects_init_from_files")
  call read_radii_variables(runObj%atoms, runObj%radii_cf, &
       & runObj%inputs%crmult, runObj%inputs%frmult, runObj%inputs%projrad)
END SUBROUTINE run_objects_associate

subroutine run_objects_system_setup(runObj, iproc, nproc, rxyz, shift, mem)
  use module_types
  use module_fragments
  use module_interfaces, only: system_initialization
  implicit none
  type(run_objects), intent(inout) :: runObj
  integer, intent(in) :: iproc, nproc
  real(gp), dimension(3,runObj%atoms%astruct%nat), intent(out) :: rxyz
  real(gp), dimension(3), intent(out) :: shift
  type(memory_estimation), intent(out) :: mem

  integer :: inputpsi, input_wf_format
  type(nonlocal_psp_descriptors) :: nlpspd
  real(wp), dimension(:), pointer :: proj
  type(system_fragment), dimension(:), pointer :: ref_frags
  character(len = *), parameter :: subname = "run_objects_estimate_memory"

  ! Copy rxyz since system_size() will shift them.
!!$  allocate(rxyz(3,runObj%atoms%astruct%nat+ndebug),stat=i_stat)
!!$  call memocc(i_stat,rxyz,'rxyz',subname)
  call vcopy(3 * runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1, rxyz(1,1), 1)

  call system_initialization(iproc, nproc, .false., inputpsi, input_wf_format, .true., &
       & runObj%inputs, runObj%atoms, rxyz, runObj%rst%KSwfn%orbs, &
       & runObj%rst%tmb%npsidim_orbs, runObj%rst%tmb%npsidim_comp, &
       & runObj%rst%tmb%orbs, runObj%rst%KSwfn%Lzd, runObj%rst%tmb%Lzd, &
       & nlpspd, runObj%rst%KSwfn%comms, shift, proj, runObj%radii_cf, &
       & ref_frags)
  call MemoryEstimator(nproc,runObj%inputs%idsx,runObj%rst%KSwfn%Lzd%Glr,&
       & runObj%rst%KSwfn%orbs%norb,runObj%rst%KSwfn%orbs%nspinor,&
       & runObj%rst%KSwfn%orbs%nkpts,nlpspd%nprojel,&
       & runObj%inputs%nspin,runObj%inputs%itrpmax,runObj%inputs%iscf,mem)

  ! De-allocations
!!$  i_all=-product(shape(rxyz))*kind(rxyz)
!!$  deallocate(rxyz,stat=i_stat)
!!$  call memocc(i_stat,i_all,'rxyz',subname)
  call deallocate_Lzd_except_Glr(runObj%rst%KSwfn%Lzd, subname)
  call deallocate_comms(runObj%rst%KSwfn%comms,subname)
  call deallocate_orbs(runObj%rst%KSwfn%orbs,subname)
  call deallocate_proj_descr(nlpspd,subname)  
END SUBROUTINE run_objects_system_setup


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
