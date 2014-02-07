!> @file
!! Manage dynamic memory allocation
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module used to manage memory allocations and de-allocations
module dynamic_memory

  !use m_profiling, except => ndebug, and=> d_nan, also=> r_nan
  use memory_profiling, except => ndebug
  use dictionaries, info_length => max_field_length
  use yaml_strings, only: yaml_toa,yaml_date_and_time_toa
  implicit none

  private 

  integer, parameter :: namelen=32          !< length of the character variables
  integer, parameter :: error_string_len=80 !< length of error string
  integer, parameter :: ndebug=0            !< size of debug parameters
  integer, parameter :: max_rank=7          !< maximum rank in fortran
  !maximum size of f_lib control variables
  integer, parameter :: max_ctrl = 5 !<maximum number of nested levels
  integer :: ictrl=0                 !<id of active control structure (<=max_ctrl)


  !> parameters for defitions of internal dictionary
  character(len=*), parameter :: arrayid='Array Id'
  character(len=*), parameter :: routineid='Allocating Routine Id'
  character(len=*), parameter :: sizeid='Size (Bytes)'
  character(len=*), parameter :: metadatadd='Address of metadata'
  character(len=*), parameter :: firstadd='Address of first element'
  character(len=*), parameter :: processid='Process Id'
  character(len=*), parameter :: subprograms='Subroutines'
  character(len=*), parameter :: no_of_calls='No. of calls'
  character(len=*), parameter :: t0_time='Time of last opening'
  character(len=*), parameter :: tot_time='Total time (s)'
  character(len=*), parameter :: prof_enabled='Profiling Enabled'

  !error codes
  integer :: ERR_ALLOCATE
  integer :: ERR_DEALLOCATE
  integer :: ERR_MEMLIMIT
  integer :: ERR_INVALID_MALLOC
  integer :: ERR_INVALID_RANK
  integer :: ERR_MALLOC_INTERNAL

  !> control structure of flib library. 
  !Contains all global variables of interest in a separate instance of f_lib
  type :: mem_ctrl 
     logical :: profile_initialized  !< global variables for initialization
     logical :: routine_opened       !< global variable (can be stored in dictionaries)
     logical :: profile_routine      !< decide whether the routine has to be profiled
     character(len=namelen) :: present_routine !< name of the active routine 
     !>dictionaries needed for profiling storage
     type(dictionary), pointer :: dict_global    !<status of the memory at higher level
     type(dictionary), pointer :: dict_routine   !<status of the memory inside the routine
     type(dictionary), pointer :: dict_calling_sequence !<profiling of the routines
     type(dictionary), pointer :: dict_codepoint !<points to where we are in the previous dictionary
  end type mem_ctrl
  
  !>global variable controlling the different instances of the calls
  !the 0 component is supposed to be unused, it is allocated to avoid segfaults
  ! if the library routines are called without initialization
  type(mem_ctrl), dimension(0:max_ctrl) :: mems

  !> Structure needed to allocate an allocatable array
  type, public :: malloc_information_all
     logical :: pin                          !< flag to control the pinning of the address
     logical :: profile                      !< activate profiling for this allocation
     logical :: put_to_zero                  !< initialize to zero after allocation
     integer :: rank                         !< rank of the array
     integer, dimension(max_rank) :: shape   !< shape of the structure 
     integer, dimension(max_rank) :: lbounds !< lower bounds
     integer, dimension(max_rank) :: ubounds !< upper bounds
     integer(kind=8) :: srcdata_add          !< physical address of source data
     character(len=namelen) :: array_id      !< label the array
     character(len=namelen) :: routine_id    !< label the routine
     
  end type malloc_information_all

  !> Structure needed to allocate an allocatable array of string of implicit length (for non-2003 compilers)
  type, public :: malloc_information_str_all
     logical :: pin                          !< flag to control the pinning of the address
     logical :: profile                      !< activate profiling for this allocation
     logical :: put_to_zero                  !< initialize to zero after allocation
     integer :: rank                         !< rank of the array
     integer :: len                          !< length of the character
     integer, dimension(max_rank) :: shape   !< shape of the structure 
     integer, dimension(max_rank) :: lbounds !< lower bounds
     integer, dimension(max_rank) :: ubounds !< upper bounds
     integer(kind=8) :: srcdata_add          !< physical address of source data
     character(len=namelen) :: array_id      !< label the array
     character(len=namelen) :: routine_id    !< label the routine
  end type malloc_information_str_all

  !> Structure needed to allocate a pointer
  type, public :: malloc_information_ptr
     logical :: ptr                          !< just to make the structures different, to see if needed
     logical :: pin                          !< flag to control the pinning of the address
     logical :: profile                      !< activate profiling for this allocation
     logical :: put_to_zero                  !< initialize to zero after allocation
     integer :: rank                         !< rank of the pointer
     integer, dimension(max_rank) :: shape   !< shape of the structure 
     integer, dimension(max_rank) :: lbounds !< lower bounds
     integer, dimension(max_rank) :: ubounds !< upper bounds
     integer(kind=8) :: srcdata_add          !< physical address of source data
     character(len=namelen) :: array_id      !< label the array
     character(len=namelen) :: routine_id    !< label the routine
  end type malloc_information_ptr

  !> Structure needed to allocate a pointer of string of implicit length (for non-2003 complilers)
  type, public :: malloc_information_str_ptr
     logical :: ptr                          !< just to make the structures different, to see if needed
     logical :: pin                          !< flag to control the pinning of the address
     logical :: profile                      !< activate profiling for this allocation
     logical :: put_to_zero                  !< initialize to zero after allocation
     integer :: rank                         !< rank of the pointer
     integer :: len                          !< length of the character
     integer, dimension(max_rank) :: shape   !< shape of the structure 
     integer, dimension(max_rank) :: lbounds !< lower bounds
     integer, dimension(max_rank) :: ubounds !< upper bounds
     integer(kind=8) :: srcdata_add          !< physical address of source data
     character(len=namelen) :: array_id      !< label the array
     character(len=namelen) :: routine_id    !< label the routine
  end type malloc_information_str_ptr

  type, public :: array_bounds
     integer :: nlow  !<lower bounds
     integer :: nhigh !<higher bounds
  end type array_bounds

  interface assignment(=)
     module procedure i1_all,i2_all,i3_all,i4_all
     module procedure l1_all,l2_all,l3_all
     module procedure d1_all,d2_all,d3_all,d4_all,d5_all,d6_all
     module procedure r1_all,r2_all,r3_all
     module procedure d1_ptr,d2_ptr,d3_ptr,d4_ptr,d5_ptr
     module procedure i1_ptr,i2_ptr
     !strings and pointers for characters
     module procedure c1_all
!     module procedure c1_ptr
  end interface

  interface operator(.to.)
     module procedure bounds
  end interface

  interface nullify_malloc_information
     module procedure nullify_malloc_information_all
     module procedure nullify_malloc_information_ptr
     module procedure nullify_malloc_information_str_all
     module procedure nullify_malloc_information_str_ptr
  end interface nullify_malloc_information

  interface f_free
     module procedure i1_all_free,i2_all_free,i3_all_free,i4_all_free
     module procedure i1_all_free_multi
     module procedure l1_all_free,l2_all_free,l3_all_free
     module procedure d1_all_free,d2_all_free,d1_all_free_multi,d3_all_free,d4_all_free,d5_all_free,d6_all_free
     module procedure r1_all_free,r2_all_free,r3_all_free
  end interface

  interface f_free_ptr
     module procedure i1_ptr_free,i2_ptr_free
     module procedure i1_ptr_free_multi
     module procedure d1_ptr_free,d2_ptr_free,d3_ptr_free,d4_ptr_free,d5_ptr_free
  end interface

!  interface f_free_str
!     module procedure c1_all_free
!  end interface f_free_str

!  interface f_free_str_ptr
!     module procedure c1_ptr_free
!  end interface f_free_str_ptr


  interface f_malloc
     module procedure f_malloc,f_malloc_simple
     module procedure f_malloc_bounds,f_malloc_bound
     !here also the procedures for the copying of arrays have to be defined
     module procedure f_malloc_i2
  end interface

  interface f_malloc0
     module procedure f_malloc0,f_malloc0_simple
     module procedure f_malloc0_bounds,f_malloc0_bound
  end interface

  interface f_malloc_ptr
     module procedure f_malloc_ptr,f_malloc_ptr_simple
     module procedure f_malloc_ptr_bounds,f_malloc_ptr_bound
     module procedure f_malloc_ptr_i2
  end interface

  interface f_malloc0_ptr
     module procedure f_malloc0_ptr,f_malloc0_ptr_simple
     module procedure f_malloc0_ptr_bounds,f_malloc0_ptr_bound
  end interface

  interface f_malloc_str
     module procedure f_malloc_str,f_malloc_str_simple
     module procedure f_malloc_str_bounds,f_malloc_str_bound
  end interface

  interface f_malloc0_str
     module procedure f_malloc0_str,f_malloc0_str_simple
     module procedure f_malloc0_str_bounds,f_malloc0_str_bound
  end interface

  interface f_malloc_str_ptr
     module procedure f_malloc_str_ptr,f_malloc_str_ptr_simple
     module procedure f_malloc_str_ptr_bounds,f_malloc_str_ptr_bound
  end interface

  interface f_malloc0_str_ptr
     module procedure f_malloc0_str_ptr,f_malloc0_str_ptr_simple
     module procedure f_malloc0_str_ptr_bounds,f_malloc0_str_ptr_bound
  end interface


  !to be verified if clock_gettime is without side-effect, otherwise the routine cannot be pure
  interface
     pure subroutine nanosec(itime)
       implicit none
       integer(kind=8), intent(out) :: itime
     end subroutine nanosec
  end interface


  !> Public routines
  public :: f_malloc,f_malloc0,f_malloc_ptr,f_malloc0_ptr,f_malloc_dump_status
  public :: f_malloc_str,f_malloc0_str,f_malloc_str_ptr,f_malloc0_str_ptr
  public :: f_free,f_free_ptr,f_free_str,f_free_str_ptr
  public :: f_routine,f_release_routine,f_malloc_set_status,f_malloc_initialize,f_malloc_finalize
  public :: f_time
  public :: assignment(=),operator(.to.)

  !for internal f_lib usage
  public :: dynamic_memory_errors

contains

  pure function f_time()
    integer(kind=8) :: f_time
    !local variables
    integer(kind=8) :: itime
    call nanosec(itime)
    f_time=itime
  end function f_time

  elemental pure function bounds(nlow,nhigh)
    implicit none
    integer, intent(in) :: nlow,nhigh
    type(array_bounds) :: bounds

    bounds%nlow=nlow
    bounds%nhigh=nhigh
  end function bounds

  pure function mem_ctrl_null() result(mem)
    type(mem_ctrl) :: mem
    call nullify_mem_ctrl(mem)
  end function mem_ctrl_null
  pure subroutine nullify_mem_ctrl(mem)
    implicit none
    type(mem_ctrl), intent(out) :: mem
    mem%profile_initialized=.false. 
    mem%routine_opened=.false.      
    mem%profile_routine=.true.
    mem%present_routine=repeat(' ',namelen)
    !>dictionaries needed for profiling storage
    nullify(mem%dict_global)
    nullify(mem%dict_routine)
    nullify(mem%dict_calling_sequence)
    nullify(mem%dict_codepoint)
  end subroutine nullify_mem_ctrl

  !pure 
  function mem_ctrl_init() result(mem)
    type(mem_ctrl) :: mem
    call nullify_mem_ctrl(mem)
    call initialize_mem_ctrl(mem)
  end function mem_ctrl_init
  !pure 
  subroutine initialize_mem_ctrl(mem)
    implicit none
    type(mem_ctrl), intent(out) :: mem
    mem%profile_initialized=.true.
    !initalize the dictionary with the allocation information
    nullify(mem%dict_routine)
    call dict_init(mem%dict_global)
    call dict_init(mem%dict_calling_sequence)
    !in principle the calling sequence starts from the main
    mem%dict_codepoint => mem%dict_calling_sequence
  end subroutine initialize_mem_ctrl

  pure subroutine nullify_malloc_information_all(m)
    implicit none
    type(malloc_information_all), intent(out) :: m
    include 'f_malloc-null-inc.f90'
  end subroutine nullify_malloc_information_all

  pure subroutine nullify_malloc_information_ptr(m)
    implicit none
    type(malloc_information_ptr), intent(out) :: m
    include 'f_malloc-null-inc.f90'
    m%ptr=.true.
  end subroutine nullify_malloc_information_ptr

  pure subroutine nullify_malloc_information_str_all(m)
    implicit none
    type(malloc_information_str_all), intent(out) :: m
    include 'f_malloc-null-inc.f90'
    m%len=0
  end subroutine nullify_malloc_information_str_all
  
  pure subroutine nullify_malloc_information_str_ptr(m)
    implicit none
    type(malloc_information_str_ptr), intent(out) :: m
    include 'f_malloc-null-inc.f90'
    m%len=0
    m%ptr=.true.
  end subroutine nullify_malloc_information_str_ptr


  !> This routine adds the corresponding subprogram name to the dictionary
  !! and prepend the dictionary to the global info dictionary
  !! if it is called more than once for the same name it has no effect
  subroutine f_routine(id,profile)
    !use yaml_output !debug
    implicit none
    logical, intent(in), optional :: profile
    character(len=*), intent(in), optional :: id
    
    !local variables
    integer :: lgt,ncalls
    integer(kind=8) :: itime

    if (f_err_raise(ictrl == 0,&
         'ERROR (f_routine): the routine f_malloc_initialize has not been called',&
         ERR_MALLOC_INTERNAL)) return

    if (.not. present(id)) return !no effect

    !take the time
    itime=f_time()

    !desactivate profile_routine if the mother routine has desactivated it
    if (present(profile)) mems(ictrl)%profile_routine=mems(ictrl)%profile_routine .and. profile

    !print *,'test',trim(mems(ictrl)%present_routine),'test2',trim(id),'ccc',mems(ictrl)%routine_opened
    !if (trim(mems(ictrl)%present_routine) /= trim(id) .or. &
    !         (trim(mems(ictrl)%present_routine) == trim(id) .and. .not. mems(ictrl)%routine_opened) ) then
    !debug
!!$    call yaml_open_map('Status before the opening of the routine')
!!$      call yaml_map('Level',ictrl)
!!$      call yaml_map('Routine opened',mems(ictrl)%routine_opened)
!!$      call yaml_map('Codepoint',trim(dict_key(mems(ictrl)%dict_codepoint)))
!!$      call yaml_open_map('Codepoint dictionary')
!!$      call yaml_dict_dump(mems(ictrl)%dict_codepoint)
!!$      call yaml_close_map()
!!$    call yaml_close_map()
    !end debug  

    if (.true.) then
       if(associated(mems(ictrl)%dict_routine)) then
          call prepend(mems(ictrl)%dict_global,mems(ictrl)%dict_routine)
          nullify(mems(ictrl)%dict_routine)
       end if
       !this means that the previous routine has not been closed
       if (mems(ictrl)%routine_opened) then
          !call open_routine(dict_codepoint)
          mems(ictrl)%dict_codepoint=>mems(ictrl)%dict_codepoint//subprograms
       end if
       mems(ictrl)%routine_opened=.true.
       !call add(dict_codepoint,trim(id))
       !see if the key existed in the codepoint
       if (has_key(mems(ictrl)%dict_codepoint,trim(id))) then
          !retrieve number of calls and increase it
          ncalls=mems(ictrl)%dict_codepoint//trim(id)//no_of_calls
          call set(mems(ictrl)%dict_codepoint//trim(id)//no_of_calls,ncalls+1)
          !write the starting point for the time
          call set(mems(ictrl)%dict_codepoint//trim(id)//t0_time,itime)
          call set(mems(ictrl)%dict_codepoint//trim(id)//prof_enabled,mems(ictrl)%profile_routine)
       else
          !create a new dictionary
          call set(mems(ictrl)%dict_codepoint//trim(id),&
               dict_new((/no_of_calls .is. yaml_toa(1), t0_time .is. yaml_toa(itime),&
                             tot_time .is. yaml_toa(0.d0,fmt='(f4.1)'), &
                             prof_enabled .is. yaml_toa(mems(ictrl)%profile_routine)/)))
       end if
       !then fix the new codepoint from this one
       mems(ictrl)%dict_codepoint=>mems(ictrl)%dict_codepoint//trim(id)

       lgt=min(len(id),namelen)
       mems(ictrl)%present_routine=repeat(' ',namelen)
       mems(ictrl)%present_routine(1:lgt)=id(1:lgt)

    end if
  end subroutine f_routine

  !> Close a previously opened routine
  subroutine f_release_routine()
    use yaml_output, only: yaml_dict_dump
    implicit none

    if (f_err_raise(ictrl == 0,&
         'ERROR (f_release_routine): the routine f_malloc_initialize has not been called',&
         ERR_MALLOC_INTERNAL)) return

    if (associated(mems(ictrl)%dict_routine)) then
       call prepend(mems(ictrl)%dict_global,mems(ictrl)%dict_routine)
       nullify(mems(ictrl)%dict_routine)
    end if
    !call yaml_map('Closing routine',trim(dict_key(dict_codepoint)))

    call close_routine(mems(ictrl)%dict_codepoint,.not. mems(ictrl)%routine_opened)!trim(dict_key(dict_codepoint)))

    if (f_err_check()) then
!!$       call yaml_warning('ERROR found!')
!!$       call f_dump_last_error()
!!$       call yaml_comment('End of ERROR')
       return
    end if
    !last_opened_routine=trim(dict_key(dict_codepoint))!repeat(' ',namelen)
    !the main program is opened until there is a subprograms keyword
    if (f_err_raise(.not. associated(mems(ictrl)%dict_codepoint%parent),'parent not associated(A)',&
         ERR_MALLOC_INTERNAL)) return

    if (dict_key(mems(ictrl)%dict_codepoint%parent) == subprograms) then
    
       mems(ictrl)%dict_codepoint=>mems(ictrl)%dict_codepoint%parent

       if (f_err_raise(.not. associated(mems(ictrl)%dict_codepoint%parent),'parent not associated(B)',&
            ERR_MALLOC_INTERNAL)) return
       mems(ictrl)%dict_codepoint=>mems(ictrl)%dict_codepoint%parent
    else !back in the main program
       mems(ictrl)%routine_opened=.false.
    end if

    mems(ictrl)%present_routine(1:len(mems(ictrl)%present_routine))=trim(dict_key(mems(ictrl)%dict_codepoint))
    if (.not. has_key(mems(ictrl)%dict_codepoint,prof_enabled)) then
       call yaml_dict_dump(mems(ictrl)%dict_codepoint)
       call f_err_throw('The key '//prof_enabled//' is not present in the codepoint',&
            err_id=ERR_MALLOC_INTERNAL)
       return
    end if

    mems(ictrl)%profile_routine=mems(ictrl)%dict_codepoint//prof_enabled! 

    !debug
!!$    call yaml_open_map('Codepoint after closing')
!!$    call yaml_map('Potential Reference Routine',trim(dict_key(mems(ictrl)%dict_codepoint)))
!!$    call yaml_dict_dump(mems(ictrl)%dict_codepoint)
!!$    call yaml_close_map()
!!$    call yaml_comment('End of release routine',hfill='=')
    !end debug
  end subroutine f_release_routine

  !>create the id of a new routine in the codepoint and points to it.
  !! works for sequences
  subroutine open_routine(dict)
    implicit none
    type(dictionary), pointer :: dict
    !local variables
    integer :: ival
    character(len=info_length) :: routinename
    type(dictionary), pointer :: dict_tmp

    !now imagine that a new routine is created
    ival=dict_len(dict)-1
    routinename=dict//ival

    !call yaml_map('The routine which has to be converted is',trim(routinename))

    call pop(dict,ival)

    dict_tmp=>dict//ival//trim(routinename)

    dict => dict_tmp
    nullify(dict_tmp)

  end subroutine open_routine

  subroutine close_routine(dict,jump_up)
    !use yaml_output !debug
    implicit none
    type(dictionary), pointer :: dict
    logical, intent(in) :: jump_up
    !character(len=*), intent(in) :: name
    !local variables
    integer(kind=8) :: itime,jtime
    real(kind=8) :: rtime
    type(dictionary), pointer :: dict_tmp

    if (f_err_raise(.not. associated(dict),'routine not associated',ERR_MALLOC_INTERNAL)) return

    itime=f_time()
    !debug
!!$    call yaml_open_map('codepoint'//trim(dict_key(dict)))
!!$    call yaml_dict_dump(dict)
!!$    call yaml_close_map()
!!$    call yaml_comment('We should jump up '//trim(yaml_toa(jump_up)),hfill='}')
    !end debug

    !update the total time, if the starting point is present
    if (has_key(dict,t0_time)) then
       jtime=dict//t0_time
       jtime=itime-jtime
       rtime=dict//tot_time
       call set(dict//tot_time,rtime+real(jtime,kind=8)*1.d-9,fmt='(1pe15.7)')
       call pop(dict,t0_time)
    else
       call f_err_throw('Key '//t0_time//&
            ' not found, most likely f_release_routine has been called too many times',&
            err_id=ERR_INVALID_MALLOC)
    end if

    !we should go up of three levels
    if (jump_up) then
       dict_tmp=>dict%parent
       if (f_err_raise(.not. associated(dict_tmp),'parent not associated(1)',&
         ERR_MALLOC_INTERNAL)) return
!       call yaml_map('Present Key 1',dict_key(dict_tmp))
       dict_tmp=>dict_tmp%parent
       if (f_err_raise(.not. associated(dict_tmp),'parent not associated(2)',&
            ERR_MALLOC_INTERNAL)) return
!       call yaml_map('Present Key 2',dict_key(dict_tmp))
       if (f_err_raise(.not. associated(dict_tmp%parent),'parent not associated(3)',&
            ERR_MALLOC_INTERNAL)) return
       dict_tmp=>dict_tmp%parent
       if (f_err_raise(.not. associated(dict_tmp%parent),'parent not associated(4)',&
            ERR_MALLOC_INTERNAL)) return
       dict=>dict_tmp%parent
    end if

  end subroutine close_routine

  !routine which is called for most of the errors of the module
  subroutine f_malloc_callback()
    use yaml_output, only: yaml_warning
    implicit none

    call yaml_warning('An error occured in dynamic memory module. Printing info')
    call f_malloc_dump_status()
    call f_err_severe()
  end subroutine f_malloc_callback

  !> Decide the error messages associated to the dynamic memory
  subroutine dynamic_memory_errors()
    use dictionaries, only: f_err_define
    implicit none
    
    call f_err_define(err_name='ERR_ALLOCATE',err_msg='Allocation error',err_id=ERR_ALLOCATE,&
         err_action='Control the order of the allocation of if the memory limit has been reached',&
         callback=f_malloc_callback)
    call f_err_define(err_name='ERR_DEALLOCATE',err_msg='Deallocation error',err_id=ERR_DEALLOCATE,&
         err_action='Control the order of the allocation of if the memory limit has been reached',&
         callback=f_malloc_callback)
    call f_err_define(err_name='ERR_MEMLIMIT',err_msg='Memory limit reached',err_id=ERR_MEMLIMIT,&
         err_action='Control the size of the arrays needed for this run with bigdft-tool program',&
         callback=f_malloc_callback)
    call f_err_define(err_name='ERR_INVALID_MALLOC',err_msg='Invalid specification of f_malloc',&
         err_id=ERR_INVALID_MALLOC,&
         err_action='Put coherent data for the memory space allocation',&
         callback=f_malloc_callback)
    call f_err_define(err_name='ERR_MALLOC_INTERNAL',err_msg='Internal error of memory profiler',&
         err_id=ERR_MALLOC_INTERNAL,&
         err_action='An invalid operation occurs, submit bug report to developers',&
         callback=f_malloc_callback)
    
  end subroutine dynamic_memory_errors

  !> opens a new instance of the dynamic memory handling
  subroutine f_malloc_initialize()
    implicit none
    
    !increase the number of active instances
    ictrl=ictrl+1
    if (f_err_raise(ictrl > max_ctrl,&
         'The number of active instances cannot exceed'//trim(yaml_toa(max_ctrl)),&
         ERR_MALLOC_INTERNAL)) return

    !extra options can be passed at the initialization
    mems(ictrl)=mem_ctrl_init()

    !initialize the memprofiling counters
    call set(mems(ictrl)%dict_global//'Timestamp of Profile initialization',&
         trim(yaml_date_and_time_toa()))
    !Process Id (used to dump)
    call set(mems(ictrl)%dict_global//processid,0)
    !start the profiling of the main program
    call f_routine(id='Main program')

    !set status of library to the initial case
    call f_malloc_set_status(memory_limit=0.e0)

  end subroutine f_malloc_initialize

  !> Initialize the library
  subroutine f_malloc_set_status(memory_limit,output_level,logfile_name,unit,iproc)
    use yaml_output, only: yaml_date_and_time_toa
    implicit none
    !Arguments
    character(len=*), intent(in), optional :: logfile_name   !< Name of the logfile
    real(kind=4), intent(in), optional :: memory_limit       !< Memory limit
    integer, intent(in), optional :: output_level            !< Level of output for memocc
                                                             !! 0 no file, 1 light, 2 full
    integer, intent(in), optional :: unit                    !< Indicate file unit for the output
    integer, intent(in), optional :: iproc                   !< Process Id (used to dump, by default one 0)

    if (f_err_raise(ictrl == 0,&
         'ERROR (f_malloc_set_status): the routine f_malloc_initialize has not been called',&
         ERR_MALLOC_INTERNAL)) return

!!$    if (.not. mems(ictrl)%profile_initialized) then
!!$       profile_initialized=.true.
!!$       !call malloc_errors()
!!$       !initalize the dictionary with the allocation information
!!$       nullify(dict_routine)
!!$       call dict_init(dict_global)
!!$       call set(dict_global//'Timestamp of Profile initialization',trim(yaml_date_and_time_toa()))
!!$       !Process Id (used to dump)
!!$       call set(dict_global//processid,0)
!!$       call dict_init(dict_calling_sequence)
!!$       !in principle the calling sequence starts from the main
!!$       dict_codepoint => dict_calling_sequence
!!$       call f_routine(id='Main program')
!!$    end if

    if (present(memory_limit)) call memocc_set_memory_limit(memory_limit)

    if (present(output_level)) call memocc_set_state(output_level)

    if (present(unit)) call memocc_set_stdout(unit)

    if (present(logfile_name)) call memocc_set_filename(logfile_name)
       
    if (present(iproc)) call set(mems(ictrl)%dict_global//processid,iproc)
  end subroutine f_malloc_set_status

  !> Finalize f_malloc (Display status)
  subroutine f_malloc_finalize(dump,process_id)
    use yaml_output, only: yaml_warning,yaml_open_map,yaml_close_map,yaml_dict_dump,yaml_get_default_stream
    implicit none
    !Arguments
    logical, intent(in), optional :: dump !< Dump always information, 
                                          !! otherwise only for Process Id == 0 and errors
    integer, intent(out), optional :: process_id !< retrieve the process_id  
                                                 !! useful to dump further information after finalization
    !local variables
    integer :: pid
    logical :: dump_status
    !integer :: unt
    
    if (f_err_raise(ictrl == 0,&
         'ERROR (f_malloc_finalize): the routine f_malloc_initialize has not been called',&
         ERR_MALLOC_INTERNAL)) return

    if (present(process_id)) process_id=-1

    !quick return if variables not associated
    if (associated(mems(ictrl)%dict_global)) then
       !put the last values in the dictionary if not freed
       if (associated(mems(ictrl)%dict_routine)) then
          call prepend(mems(ictrl)%dict_global,mems(ictrl)%dict_routine)
          nullify(mems(ictrl)%dict_routine)
       end if
       if (present(process_id)) process_id = mems(ictrl)%dict_global//processid

       if (present(dump)) then
          dump_status=.true.
       else 
          pid = mems(ictrl)%dict_global//processid
          if (pid == 0) then
             dump_status=.true.
          else
             dump_status=.false.
          end if
          !Print if error
          if (dict_size(mems(ictrl)%dict_global) == 2) dump_status=.false.
       end if
       if (dump_status) then
          call yaml_open_map('Status of the memory at finalization')
          !call yaml_dict_dump(dict_global)
          call dump_leaked_memory(mems(ictrl)%dict_global)
          call yaml_close_map()
       end if
       call dict_free(mems(ictrl)%dict_global)
       call f_release_routine() !release main
       !    call yaml_open_map('Calling sequence')
       !    call yaml_dict_dump(dict_calling_sequence)
       !    call yaml_close_map()
       call dict_free(mems(ictrl)%dict_calling_sequence)
    end if

    if (mems(ictrl)%profile_initialized) call memocc_report()

    !nullify control structure
    mems(ictrl)=mem_ctrl_null()
    !lower the level
    ictrl=ictrl-1
  end subroutine f_malloc_finalize

  !> Dump all allocations
  subroutine dump_leaked_memory(dict)
    use metadata_interfaces, only: address_toi
     use yaml_output
     implicit none
     type(dictionary), pointer, intent(in) :: dict
     !Local variables
     type(dictionary), pointer :: dict_ptr!, dict_tmp
     character(len=256) :: array_id
     dict_ptr => dict_next(dict)
     do while(associated(dict_ptr))
        if (has_key(dict_ptr,trim(arrayid))) then
           array_id = dict_ptr//arrayid
           call yaml_open_map(trim(array_id))
           call yaml_dict_dump(dict_ptr)
           call yaml_map(metadatadd,trim(dict_key(dict_ptr)))
           call yaml_close_map()
        else
           call yaml_map(trim(dict_key(dict_ptr)),dict_ptr)

!!$           call yaml_dict_dump(dict_ptr)
!!$           call yaml_close_map()
        end if
        dict_ptr=>dict_next(dict_ptr)
     end do
  end subroutine dump_leaked_memory

  subroutine f_malloc_dump_status()
    use yaml_output
    implicit none

    if (f_err_raise(ictrl == 0,&
         'ERROR (f_malloc_dump_status): the routine f_malloc_initialize has not been called',&
         ERR_MALLOC_INTERNAL)) return

    call yaml_newline()
    call yaml_open_map('Calling sequence of Main program')
      call yaml_dict_dump(mems(ictrl)%dict_calling_sequence)
    call yaml_close_map()
    if (associated(mems(ictrl)%dict_routine)) then
       call yaml_open_map('Routine dictionary')
       call dump_leaked_memory(mems(ictrl)%dict_routine)
       call yaml_close_map()
    end if
    call yaml_open_map('Global dictionary')
    call dump_leaked_memory(mems(ictrl)%dict_global)
    call yaml_close_map()

  end subroutine f_malloc_dump_status

!---routines for low-level dynamic memory handling

  !> For rank-1 arrays
  pure function f_malloc_simple(size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-simple-inc.f90'
  end function f_malloc_simple
  !> For rank-1 arrays
  pure function f_malloc0_simple(size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-simple-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_simple
  !> For rank-1 arrays
  pure function f_malloc_ptr_simple(size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-simple-inc.f90'
  end function f_malloc_ptr_simple
  !> For rank-1 arrays
  pure function f_malloc0_ptr_simple(size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-simple-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_ptr_simple
  !> For rank-1 arrays
  pure function f_malloc_str_simple(length,size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_malloc-simple-inc.f90'
    m%len=length
  end function f_malloc_str_simple
  !> For rank-1 arrays
  pure function f_malloc0_str_simple(length,size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_malloc-simple-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloc0_str_simple
  !> For rank-1 arrays
  pure function f_malloc_str_ptr_simple(length,size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloc-simple-inc.f90'
    m%len=length
  end function f_malloc_str_ptr_simple
  !> For rank-1 arrays
  pure function f_malloc0_str_ptr_simple(length,size,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloc-simple-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloc0_str_ptr_simple


  !> For rank-1 arrays, with bounds
  pure function f_malloc_bound(bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-bound-inc.f90'
  end function f_malloc_bound
  !>for rank-1 arrays, with boundaries
  pure function f_malloc0_bound(bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-bound-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_bound
  !> For rank-1 arrays
  pure function f_malloc_ptr_bound(bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-bound-inc.f90'
  end function f_malloc_ptr_bound
  !> For rank-1 arrays
  pure function f_malloc0_ptr_bound(bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-bound-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_ptr_bound
  !> For rank-1 arrays, with bounds
  pure function f_malloc_str_bound(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_malloc-bound-inc.f90'
    m%len=length
  end function f_malloc_str_bound
  !>for rank-1 arrays, with boundaries
  pure function f_malloc0_str_bound(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_malloc-bound-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloc0_str_bound
  !> For rank-1 arrays
  pure function f_malloc_str_ptr_bound(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloc-bound-inc.f90'
    m%len=length
  end function f_malloc_str_ptr_bound
  !> For rank-1 arrays
  pure function f_malloc0_str_ptr_bound(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloc-bound-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloc0_str_ptr_bound


  !> Define the allocation information for  arrays of different rank
  pure function f_malloc_bounds(bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-bounds-inc.f90'
  end function f_malloc_bounds
  pure function f_malloc0_bounds(bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-bounds-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_bounds
  !> Define the allocation information for  arrays of different rank
  pure function f_malloc_ptr_bounds(bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-bounds-inc.f90'
  end function f_malloc_ptr_bounds
  !> Define the allocation information for  arrays of different rank
  pure function f_malloc0_ptr_bounds(bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-bounds-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_ptr_bounds
  !> Define the allocation information for  arrays of different rank
  pure function f_malloc_str_bounds(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_malloc-bounds-inc.f90'
    m%len=length
  end function f_malloc_str_bounds
  pure function f_malloc0_str_bounds(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_malloc-bounds-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloc0_str_bounds
  !> Define the allocation information for  arrays of different rank
  pure function f_malloc_str_ptr_bounds(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloc-bounds-inc.f90'
    m%len=length
  end function f_malloc_str_ptr_bounds
  !> Define the allocation information for  arrays of different rank
  pure function f_malloc0_str_ptr_bounds(length,bounds,id,routine_id,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloc-bounds-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloc0_str_ptr_bounds

  !> Define the allocation information for  arrays of different rank
  function f_malloc(sizes,id,routine_id,lbounds,ubounds,profile,src) result(m)
    implicit none
    !the integer array src is here added to avoid problems in resolving the ambiguity with f_malloc_src
    integer, dimension(:), intent(in), optional :: src
    type(malloc_information_all) :: m
    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
    !local variables
    integer :: i
    include 'f_malloc-base-inc.f90'
    if (present(src)) then
       include 'f_malloc-inc.f90'
       !when src is given there is no need anymore to continue the routine
       if (present(lbounds) .or. present(ubounds) .or. present(sizes)) then
          call f_err_throw(&
               'The presence of lbounds, ubounds or sizes is forbidden whe src is present',&
               ERR_INVALID_MALLOC)
       end if
       return
    end if
    include 'f_malloc-extra-inc.f90'
  end function f_malloc
  !> define the allocation information for  arrays of different rank
  function f_malloc0(sizes,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-total-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0
  !> Define the allocation information for  arrays of different rank
  function f_malloc_ptr(sizes,id,routine_id,lbounds,ubounds,profile,src) result(m)
    implicit none
    !the integer array src is here added to avoid problems in resolving the ambiguity
    integer, dimension(:), intent(in), optional :: src
    type(malloc_information_ptr) :: m
    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
    !local variables
    integer :: i
    include 'f_malloc-base-inc.f90'
    if (present(src)) then
       include 'f_malloc-inc.f90'
       !when src is given there is no need anymore to continue the routine
       if (present(lbounds) .or. present(ubounds) .or. present(sizes)) then
          call f_err_throw(&
               'The presence of lbounds, ubounds or sizes is forbidden whe src is present',&
               ERR_INVALID_MALLOC)
       end if
       return
    end if
    include 'f_malloc-extra-inc.f90'
  end function f_malloc_ptr
  !> Define the allocation information for  arrays of different rank
  function f_malloc0_ptr(sizes,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-total-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_ptr
  !> Define the allocation information for  arrays of different rank
  function f_malloc_str(length,sizes,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_malloc-total-inc.f90'
    m%len=length
  end function f_malloc_str
  !> define the allocation information for  arrays of different rank
  function f_malloc0_str(length,sizes,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_malloc-total-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloc0_str
  !> Define the allocation information for  arrays of different rank
  function f_malloc_str_ptr(length,sizes,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloc-total-inc.f90'
    m%len=length
  end function f_malloc_str_ptr
  !> Define the allocation information for  arrays of different rank
  function f_malloc0_str_ptr(length,sizes,id,routine_id,lbounds,ubounds,profile) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloc-total-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloc0_str_ptr

  !> Define the allocation information for  arrays of different rank
  function f_malloc_i2(src,id,routine_id,profile) result(m)
    implicit none
    integer, dimension(:,:), intent(in) :: src
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_i2

  !> Define the allocation information for  arrays of different rank
  function f_malloc_ptr_i2(src,id,routine_id,profile) result(m)
    implicit none
    integer, dimension(:,:), intent(in) :: src
    type(malloc_information_ptr) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_ptr_i2

  !---Templates start here
  include 'malloc_templates-inc.f90'

end module dynamic_memory
