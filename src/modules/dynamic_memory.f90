module dynamic_memory
  use m_profiling, except => ndebug, and=> d_nan, also=> r_nan
  use dictionaries, info_length => max_field_length
  implicit none
  
  private 
  
  !> length of the character variables
  integer, parameter :: namelen=32
  !> length of error string
  integer, parameter :: error_string_len=80
  !> size of debug parameters
  integer, parameter :: ndebug=0
  !> errorcodes
  character(len=error_string_len) :: lasterror=repeat(' ',len(lasterror))
  integer, parameter :: SUCCESS                = 0
  integer, parameter :: INVALID_RANK           = -1979
  integer :: ierror=SUCCESS

  !global variables for initialization
  logical :: profile_initialized=.false.
  !dictionaries needed for profiling storage
  type(dictionary), pointer :: dict_global,dict_routine
  !global variable (can be stored in dictionaries)
  character(len=namelen) :: present_routine=repeat(' ',namelen)
  
  !parameters for defitions of internal dictionary
  character(len=*), parameter :: arrayid='Array Id'
  character(len=*), parameter :: routineid='Allocating Routine Id'
  character(len=*), parameter :: sizeid='Size (Bytes)'

  !> Structure needed to allocate an allocatable array
  type, public :: malloc_information_all
     logical :: try !<raise an exception
     integer :: rank !< rank of the array
     integer, dimension(7) :: shape 
     character(len=namelen) :: array_id !< label the array
     character(len=namelen) :: routine_id !<label the routine
  end type malloc_information_all
  
  interface assignment(=)
     module procedure i1_all
  end interface

  interface f_free
     module procedure i1_all_free
  end interface

  public :: f_malloc_set_status,f_malloc_finalize
  public :: f_malloc,f_free
  public :: assignment(=)

  contains

    function malloc_information_all_null() result(m)
      implicit none
      type(malloc_information_all) :: m
      !local variables
      integer :: i

      m%try=.false.
      m%rank=1
      m%shape=0
      do i=1,namelen
         m%array_id(i:i)=' '
         m%routine_id(i:i)=' '
      end do

    end function malloc_information_all_null

    !define the allocation information for  arrays of different rank
    function f_malloc(shape,array_id,routine_id,try) result(m)
      implicit none
      integer, dimension(:), intent(in) :: shape
      type(malloc_information_all) :: m
      logical, intent(in), optional :: try
      character(len=*), intent(in), optional :: array_id,routine_id
      !local variables
      integer :: lgt

      m=malloc_information_all_null()

      !associate the rank of the array
      m%rank=size(shape)
      m%shape(1:m%rank)=shape

      if (present(array_id)) then
         lgt=min(len(array_id),namelen)
         m%array_id(1:lgt)=array_id(1:lgt)
      end if
      if (present(routine_id)) then
         lgt=min(len(routine_id),namelen)
         m%routine_id(1:lgt)=routine_id(1:lgt)
      end if

      if(present(try)) m%try=try
    end function f_malloc

    function last_f_malloc_error(error_string) result(ierr)
      implicit none
      integer :: ierr
      character(len=*), intent(out), optional :: error_string
      !local variables
      integer :: lgt,i

      ierr=ierror
      
      if (present(error_string) .and. ierr/=SUCCESS) then
         lgt=min(len(error_string),error_string_len)
         error_string(1:lgt)=lasterror(1:lgt)
         do i=lgt+1,error_string_len
            error_string(i:i)=' '
         end do
         !clean last error
         lasterror=repeat(' ',len(lasterror))
      end if

    end function last_f_malloc_error

    !>initialize the library
    subroutine f_malloc_set_status(memory_limit,output_level,logfile_name)
      use yaml_output, only: yaml_date_and_time_toa
      implicit none
      character(len=*), intent(in), optional :: logfile_name
      real(kind=4), intent(in), optional :: memory_limit
      integer, intent(in), optional :: output_level

      if (.not. profile_initialized) then
         profile_initialized=.true.
         !initalize the dictionary with the allocation information
         call dict_init(dict_global)
         call set(dict_global/'Timestamp of Profile initialization',trim(yaml_date_and_time_toa()))
      end if
      
      if (present(memory_limit)) then
         call memocc_set_memory_limit(memory_limit)
      end if
      
    end subroutine f_malloc_set_status

    subroutine f_malloc_finalize()
      use yaml_output, only: yaml_warning,yaml_open_map,yaml_close_map
      !put the last values in the dictionalry if not freed
      if (associated(dict_routine)) then
         call yaml_warning('Not all the arrays have been freed: memory leaks are possible')
         call prepend(dict_global,dict_routine)
         call yaml_open_map('Addresses not being deallocated')
         call dictionary_print(dict_global)
         call yaml_close_map()
      end if
      call dict_free(dict_global)
      call memocc_report()
      profile_initialized=.false.
      present_routine=repeat(' ',namelen)
      
    end subroutine f_malloc_finalize

    !> use the adress of the allocated pointer to profile the deallocation
    subroutine profile_deallocation(ierror,ilsize,address)
      implicit none
      integer, intent(in) :: ierror
      integer(kind=8), intent(in) :: ilsize
      character(len=*), intent(in) :: address
      !local variables
      character(len=namelen) :: array_id,routine_id
      integer(kind=8) :: jlsize
      type(dictionary), pointer :: dict_add

      !search in the dictionaries the address
      !here an exception should be raised to search for the address
      dict_add=>dict_routine/trim(address)
      !print *,'here',trim(address),arrayid
      !call dictionary_print(dict_add)
      array_id=dict_add//arrayid
      !print *,'here2',address
      routine_id=dict_add//routineid
      !print *,'here3',address
      jlsize=dict_add//sizeid

      call memocc(ierror,-int(ilsize),trim(arrayid),trim(routineid))

      call pop(dict_routine,trim(address))
      
    end subroutine profile_deallocation


    subroutine profile_allocation(ierror,address,sizeof,m)
      implicit none
      integer, intent(in) :: ierror,sizeof
      character(len=*), intent(in) :: address
      type(malloc_information_all), intent(in) :: m
      !local variables
      integer :: i
      integer(kind=8) :: ilsize
      
      !recuperate possible error
      if (ierror /= INVALID_RANK) lasterror='allocation problem'
      
      !raise exception if exit
      if (ierror /= SUCCESS) then
         if (m%try) then
            return
         else
            write(*,*)'allocation error, exiting. Error code:',ierror
            write(*,*)'last error:',lasterror
            stop
         end if
      end if      

      !finalize the routine
      if (trim(present_routine) /= trim(m%routine_id)) then
         if (len_trim(present_routine)/=0) then
            call prepend(dict_global,dict_routine)
            present_routine=m%routine_id
         end if
         call dict_init(dict_routine)
      end if
      !size
      ilsize=int(sizeof,kind=8)
      do i=1,m%rank
         ilsize=ilsize*int(m%shape(i),kind=8)
      end do
      !create the dictionary array
      !add the array to the routine
      call set(dict_routine/trim(address),dict_array(m%routine_id,m%array_id,ilsize))

      call memocc(ierror,product(m%shape(1:m%rank))*sizeof,m%array_id,m%routine_id)
      contains
        
        function dict_array(routine_id,array_id,size)
          implicit none
          character(len=*), intent(in) :: array_id,routine_id
          integer(kind=8), intent(in) :: size !< in bytes
          type(dictionary), pointer :: dict_array

          call dict_init(dict_array)

          call set(dict_array/arrayid,trim(array_id))
          call set(dict_array/sizeid,size)
          call set(dict_array/routineid,trim(routine_id))

        end function dict_array

    end subroutine profile_allocation
    

    !define the routines which allocate the different allocatable array
    subroutine i1_all(array,m)
      implicit none
      integer, dimension(:), allocatable, intent(inout) :: array
      type(malloc_information_all), intent(in) :: m
      !local variables
      integer :: ierr,npaddim
      character(len=info_length) :: address
!      call timing(0,'AllocationProf','IR') 
      !no error if everything is OK
      ierror=SUCCESS
      if (m%rank/=size(shape(array))) then
         ierror=INVALID_RANK
         lasterror='rank not valid'
      else
         !fortran allocation
         allocate(array(m%shape(1)+ndebug),stat=ierror)
         !padding of the array in case of ndebug
         npaddim=1
         if (m%rank > 1) npaddim=product(m%shape(1:m%rank-1))
         call i_padding(npaddim,product(m%shape),array)
      end if
      !profile the array allocation
      call getaddress(array,address,len(address),ierr)
      call profile_allocation(ierror,address,kind(array),m)
!      call timing(0,'AllocationProf','RS') 
    end subroutine i1_all

    !define the routines which deallocate the different allocatable array
    subroutine i1_all_free(array)
      implicit none
      integer, dimension(:), allocatable, intent(inout) :: array
      !local variables
      integer :: ierr
      character(len=info_length) :: address
      !local variables
!!$      integer :: i_all
      integer(kind=8) :: ilsize
!      call timing(0,'AllocationProf','IR') 
      !no error if everything is OK
      ierror=SUCCESS
      !profile the array allocation
      call getaddress(array,address,len(address),ierr)
      ilsize=int(product(shape(array))*kind(array),kind=8)
      !fortran deallocation
      deallocate(array,stat=ierror)
      !hopefully only address is necessary for the deallocation
      call profile_deallocation(ierror,ilsize,address)

!!$  i_all=-product(shape(array))*kind(array)
!!$  deallocate(array,stat=ierror)
!!$  call memocc(ierror,i_all,'stuff','dosome')
!      call timing(0,'AllocationProf','RS') 
    end subroutine i1_all_free



    !!****f* ABINIT/d_nan
    !! FUNCTION
    !!   Function which specify NaN according to IEEE specifications
    !! SOURCE
    !!
    function d_nan()
      implicit none
      double precision :: d_nan
      !local variables
      double precision :: dnan
      integer, dimension(2) :: inan
      equivalence (dnan, inan)
      ! This first assignment is for big-endian machines
      inan(1) = 2147483647
      ! The second assignment is for little-endian machines
      inan(2) = 2147483647
      d_nan = dnan
    end function d_nan
    !!***

    !!****f* ABINIT/r_nan
    !! FUNCTION
    !!   Function which specify NaN according to IEEE specifications
    !! SOURCE
    !!
    function r_nan()
      implicit none
      real :: r_nan
      !local variables
      real :: rnan
      integer :: inan
      equivalence (rnan, inan)
      inan = 2147483647
      r_nan = rnan
    end function r_nan
    !!***
    subroutine sp_padding(npaddim,nstart,array)
      implicit none
      integer, intent(in) :: npaddim,nstart
      real, dimension(*) :: array
      !local variables
      integer :: i
      do i=1,npaddim*ndebug
         array(nstart+i)= r_nan() !this function is in profiling/memory.f90
      end do
    end subroutine sp_padding

    subroutine i_padding(npaddim,nstart,array)
      implicit none
      integer, intent(in) :: npaddim,nstart
      integer, dimension(*) :: array
      !local variables
      integer :: i
      do i=1,npaddim*ndebug
         array(nstart+i)= int(r_nan()) !this function is in profiling/timem.f90
      end do
    end subroutine i_padding

    subroutine l_padding(npaddim,nstart,array)
      implicit none
      integer, intent(in) :: npaddim,nstart
      logical, dimension(*) :: array
      !local variables
      integer :: i
      do i=1,npaddim*ndebug
         array(nstart+i)=.false.
      end do
    end subroutine l_padding

    subroutine c_padding(npaddim,nstart,array)
      implicit none
      integer, intent(in) :: npaddim,nstart
      character(len=20), dimension(*) :: array
      !local variables
      integer :: i
      do i=1,npaddim*ndebug
         array(nstart+i)='AAAAAAAAAAAAAAAAAAAA'
      end do
    end subroutine c_padding


end module dynamic_memory
