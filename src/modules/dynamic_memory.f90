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
  logical :: routine_closed=.false.
  character(len=namelen) :: present_routine=repeat(' ',namelen)
  
  !parameters for defitions of internal dictionary
  character(len=*), parameter :: arrayid='Array Id'
  character(len=*), parameter :: routineid='Allocating Routine Id'
  character(len=*), parameter :: sizeid='Size (Bytes)'

  !> Structure needed to allocate an allocatable array
  type, public :: malloc_information_all
     logical :: try !<raise an exception
     integer :: rank !< rank of the array
     integer, dimension(7) :: shape,lbounds,ubounds
     character(len=namelen) :: array_id !< label the array
     character(len=namelen) :: routine_id !<label the routine
  end type malloc_information_all
  
  type, public :: array_bounds
     integer :: nlow !<lower bounds
     integer :: nhigh !<higher bounds
  end type array_bounds

  interface assignment(=)
     module procedure i1_all,d1_all,d2_all,d3_all
  end interface

  interface operator(.to.)
     module procedure bounds
  end interface

  interface f_free
     module procedure i1_all_free,d1_all_free,d2_all_free
  end interface

  interface pad_with_nan
     module procedure i_padding,dp_padding,c_padding,l_padding,sp_padding,dp_padding2,dp_padding3
  end interface

  public :: f_malloc_set_status,f_malloc_finalize
  public :: f_malloc,f_free,f_malloc_routine_id
  public :: assignment(=),operator(.to.)

  contains

    function bounds(nlow,nhigh)
      implicit none
      integer, intent(in) :: nlow,nhigh
      type(array_bounds) :: bounds
      
      bounds%nlow=nlow
      bounds%nhigh=nhigh
    end function bounds


    function malloc_information_all_null() result(m)
      implicit none
      type(malloc_information_all) :: m
      !local variables
      integer :: i

      m%try=.false.
      m%rank=1
      m%shape=0
      m%lbounds=1
      m%ubounds=0
      do i=1,namelen
         m%array_id(i:i)=' '
         m%routine_id(i:i)=' '
      end do

    end function malloc_information_all_null


    !define the allocation information for  arrays of different rank
    function f_malloc(shape,id,routine_id,lbounds,ubounds,bounds,try) result(m)
      implicit none
      integer, dimension(:), intent(in), optional :: shape,lbounds,ubounds
      logical, intent(in), optional :: try
      character(len=*), intent(in), optional :: id,routine_id
      type(array_bounds), dimension(:), optional :: bounds
      type(malloc_information_all) :: m
      !local variables
      integer :: lgt,i
      m=malloc_information_all_null()
      !guess the rank
      m%rank=0

      if (present(bounds)) then
         m%rank=size(bounds)
         do i=1,m%rank
            m%lbounds(i)=bounds(i)%nlow
            m%ubounds(i)=bounds(i)%nhigh
            m%shape(i)=m%ubounds(i)-m%lbounds(i)+1
         end do
      end if

      if (present(lbounds)) then
         m%rank=size(lbounds)
         m%lbounds(1:m%rank)=lbounds
      end if
      
      if (present(shape)) then
         if (m%rank == 0) then
            m%rank=size(shape)
         else if (m%rank/=size(shape)) then
            stop 'ERROR, f_malloc: shape not conformal with lbounds'
         end if
         m%shape(1:m%rank)=shape
         do i=1,m%rank
            m%ubounds(i)=m%lbounds(i)+m%shape(i)-1
         end do
         if (present(ubounds)) then
            if (m%rank/=size(ubounds)) stop 'ERROR, f_malloc: shape not conformal with ubounds'
            do i=1,m%rank
               if (m%ubounds(i) /=ubounds(i)) stop 'ERROR, f_malloc: ubounds not conformal with shape and lbounds'
            end do
         end if
      else
         if (present(ubounds)) then
            if (m%rank == 0) then
               m%rank=size(shape)
            else if (m%rank/=size(ubounds)) then
               stop 'ERROR, f_malloc: ubounds not conformal with lbounds or shape'
            end if
            m%ubounds(1:m%rank)=ubounds
            do i=1,m%rank
               m%shape(i)=m%ubounds(i)-m%lbounds(i)+1
            end do
         else
            if (.not. present(bounds)) stop 'ERROR, f_malloc: at least shape or ubounds should be defined'
         end if
      end if

      if (present(id)) then
         lgt=min(len(id),namelen)
         m%array_id(1:lgt)=id(1:lgt)
      end if
      if (present(routine_id)) then
         lgt=min(len(routine_id),namelen)
         m%routine_id(1:lgt)=routine_id(1:lgt)
         call f_malloc_routine_id(m%routine_id)
      else
         m%routine_id=present_routine
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

    subroutine f_malloc_routine_id(routine_id)
      implicit none
      character(len=*), intent(in) :: routine_id
      !local variables
      integer :: lgt
      if (trim(present_routine) /= trim(routine_id)) then
         if(len_trim(present_routine)/=0) then
            routine_closed=.true.
         else
            !this means that we are at the initialization
            call dict_init(dict_routine)
         end if
         present_routine=repeat(' ',namelen)
         lgt=min(len(routine_id),namelen)
         present_routine(1:lgt)=routine_id(1:lgt)
      end if
    end subroutine f_malloc_routine_id

    !>initialize the library
    subroutine f_malloc_set_status(memory_limit,output_level,logfile_name,unit)
      use yaml_output, only: yaml_date_and_time_toa
      implicit none
      character(len=*), intent(in), optional :: logfile_name
      real(kind=4), intent(in), optional :: memory_limit
      integer, intent(in), optional :: output_level,unit

      if (.not. profile_initialized) then
         profile_initialized=.true.
         !initalize the dictionary with the allocation information
         call dict_init(dict_global)
         call set(dict_global//'Timestamp of Profile initialization',trim(yaml_date_and_time_toa()))
      end if
      
      if (present(memory_limit)) call memocc_set_memory_limit(memory_limit)

      if (present(output_level)) call memocc_set_state(output_level)

      if (present(unit)) call memocc_set_stdout(unit)

      if (present(logfile_name)) call memocc_set_filename(logfile_name)
      
    end subroutine f_malloc_set_status

    subroutine f_malloc_finalize()
      use yaml_output, only: yaml_warning,yaml_open_map,yaml_close_map,yaml_dict_dump,yaml_get_default_stream
      implicit none
      !local variables
      integer :: unt
      !put the last values in the dictionary if not freed
      if (associated(dict_routine)) then
         call yaml_get_default_stream(unt)
         !call yaml_stream_attributes(unit=unt)
         call yaml_warning('Not all the arrays have been freed: memory leaks are possible')
         call prepend(dict_global,dict_routine)
!      end if
!      if (.false.) then !residual memory to be defined
         call yaml_open_map('Addresses not being deallocated')
         call yaml_dict_dump(dict_global)
         call yaml_close_map()
      end if
      call dict_free(dict_global)
      call memocc_report()
      profile_initialized=.false.
      present_routine=repeat(' ',namelen)
      routine_closed=.false.
    end subroutine f_malloc_finalize

    subroutine check_for_errors(ierror,try)
      use yaml_output, only: yaml_warning,yaml_open_map,yaml_close_map,yaml_dict_dump,yaml_get_default_stream
      implicit none
      logical, intent(in) :: try
      integer, intent(in) :: ierror
      !local variables
      integer :: unt

      !recuperate possible error
      if (ierror /= INVALID_RANK) lasterror='Fortran (de)allocation problem'

      !raise exception if exit
      if (ierror /= SUCCESS) then
         if (try) then
            return
         else
            write(*,*)'(de)allocation error, exiting. Error code:',ierror
            write(*,*)'last error:',lasterror
            call yaml_get_default_stream(unt)
            call yaml_open_map('Status of the routine before exiting')
            call yaml_dict_dump(dict_routine)
            call yaml_close_map()
            stop
         end if
      end if

    end subroutine check_for_errors
    
    !> use the adress of the allocated pointer to profile the deallocation
    subroutine profile_deallocation(ierr,ilsize,address)
      use yaml_output, only: yaml_warning,yaml_open_map,yaml_close_map,yaml_dict_dump
      implicit none
      integer, intent(in) :: ierr
      integer(kind=8), intent(in) :: ilsize
      character(len=*), intent(in) :: address
      !local variables
      character(len=namelen) :: array_id,routine_id
      integer(kind=8) :: jlsize
      type(dictionary), pointer :: dict_add

      call check_for_errors(ierr,.false.)
      !search in the dictionaries the address
      dict_add=>find_key(dict_routine,trim(address))
      if (.not. associated(dict_add)) dict_add=>find_key(dict_global,trim(address))
      if (.not. associated(dict_add)) stop 'profile deallocations: address not present'
      !the global dictionary should be used instead
      array_id=dict_add//arrayid
      routine_id=dict_add//routineid
      jlsize=dict_add//sizeid

      call memocc(ierr,-int(ilsize),trim(arrayid),trim(routineid))

      call pop(dict_routine,trim(address))
      
    end subroutine profile_deallocation


    subroutine profile_allocation(ierr,address,sizeof,m)
      implicit none
      integer, intent(in) :: ierr,sizeof
      character(len=*), intent(in) :: address
      type(malloc_information_all), intent(in) :: m
      !local variables
      integer :: i
      integer(kind=8) :: ilsize
      type(dictionary), pointer :: dict_tmp
      
      !finalize the routine
      !if (trim(present_routine) /= trim(m%routine_id)) then
      if (routine_closed) then
         !if (len_trim(present_routine)/=0) then
            call prepend(dict_global,dict_routine)
      !      present_routine=m%routine_id
         !end if
         call dict_init(dict_routine)
      end if
      !size
      ilsize=int(sizeof,kind=8)
      do i=1,m%rank
         ilsize=ilsize*int(m%shape(i),kind=8)
      end do
      !create the dictionary array
      !add the array to the routine
      call dict_array(m%routine_id,m%array_id,ilsize,dict_tmp)
      call set(dict_routine//trim(address),dict_tmp)
      !call set(dict_routine//trim(address),dict_array(m%routine_id,m%array_id,ilsize))
      call check_for_errors(ierr,m%try)
      call memocc(ierr,product(m%shape(1:m%rank))*sizeof,m%array_id,m%routine_id)
      contains
        
        subroutine dict_array(routine_id,array_id,size,dict_tmp)
          implicit none
          character(len=*), intent(in) :: array_id,routine_id
          integer(kind=8), intent(in) :: size !< in bytes
          type(dictionary), pointer :: dict_tmp
          nullify(dict_tmp)
          call dict_init(dict_tmp)
          call set(dict_tmp//arrayid,trim(array_id))
          call set(dict_tmp//sizeid,size)
          call set(dict_tmp//routineid,trim(routine_id))

        end subroutine dict_array

!!$        function dict_array(routine_id,array_id,size)
!!$          implicit none
!!$          character(len=*), intent(in) :: array_id,routine_id
!!$          integer(kind=8), intent(in) :: size !< in bytes
!!$          type(dictionary), pointer :: dict_array
!!$          nullify(dict_array)
!!$print *,'test',associated(dict_array)
!!$          call dict_init(dict_array)
!!$print *,'test',associated(dict_array)
!!$          call set(dict_array//arrayid,trim(array_id))
!!$          call set(dict_array//sizeid,size)
!!$          call set(dict_array//routineid,trim(routine_id))
!!$
!!$        end function dict_array


    end subroutine profile_allocation

    function rank_is_ok(test,ref) result(ok)
      integer, intent(in) :: test,ref
      logical :: ok
      ierror=SUCCESS
      ok=test==ref
      if (.not. ok) then
         ierror=INVALID_RANK
         lasterror='rank not valid'
      end if
    end function rank_is_ok


    subroutine i1_all(array,m)
      implicit none
      integer, dimension(:), allocatable, intent(inout) :: array
      type(malloc_information_all), intent(in) :: m
      !local variables
      integer :: ierr,npaddim
      character(len=info_length) :: address
!      call timing(0,'AllocationProf','IR') 
      if (rank_is_ok(size(shape(array)),m%rank)) then
         !fortran allocation
         allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
         call pad_with_nan(array,m%rank,m%shape)
         !profile the array allocation
         call getaddress(array,address,len(address),ierr)
         call profile_allocation(ierror,address,kind(array),m)
      else
         call check_for_errors(ierror,m%try)
      end if
!      call timing(0,'AllocationProf','RS') 
    end subroutine i1_all

    subroutine d1_all(array,m)
      implicit none
      double precision, dimension(:), allocatable, intent(inout) :: array
      type(malloc_information_all), intent(in) :: m
      !local variables
      integer :: istat,ierr,npaddim
      character(len=info_length) :: address
!      call timing(0,'AllocationProf','IR') 
      if (rank_is_ok(size(shape(array)),m%rank)) then
         !fortran allocation
         allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
         call pad_with_nan(array,m%rank,m%shape)
         !profile the array allocation
         call getaddress(array,address,len(address),ierr)
         call profile_allocation(ierror,address,kind(array),m)
      else
         call check_for_errors(ierror,m%try)
      end if
!      call timing(0,'AllocationProf','RS') 
    end subroutine d1_all

    subroutine d2_all(array,m)
      implicit none
      double precision, dimension(:,:), allocatable, intent(inout) :: array
      type(malloc_information_all), intent(in) :: m
      !local variables
      integer :: istat,ierr,npaddim,i
      character(len=info_length) :: address
!      call timing(0,'AllocationProf','IR') 
      if (rank_is_ok(size(shape(array)),m%rank)) then
         !fortran allocation
         allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
         call pad_with_nan(array,m%rank,m%shape)
         !profile the array allocation
         call getaddress(array,address,len(address),ierr)
         call profile_allocation(ierror,address,kind(array),m)
      else
         call check_for_errors(ierror,m%try)
      end if
!      call timing(0,'AllocationProf','RS') 
    end subroutine d2_all

    subroutine d3_all(array,m)
      implicit none
      double precision, dimension(:,:,:), allocatable, intent(inout) :: array
      type(malloc_information_all), intent(in) :: m
      !local variables
      integer :: istat,ierr,npaddim
      character(len=info_length) :: address
!      call timing(0,'AllocationProf','IR') 
      if (rank_is_ok(size(shape(array)),m%rank)) then
         !fortran allocation
         allocate(array(m%lbounds(1):m%ubounds(1),&
              m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
         call pad_with_nan(array,m%rank,m%shape)
         !profile the array allocation
         call getaddress(array,address,len(address),ierr)
         call profile_allocation(ierror,address,kind(array),m)
      else
         call check_for_errors(ierror,m%try)
      end if
!      call timing(0,'AllocationProf','RS') 
    end subroutine d3_all


    subroutine i1_all_free(array)
      implicit none
      integer, dimension(:), allocatable, intent(inout) :: array
      include 'deallocate-inc.f90' 
!!$      !local variables
!!$      integer :: ierr
!!$      character(len=info_length) :: address
!!$      !local variables
!!$      integer :: i_all
!!$      integer(kind=8) :: ilsize
!!$!      call timing(0,'AllocationProf','IR') 
!!$      !profile the array allocation
!!$      call getaddress(array,address,len(address),ierr)
!!$      ilsize=int(product(shape(array))*kind(array),kind=8)
!!$      !fortran deallocation
!!$      deallocate(array,stat=ierror)
!!$      !hopefully only address is necessary for the deallocation
!!$      call profile_deallocation(ierror,ilsize,address)

!!$  i_all=-product(shape(array))*kind(array)
!!$  deallocate(array,stat=ierror)
!!$  call memocc(ierror,i_all,'stuff','dosome')
!      call timing(0,'AllocationProf','RS') 
    end subroutine i1_all_free


    subroutine d1_all_free(array)
      implicit none
      double precision, dimension(:), allocatable, intent(inout) :: array
      include 'deallocate-inc.f90' 
!!$      !local variables
!!$      integer :: istat,ierr
!!$      character(len=info_length) :: address
!!$      !local variables
!!$      integer :: i_all
!!$      integer(kind=8) :: ilsize
!!$!      call timing(0,'AllocationProf','IR') 
!!$      !profile the array allocation
!!$      call getaddress(array,address,len(address),ierr)
!!$      ilsize=int(product(shape(array))*kind(array),kind=8)
!!$      !fortran deallocation
!!$      deallocate(array,stat=ierror)
!!$      !hopefully only address is necessary for the deallocation
!!$      call profile_deallocation(ierror,ilsize,address)

!!$  i_all=-product(shape(array))*kind(array)
!!$  deallocate(array,stat=ierror)
!!$  call memocc(ierror,i_all,'stuff','dosome')
!      call timing(0,'AllocationProf','RS') 
    end subroutine d1_all_free

    subroutine d2_all_free(array)
      implicit none
      double precision, dimension(:,:), allocatable, intent(inout) :: array
      include 'deallocate-inc.f90' 
!!$      !local variables
!!$      integer :: istat,ierr
!!$      character(len=info_length) :: address
!!$      !local variables
!!$      integer :: i_all
!!$      integer(kind=8) :: ilsize
!!$!      call timing(0,'AllocationProf','IR') 
!!$      !profile the array allocation
!!$      call getaddress(array,address,len(address),ierr)
!!$      !here the size should be corrected with ndebug
!!$      ilsize=int(product(shape(array))*kind(array),kind=8)
!!$      !fortran deallocation
!!$      deallocate(array,stat=ierror)
!!$      !hopefully only address is necessary for the deallocation
!!$      call profile_deallocation(ierror,ilsize,address)
!!$!      call timing(0,'AllocationProf','RS') 
    end subroutine d2_all_free



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
    subroutine dp_padding(array,nrank,arr_shape)
      implicit none
      integer, intent(in) :: nrank
      integer, dimension(nrank), intent(in) :: arr_shape
      double precision, dimension(*), intent(inout) :: array
      !local variables
      integer :: i,npaddim,nstart
      npaddim=1
      do i=1,nrank-1
         npaddim=npaddim*arr_shape(i)
      end do
      nstart=npaddim*arr_shape(nrank)
      do i=1,npaddim*ndebug
         array(nstart+i)= d_nan()
      end do
    end subroutine dp_padding

    subroutine dp_padding2(array,nrank,arr_shape)
      implicit none
      integer, intent(in) :: nrank
      integer, dimension(nrank), intent(in) :: arr_shape
      double precision, dimension(arr_shape(1),*), intent(inout) :: array
      !call the padding routine as a rank-1 array
      call dp_padding(array,1,(/product(arr_shape)/))
    end subroutine dp_padding2

    subroutine dp_padding3(array,nrank,arr_shape)
      implicit none
      integer, intent(in) :: nrank
      integer, dimension(nrank), intent(in) :: arr_shape
      double precision, dimension(arr_shape(1),arr_shape(2),*), intent(inout) :: array
      !call the padding routine as a rank-1 array
      call dp_padding(array,1,(/product(arr_shape)/))
    end subroutine dp_padding3

    subroutine sp_padding(array,nrank,arr_shape)
      implicit none
      integer, intent(in) :: nrank
      integer, dimension(nrank), intent(in) :: arr_shape
      real, dimension(*), intent(inout) :: array
      !local variables
      integer :: i,npaddim,nstart
      npaddim=1
      do i=1,nrank-1
         npaddim=npaddim*arr_shape(i)
      end do
      nstart=npaddim*arr_shape(nrank)
      do i=1,npaddim*ndebug
         array(nstart+i)= r_nan()
      end do
    end subroutine sp_padding

    subroutine i_padding(array,nrank,arr_shape)
      implicit none
      integer, intent(in) :: nrank
      integer, dimension(nrank), intent(in) :: arr_shape
      integer, dimension(*), intent(inout) :: array
      !local variables
      integer :: i,npaddim,nstart
      npaddim=1
      do i=1,nrank-1
         npaddim=npaddim*arr_shape(i)
      end do
      nstart=npaddim*arr_shape(nrank)
      do i=1,npaddim*ndebug
         array(nstart+i)= int(r_nan()) !this function is in profiling/timem.f90
      end do
    end subroutine i_padding

!!$    subroutine i_padding(npaddim,nstart,array)
!!$      implicit none
!!$      integer, intent(in) :: npaddim,nstart
!!$      integer, dimension(*) :: array
!!$      !local variables
!!$      integer :: i
!!$      do i=1,npaddim*ndebug
!!$         array(nstart+i)= int(r_nan()) !this function is in profiling/timem.f90
!!$      end do
!!$    end subroutine i_padding


    subroutine l_padding(array,nrank,arr_shape)
      implicit none
      integer, intent(in) :: nrank
      integer, dimension(nrank), intent(in) :: arr_shape
      logical, dimension(*), intent(inout) :: array
      !local variables
      integer :: i,npaddim,nstart
      npaddim=1
      do i=1,nrank-1
         npaddim=npaddim*arr_shape(i)
      end do
      nstart=npaddim*arr_shape(nrank)
      do i=1,npaddim*ndebug
         array(nstart+i)= .false.
      end do
    end subroutine l_padding

    subroutine c_padding(array,nrank,arr_shape)
      implicit none
      integer, intent(in) :: nrank
      integer, dimension(nrank), intent(in) :: arr_shape
      character(len=*), dimension(*), intent(inout) :: array
      !local variables
      integer :: i,npaddim,nstart
      npaddim=1
      do i=1,nrank-1
         npaddim=npaddim*arr_shape(i)
      end do
      nstart=npaddim*arr_shape(nrank)
      do i=1,npaddim*ndebug
         array(nstart+i)=repeat('A',len(array(1)))!)'AAAAAAAAAAAAAAAAAAAA'
      end do
    end subroutine c_padding


end module dynamic_memory
