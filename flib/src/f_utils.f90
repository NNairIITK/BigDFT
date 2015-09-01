!> @file
!! Manage different low-level operations
!! like operations on external files and basic operations in memory
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Manage low-level operations (external files and basic memory operations)
module f_utils
  use dictionaries, only: f_err_throw,f_err_define, &
       & dictionary, dict_len, dict_iter, dict_next, dict_value, max_field_length
  use yaml_strings, only: yaml_toa
  use f_precisions
  !use module_razero
  implicit none

  private

  integer, private, save :: INPUT_OUTPUT_ERROR

  integer, public, save :: TCAT_INIT_TO_ZERO

  !preprocessed include file with processor-specific values
  include 'f_utils.inc' !defines recl_kind

  !> This type can be used to get strings from a file or a dictionary long string.
  type, public :: io_stream
     integer :: iunit = 0
     type(dictionary), pointer :: lstring => null()
  end type io_stream

  !> This type can be used to dump strings at bunches in a file
  type, public :: f_dump_buffer
     integer :: ipos
     character(len=1), dimension(:), pointer :: buf
  end type f_dump_buffer


  !> Interface for difference between two intrinsic types
  interface f_diff
     module procedure f_diff_i,f_diff_r,f_diff_d,f_diff_li,f_diff_l
     module procedure f_diff_d2d3,f_diff_d2d1,f_diff_d1d2,f_diff_d2,f_diff_d1
     module procedure f_diff_i2i1,f_diff_i1,f_diff_i2,f_diff_i1i2
     module procedure f_diff_li2li1,f_diff_li1,f_diff_li2,f_diff_li1li2
     module procedure f_diff_d0d1,f_diff_i0i1, f_diff_li0li1
     module procedure f_diff_c1i1,f_diff_c0i1
     module procedure f_diff_c1li1,f_diff_c0li1
  end interface f_diff

  !> Initialize to zero an array (should be called f_memset)
  interface f_zero
     module procedure zero_string
     module procedure zero_li,zero_i,zero_r,zero_d,zero_l
     !module procedure put_to_zero_simple, put_to_zero_long
     module procedure put_to_zero_double, put_to_zero_double_1, put_to_zero_double_2
     module procedure put_to_zero_double_3, put_to_zero_double_4, put_to_zero_double_5
     module procedure put_to_zero_double_6, put_to_zero_double_7
     module procedure put_to_zero_integer,put_to_zero_integer1,put_to_zero_integer2
     module procedure put_to_zero_integer3
     module procedure put_to_zero_long,put_to_zero_long1,put_to_zero_long2
     module procedure put_to_zero_long3
  end interface f_zero

  !to be verified if clock_gettime is without side-effect, otherwise the routine cannot be pure
  interface
     pure subroutine nanosec(itime)
       use f_precisions, only: f_long
       implicit none
       integer(f_long), intent(out) :: itime
     end subroutine nanosec
  end interface

  public :: f_diff,f_file_unit,f_mkdir
  public :: f_utils_errors,f_utils_recl,f_file_exists,f_close,f_zero
  public :: f_get_free_unit,f_delete_file,f_getpid,f_rewind,f_open_file
  public :: f_iostream_from_file,f_iostream_from_lstring
  public :: f_iostream_get_line,f_iostream_release,f_time,f_pause

contains
 
  subroutine f_utils_errors()

    call f_err_define('INPUT_OUTPUT_ERROR',&
         'Some of intrinsic I/O fortran routines returned an error code.',&
         INPUT_OUTPUT_ERROR,&
         err_action='Check if you have correct file system permission in I/O library or check the fortran runtime library.')

  end subroutine f_utils_errors

  pure function f_time()
    integer(f_long) :: f_time
    !local variables
    integer(f_long) :: itime
    call nanosec(itime)
    f_time=itime
  end function f_time

  !>enter in a infinite loop for sec seconds. Use cpu_time as granularity is enough
  subroutine f_pause(sec,verbose)
    implicit none
    integer, intent(in) :: sec !< seconds to be waited
    logical, intent(in), optional :: verbose !<for debugging purposes, do not eliminate
    !local variables
    logical :: verb
    integer(kind=8) :: t0,t1
    integer :: count

    verb=.false.
    if (present(verbose)) verb=verbose

    if (sec <=0) return
    t0=f_time()
    t1=t0
    !this loop has to be modified to avoid the compiler to perform too agressive optimisations
    count=0
    do while(real(t1-t0,kind=8)*1.d-9 < real(sec,kind=8))
       count=count+1
       t1=f_time()
    end do
    !this output is needed to avoid the compiler to perform too agressive optimizations
    !therefore having a infinie loop
    if (verb) print *,'Paused for '//trim(yaml_toa(sec))//' seconds, counting:'//&
         trim(yaml_toa(count))
  end subroutine f_pause

  !> gives the maximum record length allowed for a given unit
  subroutine f_utils_recl(unt,recl_max,recl)
    implicit none
    integer, intent(in) :: unt !< unit to be checked for record length
    integer, intent(in) :: recl_max !< maximum value for record length
    !> Value for the record length. This corresponds to the minimum between recl_max and the processor-dependent value
    !! provided by inquire statement
    integer, intent(out) :: recl 
    !local variables
    logical :: unit_is_open
    integer :: ierr,ierr_recl
    integer(kind=recl_kind) :: recl_file

    !in case of any error, the value is set to recl_max
    recl=recl_max
    ierr_recl=-1
    !initialize the value of recl_file
    recl_file=int(-1234567891,kind=recl_kind)
    inquire(unit=unt,opened=unit_is_open,iostat=ierr)
    if (ierr == 0 .and. .not. unit_is_open) then
       !inquire the record length for the unit
       inquire(unit=unt,recl=recl_file,iostat=ierr_recl)
    end if
    if (ierr_recl == 0) then
       recl=int(min(int(recl_max,kind=recl_kind),recl_file))
    end if
    if (recl <=0) recl=recl_max
  end subroutine f_utils_recl

  !> inquire for the existence of a file
  subroutine f_file_exists(file,exists)
    implicit none
    character(len=*), intent(in) :: file
    logical, intent(out) :: exists
    !local variables
    integer :: ierr 

    exists=.false.
    inquire(file=trim(file),exist=exists,iostat=ierr)
    if (ierr /=0) then
       call f_err_throw('Error in inquiring file='//&
         trim(file)//', iostat='//trim(yaml_toa(ierr)),&
         err_id=INPUT_OUTPUT_ERROR)
    end if
    exists = exists .and. ierr==0

  end subroutine f_file_exists

  !> call the close statement and retrieve the error
  !! do not close the file if the unit is not connected
  subroutine f_close(unit)
    implicit none
    integer, intent(in) :: unit
    !local variables
    integer :: ierr

    if (unit > 0) then
       close(unit,iostat=ierr)
       if (ierr /= 0) call f_err_throw('Error in closing unit='//&
               trim(yaml_toa(unit))//', iostat='//trim(yaml_toa(ierr)),&
               err_id=INPUT_OUTPUT_ERROR)
    end if
  end subroutine f_close

  !> Search the unit associated to a filename.
  !! the unit is -1 if the file does not exists or if the file is
  !! not connected
  subroutine f_file_unit(file,unit)
    implicit none
    character(len=*), intent(in) :: file
    integer, intent(out) :: unit
    !local variables
    logical ::  exists
    integer :: unt,ierr

    unit=-1
    call f_file_exists(file,exists)
    if (exists) then
       inquire(file=trim(file),number=unt,iostat=ierr)
       if (ierr /= 0) then
          call f_err_throw('Error in inquiring file='//&
               trim(file)//' for number, iostat='//trim(yaml_toa(ierr)),&
               err_id=INPUT_OUTPUT_ERROR)
       else
          unit=unt
       end if
    end if
  end subroutine f_file_unit

  !> Get a unit which is not opened at present
  !! start the search from the unit
  function f_get_free_unit(unit) result(unt2)
    implicit none
    !> putative free unit. Starts to search from this value
    integer, intent(in), optional :: unit
    integer :: unt2
    !local variables
    logical :: unit_is_open
    integer :: unt,ierr

    unit_is_open=.true.
    unt=7
    if (present(unit)) unt=unit
    inquire(unit=unt,opened=unit_is_open,iostat=ierr)
    do while(unit_is_open .and. ierr==0)     
       unt=unt+1
       inquire(unit=unt,opened=unit_is_open,iostat=ierr)
    end do
    if (ierr /=0) then
       call f_err_throw('Error in inquiring unit='//&
            trim(yaml_toa(unt))//', iostat='//trim(yaml_toa(ierr)),&
            err_id=INPUT_OUTPUT_ERROR)
    end if
    unt2=unt
  end function f_get_free_unit

  !> Create a directory from CWD path
  subroutine f_mkdir(dir,path)
    use f_precisions, only: f_integer
    implicit none
    character(len=*), intent(in) :: dir !<directory to be created
    character(len=*), intent(out) :: path !<path of the created directory (trailing slash added)
    !local variables
    integer(f_integer) :: ierr
    integer(f_integer) :: lin,lout

    call f_zero(path)
    lin=int(len_trim(dir),f_integer)
    lout=int(len(path),f_integer)

    call getdir(dir,lin,path,lout,ierr)
    if (ierr /= 0 .and. ierr /= 1) then
       call f_err_throw('Error in creating directory ='//&
            trim(dir)//', iostat='//trim(yaml_toa(ierr)),&
            err_id=INPUT_OUTPUT_ERROR)
    end if

  end subroutine f_mkdir

  subroutine f_delete_file(file)
    implicit none
    character(len=*), intent(in) :: file
    !local variables
    logical :: exists
    integer :: ierr,unit
    external :: delete

    call f_file_exists(trim(file),exists)
    if (exists) then
       !close the corresponding fortran unit if the file is connected to it
       call f_file_unit(trim(file),unit)
       call f_close(unit)
       !c-function in utils.c
       call delete(trim(file),len_trim(file),ierr)
       if (ierr /=0) call f_err_throw('Error in deleting file='//&
            trim(file)//'iostat='//trim(yaml_toa(ierr)),&
            err_id=INPUT_OUTPUT_ERROR)
    end if
    
  end subroutine f_delete_file

  !> get process id
  function f_getpid()
    implicit none
    integer :: f_getpid
    !local variables
    integer :: pid
    external :: getprocid

    call getprocid(pid)
    f_getpid=pid

  end function f_getpid

  !> rewind a unit
  subroutine f_rewind(unit)
    implicit none
    integer, intent(in) :: unit
    !local variables
    integer :: ierr

    rewind(unit,iostat=ierr)
    if (ierr /=0) call f_err_throw('Error in rewind unit='//&
         trim(yaml_toa(unit))//'iostat='//trim(yaml_toa(ierr)),&
         err_id=INPUT_OUTPUT_ERROR)
    
  end subroutine f_rewind

  !>tentative example of writing the data in a buffer
  subroutine f_write(unit,msg,advance,buffer)
    use f_precisions, only: cr => f_cr
    use yaml_strings
    !use dynamic_memory, only: f_memcpy
    implicit none
    integer, intent(in) :: unit
    character(len=*), intent(in) :: msg
    character(len=*), intent(in), optional :: advance
    type(f_dump_buffer), optional, intent(inout) :: buffer
    !local variables
    integer :: lpos
    character(len=3) :: adv
    character(len=len(cr)) :: crtmp
    
    adv='yes'
    if (present(advance)) call f_strcpy(src=advance,dest=adv)

    if (present(buffer)) then
       !determine the size of the input
       lpos=len(msg)
       call f_zero(crtmp)
       if (adv .eqv. 'yes') crtmp=cr 
       lpos=lpos+len_trim(cr)
       !copy the values we would like to add in the buffer
       !check if the total length is bigger than buffer size
       if (lpos+buffer%ipos > size(buffer%buf)) then
          write(unit=unit,fmt='(a)') buffer%buf(:buffer%ipos)
          buffer%ipos=1
       end if
       !copy the data 
       !call f_memcpy(n=lpos,src=msg+crtmp,dest=buffer%buf(buffer%ipos))
       buffer%ipos=buffer%ipos+lpos
    else
       !we should inquire if the unit is formatted or not
       write(unit=unit,fmt='(a)',advance=adv) msg
    end if
  end subroutine f_write
  
  !> open a filename and retrieve the unteger for the unit
  subroutine f_open_file(unit,file,status,position,action,binary)
    use yaml_strings, only: f_strcpy
    implicit none
    !> integer of the unit. On entry, it indicates the 
    !! suggested unit number. On exit, it indicates the free unit
    !! which has been used for the file opening
    integer, intent(inout) :: unit
    !> filename
    character(len=*), intent(in) :: file
    !> status
    character(len=*), intent(in), optional :: status
    !> position
    character(len=*), intent(in), optional :: position
    !> action
    character(len=*), intent(in), optional :: action
    !> if true, the file will be opened in the unformatted i/o
    !! if false or absent, the file will be opened for formatted i/o
    logical, intent(in), optional :: binary
    !local variables
    integer :: unt,ierror
    character(len=7) :: f_status
    character(len=11) :: f_form
    character(len=6) :: f_position
    character(len=9) :: f_action

    !first, determine if the file is already opened.
    call f_file_unit(file,unt)
    if (unt /= -1) then
       unit=unt
    else
       !find the first free unit
       unt=f_get_free_unit(unit)

       !useful open specifiers
       call f_strcpy(src='unknown',dest=f_status)
       if (present(status)) call f_strcpy(src=status,dest=f_status)

       call f_strcpy(src='formatted',dest=f_form)
       if (present(binary)) then
          if (binary) call f_strcpy(src='unformatted',dest=f_form)
       end if

       call f_strcpy(src='asis',dest=f_position)
       if (present(position)) call f_strcpy(src=position,dest=f_position)

       call f_strcpy(src='readwrite',dest=f_action)
       if (present(action)) call f_strcpy(src=action,dest=f_action)

       !then open the file with the given unit
       open(unit=unt,file=trim(file),status=f_status,form=f_form,&
            position=f_position,action=f_action,iostat=ierror)
       if (ierror /= 0) then
          call f_err_throw('Error in opening file='//&
               trim(file)//' with unit='//trim(yaml_toa(unt,fmt='(i0)'))//&
               ', iostat='//trim(yaml_toa(ierror)),&
               err_id=INPUT_OUTPUT_ERROR)
       else
          !when everything succeded, assign the unit
          unit=unt
       end if
    end if

  end subroutine f_open_file

  subroutine f_iostream_from_file(ios, filename)
    implicit none
    type(io_stream), intent(out) :: ios
    character(len = *), intent(in) :: filename
    !Local variables
    integer :: ierror

    ios%iunit=f_get_free_unit(742)
    open(unit=ios%iunit,file=trim(filename),status='old',iostat=ierror)
    !Check the open statement
    if (ierror /= 0) call f_err_throw('Error in opening file='//&
         trim(filename)//' iostat='//trim(yaml_toa(ierror)),&
         err_id=INPUT_OUTPUT_ERROR)
    nullify(ios%lstring)
  end subroutine f_iostream_from_file

  subroutine f_iostream_from_lstring(ios, dict)
    implicit none
    type(io_stream), intent(out) :: ios
    type(dictionary), pointer :: dict

    ios%iunit = 0
    if (dict_len(dict) < 0) call f_err_throw('Error dict is not a long string',&
         err_id=INPUT_OUTPUT_ERROR)
    ios%lstring => dict
  end subroutine f_iostream_from_lstring

  subroutine f_iostream_get_line(ios, line, eof)
    implicit none
    !Arguments
    type(io_stream), intent(inout) :: ios
    character(len=max_field_length), intent(out) :: line
    logical, optional, intent(out) :: eof
    !Local variables
    integer :: i_stat
    character(len=8) :: fmt

    if (ios%iunit > 0) then
       write(fmt, "(A,I0,A)") "(A", max_field_length, ")"
       if (present(eof)) then
          read(ios%iunit, trim(fmt), iostat = i_stat) line
       else
          read(ios%iunit, trim(fmt)) line
          i_stat = 0
       end if
       if (i_stat /= 0) then
          close(ios%iunit)
          ios%iunit = 0
       end if
       if (present(eof)) eof = (i_stat /= 0)
    else if (associated(ios%lstring)) then
       if (dict_len(ios%lstring) > 0) ios%lstring => dict_iter(ios%lstring)
       line = dict_value(ios%lstring)
       ios%lstring => dict_next(ios%lstring)
       if (present(eof)) eof = .not. associated(ios%lstring)
    else if (present(eof)) then
       eof = .true.
    end if
  end subroutine f_iostream_get_line

  subroutine f_iostream_release(ios)
    implicit none
    type(io_stream), intent(inout) :: ios

    if (ios%iunit > 0) close(ios%iunit)
    ios%iunit = 0
    nullify(ios%lstring)
  end subroutine f_iostream_release

  !>perform a difference of two objects (of similar kind)
  subroutine f_diff_i(n,a_add,b_add,diff)
    implicit none
    integer, intent(in) :: n
    integer(kind=4), intent(inout) :: a_add
    integer(kind=4), intent(inout) :: b_add
    integer(kind=4), intent(out) :: diff
    external :: diff_i
    call diff_i(n,a_add,b_add,diff)
  end subroutine f_diff_i
  subroutine f_diff_i2i1(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    integer(kind=4), dimension(:,:),   intent(in) :: a
    integer(kind=4), dimension(:), intent(in) :: b
    integer(kind=4), intent(out) :: diff
    external :: diff_i
    call diff_i(n,a(1,1),b(1),diff)
  end subroutine f_diff_i2i1
  subroutine f_diff_i2(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    integer(kind=4), dimension(:,:),   intent(in) :: a
    integer(kind=4), dimension(:,:), intent(in) :: b
    integer(kind=4), intent(out) :: diff
    external :: diff_i
    call diff_i(n,a(1,1),b(1,1),diff)
  end subroutine f_diff_i2
  subroutine f_diff_i1(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    integer(kind=4), dimension(:),   intent(in) :: a
    integer(kind=4), dimension(:), intent(in) :: b
    integer(kind=4), intent(out) :: diff
    external :: diff_i
    call diff_i(n,a(1),b(1),diff)
  end subroutine f_diff_i1
  subroutine f_diff_i1i2(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    integer(kind=4), dimension(:),   intent(in) :: a
    integer(kind=4), dimension(:,:), intent(in) :: b
    integer(kind=4), intent(out) :: diff
    external :: diff_i
    call diff_i(n,a(1),b(1,1),diff)
  end subroutine f_diff_i1i2


  subroutine f_diff_li(n,a_add,b_add,diff)
    implicit none
    integer, intent(in) :: n
    integer(kind=8), intent(inout) :: a_add
    integer(kind=8), intent(inout) :: b_add
    integer(kind=8), intent(out) :: diff
    external :: diff_li
    call diff_li(n,a_add,b_add,diff)
  end subroutine f_diff_li
  subroutine f_diff_li2li1(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    integer(kind=8), dimension(:,:),   intent(in) :: a
    integer(kind=8), dimension(:), intent(in) :: b
    integer(kind=8), intent(out) :: diff
    external :: diff_li
    call diff_li(n,a(1,1),b(1),diff)
  end subroutine f_diff_li2li1
  subroutine f_diff_li2(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    integer(kind=8), dimension(:,:),   intent(in) :: a
    integer(kind=8), dimension(:,:), intent(in) :: b
    integer(kind=8), intent(out) :: diff
    external :: diff_li
    call diff_li(n,a(1,1),b(1,1),diff)
  end subroutine f_diff_li2
  subroutine f_diff_li1(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    integer(kind=8), dimension(:),   intent(in) :: a
    integer(kind=8), dimension(:), intent(in) :: b
    integer(kind=8), intent(out) :: diff
    external :: diff_li
    call diff_li(n,a(1),b(1),diff)
  end subroutine f_diff_li1
  subroutine f_diff_li1li2(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    integer(kind=8), dimension(:),   intent(in) :: a
    integer(kind=8), dimension(:,:), intent(in) :: b
    integer(kind=8), intent(out) :: diff
    external :: diff_li
    call diff_li(n,a(1),b(1,1),diff)
  end subroutine f_diff_li1li2


  subroutine f_diff_r(n,a_add,b_add,diff)
    implicit none
    integer, intent(in) :: n
    real, intent(inout) :: a_add
    real, intent(inout) :: b_add
    real, intent(out) :: diff
    external :: diff_r
    call diff_r(n,a_add,b_add,diff)
  end subroutine f_diff_r

  subroutine f_diff_d(n,a_add,b_add,diff)
    implicit none
    integer, intent(in) :: n
    double precision, intent(inout) :: a_add
    double precision, intent(inout) :: b_add
    double precision, intent(out) :: diff
    external :: diff_d
    call diff_d(n,a_add,b_add,diff)
  end subroutine f_diff_d
  subroutine f_diff_d1(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    double precision, dimension(:),   intent(in) :: a
    double precision, dimension(:), intent(in) :: b
    double precision, intent(out) :: diff
    external :: diff_d
    call diff_d(n,a(1),b(1),diff)
  end subroutine f_diff_d1
  subroutine f_diff_d2d3(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    double precision, dimension(:,:),   intent(in) :: a
    double precision, dimension(:,:,:), intent(in) :: b
    double precision, intent(out) :: diff
    external :: diff_d
    call diff_d(n,a(1,1),b(1,1,1),diff)
  end subroutine f_diff_d2d3
  subroutine f_diff_d2d1(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    double precision, dimension(:,:),   intent(in) :: a
    double precision, dimension(:), intent(in) :: b
    double precision, intent(out) :: diff
    external :: diff_d
    call diff_d(n,a(1,1),b(1),diff)
  end subroutine f_diff_d2d1
  subroutine f_diff_d0d1(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    double precision, intent(inout) :: a
    double precision, dimension(:), intent(in) :: b
    double precision, intent(out) :: diff
    external :: diff_d
    call diff_d(n,a,b(1),diff)
  end subroutine f_diff_d0d1

  subroutine f_diff_d2(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    double precision, dimension(:,:),   intent(in) :: a
    double precision, dimension(:,:), intent(in) :: b
    double precision, intent(out) :: diff
    external :: diff_d
    call diff_d(n,a(1,1),b(1,1),diff)
  end subroutine f_diff_d2
  subroutine f_diff_d1d2(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    double precision, dimension(:),   intent(in) :: a
    double precision, dimension(:,:), intent(in) :: b
    double precision, intent(out) :: diff
    external :: diff_d
    call diff_d(n,a(1),b(1,1),diff)
  end subroutine f_diff_d1d2



  subroutine f_diff_l(n,a_add,b_add,diff)
    implicit none
    integer, intent(in) :: n
    logical, intent(inout) :: a_add
    logical, intent(inout) :: b_add
    logical, intent(out) :: diff
    external :: diff_l
    call diff_l(n,a_add,b_add,diff)
  end subroutine f_diff_l

  subroutine f_diff_c1i1(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    character, dimension(:),   intent(in) :: a
    integer(f_integer), dimension(:), intent(in) :: b
    integer(f_integer), intent(out) :: diff
    external :: diff_ci
    call diff_ci(n,a(1),b(1),diff)
  end subroutine f_diff_c1i1

  subroutine f_diff_c1li1(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    character, dimension(:),   intent(in) :: a
    integer(f_long), dimension(:), intent(in) :: b
    integer(f_long), intent(out) :: diff
    external :: diff_ci
    call diff_ci(n,a(1),b(1),diff)
  end subroutine f_diff_c1li1

  subroutine f_diff_c0i1(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    character(len=*),   intent(in) :: a
    integer(f_integer), dimension(:), intent(in) :: b
    integer(f_integer), intent(out) :: diff
    external :: diff_ci
    call diff_ci(n,a,b(1),diff)
  end subroutine f_diff_c0i1

  subroutine f_diff_c0li1(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    character(len=*),   intent(in) :: a
    integer(f_long), dimension(:), intent(in) :: b
    integer(f_long), intent(out) :: diff
    external :: diff_ci
    call diff_ci(n,a,b(1),diff)
  end subroutine f_diff_c0li1

  subroutine f_diff_li0li1(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    integer(f_long), intent(inout) :: a
    integer(f_long), dimension(:), intent(in) :: b
    integer(f_long), intent(out) :: diff
    external :: diff_li
    call diff_li(n,a,b(1),diff)
  end subroutine f_diff_li0li1

  subroutine f_diff_i0i1(n,a,b,diff)
    implicit none
    integer, intent(in) :: n
    integer(f_integer), intent(inout) :: a
    integer(f_integer), dimension(:), intent(in) :: b
    integer(f_integer), intent(out) :: diff
    external :: diff_i
    call diff_i(n,a,b(1),diff)
  end subroutine f_diff_i0i1

  pure subroutine zero_string(str)
    use yaml_strings, only: f_strcpy
    implicit none
    character(len=*), intent(out) :: str
    call f_strcpy(src=' ',dest=str)
  end subroutine zero_string

  pure subroutine zero_li(val)
    implicit none
    integer(f_long), intent(out) :: val
    val=int(0,f_long)
  end subroutine zero_li

  pure subroutine zero_i(val)
    implicit none
    integer(f_integer), intent(out) :: val
    val=0
  end subroutine zero_i

  pure subroutine zero_r(val)
    implicit none
    real, intent(out) :: val
    val=0.e0
  end subroutine zero_r

  pure subroutine zero_d(val)
    implicit none
    double precision, intent(out) :: val
    val=0.d0
  end subroutine zero_d

  pure subroutine zero_l(val)
    implicit none
    logical, intent(out) :: val
    val=.false.
  end subroutine zero_l

  subroutine put_to_zero_simple(n,da)
    implicit none
    integer, intent(in) :: n
    real :: da

    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_simple(n,da)
    call setzero(int(n,f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_simple

  subroutine put_to_zero_double(n,da)
    implicit none
    integer, intent(in) :: n
    double precision, intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero(n,da)
    call setzero(int(n,f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_double

  subroutine put_to_zero_double_1(da)
    implicit none
    double precision, dimension(:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero(size(da),da(lbound(da,1)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_double_1

  subroutine put_to_zero_double_2(da)
    implicit none
    double precision, dimension(:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero(size(da),da(lbound(da,1),lbound(da,2)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_double_2

  subroutine put_to_zero_double_3(da)
    implicit none
    double precision, dimension(:,:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO) 
    !call razero(size(da),da(lbound(da,1),lbound(da,2),lbound(da,3)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume() 
  end subroutine put_to_zero_double_3

  subroutine put_to_zero_double_4(da)
    implicit none
    double precision, dimension(:,:,:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO) 
    !call razero(size(da),da(lbound(da,1),lbound(da,2),lbound(da,3),lbound(da,4)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume() 
  end subroutine put_to_zero_double_4

  subroutine put_to_zero_double_5(da)
    implicit none
    double precision, dimension(:,:,:,:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO) 
    !call razero(size(da),da(lbound(da,1),lbound(da,2),lbound(da,3),lbound(da,4),lbound(da,5)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume() 
  end subroutine put_to_zero_double_5

  subroutine put_to_zero_double_6(da)
    implicit none
    double precision, dimension(:,:,:,:,:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO) 
    !call razero(size(da),da(lbound(da,1),lbound(da,2),lbound(da,3),lbound(da,4),lbound(da,5),lbound(da,6)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume() 
  end subroutine put_to_zero_double_6

  subroutine put_to_zero_double_7(da)
    implicit none
    double precision, dimension(:,:,:,:,:,:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO) 
    !call razero(size(da),da(lbound(da,1),lbound(da,2),lbound(da,3),lbound(da,4),lbound(da,5),lbound(da,6),lbound(da,7)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume() 
  end subroutine put_to_zero_double_7

  subroutine put_to_zero_integer(n,da)
    implicit none
    integer, intent(in) :: n
    integer(f_integer) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_integer(n,da)
    call setzero(int(n,f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_integer

  subroutine put_to_zero_integer1(da)
    implicit none
    integer(f_integer), dimension(:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_integer(size(da),da(lbound(da,1)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_integer1

  subroutine put_to_zero_integer2(da)
    implicit none
    integer(f_integer), dimension(:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_integer(size(da),da(lbound(da,1),lbound(da,2)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_integer2

  subroutine put_to_zero_integer3(da)
    implicit none
    integer(f_integer), dimension(:,:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_integer(size(da),da(lbound(da,1),lbound(da,2),lbound(da,3)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_integer3


  subroutine put_to_zero_long(n,da)
    implicit none
    integer, intent(in) :: n
    integer(f_long) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_integerlong(n,da)
    call setzero(int(n,f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_long

  subroutine put_to_zero_long1(da)
    implicit none
    integer(f_long), dimension(:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_integerlong(size(da),da(lbound(da,1)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_long1

  subroutine put_to_zero_long2(da)
    implicit none
    integer(f_long), dimension(:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_integerlong(size(da),da(lbound(da,1),lbound(da,2)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_long2

  subroutine put_to_zero_long3(da)
    implicit none
    integer(f_long), dimension(:,:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_integerlong(size(da),da(lbound(da,1),lbound(da,2),lbound(da,3)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_long3
  
end module f_utils
