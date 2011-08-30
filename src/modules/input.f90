module module_input

  use module_base

  implicit none
  private

  integer, parameter :: nmax_lines=500,max_length=100
  character(len=max_length) :: input_file,line_being_processed
  logical :: output
  integer :: iline_parsed,iline_written,iargument,ipos
  character(len=max_length), dimension(:), allocatable :: inout_lines


  interface input_var
     module procedure var_character, var_logical, var_integer, &
          & var_integer_array, var_double, var_keyword, var_ids,&
          & var_double_compulsory
  end interface

  public :: input_set_file
  public :: input_free
  public :: input_var

contains

  subroutine input_set_file(iproc, filename, exists,comment_file_usage)
    integer, intent(in) :: iproc
    character(len = *), intent(in) :: filename,comment_file_usage
    logical, intent(out) :: exists

    character(len=max_length), dimension(nmax_lines) :: lines
    integer :: i,nlines,ierror,ierr=0 !for MPIfake BCAST

    !no line parsed if the file not exists
    iline_parsed=0
    !no line has been written in the output
    iline_written=0
    !no argument has been read yet
    iargument=0
    !the present line is empty
    line_being_processed=repeat(' ',max_length)

    !check if the file is present
    inquire(file=trim(filename),exist=exists)
    if (exists) then
       !only the root processor parse the file
       if (iproc==0) then
          open(unit = 1, file = trim(filename), status = 'old')
          i = 1
          parse_file: do 
             lines(i)=repeat(' ',max_length) !initialize lines
             read(1, fmt = '(a)', iostat = ierror) lines(i)
             !eliminate leading blanks from the line
             lines(i)=adjustl(lines(i))
             if (ierror /= 0) exit parse_file
             i = i + 1
          end do parse_file
          close(1)
          nlines=i-1
       end if
       !broadcast the number of lines
       call MPI_BCAST(nlines,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       if (ierr /=0) stop 'input_file BCAST (1) '
       !broadcast all the lines
       call MPI_BCAST(lines,nmax_lines*nlines,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
       if (ierr /=0) stop 'input_file BCAST (2) '

!!$    write(0,*) "Setup input file '", trim(filename), "' with ", i - 1, "lines."
       write(input_file, "(A)") trim(filename)

       allocate(inout_lines(0:nlines)) !the 0-th line is for the description of the file
       do i=1,nlines
          inout_lines(i)=lines(i)
       end do
       !start parsing from the first line
       iline_parsed=1
    else
       !in this case the array constitute the output results
       allocate(inout_lines(0:nmax_lines))
       do i=1,nmax_lines
          inout_lines(i)=repeat(' ',max_length) !initialize lines
       end do
    end if

    !write the first line in the output
    if (exists) then
       write(inout_lines(iline_written),'(1x,3a)') '--- (file: ', trim(filename), &
            & ') -----------------------------------------'//&
            trim(comment_file_usage)
    else
       write(inout_lines(iline_written),'(1x,a)')&
            '--- (file:'//trim(filename)//'-- not present) --------------------------'//&
            trim(comment_file_usage)
    end if
    iline_written=iline_written+1

!!$    output = (iproc == 0)
!!$    ! Output
!!$    if (iproc == 0) then
!!$       write(*,*)
!!$       if (exists) then
!!$          write(*,'(1x,3a)') '--- (file: ', trim(filename), &
!!$               & ') -----------------------------------------'//&
!!$               trim(comment_file_usage)
!!$       else
!!$          write(*,'(1x,a)')&
!!$               '--- (file:'//trim(filename)//'-- not present) --------------------------'//&
!!$               trim(comment_file_usage)
!!$       end if
!!$    end if
  END SUBROUTINE input_set_file

  subroutine input_free()
    implicit none
    !local variables
    integer :: iline
    
    !add the writing of the file in the given unit
    do iline=0,iline_written
       print *,inout_lines(iline)
    end do

    if (allocated(inout_lines)) deallocate(inout_lines)
    
  END SUBROUTINE input_free

  subroutine find(name, iline, ii)
    character(len = *), intent(in) :: name
    integer, intent(out) :: iline, ii
    
    integer :: k
    !change the allocate condition, since input lines is always used now
    !if (allocated(inout_lines)) then
    if (iline_parsed /= 0) then
       do iline = 1, size(inout_lines), 1
          k = 1
          do ii = 1, len(inout_lines(iline)), 1
             if (ichar(inout_lines(iline)(ii:ii)) == ichar(name(k:k)) .or. &
                  & ichar(inout_lines(iline)(ii:ii)) == ichar(name(k:k)) + 32 .or. &
                  & ichar(inout_lines(iline)(ii:ii)) == ichar(name(k:k)) - 32) then
                k = k + 1
             else
                k = 1
             end if
             if (k == len(name) + 1) then
                return
             end if
          end do
       end do
    end if
    iline = 0
  END SUBROUTINE find

  subroutine check(ierror)
    implicit none
    integer, intent(in) :: ierror
    !local variables
    integer :: ierr

    !increment the argument at each check
    iargument=iargument+1
    if (ierror/=0) then
       write(*,'(1x,a,a,2(a,i3))')'Error while reading the file "', &
            & trim(input_file), '", line=', iline_written,' argument=', iargument
       write(*,*)inout_lines(iline_written)
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       stop
    end if

  END SUBROUTINE check

  !> Process the line needed with the default value in mind
  subroutine process_line(default,line_comment)
    implicit none
    character(len=*), intent(in), optional :: default,line_comment
    !local variables
    integer :: ierror,i,iblank

    if (present(default) .and. .not. present(line_comment)) then
       if (iargument==0) ipos=0
       !case without file, write default and continue the line
       do i=1,len_trim(default)
          inout_lines(iline_written)(i+ipos:i+ipos)=default(i:i)
       end do
       ipos=len_trim(default)+1
       inout_lines(iline_written)(ipos:ipos)=' '
    else if (present(default) .and. present(line_comment)) then
       !case without file, close the line 
       write(inout_lines(iline_written)(ipos:),'(a)')line_comment
       iline_written=iline_written+1
       iargument=0
    else if (.not. present(default) .and. .not. present(line_comment)) then
       !traditional case, the argument should be parsed one after another
       !start with the entire line
       if (iargument==0) then 
          line_being_processed=inout_lines(iline_parsed)
       else
          !search in the line the first blank
          iblank=scan(line_being_processed,' ')
          do i=1,max_length-iblank
             line_being_processed(i:i)=line_being_processed(i+iblank:i+iblank)
          end do
          do i=max_length-iblank+1,max_length
             line_being_processed(i:i)=' '
          end do
       end if
       !adjust the line to eliminate further blanks
       line_being_processed=adjustl(line_being_processed)
    else if (.not. present(default) .and. present(line_comment)) then
       !traditional case, close the line and skip to the next one
       iargument=0
       iline_parsed=iline_parsed+1
       iline_written=iline_written+1
    end if
       
  end subroutine process_line

  subroutine var_double_compulsory(var,default,ranges,exclusive,comment)
    implicit none
    character(len=*), intent(in) :: default
    real(kind=8), intent(out) :: var
    character(len=*), intent(in), optional :: comment
    real(kind=8), dimension(2), intent(in), optional :: ranges
    real(kind=8), dimension(:), intent(in), optional :: exclusive
    !local variables
    logical :: found
    integer :: ierror,ilist,ierr
    !if the file has not been opened, use the default variable 
    !then write in the output lines the default
    if (iline_parsed==0) then
       read(default,*,iostat=ierror)var
       call check(ierror)
       call process_line(default=default)
       !finalize the line if the comment is present
       if (present(comment)) then
          call process_line(default=default,line_comment=comment)
       end if
    !otherwise read the corresponding argument and check its validity
    else
       !read the argument
       call process_line()

       read(line_being_processed,fmt=*,iostat=ierror) var
       call check(ierror)

       !check the validity of the variable
       if (present(ranges)) then
          if (var < ranges(1) .or. var > ranges(2)) then
             write(*,*)' ERROR in parsing file'//input_file//'line=', iline_written,' argument=', iargument
             write(*,*)'      values should be from',ranges(1),' to ',ranges(2)
             call MPI_BARRIER(MPI_COMM_WORLD,ierr)
             stop
          end if
       else if (present(exclusive)) then
          found=.false.
          found_loop: do ilist=1,size(exclusive)
             if (var == exclusive(ilist)) then
                found=.true.
                exit found_loop
             end if
          end do found_loop
          if (.not. found) then
             write(*,*)' ERROR in parsing file'//input_file//'line=', iline_written,' argument=', iargument
             write(*,*)'      values should be in list: ',exclusive(:)
             call MPI_BARRIER(MPI_COMM_WORLD,ierr)
             stop
          end if
       end if 
      
       !increment the line if comment is present, do not touch the input file
       if (present(comment)) then
          call process_line(line_comment=comment)
       end if
    end if   
  END SUBROUTINE var_double_compulsory

  subroutine var_character(name, default, description, var)
    character(len = *), intent(in) :: name
    character(len = *), intent(in) :: default
    character(len = *), intent(in) :: description
    character(len = *), intent(out) :: var

    integer :: i, j, ierror, ierr

    write(var, "(A)") default
    call find(name, i, j)
    if (i > 0) then
       read(inout_lines(i)(j + 2:), fmt = *, iostat = ierror) var
       if (ierror/=0) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
    end if
    if (output) write(*,"(1x,a,3x,a,1x,a,t30,2a)") "|", name, var, '!', description
  END SUBROUTINE var_character

  subroutine var_logical(name, default, description, var)
    character(len = *), intent(in) :: name
    logical, intent(in) :: default
    character(len = *), intent(in) :: description
    logical, intent(out) :: var

    integer :: i, j, ierror

    var = default
    call find(name, i, j)
    if (i > 0) then
       read(inout_lines(i)(j + 2:), fmt = *, iostat = ierror) var
       if (ierror /= 0) then
          var = .true.
       end if
    end if
    if (output) write(*,"(1x,a,3x,a,1x,l1,t30,2a)") "|", name, var, '!', description
  END SUBROUTINE var_logical

  subroutine var_integer(name, default, description, var)
    character(len = *), intent(in) :: name
    integer, intent(in) :: default
    character(len = *), intent(in) :: description
    integer, intent(out) :: var

    integer :: i, j, ierror, ierr

    var = default
    call find(name, i, j)
    if (i > 0) then
       read(inout_lines(i)(j + 2:), fmt = *, iostat = ierror) var
       if (ierror/=0) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
    end if
    if (output) write(*,"(1x,a,3x,a,1x,I0,t30,2a)") "|", name, var, '!', description
  END SUBROUTINE var_integer

  subroutine var_integer_array(name, default, description, var)
    character(len = *), intent(in) :: name
    integer, intent(in) :: default(:)
    character(len = *), intent(in) :: description
    integer, intent(out) :: var(:)

    integer :: i, j, ierror, ierr

    var = default
    call find(name, i, j)
    if (i > 0) then
       read(inout_lines(i)(j + 2:), fmt = *, iostat = ierror) var
       if (ierror/=0) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
    end if
    if (output) then
       write(*,"(1x,a,3x,a,1x)", advance = "NO") "|", name
       do i = 1, size(var), 1
          write(*,"(1x,I0)", advance = "NO") var(i)
       end do
       write(*,"(t7,2a)") '!', description
    end if
  END SUBROUTINE var_integer_array

  subroutine var_double(name, default, description, var)
    character(len = *), intent(in) :: name
    double precision, intent(in) :: default
    character(len = *), intent(in) :: description
    double precision, intent(out) :: var

    integer :: i, j, ierror, ierr

    var = default
    call find(name, i, j)
    if (i > 0) then
       read(inout_lines(i)(j + 2:), fmt = *, iostat = ierror) var
       if (ierror/=0) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
    end if
    if (output) write(*,"(1x,a,3x,a,1x,es9.2,t30,2a)") "|", name, var, '!', description
  END SUBROUTINE var_double

  subroutine var_keyword(name, length, default, list, description, var)
    character(len = *), intent(in) :: name
    integer, intent(in) :: length
    character(len = length), intent(in) :: default
    character(len = length), intent(in) :: list(:)
    character(len = *), intent(in) :: description
    integer, intent(out) :: var

    integer :: i, j, ierror, ierr
    character(len = length) :: buf

    ! Set the default value to var.
    do i = 1, size(list), 1
       if (trim(default) == trim(list(i))) exit
    end do
    var = i - 1
    ! Find the keyword name in the file.
    call find(name, i, j)
    if (i > 0) then
       read(inout_lines(i)(j + 2:), fmt = *, iostat = ierror) buf
       if (ierror/=0) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
       ! Look for buf in list.
       do j = 1, size(list), 1
          if (trim(buf) == trim(list(j))) exit
       end do
       if (j > size(list)) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
       var = j - 1
    end if
    if (output) then
       write(*,"(1x,a,3x,a,1x,a,t30,3a)", advance = "NO") &
            & "|", name, list(var + 1), '!', description, " ("
       write(*,"(A)", advance = "NO") trim(list(1))
       do i = 2, size(list), 1
          write(*,"(2A)", advance = "NO") ", ", trim(list(i))
       end do
       write(*,"(A)") ")"
    end if
  END SUBROUTINE var_keyword

  subroutine var_ids(name, default, list, description, var)
    character(len = *), intent(in) :: name
    integer, intent(in) :: default
    integer, intent(in) :: list(:)
    character(len = *), intent(in) :: description
    integer, intent(out) :: var

    integer :: i, j, ierror, ierr
    
    var = default
    call find(name, i, j)
    if (i > 0) then
       read(inout_lines(i)(j + 2:), fmt = *, iostat = ierror) var
       if (ierror/=0) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
       do j = 1, size(list), 1
          if (var == list(j)) exit
       end do
       if (j > size(list)) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
    end if
    if (output) write(*,"(1x,a,3x,a,1x,I0,t30,2a)") "|", name, var, '!', description
  END SUBROUTINE var_ids

end module module_input
