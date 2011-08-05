module module_input

  use module_base

  implicit none
  private

  logical :: output
  character(len = 100) :: input_file
  character(len = 100), allocatable :: input_lines(:)

  interface input_var
     module procedure var_character, var_logical, var_integer, &
          & var_integer_array, var_double, var_keyword, var_ids
  end interface

  public :: input_set_file
  public :: input_free
  public :: input_var

contains

  subroutine input_set_file(iproc, filename, exists)
    integer, intent(in) :: iproc
    character(len = *), intent(in) :: filename
    logical, intent(out) :: exists

    character(len = 100) :: lines(500)
    integer :: i, ierror

    inquire(file=trim(filename),exist=exists)
    if (exists) then
       open(unit = 1, file = trim(filename), status = 'old')
       i = 1
       do 
          read(1, fmt = '(a)', iostat = ierror) lines(i)
          if (ierror /= 0) exit
          i = i + 1
       end do
       close(1)

!!$    write(0,*) "Setup input file '", trim(filename), "' with ", i - 1, "lines."
       write(input_file, "(A)") trim(filename)
       allocate(input_lines(i - 1))
       do i = 1, size(input_lines), 1
          input_lines(i) = lines(i)
       end do
    end if
    output = (iproc == 0)
    ! Output
    if (iproc == 0) then
       write(*,*)
       if (exists) then
          write(*,'(1x,3a)') '--- (file: ', trim(filename), &
               & ') ----------------------------------------- Performance Options'
       else
          write(*,'(1x,a)')&
               '--- (file: input.perf -- not present) -------------------------- Performance Options'
       end if
    end if
  end subroutine input_set_file

  subroutine input_free()
    if (allocated(input_lines)) deallocate(input_lines)
  end subroutine input_free

  subroutine find(name, iline, ii)
    character(len = *), intent(in) :: name
    integer, intent(out) :: iline, ii
    
    integer :: k

    if (allocated(input_lines)) then
       do iline = 1, size(input_lines), 1
          k = 1
          do ii = 1, len(input_lines(iline)), 1
             if (ichar(input_lines(iline)(ii:ii)) == ichar(name(k:k)) .or. &
                  & ichar(input_lines(iline)(ii:ii)) == ichar(name(k:k)) + 32 .or. &
                  & ichar(input_lines(iline)(ii:ii)) == ichar(name(k:k)) - 32) then
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
  end subroutine find

  subroutine var_character(name, default, description, var)
    character(len = *), intent(in) :: name
    character(len = *), intent(in) :: default
    character(len = *), intent(in) :: description
    character(len = *), intent(out) :: var

    integer :: i, j, ierror, ierr

    write(var, "(A)") default
    call find(name, i, j)
    if (i > 0) then
       read(input_lines(i)(j + 2:), fmt = *, iostat = ierror) var
       if (ierror/=0) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
    end if
    if (output) write(*,"(1x,a,3x,a,1x,a,t30,2a)") "|", name, var, '!', description
  end subroutine var_character

  subroutine var_logical(name, default, description, var)
    character(len = *), intent(in) :: name
    logical, intent(in) :: default
    character(len = *), intent(in) :: description
    logical, intent(out) :: var

    integer :: i, j, ierror

    var = default
    call find(name, i, j)
    if (i > 0) then
       read(input_lines(i)(j + 2:), fmt = *, iostat = ierror) var
       if (ierror /= 0) then
          var = .true.
       end if
    end if
    if (output) write(*,"(1x,a,3x,a,1x,l1,t30,2a)") "|", name, var, '!', description
  end subroutine var_logical

  subroutine var_integer(name, default, description, var)
    character(len = *), intent(in) :: name
    integer, intent(in) :: default
    character(len = *), intent(in) :: description
    integer, intent(out) :: var

    integer :: i, j, ierror, ierr

    var = default
    call find(name, i, j)
    if (i > 0) then
       read(input_lines(i)(j + 2:), fmt = *, iostat = ierror) var
       if (ierror/=0) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
    end if
    if (output) write(*,"(1x,a,3x,a,1x,I0,t30,2a)") "|", name, var, '!', description
  end subroutine var_integer

  subroutine var_integer_array(name, default, description, var)
    character(len = *), intent(in) :: name
    integer, intent(in) :: default(:)
    character(len = *), intent(in) :: description
    integer, intent(out) :: var(:)

    integer :: i, j, ierror, ierr

    var = default
    call find(name, i, j)
    if (i > 0) then
       read(input_lines(i)(j + 2:), fmt = *, iostat = ierror) var
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
  end subroutine var_integer_array

  subroutine var_double(name, default, description, var)
    character(len = *), intent(in) :: name
    double precision, intent(in) :: default
    character(len = *), intent(in) :: description
    double precision, intent(out) :: var

    integer :: i, j, ierror, ierr

    var = default
    call find(name, i, j)
    if (i > 0) then
       read(input_lines(i)(j + 2:), fmt = *, iostat = ierror) var
       if (ierror/=0) then
          if (output) write(*,'(1x,a,a,a,i3)')  'Error while reading the file "', &
               & trim(input_file), '", line=', i
          call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
       end if
    end if
    if (output) write(*,"(1x,a,3x,a,1x,es9.2,t30,2a)") "|", name, var, '!', description
  end subroutine var_double

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
       read(input_lines(i)(j + 2:), fmt = *, iostat = ierror) buf
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
  end subroutine var_keyword

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
       read(input_lines(i)(j + 2:), fmt = *, iostat = ierror) var
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
  end subroutine var_ids

end module module_input
