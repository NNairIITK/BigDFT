!> @file
!! Memory profiling module (flib)
!! @author
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module which contains routines to profile the code.
!! Currently only memory occupation are provided here.
!! Ideally, timab should be incorporated here.
module memory_profiling
  implicit none

  private

  !Memory profiling
  type :: memstat
     character(len=36) :: routine,array
     integer(kind=8) :: memory,peak
  end type memstat

  real :: memorylimit = 0.e0
  logical :: meminit = .false.
  integer :: ndebug =0 !desactivated here
!  integer, parameter :: mallocFile = 98
!  character(len=256) :: filename=repeat(' ',256)
!  integer :: stdout=6
  type(memstat), save :: memloc,memtot
  integer, save :: memalloc,memdealloc!,memproc = 0
  !Debug option for memocc, set in the input file
  !logical :: memdebug=.true.
!  integer :: malloc_level=2

  !> Interface for the memory allocation control, depends on ndebug
  !! Control the memory occupation by calculating the overall size of the allocated arrays
  !!   when allocating allocating an array "stuff" of dimension n in the routine "dosome":
  !!      allocate(stuff(n),stat=i_stat)
  !!      call memocc(i_stat,stuff,'stuff','dosome')
  !!   when deallocating:
  !!      i_all=-product(shape(stuff))*kind(stuff)
  !!      deallocate(stuff,stat=i_stat)
  !!      call memocc(i_stat,i_all,'stuff','dosome')
  !!   The counters are initialized once by the first allocation
  !!   and stopped with:
  !!      call memocc(0,0,'count','stop')
  !!   At the end of the calculation a short report is printed on the screen,
  !!   some information can be also written on disk following the needs
  !!
  !!   The file malloc.prc is not deleted if the final total memory is not equal
  !!   to zero.
  !!   debug (parameter)
  !!     == .true.  verbose format (useful with utils/scripts/memcheck.py)
  !!                then display a line per allocation or deallocation
  !!                a routine at the end parses the file
  !!     == .false. compact format
  interface memocc
     module procedure memocc_internal  !< Central routine to be used for deallocation 
     !mo_dp1,mo_dp2,mo_dp3,mo_dp4,mo_dp5,mo_dp6,mo_dp7,&
     !     mo_sp1,mo_sp2,mo_sp3,mo_sp4,mo_sp5,mo_sp6,mo_sp7,&
     !     mo_i1,mo_i2,mo_i3,mo_i4,mo_i5,mo_i6,mo_i7,&
     !     mo_l1,mo_l2,mo_l3,mo_l4,mo_l5,mo_l6,mo_l7,&
     !     mo_c1, mo_cmpdp1, &
  end interface

  public :: memocc!,ndebug
!  public :: memocc_set_state
!  public :: memocc_set_stdout
!  public :: memocc_set_filename
  public :: memocc_set_memory_limit
  public :: memocc_report

contains

!  !> State of malloc.prc file and of counters
!  !! The status can only be downgraded. A stop signal is produced if status is increased
!  subroutine memocc_set_state(istatus)
!    !> 0 no file malloc.prc is created, only memory allocation counters running
!    !! 1 file malloc.prc is created in a light version (only current information is written) (this version is now eliminated, equal to 0 functionality)
!    !! 2 file malloc.prc is created with full information inside (default state if not specified)
!    use yaml_output, only: yaml_warning
!    integer, intent(in) :: istatus
!    !local variable
!    logical :: linq
!    integer :: istat_del,ierr
!    if (istatus > malloc_level) then
!       !here we should replace by yaml_warning
!       !write(7,*) 'WARNING: malloc_level can be only downgraded, ignoring'
!       call yaml_warning('malloc_level can be only downgraded, ignoring')
!       return
!    end if
!
!    malloc_level = istatus
!
!    if (istatus == 2) return !the default situation
!
!    !inquire for unit opened
!    linq=.false.
!    inquire(unit=mallocFile,opened=linq,iostat=ierr)
!    linq= linq .and. ierr==0
!
!!!$    if (istatus == 1 .and. memproc==0) then 
!!!$       !clean the file situation (delete the previously existing file)
!!!$       if (linq) close(unit=mallocFile)                        
!!!$       !call delete(trim(filename),len(trim(filename)),istat_del)
!!!$       if (len_trim(filename) == 0) then
!!!$          open(unit=mallocFile,file="malloc.prc",status='unknown',action='write')
!!!$       else
!!!$          open(unit=mallocFile,file=trim(filename),status='unknown',action='write')
!!!$       end if
!!!$    end if
!
!    if (istatus <= 1 .and. memproc==0) then
!       !the file should be deleted
!       if (linq) close(unit=mallocFile)
!       !open(unit=mallocFile,file='malloc.prc',status='replace')
!       !close(unit=mallocFile)
!       linq=.false.
!       if (len_trim(filename) /= 0) then
!          inquire(file=trim(filename),exist=linq,iostat=ierr)
!          linq = linq .and. ierr ==0
!          if (linq) call delete(trim(filename),len(trim(filename)),istat_del)
!       end if
!    end if
!  end subroutine memocc_set_state
!
  subroutine memocc_set_memory_limit(limit)
    real, intent(in) :: limit

    memorylimit = limit
!    filename=repeat(' ',len(filename))
!    filename='malloc.prc'
  end subroutine memocc_set_memory_limit

!!$  subroutine memocc_set_filename(file)
!!$    character(len=*), intent(in) :: file
!!$    !local variables
!!$    integer :: ipos
!!$
!!$    ipos=min(len(trim(file)),len(filename))
!!$    filename=repeat(' ',len(filename))
!!$    filename(1:ipos)=file(1:ipos)
!!$
!!$  end subroutine memocc_set_filename

  subroutine memocc_report(dump)
    implicit none
    logical, intent(in), optional :: dump
    !local variables
    logical :: dmp
    dmp=.true.
    if (present(dump)) dmp=dump
    if (dmp) then
       call memocc(0, 0,'count','stop')
    else
       call memocc(0,-1,'count','stop')
    end if
  end subroutine memocc_report

  !> Put to zero memocc counters
  subroutine memocc_variables_init()
    memtot%memory=int(0,kind=8)
    memtot%peak=int(0,kind=8)
    memtot%routine=''
    memtot%array=''
    memalloc=0
    memdealloc=0

    memloc%routine='routine'
    memloc%array='array'
    memloc%memory=int(0,kind=8) !fake initialisation to print the first routine
    memloc%peak=int(0,kind=8)
  end subroutine memocc_variables_init

!!$  subroutine memocc_open_file()
!!$    open(unit=mallocFile,file=trim(filename),status='unknown')
!!$    !if (memdebug) then
!!$    write(mallocFile ,'(a,t40,a,t70,4(1x,a12))')&
!!$         '(Data in KB) Routine','Array name    ',&
!!$         'Array size','Total Memory'
!!$    !else
!!$    !write(mallocFile,'(a,t40,a,t70,4(1x,a12))')&
!!$    !     '(Data in KB) Routine','Peak Array    ',&
!!$    !     'Routine Mem','Routine Peak','Memory Stat.','Memory Peak'
!!$    !end if
!!$  end subroutine memocc_open_file

!!$  subroutine memocc_set_stdout(unit)
!!$    implicit none
!!$    integer, intent(in) :: unit
!!$
!!$    stdout=unit
!!$
!!$  end subroutine memocc_set_stdout


  subroutine memocc_internal(istat,isize,array,routine)
    use yaml_output
    use dictionaries!error_handling
    implicit none

    ! Arguments
    integer, intent(in) :: istat,isize
    character(len=*), intent(in) :: array,routine

    ! Local variables
    !logical :: lmpinit
    logical :: linq
    integer :: ierr,istat_del!,impinit
    character(len=256) :: message

    !print *,memproc,array,routine
    ! Initialised first
    if (.not.meminit) then
       !the mpi routines have to be eliminated
       call memocc_variables_init()
       meminit = .true.
    end if

    select case(array)
    case('count')
       if (trim(routine)=='stop') then
          if (isize /= -1) then
             call yaml_mapping_open('Memory Consumption Report')
             call yaml_map('Tot. No. of Allocations',memalloc)
             call yaml_map('Tot. No. of Deallocations',memdealloc)
             call yaml_map('Remaining Memory (B)',memtot%memory)
             call yaml_mapping_open('Memory occupation')
             call yaml_map('Peak Value (MB)',memtot%peak/int(1048576,kind=8))
             call yaml_map('for the array',trim(memtot%array))
             call yaml_map('in the routine',trim(memtot%routine))
             call yaml_mapping_close()
             call yaml_mapping_close()
          end if

       end if

    case default
       !Total counter, for all the processes
       memtot%memory=memtot%memory+int(isize,kind=8)
       if (memtot%memory > memtot%peak) then
          memtot%peak=memtot%memory
          memtot%routine=routine
          memtot%array=array
       end if
       if (isize>0) then
          memalloc=memalloc+1
       else if (isize<0) then
          memdealloc=memdealloc+1
       end if
       
       if (memorylimit /= 0.e0 .and. &
            memtot%memory > int(real(memorylimit,kind=8)*1073741824.d0,kind=8)) then !memory limit is in GB
          write(message,'(a,f7.3,a,i0,a)')&
               'Limit of ',memorylimit,' GB reached, total memory is ',memtot%memory,' B. '

          call f_err_throw(trim(message)//' Array '//trim(memtot%array)//&
               ', routine '//trim(memtot%routine),err_name='ERR_MEMLIMIT')
          return
       end if
       if (trim(memloc%routine) /= routine) then
          memloc%routine=routine
          memloc%array=array
          memloc%memory=isize
          memloc%peak=isize
       else
          memloc%memory=memloc%memory+isize
          if (memloc%memory > memloc%peak) then
             memloc%peak=memloc%memory
             memloc%array=array
          end if
       end if
    end select
  END SUBROUTINE memocc_internal


  !> Check the malloc.prc file (verbose format)
!!$  subroutine memory_malloc_check(nalloc,ndealloc)
!!$    implicit none
!!$    !Arguments
!!$    integer, intent(in) :: nalloc,ndealloc
!!$    !Local variables
!!$    if (malloc_level==2 .and. nalloc /= ndealloc) then
!!$       !Use # to be yaml compliant (is a comment in yaml)
!!$       write(*,*) &
!!$            "#Use the python script 'memcheck.py' in utils/scripts to check "//&
!!$            trim(filename)//" file"
!!$    end if
!!$  END SUBROUTINE memory_malloc_check

end module memory_profiling
