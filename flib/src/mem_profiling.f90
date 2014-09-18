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

  public !<low-level module

  !Memory profiling
  type, public :: memstat
     character(len=36) :: routine,array
     integer(kind=8) :: memory,peak
  end type memstat

  type, public :: memory_state
     integer :: memalloc !< number of allocations recorded
     integer :: memdealloc !<number of deallocations recorded
     type(memstat) :: memloc !<state of the memory in the local routine
     type(memstat) :: memtot !<global state of the memory in profiler instance
  end type memory_state

  real :: memorylimit = 0.e0 !< limit of the memory allowed, in Gb

contains

  subroutine f_set_memory_limit(limit)
    real, intent(in) :: limit
    memorylimit = limit
  end subroutine f_set_memory_limit

  !> dump variable has now become useless
  subroutine memstate_report(memstate,dump)
    use yaml_output
    implicit none
    type(memory_state), intent(in) :: memstate
    logical, intent(in), optional :: dump
    !local variables
    logical :: dmp
    dmp=.true.
    if (present(dump)) dmp=dump
    if (dmp) then
       call yaml_mapping_open('Memory Consumption Report')
       call yaml_map('Tot. No. of Allocations',memstate%memalloc)
       call yaml_map('Tot. No. of Deallocations',memstate%memdealloc)
       call yaml_map('Remaining Memory (B)',memstate%memtot%memory)
       call yaml_mapping_open('Memory occupation')
       call yaml_map('Peak Value (MB)',memstate%memtot%peak/int(1048576,kind=8))
       call yaml_map('for the array',trim(memstate%memtot%array))
       call yaml_map('in the routine',trim(memstate%memtot%routine))
       call yaml_mapping_close()
       call yaml_mapping_close()
    else
       !call memocc(0,-1,'count','stop')
    end if
  end subroutine memstate_report

  !> Put to zero memocc counters
  pure subroutine memstate_init(memstate)
    implicit none
    type(memory_state), intent(out) :: memstate

    memstate%memtot%memory=int(0,kind=8)
    memstate%memtot%peak=int(0,kind=8)
    memstate%memtot%routine=''
    memstate%memtot%array=''

    memstate%memalloc=0
    memstate%memdealloc=0

    memstate%memloc%routine='routine'
    memstate%memloc%array='array'
    memstate%memloc%memory=int(0,kind=8) !fake initialisation to print the first routine
    memstate%memloc%peak=int(0,kind=8)
  end subroutine memstate_init

  subroutine memstate_update(memstate,isize,array,routine)
    use yaml_output
    use dictionaries!error_handling
    implicit none
    ! Arguments
    type(memory_state), intent(inout) :: memstate
    integer(kind=8), intent(in) :: isize
    character(len=*), intent(in) :: array,routine
    ! Local variables
    !logical :: lmpinit
    !logical :: linq
    !integer :: ierr,istat_del,impinit
    character(len=256) :: message

    !Total counter, for all the processes
    memstate%memtot%memory=memstate%memtot%memory+isize
    if (memstate%memtot%memory > memstate%memtot%peak) then
       memstate%memtot%peak=memstate%memtot%memory
       memstate%memtot%routine=routine
       memstate%memtot%array=array
    end if
    if (isize > int(0,kind=8)) then
       memstate%memalloc=memstate%memalloc+1
    else if (isize < int(0,kind=8)) then
       memstate%memdealloc=memstate%memdealloc+1
    end if

    if (memorylimit /= 0.e0 .and. &
         memstate%memtot%memory > int(real(memorylimit,kind=8)*1073741824.d0,kind=8)) then !memory limit is in GB
       write(message,'(a,f7.3,a,i0,a)')&
            'Limit of ',memorylimit,' GB reached, total memory is ',memstate%memtot%memory,' B. '

       call f_err_throw(trim(message)//' Array '//trim(memstate%memtot%array)//&
            ', routine '//trim(memstate%memtot%routine),err_name='ERR_MEMLIMIT')
       return
    end if
    if (trim(memstate%memloc%routine) /= routine) then
       memstate%memloc%routine=routine
       memstate%memloc%array=array
       memstate%memloc%memory=isize
       memstate%memloc%peak=isize
    else
       memstate%memloc%memory=memstate%memloc%memory+isize
       if (memstate%memloc%memory > memstate%memloc%peak) then
          memstate%memloc%peak=memstate%memloc%memory
          memstate%memloc%array=array
       end if
    end if

  END SUBROUTINE memstate_update

end module memory_profiling
