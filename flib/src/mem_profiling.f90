!!****m* ABINIT/m_profiling
!! NAME
!! m_profiling
!!
!! FUNCTION
!! This module contains routines to profile the code.
!! Currently only memory occupation are provided here.
!! Ideally, timab should be incorporated here.
!!
!! COPYRIGHT
!! Copyright (C) 2010 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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
  integer, parameter :: mallocFile = 98
  character(len=256) :: filename=repeat(' ',256)
  integer :: stdout=6
  type(memstat), save :: memloc,memtot
  integer :: memalloc,memdealloc,memproc = 0
  !Debug option for memocc, set in the input file
  !logical :: memdebug=.true.
  integer :: malloc_level=2

  !interface for the memory allocation control, depends on ndebug
  interface memocc
     module procedure mo_dp1,mo_dp2,mo_dp3,mo_dp4,mo_dp5,mo_dp6,mo_dp7,&
          mo_sp1,mo_sp2,mo_sp3,mo_sp4,mo_sp5,mo_sp6,mo_sp7,&
          mo_i1,mo_i2,mo_i3,mo_i4,mo_i5,mo_i6,mo_i7,&
          mo_l1,mo_l2,mo_l3,mo_l4,mo_l5,mo_l6,mo_l7,&
          mo_c1, mo_cmpdp1, &
          memocc_internal  !central routine to be used for deallocation
  end interface

  public :: memocc,ndebug
  public :: memocc_set_state
  public :: memocc_set_stdout
  public :: memocc_set_filename
  public :: memocc_set_memory_limit
  public :: memocc_report

contains

  !> State of malloc.prc file and of counters
  !!  @param istatus 0 no file malloc.prc is created, only memory allocation counters running
  !!                 1 file malloc.prc is created in a light version (only current information is written)
  !!                 2 file malloc.prc is created with full information inside (default state if not specified)
  !! The status can only be downgraded. A stop signal is produced if status is increased
  subroutine memocc_set_state(istatus)
    integer, intent(in) :: istatus
    !local variable
    integer :: istat_del
    if (istatus > malloc_level) then
       !here we should replace by yaml_warning
       !write(7,*) 'WARNING: malloc_level can be only downgraded, ignoring'
       return
    end if

    malloc_level = istatus

    if (istatus == 2) return !the default situation

    if (istatus == 1 .and. memproc==0) then 
       !clean the file situation (delete the previously existing file)
       close(unit=mallocFile)                        
       !call delete(trim(filename),len(trim(filename)),istat_del)
       if (trim(filename) == "") then
          open(unit=mallocFile,file="malloc.prc",status='unknown',action='write')
       else
          open(unit=mallocFile,file=trim(filename),status='unknown',action='write')
       end if
    end if

    if (istatus == 0 .and. memproc==0) then
       !the file should be deleted
       close(unit=mallocFile)
       !open(unit=mallocFile,file='malloc.prc',status='replace')
       !close(unit=mallocFile)
       call delete(trim(filename),len(trim(filename)),istat_del)
    end if
  end subroutine memocc_set_state

  subroutine memocc_set_memory_limit(limit)
    real, intent(in) :: limit

    memorylimit = limit
    filename=repeat(' ',len(filename))
    filename='malloc.prc'
  end subroutine memocc_set_memory_limit

  subroutine memocc_set_filename(file)
    character(len=*), intent(in) :: file
    !local variables
    integer :: ipos

    ipos=min(len(trim(file)),len(filename))
    filename=repeat(' ',len(filename))
    filename(1:ipos)=file(1:ipos)

  end subroutine memocc_set_filename

  subroutine memocc_report()
    call memocc(0,0,'count', 'stop')
  end subroutine memocc_report

  !> Put to zero memocc counters
  subroutine memocc_variables_init()
    memtot%memory=int(0,kind=8)
    memtot%peak=int(0,kind=8)
    memalloc=0
    memdealloc=0

    memloc%routine='routine'
    memloc%array='array'
    memloc%memory=int(0,kind=8) !fake initialisation to print the first routine
    memloc%peak=int(0,kind=8)
  end subroutine memocc_variables_init

  subroutine memocc_open_file()
    open(unit=mallocFile,file=trim(filename),status='unknown')
    !if (memdebug) then
    write(mallocFile,'(a,t40,a,t70,4(1x,a12))')&
         '(Data in KB) Routine','Array name    ',&
         'Array size','Total Memory'
    !else
    !write(mallocFile,'(a,t40,a,t70,4(1x,a12))')&
    !     '(Data in KB) Routine','Peak Array    ',&
    !     'Routine Mem','Routine Peak','Memory Stat.','Memory Peak'
    !end if
  end subroutine memocc_open_file

  subroutine memocc_set_stdout(unit)
    implicit none
    integer, intent(in) :: unit

    stdout=unit

  end subroutine memocc_set_stdout


  !!****f* ABINIT/memory_occupation
  !! FUNCTION
  !! Control the memory occupation by calculating the overall size of the allocated arrays
  !! DESCRIPTION
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
  !!   memdebug (parameter)
  !!     == .true.  verbose format (useful with utils/scripts/memcheck.py)
  !!                then display a line per allocation or deallocation
  !!                a routine at the end parses the file
  !!     == .false. compact format
  !! COPYRIGHT
  !!    Copyright (C) 2007-2010, BigDFT group (Luigi Genovese)
  !!    This file is distributed under the terms of the
  !!    GNU General Public License, see ~/COPYING file
  !!    or http://www.gnu.org/copyleft/gpl.txt .
  !!    For the list of contributors, see ~/AUTHORS 
  !! SOURCE
  !!
  subroutine memocc_internal(istat,isize,array,routine)
    use yaml_output
    use dictionaries!error_handling
    implicit none

    ! Arguments
    integer, intent(in) :: istat,isize
    character(len=*), intent(in) :: array,routine

    ! Local variables
    logical :: lmpinit
    integer :: ierr,istat_del
    character(len=256) :: message

    include 'mpif.h'

    !print *,memproc,array,routine
    ! Initialised first
    if (.not.meminit) then
       !the mpi routines have to be eliminated
       call memocc_variables_init()
       !Use MPI to have the mpi rank
       call MPI_INITIALIZED(lmpinit,ierr)
       if (lmpinit) then
          call MPI_COMM_RANK(MPI_COMM_WORLD,memproc,ierr)
       else
          !no-mpi case 
          memproc=0
       end if

       !open the writing file for the root process
       if (memproc == 0 .and. malloc_level > 0) then
          if (len(trim(filename))==0) then
             filename='malloc.prc'
          end if
          call memocc_open_file()
       end if
       meminit = .true.
    end if

    select case(array)
    case('count')
       if (trim(routine)=='stop' .and. memproc==0) then
          if (malloc_level > 0) then
             if (malloc_level == 1) rewind(mallocFile)
             write(mallocFile,'(a,t40,a,t70,4(1x,i12))')&
                  trim(memloc%routine),trim(memloc%array),&
                  memloc%memory/int(1024,kind=8),memloc%peak/int(1024,kind=8),&
                  memtot%memory/int(1024,kind=8),&
                  (memtot%peak+memloc%peak-memloc%memory)/int(1024,kind=8)
             close(unit=mallocFile)
          end if
          call yaml_open_map('Memory Consumption Report')
          call yaml_map('Tot. No. of Allocations',memalloc)
          call yaml_map('Tot. No. of Deallocations',memdealloc)
          call yaml_map('Remaining Memory (B)',memtot%memory)
          call yaml_open_map('Memory occupation')
          call yaml_map('Peak Value (MB)',memtot%peak/int(1048576,kind=8))
          call yaml_map('for the array',trim(memtot%array))
          call yaml_map('in the routine',trim(memtot%routine))
          call yaml_close_map()
          call yaml_close_map()

          !here we can add a routine which open the malloc.prc file in case of some 
          !memory allocation problem, and which eliminates it for a successful run
          if (malloc_level == 1 .and. memalloc == memdealloc .and. memtot%memory==int(0,kind=8)) then
             !remove file 
             call delete(trim(filename),len(trim(filename)),istat_del)
          else
             call memory_malloc_check(memalloc,memdealloc)
          end if
          ! no need to check since the calls are wrapped
!!$         else if (trim(routine)/='stop') then
!!$            write(*,*) "memocc: ",array," ",routine
!!$            write(*,"(a,i0,a)") "Error[",memproc,"]: Use memocc and the word 'count' only with the word 'stop'."
!!$            stop
       end if

    case default
       !control of the allocation/deallocation status (to be removed once f_malloc has been inserted)
       if (istat/=0) then
          if (memproc == 0 .and. malloc_level > 0) close(unit=mallocFile)
          write(message,'(1x,a)')'subroutine '//trim(routine)//', array '//trim(array)//&
               ', error code '//trim(yaml_toa(istat))
          if (f_err_raise(isize>=0,trim(message),err_name='ERR_ALLOCATE')) return
          if (f_err_raise(isize< 0,trim(message),err_name='ERR_DEALLOCATE')) return
       end if

       !total counter, for all the processes
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
          write(message,'(1x,a,f7.3,2(a,i0),a)')&
               'Limit of ',memorylimit,' GB reached, memproc ',memproc,' total memory is ',memtot%memory,' B. '

          if (f_err_raise(.true.,trim(message)//' Array '//trim(memtot%array)//&
               ', routine '//trim(memtot%routine),err_name='ERR_MEMLIMIT')) return
       end if

       select case(memproc)
       case (0)
          if (malloc_level ==2) then
             !to be used for inspecting an array which is not deallocated
             write(mallocFile,'(a,t40,a,t70,4(1x,i12))')trim(routine),trim(array),isize,memtot%memory
          else if (malloc_level ==1) then
             !Compact format
             if (trim(memloc%routine) /= routine) then
                if (memloc%memory /= int(0,kind=8)) then
                   rewind(mallocFile)
                   write(mallocFile,'(a,t40,a,t70,4(1x,i12))')&
                        trim(memloc%routine),trim(memloc%array),&
                        memloc%memory/int(1024,kind=8),memloc%peak/int(1024,kind=8),&
                        memtot%memory/int(1024,kind=8),&
                        (memtot%memory+memloc%peak-memloc%memory)/int(1024,kind=8)
                end if
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
          end if
       case default
          return
       end select
    end select
  END SUBROUTINE memocc_internal
  !!***


  !!****f* ABINIT/memory_malloc_check
  !! FUNCTION
  !!   Check the malloc.prc file (verbose format)
  !! SOURCE
  !!
  subroutine memory_malloc_check(nalloc,ndealloc)
    implicit none
    !Arguments
    integer, intent(in) :: nalloc,ndealloc
    !Local variables
    if (malloc_level==2 .and. nalloc /= ndealloc) then
       !Use # to be yaml compliant (is a comment in yaml)
       write(*,*) &
            "#Use the python script 'memcheck.py' in utils/scripts to check"//&
            trim(filename)//" file"
    end if
  END SUBROUTINE memory_malloc_check
  !!***

  !to be desactivated
  subroutine dp_padding(npaddim,nstart,array)
    implicit none
    integer, intent(in) :: npaddim,nstart
    double precision, dimension(*) :: array
    !local variables
    integer :: i
    do i=1,npaddim*ndebug
       array(nstart+i)= 0.0d0!d_nan() !this function is in profiling/memory.f90
    end do
  end subroutine dp_padding

  subroutine sp_padding(npaddim,nstart,array)
    implicit none
    integer, intent(in) :: npaddim,nstart
    real, dimension(*) :: array
    !local variables
    integer :: i
    do i=1,npaddim*ndebug
       array(nstart+i)=  0.0e0! r_nan() !this function is in profiling/memory.f90
    end do
  end subroutine sp_padding

  subroutine i_padding(npaddim,nstart,array)
    implicit none
    integer, intent(in) :: npaddim,nstart
    integer, dimension(*) :: array
    !local variables
    integer :: i
    do i=1,npaddim*ndebug
       array(nstart+i)= 0!int(r_nan()) !this function is in profiling/timem.f90
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

  !beginning of the verbose section
  subroutine mo_dp1(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    double precision, dimension(:), intent(in) :: array
    !local variables
    integer :: ndim
    if (ndebug /=0) then
       ndim=product(shape(array))-ndebug
       call dp_padding(1,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_dp1

  subroutine mo_dp2(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    double precision, dimension(:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(2) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:1))
       ndim=product(shape(array))-ndebug*npaddim
       call dp_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_dp2

  subroutine mo_dp3(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    double precision, dimension(:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(3) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:2))
       ndim=product(shape(array))-ndebug*npaddim
       call dp_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_dp3

  subroutine mo_dp4(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    double precision, dimension(:,:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(4) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:3))
       ndim=product(shape(array))-ndebug*npaddim
       call dp_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_dp4

  subroutine mo_dp5(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    double precision, dimension(:,:,:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(5) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:4))
       ndim=product(shape(array))-ndebug*npaddim
       call dp_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_dp5

  subroutine mo_dp6(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    double precision, dimension(:,:,:,:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(6) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:5))
       ndim=product(shape(array))-ndebug*npaddim
       call dp_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_dp6

  subroutine mo_dp7(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    double precision, dimension(:,:,:,:,:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(7) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:6))
       ndim=product(shape(array))-ndebug*npaddim
       call dp_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_dp7

  subroutine mo_sp1(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    real, dimension(:), intent(in) :: array
    !local variables
    integer :: ndim
    if (ndebug /=0) then
       ndim=product(shape(array))-ndebug
       call sp_padding(1,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_sp1

  subroutine mo_sp2(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    real, dimension(:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(2) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:1))
       ndim=product(shape(array))-ndebug*npaddim
       call sp_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_sp2

  subroutine mo_sp3(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    real, dimension(:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(3) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:2))
       ndim=product(shape(array))-ndebug*npaddim
       call sp_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_sp3

  subroutine mo_sp4(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    real, dimension(:,:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(4) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:3))
       ndim=product(shape(array))-ndebug*npaddim
       call sp_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_sp4

  subroutine mo_sp5(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    real, dimension(:,:,:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(5) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:4))
       ndim=product(shape(array))-ndebug*npaddim
       call sp_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_sp5

  subroutine mo_sp6(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    real, dimension(:,:,:,:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(6) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:5))
       ndim=product(shape(array))-ndebug*npaddim
       call sp_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_sp6

  subroutine mo_sp7(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    real, dimension(:,:,:,:,:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(7) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:6))
       ndim=product(shape(array))-ndebug*npaddim
       call sp_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_sp7

  subroutine mo_i1(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    integer, dimension(:), intent(in) :: array
    !local variables
    integer :: ndim
    if (ndebug /=0) then
       ndim=product(shape(array))-ndebug
       call i_padding(1,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_i1

  subroutine mo_i2(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    integer, dimension(:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(2) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:1))
       ndim=product(shape(array))-ndebug*npaddim
       call i_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_i2

  subroutine mo_i3(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    integer, dimension(:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(3) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:2))
       ndim=product(shape(array))-ndebug*npaddim
       call i_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_i3

  subroutine mo_i4(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    integer, dimension(:,:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(4) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:3))
       ndim=product(shape(array))-ndebug*npaddim
       call i_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_i4

  subroutine mo_i5(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    integer, dimension(:,:,:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(5) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:4))
       ndim=product(shape(array))-ndebug*npaddim
       call i_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_i5

  subroutine mo_i6(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    integer, dimension(:,:,:,:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(6) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:5))
       ndim=product(shape(array))-ndebug*npaddim
       call i_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_i6

  subroutine mo_i7(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    integer, dimension(:,:,:,:,:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(7) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:6))
       ndim=product(shape(array))-ndebug*npaddim
       call i_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_i7

  subroutine mo_l1(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    logical, dimension(:), intent(in) :: array
    !local variables
    integer :: ndim
    if (ndebug /=0) then
       ndim=product(shape(array))-ndebug
       call l_padding(1,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_l1

  subroutine mo_l2(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    logical, dimension(:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(2) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:1))
       ndim=product(shape(array))-ndebug*npaddim
       call l_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_l2

  subroutine mo_l3(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    logical, dimension(:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(3) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:2))
       ndim=product(shape(array))-ndebug*npaddim
       call l_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_l3

  subroutine mo_l4(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    logical, dimension(:,:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(4) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:3))
       ndim=product(shape(array))-ndebug*npaddim
       call l_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_l4

  subroutine mo_l5(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    logical, dimension(:,:,:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(5) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:4))
       ndim=product(shape(array))-ndebug*npaddim
       call l_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_l5

  subroutine mo_l6(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    logical, dimension(:,:,:,:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(6) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:5))
       ndim=product(shape(array))-ndebug*npaddim
       call l_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_l6

  subroutine mo_l7(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    logical, dimension(:,:,:,:,:,:,:), intent(in) :: array
    !local variables
    integer :: ndim,npaddim
    integer, dimension(7) :: iashp
    if (ndebug /=0) then
       iashp=shape(array)
       npaddim=product(iashp(1:6))
       ndim=product(shape(array))-ndebug*npaddim
       call l_padding(npaddim,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_l7

  subroutine mo_c1(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    character(len=*), dimension(:), intent(in) :: array
    !local variables
    integer :: ndim
    if (ndebug /=0) then
       ndim=product(shape(array))-ndebug
       call c_padding(1,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_c1

  subroutine mo_cmpdp1(istat,array,aname,rname)
    implicit none
    character(len=*), intent(in) :: aname,rname
    integer, intent(in) :: istat
    complex(kind=8), dimension(:), intent(in) :: array
    !local variables
    integer :: ndim
    if (ndebug /=0) then
       ndim=product(shape(array))-ndebug
       !stop "I don't have this function!!!!!"
       !call cmpdp_padding(1,ndim,array)
    end if
    call memocc_internal(istat,product(shape(array))*kind(array),aname,rname)
  end subroutine mo_cmpdp1

end module memory_profiling
