!!****m* ABINIT/m_profiling_abi
!! NAME
!! m_profiling_abi
!!
!! FUNCTION
!! This module contains routines to profile the code.
!! Currently only memory occupation are provided here.
!! Ideally, timab should be incorporated here.
!!
!! COPYRIGHT
!! Copyright (C) 2010-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

  module m_profiling_abi

    use defs_basis

    implicit none

    private

    !Memory profiling
    type :: memstat_abi
       character(len=36) :: routine,array
       integer(kind=8) :: memory,peak
    end type memstat_abi

    ! Save values for memocc_abi.
#if defined DEBUG_MODE
!    integer, parameter :: ndebug = 5  !5 will not work for wavelets compiling with debug=naughty
    integer,public, parameter :: ndebug = 0

#else
    integer,public, parameter :: ndebug = 0
#endif
    real,save :: memorylimit_abi = 0.e0
    logical,save :: meminit_abi = .false.
    integer, parameter :: mallocFile = 98
    type(memstat_abi),save :: memloc_abi,memtot_abi
    integer,save :: memalloc_abi,memdealloc_abi,memproc_abi = 0
    !Debug option for memocc_abi, set in the input file
    !logical :: memdebug=.true.
    integer,save :: malloc_level_abi=2

    integer, public :: ABI_ALLOC_STAT_ABI
    integer(kind=8), public :: ABI_ALLOC_SIZE_ABI

    !interface for the memory allocation control, depends on ndebug
    interface memocc_abi
       module procedure mo_abi_dp1,mo_abi_dp2,mo_abi_dp3,mo_abi_dp4,mo_abi_dp5,mo_abi_dp6,mo_abi_dp7,&
            mo_abi_sp1,mo_abi_sp2,mo_abi_sp3,mo_abi_sp4,mo_abi_sp5,mo_abi_sp6,mo_abi_sp7,&
            mo_abi_i1,mo_abi_i2,mo_abi_i3,mo_abi_i4,mo_abi_i5,mo_abi_i6,mo_abi_i7,&
            mo_abi_l1,mo_abi_l2,mo_abi_l3,mo_abi_l4,mo_abi_l5,mo_abi_l6,mo_abi_l7,&
            mo_abi_c1, mo_abi_cmpdp1, &
            memocc_abi_internal  !central routine to be used for deallocation
    end interface

    public :: memocc_abi
    public :: memocc_abi_get_info
    public :: memocc_abi_set_state
    public :: memocc_abi_set_memory_limit
    public :: memocc_abi_report
!!***

  contains

    !!****f* m_profiling_abi/memocc_abi_set_state
    !! NAME
    !! memocc_abi_set_state
    !!
    !! FUNCTION
    !!  @param istatus 0 no file malloc.prc is created, only memory allocation counters running
    !!                 1 file malloc.prc is created in a light version (only current information is written)
    !!                 2 file malloc.prc is created with full information inside (default state if not specified)
    !! The status can only be downgraded. A stop signal is produced if status is increased
    subroutine memocc_abi_set_state(istatus)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memocc_abi_set_state'
!End of the abilint section

      integer, intent(in) :: istatus

      if (istatus > malloc_level_abi) stop 'malloc_level_abi can be only downgraded'

      malloc_level_abi = istatus
      
      if (istatus == 2) return !the default situation

      if (istatus == 1) then 
         !clean the file situation
         close(unit=mallocFile)
         open(unit=mallocFile,file='malloc.prc',status='unknown',action='write')
      end if

      if (istatus == 0) then
         !the file should be deleted
         close(unit=mallocFile)
         open(unit=mallocFile,file='malloc.prc',status='replace')
         close(unit=mallocFile)
      end if
    end subroutine memocc_abi_set_state
    !!***

    !!****f* m_profiling_abi/memocc_abi_set_memory_limit
    !! NAME
    !! memocc_abi_set_memory_limit
    !!
    !! FUNCTION
    !!  @param limit Give a memory limit above which the code will stop
    !!               properly. The unit is bytes.
    subroutine memocc_abi_set_memory_limit(limit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memocc_abi_set_memory_limit'
!End of the abilint section

      real, intent(in) :: limit

      memorylimit_abi = limit
    end subroutine memocc_abi_set_memory_limit
    !!***

    !> Put to zero memocc_abi counters
    subroutine memocc_abi_variables_init()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memocc_abi_variables_init'
!End of the abilint section

      memtot_abi%memory=int(0,kind=8)
      memtot_abi%peak=int(0,kind=8)
      memalloc_abi=0
      memdealloc_abi=0

      memloc_abi%routine='routine'
      memloc_abi%array='array'
      memloc_abi%memory=int(0,kind=8) !fake initialisation to print the first routine
      memloc_abi%peak=int(0,kind=8)
    end subroutine memocc_abi_variables_init

    subroutine memocc_abi_report()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memocc_abi_report'
!End of the abilint section

      call memocc_abi(0,int(0,kind=8),'count', 'stop')
    end subroutine memocc_abi_report

    !!****f* m_profiling_abi/memocc_abi_get_info
    !! NAME
    !! memocc_abi_get_info
    !!
    !! FUNCTION
    !!  Function that returns the number of allocations and deallocations that have
    !!  been done and the memory currently used
    !! INPUT VARIABLES
    !!
    !! OUTPUT VARIABLES
    !!  @param nalloc       number of allocations that have been done
    !!  @param ndealloc     number of deallocations that have been done
    !!  @param allocmemory  total memory used
    subroutine memocc_abi_get_info(nalloc,ndealloc,allocmemory)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memocc_abi_get_info'
!End of the abilint section

      integer(kind=8), intent(out) :: allocmemory
      integer, intent(out) :: nalloc,ndealloc

      nalloc = memalloc_abi
      ndealloc = memdealloc_abi
      allocmemory = memtot_abi%memory

    end subroutine memocc_abi_get_info
    !!***

    subroutine memocc_abi_open_file()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memocc_abi_open_file'
!End of the abilint section

      open(unit=mallocFile,file='malloc.prc',status='unknown')
      !if (memdebug) then
         write(mallocFile,'(a,t60,a,t90,4(1x,a12))')&
              '(Data in KB) Routine','Array name    ',&
              'Array size','Total Memory'
      !else
      !write(mallocFile,'(a,t60,a,t90,4(1x,a12))')&
      !     '(Data in KB) Routine','Peak Array    ',&
      !     'Routine Mem','Routine Peak','Memory Stat.','Memory Peak'
      !end if
    end subroutine memocc_abi_open_file


    !!****f* ABINIT/memory_occupation_abi
    !! FUNCTION
    !! Control the memory occupation by calculating the overall size of the allocated arrays
    !! DESCRIPTION
    !!   when allocating allocating an array "stuff" of dimension n in the routine "dosome":
    !!      allocate(stuff(n),stat=i_stat)
    !!      call memocc_abi(i_stat,product(int(shape(stuff),kind=8))*kind(stuff),'stuff','dosome')
    !!   when deallocating:
    !!      i_all=-product(int(shape(stuff),kind=8))*kind(stuff)
    !!      deallocate(stuff,stat=i_stat)
    !!      call memocc_abi(i_stat,i_all,'stuff','dosome')
    !!   The counters are initialized once by the first allocation
    !!   and stopped with:
    !!      call memocc_abi(0,0,'count','stop')
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
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine memory_occupation_abi(istat,isize,array,routine)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memory_occupation_abi'
!End of the abilint section

      implicit none

#if defined HAVE_MPI
 include 'mpif.h'
#endif

      ! Arguments
      integer, intent(in) :: istat
      integer(kind=8), intent(in) :: isize
      character(len=*), intent(in) :: array,routine

      ! Local variables
      logical :: lmpinit
      integer :: ierr

      !write(std_out,*) memproc_abi,array,routine
      ! Initialised first
      if (.not.meminit_abi) then

         call memocc_abi_variables_init()
#if defined HAVE_MPI
         !Use MPI to have the mpi rank
         call MPI_INITIALIZED(lmpinit,ierr)
         if (lmpinit) then
            call MPI_COMM_RANK(MPI_COMM_WORLD,memproc_abi,ierr)
         else
            !no-mpi case 
            memproc_abi=0
         end if
#else
         memproc_abi = 0
#endif

         !open the writing file for the root process
         if (memproc_abi == 0 .and. malloc_level_abi > 0) then
            call memocc_abi_open_file()
         end if
         meminit_abi = .true.
      end if

      select case(array)
      case('count')
         if (trim(routine)=='stop' .and. memproc_abi==0) then
            if (malloc_level_abi > 0) then
               if (malloc_level_abi == 1) rewind(mallocFile)
               write(mallocFile,'(a,t60,a,t90,4(1x,i12))')&
                    trim(memloc_abi%routine),trim(memloc_abi%array),&
                    memloc_abi%memory/int(1024,kind=8),memloc_abi%peak/int(1024,kind=8),&
                    memtot_abi%memory/int(1024,kind=8),&
                    (memtot_abi%peak+memloc_abi%peak-memloc_abi%memory)/int(1024,kind=8)
               close(unit=mallocFile)
            end if
             write(std_out,'(1x,a)')&
                  '-------------------------MEMORY CONSUMPTION REPORT-----------------------------'
             write(std_out,'(1x,2(i0,a,1x),i0)')&
                  memalloc_abi,' allocations and',memdealloc_abi,' deallocations, remaining memory(B):',&
                  memtot_abi%memory
             write(std_out,'(1x,a,i0,a)') 'memory occupation peak: ',memtot_abi%peak/int(1048576,kind=8),' MB'
             write(std_out,'(4(1x,a))') 'for the array ',trim(memtot_abi%array),&
                  'in the routine',trim(memtot_abi%routine)
               !here we can add a routine which open the malloc.prc file in case of some
               !memory allocation problem, and which eliminates it for a successful run
            if (malloc_level_abi == 1 .and. memalloc_abi == memdealloc_abi .and. memtot_abi%memory==int(0,kind=8)) then
               !remove file should be put here
               open(unit=mallocFile,file='malloc.prc',status='unknown',action='write')
               write(unit=mallocFile,fmt='()',advance='no')
               close(unit=mallocFile)
            else
               call memory_malloc_check_abi(memalloc_abi,memdealloc_abi)
            end if
         else if (trim(routine)/='stop') then
            write(std_out,*) "memocc_abi: ",array," ",routine
            write(std_out,"(a,i0,a)") "Error[",memproc_abi,"]: Use memocc_abi and the word 'count' only with the word 'stop'."
            stop
         end if

      case default
         !control of the allocation/deallocation status
         if (istat/=0) then
            if (isize>=0) then
               write(std_out,*)' subroutine ',routine,': problem of allocation of array ',array,&
                    ', error code=',istat,' exiting...'
               if (memproc_abi == 0 .and. malloc_level_abi > 0) close(unit=mallocFile)
#if defined HAVE_MPI
               call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
               stop
#endif
            else if (isize<0) then
               write(std_out,*)' subroutine ',routine,': problem of deallocation of array ',array,&
                    ', error code=',istat,' exiting...'
               if (memproc_abi == 0 .and. malloc_level_abi > 0) close(unit=mallocFile)
#if defined HAVE_MPI
               call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
               stop
#endif
            end if
         end if
         !total counter, for all the processes
         memtot_abi%memory=memtot_abi%memory+isize
         if (memtot_abi%memory > memtot_abi%peak) then
            memtot_abi%peak=memtot_abi%memory
            memtot_abi%routine=routine
            memtot_abi%array=array
         end if
         if (isize>0) then
            memalloc_abi=memalloc_abi+1
         else if (isize<0) then
            memdealloc_abi=memdealloc_abi+1
         end if

         if (memorylimit_abi /= 0.e0 .and. &
              memtot_abi%memory > int(real(memorylimit_abi,kind=8)*1073741824.d0,kind=8)) then !memory limit is in GB
            write(std_out,'(1x,a,f7.3,2(a,i0),a)')&
                 'ERROR: Memory limit of ',memorylimit_abi,&
                 ' GB reached for memproc ',memproc_abi,' : total memory is ',memtot_abi%memory,' B.'
            write(std_out,'(1x,2(a,i0))')&
                 '       this happened for array '//trim(memtot_abi%array)//' in routine '//trim(memtot_abi%routine)
#if defined HAVE_MPI
               call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
               stop
#endif
         end if

         select case(memproc_abi)
         case (0)
            if (malloc_level_abi ==2) then
               !to be used for inspecting an array which is not deallocated
               write(98,'(a,t60,a,t90,4(1x,i12))')trim(routine),trim(array),isize,memtot_abi%memory
            else if (malloc_level_abi ==1) then
               !Compact format
               if (trim(memloc_abi%routine) /= routine) then
                  if (memloc_abi%memory /= int(0,kind=8)) then
                     rewind(mallocFile)
                     write(mallocFile,'(a,t60,a,t90,4(1x,i12))')&
                          trim(memloc_abi%routine),trim(memloc_abi%array),&
                          memloc_abi%memory/int(1024,kind=8),memloc_abi%peak/int(1024,kind=8),&
                          memtot_abi%memory/int(1024,kind=8),&
                          (memtot_abi%memory+memloc_abi%peak-memloc_abi%memory)/int(1024,kind=8)
                  end if
                  memloc_abi%routine=routine
                  memloc_abi%array=array
                  memloc_abi%memory=isize
                  memloc_abi%peak=isize
               else
                  memloc_abi%memory=memloc_abi%memory+isize
                  if (memloc_abi%memory > memloc_abi%peak) then
                     memloc_abi%peak=memloc_abi%memory
                     memloc_abi%array=array
                  end if
               end if
            end if
         case default
            return
         end select
      end select
    END SUBROUTINE memory_occupation_abi
    !!***


    !!****f* m_profiling_abi/memory_malloc_check_abi
    !! NAME
    !! memory_malloc_check_abi
    !!
    !! FUNCTION
    !!   Check the malloc.prc file (verbose format)
!! PARENTS
!!      m_profiling_abi
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine memory_malloc_check_abi(nalloc,ndealloc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memory_malloc_check_abi'
!End of the abilint section

      implicit none
      !Arguments
      integer, intent(in) :: nalloc,ndealloc
      !Local variables
      if (malloc_level_abi==2 .and. nalloc /= ndealloc) then
         write(std_out,*) &
              "Use the python script 'memcheck.py' in utils/scripts to check 'malloc.prc' file"
      end if
    END SUBROUTINE memory_malloc_check_abi
    !!***


    !!****f* m_profiling_abi/d_nan_abi
    !! NAME
    !! d_nan_abi
    !!
    !! FUNCTION
    !!   Function which specify NaN according to IEEE specifications
    !! SOURCE
    !!
    function d_nan_abi()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd_nan_abi'
!End of the abilint section

      implicit none
      double precision :: d_nan_abi
      !local variables
      double precision :: dnan
      integer, dimension(2) :: inan
      equivalence (dnan, inan)
      ! This first assignment is for big-endian machines
      inan(1) = 2147483647
      ! The second assignment is for little-endian machines
      inan(2) = 2147483647
      d_nan_abi = dnan
    end function d_nan_abi
    !!***

    !!****f* m_profiling_abi/r_nan_abi
    !! NAME
    !! r_nan_abi
    !!
    !! FUNCTION
    !!   Function which specify NaN according to IEEE specifications
    !! SOURCE
    !!
    function r_nan_abi()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'r_nan_abi'
!End of the abilint section

      implicit none
      real :: r_nan_abi
      !local variables
      real :: rnan
      integer :: inan
      equivalence (rnan, inan)
      inan = 2147483647
      r_nan_abi = rnan
    end function r_nan_abi
    !!***

    !!****f* m_profiling_abi/memocc_abi_internal
    !! NAME
    !! memocc_abi_internal
    !!
    !! FUNCTION
    !!   routine used for deallocations
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine memocc_abi_internal(istat,isize,array,routine)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memocc_abi_internal'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: array,routine
      integer, intent(in) :: istat
      integer(kind=8), intent(in) :: isize
      call memory_occupation_abi(istat,isize,array,routine)
    end subroutine memocc_abi_internal
    !!***

    !!****f* m_profiling_abi/dp_padding_abi
    !! NAME
    !! dp_padding_abi
    !!
    !! FUNCTION
    !!   Pad the end of the array with NAN values.
!! PARENTS
!!      m_profiling_abi
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine dp_padding_abi(npaddim,nstart,array)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dp_padding_abi'
!End of the abilint section

      implicit none
      integer, intent(in) :: npaddim,nstart
      double precision, dimension(*) :: array
      !local variables
      integer :: i,nstart_
      nstart_=max(nstart,0)
      do i=1,npaddim*ndebug
         array(nstart_+i)= d_nan_abi()
      end do
    end subroutine dp_padding_abi
    !!***

    !!****f* m_profiling_abi/sp_padding_abi
    !! NAME
    !! sp_padding_abi
    !!
    !! FUNCTION
    !!   Pad the end of the array with NAN values.
!! PARENTS
!!      m_profiling_abi
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine sp_padding_abi(npaddim,nstart,array)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sp_padding_abi'
!End of the abilint section

      implicit none
      integer, intent(in) :: npaddim,nstart
      real, dimension(*) :: array
      !local variables
      integer :: i,nstart_
      nstart_=max(nstart,0)
      do i=1,npaddim*ndebug
         array(nstart_+i)= r_nan_abi()
      end do
    end subroutine sp_padding_abi
    !!***

    !!****f* m_profiling_abi/i_padding_abi
    !! NAME
    !! i_padding_abi
    !!
    !! FUNCTION
    !!   Pad the end of the array with NAN values.
!! PARENTS
!!      m_profiling_abi
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine i_padding_abi(npaddim,nstart,array)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'i_padding_abi'
!End of the abilint section

      implicit none
      integer, intent(in) :: npaddim,nstart
      integer, dimension(*) :: array
      !local variables
      integer :: i,nstart_
      nstart_=max(nstart,0)
      do i=1,npaddim*ndebug
         array(nstart_+i)= int(r_nan_abi())
      end do
    end subroutine i_padding_abi
    !!***

    !!****f* m_profiling_abi/l_padding_abi
    !! NAME
    !! l_padding_abi
    !!
    !! FUNCTION
    !!   Pad the end of the array with .false. values.
!! PARENTS
!!      m_profiling_abi
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine l_padding_abi(npaddim,nstart,array)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'l_padding_abi'
!End of the abilint section

      implicit none
      integer, intent(in) :: npaddim,nstart
      logical, dimension(*) :: array
      !local variables
      integer :: i,nstart_
      nstart_=max(nstart,0)
      do i=1,npaddim*ndebug
         array(nstart_+i)=.false.
      end do
    end subroutine l_padding_abi
    !!***

    !!****f* m_profiling_abi/c_padding_abi
    !! NAME
    !! c_padding_abi
    !!
    !! FUNCTION
    !!   Pad the end of the array with character values.
!! PARENTS
!!      m_profiling_abi
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine c_padding_abi(npaddim,nstart,array)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'c_padding_abi'
!End of the abilint section

      implicit none
      integer, intent(in) :: npaddim,nstart
      character(len=20), dimension(*) :: array
      !local variables
      integer :: i,nstart_
      nstart_=max(nstart,0)
      do i=1,npaddim*ndebug
         array(nstart_+i)='AAAAAAAAAAAAAAAAAAAA'
      end do
    end subroutine c_padding_abi
    !!***

    !!****f* m_profiling_abi/mo_abi_dp1
    !! NAME
    !! mo_abi_dp1
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_dp1(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_dp1'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      double precision, dimension(:), intent(in) :: array
      !local variables
      integer :: ndim
      if (ndebug /=0) then
         ndim=product(int(shape(array),kind=8))-ndebug
         call dp_padding_abi(1,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_dp1
    !!***

    !!****f* m_profiling_abi/mo_abi_dp2
    !! NAME
    !! mo_abi_dp2
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_dp2(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_dp2'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      double precision, dimension(:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(2) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:1))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call dp_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_dp2
    !!***

    !!****f* m_profiling_abi/mo_abi_dp3
    !! NAME
    !! mo_abi_dp3
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_dp3(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_dp3'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      double precision, dimension(:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(3) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:2))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call dp_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_dp3
    !!***

    !!****f* m_profiling_abi/mo_abi_dp4
    !! NAME
    !! mo_abi_dp4
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_dp4(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_dp4'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      double precision, dimension(:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(4) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:3))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call dp_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_dp4
    !!***

    !!****f* m_profiling_abi/mo_abi_dp5
    !! NAME
    !! mo_abi_dp5
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_dp5(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_dp5'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      double precision, dimension(:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(5) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:4))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call dp_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_dp5
    !!***

    !!****f* m_profiling_abi/mo_abi_dp6
    !! NAME
    !! mo_abi_dp6
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_dp6(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_dp6'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      double precision, dimension(:,:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(6) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:5))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call dp_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_dp6
    !!***

    !!****f* m_profiling_abi/mo_abi_dp7
    !! NAME
    !! mo_abi_dp7
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_dp7(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_dp7'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      double precision, dimension(:,:,:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(7) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:6))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call dp_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_dp7
    !!***

    !!****f* m_profiling_abi/mo_abi_sp1
    !! NAME
    !! mo_abi_sp1
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_sp1(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_sp1'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real, dimension(:), intent(in) :: array
      !local variables
      integer :: ndim
      if (ndebug /=0) then
         ndim=product(int(shape(array),kind=8))-ndebug
         call sp_padding_abi(1,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_sp1
    !!***

    !!****f* m_profiling_abi/mo_abi_sp2
    !! NAME
    !! mo_abi_sp2
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_sp2(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_sp2'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real, dimension(:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(2) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:1))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call sp_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_sp2
    !!***

    !!****f* m_profiling_abi/mo_abi_sp3
    !! NAME
    !! mo_abi_sp3
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_sp3(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_sp3'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real, dimension(:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(3) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:2))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call sp_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_sp3
    !!***

    !!****f* m_profiling_abi/mo_abi_sp4
    !! NAME
    !! mo_abi_sp4
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_sp4(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_sp4'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real, dimension(:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(4) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:3))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call sp_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_sp4
    !!***

    !!****f* m_profiling_abi/mo_abi_sp5
    !! NAME
    !! mo_abi_sp5
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_sp5(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_sp5'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real, dimension(:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(5) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:4))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call sp_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_sp5
    !!***

    !!****f* m_profiling_abi/mo_abi_sp6
    !! NAME
    !! mo_abi_sp6
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_sp6(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_sp6'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real, dimension(:,:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(6) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:5))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call sp_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_sp6
    !!***

    !!****f* m_profiling_abi/mo_abi_sp7
    !! NAME
    !! mo_abi_sp7
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_sp7(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_sp7'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real, dimension(:,:,:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(7) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:6))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call sp_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_sp7
    !!***

    !!****f* m_profiling_abi/mo_abi_i1
    !! NAME
    !! mo_abi_i1
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_i1(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_i1'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:), intent(in) :: array
      !local variables
      integer :: ndim
      if (ndebug /=0) then
         ndim=product(int(shape(array),kind=8))-ndebug
         call i_padding_abi(1,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_i1
    !!***

    !!****f* m_profiling_abi/mo_abi_i2
    !! NAME
    !! mo_abi_i2
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_i2(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_i2'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(2) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:1))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call i_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_i2
    !!***

    !!****f* m_profiling_abi/mo_abi_i3
    !! NAME
    !! mo_abi_i3
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_i3(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_i3'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(3) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:2))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call i_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_i3
    !!***

    !!****f* m_profiling_abi/mo_abi_i4
    !! NAME
    !! mo_abi_i4
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_i4(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_i4'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(4) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:3))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call i_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_i4
    !!***

    !!****f* m_profiling_abi/mo_abi_i5
    !! NAME
    !! mo_abi_i5
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_i5(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_i5'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(5) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:4))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call i_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_i5
    !!***

    !!****f* m_profiling_abi/mo_abi_i6
    !! NAME
    !! mo_abi_i6
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_i6(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_i6'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(6) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:5))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call i_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_i6
    !!***

    !!****f* m_profiling_abi/mo_abi_i7
    !! NAME
    !! mo_abi_i7
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_i7(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_i7'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:,:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(7) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:6))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call i_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_i7
    !!***

    !!****f* m_profiling_abi/mo_abi_l1
    !! NAME
    !! mo_abi_l1
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_l1(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_l1'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      logical, dimension(:), intent(in) :: array
      !local variables
      integer :: ndim
      if (ndebug /=0) then
         ndim=product(int(shape(array),kind=8))-ndebug
         call l_padding_abi(1,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_l1
    !!***

    !!****f* m_profiling_abi/mo_abi_l2
    !! NAME
    !! mo_abi_l2
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_l2(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_l2'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      logical, dimension(:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(2) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:1))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call l_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_l2
    !!***

    !!****f* m_profiling_abi/mo_abi_l3
    !! NAME
    !! mo_abi_l3
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_l3(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_l3'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      logical, dimension(:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(3) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:2))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call l_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_l3
    !!***

    !!****f* m_profiling_abi/mo_abi_l4
    !! NAME
    !! mo_abi_l4
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_l4(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_l4'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      logical, dimension(:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(4) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:3))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call l_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_l4
    !!***

    !!****f* m_profiling_abi/mo_abi_l5
    !! NAME
    !! mo_abi_l5
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_l5(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_l5'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      logical, dimension(:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(5) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:4))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call l_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_l5
    !!***

    !!****f* m_profiling_abi/mo_abi_l6
    !! NAME
    !! mo_abi_l6
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_l6(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_l6'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      logical, dimension(:,:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(6) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:5))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call l_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_l6
    !!***

    !!****f* m_profiling_abi/mo_abi_l7
    !! NAME
    !! mo_abi_l7
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_l7(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_l7'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      logical, dimension(:,:,:,:,:,:,:), intent(in) :: array
      !local variables
      integer :: ndim,npaddim
      integer, dimension(7) :: iashp
      if (ndebug /=0) then
         iashp=int(shape(array),kind=8)
         npaddim=product(iashp(1:6))
         ndim=product(int(shape(array),kind=8))-ndebug*npaddim
         call l_padding_abi(npaddim,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_l7
    !!***

    !!****f* m_profiling_abi/mo_abi_c1
    !! NAME
    !! mo_abi_c1
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_c1(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_c1'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      character(len=*), dimension(:), intent(in) :: array
      !local variables
      integer :: ndim
      if (ndebug /=0) then
         ndim=product(int(shape(array),kind=8))-ndebug
         call c_padding_abi(1,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_c1
    !!***

    !!****f* m_profiling_abi/mo_abi_cmpdp1
    !! NAME
    !! mo_abi_cmpdp1
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation_abi
!!
    !! SOURCE
    !!
    subroutine mo_abi_cmpdp1(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_abi_cmpdp1'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      complex(kind=8), dimension(:), intent(in) :: array
      !local variables
      integer :: ndim
      if (ndebug /=0) then
         ndim=product(int(shape(array),kind=8))-ndebug
         !TODO: to be implemented.
         !call cmpdp_padding_abi(1,ndim,array)
      end if
      call memory_occupation_abi(istat,product(int(shape(array),kind=8))*kind(array),aname,rname)
    end subroutine mo_abi_cmpdp1
    !!***

end module m_profiling_abi
