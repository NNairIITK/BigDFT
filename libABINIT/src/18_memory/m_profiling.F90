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
!! Copyright (C) 2010-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

  module m_profiling

    use defs_basis

    implicit none

    private

    !Memory profiling
    type :: memstat
       character(len=36) :: routine,array
       integer(kind=8) :: memory,peak
    end type memstat

    ! Save values for memocc.
#if defined DEBUG_MODE
!    integer, parameter :: ndebug = 5  !5 will not work for wavelets compiling with debug=naughty
    integer,public, parameter :: ndebug = 0

#else
    integer,public, parameter :: ndebug = 0
#endif
    real,save :: memorylimit = 0.e0
    logical,save :: meminit = .false.
    integer, parameter :: mallocFile = 98
    type(memstat),save :: memloc,memtot
    integer,save :: memalloc,memdealloc,memproc = 0
    !Debug option for memocc, set in the input file
    !logical :: memdebug=.true.
    integer,save :: malloc_level=2

    integer, public :: ABI_ALLOC_STAT, ABI_ALLOC_SIZE
    ! MG: These are not a good names for  public variables that are exported 
    integer, public :: sz1, sz2, sz3, sz4, sz5, sz6, sz7

    !interface for the memory allocation control, depends on ndebug
    interface memocc
       module procedure mo_dp1,mo_dp2,mo_dp3,mo_dp4,mo_dp5,mo_dp6,mo_dp7,&
            mo_sp1,mo_sp2,mo_sp3,mo_sp4,mo_sp5,mo_sp6,mo_sp7,&
            mo_i1,mo_i2,mo_i3,mo_i4,mo_i5,mo_i6,mo_i7,&
            mo_l1,mo_l2,mo_l3,mo_l4,mo_l5,mo_l6,mo_l7,&
            mo_c1, mo_cmpdp1, &
            memocc_internal  !central routine to be used for deallocation
    end interface

!    public :: memocc
    public :: memocc_get_info
    public :: memocc_set_state
    public :: memocc_set_memory_limit
    public :: memocc_report
    public :: d_nan,r_nan
!!***

  contains

    !!****f* m_profiling/memocc_set_state
    !! NAME
    !! memocc_set_state
    !!
    !! FUNCTION
    !!  @param istatus 0 no file malloc.prc is created, only memory allocation counters running
    !!                 1 file malloc.prc is created in a light version (only current information is written)
    !!                 2 file malloc.prc is created with full information inside (default state if not specified)
    !! The status can only be downgraded. A stop signal is produced if status is increased
    subroutine memocc_set_state(istatus)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memocc_set_state'
!End of the abilint section

      integer, intent(in) :: istatus

      if (istatus > malloc_level) stop 'malloc_level can be only downgraded'

      malloc_level = istatus
      
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
    end subroutine memocc_set_state
    !!***

    !!****f* m_profiling/memocc_set_memory_limit
    !! NAME
    !! memocc_set_memory_limit
    !!
    !! FUNCTION
    !!  @param limit Give a memory limit above which the code will stop
    !!               properly. The unit is bytes.
    subroutine memocc_set_memory_limit(limit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memocc_set_memory_limit'
!End of the abilint section

      real, intent(in) :: limit

      memorylimit = limit
    end subroutine memocc_set_memory_limit
    !!***

    !> Put to zero memocc counters
    subroutine memocc_variables_init()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memocc_variables_init'
!End of the abilint section

      memtot%memory=int(0,kind=8)
      memtot%peak=int(0,kind=8)
      memalloc=0
      memdealloc=0

      memloc%routine='routine'
      memloc%array='array'
      memloc%memory=int(0,kind=8) !fake initialisation to print the first routine
      memloc%peak=int(0,kind=8)
    end subroutine memocc_variables_init

    subroutine memocc_report()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memocc_report'
!End of the abilint section

      call memocc(0,0,'count', 'stop')
    end subroutine memocc_report

    !!****f* m_profiling/memocc_get_info
    !! NAME
    !! memocc_get_info
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
    subroutine memocc_get_info(nalloc,ndealloc,allocmemory)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memocc_get_info'
!End of the abilint section

      integer(kind=8), intent(out) :: allocmemory
      integer, intent(out) :: nalloc,ndealloc

      nalloc = memalloc
      ndealloc = memdealloc
      allocmemory = memtot%memory

    end subroutine memocc_get_info
    !!***

    subroutine memocc_open_file()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memocc_open_file'
!End of the abilint section

      open(unit=mallocFile,file='malloc.prc',status='unknown')
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



    !!****f* ABINIT/memory_occupation
    !! FUNCTION
    !! Control the memory occupation by calculating the overall size of the allocated arrays
    !! DESCRIPTION
    !!   when allocating allocating an array "stuff" of dimension n in the routine "dosome":
    !!      allocate(stuff(n),stat=i_stat)
    !!      call memocc(i_stat,product(shape(stuff))*kind(stuff),'stuff','dosome')
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
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine memory_occupation(istat,isize,array,routine)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memory_occupation'
!End of the abilint section

      implicit none

#if defined HAVE_MPI
 include 'mpif.h'
#endif

      ! Arguments
      integer, intent(in) :: istat,isize
      character(len=*), intent(in) :: array,routine

      ! Local variables
      logical :: lmpinit
      integer :: ierr

      !write(std_out,*) memproc,array,routine
      ! Initialised first
      if (.not.meminit) then

         call memocc_variables_init()
#if defined HAVE_MPI
         !Use MPI to have the mpi rank
         call MPI_INITIALIZED(lmpinit,ierr)
         if (lmpinit) then
            call MPI_COMM_RANK(MPI_COMM_WORLD,memproc,ierr)
         else
            !no-mpi case 
            memproc=0
         end if
#else
         memproc = 0
#endif

         !open the writing file for the root process
         if (memproc == 0 .and. malloc_level > 0) then
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
             write(std_out,'(1x,a)')&
                  '-------------------------MEMORY CONSUMPTION REPORT-----------------------------'
             write(std_out,'(1x,2(i0,a,1x),i0)')&
                  memalloc,' allocations and',memdealloc,' deallocations, remaining memory(B):',&
                  memtot%memory
             write(std_out,'(1x,a,i0,a)') 'memory occupation peak: ',memtot%peak/int(1048576,kind=8),' MB'
             write(std_out,'(4(1x,a))') 'for the array ',trim(memtot%array),&
                  'in the routine',trim(memtot%routine)
               !here we can add a routine which open the malloc.prc file in case of some
               !memory allocation problem, and which eliminates it for a successful run
            if (malloc_level == 1 .and. memalloc == memdealloc .and. memtot%memory==int(0,kind=8)) then
               !remove file should be put here
               open(unit=mallocFile,file='malloc.prc',status='unknown',action='write')
               write(unit=mallocFile,fmt='()',advance='no')
               close(unit=mallocFile)
            else
               call memory_malloc_check(memalloc,memdealloc)
            end if
         else if (trim(routine)/='stop') then
            write(std_out,*) "memocc: ",array," ",routine
            write(std_out,"(a,i0,a)") "Error[",memproc,"]: Use memocc and the word 'count' only with the word 'stop'."
            stop
         end if

      case default
         !control of the allocation/deallocation status
         if (istat/=0) then
            if (isize>=0) then
               write(std_out,*)' subroutine ',routine,': problem of allocation of array ',array,&
                    ', error code=',istat,' exiting...'
               if (memproc == 0 .and. malloc_level > 0) close(unit=mallocFile)
#if defined HAVE_MPI
               call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
               stop
#endif
            else if (isize<0) then
               write(std_out,*)' subroutine ',routine,': problem of deallocation of array ',array,&
                    ', error code=',istat,' exiting...'
               if (memproc == 0 .and. malloc_level > 0) close(unit=mallocFile)
#if defined HAVE_MPI
               call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
               stop
#endif
            end if
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
            write(std_out,'(1x,a,f7.3,2(a,i0),a)')&
                 'ERROR: Memory limit of ',memorylimit,&
                 ' GB reached for memproc ',memproc,' : total memory is ',memtot%memory,' B.'
            write(std_out,'(1x,2(a,i0))')&
                 '       this happened for array '//trim(memtot%array)//' in routine '//trim(memtot%routine)
#if defined HAVE_MPI
               call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
               stop
#endif
         end if

         select case(memproc)
         case (0)
            if (malloc_level ==2) then
               !to be used for inspecting an array which is not deallocated
               write(98,'(a,t40,a,t70,4(1x,i12))')trim(routine),trim(array),isize,memtot%memory
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
    END SUBROUTINE memory_occupation
    !!***


    !!****f* m_profiling/memory_malloc_check
    !! NAME
    !! memory_malloc_check
    !!
    !! FUNCTION
    !!   Check the malloc.prc file (verbose format)
!! PARENTS
!!      m_profiling
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine memory_malloc_check(nalloc,ndealloc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memory_malloc_check'
!End of the abilint section

      implicit none
      !Arguments
      integer, intent(in) :: nalloc,ndealloc
      !Local variables
      if (malloc_level==2 .and. nalloc /= ndealloc) then
         write(std_out,*) &
              "Use the python script 'memcheck.py' in utils/scripts to check 'malloc.prc' file"
      end if
    END SUBROUTINE memory_malloc_check
    !!***


    !!****f* m_profiling/d_nan
    !! NAME
    !! d_nan
    !!
    !! FUNCTION
    !!   Function which specify NaN according to IEEE specifications
    !! SOURCE
    !!
    function d_nan()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd_nan'
!End of the abilint section

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

    !!****f* m_profiling/r_nan
    !! NAME
    !! r_nan
    !!
    !! FUNCTION
    !!   Function which specify NaN according to IEEE specifications
    !! SOURCE
    !!
    function r_nan()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'r_nan'
!End of the abilint section

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

    !!****f* m_profiling/memocc_internal
    !! NAME
    !! memocc_internal
    !!
    !! FUNCTION
    !!   routine used for deallocations
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine memocc_internal(istat,isize,array,routine)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'memocc_internal'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: array,routine
      integer, intent(in) :: istat,isize
      call memory_occupation(istat,isize,array,routine)
    end subroutine memocc_internal
    !!***

    !!****f* m_profiling/dp_padding
    !! NAME
    !! dp_padding
    !!
    !! FUNCTION
    !!   Pad the end of the array with NAN values.
!! PARENTS
!!      m_profiling
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine dp_padding(npaddim,nstart,array)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dp_padding'
!End of the abilint section

      implicit none
      integer, intent(in) :: npaddim,nstart
      double precision, dimension(*) :: array
      !local variables
      integer :: i,nstart_
      nstart_=max(nstart,0)
      do i=1,npaddim*ndebug
         array(nstart_+i)= d_nan()
      end do
    end subroutine dp_padding
    !!***

    !!****f* m_profiling/sp_padding
    !! NAME
    !! sp_padding
    !!
    !! FUNCTION
    !!   Pad the end of the array with NAN values.
!! PARENTS
!!      m_profiling
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine sp_padding(npaddim,nstart,array)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sp_padding'
!End of the abilint section

      implicit none
      integer, intent(in) :: npaddim,nstart
      real, dimension(*) :: array
      !local variables
      integer :: i,nstart_
      nstart_=max(nstart,0)
      do i=1,npaddim*ndebug
         array(nstart_+i)= r_nan()
      end do
    end subroutine sp_padding
    !!***

    !!****f* m_profiling/i_padding
    !! NAME
    !! i_padding
    !!
    !! FUNCTION
    !!   Pad the end of the array with NAN values.
!! PARENTS
!!      m_profiling
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine i_padding(npaddim,nstart,array)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'i_padding'
!End of the abilint section

      implicit none
      integer, intent(in) :: npaddim,nstart
      integer, dimension(*) :: array
      !local variables
      integer :: i,nstart_
      nstart_=max(nstart,0)
      do i=1,npaddim*ndebug
         array(nstart_+i)= int(r_nan())
      end do
    end subroutine i_padding
    !!***

    !!****f* m_profiling/l_padding
    !! NAME
    !! l_padding
    !!
    !! FUNCTION
    !!   Pad the end of the array with .false. values.
!! PARENTS
!!      m_profiling
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine l_padding(npaddim,nstart,array)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'l_padding'
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
    end subroutine l_padding
    !!***

    !!****f* m_profiling/c_padding
    !! NAME
    !! c_padding
    !!
    !! FUNCTION
    !!   Pad the end of the array with character values.
!! PARENTS
!!      m_profiling
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine c_padding(npaddim,nstart,array)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'c_padding'
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
    end subroutine c_padding
    !!***

    !!****f* m_profiling/mo_dp1
    !! NAME
    !! mo_dp1
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_dp1(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_dp1'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_dp1
    !!***

    !!****f* m_profiling/mo_dp2
    !! NAME
    !! mo_dp2
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_dp2(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_dp2'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_dp2
    !!***

    !!****f* m_profiling/mo_dp3
    !! NAME
    !! mo_dp3
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_dp3(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_dp3'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_dp3
    !!***

    !!****f* m_profiling/mo_dp4
    !! NAME
    !! mo_dp4
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_dp4(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_dp4'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_dp4
    !!***

    !!****f* m_profiling/mo_dp5
    !! NAME
    !! mo_dp5
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_dp5(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_dp5'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_dp5
    !!***

    !!****f* m_profiling/mo_dp6
    !! NAME
    !! mo_dp6
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_dp6(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_dp6'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_dp6
    !!***

    !!****f* m_profiling/mo_dp7
    !! NAME
    !! mo_dp7
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_dp7(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_dp7'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_dp7
    !!***

    !!****f* m_profiling/mo_sp1
    !! NAME
    !! mo_sp1
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_sp1(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_sp1'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_sp1
    !!***

    !!****f* m_profiling/mo_sp2
    !! NAME
    !! mo_sp2
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_sp2(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_sp2'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_sp2
    !!***

    !!****f* m_profiling/mo_sp3
    !! NAME
    !! mo_sp3
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_sp3(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_sp3'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_sp3
    !!***

    !!****f* m_profiling/mo_sp4
    !! NAME
    !! mo_sp4
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_sp4(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_sp4'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_sp4
    !!***

    !!****f* m_profiling/mo_sp5
    !! NAME
    !! mo_sp5
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_sp5(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_sp5'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_sp5
    !!***

    !!****f* m_profiling/mo_sp6
    !! NAME
    !! mo_sp6
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_sp6(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_sp6'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_sp6
    !!***

    !!****f* m_profiling/mo_sp7
    !! NAME
    !! mo_sp7
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_sp7(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_sp7'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_sp7
    !!***

    !!****f* m_profiling/mo_i1
    !! NAME
    !! mo_i1
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_i1(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_i1'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_i1
    !!***

    !!****f* m_profiling/mo_i2
    !! NAME
    !! mo_i2
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_i2(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_i2'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_i2
    !!***

    !!****f* m_profiling/mo_i3
    !! NAME
    !! mo_i3
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_i3(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_i3'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_i3
    !!***

    !!****f* m_profiling/mo_i4
    !! NAME
    !! mo_i4
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_i4(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_i4'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_i4
    !!***

    !!****f* m_profiling/mo_i5
    !! NAME
    !! mo_i5
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_i5(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_i5'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_i5
    !!***

    !!****f* m_profiling/mo_i6
    !! NAME
    !! mo_i6
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_i6(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_i6'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_i6
    !!***

    !!****f* m_profiling/mo_i7
    !! NAME
    !! mo_i7
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_i7(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_i7'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_i7
    !!***

    !!****f* m_profiling/mo_l1
    !! NAME
    !! mo_l1
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_l1(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_l1'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_l1
    !!***

    !!****f* m_profiling/mo_l2
    !! NAME
    !! mo_l2
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_l2(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_l2'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_l2
    !!***

    !!****f* m_profiling/mo_l3
    !! NAME
    !! mo_l3
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_l3(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_l3'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_l3
    !!***

    !!****f* m_profiling/mo_l4
    !! NAME
    !! mo_l4
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_l4(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_l4'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_l4
    !!***

    !!****f* m_profiling/mo_l5
    !! NAME
    !! mo_l5
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_l5(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_l5'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_l5
    !!***

    !!****f* m_profiling/mo_l6
    !! NAME
    !! mo_l6
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_l6(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_l6'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_l6
    !!***

    !!****f* m_profiling/mo_l7
    !! NAME
    !! mo_l7
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_l7(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_l7'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_l7
    !!***

    !!****f* m_profiling/mo_c1
    !! NAME
    !! mo_c1
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_c1(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_c1'
!End of the abilint section

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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_c1
    !!***

    !!****f* m_profiling/mo_cmpdp1
    !! NAME
    !! mo_cmpdp1
    !!
    !! FUNCTION
    !!   Pad the given array with debug values and call the statistic
    !!   on allocation size.
!! PARENTS
!!
!! CHILDREN
!!      memory_occupation
!!
    !! SOURCE
    !!
    subroutine mo_cmpdp1(istat,array,aname,rname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mo_cmpdp1'
!End of the abilint section

      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      complex(kind=8), dimension(:), intent(in) :: array
      !local variables
      integer :: ndim
      if (ndebug /=0) then
         ndim=product(shape(array))-ndebug
         !TODO: to be implemented.
         !call cmpdp_padding(1,ndim,array)
      end if
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_cmpdp1
    !!***

end module m_profiling
