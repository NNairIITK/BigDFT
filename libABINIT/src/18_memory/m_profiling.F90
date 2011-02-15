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

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif
  
  module m_profiling
    implicit none

    private

    !Memory profiling
    type :: memstat
       character(len=36) :: routine,array
       integer(kind=8) :: memory,peak
    end type memstat

    ! Save values for memocc.
#if defined DEBUG_MODE
    integer, parameter :: ndebug = 5
    integer :: verbose = 2
#else
    integer, parameter :: ndebug = 0
    integer :: verbose = 0
#endif
    real :: memorylimit = 0.e0
    logical :: meminit = .false.
    type(memstat) :: memloc,memtot
    integer :: memalloc,memdealloc,memproc
    !Debug option for memocc, set in the input file
    logical :: memdebug

    !interface for the memory allocation control, depends on ndebug
    interface memocc
       module procedure mo_dp1,mo_dp2,mo_dp3,mo_dp4,mo_dp5,mo_dp6,mo_dp7,&
            mo_sp1,mo_sp2,mo_sp3,mo_sp4,mo_sp5,mo_sp6,mo_sp7,&
            mo_i1,mo_i2,mo_i3,mo_i4,mo_i5,mo_i6,mo_i7,&
            mo_l1,mo_l2,mo_l3,mo_l4,mo_l5,mo_l6,mo_l7,&
            mo_c1, &
            memocc_internal  !central routine to be used for deallocation
    end interface

    public :: ndebug
    public :: memocc
    public :: memocc_set_verbosity
    public :: memocc_set_debug
    public :: memocc_set_memory_limit

  contains

    subroutine memocc_set_verbosity(verb)
      integer, intent(in) :: verb

      verbose = verb
    end subroutine memocc_set_verbosity

    subroutine memocc_set_debug(debug)
      logical, intent(in) :: debug

      memdebug = debug
    end subroutine memocc_set_debug

    subroutine memocc_set_memory_limit(limit)
      real, intent(in) :: limit

      memorylimit = limit
    end subroutine memocc_set_memory_limit

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
    !! SOURCE
    !!
    subroutine memory_occupation(istat,isize,array,routine)

      implicit none

      ! Arguments
      integer, intent(in) :: istat,isize
      character(len=*), intent(in) :: array,routine

      ! Local variables
      integer :: ierr

      include 'mpif.h'

      !print *,memproc,array,routine

      ! Initialised first
      if (.not.meminit) then
         memtot%memory=int(0,kind=8)
         memtot%peak=int(0,kind=8)
         memalloc=0
         memdealloc=0

         memloc%routine='routine'
         memloc%array='array'
         memloc%memory=int(0,kind=8) !fake initialisation to print the first routine
         memloc%peak=int(0,kind=8)

         !Use MPI to have the mpi rank
         call MPI_INITIALIZED(memproc,ierr)
         if (memproc /= 0) then
            call MPI_COMM_RANK(MPI_COMM_WORLD,memproc,ierr)
         end if

         !open the writing file for the root process
         if (memproc == 0 .and. verbose >= 2) then
            open(unit=98,file='malloc.prc',status='unknown')
            if (memdebug) then
               write(98,'(a,t40,a,t70,4(1x,a12))')&
                    '(Data in KB) Routine','Array name    ',&
                    'Array size','Total Memory'
            else
               write(98,'(a,t40,a,t70,4(1x,a12))')&
                    '(Data in KB) Routine','Peak Array    ',&
                    'Routine Mem','Routine Peak','Memory Stat.','Memory Peak'
            end if
         end if
         meminit = .true.
      end if

      select case(array)
      case('count')
         if (trim(routine)=='stop' .and. memproc==0) then
            if (verbose >= 2) then
               write(98,'(a,t40,a,t70,4(1x,i12))')&
                    trim(memloc%routine),trim(memloc%array),&
                    memloc%memory/int(1024,kind=8),memloc%peak/int(1024,kind=8),&
                    memtot%memory/int(1024,kind=8),&
                    (memtot%peak+memloc%peak-memloc%memory)/int(1024,kind=8)
               close(unit=98)
            end if
            write(*,'(1x,a)')&
                 '-------------------------MEMORY CONSUMPTION REPORT-----------------------------'
            write(*,'(1x,2(i0,a,1x),i0)')&
                 memalloc,' allocations and',memdealloc,' deallocations, remaining memory(B):',&
                 memtot%memory
            write(*,'(1x,a,i0,a)') 'memory occupation peak: ',memtot%peak/int(1048576,kind=8),' MB'
            write(*,'(4(1x,a))') 'for the array ',trim(memtot%array),&
                 'in the routine',trim(memtot%routine)
            !here we can add a routine which open the malloc.prc file in case of some 
            !memory allocation problem, and which eliminates it for a successful run
            if (.not.memdebug .and. memalloc == memdealloc .and. memtot%memory==int(0,kind=8)) then
               !clean the malloc file
               if (memproc==0 .and. verbose >= 2) then
                  open(unit=98,file='malloc.prc',status='unknown',action='write')
                  write(unit=98,fmt='()',advance='no')
                  close(unit=98)
               end if
            else
               call memory_malloc_check(memdebug,memalloc,memdealloc)
            end if
         else if (trim(routine)/='stop') then
            write(*,*) "memocc: ",array," ",routine
            write(*,"(a,i0,a)") "Error[",memproc,"]: Use memocc and the word 'count' only with the word 'stop'."
            stop
         end if

      case default
         !control of the allocation/deallocation status
         if (istat/=0) then
            if (isize>=0) then
               write(*,*)' subroutine ',routine,': problem of allocation of array ',array,&
                    ', error code=',istat,' exiting...'
               if (memproc == 0) close(unit=98)
               call MPI_ABORT(MPI_COMM_WORLD,ierr)
            else if (isize<0) then
               write(*,*)' subroutine ',routine,': problem of deallocation of array ',array,&
                    ', error code=',istat,' exiting...'
               if (memproc == 0) close(unit=98)
               call MPI_ABORT(MPI_COMM_WORLD,ierr)
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
            write(*,'(1x,a,f7.3,2(a,i0),a)')&
                 'ERROR: Memory limit of ',memorylimit,&
                 ' GB reached for memproc ',memproc,' : total memory is ',memtot%memory,' B.'
            write(*,'(1x,2(a,i0))')&
                 '       this happened for array '//trim(memtot%array)//' in routine '//trim(memtot%routine)
               call MPI_ABORT(MPI_COMM_WORLD,ierr)
         end if

         select case(memproc)
         case (0)
            if (memdebug .and. verbose >= 2) then
               !to be used for inspecting an array which is not deallocated
               write(98,'(a,t40,a,t70,4(1x,i12))')trim(routine),trim(array),isize,memtot%memory
            else
               !Compact format
               if (trim(memloc%routine) /= routine) then
                  if (memloc%memory /= int(0,kind=8) .and. verbose >= 2) then
                     write(98,'(a,t40,a,t70,4(1x,i12))')&
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


    !!****f* ABINIT/memory_malloc_check
    !! FUNCTION
    !!   Check the malloc.prc file (verbose format)
    !! SOURCE
    !!
    subroutine memory_malloc_check(memdebug,nalloc,ndealloc)


      implicit none
      !Arguments
      logical, intent(in) :: memdebug
      integer, intent(in) :: nalloc,ndealloc
      !Local variables
      if (memdebug .and. nalloc /= ndealloc) then
         write(*,*) &
              "Use the python script 'memcheck.py' in utils/scripts to check 'malloc.prc' file"
      end if
    END SUBROUTINE memory_malloc_check
    !!***


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

    !routine used for deallocations
    subroutine memocc_internal(istat,isize,array,routine)
      implicit none
      character(len=*), intent(in) :: array,routine
      integer, intent(in) :: istat,isize
      call memory_occupation(istat,isize,array,routine) !this routine is in profiling/memory.f90
    end subroutine memocc_internal

    subroutine dp_padding(npaddim,nstart,array)
      implicit none
      integer, intent(in) :: npaddim,nstart
      double precision, dimension(*) :: array
      !local variables
      integer :: i
      do i=1,npaddim*ndebug
         array(nstart+i)= d_nan() !this function is in profiling/memory.f90
      end do
    end subroutine dp_padding

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
         array(nstart+i)= r_nan() !this function is in profiling/timem.f90
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
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
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_c1
end module m_profiling
