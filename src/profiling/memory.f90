!!****f* BigDFT/memory_occupation
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
!!   The counters are initialized with:
!!      call memocc(0,iproc,'count','start') (iproc = mpi rank, nproc=mpi size)
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
!!    Copyright (C) Luigi Genovese, CEA Grenoble, France, 2007-2010
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! SOURCE
!!
subroutine memory_occupation(istat,isize,array,routine)

  use module_base, only: memorylimit,verbose,memstat,memloc,memtot,&
       & memalloc,memdealloc,memproc,meminit

  implicit none

! Arguments
  integer, intent(in) :: istat,isize
  character(len=*), intent(in) :: array,routine

! Local variables
  logical, parameter :: memdebug=.false.

  include 'mpif.h'

  integer :: ierr

  !print *,memproc,array,routine

  ! Don't use memory profiling if not initialised.
  if (.not.meminit .and. (array /= 'count' .or. routine /= 'start')) return

  select case(array)
  case('count')
     if (routine=='start') then
        memtot%memory=int(0,kind=8)
        memtot%peak=int(0,kind=8)
        memalloc=0
        memdealloc=0

        memloc%routine='routine'
        memloc%array='array'
        memloc%memory=int(0,kind=8) !fake initialisation to print the first routine
        memloc%peak=int(0,kind=8)

        memproc=isize
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
     else if (routine=='stop' .and. memproc==0) then
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
        write(*,'(1x,2(i0,a),i0)')&
             memalloc,' allocations and ',memdealloc,' deallocations, remaining memory(B):',&
             memtot%memory
        write(*,'(1x,a,i0,a)') 'memory occupation peak: ',memtot%peak/int(1048576,kind=8),' MB'
        write(*,'(4(1x,a))') 'for the array ',trim(memtot%array),&
             'in the routine',trim(memtot%routine)
        !here we can add a routine which open the malloc.prc file in case of some 
        !memory allocation problem, and which eliminates it for a successful run
        if (.not.memdebug .and. memalloc == memdealloc .and. memtot%memory==int(0,kind=8)) then
           !clean the malloc file
           if (memproc==0 .and. verbose >= 2) then
              open(unit=98,file='malloc.prc',status='unknown')
              write(98,*)
              close(unit=98)
           end if
        else
           call memory_malloc_check(memdebug,memalloc,memdealloc)
        end if
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
end subroutine memory_occupation
!!***


!!****f* BigDFT/memory_malloc_check
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
end subroutine memory_malloc_check
!!***


!!****f* BigDFT/d_nan
!! FUNCTION
!!   Function which specify NaN according to IEEE specifications
!! SOURCE
!!
function d_nan()
  implicit none
  real(kind=8) :: d_nan
  !local variables
  real(kind=8) :: dnan
  integer, dimension(2) :: inan
  equivalence (dnan, inan)
  ! This first assignment is for big-endian machines
  inan(1) = 2147483647
  ! The second assignment is for little-endian machines
  inan(2) = 2147483647
  d_nan = dnan
end function d_nan
!!***

!!****f* BigDFT/r_nan
!! FUNCTION
!!   Function which specify NaN according to IEEE specifications
!! SOURCE
!!
function r_nan()
  implicit none
  real(kind=4) :: r_nan
  !local variables
  real(kind=4) :: rnan
  integer :: inan
  equivalence (rnan, inan)
  inan = 2147483647
  r_nan = rnan
end function r_nan
!!***
