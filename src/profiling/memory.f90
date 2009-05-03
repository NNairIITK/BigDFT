!!****f* BigDFT/memory_occupation
!! FUNCTION
!! Control the memory occupation by calculating the overall size of the allocated arrays
!! DESCRIPTION
!!   when allocating allocating an array "stuff" of dimension n in the routine "dosome"
!!     allocate(stuff(n),stat=i_stat)
!!     call memocc(i_stat,product(shape(stuff))*kind(stuff),'stuff','dosome')
!!   when deallocating 
!!     i_all=-product(shape(stuff))*kind(stuff)
!!     deallocate(stuff,stat=i_stat)
!!     call memocc(i_stat,i_all,'stuff','dosome')
!!   the counters are initialized with
!!     call memocc(0,iproc,'count','start') (iproc = mpi rank, nproc=mpi size)
!!   and stopped with
!!     call memocc(0,0,'count','stop')
!!   at the end of the calculation a short report is printed on the screen
!!   some information can be also written on disk following the needs
!!
!!   The file malloc.prc is not deleted if the final total memory is not equal
!!   to zero.
!!   memdebug (parameter)
!!     == .true.  then display a line per allocation or deallocation
!!                a routine at the end parses the file
!!     == .false. compact format
!! COPYRIGHT
!!    Copyright (C) Luigi Genovese, CEA Grenoble, France, 2007-2009
!!    This file is distributed under the terms of the
!!    GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!! SOURCE
!!
subroutine memory_occupation(istat,isize,array,routine)

  use module_base, only: memorylimit

  implicit none

! Arguments
  integer, intent(in) :: istat,isize
  character(len=*), intent(in) :: array,routine

! Local variables
  logical, parameter :: memdebug=.true.
  type :: memstat
     character(len=36) :: routine,array
     integer(kind=8) :: memory,peak
  end type memstat

  include 'mpif.h'

  type(memstat), save :: loc,tot
  integer, save :: nalloc,ndealloc,iproc
  integer :: ierr

  !print *,iproc,array,routine

  select case(array)
  case('count')
     if (routine=='start') then
        tot%memory=int(0,kind=8)
        tot%peak=int(0,kind=8)
        nalloc=0
        ndealloc=0

        loc%routine='routine'
        loc%array='array'
        loc%memory=int(0,kind=8) !fake initialisation to print the first routine
        loc%peak=int(0,kind=8)

        iproc=isize
        !open the writing file for the root process
        if (iproc == 0) then
           open(unit=98,file='malloc.prc',status='unknown')
           if (memdebug) then
              write(98,'(a,t40,a,t60,4(1x,a12))')&
                   '(Data in KB)             Routine','    Array name',&
                   'Array size','Total Memory'
           else
              write(98,'(a,t40,a,t60,4(1x,a12))')&
                   '(Data in KB)             Routine','    Peak Array',&
                   'Routine Mem','Routine Peak','Memory Stat.','Memory Peak'
           end if
        end if
     else if (routine=='stop' .and. iproc==0) then
        write(98,'(a,t40,a,t60,4(1x,i12))')&
             trim(loc%routine),trim(loc%array),&
             loc%memory/int(1024,kind=8),loc%peak/int(1024,kind=8),&
             tot%memory/int(1024,kind=8),&
             (tot%peak+loc%peak-loc%memory)/int(1024,kind=8)
        close(unit=98)
        write(*,'(1x,a)')&
             '-------------------------MEMORY CONSUMPTION REPORT-----------------------------'
        write(*,'(1x,2(i0,a),i0)')&
             nalloc,' allocations and ',ndealloc,' deallocations, remaining memory(B):',&
             tot%memory
        write(*,'(1x,a,i0,a)') 'memory occupation peak: ',tot%peak/int(1048576,kind=8),' MB'
        write(*,'(4(1x,a))') 'for the array ',trim(tot%array),&
             'in the routine',trim(tot%routine)
        !here we can add a routine which open the malloc.prc file in case of some 
        !memory allocation problem, and which eliminates it for a successful run
        if (.not.memdebug .and. nalloc == ndealloc .and. tot%memory==int(0,kind=8)) then
           !clean the malloc file
           if (iproc==0) then
              open(unit=98,file='malloc.prc',status='unknown')
              write(98,*)
              close(unit=98)
           end if
        else
           call memory_malloc_check(memdebug,nalloc,ndealloc)
        end if
     end if

  case default
     !control of the allocation/deallocation status
     if (istat/=0) then
        if (isize>=0) then
           write(*,*)' subroutine ',routine,': problem of allocation of array ',array,&
                ', error code=',istat,' exiting...'
           if (iproc == 0) close(unit=98)
           stop
        else if (isize<0) then
           write(*,*)' subroutine ',routine,': problem of deallocation of array ',array,&
                ', error code=',istat,' exiting...'
           if (iproc == 0) close(unit=98)
           stop
        end if
     end if
     !total counter, for all the processes
     tot%memory=tot%memory+int(isize,kind=8)
     if (tot%memory > tot%peak) then
        tot%peak=tot%memory
        tot%routine=routine
        tot%array=array
     end if
     if (isize>0) then
        nalloc=nalloc+1
     else if (isize<0) then
        ndealloc=ndealloc+1
     end if

     if (memorylimit /= 0.e0 .and. &
          tot%memory > int(real(memorylimit,kind=8)*1073741824.d0,kind=8)) then !memory limit is in GB
        write(*,'(1x,a,f7.3,2(a,i0),a)')&
             'ERROR: Memory limit of ',memorylimit,&
             ' GB reached for iproc ',iproc,' : total memory is ',tot%memory,' B.'
        write(*,'(1x,2(a,i0))')&
             '       this happened for array '//trim(tot%array)//' in routine '//trim(tot%routine)
        call MPI_ABORT(MPI_COMM_WORLD,ierr)
     end if

     select case(iproc)
     case (0)
        if (memdebug) then
           !to be used for inspecting an array which is not deallocated
           write(98,'(a,t40,a,t60,4(1x,i12))')trim(routine),trim(array),isize,tot%memory
        else
           !Compact format
           if (trim(loc%routine) /= routine) then
              if (loc%memory /= int(0,kind=8)) then
                 write(98,'(a,t40,a,t60,4(1x,i12))')&
                      trim(loc%routine),trim(loc%array),&
                      loc%memory/int(1024,kind=8),loc%peak/int(1024,kind=8),&
                      tot%memory/int(1024,kind=8),&
                      (tot%memory+loc%peak-loc%memory)/int(1024,kind=8)
              end if
              loc%routine=routine
              loc%array=array
              loc%memory=isize
              loc%peak=isize
           else
              loc%memory=loc%memory+isize
              if (loc%memory > loc%peak) then
                 loc%peak=loc%memory
                 loc%array=array
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
       "Use the python script 'malloc_check.py' in config/scripts to check 'malloc.prc' file"
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
