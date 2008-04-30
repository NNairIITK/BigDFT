!!****m* BigDFT/base
!! NAME
!!   base
!!
!! FUNCTION
!!  Modules which contains the low level definitions, as well as some profiling procedures
!!
!! DESCRIPTION
!!  Interfaces of:
!!
!! AUTHOR
!!    Luigi Genovese
!!
!! COPYRIGHT
!!    Copyright (C) 2008 CEA
!!
!! SOURCE
!!
module module_base

  implicit none

  interface memocc
     module procedure mo_dp1,mo_dp2,mo_dp3,mo_dp4,mo_dp5,mo_dp6,mo_dp7,&
          mo_sp1,mo_sp2,mo_sp3,mo_sp4,mo_sp5,mo_sp6,mo_sp7,&
          mo_i1,mo_i2,mo_i3,mo_i4,mo_i5,mo_i6,mo_i7,&
          memocc_internal  !central routine to be used for deallocation
  end interface

  contains
    !control the memory occupation by calculating the overall size of the allocated arrays
    !usage: 
    ! when allocating allocating an array "stuff" of dimension n in the routine "dosome"
    !  allocate(stuff(n),stat=i_stat)
    !  call memocc(i_stat,product(shape(stuff))*kind(stuff),'stuff','dosome')
    ! when deallocating 
    !  i_all=-product(shape(stuff))*kind(stuff)
    !  deallocate(stuff,stat=i_stat)
    !  call memocc(i_stat,i_all,'stuff','dosome')
    ! the counters are initialized with
    !  call memocc(0,iproc,'count','start') (iproc = mpi rank, nproc=mpi size)
    ! and stopped with
    !  call memocc(0,0,'count','stop')
    ! at the end of the calculation a short report is printed on the screen
    ! some information can be also written on disk following the needs
    !  This file is distributed under the terms of the
    !  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
    !  Copyright (C) Luigi Genovese, CEA Grenoble, France, 2007
    subroutine memory_occupation(istat,isize,array,routine)
      implicit none
      character(len=*), intent(in) :: array,routine
      integer, intent(in) :: istat,isize
      !local variables
      character(len=36) :: maxroutine,locroutine
      character(len=36) :: maxarray,locarray
      integer :: nalloc,ndealloc,locpeak,locmemory,iproc
      integer(kind=8) :: memory,maxmemory
      save :: memory,nalloc,ndealloc,maxroutine,maxarray,maxmemory
      save :: locroutine,locarray,locpeak,locmemory,iproc
      select case(array)
      case('count')
         if (routine=='start') then
            memory=int(0,kind=8)
            maxmemory=int(0,kind=8)
            nalloc=0
            ndealloc=0
            locroutine='routine'
            locarray='array'
            locmemory=0
            locpeak=0
            iproc=isize
            !open the writing file for the root process
            if (iproc == 0) then
               open(unit=98,file='malloc.prc',status='unknown')
               write(98,'(a32,a14,4(1x,a12))')&
                    '(Data in KB)             Routine','    Peak Array',&
                    'Routine Mem','Routine Peak','Memory Stat.','Memory Peak'
            end if
         else if (routine=='stop' .and. iproc==0) then
            write(98,'(a32,a14,4(1x,i12))')&
                 trim(locroutine),trim(locarray),&
                 locmemory/1024,locpeak/1024,memory/int(1024,kind=8),&
                 (memory+int(locpeak-locmemory,kind=8))/int(1024,kind=8)
            close(98)
            write(*,'(1x,a)')&
                 '-------------------------MEMORY CONSUMPTION REPORT-----------------------------'
            write(*,'(1x,2(i0,a),i0)')&
                 nalloc,' allocations and ',ndealloc,' deallocations, remaining memory(B):',memory
            write(*,'(1x,a,i0,a)') 'memory occupation peak: ',maxmemory/int(1048576,kind=8),' MB'
            write(*,'(4(1x,a))') 'for the array ',trim(maxarray),'in the routine',trim(maxroutine)
         end if

      case default
         !control of the allocation/deallocation status
         if (istat/=0) then
            if (isize>=0) then
               write(*,*)' subroutine ',routine,': problem of allocation of array ',array,&
                    'error code=',istat,' exiting...'
               stop
            else if (isize<0) then
               write(*,*)' subroutine ',routine,': problem of deallocation of array ',array,&
                    'error code=',istat,' exiting...'
               stop
            end if
         end if
         select case(iproc)
         case (0)
            !to be used for inspecting an array which is not deallocated
            !write(98,'(a32,a14,4(1x,i12))')trim(routine),trim(array),isize,memory
            if (trim(locroutine) /= routine) then
               write(98,'(a32,a14,4(1x,i12))')&
                    trim(locroutine),trim(locarray),&
                    locmemory/1024,locpeak/1024,memory/int(1024,kind=8),&
                    (memory+int(locpeak-locmemory,kind=8))/int(1024,kind=8)
               locroutine=routine
               locarray=array
               locmemory=isize
               locpeak=isize
            else
               locmemory=locmemory+isize
               if (locmemory > locpeak) then
                  locpeak=locmemory
                  locarray=array
               end if
            end if
            memory=memory+int(isize,kind=8)
            if (memory > maxmemory) then
               maxmemory=memory
               maxroutine=routine
               maxarray=array
            end if
            if (isize>0) then
               nalloc=nalloc+1
            else if (isize<0) then
               ndealloc=ndealloc+1
            end if
         case default
            return
         end select
      end select
    end subroutine memory_occupation

    subroutine memocc_internal(istat,isize,array,routine)
      implicit none
      character(len=*), intent(in) :: array,routine
      integer, intent(in) :: istat,isize
      call memory_occupation(istat,isize,array,routine)
    end subroutine memocc_internal

    !beginning of the verbose section
    subroutine mo_dp1(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=8), dimension(:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_dp1

    subroutine mo_dp2(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=8), dimension(:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_dp2

    subroutine mo_dp3(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=8), dimension(:,:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_dp3

    subroutine mo_dp4(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=8), dimension(:,:,:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_dp4

    subroutine mo_dp5(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=8), dimension(:,:,:,:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_dp5

    subroutine mo_dp6(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=8), dimension(:,:,:,:,:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_dp6

    subroutine mo_dp7(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=8), dimension(:,:,:,:,:,:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_dp7

    subroutine mo_sp1(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=4), dimension(:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_sp1

    subroutine mo_sp2(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=4), dimension(:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_sp2

    subroutine mo_sp3(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=4), dimension(:,:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_sp3

    subroutine mo_sp4(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=4), dimension(:,:,:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_sp4

    subroutine mo_sp5(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=4), dimension(:,:,:,:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_sp5

    subroutine mo_sp6(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=4), dimension(:,:,:,:,:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_sp6

    subroutine mo_sp7(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      real(kind=4), dimension(:,:,:,:,:,:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_sp7

    subroutine mo_i1(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_i1

    subroutine mo_i2(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_i2

    subroutine mo_i3(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_i3

    subroutine mo_i4(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:,:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_i4

    subroutine mo_i5(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:,:,:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_i5

    subroutine mo_i6(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:,:,:,:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_i6

    subroutine mo_i7(istat,array,aname,rname)
      implicit none
      character(len=*), intent(in) :: aname,rname
      integer, intent(in) :: istat
      integer, dimension(:,:,:,:,:,:,:), intent(in) :: array
      call memory_occupation(istat,product(shape(array))*kind(array),aname,rname)
    end subroutine mo_i7



end module module_base
!!***
