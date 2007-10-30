!!$subroutine timing(iproc,category,action)
!!$  implicit none
!!$  include 'mpif.h'
!!$  !Variables
!!$  integer, intent(in) :: iproc
!!$  character(len=14), intent(in) :: category
!!$  character(len=2), intent(in) :: action      ! possibilities: INitialize, ON, OFf, REsults
!!$  !Local variables
!!$  logical :: parallel,init
!!$  integer, parameter :: ncat=22   ! define timimg categories
!!$  integer :: i,ierr,ii
!!$  !cputime routine gives a real
!!$  real :: total,total0,time,time0
!!$  real(kind=8) :: pc,total_pc
!!$  real(kind=8) :: flops(ncat),timesum(ncat+1),timemax(ncat+1),timemin(ncat+1)
!!$  save :: time0,init,timesum,total0,parallel
!!$
!!$  character(len=14), dimension(ncat), parameter :: cats = (/ &
!!$       'ReformatWaves '    ,  &  !  Reformatting of input waves
!!$       'CrtDescriptors'    ,  &  !  Calculation of descriptor arrays
!!$       'CrtLocPot     '    ,  &  !  Calculation of local potential
!!$       'CrtProjectors '    ,  &  !  Calculation of projectors
!!$       'ApplyLocPotKin'    ,  &  !  Application of PSP, kinetic energy
!!$       'ApplyProj     '    ,  &  !  Application of nonlocal PSP
!!$       'Precondition  '    ,  &  !  Precondtioning
!!$       'Rho_comput    '    ,  &  !  Calculation of charge density (sumrho) computation
!!$       'Rho_commun    '    ,  &  !  Calculation of charge density (sumrho) communication
!!$       'Un-TransSwitch'    ,  &  !  Transposition of wavefunction, computation
!!$       'Un-TransComm  '    ,  &  !  Transposition of wavefunction, communication
!!$       'GramS_comput  '    ,  &  !  Gram Schmidt computation        
!!$       'GramS_commun  '    ,  &  !  Gram Schmidt communication
!!$       'LagrM_comput  '    ,  &  !  Lagrange Multipliers computation
!!$       'LagrM_commun  '    ,  &  !  Lagrange Multipliers communication
!!$       'Diis          '    ,  &  !  Wavefunction mixing time, steepest descent if present
!!$       'PSolv_comput  '    ,  &  !
!!$       'PSolv_commun  '    ,  &  !
!!$       'PSolvKernel   '    ,  &  !
!!$       'Exchangecorr  '    ,  &  !
!!$       'Forces        '    ,  &  !
!!$       'Tail          '    /)    !
!!$
!!$  ! write(*,*) 'ACTION=',action,'...','CATEGORY=',category,'...'
!!$  if (action.eq.'IN') then  ! INIT
!!$     call cpu_time(total0)
!!$     do i=1,ncat
!!$        flops(i)=0.d0
!!$        timesum(i)=0.d0
!!$     enddo
!!$     parallel=trim(category).eq.'parallel'
!!$     init=.false.
!!$
!!$  else if (action.eq.'RE') then ! RESULT
!!$     if (init.neqv..false.) then
!!$        print *, 'TIMING INITIALIZED BEFORE RESULTS'
!!$        stop 
!!$     endif
!!$     !   sum results over all processor
!!$     call cpu_time(total)
!!$     total=total-total0
!!$     timesum(ncat+1)=total
!!$     if (parallel) then 
!!$        call MPI_REDUCE(timesum,timemax,ncat+1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
!!$        call MPI_REDUCE(timesum,timemin,ncat+1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
!!$     else
!!$        do i=1,ncat+1
!!$           timemax(i)=timesum(i)
!!$           timemin(i)=timesum(i)
!!$        enddo
!!$     endif
!!$     total=timemax(ncat+1)
!!$
!!$     if (iproc.eq.0) then
!!$        open(unit=60,file='time.prc',status='unknown')
!!$        write(60,*) 'CATEGORY          min. TIME(sec)     max. TIME(sec)           PERCENT'
!!$        total_pc=0.d0
!!$        do i=1,ncat
!!$           pc=100.d0*timemax(i)/real(total,kind=8)
!!$           write(60,'(a14,2(10x,1pe9.2),10x,0pf8.1 )') cats(i),timemin(i),timemax(i),pc
!!$           total_pc=total_pc+pc
!!$        enddo
!!$        write(60,'(70("-"))')
!!$        write(60,'(a,10x,1pe9.2,6x,a,0pf5.1)') 'Total CPU time',total,'Total categorized percent ',total_pc
!!$     endif
!!$
!!$  else
!!$
!!$     ii=100000
!!$     do i=1,ncat
!!$        if (trim(category).eq.trim(cats(i))) then
!!$           ii=i
!!$           exit
!!$        endif
!!$     enddo
!!$     if (ii.eq.100000) then
!!$        print *, 'ACTION  ',action
!!$        stop 'TIMING CATEGORY NOT DEFINED'
!!$     end if
!!$     !   write(*,*) 'category found',ii,cats(ii)
!!$
!!$     if (action.eq.'ON') then  ! ON
!!$        if (init.neqv..false.) then
!!$           print *, cats(ii),': TIMING INITIALIZED BEFORE READ'
!!$           stop 
!!$        endif
!!$        call cpu_time(time0)
!!$        init=.true.
!!$
!!$     else if (action.eq.'OF') then  ! OFF
!!$        if (init.neqv..true.) then
!!$           print *, cats(ii), 'not initialized'
!!$           stop 
!!$        endif
!!$        call cpu_time(time)
!!$        timesum(ii)=timesum(ii)+time-time0
!!$        init=.false.
!!$     else
!!$        stop 'TIMING ACTION UNDEFINED'
!!$     endif
!!$
!!$  endif
!!$
!!$end subroutine timing


!the same timing routine but with system_clock (in case of a supported specs)
subroutine timing(iproc,category,action)
  implicit none
  include 'mpif.h'
  !Variables
  integer, intent(in) :: iproc
  character(len=14), intent(in) :: category
  character(len=2), intent(in) :: action      ! possibilities: INitialize, ON, OFf, REsults
  !Local variables
  logical :: parallel,init
  integer, parameter :: ncat=22   ! define timimg categories
  integer :: i,ierr,ii,i_all,i_stat
  integer :: istart,iend,count_rate,count_max,ielapsed
  !cputime routine gives a real
  real :: total,total0,time,time0
  real(kind=8) :: pc,total_pc
  integer, dimension(ncat+1) :: itsum
  real(kind=8), dimension(ncat+1) :: timesum
  real(kind=8), dimension(:), allocatable :: timemax,timemin
  save :: init,itsum,istart,timesum,total0,parallel

  character(len=14), dimension(ncat), parameter :: cats = (/ &
       'ReformatWaves '    ,  &  !  Reformatting of input waves
       'CrtDescriptors'    ,  &  !  Calculation of descriptor arrays
       'CrtLocPot     '    ,  &  !  Calculation of local potential
       'CrtProjectors '    ,  &  !  Calculation of projectors
       'ApplyLocPotKin'    ,  &  !  Application of PSP, kinetic energy
       'ApplyProj     '    ,  &  !  Application of nonlocal PSP
       'Precondition  '    ,  &  !  Precondtioning
       'Rho_comput    '    ,  &  !  Calculation of charge density (sumrho) computation
       'Rho_commun    '    ,  &  !  Calculation of charge density (sumrho) communication
       'Un-TransSwitch'    ,  &  !  Transposition of wavefunction, computation
       'Un-TransComm  '    ,  &  !  Transposition of wavefunction, communication
       'GramS_comput  '    ,  &  !  Gram Schmidt computation        
       'GramS_commun  '    ,  &  !  Gram Schmidt communication
       'LagrM_comput  '    ,  &  !  Lagrange Multipliers computation
       'LagrM_commun  '    ,  &  !  Lagrange Multipliers communication
       'Diis          '    ,  &  !
       'PSolv_comput  '    ,  &  !
       'PSolv_commun  '    ,  &  !
       'PSolvKernel   '    ,  &  !
       'Exchangecorr  '    ,  &  !
       'Forces        '    ,  &  !
       'Tail          '    /)    !

  ! write(*,*) 'ACTION=',action,'...','CATEGORY=',category,'...'
  if (action.eq.'IN') then  ! INIT
     !no need of using system clock for the total time (presumably more than a millisecond)
     call cpu_time(total0) 
     do i=1,ncat
        itsum(i)=0
        timesum(i)=0.d0
     enddo
     parallel=trim(category).eq.'parallel'
     init=.false.

  else if (action.eq.'RE') then ! RESULT
     if (init.neqv..false.) then
        print *, 'TIMING INITIALIZED BEFORE RESULTS'
        stop 
     endif
     !   sum results over all processor
     !total elapsed time
     call cpu_time(total)
     total=total-total0
     timesum(ncat+1)=real(total,kind=8)

     !time for other categories
     call system_clock(iend,count_rate,count_max)
     do i=1,ncat
        timesum(i)=timesum(i)+real(itsum(i),kind=8)/real(count_rate,kind=8)
     end do

     allocate(timemax(ncat+1),stat=i_stat)
     call memocc(i_stat,product(shape(timemax))*kind(timemax),'timemax','timing')
     allocate(timemin(ncat+1),stat=i_stat)
     call memocc(i_stat,product(shape(timemin))*kind(timemin),'timemin','timing')

     if (parallel) then 
        call MPI_REDUCE(timesum,timemax,ncat+1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(timesum,timemin,ncat+1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
     else
        do i=1,ncat+1
           timemax(i)=timesum(i)
           timemin(i)=timesum(i)
        enddo
     endif
     total=real(timemax(ncat+1),kind=4)

     if (iproc.eq.0) then
        open(unit=60,file='time.prc',status='unknown')
        write(60,*) 'CATEGORY          min. TIME(sec)     max. TIME(sec)           PERCENT'
        total_pc=0.d0
        do i=1,ncat
           pc=100.d0*timemax(i)/real(total,kind=8)
           write(60,'(a14,2(10x,1pe9.2),10x,0pf8.1 )') cats(i),timemin(i),timemax(i),pc
           total_pc=total_pc+pc
        enddo
        write(60,'(70("-"))')
        write(60,'(a,10x,1pe9.2,6x,a,0pf5.1)') &
             'Total CPU time',total,'Total categorized percent ',total_pc
     endif

     i_all=-product(shape(timemax))*kind(timemax)
     deallocate(timemax,stat=i_stat)
     call memocc(i_stat,i_all,'timemax','timing')
     i_all=-product(shape(timemin))*kind(timemin)
     deallocate(timemin,stat=i_stat)
     call memocc(i_stat,i_all,'timemin','timing')


  else

     !controls if the category exists
     ii=100000
     do i=1,ncat
        if (trim(category).eq.trim(cats(i))) then
           ii=i
           exit
        endif
     enddo
     if (ii.eq.100000) then
        print *, 'ACTION  ',action
        stop 'TIMING CATEGORY NOT DEFINED'
     end if

   

     if (action.eq.'ON') then  ! ON
        if (init.neqv..false.) then
           print *, cats(ii),': TIMING INITIALIZED BEFORE READ'
           stop 
        endif
        call system_clock(istart)
        init=.true.
     else if (action.eq.'OF') then  ! OFF
        if (init.neqv..true.) then
           print *, cats(ii), 'not initialized'
           stop 
        endif
        call system_clock(iend,count_rate,count_max)
        if (iend-istart < 0) then
           ielapsed=iend-istart+count_max
        else
           ielapsed=iend-istart
        end if
        itsum(ii)=itsum(ii)+ielapsed
        !resum the values of the category inside timesum if too high
        if (itsum(ii) > count_max/2) then
           timesum(ii)=timesum(ii)+real(itsum(ii),kind=8)/real(count_rate,kind=8)
           itsum(ii)=0
        end if
        init=.false.
     else
        stop 'TIMING ACTION UNDEFINED'
     endif

  endif

end subroutine timing


!control the memory occupation by calculating the overall size in bytes of the allocated arrays
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
subroutine memocc(istat,isize,array,routine)
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
              write(*,*)' subroutine ',routine,': problem of allocation of array ',array
              stop
           else if (isize<0) then
              write(*,*)' subroutine ',routine,': problem of deallocation of array ',array
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

end subroutine memocc
