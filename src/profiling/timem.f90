!the same timing routine but with system_clock (in case of a supported specs)
subroutine timing(iproc,category,action)
  implicit none
  include 'mpif.h'
  !Variables
  integer, intent(in) :: iproc
  character(len=*), intent(in) :: category
  character(len=2), intent(in) :: action      ! possibilities: INitialize, ON, OFf, REsults
  !Local variables
  logical :: parallel,init
  integer, parameter :: ncat=22   ! define timimg categories
  integer :: i,ierr,ii,i_all,i_stat,nproc
  integer :: istart,iend,count_rate,count_max,ielapsed,ncounters,itime,ittime
  !cputime routine gives a real
  !real :: total,total0,time,time0
  real(kind=8) :: pc,total_pc,total
  character(len=10), dimension(ncat) :: pcnames !names of the partial counters, to be assigned
  integer, dimension(ncat+1) :: itsum
  real(kind=8), dimension(ncat+1) :: timesum
  real(kind=8), dimension(ncat) :: pctimes !total times of the partial counters
  save :: init,itsum,istart,timesum,ittime,parallel,pcnames,pctimes,ncounters

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

  !first of all, read the time
  call system_clock(itime,count_rate,count_max)

  ! write(*,*) 'ACTION=',action,'...','CATEGORY=',category,'...'
  if (action.eq.'IN') then  ! INIT
     !!no need of using system clock for the total time (presumably more than a millisecond)
     !call cpu_time(total0)
     ittime=itime
     do i=1,ncat
        itsum(i)=0
        timesum(i)=0.d0
        pctimes(i)=0.d0
     enddo
     parallel=trim(category).eq.'parallel'
     init=.false.
     ncounters=0

  else if (action.eq.'PR') then !stop partial counters and restart from the beginning
     if (init.neqv..false.) then
        print *, 'TIMING IS INITIALIZED BEFORE PARTIAL RESULTS'
        stop 
     endif
     ncounters=ncounters+1
     if (ncounters > ncat) then
        print *, 'It is not allowed to have more partial counters that categories; ncat=',ncat
        stop
     end if
     pcnames(ncounters)=trim(category)
     if (itime-ittime < 0) then
        timesum(ncat+1)=real(itime-ittime+count_max,kind=8)/real(count_rate,kind=8)
     else
        timesum(ncat+1)=real(itime-ittime,kind=8)/real(count_rate,kind=8)
     end if
     pctimes(ncounters)=timesum(ncat+1)
     call sum_results(parallel,iproc,ncat,cats,itsum,timesum,pcnames(ncounters))
     !reset all timings
     ittime=itime
     do i=1,ncat
        itsum(i)=0
        timesum(i)=0.d0
     enddo


  else if (action.eq.'RE') then ! RESULT
     if (init.neqv..false.) then
        print *, 'TIMING IS INITIALIZED BEFORE RESULTS'
        stop 
     endif

     if (ncounters == 0) then !no partial counters selected
        if (itime-ittime < 0) then
           timesum(ncat+1)=real(itime-ittime+count_max,kind=8)/real(count_rate,kind=8)
        else
           timesum(ncat+1)=real(itime-ittime,kind=8)/real(count_rate,kind=8)
        end if

        call sum_results(parallel,iproc,ncat,cats,itsum,timesum,'ALL')
     else !consider only the results of the partial counters
        total_pc=0.d0
        do i=1,ncounters
           total_pc=total_pc+pctimes(i)
        end do
        if (parallel) then 
           !calculate total time
           call MPI_REDUCE(total_pc,total,1,&
                MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
           call MPI_REDUCE(pctimes,timesum,ncounters,&
                MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        else
           total=total_pc
           do i=1,ncounters
              timesum(i)=pctimes(i)
           end do
        end if
        if (iproc.eq.0) then
           call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
           write(60,*) 'PARTIAL COUNTER   mean TIME(sec)       PERCENT'
           total_pc=0.d0
           do i=1,ncounters
              pc=100.d0*timesum(i)/sum(timesum(1:ncounters))
              if (timesum(i) /= 0.d0) write(60,'(a14,1(10x,1pe9.2),5x,0pf8.3 )') pcnames(i),timesum(i)/real(nproc,kind=8),pc
              total_pc=total_pc+pc
           enddo
           write(60,'(70("-"))')
           write(60,'(a,10x,1pe9.2,6x,a,0pf5.1)') &
                'Total CPU time=',total,'Total categorized percent ',total_pc
        end if
     end if

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
        istart=itime
        init=.true.
     else if (action.eq.'OF') then  ! OFF
        if (init.neqv..true.) then
           print *, cats(ii), 'not initialized'
           stop 
        endif
        iend=itime
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

subroutine sum_results(parallel,iproc,ncat,cats,itsum,timesum,message)
  implicit none
  include 'mpif.h'
  character(len=*), intent(in) :: message
  logical, intent(in) :: parallel
  integer, intent(in) :: iproc,ncat
  !real, intent(in) :: total0
  character(len=14), dimension(ncat), intent(in) :: cats
  integer, dimension(ncat+1), intent(in) :: itsum
  real(kind=8), dimension(ncat+1), intent(inout) :: timesum
  !local variables
  integer :: i,iend,count_rate,count_max,ierr,nproc
  real :: total
  real(kind=8) :: total_pc,pc,totaltime,sum
  real(kind=8), dimension(ncat+1) :: timetot!timemax,timemin

  !time for other categories
  call system_clock(iend,count_rate,count_max)
  sum=0.d0
  do i=1,ncat
     timesum(i)=timesum(i)+real(itsum(i),kind=8)/real(count_rate,kind=8)
     sum=sum+timesum(i)
  end do


  !print *,'iproc, sum, tot',iproc,sum,timesum(ncat+1),count_max

  if (parallel) then 
     !old method, do not work due to the triangular inequality
     !call MPI_REDUCE(timesum,timemax,ncat+1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
     !call MPI_REDUCE(timesum,timemin,ncat+1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)

     !newmethod
     !calculate total time
     call MPI_REDUCE(timesum(ncat+1),totaltime,1,&
          MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
     call MPI_REDUCE(timesum,timetot,ncat+1,&
          MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  else
     totaltime=timesum(ncat+1)
     do i=1,ncat+1
        timetot(i)=timesum(i)
     end do
!!$     do i=1,ncat+1
!!$        timemax(i)=timesum(i)
!!$        timemin(i)=timesum(i)
!!$     enddo
  endif
  !total=real(timemax(ncat+1),kind=4)

  if (iproc.eq.0) then
     call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
     open(unit=60,file='time.prc',status='unknown')
     !write(60,*) 'CATEGORY          min. TIME(sec)     max. TIME(sec)           PERCENT'
     write(60,*) 'CATEGORY          mean TIME(sec)       PERCENT'
     total_pc=0.d0
     do i=1,ncat
        pc=100.d0*timetot(i)/timetot(ncat+1)!real(total,kind=8)
        if (timetot(i) /= 0.d0) write(60,'(a14,1(10x,1pe9.2),5x,0pf8.3 )') cats(i),timetot(i)/real(nproc,kind=8),pc
        total_pc=total_pc+pc
     enddo
     write(60,'(70("-"))')
     write(60,'(a,10x,1pe9.2,6x,a,0pf5.1)') &
          'Total CPU time for category: '//message//'=',totaltime,'Total categorized percent ',total_pc
  endif

end subroutine sum_results


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
