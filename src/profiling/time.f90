module timeData
  integer, parameter :: ncat=23   ! define timimg categories

  integer :: istart, ittime, ncounters, ncaton!, nskip
  logical :: parallel,init
  integer, dimension(ncat+1) :: itsum
  real(kind=8), dimension(ncat+1) :: timesum
  real(kind=8), dimension(ncat) :: pctimes !total times of the partial counters
  character(len=10), dimension(ncat) :: pcnames !names of the partial counters, to be assigned
end module timeData

!the same timing routine but with system_clock (in case of a supported specs)
subroutine timing(iproc,category,action)
  use timeData

  implicit none

  include 'mpif.h'
  !Variables
  integer, intent(in) :: iproc
  character(len=*), intent(in) :: category
  character(len=2), intent(in) :: action      ! possibilities: INitialize, ON, OFf, REsults
  !Local variables
  integer :: i,ierr,ii,nproc
  integer :: iend,count_rate,count_max,ielapsed,itime
  !cputime routine gives a real
  !real :: total,total0,time,time0
  real(kind=8) :: pc,total_pc,total

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
       'Tail          '    ,  &
       'Davidson      '    /)    !

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
        if (iproc == 0) then
           if (parallel) then
              call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
           else
              nproc=1
           end if
           write(60,*)
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
           write(60,*)
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
     !print *,'find category',ii,trim(category)
     if (ii.eq.100000) then
        print *, 'ACTION  ',action
        stop 'TIMING CATEGORY NOT DEFINED'
     end if

     if (action.eq.'ON') then  ! ON
        if (init.neqv..false.) then
!!!           print *, cats(ii),': TIMING INITIALIZED BEFORE READ'
!!!           stop 
           !some other category was initalized before, taking that one
           return
        endif
        istart=itime
        init=.true.
        ncaton=ii
     else if (action.eq.'OF' .and. ii==ncaton) then  ! OFF
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
        print *,action,ii,ncaton,trim(category)
        stop 'TIMING ACTION UNDEFINED'
     endif

  endif

END SUBROUTINE timing

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
!!!     do i=1,ncat+1
!!!        timemax(i)=timesum(i)
!!!        timemin(i)=timesum(i)
!!!     enddo
  endif
  !total=real(timemax(ncat+1),kind=4)

  if (iproc.eq.0) then
     if (parallel) then
        call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
     else
        nproc=1
     end if
     open(unit=60,file='time.prc',status='unknown')
     write(60,*)
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
     write(60,*)
  endif

END SUBROUTINE sum_results
