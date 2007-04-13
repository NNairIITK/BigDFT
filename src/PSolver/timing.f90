subroutine timing(iproc,category,action)
  implicit none
  include 'mpif.h'
  !Variables
  integer, intent(in) :: iproc
  character(len=14), intent(in) :: category
  character(len=2), intent(in) :: action      ! possibilities: INitialize, ON, OFf, REsults
  !Local variables
  logical :: parallel,init
  integer, parameter :: ncat=4   ! define timimg categories
  integer :: i,ierr,ii
  !cputime routine gives a real
  real :: total,total0,time,time0
  real(kind=8) :: pc,total_pc
  real(kind=8) :: flops(ncat),timesum(ncat+1),timemax(ncat+1),timemin(ncat+1)
  character(len=14) :: filename
  save :: time0,init,timesum,total0,parallel

  character(len=14), dimension(ncat), parameter :: cats = (/ &
              'PSolv_comput  '    ,  &  !
              'PSolv_commun  '    ,  &  !
              'PSolvKernel   '    ,  &  !
              'Exchangecorr  '    /)    !


! write(*,*) 'ACTION=',action,'...','CATEGORY=',category,'...'
  if (action.eq.'IN') then  ! INIT
    call cpu_time(total0)
    do i=1,ncat
      flops(i)=0.d0
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
    call cpu_time(total)
    total=total-total0
    timesum(ncat+1)=total
    if (parallel) then 
      call MPI_REDUCE(timesum,timemax,ncat+1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(timesum,timemin,ncat+1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
    else
      do i=1,ncat+1
        timemax(i)=timesum(i)
        timemin(i)=timesum(i)
      enddo
    endif
    total=timemax(ncat+1)

    if (iproc.eq.0) then

       if (parallel) then
          filename='time.par'
       else
          filename='time.ser'
       end if

      open(unit=60,file=filename,status='unknown')
      write(60,*) 'CATEGORY          min. TIME(sec)     max. TIME(sec)           PERCENT'
      total_pc=0.d0
      do i=1,ncat
        pc=100.d0*timemax(i)/real(total,kind=8)
        write(60,'(a14,2(10x,1pe9.2),10x,0pf8.1 )') cats(i),timemin(i),timemax(i),pc
        total_pc=total_pc+pc
      enddo
      write(60,'(70("-"))')
      write(60,'(a,10x,1pe9.2,6x,a,0pf5.1)') 'Total CPU time',total,'Total categorized percent ',total_pc

      close(60)

    endif
  
  else

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
!   write(*,*) 'category found',ii,cats(ii)

    if (action.eq.'ON') then  ! ON
      if (init.neqv..false.) then
        print *, cats(ii),': TIMING INITIALIZED BEFORE READ'
        stop 
      endif
      call cpu_time(time0)
      init=.true.

    else if (action.eq.'OF') then  ! OFF
      if (init.neqv..true.) then
        print *, cats(ii), 'not initialized'
        stop 
      endif
      call cpu_time(time)
      timesum(ii)=timesum(ii)+time-time0
      init=.false.
      !print *,'debug time print:,',cats(ii),' iproc=',iproc,'time',time-time0,&
      !     'time',time,'time0',time0
    else
      stop 'TIMING ACTION UNDEFINED'
  endif

endif

end subroutine timing
