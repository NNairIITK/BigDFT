
	subroutine timing(iproc,category,action)
	implicit real*8 (a-h,o-z)
	include 'mpif.h'
	character*14 category,cats
	character*2 action    ! possibilities: INitialize, ON, OFf, REsults
	logical parallel,init
	parameter(ncat=21)   ! define timimg categories
	dimension cats(ncat),flops(ncat),timesum(ncat+1),timemax(ncat+1),timemin(ncat+1)
	save cats,time0,init,timesum,total0,parallel

	data cats / 'ReformatWaves'     ,  &  !  Reformatting of input waves
                    'CrtDescriptors'    ,  &  !  Calculation of descriptor arrays
                    'CrtLocPot'         ,  &  !  Calculation of local potential
                    'CrtProjectors'     ,  &  !  Calculation of projectors
                    'ApplyLocPotKin'    ,  &  !  Application of PSP, kinetic energy
                    'ApplyProj'         ,  &  !  Application of nonlocal PSP
                    'Precondition'      ,  &  !  Precondtioning
                    'Rho_comput'        ,  &  !  Calculation of charge density (sumrho) computation
                    'Rho_commun'        ,  &  !  Calculation of charge density (sumrho) communication
                    'Un-Transall'       ,  &  !  Transposition of wavefunction, mostly communication
                    'GramS_comput'      ,  &  !  Gram Schmidt computation		
                    'GramS_commun'      ,  &  !  Gram Schmidt communication
                    'LagrM_comput'      ,  &  !  Lagrange Multipliers computation
                    'LagrM_commun'      ,  &  !  Lagrange Multipliers communication
                    'Diis'              ,  &  !
                    'PSolv_comput'      ,  &  !
                    'PSolv_commun'      ,  &  !
                    'PSolvKernel'       ,  &  !
                    'Exchangecorr'      ,  &  !
                    'Forces'            ,  &  !
                    'Tail'              /     !

! write(*,*) 'ACTION=',action,'...','CATEGORY=',category,'...'
 if (action.eq.'IN') then  ! INIT
        call cpu_time(total0)
	do i=1,ncat
	flops(i)=0.d0
  	timesum(i)=0.d0
        enddo
        parallel=trim(category).eq.'parallel'
	init=.false.

 else if (action.eq.'RE') then ! RES
!    sum results over all processor
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
	open(unit=60,file='time.prc',status='unknown')
        write(60,*) 'CATEGORY ,       min. TIME(sec) ,    max. TIME(sec),           PERCENT '
        do i=1,ncat
  	write(60,'(a14, 2(10x,e9.2),10x,f7.1 )') cats(i),timemin(i),timemax(i),100.d0*timemax(i)/total
        enddo
	endif
	
 else

        ii=100000
	do 100,i=1,ncat
	if (trim(category).eq.trim(cats(i))) then
	ii=i
	goto 200
	endif
100	continue
	print*, 'TIMING CATEGORY',category, ' NOT DEFINED'
	print*, 'ACTION  ',action
	stop 'TIMING CATEGORY NOT DEFINED'
200	continue
!        write(*,*) 'category found',ii,cats(ii)

    if (action.eq.'ON') then  ! ON
		if (init.neqv..false.) then
		print*, cats(ii),': TIMING INITIALIZED BEFORE READ'
		stop 
		endif
        call cpu_time(time0)
	init=.true.

    else if (action.eq.'OF') then  ! OFF
		if (init.neqv..true.) then
		print*, cats(ii), 'not initialized'
		stop 
		endif
        call cpu_time(time)
	timesum(ii)=timesum(ii)+time-time0
	init=.false.
    else
	stop 'TIMING ACTION UNDEFINED'
    endif

 endif

	end SUBROUTINE
