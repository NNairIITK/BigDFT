	subroutine timing(category,action,nflop,me,nproc)
	implicit real*8 (a-h,o-z)
	include 'mpif.h'
	character*6 category,cats
	character*3 action
	logical init
	parameter(ncat=18)
	dimension cats(ncat),flops(ncat),timesum(ncat+1),
     1  wrktime(ncat+1),wrksquare(ncat),wrkflops(ncat),wrkflop2(ncat)
	save cats,time0,init,timesum,flops,total0

	data cats / 'BLASL3' , 'BLASL2' , 'RWVFFT' , 'POTFFT' , 
     1              'SLFPOT' , 'ROTATN' , 'REDUCT' , 'ANALYS' ,
     1              'IGUESS' , 'TETERF' , 'PROJE1' , 'PROJE2' , 
     1              'COMOVR' , 'MONOPL' , 'FORCES' , 
     1              'DIIS00' , 'POISSO' , 'ADDNUC' /

	if (action.eq.'ZER') then
	do 11,i=1,ncat
	total0=mclock()
	flops(i)=0.d0
11	timesum(i)=0.d0
	init=.false.
	else if (action.eq.'RES') then
c sum results over all processor
	timesum(ncat+1)=mclock()-total0
        if (nproc.gt.1) then 
	call MPI_REDUCE(timesum,wrktime,ncat+1,
     1       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	else
	do 488,i=1,ncat+1
488	wrktime(i)=timesum(i)
	endif
	do 35,i=1,ncat
35	timesum(i)=timesum(i)**2
	if (nproc.gt.1) then
        call MPI_REDUCE(timesum,wrksquare,ncat,
     1       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(flops,wrkflops,ncat,
     1       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	else
	do 124,i=1,ncat
	wrksquare(i)=timesum(i)
124	wrkflops(i)=flops(i)
	endif
	do 45,i=1,ncat
45	flops(i)=flops(i)**2
	if (nproc.gt.1) then
        call MPI_REDUCE(flops,wrkflop2,ncat,
     1       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	else
	do 874,i=1,ncat
874	wrkflop2(i)=flops(i)
	endif

	if (me.eq.0) then
	open(unit=60,file='time.prc',status='unknown')
        write(60,*) 
     1  'CATEGORY , TIME(sec)+dev , PERCENT ,   FLOPS+dev ,      SPEED '
25	format(a, e10.3 , e9.2 ,2x, f4.3 ,2x, e10.3 , e9.2 ,2x, e10.3)
	sum=0.d0
	do 22,i=1,ncat
22      sum=sum+wrktime(i)
	do 20,i=1,ncat
	t1=wrktime(i)
	t2=sqrt(-t1**2+nproc*wrksquare(i))
	t3=wrktime(i)/sum
	t4=wrkflops(i)
	t5=sqrt(-t4**2+nproc*wrkflop2(i))
	t6=1.d-4*wrkflops(i)/wrktime(i)
20	write(60,25) cats(i),1.d-2*t1/nproc,1.d-2*t2/nproc,t3,
     1               t4/nproc,t5/nproc,t6
	write(60,25) 'SUM  ',1.d-2*sum/nproc
	write(60,25) 'TOTAL',1.d-2*wrktime(ncat+1)/nproc
	endif
	
	else

	do 100,i=1,ncat
	if (category.eq.cats(i)) then
	ii=i
	goto 200
	endif
100	continue
	print*, 'TIMING CATEGORY',category, ' NOT DEFINED'
	print*, 'ACTION  ',action
	stop 'TIMING CATEGORY NOT DEFINED'
200	continue

	if (action.eq.'INI') then
		if (init.neqv..false.) then
		print*, cats(ii),': TIMING INITIALIZED BEFORE READ'
		stop 
		endif
	time0=mclock()
	init=.true.
	else if (action.eq.'SUM') then
		if (init.neqv..true.) then
		print*, cats(ii), 'not initialized'
		stop 
		endif
	flops(ii)=flops(ii)+nflop
	timesum(ii)=timesum(ii)+mclock()-time0
	init=.false.
	else
	stop 'TIMING ACTION UNDEFINED'
	endif

	endif

	end

