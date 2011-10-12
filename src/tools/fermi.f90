

   implicit real*8 (a-h,o-z)
   parameter(n=100)
   dimension eval(n)


   eval(1)=0.d0
   do i=2,n
   call random_number(ts)
   eval(i)=eval(i-1)+10*ts/n
!   eval(i)=10*ts
   write(*,*) eval(i)
   enddo

   iproc=0
   melec=50
   ef=0.d0
   wf=.2d0
   call Fermilevel(iproc,n,eval,melec,ef,wf)

   end


   subroutine Fermilevel(iproc,n,eval,melec,ef,wf)
! finds  the fermi level ef for an error function distribution with a width wf
! eval are the Kohn Sham eigenvalues and melec is the total number of electrons
   implicit real*8 (a-h,o-z)
   parameter(pi=3.1415926535897932d0)
   dimension eval(n)

   factor=1.d0/(sqrt(pi)*wf)
10 electrons=0.d0
   dlectrons=0.d0
   do i=1,n
   arg=(eval(i)-ef)/wf
! next 2 line error function distribution
!   electrons=electrons+.5d0*(1.d0-derf(arg))
!   dlectrons=dlectrons-exp(-arg**2)
! next 2 line Fermi function distribution
   electrons=electrons+1.d0/(1.d0+exp(arg))
   dlectrons=dlectrons-exp(arg)/(1.d0+exp(arg))**2
   enddo
! next  line error function distribution
!   dlectrons=dlectrons*factor
! next  line Fermi function distribution
   dlectrons=dlectrons/wf

   write(*,*) electrons,ef,dlectrons
   diff=melec-electrons
   if (abs(diff).lt.1.d-12) goto 20
   corr=diff/dlectrons
   if (corr.gt.1.d0*wf) corr=1.d0*wf
   if (corr.lt.-1.d0*wf) corr=-1.d0*wf
   if (abs(dlectrons).lt.1.d-14  .and. electrons.gt.dble(melec)) corr=3.d0*wf
   if (abs(dlectrons).lt.1.d-14  .and. electrons.lt.dble(melec)) corr=-3.d0*wf
   ef=ef-corr

   goto 10
20 continue
   arg=(eval(n)-ef)/wf
!   cutoff=.5d0*(1.d0-derf(arg))
   cutoff=1.d0/(1.d0+exp(arg))
   if (iproc.eq.0) write(*,'(a,f12.5,e8.1)') 'Fermi level, Fermi distribution cut off at:',ef,cutoff
    
   end subroutine

