subroutine foe(iproc, nproc, tmb, orbs, evlow, evhigh, fscale, ef, tmprtr, ham, ovrlp, fermi, ebs)
  use module_base
  use module_types
  use module_interfaces
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(DFT_wavefunction),intent(inout) :: tmb
  type(orbitals_data),intent(in) :: orbs
  real(kind=8),intent(inout) :: evlow, evhigh, fscale, ef, tmprtr
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(in) :: ham, ovrlp
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(out) :: fermi
  real(kind=8),intent(out) :: ebs

  ! Local variables
  integer :: npl, istat, iall, iorb, jorb, lwork, info, ipl, korb,i, it
  integer,parameter :: nplx=5000
  real(8),dimension(:,:),allocatable :: cc, ovrlptemp2, hamscal, fermider, hamtemp, ovrlptemp
  real(8),dimension(:,:,:),allocatable :: penalty_ev
  real(kind=8),dimension(:),allocatable :: work, eval
  real(8) :: anoise, scale_factor, shift_value, ddot, tt, ttder, charge, avsumn, avsumder, sumn, sumnder
  logical :: restart
  character(len=*),parameter :: subname='foe'

  charge=0.d0
  do iorb=1,orbs%norb
       charge=charge+orbs%occup(iorb)
  end do


  ! Determine somehow evlow, evhigh, fscale, ev, tmprtr
  allocate(hamscal(tmb%orbs%norb,tmb%orbs%norb))
  allocate(ovrlptemp2(tmb%orbs%norb,tmb%orbs%norb))
  allocate(fermider(tmb%orbs%norb,tmb%orbs%norb))
  allocate(penalty_ev(tmb%orbs%norb,tmb%orbs%norb,2))


  lwork=10*tmb%orbs%norb
  allocate(work(lwork))
  allocate(eval(tmb%orbs%norb))
  allocate(hamtemp(tmb%orbs%norb,tmb%orbs%norb))
  allocate(ovrlptemp(tmb%orbs%norb,tmb%orbs%norb))
  !!hamtemp=ham
  !!ovrlptemp=ovrlp
  !!call dsygv(1, 'v', 'l', tmb%orbs%norb, hamtemp, tmb%orbs%norb, ovrlptemp, tmb%orbs%norb, eval, work, lwork, info)
  !!if (iproc==0) then
  !!    do iorb=1,tmb%orbs%norb
  !!        write(*,*) 'iorb, eval(iorb)', iorb,eval(iorb)
  !!    end do
  !!end if




  avsumn=0.d0
  avsumder=0.d0
  do it=1,100
  


      ! Scale the Hamiltonian such that all eigenvalues are in the intervall [-1:1]
      scale_factor=2.d0/(evhigh-evlow)
      shift_value=.5d0*(evhigh+evlow)
      do iorb=1,tmb%orbs%norb
          do jorb=1,tmb%orbs%norb
              hamscal(jorb,iorb)=scale_factor*(ham(jorb,iorb)-shift_value*ovrlp(jorb,iorb))
          end do
      end do

      ! Determine the degree of the polynomial
      npl=nint(2.0d0*(evhigh-evlow)/fscale)
      if (npl>nplx) stop 'npl>nplx'

      if (iproc==0) then
          write( *,'(1x,a,i0)') repeat('-',75 - int(log(real(it))/log(10.))) // ' FOE it=', it
          write(*,'(a,2x,i0,3es12.3,3x,i0)') 'FOE: it, evlow, evhigh, efermi, npl', it, evlow, evhigh, ef, npl
      end if


      allocate(cc(npl,3), stat=istat)
      call memocc(istat, cc, 'cc', subname)
      !cc=0.d0

      if (evlow>=0.d0) then
          stop 'ERROR: lowest eigenvalue must be negative'
      end if
      if (evhigh<=0.d0) then
          stop 'ERROR: highest eigenvalue must be positive'
      end if

      call CHEBFT(evlow, evhigh, npl, cc(1,1), ef, fscale, tmprtr)
      call CHDER(evlow, evhigh, cc(1,1), cc(1,2), npl)
      call CHEBFT2(evlow, evhigh, npl, cc(1,3))
      call evnoise(npl, cc(1,3), evlow, evhigh, anoise)
    
      !!if (iproc==0) then
      !!    call pltwght(npl,cc(1,1),cc(1,2),evlow,evhigh,ef,fscale,tmprtr)
      !!    call pltexp(anoise,npl,cc(1,3),evlow,evhigh)
      !!endif
    
    
      if (tmb%orbs%nspin==1) then
          do ipl=1,npl
              cc(ipl,1)=2.d0*cc(ipl,1)
              cc(ipl,2)=2.d0*cc(ipl,2)
              cc(ipl,3)=2.d0*cc(ipl,3)
          end do
      end if
    
    
    
      call chebyshev(iproc, nproc, npl, cc, tmb, hamscal, ovrlp, fermi, fermider, penalty_ev)

      !!if (iproc==0) then
      !!    write(*,'(a,3es17.7)') 'maxval(abs(penalty_ev(:,:,1))), maxval(abs(penalty_ev(:,:,2))), anoise', &
      !!                            maxval(abs(penalty_ev(:,:,1))), maxval(abs(penalty_ev(:,:,2))), anoise
      !!end if

      restart=.false.

      tt=maxval(abs(penalty_ev(:,:,2)))
      if (tt>anoise) then
          if (iproc==0) then
              write(*,'(a,2es12.3)') 'WARNING: lowest eigenvalue to high; penalty function, noise: ', tt, anoise
              write(*,'(a)') 'Increase magnitude by 20% and cycle'
          end if
          evlow=evlow*1.2d0
          restart=.true.
      end if
      tt=maxval(abs(penalty_ev(:,:,1)))
      if (tt>anoise) then
          if (iproc==0) then
              write(*,'(a,2es12.3)') 'WARNING: highest eigenvalue to high; penalty function, noise: ', tt, anoise
              write(*,'(a)') 'Increase magnitude by 20% and cycle'
          end if
          evhigh=evhigh*1.2d0
          restart=.true.
      end if

      iall=-product(shape(cc))*kind(cc)
      deallocate(cc, stat=istat)
      call memocc(istat, iall, 'cc', subname)

      if (restart) cycle

    
      ! Calculate the trace of the Fermi matrix and the derivative matrix. 
      sumn=0.d0
      sumnder=0.d0
      do iorb=1,tmb%orbs%norb
          do jorb=1,tmb%orbs%norb
          end do
          sumn=sumn+fermi(iorb,iorb)
          sumnder=sumnder+fermider(iorb,iorb)
      end do


      ef=ef+5.d-1*(sumn-charge)/sumnder
      !ef=ef-1.d0*(sumn-charge)/charge

      if (iproc==0) then
          write(*,'(a,2es17.8)') 'trace of the Fermi matrix, derivative matrix:', sumn, sumnder
          write(*,'(a,2es13.4)') 'charge difference, exit criterion:', sumn-charge, 1.d-3
          write(*,'(a,es17.8)') 'suggested Fermi energy for next iteration:', ef
      end if
  !!tt=0.d0
  !!do jorb=1,tmb%orbs%norb
  !!    do korb=1,tmb%orbs%norb
  !!        tt=tt+fermi(korb,jorb)*hamscal(korb,jorb)
  !!    end do
  !!end do
  !!if(iproc==0) write(*,*) 'tt',tt
  !!tt=tt*scale_factor-shift_value
  !!if(iproc==0) write(*,*) 'tt',tt

      if (abs(sumn-charge)<1.d-3) then
          exit
      end if
    

  end do

  if (iproc==0) then
      write( *,'(1x,a,i0)') repeat('-',84 - int(log(real(it))/log(10.)))
  end if

  scale_factor=1.d0/scale_factor
  shift_value=-shift_value
  ebs=0.d0
  do jorb=1,tmb%orbs%norb
      do korb=1,tmb%orbs%norb
          ebs = ebs + fermi(korb,jorb)*hamscal(korb,jorb)
      end do
  end do
  ebs=ebs*scale_factor-shift_value*sumn

  !!if (iproc==0) write(*,*) 'in FOE EBS',ebs




end subroutine foe




! Calculates chebychev expansion of fermi distribution.
! Taken from numerical receipes: press et al
subroutine chebft(A,B,N,cc,ef,fscale,tmprtr)
  implicit none
  
  ! Calling arguments
  real(kind=8),intent(in) :: A, B, ef, fscale, tmprtr
  integer,intent(in) :: n
  real(8),dimension(n),intent(out) :: cc

  ! Local variables
  integer :: k, j
  real(kind=8) :: bma, bpa, y, arg, fac, tt, erfcc
  real(kind=8),dimension(5000) :: cf
  real(kind=8),parameter :: pi=4.d0*atan(1.d0)

  if (n>5000) stop 'chebft'
  bma=0.5d0*(b-a)
  bpa=0.5d0*(b+a)
  do k=1,n
      y=cos(pi*(k-0.5D0)*(1.d0/n))
      arg=y*bma+bpa
      if (tmprtr.eq.0.d0) then
          cf(k)=.5d0*erfcc((arg-ef)*(1.d0/fscale))
      else
          cf(k)=1.d0/(1.d0+exp( (arg-ef)*(1.d0/tmprtr) ) )
      end if
  end do
  fac=2.d0/n
  do j=1,n
      tt=0.d0
      do  k=1,n
          tt=tt+cf(k)*cos((pi*(j-1))*((k-0.5d0)*(1.D0/n)))
      end do
      cc(j)=fac*tt
  end do

end subroutine chebft



! Calculates chebychev expansion of fermi distribution.
! Taken from numerical receipes: press et al
subroutine chebft2(a,b,n,cc)
  implicit none

  ! Calling arguments
  real(kind=8),intent(in) :: a, b
  integer,intent(in) :: n
  real(kind=8),dimension(n),intent(out) :: cc

  ! Local variables
  integer :: k, j
  real(kind=8),parameter :: pi=4.d0*atan(1.d0)
  real(kind=8) :: tt, y, arg, fac, bma, bpa
  real(kind=8),dimension(5000) :: cf

  if (n>5000) stop 'chebft2'
  bma=0.5d0*(b-a)
  bpa=0.5d0*(b+a)
  ! 3 gives broder safety zone than 4
  !tt=3.0d0*n/(B-A)
  tt=4.d0*n/(b-a)
  do k=1,n
      y=cos(pi*(k-0.5d0)*(1.d0/n))
      arg=y*bma+bpa
      cf(k)=exp((arg-b)*tt)
  end do
  fac=2.d0/n
  do j=1,n
      tt=0.d0
      do k=1,n
          tt=tt+cf(k)*cos((pi*(j-1))*((k-0.5d0)*(1.d0/n)))
      end do
      cc(j)=fac*tt
  end do
end subroutine chebft2


! Calculates chebychev expansion of the derivative of Fermi distribution.
subroutine chder(a,b,c,cder,n)
  implicit none

  ! Calling arguments
  real(kind=8),intent(in) :: a, b
  integer,intent(in) :: n
  real(8),dimension(n),intent(in) :: c
  real(8),dimension(n),intent(out) :: cder

  ! Local variables
  integer :: j
  real(kind=8) :: con

  cder(n)=0.d0
  cder(n-1)=2*(n-1)*c(n)
  if(n>=3)then
      do j=n-2,1,-1
        cder(j)=cder(j+2)+2*j*c(j+1)
      end do
  end if
  con=2.d0/(b-a)
  do j=1,n
      cder(j)=cder(j)*con
  end do

end subroutine chder




! determine noise level
subroutine evnoise(npl,cc,evlow,evhigh,anoise)
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: npl
  real(kind=8),dimension(npl),intent(in) :: cc
  real(kind=8),intent(in) :: evlow, evhigh
  real(kind=8),intent(out) :: anoise
  
  ! Local variables
  integer :: ic
  real(kind=8) :: fact, dist, ddx, cent, tt, x, chebev
  
  
  fact=1.d0
  dist=(fact*evhigh-fact*evlow)
  ddx=dist/(10*npl)
  cent=.5d0*(fact*evhigh+fact*evlow)
  ic=1
  tt=abs(chebev(evlow,evhigh,npl,cent,cc))
  ! Why use a real number as counter?!
  do x=ddx,.25d0*dist,ddx
      ic=ic+2
      tt=max(tt,abs(chebev(evlow,evhigh,npl,cent+x,cc)), &
         & abs(chebev(evlow,evhigh,npl,cent-x,cc)))
  end do
  anoise=2.d0*tt

end subroutine evnoise



! Calculates the error function complement with an error of less than 1.2E-7
  function erfcc(x)
  implicit none

  ! Calling arguments
  real(8),intent(in) :: x
  real(8) :: erfcc

  ! Local variables
  real(8) :: z, t

  z=abs(x)
  t=1.d0/(1.+0.5d0*z)
  erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+ &
        & t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+ &
        & t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
  if (x.lt.0.) erfcc=2.D0-erfcc

end function erfcc



!  evaluates chebychev expansion
function chebev(a,b,m,x,cc)
  implicit none
  
  ! Calling arguments
  real(kind=8),intent(in) :: a, b, x
  integer,intent(in) :: m
  real(kind=8),dimension(m),intent(in) :: cc
  real(kind=8) :: chebev
  
  ! Local variables
  integer :: j
  real(kind=8) :: d, dd, y, sv
  
  d=0.d0
  dd=0.d0
  y=2.d0*(2.d0*x-a-b)/(b-a)
  do j=m,2,-1
      sv=d
      d=y*d-dd+cc(j)
      dd=sv
  end do
  chebev= -dd + 0.5d0*(y*d+cc(1))

end function chebev
 



! plots the approximate fermi distribution
        subroutine pltwght(npl,cc,cder,evlow,evhigh,ef,fscale,tmprtr)
          implicit none

          ! Calling arguments
          integer,intent(in) :: npl
          real(kind=8),dimension(npl),intent(in) :: cc, cder
          real(kind=8),intent(in) :: evlow, evhigh, ef, fscale, tmprtr

          ! Local variables
          integer :: ic
          real(kind=8) :: ddx, x, tt, err, chebev

        open (unit=66,file='fermi',status='unknown')
!     header for my favourite plotting program
        write(66,*) ' 3'
        write(66,*) ' #LINETYPE{132}'
65        format(a,f5.2,a,i3,a)
        write(66,65) ' #TITLE{WEIGHT DISTR. for fscale=', fscale,' npl=',npl,'}'
        write(66,*) ' #XCAPT{ENERGY in eV}'
        write(66,*) ' #XYAPT{WEIGHT DISTR.}'
        write(66,*) ' #2YAXIS{2}'
        write(66,*) ' #YLOGSCALE2'
        write(66,*) ' #YCAPT2{ERROR}'
        write(66,*) ' $'
!
! plot chebechev expansion of weight distribution function
!
        ddx=(evhigh-evlow)/(10*npl)
! number of plot p[oints
        ic=0
        do x=evlow,evhigh,ddx
            ic=ic+1
        end do
! weight distribution
        write(66,*) ic
        do x=evlow,evhigh,ddx
            write(66,*) x,CHEBEV(evlow,evhigh,npl,x,cc)
        end do
! derivative
        write(66,*) ic
        do x=evlow,evhigh,ddx
            write(66,*) x,-CHEBEV(evlow,evhigh,npl,x,cder)
        end do
! error
        write(66,*) ic
        do x=evlow,evhigh,ddx
            tt=tmprtr
            if (tmprtr.eq.0.d0) tt=1.d-16
            err=CHEBEV(evlow,evhigh,npl,x,cc) -1.d0/(1.d0+exp((x-ef)/tt))
            write(66,*) x,err
        end do

        close(unit=66)
end subroutine pltwght




! plots the approximate fermi distribution
        subroutine pltexp(anoise,npl,cc,evlow,evhigh)
        implicit none

        ! Calling arguments
        integer,intent(in) :: npl
        real(kind=8),dimension(npl),intent(in) :: cc
        real(kind=8),intent(in) :: anoise, evlow, evhigh

        ! Local variables
        integer :: ic
        real(kind=8) :: fact, ddx, tt, chebev, x

        open (unit=66,file='exp',status='unknown')
!     header for my favourite plotting program
        write(66,*) ' 2'
        write(66,*) ' #LINETYPE{12}'
        write(66,*) ' #TITLE{exponential}'
        write(66,*) ' #YLOGSCALE'
        write(66,*) ' #XCAPT{ENERGY in eV}'
        write(66,*) ' $'
!
        fact=1.25d0
! plot chebechev expansion of weight distribution function
!
        ddx=(fact*evhigh-fact*evlow)/(10*npl)
! number of plot p[oints
        ic=0
        do x=fact*evlow,fact*evhigh,ddx
            ic=ic+1
        end do
! first curve
        write(66,*) ic
        do x=fact*evlow,fact*evhigh,ddx
        tt=CHEBEV(evlow,evhigh,npl,x,cc)
        if (abs(tt).lt.anoise) tt=anoise
            write(66,*) x,tt
        end do
! second curve
        write(66,*) ic
        do x=fact*evhigh,fact*evlow,-ddx
        tt=CHEBEV(evlow,evhigh,npl,x,cc)
        if (abs(tt).lt.anoise) tt=anoise
             write(66,*) fact*evhigh-(x-fact*evlow),tt
        end do

        close(unit=66)
        end subroutine pltexp
