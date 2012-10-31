subroutine foe(iproc, nproc, tmb, evlow, evhigh, fscale, ef, tmprtr, ham, ovrlp, fermi)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(DFT_wavefunction),intent(in) :: tmb
  real(kind=8),intent(inout) :: evlow, evhigh, fscale, ef, tmprtr
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(inout) :: ham, ovrlp, fermi

  ! Local variables
  integer :: npl, istat, iall, iorb, jorb, lwork, info
  integer,parameter :: nplx=200
  real(8),dimension(:,:),allocatable :: cc, hamtemp
  real(kind=8),dimension(:),allocatable :: work, eval
  real(8) :: anoise, tt1, tt2
  character(len=*),parameter :: subname='foe'

  ! Determine somehow evlow, evhigh, fscale, ev, tmprtr
  lwork=10*tmb%orbs%norb
  allocate(work(lwork))
  allocate(eval(tmb%orbs%norb))
  allocate(hamtemp(tmb%orbs%norb,tmb%orbs%norb))
  hamtemp=ham
  call dsyev('n', 'l', tmb%orbs%norb, hamtemp, tmb%orbs%norb, eval, work, lwork, info)
  evlow=eval(1)
  evhigh=eval(tmb%orbs%norb)

  if (iproc==0) then
      write(*,*) 'BEFORE: lowest eval', eval(1)
      write(*,*) 'BEFORE: highest eval', eval(tmb%orbs%norb)
  end if

  ! Determine the degree of the polynomial
  npl=nint(2.0d0*(evhigh-evlow)/fscale)
  if(iproc==0) write(*,*) 'npl',npl
  if (npl>nplx) stop 'npl>nplx'
  
  allocate(cc(npl,3), stat=istat)
  call memocc(istat, cc, 'cc', subname)

  call CHEBFT(evlow, evhigh, npl, cc(1,1), ef, fscale, tmprtr)
  call CHDER(evlow, evhigh, cc(1,1), cc(1,2), npl)
  call CHEBFT2(evlow, evhigh, npl, cc(1,3))
  call evnoise(npl, cc(1,3), evlow, evhigh, anoise)

  ! Scale the Hamiltonian such that all eigenvalues are in the intervall [-1:1]
  tt1=2.d0/(evhigh-evlow)
  tt2=.5d0*(evhigh+evlow)
  do iorb=1,tmb%orbs%norb
      do jorb=1,tmb%orbs%norb
          if (iorb==jorb) then
              ham(jorb,iorb)=tt1*(ham(jorb,iorb)-tt2)
          else
              ham(jorb,iorb)=tt1*ham(jorb,iorb)
          end if
      end do
  end do

  hamtemp=ham
  call dsyev('n', 'l', tmb%orbs%norb, hamtemp, tmb%orbs%norb, eval, work, lwork, info)

  if (iproc==0) then
     write(*,*) 'AFTER: lowest eval', eval(1)
     write(*,*) 'AFTER: highest eval', eval(tmb%orbs%norb)
 end if


  ! Unscale spectrum
  tt1=1.d0/tt1
  tt2=-tt2
  do iorb=1,tmb%orbs%norb
      eval(iorb)=tt1*eval(iorb)-tt2
  end do
  do iorb=1,tmb%orbs%norb
      do jorb=1,tmb%orbs%norb
          if (iorb==jorb) then
              fermi(jorb,iorb)=tt1*fermi(jorb,iorb)-tt2
          else
              fermi(jorb,iorb)=tt1*fermi(jorb,iorb)
          end if
      end do
  end do

  if (iproc==0) then
     write(*,*) 'AFTER UNSCALE: lowest eval', eval(1)
     write(*,*) 'AFTER UNSCALE: highest eval', eval(tmb%orbs%norb)
 end if

  iall=-product(shape(cc))*kind(cc)
  deallocate(cc, stat=istat)
  call memocc(istat, iall, 'cc', subname)



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
  real(kind=8),dimension(1000) :: cf
  real(kind=8),parameter :: pi=4.d0*atan(1.d0)

  if (n>1000) stop 'chebft'
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
  real(kind=8),dimension(1000) :: cf

  if (n>1000) stop 'chebft2'
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

