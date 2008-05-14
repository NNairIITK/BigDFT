
subroutine hit_with_kernel(x,y,z,kern_k1,kern_k2,kern_k3,n1,n2,n3,nd1,nd2,nd3,c)	
! hits the input array x with the kernel
! ((-1/2\Delta+C)_{ij})^{-1}
	implicit none
	integer,intent(in)::n1,n2,n3,nd1,nd2,nd3
	real*8,intent(in)::x(0:n1,0:n2,0:n3)! input
	real*8,intent(in)::kern_k1(0:n1)
	real*8,intent(in)::kern_k2(0:n2)
	real*8,intent(in)::kern_k3(0:n3)
	real*8,intent(in)::c

	real*8,intent(out)::y(0:n1,0:n2,0:n3)! output

	real*8::z(2,nd1,nd2,nd3,2)! work array
	real*8 tt
	integer i1,i2,i3,isign,inzee

! fft the input array x:
do i3=0,n3
	do i2=0,n2
		do i1=0,n1
			z(1,i1+1,i2+1,i3+1,1)=x(i1,i2,i3)
			z(2,i1+1,i2+1,i3+1,1)=0.d0
		enddo
	enddo
enddo

inzee=1; isign=1
call fft(n1+1,n2+1,n3+1,nd1,nd2,nd3,z,isign,inzee)

! hit the fourier transform of x with the kernel
do i3=0,n3
	do i2=0,n2
		do i1=0,n1
			tt=kern_k1(i1)+kern_k2(i2)+kern_k3(i3)+c
			z(:,i1+1,i2+1,i3+1,inzee)=z(:,i1+1,i2+1,i3+1,inzee)/tt
		enddo
	enddo
enddo

! fourier transform x back; the result is y
isign=-1
call fft(n1+1,n2+1,n3+1,nd1,nd2,nd3,z,isign,inzee)

do i3=0,n3
	do i2=0,n2
		do i1=0,n1
			y(i1,i2,i3)=z(1,i1+1,i2+1,i3+1,inzee)/((n1+1)*(n2+1)*(n3+1))
		enddo
	enddo
enddo

end subroutine hit_with_kernel




subroutine make_kernel(n1,hgrid,kern)
! construct the kernel (-1/2 d^2/dx^2)_{ij}
! at a real space grid with grid size hgrid
! and then fourier transform it to momentum space
implicit none
integer,intent(in)::n1
real*8,intent(in)::hgrid
integer,parameter::lowfil=-14,lupfil=14

real*8,intent(out)::kern(0:n1)

real*8 scale
real*8 fil(lowfil:lupfil)

!***********************************************************************************************
integer::now(7),after(7),before(7),isign=1,ic
real*8,allocatable::trig(:,:),z(:,:,:)
integer inzee,i,nd1
!***********************************************************************************************

scale=-.5d0/hgrid**2

! second derivative filters for Daubechies 16
fil(0)=   -3.5536922899131901941296809374d0*scale
fil(1)=    2.2191465938911163898794546405d0*scale
fil(2)=   -0.6156141465570069496314853949d0*scale
fil(3)=    0.2371780582153805636239247476d0*scale
fil(4)=   -0.0822663999742123340987663521d0*scale
fil(5)=    0.02207029188482255523789911295638968409d0*scale
fil(6)=   -0.409765689342633823899327051188315485d-2*scale
fil(7)=    0.45167920287502235349480037639758496d-3*scale
fil(8)=   -0.2398228524507599670405555359023135d-4*scale
fil(9)=    2.0904234952920365957922889447361d-6*scale
fil(10)=  -3.7230763047369275848791496973044d-7*scale
fil(11)=  -1.05857055496741470373494132287d-8*scale
fil(12)=  -5.813879830282540547959250667d-11*scale
fil(13)=   2.70800493626319438269856689037647576d-13*scale
fil(14)=  -6.924474940639200152025730585882d-18*scale

! construct the kernel in real space
kern=0.d0
kern(0)=fil(0)

do i=1,14
	kern(i)     =fil(i)
	kern(n1+1-i)=fil(i)
enddo


!***********************************************************************************************
! fourier transform the kernel

nd1=n1+2
allocate(trig(2,1024))
allocate(z(2,nd1,2))

inzee=1
do i=0,n1
	z(1,i+1,inzee)=kern(i)
	z(2,i+1,inzee)=0.d0
enddo

isign=1
call ctrig(n1+1,trig,after,before,now,isign,ic)
do i=1,ic
	call fftstp(1,1,nd1,1,nd1,z(1,1,inzee),z(1,1,3-inzee),trig,after(i),now(i),before(i),isign)
	inzee=3-inzee
enddo

do i=0,n1
	kern(i)=z(1,i+1,inzee)
enddo

deallocate(trig,z)
end subroutine make_kernel