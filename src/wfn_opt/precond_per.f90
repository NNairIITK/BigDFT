
!!$subroutine preconditionall(iproc,nproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hgrid,&
!!$     ncong,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,eval,&
!!$     ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi)
!!$  ! Calls the preconditioner for each orbital treated by the processor
!!$  implicit real(kind=8) (a-h,o-z)
!!$  dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
!!$  dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
!!$  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
!!$  dimension hpsi(nvctr_c+7*nvctr_f,norbp),eval(norb)
!!$  
!!$  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
!!$     
!!$     cprecr=-eval(iorb)
!!$     !       write(*,*) 'cprecr',iorb,cprecr!
!!$     call prec_fft(iorb,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
!!$          nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
!!$          ncong,cprecr,hgrid,hpsi(1,iorb-iproc*norbp))
!!$     
!!$  enddo
!!$
!!$end subroutine preconditionall

subroutine prec_fft(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     ncong,cprecr,hx,hy,hz,hpsi)
  ! Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
  ! hpsi is the right hand side on input and the solution on output
  implicit none 
  integer, intent(in) :: n1,n2,n3,ncong
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(kind=8), intent(in) :: hx,hy,hz,cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(kind=8), intent(out) ::  hpsi(nvctr_c+7*nvctr_f)
  !local variables
  integer nd1,nd2,nd3
  real(kind=8), dimension(:), allocatable :: kern_k1,kern_k2,kern_k3
  real(kind=8), dimension(:,:,:), allocatable :: x_c,y_c! in and out of Fourier preconditioning
  real(kind=8), dimension(:,:,:,:,:), allocatable::z(:,:,:,:,:) ! work array for FFT


  nd1=n1+2; nd2=n2+2; nd3=n3+2
  !these allocation statements must be changed with the memory profiling method
  allocate(kern_k1(0:n1))
  allocate(kern_k2(0:n2))
  allocate(kern_k3(0:n3))


  allocate(z(2,nd1,nd2,nd3,2)) ! work array for fft
  allocate(x_c(0:n1,0:n2,0:n3))
  allocate(y_c(0:n1,0:n2,0:n3))

  call wscal_f(nvctr_f,hpsi(nvctr_c+1),hx,hy,hz,cprecr)

  call make_kernel(n1,hx,kern_k1)
  call make_kernel(n2,hy,kern_k2)
  call make_kernel(n3,hz,kern_k3)

  call uncompress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  call  hit_with_kernel(x_c,y_c,z,kern_k1,kern_k2,kern_k3,n1,n2,n3,nd1,nd2,nd3,cprecr)

  call   compress_c(hpsi,y_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  deallocate(z,kern_k1,kern_k2,kern_k3)

  deallocate(x_c,y_c)
end subroutine prec_fft



subroutine uncompress_c(hpsi,x_c,keyg_c,keyv_c,nseg_c,nvctr_c,n1,n2,n3)
implicit none
integer,intent(in)::n1,n2,n3
integer,intent(in)::nseg_c,nvctr_c
real*8,intent(in)::hpsi(nvctr_c)
integer,intent(in)::keyg_c(2,nseg_c),keyv_c(nseg_c)

real*8,intent(out)::x_c(0:n1,0:n2,0:n3)

integer iseg,jj,j0,j1,ii,i3,i2,i0,i1,i

x_c=0.d0
do iseg=1,nseg_c
	jj=keyv_c(iseg)
	j0=keyg_c(1,iseg)
	j1=keyg_c(2,iseg)
			  
	ii=j0-1
	i3=ii/((n1+1)*(n2+1))
	ii=ii-i3*(n1+1)*(n2+1)
	i2=ii/(n1+1)
	i0=ii-i2*(n1+1)
	i1=i0+j1-j0
	
	do i=i0,i1
		x_c(i,i2,i3)=hpsi(i-i0+jj)
	enddo
enddo
end subroutine uncompress_c

subroutine compress_c(hpsi,y_c,keyg_c,keyv_c,nseg_c,nvctr_c,n1,n2,n3)
implicit none
integer,intent(in)::n1,n2,n3
integer,intent(in)::nseg_c,nvctr_c
integer,intent(in)::keyg_c(2,nseg_c),keyv_c(nseg_c)
real*8,intent(in)::y_c(0:n1,0:n2,0:n3)

real*8,intent(out)::hpsi(nvctr_c)

integer iseg,jj,j0,j1,ii,i3,i2,i0,i1,i

! coarse part
do iseg=1,nseg_c
	jj=keyv_c(iseg)
	j0=keyg_c(1,iseg)
	j1=keyg_c(2,iseg)
	ii=j0-1
	i3=ii/((n1+1)*(n2+1))
	ii=ii-i3*(n1+1)*(n2+1)
	i2=ii/(n1+1)
	i0=ii-i2*(n1+1)
	i1=i0+j1-j0
	do i=i0,i1
		hpsi(i-i0+jj)=y_c(i,i2,i3)
	enddo
enddo

end subroutine compress_c


subroutine wscal_f(mvctr_f,psi_f,hx,hy,hz,c)
! multiplies a wavefunction psi_c,psi_f (in vector form) with a scaling vector (scal)
implicit real(kind=8) (a-h,o-z)
real*8,intent(in)::c,hx,hy,hz
dimension psi_f(7,mvctr_f),scal(7),hh(3)
!WAVELET AND SCALING FUNCTION SECOND DERIVATIVE FILTERS, diagonal elements
PARAMETER(B2=24.8758460293923314D0,A2=3.55369228991319019D0)
  
hh(1)=.5d0/hx**2
hh(2)=.5d0/hy**2
hh(3)=.5d0/hz**2

scal(1)=1.d0/(b2*hh(1)+a2*hh(2)+a2*hh(3)+c)       !  2 1 1
scal(2)=1.d0/(a2*hh(1)+b2*hh(2)+a2*hh(3)+c)       !  1 2 1
scal(3)=1.d0/(b2*hh(1)+b2*hh(2)+a2*hh(3)+c)       !  2 2 1
scal(4)=1.d0/(a2*hh(1)+a2*hh(2)+b2*hh(3)+c)       !  1 1 2
scal(5)=1.d0/(b2*hh(1)+a2*hh(2)+b2*hh(3)+c)       !  2 1 2
scal(6)=1.d0/(a2*hh(1)+b2*hh(2)+b2*hh(3)+c)       !  1 2 2
scal(7)=1.d0/(b2*hh(1)+b2*hh(2)+b2*hh(3)+c)       !  2 2 2
  
do i=1,mvctr_f
	psi_f(1,i)=psi_f(1,i)*scal(1)       !  2 1 1
	psi_f(2,i)=psi_f(2,i)*scal(2)       !  1 2 1
	psi_f(3,i)=psi_f(3,i)*scal(3)       !  2 2 1
	psi_f(4,i)=psi_f(4,i)*scal(4)       !  1 1 2
	psi_f(5,i)=psi_f(5,i)*scal(5)       !  2 1 2
	psi_f(6,i)=psi_f(6,i)*scal(6)       !  1 2 2
	psi_f(7,i)=psi_f(7,i)*scal(7)       !  2 2 2
enddo

end subroutine wscal_f

