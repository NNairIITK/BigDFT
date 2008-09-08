
subroutine prec_fft_slab(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,hpsi)
  ! Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
  ! hpsi is the right hand side on input and the solution on output
  use module_base
  implicit none 
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  real(wp), dimension(nvctr_c+7*nvctr_f), intent(inout) :: hpsi
  !local variables
  character(len=*), parameter :: subname='prec_fft_slab'
  integer i_stat,i_all
  real(gp), dimension(:), allocatable :: kern_k1,kern_k3
  real(wp), dimension(:,:,:), allocatable :: x_c! in and out of Fourier preconditioning
  real(wp), allocatable::z(:,:,:,:) ! work array for FFT

  call allocate_all

  ! diagonally precondition the wavelet part  
  call wscal_f(nvctr_f,hpsi(nvctr_c+1),hx,hy,hz,cprecr)

  call make_kernel(n1,hx,kern_k1)
  call make_kernel(n3,hz,kern_k3)

  call uncompress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  !	solve the helmholtz equation for the scfunction part  
  call hit_with_kernel_slab(x_c,z,kern_k1,kern_k3,n1,n2,n3,cprecr,hy)	

  call compress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  call deallocate_all

contains

  subroutine allocate_all

    allocate(kern_k1(0:n1+ndebug),stat=i_stat)
    call memocc(i_stat,kern_k1,'kern_k1',subname)
    allocate(kern_k3(0:n3+ndebug),stat=i_stat)
    call memocc(i_stat,kern_k3,'kern_k3',subname)
    allocate(z(2,0:(n1+1)/2,0:n2,0:n3+ndebug),stat=i_stat) ! work array for fft
    call memocc(i_stat,z,'z',subname)
    allocate(x_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
    call memocc(i_stat,x_c,'x_c',subname)

  end subroutine allocate_all

  subroutine deallocate_all

    i_all=-product(shape(z))*kind(z)
    deallocate(z,stat=i_stat)
    call memocc(i_stat,i_all,'z',subname)

    i_all=-product(shape(kern_k1))*kind(kern_k1)
    deallocate(kern_k1,stat=i_stat)
    call memocc(i_stat,i_all,'kern_k1',subname)

    i_all=-product(shape(kern_k3))*kind(kern_k3)
    deallocate(kern_k3,stat=i_stat)
    call memocc(i_stat,i_all,'kern_k3',subname)

    i_all=-product(shape(x_c))*kind(x_c)
    deallocate(x_c,stat=i_stat)
    call memocc(i_stat,i_all,'x_c',subname)

  end subroutine deallocate_all
end subroutine prec_fft_slab

subroutine segment_invert(n1,n2,n3,kern_k1,kern_k3,c,zx,hgrid)
  use module_base
  implicit none
  integer,intent(in)::n1,n2,n3
  real(wp),intent(in)::kern_k1(0:n1)
  real(wp),intent(in)::kern_k3(0:n3)
  real(gp),intent(in)::c,hgrid
  real(wp),intent(inout)::zx(2,0:(n1+1)/2,0:n2,0:n3)

  real(gp) ct
  real(gp),allocatable,dimension(:,:)::ab
  real(wp) b(0:n2,2)
  !     .. Scalar Arguments ..
  INTEGER::INFO, Kd, LDAB, LDB, NRHS=2,n
  integer ipiv(n2+1)
  integer i1,i2,i3,i,j

  integer,parameter::lowfil=-14,lupfil=14
  real(gp) scale
  real(gp)::fil(lowfil:lupfil)

  scale=-.5_gp/hgrid**2

! second derivative filters for Daubechies 16
fil(0)=   -3.5536922899131901941296809374_gp*scale
fil(1)=    2.2191465938911163898794546405_gp*scale
fil(2)=   -0.6156141465570069496314853949_gp*scale
fil(3)=    0.2371780582153805636239247476_gp*scale
fil(4)=   -0.0822663999742123340987663521_gp*scale
fil(5)=    0.02207029188482255523789911295638968409_gp*scale
fil(6)=   -0.409765689342633823899327051188315485e-2_gp*scale
fil(7)=    0.45167920287502235349480037639758496e-3_gp*scale
fil(8)=   -0.2398228524507599670405555359023135e-4_gp*scale
fil(9)=    2.0904234952920365957922889447361e-6_gp*scale
fil(10)=  -3.7230763047369275848791496973044e-7_gp*scale
fil(11)=  -1.05857055496741470373494132287e-8_gp*scale
fil(12)=  -5.813879830282540547959250667e-11_gp*scale
fil(13)=   2.70800493626319438269856689037647576e-13_gp*scale
fil(14)=  -6.924474940639200152025730585882e-18_gp*scale

do i=1,14
	fil(-i)=fil(i)
enddo

fil=fil*(n1+1)*(n3+1) ! include the factor from the Fourier transform

! initialize the variables for matrix inversion
n=n2+1
kd=lupfil
ldab=kd+1!
ldb=n2+1
allocate(ab(ldab,n))

! hit the fourier transform of x with the kernel
do i3=0,n3
!	do i1=0,n1
	do i1=0,(n1+1)/2
!		ct=kern_k1(i1)+kern_k3(i3)+c
		ct=(kern_k1(i1)+kern_k3(i3)+c)*(n1+1)*(n3+1)

		! ab has to be reinitialized each time
		! since it is overwritten in the course of  dgbsv
		do j=1,n
			do i=max(1,j-lupfil),j
				ab(kd+1+i-j,j)=fil(i-j)
			enddo
			ab(kd+1,j)=fil(0)+ct
		enddo

		do i2=0,n2
			b(i2,1)=zx(1,i1,i2,i3)
			b(i2,2)=zx(2,i1,i2,i3)
		enddo

!		call DGBSV( N, KL, KU, NRHS, ab , LDAB, IPIV, b, LDB, INFO )
        call DPBSV( 'U', N, KD, NRHS, AB, LDAB, B, LDB, INFO )
		if (info.ne.0) stop 'error in matrix inversion'

		do i2=0,n2
			zx(1,i1,i2,i3)=b(i2,1)
			zx(2,i1,i2,i3)=b(i2,2)
		enddo

	enddo
enddo

deallocate(ab)
end subroutine segment_invert






