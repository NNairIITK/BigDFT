subroutine precong_slab(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     ncong,cprecr,hx,hy,hz,x)
  use module_base
  ! Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
  ! x is the right hand side on input and the solution on output
  ! This subroutine is almost similar to precong_per
  ! with minor differences: 
  !
  ! -the convolutions are in the surface BC
  ! -the work arrays psifscf and ww are slightly bigger
  ! -the work array z is smaller now
  ! -the Fourier preconditioner is in the surface version
  implicit none
  integer, intent(in) :: n1,n2,n3,ncong
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(inout) ::  x(nvctr_c+7*nvctr_f)
  ! local variables
  real(gp)::scal(0:8)
  real(wp)::rmr,rmr_new,alpha,beta
  integer i,i_stat,i_all
  real(wp),allocatable::b(:),r(:),d(:)
  real(wp),allocatable::psifscf(:),ww(:)

  integer, parameter :: lowfil=-14,lupfil=14
  real(gp), allocatable,dimension(:,:) :: af,bf,cf,ef
  integer,allocatable,dimension(:)::modul1,modul3

  call allocate_all

  call prepare_sdc_slab(n1,n2,n3,modul1,modul3,af,bf,cf,ef,hx,hy,hz)

  !	initializes the wavelet scaling coefficients	
  call wscal_init_per(scal,hx,hy,hz,cprecr)
  b=x

  !	compute the input guess x via a Fourier transform in a cubic box.
  !	Arrays psifscf and ww serve as work arrays for the Fourier
  call prec_fft_slab_fast(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
       cprecr,hx,hy,hz,x,&
       psifscf(1),psifscf(n1+2),ww(1),ww(2*((n1+1)/2+1)*(n2+1)*(n3+1)+1))

  !	call apply_hp_slab(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
  !     cprecr,hx,hy,hz,x,d,psifscf,ww) ! d:=Ax
  call apply_hp_slab_sd(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
       cprecr,hx,hy,hz,x,d,psifscf,ww,modul1,modul3,af,bf,cf,ef) ! d:=Ax
  r=b-d

  call wscal_per(nvctr_c,nvctr_f,scal,r(1),r(nvctr_c+1),d(1),d(nvctr_c+1))
  rmr=dot_product(r,d)
  do i=1,ncong 
     !		write(*,*)i,sqrt(rmr)

     !		call apply_hp_slab(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     !	     cprecr,hx,hy,hz,d,b,psifscf,ww) ! b:=Ad
     call apply_hp_slab_sd(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
          cprecr,hx,hy,hz,d,b,psifscf,ww,modul1,modul3,af,bf,cf,ef) ! b:=Ad

     alpha=rmr/dot_product(d,b)
     x=x+alpha*d
     r=r-alpha*b

     call wscal_per(nvctr_c,nvctr_f,scal,r(1),r(nvctr_c+1),b(1),b(nvctr_c+1))
     rmr_new=dot_product(r,b)

     beta=rmr_new/rmr
     d=b+beta*d
     rmr=rmr_new
  enddo

  call deallocate_all

contains
  subroutine allocate_all

    allocate(modul1(lowfil:n1+lupfil),stat=i_stat)
    call memocc(i_stat,product(shape(modul1))*kind(modul1),'modul1','precong_per')
    allocate(modul3(lowfil:n3+lupfil),stat=i_stat)
    call memocc(i_stat,product(shape(modul3))*kind(modul3),'modul3','precong_per')
    allocate(af(lowfil:lupfil,3),stat=i_stat)
    call memocc(i_stat,product(shape(af))*kind(af),'af','precong_per')
    allocate(bf(lowfil:lupfil,3),stat=i_stat)
    call memocc(i_stat,product(shape(bf))*kind(bf),'bf','precong_per')
    allocate(cf(lowfil:lupfil,3),stat=i_stat)
    call memocc(i_stat,product(shape(cf))*kind(cf),'cf','precong_per')
    allocate(ef(lowfil:lupfil,3),stat=i_stat)
    call memocc(i_stat,product(shape(ef))*kind(ef),'ef','precong_per')

    allocate(b(nvctr_c+7*nvctr_f),stat=i_stat)
    call memocc(i_stat,product(shape(b))*kind(b),'b','precong_per')
    allocate(r(nvctr_c+7*nvctr_f),stat=i_stat)
    call memocc(i_stat,product(shape(r))*kind(r),'r','precong_per')
    allocate(d(nvctr_c+7*nvctr_f),stat=i_stat)
    call memocc(i_stat,product(shape(d))*kind(d),'','precong_per')
    allocate( psifscf((2*n1+2)*(2*n2+16)*(2*n3+2)),stat=i_stat )
    call memocc(i_stat,product(shape(psifscf))*kind(psifscf),'psifscf','precong_per')
    allocate( ww((2*n1+2)*(2*n2+16)*(2*n3+2)) ,stat=i_stat)
    call memocc(i_stat,product(shape(ww))*kind(ww),'ww','precong_per')
  end subroutine allocate_all

  subroutine deallocate_all

    i_all=-product(shape(modul1))*kind(modul1)
    deallocate(modul1,stat=i_stat)
    call memocc(i_stat,i_all,'modul1','last_orthon')

    i_all=-product(shape(modul3))*kind(modul3)
    deallocate(modul3,stat=i_stat)
    call memocc(i_stat,i_all,'modul3','last_orthon')

    i_all=-product(shape(af))*kind(af)
    deallocate(af,stat=i_stat)
    call memocc(i_stat,i_all,'af','last_orthon')

    i_all=-product(shape(bf))*kind(bf)
    deallocate(bf,stat=i_stat)
    call memocc(i_stat,i_all,'bf','last_orthon')

    i_all=-product(shape(cf))*kind(cf)
    deallocate(cf,stat=i_stat)
    call memocc(i_stat,i_all,'cf','last_orthon')

    i_all=-product(shape(ef))*kind(ef)
    deallocate(ef,stat=i_stat)
    call memocc(i_stat,i_all,'ef','last_orthon')


    i_all=-product(shape(psifscf))*kind(psifscf)
    deallocate(psifscf,stat=i_stat)
    call memocc(i_stat,i_all,'psifscf','last_orthon')

    i_all=-product(shape(ww))*kind(ww)
    deallocate(ww,stat=i_stat)
    call memocc(i_stat,i_all,'ww','last_orthon')

    i_all=-product(shape(b))*kind(b)
    deallocate(b,stat=i_stat)
    call memocc(i_stat,i_all,'b','last_orthon')

    i_all=-product(shape(r))*kind(r)
    deallocate(r,stat=i_stat)
    call memocc(i_stat,i_all,'r','last_orthon')

    i_all=-product(shape(d))*kind(d)
    deallocate(d,stat=i_stat)
    call memocc(i_stat,i_all,'d','last_orthon')
  end subroutine deallocate_all
end subroutine precong_slab

subroutine apply_hp_slab(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,x,y,psifscf,ww)
  !	Applies the operator (KE+cprecr*I)*x=y
  !	array x is input, array y is output
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(in) ::  x(nvctr_c+7*nvctr_f)  
  real(wp), intent(out) ::  y(nvctr_c+7*nvctr_f)

  real(gp) hgridh(3)	
  real(wp),dimension((2*n1+2)*(2*n2+16)*(2*n3+2))::ww,psifscf

  ! x: input, ww:work
  ! psifscf: output
  call uncompress_slab(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
       x(1),x(nvctr_c+1),psifscf,ww)

  hgridh(1)=hx*.5_gp
  hgridh(2)=hy*.5_gp
  hgridh(3)=hz*.5_gp
  ! psifscf: input, ww: output
  call convolut_kinetic_slab_c(2*n1+1,2*n2+15,2*n3+1,hgridh,psifscf,ww,cprecr)

  ! ww:intput, psifscf: work
  ! y:output
  call compress_slab(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),& 
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),& 
       ww,y(1),y(nvctr_c+1),psifscf)
end subroutine apply_hp_slab


subroutine apply_hp_slab_sd(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,x,y,psifscf,ww,modul1,modul3,a,b,c,e)
  !	Applies the operator (KE+cprecr*I)*x=y
  !	array x is input, array y is output
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(in) ::  x(nvctr_c+7*nvctr_f)  
  real(wp), intent(out) ::  y(nvctr_c+7*nvctr_f)
  integer, parameter :: lowfil=-14,lupfil=14

  real(gp) hgrid(3)	
  real(wp),dimension((2*n1+2)*(2*n2+16)*(2*n3+2))::ww,psifscf

  integer,intent(in)::modul1(lowfil:n1+lupfil)
  integer,intent(in)::modul3(lowfil:n3+lupfil)
  real(gp),intent(in)::a(lowfil:lupfil,3)
  real(gp),intent(in)::b(lowfil:lupfil,3)
  real(gp),intent(in)::c(lowfil:lupfil,3)
  real(gp),intent(in)::e(lowfil:lupfil,3)
  ! x: input
  ! psifscf: output
  call uncompress_sd(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
       x(1),x(nvctr_c+1),psifscf)

  hgrid(1)=hx
  hgrid(2)=hy
  hgrid(3)=hz
  ! psifscf: input, ww: output
  !  	call convolut_kinetic_slab_c(2*n1+1,2*n2+15,2*n3+1,hgridh,psifscf,ww,cprecr)
  call convolut_kinetic_slab_sdc(n1,n2,n3,hgrid,psifscf,ww,cprecr,modul1,modul3,a,b,c,e)

  ! ww:intput
  ! y:output
  call compress_sd(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),& 
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),& 
       ww,y(1),y(nvctr_c+1))
end subroutine apply_hp_slab_sd



subroutine prec_fft_slab_fast(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,hpsi,kern_k1,kern_k3,z,x_c)
  ! Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
  ! hpsi is the right hand side on input and the solution on output
  use module_base
  implicit none 
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(inout) ::  hpsi(nvctr_c+7*nvctr_f)
  !local variables
  integer i_stat,i_all
  real(gp)::kern_k1(0:n1),kern_k3(0:n3)
  real(wp)::x_c(0:n1,0:n2,0:n3)! in and out of Fourier preconditioning
  real(wp)::z(2,0:(n1+1)/2,0:n2,0:n3)! work array for FFT

  ! diagonally precondition the wavelet part  
  call wscal_f(nvctr_f,hpsi(nvctr_c+1),hx,hy,hz,cprecr)

  call make_kernel(n1,hx,kern_k1)
  call make_kernel(n3,hz,kern_k3)

  call uncompress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  !	solve the helmholtz equation for the scfunction part  
  call hit_with_kernel_slab(x_c,z,kern_k1,kern_k3,n1,n2,n3,cprecr,hy)	

  call   compress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

end subroutine prec_fft_slab_fast




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
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(inout) ::  hpsi(nvctr_c+7*nvctr_f)
  !local variables
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

  call   compress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  call deallocate_all

contains
  subroutine allocate_all
    allocate(kern_k1(0:n1),stat=i_stat)
    call memocc(i_stat,product(shape(kern_k1))*kind(kern_k1),'kern_k1','prec_fft')
    allocate(kern_k3(0:n3),stat=i_stat)
    call memocc(i_stat,product(shape(kern_k3))*kind(kern_k3),'kern_k3','prec_fft')
    allocate(z(2,0:(n1+1)/2,0:n2,0:n3),stat=i_stat) ! work array for fft
    call memocc(i_stat,product(shape(z))*kind(z),'z','prec_fft')
    allocate(x_c(0:n1,0:n2,0:n3),stat=i_stat)
    call memocc(i_stat,product(shape(x_c))*kind(x_c),'x_c','prec_fft')
  end subroutine allocate_all
  subroutine deallocate_all
    i_all=-product(shape(z))*kind(z)
    deallocate(z,stat=i_stat)
    call memocc(i_stat,i_all,'z','last_orthon')

    i_all=-product(shape(kern_k1))*kind(kern_k1)
    deallocate(kern_k1,stat=i_stat)
    call memocc(i_stat,i_all,'kern_k1','last_orthon')

    i_all=-product(shape(kern_k3))*kind(kern_k3)
    deallocate(kern_k3,stat=i_stat)
    call memocc(i_stat,i_all,'kern_k3','last_orthon')

    i_all=-product(shape(x_c))*kind(x_c)
    deallocate(x_c,stat=i_stat)
    call memocc(i_stat,i_all,'x_c','last_orthon')

  end subroutine deallocate_all
end subroutine prec_fft_slab

subroutine segment_invert(n1,n2,n3,kern_k1,kern_k3,c,zx,hgrid)
  ! solve the discretized equation
  ! (-d^2/dy^2+ct(k1,k3)) zx(output) = zx(input)
  ! for all k1,k3 via Lapack
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






