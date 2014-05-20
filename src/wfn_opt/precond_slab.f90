!> @file
!!  Routines to precondition
!! @author
!!    Copyright (C) 2010-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
!! x is the right hand side on input and the solution on output
!! This subroutine is almost similar to precong_per
!! with minor differences: 
!! -the convolutions are in the surface BC
!! -the work arrays psifscf and ww are slightly bigger
!! -the work array z is smaller now
!! -the Fourier preconditioner is in the surface version
subroutine precong_slab(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     ncong,cprecr,hx,hy,hz,x)
  use module_base
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

  call allocate_all()

  call prepare_sdc_slab(n1,n3,modul1,modul3,af,bf,cf,ef,hx,hy,hz)

  !   initializes the wavelet scaling coefficients   
  call wscal_init_per(scal,hx,hy,hz,cprecr)
  !b=x
  call vcopy(nvctr_c+7*nvctr_f,x(1),1,b(1),1) 

  !   compute the input guess x via a Fourier transform in a cubic box.
  !   Arrays psifscf and ww serve as work arrays for the Fourier
  call prec_fft_slab_fast(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
       cprecr,hx,hy,hz,x,&
       psifscf(1),psifscf(n1+2),ww(1),ww(2*((n1+1)/2+1)*(n2+1)*(n3+1)+1))

  !   call apply_hp_slab(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
  !     cprecr,hx,hy,hz,x,d,psifscf,ww) ! d:=Ax
  call apply_hp_slab_sd(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
       cprecr,x,d,psifscf,ww,modul1,modul3,af,bf,cf,ef) ! d:=Ax
  r=b-d
 
  call wscal_per(nvctr_c,nvctr_f,scal,r(1),r(nvctr_c+1),d(1),d(nvctr_c+1))
  !rmr=dot_product(r,d)
  rmr=dot(nvctr_c+7*nvctr_f,r(1),1,d(1),1)
  do i=1,ncong 
     !write(*,*)i,rmr

     !      call apply_hp_slab(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     !        cprecr,hx,hy,hz,d,b,psifscf,ww) ! b:=Ad
     call apply_hp_slab_sd(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
          cprecr,d,b,psifscf,ww,modul1,modul3,af,bf,cf,ef) ! b:=Ad

     !alpha=rmr/dot_product(d,b)
     alpha=rmr/dot(nvctr_c+7*nvctr_f,d(1),1,b(1),1)
     x=x+alpha*d
     r=r-alpha*b

     call wscal_per(nvctr_c,nvctr_f,scal,r(1),r(nvctr_c+1),b(1),b(nvctr_c+1))
     !rmr_new=dot_product(r,b)
     rmr_new=dot(nvctr_c+7*nvctr_f,r(1),1,b(1),1)
     beta=rmr_new/rmr
     d=b+beta*d
     rmr=rmr_new
  enddo

  call deallocate_all()

contains

  subroutine allocate_all

    modul1 = f_malloc(lowfil.to.n1+lupfil,id='modul1')
    modul3 = f_malloc(lowfil.to.n3+lupfil,id='modul3')
    af = f_malloc((/ lowfil.to.lupfil, 1.to.3 /),id='af')
    bf = f_malloc((/ lowfil.to.lupfil, 1.to.3 /),id='bf')
    cf = f_malloc((/ lowfil.to.lupfil, 1.to.3 /),id='cf')
    ef = f_malloc((/ lowfil.to.lupfil, 1.to.3 /),id='ef')
    b = f_malloc(nvctr_c+7*nvctr_f,id='b')
    r = f_malloc(nvctr_c+7*nvctr_f,id='r')
    d = f_malloc(nvctr_c+7*nvctr_f,id='d')
    psifscf = f_malloc((2*n1+2)*(2*n2+16)*(2*n3+2),id='psifscf')
    ww = f_malloc((2*n1+2)*(2*n2+16)*(2*n3+2),id='ww')
  END SUBROUTINE allocate_all

  subroutine deallocate_all

    call f_free(modul1)
    call f_free(modul3)
    call f_free(af)
    call f_free(bf)
    call f_free(cf)
    call f_free(ef)
    call f_free(psifscf)
    call f_free(ww)
    call f_free(b)
    call f_free(r)
    call f_free(d)
  END SUBROUTINE deallocate_all

END SUBROUTINE precong_slab


!>   Applies the operator (KE+cprecr*I)*x=y
!!   array x is input, array y is output
subroutine apply_hp_slab(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,x,y,psifscf,ww)
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
END SUBROUTINE apply_hp_slab


!>   Applies the operator (KE+cprecr*I)*x=y
!!   array x is input, array y is output
subroutine apply_hp_slab_sd(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,x,y,psifscf,ww,modul1,modul3,a,b,c,e)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(in) ::  x(nvctr_c+7*nvctr_f)  
  real(wp), intent(out) ::  y(nvctr_c+7*nvctr_f)
  integer, parameter :: lowfil=-14,lupfil=14

  !n(c) real(gp) hgrid(3)   
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
       nseg_f,nvctr_f,keyg(1,nseg_c+min(1,nseg_f)),keyv(nseg_c+min(1,nseg_f)),   &
       x(1),x(nvctr_c+min(1,nvctr_f)),psifscf)

  !n(c) hgrid(1)=hx
  !n(c) hgrid(2)=hy
  !n(c) hgrid(3)=hz
  ! psifscf: input, ww: output
  !     call convolut_kinetic_slab_c(2*n1+1,2*n2+15,2*n3+1,hgridh,psifscf,ww,cprecr)
  call convolut_kinetic_slab_sdc(n1,n2,n3,psifscf,ww,cprecr,modul1,modul3,a,b,c,e)

  ! ww:intput
  ! y:output
  call compress_sd(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),& 
       nseg_f,nvctr_f,keyg(1,nseg_c+min(1,nseg_f)),keyv(nseg_c+min(1,nseg_f)),& 
       ww,y(1),y(nvctr_c+min(1,nvctr_f)))
END SUBROUTINE apply_hp_slab_sd


!>   Applies the operator (KE+cprecr*I)*x=y
!!   array x is input, array y is output
subroutine apply_hp_slab_sd_scal(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,x,y,psifscf,ww,modul1,modul3,a,b,c,e,scal)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp),intent(in) :: scal(0:7)
  real(gp), intent(in) :: cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(in) ::  x(nvctr_c+7*nvctr_f)  
  real(wp), intent(out) ::  y(nvctr_c+7*nvctr_f)
  integer, parameter :: lowfil=-14,lupfil=14

  !n(c) real(gp) hgrid(3)   
  real(wp),dimension((2*n1+2)*(2*n2+16)*(2*n3+2))::ww,psifscf

  integer,intent(in)::modul1(lowfil:n1+lupfil)
  integer,intent(in)::modul3(lowfil:n3+lupfil)
  real(gp),intent(in)::a(lowfil:lupfil,3)
  real(gp),intent(in)::b(lowfil:lupfil,3)
  real(gp),intent(in)::c(lowfil:lupfil,3)
  real(gp),intent(in)::e(lowfil:lupfil,3)
  ! x: input
  ! psifscf: output
  call uncompress_sd_scal(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
       nseg_f,nvctr_f,keyg(1,nseg_c+min(1,nseg_f)),keyv(nseg_c+min(1,nseg_f)),&
       x(1),x(nvctr_c+min(1,nvctr_f)),psifscf,scal)

  !n(c) hgrid(1)=hx
  !n(c) hgrid(2)=hy
  !n(c) hgrid(3)=hz
  ! psifscf: input, ww: output
  !     call convolut_kinetic_slab_c(2*n1+1,2*n2+15,2*n3+1,hgridh,psifscf,ww,cprecr)
  call convolut_kinetic_slab_sdc(n1,n2,n3,psifscf,ww,cprecr,modul1,modul3,a,b,c,e)
  ! ww:intput
  ! y:output
  call compress_sd_scal(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),& 
       nseg_f,nvctr_f,keyg(1,nseg_c+min(1,nseg_f)),keyv(nseg_c+min(1,nseg_f)),& 
       ww,y(1),y(nvctr_c+min(1,nvctr_f)),scal)

END SUBROUTINE apply_hp_slab_sd_scal


!>   Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
!!   hpsi is the right hand side on input and the solution on output
subroutine prec_fft_slab_fast(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,hpsi,kern_k1,kern_k3,z,x_c)
  use module_base
  implicit none 
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(inout) ::  hpsi(nvctr_c+7*nvctr_f)
  !local variables
  real(gp) :: kern_k1(0:n1),kern_k3(0:n3)
  real(wp) :: x_c(0:n1,0:n2,0:n3)! in and out of Fourier preconditioning
  real(wp) :: z(2,0:(n1+1)/2,0:n2,0:n3)! work array for FFT

  ! diagonally precondition the wavelet part  
  if (nvctr_f > 0) then
     call wscal_f(nvctr_f,hpsi(nvctr_c+1),hx,hy,hz,cprecr)
  end if

  call make_kernel(n1,hx,kern_k1)

  call make_kernel(n3,hz,kern_k3)

  call uncompress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  !   solve the helmholtz equation for the scfunction part  
  call hit_with_kernel_slab(x_c,z,kern_k1,kern_k3,n1,n2,n3,cprecr,hy)   

  call compress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

END SUBROUTINE prec_fft_slab_fast


!>   Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
!!   hpsi is the right hand side on input and the solution on output
subroutine prec_fft_slab(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,hpsi)
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

  !   solve the helmholtz equation for the scfunction part  
  call hit_with_kernel_slab(x_c,z,kern_k1,kern_k3,n1,n2,n3,cprecr,hy)   

  call   compress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  call deallocate_all

contains

  subroutine allocate_all
    kern_k1 = f_malloc(0.to.n1,id='kern_k1')
    kern_k3 = f_malloc(0.to.n3,id='kern_k3')
    z = f_malloc((/ 1.to.2, 0.to.(n1+1)/2, 0.to.n2, 0.to.n3 /),id='z')
    x_c = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='x_c')
  END SUBROUTINE allocate_all

  subroutine deallocate_all
    call f_free(z)
    call f_free(kern_k1)
    call f_free(kern_k3)
    call f_free(x_c)
  END SUBROUTINE deallocate_all

END SUBROUTINE prec_fft_slab


!>   Solve the discretized equation
!!   (-d^2/dy^2+ct(k1,k3)) zx(output) = zx(input)
!!   for all k1,k3 via Lapack
subroutine segment_invert(n1,n2,n3,kern_k1,kern_k3,c,zx,hgrid)
  use module_base
  implicit none
  integer,intent(in) :: n1,n2,n3
  real(wp),intent(in) :: kern_k1(0:n1)
  real(wp),intent(in) :: kern_k3(0:n3)
  real(gp),intent(in) :: c,hgrid
  real(wp),intent(inout) :: zx(2,0:(n1+1)/2,0:n2,0:n3)

  real(gp) :: ct
  real(gp),allocatable,dimension(:,:) :: ab
  !real(wp) :: b(0:n2,2)
  real(wp),allocatable,dimension(:,:) :: b
  !     .. Scalar Arguments ..
  INTEGER :: INFO, Kd, LDAB, LDB, NRHS=2,n
  integer :: i1,i2,i3,i,j

  integer,parameter :: lowfil=-14,lupfil=14
  real(gp) :: scale
  real(gp) :: fil(lowfil:lupfil)

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

  ! hit the fourier transform of x with the kernel

  !$omp parallel default(none) & 
  !$omp private (b,ab,i3,i1,i2,j,i,ct,info,i_stat,i_all) &
  !$omp shared (n1,n2,n3,zx,fil,kd,ldb,ldab,nrhs,n,c,kern_k1,kern_k3)
  !$omp critical (allocate_critical)
  ab = f_malloc((/ ldab, n /),id='ab')
  !b = f_malloc((/ 0.to.n2, 1.to.2 /),id='b')
  allocate(b(0:n2,1:2),stat=i_stat)
  call memocc(i_stat,b,'b','segment_invert')
  !$omp end critical (allocate_critical)
  !$omp do schedule(static,1)
  do i3=0,n3
     !   do i1=0,n1
     do i1=0,(n1+1)/2
        !      ct=kern_k1(i1)+kern_k3(i3)+c
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
        !      call DGBSV( N, KL, KU, NRHS, ab , LDAB, IPIV, b, LDB, INFO )
        call DPBSV( 'U', N, KD, NRHS, AB, LDAB, B, LDB, INFO )
        if (info.ne.0) stop 'error in matrix inversion'
        do i2=0,n2
           zx(1,i1,i2,i3)=b(i2,1)
           zx(2,i1,i2,i3)=b(i2,2)
        enddo
     enddo
  enddo
  !$omp end do
  !$omp critical (deallocate_critical)
  call f_free(ab)
  !call f_free(b)
  i_all = -product(shape(b))*kind(b)
  deallocate(b,stat=i_stat)
  call memocc(i_stat,i_all,'b','segment_invert')
  !$omp end critical (deallocate_critical)
  !$omp end parallel

END SUBROUTINE segment_invert
