!> @file
!!   Routines to do preconditioning on wavefunctions
!! @author
!!    Copyright (C) 2010-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
!! x is the right hand side on input and the solution on output
subroutine precong_per_hyb(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     ncong,cprecr,hx,hy,hz,x,ibyz,ibxz,ibxy)
  use module_base
  implicit none
integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,ncong
  integer,intent(in)::ibyz(2,0:n2,0:n3+ndebug),ibxz(2,0:n1,0:n3+ndebug),ibxy(2,0:n1,0:n2+ndebug)
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

  ! work arrays for adaptive wavelet data structure
  ! x_c and y_c are taken from the FFT arrays
  real(wp),allocatable::x_f(:,:,:,:)
  real(wp),allocatable,dimension(:)::x_f1,x_f2,x_f3
  real(wp),allocatable,dimension(:,:,:,:)::y_f

  ! work arrays for FFT
  real(wp), dimension(:), allocatable :: kern_k1,kern_k2,kern_k3
  real(wp), dimension(:,:,:), allocatable :: x_c! in and out of Fourier preconditioning
  real(wp), dimension(:,:,:,:,:), allocatable::z1,z3 ! work array for FFT

  integer :: nd1,nd2,nd3
  integer :: n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b
  integer :: nf

  call dimensions_fft(n1,n2,n3,nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b)

  nf=(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)

  call allocate_all

  ! initializes the wavelet scaling coefficients 
  call wscal_init_per(scal,hx,hy,hz,cprecr)
  !b=x
  call vcopy(nvctr_c+7*nvctr_f,x(1),1,b(1),1) 

  ! compute the input guess x via a Fourier transform in a cubic box.
  ! Arrays psifscf and ww serve as work arrays for the Fourier
  !        prec_fft_fast(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
  !    cprecr,hx,hy,hz,hpsi,&
  !  kern_k1,kern_k2,kern_k3,z1,z3,x_c,&
  !  nd1,nd2,nd3,n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b)
  call prec_fft_fast(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
       cprecr,hx,hy,hz,x,&
       kern_k1,kern_k2,kern_k3,z1,z3,x_c,&
       nd1,nd2,nd3,n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b)

  call apply_hp_hyb(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
       cprecr,hx,hy,hz,x,d,x_f,x_c,x_f1,x_f2,x_f3,y_f,z1,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,nf,ibyz,ibxz,ibxy)

  r=b-d

  call wscal_per(nvctr_c,nvctr_f,scal,r(1),r(nvctr_c+1),d(1),d(nvctr_c+1))
  !rmr=dot_product(r,d)
  rmr=dot(nvctr_c+7*nvctr_f,r(1),1,d(1),1)
  do i=1,ncong
     !write(*,*)i,rmr
     !  write(*,*)i,sqrt(rmr)

     call apply_hp_hyb(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
          cprecr,hx,hy,hz,d,b,x_f,x_c,x_f1,x_f2,x_f3,y_f,z1,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,nf,ibyz,ibxz,ibxy)

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
    allocate(b(nvctr_c+7*nvctr_f+ndebug),stat=i_stat)
    call memocc(i_stat,b,'b','precong_per')
    allocate(r(nvctr_c+7*nvctr_f+ndebug),stat=i_stat)
    call memocc(i_stat,r,'r','precong_per')
    allocate(d(nvctr_c+7*nvctr_f+ndebug),stat=i_stat)
    call memocc(i_stat,d,'','precong_per')

    allocate(kern_k1(0:n1+ndebug),stat=i_stat)
    call memocc(i_stat,kern_k1,'kern_k1','precong_per')
    allocate(kern_k2(0:n2+ndebug),stat=i_stat)
    call memocc(i_stat,kern_k2,'kern_k2','precong_per')
    allocate(kern_k3(0:n3+ndebug),stat=i_stat)
    call memocc(i_stat,kern_k3,'kern_k3','precong_per')
    allocate(z1(2,nd1b,nd2,nd3,2+ndebug),stat=i_stat) ! work array for fft
    call memocc(i_stat,z1,'z1','precong_per')
    allocate(z3(2,nd1,nd2,nd3f,2+ndebug),stat=i_stat) ! work array for fft
    call memocc(i_stat,z3,'z3','precong_per')
    allocate(x_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
    call memocc(i_stat,x_c,'x_c','precong_per')

    allocate(x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
    call memocc(i_stat,x_f,'x_f','precong_per')
    allocate(x_f1(nf+ndebug),stat=i_stat)
    call memocc(i_stat,x_f1,'x_f1','precong_per')
    allocate(x_f2(nf+ndebug),stat=i_stat)
    call memocc(i_stat,x_f2,'x_f2','precong_per')
    allocate(x_f3(nf+ndebug),stat=i_stat)
    call memocc(i_stat,x_f3,'x_f3','precong_per')
    allocate(y_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
    call memocc(i_stat,y_f,'y_f','precong_per')

  END SUBROUTINE allocate_all

  subroutine deallocate_all

    i_all=-product(shape(b))*kind(b)
    deallocate(b,stat=i_stat)
    call memocc(i_stat,i_all,'b','precong_per')

    i_all=-product(shape(r))*kind(r)
    deallocate(r,stat=i_stat)
    call memocc(i_stat,i_all,'r','precong_per')

    i_all=-product(shape(d))*kind(d)
    deallocate(d,stat=i_stat)
    call memocc(i_stat,i_all,'d','precong_per')

    i_all=-product(shape(z1))*kind(z1)
    deallocate(z1,stat=i_stat)
    call memocc(i_stat,i_all,'z1','prec_fft')

    i_all=-product(shape(z3))*kind(z3)
    deallocate(z3,stat=i_stat)
    call memocc(i_stat,i_all,'z3','prec_fft')

    i_all=-product(shape(kern_k1))*kind(kern_k1)
    deallocate(kern_k1,stat=i_stat)
    call memocc(i_stat,i_all,'kern_k1','prec_fft')

    i_all=-product(shape(kern_k2))*kind(kern_k2)
    deallocate(kern_k2,stat=i_stat)
    call memocc(i_stat,i_all,'kern_k2','prec_fft')

    i_all=-product(shape(kern_k3))*kind(kern_k3)
    deallocate(kern_k3,stat=i_stat)
    call memocc(i_stat,i_all,'kern_k3','prec_fft')

    i_all=-product(shape(x_c))*kind(x_c)
    deallocate(x_c,stat=i_stat)
    call memocc(i_stat,i_all,'x_c','prec_fft')

    i_all=-product(shape(x_f))*kind(x_f)
    deallocate(x_f,stat=i_stat)
    call memocc(i_stat,i_all,'x_f','precong_per')

    i_all=-product(shape(x_f1))*kind(x_f1)
    deallocate(x_f1,stat=i_stat)
    call memocc(i_stat,i_all,'x_f1','precong_per')

    i_all=-product(shape(x_f2))*kind(x_f2)
    deallocate(x_f2,stat=i_stat)
    call memocc(i_stat,i_all,'x_f2','precong_per')

    i_all=-product(shape(x_f3))*kind(x_f3)
    deallocate(x_f3,stat=i_stat)
    call memocc(i_stat,i_all,'x_f3','precong_per')

    i_all=-product(shape(y_f))*kind(y_f)
    deallocate(y_f,stat=i_stat)
    call memocc(i_stat,i_all,'y_f','precong_per')

  END SUBROUTINE deallocate_all

END SUBROUTINE precong_per_hyb


!>   Applies the operator (KE+cprecr*I)*x=y
!!   array x is input, array y is output
subroutine apply_hp_hyb(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,x,y,x_f,x_c,x_f1,x_f2,x_f3,y_f,y_c,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,nf,ibyz,ibxz,ibxy)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer,intent(in) :: ibyz(2,0:n2,0:n3+ndebug),ibxz(2,0:n1,0:n3+ndebug),ibxy(2,0:n1,0:n2+ndebug)
  integer,intent(in) :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,nf
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(in) :: x(nvctr_c+7*nvctr_f)  
  real(wp), intent(out) :: y(nvctr_c+7*nvctr_f)
  !work arrays 
  real(wp) :: x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
  real(wp) :: x_c(0:n1,0:n2,0:n3)
  real(wp) :: x_f1(nf),x_f2(nf),x_f3(nf)
  real(wp) :: y_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
  real(wp) :: y_c(0:n1,0:n2,0:n3)

  real(gp) :: hgrid(3)

  call uncompress_per_f(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
       nseg_f,nvctr_f,keyg(1,nseg_c+min(1,nseg_f)),keyv(nseg_c+min(1,nseg_f)),   &
       x(1),x(nvctr_c+min(1,nvctr_f)),x_c,x_f,x_f1,x_f2,x_f3,&
       nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)

  hgrid(1)=hx
  hgrid(2)=hy
  hgrid(3)=hz

  call convolut_kinetic_hyb_c(n1,n2,n3, &
       nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       hgrid,x_c,x_f,y_c,y_f,cprecr,x_f1,x_f2,x_f3,ibyz,ibxz,ibxy)

  call compress_per_f(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   & 
       nseg_f,nvctr_f,keyg(1,nseg_c+min(1,nseg_f)),keyv(nseg_c+min(1,nseg_f)), & 
       y_c,y_f,y(1),y(nvctr_c+min(1,nvctr_f)),nfl1,nfl2,nfl3,nfu1,nfu2,nfu3)

END SUBROUTINE apply_hp_hyb

