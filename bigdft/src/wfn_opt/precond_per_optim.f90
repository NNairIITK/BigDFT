!> @file
!!  Optimized routines to precondition wavefunctions
!! @author
!!    Copyright (C) 2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
 

!>   Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
!!   x is the right hand side on input and the solution on output
subroutine precong_per(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
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
  integer, parameter :: lowfil=-14,lupfil=14
  real(gp) :: scal(0:7),fac
  real(wp) :: rmr_old,rmr_new,alpha,beta
  integer :: i,icong
  real(wp), allocatable :: b(:),r(:),d(:)
  real(wp), allocatable :: psifscf(:),ww(:)
  integer :: nd1,nd2,nd3
  integer :: n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b
  real(gp), allocatable, dimension(:,:) :: af,bf,cf,ef
  integer, allocatable, dimension(:) :: modul1,modul2,modul3

  ! Array sizes for the real-to-complex FFT: note that n1(there)=n1(here)+1
  ! and the same for n2,n3.
  call dimensions_fft(n1,n2,n3,nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b)

  call allocate_all

  call prepare_sdc(n1,n2,n3,modul1,modul2,modul3,af,bf,cf,ef,hx,hy,hz)
  ! initializes the wavelet scaling coefficients 
  call wscal_init_per(scal,hx,hy,hz,cprecr)

  ! scale the r.h.s. that is also the scaled input guess :
  ! b'=D^{-1/2}b
  call wscal_per_self(nvctr_c,nvctr_f,scal,x(1),x(nvctr_c+1))
  !b=x
  call vcopy(nvctr_c+7*nvctr_f,x(1),1,b(1),1) 

  ! compute the input guess x via a Fourier transform in a cubic box.
  ! Arrays psifscf and ww serve as work arrays for the Fourier
  fac=1.d0/scal(0)**2
  call prec_fft_c(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
       cprecr,hx,hy,hz,x,&
       psifscf(1),psifscf(n1+2),psifscf(n1+n2+3),ww(1),ww(nd1b*nd2*nd3*4+1),&
       ww(nd1b*nd2*nd3*4+nd1*nd2*nd3f*4+1),&
       nd1,nd2,nd3,n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b,fac)

  call apply_hp_scal(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
       cprecr,x,d,psifscf,ww,modul1,modul2,modul3,af,bf,cf,ef,scal) ! d:=Ax

!!  x=d
!!  return

  r=b-d ! r=b-Ax
  d=r
  !rmr_new=dot_product(r,r)
  rmr_new=dot(nvctr_c+7*nvctr_f,r(1),1,r(1),1) !n(?) r not allocated

  do icong=1,ncong 
     !write(*,*)icong,rmr_new

     call apply_hp_scal(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
          cprecr,d,b,psifscf,ww,modul1,modul2,modul3,af,bf,cf,ef,scal) ! b:=Ad

     alpha=rmr_new/dot(nvctr_c+7*nvctr_f,d(1),1,b(1),1)

     do i=1,nvctr_c+7*nvctr_f
        x(i)=x(i)+alpha*d(i)
        r(i)=r(i)-alpha*b(i)
     enddo

     if (icong==ncong) exit

     rmr_old=rmr_new
     rmr_new=dot(nvctr_c+7*nvctr_f,r(1),1,r(1),1)

     beta=rmr_new/rmr_old
     d=r+beta*d
  enddo

  ! x=D^{-1/2}x'
  call wscal_per_self(nvctr_c,nvctr_f,scal,x(1),x(nvctr_c+1))
  ! write(30,*) x
  ! stop

  call deallocate_all

contains

  subroutine allocate_all
    modul1 = f_malloc(lowfil.to.n1+lupfil,id='modul1')
    modul2 = f_malloc(lowfil.to.n2+lupfil,id='modul2')
    modul3 = f_malloc(lowfil.to.n3+lupfil,id='modul3')
    af = f_malloc((/ lowfil.to.lupfil, 1.to.3 /),id='af')
    bf = f_malloc((/ lowfil.to.lupfil, 1.to.3 /),id='bf')
    cf = f_malloc((/ lowfil.to.lupfil, 1.to.3 /),id='cf')
    ef = f_malloc((/ lowfil.to.lupfil, 1.to.3 /),id='ef')
    b = f_malloc(nvctr_c+7*nvctr_f,id='b')
    r = f_malloc(nvctr_c+7*nvctr_f,id='r')
    d = f_malloc(nvctr_c+7*nvctr_f,id='d')
    psifscf = f_malloc((2*n1+2)*(2*n2+2)*(2*n3+2),id='psifscf')
    ww = f_malloc((2*n1+2)*(2*n2+2)*(2*n3+2),id='ww')
  END SUBROUTINE allocate_all

  subroutine deallocate_all

    call f_free(modul1)
    call f_free(modul2)
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

END SUBROUTINE precong_per


!>   Solves (KE+cprecr*I)*xx=yy by FFT in a cubic box 
!!   x_c is the right hand side on input and the solution on output
!!   This version uses work arrays kern_k1-kern_k3 and z allocated elsewhere
subroutine prec_fft_c(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,hpsi,&
     kern_k1,kern_k2,kern_k3,z1,z3,x_c,&
     nd1,nd2,nd3,n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b,fac)
  use module_base
  implicit none 
  integer, intent(in) :: n1,n2,n3
  integer,intent(in) :: nd1,nd2,nd3
  integer,intent(in) :: n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr,fac
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(inout) :: hpsi(nvctr_c+7*nvctr_f) 

  !work arrays
  real(gp):: kern_k1(0:n1),kern_k2(0:n2),kern_k3(0:n3)
  real(wp),dimension(0:n1,0:n2,0:n3):: x_c! in and out of Fourier preconditioning
  real(wp)::z1(2,nd1b,nd2,nd3,2)! work array
  real(wp)::z3(2,nd1,nd2,nd3f,2)! work array

  call make_kernel(n1,hx,kern_k1)
  call make_kernel(n2,hy,kern_k2)
  call make_kernel(n3,hz,kern_k3)

  call uncompress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  call hit_with_kernel_fac(x_c,z1,z3,kern_k1,kern_k2,kern_k3,n1+1,n2+1,n3+1,nd1,nd2,nd3,&
       n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b,cprecr,fac)

  call compress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

END SUBROUTINE prec_fft_c


!>   Solves (KE+cprecr*I)*xx=yy by FFT in a cubic box 
!!   x_c is the right hand side on input and the solution on output
!!   This version uses work arrays kern_k1-kern_k3 and z allocated elsewhere
subroutine prec_fft_fast(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,hpsi,&
     kern_k1,kern_k2,kern_k3,z1,z3,x_c,&
     nd1,nd2,nd3,n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b)
  use module_base
  implicit none 
  integer, intent(in) :: n1,n2,n3
  integer,intent(in)::nd1,nd2,nd3
  integer,intent(in)::n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(inout) ::  hpsi(nvctr_c+7*nvctr_f) 

  !work arrays
  real(gp):: kern_k1(0:n1),kern_k2(0:n2),kern_k3(0:n3)
  real(wp),dimension(0:n1,0:n2,0:n3):: x_c! in and out of Fourier preconditioning
  real(wp)::z1(2,nd1b,nd2,nd3,2)! work array
  real(wp)::z3(2,nd1,nd2,nd3f,2)! work array

  if (nvctr_f > 0) then
     call wscal_f(nvctr_f,hpsi(nvctr_c+1),hx,hy,hz,cprecr)
  end if

  call make_kernel(n1,hx,kern_k1)
  call make_kernel(n2,hy,kern_k2)
  call make_kernel(n3,hz,kern_k3)

  call uncompress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  call  hit_with_kernel(x_c,z1,z3,kern_k1,kern_k2,kern_k3,n1+1,n2+1,n3+1,nd1,nd2,nd3,&
       n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b,cprecr)

  call compress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

END SUBROUTINE prec_fft_fast


!> Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
!! hpsi is the right hand side on input and the solution on output
subroutine prec_fft(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,hpsi)
  use module_base
  implicit none 
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(inout) :: hpsi(nvctr_c+7*nvctr_f)
  !local variables
  character(len=*), parameter :: subname='prec_fft'
  integer :: nd1,nd2,nd3
  real(gp), dimension(:), allocatable :: kern_k1,kern_k2,kern_k3
  real(wp), dimension(:,:,:), allocatable :: x_c! in and out of Fourier preconditioning
  real(wp), dimension(:,:,:,:,:), allocatable :: z1,z3 ! work array for FFT
  integer :: n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b

  ! Array sizes for the real-to-complex FFT: note that n1(there)=n1(here)+1
  ! and the same for n2,n3.
  call dimensions_fft(n1,n2,n3,nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b)

  call allocate_all

  ! diagonally precondition the wavelet part  
  if (nvctr_f > 0) then
     call wscal_f(nvctr_f,hpsi(nvctr_c+1),hx,hy,hz,cprecr)
  end if

  call make_kernel(n1,hx,kern_k1)
  call make_kernel(n2,hy,kern_k2)
  call make_kernel(n3,hz,kern_k3)

  call uncompress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  ! solve the helmholtz equation for the scfunction part  
  call  hit_with_kernel(x_c,z1,z3,kern_k1,kern_k2,kern_k3,n1+1,n2+1,n3+1,nd1,nd2,nd3,&
       n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b,cprecr)

  call compress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  call deallocate_all

contains

  subroutine allocate_all
    kern_k1 = f_malloc(0.to.n1,id='kern_k1')
    kern_k2 = f_malloc(0.to.n2,id='kern_k2')
    kern_k3 = f_malloc(0.to.n3,id='kern_k3')
    z1 = f_malloc((/ 2, nd1b, nd2, nd3, 2 /),id='z1')
    z3 = f_malloc((/ 2, nd1, nd2, nd3f, 2 /),id='z3')
    x_c = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='x_c')
  END SUBROUTINE allocate_all

  subroutine deallocate_all
    call f_free(z1)
    call f_free(z3)
    call f_free(kern_k1)
    call f_free(kern_k2)
    call f_free(kern_k3)
    call f_free(x_c)
  END SUBROUTINE deallocate_all

END SUBROUTINE prec_fft


!>  Applies the operator (KE+cprecr*I)*x=y
!!  array x is input, array y is output
subroutine apply_hp(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,x,y,psifscf,ww)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(in) :: x(nvctr_c+7*nvctr_f)  
  real(wp), intent(out) :: y(nvctr_c+7*nvctr_f)

  real(gp) :: hgridh(3)
  real(wp), dimension((2*n1+2)*(2*n2+2)*(2*n3+2)) :: ww,psifscf

  call uncompress_per(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
       x(1),x(nvctr_c+1),psifscf,ww)

  hgridh(1)=hx*.5_gp
  hgridh(2)=hy*.5_gp
  hgridh(3)=hz*.5_gp
  call convolut_kinetic_per_c(2*n1+1,2*n2+1,2*n3+1,hgridh,psifscf,ww,cprecr)

  call compress_per(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   & 
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   & 
       ww,y(1),y(nvctr_c+1),psifscf)
END SUBROUTINE apply_hp


!> Applies the operator (KE+cprecr*I)*x=y
!! array x is input, array y is output
subroutine apply_hp_slab_k(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,kx,ky,kz,x,y,psifscf,ww,scal)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr,kx,ky,kz
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), dimension(0:7), intent(in) :: scal
  real(wp), intent(in) ::  x(nvctr_c+7*nvctr_f,2)  
  real(wp), intent(out) ::  y(nvctr_c+7*nvctr_f,2)
  real(wp), dimension((2*n1+2)*(2*n2+16)*(2*n3+2),2) :: ww,psifscf
  !local variables
  integer :: idx
  real(gp) :: hgridh(3)

  ! x: input, ww:work
  ! psifscf: output
  do idx=1,2
     call uncompress_slab_scal(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
          nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
          x(1,idx),x(nvctr_c+1,idx),psifscf(1,idx),ww(1,idx),scal)
  end do

  hgridh(1)=hx*.5_gp
  hgridh(2)=hy*.5_gp
  hgridh(3)=hz*.5_gp
  !transpose (to be included in the uncompression)
  call transpose_for_kpoints(2,2*n1+2,2*n2+16,2*n3+2,&
       psifscf,ww,.true.)

  ! psifscf: input, ww: output
  call convolut_kinetic_slab_c_k(2*n1+1,2*n2+15,2*n3+1,hgridh,psifscf,ww,cprecr,&
       kx,ky,kz)

  call transpose_for_kpoints(2,2*n1+2,2*n2+16,2*n3+2,&
       ww,psifscf,.false.)

  ! ww:intput, psifscf: work
  ! y:output
  do idx=1,2
     call compress_slab_scal(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),& 
          nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),& 
          ww(1,idx),y(1,idx),y(nvctr_c+1,idx),psifscf(1,idx),scal)
  end do

END SUBROUTINE apply_hp_slab_k


!> Applies the operator (KE+cprecr*I)*x=y
!! array x is input, array y is output
subroutine apply_hp_per_k(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,k1,k2,k3,x,y,psifscf,ww,scal)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr,k1,k2,k3
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), dimension(0:7), intent(in) :: scal
  real(wp), dimension(nvctr_c+7*nvctr_f,2), intent(in) :: x
  real(wp), dimension(nvctr_c+7*nvctr_f,2), intent(in) :: y
  real(wp), dimension((2*n1+2)*(2*n2+2)*(2*n3+2),2) :: ww,psifscf
  !local variables
  logical, parameter :: transpose=.false.
  integer :: idx
  real(gp) :: hgridh(3)
  
  do idx=1,2
     call uncompress_per_scal(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
          nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
          x(1,idx),x(nvctr_c+1,idx),psifscf(1,idx),ww(1,idx),scal)
  end do

  hgridh(1)=hx*.5_gp
  hgridh(2)=hy*.5_gp
  hgridh(3)=hz*.5_gp

  if (transpose) then
     !transpose (to be included in the uncompression)
     call transpose_for_kpoints(2,2*n1+2,2*n2+2,2*n3+2,&
          psifscf,ww,.true.)
     call convolut_kinetic_per_c_k(2*n1+1,2*n2+1,2*n3+1,hgridh,psifscf,ww,cprecr,k1,k2,k3)
     !call convolut_kinetic_per_c(2*n1+1,2*n2+1,2*n3+1,hgridh,psifscf,ww,cprecr)
     !this can be included in the compression
     call transpose_for_kpoints(2,2*n1+2,2*n2+2,2*n3+2,&
          ww,psifscf,.false.)
  else
     call convolut_kinetic_per_c_k_notranspose(2*n1+1,2*n2+1,2*n3+1,hgridh,psifscf,ww,cprecr,k1,k2,k3)
  end if

  do idx=1,2
     call compress_per_scal(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   & 
          nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   & 
          ww(1,idx),y(1,idx),y(nvctr_c+1,idx),psifscf(1,idx),scal)
  end do
END SUBROUTINE apply_hp_per_k


!> multiplies a wavefunction psi_c,psi_f (in vector form) with a scaling vector (scal)
subroutine wscal_f(mvctr_f,psi_f,hx,hy,hz,c)
  use module_base
  implicit none
  integer,intent(in)::mvctr_f
  real(gp),intent(in)::c,hx,hy,hz

  real(wp)::psi_f(7,mvctr_f)
  real(gp)::scal(7),hh(3)
  !WAVELET AND SCALING FUNCTION SECOND DERIVATIVE FILTERS, diagonal elements
  real(gp),PARAMETER::B2=24.8758460293923314_gp,A2=3.55369228991319019_gp

  integer i

  hh(1)=.5_gp/hx**2
  hh(2)=.5_gp/hy**2
  hh(3)=.5_gp/hz**2

  scal(1)=1._gp/(b2*hh(1)+a2*hh(2)+a2*hh(3)+c)       !  2 1 1
  scal(2)=1._gp/(a2*hh(1)+b2*hh(2)+a2*hh(3)+c)       !  1 2 1
  scal(3)=1._gp/(b2*hh(1)+b2*hh(2)+a2*hh(3)+c)       !  2 2 1
  scal(4)=1._gp/(a2*hh(1)+a2*hh(2)+b2*hh(3)+c)       !  1 1 2
  scal(5)=1._gp/(b2*hh(1)+a2*hh(2)+b2*hh(3)+c)       !  2 1 2
  scal(6)=1._gp/(a2*hh(1)+b2*hh(2)+b2*hh(3)+c)       !  1 2 2
  scal(7)=1._gp/(b2*hh(1)+b2*hh(2)+b2*hh(3)+c)       !  2 2 2

  do i=1,mvctr_f
     psi_f(1,i)=psi_f(1,i)*scal(1)       !  2 1 1
     psi_f(2,i)=psi_f(2,i)*scal(2)       !  1 2 1
     psi_f(3,i)=psi_f(3,i)*scal(3)       !  2 2 1
     psi_f(4,i)=psi_f(4,i)*scal(4)       !  1 1 2
     psi_f(5,i)=psi_f(5,i)*scal(5)       !  2 1 2
     psi_f(6,i)=psi_f(6,i)*scal(6)       !  1 2 2
     psi_f(7,i)=psi_f(7,i)*scal(7)       !  2 2 2
  enddo

END SUBROUTINE wscal_f


!> multiplies a wavefunction psi_c,psi_f (in vector form) with a scaling vector (scal)
subroutine wscal_per_self(mvctr_c,mvctr_f,scal,psi_c,psi_f)
  use module_base
  implicit none
  integer,intent(in)::mvctr_c,mvctr_f
  real(gp),intent(in)::scal(0:7)
  real(wp),intent(inout)::psi_c(mvctr_c),psi_f(7,mvctr_f)

  integer i

  do i=1,mvctr_c
     psi_c(i)=psi_c(i)*scal(0)           !  1 1 1
  enddo

  do i=1,mvctr_f
     psi_f(1,i)=psi_f(1,i)*scal(1)       !  2 1 1
     psi_f(2,i)=psi_f(2,i)*scal(2)       !  1 2 1
     psi_f(3,i)=psi_f(3,i)*scal(3)       !  2 2 1
     psi_f(4,i)=psi_f(4,i)*scal(4)       !  1 1 2
     psi_f(5,i)=psi_f(5,i)*scal(5)       !  2 1 2
     psi_f(6,i)=psi_f(6,i)*scal(6)       !  1 2 2
     psi_f(7,i)=psi_f(7,i)*scal(7)       !  2 2 2
  enddo

END SUBROUTINE wscal_per_self


!> multiplies a wavefunction psi_c,psi_f (in vector form) with a scaling vector (scal)
subroutine wscal_per(mvctr_c,mvctr_f,scal,psi_c_in,psi_f_in,psi_c_out,psi_f_out)
  use module_base
  implicit none
  integer,intent(in)::mvctr_c,mvctr_f
  real(gp),intent(in)::scal(0:7)
  real(wp),intent(in)::psi_c_in(mvctr_c),psi_f_in(7,mvctr_f)
  real(wp),intent(out)::psi_c_out(mvctr_c),psi_f_out(7,mvctr_f)

  integer i

  do i=1,mvctr_c
     psi_c_out(i)=psi_c_in(i)*scal(0)           !  1 1 1
  enddo

  do i=1,mvctr_f
     psi_f_out(1,i)=psi_f_in(1,i)*scal(1)       !  2 1 1
     psi_f_out(2,i)=psi_f_in(2,i)*scal(2)       !  1 2 1
     psi_f_out(3,i)=psi_f_in(3,i)*scal(3)       !  2 2 1
     psi_f_out(4,i)=psi_f_in(4,i)*scal(4)       !  1 1 2
     psi_f_out(5,i)=psi_f_in(5,i)*scal(5)       !  2 1 2
     psi_f_out(6,i)=psi_f_in(6,i)*scal(6)       !  1 2 2
     psi_f_out(7,i)=psi_f_in(7,i)*scal(7)       !  2 2 2
  enddo

END SUBROUTINE wscal_per


!> initialization for the array scal in the subroutine wscal_per  
subroutine wscal_init_per(scal,hx,hy,hz,c)
  use module_base
  implicit none
  real(wp), intent(in) :: c,hx,hy,hz
  real(wp), dimension(0:7), intent(out) :: scal
  !local variables
  real(wp), parameter :: b2=24.8758460293923314d0,a2=3.55369228991319019d0
  real(gp) :: hh(3)

  hh(1)=.5_wp/hx**2
  hh(2)=.5_wp/hy**2
  hh(3)=.5_wp/hz**2

  scal(0)=1._wp/sqrt(a2*hh(1)+a2*hh(2)+a2*hh(3)+c)       !  1 1 1
  scal(1)=1._wp/sqrt(b2*hh(1)+a2*hh(2)+a2*hh(3)+c)       !  2 1 1
  scal(2)=1._wp/sqrt(a2*hh(1)+b2*hh(2)+a2*hh(3)+c)       !  1 2 1
  scal(3)=1._wp/sqrt(b2*hh(1)+b2*hh(2)+a2*hh(3)+c)       !  2 2 1
  scal(4)=1._wp/sqrt(a2*hh(1)+a2*hh(2)+b2*hh(3)+c)       !  1 1 2
  scal(5)=1._wp/sqrt(b2*hh(1)+a2*hh(2)+b2*hh(3)+c)       !  2 1 2
  scal(6)=1._wp/sqrt(a2*hh(1)+b2*hh(2)+b2*hh(3)+c)       !  1 2 2
  scal(7)=1._wp/sqrt(b2*hh(1)+b2*hh(2)+b2*hh(3)+c)       !  2 2 2

END SUBROUTINE wscal_init_per
