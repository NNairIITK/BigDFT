!> @file
!!  Use integral form for Poisson solver
!! @author
!!    Copyright (c) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
program PS_Integral

  use module_base
  use module_types
  use module_interfaces
  use module_xc
  use Poisson_Solver
  use yaml_output
  use dynamic_memory

  implicit none
  integer :: n_points, n_range, n_scf, itype_scf
  integer, dimension(1:7) :: n_points_list = 0
  real(dp) :: hgrid, dx!, t0, t1
  real(dp), dimension(7) :: pgauss
  real(dp), dimension(:), allocatable :: x_scf
  real(dp), dimension(:), allocatable :: y_scf
  !real(dp), dimension(-n_range:n_range) :: work
  !real(dp), dimension(-n_range:n_range) :: kernel_scf
  !real(dp), dimension(-n_range:n_range) :: work
  real(dp), dimension(:), allocatable :: work
  !real(dp), dimension(-n_range:n_range) :: kernel_scf
  real(dp), dimension(:), allocatable :: kernel_scf
  real(dp), dimension(:,:), allocatable :: analytic_vs_naive, gaussian
  real(dp), dimension(:,:,:), allocatable :: multiple_naive, analytic_integral_result, timings
  !local variables
  character(len=*), parameter :: subname='PS_Integral'
  real(dp), parameter :: p0_ref = 1.0_dp
  integer :: i,i_stat,j,k!,n_iter,_kern
  !real(dp) :: p0_cell,p0gauss,absci,kern,moment
  character(len=64) :: chain
  !character(len=*) :: chain
  logical :: timings_switch = .false. 
  real(dp), dimension(0:2048) :: fISF

  integer :: n1,n1_old,n2,n2_old,n3,n3_old,nb1,nb2,nb3,itype,nd,i_all,nrange
  real(dp) :: hx,x,hy,y,hz,z,dy,dz,xgauss,ygauss,theta,sum,sum2,sumx,sumx2,sumy,sumy2,xmin,xmax,ymin,ymax
  real(dp), dimension(:), allocatable :: x_phi, y_phi
  real(dp), dimension(:,:,:), allocatable :: psifscf,psifscfold,psi_w,psi_w2,dx_field,dy_field
  real(dp), external :: lr_gauss

  call f_malloc_set_status(memory_limit=0.e0)

  call f_malloc_routine_id(subname)
!!$  include 'lazy_ISF_8_2048.inc'
!!$  include 'lazy_ISF_14_2048.inc'
!!$  include 'lazy_ISF_16_2048.inc'
!!$  include 'lazy_ISF_20_2048.inc'
!!$  include 'lazy_ISF_24_2048.inc'
!!$  include 'lazy_ISF_30_2048.inc'
!!$  include 'lazy_ISF_40_2048.inc'
!!$  include 'lazy_ISF_50_2048.inc'
!!$  include 'lazy_ISF_60_2048.inc'
!!$  include 'lazy_ISF_100_2048.inc'

!!$  interface 
!!$     subroutine my_analytic_integral(alpha,ntot,m,fwork,fISF,argument_nf)
!!$       use module_base
!!$       implicit none
!!$       integer, intent(in) :: ntot,m
!!$       real(dp), intent(in) :: alpha
!!$       real(dp), dimension(0:ntot), intent(inout) :: fwork
!!$       real(dp), dimension(0:2048), intent(in) :: fISF
!!$       integer, optional, intent(in) :: argument_nf  
!!$     end subroutine my_analytic_integral
!!$  end interface

  !pgauss = 1.0e21_dp
  hgrid = 1.0_dp


  call get_command_argument(1,chain)
  if(trim(chain)=='') then
     write(*,'(1x,a)')&
          'Usage: ./PS_Integral itype_scf [timings]'
     stop
  end if
  read(unit=chain,fmt=*) itype_scf

  call get_command_argument(2,chain)
  if(trim(chain)=='') then
     timings_switch = .false.
  else if (trim(chain)=='timings') then
     timings_switch = .true.
  else 
     write(*,'(1x,a)')&
          'Usage: ./PS_Integral itype_scf [timings]'
     stop
  end if

  !itype_scf = 14

  !Only itype=8,14,16,20,24,30,40,50,60,100
  call four_isf_coeffs(itype_scf,fISF)

!!$  select case(itype_scf)
!!$  case(8)
!!$     fISF => fISF8   
!!$  case(14)
!!$     fISF => fISF14
!!$  case(16)
!!$     fISF => fISF16
!!$  case(20)
!!$     fISF => fISF20    
!!$  case(24)
!!$     fISF => fISF24
!!$  case(30)
!!$     fISF => fISF30
!!$  case(40)
!!$     fISF => fISF40
!!$  case(50)
!!$     fISF => fISF50
!!$  case(60)
!!$     fISF => fISF60
!!$  case(100)
!!$     fISF => fISF100
!!$  case default
!!$     print *,"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100"
!!$     stop
!!$  end select

  n_range=2*itype_scf  

  print *,'here'


  !testing interpolate
hx=0.1
dx=1.0
hy=0.1
  dy=0.0!-0.44
  hz=0.2
  dz=0.0!0.55
n1_old=22
n1=22
n2_old=22
n2=22
  n3_old=8
  n3=8
  nb1=0
  nb2=0
  nb3=0

  itype=16
  nd=2**20

allocate(dx_field(-nb1:2*n1_old+1+nb1,-nb2:2*n2_old+1+nb2,1+ndebug),stat=i_stat)
allocate(dy_field(-nb1:2*n1_old+1+nb1,-nb2:2*n2_old+1+nb2,1+ndebug),stat=i_stat)

dx_field=dx/hx


!1d
   allocate(psifscfold(-nb1:2*n1_old+1+nb1,-nb2:2*n2_old+1+nb2,1+ndebug),stat=i_stat)
  call memocc(i_stat,psifscfold,'psifscfold',subname)
   allocate(psi_w(-nb1:2*n1+1+nb1,-nb2:2*n2_old+1+nb2,1+ndebug),stat=i_stat)
  call memocc(i_stat,psi_w,'psi_w',subname)
   allocate(psi_w2(-nb1:2*n1+1+nb1,-nb2:2*n2+1+nb2,1+ndebug),stat=i_stat)
   call memocc(i_stat,psi_w,'psi_w',subname)
   allocate(psifscf(-nb1:2*n1+1+nb1,-nb2:2*n2_old+1+nb2,1+ndebug),stat=i_stat)
  call memocc(i_stat,psifscf,'psifscf',subname)

!   allocate(psifscfold(-nb1:2*n1_old+1+nb1,-nb2:2*n2_old+1+nb2,-nb3:2*n3_old+1+nb3+ndebug),stat=i_stat)
!   call memocc(i_stat,psifscfold,'psifscfold',subname)
!   allocate(psi_w(-nb1:2*n1+1+nb1,-nb2:2*n2_old+1+nb2,-nb3:2*n3_old+1+nb3+ndebug),stat=i_stat)
!   call memocc(i_stat,psi_w,'psi_w',subname)
!   allocate(psi_w2(-nb1:2*n1+1+nb1,-nb2:2*n2+1+nb2,-nb3:2*n3_old+1+nb3+ndebug),stat=i_stat)
!   call memocc(i_stat,psi_w2,'psi_w2',subname)
!   allocate(psifscf(-nb1:2*n1+1+nb1,-nb2:2*n2+1+nb2,-nb3:2*n3+1+nb3+ndebug),stat=i_stat)
!   call memocc(i_stat,psifscf,'psifscf',subname)

  ! fill psifscfold with a gaussian
  x=-n1_old*hx!0.d0
  do i=0,2*n1_old+1
    y=-n2_old*hy!0.d0
    xgauss=lr_gauss(x,0.0d0,0.7d0)!lr_gauss(x,1.7d0,1.4d0)
    do j=0,2*n2_old+1
       z=0.d0
       !ygauss=lr_gauss(y,0.0d0,0.7d0)
       ygauss=lr_gauss(y,-0.3d0,1.4d0)
       do k=1,1!1,2*n3_old+1
          psifscfold(j,i,k)=xgauss!*ygauss!*lr_gauss(z,1.4d0,2.1d0)
theta=35.0_gp*(4.0_gp*atan(1.d0)/180.0_gp) !in degrees(converted)
dx_field(i,j,k) =x*cos(theta)-y*sin(theta)-x
!dy_field(i,j,k) =x*sin(theta)+y*cos(theta)-y
dy_field(j,i,k) =((1.0d0/cos(theta))-1.0d0)*y+tan(theta)*x
          !write(100,*) x,y,z,psifscfold(i,j,k),dx_field(i,j,k),dy_field(i,j,k)!,psifscf(i,j,k)
          z=z+hz
       end do
       !!write(100,*) ''
       y=y+hy
    end do
    !write(100,*) ''
    x=x+hx
  end do

  ! make an empty border
  do i=0,2*n1_old+1
     psifscfold(i,0,1) = 0.0_gp
     psifscfold(i,2*n2_old,1) = 0.0_gp
     psifscfold(i,2*n2_old+1,1) = 0.0_gp
     dy_field(i,0,1) = 0.0_gp
     dy_field(i,2*n2_old,1) = 0.0_gp
     dy_field(i,2*n2_old+1,1) = 0.0_gp
     dx_field(i,0,1) = 0.0_gp
     dx_field(i,2*n2_old,1) = 0.0_gp
     dx_field(i,2*n2_old+1,1) = 0.0_gp
  end do

  do j=0,2*n2_old+1
     psifscfold(0,j,1) = 0.0_gp
     psifscfold(2*n1_old,j,1) = 0.0_gp
     psifscfold(2*n1_old+1,j,1) = 0.0_gp
     dy_field(0,j,1) = 0.0_gp
     dy_field(2*n1_old,j,1) = 0.0_gp
     dy_field(2*n1_old+1,j,1) = 0.0_gp
     dx_field(0,j,1) = 0.0_gp
     dx_field(2*n1_old,j,1) = 0.0_gp
     dx_field(2*n1_old+1,j,1) = 0.0_gp
  end do
  
!!$     allocate(x_phi(0:nd+ndebug),stat=i_stat )
!!$     call memocc(i_stat,x_phi,'x_phi',subname)

  x_phi=f_malloc(bounds=(/0.to.nd/),id='x_phi')

  allocate(y_phi(0:nd+ndebug) ,stat=i_stat )
  call memocc(i_stat,y_phi,'y_phi',subname)

  print *, " scaling function for interpolation "

!    call scaling_function(itype,nd,nrange,x_phi,y_phi) 
!cut the size of the array to exclude points outside support

!    if( abs(y_phi(nd/2)-1)>1.0e-10 ) then
  do i=0,nd
    write(17,*)i,x_phi(i),y_phi(i)
  end do
  !call scaling_function4b2B(itype,nd,nrange,x_phi,y_phi)
  call my_scaling_function4b2B(itype,nd,nrange,x_phi,y_phi)  
  do i=0,nd
    write(18,*)i,x_phi(i),y_phi(i)
  end do

!      stop " wrong scaling function 4b2B: not a centered one "
!    endif

  call f_free(x_phi)

!    i_all=-product(shape(x_phi))*kind(x_phi)
!    deallocate(x_phi,stat=i_stat)
!    call memocc(i_stat,i_all,'x_phi',subname)

!call interpolate_and_transpose(h,t0,nphi,nrange,phi,ndat,nin,psi_in,nout,psi_out)

!1d
!call my_interpolate_and_transpose(dx/hx,nd,nrange,y_phi,1,&
!         (2*n1_old+2+2*nb1),psifscfold,(2*n1+2+2*nb1),psifscf)

  print *,'interpolating first dimension...'
    call my_morph_and_transpose(hx,dx_field/hx,nd,nrange,y_phi,(2*n2_old+2+2*nb2),&
         (2*n1_old+2+2*nb1),psifscfold,(2*n1+2+2*nb1),psi_w)

    call my_morph_and_transpose(hx,dx_field/hx,nd,nrange,y_phi,(2*n2_old+2+2*nb2),&
         (2*n1_old+2+2*nb1),dy_field,(2*n1+2+2*nb1),psifscf)


  ! make an empty border
  do j=0,2*n2_old+1
     psi_w(j,0,1) = 0.0_gp
     psi_w(j,2*n1_old,1) = 0.0_gp
     psi_w(j,2*n1_old+1,1) = 0.0_gp
     psifscf(j,0,1) = 0.0_gp
     psifscf(j,2*n1_old,1) = 0.0_gp
     psifscf(j,2*n1_old+1,1) = 0.0_gp
  end do

  do i=0,2*n1_old+1
     psi_w(0,i,1) = 0.0_gp
     psi_w(2*n2_old,i,1) = 0.0_gp
     psi_w(2*n2_old+1,i,1) = 0.0_gp
     psifscf(0,i,1) = 0.0_gp
     psifscf(2*n2_old,i,1) = 0.0_gp
     psifscf(2*n2_old+1,i,1) = 0.0_gp
  end do

  x=-n1_old*hx!0.d0
  do i=0,2*n1+1
    y=-n2_old*hy!0.d0
    do j=0,2*n2+1
       z=0.d0
       do k=1,1!1,2*n3+1
          write(103,*) x,y,z,psi_w(j,i,k)
          z=z+hz
       end do
       y=y+hy
    end do
    write(103,*) ''
    x=x+hx
  end do

  x=-n1_old*hx!0.d0
  do i=0,2*n1_old+1
    y=-n2_old*hy!0.d0
    do j=0,2*n2_old+1
       z=0.d0
       do k=1,1!1,2*n3_old+1
          write(100,*) x,y,z,psifscfold(i,j,k),dx_field(i,j,k),psifscf(j,i,k),dy_field(j,i,k)
          z=z+hz
       end do
       !write(100,*) ''
       y=y+hy
    end do
    write(100,*) ''
    x=x+hx
  end do



!print*,'...interpolating second dimension...'
    call my_morph_and_transpose(hy,dy_field/hy,nd,nrange,y_phi,(2*n1+2+2*nb1),&
         (2*n2_old+2+2*nb2),psi_w,(2*n2+2+2*nb2),psi_w2)

  x=-n1_old*hx!0.d0
  do i=0,2*n1+1
    y=-n2_old*hy!0.d0
    do j=0,2*n2+1
       z=0.d0
       do k=1,1!1,2*n3+1
          write(101,*) x,y,z,psi_w2(i,j,k)
          z=z+hz
       end do
       y=y+hy
    end do
    x=x+hx
  end do

deallocate(dy_field)
deallocate(dx_field)


  ! original norm
  sum=0.d0
  sum2=0.d0
  sumx=0.d0
  sumx2=0.d0
  sumy=0.d0
  sumy2=0.d0
  z=0.d0
  k=1
  xmin=-n1_old*hx
  xmax=-n1_old*hx+(2*n1_old)*hx
  ymin=-n2_old*hy
  ymax=-n2_old*hx+(2*n2_old)*hy

!ignore last rows so that we have a properly symmetric grid
  sum = sum+psifscfold(0,0,k)+psifscfold(2*n1_old,0,k)+psifscfold(2*n1_old,2*n2_old,k)+psifscfold(0,2*n2_old,k)
  sum2 = sum2+psi_w2(0,0,k)+psi_w2(2*n1_old,0,k)+psi_w2(2*n1_old,2*n2_old,k)+psi_w2(0,2*n2_old,k)
  sumx = sumx+xmin*psifscfold(0,0,k)+xmax*psifscfold(2*n1_old,0,k)&
     +xmax*psifscfold(2*n1_old,2*n2_old,k)+xmin*psifscfold(0,2*n2_old,k)
  sumx2 = sumx2+xmin*psi_w2(0,0,k)+xmax*psi_w2(2*n1_old,0,k)+xmax*psi_w2(2*n1_old,2*n2_old,k)+xmin*psi_w2(0,2*n2_old,k)
  sumy = sumy+ymin*psifscfold(0,0,k)+ymin*psifscfold(2*n1_old,0,k)&
     +ymax*psifscfold(2*n1_old,2*n2_old,k)+ymax*psifscfold(0,2*n2_old,k)
  sumy2 = sumy2+ymin*psi_w2(0,0,k)+ymin*psi_w2(2*n1_old,0,k)+ymax*psi_w2(2*n1_old,2*n2_old,k)+ymax*psi_w2(0,2*n2_old,k)
print*,sumx,xmin,xmax,psifscfold(0,0,k),psifscfold(2*n1_old,0,k),psifscfold(2*n1_old,2*n2_old,k),psifscfold(0,2*n2_old,k)
  x=-n1_old*hx+hx
  do i=1,2*n1_old-1
    sum = sum+2.0d0*(psifscfold(i,0,k)+psifscfold(i,2*n2_old,k))
    sum2 = sum2+2.0d0*(psi_w2(i,0,k)+psi_w2(i,2*n2_old,k))
    sumx = sumx+2.0d0*x*(psifscfold(i,0,k)+psifscfold(i,2*n2_old,k))
    sumx2 = sumx2+2.0d0*x*(psi_w2(i,0,k)+psi_w2(i,2*n2_old,k))
    sumy = sumy+2.0d0*(ymin*psifscfold(i,0,k)+ymax*psifscfold(i,2*n2_old,k))
    sumy2 = sumy2+2.0d0*(ymin*psi_w2(i,0,k)+ymax*psi_w2(i,2*n2_old,k))
    x=x+hx
  end do

  y=-n2_old*hx+hx
  do j=1,2*n2_old-1
    sum = sum+2.0d0*(psifscfold(0,j,k)+psifscfold(2*n1_old,j,k))
    sum2 = sum2+2.0d0*(psi_w2(0,j,k)+psi_w2(2*n1_old,j,k))
    sumx = sumx+2.0d0*(xmin*psifscfold(0,j,k)+xmax*psifscfold(2*n1_old,j,k))
    sumx2 = sumx2+2.0d0*(xmin*psi_w2(0,j,k)+xmax*psi_w2(2*n1_old,j,k))
    sumy = sumy+2.0d0*y*(psifscfold(0,j,k)+psifscfold(2*n1_old,j,k))
    sumy2 = sumy2+2.0d0*y*(psi_w2(0,j,k)+psi_w2(2*n1_old,j,k))
    y=y+hy
  end do

  x=-n1_old*hx+hx
  do i=1,2*n1_old-1
    y=-n2_old*hy+hy
    do j=1,2*n2_old-1
       sum = sum + 4.0d0*psifscfold(i,j,k)
       sum2 = sum2 + 4.0d0*psi_w2(i,j,k)
       sumx = sumx + 4.0d0*x*psifscfold(i,j,k)
       sumx2 = sumx2 + 4.0d0*x*psi_w2(i,j,k)
       sumy = sumy + 4.0d0*y*psifscfold(i,j,k)
       sumy2 = sumy2 + 4.0d0*y*psi_w2(i,j,k)
       y=y+hy
    end do
    x=x+hx
  end do

  sum=sum*0.25d0*hx*hy
  sum2=sum2*0.25d0*hx*hy
  sumx=sumx*0.25d0*hx*hy
  sumx2=sumx2*0.25d0*hx*hy
  sumy=sumy*0.25d0*hx*hy
  sumy2=sumy2*0.25d0*hx*hy

  print*,'norm in',sum,'norm out',sum2
  print*,'x in',sumx,'x out',sumx2,'y in',sumy,'y out',sumy2,'x/y in',sumx/sumy,'x/y out',sumx2/sumy2


stop


print*,'interpolating first dimension...'
    call my_interpolate_and_transpose(hx,dx/hx,nd,nrange,y_phi,(2*n3_old+2+2*nb3)*(2*n2_old+2+2*nb2),&
         (2*n1_old+2+2*nb1),psifscfold,(2*n1+2+2*nb1),psi_w)

  print *,'...interpolating second dimension...'
  call my_interpolate_and_transpose(dy/hy,nd,nrange,y_phi,(2*n3_old+2+2*nb3)*(2*n1+2+2*nb1),&
         (2*n2_old+2+2*nb2),psi_w,(2*n2+2+2*nb2),psi_w2) 

  print *,'...interpolating third dimension...'
  call my_interpolate_and_transpose(dz/hz,nd,nrange,y_phi,(2*n2+2+2*nb2)*(2*n1+2+2*nb1),&
    (2*n3_old+2+2*nb3),psi_w2,(2*n3+2+2*nb3),psifscf) 

  print *,'done'

  !print psifscf to check
  !x=0.d0
  !do i=0,2*n1+1
  !  write(101,*) x,psifscf(i,1,1)
  !  x=x+hx
  !end do

  x=0.d0
  do i=0,2*n1+1
    y=0.d0
    do j=0,2*n2+1
       z=0.d0
       do k=1,2*n3+1
          write(101,*) x,y,z,psifscf(i,j,k)
          z=z+hz
       end do
       y=y+hy
    end do
    x=x+hx
  end do

  i_all=-product(shape(psifscf))*kind(psifscf)
  deallocate(psifscf,stat=i_stat)
  call memocc(i_stat,i_all,'psifscf',subname)



  i_all=-product(shape(psifscfold))*kind(psifscfold)
  deallocate(psifscfold,stat=i_stat)
  call memocc(i_stat,i_all,'psifscfold',subname)

  i_all=-product(shape(psi_w2))*kind(psi_w2)
  deallocate(psi_w2,stat=i_stat)
  call memocc(i_stat,i_all,'psi_w2',subname)

  i_all=-product(shape(psi_w))*kind(psi_w)
  deallocate(psi_w,stat=i_stat)
  call memocc(i_stat,i_all,'psi_w',subname)

  i_all=-product(shape(y_phi))*kind(y_phi)
  deallocate(y_phi,stat=i_stat)
  call memocc(i_stat,i_all,'y_phi',subname)

  stop
  !end testing interpolate


  !Allocations
  work      =f_malloc(bounds=(/-n_range .to. n_range /),id='work')
  print *,'test'
  kernel_scf=f_malloc(bounds=(/-n_range .to. n_range /),id='kernel_scf')

  !allocate(work(-n_range:n_range), stat=i_stat)
  !allocate(kernel_scf(-n_range:n_range), stat=i_stat)


  !Convergence test for the subroutine gauss_conv_scf 
  !with 6 different values of n_points and 7 different pgauss
  multiple_naive=f_malloc(lbounds=(/1,1,-n_range/),ubounds=(/6,7,n_range/),id='multiple_naive',routine_id=subname)
  timings       =f_malloc((/2,7,7/),id='timings')
!  allocate(multiple_naive(1:6,1:7,-n_range:n_range), stat = i_stat) 
!  allocate(timings(1:2,1:7,1:7), stat = i_stat)
  timings = 0.0_dp

  n_points_list = (/ 2, 8, 32, 64, 256, 512, 0 /) 



  do i=1,6

     !n_points_list(i) = 8*2**i
     n_points = n_points_list(i)

     !Number of integration points = 2*itype_scf*n_points
     n_scf=2*itype_scf*n_points

     !Other allocations
     x_scf=f_malloc(bounds=(/0.to.n_scf/),id='x_scf')
     y_scf=f_malloc(bounds=(/0.to.n_scf/),id='y_scf')
     gaussian=f_malloc(bounds=(/1 .to. 7, 0 .to. n_scf/),id='gaussian')
     !allocate(x_scf(0:n_scf),stat=i_stat)
     !allocate(y_scf(0:n_scf),stat=i_stat)
     !allocate(gaussian(1:7,0:n_scf),stat=i_stat)

     !Build the scaling function
     call scaling_function(itype_scf,n_scf,n_range,x_scf,y_scf)

     !Step grid for the integration
     dx = real(n_range,dp)/real(n_scf,dp)

     do j = 7,7!1,7
        !pgauss(j) = real(2.0,dp)**(2*jx-3)*(8.0/real(itype_scf,dp))**2
        !pgauss(j) = 1.0e-8_dp*2.0_dp**(6*(j-1))
        pgauss(j) = 8.0e-3_dp*1.25_dp**(6*(j-1))
        !compare the gaussians together
!        n_range=200
        if (i==1) call discretize_gaussian(n_range,1.d0,pgauss(j),0.1_gp,&
             hgrid,'gaussian'//trim(adjustl(yaml_toa(j))))
     end do

!stop
!!$     if (i == 6) then
!!$        do j = 1,7
!!$           !Test: is the width of the Gaussian ok?
!!$           kern = 0.0_dp
!!$           do k=0,n_scf
!!$              gaussian(j,k) = dexp(-pgauss(j)*(hgrid*x_scf(k))**2) ! <<< For subsequent plotting purposes...
!!$              absci = x_scf(k)
!!$              absci = absci*absci*hgrid**2
!!$              kern = kern + dexp(-pgauss(j)*absci)
!!$           end do
!!$           kern = kern*dx
!!$           write(*,'(a28,1pe9.3,a7,1pe23.16)') 'Gaussian test ---> pgauss = ', pgauss(j), ' err = ', &
!!$                abs(kern-sqrt(4*datan(1.d0)/pgauss(j)))
!!$        end do
!!$        ! plot to file the Gaussian & the scaling function 
!!$        ! >>> only in the case of the highest n_points <<<
!!$        write(chain,'(a18)') 'gauss_scf_plot.out'
!!$        !write(chain,j)
!!$        open(unit=51, file = chain, status = 'unknown')
!!$        write(51, *) '# The columns correspond to x_scf, y_scf, gaussian(pgauss),'
!!$        write(51, '(a18,7(1pe11.3))') '# where pgauss =', pgauss(:)
!!$        write(51, *) ''
!!$        do k = 0, n_scf
!!$           write(51, '(9(1pe21.12e3))') x_scf(k), y_scf(k), &
!!$                gaussian(1,k), gaussian(2,k), gaussian(3,k), &
!!$                gaussian(4,k), gaussian(5,k), gaussian(6,k), gaussian(7,k)
!!$        end do
!!$        close(51) 
!!$     end if
!!$
!!$
!!$     ! !Test: is the norm of the ISF correct?
!!$     ! kern = 0.0_dp
!!$     ! do k=0,n_scf
!!$     !    kern = kern + y_scf(k)
!!$     ! end do
!!$     ! kern = kern*dx
!!$     ! write(*,*) 'ISF test      --->', pgauss, kern
!!$
!!$     do j = 1,7
!!$        !template: gauss_conv_scf(itype_scf,pgauss,hgrid,dx,n_range,n_scf,x_scf,y_scf,kernel_scf,work)
!!$        call cpu_time(t0)
!!$        if (i == 1 .and. timings_switch) then
!!$           do k = 1,1000
!!$              call my_gauss_conv_scf(itype_scf, pgauss(j),0.0_dp, hgrid, dx, n_range, n_scf, x_scf, y_scf, kernel_scf, work)
!!$           end do
!!$        else
!!$           call my_gauss_conv_scf(itype_scf, pgauss(j),0.0_dp, hgrid, dx, n_range, n_scf, x_scf, y_scf, kernel_scf, work)
!!$        end if
!!$        call cpu_time(t1)
!!$        timings(1,i,j) = t1-t0
!!$        !stores the result in multiple_naive
!!$        multiple_naive(i,j,:) = kernel_scf(:)
!!$     end do

     
   call f_free(x_scf)
   call f_free(y_scf)
   call f_free(gaussian)
   !deallocate(x_scf, stat = i_stat)
   !deallocate(y_scf, stat = i_stat)

  end do

!!$  do j = 1,7
!!$     write(chain,'(a28,i1,a4)') 'gauss_conv_scf_integral_plot', j, '.out'
!!$     !write(chain,j)
!!$     open(unit=51, file = chain, status = 'unknown')
!!$     write(51, *) '# pgauss =', pgauss(j)
!!$     write(51, *) '# 1st column    >>> k in [-n_range:n_range]'
!!$     write(51, *) '# other columns >>> ', &
!!$          'naive integrals with no. of integration points = 2*itype_scf*N, where N =', n_points_list(1:6)
!!$     write(51, *) ' '
!!$     do k = -n_range, n_range
!!$        write(51, '(i4, 6(1pe21.12e3))') k, &
!!$             multiple_naive(1,j,k), & 
!!$             multiple_naive(2,j,k), &
!!$             multiple_naive(3,j,k), &
!!$             multiple_naive(4,j,k), &
!!$             multiple_naive(5,j,k), &
!!$             multiple_naive(6,j,k)
!!$     end do
!!$     close(51)
!!$  end do
!!$
!!$
!!$  ! Evaluate the absolute error with respect to the integral computed with the highest n_points
!!$  do i = 1,5
!!$     do j = 1,7
!!$        do k = -n_range, n_range
!!$           multiple_naive(i,j,k) = abs(multiple_naive(i,j,k)-multiple_naive(6,j,k))
!!$        end do
!!$     end do
!!$  end do
!!$
!!$
!!$
!!$  ! write(chain,'(a26)') 'naive_integral_summary.out'
!!$  ! open(unit=52, file = chain, status = 'replace', position = 'append')
!!$  ! write(52, *) '# 1st column    >>> log10(pgauss)'
!!$  ! write(52, *) '# other columns >>> ', &
!!$  !      'maxval(abs("integral with 2*itype_scf*512 points" - "integral with 2*itype_scf*N points")), where N =', n_points_list(1:5)
!!$  ! write(52, *) ' '
!!$  ! do j = 1,7
!!$  !    write(52, '(1pe10.3,5(1pe20.12e3))') log10(pgauss(j)), &
!!$  !         maxval(multiple_naive(1,j,:)), &
!!$  !         maxval(multiple_naive(2,j,:)), &
!!$  !         maxval(multiple_naive(3,j,:)), &
!!$  !         maxval(multiple_naive(4,j,:)), &
!!$  !         maxval(multiple_naive(5,j,:))
!!$  !    !maxval(multiple_naive(6,j,:))
!!$  ! end do
!!$  ! close(52)
!!$
!!$
!!$  write(chain,'(a26)') 'naive_integral_summary.out'
!!$  open(unit=52, file = chain, status = 'replace', position = 'append')
!!$  write(52, *) '# 1st column    >>> N (N.B.: the no. of integration points is = 2*itype_scf*N)'
!!$  write(52, *) '# other columns >>> ', &
!!$       'maxval(abs("integral with 2*itype_scf*512 points" - "integral with 2*itype_scf*N points"))'
!!$  write(52, *) '# for pgauss =', pgauss(1:7)
!!$  write(52, *) ' '
!!$  do i = 1,5
!!$     !(i3,7(1pe20.12e3)
!!$     write(52, '(i3, 7(1pe21.12e3))') n_points_list(i), &
!!$          maxval(multiple_naive(i,1,:)), &
!!$          maxval(multiple_naive(i,2,:)), &
!!$          maxval(multiple_naive(i,3,:)), &
!!$          maxval(multiple_naive(i,4,:)), &
!!$          maxval(multiple_naive(i,5,:)), &
!!$          maxval(multiple_naive(i,6,:)), &
!!$          maxval(multiple_naive(i,7,:))
!!$     !maxval(multiple_naive(6,j,:))
!!$  end do
!!$  close(52)


  ! Convergence test for my_analytic_integral routine
  analytic_integral_result=f_malloc((/7,7,n_range+1/),lbounds=(/1,1,0/),id='analytic_integral_result')
  analytic_vs_naive=f_malloc(bounds=(/0.to.n_range,0.to.n_range/),id='analytic_vs_naive')
  !allocate(analytic_integral_result(1:7,1:7,0:n_range), stat=i_stat)
  !allocate(analytic_vs_naive(0:n_range,0:n_range), stat=i_stat)

!!$
!!$  do i = 1,7
!!$     n_points_list(i) = 2**(i+4)
!!$  end do
!!$
!!$
!!$  open(unit=54, file = 'analytic_vs_naive_summary.out', status = 'replace', position = 'append')
!!$  write(54, *) '# 1st column    >>> N, number of Fourier terms taken into account'
!!$  write(54, *) '# other columns >>> ', &
!!$       'maxval(abs("naive integral with 2*itype_scf*512 integration points" - ', &
!!$       '"analytic_integral with N Fourier terms")) for pgauss =', pgauss(:)
!!$  write(54, *) ' '

!!$
!!$  do i = 1,7 ! loop over integrals computed with different numbers of Fourier components
!!$     do j = 1,7 ! loop over different Gaussians
!!$        call cpu_time(t0)
!!$        ! template: my_analytic_integral(alpha,ntot,m,fwork,nf)
!!$        if (i == 4 .and. timings_switch) then
!!$           do k = 1, 1000
!!$              call my_analytic_integral(hgrid*sqrt(pgauss(j)),0.0_dp/hgrid,&
!!$                   n_range, itype_scf, analytic_integral_result(i,j,-n_range:n_range), fISF, n_points_list(i))
!!$           end do
!!$        else 
!!$           call my_analytic_integral(hgrid*sqrt(pgauss(j)),0.0_dp/hgrid,&
!!$                n_range, itype_scf, analytic_integral_result(i,j,-n_range:n_range), fISF, n_points_list(i))
!!$        end if
!!$        call cpu_time(t1)
!!$        timings(2,i,j) = t1-t0
!!$        ! norm infinity with respect to the naive integral computed with the highest n_points
!!$        analytic_vs_naive(i,j) = maxval(abs(analytic_integral_result(i,j,0:n_range) - & 
!!$             multiple_naive(6,j,0:n_range)))
!!$     end do
!!$
!!$     write(54,'(i5, 7(1pe21.12e3))') n_points_list(i), analytic_vs_naive(i,1),&
!!$          analytic_vs_naive(i,2),&
!!$          analytic_vs_naive(i,3),&
!!$          analytic_vs_naive(i,4),&
!!$          analytic_vs_naive(i,5),&
!!$          analytic_vs_naive(i,6),&
!!$          analytic_vs_naive(i,7)
!!$  end do
!!$
!!$  close(54)

!!$
!!$  do j = 1,7
!!$     write(chain,'(a25,i1,a4)') 'my_analytic_integral_plot', j, '.out'
!!$     open(unit=61, file = chain, status = 'unknown')
!!$     write(61, *) '# pgauss =', pgauss(j)
!!$     write(61, *) '# 1st column    >>> k in [0:n_range]'
!!$     write(61, *) '# other columns >>> ', &
!!$          'analytic_integral with N Fourier terms, where N =', n_points_list(:)
!!$     write(61, *) '# The last column contains the result from the naive method, to be used as benchmark.'
!!$     write(61, *) ' '
!!$     do k = 0, n_range
!!$        write(61, '(i3, 8(1pe21.12e3))') k, &
!!$             analytic_integral_result(1,j,k), & 
!!$             analytic_integral_result(2,j,k), &
!!$             analytic_integral_result(3,j,k), &
!!$             analytic_integral_result(4,j,k), &
!!$             analytic_integral_result(5,j,k), &
!!$             analytic_integral_result(6,j,k), &
!!$             analytic_integral_result(7,j,k), & 
!!$             multiple_naive(6,j,k) ! <<< reference
!!$     end do
!!$     close(61)
!!$  end do
!!$
!!$  !print moments of the gaussian
!!$  do j = 1,7
!!$     moment=0.0_dp
!!$     do k=0,n_scf
!!$        moment=moment+gaussian(j,k)
!!$     end do
!!$     moment = moment*dx
!!$     print '(a,4(1pe25.17))','moment',&
!!$          pgauss(j),moment,sqrt(4*datan(1.d0)/pgauss(j)),moment-sqrt(4*datan(1.d0)/pgauss(j))
!!$  end do
!!$
!!$  !print moments of the regularized gaussian
!!$  !also print collocated values
!!$  do j = 1,7
!!$     moment=0.0_dp
!!$     kern=0.0_dp
!!$     do k=-n_range,n_range
!!$        write(20+j,'(i4,4(1pe26.17e3))')&
!!$             k,dexp(-pgauss(j)*(hgrid*real(k,dp))**2),multiple_naive(6,j,k)
!!$        kern=kern+dexp(-pgauss(j)*(hgrid*real(k,dp))**2)
!!$        moment=moment+multiple_naive(6,j,k)
!!$     end do
!!$     moment = moment*hgrid
!!$     kern = kern*hgrid
!!$     print '(a,8(1pe25.17))','momentphi',&
!!$          pgauss(j),moment,kern,sqrt(4*datan(1.d0)/pgauss(j)),&
!!$          moment-sqrt(4*datan(1.d0)/pgauss(j)),kern-sqrt(4*datan(1.d0)/pgauss(j))
!!$  end do
!!$
!!$
!!$  !Timings
!!$  if (timings_switch) then
!!$     open(unit=62, file = 'timings.out', status = 'unknown')
!!$     ! write(62,*) '# 1st column >>> index running over different n_points or nf'
!!$     ! write(62,*) '# 2nd column >>> index running over Gaussian widths'
!!$     ! write(62,*) '# 3rd column >>> naive integral'
!!$     ! write(62,*) '# 4th column >>> analytic integral'
!!$     write(62,*) '# >>> Elapsed times for a bunch of 1,000 integrations <<<'
!!$     write(62,*) '# 1st column >>> index running over Gaussian widths'
!!$     write(62,*) '# 2nd column >>> naive integral with n_points = 16'
!!$     write(62,*) '# 3rd column >>> analytic integral with nf = 256'
!!$
!!$     do j = 1,7
!!$        !do i = 1,7
!!$        !   write(62, '(2(i2), 2(1pe15.7))') i, j, timings(1,i,j), timings(2,i,j)
!!$        !end do
!!$        ! timings(1,1,:) >>> n_points = 16 (already enough for achieving an acceptable accuracy)
!!$        ! timings(2,4,:) >>> nf = 256 (enough for achieving an acceptable accuracy)
!!$        write(62, '(1(i2), 2(1pe15.7))') j, timings(1,1,j), timings(2,4,j)
!!$     end do
!!$     close(62)
!!$  end if


  write(*,*) ' '
  write(*,*) 'Please give a look at the output files!'
  write(*,*) "Here's the list:"
  write(*,*) 'naive_integral_summary.out'
  write(*,*) 'analytic_vs_naive_summary.out'
  if (timings_switch) write(*,*) 'timings.out'

  write(*,*) ' '
  write(*,*) 'Moreover:'
  write(*,*) 'gauss_scf_plot.out'
  write(*,*) 'gauss_conv_scf_integral_plot(1-7).out'
  write(*,*) 'my_analytic_integral_plot(1-7).out'



  deallocate(gaussian, stat = i_stat)
  deallocate(work, stat = i_stat)
  deallocate(kernel_scf, stat = i_stat)
  deallocate(x_scf, stat = i_stat)
  deallocate(y_scf, stat = i_stat)
  deallocate(analytic_integral_result, stat=i_stat)
  deallocate(analytic_vs_naive, stat=i_stat)
  deallocate(timings, stat=i_stat)

  call f_malloc_finalize()

end program PS_Integral


function lr_gauss(x,shift,shrink)
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  real(dp) :: lr_gauss, x, shift, shrink

  lr_gauss = exp(-((x - shift) * shrink)**2)

end function lr_gauss


!> takes a Gaussian of exponent pgauss an center x0 and discretize it on a grid of size 1 in units of sqrt(0.5*[pgauss])
!! f(x)=fac*exp(-pgauss*(x-x0)**2)
subroutine discretize_gaussian(nrange,fac,pgauss,x0,hgrid,filename)
  use module_base
  use gaussians
  !not yet use dynamic_memory !let's start
  implicit none
  integer, intent(in) :: nrange
  real(gp), intent(in) :: fac,pgauss,x0,hgrid
  character(len=*), intent(in) :: filename
  !local variables
  integer, parameter :: itype_scf=16,npoints=2**6,nmoms=6
  integer :: k,n_scf,n_range,j
  real(gp) :: dx,x
  real(gp), dimension(0:nmoms,3) :: moments !<to check all methods
  real(dp), dimension(0:2048) :: fISF
  real(gp), dimension(:), allocatable :: f_i,f_phi_i,f_phi_i2,x_scf,y_scf

  !allocate the three gaussians arrays
  allocate(f_i(-nrange:nrange))
  allocate(f_phi_i(-nrange:nrange))
  allocate(f_phi_i2(-nrange:nrange))

  !Number of integration points = 2*itype_scf*n_points
  n_scf=2*itype_scf*npoints

  !Other allocations
  allocate(x_scf(0:n_scf))
  allocate(y_scf(0:n_scf))

  !Build the scaling function
  call scaling_function(itype_scf,n_scf,n_range,x_scf,y_scf)

  !Step grid for the integration
  dx = real(nrange,gp)/real(n_scf,gp)

  !first, collocate the gaussian
  do k=-nrange,nrange
     f_i(k)=fac*exp(-pgauss*(hgrid*real(k,gp)-x0)**2)
  end do

  f_phi_i2=0.0_gp
  f_phi_i=0.0_gp
  !then, calculate its values integrated with SCF
  !use the i2 array as a workspace
  call my_gauss_conv_scf(itype_scf,pgauss,x0,hgrid,dx,nrange,n_scf,x_scf,y_scf,f_phi_i,f_phi_i2)
  call dscal(2*nrange+1,fac,f_phi_i,1)
  !Only itype=8,14,16,20,24,30,40,50,60,100
  call four_isf_coeffs(itype_scf,fISF)
  f_phi_i2=0.0_gp
  !then calculate the values via the analytic method (for the moment only on the second half)
  !this routine is not useful for performances of precision and speed
  !call my_analytic_integral(hgrid*sqrt(pgauss),x0/hgrid,&
  !     nrange,itype_scf,f_phi_i2,fISF,64)
!  do k=-nrange,-1
!     f_phi_i2(k)=f_phi_i2(-k)
!  end do
  !call dscal(2*nrange+1,fac,f_phi_i2,1)
  
  !use daubechies wavelets for expressing the function
!  call gauss_to_daub(hgrid,fac,x0,sqrt(0.5/pgauss),0,&!no err, errsuc
!       2*nrange+1,n_left,n_right,c,err_norm,&         !no err_wav. nmax instead of n_intvx
!       ww,nwork,.false.)                             !added work arrays ww with dimension nwork
  !create the isf array in the gaussian module
  call initialize_real_space_conversion()
  do j=-nrange,nrange
     f_phi_i2(j)=scfdotf(j,hgrid,pgauss,x0,0)
     !print *,'j',j,f_phi_i2(j)
  end do
  !use the elemental property of the scfdotf function
  f_phi_i2=scfdotf((/(j,j=-nrange,nrange)/),hgrid,pgauss,x0,0)
  call finalize_real_space_conversion('discretize_gaussian')

  !then print the results
  moments=0.0_dp
  !open(unit=200, file = trim(filename), status = 'unknown')
  do j=0,nmoms
     open(unit=200+j)
     do k=-nrange,nrange
        x=real(k,gp)*hgrid-x0
        moments(j,1)=moments(j,1)+x**j*f_i(k)
        moments(j,2)=moments(j,2)+x**j*f_phi_i(k)
        moments(j,3)=moments(j,3)+x**j*f_phi_i2(k)
        write(200+j,'(i4,3(1pe26.17e3))')k,x**j*f_i(k),x**j*f_phi_i(k),x**j*f_phi_i2(k)
     end do
     moments(j,:) = moments(j,:)*hgrid
     print '(a,i3,8(1pe14.5))',trim(filename),j,&
          pgauss,moments(j,:),gauint0(pgauss,j),moments(j,:)-gauint0(pgauss,j)
     close(unit=200+j)
  end do

  deallocate(f_i,f_phi_i,f_phi_i2)

end subroutine discretize_gaussian


!> Write a routine which performs the integration of a function on a given grid
subroutine my_gauss_conv_scf(itype_scf,pgauss,x0,hgrid,dx,n_range,n_scf,x_scf,y_scf,kernel_scf,work)
  use module_base
  implicit none
  integer, intent(in) :: n_range,n_scf,itype_scf
  real(dp), intent(in) :: pgauss,hgrid,dx,x0
  real(dp), dimension(0:n_scf), intent(in) :: x_scf
  real(dp), dimension(0:n_scf), intent(in) :: y_scf
  real(dp), dimension(-n_range:n_range), intent(inout) :: work
  real(dp), dimension(-n_range:n_range), intent(inout) :: kernel_scf
  !local variables
  real(dp), parameter :: p0_ref = 1.0_dp
  integer :: n_iter,i_kern,i
  real(dp) :: p0_cell,p0gauss,absci,kern

  !Step grid for the integration
  !dx = real(n_range,dp)/real(n_scf,dp)

  !To have a correct integration
  p0_cell = p0_ref/(hgrid*hgrid)

  !write(*,*) 'p0_cell = ', p0_cell 

  !We calculate the number of iterations to go from pgauss to p0_ref
  n_iter = nint((log(pgauss) - log(p0_cell))/log(4.0_dp))

  write(*,*) 'n_iter = ', n_iter

  if (n_iter <= 0 .or. x0 /= 0.d0)then
     n_iter = 0
     p0gauss = pgauss
  else
     p0gauss = pgauss/4._dp**n_iter
  end if

!  write(*,*) 'p0gauss =',  p0gauss

  !Stupid integration
  !Do the integration with the exponential centered in i_kern
  kernel_scf(:) = 0.0_dp
  do i_kern=-n_range,n_range
     kern = 0.0_dp
     do i=0,n_scf
        absci = x_scf(i) - real(i_kern,dp)
        absci = x0-absci*hgrid !sign has changed here
        absci = absci*absci
        kern = kern + y_scf(i)*dexp(-p0gauss*absci)
     end do
     kernel_scf(-i_kern) = kern*dx
     !write(17,*)i_kern,kernel_scf(i_kern)
     !kernel_scf(-i_kern) = kern*dx
     !if (abs(kern) < 1.d-18) then
     !   !Too small not useful to calculate
     !   exit
     !end if
  end do

!  do i=0,n_scf
!     write(18,*)i,x_scf(i),y_scf(i)
!  end do

  !Start the iteration to go from p0gauss to pgauss
  call scf_recursion(itype_scf,n_iter,n_range,kernel_scf,work)

END SUBROUTINE my_gauss_conv_scf


!> Here alpha corresponds to sqrt(alpha) in mathematica
!! the final result is fwork(j+m)-fwork(j-m)
subroutine my_analytic_integral(alpha,x0,ntot,m,fwork,fISF,argument_nf)
  use module_base
  implicit none
  integer, intent(in) :: ntot,m
  real(dp), intent(in) :: alpha,x0 !<x0 is the deviation wrt the grid spacing
  real(dp), dimension(-ntot:ntot), intent(inout) :: fwork
  real(dp), dimension(0:2048), intent(in) :: fISF
  integer, intent(in) :: argument_nf 

  !local variables
  integer :: nf
  real(dp), parameter :: pi=3.1415926535897932384_dp
  logical :: flag,flag1,flag2
  integer :: j,q,jz
  real(dp) :: if,r1,r2,res,ypm,ymm,erfcpm,erfcmm,factor,re,ro,factorend

  !write(*,*) fISF(1000), fISF16(1000)

  !  if(present(argument_nf)) then 
  nf=argument_nf
  !  else
  !     nf = 64 ! "default value"
  !  endif

  flag=.false.
  factor=pi/real(2*m,dp)/alpha
  factorend=sqrt(pi)/alpha/real(4*m,dp)*exp((alpha*x0)**2)
  !fill work array
  !the calculation for j=0 can be separated from the rest
  !since it only contains error functions
  !the values of the loop are not symmetric in principle
  loop_nonzero: do j=-ntot,ntot
     ypm=alpha*(real(j+m,dp)+x0)
     ymm=alpha*(real(j-m,dp)+x0)
     call derfcf(erfcpm,ypm) ! erfcpm = erfc(ypm)
     call derfcf(erfcmm,ymm)

     !assume nf even
     res=0._dp
     !reso=0._dp
     do q=nf,2,-2
        !the sign of q only influences the imaginary part
        !so we can multiply by a factor of two

!!$        call wofz_mod(alpha,m,q,-j-m,r1,if,flag1)
!!$        call wofz_mod(alpha,m,q,-j+m,r2,if,flag2)
        call GcplxInt(alpha,m,q,-j-m,-x0,r1,if,flag1)
!!$        if (q==12 .and. j==15) then
!!$           print *,'here',alpha,m,q,j,x0,r1,if,flag1
!!$           stop
!!$        end if
        call GcplxInt(alpha,m,q,-j+m,-x0,r2,if,flag2)
        re=r1-r2
        flag=flag1 .or. flag2 .or. flag
        !if (flag) then
        !   print *,'here',r1,if,q,j,x0
        !   stop 
        !end if
!!$        call wofz_mod(alpha,m,q-1,-j-m,r1,if,flag1)
!!$        call wofz_mod(alpha,m,q-1,-j+m,r2,if,flag2)
        call GcplxInt(alpha,m,q-1,-j-m,-x0,r1,if,flag1)
        call GcplxInt(alpha,m,q-1,-j+m,-x0,r2,if,flag2)
        ro=r1-r2
        flag=flag1 .or. flag2 .or. flag
        !if (flag) then
        !   print *,'there',xo,y
        !   stop 
        !end if
        !write(16,'(2(i4),6(1pe15.7))')j,q,re,ro,erfcmm-erfcpm
        re=re*fISF(q)
        ro=ro*fISF(q-1)
        res=res+re-ro
     end do
     !q=0 
     !fwork(j)=derf(y)+rese-reso
     fwork(j)=erfcmm-erfcpm+2.0_dp*res!e-reso
     fwork(j)=factorend*fwork(j)
     !exit form the loop if it is close to zero
!     if (abs(fwork(j)) < 1.e-25_dp) exit loop_nonzero
     !write(17,'(i4,8(1pe15.7))')j,derf(y),erfcsgn,rese,reso,derf(y)+rese-reso,&
     !     -erfcsgn+rese-reso,erfcsgn+rese-reso
  end do loop_nonzero

  !check flag
  if (flag) then
     write (*,*)'value of alpha',alpha
     stop 'problem occurred in wofz'
  end if

  !put to zero the rest
  do jz=j+1,ntot
     fwork(jz)=0.0_dp
  end do

END SUBROUTINE my_analytic_integral



! >>> GARBAGE <<<


! open(unit=54, file = 'analytic_vs_naive_summary.out', status = 'replace', position = 'append')
! write(54, *) '# 1st column    >>> log10(pgauss)'
! write(54, *) '# other columns >>> ', &
!      'maxval(abs("naive integral with 512 points" - "analytic_integral with N Fourier terms")), where N =', n_points_list(1:7)
! write(54, *) '# Note that n_scf = 2*2*itype_scf*N'
! write(54, *) ' '


! do j = 1,7
!    ! Convergence test for my_analytic_integral
!    do i = 1,7
!       call my_analytic_integral(hgrid*sqrt(pgauss(j)),&
!            n_range, itype_scf, analytic_integral_result(i, 0:n_range), n_points_list(i))
!       !write(*,*) 'analytic_integral with nf =', 2**i, ' ---> done!'
!       analytic_vs_naive(i,j) = maxval(abs(analytic_integral_result(i,0:n_range) - & 
!            multiple_naive(6,j,0:n_range))) ! NB: norm infinity
!       !write(*,*)  analytic_integral_result(0), multiple_naive(6,j,0)
!       !analytic_vs_naive(i-5,j) = analytic_integral_result(0)
!    end do

!    write(54,'(1pe8.2e02,7(1pe20.12e3))') log10(pgauss(j)),&
!         analytic_vs_naive(1,j),&
!         analytic_vs_naive(2,j),&
!         analytic_vs_naive(3,j),&
!         analytic_vs_naive(4,j),&
!         analytic_vs_naive(5,j),&
!         analytic_vs_naive(6,j),&
!         analytic_vs_naive(7,j)

!    ! write(54,'(i1,7(1pe20.12e3))') i,&
!    !      analytic_vs_naive(i,1),&
!    !      analytic_vs_naive(i,2),&
!    !      analytic_vs_naive(i,3),&
!    !      analytic_vs_naive(i,4),&
!    !      analytic_vs_naive(i,5),&
!    !      analytic_vs_naive(i,6),&
!    !      analytic_vs_naive(i,7)



! end do

! close(54)


!> Needed to calculate the integral of a ISF with a Gaussian
!! this is exp(-y**2) W(x + i y), with x=(pi q)/(2m alpha) and y= (jm+t) alpha
subroutine GcplxInt(alpha,m,q,jm,t,u,v,flag)
  use module_base
  implicit none
  real(dp), intent(in) :: alpha,t
  integer, intent(in) :: q,m,jm
  logical, intent(out) :: flag
  real(dp), intent(out) :: u,v
  !local variables
  real(dp), parameter :: factor=1.12837916709551257388_dp,rmaxreal=0.5e+154_dp
  real(dp), parameter :: rmaxexp=708.503061461606_dp !n(c) ,rmaxgoni=3.53711887601422e+15_dp
  real(dp), parameter :: pi=3.1415926535897932384_dp
  logical :: a,b
  integer :: j,n,i,kapn,nu,np1,multiple
  real(dp) :: xabs,yabs,x,y,qrho,xquad,yquad,xsum,ysum,xaux,u1,v1,u2,v2,daux,h,qlambda,h2
  real(dp) :: rx,ry,sx,sy,tx,ty,c,w1,xabsq,xi,yi,fac,yquadmod!,tt

  flag = .false.

  fac=pi/real(2*m,dp)/alpha

  xi=fac*real(q,dp)
  yi=alpha*(real(jm,dp)+t)
  xabs = dabs(xi)
  yabs = dabs(yi)
  x    = xabs/6.3_dp
  y    = yabs/4.4_dp

!     the following if-statement protects
!     qrho = (x**2 + y**2) against overflow
  if ((xabs > rmaxreal).or.(yabs > rmaxreal)) then
     flag = .true.
     return
  end if

  qrho = x**2 + y**2

  xabsq = xabs**2
  xquad = xabsq - yabs**2
  yquad = 2*xabs*yabs !this is equal to abs((pi q (jm+t))/m)
  !must subtract a arbitrary multiple of the period (2m), chosen from q*jm
  multiple=abs(q*jm)/(2*m)
  !tt=(real(abs(q*jm)-multiple*2*m,dp)+t)/real(m,dp)*pi-real(modulo(abs(q*jm),2*m),dp)/real(m,dp)*pi
  !if (abs(tt) > 1.e-16_gp) print *,'test',q*jm,multiple,abs(q*jm)-multiple*2*m,&
  !   (real(abs(q*jm)-multiple*2*m,dp)+t)/real(m,dp)*pi,&
  !   real(modulo(abs(q*jm),2*m),dp)/real(m,dp)*pi,tt
  multiple=abs(q*jm)-multiple*2*m
  !yquad to be passed to trigonometric functions
  yquadmod = (real(multiple,dp)+real(abs(q),dp)*t)/real(m,dp)*pi
  !yquadmod = real(modulo(abs(q*jm),2*m),dp)/real(m,dp)*pi

  a = qrho < 0.085264_dp

  if (a) then
     !if (qrho < 0.085264d0) then the faddeeva-function is evaluated
     !using a power-series (abramowitz/stegun, equation (7.1.5), p.297)
     !n is the minimum number of terms needed to obtain the required
     !accuracy
     
     qrho  = (1._dp-0.85_dp*y)*sqrt(qrho)
     n     = nint(6._dp + 72._dp*qrho)
     j     = 2*n+1
     xsum  = 1.0_dp/real(j,dp)
     ysum  = 0.0_dp
     do i=n, 1, -1
        j    = j - 2
        xaux = (xsum*xquad - ysum*yquad)/real(i,dp)
        ysum = (xsum*yquad + ysum*xquad)/real(i,dp)
        xsum = xaux + 1.0_dp/real(j,dp)
     end do
     u1   = -factor*(xsum*yabs + ysum*xabs) + 1.0_dp
     v1   =  factor*(xsum*xabs - ysum*yabs)
     !MODIFICATION: the exponent is corrected, yabs disappears

     daux = fac**2
     daux = exp(-daux)
     daux = daux**(q**2)

     !daux =  exp(-xquad-yabs**2)
     u2   =  daux*dcos(yquadmod)
     v2   = -daux*dsin(yquadmod)

     u    = u1*u2 - v1*v2
     v    = u1*v2 + v1*u2

  else

     !if (qrho.gt.1.0) then w(z) is evaluated using the laplace
     !continued fraction
     !nu is the minimum number of terms needed to obtain the required
     !accuracy
     !
     !if ((qrho.gt.0.085264).and.(qrho.lt.1.0)) then w(z) is evaluated
     !by a truncated taylor expansion, where the laplace continued fraction
     !is used to calculate the derivatives of w(z)
     !kapn is the minimum number of terms in the taylor expansion needed
     !to obtain the required accuracy
     !nu is the minimum number of terms of the continued fraction needed
     !to calculate the derivatives with the required accuracy

     if (qrho > 1.0) then
        h    = 0.0_dp
        kapn = 0
        qrho = sqrt(qrho)
        nu   = int(3._dp + (1442._dp/(26._dp*qrho+77._dp)))
     else
        qrho = (1-y)*sqrt(1-qrho)
        h    = 1.88*qrho
        h2   = 2._dp*h
        kapn = nint(7._dp  + 34._dp*qrho)
        nu   = nint(16._dp + 26._dp*qrho)
     endif

     b = (h > 0.0_dp)



     if (b) qlambda = h2**kapn

     rx = 0.0_dp
     ry = 0.0_dp
     sx = 0.0_dp
     sy = 0.0_dp

     do n=nu, 0, -1
        np1 = n + 1
        tx  = yabs + h + real(np1,dp)*rx
        ty  = xabs - real(np1,dp)*ry
        c   = 0.5_dp/(tx**2 + ty**2)
        rx  = c*tx
        ry  = c*ty
        if ((b).and.(n <= kapn)) then
           tx = qlambda + sx
           sx = rx*tx - ry*sy
           sy = ry*tx + rx*sy
           qlambda = qlambda/h2
        endif
        !print *,'nu,yabs,xabs,tx,ty,rx,ry',n,nu,yabs,xabs,tx,ty,rx,ry
     end do

     if (h == 0.0_dp) then
        u = factor*rx
        v = factor*ry
     else
        u = factor*sx
        v = factor*sy
     end if

     !MODIFICATION: here the exponent is added
     daux=exp(-yabs**2)
     u=daux*u
     v=daux*v

     !still do not now what happens to v in that case
     if (yabs == 0.0_dp) u = dexp(-xabs**2)

  end if

  !
  !  evaluation of w(z) in the other quadrants
  !
  if (real(jm,gp)+t < 0.0_dp) then

     if (a) then
        !no modification is needed
        u2    = 2._dp*u2
        v2    = 2._dp*v2
     else
        xquad =  -xquad

        !the following if-statement protects 2*exp(-z**2)
        !against overflow, taking into account the modification
        if (xquad-yabs**2 > rmaxexp) then
           flag = .true.
           return
        end if
        !daux = fac**2
        !daux = exp(-daux)
        !w1 = 2.0_dp*daux**(q**2)
        !avoid floating point exceptions
        daux = real(q,kind=8)*fac
        daux = daux*daux
        daux = dexp(-daux)
        w1 = 2.0_dp*daux    
        !w1 =  2*dexp(xquad-yabs**2)
        u2  =  w1*dcos(yquadmod)
        v2  = -w1*dsin(yquadmod)
        !print *,'w1,u2,v2',yquad,yquadmod,xquad,yabs**2,xquad-yabs**2,w1,u2,v2
     end if

     u = u2 - u
     v = v2 - v
     !print *,u,v
     if (q > 0) v = -v
  else
     if (q < 0) v = -v
  end if
END SUBROUTINE GcplxInt


subroutine GcplxInt2(alpha,m,q,p,t,u,v,flag)
  use module_base
  implicit none
  real(dp), intent(in) :: alpha,t
  integer, intent(in) :: q,m,p
  logical, intent(out) :: flag
  real(dp), intent(out) :: u,v
  !local variables
  real(dp), parameter :: factor=1.12837916709551257388_dp,rmaxreal=0.5e+154_dp
  real(dp), parameter :: rmaxexp=708.503061461606_dp !n(c) ,rmaxgoni=3.53711887601422e+15_dp
  real(dp), parameter :: pi=3.1415926535897932384_dp
  logical :: a,b
  integer :: j,n,i,kapn,nu,np1,multiple
  real(dp) :: xabs,yabs,x,y,qrho,xquad,yquad,xsum,ysum,xaux,u1,v1,u2,v2,daux,h,qlambda,h2
  real(dp) :: rx,ry,sx,sy,tx,ty,c,w1,xabsq,xi,yi,fac,yquadmod

  flag = .false.

  fac=pi/real(2*m,dp)/alpha

  !here xi and yi can be eliminated
  xi=fac*real(q,dp)
  yi=alpha*(real(p,dp)+t)
  xabs = dabs(xi)
  yabs = dabs(yi)
  x    = xabs/6.3_dp
  y    = yabs/4.4_dp

!     the following if-statement protects
!     qrho = (x**2 + y**2) against overflow
  if ((xabs > rmaxreal).or.(yabs > rmaxreal)) then
     flag = .true.
     print *,'ciao'
     return
  end if

  qrho = x**2 + y**2 !the result is identical for both terms

  xabsq = xabs**2
  xquad = xabsq - yabs**2
  yquad = 2*xabs*yabs
  !yquad to be passed to trigonometric functions (here t has to be added)
  !must subtract a arbitrary multiple of the period
  multiple=(q*p)/(2*m)
  multiple=q*p-multiple*2*m
  yquadmod = (real(multiple,dp)+t)/real(m,dp)*pi
  !yquadmod = real(modulo(abs(q*(p)),2*m),dp)/real(m,dp)*pi

  a = qrho < 0.085264_dp

  if (a) then
     !if (qrho < 0.085264d0) then the faddeeva-function is evaluated
     !using a power-series (abramowitz/stegun, equation (7.1.5), p.297)
     !n is the minimum number of terms needed to obtain the required
     !accuracy
     
     qrho  = (1._dp-0.85_dp*y)*sqrt(qrho)
     n     = nint(6._dp + 72._dp*qrho)
     j     = 2*n+1
     xsum  = 1.0_dp/real(j,dp)
     ysum  = 0.0_dp
     do i=n, 1, -1
        j    = j - 2
        xaux = (xsum*xquad - ysum*yquad)/real(i,dp)
        ysum = (xsum*yquad + ysum*xquad)/real(i,dp)
        xsum = xaux + 1.0_dp/real(j,dp)
     end do
     u1   = -factor*(xsum*yabs + ysum*xabs) + 1.0_dp
     v1   =  factor*(xsum*xabs - ysum*yabs)
     !MODIFICATION: the exponent is corrected, yabs disappears

     daux = fac**2
     daux = exp(-daux)
     daux = daux**(q**2)

     !daux =  exp(-xquad-yabs**2)
     u2   =  daux*dcos(yquadmod)
     v2   = -daux*dsin(yquadmod)

     u    = u1*u2 - v1*v2
     v    = u1*v2 + v1*u2 !can be eliminated

    !changing the sign of x*y implies yquad -> -yquad

  else

     !if (qrho.gt.1.0) then w(z) is evaluated using the laplace
     !continued fraction
     !nu is the minimum number of terms needed to obtain the required
     !accuracy
     !
     !if ((qrho.gt.0.085264).and.(qrho.lt.1.0)) then w(z) is evaluated
     !by a truncated taylor expansion, where the laplace continued fraction
     !is used to calculate the derivatives of w(z)
     !kapn is the minimum number of terms in the taylor expansion needed
     !to obtain the required accuracy
     !nu is the minimum number of terms of the continued fraction needed
     !to calculate the derivatives with the required accuracy

     if (qrho > 1.0) then
        h    = 0.0_dp
        kapn = 0
        qrho = sqrt(qrho)
        nu   = int(3._dp + (1442._dp/(26._dp*qrho+77._dp)))
     else
        qrho = (1-y)*sqrt(1-qrho)
        h    = 1.88*qrho
        h2   = 2._dp*h
        kapn = nint(7._dp  + 34._dp*qrho)
        nu   = nint(16._dp + 26._dp*qrho)
     endif

     b = (h > 0.0_dp)



     if (b) qlambda = h2**kapn

     rx = 0.0_dp
     ry = 0.0_dp
     sx = 0.0_dp
     sy = 0.0_dp

     do n=nu, 0, -1
        np1 = n + 1
        tx  = yabs + h + real(np1,dp)*rx
        ty  = xabs - real(np1,dp)*ry
        c   = 0.5_dp/(tx**2 + ty**2)
        rx  = c*tx
        ry  = c*ty
        if ((b).and.(n <= kapn)) then
           tx = qlambda + sx
           sx = rx*tx - ry*sy
           sy = ry*tx + rx*sy
           qlambda = qlambda/h2
        endif
        !print *,'nu,yabs,xabs,tx,ty,rx,ry',n,nu,yabs,xabs,tx,ty,rx,ry
     end do

     if (h == 0.0_dp) then
        u = factor*rx
        v = factor*ry
     else
        u = factor*sx
        v = factor*sy
     end if

     !MODIFICATION: here the exponent is added
     daux=exp(-yabs**2)
     u=daux*u
     v=daux*v

     !still do not now what happens to v in that case
     if (yabs == 0.0_dp) u = dexp(-xabs**2)

  end if

  !
  !  evaluation of w(z) in the other quadrants
  !
  if (p < 0) then

     if (a) then
        !no modification is needed
        u2    = 2._dp*u2
        v2    = 2._dp*v2
     else
        xquad =  -xquad

        !the following if-statement protects 2*exp(-z**2)
        !against overflow, taking into account the modification
        if (xquad-yabs**2 > rmaxexp) then
           print *,'ciao2',xquad,yabs
           flag = .true.
           return
        end if
        !daux = fac**2
        !daux = exp(-daux)
        !w1 = 2.0_dp*daux**(q**2)
        !avoid floating point exceptions
        daux = real(q,kind=8)*fac
        daux = daux*daux
        daux = dexp(-daux)
        w1 = 2.0_dp*daux    
        !w1 =  2*dexp(xquad-yabs**2)
        u2  =  w1*dcos(yquadmod)
        v2  = -w1*dsin(yquadmod)
        !print *,'w1,u2,v2',yquad,yquadmod,xquad,yabs**2,xquad-yabs**2,w1,u2,v2
     end if

     u = u2 - u
     v = v2 - v
     !print *,u,v
     if (q > 0) v = -v
  else
     if (q < 0) v = -v
  end if
END SUBROUTINE GcplxInt2


!> Call the routine which performs the interpolation in each direction
subroutine my_interpolate_and_transpose(t0,nphi,nrange,phi,ndat,nin,psi_in,nout,psi_out)
 use module_base
 implicit none
 integer, intent(in) :: nphi      !< number of sampling points of the ISF function (multiple of nrange)
 integer, intent(in) :: nrange    !< extension of the ISF domain in dimensionless units (even number)
 integer, intent(in) :: nin,nout  !< sizes of the input and output array in interpolating direction
 integer, intent(in) :: ndat      !< size of the array in orthogonal directions
 real(gp), intent(in) :: t0       !< shift in the interpolating direction in grid spacing units
 real(gp), dimension(nphi), intent(in) :: phi           !< interpolating scaling function array
 real(gp), dimension(nin,ndat), intent(in) :: psi_in    !< input wavefunction psifscf
 real(gp), dimension(ndat,nout), intent(out) :: psi_out !< input wavefunction psifscf
 !local variables
 character(len=*), parameter :: subname='interpolate_and_transpose'
 integer :: i_all,i_stat,nunit,m_isf,ish,ipos,i,j,l,ms,me
 real(gp) :: dt, tt
 real(gp), dimension(:), allocatable :: shf !< shift filter

 !assume for the moment that the grid spacing is constant

 m_isf=nrange/2
 !calculate the shift filter for the given t0
 allocate(shf(-m_isf:m_isf+ndebug),stat=i_stat )
 call memocc(i_stat,shf,'shf',subname)

 !number of points for a unit displacement
 nunit=nphi/nrange 
 !this should be number between -0.5 and 0.5
 dt=t0-nint(t0)
 !evaluate the shift
 ish=nint(real(nunit,gp)*dt)

print *,'start interpolation',ish,real(nunit,gp)*dt,dt

 if (ish<=0) then
   shf(-m_isf)=0.0_gp
 else
   shf(-m_isf)=phi(ish)  
 end if 
 ipos=ish

 do i=-m_isf+1,m_isf-1 !extremes excluded
   !position of the shifted argument in the phi array
   ipos=ipos+nunit
   shf(i)=phi(ipos)  
 end do

 if (ish<=0) then
   shf(m_isf)=phi(ipos+nunit)
 else
   shf(m_isf)=0.0_gp
 end if 

 !define the shift for output results
 ish=nint(t0)

print *,'start interpolation',ish,t0,shf

do i=-m_isf,m_isf
write(104,*)i,shf(i)
end do

 !apply the interpolating filter to the output
 do j=1,ndat
   psi_out(j,:)=0.0_gp
   do i=1,nin
     !here the boundary conditions have to be considered
      tt=0.0_gp
      ms=-min(m_isf,i-1)
      me=min(m_isf,nin-i)
      do l=ms,me
         tt=tt+shf(l)*psi_in(i+l,j)
      end do
!      tt=h*tt
      
      if (i+ish > 0 .and. i+ish < nout) psi_out(j,i+ish)=tt
!	      write(102,*)i+ish,psi_out(j,i+ish)
   end do
 end do



 i_all=-product(shape(shf))*kind(shf)
 deallocate(shf,stat=i_stat)
 call memocc(i_stat,i_all,'shf',subname)

end subroutine my_interpolate_and_transpose


subroutine scaling_function4b2B(itype,nd,nrange,a,x)
   use module_base
   implicit none
   !Arguments
   !Type of interpolating functions
   integer, intent(in) :: itype
   !Number of points: must be 2**nex
   integer, intent(in) :: nd
   integer, intent(out) :: nrange
   real(kind=8), dimension(0:nd), intent(out) :: a,x
   !Local variables
   character(len=*), parameter :: subname='scaling_function4b2B'
   real(kind=8), dimension(:), allocatable :: y
   integer :: i,nt,ni,i_all,i_stat  

   !Only itype=8,14,16,20,24,30,40,50,60,100
   select case(itype)
   case(8,14,16,20,24,30,40,50,60,100)
      !O.K.
   case default
      print *,"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100"
      stop
   end select
   !!$  write(unit=*,fmt="(1x,a,i0,a)") &
   !!$       "Use interpolating scaling functions of ",itype," order"

   !Give the range of the scaling function
   !from -itype to itype
   ni=2*itype
   nrange = ni
   allocate(y(0:nd+ndebug),stat=i_stat)
   call memocc(i_stat,y,'y',subname)

   ! plot scaling function
   call zero(nd+1,x)
   call zero(nd+1,y)
   nt=ni
   x(nt/2)=1.d0
   loop1: do
      nt=2*nt
      ! write(6,*) 'nd,nt',nd,nt
      select case(itype)
      case(8)
         stop
      case(14)
         stop
      case(16)
         call back_trans_16(nd,nt,x,y)
      case(20)
         stop
      case(24)
         stop
      case(30)
         stop
      case(40)
         stop
      case(50)
         stop
      case(60)
         stop
      case(100)
         stop
      end select

      do i=0,nt-1
         x(i)=y(i)
      end do
      if (nt.eq.nd) then
         exit loop1
      end if
   end do loop1

   !open (unit=1,file='scfunction',status='unknown')
   do i=0,nd
      a(i) = real(i*ni,kind=8)/real(nd,kind=8)-(.5d0*real(ni,kind=8)-1.d0)
      !write(1,*) a(i),x(i)
   end do
   !close(1)

   i_all=-product(shape(y))*kind(y)
   deallocate(y,stat=i_stat)
   call memocc(i_stat,i_all,'y',subname)
END SUBROUTINE scaling_function4b2B


!call the routine which performs the interpolation in each direction
subroutine morph_and_transpose(h,t0_field,nphi,nrange,phi,ndat,nin,psi_in,nout,psi_out)
 use module_base
 implicit none
 integer, intent(in) :: nphi !< number of sampling points of the ISF function (multiple of nrange)
 integer, intent(in) :: nrange !< extension of the ISF domain in dimensionless units (even number)
 integer, intent(in) :: nin,nout !< sizes of the input and output array in interpolating direction
 integer, intent(in) :: ndat !< size of the array in orthogonal directions
 real(gp), intent(in) :: h !< grid spacing in the interpolating direction
 real(gp), dimension(nin,ndat), intent(in) :: t0_field !< field of shifts to be applied for each point in grid spacing units
 real(gp), dimension(nphi), intent(in) :: phi !< interpolating scaling function array
 real(gp), dimension(nin,ndat), intent(in) :: psi_in !< input wavefunction psifscf
 real(gp), dimension(ndat,nout), intent(out) :: psi_out !< input wavefunction psifscf
 !local variables
 character(len=*), parameter :: subname='my_morph_and_transpose'
 real(gp), parameter  :: tol=1.e-14_gp
 integer :: i_all,i_stat,nunit,m_isf,ish,ipos,i,j,l,ms,me,k2,k1
 real(gp) :: dt, tt, t0_l, tt2, ksh1, ksh2, k, kold, alpha,diff
 real(gp), dimension(:), allocatable :: shf !< shift filter

 !assume for the moment that the grid spacing is constant
 !call f_malloc_routine_id(subname)
 m_isf=nrange/2

 !shf=f_malloc(bounds=(/-m_isf .to. m_isf/),id='shf')

 !calculate the shift filter for the given t0
 allocate(shf(-m_isf:m_isf+ndebug),stat=i_stat )
 call memocc(i_stat,shf,'shf',subname)

 !number of points for a unit displacement
 nunit=nphi/nrange 

!print *,'start interpolation',ish,t0,shf

!do i=-m_isf,m_isf
!write(104,*)i,shf(i)
!end do

 !apply the interpolating filter to the output
 do j=1,ndat
   psi_out(j,:)=0.0_gp
do i=1,nout

kold=-1000.0_gp
find_trans: do l=1,nin
     k=real(l,gp)+t0_field(l,j)
     if (k-real(i,gp) > tol) exit find_trans
     kold=k
end do find_trans
    
!idist=k-kold

!if (idist > 1) then
 !  print *,k,kold,l,t0_field(l,j)
  ! stop 'idist'

!if (l>=nin) l=nin

!print*,i,l,k,kold
! want to use either l or l-1 to give us point i - pick closest
     if (k-real(i,gp) < -kold+real(i,gp)) then
        ksh1=k-real(i,gp)
        ksh2=-kold+real(i,gp)
        k1=l
        k2=l-1
        !if (k2==0) k2=1

        if (k2==0) then
           k2=1
           ksh2=ksh1
        end if
        if (k1==nin+1) then
           k1=nin
           ksh1=ksh2
        end if
     else
        ksh1=-kold+real(i,gp)
        ksh2=k-real(i,gp)
        k1=l-1
        k2=l
        if (k1==0) then
           k1=1
           ksh1=ksh2
        end if
        if (k2==nin+1) then
           k2=nin
           ksh2=ksh1
        end if
     end if


!print*,i,k
!print*,ksh1,ksh2,ksh1+ksh2

if (ksh1==0.0_gp .or. k1==k2) then !otherwise already have exactly on point
tt2=0.0_gp
ksh2=1
ksh1=0
end if 

    alpha=ksh2/(ksh1+ksh2)
!alpha=1.0_gp
     t0_l=alpha*t0_field(k1,j)+(1.0_gp-alpha)*t0_field(k2,j)

   !    t0_l=t0_field(k1,j)-ksh2

!   if (abs(t0_l) > 1.0_gp) then
!     dt=t0_l-nint(t0_l*(2.0_gp))/2
!else
   dt=t0_l-nint(t0_l)
   
   diff=real(i,gp)-(k1+t0_l)   

   if (abs(diff - dt) < abs(diff+dt)) dt=-dt
   !dt=-real(i,gp)+(k1+t0_l)
!end if
!     if (k2<k1.and.t0_field(k1,j)<0) dt=-dt
!     if (k2>k1.and.t0_field(k2,j)<0) dt=-dt

   !   if (abs(dt) > 0.5d

     !evaluate the shift
     ish=nint(real(nunit,gp)*dt)

     if (ish<=0) then
       shf(-m_isf)=0.0_gp
     else
       shf(-m_isf)=phi(ish)  
     end if 
     ipos=ish

     print *,'j',j,'start',k1,'shift',t0_l,'end',i,'dt',dt,'difference',real(i,gp)-(k1+t0_l)

     do l=-m_isf+1,m_isf-1 !extremes excluded
       !position of the shifted argument in the phi array
       ipos=ipos+nunit
       shf(l)=phi(ipos)  
     end do

     if (ish<=0) then
       shf(m_isf)=phi(ipos+nunit)
     else
       shf(m_isf)=0.0_gp
     end if 

     !define the shift for output results
     !ish=nint(t0_l)
     !ish=0

     !here the boundary conditions have to be considered
      tt=0.0_gp
      ms=-min(m_isf,k1-1)
      me=min(m_isf,nin-k1)
      do l=ms,me
         tt=tt+shf(l)*psi_in(k1+l,j)
      end do



!print '(a,3i4,4(1pe15.7),4i4)','i,k1,k2,ksh1,ksh2,dt,t0_l,nint(t0_l),j',i,k1,k2,ksh1,ksh2,dt,t0_l,nint(t0_l),j,&
!     nint(t0_l*2.0_gp),nint(t0_l*2.0_gp)/2 !,tt2,(abs(ksh2)*tt+abs(ksh1)*tt2)/(abs(ksh1)+abs(ksh2))
write(104,*) i,j,k1,k2,ksh1,ksh2,tt,t0_l,nint(t0_l)!+sign(t0_l,1.0_gp)*1.0_gp*epsilon(1.0_gp))
      if (i > 0 .and. i < nout) psi_out(j,i)=tt!=(abs(ksh2)*tt+abs(ksh1)*tt2)/(abs(ksh1)+abs(ksh2))
	if (i > 0 .and. i < nout)       write(102,*)i,psi_out(j,i)



end do

end do

! call f_free(shf)
! call f_malloc_free_routine()

 i_all=-product(shape(shf))*kind(shf)
 deallocate(shf,stat=i_stat)
 call memocc(i_stat,i_all,'shf',subname)

end subroutine morph_and_transpose
