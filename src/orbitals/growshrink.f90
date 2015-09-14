!> @file
!!  Synthesis convolution routines
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>  Apply synthesis wavelet transformation
subroutine comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3&
     ,w2,w1,xc,xf,y,ibyz_c,ibzxx_c,ibxxyy_c,&
     ibyz_f,ibzxx_f,ibxxyy_f)
  use module_defs, only: wp
  implicit none
  integer,intent(in)::n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c
  integer, dimension(2,0:n3,-14:2*n1+16), intent(in) :: ibzxx_c
  integer, dimension(2,-14:2*n1+16,-14:2*n2+16), intent(in) :: ibxxyy_c
  integer, dimension(2,nfl2:nfu2,nfl3:nfu3), intent(in) :: ibyz_f
  integer, dimension(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16), intent(in) :: ibzxx_f
  integer, dimension(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16), intent(in) :: ibxxyy_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: xc
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: xf
  real(wp), dimension(max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
       (2*n1+31)*(n2+1)*(n3+1))), intent(inout) :: w1 !work
  real(wp), dimension(max((n3+1)*(2*n1+31)*(2*n2+31),&
       2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))), intent(inout) :: w2 ! work
  real(wp), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(inout) :: y

  call comb_grow_c(n1,n2,n3,w1,w2,xc,y,ibyz_c,ibzxx_c,ibxxyy_c)

  call comb_grow_tree(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       w1,w2,xf,y,ibyz_f,ibzxx_f,ibxxyy_f)                

END SUBROUTINE comb_grow_all



!>   In 3d,            
!!   Applies synthesis wavelet transformation 
!!   then convolves with magic filter
!!   the size of the data is allowed to grow
!!   The input array x is not overwritten
!!   However, the output array y contains nonphysical values
!!   outside of the localization region
!!   that remain from the first comb_grow
subroutine comb_grow_c(n1,n2,n3,w1,w2,x,y,ibyz,ibzxx,ibxxyy)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz
  integer, dimension(2,0:n3,-14:2*n1+16), intent(in) :: ibzxx
  integer, dimension(2,-14:2*n1+16,-14:2*n2+16), intent(in) ::  ibxxyy
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n2,0:n3,-14:2*n1+16), intent(inout) :: w1
  real(wp), dimension(0:n3,-14:2*n1+16,-14:2*n2+16), intent(inout) :: w2
  real(wp), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(inout) :: y

  ! i1,i2,i3 -> i2,i3,I1
  call  comb_rot_grow_loc_square(n1,n2,n3,x,w1,ibyz,ibzxx) 

  ! i2,i3,I1 -> i3,I1,I2
  call  comb_rot_grow_loc_square(n2,n3,2*n1+30,w1,w2,ibzxx,ibxxyy) 

  ! i3,I1,I2  -> I1,I2,I3
  call  comb_rot_grow_loc_square_3(n3,2*n1+30,2*n2+30,w2,y,ibxxyy) 

END SUBROUTINE comb_grow_c


!>   In 3d,
!!   Applies synthesis wavelet transformation 
!!   then convolves with magic filter
!!   the size of the data is allowed to grow
!!   The input array x is not overwritten
subroutine comb_grow_tree(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3&
     ,w1,w2,x,y,ibyz,ibzxx,ibxxyy)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, dimension(2,nfl2:nfu2,nfl3:nfu3), intent(in) :: ibyz
  integer, dimension(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16), intent(in) :: ibzxx
  integer, dimension(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16), intent(in) :: ibxxyy
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x
  real(wp), dimension(4,nfl2:nfu2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16), intent(inout) :: w1
  real(wp), dimension(2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16), intent(inout) :: w2
  real(wp), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(inout) :: y
  !local variables
  !n(c) integer :: m1,m2,m3,nt

  !n(c) m1=nfu1-nfl1
  !n(c) m2=nfu2-nfl2
  !n(c) m3=nfu3-nfl3

  ! i1,i2,i3 -> i2,i3,I1
  !n(c) nt=(nfu2-nfl2+1)*(nfu3-nfl3+1)
  call comb_rot_grow_loc_1(nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,w1,ibyz,ibzxx) 

  ! i2,i3,I1 -> i3,I1,I2
  !n(c) nt=(nfu3-nfl3+1)*(2*m1+31)
  call comb_rot_grow_loc_2(nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,ibzxx,ibxxyy) 

  ! i3,I1,I2  -> I1,I2,I3: add the result to y
  call comb_rot_grow_loc_3(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w2,y,ibxxyy)

END SUBROUTINE comb_grow_tree


!>   In one dimension,    
!!   with optimised cycles
!!   Applies synthesis wavelet transformation 
!!   then convolves with magic filter
!!   the size of the data is allowed to grow
subroutine comb_rot_grow_loc_1(nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y,ib,ib2)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, dimension(2,nfl2:nfu2,nfl3:nfu3), intent(in) :: ib
  integer, dimension(2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16), intent(in) :: ib2 
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x
  real(wp), dimension(2,2,nfl2:nfu2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16), intent(inout) :: y
  !local variables
  !integer ncount0,ncount1,ncount_rate,ncount_max,nflop
  !real(kind=8) tel
  integer :: l2,l3,i,t,l1
  real(wp) y2i__11,y2i__21,y2i1_11,y2i1_21
  real(wp) y2i__12,y2i__22,y2i1_12,y2i1_22

  include 'v_17.inc'

  !    open(unit=20,file='tree.flop')

  !    nflop=0
  !    do l2=nfl2,nfu2
  !        do l3=nfl3,nfu3
  !            if (ib(2,l2,l3).ge.ib(1,l2,l3)) nflop=nflop+(ib(2,l2,l3)-ib(1,l2,l3)+1)*31*2*7
  !        enddo
  !    enddo

  !    call system_clock(ncount0,ncount_rate,ncount_max)

  !   y=0._wp

!$omp parallel default (private) shared(x,y,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)&
!$omp shared(ib,ib2,fil2)

  !$omp do 
  do l1=-14+2*nfl1,2*nfu1+16
     do l3=nfl3,nfu3
        y(:,:,ib2(1,l3,l1):ib2(2,l3,l1),l3,l1)=0._wp
     enddo
  enddo
  !$omp enddo

  !$omp barrier
  
  !$omp do
  do l3=nfl3,nfu3
     do l2=nfl2,nfu2

        if (ib(1,l2,l3).le.ib(2,l2,l3)) then

           y(1,1,l2,l3,2*ib(2,l2,l3)+16)=    fil2(16,2)*x(1,ib(2,l2,l3),l2,l3)
           y(2,1,l2,l3,2*ib(2,l2,l3)+16)=&
                fil2(16,1)*x(2,ib(2,l2,l3),l2,l3)+fil2(16,2)*x(3,ib(2,l2,l3),l2,l3)

           y(1,2,l2,l3,2*ib(2,l2,l3)+16)=&
                fil2(16,1)*x(4,ib(2,l2,l3),l2,l3)+fil2(16,2)*x(5,ib(2,l2,l3),l2,l3)
           y(2,2,l2,l3,2*ib(2,l2,l3)+16)=&
                fil2(16,1)*x(6,ib(2,l2,l3),l2,l3)+fil2(16,2)*x(7,ib(2,l2,l3),l2,l3)

           do i=ib(1,l2,l3)-7,ib(2,l2,l3)+7 
              y2i__11=0._wp
              y2i__21=0._wp
              y2i__12=0._wp
              y2i__22=0._wp

              y2i1_11=0._wp
              y2i1_21=0._wp
              y2i1_12=0._wp
              y2i1_22=0._wp
              do t=max(i-8,ib(1,l2,l3)),min(i+7,ib(2,l2,l3))
                 y2i__11=y2i__11                               +fil2(2*(i-t)  ,2)*x(1,t,l2,l3)
                 y2i__21=y2i__21+fil2(2*(i-t)  ,1)*x(2,t,l2,l3)+fil2(2*(i-t)  ,2)*x(3,t,l2,l3)
                 y2i__12=y2i__12+fil2(2*(i-t)  ,1)*x(4,t,l2,l3)+fil2(2*(i-t)  ,2)*x(5,t,l2,l3)
                 y2i__22=y2i__22+fil2(2*(i-t)  ,1)*x(6,t,l2,l3)+fil2(2*(i-t)  ,2)*x(7,t,l2,l3)

                 y2i1_11=y2i1_11                               +fil2(2*(i-t)+1,2)*x(1,t,l2,l3)
                 y2i1_21=y2i1_21+fil2(2*(i-t)+1,1)*x(2,t,l2,l3)+fil2(2*(i-t)+1,2)*x(3,t,l2,l3)
                 y2i1_12=y2i1_12+fil2(2*(i-t)+1,1)*x(4,t,l2,l3)+fil2(2*(i-t)+1,2)*x(5,t,l2,l3)
                 y2i1_22=y2i1_22+fil2(2*(i-t)+1,1)*x(6,t,l2,l3)+fil2(2*(i-t)+1,2)*x(7,t,l2,l3)
              enddo
              y(1,1,l2,l3,2*i  )=y2i__11
              y(2,1,l2,l3,2*i  )=y2i__21
              y(1,2,l2,l3,2*i  )=y2i__12
              y(2,2,l2,l3,2*i  )=y2i__22

              y(1,1,l2,l3,2*i+1)=y2i1_11
              y(2,1,l2,l3,2*i+1)=y2i1_21
              y(1,2,l2,l3,2*i+1)=y2i1_12
              y(2,2,l2,l3,2*i+1)=y2i1_22
           enddo
        endif

     enddo
  enddo
  !$omp enddo
  !$omp end parallel

  !    call system_clock(ncount1,ncount_rate,ncount_max)
  !    tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !
  !    write(20,*) tel, 1.d-6*nflop/tel
END SUBROUTINE comb_rot_grow_loc_1


!>   In one dimension,    
!!   with optimised cycles
!!   Applies synthesis wavelet transformation 
!!   then convolves with magic filter
!!   the size of the data is allowed to grow
subroutine comb_rot_grow_loc_2(nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y,ib,ib2)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, dimension(2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16), intent(in) :: ib 
  integer, dimension(2,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16), intent(in) :: ib2
  real(wp), dimension(2,2,nfl2:nfu2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16), intent(in) :: x
  real(wp), dimension(2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16), intent(inout) :: y
  !integer ncount0,ncount1,ncount_rate,ncount_max,nflop
  !real(kind=8) :: tel
  integer l1,l3,i,t,l2
  real(wp) y2i__1,y2i__2,y2i1_1,y2i1_2

  include 'v_17.inc'

  !    open(unit=20,file='tree.flop')
  !    nflop=0
  !    do l3=nfl3,nfu3
  !        do l1=-14+2*nfl1,2*nfu1+16
  !            if (ib(2,l3,l1).ge.ib(1,l3,l1)) nflop=nflop+(ib(2,l3,l1)-ib(1,l3,l1)+1)
  !        enddo
  !    enddo
  !    nflop=nflop*31*2*4

  !    call system_clock(ncount0,ncount_rate,ncount_max)

  !     y=0._wp


!$omp parallel default (private) shared(x,y,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)&
!$omp shared(ib,ib2,fil2)

  !$omp do
  do l2=-14+2*nfl2,2*nfu2+16
     do l1=-14+2*nfl1,2*nfu1+16
        y(:,ib2(1,l1,l2):ib2(2,l1,l2),l1,l2)=0._wp
     enddo
  enddo
  !$omp enddo

  !$omp barrier

  !$omp do
  do l1=-14+2*nfl1,2*nfu1+16
     do l3=nfl3,nfu3

        if (ib(1,l3,l1).le.ib(2,l3,l1)) then
           y(1,l3,l1,2*ib(2,l3,l1)+16)=&
                fil2(16,1)*x(1,1,ib(2,l3,l1),l3,l1)+fil2(16,2)*x(2,1,ib(2,l3,l1),l3,l1)
           y(2,l3,l1,2*ib(2,l3,l1)+16)=&
                fil2(16,1)*x(1,2,ib(2,l3,l1),l3,l1)+fil2(16,2)*x(2,2,ib(2,l3,l1),l3,l1)

           do i=ib(1,l3,l1)-7,ib(2,l3,l1)+7 
              y2i__1=0._wp
              y2i__2=0._wp
              y2i1_1=0._wp
              y2i1_2=0._wp
              do t=max(i-8,ib(1,l3,l1)),min(i+7,ib(2,l3,l1))
                 y2i__1=y2i__1+fil2(2*(i-t)  ,1)*x(1,1,t,l3,l1)+fil2(2*(i-t)  ,2)*x(2,1,t,l3,l1)
                 y2i__2=y2i__2+fil2(2*(i-t)  ,1)*x(1,2,t,l3,l1)+fil2(2*(i-t)  ,2)*x(2,2,t,l3,l1)
                 y2i1_1=y2i1_1+fil2(2*(i-t)+1,1)*x(1,1,t,l3,l1)+fil2(2*(i-t)+1,2)*x(2,1,t,l3,l1)
                 y2i1_2=y2i1_2+fil2(2*(i-t)+1,1)*x(1,2,t,l3,l1)+fil2(2*(i-t)+1,2)*x(2,2,t,l3,l1)
              enddo
              y(1,l3,l1,2*i  )=y2i__1
              y(2,l3,l1,2*i  )=y2i__2
              y(1,l3,l1,2*i+1)=y2i1_1
              y(2,l3,l1,2*i+1)=y2i1_2
           enddo
        endif

     enddo
  enddo
  !$omp enddo

  !$omp end parallel

  !    call system_clock(ncount1,ncount_rate,ncount_max)
  !    tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !    write(20,*) tel, 1.d-6*nflop/tel
END SUBROUTINE comb_rot_grow_loc_2


!>   In 3d,            
!!   Applies the magic filter transposed, then analysis wavelet transformation.
!!   The size of the data is forced to shrink
!!   The input array y is not overwritten
subroutine comb_shrink(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,y,&
     ibxy_c,ibzzx_c,ibyyzz_c,ibxy_f,ibzzx_f,ibyyzz_f,xc,xf)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c
  integer, dimension(2,-14:2*n3+16,0:n1), intent(in) :: ibzzx_c
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in) :: ibyyzz_c
  integer, dimension(2,nfl1:nfu1,nfl2:nfu2), intent(in) :: ibxy_f
  integer, dimension(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1), intent(in) :: ibzzx_f
  integer, dimension(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16), intent(in) :: ibyyzz_f
  real(wp), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(in) :: y
  real(wp), dimension(max(2*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31)*(nfu1-nfl1+1),&
       (2*n2+31)*(2*n3+31)*(n1+1))), intent(inout) :: w1
  real(wp), dimension(max(4*(2*(nfu3-nfl3)+31)*(nfu1-nfl1+1)*(nfu2-nfl2+1),&
       (2*n3+31)*(n1+1)*(n2+1))), intent(inout) :: w2
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: xc
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: xf
  !    perform the combined transform    
  call comb_shrink_loc_c(0,n1,0,n2,0,n3,w1,w2,y,xc,1,1,1,&
       ibxy_c,ibzzx_c,ibyyzz_c) ! for scfunctions
  !    for wavelets:
  call comb_shrink_loc_f(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,y,xf,&
       ibxy_f,ibzzx_f,ibyyzz_f)

END SUBROUTINE comb_shrink


!>   In 3d,            
!!   Applies the magic filter transposed, then analysis wavelet transformation.
!!   The output is only the l1,l2,l3 wavelet component
!!   The size of the data is forced to shrink
!!   The input array y is not overwritten
subroutine comb_shrink_loc_f(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,y,x,&
     ibxy,ibzzx,ibyyzz)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, dimension(2,nfl1:nfu1,nfl2:nfu2), intent(in) :: ibxy
  integer, dimension(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1), intent(in) :: ibzzx
  integer, dimension(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16), intent(in) :: ibyyzz
  real(wp), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(in) :: y ! input
  real(wp), dimension(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16,nfl1:nfu1), intent(inout) :: w1
  real(wp), dimension(4,-14+2*nfl3:2*nfu3+16,nfl1:nfu1,nfl2:nfu2), intent(inout) :: w2
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: x
  !local variables
  integer :: m1,m2,m3,nt

  m1=nfu1-nfl1
  m2=nfu2-nfl2
  m3=nfu3-nfl3

  ! I1,I2,I3 -> I2,I3,i1
  call comb_rot_shrink_loc_1(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,y,w1,ibyyzz)

  ! I2,I3,i1 -> I3,i1,i2
  nt=(2*m3+31)*(m1+1)
  call comb_rot_shrink_loc_2(nt,w1,w2,nfl2,nfu2,ibzzx)
  ! I3,i1,i2 -> i1,i2,i3

  nt=(m1+1)*(m2+1)
  call comb_rot_shrink_loc_3(nt,w2,x,nfl3,nfu3,ibxy)

END SUBROUTINE comb_shrink_loc_f


!>   In 3d,            
!!   Applies the magic filter transposed, then analysis wavelet transformation.
!!   The output is only the l1,l2,l3 wavelet component
!!   The size of the data is forced to shrink
!!   The input array y is not overwritten
subroutine comb_shrink_loc_c(nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,y,x,l1,l2,l3,&
     ibxy,ibzzx,ibyyzz)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,l1,l2,l3
  integer, dimension(2,nfl1:nfu1,nfl2:nfu2), intent(in) :: ibxy  
  integer, dimension(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1), intent(in) :: ibzzx 
  integer, dimension(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16), intent(in) :: ibyyzz
  real(wp), dimension(-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16),&
       intent(in) :: y!input
  real(wp), dimension(-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16,nfl1:nfu1), intent(inout) :: w1
  real(wp), dimension(-14+2*nfl3:2*nfu3+16,nfl1:nfu1,nfl2:nfu2), intent(inout) :: w2
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: x!output
  !local variables
  integer :: m1,m2,m3,nt

  m1=nfu1-nfl1
  m2=nfu2-nfl2
  m3=nfu3-nfl3

  ! I1,I2,I3 -> I2,I3,i1
  nt=(2*m2+31)*(2*m3+31)
  call comb_rot_shrink_loc(nt,y,w1,l1,nfl1,nfu1,ibyyzz)

  ! I2,I3,i1 -> I3,i1,i2
  nt=(2*m3+31)*(m1+1)
  call comb_rot_shrink_loc(nt,w1,w2,l2,nfl2,nfu2,ibzzx)

  ! I3,i1,i2 -> i1,i2,i3
  nt=(m1+1)*(m2+1)
  call comb_rot_shrink_loc(nt,w2,x,l3,nfl3,nfu3,ibxy)

END SUBROUTINE comb_shrink_loc_c


!>   In one dimension,    
!!   Applies the magic filter transposed, then analysis wavelet transformation.
!!   The size of the data is forced to shrink
subroutine comb_rot_shrink_loc_3(ndat,x,y,nfl,nfu,ib)
  use module_defs, only: wp
  implicit none
  integer, parameter :: lowfil2=-14,lupfil2=16
  integer, intent(in) :: ndat,nfl,nfu
  integer, dimension(2,ndat), intent(in) :: ib
  real(wp), dimension(2,2,lowfil2+2*nfl:2*nfu+lupfil2,ndat), intent(in) :: x
  real(wp), dimension(7,ndat,nfl:nfu), intent(inout) :: y
  !local variables
  integer :: i,j,l
  real(wp) :: ci112,ci121,ci122,ci211,ci212,ci221,ci222
  include 'v.inc'

  !nflop=0
  !open(unit=20,file='long.flop')
  !call system_clock(ncount0,ncount_rate,ncount_max)

  !$omp parallel do default(private) shared(ndat,nfl,nfu) &
  !$omp shared(x,y,ib,fil2)
  do j=1,ndat
     do i=ib(1,j),ib(2,j)
        ci112=0._wp
        ci121=0._wp
        ci122=0._wp
        ci211=0._wp
        ci212=0._wp
        ci221=0._wp
        ci222=0._wp

        !nflop=nflop+(lupfil2-lowfil2+1)*2*7
        do l=lowfil2+2*i,lupfil2+2*i
           ci112=ci112+fil2(l-2*i,2)*x(1,1,l,j)
           ci121=ci121+fil2(l-2*i,1)*x(1,2,l,j)
           ci122=ci122+fil2(l-2*i,2)*x(1,2,l,j)
           ci211=ci211+fil2(l-2*i,1)*x(2,1,l,j)
           ci212=ci212+fil2(l-2*i,2)*x(2,1,l,j)
           ci221=ci221+fil2(l-2*i,1)*x(2,2,l,j)
           ci222=ci222+fil2(l-2*i,2)*x(2,2,l,j)
        enddo
        y(1,j,i)=ci211
        y(2,j,i)=ci121    
        y(3,j,i)=ci221    
        y(4,j,i)=ci112    
        y(5,j,i)=ci212    
        y(6,j,i)=ci122    
        y(7,j,i)=ci222    
     enddo
  enddo
  !$omp end parallel do

  !call system_clock(ncount1,ncount_rate,ncount_max)
  !tel=dble(ncount1-ncount0)/dble(ncount_rate)
  ! write(20,*) tel, 1.d-6*nflop/tel

END SUBROUTINE comb_rot_shrink_loc_3
