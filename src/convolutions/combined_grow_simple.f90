!> @file
!! Non optimized version of combined synthesis + magic filters routines
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> In one dimension,    
!! with optimised cycles
!! Applies synthesis wavelet transformation 
!! then convolves with magic filter
!! then adds the result to y.
!! The size of the data is allowed to grow
subroutine  comb_rot_grow_loc_3(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y,ib)
  use module_base
  implicit none
  integer,intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,n1,n2,n3
  integer, dimension(2,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16), intent(in) :: ib
  real(wp), dimension(2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16), intent(in) :: x
  real(wp), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(out) :: y
  !local variables
!!    integer :: ncount0,ncount1,ncount_rate,ncount_max,nflop
!!    real(kind=8) :: tel
  integer :: l1,l2,i,t !n(c) l1_0, l1_1, ll1
  real(wp) :: y2i,y2i1 !n(c) y2i__0,y2i__1,y2i1_0,y2i1_1
  include 'v_17.inc'

  !    open(unit=20,file='tree.flop')

  !    nflop=0
  !    do l1=-14+2*nfl1,2*nfu1+16
  !        do l2=-14+2*nfl2,2*nfu2+16
  !            if (ib(2,l1,l2).ge.ib(1,l1,l2)) nflop=nflop+(ib(2,l1,l2)-ib(1,l1,l2)+1)*31*2*2
  !        enddo
  !    enddo

  !    call system_clock(ncount0,ncount_rate,ncount_max)

  do l2=-14+2*nfl2,2*nfu2+16
     do l1=-14+2*nfl1,2*nfu1+16
        if (ib(1,l1,l2).le.ib(2,l1,l2)) then
           y(l1,l2,2*ib(2,l1,l2)+16)=y(l1,l2,2*ib(2,l1,l2)+16)+&
                fil2(16,1)*x(1,ib(2,l1,l2),l1,l2)+fil2(16,2)*x(2,ib(2,l1,l2),l1,l2)

           do i=ib(1,l1,l2)-7,ib(2,l1,l2)+7 
              y2i=0.d0
              y2i1=0.d0
              do t=max(i-8,ib(1,l1,l2)),min(i+7,ib(2,l1,l2))
                 y2i=y2i+fil2(2*(i-t),1)*x(1,t,l1,l2)+fil2(2*(i-t),2)*x(2,t,l1,l2)
                 y2i1=y2i1+fil2(2*(i-t)+1,1)*x(1,t,l1,l2)+fil2(2*(i-t)+1,2)*x(2,t,l1,l2)
              enddo
              y(l1,l2,2*i  )=y(l1,l2,2*i  )+y2i
              y(l1,l2,2*i+1)=y(l1,l2,2*i+1)+y2i1
           enddo
        endif
     enddo
  enddo

  !    call system_clock(ncount1,ncount_rate,ncount_max)
  !    tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !
  !    write(20,*) tel, 1.d-6*nflop/tel
END SUBROUTINE comb_rot_grow_loc_3


!> In one dimension,
!! with unoptimized cycles   
!! Applies synthesis wavelet transformation 
!! then convolves with magic filter
!!  the size of the data is allowed to grow
subroutine comb_rot_grow_loc_square_1(n1,n2,n3,x,y,ib,ib2,loczero)
  use module_base
  implicit none
  logical, intent(in) :: loczero
  integer, intent(in) :: n1,n2,n3
  integer, dimension(2,0:n2,0:n3), intent(in) :: ib
  integer, dimension(2,0:n3,-14:2*n1+16), intent(in) :: ib2
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n2,0:n3,-14:2*n1+16), intent(out) :: y
  !local variables
!!    integer ncount0,ncount1,ncount2,ncount_rate,ncount_max,nflop
!!    real(kind=8) tel,t0,t1
  integer i,t,l2,l3 !n(c) ,l1
  !n(c) integer ll1,ll3,l10,l11,l30,l31,ll2,l21
  real(wp) y2i,y2i1 !n(c) y2i__11, y2i__12, y2i1_11, y2i1_12, y2i__21, y2i__22, y2i1_21, y2i1_22
  include 'v_17.inc'

  !    open(unit=10,file='zero.square')

  !    nflop=0
  !
  !    do l2=0,n2
  !        do l3=0,n3
  !            if (ib(2,l2,l3).ge.ib(1,l2,l3)) nflop=nflop+(ib(2,l2,l3)-ib(1,l2,l3)+1)*31*2
  !        enddo
  !    enddo
  !    call system_clock(ncount0,ncount_rate,ncount_max)

  !   initialize the y array with zeroes
  !   but only inside the region defined by ib2 array
  !   which is the ib array for the next step

  if (loczero) call make_loczero(n1,n2,n3,ib2,y)

  !    call system_clock(ncount1,ncount_rate,ncount_max)
  !    
  !    t0=dble(ncount1-ncount0)/dble(ncount_rate)

  !    the actual convolution

  do l3=0,n3
     do l2=0,n2

        if (ib(1,l2,l3).le.ib(2,l2,l3)) then

           y(l2,l3,2*ib(2,l2,l3)+16)=fil2(16,1)*x(ib(2,l2,l3),l2,l3)

           do i=ib(1,l2,l3)-7,ib(2,l2,l3)+7 
              y2i=0.d0
              y2i1=0.d0
              do t=max(i-8,ib(1,l2,l3)),min(i+7,ib(2,l2,l3))
                 y2i=y2i+fil2(2*(i-t),1)*x(t  ,l2,l3)
                 y2i1=y2i1+fil2(2*(i-t)+1,1)*x(t  ,l2,l3)
              enddo
              y(l2,l3,2*i)=y2i
              y(l2,l3,2*i+1)=y2i1
           enddo

        endif
     enddo
  enddo
  !    call system_clock(ncount2,ncount_rate,ncount_max)
  !    t1=dble(ncount2-ncount1)/dble(ncount_rate)
  !    tel=dble(ncount2-ncount0)/dble(ncount_rate)
  !
  !    write(10,'(3f10.3,f10.0)') t0,t1,tel, 1.d-6*nflop/tel
END SUBROUTINE comb_rot_grow_loc_square_1


!>  initialize the y array with zeroes
!!  but only inside the region defined by ib2 array
!!  which is the ib array for the next step
subroutine make_loczero(n1,n2,n3,ib2,y)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer, dimension(2,0:n3,-14:2*n1+16), intent(in) :: ib2
  real(wp), dimension(0:n2,0:n3,-14:2*n1+16), intent(out) :: y
  !local variables
  integer :: l1,l2,l3
    
        do l1=-14,2*n1+16
            do l3=0,n3
                do l2=ib2(1,l3,l1),ib2(2,l3,l1)
                    y(l2,l3,l1)=0.0_wp
                enddo
            enddo
        enddo
END SUBROUTINE make_loczero
