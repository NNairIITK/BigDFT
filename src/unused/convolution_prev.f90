!> @file
!!  Combined convolution routines
!! @deprecated
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> w2 and w1 are switched because from shrink convention we go
!! to grow convention
subroutine comb_grow_all_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3&
     ,w2,w1,xc,xf,y,ibyz_c,ibzxx_c,ibxxyy_c,&
     ibyz_f,ibzxx_f,ibxxyy_f)
  implicit none
  
  integer :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  
  real(kind=8), intent(in) :: xc(0:n1,0:n2,0:n3),  xf(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)! input
  real(kind=8), intent(out) :: y(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)! output
  
  real(kind=8) :: w1(4,nfl2:nfu2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16)!work
  !   real(kind=8) w2(0:n3,-14:2*n1+16,-14:2*n2+16) ! work
  ! real(kind=8) w2(2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16)
  real(kind=8) :: w2(max((n3+1)*(2*n1+31)*(2*n2+31),&
       2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))) ! work
  
  integer, intent(in) :: ibyz_c(2,0:n2,0:n3)
  integer, intent(in) :: ibzxx_c(2,0:n3,-14:2*n1+16)
  integer, intent(in) :: ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)
  
  integer, intent(in) :: ibyz_f(2,nfl2:nfu2,nfl3:nfu3) 
  integer, intent(in) :: ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
  integer, intent(in) :: ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)
  
  call comb_grow_c_prev(n1,n2,n3,w2,xc,y,ibyz_c,ibzxx_c,ibxxyy_c)
  
  call comb_grow_tree_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3&
       ,w1,w2,xf,y,ibyz_f,ibzxx_f,ibxxyy_f)    
  
END SUBROUTINE comb_grow_all_prev


!> In 3d,   
!! Applies synthesis wavelet transformation 
!! then convolves with magic filter
!!  the size of the data is allowed to grow
!! The input array x is not overwritten
subroutine comb_grow_tree_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3&
     ,w1,w2,x,y,ibyz,ibzxx,ibxxyy)

  implicit real(kind=8) (a-h,o-z)

  real(kind=8) :: x(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
  real(kind=8) :: w1(4,nfl2:nfu2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16)
  real(kind=8) :: w2(2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16)
  real(kind=8) :: y(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)
  integer :: ibyz(2,nfl2:nfu2,nfl3:nfu3)
  integer :: ibzxx(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
  integer :: ibxxyy(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)

  m1=nfu1-nfl1
  m2=nfu2-nfl2
  m3=nfu3-nfl3

  ! i1,i2,i3 -> i2,i3,I1
  nt=(nfu2-nfl2+1)*(nfu3-nfl3+1)
  call comb_rot_grow_loc_1_prev(nfl1,nfu1,nt,x,w1,ibyz) 
  
  ! i2,i3,I1 -> i3,I1,I2
  nt=(nfu3-nfl3+1)*(2*m1+31)
  call comb_rot_grow_loc_2_prev(nfl2,nfu2,nt,w1,w2,ibzxx) 
  
  ! i3,I1,I2  -> I1,I2,I3: add the result to y
  call comb_rot_grow_loc_3_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w2,y,ibxxyy)
  
END SUBROUTINE comb_grow_tree_prev


!> In 3d,   
!! Applies synthesis wavelet transformation 
!! then convolves with magic filter
!!  the size of the data is allowed to grow
!! The input array x is not overwritten
!! However, the output array y contains nonphysical values
!! outside of the localization region
!! that remain from the first comb_grow
subroutine comb_grow_c_prev(n1,n2,n3,ww,x,y,ibyz,ibzxx,ibxxyy)
  
  implicit real(kind=8) (a-h,o-z)
  real(kind=8) :: x(0:n1,0:n2,0:n3)
  real(kind=8) :: ww(0:n3,-14:2*n1+16,-14:2*n2+16) ! work
  real(kind=8) :: y(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)
  integer :: ibyz(2,0:n2,0:n3)
  integer :: ibzxx(2,0:n3,-14:2*n1+16)
  integer :: ibxxyy(2,-14:2*n1+16,-14:2*n2+16)
  
  ! i1,i2,i3 -> i2,i3,I1
  nt=(n2+1)*(n3+1)
  call comb_rot_grow_loc_prev(0,n1,nt,x,y,1,ibyz) 
  
  ! i2,i3,I1 -> i3,I1,I2
  nt=(n3+1)*(2*n1+31)
  call comb_rot_grow_loc_prev(0,n2,nt,y,ww,1,ibzxx) 
  ! dimension of ww: nt*(2*n2+31)=(n3+1)*(2*n1+31)*(2*n2+31)
  
  ! i3,I1,I2  -> I1,I2,I3
  nt=(2*n1+31)*(2*n2+31)
  call comb_rot_grow_loc_prev(0,n3,nt,ww,y,1,ibxxyy)
  
END SUBROUTINE comb_grow_c_prev


!> subroutine for scfunctions
!! In one dimesnion, 
!! with optimised cycles
!! Applies synthesis wavelet transformation 
!! then convolves with magic filter
!!  the size of the data is allowed to grow
subroutine comb_rot_grow_loc_prev(nfl,nfu,ndat,x,y,icf,ib)

  implicit real(kind=8) (a-h,o-z)
  integer t
  integer,parameter:: lowfil=-7,lupfil=8
  integer,parameter:: lowfil2=2*lowfil,lupfil2=2*lupfil
  dimension x(nfl:nfu,ndat),y(ndat,-14+2*nfl:2*nfu+16)
  integer ib(2,ndat)
  
  include 'v_long.inc'
  
  y=0.d0
  !open(unit=10,file='long.flop')
  
  !nflop=0
  !do l=1,ndat 
  !   if (ib(2,l).ge.ib(1,l)) nflop=nflop+(ib(2,l)-ib(1,l)+1)*31*2
  !enddo
  
  !call system_clock(ncount0,ncount_rate,ncount_max)
  
  do l=1,ndat
     
     !       loop for smaller i
     !       j=8 won't work : i-j would turn negative
     if (ib(1,l).le.ib(2,l)) then
        
        if (ib(1,l)+7<ib(2,l)-7) then
           do i=ib(1,l)-7,ib(1,l)+7
              y2i=0.d0
              y2i1=0.d0
              !     i+7=<ib(1,l)+14<ib(2,l)
              !     do t=        ib(1,l) ,min(i+7,ib(2,l))
              do t=ib(1,l),i+7
                 y2i =y2i +fil2(2*(i-t)  ,icf)*x(t,l)
                 y2i1=y2i1+fil2(2*(i-t)+1,icf)*x(t,l)
              enddo
              
              y(l,2*i)=y2i
              y(l,2*i+1)=y2i1
           enddo
           !***********************************************************************************************
           
           !   loop for ordinary i  
           
           if (ib(2,l)-ib(1,l)-16.ge.3) then
              !for the asymmetry of the daubechies filters
              !this loop can be unrolled with 3 units at most
              !otherwise the filter will go out of bounds
              do i=ib(1,l)+8,ib(2,l)-8-3,3
                 
                 y0=0.d0
                 y1=0.d0
                 y2=0.d0
                 y3=0.d0
                 y4=0.d0
                 y5=0.d0
!!!                 y6=0.d0
!!!                 y7=0.d0

                 do t=i-8,i+7+3 !this 3 must be discussed(does it depend on the loop unrolling?)
                    y0=y0+fil2(2*(i-t)+0,icf)*x(t,l)
                    y1=y1+fil2(2*(i-t)+1,icf)*x(t,l)
                    y2=y2+fil2(2*(i-t)+2,icf)*x(t,l)
                    y3=y3+fil2(2*(i-t)+3,icf)*x(t,l)
                    y4=y4+fil2(2*(i-t)+4,icf)*x(t,l)
                    y5=y5+fil2(2*(i-t)+5,icf)*x(t,l)
!!!                    y6=y6+fil2(2*(i-t)+6,icf)*x(t,l)
!!!                    y7=y7+fil2(2*(i-t)+7,icf)*x(t,l)
                 enddo

                 y(l,2*i+0)=y0
                 y(l,2*i+1)=y1
                 y(l,2*i+2)=y2
                 y(l,2*i+3)=y3
                 y(l,2*i+4)=y4
                 y(l,2*i+5)=y5
!!!                 y(l,2*i+6)=y6
!!!                 y(l,2*i+7)=y7
              enddo
              icur=i
           else
              icur=ib(1,l)+8
           endif

           !   loop for the "rigthmost of ordinary" i 
           do i=icur,ib(2,l)-8

              ! j=8  works since i is big enough
              y2i =fil2(16  ,icf)*x(i-8,l)
              y2i1=0.d0

              do t=i-7,i+7 
                 y2i =y2i +fil2(2*(i-t)  ,icf)*x(t,l)
                 y2i1=y2i1+fil2(2*(i-t)+1,icf)*x(t,l)
              enddo

              y(l,2*i)=y2i
              y(l,2*i+1)=y2i1
           enddo

           !*****************************************************************************
           !loop for the rightmost i
           do i=ib(2,l)-7,ib(2,l)+7

              y2i =fil2(16  ,icf)*x(i-8,l)
              y2i1=0.d0

              !  i+7=<ib(2,l)
              !     do t=    i-7         ,min(i+7,ib(2,l)) ! i-7>ib(1,l) since i=ib(1,l)+8,..
              do t=    i-7         ,    ib(2,l) 
                 y2i =y2i +fil2(2*(i-t)  ,icf)*x(t,l)
                 y2i1=y2i1+fil2(2*(i-t)+1,icf)*x(t,l)
              enddo

              y(l,2*i)=y2i
              y(l,2*i+1)=y2i1

           enddo

           t=ib(2,l); i=t+8;  ! to the rightmost element of y, only j=8 contributes
           y(l,2*i)=fil2(16,icf)*x(t,l)

        else!      ib(1,l)+7>=ib(2,l)-7!*****************************************************
           do i=ib(1,l)-7,ib(1,l)+7
              y2i=0.d0
              y2i1=0.d0

              ! i-7 =< ib(1,l)
              do t=        ib(1,l) ,min(i+7,ib(2,l))
                 y2i =y2i +fil2(2*(i-t)  ,icf)*x(t,l)
                 y2i1=y2i1+fil2(2*(i-t)+1,icf)*x(t,l)
              enddo

              y(l,2*i)=y2i
              y(l,2*i+1)=y2i1
           enddo
           !loop for ordinary i  
           do i=ib(1,l)+8,ib(2,l)+7

              ! j=8  works since i is big enough
              y2i =fil2(16  ,icf)*x(i-8,l)
              y2i1=0.d0

              !     do t=    i-7         ,min(i+7,ib(2,l)) 
              do t=i-7,ib(2,l) ! i-7>ib(1,l) since i=ib(1,l)+8,..
                 y2i =y2i +fil2(2*(i-t)  ,icf)*x(t,l)
                 y2i1=y2i1+fil2(2*(i-t)+1,icf)*x(t,l)
              enddo

              y(l,2*i)=y2i
              y(l,2*i+1)=y2i1
           enddo

           t=ib(2,l); i=t+8;  ! to the rightmost element of y, only j=8 contributes
           y(l,2*i)=fil2(16,icf)*x(t,l)

        endif

     endif
  enddo

  !call system_clock(ncount1,ncount_rate,ncount_max)
  !tel=dble(ncount1-ncount0)/dble(ncount_rate)

  !write(10,*) tel, 1.d-6*nflop/tel
END SUBROUTINE comb_rot_grow_loc_prev


!>   In one dimension, 
!!   with optimised cycles
!!   Applies synthesis wavelet transformation 
!!   then convolves with magic filter
!!   the size of the data is allowed to grow
subroutine  comb_rot_grow_loc_1_prev(nfl,nfu,ndat,x,y,ib)

  implicit none
!Arguments
  integer :: ndat,nfl,nfu
  real(kind=8), dimension(7,nfl:nfu,ndat) :: x
  real(kind=8), dimension(2,2,ndat,-14+2*nfl:2*nfu+16) :: y
  integer :: ib(2,ndat)
!Local variables
  integer,parameter:: lowfil=-7,lupfil=8
  integer,parameter:: lowfil2=2*lowfil,lupfil2=2*lupfil
  integer :: i,l
  real(kind=8) :: y1w11,y1w12,y1w21,y1w22,y2w11,y2w12,y2w21,y2w22
  integer :: t

  include 'v.inc'

  y=0.d0
  !open(unit=20,file='tree_unrolled.flop')

  !nflop=0
  !do l=1,ndat
  !   if (ib(2,l).ge.ib(1,l)) nflop=nflop+(ib(2,l)-ib(1,l)+1)*31*2*7
  !enddo

  !call system_clock(ncount0,ncount_rate,ncount_max)

  do l=1,ndat

     !       loop for smaller i
     !       j=8 won't work : i-j would turn negative
     if (ib(1,l).le.ib(2,l)) then
        do i=ib(1,l)-7,ib(1,l)+7

           y1w11=0.d0
           y1w21=0.d0
           y1w12=0.d0
           y1w22=0.d0

           y2w11=0.d0
           y2w21=0.d0
           y2w12=0.d0
           y2w22=0.d0

           ! i-7 =< ib(1,l)
           do t=        ib(1,l) ,min(i+7,ib(2,l))
              y1w11=y1w11+                           fil2(2*(i-t)  ,2)*x(1,t,l)
              y1w21=y1w21+fil2(2*(i-t)  ,1)*x(2,t,l)+fil2(2*(i-t)  ,2)*x(3,t,l)
              y1w12=y1w12+fil2(2*(i-t)  ,1)*x(4,t,l)+fil2(2*(i-t)  ,2)*x(5,t,l)
              y1w22=y1w22+fil2(2*(i-t)  ,1)*x(6,t,l)+fil2(2*(i-t)  ,2)*x(7,t,l)

              y2w11=y2w11+                           fil2(2*(i-t)+1,2)*x(1,t,l)
              y2w21=y2w21+fil2(2*(i-t)+1,1)*x(2,t,l)+fil2(2*(i-t)+1,2)*x(3,t,l)
              y2w12=y2w12+fil2(2*(i-t)+1,1)*x(4,t,l)+fil2(2*(i-t)+1,2)*x(5,t,l)
              y2w22=y2w22+fil2(2*(i-t)+1,1)*x(6,t,l)+fil2(2*(i-t)+1,2)*x(7,t,l)
           enddo

           y(1,1,l,2*i  )=y1w11
           y(2,1,l,2*i  )=y1w21
           y(1,2,l,2*i  )=y1w12
           y(2,2,l,2*i  )=y1w22

           y(1,1,l,2*i+1)=y2w11
           y(2,1,l,2*i+1)=y2w21
           y(1,2,l,2*i+1)=y2w12
           y(2,2,l,2*i+1)=y2w22
        enddo

        !   loop for ordinary i  
        do i=ib(1,l)+8,ib(2,l)+7

           ! j=8  works since i is big enough

           y1w11=fil2(16  ,2)*x(1,i-8,l)
           y1w21=fil2(16  ,1)*x(2,i-8,l)+fil2(16  ,2)*x(3,i-8,l)
           y1w12=fil2(16  ,1)*x(4,i-8,l)+fil2(16  ,2)*x(5,i-8,l)
           y1w22=fil2(16  ,1)*x(6,i-8,l)+fil2(16  ,2)*x(7,i-8,l)

           y2w11=0.d0
           y2w21=0.d0
           y2w12=0.d0
           y2w22=0.d0

           do t=    i-7         ,min(i+7,ib(2,l)) ! i-7>ib(1,l) since i=ib(1,l)+8,..
              y1w11=y1w11+                           fil2(2*(i-t)  ,2)*x(1,t,l)
              y1w21=y1w21+fil2(2*(i-t)  ,1)*x(2,t,l)+fil2(2*(i-t)  ,2)*x(3,t,l)
              y1w12=y1w12+fil2(2*(i-t)  ,1)*x(4,t,l)+fil2(2*(i-t)  ,2)*x(5,t,l)
              y1w22=y1w22+fil2(2*(i-t)  ,1)*x(6,t,l)+fil2(2*(i-t)  ,2)*x(7,t,l)

              y2w11=y2w11+                           fil2(2*(i-t)+1,2)*x(1,t,l)
              y2w21=y2w21+fil2(2*(i-t)+1,1)*x(2,t,l)+fil2(2*(i-t)+1,2)*x(3,t,l)
              y2w12=y2w12+fil2(2*(i-t)+1,1)*x(4,t,l)+fil2(2*(i-t)+1,2)*x(5,t,l)
              y2w22=y2w22+fil2(2*(i-t)+1,1)*x(6,t,l)+fil2(2*(i-t)+1,2)*x(7,t,l)
           enddo

           y(1,1,l,2*i  )=y1w11
           y(2,1,l,2*i  )=y1w21
           y(1,2,l,2*i  )=y1w12
           y(2,2,l,2*i  )=y1w22

           y(1,1,l,2*i+1)=y2w11
           y(2,1,l,2*i+1)=y2w21
           y(1,2,l,2*i+1)=y2w12
           y(2,2,l,2*i+1)=y2w22
        enddo

        t=ib(2,l);  i=t+8; ! to the rightmost element of y, only j=8 contributes

        y(1,1,l,2*i)=                    fil2(16,2)*x(1,t,l)
        y(2,1,l,2*i)=fil2(16,1)*x(2,t,l)+fil2(16,2)*x(3,t,l)
        y(1,2,l,2*i)=fil2(16,1)*x(4,t,l)+fil2(16,2)*x(5,t,l)
        y(2,2,l,2*i)=fil2(16,1)*x(6,t,l)+fil2(16,2)*x(7,t,l)

     endif
  enddo

  !call system_clock(ncount1,ncount_rate,ncount_max)
  !tel=dble(ncount1-ncount0)/dble(ncount_rate)

  !write(20,*) tel, 1.d-6*nflop/tel
END SUBROUTINE comb_rot_grow_loc_1_prev


!> In one dimesnion, 
!! with optimised cycles
!! Applies synthesis wavelet transformation 
!! then convolves with magic filter
!! the size of the data is allowed to grow
subroutine comb_rot_grow_loc_2_prev(nfl,nfu,ndat,x,y,ib)

  implicit real(kind=8) (a-h,o-z)
  integer t
  integer,parameter:: lowfil=-7,lupfil=8
  integer,parameter:: lowfil2=2*lowfil,lupfil2=2*lupfil
  dimension x(2,2,nfl:nfu,ndat),y(2,ndat,-14+2*nfl:2*nfu+16)
  integer ib(2,ndat)

  include 'v.inc'

  y=0.d0
  !open(unit=20,file='tree_unrolled.flop')

  nflop=0
  do l=1,ndat
     if (ib(2,l).ge.ib(1,l)) nflop=nflop+(ib(2,l)-ib(1,l)+1)*31*2*4
  enddo

  call system_clock(ncount0,ncount_rate,ncount_max)

  do l=1,ndat

     !       loop for smaller i
     !       j=8 won't work : i-j would turn negative
     if (ib(1,l).le.ib(2,l)) then
        do i=ib(1,l)-7,ib(1,l)+7

           y1w1=0.d0
           y1w2=0.d0

           y2w1=0.d0
           y2w2=0.d0

           ! i-7 =< ib(1,l)
           do t=        ib(1,l) ,min(i+7,ib(2,l))
              y1w1=y1w1+fil2(2*(i-t)  ,1)*x(1,1,t,l)+fil2(2*(i-t)  ,2)*x(2,1,t,l)
              y1w2=y1w2+fil2(2*(i-t)  ,1)*x(1,2,t,l)+fil2(2*(i-t)  ,2)*x(2,2,t,l)

              y2w1=y2w1+fil2(2*(i-t)+1,1)*x(1,1,t,l)+fil2(2*(i-t)+1,2)*x(2,1,t,l)
              y2w2=y2w2+fil2(2*(i-t)+1,1)*x(1,2,t,l)+fil2(2*(i-t)+1,2)*x(2,2,t,l)
           enddo

           y(1,l,2*i  )=y1w1
           y(2,l,2*i  )=y1w2

           y(1,l,2*i+1)=y2w1
           y(2,l,2*i+1)=y2w2
        enddo

        !   loop for ordinary i  
        do i=ib(1,l)+8,ib(2,l)+7

           ! j=8  works since i is big enough
           y1w1=fil2(16  ,1)*x(1,1,i-8,l)+fil2(16  ,2)*x(2,1,i-8,l)
           y1w2=fil2(16  ,1)*x(1,2,i-8,l)+fil2(16  ,2)*x(2,2,i-8,l)

           y2w1=0.d0
           y2w2=0.d0

           do t=    i-7         ,min(i+7,ib(2,l)) ! i-7>ib(1,l) since i=ib(1,l)+8,..
              y1w1=y1w1+fil2(2*(i-t)  ,1)*x(1,1,t,l)+fil2(2*(i-t)  ,2)*x(2,1,t,l)
              y1w2=y1w2+fil2(2*(i-t)  ,1)*x(1,2,t,l)+fil2(2*(i-t)  ,2)*x(2,2,t,l)

              y2w1=y2w1+fil2(2*(i-t)+1,1)*x(1,1,t,l)+fil2(2*(i-t)+1,2)*x(2,1,t,l)
              y2w2=y2w2+fil2(2*(i-t)+1,1)*x(1,2,t,l)+fil2(2*(i-t)+1,2)*x(2,2,t,l)
           enddo

           y(1,l,2*i  )=y1w1
           y(2,l,2*i  )=y1w2

           y(1,l,2*i+1)=y2w1
           y(2,l,2*i+1)=y2w2
        enddo

        t=ib(2,l); i=t+8;  ! to the rightmost element of y, only j=8 contributes
        y(1,l,2*i)=fil2(16,1)*x(1,1,t,l)+fil2(16,2)*x(2,1,t,l)
        y(2,l,2*i)=fil2(16,1)*x(1,2,t,l)+fil2(16,2)*x(2,2,t,l)

     endif
  enddo

  call system_clock(ncount1,ncount_rate,ncount_max)
  tel=dble(ncount1-ncount0)/dble(ncount_rate)

  !write(20,*) tel, 1.d-6*nflop/tel
END SUBROUTINE comb_rot_grow_loc_2_prev


!> In one dimension, 
!! with optimised cycles
!! Applies synthesis wavelet transformation 
!! then convolves with magic filter
!! then adds the result to y.
!! The size of the data is allowed to grow
subroutine comb_rot_grow_loc_3_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y,ib)
    implicit real(kind=8) (a-h,o-z)
    integer t
    integer,parameter:: lowfil=-7,lupfil=8
    integer,parameter:: lowfil2=2*lowfil,lupfil2=2*lupfil

    dimension  x(2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16)
    dimension  y(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)
    integer ib(2,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16)

    include 'v_long.inc'

    !open(unit=20,file='tree_unrolled.flop')

    !nflop=0

    !do l1=-14+2*nfl1,2*nfu1+16
    !    do l2=-14+2*nfl2,2*nfu2+16
    !        if (ib(2,l1,l2).ge.ib(1,l1,l2)) nflop=nflop+(ib(2,l1,l2)-ib(1,l1,l2)+1)*31*2*2
!       additional 2 because of wavelet/scfunction pair
    !    enddo
    !enddo

    !call system_clock(ncount0,ncount_rate,ncount_max)

do l1=-14+2*nfl1,2*nfu1+16
    do l2=-14+2*nfl2,2*nfu2+16
!       loop for smaller i
!       j=8 won't work : i-j would turn negative
        if (ib(1,l1,l2).le.ib(2,l1,l2)) then
            if (ib(1,l1,l2)+7<ib(2,l1,l2)-7) then
                ! the case of a long segment
                do i=ib(1,l1,l2)-7,ib(1,l1,l2)+7
                    y2i=0.d0
                    y2i1=0.d0

                    ! i-7 =< ib(1,l1,l2)
                    do t=         ib(1,l1,l2) ,min(i+7,ib(2,l1,l2))
                        y2i =y2i +fil2(2*(i-t)  ,1)*x(1,t,l1,l2)+fil2(2*(i-t)  ,2)*x(2,t,l1,l2)
                        y2i1=y2i1+fil2(2*(i-t)+1,1)*x(1,t,l1,l2)+fil2(2*(i-t)+1,2)*x(2,t,l1,l2)
                    enddo
        
                    y(l1,l2,2*i  )=y(l1,l2,2*i  )+y2i
                    y(l1,l2,2*i+1)=y(l1,l2,2*i+1)+y2i1
                enddo
    
    !            loop for ordinary i        
                if (ib(2,l1,l2)-ib(1,l1,l2)-16.ge.4) then
                    do i=ib(1,l1,l2)+8,ib(2,l1,l2)-8-4,4
        
                        y0=0.d0
                        y1=0.d0
                        y2=0.d0
                        y3=0.d0
                        y4=0.d0
                        y5=0.d0
                        y6=0.d0
                        y7=0.d0
            
                        ! i-7>ib(1,l) since i=ib(1,l)+8,..
                        do t=    i-8             ,i+7+3
                            y0=y0+fil2(2*(i-t)+0,1)*x(1,t,l1,l2)+fil2(2*(i-t)+0,2)*x(2,t,l1,l2)
                            y1=y1+fil2(2*(i-t)+1,1)*x(1,t,l1,l2)+fil2(2*(i-t)+1,2)*x(2,t,l1,l2)
                            y2=y2+fil2(2*(i-t)+2,1)*x(1,t,l1,l2)+fil2(2*(i-t)+2,2)*x(2,t,l1,l2)
                            y3=y3+fil2(2*(i-t)+3,1)*x(1,t,l1,l2)+fil2(2*(i-t)+3,2)*x(2,t,l1,l2)
                            y4=y4+fil2(2*(i-t)+4,1)*x(1,t,l1,l2)+fil2(2*(i-t)+4,2)*x(2,t,l1,l2)
                            y5=y5+fil2(2*(i-t)+5,1)*x(1,t,l1,l2)+fil2(2*(i-t)+5,2)*x(2,t,l1,l2)
                            y6=y6+fil2(2*(i-t)+6,1)*x(1,t,l1,l2)+fil2(2*(i-t)+6,2)*x(2,t,l1,l2)
                            y7=y7+fil2(2*(i-t)+7,1)*x(1,t,l1,l2)+fil2(2*(i-t)+7,2)*x(2,t,l1,l2)
                        enddo
            
                        y(l1,l2,2*i+0)=y(l1,l2,2*i+0)+y0
                        y(l1,l2,2*i+1)=y(l1,l2,2*i+1)+y1
                        y(l1,l2,2*i+2)=y(l1,l2,2*i+2)+y2
                        y(l1,l2,2*i+3)=y(l1,l2,2*i+3)+y3
                        y(l1,l2,2*i+4)=y(l1,l2,2*i+4)+y4
                        y(l1,l2,2*i+5)=y(l1,l2,2*i+5)+y5
                        y(l1,l2,2*i+6)=y(l1,l2,2*i+6)+y6
                        y(l1,l2,2*i+7)=y(l1,l2,2*i+7)+y7

                    enddo
            
                    icur=i
                else
                    icur=ib(1,l1,l2)+8
                endif

                do i=icur,ib(2,l1,l2)-8
    
                    !    j=8  works since i is big enough
                    y2i =fil2(16  ,1)*x(1,i-8,l1,l2)+fil2(16  ,2)*x(2,i-8,l1,l2)
                    y2i1=0.d0
        
                    ! i-7>ib(1,l) since i=ib(1,l)+8,..
                    do t=    i-7,i+7
                        y2i =y2i +fil2(2*(i-t)  ,1)*x(1,t,l1,l2)+fil2(2*(i-t)  ,2)*x(2,t,l1,l2)
                        y2i1=y2i1+fil2(2*(i-t)+1,1)*x(1,t,l1,l2)+fil2(2*(i-t)+1,2)*x(2,t,l1,l2)
                    enddo
        
                    y(l1,l2,2*i  )=y(l1,l2,2*i  )+y2i
                    y(l1,l2,2*i+1)=y(l1,l2,2*i+1)+y2i1
                enddo
        
                do i=ib(2,l1,l2)-7,ib(2,l1,l2)+7
    
                    y2i =fil2(16  ,1)*x(1,i-8,l1,l2)+fil2(16  ,2)*x(2,i-8,l1,l2)
                    y2i1=0.d0
        
                    do t=    i-7,ib(2,l1,l2)
                        y2i =y2i +fil2(2*(i-t)  ,1)*x(1,t,l1,l2)+fil2(2*(i-t)  ,2)*x(2,t,l1,l2)
                        y2i1=y2i1+fil2(2*(i-t)+1,1)*x(1,t,l1,l2)+fil2(2*(i-t)+1,2)*x(2,t,l1,l2)
                    enddo
        
                    y(l1,l2,2*i  )=y(l1,l2,2*i  )+y2i
                    y(l1,l2,2*i+1)=y(l1,l2,2*i+1)+y2i1
                enddo
                t=ib(2,l1,l2);    i=t+8;     ! to the rightmost element of y, only j=8 contributes
                y(l1,l2,2*i)=y(l1,l2,2*i)+fil2(16,1)*x(1,t,l1,l2)+fil2(16,2)*x(2,t,l1,l2)




            else
                ! the case of a short segment
                do i=ib(1,l1,l2)-7,ib(1,l1,l2)+7
                    y2i=0.d0
                    y2i1=0.d0

                    ! i-7 =< ib(1,l1,l2)
                    do t=         ib(1,l1,l2) ,min(i+7,ib(2,l1,l2))
                        y2i =y2i +fil2(2*(i-t)  ,1)*x(1,t,l1,l2)+fil2(2*(i-t)  ,2)*x(2,t,l1,l2)
                        y2i1=y2i1+fil2(2*(i-t)+1,1)*x(1,t,l1,l2)+fil2(2*(i-t)+1,2)*x(2,t,l1,l2)
                    enddo
        
                    y(l1,l2,2*i  )=y(l1,l2,2*i  )+y2i
                    y(l1,l2,2*i+1)=y(l1,l2,2*i+1)+y2i1
                enddo
    
    !            loop for ordinary i        
                do i=ib(1,l1,l2)+8,ib(2,l1,l2)+7
        
                    !    j=8  works since i is big enough
                    y2i =fil2(16  ,1)*x(1,i-8,l1,l2)+fil2(16  ,2)*x(2,i-8,l1,l2)
                    y2i1=0.d0
        
                    ! i-7>ib(1,l) since i=ib(1,l)+8,..
                    do t=    i-7             ,min(i+7,ib(2,l1,l2)) 
                        y2i =y2i +fil2(2*(i-t)  ,1)*x(1,t,l1,l2)+fil2(2*(i-t)  ,2)*x(2,t,l1,l2)
                        y2i1=y2i1+fil2(2*(i-t)+1,1)*x(1,t,l1,l2)+fil2(2*(i-t)+1,2)*x(2,t,l1,l2)
                    enddo
        
                    y(l1,l2,2*i  )=y(l1,l2,2*i  )+y2i
                    y(l1,l2,2*i+1)=y(l1,l2,2*i+1)+y2i1
                enddo
        
                t=ib(2,l1,l2);    i=t+8;     ! to the rightmost element of y, only j=8 contributes
                y(l1,l2,2*i)=y(l1,l2,2*i)+fil2(16,1)*x(1,t,l1,l2)+fil2(16,2)*x(2,t,l1,l2)
            endif

        endif
    enddo
enddo

    !call system_clock(ncount1,ncount_rate,ncount_max)
    !tel=dble(ncount1-ncount0)/dble(ncount_rate)

    !write(20,*) tel, 1.d-6*nflop/tel

END SUBROUTINE comb_rot_grow_loc_3_prev


!> y = (kinetic energy operator)x + (cprec*I)x 
subroutine Convolkinetic_prev(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_fc,x_f,y_c,y_f)
  implicit real(kind=8) (a-h,o-z)
  logical :: firstcall=.true. 
  integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
  dimension x_c(0:n1,0:n2,0:n3),y_c(0:n1,0:n2,0:n3)
  dimension x_fc(0:n1,0:n2,0:n3,3),x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
  dimension y_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
  dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)


  parameter(lowfil=-14,lupfil=14)
  dimension a(lowfil:lupfil),b(lowfil:lupfil),c(lowfil:lupfil),e(lowfil:lupfil)
  scale=-.5d0/hgrid**2
  !---------------------------------------------------------------------------
  ! second derivative filters for Daubechies 16
  !  <phi|D^2|phi_i>
  a(0)=   -3.5536922899131901941296809374d0*scale
  a(1)=    2.2191465938911163898794546405d0*scale
  a(2)=   -0.6156141465570069496314853949d0*scale
  a(3)=    0.2371780582153805636239247476d0*scale
  a(4)=   -0.0822663999742123340987663521d0*scale
  a(5)=    0.02207029188482255523789911295638968409d0*scale
  a(6)=   -0.409765689342633823899327051188315485d-2*scale
  a(7)=    0.45167920287502235349480037639758496d-3*scale
  a(8)=   -0.2398228524507599670405555359023135d-4*scale
  a(9)=    2.0904234952920365957922889447361d-6*scale
  a(10)=  -3.7230763047369275848791496973044d-7*scale
  a(11)=  -1.05857055496741470373494132287d-8*scale
  a(12)=  -5.813879830282540547959250667d-11*scale
  a(13)=   2.70800493626319438269856689037647576d-13*scale
  a(14)=  -6.924474940639200152025730585882d-18*scale
  do i=1,14
     a(-i)=a(i)
  enddo
  !  <phi|D^2|psi_i>
  c(-14)=     -3.869102413147656535541850057188d-18*scale
  c(-13)=      1.5130616560866154733900029272077362d-13*scale
  c(-12)=     -3.2264702314010525539061647271983988409d-11*scale
  c(-11)=     -5.96264938781402337319841002642d-9*scale
  c(-10)=     -2.1656830629214041470164889350342d-7*scale
  c(-9 )=      8.7969704055286288323596890609625d-7*scale
  c(-8 )=     -0.00001133456724516819987751818232711775d0*scale
  c(-7 )=      0.00021710795484646138591610188464622454d0*scale
  c(-6 )=     -0.0021356291838797986414312219042358542d0*scale
  c(-5 )=      0.00713761218453631422925717625758502986d0*scale
  c(-4 )=     -0.0284696165863973422636410524436931061d0*scale
  c(-3 )=      0.14327329352510759457155821037742893841d0*scale
  c(-2 )=     -0.42498050943780130143385739554118569733d0*scale
  c(-1 )=      0.65703074007121357894896358254040272157d0*scale
  c( 0 )=     -0.42081655293724308770919536332797729898d0*scale
  c( 1 )=     -0.21716117505137104371463587747283267899d0*scale
  c( 2 )=      0.63457035267892488185929915286969303251d0*scale
  c( 3 )=     -0.53298223962800395684936080758073568406d0*scale
  c( 4 )=      0.23370490631751294307619384973520033236d0*scale
  c( 5 )=     -0.05657736973328755112051544344507997075d0*scale
  c( 6 )=      0.0080872029411844780634067667008050127d0*scale
  c( 7 )=     -0.00093423623304808664741804536808932984d0*scale
  c( 8 )=      0.00005075807947289728306309081261461095d0*scale
  c( 9 )=     -4.62561497463184262755416490048242d-6*scale
  c( 10)=      6.3919128513793415587294752371778d-7*scale
  c( 11)=      1.87909235155149902916133888931d-8*scale
  c( 12)=      1.04757345962781829480207861447155543883d-10*scale
  c( 13)=     -4.84665690596158959648731537084025836d-13*scale
  c( 14)=      1.2392629629188986192855777620877d-17*scale
  !  <psi|D^2|phi_i>
  do i=-14,14
     b(i)=c(-i)
  enddo
  !<psi|D^2|psi_i>
  e(0)=   -24.875846029392331358907766562d0*scale
  e(1)=   -7.1440597663471719869313377994d0*scale
  e(2)=   -0.04251705323669172315864542163525830944d0*scale
  e(3)=   -0.26995931336279126953587091167128839196d0*scale
  e(4)=    0.08207454169225172612513390763444496516d0*scale
  e(5)=   -0.02207327034586634477996701627614752761d0*scale
  e(6)=    0.00409765642831595181639002667514310145d0*scale
  e(7)=   -0.00045167920287507774929432548999880117d0*scale
  e(8)=    0.00002398228524507599670405555359023135d0*scale
  e(9)=   -2.0904234952920365957922889447361d-6*scale
  e(10)=   3.7230763047369275848791496973044d-7*scale
  e(11)=   1.05857055496741470373494132287d-8*scale
  e(12)=   5.8138798302825405479592506674648873655d-11*scale
  e(13)=  -2.70800493626319438269856689037647576d-13*scale
  e(14)=   6.924474940639200152025730585882d-18*scale
  do i=1,14
     e(-i)=e(i)
  enddo


  if (firstcall) then

     ! (1/2) d^2/dx^2
     mflop1=0
     do i3=0,n3
        do i2=0,n2
           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
              do l=max(-i1,lowfil),min(lupfil,n1-i1)
                 mflop1=mflop1+4
              enddo
              mflop1=mflop1+3
           enddo

           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
              do l=max(-i1,lowfil),min(lupfil,n1-i1)
                 mflop1=mflop1+2
              enddo
           enddo
        enddo
     enddo

     ! + (1/2) d^2/dy^2
     mflop2=0
     do i3=0,n3
        do i1=0,n1
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
              do l=max(-i2,lowfil),min(lupfil,n2-i2)
                 mflop2=mflop2+4
              enddo
              mflop2=mflop2+2
           enddo

           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
              do l=max(-i2,lowfil),min(lupfil,n2-i2)
                 mflop2=mflop2+2
              enddo
           enddo
        enddo
     enddo


     ! + (1/2) d^2/dz^2
     mflop3=0
     do i2=0,n2
        do i1=0,n1
           do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
              do l=max(-i3,lowfil),min(lupfil,n3-i3)
                 mflop3=mflop3+4
              enddo
              mflop3=mflop3+2
           enddo

           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
              do l=max(-i3,lowfil),min(lupfil,n3-i3)
                 mflop3=mflop3+2
              enddo
           enddo
        enddo
     enddo

     ! wavelet part
     ! (1/2) d^2/dx^2
     nflop1=0
     do i3=nfl3,nfu3
        do i2=nfl2,nfu2
           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
              do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
                 nflop1=nflop1+26
              enddo
              nflop1=nflop1+17
           enddo
        enddo
     enddo

     ! + (1/2) d^2/dy^2
     nflop2=0
     do i3=nfl3,nfu3
        do i1=nfl1,nfu1
           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
              do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
                 nflop2=nflop2+26
              enddo
              nflop2=nflop2+7
           enddo
        enddo
     enddo

     ! + (1/2) d^2/dz^2
     nflop3=0
     do i2=nfl2,nfu2
        do i1=nfl1,nfu1
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
              do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
                 nflop3=nflop3+26
              enddo
              nflop3=nflop3+7
           enddo
        enddo
     enddo

     firstcall=.false.
  endif


  !---------------------------------------------------------------------------

  ! Scaling function part

  call system_clock(ncount0,ncount_rate,ncount_max)

  ! (1/2) d^2/dx^2
  do i3=0,n3
     do i2=0,n2
        do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
           t111=0.d0 ; s111=0.d0
           do l=max(-i1,lowfil),min(lupfil,n1-i1)
              t111=t111 + x_c(i1+l,i2,i3)*a(l)
              s111=s111 + x_fc(i1+l,i2,i3,1)*b(l)
           enddo
           y_c(i1,i2,i3)=t111+s111+cprecr*x_c(i1,i2,i3)
        enddo

        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t211=0.d0 
           do l=max(-i1,lowfil),min(lupfil,n1-i1)
              t211=t211 + x_c(i1+l,i2,i3)*c(l)
           enddo
           y_f(1,i1,i2,i3)=t211
        enddo
     enddo
  enddo

  call system_clock(ncount1,ncount_rate,ncount_max)
  tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'FIRST PART:x',tel,1.d-6*mflop1/tel

  ! + (1/2) d^2/dy^2
  do i3=0,n3
     do i1=0,n1
        do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
           t111=0.d0 ; s111=0.d0
           do l=max(-i2,lowfil),min(lupfil,n2-i2)
              t111=t111 + x_c(i1,i2+l,i3)*a(l)
              s111=s111 + x_fc(i1,i2+l,i3,2)*b(l)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+t111+s111
        enddo

        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t121=0.d0 
           do l=max(-i2,lowfil),min(lupfil,n2-i2)
              t121=t121 + x_c(i1,i2+l,i3)*c(l)
           enddo
           y_f(2,i1,i2,i3)=t121
        enddo
     enddo
  enddo


  call system_clock(ncount2,ncount_rate,ncount_max)
  tel=dble(ncount2-ncount1)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'FIRST PART:y',tel,1.d-6*mflop2/tel

  ! + (1/2) d^2/dz^2
  do i2=0,n2
     do i1=0,n1
        do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
           t111=0.d0 ; s111=0.d0
           do l=max(-i3,lowfil),min(lupfil,n3-i3)
              t111=t111 + x_c(i1,i2,i3+l)*a(l)
              s111=s111 + x_fc(i1,i2,i3+l,3)*b(l)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+t111+s111
        enddo

        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0.d0 
           do l=max(-i3,lowfil),min(lupfil,n3-i3)
              t112=t112 + x_c(i1,i2,i3+l)*c(l)
           enddo
           y_f(4,i1,i2,i3)=t112
        enddo
     enddo
  enddo

  call system_clock(ncount3,ncount_rate,ncount_max)
  tel=dble(ncount3-ncount2)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'FIRST PART:z',tel,1.d-6*mflop3/tel

  ! wavelet part
  ! (1/2) d^2/dx^2
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              t112=t112 + x_f(4,i1+l,i2,i3)*a(l) + x_f(5,i1+l,i2,i3)*b(l)
              t121=t121 + x_f(2,i1+l,i2,i3)*a(l) + x_f(3,i1+l,i2,i3)*b(l)
              t122=t122 + x_f(6,i1+l,i2,i3)*a(l) + x_f(7,i1+l,i2,i3)*b(l)
              t212=t212 + x_f(4,i1+l,i2,i3)*c(l) + x_f(5,i1+l,i2,i3)*e(l)
              t221=t221 + x_f(2,i1+l,i2,i3)*c(l) + x_f(3,i1+l,i2,i3)*e(l)
              t222=t222 + x_f(6,i1+l,i2,i3)*c(l) + x_f(7,i1+l,i2,i3)*e(l)
              t211=t211 + x_f(1,i1+l,i2,i3)*e(l)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112+cprecr*x_f(4,i1,i2,i3)
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121+cprecr*x_f(2,i1,i2,i3)
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211+cprecr*x_f(1,i1,i2,i3)
           y_f(6,i1,i2,i3)=t122+cprecr*x_f(6,i1,i2,i3)
           y_f(5,i1,i2,i3)=t212+cprecr*x_f(5,i1,i2,i3)
           y_f(3,i1,i2,i3)=t221+cprecr*x_f(3,i1,i2,i3)
           y_f(7,i1,i2,i3)=t222+cprecr*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo

  call system_clock(ncount4,ncount_rate,ncount_max)
  tel=dble(ncount4-ncount3)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'SECND PART:x',tel,1.d-6*nflop1/tel


  ! + (1/2) d^2/dy^2
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              t112=t112 + x_f(4,i1,i2+l,i3)*a(l) + x_f(6,i1,i2+l,i3)*b(l)
              t211=t211 + x_f(1,i1,i2+l,i3)*a(l) + x_f(3,i1,i2+l,i3)*b(l)
              t122=t122 + x_f(4,i1,i2+l,i3)*c(l) + x_f(6,i1,i2+l,i3)*e(l)
              t212=t212 + x_f(5,i1,i2+l,i3)*a(l) + x_f(7,i1,i2+l,i3)*b(l)
              t221=t221 + x_f(1,i1,i2+l,i3)*c(l) + x_f(3,i1,i2+l,i3)*e(l)
              t222=t222 + x_f(5,i1,i2+l,i3)*c(l) + x_f(7,i1,i2+l,i3)*e(l)
              t121=t121 + x_f(2,i1,i2+l,i3)*e(l)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
        enddo
     enddo
  enddo

  call system_clock(ncount5,ncount_rate,ncount_max)
  tel=dble(ncount5-ncount4)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'SECND PART:y',tel,1.d-6*nflop2/tel

  ! + (1/2) d^2/dz^2
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
              t121=t121 + x_f(2,i1,i2,i3+l)*a(l) + x_f(6,i1,i2,i3+l)*b(l)
              t211=t211 + x_f(1,i1,i2,i3+l)*a(l) + x_f(5,i1,i2,i3+l)*b(l)
              t122=t122 + x_f(2,i1,i2,i3+l)*c(l) + x_f(6,i1,i2,i3+l)*e(l)
              t212=t212 + x_f(1,i1,i2,i3+l)*c(l) + x_f(5,i1,i2,i3+l)*e(l)
              t221=t221 + x_f(3,i1,i2,i3+l)*a(l) + x_f(7,i1,i2,i3+l)*b(l)
              t222=t222 + x_f(3,i1,i2,i3+l)*c(l) + x_f(7,i1,i2,i3+l)*e(l)
              t112=t112 + x_f(4,i1,i2,i3+l)*e(l)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222

        enddo
     enddo
  enddo

  call system_clock(ncount6,ncount_rate,ncount_max)
  tel=dble(ncount6-ncount5)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'SECND PART:z',tel,1.d-6*nflop3/tel

  tel=dble(ncount6-ncount0)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'ALL   PART',  & 
  !            tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel

  return
END SUBROUTINE Convolkinetic_prev


!>  y = y+(kinetic energy operator)x 
subroutine ConvolkineticT_prev(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_fc,x_f,y_c,y_f,ekin)
  implicit real(kind=8) (a-h,o-z)
  !Arguments
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
  real(kind=8), intent(in) :: hgrid
  real(kind=8), intent(out) :: ekin
  real(kind=8) :: x_c(0:n1,0:n2,0:n3),y_c(0:n1,0:n2,0:n3)
  real(kind=8) :: x_fc(0:n1,0:n2,0:n3,3),x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
  real(kind=8) :: y_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
  integer ::  ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  integer ::  ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
  !Local variables
  logical :: firstcall=.true. 
  integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
  integer, parameter :: lowfil=-14,lupfil=14
  real(kind=8) :: a(lowfil:lupfil),b(lowfil:lupfil),c(lowfil:lupfil),e(lowfil:lupfil)
  integer :: i,i1,i2,i3,l,nb,ncount0
  integer :: ncount1,ncount2,ncount3,ncount4,ncount5,ncount6,ncount_max,ncount_rate
  real(kind=8) :: s111,t111,t112,t121,t122,t211,t212,t221,t222,tel
  scale=-.5d0/hgrid**2
  !---------------------------------------------------------------------------
  ! second derivative filters for Daubechies 16
  !  <phi|D^2|phi_i>
  a(0)=   -3.5536922899131901941296809374d0*scale
  a(1)=    2.2191465938911163898794546405d0*scale
  a(2)=   -0.6156141465570069496314853949d0*scale
  a(3)=    0.2371780582153805636239247476d0*scale
  a(4)=   -0.0822663999742123340987663521d0*scale
  a(5)=    0.02207029188482255523789911295638968409d0*scale
  a(6)=   -0.409765689342633823899327051188315485d-2*scale
  a(7)=    0.45167920287502235349480037639758496d-3*scale
  a(8)=   -0.2398228524507599670405555359023135d-4*scale
  a(9)=    2.0904234952920365957922889447361d-6*scale
  a(10)=  -3.7230763047369275848791496973044d-7*scale
  a(11)=  -1.05857055496741470373494132287d-8*scale
  a(12)=  -5.813879830282540547959250667d-11*scale
  a(13)=   2.70800493626319438269856689037647576d-13*scale
  a(14)=  -6.924474940639200152025730585882d-18*scale
  do i=1,14
     a(-i)=a(i)
  enddo
  !  <phi|D^2|psi_i>
  c(-14)=     -3.869102413147656535541850057188d-18*scale
  c(-13)=      1.5130616560866154733900029272077362d-13*scale
  c(-12)=     -3.2264702314010525539061647271983988409d-11*scale
  c(-11)=     -5.96264938781402337319841002642d-9*scale
  c(-10)=     -2.1656830629214041470164889350342d-7*scale
  c(-9 )=      8.7969704055286288323596890609625d-7*scale
  c(-8 )=     -0.00001133456724516819987751818232711775d0*scale
  c(-7 )=      0.00021710795484646138591610188464622454d0*scale
  c(-6 )=     -0.0021356291838797986414312219042358542d0*scale
  c(-5 )=      0.00713761218453631422925717625758502986d0*scale
  c(-4 )=     -0.0284696165863973422636410524436931061d0*scale
  c(-3 )=      0.14327329352510759457155821037742893841d0*scale
  c(-2 )=     -0.42498050943780130143385739554118569733d0*scale
  c(-1 )=      0.65703074007121357894896358254040272157d0*scale
  c( 0 )=     -0.42081655293724308770919536332797729898d0*scale
  c( 1 )=     -0.21716117505137104371463587747283267899d0*scale
  c( 2 )=      0.63457035267892488185929915286969303251d0*scale
  c( 3 )=     -0.53298223962800395684936080758073568406d0*scale
  c( 4 )=      0.23370490631751294307619384973520033236d0*scale
  c( 5 )=     -0.05657736973328755112051544344507997075d0*scale
  c( 6 )=      0.0080872029411844780634067667008050127d0*scale
  c( 7 )=     -0.00093423623304808664741804536808932984d0*scale
  c( 8 )=      0.00005075807947289728306309081261461095d0*scale
  c( 9 )=     -4.62561497463184262755416490048242d-6*scale
  c( 10)=      6.3919128513793415587294752371778d-7*scale
  c( 11)=      1.87909235155149902916133888931d-8*scale
  c( 12)=      1.04757345962781829480207861447155543883d-10*scale
  c( 13)=     -4.84665690596158959648731537084025836d-13*scale
  c( 14)=      1.2392629629188986192855777620877d-17*scale
  !  <psi|D^2|phi_i>
  do i=-14,14
     b(i)=c(-i)
  enddo
  !<psi|D^2|psi_i>
  e(0)=   -24.875846029392331358907766562d0*scale
  e(1)=   -7.1440597663471719869313377994d0*scale
  e(2)=   -0.04251705323669172315864542163525830944d0*scale
  e(3)=   -0.26995931336279126953587091167128839196d0*scale
  e(4)=    0.08207454169225172612513390763444496516d0*scale
  e(5)=   -0.02207327034586634477996701627614752761d0*scale
  e(6)=    0.00409765642831595181639002667514310145d0*scale
  e(7)=   -0.00045167920287507774929432548999880117d0*scale
  e(8)=    0.00002398228524507599670405555359023135d0*scale
  e(9)=   -2.0904234952920365957922889447361d-6*scale
  e(10)=   3.7230763047369275848791496973044d-7*scale
  e(11)=   1.05857055496741470373494132287d-8*scale
  e(12)=   5.8138798302825405479592506674648873655d-11*scale
  e(13)=  -2.70800493626319438269856689037647576d-13*scale
  e(14)=   6.924474940639200152025730585882d-18*scale
  do i=1,14
     e(-i)=e(i)
  enddo

  if (firstcall) then

     ! (1/2) d^2/dx^2
     mflop1=0
     do i3=0,n3
        do i2=0,n2
           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
              do l=max(-i1,lowfil),min(lupfil,n1-i1)
                 mflop1=mflop1+4
              enddo
              mflop1=mflop1+4
           enddo

           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
              do l=max(-i1,lowfil),min(lupfil,n1-i1)
                 mflop1=mflop1+3
              enddo
           enddo
        enddo
     enddo

     ! + (1/2) d^2/dy^2
     mflop2=0
     do i3=0,n3
        do i1=0,n1
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
              do l=max(-i2,lowfil),min(lupfil,n2-i2)
                 mflop2=mflop2+4
              enddo
              mflop2=mflop2+4
           enddo

           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
              do l=max(-i2,lowfil),min(lupfil,n2-i2)
                 mflop2=mflop2+3
              enddo
           enddo
        enddo
     enddo


     ! + (1/2) d^2/dz^2
     mflop3=0
     do i2=0,n2
        do i1=0,n1
           do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
              do l=max(-i3,lowfil),min(lupfil,n3-i3)
                 mflop3=mflop3+4
              enddo
              mflop3=mflop3+4
           enddo

           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
              do l=max(-i3,lowfil),min(lupfil,n3-i3)
                 mflop3=mflop3+3
              enddo
           enddo
        enddo
     enddo

     ! wavelet part
     ! (1/2) d^2/dx^2
     nflop1=0
     do i3=nfl3,nfu3
        do i2=nfl2,nfu2
           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
              do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
                 nflop1=nflop1+26
              enddo
              nflop1=nflop1+21
           enddo
        enddo
     enddo

     ! + (1/2) d^2/dy^2
     nflop2=0
     do i3=nfl3,nfu3
        do i1=nfl1,nfu1
           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
              do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
                 nflop2=nflop2+26
              enddo
              nflop2=nflop2+21
           enddo
        enddo
     enddo

     ! + (1/2) d^2/dz^2
     nflop3=0
     do i2=nfl2,nfu2
        do i1=nfl1,nfu1
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
              do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
                 nflop3=nflop3+26
              enddo
              nflop3=nflop3+21
           enddo
        enddo
     enddo

     firstcall=.false.
  endif

  !---------------------------------------------------------------------------

  ! Scaling function part

  call system_clock(ncount0,ncount_rate,ncount_max)
  ekin=0.d0

  ! (1/2) d^2/dx^2
  do i3=0,n3
     do i2=0,n2
        do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
           t111=0.d0 ; s111=0.d0
           do l=max(-i1,lowfil),min(lupfil,n1-i1)
              t111=t111 + x_c(i1+l,i2,i3)*a(l)
              s111=s111 + x_fc(i1+l,i2,i3,1)*b(l)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+(t111+s111)
           ekin=ekin+(t111+s111)*x_c(i1,i2,i3)
        enddo

        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t211=0.d0 
           do l=max(-i1,lowfil),min(lupfil,n1-i1)
              t211=t211 + x_c(i1+l,i2,i3)*c(l)
           enddo
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           ekin=ekin+t211*x_f(1,i1,i2,i3)
        enddo
     enddo
  enddo

  call system_clock(ncount1,ncount_rate,ncount_max)
  tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'T:FIRST PART:x',tel,1.d-6*mflop1/tel

  ! + (1/2) d^2/dy^2
  do i3=0,n3
     do i1=0,n1
        do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
           t111=0.d0 ; s111=0.d0
           do l=max(-i2,lowfil),min(lupfil,n2-i2)
              t111=t111 + x_c(i1,i2+l,i3)*a(l)
              s111=s111 + x_fc(i1,i2+l,i3,2)*b(l)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+(t111+s111)
           ekin=ekin+(t111+s111)*x_c(i1,i2,i3)
        enddo

        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t121=0.d0 
           do l=max(-i2,lowfil),min(lupfil,n2-i2)
              t121=t121 + x_c(i1,i2+l,i3)*c(l)
           enddo
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           ekin=ekin+t121*x_f(2,i1,i2,i3)
        enddo
     enddo
  enddo


  call system_clock(ncount2,ncount_rate,ncount_max)
  tel=dble(ncount2-ncount1)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'T:FIRST PART:y',tel,1.d-6*mflop2/tel

  ! + (1/2) d^2/dz^2
  do i2=0,n2
     do i1=0,n1
        do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
           t111=0.d0 ; s111=0.d0
           do l=max(-i3,lowfil),min(lupfil,n3-i3)
              t111=t111 + x_c(i1,i2,i3+l)*a(l)
              s111=s111 + x_fc(i1,i2,i3+l,3)*b(l)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+(t111+s111)
           ekin=ekin+(t111+s111)*x_c(i1,i2,i3)
        enddo

        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0.d0 
           do l=max(-i3,lowfil),min(lupfil,n3-i3)
              t112=t112 + x_c(i1,i2,i3+l)*c(l)
           enddo
           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           ekin=ekin+t112*x_f(4,i1,i2,i3)
        enddo
     enddo
  enddo

  call system_clock(ncount3,ncount_rate,ncount_max)
  tel=dble(ncount3-ncount2)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'T:FIRST PART:z',tel,1.d-6*mflop3/tel

  ! wavelet part
  ! (1/2) d^2/dx^2
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              t112=t112 + x_f(4,i1+l,i2,i3)*a(l) + x_f(5,i1+l,i2,i3)*b(l)
              t121=t121 + x_f(2,i1+l,i2,i3)*a(l) + x_f(3,i1+l,i2,i3)*b(l)
              t122=t122 + x_f(6,i1+l,i2,i3)*a(l) + x_f(7,i1+l,i2,i3)*b(l)
              t212=t212 + x_f(4,i1+l,i2,i3)*c(l) + x_f(5,i1+l,i2,i3)*e(l)
              t221=t221 + x_f(2,i1+l,i2,i3)*c(l) + x_f(3,i1+l,i2,i3)*e(l)
              t222=t222 + x_f(6,i1+l,i2,i3)*c(l) + x_f(7,i1+l,i2,i3)*e(l)
              t211=t211 + x_f(1,i1+l,i2,i3)*e(l)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
           ekin=ekin+t112*x_f(4,i1,i2,i3)
           ekin=ekin+t121*x_f(2,i1,i2,i3)
           ekin=ekin+t211*x_f(1,i1,i2,i3)
           ekin=ekin+t122*x_f(6,i1,i2,i3)
           ekin=ekin+t212*x_f(5,i1,i2,i3)
           ekin=ekin+t221*x_f(3,i1,i2,i3)
           ekin=ekin+t222*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo

  call system_clock(ncount4,ncount_rate,ncount_max)
  tel=dble(ncount4-ncount3)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'T:SECND PART:x',tel,1.d-6*nflop1/tel


  ! + (1/2) d^2/dy^2
  nb=16
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              t112=t112 + x_f(4,i1,i2+l,i3)*a(l) + x_f(6,i1,i2+l,i3)*b(l)
              t211=t211 + x_f(1,i1,i2+l,i3)*a(l) + x_f(3,i1,i2+l,i3)*b(l)
              t122=t122 + x_f(4,i1,i2+l,i3)*c(l) + x_f(6,i1,i2+l,i3)*e(l)
              t212=t212 + x_f(5,i1,i2+l,i3)*a(l) + x_f(7,i1,i2+l,i3)*b(l)
              t221=t221 + x_f(1,i1,i2+l,i3)*c(l) + x_f(3,i1,i2+l,i3)*e(l)
              t222=t222 + x_f(5,i1,i2+l,i3)*c(l) + x_f(7,i1,i2+l,i3)*e(l)
              t121=t121 + x_f(2,i1,i2+l,i3)*e(l)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
           ekin=ekin+t112*x_f(4,i1,i2,i3)
           ekin=ekin+t121*x_f(2,i1,i2,i3)
           ekin=ekin+t211*x_f(1,i1,i2,i3)
           ekin=ekin+t122*x_f(6,i1,i2,i3)
           ekin=ekin+t212*x_f(5,i1,i2,i3)
           ekin=ekin+t221*x_f(3,i1,i2,i3)
           ekin=ekin+t222*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo

  call system_clock(ncount5,ncount_rate,ncount_max)
  tel=dble(ncount5-ncount4)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'T:SECND PART:y',tel,1.d-6*nflop2/tel

  ! + (1/2) d^2/dz^2
  nb=16
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
              t121=t121 + x_f(2,i1,i2,i3+l)*a(l) + x_f(6,i1,i2,i3+l)*b(l)
              t211=t211 + x_f(1,i1,i2,i3+l)*a(l) + x_f(5,i1,i2,i3+l)*b(l)
              t122=t122 + x_f(2,i1,i2,i3+l)*c(l) + x_f(6,i1,i2,i3+l)*e(l)
              t212=t212 + x_f(1,i1,i2,i3+l)*c(l) + x_f(5,i1,i2,i3+l)*e(l)
              t221=t221 + x_f(3,i1,i2,i3+l)*a(l) + x_f(7,i1,i2,i3+l)*b(l)
              t222=t222 + x_f(3,i1,i2,i3+l)*c(l) + x_f(7,i1,i2,i3+l)*e(l)
              t112=t112 + x_f(4,i1,i2,i3+l)*e(l)
           enddo

           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
           ekin=ekin+t112*x_f(4,i1,i2,i3)
           ekin=ekin+t121*x_f(2,i1,i2,i3)
           ekin=ekin+t211*x_f(1,i1,i2,i3)
           ekin=ekin+t122*x_f(6,i1,i2,i3)
           ekin=ekin+t212*x_f(5,i1,i2,i3)
           ekin=ekin+t221*x_f(3,i1,i2,i3)
           ekin=ekin+t222*x_f(7,i1,i2,i3)
        enddo
     enddo
  enddo

  call system_clock(ncount6,ncount_rate,ncount_max)
  tel=dble(ncount6-ncount5)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'T:SECND PART:z',tel,1.d-6*nflop3/tel

  tel=dble(ncount6-ncount0)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'T:ALL   PART',  & 
  !            tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel

  return
END SUBROUTINE ConvolkineticT_prev


!>   y = y + (kinetic energy operator)x 
subroutine ConvolkineticP_prev(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x,y,ekin)
  implicit real(kind=8) (a-h,o-z)
  logical :: firstcall=.true. 
  integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
  dimension x(0:n1,2,0:n2,2,0:n3,2),y(0:n1,2,0:n2,2,0:n3,2)
  dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)


  parameter(lowfil=-14,lupfil=14)
  dimension a(lowfil:lupfil),b(lowfil:lupfil),c(lowfil:lupfil),e(lowfil:lupfil)
  scale=-.5d0/hgrid**2

  !---------------------------------------------------------------------------
  ! second derivative filters for Daubechies 16
  !  <phi|D^2|phi_i>
  a(0)=   -3.5536922899131901941296809374d0*scale
  a(1)=    2.2191465938911163898794546405d0*scale
  a(2)=   -0.6156141465570069496314853949d0*scale
  a(3)=    0.2371780582153805636239247476d0*scale
  a(4)=   -0.0822663999742123340987663521d0*scale
  a(5)=    0.02207029188482255523789911295638968409d0*scale
  a(6)=   -0.409765689342633823899327051188315485d-2*scale
  a(7)=    0.45167920287502235349480037639758496d-3*scale
  a(8)=   -0.2398228524507599670405555359023135d-4*scale
  a(9)=    2.0904234952920365957922889447361d-6*scale
  a(10)=  -3.7230763047369275848791496973044d-7*scale
  a(11)=  -1.05857055496741470373494132287d-8*scale
  a(12)=  -5.813879830282540547959250667d-11*scale
  a(13)=   2.70800493626319438269856689037647576d-13*scale
  a(14)=  -6.924474940639200152025730585882d-18*scale
  do i=1,14
     a(-i)=a(i)
  enddo
  !  <phi|D^2|psi_i>
  c(-14)=     -3.869102413147656535541850057188d-18*scale
  c(-13)=      1.5130616560866154733900029272077362d-13*scale
  c(-12)=     -3.2264702314010525539061647271983988409d-11*scale
  c(-11)=     -5.96264938781402337319841002642d-9*scale
  c(-10)=     -2.1656830629214041470164889350342d-7*scale
  c(-9 )=      8.7969704055286288323596890609625d-7*scale
  c(-8 )=     -0.00001133456724516819987751818232711775d0*scale
  c(-7 )=      0.00021710795484646138591610188464622454d0*scale
  c(-6 )=     -0.0021356291838797986414312219042358542d0*scale
  c(-5 )=      0.00713761218453631422925717625758502986d0*scale
  c(-4 )=     -0.0284696165863973422636410524436931061d0*scale
  c(-3 )=      0.14327329352510759457155821037742893841d0*scale
  c(-2 )=     -0.42498050943780130143385739554118569733d0*scale
  c(-1 )=      0.65703074007121357894896358254040272157d0*scale
  c( 0 )=     -0.42081655293724308770919536332797729898d0*scale
  c( 1 )=     -0.21716117505137104371463587747283267899d0*scale
  c( 2 )=      0.63457035267892488185929915286969303251d0*scale
  c( 3 )=     -0.53298223962800395684936080758073568406d0*scale
  c( 4 )=      0.23370490631751294307619384973520033236d0*scale
  c( 5 )=     -0.05657736973328755112051544344507997075d0*scale
  c( 6 )=      0.0080872029411844780634067667008050127d0*scale
  c( 7 )=     -0.00093423623304808664741804536808932984d0*scale
  c( 8 )=      0.00005075807947289728306309081261461095d0*scale
  c( 9 )=     -4.62561497463184262755416490048242d-6*scale
  c( 10)=      6.3919128513793415587294752371778d-7*scale
  c( 11)=      1.87909235155149902916133888931d-8*scale
  c( 12)=      1.04757345962781829480207861447155543883d-10*scale
  c( 13)=     -4.84665690596158959648731537084025836d-13*scale
  c( 14)=      1.2392629629188986192855777620877d-17*scale
  !  <psi|D^2|phi_i>
  do i=-14,14
     b(i)=c(-i)
  enddo
  !<psi|D^2|psi_i>
  e(0)=   -24.875846029392331358907766562d0*scale
  e(1)=   -7.1440597663471719869313377994d0*scale
  e(2)=   -0.04251705323669172315864542163525830944d0*scale
  e(3)=   -0.26995931336279126953587091167128839196d0*scale
  e(4)=    0.08207454169225172612513390763444496516d0*scale
  e(5)=   -0.02207327034586634477996701627614752761d0*scale
  e(6)=    0.00409765642831595181639002667514310145d0*scale
  e(7)=   -0.00045167920287507774929432548999880117d0*scale
  e(8)=    0.00002398228524507599670405555359023135d0*scale
  e(9)=   -2.0904234952920365957922889447361d-6*scale
  e(10)=   3.7230763047369275848791496973044d-7*scale
  e(11)=   1.05857055496741470373494132287d-8*scale
  e(12)=   5.8138798302825405479592506674648873655d-11*scale
  e(13)=  -2.70800493626319438269856689037647576d-13*scale
  e(14)=   6.924474940639200152025730585882d-18*scale
  do i=1,14
     e(-i)=e(i)
  enddo


  if (firstcall) then

     ! (1/2) d^2/dx^2
     mflop1=0
     do i3=0,n3
        do i2=0,n2
           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
              do l=max(-i1,lowfil),min(lupfil,n1-i1)
                 mflop1=mflop1+4
              enddo
              mflop1=mflop1+4
           enddo

           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
              do l=max(-i1,lowfil),min(lupfil,n1-i1)
                 mflop1=mflop1+3
              enddo
           enddo
        enddo
     enddo

     ! + (1/2) d^2/dy^2
     mflop2=0
     do i3=0,n3
        do i1=0,n1
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
              do l=max(-i2,lowfil),min(lupfil,n2-i2)
                 mflop2=mflop2+4
              enddo
              mflop2=mflop2+4
           enddo

           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
              do l=max(-i2,lowfil),min(lupfil,n2-i2)
                 mflop2=mflop2+3
              enddo
           enddo
        enddo
     enddo


     ! + (1/2) d^2/dz^2
     mflop3=0
     do i2=0,n2
        do i1=0,n1
           do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
              do l=max(-i3,lowfil),min(lupfil,n3-i3)
                 mflop3=mflop3+4
              enddo
              mflop3=mflop3+4
           enddo

           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
              do l=max(-i3,lowfil),min(lupfil,n3-i3)
                 mflop3=mflop3+3
              enddo
           enddo
        enddo
     enddo

     ! wavelet part
     ! (1/2) d^2/dx^2
     nflop1=0
     do i3=nfl3,nfu3
        do i2=nfl2,nfu2
           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
              do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
                 nflop1=nflop1+26
              enddo
              nflop1=nflop1+21
           enddo
        enddo
     enddo

     ! + (1/2) d^2/dy^2
     nflop2=0
     do i3=nfl3,nfu3
        do i1=nfl1,nfu1
           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
              do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
                 nflop2=nflop2+26
              enddo
              nflop2=nflop2+21
           enddo
        enddo
     enddo

     ! + (1/2) d^2/dz^2
     nflop3=0
     do i2=nfl2,nfu2
        do i1=nfl1,nfu1
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
              do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
                 nflop3=nflop3+26
              enddo
              nflop3=nflop3+21
           enddo
        enddo
     enddo

     firstcall=.false.
  endif


  !---------------------------------------------------------------------------

  ekin=0.d0

  ! Scaling function part

  call system_clock(ncount0,ncount_rate,ncount_max)


  ! (1/2) d^2/dx^2
  do i3=0,n3
     do i2=0,n2
        do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
           t111=0.d0 ; s111=0.d0
           do l=max(-i1,lowfil),min(lupfil,n1-i1)
              t111=t111 + x(i1+l,1,i2,1,i3,1)*a(l)
              s111=s111 + x(i1+l,2,i2,1,i3,1)*b(l)
           enddo
           y(i1,1,i2,1,i3,1)=y(i1,1,i2,1,i3,1)+(t111+s111)
           ekin=ekin+(t111+s111)*x(i1,1,i2,1,i3,1)
        enddo

        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t211=0.d0 
           do l=max(-i1,lowfil),min(lupfil,n1-i1)
              t211=t211 + x(i1+l,1,i2,1,i3,1)*c(l)
           enddo
           y(i1,2,i2,1,i3,1)=y(i1,2,i2,1,i3,1)+t211
           ekin=ekin+t211*x(i1,2,i2,1,i3,1)
        enddo
     enddo
  enddo

  call system_clock(ncount1,ncount_rate,ncount_max)
  tel=dble(ncount1-ncount0)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'P:FIRST PART:x',tel,1.d-6*mflop1/tel

  ! + (1/2) d^2/dy^2
  do i3=0,n3
     do i1=0,n1
        do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
           t111=0.d0 ; s111=0.d0
           do l=max(-i2,lowfil),min(lupfil,n2-i2)
              t111=t111 + x(i1,1,i2+l,1,i3,1)*a(l)
              s111=s111 + x(i1,1,i2+l,2,i3,1)*b(l)
           enddo
           y(i1,1,i2,1,i3,1)=y(i1,1,i2,1,i3,1)+(t111+s111)
           ekin=ekin+(t111+s111)*x(i1,1,i2,1,i3,1)
        enddo

        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t121=0.d0 
           do l=max(-i2,lowfil),min(lupfil,n2-i2)
              t121=t121 + x(i1,1,i2+l,1,i3,1)*c(l)
           enddo
           y(i1,1,i2,2,i3,1)=y(i1,1,i2,2,i3,1)+t121
           ekin=ekin+t121*x(i1,1,i2,2,i3,1)
        enddo
     enddo
  enddo


  call system_clock(ncount2,ncount_rate,ncount_max)
  tel=dble(ncount2-ncount1)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'P:FIRST PART:y',tel,1.d-6*mflop2/tel

  ! + (1/2) d^2/dz^2
  do i2=0,n2
     do i1=0,n1
        do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
           t111=0.d0 ; s111=0.d0
           do l=max(-i3,lowfil),min(lupfil,n3-i3)
              t111=t111 + x(i1,1,i2,1,i3+l,1)*a(l)
              s111=s111 + x(i1,1,i2,1,i3+l,2)*b(l)
           enddo
           y(i1,1,i2,1,i3,1)=y(i1,1,i2,1,i3,1)+(t111+s111)
           ekin=ekin+(t111+s111)*x(i1,1,i2,1,i3,1)
        enddo

        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0.d0 
           do l=max(-i3,lowfil),min(lupfil,n3-i3)
              t112=t112 + x(i1,1,i2,1,i3+l,1)*c(l)
           enddo
           y(i1,1,i2,1,i3,2)=y(i1,1,i2,1,i3,2)+t112
           ekin=ekin+t112*x(i1,1,i2,1,i3,2)
        enddo
     enddo
  enddo

  call system_clock(ncount3,ncount_rate,ncount_max)
  tel=dble(ncount3-ncount2)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'P:FIRST PART:z',tel,1.d-6*mflop3/tel

  ! wavelet part
  ! (1/2) d^2/dx^2
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              t112=t112 + x(i1+l,1,i2,1,i3,2)*a(l) + x(i1+l,2,i2,1,i3,2)*b(l)
              t121=t121 + x(i1+l,1,i2,2,i3,1)*a(l) + x(i1+l,2,i2,2,i3,1)*b(l)
              t122=t122 + x(i1+l,1,i2,2,i3,2)*a(l) + x(i1+l,2,i2,2,i3,2)*b(l)
              t212=t212 + x(i1+l,1,i2,1,i3,2)*c(l) + x(i1+l,2,i2,1,i3,2)*e(l)
              t221=t221 + x(i1+l,1,i2,2,i3,1)*c(l) + x(i1+l,2,i2,2,i3,1)*e(l)
              t222=t222 + x(i1+l,1,i2,2,i3,2)*c(l) + x(i1+l,2,i2,2,i3,2)*e(l)
              t211=t211 + x(i1+l,2,i2,1,i3,1)*e(l)
           enddo

           y(i1,1,i2,1,i3,2)=y(i1,1,i2,1,i3,2)+t112
           y(i1,1,i2,2,i3,1)=y(i1,1,i2,2,i3,1)+t121
           y(i1,2,i2,1,i3,1)=y(i1,2,i2,1,i3,1)+t211
           y(i1,1,i2,2,i3,2)=y(i1,1,i2,2,i3,2)+t122
           y(i1,2,i2,1,i3,2)=y(i1,2,i2,1,i3,2)+t212
           y(i1,2,i2,2,i3,1)=y(i1,2,i2,2,i3,1)+t221
           y(i1,2,i2,2,i3,2)=y(i1,2,i2,2,i3,2)+t222
           ekin=ekin+x(i1,1,i2,1,i3,2)*t112
           ekin=ekin+x(i1,1,i2,2,i3,1)*t121
           ekin=ekin+x(i1,2,i2,1,i3,1)*t211
           ekin=ekin+x(i1,1,i2,2,i3,2)*t122
           ekin=ekin+x(i1,2,i2,1,i3,2)*t212
           ekin=ekin+x(i1,2,i2,2,i3,1)*t221
           ekin=ekin+x(i1,2,i2,2,i3,2)*t222
        enddo
     enddo
  enddo

  call system_clock(ncount4,ncount_rate,ncount_max)
  tel=dble(ncount4-ncount3)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'P:SECND PART:x',tel,1.d-6*nflop1/tel


  ! + (1/2) d^2/dy^2
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              t112=t112 + x(i1,1,i2+l,1,i3,2)*a(l) + x(i1,1,i2+l,2,i3,2)*b(l)
              t211=t211 + x(i1,2,i2+l,1,i3,1)*a(l) + x(i1,2,i2+l,2,i3,1)*b(l)
              t122=t122 + x(i1,1,i2+l,1,i3,2)*c(l) + x(i1,1,i2+l,2,i3,2)*e(l)
              t212=t212 + x(i1,2,i2+l,1,i3,2)*a(l) + x(i1,2,i2+l,2,i3,2)*b(l)
              t221=t221 + x(i1,2,i2+l,1,i3,1)*c(l) + x(i1,2,i2+l,2,i3,1)*e(l)
              t222=t222 + x(i1,2,i2+l,1,i3,2)*c(l) + x(i1,2,i2+l,2,i3,2)*e(l)
              t121=t121 + x(i1,1,i2+l,2,i3,1)*e(l)
           enddo

           y(i1,1,i2,1,i3,2)=y(i1,1,i2,1,i3,2)+t112
           y(i1,1,i2,2,i3,1)=y(i1,1,i2,2,i3,1)+t121
           y(i1,2,i2,1,i3,1)=y(i1,2,i2,1,i3,1)+t211
           y(i1,1,i2,2,i3,2)=y(i1,1,i2,2,i3,2)+t122
           y(i1,2,i2,1,i3,2)=y(i1,2,i2,1,i3,2)+t212
           y(i1,2,i2,2,i3,1)=y(i1,2,i2,2,i3,1)+t221
           y(i1,2,i2,2,i3,2)=y(i1,2,i2,2,i3,2)+t222
           ekin=ekin+x(i1,1,i2,1,i3,2)*t112
           ekin=ekin+x(i1,1,i2,2,i3,1)*t121
           ekin=ekin+x(i1,2,i2,1,i3,1)*t211
           ekin=ekin+x(i1,1,i2,2,i3,2)*t122
           ekin=ekin+x(i1,2,i2,1,i3,2)*t212
           ekin=ekin+x(i1,2,i2,2,i3,1)*t221
           ekin=ekin+x(i1,2,i2,2,i3,2)*t222
        enddo
     enddo
  enddo

  call system_clock(ncount5,ncount_rate,ncount_max)
  tel=dble(ncount5-ncount4)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'P:SECND PART:y',tel,1.d-6*nflop2/tel

  ! + (1/2) d^2/dz^2
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
              t121=t121 + x(i1,1,i2,2,i3+l,1)*a(l) + x(i1,1,i2,2,i3+l,2)*b(l)
              t211=t211 + x(i1,2,i2,1,i3+l,1)*a(l) + x(i1,2,i2,1,i3+l,2)*b(l)
              t122=t122 + x(i1,1,i2,2,i3+l,1)*c(l) + x(i1,1,i2,2,i3+l,2)*e(l)
              t212=t212 + x(i1,2,i2,1,i3+l,1)*c(l) + x(i1,2,i2,1,i3+l,2)*e(l)
              t221=t221 + x(i1,2,i2,2,i3+l,1)*a(l) + x(i1,2,i2,2,i3+l,2)*b(l)
              t222=t222 + x(i1,2,i2,2,i3+l,1)*c(l) + x(i1,2,i2,2,i3+l,2)*e(l)
              t112=t112 + x(i1,1,i2,1,i3+l,2)*e(l)
           enddo

           y(i1,1,i2,1,i3,2)=y(i1,1,i2,1,i3,2)+t112
           y(i1,1,i2,2,i3,1)=y(i1,1,i2,2,i3,1)+t121
           y(i1,2,i2,1,i3,1)=y(i1,2,i2,1,i3,1)+t211
           y(i1,1,i2,2,i3,2)=y(i1,1,i2,2,i3,2)+t122
           y(i1,2,i2,1,i3,2)=y(i1,2,i2,1,i3,2)+t212
           y(i1,2,i2,2,i3,1)=y(i1,2,i2,2,i3,1)+t221
           y(i1,2,i2,2,i3,2)=y(i1,2,i2,2,i3,2)+t222
           ekin=ekin+x(i1,1,i2,1,i3,2)*t112
           ekin=ekin+x(i1,1,i2,2,i3,1)*t121
           ekin=ekin+x(i1,2,i2,1,i3,1)*t211
           ekin=ekin+x(i1,1,i2,2,i3,2)*t122
           ekin=ekin+x(i1,2,i2,1,i3,2)*t212
           ekin=ekin+x(i1,2,i2,2,i3,1)*t221
           ekin=ekin+x(i1,2,i2,2,i3,2)*t222

        enddo
     enddo
  enddo

  call system_clock(ncount6,ncount_rate,ncount_max)
  tel=dble(ncount6-ncount5)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'P:SECND PART:z',tel,1.d-6*nflop3/tel

  tel=dble(ncount6-ncount0)/dble(ncount_rate)
  !       write(99,'(a40,2(1x,e10.3))') 'P:ALL   PART',  & 
  !            tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel

  return
END SUBROUTINE ConvolkineticP_prev
