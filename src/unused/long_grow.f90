!> @file
!!  Combined convolution routines (unused)
!! @deprecated
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

subroutine  comb_rot_grow_loc_plus(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y,icf,ib)
  ! In one dimension, 
  ! with optimised cycles
  ! Applies synthesis wavelet transformation 
  ! then convolves with magic filter
  ! then adds the result to y.
  ! The size of the data is allowed to grow

  implicit real(kind=8) (a-h,o-z)
  integer t
  integer,parameter:: lowfil=-7,lupfil=8
  integer,parameter:: lowfil2=2*lowfil,lupfil2=2*lupfil

  dimension  x(nfl3:nfu3,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16)
  dimension  y(          -14       :2*n1  +16,-14       :2*n2  +16,-14:2*n3+16)
  integer ib(2,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16)

  include 'v.inc'

  do l1=-14+2*nfl1,2*nfu1+16
     do l2=-14+2*nfl2,2*nfu2+16
 !       loop for smaller i
 !       j=8 won't work : i-j would turn negative
        if (ib(1,l1,l2).le.ib(2,l1,l2)) then
           do i=ib(1,l1,l2)-7,ib(1,l1,l2)+7
              y2i=0.d0
              y2i1=0.d0

              ! i-7 =< ib(1,l1,l2)
              do t=      ib(1,l1,l2) ,min(i+7,ib(2,l1,l2))
                 y2i =y2i +fil2(2*(i-t)  ,icf)*x(t,l1,l2)
                 y2i1=y2i1+fil2(2*(i-t)+1,icf)*x(t,l1,l2)
              enddo

              y(l1,l2,2*i  )=y(l1,l2,2*i  )+y2i
              y(l1,l2,2*i+1)=y(l1,l2,2*i+1)+y2i1
           enddo

           !   loop for ordinary i  
           do i=ib(1,l1,l2)+8,ib(2,l1,l2)+7

              ! j=8  works since i is big enough
              y2i =fil2(16  ,icf)*x(i-8,l1,l2)
              y2i1=0.d0

              ! i-7>ib(1,l) since i=ib(1,l)+8,..
              do t=    i-7             ,min(i+7,ib(2,l1,l2)) 
                 y2i =y2i +fil2(2*(i-t)  ,icf)*x(t,l1,l2)
                 y2i1=y2i1+fil2(2*(i-t)+1,icf)*x(t,l1,l2)
              enddo

              y(l1,l2,2*i  )=y(l1,l2,2*i  )+y2i
              y(l1,l2,2*i+1)=y(l1,l2,2*i+1)+y2i1
           enddo

           t=ib(2,l1,l2); i=t+8;  ! to the rightmost element of y, only j=8 contributes
           y(l1,l2,2*i)=y(l1,l2,2*i)+fil2(16,icf)*x(t,l1,l2)

        endif
     enddo
  enddo

END SUBROUTINE comb_rot_grow_loc_plus



subroutine comb_rot_grow_loc(nfl,nfu,ndat,x,y,icf,ib)
! In one dimension, 
! with optimised cycles
! Applies synthesis wavelet transformation 
! then convolves with magic filter
!  the size of the data is allowed to grow

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
END SUBROUTINE comb_rot_grow_loc



