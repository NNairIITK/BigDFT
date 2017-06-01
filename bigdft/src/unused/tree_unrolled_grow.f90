!> @file
!!   Old convolution routines
!! @deprecated
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

subroutine  comb_rot_grow_loc_1(nfl,nfu,ndat,x,y,ib)
  ! In one dimension, 
  ! with optimised cycles
  ! Applies synthesis wavelet transformation 
  ! then convolves with magic filter
  !  the size of the data is allowed to grow

  implicit real(kind=8) (a-h,o-z)
  integer t
  integer,parameter:: lowfil=-7,lupfil=8
  integer,parameter:: lowfil2=2*lowfil,lupfil2=2*lupfil
  dimension x(7,nfl:nfu,ndat),y(2,2,ndat,-14+2*nfl:2*nfu+16)
  integer ib(2,ndat)

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
END SUBROUTINE comb_rot_grow_loc_1



subroutine  comb_rot_grow_loc_2(nfl,nfu,ndat,x,y,ib)
! In one dimension, 
! with optimised cycles
! Applies synthesis wavelet transformation 
! then convolves with magic filter
! the size of the data is allowed to grow

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
END SUBROUTINE comb_rot_grow_loc_2



subroutine comb_rot_grow_loc_3(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y,ib)
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

end



