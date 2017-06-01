!> @file
!!  Old combined convolution routines
!! @deprecated
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

subroutine comb_rot_shrink_loc_1(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y,ib)
  ! In one dimension, 
  ! Applies the magic filter transposed, then analysis wavelet transformation.
  ! The size of the data is forced to shrink
  implicit none
!Arguments 
  integer :: n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,ib(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)
  real(kind=8) :: x(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16),y(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16,nfl1:nfu1)
!Local variables 
  integer :: i,icur,j2,j3,l
  integer, parameter :: lowfil2=-14,lupfil2=16
  integer :: ncount0,ncount1,ncount_max,ncount_rate,nflop
  real(kind=8) :: c1i0,c1i1,c1i2,c1i3,c2i0,c2i1,c2i2,c2i3,ci1,ci2,tel

  include 'v_long.inc'

  nflop=0
  !open(unit=20,file='long.flop')
  call system_clock(ncount0,ncount_rate,ncount_max)

  do j2=-14+2*nfl2,2*nfu2+16
     do j3=-14+2*nfl3,2*nfu3+16
        if (ib(2,j2,j3)-ib(1,j2,j3).ge.4) then
           do i=ib(1,j2,j3),ib(2,j2,j3)-4,4
              c1i0=0.d0
              c1i1=0.d0
              c1i2=0.d0
              c1i3=0.d0

              c2i0=0.d0
              c2i1=0.d0
              c2i2=0.d0
              c2i3=0.d0

              nflop=nflop+(lupfil2-lowfil2+1)*8*2

              do l=lowfil2+2*i,lupfil2+2*i+6
                 c1i0=c1i0+fil2(l-2*(i+0),1)*x(l,j2,j3)
                 c1i1=c1i1+fil2(l-2*(i+1),1)*x(l,j2,j3)
                 c1i2=c1i2+fil2(l-2*(i+2),1)*x(l,j2,j3)
                 c1i3=c1i3+fil2(l-2*(i+3),1)*x(l,j2,j3)

                 c2i0=c2i0+fil2(l-2*(i+0),2)*x(l,j2,j3)
                 c2i1=c2i1+fil2(l-2*(i+1),2)*x(l,j2,j3)
                 c2i2=c2i2+fil2(l-2*(i+2),2)*x(l,j2,j3)
                 c2i3=c2i3+fil2(l-2*(i+3),2)*x(l,j2,j3)
              enddo

              y(1,j2,j3,i+0)=c1i0
              y(1,j2,j3,i+1)=c1i1
              y(1,j2,j3,i+2)=c1i2
              y(1,j2,j3,i+3)=c1i3

              y(2,j2,j3,i+0)=c2i0
              y(2,j2,j3,i+1)=c2i1
              y(2,j2,j3,i+2)=c2i2
              y(2,j2,j3,i+3)=c2i3
           enddo
           icur=i
        else
           icur=ib(1,j2,j3)
        endif

        do i=icur,ib(2,j2,j3)
           ci1=0.d0
           ci2=0.d0
           nflop=nflop+(lupfil2-lowfil2+1)*2*2
           do l=lowfil2+2*i,lupfil2+2*i
              ci1=ci1+fil2(l-2*i,1)*x(l,j2,j3)
              ci2=ci2+fil2(l-2*i,2)*x(l,j2,j3)
           enddo
           y(1,j2,j3,i)=ci1
           y(2,j2,j3,i)=ci2
        enddo
     enddo
  enddo

  call system_clock(ncount1,ncount_rate,ncount_max)
  tel=dble(ncount1-ncount0)/dble(ncount_rate)
  ! write(20,*) tel, 1.d-6*nflop/tel

END SUBROUTINE comb_rot_shrink_loc_1


subroutine comb_rot_shrink_loc_2(ndat,x,y,nfl,nfu,ib)
  ! In one dimension, 
  ! Applies the magic filter transposed, then analysis wavelet transformation.
  ! The size of the data is forced to shrink
  implicit real(kind=8) (a-h,o-z)
  integer,parameter:: lowfil2=-14,lupfil2=16
  dimension x(2,lowfil2+2*nfl:2*nfu+lupfil2,ndat),y(2,2,ndat,nfl:nfu)
  integer ib(2,ndat)
  include 'v_long.inc'

  nflop=0
  !open(unit=20,file='long.flop')
  call system_clock(ncount0,ncount_rate,ncount_max)

  do j=1,ndat
     if (ib(2,j)-ib(1,j).ge.2) then
        do i=ib(1,j),ib(2,j)-2,2

           c11i0=0.d0
           c12i0=0.d0
           c21i0=0.d0
           c22i0=0.d0

           c11i1=0.d0
           c12i1=0.d0
           c21i1=0.d0
           c22i1=0.d0

           nflop=nflop+(lupfil2-lowfil2+1)*2*8
           do l=lowfil2+2*i,lupfil2+2*i+2
              c11i0=c11i0+fil2(l-2*(i+0),1)*x(1,l,j)
              c12i0=c12i0+fil2(l-2*(i+0),2)*x(1,l,j)
              c21i0=c21i0+fil2(l-2*(i+0),1)*x(2,l,j)
              c22i0=c22i0+fil2(l-2*(i+0),2)*x(2,l,j)

              c11i1=c11i1+fil2(l-2*(i+1),1)*x(1,l,j)
              c12i1=c12i1+fil2(l-2*(i+1),2)*x(1,l,j)
              c21i1=c21i1+fil2(l-2*(i+1),1)*x(2,l,j)
              c22i1=c22i1+fil2(l-2*(i+1),2)*x(2,l,j)
           enddo

           y(1,1,j,i+0)=c11i0
           y(1,2,j,i+0)=c12i0
           y(2,1,j,i+0)=c21i0
           y(2,2,j,i+0)=c22i0

           y(1,1,j,i+1)=c11i1
           y(1,2,j,i+1)=c12i1
           y(2,1,j,i+1)=c21i1
           y(2,2,j,i+1)=c22i1
        enddo
        icur=i
     else
        icur=ib(1,j)
     endif

     do i=icur,ib(2,j)
        ci11=0.d0
        ci12=0.d0
        ci21=0.d0
        ci22=0.d0

        nflop=nflop+(lupfil2-lowfil2+1)*2*4
        do l=lowfil2+2*i,lupfil2+2*i
           ci11=ci11+fil2(l-2*i,1)*x(1,l,j)
           ci12=ci12+fil2(l-2*i,2)*x(1,l,j)
           ci21=ci21+fil2(l-2*i,1)*x(2,l,j)
           ci22=ci22+fil2(l-2*i,2)*x(2,l,j)
        enddo
        y(1,1,j,i)=ci11
        y(1,2,j,i)=ci12
        y(2,1,j,i)=ci21
        y(2,2,j,i)=ci22
     enddo
  enddo

  call system_clock(ncount1,ncount_rate,ncount_max)
  tel=dble(ncount1-ncount0)/dble(ncount_rate)
  ! write(20,*) tel, 1.d-6*nflop/tel

END SUBROUTINE comb_rot_shrink_loc_2



subroutine comb_rot_shrink_loc_3(ndat,x,y,nfl,nfu,ib)
  ! In one dimension, 
  ! Applies the magic filter transposed, then analysis wavelet transformation.
  ! The size of the data is forced to shrink
  implicit real(kind=8) (a-h,o-z)
  integer,parameter:: lowfil2=-14,lupfil2=16
  dimension x(2,2,lowfil2+2*nfl:2*nfu+lupfil2,ndat),y(7,ndat,nfl:nfu)
  integer ib(2,ndat)
  include 'v.inc'

  nflop=0
  !open(unit=20,file='long.flop')
  call system_clock(ncount0,ncount_rate,ncount_max)

  do j=1,ndat
     do i=ib(1,j),ib(2,j)
        ci112=0.d0
        ci121=0.d0
        ci122=0.d0
        ci211=0.d0
        ci212=0.d0
        ci221=0.d0
        ci222=0.d0

        nflop=nflop+(lupfil2-lowfil2+1)*2*7
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

  call system_clock(ncount1,ncount_rate,ncount_max)
  tel=dble(ncount1-ncount0)/dble(ncount_rate)
  ! write(20,*) tel, 1.d-6*nflop/tel

END SUBROUTINE comb_rot_shrink_loc_3
