!> @file
!!  Combined convolution routines
!! @deprecated
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

subroutine comb_rot_shrink_loc(ndat,x,y,icf,nfl,nfu,ib)
  ! In one dimension, 
  ! Applies the magic filter transposed, then analysis wavelet transformation.
  ! The size of the data is forced to shrink
  implicit real(kind=8) (a-h,o-z)
  integer,parameter:: lowfil2=-14,lupfil2=16
  dimension x(lowfil2+2*nfl:2*nfu+lupfil2,ndat),y(ndat,nfl:nfu)
  integer ib(2,ndat)
  include 'v_long.inc'

  !open(unit=10,file='longer_filter.flop')
  nflop=0
  ! count the flops:
  do j=1,ndat
     do i=ib(1,j),ib(2,j)
        do l=lowfil2+2*i,lupfil2+2*i
           nflop=nflop+2
        enddo
     enddo
  enddo

  ! the convolution itself:
  call system_clock(ncount0,ncount_rate,ncount_max)
  do j=1,ndat
     if (ib(2,j)-ib(1,j).ge.4) then
        do i=ib(1,j),ib(2,j)-4,4
           ci0=0.d0
           ci1=0.d0
           ci2=0.d0
           ci3=0.d0

           !    l=lowfil2+2*i+0
           !      ci0=ci0+fil2(l-2*(i+0),icf)*x(l,j)
           !    l=lowfil2+2*i+1
           !      ci0=ci0+fil2(l-2*(i+0),icf)*x(l,j)
           !    l=lowfil2+2*i+2
           !      ci0=ci0+fil2(l-2*(i+0),icf)*x(l,j)
           !      ci1=ci1+fil2(l-2*(i+1),icf)*x(l,j)
           !    l=lowfil2+2*i+3
           !      ci0=ci0+fil2(l-2*(i+0),icf)*x(l,j)
           !      ci1=ci1+fil2(l-2*(i+1),icf)*x(l,j)
           !    l=lowfil2+2*i+4
           !      ci0=ci0+fil2(l-2*(i+0),icf)*x(l,j)
           !      ci1=ci1+fil2(l-2*(i+1),icf)*x(l,j)
           !      ci2=ci2+fil2(l-2*(i+2),icf)*x(l,j)
           !    l=lowfil2+2*i+5
           !      ci0=ci0+fil2(l-2*(i+0),icf)*x(l,j)
           !      ci1=ci1+fil2(l-2*(i+1),icf)*x(l,j)
           !      ci2=ci2+fil2(l-2*(i+2),icf)*x(l,j)

           do l=lowfil2+2*i,lupfil2+2*i+6
              ci0=ci0+fil2(l-2*(i+0),icf)*x(l,j)
              ci1=ci1+fil2(l-2*(i+1),icf)*x(l,j)
              ci2=ci2+fil2(l-2*(i+2),icf)*x(l,j)
              ci3=ci3+fil2(l-2*(i+3),icf)*x(l,j)
           enddo

           !    l=lupfil2+2*i+1
           !      ci1=ci1+fil2(l-2*(i+1),icf)*x(l,j)
           !      ci2=ci2+fil2(l-2*(i+2),icf)*x(l,j)
           !      ci3=ci3+fil2(l-2*(i+3),icf)*x(l,j)
           !    l=lupfil2+2*i+2
           !      ci1=ci1+fil2(l-2*(i+1),icf)*x(l,j)
           !      ci2=ci2+fil2(l-2*(i+2),icf)*x(l,j)
           !      ci3=ci3+fil2(l-2*(i+3),icf)*x(l,j)
           !    l=lupfil2+2*i+3
           !      ci2=ci2+fil2(l-2*(i+2),icf)*x(l,j)
           !      ci3=ci3+fil2(l-2*(i+3),icf)*x(l,j)
           !    l=lupfil2+2*i+4
           !      ci2=ci2+fil2(l-2*(i+2),icf)*x(l,j)
           !      ci3=ci3+fil2(l-2*(i+3),icf)*x(l,j)
           !    l=lupfil2+2*i+5
           !      ci3=ci3+fil2(l-2*(i+3),icf)*x(l,j)
           !    l=lupfil2+2*i+6
           !      ci3=ci3+fil2(l-2*(i+3),icf)*x(l,j)

           y(j,i+0)=ci0
           y(j,i+1)=ci1
           y(j,i+2)=ci2
           y(j,i+3)=ci3
        enddo
        icur=i! the greatest multiple of 4 that is smaller or equal to ib(2,j)
     else
        icur=ib(1,j)
     endif

     do i=icur,ib(2,j)
        ci0=0.d0
        do l=lowfil2+2*i,lupfil2+2*i
           ci0=ci0+fil2(l-2*(i+0),icf)*x(l,j)
        enddo
        y(j,i)=ci0
     enddo
  enddo

  call system_clock(ncount1,ncount_rate,ncount_max)
  tel=dble(ncount1-ncount0)/dble(ncount_rate)

  !write(10,*) tel, 1.d-6*nflop/tel
END SUBROUTINE comb_rot_shrink_loc


subroutine comb_rot_shrink_loc_first(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y,icf,ib)
  ! In one dimension, 
  ! Applies the magic filter transposed, then analysis wavelet transformation.
  ! The size of the data is forced to shrink
  implicit real(kind=8) (a-h,o-z)
  integer,parameter:: lowfil2=-14,lupfil2=16
  real(kind=8) x(-14:2*n1+16,-14:2*n2+16,         -14:2*n3+16) ! input
  real(kind=8) y(     -14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16,nfl1:nfu1)! output
  integer ib(2, -14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)
  include 'v_long.inc'

  ! open(unit=10,file='comb_shrink.unrolled')
  nflop=0
  ! count the flops:
  do j2=-14+2*nfl2,2*nfu2+16
     do j3=-14+2*nfl3,2*nfu3+16
        do i=ib(1,j2,j3),ib(2,j2,j3)
           do l=lowfil2+2*i,lupfil2+2*i
              nflop=nflop+2
           enddo
        enddo
     enddo
  enddo

  ! the convolution itself:
  call system_clock(ncount0,ncount_rate,ncount_max)

  do j2=-14+2*nfl2,2*nfu2+16
     do j3=-14+2*nfl3,2*nfu3+16
        if (ib(2,j2,j3)-ib(1,j2,j3).ge.4) then
           do i=ib(1,j2,j3),ib(2,j2,j3)-4,4
              ci0=0.d0
              ci1=0.d0
              ci2=0.d0
              ci3=0.d0

              !      l=lowfil2+2*i+0
              !        ci0=ci0+fil2(l-2*(i+0),icf)*x(l,j2,j3)
              !      l=lowfil2+2*i+1
              !        ci0=ci0+fil2(l-2*(i+0),icf)*x(l,j2,j3)
              !      l=lowfil2+2*i+2
              !        ci0=ci0+fil2(l-2*(i+0),icf)*x(l,j2,j3)
              !        ci1=ci1+fil2(l-2*(i+1),icf)*x(l,j2,j3)
              !      l=lowfil2+2*i+3
              !        ci0=ci0+fil2(l-2*(i+0),icf)*x(l,j2,j3)
              !        ci1=ci1+fil2(l-2*(i+1),icf)*x(l,j2,j3)
              !      l=lowfil2+2*i+4
              !        ci0=ci0+fil2(l-2*(i+0),icf)*x(l,j2,j3)
              !        ci1=ci1+fil2(l-2*(i+1),icf)*x(l,j2,j3)
              !        ci2=ci2+fil2(l-2*(i+2),icf)*x(l,j2,j3)
              !      l=lowfil2+2*i+5
              !        ci0=ci0+fil2(l-2*(i+0),icf)*x(l,j2,j3)
              !        ci1=ci1+fil2(l-2*(i+1),icf)*x(l,j2,j3)
              !        ci2=ci2+fil2(l-2*(i+2),icf)*x(l,j2,j3)

              do l=lowfil2+2*i,lupfil2+2*i+6
                 ci0=ci0+fil2(l-2*(i+0),icf)*x(l,j2,j3)
                 ci1=ci1+fil2(l-2*(i+1),icf)*x(l,j2,j3)
                 ci2=ci2+fil2(l-2*(i+2),icf)*x(l,j2,j3)
                 ci3=ci3+fil2(l-2*(i+3),icf)*x(l,j2,j3)
              enddo

              !      l=lupfil2+2*i+1
              !        ci1=ci1+fil2(l-2*(i+1),icf)*x(l,j2,j3)
              !        ci2=ci2+fil2(l-2*(i+2),icf)*x(l,j2,j3)
              !        ci3=ci3+fil2(l-2*(i+3),icf)*x(l,j2,j3)
              !      l=lupfil2+2*i+2
              !        ci1=ci1+fil2(l-2*(i+1),icf)*x(l,j2,j3)
              !        ci2=ci2+fil2(l-2*(i+2),icf)*x(l,j2,j3)
              !        ci3=ci3+fil2(l-2*(i+3),icf)*x(l,j2,j3)
              !      l=lupfil2+2*i+3
              !        ci2=ci2+fil2(l-2*(i+2),icf)*x(l,j2,j3)
              !        ci3=ci3+fil2(l-2*(i+3),icf)*x(l,j2,j3)
              !      l=lupfil2+2*i+4
              !        ci2=ci2+fil2(l-2*(i+2),icf)*x(l,j2,j3)
              !        ci3=ci3+fil2(l-2*(i+3),icf)*x(l,j2,j3)
              !      l=lupfil2+2*i+5
              !        ci3=ci3+fil2(l-2*(i+3),icf)*x(l,j2,j3)
              !      l=lupfil2+2*i+6
              !        ci3=ci3+fil2(l-2*(i+3),icf)*x(l,j2,j3)

              y(j2,j3,i+0)=ci0
              y(j2,j3,i+1)=ci1
              y(j2,j3,i+2)=ci2
              y(j2,j3,i+3)=ci3
           enddo

           icur=i
        else
           icur=ib(1,j2,j3)
        endif

        do i=icur,ib(2,j2,j3)
           ci0=0.d0
           do l=lowfil2+2*i,lupfil2+2*i
              ci0=ci0+fil2(l-2*(i+0),icf)*x(l,j2,j3)
           enddo
           y(j2,j3,i)=ci0
        enddo
     enddo
  enddo

  call system_clock(ncount1,ncount_rate,ncount_max)
  tel=dble(ncount1-ncount0)/dble(ncount_rate)

  !write(10,*) tel, 1.d-6*nflop/tel

END SUBROUTINE comb_rot_shrink_loc_first

