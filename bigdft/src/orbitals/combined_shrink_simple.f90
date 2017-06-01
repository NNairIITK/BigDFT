!> @file
!!  Simple routines combining magic filters and analysis wavelet transforms
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> In one dimension,    
!! Applies the magic filter transposed, then analysis wavelet transformation.
!! The size of the data is forced to shrink
subroutine comb_rot_shrink_loc(ndat,x,y,icf,nfl,nfu,ib)
  use module_defs, only: wp
  implicit none
  integer, parameter :: lowfil2=-14,lupfil2=16
  integer, intent(in) :: ndat,icf,nfl,nfu
  integer, dimension(2,ndat), intent(in) :: ib
  real(wp), dimension(lowfil2+2*nfl:2*nfu+lupfil2,ndat), intent(in) :: x
  real(wp), dimension(ndat,nfl:nfu), intent(out) :: y
  !local variables
  integer :: j,l,i !n(c) ,icur
  real(wp) :: ci !n(c) ,ci1,ci2,ci3
  include 'v.inc'

!    open(unit=10,file='simple_shrink.flop')
!    nflop=0
!    ! count the flops:
!    do j=1,ndat
!           do i=ib(1,j),ib(2,j)
!             do l=lowfil2+2*i,lupfil2+2*i
!                nflop=nflop+2
!             enddo
!           enddo
!    enddo

    ! the convolution itself:
  !call system_clock(ncount0,ncount_rate,ncount_max)
    do j=1,ndat
       do i=ib(1,j),ib(2,j)
         ci=0.d0
         do l=lowfil2+2*i,lupfil2+2*i
           ci=ci+fil2(l-2*i,icf)*x(l,j)
         enddo
         y(j,i)=ci
       enddo
    enddo
       
!    call system_clock(ncount1,ncount_rate,ncount_max)
!    tel=dble(ncount1-ncount0)/dble(ncount_rate)

!    write(10,*) tel, 1.d-6*nflop/tel
END SUBROUTINE comb_rot_shrink_loc


!> In one dimension,    
!! Applies the magic filter transposed, then analysis wavelet transformation.
!! The size of the data is forced to shrink
subroutine comb_rot_shrink_loc_1(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y,ib)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, dimension(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16), intent(in) :: ib
  real(wp), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(in) :: x
  real(wp), dimension(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16,nfl1:nfu1), intent(out) :: y
  !local variables
  integer, parameter :: lowfil2=-14,lupfil2=16
  integer :: i,j2,j3,l !n(c) nflop,icur
  real(wp) :: ci1,ci2 !n(c) c1i0,c1i1,c1i2,c1i3,c2i0,c2i1,c2i2,c2i3
  include 'v.inc'

    !n(c) nflop=0
!    open(unit=20,file='tree_shrink.flop')
!    call system_clock(ncount0,ncount_rate,ncount_max)

    do j2=-14+2*nfl2,2*nfu2+16
        do j3=-14+2*nfl3,2*nfu3+16
           do i=ib(1,j2,j3),ib(2,j2,j3)
             ci1=0.d0
             ci2=0.d0
             !nflop=nflop+(lupfil2-lowfil2+1)*2*2
             do l=lowfil2+2*i,lupfil2+2*i
               ci1=ci1+fil2(l-2*i,1)*x(l,j2,j3)
               ci2=ci2+fil2(l-2*i,2)*x(l,j2,j3)
             enddo
             y(1,j2,j3,i)=ci1
             y(2,j2,j3,i)=ci2
           enddo
        enddo
    enddo

    !call system_clock(ncount1,ncount_rate,ncount_max)
    !tel=dble(ncount1-ncount0)/dble(ncount_rate)
    !write(20,*) tel, 1.d-6*nflop/tel

END SUBROUTINE comb_rot_shrink_loc_1


!> In one dimension,    
!! Applies the magic filter transposed, then analysis wavelet transformation.
!! The size of the data is forced to shrink
subroutine comb_rot_shrink_loc_2(ndat,x,y,nfl,nfu,ib)
  use module_defs, only: wp
  implicit none
  integer, parameter:: lowfil2=-14,lupfil2=16
  integer, intent(in) :: ndat,nfl,nfu
  integer, dimension(2,ndat), intent(in) :: ib
  real(wp), dimension(2,lowfil2+2*nfl:2*nfu+lupfil2,ndat), intent(in) :: x
  real(wp), dimension(2,2,ndat,nfl:nfu), intent(out) :: y
  !local variables
  integer :: j,i,l !n(c) nflop,icur
  real(wp) :: ci11,ci12,ci21,ci22 !n(c) c11i0,c12i0,c21i0,c22i0,c11i1,c12i1,c21i1,c22i1
  include 'v.inc'

    !n(c) nflop=0
    !open(unit=20,file='tree_shrink.flop')
    !call system_clock(ncount0,ncount_rate,ncount_max)

    do j=1,ndat
       do i=ib(1,j),ib(2,j)
         ci11=0.d0
         ci12=0.d0
         ci21=0.d0
         ci22=0.d0

         !nflop=nflop+(lupfil2-lowfil2+1)*2*4
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

    !call system_clock(ncount1,ncount_rate,ncount_max)
    !tel=dble(ncount1-ncount0)/dble(ncount_rate)
    !write(20,*) tel, 1.d-6*nflop/tel
       
END SUBROUTINE comb_rot_shrink_loc_2
