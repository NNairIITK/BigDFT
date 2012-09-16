!> @file
!!  Combined magic filter + analysis routines
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> In 3d,   
!! Applies the magic filter transposed, then analysis wavelet transformation.
!! The size of the data is forced to shrink
!! The input array y is not overwritten
subroutine comb_shrink(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,y,&
     ibxy_c,ibzzx_c,ibyyzz_c,ibxy_f,ibzzx_f,ibyyzz_f,xc,xf,ibyz_c,ibyz_f)

  implicit none
  integer :: n1,n2,n3,i1,i2,i3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(kind=8) :: y( -14:2*n1+16,-14:2*n2+16,-14:2*n3+16) ! input
  !    real(kind=8) w1( 2,           -14:2*n2+16,-14:2*n3+16,0:n1)!  work 
  ! real(kind=8) w2( 4,                       -14:2*n3+16,0:n1,0:n2)
  real(kind=8) :: w1(max(2*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31)*(nfu1-nfl1+1),&
       (2*n2+31)*(2*n3+31)*(n1+1)))
  real(kind=8) :: w2( max(4*(2*(nfu3-nfl3)+31)*(nfu1-nfl1+1)*(nfu2-nfl2+1),&
       (2*n3+31)*(n1+1)*(n2+1)))

  real(kind=8) :: xf(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),xc(0:n1,0:n2,0:n3)! output arrays

  ! boundary arrays
  integer :: ibxy_c(2,0:n1,0:n2) 
  integer :: ibzzx_c(2,-14:2*n3+16,0:n1) 
  integer :: ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16)

  integer :: ibxy_f(2,nfl1:nfu1,nfl2:nfu2)
  integer :: ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1)
  integer :: ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)

  integer, dimension(2,0:n2,0:n3) :: ibyz_c,ibyz_f

  ! perform the combined transform 
  call comb_shrink_loc_c(0,n1,0,n2,0,n3,w1,w2,y,xc,1,1,1,&
       ibxy_c,ibzzx_c,ibyyzz_c) ! for scfunctions
  ! for wavelets:

  call comb_shrink_loc_f(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,y,xf,&
       ibxy_f,ibzzx_f,ibyyzz_f)

END SUBROUTINE comb_shrink


!> In 3d,   
!! Applies the magic filter transposed, then analysis wavelet transformation.
!! The output is only the l1,l2,l3 wavelet component
!! The size of the data is forced to shrink
!! The input array y is not overwritten
subroutine comb_shrink_loc_f(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,y,x,&
     ibxy,ibzzx,ibyyzz)
  implicit real(kind=8) (a-h,o-z)
  real(kind=8) y(-14:2*n1+16,-14:2*n2+16,         -14:2*n3+16) ! input
  real(kind=8) w1(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16,nfl1:nfu1)!work
  real(kind=8) w2(4,-14+2*nfl3:2*nfu3+16,nfl1:nfu1,nfl2:nfu2)
  real(kind=8) x(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)!output

  integer ibxy(2,nfl1:nfu1,nfl2:nfu2)
  integer ibzzx(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1)
  integer ibyyzz(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)

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
  return
END SUBROUTINE comb_shrink_loc_f


!> In 3d,   
!! Applies the magic filter transposed, then analysis wavelet transformation.
!! The output is only the l1,l2,l3 wavelet component
!! The size of the data is forced to shrink
!! The input array y is not overwritten
subroutine comb_shrink_loc_c(nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,y,x,l1,l2,l3,&
     ibxy,ibzzx,ibyyzz)
  implicit none
  !Arguments
  integer, intent(in) :: l1,l2,l3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
  real(kind=8) :: y(-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)!input
  real(kind=8) :: w1(-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16,nfl1:nfu1)!work
  real(kind=8) :: w2(-14+2*nfl3:2*nfu3+16,nfl1:nfu1,nfl2:nfu2)
  real(kind=8) :: x(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)!output

  integer :: ibxy(2,nfl1:nfu1,nfl2:nfu2)
  integer :: ibzzx(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1)
  integer :: ibyyzz(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)

  !Local variables
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

  return
END SUBROUTINE comb_shrink_loc_c

! subroutine for scfunctions:
!include 'simple_shrink.f90'
!include 'unrolled_shrink.f90'
!Better to link it
!include 'long_shrink.f90'

!subroutine for wavelets:
!include 'tree.f90'
!Better to link it
!include 'tree_long_shrink.f90'
