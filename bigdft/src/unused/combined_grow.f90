!> @file
!!  Convolution routines
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!> w2 and w1 are switched because from shrink convention we go
!! to grow convention
subroutine comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3&
     ,w2,w1,xc,xf,y,ibyz_c,ibzxx_c,ibxxyy_c,&
     ibyz_f,ibzxx_f,ibxxyy_f)
  implicit none
  
  integer n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  
  real(kind=8),intent(in) :: xc(0:n1,0:n2,0:n3),  xf(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)! input
  real(kind=8),intent(out) :: y(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)! output
  
  real(kind=8) :: w1(4,nfl2:nfu2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16)!work
  !   real(kind=8) w2(0:n3,-14:2*n1+16,-14:2*n2+16) ! work
  ! real(kind=8) w2(2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16)
  real(kind=8) :: w2(max((n3+1)*(2*n1+31)*(2*n2+31),&
       2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))) ! work
  
  integer,intent(in)::ibyz_c(2,0:n2,0:n3)
  integer,intent(in)::ibzxx_c(2,0:n3,-14:2*n1+16)
  integer,intent(in)::ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)
  
  integer,intent(in)::ibyz_f(2,nfl2:nfu2,nfl3:nfu3) 
  integer,intent(in)::ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
  integer,intent(in)::ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)
  
  call comb_grow_c(n1,n2,n3,w2,xc,y,ibyz_c,ibzxx_c,ibxxyy_c)
  
  call comb_grow_tree(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3&
       ,w1,w2,xf,y,ibyz_f,ibzxx_f,ibxxyy_f)    
  
END SUBROUTINE comb_grow_all


!> In 3d,   
!! Applies synthesis wavelet transformation 
!! then convolves with magic filter
!!  the size of the data is allowed to grow
!! The input array x is not overwritten
subroutine comb_grow_tree(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3&
     ,w1,w2,x,y,ibyz,ibzxx,ibxxyy)

  implicit none

!Arguments
  integer :: n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
  real(kind=8) :: x(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3) !in, fine
  real(kind=8) :: w1(4,nfl2:nfu2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16)
  real(kind=8) :: w2(2,nfl3:nfu3,-14+2*nfl1:2*nfu1+16,-14+2*nfl2:2*nfu2+16)
  real(kind=8) :: y(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16) !out
  integer :: ibyz(2,nfl2:nfu2,nfl3:nfu3)
  integer :: ibzxx(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
  integer :: ibxxyy(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)
!Local variables
  integer :: m1,m2,m3,nt

  m1=nfu1-nfl1
  m2=nfu2-nfl2
  m3=nfu3-nfl3

  ! i1,i2,i3 -> i2,i3,I1
  nt=(nfu2-nfl2+1)*(nfu3-nfl3+1)
  call comb_rot_grow_loc_1(nfl1,nfu1,nt,x,w1,ibyz) 
  
  ! i2,i3,I1 -> i3,I1,I2
  nt=(nfu3-nfl3+1)*(2*m1+31)
  call comb_rot_grow_loc_2(nfl2,nfu2,nt,w1,w2,ibzxx) 
  
  ! i3,I1,I2  -> I1,I2,I3: add the result to y
  call comb_rot_grow_loc_3(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w2,y,ibxxyy)
  
END SUBROUTINE comb_grow_tree


!> In 3d,   
!! Applies synthesis wavelet transformation 
!! then convolves with magic filter
!!  the size of the data is allowed to grow
!! The input array x is not overwritten
!! However, the output array y contains nonphysical values
!! outside of the localization region
!! that remain from the first comb_grow
subroutine comb_grow_c(n1,n2,n3,ww,x,y,ibyz,ibzxx,ibxxyy)
  
  implicit real(kind=8) (a-h,o-z)
  real(kind=8) :: x(0:n1,0:n2,0:n3) !in
  real(kind=8) :: ww(0:n3,-14:2*n1+16,-14:2*n2+16) ! work
  real(kind=8) :: y(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16) !out
  integer :: ibyz(2,0:n2,0:n3)
  integer :: ibzxx(2,0:n3,-14:2*n1+16)
  integer :: ibxxyy(2,-14:2*n1+16,-14:2*n2+16)
  
  ! i1,i2,i3 -> i2,i3,I1
  nt=(n2+1)*(n3+1)
  call comb_rot_grow_loc(0,n1,nt,x,y,1,ibyz) 
  
  ! i2,i3,I1 -> i3,I1,I2
  nt=(n3+1)*(2*n1+31)
  call comb_rot_grow_loc(0,n2,nt,y,ww,1,ibzxx) 
  ! dimension of ww: nt*(2*n2+31)=(n3+1)*(2*n1+31)*(2*n2+31)
  
  ! i3,I1,I2  -> I1,I2,I3
  nt=(2*n1+31)*(2*n2+31)
  call comb_rot_grow_loc(0,n3,nt,ww,y,1,ibxxyy)
  
END SUBROUTINE comb_grow_c

! subroutine for scfunctions
!include 'standard_grow.f90'
!include 'unrolled_grow.f90'
!Better to link it
!include 'long_grow.f90'

! subroutine for wavelets
!include 'tree_grow.f90'
!Better to link it
!include 'tree_unrolled_grow.f90'

