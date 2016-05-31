!> @file
!! Simple convolution routines
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


subroutine comb_grow_all_hybrid(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nw1,nw2&
     ,w1,w2,xc,xf,y)
use module_defs, only: wp
implicit none
integer,intent(in)::n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nw1,nw2
real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: xc
real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: xf
real(wp), dimension(nw1), intent(inout) :: w1 !work
real(wp), dimension(nw2), intent(inout) :: w2 ! work
real(wp), dimension(0:2*n1+1,0:2*n2+1,0:2*n3+1), intent(out) :: y

call comb_grow_c_simple(n1,n2,n3,w1,w2,xc,y)

call comb_rot_grow_1(n1      ,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,xf,w1)
call comb_rot_grow_2(n1,n2   ,nfl2,nfu2,nfl3,nfu3,w1,w2)
call comb_rot_grow_3(n1,n2,n3,nfl3,nfu3,w2,y)

END SUBROUTINE comb_grow_all_hybrid


!> In one dimesnion,    
!! with optimised cycles
!! Applies synthesis wavelet transformation 
!! then convolves with magic filter
!! the size of the data is allowed to grow
subroutine comb_rot_grow_1(n1,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y)
use module_defs, only: wp
implicit none
integer, intent(in) :: n1,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x
real(wp), dimension(2,2,nfl2:nfu2,nfl3:nfu3,0:2*n1+1), intent(out) :: y
!local variables
integer :: l2,l3,i,t !n(c) ,l1
integer :: ii,ii1
real(wp) y2i__11,y2i__21,y2i1_11,y2i1_21
real(wp) y2i__12,y2i__22,y2i1_12,y2i1_22
integer::modul(-14+2*nfl1:2*nfu1+16)

include 'v_17.inc'
call fill_mod_arr(modul,-14+2*nfl1,2*nfu1+16,2*n1+2)

do l3=nfl3,nfu3
   do l2=nfl2,nfu2

      ! the parts of y that are not initialized, should be put to zero.
      do ii=0,2*nfl1-15
         y(:,:,l2,l3,ii)=0.d0
      enddo
      
      do ii=2*nfu1+17,2*n1+1
         y(:,:,l2,l3,ii)=0.d0
      enddo
      
      ii=modul(2*nfu1+16)
      y(1,1,l2,l3,ii)=fil2(16,2)*x(1,nfu1,l2,l3)
      y(2,1,l2,l3,ii)=fil2(16,1)*x(2,nfu1,l2,l3)+fil2(16,2)*x(3,nfu1,l2,l3)
      
      y(1,2,l2,l3,ii)=fil2(16,1)*x(4,nfu1,l2,l3)+fil2(16,2)*x(5,nfu1,l2,l3)
      y(2,2,l2,l3,ii)=fil2(16,1)*x(6,nfu1,l2,l3)+fil2(16,2)*x(7,nfu1,l2,l3)
   
      ! scaling function everywhere, but folded back to the box
      do i=nfl1-7,nfu1+7
         y2i__11=0.d0
         y2i__21=0.d0
         y2i__12=0.d0
         y2i__22=0.d0
         
         y2i1_11=0.d0
         y2i1_21=0.d0
         y2i1_12=0.d0
         y2i1_22=0.d0
      
         ! wavelet restricted   
         do t=max(i-8,nfl1),min(i+7,nfu1)
            y2i__11=y2i__11                                +fil2(2*(i-t) ,2)*x(1,t,l2,l3)
            y2i__21=y2i__21+fil2(2*(i-t)  ,1)*x(2,t,l2,l3)+fil2(2*(i-t)  ,2)*x(3,t,l2,l3)
            y2i__12=y2i__12+fil2(2*(i-t)  ,1)*x(4,t,l2,l3)+fil2(2*(i-t)  ,2)*x(5,t,l2,l3)
            y2i__22=y2i__22+fil2(2*(i-t)  ,1)*x(6,t,l2,l3)+fil2(2*(i-t)  ,2)*x(7,t,l2,l3)
                                                                                  
            y2i1_11=y2i1_11                               +fil2(2*(i-t)+1,2)*x(1,t,l2,l3)
            y2i1_21=y2i1_21+fil2(2*(i-t)+1,1)*x(2,t,l2,l3)+fil2(2*(i-t)+1,2)*x(3,t,l2,l3)
            y2i1_12=y2i1_12+fil2(2*(i-t)+1,1)*x(4,t,l2,l3)+fil2(2*(i-t)+1,2)*x(5,t,l2,l3)
            y2i1_22=y2i1_22+fil2(2*(i-t)+1,1)*x(6,t,l2,l3)+fil2(2*(i-t)+1,2)*x(7,t,l2,l3)
         enddo
         
         ii =modul(2*i  )
         ii1=modul(2*i+1)
         
         y(1,1,l2,l3,ii )=y2i__11
         y(2,1,l2,l3,ii )=y2i__21
         y(1,2,l2,l3,ii )=y2i__12
         y(2,2,l2,l3,ii )=y2i__22
         
         y(1,1,l2,l3,ii1)=y2i1_11
         y(2,1,l2,l3,ii1)=y2i1_21
         y(1,2,l2,l3,ii1)=y2i1_12
         y(2,2,l2,l3,ii1)=y2i1_22
      enddo
   enddo
enddo

END SUBROUTINE comb_rot_grow_1


!> In one dimesnion,    
!! with optimised cycles
!! Applies synthesis wavelet transformation 
!! then convolves with magic filter
!! the size of the data is allowed to grow
subroutine comb_rot_grow_2(n1,n2,nfl2,nfu2,nfl3,nfu3,x,y) !n(c) nfl1,nfu1 (arg:3,4)
use module_defs, only: wp
implicit none
integer,intent(in)::n1,n2
integer, intent(in) :: nfl2,nfu2,nfl3,nfu3 !n(c) nfl1,nfu1
real(wp), dimension(2,2,nfl2:nfu2,nfl3:nfu3,0:2*n1+1), intent(in) :: x
real(wp), dimension(2,nfl3:nfu3,0:2*n1+1,0:2*n2+1), intent(out) :: y
integer l1,l3,i,t !n(c) ,l2
integer :: ii,ii1
real(wp) y2i__1,y2i__2,y2i1_1,y2i1_2
integer::modul(-14+2*nfl2:2*nfu2+16)

include 'v_17.inc'
call fill_mod_arr(modul,-14+2*nfl2,2*nfu2+16,2*n2+2)

do l1=0,2*n1+1
   do l3=nfl3,nfu3
      
      do ii=0,2*nfl2-15
         y(:,l3,l1,ii)=0.d0   
      enddo         
      do ii=2*nfu2+17,2*n2+1
         y(:,l3,l1,ii)=0.d0   
      enddo
      
      ii=modul(2*nfu2+16)      
      y(1,l3,l1,ii)=fil2(16,1)*x(1,1,nfu2,l3,l1)+fil2(16,2)*x(2,1,nfu2,l3,l1)
      y(2,l3,l1,ii)=fil2(16,1)*x(1,2,nfu2,l3,l1)+fil2(16,2)*x(2,2,nfu2,l3,l1)
      
      ! scaling function everywhere, but folded
      do i=nfl2-7,nfu2+7 
         y2i__1=0.d0
         y2i__2=0.d0
         y2i1_1=0.d0
         y2i1_2=0.d0
         
         ! wavelet restricted
         do t=max(i-8,nfl2),min(i+7,nfu2)
            y2i__1=y2i__1+fil2(2*(i-t)  ,1)*x(1,1,t,l3,l1)+fil2(2*(i-t)  ,2)*x(2,1,t,l3,l1)
            y2i__2=y2i__2+fil2(2*(i-t)  ,1)*x(1,2,t,l3,l1)+fil2(2*(i-t)  ,2)*x(2,2,t,l3,l1)
            y2i1_1=y2i1_1+fil2(2*(i-t)+1,1)*x(1,1,t,l3,l1)+fil2(2*(i-t)+1,2)*x(2,1,t,l3,l1)
            y2i1_2=y2i1_2+fil2(2*(i-t)+1,1)*x(1,2,t,l3,l1)+fil2(2*(i-t)+1,2)*x(2,2,t,l3,l1)
         enddo
         
         ii =modul(2*i  )
         ii1=modul(2*i+1)
         
         y(1,l3,l1,ii )=y2i__1
         y(2,l3,l1,ii )=y2i__2
         y(1,l3,l1,ii1)=y2i1_1
         y(2,l3,l1,ii1)=y2i1_2
      enddo
   enddo
enddo

END SUBROUTINE comb_rot_grow_2



!> In one dimesnion,    
!! with optimised cycles
!! Applies synthesis wavelet transformation 
!! then convolves with magic filter
!! then adds the result to y.
!! The size of the data is allowed to grow
subroutine  comb_rot_grow_3(n1,n2,n3,nfl3,nfu3,x,y) !n(c) nfl1,nfu1,nfl2,nfu2
use module_defs, only: wp
implicit none
integer,intent(in) :: nfl3,nfu3,n1,n2,n3 !n(c) nfl1,nfu1,nfl2,nfu2
real(wp), dimension(2,nfl3:nfu3,0:2*n1+1,0:2*n2+1), intent(in) :: x
real(wp), dimension(0:2*n1+1,0:2*n2+1,0:2*n3+1), intent(out) :: y
!local variables
integer :: l1,l2,i,t !n(c) ,ll1,l1_0,l1_1
integer :: ii,ii1
real(wp) :: y2i,y2i1 !n(c) y2i__0,y2i__1,y2i1_0,y2i1_1
integer::modul(-14+2*nfl3:2*nfu3+16)
include 'v_17.inc'
call fill_mod_arr(modul,-14+2*nfl3,2*nfu3+16,2*n3+2)

do l2=0,2*n2+1
   do l1=0,2*n1+1
      ii=modul(2*nfu3+16)
      y(l1,l2,ii)=y(l1,l2,ii)+fil2(16,1)*x(1,nfu3,l1,l2)+fil2(16,2)*x(2,nfu3,l1,l2)
      
      do i=nfl3-7,nfu3+7 
         ii =modul(2*i  )
         ii1=modul(2*i+1)
         
         y2i =y(l1,l2,ii)
         y2i1=y(l1,l2,ii1)

         do t=max(i-8,nfl3),min(i+7,nfu3)
            y2i =y2i +fil2(2*(i-t)  ,1)*x(1,t,l1,l2)+fil2(2*(i-t)  ,2)*x(2,t,l1,l2)
            y2i1=y2i1+fil2(2*(i-t)+1,1)*x(1,t,l1,l2)+fil2(2*(i-t)+1,2)*x(2,t,l1,l2)
         enddo

         y(l1,l2,ii)=y2i
         y(l1,l2,ii1)=y2i1
      enddo
   enddo
enddo

END SUBROUTINE comb_rot_grow_3


!> In one dimesnion,    
!! with optimised cycles
!! Applies synthesis wavelet transformation 
!! then convolves with magic filter
!! then adds the result to y.
!! The size of the data is allowed to grow
subroutine  comb_rot_grow(n1,ndat,x,y)
use module_defs, only: wp
implicit none
integer,intent(in) :: n1,ndat
real(wp), dimension(0:n1,ndat), intent(in) :: x
real(wp), dimension(ndat,0:2*n1+1), intent(out) :: y
!local variables
integer :: l,i,t,tt
real(wp) :: y2i,y2i1
integer::modul(-8:n1+7)
include 'v_17.inc'
call fill_mod_arr(modul,-8,n1+7,n1+1)

do l=1,ndat
   do i=0,n1
      y2i =0._wp
      y2i1=0._wp
      do t=i-8,i+7
         tt=modul(t)
         y2i =y2i +fil2(2*(i-t)  ,1)*x(tt,l)
         y2i1=y2i1+fil2(2*(i-t)+1,1)*x(tt,l)
      enddo
      y(l,2*i  )=y2i
      y(l,2*i+1)=y2i1
   enddo
enddo

END SUBROUTINE comb_rot_grow


!> In 3d,            
!! Applies the magic filter transposed, then analysis wavelet transformation.
!! The size of the data is forced to shrink
!! The input array y is not overwritten
subroutine comb_shrink_hyb(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,y,xc,xf)
use module_defs, only: wp
implicit none
integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
real(wp), dimension(0:2*n1+1,0:2*n2+1,0:2*n3+1), intent(in) :: y
real(wp), dimension(max(2*(2*n2+2)*(2*n3+2)*(nfu1-nfl1+1),&
     (2*n2+2)*(2*n3+2)*(n1+1))), intent(inout) :: w1
real(wp), dimension(max(4*(2*n3+2)*(nfu1-nfl1+1)*(nfu2-nfl2+1),&
     (2*n3+2)*(n1+1)*(n2+1))), intent(inout) :: w2
real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: xc
real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: xf

integer nt

   !perform the combined transform    
   
   call comb_shrink_hyb_c(n1,n2,n3,w1,w2,y,xc)
     
   ! I1,I2,I3 -> I2,I3,i1
   call comb_rot_shrink_hyb_1(n1,n2,n3,nfl1,nfu1,y,w1)
   
   ! I2,I3,i1 -> I3,i1,i2
   nt=(2*n3+2)*(nfu1-nfl1+1)
   call comb_rot_shrink_hyb_2(nt,w1,w2,nfl2,nfu2,n2)
   
   ! I3,i1,i2 -> i1,i2,i3
   nt=(nfu1-nfl1+1)*(nfu2-nfl2+1)
   call comb_rot_shrink_hyb_3(nt,w2,xf,nfl3,nfu3,n3)
   
END SUBROUTINE comb_shrink_hyb


!> In one dimension,    
!! Applies the magic filter transposed, then analysis wavelet transformation.
!! The size of the data is forced to shrink
subroutine comb_rot_shrink_hyb_1(n1,n2,n3,nfl1,nfu1,x,y) !n(c) nfl2,nfu2,nfl3,nfu3 (arg:6,7,8,9)
use module_defs, only: wp
implicit none
integer, intent(in) :: n1,n2,n3,nfl1,nfu1 !n(c) nfl2,nfu2,nfl3,nfu3
real(wp), dimension(0:2*n1+1,0:2*n2+1,0:2*n3+1), intent(in) :: x
real(wp), dimension(2,0:2*n2+1,0:2*n3+1,nfl1:nfu1), intent(out) :: y
!local variables
integer, parameter :: lowfil2=-14,lupfil2=16
integer :: i,j2,j3,l,ll !n(c) nflop, icur
real(wp) :: ci1,ci2
integer::modul(lowfil2+2*nfl1:2*nfu1+lupfil2)
include 'v.inc'
call fill_mod_arr(modul,lowfil2+2*nfl1,2*nfu1+lupfil2,2*n1+2)

do j2=0,2*n2+1
   do j3=0,2*n3+1
      do i=nfl1,nfu1
         ci1=0.d0
         ci2=0.d0
         do l=lowfil2+2*i,lupfil2+2*i
            ! the input data are wrapped around because of periodic BC
            ll=modul(l)
            ci1=ci1+fil2(l-2*i,1)*x(ll,j2,j3)
            ci2=ci2+fil2(l-2*i,2)*x(ll,j2,j3)
         enddo
         y(1,j2,j3,i)=ci1
         y(2,j2,j3,i)=ci2
      enddo
   enddo
enddo

END SUBROUTINE comb_rot_shrink_hyb_1


!> In one dimension,    
!! Applies the magic filter transposed, then analysis wavelet transformation.
!! The size of the data is forced to shrink
subroutine comb_rot_shrink_hyb_2(ndat,x,y,nfl,nfu,n1)
use module_defs, only: wp
implicit none
integer, parameter:: lowfil2=-14,lupfil2=16
integer, intent(in) :: ndat,nfl,nfu,n1
real(wp), dimension(2,0:2*n1+1,ndat), intent(in) :: x
real(wp), dimension(2,2,ndat,nfl:nfu), intent(out) :: y
!local variables
integer :: j,i,l,ll !n(c) nflop, icur
real(wp) :: ci11,ci12,ci21,ci22
integer::modul(lowfil2+2*nfl:2*nfu+lupfil2)
include 'v.inc'
call fill_mod_arr(modul,lowfil2+2*nfl,2*nfu+lupfil2,2*n1+2)

do j=1,ndat
   do i=nfl,nfu
      ci11=0.d0
      ci12=0.d0
      ci21=0.d0
      ci22=0.d0
      
      do l=lowfil2+2*i,lupfil2+2*i
         ! the input data are wrapped around because of periodic BC
         ll=modul(l)
         ci11=ci11+fil2(l-2*i,1)*x(1,ll,j)
         ci12=ci12+fil2(l-2*i,2)*x(1,ll,j)
         ci21=ci21+fil2(l-2*i,1)*x(2,ll,j)
         ci22=ci22+fil2(l-2*i,2)*x(2,ll,j)
      enddo
      y(1,1,j,i)=ci11
      y(1,2,j,i)=ci12
      y(2,1,j,i)=ci21
      y(2,2,j,i)=ci22
   enddo
enddo

END SUBROUTINE comb_rot_shrink_hyb_2


!> In one dimension,    
!! Applies the magic filter transposed, then analysis wavelet transformation.
!! The size of the data is forced to shrink
subroutine comb_rot_shrink_hyb_3(ndat,x,y,nfl,nfu,n1)
use module_defs, only: wp
implicit none
integer, parameter :: lowfil2=-14,lupfil2=16
integer, intent(in) :: ndat,nfl,nfu,n1
real(wp), dimension(2,2,0:2*n1+1,ndat), intent(in) :: x
real(wp), dimension(7,ndat,nfl:nfu), intent(out) :: y
!local variables
integer :: i,j,l,ll
real(wp) :: ci112,ci121,ci122,ci211,ci212,ci221,ci222
integer::modul(lowfil2+2*nfl:2*nfu+lupfil2)
include 'v.inc'
call fill_mod_arr(modul,lowfil2+2*nfl,2*nfu+lupfil2,2*n1+2)

do j=1,ndat
   do i=nfl,nfu
      ci112=0._wp
      ci121=0._wp
      ci122=0._wp
      ci211=0._wp
      ci212=0._wp
      ci221=0._wp
      ci222=0._wp
      
      do l=lowfil2+2*i,lupfil2+2*i
         ! the input data are wrapped around because of periodic BC
         ll=modul(l)
         ci112=ci112+fil2(l-2*i,2)*x(1,1,ll,j)
         ci121=ci121+fil2(l-2*i,1)*x(1,2,ll,j)
         ci122=ci122+fil2(l-2*i,2)*x(1,2,ll,j)
         ci211=ci211+fil2(l-2*i,1)*x(2,1,ll,j)
         ci212=ci212+fil2(l-2*i,2)*x(2,1,ll,j)
         ci221=ci221+fil2(l-2*i,1)*x(2,2,ll,j)
         ci222=ci222+fil2(l-2*i,2)*x(2,2,ll,j)
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

END SUBROUTINE comb_rot_shrink_hyb_3


!> In one dimension,    
!! Applies the magic filter transposed, then analysis wavelet transformation.
!! The size of the data is forced to shrink
subroutine comb_rot_shrink_hyb(ndat,x,y,n1)
use module_defs, only: wp
implicit none
integer, parameter :: lowfil2=-14,lupfil2=16
integer, intent(in) :: ndat,n1
real(wp), dimension(0:2*n1+1,ndat), intent(in) :: x
real(wp), dimension(ndat,0:n1), intent(out) :: y
!local variables
integer :: j,l,i,ll !n(c) icur
real(wp) :: ci
integer::modul(lowfil2:2*n1+lupfil2)
include 'v.inc'
call fill_mod_arr(modul,lowfil2,2*n1+lupfil2,2*n1+2)

! the convolution itself:
do j=1,ndat
   do i=0,n1
      ci=0.d0
      do l=lowfil2+2*i,lupfil2+2*i
         ll=modul(l)
         ci=ci+fil2(l-2*i,1)*x(ll,j)
      enddo
      y(j,i)=ci
   enddo
enddo
   
END SUBROUTINE comb_rot_shrink_hyb
