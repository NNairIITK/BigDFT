!> @file
!!  Routines to apply synthesis and magic filters
!! @author
!!    Copyright (C) 2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!>   In one dimension,    
!!   with optimised cycles
!!   Applies synthesis wavelet transformation 
!!   then convolves with magic filter
!!   then adds the result to y.
!!   The size of the data is allowed to grow
subroutine  comb_rot_grow_ib_3(n1,n2,n3,nfl3,nfu3,x,y,ibxxyy)
   use module_defs, only: wp
   implicit none
   integer,intent(in) :: nfl3,nfu3,n1,n2,n3
   integer,intent(in)::ibxxyy(2,   0:2*n1+1,0:2*n2+1)
   real(wp), dimension(2,nfl3:nfu3,0:2*n1+1,0:2*n2+1), intent(in) :: x
   real(wp), dimension(            0:2*n1+1,0:2*n2+1,0:2*n3+1), intent(inout) :: y
   !local variables
   integer :: l1,l2,i,t
   integer :: ii,ii1
   real(wp) :: y2i,y2i1
   integer :: modul(-14+2*nfl3:2*nfu3+16)
   integer :: mfl3,mfu3
   include 'v_17.inc'
   call fill_mod_arr(modul,-14+2*nfl3,2*nfu3+16,2*n3+2)

!$omp parallel default(private)&
!$omp shared (n2,n1,ibxxyy,modul,x,y,fil2)
!$omp do
   do l2=0,2*n2+1
      do l1=0,2*n1+1
         mfl3=ibxxyy(1,l1,l2)
         mfu3=ibxxyy(2,l1,l2)

         if (mfl3.le.mfu3) then
         
            ii=modul(2*mfu3+16)
            y(l1,l2,ii)=y(l1,l2,ii)+fil2(16,1)*x(1,mfu3,l1,l2)+fil2(16,2)*x(2,mfu3,l1,l2)
            
            do i=mfl3-7,mfu3+7 
               ii =modul(2*i  )
               ii1=modul(2*i+1)
               
               y2i =y(l1,l2,ii)
               y2i1=y(l1,l2,ii1)
      
               do t=max(i-8,mfl3),min(i+7,mfu3)
                  y2i =y2i +fil2(2*(i-t)  ,1)*x(1,t,l1,l2)+fil2(2*(i-t)  ,2)*x(2,t,l1,l2)
                  y2i1=y2i1+fil2(2*(i-t)+1,1)*x(1,t,l1,l2)+fil2(2*(i-t)+1,2)*x(2,t,l1,l2)
               enddo
      
               y(l1,l2,ii)=y2i
               y(l1,l2,ii1)=y2i1
            enddo
         endif
      enddo
   enddo
!$omp enddo
!$omp end parallel

END SUBROUTINE comb_rot_grow_ib_3


!>   In one dimension,    
!!   with optimised cycles
!!   Applies synthesis wavelet transformation 
!!   then convolves with magic filter
!!   the size of the data is allowed to grow
subroutine comb_rot_grow_ib_2(n1,n2,nfl2,nfu2,nfl3,nfu3,x,y,ibzxx,ibxxyy)
use module_defs, only: wp
implicit none
integer,intent(in) :: n1,n2
integer, intent(in) :: nfl2,nfu2,nfl3,nfu3
integer,intent(in)::ibzxx(2,      nfl3:nfu3,0:2*n1+1)
integer,intent(in)::ibxxyy(2,               0:2*n1+1,0:2*n2+1)
real(wp), dimension(2,2,nfl2:nfu2,nfl3:nfu3,0:2*n1+1), intent(in) :: x
real(wp), dimension(2,            nfl3:nfu3,0:2*n1+1,0:2*n2+1), intent(out) :: y
integer :: l1,l3,i,t,l2
integer :: ii,ii1
real(wp) :: y2i__1,y2i__2,y2i1_1,y2i1_2
integer :: modul(-14+2*nfl2:2*nfu2+16)
integer :: mfl3,mfu3,mfl2,mfu2

include 'v_17.inc'
call fill_mod_arr(modul,-14+2*nfl2,2*nfu2+16,2*n2+2)

!$omp parallel default(private)&
!$omp shared (n2,n1,ibxxyy,modul,x,y,fil2,nfl3,nfu3,ibzxx)
!$omp do
do l2=0,2*n2+1
   do l1=0,2*n1+1
      mfl3=ibxxyy(1,l1,l2)
      mfu3=ibxxyy(2,l1,l2)
      y(:,mfl3:mfu3,l1,l2)=0._wp
   enddo
enddo
!$omp enddo
!$omp do
do l1=0,2*n1+1
   do l3=nfl3,nfu3
   
      mfl2=ibzxx(1,l3,l1)
      mfu2=ibzxx(2,l3,l1)
   
      if (mfl2.le.mfu2) then   
         ii=modul(2*mfu2+16)      
         y(1,l3,l1,ii)=fil2(16,1)*x(1,1,mfu2,l3,l1)+fil2(16,2)*x(2,1,mfu2,l3,l1)
         y(2,l3,l1,ii)=fil2(16,1)*x(1,2,mfu2,l3,l1)+fil2(16,2)*x(2,2,mfu2,l3,l1)
         
         ! scaling function everywhere, but folded
         do i=mfl2-7,mfu2+7 
            y2i__1=0._wp
            y2i__2=0._wp
            y2i1_1=0._wp
            y2i1_2=0._wp
            
            ! wavelet restricted
            do t=max(i-8,mfl2),min(i+7,mfu2)
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
      endif
   enddo
enddo
!$omp enddo
!$omp end parallel
END SUBROUTINE comb_rot_grow_ib_2


!>   In one dimension,    
!!   with optimised cycles
!!   Applies synthesis wavelet transformation 
!!   then convolves with magic filter
!!   the size of the data is allowed to grow
subroutine comb_rot_grow_ib_1(n1,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,x,y,ibyz,ibzxx)
use module_defs, only: wp
implicit none
integer, intent(in) :: n1,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
integer,intent(in)::ibyz(2     ,nfl2:nfu2,nfl3:nfu3)
integer,intent(in)::ibzxx(2,              nfl3:nfu3,0:2*n1+1)
real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: x
real(wp), dimension(2,2,        nfl2:nfu2,nfl3:nfu3,0:2*n1+1), intent(out) :: y
!local variables
integer :: l2,l3,i,t,l1
integer :: ii,ii1
real(wp) y2i__11,y2i__21,y2i1_11,y2i1_21
real(wp) y2i__12,y2i__22,y2i1_12,y2i1_22
integer :: modul(-14+2*nfl1:2*nfu1+16)
integer :: mfl1,mfu1,mfl2,mfu2

include 'v_17.inc'
call fill_mod_arr(modul,-14+2*nfl1,2*nfu1+16,2*n1+2)

!$omp parallel default(private)&
!$omp shared (n1,ibzxx,modul,x,y,fil2,nfl3,nfu3,nfl2,nfu2,ibyz)
!$omp do
do l1=0,2*n1+1
   do l3=nfl3,nfu3
      mfl2=ibzxx(1,l3,l1)
      mfu2=ibzxx(2,l3,l1)
      y(:,:,mfl2:mfu2,l3,l1)=0._wp
   enddo
enddo
!$omp enddo
!$omp do

do l3=nfl3,nfu3
   do l2=nfl2,nfu2
      mfl1=ibyz(1,l2,l3)
      mfu1=ibyz(2,l2,l3)

      if (mfl1.le.mfu1) then
   
         ii=modul(2*mfu1+16)
         y(1,1,l2,l3,ii)=fil2(16,2)*x(1,mfu1,l2,l3)
         y(2,1,l2,l3,ii)=fil2(16,1)*x(2,mfu1,l2,l3)+fil2(16,2)*x(3,mfu1,l2,l3)
         
         y(1,2,l2,l3,ii)=fil2(16,1)*x(4,mfu1,l2,l3)+fil2(16,2)*x(5,mfu1,l2,l3)
         y(2,2,l2,l3,ii)=fil2(16,1)*x(6,mfu1,l2,l3)+fil2(16,2)*x(7,mfu1,l2,l3)
      
         ! scaling function everywhere, but folded back to the box
         do i=mfl1-7,mfu1+7
            y2i__11=0._wp
            y2i__21=0._wp
            y2i__12=0._wp
            y2i__22=0._wp
            
            y2i1_11=0._wp
            y2i1_21=0._wp
            y2i1_12=0._wp
            y2i1_22=0._wp
         
            ! wavelet restricted   
            do t=max(i-8,mfl1),min(i+7,mfu1)
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
      endif
   enddo
enddo
!$omp enddo
!$omp end parallel
END SUBROUTINE comb_rot_grow_ib_1


!> In one dimension,    
!! with optimised cycles
!! Applies synthesis wavelet transformation 
!! then convolves with magic filter
!! then adds the result to y.
!! The size of the data is allowed to grow
subroutine comb_rot_grow(n1,ndat,x,y)
use module_defs, only: wp
implicit none
integer,intent(in) :: n1,ndat
real(wp), dimension(0:n1,ndat), intent(in) :: x
real(wp), dimension(ndat,0:2*n1+1), intent(out) :: y
!local variables
integer :: l,i,t,tt
real(wp) :: y2i,y2i1
real(wp) :: so1,so2,so3,so4,so5,so6,so7,so8
real(wp) :: se1,se2,se3,se4,se5,se6,se7,se8
integer::modul(-8:n1+7)
include 'v_17.inc'
call fill_mod_arr(modul,-8,n1+7,n1+1)

!$omp parallel default(private)&
!$omp shared (ndat,n1,modul,fil2,x,y)
!$omp do

do l=0,ndat/8-1
   do i=0,n1
      tt=modul(i-8)
      se1 =fil2(16,1)*x(tt,l*8+1)
      se2 =fil2(16,1)*x(tt,l*8+2)
      se3 =fil2(16,1)*x(tt,l*8+3)
      se4 =fil2(16,1)*x(tt,l*8+4)
      se5 =fil2(16,1)*x(tt,l*8+5)
      se6 =fil2(16,1)*x(tt,l*8+6)
      se7 =fil2(16,1)*x(tt,l*8+7)
      se8 =fil2(16,1)*x(tt,l*8+8)

      so1=0.d0
      so2=0.d0
      so3=0.d0
      so4=0.d0
      so5=0.d0
      so6=0.d0
      so7=0.d0
      so8=0.d0
      
      do t=i-7,i+7
         tt=modul(t)
         se1=se1 +fil2(2*(i-t)  ,1)*x(tt,l*8+1)
         se2=se2 +fil2(2*(i-t)  ,1)*x(tt,l*8+2)
         se3=se3 +fil2(2*(i-t)  ,1)*x(tt,l*8+3)
         se4=se4 +fil2(2*(i-t)  ,1)*x(tt,l*8+4)
         se5=se5 +fil2(2*(i-t)  ,1)*x(tt,l*8+5)
         se6=se6 +fil2(2*(i-t)  ,1)*x(tt,l*8+6)
         se7=se7 +fil2(2*(i-t)  ,1)*x(tt,l*8+7)
         se8=se8 +fil2(2*(i-t)  ,1)*x(tt,l*8+8)
         
         so1=so1+fil2(2*(i-t)+1,1)*x(tt,l*8+1)
         so2=so2+fil2(2*(i-t)+1,1)*x(tt,l*8+2)
         so3=so3+fil2(2*(i-t)+1,1)*x(tt,l*8+3)
         so4=so4+fil2(2*(i-t)+1,1)*x(tt,l*8+4)
         so5=so5+fil2(2*(i-t)+1,1)*x(tt,l*8+5)
         so6=so6+fil2(2*(i-t)+1,1)*x(tt,l*8+6)
         so7=so7+fil2(2*(i-t)+1,1)*x(tt,l*8+7)
         so8=so8+fil2(2*(i-t)+1,1)*x(tt,l*8+8)
      enddo
      y(l*8+1,2*i  )=se1
      y(l*8+2,2*i  )=se2
      y(l*8+3,2*i  )=se3
      y(l*8+4,2*i  )=se4
      y(l*8+5,2*i  )=se5
      y(l*8+6,2*i  )=se6
      y(l*8+7,2*i  )=se7
      y(l*8+8,2*i  )=se8
      
      y(l*8+1,2*i+1)=so1
      y(l*8+2,2*i+1)=so2
      y(l*8+3,2*i+1)=so3
      y(l*8+4,2*i+1)=so4
      y(l*8+5,2*i+1)=so5
      y(l*8+6,2*i+1)=so6
      y(l*8+7,2*i+1)=so7
      y(l*8+8,2*i+1)=so8
   enddo
enddo
!$omp enddo
!$omp do
do l=(ndat/8)*8+1,ndat
   do i=0,n1
      !SERIOUS DOUBTS ABOUT THIS LINE
      y2i =fil2(16,1)*x(tt,l)
      y2i1=0._wp
      do t=i-7,i+7
         tt=modul(t)
         y2i =y2i +fil2(2*(i-t)  ,1)*x(tt,l)
         y2i1=y2i1+fil2(2*(i-t)+1,1)*x(tt,l)
      enddo
      y(l,2*i  )=y2i
      y(l,2*i+1)=y2i1
   enddo
enddo
!$omp enddo
!$omp end parallel
END SUBROUTINE comb_rot_grow




!> In one dimension,    
!! Applies the magic filter transposed, then analysis wavelet transformation.
!! The size of the data is forced to shrink
subroutine comb_rot_shrink_hyb_1_ib(ndat,n1,nfl1,nfu1,x,y,ib)
use module_defs, only: wp
implicit none
integer, intent(in) :: ndat,nfl1,nfu1,n1
integer,intent(in)::ib(2,ndat)
real(wp), dimension(0:2*n1+1,ndat), intent(in) :: x
real(wp), dimension(2,ndat,nfl1:nfu1), intent(out) :: y
!local variables
integer, parameter :: lowfil2=-14,lupfil2=16
integer :: i,j,l,ll
real(wp) :: ci1,ci2
integer::modul(lowfil2+2*nfl1:2*nfu1+lupfil2)
include 'v.inc'
call fill_mod_arr(modul,lowfil2+2*nfl1,2*nfu1+lupfil2,2*n1+2)
!$omp parallel default(private)&
!$omp shared (ndat,ib,fil2,modul,x,y)
!$omp do

do j=1,ndat
   do i=ib(1,j),ib(2,j)
      ci1=0.d0
      ci2=0.d0
      do l=lowfil2+2*i,lupfil2+2*i
         ! the input data are wrapped around because of periodic BC
         ll=modul(l)
         ci1=ci1+fil2(l-2*i,1)*x(ll,j)
         ci2=ci2+fil2(l-2*i,2)*x(ll,j)
      enddo
      y(1,j,i)=ci1
      y(2,j,i)=ci2
   enddo
enddo
!$omp end do
!$omp end parallel
END SUBROUTINE comb_rot_shrink_hyb_1_ib


!> In one dimension,    
!! Applies the magic filter transposed, then analysis wavelet transformation.
!! The size of the data is forced to shrink
subroutine comb_rot_shrink_hyb_2_ib(ndat,x,y,nfl,nfu,n1,ib)
use module_defs, only: wp
implicit none
integer, parameter:: lowfil2=-14,lupfil2=16
integer, intent(in) :: ndat,nfl,nfu,n1
integer,intent(in)::ib(2,ndat)
real(wp), dimension(2,0:2*n1+1,ndat), intent(in) :: x
real(wp), dimension(2,2,ndat,nfl:nfu), intent(out) :: y
!local variables
integer :: j,i,l,ll
real(wp) :: ci11,ci12,ci21,ci22
integer::modul(lowfil2+2*nfl:2*nfu+lupfil2)
include 'v.inc'
call fill_mod_arr(modul,lowfil2+2*nfl,2*nfu+lupfil2,2*n1+2)
!$omp parallel default(private)&
!$omp shared (ndat,ib,fil2,modul,x,y)
!$omp do

do j=1,ndat
   do i=ib(1,j),ib(2,j)
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
!$omp enddo
!$omp end parallel

END SUBROUTINE comb_rot_shrink_hyb_2_ib


!> In one dimension,    
!! Applies the magic filter transposed, then analysis wavelet transformation.
!! The size of the data is forced to shrink
subroutine comb_rot_shrink_hyb_3_ib(ndat,x,y,nfl,nfu,n1,ib)
use module_defs, only: wp
implicit none
integer, parameter :: lowfil2=-14,lupfil2=16
integer, intent(in) :: ndat,nfl,nfu,n1
integer,intent(in)::ib(2,ndat)
real(wp), dimension(2,2,0:2*n1+1,ndat), intent(in) :: x
real(wp), dimension(7,ndat,nfl:nfu), intent(out) :: y
!local variables
integer :: i,j,l,ll
real(wp) :: ci112,ci121,ci122,ci211,ci212,ci221,ci222
integer::modul(lowfil2+2*nfl:2*nfu+lupfil2)
include 'v.inc'
call fill_mod_arr(modul,lowfil2+2*nfl,2*nfu+lupfil2,2*n1+2)

!$omp parallel default(private)&
!$omp shared (ndat,ib,fil2,modul,x,y)
!$omp do
do j=1,ndat
   do i=ib(1,j),ib(2,j)
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
!$omp enddo
!$omp end parallel
END SUBROUTINE comb_rot_shrink_hyb_3_ib


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
integer :: j,l,i,ll
real(wp) :: ci,ci1,ci2,ci3,ci4,ci5,ci6,ci7,ci8,ci9,ci10,ci11,ci12
integer::modul(lowfil2:2*n1+lupfil2)
include 'v.inc'
call fill_mod_arr(modul,lowfil2,2*n1+lupfil2,2*n1+2)
!$omp parallel default(private)&
!$omp shared (ndat,n1,modul,fil2,x,y)
!$omp do
! the convolution itself:
do j=0,ndat/12-1
   do i=0,n1
      ci1=0.d0
      ci2=0.d0
      ci3=0.d0
      ci4=0.d0
      ci5=0.d0
      ci6=0.d0
      ci7=0.d0
      ci8=0.d0
      ci9 =0.d0
      ci10=0.d0
      ci11=0.d0
      ci12=0.d0
      do l=lowfil2+2*i,lupfil2+2*i
         ll=modul(l)
         ci1=ci1+fil2(l-2*i,1)*x(ll,j*12+1)
         ci2=ci2+fil2(l-2*i,1)*x(ll,j*12+2)
         ci3=ci3+fil2(l-2*i,1)*x(ll,j*12+3)
         ci4=ci4+fil2(l-2*i,1)*x(ll,j*12+4)
         ci5=ci5+fil2(l-2*i,1)*x(ll,j*12+5)
         ci6=ci6+fil2(l-2*i,1)*x(ll,j*12+6)
         ci7=ci7+fil2(l-2*i,1)*x(ll,j*12+7)
         ci8=ci8+fil2(l-2*i,1)*x(ll,j*12+8)
         ci9 =ci9 +fil2(l-2*i,1)*x(ll,j*12+9 )
         ci10=ci10+fil2(l-2*i,1)*x(ll,j*12+10)
         ci11=ci11+fil2(l-2*i,1)*x(ll,j*12+11)
         ci12=ci12+fil2(l-2*i,1)*x(ll,j*12+12)
      enddo
      y(j*12+1,i)=ci1
      y(j*12+2,i)=ci2
      y(j*12+3,i)=ci3
      y(j*12+4,i)=ci4
      y(j*12+5,i)=ci5
      y(j*12+6,i)=ci6
      y(j*12+7,i)=ci7
      y(j*12+8,i)=ci8
      y(j*12+9 ,i)=ci9 
      y(j*12+10,i)=ci10
      y(j*12+11,i)=ci11
      y(j*12+12,i)=ci12
   enddo
enddo
!$omp enddo
!$omp do   
do j=(ndat/12)*12+1,ndat
   do i=0,n1
      ci=0.d0
      do l=lowfil2+2*i,lupfil2+2*i
         ll=modul(l)
         ci=ci+fil2(l-2*i,1)*x(ll,j)
      enddo
      y(j,i)=ci
   enddo
enddo
!$omp end do
!$omp end parallel
END SUBROUTINE comb_rot_shrink_hyb
