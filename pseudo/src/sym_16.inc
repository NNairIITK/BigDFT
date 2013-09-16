!> @file
!! Include file for Daubehcies-16 filter
!! @author
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!!

!> Daubechy S6
integer, parameter :: m=8,MM=M+2
real(kind=8), dimension(-mM:mM) :: cg,cht,cgt
real(kind=8), parameter :: CH(-MM:MM) =  (/ 0.D0,0.D0,0.d0, &
        -0.0033824159510050025955D0,-0.00054213233180001068935D0, &
        0.031695087811525991431D0,0.0076074873249766081919D0, &
        -0.14329423835127266284D0,-0.061273359067811077843D0, &
        0.48135965125905339159D0,0.77718575169962802862D0,0.36444189483617893676D0, &
        -0.051945838107881800736D0,-0.027219029917103486322D0, &
        0.049137179673730286787D0,0.0038087520138944894631D0, &
        -0.014952258337062199118D0,-0.00030292051472413308126D0, &
        0.0018899503327676891843D0,0.D0,0.D0 /)
integer :: i
! coefficients for wavelet transform
do i=-mM,mM
   cht(i)=0.d0
   cg(i)=0.d0
   cgt(i)=0.d0
enddo


! orthonormal wavelet
do i=-mM,mM
   cht(i)=ch(i)
enddo

! g coefficients from h coefficients
do i=-mM,mM-1
   cg(i+1)=cht(-i)*(-1)**(i+1)
   cgt(i+1)=ch(-i)*(-1)**(i+1)
enddo
