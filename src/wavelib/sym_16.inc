!> @file
!!   Daubechy S6
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

integer, parameter :: m=8,mm=m+2

real(wp) :: cg(-mm:mm),cht(-mm:mm),cgt(-mm:mm)

real(wp), parameter :: CH(-MM:MM) = (/ &
        0.0_wp,0.0_wp,0.0_wp, &
        -0.0033824159510050025955_wp,-0.00054213233180001068935_wp, &
        0.031695087811525991431_wp,0.0076074873249766081919_wp, &
        -0.14329423835127266284_wp,-0.061273359067811077843_wp, &
        0.48135965125905339159_wp,0.77718575169962802862_wp,0.36444189483617893676_wp, &
        -0.051945838107881800736_wp,-0.027219029917103486322_wp, &
        0.049137179673730286787_wp,0.0038087520138944894631_wp, &
        -0.014952258337062199118_wp,-0.00030292051472413308126_wp, &
        0.0018899503327676891843_wp,0.0_wp,0.0_wp &
        /)

! Coefficints for wavelet transform
do i=-mm,mm
   cht(i)=0.0_wp
   cg(i)=0.0_wp
   cgt(i)=0.0_wp
end do


! orthonormal wavelet
do i=-mm,mm
   cht(i)=ch(i)
end do

! g coefficients from h coefficients
do i=-mM,mM-1
   cg(i+1)=cht(-i)*(-1.0_wp)**(i+1)
   cgt(i+1)=ch(-i)*(-1.0_wp)**(i+1)
end do
