!> @file
!! Include file for Daubechies-8
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

integer, parameter :: m=4
real(kind=8) :: ch(-m:m),cg(-m:m),cht(-m:m),cgt(-m:m)

! Coefficients for wavelet transform
do i=-m,m
   ch(i)=0.d0
   cht(i)=0.d0
   cg(i)=0.d0
   cgt(i)=0.d0
end do

! Daubechy D8
ch(-3)= .2303778133088964d0
ch(-2)= .7148465705529154d0
ch(-1)= .6308807679298587d0
ch( 0)=-.0279837694168599d0
ch( 1)=-.1870348117190931d0
ch( 2)= .0308413818355607d0
ch( 3)= .0328830116668852d0
ch( 4)=-.0105974017850690d0

! orthonormal wavelet
do i=-m,m
   cht(i)=ch(i)
end do

! g coefficients from h coefficients
do i=-m,m-1
   cg(i+1)=cht(-i)*(-1.d0)**(i+1)
   cgt(i+1)=ch(-i)*(-1.d0)**(i+1)
end do
