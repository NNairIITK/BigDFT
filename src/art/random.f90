!> @file
!!   Routines to generate random numbers for ART methods
!! @author
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> Module used by ART methods (random routine)
MODULE random

  implicit none

  ! Shared variables
  integer :: idum, inext, inextp
  integer :: iff = 0
  integer, dimension(55) :: ma
end module random


!>  Random number generator (from "Numerical Recipes").
!!  Returns a uniform random deviate between 0.0 and 1.0.
!!  Set idum to any negative value to initialize or
!!  reinitialize the sequence.
real(8) function ran3()
  use RANDOM
  implicit none

  integer, parameter :: mbig=1000000000
  integer, parameter :: mseed=161803398
  integer, parameter :: mz=0
  real(8), parameter :: fac=1./mbig

  integer :: i,mj, mk, ii, k

  ! Any large mbig, and any smaller (but still large) mseed can be
  !  substituted for the above values.

  if(idum.lt.0.or.iff.eq.0)then
       iff=1
       mj=mseed-iabs(idum)
       mj=mod(mj,mbig)
       ma(55)=mj
       mk=1
       do i=1,54
         ii=mod(21*i,55)
         ma(ii)=mk
         mk=mj-mk
         if(mk.lt.mz)mk=mk+mbig
         mj=ma(ii)
       enddo
       do k=1,4
         do i=1,55
           ma(i)=ma(i)-ma(1+mod(i+30,55))
           if(ma(i).lt.mz)ma(i)=ma(i)+mbig
         enddo
       enddo
       inext=0
       inextp=31
       idum=1
  endif
  inext=inext+1
  if(inext.eq.56)inext=1
  inextp=inextp+1
  if(inextp.eq.56)inextp=1
  mj=ma(inext)-ma(inextp)
  if(mj.lt.mz)mj=mj+mbig
  ma(inext)=mj
  ran3=mj*fac
end function ran3
