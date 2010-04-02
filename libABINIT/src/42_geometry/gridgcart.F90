!{\src2tex{textfont=tt}}
!!****f* ABINIT/gridgcart
!! NAME
!! gridgcart
!!
!! FUNCTION
!! Returns the reciprocal space vectors corresponding
!! to the real space grid
!!
!! COPYRIGHT
!! Copyright (C) 2007-2010 ABINIT group (JZwanziger)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  real(dp) :: gprimd(3,3) ! primitive translation vectors in G space, units of 1/a0
!!  integer :: ngfft(18) ! fft info
!!
!! OUTPUT
!!  real(dp) :: gcart(ngfft(1),ngfft(2),ngfft(3),3) ! cartesian components of G vectors on grid
!!
!! PARENTS        ! This paragraph will be produced automatically
!!      calc_efg,respfn
!!
!! NOTES
!! this subroutine returns a three component G vector at index ijk given
!! information about the real space grid contained in ngfft(1,2,3)
!!
!! The reason for this function is that in the calc_efg.F90 code, we need to
!! compute derivatives of real space functions in reciprocal space, so
!! we are calculating sums over G_i*G_j - delta(ij)G^2, for example, where
!! G_i, G_j are components of the G vector. For the real space grid
!! [1..ngfft(1)][1..ngfft(2)][1..ngfft(3)], the grid "spacing" is Rp(1)/ngfft(1),
!! Rp(2)/ngfft(2), and Rp(3)/ngfft(3) in the three directions, where Rp(i) is the
!! length of the i'th primitive translation vector. That means that in G space,
!! the components range from -Gp(i)*ngfft(i)/2..0..+Gp(i)*ngfft(i)/2, for i =1, 2, and 3,
!! and Gp(i) the length of the i'th primitive translation vector in G space.
!!
!! In FFT algorithms, the transformed data is ordered such that zero comes first, then
!! increasing positive frequencies, then *decreasing* negative frequencies. The cases of
!! even and odd data must be distinguished, as follows.
!! If ngfft(i) is even, then components ngfft(i)/2 and -ngfft(i)/2 are identical and
!! both held in data point ngfft(i)/2 + 1 in the output vector. If ngfft(i) is odd, they
!! are offset. Example:
!!
!! Even number of points:
!! pt 1 2 3  4    5  6
!! G  0 1 2 +/-3 -2 -1
!!
!! Odd number of points
!! pt 1 2 3 4  5  6  7
!! G  0 1 2 3 -3 -2 -1
!!
!! To assign the G values for each point, the positive G values are assigned to points 1 to ngfft(i)/2 + 1
!! This means that for odd ngfft(i), the Nyquist frequencies occurs at adjacent points ngfft(i)/2+1 and
!! ngfft(i)/2 + 2, while for even ngfft(i), only the positive Nyquist G component (not
!! the negative!).
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

subroutine gridgcart(gcart,gprimd,ngfft)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
! fft and grid info
! Grid of G vectors
! primitive translation vectors in G space, units of 1/a0
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(out) :: gcart(ngfft(1),ngfft(2),ngfft(3),3)

!Local variables-------------------------------
! gred is the "reduced" g vector components, they will be dimensioned with gprimd later
!scalars
 integer :: ii,jj,kk
!arrays
 real(dp) :: gred(3)

! *************************************************************************

 do ii = 1, ngfft(1) ! looping over first dimension of grid
   if ( ii <= (1 + ngfft(1)/2) ) then ! split components by positives in the first half
     gred(1) = ii - 1
   else
     gred(1) = -ngfft(1) + ii - 1 ! and negatives in the second half
   end if
   do jj = 1, ngfft(2) ! looping over 2nd dimension
     if ( jj <= (1 + ngfft(2)/2) ) then
       gred(2) = jj - 1
     else
       gred(2) = -ngfft(2) + jj - 1
     end if
     do kk = 1, ngfft(3) ! looping over third dimension
       if ( kk <= (1 + ngfft(3)/2) ) then
         gred(3) = kk - 1
       else
         gred(3) = -ngfft(3) + kk - 1
       end if
!      now the reduced g vector is converted to cartesian components
       gcart(ii,jj,kk,1) = gred(1)*gprimd(1,1)+gred(2)*gprimd(1,2)+gred(3)*gprimd(1,3)
       gcart(ii,jj,kk,2) = gred(1)*gprimd(2,1)+gred(2)*gprimd(2,2)+gred(3)*gprimd(2,3)
       gcart(ii,jj,kk,3) = gred(1)*gprimd(3,1)+gred(2)*gprimd(3,2)+gred(3)*gprimd(3,3)
     end do ! end loop over ngfft(3)
   end do ! end loop over ngfft(2)
 end do ! end loop over ngfft(1)
end subroutine gridgcart
!!***
