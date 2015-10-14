!> @file
!!  Routines to apply the magic filter and an analysis wavelet transform
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> In 3d,            
!! Applies the magic filter transposed, then analysis wavelet transformation.
!! The size of the data is forced to shrink
!! The input array y is not overwritten
subroutine comb_shrink_hyb_c(n1,n2,n3,w1,w2,y,x)
use module_defs, only: wp
implicit none
integer,intent(in) ::n1,n2,n3 
real(wp),dimension(0:2*n1+1,0:2*n2+1,0:2*n3+1),intent(in)::y!input
real(wp),dimension         (0:2*n2+1,0:2*n3+1,0:n1), intent(inout)::w1
real(wp),dimension                  (0:2*n3+1,0:n1,0:n2), intent(inout)::w2
real(wp),dimension                           (0:n1,0:n2,0:n3), intent(out)::x!output
!local variables
integer :: nt

! I1,I2,I3 -> I2,I3,i1
nt=(2*n2+2)*(2*n3+2)
call comb_rot_shrink_hyb(nt,y,w1,n1)

! I2,I3,i1 -> I3,i1,i2
nt=(2*n3+2)*(n1+1)
call comb_rot_shrink_hyb(nt,w1,w2,n2)

! I3,i1,i2 -> i1,i2,i3
nt=(n1+1)*(n2+1)
call comb_rot_shrink_hyb(nt,w2,x,n3)

END SUBROUTINE comb_shrink_hyb_c


subroutine comb_grow_c_simple(n1,n2,n3,w1,w2,x,y)
  use module_defs, only: wp
  implicit none
  integer,intent(in)::n1,n2,n3
  real(wp),intent(in)::x(0:n1,0:n2,0:n3)   
  real(wp),dimension((2*n1+2)*(n2+1)*(n3+1)), intent(inout) :: w1 !work
  real(wp),dimension((n3+1)*(2*n1+2)*(2*n2+2)), intent(inout) :: w2 ! work
  real(wp),dimension((2*n1+2)*(2*n2+2)*(2*n3+2)), intent(out) :: y

  integer::nt

  ! i1,i2,i3 -> i2,i3,I1
  nt=(n2+1)*(n3+1)
  call  comb_rot_grow(n1,nt,x,w1) 

  ! i2,i3,I1 -> i3,I1,I2
  nt=(n3+1)*(2*n1+2)
  call  comb_rot_grow(n2,nt,w1,w2) 

  ! i3,I1,I2  -> I1,I2,I3
  nt=(2*n1+2)*(2*n2+2)
  call  comb_rot_grow(n3,nt,w2,y) 
  
END SUBROUTINE comb_grow_c_simple
