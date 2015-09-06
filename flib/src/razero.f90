!> @file
!!  Routines to initialize to zero arrays
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> set nbytes term to zero
subroutine setzero(nbytes,x)
  use f_precisions, only: f_long
  implicit none
  integer(f_long), intent(in) :: nbytes
  character, dimension(nbytes), intent(inout) :: x
  !local variables
  integer :: nthreads,ithread
  !$ integer omp_get_max_threads,omp_get_thread_num
  integer(f_long) :: nbt,it,nt,nb

  if (nbytes == 0_f_long) return
  nthreads=1
  ithread=0
  !$ nthreads=omp_get_max_threads()
  nt=int(nthreads,f_long)
  nbt=(nbytes+nt-1)/nt
  !$omp parallel private(ithread,it,nb) 
  !$ ithread=omp_get_thread_num()
  !calculate the number of elements for each thread
  it=min(nbt*ithread,nbytes-1)
  nb=min(nbytes-it,nbt)
  call memsetzero(x(it+1),nb)
  !$omp end parallel
end subroutine setzero

!> Routine initialize double precision arrays to zero
subroutine razero(n,x)
  implicit none
  !Arguments
  integer, intent(in) :: n
  double precision, dimension(n), intent(out) :: x
  !Local variables
  integer :: i
  !$ logical :: omp_in_parallel,do_omp
  !$ do_omp = n > 1024
  !$ if (do_omp) do_omp= .not. omp_in_parallel()
  !$omp parallel if (do_omp) shared(x,n) private(i)
  !$omp do
  do i=1,n
      x(i)=0.d0
  end do
  !$omp enddo
  !$omp end parallel

end subroutine razero

!> Set to zero an array x(n)
subroutine razero_simple(n,x)
  implicit none
  !Arguments
  integer, intent(in) :: n
  real(kind=4), intent(out) :: x(n)
  !Local variables
  integer :: i,m
  !$ logical :: omp_in_parallel,do_omp
  !$ do_omp = n > 1024
  !$ if (do_omp) do_omp= .not. omp_in_parallel()
  !$omp parallel if (do_omp) shared(x,n) private(i)
  !$omp do
  do i=1,n
      x(i)=0.e0
  end do
  !$omp enddo
  !$omp end parallel

  !!do i=1,n
  !!   x(i)=0.e0
  !!end do
  !m=mod(n,7)
  !if (m/=0) then
  !    do i=1,m
  !        x(i)=0.e0
  !    end do
  !    if (n<7) return
  !end if
  !m=m+1
  !do i=m,n,7
  !    x(i+0)=0.e0
  !    x(i+1)=0.e0
  !    x(i+2)=0.e0
  !    x(i+3)=0.e0
  !    x(i+4)=0.e0
  !    x(i+5)=0.e0
  !    x(i+6)=0.e0
  !end do
END SUBROUTINE razero_simple

!>   Set to zero an array x(n)
subroutine razero_integer(n,x)
  implicit none
  !Arguments
  integer, intent(in) :: n
  integer, dimension(n), intent(out) :: x
  !Local variables
  integer :: i,m
  !$ logical :: omp_in_parallel,do_omp
  !$ do_omp = n > 1024
  !$ if (do_omp) do_omp= .not. omp_in_parallel()
  !$omp parallel if (do_omp) shared(x,n) private(i)
  !$omp do
  do i=1,n
      x(i)=0
  end do
  !$omp enddo
  !$omp end parallel

  !!!!do i=1,n
  !!!!   x(i)=0
  !!!!end do
  !!m=mod(n,7)
  !!if (m/=0) then
  !!    do i=1,m
  !!        x(i)=0
  !!    end do
  !!    if (n<7) return
  !!end if
  !!m=m+1
  !!do i=m,n,7
  !!    x(i+0)=0
  !!    x(i+1)=0
  !!    x(i+2)=0
  !!    x(i+3)=0
  !!    x(i+4)=0
  !!    x(i+5)=0
  !!    x(i+6)=0
  !!end do
END SUBROUTINE razero_integer

!>   Set to zero an array x(n)
subroutine razero_integerlong(n,x)
  implicit none
  !Arguments
  integer, intent(in) :: n
  integer(kind=8), dimension(n), intent(out) :: x
  !Local variables
  integer :: i,m
  !$ logical :: omp_in_parallel,do_omp
  !$ do_omp = n > 1024
  !$ if (do_omp) do_omp= .not. omp_in_parallel()
  !$omp parallel if (do_omp) shared(x,n) private(i)
  !$omp do
  do i=1,n
      x(i)=int(0,kind=8)
  end do
  !$omp enddo
  !$omp end parallel


  !!!!do i=1,n
  !!!!   x(i)=0
  !!!!end do
  !!m=mod(n,7)
  !!if (m/=0) then
  !!    do i=1,m
  !!        x(i)=int(0,kind=8)
  !!    end do
  !!    if (n<7) return
  !!end if
  !!m=m+1
  !!do i=m,n,7
  !!    x(i+0)=int(0,kind=8)
  !!    x(i+1)=int(0,kind=8)
  !!    x(i+2)=int(0,kind=8)
  !!    x(i+3)=int(0,kind=8)
  !!    x(i+4)=int(0,kind=8)
  !!    x(i+5)=int(0,kind=8)
  !!    x(i+6)=int(0,kind=8)
  !!end do
END SUBROUTINE razero_integerlong

!!!>   Set to zero an array x(n): omp version of razero
!!subroutine omp_razero(n,x)
!!  use module_base
!!  implicit none
!!  !Arguments
!!  integer, intent(in) :: n
!!  real(kind=8), intent(out) :: x(n)
!!  !Local variables
!!  integer :: i,is
!!
!!
!!!!!$omp do
!!      do i=1,n-7,8
!!      x(i+0)=0.d0
!!      x(i+1)=0.d0
!!      x(i+2)=0.d0
!!      x(i+3)=0.d0
!!      x(i+4)=0.d0
!!      x(i+5)=0.d0
!!      x(i+6)=0.d0
!!      x(i+7)=0.d0
!!      x(i+8)=0.d0
!!      end do
!!!!!$omp enddo
!!      is=i
!!      do i=is,n
!!      x(i)=0.d0
!!      end do
!!END SUBROUTINE omp_razero


!>   Set to 10^-20 an array x(n) for exchange-correlation function of ABINIT
subroutine tenminustwenty(n,x,nproc)
  implicit none
! Arguments
  integer :: n,nproc
  real(kind=8) :: x(n)
! Local variables
  integer :: i
  do i=1,n
     x(i)=1.d-20/real(nproc,kind=8)
  end do
END SUBROUTINE tenminustwenty

!>   Set to 10^-10 an array x(n) for exchange-correlation function of ABINIT.
!!   We use 10^-10 here since the array will be squared later and we then arrive at
!!   the desired 10^-20.
subroutine tenminusten(n,x,nproc)
  implicit none
! Arguments
  integer :: n,nproc
  real(kind=8) :: x(n)
! Local variables
  integer :: i
  do i=1,n
     x(i)=1.d-10/real(nproc,kind=8)
  end do
END SUBROUTINE tenminusten


subroutine dasxpdy(n,da,dx,incx,dy,incy)
  implicit none
  integer, intent(in) :: n,incx,incy
  real(kind=8), intent(in) :: da
  real(kind=4), dimension(*), intent(in) :: dx
  real(kind=8), dimension(*), intent(inout) :: dy
  !local variables
  integer :: i,ix,iy
  
  ix=1
  iy=1
  do i=1,n
     dy(iy)=dy(iy)+da*real(dx(ix),kind=8)
     ix=ix+incx
     iy=iy+incy
  end do
end subroutine dasxpdy

subroutine dscopy(n,dx,incx,dy,incy)
  implicit none
  integer, intent(in) :: n,incx,incy
  real(kind=8), dimension(*), intent(in) :: dx
  real(kind=4), dimension(*), intent(out) :: dy
  !local variables
  integer :: i,ix,iy
  
  ix=1
  iy=1
  do i=1,n
     dy(iy)=real(dx(ix),kind=4)
     ix=ix+incx
     iy=iy+incy
  end do

end subroutine dscopy

subroutine icopy(n,dx,incx,dy,incy)
  implicit none
  integer, intent(in) :: n,incx,incy
  integer, dimension(*), intent(in) :: dx
  integer, dimension(*), intent(out) :: dy
  !local variables
  integer :: i,ix,iy
  
  ix=1
  iy=1
  do i=1,n
     dy(iy)=dx(ix)
     ix=ix+incx
     iy=iy+incy
  end do

end subroutine icopy


!> Module used in the builtin_rand function (see razero.f90)
module randomData
  implicit none

  integer, parameter :: ntab=32

  logical :: start = .true.
  integer :: iy = 0
  integer, dimension(NTAB) :: iv
end module randomData


!> Random Number generator from Numerical Recipes
!! To be used for reproducibility of the results
function builtin_rand(idum)
  use randomData, only : ntab, iy, iv, start

  implicit none

  integer, intent(inout) :: idum
  real(kind=4) :: builtin_rand
  !local variables
  integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836,ndiv=1+(im-1)/ntab
  real(kind=4), parameter :: am=1.e0/im,eps=1.2e-7,rnmx=1.-eps
  integer :: j,k

  if (start) then
     iv(:) = 0
     start = .false.
  end if
  if (idum <= 0.or. iy == 0) then
     idum=max(-idum,1)
     do j=ntab+8,1,-1
        k=idum/iq
        idum=ia*(idum-k*iq)-ir*k
        if (idum < 0) idum=idum+im
        if (j <= ntab) iv(j)=idum
     end do
     iy=iv(1)
  endif
  k=idum/iq
  idum=ia*(idum-k*iq)-ir*k
  if (idum <= 0) idum=idum+im
  j=1+iy/ndiv
  iy=iv(j)
  iv(j)=idum
  builtin_rand=min(am*iy,rnmx)
END FUNCTION builtin_rand

subroutine diff_i(n,a,b,diff)
  implicit none
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: a
  integer, dimension(n), intent(in) :: b
  integer, intent(out) :: diff
  !local variables
  integer :: i

  diff=0
  do i=1,n
     diff=max(diff,abs(a(i)-b(i)))
  end do
end subroutine diff_i
subroutine diff_li(n,a,b,diff)
  implicit none
  integer, intent(in) :: n
  integer(kind=8), dimension(n), intent(in) :: a
  integer(kind=8), dimension(n), intent(in) :: b
  integer(kind=8), intent(out) :: diff
  !local variables
  integer :: i

  diff=int(0,kind=8)
  do i=1,n
     diff=max(diff,abs(a(i)-b(i)))
  end do
end subroutine diff_li
subroutine diff_r(n,a,b,diff)
  implicit none
  integer, intent(in) :: n
  real, dimension(n), intent(in) :: a
  real, dimension(n), intent(in) :: b
  real, intent(out) :: diff
  !local variables
  integer :: i

  diff=0.0e0
  do i=1,n
     diff=max(diff,abs(a(i)-b(i)))
  end do
end subroutine diff_r
subroutine diff_d(n,a,b,diff)
  implicit none
  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: a
  double precision, dimension(n), intent(in) :: b
  double precision, intent(out) :: diff
  !local variables
  integer :: i

  diff=0.0d0
  do i=1,n
     diff=max(diff,abs(a(i)-b(i)))
  end do
end subroutine diff_d
subroutine diff_l(n,a,b,diff)
  implicit none
  integer, intent(in) :: n
  logical, dimension(n), intent(in) :: a
  logical, dimension(n), intent(in) :: b
  logical, intent(out) :: diff
  !local variables
  integer :: i

  diff=.false.
  do i=1,n
     diff=a(i) .eqv. b(i)
     if (diff) exit
  end do
end subroutine diff_l

subroutine diff_ci(n,a,b,diff)
  implicit none
  integer, intent(in) :: n
  character, dimension(n), intent(in) :: a
  integer, dimension(n), intent(in) :: b
  integer, intent(out) :: diff
  !local variables
  integer :: i

  diff=0
  do i=1,n
     diff=max(diff,abs(ichar(a(i))-b(i)))
  end do
end subroutine diff_ci

subroutine f_itoa(n,src,dest)
  implicit none
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: src
  character, dimension(n), intent(out) :: dest
  !local variables
  integer :: i
  
  do i=1,n
     dest(i)=achar(src(i))
  end do

end subroutine f_itoa

subroutine f_atoi(n,src,dest)
  implicit none
  integer, intent(in) :: n
  character, dimension(n), intent(in) :: src
  integer, dimension(n), intent(out) :: dest
  !local variables
  integer :: i
  
  do i=1,n
     dest(i)=ichar(src(i))
  end do

end subroutine f_atoi
