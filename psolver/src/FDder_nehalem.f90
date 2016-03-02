  !> @file
  !!   Module containing finite difference derivative operations
  !!
  !! @author
  !!    G. Fisicaro, L. Genovese (September 2015)
  !!    Copyright (C) 2002-2015 BigDFT group 
  !!    This file is distributed under the terms of the
  !!    GNU General Public License, see ~/COPYING file
  !!    or http://www.gnu.org/copyleft/gpl.txt .
  !!    For the list of contributors, see ~/AUTHORS 
module FDder
  use PSbase
  implicit none
  private

  public :: nabla_u_square,update_rhopol
  public :: nabla_u,div_u_i,nabla_u_and_square,nonvacuum_projection

  contains

include 'FDder_nehalem_kernels.f90'

SUBROUTINE div_u_i(geocode, n01, n02, n03, u, du, nord, hgrids, cc)
  integer, parameter :: wp=kind(1.0d0)
  character(len=1), intent(in) :: geocode !< @copydoc
  integer(kind=4), intent(in) :: n01
  integer(kind=4), intent(in) :: n02
  integer(kind=4), intent(in) :: n03
  real(kind=8), intent(in), dimension(0:n01 - (1), 0:n02 - (1), 0:n03 - (1), 0:2) :: u
  real(kind=8), intent(out), dimension(0:n01 - (1), 0:n02 - (1), 0:n03 - (1)) :: du
  integer(kind=4), intent(in) :: nord
  real(kind=8), intent(in), dimension(0:2) :: hgrids
  real(kind=8), intent(out), optional, dimension(0:n01 - (1), 0:n02 - (1), 0:n03 - (1)) :: cc
  integer(kind=4), dimension(3) :: nn
  integer(kind=4) :: bc0
  integer(kind=4) :: bc1
  integer(kind=4) :: bc2
  real(kind=8) :: a0
  real(kind=8) :: a1
  real(kind=8) :: a2
  real(kind=8), allocatable, dimension(:, :, :, :) :: tmp
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: w
  real(kind=8) :: u0
  real(kind=8) :: u1
  real(kind=8) :: u2
  nn(1) = n01
  nn(2) = n02
  nn(3) = n03
  u0 = 0.0_wp
  u1 = 0.0_wp
  u2 = 0.0_wp
  bc0 = 0
  bc1 = 0
  bc2 = 0
  a0 = (1.0) / (hgrids(0))
  a1 = (1.0) / (hgrids(1))
  a2 = (1.0) / (hgrids(2))
  if (geocode == 'F') then
    bc0 = -2
    bc2 = -2
  end if
  if (geocode /= 'P') then
    bc1 = -2
  end if
  if (present(cc)) then
    allocate( tmp(0:n01 - (1), 0:n02 - (1), 0:n03 - (1), 0:4) )
    call Poisson_broker(nord, 0, nn, bc0,  u(0, 0, 0, 0), du, a0)
    call Poisson_broker(nord, 1, nn, bc1,  u(0, 0, 0, 1), tmp(0, 0, 0, 0), a1)
    call Poisson_broker(nord, 2, nn, bc2,  u(0, 0, 0, 2), tmp(0, 0, 0, 1), a2)
    call Poisson_broker(nord, 0, nn, bc0,  u(0, 0, 0, 1), tmp(0, 0, 0, 2), a0)
    call Poisson_broker(nord, 0, nn, bc0,  u(0, 0, 0, 2), tmp(0, 0, 0, 3), a0)
    call Poisson_broker(nord, 1, nn, bc1,  u(0, 0, 0, 2), tmp(0, 0, 0, 4), a1)
!$omp parallel  default(shared) private(i1, i2, i3, u0, u1, u2)
!$omp do 
    do i3 = 0, n03 - (1), 1
      do i2 = 0, n02 - (1), 1
        do i1 = 0, n01 - (1), 5
          u0 = u(i1 + 0, i2, i3, 0)
          u1 = u(i1 + 0, i2, i3, 1)
          u2 = u(i1 + 0, i2, i3, 2)
          cc(i1 + 0, i2, i3) = ((u0) * (u0)) * (du(i1 + 0, i2, i3)) + &
((u1) * (u1)) * (tmp(i1 + 0, i2, i3, 0)) + ((u2) * (u2)) * (tmp(i1 + 0, i2, i3, 1))&
 + (2.0) * (((u0) * (u1)) * (tmp(i1 + 0, i2, i3, 2)) + &
((u0) * (u2)) * (tmp(i1 + 0, i2, i3, 3)) + ((u1) * (u2)) * (tmp(i1 + 0, i2, i3, 4)))
          du(i1 + 0, i2, i3) = du(i1 + 0, i2, i3) + tmp(i1 + 0, i2, i3, 0) + tmp(i1 + 0, i2, i3, 1)
          u0 = u(i1 + 1, i2, i3, 0)
          u1 = u(i1 + 1, i2, i3, 1)
          u2 = u(i1 + 1, i2, i3, 2)
          cc(i1 + 1, i2, i3) = ((u0) * (u0)) * (du(i1 + 1, i2, i3)) + &
((u1) * (u1)) * (tmp(i1 + 1, i2, i3, 0)) + ((u2) * (u2)) * (tmp(i1 + 1, i2, i3, 1))&
 + (2.0) * (((u0) * (u1)) * (tmp(i1 + 1, i2, i3, 2)) + &
((u0) * (u2)) * (tmp(i1 + 1, i2, i3, 3)) + ((u1) * (u2)) * (tmp(i1 + 1, i2, i3, 4)))
          du(i1 + 1, i2, i3) = du(i1 + 1, i2, i3) + tmp(i1 + 1, i2, i3, 0) + tmp(i1 + 1, i2, i3, 1)
          u0 = u(i1 + 2, i2, i3, 0)
          u1 = u(i1 + 2, i2, i3, 1)
          u2 = u(i1 + 2, i2, i3, 2)
          cc(i1 + 2, i2, i3) = ((u0) * (u0)) * (du(i1 + 2, i2, i3)) + &
((u1) * (u1)) * (tmp(i1 + 2, i2, i3, 0)) + ((u2) * (u2)) * (tmp(i1 + 2, i2, i3, 1))&
 + (2.0) * (((u0) * (u1)) * (tmp(i1 + 2, i2, i3, 2)) + &
((u0) * (u2)) * (tmp(i1 + 2, i2, i3, 3)) + ((u1) * (u2)) * (tmp(i1 + 2, i2, i3, 4)))
          du(i1 + 2, i2, i3) = du(i1 + 2, i2, i3) + tmp(i1 + 2, i2, i3, 0) + tmp(i1 + 2, i2, i3, 1)
          u0 = u(i1 + 3, i2, i3, 0)
          u1 = u(i1 + 3, i2, i3, 1)
          u2 = u(i1 + 3, i2, i3, 2)
          cc(i1 + 3, i2, i3) = ((u0) * (u0)) * (du(i1 + 3, i2, i3)) + &
((u1) * (u1)) * (tmp(i1 + 3, i2, i3, 0)) + ((u2) * (u2)) * (tmp(i1 + 3, i2, i3, 1))&
 + (2.0) * (((u0) * (u1)) * (tmp(i1 + 3, i2, i3, 2)) +&
 ((u0) * (u2)) * (tmp(i1 + 3, i2, i3, 3)) + ((u1) * (u2)) * (tmp(i1 + 3, i2, i3, 4)))
          du(i1 + 3, i2, i3) = du(i1 + 3, i2, i3) + tmp(i1 + 3, i2, i3, 0) + tmp(i1 + 3, i2, i3, 1)
          u0 = u(i1 + 4, i2, i3, 0)
          u1 = u(i1 + 4, i2, i3, 1)
          u2 = u(i1 + 4, i2, i3, 2)
          cc(i1 + 4, i2, i3) = ((u0) * (u0)) * (du(i1 + 4, i2, i3)) + &
((u1) * (u1)) * (tmp(i1 + 4, i2, i3, 0)) + ((u2) * (u2)) * (tmp(i1 + 4, i2, i3, 1))&
 + (2.0) * (((u0) * (u1)) * (tmp(i1 + 4, i2, i3, 2)) + &
((u0) * (u2)) * (tmp(i1 + 4, i2, i3, 3)) + ((u1) * (u2)) * (tmp(i1 + 4, i2, i3, 4)))
          du(i1 + 4, i2, i3) = du(i1 + 4, i2, i3) + tmp(i1 + 4, i2, i3, 0) + tmp(i1 + 4, i2, i3, 1)
        end do
      end do
    end do
!$omp end do 
!$omp end parallel 
    deallocate( tmp )
  else
    allocate( tmp(0:n01 - (1), 0:n02 - (1), 0:n03 - (1), 0:2) )
    call Poisson_broker(nord, 0, nn, bc0,  u(0, 0, 0, 0), tmp(0, 0, 0, 0), a0)
    call Poisson_broker(nord, 1, nn, bc1,  u(0, 0, 0, 1), tmp(0, 0, 0, 1), a1)
    call Poisson_broker(nord, 2, nn, bc2,  u(0, 0, 0, 2), tmp(0, 0, 0, 2), a2)
!$omp parallel  default(shared) private(i1, i2, i3)
!$omp do 
    do i3 = 0, n03 - (1), 1
      do i2 = 0, n02 - (1), 1
        do i1 = 0, n01 - (1), 1
          du(i1, i2, i3) = tmp(i1, i2, i3, 0) + tmp(i1, i2, i3, 1) + tmp(i1, i2, i3, 2)
        end do
      end do
    end do
!$omp end do 
!$omp end parallel 
  end if
END SUBROUTINE div_u_i
SUBROUTINE nabla_u_and_square(geocode, n01, n02, n03, u, du, ddu, nord, hgrids)
  integer, parameter :: wp=kind(1.0d0)
  character(len=1), intent(in) :: geocode !< @copydoc
  integer(kind=4), intent(in) :: n01
  integer(kind=4), intent(in) :: n02
  integer(kind=4), intent(in) :: n03
  real(kind=8), intent(in), dimension(0:n01 - (1), 0:n02 - (1), 0:n03 - (1)) :: u
  real(kind=8), intent(out), dimension(0:n01 - (1), 0:n02 - (1), 0:n03 - (1), 0:2) :: du
  real(kind=8), intent(out), dimension(0:n01 - (1), 0:n02 - (1), 0:n03 - (1)) :: ddu
  integer(kind=4), intent(in) :: nord
  real(kind=8), intent(in), dimension(0:2) :: hgrids
  integer(kind=4), dimension(3) :: nn
  integer(kind=4) :: bc0
  integer(kind=4) :: bc1
  integer(kind=4) :: bc2
  real(kind=8) :: a0
  real(kind=8) :: a1
  real(kind=8) :: a2
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  nn(1) = n01
  nn(2) = n02
  nn(3) = n03
  i1 = 0
  i2 = 0
  i3 = 0
  bc0 = 0
  bc1 = 0
  bc2 = 0
  if (geocode == 'F') then
    bc0 = -2
    bc2 = -2
  end if
  if (geocode /= 'P') then
    bc1 = -2
  end if
  a0 = (1.0) / (hgrids(0))
  a1 = (1.0) / (hgrids(1))
  a2 = (1.0) / (hgrids(2))
  call Poisson_broker(nord, 0, nn, bc0, u,  du(0, 0, 0, 0), a0)
  call Poisson_broker(nord, 1, nn, bc1, u,  du(0, 0, 0, 1), a1)
  call Poisson_broker(nord, 2, nn, bc2, u,  du(0, 0, 0, 2), a2)
!$omp parallel  default(shared) private(i1, i2, i3)
!$omp do 
  do i3 = 0, n03 - (1), 1
    do i2 = 0, n02 - (1), 1
      do i1 = 0, n01 - (1), 1
        ddu(i1, i2, i3) = (du(i1, i2, i3, 0)) * (du(i1, i2, i3, 0)) + &
(du(i1, i2, i3, 1)) * (du(i1, i2, i3, 1)) + &
(du(i1, i2, i3, 2)) * (du(i1, i2, i3, 2))
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE nabla_u_and_square
SUBROUTINE update_rhopol(geocode, n01, n02, n03, u, nord, hgrids, eta, dlogeps, rhopol, rhores2)
  integer, parameter :: wp=kind(1.0d0)
  character(len=1), intent(in) :: geocode !< @copydoc
  integer(kind=4), intent(in) :: n01
  integer(kind=4), intent(in) :: n02
  integer(kind=4), intent(in) :: n03
  real(kind=8), intent(in), dimension(0:n01 - (1), 0:n02 - (1), 0:n03 - (1)) :: u
  integer(kind=4) :: nord
  real(kind=8), intent(in), dimension(0:2) :: hgrids
  real(kind=8), intent(in) :: eta
  real(kind=8), intent(in), dimension(3, 0:n01 - (1), 0:n02 - (1), 0:n03 - (1)) :: dlogeps
  real(kind=8), intent(inout), dimension(0:n01 - (1), 0:n02 - (1), 0:n03 - (1)) :: rhopol
  real(kind=8), intent(out) :: rhores2
  integer(kind=4), dimension(3) :: nn
  integer(kind=4) :: bc0
  integer(kind=4) :: bc1
  integer(kind=4) :: bc2
  real(kind=8) :: a0
  real(kind=8) :: a1
  real(kind=8) :: a2
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  real(kind=8) :: res
  real(kind=8) :: rho
  real(kind=8), allocatable, dimension(:, :, :, :) :: du
  real(kind=8) :: oneo4pi
  real(kind=8) :: tmp_rhores2
  nn(1) = n01
  nn(2) = n02
  nn(3) = n03
  i1 = 0
  i2 = 0
  i3 = 0
  oneo4pi = 0.07957747154594767_wp
  bc0 = 0
  bc1 = 0
  bc2 = 0
  if (geocode == 'F') then
    bc0 = -2
    bc2 = -2
  end if
  if (geocode /= 'P') then
    bc1 = -2
  end if
  a0 = (1.0) / (hgrids(0))
  a1 = (1.0) / (hgrids(1))
  a2 = (1.0) / (hgrids(2))
  allocate( du(0:n01 - (1), 0:n02 - (1), 0:n03 - (1), 0:2) )
  call Poisson_broker(nord, 0, nn, bc0, u,  du(0, 0, 0, 0), a0)
  call Poisson_broker(nord, 1, nn, bc1, u,  du(0, 0, 0, 1), a1)
  call Poisson_broker(nord, 2, nn, bc2, u,  du(0, 0, 0, 2), a2)
  tmp_rhores2 = 0.0_wp
!$omp parallel  default(shared) private(i1, i2, i3, res, rho) reduction(+: tmp_rhores2)
!$omp do 
  do i3 = 0, n03 - (1), 1
    do i2 = 0, n02 - (1), 1
      do i1 = 0, n01 - (1), 1
        res = (dlogeps(1, i1, i2, i3)) * (du(i1, i2, i3, 0)) + &
(dlogeps(2, i1, i2, i3)) * (du(i1, i2, i3, 1)) + &
(dlogeps(3, i1, i2, i3)) * (du(i1, i2, i3, 2))
        res = (res) * (oneo4pi)
        rho = rhopol(i1, i2, i3)
        res = (res - (rho)) * (eta)
        tmp_rhores2 = tmp_rhores2 + (res) * (res)
        rhopol(i1, i2, i3) = res + rho
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
  deallocate( du )
  rhores2 = tmp_rhores2
END SUBROUTINE update_rhopol
SUBROUTINE nabla_u(geocode, n01, n02, n03, u, du, nord, hgrids)
  integer, parameter :: wp=kind(1.0d0)
  character(len=1), intent(in) :: geocode !< @copydoc
  integer(kind=4), intent(in) :: n01
  integer(kind=4), intent(in) :: n02
  integer(kind=4), intent(in) :: n03
  real(kind=8), intent(in), dimension(0:n01 - (1), 0:n02 - (1), 0:n03 - (1)) :: u
  real(kind=8), intent(out), dimension(0:n01 - (1), 0:n02 - (1), 0:n03 - (1), 0:2) :: du
  integer(kind=4), intent(in) :: nord
  real(kind=8), intent(in), dimension(0:2) :: hgrids
  integer(kind=4), dimension(3) :: nn
  integer(kind=4) :: bc0
  integer(kind=4) :: bc1
  integer(kind=4) :: bc2
  real(kind=8) :: a0
  real(kind=8) :: a1
  real(kind=8) :: a2
  nn(1) = n01
  nn(2) = n02
  nn(3) = n03
  bc0 = 0
  bc1 = 0
  bc2 = 0
  if (geocode == 'F') then
    bc0 = -2
    bc2 = -2
  end if
  if (geocode /= 'P') then
    bc1 = -2
  end if
  a0 = (1.0) / (hgrids(0))
  a1 = (1.0) / (hgrids(1))
  a2 = (1.0) / (hgrids(2))
  call Poisson_broker(nord, 0, nn, bc0, u,  du(0, 0, 0, 0), a0)
  call Poisson_broker(nord, 1, nn, bc1, u,  du(0, 0, 0, 1), a1)
  call Poisson_broker(nord, 2, nn, bc2, u,  du(0, 0, 0, 2), a2)
END SUBROUTINE nabla_u
SUBROUTINE nabla_u_epsilon(geocode, n01, n02, n03, u, du, nord, hgrids, eps)
  integer, parameter :: wp=kind(1.0d0)
  character(len=1), intent(in) :: geocode !< @copydoc
  integer(kind=4), intent(in) :: n01
  integer(kind=4), intent(in) :: n02
  integer(kind=4), intent(in) :: n03
  real(kind=8), intent(in), dimension(0:n01 - (1), 0:n02 - (1), 0:n03 - (1)) :: u
  real(kind=8), intent(out), dimension(0:n01 - (1), 0:n02 - (1), 0:n03 - (1), 0:2) :: du
  integer(kind=4), intent(in) :: nord
  real(kind=8), intent(in), dimension(0:2) :: hgrids
  real(kind=8), intent(in), dimension(0:n01 - (1), 0:n02 - (1), 0:n03 - (1)) :: eps
  integer(kind=4), dimension(3) :: nn
  integer(kind=4) :: bc0
  integer(kind=4) :: bc1
  integer(kind=4) :: bc2
  real(kind=8) :: a0
  real(kind=8) :: a1
  real(kind=8) :: a2
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  nn(1) = n01
  nn(2) = n02
  nn(3) = n03
  i1 = 0
  i2 = 0
  i3 = 0
  bc0 = 0
  bc1 = 0
  bc2 = 0
  if (geocode == 'F') then
    bc0 = -2
    bc2 = -2
  end if
  if (geocode /= 'P') then
    bc1 = -2
  end if
  a0 = (1.0) / (hgrids(0))
  a1 = (1.0) / (hgrids(1))
  a2 = (1.0) / (hgrids(2))
  call Poisson_broker(nord, 0, nn, bc0, u,  du(0, 0, 0, 0), a0)
  call Poisson_broker(nord, 1, nn, bc1, u,  du(0, 0, 0, 1), a1)
  call Poisson_broker(nord, 2, nn, bc2, u,  du(0, 0, 0, 2), a2)
!$omp parallel  default(shared) private(i1, i2, i3)
!$omp do 
  do i3 = 0, n03 - (1), 1
    do i2 = 0, n02 - (1), 1
      do i1 = 0, n01 - (1), 1
        du(i1, i2, i3, 0) = (du(i1, i2, i3, 0)) * (eps(i1, i2, i3))
        du(i1, i2, i3, 1) = (du(i1, i2, i3, 1)) * (eps(i1, i2, i3))
        du(i1, i2, i3, 2) = (du(i1, i2, i3, 2)) * (eps(i1, i2, i3))
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE nabla_u_epsilon
SUBROUTINE nabla_u_square(geocode, n01, n02, n03, u, ddu, nord, hgrids)
  integer, parameter :: wp=kind(1.0d0)
  character(len=1), intent(in) :: geocode !< @copydoc
  integer(kind=4), intent(in) :: n01
  integer(kind=4), intent(in) :: n02
  integer(kind=4), intent(in) :: n03
  real(kind=8), intent(in), dimension(0:n01 - (1), 0:n02 - (1), 0:n03 - (1)) :: u
  real(kind=8), intent(out), dimension(0:n01 - (1), 0:n02 - (1), 0:n03 - (1)) :: ddu
  integer(kind=4), intent(in) :: nord
  real(kind=8), intent(in), dimension(0:2) :: hgrids
  integer(kind=4), dimension(3) :: nn
  integer(kind=4) :: bc0
  integer(kind=4) :: bc1
  integer(kind=4) :: bc2
  real(kind=8) :: a0
  real(kind=8) :: a1
  real(kind=8) :: a2
  real(kind=8), allocatable, dimension(:, :, :, :) :: du
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  nn(1) = n01
  nn(2) = n02
  nn(3) = n03
  i1 = 0
  i2 = 0
  i3 = 0
  bc0 = 0
  bc1 = 0
  bc2 = 0
  if (geocode == 'F') then
    bc0 = -2
    bc2 = -2
  end if
  if (geocode /= 'P') then
    bc1 = -2
  end if
  a0 = (1.0) / (hgrids(0))
  a1 = (1.0) / (hgrids(1))
  a2 = (1.0) / (hgrids(2))
  allocate( du(0:n01 - (1), 0:n02 - (1), 0:n03 - (1), 0:2) )
  call Poisson_broker(nord, 0, nn, bc0, u,  du(0, 0, 0, 0), a0)
  call Poisson_broker(nord, 1, nn, bc1, u,  du(0, 0, 0, 1), a1)
  call Poisson_broker(nord, 2, nn, bc2, u,  du(0, 0, 0, 2), a2)
!$omp parallel  default(shared) private(i1, i2, i3)
!$omp do 
  do i3 = 0, n03 - (1), 1
    do i2 = 0, n02 - (1), 1
      do i1 = 0, n01 - (1), 1
        ddu(i1, i2, i3) = (du(i1, i2, i3, 0)) * (du(i1, i2, i3, 0)) + &
(du(i1, i2, i3, 1)) * (du(i1, i2, i3, 1)) + &
(du(i1, i2, i3, 2)) * (du(i1, i2, i3, 2))
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
  deallocate( du )
END SUBROUTINE nabla_u_square


    !> verify that the density is considerably zero in the region where epsilon is different from one
    subroutine nonvacuum_projection(n1,n23,rho,oneoeps,norm)
      implicit none
      integer, intent(in) :: n1,n23 !< parallelized box dimensions
      real(dp), dimension(n1,n23), intent(in) :: rho !<charge density
      !>inverse of epsilon (might also be inverse of sqrt(eps))
      real(dp), dimension(n1,n23), intent(in) :: oneoeps 
      real(dp), intent(out) :: norm !< \int of rho where epsilon /=1
      !local variables
      integer :: i1,i23
      real(dp), parameter :: tol= 5.d-1

      norm=0.0_dp
      !$omp parallel do default(shared) private(i1,i23)&
      !$omp reduction(+:norm)
      do i23=1,n23
         do i1=1,n1
            if (abs(oneoeps(i1,i23) - 1.0_dp) > tol) norm=norm+rho(i1,i23)
         end do
      end do
      !$omp end parallel do

    end subroutine nonvacuum_projection

end module FDder
