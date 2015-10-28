SUBROUTINE d_poisson2_p_10_u2_1_false_false_true(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
!$omp parallel  default(shared) private(i1, i2, tt0, tt1)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 = 0,  -(-1) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt0 = tt0 + (x(-1 + i1 - (((i1 + -1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (-0.5_wp)
      tt1 = tt1 + (x(-1 + i1 - (((i1 + -1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 1)) * (-0.5_wp)
      tt0 = tt0 + (x(0 + i1 - (((i1 + 0 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.0_wp)
      tt1 = tt1 + (x(0 + i1 - (((i1 + 0 + (n) * (2)) / (n) - (2)) * (n)), i2 + 1)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1 - (((i1 + 1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.5_wp)
      tt1 = tt1 + (x(1 + i1 - (((i1 + 1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 1)) * (0.5_wp)
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
    end do
    do i1 =  -(-1), n - (1) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt0 = tt0 + (x(-1 + i1, i2 + 0)) * (-0.5_wp)
      tt1 = tt1 + (x(-1 + i1, i2 + 1)) * (-0.5_wp)
      tt0 = tt0 + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt1 = tt1 + (x(0 + i1, i2 + 1)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1, i2 + 0)) * (0.5_wp)
      tt1 = tt1 + (x(1 + i1, i2 + 1)) * (0.5_wp)
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
    end do
    do i1 = n - (1), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt0 = tt0 + (x(-1 + i1 - (((i1 + -1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (-0.5_wp)
      tt1 = tt1 + (x(-1 + i1 - (((i1 + -1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 1)) * (-0.5_wp)
      tt0 = tt0 + (x(0 + i1 - (((i1 + 0 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.0_wp)
      tt1 = tt1 + (x(0 + i1 - (((i1 + 0 + (n) * (2)) / (n) - (2)) * (n)), i2 + 1)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1 - (((i1 + 1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.5_wp)
      tt1 = tt1 + (x(1 + i1 - (((i1 + 1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 1)) * (0.5_wp)
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 = 0,  -(-1) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(-1 + i1 - (((i1 + -1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (-0.5_wp)
      tt0 = tt0 + (x(0 + i1 - (((i1 + 0 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1 - (((i1 + 1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.5_wp)
      y(i1, i2 + 0) = tt0
    end do
    do i1 =  -(-1), n - (1) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(-1 + i1, i2 + 0)) * (-0.5_wp)
      tt0 = tt0 + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1, i2 + 0)) * (0.5_wp)
      y(i1, i2 + 0) = tt0
    end do
    do i1 = n - (1), n - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(-1 + i1 - (((i1 + -1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (-0.5_wp)
      tt0 = tt0 + (x(0 + i1 - (((i1 + 0 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1 - (((i1 + 1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.5_wp)
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_p_10_u2_1_false_false_true
SUBROUTINE d_poisson2_p_10_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_p_10_u1_1_true_false_true_cost
SUBROUTINE d_poisson2_p_10_a_u2_1_true_true_true(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
  integer(kind=4), dimension(-1 - (1):1 - (-1) - (1)) :: mod_arr
  do l = -1 - (1), 1 - (-1) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 = 0,  -(-1) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(0) = tt(0) + (x(mod_arr(-1 + i1), i2 + 0)) * (-0.5_wp)
      tt(1) = tt(1) + (x(mod_arr(-1 + i1), i2 + 1)) * (-0.5_wp)
      tt(0) = tt(0) + (x(mod_arr(0 + i1), i2 + 0)) * (0.0_wp)
      tt(1) = tt(1) + (x(mod_arr(0 + i1), i2 + 1)) * (0.0_wp)
      tt(0) = tt(0) + (x(mod_arr(1 + i1), i2 + 0)) * (0.5_wp)
      tt(1) = tt(1) + (x(mod_arr(1 + i1), i2 + 1)) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
    do i1 =  -(-1), n - (1) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.5_wp)
      tt(1) = tt(1) + (x(-1 + i1, i2 + 1)) * (-0.5_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(1) = tt(1) + (x(0 + i1, i2 + 1)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.5_wp)
      tt(1) = tt(1) + (x(1 + i1, i2 + 1)) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
    do i1 = n - (1), n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(0) = tt(0) + (x(mod_arr(-1 + i1 - (n)), i2 + 0)) * (-0.5_wp)
      tt(1) = tt(1) + (x(mod_arr(-1 + i1 - (n)), i2 + 1)) * (-0.5_wp)
      tt(0) = tt(0) + (x(mod_arr(0 + i1 - (n)), i2 + 0)) * (0.0_wp)
      tt(1) = tt(1) + (x(mod_arr(0 + i1 - (n)), i2 + 1)) * (0.0_wp)
      tt(0) = tt(0) + (x(mod_arr(1 + i1 - (n)), i2 + 0)) * (0.5_wp)
      tt(1) = tt(1) + (x(mod_arr(1 + i1 - (n)), i2 + 1)) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 = 0,  -(-1) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(mod_arr(-1 + i1), i2 + 0)) * (-0.5_wp)
      tt(0) = tt(0) + (x(mod_arr(0 + i1), i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(mod_arr(1 + i1), i2 + 0)) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-1), n - (1) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.5_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (1), n - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(mod_arr(-1 + i1 - (n)), i2 + 0)) * (-0.5_wp)
      tt(0) = tt(0) + (x(mod_arr(0 + i1 - (n)), i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(mod_arr(1 + i1 - (n)), i2 + 0)) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_p_10_a_u2_1_true_true_true
SUBROUTINE d_poisson2_p_10_a_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_p_10_a_u1_1_true_false_true_cost
SUBROUTINE d_poisson2_p_01_u4_0_false_true_false(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:3) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 = 0,  -(-1) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      do l = -1, 1, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson2_1_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson2_1_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson2_1_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson2_1_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
      y(i1 + 2, i2) = tt(2)
      y(i1 + 3, i2) = tt(3)
    end do
    do i2 =  -(-1), n - (1) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      do l = -1, 1, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson2_1_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson2_1_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson2_1_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson2_1_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
      y(i1 + 2, i2) = tt(2)
      y(i1 + 3, i2) = tt(3)
    end do
    do i2 = n - (1), n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      do l = -1, 1, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson2_1_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson2_1_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson2_1_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson2_1_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
      y(i1 + 2, i2) = tt(2)
      y(i1 + 3, i2) = tt(3)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 = 0,  -(-1) - (1), 1
      tt(0) = 0.0_wp
      do l = -1, 1, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson2_1_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 =  -(-1), n - (1) - (1), 1
      tt(0) = 0.0_wp
      do l = -1, 1, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson2_1_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 = n - (1), n - (1), 1
      tt(0) = 0.0_wp
      do l = -1, 1, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson2_1_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_p_01_u4_0_false_true_false
SUBROUTINE d_poisson2_p_01_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_p_01_u1_0_true_false_true_cost
SUBROUTINE d_poisson2_p_01_a_u5_0_false_true_true(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:4) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (5), 5
    do i2 = 0,  -(-1) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.5_wp)
      tt(1) = tt(1) + (x(i1 + 1, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.5_wp)
      tt(2) = tt(2) + (x(i1 + 2, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.5_wp)
      tt(3) = tt(3) + (x(i1 + 3, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.5_wp)
      tt(4) = tt(4) + (x(i1 + 4, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.5_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(1) = tt(1) + (x(i1 + 1, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(2) = tt(2) + (x(i1 + 2, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(3) = tt(3) + (x(i1 + 3, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(4) = tt(4) + (x(i1 + 4, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.5_wp)
      tt(1) = tt(1) + (x(i1 + 1, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.5_wp)
      tt(2) = tt(2) + (x(i1 + 2, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.5_wp)
      tt(3) = tt(3) + (x(i1 + 3, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.5_wp)
      tt(4) = tt(4) + (x(i1 + 4, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
      tt(4) = (tt(4)) * (a)
      y(i1 + 4, i2) = tt(4)
    end do
    do i2 =  -(-1), n - (1) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2)) * (-0.5_wp)
      tt(1) = tt(1) + (x(i1 + 1, -1 + i2)) * (-0.5_wp)
      tt(2) = tt(2) + (x(i1 + 2, -1 + i2)) * (-0.5_wp)
      tt(3) = tt(3) + (x(i1 + 3, -1 + i2)) * (-0.5_wp)
      tt(4) = tt(4) + (x(i1 + 4, -1 + i2)) * (-0.5_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt(1) = tt(1) + (x(i1 + 1, 0 + i2)) * (0.0_wp)
      tt(2) = tt(2) + (x(i1 + 2, 0 + i2)) * (0.0_wp)
      tt(3) = tt(3) + (x(i1 + 3, 0 + i2)) * (0.0_wp)
      tt(4) = tt(4) + (x(i1 + 4, 0 + i2)) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2)) * (0.5_wp)
      tt(1) = tt(1) + (x(i1 + 1, 1 + i2)) * (0.5_wp)
      tt(2) = tt(2) + (x(i1 + 2, 1 + i2)) * (0.5_wp)
      tt(3) = tt(3) + (x(i1 + 3, 1 + i2)) * (0.5_wp)
      tt(4) = tt(4) + (x(i1 + 4, 1 + i2)) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
      tt(4) = (tt(4)) * (a)
      y(i1 + 4, i2) = tt(4)
    end do
    do i2 = n - (1), n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.5_wp)
      tt(1) = tt(1) + (x(i1 + 1, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.5_wp)
      tt(2) = tt(2) + (x(i1 + 2, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.5_wp)
      tt(3) = tt(3) + (x(i1 + 3, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.5_wp)
      tt(4) = tt(4) + (x(i1 + 4, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.5_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(1) = tt(1) + (x(i1 + 1, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(2) = tt(2) + (x(i1 + 2, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(3) = tt(3) + (x(i1 + 3, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(4) = tt(4) + (x(i1 + 4, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.5_wp)
      tt(1) = tt(1) + (x(i1 + 1, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.5_wp)
      tt(2) = tt(2) + (x(i1 + 2, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.5_wp)
      tt(3) = tt(3) + (x(i1 + 3, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.5_wp)
      tt(4) = tt(4) + (x(i1 + 4, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
      tt(4) = (tt(4)) * (a)
      y(i1 + 4, i2) = tt(4)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (5)) * (5), ndat - (1), 1
    do i2 = 0,  -(-1) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.5_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 =  -(-1), n - (1) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2)) * (-0.5_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2)) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 = n - (1), n - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.5_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_p_01_a_u5_0_false_true_true
SUBROUTINE d_poisson2_p_01_a_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_p_01_a_u1_0_true_false_true_cost
SUBROUTINE d_poisson2_p_201_u4_0_true_true_false(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8), dimension(0:3) :: tt
  integer(kind=4), dimension(-1 - (1):1 - (-1) - (1)) :: mod_arr
  do l = -1 - (1), 1 - (-1) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 = 0,  -(-1) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        do l = -1, 1, 1
          tt(0) = tt(0) + (x(i1 + 0, mod_arr(l + i2), i3)) * (poisson2_1_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, mod_arr(l + i2), i3)) * (poisson2_1_fil(l))
          tt(2) = tt(2) + (x(i1 + 2, mod_arr(l + i2), i3)) * (poisson2_1_fil(l))
          tt(3) = tt(3) + (x(i1 + 3, mod_arr(l + i2), i3)) * (poisson2_1_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
      end do
      do i2 =  -(-1), n - (1) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        do l = -1, 1, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson2_1_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, l + i2, i3)) * (poisson2_1_fil(l))
          tt(2) = tt(2) + (x(i1 + 2, l + i2, i3)) * (poisson2_1_fil(l))
          tt(3) = tt(3) + (x(i1 + 3, l + i2, i3)) * (poisson2_1_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
      end do
      do i2 = n - (1), n - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        do l = -1, 1, 1
          tt(0) = tt(0) + (x(i1 + 0, mod_arr(l + i2 - (n)), i3)) * (poisson2_1_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, mod_arr(l + i2 - (n)), i3)) * (poisson2_1_fil(l))
          tt(2) = tt(2) + (x(i1 + 2, mod_arr(l + i2 - (n)), i3)) * (poisson2_1_fil(l))
          tt(3) = tt(3) + (x(i1 + 3, mod_arr(l + i2 - (n)), i3)) * (poisson2_1_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 = 0,  -(-1) - (1), 1
        tt(0) = 0.0_wp
        do l = -1, 1, 1
          tt(0) = tt(0) + (x(i1 + 0, mod_arr(l + i2), i3)) * (poisson2_1_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 =  -(-1), n - (1) - (1), 1
        tt(0) = 0.0_wp
        do l = -1, 1, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson2_1_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 = n - (1), n - (1), 1
        tt(0) = 0.0_wp
        do l = -1, 1, 1
          tt(0) = tt(0) + (x(i1 + 0, mod_arr(l + i2 - (n)), i3)) * (poisson2_1_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_p_201_u4_0_true_true_false
SUBROUTINE d_poisson2_p_201_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_p_201_u1_2_true_false_true_cost
SUBROUTINE d_poisson2_p_201_a_u4_0_true_true_true(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8), dimension(0:3) :: tt
  integer(kind=4), dimension(-1 - (1):1 - (-1) - (1)) :: mod_arr
  do l = -1 - (1), 1 - (-1) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 = 0,  -(-1) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(-1 + i2), i3)) * (-0.5_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(-1 + i2), i3)) * (-0.5_wp)
        tt(2) = tt(2) + (x(i1 + 2, mod_arr(-1 + i2), i3)) * (-0.5_wp)
        tt(3) = tt(3) + (x(i1 + 3, mod_arr(-1 + i2), i3)) * (-0.5_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt(2) = tt(2) + (x(i1 + 2, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt(3) = tt(3) + (x(i1 + 3, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(1 + i2), i3)) * (0.5_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(1 + i2), i3)) * (0.5_wp)
        tt(2) = tt(2) + (x(i1 + 2, mod_arr(1 + i2), i3)) * (0.5_wp)
        tt(3) = tt(3) + (x(i1 + 3, mod_arr(1 + i2), i3)) * (0.5_wp)
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
        tt(1) = (tt(1)) * (a)
        y(i1 + 1, i2, i3) = tt(1)
        tt(2) = (tt(2)) * (a)
        y(i1 + 2, i2, i3) = tt(2)
        tt(3) = (tt(3)) * (a)
        y(i1 + 3, i2, i3) = tt(3)
      end do
      do i2 =  -(-1), n - (1) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.5_wp)
        tt(1) = tt(1) + (x(i1 + 1, -1 + i2, i3)) * (-0.5_wp)
        tt(2) = tt(2) + (x(i1 + 2, -1 + i2, i3)) * (-0.5_wp)
        tt(3) = tt(3) + (x(i1 + 3, -1 + i2, i3)) * (-0.5_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(1) = tt(1) + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt(2) = tt(2) + (x(i1 + 2, 0 + i2, i3)) * (0.0_wp)
        tt(3) = tt(3) + (x(i1 + 3, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.5_wp)
        tt(1) = tt(1) + (x(i1 + 1, 1 + i2, i3)) * (0.5_wp)
        tt(2) = tt(2) + (x(i1 + 2, 1 + i2, i3)) * (0.5_wp)
        tt(3) = tt(3) + (x(i1 + 3, 1 + i2, i3)) * (0.5_wp)
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
        tt(1) = (tt(1)) * (a)
        y(i1 + 1, i2, i3) = tt(1)
        tt(2) = (tt(2)) * (a)
        y(i1 + 2, i2, i3) = tt(2)
        tt(3) = (tt(3)) * (a)
        y(i1 + 3, i2, i3) = tt(3)
      end do
      do i2 = n - (1), n - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(-1 + i2 - (n)), i3)) * (-0.5_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(-1 + i2 - (n)), i3)) * (-0.5_wp)
        tt(2) = tt(2) + (x(i1 + 2, mod_arr(-1 + i2 - (n)), i3)) * (-0.5_wp)
        tt(3) = tt(3) + (x(i1 + 3, mod_arr(-1 + i2 - (n)), i3)) * (-0.5_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt(2) = tt(2) + (x(i1 + 2, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt(3) = tt(3) + (x(i1 + 3, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(1 + i2 - (n)), i3)) * (0.5_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(1 + i2 - (n)), i3)) * (0.5_wp)
        tt(2) = tt(2) + (x(i1 + 2, mod_arr(1 + i2 - (n)), i3)) * (0.5_wp)
        tt(3) = tt(3) + (x(i1 + 3, mod_arr(1 + i2 - (n)), i3)) * (0.5_wp)
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
        tt(1) = (tt(1)) * (a)
        y(i1 + 1, i2, i3) = tt(1)
        tt(2) = (tt(2)) * (a)
        y(i1 + 2, i2, i3) = tt(2)
        tt(3) = (tt(3)) * (a)
        y(i1 + 3, i2, i3) = tt(3)
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 = 0,  -(-1) - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(-1 + i2), i3)) * (-0.5_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(1 + i2), i3)) * (0.5_wp)
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 =  -(-1), n - (1) - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.5_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.5_wp)
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 = n - (1), n - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(-1 + i2 - (n)), i3)) * (-0.5_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(1 + i2 - (n)), i3)) * (0.5_wp)
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_p_201_a_u4_0_true_true_true
SUBROUTINE d_poisson2_p_201_a_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_p_201_a_u1_2_true_false_true_cost
SUBROUTINE d_poisson2_fg_10_u5_1_false_false_false(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension( -(1):n - (-1) - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2, tt3, tt4)
!$omp do 
  do i2 = 0, ndat - (5), 5
    do i1 =  -(1),  -(-1) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt4 = 0.0_wp
      do l = max( -(i1), -1), 1, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson2_1_fil(l))
        tt1 = tt1 + (x(l + i1, i2 + 1)) * (poisson2_1_fil(l))
        tt2 = tt2 + (x(l + i1, i2 + 2)) * (poisson2_1_fil(l))
        tt3 = tt3 + (x(l + i1, i2 + 3)) * (poisson2_1_fil(l))
        tt4 = tt4 + (x(l + i1, i2 + 4)) * (poisson2_1_fil(l))
      end do
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
      y(i1, i2 + 2) = tt2
      y(i1, i2 + 3) = tt3
      y(i1, i2 + 4) = tt4
    end do
    do i1 =  -(-1), n - (1) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt4 = 0.0_wp
      do l = -1, 1, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson2_1_fil(l))
        tt1 = tt1 + (x(l + i1, i2 + 1)) * (poisson2_1_fil(l))
        tt2 = tt2 + (x(l + i1, i2 + 2)) * (poisson2_1_fil(l))
        tt3 = tt3 + (x(l + i1, i2 + 3)) * (poisson2_1_fil(l))
        tt4 = tt4 + (x(l + i1, i2 + 4)) * (poisson2_1_fil(l))
      end do
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
      y(i1, i2 + 2) = tt2
      y(i1, i2 + 3) = tt3
      y(i1, i2 + 4) = tt4
    end do
    do i1 = n - (1), n - (-1) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt4 = 0.0_wp
      do l = -1, min(1, n - (1) - (i1)), 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson2_1_fil(l))
        tt1 = tt1 + (x(l + i1, i2 + 1)) * (poisson2_1_fil(l))
        tt2 = tt2 + (x(l + i1, i2 + 2)) * (poisson2_1_fil(l))
        tt3 = tt3 + (x(l + i1, i2 + 3)) * (poisson2_1_fil(l))
        tt4 = tt4 + (x(l + i1, i2 + 4)) * (poisson2_1_fil(l))
      end do
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
      y(i1, i2 + 2) = tt2
      y(i1, i2 + 3) = tt3
      y(i1, i2 + 4) = tt4
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (5)) * (5), ndat - (1), 1
    do i1 =  -(1),  -(-1) - (1), 1
      tt0 = 0.0_wp
      do l = max( -(i1), -1), 1, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson2_1_fil(l))
      end do
      y(i1, i2 + 0) = tt0
    end do
    do i1 =  -(-1), n - (1) - (1), 1
      tt0 = 0.0_wp
      do l = -1, 1, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson2_1_fil(l))
      end do
      y(i1, i2 + 0) = tt0
    end do
    do i1 = n - (1), n - (-1) - (1), 1
      tt0 = 0.0_wp
      do l = -1, min(1, n - (1) - (i1)), 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson2_1_fil(l))
      end do
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_fg_10_u5_1_false_false_false
SUBROUTINE d_poisson2_fg_10_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_fg_10_u1_1_false_false_true_cost
SUBROUTINE d_poisson2_fg_10_a_u2_1_false_true_true(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension( -(1):n - (-1) - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 =  -(1),  -(-1) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = max( -(i1), -1), 1, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson2_1_fil(l))
        tt(1) = tt(1) + (x(l + i1, i2 + 1)) * (poisson2_1_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
    do i1 =  -(-1), n - (1) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.5_wp)
      tt(1) = tt(1) + (x(-1 + i1, i2 + 1)) * (-0.5_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(1) = tt(1) + (x(0 + i1, i2 + 1)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.5_wp)
      tt(1) = tt(1) + (x(1 + i1, i2 + 1)) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
    do i1 = n - (1), n - (-1) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -1, min(1, n - (1) - (i1)), 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson2_1_fil(l))
        tt(1) = tt(1) + (x(l + i1, i2 + 1)) * (poisson2_1_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 =  -(1),  -(-1) - (1), 1
      tt(0) = 0.0_wp
      do l = max( -(i1), -1), 1, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson2_1_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-1), n - (1) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.5_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (1), n - (-1) - (1), 1
      tt(0) = 0.0_wp
      do l = -1, min(1, n - (1) - (i1)), 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson2_1_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_fg_10_a_u2_1_false_true_true
SUBROUTINE d_poisson2_fg_10_a_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_fg_10_a_u1_1_false_false_true_cost
SUBROUTINE d_poisson2_fg_01_u2_0_false_true_false(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1),  -(1):n - (-1) - (1)) :: y
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (2), 2
    do i2 =  -(1),  -(-1) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = max( -(i2), -1), 1, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson2_1_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson2_1_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
    end do
    do i2 =  -(-1), n - (1) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -1, 1, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson2_1_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson2_1_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
    end do
    do i2 = n - (1), n - (-1) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -1, min(1, n - (1) - (i2)), 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson2_1_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson2_1_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i2 =  -(1),  -(-1) - (1), 1
      tt(0) = 0.0_wp
      do l = max( -(i2), -1), 1, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson2_1_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 =  -(-1), n - (1) - (1), 1
      tt(0) = 0.0_wp
      do l = -1, 1, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson2_1_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 = n - (1), n - (-1) - (1), 1
      tt(0) = 0.0_wp
      do l = -1, min(1, n - (1) - (i2)), 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson2_1_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_fg_01_u2_0_false_true_false
SUBROUTINE d_poisson2_fg_01_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_fg_01_u1_0_false_false_true_cost
SUBROUTINE d_poisson2_fg_01_a_u4_0_false_false_true(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1),  -(1):n - (-1) - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2, tt3)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 =  -(1),  -(-1) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = max( -(i2), -1), 1, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson2_1_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson2_1_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2)) * (poisson2_1_fil(l))
        tt3 = tt3 + (x(i1 + 3, l + i2)) * (poisson2_1_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
    do i2 =  -(-1), n - (1) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.5_wp)
      tt1 = tt1 + (x(i1 + 1, -1 + i2)) * (-0.5_wp)
      tt2 = tt2 + (x(i1 + 2, -1 + i2)) * (-0.5_wp)
      tt3 = tt3 + (x(i1 + 3, -1 + i2)) * (-0.5_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt1 = tt1 + (x(i1 + 1, 0 + i2)) * (0.0_wp)
      tt2 = tt2 + (x(i1 + 2, 0 + i2)) * (0.0_wp)
      tt3 = tt3 + (x(i1 + 3, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.5_wp)
      tt1 = tt1 + (x(i1 + 1, 1 + i2)) * (0.5_wp)
      tt2 = tt2 + (x(i1 + 2, 1 + i2)) * (0.5_wp)
      tt3 = tt3 + (x(i1 + 3, 1 + i2)) * (0.5_wp)
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
    do i2 = n - (1), n - (-1) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = -1, min(1, n - (1) - (i2)), 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson2_1_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson2_1_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2)) * (poisson2_1_fil(l))
        tt3 = tt3 + (x(i1 + 3, l + i2)) * (poisson2_1_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 =  -(1),  -(-1) - (1), 1
      tt0 = 0.0_wp
      do l = max( -(i2), -1), 1, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson2_1_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
    do i2 =  -(-1), n - (1) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.5_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.5_wp)
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
    do i2 = n - (1), n - (-1) - (1), 1
      tt0 = 0.0_wp
      do l = -1, min(1, n - (1) - (i2)), 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson2_1_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_fg_01_a_u4_0_false_false_true
SUBROUTINE d_poisson2_fg_01_a_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_fg_01_a_u1_0_false_false_true_cost
SUBROUTINE d_poisson2_fg_201_u4_0_false_true_true(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1),  -(1):n - (-1) - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8), dimension(0:3) :: tt
!$omp parallel  default(shared) private(i1, i2, i3, tt)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 =  -(1),  -(-1) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        do l = max( -(i2), -1), 1, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson2_1_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, l + i2, i3)) * (poisson2_1_fil(l))
          tt(2) = tt(2) + (x(i1 + 2, l + i2, i3)) * (poisson2_1_fil(l))
          tt(3) = tt(3) + (x(i1 + 3, l + i2, i3)) * (poisson2_1_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
      end do
      do i2 =  -(-1), n - (1) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.5_wp)
        tt(1) = tt(1) + (x(i1 + 1, -1 + i2, i3)) * (-0.5_wp)
        tt(2) = tt(2) + (x(i1 + 2, -1 + i2, i3)) * (-0.5_wp)
        tt(3) = tt(3) + (x(i1 + 3, -1 + i2, i3)) * (-0.5_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(1) = tt(1) + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt(2) = tt(2) + (x(i1 + 2, 0 + i2, i3)) * (0.0_wp)
        tt(3) = tt(3) + (x(i1 + 3, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.5_wp)
        tt(1) = tt(1) + (x(i1 + 1, 1 + i2, i3)) * (0.5_wp)
        tt(2) = tt(2) + (x(i1 + 2, 1 + i2, i3)) * (0.5_wp)
        tt(3) = tt(3) + (x(i1 + 3, 1 + i2, i3)) * (0.5_wp)
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
      end do
      do i2 = n - (1), n - (-1) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        do l = -1, min(1, n - (1) - (i2)), 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson2_1_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, l + i2, i3)) * (poisson2_1_fil(l))
          tt(2) = tt(2) + (x(i1 + 2, l + i2, i3)) * (poisson2_1_fil(l))
          tt(3) = tt(3) + (x(i1 + 3, l + i2, i3)) * (poisson2_1_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 =  -(1),  -(-1) - (1), 1
        tt(0) = 0.0_wp
        do l = max( -(i2), -1), 1, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson2_1_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 =  -(-1), n - (1) - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.5_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.5_wp)
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 = n - (1), n - (-1) - (1), 1
        tt(0) = 0.0_wp
        do l = -1, min(1, n - (1) - (i2)), 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson2_1_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_fg_201_u4_0_false_true_true
SUBROUTINE d_poisson2_fg_201_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_fg_201_u1_2_false_false_true_cost
SUBROUTINE d_poisson2_fg_201_a_u4_0_false_false_true(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1),  -(1):n - (-1) - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1, tt2, tt3)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 =  -(1),  -(-1) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        do l = max( -(i2), -1), 1, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson2_1_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson2_1_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2, i3)) * (poisson2_1_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2, i3)) * (poisson2_1_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
      do i2 =  -(-1), n - (1) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.5_wp)
        tt1 = tt1 + (x(i1 + 1, -1 + i2, i3)) * (-0.5_wp)
        tt2 = tt2 + (x(i1 + 2, -1 + i2, i3)) * (-0.5_wp)
        tt3 = tt3 + (x(i1 + 3, -1 + i2, i3)) * (-0.5_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt1 = tt1 + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt2 = tt2 + (x(i1 + 2, 0 + i2, i3)) * (0.0_wp)
        tt3 = tt3 + (x(i1 + 3, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.5_wp)
        tt1 = tt1 + (x(i1 + 1, 1 + i2, i3)) * (0.5_wp)
        tt2 = tt2 + (x(i1 + 2, 1 + i2, i3)) * (0.5_wp)
        tt3 = tt3 + (x(i1 + 3, 1 + i2, i3)) * (0.5_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
      do i2 = n - (1), n - (-1) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        do l = -1, min(1, n - (1) - (i2)), 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson2_1_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson2_1_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2, i3)) * (poisson2_1_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2, i3)) * (poisson2_1_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 =  -(1),  -(-1) - (1), 1
        tt0 = 0.0_wp
        do l = max( -(i2), -1), 1, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson2_1_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-1), n - (1) - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.5_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.5_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (1), n - (-1) - (1), 1
        tt0 = 0.0_wp
        do l = -1, min(1, n - (1) - (i2)), 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson2_1_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_fg_201_a_u4_0_false_false_true
SUBROUTINE d_poisson2_fg_201_a_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_fg_201_a_u1_2_false_false_true_cost
SUBROUTINE d_poisson2_fs_10_u1_1_false_true_true(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(-1:n + 1 - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:0) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 = 0, n - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.5_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.5_wp)
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_fs_10_u1_1_false_true_true
SUBROUTINE d_poisson2_fs_10_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_fs_10_u1_1_false_false_true_cost
SUBROUTINE d_poisson2_fs_10_a_u2_1_false_true_true(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(-1:n + 1 - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 = 0, n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.5_wp)
      tt(1) = tt(1) + (x(-1 + i1, i2 + 1)) * (-0.5_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(1) = tt(1) + (x(0 + i1, i2 + 1)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.5_wp)
      tt(1) = tt(1) + (x(1 + i1, i2 + 1)) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 = 0, n - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.5_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_fs_10_a_u2_1_false_true_true
SUBROUTINE d_poisson2_fs_10_a_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_fs_10_a_u1_1_false_false_true_cost
SUBROUTINE d_poisson2_fs_01_u5_0_false_false_false(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), -1:n + 1 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2, tt3, tt4)
!$omp do 
  do i1 = 0, ndat - (5), 5
    do i2 = 0, n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt4 = 0.0_wp
      do l = -1, 1, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson2_1_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson2_1_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2)) * (poisson2_1_fil(l))
        tt3 = tt3 + (x(i1 + 3, l + i2)) * (poisson2_1_fil(l))
        tt4 = tt4 + (x(i1 + 4, l + i2)) * (poisson2_1_fil(l))
      end do
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
      y(i1 + 3, i2) = tt3
      y(i1 + 4, i2) = tt4
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (5)) * (5), ndat - (1), 1
    do i2 = 0, n - (1), 1
      tt0 = 0.0_wp
      do l = -1, 1, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson2_1_fil(l))
      end do
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_fs_01_u5_0_false_false_false
SUBROUTINE d_poisson2_fs_01_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_fs_01_u1_0_false_false_true_cost
SUBROUTINE d_poisson2_fs_01_a_u4_0_false_false_true(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), -1:n + 1 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2, tt3)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 = 0, n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.5_wp)
      tt1 = tt1 + (x(i1 + 1, -1 + i2)) * (-0.5_wp)
      tt2 = tt2 + (x(i1 + 2, -1 + i2)) * (-0.5_wp)
      tt3 = tt3 + (x(i1 + 3, -1 + i2)) * (-0.5_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt1 = tt1 + (x(i1 + 1, 0 + i2)) * (0.0_wp)
      tt2 = tt2 + (x(i1 + 2, 0 + i2)) * (0.0_wp)
      tt3 = tt3 + (x(i1 + 3, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.5_wp)
      tt1 = tt1 + (x(i1 + 1, 1 + i2)) * (0.5_wp)
      tt2 = tt2 + (x(i1 + 2, 1 + i2)) * (0.5_wp)
      tt3 = tt3 + (x(i1 + 3, 1 + i2)) * (0.5_wp)
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 = 0, n - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.5_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.5_wp)
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_fs_01_a_u4_0_false_false_true
SUBROUTINE d_poisson2_fs_01_a_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_fs_01_a_u1_0_false_false_true_cost
SUBROUTINE d_poisson2_fs_201_u4_0_false_false_false(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), -1:n + 1 - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1, tt2, tt3)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 = 0, n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        do l = -1, 1, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson2_1_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson2_1_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2, i3)) * (poisson2_1_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2, i3)) * (poisson2_1_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
        y(i1 + 2, i2, i3) = tt2
        y(i1 + 3, i2, i3) = tt3
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 = 0, n - (1), 1
        tt0 = 0.0_wp
        do l = -1, 1, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson2_1_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_fs_201_u4_0_false_false_false
SUBROUTINE d_poisson2_fs_201_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_fs_201_u1_2_false_false_true_cost
SUBROUTINE d_poisson2_fs_201_a_u4_0_false_false_true(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), -1:n + 1 - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1, tt2, tt3)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 = 0, n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.5_wp)
        tt1 = tt1 + (x(i1 + 1, -1 + i2, i3)) * (-0.5_wp)
        tt2 = tt2 + (x(i1 + 2, -1 + i2, i3)) * (-0.5_wp)
        tt3 = tt3 + (x(i1 + 3, -1 + i2, i3)) * (-0.5_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt1 = tt1 + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt2 = tt2 + (x(i1 + 2, 0 + i2, i3)) * (0.0_wp)
        tt3 = tt3 + (x(i1 + 3, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.5_wp)
        tt1 = tt1 + (x(i1 + 1, 1 + i2, i3)) * (0.5_wp)
        tt2 = tt2 + (x(i1 + 2, 1 + i2, i3)) * (0.5_wp)
        tt3 = tt3 + (x(i1 + 3, 1 + i2, i3)) * (0.5_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 = 0, n - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.5_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.5_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_fs_201_a_u4_0_false_false_true
SUBROUTINE d_poisson2_fs_201_a_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_fs_201_a_u1_2_false_false_true_cost
SUBROUTINE d_poisson2_np_10_u2_1_true_false_false(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(0:8) :: poisson2_fil = (/ &
-1.5_wp, &
2.0_wp, &
-0.5_wp, &
-0.5_wp, &
0.0_wp, &
0.5_wp, &
0.5_wp, &
-2.0_wp, &
1.5_wp /)
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  integer(kind=4), dimension(-1 - (1):1 - (-1) - (1)) :: mod_arr
  do l = -1 - (1), 1 - (-1) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt0, tt1)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 = 0,  -(-1) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -1, 1, 1
        tt0 = tt0 + (x(l + 1, i2 + 0)) * (poisson2_fil((i1 - (1)) * (3) + 3 + l + 1))
        tt1 = tt1 + (x(l + 1, i2 + 1)) * (poisson2_fil((i1 - (1)) * (3) + 3 + l + 1))
      end do
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
    end do
    do i1 =  -(-1), n - (1) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -1, 1, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson2_1_fil(l))
        tt1 = tt1 + (x(l + i1, i2 + 1)) * (poisson2_1_fil(l))
      end do
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
    end do
    do i1 = n - (1), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -1, 1, 1
        tt0 = tt0 + (x(l + -1 + n - (1), i2 + 0)) * (poisson2_fil((i1 + 1 - (n) + 1) * (3) + 3 + l + 1))
        tt1 = tt1 + (x(l + -1 + n - (1), i2 + 1)) * (poisson2_fil((i1 + 1 - (n) + 1) * (3) + 3 + l + 1))
      end do
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 = 0,  -(-1) - (1), 1
      tt0 = 0.0_wp
      do l = -1, 1, 1
        tt0 = tt0 + (x(l + 1, i2 + 0)) * (poisson2_fil((i1 - (1)) * (3) + 3 + l + 1))
      end do
      y(i1, i2 + 0) = tt0
    end do
    do i1 =  -(-1), n - (1) - (1), 1
      tt0 = 0.0_wp
      do l = -1, 1, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson2_1_fil(l))
      end do
      y(i1, i2 + 0) = tt0
    end do
    do i1 = n - (1), n - (1), 1
      tt0 = 0.0_wp
      do l = -1, 1, 1
        tt0 = tt0 + (x(l + -1 + n - (1), i2 + 0)) * (poisson2_fil((i1 + 1 - (n) + 1) * (3) + 3 + l + 1))
      end do
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_np_10_u2_1_true_false_false
SUBROUTINE d_poisson2_np_10_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_np_10_u1_1_true_false_true_cost
SUBROUTINE d_poisson2_np_10_a_u1_1_true_true_true(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(0:8) :: poisson2_fil = (/ &
-1.5_wp, &
2.0_wp, &
-0.5_wp, &
-0.5_wp, &
0.0_wp, &
0.5_wp, &
0.5_wp, &
-2.0_wp, &
1.5_wp /)
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:0) :: tt
  integer(kind=4), dimension(-1 - (1):1 - (-1) - (1)) :: mod_arr
  do l = -1 - (1), 1 - (-1) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 = 0,  -(-1) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-1 + 1, i2 + 0)) * (poisson2_fil((i1 - (1)) * (3) + 3 + -1 + 1))
      tt(0) = tt(0) + (x(0 + 1, i2 + 0)) * (poisson2_fil((i1 - (1)) * (3) + 3 + 0 + 1))
      tt(0) = tt(0) + (x(1 + 1, i2 + 0)) * (poisson2_fil((i1 - (1)) * (3) + 3 + 1 + 1))
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-1), n - (1) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.5_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.5_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (1), n - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-1 + -1 + n - (1), i2 + 0)) * (poisson2_fil((i1 + 1 - (n) + 1) * (3) + 3 + -1 + 1))
      tt(0) = tt(0) + (x(0 + -1 + n - (1), i2 + 0)) * (poisson2_fil((i1 + 1 - (n) + 1) * (3) + 3 + 0 + 1))
      tt(0) = tt(0) + (x(1 + -1 + n - (1), i2 + 0)) * (poisson2_fil((i1 + 1 - (n) + 1) * (3) + 3 + 1 + 1))
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_np_10_a_u1_1_true_true_true
SUBROUTINE d_poisson2_np_10_a_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_np_10_a_u1_1_true_false_true_cost
SUBROUTINE d_poisson2_np_01_u4_0_false_false_true(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), parameter, dimension(0:8) :: poisson2_fil = (/ &
-1.5_wp, &
2.0_wp, &
-0.5_wp, &
-0.5_wp, &
0.0_wp, &
0.5_wp, &
0.5_wp, &
-2.0_wp, &
1.5_wp /)
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2, tt3)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 = 0,  -(-1) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -1 + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + -1 + 1))
      tt1 = tt1 + (x(i1 + 1, -1 + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + -1 + 1))
      tt2 = tt2 + (x(i1 + 2, -1 + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + -1 + 1))
      tt3 = tt3 + (x(i1 + 3, -1 + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + -1 + 1))
      tt0 = tt0 + (x(i1 + 0, 0 + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 0 + 1))
      tt1 = tt1 + (x(i1 + 1, 0 + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 0 + 1))
      tt2 = tt2 + (x(i1 + 2, 0 + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 0 + 1))
      tt3 = tt3 + (x(i1 + 3, 0 + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 0 + 1))
      tt0 = tt0 + (x(i1 + 0, 1 + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 1 + 1))
      tt1 = tt1 + (x(i1 + 1, 1 + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 1 + 1))
      tt2 = tt2 + (x(i1 + 2, 1 + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 1 + 1))
      tt3 = tt3 + (x(i1 + 3, 1 + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 1 + 1))
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
      y(i1 + 3, i2) = tt3
    end do
    do i2 =  -(-1), n - (1) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.5_wp)
      tt1 = tt1 + (x(i1 + 1, -1 + i2)) * (-0.5_wp)
      tt2 = tt2 + (x(i1 + 2, -1 + i2)) * (-0.5_wp)
      tt3 = tt3 + (x(i1 + 3, -1 + i2)) * (-0.5_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt1 = tt1 + (x(i1 + 1, 0 + i2)) * (0.0_wp)
      tt2 = tt2 + (x(i1 + 2, 0 + i2)) * (0.0_wp)
      tt3 = tt3 + (x(i1 + 3, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.5_wp)
      tt1 = tt1 + (x(i1 + 1, 1 + i2)) * (0.5_wp)
      tt2 = tt2 + (x(i1 + 2, 1 + i2)) * (0.5_wp)
      tt3 = tt3 + (x(i1 + 3, 1 + i2)) * (0.5_wp)
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
      y(i1 + 3, i2) = tt3
    end do
    do i2 = n - (1), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -1 + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + -1 + 1))
      tt1 = tt1 + (x(i1 + 1, -1 + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + -1 + 1))
      tt2 = tt2 + (x(i1 + 2, -1 + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + -1 + 1))
      tt3 = tt3 + (x(i1 + 3, -1 + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + -1 + 1))
      tt0 = tt0 + (x(i1 + 0, 0 + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 0 + 1))
      tt1 = tt1 + (x(i1 + 1, 0 + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 0 + 1))
      tt2 = tt2 + (x(i1 + 2, 0 + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 0 + 1))
      tt3 = tt3 + (x(i1 + 3, 0 + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 0 + 1))
      tt0 = tt0 + (x(i1 + 0, 1 + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 1 + 1))
      tt1 = tt1 + (x(i1 + 1, 1 + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 1 + 1))
      tt2 = tt2 + (x(i1 + 2, 1 + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 1 + 1))
      tt3 = tt3 + (x(i1 + 3, 1 + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 1 + 1))
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
      y(i1 + 3, i2) = tt3
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 = 0,  -(-1) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -1 + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + -1 + 1))
      tt0 = tt0 + (x(i1 + 0, 0 + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 0 + 1))
      tt0 = tt0 + (x(i1 + 0, 1 + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 1 + 1))
      y(i1 + 0, i2) = tt0
    end do
    do i2 =  -(-1), n - (1) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.5_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.5_wp)
      y(i1 + 0, i2) = tt0
    end do
    do i2 = n - (1), n - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -1 + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + -1 + 1))
      tt0 = tt0 + (x(i1 + 0, 0 + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 0 + 1))
      tt0 = tt0 + (x(i1 + 0, 1 + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 1 + 1))
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_np_01_u4_0_false_false_true
SUBROUTINE d_poisson2_np_01_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_np_01_u1_0_true_false_true_cost
SUBROUTINE d_poisson2_np_01_a_u4_0_true_false_false(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(0:8) :: poisson2_fil = (/ &
-1.5_wp, &
2.0_wp, &
-0.5_wp, &
-0.5_wp, &
0.0_wp, &
0.5_wp, &
0.5_wp, &
-2.0_wp, &
1.5_wp /)
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  integer(kind=4), dimension(-1 - (1):1 - (-1) - (1)) :: mod_arr
  do l = -1 - (1), 1 - (-1) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2, tt3)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 = 0,  -(-1) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = -1, 1, 1
        tt0 = tt0 + (x(i1 + 0, l + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + l + 1))
        tt1 = tt1 + (x(i1 + 1, l + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + l + 1))
        tt2 = tt2 + (x(i1 + 2, l + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + l + 1))
        tt3 = tt3 + (x(i1 + 3, l + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + l + 1))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
    do i2 =  -(-1), n - (1) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = -1, 1, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson2_1_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson2_1_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2)) * (poisson2_1_fil(l))
        tt3 = tt3 + (x(i1 + 3, l + i2)) * (poisson2_1_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
    do i2 = n - (1), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = -1, 1, 1
        tt0 = tt0 + (x(i1 + 0, l + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + l + 1))
        tt1 = tt1 + (x(i1 + 1, l + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + l + 1))
        tt2 = tt2 + (x(i1 + 2, l + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + l + 1))
        tt3 = tt3 + (x(i1 + 3, l + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + l + 1))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 = 0,  -(-1) - (1), 1
      tt0 = 0.0_wp
      do l = -1, 1, 1
        tt0 = tt0 + (x(i1 + 0, l + 1)) * (poisson2_fil((i2 - (1)) * (3) + 3 + l + 1))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
    do i2 =  -(-1), n - (1) - (1), 1
      tt0 = 0.0_wp
      do l = -1, 1, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson2_1_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
    do i2 = n - (1), n - (1), 1
      tt0 = 0.0_wp
      do l = -1, 1, 1
        tt0 = tt0 + (x(i1 + 0, l + -1 + n - (1))) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + l + 1))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_np_01_a_u4_0_true_false_false
SUBROUTINE d_poisson2_np_01_a_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_np_01_a_u1_0_true_false_true_cost
SUBROUTINE d_poisson2_np_201_u2_0_true_true_true(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(0:8) :: poisson2_fil = (/ &
-1.5_wp, &
2.0_wp, &
-0.5_wp, &
-0.5_wp, &
0.0_wp, &
0.5_wp, &
0.5_wp, &
-2.0_wp, &
1.5_wp /)
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
  integer(kind=4), dimension(-1 - (1):1 - (-1) - (1)) :: mod_arr
  do l = -1 - (1), 1 - (-1) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (2), 2
      do i2 = 0,  -(-1) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -1 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + -1 + 1))
        tt(1) = tt(1) + (x(i1 + 1, -1 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + -1 + 1))
        tt(0) = tt(0) + (x(i1 + 0, 0 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 0 + 1))
        tt(1) = tt(1) + (x(i1 + 1, 0 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 0 + 1))
        tt(0) = tt(0) + (x(i1 + 0, 1 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 1 + 1))
        tt(1) = tt(1) + (x(i1 + 1, 1 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 1 + 1))
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
      end do
      do i2 =  -(-1), n - (1) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.5_wp)
        tt(1) = tt(1) + (x(i1 + 1, -1 + i2, i3)) * (-0.5_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(1) = tt(1) + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.5_wp)
        tt(1) = tt(1) + (x(i1 + 1, 1 + i2, i3)) * (0.5_wp)
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
      end do
      do i2 = n - (1), n - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -1 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + -1 + 1))
        tt(1) = tt(1) + (x(i1 + 1, -1 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + -1 + 1))
        tt(0) = tt(0) + (x(i1 + 0, 0 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 0 + 1))
        tt(1) = tt(1) + (x(i1 + 1, 0 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 0 + 1))
        tt(0) = tt(0) + (x(i1 + 0, 1 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 1 + 1))
        tt(1) = tt(1) + (x(i1 + 1, 1 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 1 + 1))
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (2)) * (2), ndat1 - (1), 1
      do i2 = 0,  -(-1) - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -1 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + -1 + 1))
        tt(0) = tt(0) + (x(i1 + 0, 0 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 0 + 1))
        tt(0) = tt(0) + (x(i1 + 0, 1 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 1 + 1))
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 =  -(-1), n - (1) - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.5_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.5_wp)
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 = n - (1), n - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -1 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + -1 + 1))
        tt(0) = tt(0) + (x(i1 + 0, 0 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 0 + 1))
        tt(0) = tt(0) + (x(i1 + 0, 1 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 1 + 1))
        y(i1 + 0, i2, i3) = tt(0)
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_np_201_u2_0_true_true_true
SUBROUTINE d_poisson2_np_201_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_np_201_u1_2_true_false_true_cost
SUBROUTINE d_poisson2_np_201_a_u4_0_true_false_true(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -1
  integer(kind=4), parameter :: upfil = 1
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(0:8) :: poisson2_fil = (/ &
-1.5_wp, &
2.0_wp, &
-0.5_wp, &
-0.5_wp, &
0.0_wp, &
0.5_wp, &
0.5_wp, &
-2.0_wp, &
1.5_wp /)
  real(kind=8), parameter, dimension(-1:1) :: poisson2_1_fil = (/ &
-0.5_wp, &
0.0_wp, &
0.5_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  integer(kind=4), dimension(-1 - (1):1 - (-1) - (1)) :: mod_arr
  do l = -1 - (1), 1 - (-1) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1, tt2, tt3)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 = 0,  -(-1) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -1 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + -1 + 1))
        tt1 = tt1 + (x(i1 + 1, -1 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + -1 + 1))
        tt2 = tt2 + (x(i1 + 2, -1 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + -1 + 1))
        tt3 = tt3 + (x(i1 + 3, -1 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + -1 + 1))
        tt0 = tt0 + (x(i1 + 0, 0 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 0 + 1))
        tt1 = tt1 + (x(i1 + 1, 0 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 0 + 1))
        tt2 = tt2 + (x(i1 + 2, 0 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 0 + 1))
        tt3 = tt3 + (x(i1 + 3, 0 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 0 + 1))
        tt0 = tt0 + (x(i1 + 0, 1 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 1 + 1))
        tt1 = tt1 + (x(i1 + 1, 1 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 1 + 1))
        tt2 = tt2 + (x(i1 + 2, 1 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 1 + 1))
        tt3 = tt3 + (x(i1 + 3, 1 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 1 + 1))
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
      do i2 =  -(-1), n - (1) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.5_wp)
        tt1 = tt1 + (x(i1 + 1, -1 + i2, i3)) * (-0.5_wp)
        tt2 = tt2 + (x(i1 + 2, -1 + i2, i3)) * (-0.5_wp)
        tt3 = tt3 + (x(i1 + 3, -1 + i2, i3)) * (-0.5_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt1 = tt1 + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt2 = tt2 + (x(i1 + 2, 0 + i2, i3)) * (0.0_wp)
        tt3 = tt3 + (x(i1 + 3, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.5_wp)
        tt1 = tt1 + (x(i1 + 1, 1 + i2, i3)) * (0.5_wp)
        tt2 = tt2 + (x(i1 + 2, 1 + i2, i3)) * (0.5_wp)
        tt3 = tt3 + (x(i1 + 3, 1 + i2, i3)) * (0.5_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
      do i2 = n - (1), n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -1 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + -1 + 1))
        tt1 = tt1 + (x(i1 + 1, -1 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + -1 + 1))
        tt2 = tt2 + (x(i1 + 2, -1 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + -1 + 1))
        tt3 = tt3 + (x(i1 + 3, -1 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + -1 + 1))
        tt0 = tt0 + (x(i1 + 0, 0 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 0 + 1))
        tt1 = tt1 + (x(i1 + 1, 0 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 0 + 1))
        tt2 = tt2 + (x(i1 + 2, 0 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 0 + 1))
        tt3 = tt3 + (x(i1 + 3, 0 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 0 + 1))
        tt0 = tt0 + (x(i1 + 0, 1 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 1 + 1))
        tt1 = tt1 + (x(i1 + 1, 1 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 1 + 1))
        tt2 = tt2 + (x(i1 + 2, 1 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 1 + 1))
        tt3 = tt3 + (x(i1 + 3, 1 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 1 + 1))
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 = 0,  -(-1) - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -1 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + -1 + 1))
        tt0 = tt0 + (x(i1 + 0, 0 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 0 + 1))
        tt0 = tt0 + (x(i1 + 0, 1 + 1, i3)) * (poisson2_fil((i2 - (1)) * (3) + 3 + 1 + 1))
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-1), n - (1) - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.5_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.5_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (1), n - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -1 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + -1 + 1))
        tt0 = tt0 + (x(i1 + 0, 0 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 0 + 1))
        tt0 = tt0 + (x(i1 + 0, 1 + -1 + n - (1), i3)) * (poisson2_fil((i2 + 1 - (n) + 1) * (3) + 3 + 1 + 1))
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson2_np_201_a_u4_0_true_false_true
SUBROUTINE d_poisson2_np_201_a_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (3)) * (ndat_t)
END SUBROUTINE d_poisson2_np_201_a_u1_2_true_false_true_cost
SUBROUTINE d_s0s0_1d_poisson2_cost(d, idim, n, bc, x, y, a, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: d
  integer(kind=4), intent(in) :: idim
  integer(kind=4), intent(in), dimension(0:d - (1)) :: n
  integer(kind=4), intent(in) :: bc
  real(kind=8), intent(in), dimension(*) :: x
  real(kind=8), intent(out), dimension(*) :: y
  real(kind=8), intent(in) :: a
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: i
  integer(kind=4) :: ndat_left
  integer(kind=4) :: ndat_right
  integer(kind=4) :: c
  if (idim == 0) then
    ndat_right = 1
    do i = 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson2_p_10_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson2_p_10_a_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson2_fg_10_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson2_fg_10_a_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson2_fs_10_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson2_fs_10_a_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-2)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson2_np_10_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson2_np_10_a_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      end select
  else if (idim == d - (1)) then
    ndat_left = 1
    do i = 0, d - (2), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson2_p_01_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson2_p_01_a_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      case (1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson2_fg_01_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson2_fg_01_a_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      case (-1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson2_fs_01_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson2_fs_01_a_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      case (-2)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson2_np_01_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson2_np_01_a_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      end select
  else
    ndat_left = 1
    ndat_right = 1
    do i = 0, idim - (1), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    do i = idim + 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson2_p_201_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson2_p_201_a_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson2_fg_201_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson2_fg_201_a_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson2_fs_201_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson2_fs_201_a_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-2)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson2_np_201_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson2_np_201_a_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      end select
  end if
END SUBROUTINE d_s0s0_1d_poisson2_cost
SUBROUTINE d_s0s0_1d_poisson2(d, idim, n, bc, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: d
  integer(kind=4), intent(in) :: idim
  integer(kind=4), intent(in), dimension(0:d - (1)) :: n
  integer(kind=4), intent(in) :: bc
  real(kind=8), intent(in), dimension(*) :: x
  real(kind=8), intent(out), dimension(*) :: y
  real(kind=8), intent(in) :: a
  integer(kind=4) :: i
  integer(kind=4) :: ndat_left
  integer(kind=4) :: ndat_right
  if (idim == 0) then
    ndat_right = 1
    do i = 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          call d_poisson2_p_10_u2_1_false_false_true(n(idim), ndat_right, x, y)
        else
          call d_poisson2_p_10_a_u2_1_true_true_true(n(idim), ndat_right, x, y, a)
        end if
      case (1)
        if (a == 1.0_wp) then
          call d_poisson2_fg_10_u5_1_false_false_false(n(idim), ndat_right, x, y)
        else
          call d_poisson2_fg_10_a_u2_1_false_true_true(n(idim), ndat_right, x, y, a)
        end if
      case (-1)
        if (a == 1.0_wp) then
          call d_poisson2_fs_10_u1_1_false_true_true(n(idim), ndat_right, x, y)
        else
          call d_poisson2_fs_10_a_u2_1_false_true_true(n(idim), ndat_right, x, y, a)
        end if
      case (-2)
        if (a == 1.0_wp) then
          call d_poisson2_np_10_u2_1_true_false_false(n(idim), ndat_right, x, y)
        else
          call d_poisson2_np_10_a_u1_1_true_true_true(n(idim), ndat_right, x, y, a)
        end if
      end select
  else if (idim == d - (1)) then
    ndat_left = 1
    do i = 0, d - (2), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          call d_poisson2_p_01_u4_0_false_true_false(ndat_left, n(idim), x, y)
        else
          call d_poisson2_p_01_a_u5_0_false_true_true(ndat_left, n(idim), x, y, a)
        end if
      case (1)
        if (a == 1.0_wp) then
          call d_poisson2_fg_01_u2_0_false_true_false(ndat_left, n(idim), x, y)
        else
          call d_poisson2_fg_01_a_u4_0_false_false_true(ndat_left, n(idim), x, y, a)
        end if
      case (-1)
        if (a == 1.0_wp) then
          call d_poisson2_fs_01_u5_0_false_false_false(ndat_left, n(idim), x, y)
        else
          call d_poisson2_fs_01_a_u4_0_false_false_true(ndat_left, n(idim), x, y, a)
        end if
      case (-2)
        if (a == 1.0_wp) then
          call d_poisson2_np_01_u4_0_false_false_true(ndat_left, n(idim), x, y)
        else
          call d_poisson2_np_01_a_u4_0_true_false_false(ndat_left, n(idim), x, y, a)
        end if
      end select
  else
    ndat_left = 1
    ndat_right = 1
    do i = 0, idim - (1), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    do i = idim + 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          call d_poisson2_p_201_u4_0_true_true_false(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson2_p_201_a_u4_0_true_true_true(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      case (1)
        if (a == 1.0_wp) then
          call d_poisson2_fg_201_u4_0_false_true_true(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson2_fg_201_a_u4_0_false_false_true(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      case (-1)
        if (a == 1.0_wp) then
          call d_poisson2_fs_201_u4_0_false_false_false(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson2_fs_201_a_u4_0_false_false_true(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      case (-2)
        if (a == 1.0_wp) then
          call d_poisson2_np_201_u2_0_true_true_true(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson2_np_201_a_u4_0_true_false_true(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      end select
  end if
END SUBROUTINE d_s0s0_1d_poisson2
SUBROUTINE d_poisson4_p_10_u2_1_false_false_true(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
!$omp parallel  default(shared) private(i1, i2, tt0, tt1)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 = 0,  -(-2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt0 = tt0 + (x(-2 + i1 - (((i1 + -2 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.08333333333333333_wp)
      tt1 = tt1 + (x(-2 + i1 - (((i1 + -2 + (n) * (2)) / (n) - (2)) * (n)), i2 + 1)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(-1 + i1 - (((i1 + -1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (-0.6666666666666666_wp)
      tt1 = tt1 + (x(-1 + i1 - (((i1 + -1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 1)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(0 + i1 - (((i1 + 0 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.0_wp)
      tt1 = tt1 + (x(0 + i1 - (((i1 + 0 + (n) * (2)) / (n) - (2)) * (n)), i2 + 1)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1 - (((i1 + 1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.6666666666666666_wp)
      tt1 = tt1 + (x(1 + i1 - (((i1 + 1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 1)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(2 + i1 - (((i1 + 2 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (-0.08333333333333333_wp)
      tt1 = tt1 + (x(2 + i1 - (((i1 + 2 + (n) * (2)) / (n) - (2)) * (n)), i2 + 1)) * (-0.08333333333333333_wp)
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
    end do
    do i1 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt0 = tt0 + (x(-2 + i1, i2 + 0)) * (0.08333333333333333_wp)
      tt1 = tt1 + (x(-2 + i1, i2 + 1)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(-1 + i1, i2 + 0)) * (-0.6666666666666666_wp)
      tt1 = tt1 + (x(-1 + i1, i2 + 1)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt1 = tt1 + (x(0 + i1, i2 + 1)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1, i2 + 0)) * (0.6666666666666666_wp)
      tt1 = tt1 + (x(1 + i1, i2 + 1)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(2 + i1, i2 + 0)) * (-0.08333333333333333_wp)
      tt1 = tt1 + (x(2 + i1, i2 + 1)) * (-0.08333333333333333_wp)
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
    end do
    do i1 = n - (2), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt0 = tt0 + (x(-2 + i1 - (((i1 + -2 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.08333333333333333_wp)
      tt1 = tt1 + (x(-2 + i1 - (((i1 + -2 + (n) * (2)) / (n) - (2)) * (n)), i2 + 1)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(-1 + i1 - (((i1 + -1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (-0.6666666666666666_wp)
      tt1 = tt1 + (x(-1 + i1 - (((i1 + -1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 1)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(0 + i1 - (((i1 + 0 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.0_wp)
      tt1 = tt1 + (x(0 + i1 - (((i1 + 0 + (n) * (2)) / (n) - (2)) * (n)), i2 + 1)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1 - (((i1 + 1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.6666666666666666_wp)
      tt1 = tt1 + (x(1 + i1 - (((i1 + 1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 1)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(2 + i1 - (((i1 + 2 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (-0.08333333333333333_wp)
      tt1 = tt1 + (x(2 + i1 - (((i1 + 2 + (n) * (2)) / (n) - (2)) * (n)), i2 + 1)) * (-0.08333333333333333_wp)
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 = 0,  -(-2) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(-2 + i1 - (((i1 + -2 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(-1 + i1 - (((i1 + -1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(0 + i1 - (((i1 + 0 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1 - (((i1 + 1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(2 + i1 - (((i1 + 2 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (-0.08333333333333333_wp)
      y(i1, i2 + 0) = tt0
    end do
    do i1 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(-2 + i1, i2 + 0)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(-1 + i1, i2 + 0)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1, i2 + 0)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(2 + i1, i2 + 0)) * (-0.08333333333333333_wp)
      y(i1, i2 + 0) = tt0
    end do
    do i1 = n - (2), n - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(-2 + i1 - (((i1 + -2 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(-1 + i1 - (((i1 + -1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(0 + i1 - (((i1 + 0 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1 - (((i1 + 1 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(2 + i1 - (((i1 + 2 + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (-0.08333333333333333_wp)
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_p_10_u2_1_false_false_true
SUBROUTINE d_poisson4_p_10_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_p_10_u1_1_true_false_true_cost
SUBROUTINE d_poisson4_p_10_a_u3_1_true_false_true(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  integer(kind=4), dimension(-2 - (2):2 - (-2) - (1)) :: mod_arr
  do l = -2 - (2), 2 - (-2) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2)
!$omp do 
  do i2 = 0, ndat - (3), 3
    do i1 = 0,  -(-2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt0 = tt0 + (x(mod_arr(-2 + i1), i2 + 0)) * (0.08333333333333333_wp)
      tt1 = tt1 + (x(mod_arr(-2 + i1), i2 + 1)) * (0.08333333333333333_wp)
      tt2 = tt2 + (x(mod_arr(-2 + i1), i2 + 2)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(mod_arr(-1 + i1), i2 + 0)) * (-0.6666666666666666_wp)
      tt1 = tt1 + (x(mod_arr(-1 + i1), i2 + 1)) * (-0.6666666666666666_wp)
      tt2 = tt2 + (x(mod_arr(-1 + i1), i2 + 2)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(mod_arr(0 + i1), i2 + 0)) * (0.0_wp)
      tt1 = tt1 + (x(mod_arr(0 + i1), i2 + 1)) * (0.0_wp)
      tt2 = tt2 + (x(mod_arr(0 + i1), i2 + 2)) * (0.0_wp)
      tt0 = tt0 + (x(mod_arr(1 + i1), i2 + 0)) * (0.6666666666666666_wp)
      tt1 = tt1 + (x(mod_arr(1 + i1), i2 + 1)) * (0.6666666666666666_wp)
      tt2 = tt2 + (x(mod_arr(1 + i1), i2 + 2)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(mod_arr(2 + i1), i2 + 0)) * (-0.08333333333333333_wp)
      tt1 = tt1 + (x(mod_arr(2 + i1), i2 + 1)) * (-0.08333333333333333_wp)
      tt2 = tt2 + (x(mod_arr(2 + i1), i2 + 2)) * (-0.08333333333333333_wp)
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
      tt1 = (tt1) * (a)
      y(i1, i2 + 1) = tt1
      tt2 = (tt2) * (a)
      y(i1, i2 + 2) = tt2
    end do
    do i1 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt0 = tt0 + (x(-2 + i1, i2 + 0)) * (0.08333333333333333_wp)
      tt1 = tt1 + (x(-2 + i1, i2 + 1)) * (0.08333333333333333_wp)
      tt2 = tt2 + (x(-2 + i1, i2 + 2)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(-1 + i1, i2 + 0)) * (-0.6666666666666666_wp)
      tt1 = tt1 + (x(-1 + i1, i2 + 1)) * (-0.6666666666666666_wp)
      tt2 = tt2 + (x(-1 + i1, i2 + 2)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt1 = tt1 + (x(0 + i1, i2 + 1)) * (0.0_wp)
      tt2 = tt2 + (x(0 + i1, i2 + 2)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1, i2 + 0)) * (0.6666666666666666_wp)
      tt1 = tt1 + (x(1 + i1, i2 + 1)) * (0.6666666666666666_wp)
      tt2 = tt2 + (x(1 + i1, i2 + 2)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(2 + i1, i2 + 0)) * (-0.08333333333333333_wp)
      tt1 = tt1 + (x(2 + i1, i2 + 1)) * (-0.08333333333333333_wp)
      tt2 = tt2 + (x(2 + i1, i2 + 2)) * (-0.08333333333333333_wp)
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
      tt1 = (tt1) * (a)
      y(i1, i2 + 1) = tt1
      tt2 = (tt2) * (a)
      y(i1, i2 + 2) = tt2
    end do
    do i1 = n - (2), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt0 = tt0 + (x(mod_arr(-2 + i1 - (n)), i2 + 0)) * (0.08333333333333333_wp)
      tt1 = tt1 + (x(mod_arr(-2 + i1 - (n)), i2 + 1)) * (0.08333333333333333_wp)
      tt2 = tt2 + (x(mod_arr(-2 + i1 - (n)), i2 + 2)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(mod_arr(-1 + i1 - (n)), i2 + 0)) * (-0.6666666666666666_wp)
      tt1 = tt1 + (x(mod_arr(-1 + i1 - (n)), i2 + 1)) * (-0.6666666666666666_wp)
      tt2 = tt2 + (x(mod_arr(-1 + i1 - (n)), i2 + 2)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(mod_arr(0 + i1 - (n)), i2 + 0)) * (0.0_wp)
      tt1 = tt1 + (x(mod_arr(0 + i1 - (n)), i2 + 1)) * (0.0_wp)
      tt2 = tt2 + (x(mod_arr(0 + i1 - (n)), i2 + 2)) * (0.0_wp)
      tt0 = tt0 + (x(mod_arr(1 + i1 - (n)), i2 + 0)) * (0.6666666666666666_wp)
      tt1 = tt1 + (x(mod_arr(1 + i1 - (n)), i2 + 1)) * (0.6666666666666666_wp)
      tt2 = tt2 + (x(mod_arr(1 + i1 - (n)), i2 + 2)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(mod_arr(2 + i1 - (n)), i2 + 0)) * (-0.08333333333333333_wp)
      tt1 = tt1 + (x(mod_arr(2 + i1 - (n)), i2 + 1)) * (-0.08333333333333333_wp)
      tt2 = tt2 + (x(mod_arr(2 + i1 - (n)), i2 + 2)) * (-0.08333333333333333_wp)
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
      tt1 = (tt1) * (a)
      y(i1, i2 + 1) = tt1
      tt2 = (tt2) * (a)
      y(i1, i2 + 2) = tt2
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (3)) * (3), ndat - (1), 1
    do i1 = 0,  -(-2) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(mod_arr(-2 + i1), i2 + 0)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(mod_arr(-1 + i1), i2 + 0)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(mod_arr(0 + i1), i2 + 0)) * (0.0_wp)
      tt0 = tt0 + (x(mod_arr(1 + i1), i2 + 0)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(mod_arr(2 + i1), i2 + 0)) * (-0.08333333333333333_wp)
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
    end do
    do i1 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(-2 + i1, i2 + 0)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(-1 + i1, i2 + 0)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1, i2 + 0)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(2 + i1, i2 + 0)) * (-0.08333333333333333_wp)
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
    end do
    do i1 = n - (2), n - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(mod_arr(-2 + i1 - (n)), i2 + 0)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(mod_arr(-1 + i1 - (n)), i2 + 0)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(mod_arr(0 + i1 - (n)), i2 + 0)) * (0.0_wp)
      tt0 = tt0 + (x(mod_arr(1 + i1 - (n)), i2 + 0)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(mod_arr(2 + i1 - (n)), i2 + 0)) * (-0.08333333333333333_wp)
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_p_10_a_u3_1_true_false_true
SUBROUTINE d_poisson4_p_10_a_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_p_10_a_u1_1_true_false_true_cost
SUBROUTINE d_poisson4_p_01_u4_0_false_false_false(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2, tt3)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 = 0,  -(-2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = -2, 2, 1
        tt0 = tt0 + (x(i1 + 0, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson4_2_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson4_2_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson4_2_fil(l))
        tt3 = tt3 + (x(i1 + 3, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson4_2_fil(l))
      end do
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
      y(i1 + 3, i2) = tt3
    end do
    do i2 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = -2, 2, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson4_2_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson4_2_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2)) * (poisson4_2_fil(l))
        tt3 = tt3 + (x(i1 + 3, l + i2)) * (poisson4_2_fil(l))
      end do
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
      y(i1 + 3, i2) = tt3
    end do
    do i2 = n - (2), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = -2, 2, 1
        tt0 = tt0 + (x(i1 + 0, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson4_2_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson4_2_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson4_2_fil(l))
        tt3 = tt3 + (x(i1 + 3, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson4_2_fil(l))
      end do
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
      y(i1 + 3, i2) = tt3
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 = 0,  -(-2) - (1), 1
      tt0 = 0.0_wp
      do l = -2, 2, 1
        tt0 = tt0 + (x(i1 + 0, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson4_2_fil(l))
      end do
      y(i1 + 0, i2) = tt0
    end do
    do i2 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      do l = -2, 2, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson4_2_fil(l))
      end do
      y(i1 + 0, i2) = tt0
    end do
    do i2 = n - (2), n - (1), 1
      tt0 = 0.0_wp
      do l = -2, 2, 1
        tt0 = tt0 + (x(i1 + 0, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson4_2_fil(l))
      end do
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_p_01_u4_0_false_false_false
SUBROUTINE d_poisson4_p_01_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_p_01_u1_0_true_false_true_cost
SUBROUTINE d_poisson4_p_01_a_u4_0_true_false_true(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  integer(kind=4), dimension(-2 - (2):2 - (-2) - (1)) :: mod_arr
  do l = -2 - (2), 2 - (-2) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2, tt3)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 = 0,  -(-2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, mod_arr(-2 + i2))) * (0.08333333333333333_wp)
      tt1 = tt1 + (x(i1 + 1, mod_arr(-2 + i2))) * (0.08333333333333333_wp)
      tt2 = tt2 + (x(i1 + 2, mod_arr(-2 + i2))) * (0.08333333333333333_wp)
      tt3 = tt3 + (x(i1 + 3, mod_arr(-2 + i2))) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(i1 + 0, mod_arr(-1 + i2))) * (-0.6666666666666666_wp)
      tt1 = tt1 + (x(i1 + 1, mod_arr(-1 + i2))) * (-0.6666666666666666_wp)
      tt2 = tt2 + (x(i1 + 2, mod_arr(-1 + i2))) * (-0.6666666666666666_wp)
      tt3 = tt3 + (x(i1 + 3, mod_arr(-1 + i2))) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, mod_arr(0 + i2))) * (0.0_wp)
      tt1 = tt1 + (x(i1 + 1, mod_arr(0 + i2))) * (0.0_wp)
      tt2 = tt2 + (x(i1 + 2, mod_arr(0 + i2))) * (0.0_wp)
      tt3 = tt3 + (x(i1 + 3, mod_arr(0 + i2))) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, mod_arr(1 + i2))) * (0.6666666666666666_wp)
      tt1 = tt1 + (x(i1 + 1, mod_arr(1 + i2))) * (0.6666666666666666_wp)
      tt2 = tt2 + (x(i1 + 2, mod_arr(1 + i2))) * (0.6666666666666666_wp)
      tt3 = tt3 + (x(i1 + 3, mod_arr(1 + i2))) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, mod_arr(2 + i2))) * (-0.08333333333333333_wp)
      tt1 = tt1 + (x(i1 + 1, mod_arr(2 + i2))) * (-0.08333333333333333_wp)
      tt2 = tt2 + (x(i1 + 2, mod_arr(2 + i2))) * (-0.08333333333333333_wp)
      tt3 = tt3 + (x(i1 + 3, mod_arr(2 + i2))) * (-0.08333333333333333_wp)
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
    do i2 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -2 + i2)) * (0.08333333333333333_wp)
      tt1 = tt1 + (x(i1 + 1, -2 + i2)) * (0.08333333333333333_wp)
      tt2 = tt2 + (x(i1 + 2, -2 + i2)) * (0.08333333333333333_wp)
      tt3 = tt3 + (x(i1 + 3, -2 + i2)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.6666666666666666_wp)
      tt1 = tt1 + (x(i1 + 1, -1 + i2)) * (-0.6666666666666666_wp)
      tt2 = tt2 + (x(i1 + 2, -1 + i2)) * (-0.6666666666666666_wp)
      tt3 = tt3 + (x(i1 + 3, -1 + i2)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt1 = tt1 + (x(i1 + 1, 0 + i2)) * (0.0_wp)
      tt2 = tt2 + (x(i1 + 2, 0 + i2)) * (0.0_wp)
      tt3 = tt3 + (x(i1 + 3, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.6666666666666666_wp)
      tt1 = tt1 + (x(i1 + 1, 1 + i2)) * (0.6666666666666666_wp)
      tt2 = tt2 + (x(i1 + 2, 1 + i2)) * (0.6666666666666666_wp)
      tt3 = tt3 + (x(i1 + 3, 1 + i2)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, 2 + i2)) * (-0.08333333333333333_wp)
      tt1 = tt1 + (x(i1 + 1, 2 + i2)) * (-0.08333333333333333_wp)
      tt2 = tt2 + (x(i1 + 2, 2 + i2)) * (-0.08333333333333333_wp)
      tt3 = tt3 + (x(i1 + 3, 2 + i2)) * (-0.08333333333333333_wp)
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
    do i2 = n - (2), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, mod_arr(-2 + i2 - (n)))) * (0.08333333333333333_wp)
      tt1 = tt1 + (x(i1 + 1, mod_arr(-2 + i2 - (n)))) * (0.08333333333333333_wp)
      tt2 = tt2 + (x(i1 + 2, mod_arr(-2 + i2 - (n)))) * (0.08333333333333333_wp)
      tt3 = tt3 + (x(i1 + 3, mod_arr(-2 + i2 - (n)))) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(i1 + 0, mod_arr(-1 + i2 - (n)))) * (-0.6666666666666666_wp)
      tt1 = tt1 + (x(i1 + 1, mod_arr(-1 + i2 - (n)))) * (-0.6666666666666666_wp)
      tt2 = tt2 + (x(i1 + 2, mod_arr(-1 + i2 - (n)))) * (-0.6666666666666666_wp)
      tt3 = tt3 + (x(i1 + 3, mod_arr(-1 + i2 - (n)))) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, mod_arr(0 + i2 - (n)))) * (0.0_wp)
      tt1 = tt1 + (x(i1 + 1, mod_arr(0 + i2 - (n)))) * (0.0_wp)
      tt2 = tt2 + (x(i1 + 2, mod_arr(0 + i2 - (n)))) * (0.0_wp)
      tt3 = tt3 + (x(i1 + 3, mod_arr(0 + i2 - (n)))) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, mod_arr(1 + i2 - (n)))) * (0.6666666666666666_wp)
      tt1 = tt1 + (x(i1 + 1, mod_arr(1 + i2 - (n)))) * (0.6666666666666666_wp)
      tt2 = tt2 + (x(i1 + 2, mod_arr(1 + i2 - (n)))) * (0.6666666666666666_wp)
      tt3 = tt3 + (x(i1 + 3, mod_arr(1 + i2 - (n)))) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, mod_arr(2 + i2 - (n)))) * (-0.08333333333333333_wp)
      tt1 = tt1 + (x(i1 + 1, mod_arr(2 + i2 - (n)))) * (-0.08333333333333333_wp)
      tt2 = tt2 + (x(i1 + 2, mod_arr(2 + i2 - (n)))) * (-0.08333333333333333_wp)
      tt3 = tt3 + (x(i1 + 3, mod_arr(2 + i2 - (n)))) * (-0.08333333333333333_wp)
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 = 0,  -(-2) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, mod_arr(-2 + i2))) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(i1 + 0, mod_arr(-1 + i2))) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, mod_arr(0 + i2))) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, mod_arr(1 + i2))) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, mod_arr(2 + i2))) * (-0.08333333333333333_wp)
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
    do i2 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -2 + i2)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, 2 + i2)) * (-0.08333333333333333_wp)
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
    do i2 = n - (2), n - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, mod_arr(-2 + i2 - (n)))) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(i1 + 0, mod_arr(-1 + i2 - (n)))) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, mod_arr(0 + i2 - (n)))) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, mod_arr(1 + i2 - (n)))) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, mod_arr(2 + i2 - (n)))) * (-0.08333333333333333_wp)
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_p_01_a_u4_0_true_false_true
SUBROUTINE d_poisson4_p_01_a_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_p_01_a_u1_0_true_false_true_cost
SUBROUTINE d_poisson4_p_201_u2_0_true_false_true(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  integer(kind=4), dimension(-2 - (2):2 - (-2) - (1)) :: mod_arr
  do l = -2 - (2), 2 - (-2) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (2), 2
      do i2 = 0,  -(-2) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, mod_arr(-2 + i2), i3)) * (0.08333333333333333_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(-2 + i2), i3)) * (0.08333333333333333_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-1 + i2), i3)) * (-0.6666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(-1 + i2), i3)) * (-0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(1 + i2), i3)) * (0.6666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(1 + i2), i3)) * (0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(2 + i2), i3)) * (-0.08333333333333333_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(2 + i2), i3)) * (-0.08333333333333333_wp)
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
      end do
      do i2 =  -(-2), n - (2) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt1 = tt1 + (x(i1 + 1, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt1 = tt1 + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt1 = tt1 + (x(i1 + 1, 2 + i2, i3)) * (-0.08333333333333333_wp)
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
      end do
      do i2 = n - (2), n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, mod_arr(-2 + i2 - (n)), i3)) * (0.08333333333333333_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(-2 + i2 - (n)), i3)) * (0.08333333333333333_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-1 + i2 - (n)), i3)) * (-0.6666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(-1 + i2 - (n)), i3)) * (-0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(1 + i2 - (n)), i3)) * (0.6666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(1 + i2 - (n)), i3)) * (0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(2 + i2 - (n)), i3)) * (-0.08333333333333333_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(2 + i2 - (n)), i3)) * (-0.08333333333333333_wp)
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (2)) * (2), ndat1 - (1), 1
      do i2 = 0,  -(-2) - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, mod_arr(-2 + i2), i3)) * (0.08333333333333333_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-1 + i2), i3)) * (-0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(1 + i2), i3)) * (0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(2 + i2), i3)) * (-0.08333333333333333_wp)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-2), n - (2) - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, 2 + i2, i3)) * (-0.08333333333333333_wp)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (2), n - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, mod_arr(-2 + i2 - (n)), i3)) * (0.08333333333333333_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-1 + i2 - (n)), i3)) * (-0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(1 + i2 - (n)), i3)) * (0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(2 + i2 - (n)), i3)) * (-0.08333333333333333_wp)
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_p_201_u2_0_true_false_true
SUBROUTINE d_poisson4_p_201_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_p_201_u1_2_true_false_true_cost
SUBROUTINE d_poisson4_p_201_a_u4_0_true_false_true(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  integer(kind=4), dimension(-2 - (2):2 - (-2) - (1)) :: mod_arr
  do l = -2 - (2), 2 - (-2) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1, tt2, tt3)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 = 0,  -(-2) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, mod_arr(-2 + i2), i3)) * (0.08333333333333333_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(-2 + i2), i3)) * (0.08333333333333333_wp)
        tt2 = tt2 + (x(i1 + 2, mod_arr(-2 + i2), i3)) * (0.08333333333333333_wp)
        tt3 = tt3 + (x(i1 + 3, mod_arr(-2 + i2), i3)) * (0.08333333333333333_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-1 + i2), i3)) * (-0.6666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(-1 + i2), i3)) * (-0.6666666666666666_wp)
        tt2 = tt2 + (x(i1 + 2, mod_arr(-1 + i2), i3)) * (-0.6666666666666666_wp)
        tt3 = tt3 + (x(i1 + 3, mod_arr(-1 + i2), i3)) * (-0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt2 = tt2 + (x(i1 + 2, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt3 = tt3 + (x(i1 + 3, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(1 + i2), i3)) * (0.6666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(1 + i2), i3)) * (0.6666666666666666_wp)
        tt2 = tt2 + (x(i1 + 2, mod_arr(1 + i2), i3)) * (0.6666666666666666_wp)
        tt3 = tt3 + (x(i1 + 3, mod_arr(1 + i2), i3)) * (0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(2 + i2), i3)) * (-0.08333333333333333_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(2 + i2), i3)) * (-0.08333333333333333_wp)
        tt2 = tt2 + (x(i1 + 2, mod_arr(2 + i2), i3)) * (-0.08333333333333333_wp)
        tt3 = tt3 + (x(i1 + 3, mod_arr(2 + i2), i3)) * (-0.08333333333333333_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
      do i2 =  -(-2), n - (2) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt1 = tt1 + (x(i1 + 1, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt2 = tt2 + (x(i1 + 2, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt3 = tt3 + (x(i1 + 3, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt2 = tt2 + (x(i1 + 2, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt3 = tt3 + (x(i1 + 3, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt1 = tt1 + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt2 = tt2 + (x(i1 + 2, 0 + i2, i3)) * (0.0_wp)
        tt3 = tt3 + (x(i1 + 3, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt2 = tt2 + (x(i1 + 2, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt3 = tt3 + (x(i1 + 3, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt1 = tt1 + (x(i1 + 1, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt2 = tt2 + (x(i1 + 2, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt3 = tt3 + (x(i1 + 3, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
      do i2 = n - (2), n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, mod_arr(-2 + i2 - (n)), i3)) * (0.08333333333333333_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(-2 + i2 - (n)), i3)) * (0.08333333333333333_wp)
        tt2 = tt2 + (x(i1 + 2, mod_arr(-2 + i2 - (n)), i3)) * (0.08333333333333333_wp)
        tt3 = tt3 + (x(i1 + 3, mod_arr(-2 + i2 - (n)), i3)) * (0.08333333333333333_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-1 + i2 - (n)), i3)) * (-0.6666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(-1 + i2 - (n)), i3)) * (-0.6666666666666666_wp)
        tt2 = tt2 + (x(i1 + 2, mod_arr(-1 + i2 - (n)), i3)) * (-0.6666666666666666_wp)
        tt3 = tt3 + (x(i1 + 3, mod_arr(-1 + i2 - (n)), i3)) * (-0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt2 = tt2 + (x(i1 + 2, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt3 = tt3 + (x(i1 + 3, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(1 + i2 - (n)), i3)) * (0.6666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(1 + i2 - (n)), i3)) * (0.6666666666666666_wp)
        tt2 = tt2 + (x(i1 + 2, mod_arr(1 + i2 - (n)), i3)) * (0.6666666666666666_wp)
        tt3 = tt3 + (x(i1 + 3, mod_arr(1 + i2 - (n)), i3)) * (0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(2 + i2 - (n)), i3)) * (-0.08333333333333333_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(2 + i2 - (n)), i3)) * (-0.08333333333333333_wp)
        tt2 = tt2 + (x(i1 + 2, mod_arr(2 + i2 - (n)), i3)) * (-0.08333333333333333_wp)
        tt3 = tt3 + (x(i1 + 3, mod_arr(2 + i2 - (n)), i3)) * (-0.08333333333333333_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 = 0,  -(-2) - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, mod_arr(-2 + i2), i3)) * (0.08333333333333333_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-1 + i2), i3)) * (-0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(1 + i2), i3)) * (0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(2 + i2), i3)) * (-0.08333333333333333_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-2), n - (2) - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (2), n - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, mod_arr(-2 + i2 - (n)), i3)) * (0.08333333333333333_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-1 + i2 - (n)), i3)) * (-0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(1 + i2 - (n)), i3)) * (0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(2 + i2 - (n)), i3)) * (-0.08333333333333333_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_p_201_a_u4_0_true_false_true
SUBROUTINE d_poisson4_p_201_a_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_p_201_a_u1_2_true_false_true_cost
SUBROUTINE d_poisson4_fg_10_u2_1_false_false_false(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension( -(2):n - (-2) - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
!$omp parallel  default(shared) private(i1, i2, tt0, tt1)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 =  -(2),  -(-2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = max( -(i1), -2), 2, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson4_2_fil(l))
        tt1 = tt1 + (x(l + i1, i2 + 1)) * (poisson4_2_fil(l))
      end do
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
    end do
    do i1 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -2, 2, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson4_2_fil(l))
        tt1 = tt1 + (x(l + i1, i2 + 1)) * (poisson4_2_fil(l))
      end do
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
    end do
    do i1 = n - (2), n - (-2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -2, min(2, n - (1) - (i1)), 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson4_2_fil(l))
        tt1 = tt1 + (x(l + i1, i2 + 1)) * (poisson4_2_fil(l))
      end do
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 =  -(2),  -(-2) - (1), 1
      tt0 = 0.0_wp
      do l = max( -(i1), -2), 2, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson4_2_fil(l))
      end do
      y(i1, i2 + 0) = tt0
    end do
    do i1 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      do l = -2, 2, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson4_2_fil(l))
      end do
      y(i1, i2 + 0) = tt0
    end do
    do i1 = n - (2), n - (-2) - (1), 1
      tt0 = 0.0_wp
      do l = -2, min(2, n - (1) - (i1)), 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson4_2_fil(l))
      end do
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_fg_10_u2_1_false_false_false
SUBROUTINE d_poisson4_fg_10_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_fg_10_u1_1_false_false_true_cost
SUBROUTINE d_poisson4_fg_10_a_u1_1_false_false_true(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension( -(2):n - (-2) - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
!$omp parallel  default(shared) private(i1, i2, tt0)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 =  -(2),  -(-2) - (1), 1
      tt0 = 0.0_wp
      do l = max( -(i1), -2), 2, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson4_2_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
    end do
    do i1 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(-2 + i1, i2 + 0)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(-1 + i1, i2 + 0)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1, i2 + 0)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(2 + i1, i2 + 0)) * (-0.08333333333333333_wp)
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
    end do
    do i1 = n - (2), n - (-2) - (1), 1
      tt0 = 0.0_wp
      do l = -2, min(2, n - (1) - (i1)), 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson4_2_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_fg_10_a_u1_1_false_false_true
SUBROUTINE d_poisson4_fg_10_a_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_fg_10_a_u1_1_false_false_true_cost
SUBROUTINE d_poisson4_fg_01_u4_0_false_false_false(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1),  -(2):n - (-2) - (1)) :: y
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2, tt3)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 =  -(2),  -(-2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = max( -(i2), -2), 2, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson4_2_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson4_2_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2)) * (poisson4_2_fil(l))
        tt3 = tt3 + (x(i1 + 3, l + i2)) * (poisson4_2_fil(l))
      end do
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
      y(i1 + 3, i2) = tt3
    end do
    do i2 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = -2, 2, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson4_2_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson4_2_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2)) * (poisson4_2_fil(l))
        tt3 = tt3 + (x(i1 + 3, l + i2)) * (poisson4_2_fil(l))
      end do
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
      y(i1 + 3, i2) = tt3
    end do
    do i2 = n - (2), n - (-2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = -2, min(2, n - (1) - (i2)), 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson4_2_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson4_2_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2)) * (poisson4_2_fil(l))
        tt3 = tt3 + (x(i1 + 3, l + i2)) * (poisson4_2_fil(l))
      end do
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
      y(i1 + 3, i2) = tt3
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 =  -(2),  -(-2) - (1), 1
      tt0 = 0.0_wp
      do l = max( -(i2), -2), 2, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson4_2_fil(l))
      end do
      y(i1 + 0, i2) = tt0
    end do
    do i2 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      do l = -2, 2, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson4_2_fil(l))
      end do
      y(i1 + 0, i2) = tt0
    end do
    do i2 = n - (2), n - (-2) - (1), 1
      tt0 = 0.0_wp
      do l = -2, min(2, n - (1) - (i2)), 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson4_2_fil(l))
      end do
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_fg_01_u4_0_false_false_false
SUBROUTINE d_poisson4_fg_01_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_fg_01_u1_0_false_false_true_cost
SUBROUTINE d_poisson4_fg_01_a_u4_0_false_true_true(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1),  -(2):n - (-2) - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:3) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 =  -(2),  -(-2) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      do l = max( -(i2), -2), 2, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson4_2_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson4_2_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson4_2_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson4_2_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
    end do
    do i2 =  -(-2), n - (2) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -2 + i2)) * (0.08333333333333333_wp)
      tt(1) = tt(1) + (x(i1 + 1, -2 + i2)) * (0.08333333333333333_wp)
      tt(2) = tt(2) + (x(i1 + 2, -2 + i2)) * (0.08333333333333333_wp)
      tt(3) = tt(3) + (x(i1 + 3, -2 + i2)) * (0.08333333333333333_wp)
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2)) * (-0.6666666666666666_wp)
      tt(1) = tt(1) + (x(i1 + 1, -1 + i2)) * (-0.6666666666666666_wp)
      tt(2) = tt(2) + (x(i1 + 2, -1 + i2)) * (-0.6666666666666666_wp)
      tt(3) = tt(3) + (x(i1 + 3, -1 + i2)) * (-0.6666666666666666_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt(1) = tt(1) + (x(i1 + 1, 0 + i2)) * (0.0_wp)
      tt(2) = tt(2) + (x(i1 + 2, 0 + i2)) * (0.0_wp)
      tt(3) = tt(3) + (x(i1 + 3, 0 + i2)) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2)) * (0.6666666666666666_wp)
      tt(1) = tt(1) + (x(i1 + 1, 1 + i2)) * (0.6666666666666666_wp)
      tt(2) = tt(2) + (x(i1 + 2, 1 + i2)) * (0.6666666666666666_wp)
      tt(3) = tt(3) + (x(i1 + 3, 1 + i2)) * (0.6666666666666666_wp)
      tt(0) = tt(0) + (x(i1 + 0, 2 + i2)) * (-0.08333333333333333_wp)
      tt(1) = tt(1) + (x(i1 + 1, 2 + i2)) * (-0.08333333333333333_wp)
      tt(2) = tt(2) + (x(i1 + 2, 2 + i2)) * (-0.08333333333333333_wp)
      tt(3) = tt(3) + (x(i1 + 3, 2 + i2)) * (-0.08333333333333333_wp)
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
    end do
    do i2 = n - (2), n - (-2) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      do l = -2, min(2, n - (1) - (i2)), 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson4_2_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson4_2_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson4_2_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson4_2_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 =  -(2),  -(-2) - (1), 1
      tt(0) = 0.0_wp
      do l = max( -(i2), -2), 2, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson4_2_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 =  -(-2), n - (2) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -2 + i2)) * (0.08333333333333333_wp)
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2)) * (-0.6666666666666666_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2)) * (0.6666666666666666_wp)
      tt(0) = tt(0) + (x(i1 + 0, 2 + i2)) * (-0.08333333333333333_wp)
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 = n - (2), n - (-2) - (1), 1
      tt(0) = 0.0_wp
      do l = -2, min(2, n - (1) - (i2)), 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson4_2_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_fg_01_a_u4_0_false_true_true
SUBROUTINE d_poisson4_fg_01_a_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_fg_01_a_u1_0_false_false_true_cost
SUBROUTINE d_poisson4_fg_201_u2_0_false_true_false(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1),  -(2):n - (-2) - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
!$omp parallel  default(shared) private(i1, i2, i3, tt)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (2), 2
      do i2 =  -(2),  -(-2) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        do l = max( -(i2), -2), 2, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson4_2_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, l + i2, i3)) * (poisson4_2_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
      end do
      do i2 =  -(-2), n - (2) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        do l = -2, 2, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson4_2_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, l + i2, i3)) * (poisson4_2_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
      end do
      do i2 = n - (2), n - (-2) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        do l = -2, min(2, n - (1) - (i2)), 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson4_2_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, l + i2, i3)) * (poisson4_2_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (2)) * (2), ndat1 - (1), 1
      do i2 =  -(2),  -(-2) - (1), 1
        tt(0) = 0.0_wp
        do l = max( -(i2), -2), 2, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson4_2_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 =  -(-2), n - (2) - (1), 1
        tt(0) = 0.0_wp
        do l = -2, 2, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson4_2_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 = n - (2), n - (-2) - (1), 1
        tt(0) = 0.0_wp
        do l = -2, min(2, n - (1) - (i2)), 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson4_2_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_fg_201_u2_0_false_true_false
SUBROUTINE d_poisson4_fg_201_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_fg_201_u1_2_false_false_true_cost
SUBROUTINE d_poisson4_fg_201_a_u4_0_false_true_false(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1),  -(2):n - (-2) - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8), dimension(0:3) :: tt
!$omp parallel  default(shared) private(i1, i2, i3, tt)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 =  -(2),  -(-2) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        do l = max( -(i2), -2), 2, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson4_2_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, l + i2, i3)) * (poisson4_2_fil(l))
          tt(2) = tt(2) + (x(i1 + 2, l + i2, i3)) * (poisson4_2_fil(l))
          tt(3) = tt(3) + (x(i1 + 3, l + i2, i3)) * (poisson4_2_fil(l))
        end do
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
        tt(1) = (tt(1)) * (a)
        y(i1 + 1, i2, i3) = tt(1)
        tt(2) = (tt(2)) * (a)
        y(i1 + 2, i2, i3) = tt(2)
        tt(3) = (tt(3)) * (a)
        y(i1 + 3, i2, i3) = tt(3)
      end do
      do i2 =  -(-2), n - (2) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        do l = -2, 2, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson4_2_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, l + i2, i3)) * (poisson4_2_fil(l))
          tt(2) = tt(2) + (x(i1 + 2, l + i2, i3)) * (poisson4_2_fil(l))
          tt(3) = tt(3) + (x(i1 + 3, l + i2, i3)) * (poisson4_2_fil(l))
        end do
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
        tt(1) = (tt(1)) * (a)
        y(i1 + 1, i2, i3) = tt(1)
        tt(2) = (tt(2)) * (a)
        y(i1 + 2, i2, i3) = tt(2)
        tt(3) = (tt(3)) * (a)
        y(i1 + 3, i2, i3) = tt(3)
      end do
      do i2 = n - (2), n - (-2) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        do l = -2, min(2, n - (1) - (i2)), 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson4_2_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, l + i2, i3)) * (poisson4_2_fil(l))
          tt(2) = tt(2) + (x(i1 + 2, l + i2, i3)) * (poisson4_2_fil(l))
          tt(3) = tt(3) + (x(i1 + 3, l + i2, i3)) * (poisson4_2_fil(l))
        end do
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
        tt(1) = (tt(1)) * (a)
        y(i1 + 1, i2, i3) = tt(1)
        tt(2) = (tt(2)) * (a)
        y(i1 + 2, i2, i3) = tt(2)
        tt(3) = (tt(3)) * (a)
        y(i1 + 3, i2, i3) = tt(3)
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 =  -(2),  -(-2) - (1), 1
        tt(0) = 0.0_wp
        do l = max( -(i2), -2), 2, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson4_2_fil(l))
        end do
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 =  -(-2), n - (2) - (1), 1
        tt(0) = 0.0_wp
        do l = -2, 2, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson4_2_fil(l))
        end do
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 = n - (2), n - (-2) - (1), 1
        tt(0) = 0.0_wp
        do l = -2, min(2, n - (1) - (i2)), 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson4_2_fil(l))
        end do
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_fg_201_a_u4_0_false_true_false
SUBROUTINE d_poisson4_fg_201_a_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_fg_201_a_u1_2_false_false_true_cost
SUBROUTINE d_poisson4_fs_10_u3_1_false_false_false(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(-2:n + 2 - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2)
!$omp do 
  do i2 = 0, ndat - (3), 3
    do i1 = 0, n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      do l = -2, 2, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson4_2_fil(l))
        tt1 = tt1 + (x(l + i1, i2 + 1)) * (poisson4_2_fil(l))
        tt2 = tt2 + (x(l + i1, i2 + 2)) * (poisson4_2_fil(l))
      end do
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
      y(i1, i2 + 2) = tt2
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (3)) * (3), ndat - (1), 1
    do i1 = 0, n - (1), 1
      tt0 = 0.0_wp
      do l = -2, 2, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson4_2_fil(l))
      end do
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_fs_10_u3_1_false_false_false
SUBROUTINE d_poisson4_fs_10_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_fs_10_u1_1_false_false_true_cost
SUBROUTINE d_poisson4_fs_10_a_u1_1_false_false_true(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(-2:n + 2 - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
!$omp parallel  default(shared) private(i1, i2, tt0)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 = 0, n - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(-2 + i1, i2 + 0)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(-1 + i1, i2 + 0)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1, i2 + 0)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(2 + i1, i2 + 0)) * (-0.08333333333333333_wp)
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_fs_10_a_u1_1_false_false_true
SUBROUTINE d_poisson4_fs_10_a_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_fs_10_a_u1_1_false_false_true_cost
SUBROUTINE d_poisson4_fs_01_u4_0_false_true_false(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), -2:n + 2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:3) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 = 0, n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      do l = -2, 2, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson4_2_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson4_2_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson4_2_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson4_2_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
      y(i1 + 2, i2) = tt(2)
      y(i1 + 3, i2) = tt(3)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 = 0, n - (1), 1
      tt(0) = 0.0_wp
      do l = -2, 2, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson4_2_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_fs_01_u4_0_false_true_false
SUBROUTINE d_poisson4_fs_01_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_fs_01_u1_0_false_false_true_cost
SUBROUTINE d_poisson4_fs_01_a_u4_0_false_false_false(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), -2:n + 2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2, tt3)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 = 0, n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = -2, 2, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson4_2_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson4_2_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2)) * (poisson4_2_fil(l))
        tt3 = tt3 + (x(i1 + 3, l + i2)) * (poisson4_2_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 = 0, n - (1), 1
      tt0 = 0.0_wp
      do l = -2, 2, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson4_2_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_fs_01_a_u4_0_false_false_false
SUBROUTINE d_poisson4_fs_01_a_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_fs_01_a_u1_0_false_false_true_cost
SUBROUTINE d_poisson4_fs_201_u2_0_false_false_false(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), -2:n + 2 - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (2), 2
      do i2 = 0, n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        do l = -2, 2, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson4_2_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson4_2_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (2)) * (2), ndat1 - (1), 1
      do i2 = 0, n - (1), 1
        tt0 = 0.0_wp
        do l = -2, 2, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson4_2_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_fs_201_u2_0_false_false_false
SUBROUTINE d_poisson4_fs_201_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_fs_201_u1_2_false_false_true_cost
SUBROUTINE d_poisson4_fs_201_a_u2_0_false_true_true(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), -2:n + 2 - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
!$omp parallel  default(shared) private(i1, i2, i3, tt)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (2), 2
      do i2 = 0, n - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt(1) = tt(1) + (x(i1 + 1, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt(1) = tt(1) + (x(i1 + 1, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(1) = tt(1) + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt(1) = tt(1) + (x(i1 + 1, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt(0) = tt(0) + (x(i1 + 0, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt(1) = tt(1) + (x(i1 + 1, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
        tt(1) = (tt(1)) * (a)
        y(i1 + 1, i2, i3) = tt(1)
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (2)) * (2), ndat1 - (1), 1
      do i2 = 0, n - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt(0) = tt(0) + (x(i1 + 0, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_fs_201_a_u2_0_false_true_true
SUBROUTINE d_poisson4_fs_201_a_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_fs_201_a_u1_2_false_false_true_cost
SUBROUTINE d_poisson4_np_10_u5_1_false_true_false(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(0:24) :: poisson4_fil = (/ &
-2.0833333333333335_wp, &
4.0_wp, &
-3.0_wp, &
1.3333333333333333_wp, &
-0.25_wp, &
-0.25_wp, &
-0.8333333333333334_wp, &
1.5_wp, &
-0.5_wp, &
0.08333333333333333_wp, &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp, &
-0.08333333333333333_wp, &
0.5_wp, &
-1.5_wp, &
0.8333333333333334_wp, &
0.25_wp, &
0.25_wp, &
-1.3333333333333333_wp, &
3.0_wp, &
-4.0_wp, &
2.0833333333333335_wp /)
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:4) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (5), 5
    do i1 = 0,  -(-2) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      do l = -2, 2, 1
        tt(0) = tt(0) + (x(l + 2, i2 + 0)) * (poisson4_fil((i1 - (2)) * (5) + 10 + l + 2))
        tt(1) = tt(1) + (x(l + 2, i2 + 1)) * (poisson4_fil((i1 - (2)) * (5) + 10 + l + 2))
        tt(2) = tt(2) + (x(l + 2, i2 + 2)) * (poisson4_fil((i1 - (2)) * (5) + 10 + l + 2))
        tt(3) = tt(3) + (x(l + 2, i2 + 3)) * (poisson4_fil((i1 - (2)) * (5) + 10 + l + 2))
        tt(4) = tt(4) + (x(l + 2, i2 + 4)) * (poisson4_fil((i1 - (2)) * (5) + 10 + l + 2))
      end do
      y(i1, i2 + 0) = tt(0)
      y(i1, i2 + 1) = tt(1)
      y(i1, i2 + 2) = tt(2)
      y(i1, i2 + 3) = tt(3)
      y(i1, i2 + 4) = tt(4)
    end do
    do i1 =  -(-2), n - (2) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      do l = -2, 2, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson4_2_fil(l))
        tt(1) = tt(1) + (x(l + i1, i2 + 1)) * (poisson4_2_fil(l))
        tt(2) = tt(2) + (x(l + i1, i2 + 2)) * (poisson4_2_fil(l))
        tt(3) = tt(3) + (x(l + i1, i2 + 3)) * (poisson4_2_fil(l))
        tt(4) = tt(4) + (x(l + i1, i2 + 4)) * (poisson4_2_fil(l))
      end do
      y(i1, i2 + 0) = tt(0)
      y(i1, i2 + 1) = tt(1)
      y(i1, i2 + 2) = tt(2)
      y(i1, i2 + 3) = tt(3)
      y(i1, i2 + 4) = tt(4)
    end do
    do i1 = n - (2), n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      do l = -2, 2, 1
        tt(0) = tt(0) + (x(l + -2 + n - (1), i2 + 0)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + l + 2))
        tt(1) = tt(1) + (x(l + -2 + n - (1), i2 + 1)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + l + 2))
        tt(2) = tt(2) + (x(l + -2 + n - (1), i2 + 2)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + l + 2))
        tt(3) = tt(3) + (x(l + -2 + n - (1), i2 + 3)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + l + 2))
        tt(4) = tt(4) + (x(l + -2 + n - (1), i2 + 4)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + l + 2))
      end do
      y(i1, i2 + 0) = tt(0)
      y(i1, i2 + 1) = tt(1)
      y(i1, i2 + 2) = tt(2)
      y(i1, i2 + 3) = tt(3)
      y(i1, i2 + 4) = tt(4)
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (5)) * (5), ndat - (1), 1
    do i1 = 0,  -(-2) - (1), 1
      tt(0) = 0.0_wp
      do l = -2, 2, 1
        tt(0) = tt(0) + (x(l + 2, i2 + 0)) * (poisson4_fil((i1 - (2)) * (5) + 10 + l + 2))
      end do
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-2), n - (2) - (1), 1
      tt(0) = 0.0_wp
      do l = -2, 2, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson4_2_fil(l))
      end do
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (2), n - (1), 1
      tt(0) = 0.0_wp
      do l = -2, 2, 1
        tt(0) = tt(0) + (x(l + -2 + n - (1), i2 + 0)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + l + 2))
      end do
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_np_10_u5_1_false_true_false
SUBROUTINE d_poisson4_np_10_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_np_10_u1_1_true_false_true_cost
SUBROUTINE d_poisson4_np_10_a_u2_1_true_true_true(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(0:24) :: poisson4_fil = (/ &
-2.0833333333333335_wp, &
4.0_wp, &
-3.0_wp, &
1.3333333333333333_wp, &
-0.25_wp, &
-0.25_wp, &
-0.8333333333333334_wp, &
1.5_wp, &
-0.5_wp, &
0.08333333333333333_wp, &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp, &
-0.08333333333333333_wp, &
0.5_wp, &
-1.5_wp, &
0.8333333333333334_wp, &
0.25_wp, &
0.25_wp, &
-1.3333333333333333_wp, &
3.0_wp, &
-4.0_wp, &
2.0833333333333335_wp /)
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
  integer(kind=4), dimension(-2 - (2):2 - (-2) - (1)) :: mod_arr
  do l = -2 - (2), 2 - (-2) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 = 0,  -(-2) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(0) = tt(0) + (x(-2 + 2, i2 + 0)) * (poisson4_fil((i1 - (2)) * (5) + 10 + -2 + 2))
      tt(1) = tt(1) + (x(-2 + 2, i2 + 1)) * (poisson4_fil((i1 - (2)) * (5) + 10 + -2 + 2))
      tt(0) = tt(0) + (x(-1 + 2, i2 + 0)) * (poisson4_fil((i1 - (2)) * (5) + 10 + -1 + 2))
      tt(1) = tt(1) + (x(-1 + 2, i2 + 1)) * (poisson4_fil((i1 - (2)) * (5) + 10 + -1 + 2))
      tt(0) = tt(0) + (x(0 + 2, i2 + 0)) * (poisson4_fil((i1 - (2)) * (5) + 10 + 0 + 2))
      tt(1) = tt(1) + (x(0 + 2, i2 + 1)) * (poisson4_fil((i1 - (2)) * (5) + 10 + 0 + 2))
      tt(0) = tt(0) + (x(1 + 2, i2 + 0)) * (poisson4_fil((i1 - (2)) * (5) + 10 + 1 + 2))
      tt(1) = tt(1) + (x(1 + 2, i2 + 1)) * (poisson4_fil((i1 - (2)) * (5) + 10 + 1 + 2))
      tt(0) = tt(0) + (x(2 + 2, i2 + 0)) * (poisson4_fil((i1 - (2)) * (5) + 10 + 2 + 2))
      tt(1) = tt(1) + (x(2 + 2, i2 + 1)) * (poisson4_fil((i1 - (2)) * (5) + 10 + 2 + 2))
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
    do i1 =  -(-2), n - (2) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(0) = tt(0) + (x(-2 + i1, i2 + 0)) * (0.08333333333333333_wp)
      tt(1) = tt(1) + (x(-2 + i1, i2 + 1)) * (0.08333333333333333_wp)
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.6666666666666666_wp)
      tt(1) = tt(1) + (x(-1 + i1, i2 + 1)) * (-0.6666666666666666_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(1) = tt(1) + (x(0 + i1, i2 + 1)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.6666666666666666_wp)
      tt(1) = tt(1) + (x(1 + i1, i2 + 1)) * (0.6666666666666666_wp)
      tt(0) = tt(0) + (x(2 + i1, i2 + 0)) * (-0.08333333333333333_wp)
      tt(1) = tt(1) + (x(2 + i1, i2 + 1)) * (-0.08333333333333333_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
    do i1 = n - (2), n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(0) = tt(0) + (x(-2 + -2 + n - (1), i2 + 0)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
      tt(1) = tt(1) + (x(-2 + -2 + n - (1), i2 + 1)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
      tt(0) = tt(0) + (x(-1 + -2 + n - (1), i2 + 0)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
      tt(1) = tt(1) + (x(-1 + -2 + n - (1), i2 + 1)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
      tt(0) = tt(0) + (x(0 + -2 + n - (1), i2 + 0)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
      tt(1) = tt(1) + (x(0 + -2 + n - (1), i2 + 1)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
      tt(0) = tt(0) + (x(1 + -2 + n - (1), i2 + 0)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
      tt(1) = tt(1) + (x(1 + -2 + n - (1), i2 + 1)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
      tt(0) = tt(0) + (x(2 + -2 + n - (1), i2 + 0)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
      tt(1) = tt(1) + (x(2 + -2 + n - (1), i2 + 1)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 = 0,  -(-2) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-2 + 2, i2 + 0)) * (poisson4_fil((i1 - (2)) * (5) + 10 + -2 + 2))
      tt(0) = tt(0) + (x(-1 + 2, i2 + 0)) * (poisson4_fil((i1 - (2)) * (5) + 10 + -1 + 2))
      tt(0) = tt(0) + (x(0 + 2, i2 + 0)) * (poisson4_fil((i1 - (2)) * (5) + 10 + 0 + 2))
      tt(0) = tt(0) + (x(1 + 2, i2 + 0)) * (poisson4_fil((i1 - (2)) * (5) + 10 + 1 + 2))
      tt(0) = tt(0) + (x(2 + 2, i2 + 0)) * (poisson4_fil((i1 - (2)) * (5) + 10 + 2 + 2))
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-2), n - (2) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-2 + i1, i2 + 0)) * (0.08333333333333333_wp)
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.6666666666666666_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.6666666666666666_wp)
      tt(0) = tt(0) + (x(2 + i1, i2 + 0)) * (-0.08333333333333333_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (2), n - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-2 + -2 + n - (1), i2 + 0)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
      tt(0) = tt(0) + (x(-1 + -2 + n - (1), i2 + 0)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
      tt(0) = tt(0) + (x(0 + -2 + n - (1), i2 + 0)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
      tt(0) = tt(0) + (x(1 + -2 + n - (1), i2 + 0)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
      tt(0) = tt(0) + (x(2 + -2 + n - (1), i2 + 0)) * (poisson4_fil((i1 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_np_10_a_u2_1_true_true_true
SUBROUTINE d_poisson4_np_10_a_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_np_10_a_u1_1_true_false_true_cost
SUBROUTINE d_poisson4_np_01_u4_0_true_false_true(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), parameter, dimension(0:24) :: poisson4_fil = (/ &
-2.0833333333333335_wp, &
4.0_wp, &
-3.0_wp, &
1.3333333333333333_wp, &
-0.25_wp, &
-0.25_wp, &
-0.8333333333333334_wp, &
1.5_wp, &
-0.5_wp, &
0.08333333333333333_wp, &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp, &
-0.08333333333333333_wp, &
0.5_wp, &
-1.5_wp, &
0.8333333333333334_wp, &
0.25_wp, &
0.25_wp, &
-1.3333333333333333_wp, &
3.0_wp, &
-4.0_wp, &
2.0833333333333335_wp /)
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  integer(kind=4), dimension(-2 - (2):2 - (-2) - (1)) :: mod_arr
  do l = -2 - (2), 2 - (-2) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2, tt3)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 = 0,  -(-2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
      tt1 = tt1 + (x(i1 + 1, -2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
      tt2 = tt2 + (x(i1 + 2, -2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
      tt3 = tt3 + (x(i1 + 3, -2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
      tt0 = tt0 + (x(i1 + 0, -1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
      tt1 = tt1 + (x(i1 + 1, -1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
      tt2 = tt2 + (x(i1 + 2, -1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
      tt3 = tt3 + (x(i1 + 3, -1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
      tt0 = tt0 + (x(i1 + 0, 0 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
      tt1 = tt1 + (x(i1 + 1, 0 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
      tt2 = tt2 + (x(i1 + 2, 0 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
      tt3 = tt3 + (x(i1 + 3, 0 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
      tt0 = tt0 + (x(i1 + 0, 1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
      tt1 = tt1 + (x(i1 + 1, 1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
      tt2 = tt2 + (x(i1 + 2, 1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
      tt3 = tt3 + (x(i1 + 3, 1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
      tt0 = tt0 + (x(i1 + 0, 2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
      tt1 = tt1 + (x(i1 + 1, 2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
      tt2 = tt2 + (x(i1 + 2, 2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
      tt3 = tt3 + (x(i1 + 3, 2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
      y(i1 + 3, i2) = tt3
    end do
    do i2 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -2 + i2)) * (0.08333333333333333_wp)
      tt1 = tt1 + (x(i1 + 1, -2 + i2)) * (0.08333333333333333_wp)
      tt2 = tt2 + (x(i1 + 2, -2 + i2)) * (0.08333333333333333_wp)
      tt3 = tt3 + (x(i1 + 3, -2 + i2)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.6666666666666666_wp)
      tt1 = tt1 + (x(i1 + 1, -1 + i2)) * (-0.6666666666666666_wp)
      tt2 = tt2 + (x(i1 + 2, -1 + i2)) * (-0.6666666666666666_wp)
      tt3 = tt3 + (x(i1 + 3, -1 + i2)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt1 = tt1 + (x(i1 + 1, 0 + i2)) * (0.0_wp)
      tt2 = tt2 + (x(i1 + 2, 0 + i2)) * (0.0_wp)
      tt3 = tt3 + (x(i1 + 3, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.6666666666666666_wp)
      tt1 = tt1 + (x(i1 + 1, 1 + i2)) * (0.6666666666666666_wp)
      tt2 = tt2 + (x(i1 + 2, 1 + i2)) * (0.6666666666666666_wp)
      tt3 = tt3 + (x(i1 + 3, 1 + i2)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, 2 + i2)) * (-0.08333333333333333_wp)
      tt1 = tt1 + (x(i1 + 1, 2 + i2)) * (-0.08333333333333333_wp)
      tt2 = tt2 + (x(i1 + 2, 2 + i2)) * (-0.08333333333333333_wp)
      tt3 = tt3 + (x(i1 + 3, 2 + i2)) * (-0.08333333333333333_wp)
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
      y(i1 + 3, i2) = tt3
    end do
    do i2 = n - (2), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
      tt1 = tt1 + (x(i1 + 1, -2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
      tt2 = tt2 + (x(i1 + 2, -2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
      tt3 = tt3 + (x(i1 + 3, -2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
      tt0 = tt0 + (x(i1 + 0, -1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
      tt1 = tt1 + (x(i1 + 1, -1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
      tt2 = tt2 + (x(i1 + 2, -1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
      tt3 = tt3 + (x(i1 + 3, -1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
      tt0 = tt0 + (x(i1 + 0, 0 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
      tt1 = tt1 + (x(i1 + 1, 0 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
      tt2 = tt2 + (x(i1 + 2, 0 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
      tt3 = tt3 + (x(i1 + 3, 0 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
      tt0 = tt0 + (x(i1 + 0, 1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
      tt1 = tt1 + (x(i1 + 1, 1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
      tt2 = tt2 + (x(i1 + 2, 1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
      tt3 = tt3 + (x(i1 + 3, 1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
      tt0 = tt0 + (x(i1 + 0, 2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
      tt1 = tt1 + (x(i1 + 1, 2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
      tt2 = tt2 + (x(i1 + 2, 2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
      tt3 = tt3 + (x(i1 + 3, 2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
      y(i1 + 3, i2) = tt3
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 = 0,  -(-2) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
      tt0 = tt0 + (x(i1 + 0, -1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
      tt0 = tt0 + (x(i1 + 0, 0 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
      tt0 = tt0 + (x(i1 + 0, 1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
      tt0 = tt0 + (x(i1 + 0, 2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
      y(i1 + 0, i2) = tt0
    end do
    do i2 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -2 + i2)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, 2 + i2)) * (-0.08333333333333333_wp)
      y(i1 + 0, i2) = tt0
    end do
    do i2 = n - (2), n - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
      tt0 = tt0 + (x(i1 + 0, -1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
      tt0 = tt0 + (x(i1 + 0, 0 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
      tt0 = tt0 + (x(i1 + 0, 1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
      tt0 = tt0 + (x(i1 + 0, 2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_np_01_u4_0_true_false_true
SUBROUTINE d_poisson4_np_01_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_np_01_u1_0_true_false_true_cost
SUBROUTINE d_poisson4_np_01_a_u3_0_false_false_true(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(0:24) :: poisson4_fil = (/ &
-2.0833333333333335_wp, &
4.0_wp, &
-3.0_wp, &
1.3333333333333333_wp, &
-0.25_wp, &
-0.25_wp, &
-0.8333333333333334_wp, &
1.5_wp, &
-0.5_wp, &
0.08333333333333333_wp, &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp, &
-0.08333333333333333_wp, &
0.5_wp, &
-1.5_wp, &
0.8333333333333334_wp, &
0.25_wp, &
0.25_wp, &
-1.3333333333333333_wp, &
3.0_wp, &
-4.0_wp, &
2.0833333333333335_wp /)
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2)
!$omp do 
  do i1 = 0, ndat - (3), 3
    do i2 = 0,  -(-2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
      tt1 = tt1 + (x(i1 + 1, -2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
      tt2 = tt2 + (x(i1 + 2, -2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
      tt0 = tt0 + (x(i1 + 0, -1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
      tt1 = tt1 + (x(i1 + 1, -1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
      tt2 = tt2 + (x(i1 + 2, -1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
      tt0 = tt0 + (x(i1 + 0, 0 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
      tt1 = tt1 + (x(i1 + 1, 0 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
      tt2 = tt2 + (x(i1 + 2, 0 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
      tt0 = tt0 + (x(i1 + 0, 1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
      tt1 = tt1 + (x(i1 + 1, 1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
      tt2 = tt2 + (x(i1 + 2, 1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
      tt0 = tt0 + (x(i1 + 0, 2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
      tt1 = tt1 + (x(i1 + 1, 2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
      tt2 = tt2 + (x(i1 + 2, 2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
    end do
    do i2 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -2 + i2)) * (0.08333333333333333_wp)
      tt1 = tt1 + (x(i1 + 1, -2 + i2)) * (0.08333333333333333_wp)
      tt2 = tt2 + (x(i1 + 2, -2 + i2)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.6666666666666666_wp)
      tt1 = tt1 + (x(i1 + 1, -1 + i2)) * (-0.6666666666666666_wp)
      tt2 = tt2 + (x(i1 + 2, -1 + i2)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt1 = tt1 + (x(i1 + 1, 0 + i2)) * (0.0_wp)
      tt2 = tt2 + (x(i1 + 2, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.6666666666666666_wp)
      tt1 = tt1 + (x(i1 + 1, 1 + i2)) * (0.6666666666666666_wp)
      tt2 = tt2 + (x(i1 + 2, 1 + i2)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, 2 + i2)) * (-0.08333333333333333_wp)
      tt1 = tt1 + (x(i1 + 1, 2 + i2)) * (-0.08333333333333333_wp)
      tt2 = tt2 + (x(i1 + 2, 2 + i2)) * (-0.08333333333333333_wp)
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
    end do
    do i2 = n - (2), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
      tt1 = tt1 + (x(i1 + 1, -2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
      tt2 = tt2 + (x(i1 + 2, -2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
      tt0 = tt0 + (x(i1 + 0, -1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
      tt1 = tt1 + (x(i1 + 1, -1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
      tt2 = tt2 + (x(i1 + 2, -1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
      tt0 = tt0 + (x(i1 + 0, 0 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
      tt1 = tt1 + (x(i1 + 1, 0 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
      tt2 = tt2 + (x(i1 + 2, 0 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
      tt0 = tt0 + (x(i1 + 0, 1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
      tt1 = tt1 + (x(i1 + 1, 1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
      tt2 = tt2 + (x(i1 + 2, 1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
      tt0 = tt0 + (x(i1 + 0, 2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
      tt1 = tt1 + (x(i1 + 1, 2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
      tt2 = tt2 + (x(i1 + 2, 2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (3)) * (3), ndat - (1), 1
    do i2 = 0,  -(-2) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
      tt0 = tt0 + (x(i1 + 0, -1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
      tt0 = tt0 + (x(i1 + 0, 0 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
      tt0 = tt0 + (x(i1 + 0, 1 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
      tt0 = tt0 + (x(i1 + 0, 2 + 2)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
    do i2 =  -(-2), n - (2) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -2 + i2)) * (0.08333333333333333_wp)
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.6666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, 2 + i2)) * (-0.08333333333333333_wp)
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
    do i2 = n - (2), n - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
      tt0 = tt0 + (x(i1 + 0, -1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
      tt0 = tt0 + (x(i1 + 0, 0 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
      tt0 = tt0 + (x(i1 + 0, 1 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
      tt0 = tt0 + (x(i1 + 0, 2 + -2 + n - (1))) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_np_01_a_u3_0_false_false_true
SUBROUTINE d_poisson4_np_01_a_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_np_01_a_u1_0_true_false_true_cost
SUBROUTINE d_poisson4_np_201_u4_0_false_true_true(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(0:24) :: poisson4_fil = (/ &
-2.0833333333333335_wp, &
4.0_wp, &
-3.0_wp, &
1.3333333333333333_wp, &
-0.25_wp, &
-0.25_wp, &
-0.8333333333333334_wp, &
1.5_wp, &
-0.5_wp, &
0.08333333333333333_wp, &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp, &
-0.08333333333333333_wp, &
0.5_wp, &
-1.5_wp, &
0.8333333333333334_wp, &
0.25_wp, &
0.25_wp, &
-1.3333333333333333_wp, &
3.0_wp, &
-4.0_wp, &
2.0833333333333335_wp /)
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8), dimension(0:3) :: tt
!$omp parallel  default(shared) private(i1, i2, i3, tt)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 = 0,  -(-2) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
        tt(1) = tt(1) + (x(i1 + 1, -2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
        tt(2) = tt(2) + (x(i1 + 2, -2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
        tt(3) = tt(3) + (x(i1 + 3, -2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
        tt(0) = tt(0) + (x(i1 + 0, -1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
        tt(1) = tt(1) + (x(i1 + 1, -1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
        tt(2) = tt(2) + (x(i1 + 2, -1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
        tt(3) = tt(3) + (x(i1 + 3, -1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
        tt(0) = tt(0) + (x(i1 + 0, 0 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
        tt(1) = tt(1) + (x(i1 + 1, 0 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
        tt(2) = tt(2) + (x(i1 + 2, 0 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
        tt(3) = tt(3) + (x(i1 + 3, 0 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
        tt(0) = tt(0) + (x(i1 + 0, 1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
        tt(1) = tt(1) + (x(i1 + 1, 1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
        tt(2) = tt(2) + (x(i1 + 2, 1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
        tt(3) = tt(3) + (x(i1 + 3, 1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
        tt(0) = tt(0) + (x(i1 + 0, 2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
        tt(1) = tt(1) + (x(i1 + 1, 2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
        tt(2) = tt(2) + (x(i1 + 2, 2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
        tt(3) = tt(3) + (x(i1 + 3, 2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
      end do
      do i2 =  -(-2), n - (2) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt(1) = tt(1) + (x(i1 + 1, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt(2) = tt(2) + (x(i1 + 2, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt(3) = tt(3) + (x(i1 + 3, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt(1) = tt(1) + (x(i1 + 1, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt(2) = tt(2) + (x(i1 + 2, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt(3) = tt(3) + (x(i1 + 3, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(1) = tt(1) + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt(2) = tt(2) + (x(i1 + 2, 0 + i2, i3)) * (0.0_wp)
        tt(3) = tt(3) + (x(i1 + 3, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt(1) = tt(1) + (x(i1 + 1, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt(2) = tt(2) + (x(i1 + 2, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt(3) = tt(3) + (x(i1 + 3, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt(0) = tt(0) + (x(i1 + 0, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt(1) = tt(1) + (x(i1 + 1, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt(2) = tt(2) + (x(i1 + 2, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt(3) = tt(3) + (x(i1 + 3, 2 + i2, i3)) * (-0.08333333333333333_wp)
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
      end do
      do i2 = n - (2), n - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
        tt(1) = tt(1) + (x(i1 + 1, -2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
        tt(2) = tt(2) + (x(i1 + 2, -2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
        tt(3) = tt(3) + (x(i1 + 3, -2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
        tt(0) = tt(0) + (x(i1 + 0, -1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
        tt(1) = tt(1) + (x(i1 + 1, -1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
        tt(2) = tt(2) + (x(i1 + 2, -1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
        tt(3) = tt(3) + (x(i1 + 3, -1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
        tt(0) = tt(0) + (x(i1 + 0, 0 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
        tt(1) = tt(1) + (x(i1 + 1, 0 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
        tt(2) = tt(2) + (x(i1 + 2, 0 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
        tt(3) = tt(3) + (x(i1 + 3, 0 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
        tt(0) = tt(0) + (x(i1 + 0, 1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
        tt(1) = tt(1) + (x(i1 + 1, 1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
        tt(2) = tt(2) + (x(i1 + 2, 1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
        tt(3) = tt(3) + (x(i1 + 3, 1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
        tt(0) = tt(0) + (x(i1 + 0, 2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
        tt(1) = tt(1) + (x(i1 + 1, 2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
        tt(2) = tt(2) + (x(i1 + 2, 2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
        tt(3) = tt(3) + (x(i1 + 3, 2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 = 0,  -(-2) - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
        tt(0) = tt(0) + (x(i1 + 0, -1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
        tt(0) = tt(0) + (x(i1 + 0, 0 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
        tt(0) = tt(0) + (x(i1 + 0, 1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
        tt(0) = tt(0) + (x(i1 + 0, 2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 =  -(-2), n - (2) - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt(0) = tt(0) + (x(i1 + 0, 2 + i2, i3)) * (-0.08333333333333333_wp)
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 = n - (2), n - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
        tt(0) = tt(0) + (x(i1 + 0, -1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
        tt(0) = tt(0) + (x(i1 + 0, 0 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
        tt(0) = tt(0) + (x(i1 + 0, 1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
        tt(0) = tt(0) + (x(i1 + 0, 2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
        y(i1 + 0, i2, i3) = tt(0)
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_np_201_u4_0_false_true_true
SUBROUTINE d_poisson4_np_201_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_np_201_u1_2_true_false_true_cost
SUBROUTINE d_poisson4_np_201_a_u4_0_false_false_true(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -2
  integer(kind=4), parameter :: upfil = 2
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(0:24) :: poisson4_fil = (/ &
-2.0833333333333335_wp, &
4.0_wp, &
-3.0_wp, &
1.3333333333333333_wp, &
-0.25_wp, &
-0.25_wp, &
-0.8333333333333334_wp, &
1.5_wp, &
-0.5_wp, &
0.08333333333333333_wp, &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp, &
-0.08333333333333333_wp, &
0.5_wp, &
-1.5_wp, &
0.8333333333333334_wp, &
0.25_wp, &
0.25_wp, &
-1.3333333333333333_wp, &
3.0_wp, &
-4.0_wp, &
2.0833333333333335_wp /)
  real(kind=8), parameter, dimension(-2:2) :: poisson4_2_fil = (/ &
0.08333333333333333_wp, &
-0.6666666666666666_wp, &
0.0_wp, &
0.6666666666666666_wp, &
-0.08333333333333333_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1, tt2, tt3)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 = 0,  -(-2) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
        tt1 = tt1 + (x(i1 + 1, -2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
        tt2 = tt2 + (x(i1 + 2, -2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
        tt3 = tt3 + (x(i1 + 3, -2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
        tt0 = tt0 + (x(i1 + 0, -1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
        tt1 = tt1 + (x(i1 + 1, -1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
        tt2 = tt2 + (x(i1 + 2, -1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
        tt3 = tt3 + (x(i1 + 3, -1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
        tt0 = tt0 + (x(i1 + 0, 0 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
        tt1 = tt1 + (x(i1 + 1, 0 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
        tt2 = tt2 + (x(i1 + 2, 0 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
        tt3 = tt3 + (x(i1 + 3, 0 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
        tt0 = tt0 + (x(i1 + 0, 1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
        tt1 = tt1 + (x(i1 + 1, 1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
        tt2 = tt2 + (x(i1 + 2, 1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
        tt3 = tt3 + (x(i1 + 3, 1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
        tt0 = tt0 + (x(i1 + 0, 2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
        tt1 = tt1 + (x(i1 + 1, 2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
        tt2 = tt2 + (x(i1 + 2, 2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
        tt3 = tt3 + (x(i1 + 3, 2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
      do i2 =  -(-2), n - (2) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt1 = tt1 + (x(i1 + 1, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt2 = tt2 + (x(i1 + 2, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt3 = tt3 + (x(i1 + 3, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt2 = tt2 + (x(i1 + 2, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt3 = tt3 + (x(i1 + 3, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt1 = tt1 + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt2 = tt2 + (x(i1 + 2, 0 + i2, i3)) * (0.0_wp)
        tt3 = tt3 + (x(i1 + 3, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt2 = tt2 + (x(i1 + 2, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt3 = tt3 + (x(i1 + 3, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt1 = tt1 + (x(i1 + 1, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt2 = tt2 + (x(i1 + 2, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt3 = tt3 + (x(i1 + 3, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
      do i2 = n - (2), n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
        tt1 = tt1 + (x(i1 + 1, -2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
        tt2 = tt2 + (x(i1 + 2, -2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
        tt3 = tt3 + (x(i1 + 3, -2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
        tt0 = tt0 + (x(i1 + 0, -1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
        tt1 = tt1 + (x(i1 + 1, -1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
        tt2 = tt2 + (x(i1 + 2, -1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
        tt3 = tt3 + (x(i1 + 3, -1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
        tt0 = tt0 + (x(i1 + 0, 0 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
        tt1 = tt1 + (x(i1 + 1, 0 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
        tt2 = tt2 + (x(i1 + 2, 0 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
        tt3 = tt3 + (x(i1 + 3, 0 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
        tt0 = tt0 + (x(i1 + 0, 1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
        tt1 = tt1 + (x(i1 + 1, 1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
        tt2 = tt2 + (x(i1 + 2, 1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
        tt3 = tt3 + (x(i1 + 3, 1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
        tt0 = tt0 + (x(i1 + 0, 2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
        tt1 = tt1 + (x(i1 + 1, 2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
        tt2 = tt2 + (x(i1 + 2, 2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
        tt3 = tt3 + (x(i1 + 3, 2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 = 0,  -(-2) - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -2 + 2))
        tt0 = tt0 + (x(i1 + 0, -1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + -1 + 2))
        tt0 = tt0 + (x(i1 + 0, 0 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 0 + 2))
        tt0 = tt0 + (x(i1 + 0, 1 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 1 + 2))
        tt0 = tt0 + (x(i1 + 0, 2 + 2, i3)) * (poisson4_fil((i2 - (2)) * (5) + 10 + 2 + 2))
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-2), n - (2) - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -2 + i2, i3)) * (0.08333333333333333_wp)
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.6666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, 2 + i2, i3)) * (-0.08333333333333333_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (2), n - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -2 + 2))
        tt0 = tt0 + (x(i1 + 0, -1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + -1 + 2))
        tt0 = tt0 + (x(i1 + 0, 0 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 0 + 2))
        tt0 = tt0 + (x(i1 + 0, 1 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 1 + 2))
        tt0 = tt0 + (x(i1 + 0, 2 + -2 + n - (1), i3)) * (poisson4_fil((i2 + 2 - (n) + 1) * (5) + 10 + 2 + 2))
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson4_np_201_a_u4_0_false_false_true
SUBROUTINE d_poisson4_np_201_a_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (5)) * (ndat_t)
END SUBROUTINE d_poisson4_np_201_a_u1_2_true_false_true_cost
SUBROUTINE d_s0s0_1d_poisson4_cost(d, idim, n, bc, x, y, a, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: d
  integer(kind=4), intent(in) :: idim
  integer(kind=4), intent(in), dimension(0:d - (1)) :: n
  integer(kind=4), intent(in) :: bc
  real(kind=8), intent(in), dimension(*) :: x
  real(kind=8), intent(out), dimension(*) :: y
  real(kind=8), intent(in) :: a
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: i
  integer(kind=4) :: ndat_left
  integer(kind=4) :: ndat_right
  integer(kind=4) :: c
  if (idim == 0) then
    ndat_right = 1
    do i = 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson4_p_10_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson4_p_10_a_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson4_fg_10_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson4_fg_10_a_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson4_fs_10_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson4_fs_10_a_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-2)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson4_np_10_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson4_np_10_a_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      end select
  else if (idim == d - (1)) then
    ndat_left = 1
    do i = 0, d - (2), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson4_p_01_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson4_p_01_a_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      case (1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson4_fg_01_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson4_fg_01_a_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      case (-1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson4_fs_01_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson4_fs_01_a_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      case (-2)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson4_np_01_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson4_np_01_a_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      end select
  else
    ndat_left = 1
    ndat_right = 1
    do i = 0, idim - (1), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    do i = idim + 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson4_p_201_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson4_p_201_a_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson4_fg_201_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson4_fg_201_a_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson4_fs_201_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson4_fs_201_a_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-2)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson4_np_201_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson4_np_201_a_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      end select
  end if
END SUBROUTINE d_s0s0_1d_poisson4_cost
SUBROUTINE d_s0s0_1d_poisson4(d, idim, n, bc, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: d
  integer(kind=4), intent(in) :: idim
  integer(kind=4), intent(in), dimension(0:d - (1)) :: n
  integer(kind=4), intent(in) :: bc
  real(kind=8), intent(in), dimension(*) :: x
  real(kind=8), intent(out), dimension(*) :: y
  real(kind=8), intent(in) :: a
  integer(kind=4) :: i
  integer(kind=4) :: ndat_left
  integer(kind=4) :: ndat_right
  if (idim == 0) then
    ndat_right = 1
    do i = 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          call d_poisson4_p_10_u2_1_false_false_true(n(idim), ndat_right, x, y)
        else
          call d_poisson4_p_10_a_u3_1_true_false_true(n(idim), ndat_right, x, y, a)
        end if
      case (1)
        if (a == 1.0_wp) then
          call d_poisson4_fg_10_u2_1_false_false_false(n(idim), ndat_right, x, y)
        else
          call d_poisson4_fg_10_a_u1_1_false_false_true(n(idim), ndat_right, x, y, a)
        end if
      case (-1)
        if (a == 1.0_wp) then
          call d_poisson4_fs_10_u3_1_false_false_false(n(idim), ndat_right, x, y)
        else
          call d_poisson4_fs_10_a_u1_1_false_false_true(n(idim), ndat_right, x, y, a)
        end if
      case (-2)
        if (a == 1.0_wp) then
          call d_poisson4_np_10_u5_1_false_true_false(n(idim), ndat_right, x, y)
        else
          call d_poisson4_np_10_a_u2_1_true_true_true(n(idim), ndat_right, x, y, a)
        end if
      end select
  else if (idim == d - (1)) then
    ndat_left = 1
    do i = 0, d - (2), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          call d_poisson4_p_01_u4_0_false_false_false(ndat_left, n(idim), x, y)
        else
          call d_poisson4_p_01_a_u4_0_true_false_true(ndat_left, n(idim), x, y, a)
        end if
      case (1)
        if (a == 1.0_wp) then
          call d_poisson4_fg_01_u4_0_false_false_false(ndat_left, n(idim), x, y)
        else
          call d_poisson4_fg_01_a_u4_0_false_true_true(ndat_left, n(idim), x, y, a)
        end if
      case (-1)
        if (a == 1.0_wp) then
          call d_poisson4_fs_01_u4_0_false_true_false(ndat_left, n(idim), x, y)
        else
          call d_poisson4_fs_01_a_u4_0_false_false_false(ndat_left, n(idim), x, y, a)
        end if
      case (-2)
        if (a == 1.0_wp) then
          call d_poisson4_np_01_u4_0_true_false_true(ndat_left, n(idim), x, y)
        else
          call d_poisson4_np_01_a_u3_0_false_false_true(ndat_left, n(idim), x, y, a)
        end if
      end select
  else
    ndat_left = 1
    ndat_right = 1
    do i = 0, idim - (1), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    do i = idim + 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          call d_poisson4_p_201_u2_0_true_false_true(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson4_p_201_a_u4_0_true_false_true(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      case (1)
        if (a == 1.0_wp) then
          call d_poisson4_fg_201_u2_0_false_true_false(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson4_fg_201_a_u4_0_false_true_false(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      case (-1)
        if (a == 1.0_wp) then
          call d_poisson4_fs_201_u2_0_false_false_false(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson4_fs_201_a_u2_0_false_true_true(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      case (-2)
        if (a == 1.0_wp) then
          call d_poisson4_np_201_u4_0_false_true_true(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson4_np_201_a_u4_0_false_false_true(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      end select
  end if
END SUBROUTINE d_s0s0_1d_poisson4
SUBROUTINE d_poisson6_p_10_u1_1_true_false_true(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  integer(kind=4), dimension(-3 - (3):3 - (-3) - (1)) :: mod_arr
  do l = -3 - (3), 3 - (-3) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt0)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 = 0,  -(-3) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(mod_arr(-3 + i1), i2 + 0)) * (-0.016666666666666666_wp)
      tt0 = tt0 + (x(mod_arr(-2 + i1), i2 + 0)) * (0.15_wp)
      tt0 = tt0 + (x(mod_arr(-1 + i1), i2 + 0)) * (-0.75_wp)
      tt0 = tt0 + (x(mod_arr(0 + i1), i2 + 0)) * (0.0_wp)
      tt0 = tt0 + (x(mod_arr(1 + i1), i2 + 0)) * (0.75_wp)
      tt0 = tt0 + (x(mod_arr(2 + i1), i2 + 0)) * (-0.15_wp)
      tt0 = tt0 + (x(mod_arr(3 + i1), i2 + 0)) * (0.016666666666666666_wp)
      y(i1, i2 + 0) = tt0
    end do
    do i1 =  -(-3), n - (3) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(-3 + i1, i2 + 0)) * (-0.016666666666666666_wp)
      tt0 = tt0 + (x(-2 + i1, i2 + 0)) * (0.15_wp)
      tt0 = tt0 + (x(-1 + i1, i2 + 0)) * (-0.75_wp)
      tt0 = tt0 + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1, i2 + 0)) * (0.75_wp)
      tt0 = tt0 + (x(2 + i1, i2 + 0)) * (-0.15_wp)
      tt0 = tt0 + (x(3 + i1, i2 + 0)) * (0.016666666666666666_wp)
      y(i1, i2 + 0) = tt0
    end do
    do i1 = n - (3), n - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(mod_arr(-3 + i1 - (n)), i2 + 0)) * (-0.016666666666666666_wp)
      tt0 = tt0 + (x(mod_arr(-2 + i1 - (n)), i2 + 0)) * (0.15_wp)
      tt0 = tt0 + (x(mod_arr(-1 + i1 - (n)), i2 + 0)) * (-0.75_wp)
      tt0 = tt0 + (x(mod_arr(0 + i1 - (n)), i2 + 0)) * (0.0_wp)
      tt0 = tt0 + (x(mod_arr(1 + i1 - (n)), i2 + 0)) * (0.75_wp)
      tt0 = tt0 + (x(mod_arr(2 + i1 - (n)), i2 + 0)) * (-0.15_wp)
      tt0 = tt0 + (x(mod_arr(3 + i1 - (n)), i2 + 0)) * (0.016666666666666666_wp)
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_p_10_u1_1_true_false_true
SUBROUTINE d_poisson6_p_10_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_p_10_u1_1_true_false_true_cost
SUBROUTINE d_poisson6_p_10_a_u1_1_false_true_false(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:0) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 = 0,  -(-3) - (1), 1
      tt(0) = 0.0_wp
      do l = -3, 3, 1
        tt(0) = tt(0) + (x(l + i1 - (((i1 + l + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-3), n - (3) - (1), 1
      tt(0) = 0.0_wp
      do l = -3, 3, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (3), n - (1), 1
      tt(0) = 0.0_wp
      do l = -3, 3, 1
        tt(0) = tt(0) + (x(l + i1 - (((i1 + l + (n) * (2)) / (n) - (2)) * (n)), i2 + 0)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_p_10_a_u1_1_false_true_false
SUBROUTINE d_poisson6_p_10_a_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_p_10_a_u1_1_true_false_true_cost
SUBROUTINE d_poisson6_p_01_u2_0_true_false_false(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  integer(kind=4), dimension(-3 - (3):3 - (-3) - (1)) :: mod_arr
  do l = -3 - (3), 3 - (-3) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt0, tt1)
!$omp do 
  do i1 = 0, ndat - (2), 2
    do i2 = 0,  -(-3) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2))) * (poisson6_3_fil(l))
        tt1 = tt1 + (x(i1 + 1, mod_arr(l + i2))) * (poisson6_3_fil(l))
      end do
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
    end do
    do i2 =  -(-3), n - (3) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson6_3_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson6_3_fil(l))
      end do
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
    end do
    do i2 = n - (3), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2 - (n)))) * (poisson6_3_fil(l))
        tt1 = tt1 + (x(i1 + 1, mod_arr(l + i2 - (n)))) * (poisson6_3_fil(l))
      end do
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i2 = 0,  -(-3) - (1), 1
      tt0 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2))) * (poisson6_3_fil(l))
      end do
      y(i1 + 0, i2) = tt0
    end do
    do i2 =  -(-3), n - (3) - (1), 1
      tt0 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson6_3_fil(l))
      end do
      y(i1 + 0, i2) = tt0
    end do
    do i2 = n - (3), n - (1), 1
      tt0 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2 - (n)))) * (poisson6_3_fil(l))
      end do
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_p_01_u2_0_true_false_false
SUBROUTINE d_poisson6_p_01_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_p_01_u1_0_true_false_true_cost
SUBROUTINE d_poisson6_p_01_a_u4_0_false_true_true(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:3) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 = 0,  -(-3) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -3 + i2 - (((i2 + -3 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.016666666666666666_wp)
      tt(1) = tt(1) + (x(i1 + 1, -3 + i2 - (((i2 + -3 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.016666666666666666_wp)
      tt(2) = tt(2) + (x(i1 + 2, -3 + i2 - (((i2 + -3 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.016666666666666666_wp)
      tt(3) = tt(3) + (x(i1 + 3, -3 + i2 - (((i2 + -3 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.016666666666666666_wp)
      tt(0) = tt(0) + (x(i1 + 0, -2 + i2 - (((i2 + -2 + (n) * (2)) / (n) - (2)) * (n)))) * (0.15_wp)
      tt(1) = tt(1) + (x(i1 + 1, -2 + i2 - (((i2 + -2 + (n) * (2)) / (n) - (2)) * (n)))) * (0.15_wp)
      tt(2) = tt(2) + (x(i1 + 2, -2 + i2 - (((i2 + -2 + (n) * (2)) / (n) - (2)) * (n)))) * (0.15_wp)
      tt(3) = tt(3) + (x(i1 + 3, -2 + i2 - (((i2 + -2 + (n) * (2)) / (n) - (2)) * (n)))) * (0.15_wp)
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.75_wp)
      tt(1) = tt(1) + (x(i1 + 1, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.75_wp)
      tt(2) = tt(2) + (x(i1 + 2, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.75_wp)
      tt(3) = tt(3) + (x(i1 + 3, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.75_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(1) = tt(1) + (x(i1 + 1, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(2) = tt(2) + (x(i1 + 2, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(3) = tt(3) + (x(i1 + 3, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.75_wp)
      tt(1) = tt(1) + (x(i1 + 1, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.75_wp)
      tt(2) = tt(2) + (x(i1 + 2, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.75_wp)
      tt(3) = tt(3) + (x(i1 + 3, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.75_wp)
      tt(0) = tt(0) + (x(i1 + 0, 2 + i2 - (((i2 + 2 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.15_wp)
      tt(1) = tt(1) + (x(i1 + 1, 2 + i2 - (((i2 + 2 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.15_wp)
      tt(2) = tt(2) + (x(i1 + 2, 2 + i2 - (((i2 + 2 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.15_wp)
      tt(3) = tt(3) + (x(i1 + 3, 2 + i2 - (((i2 + 2 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.15_wp)
      tt(0) = tt(0) + (x(i1 + 0, 3 + i2 - (((i2 + 3 + (n) * (2)) / (n) - (2)) * (n)))) * (0.016666666666666666_wp)
      tt(1) = tt(1) + (x(i1 + 1, 3 + i2 - (((i2 + 3 + (n) * (2)) / (n) - (2)) * (n)))) * (0.016666666666666666_wp)
      tt(2) = tt(2) + (x(i1 + 2, 3 + i2 - (((i2 + 3 + (n) * (2)) / (n) - (2)) * (n)))) * (0.016666666666666666_wp)
      tt(3) = tt(3) + (x(i1 + 3, 3 + i2 - (((i2 + 3 + (n) * (2)) / (n) - (2)) * (n)))) * (0.016666666666666666_wp)
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
    end do
    do i2 =  -(-3), n - (3) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -3 + i2)) * (-0.016666666666666666_wp)
      tt(1) = tt(1) + (x(i1 + 1, -3 + i2)) * (-0.016666666666666666_wp)
      tt(2) = tt(2) + (x(i1 + 2, -3 + i2)) * (-0.016666666666666666_wp)
      tt(3) = tt(3) + (x(i1 + 3, -3 + i2)) * (-0.016666666666666666_wp)
      tt(0) = tt(0) + (x(i1 + 0, -2 + i2)) * (0.15_wp)
      tt(1) = tt(1) + (x(i1 + 1, -2 + i2)) * (0.15_wp)
      tt(2) = tt(2) + (x(i1 + 2, -2 + i2)) * (0.15_wp)
      tt(3) = tt(3) + (x(i1 + 3, -2 + i2)) * (0.15_wp)
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2)) * (-0.75_wp)
      tt(1) = tt(1) + (x(i1 + 1, -1 + i2)) * (-0.75_wp)
      tt(2) = tt(2) + (x(i1 + 2, -1 + i2)) * (-0.75_wp)
      tt(3) = tt(3) + (x(i1 + 3, -1 + i2)) * (-0.75_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt(1) = tt(1) + (x(i1 + 1, 0 + i2)) * (0.0_wp)
      tt(2) = tt(2) + (x(i1 + 2, 0 + i2)) * (0.0_wp)
      tt(3) = tt(3) + (x(i1 + 3, 0 + i2)) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2)) * (0.75_wp)
      tt(1) = tt(1) + (x(i1 + 1, 1 + i2)) * (0.75_wp)
      tt(2) = tt(2) + (x(i1 + 2, 1 + i2)) * (0.75_wp)
      tt(3) = tt(3) + (x(i1 + 3, 1 + i2)) * (0.75_wp)
      tt(0) = tt(0) + (x(i1 + 0, 2 + i2)) * (-0.15_wp)
      tt(1) = tt(1) + (x(i1 + 1, 2 + i2)) * (-0.15_wp)
      tt(2) = tt(2) + (x(i1 + 2, 2 + i2)) * (-0.15_wp)
      tt(3) = tt(3) + (x(i1 + 3, 2 + i2)) * (-0.15_wp)
      tt(0) = tt(0) + (x(i1 + 0, 3 + i2)) * (0.016666666666666666_wp)
      tt(1) = tt(1) + (x(i1 + 1, 3 + i2)) * (0.016666666666666666_wp)
      tt(2) = tt(2) + (x(i1 + 2, 3 + i2)) * (0.016666666666666666_wp)
      tt(3) = tt(3) + (x(i1 + 3, 3 + i2)) * (0.016666666666666666_wp)
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
    end do
    do i2 = n - (3), n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -3 + i2 - (((i2 + -3 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.016666666666666666_wp)
      tt(1) = tt(1) + (x(i1 + 1, -3 + i2 - (((i2 + -3 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.016666666666666666_wp)
      tt(2) = tt(2) + (x(i1 + 2, -3 + i2 - (((i2 + -3 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.016666666666666666_wp)
      tt(3) = tt(3) + (x(i1 + 3, -3 + i2 - (((i2 + -3 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.016666666666666666_wp)
      tt(0) = tt(0) + (x(i1 + 0, -2 + i2 - (((i2 + -2 + (n) * (2)) / (n) - (2)) * (n)))) * (0.15_wp)
      tt(1) = tt(1) + (x(i1 + 1, -2 + i2 - (((i2 + -2 + (n) * (2)) / (n) - (2)) * (n)))) * (0.15_wp)
      tt(2) = tt(2) + (x(i1 + 2, -2 + i2 - (((i2 + -2 + (n) * (2)) / (n) - (2)) * (n)))) * (0.15_wp)
      tt(3) = tt(3) + (x(i1 + 3, -2 + i2 - (((i2 + -2 + (n) * (2)) / (n) - (2)) * (n)))) * (0.15_wp)
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.75_wp)
      tt(1) = tt(1) + (x(i1 + 1, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.75_wp)
      tt(2) = tt(2) + (x(i1 + 2, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.75_wp)
      tt(3) = tt(3) + (x(i1 + 3, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.75_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(1) = tt(1) + (x(i1 + 1, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(2) = tt(2) + (x(i1 + 2, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(3) = tt(3) + (x(i1 + 3, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.75_wp)
      tt(1) = tt(1) + (x(i1 + 1, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.75_wp)
      tt(2) = tt(2) + (x(i1 + 2, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.75_wp)
      tt(3) = tt(3) + (x(i1 + 3, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.75_wp)
      tt(0) = tt(0) + (x(i1 + 0, 2 + i2 - (((i2 + 2 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.15_wp)
      tt(1) = tt(1) + (x(i1 + 1, 2 + i2 - (((i2 + 2 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.15_wp)
      tt(2) = tt(2) + (x(i1 + 2, 2 + i2 - (((i2 + 2 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.15_wp)
      tt(3) = tt(3) + (x(i1 + 3, 2 + i2 - (((i2 + 2 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.15_wp)
      tt(0) = tt(0) + (x(i1 + 0, 3 + i2 - (((i2 + 3 + (n) * (2)) / (n) - (2)) * (n)))) * (0.016666666666666666_wp)
      tt(1) = tt(1) + (x(i1 + 1, 3 + i2 - (((i2 + 3 + (n) * (2)) / (n) - (2)) * (n)))) * (0.016666666666666666_wp)
      tt(2) = tt(2) + (x(i1 + 2, 3 + i2 - (((i2 + 3 + (n) * (2)) / (n) - (2)) * (n)))) * (0.016666666666666666_wp)
      tt(3) = tt(3) + (x(i1 + 3, 3 + i2 - (((i2 + 3 + (n) * (2)) / (n) - (2)) * (n)))) * (0.016666666666666666_wp)
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 = 0,  -(-3) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -3 + i2 - (((i2 + -3 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.016666666666666666_wp)
      tt(0) = tt(0) + (x(i1 + 0, -2 + i2 - (((i2 + -2 + (n) * (2)) / (n) - (2)) * (n)))) * (0.15_wp)
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.75_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.75_wp)
      tt(0) = tt(0) + (x(i1 + 0, 2 + i2 - (((i2 + 2 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.15_wp)
      tt(0) = tt(0) + (x(i1 + 0, 3 + i2 - (((i2 + 3 + (n) * (2)) / (n) - (2)) * (n)))) * (0.016666666666666666_wp)
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 =  -(-3), n - (3) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -3 + i2)) * (-0.016666666666666666_wp)
      tt(0) = tt(0) + (x(i1 + 0, -2 + i2)) * (0.15_wp)
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2)) * (-0.75_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2)) * (0.75_wp)
      tt(0) = tt(0) + (x(i1 + 0, 2 + i2)) * (-0.15_wp)
      tt(0) = tt(0) + (x(i1 + 0, 3 + i2)) * (0.016666666666666666_wp)
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 = n - (3), n - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -3 + i2 - (((i2 + -3 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.016666666666666666_wp)
      tt(0) = tt(0) + (x(i1 + 0, -2 + i2 - (((i2 + -2 + (n) * (2)) / (n) - (2)) * (n)))) * (0.15_wp)
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2 - (((i2 + -1 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.75_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2 - (((i2 + 0 + (n) * (2)) / (n) - (2)) * (n)))) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2 - (((i2 + 1 + (n) * (2)) / (n) - (2)) * (n)))) * (0.75_wp)
      tt(0) = tt(0) + (x(i1 + 0, 2 + i2 - (((i2 + 2 + (n) * (2)) / (n) - (2)) * (n)))) * (-0.15_wp)
      tt(0) = tt(0) + (x(i1 + 0, 3 + i2 - (((i2 + 3 + (n) * (2)) / (n) - (2)) * (n)))) * (0.016666666666666666_wp)
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_p_01_a_u4_0_false_true_true
SUBROUTINE d_poisson6_p_01_a_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_p_01_a_u1_0_true_false_true_cost
SUBROUTINE d_poisson6_p_201_u2_0_true_true_true(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
  integer(kind=4), dimension(-3 - (3):3 - (-3) - (1)) :: mod_arr
  do l = -3 - (3), 3 - (-3) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (2), 2
      do i2 = 0,  -(-3) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(-3 + i2), i3)) * (-0.016666666666666666_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(-3 + i2), i3)) * (-0.016666666666666666_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(-2 + i2), i3)) * (0.15_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(-2 + i2), i3)) * (0.15_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(-1 + i2), i3)) * (-0.75_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(-1 + i2), i3)) * (-0.75_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(1 + i2), i3)) * (0.75_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(1 + i2), i3)) * (0.75_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(2 + i2), i3)) * (-0.15_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(2 + i2), i3)) * (-0.15_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(3 + i2), i3)) * (0.016666666666666666_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(3 + i2), i3)) * (0.016666666666666666_wp)
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
      end do
      do i2 =  -(-3), n - (3) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -3 + i2, i3)) * (-0.016666666666666666_wp)
        tt(1) = tt(1) + (x(i1 + 1, -3 + i2, i3)) * (-0.016666666666666666_wp)
        tt(0) = tt(0) + (x(i1 + 0, -2 + i2, i3)) * (0.15_wp)
        tt(1) = tt(1) + (x(i1 + 1, -2 + i2, i3)) * (0.15_wp)
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.75_wp)
        tt(1) = tt(1) + (x(i1 + 1, -1 + i2, i3)) * (-0.75_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(1) = tt(1) + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.75_wp)
        tt(1) = tt(1) + (x(i1 + 1, 1 + i2, i3)) * (0.75_wp)
        tt(0) = tt(0) + (x(i1 + 0, 2 + i2, i3)) * (-0.15_wp)
        tt(1) = tt(1) + (x(i1 + 1, 2 + i2, i3)) * (-0.15_wp)
        tt(0) = tt(0) + (x(i1 + 0, 3 + i2, i3)) * (0.016666666666666666_wp)
        tt(1) = tt(1) + (x(i1 + 1, 3 + i2, i3)) * (0.016666666666666666_wp)
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
      end do
      do i2 = n - (3), n - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(-3 + i2 - (n)), i3)) * (-0.016666666666666666_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(-3 + i2 - (n)), i3)) * (-0.016666666666666666_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(-2 + i2 - (n)), i3)) * (0.15_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(-2 + i2 - (n)), i3)) * (0.15_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(-1 + i2 - (n)), i3)) * (-0.75_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(-1 + i2 - (n)), i3)) * (-0.75_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(1 + i2 - (n)), i3)) * (0.75_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(1 + i2 - (n)), i3)) * (0.75_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(2 + i2 - (n)), i3)) * (-0.15_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(2 + i2 - (n)), i3)) * (-0.15_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(3 + i2 - (n)), i3)) * (0.016666666666666666_wp)
        tt(1) = tt(1) + (x(i1 + 1, mod_arr(3 + i2 - (n)), i3)) * (0.016666666666666666_wp)
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (2)) * (2), ndat1 - (1), 1
      do i2 = 0,  -(-3) - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(-3 + i2), i3)) * (-0.016666666666666666_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(-2 + i2), i3)) * (0.15_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(-1 + i2), i3)) * (-0.75_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(1 + i2), i3)) * (0.75_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(2 + i2), i3)) * (-0.15_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(3 + i2), i3)) * (0.016666666666666666_wp)
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 =  -(-3), n - (3) - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -3 + i2, i3)) * (-0.016666666666666666_wp)
        tt(0) = tt(0) + (x(i1 + 0, -2 + i2, i3)) * (0.15_wp)
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.75_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.75_wp)
        tt(0) = tt(0) + (x(i1 + 0, 2 + i2, i3)) * (-0.15_wp)
        tt(0) = tt(0) + (x(i1 + 0, 3 + i2, i3)) * (0.016666666666666666_wp)
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 = n - (3), n - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(-3 + i2 - (n)), i3)) * (-0.016666666666666666_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(-2 + i2 - (n)), i3)) * (0.15_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(-1 + i2 - (n)), i3)) * (-0.75_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(1 + i2 - (n)), i3)) * (0.75_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(2 + i2 - (n)), i3)) * (-0.15_wp)
        tt(0) = tt(0) + (x(i1 + 0, mod_arr(3 + i2 - (n)), i3)) * (0.016666666666666666_wp)
        y(i1 + 0, i2, i3) = tt(0)
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_p_201_u2_0_true_true_true
SUBROUTINE d_poisson6_p_201_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_p_201_u1_2_true_false_true_cost
SUBROUTINE d_poisson6_p_201_a_u4_0_false_false_false(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1, tt2, tt3)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 = 0,  -(-3) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        do l = -3, 3, 1
          tt0 = tt0 + (x(i1 + 0, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)), i3)) * (poisson6_3_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)), i3)) * (poisson6_3_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)), i3)) * (poisson6_3_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)), i3)) * (poisson6_3_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
      do i2 =  -(-3), n - (3) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        do l = -3, 3, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson6_3_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson6_3_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2, i3)) * (poisson6_3_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2, i3)) * (poisson6_3_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
      do i2 = n - (3), n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        do l = -3, 3, 1
          tt0 = tt0 + (x(i1 + 0, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)), i3)) * (poisson6_3_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)), i3)) * (poisson6_3_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)), i3)) * (poisson6_3_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)), i3)) * (poisson6_3_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 = 0,  -(-3) - (1), 1
        tt0 = 0.0_wp
        do l = -3, 3, 1
          tt0 = tt0 + (x(i1 + 0, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)), i3)) * (poisson6_3_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-3), n - (3) - (1), 1
        tt0 = 0.0_wp
        do l = -3, 3, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson6_3_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (3), n - (1), 1
        tt0 = 0.0_wp
        do l = -3, 3, 1
          tt0 = tt0 + (x(i1 + 0, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)), i3)) * (poisson6_3_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_p_201_a_u4_0_false_false_false
SUBROUTINE d_poisson6_p_201_a_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_p_201_a_u1_2_true_false_true_cost
SUBROUTINE d_poisson6_fg_10_u1_1_false_false_false(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension( -(3):n - (-3) - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
!$omp parallel  default(shared) private(i1, i2, tt0)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 =  -(3),  -(-3) - (1), 1
      tt0 = 0.0_wp
      do l = max( -(i1), -3), 3, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson6_3_fil(l))
      end do
      y(i1, i2 + 0) = tt0
    end do
    do i1 =  -(-3), n - (3) - (1), 1
      tt0 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson6_3_fil(l))
      end do
      y(i1, i2 + 0) = tt0
    end do
    do i1 = n - (3), n - (-3) - (1), 1
      tt0 = 0.0_wp
      do l = -3, min(3, n - (1) - (i1)), 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson6_3_fil(l))
      end do
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_fg_10_u1_1_false_false_false
SUBROUTINE d_poisson6_fg_10_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_fg_10_u1_1_false_false_true_cost
SUBROUTINE d_poisson6_fg_10_a_u2_1_false_true_false(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension( -(3):n - (-3) - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 =  -(3),  -(-3) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = max( -(i1), -3), 3, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson6_3_fil(l))
        tt(1) = tt(1) + (x(l + i1, i2 + 1)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
    do i1 =  -(-3), n - (3) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -3, 3, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson6_3_fil(l))
        tt(1) = tt(1) + (x(l + i1, i2 + 1)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
    do i1 = n - (3), n - (-3) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -3, min(3, n - (1) - (i1)), 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson6_3_fil(l))
        tt(1) = tt(1) + (x(l + i1, i2 + 1)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 =  -(3),  -(-3) - (1), 1
      tt(0) = 0.0_wp
      do l = max( -(i1), -3), 3, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-3), n - (3) - (1), 1
      tt(0) = 0.0_wp
      do l = -3, 3, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (3), n - (-3) - (1), 1
      tt(0) = 0.0_wp
      do l = -3, min(3, n - (1) - (i1)), 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_fg_10_a_u2_1_false_true_false
SUBROUTINE d_poisson6_fg_10_a_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_fg_10_a_u1_1_false_false_true_cost
SUBROUTINE d_poisson6_fg_01_u4_0_false_true_true(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1),  -(3):n - (-3) - (1)) :: y
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:3) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 =  -(3),  -(-3) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      do l = max( -(i2), -3), 3, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson6_3_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson6_3_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson6_3_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson6_3_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
      y(i1 + 2, i2) = tt(2)
      y(i1 + 3, i2) = tt(3)
    end do
    do i2 =  -(-3), n - (3) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -3 + i2)) * (-0.016666666666666666_wp)
      tt(1) = tt(1) + (x(i1 + 1, -3 + i2)) * (-0.016666666666666666_wp)
      tt(2) = tt(2) + (x(i1 + 2, -3 + i2)) * (-0.016666666666666666_wp)
      tt(3) = tt(3) + (x(i1 + 3, -3 + i2)) * (-0.016666666666666666_wp)
      tt(0) = tt(0) + (x(i1 + 0, -2 + i2)) * (0.15_wp)
      tt(1) = tt(1) + (x(i1 + 1, -2 + i2)) * (0.15_wp)
      tt(2) = tt(2) + (x(i1 + 2, -2 + i2)) * (0.15_wp)
      tt(3) = tt(3) + (x(i1 + 3, -2 + i2)) * (0.15_wp)
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2)) * (-0.75_wp)
      tt(1) = tt(1) + (x(i1 + 1, -1 + i2)) * (-0.75_wp)
      tt(2) = tt(2) + (x(i1 + 2, -1 + i2)) * (-0.75_wp)
      tt(3) = tt(3) + (x(i1 + 3, -1 + i2)) * (-0.75_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt(1) = tt(1) + (x(i1 + 1, 0 + i2)) * (0.0_wp)
      tt(2) = tt(2) + (x(i1 + 2, 0 + i2)) * (0.0_wp)
      tt(3) = tt(3) + (x(i1 + 3, 0 + i2)) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2)) * (0.75_wp)
      tt(1) = tt(1) + (x(i1 + 1, 1 + i2)) * (0.75_wp)
      tt(2) = tt(2) + (x(i1 + 2, 1 + i2)) * (0.75_wp)
      tt(3) = tt(3) + (x(i1 + 3, 1 + i2)) * (0.75_wp)
      tt(0) = tt(0) + (x(i1 + 0, 2 + i2)) * (-0.15_wp)
      tt(1) = tt(1) + (x(i1 + 1, 2 + i2)) * (-0.15_wp)
      tt(2) = tt(2) + (x(i1 + 2, 2 + i2)) * (-0.15_wp)
      tt(3) = tt(3) + (x(i1 + 3, 2 + i2)) * (-0.15_wp)
      tt(0) = tt(0) + (x(i1 + 0, 3 + i2)) * (0.016666666666666666_wp)
      tt(1) = tt(1) + (x(i1 + 1, 3 + i2)) * (0.016666666666666666_wp)
      tt(2) = tt(2) + (x(i1 + 2, 3 + i2)) * (0.016666666666666666_wp)
      tt(3) = tt(3) + (x(i1 + 3, 3 + i2)) * (0.016666666666666666_wp)
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
      y(i1 + 2, i2) = tt(2)
      y(i1 + 3, i2) = tt(3)
    end do
    do i2 = n - (3), n - (-3) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      do l = -3, min(3, n - (1) - (i2)), 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson6_3_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson6_3_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson6_3_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson6_3_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
      y(i1 + 2, i2) = tt(2)
      y(i1 + 3, i2) = tt(3)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 =  -(3),  -(-3) - (1), 1
      tt(0) = 0.0_wp
      do l = max( -(i2), -3), 3, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson6_3_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 =  -(-3), n - (3) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -3 + i2)) * (-0.016666666666666666_wp)
      tt(0) = tt(0) + (x(i1 + 0, -2 + i2)) * (0.15_wp)
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2)) * (-0.75_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2)) * (0.75_wp)
      tt(0) = tt(0) + (x(i1 + 0, 2 + i2)) * (-0.15_wp)
      tt(0) = tt(0) + (x(i1 + 0, 3 + i2)) * (0.016666666666666666_wp)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 = n - (3), n - (-3) - (1), 1
      tt(0) = 0.0_wp
      do l = -3, min(3, n - (1) - (i2)), 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson6_3_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_fg_01_u4_0_false_true_true
SUBROUTINE d_poisson6_fg_01_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_fg_01_u1_0_false_false_true_cost
SUBROUTINE d_poisson6_fg_01_a_u4_0_false_true_false(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1),  -(3):n - (-3) - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:3) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 =  -(3),  -(-3) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      do l = max( -(i2), -3), 3, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson6_3_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson6_3_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson6_3_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
    end do
    do i2 =  -(-3), n - (3) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      do l = -3, 3, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson6_3_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson6_3_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson6_3_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
    end do
    do i2 = n - (3), n - (-3) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      do l = -3, min(3, n - (1) - (i2)), 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson6_3_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson6_3_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson6_3_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 =  -(3),  -(-3) - (1), 1
      tt(0) = 0.0_wp
      do l = max( -(i2), -3), 3, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 =  -(-3), n - (3) - (1), 1
      tt(0) = 0.0_wp
      do l = -3, 3, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 = n - (3), n - (-3) - (1), 1
      tt(0) = 0.0_wp
      do l = -3, min(3, n - (1) - (i2)), 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_fg_01_a_u4_0_false_true_false
SUBROUTINE d_poisson6_fg_01_a_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_fg_01_a_u1_0_false_false_true_cost
SUBROUTINE d_poisson6_fg_201_u4_0_false_false_false(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1),  -(3):n - (-3) - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1, tt2, tt3)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 =  -(3),  -(-3) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        do l = max( -(i2), -3), 3, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson6_3_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson6_3_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2, i3)) * (poisson6_3_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2, i3)) * (poisson6_3_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
        y(i1 + 2, i2, i3) = tt2
        y(i1 + 3, i2, i3) = tt3
      end do
      do i2 =  -(-3), n - (3) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        do l = -3, 3, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson6_3_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson6_3_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2, i3)) * (poisson6_3_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2, i3)) * (poisson6_3_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
        y(i1 + 2, i2, i3) = tt2
        y(i1 + 3, i2, i3) = tt3
      end do
      do i2 = n - (3), n - (-3) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        do l = -3, min(3, n - (1) - (i2)), 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson6_3_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson6_3_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2, i3)) * (poisson6_3_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2, i3)) * (poisson6_3_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
        y(i1 + 2, i2, i3) = tt2
        y(i1 + 3, i2, i3) = tt3
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 =  -(3),  -(-3) - (1), 1
        tt0 = 0.0_wp
        do l = max( -(i2), -3), 3, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson6_3_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-3), n - (3) - (1), 1
        tt0 = 0.0_wp
        do l = -3, 3, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson6_3_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (3), n - (-3) - (1), 1
        tt0 = 0.0_wp
        do l = -3, min(3, n - (1) - (i2)), 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson6_3_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_fg_201_u4_0_false_false_false
SUBROUTINE d_poisson6_fg_201_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_fg_201_u1_2_false_false_true_cost
SUBROUTINE d_poisson6_fg_201_a_u2_0_false_false_true(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1),  -(3):n - (-3) - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (2), 2
      do i2 =  -(3),  -(-3) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        do l = max( -(i2), -3), 3, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson6_3_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson6_3_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
      end do
      do i2 =  -(-3), n - (3) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -3 + i2, i3)) * (-0.016666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, -3 + i2, i3)) * (-0.016666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, -2 + i2, i3)) * (0.15_wp)
        tt1 = tt1 + (x(i1 + 1, -2 + i2, i3)) * (0.15_wp)
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.75_wp)
        tt1 = tt1 + (x(i1 + 1, -1 + i2, i3)) * (-0.75_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt1 = tt1 + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.75_wp)
        tt1 = tt1 + (x(i1 + 1, 1 + i2, i3)) * (0.75_wp)
        tt0 = tt0 + (x(i1 + 0, 2 + i2, i3)) * (-0.15_wp)
        tt1 = tt1 + (x(i1 + 1, 2 + i2, i3)) * (-0.15_wp)
        tt0 = tt0 + (x(i1 + 0, 3 + i2, i3)) * (0.016666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, 3 + i2, i3)) * (0.016666666666666666_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
      end do
      do i2 = n - (3), n - (-3) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        do l = -3, min(3, n - (1) - (i2)), 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson6_3_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson6_3_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (2)) * (2), ndat1 - (1), 1
      do i2 =  -(3),  -(-3) - (1), 1
        tt0 = 0.0_wp
        do l = max( -(i2), -3), 3, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson6_3_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-3), n - (3) - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -3 + i2, i3)) * (-0.016666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, -2 + i2, i3)) * (0.15_wp)
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.75_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.75_wp)
        tt0 = tt0 + (x(i1 + 0, 2 + i2, i3)) * (-0.15_wp)
        tt0 = tt0 + (x(i1 + 0, 3 + i2, i3)) * (0.016666666666666666_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (3), n - (-3) - (1), 1
        tt0 = 0.0_wp
        do l = -3, min(3, n - (1) - (i2)), 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson6_3_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_fg_201_a_u2_0_false_false_true
SUBROUTINE d_poisson6_fg_201_a_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_fg_201_a_u1_2_false_false_true_cost
SUBROUTINE d_poisson6_fs_10_u1_1_false_true_false(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(-3:n + 3 - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:0) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 = 0, n - (1), 1
      tt(0) = 0.0_wp
      do l = -3, 3, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson6_3_fil(l))
      end do
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_fs_10_u1_1_false_true_false
SUBROUTINE d_poisson6_fs_10_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_fs_10_u1_1_false_false_true_cost
SUBROUTINE d_poisson6_fs_10_a_u2_1_false_false_true(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(-3:n + 3 - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
!$omp parallel  default(shared) private(i1, i2, tt0, tt1)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 = 0, n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt0 = tt0 + (x(-3 + i1, i2 + 0)) * (-0.016666666666666666_wp)
      tt1 = tt1 + (x(-3 + i1, i2 + 1)) * (-0.016666666666666666_wp)
      tt0 = tt0 + (x(-2 + i1, i2 + 0)) * (0.15_wp)
      tt1 = tt1 + (x(-2 + i1, i2 + 1)) * (0.15_wp)
      tt0 = tt0 + (x(-1 + i1, i2 + 0)) * (-0.75_wp)
      tt1 = tt1 + (x(-1 + i1, i2 + 1)) * (-0.75_wp)
      tt0 = tt0 + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt1 = tt1 + (x(0 + i1, i2 + 1)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1, i2 + 0)) * (0.75_wp)
      tt1 = tt1 + (x(1 + i1, i2 + 1)) * (0.75_wp)
      tt0 = tt0 + (x(2 + i1, i2 + 0)) * (-0.15_wp)
      tt1 = tt1 + (x(2 + i1, i2 + 1)) * (-0.15_wp)
      tt0 = tt0 + (x(3 + i1, i2 + 0)) * (0.016666666666666666_wp)
      tt1 = tt1 + (x(3 + i1, i2 + 1)) * (0.016666666666666666_wp)
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
      tt1 = (tt1) * (a)
      y(i1, i2 + 1) = tt1
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 = 0, n - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(-3 + i1, i2 + 0)) * (-0.016666666666666666_wp)
      tt0 = tt0 + (x(-2 + i1, i2 + 0)) * (0.15_wp)
      tt0 = tt0 + (x(-1 + i1, i2 + 0)) * (-0.75_wp)
      tt0 = tt0 + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1, i2 + 0)) * (0.75_wp)
      tt0 = tt0 + (x(2 + i1, i2 + 0)) * (-0.15_wp)
      tt0 = tt0 + (x(3 + i1, i2 + 0)) * (0.016666666666666666_wp)
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_fs_10_a_u2_1_false_false_true
SUBROUTINE d_poisson6_fs_10_a_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_fs_10_a_u1_1_false_false_true_cost
SUBROUTINE d_poisson6_fs_01_u3_0_false_false_true(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), -3:n + 3 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2)
!$omp do 
  do i1 = 0, ndat - (3), 3
    do i2 = 0, n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -3 + i2)) * (-0.016666666666666666_wp)
      tt1 = tt1 + (x(i1 + 1, -3 + i2)) * (-0.016666666666666666_wp)
      tt2 = tt2 + (x(i1 + 2, -3 + i2)) * (-0.016666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, -2 + i2)) * (0.15_wp)
      tt1 = tt1 + (x(i1 + 1, -2 + i2)) * (0.15_wp)
      tt2 = tt2 + (x(i1 + 2, -2 + i2)) * (0.15_wp)
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.75_wp)
      tt1 = tt1 + (x(i1 + 1, -1 + i2)) * (-0.75_wp)
      tt2 = tt2 + (x(i1 + 2, -1 + i2)) * (-0.75_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt1 = tt1 + (x(i1 + 1, 0 + i2)) * (0.0_wp)
      tt2 = tt2 + (x(i1 + 2, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.75_wp)
      tt1 = tt1 + (x(i1 + 1, 1 + i2)) * (0.75_wp)
      tt2 = tt2 + (x(i1 + 2, 1 + i2)) * (0.75_wp)
      tt0 = tt0 + (x(i1 + 0, 2 + i2)) * (-0.15_wp)
      tt1 = tt1 + (x(i1 + 1, 2 + i2)) * (-0.15_wp)
      tt2 = tt2 + (x(i1 + 2, 2 + i2)) * (-0.15_wp)
      tt0 = tt0 + (x(i1 + 0, 3 + i2)) * (0.016666666666666666_wp)
      tt1 = tt1 + (x(i1 + 1, 3 + i2)) * (0.016666666666666666_wp)
      tt2 = tt2 + (x(i1 + 2, 3 + i2)) * (0.016666666666666666_wp)
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (3)) * (3), ndat - (1), 1
    do i2 = 0, n - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -3 + i2)) * (-0.016666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, -2 + i2)) * (0.15_wp)
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.75_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.75_wp)
      tt0 = tt0 + (x(i1 + 0, 2 + i2)) * (-0.15_wp)
      tt0 = tt0 + (x(i1 + 0, 3 + i2)) * (0.016666666666666666_wp)
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_fs_01_u3_0_false_false_true
SUBROUTINE d_poisson6_fs_01_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_fs_01_u1_0_false_false_true_cost
SUBROUTINE d_poisson6_fs_01_a_u2_0_false_true_false(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), -3:n + 3 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (2), 2
    do i2 = 0, n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -3, 3, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson6_3_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i2 = 0, n - (1), 1
      tt(0) = 0.0_wp
      do l = -3, 3, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_fs_01_a_u2_0_false_true_false
SUBROUTINE d_poisson6_fs_01_a_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_fs_01_a_u1_0_false_false_true_cost
SUBROUTINE d_poisson6_fs_201_u2_0_false_false_false(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), -3:n + 3 - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (2), 2
      do i2 = 0, n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        do l = -3, 3, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson6_3_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson6_3_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (2)) * (2), ndat1 - (1), 1
      do i2 = 0, n - (1), 1
        tt0 = 0.0_wp
        do l = -3, 3, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson6_3_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_fs_201_u2_0_false_false_false
SUBROUTINE d_poisson6_fs_201_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_fs_201_u1_2_false_false_true_cost
SUBROUTINE d_poisson6_fs_201_a_u4_0_false_true_true(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), -3:n + 3 - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8), dimension(0:3) :: tt
!$omp parallel  default(shared) private(i1, i2, i3, tt)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 = 0, n - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -3 + i2, i3)) * (-0.016666666666666666_wp)
        tt(1) = tt(1) + (x(i1 + 1, -3 + i2, i3)) * (-0.016666666666666666_wp)
        tt(2) = tt(2) + (x(i1 + 2, -3 + i2, i3)) * (-0.016666666666666666_wp)
        tt(3) = tt(3) + (x(i1 + 3, -3 + i2, i3)) * (-0.016666666666666666_wp)
        tt(0) = tt(0) + (x(i1 + 0, -2 + i2, i3)) * (0.15_wp)
        tt(1) = tt(1) + (x(i1 + 1, -2 + i2, i3)) * (0.15_wp)
        tt(2) = tt(2) + (x(i1 + 2, -2 + i2, i3)) * (0.15_wp)
        tt(3) = tt(3) + (x(i1 + 3, -2 + i2, i3)) * (0.15_wp)
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.75_wp)
        tt(1) = tt(1) + (x(i1 + 1, -1 + i2, i3)) * (-0.75_wp)
        tt(2) = tt(2) + (x(i1 + 2, -1 + i2, i3)) * (-0.75_wp)
        tt(3) = tt(3) + (x(i1 + 3, -1 + i2, i3)) * (-0.75_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(1) = tt(1) + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt(2) = tt(2) + (x(i1 + 2, 0 + i2, i3)) * (0.0_wp)
        tt(3) = tt(3) + (x(i1 + 3, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.75_wp)
        tt(1) = tt(1) + (x(i1 + 1, 1 + i2, i3)) * (0.75_wp)
        tt(2) = tt(2) + (x(i1 + 2, 1 + i2, i3)) * (0.75_wp)
        tt(3) = tt(3) + (x(i1 + 3, 1 + i2, i3)) * (0.75_wp)
        tt(0) = tt(0) + (x(i1 + 0, 2 + i2, i3)) * (-0.15_wp)
        tt(1) = tt(1) + (x(i1 + 1, 2 + i2, i3)) * (-0.15_wp)
        tt(2) = tt(2) + (x(i1 + 2, 2 + i2, i3)) * (-0.15_wp)
        tt(3) = tt(3) + (x(i1 + 3, 2 + i2, i3)) * (-0.15_wp)
        tt(0) = tt(0) + (x(i1 + 0, 3 + i2, i3)) * (0.016666666666666666_wp)
        tt(1) = tt(1) + (x(i1 + 1, 3 + i2, i3)) * (0.016666666666666666_wp)
        tt(2) = tt(2) + (x(i1 + 2, 3 + i2, i3)) * (0.016666666666666666_wp)
        tt(3) = tt(3) + (x(i1 + 3, 3 + i2, i3)) * (0.016666666666666666_wp)
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
        tt(1) = (tt(1)) * (a)
        y(i1 + 1, i2, i3) = tt(1)
        tt(2) = (tt(2)) * (a)
        y(i1 + 2, i2, i3) = tt(2)
        tt(3) = (tt(3)) * (a)
        y(i1 + 3, i2, i3) = tt(3)
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 = 0, n - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -3 + i2, i3)) * (-0.016666666666666666_wp)
        tt(0) = tt(0) + (x(i1 + 0, -2 + i2, i3)) * (0.15_wp)
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.75_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.75_wp)
        tt(0) = tt(0) + (x(i1 + 0, 2 + i2, i3)) * (-0.15_wp)
        tt(0) = tt(0) + (x(i1 + 0, 3 + i2, i3)) * (0.016666666666666666_wp)
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_fs_201_a_u4_0_false_true_true
SUBROUTINE d_poisson6_fs_201_a_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_fs_201_a_u1_2_false_false_true_cost
SUBROUTINE d_poisson6_np_10_u2_1_false_false_false(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(0:48) :: poisson6_fil = (/ &
-2.45_wp, &
6.0_wp, &
-7.5_wp, &
6.666666666666667_wp, &
-3.75_wp, &
1.2_wp, &
-0.16666666666666666_wp, &
-0.16666666666666666_wp, &
-1.2833333333333334_wp, &
2.5_wp, &
-1.6666666666666667_wp, &
0.8333333333333334_wp, &
-0.25_wp, &
0.03333333333333333_wp, &
0.03333333333333333_wp, &
-0.4_wp, &
-0.5833333333333334_wp, &
1.3333333333333333_wp, &
-0.5_wp, &
0.13333333333333333_wp, &
-0.016666666666666666_wp, &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp, &
0.016666666666666666_wp, &
-0.13333333333333333_wp, &
0.5_wp, &
-1.3333333333333333_wp, &
0.5833333333333334_wp, &
0.4_wp, &
-0.03333333333333333_wp, &
-0.03333333333333333_wp, &
0.25_wp, &
-0.8333333333333334_wp, &
1.6666666666666667_wp, &
-2.5_wp, &
1.2833333333333334_wp, &
0.16666666666666666_wp, &
0.16666666666666666_wp, &
-1.2_wp, &
3.75_wp, &
-6.666666666666667_wp, &
7.5_wp, &
-6.0_wp, &
2.45_wp /)
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
!$omp parallel  default(shared) private(i1, i2, tt0, tt1)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 = 0,  -(-3) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(l + 3, i2 + 0)) * (poisson6_fil((i1 - (3)) * (7) + 21 + l + 3))
        tt1 = tt1 + (x(l + 3, i2 + 1)) * (poisson6_fil((i1 - (3)) * (7) + 21 + l + 3))
      end do
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
    end do
    do i1 =  -(-3), n - (3) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson6_3_fil(l))
        tt1 = tt1 + (x(l + i1, i2 + 1)) * (poisson6_3_fil(l))
      end do
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
    end do
    do i1 = n - (3), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(l + -3 + n - (1), i2 + 0)) * (poisson6_fil((i1 + 3 - (n) + 1) * (7) + 21 + l + 3))
        tt1 = tt1 + (x(l + -3 + n - (1), i2 + 1)) * (poisson6_fil((i1 + 3 - (n) + 1) * (7) + 21 + l + 3))
      end do
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 = 0,  -(-3) - (1), 1
      tt0 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(l + 3, i2 + 0)) * (poisson6_fil((i1 - (3)) * (7) + 21 + l + 3))
      end do
      y(i1, i2 + 0) = tt0
    end do
    do i1 =  -(-3), n - (3) - (1), 1
      tt0 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson6_3_fil(l))
      end do
      y(i1, i2 + 0) = tt0
    end do
    do i1 = n - (3), n - (1), 1
      tt0 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(l + -3 + n - (1), i2 + 0)) * (poisson6_fil((i1 + 3 - (n) + 1) * (7) + 21 + l + 3))
      end do
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_np_10_u2_1_false_false_false
SUBROUTINE d_poisson6_np_10_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_np_10_u1_1_true_false_true_cost
SUBROUTINE d_poisson6_np_10_a_u2_1_true_true_false(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(0:48) :: poisson6_fil = (/ &
-2.45_wp, &
6.0_wp, &
-7.5_wp, &
6.666666666666667_wp, &
-3.75_wp, &
1.2_wp, &
-0.16666666666666666_wp, &
-0.16666666666666666_wp, &
-1.2833333333333334_wp, &
2.5_wp, &
-1.6666666666666667_wp, &
0.8333333333333334_wp, &
-0.25_wp, &
0.03333333333333333_wp, &
0.03333333333333333_wp, &
-0.4_wp, &
-0.5833333333333334_wp, &
1.3333333333333333_wp, &
-0.5_wp, &
0.13333333333333333_wp, &
-0.016666666666666666_wp, &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp, &
0.016666666666666666_wp, &
-0.13333333333333333_wp, &
0.5_wp, &
-1.3333333333333333_wp, &
0.5833333333333334_wp, &
0.4_wp, &
-0.03333333333333333_wp, &
-0.03333333333333333_wp, &
0.25_wp, &
-0.8333333333333334_wp, &
1.6666666666666667_wp, &
-2.5_wp, &
1.2833333333333334_wp, &
0.16666666666666666_wp, &
0.16666666666666666_wp, &
-1.2_wp, &
3.75_wp, &
-6.666666666666667_wp, &
7.5_wp, &
-6.0_wp, &
2.45_wp /)
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
  integer(kind=4), dimension(-3 - (3):3 - (-3) - (1)) :: mod_arr
  do l = -3 - (3), 3 - (-3) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 = 0,  -(-3) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -3, 3, 1
        tt(0) = tt(0) + (x(l + 3, i2 + 0)) * (poisson6_fil((i1 - (3)) * (7) + 21 + l + 3))
        tt(1) = tt(1) + (x(l + 3, i2 + 1)) * (poisson6_fil((i1 - (3)) * (7) + 21 + l + 3))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
    do i1 =  -(-3), n - (3) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -3, 3, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson6_3_fil(l))
        tt(1) = tt(1) + (x(l + i1, i2 + 1)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
    do i1 = n - (3), n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -3, 3, 1
        tt(0) = tt(0) + (x(l + -3 + n - (1), i2 + 0)) * (poisson6_fil((i1 + 3 - (n) + 1) * (7) + 21 + l + 3))
        tt(1) = tt(1) + (x(l + -3 + n - (1), i2 + 1)) * (poisson6_fil((i1 + 3 - (n) + 1) * (7) + 21 + l + 3))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 = 0,  -(-3) - (1), 1
      tt(0) = 0.0_wp
      do l = -3, 3, 1
        tt(0) = tt(0) + (x(l + 3, i2 + 0)) * (poisson6_fil((i1 - (3)) * (7) + 21 + l + 3))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-3), n - (3) - (1), 1
      tt(0) = 0.0_wp
      do l = -3, 3, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson6_3_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (3), n - (1), 1
      tt(0) = 0.0_wp
      do l = -3, 3, 1
        tt(0) = tt(0) + (x(l + -3 + n - (1), i2 + 0)) * (poisson6_fil((i1 + 3 - (n) + 1) * (7) + 21 + l + 3))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_np_10_a_u2_1_true_true_false
SUBROUTINE d_poisson6_np_10_a_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_np_10_a_u1_1_true_false_true_cost
SUBROUTINE d_poisson6_np_01_u2_0_false_false_true(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), parameter, dimension(0:48) :: poisson6_fil = (/ &
-2.45_wp, &
6.0_wp, &
-7.5_wp, &
6.666666666666667_wp, &
-3.75_wp, &
1.2_wp, &
-0.16666666666666666_wp, &
-0.16666666666666666_wp, &
-1.2833333333333334_wp, &
2.5_wp, &
-1.6666666666666667_wp, &
0.8333333333333334_wp, &
-0.25_wp, &
0.03333333333333333_wp, &
0.03333333333333333_wp, &
-0.4_wp, &
-0.5833333333333334_wp, &
1.3333333333333333_wp, &
-0.5_wp, &
0.13333333333333333_wp, &
-0.016666666666666666_wp, &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp, &
0.016666666666666666_wp, &
-0.13333333333333333_wp, &
0.5_wp, &
-1.3333333333333333_wp, &
0.5833333333333334_wp, &
0.4_wp, &
-0.03333333333333333_wp, &
-0.03333333333333333_wp, &
0.25_wp, &
-0.8333333333333334_wp, &
1.6666666666666667_wp, &
-2.5_wp, &
1.2833333333333334_wp, &
0.16666666666666666_wp, &
0.16666666666666666_wp, &
-1.2_wp, &
3.75_wp, &
-6.666666666666667_wp, &
7.5_wp, &
-6.0_wp, &
2.45_wp /)
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
!$omp parallel  default(shared) private(i1, i2, tt0, tt1)
!$omp do 
  do i1 = 0, ndat - (2), 2
    do i2 = 0,  -(-3) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -3 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -3 + 3))
      tt1 = tt1 + (x(i1 + 1, -3 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -3 + 3))
      tt0 = tt0 + (x(i1 + 0, -2 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -2 + 3))
      tt1 = tt1 + (x(i1 + 1, -2 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -2 + 3))
      tt0 = tt0 + (x(i1 + 0, -1 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -1 + 3))
      tt1 = tt1 + (x(i1 + 1, -1 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -1 + 3))
      tt0 = tt0 + (x(i1 + 0, 0 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 0 + 3))
      tt1 = tt1 + (x(i1 + 1, 0 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 0 + 3))
      tt0 = tt0 + (x(i1 + 0, 1 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 1 + 3))
      tt1 = tt1 + (x(i1 + 1, 1 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 1 + 3))
      tt0 = tt0 + (x(i1 + 0, 2 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 2 + 3))
      tt1 = tt1 + (x(i1 + 1, 2 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 2 + 3))
      tt0 = tt0 + (x(i1 + 0, 3 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 3 + 3))
      tt1 = tt1 + (x(i1 + 1, 3 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 3 + 3))
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
    end do
    do i2 =  -(-3), n - (3) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -3 + i2)) * (-0.016666666666666666_wp)
      tt1 = tt1 + (x(i1 + 1, -3 + i2)) * (-0.016666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, -2 + i2)) * (0.15_wp)
      tt1 = tt1 + (x(i1 + 1, -2 + i2)) * (0.15_wp)
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.75_wp)
      tt1 = tt1 + (x(i1 + 1, -1 + i2)) * (-0.75_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt1 = tt1 + (x(i1 + 1, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.75_wp)
      tt1 = tt1 + (x(i1 + 1, 1 + i2)) * (0.75_wp)
      tt0 = tt0 + (x(i1 + 0, 2 + i2)) * (-0.15_wp)
      tt1 = tt1 + (x(i1 + 1, 2 + i2)) * (-0.15_wp)
      tt0 = tt0 + (x(i1 + 0, 3 + i2)) * (0.016666666666666666_wp)
      tt1 = tt1 + (x(i1 + 1, 3 + i2)) * (0.016666666666666666_wp)
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
    end do
    do i2 = n - (3), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -3 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -3 + 3))
      tt1 = tt1 + (x(i1 + 1, -3 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -3 + 3))
      tt0 = tt0 + (x(i1 + 0, -2 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -2 + 3))
      tt1 = tt1 + (x(i1 + 1, -2 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -2 + 3))
      tt0 = tt0 + (x(i1 + 0, -1 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -1 + 3))
      tt1 = tt1 + (x(i1 + 1, -1 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -1 + 3))
      tt0 = tt0 + (x(i1 + 0, 0 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 0 + 3))
      tt1 = tt1 + (x(i1 + 1, 0 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 0 + 3))
      tt0 = tt0 + (x(i1 + 0, 1 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 1 + 3))
      tt1 = tt1 + (x(i1 + 1, 1 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 1 + 3))
      tt0 = tt0 + (x(i1 + 0, 2 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 2 + 3))
      tt1 = tt1 + (x(i1 + 1, 2 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 2 + 3))
      tt0 = tt0 + (x(i1 + 0, 3 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 3 + 3))
      tt1 = tt1 + (x(i1 + 1, 3 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 3 + 3))
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i2 = 0,  -(-3) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -3 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -3 + 3))
      tt0 = tt0 + (x(i1 + 0, -2 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -2 + 3))
      tt0 = tt0 + (x(i1 + 0, -1 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -1 + 3))
      tt0 = tt0 + (x(i1 + 0, 0 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 0 + 3))
      tt0 = tt0 + (x(i1 + 0, 1 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 1 + 3))
      tt0 = tt0 + (x(i1 + 0, 2 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 2 + 3))
      tt0 = tt0 + (x(i1 + 0, 3 + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 3 + 3))
      y(i1 + 0, i2) = tt0
    end do
    do i2 =  -(-3), n - (3) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -3 + i2)) * (-0.016666666666666666_wp)
      tt0 = tt0 + (x(i1 + 0, -2 + i2)) * (0.15_wp)
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.75_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.75_wp)
      tt0 = tt0 + (x(i1 + 0, 2 + i2)) * (-0.15_wp)
      tt0 = tt0 + (x(i1 + 0, 3 + i2)) * (0.016666666666666666_wp)
      y(i1 + 0, i2) = tt0
    end do
    do i2 = n - (3), n - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -3 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -3 + 3))
      tt0 = tt0 + (x(i1 + 0, -2 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -2 + 3))
      tt0 = tt0 + (x(i1 + 0, -1 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -1 + 3))
      tt0 = tt0 + (x(i1 + 0, 0 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 0 + 3))
      tt0 = tt0 + (x(i1 + 0, 1 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 1 + 3))
      tt0 = tt0 + (x(i1 + 0, 2 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 2 + 3))
      tt0 = tt0 + (x(i1 + 0, 3 + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 3 + 3))
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_np_01_u2_0_false_false_true
SUBROUTINE d_poisson6_np_01_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_np_01_u1_0_true_false_true_cost
SUBROUTINE d_poisson6_np_01_a_u5_0_false_false_false(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(0:48) :: poisson6_fil = (/ &
-2.45_wp, &
6.0_wp, &
-7.5_wp, &
6.666666666666667_wp, &
-3.75_wp, &
1.2_wp, &
-0.16666666666666666_wp, &
-0.16666666666666666_wp, &
-1.2833333333333334_wp, &
2.5_wp, &
-1.6666666666666667_wp, &
0.8333333333333334_wp, &
-0.25_wp, &
0.03333333333333333_wp, &
0.03333333333333333_wp, &
-0.4_wp, &
-0.5833333333333334_wp, &
1.3333333333333333_wp, &
-0.5_wp, &
0.13333333333333333_wp, &
-0.016666666666666666_wp, &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp, &
0.016666666666666666_wp, &
-0.13333333333333333_wp, &
0.5_wp, &
-1.3333333333333333_wp, &
0.5833333333333334_wp, &
0.4_wp, &
-0.03333333333333333_wp, &
-0.03333333333333333_wp, &
0.25_wp, &
-0.8333333333333334_wp, &
1.6666666666666667_wp, &
-2.5_wp, &
1.2833333333333334_wp, &
0.16666666666666666_wp, &
0.16666666666666666_wp, &
-1.2_wp, &
3.75_wp, &
-6.666666666666667_wp, &
7.5_wp, &
-6.0_wp, &
2.45_wp /)
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2, tt3, tt4)
!$omp do 
  do i1 = 0, ndat - (5), 5
    do i2 = 0,  -(-3) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt4 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(i1 + 0, l + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + l + 3))
        tt1 = tt1 + (x(i1 + 1, l + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + l + 3))
        tt2 = tt2 + (x(i1 + 2, l + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + l + 3))
        tt3 = tt3 + (x(i1 + 3, l + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + l + 3))
        tt4 = tt4 + (x(i1 + 4, l + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + l + 3))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
      tt4 = (tt4) * (a)
      y(i1 + 4, i2) = tt4
    end do
    do i2 =  -(-3), n - (3) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt4 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson6_3_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson6_3_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2)) * (poisson6_3_fil(l))
        tt3 = tt3 + (x(i1 + 3, l + i2)) * (poisson6_3_fil(l))
        tt4 = tt4 + (x(i1 + 4, l + i2)) * (poisson6_3_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
      tt4 = (tt4) * (a)
      y(i1 + 4, i2) = tt4
    end do
    do i2 = n - (3), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt4 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(i1 + 0, l + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + l + 3))
        tt1 = tt1 + (x(i1 + 1, l + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + l + 3))
        tt2 = tt2 + (x(i1 + 2, l + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + l + 3))
        tt3 = tt3 + (x(i1 + 3, l + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + l + 3))
        tt4 = tt4 + (x(i1 + 4, l + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + l + 3))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
      tt4 = (tt4) * (a)
      y(i1 + 4, i2) = tt4
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (5)) * (5), ndat - (1), 1
    do i2 = 0,  -(-3) - (1), 1
      tt0 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(i1 + 0, l + 3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + l + 3))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
    do i2 =  -(-3), n - (3) - (1), 1
      tt0 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson6_3_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
    do i2 = n - (3), n - (1), 1
      tt0 = 0.0_wp
      do l = -3, 3, 1
        tt0 = tt0 + (x(i1 + 0, l + -3 + n - (1))) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + l + 3))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_np_01_a_u5_0_false_false_false
SUBROUTINE d_poisson6_np_01_a_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_np_01_a_u1_0_true_false_true_cost
SUBROUTINE d_poisson6_np_201_u4_0_true_true_false(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(0:48) :: poisson6_fil = (/ &
-2.45_wp, &
6.0_wp, &
-7.5_wp, &
6.666666666666667_wp, &
-3.75_wp, &
1.2_wp, &
-0.16666666666666666_wp, &
-0.16666666666666666_wp, &
-1.2833333333333334_wp, &
2.5_wp, &
-1.6666666666666667_wp, &
0.8333333333333334_wp, &
-0.25_wp, &
0.03333333333333333_wp, &
0.03333333333333333_wp, &
-0.4_wp, &
-0.5833333333333334_wp, &
1.3333333333333333_wp, &
-0.5_wp, &
0.13333333333333333_wp, &
-0.016666666666666666_wp, &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp, &
0.016666666666666666_wp, &
-0.13333333333333333_wp, &
0.5_wp, &
-1.3333333333333333_wp, &
0.5833333333333334_wp, &
0.4_wp, &
-0.03333333333333333_wp, &
-0.03333333333333333_wp, &
0.25_wp, &
-0.8333333333333334_wp, &
1.6666666666666667_wp, &
-2.5_wp, &
1.2833333333333334_wp, &
0.16666666666666666_wp, &
0.16666666666666666_wp, &
-1.2_wp, &
3.75_wp, &
-6.666666666666667_wp, &
7.5_wp, &
-6.0_wp, &
2.45_wp /)
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8), dimension(0:3) :: tt
  integer(kind=4), dimension(-3 - (3):3 - (-3) - (1)) :: mod_arr
  do l = -3 - (3), 3 - (-3) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 = 0,  -(-3) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        do l = -3, 3, 1
          tt(0) = tt(0) + (x(i1 + 0, l + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + l + 3))
          tt(1) = tt(1) + (x(i1 + 1, l + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + l + 3))
          tt(2) = tt(2) + (x(i1 + 2, l + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + l + 3))
          tt(3) = tt(3) + (x(i1 + 3, l + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + l + 3))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
      end do
      do i2 =  -(-3), n - (3) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        do l = -3, 3, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson6_3_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, l + i2, i3)) * (poisson6_3_fil(l))
          tt(2) = tt(2) + (x(i1 + 2, l + i2, i3)) * (poisson6_3_fil(l))
          tt(3) = tt(3) + (x(i1 + 3, l + i2, i3)) * (poisson6_3_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
      end do
      do i2 = n - (3), n - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        do l = -3, 3, 1
          tt(0) = tt(0) + (x(i1 + 0, l + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + l + 3))
          tt(1) = tt(1) + (x(i1 + 1, l + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + l + 3))
          tt(2) = tt(2) + (x(i1 + 2, l + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + l + 3))
          tt(3) = tt(3) + (x(i1 + 3, l + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + l + 3))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 = 0,  -(-3) - (1), 1
        tt(0) = 0.0_wp
        do l = -3, 3, 1
          tt(0) = tt(0) + (x(i1 + 0, l + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + l + 3))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 =  -(-3), n - (3) - (1), 1
        tt(0) = 0.0_wp
        do l = -3, 3, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson6_3_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 = n - (3), n - (1), 1
        tt(0) = 0.0_wp
        do l = -3, 3, 1
          tt(0) = tt(0) + (x(i1 + 0, l + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + l + 3))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_np_201_u4_0_true_true_false
SUBROUTINE d_poisson6_np_201_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_np_201_u1_2_true_false_true_cost
SUBROUTINE d_poisson6_np_201_a_u4_0_true_false_true(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -3
  integer(kind=4), parameter :: upfil = 3
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(0:48) :: poisson6_fil = (/ &
-2.45_wp, &
6.0_wp, &
-7.5_wp, &
6.666666666666667_wp, &
-3.75_wp, &
1.2_wp, &
-0.16666666666666666_wp, &
-0.16666666666666666_wp, &
-1.2833333333333334_wp, &
2.5_wp, &
-1.6666666666666667_wp, &
0.8333333333333334_wp, &
-0.25_wp, &
0.03333333333333333_wp, &
0.03333333333333333_wp, &
-0.4_wp, &
-0.5833333333333334_wp, &
1.3333333333333333_wp, &
-0.5_wp, &
0.13333333333333333_wp, &
-0.016666666666666666_wp, &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp, &
0.016666666666666666_wp, &
-0.13333333333333333_wp, &
0.5_wp, &
-1.3333333333333333_wp, &
0.5833333333333334_wp, &
0.4_wp, &
-0.03333333333333333_wp, &
-0.03333333333333333_wp, &
0.25_wp, &
-0.8333333333333334_wp, &
1.6666666666666667_wp, &
-2.5_wp, &
1.2833333333333334_wp, &
0.16666666666666666_wp, &
0.16666666666666666_wp, &
-1.2_wp, &
3.75_wp, &
-6.666666666666667_wp, &
7.5_wp, &
-6.0_wp, &
2.45_wp /)
  real(kind=8), parameter, dimension(-3:3) :: poisson6_3_fil = (/ &
-0.016666666666666666_wp, &
0.15_wp, &
-0.75_wp, &
0.0_wp, &
0.75_wp, &
-0.15_wp, &
0.016666666666666666_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  integer(kind=4), dimension(-3 - (3):3 - (-3) - (1)) :: mod_arr
  do l = -3 - (3), 3 - (-3) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1, tt2, tt3)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 = 0,  -(-3) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -3 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -3 + 3))
        tt1 = tt1 + (x(i1 + 1, -3 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -3 + 3))
        tt2 = tt2 + (x(i1 + 2, -3 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -3 + 3))
        tt3 = tt3 + (x(i1 + 3, -3 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -3 + 3))
        tt0 = tt0 + (x(i1 + 0, -2 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -2 + 3))
        tt1 = tt1 + (x(i1 + 1, -2 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -2 + 3))
        tt2 = tt2 + (x(i1 + 2, -2 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -2 + 3))
        tt3 = tt3 + (x(i1 + 3, -2 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -2 + 3))
        tt0 = tt0 + (x(i1 + 0, -1 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -1 + 3))
        tt1 = tt1 + (x(i1 + 1, -1 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -1 + 3))
        tt2 = tt2 + (x(i1 + 2, -1 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -1 + 3))
        tt3 = tt3 + (x(i1 + 3, -1 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -1 + 3))
        tt0 = tt0 + (x(i1 + 0, 0 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 0 + 3))
        tt1 = tt1 + (x(i1 + 1, 0 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 0 + 3))
        tt2 = tt2 + (x(i1 + 2, 0 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 0 + 3))
        tt3 = tt3 + (x(i1 + 3, 0 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 0 + 3))
        tt0 = tt0 + (x(i1 + 0, 1 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 1 + 3))
        tt1 = tt1 + (x(i1 + 1, 1 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 1 + 3))
        tt2 = tt2 + (x(i1 + 2, 1 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 1 + 3))
        tt3 = tt3 + (x(i1 + 3, 1 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 1 + 3))
        tt0 = tt0 + (x(i1 + 0, 2 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 2 + 3))
        tt1 = tt1 + (x(i1 + 1, 2 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 2 + 3))
        tt2 = tt2 + (x(i1 + 2, 2 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 2 + 3))
        tt3 = tt3 + (x(i1 + 3, 2 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 2 + 3))
        tt0 = tt0 + (x(i1 + 0, 3 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 3 + 3))
        tt1 = tt1 + (x(i1 + 1, 3 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 3 + 3))
        tt2 = tt2 + (x(i1 + 2, 3 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 3 + 3))
        tt3 = tt3 + (x(i1 + 3, 3 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 3 + 3))
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
      do i2 =  -(-3), n - (3) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -3 + i2, i3)) * (-0.016666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, -3 + i2, i3)) * (-0.016666666666666666_wp)
        tt2 = tt2 + (x(i1 + 2, -3 + i2, i3)) * (-0.016666666666666666_wp)
        tt3 = tt3 + (x(i1 + 3, -3 + i2, i3)) * (-0.016666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, -2 + i2, i3)) * (0.15_wp)
        tt1 = tt1 + (x(i1 + 1, -2 + i2, i3)) * (0.15_wp)
        tt2 = tt2 + (x(i1 + 2, -2 + i2, i3)) * (0.15_wp)
        tt3 = tt3 + (x(i1 + 3, -2 + i2, i3)) * (0.15_wp)
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.75_wp)
        tt1 = tt1 + (x(i1 + 1, -1 + i2, i3)) * (-0.75_wp)
        tt2 = tt2 + (x(i1 + 2, -1 + i2, i3)) * (-0.75_wp)
        tt3 = tt3 + (x(i1 + 3, -1 + i2, i3)) * (-0.75_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt1 = tt1 + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt2 = tt2 + (x(i1 + 2, 0 + i2, i3)) * (0.0_wp)
        tt3 = tt3 + (x(i1 + 3, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.75_wp)
        tt1 = tt1 + (x(i1 + 1, 1 + i2, i3)) * (0.75_wp)
        tt2 = tt2 + (x(i1 + 2, 1 + i2, i3)) * (0.75_wp)
        tt3 = tt3 + (x(i1 + 3, 1 + i2, i3)) * (0.75_wp)
        tt0 = tt0 + (x(i1 + 0, 2 + i2, i3)) * (-0.15_wp)
        tt1 = tt1 + (x(i1 + 1, 2 + i2, i3)) * (-0.15_wp)
        tt2 = tt2 + (x(i1 + 2, 2 + i2, i3)) * (-0.15_wp)
        tt3 = tt3 + (x(i1 + 3, 2 + i2, i3)) * (-0.15_wp)
        tt0 = tt0 + (x(i1 + 0, 3 + i2, i3)) * (0.016666666666666666_wp)
        tt1 = tt1 + (x(i1 + 1, 3 + i2, i3)) * (0.016666666666666666_wp)
        tt2 = tt2 + (x(i1 + 2, 3 + i2, i3)) * (0.016666666666666666_wp)
        tt3 = tt3 + (x(i1 + 3, 3 + i2, i3)) * (0.016666666666666666_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
      do i2 = n - (3), n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -3 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -3 + 3))
        tt1 = tt1 + (x(i1 + 1, -3 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -3 + 3))
        tt2 = tt2 + (x(i1 + 2, -3 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -3 + 3))
        tt3 = tt3 + (x(i1 + 3, -3 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -3 + 3))
        tt0 = tt0 + (x(i1 + 0, -2 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -2 + 3))
        tt1 = tt1 + (x(i1 + 1, -2 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -2 + 3))
        tt2 = tt2 + (x(i1 + 2, -2 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -2 + 3))
        tt3 = tt3 + (x(i1 + 3, -2 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -2 + 3))
        tt0 = tt0 + (x(i1 + 0, -1 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -1 + 3))
        tt1 = tt1 + (x(i1 + 1, -1 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -1 + 3))
        tt2 = tt2 + (x(i1 + 2, -1 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -1 + 3))
        tt3 = tt3 + (x(i1 + 3, -1 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -1 + 3))
        tt0 = tt0 + (x(i1 + 0, 0 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 0 + 3))
        tt1 = tt1 + (x(i1 + 1, 0 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 0 + 3))
        tt2 = tt2 + (x(i1 + 2, 0 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 0 + 3))
        tt3 = tt3 + (x(i1 + 3, 0 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 0 + 3))
        tt0 = tt0 + (x(i1 + 0, 1 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 1 + 3))
        tt1 = tt1 + (x(i1 + 1, 1 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 1 + 3))
        tt2 = tt2 + (x(i1 + 2, 1 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 1 + 3))
        tt3 = tt3 + (x(i1 + 3, 1 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 1 + 3))
        tt0 = tt0 + (x(i1 + 0, 2 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 2 + 3))
        tt1 = tt1 + (x(i1 + 1, 2 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 2 + 3))
        tt2 = tt2 + (x(i1 + 2, 2 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 2 + 3))
        tt3 = tt3 + (x(i1 + 3, 2 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 2 + 3))
        tt0 = tt0 + (x(i1 + 0, 3 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 3 + 3))
        tt1 = tt1 + (x(i1 + 1, 3 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 3 + 3))
        tt2 = tt2 + (x(i1 + 2, 3 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 3 + 3))
        tt3 = tt3 + (x(i1 + 3, 3 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 3 + 3))
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 = 0,  -(-3) - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -3 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -3 + 3))
        tt0 = tt0 + (x(i1 + 0, -2 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -2 + 3))
        tt0 = tt0 + (x(i1 + 0, -1 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + -1 + 3))
        tt0 = tt0 + (x(i1 + 0, 0 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 0 + 3))
        tt0 = tt0 + (x(i1 + 0, 1 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 1 + 3))
        tt0 = tt0 + (x(i1 + 0, 2 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 2 + 3))
        tt0 = tt0 + (x(i1 + 0, 3 + 3, i3)) * (poisson6_fil((i2 - (3)) * (7) + 21 + 3 + 3))
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-3), n - (3) - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -3 + i2, i3)) * (-0.016666666666666666_wp)
        tt0 = tt0 + (x(i1 + 0, -2 + i2, i3)) * (0.15_wp)
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.75_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.75_wp)
        tt0 = tt0 + (x(i1 + 0, 2 + i2, i3)) * (-0.15_wp)
        tt0 = tt0 + (x(i1 + 0, 3 + i2, i3)) * (0.016666666666666666_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (3), n - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -3 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -3 + 3))
        tt0 = tt0 + (x(i1 + 0, -2 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -2 + 3))
        tt0 = tt0 + (x(i1 + 0, -1 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + -1 + 3))
        tt0 = tt0 + (x(i1 + 0, 0 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 0 + 3))
        tt0 = tt0 + (x(i1 + 0, 1 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 1 + 3))
        tt0 = tt0 + (x(i1 + 0, 2 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 2 + 3))
        tt0 = tt0 + (x(i1 + 0, 3 + -3 + n - (1), i3)) * (poisson6_fil((i2 + 3 - (n) + 1) * (7) + 21 + 3 + 3))
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson6_np_201_a_u4_0_true_false_true
SUBROUTINE d_poisson6_np_201_a_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (7)) * (ndat_t)
END SUBROUTINE d_poisson6_np_201_a_u1_2_true_false_true_cost
SUBROUTINE d_s0s0_1d_poisson6_cost(d, idim, n, bc, x, y, a, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: d
  integer(kind=4), intent(in) :: idim
  integer(kind=4), intent(in), dimension(0:d - (1)) :: n
  integer(kind=4), intent(in) :: bc
  real(kind=8), intent(in), dimension(*) :: x
  real(kind=8), intent(out), dimension(*) :: y
  real(kind=8), intent(in) :: a
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: i
  integer(kind=4) :: ndat_left
  integer(kind=4) :: ndat_right
  integer(kind=4) :: c
  if (idim == 0) then
    ndat_right = 1
    do i = 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson6_p_10_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson6_p_10_a_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson6_fg_10_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson6_fg_10_a_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson6_fs_10_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson6_fs_10_a_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-2)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson6_np_10_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson6_np_10_a_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      end select
  else if (idim == d - (1)) then
    ndat_left = 1
    do i = 0, d - (2), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson6_p_01_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson6_p_01_a_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      case (1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson6_fg_01_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson6_fg_01_a_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      case (-1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson6_fs_01_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson6_fs_01_a_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      case (-2)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson6_np_01_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson6_np_01_a_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      end select
  else
    ndat_left = 1
    ndat_right = 1
    do i = 0, idim - (1), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    do i = idim + 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson6_p_201_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson6_p_201_a_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson6_fg_201_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson6_fg_201_a_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson6_fs_201_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson6_fs_201_a_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-2)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson6_np_201_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson6_np_201_a_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      end select
  end if
END SUBROUTINE d_s0s0_1d_poisson6_cost
SUBROUTINE d_s0s0_1d_poisson6(d, idim, n, bc, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: d
  integer(kind=4), intent(in) :: idim
  integer(kind=4), intent(in), dimension(0:d - (1)) :: n
  integer(kind=4), intent(in) :: bc
  real(kind=8), intent(in), dimension(*) :: x
  real(kind=8), intent(out), dimension(*) :: y
  real(kind=8), intent(in) :: a
  integer(kind=4) :: i
  integer(kind=4) :: ndat_left
  integer(kind=4) :: ndat_right
  if (idim == 0) then
    ndat_right = 1
    do i = 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          call d_poisson6_p_10_u1_1_true_false_true(n(idim), ndat_right, x, y)
        else
          call d_poisson6_p_10_a_u1_1_false_true_false(n(idim), ndat_right, x, y, a)
        end if
      case (1)
        if (a == 1.0_wp) then
          call d_poisson6_fg_10_u1_1_false_false_false(n(idim), ndat_right, x, y)
        else
          call d_poisson6_fg_10_a_u2_1_false_true_false(n(idim), ndat_right, x, y, a)
        end if
      case (-1)
        if (a == 1.0_wp) then
          call d_poisson6_fs_10_u1_1_false_true_false(n(idim), ndat_right, x, y)
        else
          call d_poisson6_fs_10_a_u2_1_false_false_true(n(idim), ndat_right, x, y, a)
        end if
      case (-2)
        if (a == 1.0_wp) then
          call d_poisson6_np_10_u2_1_false_false_false(n(idim), ndat_right, x, y)
        else
          call d_poisson6_np_10_a_u2_1_true_true_false(n(idim), ndat_right, x, y, a)
        end if
      end select
  else if (idim == d - (1)) then
    ndat_left = 1
    do i = 0, d - (2), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          call d_poisson6_p_01_u2_0_true_false_false(ndat_left, n(idim), x, y)
        else
          call d_poisson6_p_01_a_u4_0_false_true_true(ndat_left, n(idim), x, y, a)
        end if
      case (1)
        if (a == 1.0_wp) then
          call d_poisson6_fg_01_u4_0_false_true_true(ndat_left, n(idim), x, y)
        else
          call d_poisson6_fg_01_a_u4_0_false_true_false(ndat_left, n(idim), x, y, a)
        end if
      case (-1)
        if (a == 1.0_wp) then
          call d_poisson6_fs_01_u3_0_false_false_true(ndat_left, n(idim), x, y)
        else
          call d_poisson6_fs_01_a_u2_0_false_true_false(ndat_left, n(idim), x, y, a)
        end if
      case (-2)
        if (a == 1.0_wp) then
          call d_poisson6_np_01_u2_0_false_false_true(ndat_left, n(idim), x, y)
        else
          call d_poisson6_np_01_a_u5_0_false_false_false(ndat_left, n(idim), x, y, a)
        end if
      end select
  else
    ndat_left = 1
    ndat_right = 1
    do i = 0, idim - (1), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    do i = idim + 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          call d_poisson6_p_201_u2_0_true_true_true(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson6_p_201_a_u4_0_false_false_false(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      case (1)
        if (a == 1.0_wp) then
          call d_poisson6_fg_201_u4_0_false_false_false(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson6_fg_201_a_u2_0_false_false_true(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      case (-1)
        if (a == 1.0_wp) then
          call d_poisson6_fs_201_u2_0_false_false_false(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson6_fs_201_a_u4_0_false_true_true(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      case (-2)
        if (a == 1.0_wp) then
          call d_poisson6_np_201_u4_0_true_true_false(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson6_np_201_a_u4_0_true_false_true(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      end select
  end if
END SUBROUTINE d_s0s0_1d_poisson6
SUBROUTINE d_poisson8_p_10_u3_1_true_true_true(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:2) :: tt
  integer(kind=4), dimension(-4 - (4):4 - (-4) - (1)) :: mod_arr
  do l = -4 - (4), 4 - (-4) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (3), 3
    do i1 = 0,  -(-4) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(0) = tt(0) + (x(mod_arr(-4 + i1), i2 + 0)) * (0.0035714285714285713_wp)
      tt(1) = tt(1) + (x(mod_arr(-4 + i1), i2 + 1)) * (0.0035714285714285713_wp)
      tt(2) = tt(2) + (x(mod_arr(-4 + i1), i2 + 2)) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(mod_arr(-3 + i1), i2 + 0)) * (-0.0380952380952381_wp)
      tt(1) = tt(1) + (x(mod_arr(-3 + i1), i2 + 1)) * (-0.0380952380952381_wp)
      tt(2) = tt(2) + (x(mod_arr(-3 + i1), i2 + 2)) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(mod_arr(-2 + i1), i2 + 0)) * (0.2_wp)
      tt(1) = tt(1) + (x(mod_arr(-2 + i1), i2 + 1)) * (0.2_wp)
      tt(2) = tt(2) + (x(mod_arr(-2 + i1), i2 + 2)) * (0.2_wp)
      tt(0) = tt(0) + (x(mod_arr(-1 + i1), i2 + 0)) * (-0.8_wp)
      tt(1) = tt(1) + (x(mod_arr(-1 + i1), i2 + 1)) * (-0.8_wp)
      tt(2) = tt(2) + (x(mod_arr(-1 + i1), i2 + 2)) * (-0.8_wp)
      tt(0) = tt(0) + (x(mod_arr(0 + i1), i2 + 0)) * (0.0_wp)
      tt(1) = tt(1) + (x(mod_arr(0 + i1), i2 + 1)) * (0.0_wp)
      tt(2) = tt(2) + (x(mod_arr(0 + i1), i2 + 2)) * (0.0_wp)
      tt(0) = tt(0) + (x(mod_arr(1 + i1), i2 + 0)) * (0.8_wp)
      tt(1) = tt(1) + (x(mod_arr(1 + i1), i2 + 1)) * (0.8_wp)
      tt(2) = tt(2) + (x(mod_arr(1 + i1), i2 + 2)) * (0.8_wp)
      tt(0) = tt(0) + (x(mod_arr(2 + i1), i2 + 0)) * (-0.2_wp)
      tt(1) = tt(1) + (x(mod_arr(2 + i1), i2 + 1)) * (-0.2_wp)
      tt(2) = tt(2) + (x(mod_arr(2 + i1), i2 + 2)) * (-0.2_wp)
      tt(0) = tt(0) + (x(mod_arr(3 + i1), i2 + 0)) * (0.0380952380952381_wp)
      tt(1) = tt(1) + (x(mod_arr(3 + i1), i2 + 1)) * (0.0380952380952381_wp)
      tt(2) = tt(2) + (x(mod_arr(3 + i1), i2 + 2)) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(mod_arr(4 + i1), i2 + 0)) * (-0.0035714285714285713_wp)
      tt(1) = tt(1) + (x(mod_arr(4 + i1), i2 + 1)) * (-0.0035714285714285713_wp)
      tt(2) = tt(2) + (x(mod_arr(4 + i1), i2 + 2)) * (-0.0035714285714285713_wp)
      y(i1, i2 + 0) = tt(0)
      y(i1, i2 + 1) = tt(1)
      y(i1, i2 + 2) = tt(2)
    end do
    do i1 =  -(-4), n - (4) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(0) = tt(0) + (x(-4 + i1, i2 + 0)) * (0.0035714285714285713_wp)
      tt(1) = tt(1) + (x(-4 + i1, i2 + 1)) * (0.0035714285714285713_wp)
      tt(2) = tt(2) + (x(-4 + i1, i2 + 2)) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(-3 + i1, i2 + 0)) * (-0.0380952380952381_wp)
      tt(1) = tt(1) + (x(-3 + i1, i2 + 1)) * (-0.0380952380952381_wp)
      tt(2) = tt(2) + (x(-3 + i1, i2 + 2)) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(-2 + i1, i2 + 0)) * (0.2_wp)
      tt(1) = tt(1) + (x(-2 + i1, i2 + 1)) * (0.2_wp)
      tt(2) = tt(2) + (x(-2 + i1, i2 + 2)) * (0.2_wp)
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.8_wp)
      tt(1) = tt(1) + (x(-1 + i1, i2 + 1)) * (-0.8_wp)
      tt(2) = tt(2) + (x(-1 + i1, i2 + 2)) * (-0.8_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(1) = tt(1) + (x(0 + i1, i2 + 1)) * (0.0_wp)
      tt(2) = tt(2) + (x(0 + i1, i2 + 2)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.8_wp)
      tt(1) = tt(1) + (x(1 + i1, i2 + 1)) * (0.8_wp)
      tt(2) = tt(2) + (x(1 + i1, i2 + 2)) * (0.8_wp)
      tt(0) = tt(0) + (x(2 + i1, i2 + 0)) * (-0.2_wp)
      tt(1) = tt(1) + (x(2 + i1, i2 + 1)) * (-0.2_wp)
      tt(2) = tt(2) + (x(2 + i1, i2 + 2)) * (-0.2_wp)
      tt(0) = tt(0) + (x(3 + i1, i2 + 0)) * (0.0380952380952381_wp)
      tt(1) = tt(1) + (x(3 + i1, i2 + 1)) * (0.0380952380952381_wp)
      tt(2) = tt(2) + (x(3 + i1, i2 + 2)) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(4 + i1, i2 + 0)) * (-0.0035714285714285713_wp)
      tt(1) = tt(1) + (x(4 + i1, i2 + 1)) * (-0.0035714285714285713_wp)
      tt(2) = tt(2) + (x(4 + i1, i2 + 2)) * (-0.0035714285714285713_wp)
      y(i1, i2 + 0) = tt(0)
      y(i1, i2 + 1) = tt(1)
      y(i1, i2 + 2) = tt(2)
    end do
    do i1 = n - (4), n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(0) = tt(0) + (x(mod_arr(-4 + i1 - (n)), i2 + 0)) * (0.0035714285714285713_wp)
      tt(1) = tt(1) + (x(mod_arr(-4 + i1 - (n)), i2 + 1)) * (0.0035714285714285713_wp)
      tt(2) = tt(2) + (x(mod_arr(-4 + i1 - (n)), i2 + 2)) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(mod_arr(-3 + i1 - (n)), i2 + 0)) * (-0.0380952380952381_wp)
      tt(1) = tt(1) + (x(mod_arr(-3 + i1 - (n)), i2 + 1)) * (-0.0380952380952381_wp)
      tt(2) = tt(2) + (x(mod_arr(-3 + i1 - (n)), i2 + 2)) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(mod_arr(-2 + i1 - (n)), i2 + 0)) * (0.2_wp)
      tt(1) = tt(1) + (x(mod_arr(-2 + i1 - (n)), i2 + 1)) * (0.2_wp)
      tt(2) = tt(2) + (x(mod_arr(-2 + i1 - (n)), i2 + 2)) * (0.2_wp)
      tt(0) = tt(0) + (x(mod_arr(-1 + i1 - (n)), i2 + 0)) * (-0.8_wp)
      tt(1) = tt(1) + (x(mod_arr(-1 + i1 - (n)), i2 + 1)) * (-0.8_wp)
      tt(2) = tt(2) + (x(mod_arr(-1 + i1 - (n)), i2 + 2)) * (-0.8_wp)
      tt(0) = tt(0) + (x(mod_arr(0 + i1 - (n)), i2 + 0)) * (0.0_wp)
      tt(1) = tt(1) + (x(mod_arr(0 + i1 - (n)), i2 + 1)) * (0.0_wp)
      tt(2) = tt(2) + (x(mod_arr(0 + i1 - (n)), i2 + 2)) * (0.0_wp)
      tt(0) = tt(0) + (x(mod_arr(1 + i1 - (n)), i2 + 0)) * (0.8_wp)
      tt(1) = tt(1) + (x(mod_arr(1 + i1 - (n)), i2 + 1)) * (0.8_wp)
      tt(2) = tt(2) + (x(mod_arr(1 + i1 - (n)), i2 + 2)) * (0.8_wp)
      tt(0) = tt(0) + (x(mod_arr(2 + i1 - (n)), i2 + 0)) * (-0.2_wp)
      tt(1) = tt(1) + (x(mod_arr(2 + i1 - (n)), i2 + 1)) * (-0.2_wp)
      tt(2) = tt(2) + (x(mod_arr(2 + i1 - (n)), i2 + 2)) * (-0.2_wp)
      tt(0) = tt(0) + (x(mod_arr(3 + i1 - (n)), i2 + 0)) * (0.0380952380952381_wp)
      tt(1) = tt(1) + (x(mod_arr(3 + i1 - (n)), i2 + 1)) * (0.0380952380952381_wp)
      tt(2) = tt(2) + (x(mod_arr(3 + i1 - (n)), i2 + 2)) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(mod_arr(4 + i1 - (n)), i2 + 0)) * (-0.0035714285714285713_wp)
      tt(1) = tt(1) + (x(mod_arr(4 + i1 - (n)), i2 + 1)) * (-0.0035714285714285713_wp)
      tt(2) = tt(2) + (x(mod_arr(4 + i1 - (n)), i2 + 2)) * (-0.0035714285714285713_wp)
      y(i1, i2 + 0) = tt(0)
      y(i1, i2 + 1) = tt(1)
      y(i1, i2 + 2) = tt(2)
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (3)) * (3), ndat - (1), 1
    do i1 = 0,  -(-4) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(mod_arr(-4 + i1), i2 + 0)) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(mod_arr(-3 + i1), i2 + 0)) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(mod_arr(-2 + i1), i2 + 0)) * (0.2_wp)
      tt(0) = tt(0) + (x(mod_arr(-1 + i1), i2 + 0)) * (-0.8_wp)
      tt(0) = tt(0) + (x(mod_arr(0 + i1), i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(mod_arr(1 + i1), i2 + 0)) * (0.8_wp)
      tt(0) = tt(0) + (x(mod_arr(2 + i1), i2 + 0)) * (-0.2_wp)
      tt(0) = tt(0) + (x(mod_arr(3 + i1), i2 + 0)) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(mod_arr(4 + i1), i2 + 0)) * (-0.0035714285714285713_wp)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-4), n - (4) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-4 + i1, i2 + 0)) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(-3 + i1, i2 + 0)) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(-2 + i1, i2 + 0)) * (0.2_wp)
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.8_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.8_wp)
      tt(0) = tt(0) + (x(2 + i1, i2 + 0)) * (-0.2_wp)
      tt(0) = tt(0) + (x(3 + i1, i2 + 0)) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(4 + i1, i2 + 0)) * (-0.0035714285714285713_wp)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (4), n - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(mod_arr(-4 + i1 - (n)), i2 + 0)) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(mod_arr(-3 + i1 - (n)), i2 + 0)) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(mod_arr(-2 + i1 - (n)), i2 + 0)) * (0.2_wp)
      tt(0) = tt(0) + (x(mod_arr(-1 + i1 - (n)), i2 + 0)) * (-0.8_wp)
      tt(0) = tt(0) + (x(mod_arr(0 + i1 - (n)), i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(mod_arr(1 + i1 - (n)), i2 + 0)) * (0.8_wp)
      tt(0) = tt(0) + (x(mod_arr(2 + i1 - (n)), i2 + 0)) * (-0.2_wp)
      tt(0) = tt(0) + (x(mod_arr(3 + i1 - (n)), i2 + 0)) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(mod_arr(4 + i1 - (n)), i2 + 0)) * (-0.0035714285714285713_wp)
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_p_10_u3_1_true_true_true
SUBROUTINE d_poisson8_p_10_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_p_10_u1_1_true_false_true_cost
SUBROUTINE d_poisson8_p_10_a_u1_1_true_true_true(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:0) :: tt
  integer(kind=4), dimension(-4 - (4):4 - (-4) - (1)) :: mod_arr
  do l = -4 - (4), 4 - (-4) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 = 0,  -(-4) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(mod_arr(-4 + i1), i2 + 0)) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(mod_arr(-3 + i1), i2 + 0)) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(mod_arr(-2 + i1), i2 + 0)) * (0.2_wp)
      tt(0) = tt(0) + (x(mod_arr(-1 + i1), i2 + 0)) * (-0.8_wp)
      tt(0) = tt(0) + (x(mod_arr(0 + i1), i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(mod_arr(1 + i1), i2 + 0)) * (0.8_wp)
      tt(0) = tt(0) + (x(mod_arr(2 + i1), i2 + 0)) * (-0.2_wp)
      tt(0) = tt(0) + (x(mod_arr(3 + i1), i2 + 0)) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(mod_arr(4 + i1), i2 + 0)) * (-0.0035714285714285713_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-4), n - (4) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-4 + i1, i2 + 0)) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(-3 + i1, i2 + 0)) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(-2 + i1, i2 + 0)) * (0.2_wp)
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.8_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.8_wp)
      tt(0) = tt(0) + (x(2 + i1, i2 + 0)) * (-0.2_wp)
      tt(0) = tt(0) + (x(3 + i1, i2 + 0)) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(4 + i1, i2 + 0)) * (-0.0035714285714285713_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (4), n - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(mod_arr(-4 + i1 - (n)), i2 + 0)) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(mod_arr(-3 + i1 - (n)), i2 + 0)) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(mod_arr(-2 + i1 - (n)), i2 + 0)) * (0.2_wp)
      tt(0) = tt(0) + (x(mod_arr(-1 + i1 - (n)), i2 + 0)) * (-0.8_wp)
      tt(0) = tt(0) + (x(mod_arr(0 + i1 - (n)), i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(mod_arr(1 + i1 - (n)), i2 + 0)) * (0.8_wp)
      tt(0) = tt(0) + (x(mod_arr(2 + i1 - (n)), i2 + 0)) * (-0.2_wp)
      tt(0) = tt(0) + (x(mod_arr(3 + i1 - (n)), i2 + 0)) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(mod_arr(4 + i1 - (n)), i2 + 0)) * (-0.0035714285714285713_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_p_10_a_u1_1_true_true_true
SUBROUTINE d_poisson8_p_10_a_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_p_10_a_u1_1_true_false_true_cost
SUBROUTINE d_poisson8_p_01_u2_0_true_true_true(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
  integer(kind=4), dimension(-4 - (4):4 - (-4) - (1)) :: mod_arr
  do l = -4 - (4), 4 - (-4) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (2), 2
    do i2 = 0,  -(-4) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(-4 + i2))) * (0.0035714285714285713_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(-4 + i2))) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(-3 + i2))) * (-0.0380952380952381_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(-3 + i2))) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(-2 + i2))) * (0.2_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(-2 + i2))) * (0.2_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(-1 + i2))) * (-0.8_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(-1 + i2))) * (-0.8_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(0 + i2))) * (0.0_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(0 + i2))) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(1 + i2))) * (0.8_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(1 + i2))) * (0.8_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(2 + i2))) * (-0.2_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(2 + i2))) * (-0.2_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(3 + i2))) * (0.0380952380952381_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(3 + i2))) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(4 + i2))) * (-0.0035714285714285713_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(4 + i2))) * (-0.0035714285714285713_wp)
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
    end do
    do i2 =  -(-4), n - (4) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -4 + i2)) * (0.0035714285714285713_wp)
      tt(1) = tt(1) + (x(i1 + 1, -4 + i2)) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(i1 + 0, -3 + i2)) * (-0.0380952380952381_wp)
      tt(1) = tt(1) + (x(i1 + 1, -3 + i2)) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(i1 + 0, -2 + i2)) * (0.2_wp)
      tt(1) = tt(1) + (x(i1 + 1, -2 + i2)) * (0.2_wp)
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2)) * (-0.8_wp)
      tt(1) = tt(1) + (x(i1 + 1, -1 + i2)) * (-0.8_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt(1) = tt(1) + (x(i1 + 1, 0 + i2)) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2)) * (0.8_wp)
      tt(1) = tt(1) + (x(i1 + 1, 1 + i2)) * (0.8_wp)
      tt(0) = tt(0) + (x(i1 + 0, 2 + i2)) * (-0.2_wp)
      tt(1) = tt(1) + (x(i1 + 1, 2 + i2)) * (-0.2_wp)
      tt(0) = tt(0) + (x(i1 + 0, 3 + i2)) * (0.0380952380952381_wp)
      tt(1) = tt(1) + (x(i1 + 1, 3 + i2)) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(i1 + 0, 4 + i2)) * (-0.0035714285714285713_wp)
      tt(1) = tt(1) + (x(i1 + 1, 4 + i2)) * (-0.0035714285714285713_wp)
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
    end do
    do i2 = n - (4), n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(-4 + i2 - (n)))) * (0.0035714285714285713_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(-4 + i2 - (n)))) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(-3 + i2 - (n)))) * (-0.0380952380952381_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(-3 + i2 - (n)))) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(-2 + i2 - (n)))) * (0.2_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(-2 + i2 - (n)))) * (0.2_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(-1 + i2 - (n)))) * (-0.8_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(-1 + i2 - (n)))) * (-0.8_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(0 + i2 - (n)))) * (0.0_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(0 + i2 - (n)))) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(1 + i2 - (n)))) * (0.8_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(1 + i2 - (n)))) * (0.8_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(2 + i2 - (n)))) * (-0.2_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(2 + i2 - (n)))) * (-0.2_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(3 + i2 - (n)))) * (0.0380952380952381_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(3 + i2 - (n)))) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(4 + i2 - (n)))) * (-0.0035714285714285713_wp)
      tt(1) = tt(1) + (x(i1 + 1, mod_arr(4 + i2 - (n)))) * (-0.0035714285714285713_wp)
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i2 = 0,  -(-4) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(-4 + i2))) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(-3 + i2))) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(-2 + i2))) * (0.2_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(-1 + i2))) * (-0.8_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(0 + i2))) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(1 + i2))) * (0.8_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(2 + i2))) * (-0.2_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(3 + i2))) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(4 + i2))) * (-0.0035714285714285713_wp)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 =  -(-4), n - (4) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -4 + i2)) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(i1 + 0, -3 + i2)) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(i1 + 0, -2 + i2)) * (0.2_wp)
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2)) * (-0.8_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2)) * (0.8_wp)
      tt(0) = tt(0) + (x(i1 + 0, 2 + i2)) * (-0.2_wp)
      tt(0) = tt(0) + (x(i1 + 0, 3 + i2)) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(i1 + 0, 4 + i2)) * (-0.0035714285714285713_wp)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 = n - (4), n - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(-4 + i2 - (n)))) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(-3 + i2 - (n)))) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(-2 + i2 - (n)))) * (0.2_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(-1 + i2 - (n)))) * (-0.8_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(0 + i2 - (n)))) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(1 + i2 - (n)))) * (0.8_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(2 + i2 - (n)))) * (-0.2_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(3 + i2 - (n)))) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(i1 + 0, mod_arr(4 + i2 - (n)))) * (-0.0035714285714285713_wp)
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_p_01_u2_0_true_true_true
SUBROUTINE d_poisson8_p_01_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_p_01_u1_0_true_false_true_cost
SUBROUTINE d_poisson8_p_01_a_u2_0_false_true_false(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (2), 2
    do i2 = 0,  -(-4) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson8_4_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson8_4_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
    end do
    do i2 =  -(-4), n - (4) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson8_4_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson8_4_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
    end do
    do i2 = n - (4), n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson8_4_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson8_4_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i2 = 0,  -(-4) - (1), 1
      tt(0) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson8_4_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 =  -(-4), n - (4) - (1), 1
      tt(0) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson8_4_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 = n - (4), n - (1), 1
      tt(0) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2 - (((i2 + l + (n) * (2)) / (n) - (2)) * (n)))) * (poisson8_4_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_p_01_a_u2_0_false_true_false
SUBROUTINE d_poisson8_p_01_a_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_p_01_a_u1_0_true_false_true_cost
SUBROUTINE d_poisson8_p_201_u2_0_true_false_true(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  integer(kind=4), dimension(-4 - (4):4 - (-4) - (1)) :: mod_arr
  do l = -4 - (4), 4 - (-4) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (2), 2
      do i2 = 0,  -(-4) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, mod_arr(-4 + i2), i3)) * (0.0035714285714285713_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(-4 + i2), i3)) * (0.0035714285714285713_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-3 + i2), i3)) * (-0.0380952380952381_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(-3 + i2), i3)) * (-0.0380952380952381_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-2 + i2), i3)) * (0.2_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(-2 + i2), i3)) * (0.2_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-1 + i2), i3)) * (-0.8_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(-1 + i2), i3)) * (-0.8_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(1 + i2), i3)) * (0.8_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(1 + i2), i3)) * (0.8_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(2 + i2), i3)) * (-0.2_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(2 + i2), i3)) * (-0.2_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(3 + i2), i3)) * (0.0380952380952381_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(3 + i2), i3)) * (0.0380952380952381_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(4 + i2), i3)) * (-0.0035714285714285713_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(4 + i2), i3)) * (-0.0035714285714285713_wp)
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
      end do
      do i2 =  -(-4), n - (4) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -4 + i2, i3)) * (0.0035714285714285713_wp)
        tt1 = tt1 + (x(i1 + 1, -4 + i2, i3)) * (0.0035714285714285713_wp)
        tt0 = tt0 + (x(i1 + 0, -3 + i2, i3)) * (-0.0380952380952381_wp)
        tt1 = tt1 + (x(i1 + 1, -3 + i2, i3)) * (-0.0380952380952381_wp)
        tt0 = tt0 + (x(i1 + 0, -2 + i2, i3)) * (0.2_wp)
        tt1 = tt1 + (x(i1 + 1, -2 + i2, i3)) * (0.2_wp)
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.8_wp)
        tt1 = tt1 + (x(i1 + 1, -1 + i2, i3)) * (-0.8_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt1 = tt1 + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.8_wp)
        tt1 = tt1 + (x(i1 + 1, 1 + i2, i3)) * (0.8_wp)
        tt0 = tt0 + (x(i1 + 0, 2 + i2, i3)) * (-0.2_wp)
        tt1 = tt1 + (x(i1 + 1, 2 + i2, i3)) * (-0.2_wp)
        tt0 = tt0 + (x(i1 + 0, 3 + i2, i3)) * (0.0380952380952381_wp)
        tt1 = tt1 + (x(i1 + 1, 3 + i2, i3)) * (0.0380952380952381_wp)
        tt0 = tt0 + (x(i1 + 0, 4 + i2, i3)) * (-0.0035714285714285713_wp)
        tt1 = tt1 + (x(i1 + 1, 4 + i2, i3)) * (-0.0035714285714285713_wp)
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
      end do
      do i2 = n - (4), n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, mod_arr(-4 + i2 - (n)), i3)) * (0.0035714285714285713_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(-4 + i2 - (n)), i3)) * (0.0035714285714285713_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-3 + i2 - (n)), i3)) * (-0.0380952380952381_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(-3 + i2 - (n)), i3)) * (-0.0380952380952381_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-2 + i2 - (n)), i3)) * (0.2_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(-2 + i2 - (n)), i3)) * (0.2_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-1 + i2 - (n)), i3)) * (-0.8_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(-1 + i2 - (n)), i3)) * (-0.8_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(1 + i2 - (n)), i3)) * (0.8_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(1 + i2 - (n)), i3)) * (0.8_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(2 + i2 - (n)), i3)) * (-0.2_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(2 + i2 - (n)), i3)) * (-0.2_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(3 + i2 - (n)), i3)) * (0.0380952380952381_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(3 + i2 - (n)), i3)) * (0.0380952380952381_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(4 + i2 - (n)), i3)) * (-0.0035714285714285713_wp)
        tt1 = tt1 + (x(i1 + 1, mod_arr(4 + i2 - (n)), i3)) * (-0.0035714285714285713_wp)
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (2)) * (2), ndat1 - (1), 1
      do i2 = 0,  -(-4) - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, mod_arr(-4 + i2), i3)) * (0.0035714285714285713_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-3 + i2), i3)) * (-0.0380952380952381_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-2 + i2), i3)) * (0.2_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-1 + i2), i3)) * (-0.8_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(0 + i2), i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(1 + i2), i3)) * (0.8_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(2 + i2), i3)) * (-0.2_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(3 + i2), i3)) * (0.0380952380952381_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(4 + i2), i3)) * (-0.0035714285714285713_wp)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-4), n - (4) - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -4 + i2, i3)) * (0.0035714285714285713_wp)
        tt0 = tt0 + (x(i1 + 0, -3 + i2, i3)) * (-0.0380952380952381_wp)
        tt0 = tt0 + (x(i1 + 0, -2 + i2, i3)) * (0.2_wp)
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.8_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.8_wp)
        tt0 = tt0 + (x(i1 + 0, 2 + i2, i3)) * (-0.2_wp)
        tt0 = tt0 + (x(i1 + 0, 3 + i2, i3)) * (0.0380952380952381_wp)
        tt0 = tt0 + (x(i1 + 0, 4 + i2, i3)) * (-0.0035714285714285713_wp)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (4), n - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, mod_arr(-4 + i2 - (n)), i3)) * (0.0035714285714285713_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-3 + i2 - (n)), i3)) * (-0.0380952380952381_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-2 + i2 - (n)), i3)) * (0.2_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(-1 + i2 - (n)), i3)) * (-0.8_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(0 + i2 - (n)), i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(1 + i2 - (n)), i3)) * (0.8_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(2 + i2 - (n)), i3)) * (-0.2_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(3 + i2 - (n)), i3)) * (0.0380952380952381_wp)
        tt0 = tt0 + (x(i1 + 0, mod_arr(4 + i2 - (n)), i3)) * (-0.0035714285714285713_wp)
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_p_201_u2_0_true_false_true
SUBROUTINE d_poisson8_p_201_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_p_201_u1_2_true_false_true_cost
SUBROUTINE d_poisson8_p_201_a_u2_0_true_false_false(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  integer(kind=4), dimension(-4 - (4):4 - (-4) - (1)) :: mod_arr
  do l = -4 - (4), 4 - (-4) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (2), 2
      do i2 = 0,  -(-4) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2), i3)) * (poisson8_4_fil(l))
          tt1 = tt1 + (x(i1 + 1, mod_arr(l + i2), i3)) * (poisson8_4_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
      end do
      do i2 =  -(-4), n - (4) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson8_4_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson8_4_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
      end do
      do i2 = n - (4), n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2 - (n)), i3)) * (poisson8_4_fil(l))
          tt1 = tt1 + (x(i1 + 1, mod_arr(l + i2 - (n)), i3)) * (poisson8_4_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (2)) * (2), ndat1 - (1), 1
      do i2 = 0,  -(-4) - (1), 1
        tt0 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2), i3)) * (poisson8_4_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-4), n - (4) - (1), 1
        tt0 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson8_4_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (4), n - (1), 1
        tt0 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2 - (n)), i3)) * (poisson8_4_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_p_201_a_u2_0_true_false_false
SUBROUTINE d_poisson8_p_201_a_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_p_201_a_u1_2_true_false_true_cost
SUBROUTINE d_poisson8_fg_10_u1_1_false_true_true(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension( -(4):n - (-4) - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:0) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 =  -(4),  -(-4) - (1), 1
      tt(0) = 0.0_wp
      do l = max( -(i1), -4), 4, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson8_4_fil(l))
      end do
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-4), n - (4) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-4 + i1, i2 + 0)) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(-3 + i1, i2 + 0)) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(-2 + i1, i2 + 0)) * (0.2_wp)
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.8_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.8_wp)
      tt(0) = tt(0) + (x(2 + i1, i2 + 0)) * (-0.2_wp)
      tt(0) = tt(0) + (x(3 + i1, i2 + 0)) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(4 + i1, i2 + 0)) * (-0.0035714285714285713_wp)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (4), n - (-4) - (1), 1
      tt(0) = 0.0_wp
      do l = -4, min(4, n - (1) - (i1)), 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson8_4_fil(l))
      end do
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_fg_10_u1_1_false_true_true
SUBROUTINE d_poisson8_fg_10_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_fg_10_u1_1_false_false_true_cost
SUBROUTINE d_poisson8_fg_10_a_u2_1_false_true_false(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension( -(4):n - (-4) - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 =  -(4),  -(-4) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = max( -(i1), -4), 4, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson8_4_fil(l))
        tt(1) = tt(1) + (x(l + i1, i2 + 1)) * (poisson8_4_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
    do i1 =  -(-4), n - (4) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson8_4_fil(l))
        tt(1) = tt(1) + (x(l + i1, i2 + 1)) * (poisson8_4_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
    do i1 = n - (4), n - (-4) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -4, min(4, n - (1) - (i1)), 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson8_4_fil(l))
        tt(1) = tt(1) + (x(l + i1, i2 + 1)) * (poisson8_4_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 =  -(4),  -(-4) - (1), 1
      tt(0) = 0.0_wp
      do l = max( -(i1), -4), 4, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson8_4_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-4), n - (4) - (1), 1
      tt(0) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson8_4_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (4), n - (-4) - (1), 1
      tt(0) = 0.0_wp
      do l = -4, min(4, n - (1) - (i1)), 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson8_4_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_fg_10_a_u2_1_false_true_false
SUBROUTINE d_poisson8_fg_10_a_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_fg_10_a_u1_1_false_false_true_cost
SUBROUTINE d_poisson8_fg_01_u2_0_false_false_true(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1),  -(4):n - (-4) - (1)) :: y
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
!$omp parallel  default(shared) private(i1, i2, tt0, tt1)
!$omp do 
  do i1 = 0, ndat - (2), 2
    do i2 =  -(4),  -(-4) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = max( -(i2), -4), 4, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson8_4_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson8_4_fil(l))
      end do
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
    end do
    do i2 =  -(-4), n - (4) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -4 + i2)) * (0.0035714285714285713_wp)
      tt1 = tt1 + (x(i1 + 1, -4 + i2)) * (0.0035714285714285713_wp)
      tt0 = tt0 + (x(i1 + 0, -3 + i2)) * (-0.0380952380952381_wp)
      tt1 = tt1 + (x(i1 + 1, -3 + i2)) * (-0.0380952380952381_wp)
      tt0 = tt0 + (x(i1 + 0, -2 + i2)) * (0.2_wp)
      tt1 = tt1 + (x(i1 + 1, -2 + i2)) * (0.2_wp)
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.8_wp)
      tt1 = tt1 + (x(i1 + 1, -1 + i2)) * (-0.8_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt1 = tt1 + (x(i1 + 1, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.8_wp)
      tt1 = tt1 + (x(i1 + 1, 1 + i2)) * (0.8_wp)
      tt0 = tt0 + (x(i1 + 0, 2 + i2)) * (-0.2_wp)
      tt1 = tt1 + (x(i1 + 1, 2 + i2)) * (-0.2_wp)
      tt0 = tt0 + (x(i1 + 0, 3 + i2)) * (0.0380952380952381_wp)
      tt1 = tt1 + (x(i1 + 1, 3 + i2)) * (0.0380952380952381_wp)
      tt0 = tt0 + (x(i1 + 0, 4 + i2)) * (-0.0035714285714285713_wp)
      tt1 = tt1 + (x(i1 + 1, 4 + i2)) * (-0.0035714285714285713_wp)
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
    end do
    do i2 = n - (4), n - (-4) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -4, min(4, n - (1) - (i2)), 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson8_4_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson8_4_fil(l))
      end do
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i2 =  -(4),  -(-4) - (1), 1
      tt0 = 0.0_wp
      do l = max( -(i2), -4), 4, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson8_4_fil(l))
      end do
      y(i1 + 0, i2) = tt0
    end do
    do i2 =  -(-4), n - (4) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -4 + i2)) * (0.0035714285714285713_wp)
      tt0 = tt0 + (x(i1 + 0, -3 + i2)) * (-0.0380952380952381_wp)
      tt0 = tt0 + (x(i1 + 0, -2 + i2)) * (0.2_wp)
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.8_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.8_wp)
      tt0 = tt0 + (x(i1 + 0, 2 + i2)) * (-0.2_wp)
      tt0 = tt0 + (x(i1 + 0, 3 + i2)) * (0.0380952380952381_wp)
      tt0 = tt0 + (x(i1 + 0, 4 + i2)) * (-0.0035714285714285713_wp)
      y(i1 + 0, i2) = tt0
    end do
    do i2 = n - (4), n - (-4) - (1), 1
      tt0 = 0.0_wp
      do l = -4, min(4, n - (1) - (i2)), 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson8_4_fil(l))
      end do
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_fg_01_u2_0_false_false_true
SUBROUTINE d_poisson8_fg_01_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_fg_01_u1_0_false_false_true_cost
SUBROUTINE d_poisson8_fg_01_a_u3_0_false_false_false(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1),  -(4):n - (-4) - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2)
!$omp do 
  do i1 = 0, ndat - (3), 3
    do i2 =  -(4),  -(-4) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      do l = max( -(i2), -4), 4, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson8_4_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson8_4_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2)) * (poisson8_4_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
    end do
    do i2 =  -(-4), n - (4) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      do l = -4, 4, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson8_4_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson8_4_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2)) * (poisson8_4_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
    end do
    do i2 = n - (4), n - (-4) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      do l = -4, min(4, n - (1) - (i2)), 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson8_4_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson8_4_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2)) * (poisson8_4_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (3)) * (3), ndat - (1), 1
    do i2 =  -(4),  -(-4) - (1), 1
      tt0 = 0.0_wp
      do l = max( -(i2), -4), 4, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson8_4_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
    do i2 =  -(-4), n - (4) - (1), 1
      tt0 = 0.0_wp
      do l = -4, 4, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson8_4_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
    do i2 = n - (4), n - (-4) - (1), 1
      tt0 = 0.0_wp
      do l = -4, min(4, n - (1) - (i2)), 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson8_4_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_fg_01_a_u3_0_false_false_false
SUBROUTINE d_poisson8_fg_01_a_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_fg_01_a_u1_0_false_false_true_cost
SUBROUTINE d_poisson8_fg_201_u2_0_false_true_true(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1),  -(4):n - (-4) - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
!$omp parallel  default(shared) private(i1, i2, i3, tt)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (2), 2
      do i2 =  -(4),  -(-4) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        do l = max( -(i2), -4), 4, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson8_4_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, l + i2, i3)) * (poisson8_4_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
      end do
      do i2 =  -(-4), n - (4) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -4 + i2, i3)) * (0.0035714285714285713_wp)
        tt(1) = tt(1) + (x(i1 + 1, -4 + i2, i3)) * (0.0035714285714285713_wp)
        tt(0) = tt(0) + (x(i1 + 0, -3 + i2, i3)) * (-0.0380952380952381_wp)
        tt(1) = tt(1) + (x(i1 + 1, -3 + i2, i3)) * (-0.0380952380952381_wp)
        tt(0) = tt(0) + (x(i1 + 0, -2 + i2, i3)) * (0.2_wp)
        tt(1) = tt(1) + (x(i1 + 1, -2 + i2, i3)) * (0.2_wp)
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.8_wp)
        tt(1) = tt(1) + (x(i1 + 1, -1 + i2, i3)) * (-0.8_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(1) = tt(1) + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.8_wp)
        tt(1) = tt(1) + (x(i1 + 1, 1 + i2, i3)) * (0.8_wp)
        tt(0) = tt(0) + (x(i1 + 0, 2 + i2, i3)) * (-0.2_wp)
        tt(1) = tt(1) + (x(i1 + 1, 2 + i2, i3)) * (-0.2_wp)
        tt(0) = tt(0) + (x(i1 + 0, 3 + i2, i3)) * (0.0380952380952381_wp)
        tt(1) = tt(1) + (x(i1 + 1, 3 + i2, i3)) * (0.0380952380952381_wp)
        tt(0) = tt(0) + (x(i1 + 0, 4 + i2, i3)) * (-0.0035714285714285713_wp)
        tt(1) = tt(1) + (x(i1 + 1, 4 + i2, i3)) * (-0.0035714285714285713_wp)
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
      end do
      do i2 = n - (4), n - (-4) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        do l = -4, min(4, n - (1) - (i2)), 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson8_4_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, l + i2, i3)) * (poisson8_4_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (2)) * (2), ndat1 - (1), 1
      do i2 =  -(4),  -(-4) - (1), 1
        tt(0) = 0.0_wp
        do l = max( -(i2), -4), 4, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson8_4_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 =  -(-4), n - (4) - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -4 + i2, i3)) * (0.0035714285714285713_wp)
        tt(0) = tt(0) + (x(i1 + 0, -3 + i2, i3)) * (-0.0380952380952381_wp)
        tt(0) = tt(0) + (x(i1 + 0, -2 + i2, i3)) * (0.2_wp)
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.8_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.8_wp)
        tt(0) = tt(0) + (x(i1 + 0, 2 + i2, i3)) * (-0.2_wp)
        tt(0) = tt(0) + (x(i1 + 0, 3 + i2, i3)) * (0.0380952380952381_wp)
        tt(0) = tt(0) + (x(i1 + 0, 4 + i2, i3)) * (-0.0035714285714285713_wp)
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 = n - (4), n - (-4) - (1), 1
        tt(0) = 0.0_wp
        do l = -4, min(4, n - (1) - (i2)), 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson8_4_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_fg_201_u2_0_false_true_true
SUBROUTINE d_poisson8_fg_201_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_fg_201_u1_2_false_false_true_cost
SUBROUTINE d_poisson8_fg_201_a_u2_0_false_false_true(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1),  -(4):n - (-4) - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (2), 2
      do i2 =  -(4),  -(-4) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        do l = max( -(i2), -4), 4, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson8_4_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson8_4_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
      end do
      do i2 =  -(-4), n - (4) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -4 + i2, i3)) * (0.0035714285714285713_wp)
        tt1 = tt1 + (x(i1 + 1, -4 + i2, i3)) * (0.0035714285714285713_wp)
        tt0 = tt0 + (x(i1 + 0, -3 + i2, i3)) * (-0.0380952380952381_wp)
        tt1 = tt1 + (x(i1 + 1, -3 + i2, i3)) * (-0.0380952380952381_wp)
        tt0 = tt0 + (x(i1 + 0, -2 + i2, i3)) * (0.2_wp)
        tt1 = tt1 + (x(i1 + 1, -2 + i2, i3)) * (0.2_wp)
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.8_wp)
        tt1 = tt1 + (x(i1 + 1, -1 + i2, i3)) * (-0.8_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt1 = tt1 + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.8_wp)
        tt1 = tt1 + (x(i1 + 1, 1 + i2, i3)) * (0.8_wp)
        tt0 = tt0 + (x(i1 + 0, 2 + i2, i3)) * (-0.2_wp)
        tt1 = tt1 + (x(i1 + 1, 2 + i2, i3)) * (-0.2_wp)
        tt0 = tt0 + (x(i1 + 0, 3 + i2, i3)) * (0.0380952380952381_wp)
        tt1 = tt1 + (x(i1 + 1, 3 + i2, i3)) * (0.0380952380952381_wp)
        tt0 = tt0 + (x(i1 + 0, 4 + i2, i3)) * (-0.0035714285714285713_wp)
        tt1 = tt1 + (x(i1 + 1, 4 + i2, i3)) * (-0.0035714285714285713_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
      end do
      do i2 = n - (4), n - (-4) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        do l = -4, min(4, n - (1) - (i2)), 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson8_4_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson8_4_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (2)) * (2), ndat1 - (1), 1
      do i2 =  -(4),  -(-4) - (1), 1
        tt0 = 0.0_wp
        do l = max( -(i2), -4), 4, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson8_4_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-4), n - (4) - (1), 1
        tt0 = 0.0_wp
        tt0 = tt0 + (x(i1 + 0, -4 + i2, i3)) * (0.0035714285714285713_wp)
        tt0 = tt0 + (x(i1 + 0, -3 + i2, i3)) * (-0.0380952380952381_wp)
        tt0 = tt0 + (x(i1 + 0, -2 + i2, i3)) * (0.2_wp)
        tt0 = tt0 + (x(i1 + 0, -1 + i2, i3)) * (-0.8_wp)
        tt0 = tt0 + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt0 = tt0 + (x(i1 + 0, 1 + i2, i3)) * (0.8_wp)
        tt0 = tt0 + (x(i1 + 0, 2 + i2, i3)) * (-0.2_wp)
        tt0 = tt0 + (x(i1 + 0, 3 + i2, i3)) * (0.0380952380952381_wp)
        tt0 = tt0 + (x(i1 + 0, 4 + i2, i3)) * (-0.0035714285714285713_wp)
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (4), n - (-4) - (1), 1
        tt0 = 0.0_wp
        do l = -4, min(4, n - (1) - (i2)), 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson8_4_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_fg_201_a_u2_0_false_false_true
SUBROUTINE d_poisson8_fg_201_a_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_fg_201_a_u1_2_false_false_true_cost
SUBROUTINE d_poisson8_fs_10_u1_1_false_true_false(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(-4:n + 4 - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:0) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 = 0, n - (1), 1
      tt(0) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson8_4_fil(l))
      end do
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_fs_10_u1_1_false_true_false
SUBROUTINE d_poisson8_fs_10_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_fs_10_u1_1_false_false_true_cost
SUBROUTINE d_poisson8_fs_10_a_u1_1_false_false_false(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(-4:n + 4 - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
!$omp parallel  default(shared) private(i1, i2, tt0)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 = 0, n - (1), 1
      tt0 = 0.0_wp
      do l = -4, 4, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson8_4_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_fs_10_a_u1_1_false_false_false
SUBROUTINE d_poisson8_fs_10_a_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_fs_10_a_u1_1_false_false_true_cost
SUBROUTINE d_poisson8_fs_01_u2_0_false_false_false(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), -4:n + 4 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
!$omp parallel  default(shared) private(i1, i2, tt0, tt1)
!$omp do 
  do i1 = 0, ndat - (2), 2
    do i2 = 0, n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -4, 4, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson8_4_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson8_4_fil(l))
      end do
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i2 = 0, n - (1), 1
      tt0 = 0.0_wp
      do l = -4, 4, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson8_4_fil(l))
      end do
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_fs_01_u2_0_false_false_false
SUBROUTINE d_poisson8_fs_01_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_fs_01_u1_0_false_false_true_cost
SUBROUTINE d_poisson8_fs_01_a_u4_0_false_true_true(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), -4:n + 4 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:3) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 = 0, n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -4 + i2)) * (0.0035714285714285713_wp)
      tt(1) = tt(1) + (x(i1 + 1, -4 + i2)) * (0.0035714285714285713_wp)
      tt(2) = tt(2) + (x(i1 + 2, -4 + i2)) * (0.0035714285714285713_wp)
      tt(3) = tt(3) + (x(i1 + 3, -4 + i2)) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(i1 + 0, -3 + i2)) * (-0.0380952380952381_wp)
      tt(1) = tt(1) + (x(i1 + 1, -3 + i2)) * (-0.0380952380952381_wp)
      tt(2) = tt(2) + (x(i1 + 2, -3 + i2)) * (-0.0380952380952381_wp)
      tt(3) = tt(3) + (x(i1 + 3, -3 + i2)) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(i1 + 0, -2 + i2)) * (0.2_wp)
      tt(1) = tt(1) + (x(i1 + 1, -2 + i2)) * (0.2_wp)
      tt(2) = tt(2) + (x(i1 + 2, -2 + i2)) * (0.2_wp)
      tt(3) = tt(3) + (x(i1 + 3, -2 + i2)) * (0.2_wp)
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2)) * (-0.8_wp)
      tt(1) = tt(1) + (x(i1 + 1, -1 + i2)) * (-0.8_wp)
      tt(2) = tt(2) + (x(i1 + 2, -1 + i2)) * (-0.8_wp)
      tt(3) = tt(3) + (x(i1 + 3, -1 + i2)) * (-0.8_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt(1) = tt(1) + (x(i1 + 1, 0 + i2)) * (0.0_wp)
      tt(2) = tt(2) + (x(i1 + 2, 0 + i2)) * (0.0_wp)
      tt(3) = tt(3) + (x(i1 + 3, 0 + i2)) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2)) * (0.8_wp)
      tt(1) = tt(1) + (x(i1 + 1, 1 + i2)) * (0.8_wp)
      tt(2) = tt(2) + (x(i1 + 2, 1 + i2)) * (0.8_wp)
      tt(3) = tt(3) + (x(i1 + 3, 1 + i2)) * (0.8_wp)
      tt(0) = tt(0) + (x(i1 + 0, 2 + i2)) * (-0.2_wp)
      tt(1) = tt(1) + (x(i1 + 1, 2 + i2)) * (-0.2_wp)
      tt(2) = tt(2) + (x(i1 + 2, 2 + i2)) * (-0.2_wp)
      tt(3) = tt(3) + (x(i1 + 3, 2 + i2)) * (-0.2_wp)
      tt(0) = tt(0) + (x(i1 + 0, 3 + i2)) * (0.0380952380952381_wp)
      tt(1) = tt(1) + (x(i1 + 1, 3 + i2)) * (0.0380952380952381_wp)
      tt(2) = tt(2) + (x(i1 + 2, 3 + i2)) * (0.0380952380952381_wp)
      tt(3) = tt(3) + (x(i1 + 3, 3 + i2)) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(i1 + 0, 4 + i2)) * (-0.0035714285714285713_wp)
      tt(1) = tt(1) + (x(i1 + 1, 4 + i2)) * (-0.0035714285714285713_wp)
      tt(2) = tt(2) + (x(i1 + 2, 4 + i2)) * (-0.0035714285714285713_wp)
      tt(3) = tt(3) + (x(i1 + 3, 4 + i2)) * (-0.0035714285714285713_wp)
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 = 0, n - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(i1 + 0, -4 + i2)) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(i1 + 0, -3 + i2)) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(i1 + 0, -2 + i2)) * (0.2_wp)
      tt(0) = tt(0) + (x(i1 + 0, -1 + i2)) * (-0.8_wp)
      tt(0) = tt(0) + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt(0) = tt(0) + (x(i1 + 0, 1 + i2)) * (0.8_wp)
      tt(0) = tt(0) + (x(i1 + 0, 2 + i2)) * (-0.2_wp)
      tt(0) = tt(0) + (x(i1 + 0, 3 + i2)) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(i1 + 0, 4 + i2)) * (-0.0035714285714285713_wp)
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_fs_01_a_u4_0_false_true_true
SUBROUTINE d_poisson8_fs_01_a_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_fs_01_a_u1_0_false_false_true_cost
SUBROUTINE d_poisson8_fs_201_u2_0_false_false_false(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), -4:n + 4 - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (2), 2
      do i2 = 0, n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson8_4_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson8_4_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (2)) * (2), ndat1 - (1), 1
      do i2 = 0, n - (1), 1
        tt0 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson8_4_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_fs_201_u2_0_false_false_false
SUBROUTINE d_poisson8_fs_201_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_fs_201_u1_2_false_false_true_cost
SUBROUTINE d_poisson8_fs_201_a_u2_0_false_true_true(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), -4:n + 4 - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
!$omp parallel  default(shared) private(i1, i2, i3, tt)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (2), 2
      do i2 = 0, n - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -4 + i2, i3)) * (0.0035714285714285713_wp)
        tt(1) = tt(1) + (x(i1 + 1, -4 + i2, i3)) * (0.0035714285714285713_wp)
        tt(0) = tt(0) + (x(i1 + 0, -3 + i2, i3)) * (-0.0380952380952381_wp)
        tt(1) = tt(1) + (x(i1 + 1, -3 + i2, i3)) * (-0.0380952380952381_wp)
        tt(0) = tt(0) + (x(i1 + 0, -2 + i2, i3)) * (0.2_wp)
        tt(1) = tt(1) + (x(i1 + 1, -2 + i2, i3)) * (0.2_wp)
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.8_wp)
        tt(1) = tt(1) + (x(i1 + 1, -1 + i2, i3)) * (-0.8_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(1) = tt(1) + (x(i1 + 1, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.8_wp)
        tt(1) = tt(1) + (x(i1 + 1, 1 + i2, i3)) * (0.8_wp)
        tt(0) = tt(0) + (x(i1 + 0, 2 + i2, i3)) * (-0.2_wp)
        tt(1) = tt(1) + (x(i1 + 1, 2 + i2, i3)) * (-0.2_wp)
        tt(0) = tt(0) + (x(i1 + 0, 3 + i2, i3)) * (0.0380952380952381_wp)
        tt(1) = tt(1) + (x(i1 + 1, 3 + i2, i3)) * (0.0380952380952381_wp)
        tt(0) = tt(0) + (x(i1 + 0, 4 + i2, i3)) * (-0.0035714285714285713_wp)
        tt(1) = tt(1) + (x(i1 + 1, 4 + i2, i3)) * (-0.0035714285714285713_wp)
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
        tt(1) = (tt(1)) * (a)
        y(i1 + 1, i2, i3) = tt(1)
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (2)) * (2), ndat1 - (1), 1
      do i2 = 0, n - (1), 1
        tt(0) = 0.0_wp
        tt(0) = tt(0) + (x(i1 + 0, -4 + i2, i3)) * (0.0035714285714285713_wp)
        tt(0) = tt(0) + (x(i1 + 0, -3 + i2, i3)) * (-0.0380952380952381_wp)
        tt(0) = tt(0) + (x(i1 + 0, -2 + i2, i3)) * (0.2_wp)
        tt(0) = tt(0) + (x(i1 + 0, -1 + i2, i3)) * (-0.8_wp)
        tt(0) = tt(0) + (x(i1 + 0, 0 + i2, i3)) * (0.0_wp)
        tt(0) = tt(0) + (x(i1 + 0, 1 + i2, i3)) * (0.8_wp)
        tt(0) = tt(0) + (x(i1 + 0, 2 + i2, i3)) * (-0.2_wp)
        tt(0) = tt(0) + (x(i1 + 0, 3 + i2, i3)) * (0.0380952380952381_wp)
        tt(0) = tt(0) + (x(i1 + 0, 4 + i2, i3)) * (-0.0035714285714285713_wp)
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_fs_201_a_u2_0_false_true_true
SUBROUTINE d_poisson8_fs_201_a_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_fs_201_a_u1_2_false_false_true_cost
SUBROUTINE d_poisson8_np_10_u1_1_true_true_true(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(0:80) :: poisson8_fil = (/ &
-2.717857142857143_wp, &
8.0_wp, &
-14.0_wp, &
18.666666666666668_wp, &
17.5_wp, &
11.2_wp, &
-4.666666666666667_wp, &
1.1428571428571428_wp, &
-0.125_wp, &
-0.125_wp, &
-1.5928571428571427_wp, &
3.5_wp, &
-3.5_wp, &
2.9166666666666665_wp, &
-1.75_wp, &
0.7_wp, &
-0.16666666666666666_wp, &
0.017857142857142856_wp, &
0.017857142857142856_wp, &
-0.2857142857142857_wp, &
-0.95_wp, &
2.0_wp, &
-1.25_wp, &
0.6666666666666666_wp, &
-0.25_wp, &
0.05714285714285714_wp, &
-0.005952380952380952_wp, &
-0.005952380952380952_wp, &
0.07142857142857142_wp, &
-0.5_wp, &
-0.45_wp, &
1.25_wp, &
-0.5_wp, &
0.16666666666666666_wp, &
-0.03571428571428571_wp, &
0.0035714285714285713_wp, &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp, &
-0.0035714285714285713_wp, &
0.03571428571428571_wp, &
-0.16666666666666666_wp, &
0.5_wp, &
-1.25_wp, &
0.45_wp, &
0.5_wp, &
-0.07142857142857142_wp, &
0.005952380952380952_wp, &
0.005952380952380952_wp, &
-0.05714285714285714_wp, &
0.25_wp, &
-0.6666666666666666_wp, &
1.25_wp, &
-2.0_wp, &
0.95_wp, &
0.2857142857142857_wp, &
-0.017857142857142856_wp, &
-0.017857142857142856_wp, &
0.16666666666666666_wp, &
-0.7_wp, &
1.75_wp, &
-2.9166666666666665_wp, &
3.5_wp, &
-3.5_wp, &
1.5928571428571427_wp, &
0.125_wp, &
0.125_wp, &
-1.1428571428571428_wp, &
4.666666666666667_wp, &
-11.2_wp, &
17.5_wp, &
-18.666666666666668_wp, &
14.0_wp, &
-8.0_wp, &
2.717857142857143_wp /)
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:0) :: tt
  integer(kind=4), dimension(-4 - (4):4 - (-4) - (1)) :: mod_arr
  do l = -4 - (4), 4 - (-4) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 = 0,  -(-4) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-4 + 4, i2 + 0)) * (poisson8_fil((i1 - (4)) * (9) + 36 + -4 + 4))
      tt(0) = tt(0) + (x(-3 + 4, i2 + 0)) * (poisson8_fil((i1 - (4)) * (9) + 36 + -3 + 4))
      tt(0) = tt(0) + (x(-2 + 4, i2 + 0)) * (poisson8_fil((i1 - (4)) * (9) + 36 + -2 + 4))
      tt(0) = tt(0) + (x(-1 + 4, i2 + 0)) * (poisson8_fil((i1 - (4)) * (9) + 36 + -1 + 4))
      tt(0) = tt(0) + (x(0 + 4, i2 + 0)) * (poisson8_fil((i1 - (4)) * (9) + 36 + 0 + 4))
      tt(0) = tt(0) + (x(1 + 4, i2 + 0)) * (poisson8_fil((i1 - (4)) * (9) + 36 + 1 + 4))
      tt(0) = tt(0) + (x(2 + 4, i2 + 0)) * (poisson8_fil((i1 - (4)) * (9) + 36 + 2 + 4))
      tt(0) = tt(0) + (x(3 + 4, i2 + 0)) * (poisson8_fil((i1 - (4)) * (9) + 36 + 3 + 4))
      tt(0) = tt(0) + (x(4 + 4, i2 + 0)) * (poisson8_fil((i1 - (4)) * (9) + 36 + 4 + 4))
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-4), n - (4) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-4 + i1, i2 + 0)) * (0.0035714285714285713_wp)
      tt(0) = tt(0) + (x(-3 + i1, i2 + 0)) * (-0.0380952380952381_wp)
      tt(0) = tt(0) + (x(-2 + i1, i2 + 0)) * (0.2_wp)
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.8_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.8_wp)
      tt(0) = tt(0) + (x(2 + i1, i2 + 0)) * (-0.2_wp)
      tt(0) = tt(0) + (x(3 + i1, i2 + 0)) * (0.0380952380952381_wp)
      tt(0) = tt(0) + (x(4 + i1, i2 + 0)) * (-0.0035714285714285713_wp)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (4), n - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-4 + -4 + n - (1), i2 + 0)) * (poisson8_fil((i1 + 4 - (n) + 1) * (9) + 36 + -4 + 4))
      tt(0) = tt(0) + (x(-3 + -4 + n - (1), i2 + 0)) * (poisson8_fil((i1 + 4 - (n) + 1) * (9) + 36 + -3 + 4))
      tt(0) = tt(0) + (x(-2 + -4 + n - (1), i2 + 0)) * (poisson8_fil((i1 + 4 - (n) + 1) * (9) + 36 + -2 + 4))
      tt(0) = tt(0) + (x(-1 + -4 + n - (1), i2 + 0)) * (poisson8_fil((i1 + 4 - (n) + 1) * (9) + 36 + -1 + 4))
      tt(0) = tt(0) + (x(0 + -4 + n - (1), i2 + 0)) * (poisson8_fil((i1 + 4 - (n) + 1) * (9) + 36 + 0 + 4))
      tt(0) = tt(0) + (x(1 + -4 + n - (1), i2 + 0)) * (poisson8_fil((i1 + 4 - (n) + 1) * (9) + 36 + 1 + 4))
      tt(0) = tt(0) + (x(2 + -4 + n - (1), i2 + 0)) * (poisson8_fil((i1 + 4 - (n) + 1) * (9) + 36 + 2 + 4))
      tt(0) = tt(0) + (x(3 + -4 + n - (1), i2 + 0)) * (poisson8_fil((i1 + 4 - (n) + 1) * (9) + 36 + 3 + 4))
      tt(0) = tt(0) + (x(4 + -4 + n - (1), i2 + 0)) * (poisson8_fil((i1 + 4 - (n) + 1) * (9) + 36 + 4 + 4))
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_np_10_u1_1_true_true_true
SUBROUTINE d_poisson8_np_10_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_np_10_u1_1_true_false_true_cost
SUBROUTINE d_poisson8_np_10_a_u2_1_true_true_false(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(0:80) :: poisson8_fil = (/ &
-2.717857142857143_wp, &
8.0_wp, &
-14.0_wp, &
18.666666666666668_wp, &
17.5_wp, &
11.2_wp, &
-4.666666666666667_wp, &
1.1428571428571428_wp, &
-0.125_wp, &
-0.125_wp, &
-1.5928571428571427_wp, &
3.5_wp, &
-3.5_wp, &
2.9166666666666665_wp, &
-1.75_wp, &
0.7_wp, &
-0.16666666666666666_wp, &
0.017857142857142856_wp, &
0.017857142857142856_wp, &
-0.2857142857142857_wp, &
-0.95_wp, &
2.0_wp, &
-1.25_wp, &
0.6666666666666666_wp, &
-0.25_wp, &
0.05714285714285714_wp, &
-0.005952380952380952_wp, &
-0.005952380952380952_wp, &
0.07142857142857142_wp, &
-0.5_wp, &
-0.45_wp, &
1.25_wp, &
-0.5_wp, &
0.16666666666666666_wp, &
-0.03571428571428571_wp, &
0.0035714285714285713_wp, &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp, &
-0.0035714285714285713_wp, &
0.03571428571428571_wp, &
-0.16666666666666666_wp, &
0.5_wp, &
-1.25_wp, &
0.45_wp, &
0.5_wp, &
-0.07142857142857142_wp, &
0.005952380952380952_wp, &
0.005952380952380952_wp, &
-0.05714285714285714_wp, &
0.25_wp, &
-0.6666666666666666_wp, &
1.25_wp, &
-2.0_wp, &
0.95_wp, &
0.2857142857142857_wp, &
-0.017857142857142856_wp, &
-0.017857142857142856_wp, &
0.16666666666666666_wp, &
-0.7_wp, &
1.75_wp, &
-2.9166666666666665_wp, &
3.5_wp, &
-3.5_wp, &
1.5928571428571427_wp, &
0.125_wp, &
0.125_wp, &
-1.1428571428571428_wp, &
4.666666666666667_wp, &
-11.2_wp, &
17.5_wp, &
-18.666666666666668_wp, &
14.0_wp, &
-8.0_wp, &
2.717857142857143_wp /)
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
  integer(kind=4), dimension(-4 - (4):4 - (-4) - (1)) :: mod_arr
  do l = -4 - (4), 4 - (-4) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 = 0,  -(-4) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(l + 4, i2 + 0)) * (poisson8_fil((i1 - (4)) * (9) + 36 + l + 4))
        tt(1) = tt(1) + (x(l + 4, i2 + 1)) * (poisson8_fil((i1 - (4)) * (9) + 36 + l + 4))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
    do i1 =  -(-4), n - (4) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson8_4_fil(l))
        tt(1) = tt(1) + (x(l + i1, i2 + 1)) * (poisson8_4_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
    do i1 = n - (4), n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(l + -4 + n - (1), i2 + 0)) * (poisson8_fil((i1 + 4 - (n) + 1) * (9) + 36 + l + 4))
        tt(1) = tt(1) + (x(l + -4 + n - (1), i2 + 1)) * (poisson8_fil((i1 + 4 - (n) + 1) * (9) + 36 + l + 4))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1, i2 + 1) = tt(1)
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 = 0,  -(-4) - (1), 1
      tt(0) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(l + 4, i2 + 0)) * (poisson8_fil((i1 - (4)) * (9) + 36 + l + 4))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-4), n - (4) - (1), 1
      tt(0) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson8_4_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (4), n - (1), 1
      tt(0) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(l + -4 + n - (1), i2 + 0)) * (poisson8_fil((i1 + 4 - (n) + 1) * (9) + 36 + l + 4))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_np_10_a_u2_1_true_true_false
SUBROUTINE d_poisson8_np_10_a_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_np_10_a_u1_1_true_false_true_cost
SUBROUTINE d_poisson8_np_01_u2_0_false_true_false(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), parameter, dimension(0:80) :: poisson8_fil = (/ &
-2.717857142857143_wp, &
8.0_wp, &
-14.0_wp, &
18.666666666666668_wp, &
17.5_wp, &
11.2_wp, &
-4.666666666666667_wp, &
1.1428571428571428_wp, &
-0.125_wp, &
-0.125_wp, &
-1.5928571428571427_wp, &
3.5_wp, &
-3.5_wp, &
2.9166666666666665_wp, &
-1.75_wp, &
0.7_wp, &
-0.16666666666666666_wp, &
0.017857142857142856_wp, &
0.017857142857142856_wp, &
-0.2857142857142857_wp, &
-0.95_wp, &
2.0_wp, &
-1.25_wp, &
0.6666666666666666_wp, &
-0.25_wp, &
0.05714285714285714_wp, &
-0.005952380952380952_wp, &
-0.005952380952380952_wp, &
0.07142857142857142_wp, &
-0.5_wp, &
-0.45_wp, &
1.25_wp, &
-0.5_wp, &
0.16666666666666666_wp, &
-0.03571428571428571_wp, &
0.0035714285714285713_wp, &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp, &
-0.0035714285714285713_wp, &
0.03571428571428571_wp, &
-0.16666666666666666_wp, &
0.5_wp, &
-1.25_wp, &
0.45_wp, &
0.5_wp, &
-0.07142857142857142_wp, &
0.005952380952380952_wp, &
0.005952380952380952_wp, &
-0.05714285714285714_wp, &
0.25_wp, &
-0.6666666666666666_wp, &
1.25_wp, &
-2.0_wp, &
0.95_wp, &
0.2857142857142857_wp, &
-0.017857142857142856_wp, &
-0.017857142857142856_wp, &
0.16666666666666666_wp, &
-0.7_wp, &
1.75_wp, &
-2.9166666666666665_wp, &
3.5_wp, &
-3.5_wp, &
1.5928571428571427_wp, &
0.125_wp, &
0.125_wp, &
-1.1428571428571428_wp, &
4.666666666666667_wp, &
-11.2_wp, &
17.5_wp, &
-18.666666666666668_wp, &
14.0_wp, &
-8.0_wp, &
2.717857142857143_wp /)
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (2), 2
    do i2 = 0,  -(-4) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(i1 + 0, l + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + l + 4))
        tt(1) = tt(1) + (x(i1 + 1, l + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + l + 4))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
    end do
    do i2 =  -(-4), n - (4) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson8_4_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson8_4_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
    end do
    do i2 = n - (4), n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(i1 + 0, l + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + l + 4))
        tt(1) = tt(1) + (x(i1 + 1, l + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + l + 4))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i2 = 0,  -(-4) - (1), 1
      tt(0) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(i1 + 0, l + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + l + 4))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 =  -(-4), n - (4) - (1), 1
      tt(0) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson8_4_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 = n - (4), n - (1), 1
      tt(0) = 0.0_wp
      do l = -4, 4, 1
        tt(0) = tt(0) + (x(i1 + 0, l + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + l + 4))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_np_01_u2_0_false_true_false
SUBROUTINE d_poisson8_np_01_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_np_01_u1_0_true_false_true_cost
SUBROUTINE d_poisson8_np_01_a_u4_0_false_false_true(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(0:80) :: poisson8_fil = (/ &
-2.717857142857143_wp, &
8.0_wp, &
-14.0_wp, &
18.666666666666668_wp, &
17.5_wp, &
11.2_wp, &
-4.666666666666667_wp, &
1.1428571428571428_wp, &
-0.125_wp, &
-0.125_wp, &
-1.5928571428571427_wp, &
3.5_wp, &
-3.5_wp, &
2.9166666666666665_wp, &
-1.75_wp, &
0.7_wp, &
-0.16666666666666666_wp, &
0.017857142857142856_wp, &
0.017857142857142856_wp, &
-0.2857142857142857_wp, &
-0.95_wp, &
2.0_wp, &
-1.25_wp, &
0.6666666666666666_wp, &
-0.25_wp, &
0.05714285714285714_wp, &
-0.005952380952380952_wp, &
-0.005952380952380952_wp, &
0.07142857142857142_wp, &
-0.5_wp, &
-0.45_wp, &
1.25_wp, &
-0.5_wp, &
0.16666666666666666_wp, &
-0.03571428571428571_wp, &
0.0035714285714285713_wp, &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp, &
-0.0035714285714285713_wp, &
0.03571428571428571_wp, &
-0.16666666666666666_wp, &
0.5_wp, &
-1.25_wp, &
0.45_wp, &
0.5_wp, &
-0.07142857142857142_wp, &
0.005952380952380952_wp, &
0.005952380952380952_wp, &
-0.05714285714285714_wp, &
0.25_wp, &
-0.6666666666666666_wp, &
1.25_wp, &
-2.0_wp, &
0.95_wp, &
0.2857142857142857_wp, &
-0.017857142857142856_wp, &
-0.017857142857142856_wp, &
0.16666666666666666_wp, &
-0.7_wp, &
1.75_wp, &
-2.9166666666666665_wp, &
3.5_wp, &
-3.5_wp, &
1.5928571428571427_wp, &
0.125_wp, &
0.125_wp, &
-1.1428571428571428_wp, &
4.666666666666667_wp, &
-11.2_wp, &
17.5_wp, &
-18.666666666666668_wp, &
14.0_wp, &
-8.0_wp, &
2.717857142857143_wp /)
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2, tt3)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 = 0,  -(-4) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -4 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -4 + 4))
      tt1 = tt1 + (x(i1 + 1, -4 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -4 + 4))
      tt2 = tt2 + (x(i1 + 2, -4 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -4 + 4))
      tt3 = tt3 + (x(i1 + 3, -4 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -4 + 4))
      tt0 = tt0 + (x(i1 + 0, -3 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -3 + 4))
      tt1 = tt1 + (x(i1 + 1, -3 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -3 + 4))
      tt2 = tt2 + (x(i1 + 2, -3 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -3 + 4))
      tt3 = tt3 + (x(i1 + 3, -3 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -3 + 4))
      tt0 = tt0 + (x(i1 + 0, -2 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -2 + 4))
      tt1 = tt1 + (x(i1 + 1, -2 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -2 + 4))
      tt2 = tt2 + (x(i1 + 2, -2 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -2 + 4))
      tt3 = tt3 + (x(i1 + 3, -2 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -2 + 4))
      tt0 = tt0 + (x(i1 + 0, -1 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -1 + 4))
      tt1 = tt1 + (x(i1 + 1, -1 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -1 + 4))
      tt2 = tt2 + (x(i1 + 2, -1 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -1 + 4))
      tt3 = tt3 + (x(i1 + 3, -1 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -1 + 4))
      tt0 = tt0 + (x(i1 + 0, 0 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 0 + 4))
      tt1 = tt1 + (x(i1 + 1, 0 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 0 + 4))
      tt2 = tt2 + (x(i1 + 2, 0 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 0 + 4))
      tt3 = tt3 + (x(i1 + 3, 0 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 0 + 4))
      tt0 = tt0 + (x(i1 + 0, 1 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 1 + 4))
      tt1 = tt1 + (x(i1 + 1, 1 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 1 + 4))
      tt2 = tt2 + (x(i1 + 2, 1 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 1 + 4))
      tt3 = tt3 + (x(i1 + 3, 1 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 1 + 4))
      tt0 = tt0 + (x(i1 + 0, 2 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 2 + 4))
      tt1 = tt1 + (x(i1 + 1, 2 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 2 + 4))
      tt2 = tt2 + (x(i1 + 2, 2 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 2 + 4))
      tt3 = tt3 + (x(i1 + 3, 2 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 2 + 4))
      tt0 = tt0 + (x(i1 + 0, 3 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 3 + 4))
      tt1 = tt1 + (x(i1 + 1, 3 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 3 + 4))
      tt2 = tt2 + (x(i1 + 2, 3 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 3 + 4))
      tt3 = tt3 + (x(i1 + 3, 3 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 3 + 4))
      tt0 = tt0 + (x(i1 + 0, 4 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 4 + 4))
      tt1 = tt1 + (x(i1 + 1, 4 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 4 + 4))
      tt2 = tt2 + (x(i1 + 2, 4 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 4 + 4))
      tt3 = tt3 + (x(i1 + 3, 4 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 4 + 4))
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
    do i2 =  -(-4), n - (4) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -4 + i2)) * (0.0035714285714285713_wp)
      tt1 = tt1 + (x(i1 + 1, -4 + i2)) * (0.0035714285714285713_wp)
      tt2 = tt2 + (x(i1 + 2, -4 + i2)) * (0.0035714285714285713_wp)
      tt3 = tt3 + (x(i1 + 3, -4 + i2)) * (0.0035714285714285713_wp)
      tt0 = tt0 + (x(i1 + 0, -3 + i2)) * (-0.0380952380952381_wp)
      tt1 = tt1 + (x(i1 + 1, -3 + i2)) * (-0.0380952380952381_wp)
      tt2 = tt2 + (x(i1 + 2, -3 + i2)) * (-0.0380952380952381_wp)
      tt3 = tt3 + (x(i1 + 3, -3 + i2)) * (-0.0380952380952381_wp)
      tt0 = tt0 + (x(i1 + 0, -2 + i2)) * (0.2_wp)
      tt1 = tt1 + (x(i1 + 1, -2 + i2)) * (0.2_wp)
      tt2 = tt2 + (x(i1 + 2, -2 + i2)) * (0.2_wp)
      tt3 = tt3 + (x(i1 + 3, -2 + i2)) * (0.2_wp)
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.8_wp)
      tt1 = tt1 + (x(i1 + 1, -1 + i2)) * (-0.8_wp)
      tt2 = tt2 + (x(i1 + 2, -1 + i2)) * (-0.8_wp)
      tt3 = tt3 + (x(i1 + 3, -1 + i2)) * (-0.8_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt1 = tt1 + (x(i1 + 1, 0 + i2)) * (0.0_wp)
      tt2 = tt2 + (x(i1 + 2, 0 + i2)) * (0.0_wp)
      tt3 = tt3 + (x(i1 + 3, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.8_wp)
      tt1 = tt1 + (x(i1 + 1, 1 + i2)) * (0.8_wp)
      tt2 = tt2 + (x(i1 + 2, 1 + i2)) * (0.8_wp)
      tt3 = tt3 + (x(i1 + 3, 1 + i2)) * (0.8_wp)
      tt0 = tt0 + (x(i1 + 0, 2 + i2)) * (-0.2_wp)
      tt1 = tt1 + (x(i1 + 1, 2 + i2)) * (-0.2_wp)
      tt2 = tt2 + (x(i1 + 2, 2 + i2)) * (-0.2_wp)
      tt3 = tt3 + (x(i1 + 3, 2 + i2)) * (-0.2_wp)
      tt0 = tt0 + (x(i1 + 0, 3 + i2)) * (0.0380952380952381_wp)
      tt1 = tt1 + (x(i1 + 1, 3 + i2)) * (0.0380952380952381_wp)
      tt2 = tt2 + (x(i1 + 2, 3 + i2)) * (0.0380952380952381_wp)
      tt3 = tt3 + (x(i1 + 3, 3 + i2)) * (0.0380952380952381_wp)
      tt0 = tt0 + (x(i1 + 0, 4 + i2)) * (-0.0035714285714285713_wp)
      tt1 = tt1 + (x(i1 + 1, 4 + i2)) * (-0.0035714285714285713_wp)
      tt2 = tt2 + (x(i1 + 2, 4 + i2)) * (-0.0035714285714285713_wp)
      tt3 = tt3 + (x(i1 + 3, 4 + i2)) * (-0.0035714285714285713_wp)
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
    do i2 = n - (4), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -4 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -4 + 4))
      tt1 = tt1 + (x(i1 + 1, -4 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -4 + 4))
      tt2 = tt2 + (x(i1 + 2, -4 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -4 + 4))
      tt3 = tt3 + (x(i1 + 3, -4 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -4 + 4))
      tt0 = tt0 + (x(i1 + 0, -3 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -3 + 4))
      tt1 = tt1 + (x(i1 + 1, -3 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -3 + 4))
      tt2 = tt2 + (x(i1 + 2, -3 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -3 + 4))
      tt3 = tt3 + (x(i1 + 3, -3 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -3 + 4))
      tt0 = tt0 + (x(i1 + 0, -2 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -2 + 4))
      tt1 = tt1 + (x(i1 + 1, -2 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -2 + 4))
      tt2 = tt2 + (x(i1 + 2, -2 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -2 + 4))
      tt3 = tt3 + (x(i1 + 3, -2 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -2 + 4))
      tt0 = tt0 + (x(i1 + 0, -1 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -1 + 4))
      tt1 = tt1 + (x(i1 + 1, -1 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -1 + 4))
      tt2 = tt2 + (x(i1 + 2, -1 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -1 + 4))
      tt3 = tt3 + (x(i1 + 3, -1 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -1 + 4))
      tt0 = tt0 + (x(i1 + 0, 0 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 0 + 4))
      tt1 = tt1 + (x(i1 + 1, 0 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 0 + 4))
      tt2 = tt2 + (x(i1 + 2, 0 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 0 + 4))
      tt3 = tt3 + (x(i1 + 3, 0 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 0 + 4))
      tt0 = tt0 + (x(i1 + 0, 1 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 1 + 4))
      tt1 = tt1 + (x(i1 + 1, 1 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 1 + 4))
      tt2 = tt2 + (x(i1 + 2, 1 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 1 + 4))
      tt3 = tt3 + (x(i1 + 3, 1 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 1 + 4))
      tt0 = tt0 + (x(i1 + 0, 2 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 2 + 4))
      tt1 = tt1 + (x(i1 + 1, 2 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 2 + 4))
      tt2 = tt2 + (x(i1 + 2, 2 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 2 + 4))
      tt3 = tt3 + (x(i1 + 3, 2 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 2 + 4))
      tt0 = tt0 + (x(i1 + 0, 3 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 3 + 4))
      tt1 = tt1 + (x(i1 + 1, 3 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 3 + 4))
      tt2 = tt2 + (x(i1 + 2, 3 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 3 + 4))
      tt3 = tt3 + (x(i1 + 3, 3 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 3 + 4))
      tt0 = tt0 + (x(i1 + 0, 4 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 4 + 4))
      tt1 = tt1 + (x(i1 + 1, 4 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 4 + 4))
      tt2 = tt2 + (x(i1 + 2, 4 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 4 + 4))
      tt3 = tt3 + (x(i1 + 3, 4 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 4 + 4))
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 = 0,  -(-4) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -4 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -4 + 4))
      tt0 = tt0 + (x(i1 + 0, -3 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -3 + 4))
      tt0 = tt0 + (x(i1 + 0, -2 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -2 + 4))
      tt0 = tt0 + (x(i1 + 0, -1 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + -1 + 4))
      tt0 = tt0 + (x(i1 + 0, 0 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 0 + 4))
      tt0 = tt0 + (x(i1 + 0, 1 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 1 + 4))
      tt0 = tt0 + (x(i1 + 0, 2 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 2 + 4))
      tt0 = tt0 + (x(i1 + 0, 3 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 3 + 4))
      tt0 = tt0 + (x(i1 + 0, 4 + 4)) * (poisson8_fil((i2 - (4)) * (9) + 36 + 4 + 4))
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
    do i2 =  -(-4), n - (4) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -4 + i2)) * (0.0035714285714285713_wp)
      tt0 = tt0 + (x(i1 + 0, -3 + i2)) * (-0.0380952380952381_wp)
      tt0 = tt0 + (x(i1 + 0, -2 + i2)) * (0.2_wp)
      tt0 = tt0 + (x(i1 + 0, -1 + i2)) * (-0.8_wp)
      tt0 = tt0 + (x(i1 + 0, 0 + i2)) * (0.0_wp)
      tt0 = tt0 + (x(i1 + 0, 1 + i2)) * (0.8_wp)
      tt0 = tt0 + (x(i1 + 0, 2 + i2)) * (-0.2_wp)
      tt0 = tt0 + (x(i1 + 0, 3 + i2)) * (0.0380952380952381_wp)
      tt0 = tt0 + (x(i1 + 0, 4 + i2)) * (-0.0035714285714285713_wp)
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
    do i2 = n - (4), n - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(i1 + 0, -4 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -4 + 4))
      tt0 = tt0 + (x(i1 + 0, -3 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -3 + 4))
      tt0 = tt0 + (x(i1 + 0, -2 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -2 + 4))
      tt0 = tt0 + (x(i1 + 0, -1 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + -1 + 4))
      tt0 = tt0 + (x(i1 + 0, 0 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 0 + 4))
      tt0 = tt0 + (x(i1 + 0, 1 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 1 + 4))
      tt0 = tt0 + (x(i1 + 0, 2 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 2 + 4))
      tt0 = tt0 + (x(i1 + 0, 3 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 3 + 4))
      tt0 = tt0 + (x(i1 + 0, 4 + -4 + n - (1))) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + 4 + 4))
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_np_01_a_u4_0_false_false_true
SUBROUTINE d_poisson8_np_01_a_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_np_01_a_u1_0_true_false_true_cost
SUBROUTINE d_poisson8_np_201_u2_0_true_false_false(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(0:80) :: poisson8_fil = (/ &
-2.717857142857143_wp, &
8.0_wp, &
-14.0_wp, &
18.666666666666668_wp, &
17.5_wp, &
11.2_wp, &
-4.666666666666667_wp, &
1.1428571428571428_wp, &
-0.125_wp, &
-0.125_wp, &
-1.5928571428571427_wp, &
3.5_wp, &
-3.5_wp, &
2.9166666666666665_wp, &
-1.75_wp, &
0.7_wp, &
-0.16666666666666666_wp, &
0.017857142857142856_wp, &
0.017857142857142856_wp, &
-0.2857142857142857_wp, &
-0.95_wp, &
2.0_wp, &
-1.25_wp, &
0.6666666666666666_wp, &
-0.25_wp, &
0.05714285714285714_wp, &
-0.005952380952380952_wp, &
-0.005952380952380952_wp, &
0.07142857142857142_wp, &
-0.5_wp, &
-0.45_wp, &
1.25_wp, &
-0.5_wp, &
0.16666666666666666_wp, &
-0.03571428571428571_wp, &
0.0035714285714285713_wp, &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp, &
-0.0035714285714285713_wp, &
0.03571428571428571_wp, &
-0.16666666666666666_wp, &
0.5_wp, &
-1.25_wp, &
0.45_wp, &
0.5_wp, &
-0.07142857142857142_wp, &
0.005952380952380952_wp, &
0.005952380952380952_wp, &
-0.05714285714285714_wp, &
0.25_wp, &
-0.6666666666666666_wp, &
1.25_wp, &
-2.0_wp, &
0.95_wp, &
0.2857142857142857_wp, &
-0.017857142857142856_wp, &
-0.017857142857142856_wp, &
0.16666666666666666_wp, &
-0.7_wp, &
1.75_wp, &
-2.9166666666666665_wp, &
3.5_wp, &
-3.5_wp, &
1.5928571428571427_wp, &
0.125_wp, &
0.125_wp, &
-1.1428571428571428_wp, &
4.666666666666667_wp, &
-11.2_wp, &
17.5_wp, &
-18.666666666666668_wp, &
14.0_wp, &
-8.0_wp, &
2.717857142857143_wp /)
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  integer(kind=4), dimension(-4 - (4):4 - (-4) - (1)) :: mod_arr
  do l = -4 - (4), 4 - (-4) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (2), 2
      do i2 = 0,  -(-4) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, l + 4, i3)) * (poisson8_fil((i2 - (4)) * (9) + 36 + l + 4))
          tt1 = tt1 + (x(i1 + 1, l + 4, i3)) * (poisson8_fil((i2 - (4)) * (9) + 36 + l + 4))
        end do
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
      end do
      do i2 =  -(-4), n - (4) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson8_4_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson8_4_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
      end do
      do i2 = n - (4), n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, l + -4 + n - (1), i3)) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + l + 4))
          tt1 = tt1 + (x(i1 + 1, l + -4 + n - (1), i3)) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + l + 4))
        end do
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (2)) * (2), ndat1 - (1), 1
      do i2 = 0,  -(-4) - (1), 1
        tt0 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, l + 4, i3)) * (poisson8_fil((i2 - (4)) * (9) + 36 + l + 4))
        end do
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-4), n - (4) - (1), 1
        tt0 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson8_4_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (4), n - (1), 1
        tt0 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, l + -4 + n - (1), i3)) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + l + 4))
        end do
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_np_201_u2_0_true_false_false
SUBROUTINE d_poisson8_np_201_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_np_201_u1_2_true_false_true_cost
SUBROUTINE d_poisson8_np_201_a_u2_0_true_false_false(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -4
  integer(kind=4), parameter :: upfil = 4
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(0:80) :: poisson8_fil = (/ &
-2.717857142857143_wp, &
8.0_wp, &
-14.0_wp, &
18.666666666666668_wp, &
17.5_wp, &
11.2_wp, &
-4.666666666666667_wp, &
1.1428571428571428_wp, &
-0.125_wp, &
-0.125_wp, &
-1.5928571428571427_wp, &
3.5_wp, &
-3.5_wp, &
2.9166666666666665_wp, &
-1.75_wp, &
0.7_wp, &
-0.16666666666666666_wp, &
0.017857142857142856_wp, &
0.017857142857142856_wp, &
-0.2857142857142857_wp, &
-0.95_wp, &
2.0_wp, &
-1.25_wp, &
0.6666666666666666_wp, &
-0.25_wp, &
0.05714285714285714_wp, &
-0.005952380952380952_wp, &
-0.005952380952380952_wp, &
0.07142857142857142_wp, &
-0.5_wp, &
-0.45_wp, &
1.25_wp, &
-0.5_wp, &
0.16666666666666666_wp, &
-0.03571428571428571_wp, &
0.0035714285714285713_wp, &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp, &
-0.0035714285714285713_wp, &
0.03571428571428571_wp, &
-0.16666666666666666_wp, &
0.5_wp, &
-1.25_wp, &
0.45_wp, &
0.5_wp, &
-0.07142857142857142_wp, &
0.005952380952380952_wp, &
0.005952380952380952_wp, &
-0.05714285714285714_wp, &
0.25_wp, &
-0.6666666666666666_wp, &
1.25_wp, &
-2.0_wp, &
0.95_wp, &
0.2857142857142857_wp, &
-0.017857142857142856_wp, &
-0.017857142857142856_wp, &
0.16666666666666666_wp, &
-0.7_wp, &
1.75_wp, &
-2.9166666666666665_wp, &
3.5_wp, &
-3.5_wp, &
1.5928571428571427_wp, &
0.125_wp, &
0.125_wp, &
-1.1428571428571428_wp, &
4.666666666666667_wp, &
-11.2_wp, &
17.5_wp, &
-18.666666666666668_wp, &
14.0_wp, &
-8.0_wp, &
2.717857142857143_wp /)
  real(kind=8), parameter, dimension(-4:4) :: poisson8_4_fil = (/ &
0.0035714285714285713_wp, &
-0.0380952380952381_wp, &
0.2_wp, &
-0.8_wp, &
0.0_wp, &
0.8_wp, &
-0.2_wp, &
0.0380952380952381_wp, &
-0.0035714285714285713_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  integer(kind=4), dimension(-4 - (4):4 - (-4) - (1)) :: mod_arr
  do l = -4 - (4), 4 - (-4) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (2), 2
      do i2 = 0,  -(-4) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, l + 4, i3)) * (poisson8_fil((i2 - (4)) * (9) + 36 + l + 4))
          tt1 = tt1 + (x(i1 + 1, l + 4, i3)) * (poisson8_fil((i2 - (4)) * (9) + 36 + l + 4))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
      end do
      do i2 =  -(-4), n - (4) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson8_4_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson8_4_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
      end do
      do i2 = n - (4), n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, l + -4 + n - (1), i3)) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + l + 4))
          tt1 = tt1 + (x(i1 + 1, l + -4 + n - (1), i3)) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + l + 4))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (2)) * (2), ndat1 - (1), 1
      do i2 = 0,  -(-4) - (1), 1
        tt0 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, l + 4, i3)) * (poisson8_fil((i2 - (4)) * (9) + 36 + l + 4))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-4), n - (4) - (1), 1
        tt0 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson8_4_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (4), n - (1), 1
        tt0 = 0.0_wp
        do l = -4, 4, 1
          tt0 = tt0 + (x(i1 + 0, l + -4 + n - (1), i3)) * (poisson8_fil((i2 + 4 - (n) + 1) * (9) + 36 + l + 4))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson8_np_201_a_u2_0_true_false_false
SUBROUTINE d_poisson8_np_201_a_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (9)) * (ndat_t)
END SUBROUTINE d_poisson8_np_201_a_u1_2_true_false_true_cost
SUBROUTINE d_s0s0_1d_poisson8_cost(d, idim, n, bc, x, y, a, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: d
  integer(kind=4), intent(in) :: idim
  integer(kind=4), intent(in), dimension(0:d - (1)) :: n
  integer(kind=4), intent(in) :: bc
  real(kind=8), intent(in), dimension(*) :: x
  real(kind=8), intent(out), dimension(*) :: y
  real(kind=8), intent(in) :: a
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: i
  integer(kind=4) :: ndat_left
  integer(kind=4) :: ndat_right
  integer(kind=4) :: c
  if (idim == 0) then
    ndat_right = 1
    do i = 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson8_p_10_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson8_p_10_a_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson8_fg_10_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson8_fg_10_a_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson8_fs_10_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson8_fs_10_a_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-2)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson8_np_10_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson8_np_10_a_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      end select
  else if (idim == d - (1)) then
    ndat_left = 1
    do i = 0, d - (2), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson8_p_01_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson8_p_01_a_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      case (1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson8_fg_01_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson8_fg_01_a_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      case (-1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson8_fs_01_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson8_fs_01_a_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      case (-2)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson8_np_01_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson8_np_01_a_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      end select
  else
    ndat_left = 1
    ndat_right = 1
    do i = 0, idim - (1), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    do i = idim + 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson8_p_201_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson8_p_201_a_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson8_fg_201_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson8_fg_201_a_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson8_fs_201_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson8_fs_201_a_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-2)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson8_np_201_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson8_np_201_a_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      end select
  end if
END SUBROUTINE d_s0s0_1d_poisson8_cost
SUBROUTINE d_s0s0_1d_poisson8(d, idim, n, bc, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: d
  integer(kind=4), intent(in) :: idim
  integer(kind=4), intent(in), dimension(0:d - (1)) :: n
  integer(kind=4), intent(in) :: bc
  real(kind=8), intent(in), dimension(*) :: x
  real(kind=8), intent(out), dimension(*) :: y
  real(kind=8), intent(in) :: a
  integer(kind=4) :: i
  integer(kind=4) :: ndat_left
  integer(kind=4) :: ndat_right
  if (idim == 0) then
    ndat_right = 1
    do i = 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          call d_poisson8_p_10_u3_1_true_true_true(n(idim), ndat_right, x, y)
        else
          call d_poisson8_p_10_a_u1_1_true_true_true(n(idim), ndat_right, x, y, a)
        end if
      case (1)
        if (a == 1.0_wp) then
          call d_poisson8_fg_10_u1_1_false_true_true(n(idim), ndat_right, x, y)
        else
          call d_poisson8_fg_10_a_u2_1_false_true_false(n(idim), ndat_right, x, y, a)
        end if
      case (-1)
        if (a == 1.0_wp) then
          call d_poisson8_fs_10_u1_1_false_true_false(n(idim), ndat_right, x, y)
        else
          call d_poisson8_fs_10_a_u1_1_false_false_false(n(idim), ndat_right, x, y, a)
        end if
      case (-2)
        if (a == 1.0_wp) then
          call d_poisson8_np_10_u1_1_true_true_true(n(idim), ndat_right, x, y)
        else
          call d_poisson8_np_10_a_u2_1_true_true_false(n(idim), ndat_right, x, y, a)
        end if
      end select
  else if (idim == d - (1)) then
    ndat_left = 1
    do i = 0, d - (2), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          call d_poisson8_p_01_u2_0_true_true_true(ndat_left, n(idim), x, y)
        else
          call d_poisson8_p_01_a_u2_0_false_true_false(ndat_left, n(idim), x, y, a)
        end if
      case (1)
        if (a == 1.0_wp) then
          call d_poisson8_fg_01_u2_0_false_false_true(ndat_left, n(idim), x, y)
        else
          call d_poisson8_fg_01_a_u3_0_false_false_false(ndat_left, n(idim), x, y, a)
        end if
      case (-1)
        if (a == 1.0_wp) then
          call d_poisson8_fs_01_u2_0_false_false_false(ndat_left, n(idim), x, y)
        else
          call d_poisson8_fs_01_a_u4_0_false_true_true(ndat_left, n(idim), x, y, a)
        end if
      case (-2)
        if (a == 1.0_wp) then
          call d_poisson8_np_01_u2_0_false_true_false(ndat_left, n(idim), x, y)
        else
          call d_poisson8_np_01_a_u4_0_false_false_true(ndat_left, n(idim), x, y, a)
        end if
      end select
  else
    ndat_left = 1
    ndat_right = 1
    do i = 0, idim - (1), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    do i = idim + 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          call d_poisson8_p_201_u2_0_true_false_true(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson8_p_201_a_u2_0_true_false_false(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      case (1)
        if (a == 1.0_wp) then
          call d_poisson8_fg_201_u2_0_false_true_true(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson8_fg_201_a_u2_0_false_false_true(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      case (-1)
        if (a == 1.0_wp) then
          call d_poisson8_fs_201_u2_0_false_false_false(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson8_fs_201_a_u2_0_false_true_true(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      case (-2)
        if (a == 1.0_wp) then
          call d_poisson8_np_201_u2_0_true_false_false(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson8_np_201_a_u2_0_true_false_false(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      end select
  end if
END SUBROUTINE d_s0s0_1d_poisson8
SUBROUTINE d_poisson16_p_10_u2_1_true_true_false(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:1) :: tt
  integer(kind=4), dimension(-8 - (8):8 - (-8) - (1)) :: mod_arr
  do l = -8 - (8), 8 - (-8) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 = 0,  -(-8) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(mod_arr(l + i1), i2 + 0)) * (poisson16_8_fil(l))
        tt(1) = tt(1) + (x(mod_arr(l + i1), i2 + 1)) * (poisson16_8_fil(l))
      end do
      y(i1, i2 + 0) = tt(0)
      y(i1, i2 + 1) = tt(1)
    end do
    do i1 =  -(-8), n - (8) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson16_8_fil(l))
        tt(1) = tt(1) + (x(l + i1, i2 + 1)) * (poisson16_8_fil(l))
      end do
      y(i1, i2 + 0) = tt(0)
      y(i1, i2 + 1) = tt(1)
    end do
    do i1 = n - (8), n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(mod_arr(l + i1 - (n)), i2 + 0)) * (poisson16_8_fil(l))
        tt(1) = tt(1) + (x(mod_arr(l + i1 - (n)), i2 + 1)) * (poisson16_8_fil(l))
      end do
      y(i1, i2 + 0) = tt(0)
      y(i1, i2 + 1) = tt(1)
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 = 0,  -(-8) - (1), 1
      tt(0) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(mod_arr(l + i1), i2 + 0)) * (poisson16_8_fil(l))
      end do
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-8), n - (8) - (1), 1
      tt(0) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson16_8_fil(l))
      end do
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (8), n - (1), 1
      tt(0) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(mod_arr(l + i1 - (n)), i2 + 0)) * (poisson16_8_fil(l))
      end do
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_p_10_u2_1_true_true_false
SUBROUTINE d_poisson16_p_10_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_p_10_u1_1_true_false_true_cost
SUBROUTINE d_poisson16_p_10_a_u2_1_true_false_false(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  integer(kind=4), dimension(-8 - (8):8 - (-8) - (1)) :: mod_arr
  do l = -8 - (8), 8 - (-8) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt0, tt1)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 = 0,  -(-8) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(mod_arr(l + i1), i2 + 0)) * (poisson16_8_fil(l))
        tt1 = tt1 + (x(mod_arr(l + i1), i2 + 1)) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
      tt1 = (tt1) * (a)
      y(i1, i2 + 1) = tt1
    end do
    do i1 =  -(-8), n - (8) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson16_8_fil(l))
        tt1 = tt1 + (x(l + i1, i2 + 1)) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
      tt1 = (tt1) * (a)
      y(i1, i2 + 1) = tt1
    end do
    do i1 = n - (8), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(mod_arr(l + i1 - (n)), i2 + 0)) * (poisson16_8_fil(l))
        tt1 = tt1 + (x(mod_arr(l + i1 - (n)), i2 + 1)) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
      tt1 = (tt1) * (a)
      y(i1, i2 + 1) = tt1
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 = 0,  -(-8) - (1), 1
      tt0 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(mod_arr(l + i1), i2 + 0)) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
    end do
    do i1 =  -(-8), n - (8) - (1), 1
      tt0 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
    end do
    do i1 = n - (8), n - (1), 1
      tt0 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(mod_arr(l + i1 - (n)), i2 + 0)) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_p_10_a_u2_1_true_false_false
SUBROUTINE d_poisson16_p_10_a_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_p_10_a_u1_1_true_false_true_cost
SUBROUTINE d_poisson16_p_01_u4_0_true_false_false(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  integer(kind=4), dimension(-8 - (8):8 - (-8) - (1)) :: mod_arr
  do l = -8 - (8), 8 - (-8) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2, tt3)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 = 0,  -(-8) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2))) * (poisson16_8_fil(l))
        tt1 = tt1 + (x(i1 + 1, mod_arr(l + i2))) * (poisson16_8_fil(l))
        tt2 = tt2 + (x(i1 + 2, mod_arr(l + i2))) * (poisson16_8_fil(l))
        tt3 = tt3 + (x(i1 + 3, mod_arr(l + i2))) * (poisson16_8_fil(l))
      end do
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
      y(i1 + 3, i2) = tt3
    end do
    do i2 =  -(-8), n - (8) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson16_8_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2)) * (poisson16_8_fil(l))
        tt3 = tt3 + (x(i1 + 3, l + i2)) * (poisson16_8_fil(l))
      end do
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
      y(i1 + 3, i2) = tt3
    end do
    do i2 = n - (8), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2 - (n)))) * (poisson16_8_fil(l))
        tt1 = tt1 + (x(i1 + 1, mod_arr(l + i2 - (n)))) * (poisson16_8_fil(l))
        tt2 = tt2 + (x(i1 + 2, mod_arr(l + i2 - (n)))) * (poisson16_8_fil(l))
        tt3 = tt3 + (x(i1 + 3, mod_arr(l + i2 - (n)))) * (poisson16_8_fil(l))
      end do
      y(i1 + 0, i2) = tt0
      y(i1 + 1, i2) = tt1
      y(i1 + 2, i2) = tt2
      y(i1 + 3, i2) = tt3
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 = 0,  -(-8) - (1), 1
      tt0 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2))) * (poisson16_8_fil(l))
      end do
      y(i1 + 0, i2) = tt0
    end do
    do i2 =  -(-8), n - (8) - (1), 1
      tt0 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
      end do
      y(i1 + 0, i2) = tt0
    end do
    do i2 = n - (8), n - (1), 1
      tt0 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2 - (n)))) * (poisson16_8_fil(l))
      end do
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_p_01_u4_0_true_false_false
SUBROUTINE d_poisson16_p_01_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_p_01_u1_0_true_false_true_cost
SUBROUTINE d_poisson16_p_01_a_u4_0_true_false_false(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  integer(kind=4), dimension(-8 - (8):8 - (-8) - (1)) :: mod_arr
  do l = -8 - (8), 8 - (-8) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2, tt3)
!$omp do 
  do i1 = 0, ndat - (4), 4
    do i2 = 0,  -(-8) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2))) * (poisson16_8_fil(l))
        tt1 = tt1 + (x(i1 + 1, mod_arr(l + i2))) * (poisson16_8_fil(l))
        tt2 = tt2 + (x(i1 + 2, mod_arr(l + i2))) * (poisson16_8_fil(l))
        tt3 = tt3 + (x(i1 + 3, mod_arr(l + i2))) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
    do i2 =  -(-8), n - (8) - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson16_8_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2)) * (poisson16_8_fil(l))
        tt3 = tt3 + (x(i1 + 3, l + i2)) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
    do i2 = n - (8), n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2 - (n)))) * (poisson16_8_fil(l))
        tt1 = tt1 + (x(i1 + 1, mod_arr(l + i2 - (n)))) * (poisson16_8_fil(l))
        tt2 = tt2 + (x(i1 + 2, mod_arr(l + i2 - (n)))) * (poisson16_8_fil(l))
        tt3 = tt3 + (x(i1 + 3, mod_arr(l + i2 - (n)))) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (4)) * (4), ndat - (1), 1
    do i2 = 0,  -(-8) - (1), 1
      tt0 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2))) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
    do i2 =  -(-8), n - (8) - (1), 1
      tt0 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
    do i2 = n - (8), n - (1), 1
      tt0 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2 - (n)))) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_p_01_a_u4_0_true_false_false
SUBROUTINE d_poisson16_p_01_a_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_p_01_a_u1_0_true_false_true_cost
SUBROUTINE d_poisson16_p_201_u4_0_true_true_false(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8), dimension(0:3) :: tt
  integer(kind=4), dimension(-8 - (8):8 - (-8) - (1)) :: mod_arr
  do l = -8 - (8), 8 - (-8) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 = 0,  -(-8) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        do l = -8, 8, 1
          tt(0) = tt(0) + (x(i1 + 0, mod_arr(l + i2), i3)) * (poisson16_8_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, mod_arr(l + i2), i3)) * (poisson16_8_fil(l))
          tt(2) = tt(2) + (x(i1 + 2, mod_arr(l + i2), i3)) * (poisson16_8_fil(l))
          tt(3) = tt(3) + (x(i1 + 3, mod_arr(l + i2), i3)) * (poisson16_8_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
      end do
      do i2 =  -(-8), n - (8) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        do l = -8, 8, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, l + i2, i3)) * (poisson16_8_fil(l))
          tt(2) = tt(2) + (x(i1 + 2, l + i2, i3)) * (poisson16_8_fil(l))
          tt(3) = tt(3) + (x(i1 + 3, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
      end do
      do i2 = n - (8), n - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        do l = -8, 8, 1
          tt(0) = tt(0) + (x(i1 + 0, mod_arr(l + i2 - (n)), i3)) * (poisson16_8_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, mod_arr(l + i2 - (n)), i3)) * (poisson16_8_fil(l))
          tt(2) = tt(2) + (x(i1 + 2, mod_arr(l + i2 - (n)), i3)) * (poisson16_8_fil(l))
          tt(3) = tt(3) + (x(i1 + 3, mod_arr(l + i2 - (n)), i3)) * (poisson16_8_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 = 0,  -(-8) - (1), 1
        tt(0) = 0.0_wp
        do l = -8, 8, 1
          tt(0) = tt(0) + (x(i1 + 0, mod_arr(l + i2), i3)) * (poisson16_8_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 =  -(-8), n - (8) - (1), 1
        tt(0) = 0.0_wp
        do l = -8, 8, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 = n - (8), n - (1), 1
        tt(0) = 0.0_wp
        do l = -8, 8, 1
          tt(0) = tt(0) + (x(i1 + 0, mod_arr(l + i2 - (n)), i3)) * (poisson16_8_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_p_201_u4_0_true_true_false
SUBROUTINE d_poisson16_p_201_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_p_201_u1_2_true_false_true_cost
SUBROUTINE d_poisson16_p_201_a_u4_0_true_false_false(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  integer(kind=4), dimension(-8 - (8):8 - (-8) - (1)) :: mod_arr
  do l = -8 - (8), 8 - (-8) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1, tt2, tt3)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (4), 4
      do i2 = 0,  -(-8) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2), i3)) * (poisson16_8_fil(l))
          tt1 = tt1 + (x(i1 + 1, mod_arr(l + i2), i3)) * (poisson16_8_fil(l))
          tt2 = tt2 + (x(i1 + 2, mod_arr(l + i2), i3)) * (poisson16_8_fil(l))
          tt3 = tt3 + (x(i1 + 3, mod_arr(l + i2), i3)) * (poisson16_8_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
      do i2 =  -(-8), n - (8) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson16_8_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2, i3)) * (poisson16_8_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
      do i2 = n - (8), n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2 - (n)), i3)) * (poisson16_8_fil(l))
          tt1 = tt1 + (x(i1 + 1, mod_arr(l + i2 - (n)), i3)) * (poisson16_8_fil(l))
          tt2 = tt2 + (x(i1 + 2, mod_arr(l + i2 - (n)), i3)) * (poisson16_8_fil(l))
          tt3 = tt3 + (x(i1 + 3, mod_arr(l + i2 - (n)), i3)) * (poisson16_8_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (4)) * (4), ndat1 - (1), 1
      do i2 = 0,  -(-8) - (1), 1
        tt0 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2), i3)) * (poisson16_8_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-8), n - (8) - (1), 1
        tt0 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (8), n - (1), 1
        tt0 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, mod_arr(l + i2 - (n)), i3)) * (poisson16_8_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_p_201_a_u4_0_true_false_false
SUBROUTINE d_poisson16_p_201_a_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_p_201_a_u1_2_true_false_true_cost
SUBROUTINE d_poisson16_fg_10_u1_1_false_true_true(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension( -(8):n - (-8) - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:0) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 =  -(8),  -(-8) - (1), 1
      tt(0) = 0.0_wp
      do l = max( -(i1), -8), 8, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson16_8_fil(l))
      end do
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-8), n - (8) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-8 + i1, i2 + 0)) * (9.712509712509713e-06_wp)
      tt(0) = tt(0) + (x(-7 + i1, i2 + 0)) * (-0.0001776001776001776_wp)
      tt(0) = tt(0) + (x(-6 + i1, i2 + 0)) * (0.001554001554001554_wp)
      tt(0) = tt(0) + (x(-5 + i1, i2 + 0)) * (-0.008702408702408702_wp)
      tt(0) = tt(0) + (x(-4 + i1, i2 + 0)) * (0.03535353535353535_wp)
      tt(0) = tt(0) + (x(-3 + i1, i2 + 0)) * (-0.11313131313131314_wp)
      tt(0) = tt(0) + (x(-2 + i1, i2 + 0)) * (0.3111111111111111_wp)
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.8888888888888888_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.8888888888888888_wp)
      tt(0) = tt(0) + (x(2 + i1, i2 + 0)) * (-0.3111111111111111_wp)
      tt(0) = tt(0) + (x(3 + i1, i2 + 0)) * (0.11313131313131314_wp)
      tt(0) = tt(0) + (x(4 + i1, i2 + 0)) * (-0.03535353535353535_wp)
      tt(0) = tt(0) + (x(5 + i1, i2 + 0)) * (0.008702408702408702_wp)
      tt(0) = tt(0) + (x(6 + i1, i2 + 0)) * (-0.001554001554001554_wp)
      tt(0) = tt(0) + (x(7 + i1, i2 + 0)) * (0.0001776001776001776_wp)
      tt(0) = tt(0) + (x(8 + i1, i2 + 0)) * (-9.712509712509713e-06_wp)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (8), n - (-8) - (1), 1
      tt(0) = 0.0_wp
      do l = -8, min(8, n - (1) - (i1)), 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson16_8_fil(l))
      end do
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_fg_10_u1_1_false_true_true
SUBROUTINE d_poisson16_fg_10_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_fg_10_u1_1_false_false_true_cost
SUBROUTINE d_poisson16_fg_10_a_u1_1_false_true_true(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension( -(8):n - (-8) - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:0) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 =  -(8),  -(-8) - (1), 1
      tt(0) = 0.0_wp
      do l = max( -(i1), -8), 8, 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson16_8_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 =  -(-8), n - (8) - (1), 1
      tt(0) = 0.0_wp
      tt(0) = tt(0) + (x(-8 + i1, i2 + 0)) * (9.712509712509713e-06_wp)
      tt(0) = tt(0) + (x(-7 + i1, i2 + 0)) * (-0.0001776001776001776_wp)
      tt(0) = tt(0) + (x(-6 + i1, i2 + 0)) * (0.001554001554001554_wp)
      tt(0) = tt(0) + (x(-5 + i1, i2 + 0)) * (-0.008702408702408702_wp)
      tt(0) = tt(0) + (x(-4 + i1, i2 + 0)) * (0.03535353535353535_wp)
      tt(0) = tt(0) + (x(-3 + i1, i2 + 0)) * (-0.11313131313131314_wp)
      tt(0) = tt(0) + (x(-2 + i1, i2 + 0)) * (0.3111111111111111_wp)
      tt(0) = tt(0) + (x(-1 + i1, i2 + 0)) * (-0.8888888888888888_wp)
      tt(0) = tt(0) + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt(0) = tt(0) + (x(1 + i1, i2 + 0)) * (0.8888888888888888_wp)
      tt(0) = tt(0) + (x(2 + i1, i2 + 0)) * (-0.3111111111111111_wp)
      tt(0) = tt(0) + (x(3 + i1, i2 + 0)) * (0.11313131313131314_wp)
      tt(0) = tt(0) + (x(4 + i1, i2 + 0)) * (-0.03535353535353535_wp)
      tt(0) = tt(0) + (x(5 + i1, i2 + 0)) * (0.008702408702408702_wp)
      tt(0) = tt(0) + (x(6 + i1, i2 + 0)) * (-0.001554001554001554_wp)
      tt(0) = tt(0) + (x(7 + i1, i2 + 0)) * (0.0001776001776001776_wp)
      tt(0) = tt(0) + (x(8 + i1, i2 + 0)) * (-9.712509712509713e-06_wp)
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
    do i1 = n - (8), n - (-8) - (1), 1
      tt(0) = 0.0_wp
      do l = -8, min(8, n - (1) - (i1)), 1
        tt(0) = tt(0) + (x(l + i1, i2 + 0)) * (poisson16_8_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1, i2 + 0) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_fg_10_a_u1_1_false_true_true
SUBROUTINE d_poisson16_fg_10_a_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_fg_10_a_u1_1_false_false_true_cost
SUBROUTINE d_poisson16_fg_01_u5_0_false_true_false(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1),  -(8):n - (-8) - (1)) :: y
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:4) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (5), 5
    do i2 =  -(8),  -(-8) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      do l = max( -(i2), -8), 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson16_8_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson16_8_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson16_8_fil(l))
        tt(4) = tt(4) + (x(i1 + 4, l + i2)) * (poisson16_8_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
      y(i1 + 2, i2) = tt(2)
      y(i1 + 3, i2) = tt(3)
      y(i1 + 4, i2) = tt(4)
    end do
    do i2 =  -(-8), n - (8) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson16_8_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson16_8_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson16_8_fil(l))
        tt(4) = tt(4) + (x(i1 + 4, l + i2)) * (poisson16_8_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
      y(i1 + 2, i2) = tt(2)
      y(i1 + 3, i2) = tt(3)
      y(i1 + 4, i2) = tt(4)
    end do
    do i2 = n - (8), n - (-8) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      do l = -8, min(8, n - (1) - (i2)), 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson16_8_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson16_8_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson16_8_fil(l))
        tt(4) = tt(4) + (x(i1 + 4, l + i2)) * (poisson16_8_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
      y(i1 + 2, i2) = tt(2)
      y(i1 + 3, i2) = tt(3)
      y(i1 + 4, i2) = tt(4)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (5)) * (5), ndat - (1), 1
    do i2 =  -(8),  -(-8) - (1), 1
      tt(0) = 0.0_wp
      do l = max( -(i2), -8), 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 =  -(-8), n - (8) - (1), 1
      tt(0) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 = n - (8), n - (-8) - (1), 1
      tt(0) = 0.0_wp
      do l = -8, min(8, n - (1) - (i2)), 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_fg_01_u5_0_false_true_false
SUBROUTINE d_poisson16_fg_01_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_fg_01_u1_0_false_false_true_cost
SUBROUTINE d_poisson16_fg_01_a_u5_0_false_true_false(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1),  -(8):n - (-8) - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:4) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (5), 5
    do i2 =  -(8),  -(-8) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      do l = max( -(i2), -8), 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson16_8_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson16_8_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson16_8_fil(l))
        tt(4) = tt(4) + (x(i1 + 4, l + i2)) * (poisson16_8_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
      tt(4) = (tt(4)) * (a)
      y(i1 + 4, i2) = tt(4)
    end do
    do i2 =  -(-8), n - (8) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson16_8_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson16_8_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson16_8_fil(l))
        tt(4) = tt(4) + (x(i1 + 4, l + i2)) * (poisson16_8_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
      tt(4) = (tt(4)) * (a)
      y(i1 + 4, i2) = tt(4)
    end do
    do i2 = n - (8), n - (-8) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      do l = -8, min(8, n - (1) - (i2)), 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson16_8_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson16_8_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson16_8_fil(l))
        tt(4) = tt(4) + (x(i1 + 4, l + i2)) * (poisson16_8_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
      tt(4) = (tt(4)) * (a)
      y(i1 + 4, i2) = tt(4)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (5)) * (5), ndat - (1), 1
    do i2 =  -(8),  -(-8) - (1), 1
      tt(0) = 0.0_wp
      do l = max( -(i2), -8), 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 =  -(-8), n - (8) - (1), 1
      tt(0) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 = n - (8), n - (-8) - (1), 1
      tt(0) = 0.0_wp
      do l = -8, min(8, n - (1) - (i2)), 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_fg_01_a_u5_0_false_true_false
SUBROUTINE d_poisson16_fg_01_a_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_fg_01_a_u1_0_false_false_true_cost
SUBROUTINE d_poisson16_fg_201_u5_0_false_false_false(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1),  -(8):n - (-8) - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1, tt2, tt3, tt4)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (5), 5
      do i2 =  -(8),  -(-8) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt4 = 0.0_wp
        do l = max( -(i2), -8), 8, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson16_8_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2, i3)) * (poisson16_8_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2, i3)) * (poisson16_8_fil(l))
          tt4 = tt4 + (x(i1 + 4, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
        y(i1 + 2, i2, i3) = tt2
        y(i1 + 3, i2, i3) = tt3
        y(i1 + 4, i2, i3) = tt4
      end do
      do i2 =  -(-8), n - (8) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt4 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson16_8_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2, i3)) * (poisson16_8_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2, i3)) * (poisson16_8_fil(l))
          tt4 = tt4 + (x(i1 + 4, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
        y(i1 + 2, i2, i3) = tt2
        y(i1 + 3, i2, i3) = tt3
        y(i1 + 4, i2, i3) = tt4
      end do
      do i2 = n - (8), n - (-8) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt4 = 0.0_wp
        do l = -8, min(8, n - (1) - (i2)), 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson16_8_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2, i3)) * (poisson16_8_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2, i3)) * (poisson16_8_fil(l))
          tt4 = tt4 + (x(i1 + 4, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
        y(i1 + 2, i2, i3) = tt2
        y(i1 + 3, i2, i3) = tt3
        y(i1 + 4, i2, i3) = tt4
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (5)) * (5), ndat1 - (1), 1
      do i2 =  -(8),  -(-8) - (1), 1
        tt0 = 0.0_wp
        do l = max( -(i2), -8), 8, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-8), n - (8) - (1), 1
        tt0 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (8), n - (-8) - (1), 1
        tt0 = 0.0_wp
        do l = -8, min(8, n - (1) - (i2)), 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_fg_201_u5_0_false_false_false
SUBROUTINE d_poisson16_fg_201_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_fg_201_u1_2_false_false_true_cost
SUBROUTINE d_poisson16_fg_201_a_u5_0_false_false_false(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1),  -(8):n - (-8) - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1, tt2, tt3, tt4)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (5), 5
      do i2 =  -(8),  -(-8) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt4 = 0.0_wp
        do l = max( -(i2), -8), 8, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson16_8_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2, i3)) * (poisson16_8_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2, i3)) * (poisson16_8_fil(l))
          tt4 = tt4 + (x(i1 + 4, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
        tt4 = (tt4) * (a)
        y(i1 + 4, i2, i3) = tt4
      end do
      do i2 =  -(-8), n - (8) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt4 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson16_8_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2, i3)) * (poisson16_8_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2, i3)) * (poisson16_8_fil(l))
          tt4 = tt4 + (x(i1 + 4, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
        tt4 = (tt4) * (a)
        y(i1 + 4, i2, i3) = tt4
      end do
      do i2 = n - (8), n - (-8) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt4 = 0.0_wp
        do l = -8, min(8, n - (1) - (i2)), 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson16_8_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2, i3)) * (poisson16_8_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2, i3)) * (poisson16_8_fil(l))
          tt4 = tt4 + (x(i1 + 4, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
        tt4 = (tt4) * (a)
        y(i1 + 4, i2, i3) = tt4
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (5)) * (5), ndat1 - (1), 1
      do i2 =  -(8),  -(-8) - (1), 1
        tt0 = 0.0_wp
        do l = max( -(i2), -8), 8, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-8), n - (8) - (1), 1
        tt0 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (8), n - (-8) - (1), 1
        tt0 = 0.0_wp
        do l = -8, min(8, n - (1) - (i2)), 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_fg_201_a_u5_0_false_false_false
SUBROUTINE d_poisson16_fg_201_a_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_fg_201_a_u1_2_false_false_true_cost
SUBROUTINE d_poisson16_fs_10_u2_1_false_false_false(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(-8:n + 8 - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
!$omp parallel  default(shared) private(i1, i2, tt0, tt1)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 = 0, n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson16_8_fil(l))
        tt1 = tt1 + (x(l + i1, i2 + 1)) * (poisson16_8_fil(l))
      end do
      y(i1, i2 + 0) = tt0
      y(i1, i2 + 1) = tt1
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 = 0, n - (1), 1
      tt0 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson16_8_fil(l))
      end do
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_fs_10_u2_1_false_false_false
SUBROUTINE d_poisson16_fs_10_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_fs_10_u1_1_false_false_true_cost
SUBROUTINE d_poisson16_fs_10_a_u2_1_false_false_false(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(-8:n + 8 - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
!$omp parallel  default(shared) private(i1, i2, tt0, tt1)
!$omp do 
  do i2 = 0, ndat - (2), 2
    do i1 = 0, n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson16_8_fil(l))
        tt1 = tt1 + (x(l + i1, i2 + 1)) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
      tt1 = (tt1) * (a)
      y(i1, i2 + 1) = tt1
    end do
  end do
!$omp end do 
!$omp do 
  do i2 = ((ndat) / (2)) * (2), ndat - (1), 1
    do i1 = 0, n - (1), 1
      tt0 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_fs_10_a_u2_1_false_false_false
SUBROUTINE d_poisson16_fs_10_a_u1_1_false_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_fs_10_a_u1_1_false_false_true_cost
SUBROUTINE d_poisson16_fs_01_u5_0_false_true_false(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), -8:n + 8 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:4) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (5), 5
    do i2 = 0, n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson16_8_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson16_8_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson16_8_fil(l))
        tt(4) = tt(4) + (x(i1 + 4, l + i2)) * (poisson16_8_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
      y(i1 + 2, i2) = tt(2)
      y(i1 + 3, i2) = tt(3)
      y(i1 + 4, i2) = tt(4)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (5)) * (5), ndat - (1), 1
    do i2 = 0, n - (1), 1
      tt(0) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_fs_01_u5_0_false_true_false
SUBROUTINE d_poisson16_fs_01_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_fs_01_u1_0_false_false_true_cost
SUBROUTINE d_poisson16_fs_01_a_u5_0_false_false_false(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), -8:n + 8 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
!$omp parallel  default(shared) private(i1, i2, tt0, tt1, tt2, tt3, tt4)
!$omp do 
  do i1 = 0, ndat - (5), 5
    do i2 = 0, n - (1), 1
      tt0 = 0.0_wp
      tt1 = 0.0_wp
      tt2 = 0.0_wp
      tt3 = 0.0_wp
      tt4 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
        tt1 = tt1 + (x(i1 + 1, l + i2)) * (poisson16_8_fil(l))
        tt2 = tt2 + (x(i1 + 2, l + i2)) * (poisson16_8_fil(l))
        tt3 = tt3 + (x(i1 + 3, l + i2)) * (poisson16_8_fil(l))
        tt4 = tt4 + (x(i1 + 4, l + i2)) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
      tt1 = (tt1) * (a)
      y(i1 + 1, i2) = tt1
      tt2 = (tt2) * (a)
      y(i1 + 2, i2) = tt2
      tt3 = (tt3) * (a)
      y(i1 + 3, i2) = tt3
      tt4 = (tt4) * (a)
      y(i1 + 4, i2) = tt4
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (5)) * (5), ndat - (1), 1
    do i2 = 0, n - (1), 1
      tt0 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1 + 0, i2) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_fs_01_a_u5_0_false_false_false
SUBROUTINE d_poisson16_fs_01_a_u1_0_false_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_fs_01_a_u1_0_false_false_true_cost
SUBROUTINE d_poisson16_fs_201_u5_0_false_false_false(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), -8:n + 8 - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1, tt2, tt3, tt4)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (5), 5
      do i2 = 0, n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt4 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson16_8_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2, i3)) * (poisson16_8_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2, i3)) * (poisson16_8_fil(l))
          tt4 = tt4 + (x(i1 + 4, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
        y(i1 + 1, i2, i3) = tt1
        y(i1 + 2, i2, i3) = tt2
        y(i1 + 3, i2, i3) = tt3
        y(i1 + 4, i2, i3) = tt4
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (5)) * (5), ndat1 - (1), 1
      do i2 = 0, n - (1), 1
        tt0 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_fs_201_u5_0_false_false_false
SUBROUTINE d_poisson16_fs_201_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_fs_201_u1_2_false_false_true_cost
SUBROUTINE d_poisson16_fs_201_a_u5_0_false_true_false(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), -8:n + 8 - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8), dimension(0:4) :: tt
!$omp parallel  default(shared) private(i1, i2, i3, tt)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (5), 5
      do i2 = 0, n - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        tt(4) = 0.0_wp
        do l = -8, 8, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, l + i2, i3)) * (poisson16_8_fil(l))
          tt(2) = tt(2) + (x(i1 + 2, l + i2, i3)) * (poisson16_8_fil(l))
          tt(3) = tt(3) + (x(i1 + 3, l + i2, i3)) * (poisson16_8_fil(l))
          tt(4) = tt(4) + (x(i1 + 4, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
        tt(1) = (tt(1)) * (a)
        y(i1 + 1, i2, i3) = tt(1)
        tt(2) = (tt(2)) * (a)
        y(i1 + 2, i2, i3) = tt(2)
        tt(3) = (tt(3)) * (a)
        y(i1 + 3, i2, i3) = tt(3)
        tt(4) = (tt(4)) * (a)
        y(i1 + 4, i2, i3) = tt(4)
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (5)) * (5), ndat1 - (1), 1
      do i2 = 0, n - (1), 1
        tt(0) = 0.0_wp
        do l = -8, 8, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        tt(0) = (tt(0)) * (a)
        y(i1 + 0, i2, i3) = tt(0)
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_fs_201_a_u5_0_false_true_false
SUBROUTINE d_poisson16_fs_201_a_u1_2_false_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_fs_201_a_u1_2_false_false_true_cost
SUBROUTINE d_poisson16_np_10_u1_1_false_false_true(n, ndat, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), parameter, dimension(0:288) :: poisson16_fil = (/ &
-3.3807289932289932_wp, &
16.0_wp, &
-60.0_wp, &
186.66666666666666_wp, &
-455.0_wp, &
873.6_wp, &
-1334.6666666666667_wp, &
1634.2857142857142_wp, &
-1608.75_wp, &
1271.111111111111_wp, &
-800.8_wp, &
397.09090909090907_wp, &
-151.66666666666666_wp, &
43.07692307692308_wp, &
-8.571428571428571_wp, &
1.0666666666666667_wp, &
-0.0625_wp, &
-0.0625_wp, &
-2.3182289932289932_wp, &
7.5_wp, &
-17.5_wp, &
37.916666666666664_wp, &
-68.25_wp, &
100.1_wp, &
-119.16666666666667_wp, &
114.91071428571429_wp, &
-89.375_wp, &
55.611111111111114_wp, &
-27.3_wp, &
10.340909090909092_wp, &
-2.9166666666666665_wp, &
0.5769230769230769_wp, &
-0.07142857142857142_wp, &
0.004166666666666667_wp, &
0.004166666666666667_wp, &
-0.13333333333333333_wp, &
-1.7515623265623266_wp, &
4.666666666666667_wp, &
-7.583333333333333_wp, &
12.133333333333333_wp, &
-16.683333333333334_wp, &
19.066666666666666_wp, &
-17.875_wp, &
13.619047619047619_wp, &
-8.341666666666667_wp, &
4.044444444444444_wp, &
-1.5166666666666666_wp, &
0.42424242424242425_wp, &
-0.08333333333333333_wp, &
0.010256410256410256_wp, &
-0.0005952380952380953_wp, &
-0.0005952380952380953_wp, &
0.014285714285714285_wp, &
-0.21428571428571427_wp, &
-1.3468004218004217_wp, &
3.25_wp, &
-3.9_wp, &
4.766666666666667_wp, &
-5.107142857142857_wp, &
4.5964285714285715_wp, &
-3.4047619047619047_wp, &
2.0428571428571427_wp, &
-0.975_wp, &
0.3611111111111111_wp, &
-0.1_wp, &
0.01948051948051948_wp, &
-0.002380952380952381_wp, &
0.00013736263736263736_wp, &
0.00013736263736263736_wp, &
-0.0029304029304029304_wp, &
0.03296703296703297_wp, &
-0.3076923076923077_wp, &
-1.0198773448773448_wp, &
2.4_wp, &
-2.2_wp, &
2.0952380952380953_wp, &
-1.7678571428571428_wp, &
1.2571428571428571_wp, &
-0.7333333333333333_wp, &
0.34285714285714286_wp, &
-0.125_wp, &
0.03418803418803419_wp, &
-0.006593406593406593_wp, &
0.0007992007992007992_wp, &
-4.578754578754579e-05_wp, &
-4.578754578754579e-05_wp, &
0.0009157509157509158_wp, &
-0.009157509157509158_wp, &
0.0641025641025641_wp, &
-0.4166666666666667_wp, &
-0.7365440115440115_wp, &
1.8333333333333333_wp, &
-1.3095238095238095_wp, &
0.9821428571428571_wp, &
-0.6547619047619048_wp, &
0.36666666666666664_wp, &
-0.16666666666666666_wp, &
0.05952380952380952_wp, &
-0.016025641025641024_wp, &
0.0030525030525030525_wp, &
-0.0003663003663003663_wp, &
2.0812520812520813e-05_wp, &
2.0812520812520813e-05_wp, &
-0.0003996003996003996_wp, &
0.0037462537462537465_wp, &
-0.023310023310023312_wp, &
0.11363636363636363_wp, &
-0.5454545454545454_wp, &
-0.478968253968254_wp, &
1.4285714285714286_wp, &
-0.8035714285714286_wp, &
0.47619047619047616_wp, &
-0.25_wp, &
0.10909090909090909_wp, &
-0.03787878787878788_wp, &
0.00999000999000999_wp, &
-0.0018731268731268732_wp, &
0.000222000222000222_wp, &
-1.2487512487512488e-05_wp, &
-1.2487512487512488e-05_wp, &
0.0002331002331002331_wp, &
-0.002097902097902098_wp, &
0.012237762237762238_wp, &
-0.05303030303030303_wp, &
0.19090909090909092_wp, &
-0.7_wp, &
-0.2361111111111111_wp, &
1.125_wp, &
-0.5_wp, &
0.23333333333333334_wp, &
-0.09545454545454546_wp, &
0.031818181818181815_wp, &
-0.008158508158508158_wp, &
0.0014985014985014985_wp, &
-0.00017482517482517483_wp, &
9.712509712509713e-06_wp, &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp, &
-9.712509712509713e-06_wp, &
0.00017482517482517483_wp, &
-0.0014985014985014985_wp, &
0.008158508158508158_wp, &
-0.031818181818181815_wp, &
0.09545454545454546_wp, &
-0.23333333333333334_wp, &
0.5_wp, &
-1.125_wp, &
0.2361111111111111_wp, &
0.7_wp, &
-0.19090909090909092_wp, &
0.05303030303030303_wp, &
-0.012237762237762238_wp, &
0.002097902097902098_wp, &
-0.0002331002331002331_wp, &
1.2487512487512488e-05_wp, &
1.2487512487512488e-05_wp, &
-0.000222000222000222_wp, &
0.0018731268731268732_wp, &
-0.00999000999000999_wp, &
0.03787878787878788_wp, &
-0.10909090909090909_wp, &
0.25_wp, &
-0.47619047619047616_wp, &
0.8035714285714286_wp, &
-1.4285714285714286_wp, &
0.478968253968254_wp, &
0.5454545454545454_wp, &
-0.11363636363636363_wp, &
0.023310023310023312_wp, &
-0.0037462537462537465_wp, &
0.0003996003996003996_wp, &
-2.0812520812520813e-05_wp, &
-2.0812520812520813e-05_wp, &
0.0003663003663003663_wp, &
-0.0030525030525030525_wp, &
0.016025641025641024_wp, &
-0.05952380952380952_wp, &
0.16666666666666666_wp, &
-0.36666666666666664_wp, &
0.6547619047619048_wp, &
-0.9821428571428571_wp, &
1.3095238095238095_wp, &
-1.8333333333333333_wp, &
0.7365440115440115_wp, &
0.4166666666666667_wp, &
-0.0641025641025641_wp, &
0.009157509157509158_wp, &
-0.0009157509157509158_wp, &
4.578754578754579e-05_wp, &
4.578754578754579e-05_wp, &
-0.0007992007992007992_wp, &
0.006593406593406593_wp, &
-0.03418803418803419_wp, &
0.125_wp, &
-0.34285714285714286_wp, &
0.7333333333333333_wp, &
-1.2571428571428571_wp, &
1.7678571428571428_wp, &
-2.0952380952380953_wp, &
2.2_wp, &
-2.4_wp, &
1.0198773448773448_wp, &
0.3076923076923077_wp, &
-0.03296703296703297_wp, &
0.0029304029304029304_wp, &
-0.00013736263736263736_wp, &
-0.00013736263736263736_wp, &
0.002380952380952381_wp, &
-0.01948051948051948_wp, &
0.1_wp, &
-0.3611111111111111_wp, &
0.975_wp, &
-2.0428571428571427_wp, &
3.4047619047619047_wp, &
-4.5964285714285715_wp, &
5.107142857142857_wp, &
-4.766666666666667_wp, &
3.9_wp, &
-3.25_wp, &
1.3468004218004217_wp, &
0.21428571428571427_wp, &
-0.014285714285714285_wp, &
0.0005952380952380953_wp, &
0.0005952380952380953_wp, &
-0.010256410256410256_wp, &
0.08333333333333333_wp, &
-0.42424242424242425_wp, &
1.5166666666666666_wp, &
-4.044444444444444_wp, &
8.341666666666667_wp, &
-13.619047619047619_wp, &
17.875_wp, &
-19.066666666666666_wp, &
16.683333333333334_wp, &
-12.133333333333333_wp, &
7.583333333333333_wp, &
-4.666666666666667_wp, &
1.7515623265623266_wp, &
0.13333333333333333_wp, &
-0.004166666666666667_wp, &
-0.004166666666666667_wp, &
0.07142857142857142_wp, &
-0.5769230769230769_wp, &
2.9166666666666665_wp, &
-10.340909090909092_wp, &
27.3_wp, &
-55.611111111111114_wp, &
89.375_wp, &
-114.91071428571429_wp, &
119.16666666666667_wp, &
-100.1_wp, &
68.25_wp, &
-37.916666666666664_wp, &
17.5_wp, &
-7.5_wp, &
2.3182289932289932_wp, &
0.0625_wp, &
0.0625_wp, &
-1.0666666666666667_wp, &
8.571428571428571_wp, &
-43.07692307692308_wp, &
151.66666666666666_wp, &
-397.09090909090907_wp, &
800.8_wp, &
-1271.111111111111_wp, &
1608.75_wp, &
-1634.2857142857142_wp, &
1334.6666666666667_wp, &
-873.6_wp, &
455.0_wp, &
-186.66666666666666_wp, &
60.0_wp, &
-16.0_wp, &
3.3807289932289932_wp /)
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
!$omp parallel  default(shared) private(i1, i2, tt0)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 = 0,  -(-8) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(-8 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + -8 + 8))
      tt0 = tt0 + (x(-7 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + -7 + 8))
      tt0 = tt0 + (x(-6 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + -6 + 8))
      tt0 = tt0 + (x(-5 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + -5 + 8))
      tt0 = tt0 + (x(-4 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + -4 + 8))
      tt0 = tt0 + (x(-3 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + -3 + 8))
      tt0 = tt0 + (x(-2 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + -2 + 8))
      tt0 = tt0 + (x(-1 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + -1 + 8))
      tt0 = tt0 + (x(0 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + 0 + 8))
      tt0 = tt0 + (x(1 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + 1 + 8))
      tt0 = tt0 + (x(2 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + 2 + 8))
      tt0 = tt0 + (x(3 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + 3 + 8))
      tt0 = tt0 + (x(4 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + 4 + 8))
      tt0 = tt0 + (x(5 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + 5 + 8))
      tt0 = tt0 + (x(6 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + 6 + 8))
      tt0 = tt0 + (x(7 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + 7 + 8))
      tt0 = tt0 + (x(8 + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + 8 + 8))
      y(i1, i2 + 0) = tt0
    end do
    do i1 =  -(-8), n - (8) - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(-8 + i1, i2 + 0)) * (9.712509712509713e-06_wp)
      tt0 = tt0 + (x(-7 + i1, i2 + 0)) * (-0.0001776001776001776_wp)
      tt0 = tt0 + (x(-6 + i1, i2 + 0)) * (0.001554001554001554_wp)
      tt0 = tt0 + (x(-5 + i1, i2 + 0)) * (-0.008702408702408702_wp)
      tt0 = tt0 + (x(-4 + i1, i2 + 0)) * (0.03535353535353535_wp)
      tt0 = tt0 + (x(-3 + i1, i2 + 0)) * (-0.11313131313131314_wp)
      tt0 = tt0 + (x(-2 + i1, i2 + 0)) * (0.3111111111111111_wp)
      tt0 = tt0 + (x(-1 + i1, i2 + 0)) * (-0.8888888888888888_wp)
      tt0 = tt0 + (x(0 + i1, i2 + 0)) * (0.0_wp)
      tt0 = tt0 + (x(1 + i1, i2 + 0)) * (0.8888888888888888_wp)
      tt0 = tt0 + (x(2 + i1, i2 + 0)) * (-0.3111111111111111_wp)
      tt0 = tt0 + (x(3 + i1, i2 + 0)) * (0.11313131313131314_wp)
      tt0 = tt0 + (x(4 + i1, i2 + 0)) * (-0.03535353535353535_wp)
      tt0 = tt0 + (x(5 + i1, i2 + 0)) * (0.008702408702408702_wp)
      tt0 = tt0 + (x(6 + i1, i2 + 0)) * (-0.001554001554001554_wp)
      tt0 = tt0 + (x(7 + i1, i2 + 0)) * (0.0001776001776001776_wp)
      tt0 = tt0 + (x(8 + i1, i2 + 0)) * (-9.712509712509713e-06_wp)
      y(i1, i2 + 0) = tt0
    end do
    do i1 = n - (8), n - (1), 1
      tt0 = 0.0_wp
      tt0 = tt0 + (x(-8 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + -8 + 8))
      tt0 = tt0 + (x(-7 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + -7 + 8))
      tt0 = tt0 + (x(-6 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + -6 + 8))
      tt0 = tt0 + (x(-5 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + -5 + 8))
      tt0 = tt0 + (x(-4 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + -4 + 8))
      tt0 = tt0 + (x(-3 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + -3 + 8))
      tt0 = tt0 + (x(-2 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + -2 + 8))
      tt0 = tt0 + (x(-1 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + -1 + 8))
      tt0 = tt0 + (x(0 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + 0 + 8))
      tt0 = tt0 + (x(1 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + 1 + 8))
      tt0 = tt0 + (x(2 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + 2 + 8))
      tt0 = tt0 + (x(3 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + 3 + 8))
      tt0 = tt0 + (x(4 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + 4 + 8))
      tt0 = tt0 + (x(5 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + 5 + 8))
      tt0 = tt0 + (x(6 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + 6 + 8))
      tt0 = tt0 + (x(7 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + 7 + 8))
      tt0 = tt0 + (x(8 + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + 8 + 8))
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_np_10_u1_1_false_false_true
SUBROUTINE d_poisson16_np_10_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_np_10_u1_1_true_false_true_cost
SUBROUTINE d_poisson16_np_10_a_u1_1_true_false_false(n, ndat, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n - (1), 0:ndat - (1)) :: x
  real(kind=8), intent(out), dimension(0:n - (1), 0:ndat - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(0:288) :: poisson16_fil = (/ &
-3.3807289932289932_wp, &
16.0_wp, &
-60.0_wp, &
186.66666666666666_wp, &
-455.0_wp, &
873.6_wp, &
-1334.6666666666667_wp, &
1634.2857142857142_wp, &
-1608.75_wp, &
1271.111111111111_wp, &
-800.8_wp, &
397.09090909090907_wp, &
-151.66666666666666_wp, &
43.07692307692308_wp, &
-8.571428571428571_wp, &
1.0666666666666667_wp, &
-0.0625_wp, &
-0.0625_wp, &
-2.3182289932289932_wp, &
7.5_wp, &
-17.5_wp, &
37.916666666666664_wp, &
-68.25_wp, &
100.1_wp, &
-119.16666666666667_wp, &
114.91071428571429_wp, &
-89.375_wp, &
55.611111111111114_wp, &
-27.3_wp, &
10.340909090909092_wp, &
-2.9166666666666665_wp, &
0.5769230769230769_wp, &
-0.07142857142857142_wp, &
0.004166666666666667_wp, &
0.004166666666666667_wp, &
-0.13333333333333333_wp, &
-1.7515623265623266_wp, &
4.666666666666667_wp, &
-7.583333333333333_wp, &
12.133333333333333_wp, &
-16.683333333333334_wp, &
19.066666666666666_wp, &
-17.875_wp, &
13.619047619047619_wp, &
-8.341666666666667_wp, &
4.044444444444444_wp, &
-1.5166666666666666_wp, &
0.42424242424242425_wp, &
-0.08333333333333333_wp, &
0.010256410256410256_wp, &
-0.0005952380952380953_wp, &
-0.0005952380952380953_wp, &
0.014285714285714285_wp, &
-0.21428571428571427_wp, &
-1.3468004218004217_wp, &
3.25_wp, &
-3.9_wp, &
4.766666666666667_wp, &
-5.107142857142857_wp, &
4.5964285714285715_wp, &
-3.4047619047619047_wp, &
2.0428571428571427_wp, &
-0.975_wp, &
0.3611111111111111_wp, &
-0.1_wp, &
0.01948051948051948_wp, &
-0.002380952380952381_wp, &
0.00013736263736263736_wp, &
0.00013736263736263736_wp, &
-0.0029304029304029304_wp, &
0.03296703296703297_wp, &
-0.3076923076923077_wp, &
-1.0198773448773448_wp, &
2.4_wp, &
-2.2_wp, &
2.0952380952380953_wp, &
-1.7678571428571428_wp, &
1.2571428571428571_wp, &
-0.7333333333333333_wp, &
0.34285714285714286_wp, &
-0.125_wp, &
0.03418803418803419_wp, &
-0.006593406593406593_wp, &
0.0007992007992007992_wp, &
-4.578754578754579e-05_wp, &
-4.578754578754579e-05_wp, &
0.0009157509157509158_wp, &
-0.009157509157509158_wp, &
0.0641025641025641_wp, &
-0.4166666666666667_wp, &
-0.7365440115440115_wp, &
1.8333333333333333_wp, &
-1.3095238095238095_wp, &
0.9821428571428571_wp, &
-0.6547619047619048_wp, &
0.36666666666666664_wp, &
-0.16666666666666666_wp, &
0.05952380952380952_wp, &
-0.016025641025641024_wp, &
0.0030525030525030525_wp, &
-0.0003663003663003663_wp, &
2.0812520812520813e-05_wp, &
2.0812520812520813e-05_wp, &
-0.0003996003996003996_wp, &
0.0037462537462537465_wp, &
-0.023310023310023312_wp, &
0.11363636363636363_wp, &
-0.5454545454545454_wp, &
-0.478968253968254_wp, &
1.4285714285714286_wp, &
-0.8035714285714286_wp, &
0.47619047619047616_wp, &
-0.25_wp, &
0.10909090909090909_wp, &
-0.03787878787878788_wp, &
0.00999000999000999_wp, &
-0.0018731268731268732_wp, &
0.000222000222000222_wp, &
-1.2487512487512488e-05_wp, &
-1.2487512487512488e-05_wp, &
0.0002331002331002331_wp, &
-0.002097902097902098_wp, &
0.012237762237762238_wp, &
-0.05303030303030303_wp, &
0.19090909090909092_wp, &
-0.7_wp, &
-0.2361111111111111_wp, &
1.125_wp, &
-0.5_wp, &
0.23333333333333334_wp, &
-0.09545454545454546_wp, &
0.031818181818181815_wp, &
-0.008158508158508158_wp, &
0.0014985014985014985_wp, &
-0.00017482517482517483_wp, &
9.712509712509713e-06_wp, &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp, &
-9.712509712509713e-06_wp, &
0.00017482517482517483_wp, &
-0.0014985014985014985_wp, &
0.008158508158508158_wp, &
-0.031818181818181815_wp, &
0.09545454545454546_wp, &
-0.23333333333333334_wp, &
0.5_wp, &
-1.125_wp, &
0.2361111111111111_wp, &
0.7_wp, &
-0.19090909090909092_wp, &
0.05303030303030303_wp, &
-0.012237762237762238_wp, &
0.002097902097902098_wp, &
-0.0002331002331002331_wp, &
1.2487512487512488e-05_wp, &
1.2487512487512488e-05_wp, &
-0.000222000222000222_wp, &
0.0018731268731268732_wp, &
-0.00999000999000999_wp, &
0.03787878787878788_wp, &
-0.10909090909090909_wp, &
0.25_wp, &
-0.47619047619047616_wp, &
0.8035714285714286_wp, &
-1.4285714285714286_wp, &
0.478968253968254_wp, &
0.5454545454545454_wp, &
-0.11363636363636363_wp, &
0.023310023310023312_wp, &
-0.0037462537462537465_wp, &
0.0003996003996003996_wp, &
-2.0812520812520813e-05_wp, &
-2.0812520812520813e-05_wp, &
0.0003663003663003663_wp, &
-0.0030525030525030525_wp, &
0.016025641025641024_wp, &
-0.05952380952380952_wp, &
0.16666666666666666_wp, &
-0.36666666666666664_wp, &
0.6547619047619048_wp, &
-0.9821428571428571_wp, &
1.3095238095238095_wp, &
-1.8333333333333333_wp, &
0.7365440115440115_wp, &
0.4166666666666667_wp, &
-0.0641025641025641_wp, &
0.009157509157509158_wp, &
-0.0009157509157509158_wp, &
4.578754578754579e-05_wp, &
4.578754578754579e-05_wp, &
-0.0007992007992007992_wp, &
0.006593406593406593_wp, &
-0.03418803418803419_wp, &
0.125_wp, &
-0.34285714285714286_wp, &
0.7333333333333333_wp, &
-1.2571428571428571_wp, &
1.7678571428571428_wp, &
-2.0952380952380953_wp, &
2.2_wp, &
-2.4_wp, &
1.0198773448773448_wp, &
0.3076923076923077_wp, &
-0.03296703296703297_wp, &
0.0029304029304029304_wp, &
-0.00013736263736263736_wp, &
-0.00013736263736263736_wp, &
0.002380952380952381_wp, &
-0.01948051948051948_wp, &
0.1_wp, &
-0.3611111111111111_wp, &
0.975_wp, &
-2.0428571428571427_wp, &
3.4047619047619047_wp, &
-4.5964285714285715_wp, &
5.107142857142857_wp, &
-4.766666666666667_wp, &
3.9_wp, &
-3.25_wp, &
1.3468004218004217_wp, &
0.21428571428571427_wp, &
-0.014285714285714285_wp, &
0.0005952380952380953_wp, &
0.0005952380952380953_wp, &
-0.010256410256410256_wp, &
0.08333333333333333_wp, &
-0.42424242424242425_wp, &
1.5166666666666666_wp, &
-4.044444444444444_wp, &
8.341666666666667_wp, &
-13.619047619047619_wp, &
17.875_wp, &
-19.066666666666666_wp, &
16.683333333333334_wp, &
-12.133333333333333_wp, &
7.583333333333333_wp, &
-4.666666666666667_wp, &
1.7515623265623266_wp, &
0.13333333333333333_wp, &
-0.004166666666666667_wp, &
-0.004166666666666667_wp, &
0.07142857142857142_wp, &
-0.5769230769230769_wp, &
2.9166666666666665_wp, &
-10.340909090909092_wp, &
27.3_wp, &
-55.611111111111114_wp, &
89.375_wp, &
-114.91071428571429_wp, &
119.16666666666667_wp, &
-100.1_wp, &
68.25_wp, &
-37.916666666666664_wp, &
17.5_wp, &
-7.5_wp, &
2.3182289932289932_wp, &
0.0625_wp, &
0.0625_wp, &
-1.0666666666666667_wp, &
8.571428571428571_wp, &
-43.07692307692308_wp, &
151.66666666666666_wp, &
-397.09090909090907_wp, &
800.8_wp, &
-1271.111111111111_wp, &
1608.75_wp, &
-1634.2857142857142_wp, &
1334.6666666666667_wp, &
-873.6_wp, &
455.0_wp, &
-186.66666666666666_wp, &
60.0_wp, &
-16.0_wp, &
3.3807289932289932_wp /)
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8) :: tt0
  integer(kind=4), dimension(-8 - (8):8 - (-8) - (1)) :: mod_arr
  do l = -8 - (8), 8 - (-8) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, tt0)
!$omp do 
  do i2 = 0, ndat - (1), 1
    do i1 = 0,  -(-8) - (1), 1
      tt0 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(l + 8, i2 + 0)) * (poisson16_fil((i1 - (8)) * (17) + 136 + l + 8))
      end do
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
    end do
    do i1 =  -(-8), n - (8) - (1), 1
      tt0 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(l + i1, i2 + 0)) * (poisson16_8_fil(l))
      end do
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
    end do
    do i1 = n - (8), n - (1), 1
      tt0 = 0.0_wp
      do l = -8, 8, 1
        tt0 = tt0 + (x(l + -8 + n - (1), i2 + 0)) * (poisson16_fil((i1 + 8 - (n) + 1) * (17) + 136 + l + 8))
      end do
      tt0 = (tt0) * (a)
      y(i1, i2 + 0) = tt0
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_np_10_a_u1_1_true_false_false
SUBROUTINE d_poisson16_np_10_a_u1_1_true_false_true_cost(n, ndat, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_np_10_a_u1_1_true_false_true_cost
SUBROUTINE d_poisson16_np_01_u5_0_false_true_false(ndat, n, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), parameter, dimension(0:288) :: poisson16_fil = (/ &
-3.3807289932289932_wp, &
16.0_wp, &
-60.0_wp, &
186.66666666666666_wp, &
-455.0_wp, &
873.6_wp, &
-1334.6666666666667_wp, &
1634.2857142857142_wp, &
-1608.75_wp, &
1271.111111111111_wp, &
-800.8_wp, &
397.09090909090907_wp, &
-151.66666666666666_wp, &
43.07692307692308_wp, &
-8.571428571428571_wp, &
1.0666666666666667_wp, &
-0.0625_wp, &
-0.0625_wp, &
-2.3182289932289932_wp, &
7.5_wp, &
-17.5_wp, &
37.916666666666664_wp, &
-68.25_wp, &
100.1_wp, &
-119.16666666666667_wp, &
114.91071428571429_wp, &
-89.375_wp, &
55.611111111111114_wp, &
-27.3_wp, &
10.340909090909092_wp, &
-2.9166666666666665_wp, &
0.5769230769230769_wp, &
-0.07142857142857142_wp, &
0.004166666666666667_wp, &
0.004166666666666667_wp, &
-0.13333333333333333_wp, &
-1.7515623265623266_wp, &
4.666666666666667_wp, &
-7.583333333333333_wp, &
12.133333333333333_wp, &
-16.683333333333334_wp, &
19.066666666666666_wp, &
-17.875_wp, &
13.619047619047619_wp, &
-8.341666666666667_wp, &
4.044444444444444_wp, &
-1.5166666666666666_wp, &
0.42424242424242425_wp, &
-0.08333333333333333_wp, &
0.010256410256410256_wp, &
-0.0005952380952380953_wp, &
-0.0005952380952380953_wp, &
0.014285714285714285_wp, &
-0.21428571428571427_wp, &
-1.3468004218004217_wp, &
3.25_wp, &
-3.9_wp, &
4.766666666666667_wp, &
-5.107142857142857_wp, &
4.5964285714285715_wp, &
-3.4047619047619047_wp, &
2.0428571428571427_wp, &
-0.975_wp, &
0.3611111111111111_wp, &
-0.1_wp, &
0.01948051948051948_wp, &
-0.002380952380952381_wp, &
0.00013736263736263736_wp, &
0.00013736263736263736_wp, &
-0.0029304029304029304_wp, &
0.03296703296703297_wp, &
-0.3076923076923077_wp, &
-1.0198773448773448_wp, &
2.4_wp, &
-2.2_wp, &
2.0952380952380953_wp, &
-1.7678571428571428_wp, &
1.2571428571428571_wp, &
-0.7333333333333333_wp, &
0.34285714285714286_wp, &
-0.125_wp, &
0.03418803418803419_wp, &
-0.006593406593406593_wp, &
0.0007992007992007992_wp, &
-4.578754578754579e-05_wp, &
-4.578754578754579e-05_wp, &
0.0009157509157509158_wp, &
-0.009157509157509158_wp, &
0.0641025641025641_wp, &
-0.4166666666666667_wp, &
-0.7365440115440115_wp, &
1.8333333333333333_wp, &
-1.3095238095238095_wp, &
0.9821428571428571_wp, &
-0.6547619047619048_wp, &
0.36666666666666664_wp, &
-0.16666666666666666_wp, &
0.05952380952380952_wp, &
-0.016025641025641024_wp, &
0.0030525030525030525_wp, &
-0.0003663003663003663_wp, &
2.0812520812520813e-05_wp, &
2.0812520812520813e-05_wp, &
-0.0003996003996003996_wp, &
0.0037462537462537465_wp, &
-0.023310023310023312_wp, &
0.11363636363636363_wp, &
-0.5454545454545454_wp, &
-0.478968253968254_wp, &
1.4285714285714286_wp, &
-0.8035714285714286_wp, &
0.47619047619047616_wp, &
-0.25_wp, &
0.10909090909090909_wp, &
-0.03787878787878788_wp, &
0.00999000999000999_wp, &
-0.0018731268731268732_wp, &
0.000222000222000222_wp, &
-1.2487512487512488e-05_wp, &
-1.2487512487512488e-05_wp, &
0.0002331002331002331_wp, &
-0.002097902097902098_wp, &
0.012237762237762238_wp, &
-0.05303030303030303_wp, &
0.19090909090909092_wp, &
-0.7_wp, &
-0.2361111111111111_wp, &
1.125_wp, &
-0.5_wp, &
0.23333333333333334_wp, &
-0.09545454545454546_wp, &
0.031818181818181815_wp, &
-0.008158508158508158_wp, &
0.0014985014985014985_wp, &
-0.00017482517482517483_wp, &
9.712509712509713e-06_wp, &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp, &
-9.712509712509713e-06_wp, &
0.00017482517482517483_wp, &
-0.0014985014985014985_wp, &
0.008158508158508158_wp, &
-0.031818181818181815_wp, &
0.09545454545454546_wp, &
-0.23333333333333334_wp, &
0.5_wp, &
-1.125_wp, &
0.2361111111111111_wp, &
0.7_wp, &
-0.19090909090909092_wp, &
0.05303030303030303_wp, &
-0.012237762237762238_wp, &
0.002097902097902098_wp, &
-0.0002331002331002331_wp, &
1.2487512487512488e-05_wp, &
1.2487512487512488e-05_wp, &
-0.000222000222000222_wp, &
0.0018731268731268732_wp, &
-0.00999000999000999_wp, &
0.03787878787878788_wp, &
-0.10909090909090909_wp, &
0.25_wp, &
-0.47619047619047616_wp, &
0.8035714285714286_wp, &
-1.4285714285714286_wp, &
0.478968253968254_wp, &
0.5454545454545454_wp, &
-0.11363636363636363_wp, &
0.023310023310023312_wp, &
-0.0037462537462537465_wp, &
0.0003996003996003996_wp, &
-2.0812520812520813e-05_wp, &
-2.0812520812520813e-05_wp, &
0.0003663003663003663_wp, &
-0.0030525030525030525_wp, &
0.016025641025641024_wp, &
-0.05952380952380952_wp, &
0.16666666666666666_wp, &
-0.36666666666666664_wp, &
0.6547619047619048_wp, &
-0.9821428571428571_wp, &
1.3095238095238095_wp, &
-1.8333333333333333_wp, &
0.7365440115440115_wp, &
0.4166666666666667_wp, &
-0.0641025641025641_wp, &
0.009157509157509158_wp, &
-0.0009157509157509158_wp, &
4.578754578754579e-05_wp, &
4.578754578754579e-05_wp, &
-0.0007992007992007992_wp, &
0.006593406593406593_wp, &
-0.03418803418803419_wp, &
0.125_wp, &
-0.34285714285714286_wp, &
0.7333333333333333_wp, &
-1.2571428571428571_wp, &
1.7678571428571428_wp, &
-2.0952380952380953_wp, &
2.2_wp, &
-2.4_wp, &
1.0198773448773448_wp, &
0.3076923076923077_wp, &
-0.03296703296703297_wp, &
0.0029304029304029304_wp, &
-0.00013736263736263736_wp, &
-0.00013736263736263736_wp, &
0.002380952380952381_wp, &
-0.01948051948051948_wp, &
0.1_wp, &
-0.3611111111111111_wp, &
0.975_wp, &
-2.0428571428571427_wp, &
3.4047619047619047_wp, &
-4.5964285714285715_wp, &
5.107142857142857_wp, &
-4.766666666666667_wp, &
3.9_wp, &
-3.25_wp, &
1.3468004218004217_wp, &
0.21428571428571427_wp, &
-0.014285714285714285_wp, &
0.0005952380952380953_wp, &
0.0005952380952380953_wp, &
-0.010256410256410256_wp, &
0.08333333333333333_wp, &
-0.42424242424242425_wp, &
1.5166666666666666_wp, &
-4.044444444444444_wp, &
8.341666666666667_wp, &
-13.619047619047619_wp, &
17.875_wp, &
-19.066666666666666_wp, &
16.683333333333334_wp, &
-12.133333333333333_wp, &
7.583333333333333_wp, &
-4.666666666666667_wp, &
1.7515623265623266_wp, &
0.13333333333333333_wp, &
-0.004166666666666667_wp, &
-0.004166666666666667_wp, &
0.07142857142857142_wp, &
-0.5769230769230769_wp, &
2.9166666666666665_wp, &
-10.340909090909092_wp, &
27.3_wp, &
-55.611111111111114_wp, &
89.375_wp, &
-114.91071428571429_wp, &
119.16666666666667_wp, &
-100.1_wp, &
68.25_wp, &
-37.916666666666664_wp, &
17.5_wp, &
-7.5_wp, &
2.3182289932289932_wp, &
0.0625_wp, &
0.0625_wp, &
-1.0666666666666667_wp, &
8.571428571428571_wp, &
-43.07692307692308_wp, &
151.66666666666666_wp, &
-397.09090909090907_wp, &
800.8_wp, &
-1271.111111111111_wp, &
1608.75_wp, &
-1634.2857142857142_wp, &
1334.6666666666667_wp, &
-873.6_wp, &
455.0_wp, &
-186.66666666666666_wp, &
60.0_wp, &
-16.0_wp, &
3.3807289932289932_wp /)
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:4) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (5), 5
    do i2 = 0,  -(-8) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + 8)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
        tt(1) = tt(1) + (x(i1 + 1, l + 8)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
        tt(2) = tt(2) + (x(i1 + 2, l + 8)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
        tt(3) = tt(3) + (x(i1 + 3, l + 8)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
        tt(4) = tt(4) + (x(i1 + 4, l + 8)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
      y(i1 + 2, i2) = tt(2)
      y(i1 + 3, i2) = tt(3)
      y(i1 + 4, i2) = tt(4)
    end do
    do i2 =  -(-8), n - (8) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson16_8_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson16_8_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson16_8_fil(l))
        tt(4) = tt(4) + (x(i1 + 4, l + i2)) * (poisson16_8_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
      y(i1 + 2, i2) = tt(2)
      y(i1 + 3, i2) = tt(3)
      y(i1 + 4, i2) = tt(4)
    end do
    do i2 = n - (8), n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + -8 + n - (1))) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
        tt(1) = tt(1) + (x(i1 + 1, l + -8 + n - (1))) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
        tt(2) = tt(2) + (x(i1 + 2, l + -8 + n - (1))) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
        tt(3) = tt(3) + (x(i1 + 3, l + -8 + n - (1))) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
        tt(4) = tt(4) + (x(i1 + 4, l + -8 + n - (1))) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
      end do
      y(i1 + 0, i2) = tt(0)
      y(i1 + 1, i2) = tt(1)
      y(i1 + 2, i2) = tt(2)
      y(i1 + 3, i2) = tt(3)
      y(i1 + 4, i2) = tt(4)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (5)) * (5), ndat - (1), 1
    do i2 = 0,  -(-8) - (1), 1
      tt(0) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + 8)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 =  -(-8), n - (8) - (1), 1
      tt(0) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 = n - (8), n - (1), 1
      tt(0) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + -8 + n - (1))) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
      end do
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_np_01_u5_0_false_true_false
SUBROUTINE d_poisson16_np_01_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_np_01_u1_0_true_false_true_cost
SUBROUTINE d_poisson16_np_01_a_u5_0_false_true_false(ndat, n, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  real(kind=8), intent(in), dimension(0:ndat - (1), 0:n - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat - (1), 0:n - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(0:288) :: poisson16_fil = (/ &
-3.3807289932289932_wp, &
16.0_wp, &
-60.0_wp, &
186.66666666666666_wp, &
-455.0_wp, &
873.6_wp, &
-1334.6666666666667_wp, &
1634.2857142857142_wp, &
-1608.75_wp, &
1271.111111111111_wp, &
-800.8_wp, &
397.09090909090907_wp, &
-151.66666666666666_wp, &
43.07692307692308_wp, &
-8.571428571428571_wp, &
1.0666666666666667_wp, &
-0.0625_wp, &
-0.0625_wp, &
-2.3182289932289932_wp, &
7.5_wp, &
-17.5_wp, &
37.916666666666664_wp, &
-68.25_wp, &
100.1_wp, &
-119.16666666666667_wp, &
114.91071428571429_wp, &
-89.375_wp, &
55.611111111111114_wp, &
-27.3_wp, &
10.340909090909092_wp, &
-2.9166666666666665_wp, &
0.5769230769230769_wp, &
-0.07142857142857142_wp, &
0.004166666666666667_wp, &
0.004166666666666667_wp, &
-0.13333333333333333_wp, &
-1.7515623265623266_wp, &
4.666666666666667_wp, &
-7.583333333333333_wp, &
12.133333333333333_wp, &
-16.683333333333334_wp, &
19.066666666666666_wp, &
-17.875_wp, &
13.619047619047619_wp, &
-8.341666666666667_wp, &
4.044444444444444_wp, &
-1.5166666666666666_wp, &
0.42424242424242425_wp, &
-0.08333333333333333_wp, &
0.010256410256410256_wp, &
-0.0005952380952380953_wp, &
-0.0005952380952380953_wp, &
0.014285714285714285_wp, &
-0.21428571428571427_wp, &
-1.3468004218004217_wp, &
3.25_wp, &
-3.9_wp, &
4.766666666666667_wp, &
-5.107142857142857_wp, &
4.5964285714285715_wp, &
-3.4047619047619047_wp, &
2.0428571428571427_wp, &
-0.975_wp, &
0.3611111111111111_wp, &
-0.1_wp, &
0.01948051948051948_wp, &
-0.002380952380952381_wp, &
0.00013736263736263736_wp, &
0.00013736263736263736_wp, &
-0.0029304029304029304_wp, &
0.03296703296703297_wp, &
-0.3076923076923077_wp, &
-1.0198773448773448_wp, &
2.4_wp, &
-2.2_wp, &
2.0952380952380953_wp, &
-1.7678571428571428_wp, &
1.2571428571428571_wp, &
-0.7333333333333333_wp, &
0.34285714285714286_wp, &
-0.125_wp, &
0.03418803418803419_wp, &
-0.006593406593406593_wp, &
0.0007992007992007992_wp, &
-4.578754578754579e-05_wp, &
-4.578754578754579e-05_wp, &
0.0009157509157509158_wp, &
-0.009157509157509158_wp, &
0.0641025641025641_wp, &
-0.4166666666666667_wp, &
-0.7365440115440115_wp, &
1.8333333333333333_wp, &
-1.3095238095238095_wp, &
0.9821428571428571_wp, &
-0.6547619047619048_wp, &
0.36666666666666664_wp, &
-0.16666666666666666_wp, &
0.05952380952380952_wp, &
-0.016025641025641024_wp, &
0.0030525030525030525_wp, &
-0.0003663003663003663_wp, &
2.0812520812520813e-05_wp, &
2.0812520812520813e-05_wp, &
-0.0003996003996003996_wp, &
0.0037462537462537465_wp, &
-0.023310023310023312_wp, &
0.11363636363636363_wp, &
-0.5454545454545454_wp, &
-0.478968253968254_wp, &
1.4285714285714286_wp, &
-0.8035714285714286_wp, &
0.47619047619047616_wp, &
-0.25_wp, &
0.10909090909090909_wp, &
-0.03787878787878788_wp, &
0.00999000999000999_wp, &
-0.0018731268731268732_wp, &
0.000222000222000222_wp, &
-1.2487512487512488e-05_wp, &
-1.2487512487512488e-05_wp, &
0.0002331002331002331_wp, &
-0.002097902097902098_wp, &
0.012237762237762238_wp, &
-0.05303030303030303_wp, &
0.19090909090909092_wp, &
-0.7_wp, &
-0.2361111111111111_wp, &
1.125_wp, &
-0.5_wp, &
0.23333333333333334_wp, &
-0.09545454545454546_wp, &
0.031818181818181815_wp, &
-0.008158508158508158_wp, &
0.0014985014985014985_wp, &
-0.00017482517482517483_wp, &
9.712509712509713e-06_wp, &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp, &
-9.712509712509713e-06_wp, &
0.00017482517482517483_wp, &
-0.0014985014985014985_wp, &
0.008158508158508158_wp, &
-0.031818181818181815_wp, &
0.09545454545454546_wp, &
-0.23333333333333334_wp, &
0.5_wp, &
-1.125_wp, &
0.2361111111111111_wp, &
0.7_wp, &
-0.19090909090909092_wp, &
0.05303030303030303_wp, &
-0.012237762237762238_wp, &
0.002097902097902098_wp, &
-0.0002331002331002331_wp, &
1.2487512487512488e-05_wp, &
1.2487512487512488e-05_wp, &
-0.000222000222000222_wp, &
0.0018731268731268732_wp, &
-0.00999000999000999_wp, &
0.03787878787878788_wp, &
-0.10909090909090909_wp, &
0.25_wp, &
-0.47619047619047616_wp, &
0.8035714285714286_wp, &
-1.4285714285714286_wp, &
0.478968253968254_wp, &
0.5454545454545454_wp, &
-0.11363636363636363_wp, &
0.023310023310023312_wp, &
-0.0037462537462537465_wp, &
0.0003996003996003996_wp, &
-2.0812520812520813e-05_wp, &
-2.0812520812520813e-05_wp, &
0.0003663003663003663_wp, &
-0.0030525030525030525_wp, &
0.016025641025641024_wp, &
-0.05952380952380952_wp, &
0.16666666666666666_wp, &
-0.36666666666666664_wp, &
0.6547619047619048_wp, &
-0.9821428571428571_wp, &
1.3095238095238095_wp, &
-1.8333333333333333_wp, &
0.7365440115440115_wp, &
0.4166666666666667_wp, &
-0.0641025641025641_wp, &
0.009157509157509158_wp, &
-0.0009157509157509158_wp, &
4.578754578754579e-05_wp, &
4.578754578754579e-05_wp, &
-0.0007992007992007992_wp, &
0.006593406593406593_wp, &
-0.03418803418803419_wp, &
0.125_wp, &
-0.34285714285714286_wp, &
0.7333333333333333_wp, &
-1.2571428571428571_wp, &
1.7678571428571428_wp, &
-2.0952380952380953_wp, &
2.2_wp, &
-2.4_wp, &
1.0198773448773448_wp, &
0.3076923076923077_wp, &
-0.03296703296703297_wp, &
0.0029304029304029304_wp, &
-0.00013736263736263736_wp, &
-0.00013736263736263736_wp, &
0.002380952380952381_wp, &
-0.01948051948051948_wp, &
0.1_wp, &
-0.3611111111111111_wp, &
0.975_wp, &
-2.0428571428571427_wp, &
3.4047619047619047_wp, &
-4.5964285714285715_wp, &
5.107142857142857_wp, &
-4.766666666666667_wp, &
3.9_wp, &
-3.25_wp, &
1.3468004218004217_wp, &
0.21428571428571427_wp, &
-0.014285714285714285_wp, &
0.0005952380952380953_wp, &
0.0005952380952380953_wp, &
-0.010256410256410256_wp, &
0.08333333333333333_wp, &
-0.42424242424242425_wp, &
1.5166666666666666_wp, &
-4.044444444444444_wp, &
8.341666666666667_wp, &
-13.619047619047619_wp, &
17.875_wp, &
-19.066666666666666_wp, &
16.683333333333334_wp, &
-12.133333333333333_wp, &
7.583333333333333_wp, &
-4.666666666666667_wp, &
1.7515623265623266_wp, &
0.13333333333333333_wp, &
-0.004166666666666667_wp, &
-0.004166666666666667_wp, &
0.07142857142857142_wp, &
-0.5769230769230769_wp, &
2.9166666666666665_wp, &
-10.340909090909092_wp, &
27.3_wp, &
-55.611111111111114_wp, &
89.375_wp, &
-114.91071428571429_wp, &
119.16666666666667_wp, &
-100.1_wp, &
68.25_wp, &
-37.916666666666664_wp, &
17.5_wp, &
-7.5_wp, &
2.3182289932289932_wp, &
0.0625_wp, &
0.0625_wp, &
-1.0666666666666667_wp, &
8.571428571428571_wp, &
-43.07692307692308_wp, &
151.66666666666666_wp, &
-397.09090909090907_wp, &
800.8_wp, &
-1271.111111111111_wp, &
1608.75_wp, &
-1634.2857142857142_wp, &
1334.6666666666667_wp, &
-873.6_wp, &
455.0_wp, &
-186.66666666666666_wp, &
60.0_wp, &
-16.0_wp, &
3.3807289932289932_wp /)
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: l
  real(kind=8), dimension(0:4) :: tt
!$omp parallel  default(shared) private(i1, i2, tt)
!$omp do 
  do i1 = 0, ndat - (5), 5
    do i2 = 0,  -(-8) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + 8)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
        tt(1) = tt(1) + (x(i1 + 1, l + 8)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
        tt(2) = tt(2) + (x(i1 + 2, l + 8)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
        tt(3) = tt(3) + (x(i1 + 3, l + 8)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
        tt(4) = tt(4) + (x(i1 + 4, l + 8)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
      tt(4) = (tt(4)) * (a)
      y(i1 + 4, i2) = tt(4)
    end do
    do i2 =  -(-8), n - (8) - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
        tt(1) = tt(1) + (x(i1 + 1, l + i2)) * (poisson16_8_fil(l))
        tt(2) = tt(2) + (x(i1 + 2, l + i2)) * (poisson16_8_fil(l))
        tt(3) = tt(3) + (x(i1 + 3, l + i2)) * (poisson16_8_fil(l))
        tt(4) = tt(4) + (x(i1 + 4, l + i2)) * (poisson16_8_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
      tt(4) = (tt(4)) * (a)
      y(i1 + 4, i2) = tt(4)
    end do
    do i2 = n - (8), n - (1), 1
      tt(0) = 0.0_wp
      tt(1) = 0.0_wp
      tt(2) = 0.0_wp
      tt(3) = 0.0_wp
      tt(4) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + -8 + n - (1))) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
        tt(1) = tt(1) + (x(i1 + 1, l + -8 + n - (1))) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
        tt(2) = tt(2) + (x(i1 + 2, l + -8 + n - (1))) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
        tt(3) = tt(3) + (x(i1 + 3, l + -8 + n - (1))) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
        tt(4) = tt(4) + (x(i1 + 4, l + -8 + n - (1))) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
      tt(1) = (tt(1)) * (a)
      y(i1 + 1, i2) = tt(1)
      tt(2) = (tt(2)) * (a)
      y(i1 + 2, i2) = tt(2)
      tt(3) = (tt(3)) * (a)
      y(i1 + 3, i2) = tt(3)
      tt(4) = (tt(4)) * (a)
      y(i1 + 4, i2) = tt(4)
    end do
  end do
!$omp end do 
!$omp do 
  do i1 = ((ndat) / (5)) * (5), ndat - (1), 1
    do i2 = 0,  -(-8) - (1), 1
      tt(0) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + 8)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 =  -(-8), n - (8) - (1), 1
      tt(0) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + i2)) * (poisson16_8_fil(l))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
    do i2 = n - (8), n - (1), 1
      tt(0) = 0.0_wp
      do l = -8, 8, 1
        tt(0) = tt(0) + (x(i1 + 0, l + -8 + n - (1))) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
      end do
      tt(0) = (tt(0)) * (a)
      y(i1 + 0, i2) = tt(0)
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_np_01_a_u5_0_false_true_false
SUBROUTINE d_poisson16_np_01_a_u1_0_true_false_true_cost(ndat, n, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_np_01_a_u1_0_true_false_true_cost
SUBROUTINE d_poisson16_np_201_u5_0_true_true_false(ndat1, n, ndat2, x, y)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), parameter, dimension(0:288) :: poisson16_fil = (/ &
-3.3807289932289932_wp, &
16.0_wp, &
-60.0_wp, &
186.66666666666666_wp, &
-455.0_wp, &
873.6_wp, &
-1334.6666666666667_wp, &
1634.2857142857142_wp, &
-1608.75_wp, &
1271.111111111111_wp, &
-800.8_wp, &
397.09090909090907_wp, &
-151.66666666666666_wp, &
43.07692307692308_wp, &
-8.571428571428571_wp, &
1.0666666666666667_wp, &
-0.0625_wp, &
-0.0625_wp, &
-2.3182289932289932_wp, &
7.5_wp, &
-17.5_wp, &
37.916666666666664_wp, &
-68.25_wp, &
100.1_wp, &
-119.16666666666667_wp, &
114.91071428571429_wp, &
-89.375_wp, &
55.611111111111114_wp, &
-27.3_wp, &
10.340909090909092_wp, &
-2.9166666666666665_wp, &
0.5769230769230769_wp, &
-0.07142857142857142_wp, &
0.004166666666666667_wp, &
0.004166666666666667_wp, &
-0.13333333333333333_wp, &
-1.7515623265623266_wp, &
4.666666666666667_wp, &
-7.583333333333333_wp, &
12.133333333333333_wp, &
-16.683333333333334_wp, &
19.066666666666666_wp, &
-17.875_wp, &
13.619047619047619_wp, &
-8.341666666666667_wp, &
4.044444444444444_wp, &
-1.5166666666666666_wp, &
0.42424242424242425_wp, &
-0.08333333333333333_wp, &
0.010256410256410256_wp, &
-0.0005952380952380953_wp, &
-0.0005952380952380953_wp, &
0.014285714285714285_wp, &
-0.21428571428571427_wp, &
-1.3468004218004217_wp, &
3.25_wp, &
-3.9_wp, &
4.766666666666667_wp, &
-5.107142857142857_wp, &
4.5964285714285715_wp, &
-3.4047619047619047_wp, &
2.0428571428571427_wp, &
-0.975_wp, &
0.3611111111111111_wp, &
-0.1_wp, &
0.01948051948051948_wp, &
-0.002380952380952381_wp, &
0.00013736263736263736_wp, &
0.00013736263736263736_wp, &
-0.0029304029304029304_wp, &
0.03296703296703297_wp, &
-0.3076923076923077_wp, &
-1.0198773448773448_wp, &
2.4_wp, &
-2.2_wp, &
2.0952380952380953_wp, &
-1.7678571428571428_wp, &
1.2571428571428571_wp, &
-0.7333333333333333_wp, &
0.34285714285714286_wp, &
-0.125_wp, &
0.03418803418803419_wp, &
-0.006593406593406593_wp, &
0.0007992007992007992_wp, &
-4.578754578754579e-05_wp, &
-4.578754578754579e-05_wp, &
0.0009157509157509158_wp, &
-0.009157509157509158_wp, &
0.0641025641025641_wp, &
-0.4166666666666667_wp, &
-0.7365440115440115_wp, &
1.8333333333333333_wp, &
-1.3095238095238095_wp, &
0.9821428571428571_wp, &
-0.6547619047619048_wp, &
0.36666666666666664_wp, &
-0.16666666666666666_wp, &
0.05952380952380952_wp, &
-0.016025641025641024_wp, &
0.0030525030525030525_wp, &
-0.0003663003663003663_wp, &
2.0812520812520813e-05_wp, &
2.0812520812520813e-05_wp, &
-0.0003996003996003996_wp, &
0.0037462537462537465_wp, &
-0.023310023310023312_wp, &
0.11363636363636363_wp, &
-0.5454545454545454_wp, &
-0.478968253968254_wp, &
1.4285714285714286_wp, &
-0.8035714285714286_wp, &
0.47619047619047616_wp, &
-0.25_wp, &
0.10909090909090909_wp, &
-0.03787878787878788_wp, &
0.00999000999000999_wp, &
-0.0018731268731268732_wp, &
0.000222000222000222_wp, &
-1.2487512487512488e-05_wp, &
-1.2487512487512488e-05_wp, &
0.0002331002331002331_wp, &
-0.002097902097902098_wp, &
0.012237762237762238_wp, &
-0.05303030303030303_wp, &
0.19090909090909092_wp, &
-0.7_wp, &
-0.2361111111111111_wp, &
1.125_wp, &
-0.5_wp, &
0.23333333333333334_wp, &
-0.09545454545454546_wp, &
0.031818181818181815_wp, &
-0.008158508158508158_wp, &
0.0014985014985014985_wp, &
-0.00017482517482517483_wp, &
9.712509712509713e-06_wp, &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp, &
-9.712509712509713e-06_wp, &
0.00017482517482517483_wp, &
-0.0014985014985014985_wp, &
0.008158508158508158_wp, &
-0.031818181818181815_wp, &
0.09545454545454546_wp, &
-0.23333333333333334_wp, &
0.5_wp, &
-1.125_wp, &
0.2361111111111111_wp, &
0.7_wp, &
-0.19090909090909092_wp, &
0.05303030303030303_wp, &
-0.012237762237762238_wp, &
0.002097902097902098_wp, &
-0.0002331002331002331_wp, &
1.2487512487512488e-05_wp, &
1.2487512487512488e-05_wp, &
-0.000222000222000222_wp, &
0.0018731268731268732_wp, &
-0.00999000999000999_wp, &
0.03787878787878788_wp, &
-0.10909090909090909_wp, &
0.25_wp, &
-0.47619047619047616_wp, &
0.8035714285714286_wp, &
-1.4285714285714286_wp, &
0.478968253968254_wp, &
0.5454545454545454_wp, &
-0.11363636363636363_wp, &
0.023310023310023312_wp, &
-0.0037462537462537465_wp, &
0.0003996003996003996_wp, &
-2.0812520812520813e-05_wp, &
-2.0812520812520813e-05_wp, &
0.0003663003663003663_wp, &
-0.0030525030525030525_wp, &
0.016025641025641024_wp, &
-0.05952380952380952_wp, &
0.16666666666666666_wp, &
-0.36666666666666664_wp, &
0.6547619047619048_wp, &
-0.9821428571428571_wp, &
1.3095238095238095_wp, &
-1.8333333333333333_wp, &
0.7365440115440115_wp, &
0.4166666666666667_wp, &
-0.0641025641025641_wp, &
0.009157509157509158_wp, &
-0.0009157509157509158_wp, &
4.578754578754579e-05_wp, &
4.578754578754579e-05_wp, &
-0.0007992007992007992_wp, &
0.006593406593406593_wp, &
-0.03418803418803419_wp, &
0.125_wp, &
-0.34285714285714286_wp, &
0.7333333333333333_wp, &
-1.2571428571428571_wp, &
1.7678571428571428_wp, &
-2.0952380952380953_wp, &
2.2_wp, &
-2.4_wp, &
1.0198773448773448_wp, &
0.3076923076923077_wp, &
-0.03296703296703297_wp, &
0.0029304029304029304_wp, &
-0.00013736263736263736_wp, &
-0.00013736263736263736_wp, &
0.002380952380952381_wp, &
-0.01948051948051948_wp, &
0.1_wp, &
-0.3611111111111111_wp, &
0.975_wp, &
-2.0428571428571427_wp, &
3.4047619047619047_wp, &
-4.5964285714285715_wp, &
5.107142857142857_wp, &
-4.766666666666667_wp, &
3.9_wp, &
-3.25_wp, &
1.3468004218004217_wp, &
0.21428571428571427_wp, &
-0.014285714285714285_wp, &
0.0005952380952380953_wp, &
0.0005952380952380953_wp, &
-0.010256410256410256_wp, &
0.08333333333333333_wp, &
-0.42424242424242425_wp, &
1.5166666666666666_wp, &
-4.044444444444444_wp, &
8.341666666666667_wp, &
-13.619047619047619_wp, &
17.875_wp, &
-19.066666666666666_wp, &
16.683333333333334_wp, &
-12.133333333333333_wp, &
7.583333333333333_wp, &
-4.666666666666667_wp, &
1.7515623265623266_wp, &
0.13333333333333333_wp, &
-0.004166666666666667_wp, &
-0.004166666666666667_wp, &
0.07142857142857142_wp, &
-0.5769230769230769_wp, &
2.9166666666666665_wp, &
-10.340909090909092_wp, &
27.3_wp, &
-55.611111111111114_wp, &
89.375_wp, &
-114.91071428571429_wp, &
119.16666666666667_wp, &
-100.1_wp, &
68.25_wp, &
-37.916666666666664_wp, &
17.5_wp, &
-7.5_wp, &
2.3182289932289932_wp, &
0.0625_wp, &
0.0625_wp, &
-1.0666666666666667_wp, &
8.571428571428571_wp, &
-43.07692307692308_wp, &
151.66666666666666_wp, &
-397.09090909090907_wp, &
800.8_wp, &
-1271.111111111111_wp, &
1608.75_wp, &
-1634.2857142857142_wp, &
1334.6666666666667_wp, &
-873.6_wp, &
455.0_wp, &
-186.66666666666666_wp, &
60.0_wp, &
-16.0_wp, &
3.3807289932289932_wp /)
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8), dimension(0:4) :: tt
  integer(kind=4), dimension(-8 - (8):8 - (-8) - (1)) :: mod_arr
  do l = -8 - (8), 8 - (-8) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (5), 5
      do i2 = 0,  -(-8) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        tt(4) = 0.0_wp
        do l = -8, 8, 1
          tt(0) = tt(0) + (x(i1 + 0, l + 8, i3)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
          tt(1) = tt(1) + (x(i1 + 1, l + 8, i3)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
          tt(2) = tt(2) + (x(i1 + 2, l + 8, i3)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
          tt(3) = tt(3) + (x(i1 + 3, l + 8, i3)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
          tt(4) = tt(4) + (x(i1 + 4, l + 8, i3)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
        y(i1 + 4, i2, i3) = tt(4)
      end do
      do i2 =  -(-8), n - (8) - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        tt(4) = 0.0_wp
        do l = -8, 8, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
          tt(1) = tt(1) + (x(i1 + 1, l + i2, i3)) * (poisson16_8_fil(l))
          tt(2) = tt(2) + (x(i1 + 2, l + i2, i3)) * (poisson16_8_fil(l))
          tt(3) = tt(3) + (x(i1 + 3, l + i2, i3)) * (poisson16_8_fil(l))
          tt(4) = tt(4) + (x(i1 + 4, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
        y(i1 + 4, i2, i3) = tt(4)
      end do
      do i2 = n - (8), n - (1), 1
        tt(0) = 0.0_wp
        tt(1) = 0.0_wp
        tt(2) = 0.0_wp
        tt(3) = 0.0_wp
        tt(4) = 0.0_wp
        do l = -8, 8, 1
          tt(0) = tt(0) + (x(i1 + 0, l + -8 + n - (1), i3)) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
          tt(1) = tt(1) + (x(i1 + 1, l + -8 + n - (1), i3)) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
          tt(2) = tt(2) + (x(i1 + 2, l + -8 + n - (1), i3)) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
          tt(3) = tt(3) + (x(i1 + 3, l + -8 + n - (1), i3)) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
          tt(4) = tt(4) + (x(i1 + 4, l + -8 + n - (1), i3)) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
        end do
        y(i1 + 0, i2, i3) = tt(0)
        y(i1 + 1, i2, i3) = tt(1)
        y(i1 + 2, i2, i3) = tt(2)
        y(i1 + 3, i2, i3) = tt(3)
        y(i1 + 4, i2, i3) = tt(4)
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (5)) * (5), ndat1 - (1), 1
      do i2 = 0,  -(-8) - (1), 1
        tt(0) = 0.0_wp
        do l = -8, 8, 1
          tt(0) = tt(0) + (x(i1 + 0, l + 8, i3)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 =  -(-8), n - (8) - (1), 1
        tt(0) = 0.0_wp
        do l = -8, 8, 1
          tt(0) = tt(0) + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
      do i2 = n - (8), n - (1), 1
        tt(0) = 0.0_wp
        do l = -8, 8, 1
          tt(0) = tt(0) + (x(i1 + 0, l + -8 + n - (1), i3)) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
        end do
        y(i1 + 0, i2, i3) = tt(0)
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_np_201_u5_0_true_true_false
SUBROUTINE d_poisson16_np_201_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_np_201_u1_2_true_false_true_cost
SUBROUTINE d_poisson16_np_201_a_u5_0_true_false_false(ndat1, n, ndat2, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 8
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  real(kind=8), intent(in), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: x
  real(kind=8), intent(out), dimension(0:ndat1 - (1), 0:n - (1), 0:ndat2 - (1)) :: y
  real(kind=8), intent(in) :: a
  real(kind=8), parameter, dimension(0:288) :: poisson16_fil = (/ &
-3.3807289932289932_wp, &
16.0_wp, &
-60.0_wp, &
186.66666666666666_wp, &
-455.0_wp, &
873.6_wp, &
-1334.6666666666667_wp, &
1634.2857142857142_wp, &
-1608.75_wp, &
1271.111111111111_wp, &
-800.8_wp, &
397.09090909090907_wp, &
-151.66666666666666_wp, &
43.07692307692308_wp, &
-8.571428571428571_wp, &
1.0666666666666667_wp, &
-0.0625_wp, &
-0.0625_wp, &
-2.3182289932289932_wp, &
7.5_wp, &
-17.5_wp, &
37.916666666666664_wp, &
-68.25_wp, &
100.1_wp, &
-119.16666666666667_wp, &
114.91071428571429_wp, &
-89.375_wp, &
55.611111111111114_wp, &
-27.3_wp, &
10.340909090909092_wp, &
-2.9166666666666665_wp, &
0.5769230769230769_wp, &
-0.07142857142857142_wp, &
0.004166666666666667_wp, &
0.004166666666666667_wp, &
-0.13333333333333333_wp, &
-1.7515623265623266_wp, &
4.666666666666667_wp, &
-7.583333333333333_wp, &
12.133333333333333_wp, &
-16.683333333333334_wp, &
19.066666666666666_wp, &
-17.875_wp, &
13.619047619047619_wp, &
-8.341666666666667_wp, &
4.044444444444444_wp, &
-1.5166666666666666_wp, &
0.42424242424242425_wp, &
-0.08333333333333333_wp, &
0.010256410256410256_wp, &
-0.0005952380952380953_wp, &
-0.0005952380952380953_wp, &
0.014285714285714285_wp, &
-0.21428571428571427_wp, &
-1.3468004218004217_wp, &
3.25_wp, &
-3.9_wp, &
4.766666666666667_wp, &
-5.107142857142857_wp, &
4.5964285714285715_wp, &
-3.4047619047619047_wp, &
2.0428571428571427_wp, &
-0.975_wp, &
0.3611111111111111_wp, &
-0.1_wp, &
0.01948051948051948_wp, &
-0.002380952380952381_wp, &
0.00013736263736263736_wp, &
0.00013736263736263736_wp, &
-0.0029304029304029304_wp, &
0.03296703296703297_wp, &
-0.3076923076923077_wp, &
-1.0198773448773448_wp, &
2.4_wp, &
-2.2_wp, &
2.0952380952380953_wp, &
-1.7678571428571428_wp, &
1.2571428571428571_wp, &
-0.7333333333333333_wp, &
0.34285714285714286_wp, &
-0.125_wp, &
0.03418803418803419_wp, &
-0.006593406593406593_wp, &
0.0007992007992007992_wp, &
-4.578754578754579e-05_wp, &
-4.578754578754579e-05_wp, &
0.0009157509157509158_wp, &
-0.009157509157509158_wp, &
0.0641025641025641_wp, &
-0.4166666666666667_wp, &
-0.7365440115440115_wp, &
1.8333333333333333_wp, &
-1.3095238095238095_wp, &
0.9821428571428571_wp, &
-0.6547619047619048_wp, &
0.36666666666666664_wp, &
-0.16666666666666666_wp, &
0.05952380952380952_wp, &
-0.016025641025641024_wp, &
0.0030525030525030525_wp, &
-0.0003663003663003663_wp, &
2.0812520812520813e-05_wp, &
2.0812520812520813e-05_wp, &
-0.0003996003996003996_wp, &
0.0037462537462537465_wp, &
-0.023310023310023312_wp, &
0.11363636363636363_wp, &
-0.5454545454545454_wp, &
-0.478968253968254_wp, &
1.4285714285714286_wp, &
-0.8035714285714286_wp, &
0.47619047619047616_wp, &
-0.25_wp, &
0.10909090909090909_wp, &
-0.03787878787878788_wp, &
0.00999000999000999_wp, &
-0.0018731268731268732_wp, &
0.000222000222000222_wp, &
-1.2487512487512488e-05_wp, &
-1.2487512487512488e-05_wp, &
0.0002331002331002331_wp, &
-0.002097902097902098_wp, &
0.012237762237762238_wp, &
-0.05303030303030303_wp, &
0.19090909090909092_wp, &
-0.7_wp, &
-0.2361111111111111_wp, &
1.125_wp, &
-0.5_wp, &
0.23333333333333334_wp, &
-0.09545454545454546_wp, &
0.031818181818181815_wp, &
-0.008158508158508158_wp, &
0.0014985014985014985_wp, &
-0.00017482517482517483_wp, &
9.712509712509713e-06_wp, &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp, &
-9.712509712509713e-06_wp, &
0.00017482517482517483_wp, &
-0.0014985014985014985_wp, &
0.008158508158508158_wp, &
-0.031818181818181815_wp, &
0.09545454545454546_wp, &
-0.23333333333333334_wp, &
0.5_wp, &
-1.125_wp, &
0.2361111111111111_wp, &
0.7_wp, &
-0.19090909090909092_wp, &
0.05303030303030303_wp, &
-0.012237762237762238_wp, &
0.002097902097902098_wp, &
-0.0002331002331002331_wp, &
1.2487512487512488e-05_wp, &
1.2487512487512488e-05_wp, &
-0.000222000222000222_wp, &
0.0018731268731268732_wp, &
-0.00999000999000999_wp, &
0.03787878787878788_wp, &
-0.10909090909090909_wp, &
0.25_wp, &
-0.47619047619047616_wp, &
0.8035714285714286_wp, &
-1.4285714285714286_wp, &
0.478968253968254_wp, &
0.5454545454545454_wp, &
-0.11363636363636363_wp, &
0.023310023310023312_wp, &
-0.0037462537462537465_wp, &
0.0003996003996003996_wp, &
-2.0812520812520813e-05_wp, &
-2.0812520812520813e-05_wp, &
0.0003663003663003663_wp, &
-0.0030525030525030525_wp, &
0.016025641025641024_wp, &
-0.05952380952380952_wp, &
0.16666666666666666_wp, &
-0.36666666666666664_wp, &
0.6547619047619048_wp, &
-0.9821428571428571_wp, &
1.3095238095238095_wp, &
-1.8333333333333333_wp, &
0.7365440115440115_wp, &
0.4166666666666667_wp, &
-0.0641025641025641_wp, &
0.009157509157509158_wp, &
-0.0009157509157509158_wp, &
4.578754578754579e-05_wp, &
4.578754578754579e-05_wp, &
-0.0007992007992007992_wp, &
0.006593406593406593_wp, &
-0.03418803418803419_wp, &
0.125_wp, &
-0.34285714285714286_wp, &
0.7333333333333333_wp, &
-1.2571428571428571_wp, &
1.7678571428571428_wp, &
-2.0952380952380953_wp, &
2.2_wp, &
-2.4_wp, &
1.0198773448773448_wp, &
0.3076923076923077_wp, &
-0.03296703296703297_wp, &
0.0029304029304029304_wp, &
-0.00013736263736263736_wp, &
-0.00013736263736263736_wp, &
0.002380952380952381_wp, &
-0.01948051948051948_wp, &
0.1_wp, &
-0.3611111111111111_wp, &
0.975_wp, &
-2.0428571428571427_wp, &
3.4047619047619047_wp, &
-4.5964285714285715_wp, &
5.107142857142857_wp, &
-4.766666666666667_wp, &
3.9_wp, &
-3.25_wp, &
1.3468004218004217_wp, &
0.21428571428571427_wp, &
-0.014285714285714285_wp, &
0.0005952380952380953_wp, &
0.0005952380952380953_wp, &
-0.010256410256410256_wp, &
0.08333333333333333_wp, &
-0.42424242424242425_wp, &
1.5166666666666666_wp, &
-4.044444444444444_wp, &
8.341666666666667_wp, &
-13.619047619047619_wp, &
17.875_wp, &
-19.066666666666666_wp, &
16.683333333333334_wp, &
-12.133333333333333_wp, &
7.583333333333333_wp, &
-4.666666666666667_wp, &
1.7515623265623266_wp, &
0.13333333333333333_wp, &
-0.004166666666666667_wp, &
-0.004166666666666667_wp, &
0.07142857142857142_wp, &
-0.5769230769230769_wp, &
2.9166666666666665_wp, &
-10.340909090909092_wp, &
27.3_wp, &
-55.611111111111114_wp, &
89.375_wp, &
-114.91071428571429_wp, &
119.16666666666667_wp, &
-100.1_wp, &
68.25_wp, &
-37.916666666666664_wp, &
17.5_wp, &
-7.5_wp, &
2.3182289932289932_wp, &
0.0625_wp, &
0.0625_wp, &
-1.0666666666666667_wp, &
8.571428571428571_wp, &
-43.07692307692308_wp, &
151.66666666666666_wp, &
-397.09090909090907_wp, &
800.8_wp, &
-1271.111111111111_wp, &
1608.75_wp, &
-1634.2857142857142_wp, &
1334.6666666666667_wp, &
-873.6_wp, &
455.0_wp, &
-186.66666666666666_wp, &
60.0_wp, &
-16.0_wp, &
3.3807289932289932_wp /)
  real(kind=8), parameter, dimension(-8:8) :: poisson16_8_fil = (/ &
9.712509712509713e-06_wp, &
-0.0001776001776001776_wp, &
0.001554001554001554_wp, &
-0.008702408702408702_wp, &
0.03535353535353535_wp, &
-0.11313131313131314_wp, &
0.3111111111111111_wp, &
-0.8888888888888888_wp, &
0.0_wp, &
0.8888888888888888_wp, &
-0.3111111111111111_wp, &
0.11313131313131314_wp, &
-0.03535353535353535_wp, &
0.008702408702408702_wp, &
-0.001554001554001554_wp, &
0.0001776001776001776_wp, &
-9.712509712509713e-06_wp /)
  integer(kind=4) :: i1
  integer(kind=4) :: i2
  integer(kind=4) :: i3
  integer(kind=4) :: l
  real(kind=8) :: tt0
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
  integer(kind=4), dimension(-8 - (8):8 - (-8) - (1)) :: mod_arr
  do l = -8 - (8), 8 - (-8) - (1), 1
    mod_arr(l) = modulo(l, n)
  end do
!$omp parallel  default(shared) private(i1, i2, i3, tt0, tt1, tt2, tt3, tt4)
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = 0, ndat1 - (5), 5
      do i2 = 0,  -(-8) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt4 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, l + 8, i3)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
          tt1 = tt1 + (x(i1 + 1, l + 8, i3)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
          tt2 = tt2 + (x(i1 + 2, l + 8, i3)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
          tt3 = tt3 + (x(i1 + 3, l + 8, i3)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
          tt4 = tt4 + (x(i1 + 4, l + 8, i3)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
        tt4 = (tt4) * (a)
        y(i1 + 4, i2, i3) = tt4
      end do
      do i2 =  -(-8), n - (8) - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt4 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
          tt1 = tt1 + (x(i1 + 1, l + i2, i3)) * (poisson16_8_fil(l))
          tt2 = tt2 + (x(i1 + 2, l + i2, i3)) * (poisson16_8_fil(l))
          tt3 = tt3 + (x(i1 + 3, l + i2, i3)) * (poisson16_8_fil(l))
          tt4 = tt4 + (x(i1 + 4, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
        tt4 = (tt4) * (a)
        y(i1 + 4, i2, i3) = tt4
      end do
      do i2 = n - (8), n - (1), 1
        tt0 = 0.0_wp
        tt1 = 0.0_wp
        tt2 = 0.0_wp
        tt3 = 0.0_wp
        tt4 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, l + -8 + n - (1), i3)) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
          tt1 = tt1 + (x(i1 + 1, l + -8 + n - (1), i3)) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
          tt2 = tt2 + (x(i1 + 2, l + -8 + n - (1), i3)) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
          tt3 = tt3 + (x(i1 + 3, l + -8 + n - (1), i3)) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
          tt4 = tt4 + (x(i1 + 4, l + -8 + n - (1), i3)) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
        tt1 = (tt1) * (a)
        y(i1 + 1, i2, i3) = tt1
        tt2 = (tt2) * (a)
        y(i1 + 2, i2, i3) = tt2
        tt3 = (tt3) * (a)
        y(i1 + 3, i2, i3) = tt3
        tt4 = (tt4) * (a)
        y(i1 + 4, i2, i3) = tt4
      end do
    end do
  end do
!$omp end do 
!$omp do 
  do i3 = 0, ndat2 - (1), 1
    do i1 = ((ndat1) / (5)) * (5), ndat1 - (1), 1
      do i2 = 0,  -(-8) - (1), 1
        tt0 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, l + 8, i3)) * (poisson16_fil((i2 - (8)) * (17) + 136 + l + 8))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 =  -(-8), n - (8) - (1), 1
        tt0 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, l + i2, i3)) * (poisson16_8_fil(l))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
      do i2 = n - (8), n - (1), 1
        tt0 = 0.0_wp
        do l = -8, 8, 1
          tt0 = tt0 + (x(i1 + 0, l + -8 + n - (1), i3)) * (poisson16_fil((i2 + 8 - (n) + 1) * (17) + 136 + l + 8))
        end do
        tt0 = (tt0) * (a)
        y(i1 + 0, i2, i3) = tt0
      end do
    end do
  end do
!$omp end do 
!$omp end parallel 
END SUBROUTINE d_poisson16_np_201_a_u5_0_true_false_false
SUBROUTINE d_poisson16_np_201_a_u1_2_true_false_true_cost(ndat1, n, ndat2, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: ndat1
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat2
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: ndat_t
  ndat_t = 1
  ndat_t = (ndat_t) * (ndat2)
  ndat_t = (ndat_t) * (ndat1)
  cost = (((n) * (2)) * (17)) * (ndat_t)
END SUBROUTINE d_poisson16_np_201_a_u1_2_true_false_true_cost
SUBROUTINE d_s0s0_1d_poisson16_cost(d, idim, n, bc, x, y, a, cost)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: d
  integer(kind=4), intent(in) :: idim
  integer(kind=4), intent(in), dimension(0:d - (1)) :: n
  integer(kind=4), intent(in) :: bc
  real(kind=8), intent(in), dimension(*) :: x
  real(kind=8), intent(out), dimension(*) :: y
  real(kind=8), intent(in) :: a
  integer(kind=4), intent(out) :: cost
  integer(kind=4) :: i
  integer(kind=4) :: ndat_left
  integer(kind=4) :: ndat_right
  integer(kind=4) :: c
  if (idim == 0) then
    ndat_right = 1
    do i = 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson16_p_10_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson16_p_10_a_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson16_fg_10_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson16_fg_10_a_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson16_fs_10_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson16_fs_10_a_u1_1_false_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-2)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson16_np_10_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson16_np_10_a_u1_1_true_false_true_cost(n(idim), ndat_right,  c)
          cost = cost + c
        end if
      end select
  else if (idim == d - (1)) then
    ndat_left = 1
    do i = 0, d - (2), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson16_p_01_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson16_p_01_a_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      case (1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson16_fg_01_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson16_fg_01_a_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      case (-1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson16_fs_01_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson16_fs_01_a_u1_0_false_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      case (-2)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson16_np_01_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson16_np_01_a_u1_0_true_false_true_cost(ndat_left, n(idim),  c)
          cost = cost + c
        end if
      end select
  else
    ndat_left = 1
    ndat_right = 1
    do i = 0, idim - (1), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    do i = idim + 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson16_p_201_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson16_p_201_a_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson16_fg_201_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson16_fg_201_a_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-1)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson16_fs_201_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson16_fs_201_a_u1_2_false_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      case (-2)
        if (a == 1.0_wp) then
          cost = 0
          call d_poisson16_np_201_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        else
          cost = 0
          call d_poisson16_np_201_a_u1_2_true_false_true_cost(ndat_left, n(idim), ndat_right,  c)
          cost = cost + c
        end if
      end select
  end if
END SUBROUTINE d_s0s0_1d_poisson16_cost
SUBROUTINE d_s0s0_1d_poisson16(d, idim, n, bc, x, y, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: d
  integer(kind=4), intent(in) :: idim
  integer(kind=4), intent(in), dimension(0:d - (1)) :: n
  integer(kind=4), intent(in) :: bc
  real(kind=8), intent(in), dimension(*) :: x
  real(kind=8), intent(out), dimension(*) :: y
  real(kind=8), intent(in) :: a
  integer(kind=4) :: i
  integer(kind=4) :: ndat_left
  integer(kind=4) :: ndat_right
  if (idim == 0) then
    ndat_right = 1
    do i = 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          call d_poisson16_p_10_u2_1_true_true_false(n(idim), ndat_right, x, y)
        else
          call d_poisson16_p_10_a_u2_1_true_false_false(n(idim), ndat_right, x, y, a)
        end if
      case (1)
        if (a == 1.0_wp) then
          call d_poisson16_fg_10_u1_1_false_true_true(n(idim), ndat_right, x, y)
        else
          call d_poisson16_fg_10_a_u1_1_false_true_true(n(idim), ndat_right, x, y, a)
        end if
      case (-1)
        if (a == 1.0_wp) then
          call d_poisson16_fs_10_u2_1_false_false_false(n(idim), ndat_right, x, y)
        else
          call d_poisson16_fs_10_a_u2_1_false_false_false(n(idim), ndat_right, x, y, a)
        end if
      case (-2)
        if (a == 1.0_wp) then
          call d_poisson16_np_10_u1_1_false_false_true(n(idim), ndat_right, x, y)
        else
          call d_poisson16_np_10_a_u1_1_true_false_false(n(idim), ndat_right, x, y, a)
        end if
      end select
  else if (idim == d - (1)) then
    ndat_left = 1
    do i = 0, d - (2), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          call d_poisson16_p_01_u4_0_true_false_false(ndat_left, n(idim), x, y)
        else
          call d_poisson16_p_01_a_u4_0_true_false_false(ndat_left, n(idim), x, y, a)
        end if
      case (1)
        if (a == 1.0_wp) then
          call d_poisson16_fg_01_u5_0_false_true_false(ndat_left, n(idim), x, y)
        else
          call d_poisson16_fg_01_a_u5_0_false_true_false(ndat_left, n(idim), x, y, a)
        end if
      case (-1)
        if (a == 1.0_wp) then
          call d_poisson16_fs_01_u5_0_false_true_false(ndat_left, n(idim), x, y)
        else
          call d_poisson16_fs_01_a_u5_0_false_false_false(ndat_left, n(idim), x, y, a)
        end if
      case (-2)
        if (a == 1.0_wp) then
          call d_poisson16_np_01_u5_0_false_true_false(ndat_left, n(idim), x, y)
        else
          call d_poisson16_np_01_a_u5_0_false_true_false(ndat_left, n(idim), x, y, a)
        end if
      end select
  else
    ndat_left = 1
    ndat_right = 1
    do i = 0, idim - (1), 1
      ndat_left = (ndat_left) * (n(i))
    end do
    do i = idim + 1, d - (1), 1
      ndat_right = (ndat_right) * (n(i))
    end do
    select case (bc)
      case (0)
        if (a == 1.0_wp) then
          call d_poisson16_p_201_u4_0_true_true_false(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson16_p_201_a_u4_0_true_false_false(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      case (1)
        if (a == 1.0_wp) then
          call d_poisson16_fg_201_u5_0_false_false_false(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson16_fg_201_a_u5_0_false_false_false(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      case (-1)
        if (a == 1.0_wp) then
          call d_poisson16_fs_201_u5_0_false_false_false(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson16_fs_201_a_u5_0_false_true_false(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      case (-2)
        if (a == 1.0_wp) then
          call d_poisson16_np_201_u5_0_true_true_false(ndat_left, n(idim), ndat_right, x, y)
        else
          call d_poisson16_np_201_a_u5_0_true_false_false(ndat_left, n(idim), ndat_right, x, y, a)
        end if
      end select
  end if
END SUBROUTINE d_s0s0_1d_poisson16
SUBROUTINE Poisson_broker(nord, idim, nn, bc, u, du, a)
  integer, parameter :: wp=kind(1.0d0)
  integer(kind=4), intent(in) :: nord
  integer(kind=4), intent(in) :: idim
  integer(kind=4), intent(in), dimension(3) :: nn
  integer(kind=4), intent(in) :: bc
  real(kind=8), intent(in), dimension(*) :: u
  real(kind=8), intent(out), dimension(*) :: du
  real(kind=8), intent(in) :: a
  select case (nord)
    case (2)
      call d_s0s0_1d_poisson2(3, idim, nn, bc, u, du, a)
    case (4)
      call d_s0s0_1d_poisson4(3, idim, nn, bc, u, du, a)
    case (6)
      call d_s0s0_1d_poisson6(3, idim, nn, bc, u, du, a)
    case (8)
      call d_s0s0_1d_poisson8(3, idim, nn, bc, u, du, a)
    case (16)
      call d_s0s0_1d_poisson16(3, idim, nn, bc, u, du, a)
    end select
END SUBROUTINE Poisson_broker
