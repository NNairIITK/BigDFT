!> @file
!! Test of the libxc library
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Program to test the libxc library used by BigDFT
program test_xc

  use module_base
  use module_xc

  implicit none

  integer, parameter :: n_funcs = 54
  integer, dimension(n_funcs), parameter :: funcs = (/ &
       & 1, -20, 2, -1009, 3, 4, -1002, 5, -1004, 6, -6006, &
       & 7, -1012, 8, -1, 9, -1003, -1005, -1007, -1008, -1010, &
       & -1011, -1013, -1014, &
       & 11, -101130, 12, -101, 13, -160012, 14, -102130, 15, -117130, &
       & 16, -161, 17, -162, 23, -118, 26, -163, 27, -164, &
       & -102, -103, -104, -105, -106, -107, -108, -109, -110, &
       & -406 /)
!!$  integer, parameter :: n_funcs = 1
!!$  integer, dimension(n_funcs), parameter :: funcs = (/ -101130 /)
  integer :: ifunc, ixc_prev, ierr, iproc, nproc
  real(dp) :: exc_(2, n_funcs), dt_(n_funcs)
  real(dp) :: exc(2, n_funcs), dt(n_funcs)
  

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  
  exc_ = 0.d0
  dt_  = 0.d0
  do ifunc = 1, n_funcs, 1
     if (modulo(ifunc, nproc) == iproc) then
        call test(funcs(ifunc), exc_(:, ifunc), dt_(ifunc))
     end if
!!$     if (funcs(ifunc) < 0) then
!!$        call test(funcs(ifunc), exc, dt, option = XC_LIBXC)
!!$        write(*,"(1x,A,I7,3x,A,F17.8,1x,A,1x,A,F17.8,3x,A,F10.5,1x,A)") &
!!$             & "ixc = ", funcs(ifunc), "nosp = ", exc(1), "|", "scol = ", &
!!$             & exc(2), "time = ", dt, "s"
!!$     end if
  end do
  if (nproc == 1) then
     exc = exc_
     dt = dt_
  else
     call MPI_ALLREDUCE(exc_, exc, 2 * n_funcs, &
          & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(dt_, dt, n_funcs, &
          & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  end if

  if (iproc == 0) then
     ixc_prev = 1
     do ifunc = 1, n_funcs, 1
        if (ixc_prev * funcs(ifunc) > 0 .or. (ixc_prev < 0 .and. funcs(ifunc) > 0)) &
             & write(*,"(1x,A,A,A)") repeat("-", 41), "+", repeat("-", 44)
        write(*,"(1x,A,I7,3x,A,F17.8,1x,A,1x,A,F17.8,3x,A,F10.5,1x,A)") &
             & "ixc = ", funcs(ifunc), "nosp = ", exc(1, ifunc), "|", "scol = ", &
             & exc(2, ifunc), "time = ", dt(ifunc), "s"
        ixc_prev = funcs(ifunc)
     end do
     write(*,"(1x,A,A,A)") repeat("-", 41), "+", repeat("-", 44)
  end if

  call MPI_FINALIZE(ierr)

contains

  subroutine test(ixc, excs, dt, option)
    use module_base
    use module_xc

    implicit none

    integer, intent(in) :: ixc
    real(dp), intent(out) :: excs(2), dt
    integer, intent(in), optional :: option

    integer :: i, n, type
    integer, parameter :: n_rho = 100000, n_runs = 2
    real(dp), dimension(n_rho, 2) :: rho, vxc
    real(dp), dimension(n_rho, 3) :: rhogr, vxcgr
    real(dp), dimension(n_rho) :: exc
    integer :: start, end, countPerSecond

    if (present(option)) then
       type = option
    else
       if (ixc < 0) then
          type = XC_MIXED
       else
          type = XC_ABINIT
       end if
    end if
    call system_clock(count_rate = countPerSecond)
    call system_clock(start)
    do n = 1, n_runs, 1
       do i = 1, 2, 1
          call xc_init(ixc, type, i)
!!$          if (i == 1 .and. n == 1) call xc_dump()

          call gauss(rho, n_rho, i, type)
          if (xc_isgga()) call gaussgr(rhogr, rho, n_rho)
          call xc_getvxc(n_rho, exc, i, rho(1,1), vxc(1,1), rhogr, vxcgr)
          excs(i) = sum(exc)

          call xc_end()
       end do
    end do
    call system_clock(end)

    dt = real(end - start) / real(countPerSecond) / real(n_runs)
  end subroutine test

  subroutine gauss(rho, n_rho, nspin, type)
    use module_base
    use module_xc

    implicit none

    integer, intent(in) :: n_rho, nspin, type
    real(dp), intent(out) :: rho(n_rho * 2)

    integer :: j, delta
    real(dp) :: sigma

    call xc_init_rho(n_rho * 2, rho, 1)
    delta = 0
    if (nspin == 2 ) delta = int(real(n_rho) * 0.005)
    sigma = 1.d0 / (real(n_rho, dp) * 0.25d0)
    do j = 5, n_rho - 5
       if (type == XC_LIBXC .and. nspin == 2) then
          rho(2 * j - 1) = exp(-((sigma * (n_rho / 2 - j + delta)) ** 2)) / nspin / n_rho
          rho(2 * j    ) = exp(-((sigma * (n_rho / 2 - j - delta)) ** 2)) / nspin / n_rho
       else
          rho(j        ) = exp(-((sigma * (n_rho / 2 - j + delta)) ** 2)) / nspin / n_rho
          rho(j + n_rho) = exp(-((sigma * (n_rho / 2 - j - delta)) ** 2)) / nspin / n_rho
       end if
    end do
    call xc_clean_rho(n_rho * 2, rho, 1)
    rho = rho / sum(rho(1:nspin * n_rho))
  end subroutine gauss

  subroutine gaussgr(rhogr, rho, n_rho)
    implicit none

    integer, intent(in) :: n_rho
    real(dp), intent(in) :: rho(n_rho, 2)
    real(dp), intent(out) :: rhogr(n_rho, 3)

    integer :: j

    rhogr(1, :) = 0.d0
    do j = 2, n_rho - 1
       rhogr(j, 1) = (rho(j + 1, 1) - rho(j - 1, 1)) / 2. * n_rho
       rhogr(j, 2) = (rho(j + 1, 2) - rho(j - 1, 2)) / 2. * n_rho
       rhogr(j, 3) = rhogr(j, 1) + rhogr(j, 2)
       rhogr(j, 1) = rhogr(j, 1) * rhogr(j, 1)
       rhogr(j, 2) = rhogr(j, 2) * rhogr(j, 2)
       rhogr(j, 3) = rhogr(j, 3) * rhogr(j, 3)
    end do
    rhogr(n_rho, :) = 0.d0
  end subroutine gaussgr
end program test_xc
