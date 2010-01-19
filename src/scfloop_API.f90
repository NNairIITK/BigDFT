!!****m* BigDFT/scfloop_API
!! FUNCTION
!!  Self-Consistent Loop API
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2009 CEA, UNIBAS
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
module scfloop_API

  use module_base
  use module_types

  implicit none

  ! Storage of required variables for a SCF loop calculation.
  logical :: scfloop_initialised = .false.
  integer :: scfloop_nproc
  type(atoms_data), pointer :: scfloop_at
  type(input_variables), pointer :: scfloop_in
  type(restart_objects), pointer :: scfloop_rst

  public :: scfloop_init
!!!  public :: scfloop_finalise
contains

  subroutine scfloop_init(nproc_, at_, in_, rst_)
    integer, intent(in) :: nproc_
    type(atoms_data), intent(in), target :: at_
    type(input_variables), intent(in), target :: in_
    type(restart_objects), intent(in), target :: rst_

    scfloop_nproc = nproc_
    scfloop_at => at_
    scfloop_in => in_
    scfloop_rst => rst_

    scfloop_initialised = .true.
  end subroutine scfloop_init

!!!  subroutine scfloop_finalise()
!!!  end subroutine scfloop_finalise
end module scfloop_API
!!***

subroutine scfloop_main(acell, epot, fcart, grad, itime, me, natom, rprimd, xred)
  use scfloop_API
  use module_base
  use module_types
  use module_interfaces

  implicit none

  integer, intent(in) :: natom, itime, me
  real(dp), intent(out) :: epot
  real(dp), intent(in) :: acell(3)
  real(dp), intent(in) :: rprimd(3,3), xred(3,natom)
  real(dp), intent(out) :: fcart(3, natom), grad(3, natom)

  character(len=*), parameter :: subname='scfloop_main'
  integer :: infocode, i, i_stat, i_all
  real(dp) :: favg(3)
  real(dp), allocatable :: xcart(:,:)

  if (.not. scfloop_initialised) then
     write(0,*) "No previous call to scfloop_init(). On strike, refuse to work."
     stop
  end if

  if (me == 0) then
     write( *,'(1x,a,1x,i0)') &
          & 'SCFloop API, call force calculation step=', itime
  end if

  ! We transfer acell into at
  scfloop_at%alat1 = acell(1)
  scfloop_at%alat2 = rprimd(2,2)
  scfloop_at%alat3 = acell(3)

  scfloop_in%inputPsiId=1
  ! need to transform xred into xcart
  allocate(xcart(3, scfloop_at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,xcart,'xcart',subname)
  do i = 1, scfloop_at%nat, 1
     xcart(:, i) = xred(:, i) * acell(:)
  end do

  scfloop_in%inputPsiId = 1
  call call_bigdft(scfloop_nproc,me,scfloop_at,xcart,scfloop_in,epot,grad,scfloop_rst,infocode)

  ! need to transform the forces into reduced ones.
  favg(:) = real(0, dp)
  do i = 1, scfloop_at%nat, 1
     fcart(:, i) = grad(:, i)
     favg(:) = favg(:) + fcart(:, i) / real(natom, dp)
     grad(:, i) = -grad(:, i) / acell(:)
  end do
  do i = 1, scfloop_at%nat, 1
     fcart(:, i) = fcart(:, i) - favg(:)
  end do

  i_all=-product(shape(xcart))*kind(xcart)
  deallocate(xcart,stat=i_stat)
  call memocc(i_stat,i_all,'xcart',subname)
end subroutine scfloop_main

subroutine scfloop_output(acell, epot, ekin, fred, itime, me, natom, rprimd, vel, xred)
  use scfloop_API
  use module_base
  use module_types
  use module_interfaces

  implicit none

  integer, intent(in) :: natom, itime, me
  real(dp), intent(in) :: epot, ekin
  real(dp), intent(in) :: acell(3)
  real(dp), intent(in) :: rprimd(3,3), xred(3,natom)
  real(dp), intent(in) :: fred(3, natom), vel(3, natom)

  character(len=*), parameter :: subname='scfloop_output'
  character(len = 4) :: fn4
  character(len = 40) :: comment
  integer :: i, i_stat, i_all
  real :: fnrm
  real(dp), allocatable :: xcart(:,:)

  if (me /= 0) return

  fnrm = real(0, dp)
  ! need to transform xred into xcart
  allocate(xcart(3, scfloop_at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,xcart,'xcart',subname)
  do i = 1, scfloop_at%nat, 1
     xcart(:, i) = xred(:, i) * acell(:)
     fnrm = fnrm + fred(1, i) * acell(1) * fred(1, i) * acell(1) + &
          & fred(2, i) * acell(2) * fred(2, i) * acell(2) + &
          & fred(3, i) * acell(3) * fred(3, i) * acell(3)
  end do

  write(fn4,'(i4.4)') itime
  write(comment,'(a,1pe10.3)')'AB6MD:fnrm= ', sqrt(fnrm)
  call write_atomic_file('posout_'//fn4, epot + ekin, xcart, scfloop_at, trim(comment))

  i_all=-product(shape(xcart))*kind(xcart)
  deallocate(xcart,stat=i_stat)
  call memocc(i_stat,i_all,'xcart',subname)
end subroutine scfloop_output
