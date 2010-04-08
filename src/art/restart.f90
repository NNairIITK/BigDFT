!!****f* art/restart
!! FUNCTION
!!   This subroutine is used to restart an event during the activation. It uses the files:
!!     restart.dat
!!   It then continues the event where it stopped. 
!!
!! COPYRIGHT
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!! 
subroutine restart_states( istatus, ieventcurrent, iterations )

  use defs
  implicit none  

  !Arguments
  integer, intent(out) :: istatus
  integer, intent(out) :: ieventcurrent
  integer, intent(out) :: iterations 

  !Local variables
  integer :: i, ierror
  integer :: nat
  real(8), dimension(:), pointer :: dx, dy, dz   ! Pointers for reference position

  dx => direction_restart(1:NATOMS)
  dy => direction_restart(NATOMS+1:2*NATOMS)
  dz => direction_restart(2*NATOMS+1:3*NATOMS)

  open(unit=FRESTART,file=RESTARTFILE,status='old',action='read',iostat=ierror)
  read(FRESTART,*) istatus, iterations, mincounter
  read(FRESTART,*) ieventcurrent, evalf_number
  read(FRESTART,*) box, scala, boxref, scalaref
  read(FRESTART,*) total_energy, ref_energy
  read(FRESTART,*) ( type(i), i = 1, NATOMS )
  read(FRESTART,*) ( xref(i), yref(i), zref(i), i = 1, NATOMS )
  read(FRESTART,*) ( x(i), y(i), z(i), i = 1, NATOMS )
  read(FRESTART,*) ( dx(i), dy(i), dz(i), i = 1, NATOMS )
  close(FRESTART)

                                      ! Write
  if ( iproc == 0 ) then 
   write(*,*) 'BART: restart file'
   write(*,*) 'BART: box: ', box
   write(*,'(a,(1p,e17.10,0p))') ' BART: total energy: ', total_energy
   write(*,'(a,3f20.8)') ' BART: xref: ', xref(1), yref(1), zref(1)
   write(*,'(a,3f20.8)') ' BART: x   : ', x(1), y(1), z(1)
  end if

  nat = 3 * NATOMS
  call MPI_Bcast(pos,nat,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  call MPI_Bcast(posref,nat,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

END SUBROUTINE restart_states 
!!***

!!****f* restart/save_state
!! FUNCTION
!!   saves the status for a restart
!!
!! SOURCE
!!
subroutine save_state( istatus, iter, direction )

  use defs
  implicit none  

  !Arguments
  integer, intent(in) :: istatus
  integer, intent(in) :: iter
  real(8), dimension(VECSIZE), target :: direction

  !Local variables
  integer :: i, ierror
  real(8), dimension(:), pointer :: dx, dy, dz   ! Pointers for reference position

  dx => direction(1:NATOMS)
  dy => direction(NATOMS+1:2*NATOMS)
  dz => direction(2*NATOMS+1:3*NATOMS)

  open(unit=FRESTART,file=RESTARTFILE,status='unknown',action='write',iostat=ierror)
  write(FRESTART,*) istatus, iter, mincounter
  write(FRESTART,*) ievent, evalf_number
  write(FRESTART,*) box, scala, boxref, scalaref
  write(FRESTART,*) total_energy, ref_energy
  write(FRESTART,'(2x,10i6)') ( type(i), i = 1, NATOMS )
  write(FRESTART,'(2x,3f20.8)') ( xref(i), yref(i), zref(i), i = 1, NATOMS )
  write(FRESTART,'(2x,3f20.8)') ( x(i), y(i), z(i), i = 1, NATOMS )
  write(FRESTART,'(2x,3f20.8)') ( dx(i), dy(i), dz(i), i = 1, NATOMS )
  close(FRESTART)

END SUBROUTINE save_state
!!***
